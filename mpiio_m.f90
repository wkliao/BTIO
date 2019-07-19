!
!  Copyright (C) 2013, Northwestern University
!  See COPYRIGHT notice in top-level directory.
!

      module mpiio_m
      use mpi
      use header
      implicit none

      character(len=512) mpi_filenm
      integer fp, combined_btype, elemtype, combined_ftype, filetype

      contains

      !----< mpiio_setup >---------------------------------------------
      integer function mpiio_setup(io_mode)
      implicit none

      character io_mode   ! 'w' for write and 'r' for read

      ! local variables
      integer c, omode, info, err, elemsize, err_len
      integer sizes(4), starts(4), subsizes(4)
      integer(KIND=MPI_ADDRESS_KIND) lb, gsize
      integer(KIND=MPI_OFFSET_KIND) iseek
      integer, allocatable :: cell_btype(:)
      integer, allocatable :: cell_ftype(:)
      integer, allocatable :: cell_blength(:)
      integer, allocatable :: cell_disp(:)
      character(LEN=MPI_MAX_ERROR_STRING) err_string

      mpiio_setup = 1

      allocate(cell_btype(ncells))
      allocate(cell_ftype(ncells))
      allocate(cell_blength(ncells))
      allocate(cell_disp(ncells))

      mpi_filenm = trim(dir_path)//'/btio.mpi'

      call MPI_Type_contiguous(5, MPI_DOUBLE_PRECISION, elemtype, err)
      call MPI_Type_commit(elemtype, err)
      call MPI_Type_size(elemtype, elemsize, err)

      ! define an MPI derived data type for local buffer (ghosted)
      ! ncells is the number of subarrays written by each process (the square
      ! root of nproc
      do c = 1, ncells

         ! Outer array dimensions are same for every cell
         sizes(1) = IMAX+4
         sizes(2) = JMAX+4
         sizes(3) = KMAX+4

         ! 4th dimension is number of cells
         sizes(4) = ncells

         ! dimensions of local cells, can differ slightly between cells
         subsizes(1) = cell_size(1, c)
         subsizes(2) = cell_size(2, c)
         subsizes(3) = cell_size(3, c)

         ! Cell is 4th dimension, 1 cell per cell type to handle varying
         ! cell sub-array sizes
         subsizes(4) = 1

         ! local buffer MPI derived type constructors use 0-based start addresses
         starts(1) = 2
         starts(2) = 2
         starts(3) = 2
         starts(4) = c-1

         ! Create buftype for a local cell
         call MPI_Type_create_subarray(4, sizes, subsizes, starts, &
              MPI_ORDER_FORTRAN, elemtype, cell_btype(c), err)
         call MPI_Type_commit(cell_btype(c), err)

         ! block length and displacement for joining cells -
         ! 1 cell buftype per block, cell buftypes have own displacement
         ! generated from cell number (4th array dimension)
         cell_blength(c) = 1
         cell_disp(c) = 0
      enddo

      ! Create combined/concatenated  buftype for all local cells
      call MPI_Type_struct(ncells, cell_blength, cell_disp, &
                           cell_btype, combined_btype, err)
      call MPI_Type_commit(combined_btype, err)
      do c = 1, ncells
         call MPI_Type_free(cell_btype(c), err)
      enddo

      ! define an MPI derived data type for fileview
      do c = 1, ncells
         sizes(1) = grid_points(1)
         sizes(2) = grid_points(2)
         sizes(3) = grid_points(3)

         ! Size of c'th cell
         subsizes(1) = cell_size(1, c)
         subsizes(2) = cell_size(2, c)
         subsizes(3) = cell_size(3, c)

         ! Starting point in full array of c'th cell
         starts(1) = cell_low(1,c)
         starts(2) = cell_low(2,c)
         starts(3) = cell_low(3,c)

         call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
                                       MPI_ORDER_FORTRAN, &
                                       elemtype, cell_ftype(c), err)
         call MPI_Type_commit(cell_ftype(c), err)
         cell_blength(c) = 1
         cell_disp(c) = 0
      enddo

      ! concatenate into the final file type
      lb = 0
      gsize = elemsize &
            * grid_points(1) * grid_points (2) * grid_points(3)

      call MPI_Type_struct(ncells, cell_blength, cell_disp, &
                           cell_ftype, combined_ftype, err)
      call MPI_Type_commit(combined_ftype, err)
      ! expand the type size to the entire global array
      call MPI_Type_create_resized(combined_ftype, lb, gsize, &
                                   filetype, err)
      call MPI_Type_commit(filetype, err)
      call MPI_Type_free(combined_ftype, err)
      do c = 1, ncells
         call MPI_Type_free(cell_ftype(c), err)
      enddo

      deallocate(cell_btype)
      deallocate(cell_ftype)
      deallocate(cell_blength)
      deallocate(cell_disp)

      if (io_mode .EQ. 'w') then
         omode = IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE)
!         if (rank .EQ. root) &
!            call MPI_File_delete(mpi_filenm, MPI_INFO_NULL, err)
!         call MPI_Barrier(MPI_COMM_WORLD, err)
      else
         omode = MPI_MODE_RDONLY
      endif

      !      disable data sieving for write operations
      call MPI_Info_create(info, err)
      call MPI_Info_set(info, 'romio_ds_write', 'disable', err)

      call MPI_File_open(MPI_COMM_WORLD, trim(mpi_filenm), omode, &
                         info, fp, err)
      if (err .ne. MPI_SUCCESS) then
         call MPI_Error_string(err, err_string, err_len, err)
         print *, 'Error opening file: ',trim(err_string)
         mpiio_setup = 0
         return
      endif
      call MPI_Info_free(info, err)

      iseek=0  ! type of MPI_OFFSET_KIND
      call MPI_File_set_view(fp, iseek, MPI_BYTE, filetype, &
                             'native', MPI_INFO_NULL, err)
      if (err .ne. MPI_SUCCESS) then
         call MPI_Error_string(err, err_string, err_len, err)
         print *, 'Error setting file view: ',trim(err_string)
         mpiio_setup = 0
         return
      endif

      ! get MPI-IO info object
      call MPI_File_get_info(fp, info_used, err)

      end function mpiio_setup

      !----< mpiio_write >--------------------------------------------
      subroutine mpiio_write(io_method)
      implicit none
      integer io_method

      ! local variables
      integer err, err_len, mstatus(MPI_STATUS_SIZE)
      character(LEN=MPI_MAX_ERROR_STRING) err_string

      num_io = num_io +  1
      if (io_method .EQ. 0) then
          call MPI_File_write_all(fp, u, 1, combined_btype, mstatus,err)
      else
          call MPI_File_write(fp, u, 1, combined_btype, mstatus, err)
      endif
      if (err .ne. MPI_SUCCESS) then
          call MPI_Error_string(err, err_string, err_len, err)
          print *, 'Error: MPI write to file: ',trim(err_string)
          stop
      endif

      end subroutine mpiio_write

      !----< mpiio_read >--------------------------------------------
      subroutine mpiio_read(io_method)
      implicit none
      integer io_method

      ! local variables
      integer err, err_len, mstatus(MPI_STATUS_SIZE)
      character(LEN=MPI_MAX_ERROR_STRING) err_string

      num_io = num_io +  1
      if (io_method .EQ. 0) then
          call MPI_File_read_all(fp, u, 1, combined_btype, mstatus, err)
      else
          call MPI_File_read(fp, u, 1, combined_btype, mstatus, err)
      endif
      if (err .ne. MPI_SUCCESS) then
          call MPI_Error_string(err, err_string, err_len, err)
          print *, 'Error: MPI read from file: ',trim(err_string)
          stop
      endif

      end subroutine mpiio_read

      !----< mpiio_cleanup >-------------------------------------------
      subroutine mpiio_cleanup
      implicit none

      integer err

      call MPI_Type_free(combined_btype, err)
      call MPI_Type_free(filetype, err)

      call MPI_File_close(fp, err)

      end subroutine mpiio_cleanup

      end module mpiio_m
