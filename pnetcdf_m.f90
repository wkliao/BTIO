!
!  Copyright (C) 2013, Northwestern University
!  See COPYRIGHT notice in top-level directory.
!

      module pnetcdf_m
      use mpi
      use pnetcdf
      use header
      implicit none

      integer ncid, dimid(5), varid
      character(len=512) nc_filenm
      integer(KIND=MPI_OFFSET_KIND) global_five_dbl, num_dumps
      integer(KIND=MPI_OFFSET_KIND) put_size, get_size
      integer, allocatable :: buftypes(:), reqs(:), sts(:)
      logical doNonBlockingIO
      double precision t_create, t_post_w, t_wait_w
      double precision t_open,   t_post_r, t_wait_r

      private :: check

      contains

      !----< check() >-------------------------------------------------
      subroutine check(err, message)
      implicit none
      integer err
      character(len=*) message

      ! It is a good idea to check returned value for possible error
      write(6,*) trim(message), trim(nfmpi_strerror(err))
      call MPI_Abort(MPI_COMM_WORLD, -1, err)
      end subroutine check

      !----< create_buffer_type >--------------------------------------
      subroutine create_buffer_type
      implicit none

      integer c, elemtype, err
      integer sizes(3), gh_starts(3), subsizes(3)

      ! define MPI derived data type for local buffers
      ! Note dimensions of local cells, can differ slightly between cells

      call MPI_Type_contiguous(5, MPI_DOUBLE_PRECISION, elemtype, err)
      call MPI_Type_commit(elemtype, err)

      ! local array size
      sizes(1) = IMAX+4
      sizes(2) = JMAX+4
      sizes(3) = KMAX+4

      ! (ghost length = 2 on both ends)
      gh_starts(1) = 2
      gh_starts(2) = 2
      gh_starts(3) = 2

      allocate(buftypes(ncells))

      do c = 1, ncells
         subsizes(1) = cell_size(1, c)
         subsizes(2) = cell_size(2, c)
         subsizes(3) = cell_size(3, c)

         ! Create buffer data type for a local cell
         call MPI_Type_create_subarray(3, sizes, subsizes, gh_starts, &
                                       MPI_ORDER_FORTRAN, elemtype, &
                                       buftypes(c), err)
         call MPI_Type_commit(buftypes(c), err)
      enddo
      call MPI_Type_free(elemtype, err)

      end subroutine create_buffer_type

      !----< pnetcdf_setup >-------------------------------------------
      integer function pnetcdf_setup(io_mode, io_method)
      implicit none

      character io_mode   ! 'w' for write and 'r' for read
      integer   io_method ! 2 for blocking I/O and 3 for nonblocking

      ! local variables
      integer omode, info, err
      double precision t

      t_post_w = 0.0
      t_wait_w = 0.0
      t_post_r = 0.0
      t_wait_r = 0.0

      t = MPI_Wtime()

      pnetcdf_setup = 1

      if (io_method .EQ. 2) doNonBlockingIO = .FALSE.
      if (io_method .EQ. 3) doNonBlockingIO = .TRUE.

      call create_buffer_type

      allocate(reqs(ncells))
      allocate( sts(ncells))

      nc_filenm = trim(dir_path)//'/btio.nc'

      ! disable data sieving for write operations
      call MPI_Info_create(info, err)
      call MPI_Info_set(info, "nc_var_align_size", "1", err);
      call MPI_Info_set(info, "nc_in_place_swap", "enable", err);

      if (io_mode .EQ. 'w') then
         omode = IOR(NF_CLOBBER, NF_64BIT_DATA)
         err = nfmpi_create(MPI_COMM_WORLD, nc_filenm, omode, &
                            info, ncid)
         if (err .ne. NF_NOERR) call check(err, 'In nfmpi_create:')

         call pnetcdf_define
      else
         omode = NF_NOWRITE
         err = nfmpi_open(MPI_COMM_WORLD, nc_filenm, omode, info, ncid)
         if (err .ne. NF_NOERR) then
            pnetcdf_setup = 0
            if (err .EQ. NF_ENOENT) then
               write(6,*)'Error: file does not exit'
            else
               write(6,*)'Error (file open): ',trim(nfmpi_strerror(err))
            endif
            return
         endif

         call pnetcdf_inquiry
      endif

      call MPI_Info_free(info, err)

      ! get the info object used by MPI-IO library
      err = nfmpi_get_file_info(ncid, info_used)

      t = MPI_Wtime() - t

      if (io_mode .EQ. 'w') then
          t_create = t
      else
          t_open = t
      endif

      end function pnetcdf_setup

      !----< pnetcdf_define >-------------------------------------------
      subroutine pnetcdf_define
      implicit none

      integer err

      ! define dimensions

      global_five_dbl = 5
      err = nfmpi_def_dim(ncid, 'FIVE_DBL', global_five_dbl, dimid(1))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_def_dim FIVE_DBL')

      err = nfmpi_def_dim(ncid, 'X', grid_points(1), dimid(2))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_def_dim X')

      err = nfmpi_def_dim(ncid, 'Y', grid_points(2), dimid(3))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_def_dim Y')

      err = nfmpi_def_dim(ncid, 'Z', grid_points(3), dimid(4))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_def_dim Z')

      err = nfmpi_def_dim(ncid, 'NUM_DUMPS', NFMPI_UNLIMITED, dimid(5))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_def_dim NUM_DUMPS')

      err = nfmpi_def_var(ncid, 'var', NF_DOUBLE, 5, dimid, varid)
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_def_var var')

      ! exit define mode

      err = nfmpi_enddef(ncid)
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_enddef:')

      end subroutine pnetcdf_define

      !----< pnetcdf_inquiry >-------------------------------------------
      subroutine pnetcdf_inquiry
      implicit none

      integer err

      ! inquire dimension IDs and lengths

      err = nfmpi_inq_dimid(ncid, 'FIVE_DBL', dimid(1))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimid FIVE_DBL')

      err = nfmpi_inq_dimlen(ncid, dimid(1), global_five_dbl)
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimlen X')

      err = nfmpi_inq_dimid(ncid, 'X', dimid(2))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimid X')

      err = nfmpi_inq_dimlen(ncid, dimid(2), grid_points(1))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimlen X')

      err = nfmpi_inq_dimid(ncid, 'Y', dimid(3))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimid Y')

      err = nfmpi_inq_dimlen(ncid, dimid(3), grid_points(2))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimlen Y')

      err = nfmpi_inq_dimid(ncid, 'Z', dimid(4))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimid Z')

      err = nfmpi_inq_dimlen(ncid, dimid(4), grid_points(3))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimlen Z')

      err = nfmpi_inq_dimid(ncid, 'NUM_DUMPS', dimid(5))
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimid NUM_DUMPS')

      err = nfmpi_inq_dimlen(ncid, dimid(5), num_dumps)
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_dimlen NUM_DUMPS')

      ! inquire variable ID

      err = nfmpi_inq_varid(ncid, 'var', varid)
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_inq_varid')

      end subroutine pnetcdf_inquiry

      !---< pnetcdf_write >---------------------------------------------
      subroutine pnetcdf_write
      implicit none

      integer c, err
      integer(KIND=MPI_OFFSET_KIND) starts(6), counts(6), nReqs
      double precision t_start, t_end

      t_start = MPI_Wtime()

      num_io = num_io +  1

      nReqs = 1

      do c = 1, ncells
         starts(1) = 1
         starts(2) = cell_low(1,c) + 1
         starts(3) = cell_low(2,c) + 1
         starts(4) = cell_low(3,c) + 1
         starts(5) = num_io

         counts(1) = 5
         counts(2) = cell_size(1,c)
         counts(3) = cell_size(2,c)
         counts(4) = cell_size(3,c)
         counts(5) = 1

         if (doNonBlockingIO) then
             err = nfmpi_iput_vara(ncid, varid, starts, counts, &
                                   u(:,:,:,:,c), nReqs, buftypes(c), reqs(c))
             if (err .ne. NF_NOERR) call check(err, 'In nfmpi_iput_vara:')
         else
             err = nfmpi_put_vara_all(ncid, varid, starts, counts, &
                                      u(:,:,:,:,c), nReqs, buftypes(c))
             if (err .ne. NF_NOERR) call check(err, 'In nfmpi_put_vara_all:')
         endif
      enddo

      t_end = MPI_Wtime()
      t_post_w = t_post_w + t_end - t_start
      t_start = t_end

      if (doNonBlockingIO) then
          err = nfmpi_wait_all(ncid, ncells, reqs, sts)
          if (err .ne. NF_NOERR) call check(err, 'In nfmpi_wait_all:')
      endif

      do c = 1, ncells
         if (err .ne. NF_NOERR) &
            call check(sts(c), 'In nfmpi_wait_all status error: ')
      enddo

      t_end = MPI_Wtime()
      t_wait_w = t_wait_w + t_end - t_start

      end subroutine pnetcdf_write

      !---< pnetcdf_read >----------------------------------------------
      subroutine pnetcdf_read
      implicit none

      integer c, err
      integer(KIND=MPI_OFFSET_KIND) starts(6), counts(6), nReqs
      double precision t_start, t_end

      t_start = MPI_Wtime()

      num_io = num_io +  1

      nReqs = 1

      do c = 1, ncells
         starts(1) = 1
         starts(2) = cell_low(1,c) + 1
         starts(3) = cell_low(2,c) + 1
         starts(4) = cell_low(3,c) + 1
         starts(5) = num_io

         counts(1) = 5
         counts(2) = cell_size(1,c)
         counts(3) = cell_size(2,c)
         counts(4) = cell_size(3,c)
         counts(5) = 1

         if (doNonBlockingIO) then
            err = nfmpi_iget_vara(ncid, varid, starts, counts, &
                                  u(:,:,:,:,c), nReqs, buftypes(c), reqs(c))
            if (err .ne. NF_NOERR) call check(err, 'In nfmpi_iget_vara:')
         else
            err = nfmpi_get_vara_all(ncid, varid, starts, counts, &
                                     u(:,:,:,:,c), nReqs, buftypes(c))
            if (err .ne. NF_NOERR) call check(err, 'In nfmpi_get_vara_all:')
         endif
      enddo

      t_end = MPI_Wtime()
      t_post_r = t_post_r + t_end - t_start
      t_start = t_end

      if (doNonBlockingIO) then
         err = nfmpi_wait_all(ncid, ncells, reqs, sts)
         if (err .ne. NF_NOERR) call check(err, 'In nfmpi_wait_all:')
      endif

      do c = 1, ncells
         if (err .ne. NF_NOERR) &
            call check(sts(c), 'In nfmpi_wait_all status error: ')
      enddo

      t_end = MPI_Wtime()
      t_wait_r = t_wait_r + t_end - t_start

      end subroutine pnetcdf_read

      !----< pnetcdf_cleanup >-----------------------------------------
      subroutine pnetcdf_cleanup
      implicit none

      integer c, err

      err = nfmpi_inq_put_size(ncid, put_size)
      err = nfmpi_inq_get_size(ncid, get_size)

      err = nfmpi_close(ncid)
      if (err .ne. NF_NOERR) call check(err, 'In nfmpi_close: ')

      do c = 1, ncells
         call MPI_Type_free(buftypes(c), err)
      enddo

      deallocate(buftypes)
      deallocate(reqs)
      deallocate(sts)

      end subroutine pnetcdf_cleanup

      end module pnetcdf_m
