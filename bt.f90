!
!  Copyright (C) 2013, Northwestern University
!  See COPYRIGHT notice in top-level directory.
!

      !----< main >----------------------------------------------------
      program main
      use mpi
      use header
      use mpiio_m
      use pnetcdf_m
      implicit none

      character io_mode
      character(LEN=128) filename, cmd
      integer i, err, argc, iargc, fstatus, io_method
      integer*8 n3
      integer striping_factor, striping_unit
      double precision navg, t, tmax

      call MPI_Init(err)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

      ! check if number of processes is a square number
      ncells = dint(dsqrt(dble(nprocs) + 0.00001d0))

      if (nprocs .NE. ncells*ncells) then
         print*, 'Number of processes must be a square number. exit ...'
         goto 999
      endif

      root = 0

      !---------------------------------------------------------------------
      !      Root reads input file (if it exists) else takes
      !      defaults from parameters
      !---------------------------------------------------------------------
      if (rank .eq. root) then
         ! take filename from command-line argument if there is any
         call getarg(0, cmd)
         argc = IARGC()
         if (argc .GT. 1) then
            print*,'Usage: ',trim(cmd),' [filename]'
            niter = -1
            goto 777
         endif
         filename = 'inputbt.data'
         if (argc .EQ. 1) call getarg(1, filename)

         open(unit=2,file=trim(filename),status='old', iostat=fstatus)

         if (fstatus .eq. 0) then
 233        FORMAT(' Reading from input file ',A)
            write(*,233) trim(filename)
            read (2,*) io_mode
            read (2,*) io_method
            read (2,*) niter
            read (2,*) grid_points(1), grid_points(2), grid_points(3)
            read (2,'(a)') dir_path
            close(2)
         else
 234        format(' No input file inputbt.data. Exiting ...')
            write(*,234) 
            niter = -1
            goto 777
         endif
      endif

 777  call MPI_Bcast(io_mode,     1, MPI_CHARACTER, root, MPI_COMM_WORLD, err)
      call MPI_Bcast(io_method,   1, MPI_INTEGER,   root, MPI_COMM_WORLD, err)
      call MPI_Bcast(niter,       1, MPI_INTEGER,   root, MPI_COMM_WORLD, err)
      call MPI_Bcast(grid_points, 3, MPI_INTEGER8,  root, MPI_COMM_WORLD, err)
      call MPI_Bcast(dir_path,  128, MPI_CHARACTER, root, MPI_COMM_WORLD, err)

      if (niter .EQ. -1) goto 999

      call allocate_variables

      call make_set

      if (io_method .LT. 2) then ! 0: collective I/O, 1: independent I/O
         err = mpiio_setup(io_mode)
      else
         err = pnetcdf_setup(io_mode, io_method)
      endif
      if (err .EQ. 0) goto 999

      !---------------------------------------------------------------------
      !      Synchronize before placing time stamp
      !---------------------------------------------------------------------
      call MPI_Barrier(MPI_COMM_WORLD, err)
      t = MPI_Wtime()
      num_io = 0

      do i=1, niter
         if (io_mode .EQ. 'w') then
            if (io_method .LT. 2) then ! 0: collective I/O, 1: independent I/O
               call mpiio_write(io_method)
            else
               call pnetcdf_write
            endif
         else
            if (io_method .LT. 2) then ! 0: collective I/O, 1: independent I/O
               call mpiio_read(io_method)
            else
               call pnetcdf_read
            endif
         endif
      end do

      if (io_method .LT. 2) then ! 0: collective I/O, 1: independent I/O
         call mpiio_cleanup
      else
         call pnetcdf_cleanup
      endif

      call deallocate_variables

      t = MPI_Wtime() - t
      call MPI_Reduce(t, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX,  &
                      root, MPI_COMM_WORLD, err)

      if ( rank .eq. root ) then
         striping_factor = 0
         striping_unit   = 0
         call get_file_striping(info_used, striping_factor, striping_unit)
         ! call print_io_hints(info_used)

         n3 = grid_points(1)*grid_points(2)*grid_points(3)

         navg = n3 * 5 * 8       ! I/O amount per write/read
         navg = navg * num_io    ! I/O amount for all write/read
         navg = navg / 1048576.0 ! I/O amount in MB

 2000    FORMAT(A)
 2001    FORMAT(A, A)
 2002    FORMAT(A, I9)
 2003    FORMAT(A, I9, A)
 2004    FORMAT(A, F12.2)
 2005    FORMAT(A, F12.2, A)

         if (io_mode .EQ. 'w') then
            print 2000,'-- BT-IO Benchmark (write operation only) --'
         else
            print 2000,'-- BT-IO Benchmark (read  operation only) --'
         endif
         print 2002,'Number of MPI processes  : ',nprocs
         print 2002,'Global array size X      : ',grid_points(1)
         print 2002,'Global array size Y      : ',grid_points(2)
         print 2002,'Global array size Z      : ',grid_points(3)
         print 2002,'Number of I/O iterations : ',niter
         print 2005,'Totail I/O amount        : ',navg, ' MiB'
         print 2004,'Time in sec              : ',tmax

         navg = navg / tmax      ! I/O Bandwidth in MB/s
         print 2005,'I/O bandwidth            : ',navg, ' MiB/s'

         print 2000, '------------------------------------------'
         if (io_method .EQ. 0) then ! 0: collective I/O, 1: independent I/O
            print 2000,'Using MPI collective I/O method'
         elseif (io_method .EQ. 1) then
            print 2000,'Using MPI independent I/O method'
         elseif (io_method .EQ. 2) then
            print 2000,'Using Parallel netCDF blocking I/O method'
         else
            print 2000,'Using Parallel netCDF non-blocking I/O method'
         endif
         print 2001,'output file path         : ',trim(dir_path)
         print 2002,'file striping count      : ',striping_factor
         print 2003,'file striping size       : ',striping_unit,' bytes'

      endif
      if (info_used .NE. MPI_INFO_NULL) &
         call MPI_Info_free(info_used, err)

 999  call MPI_Finalize(err)

      end program main

