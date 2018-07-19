!
!  Copyright (C) 2013, Northwestern University
!  See COPYRIGHT notice in top-level directory.
!

      module header
      implicit none

      integer num_io, ncells
      integer(KIND=8) grid_points(3) ! global array size

      integer, allocatable :: cell_coord (:,:), &
                              cell_low   (:,:), &
                              cell_high  (:,:), &
                              cell_size  (:,:)

      ! cell_coord (3,ncells), cell_low (3,ncells), &
      ! cell_high  (3,ncells), cell_size(3,ncells), &
      ! start      (3,ncells), end      (3,ncells)

      integer(KIND=8) IMAX, JMAX, KMAX

      double precision, allocatable :: u(:, :, :, :, :)
      ! u(5, -2:IMAX+1, -2:JMAX+1, -2:KMAX+1, ncells)

      integer rank, nprocs, root, niter, info_used

      character(len=128) dir_path

      contains

      !----< allocate_variables >--------------------------------------
      subroutine allocate_variables
          use mpi
          implicit none

          info_used = MPI_INFO_NULL

          IMAX = (grid_points(1)/ncells) + 1
          JMAX = (grid_points(2)/ncells) + 1
          KMAX = (grid_points(3)/ncells) + 1

          allocate(cell_coord(3,ncells))
          allocate(  cell_low(3,ncells))
          allocate( cell_high(3,ncells))
          allocate( cell_size(3,ncells))

          allocate(u(5, -2:IMAX+1, -2:JMAX+1, -2:KMAX+1, ncells))

          ! initialize contents of variable u
          call RANDOM_NUMBER(u(1:5,-2:IMAX+1,-2:JMAX+1,-2:KMAX+1,1:ncells))

      end subroutine allocate_variables

      !----< deallocate_variables >------------------------------------
      subroutine deallocate_variables
          implicit none

          deallocate(u)
          deallocate(cell_coord)
          deallocate(  cell_low)
          deallocate( cell_high)
          deallocate( cell_size)
      end subroutine deallocate_variables

      end module header


