!
!  Copyright (C) 2013, Northwestern University
!  See COPYRIGHT notice in top-level directory.
!
!  $Id: make_set.f90 2176 2013-11-06 22:13:59Z wkliao $


      !----< make_set >------------------------------------------------
      subroutine make_set
      !----------------------------------------------------------------
      ! This function allocates space for a set of cells and fills the set
      ! such that communication between cells on different ranks is only
      ! nearest neighbor
      !----------------------------------------------------------------
      use mpi
      use header
      implicit none

      integer i, j, c, dir, err
      integer(KIND=8) size, excess

      ! ncells is the square root of the total number of processes
      
      !----------------------------------------------------------------
      ! determine the location of the cell at the bottom of the 3D 
      ! array of cells
      !----------------------------------------------------------------
      ! rank is the process rank assigned by MPI
      ! cell_coord is the Cardesian id
      ! 3D array is partitioned only along XY. Axis Z is not partitioned

      cell_coord(1,1) = mod(rank,ncells) 
      cell_coord(2,1) = rank/ncells 
      cell_coord(3,1) = 0

      !----------------------------------------------------------------
      ! set the cell_coords for cells in the rest of the z-layers; 
      ! this comes down to a simple linear numbering in the z-direct-
      ! ion, and to the doubly-cyclic numbering in the other dirs     
      !----------------------------------------------------------------
      do c = 2, ncells
         cell_coord(1,c) = mod(cell_coord(1,c-1)+1,ncells) 
         cell_coord(2,c) = mod(cell_coord(2,c-1)-1+ncells,ncells) 
         cell_coord(3,c) = c-1
      end do

      !----------------------------------------------------------------
      ! offset all the coordinates by 1 to adjust for Fortran arrays
      !----------------------------------------------------------------
      do dir = 1, 3
         do c = 1, ncells
            cell_coord(dir,c) = cell_coord(dir,c) + 1
         end do
      end do

      !----------------------------------------------------------------
      ! fill the predecessor and successor entries, using the indices 
      ! of the bottom cells (they are the same at each level of k 
      ! anyway) acting as if full periodicity pertains; note that ncells is
      ! added to those arguments to the mod functions that might
      ! otherwise return wrong values when using the modulo function
      !----------------------------------------------------------------
      i = cell_coord(1,1)-1
      j = cell_coord(2,1)-1

      !----------------------------------------------------------------
      ! now compute the sizes of the cells
      !----------------------------------------------------------------
      do dir= 1, 3
      !----------------------------------------------------------------
      ! set cell_coord range for each direction
      !----------------------------------------------------------------
      ! grid_points(*) = problem_size (global array size)

         size   = grid_points(dir)/ncells
         ! excess = mod(grid_points(dir),ncells)
         excess = grid_points(dir) - ((grid_points(dir) / ncells) * ncells)
         do c=1, ncells
            if (cell_coord(dir,c) .le. excess) then
               cell_size(dir,c) = size+1
               cell_low(dir,c) = (cell_coord(dir,c)-1)*(size+1)
               cell_high(dir,c) = cell_low(dir,c)+size
            else 
               cell_size(dir,c) = size
               cell_low(dir,c)  = excess*(size+1)+ &
                                  (cell_coord(dir,c)-excess-1)*size
               cell_high(dir,c) = cell_low(dir,c)+size-1
            endif
            if (cell_size(dir, c) .le. 2) then
               write(*,50)
 50            FORMAT(' Error: Cell size too small. Min size is 3')
               call MPI_Abort(MPI_COMM_WORLD, -1, err)
               stop
            endif
         end do
      end do

      end subroutine make_set

