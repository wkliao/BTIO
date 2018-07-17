!
!  Copyright (C) 2013, Northwestern University
!  See COPYRIGHT notice in top-level directory.
!

      !----< print_io_hints() >-----------------------------------------
      subroutine print_io_hints(info)
          ! print the MPI info objects to stdout
          use mpi
          implicit none
          integer, intent(in) :: info

          ! local variables
          character*(MPI_MAX_INFO_VAL) key, value
          integer                      i, nkeys, valuelen, err
          logical                      flag

 2011     FORMAT('    ',A32,' = ',A)
          call MPI_Info_get_nkeys(info, nkeys, err)
          print *, '---- MPI file info used ----'
          do i=0, nkeys-1
              key(:) = ' '
              call MPI_Info_get_nthkey(info, i, key, err)
              call MPI_Info_get(info, key, MPI_MAX_INFO_VAL, value, flag, err)
              call MPI_Info_get_valuelen(info, key, valuelen, flag, err)
              value(valuelen+1:) = ' '
              if (key(len_trim(key):len_trim(key)) .EQ. char(0)) &
                  key(len_trim(key):) = ' '
              print 2011, trim(key), trim(value)
          enddo
          print *
      end subroutine print_io_hints

      !----< get_file_striping() >-------------------------------------
      subroutine get_file_striping(info, striping_factor, striping_unit)
          ! get the file striping information from the MPI info objects
          use mpi
          implicit none
          integer, intent(in)  :: info
          integer, intent(out) :: striping_factor
          integer, intent(out) :: striping_unit

          ! local variables
          character*(MPI_MAX_INFO_VAL) key, value
          integer                      i, nkeys, valuelen, err
          logical                      flag

          call MPI_Info_get_nkeys(info, nkeys, err)
          do i=0, nkeys-1
              key(:) = ' '
              call MPI_Info_get_nthkey(info, i, key, err)
              call MPI_Info_get(info, key, MPI_MAX_INFO_VAL, value, flag, err)
              call MPI_Info_get_valuelen(info, key, valuelen, flag, err)
              value(valuelen+1:) = ' '
              if (key(len_trim(key):len_trim(key)) .EQ. char(0)) &
                  key(len_trim(key):) = ' '
              if (trim(key) .EQ. 'striping_factor') &
                  read(value, '(i10)') striping_factor
              if (trim(key) .EQ. 'striping_unit') &
                  read(value, '(i10)') striping_unit
          enddo
      end subroutine get_file_striping

