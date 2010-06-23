module pio_msg_mod
  use pio_kinds
  use pio_types
! PIO ASYNC MESSAGE TAGS
   integer, parameter, public :: pio_msg_create_file = 300
   integer, parameter, public :: pio_msg_exit = 999




contains

  subroutine pio_msg_handler(iosystem)
    use pio_types
    use mpi ! _EXTERNAL
    implicit none
    type(iosystem_desc_t) :: iosystem
    integer :: msg = 0, ierr

    do while(msg /= pio_msg_exit)
       print *,__FILE__,__LINE__, iosystem%intercomm
       call mpi_bcast(msg, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
       print *,__FILE__,__LINE__,msg
       select case(msg) 
       case (PIO_MSG_CREATE_FILE)
          call create_file_handler(iosystem, ierr)
       case (PIO_MSG_EXIT)
          print *,'Exiting'
       case default
          print *, 'Got unrecognized message ', msg, ierr
       end select   
    end do

    print *,__FILE__,__LINE__
    call mpi_finalize(ierr)
    stop

  end subroutine pio_msg_handler

end module pio_msg_mod

subroutine create_file_handler(iosystem, ierr)
  use mpi !_EXTERNAL
  use pio
  use pio_kinds
  implicit none
  type(iosystem_desc_t) :: iosystem
  integer :: ierr
  integer :: iotype, amode
  
  character(len=char_len) :: fname
  type(file_desc_t) :: file
  
  call mpi_bcast(fname, char_len, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(iotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  
  call mpi_bcast(amode, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  
  ierr= pio_createfile(iosystem, file, iotype, trim(fname), amode, 1)

  print *,__FILE__,__LINE__, file%fh

  call mpi_bcast(file%fh, 1, mpi_integer, iosystem%iomaster, iosystem%intercomm, ierr)


end subroutine create_file_handler


