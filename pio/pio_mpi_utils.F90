#define _FILE_ "pio_mpi_utils.F90"
module pio_mpi_utils
  

  use pio_support, only : piodie
  use pio_types, only : PIO_char, PIO_int, PIO_double, PIO_real
   

  implicit none
  private
  
  public :: pio_type_to_mpi_type

contains

  integer function pio_type_to_mpi_type(ptype) result(mtype)
    use pio_support
    implicit none
    include 'mpif.h'            ! _EXTERNAL
    integer, intent(in):: ptype

    select case(ptype)
      case (PIO_double)
         mtype=MPI_REAL8
      case (PIO_real)
         mtype=MPI_REAL4
      case (PIO_int)
         mtype=MPI_INTEGER
      case (PIO_char)
         mtype=MPI_CHARACTER
      case default
         call piodie( _FILE_,__LINE__, &
                      'Could not convert pio type=',ptype,' to an mpi type')
    end select

  end function pio_type_to_mpi_type

end module pio_mpi_utils
