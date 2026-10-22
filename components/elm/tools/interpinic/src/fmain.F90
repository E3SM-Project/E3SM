program fmain

! Based on earlier code by slevis which read a surface and initial file
! at one resolution and a surface file at another resolution to create
! an initial file at the latter resolution.
!
! This version only needs an input initial file and an output initial file.
! Data from the input initial file is mapped into the output initial file.

  use interpinic    , only : override_missing, interp_filei
  use interpinic_mpi_utils, only : mpi_initialize, mpi_cleanup, masterproc
  implicit none
  include 'netcdf.inc'

  character(len= 256) :: arg
  integer :: n                   !index 
  integer :: nargs               !number of arguments  
  integer :: iargc     !number of arguments function
  character(len=256) :: finidati !input initial dataset to read
  character(len=256) :: finidato !output initial dataset to create
  character(len=256) :: cmdline  !input command line
  !----------------------------------------------------

  ! Every rank parses argv -- they all receive it -- but only rank 0 prints.
  call mpi_initialize()

  finidati = ' '
  finidato = ' '

  cmdline = 'interpinic '
  nargs = iargc()
  n = 1
  do while (n <= nargs)
    arg = ' '
    call getarg (n, arg)
    n = n + 1

    select case (arg)
    case ('-i')
       call getarg (n, arg)
       n = n + 1
       finidati = trim(arg)
       cmdline = trim(cmdline) // ' -i ' // trim(arg)
    case ('-o')
      call getarg (n, arg)
      n = n + 1
      finidato = trim(arg)
      cmdline = trim(cmdline) // ' -o ' // trim(arg)
    case ('-a')
      override_missing = .false.
      cmdline = trim(cmdline) // ' -a '
    case default
       if (masterproc) write (6,*) 'Argument ', arg,' is not known'
       call usage_exit (' ')
       cmdline = trim(cmdline) // ' ' // trim(arg)
    end select
  end do
  
  if (finidati == ' ' .or. finidato == ' ') then
    call usage_exit ('Must specify all the following arguments')
  end if

  ! before calling interp_filei, check with the user that overwriting the output
  ! file is acceptable!!!

  call interp_filei (finidati, finidato, cmdline)

  call mpi_cleanup()

end program fmain


subroutine usage_exit (arg)
  ! Every rank parses the same argv, so every rank reaches here together and
  ! MPI_Finalize stays collective.  Only rank 0 prints.
  use interpinic_mpi_utils, only : masterproc, mpi_cleanup
  implicit none
  character(len=*) :: arg
  if (masterproc) then
     if (arg /= ' ') write (6,*) arg
     write (6,*) 'Usage: interpinic -i <input initial data file>  -o <output initial data file> [options]'
     write (6,*) 'options: -a = abort rather than override missing values with closest bare-soil'
     write (6,*) 'Note - the output initial data file will be overwritten with the interpolated values'
  end if
  call mpi_cleanup()
  stop 999
end subroutine
