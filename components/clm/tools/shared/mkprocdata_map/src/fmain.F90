program fmain

  use mkprocdata_map, only : mkmap
  implicit none

  character(len= 256) :: arg
  integer :: n                   !index 
  integer :: nargs               !number of arguments  
  integer, external :: iargc     !number of arguments function
  character(len=256) :: filei    !input file
  character(len=256) :: fileo    !output mapped file
  character(len=256) :: fmap     !maping file
  character(len=256) :: ftemplate !template file, containing lat & lon arrays desired in output file
  character(len=256) :: cmdline  !input command line
  integer, parameter :: inival = -999 !initial value for command-line integers
  !----------------------------------------------------

  filei = ' '
  fileo = ' '
  fmap  = ' '
  ftemplate = ' '

  cmdline = 'mkprocdata_map'
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
       filei = trim(arg)
       cmdline = trim(cmdline) // ' -i ' // trim(arg)
    case ('-o')
      call getarg (n, arg)
      n = n + 1
      fileo = trim(arg)
      cmdline = trim(cmdline) // ' -o ' // trim(arg)
    case ('-m')      
      call getarg (n, arg)
      n = n + 1
      fmap = trim(arg)
      cmdline = trim(cmdline) // ' -m ' // trim(arg)
    case ('-t')      
      call getarg (n, arg)
      n = n + 1
      ftemplate = trim(arg)
      cmdline = trim(cmdline) // ' -t ' // trim(arg)
    case default
       write (6,*) 'Argument ', arg,' is not known'
       call usage_exit (' ')
       cmdline = trim(cmdline) // ' ' // trim(arg)
    end select
  end do
  
  if (filei == ' ' .or. fileo == ' ' .or. fmap == ' ' &
      .or. ftemplate == ' ') then
    call usage_exit ('Must specify all the following arguments')
  end if

  call mkmap (filei, fileo, fmap, ftemplate)

end program fmain


subroutine usage_exit (arg)
  implicit none
  character(len=*) :: arg
  if (arg /= ' ') write (6,*) arg
  write (6,*) 'Usage: mkprocdata_map -i <input file>  -o <output file> -m <map file> -t <template file>'
  write (6,*)
  write (6,*) "The template file must contain the dimensions 'lat' and 'lon';"
  write (6,*) "these are used to determine the number of latitudes and longitudes in the output file"
  stop 1 
end subroutine
