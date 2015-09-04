program fmain
!------------------------------------------------------------------
! Purpose:
!   Usage: interpaerosols [-d] [-r rgridnl] [-v] <infile.nc> <outfile.nc>
!
!   Horizontally interpolate file containing aerosol masses (infile.nc)
!   to CAM grid (as specified by outfile.nc).  Perl script REGRID.pl
!   contained in this directory can be used to create a template
!   outfile.nc.
!
!   Sum the values in the vertical such that output arrays contain the
!   column sum from the bottom of the atmosphere to that level.
!   Convert monthly averages to mid month values.
!   Input file assumed to be MATCH ncep runs, averaged by month.
!        and backsolved to provide Mid-month values)
!     
!  Method:
!    read data from file
!    interpolate data onto CAM horizontal grid
!
!------------------------------------------------------------------
   use globals
   use netcdf
   use error_messages, only : handle_ncerr  

   implicit none

!
! Local workspace
!
   character(len=80) :: arg = ' '        ! cmd line arg
   character(len=80) :: rgridnl = ' '    ! reduced grid namelist (if applicable)
   character(len=256) :: infile = ' '    ! input file name
   character(len=256) :: outfile = ' '   ! output file name
   character(len=256) :: cmdline = ' '   ! command line

   integer :: n                          ! argument counter
   integer :: nargs                      ! number of command line arguments
   integer :: ncprec = nf90_float          ! default precision to write data
   integer :: old_mode                   ! returned from nf90_set_fill
   integer :: ret                        ! return code

   logical :: verbose = .false.          ! verbose output
   logical :: isncol = .false.
   integer iargc
   external iargc
   integer :: nxi,nyi,nxo,nyo,nz,ntime,nyof
!
! Default settings before parsing argument list
!
   nargs = iargc()
   n = 1
   cmdline = 'interpaerosols '
   do while (n <= nargs)
      arg = ' '
      call getarg (n, arg)
      n = n + 1
      
      select case (arg)
      case ('-d')
         ncprec = nf90_double
         cmdline = trim(cmdline) // ' -d'
      case ('-r')
         call getarg (n, arg)
         n = n + 1
         rgridnl = arg
         cmdline = trim(cmdline) // ' -r ' // trim(rgridnl)
      case ('-v')
         verbose = .true.
         cmdline = trim(cmdline) // ' -v'
      case default
         if (infile == ' ') then
            infile = arg
         else if (outfile == ' ') then
            outfile = arg
         else
            write (6,*) 'Argument ', arg,' is not known'
            call usage_exit (' ')
         end if
         cmdline = trim(cmdline) // ' ' // trim(arg)
      end select
   end do
  
   if (infile == ' ' .or. outfile == ' ') then
      call usage_exit ('Must enter an input file and an output file')
   end if
!
! Open input and output netcdf files
!
   call handle_ncerr( nf90_open (infile, NF90_NOWRITE, ncidi),&
     'fmain.F90:90')
   call handle_ncerr( nf90_open (outfile, NF90_WRITE, ncido),&
     'fmain.F90:92')
   ret = nf90_set_fill (ncido, NF90_NOFILL, old_mode)
   if (ret /= NF90_NOERR) then
      write(6,*)'Error calling nf90_set_fill:'
      write(6,*)nf90_strerror (ret)
      stop 999
   end if
   ! start in define mode on output file
   call handle_ncerr( nf90_redef (ncido),     &
     'fmain.F90:100')
!
! Get input and output dimensions for x, y, time
! Input file first
!
   call handle_ncerr( nf90_inq_dimid (ncidi, 'lon', londimidi),&
     'fmain.F90:106')
   call handle_ncerr( nf90_inquire_dimension (ncidi, londimidi, len=nxi),&
     'fmain.F90:108')

   call handle_ncerr( nf90_inq_dimid (ncidi, 'lat', latdimidi),&
     'fmain.F90:111')
   call handle_ncerr( nf90_inquire_dimension (ncidi, latdimidi, len=nyi),&
     'fmain.F90:113')
   
   call handle_ncerr( nf90_inq_dimid (ncidi, 'lev', levdimidi),&
     'fmain.F90:116')
   ! nz is the same for input and output
   call handle_ncerr( nf90_inquire_dimension (ncidi, levdimidi, len=nz), &
     'fmain.F90:118')

   call handle_ncerr( nf90_inq_dimid (ncidi, 'time', timedimidi),&
     'fmain.F90:121')
   call handle_ncerr( nf90_inquire_dimension (ncidi, timedimidi, len=ntime),&
     'fmain.F90:123')
!
! Ensure ntime is 12.  Reason is because mean-preserving code for the moment assumes it.
!
   if (ntime /= 12) then
      write(6,*) 'Size of input time dimension must be 12'
      stop 999
   end if
!
! Now output file
!
   
   ret = nf90_inq_dimid (ncido, 'lon', londimido)
   if(ret==nf90_EBADDIM ) then
      call handle_ncerr( nf90_inq_dimid(ncido, 'ncol',ncoldimido),&
        'fmain.F90:139')
      call handle_ncerr( nf90_inquire_dimension (ncido, ncoldimido, len=nxo),&
        'fmain.F90:141')
      nyo=nxo
      nyof=1
      isncol=.true.
   else if(ret==NF90_NOERR) then
     call handle_ncerr( nf90_inquire_dimension (ncido, londimido, len=nxo),&
       'fmain.F90:147')
     call handle_ncerr( nf90_inq_dimid (ncido, 'lat', latdimido),&
       'fmain.F90:149')
     call handle_ncerr( nf90_inquire_dimension (ncido, latdimido, len=nyo),&
       'fmain.F90:151')
     nyof=nyo
  else
    call handle_ncerr( ret,'fmain.F90:156')
   end if

!
! Assume z dimensions on output file don't exist (code will crash if they do)
! This way it will be guaranteed to match input file.
!
   call handle_ncerr( nf90_def_dim (ncido, 'lev', nz, levdimido),&
     'fmain.F90:160')
   call handle_ncerr( nf90_def_dim (ncido, 'ilev', nz+1, ilevdimido),&
     'fmain.F90:162')
!
! If a time dimension is not found on the output file, create it
!
   if (nf90_inq_dimid (ncido, 'time', timedimido) /= nf90_noerr) then
      call handle_ncerr( nf90_def_dim (ncido, 'time', nf90_unlimited, timedimido),&
        'fmain.F90:168')
   end if
!
! Add global attributes
!
   call addglobal (ncido, cmdline)
!
! Call driver code to read the input data and do the interpolations and/or copies
!
   if(isncol) then
      nyof = 1
   else
      nyof = nyo
   end if
   call driver (ncprec, rgridnl,isncol,nxi,nyi,nxo,nyo,nz,ntime,nxo,nyof)
!
! Close netcdf files.  This is crucial for output file in particular in order
! to ensure that all data get written to the file.
!
   call handle_ncerr( nf90_close (ncidi),&
     'fmain.F90:188')
   call handle_ncerr( nf90_close (ncido),&
     'fmain.F90:190')

   write(6,*)'Successfully interpolated aerosol data to file ', trim (outfile)

   stop 0
end program fmain

subroutine usage_exit (arg)
   implicit none
   character*(*) arg

   if (arg /= ' ') write (6,*) arg
   write (6,*) 'Usage: interpaerosols [-d] [-r rgridnl] [-v] infile.nc outfile.nc'
   write (6,*) '  -d: write output file in double precision'
   write (6,*) '  -r rgridnl: define output on reduced grid defined by nlon in namelist'
   write (6,*) '  -v: verbose printout'
   stop 999
end subroutine usage_exit

#ifdef UNICOSMP
subroutine getarg (n, arg)
  implicit none
  integer, intent(in) :: n
  character(len=*), intent(out) :: arg

  integer :: ilen
  integer :: ierr

  call pxfgetarg (n, arg, ilen, ierr)
  if (ierr /= 0) then
    write(6,*)'getarg: ierr not 0:', ierr
    stop 999
  end if
  return
end subroutine getarg

integer function iargc ()
  integer, external :: ipxfargc

  iargc = ipxfargc ()
  return
end function iargc
#endif
