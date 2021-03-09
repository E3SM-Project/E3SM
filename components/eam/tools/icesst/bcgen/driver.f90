program pcmdisst
!------------------------------------------------------------------------------------
!
! Purpose: Wrapper program for routines (provided by Karl Taylor of PCMDI) which modify 
!          mid-month values of SST and ice concentration to preserve monthly means upon 
!          linear time interpolation.
!
! Method: Read in a namelist containing the relevant variables.  Check for validity,
!         then call driving subroutine.
!
!------------------------------------------------------------------------------------
   implicit none

   include 'netcdf.inc'

   integer :: mon1 = -1           ! start month of period of interest plus buffer
   integer :: iyr1 = -1           ! start year of period of interest plus buffer
   integer :: monn = -1           ! end month of period of interest plus buffer
   integer :: iyrn = -1           ! end year of period of interest plus buffer
   integer :: mon1rd = -1         ! start month to read from input data
   integer :: iyr1rd = -1         ! start year to read from input data
   integer :: monnrd = -1         ! end month to read from input data
   integer :: iyrnrd = -1         ! end year to read from input data
   integer :: mon1clm = -1        ! start month to use for climatology
   integer :: iyr1clm = -1        ! start year to use for climatology
   integer :: monnclm = -1        ! end month to use for climatology
   integer :: iyrnclm = -1        ! end year to use for climatology
   integer :: nlat = -1           ! number of latitudes (e.g. 64 for typical T42)
   integer :: nlon = -1           ! number of longitudes (e.g. 128 for typical T42)
   integer :: mon1out = -1        ! start month written to output file (default mon1rd)
   integer :: iyr1out = -1        ! start year written to output file (default iyr1rd)
   integer :: monnout = -1        ! end month written to output file (default monnrd)
   integer :: iyrnout = -1        ! end year written to output file (default iyrnrd)
   integer :: nmon                ! total number of months of period of interest plus buffer

   logical :: oldttcalc = .false. ! true => bfb agreement with original code

   character(len=120) :: infil      = ' ' ! input filename.
   character(len=120) :: outfilclim = ' ' ! output climatology file
   character(len=120) :: outfilamip = ' ' ! output AMIP-style file
   character(len=256) :: string           ! temporary character variable

   character(len=256) :: arg                 ! cmd line argument
   character(len=512) :: cmdline             ! input command line
   character(len=19)  :: cur_timestamp
   character(len=1024) :: prev_history = ' ' ! history attribute from input file
   character(len=1024) :: history = ' '      ! history attribute for output files
!
! netcdf info
!
   integer :: ncidin = -1         ! input file handle
   integer :: londimid = -1       ! longitude dimension id
   integer :: latdimid = -1       ! latitude dimension id
!
! Cmd line
!
   integer, external :: iargc
   integer           :: n, nargs
!
! Namelist
!
   namelist /cntlvars/ mon1, iyr1, monn, iyrn, mon1rd, iyr1rd, monnrd, iyrnrd, &
                       mon1clm, iyr1clm, monnclm, iyrnclm, mon1out, &
                       iyr1out, monnout, iyrnout, oldttcalc
!
! Read namelist
!
   read (5,cntlvars)
!
! Check that all required input items were specified in the namelist
!
   call verify_input (mon1, 'mon1', 'start month of period of interest plus buffer')
   call verify_input (iyr1, 'iyr1', 'start year of period of interest plus buffer')
   call verify_input (monn, 'monn', 'end month of period of interest plus buffer')
   call verify_input (iyrn, 'iyrn', 'end year of period of interest plus buffer')
   call verify_input (mon1rd, 'mon1rd', 'start month of input data')
   call verify_input (iyr1rd, 'iyr1rd', 'start year of input data')
   call verify_input (monnrd, 'monnrd', 'end month of input data')
   call verify_input (iyrnrd, 'iyrnrd', 'end year of input data')
   call verify_input (mon1clm, 'mon1clm', 'start month of output climatology')
   call verify_input (iyr1clm, 'iyr1clm', 'start year of output climatology')
   call verify_input (monnclm, 'monnclm', 'end month of output climatology')
   call verify_input (iyrnclm, 'iyrnclm', 'end year of output climatology')
!
! Check that all specified input values are valid
!
   call verify_monthindx (mon1, 'mon1')
   call verify_monthindx (monn, 'monn')
   call verify_monthindx (mon1rd, 'mon1rd')
   call verify_monthindx (monnrd, 'monnrd')
   call verify_monthindx (mon1clm, 'mon1clm')
   call verify_monthindx (monnclm, 'monnclm')
!
! Check that input dates are valid with respect to each other
!
   if ((mon1clm + 12*iyr1clm) < (mon1rd + 12*iyr1rd)) then
      write(6,*)'mon1clm, iyr1clm=', mon1clm, iyr1clm
      write(6,*)'mon1rd,  iyr1rd= ', mon1rd, iyr1rd
      call err_exit ('Dates for data to be read must bracket climatology')
   end if

   if (mon1rd + 12*iyr1rd > monnrd + 12*(iyrnrd-1)) then
      write(6,*)'mon1rd, iyr1rd=', mon1rd, iyr1rd
      write(6,*)'monnrd, iyrnrd=', monnrd, iyrnrd
      call err_exit ('must read at least 1 year of data total')
   end if

   if ((monnclm + 12*iyrnclm) > (monnrd + 12*iyrnrd)) then
      write(6,*)'monnclm, iyrnclm=', monnclm, iyrnclm
      write(6,*)'monnrd,  iyrnrd= ', monnrd, iyrnrd
      call err_exit ('End date for climatology exceeds end date of data to be read')
   end if

   if (iyr1rd > iyrnrd .or. iyr1rd == iyrnrd .and. mon1rd > mon1rd) then
      write(6,*)'mon1rd, iyr1rd=', mon1rd, iyr1rd
      write(6,*)'monnrd, iyrnrd=', monnrd, iyrnrd
      call err_exit ('Start date for input data exceeds end date')
   end if

   if (iyr1 > iyrn .or. iyr1 == iyrn .and. mon1 > monn) then
      write(6,*)'mon1, iyr1=', mon1, iyr1
      write(6,*)'monn, iyrn=', monn, iyrn
      call err_exit ('Start date for output data exceeds end date')
   end if

   if (.not. (monn == mon1-1 .or. mon1 == 1 .and. monn == 12)) then
      write(6,*)'mon1, monn=', mon1, monn
      call err_exit ('Period to be treated must be an integral number of years')
   end if

   if ((mon1 + 12*iyr1) > (mon1rd + 12*iyr1rd)) then
      write(6,*)'mon1,   iyr1=  ', mon1, iyr1
      write(6,*)'mon1rd, iyr1rd=', mon1rd, iyr1rd
      call err_exit ('Start date of data to be read must not preceed period of interest')
   end if

   if ((monn + 12*iyrn) < (monnrd + 12*iyrnrd)) then
      write(6,*)'monn,   iyrn=  ', monn, iyrn
      write(6,*)'monnrd, iyrnrd=', monnrd, iyrnrd
      call err_exit ('End date of data to be read must not follow period of interest')
   end if

   if (mod((monn+12-mon1+1),12) /= 0) then
      write(6,*)'mon1, monn=',mon1, monn
      string = 'error in time specifications: only integral number of yrs allowed'
      call err_exit (string)
   end if
!
! Ensure that climatology period is as expected (1982-2001).  To enable other
! averaging periods, just comment out the following bit of code
!
   if (mon1clm /= 1  .or. iyr1clm /= 1982 .or. &
       monnclm /= 12 .or. iyrnclm /= 2001) then
      write(6,*)'Climatological averaging period is not as expected (1982-2001)'
      write(6,*)'If you REALLY want to change the averaging period, delete the'
      write(6,*)'appropriate err_exit call from driver.f90'
      call err_exit ('pcmdisst')
   end if
!
! Calculate derived variables
!
   nmon = (iyrn - iyr1)*12 + (monn - mon1 + 1)

   if (monn - mon1 + 12*(iyrn - iyr1)+1 /= nmon) then
      string = 'error in time specifications: parameter nmon must be consistent with '// &
               'first and last months specified. Check nmon, mon1, iyr1, monn, iyrn'
      call err_exit (string)
   end if

   ! parse command line arguments, saving them to be written to history attribute
   nargs = iargc ()
   n = 1
   cmdline = 'bcgen '
   do while (n <= nargs)
      arg = ' '
      call getarg (n, arg)
      n = n + 1
      select case (arg)
      case ('-i')
         call getarg (n, arg)
         n = n + 1
         infil = arg
         cmdline = trim(cmdline) // ' -i ' // trim(infil)
      case ('-c')
         call getarg (n, arg)
         n = n + 1
         outfilclim = arg
         cmdline = trim(cmdline) // ' -c ' // trim(outfilclim)
      case ('-t')
         call getarg (n, arg)
         n = n + 1
         outfilamip = arg
         cmdline = trim(cmdline) // ' -t ' // trim(outfilamip)
      case default
         write (6,*) 'Argument ', arg,' is not known'
         call usage_exit (' ')
      end select
   end do

!
! Check required input character variables
!
   if (infil == ' ') then
      call usage_exit ('input file must be specified using -i command-line option')
   end if

   if (outfilclim == ' ') then
      call usage_exit ('output climatology file must be specified using -c command-line option')
   end if
   
   if (outfilamip == ' ') then
      call usage_exit ('output time series file must be specified using -t command-line option')
   end if
!
! Set defaults for unspecified variables   
!
   if (mon1out == -1) mon1out = mon1rd
   if (iyr1out == -1) iyr1out = iyr1rd
   if (monnout == -1) monnout = monnrd
   if (iyrnout == -1) iyrnout = iyrnrd
!
! Check validity of output dates
!
   call verify_monthindx (mon1out, 'mon1out')
   call verify_monthindx (monnout, 'monnout')

   if (iyr1out > iyrnout .or. iyr1out == iyrnout .and. mon1out > mon1out) then
      call err_exit ('Start date for output exceeds end date')
   end if
!
! Calculate derived variables
!
   nmon = (iyrn - iyr1)*12 + (monn - mon1 + 1)
!
! Open input file and obtain grid size (nlon and nlat)
!
   call wrap_nf_open (infil, NF_NOWRITE, ncidin)
   call wrap_nf_inq_dimid (ncidin, 'lon', londimid)
   call wrap_nf_inq_dimid (ncidin, 'lat', latdimid)
   call wrap_nf_inq_dimlen (ncidin, londimid, nlon)
   call wrap_nf_inq_dimlen (ncidin, latdimid, nlat)
!
! Add to or define history attribute.
!
   if (nf_get_att_text (ncidin, NF_GLOBAL, 'history', prev_history) /= NF_NOERR) then
      write(6,*)'nf_get_att_text() failed for history attribute'
   end if
   
   call get_curr_timestamp(cur_timestamp)
   if (len_trim(prev_history) == 0) then
      history = cur_timestamp // ' ' // trim(cmdline)
   else
      history = trim(prev_history) // char(10) // cur_timestamp // ' ' // trim(cmdline)
   end if

   write(6,*)'Grid size is nlon,nlat=', nlon, nlat
!
! Call Karl's main program or the derivative code.  Either should give the same answers
!
   call bcgen (mon1, iyr1, monn, iyrn, mon1rd,                &
               iyr1rd, monnrd, iyrnrd, mon1clm, iyr1clm,      &
               monnclm, iyrnclm, nlat, nlon, mon1out,         &
               iyr1out, monnout, iyrnout, ncidin, outfilclim, &
               outfilamip, nmon, oldttcalc, history)
   write(6,*)'Done writing output files ', trim(outfilclim), ' and ', trim(outfilamip)
   stop 0
end program pcmdisst

subroutine verify_monthindx (indx, varname)
!------------------------------------------------------------------------------------
!
! Purpose: Check validity of input month index.
!
!------------------------------------------------------------------------------------
   implicit none
!
! Input arguments
!
   integer, intent(in) :: indx
   character(len=*), intent(in) :: varname

   if (indx < 1 .or. indx > 12) then
      write(6,*) varname, '=', indx, ' is an invalid month index'
      stop 999
   end if

   return
end subroutine verify_monthindx

subroutine err_exit (string)
!------------------------------------------------------------------------------------
!
! Purpose: Print an error message and exit
!
!------------------------------------------------------------------------------------
   implicit none
!
! Input arguments
!
   character(len=*), intent(in) :: string

   write(6,*) string
   stop 999
end subroutine err_exit

subroutine verify_input (ivar, ivarname, string)
!------------------------------------------------------------------------------------
!
! Purpose: Check that a required input variable was actually set.  If not, print an
!          error msg and exit.
!
! Method: All namelist input integer variables must have a non-negative value.  Since
!         they are initialized to -1, a check on < 0 effectively checks whether they 
!         were correctly set in the namelist.
!------------------------------------------------------------------------------------
   implicit none
!
! Input arguments
!
   integer, intent(in) :: ivar
   character(len=*), intent(in) :: ivarname
   character(len=*), intent(in) :: string

   if (ivar < 0) then
      write(6,*) ivarname, ' must be set in the namelist to a non-negative value'
      write(6,*) 'verbose description of this variable: ', trim(string)
      stop 999
   end if

   return
end subroutine verify_input

subroutine usage_exit (arg)
   implicit none
   character*(*) arg
   
   if (arg /= ' ') write (6,*) arg
   write (6,*) 'Usage: bcgen -i infile -c climfile -t tsfile < namelist_file'
   write (6,*) '  -i file: input file'
   write (6,*) '  -c file: output climatology file'
   write (6,*) '  -t file: output time series file'
   stop 999
end subroutine usage_exit

subroutine get_curr_timestamp(time)
! return timestamp formatted as "YYYY-MM-DD HH:MM:SS"   
   character(len=19), intent(out) :: time
   integer :: t(8)
   call date_and_time(values=t)
   write(time,'(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') t(1),'-',t(2),'-',t(3),' ',&
                                                         t(5),':',t(6),':',t(7)
end subroutine get_curr_timestamp
