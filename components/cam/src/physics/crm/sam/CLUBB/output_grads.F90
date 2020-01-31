!-------------------------------------------------------------------------------
! $Id: output_grads.F90 5867 2012-07-03 21:06:44Z dschanen@uwm.edu $
module output_grads


! Description:
!   This module contains structure and subroutine definitions to
!   create GrADS output data files for one dimensional arrays.
!
!   The structure type (stat_file) contains all necessay information
!   to generate a GrADS file and a list of variables to be output
!   in the data file.
!
! References:
!   None
!
! Original Author:
!   Chris Golaz, updated 2/18/2003
!-------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public  :: open_grads, write_grads

  private :: format_date, check_grads, &
    determine_time_inc

  ! Undefined value
  real( kind = core_rknd ), private, parameter :: undef = -9.99e33_core_rknd

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine open_grads( iunit, fdir, fname,  & 
                         ia, iz, z, & 
                         day, month, year, rlat, rlon, & 
                         time, dtwrite, & 
                         nvar, grads_file )
! Description:
!   Opens and initialize variable components for derived type 'grads_file'
!   If the GrADS file already exists, open_grads will overwrite it.

! References:
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only:  & 
        fstderr,  & ! Variable 
        fstdout

    use stat_file_module, only: & 
        stat_file ! Type

    use clubb_precision, only:  & 
        time_precision ! Variable

    implicit none

    ! Input Variables

    integer, intent(in) :: iunit   ! File unit being written to            [-]

    character(len=*), intent(in) ::  & 
      fdir,  & ! Directory where file is stored        [-]
      fname    ! Name of file                          [-]

    integer, intent(in) :: & 
      ia,                    & ! Lower Bound of z      [-]
      iz                       ! Upper Bound of z      [-]

    real( kind = core_rknd ), dimension(:), intent(in) :: z

    integer, intent(in) ::  & 
      day,           & ! Day of Month at Model Start    [dd]
      month,         & ! Month of Year at Model start   [mm]
      year             ! Year at Model Start            [yyyy]

    real( kind = core_rknd ), dimension(1), intent(in) :: &
      rlat, rlon ! Latitude and Longitude [Degrees N/E]

    real(kind=time_precision), intent(in) ::  & 
      time,     & ! Time since Model start          [s]
      dtwrite     ! Time interval for output        [s]

    ! Number of GrADS variables to store            [#]
    integer, intent(in) :: nvar

    ! Input/Output Variables
    type (stat_file), intent(inout) :: &
      grads_file ! File data [-]

    ! Local Variables

    integer :: k
    logical :: l_ctl, l_dat, l_error

    ! ---- Begin Code ----

    ! Define parameters for the GrADS ctl and dat files

    grads_file%iounit = iunit
    grads_file%fdir   = fdir
    grads_file%fname  = fname
    grads_file%ia     = ia
    grads_file%iz     = iz

    ! Determine if the altitudes are ascending or descending and setup the
    ! variable z accordingly.
    if ( ia <= iz ) then
      do k=1,iz-ia+1
        grads_file%z(k) = z(ia+k-1)
      end do
    else
      do k=1,ia-iz+1
        grads_file%z(k) = z(ia-k+1)
      end do
    end if

    grads_file%day   = day
    grads_file%month = month
    grads_file%year  = year

    allocate( grads_file%rlat(1), grads_file%rlon(1) )

    grads_file%rlat  = rlat
    grads_file%rlon  = rlon

    grads_file%dtwrite = dtwrite

    grads_file%nvar = nvar

    ! Check whether GrADS files already exists

    ! We don't use this feature for the single-column model.  The
    ! clubb_standalone program will simply overwrite existing data files if they
    ! exist.  The restart function will create a new GrADS file starting from
    ! the restart time in the output directory.

    ! inquire( file=trim(fdir)//trim(fname)//'.ctl', exist=l_ctl )
    ! inquire( file=trim(fdir)//trim(fname)//'.dat', exist=l_dat )

    l_ctl = .false.
    l_dat = .false.

    ! If none of the files exist, set ntimes and nrecord and
    ! to initial values and return

    if ( .not.l_ctl .and. .not.l_dat ) then

      grads_file%time = time
      grads_file%ntimes = 0
      grads_file%nrecord = 1
      return

      ! If both files exists, attempt to append to existing files

    else if ( l_ctl .and. l_dat ) then

      !  Check existing ctl file

      call check_grads( iunit, fdir, fname,  & 
                        ia, iz, & 
                        day, month, year, time, dtwrite, & 
                        nvar,  & 
                        l_error, grads_file%ntimes, grads_file%nrecord, &
                        grads_file%time )

      if ( l_error ) then
        write(unit=fstderr,fmt=*) "Error in open_grads:"
        write(unit=fstderr,fmt=*)  & 
        "Attempt to append to existing files failed"
!              call stopcode('open_grads')
        stop 'open_grads'
      end if

      return

!         If one file exists, but not the other, give up

    else
      write(unit=fstderr,fmt=*) 'Error in open_grads:'
      write(unit=fstderr,fmt=*)  & 
        "Attempt to append to existing files failed,"//  & 
        " because only one of the two GrADS files was found."
      stop "open_grads"

    end if

    return
  end subroutine open_grads

!-------------------------------------------------------------------------------
  subroutine check_grads( iunit, fdir, fname,  & 
                          ia, iz, & 
                          day, month, year, time, dtwrite, & 
                          nvar,  & 
                          l_error, ntimes, nrecord, time_grads )
! Description:
!   Given a GrADS file that already exists, this subroutine will attempt
!   to determine whether data can be safely appended to existing file.
! References:
!   None
!-------------------------------------------------------------------------------
    use stat_file_module, only: & 
        variable ! Type

    use clubb_precision, only: & 
        time_precision ! Variable

    use constants_clubb, only:  & 
        fstderr,  & ! Variable 
        fstdout,  &
        sec_per_hr, &
        sec_per_min

    implicit none

    ! Input Variables

    integer, intent(in) ::  & 
      iunit,              & ! Fortran file unit
      ia, iz,            & ! First and last level
      day, month, year,  & ! Day, month and year numbers
      nvar              ! Number of variables in the file

    character(len=*), intent(in) :: & 
      fdir, fname ! File directory and name

    real(kind=time_precision), intent(in) :: & 
      time    ! Current model time        [s]

    real(kind=time_precision), intent(in) :: & 
      dtwrite ! Time interval between writes to the file    [s]

    ! Output Variables
    logical, intent(out) ::  & 
      l_error

    integer, intent(out) ::  & 
      ntimes, nrecord

    real(kind=time_precision), intent(out) :: time_grads

    ! Local Variables
    logical :: l_done
    integer :: ierr
    character(len = 256) :: line, tmp, date, dt

    integer ::  & 
      i, nx, ny, nzmax, & 
      ihour, imin, & 
      ia_in, iz_in, ntimes_in, nvar_in, & 
      day_in, month_in, year_in

    real(kind=time_precision) :: dtwrite_in

    real( kind = core_rknd ), dimension(:), allocatable :: z_in

    type (variable), dimension(:), pointer :: var_in

!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    ! Initialize logical variables
    l_error = .false.
    l_done  = .false.

    ! Open control file
    open( unit = iunit, & 
          file = trim( fdir )//trim( fname )//'.ctl', & 
          status = 'old', iostat = ierr )
    if ( ierr < 0 ) l_done = .true.

    ! Read and process it

    read(unit=iunit,iostat=ierr,fmt='(a256)') line
    if ( ierr < 0 ) l_done = .true.

    do while ( .not. l_done )

      if ( index(line,'XDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, nx
        if ( nx /= 1 ) then
          write(unit=fstderr,fmt=*) 'Error: XDEF can only be 1'
          l_error = .true.
        end if

      else if ( index(line,'YDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, ny
        if ( ny /= 1 ) then
          write(unit=fstderr,fmt=*) "Error: YDEF can only be 1"
          l_error = .true.
        end if

      else if ( index(line,'ZDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, iz_in

        if ( index(line,'LEVELS') > 0 ) then
          ia_in = 1
          allocate( z_in(ia_in:iz_in) )
          read(unit=iunit,fmt=*) (z_in(i),i=ia_in,iz_in)
        end if

      else if ( index(line,'TDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, ntimes_in, tmp, date, dt
        read(unit=date(1:2),fmt=*) ihour
        read(unit=date(4:5),fmt=*) imin
        time_grads = real( ihour, kind=time_precision ) * sec_per_hr &
                   + real( imin, kind=time_precision ) * sec_per_min
        read(unit=date(7:8),fmt=*) day_in
        read(unit=date(12:15),fmt=*) year_in

        select case( date(9:11) )
        case( 'JAN' )
          month_in = 1
        case( 'FEB' )
          month_in = 2
        case( 'MAR' )
          month_in = 3
        case( 'APR' )
          month_in = 4
        case( 'MAY' )
          month_in = 5
        case( 'JUN' )
          month_in = 6
        case( 'JUL' )
          month_in = 7
        case( 'AUG' )
          month_in = 8
        case( 'SEP' )
          month_in = 9
        case( 'OCT' )
          month_in = 10
        case( 'NOV' )
          month_in = 11
        case( 'DEC' )
          month_in = 12
        case default
          write(unit=fstderr,fmt=*) "Unknown month: "//date(9:11)
          l_error = .true.
        end select

        read(unit=dt(1:len_trim(dt)-2),fmt=*) dtwrite_in
        dtwrite_in = dtwrite_in * sec_per_min

      else if ( index(line,'ENDVARS') > 0 ) then

        l_done = .true.

      else if ( index(line,'VARS') > 0 ) then

        read(line,*) tmp, nvar_in
        allocate( var_in(nvar_in) )
        do i=1, nvar_in
          read(unit=iunit,iostat=ierr,fmt='(a256)') line
          read(unit=line,fmt=*) var_in(i)%name, nzmax
          if ( nzmax /= iz_in ) then
            write(unit=fstderr,fmt=*)  & 
              "Error reading ", trim( var_in(i)%name )
            l_error = .true.
          end if ! nzmax /= iz_in
        end do ! 1..nvar_in
      end if

      read(unit=iunit,iostat=ierr,fmt='(a256)') line
      if ( ierr < 0 ) l_done = .true.

    end do ! while ( .not. l_done )

    close( unit=iunit )

    ! Perform some error check

    if ( abs(ia_in - iz_in) /= abs(ia - iz) ) then
      write(unit=fstderr,fmt=*) "check_grads: size mismatch"
      l_error = .true.
    end if

    if ( day_in /= day ) then
      write(unit=fstderr,fmt=*) "check_grads: day mismatch"
      l_error = .true.
    end if

    if ( month_in /= month ) then
      write(unit=fstderr,fmt=*) "check_grads: month mismatch"
      l_error = .true.
    end if

    if ( year_in /= year ) then
      write(unit=fstderr,fmt=*) "check_grads: year mismatch"
      l_error = .true.
    end if

    if ( int( time_grads ) + ntimes_in*int( dtwrite_in )  & 
         /= int( time ) ) then
      write(unit=fstderr,fmt=*) "check_grads: time mismatch"
      l_error = .true.
    end if

    if ( int( dtwrite_in ) /= int( dtwrite) ) then
      write(unit=fstderr,fmt=*) 'check_grads: dtwrite mismatch'
      l_error = .true.
    end if

    if ( nvar_in /= nvar ) then
      write(unit=fstderr,fmt=*) 'check_grads: nvar mismatch'
      l_error = .true.
    end if

    if ( l_error ) then
      write(unit=fstderr,fmt=*) "check_grads diagnostic"
      write(unit=fstderr,fmt=*) "ia      = ", ia_in, ia
      write(unit=fstderr,fmt=*) "iz      = ", iz_in, iz
      write(unit=fstderr,fmt=*) "day     = ", day_in, day
      write(unit=fstderr,fmt=*) "month   = ", month_in, month
      write(unit=fstderr,fmt=*) "year    = ", year_in, year
      write(unit=fstderr,fmt=*) "time_grads / time    = ", time_grads, time
      write(unit=fstderr,fmt=*) "dtwrite = ", dtwrite_in, dtwrite
      write(unit=fstderr,fmt=*) "nvar    = ", nvar_in, nvar
    end if

    ! Set ntimes and nrecord to append to existing files

    ntimes  = ntimes_in
    nrecord = ntimes_in * nvar_in * iz_in + 1

    deallocate( z_in )

    ! The purpose of this statement is to avoid a compiler warning
    ! for tmp
    if (tmp =="") then
    end if
    ! Joshua Fasching June 2008

    return
  end subroutine check_grads

!-------------------------------------------------------------------------------
  subroutine write_grads( grads_file )

! Description:
!   Write part of a GrADS file to data (.dat) file update control file (.ctl. 
!   Can be called as many times as necessary
! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only: & 
      fstderr ! Variable(s)

    use model_flags, only: &
      l_byteswap_io ! Variable

    use endian, only: & 
      big_endian, & ! Variable
      little_endian

    use stat_file_module, only: & 
      stat_file ! Type

    use clubb_precision, only:  & 
      time_precision, & ! Variable(s)
      core_rknd

    implicit none

    ! External
    intrinsic :: selected_real_kind

    ! Constant parameters
    integer, parameter :: &
      r4 = selected_real_kind( p=5 ) ! Specify 5 decimal digits of precision

    ! Input Variables
    type (stat_file), intent(inout) :: &
      grads_file ! Contains all information on the files to be written to

    ! Local Variables
    integer ::  & 
      i,     & ! Loop indices
      ios   ! I/O status

    character(len=15) :: date

    integer :: dtwrite_ctl ! Time increment for the ctl file
    character(len=2) :: dtwrite_units ! Units on dtwrite_ctl

    ! ---- Begin Code ----
    ! Check number of variables and write nothing if less than 1

    if ( grads_file%nvar < 1 ) return

#include "recl.inc"

    ! Output data to file
    open( unit=grads_file%iounit, & 
          file=trim( grads_file%fdir )//trim( grads_file%fname )//'.dat', & 
          form='unformatted', access='direct', & 
          recl=F_RECL*abs( grads_file%iz-grads_file%ia+1 ), & 
          status='unknown', iostat=ios )
    if ( ios /= 0 ) then
      write(unit=fstderr,fmt=*)  & 
        "write_grads: error opening binary file"
      write(unit=fstderr,fmt=*) "iostat = ", ios
      stop
    end if

    if ( grads_file%ia <= grads_file%iz ) then
      do i=1,grads_file%nvar
        write(grads_file%iounit,rec=grads_file%nrecord)  & 
          real( grads_file%var(i)%ptr(1,1,grads_file%ia:grads_file%iz), kind=r4)
        grads_file%nrecord = grads_file%nrecord + 1
      end do

    else
      do i=1, grads_file%nvar
        write(grads_file%iounit,rec=grads_file%nrecord) & 
          real( grads_file%var(i)%ptr(1,1,grads_file%ia:grads_file%iz:-1), kind=r4)
        grads_file%nrecord = grads_file%nrecord + 1
      end do

    end if ! grads_file%ia <= grads_file%iz

    close( unit=grads_file%iounit, iostat = ios )

    if ( ios /= 0 ) then
      write(unit=fstderr,fmt=*)  & 
        "write_grads: error closing binary file"
      write(unit=fstderr,fmt=*) "iostat = ", ios
      stop
    end if

    grads_file%ntimes = grads_file%ntimes + 1

    ! Write control file

    open(unit=grads_file%iounit, & 
         file=trim( grads_file%fdir )//trim( grads_file%fname )//'.ctl', & 
         status='unknown', iostat=ios)
    if ( ios > 0 ) then
      write(unit=fstderr,fmt=*)  & 
        "write_grads: error opening control file"
      write(unit=fstderr,fmt=*) "iostat = ", ios
      stop
    end if

    ! Write file header
    if ( ( big_endian .and. .not. l_byteswap_io ) &
      .or. ( little_endian .and. l_byteswap_io ) ) then
      write(unit=grads_file%iounit,fmt='(a)') 'OPTIONS BIG_ENDIAN'

    else
      write(unit=grads_file%iounit,fmt='(a)') 'OPTIONS LITTLE_ENDIAN'

    end if

    write(unit=grads_file%iounit,fmt='(a)') 'DSET ^'//trim( grads_file%fname )//'.dat'
    write(unit=grads_file%iounit,fmt='(a,e11.5)') 'UNDEF ',undef
    write(unit=grads_file%iounit,fmt='(a,f8.3,a)') 'XDEF    1 LINEAR ', grads_file%rlon, ' 1.'
    write(unit=grads_file%iounit,fmt='(a,f8.3,a)') 'YDEF    1 LINEAR ', grads_file%rlat, ' 1.'
    if ( grads_file%ia == grads_file%iz ) then
      write(unit=grads_file%iounit,fmt='(a)') 'ZDEF    1 LEVELS 0.'
    else if ( grads_file%ia < grads_file%iz ) then
      write(unit=grads_file%iounit,fmt='(a,i5,a)')  & 
        'ZDEF', abs(grads_file%iz-grads_file%ia)+1,' LEVELS '
      write(unit=grads_file%iounit,fmt='(6f13.4)')  & 
        (grads_file%z(i-grads_file%ia+1),i=grads_file%ia,grads_file%iz)
    else
      write(unit=grads_file%iounit,fmt='(a,i5,a)')  & 
        'ZDEF',abs(grads_file%iz-grads_file%ia)+1,' LEVELS '
      write(grads_file%iounit,'(6f13.4)') (grads_file%z(grads_file%ia-i+1), &
        i=grads_file%ia,grads_file%iz,-1)
    end if

    call format_date( grads_file%day, grads_file%month, grads_file%year, grads_file%time, & ! In
                      date ) ! Out

    call determine_time_inc( grads_file%dtwrite, & ! In
                             dtwrite_ctl, dtwrite_units ) ! Out

    write(unit=grads_file%iounit,fmt='(a,i6,a,a,i5,a)') 'TDEF    ', & 
      grads_file%ntimes, ' LINEAR ', date, dtwrite_ctl, dtwrite_units

    ! Variables description
    write(unit=grads_file%iounit,fmt='(a,i5)') 'VARS', grads_file%nvar

    do i=1, grads_file%nvar, 1
      write(unit=grads_file%iounit,fmt='(a,i5,a,a)') & 
        grads_file%var(i)%name(1:len_trim(grads_file%var(i)%name)), & 
        abs(grads_file%iz-grads_file%ia)+1,' 99 ', & 
        grads_file%var(i)%description(1:len_trim(grads_file%var(i)%description))
    end do

    write(unit=grads_file%iounit,fmt='(a)') 'ENDVARS'

    close( unit=grads_file%iounit, iostat=ios )
    if ( ios > 0 ) then
      write(unit=fstderr,fmt=*)  & 
        "write_grads: error closing control file"
      write(unit=fstderr,fmt=*) "iostat = ",ios
      stop
    end if

    return
  end subroutine write_grads

!---------------------------------------------------------
  subroutine format_date( day_in, month_in, year_in, time_in, &
                          date )
!
! Description: 
!   This subroutine formats the current time of the model (given in seconds 
!   since the start time) to a date format usable as GrADS output.
! References:
!   None
!---------------------------------------------------------
    use clubb_precision, only:  & 
      time_precision, & ! Variable(s)
      core_rknd

    use calendar, only:  & 
      compute_current_date ! Procedure(s)

    use calendar, only: & 
      month_names ! Variable(s)

    use constants_clubb, only: &
      sec_per_hr, & ! Variable(s)
      min_per_hr

    implicit none

    ! Input Variables
    integer, intent(in) :: & 
      day_in,   & ! Day of the Month at Model Start  [dd]
      month_in, & ! Month of the Year at Model Start [mm]
      year_in     ! Year at Model Start              [yyyy]

    real(kind=time_precision), intent(in) ::  & 
      time_in ! Time since Model Start              [s]

    ! Output Variables
    character(len=15), intent(out) ::  & 
      date ! Current Date in format 'hh:mmZddmmmyyyy'

    ! Local Variables
    integer :: iday, imonth, iyear  ! Day, month, year
    real(kind=time_precision) :: time ! time      [s]

    ! ---- Begin Code ----

    ! Copy input arguments into local variables

    iday   = day_in
    imonth = month_in
    iyear  = year_in
    time   = time_in

    call compute_current_date( day_in, month_in, &  ! In
                               year_in, & ! In
                               time_in, & ! In 
                               iday, imonth, & ! Out
                               iyear, & ! Out
                               time ) ! Out

    date = 'hh:mmZddmmmyyyy'
    write(unit=date(7:8),fmt='(i2.2)') iday
    write(unit=date(9:11),fmt='(a3)') month_names(imonth)
    write(unit=date(12:15),fmt='(i4.4)') iyear
    write(unit=date(1:2),fmt='(i2.2)') floor( time/sec_per_hr )
    write(unit=date(4:5),fmt='(i2.2)')  & 
      int( mod( nint( time ), nint(sec_per_hr) ) / nint(min_per_hr) )

    return
  end subroutine format_date

!-------------------------------------------------------------------------------
  subroutine determine_time_inc( dtwrite_sec, &
                                 dtwrite_ctl, units )
! Description:
!   Determine the units on the time increment, since GrADS only allows a 2 digit
!   time increment.
! References:
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only: &
      sec_per_day, & ! Constants
      sec_per_hr, &
      sec_per_min

    use clubb_precision, only:  & 
      time_precision ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, floor

    ! Input Variables
    real(kind=time_precision), intent(in) :: &
      dtwrite_sec ! Time increment in GrADS [s]

    ! Output Variables
    integer, intent(out) :: &
      dtwrite_ctl ! Time increment in GrADS [units vary]

    character(len=2), intent(out) :: units ! Units on dtwrite_ctl

    ! Local variables
    real(kind=time_precision) :: &
      dtwrite_min, & ! Time increment [minutes]
      dtwrite_hrs, & ! Time increment [hours]
      dtwrite_days   ! Time increment [days]

    ! ---- Begin Code ----

    ! Since GrADs can't handle a time increment of less than a minute we assume
    ! 1 minute output for an output frequency of less than a minute.
    dtwrite_min = real( floor( dtwrite_sec/sec_per_min ), kind=time_precision )
    dtwrite_min = max( 1._time_precision, dtwrite_min )

    if ( dtwrite_min <= 99._time_precision ) then
      dtwrite_ctl = int( dtwrite_min )
      units = 'mn'
    else
      dtwrite_hrs = dtwrite_sec / sec_per_hr
      if ( dtwrite_hrs <= 99._time_precision ) then
        dtwrite_ctl = int( dtwrite_hrs )
        units = 'hr'
      else
        dtwrite_days = dtwrite_sec / sec_per_day
        if ( dtwrite_days <= 99._time_precision ) then
          dtwrite_ctl = int( dtwrite_days )
          units = 'dy'
        else
          stop "Fatal error in determine_time_inc"
        end if ! dwrite_days <= 99.
      end if ! dtwrite_hrs <= 99.
    end if ! dtwrite_min <= 99.

    return
  end subroutine determine_time_inc

end module output_grads
!-------------------------------------------------------------------------------
