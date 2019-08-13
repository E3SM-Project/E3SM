! $Id: ice_calendar.F90 56 2007-03-15 14:42:35Z dbailey $
!=======================================================================
!BOP
!
! !MODULE: ice_calendar - calendar routines for managing time
!
! !DESCRIPTION:
!
! Calendar routines for managing time
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! 2006 ECH: Removed 'w' option for history; added 'h' and histfreq_n.
!           Converted to free form source (F90).
!
! !INTERFACE:
!
      module ice_calendar
!
! !USES:
!
      use ice_constants
      use ice_domain_size, only: max_nstrm
      use ice_exit, only: abort_ice
      use ice_fileunits
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: &
         days_per_year        , & ! number of days in one year
         daymo(12)            , & ! number of days in each month
         daycal(13)               ! day number at end of month

      ! 360-day year data
      integer (kind=int_kind) :: &
         daymo360(12)         , & ! number of days in each month
         daycal360(13)            ! day number at end of month
      data daymo360 /   30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30/
      data daycal360/ 0,30, 60, 90,120,150,180,210,240,270,300,330,360/

      ! 365-day year data
      integer (kind=int_kind) :: &
         daymo365(12)         , & ! number of days in each month
         daycal365(13)            ! day number at end of month
      data daymo365 /   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal365/ 0,31, 59, 90,120,151,181,212,243,273,304,334,365/

      ! 366-day year data (leap year)
      integer (kind=int_kind) :: &
         daymo366(12)         , & ! number of days in each month
         daycal366(13)            ! day number at end of month
      data daymo366 /   31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal366/ 0,31, 60, 91,121,152,182,213,244,274,305,335,366/

      integer (kind=int_kind) :: &
         istep    , & ! local step counter for time loop
         istep0   , & ! counter, number of steps taken in previous run
         istep1   , & ! counter, number of steps at current timestep
         mday     , & ! day of the month
         hour     , & ! hour of the year
         month    , & ! month number, 1 to 12
         monthp   , & ! last month
         year_init, & ! initial year
         nyr      , & ! year number
         idate    , & ! date (yyyymmdd)
         idate0   , & ! initial date (yyyymmdd)
         sec      , & ! elapsed seconds into date
         npt      , & ! total number of time steps (dt)
         stop_now     , & ! if 1, end program execution
         write_restart, & ! if 1, write restart now
         diagfreq     , & ! diagnostic output frequency (10 = once per 10 dt)
         dumpfreq_n   , & ! restart output frequency (10 = once per 10 d,m,y)
         nstreams     , & ! number of history output streams
         histfreq_n(max_nstrm) ! history output frequency 

      real (kind=dbl_kind) :: &
         dt             , & ! thermodynamics timestep (s)
         dt_thm         , & ! thermodynamics timestep (s)
         dt_dyn         , & ! dynamics/transport/ridging timestep (s)
         time           , & ! total elapsed time (s)
         time_forc      , & ! time of last forcing update (s)
         yday           , & ! day of the year
         nextsw_cday    , & ! next day for sw calculation
         tday           , & ! absolute day number
         xndt_dyn       , & ! reduced timestep for dynamics: xndt_dyn=dt/dt_dyn
         dayyr              ! number of days per year

      logical (kind=log_kind) :: &
         new_year       , & ! new year = .true.
         new_month      , & ! new month = .true.
         new_day        , & ! new day = .true.
         new_hour       , & ! new hour = .true.
         write_ic       , & ! write initial condition now
         write_history(max_nstrm) ! write history now

      character (len=1) :: &
         histfreq(max_nstrm) , & ! history output frequency, 'y','m','d','h','1'
         dumpfreq               ! restart frequency, 'y','m','d'

      character (len=char_len) :: calendar_type

      integer :: nleaps = 0    ! The number of leap days *before* the current year
                               ! This is set by ice_comp_mct and used there & here (must be consistent)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_calendar - initialize calendar variables
!
! !INTERFACE:
!
      subroutine init_calendar
!
! !DESCRIPTION:
!
! Initialize calendar variables
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         k                          , &

      istep = 0         ! local timestep number
      time=istep0*dt    ! s
#ifdef CCSMCOUPLED
      time_forc = c0    ! for coupled runs
#endif
      yday=c0           ! absolute day number
      mday=0            ! day of the month
      month=0           ! month
      nyr=0             ! year
      idate=00000101    ! date
      sec=0             ! seconds into date
      istep1 = istep0   ! number of steps at current timestep
                        ! real (dumped) or imagined (use to set calendar)
      stop_now = 0      ! end program execution if stop_now=1
      dt_thm = dt       ! convenience copy of thermodynamic timestep
      dt_dyn = dt/xndt_dyn ! dynamics et al timestep

      dayyr = real(days_per_year, kind=dbl_kind)  ! days_per_year set by ice_init


      ! determine initial date (assumes namelist year_init, istep0 unchanged)     
      sec = mod(time,secday)            ! elapsed seconds into date at
                                        ! end of dt
      tday = (time-sec)/secday + c1     ! absolute day number

      if (calendar_type /= "GREGORIAN") then 	
         nyr = int((tday-c1)/dayyr) + 1    ! year number
      endif

      ! reset the number of leap days: this is necessary to add one one
      ! the year turns from a leap-year to a non-leap year
      nleaps = leap_year_count(nyr+year_init-1)

      ! get the daycal variable, which is dependant on calendar and days_per_year
      call get_daycal(year=nyr+year_init-1,days_per_year_in=days_per_year,&
           daycal_out=daycal,daymo_out=daymo)

      ! subtract the number of days in prior years from the number of days to get the number of days this year
      yday = tday-real(nleaps,kind=dbl_kind)-real(nyr-1,kind=dbl_kind)*dayyr    ! days that have passed this year

      do k = 1, 12
        if (yday > real(daycal(k),kind=dbl_kind)) month = k
      enddo
      mday = int(yday) - daycal(month)  ! day of the month

      idate0 = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 

      end subroutine init_calendar

!=======================================================================
!BOP
!
! !IROUTINE: calendar - computes date at the end of the time step
!
! !INTERFACE:
!
      subroutine calendar(ttime)
!
! !DESCRIPTION:
!
! Determine the date at the end of the time step
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Tony Craig, NCAR
!
! !USES:
      use ice_fileunits
      use ice_communicate, only: my_task, master_task
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         ttime                          ! time variable
!
!EOP
!
      integer (kind=int_kind) :: &
         k, ns                      , &
         nyrp,mdayp,hourp           , & ! previous year, day, hour
         elapsed_days               , & ! since beginning this run
         elapsed_months             , & ! since beginning this run
         elapsed_hours                  ! since beginning this run

      ! nyrp and new_year are Handled in driver for gregorian calendar
      if (calendar_type /= "GREGORIAN") nyrp=nyr
      monthp=month
      mdayp=mday
      hourp=hour
      if (calendar_type /= "GREGORIAN") new_year=.false.
      new_month=.false.
      new_day=.false.
      new_hour=.false.
      write_history(:)=.false.
      write_restart=0

      sec = mod(ttime,secday)           ! elapsed seconds into date at
                                        ! end of dt
      tday = (ttime-sec)/secday + c1    ! absolute day number
      
      if (calendar_type /= "GREGORIAN") then 	
         nyr = int((tday-c1)/dayyr) + 1    ! year number
         if (nyr   /= nyrp)   new_year = .true.
      endif

      ! reset the number of leap days: this is necessary to add one one
      ! the year turns from a leap-year to a non-leap year
      nleaps = leap_year_count(nyr+year_init-1)

      ! get the daycal variable, depending on calendar and days_per_year
      call get_daycal(year=nyr+year_init-1,days_per_year_in=days_per_year,&
           daycal_out=daycal) 

      ! subtract the number of days in prior years from the number of days to get the number of days this year
      yday = tday-real(nleaps,kind=dbl_kind)-real(nyr-1,kind=dbl_kind)*dayyr    ! days that have passed this year

      do k = 1, 12
        if (yday > real(daycal(k),kind=dbl_kind)) month = k
      enddo
      mday = int(yday) - daycal(month)  ! day of the month 

      hour = int((ttime-dt)/c3600) + c1 ! hour  

      elapsed_months = (nyr - 1)*12 + month - 1
      elapsed_days = int(tday) - 1 
      elapsed_hours = int(ttime/3600)

      idate = (nyr+year_init-1)*10000 + month*100 + mday ! date (yyyymmdd) 

#ifndef CCSMCOUPLED
      if (istep >= npt+1)  stop_now = 1
#endif

      if (month /= monthp) new_month = .true.
      if (mday  /= mdayp)  new_day = .true.
      if (hour  /= hourp)  new_hour = .true.


      do ns = 1, nstreams
         if (histfreq(ns)=='1' .and. histfreq_n(ns)/=0) then
             if (mod(istep1, histfreq_n(ns))==0) &
                write_history(ns)=.true.
         endif
      enddo

      if (istep > 1) then

        do ns = 1, nstreams

           select case (histfreq(ns))
           case ("y", "Y")
             if (new_year  .and. histfreq_n(ns)/=0) then
                if (mod(nyr, histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("m", "M")
             if (new_month .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_months,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("d", "D")
             if (new_day  .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_days,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           case ("h", "H")
             if (new_hour  .and. histfreq_n(ns)/=0) then
                if (mod(elapsed_hours,histfreq_n(ns))==0) &
                   write_history(ns) = .true.
             endif
           end select

        enddo ! nstreams

        select case (dumpfreq)
        case ("y", "Y")
          if (new_year  .and. mod(nyr, dumpfreq_n)==0) &
                write_restart = 1
        case ("m", "M")
          if (new_month .and. mod(elapsed_months,dumpfreq_n)==0) &
                write_restart=1
        case ("d", "D")
          if (new_day   .and. mod(elapsed_days, dumpfreq_n)==0) &
                write_restart = 1
        case default
          call abort_ice('ice_calendar: Invalid dumpfreq: '//dumpfreq)
        end select
      endif

      if (my_task == master_task .and. mod(istep,diagfreq) == 0 &
                                 .and. stop_now /= 1) then
        write(nu_diag,*) ' '
        write(nu_diag,'(a7,i10,4x,a6,i10,4x,a4,i10)') &
             'istep1:', istep1, 'idate:', idate, 'sec:', sec
      endif

      end subroutine calendar

!=======================================================================
      subroutine get_daycal(year,days_per_year_in,daycal_out,daymo_out)

      ! Input/output paramters
        integer, intent(in), optional  :: year    ! year
        integer, intent(in), optional  :: days_per_year_in   ! 360 or 365
        integer, intent(out), optional :: daycal_out(13)     ! cumumulative days per month
        integer, intent(out), optional :: daymo_out(12)     ! days per month

      ! 360-day year data
      integer (kind=int_kind) :: &
         daymo360(12)         , & ! number of days in each month
         daycal360(13)            ! day number at end of month
      data daymo360 /   30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30/
      data daycal360/ 0,30, 60, 90,120,150,180,210,240,270,300,330,360/

      ! 365-day year data
      integer (kind=int_kind) :: &
         daymo365(12)         , & ! number of days in each month
         daycal365(13)            ! day number at end of month
      data daymo365 /   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daycal365/ 0,31, 59, 90,120,151,181,212,243,273,304,334,365/

      ! 366-day year data (leap year)
      integer (kind=int_kind) :: &
         daycal366(13)            ! day number at end of month
      data daycal366/ 0,31, 60, 91,121,152,182,213,244,274,305,335,366/


        if ( present(daymo_out) ) then
           if ( .not. present(days_per_year_in)) &
                call abort_ice('ice: get_daycal needs days_per_year_in to return daymo_out')
           if (days_per_year_in.eq.360) then
              daymo_out  = daymo360
           elseif (days_per_year_in.eq.365) then
              daymo_out  = daymo365
           else
              call abort_ice('ice: year must have 360 or 365 days')
           endif
        endif ! daymo_out

        if ( present(daycal_out) ) then

           ! initialize to check for non-setting
           daycal_out(:) = 0

           ! calculate from days_per_year_in
           if ( present(days_per_year_in) ) then
              if (days_per_year_in.eq.360) then
                 daycal_out = daycal360
              elseif (days_per_year_in.eq.365) then
                 daycal_out = daycal365
              else
                 call abort_ice('ice: year must have 360 or 365 days')
              endif
           endif  ! present(days_per_year_in) 


           if (calendar_type == "GREGORIAN") then 	
              if ( .not. present(year) ) &
                   call abort_ice('ice: get_daycal needs year to return daycal_out for Gregorian calendar')
              if ( is_leap_year(year) ) then
                 daycal_out = daycal366
              else
                 daycal_out = daycal365
              endif
           endif ! calendar_type GREGORIAN

           if ( daycal_out(13) .eq. 0 ) call abort_ice('ice: get_daycal failed to set daycal_out')

        endif  ! daycal_out


      end subroutine get_daycal

!=======================================================================

      logical function is_leap_year(year)
        ! returns .true. if year is a leap year

        ! Input/output paramters
        integer, intent(in) :: year

        is_leap_year = .false.
        if (mod(year,  4) == 0) is_leap_year = .true.
        if (mod(year,100) == 0) is_leap_year = .false.
        if (mod(year,400) == 0) is_leap_year = .true.

        end function is_leap_year

!=======================================================================
      integer function leap_year_count(Y)
        ! counts the number of leap years since year 1

        ! Input/output paramters
        integer, intent(in) :: Y


        if (calendar_type == "GREGORIAN") then 	
           ! count the number of leap years before Y
           if ( Y .lt. 0 ) then
              leap_year_count = 0
              write(6,*) 'WARNING: leap_year_count for year ',Y,&
                   'assumes no leap years before year 0'
           else
              leap_year_count  = ( (Y-1)/4 - (Y-1)/100 + (Y-1)/400 ) + 1
           endif
        else
           leap_year_count = 0
        endif

        ! set module variable
        nleaps = leap_year_count

        return

      end function leap_year_count

      end module ice_calendar

!=======================================================================
