module timeinfoMod

  use shr_kind_mod       , only : r8 => shr_kind_r8
  !variables needed for time steps
  implicit none
  real(r8) :: dtime_mod   = 3600.d0
  real(r8) :: dayspyr_mod = 365.d0
  logical :: first = .true.
  integer  :: year_curr = 1 , year_prev = 1
  integer  :: mon_curr  = 1 , mon_prev  = 1
  integer  :: day_curr  = 1 , day_prev  = 1
  integer  :: secs_curr = 0 , secs_prev = 0
  integer  :: nstep_mod = 0
  integer, dimension(12) :: days_per_mon = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer  :: jday_mod  = 1            ! day into year?
  real(r8) :: thiscalday_mod  = 1.0_r8  ! day number including hours
  real(r8) :: nextsw_cday_mod = 1.0_r8 !nextsw_cday = mod((nstep/(86400._r8/dtime))*1.0_r8,365._r8)+1._r8
  logical  :: end_cd_mod = .false.     ! end of current day
  logical  :: doalb = .false.
  !$acc declare copyin(dtime_mod,dayspyr_mod, year_curr, year_prev, &
  !$acc mon_curr, mon_prev, day_curr, day_prev, secs_curr, secs_prev, nstep_mod, &
  !$acc days_per_mon, jday_mod, thiscalday_mod, nextsw_cday_mod, end_cd_mod, doalb )
contains

  subroutine increment_time_vars()
    !this subroutine is meant to do all time_advancement on gpu
    !$acc routine seq
    implicit none
    integer  :: year
    integer  :: day
    integer  :: secs
    integer  :: mon
    !
    end_cd_mod = .false.
    nstep_mod = nstep_mod + 1
    year = year_curr;  mon = mon_curr;
    day  = day_curr ;  secs = secs_curr;
    !
    secs = mod(secs + int(dtime_mod), 86400)
    if(secs == 0) then !new day
      day = day + 1;
      end_cd_mod = .true.
      if ( day > days_per_mon(mon) ) then !next month no leap year
        jday_mod = jday_mod + 1 !julian day?
        day = 1
        mon = mon + 1
        if( mon > 12 ) then
          mon = 1
          year = year + 1
        end if
      end if
    end if
    !
    if(nstep_mod>1) nextsw_cday_mod = mod((nstep_mod/(86400._r8/dtime_mod))*1.0_r8,365._r8)+1._r8
    thiscalday_mod  = 1.0_r8 + nstep_mod*(1.0/24.0_r8)
    !
    !
    year_prev = year_curr; mon_prev = mon_curr;
    day_prev = day_curr  ; secs_prev = secs_curr

    year_curr = year; mon_curr = mon;
    day_curr  = day ; secs_curr = secs ;

    if (nstep_mod == 0) then
          doalb = .false.
    else if (nstep_mod == 1) then
          doalb = .false. !(abs(nextsw_cday_mod- caldayp1) < 1.e-10_r8)
    else
          doalb = (nextsw_cday_mod >= -0.5_r8)
    end if

  end subroutine



end MODULE
