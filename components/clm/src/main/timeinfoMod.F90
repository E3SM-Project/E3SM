module timeinfoMod

  use shr_kind_mod       , only : r8 => shr_kind_r8
  use ForcingUpdateMod  , only : caldaym
  !variables needed for time steps
  implicit none
  real(r8), parameter :: dtime_mod   = 3600d0
  real(r8), parameter :: dayspyr_mod = 365d0
  integer  :: year_curr = 1 , year_prev = 1
  integer  :: mon_curr  = 1 , mon_prev  = 1
  integer  :: day_curr  = 1 , day_prev  = 1
  integer  :: secs_curr = 0 , secs_prev = 0
  integer  :: nstep_mod = 0
  integer  :: jday_mod  = 1            ! day into year?
  real(r8) :: thiscalday_mod  = 1.0_r8  ! day number including hours
  real(r8) :: nextsw_cday_mod = 1.0_r8 !nextsw_cday = mod((nstep/(86400._r8/dtime))*1.0_r8,365._r8)+1._r8
  logical  :: end_cd_mod = 0           ! end of current day
  logical  :: doalb = .false.
  real(r8) :: declin, declinp1
  logical :: first = .true.
  !$acc declare create(declin   , declinp1)
  !$acc declare copyin(dtime_mod, dayspyr_mod , &
  !$acc year_curr , year_prev ,&
  !$acc mon_curr  , mon_prev  ,&
  !$acc day_curr  , day_prev  ,&
  !$acc secs_curr , secs_prev ,&
  !$acc nstep_mod       ,&
  !$acc jday_mod        ,&
  !$acc thiscalday_mod  ,&
  !$acc nextsw_cday_mod ,&
  !$acc end_cd_mod, doalb,first )                  !1   2   3   4  5   6   7    8   9  10  11  12
  integer, dimension(12) :: days_per_mon = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  !$acc declare copyin(days_per_mon)

  public :: increment_time_vars

contains
  !
  subroutine increment_time_vars()
    !this subroutine is meant to do all time_advancement on gpu
    !$acc routine seq
    implicit none
    integer  :: year
    integer  :: day
    integer  :: secs
    integer  :: mon
    !
    if(first) then
            first = .false.
            return
    end if 
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
    thiscalday_mod  = 1.0_r8 + nstep_mod*(1.0/24.0_r8)
    nextsw_cday_mod = thiscalday_mod
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
