#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module time_mod
  !------------------
  use kinds, only : real_kind
  !------------------
  implicit none
  integer,public                :: nsplit=1
  integer,public                :: nmax          ! Max number of timesteps
  integer,public                :: nEndStep      ! Number of End Step
  integer,public                :: ndays         ! Max number of days
  real (kind=real_kind), public :: tstep         ! Dynamics timestep
  real (kind=real_kind), public :: tevolve=0.    ! time evolved since start of dynamics (end of physics) - phl
  real (kind=real_kind), public :: phys_tscale=0.! Physics time scale

  real (kind=real_kind), public, parameter :: secphr = 3600.0D0 ! Timestep filter
  real (kind=real_kind), public, parameter :: secpday = 86400.0D0 ! Timestep filter

  ! smooth now in namelist
  real (kind=real_kind), public :: smooth  = 0.05D0    ! Timestep filter
  integer, parameter :: ptimelevels = 3                           ! number of time levels in the dycore

  type, public :: TimeLevel_t
     integer nm1      ! relative time level n-1
     integer n0       ! relative time level n
     integer np1      ! relative time level n+1
     integer nstep    ! time level since simulation start
     integer nstep0   ! timelevel of first complete leapfrog timestep
  end type TimeLevel_t

  ! Methods
  public :: Time_at
  public :: TimeLevel_update
  public :: TimeLevel_init
  public :: TimeLevel_Qdp

  interface TimeLevel_init
     module procedure TimeLevel_init_default
     module procedure TimeLevel_init_specific
     module procedure TimeLevel_init_copy
  end interface

contains

  function Time_at(nstep) result(tat)
    integer, intent(in) :: nstep
    real (kind=real_kind) :: tat
    tat = nstep*tstep
  end function Time_at

  subroutine TimeLevel_init_default(tl)
    type (TimeLevel_t), intent(out) :: tl
    tl%nm1   = 1
    tl%n0    = 2
    tl%np1   = 3
    tl%nstep = 0
    tl%nstep0 = 2
  end subroutine TimeLevel_init_default

  subroutine TimeLevel_init_copy(tl, tin)
    type (TimeLevel_t), intent(in) :: tin
    type (TimeLevel_t), intent(out) :: tl
    tl%nm1   = tin%nm1
    tl%n0    = tin%n0
    tl%np1   = tin%np1
    tl%nstep = tin%nstep
    tl%nstep0= tin%nstep0
  end subroutine TimeLevel_init_copy

  subroutine TimeLevel_init_specific(tl,n0,n1,n2,nstep)
    type (TimeLevel_t) :: tl
    integer, intent(in) :: n0,n1,n2,nstep
    tl%nm1= n0
    tl%n0 = n1
    tl%np1= n2
    tl%nstep= nstep
  end subroutine TimeLevel_init_specific


  !this subroutine returns the proper
  !locations for nm1 and n0 for Qdp - because
  !it only has 2 levels for storage
  subroutine TimeLevel_Qdp(tl, qsplit, n0, np1)
    type (TimeLevel_t) :: tl
    integer, intent(in) :: qsplit
    integer, intent(inout) :: n0
    integer, intent(inout), optional :: np1

    integer :: i_temp

    i_temp = tl%nstep/qsplit

    if (mod(i_temp,2)  ==0) then
       n0 = 1
       if (present(np1)) then
          np1 = 2
       endif
    else
       n0 = 2
       if (present(np1)) then 
          np1 = 1
       end if
    endif

    !print * ,'nstep = ', tl%nstep, 'qsplit= ', qsplit, 'i_temp = ', i_temp, 'n0 = ', n0

  end subroutine TimeLevel_Qdp

  subroutine TimeLevel_update(tl,uptype)
    type (TimeLevel_t) :: tl
    character(len=*)   :: uptype

    ! Local Variable

    integer :: ntmp
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    if (uptype == "leapfrog") then
       ntmp    = tl%np1
       tl%np1  = tl%nm1
       tl%nm1  = tl%n0
       tl%n0   = ntmp
    else if (uptype == "forward") then
       ntmp    = tl%np1
       tl%np1  = tl%n0
       tl%n0   = ntmp
    else 
       print *,'WARNING: TimeLevel_update called wint invalid uptype=',uptype
    end if
       
    tl%nstep = tl%nstep+1
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
!$OMP BARRIER    
#endif
  end subroutine TimeLevel_update

end module time_mod
