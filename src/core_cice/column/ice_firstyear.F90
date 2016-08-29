!  SVN:$Id: ice_firstyear.F90 1012 2015-06-26 12:34:09Z eclare $
!=======================================================================
!
! First year concentration tracer for sea ice
!
! see 
! Armour, K. C., C. M. Bitz, L. Thompson and E. C. Hunke (2011). Controls
! on Arctic sea ice from first-year and multi-year ice survivability.
! J. Climate, 24, 23782390. doi: 10.1175/2010JCLI3823.1.
!
! authors C. Bitz, University of Washington, modified from ice_age module
!
! 2012: E. Hunke adopted from CESM into CICE, changed name from ice_FY.F90
!
      module ice_firstyear

      use ice_kinds_mod
      use ice_constants_colpkg, only: secday, c0

      implicit none

      private
      public :: update_FYarea

!=======================================================================

      contains

!=======================================================================

!  Zero ice FY tracer on fixed day of year. Zeroing FY ice tracer promotes
!  ice to MY ice. Unfortunately some frazil ice may grow before the 
!  zeroing date and thus get promoted to MY ice too soon.
!  Bummer.

      subroutine update_FYarea (dt,                 &
                                nhmask,   shmask,   &
                                yday,     FYarea)

      real (kind=dbl_kind), intent(in) :: &
         dt , &                ! time step
         yday                  ! day of the year

      logical (kind=log_kind), &
         intent(in) :: &
         nhmask, shmask

      real (kind=dbl_kind), &
         intent(inout) :: &
         FYarea

      if ((yday >= 259._dbl_kind) .and. &
          (yday <  259._dbl_kind+dt/secday)) then
         if (nhmask) FYarea = c0
      endif

      if ((yday >= 75._dbl_kind) .and. &
          (yday <  75._dbl_kind+dt/secday)) then
         if (shmask) FYarea = c0
      endif

      end subroutine update_FYarea

!=======================================================================

      end module ice_firstyear

!=======================================================================
