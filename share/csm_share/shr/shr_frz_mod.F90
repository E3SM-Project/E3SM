!===============================================================================
! SVN $Id: shr_frz_mod.F90 25434 2010-11-04 22:46:24Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/release_tags/cesm1_2_x_n02_share3_130715/shr/shr_frz_mod.F90 $
!===============================================================================

MODULE shr_frz_mod

!===============================================================================
!  This is a general module that can be used for global functions.
!===============================================================================
   use shr_kind_mod, only: SHR_KIND_R8, SHR_KIND_IN, SHR_KIND_CS
   use shr_sys_mod
   use shr_const_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   IMPLICIT none

   !----------------------------------------------------------------------------
   ! PUBLIC: Interfaces and global data
   !----------------------------------------------------------------------------
   public :: shr_frz_freezetemp, shr_frz_freezetemp_init

   interface shr_frz_freezetemp
      module procedure shr_frz_freezetemp_0d
      module procedure shr_frz_freezetemp_1d
      module procedure shr_frz_freezetemp_2d
   end interface

   integer, public, parameter :: TFREEZE_OPTION_MINUS1P8 = 1
   integer, public, parameter :: TFREEZE_OPTION_LINEAR_SALT = 2
   integer, public, parameter :: TFREEZE_OPTION_MUSHY = 3

   character(SHR_KIND_CS), public :: tfreeze_option     ! option for computing freezing point
                                         ! minus1p8 is constant -1.8C
                                         ! linear_salt is linear equation
                                         ! mushy for CICE mushy-layer nonlinear equation

   private

   integer :: tfrz_option

!===============================================================================
CONTAINS
!===============================================================================

subroutine shr_frz_freezetemp_init

       !---------------------------------------------------------------
       ! Check tfreeze_option
       !---------------------------------------------------------------
       if (trim(tfreeze_option) == 'minus1p8') then
          write(s_logunit,*) ' tfreeze_option is minus1p8'
          tfrz_option = TFREEZE_OPTION_MINUS1P8
       elseif (trim(tfreeze_option) == 'linear_salt') then
          write(s_logunit,*) ' tfreeze_option is linear_salt'
          tfrz_option = TFREEZE_OPTION_LINEAR_SALT
       elseif (trim(tfreeze_option) == 'mushy') then
          write(s_logunit,*) ' tfreeze_option is mushy'
          tfrz_option = TFREEZE_OPTION_MUSHY
       else
          call shr_sys_abort('shr_frz_freezetemp_init: not a valid tfreeze_option')
       endif

end subroutine shr_frz_freezetemp_init

FUNCTION shr_frz_freezetemp_0d(s)

   !----------------------------------------------------------------------------
   !
   ! FUNCTION to return the freezing point of salt water
   !
   !--------------- Code History -----------------------------------------------
   !
   ! Original Author: David Bailey
   ! Date:            Feb, 2016
   !----------------------------------------------------------------------------

   real   (SHR_KIND_R8),intent(in) :: s ! Salinity in psu
   real   (SHR_KIND_R8) :: shr_frz_freezetemp_0d

   !----------------------------------------------------------------------------

   if (tfrz_option == TFREEZE_OPTION_MINUS1P8) then
      shr_frz_freezetemp_0d = -1.8_SHR_KIND_R8
   elseif (tfrz_option == TFREEZE_OPTION_LINEAR_SALT) then
      shr_frz_freezetemp_0d = -0.0544_SHR_KIND_R8*max(s,0.0_SHR_KIND_R8)
   elseif (tfrz_option == TFREEZE_OPTION_MUSHY) then
!     This form is the high temperature part of the liquidus relation (Assur 1958)
      shr_frz_freezetemp_0d = max(s,0.0_SHR_KIND_R8) &
         / (-18.48_SHR_KIND_R8 + (0.01848_SHR_KIND_R8*max(s,0.0_SHR_KIND_R8)))
   else
   endif

   shr_frz_freezetemp_0d = max(shr_frz_freezetemp_0d,-2.0_SHR_KIND_R8)

END FUNCTION shr_frz_freezetemp_0d

FUNCTION shr_frz_freezetemp_1d(s)

   !----------------------------------------------------------------------------
   !
   ! FUNCTION to return the freezing point of salt water
   !
   !--------------- Code History -----------------------------------------------
   !
   ! Original Author: David Bailey
   ! Date:            Feb, 2016
   !----------------------------------------------------------------------------

   real   (SHR_KIND_R8),intent(in) :: s(:) ! Salinity in psu
   real   (SHR_KIND_R8) :: shr_frz_freezetemp_1d(size(s))

   !----------------------------------------------------------------------------

   if (tfrz_option == TFREEZE_OPTION_MINUS1P8) then
      shr_frz_freezetemp_1d = -1.8_SHR_KIND_R8
   elseif (tfrz_option == TFREEZE_OPTION_LINEAR_SALT) then
      shr_frz_freezetemp_1d = -0.0544_SHR_KIND_R8*max(s,0.0_SHR_KIND_R8*s)
   elseif (tfrz_option == TFREEZE_OPTION_MUSHY) then
!     This form is the high temperature part of the liquidus relation (Assur 1958)
      shr_frz_freezetemp_1d = max(s,0.0_SHR_KIND_R8*s) &
         / (-18.48_SHR_KIND_R8 + (0.01848_SHR_KIND_R8*max(s,0.0_SHR_KIND_R8*s)))
   else
   endif

   shr_frz_freezetemp_1d = max(shr_frz_freezetemp_1d,-2.0_SHR_KIND_R8)

END FUNCTION shr_frz_freezetemp_1d

FUNCTION shr_frz_freezetemp_2d(s)

   !----------------------------------------------------------------------------
   !
   ! FUNCTION to return the freezing point of salt water
   !
   !--------------- Code History -----------------------------------------------
   !
   ! Original Author: David Bailey
   ! Date:            Feb, 2016
   !----------------------------------------------------------------------------

   real   (SHR_KIND_R8),intent(in) :: s(:,:) ! Salinity in psu
   real   (SHR_KIND_R8) :: shr_frz_freezetemp_2d(size(s,1),size(s,2))

   !----------------------------------------------------------------------------

   if (tfrz_option == TFREEZE_OPTION_MINUS1P8) then
      shr_frz_freezetemp_2d = -1.8_SHR_KIND_R8
   elseif (tfrz_option == TFREEZE_OPTION_LINEAR_SALT) then
      shr_frz_freezetemp_2d = -0.0544_SHR_KIND_R8*max(s,0.0_SHR_KIND_R8)
   elseif (tfrz_option == TFREEZE_OPTION_MUSHY) then
!     This form is the high temperature part of the liquidus relation (Assur 1958)
      shr_frz_freezetemp_2d = max(s,0.0_SHR_KIND_R8) &
         / (-18.48_SHR_KIND_R8 + (0.01848_SHR_KIND_R8*max(s,0.0_SHR_KIND_R8)))
   else
   endif

   shr_frz_freezetemp_2d = max(shr_frz_freezetemp_2d,-2.0_SHR_KIND_R8)

END FUNCTION shr_frz_freezetemp_2d
!===============================================================================

END MODULE shr_frz_mod
