!===============================================================================
! SVN $Id: shr_func_mod.F90 25434 2010-11-04 22:46:24Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/release_tags/cesm1_2_x_n02_share3_130715/shr/shr_func_mod.F90 $
!===============================================================================

MODULE shr_func_mod

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
   public :: shr_func_freezetemp

   character(SHR_KIND_CS) :: tfreeze_option     ! option for computing freezing point
                                    ! default is constant -1.8C
                                    ! 'linear_salt' for linear equation
                                    ! 'mushy' for CICE mushy-layer nonlinear equation

!===============================================================================
CONTAINS
!===============================================================================

real(SHR_KIND_R8) pure FUNCTION shr_func_freezetemp(s)

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

   real   (SHR_KIND_R8) :: salt

   !----------------------------------------------------------------------------

   salt = max(s,0.0_SHR_KIND_R8)

   if (trim(tfreeze_option) == 'linear_salt') then
      shr_func_freezetemp = -0.0544_SHR_KIND_R8*salt
   elseif(trim(tfreeze_option) == 'mushy') then
!     This form is the high temperature part of the liquidus relation (Assur 1958)
      shr_func_freezetemp = salt / (-18.48_SHR_KIND_R8 + (0.01848_SHR_KIND_R8*salt))
   else
      shr_func_freezetemp = -1.8_SHR_KIND_R8
   endif

   shr_func_freezetemp = max(shr_func_freezetemp,-2.0_SHR_KIND_R8)

END FUNCTION shr_func_freezetemp

!===============================================================================

END MODULE shr_func_mod
