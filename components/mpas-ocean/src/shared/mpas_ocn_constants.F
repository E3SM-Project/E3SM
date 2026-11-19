! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_constants
!
!> \brief MPAS ocean specific constants
!> \author Doug Jacobsen
!> \date   04/25/12
!> \details
!>  This module contains constants specific to the ocean model.
!
!-----------------------------------------------------------------------

module ocn_constants

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_kind_types

   use shr_const_mod

   implicit none
   private
   save

   public :: ocn_constants_init

   type (mpas_pool_type), public, pointer :: ocnConfigs
   type (mpas_pool_type), public, pointer :: ocnPackages

   character (len=*), public, parameter :: statsAlarmID = 'stats'
   character (len=*), public, parameter :: coupleAlarmID = 'coupling'

   ! physical constants, all from E3SM's shared constants (shr_const_mod),
   ! which take them from the PhysicalConstantsDictionary (PCD)

   real (kind=RKIND), public ::      &
      rho_fw                ,&! density of fresh water (kg/m^3)
      rho_sw                ,&! density of salt water (kg/m^3)
      cp_sw                 ,&! specific heat salt water (J/kg/K)
      cp_ice                ,&! specific heat fresh ice (J/kg/K)
      vonkar                ,&! von Karman constant
      latent_heat_fusion_mks,&! lat heat of fusion (J/kg)
      ocn_ref_salinity      ,&! ocean reference salinity (psu)
      ice_ref_salinity      ,&! sea ice reference salinity (psu)
      atm_ref_pressure      ,&! standard sea level pressure (Pa)
      molecular_weight_C    ,&! molecular weight carbon
      molecular_weight_Fe   ,&! molecular weight iron
      T0_Kelvin             ,&! zero point for Celsius
      pi                      ! pi

   !  conversion factors

   real (kind=RKIND), public :: &
      seconds_per_day  ,&! seconds per day
      days_per_second  ,&! days per second
      ppt_to_salt      ,&! salt ppt to kg/kg
      hflux_factor     ,&! heat flux (W/m^2) to temp flux (C*m/s)
      sflux_factor     ,&! salt flux (kg/m^2/s) to salt flux (psu*m/s)
      mmol_to_kg_C       ! mmol-C to kg-C

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_constants_init
!
!> \brief   Initializes the ocean constants
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details
!>  This routine sets up constants for use in the ocean model.
!
!-----------------------------------------------------------------------
   subroutine ocn_constants_init(configPool, packagePool)!{{{
       type (mpas_pool_type), pointer :: configPool
       type (mpas_pool_type), pointer :: packagePool

       ocnConfigs => configPool
       ocnPackages => packagePool

       !-----------------------------------------------------------------------
       !
       !  physical constants
       !
       !-----------------------------------------------------------------------

       T0_Kelvin = real(SHR_CONST_TKFRZ, RKIND)   ! zero point for Celsius
       pi        = real(SHR_CONST_PI, RKIND)      ! pi
       rho_sw    = real(SHR_CONST_RHOSW, RKIND)   ! density of salt water (kg/m^3)
       rho_fw    = real(SHR_CONST_RHOFW, RKIND)   ! density of fresh water (kg/m^3)
       cp_sw     = real(SHR_CONST_CPSW, RKIND)    ! specific heat salt water (J/kg/K)
       cp_ice    = real(SHR_CONST_CPICE, RKIND)   ! specific heat fresh ice (J/kg/K)
       vonkar    = real(SHR_CONST_KARMAN, RKIND)  ! von Karman constant
       latent_heat_fusion_mks = real(SHR_CONST_LATICE, RKIND)    ! lat heat of fusion (J/kg)
       ocn_ref_salinity       = real(SHR_CONST_OCN_REF_SAL, RKIND) ! (psu)
       ice_ref_salinity       = real(SHR_CONST_ICE_REF_SAL, RKIND) ! (psu)
       atm_ref_pressure       = real(SHR_CONST_PSTD, RKIND)      ! (Pa)
       molecular_weight_C     = real(SHR_CONST_MWC, RKIND)       ! molecular weight carbon
       molecular_weight_Fe    = real(SHR_CONST_MWFE, RKIND)      ! molecular weight iron

       !-----------------------------------------------------------------------
       !
       !  conversion factors
       !
       !-----------------------------------------------------------------------

       seconds_per_day = real(SHR_CONST_CDAY, RKIND) ! seconds per day
       days_per_second = 1._RKIND/seconds_per_day    ! days per second
       ppt_to_salt   = 1.e-3_RKIND        ! salt ppt to kg/kg

       ! convert heat flux (W/m^2) to temperature flux (C*m/s):
       !    temp_flux = heat_flux / (rho_sw*cp_sw)
       hflux_factor = 1.0_RKIND/(rho_sw*cp_sw)

       ! convert salt flux (kg/m^2/s) to salt flux (psu*m/s)
       sflux_factor = 1.0_RKIND

       mmol_to_kg_C  = molecular_weight_C  * 1.e-6_RKIND

   end subroutine ocn_constants_init!}}}

!***********************************************************************

end module ocn_constants

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
