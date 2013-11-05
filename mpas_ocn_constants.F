! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_constants
!
!> \brief MPAS ocean specific constants
!> \author Doug Jacobsen
!> \date   04/25/12
!> \version SVN:$Id:$
!> \details
!>  This module contains constants specific to the ocean model.
!
!-----------------------------------------------------------------------

module ocn_constants

   use mpas_kind_types

#ifdef MPAS_CESM
   use shr_kind_mod, only: SHR_KIND_R8
   use shr_const_mod
#endif

   implicit none
   save

   real (kind=RKIND), public ::      &
      rho_air               ,&! ambient air density (kg/m^3)
      rho_fw                ,&! density of fresh water (kg/m^3)
      rho_sw                ,&! density of salt water (kg/m^3)
      rho_ice               ,&! density of sea ice (kg/m^3)
      cp_sw                 ,&! specific heat salt water (J/kg/K)
      cp_air                ,&! heat capacity of air (J/kg/K)
      sound                 ,&! speed of sound (m/s)
      vonkar                ,&! von Karman constant
      emissivity            ,&!
      stefan_boltzmann      ,&! W/m^2/K^4
      latent_heat_vapor     ,&! lat heat of vaporization (erg/g)
      latent_heat_vapor_mks ,&! lat heat of vaporization (J/kg)
      latent_heat_fusion    ,&! lat heat of fusion (erg/g)
      latent_heat_fusion_mks,&! lat heat of fusion (J/kg)
      sea_ice_salinity      ,&! salinity of sea ice formed (psu)
      ocn_ref_salinity        ! ocean reference salinity (psu)

   !  conversion factors

   real (kind=RKIND), public :: &
      T0_Kelvin        ,&! zero point for Celsius
      mpercm           ,&! meters per m
      cmperm           ,&! m per meter
      salt_to_ppt      ,&! salt (kg/kg) to ppt
      ppt_to_salt      ,&! salt ppt to kg/kg
      mass_to_Sv       ,&! mass flux to Sverdrups
      heat_to_PW       ,&! heat flux to Petawatts
      salt_to_Svppt    ,&! salt flux to Sv*ppt
      salt_to_mmday    ,&! salt to water (mm/day)
      hflux_factor     ,&! heat flux (W/m^2) to temp flux (C*m/s)
      fwflux_factor    ,&! fw flux (kg/m^2/s) to salt((msu/psu)*m/s)
      salinity_factor  ,&! fw flux (kg/m^2/s) to salt flux (msu*m/s)
      sflux_factor     ,&! salt flux (kg/m^2/s) to salt flux (msu*m/s)
      fwmass_to_fwflux   ! fw flux (kg/m^2/s) to fw flux (m/s)

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_constants_init
!
!> \brief   Initializes the ocean constants
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \version SVN:$Id$
!> \details 
!>  This routine sets up constants for use in the ocean model.
!
!-----------------------------------------------------------------------

   subroutine ocn_constants_init()!{{{
       integer :: n

       !-----------------------------------------------------------------------
       !
       !  physical constants
       !  some of these constants may be over-ridden by CSM-defined
       !  constants if the CSM shared constants are available
       !
       !-----------------------------------------------------------------------

       T0_Kelvin = 273.16_RKIND               ! zero point for Celsius
       rho_air   = 1.2_RKIND                  ! ambient air density (kg/m^3)
       rho_sw    = 1.026e3_RKIND              ! density of salt water (kg/m^3)
       rho_fw    = 1.0e3_RKIND                ! avg. water density (kg/m^3)
       rho_ice   = 0.917e3_RKIND              ! density of ice (kg/m^3)
       cp_sw     = 3.996e3_RKIND              ! specific heat salt water
       cp_air    = 1005.0_RKIND               ! heat capacity of air (J/kg/K)
       sound     = 1.5e2_RKIND                ! speed of sound (m/s)
       vonkar    = 0.4_RKIND                  ! von Karman constant
       emissivity         = 1.0_RKIND         !
       stefan_boltzmann   = 567.0e-10_RKIND   !  W/m^2/K^4
       latent_heat_vapor_mks  = 2.5e6_RKIND   ! lat heat of vaporization (J/kg)
       latent_heat_fusion = 3.34e6_RKIND      ! lat heat of fusion (erg/kg)
       latent_heat_fusion_mks = 3.337e5_RKIND ! lat heat of fusion (J/kg)
       sea_ice_salinity   =  4.0_RKIND        ! (psu)
       ocn_ref_salinity   = 34.7_RKIND        ! (psu)


       !-----------------------------------------------------------------------
       !
       !  conversion factors
       !
       !-----------------------------------------------------------------------

       salt_to_ppt   = 1000._RKIND        ! salt (kg/kg) to ppt
       ppt_to_salt   = 1.e-3_RKIND        ! salt ppt to kg/kg
       mass_to_Sv    = 1.0e-12_RKIND      ! mass flux to Sverdrups
       heat_to_PW    = 4.186e-15_RKIND    ! heat flux to Petawatts
       salt_to_Svppt = 1.0e-9_RKIND       ! salt flux to Sv*ppt
       salt_to_mmday = 3.1536e+5_RKIND    ! salt to water (mm/day)

       !-----------------------------------------------------------------------
       !
       !  physical constants -- CCSM override
       !
       !-----------------------------------------------------------------------

#ifdef MPAS_CESM
       T0_Kelvin              = SHR_CONST_TKFRZ         ! zero point for Celsius
       cp_sw                  = SHR_CONST_CPSW          ! erg/g/K
       cp_air                 = SHR_CONST_CPDAIR        ! J/kg/K
       rho_air                = SHR_CONST_RHODAIR       ! kg/m^3
       rho_sw                 = SHR_CONST_RHOSW         ! kg/m^3
       rho_fw                 = SHR_CONST_RHOFW         ! kg/m^3
       rho_ice                = SHR_CONST_RHOICE        ! kg/m^3
       vonkar                 = SHR_CONST_KARMAN
       stefan_boltzmann       = SHR_CONST_STEBOL        ! W/m^2/K^4
       latent_heat_vapor_mks  = SHR_CONST_LATVAP        ! J/kg
       latent_heat_fusion_mks = SHR_CONST_LATICE        ! J/kg 
       ocn_ref_salinity       = SHR_CONST_OCN_REF_SAL   ! psu
       sea_ice_salinity       = SHR_CONST_ICE_REF_SAL   ! psu
#endif

!#ifdef ZERO_SEA_ICE_REF_SAL
!      sea_ice_salinity       = c0
!#endif

       !-----------------------------------------------------------------------
       !
       !  convert heat, solar flux (W/m^2) to temperature flux (C*cm/s):
       !  --------------------------------------------------------------
       !    heat_flux in (W/m^2) = (J/s/m^2) = 1(kg/s^3)
       !    density of seawater rho_sw in (kg/m^3)
       !    specific heat of seawater cp_sw in (erg/kg/C) = (m^2/s^2/C)
       !
       !    temp_flux          = heat_flux / (rho_sw*cp_sw)
       !    temp_flux (C*cm/s) = heat_flux (W/m^2)
       !                         * 1 (kg/s^3)/(W/m^2)
       !                         / [(rho_sw*cp_sw) (kg/m/s^2/C)]
       !
       !                       = heat_flux (W/m^2)
       !                         * hflux_factor (C*m/s)/(W/m^2)
       !
       !    ==>  hflux_factor = 1/(rho_sw*cp_sw)
       !
       !-----------------------------------------------------------------------

       hflux_factor = 1.0_RKIND/(rho_sw*cp_sw)

       !-----------------------------------------------------------------------
       !
       !  convert fresh water flux (kg/m^2/s) to virtual salt flux (msu*cm/s):
       !  --------------------------------------------------------------------
       !    ocean reference salinity in (o/oo=psu)
       !    density of freshwater rho_fw = 1.0e3 (kg/m^3)
       !    h2o_flux in (kg/m^2/s) = 1.0e2 (kg/m^2/s)
       !
       !    salt_flux            = - h2o_flux * ocn_ref_salinity / rho_fw
       !    salt_flux (msu*cm/s) = - h2o_flux (kg/m^2/s)
       !                           * ocn_ref_salinity (psu)
       !                           * 1.e-3 (msu/psu)
       !                           * 1.0e2 (kg/m^2/s)/(kg/m^2/s)
       !                           / 1.0e3 (kg/m^3)
       !                         = - h2o_flux (kg/m^2/s)
       !                           * ocn_ref_salinity (psu)
       !                           * fwflux_factor (m/s)(msu/psu)/(kg/m^2/s)
       !
       !    ==>  fwflux_factor = 1.e-6
       !
       !    salt_flux(msu*cm/s) = h2oflux(kg/m^2/s) * salinity_factor
       !
       !    ==> salinity_factor = - ocn_ref_salinity(psu) * fwflux_factor
       !
       !-----------------------------------------------------------------------

       fwflux_factor   = 1.e-6_RKIND
       salinity_factor = -ocn_ref_salinity*fwflux_factor

       !-----------------------------------------------------------------------
       !
       !  convert salt flux (kg/m^2/s) to salt flux (psu*m/s):
       !  -----------------------------------------------------
       !    density of freshwater rho_fw = 1.0e3 (kg/m^3)
       !    salt_flux_kg in (kg/m^2/s) = 1.0e2 (kg/m^2/s)
       !
       !    salt_flux            = - h2o_flux * ocn_ref_salinity / rho_fw
       !    salt_flux (msu*cm/s) = salt_flux_kg (kg/m^2/s)
       !                           * 1.0e2 (kg/m^2/s)/(kg/m^2/s)
       !                           / 1.0e3 (kg/m^3)
       !                         = salt_flux_kg (kg/m^2/s)
       !                           * sflux_factor (psu*m/s)/(kg/m^2/s)
       !
       !    ==>  sflux_factor = 1.0
       !
       !-----------------------------------------------------------------------

       sflux_factor = 1.0_RKIND

       !-----------------------------------------------------------------------
       !
       !  convert fresh water mass flux (kg/m^2/s) to fresh water flux (m/s):
       !  --------------------------------------------------------------------
       !    density of freshwater rho_fw = 1.0e3 (kg/m^3)
       !    h2o_flux in (kg/m^2/s) = 1.0e2 (kg/m^2/s)
       !
       !    fw_flux  = h2o_flux / rho_fw
       !    fw_flux (cm/s) = h2o_flux (kg/m^2/s)
       !                     * 1.0e2 (kg/m^2/s)
       !                     / 1.0e3 (kg/m^3)
       !                   = h2o_flux (kg/m^2/s)
       !                     * fwmass_to_fwflux (m/s)/(kg/m^2/s)
       !
       !    ==>  fwmass_to_fwflux = 0.1
       !
       !-----------------------------------------------------------------------

       fwmass_to_fwflux = 0.1_RKIND

   end subroutine ocn_constants_init!}}}

!***********************************************************************

end module ocn_constants

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
