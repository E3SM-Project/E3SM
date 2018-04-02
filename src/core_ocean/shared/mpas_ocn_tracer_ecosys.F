! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_ecosys
!
!> \brief MPAS ocean ecosys
!> \author Mathew Maltrud
!> \date   11/01/2015
!> \details
!>  This module contains routines for computing tracer forcing due to ecosys
!
!-----------------------------------------------------------------------

module ocn_tracer_ecosys

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use ocn_constants
   use ocn_framework_forcing

!  use BGC_mod
!  use BGC_parms
   use bgc_mod
   use bgc_parms

   implicit none

   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_ecosys_compute, &
             ocn_tracer_ecosys_surface_flux_compute,  &
             ocn_get_ecosysData,  &
             ocn_ecosys_forcing_write_restart,  &
             ocn_tracer_ecosys_init

   integer, public:: &
      numColumnsMax

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  name the necessary BGC derived types
!  all of these are defined in BGC_mod
!-----------------------------------------------------------------------

! autotroph_cnt comes from BGC_parms module
  type(autotroph_type), dimension(autotroph_cnt), public :: autotrophs
  type(BGC_indices_type)    , public :: BGC_indices
  type(BGC_input_type)      , public :: BGC_input
  type(BGC_forcing_type)    , public :: BGC_forcing
  type(BGC_output_type)     , public :: BGC_output
  type(BGC_diagnostics_type), public :: BGC_diagnostic_fields
  type(BGC_flux_diagnostics_type), public :: BGC_flux_diagnostic_fields

! hold indices in tracer pool corresponding to each eco tracer array
  type(BGC_indices_type), public :: ecosysIndices

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_ecosys_compute
!
!> \brief   computes a tracer tendency due to ecosys
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine computes a tracer tendency due to ecosys
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_ecosys_compute(activeTracers, ecosysTracers, forcingPool, nTracers, nCellsSolve, &
      latCell, maxLevelCell, nVertLevels, layerThickness, zMid, indexTemperature, indexSalinity, ecosysTracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         maxLevelCell
      real (kind=RKIND), dimension(:), intent(in) :: &
         latCell

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid, layerThickness

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         ecosysTracers
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve, nVertLevels, indexTemperature, indexSalinity

      type (mpas_pool_type), intent(inout) :: forcingPool

      !
      ! two dimensional pointers
      !
      real (kind=RKIND), dimension(:), pointer :: &
         dust_FLUX_IN, PAR_surface, shortWaveHeatFlux

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), pointer :: &
         PH_PREV_3D, PH_PREV_ALT_CO2_3D, FESEDFLUX

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        ecosysTracersTend

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: ecosysAuxiliary  ! additional forcing fields
      type (mpas_pool_type), pointer :: ecosysDiagFieldsLevel1,  &  ! diagnostics
                                        ecosysDiagFieldsLevel2,  &
                                        ecosysDiagFieldsLevel3,  &
                                        ecosysDiagFieldsLevel4,  &
                                        ecosysDiagFieldsLevel5

      logical, pointer ::  &
         config_ecosys_atm_alt_co2_use_eco

      logical, pointer ::  &
         config_ecosysTracers_diagnostic_fields_level1,  &
         config_ecosysTracers_diagnostic_fields_level2,  &
         config_ecosysTracers_diagnostic_fields_level3,  &
         config_ecosysTracers_diagnostic_fields_level4,  &
         config_ecosysTracers_diagnostic_fields_level5

      !
      ! level 1 diagnostics
      !

      real (kind=RKIND), dimension(:), pointer :: &
         ecosys_diag_photoC_TOT_zint,  &
         ecosys_diag_photoC_NO3_TOT_zint,  &
         ecosys_diag_O2_ZMIN,  &
         ecosys_diag_O2_ZMIN_DEPTH,  &
         ecosys_diag_Chl_TOT_zint_100m,  &
         ecosys_diag_Jint_Ctot,  &
         ecosys_diag_Jint_100m_Ctot,  &
         ecosys_diag_Jint_Ntot,  &
         ecosys_diag_Jint_100m_Ntot,  &
         ecosys_diag_Jint_Ptot,  &
         ecosys_diag_Jint_100m_Ptot,  &
         ecosys_diag_Jint_Sitot,  &
         ecosys_diag_Jint_100m_Sitot

      ! one dimensional arrays for each autotroph (nx,autotroph_cnt)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_photoC_zint,    &
         ecosys_diag_photoC_NO3_zint

      ! two dimensional arrays (nx,nz)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_PAR_avg,  &
         ecosys_diag_POC_FLUX_IN,  &
         ecosys_diag_CaCO3_FLUX_IN,  &
         ecosys_diag_auto_graze_TOT,  &
         ecosys_diag_zoo_loss,  &
         ecosys_diag_photoC_TOT,  &
         ecosys_diag_photoC_NO3_TOT,  &
         ecosys_diag_NITRIF,  &
         ecosys_diag_DENITRIF,  &
         ecosys_diag_calcToSed,      &
         ecosys_diag_pocToSed,       &
         ecosys_diag_pfeToSed,       &
         ecosys_diag_SedDenitrif,    &
         ecosys_diag_tot_Nfix

      !
      ! level 2 diagnostics
      !

      ! two dimensional arrays (nx,nz)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_O2_PRODUCTION,  &
         ecosys_diag_O2_CONSUMPTION,  &
         ecosys_diag_AOU,          &
         ecosys_diag_pH_3D,     &
         ecosys_diag_POC_PROD,  &
         ecosys_diag_POC_REMIN,  &
!maltrud NOT SAVED IN BGC_mod YET
         ecosys_diag_POC_ACCUM

      ! two dimensional arrays for each autotroph (nx,nz,autotroph_cnt)
      real (BGC_r8), dimension(:,:,:), pointer :: &
         ecosys_diag_N_lim,          &
         ecosys_diag_P_lim,          &
         ecosys_diag_Fe_lim,         &
         ecosys_diag_SiO3_lim,       &
         ecosys_diag_light_lim,      &
         ecosys_diag_photoC,         &
         ecosys_diag_photoC_NO3,     &
         ecosys_diag_photoFe,        &
         ecosys_diag_photoNO3,       &
         ecosys_diag_photoNH4,       &
         ecosys_diag_DOP_uptake,     &
         ecosys_diag_PO4_uptake,     &
         ecosys_diag_auto_graze,     &
         ecosys_diag_auto_loss,      &
         ecosys_diag_auto_agg,       &
         ecosys_diag_Nfix

      !
      ! level 3 diagnostics
      !

      ! one dimensional arrays (nx)
      real (kind=RKIND), dimension(:), pointer :: &
         ecosys_diag_tot_bSi_form

      ! two dimensional arrays (nx,nz)
      real (BGC_r8), dimension(:,:), pointer :: &
         ecosys_diag_SiO2_FLUX_IN,  &
         ecosys_diag_SiO2_PROD,  &
         ecosys_diag_SiO2_REMIN,  &
         ecosys_diag_dust_FLUX_IN,  &
         ecosys_diag_dust_REMIN,  &
         ecosys_diag_P_iron_FLUX_IN,  &
         ecosys_diag_P_iron_PROD,  &
         ecosys_diag_P_iron_REMIN,  &
         ecosys_diag_DOC_prod,  &
         ecosys_diag_DOC_remin,  &
         ecosys_diag_DON_prod,  &
         ecosys_diag_DON_remin,  &
         ecosys_diag_DOP_prod,  &
         ecosys_diag_DOP_remin,  &
         ecosys_diag_DOFe_prod,  &
         ecosys_diag_DOFe_remin,  &
         ecosys_diag_Fe_scavenge,  &
         ecosys_diag_Fe_scavenge_rate,  &
!maltrud NOT SAVED IN BGC_mod YET
         ecosys_diag_DONr_remin,     &
         ecosys_diag_DOPr_remin,     &
         ecosys_diag_ponToSed,       &
         ecosys_diag_popToSed,       &
         ecosys_diag_bsiToSed,       &
         ecosys_diag_dustToSed,      &
         ecosys_diag_OtherRemin

      ! two dimensional arrays for each autotroph (nx,nz,autotroph_cnt)
      real (BGC_r8), dimension(:,:,:), pointer :: &
         ecosys_diag_bSi_form

      !
      ! level 4 diagnostics
      !

      ! one dimensional arrays (nx)
      real (BGC_r8), dimension(:), pointer :: &
         ecosys_diag_tot_CaCO3_form_zint, &
         ecosys_diag_zsatcalc,            &
         ecosys_diag_zsatarag

      ! one dimensional arrays for each autotroph (nx,autotroph_cnt)
      real (BGC_r8), dimension(:,:), pointer :: &
         ecosys_diag_CaCO3_form_zint  ! diag array for CaCO3 formation vertical integral

      ! two dimensional arrays (nx,nz)
      real (BGC_r8), dimension(:,:), pointer :: &
         ecosys_diag_CaCO3_PROD,     &
         ecosys_diag_CaCO3_REMIN,    &
         ecosys_diag_tot_CaCO3_form, &
         ecosys_diag_CO3,            &
         ecosys_diag_HCO3,           &
         ecosys_diag_H2CO3,          &
         ecosys_diag_CO3_ALT_CO2,    &
         ecosys_diag_HCO3_ALT_CO2,   &
         ecosys_diag_H2CO3_ALT_CO2,  &
         ecosys_diag_pH_3D_ALT_CO2,  &
         ecosys_diag_co3_sat_calc,   &
         ecosys_diag_co3_sat_arag

      ! two dimensional arrays for each autotroph (nx,nz,autotroph_cnt)
      real (BGC_r8), dimension(:,:,:), pointer :: &
         ecosys_diag_CaCO3_form

      !
      !
      ! level 5 diagnostics
      !

      ! two dimensional arrays (nx,nz)
      real (BGC_r8), dimension(:,:), pointer :: &
         ecosys_diag_PO4_RESTORE,  &
         ecosys_diag_NO3_RESTORE,  &
         ecosys_diag_SiO3_RESTORE

      ! source/sink wants CGS units
      ! then convert back to MKS after
      real (kind=RKIND) :: convertLengthMKStoCGS = 100.0_RKIND
      real (kind=RKIND) :: convertLengthCGStoMKS = 0.01_RKIND
      real (kind=RKIND) :: zTop, zBot

      integer :: iCell, iLevel, iTracer, numColumns, column, autotroph

      err = 0

      call mpas_timer_start("ecosys source-sink")

      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_3D', PH_PREV_3D)
      call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_ALT_CO2_3D', PH_PREV_ALT_CO2_3D)
      call mpas_pool_get_array(ecosysAuxiliary, 'FESEDFLUX', FESEDFLUX)
      call mpas_pool_get_array(ecosysAuxiliary, 'dust_FLUX_IN', dust_FLUX_IN)
      call mpas_pool_get_array(ecosysAuxiliary, 'PAR_surface', PAR_surface)

      call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)

      call mpas_pool_get_config(ocnConfigs, 'config_ecosys_atm_alt_co2_use_eco', &
              config_ecosys_atm_alt_co2_use_eco)

!maltrud change to diagnostics pool at some point (needs to be passed in)

      call mpas_pool_get_config(ocnConfigs, 'config_ecosysTracers_diagnostic_fields_level1', &
              config_ecosysTracers_diagnostic_fields_level1)
      call mpas_pool_get_config(ocnConfigs, 'config_ecosysTracers_diagnostic_fields_level2', &
              config_ecosysTracers_diagnostic_fields_level2)
      call mpas_pool_get_config(ocnConfigs, 'config_ecosysTracers_diagnostic_fields_level3', &
              config_ecosysTracers_diagnostic_fields_level3)
      call mpas_pool_get_config(ocnConfigs, 'config_ecosysTracers_diagnostic_fields_level4', &
              config_ecosysTracers_diagnostic_fields_level4)
      call mpas_pool_get_config(ocnConfigs, 'config_ecosysTracers_diagnostic_fields_level5', &
              config_ecosysTracers_diagnostic_fields_level5)

      if (config_ecosysTracers_diagnostic_fields_level1) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel1', ecosysDiagFieldsLevel1)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_TOT_zint', ecosys_diag_photoC_TOT_zint)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_NO3_TOT_zint', ecosys_diag_photoC_NO3_TOT_zint)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_O2_ZMIN', ecosys_diag_O2_ZMIN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_O2_ZMIN_DEPTH', ecosys_diag_O2_ZMIN_DEPTH)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Chl_TOT_zint_100m', ecosys_diag_Chl_TOT_zint_100m)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Ctot', ecosys_diag_Jint_Ctot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_100m_Ctot', ecosys_diag_Jint_100m_Ctot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Ntot', ecosys_diag_Jint_Ntot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_100m_Ntot', ecosys_diag_Jint_100m_Ntot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Ptot', ecosys_diag_Jint_Ptot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_100m_Ptot', ecosys_diag_Jint_100m_Ptot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Sitot', ecosys_diag_Jint_Sitot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_100m_Sitot', ecosys_diag_Jint_100m_Sitot)

         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_zint', ecosys_diag_photoC_zint)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_NO3_zint', ecosys_diag_photoC_NO3_zint)

         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_PAR_avg', ecosys_diag_PAR_avg)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_POC_FLUX_IN', ecosys_diag_POC_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_CaCO3_FLUX_IN', ecosys_diag_CaCO3_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_auto_graze_TOT', ecosys_diag_auto_graze_TOT)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_zoo_loss', ecosys_diag_zoo_loss)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_TOT', ecosys_diag_photoC_TOT)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_NO3_TOT', ecosys_diag_photoC_NO3_TOT)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_NITRIF', ecosys_diag_NITRIF)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_DENITRIF', ecosys_diag_DENITRIF)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_calcToSed', ecosys_diag_calcToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_pocToSed', ecosys_diag_pocToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_pfeToSed', ecosys_diag_pfeToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_SedDenitrif', ecosys_diag_SedDenitrif)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_tot_Nfix', ecosys_diag_tot_Nfix)

      endif

      if (config_ecosysTracers_diagnostic_fields_level2) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel2', ecosysDiagFieldsLevel2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_O2_PRODUCTION', ecosys_diag_O2_PRODUCTION)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_O2_CONSUMPTION', ecosys_diag_O2_CONSUMPTION)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_AOU', ecosys_diag_AOU)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pH_3D', ecosys_diag_pH_3D)

         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POC_PROD', ecosys_diag_POC_PROD)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POC_REMIN', ecosys_diag_POC_REMIN)
!maltrud NOT SAVED IN BGC_mod YET
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POC_ACCUM', ecosys_diag_POC_ACCUM)

         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_N_lim', ecosys_diag_N_lim)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_P_lim', ecosys_diag_P_lim)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Fe_lim', ecosys_diag_Fe_lim)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_SiO3_lim', ecosys_diag_SiO3_lim)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_light_lim', ecosys_diag_light_lim)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoC', ecosys_diag_photoC)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoC_NO3', ecosys_diag_photoC_NO3)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoFe', ecosys_diag_photoFe)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoNO3', ecosys_diag_photoNO3)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoNH4', ecosys_diag_photoNH4)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_DOP_uptake', ecosys_diag_DOP_uptake)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_PO4_uptake', ecosys_diag_PO4_uptake)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_graze', ecosys_diag_auto_graze)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_loss', ecosys_diag_auto_loss)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_agg', ecosys_diag_auto_agg)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Nfix', ecosys_diag_Nfix)
      endif

      if (config_ecosysTracers_diagnostic_fields_level3) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel3', ecosysDiagFieldsLevel3)

         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_tot_bSi_form', ecosys_diag_tot_bSi_form)

         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_SiO2_FLUX_IN', ecosys_diag_SiO2_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_SiO2_PROD', ecosys_diag_SiO2_PROD)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_SiO2_REMIN', ecosys_diag_SiO2_REMIN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_dust_FLUX_IN', ecosys_diag_dust_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_dust_REMIN', ecosys_diag_dust_REMIN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_P_iron_FLUX_IN', ecosys_diag_P_iron_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_P_iron_PROD', ecosys_diag_P_iron_PROD)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_P_iron_REMIN', ecosys_diag_P_iron_REMIN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOC_prod', ecosys_diag_DOC_prod)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOC_remin', ecosys_diag_DOC_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DON_prod', ecosys_diag_DON_prod)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DON_remin', ecosys_diag_DON_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOP_prod', ecosys_diag_DOP_prod)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOP_remin', ecosys_diag_DOP_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOFe_prod', ecosys_diag_DOFe_prod)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOFe_remin', ecosys_diag_DOFe_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Fe_scavenge', ecosys_diag_Fe_scavenge)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Fe_scavenge_rate', ecosys_diag_Fe_scavenge_rate)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DONr_remin', ecosys_diag_DONr_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOPr_remin', ecosys_diag_DOPr_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_ponToSed', ecosys_diag_ponToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_popToSed', ecosys_diag_popToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_bsiToSed', ecosys_diag_bsiToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_dustToSed', ecosys_diag_dustToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_OtherRemin', ecosys_diag_OtherRemin)

         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_bSi_form', ecosys_diag_bSi_form)
      endif

      if (config_ecosysTracers_diagnostic_fields_level4) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel4', ecosysDiagFieldsLevel4)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_tot_CaCO3_form_zint', ecosys_diag_tot_CaCO3_form_zint)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_zsatcalc', ecosys_diag_zsatcalc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_zsatarag', ecosys_diag_zsatarag)

         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CaCO3_form_zint', ecosys_diag_CaCO3_form_zint)

         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CaCO3_PROD', ecosys_diag_CaCO3_PROD)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CaCO3_REMIN', ecosys_diag_CaCO3_REMIN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_tot_CaCO3_form', ecosys_diag_tot_CaCO3_form)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CaCO3_form', ecosys_diag_CaCO3_form)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CO3', ecosys_diag_CO3)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_HCO3', ecosys_diag_HCO3)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_H2CO3', ecosys_diag_H2CO3)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CO3_ALT_CO2', ecosys_diag_CO3_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_HCO3_ALT_CO2', ecosys_diag_HCO3_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_H2CO3_ALT_CO2', ecosys_diag_H2CO3_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_pH_3D_ALT_CO2', ecosys_diag_pH_3D_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_co3_sat_calc', ecosys_diag_co3_sat_calc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_co3_sat_arag', ecosys_diag_co3_sat_arag)
      endif

      if (config_ecosysTracers_diagnostic_fields_level5) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel5', ecosysDiagFieldsLevel5)
         call mpas_pool_get_array(ecosysDiagFieldsLevel5, 'ecosys_diag_PO4_RESTORE', ecosys_diag_PO4_RESTORE)
         call mpas_pool_get_array(ecosysDiagFieldsLevel5, 'ecosys_diag_NO3_RESTORE', ecosys_diag_NO3_RESTORE)
         call mpas_pool_get_array(ecosysDiagFieldsLevel5, 'ecosys_diag_SiO3_RESTORE', ecosys_diag_SiO3_RESTORE)
      endif

      numColumns = 1
      column = 1
      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve
         BGC_input%number_of_active_levels(column) = maxLevelCell(iCell)
         BGC_input%cell_latitude(column) = latCell(iCell)
! convert dust flux from kg/m2/s to g/cm2/s
         BGC_forcing%dust_FLUX_IN(column) = dust_FLUX_IN(iCell)*0.1_RKIND
         BGC_forcing%ShortWaveFlux_surface(column)  = shortWaveHeatFlux(iCell)
         zTop = 0.0_RKIND
         do iLevel=1,maxLevelCell(iCell)
            BGC_input%PotentialTemperature(iLevel,column) = activeTracers(indexTemperature,iLevel,iCell)
            BGC_input%Salinity(iLevel,column) = activeTracers(indexSalinity,iLevel,iCell)
            BGC_input%cell_center_depth(iLevel,column) = -1.0_RKIND*zMid(iLevel,iCell)*convertLengthMKStoCGS
            BGC_input%cell_thickness(iLevel,column)    = layerThickness(iLevel,iCell)*convertLengthMKStoCGS
            zBot = zTop + layerThickness(iLevel,iCell)
            BGC_input%cell_bottom_depth(iLevel,column) = zBot*convertLengthMKStoCGS
            zTop = zBot

            BGC_output%PH_PREV_3D(iLevel,column) = PH_PREV_3D(iLevel,iCell)
            BGC_output%PH_PREV_ALT_CO2_3D(iLevel,column) = PH_PREV_ALT_CO2_3D(iLevel,iCell)

!maltrud increase FESEDFLUX by factor of 5
!           BGC_forcing%FESEDFLUX(iLevel,column) = FESEDFLUX(iLevel,iCell)*convertLengthMKStoCGS
            BGC_forcing%FESEDFLUX(iLevel,column) = FESEDFLUX(iLevel,iCell)*convertLengthMKStoCGS*5.0_RKIND
            BGC_forcing%NUTR_RESTORE_RTAU(iLevel,column) = 0.0_RKIND
            BGC_forcing%NO3_CLIM(iLevel,column) = 0.0_RKIND
            BGC_forcing%PO4_CLIM(iLevel,column) = 0.0_RKIND
            BGC_forcing%SiO3_CLIM(iLevel,column) = 0.0_RKIND

            BGC_input%BGC_tracers(iLevel,column,BGC_indices%po4_ind)         =  &
               ecosysTracers(ecosysIndices%po4_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%no3_ind)         =  &
               ecosysTracers(ecosysIndices%no3_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%sio3_ind)        =  &
               ecosysTracers(ecosysIndices%sio3_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%nh4_ind)         =  &
               ecosysTracers(ecosysIndices%nh4_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%fe_ind)          =  &
               ecosysTracers(ecosysIndices%fe_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%o2_ind)          =  &
               ecosysTracers(ecosysIndices%o2_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%dic_ind)         =  &
               ecosysTracers(ecosysIndices%dic_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%dic_alt_co2_ind) =  &
               ecosysTracers(ecosysIndices%dic_alt_co2_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%alk_ind)         =  &
               ecosysTracers(ecosysIndices%alk_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%doc_ind)         =  &
               ecosysTracers(ecosysIndices%doc_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%don_ind)         =  &
               ecosysTracers(ecosysIndices%don_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%dofe_ind)        =  &
               ecosysTracers(ecosysIndices%dofe_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%dop_ind)         =  &
               ecosysTracers(ecosysIndices%dop_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%donr_ind)        =  &
               ecosysTracers(ecosysIndices%donr_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%dopr_ind)        =  &
               ecosysTracers(ecosysIndices%dopr_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%zooC_ind)        =  &
               ecosysTracers(ecosysIndices%zooC_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%spC_ind)         =  &
               ecosysTracers(ecosysIndices%spC_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%spChl_ind)       =  &
               ecosysTracers(ecosysIndices%spChl_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%spFe_ind)        =  &
               ecosysTracers(ecosysIndices%spFe_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%spCaCO3_ind)     =  &
               ecosysTracers(ecosysIndices%spCaCO3_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%diatC_ind)       =  &
               ecosysTracers(ecosysIndices%diatC_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%diatChl_ind)     =  &
               ecosysTracers(ecosysIndices%diatChl_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%diatFe_ind)      =  &
               ecosysTracers(ecosysIndices%diatFe_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%diatSi_ind)      =  &
               ecosysTracers(ecosysIndices%diatSi_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%phaeoC_ind)      =  &
               ecosysTracers(ecosysIndices%phaeoC_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%phaeoChl_ind)    =  &
               ecosysTracers(ecosysIndices%phaeoChl_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%phaeoFe_ind)     =  &
               ecosysTracers(ecosysIndices%phaeoFe_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%diazC_ind)       =  &
               ecosysTracers(ecosysIndices%diazC_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%diazChl_ind)     =  &
               ecosysTracers(ecosysIndices%diazChl_ind,iLevel,iCell)
            BGC_input%BGC_tracers(iLevel,column,BGC_indices%diazFe_ind)      =  &
               ecosysTracers(ecosysIndices%diazFe_ind,iLevel,iCell)

         enddo  !  iLevel

         call BGC_SourceSink(autotrophs, BGC_indices, BGC_input, BGC_forcing,   &
                             BGC_output, BGC_diagnostic_fields, nVertLevels,   &
                             numColumnsMax, numColumns, config_ecosys_atm_alt_co2_use_eco)

         do iLevel=1,maxLevelCell(iCell)

            ecosysTracersTend(ecosysIndices%po4_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%po4_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%po4_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%no3_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%no3_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%no3_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%sio3_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%sio3_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%sio3_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%nh4_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%nh4_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%nh4_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%fe_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%fe_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%fe_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%o2_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%o2_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%o2_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%dic_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%dic_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%dic_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%dic_alt_co2_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%dic_alt_co2_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%dic_alt_co2_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%alk_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%alk_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%alk_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%doc_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%doc_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%doc_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%don_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%don_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%don_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%dofe_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%dofe_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%dofe_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%dop_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%dop_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%dop_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%dopr_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%dopr_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%dopr_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%donr_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%donr_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%donr_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%zooC_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%zooC_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%zooC_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%spC_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%spC_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%spC_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%spChl_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%spChl_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%spChl_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%spFe_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%spFe_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%spFe_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%spCaCO3_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%spCaCO3_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%spCaCO3_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%diatC_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%diatC_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%diatC_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%diatChl_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%diatChl_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%diatChl_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%diatFe_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%diatFe_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%diatFe_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%diatSi_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%diatSi_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%diatSi_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%diazC_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%diazC_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%diazC_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%diazChl_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%diazChl_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%diazChl_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%diazFe_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%diazFe_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%diazFe_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%phaeoC_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%phaeoC_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%phaeoC_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%phaeoChl_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%phaeoChl_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%phaeoChl_ind)*layerThickness(iLevel,iCell)
            ecosysTracersTend(ecosysIndices%phaeoFe_ind,iLevel,iCell) =  &
               ecosysTracersTend(ecosysIndices%phaeoFe_ind,iLevel,iCell)   &
               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%phaeoFe_ind)*layerThickness(iLevel,iCell)

            PH_PREV_3D(iLevel,iCell)         = BGC_output%PH_PREV_3D(iLevel,column)
            PH_PREV_ALT_CO2_3D(iLevel,iCell) = BGC_output%PH_PREV_ALT_CO2_3D(iLevel,column)

            !
            ! level 1 diagnostics
            !

            if (config_ecosysTracers_diagnostic_fields_level1) then
               ecosys_diag_PAR_avg(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_PAR_avg(iLevel,column)
               ecosys_diag_POC_FLUX_IN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_POC_FLUX_IN(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_CaCO3_FLUX_IN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_CaCO3_FLUX_IN(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_auto_graze_TOT(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_auto_graze_TOT(iLevel,column)
               ecosys_diag_zoo_loss(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_zoo_loss(iLevel,column)
               ecosys_diag_photoC_TOT(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_photoC_TOT(iLevel,column)
               ecosys_diag_photoC_NO3_TOT(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_photoC_NO3_TOT(iLevel,column)
               ecosys_diag_NITRIF(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_NITRIF(iLevel,column)
               ecosys_diag_DENITRIF(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DENITRIF(iLevel,column)
               ecosys_diag_calcToSed(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_calcToSed(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_pocToSed(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_pocToSed(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_pfeToSed(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_pfeToSed(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_SedDenitrif(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_SedDenitrif(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_tot_Nfix(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_tot_Nfix(iLevel,column)
            endif  !  config_ecosysTracers_diagnostic_fields_level1

            !
            ! level 2 diagnostics
            !
            if (config_ecosysTracers_diagnostic_fields_level2) then
               ecosys_diag_O2_PRODUCTION(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_O2_PRODUCTION(iLevel,column)
               ecosys_diag_O2_CONSUMPTION(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_O2_CONSUMPTION(iLevel,column)
               ecosys_diag_AOU(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_AOU(iLevel,column)
               ecosys_diag_pH_3D(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_pH_3D(iLevel,column)
               ecosys_diag_POC_PROD(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_POC_PROD(iLevel,column)
               ecosys_diag_POC_REMIN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_POC_REMIN(iLevel,column)
!maltrud NOT SAVED IN BGC_mod YET
!              ecosys_diag_POC_ACCUM(iLevel,iCell) =  &
!                          BGC_diagnostic_fields%diag_POC_ACCUM(iLevel,column)
               ecosys_diag_POC_ACCUM(iLevel,iCell) = 0.0_RKIND

               do autotroph = 1, autotroph_cnt
                  ecosys_diag_N_lim(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_N_lim(iLevel,column,autotroph)
                  ecosys_diag_P_lim(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_P_lim(iLevel,column,autotroph)
                  ecosys_diag_Fe_lim(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_Fe_lim(iLevel,column,autotroph)
                  ecosys_diag_SiO3_lim(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_SiO3_lim(iLevel,column,autotroph)
                  ecosys_diag_light_lim(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_light_lim(iLevel,column,autotroph)
                  ecosys_diag_photoC(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_photoC(iLevel,column,autotroph)
                  ecosys_diag_photoC_NO3(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_photoC_NO3(iLevel,column,autotroph)
                  ecosys_diag_photoFe(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_photoFe(iLevel,column,autotroph)
                  ecosys_diag_photoNO3(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_photoNO3(iLevel,column,autotroph)
                  ecosys_diag_photoNH4(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_photoNH4(iLevel,column,autotroph)
                  ecosys_diag_DOP_uptake(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOP_uptake(iLevel,column,autotroph)
                  ecosys_diag_PO4_uptake(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_PO4_uptake(iLevel,column,autotroph)
                  ecosys_diag_auto_graze(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_auto_graze(iLevel,column,autotroph)
                  ecosys_diag_auto_loss(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_auto_loss(iLevel,column,autotroph)
                  ecosys_diag_auto_agg(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_auto_agg(iLevel,column,autotroph)
                  ecosys_diag_Nfix(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_Nfix(iLevel,column,autotroph)
               enddo  !  autotroph loop

            endif  !  config_ecosysTracers_diagnostic_fields_level2

            !
            ! level 3 diagnostics
            !
            if (config_ecosysTracers_diagnostic_fields_level3) then
               ecosys_diag_SiO2_FLUX_IN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_SiO2_FLUX_IN(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_SiO2_PROD(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_SiO2_PROD(iLevel,column)
               ecosys_diag_SiO2_REMIN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_SiO2_REMIN(iLevel,column)
               ecosys_diag_dust_FLUX_IN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_dust_FLUX_IN(iLevel,column)
               ecosys_diag_dust_REMIN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_dust_REMIN(iLevel,column)
               ecosys_diag_P_iron_FLUX_IN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_P_iron_FLUX_IN(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_P_iron_PROD(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_P_iron_PROD(iLevel,column)
               ecosys_diag_P_iron_REMIN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_P_iron_REMIN(iLevel,column)
               ecosys_diag_DOC_prod(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOC_prod(iLevel,column)
               ecosys_diag_DOC_remin(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOC_remin(iLevel,column)
               ecosys_diag_DON_prod(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DON_prod(iLevel,column)
               ecosys_diag_DON_remin(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DON_remin(iLevel,column)
               ecosys_diag_DOP_prod(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOP_prod(iLevel,column)
               ecosys_diag_DOP_remin(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOP_remin(iLevel,column)
               ecosys_diag_DOFe_prod(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOFe_prod(iLevel,column)
               ecosys_diag_DOFe_remin(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOFe_remin(iLevel,column)
               ecosys_diag_Fe_scavenge(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_Fe_scavenge(iLevel,column)
               ecosys_diag_Fe_scavenge_rate(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_Fe_scavenge_rate(iLevel,column)
               ecosys_diag_DONr_remin(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DONr_remin(iLevel,column)
               ecosys_diag_DOPr_remin(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_DOPr_remin(iLevel,column)
               ecosys_diag_ponToSed(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_ponToSed(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_popToSed(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_popToSed(iLevel,column)*convertLengthCGStoMKS
               ecosys_diag_bsiToSed(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_bsiToSed(iLevel,column)*convertLengthCGStoMKS
! convert units from g/cm2/s to kg/m2/s
               ecosys_diag_dustToSed(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_dustToSed(iLevel,column)*1.e-7_RKIND
               ecosys_diag_OtherRemin(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_OtherRemin(iLevel,column)*convertLengthCGStoMKS

               do autotroph = 1, autotroph_cnt
                  ecosys_diag_bSi_form(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_bSi_form(iLevel,column,autotroph)
               enddo  !  autotroph loop

            endif  !  config_ecosysTracers_diagnostic_fields_level3

            !
            ! level 4 diagnostics
            !
            if (config_ecosysTracers_diagnostic_fields_level4) then
               ecosys_diag_CaCO3_PROD(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_CaCO3_PROD(iLevel,column)
               ecosys_diag_CaCO3_REMIN(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_CaCO3_REMIN(iLevel,column)
               ecosys_diag_tot_CaCO3_form(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_tot_CaCO3_form(iLevel,column)
               ecosys_diag_CO3(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_CO3(iLevel,column)
               ecosys_diag_HCO3(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_HCO3(iLevel,column)
               ecosys_diag_H2CO3(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_H2CO3(iLevel,column)
               ecosys_diag_CO3_ALT_CO2(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_CO3_ALT_CO2(iLevel,column)
               ecosys_diag_HCO3_ALT_CO2(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_HCO3_ALT_CO2(iLevel,column)
               ecosys_diag_H2CO3_ALT_CO2(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_H2CO3_ALT_CO2(iLevel,column)
               ecosys_diag_pH_3D_ALT_CO2(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_pH_3D_ALT_CO2(iLevel,column)
               ecosys_diag_co3_sat_calc(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_co3_sat_calc(iLevel,column)
               ecosys_diag_co3_sat_arag(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_co3_sat_arag(iLevel,column)

               do autotroph = 1, autotroph_cnt
                  ecosys_diag_CaCO3_form(autotroph,iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_CaCO3_form(iLevel,column,autotroph)
               enddo
            endif  !  config_ecosysTracers_diagnostic_fields_level4

            !
            ! level 5 diagnostics
            !

            if (config_ecosysTracers_diagnostic_fields_level5) then
                  ecosys_diag_PO4_RESTORE(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_PO4_RESTORE(iLevel,column)
                     ecosys_diag_NO3_RESTORE(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_NO3_RESTORE(iLevel,column)
                     ecosys_diag_SiO3_RESTORE(iLevel,iCell) =  &
                           BGC_diagnostic_fields%diag_SiO3_RESTORE(iLevel,column)
            endif  !  config_ecosysTracers_diagnostic_fields_level5

         enddo  !  iLevel

         !
         ! level 1 diagnostics with no depth dependence
         !
         if (config_ecosysTracers_diagnostic_fields_level1) then
            ecosys_diag_photoC_TOT_zint(iCell) =  &
                        BGC_diagnostic_fields%diag_photoC_TOT_zint(column)*convertLengthCGStoMKS
            ecosys_diag_photoC_NO3_TOT_zint(iCell) =  &
                        BGC_diagnostic_fields%diag_photoC_NO3_TOT_zint(column)*convertLengthCGStoMKS
            ecosys_diag_O2_ZMIN(iCell) =  &
                        BGC_diagnostic_fields%diag_O2_ZMIN(column)
            ecosys_diag_O2_ZMIN_DEPTH(iCell) =  &
                        BGC_diagnostic_fields%diag_O2_ZMIN_DEPTH(column)*convertLengthCGStoMKS
            ecosys_diag_Chl_TOT_zint_100m(iCell) =  &
                        BGC_diagnostic_fields%diag_Chl_TOT_zint_100m(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_Ctot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_Ctot(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_100m_Ctot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_100m_Ctot(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_Ntot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_Ntot(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_100m_Ntot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_100m_Ntot(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_Ptot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_Ptot(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_100m_Ptot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_100m_Ptot(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_Sitot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_Sitot(column)*convertLengthCGStoMKS
            ecosys_diag_Jint_100m_Sitot(iCell) =  &
                        BGC_diagnostic_fields%diag_Jint_100m_Sitot(column)*convertLengthCGStoMKS

            do autotroph = 1, autotroph_cnt
               ecosys_diag_photoC_zint(autotroph,iCell) =  &
                        BGC_diagnostic_fields%diag_photoC_zint(column,autotroph)*convertLengthCGStoMKS
               ecosys_diag_photoC_NO3_zint(autotroph,iCell) =  &
                        BGC_diagnostic_fields%diag_photoC_NO3_zint(column,autotroph)*convertLengthCGStoMKS
            enddo

         endif  !  config_ecosysTracers_diagnostic_fields_level1

         !
         ! level 3 diagnostics with no depth dependence
         !
         if (config_ecosysTracers_diagnostic_fields_level3) then
            ecosys_diag_tot_bSi_form(iCell) =  &
                           BGC_diagnostic_fields%diag_tot_bSi_form(column)
         endif  !  config_ecosysTracers_diagnostic_fields_level3

         !
         ! level 4 diagnostics with no depth dependence
         !
         if (config_ecosysTracers_diagnostic_fields_level4) then
            ecosys_diag_tot_CaCO3_form_zint(iCell) =  &
                           BGC_diagnostic_fields%diag_tot_CaCO3_form_zint(column)*convertLengthCGStoMKS
            ecosys_diag_zsatcalc(iCell) =  &
                           BGC_diagnostic_fields%diag_zsatcalc(column)*convertLengthCGStoMKS
            ecosys_diag_zsatarag(iCell) =  &
                           BGC_diagnostic_fields%diag_zsatarag(column)*convertLengthCGStoMKS

            do autotroph = 1, autotroph_cnt
               ecosys_diag_CaCO3_form_zint(autotroph,iCell) =  &
                           BGC_diagnostic_fields%diag_CaCO3_form_zint(column,autotroph)*convertLengthCGStoMKS
            enddo

         endif  !  config_ecosysTracers_diagnostic_fields_level4

      enddo  !  iCell

      call mpas_timer_stop("ecosys source-sink")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_ecosys_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_ecosys_surface_flux_compute
!
!> \brief   computes a tracer tendency due to ecosys
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine computes a tracer tendency due to ecosys
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_ecosys_surface_flux_compute(activeTracers, ecosysTracers, forcingPool,  &
      nTracers, nCellsSolve, zMid, indexTemperature, indexSalinity, ecosysSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         ecosysSurfaceFlux

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         ecosysTracers
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve, indexTemperature, indexSalinity

      type (mpas_pool_type), intent(inout) :: forcingPool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: ecosysAuxiliary,  &  ! additional forcing fields
                                        ecosysSeaIceCoupling

      integer :: numColumns, column, iCell, iTracer, iLevelSurface

! input flux components in forcing pool
      real (kind=RKIND), dimension(:), pointer :: &
         atmosphericPressure,   &
         iceFraction,          &
         landIceFraction

! input flux components in ecosysAuxiliary
      real (kind=RKIND), dimension(:), pointer :: &
         windSpeedSquared10m,  &
         depositionFluxNO3,    &
         depositionFluxNH4,    &
         IRON_FLUX_IN,         &
         riverFluxNO3,         &
         riverFluxPO4,         &
         riverFluxDON,         &
         riverFluxDOP,         &
         riverFluxSiO3,        &
         riverFluxFe,          &
         riverFluxDIC,         &
         riverFluxALK,         &
         riverFluxDOC,         &
         atmosphericCO2,       &
         atmosphericCO2_ALT_CO2

! input flux components in ecosysSeaIceCoupling
      real (kind=RKIND), dimension(:), pointer :: &
         iceFluxDIC,           &
         iceFluxDON,           &
         iceFluxNO3,           &
         iceFluxSiO3,          &
         iceFluxNH4,           &
         iceFluxDOCr,          &
         iceFluxFeDissolved

   real (kind=RKIND), dimension(:,:), pointer :: iceFluxPhytoC, &
                                                 iceFluxDOC


! specific output fluxes
      real (kind=RKIND), dimension(:), pointer :: &
         CO2_gas_flux,             &
         CO2_alt_gas_flux

! input/output terms
      real (kind=RKIND), dimension(:), pointer :: &
         PH_PREV,              &
         PH_PREV_ALT_CO2

      ! flux routine wants atmospheres
      real (kind=RKIND) :: PascalsToAtmospheres = 1.0_RKIND/101.325e+3_RKIND
      ! flux routine wants CGS units
      ! then convert back to MKS after
      real (kind=RKIND) :: convertLengthMKStoCGS = 100.0_RKIND
      real (kind=RKIND) :: convertLengthCGStoMKS = 0.01_RKIND
      real (kind=RKIND) :: convertLengthSquaredMKStoCGS = 1.0e4_RKIND

      logical, pointer :: config_use_ecosysTracers_sea_ice_coupling

      !
      ! level 2 diagnostics
      !

      type (mpas_pool_type), pointer ::  &
         ecosysDiagFieldsLevel2

      logical, pointer ::  &
         config_ecosysTracers_diagnostic_fields_level2

      real (kind=RKIND), dimension(:), pointer :: &
         ecosys_diag_pistonVel_O2,         &
         ecosys_diag_pistonVel_CO2,        &
         ecosys_diag_Schmidt_O2,           &
         ecosys_diag_Schmidt_CO2,          &
         ecosys_diag_O2_saturation,        &
         ecosys_diag_xkw,                  &
         ecosys_diag_CO2star,              &
         ecosys_diag_dCO2star,             &
         ecosys_diag_pCO2surface,          &
         ecosys_diag_dpCO2,                &
         ecosys_diag_CO2star_ALT_CO2,      &
         ecosys_diag_dCO2star_ALT_CO2,     &
         ecosys_diag_pCO2surface_ALT_CO2,  &
         ecosys_diag_dpCO2_ALT_CO2





      call mpas_timer_start("ecosys surface flux")
      err = 0

      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
      call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)

      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_array(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10m)
      call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV', PH_PREV)
      call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_ALT_CO2', PH_PREV_ALT_CO2)
      call mpas_pool_get_array(ecosysAuxiliary, 'depositionFluxNO3', depositionFluxNO3)
      call mpas_pool_get_array(ecosysAuxiliary, 'depositionFluxNH4', depositionFluxNH4)
      call mpas_pool_get_array(ecosysAuxiliary, 'IRON_FLUX_IN', IRON_FLUX_IN)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxNO3', riverFluxNO3)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxPO4', riverFluxPO4)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDON', riverFluxDON)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOP', riverFluxDOP)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxFe', riverFluxFe)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDIC', riverFluxDIC)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxALK', riverFluxALK)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOC', riverFluxDOC)
      call mpas_pool_get_array(ecosysAuxiliary, 'atmosphericCO2', atmosphericCO2)
      call mpas_pool_get_array(ecosysAuxiliary, 'atmosphericCO2_ALT_CO2', atmosphericCO2_ALT_CO2)

      call mpas_pool_get_array(ecosysAuxiliary, 'CO2_gas_flux', CO2_gas_flux)
      call mpas_pool_get_array(ecosysAuxiliary, 'CO2_alt_gas_flux', CO2_alt_gas_flux)

      call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers_sea_ice_coupling',  &
                                             config_use_ecosysTracers_sea_ice_coupling)

      if (config_use_ecosysTracers_sea_ice_coupling) then

         call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxPhytoC', iceFluxPhytoC)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxDIC', iceFluxDIC)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxNO3', iceFluxNO3)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxSiO3', iceFluxSiO3)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxNH4', iceFluxNH4)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxDOCr', iceFluxDOCr)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxDOC', iceFluxDOC)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxDON', iceFluxDON)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'iceFluxFeDissolved', iceFluxFeDissolved)

      endif

      call mpas_pool_get_config(ocnConfigs, 'config_ecosysTracers_diagnostic_fields_level2', &
              config_ecosysTracers_diagnostic_fields_level2)

      if (config_ecosysTracers_diagnostic_fields_level2) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel2', ecosysDiagFieldsLevel2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pistonVel_O2', ecosys_diag_pistonVel_O2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pistonVel_CO2', ecosys_diag_pistonVel_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Schmidt_O2', ecosys_diag_Schmidt_O2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Schmidt_CO2', ecosys_diag_Schmidt_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_O2_saturation', ecosys_diag_O2_saturation)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_xkw', ecosys_diag_xkw)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_CO2star', ecosys_diag_CO2star)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dCO2star', ecosys_diag_dCO2star)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pCO2surface', ecosys_diag_pCO2surface)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dpCO2', ecosys_diag_dpCO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_CO2star_ALT_CO2', ecosys_diag_CO2star_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dCO2star_ALT_CO2', ecosys_diag_dCO2star_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pCO2surface_ALT_CO2', ecosys_diag_pCO2surface_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dpCO2_ALT_CO2', ecosys_diag_dpCO2_ALT_CO2)
      endif

      BGC_forcing%lcalc_O2_gas_flux  = .true.
      BGC_forcing%lcalc_CO2_gas_flux = .true.

      numColumns = 1
      column = 1
      iLevelSurface = 1
      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve

!        BGC_forcing%surfacePressure(column) = atmosphericPressure(iCell)*PascalsToAtmospheres
         BGC_forcing%surfacePressure(column) = 1.0_RKIND
         BGC_forcing%iceFraction(column) = iceFraction(iCell)
!maltrud assume for now that if there is any land ice, it is all land ice
         if (associated(landIceFraction)) then
            if (landIceFraction(iCell) > 0.0_RKIND) BGC_forcing%iceFraction(column) = 1.0_RKIND
         endif
         BGC_forcing%windSpeedSquared10m(column) = windSpeedSquared10m(iCell)*convertLengthSquaredMKStoCGS
         BGC_forcing%atmCO2(column) = atmosphericCO2(iCell)
         BGC_forcing%atmCO2_ALT_CO2(column) = atmosphericCO2_ALT_CO2(iCell)
         BGC_forcing%surface_pH(column) = PH_PREV(iCell)
         BGC_forcing%surface_pH_alt_co2(column) = PH_PREV_ALT_CO2(iCell)
         BGC_forcing%surfaceDepth(column) = -1.0_RKIND*zMid(iLevelSurface,iCell)
         BGC_forcing%SST(column) = activeTracers(indexTemperature,iLevelSurface,iCell)
         BGC_forcing%SSS(column) = activeTracers(indexSalinity,iLevelSurface,iCell)

         BGC_input%BGC_tracers(1,column,BGC_indices%dic_ind) = ecosysTracers(ecosysIndices%dic_ind,1,iCell)
         BGC_input%BGC_tracers(1,column,BGC_indices%dic_alt_co2_ind) = ecosysTracers(ecosysIndices%dic_alt_co2_ind,1,iCell)
         BGC_input%BGC_tracers(1,column,BGC_indices%alk_ind) = ecosysTracers(ecosysIndices%alk_ind,1,iCell)
         BGC_input%BGC_tracers(1,column,BGC_indices%po4_ind) = ecosysTracers(ecosysIndices%po4_ind,1,iCell)
         BGC_input%BGC_tracers(1,column,BGC_indices%sio3_ind) = ecosysTracers(ecosysIndices%sio3_ind,1,iCell)
         BGC_input%BGC_tracers(1,column,BGC_indices%o2_ind) = ecosysTracers(ecosysIndices%o2_ind,1,iCell)

! NOTE pass in total Fe and mult by parm_Fe_bioavail inside the flux routine
!             divide river Fe by bioavail since it is already the available to make it total

         BGC_forcing%depositionFlux(column,BGC_indices%no3_ind) = depositionFluxNO3(iCell)*convertLengthMKStoCGS
         BGC_forcing%depositionFlux(column,BGC_indices%nh4_ind) = depositionFluxNH4(iCell)*convertLengthMKStoCGS
         BGC_forcing%depositionFlux(column,BGC_indices%fe_ind)  = IRON_FLUX_IN(iCell)*convertLengthMKStoCGS

         BGC_forcing%riverFlux(column,BGC_indices%no3_ind)         = riverFluxNO3(iCell)*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%po4_ind)         = riverFluxPO4(iCell)*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%don_ind)         = riverFluxDON(iCell) * 0.9_BGC_r8*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%donr_ind)        = riverFluxDON(iCell) * 0.1_BGC_r8*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%dop_ind)         = riverFluxDOP(iCell) * 0.975_BGC_r8*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%dopr_ind)        = riverFluxDOP(iCell) * 0.025_BGC_r8*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%sio3_ind)        = riverFluxSiO3(iCell)*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%fe_ind)          = riverFluxFe(iCell)*convertLengthMKStoCGS / parm_Fe_bioavail
         BGC_forcing%riverFlux(column,BGC_indices%dic_ind)         = riverFluxDIC(iCell)*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%dic_alt_co2_ind) = riverFluxDIC(iCell)*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%alk_ind)         = riverFluxALK(iCell)*convertLengthMKStoCGS
         BGC_forcing%riverFlux(column,BGC_indices%doc_ind)         = riverFluxDOC(iCell)*convertLengthMKStoCGS

         if (config_use_ecosysTracers_sea_ice_coupling) then

!maltrud sea ice fluxes are already correct MPAS units, so change to POP units
            BGC_forcing%seaIceFlux(column,BGC_indices%no3_ind)         = iceFluxNO3(iCell)*convertLengthMKStoCGS
            BGC_forcing%seaIceFlux(column,BGC_indices%don_ind)         = iceFluxDON(iCell)*convertLengthMKStoCGS
!maltrud sea ice BGC doesnt care about C or N--this is just a placeholder--so assume N units
            BGC_forcing%seaIceFlux(column,BGC_indices%donr_ind)        = iceFluxDOCr(iCell)*convertLengthMKStoCGS
            BGC_forcing%seaIceFlux(column,BGC_indices%sio3_ind)        = iceFluxSiO3(iCell)*convertLengthMKStoCGS
            BGC_forcing%seaIceFlux(column,BGC_indices%nh4_ind)         = iceFluxNH4(iCell)*convertLengthMKStoCGS
            BGC_forcing%seaIceFlux(column,BGC_indices%fe_ind)          = iceFluxFeDissolved(iCell)*convertLengthMKStoCGS / parm_Fe_bioavail
            BGC_forcing%seaIceFlux(column,BGC_indices%dic_ind)         = iceFluxDIC(iCell)*convertLengthMKStoCGS
            BGC_forcing%seaIceFlux(column,BGC_indices%dic_alt_co2_ind) = iceFluxDIC(iCell)*convertLengthMKStoCGS

            BGC_forcing%seaIceFlux(column,BGC_indices%doc_ind)         = iceFluxDOC(1,iCell)*convertLengthMKStoCGS + &
                                                                         iceFluxDOC(2,iCell)*convertLengthMKStoCGS

            BGC_forcing%seaIceFlux(column,BGC_indices%diatC_ind)       = iceFluxPhytoC(1,iCell)*convertLengthMKStoCGS
            BGC_forcing%seaIceFlux(column,BGC_indices%spC_ind)         = iceFluxPhytoC(2,iCell)*convertLengthMKStoCGS
            BGC_forcing%seaIceFlux(column,BGC_indices%phaeoC_ind)      = iceFluxPhytoC(3,iCell)*convertLengthMKStoCGS

         endif

         call BGC_SurfaceFluxes(BGC_indices, BGC_input, BGC_forcing,   &
                                BGC_flux_diagnostic_fields,   &
                                numColumnsMax, column)

         PH_PREV(iCell)         = BGC_forcing%surface_pH(column)
         PH_PREV_ALT_CO2(iCell) = BGC_forcing%surface_pH_alt_co2(column)

         ecosysSurfaceFlux(ecosysIndices%no3_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%no3_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%po4_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%po4_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%sio3_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%sio3_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%nh4_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%nh4_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%don_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%don_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%donr_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%donr_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%dop_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%dop_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%dopr_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%dopr_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%fe_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%fe_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%alk_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%alk_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%doc_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%doc_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%o2_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%o2_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%dic_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%dic_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%dic_alt_co2_ind,iCell) =   &
            BGC_forcing%netFlux(column,BGC_indices%dic_alt_co2_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%diatC_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%diatC_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%spC_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%spC_ind)*convertLengthCGStoMKS
         ecosysSurfaceFlux(ecosysIndices%phaeoC_ind,iCell) = BGC_forcing%netFlux(column,BGC_indices%phaeoC_ind)*convertLengthCGStoMKS

!explicitly set the rest to 0
         ecosysSurfaceFlux(ecosysIndices%dofe_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%zooC_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%spChl_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%spFe_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%spCaCO3_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%diatChl_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%diatFe_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%diatSi_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%diazC_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%diazChl_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%diazFe_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%phaeoChl_ind,iCell) = 0.0_RKIND
         ecosysSurfaceFlux(ecosysIndices%phaeoFe_ind,iCell) = 0.0_RKIND

         CO2_gas_flux(iCell)     = BGC_forcing%gasFlux(column,BGC_indices%dic_ind)*convertLengthCGStoMKS
         CO2_alt_gas_flux(iCell) = BGC_forcing%gasFlux(column,BGC_indices%dic_alt_co2_ind)*convertLengthCGStoMKS

         if (config_ecosysTracers_diagnostic_fields_level2) then
            ecosys_diag_pistonVel_O2(iCell) = BGC_flux_diagnostic_fields%pistonVel_O2(column)
            ecosys_diag_pistonVel_CO2(iCell) = BGC_flux_diagnostic_fields%pistonVel_CO2(column)
            ecosys_diag_Schmidt_O2(iCell) = BGC_flux_diagnostic_fields%SCHMIDT_O2(column)
            ecosys_diag_Schmidt_CO2(iCell) = BGC_flux_diagnostic_fields%SCHMIDT_CO2(column)
            ecosys_diag_O2_saturation(iCell) = BGC_flux_diagnostic_fields%O2SAT(column)
            ecosys_diag_xkw(iCell) = BGC_flux_diagnostic_fields%xkw(column)
            ecosys_diag_CO2star(iCell) = BGC_flux_diagnostic_fields%co2star(column)
            ecosys_diag_dCO2star(iCell) = BGC_flux_diagnostic_fields%dco2star(column)
            ecosys_diag_pCO2surface(iCell) = BGC_flux_diagnostic_fields%pco2surf(column)
            ecosys_diag_dpCO2(iCell) = BGC_flux_diagnostic_fields%dpco2(column)
            ecosys_diag_CO2star_ALT_CO2(iCell) = BGC_flux_diagnostic_fields%co2star_alt_co2(column)
            ecosys_diag_dCO2star_ALT_CO2(iCell) = BGC_flux_diagnostic_fields%dco2star_alt_co2(column)
            ecosys_diag_pCO2surface_ALT_CO2(iCell) = BGC_flux_diagnostic_fields%pco2surf_alt_co2(column)
            ecosys_diag_dpCO2_ALT_CO2(iCell) = BGC_flux_diagnostic_fields%dpco2_alt_co2(column)
         endif

      enddo  !  iCell

      call mpas_timer_stop("ecosys surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_ecosys_surface_flux_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_ecosys_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_ecosys_init(domain,err)!{{{

!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: ecosysAuxiliary
      type (mpas_pool_type), pointer :: ecosysMonthlyForcing

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        ecosysTracers

! input flux components in ecosysAuxiliary
      real (kind=RKIND), dimension(:), pointer :: &
         depositionFluxNO3,    &
         depositionFluxNH4,    &
         IRON_FLUX_IN,         &
         dust_FLUX_IN,         &
         riverFluxNO3,         &
         riverFluxPO4,         &
         riverFluxDON,         &
         riverFluxDOP,         &
         riverFluxSiO3,        &
         riverFluxFe,          &
         riverFluxDIC,         &
         riverFluxALK,         &
         riverFluxDOC

! input flux components in ecosysMonthlyForcing
      real (kind=RKIND), dimension(:), pointer :: &
         depositionFluzNO3,    &
         depositionFluzNH4,    &
         IRON_FLUZ_IN,         &
         dust_FLUZ_IN,         &
         riverFluzNO3,         &
         riverFluzPO4,         &
         riverFluzDON,         &
         riverFluzDOP,         &
         riverFluzSiO3,        &
         riverFluzFe,          &
         riverFluzDIC,         &
         riverFluzALK,         &
         riverFluzDOC

      ! scalars
      integer :: nTracers, numColumnsMax

      ! scalar pointers
      integer, pointer :: nVertLevels, index_dummy

      character(len=strKIND) :: &
         forcingIntervalMonthly,  &
         forcingReferenceTimeMonthly

      logical, pointer :: &
         config_do_restart


      !
      ! get tracers pools
      !

      err = 0

      !
      ! Get tracer group so we can get the number of tracers in it
      !

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, 1)

      if (associated(ecosysTracers)) then

      nTracers = size(ecosysTracers, dim=1)
      if (BGC_tracer_cnt /= nTracers) then
         err = 1
         return
      endif

      !
      ! pull nVertLevels out of the mesh structure
      !

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

!-----------------------------------------------------------------------
!  initialize ecosystem parameters
!-----------------------------------------------------------------------

   allocate( BGC_indices%short_name(BGC_tracer_cnt) )
   allocate( BGC_indices%long_name(BGC_tracer_cnt) )
   allocate( BGC_indices%units(BGC_tracer_cnt) )

! no need to allocate the above fields for ecosysIndices (?)

!-----------------------------------------------------------------------
!  sets most of BGC parameters
!  sets namelist defaults
!  sets autotroph sp_ind, diat_ind, diaz_ind, phaeo_ind (swang)
!-----------------------------------------------------------------------

   call BGC_parms_init(BGC_indices, autotrophs)

! modify autotroph values here....
!  for example to change sp_kFe
!  autotrophs(BGC_indices%sp_ind)%kFe = 0.05e-3_BGC_r8

!maltrud how to handle this?
   T0_Kelvin_BGC = T0_Kelvin

      !
      ! for now only do 1 column at a time
      !
      numColumnsMax = 1

      BGC_indices%po4_ind         =  1
      BGC_indices%no3_ind         =  2
      BGC_indices%sio3_ind        =  3
      BGC_indices%nh4_ind         =  4
      BGC_indices%fe_ind          =  5
      BGC_indices%o2_ind          =  6
      BGC_indices%dic_ind         =  7
      BGC_indices%dic_alt_co2_ind =  8
      BGC_indices%alk_ind         =  9
      BGC_indices%doc_ind         = 10
      BGC_indices%don_ind         = 11
      BGC_indices%dofe_ind        = 12
      BGC_indices%dop_ind         = 13
      BGC_indices%dopr_ind        = 14
      BGC_indices%donr_ind        = 15
      BGC_indices%zooC_ind        = 16
      BGC_indices%spChl_ind       = 17
      BGC_indices%spC_ind         = 18
      BGC_indices%spFe_ind        = 19
      BGC_indices%spCaCO3_ind     = 20
      BGC_indices%diatChl_ind     = 21
      BGC_indices%diatC_ind       = 22
      BGC_indices%diatFe_ind      = 23
      BGC_indices%diatSi_ind      = 24
      BGC_indices%diazChl_ind     = 25
      BGC_indices%diazC_ind       = 26
      BGC_indices%diazFe_ind      = 27
      BGC_indices%phaeoChl_ind    = 28
      BGC_indices%phaeoC_ind      = 29
      BGC_indices%phaeoFe_ind     = 30

      call mpas_pool_get_dimension(tracersPool, 'index_PO4',      index_dummy)
      ecosysIndices%po4_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_NO3',      index_dummy)
      ecosysIndices%no3_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_SiO3',     index_dummy)
      ecosysIndices%sio3_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_NH4',      index_dummy)
      ecosysIndices%nh4_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_Fe',       index_dummy)
      ecosysIndices%fe_ind              = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_O2',       index_dummy)
      ecosysIndices%o2_ind              = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DIC',      index_dummy)
      ecosysIndices%dic_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DIC_ALT_CO2', index_dummy)
      ecosysIndices%dic_alt_co2_ind     = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_ALK',      index_dummy)
      ecosysIndices%alk_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOC',      index_dummy)
      ecosysIndices%doc_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DON',      index_dummy)
      ecosysIndices%don_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOFe',     index_dummy)
      ecosysIndices%dofe_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOP',      index_dummy)
      ecosysIndices%dop_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOPr',     index_dummy)
      ecosysIndices%dopr_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DONr',     index_dummy)
      ecosysIndices%donr_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_zooC',     index_dummy)
      ecosysIndices%zooC_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spChl',    index_dummy)
      ecosysIndices%spChl_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spC',      index_dummy)
      ecosysIndices%spC_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spFe',     index_dummy)
      ecosysIndices%spFe_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spCaCO3',  index_dummy)
      ecosysIndices%spCaCO3_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatChl',  index_dummy)
      ecosysIndices%diatChl_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatC',    index_dummy)
      ecosysIndices%diatC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatFe',   index_dummy)
      ecosysIndices%diatFe_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatSi',   index_dummy)
      ecosysIndices%diatSi_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazChl',  index_dummy)
      ecosysIndices%diazChl_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazC',    index_dummy)
      ecosysIndices%diazC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazFe',   index_dummy)
      ecosysIndices%diazFe_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoChl', index_dummy)
      ecosysIndices%phaeoChl_ind        = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoC',   index_dummy)
      ecosysIndices%phaeoC_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoFe',  index_dummy)
      ecosysIndices%phaeoFe_ind         = index_dummy

! BGC_init sets short and long names, units in BGC_indices
! also sets autotroph indices within the autotroph derived type

      call BGC_init(BGC_indices, autotrophs)

!NOTES:

!also check short_name with mpas variable name

!-----------------------------------------------------------------------
!  allocate input, forcing, diagnostic arrays
!-----------------------------------------------------------------------

      allocate ( BGC_input%BGC_tracers(nVertLevels, numColumnsMax, BGC_tracer_cnt) )
      allocate ( BGC_input%PotentialTemperature(nVertLevels, numColumnsMax) )
      allocate ( BGC_input%Salinity(nVertLevels, numColumnsMax) )
      allocate ( BGC_input%cell_center_depth(nVertLevels, numColumnsMax) )
      allocate ( BGC_input%cell_thickness(nVertLevels, numColumnsMax) )
      allocate ( BGC_input%cell_bottom_depth(nVertLevels, numColumnsMax) )
      allocate ( BGC_input%number_of_active_levels(numColumnsMax) )
      allocate ( BGC_input%cell_latitude(numColumnsMax) )

      allocate ( BGC_forcing%FESEDFLUX(nVertLevels, numColumnsMax) )
      allocate ( BGC_forcing%NUTR_RESTORE_RTAU(nVertLevels, numColumnsMax) )
      allocate ( BGC_forcing%NO3_CLIM(nVertLevels, numColumnsMax) )
      allocate ( BGC_forcing%PO4_CLIM(nVertLevels, numColumnsMax) )
      allocate ( BGC_forcing%SiO3_CLIM(nVertLevels, numColumnsMax) )

      allocate ( BGC_forcing%dust_FLUX_IN(numColumnsMax) )
      allocate ( BGC_forcing%ShortWaveFlux_surface(numColumnsMax) )
      allocate ( BGC_forcing%surfacePressure(numColumnsMax) )
      allocate ( BGC_forcing%iceFraction(numColumnsMax) )
      allocate ( BGC_forcing%windSpeedSquared10m(numColumnsMax) )
      allocate ( BGC_forcing%atmCO2(numColumnsMax) )
      allocate ( BGC_forcing%atmCO2_ALT_CO2(numColumnsMax) )
      allocate ( BGC_forcing%surface_pH(numColumnsMax) )
      allocate ( BGC_forcing%surface_pH_alt_co2(numColumnsMax) )
      allocate ( BGC_forcing%surfaceDepth(numColumnsMax) )
      allocate ( BGC_forcing%SST(numColumnsMax) )
      allocate ( BGC_forcing%SSS(numColumnsMax) )

      allocate ( BGC_forcing%depositionFlux(numColumnsMax, BGC_tracer_cnt) )
      allocate ( BGC_forcing%riverFlux(numColumnsMax, BGC_tracer_cnt) )
      allocate ( BGC_forcing%gasFlux(numColumnsMax, BGC_tracer_cnt) )
      allocate ( BGC_forcing%seaIceFlux(numColumnsMax, BGC_tracer_cnt) )
      allocate ( BGC_forcing%netFlux(numColumnsMax, BGC_tracer_cnt) )
      BGC_forcing%depositionFlux = 0.0_RKIND
      BGC_forcing%riverFlux      = 0.0_RKIND
      BGC_forcing%gasFlux        = 0.0_RKIND
      BGC_forcing%seaIceFlux     = 0.0_RKIND
      BGC_forcing%netFlux        = 0.0_RKIND

      allocate ( BGC_output%BGC_tendencies(nVertLevels, numColumnsMax, BGC_tracer_cnt) )
      allocate ( BGC_output%PH_PREV_3D(nVertLevels, numColumnsMax) )
      allocate ( BGC_output%PH_PREV_ALT_CO2_3D(nVertLevels, numColumnsMax) )

        !---------------------------------------------------------------------------
        !   allocate flux diagnostic output fields
        !---------------------------------------------------------------------------

      allocate (BGC_flux_diagnostic_fields%pistonVel_O2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%pistonVel_CO2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%SCHMIDT_O2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%SCHMIDT_CO2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%O2SAT(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%xkw(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%co2star(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%dco2star(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%pco2surf(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%dpco2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%co2star_alt_co2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%dco2star_alt_co2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%pco2surf_alt_co2(numColumnsMax) )
      allocate (BGC_flux_diagnostic_fields%dpco2_alt_co2(numColumnsMax) )

      !---------------------------------------------------------------------------
      !   allocate diagnostic output fields
      !---------------------------------------------------------------------------

  ! 3D stuff
      allocate (BGC_diagnostic_fields%diag_tot_Nfix(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_O2_PRODUCTION(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_O2_CONSUMPTION(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_AOU(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_PO4_RESTORE(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_NO3_RESTORE(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_SiO3_RESTORE(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_PAR_avg(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_POC_FLUX_IN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_POC_PROD(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_POC_REMIN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_POC_ACCUM(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_CaCO3_FLUX_IN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_CaCO3_PROD(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_CaCO3_REMIN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_SiO2_FLUX_IN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_SiO2_PROD(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_SiO2_REMIN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_dust_FLUX_IN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_dust_REMIN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_P_iron_FLUX_IN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_P_iron_PROD(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_P_iron_REMIN(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_auto_graze_TOT(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_zoo_loss(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_photoC_TOT(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_photoC_NO3_TOT(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DOC_prod(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DOC_remin(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DON_prod(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DON_remin(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DOFe_prod(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DOFe_remin(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DOP_prod(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DOP_remin(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Fe_scavenge(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Fe_scavenge_rate(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_NITRIF(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DENITRIF(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DONr_remin(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_DOPr_remin(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_CO3(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_HCO3(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_H2CO3(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_pH_3D(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_CO3_ALT_CO2(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_HCO3_ALT_CO2(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_H2CO3_ALT_CO2(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_pH_3D_ALT_CO2(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_co3_sat_calc(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_co3_sat_arag(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_calcToSed(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_pocToSed(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_ponToSed(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_popToSed(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_bsiToSed(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_dustToSed(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_pfeToSed(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_SedDenitrif(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_OtherRemin(nVertLevels, numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_tot_CaCO3_form(nVertLevels, numColumnsMax) )

! 3D stuff for each autotroph
      allocate (BGC_diagnostic_fields%diag_N_lim(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_P_lim(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_Fe_lim(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_SiO3_lim(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_light_lim(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_photoC(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_photoC_NO3(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_photoFe(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_photoNO3(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_photoNH4(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_DOP_uptake(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_PO4_uptake(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_auto_graze(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_auto_loss(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_auto_agg(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_bSi_form(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_CaCO3_form(nVertLevels, numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_Nfix(nVertLevels, numColumnsMax, autotroph_cnt) )

! 2D stuff for each autotroph
      allocate (BGC_diagnostic_fields%diag_photoC_zint(numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_photoC_NO3_zint(numColumnsMax, autotroph_cnt) )
      allocate (BGC_diagnostic_fields%diag_CaCO3_form_zint(numColumnsMax, autotroph_cnt) )

! 2D vertical integrals for photoC
      allocate (BGC_diagnostic_fields%diag_photoC_TOT_zint(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_photoC_NO3_TOT_zint(numColumnsMax) )

! 2D vertical integrals for nutrients
      allocate (BGC_diagnostic_fields%diag_Chl_TOT_zint_100m(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_Ctot(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_100m_Ctot(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_Ntot(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_100m_Ntot(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_Ptot(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_100m_Ptot(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_Sitot(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_Jint_100m_Sitot(numColumnsMax) )

! 2D stuff
      allocate (BGC_diagnostic_fields%diag_tot_bSi_form(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_tot_CaCO3_form_zint(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_zsatcalc(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_zsatarag(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_O2_ZMIN(numColumnsMax) )
      allocate (BGC_diagnostic_fields%diag_O2_ZMIN_DEPTH(numColumnsMax) )

! initialize monthly forcing to be read from file

      forcingIntervalMonthly = "0000-01-00_00:00:00"
      forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

      call MPAS_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

      call MPAS_forcing_init_group( forcingGroupHead,  &
                "ecosysMonthlyClimatology", &
                domain, &
                '0000-01-01_00:00:00', &
                '0000-01-01_00:00:00', &
                '0001-00-00_00:00:00', &
                config_do_restart)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_array(ecosysAuxiliary, 'depositionFluxNO3', depositionFluxNO3)
      call mpas_pool_get_array(ecosysAuxiliary, 'depositionFluxNH4', depositionFluxNH4)
      call mpas_pool_get_array(ecosysAuxiliary, 'IRON_FLUX_IN', IRON_FLUX_IN)
      call mpas_pool_get_array(ecosysAuxiliary, 'dust_FLUX_IN', dust_FLUX_IN)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxNO3', riverFluxNO3)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxPO4', riverFluxPO4)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDON', riverFluxDON)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOP', riverFluxDOP)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxFe', riverFluxFe)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDIC', riverFluxDIC)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxALK', riverFluxALK)
      call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOC', riverFluxDOC)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'ecosysMonthlyForcing', ecosysMonthlyForcing)

      call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNO3', depositionFluzNO3)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNH4', depositionFluzNH4)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'IRON_FLUZ_IN', IRON_FLUZ_IN)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'dust_FLUZ_IN', dust_FLUZ_IN)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzNO3', riverFluzNO3)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzPO4', riverFluzPO4)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDON', riverFluzDON)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOP', riverFluzDOP)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzSiO3', riverFluzSiO3)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzFe', riverFluzFe)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDIC', riverFluzDIC)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzALK', riverFluzALK)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOC', riverFluzDOC)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'depositionFluzNO3', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'depositionFluzNO3',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'depositionFluzNH4', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'depositionFluzNH4',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'IRON_FLUZ_IN', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'IRON_FLUZ_IN',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'dust_FLUZ_IN', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'dust_FLUZ_IN',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzNO3', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzNO3',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzPO4', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzPO4',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzDON', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzDON',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzDOP', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzDOP',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzSiO3', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzSiO3',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzFe', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzFe',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzDIC', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzDIC',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzALK', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzALK',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'ecosysMonthlyClimatology', &
                'riverFluzDOC', &
                'ecosys_monthly_flux', &
                'ecosysMonthlyForcing',  &
                'riverFluzDOC',  &
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

      call MPAS_forcing_init_field_data( forcingGroupHead, &
        'ecosysMonthlyClimatology',  &
         domain % streamManager,  &
         config_do_restart, &
         .false.)

      end if  !  associated(ecosysTracers)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_ecosys_init!}}}

!***********************************************************************

!***********************************************************************
!
!  routine get_ecosysData
!
!> \brief   retrieve data needed to compute ecosys deposition and runoff inputs
!> \author  Mathew Maltrud
!> \date    03/07/16
!> \details
!>  This routine calls mpas_forcing routines to acquire needed ecosys forcing data and interpolates
!>    between time levels.  directly copied from ocn_get_shortWaveData.
!
!-----------------------------------------------------------------------

    subroutine ocn_get_ecosysData( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepEco
        type (MPAS_clock_type) :: simulationClock

        logical, intent(in) :: firstTimeStep
        character(len=strKind), pointer :: config_dt
        real(kind=RKIND) :: dt

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: ecosysAuxiliary
        type (mpas_pool_type), pointer :: ecosysMonthlyForcing

        real (kind=RKIND), dimension(:), pointer :: &
         depositionFluxNO3,    &
         depositionFluxNH4,    &
         IRON_FLUX_IN,         &
         dust_FLUX_IN,         &
         riverFluxNO3,         &
         riverFluxPO4,         &
         riverFluxDON,         &
         riverFluxDOP,         &
         riverFluxSiO3,        &
         riverFluxFe,          &
         riverFluxDIC,         &
         riverFluxALK,         &
         riverFluxDOC

! input flux components in ecosysMonthlyForcing
        real (kind=RKIND), dimension(:), pointer :: &
         depositionFluzNO3,    &
         depositionFluzNH4,    &
         IRON_FLUZ_IN,         &
         dust_FLUZ_IN,         &
         riverFluzNO3,         &
         riverFluzPO4,         &
         riverFluzDON,         &
         riverFluzDOP,         &
         riverFluzSiO3,        &
         riverFluzFe,          &
         riverFluzDIC,         &
         riverFluzALK,         &
         riverFluzDOC

        integer, pointer :: nCells
        integer :: iCell

        logical, pointer :: config_do_restart

        call MPAS_pool_get_config(domain%configs, 'config_dt', config_dt)
        call MPAS_pool_get_config(domain%configs, 'config_do_restart', config_do_restart)

        call mpas_set_timeInterval(timeStepEco,timeString=config_dt)
        call mpas_get_timeInterval(timeStepEco,dt=dt)

!maltrud debug
        if (firstTimeStep .and. config_do_restart) then
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'ecosysMonthlyClimatology', streamManager, 0.0_RKIND)
        else
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'ecosysMonthlyClimatology', streamManager, dt)
        endif

        call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'ecosysMonthlyForcing', ecosysMonthlyForcing)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(ecosysAuxiliary, 'depositionFluxNO3', depositionFluxNO3)
        call mpas_pool_get_array(ecosysAuxiliary, 'depositionFluxNH4', depositionFluxNH4)
        call mpas_pool_get_array(ecosysAuxiliary, 'IRON_FLUX_IN', IRON_FLUX_IN)
        call mpas_pool_get_array(ecosysAuxiliary, 'dust_FLUX_IN', dust_FLUX_IN)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxNO3', riverFluxNO3)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxPO4', riverFluxPO4)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDON', riverFluxDON)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOP', riverFluxDOP)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxFe', riverFluxFe)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDIC', riverFluxDIC)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxALK', riverFluxALK)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOC', riverFluxDOC)

        call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNO3', depositionFluzNO3)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNH4', depositionFluzNH4)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'IRON_FLUZ_IN', IRON_FLUZ_IN)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'dust_FLUZ_IN', dust_FLUZ_IN)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzNO3', riverFluzNO3)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzPO4', riverFluzPO4)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDON', riverFluzDON)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOP', riverFluzDOP)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzSiO3', riverFluzSiO3)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzFe', riverFluzFe)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDIC', riverFluzDIC)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzALK', riverFluzALK)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOC', riverFluzDOC)

        do iCell = 1, nCells
           depositionFluxNO3(iCell) = depositionFluzNO3(iCell)
           depositionFluxNH4(iCell) = depositionFluzNH4(iCell)
           IRON_FLUX_IN(iCell)      = IRON_FLUZ_IN(iCell)
           dust_FLUX_IN(iCell)      = dust_FLUZ_IN(iCell)
           riverFluxNO3(iCell)      = riverFluzNO3(iCell)
           riverFluxPO4(iCell)      = riverFluzPO4(iCell)
           riverFluxDON(iCell)      = riverFluzDON(iCell)
           riverFluxDOP(iCell)      = riverFluzDOP(iCell)
           riverFluxSiO3(iCell)     = riverFluzSiO3(iCell)
           riverFluxFe(iCell)       = riverFluzFe(iCell)
           riverFluxDIC(iCell)      = riverFluzDIC(iCell)
           riverFluxALK(iCell)      = riverFluzALK(iCell)
           riverFluxDOC(iCell)      = riverFluzDOC(iCell)
        enddo

    end subroutine ocn_get_ecosysData!}}}

!***********************************************************************
!
!  routine ocn_ecosys_forcing_write_restart
!
!> \brief   writes restart timestamp for ecosys data to be read in on future restart
!> \author  Mathew Maltrud
!> \date    03/07/2016

!
!-----------------------------------------------------------------------

   subroutine ocn_ecosys_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_ecosys_forcing_write_restart!}}}

end module ocn_tracer_ecosys

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
