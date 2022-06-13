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
   use mpas_log
   use ocn_constants
   use ocn_config
   use ocn_framework_forcing
   use ocn_diagnostics_variables

   use marbl_logging             , only : marbl_log_type
   use marbl_interface           , only : marbl_interface_class

   implicit none

   type(marbl_interface_class) :: marbl_instances(1)
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

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer, parameter ::  &
      ecosys_tracer_cnt = 32,  &
      autotroph_cnt     = 3,   &
      zooplankton_cnt   = 1

   integer ::  &
      ecosys_interior_diagnostics_cnt, &
      ecosys_surface_forcing_diagnostics_cnt

! Bioavailability is already included in surface Fe flux
   real (kind=RKIND) :: &
      parm_Fe_bioavail = 1.0_RKIND

   integer, public:: &
      numColumnsMax

   integer, public ::         &
      po4_ind_MPAS,           &
      no3_ind_MPAS,           &
      sio3_ind_MPAS,          &
      nh4_ind_MPAS,           &
      fe_ind_MPAS,            &
      lig_ind_MPAS,           &
      o2_ind_MPAS,            &
      dic_ind_MPAS,           &
      dic_alt_co2_ind_MPAS,   &
      alk_ind_MPAS,           &
      alk_alt_co2_ind_MPAS,   &
      doc_ind_MPAS,           &
      don_ind_MPAS,           &
      dop_ind_MPAS,           &
      dopr_ind_MPAS,          &
      donr_ind_MPAS,          &
      docr_ind_MPAS,          &
      zooC_ind_MPAS,          &
      spChl_ind_MPAS,         &
      spC_ind_MPAS,           &
      spFe_ind_MPAS,          &
      spP_ind_MPAS,           &
      spCaCO3_ind_MPAS,       &
      diatChl_ind_MPAS,       &
      diatC_ind_MPAS,         &
      diatFe_ind_MPAS,        &
      diatP_ind_MPAS,         &
      diatSi_ind_MPAS,        &
      diazChl_ind_MPAS,       &
      diazC_ind_MPAS,         &
      diazFe_ind_MPAS,        &
      diazP_ind_MPAS

      integer, allocatable, dimension(:) :: &
         ecosysIndexMARBLtoMPASO,  &
         ecosysIndexMPASOtoMARBL

      integer ::  &
         dustFlux_indexInMARBL,          &
         shortwave_indexInMARBL,         &
         PAR_col_frac_indexInMARBL,      &
         temperature_indexInMARBL,       &
         salinity_indexInMARBL,          &
         pressure_indexInMARBL,          &
         fesedFlux_indexInMARBL,         &
         u10sqr_indexInMARBL,            &
         ifrac_indexInMARBL,             &
         sst_indexInMARBL,               &
         sss_indexInMARBL,               &
         atmPressure_indexInMARBL,       &
         atmCO2_indexInMARBL,            &
         atmAltCO2_indexInMARBL,         &
         dustFluxSurface_indexInMARBL,   &
         ironFlux_indexInMARBL,          &
         noxFlux_indexInMARBL,           &
         nhyFlux_indexInMARBL,           &
         co2_flux_indexInMARBL,          &
         co2_flux_alt_co2_indexInMARBL,  &
         totalChl_indexInMARBL,          &
         surface_pH_indexInMARBL,        &
         surface_pH_altCO2_indexInMARBL, &
         interior_pH_indexInMARBL,       &
         interior_pH_altCO2_indexInMARBL

      ! indices for level 1 diagnostics, no depth dependence (nx)
      integer :: &
         ecosys_diag_photoC_TOT_zint_indexInMARBL,  &
         ecosys_diag_photoC_NO3_TOT_zint_indexInMARBL,  &
         ecosys_diag_O2_ZMIN_indexInMARBL,  &
         ecosys_diag_O2_ZMIN_DEPTH_indexInMARBL,  &
         ecosys_diag_calcToSed_indexInMARBL,      &
         ecosys_diag_calcToSed_ALT_CO2_indexInMARBL,      &
         ecosys_diag_pocToSed_indexInMARBL,       &
         ecosys_diag_pfeToSed_indexInMARBL,       &
         ecosys_diag_SedDenitrif_indexInMARBL,    &
         ecosys_diag_Jint_Ctot_indexInMARBL,  &
         ecosys_diag_Jint_Ntot_indexInMARBL,  &
         ecosys_diag_Jint_Ptot_indexInMARBL,  &
         ecosys_diag_Jint_Sitot_indexInMARBL,  &
         ecosys_diag_Jint_Fetot_indexInMARBL
      logical :: &
         ecosys_diag_photoC_TOT_zint_allocatedInMARBL = .false.,  &
         ecosys_diag_photoC_NO3_TOT_zint_allocatedInMARBL = .false.,  &
         ecosys_diag_O2_ZMIN_allocatedInMARBL = .false.,  &
         ecosys_diag_O2_ZMIN_DEPTH_allocatedInMARBL = .false.,  &
         ecosys_diag_calcToSed_allocatedInMARBL = .false.,      &
         ecosys_diag_calcToSed_ALT_CO2_allocatedInMARBL = .false.,      &
         ecosys_diag_pocToSed_allocatedInMARBL = .false.,       &
         ecosys_diag_pfeToSed_allocatedInMARBL = .false.,       &
         ecosys_diag_SedDenitrif_allocatedInMARBL = .false.,    &
         ecosys_diag_Jint_Ctot_allocatedInMARBL = .false.,  &
         ecosys_diag_Jint_Ntot_allocatedInMARBL = .false.,  &
         ecosys_diag_Jint_Ptot_allocatedInMARBL = .false.,  &
         ecosys_diag_Jint_Sitot_allocatedInMARBL = .false.,  &
         ecosys_diag_Jint_Fetot_allocatedInMARBL = .false.

      ! indices for level 1 diagnostics, each autotroph, no depth dependence (nx,autotroph_cnt)
      integer, dimension(autotroph_cnt) :: &
         ecosys_diag_photoC_zint_indexInMARBL,    &
         ecosys_diag_photoC_NO3_zint_indexInMARBL
      logical, dimension(autotroph_cnt) :: &
         ecosys_diag_photoC_zint_allocatedInMARBL = .false.,    &
         ecosys_diag_photoC_NO3_zint_allocatedInMARBL = .false.

      ! indices for level 1 diagnostics, full depth dependence (nx,nz)
      integer :: &
         ecosys_diag_PAR_avg_indexInMARBL,  &
         ecosys_diag_POC_FLUX_IN_indexInMARBL,  &
         ecosys_diag_CaCO3_FLUX_IN_indexInMARBL,  &
         ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_indexInMARBL,  &
         ecosys_diag_graze_auto_TOT_indexInMARBL,  &
         ecosys_diag_zoo_loss_indexInMARBL,  &
         ecosys_diag_zoo_loss_poc_indexInMARBL,  &
         ecosys_diag_zoo_loss_doc_indexInMARBL,  &
         ecosys_diag_graze_zoo_indexInMARBL,  &
         ecosys_diag_graze_zoo_poc_indexInMARBL,  &
         ecosys_diag_graze_zoo_doc_indexInMARBL,  &
         ecosys_diag_graze_zoo_zoo_indexInMARBL,  &
         ecosys_diag_x_graze_zoo_indexInMARBL,  &
         ecosys_diag_photoC_TOT_indexInMARBL,  &
         ecosys_diag_photoC_NO3_TOT_indexInMARBL,  &
         ecosys_diag_NITRIF_indexInMARBL,  &
         ecosys_diag_DENITRIF_indexInMARBL,  &
         ecosys_diag_tot_Nfix_indexInMARBL
      logical :: &
         ecosys_diag_PAR_avg_allocatedInMARBL = .false.,  &
         ecosys_diag_POC_FLUX_IN_allocatedInMARBL = .false.,  &
         ecosys_diag_CaCO3_FLUX_IN_allocatedInMARBL = .false.,  &
         ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_allocatedInMARBL = .false.,  &
         ecosys_diag_graze_auto_TOT_allocatedInMARBL = .false.,  &
         ecosys_diag_zoo_loss_allocatedInMARBL = .false.,  &
         ecosys_diag_zoo_loss_poc_allocatedInMARBL = .false.,  &
         ecosys_diag_zoo_loss_doc_allocatedInMARBL = .false.,  &
         ecosys_diag_graze_zoo_allocatedInMARBL = .false.,  &
         ecosys_diag_graze_zoo_poc_allocatedInMARBL = .false.,  &
         ecosys_diag_graze_zoo_doc_allocatedInMARBL = .false.,  &
         ecosys_diag_graze_zoo_zoo_allocatedInMARBL = .false.,  &
         ecosys_diag_x_graze_zoo_allocatedInMARBL = .false.,  &
         ecosys_diag_photoC_TOT_allocatedInMARBL = .false.,  &
         ecosys_diag_photoC_NO3_TOT_allocatedInMARBL = .false.,  &
         ecosys_diag_NITRIF_allocatedInMARBL = .false.,  &
         ecosys_diag_DENITRIF_allocatedInMARBL = .false.,  &
         ecosys_diag_tot_Nfix_allocatedInMARBL = .false.

      ! indices for level 2 diagnostics, surface fields (nx)
      integer :: &
         ecosys_diag_pistonVel_O2_indexInMARBL,  &
         ecosys_diag_pistonVel_CO2_indexInMARBL,  &
         ecosys_diag_Schmidt_O2_indexInMARBL,  &
         ecosys_diag_Schmidt_CO2_indexInMARBL,  &
         ecosys_diag_O2_saturation_indexInMARBL,  &
         ecosys_diag_xkw_indexInMARBL,  &
         ecosys_diag_u10sqr_indexInMARBL,  &
         ecosys_diag_CO2star_indexInMARBL,  &
         ecosys_diag_dCO2star_indexInMARBL,  &
         ecosys_diag_pCO2surface_indexInMARBL,  &
         ecosys_diag_dpCO2_indexInMARBL,  &
         ecosys_diag_CO2star_ALT_CO2_indexInMARBL,  &
         ecosys_diag_dCO2star_ALT_CO2_indexInMARBL,  &
         ecosys_diag_pCO2surface_ALT_CO2_indexInMARBL,  &
         ecosys_diag_dpCO2_ALT_CO2_indexInMARBL,        &
         ecosys_diag_atmPressure_indexInMARBL,            &
         ecosys_diag_ifrac_indexInMARBL
      logical :: &
         ecosys_diag_pistonVel_O2_allocatedInMARBL = .false.,  &
         ecosys_diag_pistonVel_CO2_allocatedInMARBL = .false.,  &
         ecosys_diag_Schmidt_O2_allocatedInMARBL = .false.,  &
         ecosys_diag_Schmidt_CO2_allocatedInMARBL = .false.,  &
         ecosys_diag_O2_saturation_allocatedInMARBL = .false.,  &
         ecosys_diag_xkw_allocatedInMARBL = .false.,  &
         ecosys_diag_u10sqr_allocatedInMARBL = .false.,  &
         ecosys_diag_ifrac_allocatedInMARBL = .false.,  &
         ecosys_diag_atmPressure_allocatedInMARBL = .false.,  &
         ecosys_diag_CO2star_allocatedInMARBL = .false.,  &
         ecosys_diag_dCO2star_allocatedInMARBL = .false.,  &
         ecosys_diag_pCO2surface_allocatedInMARBL = .false.,  &
         ecosys_diag_dpCO2_allocatedInMARBL = .false.,  &
         ecosys_diag_CO2star_ALT_CO2_allocatedInMARBL = .false.,  &
         ecosys_diag_dCO2star_ALT_CO2_allocatedInMARBL = .false.,  &
         ecosys_diag_pCO2surface_ALT_CO2_allocatedInMARBL = .false.,  &
         ecosys_diag_dpCO2_ALT_CO2_allocatedInMARBL = .false.

      ! indices for level 2 diagnostics, each autotroph, no depth dependence (nx,autotroph_cnt)
      integer, dimension(autotroph_cnt) :: &
         ecosys_diag_N_lim_surf_indexInMARBL,          &
         ecosys_diag_P_lim_surf_indexInMARBL,          &
         ecosys_diag_Fe_lim_surf_indexInMARBL,         &
         ecosys_diag_SiO3_lim_surf_indexInMARBL,       &
         ecosys_diag_light_lim_surf_indexInMARBL
      logical, dimension(autotroph_cnt) :: &
         ecosys_diag_N_lim_surf_allocatedInMARBL = .false.,          &
         ecosys_diag_P_lim_surf_allocatedInMARBL = .false.,          &
         ecosys_diag_Fe_lim_surf_allocatedInMARBL = .false.,         &
         ecosys_diag_SiO3_lim_surf_allocatedInMARBL = .false.,       &
         ecosys_diag_light_lim_surf_allocatedInMARBL = .false.

      ! indices for level 2 diagnostics, full depth dependence (nx,nz)
      integer :: &
         ecosys_diag_O2_PRODUCTION_indexInMARBL,  &
         ecosys_diag_O2_CONSUMPTION_indexInMARBL,  &
         ecosys_diag_AOU_indexInMARBL,          &
         ecosys_diag_pH_3D_indexInMARBL,     &
         ecosys_diag_POC_PROD_indexInMARBL,  &
         ecosys_diag_POC_REMIN_DOCr_indexInMARBL,  &
         ecosys_diag_POC_REMIN_DIC_indexInMARBL,  &
         ecosys_diag_POP_FLUX_IN_indexInMARBL,  &
         ecosys_diag_POP_PROD_indexInMARBL,  &
         ecosys_diag_POP_REMIN_DOPr_indexInMARBL,  &
         ecosys_diag_POP_REMIN_PO4_indexInMARBL,  &
         ecosys_diag_PON_REMIN_DONr_indexInMARBL, &
         ecosys_diag_PON_REMIN_NH4_indexInMARBL
      logical :: &
         ecosys_diag_O2_PRODUCTION_allocatedInMARBL = .false.,  &
         ecosys_diag_O2_CONSUMPTION_allocatedInMARBL = .false.,  &
         ecosys_diag_AOU_allocatedInMARBL = .false.,          &
         ecosys_diag_pH_3D_allocatedInMARBL = .false.,     &
         ecosys_diag_POC_PROD_allocatedInMARBL = .false.,  &
         ecosys_diag_POC_REMIN_DOCr_allocatedInMARBL = .false.,  &
         ecosys_diag_POC_REMIN_DIC_allocatedInMARBL = .false.,  &
         ecosys_diag_POP_FLUX_IN_allocatedInMARBL = .false.,  &
         ecosys_diag_POP_PROD_allocatedInMARBL = .false.,  &
         ecosys_diag_POP_REMIN_DOPr_allocatedInMARBL = .false.,  &
         ecosys_diag_POP_REMIN_PO4_allocatedInMARBL = .false.,  &
         ecosys_diag_PON_REMIN_DONr_allocatedInMARBL = .false., &
         ecosys_diag_PON_REMIN_NH4_allocatedInMARBL = .false.

      ! indices for level 2 diagnostics, each autotroph, full depth dependence (nx,nz,autotroph_cnt)
      integer, dimension(autotroph_cnt) :: &
         ecosys_diag_Qp_indexInMARBL,             &
         ecosys_diag_photoC_indexInMARBL,         &
         ecosys_diag_photoC_NO3_indexInMARBL,     &
         ecosys_diag_photoFe_indexInMARBL,        &
         ecosys_diag_photoNO3_indexInMARBL,       &
         ecosys_diag_photoNH4_indexInMARBL,       &
         ecosys_diag_DOP_uptake_indexInMARBL,     &
         ecosys_diag_PO4_uptake_indexInMARBL,     &
         ecosys_diag_auto_graze_indexInMARBL,     &
         ecosys_diag_auto_graze_poc_indexInMARBL,     &
         ecosys_diag_auto_graze_doc_indexInMARBL,     &
         ecosys_diag_auto_graze_zoo_indexInMARBL,     &
         ecosys_diag_auto_loss_indexInMARBL,          &
         ecosys_diag_auto_loss_poc_indexInMARBL,      &
         ecosys_diag_auto_loss_doc_indexInMARBL,      &
         ecosys_diag_auto_agg_indexInMARBL,           &
         ecosys_diag_Nfix_indexInMARBL
      logical, dimension(autotroph_cnt) :: &
         ecosys_diag_Qp_allocatedInMARBL = .false.,             &
         ecosys_diag_photoC_allocatedInMARBL = .false.,         &
         ecosys_diag_photoC_NO3_allocatedInMARBL = .false.,     &
         ecosys_diag_photoFe_allocatedInMARBL = .false.,        &
         ecosys_diag_photoNO3_allocatedInMARBL = .false.,       &
         ecosys_diag_photoNH4_allocatedInMARBL = .false.,       &
         ecosys_diag_DOP_uptake_allocatedInMARBL = .false.,     &
         ecosys_diag_PO4_uptake_allocatedInMARBL = .false.,     &
         ecosys_diag_auto_graze_allocatedInMARBL = .false.,     &
         ecosys_diag_auto_graze_poc_allocatedInMARBL = .false.,     &
         ecosys_diag_auto_graze_doc_allocatedInMARBL = .false.,     &
         ecosys_diag_auto_graze_zoo_allocatedInMARBL = .false.,     &
         ecosys_diag_auto_loss_allocatedInMARBL = .false.,          &
         ecosys_diag_auto_loss_poc_allocatedInMARBL = .false.,      &
         ecosys_diag_auto_loss_doc_allocatedInMARBL = .false.,      &
         ecosys_diag_auto_agg_allocatedInMARBL = .false.,           &
         ecosys_diag_Nfix_allocatedInMARBL = .false.

      ! indices for level 3 diagnostics, no depth dependence (nx)
      integer :: &
         ecosys_diag_ponToSed_indexInMARBL,       &
         ecosys_diag_popToSed_indexInMARBL,       &
         ecosys_diag_bsiToSed_indexInMARBL,       &
         ecosys_diag_dustToSed_indexInMARBL,      &
         ecosys_diag_OtherRemin_indexInMARBL
      logical :: &
         ecosys_diag_ponToSed_allocatedInMARBL = .false.,       &
         ecosys_diag_popToSed_allocatedInMARBL = .false.,       &
         ecosys_diag_bsiToSed_allocatedInMARBL = .false.,       &
         ecosys_diag_dustToSed_allocatedInMARBL = .false.,      &
         ecosys_diag_OtherRemin_allocatedInMARBL = .false.

      ! indices for level 3 diagnostics, full depth dependence (nx,nz)
      integer :: &
         ecosys_diag_SiO2_FLUX_IN_indexInMARBL,  &
         ecosys_diag_SiO2_PROD_indexInMARBL,  &
         ecosys_diag_SiO2_REMIN_indexInMARBL,  &
         ecosys_diag_dust_FLUX_IN_indexInMARBL,  &
         ecosys_diag_dust_REMIN_indexInMARBL,  &
         ecosys_diag_P_iron_FLUX_IN_indexInMARBL,  &
         ecosys_diag_P_iron_PROD_indexInMARBL,  &
         ecosys_diag_P_iron_REMIN_indexInMARBL,  &
         ecosys_diag_DOC_prod_indexInMARBL,  &
         ecosys_diag_DOC_remin_indexInMARBL,  &
         ecosys_diag_DOCr_remin_indexInMARBL,  &
         ecosys_diag_DON_prod_indexInMARBL,  &
         ecosys_diag_DON_remin_indexInMARBL,  &
         ecosys_diag_DOP_prod_indexInMARBL,  &
         ecosys_diag_DOP_remin_indexInMARBL,  &
         ecosys_diag_Fe_scavenge_indexInMARBL,  &
         ecosys_diag_Fe_scavenge_rate_indexInMARBL,  &
         ecosys_diag_Lig_prod_indexInMARBL,  &
         ecosys_diag_Lig_loss_indexInMARBL,  &
         ecosys_diag_Lig_scavenge_indexInMARBL,  &
         ecosys_diag_Fefree_indexInMARBL,  &
         ecosys_diag_Lig_photochem_indexInMARBL,  &
         ecosys_diag_Lig_deg_indexInMARBL,  &
         ecosys_diag_DONr_remin_indexInMARBL,     &
         ecosys_diag_DOPr_remin_indexInMARBL,     &
         ecosys_diag_tot_bSi_form_indexInMARBL
      logical :: &
         ecosys_diag_SiO2_FLUX_IN_allocatedInMARBL = .false.,  &
         ecosys_diag_SiO2_PROD_allocatedInMARBL = .false.,  &
         ecosys_diag_SiO2_REMIN_allocatedInMARBL = .false.,  &
         ecosys_diag_dust_FLUX_IN_allocatedInMARBL = .false.,  &
         ecosys_diag_dust_REMIN_allocatedInMARBL = .false.,  &
         ecosys_diag_P_iron_FLUX_IN_allocatedInMARBL = .false.,  &
         ecosys_diag_P_iron_PROD_allocatedInMARBL = .false.,  &
         ecosys_diag_P_iron_REMIN_allocatedInMARBL = .false.,  &
         ecosys_diag_DOC_prod_allocatedInMARBL = .false.,  &
         ecosys_diag_DOC_remin_allocatedInMARBL = .false.,  &
         ecosys_diag_DOCr_remin_allocatedInMARBL = .false.,  &
         ecosys_diag_DON_prod_allocatedInMARBL = .false.,  &
         ecosys_diag_DON_remin_allocatedInMARBL = .false.,  &
         ecosys_diag_DOP_prod_allocatedInMARBL = .false.,  &
         ecosys_diag_DOP_remin_allocatedInMARBL = .false.,  &
         ecosys_diag_Fe_scavenge_allocatedInMARBL = .false.,  &
         ecosys_diag_Fe_scavenge_rate_allocatedInMARBL = .false.,  &
         ecosys_diag_Lig_prod_allocatedInMARBL = .false.,  &
         ecosys_diag_Lig_loss_allocatedInMARBL = .false.,  &
         ecosys_diag_Lig_scavenge_allocatedInMARBL = .false.,  &
         ecosys_diag_Fefree_allocatedInMARBL = .false.,  &
         ecosys_diag_Lig_photochem_allocatedInMARBL = .false.,  &
         ecosys_diag_Lig_deg_allocatedInMARBL = .false.,  &
         ecosys_diag_DONr_remin_allocatedInMARBL = .false.,     &
         ecosys_diag_DOPr_remin_allocatedInMARBL = .false.,     &
         ecosys_diag_tot_bSi_form_allocatedInMARBL = .false.

      ! indices for level 3 diagnostics, each autotroph, full depth dependence (nx,nz,autotroph_cnt)
      integer, dimension(autotroph_cnt) :: &
         ecosys_diag_bSi_form_indexInMARBL
      logical, dimension(autotroph_cnt) :: &
         ecosys_diag_bSi_form_allocatedInMARBL = .false.

      ! indices for level 4 diagnostics, no depth dependence (nx)
      integer :: &
         ecosys_diag_tot_CaCO3_form_zint_indexInMARBL, &
         ecosys_diag_zsatcalc_indexInMARBL,            &
         ecosys_diag_zsatarag_indexInMARBL
      logical :: &
         ecosys_diag_tot_CaCO3_form_zint_allocatedInMARBL, &
         ecosys_diag_zsatcalc_allocatedInMARBL,            &
         ecosys_diag_zsatarag_allocatedInMARBL

      ! indices for level 4 diagnostics, each autotroph, no depth dependence (nx,autotroph_cnt)
      integer, dimension(autotroph_cnt) :: &
         ecosys_diag_CaCO3_form_zint_indexInMARBL
      logical, dimension(autotroph_cnt) :: &
         ecosys_diag_CaCO3_form_zint_allocatedInMARBL = .false.

      ! indices for level 4 diagnostics, full depth dependence (nx,nz)
      integer :: &
         ecosys_diag_CaCO3_PROD_indexInMARBL,     &
         ecosys_diag_CaCO3_REMIN_indexInMARBL,    &
         ecosys_diag_CaCO3_ALT_CO2_PROD_indexInMARBL,     &
         ecosys_diag_CaCO3_ALT_CO2_REMIN_indexInMARBL,    &
         ecosys_diag_tot_CaCO3_form_indexInMARBL, &
         ecosys_diag_CO3_indexInMARBL,            &
         ecosys_diag_HCO3_indexInMARBL,           &
         ecosys_diag_H2CO3_indexInMARBL,          &
         ecosys_diag_CO3_ALT_CO2_indexInMARBL,    &
         ecosys_diag_HCO3_ALT_CO2_indexInMARBL,   &
         ecosys_diag_H2CO3_ALT_CO2_indexInMARBL,  &
         ecosys_diag_pH_3D_ALT_CO2_indexInMARBL,  &
         ecosys_diag_co3_sat_calc_indexInMARBL,   &
         ecosys_diag_co3_sat_arag_indexInMARBL
      logical :: &
         ecosys_diag_CaCO3_PROD_allocatedInMARBL = .false.,     &
         ecosys_diag_CaCO3_REMIN_allocatedInMARBL = .false.,    &
         ecosys_diag_CaCO3_ALT_CO2_PROD_allocatedInMARBL = .false.,     &
         ecosys_diag_CaCO3_ALT_CO2_REMIN_allocatedInMARBL = .false.,    &
         ecosys_diag_tot_CaCO3_form_allocatedInMARBL = .false., &
         ecosys_diag_CO3_allocatedInMARBL = .false.,            &
         ecosys_diag_HCO3_allocatedInMARBL = .false.,           &
         ecosys_diag_H2CO3_allocatedInMARBL = .false.,          &
         ecosys_diag_CO3_ALT_CO2_allocatedInMARBL = .false.,    &
         ecosys_diag_HCO3_ALT_CO2_allocatedInMARBL = .false.,   &
         ecosys_diag_H2CO3_ALT_CO2_allocatedInMARBL = .false.,  &
         ecosys_diag_pH_3D_ALT_CO2_allocatedInMARBL = .false.,  &
         ecosys_diag_co3_sat_calc_allocatedInMARBL = .false.,   &
         ecosys_diag_co3_sat_arag_allocatedInMARBL = .false.

      ! indices for level 4 diagnostics, each autotroph, full depth dependence (nx,nz,autotroph_cnt)
      integer, dimension(autotroph_cnt) :: &
         ecosys_diag_CaCO3_form_indexInMARBL
      logical, dimension(autotroph_cnt) :: &
         ecosys_diag_CaCO3_form_allocatedInMARBL = .false.

      ! indices for level 5 diagnostics, full depth dependence (nx,nz)
      integer :: &
         ecosys_diag_PO4_RESTORE_indexInMARBL,  &
         ecosys_diag_NO3_RESTORE_indexInMARBL,  &
         ecosys_diag_SiO3_RESTORE_indexInMARBL
      logical :: &
         ecosys_diag_PO4_RESTORE_allocatedInMARBL = .false.,  &
         ecosys_diag_NO3_RESTORE_allocatedInMARBL = .false.,  &
         ecosys_diag_SiO3_RESTORE_allocatedInMARBL = .false.

      real (kind=RKIND) :: ecosysMissingValue = -999.

      real (kind=RKIND) :: &
         PascalsToAtmospheres = 1.0_RKIND/101.325e+3_RKIND

      ! flux routine wants CGS units
      ! then convert back to MKS after
      real (kind=RKIND) :: convertLengthMKStoCGS = 100.0_RKIND
      real (kind=RKIND) :: convertLengthCGStoMKS = 0.01_RKIND
      real (kind=RKIND) :: convertLengthSquaredMKStoCGS = 1.0e4_RKIND
      real (kind=RKIND) :: convertLengthSquaredCGStoMKS = 1.0e-4_RKIND

      logical :: verbose = .false.
      character(strKIND) :: logMessage

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
      latCell, minLevelCell, maxLevelCell, nVertLevels, layerThickness, zMid, indexTemperature, indexSalinity,  &
      ecosysTracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! one dimensional arrays
      integer, dimension(:), intent(in) :: &
         minLevelCell, maxLevelCell
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

      !
      ! level 1 diagnostics
      !

      real (kind=RKIND), dimension(:), pointer :: &
         ecosys_diag_photoC_TOT_zint,  &
         ecosys_diag_photoC_NO3_TOT_zint,  &
         ecosys_diag_O2_ZMIN,  &
         ecosys_diag_O2_ZMIN_DEPTH,  &
         ecosys_diag_Chl_TOT_zint_100m,  &
         ecosys_diag_calcToSed,      &
         ecosys_diag_calcToSed_ALT_CO2,      &
         ecosys_diag_pocToSed,       &
         ecosys_diag_pfeToSed,       &
         ecosys_diag_SedDenitrif,    &
         ecosys_diag_Jint_Ctot,  &
         ecosys_diag_Jint_Ntot,  &
         ecosys_diag_Jint_Ptot,  &
         ecosys_diag_Jint_Sitot,  &
         ecosys_diag_Jint_Fetot

      ! one dimensional arrays for each autotroph (nx,autotroph_cnt)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_photoC_zint,    &
         ecosys_diag_photoC_NO3_zint

      ! two dimensional arrays (nx,nz)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_PAR_avg,  &
         ecosys_diag_POC_FLUX_IN,  &
         ecosys_diag_CaCO3_FLUX_IN,  &
         ecosys_diag_CaCO3_ALT_CO2_FLUX_IN,  &
         ecosys_diag_graze_auto_TOT,  &
         ecosys_diag_zoo_loss,  &
         ecosys_diag_zoo_loss_poc,  &
         ecosys_diag_zoo_loss_doc,  &
         ecosys_diag_graze_zoo,  &
         ecosys_diag_graze_zoo_poc,  &
         ecosys_diag_graze_zoo_doc,  &
         ecosys_diag_graze_zoo_zoo,  &
         ecosys_diag_x_graze_zoo,  &
         ecosys_diag_photoC_TOT,  &
         ecosys_diag_photoC_NO3_TOT,  &
         ecosys_diag_NITRIF,  &
         ecosys_diag_DENITRIF,  &
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
         ecosys_diag_POC_REMIN_DOCr,  &
         ecosys_diag_POC_REMIN_DIC,  &
         ecosys_diag_POP_FLUX_IN,  &
         ecosys_diag_POP_PROD,  &
         ecosys_diag_POP_REMIN_DOPr,  &
         ecosys_diag_POP_REMIN_PO4,  &
         ecosys_diag_PON_REMIN_DONr, &
         ecosys_diag_PON_REMIN_NH4

      ! one dimensional arrays for each autotroph (nx,autotroph_cnt)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_N_lim_surf,          &
         ecosys_diag_P_lim_surf,          &
         ecosys_diag_Fe_lim_surf,         &
         ecosys_diag_SiO3_lim_surf,       &
         ecosys_diag_light_lim_surf

      ! two dimensional arrays for each autotroph (nx,nz,autotroph_cnt)
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         ecosys_diag_Qp,             &
         ecosys_diag_photoC,         &
         ecosys_diag_photoC_NO3,     &
         ecosys_diag_photoFe,        &
         ecosys_diag_photoNO3,       &
         ecosys_diag_photoNH4,       &
         ecosys_diag_DOP_uptake,     &
         ecosys_diag_PO4_uptake,     &
         ecosys_diag_auto_graze,     &
         ecosys_diag_auto_graze_poc,     &
         ecosys_diag_auto_graze_doc,     &
         ecosys_diag_auto_graze_zoo,     &
         ecosys_diag_auto_loss,          &
         ecosys_diag_auto_loss_poc,      &
         ecosys_diag_auto_loss_doc,      &
         ecosys_diag_auto_agg,           &
         ecosys_diag_Nfix

      !
      ! level 3 diagnostics
      !

      ! one dimensional arrays (nx)
      real (kind=RKIND), dimension(:), pointer :: &
         ecosys_diag_ponToSed,       &
         ecosys_diag_popToSed,       &
         ecosys_diag_bsiToSed,       &
         ecosys_diag_dustToSed,      &
         ecosys_diag_OtherRemin

      ! two dimensional arrays (nx,nz)
      real (kind=RKIND), dimension(:,:), pointer :: &
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
         ecosys_diag_DOCr_remin,  &
         ecosys_diag_DON_prod,  &
         ecosys_diag_DON_remin,  &
         ecosys_diag_DOP_prod,  &
         ecosys_diag_DOP_remin,  &
         ecosys_diag_Fe_scavenge,  &
         ecosys_diag_Fe_scavenge_rate,  &
         ecosys_diag_Lig_prod,       &
         ecosys_diag_Lig_loss,       &
         ecosys_diag_Lig_scavenge,   &
         ecosys_diag_Fefree,         &
         ecosys_diag_Lig_photochem,  &
         ecosys_diag_Lig_deg,        &
         ecosys_diag_DONr_remin,     &
         ecosys_diag_DOPr_remin,     &
         ecosys_diag_tot_bSi_form

      ! two dimensional arrays for each autotroph (nx,nz,autotroph_cnt)
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         ecosys_diag_bSi_form

      !
      ! level 4 diagnostics
      !

      ! one dimensional arrays (nx)
      real (kind=RKIND), dimension(:), pointer :: &
         ecosys_diag_tot_CaCO3_form_zint, &
         ecosys_diag_zsatcalc,            &
         ecosys_diag_zsatarag

      ! one dimensional arrays for each autotroph (nx,autotroph_cnt)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_CaCO3_form_zint  ! diag array for CaCO3 formation vertical integral

      ! two dimensional arrays (nx,nz)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_CaCO3_PROD,     &
         ecosys_diag_CaCO3_REMIN,    &
         ecosys_diag_CaCO3_ALT_CO2_PROD,     &
         ecosys_diag_CaCO3_ALT_CO2_REMIN,    &
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
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         ecosys_diag_CaCO3_form

      !
      !
      ! level 5 diagnostics
      !

      ! two dimensional arrays (nx,nz)
      real (kind=RKIND), dimension(:,:), pointer :: &
         ecosys_diag_PO4_RESTORE,  &
         ecosys_diag_NO3_RESTORE,  &
         ecosys_diag_SiO3_RESTORE

      real (kind=RKIND) :: zTop, zBot

      integer :: iCell, iLevel, numColumns, column, autotroph, n

      err = 0

      call mpas_timer_start("ecosys source-sink")

      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_3D', PH_PREV_3D)
      call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_ALT_CO2_3D', PH_PREV_ALT_CO2_3D)
      call mpas_pool_get_array(ecosysAuxiliary, 'FESEDFLUX', FESEDFLUX)
      call mpas_pool_get_array(ecosysAuxiliary, 'dust_FLUX_IN', dust_FLUX_IN)
      call mpas_pool_get_array(ecosysAuxiliary, 'PAR_surface', PAR_surface)

      call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)

      if (config_ecosysTracers_diagnostic_fields_level1) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel1', ecosysDiagFieldsLevel1)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_TOT_zint', ecosys_diag_photoC_TOT_zint)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_NO3_TOT_zint', ecosys_diag_photoC_NO3_TOT_zint)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_O2_ZMIN', ecosys_diag_O2_ZMIN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_O2_ZMIN_DEPTH', ecosys_diag_O2_ZMIN_DEPTH)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Ctot', ecosys_diag_Jint_Ctot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Ntot', ecosys_diag_Jint_Ntot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Ptot', ecosys_diag_Jint_Ptot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Sitot', ecosys_diag_Jint_Sitot)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_Jint_Fetot', ecosys_diag_Jint_Fetot)

         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_zint', ecosys_diag_photoC_zint)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_NO3_zint', ecosys_diag_photoC_NO3_zint)

         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_PAR_avg', ecosys_diag_PAR_avg)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_POC_FLUX_IN', ecosys_diag_POC_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_CaCO3_FLUX_IN', ecosys_diag_CaCO3_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_CaCO3_ALT_CO2_FLUX_IN',  &
            ecosys_diag_CaCO3_ALT_CO2_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_graze_auto_TOT', ecosys_diag_graze_auto_TOT)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_zoo_loss', ecosys_diag_zoo_loss)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_zoo_loss_poc', ecosys_diag_zoo_loss_poc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_zoo_loss_doc', ecosys_diag_zoo_loss_doc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_graze_zoo', ecosys_diag_graze_zoo)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_graze_zoo_poc', ecosys_diag_graze_zoo_poc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_graze_zoo_doc', ecosys_diag_graze_zoo_doc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_graze_zoo_zoo', ecosys_diag_graze_zoo_zoo)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_x_graze_zoo', ecosys_diag_x_graze_zoo)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_TOT', ecosys_diag_photoC_TOT)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_photoC_NO3_TOT', ecosys_diag_photoC_NO3_TOT)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_NITRIF', ecosys_diag_NITRIF)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_DENITRIF', ecosys_diag_DENITRIF)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_calcToSed', ecosys_diag_calcToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel1, 'ecosys_diag_calcToSed_ALT_CO2',  &
            ecosys_diag_calcToSed_ALT_CO2)
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
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POC_REMIN_DOCr', ecosys_diag_POC_REMIN_DOCr)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POC_REMIN_DIC', ecosys_diag_POC_REMIN_DIC)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POP_FLUX_IN', ecosys_diag_POP_FLUX_IN)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POP_PROD', ecosys_diag_POP_PROD)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POP_REMIN_DOPr', ecosys_diag_POP_REMIN_DOPr)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_POP_REMIN_PO4', ecosys_diag_POP_REMIN_PO4)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_PON_REMIN_DONr', ecosys_diag_PON_REMIN_DONr)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_PON_REMIN_NH4', ecosys_diag_PON_REMIN_NH4)

         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_N_lim_surf', ecosys_diag_N_lim_surf)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_P_lim_surf', ecosys_diag_P_lim_surf)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Fe_lim_surf', ecosys_diag_Fe_lim_surf)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_SiO3_lim_surf', ecosys_diag_SiO3_lim_surf)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_light_lim_surf', ecosys_diag_light_lim_surf)

         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Qp', ecosys_diag_Qp)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoC', ecosys_diag_photoC)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoC_NO3', ecosys_diag_photoC_NO3)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoFe', ecosys_diag_photoFe)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoNO3', ecosys_diag_photoNO3)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_photoNH4', ecosys_diag_photoNH4)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_DOP_uptake', ecosys_diag_DOP_uptake)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_PO4_uptake', ecosys_diag_PO4_uptake)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_graze', ecosys_diag_auto_graze)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_graze_poc', ecosys_diag_auto_graze_poc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_graze_doc', ecosys_diag_auto_graze_doc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_graze_zoo', ecosys_diag_auto_graze_zoo)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_loss', ecosys_diag_auto_loss)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_loss_poc', ecosys_diag_auto_loss_poc)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_auto_loss_doc', ecosys_diag_auto_loss_doc)
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
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOCr_remin', ecosys_diag_DOCr_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DON_prod', ecosys_diag_DON_prod)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DON_remin', ecosys_diag_DON_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DONr_remin', ecosys_diag_DONr_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOP_prod', ecosys_diag_DOP_prod)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOP_remin', ecosys_diag_DOP_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_DOPr_remin', ecosys_diag_DOPr_remin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Fe_scavenge', ecosys_diag_Fe_scavenge)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Fe_scavenge_rate', ecosys_diag_Fe_scavenge_rate)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_ponToSed', ecosys_diag_ponToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_popToSed', ecosys_diag_popToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_bsiToSed', ecosys_diag_bsiToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_dustToSed', ecosys_diag_dustToSed)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_OtherRemin', ecosys_diag_OtherRemin)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Lig_prod', ecosys_diag_Lig_prod)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Lig_loss', ecosys_diag_Lig_loss)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Lig_scavenge', ecosys_diag_Lig_scavenge)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Fefree', ecosys_diag_Fefree)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Lig_photochem', ecosys_diag_Lig_photochem)
         call mpas_pool_get_array(ecosysDiagFieldsLevel3, 'ecosys_diag_Lig_deg', ecosys_diag_Lig_deg)

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
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CaCO3_ALT_CO2_PROD',  &
            ecosys_diag_CaCO3_ALT_CO2_PROD)
         call mpas_pool_get_array(ecosysDiagFieldsLevel4, 'ecosys_diag_CaCO3_ALT_CO2_REMIN',  &
            ecosys_diag_CaCO3_ALT_CO2_REMIN)
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

         marbl_instances(1)%domain%kmt = maxLevelCell(iCell)
         marbl_instances(1)%interior_tendency_forcings(dustFlux_indexInMARBL)%field_0d(1)  = dust_FLUX_IN(iCell)
! NOTE shortwave has a second (subcolumn) index, so use field1d instead of field0d
!  also we need to define the PAR column fraction (always 1 for now)
         marbl_instances(1)%interior_tendency_forcings(shortwave_indexInMARBL)%field_1d(1,1) = shortWaveHeatFlux(iCell)
         if (PAR_col_frac_indexInMARBL > 0) &
         marbl_instances(1)%interior_tendency_forcings(PAR_col_frac_indexInMARBL)%field_1d(1,1) = 1.0_RKIND

! make sure full column is free of garbage values--even below the bottom
         marbl_instances(1)%interior_tendency_forcings(temperature_indexInMARBL)%field_1d(1,:) = 0.0_RKIND
         marbl_instances(1)%interior_tendency_forcings(salinity_indexInMARBL)%field_1d(1,:) = 0.0_RKIND
         marbl_instances(1)%interior_tendency_forcings(pressure_indexInMARBL)%field_1d(1,:) = 1.0_RKIND
         marbl_instances(1)%interior_tendency_forcings(fesedFlux_indexInMARBL)%field_1d(1,:) = 0.0_RKIND
         marbl_instances(1)%domain%delta_z(:) = 100.0_RKIND
         marbl_instances(1)%domain%zw(:) = 100.0_RKIND
         marbl_instances(1)%interior_tendency_saved_state%state(interior_pH_indexInMARBL)%field_3d(:,1) = 8.0_RKIND
         marbl_instances(1)%interior_tendency_saved_state%state(interior_pH_altCO2_indexInMARBL)%field_3d(:,1) = 8.0_RKIND
         do n = 1, ecosys_tracer_cnt
            marbl_instances(1)%tracers(n, :) = 0.0_RKIND
         end do

! NOTE we will likely have to modify logic for under ice shelves
! for now ASSUME NO ICE SHELF CAVITIES
         zTop = 0.0_RKIND
         do iLevel=minLevelCell(iCell),maxLevelCell(iCell)
            marbl_instances(1)%interior_tendency_forcings(temperature_indexInMARBL)%field_1d(1,iLevel) = &
               activeTracers(indexTemperature,iLevel,iCell)
            marbl_instances(1)%interior_tendency_forcings(salinity_indexInMARBL)%field_1d(1,iLevel) = &
               activeTracers(indexSalinity,iLevel,iCell)
            marbl_instances(1)%interior_tendency_forcings(pressure_indexInMARBL)%field_1d(1,iLevel) = &
               pressure(iLevel,iCell)*PascalsToAtmospheres
            marbl_instances(1)%interior_tendency_forcings(fesedFlux_indexInMARBL)%field_1d(1,iLevel) = &
               FESEDFLUX(iLevel,iCell)

            marbl_instances(1)%domain%zt(iLevel) = -1.0_RKIND*zMid(iLevel,iCell)*convertLengthMKStoCGS

            marbl_instances(1)%domain%delta_z(iLevel) = layerThickness(iLevel,iCell)*convertLengthMKStoCGS
            zBot = zTop + layerThickness(iLevel,iCell)
            marbl_instances(1)%domain%zw(iLevel) = zBot*convertLengthMKStoCGS
            zTop = zBot

! load pH fields into saved state
            marbl_instances(1)%interior_tendency_saved_state%state(interior_pH_indexInMARBL)%field_3d(iLevel,1) = &
               PH_PREV_3D(iLevel,iCell)
            marbl_instances(1)%interior_tendency_saved_state%state(interior_pH_altCO2_indexInMARBL)%field_3d(iLevel,1) = &
               PH_PREV_ALT_CO2_3D(iLevel,iCell)

! load BGC fields
            do n = 1, ecosys_tracer_cnt
               marbl_instances(1)%tracers(n, iLevel) = ecosysTracers(ecosysIndexMARBLtoMPASO(n),iLevel,iCell)
            end do

         enddo  !  iLevel

! compute source-sink terms
         call marbl_instances(1)%interior_tendency_compute()

         if (marbl_instances(1)%StatusLog%labort_marbl) &
           call mpas_log_write("WARNING: problem in MARBL interior_tendency_compute", MPAS_LOG_CRIT)
        call marbl_instances(1)%StatusLog%erase()

         do iLevel=minLevelCell(iCell),maxLevelCell(iCell)

! load tendencies
            do n = 1, ecosys_tracer_cnt
               ecosysTracersTend(ecosysIndexMARBLtoMPASO(n),iLevel,iCell) =  &
                  ecosysTracersTend(ecosysIndexMARBLtoMPASO(n),iLevel,iCell) &
                  + marbl_instances(1)%interior_tendencies(n,iLevel)*layerThickness(iLevel,iCell)
            end do

! load pH
            PH_PREV_3D(iLevel,iCell) = &
               marbl_instances(1)%interior_tendency_saved_state%state(interior_pH_indexInMARBL)%field_3d(iLevel,1)
            PH_PREV_ALT_CO2_3D(iLevel,iCell) = &
               marbl_instances(1)%interior_tendency_saved_state%state(interior_pH_altCO2_indexInMARBL)%field_3d(iLevel,1)

            !
            ! level 1 diagnostics, full depth
            !

            if (config_ecosysTracers_diagnostic_fields_level1) then

              ecosys_diag_PAR_avg(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_PAR_avg_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_PAR_avg_allocatedInMARBL)

              ecosys_diag_POC_FLUX_IN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POC_FLUX_IN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POC_FLUX_IN_allocatedInMARBL)

              ecosys_diag_CaCO3_FLUX_IN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CaCO3_FLUX_IN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CaCO3_FLUX_IN_allocatedInMARBL)

              ecosys_diag_CaCO3_ALT_CO2_FLUX_IN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_allocatedInMARBL)

              ecosys_diag_graze_auto_TOT(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_graze_auto_TOT_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_graze_auto_TOT_allocatedInMARBL)

              ecosys_diag_zoo_loss(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_zoo_loss_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_zoo_loss_allocatedInMARBL)

              ecosys_diag_zoo_loss_poc(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_zoo_loss_poc_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_zoo_loss_poc_allocatedInMARBL)

              ecosys_diag_zoo_loss_doc(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_zoo_loss_doc_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_zoo_loss_doc_allocatedInMARBL)

              ecosys_diag_graze_zoo(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_graze_zoo_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_graze_zoo_allocatedInMARBL)

              ecosys_diag_graze_zoo_poc(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_graze_zoo_poc_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_graze_zoo_poc_allocatedInMARBL)

              ecosys_diag_graze_zoo_doc(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_graze_zoo_doc_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_graze_zoo_doc_allocatedInMARBL)

              ecosys_diag_graze_zoo_zoo(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_graze_zoo_zoo_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_graze_zoo_zoo_allocatedInMARBL)

              ecosys_diag_x_graze_zoo(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_x_graze_zoo_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_x_graze_zoo_allocatedInMARBL)

              ecosys_diag_photoC_TOT(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_photoC_TOT_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_photoC_TOT_allocatedInMARBL)

              ecosys_diag_photoC_NO3_TOT(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_photoC_NO3_TOT_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_photoC_NO3_TOT_allocatedInMARBL)

              ecosys_diag_NITRIF(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_NITRIF_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_NITRIF_allocatedInMARBL)

              ecosys_diag_DENITRIF(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DENITRIF_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DENITRIF_allocatedInMARBL)

              ecosys_diag_tot_Nfix(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_tot_Nfix_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_tot_Nfix_allocatedInMARBL)

            endif  !  config_ecosysTracers_diagnostic_fields_level1

            !
            ! level 2 diagnostics, full depth
            !
            if (config_ecosysTracers_diagnostic_fields_level2) then

              ecosys_diag_O2_PRODUCTION(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_O2_PRODUCTION_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_O2_PRODUCTION_allocatedInMARBL)

              ecosys_diag_O2_CONSUMPTION(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_O2_CONSUMPTION_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_O2_CONSUMPTION_allocatedInMARBL)

              ecosys_diag_AOU(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_AOU_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_AOU_allocatedInMARBL)

              ecosys_diag_pH_3D(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_pH_3D_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_pH_3D_allocatedInMARBL)

              ecosys_diag_POC_PROD(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POC_PROD_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POC_PROD_allocatedInMARBL)

              ecosys_diag_POC_REMIN_DOCr(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POC_REMIN_DOCr_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POC_REMIN_DOCr_allocatedInMARBL)

              ecosys_diag_POC_REMIN_DIC(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POC_REMIN_DIC_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POC_REMIN_DIC_allocatedInMARBL)

              ecosys_diag_POP_FLUX_IN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POP_FLUX_IN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POP_FLUX_IN_allocatedInMARBL)

              ecosys_diag_POP_PROD(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POP_PROD_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POP_PROD_allocatedInMARBL)

              ecosys_diag_POP_REMIN_DOPr(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POP_REMIN_DOPr_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POP_REMIN_DOPr_allocatedInMARBL)

              ecosys_diag_POP_REMIN_PO4(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_POP_REMIN_PO4_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_POP_REMIN_PO4_allocatedInMARBL)

              ecosys_diag_PON_REMIN_DONr(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_PON_REMIN_DONr_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_PON_REMIN_DONr_allocatedInMARBL)

              ecosys_diag_PON_REMIN_NH4(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_PON_REMIN_NH4_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_PON_REMIN_NH4_allocatedInMARBL)

              do autotroph = 1, autotroph_cnt

                ecosys_diag_Qp(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_Qp_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_Qp_allocatedInMARBL(autotroph))

                ecosys_diag_photoC(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_photoC_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_photoC_allocatedInMARBL(autotroph))

                ecosys_diag_photoC_NO3(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_photoC_NO3_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_photoC_NO3_allocatedInMARBL(autotroph))

                ecosys_diag_photoFe(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_photoFe_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_photoFe_allocatedInMARBL(autotroph))

                ecosys_diag_photoNO3(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_photoNO3_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_photoNO3_allocatedInMARBL(autotroph))

                ecosys_diag_photoNH4(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_photoNH4_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_photoNH4_allocatedInMARBL(autotroph))

                ecosys_diag_DOP_uptake(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_DOP_uptake_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_DOP_uptake_allocatedInMARBL(autotroph))

                ecosys_diag_PO4_uptake(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_PO4_uptake_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_PO4_uptake_allocatedInMARBL(autotroph))

                ecosys_diag_auto_graze(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_graze_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_graze_allocatedInMARBL(autotroph))

                ecosys_diag_auto_graze_poc(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_graze_poc_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_graze_poc_allocatedInMARBL(autotroph))

                ecosys_diag_auto_graze_doc(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_graze_doc_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_graze_doc_allocatedInMARBL(autotroph))

                ecosys_diag_auto_graze_zoo(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_graze_zoo_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_graze_zoo_allocatedInMARBL(autotroph))

                ecosys_diag_auto_loss(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_loss_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_loss_allocatedInMARBL(autotroph))

                ecosys_diag_auto_loss_poc(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_loss_poc_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_loss_poc_allocatedInMARBL(autotroph))

                ecosys_diag_auto_loss_doc(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_loss_doc_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_loss_doc_allocatedInMARBL(autotroph))

                ecosys_diag_auto_agg(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_auto_agg_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_auto_agg_allocatedInMARBL(autotroph))

! Nfix only defined for N-fixers
                if (ecosys_diag_Nfix_indexInMARBL(autotroph) > 0) then
                ecosys_diag_Nfix(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_Nfix_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_Nfix_allocatedInMARBL(autotroph))
                endif

              end do  !  autotroph

            endif  !  config_ecosysTracers_diagnostic_fields_level2

            !
            ! level 3 diagnostics, full depth
            !

            if (config_ecosysTracers_diagnostic_fields_level3) then

              ecosys_diag_SiO2_FLUX_IN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_SiO2_FLUX_IN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_SiO2_FLUX_IN_allocatedInMARBL)

              ecosys_diag_SiO2_PROD(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_SiO2_PROD_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_SiO2_PROD_allocatedInMARBL)

              ecosys_diag_SiO2_REMIN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_SiO2_REMIN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_SiO2_REMIN_allocatedInMARBL)

              ecosys_diag_dust_FLUX_IN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_dust_FLUX_IN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_dust_FLUX_IN_allocatedInMARBL)

              ecosys_diag_dust_REMIN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_dust_REMIN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_dust_REMIN_allocatedInMARBL)

              ecosys_diag_P_iron_FLUX_IN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_P_iron_FLUX_IN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_P_iron_FLUX_IN_allocatedInMARBL)

              ecosys_diag_P_iron_PROD(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_P_iron_PROD_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_P_iron_PROD_allocatedInMARBL)

              ecosys_diag_P_iron_REMIN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_P_iron_REMIN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_P_iron_REMIN_allocatedInMARBL)

              ecosys_diag_DOC_prod(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DOC_prod_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DOC_prod_allocatedInMARBL)

              ecosys_diag_DOC_remin(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DOC_remin_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DOC_remin_allocatedInMARBL)

              ecosys_diag_DOCr_remin(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DOCr_remin_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DOCr_remin_allocatedInMARBL)

              ecosys_diag_DON_prod(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DON_prod_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DON_prod_allocatedInMARBL)

              ecosys_diag_DON_remin(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DON_remin_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DON_remin_allocatedInMARBL)

              ecosys_diag_DONr_remin(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DONr_remin_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DONr_remin_allocatedInMARBL)

              ecosys_diag_DOP_prod(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DOP_prod_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DOP_prod_allocatedInMARBL)

              ecosys_diag_DOP_remin(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DOP_remin_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DOP_remin_allocatedInMARBL)

              ecosys_diag_DOPr_remin(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_DOPr_remin_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_DOPr_remin_allocatedInMARBL)

              ecosys_diag_Fe_scavenge(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Fe_scavenge_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Fe_scavenge_allocatedInMARBL)

              ecosys_diag_Fe_scavenge_rate(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Fe_scavenge_rate_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Fe_scavenge_rate_allocatedInMARBL)

              ecosys_diag_Lig_prod(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Lig_prod_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Lig_prod_allocatedInMARBL)

              ecosys_diag_Lig_loss(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Lig_loss_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Lig_loss_allocatedInMARBL)

              ecosys_diag_Lig_scavenge(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Lig_scavenge_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Lig_scavenge_allocatedInMARBL)

              ecosys_diag_Fefree(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Fefree_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Fefree_allocatedInMARBL)

              ecosys_diag_Lig_photochem(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Lig_photochem_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Lig_photochem_allocatedInMARBL)

              ecosys_diag_Lig_deg(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_Lig_deg_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_Lig_deg_allocatedInMARBL)

              ecosys_diag_tot_bSi_form(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_tot_bSi_form_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_tot_bSi_form_allocatedInMARBL)

              do autotroph = 1, autotroph_cnt
! bSi_form only defined for silicifiers
                if (ecosys_diag_bSi_form_indexInMARBL(autotroph) > 0) then
                ecosys_diag_bSi_form(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_bSi_form_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_bSi_form_allocatedInMARBL(autotroph))
                endif
              enddo  !  autotroph loop

            endif  !  config_ecosysTracers_diagnostic_fields_level3

            !
            ! level 4 diagnostics
            !
            if (config_ecosysTracers_diagnostic_fields_level4) then

              ecosys_diag_CaCO3_PROD(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CaCO3_PROD_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CaCO3_PROD_allocatedInMARBL)

              ecosys_diag_CaCO3_REMIN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CaCO3_REMIN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CaCO3_REMIN_allocatedInMARBL)

              ecosys_diag_CaCO3_ALT_CO2_PROD(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CaCO3_ALT_CO2_PROD_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CaCO3_ALT_CO2_PROD_allocatedInMARBL)

              ecosys_diag_CaCO3_ALT_CO2_REMIN(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CaCO3_ALT_CO2_REMIN_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CaCO3_ALT_CO2_REMIN_allocatedInMARBL)

              ecosys_diag_tot_CaCO3_form(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_tot_CaCO3_form_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_tot_CaCO3_form_allocatedInMARBL)

              ecosys_diag_CO3(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CO3_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CO3_allocatedInMARBL)

              ecosys_diag_HCO3(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_HCO3_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_HCO3_allocatedInMARBL)

              ecosys_diag_H2CO3(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_H2CO3_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_H2CO3_allocatedInMARBL)

              ecosys_diag_CO3_ALT_CO2(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CO3_ALT_CO2_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_CO3_ALT_CO2_allocatedInMARBL)

              ecosys_diag_HCO3_ALT_CO2(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_HCO3_ALT_CO2_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_HCO3_ALT_CO2_allocatedInMARBL)

              ecosys_diag_H2CO3_ALT_CO2(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_H2CO3_ALT_CO2_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_H2CO3_ALT_CO2_allocatedInMARBL)

              ecosys_diag_pH_3D_ALT_CO2(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_pH_3D_ALT_CO2_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_pH_3D_ALT_CO2_allocatedInMARBL)

              ecosys_diag_co3_sat_calc(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_co3_sat_calc_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_co3_sat_calc_allocatedInMARBL)

              ecosys_diag_co3_sat_arag(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_co3_sat_arag_indexInMARBL)%field_3d(iLevel,column),  &
                ecosysMissingValue, ecosys_diag_co3_sat_arag_allocatedInMARBL)

              do autotroph = 1, autotroph_cnt
! CaCO3_form only defined for calcifiers
                if (ecosys_diag_CaCO3_form_indexInMARBL(autotroph) > 0) then
                ecosys_diag_CaCO3_form(autotroph,iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                  diags(ecosys_diag_CaCO3_form_indexInMARBL(autotroph))%field_3d(iLevel,column),  &
                  ecosysMissingValue, ecosys_diag_CaCO3_form_allocatedInMARBL(autotroph))
                endif
              enddo  !  autotroph

            endif  !  config_ecosysTracers_diagnostic_fields_level4

            !
            ! level 5 diagnostics, full depth
            !

            if (config_ecosysTracers_diagnostic_fields_level5) then

! nutrient restoring is currently not implemented

!             ecosys_diag_PO4_RESTORE(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
!               diags(ecosys_diag_PO4_RESTORE_indexInMARBL)%field_3d(iLevel,column),  &
!               ecosysMissingValue, ecosys_diag_PO4_RESTORE_allocatedInMARBL)

!             ecosys_diag_NO3_RESTORE(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
!               diags(ecosys_diag_NO3_RESTORE_indexInMARBL)%field_3d(iLevel,column),  &
!               ecosysMissingValue, ecosys_diag_NO3_RESTORE_allocatedInMARBL)

!             ecosys_diag_SiO3_RESTORE(iLevel,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
!               diags(ecosys_diag_SiO3_RESTORE_indexInMARBL)%field_3d(iLevel,column),  &
!               ecosysMissingValue, ecosys_diag_SiO3_RESTORE_allocatedInMARBL)

            endif  !  config_ecosysTracers_diagnostic_fields_level5

         enddo  !  iLevel

         !
         ! level 1 diagnostics with no depth dependence
         !
         if (config_ecosysTracers_diagnostic_fields_level1) then

           ecosys_diag_photoC_TOT_zint(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_photoC_TOT_zint_indexInMARBL)%field_2d(column)*convertLengthCGStoMKS,  &
             ecosysMissingValue, ecosys_diag_photoC_TOT_zint_allocatedInMARBL)

           ecosys_diag_photoC_NO3_TOT_zint(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_photoC_NO3_TOT_zint_indexInMARBL)%field_2d(column)*convertLengthCGStoMKS,  &
             ecosysMissingValue, ecosys_diag_photoC_NO3_TOT_zint_allocatedInMARBL)

           ecosys_diag_O2_ZMIN(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_O2_ZMIN_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_O2_ZMIN_allocatedInMARBL)

           ecosys_diag_O2_ZMIN_DEPTH(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_O2_ZMIN_DEPTH_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_O2_ZMIN_DEPTH_allocatedInMARBL)

           ecosys_diag_calcToSed(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_calcToSed_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_calcToSed_allocatedInMARBL)

           ecosys_diag_calcToSed_ALT_CO2(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_calcToSed_ALT_CO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_calcToSed_ALT_CO2_allocatedInMARBL)

           ecosys_diag_pocToSed(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_pocToSed_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_pocToSed_allocatedInMARBL)

           ecosys_diag_pfeToSed(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_pfeToSed_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_pfeToSed_allocatedInMARBL)

           ecosys_diag_SedDenitrif(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_SedDenitrif_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_SedDenitrif_allocatedInMARBL)

           ecosys_diag_Jint_Ctot(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_Jint_Ctot_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_Jint_Ctot_allocatedInMARBL)

           ecosys_diag_Jint_Ptot(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_Jint_Ptot_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_Jint_Ptot_allocatedInMARBL)

           ecosys_diag_Jint_Sitot(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_Jint_Sitot_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_Jint_Sitot_allocatedInMARBL)

           ecosys_diag_Jint_Fetot(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_Jint_Fetot_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_Jint_Fetot_allocatedInMARBL)

           do autotroph = 1, autotroph_cnt

             ecosys_diag_photoC_zint(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
               diags(ecosys_diag_photoC_zint_indexInMARBL(autotroph))%field_2d(column)*convertLengthCGStoMKS,  &
               ecosysMissingValue, ecosys_diag_photoC_zint_allocatedInMARBL(autotroph))

             ecosys_diag_photoC_NO3_zint(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
               diags(ecosys_diag_photoC_NO3_zint_indexInMARBL(autotroph))%field_2d(column)*convertLengthCGStoMKS,  &
               ecosysMissingValue, ecosys_diag_photoC_NO3_zint_allocatedInMARBL(autotroph))

           enddo  ! autotroph

         endif  !  config_ecosysTracers_diagnostic_fields_level1

         !
         ! level 2 diagnostics with no depth dependence
         !
         if (config_ecosysTracers_diagnostic_fields_level2) then

           do autotroph = 1, autotroph_cnt

             ecosys_diag_N_lim_surf(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
               diags(ecosys_diag_N_lim_surf_indexInMARBL(autotroph))%field_2d(column),  &
               ecosysMissingValue, ecosys_diag_N_lim_surf_allocatedInMARBL(autotroph))

             ecosys_diag_P_lim_surf(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
               diags(ecosys_diag_P_lim_surf_indexInMARBL(autotroph))%field_2d(column),  &
               ecosysMissingValue, ecosys_diag_P_lim_surf_allocatedInMARBL(autotroph))

             ecosys_diag_Fe_lim_surf(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
               diags(ecosys_diag_Fe_lim_surf_indexInMARBL(autotroph))%field_2d(column),  &
               ecosysMissingValue, ecosys_diag_Fe_lim_surf_allocatedInMARBL(autotroph))

             if(ecosys_diag_SiO3_lim_surf_indexInMARBL(autotroph) > 0)  &
               ecosys_diag_SiO3_lim_surf(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
               diags(ecosys_diag_SiO3_lim_surf_indexInMARBL(autotroph))%field_2d(column),  &
               ecosysMissingValue, ecosys_diag_SiO3_lim_surf_allocatedInMARBL(autotroph))

             ecosys_diag_light_lim_surf(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
               diags(ecosys_diag_light_lim_surf_indexInMARBL(autotroph))%field_2d(column),  &
               ecosysMissingValue, ecosys_diag_light_lim_surf_allocatedInMARBL(autotroph))

           enddo  ! autotroph

         endif  !  config_ecosysTracers_diagnostic_fields_level2

         !
         ! level 3 diagnostics with no depth dependence
         !
         if (config_ecosysTracers_diagnostic_fields_level3) then
           ecosys_diag_ponToSed(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_ponToSed_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_ponToSed_allocatedInMARBL)

           ecosys_diag_popToSed(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_popToSed_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_popToSed_allocatedInMARBL)

           ecosys_diag_bsiToSed(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_bsiToSed_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_bsiToSed_allocatedInMARBL)

           ecosys_diag_dustToSed(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_dustToSed_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_dustToSed_allocatedInMARBL)

           ecosys_diag_OtherRemin(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
             diags(ecosys_diag_OtherRemin_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_OtherRemin_allocatedInMARBL)
         endif  !  config_ecosysTracers_diagnostic_fields_level3

         !
         ! level 4 diagnostics with no depth dependence
         !
         if (config_ecosysTracers_diagnostic_fields_level4) then

            ecosys_diag_tot_CaCO3_form_zint(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
              diags(ecosys_diag_tot_CaCO3_form_zint_indexInMARBL)%field_2d(column)*convertLengthCGStoMKS,  &
              ecosysMissingValue, ecosys_diag_tot_CaCO3_form_zint_allocatedInMARBL)

            ecosys_diag_zsatcalc(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
              diags(ecosys_diag_zsatcalc_indexInMARBL)%field_2d(column),  &
              ecosysMissingValue, ecosys_diag_zsatcalc_allocatedInMARBL)

            ecosys_diag_zsatarag(iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
              diags(ecosys_diag_zsatarag_indexInMARBL)%field_2d(column),  &
              ecosysMissingValue, ecosys_diag_zsatarag_allocatedInMARBL)

            do autotroph = 1, autotroph_cnt
! CaCO3_form_zint only defined for calcifiers
                if (ecosys_diag_CaCO3_form_zint_indexInMARBL(autotroph) > 0) then
              ecosys_diag_CaCO3_form_zint(autotroph,iCell) = merge(marbl_instances(1)%interior_tendency_diags%  &
                diags(ecosys_diag_CaCO3_form_zint_indexInMARBL(autotroph))%field_2d(column)*convertLengthCGStoMKS,  &
                ecosysMissingValue, ecosys_diag_CaCO3_form_zint_allocatedInMARBL(autotroph))
                endif
            enddo  !  autotroph

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
      nTracers, nCellsSolve, zMid, minLevelCell, indexTemperature, indexSalinity, ecosysSurfaceFlux, err)!{{{

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
      integer, dimension(:), intent(in) :: minLevelCell

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

      integer :: numColumns, column, iCell, iLevelSurface, n

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
         dust_FLUX_IN,         &
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
         CO2_alt_gas_flux,         &
         total_Chl

! input/output terms
      real (kind=RKIND), dimension(:), pointer :: &
         PH_PREV,              &
         PH_PREV_ALT_CO2

      !
      ! level 2 diagnostics
      !

      type (mpas_pool_type), pointer ::  &
         ecosysDiagFieldsLevel2

      real (kind=RKIND), dimension(:), pointer :: &
         ecosys_diag_pistonVel_O2,         &
         ecosys_diag_pistonVel_CO2,        &
         ecosys_diag_Schmidt_O2,           &
         ecosys_diag_Schmidt_CO2,          &
         ecosys_diag_O2_saturation,        &
         ecosys_diag_xkw,                  &
         ecosys_diag_u10sqr,               &
         ecosys_diag_CO2star,              &
         ecosys_diag_dCO2star,             &
         ecosys_diag_pCO2surface,          &
         ecosys_diag_dpCO2,                &
         ecosys_diag_CO2star_ALT_CO2,      &
         ecosys_diag_dCO2star_ALT_CO2,     &
         ecosys_diag_pCO2surface_ALT_CO2,  &
         ecosys_diag_dpCO2_ALT_CO2,        &
         ecosys_diag_atmPressure,          &
         ecosys_diag_ifrac

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
      call mpas_pool_get_array(ecosysAuxiliary, 'atmosphericCO2', atmosphericCO2)
      call mpas_pool_get_array(ecosysAuxiliary, 'atmosphericCO2_ALT_CO2', atmosphericCO2_ALT_CO2)

      call mpas_pool_get_array(ecosysAuxiliary, 'CO2_gas_flux', CO2_gas_flux)
      call mpas_pool_get_array(ecosysAuxiliary, 'CO2_alt_gas_flux', CO2_alt_gas_flux)
      call mpas_pool_get_array(ecosysAuxiliary, 'total_Chl', total_Chl)

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

      if (config_ecosysTracers_diagnostic_fields_level2) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysDiagFieldsLevel2', ecosysDiagFieldsLevel2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pistonVel_O2', ecosys_diag_pistonVel_O2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pistonVel_CO2', ecosys_diag_pistonVel_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Schmidt_O2', ecosys_diag_Schmidt_O2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_Schmidt_CO2', ecosys_diag_Schmidt_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_O2_saturation', ecosys_diag_O2_saturation)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_xkw', ecosys_diag_xkw)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_u10sqr', ecosys_diag_u10sqr)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_ifrac', ecosys_diag_ifrac)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_atmPressure', ecosys_diag_atmPressure)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_CO2star', ecosys_diag_CO2star)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dCO2star', ecosys_diag_dCO2star)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pCO2surface', ecosys_diag_pCO2surface)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dpCO2', ecosys_diag_dpCO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_CO2star_ALT_CO2', ecosys_diag_CO2star_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dCO2star_ALT_CO2', ecosys_diag_dCO2star_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_pCO2surface_ALT_CO2', ecosys_diag_pCO2surface_ALT_CO2)
         call mpas_pool_get_array(ecosysDiagFieldsLevel2, 'ecosys_diag_dpCO2_ALT_CO2', ecosys_diag_dpCO2_ALT_CO2)
      endif

!  for surface fluxes, marbl does all surface points in 1 call
      numColumns = 1
      column = 1
      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve

         iLevelSurface = minLevelCell(iCell)

! load pH into saved state
         marbl_instances(1)%surface_flux_saved_state%state(surface_pH_indexInMARBL)%field_2d(iCell) = &
            PH_PREV(iCell)
         marbl_instances(1)%surface_flux_saved_state%state(surface_pH_altCO2_indexInMARBL)%field_2d(iCell) = &
            PH_PREV_ALT_CO2(iCell)

! load surface forcing fields
         marbl_instances(1)%surface_flux_forcings(atmPressure_indexInMARBL)%field_0d(iCell) = &
            atmosphericPressure(iCell)*PascalsToAtmospheres
         marbl_instances(1)%surface_flux_forcings(ifrac_indexInMARBL)%field_0d(iCell) = &
            iceFraction(iCell)
         marbl_instances(1)%surface_flux_forcings(u10sqr_indexInMARBL)%field_0d(iCell) = &
            windSpeedSquared10m(iCell)*convertLengthSquaredMKStoCGS
         marbl_instances(1)%surface_flux_forcings(atmCO2_indexInMARBL)%field_0d(iCell) = &
            atmosphericCO2(iCell)
         marbl_instances(1)%surface_flux_forcings(atmAltCO2_indexInMARBL)%field_0d(iCell) = &
            atmosphericCO2_ALT_CO2(iCell)
         marbl_instances(1)%surface_flux_forcings(sst_indexInMARBL)%field_0d(iCell) = &
            activeTracers(indexTemperature,iLevelSurface,iCell)
         marbl_instances(1)%surface_flux_forcings(sss_indexInMARBL)%field_0d(iCell) = &
            activeTracers(indexSalinity,iLevelSurface,iCell)
! NOTE pass in total Fe and mult by parm_Fe_bioavail inside the flux routine
!  divide river Fe by bioavail since it is already the available to make it total
         marbl_instances(1)%surface_flux_forcings(ironFlux_indexInMARBL)%field_0d(iCell) = &
            IRON_FLUX_IN(iCell)
         marbl_instances(1)%surface_flux_forcings(noxFlux_indexInMARBL)%field_0d(iCell) = &
            depositionFluxNO3(iCell)
         marbl_instances(1)%surface_flux_forcings(nhyFlux_indexInMARBL)%field_0d(iCell) = &
            depositionFluxNH4(iCell)
         marbl_instances(1)%surface_flux_forcings(dustFluxSurface_indexInMARBL)%field_0d(iCell) = &
            dust_FLUX_IN(iCell)

! load surface values of all BGC prognostic variables
         do n = 1,ecosys_tracer_cnt
            marbl_instances(1)%tracers_at_surface(iCell,n) = &
                  ecosysTracers(ecosysIndexMARBLtoMPASO(n),iLevelSurface,iCell)
         end do

      end do   !  iCell

! add in other fluxes after call to surface forcing
! marbl flux routine calculates:
!  gas fluxes of dic, dic_alt, o2
!  alk, alk_alt fluxes related to no3, nh4 deposition
!  can modify nh4 fluxes if the ocean EMITS them (not turned on yet in mpas-ocean)
!  copies no3, nh4, iron deposition fluxes into output flux arrays
!  uses dust flux to calculate po4 and sio3 deposition fluxes

      call marbl_instances(1)%surface_flux_compute()

      if (marbl_instances(1)%StatusLog%labort_marbl) &
        call mpas_log_write("WARNING: problem in MARBL surface_flux_compute", MPAS_LOG_CRIT)
        call marbl_instances(1)%StatusLog%erase()

! copy marbl output into mpas arrays

      column = 1
      do iCell=1,nCellsSolve

! copy updated surface fluxes
         do n = 1, ecosys_tracer_cnt
            ecosysSurfaceFlux(ecosysIndexMARBLtoMPASO(n),iCell) =  &
               marbl_instances(1)%surface_fluxes(iCell,n)
         end do

! gas fluxes (co2, o2) come back in CGS units
         ecosysSurfaceFlux(dic_ind_MPAS,iCell) = ecosysSurfaceFlux(dic_ind_MPAS,iCell)*convertLengthCGStoMKS
         ecosysSurfaceFlux(dic_alt_co2_ind_MPAS,iCell) = ecosysSurfaceFlux(dic_alt_co2_ind_MPAS,iCell)*convertLengthCGStoMKS
         ecosysSurfaceFlux(o2_ind_MPAS,iCell) = ecosysSurfaceFlux(o2_ind_MPAS,iCell)*convertLengthCGStoMKS

! po4 and sio3 deposition fluxes come back in CGS units because they are proportional to dust_flux_in
         ecosysSurfaceFlux(po4_ind_MPAS,iCell) = ecosysSurfaceFlux(po4_ind_MPAS,iCell)*convertLengthCGStoMKS
         ecosysSurfaceFlux(sio3_ind_MPAS,iCell) = ecosysSurfaceFlux(sio3_ind_MPAS,iCell)*convertLengthCGStoMKS

! copy updated pH
         PH_PREV(iCell) = marbl_instances(1)%surface_flux_saved_state%state(surface_pH_indexInMARBL)%field_2d(iCell)
         PH_PREV_ALT_CO2(iCell) = &
            marbl_instances(1)%surface_flux_saved_state%state(surface_pH_altCO2_indexInMARBL)%field_2d(iCell)

! copy updated CO2 flux and total Chl
! gas fluxes come back in CGS units
         CO2_gas_flux(iCell) = marbl_instances(1)%surface_flux_output%sfo(co2_flux_indexInMARBL)%forcing_field(iCell)*convertLengthCGStoMKS

! marbl doesnt return alt_co2_flux anymore !?!?
!        CO2_alt_gas_flux(iCell) = marbl_instances(1)%surface_flux_output%sfo(co2_flux_alt_co2_indexInMARBL)%forcing_field(iCell)
         CO2_alt_gas_flux(iCell) = 0.0_RKIND
         total_Chl(iCell) = marbl_instances(1)%surface_flux_output%sfo(totalChl_indexInMARBL)%forcing_field(iCell)

! add contribution from river fluxes
         ecosysSurfaceFlux(no3_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(no3_ind_MPAS,iCell) + riverFluxNO3(iCell)
         ecosysSurfaceFlux(po4_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(po4_ind_MPAS,iCell) + riverFluxPO4(iCell)
         ecosysSurfaceFlux(don_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(don_ind_MPAS,iCell) + riverFluxDON(iCell) * 0.9_RKIND
         ecosysSurfaceFlux(donr_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(donr_ind_MPAS,iCell) + riverFluxDON(iCell) * 0.1_RKIND
         ecosysSurfaceFlux(dop_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(dop_ind_MPAS,iCell) + riverFluxDOP(iCell) * 0.975_RKIND
         ecosysSurfaceFlux(dopr_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(dopr_ind_MPAS,iCell) + riverFluxDOP(iCell) * 0.025_RKIND
         ecosysSurfaceFlux(sio3_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(sio3_ind_MPAS,iCell) + riverFluxSiO3(iCell)
         ecosysSurfaceFlux(fe_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(fe_ind_MPAS,iCell) + riverFluxFe(iCell) / parm_Fe_bioavail
         ecosysSurfaceFlux(dic_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(dic_ind_MPAS,iCell) + riverFluxDIC(iCell)
         ecosysSurfaceFlux(dic_alt_co2_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(dic_alt_co2_ind_MPAS,iCell) + riverFluxDIC(iCell)
         ecosysSurfaceFlux(alk_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(alk_ind_MPAS,iCell) + riverFluxALK(iCell)
         ecosysSurfaceFlux(alk_alt_co2_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(alk_alt_co2_ind_MPAS,iCell) + riverFluxALK(iCell)
         ecosysSurfaceFlux(doc_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(doc_ind_MPAS,iCell) + riverFluxDOC(iCell) * 0.8_RKIND
         ecosysSurfaceFlux(docr_ind_MPAS,iCell) =  &
            ecosysSurfaceFlux(docr_ind_MPAS,iCell) + riverFluxDOC(iCell) * 0.2_RKIND

! now add contribution from sea ice fluxes if available (assume MKS)
         if (config_use_ecosysTracers_sea_ice_coupling) then
            ecosysSurfaceFlux(no3_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(no3_ind_MPAS,iCell) + iceFluxNO3(iCell)
            ecosysSurfaceFlux(don_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(don_ind_MPAS,iCell) + iceFluxDON(iCell)
            ecosysSurfaceFlux(nh4_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(nh4_ind_MPAS,iCell) + iceFluxNH4(iCell)
            ecosysSurfaceFlux(sio3_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(sio3_ind_MPAS,iCell) + iceFluxSiO3(iCell)
            ecosysSurfaceFlux(fe_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(fe_ind_MPAS,iCell) + iceFluxFeDissolved(iCell)/parm_Fe_bioavail
            ecosysSurfaceFlux(dic_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(dic_ind_MPAS,iCell) + iceFluxDIC(iCell)
            ecosysSurfaceFlux(dic_alt_co2_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(dic_alt_co2_ind_MPAS,iCell) + iceFluxDIC(iCell)
            ecosysSurfaceFlux(doc_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(doc_ind_MPAS,iCell) + iceFluxDOC(1,iCell) + iceFluxDOC(2,iCell) &
               + iceFluxDOC(3,iCell)
! this is just a placeholder
            ecosysSurfaceFlux(donr_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(donr_ind_MPAS,iCell) + iceFluxDOCr(iCell)

            ecosysSurfaceFlux(diatC_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(diatC_ind_MPAS,iCell) + iceFluxPhytoC(1,iCell)
            ecosysSurfaceFlux(spC_ind_MPAS,iCell) =  &
               ecosysSurfaceFlux(spC_ind_MPAS,iCell) + iceFluxPhytoC(2,iCell)
         endif

!explicitly set the rest of the fluxes to 0
         ecosysSurfaceFlux(lig_ind_MPAS,iCell)      = 0.0_RKIND
         ecosysSurfaceFlux(zooC_ind_MPAS,iCell)     = 0.0_RKIND
         ecosysSurfaceFlux(spChl_ind_MPAS,iCell)    = 0.0_RKIND
         ecosysSurfaceFlux(spFe_ind_MPAS,iCell)     = 0.0_RKIND
         ecosysSurfaceFlux(spCaCO3_ind_MPAS,iCell)  = 0.0_RKIND
         ecosysSurfaceFlux(spP_ind_MPAS,iCell)      = 0.0_RKIND
         ecosysSurfaceFlux(diatChl_ind_MPAS,iCell)  = 0.0_RKIND
         ecosysSurfaceFlux(diatFe_ind_MPAS,iCell)   = 0.0_RKIND
         ecosysSurfaceFlux(diatSi_ind_MPAS,iCell)   = 0.0_RKIND
         ecosysSurfaceFlux(diatP_ind_MPAS,iCell)    = 0.0_RKIND
         ecosysSurfaceFlux(diazC_ind_MPAS,iCell)    = 0.0_RKIND
         ecosysSurfaceFlux(diazChl_ind_MPAS,iCell)  = 0.0_RKIND
         ecosysSurfaceFlux(diazFe_ind_MPAS,iCell)   = 0.0_RKIND
         ecosysSurfaceFlux(diazP_ind_MPAS,iCell)    = 0.0_RKIND

! level 2 diagnostics of surface fields

         if (config_ecosysTracers_diagnostic_fields_level2) then

           ecosys_diag_pistonVel_O2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_pistonVel_O2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_pistonVel_O2_allocatedInMARBL)

           ecosys_diag_pistonVel_CO2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_pistonVel_CO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_pistonVel_CO2_allocatedInMARBL)

           ecosys_diag_Schmidt_O2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_Schmidt_O2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_Schmidt_O2_allocatedInMARBL)

           ecosys_diag_Schmidt_CO2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_Schmidt_CO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_Schmidt_CO2_allocatedInMARBL)

           ecosys_diag_O2_saturation(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_O2_saturation_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_O2_saturation_allocatedInMARBL)

           ecosys_diag_xkw(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_xkw_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_xkw_allocatedInMARBL)

! no u10sqr diagnostic in marbl so use forcing field
           ecosys_diag_u10sqr(iCell) = marbl_instances(1)%surface_flux_forcings(u10sqr_indexInMARBL)%field_0d(column)*convertLengthSquaredCGStoMKS

           ecosys_diag_ifrac(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_ifrac_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_ifrac_allocatedInMARBL)

! convert back to Pa
           ecosys_diag_atmPressure(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_atmPressure_indexInMARBL)%field_2d(column)/PascalsToAtmospheres,  &
             ecosysMissingValue, ecosys_diag_atmPressure_allocatedInMARBL)

           ecosys_diag_CO2star(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_CO2star_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_CO2star_allocatedInMARBL)

           ecosys_diag_dCO2star(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_dCO2star_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_dCO2star_allocatedInMARBL)

           ecosys_diag_pCO2surface(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_pCO2surface_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_pCO2surface_allocatedInMARBL)

           ecosys_diag_dpCO2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_dpCO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_dpCO2_allocatedInMARBL)

           ecosys_diag_CO2star_ALT_CO2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_CO2star_ALT_CO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_CO2star_ALT_CO2_allocatedInMARBL)

           ecosys_diag_dCO2star_ALT_CO2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_dCO2star_ALT_CO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_dCO2star_ALT_CO2_allocatedInMARBL)

           ecosys_diag_pCO2surface_ALT_CO2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_pCO2surface_ALT_CO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_pCO2surface_ALT_CO2_allocatedInMARBL)

           ecosys_diag_dpCO2_ALT_CO2(iCell) = merge(marbl_instances(1)%surface_flux_diags%  &
             diags(ecosys_diag_dpCO2_ALT_CO2_indexInMARBL)%field_2d(column),  &
             ecosysMissingValue, ecosys_diag_dpCO2_ALT_CO2_allocatedInMARBL)

         endif  !  level 2 diagnostics

         column = column + 1
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

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: ecosysAuxiliary
      type (mpas_pool_type), pointer :: ecosysMonthlyForcing
      type (mpas_pool_type), pointer :: verticalMeshPool

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        ecosysTracers

      real (kind=RKIND), dimension(:), pointer :: &
         refZMid, refLayerThickness

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
      integer :: nTracers, numColumnsMax, iLevel, n

      ! scalar pointers
      integer, pointer :: nVertLevels, index_dummy, nCellsSolve

      character(len=strKIND) :: &
         forcingIntervalMonthly,  &
         forcingReferenceTimeMonthly,  &
         fieldIndexNameForMPAS,  &
         fieldNameForMPAS

! maltrud marbl char_len = 384
      character(len=384) :: &
         fieldNameInMARBL

      logical :: foundIndex

      real (kind=RKIND) :: &
         zTop, zBot

      real (kind=RKIND), allocatable, dimension(:) :: &
         zt, zw, dz

      real (kind=RKIND) :: &
         convertLengthScale = 100.0_RKIND

      !
      ! get tracers pools
      !

      err = 0

      if (config_use_ecosysTracers) then

      !
      ! Get tracer group so we can get the number of tracers in it
      !

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, 1)

    !--------------------------------------------------------------------
    !  MARBL setup is 3 steps:
    !  1) Configure (set variables that affect tracer count / other parms)
    !  2) Initialize ("lock" config vars so the aren't changed during init
    !     or in the time loop; write config vars to log; set parameters)
    !  3) Complete setup ("lock" parmameters so they aren't changed during
    !     time loop; write parameters to log)
    !
    !  can set up saved state, tracers, and forcing fields after (3)
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    !  get necessary vertical mesh quantities
    !  for init, just use the reference vertical grid
    !--------------------------------------------------------------------

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
      call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

      allocate(zt(nVertLevels), zw(nVertLevels), dz(nVertLevels))
      zTop = 0.0_RKIND
      do iLevel=1, nVertLevels
         zt(iLevel) = -1.0_RKIND*refZMid(iLevel)*convertLengthScale
         dz(iLevel) = refLayerThickness(iLevel)*convertLengthScale
         zBot = zTop + refLayerThickness(iLevel)
         zw(iLevel) = zBot*convertLengthScale
         zTop = zBot
      enddo


    !--------------------------------------------------------------------
    !  Configure marbl
    !--------------------------------------------------------------------

! NOTE below didnt work so hardwired in marbl_settings_mod.F90
!    call marbl_instances%put_setting('ciso_on', .false.)
!    call marbl_instances%put_setting('ciso_lsource_sink', .false.)

    !--------------------------------------------------------------------
    !  Initialize marbl
    !--------------------------------------------------------------------

     call marbl_instances(1)%init(                      &
          gcm_num_levels = nVertLevels,                 &
          gcm_num_PAR_subcols = 1,                      &
          gcm_num_elements_surface_flux = nCellsSolve,  &
          gcm_delta_z = dz,                             &
          gcm_zw = zw,                                  &
          gcm_zt = zt,                                  &
          lgcm_has_global_ops = .false.)

     if (marbl_instances(1)%StatusLog%labort_marbl) &
        call mpas_log_write("WARNING: problem in MARBL initialization", MPAS_LOG_CRIT)
        call marbl_instances(1)%StatusLog%erase()

     ! Make sure MARBL tracer count lines up with what MPASO expects
      nTracers = size(ecosysTracers, dim=1)

! NOTE we dont actually have a value for marbl_actual_tracer_cnt
!     if (marbl_actual_tracer_cnt /= nTracers) then
!      write(log_message,"(A,I0,A,I0)") 'MARBL is computing tendencies for ', &
!                 marbl_actual_tracer_cnt, ' tracers, but MPAS-OCEAN is expecting ', &
!                 ecosys_tracer_cnt
!        err = 1
!        return
!     endif

      deallocate(zt, zw, dz)

    !--------------------------------------------------------------------
    !  make connections between marbl and mpaso
    !--------------------------------------------------------------------

! figure out correspondence of marbl to mpaso specific indices
     allocate(ecosysIndexMARBLtoMPASO(ecosys_tracer_cnt))
     allocate(ecosysIndexMPASOtoMARBL(ecosys_tracer_cnt))
     ecosysIndexMARBLtoMPASO = -999
     ecosysIndexMPASOtoMARBL = -999

     do n = 1, ecosys_tracer_cnt

        fieldNameInMARBL = marbl_instances(1)%tracer_metadata(n)%short_name
        fieldNameForMPAS = trim(fieldNameInMARBL)
        fieldIndexNameForMPAS = 'index_'//trim(fieldNameInMARBL)
        logMessage = ' looking for '//trim(fieldIndexNameForMPAS)
        if (verbose) call mpas_log_write(trim(logMessage))
        call mpas_pool_get_dimension(tracersPool, trim(fieldIndexNameForMPAS), index_dummy)
        if (index_dummy > 0) then
           ecosysIndexMARBLtoMPASO(n) = index_dummy
           ecosysIndexMPASOtoMARBL(index_dummy) = n
           logMessage = ' found '//trim(fieldIndexNameForMPAS)
           if (verbose) call mpas_log_write(trim(logMessage))
!          write(*,*)' MARBL index = ', n
!          write(*,*)' MPASO index = ', ecosysIndexMARBLtoMPASO(n)
           call ocn_ecosys_assign_MPAS_indices_prognostic(fieldNameForMPAS, index_dummy)
        endif
    end do

! make sure all field indices were found
    do n = 1, ecosys_tracer_cnt
       if (ecosysIndexMARBLtoMPASO(n) == -999) call mpas_log_write(  &
          "' could not find MARBL index for '//trim(marbl_instances(1)%tracer_metadata(n)%short_name)",MPAS_LOG_CRIT)
    enddo
    if (any(ecosysIndexMARBLtoMPASO <= 0)) then
       err = 1
       return
    endif

! set indices for surface forcing
    u10sqr_indexInMARBL          = marbl_instances(1)%surface_flux_forcing_ind%u10_sqr_id
    ifrac_indexInMARBL           = marbl_instances(1)%surface_flux_forcing_ind%ifrac_id
    sst_indexInMARBL             = marbl_instances(1)%surface_flux_forcing_ind%sst_id
    sss_indexInMARBL             = marbl_instances(1)%surface_flux_forcing_ind%sss_id
    atmPressure_indexInMARBL     = marbl_instances(1)%surface_flux_forcing_ind%atm_pressure_id
    atmCO2_indexInMARBL          = marbl_instances(1)%surface_flux_forcing_ind%xco2_id
    atmAltCO2_indexInMARBL       = marbl_instances(1)%surface_flux_forcing_ind%xco2_alt_co2_id
    dustFluxSurface_indexInMARBL = marbl_instances(1)%surface_flux_forcing_ind%dust_flux_id
    ironFlux_indexInMARBL        = marbl_instances(1)%surface_flux_forcing_ind%iron_flux_id
    noxFlux_indexInMARBL         = marbl_instances(1)%surface_flux_forcing_ind%nox_flux_id
    nhyFlux_indexInMARBL         = marbl_instances(1)%surface_flux_forcing_ind%nhy_flux_id

! set indices for interior forcing
    dustFlux_indexInMARBL    = marbl_instances(1)%interior_tendency_forcing_ind%dustflux_id
    shortwave_indexInMARBL   = marbl_instances(1)%interior_tendency_forcing_ind%surf_shortwave_id
    PAR_col_frac_indexInMARBL = marbl_instances(1)%interior_tendency_forcing_ind%PAR_col_frac_id
    temperature_indexInMARBL = marbl_instances(1)%interior_tendency_forcing_ind%potemp_id
    salinity_indexInMARBL    = marbl_instances(1)%interior_tendency_forcing_ind%salinity_id
    pressure_indexInMARBL    = marbl_instances(1)%interior_tendency_forcing_ind%pressure_id
    fesedFlux_indexInMARBL   = marbl_instances(1)%interior_tendency_forcing_ind%fesedflux_id

    call ocn_ecosys_init_indices

! search for interior diagnostic field indices
    ecosys_interior_diagnostics_cnt = size(marbl_instances(1)%interior_tendency_diags%diags)
    do n = 1, ecosys_interior_diagnostics_cnt
       fieldNameInMARBL = marbl_instances(1)%interior_tendency_diags%diags(n)%short_name
       fieldNameForMPAS = trim(fieldNameInMARBL)
       index_dummy = n
       foundIndex = .false.
       call ocn_ecosys_assign_MPAS_indices_diagnostic(fieldNameForMPAS, index_dummy, foundIndex)
       logMessage = ' could not find MARBL interior diagnostic index for '//trim(fieldNameForMPAS)
       if (.not. foundIndex) call mpas_log_write(trim(logMessage))
    end do

! search for surface diagnostic field indices
    ecosys_surface_forcing_diagnostics_cnt = size(marbl_instances(1)%surface_flux_diags%diags)
    do n = 1, ecosys_surface_forcing_diagnostics_cnt
       fieldNameInMARBL = marbl_instances(1)%surface_flux_diags%diags(n)%short_name
       fieldNameForMPAS = trim(fieldNameInMARBL)
       index_dummy = n
       foundIndex = .false.
       call ocn_ecosys_assign_MPAS_indices_diagnostic(fieldNameForMPAS, index_dummy, foundIndex)
       logMessage = ' could not find MARBL surface diagnostic index for '//trim(fieldNameForMPAS)
       if (.not. foundIndex) call mpas_log_write(trim(logMessage))
    end do

! no marbl diagnostic for u10sqr so use the forcing index instead
    ecosys_diag_u10sqr_indexInMARBL = u10sqr_indexInMARBL

    if (verbose) call ocn_ecosys_echo_indices

      !
      ! for now only do 1 column at a time
      !
    numColumnsMax = 1

      ! identify saved state indices as defined by MARBL
    surface_pH_indexInMARBL = marbl_instances(1)%surf_state_ind%ph_surf
    surface_pH_altCO2_indexInMARBL = marbl_instances(1)%surf_state_ind%ph_alt_co2_surf
    interior_pH_indexInMARBL = marbl_instances(1)%interior_state_ind%ph_col
    interior_pH_altCO2_indexInMARBL = marbl_instances(1)%interior_state_ind%ph_alt_co2_col

      ! tell MARBL we want CO2 flux and total Chl field for (short-wave absorption) as additional outputs

    call marbl_instances(1)%surface_flux_output%add_sfo(   &
          num_elements = nCellsSolve,                      &
          field_name   = "flux_co2",                       &
          sfo_id       = co2_flux_indexInMARBL,            &
          marbl_status_log = marbl_instances(1)%StatusLog)

    if (marbl_instances(1)%StatusLog%labort_marbl) &
        call mpas_log_write("WARNING: problem in MARBL flux_co2 initialization", MPAS_LOG_CRIT)
        call marbl_instances(1)%StatusLog%erase()

! marbl doesnt return flux_alt_co2 anymore ???!!!
!      call marbl_instances(1)%surface_flux_output%add_sfo(  &
!             num_elements = nCellsSolve,                    &
!             field_name   = "flux_alt_co2",                 &
!             sfo_id       = co2_flux_alt_co2_indexInMARBL,  &
!             marbl_status_log = marbl_instances(1)%StatusLog)

    call marbl_instances(1)%surface_flux_output%add_sfo(  &
          num_elements = nCellsSolve,                    &
          field_name   = "totalChl",                     &
          sfo_id       = totalChl_indexInMARBL,          &
          marbl_status_log = marbl_instances(1)%StatusLog)

      if (marbl_instances(1)%StatusLog%labort_marbl) &
        call mpas_log_write("WARNING: problem in MARBL totalChl initialization", MPAS_LOG_CRIT)
        call marbl_instances(1)%StatusLog%erase()

! initialize monthly forcing to be read from file

      forcingIntervalMonthly = "0000-01-00_00:00:00"
      forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

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
      if (.not. config_use_ecosysTracers_river_inputs_from_coupler) then
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxNO3', riverFluxNO3)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxPO4', riverFluxPO4)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDON', riverFluxDON)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOP', riverFluxDOP)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxFe', riverFluxFe)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDIC', riverFluxDIC)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxALK', riverFluxALK)
        call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOC', riverFluxDOC)
      endif

      call mpas_pool_get_subpool(domain % blocklist % structs, 'ecosysMonthlyForcing', ecosysMonthlyForcing)

      call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNO3', depositionFluzNO3)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNH4', depositionFluzNH4)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'IRON_FLUZ_IN', IRON_FLUZ_IN)
      call mpas_pool_get_array(ecosysMonthlyForcing, 'dust_FLUZ_IN', dust_FLUZ_IN)
      if (.not. config_use_ecosysTracers_river_inputs_from_coupler) then
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzNO3', riverFluzNO3)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzPO4', riverFluzPO4)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDON', riverFluzDON)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOP', riverFluzDOP)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzSiO3', riverFluzSiO3)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzFe', riverFluzFe)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDIC', riverFluzDIC)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzALK', riverFluzALK)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOC', riverFluzDOC)
      endif

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

! only read in river fluxes if not getting them from coupler
      if (.not. config_use_ecosysTracers_river_inputs_from_coupler) then

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

      endif  ! config_use_ecosysTracers_river_inputs_from_coupler

      call MPAS_forcing_init_field_data( forcingGroupHead, &
        'ecosysMonthlyClimatology',  &
         domain % streamManager,  &
         config_do_restart, &
         .false.)

      end if  !  config_use_ecosysTracers)

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

        call mpas_set_timeInterval(timeStepEco,timeString=config_dt)
        call mpas_get_timeInterval(timeStepEco,dt=dt)

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
        if (.not. config_use_ecosysTracers_river_inputs_from_coupler) then
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxNO3', riverFluxNO3)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxPO4', riverFluxPO4)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDON', riverFluxDON)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOP', riverFluxDOP)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxFe', riverFluxFe)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDIC', riverFluxDIC)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxALK', riverFluxALK)
          call mpas_pool_get_array(ecosysAuxiliary, 'riverFluxDOC', riverFluxDOC)
        endif

        call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNO3', depositionFluzNO3)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'depositionFluzNH4', depositionFluzNH4)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'IRON_FLUZ_IN', IRON_FLUZ_IN)
        call mpas_pool_get_array(ecosysMonthlyForcing, 'dust_FLUZ_IN', dust_FLUZ_IN)
        if (.not. config_use_ecosysTracers_river_inputs_from_coupler) then
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzNO3', riverFluzNO3)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzPO4', riverFluzPO4)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDON', riverFluzDON)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOP', riverFluzDOP)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzSiO3', riverFluzSiO3)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzFe', riverFluzFe)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDIC', riverFluzDIC)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzALK', riverFluzALK)
          call mpas_pool_get_array(ecosysMonthlyForcing, 'riverFluzDOC', riverFluzDOC)
        endif

        do iCell = 1, nCells
           depositionFluxNO3(iCell) = depositionFluzNO3(iCell)
           depositionFluxNH4(iCell) = depositionFluzNH4(iCell)
           IRON_FLUX_IN(iCell)      = IRON_FLUZ_IN(iCell)
           dust_FLUX_IN(iCell)      = dust_FLUZ_IN(iCell)
        enddo

        if (.not. config_use_ecosysTracers_river_inputs_from_coupler) then
          do iCell = 1, nCells
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
        endif

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

!***********************************************************************
!
!  routine ocn_ecosys_assign_MPAS_indices_prognostic
!
!> \brief   find the specific field index corresponding to what is assigned in the ecosysTracers structure
!> \author  Mathew Maltrud
!> \date    07/20/2017

!
!-----------------------------------------------------------------------

   subroutine ocn_ecosys_assign_MPAS_indices_prognostic(field_name, index_dummy)!{{{

      character(len=strKIND), intent(in) :: field_name
      integer, intent(inout) :: index_dummy

      if (trim(field_name) == 'PO4')      po4_ind_MPAS            = index_dummy
      if (trim(field_name) == 'NO3')      no3_ind_MPAS            = index_dummy
      if (trim(field_name) == 'SiO3')     sio3_ind_MPAS           = index_dummy
      if (trim(field_name) == 'NH4')      nh4_ind_MPAS            = index_dummy
      if (trim(field_name) == 'Fe')       fe_ind_MPAS             = index_dummy
      if (trim(field_name) == 'Lig')      lig_ind_MPAS            = index_dummy
      if (trim(field_name) == 'O2')       o2_ind_MPAS             = index_dummy
      if (trim(field_name) == 'DIC')      dic_ind_MPAS            = index_dummy
      if (trim(field_name) == 'DIC_ALT_CO2') dic_alt_co2_ind_MPAS = index_dummy
      if (trim(field_name) == 'ALK')      alk_ind_MPAS            = index_dummy
      if (trim(field_name) == 'ALK_ALT_CO2') alk_alt_co2_ind_MPAS = index_dummy
      if (trim(field_name) == 'DOC')      doc_ind_MPAS            = index_dummy
      if (trim(field_name) == 'DON')      don_ind_MPAS            = index_dummy
      if (trim(field_name) == 'DOP')      dop_ind_MPAS            = index_dummy
      if (trim(field_name) == 'DOPr')     dopr_ind_MPAS           = index_dummy
      if (trim(field_name) == 'DONr')     donr_ind_MPAS           = index_dummy
      if (trim(field_name) == 'DOCr')     docr_ind_MPAS           = index_dummy
      if (trim(field_name) == 'zooC')     zooC_ind_MPAS           = index_dummy
      if (trim(field_name) == 'spChl')    spChl_ind_MPAS          = index_dummy
      if (trim(field_name) == 'spC')      spC_ind_MPAS            = index_dummy
      if (trim(field_name) == 'spFe')     spFe_ind_MPAS           = index_dummy
      if (trim(field_name) == 'spP')      spP_ind_MPAS            = index_dummy
      if (trim(field_name) == 'spCaCO3')  spCaCO3_ind_MPAS        = index_dummy
      if (trim(field_name) == 'diatChl')  diatChl_ind_MPAS        = index_dummy
      if (trim(field_name) == 'diatC')    diatC_ind_MPAS          = index_dummy
      if (trim(field_name) == 'diatFe')   diatFe_ind_MPAS         = index_dummy
      if (trim(field_name) == 'diatP')    diatP_ind_MPAS          = index_dummy
      if (trim(field_name) == 'diatSi')   diatSi_ind_MPAS         = index_dummy
      if (trim(field_name) == 'diazChl')  diazChl_ind_MPAS        = index_dummy
      if (trim(field_name) == 'diazC')    diazC_ind_MPAS          = index_dummy
      if (trim(field_name) == 'diazFe')   diazFe_ind_MPAS         = index_dummy
      if (trim(field_name) == 'diazP')    diazP_ind_MPAS          = index_dummy

    end subroutine ocn_ecosys_assign_MPAS_indices_prognostic!}}}

!***********************************************************************
!
!  routine ocn_ecosys_assign_MPAS_indices_diagnostic
!
!> \brief   find the specific field index corresponding to what is assigned in the ecosysTracers structure
!> \author  Mathew Maltrud
!> \date    07/20/2017

!
!-----------------------------------------------------------------------

   subroutine ocn_ecosys_assign_MPAS_indices_diagnostic(field_name, indexInMARBL, found)!{{{

      character(strKIND), intent(in) :: field_name
      integer, intent(inout) :: indexInMARBL
      logical, intent(inout) :: found

! level 1 diagnostics, no depth dependence

      if (trim(field_name) == 'photoC_TOT_zint') then
         ecosys_diag_photoC_TOT_zint_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_TOT_zint_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_TOT_zint') then
         ecosys_diag_photoC_NO3_TOT_zint_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_NO3_TOT_zint_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'O2_ZMIN') then
         ecosys_diag_O2_ZMIN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_O2_ZMIN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'O2_ZMIN_DEPTH') then
         ecosys_diag_O2_ZMIN_DEPTH_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_O2_ZMIN_DEPTH_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'calcToSed') then
         ecosys_diag_calcToSed_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_calcToSed_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'calcToSed_ALT_CO2') then
         ecosys_diag_calcToSed_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_calcToSed_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'pocToSed') then
         ecosys_diag_pocToSed_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_pocToSed_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'pfeToSed') then
         ecosys_diag_pfeToSed_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_pfeToSed_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'SedDenitrif') then
         ecosys_diag_SedDenitrif_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_SedDenitrif_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Jint_Ctot') then
         ecosys_diag_Jint_Ctot_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Jint_Ctot_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Jint_Ntot') then
         ecosys_diag_Jint_Ntot_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Jint_Ntot_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Jint_Ptot') then
         ecosys_diag_Jint_Ptot_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Jint_Ptot_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Jint_Sitot') then
         ecosys_diag_Jint_Sitot_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Jint_Sitot_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Jint_Fetot') then
         ecosys_diag_Jint_Fetot_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Jint_Fetot_allocatedInMARBL = .true.
         found = .true.
      endif

! small
      if (trim(field_name) == 'photoC_sp_zint') then
         ecosys_diag_photoC_zint_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_zint_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_sp_zint') then
         ecosys_diag_photoC_NO3_zint_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_NO3_zint_allocatedInMARBL(1) = .true.
         found = .true.
      endif
! diat
      if (trim(field_name) == 'photoC_diat_zint') then
         ecosys_diag_photoC_zint_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_zint_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_diat_zint') then
         ecosys_diag_photoC_NO3_zint_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_NO3_zint_allocatedInMARBL(2) = .true.
         found = .true.
      endif
! diaz
      if (trim(field_name) == 'photoC_diaz_zint') then
         ecosys_diag_photoC_zint_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_zint_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_diaz_zint') then
         ecosys_diag_photoC_NO3_zint_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_photoC_NO3_zint_allocatedInMARBL(3) = .true.
         found = .true.
      endif

! level 1 diagnostics, full depth

      if (trim(field_name) == 'PAR_avg') then
         ecosys_diag_PAR_avg_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_PAR_avg_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POC_FLUX_IN') then
         ecosys_diag_POC_FLUX_IN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POC_FLUX_IN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CaCO3_FLUX_IN') then
         ecosys_diag_CaCO3_FLUX_IN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_FLUX_IN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CaCO3_ALT_CO2_FLUX_IN') then
         ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_auto_TOT') then
         ecosys_diag_graze_auto_TOT_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_graze_auto_TOT_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'zoo_loss') then
         ecosys_diag_zoo_loss_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_zoo_loss_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'zoo_loss_poc') then
         ecosys_diag_zoo_loss_poc_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_zoo_loss_poc_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'zoo_loss_doc') then
         ecosys_diag_zoo_loss_doc_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_zoo_loss_doc_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_zoo') then
         ecosys_diag_graze_zoo_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_graze_zoo_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_zoo_poc') then
         ecosys_diag_graze_zoo_poc_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_graze_zoo_poc_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_zoo_doc') then
         ecosys_diag_graze_zoo_doc_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_graze_zoo_doc_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_zoo_zoo') then
         ecosys_diag_graze_zoo_zoo_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_graze_zoo_zoo_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'x_graze_zoo') then
         ecosys_diag_x_graze_zoo_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_x_graze_zoo_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_TOT') then
         ecosys_diag_photoC_TOT_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_TOT_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_TOT') then
         ecosys_diag_photoC_NO3_TOT_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_NO3_TOT_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'NITRIF') then
         ecosys_diag_NITRIF_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_NITRIF_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DENITRIF') then
         ecosys_diag_DENITRIF_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DENITRIF_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Nfix') then
         ecosys_diag_tot_Nfix_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_tot_Nfix_allocatedInMARBL = .true.
         found = .true.
      endif

! level 2 diagnostics, no depth dependence

      if (trim(field_name) == 'PV_O2') then
         ecosys_diag_pistonVel_O2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_pistonVel_O2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'PV_CO2') then
         ecosys_diag_pistonVel_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_pistonVel_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'SCHMIDT_O2') then
         ecosys_diag_Schmidt_O2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Schmidt_O2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'SCHMIDT_CO2') then
         ecosys_diag_Schmidt_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Schmidt_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'O2SAT') then
         ecosys_diag_O2_saturation_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_O2_saturation_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'ECOSYS_XKW') then
         ecosys_diag_xkw_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_xkw_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CO2STAR') then
         ecosys_diag_CO2star_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_CO2star_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DCO2STAR') then
         ecosys_diag_dCO2star_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_dCO2star_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'pCO2SURF') then
         ecosys_diag_pCO2surface_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_pCO2surface_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DpCO2') then
         ecosys_diag_dpCO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_dpCO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CO2STAR_ALT_CO2') then
         ecosys_diag_CO2star_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_CO2star_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DCO2STAR_ALT_CO2') then
         ecosys_diag_dCO2star_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_dCO2star_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'pCO2SURF_ALT_CO2') then
         ecosys_diag_pCO2surface_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_pCO2surface_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DpCO2_ALT_CO2') then
         ecosys_diag_dpCO2_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_dpCO2_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif

      if (trim(field_name) == 'ECOSYS_IFRAC') then
         ecosys_diag_ifrac_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_ifrac_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'ECOSYS_ATM_PRESS') then
         ecosys_diag_atmPressure_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%surface_flux_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_atmPressure_allocatedInMARBL = .true.
         found = .true.
      endif

! level 2 diagnostics, full depth

      if (trim(field_name) == 'O2_PRODUCTION') then
         ecosys_diag_O2_PRODUCTION_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_O2_PRODUCTION_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'O2_CONSUMPTION') then
         ecosys_diag_O2_CONSUMPTION_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_O2_CONSUMPTION_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'AOU') then
         ecosys_diag_AOU_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_AOU_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'pH_3D') then
         ecosys_diag_pH_3D_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_pH_3D_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POC_PROD') then
         ecosys_diag_POC_PROD_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POC_PROD_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POC_REMIN_DOCr') then
         ecosys_diag_POC_REMIN_DOCr_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POC_REMIN_DOCr_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POC_REMIN_DIC') then
         ecosys_diag_POC_REMIN_DIC_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POC_REMIN_DIC_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POP_FLUX_IN') then
         ecosys_diag_POP_FLUX_IN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POP_FLUX_IN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POP_PROD') then
         ecosys_diag_POP_PROD_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POP_PROD_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POP_REMIN_DOPr') then
         ecosys_diag_POP_REMIN_DOPr_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POP_REMIN_DOPr_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'POP_REMIN_PO4') then
         ecosys_diag_POP_REMIN_PO4_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_POP_REMIN_PO4_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'PON_REMIN_DONr') then
         ecosys_diag_PON_REMIN_DONr_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_PON_REMIN_DONr_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'PON_REMIN_NH4') then
         ecosys_diag_PON_REMIN_NH4_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_PON_REMIN_NH4_allocatedInMARBL = .true.
         found = .true.
      endif

! small
      if (trim(field_name) == 'sp_Qp') then
         ecosys_diag_Qp_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Qp_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_N_lim_surf') then
         ecosys_diag_N_lim_surf_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_N_lim_surf_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_P_lim_surf') then
         ecosys_diag_P_lim_surf_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_P_lim_surf_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_Fe_lim_surf') then
         ecosys_diag_Fe_lim_surf_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Fe_lim_surf_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_SiO3_lim_surf') then
         ecosys_diag_SiO3_lim_surf_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_SiO3_lim_surf_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_light_lim_surf') then
         ecosys_diag_light_lim_surf_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_light_lim_surf_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_sp') then
         ecosys_diag_photoC_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_sp') then
         ecosys_diag_photoC_NO3_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_NO3_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoFe_sp') then
         ecosys_diag_photoFe_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoFe_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoNO3_sp') then
         ecosys_diag_photoNO3_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoNO3_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoNH4_sp') then
         ecosys_diag_photoNH4_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoNH4_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOP_sp_uptake') then
         ecosys_diag_DOP_uptake_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOP_uptake_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'PO4_sp_uptake') then
         ecosys_diag_PO4_uptake_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_PO4_uptake_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_sp') then
         ecosys_diag_auto_graze_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_sp_poc') then
         ecosys_diag_auto_graze_poc_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_poc_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_sp_doc') then
         ecosys_diag_auto_graze_doc_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_doc_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_sp_zoo') then
         ecosys_diag_auto_graze_zoo_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_zoo_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_loss') then
         ecosys_diag_auto_loss_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_loss_poc') then
         ecosys_diag_auto_loss_poc_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_poc_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_loss_doc') then
         ecosys_diag_auto_loss_doc_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_doc_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_agg') then
         ecosys_diag_auto_agg_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_agg_allocatedInMARBL(1) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'sp_Nfix') then
         ecosys_diag_Nfix_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Nfix_allocatedInMARBL(1) = .true.
         found = .true.
      endif
! diat
      if (trim(field_name) == 'diat_Qp') then
         ecosys_diag_Qp_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Qp_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_N_lim_surf') then
         ecosys_diag_N_lim_surf_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_N_lim_surf_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_P_lim_surf') then
         ecosys_diag_P_lim_surf_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_P_lim_surf_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_Fe_lim_surf') then
         ecosys_diag_Fe_lim_surf_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Fe_lim_surf_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_SiO3_lim_surf') then
         ecosys_diag_SiO3_lim_surf_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_SiO3_lim_surf_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_light_lim_surf') then
         ecosys_diag_light_lim_surf_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_light_lim_surf_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_diat') then
         ecosys_diag_photoC_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_diat') then
         ecosys_diag_photoC_NO3_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_NO3_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoFe_diat') then
         ecosys_diag_photoFe_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoFe_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoNO3_diat') then
         ecosys_diag_photoNO3_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoNO3_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoNH4_diat') then
         ecosys_diag_photoNH4_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoNH4_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOP_diat_uptake') then
         ecosys_diag_DOP_uptake_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOP_uptake_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'PO4_diat_uptake') then
         ecosys_diag_PO4_uptake_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_PO4_uptake_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diat') then
         ecosys_diag_auto_graze_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diat_poc') then
         ecosys_diag_auto_graze_poc_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_poc_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diat_doc') then
         ecosys_diag_auto_graze_doc_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_doc_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diat_zoo') then
         ecosys_diag_auto_graze_zoo_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_zoo_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_loss') then
         ecosys_diag_auto_loss_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_loss_poc') then
         ecosys_diag_auto_loss_poc_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_poc_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_loss_doc') then
         ecosys_diag_auto_loss_doc_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_doc_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_agg') then
         ecosys_diag_auto_agg_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_agg_allocatedInMARBL(2) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diat_Nfix') then
         ecosys_diag_Nfix_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Nfix_allocatedInMARBL(2) = .true.
         found = .true.
      endif
! diaz
      if (trim(field_name) == 'diaz_Qp') then
         ecosys_diag_Qp_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Qp_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_N_lim_surf') then
         ecosys_diag_N_lim_surf_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_N_lim_surf_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_P_lim_surf') then
         ecosys_diag_P_lim_surf_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_P_lim_surf_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_Fe_lim_surf') then
         ecosys_diag_Fe_lim_surf_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_Fe_lim_surf_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_SiO3_lim_surf') then
         ecosys_diag_SiO3_lim_surf_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_SiO3_lim_surf_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_light_lim_surf') then
         ecosys_diag_light_lim_surf_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_light_lim_surf_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_diaz') then
         ecosys_diag_photoC_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoC_NO3_diaz') then
         ecosys_diag_photoC_NO3_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoC_NO3_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoFe_diaz') then
         ecosys_diag_photoFe_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoFe_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoNO3_diaz') then
         ecosys_diag_photoNO3_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoNO3_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'photoNH4_diaz') then
         ecosys_diag_photoNH4_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_photoNH4_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOP_diaz_uptake') then
         ecosys_diag_DOP_uptake_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOP_uptake_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'PO4_diaz_uptake') then
         ecosys_diag_PO4_uptake_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_PO4_uptake_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diaz') then
         ecosys_diag_auto_graze_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diaz_poc') then
         ecosys_diag_auto_graze_poc_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_poc_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diaz_doc') then
         ecosys_diag_auto_graze_doc_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_doc_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'graze_diaz_zoo') then
         ecosys_diag_auto_graze_zoo_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_graze_zoo_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_loss') then
         ecosys_diag_auto_loss_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_loss_poc') then
         ecosys_diag_auto_loss_poc_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_poc_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_loss_doc') then
         ecosys_diag_auto_loss_doc_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_loss_doc_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_agg') then
         ecosys_diag_auto_agg_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_auto_agg_allocatedInMARBL(3) = .true.
         found = .true.
      endif
      if (trim(field_name) == 'diaz_Nfix') then
         ecosys_diag_Nfix_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Nfix_allocatedInMARBL(3) = .true.
         found = .true.
      endif

! level 3 diagnostics, no depth dependence

      if (trim(field_name) == 'ponToSed') then
         ecosys_diag_ponToSed_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_ponToSed_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'popToSed') then
         ecosys_diag_popToSed_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_popToSed_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'bsiToSed') then
         ecosys_diag_bsiToSed_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_bsiToSed_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'dustToSed') then
         ecosys_diag_dustToSed_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_dustToSed_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'OtherRemin') then
         ecosys_diag_OtherRemin_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_OtherRemin_allocatedInMARBL = .true.
         found = .true.
      endif

! level 3 diagnostics, full depth

      if (trim(field_name) == 'SiO2_FLUX_IN') then
         ecosys_diag_SiO2_FLUX_IN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_SiO2_FLUX_IN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'SiO2_PROD') then
         ecosys_diag_SiO2_PROD_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_SiO2_PROD_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'SiO2_REMIN') then
         ecosys_diag_SiO2_REMIN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_SiO2_REMIN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'dust_FLUX_IN') then
         ecosys_diag_dust_FLUX_IN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_dust_FLUX_IN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'dust_REMIN') then
         ecosys_diag_dust_REMIN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_dust_REMIN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'P_iron_FLUX_IN') then
         ecosys_diag_P_iron_FLUX_IN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_P_iron_FLUX_IN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'P_iron_PROD') then
         ecosys_diag_P_iron_PROD_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_P_iron_PROD_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'P_iron_REMIN') then
         ecosys_diag_P_iron_REMIN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_P_iron_REMIN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOC_prod') then
         ecosys_diag_DOC_prod_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOC_prod_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOC_remin') then
         ecosys_diag_DOC_remin_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOC_remin_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOCr_remin') then
         ecosys_diag_DOCr_remin_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOCr_remin_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DON_prod') then
         ecosys_diag_DON_prod_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DON_prod_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DON_remin') then
         ecosys_diag_DON_remin_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DON_remin_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOP_prod') then
         ecosys_diag_DOP_prod_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOP_prod_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOP_remin') then
         ecosys_diag_DOP_remin_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOP_remin_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Fe_scavenge') then
         ecosys_diag_Fe_scavenge_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Fe_scavenge_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Fe_scavenge_rate') then
         ecosys_diag_Fe_scavenge_rate_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Fe_scavenge_rate_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Lig_prod') then
         ecosys_diag_Lig_prod_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Lig_prod_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Lig_loss') then
         ecosys_diag_Lig_loss_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Lig_loss_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Lig_scavenge') then
         ecosys_diag_Lig_scavenge_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Lig_scavenge_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Fefree') then
         ecosys_diag_Fefree_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Fefree_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Lig_photochem') then
         ecosys_diag_Lig_photochem_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Lig_photochem_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'Lig_deg') then
         ecosys_diag_Lig_deg_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_Lig_deg_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DONr_remin') then
         ecosys_diag_DONr_remin_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DONr_remin_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'DOPr_remin') then
         ecosys_diag_DOPr_remin_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_DOPr_remin_allocatedInMARBL = .true.
         found = .true.
      endif

      if (trim(field_name) == 'bSi_form') then
         ecosys_diag_tot_bSi_form_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_tot_bSi_form_allocatedInMARBL = .true.
         found = .true.
      endif

! small
      if (trim(field_name) == 'sp_bSi_form') then
         ecosys_diag_bSi_form_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_bSi_form_allocatedInMARBL(1) = .true.
         found = .true.
      endif
! diat
      if (trim(field_name) == 'diat_bSi_form') then
         ecosys_diag_bSi_form_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_bSi_form_allocatedInMARBL(2) = .true.
         found = .true.
      endif
! diaz
      if (trim(field_name) == 'diaz_bSi_form') then
         ecosys_diag_bSi_form_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_bSi_form_allocatedInMARBL(3) = .true.
         found = .true.
      endif

! level 4 diagnostics, no depth dependence

      if (trim(field_name) == 'CaCO3_form_zint') then
         ecosys_diag_tot_CaCO3_form_zint_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_tot_CaCO3_form_zint_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'zsatcalc') then
         ecosys_diag_zsatcalc_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_zsatcalc_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'zsatarag') then
         ecosys_diag_zsatarag_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_zsatarag_allocatedInMARBL = .true.
         found = .true.
      endif

! level 4 diagnostics, no depth dependence, each autotroph

! small
      if (trim(field_name) == 'sp_CaCO3_form_zint') then
         ecosys_diag_CaCO3_form_zint_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_CaCO3_form_zint_allocatedInMARBL(1) = .true.
         found = .true.
      endif
! diat
      if (trim(field_name) == 'diat_CaCO3_form_zint') then
         ecosys_diag_CaCO3_form_zint_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_CaCO3_form_zint_allocatedInMARBL(2) = .true.
         found = .true.
      endif
! diaz
      if (trim(field_name) == 'diaz_CaCO3_form_zint') then
         ecosys_diag_CaCO3_form_zint_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_2d) ) &
            ecosys_diag_CaCO3_form_zint_allocatedInMARBL(3) = .true.
         found = .true.
      endif

! level 4 diagnostics, full depth

      if (trim(field_name) == 'CaCO3_PROD') then
         ecosys_diag_CaCO3_PROD_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_PROD_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CaCO3_REMIN') then
         ecosys_diag_CaCO3_REMIN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_REMIN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CaCO3_ALT_CO2_PROD') then
         ecosys_diag_CaCO3_ALT_CO2_PROD_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_ALT_CO2_PROD_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CaCO3_ALT_CO2_REMIN') then
         ecosys_diag_CaCO3_ALT_CO2_REMIN_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_ALT_CO2_REMIN_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CaCO3_form') then
         ecosys_diag_tot_CaCO3_form_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_tot_CaCO3_form_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CO3') then
         ecosys_diag_CO3_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CO3_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'HCO3') then
         ecosys_diag_HCO3_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_HCO3_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'H2CO3') then
         ecosys_diag_H2CO3_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_H2CO3_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'CO3_ALT_CO2') then
         ecosys_diag_CO3_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CO3_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'HCO3_ALT_CO2') then
         ecosys_diag_HCO3_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_HCO3_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'H2CO3_ALT_CO2') then
         ecosys_diag_H2CO3_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_HCO3_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'pH_3D_ALT_CO2') then
         ecosys_diag_pH_3D_ALT_CO2_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_pH_3D_ALT_CO2_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'co3_sat_calc') then
         ecosys_diag_co3_sat_calc_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_co3_sat_calc_allocatedInMARBL = .true.
         found = .true.
      endif
      if (trim(field_name) == 'co3_sat_arag') then
         ecosys_diag_co3_sat_arag_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_co3_sat_calc_allocatedInMARBL = .true.
         found = .true.
      endif

! level 4 diagnostics, full depth, each autotroph

! small
      if (trim(field_name) == 'sp_CaCO3_form') then
         ecosys_diag_CaCO3_form_indexInMARBL(1) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_form_allocatedInMARBL(1) = .true.
         found = .true.
      endif
! diat
      if (trim(field_name) == 'diat_CaCO3_form') then
         ecosys_diag_CaCO3_form_indexInMARBL(2) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_form_allocatedInMARBL(2) = .true.
         found = .true.
      endif
! diaz
      if (trim(field_name) == 'diaz_CaCO3_form') then
         ecosys_diag_CaCO3_form_indexInMARBL(3) = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_CaCO3_form_allocatedInMARBL(3) = .true.
         found = .true.
      endif

! level 5 diagnostics, full depth

      if (trim(field_name) == 'PO4_RESTORE') then
         ecosys_diag_PO4_RESTORE_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_PO4_RESTORE_allocatedInMARBL = .true.
         found = .true.
      endif

      if (trim(field_name) == 'NO3_RESTORE') then
         ecosys_diag_NO3_RESTORE_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_NO3_RESTORE_allocatedInMARBL = .true.
         found = .true.
      endif

      if (trim(field_name) == 'SiO3_RESTORE') then
         ecosys_diag_SiO3_RESTORE_indexInMARBL = indexInMARBL
         if ( allocated( marbl_instances(1)%interior_tendency_diags%diags(indexInMARBL)%field_3d) ) &
            ecosys_diag_SiO3_RESTORE_allocatedInMARBL = .true.
         found = .true.
      endif

    end subroutine ocn_ecosys_assign_MPAS_indices_diagnostic!}}}

!***********************************************************************
!
!  routine ocn_ecosys_init_indices
!
!-----------------------------------------------------------------------

   subroutine ocn_ecosys_init_indices!{{{

! initialize indices so we can check to see if any have not been defined

   ecosys_diag_photoC_TOT_zint_indexInMARBL = -999
   ecosys_diag_photoC_NO3_TOT_zint_indexInMARBL = -999
   ecosys_diag_O2_ZMIN_indexInMARBL = -999
   ecosys_diag_O2_ZMIN_DEPTH_indexInMARBL = -999
   ecosys_diag_calcToSed_indexInMARBL = -999
   ecosys_diag_calcToSed_ALT_CO2_indexInMARBL = -999
   ecosys_diag_pocToSed_indexInMARBL = -999
   ecosys_diag_pfeToSed_indexInMARBL = -999
   ecosys_diag_SedDenitrif_indexInMARBL = -999
   ecosys_diag_Jint_Ctot_indexInMARBL = -999
   ecosys_diag_Jint_Ntot_indexInMARBL = -999
   ecosys_diag_Jint_Ptot_indexInMARBL = -999
   ecosys_diag_Jint_Sitot_indexInMARBL = -999
   ecosys_diag_Jint_Fetot_indexInMARBL = -999

   ecosys_diag_photoC_zint_indexInMARBL(1) = -999
   ecosys_diag_photoC_NO3_zint_indexInMARBL(1) = -999
   ecosys_diag_photoC_zint_indexInMARBL(2) = -999
   ecosys_diag_photoC_NO3_zint_indexInMARBL(2) = -999
   ecosys_diag_photoC_zint_indexInMARBL(3) = -999
   ecosys_diag_photoC_NO3_zint_indexInMARBL(3) = -999

   ecosys_diag_PAR_avg_indexInMARBL = -999
   ecosys_diag_POC_FLUX_IN_indexInMARBL = -999
   ecosys_diag_CaCO3_FLUX_IN_indexInMARBL = -999
   ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_indexInMARBL = -999
   ecosys_diag_graze_auto_TOT_indexInMARBL = -999
   ecosys_diag_zoo_loss_indexInMARBL = -999
   ecosys_diag_zoo_loss_poc_indexInMARBL = -999
   ecosys_diag_zoo_loss_doc_indexInMARBL = -999
   ecosys_diag_graze_zoo_indexInMARBL = -999
   ecosys_diag_graze_zoo_poc_indexInMARBL = -999
   ecosys_diag_graze_zoo_doc_indexInMARBL = -999
   ecosys_diag_graze_zoo_zoo_indexInMARBL = -999
   ecosys_diag_x_graze_zoo_indexInMARBL = -999
   ecosys_diag_photoC_TOT_indexInMARBL = -999
   ecosys_diag_photoC_NO3_TOT_indexInMARBL = -999
   ecosys_diag_NITRIF_indexInMARBL = -999
   ecosys_diag_DENITRIF_indexInMARBL = -999
   ecosys_diag_tot_Nfix_indexInMARBL = -999

   ecosys_diag_pistonVel_O2_indexInMARBL = -999
   ecosys_diag_pistonVel_CO2_indexInMARBL = -999
   ecosys_diag_Schmidt_O2_indexInMARBL = -999
   ecosys_diag_Schmidt_CO2_indexInMARBL = -999
   ecosys_diag_O2_saturation_indexInMARBL = -999
   ecosys_diag_xkw_indexInMARBL = -999
   ecosys_diag_CO2star_indexInMARBL = -999
   ecosys_diag_dCO2star_indexInMARBL = -999
   ecosys_diag_pCO2surface_indexInMARBL = -999
   ecosys_diag_dpCO2_indexInMARBL = -999
   ecosys_diag_CO2star_ALT_CO2_indexInMARBL = -999
   ecosys_diag_dCO2star_ALT_CO2_indexInMARBL = -999
   ecosys_diag_pCO2surface_ALT_CO2_indexInMARBL = -999
   ecosys_diag_dpCO2_ALT_CO2_indexInMARBL = -999

   ecosys_diag_ifrac_indexInMARBL = -999
   ecosys_diag_atmPressure_indexInMARBL = -999

   ecosys_diag_O2_PRODUCTION_indexInMARBL = -999
   ecosys_diag_O2_CONSUMPTION_indexInMARBL = -999
   ecosys_diag_AOU_indexInMARBL = -999
   ecosys_diag_pH_3D_indexInMARBL = -999
   ecosys_diag_POC_PROD_indexInMARBL = -999
   ecosys_diag_POC_REMIN_DOCr_indexInMARBL = -999
   ecosys_diag_POC_REMIN_DIC_indexInMARBL = -999
   ecosys_diag_POP_FLUX_IN_indexInMARBL = -999
   ecosys_diag_POP_PROD_indexInMARBL = -999
   ecosys_diag_POP_REMIN_DOPr_indexInMARBL = -999
   ecosys_diag_POP_REMIN_PO4_indexInMARBL = -999
   ecosys_diag_PON_REMIN_DONr_indexInMARBL = -999
   ecosys_diag_PON_REMIN_NH4_indexInMARBL = -999

   ecosys_diag_Qp_indexInMARBL(1) = -999
   ecosys_diag_N_lim_surf_indexInMARBL(1) = -999
   ecosys_diag_P_lim_surf_indexInMARBL(1) = -999
   ecosys_diag_Fe_lim_surf_indexInMARBL(1) = -999
   ecosys_diag_SiO3_lim_surf_indexInMARBL(1) = -999
   ecosys_diag_light_lim_surf_indexInMARBL(1) = -999
   ecosys_diag_photoC_indexInMARBL(1) = -999
   ecosys_diag_photoC_NO3_indexInMARBL(1) = -999
   ecosys_diag_photoFe_indexInMARBL(1) = -999
   ecosys_diag_photoNO3_indexInMARBL(1) = -999
   ecosys_diag_photoNH4_indexInMARBL(1) = -999
   ecosys_diag_DOP_uptake_indexInMARBL(1) = -999
   ecosys_diag_PO4_uptake_indexInMARBL(1) = -999
   ecosys_diag_auto_graze_indexInMARBL(1) = -999
   ecosys_diag_auto_graze_poc_indexInMARBL(1) = -999
   ecosys_diag_auto_graze_doc_indexInMARBL(1) = -999
   ecosys_diag_auto_graze_zoo_indexInMARBL(1) = -999
   ecosys_diag_auto_loss_indexInMARBL(1) = -999
   ecosys_diag_auto_loss_poc_indexInMARBL(1) = -999
   ecosys_diag_auto_loss_doc_indexInMARBL(1) = -999
   ecosys_diag_auto_agg_indexInMARBL(1) = -999
   ecosys_diag_Nfix_indexInMARBL(2) = -999
   ecosys_diag_Qp_indexInMARBL(2) = -999
   ecosys_diag_N_lim_surf_indexInMARBL(2) = -999
   ecosys_diag_P_lim_surf_indexInMARBL(2) = -999
   ecosys_diag_Fe_lim_surf_indexInMARBL(2) = -999
   ecosys_diag_SiO3_lim_surf_indexInMARBL(2) = -999
   ecosys_diag_light_lim_surf_indexInMARBL(2) = -999
   ecosys_diag_photoC_indexInMARBL(2) = -999
   ecosys_diag_photoC_NO3_indexInMARBL(2) = -999
   ecosys_diag_photoFe_indexInMARBL(2) = -999
   ecosys_diag_photoNO3_indexInMARBL(2) = -999
   ecosys_diag_photoNH4_indexInMARBL(2) = -999
   ecosys_diag_DOP_uptake_indexInMARBL(2) = -999
   ecosys_diag_PO4_uptake_indexInMARBL(2) = -999
   ecosys_diag_auto_graze_indexInMARBL(2) = -999
   ecosys_diag_auto_graze_poc_indexInMARBL(2) = -999
   ecosys_diag_auto_graze_doc_indexInMARBL(2) = -999
   ecosys_diag_auto_graze_zoo_indexInMARBL(2) = -999
   ecosys_diag_auto_loss_indexInMARBL(2) = -999
   ecosys_diag_auto_loss_poc_indexInMARBL(2) = -999
   ecosys_diag_auto_loss_doc_indexInMARBL(2) = -999
   ecosys_diag_auto_agg_indexInMARBL(2) = -999
   ecosys_diag_Nfix_indexInMARBL(2) = -999
   ecosys_diag_Qp_indexInMARBL(3) = -999
   ecosys_diag_N_lim_surf_indexInMARBL(3) = -999
   ecosys_diag_P_lim_surf_indexInMARBL(3) = -999
   ecosys_diag_Fe_lim_surf_indexInMARBL(3) = -999
   ecosys_diag_SiO3_lim_surf_indexInMARBL(3) = -999
   ecosys_diag_light_lim_surf_indexInMARBL(3) = -999
   ecosys_diag_photoC_indexInMARBL(3) = -999
   ecosys_diag_photoC_NO3_indexInMARBL(3) = -999
   ecosys_diag_photoFe_indexInMARBL(3) = -999
   ecosys_diag_photoNO3_indexInMARBL(3) = -999
   ecosys_diag_photoNH4_indexInMARBL(3) = -999
   ecosys_diag_DOP_uptake_indexInMARBL(3) = -999
   ecosys_diag_PO4_uptake_indexInMARBL(3) = -999
   ecosys_diag_auto_graze_indexInMARBL(3) = -999
   ecosys_diag_auto_graze_poc_indexInMARBL(3) = -999
   ecosys_diag_auto_graze_doc_indexInMARBL(3) = -999
   ecosys_diag_auto_graze_zoo_indexInMARBL(3) = -999
   ecosys_diag_auto_loss_indexInMARBL(3) = -999
   ecosys_diag_auto_loss_poc_indexInMARBL(3) = -999
   ecosys_diag_auto_loss_doc_indexInMARBL(3) = -999
   ecosys_diag_auto_agg_indexInMARBL(3) = -999
   ecosys_diag_Nfix_indexInMARBL(3) = -999

   ecosys_diag_ponToSed_indexInMARBL = -999
   ecosys_diag_popToSed_indexInMARBL = -999
   ecosys_diag_bsiToSed_indexInMARBL = -999
   ecosys_diag_dustToSed_indexInMARBL = -999
   ecosys_diag_OtherRemin_indexInMARBL = -999
   ecosys_diag_tot_bSi_form_indexInMARBL = -999

   ecosys_diag_SiO2_FLUX_IN_indexInMARBL = -999
   ecosys_diag_SiO2_PROD_indexInMARBL = -999
   ecosys_diag_SiO2_REMIN_indexInMARBL = -999
   ecosys_diag_dust_FLUX_IN_indexInMARBL = -999
   ecosys_diag_dust_REMIN_indexInMARBL = -999
   ecosys_diag_P_iron_FLUX_IN_indexInMARBL = -999
   ecosys_diag_P_iron_PROD_indexInMARBL = -999
   ecosys_diag_P_iron_REMIN_indexInMARBL = -999
   ecosys_diag_DOC_prod_indexInMARBL = -999
   ecosys_diag_DOC_remin_indexInMARBL = -999
   ecosys_diag_DOCr_remin_indexInMARBL = -999
   ecosys_diag_DON_prod_indexInMARBL = -999
   ecosys_diag_DON_remin_indexInMARBL = -999
   ecosys_diag_DOP_prod_indexInMARBL = -999
   ecosys_diag_DOP_remin_indexInMARBL = -999
   ecosys_diag_Fe_scavenge_indexInMARBL = -999
   ecosys_diag_Fe_scavenge_rate_indexInMARBL = -999
   ecosys_diag_Lig_prod_indexInMARBL = -999
   ecosys_diag_Lig_loss_indexInMARBL = -999
   ecosys_diag_Lig_scavenge_indexInMARBL = -999
   ecosys_diag_Fefree_indexInMARBL = -999
   ecosys_diag_Lig_photochem_indexInMARBL = -999
   ecosys_diag_Lig_deg_indexInMARBL = -999
   ecosys_diag_DONr_remin_indexInMARBL = -999
   ecosys_diag_DOPr_remin_indexInMARBL = -999

   ecosys_diag_bSi_form_indexInMARBL(1) = -999
   ecosys_diag_bSi_form_indexInMARBL(2) = -999
   ecosys_diag_bSi_form_indexInMARBL(3) = -999

   ecosys_diag_tot_CaCO3_form_zint_indexInMARBL = -999
   ecosys_diag_zsatcalc_indexInMARBL = -999
   ecosys_diag_zsatarag_indexInMARBL = -999

   ecosys_diag_CaCO3_form_zint_indexInMARBL(1) = -999
   ecosys_diag_CaCO3_form_zint_indexInMARBL(2) = -999
   ecosys_diag_CaCO3_form_zint_indexInMARBL(3) = -999

   ecosys_diag_CaCO3_PROD_indexInMARBL = -999
   ecosys_diag_CaCO3_REMIN_indexInMARBL = -999
   ecosys_diag_CaCO3_ALT_CO2_PROD_indexInMARBL = -999
   ecosys_diag_CaCO3_ALT_CO2_REMIN_indexInMARBL = -999
   ecosys_diag_tot_CaCO3_form_indexInMARBL = -999
   ecosys_diag_CO3_indexInMARBL = -999
   ecosys_diag_HCO3_indexInMARBL = -999
   ecosys_diag_H2CO3_indexInMARBL = -999
   ecosys_diag_CO3_ALT_CO2_indexInMARBL = -999
   ecosys_diag_HCO3_ALT_CO2_indexInMARBL = -999
   ecosys_diag_H2CO3_ALT_CO2_indexInMARBL = -999
   ecosys_diag_pH_3D_ALT_CO2_indexInMARBL = -999
   ecosys_diag_co3_sat_calc_indexInMARBL = -999
   ecosys_diag_co3_sat_arag_indexInMARBL = -999

   ecosys_diag_CaCO3_form_indexInMARBL(1) = -999
   ecosys_diag_CaCO3_form_indexInMARBL(2) = -999
   ecosys_diag_CaCO3_form_indexInMARBL(3) = -999

    end subroutine ocn_ecosys_init_indices!}}}

!***********************************************************************
!
!  routine ocn_ecosys_echo_indices
!
!-----------------------------------------------------------------------

   subroutine ocn_ecosys_echo_indices!{{{

! check for any indices that have not been defined
! NOTE this writes to ALL processors

 write(*,*) 'ecosys_diag_photoC_TOT_zint_indexInMARBL = ', 	ecosys_diag_photoC_TOT_zint_indexInMARBL
 write(*,*) 'ecosys_diag_photoC_NO3_TOT_zint_indexInMARBL = ', 	ecosys_diag_photoC_NO3_TOT_zint_indexInMARBL
 write(*,*) 'ecosys_diag_O2_ZMIN_indexInMARBL = ', 	ecosys_diag_O2_ZMIN_indexInMARBL
 write(*,*) 'ecosys_diag_O2_ZMIN_DEPTH_indexInMARBL = ', 	ecosys_diag_O2_ZMIN_DEPTH_indexInMARBL
 write(*,*) 'ecosys_diag_Jint_Ctot_indexInMARBL = ', 	ecosys_diag_Jint_Ctot_indexInMARBL
 write(*,*) 'ecosys_diag_Jint_Ntot_indexInMARBL = ', 	ecosys_diag_Jint_Ntot_indexInMARBL
 write(*,*) 'ecosys_diag_Jint_Ptot_indexInMARBL = ', 	ecosys_diag_Jint_Ptot_indexInMARBL
 write(*,*) 'ecosys_diag_Jint_Sitot_indexInMARBL = ', 	ecosys_diag_Jint_Sitot_indexInMARBL
 write(*,*) 'ecosys_diag_Jint_Fetot_indexInMARBL = ', 	ecosys_diag_Jint_Fetot_indexInMARBL
	
 write(*,*) 'ecosys_diag_photoC_zint_indexInMARBL(1) = ', 	ecosys_diag_photoC_zint_indexInMARBL(1)
 write(*,*) 'ecosys_diag_photoC_NO3_zint_indexInMARBL(1) = ', 	ecosys_diag_photoC_NO3_zint_indexInMARBL(1)
 write(*,*) 'ecosys_diag_photoC_zint_indexInMARBL(2) = ', 	ecosys_diag_photoC_zint_indexInMARBL(2)
 write(*,*) 'ecosys_diag_photoC_NO3_zint_indexInMARBL(2) = ', 	ecosys_diag_photoC_NO3_zint_indexInMARBL(2)
 write(*,*) 'ecosys_diag_photoC_zint_indexInMARBL(3) = ', 	ecosys_diag_photoC_zint_indexInMARBL(3)
 write(*,*) 'ecosys_diag_photoC_NO3_zint_indexInMARBL(3) = ', 	ecosys_diag_photoC_NO3_zint_indexInMARBL(3)
	
 write(*,*) 'ecosys_diag_PAR_avg_indexInMARBL = ', 	ecosys_diag_PAR_avg_indexInMARBL
 write(*,*) 'ecosys_diag_POC_FLUX_IN_indexInMARBL = ', 	ecosys_diag_POC_FLUX_IN_indexInMARBL
 write(*,*) 'ecosys_diag_CaCO3_FLUX_IN_indexInMARBL = ', 	ecosys_diag_CaCO3_FLUX_IN_indexInMARBL
 write(*,*) 'ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_indexInMARBL = ', 	ecosys_diag_CaCO3_ALT_CO2_FLUX_IN_indexInMARBL
 write(*,*) 'ecosys_diag_graze_auto_TOT_indexInMARBL = ', 	ecosys_diag_graze_auto_TOT_indexInMARBL
 write(*,*) 'ecosys_diag_zoo_loss_indexInMARBL = ', 	ecosys_diag_zoo_loss_indexInMARBL
 write(*,*) 'ecosys_diag_zoo_loss_poc_indexInMARBL = ', 	ecosys_diag_zoo_loss_poc_indexInMARBL
 write(*,*) 'ecosys_diag_zoo_loss_doc_indexInMARBL = ', 	ecosys_diag_zoo_loss_doc_indexInMARBL
 write(*,*) 'ecosys_diag_graze_zoo_indexInMARBL = ', 	ecosys_diag_graze_zoo_indexInMARBL
 write(*,*) 'ecosys_diag_graze_zoo_poc_indexInMARBL = ', 	ecosys_diag_graze_zoo_poc_indexInMARBL
 write(*,*) 'ecosys_diag_graze_zoo_doc_indexInMARBL = ', 	ecosys_diag_graze_zoo_doc_indexInMARBL
 write(*,*) 'ecosys_diag_graze_zoo_zoo_indexInMARBL = ', 	ecosys_diag_graze_zoo_zoo_indexInMARBL
 write(*,*) 'ecosys_diag_x_graze_zoo_indexInMARBL = ', 	ecosys_diag_x_graze_zoo_indexInMARBL
 write(*,*) 'ecosys_diag_photoC_TOT_indexInMARBL = ', 	ecosys_diag_photoC_TOT_indexInMARBL
 write(*,*) 'ecosys_diag_photoC_NO3_TOT_indexInMARBL = ', 	ecosys_diag_photoC_NO3_TOT_indexInMARBL
 write(*,*) 'ecosys_diag_NITRIF_indexInMARBL = ', 	ecosys_diag_NITRIF_indexInMARBL
 write(*,*) 'ecosys_diag_DENITRIF_indexInMARBL = ', 	ecosys_diag_DENITRIF_indexInMARBL
 write(*,*) 'ecosys_diag_calcToSed_indexInMARBL = ', 	ecosys_diag_calcToSed_indexInMARBL
 write(*,*) 'ecosys_diag_calcToSed_ALT_CO2_indexInMARBL = ', 	ecosys_diag_calcToSed_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_pocToSed_indexInMARBL = ', 	ecosys_diag_pocToSed_indexInMARBL
 write(*,*) 'ecosys_diag_pfeToSed_indexInMARBL = ', 	ecosys_diag_pfeToSed_indexInMARBL
 write(*,*) 'ecosys_diag_SedDenitrif_indexInMARBL = ', 	ecosys_diag_SedDenitrif_indexInMARBL
 write(*,*) 'ecosys_diag_tot_Nfix_indexInMARBL = ', 	ecosys_diag_tot_Nfix_indexInMARBL
	
 write(*,*) 'ecosys_diag_pistonVel_O2_indexInMARBL = ', 	ecosys_diag_pistonVel_O2_indexInMARBL
 write(*,*) 'ecosys_diag_pistonVel_CO2_indexInMARBL = ', 	ecosys_diag_pistonVel_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_Schmidt_O2_indexInMARBL = ', 	ecosys_diag_Schmidt_O2_indexInMARBL
 write(*,*) 'ecosys_diag_Schmidt_CO2_indexInMARBL = ', 	ecosys_diag_Schmidt_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_O2_saturation_indexInMARBL = ', 	ecosys_diag_O2_saturation_indexInMARBL
 write(*,*) 'ecosys_diag_xkw_indexInMARBL = ', 	ecosys_diag_xkw_indexInMARBL
 write(*,*) 'ecosys_diag_CO2star_indexInMARBL = ', 	ecosys_diag_CO2star_indexInMARBL
 write(*,*) 'ecosys_diag_dCO2star_indexInMARBL = ', 	ecosys_diag_dCO2star_indexInMARBL
 write(*,*) 'ecosys_diag_pCO2surface_indexInMARBL = ', 	ecosys_diag_pCO2surface_indexInMARBL
 write(*,*) 'ecosys_diag_dpCO2_indexInMARBL = ', 	ecosys_diag_dpCO2_indexInMARBL
 write(*,*) 'ecosys_diag_CO2star_ALT_CO2_indexInMARBL = ', 	ecosys_diag_CO2star_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_dCO2star_ALT_CO2_indexInMARBL = ', 	ecosys_diag_dCO2star_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_pCO2surface_ALT_CO2_indexInMARBL = ', 	ecosys_diag_pCO2surface_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_dpCO2_ALT_CO2_indexInMARBL = ', 	ecosys_diag_dpCO2_ALT_CO2_indexInMARBL

 write(*,*) 'ecosys_diag_ifrac_indexInMARBL = ',  ecosys_diag_ifrac_indexInMARBL
 write(*,*) 'ecosys_diag_atmPressure_indexInMARBL = ', ecosys_diag_atmPressure_indexInMARBL
	
 write(*,*) 'ecosys_diag_O2_PRODUCTION_indexInMARBL = ', 	ecosys_diag_O2_PRODUCTION_indexInMARBL
 write(*,*) 'ecosys_diag_O2_CONSUMPTION_indexInMARBL = ', 	ecosys_diag_O2_CONSUMPTION_indexInMARBL
 write(*,*) 'ecosys_diag_AOU_indexInMARBL = ', 	ecosys_diag_AOU_indexInMARBL
 write(*,*) 'ecosys_diag_pH_3D_indexInMARBL = ', 	ecosys_diag_pH_3D_indexInMARBL
 write(*,*) 'ecosys_diag_POC_PROD_indexInMARBL = ', 	ecosys_diag_POC_PROD_indexInMARBL
 write(*,*) 'ecosys_diag_POC_REMIN_DOCr_indexInMARBL = ', 	ecosys_diag_POC_REMIN_DOCr_indexInMARBL
 write(*,*) 'ecosys_diag_POC_REMIN_DIC_indexInMARBL = ', 	ecosys_diag_POC_REMIN_DIC_indexInMARBL
 write(*,*) 'ecosys_diag_POP_FLUX_IN_indexInMARBL = ', 	ecosys_diag_POP_FLUX_IN_indexInMARBL
 write(*,*) 'ecosys_diag_POP_PROD_indexInMARBL = ', 	ecosys_diag_POP_PROD_indexInMARBL
 write(*,*) 'ecosys_diag_POP_REMIN_DOPr_indexInMARBL = ', 	ecosys_diag_POP_REMIN_DOPr_indexInMARBL
 write(*,*) 'ecosys_diag_POP_REMIN_PO4_indexInMARBL = ', 	ecosys_diag_POP_REMIN_PO4_indexInMARBL
 write(*,*) 'ecosys_diag_PON_REMIN_DONr_indexInMARBL = ', 	ecosys_diag_PON_REMIN_DONr_indexInMARBL
 write(*,*) 'ecosys_diag_PON_REMIN_NH4_indexInMARBL = ', 	ecosys_diag_PON_REMIN_NH4_indexInMARBL
	
 write(*,*) 'ecosys_diag_Qp_indexInMARBL(1) = ', 	ecosys_diag_Qp_indexInMARBL(1)
 write(*,*) 'ecosys_diag_N_lim_surf_indexInMARBL(1) = ', 	ecosys_diag_N_lim_surf_indexInMARBL(1)
 write(*,*) 'ecosys_diag_P_lim_surf_indexInMARBL(1) = ', 	ecosys_diag_P_lim_surf_indexInMARBL(1)
 write(*,*) 'ecosys_diag_Fe_lim_surf_indexInMARBL(1) = ', 	ecosys_diag_Fe_lim_surf_indexInMARBL(1)
 write(*,*) 'ecosys_diag_SiO3_lim_surf_indexInMARBL(1) = ', 	ecosys_diag_SiO3_lim_surf_indexInMARBL(1)
 write(*,*) 'ecosys_diag_light_lim_surf_indexInMARBL(1) = ', 	ecosys_diag_light_lim_surf_indexInMARBL(1)
 write(*,*) 'ecosys_diag_photoC_indexInMARBL(1) = ', 	ecosys_diag_photoC_indexInMARBL(1)
 write(*,*) 'ecosys_diag_photoC_NO3_indexInMARBL(1) = ', 	ecosys_diag_photoC_NO3_indexInMARBL(1)
 write(*,*) 'ecosys_diag_photoFe_indexInMARBL(1) = ', 	ecosys_diag_photoFe_indexInMARBL(1)
 write(*,*) 'ecosys_diag_photoNO3_indexInMARBL(1) = ', 	ecosys_diag_photoNO3_indexInMARBL(1)
 write(*,*) 'ecosys_diag_photoNH4_indexInMARBL(1) = ', 	ecosys_diag_photoNH4_indexInMARBL(1)
 write(*,*) 'ecosys_diag_DOP_uptake_indexInMARBL(1) = ', 	ecosys_diag_DOP_uptake_indexInMARBL(1)
 write(*,*) 'ecosys_diag_PO4_uptake_indexInMARBL(1) = ', 	ecosys_diag_PO4_uptake_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_graze_indexInMARBL(1) = ', 	ecosys_diag_auto_graze_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_graze_poc_indexInMARBL(1) = ', 	ecosys_diag_auto_graze_poc_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_graze_doc_indexInMARBL(1) = ', 	ecosys_diag_auto_graze_doc_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_graze_zoo_indexInMARBL(1) = ', 	ecosys_diag_auto_graze_zoo_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_loss_indexInMARBL(1) = ', 	ecosys_diag_auto_loss_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_loss_poc_indexInMARBL(1) = ', 	ecosys_diag_auto_loss_poc_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_loss_doc_indexInMARBL(1) = ', 	ecosys_diag_auto_loss_doc_indexInMARBL(1)
 write(*,*) 'ecosys_diag_auto_agg_indexInMARBL(1) = ', 	ecosys_diag_auto_agg_indexInMARBL(1)
 write(*,*) 'ecosys_diag_Nfix_indexInMARBL(2) = ', 	ecosys_diag_Nfix_indexInMARBL(2)
 write(*,*) 'ecosys_diag_Qp_indexInMARBL(2) = ', 	ecosys_diag_Qp_indexInMARBL(2)
 write(*,*) 'ecosys_diag_N_lim_surf_indexInMARBL(2) = ', 	ecosys_diag_N_lim_surf_indexInMARBL(2)
 write(*,*) 'ecosys_diag_P_lim_surf_indexInMARBL(2) = ', 	ecosys_diag_P_lim_surf_indexInMARBL(2)
 write(*,*) 'ecosys_diag_Fe_lim_surf_indexInMARBL(2) = ', 	ecosys_diag_Fe_lim_surf_indexInMARBL(2)
 write(*,*) 'ecosys_diag_SiO3_lim_surf_indexInMARBL(2) = ', 	ecosys_diag_SiO3_lim_surf_indexInMARBL(2)
 write(*,*) 'ecosys_diag_light_lim_surf_indexInMARBL(2) = ', 	ecosys_diag_light_lim_surf_indexInMARBL(2)
 write(*,*) 'ecosys_diag_photoC_indexInMARBL(2) = ', 	ecosys_diag_photoC_indexInMARBL(2)
 write(*,*) 'ecosys_diag_photoC_NO3_indexInMARBL(2) = ', 	ecosys_diag_photoC_NO3_indexInMARBL(2)
 write(*,*) 'ecosys_diag_photoFe_indexInMARBL(2) = ', 	ecosys_diag_photoFe_indexInMARBL(2)
 write(*,*) 'ecosys_diag_photoNO3_indexInMARBL(2) = ', 	ecosys_diag_photoNO3_indexInMARBL(2)
 write(*,*) 'ecosys_diag_photoNH4_indexInMARBL(2) = ', 	ecosys_diag_photoNH4_indexInMARBL(2)
 write(*,*) 'ecosys_diag_DOP_uptake_indexInMARBL(2) = ', 	ecosys_diag_DOP_uptake_indexInMARBL(2)
 write(*,*) 'ecosys_diag_PO4_uptake_indexInMARBL(2) = ', 	ecosys_diag_PO4_uptake_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_graze_indexInMARBL(2) = ', 	ecosys_diag_auto_graze_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_graze_poc_indexInMARBL(2) = ', 	ecosys_diag_auto_graze_poc_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_graze_doc_indexInMARBL(2) = ', 	ecosys_diag_auto_graze_doc_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_graze_zoo_indexInMARBL(2) = ', 	ecosys_diag_auto_graze_zoo_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_loss_indexInMARBL(2) = ', 	ecosys_diag_auto_loss_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_loss_poc_indexInMARBL(2) = ', 	ecosys_diag_auto_loss_poc_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_loss_doc_indexInMARBL(2) = ', 	ecosys_diag_auto_loss_doc_indexInMARBL(2)
 write(*,*) 'ecosys_diag_auto_agg_indexInMARBL(2) = ', 	ecosys_diag_auto_agg_indexInMARBL(2)
 write(*,*) 'ecosys_diag_Nfix_indexInMARBL(2) = ', 	ecosys_diag_Nfix_indexInMARBL(2)
 write(*,*) 'ecosys_diag_Qp_indexInMARBL(3) = ', 	ecosys_diag_Qp_indexInMARBL(3)
 write(*,*) 'ecosys_diag_N_lim_surf_indexInMARBL(3) = ', 	ecosys_diag_N_lim_surf_indexInMARBL(3)
 write(*,*) 'ecosys_diag_P_lim_surf_indexInMARBL(3) = ', 	ecosys_diag_P_lim_surf_indexInMARBL(3)
 write(*,*) 'ecosys_diag_Fe_lim_surf_indexInMARBL(3) = ', 	ecosys_diag_Fe_lim_surf_indexInMARBL(3)
 write(*,*) 'ecosys_diag_SiO3_lim_surf_indexInMARBL(3) = ', 	ecosys_diag_SiO3_lim_surf_indexInMARBL(3)
 write(*,*) 'ecosys_diag_light_lim_surf_indexInMARBL(3) = ', 	ecosys_diag_light_lim_surf_indexInMARBL(3)
 write(*,*) 'ecosys_diag_photoC_indexInMARBL(3) = ', 	ecosys_diag_photoC_indexInMARBL(3)
 write(*,*) 'ecosys_diag_photoC_NO3_indexInMARBL(3) = ', 	ecosys_diag_photoC_NO3_indexInMARBL(3)
 write(*,*) 'ecosys_diag_photoFe_indexInMARBL(3) = ', 	ecosys_diag_photoFe_indexInMARBL(3)
 write(*,*) 'ecosys_diag_photoNO3_indexInMARBL(3) = ', 	ecosys_diag_photoNO3_indexInMARBL(3)
 write(*,*) 'ecosys_diag_photoNH4_indexInMARBL(3) = ', 	ecosys_diag_photoNH4_indexInMARBL(3)
 write(*,*) 'ecosys_diag_DOP_uptake_indexInMARBL(3) = ', 	ecosys_diag_DOP_uptake_indexInMARBL(3)
 write(*,*) 'ecosys_diag_PO4_uptake_indexInMARBL(3) = ', 	ecosys_diag_PO4_uptake_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_graze_indexInMARBL(3) = ', 	ecosys_diag_auto_graze_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_graze_poc_indexInMARBL(3) = ', 	ecosys_diag_auto_graze_poc_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_graze_doc_indexInMARBL(3) = ', 	ecosys_diag_auto_graze_doc_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_graze_zoo_indexInMARBL(3) = ', 	ecosys_diag_auto_graze_zoo_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_loss_indexInMARBL(3) = ', 	ecosys_diag_auto_loss_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_loss_poc_indexInMARBL(3) = ', 	ecosys_diag_auto_loss_poc_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_loss_doc_indexInMARBL(3) = ', 	ecosys_diag_auto_loss_doc_indexInMARBL(3)
 write(*,*) 'ecosys_diag_auto_agg_indexInMARBL(3) = ', 	ecosys_diag_auto_agg_indexInMARBL(3)
 write(*,*) 'ecosys_diag_Nfix_indexInMARBL(3) = ', 	ecosys_diag_Nfix_indexInMARBL(3)
	
 write(*,*) 'ecosys_diag_tot_bSi_form_indexInMARBL = ', 	ecosys_diag_tot_bSi_form_indexInMARBL
	
 write(*,*) 'ecosys_diag_SiO2_FLUX_IN_indexInMARBL = ', 	ecosys_diag_SiO2_FLUX_IN_indexInMARBL
 write(*,*) 'ecosys_diag_SiO2_PROD_indexInMARBL = ', 	ecosys_diag_SiO2_PROD_indexInMARBL
 write(*,*) 'ecosys_diag_SiO2_REMIN_indexInMARBL = ', 	ecosys_diag_SiO2_REMIN_indexInMARBL
 write(*,*) 'ecosys_diag_dust_FLUX_IN_indexInMARBL = ', 	ecosys_diag_dust_FLUX_IN_indexInMARBL
 write(*,*) 'ecosys_diag_dust_REMIN_indexInMARBL = ', 	ecosys_diag_dust_REMIN_indexInMARBL
 write(*,*) 'ecosys_diag_P_iron_FLUX_IN_indexInMARBL = ', 	ecosys_diag_P_iron_FLUX_IN_indexInMARBL
 write(*,*) 'ecosys_diag_P_iron_PROD_indexInMARBL = ', 	ecosys_diag_P_iron_PROD_indexInMARBL
 write(*,*) 'ecosys_diag_P_iron_REMIN_indexInMARBL = ', 	ecosys_diag_P_iron_REMIN_indexInMARBL
 write(*,*) 'ecosys_diag_DOC_prod_indexInMARBL = ', 	ecosys_diag_DOC_prod_indexInMARBL
 write(*,*) 'ecosys_diag_DOC_remin_indexInMARBL = ', 	ecosys_diag_DOC_remin_indexInMARBL
 write(*,*) 'ecosys_diag_DOCr_remin_indexInMARBL = ', 	ecosys_diag_DOCr_remin_indexInMARBL
 write(*,*) 'ecosys_diag_DON_prod_indexInMARBL = ', 	ecosys_diag_DON_prod_indexInMARBL
 write(*,*) 'ecosys_diag_DON_remin_indexInMARBL = ', 	ecosys_diag_DON_remin_indexInMARBL
 write(*,*) 'ecosys_diag_DOP_prod_indexInMARBL = ', 	ecosys_diag_DOP_prod_indexInMARBL
 write(*,*) 'ecosys_diag_DOP_remin_indexInMARBL = ', 	ecosys_diag_DOP_remin_indexInMARBL
 write(*,*) 'ecosys_diag_Fe_scavenge_indexInMARBL = ', 	ecosys_diag_Fe_scavenge_indexInMARBL
 write(*,*) 'ecosys_diag_Fe_scavenge_rate_indexInMARBL = ', 	ecosys_diag_Fe_scavenge_rate_indexInMARBL
 write(*,*) 'ecosys_diag_Lig_prod_indexInMARBL = ', 	ecosys_diag_Lig_prod_indexInMARBL
 write(*,*) 'ecosys_diag_Lig_loss_indexInMARBL = ', 	ecosys_diag_Lig_loss_indexInMARBL
 write(*,*) 'ecosys_diag_Lig_scavenge_indexInMARBL = ', 	ecosys_diag_Lig_scavenge_indexInMARBL
 write(*,*) 'ecosys_diag_Fefree_indexInMARBL = ', 	ecosys_diag_Fefree_indexInMARBL
 write(*,*) 'ecosys_diag_Lig_photochem_indexInMARBL = ', 	ecosys_diag_Lig_photochem_indexInMARBL
 write(*,*) 'ecosys_diag_Lig_deg_indexInMARBL = ', 	ecosys_diag_Lig_deg_indexInMARBL
 write(*,*) 'ecosys_diag_DONr_remin_indexInMARBL = ', 	ecosys_diag_DONr_remin_indexInMARBL
 write(*,*) 'ecosys_diag_DOPr_remin_indexInMARBL = ', 	ecosys_diag_DOPr_remin_indexInMARBL
 write(*,*) 'ecosys_diag_ponToSed_indexInMARBL = ', 	ecosys_diag_ponToSed_indexInMARBL
 write(*,*) 'ecosys_diag_popToSed_indexInMARBL = ', 	ecosys_diag_popToSed_indexInMARBL
 write(*,*) 'ecosys_diag_bsiToSed_indexInMARBL = ', 	ecosys_diag_bsiToSed_indexInMARBL
 write(*,*) 'ecosys_diag_dustToSed_indexInMARBL = ', 	ecosys_diag_dustToSed_indexInMARBL
 write(*,*) 'ecosys_diag_OtherRemin_indexInMARBL = ', 	ecosys_diag_OtherRemin_indexInMARBL
	
 write(*,*) 'ecosys_diag_bSi_form_indexInMARBL(1) = ', 	ecosys_diag_bSi_form_indexInMARBL(1)
 write(*,*) 'ecosys_diag_bSi_form_indexInMARBL(2) = ', 	ecosys_diag_bSi_form_indexInMARBL(2)
 write(*,*) 'ecosys_diag_bSi_form_indexInMARBL(3) = ', 	ecosys_diag_bSi_form_indexInMARBL(3)
	
 write(*,*) 'ecosys_diag_tot_CaCO3_form_zint_indexInMARBL = ', 	ecosys_diag_tot_CaCO3_form_zint_indexInMARBL
 write(*,*) 'ecosys_diag_zsatcalc_indexInMARBL = ', 	ecosys_diag_zsatcalc_indexInMARBL
 write(*,*) 'ecosys_diag_zsatarag_indexInMARBL = ', 	ecosys_diag_zsatarag_indexInMARBL
	
 write(*,*) 'ecosys_diag_CaCO3_form_zint_indexInMARBL(1) = ', 	ecosys_diag_CaCO3_form_zint_indexInMARBL(1)
 write(*,*) 'ecosys_diag_CaCO3_form_zint_indexInMARBL(2) = ', 	ecosys_diag_CaCO3_form_zint_indexInMARBL(2)
 write(*,*) 'ecosys_diag_CaCO3_form_zint_indexInMARBL(3) = ', 	ecosys_diag_CaCO3_form_zint_indexInMARBL(3)
	
 write(*,*) 'ecosys_diag_CaCO3_PROD_indexInMARBL = ', 	ecosys_diag_CaCO3_PROD_indexInMARBL
 write(*,*) 'ecosys_diag_CaCO3_REMIN_indexInMARBL = ', 	ecosys_diag_CaCO3_REMIN_indexInMARBL
 write(*,*) 'ecosys_diag_CaCO3_ALT_CO2_PROD_indexInMARBL = ', 	ecosys_diag_CaCO3_ALT_CO2_PROD_indexInMARBL
 write(*,*) 'ecosys_diag_CaCO3_ALT_CO2_REMIN_indexInMARBL = ', 	ecosys_diag_CaCO3_ALT_CO2_REMIN_indexInMARBL
 write(*,*) 'ecosys_diag_tot_CaCO3_form_indexInMARBL = ', 	ecosys_diag_tot_CaCO3_form_indexInMARBL
 write(*,*) 'ecosys_diag_CO3_indexInMARBL = ', 	ecosys_diag_CO3_indexInMARBL
 write(*,*) 'ecosys_diag_HCO3_indexInMARBL = ', 	ecosys_diag_HCO3_indexInMARBL
 write(*,*) 'ecosys_diag_H2CO3_indexInMARBL = ', 	ecosys_diag_H2CO3_indexInMARBL
 write(*,*) 'ecosys_diag_CO3_ALT_CO2_indexInMARBL = ', 	ecosys_diag_CO3_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_HCO3_ALT_CO2_indexInMARBL = ', 	ecosys_diag_HCO3_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_H2CO3_ALT_CO2_indexInMARBL = ', 	ecosys_diag_H2CO3_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_pH_3D_ALT_CO2_indexInMARBL = ', 	ecosys_diag_pH_3D_ALT_CO2_indexInMARBL
 write(*,*) 'ecosys_diag_co3_sat_calc_indexInMARBL = ', 	ecosys_diag_co3_sat_calc_indexInMARBL
 write(*,*) 'ecosys_diag_co3_sat_arag_indexInMARBL = ', 	ecosys_diag_co3_sat_arag_indexInMARBL
	
 write(*,*) 'ecosys_diag_CaCO3_form_indexInMARBL(1) = ', 	ecosys_diag_CaCO3_form_indexInMARBL(1)
 write(*,*) 'ecosys_diag_CaCO3_form_indexInMARBL(2) = ', 	ecosys_diag_CaCO3_form_indexInMARBL(2)
 write(*,*) 'ecosys_diag_CaCO3_form_indexInMARBL(3) = ', 	ecosys_diag_CaCO3_form_indexInMARBL(3)

    end subroutine ocn_ecosys_echo_indices!}}}

end module ocn_tracer_ecosys

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
