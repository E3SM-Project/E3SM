module clm_bgc_interface_data
!!=================================================================================================
! CLM BioGeoChemistry (BGC) Interface: Data Type (Variables)
!
! Created by wgs @ ORNL
!
! date: 8/25/2015
!!=================================================================================================
  !! USES:
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar            , only : nlevsno, nlevgrnd
  use clm_varpar            , only : nlevdecomp_full, ndecomp_pools,  ndecomp_cascade_transitions
  use clm_varcon            , only : spval
  use decompMod             , only : bounds_type

  implicit none
!  save
  private

  type, public :: clm_bgc_interface_data_type

     ! clm_varpar
     integer                    :: nlevdecomp                               ! num of CLM soil layers that are mapped to/from PFLOTRAN
     integer                    :: ndecomp_pools                            ! num of decomposition pools

     ! decomp_cascade_con
     logical, pointer           :: floating_cn_ratio                (:)     ! TRUE => pool has fixed C:N ratio
     logical, pointer           :: floating_cp_ratio                (:)     ! TRUE => pool has fixed C:P ratio
     character(len=8), pointer  :: decomp_pool_name                 (:)     ! name of pool
     real(r8), pointer          :: initial_cn_ratio                 (:)     ! c:n ratio for initialization of pools
     real(r8), pointer          :: initial_cp_ratio                 (:)     ! c:n ratio for initialization of pools

     ! col:
     real(r8), pointer :: z                                         (:,:)   ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: dz                                        (:,:)   ! layer thickness (m)  (-nlevsno+1:nlevgrnd)

     ! soilstate_vars:
     real(r8), pointer :: bd_col                                    (:,:)   ! col bulk density of dry soil material [kg/m^3] (CN)
     real(r8), pointer :: hksat_col                                 (:,:)   ! col hydraulic conductivity at saturation (mm H2O /s)
     real(r8), pointer :: bsw_col                                   (:,:)   ! col Clapp and Hornberger "b" (nlevgrnd)
     real(r8), pointer :: watsat_col                                (:,:)   ! col volumetric soil water at saturation (porosity)
     real(r8), pointer :: sucsat_col                                (:,:)   ! col minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: watfc_col                                 (:,:)   ! col volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: porosity_col                              (:,:)   ! col soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity_col                          (:,:)   ! col effective porosity = porosity - vol_ice (nlevgrnd)
     real(r8), pointer :: cellorg_col                               (:,:)   ! col organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: soilpsi_col                               (:,:)   ! col soil water potential in each soil layer (MPa) (CN)
     real(r8), pointer :: rootfr_col                                (:,:)   ! col fraction of roots in each soil layer (nlevgrnd)

     ! waterstate_vars:
     real(r8), pointer :: h2osoi_liq_col                            (:,:)   ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: h2osoi_ice_col                            (:,:)   ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: frac_sno_eff_col                          (:)     ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_h2osfc_col                           (:)     ! col fractional area with surface water greater than zero
     real(r8), pointer :: h2osoi_vol_col                            (:,:)   ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)

     ! temperature_vars:
     real(r8), pointer :: t_soisno_col                              (:,:)   ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: t_grnd_col                                (:)     ! col ground temperature (Kelvin)

     ! canopystate_vars
     integer  , pointer :: alt_indx_col                             (:)     ! col current depth of thaw

     ! ch4
     real(r8), pointer  :: finundated_col                           (:)     ! col fractional inundated area (excluding dedicated wetland cols)
     real(r8), pointer  :: o2stress_unsat_col                       (:,:)   ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer  :: o2stress_sat_col                         (:,:)   ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer  :: conc_o2_sat_col                          (:,:)   ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer  :: conc_o2_unsat_col                        (:,:)   ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer  :: o2_decomp_depth_sat_col                  (:,:)   ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer  :: o2_decomp_depth_unsat_col                (:,:)   ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)

     ! carbonstate_vars:
     real(r8), pointer :: decomp_cpools_vr_col                      (:,:,:) ! col (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools

     ! nitrogenstate_vars:
     real(r8), pointer :: decomp_npools_vr_col                      (:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: sminn_vr_col                              (:,:)   ! col (gN/m3) vertically-resolved soil mineral N
     real(r8), pointer :: smin_no3_vr_col                           (:,:)   ! col (gN/m3) vertically-resolved soil mineral NO3
     real(r8), pointer :: smin_nh4_vr_col                           (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4
     real(r8), pointer :: smin_nh4sorb_vr_col                       (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4 absorbed

     ! phosphorusstate_vars:
     real(r8), pointer :: decomp_ppools_vr_col                      (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
     real(r8), pointer :: solutionp_vr_col                          (:,:)   ! col (gP/m3) vertically-resolved soil solution P
     real(r8), pointer :: labilep_vr_col                            (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
     real(r8), pointer :: secondp_vr_col                            (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
     real(r8), pointer :: occlp_vr_col                              (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
     real(r8), pointer :: primp_vr_col                              (:,:)   ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: sminp_vr_col                              (:,:)   ! col (gP/m3) vertically-resolved soil parimary mineral P


     ! plant NP demand
     real(r8), pointer :: plant_ndemand_col                         (:)     ! col N flux required to support initial GPP (gN/m3/s)
     real(r8), pointer :: plant_pdemand_col                         (:)     ! col P flux required to support initial GPP (gP/m3/s)
     real(r8), pointer :: plant_ndemand_vr_col                      (:,:)   ! col vertically-resolved N flux required to support initial GPP (gN/m3/s)
     real(r8), pointer :: plant_pdemand_vr_col                      (:,:)   ! col vertically-resolved P flux required to support initial GPP (gP/m3/s)

     ! decomposition flux:
     real(r8), pointer :: decomp_cpools_sourcesink_col              (:,:,:) ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)
     real(r8), pointer :: decomp_npools_sourcesink_col              (:,:,:) ! col (gN/m3) change in decomposing n pools
     real(r8), pointer :: decomp_ppools_sourcesink_col              (:,:,:) ! col (gP/m3) change in decomposing n pools

     real(r8), pointer :: decomp_cascade_ctransfer_vr_col           (:,:,:) ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
     real(r8), pointer :: decomp_cascade_ntransfer_vr_col           (:,:,:) ! col vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_ptransfer_vr_col           (:,:,:) ! col vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_vr_col          (:,:,:) ! col vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_sminp_flux_vr_col          (:,:,:) ! col vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
     real(r8), pointer :: sminn_to_denit_decomp_cascade_vr_col      (:,:,:) ! col vertically-resolved denitrification along decomp cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_hr_vr_col                  (:,:,:) ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)

     ! mineralization / immobilization / uptake flux
     real(r8), pointer :: gross_nmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of N mineralization (gN/m3/s)
     real(r8), pointer :: net_nmin_vr_col                           (:,:)   ! col vertically-resolved net rate of N mineralization (gN/m3/s)

     real(r8), pointer :: potential_immob_vr_col                    (:,:)   ! col vertically-resolved potential N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_vr_col                       (:,:)   ! col vertically-resolved actual N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_no3_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
     real(r8), pointer :: actual_immob_nh4_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)
     real(r8), pointer :: supplement_to_sminn_vr_col                (:,:)   ! col vertically-resolved supplemental N supply (gN/m3/s)

     real(r8), pointer :: sminn_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral N (gN/m3/s)
     real(r8), pointer :: smin_no3_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
     real(r8), pointer :: smin_nh4_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)

     real(r8), pointer :: potential_immob_col                       (:)     ! col vert-int (diagnostic) potential N immobilization (gN/m2/s)
     real(r8), pointer :: actual_immob_col                          (:)     ! col vert-int (diagnostic) actual N immobilization (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_col                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)

     real(r8), pointer :: potential_immob_p_vr_col                  (:,:)   ! col vertically-resolved potential P immobilization (gP/m3/s) at each level
     real(r8), pointer :: actual_immob_p_vr_col                     (:,:)   ! col vertically-resolved actual P immobilization (gP/m3/s) at each level
     real(r8), pointer :: sminp_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral P (gP/m3/s)
     real(r8), pointer :: supplement_to_sminp_vr_col                (:,:)   ! col vertically-resolved supplemental P supply (gP/m3/s)
     
     real(r8), pointer :: gross_pmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of P mineralization (gP/m3/s)
     real(r8), pointer :: net_pmin_vr_col                           (:,:)   ! col vertically-resolved net rate of P mineralization (gP/m3/s)

     real(r8), pointer :: potential_immob_p_col                     (:)     ! col vert-int (diagnostic) potential P immobilization (gP/m2/s)
     real(r8), pointer :: actual_immob_p_col                        (:)     ! col vert-int (diagnostic) actual P immobilization (gP/m2/s)
     real(r8), pointer :: sminp_to_plant_col                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral P (gP/m2/s)
     
     ! nitrification / denitrification flux:
     real(r8), pointer :: f_nit_vr_col                              (:,:)   ! col (gN/m3/s) soil nitrification flux
     real(r8), pointer :: f_denit_vr_col                            (:,:)   ! col (gN/m3/s) soil denitrification flux
     real(r8), pointer :: pot_f_nit_vr_col                          (:,:)   ! col (gN/m3/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_vr_col                        (:,:)   ! col (gN/m3/s) potential soil denitrification flux
     real(r8), pointer :: n2_n2o_ratio_denit_vr_col                 (:,:)   ! col ratio of N2 to N2O production by denitrification [gN/gN]
     real(r8), pointer :: f_n2o_denit_vr_col                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_nit_vr_col                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]

     real(r8), pointer :: sminn_to_denit_excess_vr_col              (:,:)   ! col vertically-resolved denitrification from excess mineral N pool (gN/m3/s)

     ! inorganic P transformation
     real(r8), pointer :: primp_to_labilep_vr_col                   (:,:)   ! col (gP/m3/s) flux of P from primary mineral to labile
     real(r8), pointer :: labilep_to_secondp_vr_col                 (:,:)   ! col (gP/m3/s) flux of labile P to secondary mineral P
     real(r8), pointer :: secondp_to_labilep_vr_col                 (:,:)   ! col (gP/m3/s) flux of the desorption of secondary mineral P to labile P
     real(r8), pointer :: secondp_to_occlp_vr_col                   (:,:)   ! col (gP/m3/s) flux of the occlusion of secondary P to occluded P

     ! gases
     real(r8), pointer :: hr_vr_col                                 (:,:)   ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
     
     real(r8), pointer :: f_co2_soil_vr_col                         (:,:)   ! total vertically-resolved soil-atm. CO2 exchange (gC/m3/s)
     real(r8), pointer :: f_n2o_soil_vr_col                         (:,:)   ! col flux of N2o from soil-N processes [gN/m^3/s]
     real(r8), pointer :: f_n2_soil_vr_col                          (:,:)   ! col flux of N2 from soil-N processes [gN/m^3/s]

     real(r8), pointer :: phr_vr_col                                (:,:)   ! potential hr (not N-limited) (gC/m3/s)
     real(r8), pointer :: fphr_col                                  (:,:)   ! fraction of potential heterotrophic respiration

     ! fpi, fpg
     real(r8) , pointer :: fpi_vr_col                               (:,:)   ! col fraction of potential immobilization (no units)
     real(r8) , pointer :: fpi_col                                  (:)     ! col fraction of potential immobilization (no units)
     real(r8),  pointer :: fpg_col                                  (:)     ! col fraction of potential gpp (no units)

     real(r8) , pointer :: fpi_p_vr_col                             (:,:)   ! col fraction of potential immobilization (no units)
     real(r8) , pointer :: fpi_p_col                                (:)     ! col fraction of potential immobilization (no units)
     real(r8),  pointer :: fpg_p_col                                (:)     ! col fraction of potential gpp (no units)

     !!------------------------------------------------------------------------------------------
     !! pflotran variables: BEGIN
     !!------------------------------------------------------------------------------------------
     ! bgc rates/fluxes (previous time-step) to decomposition pools
     real(r8), pointer :: externalc_to_decomp_cpools_col            (:,:,:) ! col (gC/m3/s) net C fluxes associated with litter/som-adding/removal to decomp pools
     real(r8), pointer :: externaln_to_decomp_npools_col            (:,:,:) ! col (gN/m3/s) net N fluxes associated with litter/som-adding/removal to decomp pools
     real(r8), pointer :: externalp_to_decomp_ppools_col            (:,:,:) ! col (gP/m3/s) net P fluxes associated with litter/som-adding/removal to decomp pools

     ! bgc rates/fluxes (previous time-step) to nh4 / no3
     real(r8), pointer :: externaln_to_nh4_col                      (:,:)   ! col (gN/m3/s) net N fluxes to nh4 pool: deposition + fertilization + supplement + nfix + soyfixn
     real(r8), pointer :: externaln_to_no3_col                      (:,:)   ! col (gN/m3/s) net N fluxes to no3 pool: deposition + fertilization + supplement
     real(r8), pointer :: externaln_to_sminn_col                    (:,:)   ! col (gN/m3/s) net N fluxes to sminn pool: deposition + fertilization + supplement + nfix + soyfixn
     real(r8), pointer :: smin_no3_leached_vr_col                   (:,:)   ! col vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
     real(r8), pointer :: smin_no3_runoff_vr_col                    (:,:)   ! col vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
     real(r8), pointer :: no3_net_transport_vr_col                  (:,:)   ! col net NO3 transport associated with runoff/leaching (gN/m3/s)

     ! bgc rates/fluxes (previous time-step) to mineral P
     real(r8), pointer :: externalp_to_primp_col                    (:,:)   ! pdep_to_sminp_col                         (:)     ! col atmospheric P deposition to soil mineral P (gP/m2/s)
     real(r8), pointer :: externalp_to_labilep_col                  (:,:)   ! fert_p_to_sminp_col                       (:)     ! col fertilizer P to soil mineral P (gP/m2/s)
     real(r8), pointer :: externalp_to_solutionp_col                (:,:)   ! supplement_to_sminp_vr_col                (:,:)   ! col vertically-resolved supplemental P supply (gP/m3/s)
     real(r8), pointer :: sminp_leached_vr_col                      (:,:)   ! col vertically-resolved soil mineral P pool loss to leaching (gP/m3/s)
     real(r8), pointer :: sminp_net_transport_vr_col                (:,:)   ! sminp_leached_vr_col                      (:,:) col net sminp transport associated with runoff/leaching (gN/m3/s)

     ! gases:
     real(r8), pointer :: f_ngas_decomp_vr_col                      (:,:)   ! col vertically-resolved N emission from excess mineral N pool due to mineralization (gN/m3/s)
     real(r8), pointer :: f_ngas_nitri_vr_col                       (:,:)   ! col vertically-resolved N emission from nitrification (gN/m3/s)
     real(r8), pointer :: f_ngas_denit_vr_col                       (:,:)   ! col vertically-resolved N emission from denitrification (gN/m3/s)

     ! waterflux_vars:
     real(r8), pointer :: qflx_top_soil_col                         (:)     ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_sub_snow_col                         (:)     ! col sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_evap_soi_col                         (:)     ! col soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_ev_h2osfc_col                        (:)     ! col evaporation heat flux from soil         (W/m**2) [+ to atm]
     real(r8), pointer :: qflx_tran_veg_col                         (:)     ! col vegetation transpiration (mm H2O/s) (+ = to atm)

     ! energyflux_vars:
     real(r8), pointer :: htvp_col                                  (:)     ! latent heat of vapor of water (or sublimation) [j/kg]
     real(r8), pointer :: eflx_bot_col                              (:)     ! col heat flux from beneath the soil or ice column (W/m**2)
     real(r8), pointer :: eflx_gnet_col                             (:)     ! col net heat flux into ground  (W/m**2)
     real(r8), pointer :: eflx_soil_grnd_col                        (:)     ! col soil heat flux (W/m**2) [+ = into soil]

     ! atm2lnd:
     real(r8), pointer :: forc_pbot_not_downscaled_grc              (:)     ! not downscaled atm pressure (Pa)
     real(r8), pointer :: forc_pco2_grc                             (:)     ! CO2 partial pressure (Pa)
     real(r8), pointer :: forc_pch4_grc                             (:)     ! CH4 partial pressure (Pa)

     !!------------------------------------------------------------------------------------------
     !! pflotran variables: END
     !!------------------------------------------------------------------------------------------


!!-------------------------------------------------------------------------------------------------
  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type clm_bgc_interface_data_type
!!-------------------------------------------------------------------------------------------------

!!  type(clm_bgc_interface_data_type) , public, target , save :: clm_bgc_data
  
contains


!!-------------------------------------------------------------------------------------------------
  subroutine Init(this, bounds)
     class(clm_bgc_interface_data_type) :: this
     type(bounds_type), intent(in)      :: bounds

     call this%InitAllocate (bounds)
  end subroutine Init
!!-------------------------------------------------------------------------------------------------

!!-------------------------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
    !! USES


    !! ARGUMENTS:
    real(r8) :: ival  = 0.0_r8  ! initial value
    class(clm_bgc_interface_data_type) :: this
    type(bounds_type), intent(in) :: bounds

    !! LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !------------------------------------------------------------------------
    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    ! decomp_cascade_type
    allocate(this%decomp_pool_name      (1:ndecomp_pools))                  ; this%decomp_pool_name         (:)   = ''
    allocate(this%floating_cn_ratio     (1:ndecomp_pools))                  ; this%floating_cn_ratio        (:)   = .false.
    allocate(this%floating_cp_ratio     (1:ndecomp_pools))                  ; this%floating_cp_ratio        (:)   = .false.
    allocate(this%initial_cn_ratio      (0:ndecomp_pools))                  ; this%initial_cn_ratio         (:)   = nan
    allocate(this%initial_cp_ratio      (0:ndecomp_pools))                  ; this%initial_cp_ratio         (:)   = nan

    ! col:
    allocate(this%z                     (begc:endc,-nlevsno+1:nlevgrnd))    ; this%z                    (:,:) = nan
    allocate(this%dz                    (begc:endc,-nlevsno+1:nlevgrnd))    ; this%dz                   (:,:) = nan
    ! soilstate_vars:
    allocate(this%bd_col                (begc:endc,nlevgrnd))               ; this%bd_col               (:,:) = nan
    allocate(this%hksat_col             (begc:endc,nlevgrnd))               ; this%hksat_col            (:,:) = spval
    allocate(this%bsw_col               (begc:endc,nlevgrnd))               ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col            (begc:endc,nlevgrnd))               ; this%watsat_col           (:,:) = nan
    allocate(this%sucsat_col            (begc:endc,nlevgrnd))               ; this%sucsat_col           (:,:) = spval
    allocate(this%watfc_col             (begc:endc,nlevgrnd))               ; this%watfc_col            (:,:) = nan
    allocate(this%porosity_col          (begc:endc,nlevgrnd))               ; this%porosity_col         (:,:) = spval
    allocate(this%eff_porosity_col      (begc:endc,nlevgrnd))               ; this%eff_porosity_col     (:,:) = spval
    allocate(this%cellorg_col           (begc:endc,nlevgrnd))               ; this%cellorg_col          (:,:) = nan
    allocate(this%soilpsi_col           (begc:endc,nlevgrnd))               ; this%soilpsi_col          (:,:) = nan
    allocate(this%rootfr_col            (begc:endc,1:nlevgrnd))             ; this%rootfr_col           (:,:) = nan

    ! waterstate_vars:
    allocate(this%h2osoi_liq_col        (begc:endc,-nlevsno+1:nlevgrnd))    ; this%h2osoi_liq_col       (:,:) = nan
    allocate(this%h2osoi_ice_col        (begc:endc,-nlevsno+1:nlevgrnd))    ; this%h2osoi_ice_col       (:,:) = nan
    allocate(this%frac_sno_eff_col      (begc:endc))                        ; this%frac_sno_eff_col     (:)   = nan
    allocate(this%frac_h2osfc_col       (begc:endc))                        ; this%frac_h2osfc_col      (:)   = nan
    allocate(this%h2osoi_vol_col        (begc:endc, 1:nlevgrnd))            ; this%h2osoi_vol_col       (:,:) = nan

    ! temperature_vars:
    allocate(this%t_soisno_col          (begc:endc,-nlevsno+1:nlevgrnd))    ; this%t_soisno_col         (:,:) = nan
    allocate(this%t_grnd_col            (begc:endc))                        ; this%t_grnd_col           (:)   = nan

    ! canopystate_vars
    allocate(this%alt_indx_col                  (begc:endc))                ; this%alt_indx_col                 (:)   = huge(1)

    ! ch4
    allocate(this%finundated_col                (begc:endc))                ; this%finundated_col               (:)   = nan
    allocate(this%o2stress_unsat_col            (begc:endc,1:nlevgrnd))     ; this%o2stress_unsat_col           (:,:) = nan
    allocate(this%o2stress_sat_col              (begc:endc,1:nlevgrnd))     ; this%o2stress_sat_col             (:,:) = nan
    allocate(this%conc_o2_sat_col               (begc:endc,1:nlevgrnd))     ; this%conc_o2_sat_col              (:,:) = nan
    allocate(this%conc_o2_unsat_col             (begc:endc,1:nlevgrnd))     ; this%conc_o2_unsat_col            (:,:) = nan
    allocate(this%o2_decomp_depth_sat_col       (begc:endc,1:nlevgrnd))     ; this%o2_decomp_depth_sat_col      (:,:) = nan
    allocate(this%o2_decomp_depth_unsat_col     (begc:endc,1:nlevgrnd))     ; this%o2_decomp_depth_unsat_col    (:,:) = nan

    ! carbonstate_vars:
    allocate(this%decomp_cpools_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_cpools_vr_col(:,:,:)= ival
    allocate(this%decomp_npools_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_npools_vr_col(:,:,:)= ival
    allocate(this%decomp_ppools_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_ppools_vr_col(:,:,:)= ival

    ! nitrogenstate_vars:
    allocate(this%sminn_vr_col          (begc:endc,1:nlevdecomp_full))      ; this%sminn_vr_col             (:,:) = ival
    allocate(this%smin_no3_vr_col       (begc:endc,1:nlevdecomp_full))      ; this%smin_no3_vr_col          (:,:) = ival
    allocate(this%smin_nh4_vr_col       (begc:endc,1:nlevdecomp_full))      ; this%smin_nh4_vr_col          (:,:) = ival
    allocate(this%smin_nh4sorb_vr_col   (begc:endc,1:nlevdecomp_full))      ; this%smin_nh4sorb_vr_col      (:,:) = ival

    ! phosphorusstate_vars:
    allocate(this%solutionp_vr_col      (begc:endc,1:nlevdecomp_full))      ; this%solutionp_vr_col         (:,:) = ival
    allocate(this%labilep_vr_col        (begc:endc,1:nlevdecomp_full))      ; this%labilep_vr_col           (:,:) = ival
    allocate(this%secondp_vr_col        (begc:endc,1:nlevdecomp_full))      ; this%secondp_vr_col           (:,:) = ival
    allocate(this%occlp_vr_col          (begc:endc,1:nlevdecomp_full))      ; this%occlp_vr_col             (:,:) = ival
    allocate(this%primp_vr_col          (begc:endc,1:nlevdecomp_full))      ; this%primp_vr_col             (:,:) = ival
    allocate(this%sminp_vr_col          (begc:endc,1:nlevdecomp_full))      ; this%sminp_vr_col             (:,:) = ival

    ! plant N/P demand
    allocate(this%plant_ndemand_col         (begc:endc))                    ; this%plant_ndemand_col            (:)    = ival
    allocate(this%plant_pdemand_col         (begc:endc))                    ; this%plant_pdemand_col            (:)    = ival
    allocate(this%plant_ndemand_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%plant_ndemand_vr_col         (:,:)  = ival
    allocate(this%plant_pdemand_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%plant_pdemand_vr_col         (:,:)  = ival

    ! decomposition flux:
    allocate(this%decomp_cpools_sourcesink_col      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))              ; this%decomp_cpools_sourcesink_col     (:,:,:) = ival
    allocate(this%decomp_npools_sourcesink_col      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))              ; this%decomp_npools_sourcesink_col     (:,:,:) = ival
    allocate(this%decomp_ppools_sourcesink_col      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))              ; this%decomp_ppools_sourcesink_col     (:,:,:) = ival

    allocate(this%decomp_cascade_ctransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_ctransfer_vr_col  (:,:,:) = ival
    allocate(this%decomp_cascade_ntransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_ntransfer_vr_col  (:,:,:) = ival
    allocate(this%decomp_cascade_ptransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_ptransfer_vr_col  (:,:,:) = ival

    allocate(this%decomp_cascade_sminn_flux_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_sminn_flux_vr_col (:,:,:) = ival
    allocate(this%decomp_cascade_sminp_flux_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_sminp_flux_vr_col (:,:,:) = ival
    allocate(this%sminn_to_denit_decomp_cascade_vr_col (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions )); this%sminn_to_denit_decomp_cascade_vr_col (:,:,:) = ival

    allocate(this%decomp_cascade_hr_vr_col          (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_hr_vr_col         (:,:,:) = ival

    ! mineralization / immobilization / uptake fluxes
    allocate(this%gross_nmin_vr_col         (begc:endc,1:nlevdecomp_full))  ; this%gross_nmin_vr_col            (:,:) = ival
    allocate(this%net_nmin_vr_col           (begc:endc,1:nlevdecomp_full))  ; this%net_nmin_vr_col              (:,:) = ival

    allocate(this%potential_immob_vr_col    (begc:endc,1:nlevdecomp_full))  ; this%potential_immob_vr_col       (:,:) = ival
    allocate(this%actual_immob_vr_col       (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_vr_col          (:,:) = ival
    allocate(this%actual_immob_no3_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_no3_vr_col      (:,:) = ival
    allocate(this%actual_immob_nh4_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_nh4_vr_col      (:,:) = ival

    allocate(this%sminn_to_plant_vr_col     (begc:endc,1:nlevdecomp_full))  ; this%sminn_to_plant_vr_col        (:,:) = ival
    allocate(this%smin_no3_to_plant_vr_col  (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_to_plant_vr_col     (:,:) = ival
    allocate(this%smin_nh4_to_plant_vr_col  (begc:endc,1:nlevdecomp_full))  ; this%smin_nh4_to_plant_vr_col     (:,:) = ival
    allocate(this%supplement_to_sminn_vr_col(begc:endc,1:nlevdecomp_full))  ; this%supplement_to_sminn_vr_col   (:,:) = ival

    allocate(this%sminn_to_plant_col        (begc:endc))                    ; this%sminn_to_plant_col           (:)   = ival
    allocate(this%potential_immob_col       (begc:endc))                    ; this%potential_immob_col          (:)   = ival
    allocate(this%actual_immob_col          (begc:endc))                    ; this%actual_immob_col             (:)   = ival

    allocate(this%gross_pmin_vr_col         (begc:endc,1:nlevdecomp_full))  ; this%gross_pmin_vr_col            (:,:) = ival
    allocate(this%net_pmin_vr_col           (begc:endc,1:nlevdecomp_full))  ; this%net_pmin_vr_col              (:,:) = ival

    allocate(this%potential_immob_p_vr_col  (begc:endc,1:nlevdecomp_full))  ; this%potential_immob_p_vr_col     (:,:) = ival
    allocate(this%actual_immob_p_vr_col     (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_p_vr_col        (:,:) = ival
    allocate(this%sminp_to_plant_vr_col     (begc:endc,1:nlevdecomp_full))  ; this%sminp_to_plant_vr_col        (:,:) = ival
    allocate(this%supplement_to_sminp_vr_col(begc:endc,1:nlevdecomp_full))  ; this%supplement_to_sminp_vr_col   (:,:) = ival

    allocate(this%sminp_to_plant_col        (begc:endc))                    ; this%sminp_to_plant_col           (:)   = ival
    allocate(this%potential_immob_p_col     (begc:endc))                    ; this%potential_immob_p_col        (:)   = ival
    allocate(this%actual_immob_p_col        (begc:endc))                    ; this%actual_immob_p_col           (:)   = ival

    ! nitrification / denitrification flux
    allocate(this%f_nit_vr_col              (begc:endc,1:nlevdecomp_full))  ; this%f_nit_vr_col                 (:,:) = ival
    allocate(this%f_denit_vr_col            (begc:endc,1:nlevdecomp_full))  ; this%f_denit_vr_col               (:,:) = ival
    allocate(this%pot_f_nit_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%pot_f_nit_vr_col             (:,:) = nan
    allocate(this%pot_f_denit_vr_col        (begc:endc,1:nlevdecomp_full))  ; this%pot_f_denit_vr_col           (:,:) = nan
    allocate(this%n2_n2o_ratio_denit_vr_col (begc:endc,1:nlevdecomp_full))  ; this%n2_n2o_ratio_denit_vr_col    (:,:) = ival
    allocate(this%f_n2o_denit_vr_col        (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_denit_vr_col           (:,:) = ival
    allocate(this%f_n2o_nit_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_nit_vr_col             (:,:) = ival

    allocate(this%sminn_to_denit_excess_vr_col(begc:endc,1:nlevdecomp_full)); this%sminn_to_denit_excess_vr_col (:,:)   = ival

    ! inorganic P transformation
    allocate(this%primp_to_labilep_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%primp_to_labilep_vr_col      (:,:) = ival
    allocate(this%labilep_to_secondp_vr_col (begc:endc,1:nlevdecomp_full))  ; this%labilep_to_secondp_vr_col    (:,:) = ival
    allocate(this%secondp_to_labilep_vr_col (begc:endc,1:nlevdecomp_full))  ; this%secondp_to_labilep_vr_col    (:,:) = ival
    allocate(this%secondp_to_occlp_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%secondp_to_occlp_vr_col      (:,:) = ival

    ! gases
    allocate(this%hr_vr_col                 (begc:endc,1:nlevdecomp_full))  ; this%hr_vr_col                    (:,:) = ival

    allocate(this%f_co2_soil_vr_col         (begc:endc,1:nlevdecomp_full))  ; this%f_co2_soil_vr_col            (:,:) = ival
    allocate(this%f_n2o_soil_vr_col         (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_soil_vr_col            (:,:) = ival
    allocate(this%f_n2_soil_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%f_n2_soil_vr_col             (:,:) = ival
    
    allocate(this%phr_vr_col                (begc:endc,1:nlevdecomp_full))  ; this%phr_vr_col                   (:,:) = nan
    allocate(this%fphr_col                  (begc:endc,1:nlevgrnd))         ; this%fphr_col                     (:,:) = nan

    ! fpi, fpg
    allocate(this%fpi_vr_col                (begc:endc,1:nlevdecomp_full))  ; this%fpi_vr_col                   (:,:) = nan
    allocate(this%fpi_col                   (begc:endc))                    ; this%fpi_col                      (:)   = nan
    allocate(this%fpg_col                   (begc:endc))                    ; this%fpg_col                      (:)   = nan
    allocate(this%fpi_p_vr_col              (begc:endc,1:nlevdecomp_full))  ; this%fpi_p_vr_col                 (:,:) = nan
    allocate(this%fpi_p_col                 (begc:endc))                    ; this%fpi_p_col                    (:)   = nan
    allocate(this%fpg_p_col                 (begc:endc))                    ; this%fpg_p_col                    (:)   = nan

    !!------------------------------------------------------------------------------------------
    !! pflotran variables: BEGIN
    !!------------------------------------------------------------------------------------------
    ! bgc rates/fluxes to decomposition pools
    allocate(this%externalc_to_decomp_cpools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%externalc_to_decomp_cpools_col(:,:,:) = spval
    allocate(this%externaln_to_decomp_npools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%externaln_to_decomp_npools_col(:,:,:) = spval
    allocate(this%externalp_to_decomp_ppools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%externalp_to_decomp_ppools_col(:,:,:) = spval

    ! bgc rates/fluxes to nh4 / no3
    allocate(this%externaln_to_nh4_col      (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_nh4_col         (:,:) = spval
    allocate(this%externaln_to_no3_col      (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_no3_col         (:,:) = spval
    allocate(this%externaln_to_sminn_col    (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_sminn_col       (:,:) = spval
    allocate(this%smin_no3_leached_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_leached_vr_col      (:,:) = ival
    allocate(this%smin_no3_runoff_vr_col    (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_runoff_vr_col       (:,:) = ival
    allocate(this%no3_net_transport_vr_col  (begc:endc,1:nlevdecomp_full))  ; this%no3_net_transport_vr_col     (:,:) = spval

    ! bgc rates/fluxes to mineral P
    allocate(this%externalp_to_primp_col    (begc:endc,1:nlevdecomp_full))  ; this%externalp_to_primp_col       (:,:) = spval
    allocate(this%externalp_to_labilep_col  (begc:endc,1:nlevdecomp_full))  ; this%externalp_to_labilep_col     (:,:) = spval
    allocate(this%externalp_to_solutionp_col(begc:endc,1:nlevdecomp_full))  ; this%externalp_to_solutionp_col   (:,:) = spval
    allocate(this%sminp_leached_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%sminp_leached_vr_col         (:,:) = ival
    allocate(this%sminp_net_transport_vr_col(begc:endc, 1:nlevdecomp_full)) ; this%sminp_net_transport_vr_col   (:,:) = spval

    ! gases:
    allocate(this%f_ngas_decomp_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_decomp_vr_col         (:,:) = ival
    allocate(this%f_ngas_nitri_vr_col       (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_nitri_vr_col          (:,:) = ival
    allocate(this%f_ngas_denit_vr_col       (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_denit_vr_col          (:,:) = ival

    ! waterflux_vars:
    allocate(this%qflx_top_soil_col     (begc:endc))                        ; this%qflx_top_soil_col    (:)   = ival
    allocate(this%qflx_ev_h2osfc_col    (begc:endc))                        ; this%qflx_ev_h2osfc_col   (:)   = ival
    allocate(this%qflx_evap_soi_col     (begc:endc))                        ; this%qflx_evap_soi_col    (:)   = ival
    allocate(this%qflx_sub_snow_col     (begc:endc))                        ; this%qflx_sub_snow_col    (:)   = ival
    allocate(this%qflx_tran_veg_col     (begc:endc))                        ; this%qflx_tran_veg_col    (:)   = ival

    ! energyflux_vars:
    allocate( this%htvp_col             (begc:endc))                        ; this%htvp_col             (:)   = ival
    allocate( this%eflx_bot_col         (begc:endc))                        ; this%eflx_bot_col         (:)   = ival
    allocate( this%eflx_gnet_col        (begc:endc))                        ; this%eflx_bot_col         (:)   = ival
    allocate( this%eflx_soil_grnd_col   (begc:endc))                        ; this%eflx_soil_grnd_col   (:)   = ival

    ! atm2lnd:
    allocate(this%forc_pbot_not_downscaled_grc  (begg:endg))                ; this%forc_pbot_not_downscaled_grc  (:)   = ival
    allocate(this%forc_pco2_grc                 (begg:endg))                ; this%forc_pco2_grc                 (:)   = ival
    allocate(this%forc_pch4_grc                 (begg:endg))                ; this%forc_pch4_grc                 (:)   = ival
    !!------------------------------------------------------------------------------------------
    !! pflotran variables: END
    !!------------------------------------------------------------------------------------------

  end subroutine InitAllocate
!!-------------------------------------------------------------------------------------------------

end module clm_bgc_interface_data
