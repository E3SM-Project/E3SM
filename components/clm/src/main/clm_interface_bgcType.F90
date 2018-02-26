module clm_interface_bgcType
!=================================================================================================
! CLM BioGeoChemistry (BGC) Interface: Data Type (Variables)
! created: 8/25/2015
! update: 9/16/2016, 2/2/2017, May-2017, June-2017
! update: 5/13/2019 (all are IMPLICITLY column-based coupling,
!        esp. after ELM v2 data-structure modification @3/12/2019, commit 1bf22e32d)
!=================================================================================================
  ! USES:
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)

  implicit none
!
  private

  type, public :: clm_interface_bgc_datatype
     ! coupling level (grid/column/patch):
     integer :: cpl_level = 2                                ! coupling level: 1 - grid, 2 - column (default), 3 - patch

     ! clm_varpar
     integer                    :: nlevdecomp_full                          ! num of CLM soil layers that are mapped to/from PFLOTRAN
     integer                    :: ndecomp_pools                            ! num of decomposition pools

     ! decomp_cascade_con
     logical, pointer           :: floating_cn_ratio                (:)     ! TRUE => pool has fixed C:N ratio
     logical, pointer           :: floating_cp_ratio                (:)     ! TRUE => pool has fixed C:P ratio
     character(len=8), pointer  :: decomp_pool_name                 (:)     ! name of pool
     real(r8), pointer          :: initial_cn_ratio                 (:)     ! c:n ratio for initialization of pools
     real(r8), pointer          :: initial_cp_ratio                 (:)     ! c:p ratio for initialization of pools

     ! ch4
     real(r8), pointer  :: finundated                           (:)     ! col fractional inundated area (excluding dedicated wetland cols)
     real(r8), pointer  :: o2stress_unsat                       (:,:)   ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer  :: o2stress_sat                         (:,:)   ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer  :: conc_o2_sat                          (:,:)   ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer  :: conc_o2_unsat                        (:,:)   ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer  :: o2_decomp_depth_sat                  (:,:)   ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer  :: o2_decomp_depth_unsat                (:,:)   ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)

     ! cnstate_vars:
     real(r8) , pointer :: rf_decomp_cascade                    (:,:,:) ! col respired fraction in decomposition step (frac)
     real(r8) , pointer :: pathfrac_decomp_cascade              (:,:,:) ! col what fraction of C leaving a given pool passes through a given

     ! carbonstate_vars:
     real(r8), pointer :: decomp_cpools_vr                      (:,:,:) ! col (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools

     ! nitrogenstate_vars:
     real(r8), pointer :: decomp_npools_vr                      (:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: sminn_vr                              (:,:)   ! col (gN/m3) vertically-resolved soil mineral N
     real(r8), pointer :: smin_no3_vr                           (:,:)   ! col (gN/m3) vertically-resolved soil mineral NO3
     real(r8), pointer :: smin_nh4_vr                           (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4
     real(r8), pointer :: smin_nh4sorb_vr                       (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4 absorbed

     ! phosphorusstate_vars:
     real(r8), pointer :: decomp_ppools_vr                      (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
     real(r8), pointer :: solutionp_vr                          (:,:)   ! col (gP/m3) vertically-resolved soil solution P
     real(r8), pointer :: labilep_vr                            (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
     real(r8), pointer :: secondp_vr                            (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
     real(r8), pointer :: occlp_vr                              (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
     real(r8), pointer :: primp_vr                              (:,:)   ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: sminp_vr                              (:,:)   ! col (gP/m3) vertically-resolved soil parimary mineral P


     ! plant NP demand
     real(r8), pointer :: plant_ndemand                         (:)     ! col N flux required to support initial GPP (gN/m3/s)
     real(r8), pointer :: plant_pdemand                         (:)     ! col P flux required to support initial GPP (gP/m3/s)
     real(r8), pointer :: plant_ndemand_vr                      (:,:)   ! col vertically-resolved N flux required to support initial GPP (gN/m3/s)
     real(r8), pointer :: plant_pdemand_vr                      (:,:)   ! col vertically-resolved P flux required to support initial GPP (gP/m3/s)

     ! decomposition flux:
     real(r8), pointer :: decomp_cpools_sourcesink              (:,:,:) ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)
     real(r8), pointer :: decomp_npools_sourcesink              (:,:,:) ! col (gN/m3) change in decomposing n pools
     real(r8), pointer :: decomp_ppools_sourcesink              (:,:,:) ! col (gP/m3) change in decomposing n pools

     real(r8), pointer :: decomp_cascade_ctransfer_vr           (:,:,:) ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
     real(r8), pointer :: decomp_cascade_ntransfer_vr           (:,:,:) ! col vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_ptransfer_vr           (:,:,:) ! col vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_vr          (:,:,:) ! col vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_sminp_flux_vr          (:,:,:) ! col vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
     real(r8), pointer :: sminn_to_denit_decomp_cascade_vr      (:,:,:) ! col vertically-resolved denitrification along decomp cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_hr_vr                  (:,:,:) ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)

     real(r8), pointer :: o_scalar                              (:,:)   ! fraction by which decomposition is limited by anoxia
     real(r8), pointer :: w_scalar                              (:,:)   ! fraction by which decomposition is limited by moisture availability
     real(r8), pointer :: t_scalar                              (:,:)   ! fraction by which decomposition is limited by temperature

     ! mineralization / immobilization / uptake flux
     real(r8), pointer :: gross_nmin_vr                         (:,:)   ! col vertically-resolved gross rate of N mineralization (gN/m3/s)
     real(r8), pointer :: net_nmin_vr                           (:,:)   ! col vertically-resolved net rate of N mineralization (gN/m3/s)

     real(r8), pointer :: potential_immob_vr                    (:,:)   ! col vertically-resolved potential N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_vr                       (:,:)   ! col vertically-resolved actual N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_no3_vr                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
     real(r8), pointer :: actual_immob_nh4_vr                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)
     real(r8), pointer :: supplement_to_sminn_vr                (:,:)   ! col vertically-resolved supplemental N supply (gN/m3/s)

     real(r8), pointer :: sminn_to_plant_vr                     (:,:)   ! col vertically-resolved plant uptake of soil mineral N (gN/m3/s)
     real(r8), pointer :: smin_no3_to_plant_vr                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
     real(r8), pointer :: smin_nh4_to_plant_vr                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)

     real(r8), pointer :: potential_immob                       (:)     ! col vert-int (diagnostic) potential N immobilization (gN/m2/s)
     real(r8), pointer :: actual_immob                          (:)     ! col vert-int (diagnostic) actual N immobilization (gN/m2/s)
     real(r8), pointer :: sminn_to_plant                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)

     real(r8), pointer :: potential_immob_p_vr                  (:,:)   ! col vertically-resolved potential P immobilization (gP/m3/s) at each level
     real(r8), pointer :: actual_immob_p_vr                     (:,:)   ! col vertically-resolved actual P immobilization (gP/m3/s) at each level
     real(r8), pointer :: sminp_to_plant_vr                     (:,:)   ! col vertically-resolved plant uptake of soil mineral P (gP/m3/s)
     real(r8), pointer :: supplement_to_sminp_vr                (:,:)   ! col vertically-resolved supplemental P supply (gP/m3/s)
     
     real(r8), pointer :: gross_pmin_vr                         (:,:)   ! col vertically-resolved gross rate of P mineralization (gP/m3/s)
     real(r8), pointer :: net_pmin_vr                           (:,:)   ! col vertically-resolved net rate of P mineralization (gP/m3/s)

     real(r8), pointer :: potential_immob_p                     (:)     ! col vert-int (diagnostic) potential P immobilization (gP/m2/s)
     real(r8), pointer :: actual_immob_p                        (:)     ! col vert-int (diagnostic) actual P immobilization (gP/m2/s)
     real(r8), pointer :: sminp_to_plant                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral P (gP/m2/s)
     
     ! nitrification / denitrification flux:
     real(r8), pointer :: f_nit_vr                              (:,:)   ! col (gN/m3/s) soil nitrification flux
     real(r8), pointer :: f_denit_vr                            (:,:)   ! col (gN/m3/s) soil denitrification flux
     real(r8), pointer :: pot_f_nit_vr                          (:,:)   ! col (gN/m3/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_vr                        (:,:)   ! col (gN/m3/s) potential soil denitrification flux
     real(r8), pointer :: n2_n2o_ratio_denit_vr                 (:,:)   ! col ratio of N2 to N2O production by denitrification [gN/gN]
     real(r8), pointer :: f_n2o_denit_vr                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_nit_vr                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]

     real(r8), pointer :: sminn_to_denit_excess_vr              (:,:)   ! col vertically-resolved denitrification from excess mineral N pool (gN/m3/s)

     ! inorganic P transformation
     real(r8), pointer :: primp_to_labilep_vr                   (:,:)   ! col (gP/m3/s) flux of P from primary mineral to labile
     real(r8), pointer :: labilep_to_secondp_vr                 (:,:)   ! col (gP/m3/s) flux of labile P to secondary mineral P
     real(r8), pointer :: secondp_to_labilep_vr                 (:,:)   ! col (gP/m3/s) flux of the desorption of secondary mineral P to labile P
     real(r8), pointer :: secondp_to_occlp_vr                   (:,:)   ! col (gP/m3/s) flux of the occlusion of secondary P to occluded P

     ! gases
     real(r8), pointer :: hr_vr                                 (:,:)   ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
     
     real(r8), pointer :: f_co2_soil_vr                         (:,:)   ! total vertically-resolved soil-atm. CO2 exchange (gC/m3/s)
     real(r8), pointer :: f_n2o_soil_vr                         (:,:)   ! col flux of N2o from soil-N processes [gN/m^3/s]
     real(r8), pointer :: f_n2_soil_vr                          (:,:)   ! col flux of N2 from soil-N processes [gN/m^3/s]

     real(r8), pointer :: phr_vr                                (:,:)   ! potential hr (not N-limited) (gC/m3/s)
     real(r8), pointer :: fphr                                  (:,:)   ! fraction of potential heterotrophic respiration

     ! fpi, fpg
     real(r8) , pointer :: fpi_vr                               (:,:)   ! col fraction of potential immobilization (no units)
     real(r8) , pointer :: fpi                                  (:)     ! col fraction of potential immobilization (no units)
     real(r8),  pointer :: fpg                                  (:)     ! col fraction of potential gpp (no units)

     real(r8) , pointer :: fpi_p_vr                             (:,:)   ! col fraction of potential immobilization (no units)
     real(r8) , pointer :: fpi_p                                (:)     ! col fraction of potential immobilization (no units)
     real(r8),  pointer :: fpg_p                                (:)     ! col fraction of potential gpp (no units)

     !------------------------------------------------------------------------------------------
     ! pflotran variables: BEGIN
     !------------------------------------------------------------------------------------------
     ! bgc rates/fluxes (previous time-step) to decomposition pools
     real(r8), pointer :: externalc_to_decomp_cpools            (:,:,:) ! col (gC/m3/s) net C fluxes associated with litter/som-adding/removal to decomp pools
     real(r8), pointer :: externaln_to_decomp_npools            (:,:,:) ! col (gN/m3/s) net N fluxes associated with litter/som-adding/removal to decomp pools
     real(r8), pointer :: externalp_to_decomp_ppools            (:,:,:) ! col (gP/m3/s) net P fluxes associated with litter/som-adding/removal to decomp pools
     real(r8), pointer :: decomp_k_pools                            (:)     ! rate constant for each decomposition pool (1./sec)
     real(r8), pointer :: sitefactor_kd_vr                      (:,:)   ! a site factor for adjusting rate constant of all decomposition pools (-) (c,j)
     real(r8), pointer :: adfactor_kd_pools                         (:)     ! a speed-up factor for adjusting rate constant of individual decomposition pool (-) (k)
     
     ! bgc rates/fluxes (previous time-step) to nh4 / no3
     real(r8), pointer :: externaln_to_nh4                      (:,:)   ! col (gN/m3/s) net N fluxes to nh4 pool: deposition + fertilization + supplement + nfix + soyfixn
     real(r8), pointer :: externaln_to_no3                      (:,:)   ! col (gN/m3/s) net N fluxes to no3 pool: deposition + fertilization + supplement
     real(r8), pointer :: externaln_to_sminn                    (:,:)   ! col (gN/m3/s) net N fluxes to sminn pool: deposition + fertilization + supplement + nfix + soyfixn
     real(r8), pointer :: smin_no3_leached_vr                   (:,:)   ! col vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
     real(r8), pointer :: smin_no3_runoff_vr                    (:,:)   ! col vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)

     ! bgc rates/fluxes (previous time-step) to mineral P
     real(r8), pointer :: externalp_to_primp                    (:,:)   ! pdep_to_sminp                         (:)     ! col atmospheric P deposition to soil mineral P (gP/m2/s)
     real(r8), pointer :: externalp_to_labilep                  (:,:)   ! fert_p_to_sminp                       (:)     ! col fertilizer P to soil mineral P (gP/m2/s)
     real(r8), pointer :: externalp_to_solutionp                (:,:)   ! supplement_to_sminp_vr                (:,:)   ! col vertically-resolved supplemental P supply (gP/m3/s)
     real(r8), pointer :: sminp_leached_vr                      (:,:)   ! col vertically-resolved soil mineral P pool loss to leaching (gP/m3/s)
     real(r8), pointer :: sminp_net_transport_vr                (:,:)   ! sminp_leached_vr                      (:,:) col net sminp transport associated with runoff/leaching (gN/m3/s)

     ! gases:
     real(r8), pointer :: f_ngas_decomp_vr                      (:,:)   ! col vertically-resolved N emission from excess mineral N pool due to mineralization (gN/m3/s)
     real(r8), pointer :: f_ngas_nitri_vr                       (:,:)   ! col vertically-resolved N emission from nitrification (gN/m3/s)
     real(r8), pointer :: f_ngas_denit_vr                       (:,:)   ! col vertically-resolved N emission from denitrification (gN/m3/s)

     ! aq. phases:
     real(r8), pointer :: no3_net_transport_vr                  (:,:)   ! col net NO3 transport associated with runoff/leaching (gN/m3/s) - also store PF's N transport inc. diffusion at current time-step
     real(r8), pointer :: nh4_net_transport_vr                  (:,:)   ! col net NH4 transport associated with runoff/leaching (gN/m3/s) - also store PF's N transport inc. diffusion at current time-step

     ! atm2lnd:
     real(r8), pointer :: forc_pco2                             (:)     ! CO2 partial pressure (Pa)
     real(r8), pointer :: forc_pch4                             (:)     ! CH4 partial pressure (Pa)

     ! mass balance check:
     ! summary of layer 1:nlevdecomp_full
     real(r8), pointer :: soil_begcb                            (:)     ! soil organic carbon mass, beginning of time step (gC/m**2)
     real(r8), pointer :: soil_begnb                            (:)     ! soil nitrogen mass, beginning of time step (gN/m**2)
     real(r8), pointer :: soil_begnb_org                        (:)     ! soil organic nitrogen mass, beginning of time step (gN/m**2)
     real(r8), pointer :: soil_begnb_min                        (:)     ! soil mineral nitrogen mass, beginning of time step (gN/m**2) = no3 + nh4 + nh4sorb

     !------------------------------------------------------------------------------------------
     ! pflotran variables: END
     !------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type clm_interface_bgc_datatype
!-------------------------------------------------------------------------------------------------

  
contains


!-------------------------------------------------------------------------------------------------
  subroutine Init(this, bounds)
     use decompMod               , only : bounds_type
     class(clm_interface_bgc_datatype) :: this
     type(bounds_type), intent(in)      :: bounds

     call this%InitAllocate (bounds)
  end subroutine Init
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
    ! USES
    use clm_varpar            , only : nlevsno, nlevgrnd
    use clm_varpar            , only : nlevdecomp_full, ndecomp_pools,  ndecomp_cascade_transitions
    use clm_varcon            , only : spval
    use decompMod             , only : bounds_type

    ! ARGUMENTS:
    real(r8) :: ival  = 0.0_r8  ! initial value
    class(clm_interface_bgc_datatype)   :: this
    type(bounds_type), intent(in)       :: bounds

    ! LOCAL VARIABLES:
    integer  :: begc, endc
    !------------------------------------------------------------------------
    begc = bounds%begc; endc= bounds%endc
    if (this%cpl_level==1) then
      begc = bounds%begg; endc= bounds%endg
    elseif (this%cpl_level==3) then
      begc = bounds%begp; endc= bounds%endp
    endif

    ! decomp_cascade_type
    allocate(this%decomp_pool_name      (1:ndecomp_pools))                  ; this%decomp_pool_name         (:)   = ''
    allocate(this%floating_cn_ratio     (1:ndecomp_pools))                  ; this%floating_cn_ratio        (:)   = .false.
    allocate(this%floating_cp_ratio     (1:ndecomp_pools))                  ; this%floating_cp_ratio        (:)   = .false.
    allocate(this%initial_cn_ratio      (0:ndecomp_pools))                  ; this%initial_cn_ratio         (:)   = nan
    allocate(this%initial_cp_ratio      (0:ndecomp_pools))                  ; this%initial_cp_ratio         (:)   = nan

    ! ch4
    allocate(this%finundated                (begc:endc))                ; this%finundated               (:)   = nan
    allocate(this%o2stress_unsat            (begc:endc,1:nlevgrnd))     ; this%o2stress_unsat           (:,:) = nan
    allocate(this%o2stress_sat              (begc:endc,1:nlevgrnd))     ; this%o2stress_sat             (:,:) = nan
    allocate(this%conc_o2_sat               (begc:endc,1:nlevgrnd))     ; this%conc_o2_sat              (:,:) = nan
    allocate(this%conc_o2_unsat             (begc:endc,1:nlevgrnd))     ; this%conc_o2_unsat            (:,:) = nan
    allocate(this%o2_decomp_depth_sat       (begc:endc,1:nlevgrnd))     ; this%o2_decomp_depth_sat      (:,:) = nan
    allocate(this%o2_decomp_depth_unsat     (begc:endc,1:nlevgrnd))     ; this%o2_decomp_depth_unsat    (:,:) = nan

    ! cnstate_vars:
    allocate(this%rf_decomp_cascade(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));
    this%rf_decomp_cascade(:,:,:) = nan

    allocate(this%pathfrac_decomp_cascade(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));
    this%pathfrac_decomp_cascade(:,:,:) = nan

    ! carbonstate_vars:
    allocate(this%decomp_cpools_vr  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_cpools_vr(:,:,:)= ival
    allocate(this%decomp_npools_vr  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_npools_vr(:,:,:)= ival
    allocate(this%decomp_ppools_vr  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_ppools_vr(:,:,:)= ival

    ! nitrogenstate_vars:
    allocate(this%sminn_vr          (begc:endc,1:nlevdecomp_full))      ; this%sminn_vr             (:,:) = ival
    allocate(this%smin_no3_vr       (begc:endc,1:nlevdecomp_full))      ; this%smin_no3_vr          (:,:) = ival
    allocate(this%smin_nh4_vr       (begc:endc,1:nlevdecomp_full))      ; this%smin_nh4_vr          (:,:) = ival
    allocate(this%smin_nh4sorb_vr   (begc:endc,1:nlevdecomp_full))      ; this%smin_nh4sorb_vr      (:,:) = ival

    ! phosphorusstate_vars:
    allocate(this%solutionp_vr      (begc:endc,1:nlevdecomp_full))      ; this%solutionp_vr         (:,:) = ival
    allocate(this%labilep_vr        (begc:endc,1:nlevdecomp_full))      ; this%labilep_vr           (:,:) = ival
    allocate(this%secondp_vr        (begc:endc,1:nlevdecomp_full))      ; this%secondp_vr           (:,:) = ival
    allocate(this%occlp_vr          (begc:endc,1:nlevdecomp_full))      ; this%occlp_vr             (:,:) = ival
    allocate(this%primp_vr          (begc:endc,1:nlevdecomp_full))      ; this%primp_vr             (:,:) = ival
    allocate(this%sminp_vr          (begc:endc,1:nlevdecomp_full))      ; this%sminp_vr             (:,:) = ival

    ! plant N/P demand
    allocate(this%plant_ndemand         (begc:endc))                    ; this%plant_ndemand            (:)    = ival
    allocate(this%plant_pdemand         (begc:endc))                    ; this%plant_pdemand            (:)    = ival
    allocate(this%plant_ndemand_vr      (begc:endc,1:nlevdecomp_full))  ; this%plant_ndemand_vr         (:,:)  = ival
    allocate(this%plant_pdemand_vr      (begc:endc,1:nlevdecomp_full))  ; this%plant_pdemand_vr         (:,:)  = ival

    ! decomposition flux:
    allocate(this%decomp_cpools_sourcesink      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_cpools_sourcesink     (:,:,:) = ival
    allocate(this%decomp_npools_sourcesink      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_npools_sourcesink     (:,:,:) = ival
    allocate(this%decomp_ppools_sourcesink      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_ppools_sourcesink     (:,:,:) = ival

    allocate(this%decomp_cascade_ctransfer_vr   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    this%decomp_cascade_ctransfer_vr  (:,:,:) = ival
    allocate(this%decomp_cascade_ntransfer_vr   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    this%decomp_cascade_ntransfer_vr  (:,:,:) = ival
    allocate(this%decomp_cascade_ptransfer_vr   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    this%decomp_cascade_ptransfer_vr  (:,:,:) = ival

    allocate(this%decomp_cascade_sminn_flux_vr  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    this%decomp_cascade_sminn_flux_vr (:,:,:) = ival
    allocate(this%decomp_cascade_sminp_flux_vr  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    this%decomp_cascade_sminp_flux_vr (:,:,:) = ival
    allocate(this%sminn_to_denit_decomp_cascade_vr (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    this%sminn_to_denit_decomp_cascade_vr (:,:,:) = ival

    allocate(this%decomp_cascade_hr_vr          (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    this%decomp_cascade_hr_vr         (:,:,:) = ival

    allocate(this%t_scalar                      (begc:endc,1:nlevdecomp_full)); this%t_scalar (:,:)=spval
    allocate(this%w_scalar                      (begc:endc,1:nlevdecomp_full)); this%w_scalar (:,:)=spval
    allocate(this%o_scalar                      (begc:endc,1:nlevdecomp_full)); this%o_scalar (:,:)=spval


    ! mineralization / immobilization / uptake fluxes
    allocate(this%gross_nmin_vr         (begc:endc,1:nlevdecomp_full))  ; this%gross_nmin_vr            (:,:) = ival
    allocate(this%net_nmin_vr           (begc:endc,1:nlevdecomp_full))  ; this%net_nmin_vr              (:,:) = ival

    allocate(this%potential_immob_vr    (begc:endc,1:nlevdecomp_full))  ; this%potential_immob_vr       (:,:) = ival
    allocate(this%actual_immob_vr       (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_vr          (:,:) = ival
    allocate(this%actual_immob_no3_vr   (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_no3_vr      (:,:) = ival
    allocate(this%actual_immob_nh4_vr   (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_nh4_vr      (:,:) = ival

    allocate(this%sminn_to_plant_vr     (begc:endc,1:nlevdecomp_full))  ; this%sminn_to_plant_vr        (:,:) = ival
    allocate(this%smin_no3_to_plant_vr  (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_to_plant_vr     (:,:) = ival
    allocate(this%smin_nh4_to_plant_vr  (begc:endc,1:nlevdecomp_full))  ; this%smin_nh4_to_plant_vr     (:,:) = ival
    allocate(this%supplement_to_sminn_vr(begc:endc,1:nlevdecomp_full))  ; this%supplement_to_sminn_vr   (:,:) = ival

    allocate(this%sminn_to_plant        (begc:endc))                    ; this%sminn_to_plant           (:)   = ival
    allocate(this%potential_immob       (begc:endc))                    ; this%potential_immob          (:)   = ival
    allocate(this%actual_immob          (begc:endc))                    ; this%actual_immob             (:)   = ival

    allocate(this%gross_pmin_vr         (begc:endc,1:nlevdecomp_full))  ; this%gross_pmin_vr            (:,:) = ival
    allocate(this%net_pmin_vr           (begc:endc,1:nlevdecomp_full))  ; this%net_pmin_vr              (:,:) = ival

    allocate(this%potential_immob_p_vr  (begc:endc,1:nlevdecomp_full))  ; this%potential_immob_p_vr     (:,:) = ival
    allocate(this%actual_immob_p_vr     (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_p_vr        (:,:) = ival
    allocate(this%sminp_to_plant_vr     (begc:endc,1:nlevdecomp_full))  ; this%sminp_to_plant_vr        (:,:) = ival
    allocate(this%supplement_to_sminp_vr(begc:endc,1:nlevdecomp_full))  ; this%supplement_to_sminp_vr   (:,:) = ival

    allocate(this%sminp_to_plant        (begc:endc))                    ; this%sminp_to_plant           (:)   = ival
    allocate(this%potential_immob_p     (begc:endc))                    ; this%potential_immob_p        (:)   = ival
    allocate(this%actual_immob_p        (begc:endc))                    ; this%actual_immob_p           (:)   = ival

    ! nitrification / denitrification flux
    allocate(this%f_nit_vr              (begc:endc,1:nlevdecomp_full))  ; this%f_nit_vr                 (:,:) = ival
    allocate(this%f_denit_vr            (begc:endc,1:nlevdecomp_full))  ; this%f_denit_vr               (:,:) = ival
    allocate(this%pot_f_nit_vr          (begc:endc,1:nlevdecomp_full))  ; this%pot_f_nit_vr             (:,:) = nan
    allocate(this%pot_f_denit_vr        (begc:endc,1:nlevdecomp_full))  ; this%pot_f_denit_vr           (:,:) = nan
    allocate(this%n2_n2o_ratio_denit_vr (begc:endc,1:nlevdecomp_full))  ; this%n2_n2o_ratio_denit_vr    (:,:) = ival
    allocate(this%f_n2o_denit_vr        (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_denit_vr           (:,:) = ival
    allocate(this%f_n2o_nit_vr          (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_nit_vr             (:,:) = ival

    allocate(this%sminn_to_denit_excess_vr(begc:endc,1:nlevdecomp_full)); this%sminn_to_denit_excess_vr (:,:)   = ival

    ! inorganic P transformation
    allocate(this%primp_to_labilep_vr   (begc:endc,1:nlevdecomp_full))  ; this%primp_to_labilep_vr      (:,:) = ival
    allocate(this%labilep_to_secondp_vr (begc:endc,1:nlevdecomp_full))  ; this%labilep_to_secondp_vr    (:,:) = ival
    allocate(this%secondp_to_labilep_vr (begc:endc,1:nlevdecomp_full))  ; this%secondp_to_labilep_vr    (:,:) = ival
    allocate(this%secondp_to_occlp_vr   (begc:endc,1:nlevdecomp_full))  ; this%secondp_to_occlp_vr      (:,:) = ival

    ! gases
    allocate(this%hr_vr                 (begc:endc,1:nlevdecomp_full))  ; this%hr_vr                    (:,:) = ival

    allocate(this%f_co2_soil_vr         (begc:endc,1:nlevdecomp_full))  ; this%f_co2_soil_vr            (:,:) = ival
    allocate(this%f_n2o_soil_vr         (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_soil_vr            (:,:) = ival
    allocate(this%f_n2_soil_vr          (begc:endc,1:nlevdecomp_full))  ; this%f_n2_soil_vr             (:,:) = ival
    
    allocate(this%phr_vr                (begc:endc,1:nlevdecomp_full))  ; this%phr_vr                   (:,:) = nan
    allocate(this%fphr                  (begc:endc,1:nlevgrnd))         ; this%fphr                     (:,:) = nan

    ! fpi, fpg
    allocate(this%fpi_vr                (begc:endc,1:nlevdecomp_full))  ; this%fpi_vr                   (:,:) = nan
    allocate(this%fpi                   (begc:endc))                    ; this%fpi                      (:)   = nan
    allocate(this%fpg                   (begc:endc))                    ; this%fpg                      (:)   = nan
    allocate(this%fpi_p_vr              (begc:endc,1:nlevdecomp_full))  ; this%fpi_p_vr                 (:,:) = nan
    allocate(this%fpi_p                 (begc:endc))                    ; this%fpi_p                    (:)   = nan
    allocate(this%fpg_p                 (begc:endc))                    ; this%fpg_p                    (:)   = nan

    !------------------------------------------------------------------------------------------
    ! pflotran variables: BEGIN
    !------------------------------------------------------------------------------------------
    ! bgc rates/fluxes to decomposition pools
    allocate(this%externalc_to_decomp_cpools(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%externalc_to_decomp_cpools(:,:,:) = spval
    allocate(this%externaln_to_decomp_npools(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%externaln_to_decomp_npools(:,:,:) = spval
    allocate(this%externalp_to_decomp_ppools(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%externalp_to_decomp_ppools(:,:,:) = spval
    allocate(this%decomp_k_pools                (1:ndecomp_pools))                            ; 
    this%decomp_k_pools                (:)     = spval
    allocate(this%adfactor_kd_pools             (1:ndecomp_pools))                            ; 
    this%adfactor_kd_pools             (:)     = spval

    allocate(this%sitefactor_kd_vr      (begc:endc,1:nlevdecomp_full))  ; this%sitefactor_kd_vr         (:,:) = spval

    ! bgc rates/fluxes to nh4 / no3
    allocate(this%externaln_to_nh4      (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_nh4         (:,:) = spval
    allocate(this%externaln_to_no3      (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_no3         (:,:) = spval
    allocate(this%externaln_to_sminn    (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_sminn       (:,:) = spval
    allocate(this%smin_no3_leached_vr   (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_leached_vr      (:,:) = ival
    allocate(this%smin_no3_runoff_vr    (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_runoff_vr       (:,:) = ival
    allocate(this%no3_net_transport_vr  (begc:endc,1:nlevdecomp_full))  ; this%no3_net_transport_vr     (:,:) = spval
    allocate(this%nh4_net_transport_vr  (begc:endc,1:nlevdecomp_full))  ; this%nh4_net_transport_vr     (:,:) = spval

    ! bgc rates/fluxes to mineral P
    allocate(this%externalp_to_primp    (begc:endc,1:nlevdecomp_full))  ; this%externalp_to_primp       (:,:) = spval
    allocate(this%externalp_to_labilep  (begc:endc,1:nlevdecomp_full))  ; this%externalp_to_labilep     (:,:) = spval
    allocate(this%externalp_to_solutionp(begc:endc,1:nlevdecomp_full))  ; this%externalp_to_solutionp   (:,:) = spval
    allocate(this%sminp_leached_vr      (begc:endc,1:nlevdecomp_full))  ; this%sminp_leached_vr         (:,:) = ival
    allocate(this%sminp_net_transport_vr(begc:endc, 1:nlevdecomp_full)) ; this%sminp_net_transport_vr   (:,:) = spval

    ! gases:
    allocate(this%f_ngas_decomp_vr      (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_decomp_vr         (:,:) = ival
    allocate(this%f_ngas_nitri_vr       (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_nitri_vr          (:,:) = ival
    allocate(this%f_ngas_denit_vr       (begc:endc,1:nlevdecomp_full))  ; this%f_ngas_denit_vr          (:,:) = ival

    ! atm2lnd:
    allocate(this%forc_pco2                 (begc:endc))                ; this%forc_pco2                 (:)   = ival
    allocate(this%forc_pch4                 (begc:endc))                ; this%forc_pch4                 (:)   = ival

    ! mass balance check
    allocate(this%soil_begcb                (begc:endc))                ; this%soil_begcb                (:)   = ival
    allocate(this%soil_begnb                (begc:endc))                ; this%soil_begnb                (:)   = ival
    allocate(this%soil_begnb_org            (begc:endc))                ; this%soil_begnb_org            (:)   = ival
    allocate(this%soil_begnb_min            (begc:endc))                ; this%soil_begnb_min            (:)   = ival

    !------------------------------------------------------------------------------------------
    ! pflotran variables: END
    !------------------------------------------------------------------------------------------

  end subroutine InitAllocate
!-------------------------------------------------------------------------------------------------

end module clm_interface_bgcType
