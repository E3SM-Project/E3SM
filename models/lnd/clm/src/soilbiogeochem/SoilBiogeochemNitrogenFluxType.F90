module SoilBiogeochemNitrogenFluxType

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp
  use decompMod                          , only : bounds_type
  use clm_varctl                         , only : use_nitrif_denitrif, use_vertsoilc
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use abortutils                         , only : endrun
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: SoilBiogeochem_nitrogenflux_type

     ! deposition fluxes
     real(r8), pointer :: ndep_to_sminn_col                         (:)     ! col atmospheric N deposition to soil mineral N (gN/m2/s)
     real(r8), pointer :: nfix_to_sminn_col                         (:)     ! col symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s) 
     real(r8), pointer :: fert_to_sminn_col                         (:)     ! col fertilizer N to soil mineral N (gN/m2/s)
     real(r8), pointer :: soyfixn_to_sminn_col                      (:)     ! col soybean fixation to soil mineral N (gN/m2/s)

     ! decomposition fluxes
     real(r8), pointer :: decomp_cascade_ntransfer_vr_col           (:,:,:) ! col vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_ntransfer_col              (:,:)   ! col vert-int (diagnostic) transfer of N from donor to receiver pool along decomp. cascade (gN/m2/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_vr_col          (:,:,:) ! col vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_col             (:,:)   ! col vert-int (diagnostic) mineral N flux for transition along decomposition cascade (gN/m2/s)

     ! Used to update concentrations concurrently with vertical transport
     ! vertically-resolved immobilization fluxes
     real(r8), pointer :: potential_immob_vr_col                    (:,:)   ! col vertically-resolved potential N immobilization (gN/m3/s) at each level
     real(r8), pointer :: potential_immob_col                       (:)     ! col vert-int (diagnostic) potential N immobilization (gN/m2/s)
     real(r8), pointer :: actual_immob_vr_col                       (:,:)   ! col vertically-resolved actual N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_col                          (:)     ! col vert-int (diagnostic) actual N immobilization (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral N (gN/m3/s)
     real(r8), pointer :: sminn_to_plant_col                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
     real(r8), pointer :: supplement_to_sminn_vr_col                (:,:)   ! col vertically-resolved supplemental N supply (gN/m3/s)
     real(r8), pointer :: supplement_to_sminn_col                   (:)     ! col vert-int (diagnostic) supplemental N supply (gN/m2/s)
     real(r8), pointer :: gross_nmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of N mineralization (gN/m3/s)
     real(r8), pointer :: gross_nmin_col                            (:)     ! col vert-int (diagnostic) gross rate of N mineralization (gN/m2/s)
     real(r8), pointer :: net_nmin_vr_col                           (:,:)   ! col vertically-resolved net rate of N mineralization (gN/m3/s)
     real(r8), pointer :: net_nmin_col                              (:)     ! col vert-int (diagnostic) net rate of N mineralization (gN/m2/s)

     ! ---------- NITRIF_DENITRIF  ---------------------

     ! nitrification / denitrification fluxes
     real(r8), pointer :: f_nit_vr_col                              (:,:)   ! col (gN/m3/s) soil nitrification flux
     real(r8), pointer :: f_denit_vr_col                            (:,:)   ! col (gN/m3/s) soil denitrification flux
     real(r8), pointer :: f_nit_col                                 (:)     ! col (gN/m2/s) soil nitrification flux
     real(r8), pointer :: f_denit_col                               (:)     ! col (gN/m2/s) soil denitrification flux

     real(r8), pointer :: pot_f_nit_vr_col                          (:,:)   ! col (gN/m3/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_vr_col                        (:,:)   ! col (gN/m3/s) potential soil denitrification flux
     real(r8), pointer :: pot_f_nit_col                             (:)     ! col (gN/m2/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_col                           (:)     ! col (gN/m2/s) potential soil denitrification flux
     real(r8), pointer :: n2_n2o_ratio_denit_vr_col                 (:,:)   ! col ratio of N2 to N2O production by denitrification [gN/gN]
     real(r8), pointer :: f_n2o_denit_vr_col                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_denit_col                           (:)     ! col flux of N2o from denitrification [gN/m^2/s]
     real(r8), pointer :: f_n2o_nit_vr_col                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_nit_col                             (:)     ! col flux of N2o from nitrification [gN/m^2/s]

     ! immobilization / uptake fluxes
     real(r8), pointer :: actual_immob_no3_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
     real(r8), pointer :: actual_immob_nh4_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)
     real(r8), pointer :: smin_no3_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
     real(r8), pointer :: smin_nh4_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)
     real(r8), pointer :: actual_immob_no3_col                      (:)     ! col actual immobilization of NO3 (gN/m2/s)
     real(r8), pointer :: actual_immob_nh4_col                      (:)     ! col actual immobilization of NH4 (gN/m2/s)
     real(r8), pointer :: smin_no3_to_plant_col                     (:)     ! col plant uptake of soil NO3 (gN/m2/s)
     real(r8), pointer :: smin_nh4_to_plant_col                     (:)     ! col plant uptake of soil Nh4 (gN/m2/s)

     ! leaching fluxes
     real(r8), pointer :: smin_no3_leached_vr_col                   (:,:)   ! col vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
     real(r8), pointer :: smin_no3_leached_col                      (:)     ! col soil mineral NO3 pool loss to leaching (gN/m2/s)
     real(r8), pointer :: smin_no3_runoff_vr_col                    (:,:)   ! col vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
     real(r8), pointer :: smin_no3_runoff_col                       (:)     ! col soil mineral NO3 pool loss to runoff (gN/m2/s)

     ! nitrification /denitrification diagnostic quantities
     real(r8), pointer :: smin_no3_massdens_vr_col                  (:,:)   ! col (ugN / g soil) soil nitrate concentration
     real(r8), pointer :: soil_bulkdensity_col                      (:,:)   ! col (kg soil / m3) bulk density of soil
     real(r8), pointer :: k_nitr_t_vr_col                           (:,:)
     real(r8), pointer :: k_nitr_ph_vr_col                          (:,:)
     real(r8), pointer :: k_nitr_h2o_vr_col                         (:,:)
     real(r8), pointer :: k_nitr_vr_col                             (:,:)
     real(r8), pointer :: wfps_vr_col                               (:,:)
     real(r8), pointer :: fmax_denit_carbonsubstrate_vr_col         (:,:)
     real(r8), pointer :: fmax_denit_nitrate_vr_col                 (:,:)
     real(r8), pointer :: f_denit_base_vr_col                       (:,:)   ! col nitrification and denitrification fluxes
     real(r8), pointer :: diffus_col                                (:,:)   ! col diffusivity (m2/s)
     real(r8), pointer :: ratio_k1_col                              (:,:)
     real(r8), pointer :: ratio_no3_co2_col                         (:,:)
     real(r8), pointer :: soil_co2_prod_col                         (:,:)
     real(r8), pointer :: fr_WFPS_col                               (:,:)

     real(r8), pointer :: r_psi_col                                 (:,:)
     real(r8), pointer :: anaerobic_frac_col                        (:,:)

     !----------- no NITRIF_DENITRIF--------------

     ! denitrification fluxes
     real(r8), pointer :: sminn_to_denit_decomp_cascade_vr_col      (:,:,:) ! col vertically-resolved denitrification along decomp cascade (gN/m3/s) 
     real(r8), pointer :: sminn_to_denit_decomp_cascade_col         (:,:)   ! col vertically-integrated (diagnostic) denitrification along decomp cascade (gN/m2/s) 
     real(r8), pointer :: sminn_to_denit_excess_vr_col              (:,:)   ! col vertically-resolved denitrification from excess mineral N pool (gN/m3/s)
     real(r8), pointer :: sminn_to_denit_excess_col                 (:)     ! col vertically-integrated (diagnostic) denitrification from excess mineral N pool (gN/m2/s)

     ! leaching fluxes
     real(r8), pointer :: sminn_leached_vr_col                      (:,:)   ! col vertically-resolved soil mineral N pool loss to leaching (gN/m3/s)
     real(r8), pointer :: sminn_leached_col                         (:)     ! col soil mineral N pool loss to leaching (gN/m2/s)

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: denit_col                                 (:)     ! col total rate of denitrification (gN/m2/s)
     real(r8), pointer :: ninputs_col                               (:)     ! col column-level N inputs (gN/m2/s)
     real(r8), pointer :: noutputs_col                              (:)     ! col column-level N outputs (gN/m2/s)
     real(r8), pointer :: som_n_leached_col                         (:)     ! col total SOM N loss from vertical transport (gN/m^2/s)
     real(r8), pointer :: decomp_npools_leached_col                 (:,:)   ! col N loss from vertical transport from each decomposing N pool (gN/m^2/s)
     real(r8), pointer :: decomp_npools_transport_tendency_col      (:,:,:) ! col N tendency due to vertical transport in decomposing N pools (gN/m^3/s)

     ! all n pools involved in decomposition
     real(r8), pointer :: decomp_npools_sourcesink_col              (:,:,:) ! col (gN/m3) change in decomposing n pools 
                                                                            ! (sum of all additions and subtractions from stateupdate1).  

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: Summary
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory
     procedure , private :: InitCold

  end type SoilBiogeochem_nitrogenflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilbiogeochem_nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate (bounds)
    call this%InitHistory (bounds)
    call this%InitCold (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize nitrogen flux
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenflux_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begc,endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    allocate(this%ndep_to_sminn_col                 (begc:endc))                   ; this%ndep_to_sminn_col          (:)   = nan
    allocate(this%nfix_to_sminn_col                 (begc:endc))                   ; this%nfix_to_sminn_col          (:)   = nan
    allocate(this%fert_to_sminn_col                 (begc:endc))                   ; this%fert_to_sminn_col          (:)   = nan
    allocate(this%soyfixn_to_sminn_col              (begc:endc))                   ; this%soyfixn_to_sminn_col       (:)   = nan
    allocate(this%sminn_to_plant_col                (begc:endc))                   ; this%sminn_to_plant_col         (:)   = nan
    allocate(this%potential_immob_col               (begc:endc))                   ; this%potential_immob_col        (:)   = nan
    allocate(this%actual_immob_col                  (begc:endc))                   ; this%actual_immob_col           (:)   = nan
    allocate(this%gross_nmin_col                    (begc:endc))                   ; this%gross_nmin_col             (:)   = nan
    allocate(this%net_nmin_col                      (begc:endc))                   ; this%net_nmin_col               (:)   = nan
    allocate(this%denit_col                         (begc:endc))                   ; this%denit_col                  (:)   = nan
    allocate(this%supplement_to_sminn_col           (begc:endc))                   ; this%supplement_to_sminn_col    (:)   = nan
    allocate(this%ninputs_col                       (begc:endc))                   ; this%ninputs_col                (:)   = nan
    allocate(this%noutputs_col                      (begc:endc))                   ; this%noutputs_col               (:)   = nan
    allocate(this%som_n_leached_col                 (begc:endc))                   ; this%som_n_leached_col          (:)   = nan

    allocate(this%r_psi_col                         (begc:endc,1:nlevdecomp_full)) ; this%r_psi_col                  (:,:) = spval
    allocate(this%anaerobic_frac_col                (begc:endc,1:nlevdecomp_full)) ; this%anaerobic_frac_col         (:,:) = spval
    allocate(this%potential_immob_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%potential_immob_vr_col     (:,:) = nan
    allocate(this%actual_immob_vr_col               (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_vr_col        (:,:) = nan
    allocate(this%sminn_to_plant_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%sminn_to_plant_vr_col      (:,:) = nan
    allocate(this%supplement_to_sminn_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminn_vr_col (:,:) = nan
    allocate(this%gross_nmin_vr_col                 (begc:endc,1:nlevdecomp_full)) ; this%gross_nmin_vr_col          (:,:) = nan
    allocate(this%net_nmin_vr_col                   (begc:endc,1:nlevdecomp_full)) ; this%net_nmin_vr_col            (:,:) = nan

    allocate(this%f_nit_vr_col                      (begc:endc,1:nlevdecomp_full)) ; this%f_nit_vr_col               (:,:) = nan
    allocate(this%f_denit_vr_col                    (begc:endc,1:nlevdecomp_full)) ; this%f_denit_vr_col             (:,:) = nan
    allocate(this%smin_no3_leached_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_leached_vr_col    (:,:) = nan
    allocate(this%smin_no3_leached_col              (begc:endc))                   ; this%smin_no3_leached_col       (:)   = nan
    allocate(this%smin_no3_runoff_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_runoff_vr_col     (:,:) = nan
    allocate(this%smin_no3_runoff_col               (begc:endc))                   ; this%smin_no3_runoff_col        (:)   = nan
    allocate(this%pot_f_nit_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%pot_f_nit_vr_col           (:,:) = nan
    allocate(this%pot_f_nit_col                     (begc:endc))                   ; this%pot_f_nit_col              (:)   = nan
    allocate(this%pot_f_denit_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%pot_f_denit_vr_col         (:,:) = nan
    allocate(this%pot_f_denit_col                   (begc:endc))                   ; this%pot_f_denit_col            (:)   = nan
    allocate(this%actual_immob_no3_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_no3_vr_col    (:,:) = nan
    allocate(this%actual_immob_nh4_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_nh4_vr_col    (:,:) = nan
    allocate(this%smin_no3_to_plant_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_to_plant_vr_col   (:,:) = nan
    allocate(this%smin_nh4_to_plant_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_to_plant_vr_col   (:,:) = nan
    allocate(this%f_nit_col                         (begc:endc))                   ; this%f_nit_col                  (:)   = nan
    allocate(this%f_denit_col                       (begc:endc))                   ; this%f_denit_col                (:)   = nan
    allocate(this%n2_n2o_ratio_denit_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%n2_n2o_ratio_denit_vr_col  (:,:) = nan
    allocate(this%f_n2o_denit_col                   (begc:endc))                   ; this%f_n2o_denit_col            (:)   = nan
    allocate(this%f_n2o_denit_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_denit_vr_col         (:,:) = nan
    allocate(this%f_n2o_nit_col                     (begc:endc))                   ; this%f_n2o_nit_col              (:)   = nan
    allocate(this%f_n2o_nit_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_nit_vr_col           (:,:) = nan

    allocate(this%smin_no3_massdens_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_massdens_vr_col   (:,:) = nan
    allocate(this%soil_bulkdensity_col              (begc:endc,1:nlevdecomp_full)) ; this%soil_bulkdensity_col       (:,:) = nan
    allocate(this%k_nitr_t_vr_col                   (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_t_vr_col            (:,:) = nan
    allocate(this%k_nitr_ph_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_ph_vr_col           (:,:) = nan
    allocate(this%k_nitr_h2o_vr_col                 (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_h2o_vr_col          (:,:) = nan
    allocate(this%k_nitr_vr_col                     (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_vr_col              (:,:) = nan
    allocate(this%wfps_vr_col                       (begc:endc,1:nlevdecomp_full)) ; this%wfps_vr_col                (:,:) = nan
    allocate(this%f_denit_base_vr_col               (begc:endc,1:nlevdecomp_full)) ; this%f_denit_base_vr_col        (:,:) = nan
    allocate(this%diffus_col                        (begc:endc,1:nlevdecomp_full)) ; this%diffus_col                 (:,:) = spval
    allocate(this%ratio_k1_col                      (begc:endc,1:nlevdecomp_full)) ; this%ratio_k1_col               (:,:) = nan
    allocate(this%ratio_no3_co2_col                 (begc:endc,1:nlevdecomp_full)) ; this%ratio_no3_co2_col          (:,:) = spval
    allocate(this%soil_co2_prod_col                 (begc:endc,1:nlevdecomp_full)) ; this%soil_co2_prod_col          (:,:) = nan
    allocate(this%fr_WFPS_col                       (begc:endc,1:nlevdecomp_full)) ; this%fr_WFPS_col                (:,:) = spval

    allocate(this%fmax_denit_carbonsubstrate_vr_col (begc:endc,1:nlevdecomp_full)) ; 
    this%fmax_denit_carbonsubstrate_vr_col (:,:) = nan
    allocate(this%fmax_denit_nitrate_vr_col         (begc:endc,1:nlevdecomp_full)) ; 
    this%fmax_denit_nitrate_vr_col         (:,:) = nan

    allocate(this%decomp_cascade_ntransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%decomp_cascade_sminn_flux_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%decomp_cascade_ntransfer_col      (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%decomp_cascade_sminn_flux_col     (begc:endc,1:ndecomp_cascade_transitions                   ))

    this%decomp_cascade_ntransfer_vr_col  (:,:,:) = nan
    this%decomp_cascade_sminn_flux_vr_col (:,:,:) = nan
    this%decomp_cascade_ntransfer_col     (:,:)   = nan
    this%decomp_cascade_sminn_flux_col    (:,:)   = nan

    allocate(this%sminn_to_denit_decomp_cascade_vr_col (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%sminn_to_denit_decomp_cascade_col    (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%sminn_to_denit_excess_vr_col         (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminn_to_denit_excess_col            (begc:endc                                                 ))
    allocate(this%sminn_leached_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminn_leached_col                    (begc:endc                                                 ))
    allocate(this%decomp_npools_leached_col            (begc:endc,1:ndecomp_pools                                 ))
    allocate(this%decomp_npools_transport_tendency_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))

    this%sminn_to_denit_decomp_cascade_vr_col (:,:,:) = nan
    this%sminn_to_denit_decomp_cascade_col    (:,:)   = nan
    this%sminn_to_denit_excess_vr_col         (:,:)   = nan
    this%sminn_to_denit_excess_col            (:)     = nan
    this%sminn_leached_vr_col                 (:,:)   = nan
    this%sminn_leached_col                    (:)     = nan
    this%decomp_npools_leached_col            (:,:)   = nan
    this%decomp_npools_transport_tendency_col (:,:,:) = nan  

    allocate(this%decomp_npools_sourcesink_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_npools_sourcesink_col (:,:,:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d, hist_addfld_decomp
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer        :: k,l
    integer        :: begc, endc
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    !-------------------------------
    ! N flux variables - native to column
    !-------------------------------

    this%ndep_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='NDEP_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='atmospheric N deposition to soil mineral N', &
         ptr_col=this%ndep_to_sminn_col)

    this%nfix_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='NFIX_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='symbiotic/asymbiotic N fixation to soil mineral N', &
         ptr_col=this%nfix_to_sminn_col)

    do l = 1, ndecomp_cascade_transitions
       ! vertically integrated fluxes
       !-- mineralization/immobilization fluxes (none from CWD)
       if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
          this%decomp_cascade_sminn_flux_col(begc:endc,l) = spval
          data1dptr => this%decomp_cascade_sminn_flux_col(:,l)
          if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
             fieldname = 'SMINN_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N_'//&
                  trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))
             longname =  'mineral N flux for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                  'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
          else
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                  //'N_TO_SMINN'
             longname =  'mineral N flux for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
          endif
          call hist_addfld1d (fname=fieldname, units='gN/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       end if

       !-- transfer fluxes (none from terminal pool, if present)
       if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
          this%decomp_cascade_ntransfer_col(begc:endc,l) = spval
          data1dptr => this%decomp_cascade_ntransfer_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_'//&
               trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N'
          longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
               ' N to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
          call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       end if

       ! vertically resolved fluxes
       if ( nlevdecomp_full > 1 ) then
          !-- mineralization/immobilization fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             this%decomp_cascade_sminn_flux_vr_col(begc:endc,:,l) = spval
             data2dptr => this%decomp_cascade_sminn_flux_vr_col(:,:,l)
             if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                fieldname = 'SMINN_TO_'&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N_'//&
                     trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))//trim(vr_suffix)
                longname =  'mineral N flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                     'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             else
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                     //'N_TO_SMINN'//trim(vr_suffix)
                longname =  'mineral N flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
             endif
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
             this%decomp_cascade_ntransfer_vr_col(begc:endc,:,l) = spval
             data2dptr => this%decomp_cascade_ntransfer_vr_col(:,:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                  //'N'//trim(vr_suffix)
             longname =  'decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' N to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

       endif
    end do

    this%denit_col(begc:endc) = spval
    call hist_addfld1d (fname='DENIT', units='gN/m^2/s', &
         avgflag='A', long_name='total rate of denitrification', &
         ptr_col=this%denit_col)

    this%som_n_leached_col(begc:endc) = spval
    call hist_addfld1d (fname='SOM_N_LEACHED', units='gN/m^2/s', &
         avgflag='A', long_name='total flux of N from SOM pools due to leaching', &
         ptr_col=this%som_n_leached_col, default='inactive')

    do k = 1, ndecomp_pools
       if ( .not. decomp_cascade_con%is_cwd(k) ) then
          this%decomp_npools_leached_col(begc:endc,k) = spval
          data1dptr => this%decomp_npools_leached_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TO_LEACHING'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N leaching loss'
          call hist_addfld1d (fname=fieldname, units='gN/m^2/s', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          this%decomp_npools_transport_tendency_col(begc:endc,:,k) = spval
          data2dptr => this%decomp_npools_transport_tendency_col(:,:,k)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TNDNCY_VERT_TRANSPORT'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N tendency due to vertical transport'
          call hist_addfld_decomp (fname=fieldname, units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       end if
    end do

    if (.not. use_nitrif_denitrif) then
       do l = 1, ndecomp_cascade_transitions
          !-- denitrification fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             this%sminn_to_denit_decomp_cascade_col(begc:endc,l) = spval
             data1dptr => this%sminn_to_denit_decomp_cascade_col(:,l)
             fieldname = 'SMINN_TO_DENIT_'//trim(decomp_cascade_con%cascade_step_name(l))
             longname =  'denitrification for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr)
          endif

          if ( nlevdecomp_full > 1 ) then       
             !-- denitrification fluxes (none from CWD)
             if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
                this%sminn_to_denit_decomp_cascade_vr_col(begc:endc,:,l) = spval
                data2dptr => this%sminn_to_denit_decomp_cascade_vr_col(:,:,l)
                fieldname = 'SMINN_TO_DENIT_'//trim(decomp_cascade_con%cascade_step_name(l))//trim(vr_suffix)
                longname =  'denitrification for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                     'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
                call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                     avgflag='A', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             endif
          endif
       end do
    end if

    if (.not. use_nitrif_denitrif) then
       this%sminn_to_denit_excess_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINN_TO_DENIT_EXCESS', units='gN/m^2/s',  &
            avgflag='A', long_name='denitrification from excess mineral N pool', &
            ptr_col=this%sminn_to_denit_excess_col, default='inactive')
    end if

    if (.not. use_nitrif_denitrif) then
       this%sminn_leached_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINN_LEACHED', units='gN/m^2/s',   &
            avgflag='A', long_name='soil mineral N pool loss to leaching', &
            ptr_col=this%sminn_leached_col)
    end if

    if (.not. use_nitrif_denitrif) then
       if ( nlevdecomp_full > 1 ) then  
          this%sminn_to_denit_excess_vr_col(begc:endc,:) = spval
          call hist_addfld_decomp (fname='SMINN_TO_DENIT_EXCESS'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name='denitrification from excess mineral N pool', &
               ptr_col=this%sminn_to_denit_excess_vr_col, default='inactive')   

          this%sminn_leached_vr_col(begc:endc,:) = spval
          call hist_addfld_decomp (fname='SMINN_LEACHED'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name='soil mineral N pool loss to leaching', &
               ptr_col=this%sminn_leached_vr_col, default='inactive')
       endif
    end if

    if (use_nitrif_denitrif) then
       this%f_nit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_NIT', units='gN/m^2/s',  &
            avgflag='A', long_name='nitrification flux', &
            ptr_col=this%f_nit_col)
    end if

    if (use_nitrif_denitrif) then
       this%f_denit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_DENIT', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification flux', &
            ptr_col=this%f_denit_col)
    end if

    if (use_nitrif_denitrif) then
       this%pot_f_nit_col(begc:endc) = spval
       call hist_addfld1d (fname='POT_F_NIT', units='gN/m^2/s', &
            avgflag='A', long_name='potential nitrification flux', &
            ptr_col=this%pot_f_nit_col)
    end if

    if (use_nitrif_denitrif) then
       this%pot_f_denit_col(begc:endc) = spval
       call hist_addfld1d (fname='POT_F_DENIT', units='gN/m^2/s', &
            avgflag='A', long_name='potential denitrification flux', &
            ptr_col=this%pot_f_denit_col)
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_leached_col(begc:endc) = spval
       call hist_addfld1d (fname='SMIN_NO3_LEACHED', units='gN/m^2/s', &
            avgflag='A', long_name='soil NO3 pool loss to leaching', &
            ptr_col=this%smin_no3_leached_col)
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_runoff_col(begc:endc) = spval
       call hist_addfld1d (fname='SMIN_NO3_RUNOFF', units='gN/m^2/s', &
            avgflag='A', long_name='soil NO3 pool loss to runoff', &
            ptr_col=this%smin_no3_runoff_col)
    end if
       
    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then 
       this%f_nit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='F_NIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='nitrification flux', &
            ptr_col=this%f_nit_vr_col)
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then 
       this%f_denit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='F_DENIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='denitrification flux', &
            ptr_col=this%f_denit_vr_col)
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then 
       this%pot_f_nit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POT_F_NIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='potential nitrification flux', &
            ptr_col=this%pot_f_nit_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then 
       this%pot_f_denit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POT_F_DENIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='potential denitrification flux', &
            ptr_col=this%pot_f_denit_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then 
       this%smin_no3_leached_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_LEACHED'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='soil NO3 pool loss to leaching', &
            ptr_col=this%smin_no3_leached_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then 
       this%smin_no3_runoff_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_RUNOFF'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='soil NO3 pool loss to runoff', &
            ptr_col=this%smin_no3_runoff_vr_col, default='inactive')
    endif

    if (use_nitrif_denitrif) then
       this%n2_n2o_ratio_denit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='n2_n2o_ratio_denit', units='gN/gN', type2d='levdcmp', &
            avgflag='A', long_name='n2_n2o_ratio_denit', &
            ptr_col=this%n2_n2o_ratio_denit_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%actual_immob_no3_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB_NO3', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='immobilization of NO3', &
            ptr_col=this%actual_immob_no3_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%actual_immob_nh4_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB_NH4', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='immobilization of NH4', &
            ptr_col=this%actual_immob_nh4_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_TO_PLANT', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of NO3', &
            ptr_col=this%smin_no3_to_plant_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%smin_nh4_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NH4_TO_PLANT', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of NH4', &
            ptr_col=this%smin_nh4_to_plant_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_massdens_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_MASSDENS', units='ugN/cm^3 soil', type2d='levdcmp', &
            avgflag='A', long_name='SMIN_NO3_MASSDENS', &
            ptr_col=this%smin_no3_massdens_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_t_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR_T', units='unitless', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR_T', &
            ptr_col=this%k_nitr_t_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_ph_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR_PH', units='unitless', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR_PH', &
            ptr_col=this%k_nitr_ph_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_h2o_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR_H2O', units='unitless', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR_H2O', &
            ptr_col=this%k_nitr_h2o_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR', units='1/s', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR', &
            ptr_col=this%k_nitr_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%wfps_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='WFPS', units='percent', type2d='levdcmp', &
            avgflag='A', long_name='WFPS', &
            ptr_col=this%wfps_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%fmax_denit_carbonsubstrate_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='FMAX_DENIT_CARBONSUBSTRATE', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='FMAX_DENIT_CARBONSUBSTRATE', &
            ptr_col=this%fmax_denit_carbonsubstrate_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%fmax_denit_nitrate_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='FMAX_DENIT_NITRATE', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='FMAX_DENIT_NITRATE', &
            ptr_col=this%fmax_denit_nitrate_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%f_denit_base_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='F_DENIT_BASE', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='F_DENIT_BASE', &
            ptr_col=this%f_denit_base_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%diffus_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='diffus', units='m^2/s', type2d='levdcmp', &
            avgflag='A', long_name='diffusivity', &
            ptr_col=this%diffus_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%ratio_k1_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ratio_k1', units='none', type2d='levdcmp', &
            avgflag='A', long_name='ratio_k1', &
            ptr_col=this%ratio_k1_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%ratio_no3_co2_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ratio_no3_co2', units='ratio', type2d='levdcmp', &
            avgflag='A', long_name='ratio_no3_co2', &
            ptr_col=this%ratio_no3_co2_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%soil_co2_prod_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='soil_co2_prod', units='ug C / g soil / day', type2d='levdcmp', &
            avgflag='A', long_name='soil_co2_prod', &
            ptr_col=this%soil_co2_prod_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%fr_WFPS_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='fr_WFPS', units='fraction', type2d='levdcmp', &
            avgflag='A', long_name='fr_WFPS', &
            ptr_col=this%fr_WFPS_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%soil_bulkdensity_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='soil_bulkdensity', units='kg/m3', type2d='levdcmp', &
            avgflag='A', long_name='soil_bulkdensity', &
            ptr_col=this%soil_bulkdensity_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%anaerobic_frac_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='anaerobic_frac', units='m3/m3', type2d='levdcmp', &
            avgflag='A', long_name='anaerobic_frac', &
            ptr_col=this%anaerobic_frac_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%r_psi_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='r_psi', units='m', type2d='levdcmp', &
            avgflag='A', long_name='r_psi', &
            ptr_col=this%r_psi_col, default='inactive')
    end if


    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%potential_immob_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POTENTIAL_IMMOB'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='potential N immobilization', &
            ptr_col=this%potential_immob_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%actual_immob_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='actual N immobilization', &
            ptr_col=this%actual_immob_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%sminn_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINN_TO_PLANT'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of soil mineral N', &
            ptr_col=this%sminn_to_plant_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%supplement_to_sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SUPPLEMENT_TO_SMINN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='supplemental N supply', &
            ptr_col=this%supplement_to_sminn_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%gross_nmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='GROSS_NMIN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='gross rate of N mineralization', &
            ptr_col=this%gross_nmin_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%net_nmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='NET_NMIN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='net rate of N mineralization', &
            ptr_col=this%net_nmin_vr_col, default='inactive')
    end if

    this%potential_immob_col(begc:endc) = spval
    call hist_addfld1d (fname='POTENTIAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='potential N immobilization', &
         ptr_col=this%potential_immob_col)

    this%actual_immob_col(begc:endc) = spval
    call hist_addfld1d (fname='ACTUAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='actual N immobilization', &
         ptr_col=this%actual_immob_col)

    this%sminn_to_plant_col(begc:endc) = spval
    call hist_addfld1d (fname='SMINN_TO_PLANT', units='gN/m^2/s', &
         avgflag='A', long_name='plant uptake of soil mineral N', &
         ptr_col=this%sminn_to_plant_col)

    this%supplement_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='SUPPLEMENT_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='supplemental N supply', &
         ptr_col=this%supplement_to_sminn_col)

    this%gross_nmin_col(begc:endc) = spval
    call hist_addfld1d (fname='GROSS_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='gross rate of N mineralization', &
         ptr_col=this%gross_nmin_col)

    this%net_nmin_col(begc:endc) = spval
    call hist_addfld1d (fname='NET_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='net rate of N mineralization', &
         ptr_col=this%net_nmin_col)

    if (use_nitrif_denitrif) then
       this%f_n2o_nit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_N2O_NIT', units='gN/m^2/s', &
            avgflag='A', long_name='nitrification N2O flux', &
            ptr_col=this%f_n2o_nit_col)

       this%f_n2o_denit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_N2O_DENIT', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification N2O flux', &
            ptr_col=this%f_n2o_denit_col)
    end if

    if (crop_prog) then
       this%fert_to_sminn_col(begc:endc) = spval
       call hist_addfld1d (fname='FERT_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='fertilizer to soil mineral N', &
            ptr_col=this%fert_to_sminn_col)
    end if

    if (crop_prog) then
       this%soyfixn_to_sminn_col(begc:endc) = spval
       call hist_addfld1d (fname='SOYFIXN_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='Soybean fixation to soil mineral N', &
            ptr_col=this%soyfixn_to_sminn_col)
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use clm_varpar      , only : crop_prog
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: c,l
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    !---------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    !-----------------------------------------------
    ! initialize nitrogen flux variables
    !-----------------------------------------------

    call this%SetValues (&
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart (this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use clm_varpar, only : crop_prog
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenflux_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    if (use_nitrif_denitrif) then
       ! pot_f_nit_vr
       if (use_vertsoilc) then
          ptr2d => this%pot_f_nit_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='pot_f_nit_vr_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='potential soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%pot_f_nit_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='pot_f_nit_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: pot_f_nit_vr'//' is required on an initialization dataset' )
       end if
    end if

    if (use_nitrif_denitrif) then
       ! f_nit_vr
       if (use_vertsoilc) then
          ptr2d => this%f_nit_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='f_nit_vr_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%f_nit_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='f_nit_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: f_nit_vr'//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenflux_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          if (.not. use_nitrif_denitrif) then
             this%sminn_to_denit_excess_vr_col(i,j)      = value_column
             this%sminn_leached_vr_col(i,j)              = value_column
          else
             this%f_nit_vr_col(i,j)                      = value_column
             this%f_denit_vr_col(i,j)                    = value_column
             this%smin_no3_leached_vr_col(i,j)           = value_column
             this%smin_no3_runoff_vr_col(i,j)            = value_column 
             this%n2_n2o_ratio_denit_vr_col(i,j)         = value_column
             this%pot_f_nit_vr_col(i,j)                  = value_column
             this%pot_f_denit_vr_col(i,j)                = value_column
             this%actual_immob_no3_vr_col(i,j)           = value_column
             this%actual_immob_nh4_vr_col(i,j)           = value_column
             this%smin_no3_to_plant_vr_col(i,j)          = value_column
             this%smin_nh4_to_plant_vr_col(i,j)          = value_column
             this%f_n2o_denit_vr_col(i,j)                = value_column
             this%f_n2o_nit_vr_col(i,j)                  = value_column

             this%smin_no3_massdens_vr_col(i,j)          = value_column
             this%k_nitr_t_vr_col(i,j)                   = value_column
             this%k_nitr_ph_vr_col(i,j)                  = value_column
             this%k_nitr_h2o_vr_col(i,j)                 = value_column
             this%k_nitr_vr_col(i,j)                     = value_column 
             this%wfps_vr_col(i,j)                       = value_column 
             this%fmax_denit_carbonsubstrate_vr_col(i,j) = value_column 
             this%fmax_denit_nitrate_vr_col(i,j)         = value_column 
             this%f_denit_base_vr_col(i,j)               = value_column

             this%diffus_col(i,j)                        = value_column
             this%ratio_k1_col(i,j)                      = value_column
             this%ratio_no3_co2_col(i,j)                 = value_column
             this%soil_co2_prod_col(i,j)                 = value_column
             this%fr_WFPS_col(i,j)                       = value_column
             this%soil_bulkdensity_col(i,j)              = value_column

             this%r_psi_col(i,j)                         = value_column
             this%anaerobic_frac_col(i,j)                = value_column
          end if
          this%potential_immob_vr_col(i,j)               = value_column
          this%actual_immob_vr_col(i,j)                  = value_column
          this%sminn_to_plant_vr_col(i,j)                = value_column
          this%supplement_to_sminn_vr_col(i,j)           = value_column
          this%gross_nmin_vr_col(i,j)                    = value_column
          this%net_nmin_vr_col(i,j)                      = value_column
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)

       this%ndep_to_sminn_col(i)             = value_column
       this%nfix_to_sminn_col(i)             = value_column
       this%fert_to_sminn_col(i)             = value_column
       this%soyfixn_to_sminn_col(i)          = value_column
       this%potential_immob_col(i)           = value_column
       this%actual_immob_col(i)              = value_column
       this%sminn_to_plant_col(i)            = value_column
       this%supplement_to_sminn_col(i)       = value_column
       this%gross_nmin_col(i)                = value_column
       this%net_nmin_col(i)                  = value_column
       this%denit_col(i)                     = value_column
       if (use_nitrif_denitrif) then
          this%f_nit_col(i)                  = value_column
          this%pot_f_nit_col(i)              = value_column
          this%f_denit_col(i)                = value_column
          this%pot_f_denit_col(i)            = value_column
          this%f_n2o_denit_col(i)            = value_column
          this%f_n2o_nit_col(i)              = value_column
          this%smin_no3_leached_col(i)       = value_column
          this%smin_no3_runoff_col(i)        = value_column
       else
          this%sminn_to_denit_excess_col(i)  = value_column
          this%sminn_leached_col(i)          = value_column
       end if
       this%ninputs_col(i)                   = value_column
       this%noutputs_col(i)                  = value_column
       this%som_n_leached_col(i)             = value_column
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_npools_leached_col(i,k) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_npools_transport_tendency_col(i,j,k) = value_column
          end do
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cascade_ntransfer_col(i,l) = value_column
          this%decomp_cascade_sminn_flux_col(i,l) = value_column
          if (.not. use_nitrif_denitrif) then
             this%sminn_to_denit_decomp_cascade_col(i,l) = value_column
          end if
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_ntransfer_vr_col(i,j,l) = value_column
             this%decomp_cascade_sminn_flux_vr_col(i,j,l) = value_column
             if (.not. use_nitrif_denitrif) then
                this%sminn_to_denit_decomp_cascade_vr_col(i,j,l) = value_column
             end if
          end do
       end do
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_npools_sourcesink_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc)
    !
    ! !USES:
    use clm_varpar , only: nlevdecomp, ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl , only: use_nitrif_denitrif
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_nitrogenflux_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l   ! indices
    integer  :: fc        ! filter indices
    !-----------------------------------------------------------------------

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%denit_col(c) = 0._r8
       this%supplement_to_sminn_col(c) = 0._r8
       this%som_n_leached_col(c)       = 0._r8
    end do

    ! vertically integrate decomposing N cascade fluxes and soil mineral N fluxes associated with decomposition cascade
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             this%decomp_cascade_ntransfer_col(c,k) = &
                  this%decomp_cascade_ntransfer_col(c,k) + &
                  this%decomp_cascade_ntransfer_vr_col(c,j,k) * dzsoi_decomp(j) 

             this%decomp_cascade_sminn_flux_col(c,k) = &
                  this%decomp_cascade_sminn_flux_col(c,k) + &
                  this%decomp_cascade_sminn_flux_vr_col(c,j,k) * dzsoi_decomp(j) 
          end do
       end do
    end do

    if (.not. use_nitrif_denitrif) then

       ! vertically integrate each denitrification flux
       do l = 1, ndecomp_cascade_transitions
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%sminn_to_denit_decomp_cascade_col(c,l) = &
                     this%sminn_to_denit_decomp_cascade_col(c,l) + &
                     this%sminn_to_denit_decomp_cascade_vr_col(c,j,l) * dzsoi_decomp(j)
             end do
          end do
       end do

       ! vertically integrate bulk denitrification and  leaching flux
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%sminn_to_denit_excess_col(c) = &
                  this%sminn_to_denit_excess_col(c) + &
                  this%sminn_to_denit_excess_vr_col(c,j) * dzsoi_decomp(j)

             this%sminn_leached_col(c) = &
                  this%sminn_leached_col(c) + &
                  this%sminn_leached_vr_col(c,j) * dzsoi_decomp(j)
          end do
       end do

       ! total N denitrification (DENIT)
       do l = 1, ndecomp_cascade_transitions
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%denit_col(c) = &
                  this%denit_col(c) + &
                  this%sminn_to_denit_decomp_cascade_col(c,l)
          end do
       end do

       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%denit_col(c) =  &
               this%denit_col(c) + &
               this%sminn_to_denit_excess_col(c)
       end do

    else

       ! vertically integrate NO3 NH4 N2O fluxes and pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             ! nitrification and denitrification fluxes
             this%f_nit_col(c) = &
                  this%f_nit_col(c) + &
                  this%f_nit_vr_col(c,j) * dzsoi_decomp(j)

             this%f_denit_col(c) = &
                  this%f_denit_col(c) + &
                  this%f_denit_vr_col(c,j) * dzsoi_decomp(j)

             this%pot_f_nit_col(c) = &
                  this%pot_f_nit_col(c) + &
                  this%pot_f_nit_vr_col(c,j) * dzsoi_decomp(j)

             this%pot_f_denit_col(c) = &
                  this%pot_f_denit_col(c) + &
                  this%pot_f_denit_vr_col(c,j) * dzsoi_decomp(j)

             this%f_n2o_nit_col(c) = &
                  this%f_n2o_nit_col(c) + &
                  this%f_n2o_nit_vr_col(c,j) * dzsoi_decomp(j)

             this%f_n2o_denit_col(c) = &
                  this%f_n2o_denit_col(c) + &
                  this%f_n2o_denit_vr_col(c,j) * dzsoi_decomp(j)

             ! leaching/runoff flux
             this%smin_no3_leached_col(c) = &
                  this%smin_no3_leached_col(c) + &
                  this%smin_no3_leached_vr_col(c,j) * dzsoi_decomp(j)

             this%smin_no3_runoff_col(c) = &
                  this%smin_no3_runoff_col(c) + &
                  this%smin_no3_runoff_vr_col(c,j) * dzsoi_decomp(j)

          end do
       end do

       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%denit_col(c) = this%f_denit_col(c)
       end do

    end if

    ! supplementary N supplement_to_sminn
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%supplement_to_sminn_col(c) = &
               this%supplement_to_sminn_col(c) + &
               this%supplement_to_sminn_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_npools_leached_col(c,l) = 0._r8
       end do

       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_npools_leached_col(c,l) = &
                  this%decomp_npools_leached_col(c,l) + &
                  this%decomp_npools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do

       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%som_n_leached_col(c) = &
               this%som_n_leached_col(c) + &
               this%decomp_npools_leached_col(c,l)
       end do
    end do

  end subroutine Summary

end module soilbiogeochemNitrogenFluxType

