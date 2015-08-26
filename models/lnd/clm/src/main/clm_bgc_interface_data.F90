module clm_bgc_interface_data
!----------------------------------------------------------------------------------------------
! CLM soil bgc interface
! authors: Gangsheng WANG
!
!          Climate Change Science Institute & Environmental Science Division
!          Oak Ridge National Laboratory
!
! date: 8/25/2015

! NOTES for convenience:
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

     ! Time invariant data:

     ! num of CLM soil layers that are mapped to/from PFLOTRAN
     integer :: nlevdecomp

     ! Soil BGC decomposing pools
     integer :: ndecomp_pools
     logical, pointer :: floating_cn_ratio(:)           ! TRUE => pool has fixed C:N ratio
     character(len=8), pointer :: decomp_pool_name(:)   ! name of pool

     ! (1) hydraulic properties
     ! (1.1) col:
     real(r8), pointer :: z                    (:,:) ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: dz                   (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd)
     ! (1.2) soilstate_vars:
     real(r8), pointer :: bd_col               (:,:) ! col bulk density of dry soil material [kg/m^3] (CN)
     real(r8), pointer :: hksat_col            (:,:) ! col hydraulic conductivity at saturation (mm H2O /s)
     real(r8), pointer :: bsw_col              (:,:) ! col Clapp and Hornberger "b" (nlevgrnd)
     real(r8), pointer :: watsat_col           (:,:) ! col volumetric soil water at saturation (porosity)
     real(r8), pointer :: sucsat_col           (:,:) ! col minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: watfc_col            (:,:) ! col volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: porosity_col         (:,:) ! col soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity_col     (:,:) ! col effective porosity = porosity - vol_ice (nlevgrnd)

     ! (2) soil thermohydrology
     ! (2.1) soilstate_vars:
     real(r8), pointer :: soilpsi_col          (:,:) ! col soil water potential in each soil layer (MPa) (CN)
     ! (2.2) waterstate_vars:
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     ! (2.3) temperature_vars:
     real(r8), pointer :: t_soisno_col             (:,:) ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

     ! (3) bgc state variables:
     ! (3.1) carbonstate_vars:
     real(r8), pointer :: decomp_cpools_vr_col    (:,:,:)  ! col (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     ! (3.2) nitrogenstate_vars:
     real(r8), pointer :: decomp_npools_vr_col         (:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: sminn_vr_col                 (:,:)   ! col (gN/m3) vertically-resolved soil mineral N
     real(r8), pointer :: smin_no3_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NO3
     real(r8), pointer :: smin_nh4_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4
     real(r8), pointer :: smin_nh4sorb_vr_col          (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4 absorbed
     ! (3.3) phosphorusstate_vars:
     real(r8), pointer :: decomp_ppools_vr_col         (:,:,:)     ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
     real(r8), pointer :: solutionp_vr_col             (:,:)       ! col (gP/m3) vertically-resolved soil solution P
     real(r8), pointer :: labilep_vr_col               (:,:)       ! col (gP/m3) vertically-resolved soil labile mineral P
     real(r8), pointer :: secondp_vr_col               (:,:)       ! col (gP/m3) vertically-resolved soil secondary mineralP
     real(r8), pointer :: occlp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil occluded mineral P
     real(r8), pointer :: primp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: sminp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil parimary mineral P

     ! (4) bgc rates (fluxes), pass from clm to interface
     ! (4.1) to decomposition pools
     real(r8), pointer :: externalc_to_decomp_cpools_col            (:,:,:) ! col (gC/m3/s) net C fluxes associated with litter/som-adding/removal to decomp pools
     real(r8), pointer :: externaln_to_decomp_npools_col            (:,:,:) ! col (gN/m3/s) net N fluxes associated with litter/som-adding/removal to decomp pools
     real(r8), pointer :: externalp_to_decomp_ppools_col            (:,:,:) ! col (gP/m3/s) net P fluxes associated with litter/som-adding/removal to decomp pools
     ! (4.2) to nh4 / no3
     real(r8), pointer :: externaln_to_nh4_col            (:,:) ! col (gN/m3/s) net N fluxes to nh4 pool: deposition + fertilization + supplement + nfix + soyfixn
     real(r8), pointer :: externaln_to_no3_col            (:,:) ! col (gN/m3/s) net N fluxes to no3 pool: deposition + fertilization + supplement
     real(r8), pointer :: smin_no3_leached_vr_col         (:,:)   ! col vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
     real(r8), pointer :: smin_no3_runoff_vr_col          (:,:)   ! col vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
     real(r8), pointer :: no3_net_transport_vr_col        (:,:)   ! col net NO3 transport associated with runoff/leaching (gN/m3/s)
     ! (4.3) to mineral P
     real(r8), pointer :: externalp_to_primp_col          (:,:) ! pdep_to_sminp_col                         (:)     ! col atmospheric P deposition to soil mineral P (gP/m2/s)
     real(r8), pointer :: externalp_to_labilep_col        (:,:) ! fert_p_to_sminp_col                       (:)     ! col fertilizer P to soil mineral P (gP/m2/s)
     real(r8), pointer :: externalp_to_solutionp_col      (:,:) ! supplement_to_sminp_vr_col                (:,:)   ! col vertically-resolved supplemental P supply (gP/m3/s)
     real(r8), pointer :: sminp_leached_vr_col            (:,:)   ! col vertically-resolved soil mineral P pool loss to leaching (gP/m3/s)
     real(r8), pointer :: sminp_net_transport_vr_col      (:,:) ! sminp_leached_vr_col                      (:,:) col net sminp transport associated with runoff/leaching (gN/m3/s)
     ! (4.4) plant demand
     real(r8), pointer :: plant_ndemand_vr_col                      (:,:)   ! col vertically-resolved N flux required to support initial GPP (gN/m3/s)
     real(r8), pointer :: plant_pdemand_vr_col                      (:,:)   ! col vertically-resolved P flux required to support initial GPP (gP/m3/s)

     ! (5) bgc fluxes, return from interface to clm
     ! (5.1) decomposition
     real(r8), pointer :: decomp_cpools_sourcesink_col              (:,:,:) ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)
     real(r8), pointer :: decomp_npools_sourcesink_col              (:,:,:) ! col (gN/m3) change in decomposing n pools
     real(r8), pointer :: decomp_ppools_sourcesink_col              (:,:,:) ! col (gP/m3) change in decomposing n pools

     real(r8), pointer :: decomp_cascade_ctransfer_vr_col           (:,:,:) ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
     real(r8), pointer :: decomp_cascade_ntransfer_vr_col           (:,:,:) ! col vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_ptransfer_vr_col           (:,:,:) ! col vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)

     ! (5.2) mineralization / immobilization / uptake fluxes
     real(r8), pointer :: gross_nmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of N mineralization (gN/m3/s)
     real(r8), pointer :: net_nmin_vr_col                           (:,:)   ! col vertically-resolved net rate of N mineralization (gN/m3/s)

     real(r8), pointer :: potential_immob_vr_col                    (:,:)   ! col vertically-resolved potential N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_vr_col                       (:,:)   ! col vertically-resolved actual N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_no3_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
     real(r8), pointer :: actual_immob_nh4_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)

     real(r8), pointer :: sminn_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral N (gN/m3/s)
     real(r8), pointer :: smin_no3_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
     real(r8), pointer :: smin_nh4_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)

     real(r8), pointer :: gross_pmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of P mineralization (gP/m3/s)
     real(r8), pointer :: net_pmin_vr_col                           (:,:)   ! col vertically-resolved net rate of P mineralization (gP/m3/s)

     real(r8), pointer :: potential_immob_p_vr_col                  (:,:)   ! col vertically-resolved potential P immobilization (gP/m3/s) at each level
     real(r8), pointer :: actual_immob_p_vr_col                     (:,:)   ! col vertically-resolved actual P immobilization (gP/m3/s) at each level
     real(r8), pointer :: sminp_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral P (gP/m3/s)

     ! (5.3) nitrification / denitrification
     real(r8), pointer :: f_nit_vr_col                              (:,:)   ! col (gN/m3/s) soil nitrification flux
     real(r8), pointer :: f_denit_vr_col                            (:,:)   ! col (gN/m3/s) soil denitrification flux

     ! (5.4) inorganic P transformation
     real(r8), pointer :: primp_to_labilep_vr_col                     (:,:)   ! col (gP/m3/s) flux of P from primary mineral to labile
     real(r8), pointer :: labilep_to_secondp_vr_col                   (:,:)   ! col (gP/m3/s) flux of labile P to secondary mineral P
     real(r8), pointer :: secondp_to_labilep_vr_col                   (:,:)   ! col (gP/m3/s) flux of the desorption of secondary mineral P to labile P
     real(r8), pointer :: secondp_to_occlp_vr_col                     (:,:)   ! col (gP/m3/s) flux of the occlusion of secondary P to occluded P

     ! (5.5) gases
     real(r8), pointer :: decomp_cascade_hr_vr_col                  (:,:,:) ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
     real(r8), pointer :: hr_vr_col                                 (:,:)   ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)

     real(r8), pointer :: f_n2o_denit_vr_col                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_nit_vr_col                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]

  !---------------------------------------------------------------
  contains
     procedure , public  :: Init
!     procedure , public  :: SetValues
     procedure , private :: InitAllocate
  end type clm_bgc_interface_data_type
  !---------------------------------------------------------------

!  type(clm_bgc_interface_data_type) , public, target , save :: clm_bgc_data
  
contains

! ************************************************************************** !
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
     class(clm_bgc_interface_data_type) :: this
     type(bounds_type), intent(in) :: bounds

     call this%InitAllocate ( bounds)
  end subroutine Init
  !------------------------------------------------------------------------

!  subroutine SetValues(this, bounds)
!     class(clm_bgc_interface_data_type) :: this
!     type(bounds_type), intent(in) :: bounds
!  end subroutine SetValues
  !------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
    !! USES


    !! ARGUMENTS:
    class(clm_bgc_interface_data_type) :: this
    type(bounds_type), intent(in) :: bounds

    !! LOCAL VARIABLES:
    integer           :: begc,endc
    !------------------------------------------------------------------------
    begc = bounds%begc; endc = bounds%endc

    allocate(this%z                     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%z           (:,:) = nan
    allocate(this%dz                    (begc:endc,-nlevsno+1:nlevgrnd)) ; this%dz          (:,:) = nan

    allocate(this%bd_col                (begc:endc,nlevgrnd))            ; this%bd_col               (:,:) = nan
    allocate(this%hksat_col             (begc:endc,nlevgrnd))            ; this%hksat_col            (:,:) = spval
    allocate(this%bsw_col               (begc:endc,nlevgrnd))            ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col            (begc:endc,nlevgrnd))            ; this%watsat_col           (:,:) = nan
    allocate(this%sucsat_col            (begc:endc,nlevgrnd))            ; this%sucsat_col           (:,:) = spval
    allocate(this%watfc_col             (begc:endc,nlevgrnd))            ; this%watfc_col            (:,:) = nan
    allocate(this%porosity_col          (begc:endc,nlevgrnd))              ; this%porosity_col         (:,:) = spval
    allocate(this%eff_porosity_col      (begc:endc,nlevgrnd))            ; this%eff_porosity_col     (:,:) = spval

    allocate(this%soilpsi_col           (begc:endc,nlevgrnd))            ; this%soilpsi_col          (:,:) = nan

    allocate(this%h2osoi_liq_col        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_col         (:,:) = nan
    allocate(this%h2osoi_ice_col        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_col         (:,:) = nan

    allocate(this%t_soisno_col          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_soisno_col             (:,:) = nan

    allocate(this%decomp_cpools_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_cpools_vr_col(:,:,:)= nan
    allocate(this%decomp_npools_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_npools_vr_col(:,:,:)= nan
    allocate(this%decomp_ppools_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_pools));  this%decomp_ppools_vr_col(:,:,:)= nan

    allocate(this%sminn_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%sminn_vr_col             (:,:) = nan
    allocate(this%smin_no3_vr_col       (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_vr_col          (:,:) = nan
    allocate(this%smin_nh4_vr_col       (begc:endc,1:nlevdecomp_full))  ; this%smin_nh4_vr_col          (:,:) = nan
    allocate(this%smin_nh4sorb_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%smin_nh4sorb_vr_col      (:,:) = nan

    allocate(this%solutionp_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%solutionp_vr_col         (:,:) = nan
    allocate(this%labilep_vr_col        (begc:endc,1:nlevdecomp_full))  ; this%labilep_vr_col           (:,:) = nan
    allocate(this%secondp_vr_col        (begc:endc,1:nlevdecomp_full))  ; this%secondp_vr_col           (:,:) = nan
    allocate(this%occlp_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%occlp_vr_col             (:,:) = nan
    allocate(this%primp_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%primp_vr_col             (:,:) = nan
    allocate(this%sminp_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%sminp_vr_col             (:,:) = nan

    allocate(this%externalc_to_decomp_cpools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%externalc_to_decomp_cpools_col(:,:,:) = spval
    allocate(this%externaln_to_decomp_npools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%externaln_to_decomp_npools_col(:,:,:) = spval
    allocate(this%externalp_to_decomp_ppools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%externalp_to_decomp_ppools_col(:,:,:) = spval

    allocate(this%externaln_to_nh4_col      (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_nh4_col         (:,:) = spval
    allocate(this%externaln_to_no3_col      (begc:endc,1:nlevdecomp_full))  ; this%externaln_to_no3_col         (:,:) = spval
    allocate(this%smin_no3_leached_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_leached_vr_col      (:,:) = nan
    allocate(this%smin_no3_runoff_vr_col    (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_runoff_vr_col       (:,:) = nan
    allocate(this%no3_net_transport_vr_col  (begc:endc, 1:nlevdecomp_full)) ; this%no3_net_transport_vr_col     (:,:) = spval

    allocate(this%externalp_to_primp_col    (begc:endc,1:nlevdecomp_full))  ; this%externalp_to_primp_col       (:,:) = spval
    allocate(this%externalp_to_labilep_col  (begc:endc,1:nlevdecomp_full))  ; this%externalp_to_labilep_col     (:,:) = spval
    allocate(this%externalp_to_solutionp_col(begc:endc,1:nlevdecomp_full))  ; this%externalp_to_solutionp_col   (:,:) = spval
    allocate(this%sminp_leached_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%sminp_leached_vr_col         (:,:) = nan
    allocate(this%sminp_net_transport_vr_col(begc:endc, 1:nlevdecomp_full)) ; this%sminp_net_transport_vr_col   (:,:) = spval

    allocate(this%plant_ndemand_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%plant_ndemand_vr_col         (:,:)  = nan
    allocate(this%plant_pdemand_vr_col      (begc:endc,1:nlevdecomp_full))  ; this%plant_pdemand_vr_col         (:,:)  = nan

    allocate(this%decomp_cpools_sourcesink_col      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))              ; this%decomp_cpools_sourcesink_col(:,:,:)= nan
    allocate(this%decomp_npools_sourcesink_col      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))              ; this%decomp_npools_sourcesink_col(:,:,:)= nan
    allocate(this%decomp_ppools_sourcesink_col      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))              ; this%decomp_ppools_sourcesink_col(:,:,:)= nan

    allocate(this%decomp_cascade_ctransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_ctransfer_vr_col(:,:,:)= nan
    allocate(this%decomp_cascade_ntransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_ntransfer_vr_col(:,:,:)= nan
    allocate(this%decomp_cascade_ptransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); this%decomp_cascade_ptransfer_vr_col(:,:,:)= nan

    allocate(this%gross_nmin_vr_col         (begc:endc,1:nlevdecomp_full))  ; this%gross_nmin_vr_col          (:,:) = nan
    allocate(this%net_nmin_vr_col           (begc:endc,1:nlevdecomp_full))  ; this%net_nmin_vr_col            (:,:) = nan

    allocate(this%potential_immob_vr_col    (begc:endc,1:nlevdecomp_full))  ; this%potential_immob_vr_col     (:,:) = nan
    allocate(this%actual_immob_vr_col       (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_vr_col        (:,:) = nan
    allocate(this%actual_immob_no3_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_no3_vr_col    (:,:) = nan
    allocate(this%actual_immob_nh4_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_nh4_vr_col    (:,:) = nan

    allocate(this%sminn_to_plant_vr_col     (begc:endc,1:nlevdecomp_full))  ; this%sminn_to_plant_vr_col      (:,:) = nan
    allocate(this%smin_no3_to_plant_vr_col  (begc:endc,1:nlevdecomp_full))  ; this%smin_no3_to_plant_vr_col   (:,:) = nan
    allocate(this%smin_nh4_to_plant_vr_col  (begc:endc,1:nlevdecomp_full))  ; this%smin_nh4_to_plant_vr_col   (:,:) = nan

    allocate(this%gross_pmin_vr_col         (begc:endc,1:nlevdecomp_full))  ; this%gross_pmin_vr_col          (:,:) = nan
    allocate(this%net_pmin_vr_col           (begc:endc,1:nlevdecomp_full))  ; this%net_pmin_vr_col            (:,:) = nan

    allocate(this%potential_immob_p_vr_col  (begc:endc,1:nlevdecomp_full))  ; this%potential_immob_p_vr_col   (:,:) = nan
    allocate(this%actual_immob_p_vr_col     (begc:endc,1:nlevdecomp_full))  ; this%actual_immob_p_vr_col      (:,:) = nan
    allocate(this%sminp_to_plant_vr_col     (begc:endc,1:nlevdecomp_full))  ; this%sminp_to_plant_vr_col      (:,:) = nan

    allocate(this%f_nit_vr_col              (begc:endc,1:nlevdecomp_full))  ; this%f_nit_vr_col               (:,:) = nan
    allocate(this%f_denit_vr_col            (begc:endc,1:nlevdecomp_full))  ; this%f_denit_vr_col             (:,:) = nan

    allocate(this%primp_to_labilep_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%primp_to_labilep_vr_col    (:,:) = nan
    allocate(this%labilep_to_secondp_vr_col (begc:endc,1:nlevdecomp_full))  ; this%labilep_to_secondp_vr_col  (:,:) = nan
    allocate(this%secondp_to_labilep_vr_col (begc:endc,1:nlevdecomp_full))  ; this%secondp_to_labilep_vr_col  (:,:) = nan
    allocate(this%secondp_to_occlp_vr_col   (begc:endc,1:nlevdecomp_full))  ; this%secondp_to_occlp_vr_col    (:,:) = nan

    allocate(this%decomp_cascade_hr_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))    ; this%decomp_cascade_hr_vr_col(:,:,:)= spval
    allocate(this%hr_vr_col                 (begc:endc,1:nlevdecomp_full))  ; this%hr_vr_col                  (:,:) = nan
    allocate(this%f_n2o_denit_vr_col        (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_denit_vr_col         (:,:) = nan
    allocate(this%f_n2o_nit_vr_col          (begc:endc,1:nlevdecomp_full))  ; this%f_n2o_nit_vr_col           (:,:) = nan

  end subroutine InitAllocate


end module clm_bgc_interface_data
