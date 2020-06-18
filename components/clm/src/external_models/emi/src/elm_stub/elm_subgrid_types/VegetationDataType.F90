module VegetationDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Vegetation data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval, spval
  use clm_varctl     , only : use_fates
  use clm_varpar      , only : nlevdecomp, nlevdecomp_full
  use clm_varpar      , only : nlevsno, nlevgrnd, nlevlak, nlevurb, nlevcan, crop_prog
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_state
    ! temperature variables
    real(r8), pointer :: t_veg              (:) => null() ! vegetation temperature (K)
    real(r8), pointer :: t_ref2m            (:) => null() ! 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_r          (:) => null() ! rural 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_u          (:) => null() ! urban 2 m height surface air temperature (K)
    ! temperature summary and accumulator variables
    real(r8), pointer :: t_a10              (:) => null() ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_a10min           (:) => null() ! 10-day running mean of min 2-m temperature
    real(r8), pointer :: t_a5min            (:) => null() ! 5-day running mean of min 2-m temperature
    real(r8), pointer :: t_ref2m_min        (:) => null() ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_r      (:) => null() ! daily minimum of average 2 m height surface air temperature - rural(K)
    real(r8), pointer :: t_ref2m_min_u      (:) => null() ! daily minimum of average 2 m height surface air temperature - urban (K)
    real(r8), pointer :: t_ref2m_max        (:) => null() ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max_r      (:) => null() ! daily maximum of average 2 m height surface air temperature - rural(K)
    real(r8), pointer :: t_ref2m_max_u      (:) => null() ! daily maximum of average 2 m height surface air temperature - urban (K)
    real(r8), pointer :: t_ref2m_min_inst   (:) => null() ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_min_inst_r (:) => null() ! instantaneous daily min of average 2 m height surface air temp - rural (K)
    real(r8), pointer :: t_ref2m_min_inst_u (:) => null() ! instantaneous daily min of average 2 m height surface air temp - urban (K)
    real(r8), pointer :: t_ref2m_max_inst   (:) => null() ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst_r (:) => null() ! instantaneous daily max of average 2 m height surface air temp - rural (K)
    real(r8), pointer :: t_ref2m_max_inst_u (:) => null() ! instantaneous daily max of average 2 m height surface air temp - urban (K)
    real(r8), pointer :: t_veg24            (:) => null() ! 24hr average vegetation temperature (K)
    real(r8), pointer :: t_veg240           (:) => null() ! 240hr average vegetation temperature (K)
    real(r8), pointer :: gdd0               (:) => null() ! growing degree-days base  0C from planting  (ddays)
    real(r8), pointer :: gdd8               (:) => null() ! growing degree-days base  8C from planting  (ddays)
    real(r8), pointer :: gdd10              (:) => null() ! growing degree-days base 10C from planting  (ddays)
    real(r8), pointer :: gdd020             (:) => null() ! 20-year average of gdd0                     (ddays)
    real(r8), pointer :: gdd820             (:) => null() ! 20-year average of gdd8                     (ddays)
    real(r8), pointer :: gdd1020            (:) => null() ! 20-year average of gdd10                    (ddays)
    ! temperature-related variables
    real(r8), pointer :: thm                (:) => null() ! intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
    real(r8), pointer :: emv                (:) => null() ! vegetation emissivity (unitless)

    contains
    procedure, public :: Init    => veg_es_init
  end type vegetation_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_state
    ! Note: All units given as kg in this data type imply kg of H2O
    real(r8), pointer :: h2ocan       (:) => null() ! canopy water (kg/m2)
    real(r8), pointer :: q_ref2m      (:) => null() ! 2 m height surface specific humidity (kg H2O/kg moist air)
    real(r8), pointer :: rh_ref2m     (:) => null() ! 2 m height surface relative humidity (%)
    real(r8), pointer :: rh_ref2m_r   (:) => null() ! 2 m height surface relative humidity - rural (%)
    real(r8), pointer :: rh_ref2m_u   (:) => null() ! 2 m height surface relative humidity - urban (%)
    real(r8), pointer :: rh_af        (:) => null() ! fractional humidity of canopy air (dimensionless)
    real(r8), pointer :: fwet         (:) => null() ! canopy fraction that is wet (0 to 1)
    real(r8), pointer :: fdry         (:) => null() ! canopy fraction of foliage that is green and dry (0 to 1)
    real(r8), pointer :: begwb        (:) => null() ! water mass begining of the time step
    real(r8), pointer :: endwb        (:) => null() ! water mass end of the time step
    real(r8), pointer :: errh2o       (:) => null() ! water conservation error (mm H2O)
  contains
    procedure, public :: Init    => veg_ws_init
  end type vegetation_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_state
    integer           :: species                          ! c12, c13, c14
    real(r8), pointer :: leafc              (:) => null() ! (gC/m2) leaf C
    real(r8), pointer :: leafc_storage      (:) => null() ! (gC/m2) leaf C storage
    real(r8), pointer :: leafc_xfer         (:) => null() ! (gC/m2) leaf C transfer
    real(r8), pointer :: frootc             (:) => null() ! (gC/m2) fine root C
    real(r8), pointer :: frootc_storage     (:) => null() ! (gC/m2) fine root C storage
    real(r8), pointer :: frootc_xfer        (:) => null() ! (gC/m2) fine root C transfer
    real(r8), pointer :: livestemc          (:) => null() ! (gC/m2) live stem C
    real(r8), pointer :: livestemc_storage  (:) => null() ! (gC/m2) live stem C storage
    real(r8), pointer :: livestemc_xfer     (:) => null() ! (gC/m2) live stem C transfer
    real(r8), pointer :: deadstemc          (:) => null() ! (gC/m2) dead stem C
    real(r8), pointer :: deadstemc_storage  (:) => null() ! (gC/m2) dead stem C storage
    real(r8), pointer :: deadstemc_xfer     (:) => null() ! (gC/m2) dead stem C transfer
    real(r8), pointer :: livecrootc         (:) => null() ! (gC/m2) live coarse root C
    real(r8), pointer :: livecrootc_storage (:) => null() ! (gC/m2) live coarse root C storage
    real(r8), pointer :: livecrootc_xfer    (:) => null() ! (gC/m2) live coarse root C transfer
    real(r8), pointer :: deadcrootc         (:) => null() ! (gC/m2) dead coarse root C
    real(r8), pointer :: deadcrootc_storage (:) => null() ! (gC/m2) dead coarse root C storage
    real(r8), pointer :: deadcrootc_xfer    (:) => null() ! (gC/m2) dead coarse root C transfer
    real(r8), pointer :: gresp_storage      (:) => null() ! (gC/m2) growth respiration storage
    real(r8), pointer :: gresp_xfer         (:) => null() ! (gC/m2) growth respiration transfer
    real(r8), pointer :: cpool              (:) => null() ! (gC/m2) temporary photosynthate C pool
    real(r8), pointer :: xsmrpool           (:) => null() ! (gC/m2) abstract C pool to meet excess MR demand
    real(r8), pointer :: grainc             (:) => null() ! (gC/m2) grain C (crop model)
    real(r8), pointer :: grainc_storage     (:) => null() ! (gC/m2) grain C storage (crop model)
    real(r8), pointer :: grainc_xfer        (:) => null() ! (gC/m2) grain C transfer (crop model)
    real(r8), pointer :: ctrunc             (:) => null() ! (gC/m2) patch-level sink for C truncation
    real(r8), pointer :: woodc              (:) => null() ! (gC/m2) wood C
    real(r8), pointer :: leafcmax           (:) => null() ! (gC/m2) ann max leaf C
    real(r8), pointer :: cropseedc_deficit  (:) => null() ! (gC/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
    real(r8), pointer :: dispvegc           (:) => null() ! (gC/m2) displayed veg carbon, excluding storage and cpool
    real(r8), pointer :: storvegc           (:) => null() ! (gC/m2) stored vegetation carbon, excluding cpool
    real(r8), pointer :: totvegc            (:) => null() ! (gC/m2) total vegetation carbon, excluding cpool
    real(r8), pointer :: totpftc            (:) => null() ! (gC/m2) total patch-level carbon, including cpool
    real(r8), pointer :: totvegc_abg        (:) => null() ! (gC/m2) total above vegetation carbon, excluding cpool
    real(r8), pointer :: begcb              (:) => null() ! patch carbon mass, beginning of time step (gC/m**2)
    real(r8), pointer :: endcb              (:) => null() ! patch carbon mass, end of time step (gC/m**2)
    real(r8), pointer :: errcb              (:) => null() ! patch carbon balance error for the timestep (gC/m**2)
  contains
    procedure, public :: Init     => veg_cs_init
  end type vegetation_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_state
    real(r8), pointer :: leafn                  (:)   => null()  ! (gN/m2) leaf N 
    real(r8), pointer :: leafn_storage          (:)   => null()  ! (gN/m2) leaf N storage
    real(r8), pointer :: leafn_xfer             (:)   => null()  ! (gN/m2) leaf N transfer
    real(r8), pointer :: frootn                 (:)   => null()  ! (gN/m2) fine root N
    real(r8), pointer :: frootn_storage         (:)   => null()  ! (gN/m2) fine root N storage
    real(r8), pointer :: frootn_xfer            (:)   => null()  ! (gN/m2) fine root N transfer
    real(r8), pointer :: livestemn              (:)   => null()  ! (gN/m2) live stem N
    real(r8), pointer :: livestemn_storage      (:)   => null()  ! (gN/m2) live stem N storage
    real(r8), pointer :: livestemn_xfer         (:)   => null()  ! (gN/m2) live stem N transfer
    real(r8), pointer :: deadstemn              (:)   => null()  ! (gN/m2) dead stem N
    real(r8), pointer :: deadstemn_storage      (:)   => null()  ! (gN/m2) dead stem N storage
    real(r8), pointer :: deadstemn_xfer         (:)   => null()  ! (gN/m2) dead stem N transfer
    real(r8), pointer :: livecrootn             (:)   => null()  ! (gN/m2) live coarse root N
    real(r8), pointer :: livecrootn_storage     (:)   => null()  ! (gN/m2) live coarse root N storage
    real(r8), pointer :: livecrootn_xfer        (:)   => null()  ! (gN/m2) live coarse root N transfer
    real(r8), pointer :: deadcrootn             (:)   => null()  ! (gN/m2) dead coarse root N
    real(r8), pointer :: deadcrootn_storage     (:)   => null()  ! (gN/m2) dead coarse root N storage
    real(r8), pointer :: deadcrootn_xfer        (:)   => null()  ! (gN/m2) dead coarse root N transfer
    real(r8), pointer :: retransn               (:)   => null()  ! (gN/m2) plant pool of retranslocated N
    real(r8), pointer :: npool                  (:)   => null()  ! (gN/m2) temporary plant N pool
    real(r8), pointer :: ntrunc                 (:)   => null()  ! (gN/m2) pft-level sink for N truncation
    real(r8), pointer :: plant_n_buffer         (:)   => null()  ! (gN/m2) pft-level abstract N storage
    real(r8), pointer :: grainn                 (:)   => null()  ! (gN/m2) grain N (crop)
    real(r8), pointer :: grainn_storage         (:)   => null()  ! (gN/m2) grain N storage (crop)
    real(r8), pointer :: grainn_xfer            (:)   => null()  ! (gN/m2) grain N transfer (crop)
    real(r8), pointer :: cropseedn_deficit      (:)   => null()  ! (gN/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid     
    real(r8), pointer :: dispvegn               (:)   => null()  ! (gN/m2) displayed veg nitrogen, excluding storage
    real(r8), pointer :: storvegn               (:)   => null()  ! (gN/m2) stored vegetation nitrogen
    real(r8), pointer :: totvegn                (:)   => null()  ! (gN/m2) total vegetation nitrogen
    real(r8), pointer :: totpftn                (:)   => null()  ! (gN/m2) total pft-level nitrogen
    real(r8), pointer :: begnb                  (:)   => null()  ! (gN/m2) nitrogen mass, beginning of time step 
    real(r8), pointer :: endnb                  (:)   => null()  ! (gN/m2) nitrogen mass, end of time step 
    real(r8), pointer :: errnb                  (:)   => null()  ! (gN/m2) nitrogen balance error for the timestep 
    ! variables supporting dynamic C/N/P allocation cost/benefit analysis
    real(r8), pointer :: npimbalance            (:)   => null()
    real(r8), pointer :: pnup_pfrootc           (:)   => null()
    real(r8), pointer :: ppup_pfrootc           (:)   => null()
    real(r8), pointer :: ptlai_pleafc           (:)   => null()
    real(r8), pointer :: ppsnsun_ptlai          (:)   => null()
    real(r8), pointer :: ppsnsun_pleafn         (:)   => null()
    real(r8), pointer :: ppsnsun_pleafp         (:)   => null()
    real(r8), pointer :: plmrsun_ptlai          (:)   => null()
    real(r8), pointer :: plmrsun_pleafn         (:)   => null()
    real(r8), pointer :: plaisun_ptlai          (:)   => null()
    real(r8), pointer :: ppsnsha_ptlai          (:)   => null()
    real(r8), pointer :: ppsnsha_pleafn         (:)   => null()
    real(r8), pointer :: ppsnsha_pleafp         (:)   => null()
    real(r8), pointer :: plmrsha_ptlai          (:)   => null()
    real(r8), pointer :: plmrsha_pleafn         (:)   => null()
    real(r8), pointer :: plaisha_ptlai          (:)   => null()
    real(r8), pointer :: benefit_pgpp_pleafc    (:)   => null()  ! partial gpp / partial leaf carbon (used by symbiotic n2 fixation and dynamic allocation)
    real(r8), pointer :: benefit_pgpp_pleafn    (:)   => null()  ! partial gpp / partial leaf nitrogen (used by phosphatase activity and dynamic allocation)
    real(r8), pointer :: benefit_pgpp_pleafp    (:)   => null()  ! partial gpp / partial leaf phosphorus (used by phosphatase activity and dynamic allocation)
    real(r8), pointer :: cost_pgpp_pfrootc      (:)   => null()  ! partial gpp /  partial fine root carbon (used by dynamic allocation)
    real(r8), pointer :: cost_plmr_pleafc       (:)   => null()  ! partial maintenance respiration /  partial leaf carbon (used by dynamic allocation)
    real(r8), pointer :: cost_plmr_pleafn       (:)   => null()  ! partial maintenance respiration /  partial leaf nitrogen (used by dynamic allocation)
    ! variables supporting multi-layer canopy
    real(r8), pointer :: ppsn_ptlai_z           (:,:) => null()
    real(r8), pointer :: ppsn_pleafn_z          (:,:) => null()
    real(r8), pointer :: ppsn_pleafp_z          (:,:) => null()
    real(r8), pointer :: ppsn_ptlai_z_vcmax     (:,:) => null()
    real(r8), pointer :: ppsn_pleafn_z_vcmax    (:,:) => null()
    real(r8), pointer :: ppsn_pleafp_z_vcmax    (:,:) => null()
    real(r8), pointer :: ppsn_ptlai_z_jmax      (:,:) => null()
    real(r8), pointer :: ppsn_pleafn_z_jmax     (:,:) => null()
    real(r8), pointer :: ppsn_pleafp_z_jmax     (:,:) => null()
    real(r8), pointer :: ppsn_ptlai_z_tpu       (:,:) => null()
    real(r8), pointer :: ppsn_pleafn_z_tpu      (:,:) => null()
    real(r8), pointer :: ppsn_pleafp_z_tpu      (:,:) => null()
    real(r8), pointer :: plmr_ptlai_z           (:,:) => null()
    real(r8), pointer :: plmr_pleafn_z          (:,:) => null()
     
  contains
    procedure, public :: Init      => veg_ns_init
  end type vegetation_nitrogen_state
 
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_state
    real(r8), pointer :: grainp                 (:)     ! (gP/m2) grain P (crop)
    real(r8), pointer :: grainp_storage         (:)     ! (gP/m2) grain P storage (crop)
    real(r8), pointer :: grainp_xfer            (:)     ! (gP/m2) grain P transfer (crop)
    real(r8), pointer :: leafp                  (:)     ! (gP/m2) leaf P 
    real(r8), pointer :: leafp_storage          (:)     ! (gP/m2) leaf P storage
    real(r8), pointer :: leafp_xfer             (:)     ! (gP/m2) leaf P transfer
    real(r8), pointer :: frootp                 (:)     ! (gP/m2) fine root P
    real(r8), pointer :: frootp_storage         (:)     ! (gP/m2) fine root P storage
    real(r8), pointer :: frootp_xfer            (:)     ! (gP/m2) fine root P transfer
    real(r8), pointer :: livestemp              (:)     ! (gP/m2) live stem P
    real(r8), pointer :: livestemp_storage      (:)     ! (gP/m2) live stem P storage
    real(r8), pointer :: livestemp_xfer         (:)     ! (gP/m2) live stem P transfer
    real(r8), pointer :: deadstemp              (:)     ! (gP/m2) dead stem P
    real(r8), pointer :: deadstemp_storage      (:)     ! (gP/m2) dead stem P storage
    real(r8), pointer :: deadstemp_xfer         (:)     ! (gP/m2) dead stem P transfer
    real(r8), pointer :: livecrootp             (:)     ! (gP/m2) live coarse root P
    real(r8), pointer :: livecrootp_storage     (:)     ! (gP/m2) live coarse root P storage
    real(r8), pointer :: livecrootp_xfer        (:)     ! (gP/m2) live coarse root P transfer
    real(r8), pointer :: deadcrootp             (:)     ! (gP/m2) dead coarse root P
    real(r8), pointer :: deadcrootp_storage     (:)     ! (gP/m2) dead coarse root P storage
    real(r8), pointer :: deadcrootp_xfer        (:)     ! (gP/m2) dead coarse root P transfer
    real(r8), pointer :: retransp               (:)     ! (gP/m2) plant pool of retranslocated P
    real(r8), pointer :: ppool                  (:)     ! (gP/m2) temporary plant P pool
    real(r8), pointer :: ptrunc                 (:)     ! (gP/m2) pft-level sink for P truncation
    real(r8), pointer :: plant_p_buffer         (:)     ! (gP/m2) pft-level abstract p storage
    real(r8), pointer :: cropseedp_deficit      (:)     ! (gP/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
    real(r8), pointer :: dispvegp               (:)     ! (gP/m2) displayed veg phosphorus, excluding storage
    real(r8), pointer :: storvegp               (:)     ! (gP/m2) stored vegetation phosphorus
    real(r8), pointer :: totvegp                (:)     ! (gP/m2) total vegetation phosphorus
    real(r8), pointer :: totpftp                (:)     ! (gP/m2) total pft-level phosphorus
    real(r8), pointer :: begpb                  (:)     ! phosphorus mass, beginning of time step (gP/m**2)
    real(r8), pointer :: endpb                  (:)     ! phosphorus mass, end of time step (gP/m**2)
    real(r8), pointer :: errpb                  (:)     ! phosphorus balance error for the timestep (gP/m**2)
  contains
    procedure, public :: Init      => veg_ps_init
  end type vegetation_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_flux
    real(r8), pointer :: eflx_sh_grnd      (:) => null() ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_veg       (:) => null() ! sensible heat flux from leaves (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_snow      (:) => null() ! sensible heat flux from snow (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_soil      (:) => null() ! sensible heat flux from soil  (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_h2osfc    (:) => null() ! sensible heat flux from surface water (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot       (:) => null() ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_u     (:) => null() ! urban total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_r     (:) => null() ! rural total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot       (:) => null() ! total latent heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot_u     (:) => null() ! urban total latent heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot_r     (:) => null() ! rural total latent heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vegt      (:) => null() ! transpiration heat flux from veg (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vege      (:) => null() ! evaporation heat flux from veg (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_grnd      (:) => null() ! evaporation heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_soil_grnd    (:) => null() ! soil heat flux (W/m**2) [+ = into soil] 
    real(r8), pointer :: eflx_soil_grnd_u  (:) => null() ! urban soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_soil_grnd_r  (:) => null() ! rural soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_lwrad_net    (:) => null() ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_net_r  (:) => null() ! rural net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_net_u  (:) => null() ! urban net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_out    (:) => null() ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_out_r  (:) => null() ! rural emitted infrared (longwave) rad (W/m**2)
    real(r8), pointer :: eflx_lwrad_out_u  (:) => null() ! urban emitted infrared (longwave) rad (W/m**2)
    real(r8), pointer :: eflx_gnet         (:) => null() ! net heat flux into ground  (W/m**2)
    real(r8), pointer :: eflx_grnd_lake    (:) => null() ! net heat flux into lake / snow surface, excluding light transmission (W/m**2)
    real(r8), pointer :: eflx_anthro       (:) => null() ! total anthropogenic heat flux (W/m**2)
    real(r8), pointer :: eflx_traffic      (:) => null() ! traffic sensible heat flux (W/m2)
    real(r8), pointer :: eflx_wasteheat    (:) => null() ! sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac (:) => null() ! sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: dlrad             (:) => null() ! downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad             (:) => null() ! upward longwave radiation above the canopy [W/m2]
    ! Wind Stress
    real(r8), pointer :: taux              (:) => null() ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy              (:) => null() ! wind (shear) stress: n-s (kg/m/s**2)
    ! Derivatives of energy fluxes
    real(r8), pointer :: dgnetdT           (:) => null() ! derivative of net ground heat flux wrt soil temp  (W/m**2 K)
    real(r8), pointer :: netrad            (:) => null() ! col net radiation (W/m**2) [+ = to sfc]
    real(r8), pointer :: cgrnd             (:) => null() ! col deriv. of soil energy flux wrt to soil temp [W/m2/k]
    real(r8), pointer :: cgrndl            (:) => null() ! col deriv. of soil latent heat flux wrt soil temp  [W/m**2/k]
    real(r8), pointer :: cgrnds            (:) => null() ! col deriv. of soil sensible heat flux wrt soil temp [W/m2/k]
    ! Balance Checks
    real(r8), pointer :: errsoi            (:) => null() ! soil/lake energy conservation error   (W/m**2)
    real(r8), pointer :: errseb            (:) => null() ! surface energy conservation error     (W/m**2)
    real(r8), pointer :: errsol            (:) => null() ! solar radiation conservation error    (W/m**2)
    real(r8), pointer :: errlon            (:) => null() ! longwave radiation conservation error (W/m**2)

  contains
    procedure, public :: Init    => veg_ef_init
  end type vegetation_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_flux
    real(r8), pointer :: qflx_prec_grnd     (:)   => null() ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd     (:)   => null() ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd     (:)   => null() ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow      (:)   => null() ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_soi      (:)   => null() ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_veg      (:)   => null() ! vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can      (:)   => null() ! evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot      (:)   => null() ! pft_qflx_evap_soi + pft_qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd     (:)   => null() ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq     (:)   => null() ! excess rainfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_ice     (:)   => null() ! excess snowfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_tran_veg      (:)   => null() ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_dew_snow      (:)   => null() ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd      (:)   => null() ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_prec_intr     (:)   => null() ! interception of precipitation [mm/s]
    real(r8), pointer :: qflx_ev_snow       (:)   => null() ! evaporation heat flux from snow       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_soil       (:)   => null() ! evaporation heat flux from soil       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_h2osfc     (:)   => null() ! evaporation heat flux from soil       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_rootsoi_frac  (:,:) => null() !  
    real(r8), pointer :: qflx_irrig         (:)   => null() ! irrigation flux (mm H2O/s)
    real(r8), pointer :: irrig_rate         (:)   => null() ! current irrigation rate [mm/s]
    integer , pointer :: n_irrig_steps_left (:)   => null() ! number of time steps for which we still need to irrigate today (if 0, ignore)

  contains
    procedure, public :: Init    => veg_wf_init
  end type vegetation_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_flux
    ! gap mortality fluxes
    real(r8), pointer :: m_leafc_to_litter                   (:) => null()    ! leaf C mortality (gC/m2/s)
    real(r8), pointer :: m_leafc_storage_to_litter           (:) => null()    ! leaf C storage mortality (gC/m2/s)
    real(r8), pointer :: m_leafc_xfer_to_litter              (:) => null()    ! leaf C transfer mortality (gC/m2/s)
    real(r8), pointer :: m_frootc_to_litter                  (:) => null()    ! fine root C mortality (gC/m2/s)
    real(r8), pointer :: m_frootc_storage_to_litter          (:) => null()    ! fine root C storage mortality (gC/m2/s)
    real(r8), pointer :: m_frootc_xfer_to_litter             (:) => null()    ! fine root C transfer mortality (gC/m2/s)
    real(r8), pointer :: m_livestemc_to_litter               (:) => null()    ! live stem C mortality (gC/m2/s)
    real(r8), pointer :: m_livestemc_storage_to_litter       (:) => null()    ! live stem C storage mortality (gC/m2/s)
    real(r8), pointer :: m_livestemc_xfer_to_litter          (:) => null()    ! live stem C transfer mortality (gC/m2/s)
    real(r8), pointer :: m_deadstemc_to_litter               (:) => null()    ! dead stem C mortality (gC/m2/s)
    real(r8), pointer :: m_deadstemc_storage_to_litter       (:) => null()    ! dead stem C storage mortality (gC/m2/s)
    real(r8), pointer :: m_deadstemc_xfer_to_litter          (:) => null()    ! dead stem C transfer mortality (gC/m2/s)
    real(r8), pointer :: m_livecrootc_to_litter              (:) => null()    ! live coarse root C mortality (gC/m2/s)
    real(r8), pointer :: m_livecrootc_storage_to_litter      (:) => null()    ! live coarse root C storage mortality (gC/m2/s)
    real(r8), pointer :: m_livecrootc_xfer_to_litter         (:) => null()    ! live coarse root C transfer mortality (gC/m2/s)
    real(r8), pointer :: m_deadcrootc_to_litter              (:) => null()    ! dead coarse root C mortality (gC/m2/s)
    real(r8), pointer :: m_deadcrootc_storage_to_litter      (:) => null()    ! dead coarse root C storage mortality (gC/m2/s)
    real(r8), pointer :: m_deadcrootc_xfer_to_litter         (:) => null()    ! dead coarse root C transfer mortality (gC/m2/s)
    real(r8), pointer :: m_gresp_storage_to_litter           (:) => null()    ! growth respiration storage mortality (gC/m2/s)
    real(r8), pointer :: m_gresp_xfer_to_litter              (:) => null()    ! growth respiration transfer mortality (gC/m2/s)
    real(r8), pointer :: m_cpool_to_litter                   (:) => null()    ! plant storage C pool to litter (gC/m2/s)
                                                                 
    ! harvest mortality fluxes                                   
    real(r8), pointer :: hrv_leafc_to_litter                 (:) => null()    ! leaf C harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_leafc_storage_to_litter         (:) => null()    ! leaf C storage harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_leafc_xfer_to_litter            (:) => null()    ! leaf C transfer harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_frootc_to_litter                (:) => null()    ! fine root C harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_frootc_storage_to_litter        (:) => null()    ! fine root C storage harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_frootc_xfer_to_litter           (:) => null()    ! fine root C transfer harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_livestemc_to_litter             (:) => null()    ! live stem C harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_livestemc_storage_to_litter     (:) => null()    ! live stem C storage harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_livestemc_xfer_to_litter        (:) => null()    ! live stem C transfer harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_deadstemc_to_prod10c            (:) => null()    ! dead stem C harvest to 10-year product pool (gC/m2/s)
    real(r8), pointer :: hrv_deadstemc_to_prod100c           (:) => null()    ! dead stem C harvest to 100-year product pool (gC/m2/s)
    real(r8), pointer :: hrv_deadstemc_storage_to_litter     (:) => null()    ! dead stem C storage harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_deadstemc_xfer_to_litter        (:) => null()    ! dead stem C transfer harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_livecrootc_to_litter            (:) => null()    ! live coarse root C harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_livecrootc_storage_to_litter    (:) => null()    ! live coarse root C storage harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_livecrootc_xfer_to_litter       (:) => null()    ! live coarse root C transfer harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_deadcrootc_to_litter            (:) => null()    ! dead coarse root C harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_deadcrootc_storage_to_litter    (:) => null()    ! dead coarse root C storage harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_deadcrootc_xfer_to_litter       (:) => null()    ! dead coarse root C transfer harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_gresp_storage_to_litter         (:) => null()    ! growth respiration storage harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_gresp_xfer_to_litter            (:) => null()    ! growth respiration transfer harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_xsmrpool_to_atm                 (:) => null()    ! excess MR pool harvest mortality (gC/m2/s)
    real(r8), pointer :: hrv_cpool_to_litter                 (:) => null()    ! Harvest cpool to litter (gC/m2/s)
    ! crop harvest                                               
    real(r8), pointer :: hrv_leafc_to_prod1c                 (:) => null()    ! crop leafc harvested (gC/m2/s)
    real(r8), pointer :: hrv_livestemc_to_prod1c             (:) => null()    ! crop stemc harvested (gC/m2/s)
    real(r8), pointer :: hrv_grainc_to_prod1c                (:) => null()    ! crop grain harvested (gC/m2/s)
    real(r8), pointer :: hrv_cropc_to_prod1c                 (:) => null()    ! total amount of crop C harvested (gC/m2/s)
                                                                
    ! fire C fluxes                                             
    real(r8), pointer :: m_leafc_to_fire                     (:) => null()    ! (gC/m2/s) fire C emissions from leafc 
    real(r8), pointer :: m_leafc_storage_to_fire             (:) => null()    ! (gC/m2/s) fire C emissions from leafc_storage             
    real(r8), pointer :: m_leafc_xfer_to_fire                (:) => null()    ! (gC/m2/s) fire C emissions from leafc_xfer
    real(r8), pointer :: m_livestemc_to_fire                 (:) => null()    ! (gC/m2/s) fire C emissions from livestemc
    real(r8), pointer :: m_livestemc_storage_to_fire         (:) => null()    ! (gC/m2/s) fire C emissions from livestemc_storage       
    real(r8), pointer :: m_livestemc_xfer_to_fire            (:) => null()    ! (gC/m2/s) fire C emissions from livestemc_xfer
    real(r8), pointer :: m_deadstemc_to_fire                 (:) => null()    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(r8), pointer :: m_deadstemc_storage_to_fire         (:) => null()    ! (gC/m2/s) fire C emissions from deadstemc_storage         
    real(r8), pointer :: m_deadstemc_xfer_to_fire            (:) => null()    ! (gC/m2/s) fire C emissions from deadstemc_xfer
    real(r8), pointer :: m_frootc_to_fire                    (:) => null()    ! (gC/m2/s) fire C emissions from frootc
    real(r8), pointer :: m_frootc_storage_to_fire            (:) => null()    ! (gC/m2/s) fire C emissions from frootc_storage
    real(r8), pointer :: m_frootc_xfer_to_fire               (:) => null()    ! (gC/m2/s) fire C emissions from frootc_xfer
    real(r8), pointer :: m_livecrootc_to_fire                (:) => null()    ! (gC/m2/s) fire C emissions from livecrootc
    real(r8), pointer :: m_livecrootc_storage_to_fire        (:) => null()    ! (gC/m2/s) fire C emissions from livecrootc_storage     
    real(r8), pointer :: m_livecrootc_xfer_to_fire           (:) => null()    ! (gC/m2/s) fire C emissions from livecrootc_xfer
    real(r8), pointer :: m_deadcrootc_to_fire                (:) => null()    ! (gC/m2/s) fire C emissions from deadcrootc
    real(r8), pointer :: m_deadcrootc_storage_to_fire        (:) => null()    ! (gC/m2/s) fire C emissions from deadcrootc_storage 
    real(r8), pointer :: m_deadcrootc_xfer_to_fire           (:) => null()    ! (gC/m2/s) fire C emissions from deadcrootc_xfer
    real(r8), pointer :: m_gresp_storage_to_fire             (:) => null()    ! (gC/m2/s) fire C emissions from gresp_storage 
    real(r8), pointer :: m_gresp_xfer_to_fire                (:) => null()    ! (gC/m2/s) fire C emissions from gresp_xfer
    real(r8), pointer :: m_cpool_to_fire                     (:) => null()    ! (gC/m2/s) fire C emissions from cpool
    real(r8), pointer :: m_leafc_to_litter_fire              (:) => null()    ! (gC/m2/s) from leafc to litter c due to fire
    real(r8), pointer :: m_leafc_storage_to_litter_fire      (:) => null()    ! (gC/m2/s) from leafc_storage to litter C  due to fire               
    real(r8), pointer :: m_leafc_xfer_to_litter_fire         (:) => null()    ! (gC/m2/s) from leafc_xfer to litter C  due to fire               
    real(r8), pointer :: m_livestemc_to_litter_fire          (:) => null()    ! (gC/m2/s) from livestemc to litter C  due to fire               
    real(r8), pointer :: m_livestemc_storage_to_litter_fire  (:) => null()    ! (gC/m2/s) from livestemc_storage to litter C due to fire      
    real(r8), pointer :: m_livestemc_xfer_to_litter_fire     (:) => null()    ! (gC/m2/s) from livestemc_xfer to litter C due to fire      
    real(r8), pointer :: m_livestemc_to_deadstemc_fire       (:) => null()    ! (gC/m2/s) from livestemc to deadstemc due to fire       
    real(r8), pointer :: m_deadstemc_to_litter_fire          (:) => null()    ! (gC/m2/s) from deadstemc to litter C due to fire      
    real(r8), pointer :: m_deadstemc_storage_to_litter_fire  (:) => null()    ! (gC/m2/s) from deadstemc_storage to litter C due to fire               
    real(r8), pointer :: m_deadstemc_xfer_to_litter_fire     (:) => null()    ! (gC/m2/s) from deadstemc_xfer to litter C due to fire               
    real(r8), pointer :: m_frootc_to_litter_fire             (:) => null()    ! (gC/m2/s) from frootc to litter C due to fire               
    real(r8), pointer :: m_frootc_storage_to_litter_fire     (:) => null()    ! (gC/m2/s) from frootc_storage to litter C due to fire               
    real(r8), pointer :: m_frootc_xfer_to_litter_fire        (:) => null()    ! (gC/m2/s) from frootc_xfer to litter C due to fire               
    real(r8), pointer :: m_livecrootc_to_litter_fire         (:) => null()    ! (gC/m2/s) from livecrootc to litter C due to fire                     
    real(r8), pointer :: m_livecrootc_storage_to_litter_fire (:) => null()    ! (gC/m2/s) from livecrootc_storage to litter C due to fire                     
    real(r8), pointer :: m_livecrootc_xfer_to_litter_fire    (:) => null()    ! (gC/m2/s) from livecrootc_xfer to litter C due to fire                     
    real(r8), pointer :: m_livecrootc_to_deadcrootc_fire     (:) => null()    ! (gC/m2/s) from livecrootc to deadstemc due to fire        
    real(r8), pointer :: m_deadcrootc_to_litter_fire         (:) => null()    ! (gC/m2/s) from deadcrootc to litter C due to fire                       
    real(r8), pointer :: m_deadcrootc_storage_to_litter_fire (:) => null()    ! (gC/m2/s) from deadcrootc_storage to litter C due to fire                       
    real(r8), pointer :: m_deadcrootc_xfer_to_litter_fire    (:) => null()    ! (gC/m2/s) from deadcrootc_xfer to litter C due to fire                       
    real(r8), pointer :: m_gresp_storage_to_litter_fire      (:) => null()    ! (gC/m2/s) from gresp_storage to litter C due to fire                       
    real(r8), pointer :: m_gresp_xfer_to_litter_fire         (:) => null()    ! (gC/m2/s) from gresp_xfer to litter C due to fire                       
    real(r8), pointer :: m_cpool_to_litter_fire              (:) => null()    ! (gC/m2/s) from cpool to litter C due to fire              
                                                                 
    ! phenology fluxes from transfer pools                       
    real(r8), pointer :: grainc_xfer_to_grainc               (:) => null()    ! grain C growth from storage for prognostic crop(gC/m2/s)
    real(r8), pointer :: leafc_xfer_to_leafc                 (:) => null()    ! leaf C growth from storage (gC/m2/s)
    real(r8), pointer :: frootc_xfer_to_frootc               (:) => null()    ! fine root C growth from storage (gC/m2/s)
    real(r8), pointer :: livestemc_xfer_to_livestemc         (:) => null()    ! live stem C growth from storage (gC/m2/s)
    real(r8), pointer :: deadstemc_xfer_to_deadstemc         (:) => null()    ! dead stem C growth from storage (gC/m2/s)
    real(r8), pointer :: livecrootc_xfer_to_livecrootc       (:) => null()    ! live coarse root C growth from storage (gC/m2/s)
    real(r8), pointer :: deadcrootc_xfer_to_deadcrootc       (:) => null()    ! dead coarse root C growth from storage (gC/m2/s)
                                                                 
    ! leaf and fine root litterfall fluxes                          
    real(r8), pointer :: leafc_to_litter                     (:) => null()    ! leaf C litterfall (gC/m2/s)
    real(r8), pointer :: frootc_to_litter                    (:) => null()    ! fine root C litterfall (gC/m2/s)
    real(r8), pointer :: livestemc_to_litter                 (:) => null()    ! live stem C litterfall (gC/m2/s)
    real(r8), pointer :: grainc_to_food                      (:) => null()    ! grain C to food for prognostic crop(gC/m2/s)
                                                                 
    ! maintenance respiration fluxes                             
    real(r8), pointer :: leaf_mr                             (:) => null()    ! leaf maintenance respiration (gC/m2/s)
    real(r8), pointer :: froot_mr                            (:) => null()    ! fine root maintenance respiration (gC/m2/s)
    real(r8), pointer :: livestem_mr                         (:) => null()    ! live stem maintenance respiration (gC/m2/s)
    real(r8), pointer :: livecroot_mr                        (:) => null()    ! live coarse root maintenance respiration (gC/m2/s)
    real(r8), pointer :: grain_mr                            (:) => null()    ! crop grain or organs maint. respiration (gC/m2/s)
    real(r8), pointer :: leaf_curmr                          (:) => null()    ! leaf maintenance respiration from current GPP (gC/m2/s)
    real(r8), pointer :: froot_curmr                         (:) => null()    ! fine root maintenance respiration from current GPP (gC/m2/s)
    real(r8), pointer :: livestem_curmr                      (:) => null()    ! live stem maintenance respiration from current GPP (gC/m2/s)
    real(r8), pointer :: livecroot_curmr                     (:) => null()    ! live coarse root maintenance respiration from current GPP (gC/m2/s)
    real(r8), pointer :: grain_curmr                         (:) => null()    ! crop grain or organs maint. respiration from current GPP (gC/m2/s)
    real(r8), pointer :: leaf_xsmr                           (:) => null()    ! leaf maintenance respiration from storage (gC/m2/s)
    real(r8), pointer :: froot_xsmr                          (:) => null()    ! fine root maintenance respiration from storage (gC/m2/s)
    real(r8), pointer :: livestem_xsmr                       (:) => null()    ! live stem maintenance respiration from storage (gC/m2/s)
    real(r8), pointer :: livecroot_xsmr                      (:) => null()    ! live coarse root maintenance respiration from storage (gC/m2/s)
    real(r8), pointer :: grain_xsmr                          (:) => null()    ! crop grain or organs maint. respiration from storage (gC/m2/s)
    !turnover of excess carbon                                   
    real(r8), pointer :: xr                                  (:) => null()    ! respiration from excess carbon cpool (gC/m2/s)
                                                                 
    ! photosynthesis fluxes                                      
    real(r8), pointer :: psnsun_to_cpool                     (:) => null()    ! C fixation from sunlit canopy (gC/m2/s)
    real(r8), pointer :: psnshade_to_cpool                   (:) => null()    ! C fixation from shaded canopy (gC/m2/s)
                                                                 
    ! allocation fluxes, from current GPP                        
    real(r8), pointer :: cpool_to_xsmrpool                   (:) => null()    ! allocation to maintenance respiration storage pool (gC/m2/s)
    real(r8), pointer :: cpool_to_grainc                     (:) => null()    ! allocation to grain C for prognostic crop(gC/m2/s)
    real(r8), pointer :: cpool_to_grainc_storage             (:) => null()    ! allocation to grain C storage for prognostic crop(gC/m2/s)
    real(r8), pointer :: cpool_to_leafc                      (:) => null()    ! allocation to leaf C (gC/m2/s)
    real(r8), pointer :: cpool_to_leafc_storage              (:) => null()    ! allocation to leaf C storage (gC/m2/s)
    real(r8), pointer :: cpool_to_frootc                     (:) => null()    ! allocation to fine root C (gC/m2/s)
    real(r8), pointer :: cpool_to_frootc_storage             (:) => null()    ! allocation to fine root C storage (gC/m2/s)
    real(r8), pointer :: cpool_to_livestemc                  (:) => null()    ! allocation to live stem C (gC/m2/s)
    real(r8), pointer :: cpool_to_livestemc_storage          (:) => null()    ! allocation to live stem C storage (gC/m2/s)
    real(r8), pointer :: cpool_to_deadstemc                  (:) => null()    ! allocation to dead stem C (gC/m2/s)
    real(r8), pointer :: cpool_to_deadstemc_storage          (:) => null()    ! allocation to dead stem C storage (gC/m2/s)
    real(r8), pointer :: cpool_to_livecrootc                 (:) => null()    ! allocation to live coarse root C (gC/m2/s)
    real(r8), pointer :: cpool_to_livecrootc_storage         (:) => null()    ! allocation to live coarse root C storage (gC/m2/s)
    real(r8), pointer :: cpool_to_deadcrootc                 (:) => null()    ! allocation to dead coarse root C (gC/m2/s)
    real(r8), pointer :: cpool_to_deadcrootc_storage         (:) => null()    ! allocation to dead coarse root C storage (gC/m2/s)
    real(r8), pointer :: cpool_to_gresp_storage              (:) => null()    ! allocation to growth respiration storage (gC/m2/s)
                                                                 
    ! growth respiration fluxes                                  
    real(r8), pointer :: xsmrpool_to_atm                     (:) => null()    ! excess MR pool harvest mortality (gC/m2/s)
    real(r8), pointer :: cpool_leaf_gr                       (:) => null()    ! leaf growth respiration (gC/m2/s)
    real(r8), pointer :: cpool_leaf_storage_gr               (:) => null()    ! leaf growth respiration to storage (gC/m2/s)
    real(r8), pointer :: transfer_leaf_gr                    (:) => null()    ! leaf growth respiration from storage (gC/m2/s)
    real(r8), pointer :: cpool_froot_gr                      (:) => null()    ! fine root growth respiration (gC/m2/s)
    real(r8), pointer :: cpool_froot_storage_gr              (:) => null()    ! fine root  growth respiration to storage (gC/m2/s)
    real(r8), pointer :: transfer_froot_gr                   (:) => null()    ! fine root  growth respiration from storage (gC/m2/s)
    real(r8), pointer :: cpool_livestem_gr                   (:) => null()    ! live stem growth respiration (gC/m2/s)
    real(r8), pointer :: cpool_livestem_storage_gr           (:) => null()    ! live stem growth respiration to storage (gC/m2/s)
    real(r8), pointer :: transfer_livestem_gr                (:) => null()    ! live stem growth respiration from storage (gC/m2/s)
    real(r8), pointer :: cpool_deadstem_gr                   (:) => null()    ! dead stem growth respiration (gC/m2/s)
    real(r8), pointer :: cpool_deadstem_storage_gr           (:) => null()    ! dead stem growth respiration to storage (gC/m2/s)
    real(r8), pointer :: transfer_deadstem_gr                (:) => null()    ! dead stem growth respiration from storage (gC/m2/s)
    real(r8), pointer :: cpool_livecroot_gr                  (:) => null()    ! live coarse root growth respiration (gC/m2/s)
    real(r8), pointer :: cpool_livecroot_storage_gr          (:) => null()    ! live coarse root growth respiration to storage (gC/m2/s)
    real(r8), pointer :: transfer_livecroot_gr               (:) => null()    ! live coarse root growth respiration from storage (gC/m2/s)
    real(r8), pointer :: cpool_deadcroot_gr                  (:) => null()    ! dead coarse root growth respiration (gC/m2/s)
    real(r8), pointer :: cpool_deadcroot_storage_gr          (:) => null()    ! dead coarse root growth respiration to storage (gC/m2/s)
    real(r8), pointer :: transfer_deadcroot_gr               (:) => null()    ! dead coarse root growth respiration from storage (gC/m2/s)
                                                                 
    ! growth respiration for prognostic crop model               
    real(r8), pointer :: cpool_grain_gr                      (:) => null()    ! grain growth respiration (gC/m2/s)
    real(r8), pointer :: cpool_grain_storage_gr              (:) => null()    ! grain growth respiration to storage (gC/m2/s)
    real(r8), pointer :: transfer_grain_gr                   (:) => null()    ! grain growth respiration from storage (gC/m2/s)
                                                                 
    ! annual turnover of storage to transfer pools               
    real(r8), pointer :: grainc_storage_to_xfer              (:) => null()    ! grain C shift storage to transfer for prognostic crop model (gC/m2/s)
    real(r8), pointer :: leafc_storage_to_xfer               (:) => null()    ! leaf C shift storage to transfer (gC/m2/s)
    real(r8), pointer :: frootc_storage_to_xfer              (:) => null()    ! fine root C shift storage to transfer (gC/m2/s)
    real(r8), pointer :: livestemc_storage_to_xfer           (:) => null()    ! live stem C shift storage to transfer (gC/m2/s)
    real(r8), pointer :: deadstemc_storage_to_xfer           (:) => null()    ! dead stem C shift storage to transfer (gC/m2/s)
    real(r8), pointer :: livecrootc_storage_to_xfer          (:) => null()    ! live coarse root C shift storage to transfer (gC/m2/s)
    real(r8), pointer :: deadcrootc_storage_to_xfer          (:) => null()    ! dead coarse root C shift storage to transfer (gC/m2/s)
    real(r8), pointer :: gresp_storage_to_xfer               (:) => null()    ! growth respiration shift storage to transfer (gC/m2/s)
                                                                 
    ! turnover of livewood to deadwood                           
    real(r8), pointer :: livestemc_to_deadstemc              (:) => null()    ! live stem C turnover (gC/m2/s)
    real(r8), pointer :: livecrootc_to_deadcrootc            (:) => null()    ! live coarse root C turnover (gC/m2/s)
                                                                 
    ! summary (diagnostic) flux variables, not involved in mass balance
    real(r8), pointer :: gpp                                 (:) => null()    ! (gC/m2/s) gross primary production 
    real(r8), pointer :: gpp_before_downreg                  (:) => null()    ! (gC/m2/s) gross primary production before down regulation
    real(r8), pointer :: mr                                  (:) => null()    ! (gC/m2/s) maintenance respiration
    real(r8), pointer :: current_gr                          (:) => null()    ! (gC/m2/s) growth resp for new growth displayed in this timestep
    real(r8), pointer :: transfer_gr                         (:) => null()    ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
    real(r8), pointer :: storage_gr                          (:) => null()    ! (gC/m2/s) growth resp for growth sent to storage for later display
    real(r8), pointer :: gr                                  (:) => null()    ! (gC/m2/s) total growth respiration
    real(r8), pointer :: ar                                  (:) => null()    ! (gC/m2/s) autotrophic respiration (MR + GR)
    real(r8), pointer :: rr                                  (:) => null()    ! (gC/m2/s) root respiration (fine root MR + total root GR)
    real(r8), pointer :: npp                                 (:) => null()    ! (gC/m2/s) net primary production
    real(r8), pointer :: agnpp                               (:) => null()    ! (gC/m2/s) aboveground NPP
    real(r8), pointer :: bgnpp                               (:) => null()    ! (gC/m2/s) belowground NPP
    real(r8), pointer :: litfall                             (:) => null()    ! (gC/m2/s) litterfall (leaves and fine roots)
    real(r8), pointer :: vegfire                             (:) => null()    ! (gC/m2/s) patch-level fire loss (obsolete, mark for removal)
    real(r8), pointer :: wood_harvestc                       (:) => null()    ! (gC/m2/s) patch-level wood harvest (to product pools)
    real(r8), pointer :: cinputs                             (:) => null()    ! (gC/m2/s) patch-level carbon inputs (for balance checking)
    real(r8), pointer :: coutputs                            (:) => null()    ! (gC/m2/s) patch-level carbon outputs (for balance checking)
                                                                 
    real(r8), pointer :: plant_calloc                        (:) => null()    ! total allocated C flux (gC/m2/s)
    real(r8), pointer :: excess_cflux                        (:) => null()    ! C flux not allocated due to downregulation (gC/m2/s)
    real(r8), pointer :: prev_leafc_to_litter                (:) => null()    ! previous timestep leaf C litterfall flux (gC/m2/s)
    real(r8), pointer :: prev_frootc_to_litter               (:) => null()    ! previous timestep froot C litterfall flux (gC/m2/s)
    real(r8), pointer :: availc                              (:) => null()    ! C flux available for allocation (gC/m2/s)
    real(r8), pointer :: xsmrpool_recover                    (:) => null()    ! C flux assigned to recovery of negative cpool (gC/m2/s)
    real(r8), pointer :: xsmrpool_c13ratio                   (:) => null()    ! C13/C(12+13) ratio for xsmrpool (proportion)
    real(r8), pointer :: xsmrpool_turnover                   (:) => null()    ! xsmrpool flux to atmosphere due to turnover
                                                                 
    ! CN: CLAMP summary (diagnostic) variables, not involved in mass balance
    real(r8), pointer :: frootc_alloc                        (:) => null()    ! (gC/m2/s) patch-level fine root C alloc
    real(r8), pointer :: frootc_loss                         (:) => null()    ! (gC/m2/s) patch-level fine root C loss
    real(r8), pointer :: leafc_alloc                         (:) => null()    ! (gC/m2/s) patch-level leaf C alloc
    real(r8), pointer :: leafc_loss                          (:) => null()    ! (gC/m2/s) patch-level leaf C loss
    real(r8), pointer :: woodc_alloc                         (:) => null()    ! (gC/m2/s) patch-level wood C alloc
    real(r8), pointer :: woodc_loss                          (:) => null()    ! (gC/m2/s) patch-level wood C loss

    ! fire code                                                 
    real(r8), pointer :: fire_closs                          (:) => null()    ! (gC/m2/s) total patch-level fire C loss 
                                                                 
    ! crop fluxes
    real(r8), pointer :: crop_seedc_to_leaf                  (:) => null()    ! (gC/m2/s) seed source to leaf, for crops

    ! CN dynamic landcover fluxes
    real(r8), pointer :: dwt_seedc_to_leaf                   (:) => null()    ! (gC/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_seedc_to_deadstem               (:) => null()    ! (gC/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_conv_cflux                      (:) => null()    ! (gC/m2/s) conversion C flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_prod10c_gain                    (:) => null()    ! (gC/m2/s) addition to 10-yr wood product pool; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_prod100c_gain                   (:) => null()    ! (gC/m2/s) addition to 100-yr wood product pool; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_crop_productc_gain              (:) => null()    ! (gC/m2/s) addition to crop product pools from landcover change; although this is a patch-level flux, it is expressed per unit GRIDCELL area
                                                             
    ! Debug, Temporary, and annual sums                              
    real(r8), pointer :: tempsum_npp                         (:) => null()    ! patch temporary annual sum of NPP (gC/m2/yr)
    real(r8), pointer :: annsum_npp                          (:) => null()    ! patch annual sum of NPP (gC/m2/yr)
    real(r8), pointer :: annavg_agnpp                        (:) => null()    ! (gC/m2/s) annual average aboveground NPP
    real(r8), pointer :: annavg_bgnpp                        (:) => null()    ! (gC/m2/s) annual average belowground NPP
    real(r8), pointer :: tempavg_agnpp                       (:) => null()    ! (gC/m2/s) temp. average aboveground NPP
    real(r8), pointer :: tempavg_bgnpp                       (:) => null()    ! (gC/m2/s) temp. average belowground NPP
    real(r8), pointer :: allocation_leaf 		                (:) => null()    ! check allocation to leaf for dynamic allocation scheme
    real(r8), pointer :: allocation_stem 		                (:) => null()    ! check allocation to stem for dynamic allocation scheme
    real(r8), pointer :: allocation_froot 		             (:) => null()    ! check allocation to fine root for dynamic allocation scheme
                                                                 
    ! For comparison with RAINFOR wood productivity data         
    real(r8), pointer :: agwdnpp                             (:) => null()    !(gC/m2/s) aboveground NPP

  contains
    procedure, public :: Init       => veg_cf_init
  end type vegetation_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_flux
    ! gap mortality fluxes
    real(r8), pointer :: m_leafn_to_litter                   (:)   => null()  ! leaf N mortality (gN/m2/s)
    real(r8), pointer :: m_frootn_to_litter                  (:)   => null()  ! fine root N mortality (gN/m2/s)
    real(r8), pointer :: m_leafn_storage_to_litter           (:)   => null()  ! leaf N storage mortality (gN/m2/s)
    real(r8), pointer :: m_frootn_storage_to_litter          (:)   => null()  ! fine root N storage mortality (gN/m2/s)
    real(r8), pointer :: m_livestemn_storage_to_litter       (:)   => null()  ! live stem N storage mortality (gN/m2/s)
    real(r8), pointer :: m_deadstemn_storage_to_litter       (:)   => null()  ! dead stem N storage mortality (gN/m2/s)
    real(r8), pointer :: m_livecrootn_storage_to_litter      (:)   => null()  ! live coarse root N storage mortality (gN/m2/s)
    real(r8), pointer :: m_deadcrootn_storage_to_litter      (:)   => null()  ! dead coarse root N storage mortality (gN/m2/s)
    real(r8), pointer :: m_leafn_xfer_to_litter              (:)   => null()  ! leaf N transfer mortality (gN/m2/s)
    real(r8), pointer :: m_frootn_xfer_to_litter             (:)   => null()  ! fine root N transfer mortality (gN/m2/s)
    real(r8), pointer :: m_livestemn_xfer_to_litter          (:)   => null()  ! live stem N transfer mortality (gN/m2/s)
    real(r8), pointer :: m_deadstemn_xfer_to_litter          (:)   => null()  ! dead stem N transfer mortality (gN/m2/s)
    real(r8), pointer :: m_livecrootn_xfer_to_litter         (:)   => null()  ! live coarse root N transfer mortality (gN/m2/s)
    real(r8), pointer :: m_deadcrootn_xfer_to_litter         (:)   => null()  ! dead coarse root N transfer mortality (gN/m2/s)
    real(r8), pointer :: m_livestemn_to_litter               (:)   => null()  ! live stem N mortality (gN/m2/s)
    real(r8), pointer :: m_deadstemn_to_litter               (:)   => null()  ! dead stem N mortality (gN/m2/s)
    real(r8), pointer :: m_livecrootn_to_litter              (:)   => null()  ! live coarse root N mortality (gN/m2/s)
    real(r8), pointer :: m_deadcrootn_to_litter              (:)   => null()  ! dead coarse root N mortality (gN/m2/s)
    real(r8), pointer :: m_retransn_to_litter                (:)   => null()  ! retranslocated N pool mortality (gN/m2/s)
    real(r8), pointer :: m_npool_to_litter                   (:)   => null()  ! npool mortality (gN/m2/s)
    ! harvest fluxes                                                   
    real(r8), pointer :: hrv_leafn_to_litter                 (:)   => null()  ! leaf N harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_frootn_to_litter                (:)   => null()  ! fine root N harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_leafn_storage_to_litter         (:)   => null()  ! leaf N storage harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_frootn_storage_to_litter        (:)   => null()  ! fine root N storage harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_livestemn_storage_to_litter     (:)   => null()  ! live stem N storage harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_deadstemn_storage_to_litter     (:)   => null()  ! dead stem N storage harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_livecrootn_storage_to_litter    (:)   => null()  ! live coarse root N storage harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_deadcrootn_storage_to_litter    (:)   => null()  ! dead coarse root N storage harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_leafn_xfer_to_litter            (:)   => null()  ! leaf N transfer harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_frootn_xfer_to_litter           (:)   => null()  ! fine root N transfer harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_livestemn_xfer_to_litter        (:)   => null()  ! live stem N transfer harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_deadstemn_xfer_to_litter        (:)   => null()  ! dead stem N transfer harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_livecrootn_xfer_to_litter       (:)   => null()  ! live coarse root N transfer harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_deadcrootn_xfer_to_litter       (:)   => null()  ! dead coarse root N transfer harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_livestemn_to_litter             (:)   => null()  ! live stem N harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_deadstemn_to_prod10n            (:)   => null()  ! dead stem N harvest to 10-year product pool (gN/m2/s)
    real(r8), pointer :: hrv_deadstemn_to_prod100n           (:)   => null()  ! dead stem N harvest to 100-year product pool (gN/m2/s)
    real(r8), pointer :: hrv_livecrootn_to_litter            (:)   => null()  ! live coarse root N harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_deadcrootn_to_litter            (:)   => null()  ! dead coarse root N harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_retransn_to_litter              (:)   => null()  ! retranslocated N pool harvest mortality (gN/m2/s)
    real(r8), pointer :: hrv_npool_to_litter                 (:)   => null()  ! npool to harvest mortalty (gN/m2/s)
    ! crop harvest
    real(r8), pointer :: hrv_leafn_to_prod1n                 (:)   => null()  ! crop leaf N harvested (gN/m2/s)
    real(r8), pointer :: hrv_livestemn_to_prod1n             (:)   => null()  ! crop stem N harvested (gN/m2/s)
    real(r8), pointer :: hrv_grainn_to_prod1n                (:)   => null()  ! crop grain N harvested (gN/m2/s)
    real(r8), pointer :: hrv_cropn_to_prod1n                 (:)   => null()  ! total amount of crop N harvested (gN/m2/s)
    ! fire N fluxes 
    real(r8), pointer :: m_leafn_to_fire                     (:)   => null()  ! (gN/m2/s) fire N emissions from leafn            
    real(r8), pointer :: m_leafn_storage_to_fire             (:)   => null()  ! (gN/m2/s) fire N emissions from leafn_storage            
    real(r8), pointer :: m_leafn_xfer_to_fire                (:)   => null()  ! (gN/m2/s) fire N emissions from leafn_xfer     
    real(r8), pointer :: m_livestemn_to_fire                 (:)   => null()  ! (gN/m2/s) fire N emissions from livestemn 
    real(r8), pointer :: m_livestemn_storage_to_fire         (:)   => null()  ! (gN/m2/s) fire N emissions from livestemn_storage      
    real(r8), pointer :: m_livestemn_xfer_to_fire            (:)   => null()  ! (gN/m2/s) fire N emissions from livestemn_xfer
    real(r8), pointer :: m_deadstemn_to_fire                 (:)   => null()  ! (gN/m2/s) fire N emissions from deadstemn
    real(r8), pointer :: m_deadstemn_storage_to_fire         (:)   => null()  ! (gN/m2/s) fire N emissions from deadstemn_storage         
    real(r8), pointer :: m_deadstemn_xfer_to_fire            (:)   => null()  ! (gN/m2/s) fire N emissions from deadstemn_xfer
    real(r8), pointer :: m_frootn_to_fire                    (:)   => null()  ! (gN/m2/s) fire N emissions from frootn
    real(r8), pointer :: m_frootn_storage_to_fire            (:)   => null()  ! (gN/m2/s) fire N emissions from frootn_storage
    real(r8), pointer :: m_frootn_xfer_to_fire               (:)   => null()  ! (gN/m2/s) fire N emissions from frootn_xfer
    real(r8), pointer :: m_livecrootn_to_fire                (:)   => null()  ! (gN/m2/s) fire N emissions from m_livecrootn_to_fire
    real(r8), pointer :: m_livecrootn_storage_to_fire        (:)   => null()  ! (gN/m2/s) fire N emissions from livecrootn_storage     
    real(r8), pointer :: m_livecrootn_xfer_to_fire           (:)   => null()  ! (gN/m2/s) fire N emissions from livecrootn_xfer
    real(r8), pointer :: m_deadcrootn_to_fire                (:)   => null()  ! (gN/m2/s) fire N emissions from deadcrootn
    real(r8), pointer :: m_deadcrootn_storage_to_fire        (:)   => null()  ! (gN/m2/s) fire N emissions from deadcrootn_storage  
    real(r8), pointer :: m_deadcrootn_xfer_to_fire           (:)   => null()  ! (gN/m2/s) fire N emissions from deadcrootn_xfer
    real(r8), pointer :: m_retransn_to_fire                  (:)   => null()  ! (gN/m2/s) fire N emissions from retransn
    real(r8), pointer :: m_npool_to_fire                     (:)   => null()  ! (gN/m2/s) fire N emissions from npool
    real(r8), pointer :: m_leafn_to_litter_fire              (:)   => null()  ! (gN/m2/s) from leafn to litter N  due to fire               
    real(r8), pointer :: m_leafn_storage_to_litter_fire      (:)   => null()  ! (gN/m2/s) from leafn_storage to litter N  due to fire                              
    real(r8), pointer :: m_leafn_xfer_to_litter_fire         (:)   => null()  ! (gN/m2/s) from leafn_xfer to litter N  due to fire                              
    real(r8), pointer :: m_livestemn_to_litter_fire          (:)   => null()  ! (gN/m2/s) from livestemn to litter N  due to fire                              
    real(r8), pointer :: m_livestemn_storage_to_litter_fire  (:)   => null()  ! (gN/m2/s) from livestemn_storage to litter N  due to fire                                     
    real(r8), pointer :: m_livestemn_xfer_to_litter_fire     (:)   => null()  ! (gN/m2/s) from livestemn_xfer to litter N  due to fire                                     
    real(r8), pointer :: m_livestemn_to_deadstemn_fire       (:)   => null()  ! (gN/m2/s) from livestemn to deadstemn N  due to fire                                     
    real(r8), pointer :: m_deadstemn_to_litter_fire          (:)   => null()  ! (gN/m2/s) from deadstemn to litter N  due to fire                                     
    real(r8), pointer :: m_deadstemn_storage_to_litter_fire  (:)   => null()  ! (gN/m2/s) from deadstemn_storage to litter N  due to fire                                               
    real(r8), pointer :: m_deadstemn_xfer_to_litter_fire     (:)   => null()  ! (gN/m2/s) from deadstemn_xfer to litter N  due to fire                                               
    real(r8), pointer :: m_frootn_to_litter_fire             (:)   => null()  ! (gN/m2/s) from frootn to litter N  due to fire                                               
    real(r8), pointer :: m_frootn_storage_to_litter_fire     (:)   => null()  ! (gN/m2/s) from frootn_storage to litter N  due to fire                                               
    real(r8), pointer :: m_frootn_xfer_to_litter_fire        (:)   => null()  ! (gN/m2/s) from frootn_xfer to litter N  due to fire                                               
    real(r8), pointer :: m_livecrootn_to_litter_fire         (:)   => null()  ! (gN/m2/s) from livecrootn to litter N  due to fire                                               
    real(r8), pointer :: m_livecrootn_storage_to_litter_fire (:)   => null()  ! (gN/m2/s) from livecrootn_storage to litter N  due to fire                                                     
    real(r8), pointer :: m_livecrootn_xfer_to_litter_fire    (:)   => null()  ! (gN/m2/s) from livecrootn_xfer to litter N  due to fire                                                     
    real(r8), pointer :: m_livecrootn_to_deadcrootn_fire     (:)   => null()  ! (gN/m2/s) from livecrootn_xfer to deadcrootn due to fire                                                     
    real(r8), pointer :: m_deadcrootn_to_litter_fire         (:)   => null()  ! (gN/m2/s) from deadcrootn to deadcrootn due to fire                                                       
    real(r8), pointer :: m_deadcrootn_storage_to_litter_fire (:)   => null()  ! (gN/m2/s) from deadcrootn_storage to deadcrootn due to fire                                                        
    real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire    (:)   => null()  ! (gN/m2/s) from deadcrootn_xfer to deadcrootn due to fire                                                         
    real(r8), pointer :: m_retransn_to_litter_fire           (:)   => null()  ! (gN/m2/s) from retransn to deadcrootn due to fire                                                         
    real(r8), pointer :: m_npool_to_litter_fire              (:)   => null()  ! (gN/m2/s) from npool to litter due to fire
    real(r8), pointer :: fire_nloss                          (:)   => null()  ! total pft-level fire N loss (gN/m2/s) 
    ! phenology fluxes from transfer pool
    real(r8), pointer :: grainn_xfer_to_grainn               (:)   => null()  ! grain N growth from storage for prognostic crop model (gN/m2/s)
    real(r8), pointer :: leafn_xfer_to_leafn                 (:)   => null()  ! leaf N growth from storage (gN/m2/s)
    real(r8), pointer :: frootn_xfer_to_frootn               (:)   => null()  ! fine root N growth from storage (gN/m2/s)
    real(r8), pointer :: livestemn_xfer_to_livestemn         (:)   => null()  ! live stem N growth from storage (gN/m2/s)
    real(r8), pointer :: deadstemn_xfer_to_deadstemn         (:)   => null()  ! dead stem N growth from storage (gN/m2/s)
    real(r8), pointer :: livecrootn_xfer_to_livecrootn       (:)   => null()  ! live coarse root N growth from storage (gN/m2/s)
    real(r8), pointer :: deadcrootn_xfer_to_deadcrootn       (:)   => null()  ! dead coarse root N growth from storage (gN/m2/s)
    ! litterfall fluxes
    real(r8), pointer :: livestemn_to_litter                 (:)   => null()  ! livestem N to litter (gN/m2/s)
    real(r8), pointer :: grainn_to_food                      (:)   => null()  ! grain N to food for prognostic crop (gN/m2/s)
    real(r8), pointer :: leafn_to_litter                     (:)   => null()  ! leaf N litterfall (gN/m2/s)
    real(r8), pointer :: leafn_to_retransn                   (:)   => null()  ! leaf N to retranslocated N pool (gN/m2/s)
    real(r8), pointer :: frootn_to_retransn                  (:)   => null()  ! fine root N to retranslocated N pool (gN/m2/s)
    real(r8), pointer :: frootn_to_litter                    (:)   => null()  ! fine root N litterfall (gN/m2/s)
    ! allocation fluxes
    real(r8), pointer :: retransn_to_npool                   (:)   => null()  ! deployment of retranslocated N (gN/m2/s)       
    real(r8), pointer :: sminn_to_npool                      (:)   => null()  ! deployment of soil mineral N uptake (gN/m2/s)
    real(r8), pointer :: npool_to_grainn                     (:)   => null()  ! allocation to grain N for prognostic crop (gN/m2/s)
    real(r8), pointer :: npool_to_grainn_storage             (:)   => null()  ! allocation to grain N storage for prognostic crop (gN/m2/s)
    real(r8), pointer :: npool_to_leafn                      (:)   => null()  ! allocation to leaf N (gN/m2/s)
    real(r8), pointer :: npool_to_leafn_storage              (:)   => null()  ! allocation to leaf N storage (gN/m2/s)
    real(r8), pointer :: npool_to_frootn                     (:)   => null()  ! allocation to fine root N (gN/m2/s)
    real(r8), pointer :: npool_to_frootn_storage             (:)   => null()  ! allocation to fine root N storage (gN/m2/s)
    real(r8), pointer :: npool_to_livestemn                  (:)   => null()  ! allocation to live stem N (gN/m2/s)
    real(r8), pointer :: npool_to_livestemn_storage          (:)   => null()  ! allocation to live stem N storage (gN/m2/s)
    real(r8), pointer :: npool_to_deadstemn                  (:)   => null()  ! allocation to dead stem N (gN/m2/s)
    real(r8), pointer :: npool_to_deadstemn_storage          (:)   => null()  ! allocation to dead stem N storage (gN/m2/s)
    real(r8), pointer :: npool_to_livecrootn                 (:)   => null()  ! allocation to live coarse root N (gN/m2/s)
    real(r8), pointer :: npool_to_livecrootn_storage         (:)   => null()  ! allocation to live coarse root N storage (gN/m2/s)
    real(r8), pointer :: npool_to_deadcrootn                 (:)   => null()  ! allocation to dead coarse root N (gN/m2/s)
    real(r8), pointer :: npool_to_deadcrootn_storage         (:)   => null()  ! allocation to dead coarse root N storage (gN/m2/s)
    ! annual turnover of storage to transfer pools           
    real(r8), pointer :: grainn_storage_to_xfer              (:)   => null()  ! grain N shift storage to transfer for prognostic crop (gN/m2/s)
    real(r8), pointer :: leafn_storage_to_xfer               (:)   => null()  ! leaf N shift storage to transfer (gN/m2/s)
    real(r8), pointer :: frootn_storage_to_xfer              (:)   => null()  ! fine root N shift storage to transfer (gN/m2/s)
    real(r8), pointer :: livestemn_storage_to_xfer           (:)   => null()  ! live stem N shift storage to transfer (gN/m2/s)
    real(r8), pointer :: deadstemn_storage_to_xfer           (:)   => null()  ! dead stem N shift storage to transfer (gN/m2/s)
    real(r8), pointer :: livecrootn_storage_to_xfer          (:)   => null()  ! live coarse root N shift storage to transfer (gN/m2/s)
    real(r8), pointer :: deadcrootn_storage_to_xfer          (:)   => null()  ! dead coarse root N shift storage to transfer (gN/m2/s)
    real(r8), pointer :: fert                                (:)   => null()  ! applied fertilizer (gN/m2/s)
    real(r8), pointer :: fert_counter                        (:)   => null()  ! >0 fertilize; <=0 not
    real(r8), pointer :: soyfixn                             (:)   => null()  ! soybean fixed N (gN/m2/s)
    ! turnover of livewood to deadwood, with retranslocation 
    real(r8), pointer :: livestemn_to_deadstemn              (:)   => null()  ! live stem N turnover (gN/m2/s)
    real(r8), pointer :: livestemn_to_retransn               (:)   => null()  ! live stem N to retranslocated N pool (gN/m2/s)
    real(r8), pointer :: livecrootn_to_deadcrootn            (:)   => null()  ! live coarse root N turnover (gN/m2/s)
    real(r8), pointer :: livecrootn_to_retransn              (:)   => null()  ! live coarse root N to retranslocated N pool (gN/m2/s)
    ! summary (diagnostic) flux variables, not involved in mass balance
    real(r8), pointer :: ndeploy                             (:)   => null()  ! total N deployed to growth and storage (gN/m2/s)
    real(r8), pointer :: ninputs                             (:)   => null()  ! total N inputs to pft-level (gN/m2/s)
    real(r8), pointer :: noutputs                            (:)   => null()  ! total N outputs from pft-level (gN/m2/s)
    real(r8), pointer :: wood_harvestn                       (:)   => null()  ! total N losses to wood product pools (gN/m2/s)
    ! deposition fluxes
    real(r8), pointer :: nfix_to_plantn                      (:)   => null()  ! nitrogen fixation goes to plant
    ! dynamic landcover fluxes
    real(r8), pointer :: crop_seedn_to_leaf                  (:)   => null()  ! (gN/m2/s) seed source to leaf, for crops
    real(r8), pointer :: dwt_seedn_to_leaf                   (:)   => null()  ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_seedn_to_deadstem               (:)   => null()  ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_conv_nflux                      (:)   => null()  ! (gN/m2/s) conversion N flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_prod10n_gain                    (:)   => null()  ! (gN/m2/s) addition to 10-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_prod100n_gain                   (:)   => null()  ! (gN/m2/s) addition to 100-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_crop_productn_gain              (:)   => null()  ! (gN/m2/s) addition to crop product pool from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_seedn_to_npool                  (:)   => null()  ! (gN/m2/s) seed source to PFT level
    ! Misc
    real(r8), pointer :: plant_ndemand                       (:)   => null()  ! N flux required to support initial GPP (gN/m2/s)
    real(r8), pointer :: avail_retransn                      (:)   => null()  ! N flux available from retranslocation pool (gN/m2/s)
    real(r8), pointer :: plant_nalloc                        (:)   => null()  ! total allocated N flux (gN/m2/s)
    real(r8), pointer :: smin_no3_to_plant                   (:)   => null()  ! pft level plant uptake of soil NO3 (gN/m2/s) BGC mode
    real(r8), pointer :: smin_nh4_to_plant                   (:)   => null()  ! pft level plant uptake of soil Nh4 (gN/m2/s) BGC mode
    real(r8), pointer :: sminn_to_plant                      (:)   => null()  ! pft level plant uptake of soil N (gN/m2/s) CN mode
    real(r8), pointer :: plant_nh4demand_vr                  (:,:) => null()  ! pft-level plant NH4 demand BGC mode
    real(r8), pointer :: plant_no3demand_vr                  (:,:) => null()  ! pft-level plant NO3 demand BGC mode
    real(r8), pointer :: plant_ndemand_vr                    (:,:) => null()  ! pft-level plant N demand CN mode
    real(r8), pointer :: prev_leafn_to_litter                (:)   => null()  ! previous timestep leaf N litterfall flux (gN/m2/s)
    real(r8), pointer :: prev_frootn_to_litter               (:)   => null()  ! previous timestep froot N litterfall flux (gN/m2/s)
    real(r8), pointer :: supplement_to_plantn                (:)   => null()  ! supplementary N flux for plant
    real(r8), pointer :: gap_nloss_litter                    (:)   => null()  ! total nloss from veg to litter pool due to gap mortality
    real(r8), pointer :: fire_nloss_litter                   (:)   => null()  ! total nloss from veg to litter pool due to fire
    real(r8), pointer :: hrv_nloss_litter                    (:)   => null()  ! total nloss from veg to litter pool due to harvest mortality
    real(r8), pointer :: sen_nloss_litter                    (:)   => null()  ! total nloss from veg to litter pool due to senescence
  contains
    procedure, public :: Init      => veg_nf_init
  end type vegetation_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_flux
    real(r8), pointer :: m_leafp_to_litter                   (:)     ! leaf P mortality (gP/m2/s)
    real(r8), pointer :: m_frootp_to_litter                  (:)     ! fine root P mortality (gP/m2/s)
    real(r8), pointer :: m_leafp_storage_to_litter           (:)     ! leaf P storage mortality (gP/m2/s)
    real(r8), pointer :: m_frootp_storage_to_litter          (:)     ! fine root P storage mortality (gP/m2/s)
    real(r8), pointer :: m_livestemp_storage_to_litter       (:)     ! live stem P storage mortality (gP/m2/s)
    real(r8), pointer :: m_deadstemp_storage_to_litter       (:)     ! dead stem P storage mortality (gP/m2/s)
    real(r8), pointer :: m_livecrootp_storage_to_litter      (:)     ! live coarse root P storage mortality (gP/m2/s)
    real(r8), pointer :: m_deadcrootp_storage_to_litter      (:)     ! dead coarse root P storage mortality (gP/m2/s)
    real(r8), pointer :: m_leafp_xfer_to_litter              (:)     ! leaf P transfer mortality (gP/m2/s)
    real(r8), pointer :: m_frootp_xfer_to_litter             (:)     ! fine root P transfer mortality (gP/m2/s)
    real(r8), pointer :: m_livestemp_xfer_to_litter          (:)     ! live stem P transfer mortality (gP/m2/s)
    real(r8), pointer :: m_deadstemp_xfer_to_litter          (:)     ! dead stem P transfer mortality (gP/m2/s)
    real(r8), pointer :: m_livecrootp_xfer_to_litter         (:)     ! live coarse root P transfer mortality (gP/m2/s)
    real(r8), pointer :: m_deadcrootp_xfer_to_litter         (:)     ! dead coarse root P transfer mortality (gP/m2/s)
    real(r8), pointer :: m_livestemp_to_litter               (:)     ! live stem P mortality (gP/m2/s)
    real(r8), pointer :: m_deadstemp_to_litter               (:)     ! dead stem P mortality (gP/m2/s)
    real(r8), pointer :: m_livecrootp_to_litter              (:)     ! live coarse root P mortality (gP/m2/s)
    real(r8), pointer :: m_deadcrootp_to_litter              (:)     ! dead coarse root P mortality (gP/m2/s)
    real(r8), pointer :: m_retransp_to_litter                (:)     ! retranslocated P pool mortality (gP/m2/s)
    real(r8), pointer :: m_ppool_to_litter                   (:)     ! storage P pool mortality (gP/m2/s)
    real(r8), pointer :: hrv_leafp_to_litter                 (:)     ! leaf P harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_frootp_to_litter                (:)     ! fine root P harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_leafp_storage_to_litter         (:)     ! leaf P storage harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_frootp_storage_to_litter        (:)     ! fine root P storage harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_livestemp_storage_to_litter     (:)     ! live stem P storage harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_deadstemp_storage_to_litter     (:)     ! dead stem P storage harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_livecrootp_storage_to_litter    (:)     ! live coarse root P storage harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_deadcrootp_storage_to_litter    (:)     ! dead coarse root P storage harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_leafp_xfer_to_litter            (:)     ! leaf P transfer harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_frootp_xfer_to_litter           (:)     ! fine root P transfer harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_livestemp_xfer_to_litter        (:)     ! live stem P transfer harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_deadstemp_xfer_to_litter        (:)     ! dead stem P transfer harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_livecrootp_xfer_to_litter       (:)     ! live coarse root P transfer harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_deadcrootp_xfer_to_litter       (:)     ! dead coarse root P transfer harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_livestemp_to_litter             (:)     ! live stem P harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_deadstemp_to_prod10p            (:)     ! dead stem P harvest to 10-year product pool (gP/m2/s)
    real(r8), pointer :: hrv_deadstemp_to_prod100p           (:)     ! dead stem P harvest to 100-year product pool (gP/m2/s)
    real(r8), pointer :: hrv_livecrootp_to_litter            (:)     ! live coarse root P harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_deadcrootp_to_litter            (:)     ! dead coarse root P harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_retransp_to_litter              (:)     ! retranslocated P pool harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_ppool_to_litter                 (:)     ! retranslocated P pool harvest mortality (gP/m2/s)
    real(r8), pointer :: hrv_leafp_to_prod1p                 (:)     ! crop leafp harvested (gP/m2/s)
    real(r8), pointer :: hrv_livestemp_to_prod1p             (:)     ! crop stemp harvested (gP/m2/s)
    real(r8), pointer :: hrv_grainp_to_prod1p                (:)     ! crop grain harvested (gP/m2/s)
    real(r8), pointer :: hrv_cropp_to_prod1p                 (:)     ! total amount of crop P harvested (gP/m2/s)
    real(r8), pointer :: m_leafp_to_fire                     (:)     ! (gP/m2/s) fire P emissions from leafp 
    real(r8), pointer :: m_leafp_storage_to_fire             (:)     ! (gP/m2/s) fire P emissions from leafp_storage            
    real(r8), pointer :: m_leafp_xfer_to_fire                (:)     ! (gP/m2/s) fire P emissions from leafp_xfer     
    real(r8), pointer :: m_livestemp_to_fire                 (:)     ! (gP/m2/s) fire P emissions from livestemp 
    real(r8), pointer :: m_livestemp_storage_to_fire         (:)     ! (gP/m2/s) fire P emissions from livestemp_storage      
    real(r8), pointer :: m_livestemp_xfer_to_fire            (:)     ! (gP/m2/s) fire P emissions from livestemp_xfer
    real(r8), pointer :: m_deadstemp_to_fire                 (:)     ! (gP/m2/s) fire P emissions from deadstemp
    real(r8), pointer :: m_deadstemp_storage_to_fire         (:)     ! (gP/m2/s) fire P emissions from deadstemp_storage         
    real(r8), pointer :: m_deadstemp_xfer_to_fire            (:)     ! (gP/m2/s) fire P emissions from deadstemp_xfer
    real(r8), pointer :: m_frootp_to_fire                    (:)     ! (gP/m2/s) fire P emissions from frootp
    real(r8), pointer :: m_frootp_storage_to_fire            (:)     ! (gP/m2/s) fire P emissions from frootp_storage
    real(r8), pointer :: m_frootp_xfer_to_fire               (:)     ! (gP/m2/s) fire P emissions from frootp_xfer
    real(r8), pointer :: m_livecrootp_to_fire                (:)     ! (gP/m2/s) fire P emissions from m_livecrootp_to_fire
    real(r8), pointer :: m_livecrootp_storage_to_fire        (:)     ! (gP/m2/s) fire P emissions from livecrootp_storage     
    real(r8), pointer :: m_livecrootp_xfer_to_fire           (:)     ! (gP/m2/s) fire P emissions from livecrootp_xfer
    real(r8), pointer :: m_deadcrootp_to_fire                (:)     ! (gP/m2/s) fire P emissions from deadcrootp
    real(r8), pointer :: m_deadcrootp_storage_to_fire        (:)     ! (gP/m2/s) fire P emissions from deadcrootp_storage  
    real(r8), pointer :: m_deadcrootp_xfer_to_fire           (:)     ! (gP/m2/s) fire P emissions from deadcrootp_xfer
    real(r8), pointer :: m_retransp_to_fire                  (:)     ! (gP/m2/s) fire P emissions from retransp
    real(r8), pointer :: m_ppool_to_fire                     (:)     ! (gP/m2/s) fire P emissions from ppool
    real(r8), pointer :: m_leafp_to_litter_fire              (:)     ! (gP/m2/s) from leafp to litter P  due to fire               
    real(r8), pointer :: m_leafp_storage_to_litter_fire      (:)     ! (gP/m2/s) from leafp_storage to litter P  due to fire                              
    real(r8), pointer :: m_leafp_xfer_to_litter_fire         (:)     ! (gP/m2/s) from leafp_xfer to litter P  due to fire                              
    real(r8), pointer :: m_livestemp_to_litter_fire          (:)     ! (gP/m2/s) from livestemp to litter P  due to fire                              
    real(r8), pointer :: m_livestemp_storage_to_litter_fire  (:)     ! (gP/m2/s) from livestemp_storage to litter P  due to fire                                     
    real(r8), pointer :: m_livestemp_xfer_to_litter_fire     (:)     ! (gP/m2/s) from livestemp_xfer to litter P  due to fire                                     
    real(r8), pointer :: m_livestemp_to_deadstemp_fire       (:)     ! (gP/m2/s) from livestemp to deadstemp P  due to fire                                     
    real(r8), pointer :: m_deadstemp_to_litter_fire          (:)     ! (gP/m2/s) from deadstemp to litter P  due to fire                                     
    real(r8), pointer :: m_deadstemp_storage_to_litter_fire  (:)     ! (gP/m2/s) from deadstemp_storage to litter P  due to fire                                               
    real(r8), pointer :: m_deadstemp_xfer_to_litter_fire     (:)     ! (gP/m2/s) from deadstemp_xfer to litter P  due to fire                                               
    real(r8), pointer :: m_frootp_to_litter_fire             (:)     ! (gP/m2/s) from frootp to litter P  due to fire                                               
    real(r8), pointer :: m_frootp_storage_to_litter_fire     (:)     ! (gP/m2/s) from frootp_storage to litter P  due to fire                                               
    real(r8), pointer :: m_frootp_xfer_to_litter_fire        (:)     ! (gP/m2/s) from frootp_xfer to litter P  due to fire                                               
    real(r8), pointer :: m_livecrootp_to_litter_fire         (:)     ! (gP/m2/s) from livecrootp to litter P  due to fire                                               
    real(r8), pointer :: m_livecrootp_storage_to_litter_fire (:)     ! (gP/m2/s) from livecrootp_storage to litter P  due to fire                                                     
    real(r8), pointer :: m_livecrootp_xfer_to_litter_fire    (:)     ! (gP/m2/s) from livecrootp_xfer to litter P  due to fire                                                     
    real(r8), pointer :: m_livecrootp_to_deadcrootp_fire     (:)     ! (gP/m2/s) from livecrootp_xfer to deadcrootp due to fire                                                     
    real(r8), pointer :: m_deadcrootp_to_litter_fire         (:)     ! (gP/m2/s) from deadcrootp to deadcrootp due to fire                                                       
    real(r8), pointer :: m_deadcrootp_storage_to_litter_fire (:)     ! (gP/m2/s) from deadcrootp_storage to deadcrootp due to fire                                                        
    real(r8), pointer :: m_deadcrootp_xfer_to_litter_fire    (:)     ! (gP/m2/s) from deadcrootp_xfer to deadcrootp due to fire 
    real(r8), pointer :: m_retransp_to_litter_fire           (:)     ! (gP/m2/s) from retransp to deadcrootp due to fire                                                               
    real(r8), pointer :: m_ppool_to_litter_fire              (:)     ! (gP/m2/s) from ppool to deadcrootp due to fire                                                         
    real(r8), pointer :: fire_ploss                          (:)     ! total pft-level fire P loss (gP/m2/s) 
    real(r8), pointer :: grainp_xfer_to_grainp               (:)     ! grain P growth from storage for prognostic crop model (gP/m2/s)
    real(r8), pointer :: leafp_xfer_to_leafp                 (:)     ! leaf P growth from storage (gP/m2/s)
    real(r8), pointer :: frootp_xfer_to_frootp               (:)     ! fine root P growth from storage (gP/m2/s)
    real(r8), pointer :: livestemp_xfer_to_livestemp         (:)     ! live stem P growth from storage (gP/m2/s)
    real(r8), pointer :: deadstemp_xfer_to_deadstemp         (:)     ! dead stem P growth from storage (gP/m2/s)
    real(r8), pointer :: livecrootp_xfer_to_livecrootp       (:)     ! live coarse root P growth from storage (gP/m2/s)
    real(r8), pointer :: deadcrootp_xfer_to_deadcrootp       (:)     ! dead coarse root P growth from storage (gP/m2/s)
    real(r8), pointer :: livestemp_to_litter                 (:)     ! livestem P to litter (gP/m2/s)
    real(r8), pointer :: grainp_to_food                      (:)     ! grain P to food for prognostic crop (gP/m2/s)
    real(r8), pointer :: leafp_to_litter                     (:)     ! leaf P litterfall (gP/m2/s)
    real(r8), pointer :: leafp_to_retransp                   (:)     ! leaf P to retranslocated P pool (gP/m2/s)
    real(r8), pointer :: frootp_to_retransp                  (:)     ! fine root P to retranslocated P pool (gP/m2/s)
    real(r8), pointer :: frootp_to_litter                    (:)     ! fine root P litterfall (gP/m2/s)
    real(r8), pointer :: retransp_to_ppool                   (:)     ! deployment of retranslocated P (gP/m2/s)       
    real(r8), pointer :: sminp_to_ppool                      (:)     ! deployment of soil mineral P uptake (gP/m2/s)
    real(r8), pointer :: ppool_to_grainp                     (:)     ! allocation to grain P for prognostic crop (gP/m2/s)
    real(r8), pointer :: ppool_to_grainp_storage             (:)     ! allocation to grain P storage for prognostic crop (gP/m2/s)
    real(r8), pointer :: ppool_to_leafp                      (:)     ! allocation to leaf P (gP/m2/s)
    real(r8), pointer :: ppool_to_leafp_storage              (:)     ! allocation to leaf P storage (gP/m2/s)
    real(r8), pointer :: ppool_to_frootp                     (:)     ! allocation to fine root P (gP/m2/s)
    real(r8), pointer :: ppool_to_frootp_storage             (:)     ! allocation to fine root P storage (gP/m2/s)
    real(r8), pointer :: ppool_to_livestemp                  (:)     ! allocation to live stem P (gP/m2/s)
    real(r8), pointer :: ppool_to_livestemp_storage          (:)     ! allocation to live stem P storage (gP/m2/s)
    real(r8), pointer :: ppool_to_deadstemp                  (:)     ! allocation to dead stem P (gP/m2/s)
    real(r8), pointer :: ppool_to_deadstemp_storage          (:)     ! allocation to dead stem P storage (gP/m2/s)
    real(r8), pointer :: ppool_to_livecrootp                 (:)     ! allocation to live coarse root P (gP/m2/s)
    real(r8), pointer :: ppool_to_livecrootp_storage         (:)     ! allocation to live coarse root P storage (gP/m2/s)
    real(r8), pointer :: ppool_to_deadcrootp                 (:)     ! allocation to dead coarse root P (gP/m2/s)
    real(r8), pointer :: ppool_to_deadcrootp_storage         (:)     ! allocation to dead coarse root P storage (gP/m2/s)
    real(r8), pointer :: grainp_storage_to_xfer              (:)     ! grain P shift storage to transfer for prognostic crop (gP/m2/s)
    real(r8), pointer :: leafp_storage_to_xfer               (:)     ! leaf P shift storage to transfer (gP/m2/s)
    real(r8), pointer :: frootp_storage_to_xfer              (:)     ! fine root P shift storage to transfer (gP/m2/s)
    real(r8), pointer :: livestemp_storage_to_xfer           (:)     ! live stem P shift storage to transfer (gP/m2/s)
    real(r8), pointer :: deadstemp_storage_to_xfer           (:)     ! dead stem P shift storage to transfer (gP/m2/s)
    real(r8), pointer :: livecrootp_storage_to_xfer          (:)     ! live coarse root P shift storage to transfer (gP/m2/s)
    real(r8), pointer :: deadcrootp_storage_to_xfer          (:)     ! dead coarse root P shift storage to transfer (gP/m2/s)
    real(r8), pointer :: fert_p                              (:)     ! applied fertilizer (gP/m2/s)
    real(r8), pointer :: fert_p_counter                      (:)     ! >0 fertilize; <=0 not
    real(r8), pointer :: livestemp_to_deadstemp              (:)     ! live stem P turnover (gP/m2/s)
    real(r8), pointer :: livestemp_to_retransp               (:)     ! live stem P to retranslocated P pool (gP/m2/s)
    real(r8), pointer :: livecrootp_to_deadcrootp            (:)     ! live coarse root P turnover (gP/m2/s)
    real(r8), pointer :: livecrootp_to_retransp              (:)     ! live coarse root P to retranslocated P pool (gP/m2/s)
    real(r8), pointer :: pdeploy                             (:)     ! total P deployed to growth and storage (gP/m2/s)
    real(r8), pointer :: pinputs                             (:)     ! total P inputs to pft-level (gP/m2/s)
    real(r8), pointer :: poutputs                            (:)     ! total P outputs from pft-level (gP/m2/s)
    real(r8), pointer :: wood_harvestp                       (:)     ! total P losses to wood product pools (gP/m2/s)
    real(r8), pointer :: biochem_pmin_to_plant               (:)     ! biochemical P mineralization directly goes to plant (gP/m2/s)
    real(r8), pointer :: crop_seedp_to_leaf                  (:)     ! (gP/m2/s) seed source to leaf, for crops
    real(r8), pointer :: dwt_seedp_to_leaf                   (:)     ! (gP/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_seedp_to_deadstem               (:)     ! (gP/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_conv_pflux                      (:)     ! (gP/m2/s) conversion N flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_prod10p_gain                    (:)     ! (gP/m2/s) addition to 10-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_prod100p_gain                   (:)     ! (gP/m2/s) addition to 100-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_crop_productp_gain              (:)     ! (gP/m2/s) addition to crop product pool from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), pointer :: dwt_seedp_to_ppool                  (:)     ! (gP/m2/s) seed source to PFT-level
    real(r8), pointer :: plant_pdemand                       (:)     ! P flux required to support initial GPP (gP/m2/s)
    real(r8), pointer :: avail_retransp                      (:)     ! P flux available from retranslocation pool (gP/m2/s)
    real(r8), pointer :: plant_palloc                        (:)     ! total allocated P flux (gP/m2/s)
    real(r8), pointer :: sminp_to_plant                      (:)     ! plant p uptake (gP/m2/s)
    real(r8), pointer :: plant_pdemand_vr                    (:,:)   ! plant P demand
    real(r8), pointer :: prev_leafp_to_litter                (:)     ! previous timestep leaf P litterfall flux (gP/m2/s)
    real(r8), pointer :: prev_frootp_to_litter               (:)     ! previous timestep froot P litterfall flux (gP/m2/s)
    real(r8), pointer :: supplement_to_plantp                (:)     ! supplementary P flux for plant 
    real(r8), pointer :: gap_ploss_litter                    (:)     ! total ploss from veg to litter pool due to gap mortality
    real(r8), pointer :: fire_ploss_litter                   (:)     ! total ploss from veg to litter pool due to fire
    real(r8), pointer :: hrv_ploss_litter                    (:)     ! total ploss from veg to litter pool due to harvest mortality
    real(r8), pointer :: sen_ploss_litter                    (:)     ! total ploss from veg to litter pool due to senescence
  contains
    procedure, public :: Init      => veg_pf_init
  end type vegetation_phosphorus_flux

   
  !-----------------------------------------------------------------------
  ! declare the public instances of vegetation-level data types
  !-----------------------------------------------------------------------
  type(vegetation_energy_state)          , public, target :: veg_es     ! vegetation energy state
  type(vegetation_water_state)           , public, target :: veg_ws     ! vegetation water state
  type(vegetation_carbon_state)          , public, target :: veg_cs     ! vegetation carbon state
  type(vegetation_carbon_state)          , public, target :: c13_veg_cs ! vegetation carbon state (C13)
  type(vegetation_carbon_state)          , public, target :: c14_veg_cs ! vegetation carbon state (C14)
  type(vegetation_nitrogen_state)        , public, target :: veg_ns     ! vegetation nitrogen state
  type(vegetation_phosphorus_state)      , public, target :: veg_ps     ! vegetation phosphorus state
  type(vegetation_energy_flux)           , public, target :: veg_ef     ! vegetation energy flux
  type(vegetation_water_flux)            , public, target :: veg_wf     ! vegetation water flux
  type(vegetation_carbon_flux)           , public, target :: veg_cf     ! vegetation carbon flux
  type(vegetation_carbon_flux)           , public, target :: c13_veg_cf ! vegetation carbon flux (C13)
  type(vegetation_carbon_flux)           , public, target :: c14_veg_cf ! vegetation carbon flux (C14)
  type(vegetation_nitrogen_flux)         , public, target :: veg_nf     ! vegetation nitrogen flux
  type(vegetation_phosphorus_flux)       , public, target :: veg_pf     ! vegetation phosphorus flux

  !------------------------------------------------------------------------
  
  contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy state data structure
  !------------------------------------------------------------------------
  subroutine veg_es_init(this, begp, endp)
    !
    ! !USES:
    use clm_varctl     , only : use_vancouver, use_mexicocity
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,j                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_es
    !-----------------------------------------------------------------------
    allocate(this%t_veg              (begp:endp))                   ; this%t_veg              (:)   = nan
    allocate(this%t_ref2m            (begp:endp))                   ; this%t_ref2m            (:)   = nan
    allocate(this%t_ref2m_r          (begp:endp))                   ; this%t_ref2m_r          (:)   = nan
    allocate(this%t_ref2m_u          (begp:endp))                   ; this%t_ref2m_u          (:)   = nan
    allocate(this%t_a10              (begp:endp))                   ; this%t_a10              (:)   = nan
    allocate(this%t_a10min           (begp:endp))                   ; this%t_a10min           (:)   = nan
    allocate(this%t_a5min            (begp:endp))                   ; this%t_a5min            (:)   = nan
    allocate(this%t_ref2m_min        (begp:endp))                   ; this%t_ref2m_min        (:)   = nan
    allocate(this%t_ref2m_min_r      (begp:endp))                   ; this%t_ref2m_min_r      (:)   = nan
    allocate(this%t_ref2m_min_u      (begp:endp))                   ; this%t_ref2m_min_u      (:)   = nan
    allocate(this%t_ref2m_max        (begp:endp))                   ; this%t_ref2m_max        (:)   = nan
    allocate(this%t_ref2m_max_r      (begp:endp))                   ; this%t_ref2m_max_r      (:)   = nan
    allocate(this%t_ref2m_max_u      (begp:endp))                   ; this%t_ref2m_max_u      (:)   = nan
    allocate(this%t_ref2m_min_inst   (begp:endp))                   ; this%t_ref2m_min_inst   (:)   = nan
    allocate(this%t_ref2m_min_inst_r (begp:endp))                   ; this%t_ref2m_min_inst_r (:)   = nan
    allocate(this%t_ref2m_min_inst_u (begp:endp))                   ; this%t_ref2m_min_inst_u (:)   = nan
    allocate(this%t_ref2m_max_inst   (begp:endp))                   ; this%t_ref2m_max_inst   (:)   = nan
    allocate(this%t_ref2m_max_inst_r (begp:endp))                   ; this%t_ref2m_max_inst_r (:)   = nan
    allocate(this%t_ref2m_max_inst_u (begp:endp))                   ; this%t_ref2m_max_inst_u (:)   = nan
    allocate(this%t_veg24            (begp:endp))                   ; this%t_veg24            (:)   = nan
    allocate(this%t_veg240           (begp:endp))                   ; this%t_veg240           (:)   = nan
    allocate(this%gdd0               (begp:endp))                   ; this%gdd0               (:)   = spval
    allocate(this%gdd8               (begp:endp))                   ; this%gdd8               (:)   = spval
    allocate(this%gdd10              (begp:endp))                   ; this%gdd10              (:)   = spval
    allocate(this%gdd020             (begp:endp))                   ; this%gdd020             (:)   = spval
    allocate(this%gdd820             (begp:endp))                   ; this%gdd820             (:)   = spval
    allocate(this%gdd1020            (begp:endp))                   ; this%gdd1020            (:)   = spval
    allocate(this%thm                (begp:endp))                   ; this%thm                (:)   = nan
    allocate(this%emv                (begp:endp))                   ; this%emv                (:)   = nan

  end subroutine veg_es_init
    

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water state data structure
  !------------------------------------------------------------------------
  subroutine veg_ws_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer             :: p
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of veg_ws
    !-----------------------------------------------------------------------
    allocate(this%h2ocan              (begp:endp))          ; this%h2ocan            (:) = nan
    allocate(this%q_ref2m             (begp:endp))          ; this%q_ref2m           (:) = nan
    allocate(this%rh_ref2m            (begp:endp))          ; this%rh_ref2m          (:) = nan
    allocate(this%rh_ref2m_r          (begp:endp))          ; this%rh_ref2m_r        (:) = nan
    allocate(this%rh_ref2m_u          (begp:endp))          ; this%rh_ref2m_u        (:) = nan
    allocate(this%rh_af               (begp:endp))          ; this%rh_af             (:) = nan
    allocate(this%fwet                (begp:endp))          ; this%fwet              (:) = nan
    allocate(this%fdry                (begp:endp))          ; this%fdry              (:) = nan
    allocate(this%begwb               (begp:endp))          ; this%begwb             (:) = nan
    allocate(this%endwb               (begp:endp))          ; this%endwb             (:) = nan
    allocate(this%errh2o              (begp:endp))          ; this%errh2o            (:) = nan

  end subroutine veg_ws_init
    
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon state data structure
  !------------------------------------------------------------------------
  subroutine veg_cs_init(this, begp, endp, carbon_type, ratio)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    integer          , intent(in) :: begp,endp
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    real(r8)         , intent(in) :: ratio
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,i,j,k
    integer :: fp                       ! filter index
    integer :: num_special_veg          ! number of good values in special_veg filter
    integer :: special_veg(endp-begp+1) ! special landunit filter - veg
    real(r8):: value_veg                ! used to reset state variables on special landunits
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_cs
    !-----------------------------------------------------------------------
    if ( .not. use_fates ) then
       allocate(this%leafc              (begp :endp))   ;  this%leafc              (:)   = nan
       allocate(this%leafc_storage      (begp :endp))   ;  this%leafc_storage      (:)   = nan
       allocate(this%leafc_xfer         (begp :endp))   ;  this%leafc_xfer         (:)   = nan
       allocate(this%frootc             (begp :endp))   ;  this%frootc             (:)   = nan
       allocate(this%frootc_storage     (begp :endp))   ;  this%frootc_storage     (:)   = nan
       allocate(this%frootc_xfer        (begp :endp))   ;  this%frootc_xfer        (:)   = nan
       allocate(this%livestemc          (begp :endp))   ;  this%livestemc          (:)   = nan
       allocate(this%livestemc_storage  (begp :endp))   ;  this%livestemc_storage  (:)   = nan
       allocate(this%livestemc_xfer     (begp :endp))   ;  this%livestemc_xfer     (:)   = nan
       allocate(this%deadstemc          (begp :endp))   ;  this%deadstemc          (:)   = nan
       allocate(this%deadstemc_storage  (begp :endp))   ;  this%deadstemc_storage  (:)   = nan
       allocate(this%deadstemc_xfer     (begp :endp))   ;  this%deadstemc_xfer     (:)   = nan
       allocate(this%livecrootc         (begp :endp))   ;  this%livecrootc         (:)   = nan
       allocate(this%livecrootc_storage (begp :endp))   ;  this%livecrootc_storage (:)   = nan
       allocate(this%livecrootc_xfer    (begp :endp))   ;  this%livecrootc_xfer    (:)   = nan
       allocate(this%deadcrootc         (begp :endp))   ;  this%deadcrootc         (:)   = nan
       allocate(this%deadcrootc_storage (begp :endp))   ;  this%deadcrootc_storage (:)   = nan
       allocate(this%deadcrootc_xfer    (begp :endp))   ;  this%deadcrootc_xfer    (:)   = nan
       allocate(this%gresp_storage      (begp :endp))   ;  this%gresp_storage      (:)   = nan
       allocate(this%gresp_xfer         (begp :endp))   ;  this%gresp_xfer         (:)   = nan
       allocate(this%cpool              (begp :endp))   ;  this%cpool              (:)   = nan
       allocate(this%xsmrpool           (begp :endp))   ;  this%xsmrpool           (:)   = nan
       allocate(this%ctrunc             (begp :endp))   ;  this%ctrunc             (:)   = nan
       allocate(this%dispvegc           (begp :endp))   ;  this%dispvegc           (:)   = nan
       allocate(this%storvegc           (begp :endp))   ;  this%storvegc           (:)   = nan
       allocate(this%totvegc            (begp :endp))   ;  this%totvegc            (:)   = nan
       allocate(this%totpftc            (begp :endp))   ;  this%totpftc            (:)   = nan
       allocate(this%leafcmax           (begp :endp))   ;  this%leafcmax           (:)   = nan
       allocate(this%grainc             (begp :endp))   ;  this%grainc             (:)   = nan
       allocate(this%grainc_storage     (begp :endp))   ;  this%grainc_storage     (:)   = nan
       allocate(this%grainc_xfer        (begp :endp))   ;  this%grainc_xfer        (:)   = nan
       allocate(this%woodc              (begp :endp))   ;  this%woodc              (:)   = nan
       allocate(this%totvegc_abg        (begp :endp))   ;  this%totvegc_abg        (:)   = nan
    endif  !  not use_fates

    allocate(this%begcb              (begp :endp))   ;  this%begcb              (:) = nan
    allocate(this%endcb              (begp :endp))   ;  this%endcb              (:) = nan
    allocate(this%errcb              (begp :endp))   ;  this%errcb              (:) = nan
    allocate(this%cropseedc_deficit  (begp :endp))   ;  this%cropseedc_deficit  (:) = nan

    end subroutine veg_cs_init
    
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine veg_ns_init(this, begp, endp, veg_cs)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    integer, intent(in)              :: begp,endp
    type(vegetation_carbon_state)    :: veg_cs ! used with C:N ratios on cold start
    !
    ! !LOCAL VARIABLES:
    integer :: fp,l,c,p                                 ! indices
    integer :: num_special_veg                          ! number of good values in special_patch filter
    integer :: special_veg (endp-begp+1)  ! special landunit filter - patches
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_ns
    !-----------------------------------------------------------------------
    allocate(this%leafn                  (begp:endp))           ; this%leafn               (:)   = nan
    allocate(this%leafn_storage          (begp:endp))           ; this%leafn_storage       (:)   = nan
    allocate(this%leafn_xfer             (begp:endp))           ; this%leafn_xfer          (:)   = nan
    allocate(this%frootn                 (begp:endp))           ; this%frootn              (:)   = nan
    allocate(this%frootn_storage         (begp:endp))           ; this%frootn_storage      (:)   = nan
    allocate(this%frootn_xfer            (begp:endp))           ; this%frootn_xfer         (:)   = nan
    allocate(this%livestemn              (begp:endp))           ; this%livestemn           (:)   = nan
    allocate(this%livestemn_storage      (begp:endp))           ; this%livestemn_storage   (:)   = nan
    allocate(this%livestemn_xfer         (begp:endp))           ; this%livestemn_xfer      (:)   = nan
    allocate(this%deadstemn              (begp:endp))           ; this%deadstemn           (:)   = nan
    allocate(this%deadstemn_storage      (begp:endp))           ; this%deadstemn_storage   (:)   = nan
    allocate(this%deadstemn_xfer         (begp:endp))           ; this%deadstemn_xfer      (:)   = nan
    allocate(this%livecrootn             (begp:endp))           ; this%livecrootn          (:)   = nan
    allocate(this%livecrootn_storage     (begp:endp))           ; this%livecrootn_storage  (:)   = nan
    allocate(this%livecrootn_xfer        (begp:endp))           ; this%livecrootn_xfer     (:)   = nan
    allocate(this%deadcrootn             (begp:endp))           ; this%deadcrootn          (:)   = nan
    allocate(this%deadcrootn_storage     (begp:endp))           ; this%deadcrootn_storage  (:)   = nan
    allocate(this%deadcrootn_xfer        (begp:endp))           ; this%deadcrootn_xfer     (:)   = nan
    allocate(this%retransn               (begp:endp))           ; this%retransn            (:)   = nan
    allocate(this%npool                  (begp:endp))           ; this%npool               (:)   = nan
    allocate(this%ntrunc                 (begp:endp))           ; this%ntrunc              (:)   = nan
    allocate(this%plant_n_buffer         (begp:endp))           ; this%plant_n_buffer      (:)   = nan
    allocate(this%grainn                 (begp:endp))           ; this%grainn              (:)   = nan
    allocate(this%grainn_storage         (begp:endp))           ; this%grainn_storage      (:)   = nan
    allocate(this%grainn_xfer            (begp:endp))           ; this%grainn_xfer         (:)   = nan
    allocate(this%cropseedn_deficit      (begp:endp))           ; this%cropseedn_deficit   (:)   = nan
    allocate(this%dispvegn               (begp:endp))           ; this%dispvegn            (:)   = nan
    allocate(this%storvegn               (begp:endp))           ; this%storvegn            (:)   = nan
    allocate(this%totvegn                (begp:endp))           ; this%totvegn             (:)   = nan
    allocate(this%totpftn                (begp:endp))           ; this%totpftn             (:)   = nan
    allocate(this%begnb                  (begp:endp))           ; this%begnb               (:)   = nan
    allocate(this%endnb                  (begp:endp))           ; this%endnb               (:)   = nan
    allocate(this%errnb                  (begp:endp))           ; this%errnb               (:)   = nan
    allocate(this%npimbalance            (begp:endp))           ; this%npimbalance         (:)   = nan
    allocate(this%pnup_pfrootc           (begp:endp))           ; this%pnup_pfrootc        (:)   = nan
    allocate(this%ppup_pfrootc           (begp:endp))           ; this%ppup_pfrootc        (:)   = nan
    allocate(this%ptlai_pleafc           (begp:endp))           ; this%ptlai_pleafc        (:)   = nan
    allocate(this%ppsnsun_ptlai          (begp:endp))           ; this%ppsnsun_ptlai       (:)   = nan
    allocate(this%ppsnsun_pleafn         (begp:endp))           ; this%ppsnsun_pleafn      (:)   = nan
    allocate(this%ppsnsun_pleafp         (begp:endp))           ; this%ppsnsun_pleafp      (:)   = nan
    allocate(this%plmrsun_ptlai          (begp:endp))           ; this%plmrsun_ptlai       (:)   = nan
    allocate(this%plmrsun_pleafn         (begp:endp))           ; this%plmrsun_pleafn      (:)   = nan
    allocate(this%plaisun_ptlai          (begp:endp))           ; this%plaisun_ptlai       (:)   = nan
    allocate(this%ppsnsha_ptlai          (begp:endp))           ; this%ppsnsha_ptlai       (:)   = nan
    allocate(this%ppsnsha_pleafn         (begp:endp))           ; this%ppsnsha_pleafn      (:)   = nan
    allocate(this%ppsnsha_pleafp         (begp:endp))           ; this%ppsnsha_pleafp      (:)   = nan
    allocate(this%plmrsha_ptlai          (begp:endp))           ; this%plmrsha_ptlai       (:)   = nan
    allocate(this%plmrsha_pleafn         (begp:endp))           ; this%plmrsha_pleafn      (:)   = nan
    allocate(this%plaisha_ptlai          (begp:endp))           ; this%plaisha_ptlai       (:)   = nan
    allocate(this%benefit_pgpp_pleafc    (begp:endp))           ; this%benefit_pgpp_pleafc (:)   = nan
    allocate(this%benefit_pgpp_pleafn    (begp:endp))           ; this%benefit_pgpp_pleafn (:)   = nan
    allocate(this%benefit_pgpp_pleafp    (begp:endp))           ; this%benefit_pgpp_pleafp (:)   = nan
    allocate(this%cost_pgpp_pfrootc      (begp:endp))           ; this%cost_pgpp_pfrootc   (:)   = nan
    allocate(this%cost_plmr_pleafc       (begp:endp))           ; this%cost_plmr_pleafc    (:)   = nan
    allocate(this%cost_plmr_pleafn       (begp:endp))           ; this%cost_plmr_pleafn    (:)   = nan
    allocate(this%ppsn_ptlai_z           (begp:endp,1:nlevcan)) ; this%ppsn_ptlai_z        (:,:) = nan
    allocate(this%ppsn_pleafn_z          (begp:endp,1:nlevcan)) ; this%ppsn_pleafn_z       (:,:) = nan
    allocate(this%ppsn_pleafp_z          (begp:endp,1:nlevcan)) ; this%ppsn_pleafp_z       (:,:) = nan
    allocate(this%ppsn_ptlai_z_vcmax     (begp:endp,1:nlevcan)) ; this%ppsn_ptlai_z_vcmax  (:,:) = nan
    allocate(this%ppsn_pleafn_z_vcmax    (begp:endp,1:nlevcan)) ; this%ppsn_pleafn_z_vcmax (:,:) = nan
    allocate(this%ppsn_pleafp_z_vcmax    (begp:endp,1:nlevcan)) ; this%ppsn_pleafp_z_vcmax (:,:) = nan
    allocate(this%ppsn_ptlai_z_jmax      (begp:endp,1:nlevcan)) ; this%ppsn_ptlai_z_jmax   (:,:) = nan
    allocate(this%ppsn_pleafn_z_jmax     (begp:endp,1:nlevcan)) ; this%ppsn_pleafn_z_jmax  (:,:) = nan
    allocate(this%ppsn_pleafp_z_jmax     (begp:endp,1:nlevcan)) ; this%ppsn_pleafp_z_jmax  (:,:) = nan
    allocate(this%ppsn_ptlai_z_tpu       (begp:endp,1:nlevcan)) ; this%ppsn_ptlai_z_tpu    (:,:) = nan
    allocate(this%ppsn_pleafn_z_tpu      (begp:endp,1:nlevcan)) ; this%ppsn_pleafn_z_tpu   (:,:) = nan
    allocate(this%ppsn_pleafp_z_tpu      (begp:endp,1:nlevcan)) ; this%ppsn_pleafp_z_tpu   (:,:) = nan
    allocate(this%plmr_ptlai_z           (begp:endp,1:nlevcan)) ; this%plmr_ptlai_z        (:,:) = nan
    allocate(this%plmr_pleafn_z          (begp:endp,1:nlevcan)) ; this%plmr_pleafn_z       (:,:) = nan

  end subroutine veg_ns_init
    
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine veg_ps_init(this, begp, endp, veg_cs)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state)        :: this
    integer, intent(in)                       :: begp,endp
    type(vegetation_carbon_state), intent(in) :: veg_cs
    !
    ! !LOCAL VARIABLES:
    integer :: fp,l,p                      ! indices
    integer :: num_special_patch           ! number of good values in special_patch filter
    integer :: special_patch (endp-begp+1) ! special landunit filter - patches
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_ps
    !-----------------------------------------------------------------------
    allocate(this%grainp             (begp:endp)) ; this%grainp             (:) = nan     
    allocate(this%grainp_storage     (begp:endp)) ; this%grainp_storage     (:) = nan
    allocate(this%grainp_xfer        (begp:endp)) ; this%grainp_xfer        (:) = nan     
    allocate(this%leafp              (begp:endp)) ; this%leafp              (:) = nan
    allocate(this%leafp_storage      (begp:endp)) ; this%leafp_storage      (:) = nan     
    allocate(this%leafp_xfer         (begp:endp)) ; this%leafp_xfer         (:) = nan     
    allocate(this%frootp             (begp:endp)) ; this%frootp             (:) = nan
    allocate(this%frootp_storage     (begp:endp)) ; this%frootp_storage     (:) = nan     
    allocate(this%frootp_xfer        (begp:endp)) ; this%frootp_xfer        (:) = nan     
    allocate(this%livestemp          (begp:endp)) ; this%livestemp          (:) = nan
    allocate(this%livestemp_storage  (begp:endp)) ; this%livestemp_storage  (:) = nan
    allocate(this%livestemp_xfer     (begp:endp)) ; this%livestemp_xfer     (:) = nan
    allocate(this%deadstemp          (begp:endp)) ; this%deadstemp          (:) = nan
    allocate(this%deadstemp_storage  (begp:endp)) ; this%deadstemp_storage  (:) = nan
    allocate(this%deadstemp_xfer     (begp:endp)) ; this%deadstemp_xfer     (:) = nan
    allocate(this%livecrootp         (begp:endp)) ; this%livecrootp         (:) = nan
    allocate(this%livecrootp_storage (begp:endp)) ; this%livecrootp_storage (:) = nan
    allocate(this%livecrootp_xfer    (begp:endp)) ; this%livecrootp_xfer    (:) = nan
    allocate(this%deadcrootp         (begp:endp)) ; this%deadcrootp         (:) = nan
    allocate(this%deadcrootp_storage (begp:endp)) ; this%deadcrootp_storage (:) = nan
    allocate(this%deadcrootp_xfer    (begp:endp)) ; this%deadcrootp_xfer    (:) = nan
    allocate(this%retransp           (begp:endp)) ; this%retransp           (:) = nan
    allocate(this%ppool              (begp:endp)) ; this%ppool              (:) = nan
    allocate(this%ptrunc             (begp:endp)) ; this%ptrunc             (:) = nan
    allocate(this%dispvegp           (begp:endp)) ; this%dispvegp           (:) = nan
    allocate(this%storvegp           (begp:endp)) ; this%storvegp           (:) = nan
    allocate(this%totvegp            (begp:endp)) ; this%totvegp            (:) = nan
    allocate(this%totpftp            (begp:endp)) ; this%totpftp            (:) = nan
    allocate(this%plant_p_buffer     (begp:endp)) ; this%plant_p_buffer     (:) = nan
    allocate(this%cropseedp_deficit  (begp:endp)) ; this%cropseedp_deficit  (:) = nan
    allocate(this%begpb              (begp:endp)) ; this%begpb              (:) = nan
    allocate(this%endpb              (begp:endp)) ; this%endpb              (:) = nan
    allocate(this%errpb              (begp:endp)) ; this%errpb              (:) = nan
    

  end subroutine veg_ps_init

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy flux data structure
  !------------------------------------------------------------------------
  subroutine veg_ef_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer  :: l,c,p
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of veg_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_sh_grnd        (begp:endp))   ; this%eflx_sh_grnd       (:)   = nan
    allocate(this%eflx_sh_veg         (begp:endp))   ; this%eflx_sh_veg        (:)   = nan
    allocate(this%eflx_sh_snow        (begp:endp))   ; this%eflx_sh_snow       (:)   = nan
    allocate(this%eflx_sh_soil        (begp:endp))   ; this%eflx_sh_soil       (:)   = nan
    allocate(this%eflx_sh_h2osfc      (begp:endp))   ; this%eflx_sh_h2osfc     (:)   = nan
    allocate(this%eflx_sh_tot         (begp:endp))   ; this%eflx_sh_tot        (:)   = nan
    allocate(this%eflx_sh_tot_u       (begp:endp))   ; this%eflx_sh_tot_u      (:)   = nan
    allocate(this%eflx_sh_tot_r       (begp:endp))   ; this%eflx_sh_tot_r      (:)   = nan
    allocate(this%eflx_lh_tot         (begp:endp))   ; this%eflx_lh_tot        (:)   = nan
    allocate(this%eflx_lh_tot_u       (begp:endp))   ; this%eflx_lh_tot_u      (:)   = nan
    allocate(this%eflx_lh_tot_r       (begp:endp))   ; this%eflx_lh_tot_r      (:)   = nan
    allocate(this%eflx_lh_vegt        (begp:endp))   ; this%eflx_lh_vegt       (:)   = nan
    allocate(this%eflx_lh_vege        (begp:endp))   ; this%eflx_lh_vege       (:)   = nan
    allocate(this%eflx_lh_grnd        (begp:endp))   ; this%eflx_lh_grnd       (:)   = nan
    allocate(this%eflx_soil_grnd      (begp:endp))   ; this%eflx_soil_grnd     (:)   = nan
    allocate(this%eflx_soil_grnd_u    (begp:endp))   ; this%eflx_soil_grnd_u   (:)   = nan
    allocate(this%eflx_soil_grnd_r    (begp:endp))   ; this%eflx_soil_grnd_r   (:)   = nan
    allocate(this%eflx_lwrad_net      (begp:endp))   ; this%eflx_lwrad_net     (:)   = nan
    allocate(this%eflx_lwrad_net_r    (begp:endp))   ; this%eflx_lwrad_net_r   (:)   = nan
    allocate(this%eflx_lwrad_net_u    (begp:endp))   ; this%eflx_lwrad_net_u   (:)   = nan
    allocate(this%eflx_lwrad_out      (begp:endp))   ; this%eflx_lwrad_out     (:)   = nan
    allocate(this%eflx_lwrad_out_r    (begp:endp))   ; this%eflx_lwrad_out_r   (:)   = nan
    allocate(this%eflx_lwrad_out_u    (begp:endp))   ; this%eflx_lwrad_out_u   (:)   = nan
    allocate(this%eflx_gnet           (begp:endp))   ; this%eflx_gnet          (:)   = nan
    allocate(this%eflx_grnd_lake      (begp:endp))   ; this%eflx_grnd_lake     (:)   = nan
    allocate(this%eflx_anthro         (begp:endp))   ; this%eflx_anthro        (:)   = nan
    allocate(this%eflx_traffic        (begp:endp))   ; this%eflx_traffic       (:)   = nan
    allocate(this%eflx_wasteheat      (begp:endp))   ; this%eflx_wasteheat     (:)   = nan
    allocate(this%eflx_heat_from_ac   (begp:endp))   ; this%eflx_heat_from_ac  (:)   = nan
    allocate(this%dlrad               (begp:endp))   ; this%dlrad              (:)   = nan
    allocate(this%ulrad               (begp:endp))   ; this%ulrad              (:)   = nan
    allocate(this%taux                (begp:endp))   ; this%taux               (:)   = nan
    allocate(this%tauy                (begp:endp))   ; this%tauy               (:)   = nan
    allocate(this%dgnetdT             (begp:endp))   ; this%dgnetdT            (:)   = nan
    allocate(this%netrad              (begp:endp))   ; this%netrad             (:)   = nan
    allocate(this%cgrnd               (begp:endp))   ; this%cgrnd              (:)   = nan
    allocate(this%cgrndl              (begp:endp))   ; this%cgrndl             (:)   = nan
    allocate(this%cgrnds              (begp:endp))   ; this%cgrnds             (:)   = nan
    allocate(this%errsoi              (begp:endp))   ; this%errsoi             (:)   = nan
    allocate(this%errseb              (begp:endp))   ; this%errseb             (:)   = nan
    allocate(this%errsol              (begp:endp))   ; this%errsol             (:)   = nan
    allocate(this%errlon              (begp:endp))   ; this%errlon             (:)   = nan

  end subroutine veg_ef_init
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water flux data structure
  !------------------------------------------------------------------------
  subroutine veg_wf_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,l                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_wf
    !-----------------------------------------------------------------------
    allocate(this%qflx_prec_grnd         (begp:endp))             ; this%qflx_prec_grnd       (:)   = nan
    allocate(this%qflx_rain_grnd         (begp:endp))             ; this%qflx_rain_grnd       (:)   = nan
    allocate(this%qflx_snow_grnd         (begp:endp))             ; this%qflx_snow_grnd       (:)   = nan
    allocate(this%qflx_sub_snow          (begp:endp))             ; this%qflx_sub_snow        (:)   = 0._r8
    allocate(this%qflx_evap_soi          (begp:endp))             ; this%qflx_evap_soi        (:)   = nan
    allocate(this%qflx_evap_veg          (begp:endp))             ; this%qflx_evap_veg        (:)   = nan
    allocate(this%qflx_evap_can          (begp:endp))             ; this%qflx_evap_can        (:)   = nan
    allocate(this%qflx_evap_tot          (begp:endp))             ; this%qflx_evap_tot        (:)   = nan
    allocate(this%qflx_evap_grnd         (begp:endp))             ; this%qflx_evap_grnd       (:)   = nan
    allocate(this%qflx_snwcp_liq         (begp:endp))             ; this%qflx_snwcp_liq       (:)   = nan
    allocate(this%qflx_snwcp_ice         (begp:endp))             ; this%qflx_snwcp_ice       (:)   = nan
    allocate(this%qflx_tran_veg          (begp:endp))             ; this%qflx_tran_veg        (:)   = nan
    allocate(this%qflx_dew_snow          (begp:endp))             ; this%qflx_dew_snow        (:)   = nan
    allocate(this%qflx_dew_grnd          (begp:endp))             ; this%qflx_dew_grnd        (:)   = nan
    allocate(this%qflx_prec_intr         (begp:endp))             ; this%qflx_prec_intr       (:)   = nan
    allocate(this%qflx_ev_snow           (begp:endp))             ; this%qflx_ev_snow         (:)   = nan
    allocate(this%qflx_ev_soil           (begp:endp))             ; this%qflx_ev_soil         (:)   = nan
    allocate(this%qflx_ev_h2osfc         (begp:endp))             ; this%qflx_ev_h2osfc       (:)   = nan
    allocate(this%qflx_rootsoi_frac      (begp:endp,1:nlevgrnd))  ; this%qflx_rootsoi_frac    (:,:) = nan
    allocate(this%qflx_irrig             (begp:endp))             ; this%qflx_irrig           (:)   = nan
    allocate(this%irrig_rate             (begp:endp))             ; this%irrig_rate           (:)   = nan
    allocate(this%n_irrig_steps_left     (begp:endp))             ; this%n_irrig_steps_left   (:)   = 0
    
  end subroutine veg_wf_init
    
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon flux data structure
  !------------------------------------------------------------------------
  subroutine veg_cf_init(this, begp, endp, carbon_type)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    integer, intent(in)           :: begp,endp
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,j                     ! indices
    integer :: num_special_patch           ! number of good values in special_patch filter
    integer :: special_patch(endp-begp+1)  ! special landunit filter - patches
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_cf
    !-----------------------------------------------------------------------
    if (.not. use_fates) then
       allocate(this%m_leafc_to_litter                   (begp:endp)) ;    this%m_leafc_to_litter                    (:) = nan     
       allocate(this%m_leafc_storage_to_litter           (begp:endp)) ;    this%m_leafc_storage_to_litter            (:) = nan
       allocate(this%m_leafc_xfer_to_litter              (begp:endp)) ;    this%m_leafc_xfer_to_litter               (:) = nan
       allocate(this%m_frootc_to_litter                  (begp:endp)) ;    this%m_frootc_to_litter                   (:) = nan
       allocate(this%m_frootc_storage_to_litter          (begp:endp)) ;    this%m_frootc_storage_to_litter           (:) = nan
       allocate(this%m_frootc_xfer_to_litter             (begp:endp)) ;    this%m_frootc_xfer_to_litter              (:) = nan
       allocate(this%m_livestemc_to_litter               (begp:endp)) ;    this%m_livestemc_to_litter                (:) = nan
       allocate(this%m_livestemc_storage_to_litter       (begp:endp)) ;    this%m_livestemc_storage_to_litter        (:) = nan
       allocate(this%m_livestemc_xfer_to_litter          (begp:endp)) ;    this%m_livestemc_xfer_to_litter           (:) = nan
       allocate(this%m_deadstemc_to_litter               (begp:endp)) ;    this%m_deadstemc_to_litter                (:) = nan
       allocate(this%m_deadstemc_storage_to_litter       (begp:endp)) ;    this%m_deadstemc_storage_to_litter        (:) = nan
       allocate(this%m_deadstemc_xfer_to_litter          (begp:endp)) ;    this%m_deadstemc_xfer_to_litter           (:) = nan
       allocate(this%m_livecrootc_to_litter              (begp:endp)) ;    this%m_livecrootc_to_litter               (:) = nan
       allocate(this%m_livecrootc_storage_to_litter      (begp:endp)) ;    this%m_livecrootc_storage_to_litter       (:) = nan
       allocate(this%m_livecrootc_xfer_to_litter         (begp:endp)) ;    this%m_livecrootc_xfer_to_litter          (:) = nan
       allocate(this%m_deadcrootc_to_litter              (begp:endp)) ;    this%m_deadcrootc_to_litter               (:) = nan
       allocate(this%m_deadcrootc_storage_to_litter      (begp:endp)) ;    this%m_deadcrootc_storage_to_litter       (:) = nan
       allocate(this%m_deadcrootc_xfer_to_litter         (begp:endp)) ;    this%m_deadcrootc_xfer_to_litter          (:) = nan
       allocate(this%m_gresp_storage_to_litter           (begp:endp)) ;    this%m_gresp_storage_to_litter            (:) = nan
       allocate(this%m_gresp_xfer_to_litter              (begp:endp)) ;    this%m_gresp_xfer_to_litter               (:) = nan
       allocate(this%m_cpool_to_litter                   (begp:endp)) ;    this%m_cpool_to_litter                    (:) = nan
       allocate(this%hrv_leafc_to_litter                 (begp:endp)) ;    this%hrv_leafc_to_litter                  (:) = nan
       allocate(this%hrv_leafc_storage_to_litter         (begp:endp)) ;    this%hrv_leafc_storage_to_litter          (:) = nan
       allocate(this%hrv_leafc_xfer_to_litter            (begp:endp)) ;    this%hrv_leafc_xfer_to_litter             (:) = nan
       allocate(this%hrv_frootc_to_litter                (begp:endp)) ;    this%hrv_frootc_to_litter                 (:) = nan
       allocate(this%hrv_frootc_storage_to_litter        (begp:endp)) ;    this%hrv_frootc_storage_to_litter         (:) = nan
       allocate(this%hrv_frootc_xfer_to_litter           (begp:endp)) ;    this%hrv_frootc_xfer_to_litter            (:) = nan
       allocate(this%hrv_livestemc_to_litter             (begp:endp)) ;    this%hrv_livestemc_to_litter              (:) = nan
       allocate(this%hrv_livestemc_storage_to_litter     (begp:endp)) ;    this%hrv_livestemc_storage_to_litter      (:) = nan
       allocate(this%hrv_livestemc_xfer_to_litter        (begp:endp)) ;    this%hrv_livestemc_xfer_to_litter         (:) = nan
       allocate(this%hrv_deadstemc_to_prod10c            (begp:endp)) ;    this%hrv_deadstemc_to_prod10c             (:) = nan
       allocate(this%hrv_deadstemc_to_prod100c           (begp:endp)) ;    this%hrv_deadstemc_to_prod100c            (:) = nan
       allocate(this%hrv_deadstemc_storage_to_litter     (begp:endp)) ;    this%hrv_deadstemc_storage_to_litter      (:) = nan
       allocate(this%hrv_deadstemc_xfer_to_litter        (begp:endp)) ;    this%hrv_deadstemc_xfer_to_litter         (:) = nan
       allocate(this%hrv_livecrootc_to_litter            (begp:endp)) ;    this%hrv_livecrootc_to_litter             (:) = nan
       allocate(this%hrv_livecrootc_storage_to_litter    (begp:endp)) ;    this%hrv_livecrootc_storage_to_litter     (:) = nan
       allocate(this%hrv_livecrootc_xfer_to_litter       (begp:endp)) ;    this%hrv_livecrootc_xfer_to_litter        (:) = nan
       allocate(this%hrv_deadcrootc_to_litter            (begp:endp)) ;    this%hrv_deadcrootc_to_litter             (:) = nan
       allocate(this%hrv_deadcrootc_storage_to_litter    (begp:endp)) ;    this%hrv_deadcrootc_storage_to_litter     (:) = nan
       allocate(this%hrv_deadcrootc_xfer_to_litter       (begp:endp)) ;    this%hrv_deadcrootc_xfer_to_litter        (:) = nan
       allocate(this%hrv_gresp_storage_to_litter         (begp:endp)) ;    this%hrv_gresp_storage_to_litter          (:) = nan
       allocate(this%hrv_gresp_xfer_to_litter            (begp:endp)) ;    this%hrv_gresp_xfer_to_litter             (:) = nan
       allocate(this%hrv_xsmrpool_to_atm                 (begp:endp)) ;    this%hrv_xsmrpool_to_atm                  (:) = nan
       allocate(this%hrv_cpool_to_litter                 (begp:endp)) ;    this%hrv_cpool_to_litter                  (:) = nan
       allocate(this%hrv_leafc_to_prod1c                 (begp:endp)) ;    this%hrv_leafc_to_prod1c                  (:) = nan
       allocate(this%hrv_livestemc_to_prod1c             (begp:endp)) ;    this%hrv_livestemc_to_prod1c              (:) = nan
       allocate(this%hrv_grainc_to_prod1c                (begp:endp)) ;    this%hrv_grainc_to_prod1c                 (:) = nan
       allocate(this%hrv_cropc_to_prod1c                 (begp:endp)) ;    this%hrv_cropc_to_prod1c                  (:) = nan
       allocate(this%m_leafc_to_fire                     (begp:endp)) ;    this%m_leafc_to_fire                      (:) = nan
       allocate(this%m_leafc_storage_to_fire             (begp:endp)) ;    this%m_leafc_storage_to_fire              (:) = nan
       allocate(this%m_leafc_xfer_to_fire                (begp:endp)) ;    this%m_leafc_xfer_to_fire                 (:) = nan
       allocate(this%m_livestemc_to_fire                 (begp:endp)) ;    this%m_livestemc_to_fire                  (:) = nan
       allocate(this%m_livestemc_storage_to_fire         (begp:endp)) ;    this%m_livestemc_storage_to_fire          (:) = nan
       allocate(this%m_livestemc_xfer_to_fire            (begp:endp)) ;    this%m_livestemc_xfer_to_fire             (:) = nan
       allocate(this%m_deadstemc_to_fire                 (begp:endp)) ;    this%m_deadstemc_to_fire                  (:) = nan
       allocate(this%m_deadstemc_storage_to_fire         (begp:endp)) ;    this%m_deadstemc_storage_to_fire          (:) = nan
       allocate(this%m_deadstemc_xfer_to_fire            (begp:endp)) ;    this%m_deadstemc_xfer_to_fire             (:) = nan
       allocate(this%m_frootc_to_fire                    (begp:endp)) ;    this%m_frootc_to_fire                     (:) = nan
       allocate(this%m_frootc_storage_to_fire            (begp:endp)) ;    this%m_frootc_storage_to_fire             (:) = nan
       allocate(this%m_frootc_xfer_to_fire               (begp:endp)) ;    this%m_frootc_xfer_to_fire                (:) = nan
       allocate(this%m_livecrootc_to_fire                (begp:endp)) ;    this%m_livecrootc_to_fire                 (:) = nan
       allocate(this%m_livecrootc_storage_to_fire        (begp:endp)) ;    this%m_livecrootc_storage_to_fire         (:) = nan
       allocate(this%m_livecrootc_xfer_to_fire           (begp:endp)) ;    this%m_livecrootc_xfer_to_fire            (:) = nan
       allocate(this%m_deadcrootc_to_fire                (begp:endp)) ;    this%m_deadcrootc_to_fire                 (:) = nan
       allocate(this%m_deadcrootc_storage_to_fire        (begp:endp)) ;    this%m_deadcrootc_storage_to_fire         (:) = nan
       allocate(this%m_deadcrootc_xfer_to_fire           (begp:endp)) ;    this%m_deadcrootc_xfer_to_fire            (:) = nan
       allocate(this%m_gresp_storage_to_fire             (begp:endp)) ;    this%m_gresp_storage_to_fire              (:) = nan
       allocate(this%m_gresp_xfer_to_fire                (begp:endp)) ;    this%m_gresp_xfer_to_fire                 (:) = nan
       allocate(this%m_cpool_to_fire                     (begp:endp)) ;    this%m_cpool_to_fire                      (:) = nan
       allocate(this%m_leafc_to_litter_fire              (begp:endp)) ;    this%m_leafc_to_litter_fire               (:) = nan
       allocate(this%m_leafc_storage_to_litter_fire      (begp:endp)) ;    this%m_leafc_storage_to_litter_fire       (:) = nan
       allocate(this%m_leafc_xfer_to_litter_fire         (begp:endp)) ;    this%m_leafc_xfer_to_litter_fire          (:) = nan
       allocate(this%m_livestemc_to_litter_fire          (begp:endp)) ;    this%m_livestemc_to_litter_fire           (:) = nan
       allocate(this%m_livestemc_storage_to_litter_fire  (begp:endp)) ;    this%m_livestemc_storage_to_litter_fire   (:) = nan
       allocate(this%m_livestemc_xfer_to_litter_fire     (begp:endp)) ;    this%m_livestemc_xfer_to_litter_fire      (:) = nan
       allocate(this%m_livestemc_to_deadstemc_fire       (begp:endp)) ;    this%m_livestemc_to_deadstemc_fire        (:) = nan
       allocate(this%m_deadstemc_to_litter_fire          (begp:endp)) ;    this%m_deadstemc_to_litter_fire           (:) = nan
       allocate(this%m_deadstemc_storage_to_litter_fire  (begp:endp)) ;    this%m_deadstemc_storage_to_litter_fire   (:) = nan
       allocate(this%m_deadstemc_xfer_to_litter_fire     (begp:endp)) ;    this%m_deadstemc_xfer_to_litter_fire      (:) = nan
       allocate(this%m_frootc_to_litter_fire             (begp:endp)) ;    this%m_frootc_to_litter_fire              (:) = nan
       allocate(this%m_frootc_storage_to_litter_fire     (begp:endp)) ;    this%m_frootc_storage_to_litter_fire      (:) = nan
       allocate(this%m_frootc_xfer_to_litter_fire        (begp:endp)) ;    this%m_frootc_xfer_to_litter_fire         (:) = nan
       allocate(this%m_livecrootc_to_litter_fire         (begp:endp)) ;    this%m_livecrootc_to_litter_fire          (:) = nan
       allocate(this%m_livecrootc_storage_to_litter_fire (begp:endp)) ;    this%m_livecrootc_storage_to_litter_fire  (:) = nan
       allocate(this%m_livecrootc_xfer_to_litter_fire    (begp:endp)) ;    this%m_livecrootc_xfer_to_litter_fire     (:) = nan
       allocate(this%m_livecrootc_to_deadcrootc_fire     (begp:endp)) ;    this%m_livecrootc_to_deadcrootc_fire      (:) = nan
       allocate(this%m_deadcrootc_to_litter_fire         (begp:endp)) ;    this%m_deadcrootc_to_litter_fire          (:) = nan
       allocate(this%m_deadcrootc_storage_to_litter_fire (begp:endp)) ;    this%m_deadcrootc_storage_to_litter_fire  (:) = nan
       allocate(this%m_deadcrootc_xfer_to_litter_fire    (begp:endp)) ;    this%m_deadcrootc_xfer_to_litter_fire     (:) = nan
       allocate(this%m_gresp_storage_to_litter_fire      (begp:endp)) ;    this%m_gresp_storage_to_litter_fire       (:) = nan
       allocate(this%m_gresp_xfer_to_litter_fire         (begp:endp)) ;    this%m_gresp_xfer_to_litter_fire          (:) = nan
       allocate(this%m_cpool_to_litter_fire              (begp:endp)) ;    this%m_cpool_to_litter_fire               (:) = nan
       allocate(this%grainc_xfer_to_grainc               (begp:endp)) ;    this%grainc_xfer_to_grainc                (:) = nan
       allocate(this%leafc_xfer_to_leafc                 (begp:endp)) ;    this%leafc_xfer_to_leafc                  (:) = nan
       allocate(this%frootc_xfer_to_frootc               (begp:endp)) ;    this%frootc_xfer_to_frootc                (:) = nan
       allocate(this%livestemc_xfer_to_livestemc         (begp:endp)) ;    this%livestemc_xfer_to_livestemc          (:) = nan
       allocate(this%deadstemc_xfer_to_deadstemc         (begp:endp)) ;    this%deadstemc_xfer_to_deadstemc          (:) = nan
       allocate(this%livecrootc_xfer_to_livecrootc       (begp:endp)) ;    this%livecrootc_xfer_to_livecrootc        (:) = nan
       allocate(this%deadcrootc_xfer_to_deadcrootc       (begp:endp)) ;    this%deadcrootc_xfer_to_deadcrootc        (:) = nan
       allocate(this%leafc_to_litter                     (begp:endp)) ;    this%leafc_to_litter                      (:) = nan
       allocate(this%frootc_to_litter                    (begp:endp)) ;    this%frootc_to_litter                     (:) = nan
       allocate(this%livestemc_to_litter                 (begp:endp)) ;    this%livestemc_to_litter                  (:) = nan
       allocate(this%grainc_to_food                      (begp:endp)) ;    this%grainc_to_food                       (:) = nan
       allocate(this%leaf_mr                             (begp:endp)) ;    this%leaf_mr                              (:) = nan
       allocate(this%froot_mr                            (begp:endp)) ;    this%froot_mr                             (:) = nan
       allocate(this%livestem_mr                         (begp:endp)) ;    this%livestem_mr                          (:) = nan
       allocate(this%livecroot_mr                        (begp:endp)) ;    this%livecroot_mr                         (:) = nan
       allocate(this%grain_mr                            (begp:endp)) ;    this%grain_mr                             (:) = nan
       allocate(this%leaf_curmr                          (begp:endp)) ;    this%leaf_curmr                           (:) = nan
       allocate(this%froot_curmr                         (begp:endp)) ;    this%froot_curmr                          (:) = nan
       allocate(this%livestem_curmr                      (begp:endp)) ;    this%livestem_curmr                       (:) = nan
       allocate(this%livecroot_curmr                     (begp:endp)) ;    this%livecroot_curmr                      (:) = nan
       allocate(this%grain_curmr                         (begp:endp)) ;    this%grain_curmr                          (:) = nan
       allocate(this%leaf_xsmr                           (begp:endp)) ;    this%leaf_xsmr                            (:) = nan
       allocate(this%froot_xsmr                          (begp:endp)) ;    this%froot_xsmr                           (:) = nan
       allocate(this%livestem_xsmr                       (begp:endp)) ;    this%livestem_xsmr                        (:) = nan
       allocate(this%livecroot_xsmr                      (begp:endp)) ;    this%livecroot_xsmr                       (:) = nan
       allocate(this%grain_xsmr                          (begp:endp)) ;    this%grain_xsmr                           (:) = nan
       allocate(this%xr                                  (begp:endp)) ;    this%xr                                   (:) = nan
       allocate(this%psnsun_to_cpool                     (begp:endp)) ;    this%psnsun_to_cpool                      (:) = nan
       allocate(this%psnshade_to_cpool                   (begp:endp)) ;    this%psnshade_to_cpool                    (:) = nan
       allocate(this%cpool_to_xsmrpool                   (begp:endp)) ;    this%cpool_to_xsmrpool                    (:) = nan
       allocate(this%cpool_to_grainc                     (begp:endp)) ;    this%cpool_to_grainc                      (:) = nan
       allocate(this%cpool_to_grainc_storage             (begp:endp)) ;    this%cpool_to_grainc_storage              (:) = nan
       allocate(this%cpool_to_leafc                      (begp:endp)) ;    this%cpool_to_leafc                       (:) = nan
       allocate(this%cpool_to_leafc_storage              (begp:endp)) ;    this%cpool_to_leafc_storage               (:) = nan
       allocate(this%cpool_to_frootc                     (begp:endp)) ;    this%cpool_to_frootc                      (:) = nan
       allocate(this%cpool_to_frootc_storage             (begp:endp)) ;    this%cpool_to_frootc_storage              (:) = nan
       allocate(this%cpool_to_livestemc                  (begp:endp)) ;    this%cpool_to_livestemc                   (:) = nan
       allocate(this%cpool_to_livestemc_storage          (begp:endp)) ;    this%cpool_to_livestemc_storage           (:) = nan
       allocate(this%cpool_to_deadstemc                  (begp:endp)) ;    this%cpool_to_deadstemc                   (:) = nan
       allocate(this%cpool_to_deadstemc_storage          (begp:endp)) ;    this%cpool_to_deadstemc_storage           (:) = nan
       allocate(this%cpool_to_livecrootc                 (begp:endp)) ;    this%cpool_to_livecrootc                  (:) = nan
       allocate(this%cpool_to_livecrootc_storage         (begp:endp)) ;    this%cpool_to_livecrootc_storage          (:) = nan
       allocate(this%cpool_to_deadcrootc                 (begp:endp)) ;    this%cpool_to_deadcrootc                  (:) = nan
       allocate(this%cpool_to_deadcrootc_storage         (begp:endp)) ;    this%cpool_to_deadcrootc_storage          (:) = nan
       allocate(this%cpool_to_gresp_storage              (begp:endp)) ;    this%cpool_to_gresp_storage               (:) = nan
       allocate(this%xsmrpool_to_atm                     (begp:endp)) ;    this%xsmrpool_to_atm                      (:) = nan
       allocate(this%cpool_leaf_gr                       (begp:endp)) ;    this%cpool_leaf_gr                        (:) = nan
       allocate(this%cpool_leaf_storage_gr               (begp:endp)) ;    this%cpool_leaf_storage_gr                (:) = nan
       allocate(this%transfer_leaf_gr                    (begp:endp)) ;    this%transfer_leaf_gr                     (:) = nan
       allocate(this%cpool_froot_gr                      (begp:endp)) ;    this%cpool_froot_gr                       (:) = nan
       allocate(this%cpool_froot_storage_gr              (begp:endp)) ;    this%cpool_froot_storage_gr               (:) = nan
       allocate(this%transfer_froot_gr                   (begp:endp)) ;    this%transfer_froot_gr                    (:) = nan
       allocate(this%cpool_livestem_gr                   (begp:endp)) ;    this%cpool_livestem_gr                    (:) = nan
       allocate(this%cpool_livestem_storage_gr           (begp:endp)) ;    this%cpool_livestem_storage_gr            (:) = nan
       allocate(this%transfer_livestem_gr                (begp:endp)) ;    this%transfer_livestem_gr                 (:) = nan
       allocate(this%cpool_deadstem_gr                   (begp:endp)) ;    this%cpool_deadstem_gr                    (:) = nan
       allocate(this%cpool_deadstem_storage_gr           (begp:endp)) ;    this%cpool_deadstem_storage_gr            (:) = nan
       allocate(this%transfer_deadstem_gr                (begp:endp)) ;    this%transfer_deadstem_gr                 (:) = nan
       allocate(this%cpool_livecroot_gr                  (begp:endp)) ;    this%cpool_livecroot_gr                   (:) = nan
       allocate(this%cpool_livecroot_storage_gr          (begp:endp)) ;    this%cpool_livecroot_storage_gr           (:) = nan
       allocate(this%transfer_livecroot_gr               (begp:endp)) ;    this%transfer_livecroot_gr                (:) = nan
       allocate(this%cpool_deadcroot_gr                  (begp:endp)) ;    this%cpool_deadcroot_gr                   (:) = nan
       allocate(this%cpool_deadcroot_storage_gr          (begp:endp)) ;    this%cpool_deadcroot_storage_gr           (:) = nan
       allocate(this%transfer_deadcroot_gr               (begp:endp)) ;    this%transfer_deadcroot_gr                (:) = nan
       allocate(this%cpool_grain_gr                      (begp:endp)) ;    this%cpool_grain_gr                       (:) = nan
       allocate(this%cpool_grain_storage_gr              (begp:endp)) ;    this%cpool_grain_storage_gr               (:) = nan
       allocate(this%transfer_grain_gr                   (begp:endp)) ;    this%transfer_grain_gr                    (:) = nan
       allocate(this%grainc_storage_to_xfer              (begp:endp)) ;    this%grainc_storage_to_xfer               (:) = nan
       allocate(this%leafc_storage_to_xfer               (begp:endp)) ;    this%leafc_storage_to_xfer                (:) = nan
       allocate(this%frootc_storage_to_xfer              (begp:endp)) ;    this%frootc_storage_to_xfer               (:) = nan
       allocate(this%livestemc_storage_to_xfer           (begp:endp)) ;    this%livestemc_storage_to_xfer            (:) = nan
       allocate(this%deadstemc_storage_to_xfer           (begp:endp)) ;    this%deadstemc_storage_to_xfer            (:) = nan
       allocate(this%livecrootc_storage_to_xfer          (begp:endp)) ;    this%livecrootc_storage_to_xfer           (:) = nan
       allocate(this%deadcrootc_storage_to_xfer          (begp:endp)) ;    this%deadcrootc_storage_to_xfer           (:) = nan
       allocate(this%gresp_storage_to_xfer               (begp:endp)) ;    this%gresp_storage_to_xfer                (:) = nan
       allocate(this%livestemc_to_deadstemc              (begp:endp)) ;    this%livestemc_to_deadstemc               (:) = nan
       allocate(this%livecrootc_to_deadcrootc            (begp:endp)) ;    this%livecrootc_to_deadcrootc             (:) = nan
       allocate(this%gpp                                 (begp:endp)) ;    this%gpp                                  (:) = nan
       allocate(this%gpp_before_downreg                  (begp:endp)) ;    this%gpp_before_downreg                   (:) = nan
       allocate(this%mr                                  (begp:endp)) ;    this%mr                                   (:) = nan
       allocate(this%current_gr                          (begp:endp)) ;    this%current_gr                           (:) = nan
       allocate(this%transfer_gr                         (begp:endp)) ;    this%transfer_gr                          (:) = nan
       allocate(this%storage_gr                          (begp:endp)) ;    this%storage_gr                           (:) = nan
       allocate(this%gr                                  (begp:endp)) ;    this%gr                                   (:) = nan
       allocate(this%ar                                  (begp:endp)) ;    this%ar                                   (:) = nan
       allocate(this%rr                                  (begp:endp)) ;    this%rr                                   (:) = nan
       allocate(this%npp                                 (begp:endp)) ;    this%npp                                  (:) = nan
       allocate(this%agnpp                               (begp:endp)) ;    this%agnpp                                (:) = nan
       allocate(this%bgnpp                               (begp:endp)) ;    this%bgnpp                                (:) = nan
       allocate(this%litfall                             (begp:endp)) ;    this%litfall                              (:) = nan
       allocate(this%vegfire                             (begp:endp)) ;    this%vegfire                              (:) = nan
       allocate(this%wood_harvestc                       (begp:endp)) ;    this%wood_harvestc                        (:) = nan
       allocate(this%cinputs                             (begp:endp)) ;    this%cinputs                              (:) = nan
       allocate(this%coutputs                            (begp:endp)) ;    this%coutputs                             (:) = nan
       allocate(this%plant_calloc                        (begp:endp)) ;    this%plant_calloc                         (:) = nan
       allocate(this%excess_cflux                        (begp:endp)) ;    this%excess_cflux                         (:) = nan
       allocate(this%prev_leafc_to_litter                (begp:endp)) ;    this%prev_leafc_to_litter                 (:) = nan
       allocate(this%prev_frootc_to_litter               (begp:endp)) ;    this%prev_frootc_to_litter                (:) = nan
       allocate(this%availc                              (begp:endp)) ;    this%availc                               (:) = nan
       allocate(this%xsmrpool_recover                    (begp:endp)) ;    this%xsmrpool_recover                     (:) = nan
       allocate(this%xsmrpool_c13ratio                   (begp:endp)) ;    this%xsmrpool_c13ratio                    (:) = nan
       allocate(this%xsmrpool_turnover                   (begp:endp)) ;    this%xsmrpool_turnover                    (:) = nan
       allocate(this%frootc_alloc                        (begp:endp)) ;    this%frootc_alloc                         (:) = nan
       allocate(this%frootc_loss                         (begp:endp)) ;    this%frootc_loss                          (:) = nan
       allocate(this%leafc_alloc                         (begp:endp)) ;    this%leafc_alloc                          (:) = nan
       allocate(this%leafc_loss                          (begp:endp)) ;    this%leafc_loss                           (:) = nan
       allocate(this%woodc_alloc                         (begp:endp)) ;    this%woodc_alloc                          (:) = nan
       allocate(this%woodc_loss                          (begp:endp)) ;    this%woodc_loss                           (:) = nan
       allocate(this%fire_closs                          (begp:endp)) ;    this%fire_closs                           (:) = nan
       allocate(this%crop_seedc_to_leaf                  (begp:endp)) ;    this%crop_seedc_to_leaf                   (:) = nan
    end if ! .not use fates

    allocate(this%dwt_seedc_to_leaf                   (begp:endp)) ;    this%dwt_seedc_to_leaf                    (:) = nan
    allocate(this%dwt_seedc_to_deadstem               (begp:endp)) ;    this%dwt_seedc_to_deadstem                (:) = nan
    allocate(this%dwt_conv_cflux                      (begp:endp)) ;    this%dwt_conv_cflux                       (:) = nan
    allocate(this%dwt_prod10c_gain                    (begp:endp)) ;    this%dwt_prod10c_gain                     (:) = nan
    allocate(this%dwt_prod100c_gain                   (begp:endp)) ;    this%dwt_prod100c_gain                    (:) = nan
    allocate(this%dwt_crop_productc_gain              (begp:endp)) ;    this%dwt_crop_productc_gain               (:) = nan
    allocate(this%tempsum_npp                         (begp:endp)) ;    this%tempsum_npp                          (:) = nan
    allocate(this%annsum_npp                          (begp:endp)) ;    this%annsum_npp                           (:) = nan
    allocate(this%annavg_agnpp                        (begp:endp)) ;    this%annavg_agnpp                         (:) = spval
    allocate(this%annavg_bgnpp                        (begp:endp)) ;    this%annavg_bgnpp                         (:) = spval
    allocate(this%tempavg_agnpp                       (begp:endp)) ;    this%tempavg_agnpp                        (:) = spval
    allocate(this%tempavg_bgnpp                       (begp:endp)) ;    this%tempavg_bgnpp                        (:) = spval
    allocate(this%agwdnpp                             (begp:endp)) ;    this%agwdnpp                              (:) = nan
    allocate(this%allocation_leaf                     (begp:endp)) ;    this%allocation_leaf                      (:) = nan
    allocate(this%allocation_stem                     (begp:endp)) ;    this%allocation_stem                      (:) = nan
    allocate(this%allocation_froot                    (begp:endp)) ;    this%allocation_froot                     (:) = nan

  end subroutine veg_cf_init
    
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine veg_nf_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    integer :: fp                                        ! filter indices
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_patch(endp-begp+1)  ! special landunit filter - patches
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_nf
    !-----------------------------------------------------------------------
    allocate(this%m_leafn_to_litter                   (begp:endp)) ; this%m_leafn_to_litter                   (:) = nan
    allocate(this%m_frootn_to_litter                  (begp:endp)) ; this%m_frootn_to_litter                  (:) = nan
    allocate(this%m_leafn_storage_to_litter           (begp:endp)) ; this%m_leafn_storage_to_litter           (:) = nan
    allocate(this%m_frootn_storage_to_litter          (begp:endp)) ; this%m_frootn_storage_to_litter          (:) = nan
    allocate(this%m_livestemn_storage_to_litter       (begp:endp)) ; this%m_livestemn_storage_to_litter       (:) = nan
    allocate(this%m_deadstemn_storage_to_litter       (begp:endp)) ; this%m_deadstemn_storage_to_litter       (:) = nan
    allocate(this%m_livecrootn_storage_to_litter      (begp:endp)) ; this%m_livecrootn_storage_to_litter      (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter      (begp:endp)) ; this%m_deadcrootn_storage_to_litter      (:) = nan
    allocate(this%m_leafn_xfer_to_litter              (begp:endp)) ; this%m_leafn_xfer_to_litter              (:) = nan
    allocate(this%m_frootn_xfer_to_litter             (begp:endp)) ; this%m_frootn_xfer_to_litter             (:) = nan
    allocate(this%m_livestemn_xfer_to_litter          (begp:endp)) ; this%m_livestemn_xfer_to_litter          (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter          (begp:endp)) ; this%m_deadstemn_xfer_to_litter          (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter         (begp:endp)) ; this%m_livecrootn_xfer_to_litter         (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter         (begp:endp)) ; this%m_deadcrootn_xfer_to_litter         (:) = nan
    allocate(this%m_livestemn_to_litter               (begp:endp)) ; this%m_livestemn_to_litter               (:) = nan
    allocate(this%m_deadstemn_to_litter               (begp:endp)) ; this%m_deadstemn_to_litter               (:) = nan
    allocate(this%m_livecrootn_to_litter              (begp:endp)) ; this%m_livecrootn_to_litter              (:) = nan
    allocate(this%m_deadcrootn_to_litter              (begp:endp)) ; this%m_deadcrootn_to_litter              (:) = nan
    allocate(this%m_retransn_to_litter                (begp:endp)) ; this%m_retransn_to_litter                (:) = nan
    allocate(this%m_npool_to_litter                   (begp:endp)) ; this%m_npool_to_litter                   (:) = nan
    allocate(this%hrv_leafn_to_litter                 (begp:endp)) ; this%hrv_leafn_to_litter                 (:) = nan
    allocate(this%hrv_frootn_to_litter                (begp:endp)) ; this%hrv_frootn_to_litter                (:) = nan
    allocate(this%hrv_leafn_storage_to_litter         (begp:endp)) ; this%hrv_leafn_storage_to_litter         (:) = nan
    allocate(this%hrv_frootn_storage_to_litter        (begp:endp)) ; this%hrv_frootn_storage_to_litter        (:) = nan
    allocate(this%hrv_livestemn_storage_to_litter     (begp:endp)) ; this%hrv_livestemn_storage_to_litter     (:) = nan
    allocate(this%hrv_deadstemn_storage_to_litter     (begp:endp)) ; this%hrv_deadstemn_storage_to_litter     (:) = nan
    allocate(this%hrv_livecrootn_storage_to_litter    (begp:endp)) ; this%hrv_livecrootn_storage_to_litter    (:) = nan
    allocate(this%hrv_deadcrootn_storage_to_litter    (begp:endp)) ; this%hrv_deadcrootn_storage_to_litter    (:) = nan
    allocate(this%hrv_leafn_xfer_to_litter            (begp:endp)) ; this%hrv_leafn_xfer_to_litter            (:) = nan
    allocate(this%hrv_frootn_xfer_to_litter           (begp:endp)) ; this%hrv_frootn_xfer_to_litter           (:) = nan
    allocate(this%hrv_livestemn_xfer_to_litter        (begp:endp)) ; this%hrv_livestemn_xfer_to_litter        (:) = nan
    allocate(this%hrv_deadstemn_xfer_to_litter        (begp:endp)) ; this%hrv_deadstemn_xfer_to_litter        (:) = nan
    allocate(this%hrv_livecrootn_xfer_to_litter       (begp:endp)) ; this%hrv_livecrootn_xfer_to_litter       (:) = nan
    allocate(this%hrv_deadcrootn_xfer_to_litter       (begp:endp)) ; this%hrv_deadcrootn_xfer_to_litter       (:) = nan
    allocate(this%hrv_livestemn_to_litter             (begp:endp)) ; this%hrv_livestemn_to_litter             (:) = nan
    allocate(this%hrv_deadstemn_to_prod10n            (begp:endp)) ; this%hrv_deadstemn_to_prod10n            (:) = nan
    allocate(this%hrv_deadstemn_to_prod100n           (begp:endp)) ; this%hrv_deadstemn_to_prod100n           (:) = nan
    allocate(this%hrv_leafn_to_prod1n                 (begp:endp)) ; this%hrv_leafn_to_prod1n                 (:) = nan
    allocate(this%hrv_livestemn_to_prod1n             (begp:endp)) ; this%hrv_livestemn_to_prod1n             (:) = nan
    allocate(this%hrv_grainn_to_prod1n                (begp:endp)) ; this%hrv_grainn_to_prod1n                (:) = nan
    allocate(this%hrv_cropn_to_prod1n                 (begp:endp)) ; this%hrv_cropn_to_prod1n                 (:) = nan
    allocate(this%hrv_livecrootn_to_litter            (begp:endp)) ; this%hrv_livecrootn_to_litter            (:) = nan
    allocate(this%hrv_deadcrootn_to_litter            (begp:endp)) ; this%hrv_deadcrootn_to_litter            (:) = nan
    allocate(this%hrv_retransn_to_litter              (begp:endp)) ; this%hrv_retransn_to_litter              (:) = nan
    allocate(this%hrv_npool_to_litter                 (begp:endp)) ; this%hrv_npool_to_litter                 (:) = nan
    allocate(this%m_leafn_to_fire                     (begp:endp)) ; this%m_leafn_to_fire                     (:) = nan
    allocate(this%m_leafn_storage_to_fire             (begp:endp)) ; this%m_leafn_storage_to_fire             (:) = nan
    allocate(this%m_leafn_xfer_to_fire                (begp:endp)) ; this%m_leafn_xfer_to_fire                (:) = nan
    allocate(this%m_livestemn_to_fire                 (begp:endp)) ; this%m_livestemn_to_fire                 (:) = nan
    allocate(this%m_livestemn_storage_to_fire         (begp:endp)) ; this%m_livestemn_storage_to_fire         (:) = nan
    allocate(this%m_livestemn_xfer_to_fire            (begp:endp)) ; this%m_livestemn_xfer_to_fire            (:) = nan
    allocate(this%m_deadstemn_to_fire                 (begp:endp)) ; this%m_deadstemn_to_fire                 (:) = nan
    allocate(this%m_deadstemn_storage_to_fire         (begp:endp)) ; this%m_deadstemn_storage_to_fire         (:) = nan
    allocate(this%m_deadstemn_xfer_to_fire            (begp:endp)) ; this%m_deadstemn_xfer_to_fire            (:) = nan
    allocate(this%m_frootn_to_fire                    (begp:endp)) ; this%m_frootn_to_fire                    (:) = nan
    allocate(this%m_frootn_storage_to_fire            (begp:endp)) ; this%m_frootn_storage_to_fire            (:) = nan
    allocate(this%m_frootn_xfer_to_fire               (begp:endp)) ; this%m_frootn_xfer_to_fire               (:) = nan
    allocate(this%m_livecrootn_to_fire                (begp:endp)) ; this%m_livecrootn_to_fire                (:) = nan     
    allocate(this%m_livecrootn_storage_to_fire        (begp:endp)) ; this%m_livecrootn_storage_to_fire        (:) = nan
    allocate(this%m_livecrootn_xfer_to_fire           (begp:endp)) ; this%m_livecrootn_xfer_to_fire           (:) = nan
    allocate(this%m_deadcrootn_to_fire                (begp:endp)) ; this%m_deadcrootn_to_fire                (:) = nan
    allocate(this%m_deadcrootn_storage_to_fire        (begp:endp)) ; this%m_deadcrootn_storage_to_fire        (:) = nan
    allocate(this%m_deadcrootn_xfer_to_fire           (begp:endp)) ; this%m_deadcrootn_xfer_to_fire           (:) = nan
    allocate(this%m_retransn_to_fire                  (begp:endp)) ; this%m_retransn_to_fire                  (:) = nan
    allocate(this%m_npool_to_fire                     (begp:endp)) ; this%m_npool_to_fire                     (:) = nan
    allocate(this%m_leafn_to_litter_fire              (begp:endp)) ; this%m_leafn_to_litter_fire              (:) = nan
    allocate(this%m_leafn_storage_to_litter_fire      (begp:endp)) ; this%m_leafn_storage_to_litter_fire      (:) = nan
    allocate(this%m_leafn_xfer_to_litter_fire         (begp:endp)) ; this%m_leafn_xfer_to_litter_fire         (:) = nan
    allocate(this%m_livestemn_to_litter_fire          (begp:endp)) ; this%m_livestemn_to_litter_fire          (:) = nan
    allocate(this%m_livestemn_storage_to_litter_fire  (begp:endp)) ; this%m_livestemn_storage_to_litter_fire  (:) = nan
    allocate(this%m_livestemn_xfer_to_litter_fire     (begp:endp)) ; this%m_livestemn_xfer_to_litter_fire     (:) = nan
    allocate(this%m_livestemn_to_deadstemn_fire       (begp:endp)) ; this%m_livestemn_to_deadstemn_fire       (:) = nan
    allocate(this%m_deadstemn_to_litter_fire          (begp:endp)) ; this%m_deadstemn_to_litter_fire          (:) = nan
    allocate(this%m_deadstemn_storage_to_litter_fire  (begp:endp)) ; this%m_deadstemn_storage_to_litter_fire  (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter_fire     (begp:endp)) ; this%m_deadstemn_xfer_to_litter_fire     (:) = nan
    allocate(this%m_frootn_to_litter_fire             (begp:endp)) ; this%m_frootn_to_litter_fire             (:) = nan
    allocate(this%m_frootn_storage_to_litter_fire     (begp:endp)) ; this%m_frootn_storage_to_litter_fire     (:) = nan
    allocate(this%m_frootn_xfer_to_litter_fire        (begp:endp)) ; this%m_frootn_xfer_to_litter_fire        (:) = nan
    allocate(this%m_livecrootn_to_litter_fire         (begp:endp)) ; this%m_livecrootn_to_litter_fire         (:) = nan
    allocate(this%m_livecrootn_storage_to_litter_fire (begp:endp)) ; this%m_livecrootn_storage_to_litter_fire (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter_fire    (begp:endp)) ; this%m_livecrootn_xfer_to_litter_fire    (:) = nan
    allocate(this%m_livecrootn_to_deadcrootn_fire     (begp:endp)) ; this%m_livecrootn_to_deadcrootn_fire     (:) = nan
    allocate(this%m_deadcrootn_to_litter_fire         (begp:endp)) ; this%m_deadcrootn_to_litter_fire         (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter_fire (begp:endp)) ; this%m_deadcrootn_storage_to_litter_fire (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter_fire    (begp:endp)) ; this%m_deadcrootn_xfer_to_litter_fire    (:) = nan
    allocate(this%m_retransn_to_litter_fire           (begp:endp)) ; this%m_retransn_to_litter_fire           (:) = nan
    allocate(this%m_npool_to_litter_fire              (begp:endp)) ; this%m_npool_to_litter_fire              (:) = nan
    allocate(this%leafn_xfer_to_leafn                 (begp:endp)) ; this%leafn_xfer_to_leafn                 (:) = nan
    allocate(this%frootn_xfer_to_frootn               (begp:endp)) ; this%frootn_xfer_to_frootn               (:) = nan
    allocate(this%livestemn_xfer_to_livestemn         (begp:endp)) ; this%livestemn_xfer_to_livestemn         (:) = nan
    allocate(this%deadstemn_xfer_to_deadstemn         (begp:endp)) ; this%deadstemn_xfer_to_deadstemn         (:) = nan
    allocate(this%livecrootn_xfer_to_livecrootn       (begp:endp)) ; this%livecrootn_xfer_to_livecrootn       (:) = nan
    allocate(this%deadcrootn_xfer_to_deadcrootn       (begp:endp)) ; this%deadcrootn_xfer_to_deadcrootn       (:) = nan
    allocate(this%leafn_to_litter                     (begp:endp)) ; this%leafn_to_litter                     (:) = nan
    allocate(this%leafn_to_retransn                   (begp:endp)) ; this%leafn_to_retransn                   (:) = nan
    allocate(this%frootn_to_retransn                  (begp:endp)) ; this%frootn_to_retransn                  (:) = nan
    allocate(this%frootn_to_litter                    (begp:endp)) ; this%frootn_to_litter                    (:) = nan
    allocate(this%retransn_to_npool                   (begp:endp)) ; this%retransn_to_npool                   (:) = nan
    allocate(this%sminn_to_npool                      (begp:endp)) ; this%sminn_to_npool                      (:) = nan
    allocate(this%npool_to_leafn                      (begp:endp)) ; this%npool_to_leafn                      (:) = nan
    allocate(this%npool_to_leafn_storage              (begp:endp)) ; this%npool_to_leafn_storage              (:) = nan
    allocate(this%npool_to_frootn                     (begp:endp)) ; this%npool_to_frootn                     (:) = nan
    allocate(this%npool_to_frootn_storage             (begp:endp)) ; this%npool_to_frootn_storage             (:) = nan
    allocate(this%npool_to_livestemn                  (begp:endp)) ; this%npool_to_livestemn                  (:) = nan
    allocate(this%npool_to_livestemn_storage          (begp:endp)) ; this%npool_to_livestemn_storage          (:) = nan
    allocate(this%npool_to_deadstemn                  (begp:endp)) ; this%npool_to_deadstemn                  (:) = nan
    allocate(this%npool_to_deadstemn_storage          (begp:endp)) ; this%npool_to_deadstemn_storage          (:) = nan
    allocate(this%npool_to_livecrootn                 (begp:endp)) ; this%npool_to_livecrootn                 (:) = nan
    allocate(this%npool_to_livecrootn_storage         (begp:endp)) ; this%npool_to_livecrootn_storage         (:) = nan
    allocate(this%npool_to_deadcrootn                 (begp:endp)) ; this%npool_to_deadcrootn                 (:) = nan
    allocate(this%npool_to_deadcrootn_storage         (begp:endp)) ; this%npool_to_deadcrootn_storage         (:) = nan
    allocate(this%leafn_storage_to_xfer               (begp:endp)) ; this%leafn_storage_to_xfer               (:) = nan
    allocate(this%frootn_storage_to_xfer              (begp:endp)) ; this%frootn_storage_to_xfer              (:) = nan
    allocate(this%livestemn_storage_to_xfer           (begp:endp)) ; this%livestemn_storage_to_xfer           (:) = nan
    allocate(this%deadstemn_storage_to_xfer           (begp:endp)) ; this%deadstemn_storage_to_xfer           (:) = nan
    allocate(this%livecrootn_storage_to_xfer          (begp:endp)) ; this%livecrootn_storage_to_xfer          (:) = nan
    allocate(this%deadcrootn_storage_to_xfer          (begp:endp)) ; this%deadcrootn_storage_to_xfer          (:) = nan
    allocate(this%livestemn_to_deadstemn              (begp:endp)) ; this%livestemn_to_deadstemn              (:) = nan
    allocate(this%livestemn_to_retransn               (begp:endp)) ; this%livestemn_to_retransn               (:) = nan
    allocate(this%livecrootn_to_deadcrootn            (begp:endp)) ; this%livecrootn_to_deadcrootn            (:) = nan
    allocate(this%livecrootn_to_retransn              (begp:endp)) ; this%livecrootn_to_retransn              (:) = nan
    allocate(this%ndeploy                             (begp:endp)) ; this%ndeploy                             (:) = nan
    allocate(this%ninputs                             (begp:endp)) ; this%ninputs                             (:) = nan
    allocate(this%noutputs                            (begp:endp)) ; this%noutputs                            (:) = nan
    allocate(this%wood_harvestn                       (begp:endp)) ; this%wood_harvestn                       (:) = nan
    allocate(this%fire_nloss                          (begp:endp)) ; this%fire_nloss                          (:) = nan
    allocate(this%npool_to_grainn                     (begp:endp)) ; this%npool_to_grainn                     (:) = nan
    allocate(this%npool_to_grainn_storage             (begp:endp)) ; this%npool_to_grainn_storage             (:) = nan
    allocate(this%livestemn_to_litter                 (begp:endp)) ; this%livestemn_to_litter                 (:) = nan
    allocate(this%grainn_to_food                      (begp:endp)) ; this%grainn_to_food                      (:) = nan
    allocate(this%grainn_xfer_to_grainn               (begp:endp)) ; this%grainn_xfer_to_grainn               (:) = nan
    allocate(this%grainn_storage_to_xfer              (begp:endp)) ; this%grainn_storage_to_xfer              (:) = nan
    allocate(this%fert                                (begp:endp)) ; this%fert                                (:) = nan
    allocate(this%fert_counter                        (begp:endp)) ; this%fert_counter                        (:) = nan
    allocate(this%soyfixn                             (begp:endp)) ; this%soyfixn                             (:) = nan
    allocate(this%nfix_to_plantn                      (begp:endp)) ; this%nfix_to_plantn                      (:) = nan
    allocate(this%crop_seedn_to_leaf                  (begp:endp)) ; this%crop_seedn_to_leaf                  (:) = nan
    allocate(this%dwt_seedn_to_leaf                   (begp:endp)) ; this%dwt_seedn_to_leaf                   (:) = nan
    allocate(this%dwt_seedn_to_deadstem               (begp:endp)) ; this%dwt_seedn_to_deadstem               (:) = nan
    allocate(this%dwt_conv_nflux                      (begp:endp)) ; this%dwt_conv_nflux                      (:) = nan
    allocate(this%dwt_prod10n_gain                    (begp:endp)) ; this%dwt_prod10n_gain                    (:) = nan
    allocate(this%dwt_prod100n_gain                   (begp:endp)) ; this%dwt_prod100n_gain                   (:) = nan
    allocate(this%dwt_crop_productn_gain              (begp:endp)) ; this%dwt_crop_productn_gain              (:) = nan
    allocate(this%dwt_seedn_to_npool                  (begp:endp)) ; this%dwt_seedn_to_npool                  (:) = nan
    allocate(this%plant_ndemand                       (begp:endp)) ; this%plant_ndemand                       (:) = nan
    allocate(this%avail_retransn                      (begp:endp)) ; this%avail_retransn                      (:) = nan
    allocate(this%plant_nalloc                        (begp:endp)) ; this%plant_nalloc                        (:) = nan
    allocate(this%smin_no3_to_plant                   (begp:endp)) ; this%smin_no3_to_plant                   (:) = nan
    allocate(this%smin_nh4_to_plant                   (begp:endp)) ; this%smin_nh4_to_plant                   (:) = nan
    allocate(this%sminn_to_plant                      (begp:endp)) ; this%sminn_to_plant                      (:) = nan
    allocate(this%plant_nh4demand_vr                  (begp:endp,1:nlevdecomp)); this%plant_nh4demand_vr    (:,:) = nan
    allocate(this%plant_no3demand_vr                  (begp:endp,1:nlevdecomp)); this%plant_no3demand_vr    (:,:) = nan
    allocate(this%plant_ndemand_vr                    (begp:endp,1:nlevdecomp)); this%plant_ndemand_vr      (:,:) = nan
    allocate(this%prev_leafn_to_litter                (begp:endp)) ; this%prev_leafn_to_litter                (:) = nan
    allocate(this%prev_frootn_to_litter               (begp:endp)) ; this%prev_frootn_to_litter               (:) = nan
    allocate(this%supplement_to_plantn                (begp:endp)) ; this%supplement_to_plantn                (:) = 0.d0
    allocate(this%gap_nloss_litter                    (begp:endp)) ; this%gap_nloss_litter                    (:) = nan
    allocate(this%fire_nloss_litter                   (begp:endp)) ; this%fire_nloss_litter                   (:) = nan
    allocate(this%hrv_nloss_litter                    (begp:endp)) ; this%hrv_nloss_litter                    (:) = nan
    allocate(this%sen_nloss_litter                    (begp:endp)) ; this%sen_nloss_litter                    (:) = nan
    
  end subroutine veg_nf_init

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine veg_pf_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer :: p,l                         ! indices
    integer :: fp                          ! filter indices
    integer :: num_special_patch           ! number of good values in special_patch filter
    integer :: special_patch(endp-begp+1)  ! special landunit filter - patches
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_pf
    !-----------------------------------------------------------------------
    allocate(this%m_leafp_to_litter                   (begp:endp)) ; this%m_leafp_to_litter                   (:) = nan
    allocate(this%m_frootp_to_litter                  (begp:endp)) ; this%m_frootp_to_litter                  (:) = nan
    allocate(this%m_leafp_storage_to_litter           (begp:endp)) ; this%m_leafp_storage_to_litter           (:) = nan
    allocate(this%m_frootp_storage_to_litter          (begp:endp)) ; this%m_frootp_storage_to_litter          (:) = nan
    allocate(this%m_livestemp_storage_to_litter       (begp:endp)) ; this%m_livestemp_storage_to_litter       (:) = nan
    allocate(this%m_deadstemp_storage_to_litter       (begp:endp)) ; this%m_deadstemp_storage_to_litter       (:) = nan
    allocate(this%m_livecrootp_storage_to_litter      (begp:endp)) ; this%m_livecrootp_storage_to_litter      (:) = nan
    allocate(this%m_deadcrootp_storage_to_litter      (begp:endp)) ; this%m_deadcrootp_storage_to_litter      (:) = nan
    allocate(this%m_leafp_xfer_to_litter              (begp:endp)) ; this%m_leafp_xfer_to_litter              (:) = nan
    allocate(this%m_frootp_xfer_to_litter             (begp:endp)) ; this%m_frootp_xfer_to_litter             (:) = nan
    allocate(this%m_livestemp_xfer_to_litter          (begp:endp)) ; this%m_livestemp_xfer_to_litter          (:) = nan
    allocate(this%m_deadstemp_xfer_to_litter          (begp:endp)) ; this%m_deadstemp_xfer_to_litter          (:) = nan
    allocate(this%m_livecrootp_xfer_to_litter         (begp:endp)) ; this%m_livecrootp_xfer_to_litter         (:) = nan
    allocate(this%m_deadcrootp_xfer_to_litter         (begp:endp)) ; this%m_deadcrootp_xfer_to_litter         (:) = nan
    allocate(this%m_livestemp_to_litter               (begp:endp)) ; this%m_livestemp_to_litter               (:) = nan
    allocate(this%m_deadstemp_to_litter               (begp:endp)) ; this%m_deadstemp_to_litter               (:) = nan
    allocate(this%m_livecrootp_to_litter              (begp:endp)) ; this%m_livecrootp_to_litter              (:) = nan
    allocate(this%m_deadcrootp_to_litter              (begp:endp)) ; this%m_deadcrootp_to_litter              (:) = nan
    allocate(this%m_retransp_to_litter                (begp:endp)) ; this%m_retransp_to_litter                (:) = nan
    allocate(this%m_ppool_to_litter                   (begp:endp)) ; this%m_ppool_to_litter                   (:) = nan
    allocate(this%hrv_leafp_to_litter                 (begp:endp)) ; this%hrv_leafp_to_litter                 (:) = nan
    allocate(this%hrv_frootp_to_litter                (begp:endp)) ; this%hrv_frootp_to_litter                (:) = nan
    allocate(this%hrv_leafp_storage_to_litter         (begp:endp)) ; this%hrv_leafp_storage_to_litter         (:) = nan
    allocate(this%hrv_frootp_storage_to_litter        (begp:endp)) ; this%hrv_frootp_storage_to_litter        (:) = nan
    allocate(this%hrv_livestemp_storage_to_litter     (begp:endp)) ; this%hrv_livestemp_storage_to_litter     (:) = nan
    allocate(this%hrv_deadstemp_storage_to_litter     (begp:endp)) ; this%hrv_deadstemp_storage_to_litter     (:) = nan
    allocate(this%hrv_livecrootp_storage_to_litter    (begp:endp)) ; this%hrv_livecrootp_storage_to_litter    (:) = nan
    allocate(this%hrv_deadcrootp_storage_to_litter    (begp:endp)) ; this%hrv_deadcrootp_storage_to_litter    (:) = nan
    allocate(this%hrv_leafp_xfer_to_litter            (begp:endp)) ; this%hrv_leafp_xfer_to_litter            (:) = nan
    allocate(this%hrv_frootp_xfer_to_litter           (begp:endp)) ; this%hrv_frootp_xfer_to_litter           (:) = nan
    allocate(this%hrv_livestemp_xfer_to_litter        (begp:endp)) ; this%hrv_livestemp_xfer_to_litter        (:) = nan
    allocate(this%hrv_deadstemp_xfer_to_litter        (begp:endp)) ; this%hrv_deadstemp_xfer_to_litter        (:) = nan
    allocate(this%hrv_livecrootp_xfer_to_litter       (begp:endp)) ; this%hrv_livecrootp_xfer_to_litter       (:) = nan
    allocate(this%hrv_deadcrootp_xfer_to_litter       (begp:endp)) ; this%hrv_deadcrootp_xfer_to_litter       (:) = nan
    allocate(this%hrv_livestemp_to_litter             (begp:endp)) ; this%hrv_livestemp_to_litter             (:) = nan
    allocate(this%hrv_deadstemp_to_prod10p            (begp:endp)) ; this%hrv_deadstemp_to_prod10p            (:) = nan
    allocate(this%hrv_deadstemp_to_prod100p           (begp:endp)) ; this%hrv_deadstemp_to_prod100p           (:) = nan
    allocate(this%hrv_leafp_to_prod1p                 (begp:endp)) ; this%hrv_leafp_to_prod1p                 (:) = nan
    allocate(this%hrv_livestemp_to_prod1p             (begp:endp)) ; this%hrv_livestemp_to_prod1p             (:) = nan
    allocate(this%hrv_grainp_to_prod1p                (begp:endp)) ; this%hrv_grainp_to_prod1p                (:) = nan
    allocate(this%hrv_cropp_to_prod1p                 (begp:endp)) ; this%hrv_cropp_to_prod1p                 (:) = nan
    allocate(this%hrv_livecrootp_to_litter            (begp:endp)) ; this%hrv_livecrootp_to_litter            (:) = nan
    allocate(this%hrv_deadcrootp_to_litter            (begp:endp)) ; this%hrv_deadcrootp_to_litter            (:) = nan
    allocate(this%hrv_retransp_to_litter              (begp:endp)) ; this%hrv_retransp_to_litter              (:) = nan
    allocate(this%hrv_ppool_to_litter                 (begp:endp)) ; this%hrv_ppool_to_litter                 (:) = nan
    allocate(this%m_leafp_to_fire                     (begp:endp)) ; this%m_leafp_to_fire                     (:) = nan
    allocate(this%m_leafp_storage_to_fire             (begp:endp)) ; this%m_leafp_storage_to_fire             (:) = nan
    allocate(this%m_leafp_xfer_to_fire                (begp:endp)) ; this%m_leafp_xfer_to_fire                (:) = nan
    allocate(this%m_livestemp_to_fire                 (begp:endp)) ; this%m_livestemp_to_fire                 (:) = nan
    allocate(this%m_livestemp_storage_to_fire         (begp:endp)) ; this%m_livestemp_storage_to_fire         (:) = nan
    allocate(this%m_livestemp_xfer_to_fire            (begp:endp)) ; this%m_livestemp_xfer_to_fire            (:) = nan
    allocate(this%m_deadstemp_to_fire                 (begp:endp)) ; this%m_deadstemp_to_fire                 (:) = nan
    allocate(this%m_deadstemp_storage_to_fire         (begp:endp)) ; this%m_deadstemp_storage_to_fire         (:) = nan
    allocate(this%m_deadstemp_xfer_to_fire            (begp:endp)) ; this%m_deadstemp_xfer_to_fire            (:) = nan
    allocate(this%m_frootp_to_fire                    (begp:endp)) ; this%m_frootp_to_fire                    (:) = nan
    allocate(this%m_frootp_storage_to_fire            (begp:endp)) ; this%m_frootp_storage_to_fire            (:) = nan
    allocate(this%m_frootp_xfer_to_fire               (begp:endp)) ; this%m_frootp_xfer_to_fire               (:) = nan
    allocate(this%m_livecrootp_to_fire                (begp:endp)) ; this%m_livecrootp_to_fire                (:) = nan    
    allocate(this%m_livecrootp_storage_to_fire        (begp:endp)) ; this%m_livecrootp_storage_to_fire        (:) = nan
    allocate(this%m_livecrootp_xfer_to_fire           (begp:endp)) ; this%m_livecrootp_xfer_to_fire           (:) = nan
    allocate(this%m_deadcrootp_to_fire                (begp:endp)) ; this%m_deadcrootp_to_fire                (:) = nan
    allocate(this%m_deadcrootp_storage_to_fire        (begp:endp)) ; this%m_deadcrootp_storage_to_fire        (:) = nan
    allocate(this%m_deadcrootp_xfer_to_fire           (begp:endp)) ; this%m_deadcrootp_xfer_to_fire           (:) = nan
    allocate(this%m_retransp_to_fire                  (begp:endp)) ; this%m_retransp_to_fire                  (:) = nan
    allocate(this%m_ppool_to_fire                     (begp:endp)) ; this%m_ppool_to_fire                     (:) = nan
    allocate(this%m_leafp_to_litter_fire              (begp:endp)) ; this%m_leafp_to_litter_fire              (:) = nan
    allocate(this%m_leafp_storage_to_litter_fire      (begp:endp)) ; this%m_leafp_storage_to_litter_fire      (:) = nan
    allocate(this%m_leafp_xfer_to_litter_fire         (begp:endp)) ; this%m_leafp_xfer_to_litter_fire         (:) = nan
    allocate(this%m_livestemp_to_litter_fire          (begp:endp)) ; this%m_livestemp_to_litter_fire          (:) = nan
    allocate(this%m_livestemp_storage_to_litter_fire  (begp:endp)) ; this%m_livestemp_storage_to_litter_fire  (:) = nan
    allocate(this%m_livestemp_xfer_to_litter_fire     (begp:endp)) ; this%m_livestemp_xfer_to_litter_fire     (:) = nan
    allocate(this%m_livestemp_to_deadstemp_fire       (begp:endp)) ; this%m_livestemp_to_deadstemp_fire       (:) = nan
    allocate(this%m_deadstemp_to_litter_fire          (begp:endp)) ; this%m_deadstemp_to_litter_fire          (:) = nan
    allocate(this%m_deadstemp_storage_to_litter_fire  (begp:endp)) ; this%m_deadstemp_storage_to_litter_fire  (:) = nan
    allocate(this%m_deadstemp_xfer_to_litter_fire     (begp:endp)) ; this%m_deadstemp_xfer_to_litter_fire     (:) = nan
    allocate(this%m_frootp_to_litter_fire             (begp:endp)) ; this%m_frootp_to_litter_fire             (:) = nan
    allocate(this%m_frootp_storage_to_litter_fire     (begp:endp)) ; this%m_frootp_storage_to_litter_fire     (:) = nan
    allocate(this%m_frootp_xfer_to_litter_fire        (begp:endp)) ; this%m_frootp_xfer_to_litter_fire        (:) = nan
    allocate(this%m_livecrootp_to_litter_fire         (begp:endp)) ; this%m_livecrootp_to_litter_fire         (:) = nan
    allocate(this%m_livecrootp_storage_to_litter_fire (begp:endp)) ; this%m_livecrootp_storage_to_litter_fire (:) = nan
    allocate(this%m_livecrootp_xfer_to_litter_fire    (begp:endp)) ; this%m_livecrootp_xfer_to_litter_fire    (:) = nan
    allocate(this%m_livecrootp_to_deadcrootp_fire     (begp:endp)) ; this%m_livecrootp_to_deadcrootp_fire     (:) = nan
    allocate(this%m_deadcrootp_to_litter_fire         (begp:endp)) ; this%m_deadcrootp_to_litter_fire         (:) = nan
    allocate(this%m_deadcrootp_storage_to_litter_fire (begp:endp)) ; this%m_deadcrootp_storage_to_litter_fire (:) = nan
    allocate(this%m_deadcrootp_xfer_to_litter_fire    (begp:endp)) ; this%m_deadcrootp_xfer_to_litter_fire    (:) = nan
    allocate(this%m_retransp_to_litter_fire           (begp:endp)) ; this%m_retransp_to_litter_fire           (:) = nan
    allocate(this%m_ppool_to_litter_fire              (begp:endp)) ; this%m_ppool_to_litter_fire              (:) = nan
    allocate(this%leafp_xfer_to_leafp                 (begp:endp)) ; this%leafp_xfer_to_leafp                 (:) = nan
    allocate(this%frootp_xfer_to_frootp               (begp:endp)) ; this%frootp_xfer_to_frootp               (:) = nan
    allocate(this%livestemp_xfer_to_livestemp         (begp:endp)) ; this%livestemp_xfer_to_livestemp         (:) = nan
    allocate(this%deadstemp_xfer_to_deadstemp         (begp:endp)) ; this%deadstemp_xfer_to_deadstemp         (:) = nan
    allocate(this%livecrootp_xfer_to_livecrootp       (begp:endp)) ; this%livecrootp_xfer_to_livecrootp       (:) = nan
    allocate(this%deadcrootp_xfer_to_deadcrootp       (begp:endp)) ; this%deadcrootp_xfer_to_deadcrootp       (:) = nan
    allocate(this%leafp_to_litter                     (begp:endp)) ; this%leafp_to_litter                     (:) = nan
    allocate(this%leafp_to_retransp                   (begp:endp)) ; this%leafp_to_retransp                   (:) = nan
    allocate(this%frootp_to_retransp                  (begp:endp)) ; this%frootp_to_retransp                  (:) = nan
    allocate(this%frootp_to_litter                    (begp:endp)) ; this%frootp_to_litter                    (:) = nan
    allocate(this%retransp_to_ppool                   (begp:endp)) ; this%retransp_to_ppool                   (:) = nan
    allocate(this%sminp_to_ppool                      (begp:endp)) ; this%sminp_to_ppool                      (:) = nan
    allocate(this%biochem_pmin_to_plant               (begp:endp)) ; this%biochem_pmin_to_plant               (:) = nan
    allocate(this%ppool_to_leafp                      (begp:endp)) ; this%ppool_to_leafp                      (:) = nan
    allocate(this%ppool_to_leafp_storage              (begp:endp)) ; this%ppool_to_leafp_storage              (:) = nan
    allocate(this%ppool_to_frootp                     (begp:endp)) ; this%ppool_to_frootp                     (:) = nan
    allocate(this%ppool_to_frootp_storage             (begp:endp)) ; this%ppool_to_frootp_storage             (:) = nan
    allocate(this%ppool_to_livestemp                  (begp:endp)) ; this%ppool_to_livestemp                  (:) = nan
    allocate(this%ppool_to_livestemp_storage          (begp:endp)) ; this%ppool_to_livestemp_storage          (:) = nan
    allocate(this%ppool_to_deadstemp                  (begp:endp)) ; this%ppool_to_deadstemp                  (:) = nan
    allocate(this%ppool_to_deadstemp_storage          (begp:endp)) ; this%ppool_to_deadstemp_storage          (:) = nan
    allocate(this%ppool_to_livecrootp                 (begp:endp)) ; this%ppool_to_livecrootp                 (:) = nan
    allocate(this%ppool_to_livecrootp_storage         (begp:endp)) ; this%ppool_to_livecrootp_storage         (:) = nan
    allocate(this%ppool_to_deadcrootp                 (begp:endp)) ; this%ppool_to_deadcrootp                 (:) = nan
    allocate(this%ppool_to_deadcrootp_storage         (begp:endp)) ; this%ppool_to_deadcrootp_storage         (:) = nan
    allocate(this%leafp_storage_to_xfer               (begp:endp)) ; this%leafp_storage_to_xfer               (:) = nan
    allocate(this%frootp_storage_to_xfer              (begp:endp)) ; this%frootp_storage_to_xfer              (:) = nan
    allocate(this%livestemp_storage_to_xfer           (begp:endp)) ; this%livestemp_storage_to_xfer           (:) = nan
    allocate(this%deadstemp_storage_to_xfer           (begp:endp)) ; this%deadstemp_storage_to_xfer           (:) = nan
    allocate(this%livecrootp_storage_to_xfer          (begp:endp)) ; this%livecrootp_storage_to_xfer          (:) = nan
    allocate(this%deadcrootp_storage_to_xfer          (begp:endp)) ; this%deadcrootp_storage_to_xfer          (:) = nan
    allocate(this%livestemp_to_deadstemp              (begp:endp)) ; this%livestemp_to_deadstemp              (:) = nan
    allocate(this%livestemp_to_retransp               (begp:endp)) ; this%livestemp_to_retransp               (:) = nan
    allocate(this%livecrootp_to_deadcrootp            (begp:endp)) ; this%livecrootp_to_deadcrootp            (:) = nan
    allocate(this%livecrootp_to_retransp              (begp:endp)) ; this%livecrootp_to_retransp              (:) = nan
    allocate(this%pdeploy                             (begp:endp)) ; this%pdeploy                             (:) = nan
    allocate(this%pinputs                             (begp:endp)) ; this%pinputs                             (:) = nan
    allocate(this%poutputs                            (begp:endp)) ; this%poutputs                            (:) = nan
    allocate(this%wood_harvestp                       (begp:endp)) ; this%wood_harvestp                       (:) = nan
    allocate(this%fire_ploss                          (begp:endp)) ; this%fire_ploss                          (:) = nan
    allocate(this%ppool_to_grainp                     (begp:endp)) ; this%ppool_to_grainp                     (:) = nan
    allocate(this%ppool_to_grainp_storage             (begp:endp)) ; this%ppool_to_grainp_storage             (:) = nan
    allocate(this%livestemp_to_litter                 (begp:endp)) ; this%livestemp_to_litter                 (:) = nan
    allocate(this%grainp_to_food                      (begp:endp)) ; this%grainp_to_food                      (:) = nan
    allocate(this%grainp_xfer_to_grainp               (begp:endp)) ; this%grainp_xfer_to_grainp               (:) = nan
    allocate(this%grainp_storage_to_xfer              (begp:endp)) ; this%grainp_storage_to_xfer              (:) = nan
    allocate(this%fert_p                              (begp:endp)) ; this%fert_p                              (:) = nan
    allocate(this%fert_p_counter                      (begp:endp)) ; this%fert_p_counter                      (:) = nan
    allocate(this%crop_seedp_to_leaf                  (begp:endp)) ; this%crop_seedp_to_leaf                  (:) = nan
    allocate(this%dwt_seedp_to_leaf                   (begp:endp)) ; this%dwt_seedp_to_leaf                   (:) = nan
    allocate(this%dwt_seedp_to_deadstem               (begp:endp)) ; this%dwt_seedp_to_deadstem               (:) = nan
    allocate(this%dwt_conv_pflux                      (begp:endp)) ; this%dwt_conv_pflux                      (:) = nan
    allocate(this%dwt_prod10p_gain                    (begp:endp)) ; this%dwt_prod10p_gain                    (:) = nan
    allocate(this%dwt_prod100p_gain                   (begp:endp)) ; this%dwt_prod100p_gain                   (:) = nan
    allocate(this%dwt_crop_productp_gain              (begp:endp)) ; this%dwt_crop_productp_gain              (:) = nan
    allocate(this%dwt_seedp_to_ppool                  (begp:endp)) ; this%dwt_seedp_to_ppool                  (:) = nan
    allocate(this%plant_pdemand                       (begp:endp)) ; this%plant_pdemand                       (:) = nan
    allocate(this%avail_retransp                      (begp:endp)) ; this%avail_retransp                      (:) = nan
    allocate(this%plant_palloc                        (begp:endp)) ; this%plant_palloc                        (:) = nan
    allocate(this%sminp_to_plant                      (begp:endp)) ; this%sminp_to_plant                      (:) = nan
    allocate(this%plant_pdemand_vr                    (begp:endp,1:nlevdecomp_full )) ; this%plant_pdemand_vr (:,:) = nan
    allocate(this%prev_leafp_to_litter                (begp:endp)) ; this%prev_leafp_to_litter                (:) = nan
    allocate(this%prev_frootp_to_litter               (begp:endp)) ; this%prev_frootp_to_litter               (:) = nan
    allocate(this%supplement_to_plantp                (begp:endp)) ; this%supplement_to_plantp                (:) = 0.d0
    allocate(this%gap_ploss_litter                    (begp:endp)) ; this%gap_ploss_litter                    (:) = nan
    allocate(this%fire_ploss_litter                   (begp:endp)) ; this%fire_ploss_litter                   (:) = nan
    allocate(this%hrv_ploss_litter                    (begp:endp)) ; this%hrv_ploss_litter                    (:) = nan
    allocate(this%sen_ploss_litter                    (begp:endp)) ; this%sen_ploss_litter                    (:) = nan

  end subroutine veg_pf_init

end module vegetationDataType
