module VegetationDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Vegetation data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=),isnan => shr_infnan_isnan
  use shr_const_mod   , only : SHR_CONST_PDB
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use spmdMod         , only : masterproc
  use abortutils      , only : endrun
  use clm_time_manager, only : is_restart, get_nstep
  use elm_varpar      , only : nlevsno, nlevgrnd, nlevlak, nlevurb, nlevcan, crop_prog
  use elm_varpar      , only : nlevdecomp, nlevdecomp_full
  use elm_varcon      , only : spval, ispval, sb
  use elm_varcon      , only : c13ratio, c14ratio
  use landunit_varcon , only : istsoil, istcrop
  use pftvarcon       , only : npcropmin, noveg, nstor
  use elm_varctl      , only : iulog, use_cn, spinup_state, spinup_mortality_factor, use_fates  
  use elm_varctl      , only : nu_com, use_crop, use_c13
  use elm_varctl      , only : use_lch4, use_betr
  use histFileMod     , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ncdio_pio       , only : file_desc_t, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
  use decompMod       , only : bounds_type, get_proc_global
  use subgridAveMod   , only : p2c
  use restUtilMod
  use CNStateType     , only: cnstate_type
  use SpeciesMod              , only : species_from_string
  use VegetationType            , only : veg_pp
  use VegetationPropertiesType  , only : veg_vp
  use LandunitType              , only : lun_pp
  use GridcellType              , only : grc_pp
  use ColumnDataType            , only : col_es
  use ColumnDataType            , only : column_carbon_state, column_carbon_flux
  use ColumnDataType            , only : column_nitrogen_state, column_nitrogen_flux
  use ColumnDataType            , only : column_phosphorus_state, column_phosphorus_flux
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
    procedure, public :: Restart => veg_es_restart
    procedure, public :: Clean   => veg_es_clean
    procedure, public :: InitAccBuffer => init_acc_buffer_veg_es
    procedure, public :: InitAccVars   => init_acc_vars_veg_es
    procedure, public :: UpdateAccVars => update_acc_vars_veg_es
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
    procedure, public :: Restart => veg_ws_restart
    procedure, public :: Clean   => veg_ws_clean
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
    procedure, public :: Restart  => veg_cs_restart
    procedure, public :: Summary  => veg_cs_summary
    procedure, public :: ZeroDwt  => veg_cs_zerodwt
    procedure, public :: Clean    => veg_cs_clean
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
    procedure, public :: Restart   => veg_ns_restart
    procedure, public :: Summary   => veg_ns_summary
    procedure, public :: SetValues => veg_ns_setvalues
    procedure, public :: ZeroDwt   => veg_ns_zerodwt
    procedure, public :: Clean     => veg_ns_clean
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
    procedure, public :: Restart   => veg_ps_restart
    procedure, public :: SetValues => veg_ps_setvalues
    procedure, public :: ZeroDWT   => veg_ps_zerodwt
    procedure, public :: Summary   => veg_ps_summary
    procedure, public :: Clean     => veg_ps_clean
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
    procedure, public :: Restart => veg_ef_restart
    procedure, public :: Clean   => veg_ef_clean
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
    real(r8), pointer :: qflx_dirct_rain    (:)   => null() ! direct through rainfall [mm/s]
    real(r8), pointer :: qflx_leafdrip      (:)   => null() ! leaf rain drip [mm/s]
    real(r8), pointer :: qflx_ev_snow       (:)   => null() ! evaporation heat flux from snow       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_soil       (:)   => null() ! evaporation heat flux from soil       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_h2osfc     (:)   => null() ! evaporation heat flux from soil       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_rootsoi_frac  (:,:) => null() !  
   
    real(r8), pointer :: irrig_rate               (:)   => null() ! current irrigation rate [mm/s]
    real(r8), pointer :: qflx_irrig_patch         (:)   => null()   ! patch irrigation flux (mm H2O/s)
    real(r8), pointer :: qflx_real_irrig_patch    (:)   => null()   ! patch real irrigation flux (mm H2O/s) 
    real(r8), pointer :: qflx_grnd_irrig_patch    (:)   => null()   ! groundwater irrigation (mm H2O/s) 
    real(r8), pointer :: qflx_surf_irrig_patch    (:)   => null()   ! surface water irrigation(mm H2O/s) 
    real(r8), pointer :: qflx_supply_patch        (:)   => null()   ! patch supply flux (mm H2O/s) 

    real(r8), pointer :: qflx_over_supply_patch   (:)   => null()   ! over supplied irrigation
    integer , pointer :: n_irrig_steps_left (:)   => null() ! number of time steps for which we still need to irrigate today (if 0, ignore)

  contains
    procedure, public :: Init    => veg_wf_init
    procedure, public :: Restart => veg_wf_restart
    procedure, public :: Clean   => veg_wf_clean
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
    procedure, public :: Restart    => veg_cf_restart
    procedure, public :: Summary    => veg_cf_summary
    procedure, public :: SummaryRR  => veg_cf_summary_rr      ! Root respiration summary
    procedure, public :: SummaryCH4 => veg_cf_summary_for_ch4 ! Summary for CH4 model
    procedure, public :: SetValues  => veg_cf_setvalues
    procedure, public :: Clean      => veg_cf_clean
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
    procedure, public :: Restart   => veg_nf_restart
    procedure, public :: SetValues => veg_nf_setvalues
    procedure, public :: Summary   => veg_nf_summary
    procedure, public :: Clean     => veg_nf_clean
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
    procedure, public :: Restart   => veg_pf_restart
    procedure, public :: SetValues => veg_pf_setvalues
    procedure, public :: Summary   => veg_pf_summary
    procedure, public :: Clean     => veg_pf_clean
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
    use elm_varctl     , only : use_vancouver, use_mexicocity
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

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_es
    !-----------------------------------------------------------------------
    this%t_veg(begp:endp) = spval
    call hist_addfld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_patch=this%t_veg)

    this%t_ref2m(begp:endp) = spval
    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_patch=this%t_ref2m)

    this%t_ref2m_r(begp:endp) = spval
    call hist_addfld1d (fname='TSA_R', units='K',  &
         avgflag='A', long_name='Rural 2m air temperature', &
         ptr_patch=this%t_ref2m_r, set_spec=spval)

    this%t_ref2m_u(begp:endp) = spval
    call hist_addfld1d (fname='TSA_U', units='K',  &
         avgflag='A', long_name='Urban 2m air temperature', &
         ptr_patch=this%t_ref2m_u, set_nourb=spval)

    this%t_a10(begp:endp) = spval
    call hist_addfld1d (fname='T10', units='K',  &
         avgflag='A', long_name='10-day running mean of 2-m temperature', &
         ptr_patch=this%t_a10, default='inactive')

    if (use_cn .and. crop_prog )then
       this%t_a10min(begp:endp) = spval
       call hist_addfld1d (fname='A10TMIN', units='K',  &
            avgflag='A', long_name='10-day running mean of min 2-m temperature', &
            ptr_patch=this%t_a10min, default='inactive')
    end if

    if (use_cn .and.  crop_prog )then
       this%t_a5min(begp:endp) = spval
       call hist_addfld1d (fname='A5TMIN', units='K',  &
            avgflag='A', long_name='5-day running mean of min 2-m temperature', &
            ptr_patch=this%t_a5min, default='inactive')
    end if

    this%t_ref2m_min(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min)

    this%t_ref2m_max(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max)

    this%t_ref2m_min_r(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV_R', units='K',  &
         avgflag='A', long_name='Rural daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_r, set_spec=spval)

    this%t_ref2m_max_r(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV_R', units='K',  &
         avgflag='A', long_name='Rural daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_r, set_spec=spval)

    this%t_ref2m_min_u(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV_U', units='K',  &
         avgflag='A', long_name='Urban daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_u, set_nourb=spval)

    this%t_ref2m_max_u(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV_U', units='K',  &
         avgflag='A', long_name='Urban daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_u, set_nourb=spval)

    this%t_veg24(begp:endp) = spval
    call hist_addfld1d (fname='TV24', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 24hrs)', &
         ptr_patch=this%t_veg24, default='inactive')

    this%t_veg240(begp:endp)  = spval
    call hist_addfld1d (fname='TV240', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 240hrs)', &
         ptr_patch=this%t_veg240, default='inactive')

    if (crop_prog) then
       this%gdd0(begp:endp) = spval
       call hist_addfld1d (fname='GDD0', units='ddays', &
            avgflag='A', long_name='Growing degree days base  0C from planting', &
            ptr_patch=this%gdd0, default='inactive')

       this%gdd8(begp:endp) = spval
       call hist_addfld1d (fname='GDD8', units='ddays', &
            avgflag='A', long_name='Growing degree days base  8C from planting', &
            ptr_patch=this%gdd8, default='inactive')

       this%gdd10(begp:endp) = spval
       call hist_addfld1d (fname='GDD10', units='ddays', &
            avgflag='A', long_name='Growing degree days base 10C from planting', &
            ptr_patch=this%gdd10, default='inactive')

       this%gdd020(begp:endp) = spval
       call hist_addfld1d (fname='GDD020', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base  0C from planting', &
            ptr_patch=this%gdd020, default='inactive')

       this%gdd820(begp:endp) = spval
       call hist_addfld1d (fname='GDD820', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base  8C from planting', &
            ptr_patch=this%gdd820, default='inactive')

       this%gdd1020(begp:endp) = spval
       call hist_addfld1d (fname='GDD1020', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base 10C from planting', &
            ptr_patch=this%gdd1020, default='inactive')
    end if
    
    if (use_cn ) then
       this%emv(begp:endp) = spval
       call hist_addfld1d (fname='EMV', units='proportion', &
            avgflag='A', long_name='vegetation emissivity', &
            ptr_patch=this%emv, default='inactive')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_es
    !-----------------------------------------------------------------------
    
    ! Set t_veg, t_ref2m, t_ref2m_u and tref2m_r 

    do p = begp, endp
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)

       if (use_vancouver) then
          this%t_veg(p)   = 297.56
       else if (use_mexicocity) then
          this%t_veg(p)   = 289.46
       else
          this%t_veg(p)   = 283._r8
       end if

       if (use_vancouver) then
          this%t_ref2m(p) = 297.56
       else if (use_mexicocity) then
          this%t_ref2m(p) = 289.46
       else
          this%t_ref2m(p) = 283._r8
       end if

       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%t_ref2m_u(p) = 297.56
          else if (use_mexicocity) then
             this%t_ref2m_u(p) = 289.46
          else
             this%t_ref2m_u(p) = 283._r8
          end if
       else 
          if (.not. lun_pp%ifspecial(l)) then 
             if (use_vancouver) then
                this%t_ref2m_r(p) = 297.56
             else if (use_mexicocity) then
                this%t_ref2m_r(p) = 289.46
             else
                this%t_ref2m_r(p) = 283._r8
             end if
          else 
             this%t_ref2m_r(p) = spval
          end if
       end if

    end do ! veg loop

  end subroutine veg_es_init
    
  !------------------------------------------------------------------------
  subroutine veg_es_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='T_VEG', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='vegetation temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_veg)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Rural 2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_U', xtype=ncd_double, dim1name='pft',                      &
         long_name='Urban 2m height surface air temperature', units='K',                                              &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_U', xtype=ncd_double, dim1name='pft',                  &
         long_name='urban daily minimum of average 2 m height surface air temperature (K)', units='K',                &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_U', xtype=ncd_double, dim1name='pft',                  &
         long_name='urban daily maximum of average 2 m height surface air temperature (K)', units='K',                &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_U', xtype=ncd_double, dim1name='pft',             &
         long_name='urban instantaneous daily min of average 2 m height surface air temp (K)', units='K',             &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_U', xtype=ncd_double,  dim1name='pft',            &
         long_name='urban instantaneous daily max of average 2 m height surface air temp (K)', units='K',             &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_u)

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='gdd1020', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 10C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd1020)

       call restartvar(ncid=ncid, flag=flag,  varname='gdd820', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 8C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd820)

       call restartvar(ncid=ncid, flag=flag,  varname='gdd020', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 0C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd020)
    end if


  end subroutine veg_es_restart

  !------------------------------------------------------------------------
  subroutine veg_es_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    !------------------------------------------------------------------------
  end subroutine veg_es_clean
  
  !-----------------------------------------------------------------------
  subroutine init_acc_buffer_veg_es (this, bounds)
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for accumulated fields for vegetation energy state
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use accumulMod       , only : init_accum_field
    use clm_time_manager , only : get_step_size
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES: 
    real(r8) :: dtime
    integer, parameter :: not_used = huge(1)
    !---------------------------------------------------------------------

    dtime = get_step_size()

    ! The following is a running mean. The accumulation period is set to -10 for a 10-day running mean.
    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20._r8)

    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='TREFAV_U', units='K', &
         desc='average over an hour of urban 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='TREFAV_R', units='K', &
         desc='average over an hour of rural 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%t_veg24(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='T_VEG24', units='K',                                              &
         desc='24hr average of vegetation temperature',  accum_type='runmean', accum_period=-1,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%t_veg240(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='T_VEG240', units='K',                                             &
         desc='240hr average of vegetation temperature',  accum_type='runmean', accum_period=-10,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    if ( crop_prog )then
       call init_accum_field (name='TDM10', units='K', &
            desc='10-day running mean of min 2-m temperature', accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       call init_accum_field (name='TDM5', units='K', &
            desc='5-day running mean of min 2-m temperature', accum_type='runmean', accum_period=-5, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       ! All GDD summations are relative to the planting date (Kucharik & Brye 2003)
       call init_accum_field (name='GDD0', units='K', &
            desc='growing degree-days base 0C from planting', accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD8', units='K', &
            desc='growing degree-days base 8C from planting', accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD10', units='K', &
            desc='growing degree-days base 10C from planting', accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)
    end if
    
  end subroutine init_acc_buffer_veg_es

  !-----------------------------------------------------------------------
  subroutine init_acc_vars_veg_es(this, bounds)
    ! !DESCRIPTION:
    ! Initialize variables associated with vegetation energy state
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    use elm_varctl       , only : nsrest, nsrStartup
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level pft field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="extract_accum_hist allocation error for rbufslp"//&
            errMsg(__FILE__, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('T10', rbufslp, nstep)
    this%t_a10(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('T_VEG24', rbufslp, nstep)
    this%t_veg24(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('T_VEG240', rbufslp, nstep)
    this%t_veg240(begp:endp) = rbufslp(begp:endp)

    if (crop_prog) then
       call extract_accum_field ('TDM10', rbufslp, nstep) 
       this%t_a10min(begp:endp)= rbufslp(begp:endp)

       call extract_accum_field ('TDM5', rbufslp, nstep) 
       this%t_a5min(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD0', rbufslp, nstep)
       this%gdd0(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD8', rbufslp, nstep) ;
       this%gdd8(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD10', rbufslp, nstep) 
       this%gdd10(begp:endp) = rbufslp(begp:endp)

    end if

    ! Initialize variables that are to be time accumulated
    ! Initialize 2m ref temperature max and min values

    if (nsrest == nsrStartup) then 
       this%t_ref2m_max(begp:endp)        =  spval
       this%t_ref2m_max_r(begp:endp)      =  spval
       this%t_ref2m_max_u(begp:endp)      =  spval

       this%t_ref2m_min(begp:endp)        =  spval
       this%t_ref2m_min_r(begp:endp)      =  spval
       this%t_ref2m_min_u(begp:endp)      =  spval

       this%t_ref2m_max_inst(begp:endp)   = -spval
       this%t_ref2m_max_inst_r(begp:endp) = -spval
       this%t_ref2m_max_inst_u(begp:endp) = -spval

       this%t_ref2m_min_inst(begp:endp)   =  spval
       this%t_ref2m_min_inst_r(begp:endp) =  spval
       this%t_ref2m_min_inst_u(begp:endp) =  spval
    end if

    deallocate(rbufslp)
    
  end subroutine init_acc_vars_veg_es

  !-----------------------------------------------------------------------
  subroutine update_acc_vars_veg_es (this, bounds)
    !
    ! USES
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod       , only : update_accum_field, extract_accum_field, accumResetVal
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state)    :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: m,g,l,c,p                 ! indices
    integer :: ier                       ! error status
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    logical :: end_cd                    ! temporary for is_end_curr_day() value
    integer :: begp, endp
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    
    ! fill the temporary variable
    do p = begp,endp
       rbufslp(p) = this%t_veg(p)
    end do
    
    call update_accum_field  ('T10', this%t_ref2m, nstep)
    call extract_accum_field ('T10', this%t_a10  , nstep)
    call update_accum_field  ('T_VEG24' , rbufslp       , nstep)
    call extract_accum_field ('T_VEG24' , this%t_veg24  , nstep)
    call update_accum_field  ('T_VEG240', rbufslp       , nstep)
    call extract_accum_field ('T_VEG240', this%t_veg240 , nstep)


    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', this%t_ref2m, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)
    end_cd = is_end_curr_day()
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst(p) = max(rbufslp(p), this%t_ref2m_max_inst(p))
          this%t_ref2m_min_inst(p) = min(rbufslp(p), this%t_ref2m_min_inst(p))
       endif
       if (end_cd) then
          this%t_ref2m_max(p) = this%t_ref2m_max_inst(p)
          this%t_ref2m_min(p) = this%t_ref2m_min_inst(p)
          this%t_ref2m_max_inst(p) = -spval
          this%t_ref2m_min_inst(p) =  spval
       else if (secs == int(dtime)) then
          this%t_ref2m_max(p) = spval
          this%t_ref2m_min(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_U - hourly average urban 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_U', this%t_ref2m_u, nstep)
    call extract_accum_field ('TREFAV_U', rbufslp, nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst_u(p) = max(rbufslp(p), this%t_ref2m_max_inst_u(p))
          this%t_ref2m_min_inst_u(p) = min(rbufslp(p), this%t_ref2m_min_inst_u(p))
       endif
       if (end_cd) then
         if (lun_pp%urbpoi(l)) then
          this%t_ref2m_max_u(p) = this%t_ref2m_max_inst_u(p)
          this%t_ref2m_min_u(p) = this%t_ref2m_min_inst_u(p)
          this%t_ref2m_max_inst_u(p) = -spval
          this%t_ref2m_min_inst_u(p) =  spval
         end if
       else if (secs == int(dtime)) then
          this%t_ref2m_max_u(p) = spval
          this%t_ref2m_min_u(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_R - hourly average rural 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_R', this%t_ref2m_r, nstep)
    call extract_accum_field ('TREFAV_R', rbufslp, nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst_r(p) = max(rbufslp(p), this%t_ref2m_max_inst_r(p))
          this%t_ref2m_min_inst_r(p) = min(rbufslp(p), this%t_ref2m_min_inst_r(p))
       endif
       if (end_cd) then
         if (.not.(lun_pp%ifspecial(l))) then
          this%t_ref2m_max_r(p) = this%t_ref2m_max_inst_r(p)
          this%t_ref2m_min_r(p) = this%t_ref2m_min_inst_r(p)
          this%t_ref2m_max_inst_r(p) = -spval
          this%t_ref2m_min_inst_r(p) =  spval
         end if
       else if (secs == int(dtime)) then
          this%t_ref2m_max_r(p) = spval
          this%t_ref2m_min_r(p) = spval
       endif
    end do

    if ( crop_prog )then
       ! Accumulate and extract TDM10

       do p = begp,endp
          rbufslp(p) = min(this%t_ref2m_min(p),this%t_ref2m_min_inst(p)) 
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                                     !'min_inst' not initialized?
       call update_accum_field  ('TDM10', rbufslp, nstep)
       call extract_accum_field ('TDM10', this%t_a10min, nstep)

       ! Accumulate and extract TDM5

       do p = begp,endp
          rbufslp(p) = min(this%t_ref2m_min(p),this%t_ref2m_min_inst(p)) 
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                         !'min_inst' not initialized?
       call update_accum_field  ('TDM5', rbufslp, nstep)
       call extract_accum_field ('TDM5', this%t_a5min, nstep)

       ! Accumulate and extract GDD0

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(26._r8, this%t_ref2m(p)-SHR_CONST_TKFRZ)) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD0', rbufslp, nstep)
       call extract_accum_field ('GDD0', this%gdd0, nstep)

       ! Accumulate and extract GDD8

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  this%t_ref2m(p)-(SHR_CONST_TKFRZ + 8._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD8', rbufslp, nstep)
       call extract_accum_field ('GDD8', this%gdd8, nstep)

       ! Accumulate and extract GDD10

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  this%t_ref2m(p)-(SHR_CONST_TKFRZ + 10._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD10', rbufslp, nstep)
       call extract_accum_field ('GDD10', this%gdd10, nstep)
    end if

    deallocate(rbufslp)

  end subroutine update_acc_vars_veg_es

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

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_ws
    !-----------------------------------------------------------------------
    this%h2ocan(begp:endp) = spval 
    call hist_addfld1d (fname='H2OCAN', units='mm',  &
         avgflag='A', long_name='intercepted water', &
         ptr_patch=this%h2ocan, set_lake=0._r8)

    this%q_ref2m(begp:endp) = spval
    call hist_addfld1d (fname='Q2M', units='kg/kg',  &
         avgflag='A', long_name='2m specific humidity', &
         ptr_patch=this%q_ref2m)

    this%rh_ref2m(begp:endp) = spval
    call hist_addfld1d (fname='RH2M', units='%',  &
         avgflag='A', long_name='2m relative humidity', &
         ptr_patch=this%rh_ref2m)

    this%rh_ref2m_r(begp:endp) = spval
    call hist_addfld1d (fname='RH2M_R', units='%',  &
         avgflag='A', long_name='Rural 2m specific humidity', &
         ptr_patch=this%rh_ref2m_r, set_spec=spval)

    this%rh_ref2m_u(begp:endp) = spval
    call hist_addfld1d (fname='RH2M_U', units='%',  &
         avgflag='A', long_name='Urban 2m relative humidity', &
         ptr_patch=this%rh_ref2m_u, set_nourb=spval)

    this%rh_af(begp:endp) = spval
    call hist_addfld1d (fname='RHAF', units='fraction', &
         avgflag='A', long_name='fractional humidity of canopy air', &
         ptr_patch=this%rh_af, set_spec=spval, default='inactive')

    if (use_cn) then
       this%fwet(begp:endp) = spval
       call hist_addfld1d (fname='FWET', units='proportion', &
            avgflag='A', long_name='fraction of canopy that is wet', &
            ptr_patch=this%fwet, default='inactive')
    end if

    if (use_cn) then
       this%fdry(begp:endp) = spval
       call hist_addfld1d (fname='FDRY', units='proportion', &
            avgflag='A', long_name='fraction of foliage that is green and dry', &
            ptr_patch=this%fdry, default='inactive')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_ws
    !-----------------------------------------------------------------------
    do p = begp,endp
       this%h2ocan(begp:endp) = 0._r8
       this%fwet(begp:endp)   = 0._r8
       this%fdry(begp:endp)   = 0._r8
    end do
    
  end subroutine veg_ws_init
    
  !------------------------------------------------------------------------
  subroutine veg_ws_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_water_state)    :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='H2OCAN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='canopy water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2ocan)

    call restartvar(ncid=ncid, flag=flag, varname='FWET', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='fraction of canopy that is wet (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fwet)
   
  end subroutine veg_ws_restart
  
  !------------------------------------------------------------------------
  subroutine veg_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ws_clean
  
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
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_cs
    !-----------------------------------------------------------------------

    if (use_fates) then
       ! no veg-level carbon state history fields defined by host model
    
    else if (carbon_type == 'c12') then
       this%leafc(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
             avgflag='A', long_name='leaf C', &
             ptr_patch=this%leafc)

       this%leafc_storage(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='leaf C storage', &
             ptr_patch=this%leafc_storage, default='inactive')

       this%leafc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
             avgflag='A', long_name='leaf C transfer', &
             ptr_patch=this%leafc_xfer, default='inactive')

       this%frootc(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
             avgflag='A', long_name='fine root C', &
             ptr_patch=this%frootc)

       this%frootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='fine root C storage', &
             ptr_patch=this%frootc_storage, default='inactive')

       this%frootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='fine root C transfer', &
             ptr_patch=this%frootc_xfer, default='inactive')

       this%livestemc(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
             avgflag='A', long_name='live stem C', &
             ptr_patch=this%livestemc)

       this%livestemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='live stem C storage', &
             ptr_patch=this%livestemc_storage, default='inactive')

       this%livestemc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
             avgflag='A', long_name='live stem C transfer', &
             ptr_patch=this%livestemc_xfer, default='inactive')

       this%deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
             avgflag='A', long_name='dead stem C', &
             ptr_patch=this%deadstemc)

       this%deadstemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='dead stem C storage', &
             ptr_patch=this%deadstemc_storage, default='inactive')

       this%deadstemc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
             avgflag='A', long_name='dead stem C transfer', &
             ptr_patch=this%deadstemc_xfer, default='inactive')

       this%livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C', &
             ptr_patch=this%livecrootc)

       this%livecrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C storage', &
             ptr_patch=this%livecrootc_storage, default='inactive')

       this%livecrootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C transfer', &
             ptr_patch=this%livecrootc_xfer, default='inactive')

       this%deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C', &
             ptr_patch=this%deadcrootc)

       this%deadcrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C storage', &
             ptr_patch=this%deadcrootc_storage,  default='inactive')

       this%deadcrootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C transfer', &
             ptr_patch=this%deadcrootc_xfer, default='inactive')

       this%gresp_storage(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='growth respiration storage', &
             ptr_patch=this%gresp_storage, default='inactive')

       this%gresp_xfer(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
             avgflag='A', long_name='growth respiration transfer', &
             ptr_patch=this%gresp_xfer, default='inactive')

       this%cpool(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
             avgflag='A', long_name='temporary photosynthate C pool', &
             ptr_patch=this%cpool)

       this%xsmrpool(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
             avgflag='A', long_name='temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool, default='active')

       if (crop_prog) then
          this%grainc(begp:endp) = spval
          call hist_addfld1d (fname='GRAINC', units='gC/m^2', &
                avgflag='A', long_name='grain C', &
                ptr_patch=this%grainc, default='inactive')

          this%cropseedc_deficit(begp:endp) = spval
          call hist_addfld1d (fname='CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit)
       end if
             
       this%ctrunc(begp:endp) = spval
       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
             avgflag='A', long_name='patch-level sink for C truncation', &
             ptr_patch=this%ctrunc, default='inactive')

       this%woodc(begp:endp) = spval
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
             avgflag='A', long_name='wood C', &
             ptr_patch=this%woodc)

       this%dispvegc(begp:endp) = spval
       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
             avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispvegc)

       this%storvegc(begp:endp) = spval
       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
             avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storvegc)

       this%totvegc(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
             avgflag='A', long_name='total vegetation carbon, excluding cpool', &
             ptr_patch=this%totvegc)

       this%totpftc(begp:endp) = spval
       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
             avgflag='A', long_name='total patch-level carbon, including cpool', &
             ptr_patch=this%totpftc)

       this%totvegc_abg(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC_ABG', units='gC/m^2', &
            avgflag='A', long_name='total aboveground vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_abg)
       
       ! end of c12 block

    else if ( carbon_type == 'c13' ) then
       this%leafc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C', &
             ptr_patch=this%leafc)

       this%leafc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C storage', &
             ptr_patch=this%leafc_storage, default='inactive')

       this%leafc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C transfer', &
             ptr_patch=this%leafc_xfer, default='inactive')

       this%frootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C', &
             ptr_patch=this%frootc)

       this%frootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C storage', &
             ptr_patch=this%frootc_storage, default='inactive')

       this%frootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C transfer', &
             ptr_patch=this%frootc_xfer, default='inactive')

       this%livestemc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C', &
             ptr_patch=this%livestemc)

       this%livestemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C storage', &
             ptr_patch=this%livestemc_storage, default='inactive')

       this%livestemc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C transfer', &
             ptr_patch=this%livestemc_xfer, default='inactive')

       this%deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C', &
             ptr_patch=this%deadstemc)

       this%deadstemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C storage', &
             ptr_patch=this%deadstemc_storage, default='inactive')

       this%deadstemc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C transfer', &
             ptr_patch=this%deadstemc_xfer, default='inactive')

       this%livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C', &
             ptr_patch=this%livecrootc)

       this%livecrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C storage', &
             ptr_patch=this%livecrootc_storage, default='inactive')

       this%livecrootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C transfer', &
             ptr_patch=this%livecrootc_xfer, default='inactive')

       this%deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C', &
             ptr_patch=this%deadcrootc)

       this%deadcrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C storage', &
             ptr_patch=this%deadcrootc_storage,  default='inactive')

       this%deadcrootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C transfer', &
             ptr_patch=this%deadcrootc_xfer, default='inactive')

       this%gresp_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 growth respiration storage', &
             ptr_patch=this%gresp_storage, default='inactive')

       this%gresp_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 growth respiration transfer', &
             ptr_patch=this%gresp_xfer, default='inactive')

       this%cpool(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
             avgflag='A', long_name='C13 temporary photosynthate C pool', &
             ptr_patch=this%cpool)

       this%xsmrpool(begp:endp) = spval
       call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
             avgflag='A', long_name='C13 temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool)

       this%ctrunc(begp:endp) = spval
       call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
             avgflag='A', long_name='C13 patch-level sink for C truncation', &
             ptr_patch=this%ctrunc)

       this%dispvegc(begp:endp) = spval
       call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispvegc)

       this%storvegc(begp:endp) = spval
       call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storvegc)

       this%totvegc(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
             ptr_patch=this%totvegc)

       this%totpftc(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total patch-level carbon, including cpool', &
             ptr_patch=this%totpftc)

       if (use_crop) then
          this%grainc(begp:endp) = spval
          call hist_addfld1d (fname='C13_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C13 grain C (does not equal yield)', &
               ptr_patch=this%grainc)

          this%cropseedc_deficit(begp:endp) = spval
          call hist_addfld1d (fname='C13_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C13 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit)
       end if
       ! end of c13 block

    else if ( carbon_type == 'c14' ) then

       this%leafc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C', &
             ptr_patch=this%leafc)

       this%leafc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C storage', &
             ptr_patch=this%leafc_storage, default='inactive')

       this%leafc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C transfer', &
             ptr_patch=this%leafc_xfer, default='inactive')

       this%frootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C', &
             ptr_patch=this%frootc)

       this%frootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C storage', &
             ptr_patch=this%frootc_storage, default='inactive')

       this%frootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C transfer', &
             ptr_patch=this%frootc_xfer, default='inactive')

       this%livestemc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C', &
             ptr_patch=this%livestemc)

       this%livestemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C storage', &
             ptr_patch=this%livestemc_storage, default='inactive')

       this%livestemc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C transfer', &
             ptr_patch=this%livestemc_xfer, default='inactive')

       this%deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C', &
             ptr_patch=this%deadstemc)

       this%deadstemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C storage', &
             ptr_patch=this%deadstemc_storage, default='inactive')

       this%deadstemc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C transfer', &
             ptr_patch=this%deadstemc_xfer, default='inactive')

       this%livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C', &
             ptr_patch=this%livecrootc)

       this%livecrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C storage', &
             ptr_patch=this%livecrootc_storage, default='inactive')

       this%livecrootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C transfer', &
             ptr_patch=this%livecrootc_xfer, default='inactive')

       this%deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C', &
             ptr_patch=this%deadcrootc)

       this%deadcrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C storage', &
             ptr_patch=this%deadcrootc_storage,  default='inactive')

       this%deadcrootc_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C transfer', &
             ptr_patch=this%deadcrootc_xfer, default='inactive')

       this%gresp_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 growth respiration storage', &
             ptr_patch=this%gresp_storage, default='inactive')

       this%gresp_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 growth respiration transfer', &
             ptr_patch=this%gresp_xfer, default='inactive')

       this%cpool(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL', units='gC14/m^2', &
             avgflag='A', long_name='C14 temporary photosynthate C pool', &
             ptr_patch=this%cpool)

       this%xsmrpool(begp:endp) = spval
       call hist_addfld1d (fname='C14_XSMRPOOL', units='gC14/m^2', &
             avgflag='A', long_name='C14 temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool)

       this%ctrunc(begp:endp) = spval
       call hist_addfld1d (fname='C14_PFT_CTRUNC', units='gC14/m^2', &
             avgflag='A', long_name='C14 patch-level sink for C truncation', &
             ptr_patch=this%ctrunc)

       this%dispvegc(begp:endp) = spval
       call hist_addfld1d (fname='C14_DISPVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispvegc)

       this%storvegc(begp:endp) = spval
       call hist_addfld1d (fname='C14_STORVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storvegc)

       this%totvegc(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total vegetation carbon, excluding cpool', &
             ptr_patch=this%totvegc)

       this%totpftc(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTPFTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total patch-level carbon, including cpool', &
             ptr_patch=this%totpftc)
    
       if (use_crop) then
          this%grainc(begp:endp) = spval
          call hist_addfld1d (fname='C14_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C14 grain C (does not equal yield)', &
               ptr_patch=this%grainc)
          
          this%cropseedc_deficit(begp:endp) = spval
          call hist_addfld1d (fname='C14_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C14 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit)
       end if
       ! end of c14 block
    
    endif  ! fates or c12 or c13 or c14
    
    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_cs
    !-----------------------------------------------------------------------

    this%species = species_from_string(carbon_type)

    if ( .not. use_fates ) then
       do p = begp,endp

          this%leafcmax(p) = 0._r8

          l = veg_pp%landunit(p)
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

             if (veg_pp%itype(p) == noveg) then
                this%leafc(p)         = 0._r8
                this%leafc_storage(p) = 0._r8
             else
                if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                   this%leafc(p)         = 1._r8 * ratio
                   this%leafc_storage(p) = 0._r8
                else if (veg_pp%itype(p) >= npcropmin) then ! prognostic crop types
                   this%leafc(p) = 0._r8
                   this%leafc_storage(p) = 0._r8
                else
                   this%leafc(p) = 0._r8
                   this%leafc_storage(p) = 1._r8 * ratio
                end if
             end if
             this%leafc_xfer(p) = 0._r8

             this%frootc(p)            = 0._r8 
             this%frootc_storage(p)    = 0._r8 
             this%frootc_xfer(p)       = 0._r8 

             this%livestemc(p)         = 0._r8 
             this%livestemc_storage(p) = 0._r8 
             this%livestemc_xfer(p)    = 0._r8 

             if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
                this%deadstemc(p) = 0.1_r8 * ratio
             else
                this%deadstemc(p) = 0._r8 
             end if
             this%deadstemc_storage(p)  = 0._r8 
             this%deadstemc_xfer(p)     = 0._r8

             if (nu_com .ne. 'RD') then
                ! ECA competition calculate root NP uptake as a function of fine root biomass
                ! better to initialize root CNP pools with a non-zero value
                if (veg_pp%itype(p) .ne. noveg) then
                   if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                      this%leafc(p) = 20._r8 * ratio
                      this%leafc_storage(p) = 0._r8
                      this%frootc(p) = 20._r8 * ratio
                      this%frootc_storage(p) = 0._r8
                   else
                      this%leafc(p) = 0._r8 
                      this%leafc_storage(p) = 20._r8 * ratio
                      this%frootc(p) = 0._r8
                      this%frootc_storage(p) = 20._r8 * ratio
                   end if
                end if
             end if

             this%livecrootc(p)         = 0._r8 
             this%livecrootc_storage(p) = 0._r8 
             this%livecrootc_xfer(p)    = 0._r8 

             this%deadcrootc(p)         = 0._r8 
             this%deadcrootc_storage(p) = 0._r8 
             this%deadcrootc_xfer(p)    = 0._r8 

             this%gresp_storage(p)      = 0._r8 
             this%gresp_xfer(p)         = 0._r8 

             this%cpool(p)              = 0._r8 
             this%xsmrpool(p)           = 0._r8 
             this%ctrunc(p)             = 0._r8 
             this%dispvegc(p)           = 0._r8 
             this%storvegc(p)           = 0._r8 
             this%totpftc(p)            = 0._r8 
             this%woodc(p)              = 0._r8

             if ( crop_prog )then
                this%grainc(p)            = 0._r8 
                this%grainc_storage(p)    = 0._r8 
                this%grainc_xfer(p)       = 0._r8 
             end if
             this%cropseedc_deficit(p) = 0._r8
             ! calculate totvegc explicitly so that it is available for the isotope 
             ! code on the first time step.

             this%totvegc(p) = &
                  this%leafc(p)              + &
                  this%leafc_storage(p)      + &
                  this%leafc_xfer(p)         + &
                  this%frootc(p)             + &
                  this%frootc_storage(p)     + &
                  this%frootc_xfer(p)        + &
                  this%livestemc(p)          + &
                  this%livestemc_storage(p)  + &
                  this%livestemc_xfer(p)     + &
                  this%deadstemc(p)          + &
                  this%deadstemc_storage(p)  + &
                  this%deadstemc_xfer(p)     + &
                  this%livecrootc(p)         + &
                  this%livecrootc_storage(p) + &
                  this%livecrootc_xfer(p)    + &
                  this%deadcrootc(p)         + &
                  this%deadcrootc_storage(p) + &
                  this%deadcrootc_xfer(p)    + &
                  this%gresp_storage(p)      + &
                  this%gresp_xfer(p)         + &
                  this%cpool(p)

             if ( crop_prog )then
                this%totvegc(p) =  this%totvegc(p) + &
                     this%grainc(p)                            + &
                     this%grainc_storage(p)                    + &
                     this%grainc_xfer(p)
             end if
          endif ! is soil or crop
       end do ! begp:endp

       ! Set temporary veg filter for special landunits
       num_special_veg = 0
       do p = begp,endp
          l = veg_pp%landunit(p)
          if (lun_pp%ifspecial(l)) then
             num_special_veg = num_special_veg + 1
             special_veg(num_special_veg) = p
          end if
       end do
       ! reset state variables for special landunits
       value_veg = 0._r8
       do fp = 1,num_special_veg
          p = special_veg(fp)

          this%leafc(p)                = value_veg
          this%leafc_storage(p)        = value_veg
          this%leafc_xfer(p)           = value_veg
          this%frootc(p)               = value_veg
          this%frootc_storage(p)       = value_veg
          this%frootc_xfer(p)          = value_veg
          this%livestemc(p)            = value_veg
          this%livestemc_storage(p)    = value_veg
          this%livestemc_xfer(p)       = value_veg
          this%deadstemc(p)            = value_veg
          this%deadstemc_storage(p)    = value_veg
          this%deadstemc_xfer(p)       = value_veg
          this%livecrootc(p)           = value_veg
          this%livecrootc_storage(p)   = value_veg
          this%livecrootc_xfer(p)      = value_veg
          this%deadcrootc(p)           = value_veg
          this%deadcrootc_storage(p)   = value_veg
          this%deadcrootc_xfer(p)      = value_veg
          this%gresp_storage(p)        = value_veg
          this%gresp_xfer(p)           = value_veg
          this%cpool(p)                = value_veg
          this%xsmrpool(p)             = value_veg
          this%ctrunc(p)               = value_veg
          this%dispvegc(p)             = value_veg
          this%storvegc(p)             = value_veg
          this%totvegc(p)              = value_veg
          this%totpftc(p)              = value_veg
          this%woodc(p)                = value_veg
          this%totvegc_abg(p)          = value_veg
       
          if ( crop_prog ) then
             this%grainc(p)            = value_veg
             this%grainc_storage(p)    = value_veg
             this%grainc_xfer(p)       = value_veg
             this%cropseedc_deficit(p) = value_veg
          end if
       end do
    endif ! .not. use_fates

    end subroutine veg_cs_init
    
  !------------------------------------------------------------------------
  subroutine veg_cs_restart ( this,  bounds, ncid, flag, carbon_type, c12_veg_cs, cnstate_vars)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation carbon state information to/from restart file.
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state)   :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    character(len=3) , intent(in)    :: carbon_type ! 'c12' or 'c13' or 'c14'
    type (vegetation_carbon_state) , intent(in), optional :: c12_veg_cs 
    type (cnstate_type)            , intent(in)           :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar    ! determine if variable is on initial file
    character(len=128) :: varname    ! temporary
    integer            :: i,l,c,p      ! indices
    real(r8)           :: c3_del13c  ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8)           :: c4_del13c  ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8)           :: c3_r1      ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8)           :: c4_r1      ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8)           :: c3_r2      ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8)           :: c4_r2      ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8)           :: m_veg      ! multiplier for the exit_spinup code
    integer            :: idata
    logical            :: exit_spinup  = .false.
    logical            :: enter_spinup = .false.
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !-----------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_veg_cs)) then
          call endrun(msg=' ERROR: for C14 must pass in c12_veg_cs as argument' //&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    if ( .not. use_fates ) then
       !--------------------------------
       ! C12 vegetation carbon state variables
       !--------------------------------
       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='leafc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_storage) 
          
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_storage) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer) 

          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage) 

          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer) 

          call restartvar(ncid=ncid, flag=flag, varname='cpool', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%cpool) 

          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool) 

          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%ctrunc) 

          call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc) 

          if (crop_prog) then
             call restartvar(ncid=ncid, flag=flag,  varname='grainc', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grainc)

             call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C storage', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grainc_storage)

             call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C transfer', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer)

             call restartvar(ncid=ncid, flag=flag, varname='cropseedc_deficit', xtype=ncd_double,  &
                  dim1name='pft', long_name='pool for seeding new crop growth', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit)
          end if ! crop_prog

       end if  ! c12

       !--------------------------------
       ! C13 vegetation carbon state variables 
       !--------------------------------
       if ( carbon_type == 'c13')  then
          if ( .not. is_restart() .and. get_nstep() == 1 ) then
             c3_del13c = -28._r8
             c4_del13c = -13._r8
             c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
             c3_r2 = c3_r1/(1._r8 + c3_r1)
             c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
             c4_r2 = c4_r1/(1._r8 + c4_r1)

             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%grainc(i)            = c12_veg_cs%grainc(i)         * c3_r2
                   this%grainc_storage(i)    = c12_veg_cs%grainc_storage(i) * c3_r2
                   this%grainc_xfer(i)       = c12_veg_cs%grainc_xfer(i)    * c3_r2
                   this%dispvegc(i)          = c12_veg_cs%dispvegc(i)       * c3_r2
                   this%storvegc(i)          = c12_veg_cs%storvegc(i)       * c3_r2
                   this%totvegc(i)           = c12_veg_cs%totvegc(i)        * c3_r2
                   this%totpftc(i)           = c12_veg_cs%totpftc(i)        * c3_r2
                   this%woodc(i)             = c12_veg_cs%woodc(i)          * c3_r2
                else
                   this%grainc(i)            = c12_veg_cs%grainc(i)         * c4_r2
                   this%grainc_storage(i)    = c12_veg_cs%grainc_storage(i) * c4_r2
                   this%grainc_xfer(i)       = c12_veg_cs%grainc_xfer(i)    * c4_r2
                   this%dispvegc(i)          = c12_veg_cs%dispvegc(i)       * c4_r2
                   this%storvegc(i)          = c12_veg_cs%storvegc(i)       * c4_r2
                   this%totvegc(i)           = c12_veg_cs%totvegc(i)        * c4_r2
                   this%totpftc(i)           = c12_veg_cs%totpftc(i)        * c4_r2
                   this%woodc(i)             = c12_veg_cs%woodc(i)          * c4_r2
                end if
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc)
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leafc(i) = c12_veg_cs%leafc(i) * c3_r2
                else
                   this%leafc(i) = c12_veg_cs%leafc(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leafc_storage(i) = c12_veg_cs%leafc_storage(i) * c3_r2
                else
                   this%leafc_storage(i) = c12_veg_cs%leafc_storage(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leafc_xfer(i) = c12_veg_cs%leafc_xfer(i) * c3_r2
                else
                   this%leafc_xfer(i) = c12_veg_cs%leafc_xfer(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%frootc(i) = c12_veg_cs%frootc(i) * c3_r2
                else
                   this%frootc(i) = c12_veg_cs%frootc(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%frootc_storage(i) = c12_veg_cs%frootc_storage(i) * c3_r2
                else
                   this%frootc_storage(i) = c12_veg_cs%frootc_storage(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%frootc_xfer(i) = c12_veg_cs%frootc_xfer(i) * c3_r2
                else
                   this%frootc_xfer(i) = c12_veg_cs%frootc_xfer(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestemc(i) = c12_veg_cs%livestemc(i) * c3_r2
                else
                   this%livestemc(i) = c12_veg_cs%livestemc(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestemc_storage(i) = c12_veg_cs%livestemc_storage(i) * c3_r2
                else
                   this%livestemc_storage(i) = c12_veg_cs%livestemc_storage(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestemc_xfer(i) = c12_veg_cs%livestemc_xfer(i) * c3_r2
                else
                   this%livestemc_xfer(i) = c12_veg_cs%livestemc_xfer(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstemc(i) = c12_veg_cs%deadstemc(i) * c3_r2
                else
                   this%deadstemc(i) = c12_veg_cs%deadstemc(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstemc_storage(i) = c12_veg_cs%deadstemc_storage(i) * c3_r2
                else
                   this%deadstemc_storage(i) = c12_veg_cs%deadstemc_storage(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstemc_xfer(i) = c12_veg_cs%deadstemc_xfer(i) * c3_r2
                else
                   this%deadstemc_xfer(i) = c12_veg_cs%deadstemc_xfer(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecrootc(i) = c12_veg_cs%livecrootc(i) * c3_r2
                else
                   this%livecrootc(i) = c12_veg_cs%livecrootc(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecrootc_storage(i) = c12_veg_cs%livecrootc_storage(i) * c3_r2
                else
                   this%livecrootc_storage(i) = c12_veg_cs%livecrootc_storage(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecrootc_xfer(i) = c12_veg_cs%livecrootc_xfer(i) * c3_r2
                else
                   this%livecrootc_xfer(i) = c12_veg_cs%livecrootc_xfer(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcrootc(i) = c12_veg_cs%deadcrootc(i) * c3_r2
                else
                   this%deadcrootc(i) = c12_veg_cs%deadcrootc(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcrootc_storage(i) = c12_veg_cs%deadcrootc_storage(i) * c3_r2
                else
                   this%deadcrootc_storage(i) = c12_veg_cs%deadcrootc_storage(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcrootc_xfer(i) = c12_veg_cs%deadcrootc_xfer(i) * c3_r2
                else
                   this%deadcrootc_xfer(i) = c12_veg_cs%deadcrootc_xfer(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%gresp_storage(i) = c12_veg_cs%gresp_storage(i) * c3_r2
                else
                   this%gresp_storage(i) = c12_veg_cs%gresp_storage(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%gresp_xfer(i) = c12_veg_cs%gresp_xfer(i) * c3_r2
                else
                   this%gresp_xfer(i) = c12_veg_cs%gresp_xfer(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='cpool_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%cpool) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%cpool with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%cpool(i) = c12_veg_cs%cpool(i) * c3_r2
                else
                   this%cpool(i) = c12_veg_cs%cpool(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%xsmrpool with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%xsmrpool(i) = c12_veg_cs%xsmrpool(i) * c3_r2
                else
                   this%xsmrpool(i) = c12_veg_cs%xsmrpool(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%ctrunc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%ctrunc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%ctrunc(i) = c12_veg_cs%ctrunc(i) * c3_r2
                else
                   this%ctrunc(i) = c12_veg_cs%ctrunc(i) * c4_r2
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing carbonstate_vars %totvegc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%totvegc(i) = c12_veg_cs%totvegc(i) * c3_r2
                else
                   this%totvegc(i) = c12_veg_cs%totvegc(i) * c4_r2
                endif
             end do
          end if
          
       endif ! C13 block

       !--------------------------------
       ! C14 pft carbon state variables 
       !--------------------------------
       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leafc(i) /= spval .and. &
                     .not. isnan(this%leafc(i)) ) then
                   this%leafc(i) = c12_veg_cs%leafc(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_storage with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leafc_storage(i) /= spval .and. &
                     .not. isnan(this%leafc_storage(i)) ) then
                   this%leafc_storage(i) = c12_veg_cs%leafc_storage(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_xfer with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leafc_xfer(i) /= spval .and. .not. isnan(this%leafc_xfer(i)) ) then
                   this%leafc_xfer(i) = c12_veg_cs%leafc_xfer(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%frootc(i) /= spval .and. &
                     .not. isnan(this%frootc(i)) ) then
                   this%frootc(i) = c12_veg_cs%frootc(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_storage with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%frootc_storage(i) /= spval .and. &
                     .not. isnan(this%frootc_storage(i)) ) then
                   this%frootc_storage(i) = c12_veg_cs%frootc_storage(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_xfer with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%frootc_xfer(i) /= spval .and. &
                     .not. isnan(this%frootc_xfer(i)) ) then
                   this%frootc_xfer(i) = c12_veg_cs%frootc_xfer(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestemc(i) /= spval .and. .not. isnan(this%livestemc(i)) ) then
                   this%livestemc(i) = c12_veg_cs%livestemc(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_storage with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestemc_storage(i) /= spval .and. .not. isnan(this%livestemc_storage(i)) ) then
                   this%livestemc_storage(i) = c12_veg_cs%livestemc_storage(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestemc_xfer(i) /= spval .and. .not. isnan(this%livestemc_xfer(i)) ) then
                   this%livestemc_xfer(i) = c12_veg_cs%livestemc_xfer(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstemc(i) /= spval .and. .not. isnan(this%deadstemc(i)) ) then
                   this%deadstemc(i) = c12_veg_cs%deadstemc(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstemc_storage(i) /= spval .and. .not. isnan(this%deadstemc_storage(i)) ) then
                   this%deadstemc_storage(i) = c12_veg_cs%deadstemc_storage(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstemc_xfer(i) /= spval .and. .not. isnan(this%deadstemc_xfer(i)) ) then
                   this%deadstemc_xfer(i) = c12_veg_cs%deadstemc_xfer(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecrootc(i) /= spval .and. .not. isnan(this%livecrootc(i)) ) then
                   this%livecrootc(i) = c12_veg_cs%livecrootc(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecrootc_storage(i) /= spval .and. .not. isnan(this%livecrootc_storage(i)) ) then
                   this%livecrootc_storage(i) = c12_veg_cs%livecrootc_storage(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecrootc_xfer(i) /= spval .and. .not. isnan(this%livecrootc_xfer(i)) ) then
                   this%livecrootc_xfer(i) = c12_veg_cs%livecrootc_xfer(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcrootc(i) /= spval .and. .not. isnan(this%deadcrootc(i)) ) then
                   this%deadcrootc(i) = c12_veg_cs%deadcrootc(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcrootc_storage(i) /= spval .and. .not. isnan(this%deadcrootc_storage(i)) ) then
                   this%deadcrootc_storage(i) = c12_veg_cs%deadcrootc_storage(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog) 'initializing this%deadcrootc_xfer with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcrootc_xfer(i) /= spval .and. .not. isnan(this%deadcrootc_xfer(i)) ) then
                   this%deadcrootc_xfer(i) = c12_veg_cs%deadcrootc_xfer(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_storage with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%gresp_storage(i) /= spval .and. .not. isnan(this%gresp_storage(i)) ) then
                   this%gresp_storage(i) = c12_veg_cs%gresp_storage(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_xfer with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%gresp_xfer(i) /= spval .and. .not. isnan(this%gresp_xfer(i)) ) then
                   this%gresp_xfer(i) = c12_veg_cs%gresp_xfer(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='cpool_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%cpool) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%cpool with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%cpool(i) /= spval .and. .not. isnan(this%cpool(i)) ) then
                   this%cpool(i) = c12_veg_cs%cpool(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%xsmrpool with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%xsmrpool(i) /= spval .and. .not. isnan(this%xsmrpool(i)) ) then
                   this%xsmrpool(i) = c12_veg_cs%xsmrpool(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%ctrunc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%ctrunc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%ctrunc(i) /= spval .and. .not. isnan(this%ctrunc(i)) ) then
                   this%ctrunc(i) = c12_veg_cs%ctrunc(i) * c14ratio
                endif
             end do
          end if

          call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%totvegc with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%totvegc(i) /= spval .and. .not. isnan(this%totvegc(i)) ) then
                   this%totvegc(i) = c12_veg_cs%totvegc(i) * c14ratio
                endif
             end do
          end if
          
       endif  ! C14 block

    endif  ! .not. use_fates

    !--------------------------------
    ! Spinup state operations happen on restart write/read
    !--------------------------------
    ! the spinup_state variable is being written by the column-level carbon state restart
    ! routine, so only need to handle the reading part here    
    if (carbon_type == 'c12'  .or. carbon_type == 'c14') then
        if (flag == 'read') then
           call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
                long_name='Spinup state of the model that wrote this restart file: ' &
                // ' 0 = normal model mode, 1 = AD spinup', units='', &
                interpinic_flag='copy', readvar=readvar,  data=idata)
           if (readvar) then
              restart_file_spinup_state = idata
           else
              ! assume, for sake of backwards compatibility, that if spinup_state is not in 
              ! the restart file then current model state is the same as prior model state
              restart_file_spinup_state = spinup_state
              if ( masterproc ) then
                 write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                      // ' on spinup state used to generate the restart file. '
                 write(iulog,*) '   Assuming the same as current setting: ', spinup_state
              end if
           end if
        end if

        ! now compare the model and restart file spinup states, and either take the 
        ! model into spinup mode or out of it if they are not identical
        ! taking model out of spinup mode requires multiplying each pool 
        ! by the associated AD factor.
        ! putting model into spinup mode requires dividing each pool 
        ! by the associated AD factor.
        ! only allow this to occur on first timestep of model run.
        
        if (flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
           if (spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
              if ( masterproc ) write(iulog,*) ' veg_cs_restart: taking deadwood pools out of AD spinup mode'
              exit_spinup = .true.
           else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
              if ( masterproc ) write(iulog,*) ' veg_cs_restart: taking deadwood pools into AD spinup mode'
              enter_spinup = .true.
           else
              call endrun(msg=' veg_cs_restart: error in entering/exiting spinup.  spinup_state ' &
                   // ' != restart_file_spinup_state, but do not know what to do'//&
                   errMsg(__FILE__, __LINE__))
           end if
           if (get_nstep() >= 2) then
              call endrun(msg=' veg_cs_restart: error in entering/exiting spinup - should occur only when nstep = 1'//&
                   errMsg(__FILE__, __LINE__))
           endif
           do i = bounds%begp, bounds%endp
              if (exit_spinup) then 
                 m_veg = spinup_mortality_factor
              else if (enter_spinup) then 
                 m_veg = 1._r8 / spinup_mortality_factor
              end if
              this%deadstemc(i)  = this%deadstemc(i) * m_veg
              this%deadcrootc(i) = this%deadcrootc(i) * m_veg
           end do
        end if ! read
     end if ! c12 or c14 (PET: why not c13?) 

  end subroutine veg_cs_restart

  !-----------------------------------------------------------------------
  subroutine veg_cs_summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_cs)
    !
    ! !DESCRIPTION:
    ! Vegetation-level carbon state summary calculations
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state)                :: this
    type(bounds_type)         , intent(in)    :: bounds          
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                   , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                   , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type (column_carbon_state), intent(inout) :: col_cs          ! column-level state for p2c
    
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p             ! indices
    integer  :: fp              ! filter indices
    real(r8) :: maxdepth        ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    if (use_fates) return

    do fp = 1,num_soilp
       p = filter_soilp(fp)
       ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
       this%dispvegc(p) =        &
            this%leafc(p)      + &
            this%frootc(p)     + &
            this%livestemc(p)  + &
            this%deadstemc(p)  + &
            this%livecrootc(p) + &
            this%deadcrootc(p)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
       this%storvegc(p) =                &
            this%cpool(p)              + &
            this%leafc_storage(p)      + &
            this%frootc_storage(p)     + &
            this%livestemc_storage(p)  + &
            this%deadstemc_storage(p)  + &
            this%livecrootc_storage(p) + &
            this%deadcrootc_storage(p) + &
            this%leafc_xfer(p)         + &
            this%frootc_xfer(p)        + &
            this%livestemc_xfer(p)     + &
            this%deadstemc_xfer(p)     + &
            this%livecrootc_xfer(p)    + &
            this%deadcrootc_xfer(p)    + &
            this%gresp_storage(p)      + &
            this%gresp_xfer(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%storvegc(p) =            &
               this%storvegc(p)       + &
               this%grainc_storage(p) + &
               this%grainc_xfer(p)

          this%dispvegc(p) =            &
               this%dispvegc(p)       + &
               this%grainc(p)
       end if

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totvegc(p) = &
            this%dispvegc(p) + &
            this%storvegc(p)

       ! total pft-level carbon, including xsmrpool, ctrunc
       this%totpftc(p) = &
            this%totvegc(p) + &
            this%xsmrpool(p) + &
            this%ctrunc(p)

       ! (WOODC) - wood C
       this%woodc(p) = &
            this%deadstemc(p)    + &
            this%livestemc(p)    + &
            this%deadcrootc(p)   + &
            this%livecrootc(p)

       this%totvegc_abg(p) = &
               this%leafc(p)              + &
               this%leafc_storage(p)      + &
               this%leafc_xfer(p)         + &
               this%livestemc(p)          + &
               this%livestemc_storage(p)  + &
               this%livestemc_xfer(p)     + &
               this%deadstemc(p)          + &
               this%deadstemc_storage(p)  + &
               this%deadstemc_xfer(p)
    end do ! filtered veg list
    
    ! a few vegetation-to-column summaries
    call p2c(bounds, num_soilc, filter_soilc, &
         this%totpftc(bounds%begp:bounds%endp), &
         col_cs%totpftc(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegc(bounds%begp:bounds%endp), &
         col_cs%totvegc(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegc_abg(bounds%begp:bounds%endp), &
         col_cs%totvegc_abg(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%cropseedc_deficit(bounds%begp:bounds%endp), &
         col_cs%cropseedc_deficit(bounds%begc:bounds%endc))

  end subroutine veg_cs_summary
    
  !-----------------------------------------------------------------------
  subroutine veg_cs_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegc(p)   = 0._r8
       this%storvegc(p)   = 0._r8
       this%totpftc(p)    = 0._r8
    end do

  end subroutine veg_cs_zerodwt
  !------------------------------------------------------------------------
  subroutine veg_cs_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_cs_clean
  
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

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_ns
    !-----------------------------------------------------------------------
    if (crop_prog) then
       this%grainn(begp:endp) = spval
       call hist_addfld1d (fname='GRAINN', units='gN/m^2', &
            avgflag='A', long_name='grain N', &
            ptr_patch=this%grainn, default='inactive')

       this%cropseedn_deficit(begp:endp) = spval
       call hist_addfld1d (fname='CROPSEEDN_DEFICIT', units='gN/m^2', &
            avgflag='A', long_name='N used for crop seed that needs to be repaid', &
            ptr_patch=this%cropseedn_deficit)
    end if

    this%leafn(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN', units='gN/m^2', &
         avgflag='A', long_name='leaf N', &
         ptr_patch=this%leafn)

    this%leafn_storage(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='leaf N storage', &
         ptr_patch=this%leafn_storage, default='inactive')

    this%leafn_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_XFER', units='gN/m^2', &
         avgflag='A', long_name='leaf N transfer', &
         ptr_patch=this%leafn_xfer, default='inactive')

    this%frootn(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN', units='gN/m^2', &
         avgflag='A', long_name='fine root N', &
         ptr_patch=this%frootn)

    this%frootn_storage(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='fine root N storage', &
         ptr_patch=this%frootn_storage, default='inactive')

    this%frootn_xfer(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='fine root N transfer', &
         ptr_patch=this%frootn_xfer, default='inactive')

    this%livestemn(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN', units='gN/m^2', &
         avgflag='A', long_name='live stem N', &
         ptr_patch=this%livestemn)

    this%livestemn_storage(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live stem N storage', &
         ptr_patch=this%livestemn_storage, default='inactive')

    this%livestemn_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live stem N transfer', &
         ptr_patch=this%livestemn_xfer, default='inactive')

    this%deadstemn(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN', units='gN/m^2', &
         avgflag='A', long_name='dead stem N', &
         ptr_patch=this%deadstemn)

    this%deadstemn_storage(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead stem N storage', &
         ptr_patch=this%deadstemn_storage, default='inactive')

    this%deadstemn_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead stem N transfer', &
         ptr_patch=this%deadstemn_xfer, default='inactive')

    this%livecrootn(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N', &
         ptr_patch=this%livecrootn)

    this%livecrootn_storage(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N storage', &
         ptr_patch=this%livecrootn_storage, default='inactive')

    this%livecrootn_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N transfer', &
         ptr_patch=this%livecrootn_xfer, default='inactive')

    this%deadcrootn(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N', &
         ptr_patch=this%deadcrootn)

    this%deadcrootn_storage(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N storage', &
         ptr_patch=this%deadcrootn_storage, default='inactive')

    this%deadcrootn_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N transfer', &
         ptr_patch=this%deadcrootn_xfer, default='inactive')

    this%retransn(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='plant pool of retranslocated N', &
         ptr_patch=this%retransn)

    this%npool(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL', units='gN/m^2', &
         avgflag='A', long_name='temporary plant N pool', &
         ptr_patch=this%npool, default='inactive')

    this%ntrunc(begp:endp) = spval
    call hist_addfld1d (fname='PFT_NTRUNC', units='gN/m^2', &
         avgflag='A', long_name='pft-level sink for N truncation', &
         ptr_patch=this%ntrunc, default='inactive')

    this%dispvegn(begp:endp) = spval
    call hist_addfld1d (fname='DISPVEGN', units='gN/m^2', &
         avgflag='A', long_name='displayed vegetation nitrogen', &
         ptr_patch=this%dispvegn)

    this%storvegn(begp:endp) = spval
    call hist_addfld1d (fname='STORVEGN', units='gN/m^2', &
         avgflag='A', long_name='stored vegetation nitrogen', &
         ptr_patch=this%storvegn)

    this%totvegn(begp:endp) = spval
    call hist_addfld1d (fname='TOTVEGN', units='gN/m^2', &
         avgflag='A', long_name='total vegetation nitrogen', &
         ptr_patch=this%totvegn)

    this%totpftn(begp:endp) = spval
    call hist_addfld1d (fname='TOTPFTN', units='gN/m^2', &
         avgflag='A', long_name='total PFT-level nitrogen', &
         ptr_patch=this%totpftn)

    this%npimbalance(begp:endp) = spval
    call hist_addfld1d (fname='leaf_npimbalance', units='gN/gP', &
         avgflag='A', long_name='leaf np imbalance partial C partial P/partial C partial N', &
         ptr_patch=this%npimbalance)

    ! Note: this is a patch-level variable, reported to the column-level
    this%plant_n_buffer(begp:endp) = spval
    call hist_addfld1d (fname='PLANTN_BUFFER', units='gN/m^2', &
         avgflag='A', long_name='plant nitrogen stored as buffer', &
         ptr_col=this%plant_n_buffer,default='inactive')
    
         
    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_ns
    !-----------------------------------------------------------------------

    do p = begp,endp

       l = veg_pp%landunit(p)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then       
          if (veg_pp%itype(p) == noveg) then
             this%leafn(p) = 0._r8
             this%leafn_storage(p) = 0._r8
          else
             this%leafn(p)         = veg_cs%leafc(p)         / veg_vp%leafcn(veg_pp%itype(p))
             this%leafn_storage(p) = veg_cs%leafc_storage(p) / veg_vp%leafcn(veg_pp%itype(p))
          end if

          this%leafn_xfer(p)        = 0._r8
          if ( crop_prog )then
             this%grainn(p)            = 0._r8
             this%grainn_storage(p)    = 0._r8
             this%grainn_xfer(p)       = 0._r8
          end if
          this%cropseedn_deficit(p) = 0._r8
          this%frootn(p)            = 0._r8
          this%frootn_storage(p)    = 0._r8
          this%frootn_xfer(p)       = 0._r8
          this%livestemn(p)         = 0._r8
          this%livestemn_storage(p) = 0._r8
          this%livestemn_xfer(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
             this%deadstemn(p) = veg_cs%deadstemc(p) / veg_vp%deadwdcn(veg_pp%itype(p))
          else
             this%deadstemn(p) = 0._r8
          end if
          
          if (nu_com .ne. 'RD') then
              ! ECA competition calculate root NP uptake as a function of fine root biomass
              ! better to initialize root CNP pools with a non-zero value
              if (veg_pp%itype(p) .ne. noveg) then
                 this%frootn(p) = veg_cs%frootc(p) / veg_vp%frootcn(veg_pp%itype(p))
                 this%frootn_storage(p) = veg_cs%frootc_storage(p) / veg_vp%frootcn(veg_pp%itype(p))
              end if
          end if

          this%deadstemn_storage(p)  = 0._r8
          this%deadstemn_xfer(p)     = 0._r8
          this%livecrootn(p)         = 0._r8
          this%livecrootn_storage(p) = 0._r8
          this%livecrootn_xfer(p)    = 0._r8
          this%deadcrootn(p)         = 0._r8
          this%deadcrootn_storage(p) = 0._r8
          this%deadcrootn_xfer(p)    = 0._r8
          this%retransn(p)           = 0._r8
          this%npool(p)              = 0._r8
          if (nstor(veg_pp%itype(p)) .gt. 1e-6_r8) then 
              this%npool(p)          = 10.0_r8
          end if
          this%ntrunc(p)             = 0._r8
          this%dispvegn(p)           = 0._r8
          this%storvegn(p)           = 0._r8
          this%totvegn(p)            = 0._r8
          this%totpftn(p)            = 0._r8          
          this%plant_n_buffer(p)     = 1._r8
       end if

       this%npimbalance(p) = 0.0_r8
       this%pnup_pfrootc(p) = 0.0_r8 
       this%benefit_pgpp_pleafc(p) = 0.0_r8   
    end do

    ! Set filter for special landunits, set values
    num_special_veg = 0
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_veg = num_special_veg + 1
          special_veg(num_special_veg) = p
       end if
    end do
    call this%SetValues (num_veg=num_special_veg, filter_veg=special_veg, value_veg=0._r8)
    
  end subroutine veg_ns_init
    
  !-----------------------------------------------------------------------
  subroutine veg_ns_restart ( this,  bounds, ncid, flag)
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for vegetation nitrogen state variables
    !
    ! !ARGUMENTS:
    class (vegetation_nitrogen_state)          :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: m_veg          ! multiplier for the exit_spinup code
    character(len=128) :: varname    ! temporary
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='leafn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='retransn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransn) 

    call restartvar(ncid=ncid, flag=flag, varname='npool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%npool) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ntrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ntrunc) 

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainn', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N storage', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_storage)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N transfer', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_xfer)
    end if
    
    call restartvar(ncid=ncid, flag=flag,  varname='npimbalance_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='npimbalance_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%npimbalance)
     call restartvar(ncid=ncid, flag=flag,  varname='pnup_pfrootc_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='pnup_pfrootc_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%pnup_pfrootc)
    call restartvar(ncid=ncid, flag=flag,  varname='benefit_pgpp_pleafc_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='benefit_pgpp_pleafc_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%benefit_pgpp_pleafc)

    !--------------------------------
    ! Spinup state operations happen on restart write/read
    !--------------------------------
    ! the spinup_state variable is being written by the column-level carbon state restart
    ! routine, so only need to handle the reading part here    
     if (flag == 'read') then
        call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
             long_name='Spinup state of the model that wrote this restart file: ' &
             // ' 0 = normal model mode, 1 = AD spinup', units='', &
             interpinic_flag='copy', readvar=readvar,  data=idata)
        if (readvar) then
           restart_file_spinup_state = idata
        else
           ! assume, for sake of backwards compatibility, that if spinup_state is not in 
           ! the restart file then current model state is the same as prior model state
           restart_file_spinup_state = spinup_state
           if ( masterproc ) then
              write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                   // ' on spinup state used to generate the restart file. '
              write(iulog,*) '   Assuming the same as current setting: ', spinup_state
           end if
        end if
     end if ! read

     ! now compare the model and restart file spinup states, and either take the 
     ! model into spinup mode or out of it if they are not identical.
     ! taking model out of spinup mode requires multiplying each pool 
     ! by the associated AD factor.
     ! putting model into spinup mode requires dividing each pool 
     ! by the associated AD factor.
     ! only allow this to occur on first timestep of model run.
     
     if (flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
        if (spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
           if ( masterproc ) write(iulog,*) ' veg_cs_restart: taking deadwood pools out of AD spinup mode'
           exit_spinup = .true.
        else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
           if ( masterproc ) write(iulog,*) ' veg_cs_restart: taking deadwood pools into AD spinup mode'
           enter_spinup = .true.
        else
           call endrun(msg=' veg_ns_restart: error in entering/exiting spinup.  spinup_state ' &
                // ' != restart_file_spinup_state, but do not know what to do'//&
                errMsg(__FILE__, __LINE__))
        end if
        if (get_nstep() >= 2) then
           call endrun(msg=' veg_ns_restart: error in entering/exiting spinup - should occur only when nstep = 1'//&
                errMsg(__FILE__, __LINE__))
        endif
        do i = bounds%begp, bounds%endp
           if (exit_spinup) then 
              m_veg = spinup_mortality_factor
           else if (enter_spinup) then 
              m_veg = 1._r8 / spinup_mortality_factor
           end if
           this%deadstemn(i)  = this%deadstemn(i) * m_veg
           this%deadcrootn(i) = this%deadcrootn(i) * m_veg
        end do
     end if ! read

  end subroutine veg_ns_restart
  
  !-----------------------------------------------------------------------
  subroutine veg_ns_summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ns)
    !
    ! !DESCRIPTION:
    ! Vegetation-level nitrogen state summary calculations
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state)            :: this
    type(bounds_type)           , intent(in)    :: bounds          
    integer                     , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                     , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type (column_nitrogen_state), intent(inout) :: col_ns          ! column-level nitrogen state for p2c
    
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p             ! indices
    integer  :: fp              ! filter indices
    !-----------------------------------------------------------------------
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
       this%dispvegn(p) = &
            this%leafn(p)      + &
            this%frootn(p)     + &
            this%livestemn(p)  + &
            this%deadstemn(p)  + &
            this%livecrootn(p) + &
            this%deadcrootn(p)
       
      ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
      this%storvegn(p) = &
           this%leafn_storage(p)      + &
           this%frootn_storage(p)     + &
           this%livestemn_storage(p)  + &
           this%deadstemn_storage(p)  + &
           this%livecrootn_storage(p) + &
           this%deadcrootn_storage(p) + &
           this%leafn_xfer(p)         + &
           this%frootn_xfer(p)        + &
           this%livestemn_xfer(p)     + &
           this%deadstemn_xfer(p)     + &
           this%livecrootn_xfer(p)    + &
           this%deadcrootn_xfer(p)    + &
           this%npool(p)              + &
           this%retransn(p)

      if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
         this%dispvegn(p) = &
              this%dispvegn(p) + &
              this%grainn(p)

         this%storvegn(p) = &
              this%storvegn(p) + &
              this%grainn_storage(p)     + &
              this%grainn_xfer(p)
      end if

      ! total vegetation nitrogen (TOTVEGN)
      this%totvegn(p) = &
           this%dispvegn(p) + &
           this%storvegn(p)

      ! total pft-level carbon (add pft_ntrunc)
      this%totpftn(p) = &
           this%totvegn(p) + &
           this%ntrunc(p)
   end do ! filtered veg loop

   call p2c(bounds, num_soilc, filter_soilc, &
        this%plant_n_buffer(bounds%begp:bounds%endp), &
        col_ns%plant_n_buffer(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totvegn(bounds%begp:bounds%endp), &
        col_ns%totvegn(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totpftn(bounds%begp:bounds%endp), &
        col_ns%totpftn(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%cropseedn_deficit(bounds%begp:bounds%endp), &
        col_ns%cropseedn_deficit(bounds%begc:bounds%endc))
  
  end subroutine veg_ns_summary

  !-----------------------------------------------------------------------
  subroutine veg_ns_setvalues ( this, num_veg, filter_veg, value_veg)
    !
    ! !DESCRIPTION:
    ! Set nitrogen state variables
    !
    ! !ARGUMENTS:
    class (vegetation_nitrogen_state) :: this
    integer , intent(in) :: num_veg
    integer , intent(in) :: filter_veg(:)
    real(r8), intent(in) :: value_veg
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_veg
       i = filter_veg(fi)

       this%leafn(i)              = value_veg
       this%leafn_storage(i)      = value_veg
       this%leafn_xfer(i)         = value_veg
       this%frootn(i)             = value_veg
       this%frootn_storage(i)     = value_veg
       this%frootn_xfer(i)        = value_veg
       this%livestemn(i)          = value_veg
       this%livestemn_storage(i)  = value_veg
       this%livestemn_xfer(i)     = value_veg
       this%deadstemn(i)          = value_veg
       this%deadstemn_storage(i)  = value_veg
       this%deadstemn_xfer(i)     = value_veg
       this%livecrootn(i)         = value_veg
       this%livecrootn_storage(i) = value_veg
       this%livecrootn_xfer(i)    = value_veg
       this%deadcrootn(i)         = value_veg
       this%deadcrootn_storage(i) = value_veg
       this%deadcrootn_xfer(i)    = value_veg
       this%retransn(i)           = value_veg
       this%npool(i)              = value_veg
       this%ntrunc(i)             = value_veg
       this%dispvegn(i)           = value_veg
       this%storvegn(i)           = value_veg
       this%totvegn(i)            = value_veg
       this%totpftn(i)            = value_veg
    end do

    if ( crop_prog )then
       do fi = 1,num_veg
          i = filter_veg(fi)
          this%grainn(i)            = value_veg
          this%grainn_storage(i)    = value_veg
          this%grainn_xfer(i)       = value_veg
          this%cropseedn_deficit(i) = value_veg
       end do
    end if

  end subroutine veg_ns_setvalues

  !-----------------------------------------------------------------------
  subroutine veg_ns_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegn(p) = 0._r8
       this%storvegn(p) = 0._r8
       this%totvegn(p)  = 0._r8
       this%totpftn(p)  = 0._r8
    end do

  end subroutine veg_ns_zerodwt
  
  !------------------------------------------------------------------------
  subroutine veg_ns_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ns_clean
  
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
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_ps
    !-----------------------------------------------------------------------
    if (crop_prog) then
       this%grainp(begp:endp) = spval
       call hist_addfld1d (fname='GRAINP', units='gP/m^2', &
            avgflag='A', long_name='grain P', &
            ptr_patch=this%grainp, default='inactive')

       this%cropseedp_deficit(begp:endp) = spval
       call hist_addfld1d (fname='CROPSEEDP_DEFICIT', units='gP/m^2', &
            avgflag='A', long_name='P used for crop seed that needs to be repaid', &
            ptr_patch=this%cropseedp_deficit)
    end if

    this%leafp(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP', units='gP/m^2', &
         avgflag='A', long_name='leaf P', &
         ptr_patch=this%leafp)

    this%leafp_storage(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='leaf P storage', &
         ptr_patch=this%leafp_storage, default='inactive')

    this%leafp_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_XFER', units='gP/m^2', &
         avgflag='A', long_name='leaf P transfer', &
         ptr_patch=this%leafp_xfer, default='inactive')

    this%frootp(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP', units='gP/m^2', &
         avgflag='A', long_name='fine root P', &
         ptr_patch=this%frootp)

    this%frootp_storage(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='fine root P storage', &
         ptr_patch=this%frootp_storage, default='inactive')

    this%frootp_xfer(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='fine root P transfer', &
         ptr_patch=this%frootp_xfer, default='inactive')

    this%livestemp(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP', units='gP/m^2', &
         avgflag='A', long_name='live stem P', &
         ptr_patch=this%livestemp)

    this%livestemp_storage(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='live stem P storage', &
         ptr_patch=this%livestemp_storage, default='inactive')

    this%livestemp_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_XFER', units='gP/m^2', &
         avgflag='A', long_name='live stem P transfer', &
         ptr_patch=this%livestemp_xfer, default='inactive')

    this%deadstemp(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP', units='gP/m^2', &
         avgflag='A', long_name='dead stem P', &
         ptr_patch=this%deadstemp)

    this%deadstemp_storage(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='dead stem P storage', &
         ptr_patch=this%deadstemp_storage, default='inactive')

    this%deadstemp_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_XFER', units='gP/m^2', &
         avgflag='A', long_name='dead stem P transfer', &
         ptr_patch=this%deadstemp_xfer, default='inactive')

    this%livecrootp(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P', &
         ptr_patch=this%livecrootp)

    this%livecrootp_storage(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P storage', &
         ptr_patch=this%livecrootp_storage, default='inactive')

    this%livecrootp_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P transfer', &
         ptr_patch=this%livecrootp_xfer, default='inactive')

    this%deadcrootp(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P', &
         ptr_patch=this%deadcrootp)

    this%deadcrootp_storage(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P storage', &
         ptr_patch=this%deadcrootp_storage, default='inactive')

    this%deadcrootp_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P transfer', &
         ptr_patch=this%deadcrootp_xfer, default='inactive')

    this%retransp(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSP', units='gP/m^2', &
         avgflag='A', long_name='plant pool of retranslocated P', &
         ptr_patch=this%retransp)

    this%ppool(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL', units='gP/m^2', &
         avgflag='A', long_name='temporary plant P pool', &
         ptr_patch=this%ppool, default='inactive')

    this%ptrunc(begp:endp) = spval
    call hist_addfld1d (fname='PFT_PTRUNC', units='gP/m^2', &
         avgflag='A', long_name='pft-level sink for P truncation', &
         ptr_patch=this%ptrunc, default='inactive')

    this%dispvegp(begp:endp) = spval
    call hist_addfld1d (fname='DISPVEGP', units='gP/m^2', &
         avgflag='A', long_name='displayed vegetation phosphorus', &
         ptr_patch=this%dispvegp)

    this%storvegp(begp:endp) = spval
    call hist_addfld1d (fname='STORVEGP', units='gP/m^2', &
         avgflag='A', long_name='stored vegetation phosphorus', &
         ptr_patch=this%storvegp)

    this%totvegp(begp:endp) = spval
    call hist_addfld1d (fname='TOTVEGP', units='gP/m^2', &
         avgflag='A', long_name='total vegetation phosphorus', &
         ptr_patch=this%totvegp)

    this%totpftp(begp:endp) = spval
    call hist_addfld1d (fname='TOTPFTP', units='gP/m^2', &
         avgflag='A', long_name='total PFT-level phosphorus', &
         ptr_patch=this%totpftp)

    this%plant_p_buffer(begp:endp) = spval
    call hist_addfld1d (fname='PLANTP_BUFFER', units='gP/m^2', &
            avgflag='A', long_name='plant phosphorus stored as buffer', &
            ptr_col=this%plant_p_buffer,default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_ps
    !------------------------------------------------------------------------
    num_special_patch = 0
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do
    
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          if (veg_pp%itype(p) == noveg) then
             this%leafp(p) = 0._r8
             this%leafp_storage(p) = 0._r8
          else
             this%leafp(p)         = veg_cs%leafc(p)  / veg_vp%leafcp(veg_pp%itype(p))
             this%leafp_storage(p) = veg_cs%leafc_storage(p) / veg_vp%leafcp(veg_pp%itype(p))
          end if

          this%leafp_xfer(p)        = 0._r8
          if ( crop_prog )then
             this%grainp(p)            = 0._r8
             this%grainp_storage(p)    = 0._r8
             this%grainp_xfer(p)       = 0._r8
          end if
          this%cropseedp_deficit(p) = 0._r8
          this%frootp(p)            = 0._r8
          this%frootp_storage(p)    = 0._r8
          this%frootp_xfer(p)       = 0._r8
          this%livestemp(p)         = 0._r8
          this%livestemp_storage(p) = 0._r8
          this%livestemp_xfer(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
             this%deadstemp(p) = veg_cs%deadstemc(p) / veg_vp%deadwdcp(veg_pp%itype(p))
          else
             this%deadstemp(p) = 0._r8
          end if
          
          if (nu_com .ne. 'RD') then
              ! ECA competition calculate root NP uptake as a function of fine root biomass
              ! better to initialize root CNP pools with a non-zero value
              if (veg_pp%itype(p) .ne. noveg) then
                 this%frootp(p) = veg_cs%frootc(p) / veg_vp%frootcp(veg_pp%itype(p))
                 this%frootp_storage(p) = veg_cs%frootc_storage(p) / veg_vp%frootcp(veg_pp%itype(p))
              end if
          end if
           
          this%deadstemp_storage(p)  = 0._r8
          this%deadstemp_xfer(p)     = 0._r8
          this%livecrootp(p)         = 0._r8
          this%livecrootp_storage(p) = 0._r8
          this%livecrootp_xfer(p)    = 0._r8
          this%deadcrootp(p)         = 0._r8
          this%deadcrootp_storage(p) = 0._r8
          this%deadcrootp_xfer(p)    = 0._r8
          this%retransp(p)           = 0._r8
          this%ppool(p)              = 0._r8
          if (nstor(veg_pp%itype(p)) .gt. 1e-6_r8) then
              this%ppool(p)          = 1.0_r8
          end if
          this%ptrunc(p)             = 0._r8
          this%dispvegp(p)           = 0._r8
          this%storvegp(p)           = 0._r8
          this%totvegp(p)            = 0._r8
          this%totpftp(p)            = 0._r8
       end if
       this%plant_p_buffer(p)= 1.e-4_r8 
    end do
    
    call this%SetValues (num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8)
    
  end subroutine veg_ps_init
    
  !-----------------------------------------------------------------------
  subroutine veg_ps_restart ( this,  bounds, ncid, flag)
    !
    ! !DESCRIPTION: 
    ! Read/write vegetation-level phosphorus state restart data 
    !
    ! !ARGUMENTS:
    class (vegetation_phosphorus_state)        :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c,a,b,d
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: m, m_veg         ! multiplier for the exit_spinup code
    real(r8), pointer  :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    integer            :: itemp      ! temporary 
    integer , pointer  :: iptemp(:)  ! pointer to memory to be allocated
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    real(r8)           :: smax_c, ks_sorption_c
    real(r8)           :: rootfr(1:nlevdecomp)
    real(r8)           :: pinit_prof(1:nlevdecomp)
    real(r8)           :: rootfr_tot
    !------------------------------------------------------------------------
    call restartvar(ncid=ncid, flag=flag, varname='leafp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafp) 

    call restartvar(ncid=ncid, flag=flag, varname='leafp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafp_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='leafp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafp_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootp) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootp_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootp_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemp) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemp_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemp_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemp) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemp_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemp_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootp) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootp_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootp_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootp) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootp_storage) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootp_xfer) 

    call restartvar(ncid=ncid, flag=flag, varname='retransp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransp) 

    call restartvar(ncid=ncid, flag=flag, varname='ppool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ppool) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ptrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ptrunc) 

    call restartvar(ncid=ncid, flag=flag, varname='plant_p_buffer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_p_buffer)

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainp', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp)

       call restartvar(ncid=ncid, flag=flag,  varname='grainp_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P storage', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_storage)

       call restartvar(ncid=ncid, flag=flag,  varname='grainp_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P transfer', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_xfer)
    end if

    !--------------------------------
    ! Spinup state
    !--------------------------------
    ! Do nothing for write
    ! Note that the call to write spinup_state out was done in CNCarbonStateType and
    ! cannot be called again because it will try to define the variable twice
    ! when the flag below is set to define
    if (flag == 'read') then
       call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
            long_name='Spinup state of the model that wrote this restart file: ' &
            // ' 0 = normal model mode, 1 = AD spinup', units='', &
            interpinic_flag='copy', readvar=readvar,  data=idata)
       if (readvar) then
          restart_file_spinup_state = idata
       else
          ! assume, for sake of backwards compatibility, that if spinup_state is not in 
          ! the restart file then current model state is the same as prior model state
          restart_file_spinup_state = spinup_state
          if ( masterproc ) then
             write(iulog,*) ' WARNING!  Restart file does not contain info ' &
                  // ' on spinup state used to generate the restart file. '
             write(iulog,*) '   Assuming the same as current setting: ', spinup_state
          end if
       end if
    end if
    if (flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
       if (spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
          if ( masterproc ) write(iulog,*) ' Phosphorus State Restart: taking pools out of AD spinup mode'
          exit_spinup = .true.
       else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
          if ( masterproc ) write(iulog,*) ' Phosphorus State Restart: taking pools into AD spinup mode'
          enter_spinup = .true.
       else
          call endrun(msg=' Error in entering/exiting spinup.  spinup_state ' &
               // ' != restart_file_spinup_state, but do not know what to do'//&
               errMsg(__FILE__, __LINE__))
       end if
       if (get_nstep() >= 2) then
          call endrun(msg=' Error in entering/exiting spinup - should occur only when nstep = 1'//&
               errMsg(__FILE__, __LINE__))
       endif
       do i = bounds%begp, bounds%endp
          if (exit_spinup) then
             m_veg = spinup_mortality_factor
          else if (enter_spinup) then
             m_veg = 1._r8 / spinup_mortality_factor
          end if
          this%deadstemp(i)  = this%deadstemp(i) * m_veg
          this%deadcrootp(i) = this%deadcrootp(i) * m_veg
          if (nu_com == 'RD' .and. exit_spinup) then
             !Initialize plant P storage pool when exiting spinup from CN only mode 
             if (this%ppool(i) .lt. this%leafp(i)) then 
                this%ppool(i) = this%leafp(i)
             else if (this%ppool(i) .lt. this%leafp_storage(i)) then 
                this%ppool(i) = this%leafp_storage(i)
             end if
          end if
       end do
    end if
  
  end subroutine veg_ps_restart
  
  !-----------------------------------------------------------------------
  subroutine veg_ps_setvalues ( this, num_patch, filter_patch, value_patch)
    !
    ! !DESCRIPTION:
    ! Set phosphorus state variables, column-level
    !
    ! !ARGUMENTS:
    class (vegetation_phosphorus_state) :: this
    integer , intent(in)                :: num_patch
    integer , intent(in)                :: filter_patch(:)
    real(r8), intent(in)                :: value_patch
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k      ! indices
    !------------------------------------------------------------------------
    do fi = 1,num_patch
       i = filter_patch(fi)

       this%leafp(i)              = value_patch
       this%leafp_storage(i)      = value_patch
       this%leafp_xfer(i)         = value_patch
       this%frootp(i)             = value_patch
       this%frootp_storage(i)     = value_patch
       this%frootp_xfer(i)        = value_patch
       this%livestemp(i)          = value_patch
       this%livestemp_storage(i)  = value_patch
       this%livestemp_xfer(i)     = value_patch
       this%deadstemp(i)          = value_patch
       this%deadstemp_storage(i)  = value_patch
       this%deadstemp_xfer(i)     = value_patch
       this%livecrootp(i)         = value_patch
       this%livecrootp_storage(i) = value_patch
       this%livecrootp_xfer(i)    = value_patch
       this%deadcrootp(i)         = value_patch
       this%deadcrootp_storage(i) = value_patch
       this%deadcrootp_xfer(i)    = value_patch
       this%retransp(i)           = value_patch
       this%ppool(i)              = value_patch
       this%ptrunc(i)             = value_patch
       this%dispvegp(i)           = value_patch
       this%storvegp(i)           = value_patch
       this%totvegp(i)            = value_patch
       this%totpftp(i)            = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grainp(i)            = value_patch
          this%grainp_storage(i)    = value_patch
          this%grainp_xfer(i)       = value_patch
          this%cropseedp_deficit(i) = value_patch
       end do
    end if

  end subroutine veg_ps_setvalues
  
  !-----------------------------------------------------------------------
  subroutine veg_ps_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    type(bounds_type), intent(in)      :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegp(p) = 0._r8
       this%storvegp(p) = 0._r8
       this%totvegp(p)  = 0._r8
       this%totpftp(p)  = 0._r8
    end do

  end subroutine veg_ps_zerodwt

  !-----------------------------------------------------------------------
  subroutine veg_ps_summary (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ps)
    !
    ! !ARGUMENTS:
    class (vegetation_phosphorus_state) :: this
    type(bounds_type) , intent(in)      :: bounds  
    integer           , intent(in)      :: num_soilc       ! number of soil columns in filter
    integer           , intent(in)      :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in)      :: num_soilp       ! number of soil patches in filter
    integer           , intent(in)      :: filter_soilp(:) ! filter for soil patches
    type(column_phosphorus_state), intent(inout) :: col_ps
    !
    ! !LOCAL VARIABLES:
    integer  :: p        ! indices
    integer  :: fp       ! lake filter indices
    !-----------------------------------------------------------------------
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation phosphorus, excluding storage (DISPVEGN)
       this%dispvegp(p) = &
            this%leafp(p)      + &
            this%frootp(p)     + &
            this%livestemp(p)  + &
            this%deadstemp(p)  + &
            this%livecrootp(p) + &
            this%deadcrootp(p)
       
      ! stored vegetation phosphorus, including retranslocated N pool (STORVEGN)
      this%storvegp(p) = &
           this%leafp_storage(p)      + &
           this%frootp_storage(p)     + &
           this%livestemp_storage(p)  + &
           this%deadstemp_storage(p)  + &
           this%livecrootp_storage(p) + &
           this%deadcrootp_storage(p) + &
           this%leafp_xfer(p)         + &
           this%frootp_xfer(p)        + &
           this%livestemp_xfer(p)     + &
           this%deadstemp_xfer(p)     + &
           this%livecrootp_xfer(p)    + &
           this%deadcrootp_xfer(p)    + &
           this%ppool(p)              + &
           this%retransp(p)

      if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
         this%dispvegp(p) = &
              this%dispvegp(p) + &
              this%grainp(p)

         this%storvegp(p) = &
              this%storvegp(p) + &
              this%grainp_storage(p)     + &
              this%grainp_xfer(p)
      end if

      ! total vegetation phosphorus (TOTVEGN)
      this%totvegp(p) = &
           this%dispvegp(p) + &
           this%storvegp(p)

      ! total pft-level carbon (add pft_ntrunc)
      this%totpftp(p) = &
           this%totvegp(p) + &
           this%ptrunc(p)

   end do

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totvegp(bounds%begp:bounds%endp), &
        col_ps%totvegp(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totpftp(bounds%begp:bounds%endp), &
        col_ps%totpftp(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%cropseedp_deficit(bounds%begp:bounds%endp), &
        col_ps%cropseedp_deficit(bounds%begc:bounds%endc))

  end subroutine veg_ps_summary
  
  !------------------------------------------------------------------------
  subroutine veg_ps_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ps_clean

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

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_ef
    !-----------------------------------------------------------------------

    this%eflx_lwrad_net(begp:endp) = spval
    call hist_addfld1d (fname='FIRA', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net, c2l_scale_type='urbanf')

    this%eflx_lwrad_net_r(begp:endp) = spval 
    call hist_addfld1d (fname='FIRA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_r, set_spec=spval)

    this%eflx_lwrad_out(begp:endp) = spval 
    call hist_addfld1d (fname='FIRE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out, c2l_scale_type='urbanf')

    this%eflx_lwrad_out(begp:endp) = spval 
    call hist_addfld1d (fname='LWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling longwave radiation', &
         ptr_patch=this%eflx_lwrad_out, c2l_scale_type='urbanf', default='inactive')

    this%eflx_lwrad_out_r(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_R', units='W/m^2',  &
         avgflag='A', long_name='Rural emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_r, set_spec=spval)

    this%eflx_lh_vegt(begp:endp) = spval
    call hist_addfld1d (fname='FCTR', units='W/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_patch=this%eflx_lh_vegt, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_vege(begp:endp) = spval
    call hist_addfld1d (fname='FCEV', units='W/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_patch=this%eflx_lh_vege, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_grnd(begp:endp) = spval
    call hist_addfld1d (fname='FGEV', units='W/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_patch=this%eflx_lh_grnd, c2l_scale_type='urbanf') 

    this%eflx_sh_tot(begp:endp) = spval
    call hist_addfld1d (fname='FSH_NODYNLNDUSE', units='W/m^2',  &
         avgflag='A', long_name='sensible heat not including correction for land use change', &
         ptr_patch=this%eflx_sh_tot, c2l_scale_type='urbanf')

    this%eflx_sh_tot_r(begp:endp) = spval
    call hist_addfld1d (fname='FSH_R', units='W/m^2',  &
         avgflag='A', long_name='Rural sensible heat', &
         ptr_patch=this%eflx_sh_tot_r, set_spec=spval)

    this%eflx_sh_tot(begp:endp) = spval
    call hist_addfld1d (fname='Qh', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_patch=this%eflx_sh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot(begp:endp) = spval
    call hist_addfld1d (fname='Qle', units='W/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_patch=this%eflx_lh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
         avgflag='A', long_name='total latent heat flux [+ to atm]', &
         ptr_patch=this%eflx_lh_tot, c2l_scale_type='urbanf')

    this%eflx_lh_tot_r(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_R', units='W/m^2',  &
         avgflag='A', long_name='Rural total evaporation', &
         ptr_patch=this%eflx_lh_tot_r, set_spec=spval)

    this%eflx_soil_grnd(begp:endp) = spval
    call hist_addfld1d (fname='Qstor', units='W/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_patch=this%eflx_soil_grnd, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_sh_veg(begp:endp) = spval
    call hist_addfld1d (fname='FSH_V', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_patch=this%eflx_sh_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_sh_grnd(begp:endp) = spval
    call hist_addfld1d (fname='FSH_G', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_patch=this%eflx_sh_grnd, c2l_scale_type='urbanf')

    this%eflx_soil_grnd(begp:endp) = spval
    call hist_addfld1d (fname='FGR', units='W/m^2',  &
         avgflag='A', long_name='heat flux into soil/snow including snow melt and lake / snow light transmission', &
         ptr_patch=this%eflx_soil_grnd, c2l_scale_type='urbanf')

    this%eflx_soil_grnd_r(begp:endp) = spval
    call hist_addfld1d (fname='FGR_R', units='W/m^2',  &
         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt and snow light transmission', &
         ptr_patch=this%eflx_soil_grnd_r, set_spec=spval)

    this%eflx_lwrad_net_u(begp:endp) = spval
    call hist_addfld1d (fname='FIRA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_soil_grnd(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
         avgflag='A', long_name='soil heat flux [+ into soil]', &
         ptr_patch=this%eflx_soil_grnd, default='inactive', c2l_scale_type='urbanf')

    this%eflx_lwrad_out_u(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_U', units='W/m^2',  &
         avgflag='A', long_name='Urban emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_sh_tot_u(begp:endp) = spval
    call hist_addfld1d (fname='FSH_U', units='W/m^2',  &
         avgflag='A', long_name='Urban sensible heat', &
         ptr_patch=this%eflx_sh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_lh_tot_u(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_U', units='W/m^2',  &
         avgflag='A', long_name='Urban total evaporation', &
         ptr_patch=this%eflx_lh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_soil_grnd_u(begp:endp) = spval
    call hist_addfld1d (fname='FGR_U', units='W/m^2',  &
         avgflag='A', long_name='Urban heat flux into soil/snow including snow melt', &
         ptr_patch=this%eflx_soil_grnd_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%netrad(begp:endp) = spval
    call hist_addfld1d (fname='Rnet', units='W/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_patch=this%netrad, c2l_scale_type='urbanf', &
         default='inactive')

    if (use_cn) then
       this%dlrad(begp:endp) = spval
       call hist_addfld1d (fname='DLRAD', units='W/m^2', &
            avgflag='A', long_name='downward longwave radiation below the canopy', &
            ptr_patch=this%dlrad, default='inactive', c2l_scale_type='urbanf')

       this%ulrad(begp:endp) = spval
       call hist_addfld1d (fname='ULRAD', units='W/m^2', &
            avgflag='A', long_name='upward longwave radiation above the canopy', &
            ptr_patch=this%ulrad, default='inactive', c2l_scale_type='urbanf')

       this%cgrnd(begp:endp) = spval
       call hist_addfld1d (fname='CGRND', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil energy flux wrt to soil temp', &
            ptr_patch=this%cgrnd, default='inactive', c2l_scale_type='urbanf')

       this%cgrndl(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDL', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil latent heat flux wrt soil temp', &
            ptr_patch=this%cgrndl, default='inactive', c2l_scale_type='urbanf')

       this%cgrnds(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDS', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil sensible heat flux wrt soil temp', &
            ptr_patch=this%cgrnds, default='inactive', c2l_scale_type='urbanf')

       this%eflx_gnet(begp:endp) = spval
       call hist_addfld1d (fname='EFLX_GNET', units='W/m^2', &
            avgflag='A', long_name='net heat flux into ground', &
            ptr_patch=this%eflx_gnet, default='inactive', c2l_scale_type='urbanf')
    end if ! use_cn

    this%eflx_grnd_lake(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_GRND_LAKE', units='W/m^2', &
         avgflag='A', long_name='net heat flux into lake/snow surface, excluding light transmission', &
         ptr_patch=this%eflx_grnd_lake, set_nolake=spval)

    this%dgnetdT(begp:endp) = spval
    call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
         avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
         ptr_patch=this%dgnetdT, default='inactive', c2l_scale_type='urbanf')

    this%eflx_traffic(begp:endp) = spval
    call hist_addfld1d (fname='TRAFFICFLUX', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from urban traffic', &
         ptr_patch=this%eflx_traffic, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    this%eflx_wasteheat(begp:endp) = spval
    call hist_addfld1d (fname='WASTEHEAT', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from heating/cooling sources of urban waste heat', &
         ptr_patch=this%eflx_wasteheat, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_heat_from_ac(begp:endp) = spval
    call hist_addfld1d (fname='HEAT_FROM_AC', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux put into canyon due to heat removed from air conditioning', &
         ptr_patch=this%eflx_heat_from_ac, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_anthro(begp:endp) = spval
    call hist_addfld1d (fname='Qanth', units='W/m^2',  &
         avgflag='A', long_name='anthropogenic heat flux', &
         ptr_patch=this%eflx_anthro, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    this%taux(begp:endp) = spval
    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_patch=this%taux)

    this%tauy(begp:endp) = spval
    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_patch=this%tauy)

    this%errseb(begp:endp) = spval
    call hist_addfld1d (fname='ERRSEB',  units='W/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_patch=this%errseb)

    this%errsol(begp:endp) = spval
    call hist_addfld1d (fname='ERRSOL',  units='W/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_patch=this%errsol, set_urb=spval)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_ef
    !-----------------------------------------------------------------------
    
     do p = begp, endp 
        c = veg_pp%column(p)
        l = veg_pp%landunit(p)

        if (.not. lun_pp%urbpoi(l)) then ! non-urban
           this%eflx_lwrad_net_u(p)  = spval
           this%eflx_lwrad_out_u(p)  = spval
           this%eflx_lh_tot_u(p)     = spval
           this%eflx_sh_tot_u(p)     = spval
           this%eflx_soil_grnd_u(p)  = spval
           this%eflx_wasteheat(p)    = 0._r8
           this%eflx_heat_from_ac(p) = 0._r8
           this%eflx_traffic(p)      = 0._r8
           this%eflx_anthro(p)       = 0._r8
        end if

        this%eflx_lwrad_out(p) = sb * (col_es%t_grnd(c))**4
     end do
    
  end subroutine veg_ef_init
    
  !------------------------------------------------------------------------
  subroutine veg_ef_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation energy flux information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  & 
         dim1name='pft', &
         long_name='emitted infrared (longwave) radiation', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_lwrad_out)

  end subroutine veg_ef_restart

  !------------------------------------------------------------------------
  subroutine veg_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ef_clean
  
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
    allocate(this%qflx_dirct_rain        (begp:endp))             ; this%qflx_dirct_rain      (:)   = nan
    allocate(this%qflx_leafdrip          (begp:endp))             ; this%qflx_leafdrip        (:)   = nan
    allocate(this%qflx_ev_snow           (begp:endp))             ; this%qflx_ev_snow         (:)   = nan
    allocate(this%qflx_ev_soil           (begp:endp))             ; this%qflx_ev_soil         (:)   = nan
    allocate(this%qflx_ev_h2osfc         (begp:endp))             ; this%qflx_ev_h2osfc       (:)   = nan
    allocate(this%qflx_rootsoi_frac      (begp:endp,1:nlevgrnd))  ; this%qflx_rootsoi_frac    (:,:) = nan
    
    allocate(this%irrig_rate               (begp:endp))              ; this%irrig_rate               (:)   = nan
    allocate(this%qflx_irrig_patch         (begp:endp))              ; this%qflx_irrig_patch         (:)   = nan
    allocate(this%qflx_real_irrig_patch    (begp:endp))              ; this%qflx_real_irrig_patch    (:)   = nan
    allocate(this%qflx_grnd_irrig_patch    (begp:endp))              ; this%qflx_grnd_irrig_patch    (:)   = nan
    allocate(this%qflx_surf_irrig_patch    (begp:endp))              ; this%qflx_surf_irrig_patch    (:)   = nan
    allocate(this%qflx_supply_patch        (begp:endp))              ; this%qflx_supply_patch        (:)   = nan
    allocate(this%qflx_over_supply_patch   (begp:endp))              ; this%qflx_over_supply_patch   (:)   = nan 
    allocate(this%n_irrig_steps_left       (begp:endp))              ; this%n_irrig_steps_left       (:)   = 0
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_wf
    !-----------------------------------------------------------------------
    this%qflx_irrig_patch(begp:endp) = spval
    call hist_addfld1d (fname='QIRRIG_ORIG', units='mm/s', &
         avgflag='A', long_name='Original total irrigation water demand (surface + ground)', &
         ptr_patch=this%qflx_irrig_patch)
		 
    this%qflx_real_irrig_patch(begp:endp) = spval     ! real irrig
    call hist_addfld1d (fname='QIRRIG_REAL', units='mm/s', &
         avgflag='A', long_name='actual water added through irrigation (surface + ground)', &
         ptr_patch=this%qflx_real_irrig_patch)
         
    this%qflx_surf_irrig_patch(begp:endp) = spval    
    call hist_addfld1d (fname='QIRRIG_SURF', units='mm/s', &
         avgflag='A', long_name='Surface water irrigation', &
         ptr_patch=this%qflx_surf_irrig_patch)
    
    this%qflx_grnd_irrig_patch(begp:endp) = spval   
    call hist_addfld1d (fname='QIRRIG_GRND', units='mm/s', &
         avgflag='A', long_name='Groundwater irrigation', &
         ptr_patch=this%qflx_grnd_irrig_patch)	 

    this%qflx_prec_intr(begp:endp) = spval
    call hist_addfld1d (fname='QINTR', units='mm/s',  &
         avgflag='A', long_name='interception', &
         ptr_patch=this%qflx_prec_intr, set_lake=0._r8)

    this%qflx_dirct_rain(begp:endp) = spval
    call hist_addfld1d (fname='QWTRGH', units='mm/s',  &
         avgflag='A', long_name='direct rain throughfall', &
         ptr_patch=this%qflx_dirct_rain, c2l_scale_type='urbanf', default='inactive')

    this%qflx_leafdrip(begp:endp) = spval
    call hist_addfld1d (fname='QWDRIP', units='mm/s',  &
         avgflag='A', long_name='leaf rain drip', &
         ptr_patch=this%qflx_leafdrip, c2l_scale_type='urbanf', default='inactive')

    this%qflx_prec_grnd(begp:endp) = spval
    call hist_addfld1d (fname='QDRIP', units='mm/s',  &
         avgflag='A', long_name='throughfall', &
         ptr_patch=this%qflx_prec_grnd, c2l_scale_type='urbanf')

    this%qflx_evap_soi(begp:endp) = spval
    call hist_addfld1d (fname='QSOIL', units='mm/s',  &
         avgflag='A', long_name= 'Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)', &
         ptr_patch=this%qflx_evap_soi, c2l_scale_type='urbanf')

    this%qflx_evap_can(begp:endp) = spval
    call hist_addfld1d (fname='QVEGE', units='mm/s',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_patch=this%qflx_evap_can, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_tran_veg(begp:endp) = spval
    call hist_addfld1d (fname='QVEGT', units='mm/s',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_patch=this%qflx_tran_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snwcp_liq(begp:endp) = spval
    call hist_addfld1d (fname='QSNWCPLIQ', units='mm H2O/s', &
         avgflag='A', long_name='excess rainfall due to snow capping', &
         ptr_patch=this%qflx_snwcp_liq, c2l_scale_type='urbanf', default='inactive')

    if (use_cn) then
       this%qflx_rain_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_RAIN_GRND', units='mm H2O/s', &
            avgflag='A', long_name='rain on ground after interception', &
            ptr_patch=this%qflx_rain_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_snow_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_SNOW_GRND', units='mm H2O/s', &
            avgflag='A', long_name='snow on ground after interception', &
            ptr_patch=this%qflx_snow_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_evap_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_EVAP_GRND', units='mm H2O/s', &
            avgflag='A', long_name='ground surface evaporation', &
            ptr_patch=this%qflx_evap_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_evap_veg(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_EVAP_VEG', units='mm H2O/s', &
            avgflag='A', long_name='vegetation evaporation', &
            ptr_patch=this%qflx_evap_veg, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_evap_tot(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_EVAP_TOT', units='mm H2O/s', &
            avgflag='A', long_name='qflx_evap_soi + qflx_evap_can + qflx_tran_veg', &
            ptr_patch=this%qflx_evap_tot, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_dew_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_DEW_GRND', units='mm H2O/s', &
            avgflag='A', long_name='ground surface dew formation', &
            ptr_patch=this%qflx_dew_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_sub_snow(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_SUB_SNOW', units='mm H2O/s', &
            avgflag='A', long_name='sublimation rate from snow pack', &
            ptr_patch=this%qflx_sub_snow, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_dew_snow(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_DEW_SNOW', units='mm H2O/s', &
            avgflag='A', long_name='surface dew added to snow pacK', &
            ptr_patch=this%qflx_dew_snow, default='inactive', c2l_scale_type='urbanf')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_wf
    !-----------------------------------------------------------------------
    this%qflx_evap_grnd(begp:endp) = 0.0_r8
    this%qflx_dew_grnd (begp:endp) = 0.0_r8
    this%qflx_dew_snow (begp:endp) = 0.0_r8

    do p = begp, endp
       l = veg_pp%landunit(p)
       
       if (lun_pp%itype(l)==istsoil) then
          this%n_irrig_steps_left(p) = 0
          this%irrig_rate(p)         = 0.0_r8
       end if
    end do
    
  end subroutine veg_wf_init
    
  !------------------------------------------------------------------------
  subroutine veg_wf_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation water flux information to/from restart file.
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux)     :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: dimlen       ! dimension length
    integer :: nump_global  ! total number of pfts, globally
    integer :: err_code     ! error code
    logical :: do_io
    !-----------------------------------------------------------------------

    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'n_irrig_steps_left', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='n_irrig_steps_left', xtype=ncd_int,  &
            dim1name='pft', &
            long_name='number of irrigation time steps left', units='#', &
            interpinic_flag='interp', readvar=readvar, data=this%n_irrig_steps_left)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%n_irrig_steps_left = 0
    end if

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'irrig_rate', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='irrig_rate', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='irrigation rate', units='mm/s', &
            interpinic_flag='interp', readvar=readvar, data=this%irrig_rate)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%irrig_rate = 0.0_r8
    end if

  end subroutine veg_wf_restart

  !------------------------------------------------------------------------
  subroutine veg_wf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_wf_clean
  
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
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_cf
    !-----------------------------------------------------------------------
    if (use_fates) then
       ! no veg-level carbon flux history fields defined by host model

    else if (carbon_type == 'c12') then
       if (crop_prog) then
          this%grainc_to_food(begp:endp) = spval
          call hist_addfld1d (fname='GRAINC_TO_FOOD', units='gC/m^2/s', &
               avgflag='A', long_name='grain C to food', &
               ptr_patch=this%grainc_to_food, default='inactive')
       end if

       this%woodc_alloc(begp:endp) = spval
       call hist_addfld1d (fname='WOODC_ALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='wood C eallocation', &
            ptr_patch=this%woodc_alloc)

       this%woodc_loss(begp:endp) = spval
       call hist_addfld1d (fname='WOODC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='wood C loss', &
            ptr_patch=this%woodc_loss)

       this%leafc_loss(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C loss', &
            ptr_patch=this%leafc_loss)

       this%leafc_alloc(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_ALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C allocation', &
            ptr_patch=this%leafc_alloc)

       this%frootc_loss(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_LOSS', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C loss', &
            ptr_patch=this%frootc_loss)

       this%frootc_alloc(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_ALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C allocation', &
            ptr_patch=this%frootc_alloc)

       this%m_leafc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C mortality', &
            ptr_patch=this%m_leafc_to_litter, default='inactive')

       this%m_frootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C mortality', &
            ptr_patch=this%m_frootc_to_litter, default='inactive')

       this%m_leafc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C storage mortality', &
            ptr_patch=this%m_leafc_storage_to_litter, default='inactive')

       this%m_frootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C storage mortality', &
            ptr_patch=this%m_frootc_storage_to_litter, default='inactive')

       this%m_livestemc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C storage mortality', &
            ptr_patch=this%m_livestemc_storage_to_litter, default='inactive')

       this%m_deadstemc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C storage mortality', &
            ptr_patch=this%m_deadstemc_storage_to_litter, default='inactive')

       this%m_livecrootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C storage mortality', &
            ptr_patch=this%m_livecrootc_storage_to_litter, default='inactive')

       this%m_deadcrootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C storage mortality', &
            ptr_patch=this%m_deadcrootc_storage_to_litter, default='inactive')

       this%m_leafc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C transfer mortality', &
            ptr_patch=this%m_leafc_xfer_to_litter, default='inactive')

       this%m_frootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C transfer mortality', &
            ptr_patch=this%m_frootc_xfer_to_litter, default='inactive')

       this%m_livestemc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C transfer mortality', &
            ptr_patch=this%m_livestemc_xfer_to_litter, default='inactive')

       this%m_deadstemc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C transfer mortality', &
            ptr_patch=this%m_deadstemc_xfer_to_litter, default='inactive')

       this%m_livecrootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C transfer mortality', &
            ptr_patch=this%m_livecrootc_xfer_to_litter, default='inactive')

       this%m_deadcrootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C transfer mortality', &
            ptr_patch=this%m_deadcrootc_xfer_to_litter, default='inactive')

       this%m_livestemc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C mortality', &
            ptr_patch=this%m_livestemc_to_litter, default='inactive')

       this%m_deadstemc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C mortality', &
            ptr_patch=this%m_deadstemc_to_litter, default='inactive')

       this%m_livecrootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVECROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C mortality', &
            ptr_patch=this%m_livecrootc_to_litter, default='inactive')

       this%m_deadcrootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADCROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C mortality', &
            ptr_patch=this%m_deadcrootc_to_litter, default='inactive')

       this%m_gresp_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration storage mortality', &
            ptr_patch=this%m_gresp_storage_to_litter, default='inactive')

       this%m_gresp_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration transfer mortality', &
            ptr_patch=this%m_gresp_xfer_to_litter, default='inactive')

       this%m_leafc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C fire loss', &
            ptr_patch=this%m_leafc_to_fire, default='inactive')

       this%m_leafc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C storage fire loss', &
            ptr_patch=this%m_leafc_storage_to_fire, default='inactive')

       this%m_leafc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C transfer fire loss', &
            ptr_patch=this%m_leafc_xfer_to_fire, default='inactive')

       this%m_livestemc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C fire loss', &
            ptr_patch=this%m_livestemc_to_fire, default='inactive')

       this%m_livestemc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C storage fire loss', &
            ptr_patch=this%m_livestemc_storage_to_fire, default='inactive')

       this%m_livestemc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C transfer fire loss', &
            ptr_patch=this%m_livestemc_xfer_to_fire, default='inactive')

       this%m_deadstemc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C fire loss', &
            ptr_patch=this%m_deadstemc_to_fire, default='inactive')

       this%m_deadstemc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C storage fire loss', &
            ptr_patch=this%m_deadstemc_storage_to_fire, default='inactive')

       this%m_deadstemc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C transfer fire loss', &
            ptr_patch=this%m_deadstemc_xfer_to_fire, default='inactive')

       this%m_frootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C fire loss', &
            ptr_patch=this%m_frootc_to_fire, default='inactive')

       this%m_frootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C storage fire loss', &
            ptr_patch=this%m_frootc_storage_to_fire, default='inactive')

       this%m_frootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C transfer fire loss', &
            ptr_patch=this%m_frootc_xfer_to_fire, default='inactive')

       this%m_livecrootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVEROOTC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live root C fire loss', &
            ptr_patch=this%m_livecrootc_to_fire, default='inactive')

       this%m_livecrootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVEROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live root C storage fire loss', &
            ptr_patch=this%m_livecrootc_storage_to_fire, default='inactive')

       this%m_livecrootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVEROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live root C transfer fire loss', &
            ptr_patch=this%m_livecrootc_xfer_to_fire, default='inactive')

       this%m_deadcrootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADROOTC_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead root C fire loss', &
            ptr_patch=this%m_deadcrootc_to_fire, default='inactive')

       this%m_deadcrootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead root C storage fire loss', &
            ptr_patch=this%m_deadcrootc_storage_to_fire, default='inactive')

       this%m_deadcrootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead root C transfer fire loss', &
            ptr_patch=this%m_deadcrootc_xfer_to_fire, default='inactive')

       this%m_gresp_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_GRESP_STORAGE_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration storage fire loss', &
            ptr_patch=this%m_gresp_storage_to_fire, default='inactive')

       this%m_gresp_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_GRESP_XFER_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration transfer fire loss', &
            ptr_patch=this%m_gresp_xfer_to_fire, default='inactive')

       this%m_cpool_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_CPOOL_TO_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='cpool fire loss', &
            ptr_patch=this%m_cpool_to_fire, default='inactive')

       this%m_leafc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C fire mortality to litter', &
            ptr_patch=this%m_leafc_to_litter_fire, default='inactive')

       ! add by F. Li and S. Levis
       this%m_leafc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C fire mortality to litter', &
            ptr_patch=this%m_leafc_storage_to_litter_fire, default='inactive')

       this%m_leafc_xfer_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C transfer fire mortality to litter', &
            ptr_patch=this%m_leafc_xfer_to_litter_fire, default='inactive')

       this%m_livestemc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C fire mortality to litter', &
            ptr_patch=this%m_livestemc_to_litter_fire, default='inactive')

       this%m_livestemc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C storage fire mortality to litter', &
            ptr_patch=this%m_livestemc_storage_to_litter_fire, default='inactive')

       this%m_livestemc_xfer_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C transfer fire mortality to litter', &
            ptr_patch=this%m_livestemc_xfer_to_litter_fire, default='inactive')

       this%m_livestemc_to_deadstemc_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVESTEMC_TO_DEADSTEMC_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C fire mortality to dead stem C', &
            ptr_patch=this%m_livestemc_to_deadstemc_fire, default='inactive')

       this%m_deadstemc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C fire mortality to litter', &
            ptr_patch=this%m_deadstemc_to_litter_fire, default='inactive')

       this%m_deadstemc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C storage fire mortality to litter', &
            ptr_patch=this%m_deadstemc_storage_to_litter_fire, default='inactive')

       this%m_deadstemc_xfer_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C transfer fire mortality to litter', &
            ptr_patch=this%m_deadstemc_xfer_to_litter_fire, default='inactive')

       this%m_frootc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C fire mortality to litter', &
            ptr_patch=this%m_frootc_to_litter_fire, default='inactive')

       this%m_frootc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C storage fire mortality to litter', &
            ptr_patch=this%m_frootc_storage_to_litter_fire, default='inactive')

       this%m_frootc_xfer_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C transfer fire mortality to litter', &
            ptr_patch=this%m_frootc_xfer_to_litter_fire, default='inactive')

       this%m_livecrootc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVEROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live root C fire mortality to litter', &
            ptr_patch=this%m_livecrootc_to_litter_fire, default='inactive')

       this%m_livecrootc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVEROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live root C storage fire mortality to litter', &
            ptr_patch=this%m_livecrootc_storage_to_litter_fire, default='inactive')

       this%m_livecrootc_xfer_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVEROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live root C transfer fire mortality to litter', &
            ptr_patch=this%m_livecrootc_xfer_to_litter_fire, default='inactive')

       this%m_livecrootc_to_deadcrootc_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVEROOTC_TO_DEADROOTC_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live root C fire mortality to dead root C', &
            ptr_patch=this%m_livecrootc_to_deadcrootc_fire, default='inactive')

       this%m_deadcrootc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead root C fire mortality to litter', &
            ptr_patch=this%m_deadcrootc_to_litter_fire, default='inactive')

       this%m_deadcrootc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead root C storage fire mortality to litter', &
            ptr_patch=this%m_deadcrootc_storage_to_litter_fire, default='inactive')

       this%m_deadcrootc_xfer_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead root C transfer fire mortality to litter', &
            ptr_patch=this%m_deadcrootc_xfer_to_litter_fire, default='inactive')

       this%m_livecrootc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C fire mortality to litter', &
            ptr_patch=this%m_livecrootc_storage_to_litter_fire, default='inactive')

       this%m_deadcrootc_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C storage fire mortality to litter', &
            ptr_patch=this%m_deadcrootc_storage_to_litter_fire,  default='inactive')

       this%m_gresp_storage_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration storage fire mortality to litter', &
            ptr_patch=this%m_gresp_storage_to_litter_fire, default='inactive')

       this%m_gresp_xfer_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration transfer fire mortality to litter', &
            ptr_patch=this%m_gresp_xfer_to_litter_fire, default='inactive')   

      this%m_cpool_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='M_CPOOL_TO_LITTER_FIRE', units='gC/m^2/s', &
            avgflag='A', long_name='cpool fire mortality to litter', &
            ptr_patch=this%m_cpool_to_litter_fire, default='inactive')

       this%m_cpool_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='M_CPOOL_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='cpool mortality to litter', &
            ptr_patch=this%m_cpool_to_litter, default='inactive')

       this%leafc_xfer_to_leafc(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER_TO_LEAFC', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C growth from storage', &
            ptr_patch=this%leafc_xfer_to_leafc, default='inactive')

       this%frootc_xfer_to_frootc(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER_TO_FROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C growth from storage', &
            ptr_patch=this%frootc_xfer_to_frootc, default='inactive')

       this%livestemc_xfer_to_livestemc(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER_TO_LIVESTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C growth from storage', &
            ptr_patch=this%livestemc_xfer_to_livestemc, default='inactive')

       this%deadstemc_xfer_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER_TO_DEADSTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C growth from storage', &
            ptr_patch=this%deadstemc_xfer_to_deadstemc, default='inactive')

       this%livecrootc_xfer_to_livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER_TO_LIVECROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C growth from storage', &
            ptr_patch=this%livecrootc_xfer_to_livecrootc, default='inactive')

       this%deadcrootc_xfer_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER_TO_DEADCROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C growth from storage', &
            ptr_patch=this%deadcrootc_xfer_to_deadcrootc, default='inactive')

       this%leafc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C litterfall', &
            ptr_patch=this%leafc_to_litter, default='active')

       this%frootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C litterfall', &
            ptr_patch=this%frootc_to_litter, default='inactive')

       this%leaf_mr(begp:endp) = spval
       call hist_addfld1d (fname='LEAF_MR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf maintenance respiration', &
            ptr_patch=this%leaf_mr)

       this%froot_mr(begp:endp) = spval
       call hist_addfld1d (fname='FROOT_MR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root maintenance respiration', &
            ptr_patch=this%froot_mr, default='inactive')

       this%livestem_mr(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEM_MR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem maintenance respiration', &
            ptr_patch=this%livestem_mr, default='inactive')

       this%livecroot_mr(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOT_MR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root maintenance respiration', &
            ptr_patch=this%livecroot_mr, default='inactive')

       this%psnsun_to_cpool(begp:endp) = spval
       call hist_addfld1d (fname='PSNSUN_TO_CPOOL', units='gC/m^2/s', &
            avgflag='A', long_name='C fixation from sunlit canopy', &
            ptr_patch=this%psnsun_to_cpool)

       this%psnshade_to_cpool(begp:endp) = spval
       call hist_addfld1d (fname='PSNSHADE_TO_CPOOL', units='gC/m^2/s', &
            avgflag='A', long_name='C fixation from shaded canopy', &
            ptr_patch=this%psnshade_to_cpool)

       this%cpool_to_leafc(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_LEAFC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to leaf C', &
            ptr_patch=this%cpool_to_leafc, default='inactive')

       this%cpool_to_leafc_storage(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_LEAFC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to leaf C storage', &
            ptr_patch=this%cpool_to_leafc_storage, default='inactive')

       this%cpool_to_frootc(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_FROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to fine root C', &
            ptr_patch=this%cpool_to_frootc, default='inactive')

       this%cpool_to_frootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_FROOTC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to fine root C storage', &
            ptr_patch=this%cpool_to_frootc_storage, default='inactive')

       this%cpool_to_livestemc(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live stem C', &
            ptr_patch=this%cpool_to_livestemc, default='inactive')

       this%cpool_to_livestemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live stem C storage', &
            ptr_patch=this%cpool_to_livestemc_storage, default='inactive')

       this%cpool_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead stem C', &
            ptr_patch=this%cpool_to_deadstemc, default='inactive')

       this%cpool_to_deadstemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead stem C storage', &
            ptr_patch=this%cpool_to_deadstemc_storage, default='inactive')

       this%cpool_to_livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live coarse root C', &
            ptr_patch=this%cpool_to_livecrootc, default='inactive')

       this%cpool_to_livecrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to live coarse root C storage', &
            ptr_patch=this%cpool_to_livecrootc_storage, default='inactive')

       this%cpool_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead coarse root C', &
            ptr_patch=this%cpool_to_deadcrootc, default='inactive')

       this%cpool_to_deadcrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to dead coarse root C storage', &
            ptr_patch=this%cpool_to_deadcrootc_storage, default='inactive')

       this%cpool_to_gresp_storage(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_TO_GRESP_STORAGE', units='gC/m^2/s', &
            avgflag='A', long_name='allocation to growth respiration storage', &
            ptr_patch=this%cpool_to_gresp_storage, default='inactive')

       this%cpool_leaf_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_LEAF_GR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf growth respiration', &
            ptr_patch=this%cpool_leaf_gr, default='inactive')

       this%cpool_leaf_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_LEAF_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf growth respiration to storage', &
            ptr_patch=this%cpool_leaf_storage_gr, default='inactive')

       this%transfer_leaf_gr(begp:endp) = spval
       call hist_addfld1d (fname='TRANSFER_LEAF_GR', units='gC/m^2/s', &
            avgflag='A', long_name='leaf growth respiration from storage', &
            ptr_patch=this%transfer_leaf_gr, default='inactive')

       this%cpool_froot_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_FROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root growth respiration', &
            ptr_patch=this%cpool_froot_gr, default='inactive')

       this%cpool_froot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_FROOT_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root  growth respiration to storage', &
            ptr_patch=this%cpool_froot_storage_gr, default='inactive')

       this%transfer_froot_gr(begp:endp) = spval
       call hist_addfld1d (fname='TRANSFER_FROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='fine root  growth respiration from storage', &
            ptr_patch=this%transfer_froot_gr, default='inactive')

       this%cpool_livestem_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_LIVESTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem growth respiration', &
            ptr_patch=this%cpool_livestem_gr, default='inactive')

       this%cpool_livestem_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_LIVESTEM_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem growth respiration to storage', &
            ptr_patch=this%cpool_livestem_storage_gr, default='inactive')

       this%transfer_livestem_gr(begp:endp) = spval
       call hist_addfld1d (fname='TRANSFER_LIVESTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live stem growth respiration from storage', &
            ptr_patch=this%transfer_livestem_gr, default='inactive')

       this%cpool_deadstem_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_DEADSTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem growth respiration', &
            ptr_patch=this%cpool_deadstem_gr, default='inactive')

       this%cpool_deadstem_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_DEADSTEM_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem growth respiration to storage', &
            ptr_patch=this%cpool_deadstem_storage_gr, default='inactive')

       this%transfer_deadstem_gr(begp:endp) = spval
       call hist_addfld1d (fname='TRANSFER_DEADSTEM_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem growth respiration from storage', &
            ptr_patch=this%transfer_deadstem_gr, default='inactive')

       this%cpool_livecroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_LIVECROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root growth respiration', &
            ptr_patch=this%cpool_livecroot_gr, default='inactive')

       this%cpool_livecroot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_LIVECROOT_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root growth respiration to storage', &
            ptr_patch=this%cpool_livecroot_storage_gr, default='inactive')

       this%transfer_livecroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='TRANSFER_LIVECROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root growth respiration from storage', &
            ptr_patch=this%transfer_livecroot_gr, default='inactive')

       this%cpool_deadcroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_DEADCROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root growth respiration', &
            ptr_patch=this%cpool_deadcroot_gr, default='inactive')

       this%cpool_deadcroot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL_DEADCROOT_STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root growth respiration to storage', &
            ptr_patch=this%cpool_deadcroot_storage_gr, default='inactive')

       this%transfer_deadcroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='TRANSFER_DEADCROOT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root growth respiration from storage', &
            ptr_patch=this%transfer_deadcroot_gr, default='inactive')

       this%leafc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='leaf C shift storage to transfer', &
            ptr_patch=this%leafc_storage_to_xfer, default='inactive')

       this%frootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='fine root C shift storage to transfer', &
            ptr_patch=this%frootc_storage_to_xfer, default='inactive')

       this%livestemc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C shift storage to transfer', &
            ptr_patch=this%livestemc_storage_to_xfer, default='inactive')

       this%deadstemc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem C shift storage to transfer', &
            ptr_patch=this%deadstemc_storage_to_xfer, default='inactive')

       this%livecrootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C shift storage to transfer', &
            ptr_patch=this%livecrootc_storage_to_xfer, default='inactive')

       this%deadcrootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='dead coarse root C shift storage to transfer', &
            ptr_patch=this%deadcrootc_storage_to_xfer, default='inactive')

       this%gresp_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE_TO_XFER', units='gC/m^2/s', &
            avgflag='A', long_name='growth respiration shift storage to transfer', &
            ptr_patch=this%gresp_storage_to_xfer, default='inactive')

       this%livestemc_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_TO_DEADSTEMC', units='gC/m^2/s', &
            avgflag='A', long_name='live stem C turnover', &
            ptr_patch=this%livestemc_to_deadstemc, default='inactive')

       this%livecrootc_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_TO_DEADCROOTC', units='gC/m^2/s', &
            avgflag='A', long_name='live coarse root C turnover', &
            ptr_patch=this%livecrootc_to_deadcrootc, default='inactive')

       this%gpp(begp:endp) = spval
       call hist_addfld1d (fname='GPP', units='gC/m^2/s', &
            avgflag='A', long_name='gross primary production', &
            ptr_patch=this%gpp)

       this%gpp_before_downreg(begp:endp) = spval
       call hist_addfld1d (fname='INIT_GPP', units='gC/m^2/s', &
            avgflag='A', long_name='GPP flux before downregulation', &
            ptr_patch=this%gpp_before_downreg, default='inactive')

       this%mr(begp:endp) = spval
       call hist_addfld1d (fname='MR', units='gC/m^2/s', &
            avgflag='A', long_name='maintenance respiration', &
            ptr_patch=this%mr)

       this%current_gr(begp:endp) = spval
       call hist_addfld1d (fname='CURRENT_GR', units='gC/m^2/s', &
            avgflag='A', long_name='growth resp for new growth displayed in this timestep', &
            ptr_patch=this%current_gr, default='inactive')

       this%transfer_gr(begp:endp) = spval
       call hist_addfld1d (fname='TRANSFER_GR', units='gC/m^2/s', &
            avgflag='A', long_name='growth resp for transfer growth displayed in this timestep', &
            ptr_patch=this%transfer_gr, default='inactive')

       this%storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='STORAGE_GR', units='gC/m^2/s', &
            avgflag='A', long_name='growth resp for growth sent to storage for later display', &
            ptr_patch=this%storage_gr, default='inactive')

       this%gr(begp:endp) = spval
       call hist_addfld1d (fname='GR', units='gC/m^2/s', &
            avgflag='A', long_name='total growth respiration', &
            ptr_patch=this%gr)

       this%xr(begp:endp) = spval
       call hist_addfld1d (fname='XR', units='gC/m^2/s', &
            avgflag='A', long_name='total excess respiration', &
            ptr_patch=this%xr)

       this%ar(begp:endp) = spval
       call hist_addfld1d (fname='AR', units='gC/m^2/s', &
            avgflag='A', long_name='autotrophic respiration (MR + GR)', &
            ptr_patch=this%ar)

       this%rr(begp:endp) = spval
       call hist_addfld1d (fname='RR', units='gC/m^2/s', &
            avgflag='A', long_name='root respiration (fine root MR + total root GR)', &
            ptr_patch=this%rr)

       this%npp(begp:endp) = spval
       call hist_addfld1d (fname='NPP', units='gC/m^2/s', &
            avgflag='A', long_name='net primary production', &
            ptr_patch=this%npp)

       this%agnpp(begp:endp) = spval
       call hist_addfld1d (fname='AGNPP', units='gC/m^2/s', &
            avgflag='A', long_name='aboveground NPP', &
            ptr_patch=this%agnpp)

       this%bgnpp(begp:endp) = spval
       call hist_addfld1d (fname='BGNPP', units='gC/m^2/s', &
            avgflag='A', long_name='belowground NPP', &
            ptr_patch=this%bgnpp)

       this%agwdnpp(begp:endp) = spval
       call hist_addfld1d (fname='AGWDNPP', units='gC/m^2/s', &
            avgflag='A', long_name='aboveground wood NPP', &
            ptr_patch=this%agwdnpp)

       this%litfall(begp:endp) = spval
       call hist_addfld1d (fname='LITFALL', units='gC/m^2/s', &
            avgflag='A', long_name='litterfall (leaves and fine roots)', &
            ptr_patch=this%litfall)

       this%vegfire(begp:endp) = spval
       call hist_addfld1d (fname='VEGFIRE', units='gC/m^2/s', &
            avgflag='A', long_name='patch-level fire loss', &
            ptr_patch=this%vegfire, default='inactive')

       this%wood_harvestc(begp:endp) = spval
       call hist_addfld1d (fname='WOOD_HARVESTC', units='gC/m^2/s', &
            avgflag='A', long_name='wood harvest carbon (to product pools)', &
            ptr_patch=this%wood_harvestc)

       this%fire_closs(begp:endp) = spval
       call hist_addfld1d (fname='PFT_FIRE_CLOSS', units='gC/m^2/s', &
            avgflag='A', long_name='total patch-level fire C loss for non-peat fires outside land-type converted region', &
            ptr_patch=this%fire_closs)

       this%availc(begp:endp) = spval
       call hist_addfld1d (fname='AVAILC', units='gC/m^2/s', &
            avgflag='A', long_name='C flux available for allocation', &
            ptr_patch=this%availc, default='active')

       this%plant_calloc(begp:endp) = spval
       call hist_addfld1d (fname='PLANT_CALLOC', units='gC/m^2/s', &
            avgflag='A', long_name='total allocated C flux', &
            ptr_patch=this%plant_calloc, default='active')

       this%excess_cflux(begp:endp) = spval
       call hist_addfld1d (fname='EXCESS_CFLUX', units='gC/m^2/s', &
            avgflag='A', long_name='C flux not allocated due to downregulation', &
            ptr_patch=this%excess_cflux, default='inactive')

       this%prev_leafc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='PREV_LEAFC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='previous timestep leaf C litterfall flux', &
            ptr_patch=this%prev_leafc_to_litter, default='inactive')

       this%prev_frootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='PREV_FROOTC_TO_LITTER', units='gC/m^2/s', &
            avgflag='A', long_name='previous timestep froot C litterfall flux', &
            ptr_patch=this%prev_frootc_to_litter, default='inactive')

       this%xsmrpool_recover(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL_RECOVER', units='gC/m^2/s', &
            avgflag='A', long_name='C flux assigned to recovery of negative xsmrpool', &
            ptr_patch=this%xsmrpool_recover, default='inactive')

       if (nu_com .ne. 'RD' ) then
          this%allocation_leaf(begp:endp) = spval
          call hist_addfld1d (fname='allocation_leaf', units='', &
               avgflag='A', long_name='fraction of availc allocated to leaf', &
               ptr_patch=this%allocation_leaf)
          this%allocation_stem(begp:endp) = spval
          call hist_addfld1d (fname='allocation_stem', units='', &
               avgflag='A', long_name='fraction of availc allocated to stem', &
               ptr_patch=this%allocation_stem)
          this%allocation_froot(begp:endp) = spval
          call hist_addfld1d (fname='allocation_froot', units='', &
               avgflag='A', long_name='fraction of availc allocated to fine root', &
               ptr_patch=this%allocation_froot)
       end if

       this%crop_seedc_to_leaf(begp:endp) = spval
       call hist_addfld1d (fname='CROP_SEEDC_TO_LEAF', units='gC/m^2/s', &
            avgflag='A', long_name='crop seed source to leaf', &
            ptr_patch=this%crop_seedc_to_leaf, default='inactive')

       this%dwt_seedc_to_leaf(begp:endp) = spval
       call hist_addfld1d (fname='DWT_SEEDC_TO_LEAF_PATCH', units='gC/m^2/s', &
            avgflag='A', &
            long_name='patch-level seed source to patch-level leaf ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_seedc_to_leaf, default='inactive')

       this%dwt_seedc_to_deadstem(begp:endp) = spval
       call hist_addfld1d (fname='DWT_SEEDC_TO_DEADSTEM_PATCH', units='gC/m^2/s', &
            avgflag='A', &
            long_name='patch-level seed source to patch-level deadstem ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_seedc_to_deadstem, default='inactive')

       this%dwt_conv_cflux(begp:endp) = spval
       call hist_addfld1d (fname='DWT_CONV_CFLUX_PATCH', units='gC/m^2/s', &
            avgflag='A', &
            long_name='patch-level conversion C flux (immediate loss to atm) ' // &
            '(0 at all times except first timestep of year) ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_conv_cflux, default='inactive')

       this%dwt_prod10c_gain(begp:endp) = spval
       call hist_addfld1d (fname='DWT_PROD10C_GAIN_PATCH', units='gC/m^2/s', &
            avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
            ptr_patch=this%dwt_prod10c_gain, default='inactive')

       this%dwt_prod100c_gain(begp:endp) = spval
       call hist_addfld1d (fname='DWT_PROD100C_GAIN_PATCH', units='gC/m^2/s', &
            avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
            ptr_patch=this%dwt_prod100c_gain, default='inactive')

       this%annsum_npp(begp:endp) = spval
       call hist_addfld1d (fname='ANNSUM_NPP', units='gC/m^2/yr', &
            avgflag='A', long_name='annual sum of NPP', &
            ptr_patch=this%annsum_npp, default='inactive')

       ! end of C12 block

    else if ( carbon_type == 'c13') then
       this%m_leafc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LEAFC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C mortality', &
            ptr_patch=this%m_leafc_to_litter, default='inactive')

       this%m_frootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_FROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C mortality', &
            ptr_patch=this%m_frootc_to_litter, default='inactive')

       this%m_leafc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C storage mortality', &
            ptr_patch=this%m_leafc_storage_to_litter, default='inactive')

       this%m_frootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C storage mortality', &
            ptr_patch=this%m_frootc_storage_to_litter, default='inactive')

       this%m_livestemc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C storage mortality', &
            ptr_patch=this%m_livestemc_storage_to_litter, default='inactive')

       this%m_deadstemc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C storage mortality', &
            ptr_patch=this%m_deadstemc_storage_to_litter, default='inactive')

       this%m_livecrootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C storage mortality', &
            ptr_patch=this%m_livecrootc_storage_to_litter, default='inactive')

       this%m_deadcrootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C storage mortality', &
            ptr_patch=this%m_deadcrootc_storage_to_litter, default='inactive')

       this%m_leafc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C transfer mortality', &
            ptr_patch=this%m_leafc_xfer_to_litter, default='inactive')

       this%m_frootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C transfer mortality', &
            ptr_patch=this%m_frootc_xfer_to_litter, default='inactive')

       this%m_livestemc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C transfer mortality', &
            ptr_patch=this%m_livestemc_xfer_to_litter, default='inactive')

       this%m_deadstemc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C transfer mortality', &
            ptr_patch=this%m_deadstemc_xfer_to_litter, default='inactive')

       this%m_livecrootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C transfer mortality', &
            ptr_patch=this%m_livecrootc_xfer_to_litter, default='inactive')

       this%m_deadcrootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C transfer mortality', &
            ptr_patch=this%m_deadcrootc_xfer_to_litter, default='inactive')

       this%m_livestemc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C mortality', &
            ptr_patch=this%m_livestemc_to_litter, default='inactive')

       this%m_deadstemc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C mortality', &
            ptr_patch=this%m_deadstemc_to_litter, default='inactive')

       this%m_livecrootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C mortality', &
            ptr_patch=this%m_livecrootc_to_litter, default='inactive')

       this%m_deadcrootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C mortality', &
            ptr_patch=this%m_deadcrootc_to_litter, default='inactive')

       this%m_gresp_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration storage mortality', &
            ptr_patch=this%m_gresp_storage_to_litter, default='inactive')

       this%m_gresp_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration transfer mortality', &
            ptr_patch=this%m_gresp_xfer_to_litter, default='inactive')

       this%m_leafc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LEAFC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C fire loss', &
            ptr_patch=this%m_leafc_to_fire, default='inactive')

       this%m_frootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_FROOTC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C fire loss', &
            ptr_patch=this%m_frootc_to_fire, default='inactive')

       this%m_leafc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C storage fire loss', &
            ptr_patch=this%m_leafc_storage_to_fire, default='inactive')

       this%m_frootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C storage fire loss', &
            ptr_patch=this%m_frootc_storage_to_fire, default='inactive')

       this%m_livestemc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C storage fire loss', &
            ptr_patch=this%m_livestemc_storage_to_fire, default='inactive')

       this%m_deadstemc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C storage fire loss', &
            ptr_patch=this%m_deadstemc_storage_to_fire, default='inactive')

       this%m_livecrootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C storage fire loss', &
            ptr_patch=this%m_livecrootc_storage_to_fire, default='inactive')

       this%m_deadcrootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C storage fire loss', &
            ptr_patch=this%m_deadcrootc_storage_to_fire,  default='inactive')

       this%m_leafc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C transfer fire loss', &
            ptr_patch=this%m_leafc_xfer_to_fire, default='inactive')

       this%m_frootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C transfer fire loss', &
            ptr_patch=this%m_frootc_xfer_to_fire, default='inactive')

       this%m_livestemc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C transfer fire loss', &
            ptr_patch=this%m_livestemc_xfer_to_fire, default='inactive')

       this%m_deadstemc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C transfer fire loss', &
            ptr_patch=this%m_deadstemc_xfer_to_fire, default='inactive')

       this%m_livecrootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C transfer fire loss', &
            ptr_patch=this%m_livecrootc_xfer_to_fire, default='inactive')

       this%m_deadcrootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C transfer fire loss', &
            ptr_patch=this%m_deadcrootc_xfer_to_fire, default='inactive')

       this%m_livestemc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C fire loss', &
            ptr_patch=this%m_livestemc_to_fire, default='inactive')

       this%m_deadstemc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C fire loss', &
            ptr_patch=this%m_deadstemc_to_fire, default='inactive')

       this%m_deadstemc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C fire mortality to litter', &
            ptr_patch=this%m_deadstemc_to_litter_fire, default='inactive')

       this%m_livecrootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C fire loss', &
            ptr_patch=this%m_livecrootc_to_fire, default='inactive')

       this%m_deadcrootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C fire loss', &
            ptr_patch=this%m_deadcrootc_to_fire, default='inactive')

       this%m_deadcrootc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C fire mortality to litter', &
            ptr_patch=this%m_deadcrootc_to_litter_fire, default='inactive')

       this%m_gresp_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration storage fire loss', &
            ptr_patch=this%m_gresp_storage_to_fire, default='inactive')

       this%m_gresp_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_FIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration transfer fire loss', &
            ptr_patch=this%m_gresp_xfer_to_fire, default='inactive')

       this%leafc_xfer_to_leafc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_XFER_TO_LEAFC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C growth from storage', &
            ptr_patch=this%leafc_xfer_to_leafc, default='inactive')

       this%frootc_xfer_to_frootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_XFER_TO_FROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C growth from storage', &
            ptr_patch=this%frootc_xfer_to_frootc, default='inactive')

       this%livestemc_xfer_to_livestemc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C growth from storage', &
            ptr_patch=this%livestemc_xfer_to_livestemc, default='inactive')

       this%deadstemc_xfer_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C growth from storage', &
            ptr_patch=this%deadstemc_xfer_to_deadstemc, default='inactive')

       this%livecrootc_xfer_to_livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C growth from storage', &
            ptr_patch=this%livecrootc_xfer_to_livecrootc, default='inactive')

       this%deadcrootc_xfer_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C growth from storage', &
            ptr_patch=this%deadcrootc_xfer_to_deadcrootc, default='inactive')

       this%leafc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C litterfall', &
            ptr_patch=this%leafc_to_litter, default='inactive')

       this%frootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_TO_LITTER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C litterfall', &
            ptr_patch=this%frootc_to_litter, default='inactive')

       this%leaf_mr(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAF_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf maintenance respiration', &
            ptr_patch=this%leaf_mr, default='inactive')

       this%froot_mr(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOT_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root maintenance respiration', &
            ptr_patch=this%froot_mr, default='inactive')

       this%livestem_mr(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEM_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem maintenance respiration', &
            ptr_patch=this%livestem_mr, default='inactive')

       this%livecroot_mr(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOT_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root maintenance respiration', &
            ptr_patch=this%livecroot_mr, default='inactive')

       this%psnsun_to_cpool(begp:endp) = spval
       call hist_addfld1d (fname='C13_PSNSUN_TO_CPOOL', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 C fixation from sunlit canopy', &
            ptr_patch=this%psnsun_to_cpool)

       this%psnshade_to_cpool(begp:endp) = spval
       call hist_addfld1d (fname='C13_PSNSHADE_TO_CPOOL', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 C fixation from shaded canopy', &
            ptr_patch=this%psnshade_to_cpool)

       this%cpool_to_leafc(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to leaf C', &
            ptr_patch=this%cpool_to_leafc, default='inactive')

       this%cpool_to_leafc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to leaf C storage', &
            ptr_patch=this%cpool_to_leafc_storage, default='inactive')

       this%cpool_to_frootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to fine root C', &
            ptr_patch=this%cpool_to_frootc, default='inactive')

       this%cpool_to_frootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to fine root C storage', &
            ptr_patch=this%cpool_to_frootc_storage, default='inactive')

       this%cpool_to_livestemc(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live stem C', &
            ptr_patch=this%cpool_to_livestemc, default='inactive')

       this%cpool_to_livestemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live stem C storage', &
            ptr_patch=this%cpool_to_livestemc_storage, default='inactive')

       this%cpool_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead stem C', &
            ptr_patch=this%cpool_to_deadstemc, default='inactive')

       this%cpool_to_deadstemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead stem C storage', &
            ptr_patch=this%cpool_to_deadstemc_storage, default='inactive')

       this%cpool_to_livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live coarse root C', &
            ptr_patch=this%cpool_to_livecrootc, default='inactive')

       this%cpool_to_livecrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to live coarse root C storage', &
            ptr_patch=this%cpool_to_livecrootc_storage, default='inactive')

       this%cpool_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead coarse root C', &
            ptr_patch=this%cpool_to_deadcrootc, default='inactive')

       this%cpool_to_deadcrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to dead coarse root C storage', &
            ptr_patch=this%cpool_to_deadcrootc_storage, default='inactive')

       this%cpool_to_gresp_storage(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_TO_GRESP_STORAGE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 allocation to growth respiration storage', &
            ptr_patch=this%cpool_to_gresp_storage, default='inactive')

       this%cpool_leaf_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_LEAF_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf growth respiration', &
            ptr_patch=this%cpool_leaf_gr, default='inactive')

       this%cpool_leaf_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_LEAF_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf growth respiration to storage', &
            ptr_patch=this%cpool_leaf_storage_gr, default='inactive')

       this%transfer_leaf_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_TRANSFER_LEAF_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf growth respiration from storage', &
            ptr_patch=this%transfer_leaf_gr, default='inactive')

       this%cpool_froot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_FROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root growth respiration', &
            ptr_patch=this%cpool_froot_gr, default='inactive')

       this%cpool_froot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_FROOT_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root  growth respiration to storage', &
            ptr_patch=this%cpool_froot_storage_gr, default='inactive')

       this%transfer_froot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_TRANSFER_FROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root  growth respiration from storage', &
            ptr_patch=this%transfer_froot_gr, default='inactive')

       this%cpool_livestem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem growth respiration', &
            ptr_patch=this%cpool_livestem_gr, default='inactive')

       this%cpool_livestem_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem growth respiration to storage', &
            ptr_patch=this%cpool_livestem_storage_gr, default='inactive')

       this%transfer_livestem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_TRANSFER_LIVESTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem growth respiration from storage', &
            ptr_patch=this%transfer_livestem_gr, default='inactive')

       this%cpool_deadstem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem growth respiration', &
            ptr_patch=this%cpool_deadstem_gr, default='inactive')

       this%cpool_deadstem_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem growth respiration to storage', &
            ptr_patch=this%cpool_deadstem_storage_gr, default='inactive')

       this%transfer_deadstem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_TRANSFER_DEADSTEM_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem growth respiration from storage', &
            ptr_patch=this%transfer_deadstem_gr, default='inactive')

       this%cpool_livecroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root growth respiration', &
            ptr_patch=this%cpool_livecroot_gr, default='inactive')

       this%cpool_livecroot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root growth respiration to storage', &
            ptr_patch=this%cpool_livecroot_storage_gr, default='inactive')

       this%transfer_livecroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_TRANSFER_LIVECROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root growth respiration from storage', &
            ptr_patch=this%transfer_livecroot_gr, default='inactive')

       this%cpool_deadcroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root growth respiration', &
            ptr_patch=this%cpool_deadcroot_gr, default='inactive')

       this%cpool_deadcroot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root growth respiration to storage', &
            ptr_patch=this%cpool_deadcroot_storage_gr, default='inactive')

       this%transfer_deadcroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_TRANSFER_DEADCROOT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root growth respiration from storage', &
            ptr_patch=this%transfer_deadcroot_gr, default='inactive')

       this%leafc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 leaf C shift storage to transfer', &
            ptr_patch=this%leafc_storage_to_xfer, default='inactive')

       this%frootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 fine root C shift storage to transfer', &
            ptr_patch=this%frootc_storage_to_xfer, default='inactive')

       this%livestemc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C shift storage to transfer', &
            ptr_patch=this%livestemc_storage_to_xfer, default='inactive')

       this%deadstemc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem C shift storage to transfer', &
            ptr_patch=this%deadstemc_storage_to_xfer, default='inactive')

       this%livecrootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C shift storage to transfer', &
            ptr_patch=this%livecrootc_storage_to_xfer, default='inactive')

       this%deadcrootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead coarse root C shift storage to transfer', &
            ptr_patch=this%deadcrootc_storage_to_xfer, default='inactive')

       this%gresp_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_STORAGE_TO_XFER', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth respiration shift storage to transfer', &
            ptr_patch=this%gresp_storage_to_xfer, default='inactive')

       this%livestemc_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_TO_DEADSTEMC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live stem C turnover', &
            ptr_patch=this%livestemc_to_deadstemc, default='inactive')

       this%livecrootc_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_TO_DEADCROOTC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 live coarse root C turnover', &
            ptr_patch=this%livecrootc_to_deadcrootc, default='inactive')

       this%gpp(begp:endp) = spval
       call hist_addfld1d (fname='C13_GPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 gross primary production', &
            ptr_patch=this%gpp)

       this%mr(begp:endp) = spval
       call hist_addfld1d (fname='C13_MR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 maintenance respiration', &
            ptr_patch=this%mr)

       this%current_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_CURRENT_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth resp for new growth displayed in this timestep', &
            ptr_patch=this%current_gr, default='inactive')

       this%transfer_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_TRANSFER_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth resp for transfer growth displayed in this timestep', &
            ptr_patch=this%transfer_gr, default='inactive')

       this%storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_STORAGE_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 growth resp for growth sent to storage for later display', &
            ptr_patch=this%storage_gr, default='inactive')

       this%gr(begp:endp) = spval
       call hist_addfld1d (fname='C13_GR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total growth respiration', &
            ptr_patch=this%gr)

       this%ar(begp:endp) = spval
       call hist_addfld1d (fname='C13_AR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 autotrophic respiration (MR + GR)', &
            ptr_patch=this%ar)

       this%rr(begp:endp) = spval
       call hist_addfld1d (fname='C13_RR', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 root respiration (fine root MR + total root GR)', &
            ptr_patch=this%rr)

       this%npp(begp:endp) = spval
       call hist_addfld1d (fname='C13_NPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 net primary production', &
            ptr_patch=this%npp)

       this%agnpp(begp:endp) = spval
       call hist_addfld1d (fname='C13_AGNPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 aboveground NPP', &
            ptr_patch=this%agnpp)

       this%bgnpp(begp:endp) = spval
       call hist_addfld1d (fname='C13_BGNPP', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 belowground NPP', &
            ptr_patch=this%bgnpp)

       this%litfall(begp:endp) = spval
       call hist_addfld1d (fname='C13_LITFALL', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 litterfall (leaves and fine roots)', &
            ptr_patch=this%litfall, default='inactive')

       this%vegfire(begp:endp) = spval
       call hist_addfld1d (fname='C13_VEGFIRE', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 patch-level fire loss', &
            ptr_patch=this%vegfire, default='inactive')

       this%fire_closs(begp:endp) = spval
       call hist_addfld1d (fname='C13_PFT_FIRE_CLOSS', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 total patch-level fire C loss', &
            ptr_patch=this%fire_closs)

       this%crop_seedc_to_leaf(begp:endp) = spval
       call hist_addfld1d (fname='C13_CROP_SEEDC_TO_LEAF', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 crop seed source to leaf', &
            ptr_patch=this%crop_seedc_to_leaf, default='inactive')

       this%dwt_conv_cflux(begp:endp) = spval
       call hist_addfld1d (fname='C13_DWT_CONV_CFLUX_PATCH', units='gC13/m^2/s', &
            avgflag='A', &
            long_name='C13 patch-level conversion C flux (immediate loss to atm) ' // &
            '(0 at all times except first timestep of year) ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_conv_cflux, default='inactive')

       this%dwt_prod10c_gain(begp:endp) = spval
       call hist_addfld1d (fname='C13_DWT_PROD10C_GAIN_PATCH', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 landcover change-driven addition to 10-yr wood product pool', &
            ptr_patch=this%dwt_prod10c_gain, default='inactive')

       this%dwt_prod100c_gain(begp:endp) = spval
       call hist_addfld1d (fname='C13_DWT_PROD100C_GAIN_PATCH', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 landcover change-driven addition to 100-yr wood product pool', &
            ptr_patch=this%dwt_prod100c_gain, default='inactive')
            
       this%dwt_seedc_to_leaf(begp:endp) = spval
       call hist_addfld1d (fname='C13_DWT_SEEDC_TO_LEAF_PATCH', units='gC13/m^2/s', &
            avgflag='A', &
            long_name='patch-level C13 seed source to patch-level leaf ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_seedc_to_leaf, default='inactive')

       this%dwt_seedc_to_deadstem(begp:endp) = spval
       call hist_addfld1d (fname='C13_DWT_SEEDC_TO_DEADSTEM_PATCH', units='gC13/m^2/s', &
            avgflag='A', &
            long_name='patch-level C13 seed source to patch-level deadstem ' // &
           '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_seedc_to_deadstem, default='inactive')

       this%xsmrpool_c13ratio(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL_C13RATIO', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) ratio for xsmrpool', &
            ptr_patch=this%xsmrpool_c13ratio, default='inactive')
       
       ! end of C13 block     
    
    else if ( carbon_type == 'c14' ) then
       this%m_leafc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LEAFC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C mortality', &
            ptr_patch=this%m_leafc_to_litter, default='inactive')

       this%m_frootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_FROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C mortality', &
            ptr_patch=this%m_frootc_to_litter, default='inactive')

       this%m_leafc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LEAFC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C storage mortality', &
            ptr_patch=this%m_leafc_storage_to_litter, default='inactive')

       this%m_frootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_FROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C storage mortality', &
            ptr_patch=this%m_frootc_storage_to_litter, default='inactive')

       this%m_livestemc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C storage mortality', &
            ptr_patch=this%m_livestemc_storage_to_litter, default='inactive')

       this%m_deadstemc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C storage mortality', &
            ptr_patch=this%m_deadstemc_storage_to_litter, default='inactive')

       this%m_livecrootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C storage mortality', &
            ptr_patch=this%m_livecrootc_storage_to_litter, default='inactive')

       this%m_deadcrootc_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C storage mortality', &
            ptr_patch=this%m_deadcrootc_storage_to_litter, default='inactive')

       this%m_leafc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LEAFC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C transfer mortality', &
            ptr_patch=this%m_leafc_xfer_to_litter, default='inactive')

       this%m_frootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_FROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C transfer mortality', &
            ptr_patch=this%m_frootc_xfer_to_litter, default='inactive')

       this%m_livestemc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVESTEMC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C transfer mortality', &
            ptr_patch=this%m_livestemc_xfer_to_litter, default='inactive')

       this%m_deadstemc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADSTEMC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C transfer mortality', &
            ptr_patch=this%m_deadstemc_xfer_to_litter, default='inactive')

       this%m_livecrootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVECROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C transfer mortality', &
            ptr_patch=this%m_livecrootc_xfer_to_litter, default='inactive')

       this%m_deadcrootc_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADCROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C transfer mortality', &
            ptr_patch=this%m_deadcrootc_xfer_to_litter, default='inactive')

       this%m_livestemc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVESTEMC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C mortality', &
            ptr_patch=this%m_livestemc_to_litter, default='inactive')

       this%m_deadstemc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C mortality', &
            ptr_patch=this%m_deadstemc_to_litter, default='inactive')

       this%m_livecrootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVECROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C mortality', &
            ptr_patch=this%m_livecrootc_to_litter, default='inactive')

       this%m_deadcrootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C mortality', &
            ptr_patch=this%m_deadcrootc_to_litter, default='inactive')

       this%m_gresp_storage_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_GRESP_STORAGE_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration storage mortality', &
            ptr_patch=this%m_gresp_storage_to_litter, default='inactive')

       this%m_gresp_xfer_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_GRESP_XFER_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration transfer mortality', &
            ptr_patch=this%m_gresp_xfer_to_litter, default='inactive')

       this%m_leafc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LEAFC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C fire loss', &
            ptr_patch=this%m_leafc_to_fire, default='inactive')

       this%m_frootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_FROOTC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C fire loss', &
            ptr_patch=this%m_frootc_to_fire, default='inactive')

       this%m_leafc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LEAFC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C storage fire loss', &
            ptr_patch=this%m_leafc_storage_to_fire, default='inactive')

       this%m_frootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_FROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C storage fire loss', &
            ptr_patch=this%m_frootc_storage_to_fire, default='inactive')

       this%m_livestemc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C storage fire loss', &
            ptr_patch=this%m_livestemc_storage_to_fire, default='inactive')

       this%m_deadstemc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C storage fire loss', &
            ptr_patch=this%m_deadstemc_storage_to_fire, default='inactive')

       this%m_livecrootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C storage fire loss', &
            ptr_patch=this%m_livecrootc_storage_to_fire, default='inactive')

       this%m_deadcrootc_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C storage fire loss', &
            ptr_patch=this%m_deadcrootc_storage_to_fire,  default='inactive')

       this%m_leafc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LEAFC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C transfer fire loss', &
            ptr_patch=this%m_leafc_xfer_to_fire, default='inactive')

       this%m_frootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_FROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C transfer fire loss', &
            ptr_patch=this%m_frootc_xfer_to_fire, default='inactive')

       this%m_livestemc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVESTEMC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C transfer fire loss', &
            ptr_patch=this%m_livestemc_xfer_to_fire, default='inactive')

       this%m_deadstemc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADSTEMC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C transfer fire loss', &
            ptr_patch=this%m_deadstemc_xfer_to_fire, default='inactive')

       this%m_livecrootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVECROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C transfer fire loss', &
            ptr_patch=this%m_livecrootc_xfer_to_fire, default='inactive')

       this%m_deadcrootc_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADCROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C transfer fire loss', &
            ptr_patch=this%m_deadcrootc_xfer_to_fire, default='inactive')

       this%m_livestemc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVESTEMC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C fire loss', &
            ptr_patch=this%m_livestemc_to_fire, default='inactive')

       this%m_deadstemc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C fire loss', &
            ptr_patch=this%m_deadstemc_to_fire, default='inactive')

       this%m_deadstemc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_LITTER_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C fire mortality to litter', &
            ptr_patch=this%m_deadstemc_to_litter_fire, default='inactive')

       this%m_livecrootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_LIVECROOTC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C fire loss', &
            ptr_patch=this%m_livecrootc_to_fire, default='inactive')

       this%m_deadcrootc_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C fire loss', &
            ptr_patch=this%m_deadcrootc_to_fire, default='inactive')

       this%m_deadcrootc_to_litter_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_LITTER_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C fire mortality to litter', &
            ptr_patch=this%m_deadcrootc_to_litter_fire, default='inactive')

       this%m_gresp_storage_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_GRESP_STORAGE_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration storage fire loss', &
            ptr_patch=this%m_gresp_storage_to_fire, default='inactive')

       this%m_gresp_xfer_to_fire(begp:endp) = spval
       call hist_addfld1d (fname='C14_M_GRESP_XFER_TO_FIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration transfer fire loss', &
            ptr_patch=this%m_gresp_xfer_to_fire, default='inactive')

       this%leafc_xfer_to_leafc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_XFER_TO_LEAFC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C growth from storage', &
            ptr_patch=this%leafc_xfer_to_leafc, default='inactive')

       this%frootc_xfer_to_frootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_XFER_TO_FROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C growth from storage', &
            ptr_patch=this%frootc_xfer_to_frootc, default='inactive')

       this%livestemc_xfer_to_livestemc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C growth from storage', &
            ptr_patch=this%livestemc_xfer_to_livestemc, default='inactive')

       this%deadstemc_xfer_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C growth from storage', &
            ptr_patch=this%deadstemc_xfer_to_deadstemc, default='inactive')

       this%livecrootc_xfer_to_livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C growth from storage', &
            ptr_patch=this%livecrootc_xfer_to_livecrootc, default='inactive')

       this%deadcrootc_xfer_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C growth from storage', &
            ptr_patch=this%deadcrootc_xfer_to_deadcrootc, default='inactive')

       this%leafc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C litterfall', &
            ptr_patch=this%leafc_to_litter, default='inactive')

       this%frootc_to_litter(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_TO_LITTER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C litterfall', &
            ptr_patch=this%frootc_to_litter, default='inactive')

       this%leaf_mr(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAF_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf maintenance respiration', &
            ptr_patch=this%leaf_mr, default='inactive')

       this%froot_mr(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOT_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root maintenance respiration', &
            ptr_patch=this%froot_mr, default='inactive')

       this%livestem_mr(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEM_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem maintenance respiration', &
            ptr_patch=this%livestem_mr, default='inactive')

       this%livecroot_mr(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOT_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root maintenance respiration', &
            ptr_patch=this%livecroot_mr, default='inactive')

       this%psnsun_to_cpool(begp:endp) = spval
       call hist_addfld1d (fname='C14_PSNSUN_TO_CPOOL', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 C fixation from sunlit canopy', &
            ptr_patch=this%psnsun_to_cpool)

       this%psnshade_to_cpool(begp:endp) = spval
       call hist_addfld1d (fname='C14_PSNSHADE_TO_CPOOL', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 C fixation from shaded canopy', &
            ptr_patch=this%psnshade_to_cpool)

       this%cpool_to_leafc(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_LEAFC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to leaf C', &
            ptr_patch=this%cpool_to_leafc, default='inactive')

       this%cpool_to_leafc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_LEAFC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to leaf C storage', &
            ptr_patch=this%cpool_to_leafc_storage, default='inactive')

       this%cpool_to_frootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_FROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to fine root C', &
            ptr_patch=this%cpool_to_frootc, default='inactive')

       this%cpool_to_frootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_FROOTC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to fine root C storage', &
            ptr_patch=this%cpool_to_frootc_storage, default='inactive')

       this%cpool_to_livestemc(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVESTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live stem C', &
            ptr_patch=this%cpool_to_livestemc, default='inactive')

       this%cpool_to_livestemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVESTEMC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live stem C storage', &
            ptr_patch=this%cpool_to_livestemc_storage, default='inactive')

       this%cpool_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADSTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead stem C', &
            ptr_patch=this%cpool_to_deadstemc, default='inactive')

       this%cpool_to_deadstemc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADSTEMC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead stem C storage', &
            ptr_patch=this%cpool_to_deadstemc_storage, default='inactive')

       this%cpool_to_livecrootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVECROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live coarse root C', &
            ptr_patch=this%cpool_to_livecrootc, default='inactive')

       this%cpool_to_livecrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_LIVECROOTC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to live coarse root C storage', &
            ptr_patch=this%cpool_to_livecrootc_storage, default='inactive')

       this%cpool_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADCROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead coarse root C', &
            ptr_patch=this%cpool_to_deadcrootc, default='inactive')

       this%cpool_to_deadcrootc_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_DEADCROOTC_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to dead coarse root C storage', &
            ptr_patch=this%cpool_to_deadcrootc_storage, default='inactive')

       this%cpool_to_gresp_storage(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_TO_GRESP_STORAGE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 allocation to growth respiration storage', &
            ptr_patch=this%cpool_to_gresp_storage, default='inactive')

       this%cpool_leaf_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_LEAF_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf growth respiration', &
            ptr_patch=this%cpool_leaf_gr, default='inactive')

       this%cpool_leaf_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_LEAF_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf growth respiration to storage', &
            ptr_patch=this%cpool_leaf_storage_gr, default='inactive')

       this%transfer_leaf_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_TRANSFER_LEAF_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf growth respiration from storage', &
            ptr_patch=this%transfer_leaf_gr, default='inactive')

       this%cpool_froot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_FROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root growth respiration', &
            ptr_patch=this%cpool_froot_gr, default='inactive')

       this%cpool_froot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_FROOT_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root  growth respiration to storage', &
            ptr_patch=this%cpool_froot_storage_gr, default='inactive')

       this%transfer_froot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_TRANSFER_FROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root  growth respiration from storage', &
            ptr_patch=this%transfer_froot_gr, default='inactive')

       this%cpool_livestem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_LIVESTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem growth respiration', &
            ptr_patch=this%cpool_livestem_gr, default='inactive')

       this%cpool_livestem_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_LIVESTEM_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem growth respiration to storage', &
            ptr_patch=this%cpool_livestem_storage_gr, default='inactive')

       this%transfer_livestem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_TRANSFER_LIVESTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem growth respiration from storage', &
            ptr_patch=this%transfer_livestem_gr, default='inactive')

       this%cpool_deadstem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_DEADSTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem growth respiration', &
            ptr_patch=this%cpool_deadstem_gr, default='inactive')

       this%cpool_deadstem_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_DEADSTEM_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem growth respiration to storage', &
            ptr_patch=this%cpool_deadstem_storage_gr, default='inactive')

       this%transfer_deadstem_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_TRANSFER_DEADSTEM_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem growth respiration from storage', &
            ptr_patch=this%transfer_deadstem_gr, default='inactive')

       this%cpool_livecroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_LIVECROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root growth respiration', &
            ptr_patch=this%cpool_livecroot_gr, default='inactive')

       this%cpool_livecroot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_LIVECROOT_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root growth respiration to storage', &
            ptr_patch=this%cpool_livecroot_storage_gr, default='inactive')

       this%transfer_livecroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_TRANSFER_LIVECROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root growth respiration from storage', &
            ptr_patch=this%transfer_livecroot_gr, default='inactive')

       this%cpool_deadcroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_DEADCROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root growth respiration', &
            ptr_patch=this%cpool_deadcroot_gr, default='inactive')

       this%cpool_deadcroot_storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL_DEADCROOT_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root growth respiration to storage', &
            ptr_patch=this%cpool_deadcroot_storage_gr, default='inactive')

       this%transfer_deadcroot_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_TRANSFER_DEADCROOT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root growth respiration from storage', &
            ptr_patch=this%transfer_deadcroot_gr, default='inactive')

       this%leafc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 leaf C shift storage to transfer', &
            ptr_patch=this%leafc_storage_to_xfer, default='inactive')

       this%frootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 fine root C shift storage to transfer', &
            ptr_patch=this%frootc_storage_to_xfer, default='inactive')

       this%livestemc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C shift storage to transfer', &
            ptr_patch=this%livestemc_storage_to_xfer, default='inactive')

       this%deadstemc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem C shift storage to transfer', &
            ptr_patch=this%deadstemc_storage_to_xfer, default='inactive')

       this%livecrootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C shift storage to transfer', &
            ptr_patch=this%livecrootc_storage_to_xfer, default='inactive')

       this%deadcrootc_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead coarse root C shift storage to transfer', &
            ptr_patch=this%deadcrootc_storage_to_xfer, default='inactive')

       this%gresp_storage_to_xfer(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_STORAGE_TO_XFER', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth respiration shift storage to transfer', &
            ptr_patch=this%gresp_storage_to_xfer, default='inactive')

       this%livestemc_to_deadstemc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_TO_DEADSTEMC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live stem C turnover', &
            ptr_patch=this%livestemc_to_deadstemc, default='inactive')

       this%livecrootc_to_deadcrootc(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_TO_DEADCROOTC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 live coarse root C turnover', &
            ptr_patch=this%livecrootc_to_deadcrootc, default='inactive')

       this%gpp(begp:endp) = spval
       call hist_addfld1d (fname='C14_GPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 gross primary production', &
            ptr_patch=this%gpp)

       this%mr(begp:endp) = spval
       call hist_addfld1d (fname='C14_MR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 maintenance respiration', &
            ptr_patch=this%mr)

       this%current_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_CURRENT_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth resp for new growth displayed in this timestep', &
            ptr_patch=this%current_gr, default='inactive')

       this%transfer_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_TRANSFER_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth resp for transfer growth displayed in this timestep', &
            ptr_patch=this%transfer_gr, default='inactive')

       this%storage_gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_STORAGE_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 growth resp for growth sent to storage for later display', &
            ptr_patch=this%storage_gr, default='inactive')

       this%gr(begp:endp) = spval
       call hist_addfld1d (fname='C14_GR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total growth respiration', &
            ptr_patch=this%gr)

       this%ar(begp:endp) = spval
       call hist_addfld1d (fname='C14_AR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 autotrophic respiration (MR + GR)', &
            ptr_patch=this%ar)

       this%rr(begp:endp) = spval
       call hist_addfld1d (fname='C14_RR', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 root respiration (fine root MR + total root GR)', &
            ptr_patch=this%rr)

       this%npp(begp:endp) = spval
       call hist_addfld1d (fname='C14_NPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 net primary production', &
            ptr_patch=this%npp)

       this%agnpp(begp:endp) = spval
       call hist_addfld1d (fname='C14_AGNPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 aboveground NPP', &
            ptr_patch=this%agnpp)

       this%bgnpp(begp:endp) = spval
       call hist_addfld1d (fname='C14_BGNPP', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 belowground NPP', &
            ptr_patch=this%bgnpp)

       this%litfall(begp:endp) = spval
       call hist_addfld1d (fname='C14_LITFALL', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 litterfall (leaves and fine roots)', &
            ptr_patch=this%litfall, default='inactive')

       this%vegfire(begp:endp) = spval
       call hist_addfld1d (fname='C14_VEGFIRE', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 patch-level fire loss', &
            ptr_patch=this%vegfire, default='inactive')

       this%fire_closs(begp:endp) = spval
       call hist_addfld1d (fname='C14_PFT_FIRE_CLOSS', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 total patch-level fire C loss', &
            ptr_patch=this%fire_closs)
    
       this%crop_seedc_to_leaf(begp:endp) = spval
       call hist_addfld1d (fname='C14_CROP_SEEDC_TO_LEAF', units='gC13/m^2/s', &
            avgflag='A', long_name='C14 crop seed source to leaf', &
            ptr_patch=this%crop_seedc_to_leaf, default='inactive')

       this%dwt_conv_cflux(begp:endp) = spval
       call hist_addfld1d (fname='C14_DWT_CONV_CFLUX_PATCH', units='gC14/m^2/s', &
            avgflag='A', &
            long_name='C14 patch-level conversion C flux (immediate loss to atm) ' // &
            '(0 at all times except first timestep of year) ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_conv_cflux, default='inactive')

       this%dwt_prod10c_gain(begp:endp) = spval
       call hist_addfld1d (fname='C14_DWT_PROD10C_GAIN_PATCH', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 landcover change-driven addition to 10-yr wood product pool', &
            ptr_patch=this%dwt_prod10c_gain, default='inactive')

       this%dwt_prod100c_gain(begp:endp) = spval
       call hist_addfld1d (fname='C14_DWT_PROD100C_GAIN_PATCH', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 landcover change-driven addition to 100-yr wood product pool', &
            ptr_patch=this%dwt_prod100c_gain, default='inactive')

       this%dwt_seedc_to_leaf(begp:endp) = spval
       call hist_addfld1d (fname='C14_DWT_SEEDC_TO_LEAF_PATCH', units='gC14/m^2/s', &
            avgflag='A', &
            long_name='patch-level C14 seed source to patch-level leaf ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_seedc_to_leaf, default='inactive')

       this%dwt_seedc_to_deadstem(begp:endp) = spval
       call hist_addfld1d (fname='C14_DWT_SEEDC_TO_DEADSTEM_PATCH', units='gC14/m^2/s', &
            avgflag='A', &
            long_name='patch-level C14 seed source to patch-level deadstem ' // &
            '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
            ptr_patch=this%dwt_seedc_to_deadstem, default='inactive')

       ! end of C14 block
    endif

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_cf
    !-----------------------------------------------------------------------
    if (.not.use_fates) then
       do p = begp,endp
          l = veg_pp%landunit(p)

          this%gpp(p)                      = 0._r8
          this%gpp_before_downreg(p)       = 0._r8

          if (lun_pp%ifspecial(l)) then
             this%tempsum_npp(p)           = spval
             this%annsum_npp(p)            = spval
             this%availc(p)                = spval
             this%xsmrpool_recover(p)      = spval
             this%excess_cflux(p)          = spval
             this%plant_calloc(p)          = spval
             this%prev_leafc_to_litter(p)  = spval
             this%prev_frootc_to_litter(p) = spval
             if ( use_c13 ) then
                this%xsmrpool_c13ratio(p)  = spval
             endif
          end if
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
             this%tempsum_npp(p)           = 0._r8
             this%annsum_npp(p)            = 0._r8
             this%availc(p)                = 0._r8
             this%xsmrpool_recover(p)      = 0._r8
             this%excess_cflux(p)          = 0._r8
             this%prev_leafc_to_litter(p)  = 0._r8
             this%prev_frootc_to_litter(p) = 0._r8
             this%plant_calloc(p)          = 0._r8
          end if
       end do

    end if !(.not.use_fates)
       
    ! Set special patch filters, call SetValue
    num_special_patch = 0
    do p = begp,endp
       l = veg_pp%landunit(p)

       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do
    
    call this%SetValues (num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8)
   
  end subroutine veg_cf_init
    
  !-----------------------------------------------------------------------
  subroutine veg_cf_restart ( this, bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write restart data for vegetation carbon fluxes
    !
    ! !ARGUMENTS:
    class (vegetation_carbon_flux)    :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !------------------------------------------------------------------------

    ! -------------------------------------------
    ! None of these restarts are needed for FATES
    ! -------------------------------------------
    if (use_fates) return
    
    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer_to_grainc', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain C growth from storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer_to_grainc)

       call restartvar(ncid=ncid, flag=flag,  varname='livestemc_to_litter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='live stem C litterfall', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_to_litter)

       call restartvar(ncid=ncid, flag=flag,  varname='grainc_to_food', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain C to food', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_to_food)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_to_grainc', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain C', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_to_grainc)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_to_grainc_storage', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain C storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_to_grainc_storage)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_grain_gr', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain growth respiration', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_grain_gr)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_grain_storage_gr', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain growth respiration to storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_grain_storage_gr)

       call restartvar(ncid=ncid, flag=flag,  varname='transfer_grain_gr', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain growth respiration from storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%transfer_grain_gr)

       call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain C shift storage to transfer', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_storage_to_xfer)
    end if

    if (use_lch4 .or. use_betr) then
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Temp. Average AGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%tempavg_agnpp)
       
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Temp. Average BGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%tempavg_bgnpp)
       
       call restartvar(ncid=ncid, flag=flag, varname='annavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Ann. Average AGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%annavg_agnpp)
       
       call restartvar(ncid=ncid, flag=flag, varname='annavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Ann. Average BGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%annavg_bgnpp)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='gpp_pepv', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%gpp_before_downreg) 

    call restartvar(ncid=ncid, flag=flag, varname='availc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%availc) 

    call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_recover', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_recover) 

    call restartvar(ncid=ncid, flag=flag, varname='plant_calloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_calloc) 

    call restartvar(ncid=ncid, flag=flag, varname='excess_cflux', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%excess_cflux) 

    call restartvar(ncid=ncid, flag=flag, varname='prev_leafc_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prev_leafc_to_litter) 

    call restartvar(ncid=ncid, flag=flag, varname='prev_frootc_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prev_frootc_to_litter) 

    call restartvar(ncid=ncid, flag=flag, varname='tempsum_npp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempsum_npp) 
 
    call restartvar(ncid=ncid, flag=flag, varname='annsum_npp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_npp) 

  end subroutine veg_cf_restart
  
  !-----------------------------------------------------------------------
  subroutine veg_cf_summary(this, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, isotope, col_cf_input)
    !
    ! !DESCRIPTION:
    ! patch-level carbon flux summary calculations
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux)                 :: this
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    character(len=*)       , intent(in)    :: isotope   
    type(column_carbon_flux), intent(inout):: col_cf_input    ! receives p2c output
    !
    ! !LOCAL VARIABLES:
    integer  :: p,j,k,l       ! indices
    integer  :: fp            ! lake filter indices
    !-----------------------------------------------------------------------
    
    if (use_fates) return

    ! patch loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! calculate pft-level summary carbon fluxes and states

       ! gross primary production (GPP)
       this%gpp(p) = &
            this%psnsun_to_cpool(p) + &
            this%psnshade_to_cpool(p)

       ! maintenance respiration (MR)
       this%leaf_mr(p)      = this%leaf_curmr(p)      + this%leaf_xsmr(p)
       this%froot_mr(p)     = this%froot_curmr(p)     + this%froot_xsmr(p)
       this%livestem_mr(p)  = this%livestem_curmr(p)  + this%livestem_xsmr(p)
       this%livecroot_mr(p) = this%livecroot_curmr(p) + this%livecroot_xsmr(p)

       this%mr(p)  = &
            this%leaf_mr(p)     + &
            this%froot_mr(p)    + &
            this%livestem_mr(p) + &
            this%livecroot_mr(p)

       ! growth respiration (GR)
       ! current GR is respired this time step for new growth displayed in this timestep
       this%current_gr(p) = &
            this%cpool_leaf_gr(p)      + &
            this%cpool_froot_gr(p)     + &
            this%cpool_livestem_gr(p)  + &
            this%cpool_deadstem_gr(p)  + &
            this%cpool_livecroot_gr(p) + &
            this%cpool_deadcroot_gr(p)

       ! transfer GR is respired this time step for transfer growth displayed in this timestep
       this%transfer_gr(p) = &
            this%transfer_leaf_gr(p)      + &
            this%transfer_froot_gr(p)     + &
            this%transfer_livestem_gr(p)  + &
            this%transfer_deadstem_gr(p)  + &
            this%transfer_livecroot_gr(p) + &
            this%transfer_deadcroot_gr(p)

       ! storage GR is respired this time step for growth sent to storage for later display
       this%storage_gr(p) = &
            this%cpool_leaf_storage_gr(p)      + &
            this%cpool_froot_storage_gr(p)     + &
            this%cpool_livestem_storage_gr(p)  + &
            this%cpool_deadstem_storage_gr(p)  + &
            this%cpool_livecroot_storage_gr(p) + &
            this%cpool_deadcroot_storage_gr(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%mr(p) = &
               this%mr(p) + &
               this%grain_mr(p)

          this%current_gr(p) = &
               this%current_gr(p) + &
               this%cpool_grain_gr(p)

          this%transfer_gr(p) = &
               this%transfer_gr(p) + &
               this%transfer_grain_gr(p)

          this%storage_gr(p) = &
               this%storage_gr(p) + &
               this%cpool_grain_storage_gr(p)
       end if

       ! GR is the sum of current + transfer + storage GR
       this%gr(p) = &
            this%current_gr(p)  + &
            this%transfer_gr(p) + &
            this%storage_gr(p)

       ! autotrophic respiration (AR)
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%ar(p) = &
               this%mr(p) + &
               this%gr(p) + &
               this%xr(p) + &
               this%xsmrpool_to_atm(p) ! xsmr... is -ve (slevis)
          if (nu_com .ne. 'RD' ) then
             this%ar(p) = this%ar(p) + &
                  this%xsmrpool_turnover(p)
          end if
       else
          this%ar(p) = &
               this%mr(p) + &
               this%gr(p) + &
               this%xr(p)
          if (nu_com .ne. 'RD' ) then
             this%ar(p) = this%ar(p) + &
                  this%xsmrpool_turnover(p)
          end if
       end if



       ! net primary production (NPP)
       this%npp(p) = &
            this%gpp(p) - &
            this%ar(p)

       ! update the annual NPP accumulator, for use in allocation code 
       if (trim(isotope) == 'bulk') then      
          this%tempsum_npp(p) = &
               this%tempsum_npp(p) + &
               this%npp(p)
       end if

       ! litterfall (LITFALL)

       this%litfall(p) = &
            this%leafc_to_litter(p)                     + &
            this%frootc_to_litter(p)                    + &
            this%m_leafc_to_litter(p)                   + &
            this%m_leafc_storage_to_litter(p)           + &
            this%m_leafc_xfer_to_litter(p)              + &
            this%m_frootc_to_litter(p)                  + &
            this%m_frootc_storage_to_litter(p)          + &
            this%m_frootc_xfer_to_litter(p)             + &
            this%m_livestemc_to_litter(p)               + &
            this%m_livestemc_storage_to_litter(p)       + &
            this%m_livestemc_xfer_to_litter(p)          + &
            this%m_deadstemc_to_litter(p)               + &
            this%m_deadstemc_storage_to_litter(p)       + &
            this%m_deadstemc_xfer_to_litter(p)          + &
            this%m_livecrootc_to_litter(p)              + &
            this%m_livecrootc_storage_to_litter(p)      + &
            this%m_livecrootc_xfer_to_litter(p)         + &
            this%m_deadcrootc_to_litter(p)              + &
            this%m_deadcrootc_storage_to_litter(p)      + &
            this%m_deadcrootc_xfer_to_litter(p)         + &
            this%m_gresp_storage_to_litter(p)           + &
            this%m_gresp_xfer_to_litter(p)              + &
            
            this%m_leafc_to_litter_fire(p)              + &
            this%m_leafc_storage_to_litter_fire(p)      + &
            this%m_leafc_xfer_to_litter_fire(p)         + &
            this%m_livestemc_to_litter_fire(p)          + &
            this%m_livestemc_storage_to_litter_fire(p)  + &
            this%m_livestemc_xfer_to_litter_fire(p)     + &
            this%m_deadstemc_to_litter_fire(p)          + &
            this%m_deadstemc_storage_to_litter_fire(p)  + &
            this%m_deadstemc_xfer_to_litter_fire(p)     + &
            this%m_frootc_to_litter_fire(p)             + &
            this%m_frootc_storage_to_litter_fire(p)     + &
            this%m_frootc_xfer_to_litter_fire(p)        + &
            this%m_livecrootc_to_litter_fire(p)         + &
            this%m_livecrootc_storage_to_litter_fire(p) + &
            this%m_livecrootc_xfer_to_litter_fire(p)    + &
            this%m_deadcrootc_to_litter_fire(p)         + &
            this%m_deadcrootc_storage_to_litter_fire(p) + &
            this%m_deadcrootc_xfer_to_litter_fire(p)    + &
            this%m_gresp_storage_to_litter_fire(p)      + &
            this%m_gresp_xfer_to_litter_fire(p)         + &
            
            this%hrv_leafc_to_litter(p)                 + &
            this%hrv_leafc_storage_to_litter(p)         + &
            this%hrv_leafc_xfer_to_litter(p)            + &
            this%hrv_frootc_to_litter(p)                + &
            this%hrv_frootc_storage_to_litter(p)        + &
            this%hrv_frootc_xfer_to_litter(p)           + &
            this%hrv_livestemc_to_litter(p)             + &
            this%hrv_livestemc_storage_to_litter(p)     + &
            this%hrv_livestemc_xfer_to_litter(p)        + &
            this%hrv_deadstemc_storage_to_litter(p)     + &
            this%hrv_deadstemc_xfer_to_litter(p)        + &
            this%hrv_livecrootc_to_litter(p)            + &
            this%hrv_livecrootc_storage_to_litter(p)    + &
            this%hrv_livecrootc_xfer_to_litter(p)       + &
            this%hrv_deadcrootc_to_litter(p)            + &
            this%hrv_deadcrootc_storage_to_litter(p)    + &
            this%hrv_deadcrootc_xfer_to_litter(p)       + &
            this%hrv_gresp_storage_to_litter(p)         + &
            this%hrv_gresp_xfer_to_litter(p)            + &
            this%hrv_cpool_to_litter(p)

       ! patch-level fire losses (VEGFIRE)
       this%vegfire(p) = 0._r8

       ! patch-level wood harvest
       this%wood_harvestc(p) = &
            this%hrv_deadstemc_to_prod10c(p) + &
            this%hrv_deadstemc_to_prod100c(p)
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%wood_harvestc(p) = &
               this%wood_harvestc(p) + &
               this%hrv_cropc_to_prod1c(p)
       end if

       ! patch-level carbon losses to fire changed by F. Li and S. Levis
       this%fire_closs(p) = &
            this%m_leafc_to_fire(p)                + &
            this%m_leafc_storage_to_fire(p)        + &
            this%m_leafc_xfer_to_fire(p)           + &
            this%m_frootc_to_fire(p)               + &
            this%m_frootc_storage_to_fire(p)       + &
            this%m_frootc_xfer_to_fire(p)          + &
            this%m_livestemc_to_fire(p)            + &
            this%m_livestemc_storage_to_fire(p)    + &
            this%m_livestemc_xfer_to_fire(p)       + &
            this%m_deadstemc_to_fire(p)            + &
            this%m_deadstemc_storage_to_fire(p)    + &
            this%m_deadstemc_xfer_to_fire(p)       + &
            this%m_livecrootc_to_fire(p)           + &
            this%m_livecrootc_storage_to_fire(p)   + &
            this%m_livecrootc_xfer_to_fire(p)      + &
            this%m_deadcrootc_to_fire(p)           + &
            this%m_deadcrootc_storage_to_fire(p)   + &
            this%m_deadcrootc_xfer_to_fire(p)      + &
            this%m_gresp_storage_to_fire(p)        + &
            this%m_gresp_xfer_to_fire(p)           + &
            this%m_cpool_to_fire(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%litfall(p) =                  &
               this%litfall(p)             + &
               this%livestemc_to_litter(p) + &
               this%grainc_to_food(p)
       end if

       ! new summary variables for CLAMP

       ! (FROOTC_ALLOC) - fine root C allocation
       this%frootc_alloc(p) = &
            this%frootc_xfer_to_frootc(p)    + &
            this%cpool_to_frootc(p)     

       ! (FROOTC_LOSS) - fine root C loss changed by F. Li and S. Levis
       this%frootc_loss(p) = &
            this%m_frootc_to_litter(p)       + &
            this%m_frootc_to_fire(p)         + &
            this%m_frootc_to_litter_fire(p)  + &
            this%hrv_frootc_to_litter(p)     + &
            this%frootc_to_litter(p)

       ! (LEAFC_ALLOC) - leaf C allocation
       this%leafc_alloc(p) = &
            this%leafc_xfer_to_leafc(p)    + &
            this%cpool_to_leafc(p)     

       ! (LEAFC_LOSS) - leaf C loss changed by F. Li and S. Levis
       this%leafc_loss(p) = &
            this%m_leafc_to_litter(p)      + &
            this%m_leafc_to_fire(p)        + &
            this%m_leafc_to_litter_fire(p) + &
            this%hrv_leafc_to_litter(p)    + &
            this%leafc_to_litter(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%leafc_loss(p) = &
               this%leafc_loss(p) + &
               this%hrv_leafc_to_prod1c(p)
       end if


       ! (WOODC_ALLOC) - wood C allocation
       this%woodc_alloc(p) = &
            this%livestemc_xfer_to_livestemc(p)   + &
            this%deadstemc_xfer_to_deadstemc(p)   + &
            this%livecrootc_xfer_to_livecrootc(p) + &
            this%deadcrootc_xfer_to_deadcrootc(p) + &
            this%cpool_to_livestemc(p)            + &
            this%cpool_to_deadstemc(p)            + &
            this%cpool_to_livecrootc(p)           + &
            this%cpool_to_deadcrootc(p)

       ! (WOODC_LOSS) - wood C loss
       this%woodc_loss(p) = &
            this%m_livestemc_to_litter(p)            + &
            this%m_deadstemc_to_litter(p)            + &
            this%m_livecrootc_to_litter(p)           + &
            this%m_deadcrootc_to_litter(p)           + &
            this%m_livestemc_to_fire(p)              + &
            this%m_deadstemc_to_fire(p)              + &
            this%m_livecrootc_to_fire(p)             + &
            this%m_deadcrootc_to_fire(p)             + &
            this%hrv_livestemc_to_litter(p)          + &
            this%hrv_livestemc_storage_to_litter(p)  + &
            this%hrv_livestemc_xfer_to_litter(p)     + &
            this%hrv_deadstemc_to_prod10c(p)         + &
            this%hrv_deadstemc_to_prod100c(p)        + &
            this%hrv_deadstemc_storage_to_litter(p)  + &
            this%hrv_deadstemc_xfer_to_litter(p)     + &
            this%hrv_livecrootc_to_litter(p)         + &
            this%hrv_livecrootc_storage_to_litter(p) + &
            this%hrv_livecrootc_xfer_to_litter(p)    + &
            this%hrv_deadcrootc_to_litter(p)         + &
            this%hrv_deadcrootc_storage_to_litter(p) + &
            this%hrv_deadcrootc_xfer_to_litter(p)   
       ! putting the harvested crop stem and grain in the wood loss bdrewniak
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%woodc_loss(p) = &
               this%woodc_loss(p) + &
               this%hrv_grainc_to_prod1c(p) + &
               this%hrv_livestemc_to_prod1c(p)
       end if
    end do  ! end of patches loop

    ! use p2c routine to get selected column-average patch-level fluxes and states
    call p2c(bounds, num_soilc, filter_soilc, &
         this%gpp(bounds%begp:bounds%endp), &
         col_cf_input%gpp(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%ar(bounds%begp:bounds%endp), &
         col_cf_input%ar(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%npp(bounds%begp:bounds%endp), &
         col_cf_input%npp(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%vegfire(bounds%begp:bounds%endp), &
         col_cf_input%vegfire(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%wood_harvestc(bounds%begp:bounds%endp), &
         col_cf_input%wood_harvestc(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_closs(bounds%begp:bounds%endp), &
         col_cf_input%fire_closs_p2c(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%litfall(bounds%begp:bounds%endp), &
         col_cf_input%litfall(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%hrv_xsmrpool_to_atm(bounds%begp:bounds%endp), &
         col_cf_input%hrv_xsmrpool_to_atm(bounds%begc:bounds%endc))
    
  end subroutine veg_cf_summary 
  
  !------------------------------------------------------------  
  subroutine veg_cf_summary_rr(this, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf_input)
    !
    ! !DESCRIPTION:
    ! summarize root respiration
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this  
    type(bounds_type), intent(in) :: bounds  
    integer, intent(in) :: num_soilp
    integer, intent(in) :: filter_soilp(:)
    integer, intent(in) :: num_soilc
    integer, intent(in) :: filter_soilc(:)
    type(column_carbon_flux), intent(inout) :: col_cf_input
    !
    ! !LOCAL VARIABLES
    integer :: fp, p
    !------------------------------------------------------------  

    do fp = 1,num_soilp
      p = filter_soilp(fp)  
      ! root respiration (RR)
      this%rr(p) = &
      this%froot_mr(p) + &
      this%cpool_froot_gr(p) + &
      this%cpool_livecroot_gr(p) + &
      this%cpool_deadcroot_gr(p) + &
      this%transfer_froot_gr(p) + &
      this%transfer_livecroot_gr(p) + &
      this%transfer_deadcroot_gr(p) + &
      this%cpool_froot_storage_gr(p) + &
      this%cpool_livecroot_storage_gr(p) + &
      this%cpool_deadcroot_storage_gr(p)
    enddo  
      call p2c(bounds, num_soilc, filter_soilc, &
           this%rr(bounds%begp:bounds%endp), &
           col_cf_input%rr(bounds%begc:bounds%endc))
           
  end subroutine veg_cf_summary_rr

  !------------------------------------------------------------  
  subroutine veg_cf_summary_for_ch4( this, bounds, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! summarize vegetation-level fluxes for methane calculation
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this  
    type(bounds_type), intent(in) :: bounds  
    integer, intent(in) :: num_soilp
    integer, intent(in) :: filter_soilp(:)
    !
    ! !LOCAL VARIABLES
    integer :: fp, p
    !------------------------------------------------------------  

    do fp = 1,num_soilp
       p = filter_soilp(fp)   

       ! aboveground NPP: leaf, live stem, dead stem (AGNPP)
       ! This is supposed to correspond as closely as possible to
       ! field measurements of AGNPP, so it ignores the storage pools
       ! and only treats the fluxes into displayed pools.

       this%agnpp(p) = &
            this%cpool_to_leafc(p)                  + &
            this%leafc_xfer_to_leafc(p)             + &
            this%cpool_to_livestemc(p)              + &
            this%livestemc_xfer_to_livestemc(p)     + &
            this%cpool_to_deadstemc(p)              + &
            this%deadstemc_xfer_to_deadstemc(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%agnpp(p) =                    &
               this%agnpp(p)               + &
               this%cpool_to_grainc(p)     + &
               this%grainc_xfer_to_grainc(p)
       endif

       ! belowground NPP: fine root, live coarse root, dead coarse root (BGNPP)
       ! This is supposed to correspond as closely as possible to
       ! field measurements of BGNPP, so it ignores the storage pools
       ! and only treats the fluxes into displayed pools.
       this%bgnpp(p) = &
            this%cpool_to_frootc(p)                   + &
            this%frootc_xfer_to_frootc(p)             + &
            this%cpool_to_livecrootc(p)               + &
            this%livecrootc_xfer_to_livecrootc(p)     + &
            this%cpool_to_deadcrootc(p)               + &
            this%deadcrootc_xfer_to_deadcrootc(p)

       this%agwdnpp(p) = &
            this%cpool_to_livestemc(p)              + &
            this%livestemc_xfer_to_livestemc(p)     + &
            this%cpool_to_deadstemc(p)              + &
            this%deadstemc_xfer_to_deadstemc(p)

    end do
  
  end subroutine veg_cf_summary_for_ch4
  
  !-----------------------------------------------------------------------
  subroutine veg_cf_setvalues ( this, num_patch, filter_patch, value_patch)
    !
    ! !DESCRIPTION:
    ! Set vegetation-level carbon fluxes
    !
    ! !ARGUMENTS:
    class (vegetation_carbon_flux) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    !------------------------------------------------------------------------

    if(.not.use_fates) then
       do fi = 1,num_patch
          i = filter_patch(fi)

          this%m_leafc_to_litter(i)                   = value_patch
          this%m_frootc_to_litter(i)                  = value_patch
          this%m_leafc_storage_to_litter(i)           = value_patch
          this%m_frootc_storage_to_litter(i)          = value_patch
          this%m_livestemc_storage_to_litter(i)       = value_patch
          this%m_deadstemc_storage_to_litter(i)       = value_patch
          this%m_livecrootc_storage_to_litter(i)      = value_patch
          this%m_deadcrootc_storage_to_litter(i)      = value_patch
          this%m_leafc_xfer_to_litter(i)              = value_patch
          this%m_frootc_xfer_to_litter(i)             = value_patch
          this%m_livestemc_xfer_to_litter(i)          = value_patch
          this%m_deadstemc_xfer_to_litter(i)          = value_patch
          this%m_livecrootc_xfer_to_litter(i)         = value_patch
          this%m_deadcrootc_xfer_to_litter(i)         = value_patch
          this%m_livestemc_to_litter(i)               = value_patch
          this%m_deadstemc_to_litter(i)               = value_patch
          this%m_livecrootc_to_litter(i)              = value_patch
          this%m_deadcrootc_to_litter(i)              = value_patch
          this%m_gresp_storage_to_litter(i)           = value_patch
          this%m_gresp_xfer_to_litter(i)              = value_patch
          this%m_cpool_to_litter(i)                   = value_patch
          this%hrv_leafc_to_litter(i)                 = value_patch             
          this%hrv_leafc_storage_to_litter(i)         = value_patch     
          this%hrv_leafc_xfer_to_litter(i)            = value_patch        
          this%hrv_frootc_to_litter(i)                = value_patch            
          this%hrv_frootc_storage_to_litter(i)        = value_patch    
          this%hrv_frootc_xfer_to_litter(i)           = value_patch       
          this%hrv_livestemc_to_litter(i)             = value_patch         
          this%hrv_livestemc_storage_to_litter(i)     = value_patch 
          this%hrv_livestemc_xfer_to_litter(i)        = value_patch    
          this%hrv_deadstemc_to_prod10c(i)            = value_patch        
          this%hrv_deadstemc_to_prod100c(i)           = value_patch       
          this%hrv_leafc_to_prod1c(i)                 = value_patch
          this%hrv_livestemc_to_prod1c(i)             = value_patch
          this%hrv_grainc_to_prod1c(i)                = value_patch
          this%hrv_cropc_to_prod1c(i)                 = value_patch
          this%hrv_deadstemc_storage_to_litter(i)     = value_patch 
          this%hrv_deadstemc_xfer_to_litter(i)        = value_patch    
          this%hrv_livecrootc_to_litter(i)            = value_patch        
          this%hrv_livecrootc_storage_to_litter(i)    = value_patch
          this%hrv_livecrootc_xfer_to_litter(i)       = value_patch   
          this%hrv_deadcrootc_to_litter(i)            = value_patch        
          this%hrv_deadcrootc_storage_to_litter(i)    = value_patch
          this%hrv_deadcrootc_xfer_to_litter(i)       = value_patch   
          this%hrv_gresp_storage_to_litter(i)         = value_patch     
          this%hrv_gresp_xfer_to_litter(i)            = value_patch        
          this%hrv_xsmrpool_to_atm(i)                 = value_patch
          this%hrv_cpool_to_litter(i)                 = value_patch

          this%m_leafc_to_fire(i)                     = value_patch
          this%m_leafc_storage_to_fire(i)             = value_patch
          this%m_leafc_xfer_to_fire(i)                = value_patch
          this%m_livestemc_to_fire(i)                 = value_patch
          this%m_livestemc_storage_to_fire(i)         = value_patch
          this%m_livestemc_xfer_to_fire(i)            = value_patch
          this%m_deadstemc_to_fire(i)                 = value_patch
          this%m_deadstemc_storage_to_fire(i)         = value_patch
          this%m_deadstemc_xfer_to_fire(i)            = value_patch
          this%m_frootc_to_fire(i)                    = value_patch
          this%m_frootc_storage_to_fire(i)            = value_patch
          this%m_frootc_xfer_to_fire(i)               = value_patch
          this%m_livecrootc_to_fire(i)                = value_patch
          this%m_livecrootc_storage_to_fire(i)        = value_patch
          this%m_livecrootc_xfer_to_fire(i)           = value_patch
          this%m_deadcrootc_to_fire(i)                = value_patch
          this%m_deadcrootc_storage_to_fire(i)        = value_patch
          this%m_deadcrootc_xfer_to_fire(i)           = value_patch
          this%m_gresp_storage_to_fire(i)             = value_patch
          this%m_gresp_xfer_to_fire(i)                = value_patch
          this%m_cpool_to_fire(i)                     = value_patch

          this%m_leafc_to_litter_fire(i)              = value_patch
          this%m_leafc_storage_to_litter_fire(i)      = value_patch
          this%m_leafc_xfer_to_litter_fire(i)         = value_patch
          this%m_livestemc_to_litter_fire(i)          = value_patch
          this%m_livestemc_storage_to_litter_fire(i)  = value_patch
          this%m_livestemc_xfer_to_litter_fire(i)     = value_patch
          this%m_livestemc_to_deadstemc_fire(i)       = value_patch
          this%m_deadstemc_to_litter_fire(i)          = value_patch
          this%m_deadstemc_storage_to_litter_fire(i)  = value_patch
          this%m_deadstemc_xfer_to_litter_fire(i)     = value_patch
          this%m_frootc_to_litter_fire(i)             = value_patch
          this%m_frootc_storage_to_litter_fire(i)     = value_patch
          this%m_frootc_xfer_to_litter_fire(i)        = value_patch
          this%m_livecrootc_to_litter_fire(i)         = value_patch
          this%m_livecrootc_storage_to_litter_fire(i) = value_patch
          this%m_livecrootc_xfer_to_litter_fire(i)    = value_patch
          this%m_livecrootc_to_deadcrootc_fire(i)     = value_patch
          this%m_deadcrootc_to_litter_fire(i)         = value_patch
          this%m_deadcrootc_storage_to_litter_fire(i) = value_patch
          this%m_deadcrootc_xfer_to_litter_fire(i)    = value_patch
          this%m_gresp_storage_to_litter_fire(i)      = value_patch
          this%m_gresp_xfer_to_litter_fire(i)         = value_patch
          this%m_cpool_to_litter_fire(i)              = value_patch

          this%leafc_xfer_to_leafc(i)                 = value_patch
          this%frootc_xfer_to_frootc(i)               = value_patch
          this%livestemc_xfer_to_livestemc(i)         = value_patch
          this%deadstemc_xfer_to_deadstemc(i)         = value_patch
          this%livecrootc_xfer_to_livecrootc(i)       = value_patch
          this%deadcrootc_xfer_to_deadcrootc(i)       = value_patch
          this%leafc_to_litter(i)                     = value_patch
          this%frootc_to_litter(i)                    = value_patch
          this%leaf_mr(i)                             = value_patch
          this%froot_mr(i)                            = value_patch
          this%livestem_mr(i)                         = value_patch
          this%livecroot_mr(i)                        = value_patch
          this%grain_mr(i)                            = value_patch
          this%leaf_curmr(i)                          = value_patch
          this%froot_curmr(i)                         = value_patch
          this%livestem_curmr(i)                      = value_patch
          this%livecroot_curmr(i)                     = value_patch
          this%grain_curmr(i)                         = value_patch
          this%leaf_xsmr(i)                           = value_patch
          this%froot_xsmr(i)                          = value_patch
          this%livestem_xsmr(i)                       = value_patch
          this%livecroot_xsmr(i)                      = value_patch
          this%grain_xsmr(i)                          = value_patch
          this%xr(i)                                  = value_patch
          this%psnsun_to_cpool(i)                     = value_patch
          this%psnshade_to_cpool(i)                   = value_patch
          this%cpool_to_xsmrpool(i)                   = value_patch
          this%cpool_to_leafc(i)                      = value_patch
          this%cpool_to_leafc_storage(i)              = value_patch
          this%cpool_to_frootc(i)                     = value_patch
          this%cpool_to_frootc_storage(i)             = value_patch
          this%cpool_to_livestemc(i)                  = value_patch
          this%cpool_to_livestemc_storage(i)          = value_patch
          this%cpool_to_deadstemc(i)                  = value_patch
          this%cpool_to_deadstemc_storage(i)          = value_patch
          this%cpool_to_livecrootc(i)                 = value_patch
          this%cpool_to_livecrootc_storage(i)         = value_patch
          this%cpool_to_deadcrootc(i)                 = value_patch
          this%cpool_to_deadcrootc_storage(i)         = value_patch
          this%cpool_to_gresp_storage(i)              = value_patch
          this%cpool_leaf_gr(i)                       = value_patch
          this%cpool_leaf_storage_gr(i)               = value_patch
          this%transfer_leaf_gr(i)                    = value_patch
          this%cpool_froot_gr(i)                      = value_patch
          this%cpool_froot_storage_gr(i)              = value_patch
          this%transfer_froot_gr(i)                   = value_patch
          this%cpool_livestem_gr(i)                   = value_patch
          this%cpool_livestem_storage_gr(i)           = value_patch
          this%transfer_livestem_gr(i)                = value_patch
          this%cpool_deadstem_gr(i)                   = value_patch
          this%cpool_deadstem_storage_gr(i)           = value_patch
          this%transfer_deadstem_gr(i)                = value_patch
          this%cpool_livecroot_gr(i)                  = value_patch
          this%cpool_livecroot_storage_gr(i)          = value_patch
          this%transfer_livecroot_gr(i)               = value_patch
          this%cpool_deadcroot_gr(i)                  = value_patch
          this%cpool_deadcroot_storage_gr(i)          = value_patch
          this%transfer_deadcroot_gr(i)               = value_patch
          this%leafc_storage_to_xfer(i)               = value_patch
          this%frootc_storage_to_xfer(i)              = value_patch
          this%livestemc_storage_to_xfer(i)           = value_patch
          this%deadstemc_storage_to_xfer(i)           = value_patch
          this%livecrootc_storage_to_xfer(i)          = value_patch
          this%deadcrootc_storage_to_xfer(i)          = value_patch
          this%gresp_storage_to_xfer(i)               = value_patch
          this%livestemc_to_deadstemc(i)              = value_patch
          this%livecrootc_to_deadcrootc(i)            = value_patch
          this%gpp(i)                                 = value_patch
          this%gpp_before_downreg(i)                  = value_patch
          this%mr(i)                                  = value_patch
          this%current_gr(i)                          = value_patch
          this%transfer_gr(i)                         = value_patch
          this%storage_gr(i)                          = value_patch
          this%gr(i)                                  = value_patch
          this%ar(i)                                  = value_patch
          this%rr(i)                                  = value_patch
          this%npp(i)                                 = value_patch 
          this%agnpp(i)                               = value_patch
          this%bgnpp(i)                               = value_patch
          this%agwdnpp(i)                               = value_patch
          this%litfall(i)                             = value_patch
          this%vegfire(i)                             = value_patch
          this%wood_harvestc(i)                       = value_patch
          this%cinputs(i)                             = value_patch
          this%coutputs(i)                            = value_patch
          this%fire_closs(i)                          = value_patch
          this%frootc_alloc(i)                        = value_patch
          this%frootc_loss(i)                         = value_patch
          this%leafc_alloc(i)                         = value_patch
          this%leafc_loss(i)                          = value_patch
          this%woodc_alloc(i)                         = value_patch
          this%woodc_loss(i)                          = value_patch
          this%xsmrpool_turnover(i)                   = value_patch
       end do
    end if !(.not.use_fates)

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%xsmrpool_to_atm(i)         = value_patch
          this%livestemc_to_litter(i)     = value_patch
          this%grainc_to_food(i)          = value_patch
          this%grainc_xfer_to_grainc(i)   = value_patch
          this%cpool_to_grainc(i)         = value_patch
          this%cpool_to_grainc_storage(i) = value_patch
          this%cpool_grain_gr(i)          = value_patch
          this%cpool_grain_storage_gr(i)  = value_patch
          this%transfer_grain_gr(i)       = value_patch
          this%grainc_storage_to_xfer(i)  = value_patch
          this%crop_seedc_to_leaf(i)      = value_patch
       end do
    end if
  
  end subroutine veg_cf_setvalues 
  
  !------------------------------------------------------------------------
  subroutine veg_cf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_cf_clean
  
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

    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_nf
    !-----------------------------------------------------------------------
    ! add suffix if number of soil decomposition depths is greater than 1

    this%m_leafn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality', &
         ptr_patch=this%m_leafn_to_litter, default='inactive')

    this%m_frootn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality', &
         ptr_patch=this%m_frootn_to_litter, default='inactive')

    this%m_leafn_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage mortality', &
         ptr_patch=this%m_leafn_storage_to_litter, default='inactive')

    this%m_frootn_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage mortality', &
         ptr_patch=this%m_frootn_storage_to_litter, default='inactive')

    this%m_livestemn_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage mortality', &
         ptr_patch=this%m_livestemn_storage_to_litter, default='inactive')

    this%m_deadstemn_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage mortality', &
         ptr_patch=this%m_deadstemn_storage_to_litter, default='inactive')

    this%m_livecrootn_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage mortality', &
         ptr_patch=this%m_livecrootn_storage_to_litter, default='inactive')

    this%m_deadcrootn_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage mortality', &
         ptr_patch=this%m_deadcrootn_storage_to_litter, default='inactive')

    this%m_leafn_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer mortality', &
         ptr_patch=this%m_leafn_xfer_to_litter, default='inactive')

    this%m_frootn_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer mortality', &
         ptr_patch=this%m_frootn_xfer_to_litter, default='inactive')

    this%m_livestemn_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer mortality', &
         ptr_patch=this%m_livestemn_xfer_to_litter, default='inactive')

    this%m_deadstemn_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer mortality', &
         ptr_patch=this%m_deadstemn_xfer_to_litter, default='inactive')

    this%m_livecrootn_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer mortality', &
         ptr_patch=this%m_livecrootn_xfer_to_litter, default='inactive')

    this%m_deadcrootn_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer mortality', &
         ptr_patch=this%m_deadcrootn_xfer_to_litter, default='inactive')

    this%m_livestemn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N mortality', &
         ptr_patch=this%m_livestemn_to_litter, default='inactive')

    this%m_deadstemn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N mortality', &
         ptr_patch=this%m_deadstemn_to_litter, default='inactive')

    this%m_livecrootn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N mortality', &
         ptr_patch=this%m_livecrootn_to_litter, default='inactive')

    this%m_deadcrootn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N mortality', &
         ptr_patch=this%m_deadcrootn_to_litter, default='inactive')

    this%m_retransn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool mortality', &
         ptr_patch=this%m_retransn_to_litter, default='inactive')

    this%m_npool_to_litter_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_NPOOL_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='N pool mortality due to fire', &
         ptr_patch=this%m_npool_to_litter_fire, default='inactive')

    this%m_npool_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_NPOOL_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool mortality', &
         ptr_patch=this%m_npool_to_litter, default='inactive')

    this%m_leafn_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N fire loss', &
         ptr_patch=this%m_leafn_to_fire, default='inactive')

    this%m_frootn_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N fire loss ', &
         ptr_patch=this%m_frootn_to_fire, default='inactive')

    this%m_leafn_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage fire loss', &
         ptr_patch=this%m_leafn_storage_to_fire, default='inactive')

    this%m_frootn_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage fire loss', &
         ptr_patch=this%m_frootn_storage_to_fire, default='inactive')

    this%m_livestemn_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage fire loss', &
         ptr_patch=this%m_livestemn_storage_to_fire, default='inactive')

    this%m_deadstemn_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage fire loss', &
         ptr_patch=this%m_deadstemn_storage_to_fire, default='inactive')

    this%m_livecrootn_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage fire loss', &
         ptr_patch=this%m_livecrootn_storage_to_fire, default='inactive')

    this%m_deadcrootn_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage fire loss', &
         ptr_patch=this%m_deadcrootn_storage_to_fire, default='inactive')

    this%m_leafn_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer fire loss', &
         ptr_patch=this%m_leafn_xfer_to_fire, default='inactive')

    this%m_frootn_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer fire loss', &
         ptr_patch=this%m_frootn_xfer_to_fire, default='inactive')

    this%m_livestemn_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer fire loss', &
         ptr_patch=this%m_livestemn_xfer_to_fire, default='inactive')

    this%m_deadstemn_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer fire loss', &
         ptr_patch=this%m_deadstemn_xfer_to_fire, default='inactive')

    this%m_livecrootn_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer fire loss', &
         ptr_patch=this%m_livecrootn_xfer_to_fire, default='inactive')

    this%m_deadcrootn_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer fire loss', &
         ptr_patch=this%m_deadcrootn_xfer_to_fire, default='inactive')

    this%m_livestemn_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N fire loss', &
         ptr_patch=this%m_livestemn_to_fire, default='inactive')

    this%m_deadstemn_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire loss', &
         ptr_patch=this%m_deadstemn_to_fire, default='inactive')

    this%m_deadstemn_to_litter_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire mortality to litter', &
         ptr_patch=this%m_deadstemn_to_litter_fire, default='inactive')

    this%m_livecrootn_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N fire loss', &
         ptr_patch=this%m_livecrootn_to_fire, default='inactive')

    this%m_deadcrootn_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire loss', &
         ptr_patch=this%m_deadcrootn_to_fire, default='inactive')

    this%m_deadcrootn_to_litter_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire mortality to litter', &
         ptr_patch=this%m_deadcrootn_to_litter_fire, default='inactive')

    this%m_retransn_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool fire loss', &
         ptr_patch=this%m_retransn_to_fire, default='inactive')

    this%m_npool_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_NPOOL_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool fire loss', &
         ptr_patch=this%m_npool_to_fire, default='inactive')

    this%leafn_xfer_to_leafn(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_XFER_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N growth from storage', &
         ptr_patch=this%leafn_xfer_to_leafn, default='inactive')

    this%frootn_xfer_to_frootn(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_XFER_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N growth from storage', &
         ptr_patch=this%frootn_xfer_to_frootn, default='inactive')

    this%livestemn_xfer_to_livestemn(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_XFER_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N growth from storage', &
         ptr_patch=this%livestemn_xfer_to_livestemn, default='inactive')

    this%deadstemn_xfer_to_deadstemn(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_XFER_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N growth from storage', &
         ptr_patch=this%deadstemn_xfer_to_deadstemn, default='inactive')

    this%livecrootn_xfer_to_livecrootn(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_XFER_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N growth from storage', &
         ptr_patch=this%livecrootn_xfer_to_livecrootn, default='inactive')

    this%deadcrootn_xfer_to_deadcrootn(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_XFER_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N growth from storage', &
         ptr_patch=this%deadcrootn_xfer_to_deadcrootn, default='inactive')

    this%leafn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall', &
         ptr_patch=this%leafn_to_litter, default='inactive')

    this%leafn_to_retransn(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N to retranslocated N pool', &
         ptr_patch=this%leafn_to_retransn, default='inactive')

    this%frootn_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall', &
         ptr_patch=this%frootn_to_litter, default='inactive')

    this%retransn_to_npool(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated N', &
         ptr_patch=this%retransn_to_npool)

    this%sminn_to_npool(begp:endp) = spval
    call hist_addfld1d (fname='SMINN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of soil mineral N uptake', &
         ptr_patch=this%sminn_to_npool)

    this%npool_to_leafn(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N', &
         ptr_patch=this%npool_to_leafn, default='inactive')

    this%npool_to_leafn_storage(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LEAFN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N storage', &
         ptr_patch=this%npool_to_leafn_storage, default='inactive')

    this%npool_to_frootn(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N', &
         ptr_patch=this%npool_to_frootn, default='inactive')

    this%npool_to_frootn_storage(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_FROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N storage', &
         ptr_patch=this%npool_to_frootn_storage, default='inactive')

    this%npool_to_livestemn(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N', &
         ptr_patch=this%npool_to_livestemn, default='inactive')

    this%npool_to_livestemn_storage(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N storage', &
         ptr_patch=this%npool_to_livestemn_storage, default='inactive')

    this%npool_to_deadstemn(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N', &
         ptr_patch=this%npool_to_deadstemn, default='inactive')

    this%npool_to_deadstemn_storage(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N storage', &
         ptr_patch=this%npool_to_deadstemn_storage, default='inactive')

    this%npool_to_livecrootn(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N', &
         ptr_patch=this%npool_to_livecrootn, default='inactive')

    this%npool_to_livecrootn_storage(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N storage', &
         ptr_patch=this%npool_to_livecrootn_storage, default='inactive')

    this%npool_to_deadcrootn(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N', &
         ptr_patch=this%npool_to_deadcrootn, default='inactive')

    this%npool_to_deadcrootn_storage(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N storage', &
         ptr_patch=this%npool_to_deadcrootn_storage, default='inactive')

    this%leafn_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N shift storage to transfer', &
         ptr_patch=this%leafn_storage_to_xfer, default='inactive')

    this%frootn_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N shift storage to transfer', &
         ptr_patch=this%frootn_storage_to_xfer, default='inactive')

    this%livestemn_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N shift storage to transfer', &
         ptr_patch=this%livestemn_storage_to_xfer, default='inactive')

    this%deadstemn_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N shift storage to transfer', &
         ptr_patch=this%deadstemn_storage_to_xfer, default='inactive')

    this%livecrootn_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N shift storage to transfer', &
         ptr_patch=this%livecrootn_storage_to_xfer, default='inactive')

    this%deadcrootn_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N shift storage to transfer', &
         ptr_patch=this%deadcrootn_storage_to_xfer, default='inactive')

    this%livestemn_to_deadstemn(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N turnover', &
         ptr_patch=this%livestemn_to_deadstemn, default='inactive')

    this%livestemn_to_retransn(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N to retranslocated N pool', &
         ptr_patch=this%livestemn_to_retransn, default='inactive')

    this%livecrootn_to_deadcrootn(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N turnover', &
         ptr_patch=this%livecrootn_to_deadcrootn, default='inactive')

    this%livecrootn_to_retransn(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N to retranslocated N pool', &
         ptr_patch=this%livecrootn_to_retransn, default='inactive')

    this%ndeploy(begp:endp) = spval
    call hist_addfld1d (fname='NDEPLOY', units='gN/m^2/s', &
         avgflag='A', long_name='total N deployed in new growth', &
         ptr_patch=this%ndeploy)

    this%wood_harvestn(begp:endp) = spval
    call hist_addfld1d (fname='WOOD_HARVESTN', units='gN/m^2/s', &
         avgflag='A', long_name='wood harvest N (to product pools)', &
         ptr_patch=this%wood_harvestn)

    this%fire_nloss(begp:endp) = spval
    call hist_addfld1d (fname='PFT_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total pft-level fire N loss', &
         ptr_patch=this%fire_nloss)

    if (crop_prog) then
       this%fert(begp:endp) = spval
       call hist_addfld1d (fname='FERT', units='gN/m^2/s', &
            avgflag='A', long_name='fertilizer added', &
            ptr_patch=this%fert)
    end if

    if (crop_prog) then
       this%soyfixn(begp:endp) = spval
       call hist_addfld1d (fname='SOYFIXN', units='gN/m^2/s', &
            avgflag='A', long_name='soybean fixation', &
            ptr_patch=this%soyfixn)
    end if

    if (crop_prog) then
       this%fert_counter(begp:endp) = spval
       call hist_addfld1d (fname='FERT_COUNTER', units='seconds', &
            avgflag='A', long_name='time left to fertilize', &
            ptr_patch=this%fert_counter)
    end if

    this%crop_seedn_to_leaf(begp:endp) = spval
    call hist_addfld1d (fname='CROP_SEEDN_TO_LEAF', units='gN/m^2/s', &
         avgflag='A', long_name='crop seed source to leaf', &
         ptr_patch=this%crop_seedn_to_leaf, default='inactive')

    this%dwt_seedn_to_leaf(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF_PATCH', units='gN/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level leaf ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_seedn_to_leaf, default='inactive')

    this%dwt_seedn_to_deadstem(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM_PATCH', units='gN/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level deadstem ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_seedn_to_deadstem, default='inactive')

    this%dwt_conv_nflux(begp:endp) = spval
    call hist_addfld1d (fname='DWT_CONV_NFLUX_PATCH', units='gN/m^2/s', &
         avgflag='A', &
         long_name='patch-level conversion C flux (immediate loss to atm) ' // &
         '(0 at all times except first timestep of year) ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_conv_nflux, default='inactive')

    this%dwt_seedn_to_npool(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_NPOOL_PATCH', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level npool', &
         ptr_patch=this%dwt_seedn_to_npool, default='inactive')

    this%plant_ndemand(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_NDEMAND', units='gN/m^2/s', &
         avgflag='A', long_name='N flux required to support initial GPP', &
         ptr_patch=this%plant_ndemand)

    this%avail_retransn(begp:endp) = spval
    call hist_addfld1d (fname='AVAIL_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='N flux available from retranslocation pool', &
         ptr_patch=this%avail_retransn, default='inactive')

    this%plant_nalloc(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_NALLOC', units='gN/m^2/s', &
         avgflag='A', long_name='total allocated N flux', &
         ptr_patch=this%plant_nalloc, default='inactive')

    this%gap_nloss_litter(begp:endp) = spval
    call hist_addfld1d (fname='GAP_NLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total nloss from veg to litter due to gap mortality', &
         ptr_patch=this%gap_nloss_litter, default='inactive')

    this%fire_nloss_litter(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_NLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total nloss from veg to litter due to fire mortality', &
         ptr_patch=this%fire_nloss_litter, default='inactive')

    this%hrv_nloss_litter(begp:endp) = spval
    call hist_addfld1d (fname='HRV_NLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total nloss from veg to litter due to harvest mortality', &
         ptr_patch=this%hrv_nloss_litter, default='inactive')

    this%sen_nloss_litter(begp:endp) = spval
    call hist_addfld1d (fname='SEN_NLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total nloss from veg to litter pool due to senescence', &
         ptr_patch=this%sen_nloss_litter, default='inactive')

    ! Note: the following are vegetation-level variables, reported to the column-level
    this%dwt_prod10n_gain(begp:endp) = spval
    call hist_addfld1d (fname='DWT_PROD10N_GAIN_PATCH', units='gN/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_col=this%dwt_prod10n_gain, default='inactive')

    this%dwt_prod100n_gain(begp:endp) = spval
    call hist_addfld1d (fname='DWT_PROD100N_GAIN_PATCH', units='gN/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_col=this%dwt_prod100n_gain, default='inactive')

         
    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_nf
    !-----------------------------------------------------------------------
    num_special_patch = 0
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    do p = begp,endp
       l = veg_pp%landunit(p)
       
       this%prev_leafn_to_litter(p)  = 0._r8 
       this%prev_frootn_to_litter(p) = 0._r8 
       
       if ( crop_prog )then
          this%fert_counter(p)  = spval
          this%fert(p)          = 0._r8 
          this%soyfixn(p)       = 0._r8 
       end if

       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%fert_counter(p)  = 0._r8
       end if

       if (lun_pp%ifspecial(l)) then
          this%plant_ndemand(p)  = spval
          this%avail_retransn(p) = spval
          this%plant_nalloc(p)   = spval
       end if
    end do

    call this%SetValues (num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8)
    
  end subroutine veg_nf_init
    
  !-----------------------------------------------------------------------
  subroutine veg_nf_restart (this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write restart data for vegetation-level nitrogen fluxes
    !
    ! !ARGUMENTS:
    class (vegetation_nitrogen_flux)  :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !------------------------------------------------------------------------
    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag, varname='fert_counter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_counter)

       call restartvar(ncid=ncid, flag=flag, varname='fert', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer_to_grainn', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain N growth from storage', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_xfer_to_grainn)

       call restartvar(ncid=ncid, flag=flag,  varname='livestemn_to_litter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='livestem N to litter', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemn_to_litter)
    
       call restartvar(ncid=ncid, flag=flag,  varname='grainn_to_food', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain N to food', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_to_food)
    
       call restartvar(ncid=ncid, flag=flag,  varname='npool_to_grainn', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain N', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%npool_to_grainn)
    
       call restartvar(ncid=ncid, flag=flag,  varname='npool_to_grainn_storage', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain N storage', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%npool_to_grainn_storage)
    
       call restartvar(ncid=ncid, flag=flag, varname='grainn_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain N shift storage to transfer', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_storage_to_xfer)
    end if ! crop_prog

    call restartvar(ncid=ncid, flag=flag, varname='plant_ndemand', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_ndemand) 

    call restartvar(ncid=ncid, flag=flag, varname='avail_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%avail_retransn) 

    call restartvar(ncid=ncid, flag=flag, varname='plant_nalloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_nalloc) 

  end subroutine veg_nf_restart  
  
  !-----------------------------------------------------------------------
  subroutine veg_nf_setvalues ( this, num_patch, filter_patch, value_patch)
    !
    ! !DESCRIPTION:
    ! Set vegetation-level nitrogen fluxes
    !
    ! !ARGUMENTS:
    class (vegetation_nitrogen_flux) :: this
    integer , intent(in)             :: num_patch
    integer , intent(in)             :: filter_patch(:)
    real(r8), intent(in)             :: value_patch
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    !------------------------------------------------------------------------
  
    do fi = 1,num_patch
       i=filter_patch(fi)

       this%m_leafn_to_litter(i)                   = value_patch
       this%m_frootn_to_litter(i)                  = value_patch
       this%m_leafn_storage_to_litter(i)           = value_patch
       this%m_frootn_storage_to_litter(i)          = value_patch
       this%m_livestemn_storage_to_litter(i)       = value_patch
       this%m_deadstemn_storage_to_litter(i)       = value_patch
       this%m_livecrootn_storage_to_litter(i)      = value_patch
       this%m_deadcrootn_storage_to_litter(i)      = value_patch
       this%m_leafn_xfer_to_litter(i)              = value_patch
       this%m_frootn_xfer_to_litter(i)             = value_patch
       this%m_livestemn_xfer_to_litter(i)          = value_patch
       this%m_deadstemn_xfer_to_litter(i)          = value_patch
       this%m_livecrootn_xfer_to_litter(i)         = value_patch
       this%m_deadcrootn_xfer_to_litter(i)         = value_patch
       this%m_livestemn_to_litter(i)               = value_patch
       this%m_deadstemn_to_litter(i)               = value_patch
       this%m_livecrootn_to_litter(i)              = value_patch
       this%m_deadcrootn_to_litter(i)              = value_patch
       this%m_retransn_to_litter(i)                = value_patch
       this%m_npool_to_litter(i)                   = value_patch
       this%hrv_leafn_to_litter(i)                 = value_patch             
       this%hrv_frootn_to_litter(i)                = value_patch            
       this%hrv_leafn_storage_to_litter(i)         = value_patch     
       this%hrv_frootn_storage_to_litter(i)        = value_patch    
       this%hrv_livestemn_storage_to_litter(i)     = value_patch 
       this%hrv_deadstemn_storage_to_litter(i)     = value_patch 
       this%hrv_livecrootn_storage_to_litter(i)    = value_patch
       this%hrv_deadcrootn_storage_to_litter(i)    = value_patch
       this%hrv_leafn_xfer_to_litter(i)            = value_patch        
       this%hrv_frootn_xfer_to_litter(i)           = value_patch       
       this%hrv_livestemn_xfer_to_litter(i)        = value_patch    
       this%hrv_deadstemn_xfer_to_litter(i)        = value_patch    
       this%hrv_livecrootn_xfer_to_litter(i)       = value_patch   
       this%hrv_deadcrootn_xfer_to_litter(i)       = value_patch   
       this%hrv_livestemn_to_litter(i)             = value_patch         
       this%hrv_deadstemn_to_prod10n(i)            = value_patch        
       this%hrv_deadstemn_to_prod100n(i)           = value_patch       
       this%hrv_leafn_to_prod1n(i)                 = value_patch
       this%hrv_livestemn_to_prod1n(i)             = value_patch
       this%hrv_grainn_to_prod1n(i)                = value_patch
       this%hrv_cropn_to_prod1n(i)                 = value_patch
       this%hrv_livecrootn_to_litter(i)            = value_patch        
       this%hrv_deadcrootn_to_litter(i)            = value_patch        
       this%hrv_retransn_to_litter(i)              = value_patch    
       this%hrv_npool_to_litter(i)                 = value_patch

       this%m_leafn_to_fire(i)                     = value_patch
       this%m_leafn_storage_to_fire(i)             = value_patch
       this%m_leafn_xfer_to_fire(i)                = value_patch
       this%m_livestemn_to_fire(i)                 = value_patch
       this%m_livestemn_storage_to_fire(i)         = value_patch
       this%m_livestemn_xfer_to_fire(i)            = value_patch
       this%m_deadstemn_to_fire(i)                 = value_patch
       this%m_deadstemn_storage_to_fire(i)         = value_patch
       this%m_deadstemn_xfer_to_fire(i)            = value_patch
       this%m_frootn_to_fire(i)                    = value_patch
       this%m_frootn_storage_to_fire(i)            = value_patch
       this%m_frootn_xfer_to_fire(i)               = value_patch
       this%m_livecrootn_to_fire(i)                = value_patch
       this%m_livecrootn_storage_to_fire(i)        = value_patch
       this%m_livecrootn_xfer_to_fire(i)           = value_patch
       this%m_deadcrootn_to_fire(i)                = value_patch
       this%m_deadcrootn_storage_to_fire(i)        = value_patch
       this%m_deadcrootn_xfer_to_fire(i)           = value_patch
       this%m_retransn_to_fire(i)                  = value_patch
       this%m_npool_to_fire(i)                     = value_patch

       this%m_leafn_to_litter_fire(i)              = value_patch
       this%m_leafn_storage_to_litter_fire(i)      = value_patch
       this%m_leafn_xfer_to_litter_fire(i)         = value_patch
       this%m_livestemn_to_litter_fire(i)          = value_patch
       this%m_livestemn_storage_to_litter_fire(i)  = value_patch
       this%m_livestemn_xfer_to_litter_fire(i)     = value_patch
       this%m_livestemn_to_deadstemn_fire(i)       = value_patch
       this%m_deadstemn_to_litter_fire(i)          = value_patch
       this%m_deadstemn_storage_to_litter_fire(i)  = value_patch
       this%m_deadstemn_xfer_to_litter_fire(i)     = value_patch
       this%m_frootn_to_litter_fire(i)             = value_patch
       this%m_frootn_storage_to_litter_fire(i)     = value_patch
       this%m_frootn_xfer_to_litter_fire(i)        = value_patch
       this%m_livecrootn_to_litter_fire(i)         = value_patch
       this%m_livecrootn_storage_to_litter_fire(i) = value_patch
       this%m_livecrootn_xfer_to_litter_fire(i)    = value_patch
       this%m_livecrootn_to_deadcrootn_fire(i)     = value_patch
       this%m_deadcrootn_to_litter_fire(i)         = value_patch
       this%m_deadcrootn_storage_to_litter_fire(i) = value_patch
       this%m_deadcrootn_xfer_to_litter_fire(i)    = value_patch
       this%m_retransn_to_litter_fire(i)           = value_patch
       this%m_npool_to_litter_fire(i)              = value_patch

       this%leafn_xfer_to_leafn(i)                 = value_patch
       this%frootn_xfer_to_frootn(i)               = value_patch
       this%livestemn_xfer_to_livestemn(i)         = value_patch
       this%deadstemn_xfer_to_deadstemn(i)         = value_patch
       this%livecrootn_xfer_to_livecrootn(i)       = value_patch
       this%deadcrootn_xfer_to_deadcrootn(i)       = value_patch
       this%leafn_to_litter(i)                     = value_patch
       this%leafn_to_retransn(i)                   = value_patch
       this%frootn_to_litter(i)                    = value_patch
       this%retransn_to_npool(i)                   = value_patch
       this%sminn_to_npool(i)                      = value_patch
       this%npool_to_leafn(i)                      = value_patch
       this%npool_to_leafn_storage(i)              = value_patch
       this%npool_to_frootn(i)                     = value_patch
       this%npool_to_frootn_storage(i)             = value_patch
       this%npool_to_livestemn(i)                  = value_patch
       this%npool_to_livestemn_storage(i)          = value_patch
       this%npool_to_deadstemn(i)                  = value_patch
       this%npool_to_deadstemn_storage(i)          = value_patch
       this%npool_to_livecrootn(i)                 = value_patch
       this%npool_to_livecrootn_storage(i)         = value_patch
       this%npool_to_deadcrootn(i)                 = value_patch
       this%npool_to_deadcrootn_storage(i)         = value_patch
       this%leafn_storage_to_xfer(i)               = value_patch
       this%frootn_storage_to_xfer(i)              = value_patch
       this%livestemn_storage_to_xfer(i)           = value_patch
       this%deadstemn_storage_to_xfer(i)           = value_patch
       this%livecrootn_storage_to_xfer(i)          = value_patch
       this%deadcrootn_storage_to_xfer(i)          = value_patch
       this%livestemn_to_deadstemn(i)              = value_patch
       this%livestemn_to_retransn(i)               = value_patch
       this%livecrootn_to_deadcrootn(i)            = value_patch
       this%livecrootn_to_retransn(i)              = value_patch
       this%ndeploy(i)                             = value_patch
       this%wood_harvestn(i)                       = value_patch
       this%fire_nloss(i)                          = value_patch
       this%nfix_to_plantn(i)                      = value_patch
       this%gap_nloss_litter(i)                    = value_patch
       this%fire_nloss_litter(i)                   = value_patch
       this%hrv_nloss_litter(i)                    = value_patch
       this%sen_nloss_litter(i)                    = value_patch
       this%crop_seedn_to_leaf(i)                  = value_patch
       this%livestemn_to_litter(i)                 = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grainn_to_food(i)                   = value_patch
          this%grainn_xfer_to_grainn(i)            = value_patch
          this%npool_to_grainn(i)                  = value_patch
          this%npool_to_grainn_storage(i)          = value_patch
          this%grainn_storage_to_xfer(i)           = value_patch
          this%soyfixn(i)                          = value_patch
          this%frootn_to_retransn(i)               = value_patch
       end do
    end if
  
  end subroutine veg_nf_setvalues  
  
  !-----------------------------------------------------------------------
  subroutine veg_nf_summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_nf)
    !
    ! !DESCRIPTION:
    ! Vegetation-level nitrogen flux summary calculations
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux)             :: this
    type(bounds_type)           , intent(in)    :: bounds          
    integer                     , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                     , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type (column_nitrogen_flux), intent(inout)  :: col_nf          ! column-level nitrogen state for p2c
    
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p             ! indices
    integer  :: fp              ! filter indices
    !-----------------------------------------------------------------------
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! total N deployment (from sminn and retranslocated N pool) (NDEPLOY)
       this%ndeploy(p) = &
            this%sminn_to_npool(p) + &
            this%retransn_to_npool(p)

       ! pft-level wood harvest
       this%wood_harvestn(p) = &
            this%hrv_deadstemn_to_prod10n(p) + &
            this%hrv_deadstemn_to_prod100n(p)
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            this%wood_harvestn(p) = &
            this%wood_harvestn(p) + &
            this%hrv_cropn_to_prod1n(p)
       end if

       ! total pft-level fire N losses
       this%fire_nloss(p) = &
            this%m_leafn_to_fire(p)               + &
            this%m_leafn_storage_to_fire(p)       + &
            this%m_leafn_xfer_to_fire(p)          + &
            this%m_frootn_to_fire(p)              + &
            this%m_frootn_storage_to_fire(p)      + &
            this%m_frootn_xfer_to_fire(p)         + &
            this%m_livestemn_to_fire(p)           + &
            this%m_livestemn_storage_to_fire(p)   + &
            this%m_livestemn_xfer_to_fire(p)      + &
            this%m_deadstemn_to_fire(p)           + &
            this%m_deadstemn_storage_to_fire(p)   + &
            this%m_deadstemn_xfer_to_fire(p)      + &
            this%m_livecrootn_to_fire(p)          + &
            this%m_livecrootn_storage_to_fire(p)  + &
            this%m_livecrootn_xfer_to_fire(p)     + &
            this%m_deadcrootn_to_fire(p)          + &
            this%m_deadcrootn_storage_to_fire(p)  + &
            this%m_deadcrootn_xfer_to_fire(p)     + &
            this%m_retransn_to_fire(p)            + &
            this%m_npool_to_fire(p)

      this%gap_nloss_litter(p) = &
           this%m_leafn_to_litter(p)              + &
           this%m_leafn_storage_to_litter(p)      + &
           this%m_leafn_xfer_to_litter(p)         + &
           this%m_frootn_to_litter(p)             + &
           this%m_frootn_storage_to_litter(p)     + &
           this%m_frootn_xfer_to_litter(p)        + &
           this%m_livestemn_to_litter(p)          + &
           this%m_livestemn_storage_to_litter(p)  + &
           this%m_livestemn_xfer_to_litter(p)     + &
           this%m_deadstemn_to_litter(p)          + &
           this%m_deadstemn_storage_to_litter(p)  + &
           this%m_deadstemn_xfer_to_litter(p)     + &
           this%m_livecrootn_to_litter(p)         + &
           this%m_livecrootn_storage_to_litter(p) + &
           this%m_livecrootn_xfer_to_litter(p)    + &
           this%m_deadcrootn_to_litter(p)         + &
           this%m_deadcrootn_storage_to_litter(p) + &
           this%m_deadcrootn_xfer_to_litter(p)    + &
           this%m_retransn_to_litter(p)           + &
           this%m_npool_to_litter(p)

      this%fire_nloss_litter(p) = &
           this%m_deadstemn_to_litter_fire(p)     + &
           this%m_deadcrootn_to_litter_fire(p)    + &
           this%m_retransn_to_litter_fire(p)      + &
           this%m_npool_to_litter_fire(p)         + &
           this%m_leafn_to_litter_fire(p)         + &
           this%m_frootn_to_litter_fire(p)        + &
           this%m_livestemn_to_litter_fire(p)     + &
           this%m_livecrootn_to_litter_fire(p)    + &
           this%m_leafn_storage_to_litter_fire(p) + &
           this%m_frootn_storage_to_litter_fire(p)       + &
           this%m_livestemn_storage_to_litter_fire(p)    + &
           this%m_deadstemn_storage_to_litter_fire(p)    + &
           this%m_livecrootn_storage_to_litter_fire(p)   + &
           this%m_deadcrootn_storage_to_litter_fire(p)   + &
           this%m_leafn_xfer_to_litter_fire(p)           + &
           this%m_frootn_xfer_to_litter_fire(p)          + &
           this%m_livestemn_xfer_to_litter_fire(p)       + &
           this%m_deadstemn_xfer_to_litter_fire(p)       + &
           this%m_livecrootn_xfer_to_litter_fire(p)      + &
           this%m_deadcrootn_xfer_to_litter_fire(p)

      this%hrv_nloss_litter(p) = &
           this%hrv_retransn_to_litter(p)          + &
           this%hrv_npool_to_litter(p)             + &
           this%hrv_leafn_to_litter(p)             + &
           this%hrv_leafn_storage_to_litter(p)     + &
           this%hrv_leafn_xfer_to_litter(p)        + &
           this%hrv_frootn_to_litter(p)            + &
           this%hrv_frootn_storage_to_litter(p)    + &
           this%hrv_frootn_xfer_to_litter(p)       + &
           this%hrv_livestemn_to_litter(p)         + &
           this%hrv_livestemn_storage_to_litter(p) + &
           this%hrv_livestemn_xfer_to_litter(p)    + &
           this%hrv_deadstemn_storage_to_litter(p) + &
           this%hrv_deadstemn_xfer_to_litter(p)    + &
           this%hrv_livecrootn_to_litter(p)        + &
           this%hrv_livecrootn_storage_to_litter(p)+ &
           this%hrv_livecrootn_xfer_to_litter(p)   + &
           this%hrv_deadcrootn_to_litter(p)        + &
           this%hrv_deadcrootn_storage_to_litter(p)+ &
           this%hrv_deadcrootn_xfer_to_litter(p)
      if (crop_prog) then
         this%sen_nloss_litter(p) = &
             this%livestemn_to_litter(p)            + &
             this%leafn_to_litter(p)                + &
             this%frootn_to_litter(p)
      else
         this%sen_nloss_litter(p) = &
             this%leafn_to_litter(p)                + &
             this%frootn_to_litter(p)
      end if

    end do

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_nloss(bounds%begp:bounds%endp), &
         col_nf%fire_nloss_p2c(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%wood_harvestn(bounds%begp:bounds%endp), &
         col_nf%wood_harvestn(bounds%begc:bounds%endc))

  end subroutine veg_nf_summary  
  
  !------------------------------------------------------------------------
  subroutine veg_nf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_nf_clean
  
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
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_pf
    !-----------------------------------------------------------------------
    this%m_leafp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P mortality', &
         ptr_patch=this%m_leafp_to_litter, default='inactive')

    this%m_frootp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P mortality', &
         ptr_patch=this%m_frootp_to_litter, default='inactive')

    this%m_leafp_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P storage mortality', &
         ptr_patch=this%m_leafp_storage_to_litter, default='inactive')

    this%m_frootp_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P storage mortality', &
         ptr_patch=this%m_frootp_storage_to_litter, default='inactive')

    this%m_livestemp_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P storage mortality', &
         ptr_patch=this%m_livestemp_storage_to_litter, default='inactive')

    this%m_deadstemp_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P storage mortality', &
         ptr_patch=this%m_deadstemp_storage_to_litter, default='inactive')

    this%m_livecrootp_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P storage mortality', &
         ptr_patch=this%m_livecrootp_storage_to_litter, default='inactive')

    this%m_deadcrootp_storage_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P storage mortality', &
         ptr_patch=this%m_deadcrootp_storage_to_litter, default='inactive')

    this%m_leafp_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P transfer mortality', &
         ptr_patch=this%m_leafp_xfer_to_litter, default='inactive')

    this%m_frootp_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P transfer mortality', &
         ptr_patch=this%m_frootp_xfer_to_litter, default='inactive')

    this%m_livestemp_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P transfer mortality', &
         ptr_patch=this%m_livestemp_xfer_to_litter, default='inactive')

    this%m_deadstemp_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P transfer mortality', &
         ptr_patch=this%m_deadstemp_xfer_to_litter, default='inactive')

    this%m_livecrootp_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P transfer mortality', &
         ptr_patch=this%m_livecrootp_xfer_to_litter, default='inactive')

    this%m_deadcrootp_xfer_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P transfer mortality', &
         ptr_patch=this%m_deadcrootp_xfer_to_litter, default='inactive')

    this%m_livestemp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P mortality', &
         ptr_patch=this%m_livestemp_to_litter, default='inactive')

    this%m_deadstemp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P mortality', &
         ptr_patch=this%m_deadstemp_to_litter, default='inactive')

    this%m_livecrootp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P mortality', &
         ptr_patch=this%m_livecrootp_to_litter, default='inactive')

    this%m_deadcrootp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P mortality', &
         ptr_patch=this%m_deadcrootp_to_litter, default='inactive')

    this%m_retransp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='retranslocated P pool mortality', &
         ptr_patch=this%m_retransp_to_litter, default='inactive')

    this%m_ppool_to_litter_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_PPOOL_TO_LITTER_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='P pool mortality due to fire', &
         ptr_patch=this%m_ppool_to_litter_fire, default='inactive')

    this%m_ppool_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='M_PPOOL_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='Storage P pool mortality', &
         ptr_patch=this%m_ppool_to_litter, default='inactive')

    this%m_leafp_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P fire loss', &
         ptr_patch=this%m_leafp_to_fire, default='inactive')

    this%m_frootp_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P fire loss ', &
         ptr_patch=this%m_frootp_to_fire, default='inactive')

    this%m_leafp_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P storage fire loss', &
         ptr_patch=this%m_leafp_storage_to_fire, default='inactive')

    this%m_frootp_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P storage fire loss', &
         ptr_patch=this%m_frootp_storage_to_fire, default='inactive')

    this%m_livestemp_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P storage fire loss', &
         ptr_patch=this%m_livestemp_storage_to_fire, default='inactive')

    this%m_deadstemp_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P storage fire loss', &
         ptr_patch=this%m_deadstemp_storage_to_fire, default='inactive')

    this%m_livecrootp_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P storage fire loss', &
         ptr_patch=this%m_livecrootp_storage_to_fire, default='inactive')

    this%m_deadcrootp_storage_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P storage fire loss', &
         ptr_patch=this%m_deadcrootp_storage_to_fire, default='inactive')

    this%m_leafp_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P transfer fire loss', &
         ptr_patch=this%m_leafp_xfer_to_fire, default='inactive')

    this%m_frootp_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P transfer fire loss', &
         ptr_patch=this%m_frootp_xfer_to_fire, default='inactive')

    this%m_livestemp_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P transfer fire loss', &
         ptr_patch=this%m_livestemp_xfer_to_fire, default='inactive')

    this%m_deadstemp_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P transfer fire loss', &
         ptr_patch=this%m_deadstemp_xfer_to_fire, default='inactive')

    this%m_livecrootp_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P transfer fire loss', &
         ptr_patch=this%m_livecrootp_xfer_to_fire, default='inactive')

    this%m_deadcrootp_xfer_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P transfer fire loss', &
         ptr_patch=this%m_deadcrootp_xfer_to_fire, default='inactive')

    this%m_livestemp_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P fire loss', &
         ptr_patch=this%m_livestemp_to_fire, default='inactive')

    this%m_deadstemp_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P fire loss', &
         ptr_patch=this%m_deadstemp_to_fire, default='inactive')

    this%m_deadstemp_to_litter_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_TO_LITTER_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P fire mortality to litter', &
         ptr_patch=this%m_deadstemp_to_litter_fire, default='inactive')

    this%m_livecrootp_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P fire loss', &
         ptr_patch=this%m_livecrootp_to_fire, default='inactive')

    this%m_deadcrootp_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P fire loss', &
         ptr_patch=this%m_deadcrootp_to_fire, default='inactive')

    this%m_deadcrootp_to_litter_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_TO_LITTER_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P fire mortality to litter', &
         ptr_patch=this%m_deadcrootp_to_litter_fire, default='inactive')

    this%m_retransp_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='retranslocated P pool fire loss', &
         ptr_patch=this%m_retransp_to_fire, default='inactive')

    this%m_ppool_to_fire(begp:endp) = spval
    call hist_addfld1d (fname='M_PPOOL_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='Storage P pool fire loss', &
         ptr_patch=this%m_ppool_to_fire, default='inactive')

    this%leafp_xfer_to_leafp(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_XFER_TO_LEAFP', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P growth from storage', &
         ptr_patch=this%leafp_xfer_to_leafp, default='inactive')

    this%frootp_xfer_to_frootp(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_XFER_TO_FROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P growth from storage', &
         ptr_patch=this%frootp_xfer_to_frootp, default='inactive')

    this%livestemp_xfer_to_livestemp(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_XFER_TO_LIVESTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P growth from storage', &
         ptr_patch=this%livestemp_xfer_to_livestemp, default='inactive')

    this%deadstemp_xfer_to_deadstemp(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_XFER_TO_DEADSTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P growth from storage', &
         ptr_patch=this%deadstemp_xfer_to_deadstemp, default='inactive')

    this%livecrootp_xfer_to_livecrootp(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_XFER_TO_LIVECROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P growth from storage', &
         ptr_patch=this%livecrootp_xfer_to_livecrootp, default='inactive')

    this%deadcrootp_xfer_to_deadcrootp(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_XFER_TO_DEADCROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P growth from storage', &
         ptr_patch=this%deadcrootp_xfer_to_deadcrootp, default='inactive')

    this%leafp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P litterfall', &
         ptr_patch=this%leafp_to_litter, default='inactive')

    this%leafp_to_retransp(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_TO_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P to retranslocated P pool', &
         ptr_patch=this%leafp_to_retransp, default='inactive')

    this%frootp_to_litter(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P litterfall', &
         ptr_patch=this%frootp_to_litter, default='inactive')

    this%retransp_to_ppool(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSP_TO_PPOOL', units='gP/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated P', &
         ptr_patch=this%retransp_to_ppool)

    this%sminp_to_ppool(begp:endp) = spval
    call hist_addfld1d (fname='SMINP_TO_PPOOL', units='gP/m^2/s', &
         avgflag='A', long_name='deployment of soil mineral P uptake', &
         ptr_patch=this%sminp_to_ppool)

    this%ppool_to_leafp(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LEAFP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to leaf P', &
         ptr_patch=this%ppool_to_leafp, default='inactive')

    this%ppool_to_leafp_storage(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LEAFP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to leaf P storage', &
         ptr_patch=this%ppool_to_leafp_storage, default='inactive')

    this%ppool_to_frootp(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_FROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to fine root P', &
         ptr_patch=this%ppool_to_frootp, default='inactive')

    this%ppool_to_frootp_storage(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_FROOTP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to fine root P storage', &
         ptr_patch=this%ppool_to_frootp_storage, default='inactive')

    this%ppool_to_livestemp(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVESTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live stem P', &
         ptr_patch=this%ppool_to_livestemp, default='inactive')

    this%ppool_to_livestemp_storage(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVESTEMP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live stem P storage', &
         ptr_patch=this%ppool_to_livestemp_storage, default='inactive')

    this%ppool_to_deadstemp(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADSTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead stem P', &
         ptr_patch=this%ppool_to_deadstemp, default='inactive')

    this%ppool_to_deadstemp_storage(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADSTEMP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead stem P storage', &
         ptr_patch=this%ppool_to_deadstemp_storage, default='inactive')

    this%ppool_to_livecrootp(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVECROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root P', &
         ptr_patch=this%ppool_to_livecrootp, default='inactive')

    this%ppool_to_livecrootp_storage(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVECROOTP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root P storage', &
         ptr_patch=this%ppool_to_livecrootp_storage, default='inactive')

    this%ppool_to_deadcrootp(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADCROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root P', &
         ptr_patch=this%ppool_to_deadcrootp, default='inactive')

    this%ppool_to_deadcrootp_storage(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADCROOTP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root P storage', &
         ptr_patch=this%ppool_to_deadcrootp_storage, default='inactive')

    this%leafp_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P shift storage to transfer', &
         ptr_patch=this%leafp_storage_to_xfer, default='inactive')

    this%frootp_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P shift storage to transfer', &
         ptr_patch=this%frootp_storage_to_xfer, default='inactive')

    this%livestemp_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P shift storage to transfer', &
         ptr_patch=this%livestemp_storage_to_xfer, default='inactive')

    this%deadstemp_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P shift storage to transfer', &
         ptr_patch=this%deadstemp_storage_to_xfer, default='inactive')

    this%livecrootp_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P shift storage to transfer', &
         ptr_patch=this%livecrootp_storage_to_xfer, default='inactive')

    this%deadcrootp_storage_to_xfer(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P shift storage to transfer', &
         ptr_patch=this%deadcrootp_storage_to_xfer, default='inactive')

    this%livestemp_to_deadstemp(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_TO_DEADSTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P turnover', &
         ptr_patch=this%livestemp_to_deadstemp, default='inactive')

    this%livestemp_to_retransp(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_TO_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P to retranslocated P pool', &
         ptr_patch=this%livestemp_to_retransp, default='inactive')

    this%livecrootp_to_deadcrootp(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_TO_DEADCROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P turnover', &
         ptr_patch=this%livecrootp_to_deadcrootp, default='inactive')

    this%livecrootp_to_retransp(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_TO_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P to retranslocated N pool', &
         ptr_patch=this%livecrootp_to_retransp, default='inactive')

    this%pdeploy(begp:endp) = spval
    call hist_addfld1d (fname='PDEPLOY', units='gP/m^2/s', &
         avgflag='A', long_name='total P deployed in new growth', &
         ptr_patch=this%pdeploy)

    this%wood_harvestp(begp:endp) = spval
    call hist_addfld1d (fname='WOOD_HARVESTP', units='gP/m^2/s', &
         avgflag='A', long_name='wood harvest P (to product pools)', &
         ptr_patch=this%wood_harvestp, default='inactive')

    this%fire_ploss(begp:endp) = spval
    call hist_addfld1d (fname='PFT_FIRE_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total pft-level fire P loss', &
         ptr_patch=this%fire_ploss, default='inactive')

    this%gap_ploss_litter(begp:endp) = spval
    call hist_addfld1d (fname='GAP_PLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total ploss from veg to litter due to gap mortality', &
         ptr_patch=this%gap_ploss_litter, default='inactive')

    this%fire_ploss_litter(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_PLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total ploss from veg to litter due to fire mortality', &
         ptr_patch=this%fire_ploss_litter, default='inactive')
    
    this%hrv_ploss_litter(begp:endp) = spval
    call hist_addfld1d (fname='HRV_PLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total ploss from veg to litter due to harvest mortality', &
         ptr_patch=this%hrv_ploss_litter, default='inactive')
    
    this%sen_ploss_litter(begp:endp) = spval
    call hist_addfld1d (fname='SEN_PLOSS_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='total ploss from veg to litter pool due to senescence', &
         ptr_patch=this%sen_ploss_litter, default='inactive')

    if (crop_prog) then
       this%fert_p(begp:endp) = spval
       call hist_addfld1d (fname='FERT_P', units='gP/m^2/s', &
            avgflag='A', long_name='fertilizer P added', &
            ptr_patch=this%fert_p)
    end if

    if (crop_prog) then
       this%fert_p_counter(begp:endp) = spval
       call hist_addfld1d (fname='FERT_P_COUNTER', units='seconds', &
            avgflag='A', long_name='time left to fertilize', &
            ptr_patch=this%fert_p_counter)
    end if

    this%crop_seedp_to_leaf(begp:endp) = spval
    call hist_addfld1d (fname='CROP_SEEDP_TO_LEAF', units='gP/m^2/s', &
         avgflag='A', long_name='crop seed source to leaf', &
         ptr_patch=this%crop_seedp_to_leaf, default='inactive')

    this%dwt_seedp_to_leaf(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_LEAF_PATCH', units='gP/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level leaf ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_seedp_to_leaf, default='inactive')

    this%dwt_seedp_to_deadstem(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_DEADSTEM_PATCH', units='gP/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level deadstem ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_seedp_to_deadstem, default='inactive')

    this%dwt_conv_pflux(begp:endp) = spval
    call hist_addfld1d (fname='DWT_CONV_PFLUX_PATCH', units='gP/m^2/s', &
         avgflag='A', &
         long_name='patch-level conversion C flux (immediate loss to atm) ' // &
         '(0 at all times except first timestep of year) ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_conv_pflux, default='inactive')

    this%dwt_prod10p_gain(begp:endp) = spval
    call hist_addfld1d (fname='DWT_PROD10P_GAIN_PATCH', units='gP/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_col=this%dwt_prod10p_gain, default='inactive')

    this%dwt_prod100p_gain(begp:endp) = spval
    call hist_addfld1d (fname='DWT_PROD100P_GAIN_PATCH', units='gP/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_col=this%dwt_prod100p_gain, default='inactive')

    this%dwt_seedp_to_ppool(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_PPOOL_PATCH', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level', &
         ptr_patch=this%dwt_seedp_to_ppool, default='inactive')

    this%plant_pdemand(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_PDEMAND', units='gP/m^2/s', &
         avgflag='A', long_name='P flux required to support initial GPP', &
         ptr_patch=this%plant_pdemand)

    this%avail_retransp(begp:endp) = spval
    call hist_addfld1d (fname='AVAIL_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='P flux available from retranslocation pool', &
         ptr_patch=this%avail_retransp, default='active')

    this%plant_palloc(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_PALLOC', units='gP/m^2/s', &
         avgflag='A', long_name='total allocated P flux', &
         ptr_patch=this%plant_palloc, default='active')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_pf
    !------------------------------------------------------------------------
    num_special_patch = 0
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do
    do p = begp,endp
       l = veg_pp%landunit(p)

       this%prev_leafp_to_litter (p)  = 0._r8 
       this%prev_frootp_to_litter(p)  = 0._r8 
     
       if ( crop_prog )then
          this%fert_p_counter(p)  = spval
          this%fert_p(p)          = 0._r8 
       end if

       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%fert_p_counter(p)  = 0._r8
       end if

       if (lun_pp%ifspecial(l)) then
          this%plant_pdemand(p)  = spval
          this%avail_retransp(p) = spval
          this%plant_palloc(p)   = spval
       end if
    end do
    
    call this%SetValues (num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8)

  end subroutine veg_pf_init
    
  !-----------------------------------------------------------------------
  subroutine veg_pf_restart (this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write restart data for vegetation-level phosphorus fluxes
    !
    ! !ARGUMENTS:
    class (vegetation_phosphorus_flux) :: this
    type(bounds_type) , intent(in)     :: bounds 
    type(file_desc_t) , intent(inout)  :: ncid   ! netcdf id
    character(len=*)  , intent(in)     :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !------------------------------------------------------------------------

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag, varname='fert_p_counter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_p_counter)

       call restartvar(ncid=ncid, flag=flag, varname='fert_p', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_p)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainp_xfer_to_grainp', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain P growth from storage', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_xfer_to_grainp)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='livestemp_to_litter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='livestem P to litter', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemp_to_litter)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainp_to_food', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain P to food', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_to_food)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='ppool_to_grainp', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain P', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%ppool_to_grainp)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='ppool_to_grainp_storage', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain P storage', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%ppool_to_grainp_storage)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag, varname='grainp_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain P shift storage to transfer', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_storage_to_xfer)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='plant_pdemand', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_pdemand) 

    call restartvar(ncid=ncid, flag=flag, varname='avail_retransp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%avail_retransp) 

    call restartvar(ncid=ncid, flag=flag, varname='plant_palloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_palloc) 
  
  end subroutine veg_pf_restart
  
  !-----------------------------------------------------------------------
  subroutine veg_pf_setvalues ( this, num_patch, filter_patch, value_patch)
    !
    ! !DESCRIPTION:
    ! Set phosphorus flux variables
    !
    ! !ARGUMENTS:
    class (vegetation_phosphorus_flux) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    !------------------------------------------------------------------------
    do fi = 1,num_patch
       i=filter_patch(fi)

       this%m_leafp_to_litter(i)                   = value_patch
       this%m_frootp_to_litter(i)                  = value_patch
       this%m_leafp_storage_to_litter(i)           = value_patch
       this%m_frootp_storage_to_litter(i)          = value_patch
       this%m_livestemp_storage_to_litter(i)       = value_patch
       this%m_deadstemp_storage_to_litter(i)       = value_patch
       this%m_livecrootp_storage_to_litter(i)      = value_patch
       this%m_deadcrootp_storage_to_litter(i)      = value_patch
       this%m_leafp_xfer_to_litter(i)              = value_patch
       this%m_frootp_xfer_to_litter(i)             = value_patch
       this%m_livestemp_xfer_to_litter(i)          = value_patch
       this%m_deadstemp_xfer_to_litter(i)          = value_patch
       this%m_livecrootp_xfer_to_litter(i)         = value_patch
       this%m_deadcrootp_xfer_to_litter(i)         = value_patch
       this%m_livestemp_to_litter(i)               = value_patch
       this%m_deadstemp_to_litter(i)               = value_patch
       this%m_livecrootp_to_litter(i)              = value_patch
       this%m_deadcrootp_to_litter(i)              = value_patch
       this%m_retransp_to_litter(i)                = value_patch
       this%m_ppool_to_litter(i)                   = value_patch
       this%hrv_leafp_to_litter(i)                 = value_patch             
       this%hrv_frootp_to_litter(i)                = value_patch            
       this%hrv_leafp_storage_to_litter(i)         = value_patch     
       this%hrv_frootp_storage_to_litter(i)        = value_patch    
       this%hrv_livestemp_storage_to_litter(i)     = value_patch 
       this%hrv_deadstemp_storage_to_litter(i)     = value_patch 
       this%hrv_livecrootp_storage_to_litter(i)    = value_patch
       this%hrv_deadcrootp_storage_to_litter(i)    = value_patch
       this%hrv_leafp_xfer_to_litter(i)            = value_patch        
       this%hrv_frootp_xfer_to_litter(i)           = value_patch       
       this%hrv_livestemp_xfer_to_litter(i)        = value_patch    
       this%hrv_deadstemp_xfer_to_litter(i)        = value_patch    
       this%hrv_livecrootp_xfer_to_litter(i)       = value_patch   
       this%hrv_deadcrootp_xfer_to_litter(i)       = value_patch   
       this%hrv_livestemp_to_litter(i)             = value_patch         
       this%hrv_deadstemp_to_prod10p(i)            = value_patch        
       this%hrv_deadstemp_to_prod100p(i)           = value_patch       
       this%hrv_leafp_to_prod1p(i)                 = value_patch
       this%hrv_livestemp_to_prod1p(i)             = value_patch
       this%hrv_grainp_to_prod1p(i)                = value_patch
       this%hrv_cropp_to_prod1p(i)                 = value_patch
       this%hrv_livecrootp_to_litter(i)            = value_patch        
       this%hrv_deadcrootp_to_litter(i)            = value_patch        
       this%hrv_retransp_to_litter(i)              = value_patch    
       this%hrv_ppool_to_litter(i)                 = value_patch

       this%m_leafp_to_fire(i)                     = value_patch
       this%m_leafp_storage_to_fire(i)             = value_patch
       this%m_leafp_xfer_to_fire(i)                = value_patch
       this%m_livestemp_to_fire(i)                 = value_patch
       this%m_livestemp_storage_to_fire(i)         = value_patch
       this%m_livestemp_xfer_to_fire(i)            = value_patch
       this%m_deadstemp_to_fire(i)                 = value_patch
       this%m_deadstemp_storage_to_fire(i)         = value_patch
       this%m_deadstemp_xfer_to_fire(i)            = value_patch
       this%m_frootp_to_fire(i)                    = value_patch
       this%m_frootp_storage_to_fire(i)            = value_patch
       this%m_frootp_xfer_to_fire(i)               = value_patch
       this%m_livecrootp_to_fire(i)                = value_patch
       this%m_livecrootp_storage_to_fire(i)        = value_patch
       this%m_livecrootp_xfer_to_fire(i)           = value_patch
       this%m_deadcrootp_to_fire(i)                = value_patch
       this%m_deadcrootp_storage_to_fire(i)        = value_patch
       this%m_deadcrootp_xfer_to_fire(i)           = value_patch
       this%m_retransp_to_fire(i)                  = value_patch
       this%m_ppool_to_fire(i)                     = value_patch

       this%m_leafp_to_litter_fire(i)              = value_patch
       this%m_leafp_storage_to_litter_fire(i)      = value_patch
       this%m_leafp_xfer_to_litter_fire(i)         = value_patch
       this%m_livestemp_to_litter_fire(i)          = value_patch
       this%m_livestemp_storage_to_litter_fire(i)  = value_patch
       this%m_livestemp_xfer_to_litter_fire(i)     = value_patch
       this%m_livestemp_to_deadstemp_fire(i)       = value_patch
       this%m_deadstemp_to_litter_fire(i)          = value_patch
       this%m_deadstemp_storage_to_litter_fire(i)  = value_patch
       this%m_deadstemp_xfer_to_litter_fire(i)     = value_patch
       this%m_frootp_to_litter_fire(i)             = value_patch
       this%m_frootp_storage_to_litter_fire(i)     = value_patch
       this%m_frootp_xfer_to_litter_fire(i)        = value_patch
       this%m_livecrootp_to_litter_fire(i)         = value_patch
       this%m_livecrootp_storage_to_litter_fire(i) = value_patch
       this%m_livecrootp_xfer_to_litter_fire(i)    = value_patch
       this%m_livecrootp_to_deadcrootp_fire(i)     = value_patch
       this%m_deadcrootp_to_litter_fire(i)         = value_patch
       this%m_deadcrootp_storage_to_litter_fire(i) = value_patch
       this%m_deadcrootp_xfer_to_litter_fire(i)    = value_patch
       this%m_retransp_to_litter_fire(i)           = value_patch
       this%m_ppool_to_litter_fire(i)              = value_patch

       this%leafp_xfer_to_leafp(i)                 = value_patch
       this%frootp_xfer_to_frootp(i)               = value_patch
       this%livestemp_xfer_to_livestemp(i)         = value_patch
       this%deadstemp_xfer_to_deadstemp(i)         = value_patch
       this%livecrootp_xfer_to_livecrootp(i)       = value_patch
       this%deadcrootp_xfer_to_deadcrootp(i)       = value_patch
       this%leafp_to_litter(i)                     = value_patch
       this%leafp_to_retransp(i)                   = value_patch
       this%frootp_to_litter(i)                    = value_patch
       this%retransp_to_ppool(i)                   = value_patch
       this%sminp_to_ppool(i)                      = value_patch
       this%ppool_to_leafp(i)                      = value_patch
       this%ppool_to_leafp_storage(i)              = value_patch
       this%ppool_to_frootp(i)                     = value_patch
       this%ppool_to_frootp_storage(i)             = value_patch
       this%ppool_to_livestemp(i)                  = value_patch
       this%ppool_to_livestemp_storage(i)          = value_patch
       this%ppool_to_deadstemp(i)                  = value_patch
       this%ppool_to_deadstemp_storage(i)          = value_patch
       this%ppool_to_livecrootp(i)                 = value_patch
       this%ppool_to_livecrootp_storage(i)         = value_patch
       this%ppool_to_deadcrootp(i)                 = value_patch
       this%ppool_to_deadcrootp_storage(i)         = value_patch
       this%leafp_storage_to_xfer(i)               = value_patch
       this%frootp_storage_to_xfer(i)              = value_patch
       this%livestemp_storage_to_xfer(i)           = value_patch
       this%deadstemp_storage_to_xfer(i)           = value_patch
       this%livecrootp_storage_to_xfer(i)          = value_patch
       this%deadcrootp_storage_to_xfer(i)          = value_patch
       this%livestemp_to_deadstemp(i)              = value_patch
       this%livestemp_to_retransp(i)               = value_patch
       this%livecrootp_to_deadcrootp(i)            = value_patch
       this%livecrootp_to_retransp(i)              = value_patch
       this%pdeploy(i)                             = value_patch
       this%wood_harvestp(i)                       = value_patch
       this%fire_ploss(i)                          = value_patch
       this%biochem_pmin_to_plant(i)               = value_patch
       this%gap_ploss_litter(i)                    = value_patch
       this%fire_ploss_litter(i)                   = value_patch
       this%hrv_ploss_litter(i)                    = value_patch
       this%sen_ploss_litter(i)                    = value_patch
       this%livestemp_to_litter(i)                 = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grainp_to_food(i)                   = value_patch
          this%grainp_xfer_to_grainp(i)            = value_patch
          this%ppool_to_grainp(i)                  = value_patch
          this%ppool_to_grainp_storage(i)          = value_patch
          this%grainp_storage_to_xfer(i)           = value_patch
          this%frootp_to_retransp(i)               = value_patch
          this%crop_seedp_to_leaf(i)               = value_patch
       end do
    end if
  
  end subroutine veg_pf_setvalues
  
  !-----------------------------------------------------------------------
  subroutine veg_pf_summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_pf)
    !
    ! !ARGUMENTS:
    class (vegetation_phosphorus_flux) :: this
    type(bounds_type) , intent(in)     :: bounds  
    integer           , intent(in)     :: num_soilc       ! number of soil columns in filter
    integer           , intent(in)     :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in)     :: num_soilp       ! number of soil patches in filter
    integer           , intent(in)     :: filter_soilp(:) ! filter for soil patches
    type(column_phosphorus_flux), intent(inout) :: col_pf
    !
    ! !LOCAL VARIABLES:
    integer  :: p, fp   ! indices
    !-----------------------------------------------------------------------
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! total P deployment (from sminn and retranslocated P pool) (PDEPLOY)
       this%pdeploy(p) = &
            this%sminp_to_ppool(p) + &
            this%retransp_to_ppool(p)

       ! pft-level wood harvest
       this%wood_harvestp(p) = &
            this%hrv_deadstemp_to_prod10p(p) + &
            this%hrv_deadstemp_to_prod100p(p)
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            this%wood_harvestp(p) = &
            this%wood_harvestp(p) + &
            this%hrv_cropp_to_prod1p(p)
       end if

       ! total pft-level fire P losses
       this%fire_ploss(p) = &
            this%m_leafp_to_fire(p)               + &
            this%m_leafp_storage_to_fire(p)       + &
            this%m_leafp_xfer_to_fire(p)          + &
            this%m_frootp_to_fire(p)              + &
            this%m_frootp_storage_to_fire(p)      + &
            this%m_frootp_xfer_to_fire(p)         + &
            this%m_livestemp_to_fire(p)           + &
            this%m_livestemp_storage_to_fire(p)   + &
            this%m_livestemp_xfer_to_fire(p)      + &
            this%m_deadstemp_to_fire(p)           + &
            this%m_deadstemp_storage_to_fire(p)   + &
            this%m_deadstemp_xfer_to_fire(p)      + &
            this%m_livecrootp_to_fire(p)          + &
            this%m_livecrootp_storage_to_fire(p)  + &
            this%m_livecrootp_xfer_to_fire(p)     + &
            this%m_deadcrootp_to_fire(p)          + &
            this%m_deadcrootp_storage_to_fire(p)  + &
            this%m_deadcrootp_xfer_to_fire(p)     + &
            this%m_retransp_to_fire(p)            + &
            this%m_ppool_to_fire(p)

      this%gap_ploss_litter(p) = &
           this%m_leafp_to_litter(p)              + &
           this%m_leafp_storage_to_litter(p)      + &
           this%m_leafp_xfer_to_litter(p)         + &
           this%m_frootp_to_litter(p)             + &
           this%m_frootp_storage_to_litter(p)     + &
           this%m_frootp_xfer_to_litter(p)        + &
           this%m_livestemp_to_litter(p)          + &
           this%m_livestemp_storage_to_litter(p)  + &
           this%m_livestemp_xfer_to_litter(p)     + &
           this%m_deadstemp_to_litter(p)          + &
           this%m_deadstemp_storage_to_litter(p)  + &
           this%m_deadstemp_xfer_to_litter(p)     + &
           this%m_livecrootp_to_litter(p)         + &
           this%m_livecrootp_storage_to_litter(p) + &
           this%m_livecrootp_xfer_to_litter(p)    + &
           this%m_deadcrootp_to_litter(p)         + &
           this%m_deadcrootp_storage_to_litter(p) + &
           this%m_deadcrootp_xfer_to_litter(p)    + &
           this%m_retransp_to_litter(p)           + &
           this%m_ppool_to_litter(p)

      this%fire_ploss_litter(p) = &
           this%m_deadstemp_to_litter_fire(p)     + &
           this%m_deadcrootp_to_litter_fire(p)    + &
           this%m_retransp_to_litter_fire(p)      + &
           this%m_ppool_to_litter_fire(p)         + &
           this%m_leafp_to_litter_fire(p)         + &
           this%m_frootp_to_litter_fire(p)        + &
           this%m_livestemp_to_litter_fire(p)     + &
           this%m_livecrootp_to_litter_fire(p)    + &
           this%m_leafp_storage_to_litter_fire(p) + &
           this%m_frootp_storage_to_litter_fire(p)       + &
           this%m_livestemp_storage_to_litter_fire(p)    + &
           this%m_deadstemp_storage_to_litter_fire(p)    + &
           this%m_livecrootp_storage_to_litter_fire(p)   + &
           this%m_deadcrootp_storage_to_litter_fire(p)   + &
           this%m_leafp_xfer_to_litter_fire(p)           + &
           this%m_frootp_xfer_to_litter_fire(p)          + &
           this%m_livestemp_xfer_to_litter_fire(p)       + &
           this%m_deadstemp_xfer_to_litter_fire(p)       + &
           this%m_livecrootp_xfer_to_litter_fire(p)      + &
           this%m_deadcrootp_xfer_to_litter_fire(p)

      this%hrv_ploss_litter(p) = &
           this%hrv_retransp_to_litter(p)         + &
           this%hrv_ppool_to_litter(p)            + &
           this%hrv_leafp_to_litter(p)            + &
           this%hrv_leafp_storage_to_litter(p)    + &
           this%hrv_leafp_xfer_to_litter(p)       + &
           this%hrv_frootp_to_litter(p)           + &
           this%hrv_frootp_storage_to_litter(p)   + &
           this%hrv_frootp_xfer_to_litter(p)      + &
           this%hrv_livestemp_to_litter(p)        + &
           this%hrv_livestemp_storage_to_litter(p)+ &
           this%hrv_livestemp_xfer_to_litter(p)   + &
           this%hrv_deadstemp_storage_to_litter(p)+ &
           this%hrv_deadstemp_xfer_to_litter(p)   + &
           this%hrv_livecrootp_to_litter(p)       + &
           this%hrv_livecrootp_storage_to_litter(p)+ &
           this%hrv_livecrootp_xfer_to_litter(p)  + &
           this%hrv_deadcrootp_to_litter(p)       + &
           this%hrv_deadcrootp_storage_to_litter(p)+ &
           this%hrv_deadcrootp_xfer_to_litter(p)

      if (crop_prog) then
         this%sen_ploss_litter(p) = &
             this%livestemp_to_litter(p)            + &
             this%leafp_to_litter(p)                + &
             this%frootp_to_litter(p)
      else
         this%sen_ploss_litter(p) = &
             this%leafp_to_litter(p)                + &
             this%frootp_to_litter(p)
      end if
    end do

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_ploss(bounds%begp:bounds%endp), &
         col_pf%fire_ploss_p2c(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%wood_harvestp(bounds%begp:bounds%endp), &
         col_pf%wood_harvestp(bounds%begc:bounds%endc))
  
  end subroutine veg_pf_summary
  
  !------------------------------------------------------------------------
  subroutine veg_pf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    !-----------------------------------------------------------------------
   
  end subroutine veg_pf_clean
  
    !------------------------------------------------------------------------
    
end module vegetationDataType
