module zm_conv_intr
   !----------------------------------------------------------------------------
   ! 
   ! Interface to the Zhang-McFarlane deep convection scheme
   !
   ! Author: D.B. Coleman
   ! January 2010 modified by J. Kay to add COSP simulator fields to pbuf
   ! July 2015 B. Singh Added code for unified convective trasport
   !----------------------------------------------------------------------------
   use shr_kind_mod,          only: r8=>shr_kind_r8
   use spmd_utils,            only: masterproc
   use perf_mod,              only: t_startf, t_stopf
   use cam_abortutils,        only: endrun
   use cam_history,           only: outfld, addfld, horiz_only, add_default
   use cam_logfile,           only: iulog
   use physconst,             only: pi, cpair
   use ppgrid,                only: pver, pcols, pverp, begchunk, endchunk
   use rad_constituents,      only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr
   use rad_constituents,      only: rad_cnst_get_aer_props, rad_cnst_get_mode_props
   use ndrop_bam,             only: ndrop_bam_init
   use zm_conv,               only: zm_conv_evap, zm_convr, convtran, momtran, trigdcape_ull, trig_dcape_only
   use zm_conv,               only: MCSP, MCSP_heat_coeff, MCSP_moisture_coeff, MCSP_uwind_coeff, MCSP_vwind_coeff
   use zm_conv,               only: zm_microp
   use zm_microphysics,       only: zm_aero_t
   use zm_microphysics_state, only: zm_microp_st

   implicit none
   private
   save

   ! public methods
   public :: zm_conv_register  ! register fields in physics buffer
   public :: zm_conv_init      ! initialize donner_deep module
   public :: zm_conv_tend      ! return tendencies
   public :: zm_conv_tend_2    ! return tendencies

   ! physics buffer field indices
   integer :: dp_flxprc_idx    ! deep conv flux of precipitation from deep convection (kg/m2/s)
   integer :: dp_flxsnw_idx    ! deep conv flux of snow from deep convection (kg/m2/s)
   integer :: dp_cldliq_idx    ! deep conv cloud liq water (kg/kg)
   integer :: dp_cldice_idx    ! deep conv cloud ice water (kg/kg)
   integer :: dlfzm_idx        ! detrained convective cloud water mixing ratio
   integer :: difzm_idx        ! detrained convective cloud ice mixing ratio
   integer :: dsfzm_idx        ! detrained convective snow mixing ratio
   integer :: dnlfzm_idx       ! detrained convective cloud water num concen
   integer :: dnifzm_idx       ! detrained convective cloud ice num concen
   integer :: dnsfzm_idx       ! detrained convective snow num concen
   integer :: prec_dp_idx      ! total surface precipitation rate from deep conv
   integer :: snow_dp_idx      ! frozen surface precipitation rate from deep conv
   integer :: wuc_idx          ! vertical velocity in deep convection
   integer :: t_star_idx       ! DCAPE temperature from previous time step
   integer :: q_star_idx       ! DCAPE water vapor from previous time step
   integer :: cld_idx          ! cloud fraction
   integer :: icwmrdp_idx      ! in-cloud water mixing ratio
   integer :: rprddp_idx       ! rain production
   integer :: fracis_idx       ! fraction of transported species that are insoluble
   integer :: nevapr_dpcu_idx  ! evaporation of deep conv precipitation
   integer :: dgnum_idx        ! dry aerosol mode diameter
   integer :: lambdadpcu_idx   ! slope of cloud liquid size distribution
   integer :: mudpcu_idx       ! width parameter of droplet size distr
   integer :: icimrdp_idx      ! in-cloud ice mixing ratio

   ! other private module data
   logical :: convproc_do_aer 
   logical :: convproc_do_gas 
   logical :: clim_modal_aero
   logical :: old_snow  = .true.   ! flag to use old estimate of snow production in zm_conv_evap
                                   ! set false to use snow production from zm microphysics
   integer :: nmodes
   integer :: nbulk
   type(zm_aero_t), allocatable :: aero(:) ! object contains aerosol information
   
   real(r8), parameter :: ZM_upper_limit_pref   = 40e2_r8  ! pressure limit above which deep convection is not allowed [Pa]
   real(r8), parameter :: MCSP_storm_speed_pref = 600e2_r8 ! pressure level for winds in MCSP calculation [Pa]
   real(r8), parameter :: MCSP_conv_depth_min   = 700e2_r8 ! pressure thickness of convective heating [Pa]
   real(r8), parameter :: MCSP_shear_min        = 3.0_r8   ! min shear value for MCSP to be active
   real(r8), parameter :: MCSP_shear_max        = 200.0_r8 ! max shear value for MCSP to be active

contains

!===================================================================================================

subroutine zm_conv_register
   !----------------------------------------------------------------------------
   ! Purpose: register fields with the physics buffer
   !----------------------------------------------------------------------------
   use physics_buffer, only : pbuf_add_field, dtype_r8
   use misc_diagnostics,only: dcape_diags_register
   implicit none

   integer idx

   ! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx)

   ! Flux of snow from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx)

   ! deep conv cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)

   ! deep conv cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)

   ! vertical velocity (m/s)
   call pbuf_add_field('WUC','global',dtype_r8,(/pcols,pver/), wuc_idx)

   ! previous time step data for DCAPE calculation
   if (trigdcape_ull .or. trig_dcape_only) then
      call pbuf_add_field('T_STAR','global',dtype_r8,(/pcols,pver/), t_star_idx)
      call pbuf_add_field('Q_STAR','global',dtype_r8,(/pcols,pver/), q_star_idx)
   end if

   ! detrained convective cloud water mixing ratio.
   call pbuf_add_field('DLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dlfzm_idx)
   ! detrained convective cloud ice mixing ratio.
   call pbuf_add_field('DIFZM', 'physpkg', dtype_r8, (/pcols,pver/), difzm_idx)

   ! Only add the number conc fields if the microphysics is active.
   if (zm_microp) then
      ! detrained convective cloud water num concen.
      call pbuf_add_field('DNLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnlfzm_idx)
      ! detrained convective cloud ice num concen.
      call pbuf_add_field('DNIFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnifzm_idx)
      ! detrained convective snow num concen.
      call pbuf_add_field('DNSFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnsfzm_idx)
      ! detrained convective snow mixing ratio.
      call pbuf_add_field('DSFZM',  'physpkg', dtype_r8, (/pcols,pver/), dsfzm_idx)
   end if

   ! Register variables for dCAPE diagnosis and decomposition
   call dcape_diags_register( pcols )

end subroutine zm_conv_register

!===================================================================================================

subroutine zm_conv_init(pref_edge)
   !----------------------------------------------------------------------------
   ! Purpose:  declare output fields, initialize variables needed by convection
   !----------------------------------------------------------------------------
   use zm_conv,            only: zm_convi
   use pmgrid,             only: plev,plevp
   use spmd_utils,         only: masterproc
   use error_messages,     only: alloc_err
   use phys_control,       only: phys_deepconv_pbl, phys_getopts
   use physics_buffer,     only: pbuf_get_index
   use rad_constituents,   only: rad_cnst_get_info
   use zm_microphysics,    only: zm_mphyi

   implicit none

   real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces

   logical :: no_deep_pbl                 ! if true, no deep convection in PBL
   integer :: limcnv                      ! top interface level limit for convection
   logical :: history_budget              ! output tendencies and state variables for 
                                          ! temperature, water vapor, cloud ice/liq budgets
   integer :: history_budget_histfile_num ! output history file number for budget fields
   integer i, k, istat
   character(len=*), parameter :: routine = 'zm_conv_init'

   ! Allocate the basic aero structure outside the zmconv_microp logical
   ! This allows the aero structure to be passed
   ! Note that all of the arrays inside this structure are conditionally allocated
   allocate(aero(begchunk:endchunk))

   ! Register fields with the output buffer

   call addfld('PRECZ',        horiz_only, 'A', 'm/s',      'ZM total precipitation rate')
   call addfld('ZMDT',         (/ 'lev'/), 'A', 'K/s',      'ZM T tendency')
   call addfld('ZMDQ',         (/ 'lev'/), 'A', 'kg/kg/s',  'ZM Q tendency')
   call addfld('ZMDICE',       (/ 'lev'/), 'A', 'kg/kg/s',  'ZM cloud ice tendency')
   call addfld('ZMDLIQ',       (/ 'lev'/), 'A', 'kg/kg/s',  'ZM cloud liq tendency')
   call addfld('EVAPTZM',      (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from evaporation/snow prod')
   call addfld('FZSNTZM',      (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from rain to snow conversion')
   call addfld('EVSNTZM',      (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from snow to rain prod')
   call addfld('EVAPQZM',      (/ 'lev'/), 'A', 'kg/kg/s',  'ZM Q tendency from evaporation')
   call addfld('ZMFLXPRC',     (/'ilev'/), 'A', 'kg/m2/s',  'ZM flux of precipitation')
   call addfld('ZMFLXSNW',     (/'ilev'/), 'A', 'kg/m2/s',  'ZM flux of snow')
   call addfld('ZMNTPRPD',     (/ 'lev'/), 'A', 'kg/kg/s',  'ZM net precipitation production')
   call addfld('ZMNTSNPD',     (/ 'lev'/), 'A', 'kg/kg/s',  'ZM net snow production')
   call addfld('ZMEIHEAT',     (/ 'lev'/), 'A', 'W/kg',     'ZM heating by ice and evaporation')
   call addfld('CMFMCDZM',     (/'ilev'/), 'A', 'kg/m2/s',  'ZM convection mass flux')
   call addfld('PRECCDZM',     horiz_only, 'A', 'm/s',      'ZM convective precipitation rate')
   call addfld('PCONVB',       horiz_only, 'A', 'Pa',       'ZM convection base pressure')
   call addfld('PCONVT',       horiz_only, 'A', 'Pa',       'ZM convection top  pressure')
   call addfld('MAXI',         horiz_only, 'A', 'level',    'ZM model level of launching parcel')
   call addfld('CAPE_ZM',      horiz_only, 'A', 'J/kg',     'ZM convectively available potential energy')
   call addfld('DCAPE',        horiz_only, 'A', 'J/kg',     'ZM rate of change of CAPE')
   call addfld('FREQZM',       horiz_only, 'A', 'fraction', 'ZM fractional occurrence of convection')
   call addfld('ZMMTT',        (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from convective momentum transport')
   call addfld('ZMMTU',        (/ 'lev'/), 'A', 'm/s2',     'ZM U tendency from convective momentum transport')
   call addfld('ZMMTV',        (/ 'lev'/), 'A', 'm/s2',     'ZM V tendency from convective momentum transport')
   call addfld('ZMMU',         (/ 'lev'/), 'A', 'kg/m2/s',  'ZM convection updraft mass flux')
   call addfld('ZMMD',         (/ 'lev'/), 'A', 'kg/m2/s',  'ZM convection downdraft mass flux')
   call addfld('ZMUPGU',       (/ 'lev'/), 'A', 'm/s2',     'ZM zonal force from updraft pressure gradient term')
   call addfld('ZMUPGD',       (/ 'lev'/), 'A', 'm/s2',     'ZM zonal force from downdraft pressure gradient term')
   call addfld('ZMVPGU',       (/ 'lev'/), 'A', 'm/s2',     'ZM meridional force from updraft pressure gradient term')
   call addfld('ZMVPGD',       (/ 'lev'/), 'A', 'm/s2',     'ZM merdional force from downdraft pressure gradient term')
   call addfld('ZMICUU',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud U updrafts')
   call addfld('ZMICUD',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud U downdrafts')
   call addfld('ZMICVU',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud V updrafts')
   call addfld('ZMICVD',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud V downdrafts')

   if (MCSP) then 
      call addfld('MCSP_DT',   (/ 'lev'/), 'A', 'K/s',      'MCSP T tendency')
      call addfld('MCSP_freq', horiz_only, 'A', '1',        'MCSP frequency of activation')
      call addfld('MCSP_DU',   (/ 'lev'/), 'A', 'm/s/day',  'MCSP U tendency')
      call addfld('MCSP_DV',   (/ 'lev'/), 'A', 'm/s/day',  'MCSP V tendency')
      call addfld('MCSP_shear',horiz_only, 'A', 'm/s',      'MCSP vertical zonal wind shear')
      call addfld('ZM_depth',  horiz_only, 'A', 'Pa',       'ZM convection depth')
   end if

   if (zm_microp) then

      call addfld ('CLDLIQZM',(/ 'lev' /), 'A', 'g/m3',     'ZM cloud liq water')
      call addfld ('CLDICEZM',(/ 'lev' /), 'A', 'g/m3',     'ZM cloud ice water')
      call addfld ('CLIQSNUM',(/ 'lev' /), 'A', '1',        'ZM cloud liq water sample number')
      call addfld ('CICESNUM',(/ 'lev' /), 'A', '1',        'ZM cloud ice water sample number')
      call addfld ('QRAINZM' ,(/ 'lev' /), 'A', 'g/m3',     'ZM rain water')
      call addfld ('QSNOWZM' ,(/ 'lev' /), 'A', 'g/m3',     'ZM snow')
      call addfld ('QGRAPZM' ,(/ 'lev' /), 'A', 'g/m3',     'ZM graupel')
      call addfld ('CRAINNUM',(/ 'lev' /), 'A', '1',        'ZM cloud rain water sample number')
      call addfld ('CSNOWNUM',(/ 'lev' /), 'A', '1',        'ZM cloud snow sample number')
      call addfld ('CGRAPNUM',(/ 'lev' /), 'A', '1',        'ZM cloud graupel sample number')

      call addfld ('DIFZM',   (/ 'lev' /), 'A', 'kg/kg/s ', 'ZM detrained ice water')
      call addfld ('DLFZM',   (/ 'lev' /), 'A', 'kg/kg/s ', 'ZM detrained liq water')
      call addfld ('DNIFZM',  (/ 'lev' /), 'A', '1/kg/s ',  'ZM detrained ice water num concen')
      call addfld ('DNLFZM',  (/ 'lev' /), 'A', '1/kg/s ',  'ZM detrained liquid water num concen')
      call addfld ('WUZM',    (/ 'lev' /), 'A', 'm/s',      'ZM vertical velocity')
      call addfld ('WUZMSNUM',(/ 'lev' /), 'A', '1',        'ZM vertical velocity sample number')

      call addfld ('QNLZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud liq water number concen')
      call addfld ('QNIZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud ice number concen')
      call addfld ('QNRZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud rain water number concen')
      call addfld ('QNSZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud snow number concen')
      call addfld ('QNGZM',   (/ 'lev' /), 'A', '1/m3',     'ZM cloud graupel number concen')

      call addfld ('FRZZM',   (/ 'lev' /), 'A', 'K/s',      'ZM heating tendency due to freezing')

      call addfld ('AUTOL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to autoconversion of droplets to rain')
      call addfld ('ACCRL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to accretion of droplets by rain')
      call addfld ('BERGN_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to Bergeron process')
      call addfld ('FHTIM_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to immersion freezing')
      call addfld ('FHTCT_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to contact freezing')
      call addfld ('FHML_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to homogeneous freezing of droplet')
      call addfld ('HMPI_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to HM process')
      call addfld ('ACCSL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to accretion of droplet by snow')
      call addfld ('DLF_M',   (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to detrainment of droplet')
      call addfld ('COND_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to condensation')

      call addfld ('AUTOL_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to autoconversion of droplets to rain')
      call addfld ('ACCRL_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to accretion of droplets by rain')
      call addfld ('BERGN_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to Bergeron process')
      call addfld ('FHTIM_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to immersion freezing')
      call addfld ('FHTCT_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to contact freezing')
      call addfld ('FHML_N',  (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to homogeneous freezing of droplet')
      call addfld ('ACCSL_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to accretion of droplet by snow')
      call addfld ('ACTIV_N', (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to droplets activation')
      call addfld ('DLF_N',   (/ 'lev' /), 'A', '1/kg/m',   'ZM num tendency due to detrainment of droplet')

      call addfld ('AUTOI_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to autoconversion of ice to snow')
      call addfld ('ACCSI_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to accretion of ice by snow')
      call addfld ('DIF_M',   (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to detrainment of cloud ice')
      call addfld ('DEPOS_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to deposition')

      call addfld ('NUCLI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to ice nucleation')
      call addfld ('AUTOI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to autoconversion of ice to snow')
      call addfld ('ACCSI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to accretion of ice by snow')
      call addfld ('HMPI_N',  (/ 'lev' /), 'A', '1/kg/s' ,  'ZM num tendency due to HM process')
      call addfld ('DIF_N',   (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to detrainment of cloud ice')
      call addfld ('TRSPC_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of droplets due to convective transport')
      call addfld ('TRSPC_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of droplets due to convective transport')
      call addfld ('TRSPI_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of ice crystal due to convective transport')
      call addfld ('TRSPI_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of ice crystal due to convective transport')

      call addfld ('ACCGR_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to collection of rain by graupel')
      call addfld ('ACCGL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to collection of droplets by graupel')
      call addfld ('ACCGSL_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of droplets by snow')
      call addfld ('ACCGSR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of rain by snow')
      call addfld ('ACCGIR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of rain by ice')
      call addfld ('ACCGRI_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of ice by rain')
      call addfld ('ACCGRS_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel due to collection of snow by rain')

      call addfld ('ACCGSL_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel due to collection of droplets by snow')
      call addfld ('ACCGSR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel due to collection of rain by snow')
      call addfld ('ACCGIR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel due to collection of rain by ice')

      call addfld ('ACCSRI_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of snow due to collection of ice by rain')
      call addfld ('ACCIGL_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of ice mult(splintering) due to acc droplets by graupel')
      call addfld ('ACCIGR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of ice mult(splintering) due to acc rain by graupel')
      call addfld ('ACCSIR_M',(/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of snow due to collection of rain by ice')

      call addfld ('ACCIGL_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of ice mult(splintering) due to acc droplets by graupel')
      call addfld ('ACCIGR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of ice mult(splintering) due to acc rain by graupel')
      call addfld ('ACCSIR_N',(/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of snow due to collection of rain by ice')
      call addfld ('ACCGL_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to collection of droplets by graupel')
      call addfld ('ACCGR_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency due to collection of rain by graupel')
      call addfld ('ACCIL_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of cloud ice due to collection of droplet by cloud ice')
      call addfld ('ACCIL_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of cloud ice due to collection of droplet by cloud ice')

      call addfld ('FALLR_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of rain fallout')
      call addfld ('FALLS_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of snow fallout')
      call addfld ('FALLG_M', (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency of graupel fallout')
      call addfld ('FALLR_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of rain fallout')
      call addfld ('FALLS_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of snow fallout')
      call addfld ('FALLG_N', (/ 'lev' /), 'A', '1/kg/m' ,  'ZM num tendency of graupel fallout')
      call addfld ('FHMR_M',  (/ 'lev' /), 'A', 'kg/kg/m',  'ZM mass tendency due to homogeneous freezing of rain')

      call addfld ('PRECZ_SN',horiz_only , 'A', '#',        'ZM sample num of convective precipitation rate')

      call add_default ('CLDLIQZM', 1, ' ')
      call add_default ('CLDICEZM', 1, ' ')
      call add_default ('CLIQSNUM', 1, ' ')
      call add_default ('CICESNUM', 1, ' ')
      call add_default ('DIFZM',    1, ' ')
      call add_default ('DLFZM',    1, ' ')
      call add_default ('DNIFZM',   1, ' ')
      call add_default ('DNLFZM',   1, ' ')
      call add_default ('WUZM',     1, ' ')
      call add_default ('QRAINZM',  1, ' ')
      call add_default ('QSNOWZM',  1, ' ')
      call add_default ('QGRAPZM',  1, ' ')
      call add_default ('CRAINNUM', 1, ' ')
      call add_default ('CSNOWNUM', 1, ' ')
      call add_default ('CGRAPNUM', 1, ' ')
      call add_default ('QNLZM',    1, ' ')
      call add_default ('QNIZM',    1, ' ')
      call add_default ('QNRZM',    1, ' ')
      call add_default ('QNSZM',    1, ' ')
      call add_default ('QNGZM',    1, ' ')
      call add_default ('FRZZM',    1, ' ')

   end if

   call phys_getopts( history_budget_out = history_budget, &
                      history_budget_histfile_num_out = history_budget_histfile_num, &
                      convproc_do_aer_out = convproc_do_aer, &
                      convproc_do_gas_out = convproc_do_gas)
   
   ! Determine whether modal aerosols are used
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   if ( history_budget ) then
      call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
      call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
      call add_default('ZMDT     ', history_budget_histfile_num, ' ')
      call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
      call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
      call add_default('ZMDICE   ', history_budget_histfile_num, ' ')
      call add_default('ZMMTT    ', history_budget_histfile_num, ' ')
   end if

   ! Limit deep convection to regions below ZM_upper_limit_pref
   limcnv = 0 ! initialize to null value to check against below
   if (pref_edge(1) >= ZM_upper_limit_pref) then
      limcnv = 1
   else
      do k = 1,plev
         if (pref_edge(k)   <  ZM_upper_limit_pref .and. &
             pref_edge(k+1) >= ZM_upper_limit_pref) then
            limcnv = k
            exit
         end if
      end do
      if ( limcnv == 0 ) limcnv = plevp
   end if
    
   if (masterproc) write(iulog,*)'ZM_CONV_INIT: Deep convection will be capped at ', &
                                 'intfc ',limcnv,' which is ',pref_edge(limcnv),' pascals'

   no_deep_pbl = phys_deepconv_pbl()
   call zm_convi( limcnv, no_deep_pbl_in=no_deep_pbl )

   dp_flxprc_idx   = pbuf_get_index('DP_FLXPRC')
   dp_flxsnw_idx   = pbuf_get_index('DP_FLXSNW')
   dp_cldliq_idx   = pbuf_get_index('DP_CLDLIQ')
   dp_cldice_idx   = pbuf_get_index('DP_CLDICE')
   cld_idx         = pbuf_get_index('CLD')
   icwmrdp_idx     = pbuf_get_index('ICWMRDP')
   rprddp_idx      = pbuf_get_index('RPRDDP')
   fracis_idx      = pbuf_get_index('FRACIS')
   nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
   prec_dp_idx     = pbuf_get_index('PREC_DP')
   snow_dp_idx     = pbuf_get_index('SNOW_DP')
   wuc_idx         = pbuf_get_index('WUC')
   lambdadpcu_idx  = pbuf_get_index('LAMBDADPCU')
   mudpcu_idx      = pbuf_get_index('MUDPCU')
   icimrdp_idx     = pbuf_get_index('ICIMRDP')

   ! Initialization for the microphysics
   if (zm_microp) then

      call zm_mphyi()

      ! use old estimate of snow production in zm_conv_evap
      old_snow = .false. 

      ! Initialize the aerosol object with data from the modes/species
      ! affecting climate, i.e., the list index is hardcoded to 0
      call rad_cnst_get_info( 0, nmodes=nmodes, naero=nbulk )

      do i = begchunk, endchunk
         call zm_aero_init(nmodes, nbulk, aero(i))
      end do

      if (nmodes > 0) then
         dgnum_idx = pbuf_get_index('DGNUM')
      else if (nbulk > 0) then
         call ndrop_bam_init()
      end if

   end if ! zmconv_microp

   !----------------------------------------------------------------------------
   contains
   subroutine zm_aero_init(nmodes, nbulk, aero)
      ! Initialize the zm_aero_t object for modal aerosols
      integer,         intent(in)  :: nmodes
      integer,         intent(in)  :: nbulk
      type(zm_aero_t), intent(out) :: aero
      integer :: iaer, l, m
      integer :: nspecmx   ! max number of species in a mode
      character(len=20), allocatable :: aername(:)
      character(len=32) :: str32
      real(r8) :: sigmag, dgnumlo, dgnumhi
      real(r8) :: alnsg
      !-------------------------------------------------------------------------
      aero%nmodes = nmodes
      aero%nbulk  = nbulk

      if (nmodes > 0) then
         ! Initialize the modal aerosol information
         aero%scheme = 'modal'

         ! Get number of species in each mode, and find max.
         allocate(aero%nspec(aero%nmodes))
         nspecmx = 0
         do m = 1, aero%nmodes
            call rad_cnst_get_info(0, m, nspec=aero%nspec(m), mode_type=str32)
            nspecmx = max(nspecmx, aero%nspec(m))
            ! save mode index for specified mode types
            select case (trim(str32))
            case ('accum')
               aero%mode_accum_idx = m
            case ('aitken')
               aero%mode_aitken_idx = m
            case ('coarse')
               aero%mode_coarse_idx = m
            end select
         end do

         ! Check that required mode types were found
         if (aero%mode_accum_idx == -1 .or. &
             aero%mode_aitken_idx == -1 .or. &
             aero%mode_coarse_idx == -1) then
            write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
               aero%mode_accum_idx, aero%mode_aitken_idx, aero%mode_coarse_idx
            call endrun(routine//': ERROR required mode type not found')
         end if

         ! find indices for the dust and seasalt species in the coarse mode
         do l = 1, aero%nspec(aero%mode_coarse_idx)
            call rad_cnst_get_info(0, aero%mode_coarse_idx, l, spec_type=str32)
            select case (trim(str32))
            case ('dust')
               aero%coarse_dust_idx = l
            case ('seasalt')
               aero%coarse_nacl_idx = l
            end select
         end do

         ! Check that required modal species types were found
         if (aero%coarse_dust_idx == -1 .or. &
             aero%coarse_nacl_idx == -1) then
            write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
               aero%coarse_dust_idx, aero%coarse_nacl_idx
            call endrun(routine//': ERROR required mode-species type not found')
         end if

         allocate( &
            aero%num_a(nmodes), &
            aero%mmr_a(nspecmx,nmodes), &
            aero%numg_a(pcols,pver,nmodes), &
            aero%mmrg_a(pcols,pver,nspecmx,nmodes), &
            aero%voltonumblo(nmodes), &
            aero%voltonumbhi(nmodes), &
            aero%specdens(nspecmx,nmodes), &
            aero%spechygro(nspecmx,nmodes), &
            aero%dgnum(nmodes), &
            aero%dgnumg(pcols,pver,nmodes) )

         do m = 1, nmodes

            ! Properties of modes
            call rad_cnst_get_mode_props( 0, m, sigmag=sigmag, dgnumlo=dgnumlo, dgnumhi=dgnumhi )

            alnsg               = log(sigmag)
            aero%voltonumblo(m) = 1 / ( (pi/6.0_r8)*(dgnumlo**3)*exp(4.5_r8*alnsg**2) )
            aero%voltonumbhi(m) = 1 / ( (pi/6.0_r8)*(dgnumhi**3)*exp(4.5_r8*alnsg**2) )

            ! save sigmag of aitken mode
            if (m == aero%mode_aitken_idx) aero%sigmag_aitken = sigmag

            ! Properties of modal species
            do l = 1, aero%nspec(m)
               call rad_cnst_get_aer_props(0, m, l, &
                  density_aer = aero%specdens(l,m), &
                  hygro_aer   = aero%spechygro(l,m))
            end do

         end do

      else if (nbulk > 0) then

         aero%scheme = 'bulk'

         ! Props needed for BAM number concentration calcs.
         allocate( &
            aername(nbulk),                   &
            aero%num_to_mass_aer(nbulk),      &
            aero%mmr_bulk(nbulk),             &
            aero%mmrg_bulk(pcols,plev,nbulk)  )

         do iaer = 1, aero%nbulk
            call rad_cnst_get_aer_props(0, iaer, &
               aername         = aername(iaer), &
               num_to_mass_aer = aero%num_to_mass_aer(iaer) )
            ! Look for sulfate aerosol in this list (Bulk aerosol only)
            if (trim(aername(iaer)) == 'SULFATE') aero%idxsul   = iaer
            if (trim(aername(iaer)) == 'DUST1')   aero%idxdst1  = iaer
            if (trim(aername(iaer)) == 'DUST2')   aero%idxdst2  = iaer
            if (trim(aername(iaer)) == 'DUST3')   aero%idxdst3  = iaer
            if (trim(aername(iaer)) == 'DUST4')   aero%idxdst4  = iaer
            if (trim(aername(iaer)) == 'BCPHI')   aero%idxbcphi = iaer
         end do

      end if

   end subroutine zm_aero_init
   !----------------------------------------------------------------------------

end subroutine zm_conv_init

!===================================================================================================

subroutine zm_conv_tend(pblh, mcon, cme, tpert, dlftot, pflx, zdu, &
                        rliq, rice, ztodt, jctop, jcbot, &
                        state, ptend_all, landfrac, pbuf, mu, eu, &
                        du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath )
   !----------------------------------------------------------------------------
   use physics_types,      only: physics_state, physics_ptend
   use physics_types,      only: physics_ptend_init
   use physics_update_mod, only: physics_update
   use physics_types,      only: physics_state_copy, physics_state_dealloc
   use physics_types,      only: physics_ptend_sum, physics_ptend_dealloc
   use phys_grid,          only: get_lat_p, get_lon_p
   use time_manager,       only: get_nstep, is_first_step
   use physics_buffer,     only: pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,       only: pcnst, cnst_get_ind, cnst_is_convtran1
   use physconst,          only: gravit
   use time_manager,       only: get_curr_date
   use interpolate_data,   only: vertinterp
   !----------------------------------------------------------------------------
   ! Arguments
   type(physics_state),target,       intent(in)  :: state      ! Physics state variables
   type(physics_ptend),              intent(out) :: ptend_all  ! individual parameterization tendencies
   type(physics_buffer_desc),        pointer     :: pbuf(:)    ! physics buffer
   real(r8),                         intent(in)  :: ztodt      ! 2 delta t (model time increment)
   real(r8), dimension(pcols),       intent(in)  :: pblh       ! Planetary boundary layer height
   real(r8), dimension(pcols),       intent(in)  :: tpert      ! Thermal temperature excess
   real(r8), dimension(pcols),       intent(in)  :: landfrac   ! Land fraction
   real(r8), dimension(pcols,pverp), intent(out) :: mcon       ! Convective mass flux--m sub c
   real(r8), dimension(pcols,pver ), intent(out) :: dlftot     ! scattrd version of the detraining cld h2o tend
   real(r8), dimension(pcols,pverp), intent(out) :: pflx       ! scattered precip flux at each level
   real(r8), dimension(pcols,pver ), intent(out) :: cme        ! cmf condensation - evaporation
   real(r8), dimension(pcols,pver ), intent(out) :: zdu        ! detraining mass flux
   real(r8), dimension(pcols),       intent(out) :: rliq       ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), dimension(pcols),       intent(out) :: rice       ! reserved ice (not yet in cldice) for energy integrals
   real(r8), dimension(pcols,pver ), intent(out) :: mu         ! upward cloud mass flux
   real(r8), dimension(pcols,pver ), intent(out) :: eu         ! entrainment in updraft
   real(r8), dimension(pcols,pver ), intent(out) :: du         ! detrainment in updraft
   real(r8), dimension(pcols,pver ), intent(out) :: md         ! downward cloud mass flux
   real(r8), dimension(pcols,pver ), intent(out) :: ed         ! entrainment in downdraft
   real(r8), dimension(pcols,pver ), intent(out) :: dp         ! layer thickness
   real(r8), dimension(pcols),       intent(out) :: dsubcld    ! wg layer thickness in mbs (between upper/lower interface)
   integer,  dimension(pcols),       intent(out) :: jt         ! wg layer thickness in mbs between lcl and maxi
   integer,  dimension(pcols),       intent(out) :: maxg       ! wg top  level index of deep cumulus convection
   integer,  dimension(pcols),       intent(out) :: ideep      ! wg gathered values of maxi
   integer,                          intent(out) :: lengath    ! gathered points vs longitude index
   !----------------------------------------------------------------------------
   ! Local variables
   type(zm_microp_st) :: microp_st     ! ZM microphysics data structure

   integer :: i,k,l,m                  ! loop iterators
   integer :: nstep                    ! model time step number
   integer :: ixcldice, ixcldliq       ! constituent indices for cloud liquid and ice water
   integer :: lchnk                    ! chunk identifier
   integer :: ncol                     ! number of atmospheric columns
   integer :: itim_old                 ! for physics buffer fields
   logical :: l_windt(2)               ! flags for ptend initialization
   logical :: lq(pcnst)                ! flags for ptend initialization

   real(r8), dimension(pcols,pver) :: ftem            ! Temporary workspace for outfld variables
   real(r8), dimension(pcols,pver) :: ntprprd         ! evap outfld: net precip production in layer
   real(r8), dimension(pcols,pver) :: ntsnprd         ! evap outfld: net snow production in layer
   real(r8), dimension(pcols,pver) :: tend_s_snwprd   ! Heating rate of snow production
   real(r8), dimension(pcols,pver) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
   real(r8), dimension(pcols,pver) :: fake_dpdry      ! used in convtran call

   ! physics types
   type(physics_state)        :: state1               ! copy of state for evaporation
   type(physics_ptend),target :: ptend_loc            ! output tendencies

   ! physics buffer fields
   real(r8), pointer, dimension(:)     :: prec        ! total precipitation
   real(r8), pointer, dimension(:)     :: snow        ! snow from ZM convection 
   real(r8), pointer, dimension(:,:)   :: cld         ! cloud fraction
   real(r8), pointer, dimension(:,:)   :: ql          ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:)   :: rprd        ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis      ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:)   :: evapcdp     ! evaporation of precipitation
   real(r8), pointer, dimension(:,:)   :: flxprec     ! convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:)   :: flxsnow     ! convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:)   :: dp_cldliq   ! cloud liq water
   real(r8), pointer, dimension(:,:)   :: dp_cldice   ! cloud ice water
   real(r8), pointer, dimension(:,:)   :: wuc         ! vertical velocity

   ! DCAPE-ULL
   real(r8), pointer, dimension(:,:) :: t_star        ! DCAPE T from time step n-1
   real(r8), pointer, dimension(:,:) :: q_star        ! DCAPE q from time step n-1
   real(r8), dimension(pcols)        :: dcape         ! DCAPE cape change
   real(r8), dimension(pcols)        :: maxgsav       ! DCAPE tmp array for recording and outfld to MAXI

   real(r8), pointer, dimension(:,:) :: dlf           ! detrained convective cloud water mixing ratio
   real(r8), pointer, dimension(:,:) :: dif           ! detrained convective cloud ice mixing ratio
   real(r8), pointer, dimension(:,:) :: dsf           ! detrained convective snow mixing ratio
   real(r8), pointer, dimension(:,:) :: dnlf          ! detrained convective cloud water num concen
   real(r8), pointer, dimension(:,:) :: dnif          ! detrained convective cloud ice num concen
   real(r8), pointer, dimension(:,:) :: dnsf          ! detrained convective snow num concen
   real(r8), pointer, dimension(:,:) :: lambdadpcu    ! slope of cloud liquid size distr
   real(r8), pointer, dimension(:,:) :: mudpcu        ! width parameter of droplet size distr
   real(r8), pointer, dimension(:,:) :: qi            ! grid slice of cloud ice

   integer , dimension(pcols) :: jctop                ! row of top-of-deep-convection indices passed out
   integer , dimension(pcols) :: jcbot                ! row of base of cloud indices passed out
   real(r8), dimension(pcols) :: pcont                ! convection top pressure
   real(r8), dimension(pcols) :: pconb                ! convection base pressure
   real(r8), dimension(pcols) :: freqzm               ! fractional occurrence of ZM convection

   ! history output fields
   real(r8), dimension(pcols)      :: cape            ! convective available potential energy
   real(r8), dimension(pcols,pver) :: mu_out          ! updraft mass flux for output
   real(r8), dimension(pcols,pver) :: md_out          ! downdraft mass flux for output

   ! used in momentum transport calculations
   real(r8), dimension(pcols,pver,2) :: winds
   real(r8), dimension(pcols,pver,2) :: wind_tends
   real(r8), dimension(pcols,pver,2) :: pguall
   real(r8), dimension(pcols,pver,2) :: pgdall
   real(r8), dimension(pcols,pver,2) :: icwu
   real(r8), dimension(pcols,pver,2) :: icwd
   real(r8), dimension(pcols,pver)   :: seten

   real(r8), dimension(pcols,pver) :: sprd
   real(r8), dimension(pcols,pver) :: frz
   real(r8), dimension(pcols)      :: precz_snum

   ! MCSP
   logical  :: doslop
   logical  :: doslop_heat
   logical  :: doslop_moisture
   logical  :: doslop_uwind
   logical  :: doslop_vwind
   real(r8) :: alpha2, alpha_moisture, alphau, alphav
   real(r8) :: mcsp_top, mcsp_bot
   real(r8) :: dpg
   real(r8) :: Qcq_adjust
   real(r8), dimension(pcols)      :: Q_dis
   real(r8), dimension(pcols)      :: Qq_dis
   real(r8), dimension(pcols,pver) :: Qm
   real(r8), dimension(pcols,pver) :: Qmq
   real(r8), dimension(pcols,pver) :: Qmu
   real(r8), dimension(pcols,pver) :: Qmv
   real(r8), dimension(pcols)      :: Qm_int_end
   real(r8), dimension(pcols)      :: Qmq_int_end
   real(r8), dimension(pcols)      :: Pa_int_end
   real(r8), dimension(pcols)      :: Qs_zmconv
   real(r8), dimension(pcols)      :: Qv_zmconv
   real(r8), dimension(pcols)      :: MCSP_freq
   real(r8), dimension(pcols,pver) :: MCSP_DT
   real(r8), dimension(pcols)      :: ZM_depth
   real(r8), dimension(pcols)      :: MCSP_shear
   real(r8), dimension(pcols)      :: du600
   real(r8), dimension(pcols)      :: dv600

   !----------------------------------------------------------------------------

   if (zm_microp) then
      allocate( &
         microp_st%wu(pcols,pver),       & ! vertical velocity
         microp_st%qliq(pcols,pver),     & ! convective cloud liquid water.
         microp_st%qice(pcols,pver),     & ! convective cloud ice.
         microp_st%qrain(pcols,pver),    & ! convective rain water.
         microp_st%qsnow(pcols,pver),    & ! convective snow.
         microp_st%qgraupel(pcols,pver), & ! convective graupel
         microp_st%qnl(pcols,pver),      & ! convective cloud liquid water num concen.
         microp_st%qni(pcols,pver),      & ! convective cloud ice num concen.
         microp_st%qnr(pcols,pver),      & ! convective rain water num concen.
         microp_st%qns(pcols,pver),      & ! convective snow num concen.
         microp_st%qng(pcols,pver),      & ! convective graupel num concen.
         microp_st%autolm(pcols,pver),   & ! mass tendency due to autoconversion of droplets to rain
         microp_st%accrlm(pcols,pver),   & ! mass tendency due to accretion of droplets by rain
         microp_st%bergnm(pcols,pver),   & ! mass tendency due to Bergeron process
         microp_st%fhtimm(pcols,pver),   & ! mass tendency due to immersion freezing
         microp_st%fhtctm(pcols,pver),   & ! mass tendency due to contact freezing
         microp_st%fhmlm (pcols,pver),   & ! mass tendency due to homogeneous freezing
         microp_st%hmpim (pcols,pver),   & ! mass tendency due to HM process
         microp_st%accslm(pcols,pver),   & ! mass tendency due to accretion of droplets by snow
         microp_st%dlfm  (pcols,pver),   & ! mass tendency due to detrainment of droplet
         microp_st%autoln(pcols,pver),   & ! num tendency due to autoconversion of droplets to rain
         microp_st%accrln(pcols,pver),   & ! num tendency due to accretion of droplets by rain
         microp_st%bergnn(pcols,pver),   & ! num tendency due to Bergeron process
         microp_st%fhtimn(pcols,pver),   & ! num tendency due to immersion freezing
         microp_st%fhtctn(pcols,pver),   & ! num tendency due to contact freezing
         microp_st%fhmln (pcols,pver),   & ! num tendency due to homogeneous freezing
         microp_st%accsln(pcols,pver),   & ! num tendency due to accretion of droplets by snow
         microp_st%activn(pcols,pver),   & ! num tendency due to droplets activation
         microp_st%dlfn  (pcols,pver),   & ! num tendency due to detrainment of droplet
         microp_st%autoim(pcols,pver),   & ! mass tendency due to autoconversion of cloud ice to snow
         microp_st%accsim(pcols,pver),   & ! mass tendency due to accretion of cloud ice by snow
         microp_st%difm  (pcols,pver),   & ! mass tendency due to detrainment of cloud ice
         microp_st%nuclin(pcols,pver),   & ! num tendency due to ice nucleation
         microp_st%autoin(pcols,pver),   & ! num tendency due to autoconversion of cloud ice to snow
         microp_st%accsin(pcols,pver),   & ! num tendency due to accretion of cloud ice by snow
         microp_st%hmpin (pcols,pver),   & ! num tendency due to HM process
         microp_st%difn  (pcols,pver),   & ! num tendency due to detrainment of cloud ice
         microp_st%cmel  (pcols,pver),   & ! mass tendency due to condensation
         microp_st%cmei  (pcols,pver),   & ! mass tendency due to deposition
         microp_st%trspcm(pcols,pver),   & ! LWC tendency due to convective transport
         microp_st%trspcn(pcols,pver),   & ! droplet num tendency due to convective transport
         microp_st%trspim(pcols,pver),   & ! IWC tendency due to convective transport
         microp_st%trspin(pcols,pver),   & ! ice crystal num tendency due to convective transport
         microp_st%accgrm(pcols,pver),   & ! mass tendency due to collection of rain by graupel
         microp_st%accglm(pcols,pver),   & ! mass tendency due to collection of droplets by graupel
         microp_st%accgslm(pcols,pver),  & ! mass tendency of graupel due to collection of droplets by snow
         microp_st%accgsrm(pcols,pver),  & ! mass tendency of graupel due to collection of rain by snow
         microp_st%accgirm(pcols,pver),  & ! mass tendency of graupel due to collection of rain by ice
         microp_st%accgrim(pcols,pver),  & ! mass tendency of graupel due to collection of ice by rain
         microp_st%accgrsm(pcols,pver),  & ! mass tendency due to collection of snow by rain
         microp_st%accgsln(pcols,pver),  & ! num tendency of graupel due to collection of droplets by snow
         microp_st%accgsrn(pcols,pver),  & ! num tendency of graupel due to collection of rain by snow
         microp_st%accgirn(pcols,pver),  & ! num tendency of graupel due to collection of rain by ice   
         microp_st%accsrim(pcols,pver),  & ! mass tendency of snow due to collection of ice by rain
         microp_st%acciglm(pcols,pver),  & ! mass tendency of ice mult(splintering) due to acc droplets by graupel
         microp_st%accigrm(pcols,pver),  & ! mass tendency of ice mult(splintering) due to acc rain by graupel
         microp_st%accsirm(pcols,pver),  & ! mass tendency of snow due to collection of rain by ice
         microp_st%accigln(pcols,pver),  & ! num tendency of ice mult(splintering) due to acc droplets by graupel
         microp_st%accigrn(pcols,pver),  & ! num tendency of ice mult(splintering) due to acc rain by graupel
         microp_st%accsirn(pcols,pver),  & ! num tendency of snow due to collection of rain by ice
         microp_st%accgln(pcols,pver),   & ! num tendency due to collection of droplets by graupel
         microp_st%accgrn(pcols,pver),   & ! num tendency due to collection of rain by graupel
         microp_st%accilm(pcols,pver),   & ! mass tendency of cloud ice due to collection of droplet by cloud ice
         microp_st%acciln(pcols,pver),   & ! number conc tendency of cloud ice due to collection of droplet by cloud ice
         microp_st%fallrm(pcols,pver),   & ! mass tendency of rain fallout
         microp_st%fallsm(pcols,pver),   & ! mass tendency of snow fallout
         microp_st%fallgm(pcols,pver),   & ! mass tendency of graupel fallout
         microp_st%fallrn(pcols,pver),   & ! num tendency of rain fallout
         microp_st%fallsn(pcols,pver),   & ! num tendency of snow fallout
         microp_st%fallgn(pcols,pver),   & ! num tendency of graupel fallout
         microp_st%fhmrm (pcols,pver)    ) ! mass tendency due to homogeneous freezing of rain
   end if

   doslop          = .false.
   doslop_heat     = .false.
   doslop_moisture = .false.
   doslop_uwind    = .false.
   doslop_vwind    = .false.
   alphau          = 0
   alphav          = 0
   if ( MCSP ) then
      if ( MCSP_heat_coeff > 0 ) then
         doslop_heat = .true.
         alpha2 = MCSP_heat_coeff
      end if
      if ( MCSP_moisture_coeff > 0 ) then
         doslop_moisture = .true.
         alpha_moisture = MCSP_moisture_coeff
      end if
      if ( MCSP_uwind_coeff > 0 ) then
         doslop_uwind = .true.
         alphau = MCSP_uwind_coeff
      end if
      if ( MCSP_vwind_coeff > 0 ) then
         doslop_vwind = .true.
         alphav = MCSP_vwind_coeff
      end if
   end if

   if (doslop_heat .or. doslop_moisture .or. doslop_uwind .or. doslop_vwind) then
      doslop = .true.
   end if

   !----------------------------------------------------------------------------
   ! Initialize stuff

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem      (1:ncol,1:pver)   = 0
   mu_out    (1:ncol,1:pver)   = 0
   md_out    (1:ncol,1:pver)   = 0
   dlftot    (1:ncol,1:pver)   = 0
   wind_tends(1:ncol,1:pver,:) = 0

   call physics_state_copy(state,state1) ! make local copy of state

   lq(:) = .FALSE.
   lq(1) = .TRUE.
 
   call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq, lu=doslop_uwind, lv=doslop_vwind)

   !----------------------------------------------------------------------------
   ! Associate pointers with physics buffer fields

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,           cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, fracis_idx,        fracis, start=(/1,1,1/),        kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, icwmrdp_idx,       ql)
   call pbuf_get_field(pbuf, rprddp_idx,        rprd)
   call pbuf_get_field(pbuf, nevapr_dpcu_idx,   evapcdp)
   call pbuf_get_field(pbuf, prec_dp_idx,       prec)
   call pbuf_get_field(pbuf, snow_dp_idx,       snow)
   call pbuf_get_field(pbuf, icimrdp_idx,       qi)
   call pbuf_get_field(pbuf, dlfzm_idx,         dlf)
   call pbuf_get_field(pbuf, difzm_idx,         dif)

   ! DCAPE-ULL
   if (trigdcape_ull .or. trig_dcape_only) then
      call pbuf_get_field(pbuf, t_star_idx, t_star)
      call pbuf_get_field(pbuf, q_star_idx, q_star)
      if ( is_first_step()) then
         q_star(1:ncol,1:pver) = state%q(1:ncol,1:pver,1)
         t_star(1:ncol,1:pver) = state%t(1:ncol,1:pver)
      end if
   end if

   if (zm_microp) then
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlf)
      call pbuf_get_field(pbuf, dnifzm_idx, dnif)
      call pbuf_get_field(pbuf, dsfzm_idx,  dsf)
      call pbuf_get_field(pbuf, dnsfzm_idx, dnsf)
      call pbuf_get_field(pbuf, wuc_idx,    wuc)
   else
      allocate(dnlf(pcols,pver), &
               dnif(pcols,pver), &
               dsf(pcols,pver),  &
               dnsf(pcols,pver), &
               wuc(pcols,pver)   )
   end if
   wuc(1:pcols,1:pver) = 0

   call pbuf_get_field(pbuf, lambdadpcu_idx, lambdadpcu)
   call pbuf_get_field(pbuf, mudpcu_idx,     mudpcu)

   if (zm_microp) then
      if (nmodes > 0) then

         ! Associate pointers with the modes and species that affect the climate (list 0)
         do m = 1, nmodes
            call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, aero(lchnk)%num_a(m)%val)
            call pbuf_get_field(pbuf, dgnum_idx, aero(lchnk)%dgnum(m)%val, start=(/1,1,m/), kount=(/pcols,pver,1/))
            do l = 1, aero(lchnk)%nspec(m)
               call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, aero(lchnk)%mmr_a(l,m)%val)
            end do
         end do

      else if (nbulk > 0) then

         ! Associate pointers with the bulk aerosols that affect the climate (list 0)
         do m = 1, nbulk
            call rad_cnst_get_aer_mmr(0, m, state, pbuf, aero(lchnk)%mmr_bulk(m)%val)
         end do

      end if
   end if

   ! Call the primary Zhang-McFarlane convection parameterization
   call t_startf ('zm_convr')
   call zm_convr( lchnk, ncol, state%t, state%q(:,:,1), prec, jctop, jcbot, &
                  pblh, state%zm, state%phis, state%zi, ptend_loc%q(:,:,1), &
                  ptend_loc%s, state%pmid, state%pint, state%pdel, state%omega, &
                  0.5*ztodt, mcon, cme, cape, tpert, dlf, pflx, zdu, rprd, mu, md, du, eu, ed, &
                  dp, dsubcld, jt, maxg, ideep, lengath, ql, rliq, landfrac, &
                  t_star, q_star, dcape, &  
                  aero(lchnk), qi, dif, dnlf, dnif, dsf, dnsf, sprd, rice, frz, mudpcu, &
                  lambdadpcu, microp_st, wuc )
   call t_stopf ('zm_convr')

   if (zm_microp) then
      dlftot(1:ncol,1:pver) = dlf(1:ncol,1:pver) + dif(1:ncol,1:pver) + dsf(1:ncol,1:pver)
   else
      dlftot(1:ncol,1:pver) = dlf(1:ncol,1:pver)
   end if

   !----------------------------------------------------------------------------
   ! Begin MCSP parameterization here
   !----------------------------------------------------------------------------
   if (doslop) then
      du600              = 0
      dv600              = 0
      ZM_depth           = 0
      MCSP_shear         = 0
      Qs_zmconv (1:ncol) = 0
      Qv_zmconv (1:ncol) = 0
      Pa_int_end(1:ncol) = 0
      Qm (1:ncol,1:pver) = 0
      Qmq(1:ncol,1:pver) = 0
      Qmu(1:ncol,1:pver) = 0
      Qmv(1:ncol,1:pver) = 0

      call vertinterp(ncol, pcols, pver, state%pmid, MCSP_storm_speed_pref, state%u,du600)
      call vertinterp(ncol, pcols, pver, state%pmid, MCSP_storm_speed_pref, state%v,dv600)

      do i = 1,ncol
         if (state%pmid(i,pver).gt.MCSP_storm_speed_pref) then
            du600(i) = du600(i)-state%u(i,pver)
            dv600(i) = dv600(i)-state%v(i,pver)
         else
            du600(i) = -999
            dv600(i) = -999
         end if
      end do

      MCSP_shear = du600
      do i = 1,ncol
         if ( jctop(i).ne.pver ) ZM_depth(i) = state%pint(i,pver+1) - state%pmid(i,jctop(i))
      end do

      ! Define parameters
      do i = 1,ncol
         do k = jctop(i),pver
            Qs_zmconv(i) = Qs_zmconv(i) + ptend_loc%s(i,k) * state%pdel(i,k)
            Qv_zmconv(i) = Qv_zmconv(i) + ptend_loc%q(i,k,1) * state%pdel(i,k)
            Pa_int_end(i) = Pa_int_end(i) + state%pdel(i,k)
         end do
      end do

      do i = 1,ncol
         if (jctop(i) .ne. pver) then
            Qs_zmconv(i) = Qs_zmconv(i) / Pa_int_end(i)
            Qv_zmconv(i) = Qv_zmconv(i) / Pa_int_end(i)
         else
            Qs_zmconv(i) = 0
            Qv_zmconv(i) = 0
         end if
      end do

      do i = 1,ncol
         Qm_int_end(i)  = 0
         Qmq_int_end(i) = 0
         Pa_int_end(i)  = 0
         Q_dis(i)       = 0
         Qq_dis(i)      = 0

         if ( (state%pint(i,pver+1)-state%pmid(i,jctop(i))) >= MCSP_conv_depth_min ) then
            if ( abs(du600(i)).ge.MCSP_shear_min .and. &
                 abs(du600(i)).lt.MCSP_shear_max .and. &
                 Qs_zmconv(i).gt.0 ) then
               do k = jctop(i),pver
                  mcsp_top = state%pint(i,pver+1) - state%pmid(i,k)
                  mcsp_bot = state%pint(i,pver+1) - state%pmid(i,jctop(i))

                  Qm(i,k)  = -1 * Qs_zmconv(i) * alpha2 * sin(2.0_r8*pi*(mcsp_top/mcsp_bot))
                  Qmq(i,k) = -1 * Qv_zmconv(i) * alpha2 * sin(2.0_r8*pi*(mcsp_top/mcsp_bot))
                  Qmq(i,k) = Qm(i,k)/2500000.0_r8/4.0_r8

                  Qmu(i,k) = alphau * (cos(pi*(mcsp_top/mcsp_bot)))
                  Qmv(i,k) = alphav * (cos(pi*(mcsp_top/mcsp_bot)))

                  dpg = state%pdel(i,k)/gravit

                  Qm_int_end(i)  = Qm_int_end(i) + Qm(i,k) * dpg
                  Qm_int_end(i)  = Qm_int_end(i) + (2.0_r8*Qmu(i,k)*ztodt*state%u(i,k)+ &
                                                   Qmu(i,k)*Qmu(i,k)*ztodt*ztodt)/2.0_r8 * dpg/ztodt
                  Qm_int_end(i)  = Qm_int_end(i) + (2.0_r8*Qmv(i,k)*ztodt*state%v(i,k)+ &
                                                   Qmv(i,k)*Qmv(i,k)*ztodt*ztodt)/2.0_r8 * dpg/ztodt
                  Qmq_int_end(i) = Qmq_int_end(i) + Qmq(i,k) * dpg
                  Pa_int_end(i)  = Pa_int_end(i) + state%pdel(i,k)
               end do
               Q_dis(i)  = Qm_int_end(i)  / Pa_int_end(i)
               Qq_dis(i) = Qmq_int_end(i) / Pa_int_end(i)
            end if
         end if
      end do

      MCSP_DT   = 0
      MCSP_freq = 0

      do i = 1,ncol
         do k = jctop(i),pver
            Qcq_adjust = ptend_loc%q(i,k,1) - Qq_dis(i) * gravit

            ptend_loc%s(i,k) = ptend_loc%s(i,k) - Q_dis(i)*gravit ! energy fixer

            MCSP_DT(i,k) = -Q_dis(i)*gravit+Qm(i,k)
            if (abs(Qm(i,k)).gt.0 .and. abs(Qmu(i,k)).gt.0) MCSP_freq(i) = 1

            if (doslop_heat)     ptend_loc%s(i,k)   = ptend_loc%s(i,k) + Qm(i,k)
            if (doslop_moisture) ptend_loc%q(i,k,1) = Qcq_adjust + Qmq(i,k)
            if (doslop_uwind)    ptend_loc%u(i,k)   = Qmu(i,k)
            if (doslop_vwind)    ptend_loc%v(i,k)   = Qmv(i,k)
         end do
      end do

      MCSP_DT(1:ncol,1:pver) = MCSP_DT(1:ncol,1:pver)/cpair
      call outfld('MCSP_DT    ',MCSP_DT,   pcols, lchnk )
      call outfld('MCSP_freq  ',MCSP_freq, pcols, lchnk )
      if (doslop_uwind) call outfld('MCSP_DU    ',ptend_loc%u*86400.0_r8, pcols, lchnk )
      if (doslop_vwind) call outfld('MCSP_DV    ',ptend_loc%v*86400.0_r8, pcols, lchnk )
      call outfld('ZM_depth   ',ZM_depth,   pcols, lchnk )
      call outfld('MCSP_shear ',MCSP_shear, pcols, lchnk )

   end if
   !----------------------------------------------------------------------------
   ! End of MCSP parameterization calculations
   !----------------------------------------------------------------------------

   call outfld('DCAPE',  dcape, pcols, lchnk )
   call outfld('CAPE_ZM', cape, pcols, lchnk )

   ! Output fractional occurrence of ZM convection
   freqzm(:) = 0
   do i = 1,lengath
      freqzm(ideep(i)) = 1.0_r8
   end do
   call outfld('FREQZM  ', freqzm, pcols, lchnk )

   ! Convert mass flux from reported mb/s to kg/m^2/s
   mcon(1:ncol,1:pver) = mcon(1:ncol,1:pver) * 100.0_r8/gravit

   call outfld('CMFMCDZM', mcon, pcols, lchnk )

   ! Store upward and downward mass fluxes in un-gathered arrays + convert from mb/s to kg/m^2/s
   do i = 1,lengath
      do k = 1,pver
         mu_out(ideep(i),k) = mu(i,k) * 100.0_r8/gravit
         md_out(ideep(i),k) = md(i,k) * 100.0_r8/gravit
      end do
   end do

   if (convproc_do_aer .or. convproc_do_gas) then 
      call outfld('ZMMU', mu_out, pcols, lchnk )
      call outfld('ZMMD', md_out, pcols, lchnk )
   else
      call outfld('ZMMU', mu_out(1,1), pcols, lchnk )
      call outfld('ZMMD', md_out(1,1), pcols, lchnk )
   end if

   ftem(1:ncol,1:pver) = ptend_loc%s(1:ncol,1:pver)/cpair
   call outfld('ZMDT    ', ftem, pcols, lchnk )
   call outfld('ZMDQ    ', ptend_loc%q(1,1,1), pcols, lchnk )

   if (zm_microp) call zm_conv_micro_outfld( microp_st, dlf, dif, dnlf, dnif, frz, lchnk, ncol )

   maxgsav(1:ncol) = 0 ! zero if no convection. true mean to be MAXI/FREQZM
   pcont(1:ncol) = state%ps(1:ncol)
   pconb(1:ncol) = state%ps(1:ncol)
   do i = 1,lengath
      if (maxg(i).gt.jt(i)) then
         pcont(ideep(i)) = state%pmid(ideep(i),jt(i))  ! gathered array (or jctop ungathered)
         pconb(ideep(i)) = state%pmid(ideep(i),maxg(i))! gathered array
         maxgsav(ideep(i)) = dble(maxg(i))             ! gathered array for launching level
      end if
   end do
   call outfld('PCONVT  ', pcont, pcols, lchnk )
   call outfld('PCONVB  ', pconb, pcols, lchnk )
   call outfld('MAXI  ', maxgsav, pcols, lchnk )

   call physics_ptend_init(ptend_all, state%psetcols, 'zm_conv_tend')

   ! add tendency from this process to tendencies from other processes
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   ! update physics state type state1 with ptend_loc 
   call physics_update(state1, ptend_loc, ztodt)

   ! initialize ptend for next process
   lq(:) = .FALSE.
   lq(1) = .TRUE.
   call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

   ! Determine the phase of the precipitation produced and add latent heat of fusion
   ! Evaporate some of the precip directly into the environment
   ! Allow this to use the updated copy of state (state1) and the fresh ptend_loc type
   ! heating and specific humidity tendencies produced

   call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
   call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
   call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
   call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
   dp_cldliq(1:ncol,1:pver) = 0
   dp_cldice(1:ncol,1:pver) = 0

   call t_startf ('zm_conv_evap')
   call zm_conv_evap(state1%ncol, state1%lchnk, &
                     state1%t, state1%pmid, state1%pdel, state1%q(1:pcols,1:pver,1), &
                     ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
                     rprd, cld, ztodt, prec, snow, ntprprd, ntsnprd, flxprec, flxsnow, sprd, old_snow)
   call t_stopf ('zm_conv_evap')

   evapcdp(1:ncol,1:pver) = ptend_loc%q(1:ncol,1:pver,1)

   ! Write out variables from zm_conv_evap
   ftem(1:ncol,1:pver) = ptend_loc%s(1:ncol,1:pver)/cpair
   call outfld('EVAPTZM ', ftem, pcols, lchnk )
   ftem(1:ncol,1:pver) = tend_s_snwprd  (1:ncol,1:pver)/cpair
   call outfld('FZSNTZM ', ftem, pcols, lchnk )
   ftem(1:ncol,1:pver) = tend_s_snwevmlt(1:ncol,1:pver)/cpair
   call outfld('EVSNTZM ', ftem, pcols, lchnk )

   call outfld('EVAPQZM ', ptend_loc%q(1,1,1), pcols, lchnk )
   call outfld('ZMFLXPRC', flxprec,            pcols, lchnk )
   call outfld('ZMFLXSNW', flxsnow,            pcols, lchnk )
   call outfld('ZMNTPRPD', ntprprd,            pcols, lchnk )
   call outfld('ZMNTSNPD', ntsnprd,            pcols, lchnk )
   call outfld('ZMEIHEAT', ptend_loc%s,        pcols, lchnk )
   call outfld('PRECCDZM', prec,               pcols, lchnk )
   call outfld('PRECZ   ', prec,               pcols, lchnk )

   if (zm_microp) then
      do i = 1,ncol
         if (prec(i) .gt. 0) then
            precz_snum(i) = 1
         else
            precz_snum(i) = 0
         end if
      end do
      call outfld('PRECZ_SN', precz_snum, pcols, lchnk )
   end if

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   ! update physics state type state1 with ptend_loc
   call physics_update(state1, ptend_loc, ztodt)

   ! Momentum Transport
   call physics_ptend_init(ptend_loc, state1%psetcols, 'momtran', ls=.true., lu=.true., lv=.true.)

   winds(1:ncol,1:pver,1) = state1%u(1:ncol,1:pver)
   winds(1:ncol,1:pver,2) = state1%v(1:ncol,1:pver)

   l_windt(1) = .true.
   l_windt(2) = .true.

   call t_startf ('momtran')
   call momtran(lchnk, ncol, l_windt, winds, 2, mu, md, du, eu, ed, dp, dsubcld, &
                jt, maxg, ideep, 1, lengath, nstep, wind_tends, pguall, pgdall, &
                icwu, icwd, ztodt, seten )
   call t_stopf ('momtran')

   ptend_loc%u(1:ncol,1:pver) = wind_tends(1:ncol,1:pver,1)
   ptend_loc%v(1:ncol,1:pver) = wind_tends(1:ncol,1:pver,2)
   ptend_loc%s(1:ncol,1:pver) = seten     (1:ncol,1:pver)

   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   ! update physics state type state1 with ptend_loc
   call physics_update(state1, ptend_loc, ztodt)

   ftem(1:ncol,1:pver) = seten(1:ncol,1:pver)/cpair
   call outfld('ZMMTT', ftem             , pcols, lchnk )
   call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk )
   call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk )

   ! Output apparent force from  pressure gradient
   call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk )
   call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk )
   call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk )
   call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk )

   ! Output in-cloud winds
   call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk )
   call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk )
   call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk )
   call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk )

   ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)
   call physics_ptend_init(ptend_loc, state1%psetcols, 'convtran1', lq=lq)

   ! dpdry is not used in next convtran call since cloud liq/ice mixing ratios are moist
   fake_dpdry(1:ncol,1:pver) = 0

   call t_startf ('convtran1')
   call convtran( lchnk, ptend_loc%lq, state1%q, pcnst, mu, md, &
                  du, eu, ed, dp, dsubcld, jt, maxg, ideep, 1, lengath, &
                  nstep, fracis, ptend_loc%q, fake_dpdry, ztodt)  
   call t_stopf ('convtran1')

   call outfld('ZMDICE ', ptend_loc%q(1,1,ixcldice), pcols, lchnk )
   call outfld('ZMDLIQ ', ptend_loc%q(1,1,ixcldliq), pcols, lchnk )

   ! add tendency from this process to tendency from other processes
   call physics_ptend_sum( ptend_loc, ptend_all, ncol )

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)

   if (zm_microp) then
      deallocate( &
         microp_st%wu,       &
         microp_st%qliq,     &
         microp_st%qice,     &
         microp_st%qrain,    &
         microp_st%qsnow,    &
         microp_st%qgraupel, &
         microp_st%qnl,      &
         microp_st%qni,      &
         microp_st%qnr,      &
         microp_st%qns,      &
         microp_st%qng,      &
         microp_st%autolm,   &
         microp_st%accrlm,   &
         microp_st%bergnm,   &
         microp_st%fhtimm,   &
         microp_st%fhtctm,   &
         microp_st%fhmlm ,   &
         microp_st%hmpim ,   &
         microp_st%accslm,   &
         microp_st%dlfm  ,   &
         microp_st%autoln,   &
         microp_st%accrln,   &
         microp_st%bergnn,   &
         microp_st%fhtimn,   &
         microp_st%fhtctn,   &
         microp_st%fhmln ,   &
         microp_st%accsln,   &
         microp_st%activn,   &
         microp_st%dlfn  ,   &
         microp_st%autoim,   &
         microp_st%accsim,   &
         microp_st%difm  ,   &
         microp_st%nuclin,   &
         microp_st%autoin,   &
         microp_st%accsin,   &
         microp_st%hmpin ,   &
         microp_st%difn  ,   &
         microp_st%cmel  ,   &
         microp_st%cmei  ,   &
         microp_st%trspcm,   &
         microp_st%trspcn,   &
         microp_st%trspim,   &
         microp_st%trspin,   &
         microp_st%accgrm,   &
         microp_st%accglm,   &
         microp_st%accgslm,  &
         microp_st%accgsrm,  &
         microp_st%accgirm,  &
         microp_st%accgrim,  &
         microp_st%accgrsm,  &
         microp_st%accgsln,  &
         microp_st%accgsrn,  &
         microp_st%accgirn,  &
         microp_st%accsrim,  &
         microp_st%acciglm,  &
         microp_st%accigrm,  &
         microp_st%accsirm,  &
         microp_st%accigln,  &
         microp_st%accigrn,  &
         microp_st%accsirn,  &
         microp_st%accgln,   &
         microp_st%accgrn,   &
         microp_st%accilm,   &
         microp_st%acciln,   &
         microp_st%fallrm,   &
         microp_st%fallsm,   &
         microp_st%fallgm,   &
         microp_st%fallrn,   &
         microp_st%fallsn,   &
         microp_st%fallgn,   &
         microp_st%fhmrm     )
   else
      deallocate(dnlf, dnif, dsf, dnsf)
   end if

end subroutine zm_conv_tend

!===================================================================================================

subroutine zm_conv_tend_2( state,  ptend,  ztodt, pbuf, mu, eu, du, md, ed, dp, &
                           dsubcld, jt, maxg, ideep, lengath, species_class)
   !----------------------------------------------------------------------------
   use physics_types,      only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,       only: get_nstep
   use physics_buffer,     only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,       only: pcnst, cnst_get_ind, cnst_is_convtran1
   use error_messages,     only: alloc_err
   use physconst,          only: spec_class_aerosol, spec_class_gas
   !----------------------------------------------------------------------------
   ! Arguments
   type(physics_state),             intent(in ):: state        ! Physics state variables
   type(physics_ptend),             intent(out):: ptend        ! indivdual parameterization tendencies
   type(physics_buffer_desc),          pointer :: pbuf(:)      ! physics buffer
   real(r8),                        intent(in) :: ztodt        ! 2 delta t (model time increment)
   real(r8), dimension(pcols,pver), intent(in) :: mu           ! upward cloud mass flux
   real(r8), dimension(pcols,pver), intent(in) :: eu           ! entrainment in updraft
   real(r8), dimension(pcols,pver), intent(in) :: du           ! detrainment in updraft
   real(r8), dimension(pcols,pver), intent(in) :: md           ! downward cloud mass flux
   real(r8), dimension(pcols,pver), intent(in) :: ed           ! entrainment in downdraft
   real(r8), dimension(pcols,pver), intent(in) :: dp           ! layer thickness
   real(r8), dimension(pcols),      intent(in) :: dsubcld      ! wg layer thickness in mbs (between interfaces)
   integer,  dimension(pcols),      intent(in) :: jt           ! wg layer thickness in mbs between lcl and maxi
   integer,  dimension(pcols),      intent(in) :: maxg         ! wg top  level index of deep cumulus convection
   integer,  dimension(pcols),      intent(in) :: ideep        ! wg gathered values of maxi
   integer,                         intent(in) :: lengath      ! number of gathered columns
   integer,  dimension(:),          intent(in) :: species_class! constituent tracer type
   !----------------------------------------------------------------------------
   ! Local variables
   integer :: i, lchnk, ncol, istat, m
   integer :: nstep
   logical :: lq(pcnst)
   integer :: ifld
   real(r8), dimension(pcols,pver) :: dpdry
   real(r8), pointer, dimension(:,:,:) :: fracis ! fraction of transported species that are insoluble
   !----------------------------------------------------------------------------
   ! Initialize

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   lq(:) = .FALSE.
   lq(:) = .not. cnst_is_convtran1(:)
   call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=lq )

   ! Associate pointers with physics buffer fields
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

   ! Transport all constituents except cloud water and ice
   if ( (convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero ) then
      do m = 1, pcnst
         if ( (species_class(m)==spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m)==spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
         end if
      end do
   end if

   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran for tracers of dry mixing ratio type
      dpdry(1:ncol,1:pver) = 0
      do i = 1,lengath
         dpdry(i,1:pver) = state%pdeldry(ideep(i),1:pver)/100_r8
      end do
      call t_startf ('convtran2')
      call convtran (lchnk, ptend%lq, state%q, pcnst,  mu, md, &
                     du, eu, ed, dp, dsubcld,  jt, maxg, ideep, 1, lengath, &
                     nstep, fracis, ptend%q, dpdry, ztodt)
      call t_stopf ('convtran2')
   end if

   if ( (convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero ) then
      do m = 1, pcnst
         if ( (species_class(m)==spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m)==spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
            ptend%q(1:ncol,1:pver,m) = 0
         end if
      end do
   end if

end subroutine zm_conv_tend_2

!===================================================================================================

subroutine zm_conv_micro_outfld(microp_st, dlf, dif, dnlf, dnif, frz, lchnk, ncol)
   !----------------------------------------------------------------------------
   ! Arguments
   type(zm_microp_st),intent(in) :: microp_st   ! ZM microphysics data structure
   real(r8),          intent(in) :: dlf(:,:)    ! detrainment of conv cld liq water mixing ratio
   real(r8),          intent(in) :: dif(:,:)    ! detrainment of conv cld ice mixing ratio
   real(r8),          intent(in) :: dnlf(:,:)   ! detrainment of conv cld liq water num concen
   real(r8),          intent(in) :: dnif(:,:)   ! detrainment of conv cld ice num concen
   real(r8),          intent(in) :: frz(:,:)    ! heating rate due to freezing 
   integer,           intent(in) :: lchnk       ! chunk identifier
   integer,           intent(in) :: ncol        ! number of columns in chunk
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k
   real(r8) :: cice_snum(pcols,pver)            ! convective cloud ice sample number
   real(r8) :: cliq_snum(pcols,pver)            ! convective cloud liquid sample number
   real(r8) :: crain_snum(pcols,pver)           ! convective rain water sample number
   real(r8) :: csnow_snum(pcols,pver)           ! convective snow sample number
   real(r8) :: cgraupel_snum(pcols,pver)        ! convective graupel sample number
   real(r8) :: wu_snum(pcols,pver)              ! vertical velocity sample number
   !----------------------------------------------------------------------------
   do k = 1,pver
      do i = 1,ncol
         if (microp_st%qice(i,k)     >  0) cice_snum(i,k)     = 1
         if (microp_st%qice(i,k)     <= 0) cice_snum(i,k)     = 0
         if (microp_st%qliq(i,k)     >  0) cliq_snum(i,k)     = 1
         if (microp_st%qliq(i,k)     <= 0) cliq_snum(i,k)     = 0
         if (microp_st%qsnow(i,k)    >  0) csnow_snum(i,k)    = 1
         if (microp_st%qsnow(i,k)    <= 0) csnow_snum(i,k)    = 0
         if (microp_st%qrain(i,k)    >  0) crain_snum(i,k)    = 1
         if (microp_st%qrain(i,k)    <= 0) crain_snum(i,k)    = 0
         if (microp_st%qgraupel(i,k) >  0) cgraupel_snum(i,k) = 1
         if (microp_st%qgraupel(i,k) <= 0) cgraupel_snum(i,k) = 0
         if (microp_st%wu(i,k)       >  0) wu_snum(i,k)       = 1
         if (microp_st%wu(i,k)       <= 0) wu_snum(i,k)       = 0
      end do
   end do

   call outfld('CLIQSNUM',cliq_snum          , pcols, lchnk )
   call outfld('CICESNUM',cice_snum          , pcols, lchnk )
   call outfld('CRAINNUM',crain_snum         , pcols, lchnk )
   call outfld('CSNOWNUM',csnow_snum         , pcols, lchnk )
   call outfld('CGRAPNUM',cgraupel_snum      , pcols, lchnk )
   call outfld('WUZMSNUM',wu_snum            , pcols, lchnk )

   call outfld('DIFZM'   ,dif                , pcols, lchnk )
   call outfld('DLFZM'   ,dlf                , pcols, lchnk )
   call outfld('DNIFZM'  ,dnif               , pcols, lchnk )
   call outfld('DNLFZM'  ,dnlf               , pcols, lchnk )
   call outfld('FRZZM'   ,frz                , pcols, lchnk )

   call outfld('WUZM'    ,microp_st%wu       , pcols, lchnk )

   call outfld('CLDLIQZM',microp_st%qliq     , pcols, lchnk )
   call outfld('CLDICEZM',microp_st%qice     , pcols, lchnk )
   call outfld('QRAINZM' ,microp_st%qrain    , pcols, lchnk )
   call outfld('QSNOWZM' ,microp_st%qsnow    , pcols, lchnk )
   call outfld('QGRAPZM' ,microp_st%qgraupel , pcols, lchnk )

   call outfld('QNLZM'   ,microp_st%qnl      , pcols, lchnk )
   call outfld('QNIZM'   ,microp_st%qni      , pcols, lchnk )
   call outfld('QNRZM'   ,microp_st%qnr      , pcols, lchnk )
   call outfld('QNSZM'   ,microp_st%qns      , pcols, lchnk )
   call outfld('QNGZM'   ,microp_st%qng      , pcols, lchnk )

   call outfld('AUTOL_M' ,microp_st%autolm   , pcols, lchnk )
   call outfld('ACCRL_M' ,microp_st%accrlm   , pcols, lchnk )
   call outfld('BERGN_M' ,microp_st%bergnm   , pcols, lchnk )
   call outfld('FHTIM_M' ,microp_st%fhtimm   , pcols, lchnk )
   call outfld('FHTCT_M' ,microp_st%fhtctm   , pcols, lchnk )
   call outfld('FHML_M'  ,microp_st%fhmlm    , pcols, lchnk )
   call outfld('HMPI_M'  ,microp_st%hmpim    , pcols, lchnk )
   call outfld('ACCSL_M' ,microp_st%accslm   , pcols, lchnk )
   call outfld('DLF_M'   ,microp_st%dlfm     , pcols, lchnk )

   call outfld('AUTOL_N' ,microp_st%autoln   , pcols, lchnk )
   call outfld('ACCRL_N' ,microp_st%accrln   , pcols, lchnk )
   call outfld('BERGN_N' ,microp_st%bergnn   , pcols, lchnk )
   call outfld('FHTIM_N' ,microp_st%fhtimn   , pcols, lchnk )
   call outfld('FHTCT_N' ,microp_st%fhtctn   , pcols, lchnk )
   call outfld('FHML_N'  ,microp_st%fhmln    , pcols, lchnk )
   call outfld('ACCSL_N' ,microp_st%accsln   , pcols, lchnk )
   call outfld('ACTIV_N' ,microp_st%activn   , pcols, lchnk )
   call outfld('DLF_N'   ,microp_st%dlfn     , pcols, lchnk )
   call outfld('AUTOI_M' ,microp_st%autoim   , pcols, lchnk )
   call outfld('ACCSI_M' ,microp_st%accsim   , pcols, lchnk )
   call outfld('DIF_M'   ,microp_st%difm     , pcols, lchnk )
   call outfld('NUCLI_N' ,microp_st%nuclin   , pcols, lchnk )
   call outfld('AUTOI_N' ,microp_st%autoin   , pcols, lchnk )
   call outfld('ACCSI_N' ,microp_st%accsin   , pcols, lchnk )
   call outfld('HMPI_N'  ,microp_st%hmpin    , pcols, lchnk )
   call outfld('DIF_N'   ,microp_st%difn     , pcols, lchnk )
   call outfld('COND_M'  ,microp_st%cmel     , pcols, lchnk )
   call outfld('DEPOS_M' ,microp_st%cmei     , pcols, lchnk )

   call outfld('TRSPC_M' ,microp_st%trspcm   , pcols, lchnk )
   call outfld('TRSPC_N' ,microp_st%trspcn   , pcols, lchnk )
   call outfld('TRSPI_M' ,microp_st%trspim   , pcols, lchnk )
   call outfld('TRSPI_N' ,microp_st%trspin   , pcols, lchnk )

   call outfld('ACCGR_M' ,microp_st%accgrm   , pcols, lchnk )
   call outfld('ACCGL_M' ,microp_st%accglm   , pcols, lchnk )
   call outfld('ACCGSL_M',microp_st%accgslm  , pcols, lchnk )
   call outfld('ACCGSR_M',microp_st%accgsrm  , pcols, lchnk )
   call outfld('ACCGIR_M',microp_st%accgirm  , pcols, lchnk )
   call outfld('ACCGRI_M',microp_st%accgrim  , pcols, lchnk )
   call outfld('ACCGRS_M',microp_st%accgrsm  , pcols, lchnk )

   call outfld('ACCGSL_N',microp_st%accgsln  , pcols, lchnk )
   call outfld('ACCGSR_N',microp_st%accgsrn  , pcols, lchnk )
   call outfld('ACCGIR_N',microp_st%accgirn  , pcols, lchnk )

   call outfld('ACCSRI_M',microp_st%accsrim  , pcols, lchnk )
   call outfld('ACCIGL_M',microp_st%acciglm  , pcols, lchnk )
   call outfld('ACCIGR_M',microp_st%accigrm  , pcols, lchnk )
   call outfld('ACCSIR_M',microp_st%accsirm  , pcols, lchnk )

   call outfld('ACCIGL_N',microp_st%accigln  , pcols, lchnk )
   call outfld('ACCIGR_N',microp_st%accigrn  , pcols, lchnk )
   call outfld('ACCSIR_N',microp_st%accsirn  , pcols, lchnk )
   call outfld('ACCGL_N' ,microp_st%accgln   , pcols, lchnk )
   call outfld('ACCGR_N' ,microp_st%accgrn   , pcols, lchnk )

   call outfld('ACCIL_M' ,microp_st%accilm   , pcols, lchnk )
   call outfld('ACCIL_N' ,microp_st%acciln   , pcols, lchnk )

   call outfld('FALLR_M' ,microp_st%fallrm   , pcols, lchnk )
   call outfld('FALLS_M' ,microp_st%fallsm   , pcols, lchnk )
   call outfld('FALLG_M' ,microp_st%fallgm   , pcols, lchnk )
   call outfld('FALLR_N' ,microp_st%fallrn   , pcols, lchnk )
   call outfld('FALLS_N' ,microp_st%fallsn   , pcols, lchnk )
   call outfld('FALLG_N' ,microp_st%fallgn   , pcols, lchnk )

   call outfld('FHMR_M'  ,microp_st%fhmrm    , pcols, lchnk )

end subroutine zm_conv_micro_outfld

!===================================================================================================

end module zm_conv_intr
