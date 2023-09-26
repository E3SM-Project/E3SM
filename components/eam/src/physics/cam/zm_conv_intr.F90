


module zm_conv_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Zhang-McFarlane deep convection scheme
!
! Author: D.B. Coleman
! January 2010 modified by J. Kay to add COSP simulator fields to physics buffer
! July 2015 B. Singh Added code for unified convective trasport
! March 2023 X. Song added code for the initialization of convective microphysics
!            and aerosol object, and the output of microphysical properties and tendencies
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv,      only: zm_conv_evap, zm_convr, convtran, momtran, trigdcape_ull, trig_dcape_only
   use zm_conv,      only: zm_microp
   use zm_microphysics,  only: zm_aero_t, zm_microp_st
   use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
                               rad_cnst_get_aer_props, rad_cnst_get_mode_props  

   use ndrop_bam,        only: ndrop_bam_init
   use cam_abortutils,   only: endrun
   use physconst,        only: pi

   use spmd_utils,       only: masterproc
   use cam_history,  only: outfld, addfld, horiz_only, add_default
   use perf_mod
   use cam_logfile,  only: iulog 
   use zm_conv,      only: MCSP, MCSP_heat_coeff, MCSP_moisture_coeff, MCSP_uwind_coeff, MCSP_vwind_coeff

   implicit none
   private
   save

   ! Public methods

   public ::&
      zm_conv_register,           &! register fields in physics buffer
      zm_conv_init,               &! initialize donner_deep module
      zm_conv_tend,               &! return tendencies
      zm_conv_tend_2               ! return tendencies

   ! Private module data
   integer ::& ! indices for fields in the physics buffer
      dp_flxprc_idx, &
      dp_flxsnw_idx, &
      dp_cldliq_idx, &
      dp_cldice_idx, &
      dlfzm_idx,     &     ! detrained convective cloud water mixing ratio.
      difzm_idx,     &     ! detrained convective cloud ice mixing ratio.
      dsfzm_idx,     &     ! detrained convective snow mixing ratio. 
      dnlfzm_idx,    &     ! detrained convective cloud water num concen.
      dnifzm_idx,    &     ! detrained convective cloud ice num concen.
      dnsfzm_idx,    &     ! detrained convective snow num concen.
      prec_dp_idx,   &
      snow_dp_idx,   &
      wuc_idx

! DCAPE-ULL
   integer :: t_star_idx       !t_star index in physics buffer
   integer :: q_star_idx       !q_star index in physics buffer

!  indices for fields in the physics buffer
   integer  ::    cld_idx          = 0    
   integer  ::    icwmrdp_idx      = 0     
   integer  ::    rprddp_idx       = 0    
   integer  ::    fracis_idx       = 0   
   integer  ::    nevapr_dpcu_idx  = 0    
   logical  ::    convproc_do_aer 
   logical  ::    convproc_do_gas 
   logical  ::    clim_modal_aero

   integer  ::    dgnum_idx        = 0
   integer  ::    lambdadpcu_idx   = 0
   integer  ::    mudpcu_idx       = 0
   integer  ::    icimrdp_idx      = 0


   logical :: old_snow  = .true.   ! set true to use old estimate of snow production in zm_conv_evap
                                   ! set false to use snow production from zm
                                   ! microphysics
   integer :: nmodes
   integer :: nbulk

   type(zm_aero_t), allocatable :: aero(:)   ! object contains information about the aerosols

!=========================================================================================
contains
!=========================================================================================

subroutine zm_conv_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use physics_buffer, only : pbuf_add_field, dtype_r8
  use misc_diagnostics,only: dcape_diags_register

  implicit none

  integer idx

! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx) 

! Flux of snow from deep convection (kg/m2/s) 
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx) 

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)  

! deep gbm cloud liquid water (kg/kg)    
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)  

! vertical velocity (m/s)
   call pbuf_add_field('WUC','global',dtype_r8,(/pcols,pver/), wuc_idx)

! DCAPE-UPL

   if (trigdcape_ull .or. trig_dcape_only) then
    ! temperature from physics in n-1 time step
    call pbuf_add_field('T_STAR','global',dtype_r8,(/pcols,pver/), t_star_idx)
    ! moisturetendency from physics in n-1 time step
    call pbuf_add_field('Q_STAR','global',dtype_r8,(/pcols,pver/), q_star_idx)
   endif


   ! detrained convective cloud water mixing ratio.
   call pbuf_add_field('DLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dlfzm_idx)
   ! detrained convective cloud ice mixing ratio.
   call pbuf_add_field('DIFZM', 'physpkg', dtype_r8, (/pcols,pver/), difzm_idx)

   if (zm_microp) then
      ! Only add the number conc fields if the microphysics is active.

      ! detrained convective cloud water num concen.
      call pbuf_add_field('DNLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnlfzm_idx)
      ! detrained convective cloud ice num concen.
      call pbuf_add_field('DNIFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnifzm_idx)
      ! detrained convective snow num concen.
      call pbuf_add_field('DNSFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnsfzm_idx)
      ! detrained convective snow mixing ratio.
      call pbuf_add_field('DSFZM', 'physpkg', dtype_r8, (/pcols,pver/), dsfzm_idx)
       
   end if

! Variables for dCAPE diagnosis and decomposition

   call dcape_diags_register( pcols )

end subroutine zm_conv_register

!=========================================================================================

subroutine zm_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: outfld, addfld, horiz_only, add_default
  use ppgrid,         only: pcols, pver
  use zm_conv,        only: zm_convi
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use error_messages, only: alloc_err	
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index
  use rad_constituents, only: rad_cnst_get_info 
  use zm_microphysics, only: zm_mphyi


  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces


  logical :: no_deep_pbl    ! if true, no deep convection in PBL
  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

!  integer :: nmodes 

  ! Aerosols
  integer :: i
  character(len=*), parameter :: routine = 'zm_conv_init'

! Allocate the basic aero structure outside the zmconv_microp logical
! This allows the aero structure to be passed
! Note that all of the arrays inside this structure are conditionally allocated

  allocate(aero(begchunk:endchunk))

! 
! Register fields with the output buffer
!


    call addfld ('PRECZ',horiz_only,    'A','m/s','total precipitation from ZM convection')
    call addfld ('ZMDT',(/ 'lev' /), 'A','K/s','T tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDQ',(/ 'lev' /), 'A','kg/kg/s','Q tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDICE',(/ 'lev' /), 'A','kg/kg/s','Cloud ice tendency - Zhang-McFarlane convection')
    call addfld ('ZMDLIQ',(/ 'lev' /), 'A','kg/kg/s','Cloud liq tendency - Zhang-McFarlane convection')
    call addfld ('EVAPTZM',(/ 'lev' /), 'A','K/s','T tendency - Evaporation/snow prod from Zhang convection')
    call addfld ('FZSNTZM',(/ 'lev' /), 'A','K/s','T tendency - Rain to snow conversion from Zhang convection')
    call addfld ('EVSNTZM',(/ 'lev' /), 'A','K/s','T tendency - Snow to rain prod from Zhang convection')
    call addfld ('EVAPQZM',(/ 'lev' /), 'A','kg/kg/s','Q tendency - Evaporation from Zhang-McFarlane moist convection')
    
    call addfld ('ZMFLXPRC',(/ 'ilev' /), 'A','kg/m2/s','Flux of precipitation from ZM convection'       )
    call addfld ('ZMFLXSNW',(/ 'ilev' /), 'A','kg/m2/s','Flux of snow from ZM convection'                )
    call addfld ('ZMNTPRPD',(/ 'lev' /) , 'A','kg/kg/s','Net precipitation production from ZM convection')
    call addfld ('ZMNTSNPD',(/ 'lev' /) , 'A','kg/kg/s','Net snow production from ZM convection'         )
    call addfld ('ZMEIHEAT',(/ 'lev' /) , 'A','W/kg'    ,'Heating by ice and evaporation in ZM convection')
    
    call addfld ('CMFMCDZM',(/ 'ilev' /),'A','kg/m2/s','Convection mass flux from ZM deep ')
    call addfld ('PRECCDZM',horiz_only,    'A','m/s','Convective precipitation rate from ZM deep')

    call addfld ('PCONVB',horiz_only , 'A','Pa'    ,'convection base pressure')
    call addfld ('PCONVT',horiz_only , 'A','Pa'    ,'convection top  pressure')

    call addfld ('MAXI',horiz_only , 'A','level'    ,'model level of launching parcel')

    call addfld ('CAPE_ZM',horiz_only, 'A',   'J/kg', 'Convectively available potential energy')
    call addfld ('DCAPE',  horiz_only, 'A',   'J/kg', 'change rate of Convectively available potential energy')
    call addfld ('FREQZM',horiz_only  ,'A','fraction', 'Fractional occurance of ZM convection') 

    call addfld ('ZMMTT',     (/ 'lev' /), 'A', 'K/s', 'T tendency - ZM convective momentum transport')
    call addfld ('ZMMTU',    (/ 'lev' /), 'A',  'm/s2', 'U tendency - ZM convective momentum transport')
    call addfld ('ZMMTV',    (/ 'lev' /), 'A',  'm/s2', 'V tendency - ZM convective momentum transport')

    call addfld ('ZMMU', (/ 'lev' /), 'A',   'kg/m2/s', 'ZM convection updraft mass flux')
    call addfld ('ZMMD', (/ 'lev' /), 'A',   'kg/m2/s', 'ZM convection downdraft mass flux')

    call addfld ('ZMUPGU',    (/ 'lev' /), 'A', 'm/s2', 'zonal force from ZM updraft pressure gradient term')
    call addfld ('ZMUPGD',    (/ 'lev' /), 'A', 'm/s2', 'zonal force from ZM downdraft pressure gradient term')
    call addfld ('ZMVPGU',    (/ 'lev' /), 'A', 'm/s2', 'meridional force from ZM updraft pressure gradient term')
    call addfld ('ZMVPGD',    (/ 'lev' /), 'A', 'm/s2', 'merdional force from ZM downdraft pressure gradient term')

    call addfld ('ZMICUU',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud U updrafts')
    call addfld ('ZMICUD',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud U downdrafts')
    call addfld ('ZMICVU',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud V updrafts')
    call addfld ('ZMICVD',     (/ 'lev' /), 'A', 'm/s', 'ZM in-cloud V downdrafts')

    if (MCSP) then 
       call addfld ('MCSP_DT',(/ 'lev' /), 'A','K/s','T tedency due to MCSP')
       call addfld ('MCSP_freq',horiz_only, 'A','1','frequency of MCSP activated')
       call addfld ('MCSP_DU',(/ 'lev' /), 'A','m/s/day','U tedency due to MCSP')
       call addfld ('MCSP_DV',(/ 'lev' /), 'A','m/s/day','V tedency due to MCSP')
       call addfld ('ZM_freq',horiz_only, 'A','1','frequency for ZM to be activated')
       call addfld ('ZM_depth',horiz_only,'A','Pa','ZM depth')
       call addfld ('MCSP_shear',horiz_only,'A','m/s','vertical zonal wind shear')
    end if



    if (zm_microp) then

       call addfld ('CLDLIQZM',(/ 'lev' /), 'A', 'g/m3', 'Cloud liquid water - ZM convection')
       call addfld ('CLDICEZM',(/ 'lev' /), 'A', 'g/m3', 'Cloud ice water - ZM convection')
       call addfld ('CLIQSNUM',(/ 'lev' /), 'A', '1'   , 'Cloud liquid water sample number - ZM convection')
       call addfld ('CICESNUM',(/ 'lev' /), 'A', '1'   , 'Cloud ice water sample number - ZM convection')
       call addfld ('QRAINZM' ,(/ 'lev' /), 'A', 'g/m3', 'rain water - ZM convection')
       call addfld ('QSNOWZM' ,(/ 'lev' /), 'A', 'g/m3', 'snow - ZM convection')
       call addfld ('QGRAPZM' ,(/ 'lev' /), 'A', 'g/m3', 'graupel - ZM convection')
       call addfld ('CRAINNUM',(/ 'lev' /), 'A', '1'   , 'Cloud rain water sample number - ZM convection')
       call addfld ('CSNOWNUM',(/ 'lev' /), 'A', '1'   , 'Cloud snow sample number - ZM convection')
       call addfld ('CGRAPNUM',(/ 'lev' /), 'A', '1'   , 'Cloud graupel sample number -ZM convection')

       call addfld ('DIFZM',(/ 'lev' /), 'A','kg/kg/s ', 'Detrained ice water from ZM convection')
       call addfld ('DLFZM',(/ 'lev' /), 'A','kg/kg/s ', 'Detrained liquid water from ZM convection')
       call addfld ('DNIFZM',(/ 'lev' /), 'A','1/kg/s ', 'Detrained ice water num concen from ZM convection')
       call addfld ('DNLFZM',(/ 'lev' /), 'A','1/kg/s ', 'Detrained liquid water num concen from ZM convection')
       call addfld ('WUZM',(/ 'lev' /)  , 'A','m/s'    , 'vertical velocity - ZM convection')
       call addfld ('WUZMSNUM',(/ 'lev' /),'A','1'     , 'vertical velocity sample number - ZM convection')

       call addfld ('QNLZM',(/ 'lev' /), 'A','1/m3'    , 'Cloud liquid water number concen - ZM convection')
       call addfld ('QNIZM',(/ 'lev' /), 'A','1/m3'    , 'Cloud ice number concen - ZM convection')
       call addfld ('QNRZM',(/ 'lev' /), 'A','1/m3'    , 'Cloud rain water number concen - ZM convection')
       call addfld ('QNSZM',(/ 'lev' /), 'A','1/m3'    , 'Cloud snow number concen - ZM convection')
       call addfld ('QNGZM',(/ 'lev' /), 'A','1/m3'    , 'Cloud graupel number concen -ZM convection')

       call addfld ('FRZZM',(/ 'lev' /), 'A','K/s'     , 'heating tendency due to freezing - ZM convection')

       call addfld ('AUTOL_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to autoconversion of droplets to rain')
       call addfld ('ACCRL_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to accretion of droplets by rain')
       call addfld ('BERGN_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to Bergeron process')
       call addfld ('FHTIM_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to immersion freezing')
       call addfld ('FHTCT_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to contact freezing')
       call addfld ('FHML_M'  ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to homogeneous freezing of droplet')
       call addfld ('HMPI_M'  ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to HM process')
       call addfld ('ACCSL_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to accretion of droplet by snow')
       call addfld ('DLF_M'   ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to detrainment of droplet')
       call addfld ('COND_M'  ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to condensation')

       call addfld ('AUTOL_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to autoconversion of droplets to rain')
       call addfld ('ACCRL_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to accretion of droplets by rain')
       call addfld ('BERGN_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to Bergeron process')
       call addfld ('FHTIM_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to immersion freezing')
       call addfld ('FHTCT_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to contact freezing')
       call addfld ('FHML_N'  ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to homogeneous freezing of droplet')
       call addfld ('ACCSL_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to accretion of droplet by snow')
       call addfld ('ACTIV_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to droplets activation')
       call addfld ('DLF_N'   ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to detrainment of droplet')

       call addfld ('AUTOI_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to autoconversion of ice to snow')
       call addfld ('ACCSI_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to accretion of ice by snow')
       call addfld ('DIF_M'   ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to detrainment of cloud ice')
       call addfld ('DEPOS_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to deposition')

       call addfld ('NUCLI_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to ice nucleation')
       call addfld ('AUTOI_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to autoconversion of ice to snow')
       call addfld ('ACCSI_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to accretion of ice by snow')
       call addfld ('HMPI_N'  ,(/ 'lev' /), 'A', '1/kg/s' , 'num tendency due to HM process')
       call addfld ('DIF_N'   ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to detrainment of cloud ice')
       call addfld ('TRSPC_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of droplets due to convective transport')
       call addfld ('TRSPC_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of droplets due to convective transport')
       call addfld ('TRSPI_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of ice crystal due to convective transport')
       call addfld ('TRSPI_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of ice crystal due to convective transport')

       call addfld ('ACCGR_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to collection of rain by graupel')
       call addfld ('ACCGL_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to collection of droplets by graupel')
       call addfld ('ACCGSL_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of graupel due to collection of droplets by snow')
       call addfld ('ACCGSR_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of graupel due to collection of rain by snow')
       call addfld ('ACCGIR_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of graupel due to collection of rain by ice')
       call addfld ('ACCGRI_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of graupel due to collection of ice by rain')
       call addfld ('ACCGRS_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of graupel due to collection of snow by rain')

       call addfld ('ACCGSL_N',(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of graupel due to collection of droplets by snow')
       call addfld ('ACCGSR_N',(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of graupel due to collection of rain by snow')
       call addfld ('ACCGIR_N',(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of graupel due to collection of rain by ice')

       call addfld ('ACCSRI_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of snow due to collection of ice by rain')
       call addfld ('ACCIGL_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of ice mult(splintering) due to acc droplets by graupel')
       call addfld ('ACCIGR_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of ice mult(splintering) due to acc rain by graupel')
       call addfld ('ACCSIR_M',(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of snow due to collection of rain by ice')

       call addfld ('ACCIGL_N',(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of ice mult(splintering) due to acc droplets by graupel')
       call addfld ('ACCIGR_N',(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of ice mult(splintering) due to acc rain by graupel')
       call addfld ('ACCSIR_N',(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of snow due to collection of rain by ice')
       call addfld ('ACCGL_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to collection of droplets by graupel')
       call addfld ('ACCGR_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency due to collection of rain by graupel')
       call addfld ('ACCIL_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of cloud ice due to collection of droplet by cloud ice')
       call addfld ('ACCIL_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of cloud ice due to collection of droplet by cloud ice')

       call addfld ('FALLR_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of rain fallout')
       call addfld ('FALLS_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of snow fallout')
       call addfld ('FALLG_M' ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency of graupel fallout')
       call addfld ('FALLR_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of rain fallout')
       call addfld ('FALLS_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of snow fallout')
       call addfld ('FALLG_N' ,(/ 'lev' /), 'A', '1/kg/m' , 'num tendency of graupel fallout')
       call addfld ('FHMR_M'  ,(/ 'lev' /), 'A', 'kg/kg/m', 'mass tendency due to homogeneous freezing of rain')

       call addfld ('PRECZ_SN',horiz_only , 'A', '#'      , 'sample num of convective precipitation rate from ZM deep')

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
       call add_default ('FRZZM',   1, ' ')

    end if


    
    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       convproc_do_aer_out = convproc_do_aer, & 
                       convproc_do_gas_out = convproc_do_gas)   
    ! Determine whether its a 'modal' aerosol simulation  or not
    call rad_cnst_get_info(0, nmodes=nmodes)
    clim_modal_aero = (nmodes > 0)

    if ( history_budget ) then
       call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
       call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
       call add_default('ZMDT     ', history_budget_histfile_num, ' ')
       call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
       call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
       call add_default('ZMDICE   ', history_budget_histfile_num, ' ')

       if( cam_physpkg_is('cam4') .or. cam_physpkg_is('default') ) then
          call add_default('ZMMTT    ', history_budget_histfile_num, ' ')
       end if

    end if
!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
    limcnv = 0   ! null value to check against below
    if (pref_edge(1) >= 4.e3_r8) then
       limcnv = 1
    else
       do k=1,plev
          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = plevp
    end if
    
    if (masterproc) then
       write(iulog,*)'ZM_CONV_INIT: Deep convection will be capped at intfc ',limcnv, &
            ' which is ',pref_edge(limcnv),' pascals'
    end if
        
    no_deep_pbl = phys_deepconv_pbl()
    call zm_convi(limcnv,no_deep_pbl_in = no_deep_pbl)

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

       ! if zmconv_microp not enabled then use old estimate of snow production in
       ! zm_conv_evap
       old_snow = .false. 

       ! Initialize the aerosol object with data from the modes/species
       ! affecting climate,
       ! i.e., the list index is hardcoded to 0.

       call rad_cnst_get_info(0, nmodes=nmodes, naero=nbulk)


       do i = begchunk, endchunk
          call zm_aero_init(nmodes, nbulk, aero(i))
       end do

       if (nmodes > 0) then

          dgnum_idx = pbuf_get_index('DGNUM')

       else if (nbulk > 0) then

          ! This call is needed to allow running the ZM microphysics with the
          ! cam4 physics package.
          call ndrop_bam_init()

       end if

     end if ! zmconv_microp

    !-------------------------------------------------------------------------------------
    contains
    !-------------------------------------------------------------------------------------

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
       !----------------------------------------------------------------------------------

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
          if (aero%mode_accum_idx == -1 .or. aero%mode_aitken_idx == -1 .or. aero%mode_coarse_idx == -1) then
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
          ! Check that required modal specie types were found
          if (aero%coarse_dust_idx == -1 .or. aero%coarse_nacl_idx == -1) then
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
             call rad_cnst_get_mode_props(0, m, &
                sigmag=sigmag, dgnumlo=dgnumlo, dgnumhi=dgnumhi)

             alnsg               = log(sigmag)
             aero%voltonumblo(m) = 1._r8 / ( (pi/6._r8)*(dgnumlo**3._r8)*exp(4.5_r8*alnsg**2._r8) )
             aero%voltonumbhi(m) = 1._r8 / ( (pi/6._r8)*(dgnumhi**3._r8)*exp(4.5_r8*alnsg**2._r8) )

             ! save sigmag of aitken mode
             if (m == aero%mode_aitken_idx) aero%sigmag_aitken = sigmag

             ! Properties of modal species
             do l = 1, aero%nspec(m)
                call rad_cnst_get_aer_props(0, m, l, density_aer=aero%specdens(l,m), &
                   hygro_aer=aero%spechygro(l,m))
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
             if (trim(aername(iaer)) == 'SULFATE') aero%idxsul = iaer
             if (trim(aername(iaer)) == 'DUST1')   aero%idxdst1 = iaer
             if (trim(aername(iaer)) == 'DUST2')   aero%idxdst2 = iaer
             if (trim(aername(iaer)) == 'DUST3')   aero%idxdst3 = iaer
             if (trim(aername(iaer)) == 'DUST4')   aero%idxdst4 = iaer
             if (trim(aername(iaer)) == 'BCPHI')   aero%idxbcphi = iaer
          end do

       end if

    end subroutine zm_aero_init

end subroutine zm_conv_init
!=========================================================================================
!subroutine zm_conv_tend(state, ptend, tdt)

subroutine zm_conv_tend(pblh    ,mcon    ,cme     , &
     tpert   ,dlftot  ,pflx    ,zdu      , &
     rliq    ,rice    ,&
     ztodt   , &
     jctop   ,jcbot , &
     state   ,ptend_all   ,landfrac,  pbuf, mu, eu, &
     du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath) 

   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init
   use physics_update_mod, only: physics_update
   use physics_types, only: physics_state_copy, physics_state_dealloc
   use physics_types, only: physics_ptend_sum, physics_ptend_dealloc

   use phys_grid,     only: get_lat_p, get_lon_p
   use time_manager,  only: get_nstep, is_first_step
   use physics_buffer, only : pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use physconst,     only: gravit
   use phys_control,  only: cam_physpkg_is
   use time_manager,       only: get_curr_date
   use interpolate_data, only: vertinterp


   ! Arguments
   type(physics_state), target, intent(in ) :: state      ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend_all      ! individual parameterization tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblh(pcols)                 ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: landfrac(pcols)             ! RBN - Landfrac 

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: dlftot(pcols,pver) ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out) :: rice(pcols) ! reserved ice (not yet in cldice) for energy integrals
   real(r8), intent(out):: mu(pcols,pver) 
   real(r8), intent(out):: eu(pcols,pver) 
   real(r8), intent(out):: du(pcols,pver) 
   real(r8), intent(out):: md(pcols,pver) 
   real(r8), intent(out):: ed(pcols,pver) 
   real(r8), intent(out):: dp(pcols,pver) 
   
   ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out):: dsubcld(pcols) 
   
   ! wg layer thickness in mbs between lcl and maxi.    
   integer, intent(out) :: jt(pcols)   
   
   ! wg top  level index of deep cumulus convection.
   integer, intent(out) :: maxg(pcols) 
   
   ! wg gathered values of maxi.
   integer, intent(out) :: ideep(pcols)
   
   ! w holds position of gathered points vs longitude index   
   integer, intent(out)  :: lengath
   ! Local variables

    type(zm_microp_st)        :: microp_st 

   integer :: i,k,l,m
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer :: itim_old                ! for physics buffer fields

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8) :: tend_s_snwprd  (pcols,pver) ! Heating rate of snow production
   real(r8) :: tend_s_snwevmlt(pcols,pver) ! Heating rate of evap/melting of snow
   real(r8) :: fake_dpdry(pcols,pver) ! used in convtran call

   ! physics types
   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
   type(physics_ptend),target :: ptend_loc     ! package tendencies
   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow from ZM convection 
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation
   real(r8), pointer, dimension(:,:) :: flxprec      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: flxsnow      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: dp_cldliq
   real(r8), pointer, dimension(:,:) :: dp_cldice
   real(r8), pointer, dimension(:,:) :: wuc

   ! DCAPE-ULL
   real(r8), pointer, dimension(:,:) :: t_star ! intermediate T between n and n-1 time step
   real(r8), pointer, dimension(:,:) :: q_star ! intermediate q between n and n-1 time step
   real(r8) :: dcape(pcols)                    ! dynamical cape
   real(r8) :: maxgsav(pcols)                  ! tmp array for recording and outfld to MAXI

   real(r8), pointer :: dlf(:,:)    ! detrained convective cloud water mixing ratio.
   real(r8), pointer :: dif(:,:)    ! detrained convective cloud ice mixing ratio.
   real(r8), pointer :: dsf(:,:)    ! detrained convective snow mixing ratio.
   real(r8), pointer :: dnlf(:,:)   ! detrained convective cloud water num concen.
   real(r8), pointer :: dnif(:,:)   ! detrained convective cloud ice num concen.
   real(r8), pointer :: dnsf(:,:)   ! detrained convective snow num concen.

   real(r8), pointer :: lambdadpcu(:,:) ! slope of cloud liquid size distr
   real(r8), pointer :: mudpcu(:,:)     ! width parameter of droplet size distr
   real(r8), pointer :: qi(:,:)         ! wg grid slice of cloud ice.

   integer :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   integer :: jcbot(pcols)  ! o row of base of cloud indices passed out.

   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols)

   ! history output fields
   real(r8) :: cape(pcols)        ! w  convective available potential energy.
   real(r8) :: mu_out(pcols,pver)
   real(r8) :: md_out(pcols,pver)


   ! used in momentum transport calculation
   real(r8) :: winds(pcols, pver, 2)
   real(r8) :: wind_tends(pcols, pver, 2)
   real(r8) :: pguall(pcols, pver, 2)
   real(r8) :: pgdall(pcols, pver, 2)
   real(r8) :: icwu(pcols,pver, 2)
   real(r8) :: icwd(pcols,pver, 2)
   real(r8) :: seten(pcols, pver)
   logical  :: l_windt(2)
   real(r8) :: tfinal1, tfinal2
   integer  :: ii

   logical  :: lq(pcnst)

   real(r8) :: alpha2, alpha_moisture, alphau, alphav, top, bottom
   real(r8) :: Q_dis(pcols), Qc_adjust, Qq_dis(pcols), Qcq_adjust
   real(r8) :: Qm(pcols,pver), Qmq(pcols,pver), Qmu(pcols,pver), Qmv(pcols,pver)
   real(r8) :: Qc_int_start(pcols), Qcq_int_start(pcols)
   real(r8) :: Qm_int_end(pcols), Qmq_int_end(pcols), Pa_int_end(pcols)
   real(r8) :: Qs_zmconv(pcols), Qv_zmconv(pcols)

   real(r8) :: MCSP_freq(pcols), MCSP_DT(pcols,pver)
   real(r8) :: ZM_depth(pcols), MCSP_shear(pcols)
   real(r8) :: du600(pcols), dv600(pcols)
   integer :: iyr, imon, iday, isec

   logical  :: doslop
   logical  :: doslop_heat
   logical  :: doslop_moisture
   logical  :: doslop_uwind
   logical  :: doslop_vwind


   real(r8) :: sprd(pcols,pver)
   real(r8) :: frz(pcols,pver)
   real(r8)  precz_snum(pcols)


   if (zm_microp) then
     allocate( &
       microp_st%wu(pcols,pver),      &  ! vertical velocity
       microp_st%qliq(pcols,pver),    &  ! convective cloud liquid water.
       microp_st%qice(pcols,pver),    &  ! convective cloud ice.
       microp_st%qrain(pcols,pver),   &  ! convective rain water.
       microp_st%qsnow(pcols,pver),   &  ! convective snow.
       microp_st%qgraupel(pcols,pver),&  ! convective graupel
       microp_st%qnl(pcols,pver),     &  ! convective cloud liquid water num concen.
       microp_st%qni(pcols,pver),     &  ! convective cloud ice num concen.
       microp_st%qnr(pcols,pver),     &  ! convective rain water num concen.
       microp_st%qns(pcols,pver),     &  ! convective snow num concen.
       microp_st%qng(pcols,pver),     &  ! convective graupel num concen.
       microp_st%autolm(pcols,pver),  &  !mass tendency due to autoconversion of droplets to rain
       microp_st%accrlm(pcols,pver),  &  !mass tendency due to accretion of droplets by rain
       microp_st%bergnm(pcols,pver),  &  !mass tendency due to Bergeron process
       microp_st%fhtimm(pcols,pver),  &  !mass tendency due to immersion freezing
       microp_st%fhtctm(pcols,pver),  &  !mass tendency due to contact freezing
       microp_st%fhmlm (pcols,pver),  &  !mass tendency due to homogeneous freezing
       microp_st%hmpim (pcols,pver),  &  !mass tendency due to HM process
       microp_st%accslm(pcols,pver),  &  !mass tendency due to accretion of droplets by snow
       microp_st%dlfm  (pcols,pver),  &  !mass tendency due to detrainment of droplet
       microp_st%autoln(pcols,pver),  &  !num tendency due to autoconversion of droplets to rain
       microp_st%accrln(pcols,pver),  &  !num tendency due to accretion of droplets by rain
       microp_st%bergnn(pcols,pver),  &  !num tendency due to Bergeron process
       microp_st%fhtimn(pcols,pver),  &  !num tendency due to immersion freezing
       microp_st%fhtctn(pcols,pver),  &  !num tendency due to contact freezing
       microp_st%fhmln (pcols,pver),  &  !num tendency due to homogeneous freezing
       microp_st%accsln(pcols,pver),  &  !num tendency due to accretion of droplets by snow
       microp_st%activn(pcols,pver),  &  !num tendency due to droplets activation
       microp_st%dlfn  (pcols,pver),  &  !num tendency due to detrainment of droplet
       microp_st%autoim(pcols,pver),  &  !mass tendency due to autoconversion of cloud ice to snow
       microp_st%accsim(pcols,pver),  &  !mass tendency due to accretion of cloud ice by snow
       microp_st%difm  (pcols,pver),  &  !mass tendency due to detrainment of cloud ice
       microp_st%nuclin(pcols,pver),  &  !num tendency due to ice nucleation
       microp_st%autoin(pcols,pver),  &  !num tendency due to autoconversion of cloud ice to snow
       microp_st%accsin(pcols,pver),  &  !num tendency due to accretion of cloud ice by snow
       microp_st%hmpin (pcols,pver),  &  !num tendency due to HM process
       microp_st%difn  (pcols,pver),  &  !num tendency due to detrainment of cloud ice
       microp_st%cmel  (pcols,pver),  &  !mass tendency due to condensation
       microp_st%cmei  (pcols,pver),  &  !mass tendency due to deposition
       microp_st%trspcm(pcols,pver),  &  !LWC tendency due to convective transport
       microp_st%trspcn(pcols,pver),  &  !droplet num tendency due to convective transport
       microp_st%trspim(pcols,pver),  &  !IWC tendency due to convective transport
       microp_st%trspin(pcols,pver),  &  !ice crystal num tendency due to convective transport
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
       microp_st%fallrm(pcols, pver),  & ! mass tendency of rain fallout
       microp_st%fallsm(pcols, pver),  & ! mass tendency of snow fallout
       microp_st%fallgm(pcols, pver),  & ! mass tendency of graupel fallout
       microp_st%fallrn(pcols, pver),  & ! num tendency of rain fallout
       microp_st%fallsn(pcols, pver),  & ! num tendency of snow fallout
       microp_st%fallgn(pcols, pver),  & ! num tendency of graupel fallout
       microp_st%fhmrm (pcols,pver) )    ! mass tendency due to homogeneous freezing of rain
    end if


   doslop          = .false.
   doslop_heat     = .false.
   doslop_moisture = .false.
   doslop_uwind    = .false.
   doslop_vwind    = .false.
   alphau = 0.0_r8
   alphav = 0.0_r8
   if( MCSP ) then
        if( MCSP_heat_coeff > 0._r8 ) then
                doslop_heat = .true.
                alpha2 = MCSP_heat_coeff
        end if
        if( MCSP_moisture_coeff > 0._r8 ) then
                doslop_moisture = .true.
                alpha_moisture = MCSP_moisture_coeff
        end if
        if( MCSP_uwind_coeff > 0._r8 ) then
                doslop_uwind = .true.
                alphau = MCSP_uwind_coeff
        end if
        if( MCSP_vwind_coeff > 0._r8 ) then
                doslop_vwind = .true.
                alphav = MCSP_vwind_coeff
        end if
    end if

   if (doslop_heat .or. doslop_moisture .or. doslop_uwind .or. doslop_vwind) then
      doslop = .true.
   endif


   !----------------------------------------------------------------------

   ! initialize
   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem = 0._r8   
   mu_out(:,:) = 0._r8
   md_out(:,:) = 0._r8
   dlftot(:,:) = 0._r8
   wind_tends(:ncol,:pver,:) = 0.0_r8

   call physics_state_copy(state,state1)             ! copy state to local state1.

   lq(:) = .FALSE.
   lq(1) = .TRUE.
 
   if (doslop) then
        call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lu = doslop_uwind, lv = doslop_vwind, lq=lq)! initialize local ptend type
   else
        call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq)! initialize local ptend type
   end if

!
! Associate pointers with physics buffer fields
!
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, icwmrdp_idx,     ql )
   call pbuf_get_field(pbuf, rprddp_idx,      rprd )
   call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/),    kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
   call pbuf_get_field(pbuf, prec_dp_idx,     prec )
   call pbuf_get_field(pbuf, snow_dp_idx,     snow )


! DCAPE-ULL
   if(trigdcape_ull .or. trig_dcape_only)then
     call pbuf_get_field(pbuf, t_star_idx,     t_star)
     call pbuf_get_field(pbuf, q_star_idx,     q_star)
     if ( is_first_step()) then
       q_star(:ncol,:pver) = state%q(:ncol,:pver,1)
       t_star(:ncol,:pver) = state%t(:ncol,:pver)
     end if
   end if

   call pbuf_get_field(pbuf, icimrdp_idx,     qi )
   call pbuf_get_field(pbuf, dlfzm_idx,  dlf)
   call pbuf_get_field(pbuf, difzm_idx,  dif)

   if (zm_microp) then
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlf)
      call pbuf_get_field(pbuf, dnifzm_idx, dnif)
      call pbuf_get_field(pbuf, dsfzm_idx,  dsf)
      call pbuf_get_field(pbuf, dnsfzm_idx, dnsf)
      call pbuf_get_field(pbuf, wuc_idx, wuc  )
   else
      allocate(dnlf(pcols,pver), dnif(pcols,pver),dsf(pcols,pver), dnsf(pcols,pver))
      allocate(wuc(pcols,pver))
   end if
   wuc(:,:) = 0._r8

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




!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')
   call zm_convr(   lchnk   ,ncol    , &
                    state%t       ,state%q(:,:,1)     ,prec    ,jctop   ,jcbot   , &
                    pblh    ,state%zm      ,state%phis    ,state%zi      ,ptend_loc%q(:,:,1)    , &
                    ptend_loc%s    ,state%pmid     ,state%pint    ,state%pdel     ,state%omega  , &
                    .5_r8*ztodt    ,mcon    ,cme     , cape,      &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu,md,du,eu,ed      , &
                    dp ,dsubcld ,jt,maxg,ideep   , &
                    lengath ,ql      ,rliq  ,landfrac,  &
                    t_star, q_star, dcape, &  
                    aero(lchnk), qi, dif, dnlf, dnif, dsf, dnsf, sprd, rice, frz, mudpcu, &
                    lambdadpcu,  microp_st, wuc)

   if (zm_microp) then
     dlftot(:ncol,:pver) = dlf(:ncol,:pver) + dif(:ncol,:pver) + dsf(:ncol,:pver)
   else
     dlftot(:ncol,:pver) = dlf(:ncol,:pver)
   end if

   
   call t_stopf ('zm_convr')


   if (doslop) then
  !   
  !  Begin MCSP parameterization here 
  !
     du600 = 0._r8
     dv600 = 0._r8
     ZM_depth = 0._r8
     MCSP_shear = 0._r8

     call vertinterp(ncol, pcols, pver, state%pmid, 60000._r8, state%u,du600)
     call vertinterp(ncol, pcols, pver, state%pmid, 60000._r8, state%v,dv600)

     do i=1,ncol
        if(state%pmid(i,pver).gt.60000._r8) then
          du600(i) = du600(i)-state%u(i,pver)
          dv600(i) = dv600(i)-state%v(i,pver)
        else
          du600(i) = -999._r8
          dv600(i) = -999._r8
        end if
     end do

     MCSP_shear = du600
     do i=1,ncol
        if( jctop(i).ne.pver ) ZM_depth(i) = state%pint(i,pver+1)-state%pmid(i,jctop(i))
     end do

     ! Define parameters

     Qs_zmconv(:) = 0._r8
     Qv_zmconv(:) = 0._r8
     Pa_int_end(:) = 0._r8
     do i = 1,ncol
       do k = jctop(i),pver

         Qs_zmconv(i) = Qs_zmconv(i) + ptend_loc%s(i,k) * state%pdel(i,k)
         Qv_zmconv(i) = Qv_zmconv(i) + ptend_loc%q(i,k,1) * state%pdel(i,k)
         Pa_int_end(i) = Pa_int_end(i) + state%pdel(i,k)

       enddo
     enddo

     do i = 1,ncol
       if (jctop(i) .ne. pver) then
         Qs_zmconv(i) = Qs_zmconv(i) / Pa_int_end(i)
         Qv_zmconv(i) = Qv_zmconv(i) / Pa_int_end(i)
       else
         Qs_zmconv(i) = 0._r8
         Qv_zmconv(i) = 0._r8
       endif
     enddo

     Qm(:,:) = 0.0_r8
     Qmq(:,:) = 0.0_r8
     Qmu(:,:) = 0.0_r8
     Qmv(:,:) = 0.0_r8

     do i = 1,ncol

      Qm_int_end(i) = 0._r8
      Qmq_int_end(i) = 0._r8
      Pa_int_end(i) = 0._r8

      Q_dis(i) = 0._r8
      Qq_dis(i) = 0._r8

      if( (state%pint(i,pver+1)-state%pmid(i,jctop(i))) >= 70000._r8 ) then
      if(abs(du600(i)).ge.3._r8.and.abs(du600(i)).lt.200.and.Qs_zmconv(i).gt.0._r8) then

       do k = jctop(i),pver

           top = state%pint(i,pver+1) - state%pmid(i,k)
           bottom = state%pint(i,pver+1) - state%pmid(i,jctop(i))

           Qm(i,k) = -1._r8 * Qs_zmconv(i) * alpha2 * sin(2.0_r8*pi*(top/bottom))
           Qmq(i,k) = -1._r8 * Qv_zmconv(i) * alpha2 * sin(2.0_r8*pi*(top/bottom))
           Qmq(i,k) = Qm(i,k)/2500000._r8/4._r8


           Qmu(i,k) = alphau * (cos(pi*(top/bottom)))

           Qmv(i,k) = alphav * (cos(pi*(top/bottom)))

           Qm_int_end(i) = Qm_int_end(i) + Qm(i,k) * (state%pdel(i,k)/gravit)
           Qm_int_end(i) = Qm_int_end(i) + (2._r8*Qmu(i,k)*ztodt*state%u(i,k)+ &
                                            Qmu(i,k)*Qmu(i,k)*ztodt*ztodt)/2._r8 * (state%pdel(i,k)/gravit)/ztodt
           Qm_int_end(i) = Qm_int_end(i) + (2._r8*Qmv(i,k)*ztodt*state%v(i,k)+ &
                                            Qmv(i,k)*Qmv(i,k)*ztodt*ztodt)/2._r8 * (state%pdel(i,k)/gravit)/ztodt
           Qmq_int_end(i) = Qmq_int_end(i) + Qmq(i,k) * (state%pdel(i,k)/gravit)
           Pa_int_end(i) = Pa_int_end(i) + state%pdel(i,k)
        end do

        Q_dis(i) = Qm_int_end(i) / Pa_int_end(i)
        Qq_dis(i) = Qmq_int_end(i) / Pa_int_end(i)

       endif
       endif
     enddo

     MCSP_DT = 0._r8
     MCSP_freq = 0._r8

     do i = 1,ncol
       do k = jctop(i),pver
         Qcq_adjust = ptend_loc%q(i,k,1) - Qq_dis(i) * gravit

         ptend_loc%s(i,k) = ptend_loc%s(i,k) - Q_dis(i)*gravit ! energy fixer

         MCSP_DT(i,k) = -Q_dis(i)*gravit+Qm(i,k)
         if(abs(Qm(i,k)).gt.0._r8 .and. abs(Qmu(i,k)).gt.0._r8) MCSP_freq(i) = 1._r8

         if (doslop_heat) ptend_loc%s(i,k) = ptend_loc%s(i,k) + Qm(i,k)
         if (doslop_moisture) ptend_loc%q(i,k,1) = Qcq_adjust + Qmq(i,k)
         if (doslop_uwind) ptend_loc%u(i,k) = Qmu(i,k)
         if (doslop_vwind) ptend_loc%v(i,k) = Qmv(i,k)
       enddo
     enddo

   
   !   
   ! End the MCSP parameterization here 
   !   
      

      MCSP_DT(:ncol,:pver) = MCSP_DT(:ncol,:pver)/cpair
      call outfld('MCSP_DT    ',MCSP_DT           ,pcols   ,lchnk   )
      call outfld('MCSP_freq    ',MCSP_freq           ,pcols   ,lchnk   )
      if (doslop_uwind) call outfld('MCSP_DU    ',ptend_loc%u*86400._r8           ,pcols   ,lchnk   )
      if (doslop_vwind) call outfld('MCSP_DV    ',ptend_loc%v*86400._r8           ,pcols   ,lchnk   )
      call outfld('ZM_depth   ',ZM_depth                        ,pcols   ,lchnk   )
      call outfld('MCSP_shear ',MCSP_shear                      ,pcols   ,lchnk   )

   end if

   call outfld('DCAPE', dcape, pcols, lchnk)
   call outfld('CAPE_ZM', cape, pcols, lchnk)        ! RBN - CAPE output
!
! Output fractional occurance of ZM convection
!
   freqzm(:) = 0._r8
   do i = 1,lengath
      freqzm(ideep(i)) = 1.0_r8
   end do
   call outfld('FREQZM  ',freqzm          ,pcols   ,lchnk   )
!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   mcon(:ncol,:pver) = mcon(:ncol,:pver) * 100._r8/gravit

   call outfld('CMFMCDZM', mcon, pcols, lchnk)

   ! Store upward and downward mass fluxes in un-gathered arrays
   ! + convert from mb/s to kg/m^2/s
   do i=1,lengath
      do k=1,pver
         ii = ideep(i)
         mu_out(ii,k) = mu(i,k) * 100._r8/gravit
         md_out(ii,k) = md(i,k) * 100._r8/gravit
      end do
   end do


   if(convproc_do_aer .or. convproc_do_gas) then 
      call outfld('ZMMU', mu_out,      pcols, lchnk)
      call outfld('ZMMD', md_out,      pcols, lchnk)
   else
      call outfld('ZMMU', mu_out(1,1), pcols, lchnk)
      call outfld('ZMMD', md_out(1,1), pcols, lchnk)
   endif

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )

   if (zm_microp) call zm_conv_micro_outfld(microp_st, dlf, dif, dnlf, dnif, frz, lchnk, ncol)


   maxgsav(:) = 0._r8 ! zero if no convection. true mean to be MAXI/FREQZM
   pcont(:ncol) = state%ps(:ncol)
   pconb(:ncol) = state%ps(:ncol)
   do i = 1,lengath
       if (maxg(i).gt.jt(i)) then
          pcont(ideep(i)) = state%pmid(ideep(i),jt(i))  ! gathered array (or jctop ungathered)
          pconb(ideep(i)) = state%pmid(ideep(i),maxg(i))! gathered array
          maxgsav(ideep(i)) = dble(maxg(i))             ! gathered array for launching level
       endif
       !     write(iulog,*) ' pcont, pconb ', pcont(i), pconb(i), cnt(i), cnb(i)
   end do
   call outfld('PCONVT  ',pcont          ,pcols   ,lchnk   )
   call outfld('PCONVB  ',pconb          ,pcols   ,lchnk   )
   call outfld('MAXI  ',maxgsav          ,pcols   ,lchnk   )


  call physics_ptend_init(ptend_all, state%psetcols, 'zm_conv_tend')

  ! add tendency from this process to tendencies from other processes
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)

  ! initialize ptend for next process
  lq(:) = .FALSE.
  lq(1) = .TRUE.
  call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
!

    call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
    call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
    dp_cldliq(:ncol,:) = 0._r8
    dp_cldice(:ncol,:) = 0._r8

    call t_startf ('zm_conv_evap')
    call zm_conv_evap(state1%ncol,state1%lchnk, &
         state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
         ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
         rprd, cld, ztodt, &
         prec, snow, ntprprd, ntsnprd , flxprec, flxsnow, sprd, old_snow)
    call t_stopf ('zm_conv_evap')

    evapcdp(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)
!
! Write out variables from zm_conv_evap
!
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('EVAPTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwprd  (:ncol,:pver)/cpair
   call outfld('FZSNTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver)/cpair
   call outfld('EVSNTZM ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQZM ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('ZMFLXPRC', flxprec, pcols, lchnk)
   call outfld('ZMFLXSNW', flxsnow, pcols, lchnk)
   call outfld('ZMNTPRPD', ntprprd, pcols, lchnk)
   call outfld('ZMNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('ZMEIHEAT', ptend_loc%s, pcols, lchnk)
   call outfld('PRECCDZM   ',prec,  pcols   ,lchnk   )



   call outfld('PRECZ   ', prec   , pcols, lchnk)

   if (zm_microp) then
      do i = 1,ncol
         if (prec(i) .gt. 0.0_r8) then
            precz_snum(i) = 1.0_r8
         else
            precz_snum(i) = 0.0_r8
         end if
      end do
      call outfld('PRECZ_SN', precz_snum, pcols, lchnk)
   end if

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)


  ! Momentum Transport (non-cam3 physics)

  if ( .not. cam_physpkg_is('cam3')) then

     call physics_ptend_init(ptend_loc, state1%psetcols, 'momtran', ls=.true., lu=.true., lv=.true.)

     winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
     winds(:ncol,:pver,2) = state1%v(:ncol,:pver)
   
     l_windt(1) = .true.
     l_windt(2) = .true.

     call t_startf ('momtran')
     call momtran (lchnk, ncol,                                        &
                   l_windt,winds, 2,  mu(1,1), md(1,1),   &
                   du(1,1), eu(1,1), ed(1,1), dp(1,1), dsubcld(1),  &
                   jt(1),maxg(1), ideep(1), 1, lengath,  &
                   nstep,  wind_tends, pguall, pgdall, icwu, icwd, ztodt, seten )  
     call t_stopf ('momtran')

     ptend_loc%u(:ncol,:pver) = wind_tends(:ncol,:pver,1)
     ptend_loc%v(:ncol,:pver) = wind_tends(:ncol,:pver,2)
     ptend_loc%s(:ncol,:pver) = seten(:ncol,:pver)  

     call physics_ptend_sum(ptend_loc,ptend_all, ncol)

     ! update physics state type state1 with ptend_loc 
     call physics_update(state1, ptend_loc, ztodt)

     ftem(:ncol,:pver) = seten(:ncol,:pver)/cpair
     call outfld('ZMMTT', ftem             , pcols, lchnk)
     call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk)
     call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk)
   
     ! Output apparent force from  pressure gradient
     call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk)
     call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk)
     call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk)
     call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk)

     ! Output in-cloud winds
     call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk)
     call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk)
     call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk)
     call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk)

   end if

   ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)
   call physics_ptend_init(ptend_loc, state1%psetcols, 'convtran1', lq=lq)


   ! dpdry is not used in this call to convtran since the cloud liquid and ice mixing
   ! ratios are moist
   fake_dpdry(:,:) = 0._r8

   call t_startf ('convtran1')
   call convtran (lchnk,                                        &
                  ptend_loc%lq,state1%q, pcnst,  mu, md,   &
                  du, eu, ed, dp, dsubcld,  &
                  jt,maxg, ideep, 1, lengath,  &
                  nstep,   fracis,  ptend_loc%q, fake_dpdry, ztodt)  
   call t_stopf ('convtran1')

   call outfld('ZMDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('ZMDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)

   if (zm_microp) then
     deallocate( &
       microp_st%wu,      &  
       microp_st%qliq,    &  
       microp_st%qice,    &  
       microp_st%qrain,   &  
       microp_st%qsnow,   &  
       microp_st%qgraupel,&  
       microp_st%qnl,     &  
       microp_st%qni,     &  
       microp_st%qnr,     &  
       microp_st%qns,     &  
       microp_st%qng,     &  
       microp_st%autolm,  &  
       microp_st%accrlm,  &  
       microp_st%bergnm,  &  
       microp_st%fhtimm,  &  
       microp_st%fhtctm,  &  
       microp_st%fhmlm ,  &  
       microp_st%hmpim ,  &  
       microp_st%accslm,  &  
       microp_st%dlfm  ,  &  
       microp_st%autoln,  &  
       microp_st%accrln,  &  
       microp_st%bergnn,  &  
       microp_st%fhtimn,  &  
       microp_st%fhtctn,  &  
       microp_st%fhmln ,  &  
       microp_st%accsln,  &  
       microp_st%activn,  &  
       microp_st%dlfn  ,  &  
       microp_st%autoim,  &  
       microp_st%accsim,  &  
       microp_st%difm  ,  &  
       microp_st%nuclin,  &  
       microp_st%autoin,  &  
       microp_st%accsin,  &  
       microp_st%hmpin ,  &  
       microp_st%difn  ,  &  
       microp_st%cmel  ,  &  
       microp_st%cmei  ,  &   
       microp_st%trspcm,  &  
       microp_st%trspcn,  &  
       microp_st%trspim,  &  
       microp_st%trspin,  &  
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
       microp_st%fallrm,  & 
       microp_st%fallsm,  & 
       microp_st%fallgm,  & 
       microp_st%fallrn,  & 
       microp_st%fallsn,  & 
       microp_st%fallgn,  & 
       microp_st%fhmrm )
    else
       deallocate(dnlf, dnif, dsf, dnsf)    
    end if
end subroutine zm_conv_tend
!=========================================================================================


subroutine zm_conv_tend_2( state,  ptend,  ztodt, pbuf,mu, eu, &
     du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath, species_class) 

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use error_messages, only: alloc_err
   use physconst,      only: spec_class_aerosol, spec_class_gas 
 
! Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend          ! indivdual parameterization tendencies
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in):: mu(pcols,pver) 
   real(r8), intent(in):: eu(pcols,pver) 
   real(r8), intent(in):: du(pcols,pver) 
   real(r8), intent(in):: md(pcols,pver) 
   real(r8), intent(in):: ed(pcols,pver) 
   real(r8), intent(in):: dp(pcols,pver) 
   
   ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(in):: dsubcld(pcols) 
   
   ! wg layer thickness in mbs between lcl and maxi.    
   integer, intent(in) :: jt(pcols)   
   
   ! wg top  level index of deep cumulus convection.
   integer, intent(in) :: maxg(pcols) 
   
   ! wg gathered values of maxi.
   integer, intent(in) :: ideep(pcols)
   
   ! w holds position of gathered points vs longitude index 
   integer, intent(in)  :: lengath

   integer, intent(in) :: species_class(:)

! Local variables
   integer :: i, lchnk, istat, m 
   integer :: nstep
   real(r8), dimension(pcols,pver) :: dpdry

! physics buffer fields 
   integer ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   logical   :: lq(pcnst)

!
! Initialize
!
  lq(:) = .FALSE.
  lq(:) = .not. cnst_is_convtran1(:)
  call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=lq )

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

!
! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk

   nstep = get_nstep()
   if((convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero) then
      do m = 1, pcnst
         if ( (species_class(m) == spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m) == spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
         end if
      enddo
   endif


   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1,lengath
         dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
      end do

      call t_startf ('convtran2')
      call convtran (lchnk,                                        &
                     ptend%lq,state%q, pcnst,  mu, md,   &
                     du, eu, ed, dp, dsubcld,  &
                     jt,maxg,ideep, 1, lengath,  &
                     nstep,   fracis,  ptend%q, dpdry,  ztodt)  
      call t_stopf ('convtran2')
   end if

   if((convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero) then
      do m = 1, pcnst
         if ( (species_class(m) == spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m) == spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
            ptend%q(:,:,m) = 0.0_r8
         end if
      enddo
   endif

end subroutine zm_conv_tend_2


subroutine zm_conv_micro_outfld(microp_st, dlf, dif, dnlf, dnif, frz, lchnk, ncol)

   use cam_history,   only: outfld

   type(zm_microp_st),intent(in)  :: microp_st
   real(r8), intent(in) :: dlf(:,:)    ! detrainment of convective cloud liquid water mixing ratio.
   real(r8), intent(in) :: dif(:,:)    ! detrainment of convective cloud ice mixing ratio. 
   real(r8), intent(in) :: dnlf(:,:)   ! detrainment of convective cloud liquid water num concen.
   real(r8), intent(in) :: dnif(:,:)   ! detrainment of convective cloud ice num concen.
   real(r8), intent(in) :: frz(:,:)    ! heating rate due to freezing 
   integer, intent(in)         :: lchnk
   integer, intent(in)         :: ncol

   integer :: i,k

   real(r8) :: cice_snum(pcols,pver)      ! convective cloud ice sample number.
   real(r8) :: cliq_snum(pcols,pver)      ! convective cloud liquid sample number.
   real(r8) :: crain_snum(pcols,pver)     ! convective rain water sample number.
   real(r8) :: csnow_snum(pcols,pver)     ! convective snow sample number.
   real(r8) :: cgraupel_snum(pcols,pver)  ! convective graupel sample number.
   real(r8) :: wu_snum(pcols,pver)        ! vertical velocity sample number

       do k = 1,pver
          do i = 1,ncol
             if (microp_st%qice(i,k) .gt. 0.0_r8) then
                cice_snum(i,k) = 1.0_r8
             else
                cice_snum(i,k) = 0.0_r8
             end if
             if (microp_st%qliq(i,k) .gt. 0.0_r8) then
                cliq_snum(i,k) = 1.0_r8
             else
                cliq_snum(i,k) = 0.0_r8
             end if
             if (microp_st%qsnow(i,k) .gt. 0.0_r8) then
                csnow_snum(i,k) = 1.0_r8
             else
                csnow_snum(i,k) = 0.0_r8
             end if
             if (microp_st%qrain(i,k) .gt. 0.0_r8) then
                crain_snum(i,k) = 1.0_r8
             else
                crain_snum(i,k) = 0.0_r8
             end if
             if (microp_st%qgraupel(i,k) .gt. 0.0_r8) then
                cgraupel_snum(i,k) = 1.0_r8
             else
                cgraupel_snum(i,k) = 0.0_r8
             end if
             if (microp_st%wu(i,k) .gt. 0.0_r8) then
                wu_snum(i,k) = 1.0_r8
             else
                wu_snum(i,k) = 0.0_r8
             end if

          end do
       end do

      call outfld('CLIQSNUM',cliq_snum      ,pcols, lchnk)
      call outfld('CICESNUM',cice_snum      ,pcols, lchnk)
      call outfld('CRAINNUM',crain_snum     ,pcols, lchnk)
      call outfld('CSNOWNUM',csnow_snum     ,pcols, lchnk)
      call outfld('CGRAPNUM',cgraupel_snum  ,pcols, lchnk)
      call outfld('WUZMSNUM',wu_snum        ,pcols, lchnk)

      call outfld('DIFZM'   ,dif            ,pcols, lchnk)
      call outfld('DLFZM'   ,dlf            ,pcols, lchnk)
      call outfld('DNIFZM'  ,dnif           ,pcols, lchnk)
      call outfld('DNLFZM'  ,dnlf           ,pcols, lchnk)
      call outfld('FRZZM'   ,frz            ,pcols, lchnk)

      call outfld('WUZM'    ,microp_st%wu   ,pcols, lchnk)
!      call outfld('FRZZM'   ,microp_st%frz  ,pcols, lchnk)

      call outfld('CLDLIQZM',microp_st%qliq ,pcols, lchnk)
      call outfld('CLDICEZM',microp_st%qice ,pcols, lchnk)
      call outfld('QRAINZM' ,microp_st%qrain          ,pcols, lchnk)
      call outfld('QSNOWZM' ,microp_st%qsnow          ,pcols, lchnk)
      call outfld('QGRAPZM' ,microp_st%qgraupel       ,pcols, lchnk)
      
      call outfld('QNLZM'   ,microp_st%qnl            ,pcols, lchnk)
      call outfld('QNIZM'   ,microp_st%qni            ,pcols, lchnk)
      call outfld('QNRZM'   ,microp_st%qnr            ,pcols, lchnk)
      call outfld('QNSZM'   ,microp_st%qns            ,pcols, lchnk)
      call outfld('QNGZM'   ,microp_st%qng            ,pcols, lchnk)

      call outfld('AUTOL_M' ,microp_st%autolm         ,pcols, lchnk)
      call outfld('ACCRL_M' ,microp_st%accrlm         ,pcols, lchnk)
      call outfld('BERGN_M' ,microp_st%bergnm         ,pcols, lchnk)
      call outfld('FHTIM_M' ,microp_st%fhtimm         ,pcols, lchnk)
      call outfld('FHTCT_M' ,microp_st%fhtctm         ,pcols, lchnk)
      call outfld('FHML_M'  ,microp_st%fhmlm          ,pcols, lchnk)
      call outfld('HMPI_M'  ,microp_st%hmpim          ,pcols, lchnk)
      call outfld('ACCSL_M' ,microp_st%accslm         ,pcols, lchnk)
      call outfld('DLF_M'   ,microp_st%dlfm           ,pcols, lchnk)  

      call outfld('AUTOL_N' ,microp_st%autoln         ,pcols, lchnk)
      call outfld('ACCRL_N' ,microp_st%accrln         ,pcols, lchnk)
      call outfld('BERGN_N' ,microp_st%bergnn         ,pcols, lchnk)
      call outfld('FHTIM_N' ,microp_st%fhtimn         ,pcols, lchnk)
      call outfld('FHTCT_N' ,microp_st%fhtctn         ,pcols, lchnk)
      call outfld('FHML_N'  ,microp_st%fhmln          ,pcols, lchnk)
      call outfld('ACCSL_N' ,microp_st%accsln         ,pcols, lchnk)
      call outfld('ACTIV_N' ,microp_st%activn         ,pcols, lchnk)
      call outfld('DLF_N'   ,microp_st%dlfn           ,pcols, lchnk)
      call outfld('AUTOI_M' ,microp_st%autoim         ,pcols, lchnk)
      call outfld('ACCSI_M' ,microp_st%accsim         ,pcols, lchnk)
      call outfld('DIF_M'   ,microp_st%difm           ,pcols, lchnk)
      call outfld('NUCLI_N' ,microp_st%nuclin         ,pcols, lchnk)
      call outfld('AUTOI_N' ,microp_st%autoin         ,pcols, lchnk)
      call outfld('ACCSI_N' ,microp_st%accsin         ,pcols, lchnk)
      call outfld('HMPI_N'  ,microp_st%hmpin          ,pcols, lchnk)
      call outfld('DIF_N'   ,microp_st%difn           ,pcols, lchnk)
      call outfld('COND_M'  ,microp_st%cmel           ,pcols, lchnk)
      call outfld('DEPOS_M' ,microp_st%cmei           ,pcols, lchnk)

      call outfld('TRSPC_M' ,microp_st%trspcm         ,pcols, lchnk)
      call outfld('TRSPC_N' ,microp_st%trspcn         ,pcols, lchnk)
      call outfld('TRSPI_M' ,microp_st%trspim         ,pcols, lchnk)
      call outfld('TRSPI_N' ,microp_st%trspin         ,pcols, lchnk)

      call outfld('ACCGR_M' ,microp_st%accgrm         ,pcols, lchnk)
      call outfld('ACCGL_M' ,microp_st%accglm         ,pcols, lchnk)
      call outfld('ACCGSL_M',microp_st%accgslm       ,pcols, lchnk)
      call outfld('ACCGSR_M',microp_st%accgsrm       ,pcols, lchnk)
      call outfld('ACCGIR_M',microp_st%accgirm       ,pcols, lchnk)
      call outfld('ACCGRI_M',microp_st%accgrim       ,pcols, lchnk)
      call outfld('ACCGRS_M',microp_st%accgrsm       ,pcols, lchnk)

      call outfld('ACCGSL_N',microp_st%accgsln       ,pcols, lchnk)
      call outfld('ACCGSR_N',microp_st%accgsrn       ,pcols, lchnk)
      call outfld('ACCGIR_N',microp_st%accgirn       ,pcols, lchnk)

      call outfld('ACCSRI_M',microp_st%accsrim       ,pcols, lchnk)
      call outfld('ACCIGL_M',microp_st%acciglm       ,pcols, lchnk)
      call outfld('ACCIGR_M',microp_st%accigrm       ,pcols, lchnk)
      call outfld('ACCSIR_M',microp_st%accsirm       ,pcols, lchnk)

      call outfld('ACCIGL_N',microp_st%accigln       ,pcols, lchnk)
      call outfld('ACCIGR_N',microp_st%accigrn       ,pcols, lchnk)
      call outfld('ACCSIR_N',microp_st%accsirn       ,pcols, lchnk)
      call outfld('ACCGL_N' ,microp_st%accgln        ,pcols, lchnk)
      call outfld('ACCGR_N' ,microp_st%accgrn        ,pcols, lchnk)

      call outfld('ACCIL_M' ,microp_st%accilm        ,pcols, lchnk)
      call outfld('ACCIL_N' ,microp_st%acciln        ,pcols, lchnk)

      call outfld('FALLR_M' ,microp_st%fallrm        ,pcols, lchnk)
      call outfld('FALLS_M' ,microp_st%fallsm        ,pcols, lchnk)
      call outfld('FALLG_M' ,microp_st%fallgm        ,pcols, lchnk)
      call outfld('FALLR_N' ,microp_st%fallrn        ,pcols, lchnk)
      call outfld('FALLS_N' ,microp_st%fallsn        ,pcols, lchnk)
      call outfld('FALLG_N' ,microp_st%fallgn        ,pcols, lchnk)

      call outfld('FHMR_M'  ,microp_st%fhmrm         ,pcols, lchnk)

 end subroutine zm_conv_micro_outfld
!=========================================================================================



end module zm_conv_intr
