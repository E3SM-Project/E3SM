! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_vars
 use mpas_kind_types
 
 implicit none
 public
 save


!>\brief contains all local variables and arrays used in the physics parameterizations.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> add-ons and modifications:
!> --------------------------
!>    * added the variables sf_surface_physics,alswvisdir_p,alswvisdif_p,alswnirdir_p,alswnirdif_p,
!>      swvisdir_p,swvisdif_p,swnirdir_p,and swnirdif_p to upgrade the RRTMG short wave radiation
!>      code to WRF version 3.4.1. see definition of each individual variables below.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-03-11.
!>    * removed call to the updated Kain-Fritsch convection scheme.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-05-29.
!>    * added the arrays o3clim_p for implementation of monthly-varying climatological ozone in the
!>      long wave and short wave RRTMG radiation codes.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-07-08.
!>    * corrected definition of local variable dx_p.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2013-08-23.
!>    * renamed local variable conv_deep_scheme to convection_scheme.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-09-18.


!==================================================================================================
!list of physics parameterizations:
!==================================================================================================

 character(len=StrKIND),public:: microp_scheme
 character(len=StrKIND),public:: convection_scheme
 character(len=StrKIND),public:: gwdo_scheme
 character(len=StrKIND),public:: lsm_scheme
 character(len=StrKIND),public:: pbl_scheme
 character(len=StrKIND),public:: radt_cld_scheme
 character(len=StrKIND),public:: radt_lw_scheme
 character(len=StrKIND),public:: radt_sw_scheme
 character(len=StrKIND),public:: sfclayer_scheme

!==================================================================================================
!wrf-variables:these variables are needed to keep calls to different physics parameterizations
!as in wrf model.
!==================================================================================================

 logical:: l_radtlw                   !controls call to longwave radiation parameterization.
 logical:: l_radtsw                   !controls call to shortwave radiation parameterization.
 logical:: l_conv                     !controls call to convective parameterization.
 logical:: l_camlw                    !controls when to save local CAM LW abs and ems arrays.
 logical:: l_diags                    !controls when to calculate physics diagnostics.
 logical:: l_acrain                   !when .true., limit to accumulated rain is applied.
 logical:: l_acradt                   !when .true., limit to lw and sw radiation is applied.

 integer,public:: ids,ide,jds,jde,kds,kde
 integer,public:: ims,ime,jms,jme,kms,kme
 integer,public:: its,ite,jts,jte,kts,kte
 integer,public:: iall
 integer,public:: n_microp

 integer,public:: num_months          !number of months                                         [-]

 real(kind=RKIND),public:: dt_dyn     !time-step for dynamics
 real(kind=RKIND),public:: dt_microp  !time-step for cloud microphysics parameterization.
 real(kind=RKIND),public:: dt_radtlw  !time-step for longwave radiation parameterization      [mns]
 real(kind=RKIND),public:: dt_radtsw  !time-step for shortwave radiation parameterization     [mns]
 
 real(kind=RKIND),public:: xice_threshold

 real(kind=RKIND),dimension(:,:),allocatable:: &
    area_p             !grid cell area                                                         [m2]

!... arrays related to surface:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    psfc_p,           &!surface pressure                                                       [Pa]
    ptop_p             !model-top pressure                                                     [Pa]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    fzm_p,            &!weight for interpolation to w points                                    [-]
    fzp_p              !weight for interpolation to w points                                    [-]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
!... arrays related to u- and v-velocities interpolated to theta points:
    u_p,              &!u-velocity interpolated to theta points                               [m/s]
    v_p                !v-velocity interpolated to theta points                               [m/s]
    
!... arrays related to vertical sounding:
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    zz_p,             &!
    pres_p,           &!pressure                                                               [Pa]
    pi_p,             &!(p_phy/p0)**(r_d/cp)                                                    [-]
    z_p,              &!height of layer                                                         [m]
    zmid_p,           &!height of middle of layer                                               [m]
    dz_p,             &!layer thickness                                                         [m]
    t_p,              &!temperature                                                             [K]
    th_p,             &!potential temperature                                                   [K]
    al_p,             &!inverse of air density                                              [m3/kg]
    rho_p,            &!air density                                                         [kg/m3]
    rh_p               !relative humidity                                                       [-]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    qv_p,             &!water vapor mixing ratio                                            [kg/kg]
    qc_p,             &!cloud water mixing ratio                                            [kg/kg]
    qr_p,             &!rain mixing ratio                                                   [kg/kg]
    qi_p,             &!cloud ice mixing ratio                                              [kg/kg]
    qs_p,             &!snow mixing ratio                                                   [kg/kg]
    qg_p               !graupel mixing ratio                                                [kg/kg]
 
!... arrays located at w (vertical velocity) points, or at interface between layers:
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    w_p,              &!vertical velocity                                                     [m/s]
    pres2_p,          &!pressure                                                               [Pa]
    t2_p               !temperature                                                             [K]

!... arrays used for calculating the hydrostatic pressure and exner function:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    psfc_hyd_p,       &!surface pressure                                                       [Pa]
    psfc_hydd_p        !"dry" surface pressure                                                 [Pa]
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    pres_hyd_p,       &!pressure located at theta levels                                       [Pa]
    pres_hydd_p,      &!"dry" pressure located at theta levels                                 [Pa]
    pres2_hyd_p,      &!pressure located at w-velocity levels                                  [Pa]
    pres2_hydd_p,     &!"dry" pressure located at w-velocity levels                            [Pa]
    znu_hyd_p          !(pres_hyd_p / P0) needed in the Tiedtke convection scheme              [Pa]

!=============================================================================================
!... variables related to ozone climatlogy:
!=============================================================================================

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    o3clim_p           !climatological ozone volume mixing ratio                         [???]

!==================================================================================================
!... variables and arrays related to parameterization of cloud microphysics:
!    warm_phase: logical that determines if we want to run warm-phase cloud microphysics only.
!          If set to false, cold-phase cloud microphysics is active. In MPAS, we always assume
!          that the ice phase is included (except for the Kessler scheme which includes water
!          clouds only.

!    f_qv,f_qc,f_qr,f_qi,f_qs,f_qg: These logicals were initially defined in WRF to determine
!          which kind of hydrometeors are present. Here, we assume that all six water species
!          are present, even if their mixing ratios and number concentrations are zero.

!==================================================================================================

 logical,parameter:: &
    warm_rain=.false.  !warm-phase cloud microphysics only (used in WRF).

 logical,parameter:: &
    f_qv = .true.,    &!
    f_qc = .true.,    &!
    f_qr = .true.,    &!
    f_qi = .true.,    &!
    f_qs = .true.,    &!
    f_qg = .true.      !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    f_ice,            &!fraction of cloud ice (used in WRF only).
    f_rain             !fraction of rain (used in WRF only).

 real(kind=RKIND),dimension(:,:),allocatable:: &
    rainnc_p,         &!
    rainncv_p,        &!
    snownc_p,         &!
    snowncv_p,        &!
    graupelnc_p,      &!
    graupelncv_p,     &!
    sr_p

!==================================================================================================
!... variables and arrays related to parameterization of convection:
!==================================================================================================
 integer,public:: n_cu
 real(kind=RKIND),public:: dt_cu

 logical,dimension(:,:),allocatable:: &
 	cu_act_flag
 real(kind=RKIND),dimension(:,:),allocatable::   &
    rainc_p,          &!
    raincv_p,         &!
    pratec_p           !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    rthcuten_p,       &!
    rqvcuten_p,       &!
    rqccuten_p,       &!
    rqicuten_p         !

!... kain fritsch specific arrays:
 real(kind=RKIND),dimension(:,:),allocatable::   &
    cubot_p,          &!lowest convective level                                                 [-]
    cutop_p,          &!highest convective level                                                [-]
    nca_p              !counter for cloud relaxation time                                       [-]
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    w0avg_p          !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    rqrcuten_p,       &!
    rqscuten_p         !

!... tiedtke specific arrays:
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    znu_p,            &!
    rqvdynten_p,      &!
    rqvdynblten_p      !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    rucuten_p,        &!
    rvcuten_p          !

!==================================================================================================
!... variables and arrays related to parameterization of pbl:
!==================================================================================================

 integer,dimension(:,:),allocatable:: &
    kpbl_p             !index of PBL top                                                        [-]

 real(kind=RKIND),public:: dt_pbl

 real(kind=RKIND),dimension(:,:),allocatable:: &
    ctopo_p,          &!correction to topography                                                [-]
    ctopo2_p,         &!correction to topography 2                                              [-]                  
    hpbl_p             !PBL height                                                              [m]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    exch_p             !exchange coefficient                                                    [-]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    rublten_p,        &!
    rvblten_p,        &!
    rthblten_p,       &!
    rqvblten_p,       &!
    rqcblten_p,       &!
    rqiblten_p         !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    kzh_p,            &!
    kzm_p,            &!
    kzq_p              !

!==================================================================================================
!... variables and arrays related to parameterization of gravity wave drag over orography:
!==================================================================================================

 real(kind=RKIND),dimension(:,:),allocatable:: &
    var2d_p,          &!orographic variance                                                    [m2]
    con_p,            &!orographic convexity                                                   [m2]
    oa1_p,            &!orographic direction asymmetry function                                 [-]
    oa2_p,            &!orographic direction asymmetry function                                 [-]
    oa3_p,            &!orographic direction asymmetry function                                 [-]
    oa4_p,            &!orographic direction asymmetry function                                 [-]
    ol1_p,            &!orographic direction asymmetry function                                 [-]
    ol2_p,            &!orographic direction asymmetry function                                 [-]
    ol3_p,            &!orographic direction asymmetry function                                 [-]
    ol4_p              !orographic direction asymmetry function                                 [-]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    dx_p               !mean distance between cell centers                                      [m]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    dusfcg_p,         &!vertically-integrated gwdo u-stress                              [Pa m s-1]
    dvsfcg_p           !vertically-integrated gwdo v -stress                             [Pa m s-1]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    dtaux3d_p,        &!gravity wave drag over orography u-stress                           [m s-1]
    dtauy3d_p          !gravity wave drag over orography u-stress                           [m s-1]

!==================================================================================================
!... variables and arrays related to parameterization of surface layer:
!==================================================================================================
 real(kind=RKIND),dimension(:,:),allocatable:: &
    br_p,             &!bulk richardson number                                                  [-]
    cd_p,             &!momentum exchange coeff at 10 meters                                    [?]
    cda_p,            &!momentum exchange coeff at the lowest model level                       [?]
    cpm_p,            &!
    chs_p,            &!
    chs2_p,           &!
    ck_p,             &!enthalpy exchange coeff at 10 meters                                    [?]
    cka_p,            &!enthalpy exchange coeff at the lowest model level                       [?]
    cqs2_p,           &!
    gz1oz0_p,         &!log of z1 over z0                                                       [-]
    flhc_p,           &!exchange coefficient for heat                                           [-]
    flqc_p,           &!exchange coefficient for moisture                                       [-]
    hfx_p,            &!upward heat flux at the surface                                      [W/m2]
    fh_p,             &!integrated stability function for heat                                  [-]
    fm_p,             &!integrated stability function for momentum                              [-]             
    lh_p,             &!latent heat flux at the surface                                      [W/m2]
    mavail_p,         &!surface moisture availability                                           [-]
    mol_p,            &!T* in similarity theory                                                 [K]
    pblh_p,           &!PBL height                                                              [m]
    psih_p,           &!similarity theory for heat                                              [-]
    psim_p,           &!similarity theory for momentum                                          [-]
    q2_p,             &!specific humidity at 2m                                             [kg/kg]
    qfx_p,            &!upward moisture flux at the surface                               [kg/m2/s]
    qgh_p,            &!
    qsfc_p,           &!specific humidity at lower boundary                                 [kg/kg]
    regime_p,         &!flag indicating PBL regime (stable_p,unstable_p,etc...)                 [-]
    rmol_p,           &!1 / Monin Ob length                                                     [-]
    t2m_p,            &!temperature at 2m                                                       [K]
    th2m_p,           &!potential temperature at 2m                                             [K]
    u10_p,            &!u at 10 m                                                             [m/s]
    ust_p,            &!u* in similarity theory                                               [m/s]
    ustm_p,           &!u* in similarity theory       without vconv correction                [m/s]
    v10_p,            &!v at 10 m                                                             [m/s]
    wspd_p,           &!wind speed                                                            [m/s]
    znt_p,            &!time-varying roughness length                                           [m]
    zol_p              !


!==================================================================================================
!... variables and arrays related to parameterization of short-wave radiation:
!==================================================================================================

 real(kind=RKIND):: &
    declin,           &!solar declination                                                       [-]
    solcon             !solar constant                                                      [W m-2]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    coszr_p,          &!cosine of the solar zenith angle                                        [-]
    gsw_p,            &!net shortwave flux at surface                                       [W m-2]
    swcf_p,           &!shortwave cloud forcing at top-of-atmosphere                        [W m-2]
    swdnb_p,          &!all-sky downwelling shortwave flux at bottom-of-atmosphere          [J m-2]
    swdnbc_p,         &!clear-sky downwelling shortwave flux at bottom-of-atmosphere        [J m-2]
    swdnt_p,          &!all-sky downwelling shortwave flux at top-of-atmosphere             [J m-2]
    swdntc_p,         &!clear-sky downwelling shortwave flux at top-of-atmosphere           [J m-2]
    swupb_p,          &!all-sky upwelling shortwave flux at bottom-of-atmosphere            [J m-2]
    swupbc_p,         &!clear-sky upwelling shortwave flux at bottom-of-atmosphere          [J m-2]
    swupt_p,          &!all-sky upwelling shortwave flux at top-of-atmosphere               [J m-2]
    swuptc_p           !clear-sky upwelling shortwave flux at top-of-atmosphere             [J m-2]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    swvisdir_p,       &!visible direct downward flux                                        [W m-2]
    swvisdif_p,       &!visible diffuse downward flux                                       [W m-2]
    swnirdir_p,       &!near-IR direct downward flux                                        [W m-2]
    swnirdif_p         !near-IR diffuse downward flux                                       [W m-2]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    swdnflx_p,        &!
    swdnflxc_p,       &!
    swupflx_p,        &!
    swupflxc_p         !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    rthratensw_p       !uncoupled theta tendency due to shortwave radiation                 [K s-1]

!==================================================================================================
!... variables and arrays related to parameterization of long-wave radiation:
!==================================================================================================

 integer,dimension(:,:),allocatable:: &
    nlrad_p            !number of layers added above the model top                              [-]
 real(kind=RKIND),dimension(:,:),allocatable:: &
    plrad_p            !pressure at model_top                                                  [Pa]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    glw_p,            &!net longwave flux at surface                                        [W m-2]
    lwcf_p,           &!longwave cloud forcing at top-of-atmosphere                         [W m-2]
    lwdnb_p,          &!all-sky downwelling longwave flux at bottom-of-atmosphere           [J m-2]
    lwdnbc_p,         &!clear-sky downwelling longwave flux at bottom-of-atmosphere         [J m-2]
    lwdnt_p,          &!all-sky downwelling longwave flux at top-of-atmosphere              [J m-2]
    lwdntc_p,         &!clear-sky downwelling longwave flux at top-of-atmosphere            [J m-2]
    lwupb_p,          &!all-sky upwelling longwave flux at bottom-of-atmosphere             [J m-2]
    lwupbc_p,         &!clear-sky upwelling longwave flux at bottom-of-atmosphere           [J m-2]
    lwupt_p,          &!all-sky upwelling longwave flux at top-of-atmosphere                [J m-2]
    lwuptc_p,         &!clear-sky upwelling longwave flux at top-of-atmosphere              [J m-2]
    olrtoa_p           !outgoing longwave radiation at top-of-the-atmosphere                [W m-2]

  real(kind=RKIND),dimension(:,:,:),allocatable:: &
    lwdnflx_p,        &!
    lwdnflxc_p,       &!
    lwupflx_p,        &!
    lwupflxc_p         !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    rthratenlw_p       !uncoupled theta tendency due to longwave radiation                  [K s-1]

!==================================================================================================
!... variables and arrays related to parameterization of long- and short-wave radiation needed
!    only by the "CAM" radiation codes:
!==================================================================================================

 logical:: doabsems

 integer:: cam_abs_dim1
 integer:: cam_abs_dim2
 integer:: num_moist
 integer:: num_aerosols
 integer:: num_aerlevels
 integer:: num_oznlevels

 real(kind=RKIND),dimension(:),allocatable:: &
    pin_p,            &!pressure levels for ozone concentration                                [Pa]
    m_hybi_p           !hybrid levels for aerosols                                             [-]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    m_psn_p,          &!
    m_psp_p            !

 real(kind=RKIND),dimension(:,:,:,:),allocatable:: &
    aerosolcn_p,      &!
    aerosolcp_p        !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    emstot_p,         &!total emissivity                                                        [-]
    cemiss_p,         &!cloud emissivity for ISCCP                                              [-]
    taucldc_p,        &!cloud water optical depth for ISCCP                                     [-]
    taucldi_p          !cloud ice optical depth for ISCCP                                       [-]    

 real(kind=RKIND),dimension(:,:,:,:),allocatable:: &
    abstot_p,         &!total layer absorptivity                                                [-]
    absnxt_p,         &!total nearest layer absorptivity                                        [-]
    ozmixm_p           !ozone mixing ratio.

!==================================================================================================
!.. variables and arrays related to cloudiness:
!==================================================================================================

 integer,parameter:: &
    icloud= 1          !used in WRF only.

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    cldfrac_p          !cloud fraction                                                          [-]

!==================================================================================================
!.. variables and arrays related to land-surface parameterization:
!==================================================================================================

 integer,public:: &
    sf_surface_physics !used to define the land surface scheme by a number instead of name. It
                       !is only needed in module_ra_rrtmg_sw.F to define the spectral surface
                       !albedos as functions of the land surface scheme.

 integer,public:: &
    num_soils          !number of soil layers                                                   [-]
    
 integer,dimension(:,:),allocatable:: &
    isltyp_p,         &!dominant soil type category                                             [-]
    ivgtyp_p           !dominant vegetation category                                            [-]

 real(kind=RKIND),dimension(:),allocatable:: &
    dzs_p              !thickness of soil layers                                                [m]
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    smcrel_p,         &!soil moisture threshold below which transpiration starts to stress      [-]
    sh2o_p,           &!unfrozen soil moisture content                        [volumetric fraction]
    smois_p,          &!soil moisture                                         [volumetric fraction]
    tslb_p             !soil temperature                                                        [K]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    acsnom_p,         &!accumulated melted snow                                            [kg m-2]
    acsnow_p,         &!accumulated snow                                                   [kg m-2]
    canwat_p,         &!canopy water                                                       [kg m-2]
    chklowq_p,        &!surface saturation flag                                                 [-]
    grdflx_p,         &!ground heat flux                                                    [W m-2]
    lai_p,            &!leaf area index                                                         [-]
    noahres_p,        &!residual of the noah land-surface scheme energy budget              [W m-2]
    potevp_p,         &!potential evaporation                                               [W m-2]
    qz0_p,            &!specific humidity at znt                                          [kg kg-1]
    rainbl_p,         &!
    sfcrunoff_p,      &!surface runoff                                                      [m s-1]
    shdmin_p,         &!minimum areal fractional coverage of annual green vegetation            [-]
    shdmax_p,         &!maximum areal fractional coverage of annual green vegetation            [-]
    smstav_p,         &!moisture availability                                                   [-]
    smstot_p,         &!total moisture                                                     [m3 m-3]
    snopcx_p,         &!snow phase change heat flux                                         [W m-2]
    snotime_p,        &!
    snowc_p,          &!snow water equivalent                                              [kg m-2]
    snowh_p,          &!physical snow depth                                                     [m]
    swdown_p,         &!downward shortwave flux at the surface                              [W m-2]
    udrunoff_p,       &!sub-surface runoff                                                  [m s-1]
    tmn_p,            &!soil temperature at lower boundary                                      [K]
    vegfra_p,         &!vegetation fraction                                                     [-]
    z0_p               !background roughness length                                             [m]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    alswvisdir_p,     &!direct-beam surface albedo in visible spectrum                          [-]
    alswvisdif_p,     &!diffuse-beam surface albedo in visible spectrum                         [-]
    alswnirdir_p,     &!direct-beam surface albedo in near-IR spectrum                          [-]
    alswnirdif_p       !diffuse-beam surface albedo in near-IR spectrum                         [-]

!==================================================================================================
!.. variables and arrays related to surface characteristics:
!==================================================================================================

 real(kind=RKIND),dimension(:,:),allocatable:: &
    xlat_p,           &!longitude, west is negative                                       [degrees]
    xlon_p             !latitude, south is negative                                       [degrees]

 real(kind=RKIND),dimension(:,:),allocatable:: &
    sfc_albedo_p,     &!surface albedo                                                          [-]
    sfc_albbck_p,     &!surface background albedo                                               [-]
    sfc_emibck_p,     &!land surface background emissivity                                      [-]
    sfc_emiss_p,      &!land surface emissivity                                                 [-]
    snoalb_p,         &!annual max snow albedo                                                  [-]
    snow_p,           &!snow water equivalent                                              [kg m-2]
    tsk_p,            &!surface-skin temperature                                                [K]
    xice_p,           &!ice mask                                                                [-]
    xland_p            !land mask    (1 for land; 2 for water)                                  [-]

!==================================================================================================
 end module mpas_atmphys_vars
!==================================================================================================
