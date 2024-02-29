! COSP interface to unwrap derived types for interoperability with C++
module cosp_c2f
  use iso_c_binding
  use cosp_kinds,          only: wp
  use mod_cosp_config,     only: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE,SR_BINS,    &
                                 N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,               &
                                 cloudsat_DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,         &
                                 CFODD_NDBZE,      CFODD_NICOD,                           &
                                 CFODD_BNDRE,      CFODD_NCLASS,                          &
                                 CFODD_DBZE_MIN,   CFODD_DBZE_MAX,                        &
                                 CFODD_ICOD_MIN,   CFODD_ICOD_MAX,                        &
                                 CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH,                      &
                                 WR_NREGIME,                                              &
                                 numMODISTauBins,numMODISPresBins,                        &
                                 numMODISReffIceBins,numMODISReffLiqBins,                 &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins,         &
                                 ntau,modis_histTau,tau_binBounds,                        &
                                 modis_histTauEdges,tau_binEdges,                         &
                                 modis_histTauCenters,tau_binCenters,ntauV1p4,            &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4, tau_binCentersV1p4,  &
                                 grLidar532_histBsct,atlid_histBsct,vgrid_zu,vgrid_zl,    &
                                 Nlvgrid, vgrid_z,cloudsat_preclvl
  use cosp_phys_constants, only: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
  use mod_quickbeam_optics,only: size_distribution,hydro_class_init,quickbeam_optics,     &
                                 quickbeam_optics_init,gases
  use quickbeam,           only: radar_cfg
  use mod_cosp,            only: cosp_init,cosp_optical_inputs,cosp_column_inputs,        &
                                 cosp_outputs,cosp_cleanup,cosp_simulator
  use mod_rng,             only: rng_state, init_rng
  use mod_scops,           only: scops
  use mod_prec_scops,      only: prec_scops
  use mod_cosp_utils,      only: cosp_precip_mxratio
  use cosp_optics,         only: cosp_simulator_optics,lidar_optics,modis_optics,         &
                                 modis_optics_partition
  use mod_cosp_stats,      only: cosp_change_vertical_grid

  implicit none

  public :: cosp_c2f_init, cosp_c2f_run, cosp_c2f_final

  ! Local variables; control what runs and what does not
  logical :: &
       lsingle     = .false.,  & ! True if using MMF_v3_single_moment CLOUDSAT microphysical scheme (default)
       ldouble     = .true., & ! True if using MMF_v3.5_two_moment CLOUDSAT microphysical scheme
       lisccp      = .true. , & ! Local on/off switch for simulators (used by initialization)
       lmodis      = .false., & !
       lmisr       = .false., & !
       lcalipso    = .false., & !
       lgrLidar532 = .false., & !
       latlid      = .false., & !
       lcloudsat   = .false., & !
       lrttov      = .false., & !
       lparasol    = .false.    !

  ! Logicals to control output
  ! Hard-code logicals for now; these need to be consistent with EAMxx outputs anyways
  logical :: &
         Lpctisccp           = .false., & ! ISCCP mean cloud top pressure
         Lclisccp            = .true. , & ! ISCCP cloud area fraction
         Lboxptopisccp       = .false., & ! ISCCP CTP in each column
         Lboxtauisccp        = .false., & ! ISCCP optical epth in each column
         Ltauisccp           = .false., & ! ISCCP mean optical depth
         Lcltisccp           = .true. , & ! ISCCP total cloud fraction
         Lmeantbisccp        = .false., & ! ISCCP mean all-sky 10.5micron brightness temperature
         Lmeantbclrisccp     = .false., & ! ISCCP mean clear-sky 10.5micron brightness temperature
         Lalbisccp           = .false., & ! ISCCP mean cloud albedo         
         LclMISR             = .false., & ! MISR cloud fraction
         Lcltmodis           = .false., & ! MODIS total cloud fraction
         Lclwmodis           = .false., & ! MODIS liquid cloud fraction
         Lclimodis           = .false., & ! MODIS ice cloud fraction
         Lclhmodis           = .false., & ! MODIS high-level cloud fraction
         Lclmmodis           = .false., & ! MODIS mid-level cloud fraction
         Lcllmodis           = .false., & ! MODIS low-level cloud fraction
         Ltautmodis          = .false., & ! MODIS total cloud optical thicknes
         Ltauwmodis          = .false., & ! MODIS liquid optical thickness
         Ltauimodis          = .false., & ! MODIS ice optical thickness
         Ltautlogmodis       = .false., & ! MODIS total cloud optical thickness (log10 mean)
         Ltauwlogmodis       = .false., & ! MODIS liquid optical thickness (log10 mean)
         Ltauilogmodis       = .false., & ! MODIS ice optical thickness (log10 mean)
         Lreffclwmodis       = .false., & ! MODIS liquid cloud particle size
         Lreffclimodis       = .false., & ! MODIS ice particle size
         Lpctmodis           = .false., & ! MODIS cloud top pressure
         Llwpmodis           = .false., & ! MODIS cloud liquid water path
         Liwpmodis           = .false., & ! MODIS cloud ice water path
         Lclmodis            = .false., & ! MODIS cloud area fraction
         Latb532             = .false., & ! CALIPSO attenuated total backscatter (532nm)
         Latb532gr           = .false., & ! GROUND LIDAR @ 532NM attenuated total backscatter (532nm)
         Latb355             = .false., & ! ATLID attenuated total backscatter (355nm)
         LlidarBetaMol532    = .false., & ! CALIPSO molecular backscatter (532nm)         
         LlidarBetaMol532gr  = .false., & ! GROUND LIDAR @ 532NM molecular backscatter (532nm)
         LlidarBetaMol355    = .false., & ! ATLID molecular backscatter (355nm) 
         LcfadLidarsr532     = .false., & ! CALIPSO scattering ratio CFAD
         LcfadLidarsr532gr   = .false., & ! GROUND LIDAR @ 532NM scattering ratio CFAD  
         LcfadLidarsr355     = .false., & ! ATLID scattering ratio CFAD 
         Lclcalipso2         = .false., & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso          = .false., & ! CALIPSO cloud area fraction
         LclgrLidar532       = .false., & ! GROUND LIDAR @ 532NM cloud area fraction 
         Lclatlid            = .false., & ! ATLID cloud area fraction 
         Lclhcalipso         = .false., & ! CALIPSO high-level cloud fraction
         Lcllcalipso         = .false., & ! CALIPSO low-level cloud fraction
         Lclmcalipso         = .false., & ! CALIPSO mid-level cloud fraction
         Lcltcalipso         = .false., & ! CALIPSO total cloud fraction
         LclhgrLidar532      = .false., & ! GROUND LIDAR @ 532NM high-level cloud fraction 
         LcllgrLidar532      = .false., & ! GROUND LIDAR @ 532NM low-level cloud fraction 
         LclmgrLidar532      = .false., & ! GROUND LIDAR @ 532NM mid-level cloud fraction
         LcltgrLidar532      = .false., & ! GROUND LIDAR @ 532NM total cloud fraction
         Lclhatlid           = .false., & ! ATLID high-level cloud fraction
         Lcllatlid           = .false., & ! ATLID low-level cloud fraction  
         Lclmatlid           = .false., & ! ATLID mid-level cloud fraction 
         Lcltatlid           = .false., & ! ATLID total cloud fraction
         Lcltlidarradar      = .false., & ! CALIPSO-CLOUDSAT total cloud fraction
         Lcloudsat_tcc       = .false., & !
         Lcloudsat_tcc2      = .false., & !
         Lclcalipsoliq       = .false., & ! CALIPSO liquid cloud area fraction
         Lclcalipsoice       = .false., & ! CALIPSO ice cloud area fraction 
         Lclcalipsoun        = .false., & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmp       = .false., & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmpliq    = .false., & ! CALIPSO liquid cloud area fraction
         Lclcalipsotmpice    = .false., & ! CALIPSO ice cloud area fraction
         Lclcalipsotmpun     = .false., & ! CALIPSO undetected cloud area fraction
         Lcltcalipsoliq      = .false., & ! CALIPSO liquid total cloud fraction
         Lcltcalipsoice      = .false., & ! CALIPSO ice total cloud fraction
         Lcltcalipsoun       = .false., & ! CALIPSO undetected total cloud fraction
         Lclhcalipsoliq      = .false., & ! CALIPSO high-level liquid cloud fraction
         Lclhcalipsoice      = .false., & ! CALIPSO high-level ice cloud fraction
         Lclhcalipsoun       = .false., & ! CALIPSO high-level undetected cloud fraction
         Lclmcalipsoliq      = .false., & ! CALIPSO mid-level liquid cloud fraction
         Lclmcalipsoice      = .false., & ! CALIPSO mid-level ice cloud fraction
         Lclmcalipsoun       = .false., & ! CALIPSO mid-level undetected cloud fraction
         Lcllcalipsoliq      = .false., & ! CALIPSO low-level liquid cloud fraction
         Lcllcalipsoice      = .false., & ! CALIPSO low-level ice cloud fraction
         Lcllcalipsoun       = .false., & ! CALIPSO low-level undetected cloud fraction
         Lclopaquecalipso    = .false., & ! CALIPSO opaque cloud cover (2D Map)
         Lclthincalipso      = .false., & ! CALIPSO thin cloud cover (2D Map)
         Lclzopaquecalipso   = .false., & ! CALIPSO z_opaque altitude (opaque clouds only, 2D Map)
         Lclcalipsoopaque    = .false., & ! CALIPSO opaque cloud profiles 3D fraction 
         Lclcalipsothin      = .false., & ! CALIPSO thin cloud profiles 3D fraction 
         Lclcalipsozopaque   = .false., & ! CALIPSO z_opaque 3D fraction 
         Lclcalipsoopacity   = .false., & ! CALIPSO opacity 3D fraction 
         Lclopaquetemp       = .false., & ! CALIPSO opaque cloud temperature 
         Lclthintemp         = .false., & ! CALIPSO thin cloud temperature
         Lclzopaquetemp      = .false., & ! CALIPSO z_opaque temperature  
         Lclopaquemeanz      = .false., & ! CALIPSO opaque cloud altitude  
         Lclthinmeanz        = .false., & ! CALIPSO thin cloud altitude 
         Lclthinemis         = .false., & ! CALIPSO thin cloud emissivity
         Lclopaquemeanzse    = .false., & ! CALIPSO opaque cloud altitude with respect to SE 
         Lclthinmeanzse      = .false., & ! CALIPSO thin cloud altitude with respect to SE
         Lclzopaquecalipsose = .false., & ! CALIPSO z_opaque altitude with respect to SE
         LcfadDbze94         = .false., & ! CLOUDSAT radar reflectivity CFAD
         Ldbze94             = .false., & ! CLOUDSAT radar reflectivity
         LparasolRefl        = .false., & ! PARASOL reflectance
         Ltbrttov            = .false., & ! RTTOV mean clear-sky brightness temperature
         Lptradarflag0       = .false., & ! CLOUDSAT 
         Lptradarflag1       = .false., & ! CLOUDSAT 
         Lptradarflag2       = .false., & ! CLOUDSAT 
         Lptradarflag3       = .false., & ! CLOUDSAT 
         Lptradarflag4       = .false., & ! CLOUDSAT 
         Lptradarflag5       = .false., & ! CLOUDSAT 
         Lptradarflag6       = .false., & ! CLOUDSAT 
         Lptradarflag7       = .false., & ! CLOUDSAT 
         Lptradarflag8       = .false., & ! CLOUDSAT 
         Lptradarflag9       = .false., & ! CLOUDSAT 
         Lradarpia           = .false., & ! CLOUDSAT 
         Lwr_occfreq         = .false., & ! CloudSat+MODIS joint diagnostics
         Lcfodd              = .false.    ! CloudSat+MODIS joint diagnostics
 
  ! Input namelist fields (hard-code these)
  integer, parameter :: rttov_Nchannels = 3
  real(wp), dimension(rttov_Nchannels) :: rttov_surfem = (/0.0, 0.0, 0.0/)
  integer ::                          & !
       !Nlvgrid = 40,                  & ! Number of vertical levels for statistical outputs (USE_VGRID=.true.)
       surface_radar = 0,             & ! surface=1/spaceborne=0
       cloudsat_use_gas_abs = 1,      & ! Include gaseous absorption (1=yes/0=no)
       cloudsat_do_ray = 0,           & ! Calculate output Rayleigh (1=yes/0=no)
       lidar_ice_type = 0,            & ! Ice particle shape in lidar calculations (0=ice-spheres/1=ice-non-spherical)
       overlap = 3,                   & ! Overlap type: 1=max, 2=rand, 3=max/rand
       isccp_topheight = 1,           & ! ISCCP cloud top height
       isccp_topheight_direction = 2, & ! ISCCP cloud top height direction
       rttov_platform = 1,            & ! RTTOV: Satellite platform
       rttov_satellite = 15,          & ! RTTOV: Satellite
       rttov_instrument = 5,          & ! RTTOV: Instrument
       rttov_channels(rttov_Nchannels) = (/1, 2, 3/)  ! RTTOV: Number of channels to be computed
  real(wp) ::                         & !
       cloudsat_radar_freq = 94.0,    & ! CloudSat radar frequency (GHz)
       cloudsat_k2 = -1,              & ! |K|^2, -1=use frequency dependent default
       rttov_ZenAng = 50.0,           & ! RTTOV: Satellite Zenith Angle
       co2 = 5.241e-04,               & ! CO2 mixing ratio
       ch4 = 9.139e-07,               & ! CH4 mixing ratio
       n2o = 4.665e-07,               & ! n2o mixing ratio
       co  = 2.098e-07                  ! co mixing ratio
  logical ::                          & !
       use_vgrid = .true.,            & ! Use fixed vertical grid for outputs?
       csat_vgrid = .true.              ! CloudSat vertical grid?

  ! These only need to be allocated once, so we let them persist as module data
  type(size_distribution) :: sd
  type(radar_cfg) :: rcfg_cloudsat
  type(cosp_outputs) :: cospOUT
  type(cosp_optical_inputs) :: cospIN
  type(cosp_column_inputs) :: cospstateIn


  ! Indices to address arrays of LS and CONV hydrometeors
  integer,parameter :: &
       I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
       I_LSCICE = 2, & ! Large-scale (stratiform) ice
       I_LSRAIN = 3, & ! Large-scale (stratiform) rain
       I_LSSNOW = 4, & ! Large-scale (stratiform) snow
       I_CVCLIQ = 5, & ! Convective liquid
       I_CVCICE = 6, & ! Convective ice
       I_CVRAIN = 7, & ! Convective rain
       I_CVSNOW = 8, & ! Convective snow
       I_LSGRPL = 9    ! Large-scale (stratiform) groupel

  ! Stratiform and convective clouds in frac_out (scops output).
  integer, parameter :: &
       I_LSC = 1, & ! Large-scale clouds
       I_CVC = 2    ! Convective clouds

  ! Microphysical settings for the precipitation flux to mixing ratio conversion
  real(wp),parameter,dimension(N_HYDRO) :: &
                 ! LSL   LSI      LSR       LSS   CVL  CVI      CVR       CVS       LSG
       N_ax    = (/-1., -1.,     8.e6,     3.e6, -1., -1.,     8.e6,     3.e6,     4.e6/),&
       N_bx    = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
       alpha_x = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
       c_x     = (/-1., -1.,    842.0,     4.84, -1., -1.,    842.0,     4.84,     94.5/),&
       d_x     = (/-1., -1.,      0.8,     0.25, -1., -1.,      0.8,     0.25,      0.5/),&
       g_x     = (/-1., -1.,      0.5,      0.5, -1., -1.,      0.5,      0.5,      0.5/),&
       a_x     = (/-1., -1.,    524.0,    52.36, -1., -1.,    524.0,    52.36,   209.44/),&
       b_x     = (/-1., -1.,      3.0,      3.0, -1., -1.,      3.0,      3.0,      3.0/),&
       gamma_1 = (/-1., -1., 17.83725, 8.284701, -1., -1., 17.83725, 8.284701, 11.63230/),&
       gamma_2 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/),&
       gamma_3 = (/-1., -1.,      2.0,      2.0, -1., -1.,      2.0,      2.0,      2.0/),&
       gamma_4 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/)

  character(len=64) :: cloudsat_micro_scheme = 'MMF_v3.5_two_moment'

contains

  subroutine cosp_c2f_init(npoints, ncolumns, nlevels) bind(c, name='cosp_c2f_init')
    integer(kind=c_int), value, intent(in) :: npoints, ncolumns, nlevels
    ! Initialize/allocate COSP input and output derived types
    nlvgrid = 40
    call construct_cospIN(npoints,ncolumns,nlevels,cospIN)
    call construct_cospstatein(npoints,nlevels,rttov_nchannels,cospstateIN)
    call construct_cosp_outputs(npoints, ncolumns, nlevels, nlvgrid, rttov_nchannels, cospOUT)

    ! Initialize quickbeam_optics, also if two-moment radar microphysics scheme is wanted...
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
       ldouble = .true.
       lsingle = .false.
    endif

    if (Lcloudsat) then
       call quickbeam_optics_init()

       ! Initialize the distributional parameters for hydrometeors in radar simulator
       call hydro_class_init(lsingle,ldouble,sd)
    end if

    ! Initialize COSP simulator
    call cosp_init(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, LgrLidar532, Latlid,        &
         Lparasol, Lrttov,                                                                 &
         cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,                           &
         cloudsat_do_ray, isccp_topheight, isccp_topheight_direction, surface_radar,       &
         rcfg_cloudsat, use_vgrid, csat_vgrid, Nlvgrid, Nlevels, cloudsat_micro_scheme)
  end subroutine cosp_c2f_init 

  subroutine cosp_c2f_run(npoints, ncolumns, nlevels, ntau, nctp, &
       emsfc_lw, sunlit, skt, T_mid, p_mid, p_int, qv, &
       cldfrac, reff_qc, reff_qi, dtau067, dtau105, isccp_cldtot, isccp_ctptau &
       ) bind(C, name='cosp_c2f_run')
    integer(kind=c_int), value, intent(in) :: npoints, ncolumns, nlevels, ntau, nctp
    real(kind=c_double), value, intent(in) :: emsfc_lw
    real(kind=c_double), intent(in), dimension(npoints) :: sunlit, skt
    real(kind=c_double), intent(in), dimension(npoints,nlevels) :: T_mid, p_mid, qv, cldfrac, reff_qc, reff_qi, dtau067, dtau105
    real(kind=c_double), intent(in), dimension(npoints,nlevels+1) :: p_int
    real(kind=c_double), intent(inout), dimension(npoints) :: isccp_cldtot
    real(kind=c_double), intent(inout), dimension(npoints,ntau,nctp) :: isccp_ctptau
    ! Takes normal arrays as input and populates COSP derived types
    character(len=256),dimension(100) :: cosp_status
    integer :: nptsperit
    integer :: start_idx
    integer :: end_idx

    ! Locals for subsample_and_optics
    ! TODO: trim this down
    real(wp), dimension(nPoints,nLevels) :: tca,cca,mr_lsliq,mr_lsice,mr_ccliq,   &
         mr_ccice,dtau_c,dtau_s,dem_c,dem_s,mr_lsrain,mr_lssnow,mr_lsgrpl,mr_ccrain,&
         mr_ccsnow
    real(wp), dimension(nPoints,nLevels,N_HYDRO) :: reff
 
    nptsperit = npoints

    ! In-cloud values are assumed. If ncolumns = 1, then convert in-cloud values to gridbox
    if (ncolumns == 1) then
       tca(:npoints,:nlevels) = cldfrac(:npoints,:nlevels)
       cca(:npoints,:nlevels) = 0
       mr_lsliq(:npoints,:nlevels) = 0
       mr_ccliq(:npoints,:nlevels) = 0
       mr_lsice(:npoints,:nlevels) = 0
       mr_ccice(:npoints,:nlevels) = 0
       dtau_c(:npoints,:nlevels) = 0
       dtau_s(:npoints,:nlevels) = cldfrac(:npoints,:nlevels) * dtau067(:npoints,:nlevels)
       dem_c (:npoints,:nlevels) = 0
       dem_s (:npoints,:nlevels) = 1._wp - exp(-cldfrac(:npoints,:nlevels) * dtau105(:npoints,:nlevels))
       mr_lsrain(:npoints,:nlevels) = 0
       mr_ccrain(:npoints,:nlevels) = 0
       mr_lssnow(:npoints,:nlevels) = 0
       mr_lssnow(:npoints,:nlevels) = 0
       mr_ccsnow(:npoints,:nlevels) = 0
       mr_lsgrpl(:npoints,:nlevels) = 0
       reff = 0  ! FIXME
    else
       tca(:npoints,:nlevels) = cldfrac(:npoints,:nlevels)
       cca(:npoints,:nlevels) = 0
       mr_lsliq(:npoints,:nlevels) = 0
       mr_ccliq(:npoints,:nlevels) = 0
       mr_lsice(:npoints,:nlevels) = 0
       mr_ccice(:npoints,:nlevels) = 0
       dtau_c(:npoints,:nlevels) = 0
       dtau_s(:npoints,:nlevels) = dtau067(:npoints,:nlevels)
       dem_c (:npoints,:nlevels) = 0
       dem_s (:npoints,:nlevels) = 1._wp - exp(-dtau105(:npoints,:nlevels))
       mr_lsrain(:npoints,:nlevels) = 0
       mr_ccrain(:npoints,:nlevels) = 0
       mr_lssnow(:npoints,:nlevels) = 0
       mr_lssnow(:npoints,:nlevels) = 0
       mr_ccsnow(:npoints,:nlevels) = 0
       mr_lsgrpl(:npoints,:nlevels) = 0
       reff = 0  ! FIXME
    end if

    start_idx = 1
    end_idx = npoints
    ! Translate arrays to derived types
    cospIN%emsfc_lw         = emsfc_lw
    cospIN%rcfg_cloudsat    = rcfg_cloudsat
!   cospstateIN%hgt_matrix  = zlev(start_idx:end_idx,Nlevels:1:-1) ! km
    cospstateIN%sunlit      = sunlit(start_idx:end_idx)            ! 0-1
    cospstateIN%skt         = skt(start_idx:end_idx)               ! K
!   cospstateIN%surfelev    = surfelev(start_idx:end_idx)          ! m
!   cospstateIN%land        = landmask(start_idx:end_idx)          ! 0-1 (*note* model specific)
    cospstateIN%qv          = qv(start_idx:end_idx,1:Nlevels)   ! kg/kg
    cospstateIN%at          = T_mid(start_idx:end_idx,1:Nlevels) !Nlevels:1:-1)    ! K
    cospstateIN%pfull       = p_mid(start_idx:end_idx,1:Nlevels) !Nlevels:1:-1)    ! Pa
    ! Pressure at interface (nlevels+1). Set uppermost interface to 0.
    !cospstateIN%phalf(:,2:Nlevels+1) = p_int(start_idx:end_idx,Nlevels:1:-1)   ! Pa
    cospstateIN%phalf(:,1:Nlevels+1) = p_int(start_idx:end_idx,1:Nlevels+1)   ! Pa
!   ! Height of bottom interfaces of model layers (nlevels).
!   ! cospstateIN%hgt_matrix_half(:,1) contains the bottom of the top layer.
!   ! cospstateIN%hgt_matrix_half(:,Nlevels) contains the bottom of the surface layer.
!   cospstateIN%hgt_matrix_half(:,1:Nlevels) = zlev_half(start_idx:end_idx,Nlevels:1:-1) ! km

    ! Generate subcolumns and compute optical inputs.
    call subsample_and_optics(nPoints, nLevels, nColumns, N_HYDRO, overlap, use_vgrid, &
       lidar_ice_type, sd, tca, cca, mr_lsrain, mr_lssnow,  &
       mr_lsgrpl, mr_ccrain, mr_ccsnow, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,       &
       reff, dtau_c, dtau_s, dem_c, dem_s, cospstateIN, cospIN)

    ! Call cosp
    cosp_status = cosp_simulator(cospIN, cospstateIN, cospOUT, start_idx, end_idx, .false.)

    ! Translate derived types to output arrays
    isccp_cldtot(:npoints) = cospOUT%isccp_totalcldarea(:npoints)
    isccp_ctptau(:npoints,:,:) = cospOUT%isccp_fq(:npoints,:,:)

  end subroutine cosp_c2f_run

  subroutine cosp_c2f_final() bind(C, name='cosp_c2f_final')
    call destroy_cospIN(cospIN)
    call destroy_cospstateIN(cospstateIN)
    call destroy_cosp_outputs(cospOUT)
  end subroutine cosp_c2f_final

  ! These are mostly copied from cosp_test.f90, and are pretty obnoxious
  ! TODO: clean this up!

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(npoints,ncolumns,nlevels,y)
    ! Inputs
    integer,intent(in) :: &
         npoints,  & ! Number of horizontal gridpoints
         ncolumns, & ! Number of subcolumns
         nlevels     ! Number of vertical levels
    ! Outputs 
    type(cosp_optical_inputs),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Npart    = 4
    y%Nrefl    = PARASOL_NREFL
    allocate(y%frac_out(npoints,       ncolumns,nlevels))

    if (Lmodis .or. Lmisr .or. Lisccp) then
       allocate(y%tau_067(npoints,        ncolumns,nlevels),&
                y%emiss_11(npoints,       ncolumns,nlevels))
    endif
    if (Lcalipso) then
       allocate(y%betatot_calipso(npoints,        ncolumns,nlevels),&
                y%betatot_ice_calipso(npoints,    ncolumns,nlevels),&
                y%betatot_liq_calipso(npoints,    ncolumns,nlevels),&
                y%tautot_calipso(npoints,         ncolumns,nlevels),&
                y%tautot_ice_calipso(npoints,     ncolumns,nlevels),&
                y%tautot_liq_calipso(npoints,     ncolumns,nlevels),&
                y%beta_mol_calipso(npoints,                nlevels),&
                y%tau_mol_calipso(npoints,                 nlevels),&
                y%tautot_S_ice(npoints,   ncolumns        ),&
                y%tautot_S_liq(npoints,   ncolumns        ))
    endif

    if (LgrLidar532) then
       allocate(y%beta_mol_grLidar532(npoints,          nlevels),& 
                y%betatot_grLidar532(npoints,  ncolumns,nlevels),& 
                y%tau_mol_grLidar532(npoints,           nlevels),& 
                y%tautot_grLidar532(npoints,   ncolumns,nlevels)) 
    endif

    if (Latlid) then
       allocate(y%beta_mol_atlid(npoints,             nlevels),& 
                y%betatot_atlid(npoints,     ncolumns,nlevels),& 
                y%tau_mol_atlid(npoints,              nlevels),& 
                y%tautot_atlid(npoints,      ncolumns,nlevels))
    endif 

    if (Lcloudsat) then
       allocate(y%z_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%kr_vol_cloudsat(npoints, ncolumns,nlevels),&
                y%g_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%fracPrecipIce(npoints,   ncolumns))
    endif
    if (Lmodis) then
       allocate(y%fracLiq(npoints,        ncolumns,nlevels),&
                y%asym(npoints,           ncolumns,nlevels),&
                y%ss_alb(npoints,         ncolumns,nlevels))
    endif
    

  end subroutine construct_cospIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  subroutine construct_cospstateIN(npoints,nlevels,nchan,y)
    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal gridpoints
         nlevels, & ! Number of vertical levels
         nchan      ! Number of channels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y         
    
    allocate(y%sunlit(npoints),y%skt(npoints),y%land(npoints),y%at(npoints,nlevels),     &
             y%pfull(npoints,nlevels),y%phalf(npoints,nlevels+1),y%qv(npoints,nlevels),  &
             y%o3(npoints,nlevels),y%hgt_matrix(npoints,nlevels),y%u_sfc(npoints),       &
             y%v_sfc(npoints),y%lat(npoints),y%lon(nPoints),y%emis_sfc(nchan),           &
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),&
             y%fl_snow(nPoints,nLevels),y%fl_rain(nPoints,nLevels),y%seaice(npoints),    &
             y%tca(nPoints,nLevels),y%hgt_matrix_half(npoints,nlevels))

  end subroutine construct_cospstateIN

  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  ! ######################################################################################  
  ! TODO: This is WAY too many dummy arguments! These can just be defined at module scope I think
  ! and then initialized once at init
  subroutine construct_cosp_outputs(Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
    ! Inputs
    integer,intent(in) :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
          Nlvgrid,         & ! Number of levels in L3 stats computation
          Nchan              ! Number of RTTOV channels  
          
    ! Outputs
    type(cosp_outputs),intent(out) :: x  ! COSP output structure  
   
    ! ISCCP simulator outputs
    if (Lboxtauisccp)    allocate(x%isccp_boxtau(Npoints,Ncolumns)) 
    if (Lboxptopisccp)   allocate(x%isccp_boxptop(Npoints,Ncolumns))
    if (Lclisccp)        allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
    if (Lcltisccp)       allocate(x%isccp_totalcldarea(Npoints))
    if (Lpctisccp)       allocate(x%isccp_meanptop(Npoints))
    if (Ltauisccp)       allocate(x%isccp_meantaucld(Npoints))
    if (Lmeantbisccp)    allocate(x%isccp_meantb(Npoints))
    if (Lmeantbclrisccp) allocate(x%isccp_meantbclr(Npoints))
    if (Lalbisccp)       allocate(x%isccp_meanalbedocld(Npoints))

    ! MISR simulator
    if (LclMISR) then 
       allocate(x%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins))
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.
       allocate(x%misr_dist_model_layertops(Npoints,numMISRHgtBins))
       allocate(x%misr_meanztop(Npoints))
       allocate(x%misr_cldarea(Npoints))    
    endif
    
    ! MODIS simulator
    if (Lcltmodis)     allocate(x%modis_Cloud_Fraction_Total_Mean(Npoints))
    if (Lclwmodis)     allocate(x%modis_Cloud_Fraction_Water_Mean(Npoints))
    if (Lclimodis)     allocate(x%modis_Cloud_Fraction_Ice_Mean(Npoints))
    if (Lclhmodis)     allocate(x%modis_Cloud_Fraction_High_Mean(Npoints))
    if (Lclmmodis)     allocate(x%modis_Cloud_Fraction_Mid_Mean(Npoints))
    if (Lcllmodis)     allocate(x%modis_Cloud_Fraction_Low_Mean(Npoints))
    if (Ltautmodis)    allocate(x%modis_Optical_Thickness_Total_Mean(Npoints))
    if (Ltauwmodis)    allocate(x%modis_Optical_Thickness_Water_Mean(Npoints))
    if (Ltauimodis)    allocate(x%modis_Optical_Thickness_Ice_Mean(Npoints))
    if (Ltautlogmodis) allocate(x%modis_Optical_Thickness_Total_LogMean(Npoints))
    if (Ltauwlogmodis) allocate(x%modis_Optical_Thickness_Water_LogMean(Npoints))
    if (Ltauilogmodis) allocate(x%modis_Optical_Thickness_Ice_LogMean(Npoints))
    if (Lreffclwmodis) allocate(x%modis_Cloud_Particle_Size_Water_Mean(Npoints))
    if (Lreffclimodis) allocate(x%modis_Cloud_Particle_Size_Ice_Mean(Npoints))
    if (Lpctmodis)     allocate(x%modis_Cloud_Top_Pressure_Total_Mean(Npoints))
    if (Llwpmodis)     allocate(x%modis_Liquid_Water_Path_Mean(Npoints))
    if (Liwpmodis)     allocate(x%modis_Ice_Water_Path_Mean(Npoints))
    if (Lclmodis) then
        allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins,numMODISPresBins))
        allocate(x%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins,numMODISReffLiqBins))   
        allocate(x%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins,numMODISReffIceBins))
    endif
    
    ! LIDAR simulator
    if (LlidarBetaMol532) allocate(x%calipso_beta_mol(Npoints,Nlevels))
    if (Latb532)          allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532)  then
        allocate(x%calipso_srbval(SR_BINS+1))
        allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid))
        allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))  
    endif
    if (Lclcalipso)       allocate(x%calipso_lidarcld(Npoints,Nlvgrid))
    if (Lclhcalipso .or. Lclmcalipso .or. Lcllcalipso .or. Lcltcalipso) then
        allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))
    endif   
    if (Lclcalipsoice .or. Lclcalipsoliq .or. Lclcalipsoun) then
        allocate(x%calipso_lidarcldphase(Npoints,Nlvgrid,6))
    endif
    if (Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsoice .or. Lclcalipsotmpun .or. Lclcalipsotmpice) then
        allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
    endif
    if (Lcllcalipsoice .or. Lclmcalipsoice .or. Lclhcalipsoice .or.                   &
        Lcltcalipsoice .or. Lcllcalipsoliq .or. Lclmcalipsoliq .or.                   &
        Lclhcalipsoliq .or. Lcltcalipsoliq .or. Lcllcalipsoun  .or.                   &
        Lclmcalipsoun  .or. Lclhcalipsoun  .or. Lcltcalipsoun) then
        allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))     
    endif
    if (Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso) then
        allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
    endif 
    if (Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp) then 
        allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))  
    endif
    if (Lclopaquemeanz .or. Lclthinmeanz) then 
        allocate(x%calipso_cldtypemeanz(Npoints,2))
    endif 
    if (Lclopaquemeanzse .or. Lclthinmeanzse .or. Lclzopaquecalipsose) then 
        allocate(x%calipso_cldtypemeanzse(Npoints,3)) 
    endif 
    if (Lclthinemis) then 
        allocate(x%calipso_cldthinemis(Npoints))
    endif
    if (Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or. Lclcalipsoopacity) then 
        allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))
    endif
    ! These 2 outputs are part of the calipso output type, but are not controlled by an 
    ! logical switch in the output namelist, so if all other fields are on, then allocate
    if (LlidarBetaMol532 .or. Latb532        .or. LcfadLidarsr532 .or. Lclcalipso  .or.  &
        Lclcalipsoice    .or. Lclcalipsoliq  .or. Lclcalipsoun    .or. Lclcalipso2 .or.  &
        Lclhcalipso      .or. Lclmcalipso    .or. Lcllcalipso     .or. Lcltcalipso .or.  &
        Lclcalipsotmp    .or. Lclcalipsoice  .or. Lclcalipsotmpun .or.                   &
        Lclcalipsotmpliq .or. Lcllcalipsoice .or. Lclmcalipsoice  .or.                   &
        Lclhcalipsoice   .or. Lcltcalipsoice .or. Lcllcalipsoliq  .or.                   &
        Lclmcalipsoliq   .or. Lclhcalipsoliq .or. Lcltcalipsoliq  .or.                   &
        Lcllcalipsoun    .or. Lclmcalipsoun  .or. Lclhcalipsoun   .or. Lcltcalipsoun) then
       allocate(x%calipso_tau_tot(Npoints,Ncolumns,Nlevels))       
       allocate(x%calipso_temp_tot(Npoints,Nlevels))               
    endif

    ! GROUND LIDAR @ 532NM simulator
    if (LlidarBetaMol532gr) allocate(x%grLidar532_beta_mol(Npoints,Nlevels))
    if (Latb532gr)          allocate(x%grLidar532_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532gr) then 
        allocate(x%grLidar532_srbval(SR_BINS+1)) 
        allocate(x%grLidar532_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif
    if (LclgrLidar532)     allocate(x%grLidar532_lidarcld(Npoints,Nlvgrid)) 
    if (LclhgrLidar532 .or. LclmgrLidar532 .or. LcllgrLidar532 .or. LcltgrLidar532) then
        allocate(x%grLidar532_cldlayer(Npoints,LIDAR_NCAT))  
    endif
      
    ! ATLID simulator
    if (LlidarBetaMol355) allocate(x%atlid_beta_mol(Npoints,Nlevels))
    if (Latb355)          allocate(x%atlid_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr355) then
        allocate(x%atlid_srbval(SR_BINS+1)) 
        allocate(x%atlid_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif 
    if (Lclatlid)     allocate(x%atlid_lidarcld(Npoints,Nlvgrid)) 
    if (Lclhatlid .or. Lclmatlid .or. Lcllatlid .or. Lcltatlid) then
        allocate(x%atlid_cldlayer(Npoints,LIDAR_NCAT)) 
    endif
      
    ! PARASOL
    if (Lparasolrefl) then
        allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
        allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif 

    ! Cloudsat simulator
    if (Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,cloudsat_DBZE_BINS,Nlvgrid))
    if (Lptradarflag0 .or. Lptradarflag1 .or. Lptradarflag2 .or. Lptradarflag3 .or. &
        Lptradarflag4 .or. Lptradarflag5 .or. Lptradarflag6 .or. Lptradarflag7 .or. &
        Lptradarflag8 .or. Lptradarflag9) then
       allocate(x%cloudsat_precip_cover(Npoints,cloudsat_DBZE_BINS))
    endif
    if (Lradarpia) allocate(x%cloudsat_pia(Npoints))

    ! Combined CALIPSO/CLOUDSAT fields
    if (Lclcalipso2)    allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
    if (Lcltlidarradar) allocate(x%radar_lidar_tcc(Npoints))
    if (Lcloudsat_tcc) allocate(x%cloudsat_tcc(Npoints))
    if (Lcloudsat_tcc2) allocate(x%cloudsat_tcc2(Npoints))
            
    ! RTTOV
    if (Ltbrttov) allocate(x%rttov_tbs(Npoints,Nchan))

    ! Joint MODIS/CloudSat Statistics
    !if (Lwr_occfreq)  allocate(x%wr_occfreq_ntotal(Npoints,WR_NREGIME))
    !if (Lcfodd)       allocate(x%cfodd_ntotal(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS))

  end subroutine construct_cosp_outputs
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y
    if (allocated(y%tau_067))             deallocate(y%tau_067)
    if (allocated(y%emiss_11))            deallocate(y%emiss_11)
    if (allocated(y%frac_out))            deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))        deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))        deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))                deallocate(y%asym)
    if (allocated(y%ss_alb))              deallocate(y%ss_alb)
    if (allocated(y%fracLiq))             deallocate(y%fracLiq)
    if (allocated(y%beta_mol_grLidar532)) deallocate(y%beta_mol_grLidar532)
    if (allocated(y%betatot_grLidar532))  deallocate(y%betatot_grLidar532)
    if (allocated(y%tau_mol_grLidar532))  deallocate(y%tau_mol_grLidar532)
    if (allocated(y%tautot_grLidar532))   deallocate(y%tautot_grLidar532)
    if (allocated(y%beta_mol_atlid))      deallocate(y%beta_mol_atlid) 
    if (allocated(y%betatot_atlid))       deallocate(y%betatot_atlid) 
    if (allocated(y%tau_mol_atlid))       deallocate(y%tau_mol_atlid) 
    if (allocated(y%tautot_atlid))        deallocate(y%tautot_atlid)
    if (allocated(y%fracPrecipIce))      deallocate(y%fracPrecipIce)
  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)    
    if (allocated(y%surfelev))        deallocate(y%surfelev)
    
  end subroutine destroy_cospstateIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cosp_outputs(y)
     type(cosp_outputs),intent(inout) :: y

     ! Deallocate and nullify
     if (associated(y%calipso_beta_mol))          then
        deallocate(y%calipso_beta_mol)
        nullify(y%calipso_beta_mol)
     endif
     if (associated(y%calipso_temp_tot))          then
        deallocate(y%calipso_temp_tot)
        nullify(y%calipso_temp_tot)     
     endif
     if (associated(y%calipso_betaperp_tot))      then
        deallocate(y%calipso_betaperp_tot)
        nullify(y%calipso_betaperp_tot)     
     endif
     if (associated(y%calipso_beta_tot))          then
        deallocate(y%calipso_beta_tot)    
        nullify(y%calipso_beta_tot)     
     endif
     if (associated(y%calipso_tau_tot))           then
        deallocate(y%calipso_tau_tot) 
        nullify(y%calipso_tau_tot)     
     endif
     if (associated(y%calipso_lidarcldphase))     then
        deallocate(y%calipso_lidarcldphase)
        nullify(y%calipso_lidarcldphase)     
     endif
     if (associated(y%calipso_lidarcldtype))     then
        deallocate(y%calipso_lidarcldtype) 
        nullify(y%calipso_lidarcldtype) 
     endif
     if (associated(y%calipso_cldlayerphase))     then
        deallocate(y%calipso_cldlayerphase)
        nullify(y%calipso_cldlayerphase)     
     endif
     if (associated(y%calipso_lidarcldtmp))       then
        deallocate(y%calipso_lidarcldtmp)
        nullify(y%calipso_lidarcldtmp)     
     endif
     if (associated(y%calipso_cldlayer))          then
        deallocate(y%calipso_cldlayer)
        nullify(y%calipso_cldlayer)     
     endif
     if (associated(y%calipso_cldtype))          then
        deallocate(y%calipso_cldtype)
        nullify(y%calipso_cldtype)  
     endif  
     if (associated(y%calipso_cldtypetemp))      then
        deallocate(y%calipso_cldtypetemp) 
        nullify(y%calipso_cldtypetemp) 
     endif 
     if (associated(y%calipso_cldtypemeanz))     then
        deallocate(y%calipso_cldtypemeanz)
        nullify(y%calipso_cldtypemeanz)  
     endif  
     if (associated(y%calipso_cldtypemeanzse))   then 
        deallocate(y%calipso_cldtypemeanzse) 
        nullify(y%calipso_cldtypemeanzse)  
     endif 
     if (associated(y%calipso_cldthinemis))      then
        deallocate(y%calipso_cldthinemis)
        nullify(y%calipso_cldthinemis)
     endif 
     if (associated(y%calipso_lidarcld))         then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)     
     endif
     if (associated(y%calipso_srbval))            then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)     
     endif
     if (associated(y%calipso_cfad_sr))          then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)     
     endif
     if (associated(y%grLidar532_beta_mol))     then
        deallocate(y%grLidar532_beta_mol)
        nullify(y%grLidar532_beta_mol)  
     endif
     if (associated(y%grLidar532_beta_tot))     then 
        deallocate(y%grLidar532_beta_tot) 
        nullify(y%grLidar532_beta_tot)
     endif 
     if (associated(y%grLidar532_cldlayer))     then 
        deallocate(y%grLidar532_cldlayer) 
        nullify(y%grLidar532_cldlayer) 
     endif 
     if (associated(y%grLidar532_lidarcld))     then 
        deallocate(y%grLidar532_lidarcld) 
        nullify(y%grLidar532_lidarcld) 
     endif 
     if (associated(y%grLidar532_cfad_sr))      then 
        deallocate(y%grLidar532_cfad_sr)  
        nullify(y%grLidar532_cfad_sr) 
     endif 
     if (associated(y%grLidar532_srbval))       then
        deallocate(y%grLidar532_srbval) 
        nullify(y%grLidar532_srbval) 
     endif  
     if (associated(y%atlid_beta_mol))           then
        deallocate(y%atlid_beta_mol) 
        nullify(y%atlid_beta_mol) 
     endif 
     if (associated(y%atlid_beta_tot))           then
        deallocate(y%atlid_beta_tot) 
        nullify(y%atlid_beta_tot) 
     endif
     if (associated(y%atlid_cldlayer))           then 
        deallocate(y%atlid_cldlayer)  
        nullify(y%atlid_cldlayer) 
     endif 
     if (associated(y%atlid_lidarcld))           then 
        deallocate(y%atlid_lidarcld)  
        nullify(y%atlid_lidarcld) 
     endif  
     if (associated(y%atlid_cfad_sr))            then
        deallocate(y%atlid_cfad_sr) 
        nullify(y%atlid_cfad_sr) 
     endif
     if (associated(y%atlid_srbval))             then 
        deallocate(y%atlid_srbval) 
        nullify(y%atlid_srbval) 
     endif
     if (associated(y%parasolPix_refl))           then
        deallocate(y%parasolPix_refl)
        nullify(y%parasolPix_refl)     
     endif
     if (associated(y%parasolGrid_refl))          then
        deallocate(y%parasolGrid_refl) 
        nullify(y%parasolGrid_refl)     
     endif
     if (associated(y%cloudsat_Ze_tot))           then
        deallocate(y%cloudsat_Ze_tot) 
        nullify(y%cloudsat_Ze_tot)  
     endif
     if (associated(y%cloudsat_cfad_ze))          then
        deallocate(y%cloudsat_cfad_ze)
        nullify(y%cloudsat_cfad_ze)     
     endif
     if (associated(y%cloudsat_precip_cover))     then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     endif
     if (associated(y%cloudsat_pia))              then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc) 
        nullify(y%cloudsat_tcc)  
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2) 
        nullify(y%cloudsat_tcc2)  
     endif
     if (associated(y%radar_lidar_tcc))           then
        deallocate(y%radar_lidar_tcc) 
        nullify(y%radar_lidar_tcc)  
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc) 
        nullify(y%cloudsat_tcc)  
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2) 
        nullify(y%cloudsat_tcc2)  
     endif
     if (associated(y%lidar_only_freq_cloud))     then
        deallocate(y%lidar_only_freq_cloud)
        nullify(y%lidar_only_freq_cloud)     
     endif
     if (associated(y%isccp_totalcldarea))        then
        deallocate(y%isccp_totalcldarea) 
        nullify(y%isccp_totalcldarea)  
     endif
     if (associated(y%isccp_meantb))              then
        deallocate(y%isccp_meantb) 
        nullify(y%isccp_meantb)     
     endif
     if (associated(y%isccp_meantbclr))           then
        deallocate(y%isccp_meantbclr)
        nullify(y%isccp_meantbclr)  
     endif
     if (associated(y%isccp_meanptop))            then
        deallocate(y%isccp_meanptop)
        nullify(y%isccp_meanptop)     
     endif
     if (associated(y%isccp_meantaucld))          then
        deallocate(y%isccp_meantaucld) 
        nullify(y%isccp_meantaucld)       
     endif
     if (associated(y%isccp_meanalbedocld))       then
        deallocate(y%isccp_meanalbedocld)
        nullify(y%isccp_meanalbedocld)     
     endif
     if (associated(y%isccp_boxtau))              then
        deallocate(y%isccp_boxtau)
        nullify(y%isccp_boxtau)       
     endif
     if (associated(y%isccp_boxptop))             then
        deallocate(y%isccp_boxptop)
        nullify(y%isccp_boxptop)     
     endif
     if (associated(y%isccp_fq))                  then
        deallocate(y%isccp_fq)
        nullify(y%isccp_fq)       
     endif
     if (associated(y%misr_fq))                   then
        deallocate(y%misr_fq) 
        nullify(y%misr_fq)     
     endif
     if (associated(y%misr_dist_model_layertops)) then
        deallocate(y%misr_dist_model_layertops)
        nullify(y%misr_dist_model_layertops)       
     endif
     if (associated(y%misr_meanztop))             then
        deallocate(y%misr_meanztop)
        nullify(y%misr_meanztop)     
     endif
     if (associated(y%misr_cldarea))              then
        deallocate(y%misr_cldarea)
        nullify(y%misr_cldarea)      
     endif
     if (associated(y%rttov_tbs))                 then
        deallocate(y%rttov_tbs)
        nullify(y%rttov_tbs)     
     endif
     if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Total_Mean)       
        nullify(y%modis_Cloud_Fraction_Total_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)     
        nullify(y%modis_Cloud_Fraction_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Water_Mean)           
        nullify(y%modis_Cloud_Fraction_Water_Mean)           
     endif
     if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
        deallocate(y%modis_Cloud_Fraction_High_Mean)     
        nullify(y%modis_Cloud_Fraction_High_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)       
        nullify(y%modis_Cloud_Fraction_Mid_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Low_Mean)     
        nullify(y%modis_Cloud_Fraction_Low_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Total_Mean)  
        nullify(y%modis_Optical_Thickness_Total_Mean)  
     endif
     if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Water_Mean)     
        nullify(y%modis_Optical_Thickness_Water_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
        deallocate(y%modis_Optical_Thickness_Ice_Mean)       
        nullify(y%modis_Optical_Thickness_Ice_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Total_LogMean)    
        nullify(y%modis_Optical_Thickness_Total_LogMean)    
     endif
     if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Water_LogMean)     
        nullify(y%modis_Optical_Thickness_Water_LogMean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)     
        nullify(y%modis_Optical_Thickness_Ice_LogMean)     
     endif
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)       
        nullify(y%modis_Cloud_Particle_Size_Water_Mean)       
     endif
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)     
        nullify(y%modis_Cloud_Particle_Size_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)           
        nullify(y%modis_Cloud_Top_Pressure_Total_Mean)           
     endif
     if (associated(y%modis_Liquid_Water_Path_Mean))                         then
        deallocate(y%modis_Liquid_Water_Path_Mean)     
        nullify(y%modis_Liquid_Water_Path_Mean)     
     endif
     if (associated(y%modis_Ice_Water_Path_Mean))                            then
        deallocate(y%modis_Ice_Water_Path_Mean)       
        nullify(y%modis_Ice_Water_Path_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
        nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
        nullify(y%modis_Optical_thickness_vs_ReffLIQ)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        nullify(y%modis_Optical_thickness_vs_ReffICE)
     endif
     !if (associated(y%cfodd_ntotal)) then
     !   deallocate(y%cfodd_ntotal)
     !   nullify(y%cfodd_ntotal)
     !endif
     !if (associated(y%wr_occfreq_ntotal)) then
     !   deallocate(y%wr_occfreq_ntotal)
     !   nullify(y%wr_occfreq_ntotal)
     !endif

  end subroutine destroy_cosp_outputs
 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  ! SUBROUTINE subsample_and_optics
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine subsample_and_optics(nPoints, nLevels, nColumns, nHydro, overlap, use_vgrid,   &
       lidar_ice_type, sd, tca, cca, fl_lsrainIN, fl_lssnowIN,    &
       fl_lsgrplIN, fl_ccrainIN, fl_ccsnowIN, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,       &
       reffIN, dtau_c, dtau_s, dem_c, dem_s, cospstateIN, cospIN)
    ! Inputs
    integer,intent(in) :: nPoints, nLevels, nColumns, nHydro, overlap, lidar_ice_type
    real(wp),intent(in),dimension(nPoints,nLevels) :: tca,cca,mr_lsliq,mr_lsice,mr_ccliq,   &
         mr_ccice,dtau_c,dtau_s,dem_c,dem_s,fl_lsrainIN,fl_lssnowIN,fl_lsgrplIN,fl_ccrainIN,&
         fl_ccsnowIN
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: reffIN
    logical,intent(in) :: use_vgrid ! .false.: outputs on model levels
                                    ! .true.:  outputs on evenly-spaced vertical levels.
    type(size_distribution),intent(inout) :: sd
    
    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN

    ! Local variables
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable :: seed
    integer,dimension(:),allocatable :: cloudsat_preclvl_index
    integer :: i,j,k
    real(wp) :: zstep
    real(wp),dimension(:,:), allocatable :: &
         ls_p_rate, cv_p_rate, frac_ls, frac_cv, prec_ls, prec_cv,g_vol
    real(wp),dimension(:,:,:),  allocatable :: &
         frac_prec, MODIS_cloudWater, MODIS_cloudIce, fracPrecipIce, fracPrecipIce_statGrid,&
         MODIS_watersize,MODIS_iceSize, MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable :: &
         mr_hydro, Reff, Np
    real(wp),dimension(nPoints,nLevels) :: &
         column_frac_out, column_prec_out, fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow
    real(wp),dimension(nPoints,nColumns,Nlvgrid) :: tempOut
    logical :: cmpGases=.true.


    if (Ncolumns .gt. 1) then
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! RNG used for subcolumn generation
       allocate(rngs(nPoints),seed(nPoints))
       seed(:)=0
       seed = int(cospstateIN%phalf(:,Nlevels+1))  ! In case of NPoints=1
       ! *NOTE* Chunking will change the seed
       if (NPoints .gt. 1) seed=int((cospstateIN%phalf(:,Nlevels+1)-minval(cospstateIN%phalf(:,Nlevels+1)))/      &
            (maxval(cospstateIN%phalf(:,Nlevels+1))-minval(cospstateIN%phalf(:,Nlevels+1)))*100000) + 1
       call init_rng(rngs, seed)
      
       ! Call scops
       call scops(NPoints,Nlevels,Ncolumns,rngs,tca,cca,overlap,cospIN%frac_out,0)
       deallocate(seed,rngs)
       
       ! Sum up precipitation rates
       allocate(ls_p_rate(nPoints,nLevels),cv_p_rate(nPoints,Nlevels))
       ls_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow) + mixing_ratio (groupel)
       cv_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow)
       
       ! Call PREC_SCOPS
       allocate(frac_prec(nPoints,nColumns,nLevels))
       call prec_scops(nPoints,nLevels,nColumns,ls_p_rate,cv_p_rate,cospIN%frac_out,frac_prec)
       deallocate(ls_p_rate,cv_p_rate)
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute fraction in each gridbox for precipitation  and cloud type.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Allocate
       allocate(frac_ls(nPoints,nLevels),prec_ls(nPoints,nLevels),                       &
                frac_cv(nPoints,nLevels),prec_cv(nPoints,nLevels))
       
       ! Initialize
       frac_ls(1:nPoints,1:nLevels) = 0._wp
       prec_ls(1:nPoints,1:nLevels) = 0._wp
       frac_cv(1:nPoints,1:nLevels) = 0._wp
       prec_cv(1:nPoints,1:nLevels) = 0._wp
       do j=1,nPoints
          do k=1,nLevels
             do i=1,nColumns
                if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
                if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 1)  prec_ls(j,k) = prec_ls(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 2)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)  prec_ls(j,k) = prec_ls(j,k)+1._wp
             enddo
             frac_ls(j,k)=frac_ls(j,k)/nColumns
             frac_cv(j,k)=frac_cv(j,k)/nColumns
             prec_ls(j,k)=prec_ls(j,k)/nColumns
             prec_cv(j,k)=prec_cv(j,k)/nColumns
          enddo
       enddo       

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Assign gridmean mixing-ratios (mr_XXXXX), effective radius (ReffIN) and number
       ! concentration (not defined) to appropriate sub-column. Here we are using scops. 
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(mr_hydro(nPoints,nColumns,nLevels,nHydro),                               &
                Reff(nPoints,nColumns,nLevels,nHydro),                                   &
                Np(nPoints,nColumns,nLevels,nHydro))

       ! Initialize
       mr_hydro(:,:,:,:) = 0._wp
       Reff(:,:,:,:)     = 0._wp
       Np(:,:,:,:)       = 0._wp
       do k=1,nColumns
          ! Subcolumn cloud fraction
          column_frac_out = cospIN%frac_out(:,k,:)
               
          ! LS clouds
          where (column_frac_out == I_LSC)
             mr_hydro(:,k,:,I_LSCLIQ) = mr_lsliq
             mr_hydro(:,k,:,I_LSCICE) = mr_lsice
             Reff(:,k,:,I_LSCLIQ)     = ReffIN(:,:,I_LSCLIQ)
             Reff(:,k,:,I_LSCICE)     = ReffIN(:,:,I_LSCICE)
          ! CONV clouds   
          elsewhere (column_frac_out == I_CVC)
             mr_hydro(:,k,:,I_CVCLIQ) = mr_ccliq
             mr_hydro(:,k,:,I_CVCICE) = mr_ccice
             Reff(:,k,:,I_CVCLIQ)     = ReffIN(:,:,I_CVCLIQ)
             Reff(:,k,:,I_CVCICE)     = ReffIN(:,:,I_CVCICE)
          end where
          
          ! Subcolumn precipitation
          column_prec_out = frac_prec(:,k,:)

          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3) )
             Reff(:,k,:,I_LSRAIN) = ReffIN(:,:,I_LSRAIN)
             Reff(:,k,:,I_LSSNOW) = ReffIN(:,:,I_LSSNOW)
             Reff(:,k,:,I_LSGRPL) = ReffIN(:,:,I_LSGRPL)
             ! CONV precipitation   
          elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
             Reff(:,k,:,I_CVRAIN) = ReffIN(:,:,I_CVRAIN)
             Reff(:,k,:,I_CVSNOW) = ReffIN(:,:,I_CVSNOW)
          end where
       enddo
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert the subcolumn mixing ratio and precipitation fluxes from gridbox mean
       ! values to fraction-based values. 
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Initialize
       fl_lsrain(:,:) = 0._wp
       fl_lssnow(:,:) = 0._wp
       fl_lsgrpl(:,:) = 0._wp
       fl_ccrain(:,:) = 0._wp
       fl_ccsnow(:,:) = 0._wp
       do k=1,nLevels
          do j=1,nPoints
             ! In-cloud mixing ratios.
             if (frac_ls(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
             endif
             if (frac_cv(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
             endif
             ! Precipitation
             if (prec_ls(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_LSRAIN) = mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                mr_hydro(j,:,k,I_LSSNOW) = mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                mr_hydro(j,:,k,I_LSGRPL) = mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
             endif
             if (prec_cv(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_CVRAIN) = mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                mr_hydro(j,:,k,I_CVSNOW) = mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
             endif
          enddo
       enddo
       deallocate(frac_ls,prec_ls,frac_cv,prec_cv)

    else
       cospIN%frac_out(:,:,:) = 1  
       allocate(mr_hydro(nPoints,1,nLevels,nHydro),Reff(nPoints,1,nLevels,nHydro),       &
                Np(nPoints,1,nLevels,nHydro))
       mr_hydro(:,1,:,I_LSCLIQ) = mr_lsliq
       mr_hydro(:,1,:,I_LSCICE) = mr_lsice
       mr_hydro(:,1,:,I_CVCLIQ) = mr_ccliq
       mr_hydro(:,1,:,I_CVCICE) = mr_ccice
       Reff(:,1,:,:)            = ReffIN
    endif
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 11 micron emissivity
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dem_c,dem_s,    &
                                  cospIN%emiss_11)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 0.67 micron optical depth
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp .or. Lmisr .or. Lmodis) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dtau_c,dtau_s,  &
                                  cospIN%tau_067)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LIDAR Polarized optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lcalipso) then
       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 532, .false.,      &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
            cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_calipso,                    &
            cospIN%betatot_calipso, cospIN%tau_mol_calipso, cospIN%tautot_calipso,         &
            cospIN%tautot_S_liq, cospIN%tautot_S_ice, cospIN%betatot_ice_calipso,          &
            cospIN%betatot_liq_calipso, cospIN%tautot_ice_calipso, cospIN%tautot_liq_calipso)
    endif

    if (LgrLidar532) then
       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 532, .true.,       &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
            cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_grLidar532,                 &
            cospIN%betatot_grLidar532, cospIN%tau_mol_grLidar532, cospIN%tautot_grLidar532)
    endif
    
    if (Latlid) then
       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 355, .false.,      &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
            cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_atlid, cospIN%betatot_atlid,&
            cospIN%tau_mol_atlid, cospIN%tautot_atlid)
    endif
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (lcloudsat) then

       ! Compute gaseous absorption (assume identical for each subcolun)
       allocate(g_vol(nPoints,nLevels))
       g_vol(:,:)=0._wp
       do i=1,nPoints
          do j=1,nLevels
             if (rcfg_cloudsat%use_gas_abs == 1 .or. (rcfg_cloudsat%use_gas_abs == 2 .and. j .eq. 1)) then
                g_vol(i,j) = gases(cospstateIN%pfull(i,j), cospstateIN%at(i,j),cospstateIN%qv(i,j),rcfg_cloudsat%freq)
             endif
             cospIN%g_vol_cloudsat(i,:,j)=g_vol(i,j)
          end do
       end do
       
       ! Loop over all subcolumns
       allocate(fracPrecipIce(nPoints,nColumns,nLevels))
       fracPrecipIce(:,:,:) = 0._wp
       do k=1,nColumns
          call quickbeam_optics(sd, rcfg_cloudsat, nPoints, nLevels, R_UNDEF,  &
               mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,&
               Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,          &
               cospstateIN%qv, cospIN%z_vol_cloudsat(1:nPoints,k,:),           &
               cospIN%kr_vol_cloudsat(1:nPoints,k,:))
          
          ! At each model level, what fraction of the precipitation is frozen?
          where(mr_hydro(:,k,:,I_LSRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_LSSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_CVRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_CVSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_LSGRPL) .gt. 0)
             fracPrecipIce(:,k,:) = (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + &
                  mr_hydro(:,k,:,I_LSGRPL)) / &
                  (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + mr_hydro(:,k,:,I_LSGRPL) + &
                  mr_hydro(:,k,:,I_LSRAIN)  + mr_hydro(:,k,:,I_CVRAIN))
          elsewhere
             fracPrecipIce(:,k,:) = 0._wp
          endwhere
       enddo
         
       ! Regrid frozen fraction to Cloudsat/Calipso statistical grid
       if (use_vgrid) then
         allocate(fracPrecipIce_statGrid(nPoints,nColumns,Nlvgrid))
         fracPrecipIce_statGrid(:,:,:) = 0._wp
         call cosp_change_vertical_grid(Npoints, Ncolumns, Nlevels, cospstateIN%hgt_matrix(:,Nlevels:1:-1), &
              cospstateIN%hgt_matrix_half(:,Nlevels:1:-1), fracPrecipIce(:,:,Nlevels:1:-1), Nlvgrid,  &
              vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1), fracPrecipIce_statGrid(:,:,Nlvgrid:1:-1))

         ! Find proper layer above de surface elevation to compute precip flags in Cloudsat/Calipso statistical grid
         allocate(cloudsat_preclvl_index(nPoints))
         cloudsat_preclvl_index(:) = 0._wp
         ! Compute the zstep distance between two atmopsheric layers
         zstep = vgrid_zl(1)-vgrid_zl(2)
         ! Computing altitude index for precip flags calculation (one layer above surfelev layer)
         cloudsat_preclvl_index(:) = cloudsat_preclvl - floor( cospstateIN%surfelev(:)/zstep )

         ! For near-surface diagnostics, we only need the frozen fraction at one layer.
         do i=1,nPoints
           cospIN%fracPrecipIce(i,:) = fracPrecipIce_statGrid(i,:,cloudsat_preclvl_index(i))
         enddo
         deallocate(cloudsat_preclvl_index)
         deallocate(fracPrecipIce_statGrid)
       endif

    endif
   
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lmodis) then
       allocate(MODIS_cloudWater(nPoints,nColumns,nLevels),                                &
                MODIS_cloudIce(nPoints,nColumns,nLevels),                                  &
                MODIS_waterSize(nPoints,nColumns,nLevels),                                 &
                MODIS_iceSize(nPoints,nColumns,nLevels),                                   &
                MODIS_opticalThicknessLiq(nPoints,nColumns,nLevels),                       &
                MODIS_opticalThicknessIce(nPoints,nColumns,nLevels))
       ! Cloud water
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater)
       ! Cloud ice
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            mr_hydro(:,:,:,I_CVCICE),mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce)  
       ! Water droplet size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            Reff(:,:,:,I_CVCLIQ),Reff(:,:,:,I_LSCLIQ),MODIS_waterSize)
       ! Ice crystal size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            Reff(:,:,:,I_CVCICE),Reff(:,:,:,I_LSCICE),MODIS_iceSize)
       
       ! Partition optical thickness into liquid and ice parts
       call modis_optics_partition(nPoints, nLevels, nColumns, MODIS_cloudWater,           &
            MODIS_cloudIce, MODIS_waterSize, MODIS_iceSize, cospIN%tau_067,                &
            MODIS_opticalThicknessLiq, MODIS_opticalThicknessIce)
       
       ! Compute assymetry parameter and single scattering albedo 
       call modis_optics(nPoints, nLevels, nColumns, MODIS_opticalThicknessLiq,            &
            MODIS_waterSize*1.0e6_wp, MODIS_opticalThicknessIce,                           &
            MODIS_iceSize*1.0e6_wp, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)
       
       ! Deallocate memory
       deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,           &
            MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce,mr_hydro,                  &
            Np,Reff)
    endif
  end subroutine subsample_and_optics
  

  
end module cosp_c2f
