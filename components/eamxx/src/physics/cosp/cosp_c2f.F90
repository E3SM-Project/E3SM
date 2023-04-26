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
       lsingle     = .true.,  & ! True if using MMF_v3_single_moment CLOUDSAT microphysical scheme (default)
       ldouble     = .false., & ! True if using MMF_v3.5_two_moment CLOUDSAT microphysical scheme
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
         Lclisccp            = .false., & ! ISCCP cloud area fraction
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
 
  ! These only need to be allocated once, so we let them persist as module data
  type(size_distribution) :: sd
  type(radar_cfg) :: rcfg_cloudsat
  type(cosp_outputs) :: cospOUT
  type(cosp_optical_inputs) :: cospIN
  type(cosp_column_inputs) :: cospstateIn

  integer, parameter :: rttov_Nchannels = 1

contains

  subroutine cosp_c2f_init(npoints, ncolumns, nlevels) bind(c, name='cosp_c3f_init')
    integer(kind=c_int), intent(in) :: npoints, ncolumns, nlevels
    ! Initialize/allocate COSP input and output derived types
    call construct_cospIN(npoints,ncolumns,nlevels,cospIN)
    call construct_cospstatein(npoints,nlevels,rttov_nchannels,cospstatein)
    call construct_cosp_outputs(npoints, ncolumns, nlevels, nlvgrid, rttov_nchannels, cospout)
  end subroutine cosp_c2f_init 

  subroutine cosp_c2f_run() bind(C, name='cosp_c2f_run')
    ! Takes normal arrays as input and populates COSP derived types
    character(len=256),dimension(100) :: cosp_status
    integer :: start_idx
    integer :: end_idx
    ! Translate arrays to derived types
    ! Call cosp
    !cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx, end_idx, .false.)
    ! Translate derived types to output arrays
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
  
end module cosp_c2f
