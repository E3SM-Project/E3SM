#define _IDX321(i, j, k, nx, ny, nz) (nx * (ny * (k - 1) + (j - 1)) + i)
module cospsimulator_intr
  ! ######################################################################################
  ! Purpose: CAM interface to
  !         Name:         CFMIP Observational Simulator Package Version 2 (COSP2)
  !         What:         Simulate ISCCP/CloudSat/CALIPSO/MISR/MODIS cloud products from 
  !                       GCM inputs
  !         Version:      v2.1.4 (August 2019)
  !         Authors:      Dustin Swales (dustin.swales@noaa.gov)
  !                       Ben Hillman (bhillma@sandia.gov)
  !
  ! Modifications:
  !
  ! ######################################################################################
  use shr_kind_mod,         only: r8 => shr_kind_r8
  use spmd_utils,           only: masterproc
  use ppgrid,               only: pcols, pver, pverp, begchunk, endchunk
  use perf_mod,             only: t_startf, t_stopf
  use cam_abortutils,       only: endrun
  use phys_control,         only: cam_physpkg_is, phys_getopts
  use cam_logfile,          only: iulog
#ifdef USE_COSP
  use quickbeam,            only: radar_cfg
  use mod_quickbeam_optics, only: size_distribution
  use mod_cosp,             only: cosp_outputs,cosp_optical_inputs,cosp_column_inputs
  use mod_cosp_config,      only: pres_binCenters, pres_binEdges, tau_binCenters,      &
       tau_binEdges, cloudsat_binCenters, cloudsat_binEdges, calipso_binCenters,       &
       calipso_binEdges, misr_histHgtCenters, misr_histHgtEdges,  PARASOL_SZA,         &
       R_UNDEF, PARASOL_NREFL, LIDAR_NCAT,SR_BINS, N_HYDRO, RTTOV_MAX_CHANNELS,        &
       numMISRHgtBins, CLOUDSAT_DBZE_BINS, LIDAR_NTEMP, calipso_histBsct,              &
       numMODISTauBins, numMODISPresBins, numMODISReffIceBins, numMODISReffLiqBins,    &
       numISCCPTauBins, numISCCPPresBins, numMISRTauBins, reffICE_binEdges,            &
       reffICE_binCenters, reffLIQ_binEdges, reffLIQ_binCenters, LIDAR_NTYPE,          &
       nCloudsatPrecipClass, &
       nsza_cosp         => PARASOL_NREFL,       &
       nprs_cosp         => npres,               &
       ntau_cosp         => ntau,                &
       ntau_cosp_modis   => ntau,                &
       nsr_cosp          => SR_BINS,             &
       nhtmisr_cosp      => numMISRHgtBins,      &
       nhydro            => N_HYDRO, &
       cloudsat_preclvl
    use mod_cosp_stats,       only: cosp_change_vertical_grid
#endif
  use crmdims, only: crm_nx_rad, crm_ny_rad, crm_nz

  implicit none
  private
  save
   
  ! Public functions/subroutines
  public :: &
       cospsimulator_intr_readnl,  &
       cospsimulator_intr_register,&
       cospsimulator_intr_init,    &
       cospsimulator_intr_run

  ! ######################################################################################
  ! Public declarations
  ! ######################################################################################
  ! Whether to do COSP calcs and I/O, default is false. If docosp is specified in 
  ! the atm_in namelist, this value is overwritten and cosp is run
  logical, public :: docosp = .false.  

  ! Frequency at which cosp is called, every cosp_nradsteps radiation timestep
  integer, public :: cosp_nradsteps = 1 ! CAM namelist variable default, not in COSP namelist
  
#ifdef USE_COSP

  ! ######################################################################################  
  ! Local declarations
  ! ######################################################################################
  interface packed_average
     module procedure packed_average1d, packed_average2d, packed_average3d
  end interface

  integer, parameter :: &
       nhtml_cosp = pver  ! Mumber of model levels is pver
  integer ::  &
       nscol_cosp,  &     ! Number of subcolumns, use namelist input Ncolumns to set.
       nht_cosp           ! Number of height for COSP radar and calipso simulator outputs.  
                          !  *set to 40 if csat_vgrid=.true., else set to Nlr*
  
  ! ######################################################################################
  ! Bin-boundaries for mixed dimensions. Calculated in cospsetupvales OR in cosp_config.F90
  ! ######################################################################################
  real(r8), target :: prsmid_cosp(nprs_cosp)            ! pressure midpoints of COSP ISCCP output
  real(r8), target :: prslim_cosp(2,nprs_cosp)
  real(r8), target :: taumid_cosp(ntau_cosp)            ! optical depth midpoints of COSP ISCCP output
  real(r8), target :: taulim_cosp(2,ntau_cosp)
  real(r8), target :: srmid_cosp(nsr_cosp)              ! sr midpoints of COSP lidar output     
  real(r8), target :: srlim_cosp(2,nsr_cosp)
  real(r8), target :: sza_cosp(nsza_cosp)
  real(r8), target :: dbzemid_cosp(CLOUDSAT_DBZE_BINS)          ! dbze midpoints of COSP radar output
  real(r8), target :: dbzelim_cosp(2,CLOUDSAT_DBZE_BINS)
  real(r8), target :: htmisrmid_cosp(nhtmisr_cosp)      ! htmisr midpoints of COSP misr simulator output
  real(r8), target :: htmisrlim_cosp(2,nhtmisr_cosp)
  real(r8), target :: taumid_cosp_modis(ntau_cosp_modis)! optical depth midpoints of COSP MODIS output
  real(r8), target :: taulim_cosp_modis(2,ntau_cosp_modis)
  real(r8), target :: reffICE_binEdges_cosp(2,numMODISReffIceBins)
  real(r8), target :: reffLIQ_binEdges_cosp(2,numMODISReffLiqBins)
  real(r8), target :: reffICE_binCenters_cosp(numMODISReffIceBins)
  real(r8), target :: reffLIQ_binCenters_cosp(numMODISReffLiqBins)

  real(r8) :: htmlmid_cosp(nhtml_cosp)                     ! Model level height midpoints for output
  real(r8),allocatable, target :: htlim_cosp(:,:)          ! height limits for COSP outputs (nht_cosp+1)
  real(r8),allocatable, target :: htmid_cosp(:)            ! height midpoints of COSP radar/lidar output (nht_cosp)
  integer, allocatable, target :: scol_cosp(:)             ! sub-column number (nscol_cosp)

  ! ######################################################################################
  ! Default namelists
  ! The CAM and COSP namelists defaults are set below.  Some of the COSP namelist 
  ! variables are part of the CAM namelist - they all begin with "cosp_" to keep their 
  ! names specific to COSP. I set their CAM namelist defaults here, not in namelist_defaults_cam.xml
  !  Variables identified as namelist variables are defined in
  !  ../models/atm/cam/bld/namelist_files/namelist_definition.xml
  ! ######################################################################################
  ! CAM
  logical :: cosp_amwg             = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_lite             = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_passive          = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_active           = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_isccp            = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_lradar_sim       = .false. ! CAM namelist variable default
  logical :: cosp_llidar_sim     = .false. ! CAM namelist variable default
  logical :: cosp_lisccp_sim       = .false. ! CAM namelist variable default
  logical :: cosp_lmisr_sim        = .false. ! CAM namelist variable default
  logical :: cosp_lmodis_sim       = .false. ! CAM namelist variable default
  logical :: cosp_histfile_aux     = .false. ! CAM namelist variable default
  logical :: cosp_lfrac_out        = .false. ! CAM namelist variable default
  integer :: cosp_ncolumns         = 50      ! CAM namelist variable default
  integer :: cosp_histfile_num     =1        ! CAM namelist variable default, not in COSP namelist 
  integer :: cosp_histfile_aux_num =-1       ! CAM namelist variable default, not in COSP namelist
  
  ! COSP
  logical :: lradar_sim       = .false.      ! COSP namelist variable, can be changed from default by CAM namelist
  logical :: llidar_sim     = .false.      ! 
  logical :: lparasol_sim     = .false.      ! 
  logical :: lgrLidar532      = .false.      !
  logical :: latlid           = .false.      !
  logical :: lisccp_sim       = .false.      ! ""
  logical :: lmisr_sim        = .false.      ! ""
  logical :: lmodis_sim       = .false.      ! ""
  logical :: lrttov_sim       = .false.      ! not running rttov, always set to .false.
  logical :: lfrac_out        = .false.      ! COSP namelist variable, can be changed from default by CAM namelist

  ! ######################################################################################  
  ! COSP parameters
  ! ######################################################################################
  ! Note: Unless otherwise specified, these are parameters that cannot be set by the CAM namelist.
  integer, parameter :: Npoints_it = 10000       ! Max # gridpoints to be processed in one iteration (10,000)
  integer :: ncolumns = 50                       ! Number of subcolumns in SCOPS (50), can be changed from default by CAM namelist
  integer :: nlr = 40                            ! Number of levels in statistical outputs 
                                                 ! (only used if USE_VGRID=.true.)  (40)
  logical :: use_vgrid = .true.                  ! Use fixed vertical grid for outputs? 
                                                 ! (if .true. then define # of levels with nlr)  (.true.)
  logical :: csat_vgrid = .true.                 ! CloudSat vertical grid? 
                                                 ! (if .true. then the CloudSat standard grid is used.
                                                 ! If set, overides use_vgrid.) (.true.)
  ! namelist variables for COSP input related to radar simulator
  real(r8) :: radar_freq = 94.0_r8               ! CloudSat radar frequency (GHz) (94.0)
  integer :: surface_radar = 0                   ! surface=1, spaceborne=0 (0)
  integer :: use_mie_tables = 0                  ! use a precomputed lookup table? yes=1,no=0 (0)
  integer :: use_gas_abs = 1                     ! include gaseous absorption? yes=1,no=0 (1)
  integer :: do_ray = 0                          ! calculate/output Rayleigh refl=1, not=0 (0)
  integer :: melt_lay = 0                        ! melting layer model off=0, on=1 (0)
  real(r8) :: k2 = -1                            ! |K|^2, -1=use frequency dependent default (-1)
  ! namelist variables for COSP input related to lidar simulator
  integer, parameter :: Nprmts_max_hydro = 12    ! Max # params for hydrometeor size distributions (12)
  integer, parameter :: Naero = 1                ! Number of aerosol species (Not used) (1)
  integer, parameter :: Nprmts_max_aero = 1      ! Max # params for aerosol size distributions (not used) (1)
  integer :: lidar_ice_type = 0                  ! Ice particle shape in lidar calculations
                                                 ! (0=ice-spheres ; 1=ice-non-spherical) (0)
  integer, parameter :: overlap = 3              ! overlap type: 1=max, 2=rand, 3=max/rand (3)

  ! namelist variables for COSP input related to ISCCP simulator
  integer :: isccp_topheight = 1                 ! 1 = adjust top height using both a computed infrared
                                                 ! brightness temperature and the visible
                                                 ! optical depth to adjust cloud top pressure.
                                                 ! Note that this calculation is most appropriate to compare
                                                 ! to ISCCP data during sunlit hours.
                                                 ! 2 = do not adjust top height, that is cloud top pressure
                                                 ! is the actual cloud top pressure in the model
                                                 ! 3 = adjust top height using only the computed infrared
                                                 ! brightness temperature. Note that this calculation is most
                                                 ! appropriate to compare to ISCCP IR only algortihm (i.e.
                                                 ! you can compare to nighttime ISCCP data with this option) (1)
  integer :: isccp_topheight_direction = 2       ! direction for finding atmosphere pressure level with
                                                 ! interpolated temperature equal to the radiance
                                                 ! determined cloud-top temperature
                                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                                 ! 1 = default setting in COSP v1.1, matches all versions of
                                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                                 ! 2 = default setting in COSP v1.3. default since V4.0 of ISCCP simulator

  ! ######################################################################################
  ! Other variables
  ! ######################################################################################  
  ! pbuf indices
  integer :: cld_idx, concld_idx, lsreffrain_idx, lsreffsnow_idx, cvreffliq_idx
  integer :: cvreffice_idx, dpcldliq_idx, dpcldice_idx
  integer :: shcldliq_idx, shcldice_idx, shcldliq1_idx, shcldice1_idx, dpflxprc_idx
  integer :: dpflxsnw_idx, shflxprc_idx, shflxsnw_idx, lsflxprc_idx, lsflxsnw_idx
  integer :: rei_idx, rel_idx
  
  ! ######################################################################################
  ! Declarations specific to COSP2
  ! ######################################################################################
  type(radar_cfg)              :: rcfg_cloudsat ! Radar configuration (Cloudsat)
  type(radar_cfg), allocatable :: rcfg_cs(:)    ! chunked version of rcfg_cloudsat
  type(size_distribution)              :: sd       ! Size distribution used by radar simulator
  type(size_distribution), allocatable :: sd_cs(:) ! chunked version of sd
  character(len=64)         :: cloudsat_micro_scheme = 'MMF_v3.5_single_moment'
  
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
  real(r8),parameter,dimension(nhydro) :: &
                 !  LSL     LSI         LSR         LSS       CVL     CVI        CVR          CVS          LSG
       N_ax    = (/-1._r8, -1._r8,     8.e6_r8,     3.e6_r8, -1._r8, -1._r8,     8.e6_r8,     3.e6_r8,     4.e6_r8/),&
       N_bx    = (/-1._r8, -1._r8,      0.0_r8,      0.0_r8, -1._r8, -1._r8,      0.0_r8,      0.0_r8,      0.0_r8/),&
       alpha_x = (/-1._r8, -1._r8,      0.0_r8,      0.0_r8, -1._r8, -1._r8,      0.0_r8,      0.0_r8,      0.0_r8/),&
       c_x     = (/-1._r8, -1._r8,    842.0_r8,     4.84_r8, -1._r8, -1._r8,    842.0_r8,     4.84_r8,     94.5_r8/),&
       d_x     = (/-1._r8, -1._r8,      0.8_r8,     0.25_r8, -1._r8, -1._r8,      0.8_r8,     0.25_r8,      0.5_r8/),&
       g_x     = (/-1._r8, -1._r8,      0.5_r8,      0.5_r8, -1._r8, -1._r8,      0.5_r8,      0.5_r8,      0.5_r8/),&
       a_x     = (/-1._r8, -1._r8,    524.0_r8,    52.36_r8, -1._r8, -1._r8,    524.0_r8,    52.36_r8,   209.44_r8/),&
       b_x     = (/-1._r8, -1._r8,      3.0_r8,      3.0_r8, -1._r8, -1._r8,      3.0_r8,      3.0_r8,      3.0_r8/),&
       gamma_1 = (/-1._r8, -1._r8, 17.83725_r8, 8.284701_r8, -1._r8, -1._r8, 17.83725_r8, 8.284701_r8, 11.63230_r8/),&
       gamma_2 = (/-1._r8, -1._r8,      6.0_r8,      6.0_r8, -1._r8, -1._r8,      6.0_r8,      6.0_r8,      6.0_r8/),&
       gamma_3 = (/-1._r8, -1._r8,      2.0_r8,      2.0_r8, -1._r8, -1._r8,      2.0_r8,      2.0_r8,      2.0_r8/),&
       gamma_4 = (/-1._r8, -1._r8,      6.0_r8,      6.0_r8, -1._r8, -1._r8,      6.0_r8,      6.0_r8,      6.0_r8/)       
#endif

  logical :: use_MMF

CONTAINS

  ! ######################################################################################
  ! SUBROUTINE setcosp2values
  ! ######################################################################################
#ifdef USE_COSP
  subroutine setcosp2values(Nlr_in,use_vgrid_in,csat_vgrid_in,Ncolumns_in,cosp_nradsteps_in)
    use mod_cosp,             only: cosp_init 
    use mod_cosp_config,      only: vgrid_zl, vgrid_zu, vgrid_z
    use mod_quickbeam_optics, only: hydro_class_init, quickbeam_optics_init
    ! Inputs
    integer, intent(in) :: Nlr_in             ! Number of vertical levels for CALIPSO and Cloudsat products
    integer, intent(in) :: Ncolumns_in        ! Number of sub-columns
    integer, intent(in) :: cosp_nradsteps_in  ! How often to call COSP?
    logical, intent(in) :: use_vgrid_in       ! Logical switch to use interpolated, to Nlr_in, grid for CALIPSO and Cloudsat
    logical, intent(in) :: csat_vgrid_in      !
    
    ! Local
    logical :: ldouble=.false.
    logical :: lsingle=.true. ! Default is to use single moment
    integer :: i,k

    prsmid_cosp  = pres_binCenters
    prslim_cosp  = pres_binEdges
    taumid_cosp  = tau_binCenters
    taulim_cosp  = tau_binEdges
    srmid_cosp   = calipso_binCenters
    srlim_cosp   = calipso_binEdges
    sza_cosp     = parasol_sza
    dbzemid_cosp = cloudsat_binCenters
    dbzelim_cosp = cloudsat_binEdges
    htmisrmid_cosp = misr_histHgtCenters
    htmisrlim_cosp = misr_histHgtEdges
    taumid_cosp_modis = tau_binCenters
    taulim_cosp_modis = tau_binEdges
    reffICE_binCenters_cosp = reffICE_binCenters
    reffICE_binEdges_cosp   = reffICE_binEdges
    reffLIQ_binCenters_cosp = reffLIQ_binCenters
    reffLIQ_binEdges_cosp   = reffLIQ_binEdges
                                  
    ! Initialize the distributional parameters for hydrometeors in radar simulator. In COSPv1.4, this was declared in
    ! cosp_defs.f.
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
       ldouble = .true. 
       lsingle = .false.
    end if
    call hydro_class_init(lsingle,ldouble,sd)
    call quickbeam_optics_init()

    ! DS2017: The setting up of the vertical grid for regridding the CALIPSO and Cloudsat products is 
    !         now donein cosp_init, but these fields are stored in cosp_config.F90.
    !         Additionally all static fields used by the individual simulators are set up by calls
    !         to _init functions in cosp_init.
    ! DS2019: Add logicals, default=.false., for new Lidar simuldators (Earthcare (atlid) and ground-based
    !         lidar at 532nm)
    call COSP_INIT(Lisccp_sim, Lmodis_sim, Lmisr_sim, Lradar_sim, Llidar_sim, LgrLidar532, &
         Latlid, Lparasol_sim, Lrttov_sim, radar_freq, k2, use_gas_abs, do_ray,              &
         isccp_topheight, isccp_topheight_direction, surface_radar, rcfg_cloudsat,           &
         use_vgrid_in, csat_vgrid_in, Nlr_in, pver, cloudsat_micro_scheme)

    ! Set number of sub-columns, from namelist
    nscol_cosp = Ncolumns_in
    
    if (use_vgrid_in) then  ! using fixed vertical grid
       if (csat_vgrid_in)       then
          nht_cosp = 40
       else
          nht_cosp = Nlr_in
       end if
    end if
    
    ! Set COSP call frequency, from namelist.
    cosp_nradsteps = cosp_nradsteps_in
    
    ! DJS2017: In COSP2, most of the bin boundaries, centers, and edges are declared in src/cosp_config.F90.
    !          Above I just assign them accordingly in the USE statement. Other bin bounds needed by CAM 
    !          are calculated here.
    ! Allocate
    allocate(htlim_cosp(2,nht_cosp),htmid_cosp(nht_cosp),scol_cosp(nscol_cosp))
    
    htmid_cosp      = vgrid_z
    htlim_cosp(1,:) = vgrid_zu
    htlim_cosp(2,:) = vgrid_zl

    scol_cosp(:) = (/(k,k=1,nscol_cosp)/)
    
    !  Just using an index here, model height is a prognostic variable
    htmlmid_cosp(:) = (/(k,k=1,nhtml_cosp)/)
    
  end subroutine setcosp2values
#endif    

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_readnl
  ! 
  ! PURPOSE: to read namelist variables and run setcospvalues subroutine.note: cldfrc_readnl
  ! is a good template in cloud_fraction.F90. Make sure that this routine is reading in a 
  ! namelist. models/atm/cam/bld/build-namelist is the perl script to check.
  ! ######################################################################################
  subroutine cospsimulator_intr_readnl(nlfile)
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
#ifdef SPMD
    use mpishorthand,    only: mpicom, mpilog, mpiint, mpichar
#endif

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input  (nlfile=atm_in)

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cospsimulator_intr_readnl'

#ifdef USE_COSP
    ! this list should include any variable that you might want to include in the namelist
    ! philosophy is to not include COSP output flags but just important COSP settings and cfmip controls. 
    namelist /cospsimulator_nl/ docosp, cosp_active, cosp_amwg, &
         cosp_histfile_num, cosp_histfile_aux, cosp_histfile_aux_num, cosp_isccp, cosp_lfrac_out, &
         cosp_lite, cosp_lradar_sim, cosp_llidar_sim, cosp_lisccp_sim,  cosp_lmisr_sim, cosp_lmodis_sim, cosp_ncolumns, &
         cosp_nradsteps, cosp_passive
    
    ! read in the namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'cospsimulator_nl', status=ierr)   
       if (ierr == 0) then
          read(unitn, cospsimulator_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
    
#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(docosp,               1,  mpilog, 0, mpicom)
    call mpibcast(cosp_amwg,            1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lite,            1,  mpilog, 0, mpicom)
    call mpibcast(cosp_passive,         1,  mpilog, 0, mpicom)
    call mpibcast(cosp_active,          1,  mpilog, 0, mpicom)
    call mpibcast(cosp_isccp,           1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lfrac_out,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lradar_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_llidar_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lisccp_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lmisr_sim,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lmodis_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_ncolumns,        1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_num,    1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_aux_num,1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_aux,    1,  mpilog, 0, mpicom)
    call mpibcast(cosp_nradsteps,       1,  mpiint, 0, mpicom)
#endif   
    
    if (cosp_lfrac_out) then
       lfrac_out = .true.
    end if
    if (cosp_lradar_sim) then
       lradar_sim = .true.
    end if
    if (cosp_llidar_sim) then
       llidar_sim = .true.
       lparasol_sim = .true.
    end if
    if (cosp_lisccp_sim) then
       lisccp_sim = .true.
    end if
    if (cosp_lmisr_sim) then
       lmisr_sim = .true.
    end if
    if (cosp_lmodis_sim) then
       lmodis_sim = .true.
    end if
    
    if (cosp_histfile_aux .and. cosp_histfile_aux_num == -1) then
       cosp_histfile_aux_num = cosp_histfile_num
    end if
    
    if (cosp_lite) then
       llidar_sim = .true.
       lparasol_sim = .true.
       lisccp_sim = .true.
       lmisr_sim = .true.
       lmodis_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    if (cosp_passive) then
       lisccp_sim = .true.
       lmisr_sim = .true.
       lmodis_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    if (cosp_active) then
       lradar_sim = .true.
       llidar_sim = .true.
       lparasol_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    if (cosp_isccp) then
       lisccp_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    ! if no simulators are turned on at all and docosp is, set cosp_amwg = .true.
    if((docosp) .and. (.not.lradar_sim) .and. (.not.llidar_sim) .and. (.not.lisccp_sim) .and. &
         (.not.lmisr_sim) .and. (.not.lmodis_sim)) then
       cosp_amwg = .true.
    end if
    if (cosp_amwg) then
       lradar_sim = .true.
       llidar_sim = .true.
       lparasol_sim = .true.
       lisccp_sim = .true.
       lmisr_sim = .true.
       lmodis_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    ! reset COSP namelist variables based on input from cam namelist variables
    if (cosp_ncolumns .ne. ncolumns) then
       ncolumns = cosp_ncolumns
    end if
        
    ! *NOTE* COSP is configured in CAM such that if a simulator is requested, all diagnostics
    ! are output. So no need turn on/aff outputs if simulator is requested.

    ! Set vertical coordinate, subcolumn, and calculation frequency cosp options based on namelist inputs
    call setcosp2values(nlr,use_vgrid,csat_vgrid,ncolumns,cosp_nradsteps)

    if (masterproc) then
       if (docosp) then 
          write(iulog,*)'COSP configuration:'
          write(iulog,*)'  Number of COSP subcolumns                = ', cosp_ncolumns
          write(iulog,*)'  Frequency at which cosp is called        = ', cosp_nradsteps
          write(iulog,*)'  Enable radar simulator                   = ', lradar_sim
          write(iulog,*)'  Enable calipso simulator                 = ', llidar_sim
          write(iulog,*)'  Enable ISCCP simulator                   = ', lisccp_sim
          write(iulog,*)'  Enable MISR simulator                    = ', lmisr_sim
          write(iulog,*)'  Enable MODIS simulator                   = ', lmodis_sim
          write(iulog,*)'  RADAR_SIM microphysics scheme            = ', trim(cloudsat_micro_scheme)
          write(iulog,*)'  Write COSP output to history file        = ', cosp_histfile_num
          write(iulog,*)'  Write COSP input fields                  = ', cosp_histfile_aux
          write(iulog,*)'  Write COSP input fields to history file  = ', cosp_histfile_aux_num
          write(iulog,*)'  Write COSP subcolumn fields              = ', cosp_lfrac_out
       else
          write(iulog,*)'COSP not enabled'
       end if
    end if
#endif
  end subroutine cospsimulator_intr_readnl

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_register
  ! ######################################################################################
  subroutine cospsimulator_intr_register()

    use cam_history_support, only: add_hist_coord
    
#ifdef USE_COSP
    ! register non-standard variable dimensions
    if (lisccp_sim .or. lmodis_sim) then
       call add_hist_coord('cosp_prs', nprs_cosp, 'COSP Mean ISCCP pressure',  &
            'Pa', prsmid_cosp, bounds_name='cosp_prs_bnds', bounds=prslim_cosp)
    end if
    
    if (lisccp_sim .or. lmisr_sim) then
       call add_hist_coord('cosp_tau', ntau_cosp,                              &
            'COSP Mean ISCCP optical depth', '1', taumid_cosp,                 &
            bounds_name='cosp_tau_bnds', bounds=taulim_cosp)
    end if
    
    if (lisccp_sim .or. llidar_sim .or. lradar_sim .or. lmisr_sim) then
       call add_hist_coord('cosp_scol', nscol_cosp, 'COSP subcolumn',          &
            values=scol_cosp)
    end if
    
    if (llidar_sim .or. lradar_sim) then
       call add_hist_coord('cosp_ht', nht_cosp,                                &
            'COSP Mean Height for calipso and radar simulator outputs', 'm',     &
            htmid_cosp, bounds_name='cosp_ht_bnds', bounds=htlim_cosp,         &
            vertical_coord=.true.)
    end if
    
    if (llidar_sim) then
       call add_hist_coord('cosp_sr', nsr_cosp,                                &
            'COSP Mean Scattering Ratio for calipso simulator CFAD output', '1', &
            srmid_cosp, bounds_name='cosp_sr_bnds', bounds=srlim_cosp)
    end if
    
    if (llidar_sim) then
       call add_hist_coord('cosp_sza', nsza_cosp, 'COSP Parasol SZA',          &
            'degrees', sza_cosp)
    end if
    
    if (lradar_sim) then
       call add_hist_coord('cosp_dbze', CLOUDSAT_DBZE_BINS,                            &
            'COSP Mean dBZe for radar simulator CFAD output', 'dBZ',           &
            dbzemid_cosp, bounds_name='cosp_dbze_bnds', bounds=dbzelim_cosp)
    end if
    
    if (lmisr_sim) then
       call add_hist_coord('cosp_htmisr', nhtmisr_cosp, 'COSP MISR height', &
            'm', htmisrmid_cosp,                                           &
            bounds_name='cosp_htmisr_bnds', bounds=htmisrlim_cosp)
    end if
    
    if (lmodis_sim) then
       call add_hist_coord('cosp_tau_modis', ntau_cosp_modis,                  &
            'COSP Mean MODIS optical depth', '1', taumid_cosp_modis,           &
            bounds_name='cosp_tau_modis_bnds', bounds=taulim_cosp_modis)
       call add_hist_coord('cosp_reffice',numMODISReffIceBins,                 &
            'COSP Mean MODIS effective radius (ice)', 'm', reffICE_binCenters_cosp, &
            bounds_name='cosp_reffice_bnds',bounds=reffICE_binEdges_cosp)
       call add_hist_coord('cosp_reffliq',numMODISReffLiqBins,                 &
            'COSP Mean MODIS effective radius (liquid)', 'm', reffLIQ_binCenters_cosp, &
            bounds_name='cosp_reffliq_bnds',bounds=reffLIQ_binEdges_cosp)      
    end if
    
#endif
  end subroutine cospsimulator_intr_register
  
  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_init
  ! ######################################################################################
  subroutine cospsimulator_intr_init()

#ifdef USE_COSP     

    use cam_history,         only: addfld, add_default, horiz_only
#ifdef SPMD
    use mpishorthand,        only : mpir8, mpiint, mpicom
#endif
    use netcdf,              only : nf90_open, nf90_inq_varid, nf90_get_var, nf90_close, nf90_nowrite
    use error_messages,      only : handle_ncerr, alloc_err
    
    use physics_buffer,  only: pbuf_get_index

    use mod_cosp_config,  only : R_UNDEF    
    
    integer :: errcode
    integer :: i

    ! Are we using MMF?
    call phys_getopts(use_MMF_out=use_MMF)
    
    ! ISCCP OUTPUTS
    if (lisccp_sim) then
       ! addfld calls for all
       call addfld('FISCCP1_COSP',(/'cosp_tau','cosp_prs'/),'A','percent', &
            'Grid-box fraction covered by each ISCCP D level cloud type',&
            flag_xyfill=.true., fill_value=R_UNDEF)
       
       !*cfMon,cfDa* tclisccp (time,profile), CFMIP wants "gridbox mean cloud cover from ISCCP"
       call addfld('CLDTOT_ISCCP', horiz_only,'A','percent', &
            'Total Cloud Fraction Calculated by the ISCCP Simulator ',flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfMon,cfDa* albisccp (time,profile)
       ! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
       call addfld('MEANCLDALB_ISCCP',horiz_only,'A','1','Mean cloud albedo*CLDTOT_ISCCP',flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfMon,cfDa* ctpisccp (time,profile)
       ! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
       call addfld('MEANPTOP_ISCCP',horiz_only,'A','Pa','Mean cloud top pressure*CLDTOT_ISCCP',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       ! tauisccp (time,profile)
       ! For averaging, weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
       call addfld ('MEANTAU_ISCCP',horiz_only,'A','1','Mean optical thickness*CLDTOT_ISCCP',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       ! meantbisccp (time,profile), at 10.5 um
       call addfld ('MEANTB_ISCCP',horiz_only,'A','K','Mean Infrared Tb from ISCCP simulator',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       ! meantbclrisccp (time,profile)
       call addfld ('MEANTBCLR_ISCCP',horiz_only,'A','K','Mean Clear-sky Infrared Tb from ISCCP simulator',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! boxtauisccp (time,column,profile)
       call addfld ('TAU_ISCCP',(/'cosp_scol'/),'I','1','Optical Depth in each Subcolumn',flag_xyfill=.true., fill_value=R_UNDEF)
       ! boxptopisccp (time,column,profile)
       call addfld ('CLDPTOP_ISCCP',(/'cosp_scol'/),'I','Pa','Cloud Top Pressure in each Subcolumn',  &
            flag_xyfill=.true., fill_value=R_UNDEF)

       ! add all isccp outputs to the history file specified by the CAM namelist variable cosp_histfile_num
       call add_default ('FISCCP1_COSP',cosp_histfile_num,' ')
       call add_default ('CLDTOT_ISCCP',cosp_histfile_num,' ')
       call add_default ('MEANCLDALB_ISCCP',cosp_histfile_num,' ')
       call add_default ('MEANPTOP_ISCCP',cosp_histfile_num,' ')
       call add_default ('MEANTAU_ISCCP',cosp_histfile_num,' ')
       call add_default ('MEANTB_ISCCP',cosp_histfile_num,' ')
       call add_default ('MEANTBCLR_ISCCP',cosp_histfile_num,' ')
      
    end if

    ! CALIPSO SIMULATOR OUTPUTS
    if (llidar_sim) then
       ! addfld calls for all
       !*cfMon,cfOff,cfDa,cf3hr* cllcalipso (time,profile)
       call addfld('CLDLOW_CAL',horiz_only,'A','percent','Calipso Low-level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfMon,cfOff,cfDa,cf3hr* clmcalipso (time,profile)
       call addfld('CLDMED_CAL',horiz_only,'A','percent','Calipso Mid-level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfMon,cfOff,cfDa,cf3hr* clhcalipso (time,profile)
       call addfld('CLDHGH_CAL',horiz_only,'A','percent','Calipso High-level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfMon,cfOff,cfDa,cf3hr* cltcalipso (time,profile)
       call addfld('CLDTOT_CAL',horiz_only,'A','percent','Calipso Total Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfMon,cfOff,cfDa,cf3hr* clcalipso (time,height,profile)
       call addfld('CLD_CAL',(/'cosp_ht'/),'A','percent','Calipso Cloud Fraction (532 nm)', flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfMon,cfOff,cfDa,cf3hr* parasol_refl (time,sza,profile)
       call addfld ('RFL_PARASOL',(/'cosp_sza'/),'A','1','PARASOL-like mono-directional reflectance ',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       !*cfOff,cf3hr* cfad_calipsosr532 (time,height,scat_ratio,profile), %11%, default is 40 vert levs, 15 SR  bins
       call addfld('CFAD_SR532_CAL',(/'cosp_sr','cosp_ht'/),'A','1',                                    &
            'Calipso Scattering Ratio CFAD (532 nm)',                                                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MOL532_CAL',(/'lev'/),'A','m-1sr-1','Calipso Molecular Backscatter (532 nm) ',              &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ATB532_CAL',(/'cosp_scol','lev      '/),'I','no_unit_log10(x)',                           &
            'Calipso Attenuated Total Backscatter (532 nm) in each Subcolumn',                        &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclcalipsoliq (time,alt40,loc)
       call addfld('CLD_CAL_LIQ', (/'cosp_ht'/), 'A','percent', 'Calipso Liquid Cloud Fraction',                 &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclcalipsoice (time,alt40,loc)
       call addfld('CLD_CAL_ICE', (/'cosp_ht'/), 'A','percent', 'Calipso Ice Cloud Fraction',                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclcalipsoun (time,alt40,loc)
       call addfld('CLD_CAL_UN', (/'cosp_ht'/),'A','percent', 'Calipso Undefined-Phase Cloud Fraction',          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclcalipsotmp (time,alt40,loc)
       call addfld('CLD_CAL_TMP', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',        &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclcalipsotmpliq (time,alt40,loc)
       call addfld('CLD_CAL_TMPLIQ', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclcalipsotmpice (time,alt40,loc)
       call addfld('CLD_CAL_TMPICE', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclcalipsotmpun (time,alt40,loc)
       call addfld('CLD_CAL_TMPUN', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',      &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lcltcalipsoice (time,loc)
       call addfld('CLDTOT_CAL_ICE', horiz_only,'A','percent','Calipso Total Ice Cloud Fraction',    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lcltcalipsoliq (time,loc)
       call addfld('CLDTOT_CAL_LIQ', horiz_only,'A','percent','Calipso Total Liquid Cloud Fraction', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lcltcalipsoun (time,loc)
       call addfld('CLDTOT_CAL_UN',horiz_only,'A','percent','Calipso Total Undefined-Phase Cloud Fraction',      &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclhcalipsoice (time,loc)
       call addfld('CLDHGH_CAL_ICE',horiz_only,'A','percent','Calipso High-level Ice Cloud Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclhcalipsoliq (time,loc)
       call addfld('CLDHGH_CAL_LIQ',horiz_only,'A','percent','Calipso High-level Liquid Cloud Fraction',         &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclhcalipsoun (time,loc)
       call addfld('CLDHGH_CAL_UN',horiz_only,'A','percent','Calipso High-level Undefined-Phase Cloud Fraction', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclmcalipsoice (time,loc)
       call addfld('CLDMED_CAL_ICE',horiz_only,'A','percent','Calipso Mid-level Ice Cloud Fraction',             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclmcalipsoliq (time,loc)
       call addfld('CLDMED_CAL_LIQ',horiz_only,'A','percent','Calipso Mid-level Liquid Cloud Fraction',          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lclmcalipsoun (time,loc)
       call addfld('CLDMED_CAL_UN',horiz_only,'A','percent','Calipso Mid-level Undefined-Phase Cloud Fraction',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lcllcalipsoice (time,loc)
       call addfld('CLDLOW_CAL_ICE',horiz_only,'A','percent','Calipso Low-level Ice Cloud Fraction',             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lcllcalipsoliq (time,loc)
       call addfld('CLDLOW_CAL_LIQ',horiz_only,'A','percent','Calipso Low-level Liquid Cloud Fraction',          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! lcllcalipsoun (time,loc) !+cosp1.4
       call addfld('CLDLOW_CAL_UN',horiz_only,'A','percent','Calipso Low-level Undefined-Phase Cloud Fraction',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
    
       ! add_default calls for CFMIP experiments or else all fields are added to history file
       !     except those with sub-column dimension/experimental variables
       ! add all calipso outputs to the history file specified by the CAM namelist variable cosp_histfile_num
       call add_default ('CLDLOW_CAL',cosp_histfile_num,' ')
       call add_default ('CLDMED_CAL',cosp_histfile_num,' ')
       call add_default ('CLDHGH_CAL',cosp_histfile_num,' ')
       call add_default ('CLDTOT_CAL',cosp_histfile_num,' ')
       call add_default ('CLD_CAL',cosp_histfile_num,' ')
       call add_default ('RFL_PARASOL',cosp_histfile_num,' ')
       call add_default ('CFAD_SR532_CAL',cosp_histfile_num,' ')
       call add_default ('CLD_CAL_LIQ',cosp_histfile_num,' ')  !+COSP1.4
       call add_default ('CLD_CAL_ICE',cosp_histfile_num,' ')
       call add_default ('CLD_CAL_UN',cosp_histfile_num,' ')
       call add_default ('CLDTOT_CAL_ICE',cosp_histfile_num,' ')
       call add_default ('CLDTOT_CAL_LIQ',cosp_histfile_num,' ')
       call add_default ('CLDTOT_CAL_UN',cosp_histfile_num,' ')
       call add_default ('CLDHGH_CAL_ICE',cosp_histfile_num,' ')
       call add_default ('CLDHGH_CAL_LIQ',cosp_histfile_num,' ')
       call add_default ('CLDHGH_CAL_UN',cosp_histfile_num,' ')
       call add_default ('CLDMED_CAL_ICE',cosp_histfile_num,' ')
       call add_default ('CLDMED_CAL_LIQ',cosp_histfile_num,' ')
       call add_default ('CLDMED_CAL_UN',cosp_histfile_num,' ')
       call add_default ('CLDLOW_CAL_ICE',cosp_histfile_num,' ')
       call add_default ('CLDLOW_CAL_LIQ',cosp_histfile_num,' ')
       call add_default ('CLDLOW_CAL_UN',cosp_histfile_num,' ')
       if ((.not.cosp_amwg) .and. (.not.cosp_lite) .and. (.not.cosp_passive) .and. (.not.cosp_active) &
            .and. (.not.cosp_isccp)) then
          call add_default ('MOL532_CAL',cosp_histfile_num,' ')
       end if
    end if

    ! RADAR SIMULATOR OUTPUTS
    if (lradar_sim) then

       allocate(sd_cs(begchunk:endchunk), rcfg_cs(begchunk:endchunk))
       do i = begchunk, endchunk
          sd_cs(i)   = sd
          rcfg_cs(i) = rcfg_cloudsat
       end do

       ! addfld calls
       call addfld('CFAD_DBZE94_CS',(/'cosp_dbze','cosp_ht  '/),'A','1',&
            'Radar Reflectivity Factor CFAD (94 GHz)',&
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLD_CAL_NOTCS',(/'cosp_ht'/),'A','percent','Cloud occurrence seen by CALIPSO but not CloudSat ',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLDTOT_CALCS',horiz_only,'A','percent',' Calipso and Radar Total Cloud Fraction ',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       call addfld ('CLDTOT_CS',horiz_only,'A','percent',' Radar total cloud amount ',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLDTOT_CS2',horiz_only,'A','percent', &
            ' Radar total cloud amount without the data for the first kilometer above surface ', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('DBZE_CS',(/'cosp_scol','lev      '/),'I','dBZe',' Radar dBZe (94 GHz) in each Subcolumn',&
            flag_xyfill=.true., fill_value=R_UNDEF)

       ! Cloudsat near-sfc precipitation diagnostics
       call addfld('CS_NOPRECIP',  horiz_only, 'A', '1',    'CloudSat No Rain Fraction',                   flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINPOSS',  horiz_only, 'A', '1',    'Cloudsat Rain Possible Fraction',             flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINPROB',  horiz_only, 'A', '1',    'CloudSat Rain Probable Fraction',             flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINCERT',  horiz_only, 'A', '1',    'CloudSat Rain Certain Fraction',              flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_SNOWPOSS',  horiz_only, 'A', '1',    'CloudSat Snow Possible Fraction',             flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_SNOWCERT',  horiz_only, 'A', '1',    'CloudSat Snow Certain Fraction',              flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_MIXPOSS',   horiz_only, 'A', '1',    'CloudSat Mixed Possible Fraction',            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_MIXCERT',   horiz_only, 'A', '1',    'CloudSat Mixed Certain Fraction',             flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINHARD',  horiz_only, 'A', '1',    'CloudSat Heavy Rain Fraction',                flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_UN',        horiz_only, 'A', '1',    'CloudSat Unclassified Precipitation Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_PIA',       horiz_only, 'A', 'dBZ',  'CloudSat Radar Path Integrated Attenuation',  flag_xyfill=.true., fill_value=R_UNDEF)
       ! Associated CAM microphysics
       !call addfld('CAM_MP_CVRAIN',horiz_only, 'A', 'kg/kg','CAM Microphysics Convective Rain',            flag_xyfill=.true., fill_value=R_UNDEF)
       !call addfld('CAM_MP_CVSNOW',horiz_only, 'A', 'kg/kg','CAM Microphysics Convective Snow',            flag_xyfill=.true., fill_value=R_UNDEF)
       !call addfld('CAM_MP_LSRAIN',horiz_only, 'A', 'kg/kg','CAM Microphysics Large-Scale Rain',           flag_xyfill=.true., fill_value=R_UNDEF)
       !call addfld('CAM_MP_LSSNOW',horiz_only, 'A', 'kg/kg','CAM Microphysics Large-Scale Snow',           flag_xyfill=.true., fill_value=R_UNDEF)
       !call addfld('CAM_MP_LSGRPL',horiz_only, 'A', 'kg/kg','CAM Microphysics Large-Scale Graupel',        flag_xyfill=.true., fill_value=R_UNDEF)

       ! add_default calls for CFMIP experiments or else all fields are added to history file except those with sub-column dimension
       ! add all radar outputs to the history file specified by the CAM namelist variable cosp_histfile_num
       call add_default ('CFAD_DBZE94_CS',cosp_histfile_num,' ')
       call add_default ('CLD_CAL_NOTCS', cosp_histfile_num,' ')
       call add_default ('CLDTOT_CALCS',  cosp_histfile_num,' ')
       call add_default ('CLDTOT_CS',     cosp_histfile_num,' ')
       call add_default ('CLDTOT_CS2',    cosp_histfile_num,' ')
       call add_default ('CS_NOPRECIP',   cosp_histfile_num,' ')
       call add_default ('CS_RAINPOSS',   cosp_histfile_num,' ')
       call add_default ('CS_RAINPROB',   cosp_histfile_num,' ')
       call add_default ('CS_RAINCERT',   cosp_histfile_num,' ')
       call add_default ('CS_SNOWPOSS',   cosp_histfile_num,' ')
       call add_default ('CS_SNOWCERT',   cosp_histfile_num,' ')
       call add_default ('CS_MIXPOSS',    cosp_histfile_num,' ')
       call add_default ('CS_MIXCERT',    cosp_histfile_num,' ')
       call add_default ('CS_RAINHARD',   cosp_histfile_num,' ')
       call add_default ('CS_UN',         cosp_histfile_num,' ')
       call add_default ('CS_PIA',        cosp_histfile_num,' ')
    end if
    
    ! MISR SIMULATOR OUTPUTS
    if (lmisr_sim) then
       call addfld ('CLD_MISR',(/'cosp_tau   ','cosp_htmisr'/),'A','percent','Cloud Fraction from MISR Simulator',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call add_default ('CLD_MISR',cosp_histfile_num,' ')
    end if

    ! MODIS OUTPUT
    if (lmodis_sim) then
       call addfld ('CLTMODIS',horiz_only,'A','%','MODIS Total Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLWMODIS',horiz_only,'A','%','MODIS Liquid Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLIMODIS',horiz_only,'A','%','MODIS Ice Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLHMODIS',horiz_only,'A','%','MODIS High Level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLMMODIS',horiz_only,'A','%','MODIS Mid Level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLLMODIS',horiz_only,'A','%','MODIS Low Level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUTMODIS',horiz_only,'A','1','MODIS Total Cloud Optical Thickness*CLTMODIS',                  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUWMODIS',horiz_only,'A','1','MODIS Liquid Cloud Optical Thickness*CLWMODIS',                 &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUIMODIS',horiz_only,'A','1','MODIS Ice Cloud Optical Thickness*CLIMODIS',                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUTLOGMODIS',horiz_only,'A','1','MODIS Total Cloud Optical Thickness (Log10 Mean)*CLTMODIS',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUWLOGMODIS',horiz_only,'A','1','MODIS Liquid Cloud Optical Thickness (Log10 Mean)*CLWMODIS', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUILOGMODIS',horiz_only,'A','1','MODIS Ice Cloud Optical Thickness (Log10 Mean)*CLIMODIS',    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('REFFCLWMODIS',horiz_only,'A','m','MODIS Liquid Cloud Particle Size*CLWMODIS',                  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('REFFCLIMODIS',horiz_only,'A','m','MODIS Ice Cloud Particle Size*CLIMODIS',                     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('PCTMODIS',horiz_only,'A','Pa','MODIS Cloud Top Pressure*CLTMODIS',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('LWPMODIS',horiz_only,'A','kg m-2','MODIS Cloud Liquid Water Path*CLWMODIS',                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('IWPMODIS',horiz_only,'A','kg m-2','MODIS Cloud Ice Water Path*CLIMODIS',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLMODIS',(/'cosp_tau_modis','cosp_prs      '/),'A','%','MODIS Cloud Area Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLRIMODIS',(/'cosp_tau_modis','cosp_reffice  '/),'A','%','MODIS Cloud Area Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLRLMODIS',(/'cosp_tau_modis','cosp_reffliq  '/),'A','%','MODIS Cloud Area Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       
       ! add MODIS output to history file specified by the CAM namelist variable cosp_histfile_num
       call add_default ('CLTMODIS',cosp_histfile_num,' ')
       call add_default ('CLWMODIS',cosp_histfile_num,' ')
       call add_default ('CLIMODIS',cosp_histfile_num,' ')
       call add_default ('CLHMODIS',cosp_histfile_num,' ')
       call add_default ('CLMMODIS',cosp_histfile_num,' ')
       call add_default ('CLLMODIS',cosp_histfile_num,' ')
       call add_default ('TAUTMODIS',cosp_histfile_num,' ')
       call add_default ('TAUWMODIS',cosp_histfile_num,' ')
       call add_default ('TAUIMODIS',cosp_histfile_num,' ')
       call add_default ('TAUTLOGMODIS',cosp_histfile_num,' ')
       call add_default ('TAUWLOGMODIS',cosp_histfile_num,' ')
       call add_default ('TAUILOGMODIS',cosp_histfile_num,' ')
       call add_default ('REFFCLWMODIS',cosp_histfile_num,' ')
       call add_default ('REFFCLIMODIS',cosp_histfile_num,' ')
       call add_default ('PCTMODIS',cosp_histfile_num,' ')
       call add_default ('LWPMODIS',cosp_histfile_num,' ')
       call add_default ('IWPMODIS',cosp_histfile_num,' ')
       call add_default ('CLMODIS',cosp_histfile_num,' ')
       call add_default ('CLRIMODIS',cosp_histfile_num,' ')
       call add_default ('CLRLMODIS',cosp_histfile_num,' ')
    end if
    
    ! SUB-COLUMN OUTPUT
    if (lfrac_out) then
       ! frac_out (time,height_mlev,column,profile)
       call addfld ('SCOPS_OUT',(/'cosp_scol','lev      '/),'I','0=nocld,1=strcld,2=cnvcld','SCOPS Subcolumn output', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call add_default ('SCOPS_OUT',cosp_histfile_num,' ')
       ! save sub-column outputs from ISCCP if ISCCP is run
       if (lisccp_sim) then
          call add_default ('TAU_ISCCP',cosp_histfile_num,' ')
          call add_default ('CLDPTOP_ISCCP',cosp_histfile_num,' ')
       end if
       ! save sub-column outputs from calipso if calipso is run
       if (llidar_sim) then
          call add_default ('ATB532_CAL',cosp_histfile_num,' ')
       end if
       ! save sub-column outputs from radar if radar is run
       if (lradar_sim) then
          call add_default ('DBZE_CS',cosp_histfile_num,' ')
       end if
    end if
    
    ! ADDFLD, ADD_DEFAULT, OUTFLD CALLS FOR COSP OUTPUTS IF RUNNING COSP OFF-LINE
    ! Note: A suggestion was to add all of the CAM variables needed to add to make it possible to run COSP off-line
    ! These fields are available and can be called from the namelist though.
    if (cosp_histfile_aux) then
       call addfld ('PS_COSP',         horiz_only,            'I','Pa',     'PS_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TS_COSP',         horiz_only,            'I','K',      'TS_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('P_COSP',          (/            'lev'/), 'I','Pa',     'P_COSP',                             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('PH_COSP',         (/            'lev'/), 'I','Pa',     'PH_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ZLEV_COSP',       (/            'lev'/), 'I','m',      'ZLEV_COSP',                          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ZLEV_HALF_COSP',  (/            'lev'/), 'I','m',      'ZLEV_HALF_COSP',                     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('T_COSP',          (/            'lev'/), 'I','K',      'T_COSP',                             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('RH_COSP',         (/            'lev'/), 'I','percent','RH_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('Q_COSP',          (/            'lev'/), 'I','kg/kg',  'Q_COSP',                             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAU_067',         (/'cosp_scol','lev      '/), 'I','1',      'Subcolumn 0.67micron optical depth', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('EMISS_11',        (/'cosp_scol','lev      '/), 'I','1',      'Subcolumn 11micron emissivity',      &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_fracliq',   (/'cosp_scol','lev      '/), 'I','1',      'Fraction of tau from liquid water',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_asym',      (/'cosp_scol','lev      '/), 'I','1',      'Assymetry parameter (MODIS)',        &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_ssa',       (/'cosp_scol','lev      '/), 'I','1',      'Single-scattering albedo (MODIS)',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_betatot',     (/'cosp_scol','lev      '/), 'I','1',      'Backscatter coefficient (CALIPSO)',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_betatot_ice', (/'cosp_scol','lev      '/), 'I','1',      'Backscatter coefficient (CALIPSO)',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_betatot_liq', (/'cosp_scol','lev      '/), 'I','1',      'Backscatter coefficient (CALIPSO)',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_tautot',      (/'cosp_scol','lev      '/), 'I','1',      'Vertically integrated ptical-depth (CALIPSO)', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_tautot_ice',  (/'cosp_scol','lev      '/), 'I','1',      'Vertically integrated ptical-depth (CALIPSO)', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_tautot_liq',  (/'cosp_scol','lev      '/), 'I','1',      'Vertically integrated ptical-depth (CALIPSO)', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_z_vol',        (/'cosp_scol','lev      '/), 'I','1',      'Effective reflectivity factor (CLOUDSAT)',     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_kr_vol',       (/'cosp_scol','lev      '/), 'I','1',      'Attenuation coefficient (hydro) (CLOUDSAT)',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_g_vol',        (/'cosp_scol','lev      '/), 'I','1',      'Attenuation coefficient (gases) (CLOUDSAT)',   &
            flag_xyfill=.true., fill_value=R_UNDEF)

       call add_default ('PS_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('TS_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('P_COSP',          cosp_histfile_aux_num,' ')
       call add_default ('PH_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('ZLEV_COSP',       cosp_histfile_aux_num,' ')
       call add_default ('ZLEV_HALF_COSP',  cosp_histfile_aux_num,' ')
       call add_default ('T_COSP',          cosp_histfile_aux_num,' ')
       call add_default ('RH_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('TAU_067',         cosp_histfile_aux_num,' ')
       call add_default ('EMISS_11',        cosp_histfile_aux_num,' ')
       call add_default ('MODIS_fracliq',   cosp_histfile_aux_num,' ')
       call add_default ('MODIS_asym',      cosp_histfile_aux_num,' ')
       call add_default ('MODIS_ssa',       cosp_histfile_aux_num,' ')
       call add_default ('CAL_betatot',     cosp_histfile_aux_num,' ')
       call add_default ('CAL_betatot_ice', cosp_histfile_aux_num,' ')
       call add_default ('CAL_betatot_liq', cosp_histfile_aux_num,' ')
       call add_default ('CAL_tautot',      cosp_histfile_aux_num,' ')
       call add_default ('CAL_tautot_ice',  cosp_histfile_aux_num,' ')
       call add_default ('CAL_tautot_liq',  cosp_histfile_aux_num,' ')
       call add_default ('CS_z_vol',        cosp_histfile_aux_num,' ')
       call add_default ('CS_kr_vol',       cosp_histfile_aux_num,' ')
       call add_default ('CS_g_vol',        cosp_histfile_aux_num,' ')
    end if
    
    ! These fields may or may not exist in pbuf, depending on physics
    ! configuration. The extra particle sizes and precipitation fluxes will not
    ! exist for MMF, for example. Rather than add logic here that we would need
    ! to keep consistent with the accesses in cospsimulator_intr_run, pass an
    ! error code here to allow for the possibility that these do not exist. If
    ! a field does not exist in pbuf, the resulting index will be zero, and pbuf
    ! will abort the run if we try to access these during runtime.
    rei_idx        = pbuf_get_index('REI', errcode=errcode)
    rel_idx        = pbuf_get_index('REL', errcode=errcode)
    cld_idx        = pbuf_get_index('CLD', errcode=errcode)
    concld_idx     = pbuf_get_index('CONCLD', errcode=errcode)
    lsreffrain_idx = pbuf_get_index('LS_REFFRAIN', errcode=errcode)
    lsreffsnow_idx = pbuf_get_index('LS_REFFSNOW', errcode=errcode)
    cvreffliq_idx  = pbuf_get_index('CV_REFFLIQ', errcode=errcode)
    cvreffice_idx  = pbuf_get_index('CV_REFFICE', errcode=errcode)
    dpcldliq_idx   = pbuf_get_index('DP_CLDLIQ', errcode=errcode)
    dpcldice_idx   = pbuf_get_index('DP_CLDICE', errcode=errcode)
    shcldliq_idx   = pbuf_get_index('SH_CLDLIQ', errcode=errcode)
    shcldice_idx   = pbuf_get_index('SH_CLDICE', errcode=errcode)
    shcldliq1_idx  = pbuf_get_index('SH_CLDLIQ1', errcode=errcode)
    shcldice1_idx  = pbuf_get_index('SH_CLDICE1', errcode=errcode)
    dpflxprc_idx   = pbuf_get_index('DP_FLXPRC', errcode=errcode)
    dpflxsnw_idx   = pbuf_get_index('DP_FLXSNW', errcode=errcode)
    shflxprc_idx   = pbuf_get_index('SH_FLXPRC', errcode=errcode)
    shflxsnw_idx   = pbuf_get_index('SH_FLXSNW', errcode=errcode)
    lsflxprc_idx   = pbuf_get_index('LS_FLXPRC', errcode=errcode)
    lsflxsnw_idx   = pbuf_get_index('LS_FLXSNW', errcode=errcode)
    
#endif    
  end subroutine cospsimulator_intr_init

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_run
  ! ######################################################################################
  subroutine cospsimulator_intr_run(state,pbuf, cam_in,emis,coszrs,cld_swtau,snow_tau,snow_emis)    
    use physics_types,        only: physics_state
    use physics_buffer,       only: physics_buffer_desc
    use camsrfexch,           only: cam_in_t
    use interpolate_data,     only: lininterp_init,lininterp,lininterp_finish,interp_type
#ifdef USE_COSP
    use mod_cosp_config,      only: R_UNDEF,parasol_nrefl, Nlvgrid, vgrid_zl, vgrid_zu
    use mod_cosp,             only: cosp_simulator
#endif

    ! Inputs
    type(physics_state), intent(in),target  :: state
    type(physics_buffer_desc),      pointer :: pbuf(:)
    type(cam_in_t),      intent(in)         :: cam_in
    real(r8), intent(in) :: emis(pcols,pver)               ! cloud longwave emissivity
    real(r8), intent(in) :: coszrs(pcols)                  ! cosine solar zenith angle (to tell if day or night)
    real(r8), intent(in) :: cld_swtau(pcols,pver)          ! cld_swtau, read in using this variable
    real(r8), intent(in),optional :: snow_tau(pcols,pver)  ! grid-box mean SW snow optical depth
    real(r8), intent(in),optional :: snow_emis(pcols,pver) ! grid-box mean LW snow optical depth

#ifdef USE_COSP
    ! Local variables
    integer :: npoints ! Number of points to use in COSP call
    integer :: ncol    ! number of active atmospheric columns
    
    ! COSP-related local vars
    type(cosp_outputs)        :: cospOUT, cospOUTave  ! COSP simulator outputs
    type(cosp_optical_inputs) :: cospIN, cospINave    ! COSP optical (or derived?) fields needed by simulators
    type(cosp_column_inputs)  :: cospstateIN          ! COSP model fields needed by simulators
    
    ! COSP error handling
    character(len=256),dimension(100) :: cosp_status
    integer :: nerror, i


    ! Number of columns in this physics chunk
    ncol = state%ncol

    ! Number of points to use in COSP call; could be packed GCM/CRM if using MMF
    if (use_MMF) then
       npoints = ncol * crm_nx_rad * crm_ny_rad
    else
       npoints = ncol
    end if

    ! Construct COSP output derived type.
    call t_startf('cosp_construct_cosp_outputs')
    call construct_cosp_outputs(npoints,nscol_cosp,pver,Nlvgrid,0,cospOUT)
    call t_stopf('cosp_construct_cosp_outputs')

    ! Model state inputs
    call t_startf('cosp_construct_cospstateIN')
    call construct_cospstateIN(npoints,pver,0,cospstateIN)      
    call t_stopf('cosp_construct_cospstateIN')
    
    call t_startf('cosp_populate_cosp_gridbox')
    call populate_cosp_gridbox(state, pbuf, cam_in, coszrs, cospstateIN)
    call t_stopf('cosp_populate_cosp_gridbox')

    ! Optical inputs
    call t_startf('cosp_construct_cospIN')
    call construct_cospIN(npoints,nscol_cosp,pver,cospIN)
    call t_stopf('cosp_construct_cospIN')

    call t_startf('cosp_subsample_and_optics')
    call populate_cosp_subcol(npoints, state, pbuf, emis, cld_swtau, cospstateIN, cospIN, snow_tau, snow_emis)    
    call t_stopf('cosp_subsample_and_optics')

    ! Call COSP and check status
    call t_startf('cosp_simulator')
    cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx=1, stop_idx=npoints,debug=.false.)
    nerror = 0
    do i = 1, ubound(cosp_status, 1)
       if (len_trim(cosp_status(i)) > 0) then
          write(iulog,*) "cosp_simulator: ERROR: "//trim(cosp_status(i))
          nerror = nerror + 1
       end if
    end do
    if (nerror > 0) then
       call endrun('cospsimulator_intr_run: error return from cosp_simulator')
    end if
    call t_stopf('cosp_simulator')
  
    ! Write COSP inputs to output file for offline use.
    if (cosp_histfile_aux) then
       call t_startf("cosp_histfile_aux")
       call cosp_histfile_aux_out(state, cospstateIN, cospIN)
       call t_stopf("cosp_histfile_aux")
    end if

    ! Set dark-scenes to fill value. Only done for passive simulators
    ! TODO: revisit this! We should NOT have to do this here! The simulators
    ! should be masking night values for us.
    call t_startf('cosp_remask_passive')
    call cosp_remask_passive(cospstateIN, cospOUT)
    call t_stopf('cosp_remask_passive')

    if (use_MMF) then
       call t_startf('cosp_write_outputs')
       ! Allocate space for unpacked outputs
       call construct_cospIN(ncol,nscol_cosp,pver,cospINave)
       call construct_cosp_outputs(ncol,nscol_cosp,pver,Nlvgrid,0,cospOUTave)
       ! Average packed COSP outputs
       call cosp_unpack_outputs(ncol, crm_nx_rad, crm_ny_rad, cospIN, cospOUT, cospINave, cospOUTave)
       ! Write COSP outputs to history files
       call cosp_write_outputs(state, cospINave, cospOUTave)
       ! Free space for unpacked outputs
       call destroy_cospIN(cospINave)
       call destroy_cosp_outputs(cospOUTave)
       call t_stopf('cosp_write_outputs')
    else
       ! Write COSP outputs to history files
       call t_startf('cosp_write_outputs')
       call cosp_write_outputs(state, cospIN, cospOUT)
       call t_stopf('cosp_write_outputs')
    end if

    ! Clean up
    call t_startf('cosp_finalize')
    call destroy_cospIN(cospIN)
    call destroy_cospstateIN(cospstateIN)
    call destroy_cosp_outputs(cospOUT) 
    call t_stopf('cosp_finalize')
#endif /* USE_COSP */
  end subroutine cospsimulator_intr_run


#ifdef USE_COSP
  subroutine populate_cosp_gridbox(state, pbuf, cam_in, coszrs, cospstateIN)
    use physics_types,        only: physics_state
    use physics_buffer,       only: physics_buffer_desc
    use camsrfexch,           only: cam_in_t
    use rad_constituents,     only: rad_cnst_get_gas
    use physconst,            only: pi, gravit
    type(physics_state)       , intent(in),target  :: state
    type(physics_buffer_desc) , pointer            :: pbuf(:)
    type(cam_in_t)            , intent(in)         :: cam_in
    real(r8)                  , intent(in)         :: coszrs(:)  ! cosine solar zenith angle (to tell if day or night)
    type(cosp_column_inputs), intent(inout)        :: cospstateIN  ! COSP model fields needed by simulators
    ! CAM pointers to get variables from radiation interface (get from rad_cnst_get_gas)
    real(r8), pointer, dimension(:,:) :: q               ! specific humidity (kg/kg)
    real(r8), pointer, dimension(:,:) :: o3              ! Mass mixing ratio 03
    ! Other local variables
    integer :: npoints, ncol, i, j, k, ix, iy

    ! ncol is number of points in model, npoints is number of points used in
    ! COSP. For MMF, npoints will be equal to ncol * crm_nx_rad * crm_ny_rad
    ncol = state%ncol
    npoints = cospstateIN%npoints

    ! radiative constituent interface variables:
    ! specific humidity (q), 03, CH4,C02, N20 mass mixing ratio
    ! Note: these all have dims (pcol,pver) but the values don't change much for the well-mixed gases.
    call rad_cnst_get_gas(0,'H2O', state, pbuf,  q)                     
    call rad_cnst_get_gas(0,'O3',  state, pbuf,  o3)

    ! Populate 2D fields
    do iy = 1,crm_ny_rad
      do ix = 1,crm_nx_rad
        do i = 1,ncol
          j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
          ! Coordinate variables
          cospstateIN%lat(j)                      = state%lat(i)*180._r8/ pi
          cospstateIN%lon(j)                      = state%lon(i)*180._r8 / pi
          ! Compute sunlit flag
          if (coszrs(i) > 0.0_r8) then
            cospstateIN%sunlit(j) = 1
          else
            cospstateIN%sunlit(j) = 0
          end if
          ! Compute land mask
          if (cam_in%landfrac(i) > 0.01_r8) then
            cospstateIN%land(j) = 1
          else
            cospstateIN%land(j) = 0
          end if
          cospstateIN%surfelev(j) = state%zi(i,pver+1) + state%phis(i) / gravit
          cospstateIN%skt(j) = cam_in%ts(i)
        end do
      end do
    end do

    ! Populate 3D fields
    do k = 1,pver
      do iy = 1,crm_ny_rad
        do ix = 1,crm_nx_rad
          do i = 1,ncol
            j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
            cospstateIN%at(j,k)                        = state%t(i,k) 
            cospstateIN%qv(j,k)                        = q(i,k)
            cospstateIN%o3(j,k)                        = o3(i,k)  
            cospstateIN%pfull(j,k)                     = state%pmid(i,k)
            cospstateIN%phalf(j,k) = state%pint(i,k) !0._r8
            cospstateIN%phalf(j,k+1) = state%pint(i,k+1) !0._r8
            ! add surface height (surface geopotential/gravity) to convert CAM heights
            ! based on geopotential above surface into height above sea level
            cospstateIN%hgt_matrix(j,k)          = state%zm(i,k) + state%phis(i) / gravit
            cospstateIN%hgt_matrix_half(j,k)     = state%zi(i,k+1) + state%phis(i) / gravit
          end do
        end do
      end do
    end do
  end subroutine populate_cosp_gridbox
#endif /* USE_COSP */


#ifdef USE_COSP
  subroutine populate_cosp_subcol(npoints, state, pbuf, emis, cld_swtau, cospstateIN, cospIN, snow_tau, snow_emis)    
    use physics_types,        only: physics_state
    use physics_buffer,       only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx, pbuf_get_index
    use constituents,         only: cnst_get_ind
    use rad_constituents,     only: rad_cnst_get_gas
    use interpolate_data,     only: lininterp_init,lininterp,lininterp_finish,interp_type
    use physconst,            only: pi, gravit
    use cam_history,          only: hist_fld_col_active 
    use cam_history_support,  only: max_fieldname_len
    use cmparray_mod,         only: CmpDayNite, ExpDayNite
    use mod_cosp_config,      only: R_UNDEF,parasol_nrefl, Nlvgrid, vgrid_zl, vgrid_zu
    use mod_cosp,             only: cosp_simulator
    use mod_quickbeam_optics, only: size_distribution

    ! ######################################################################################
    ! Inputs
    ! ######################################################################################
    integer, intent(in) :: npoints
    type(physics_state), intent(in),target  :: state
    type(physics_buffer_desc),      pointer :: pbuf(:)
    real(r8), intent(in) :: emis(pcols,pver)                  ! cloud longwave emissivity
    real(r8), intent(in) :: cld_swtau(pcols,pver)             ! RRTM cld_swtau, read in using this variable
    real(r8), intent(in),optional :: snow_tau(pcols,pver)  ! RRTM grid-box mean SW snow optical depth, used for CAM5 simulations 
    real(r8), intent(in),optional :: snow_emis(pcols,pver) ! RRTM grid-box mean LW snow optical depth, used for CAM5 simulations 

    ! ######################################################################################
    ! Local variables
    ! ######################################################################################
    integer :: lchnk                             ! chunk identifier
    integer :: ncol                              ! number of active atmospheric columns
    integer :: i, j, k, ix, iy, iz, itim_old
    
    ! Microphysics variables
    integer :: ixcldliq                                   ! cloud liquid amount index for state%q
    integer :: ixcldice                                   ! cloud ice amount index
    
    ! COSP-related local vars
    type(cosp_optical_inputs), intent(inout) :: cospIN                   ! COSP optical (or derived?) fields needed by simulators
    type(cosp_column_inputs), intent(inout)  :: cospstateIN              ! COSP model fields needed by simulators
    
    ! COSP input variables that depend on CAM
    logical :: use_reff = .true.                          ! True if effective radius to be used by radar simulator 
                                                          ! (always used by lidar)
    logical :: use_precipitation_fluxes = .true.          ! True if precipitation fluxes are input to the algorithm 
    real(r8), parameter :: emsfc_lw = 0.99_r8             ! longwave emissivity of surface at 10.5 microns 
    
    ! Local vars related to calculations to go from CAM input to COSP input
    ! cosp convective value includes both deep and shallow convection
    real(r8) :: mr_lsliq(npoints,pver)                     ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
    real(r8) :: mr_lsice(npoints,pver)                     ! mixing_ratio_large_scale_cloud_ice (kg/kg)
    real(r8) :: mr_ccliq(npoints,pver)                     ! mixing_ratio_convective_cloud_liquid (kg/kg)
    real(r8) :: mr_ccice(npoints,pver)                     ! mixing_ratio_convective_cloud_ice (kg/kg)
    real(r8) :: rain_cv(npoints,pverp)                     ! interface flux_convective_cloud_rain (kg m^-2 s^-1)
    real(r8) :: snow_cv(npoints,pverp)                     ! interface flux_convective_cloud_snow (kg m^-2 s^-1)
    real(r8) :: rain_cv_interp(npoints,pver)               ! midpoint flux_convective_cloud_rain (kg m^-2 s^-1)
    real(r8) :: snow_cv_interp(npoints,pver)               ! midpoint flux_convective_cloud_snow (kg m^-2 s^-1)
    real(r8) :: grpl_ls_interp(npoints,pver)               ! midpoint ls grp flux, should be 0
    real(r8) :: rain_ls_interp(npoints,pver)               ! midpoint ls rain flux (kg m^-2 s^-1)
    real(r8) :: snow_ls_interp(npoints,pver)               ! midpoint ls snow flux
    real(r8) :: reff_cosp(npoints,pver,nhydro)             ! effective radius for cosp input
    real(r8) :: dtau_s(npoints,pver)                       ! Optical depth of stratiform cloud at 0.67 um
    real(r8) :: dtau_c(npoints,pver)                       ! Optical depth of convective cloud at 0.67 um
    real(r8) :: dtau_s_snow(npoints,pver)                  ! Grid-box mean Optical depth of stratiform snow at 0.67 um
    real(r8) :: dem_s(npoints,pver)                        ! Longwave emis of stratiform cloud at 10.5 um
    real(r8) :: dem_c(npoints,pver)                        ! Longwave emis of convective cloud at 10.5 um
    real(r8) :: dem_s_snow(npoints,pver)                   ! Grid-box mean Optical depth of stratiform snow at 10.5 um
    real(r8) :: psfc(npoints)                              ! Surface pressure
    
    ! CAM pointers to get variables from radiation interface (get from rad_cnst_get_gas)
    real(r8), pointer, dimension(:,:) :: q               ! specific humidity (kg/kg)
    real(r8), pointer, dimension(:,:) :: o3              ! Mass mixing ratio 03
    
    ! CAM pointers to get variables from the physics buffer
    real(r8), pointer, dimension(:,:) :: cld             ! cloud fraction, tca - total_cloud_amount (0-1)
    real(r8), pointer, dimension(:,:) :: concld          ! concld fraction, cca - convective_cloud_amount (0-1)
    real(r8), pointer, dimension(:,:) :: rel             ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei             ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: ls_reffrain     ! rain effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: ls_reffsnow     ! snow effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cv_reffliq      ! convective cld liq effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: cv_reffice      ! convective cld ice effective drop size (microns)
    
    ! precip flux pointers (use for cam4 or cam5)
    ! Added pointers;  pbuff in zm_conv_intr.F90, calc in zm_conv.F90 
    real(r8), pointer, dimension(:,:) :: dp_flxprc       ! deep interface gbm flux_convective_cloud_rain+snow (kg m^-2 s^-1)
    real(r8), pointer, dimension(:,:) :: dp_flxsnw       ! deep interface gbm flux_convective_cloud_snow (kg m^-2 s^-1) 
    ! More pointers;  pbuf in convect_shallow.F90, calc in hk_conv.F90/convect_shallow.F90 (CAM4), uwshcu.F90 (CAM5)
    real(r8), pointer, dimension(:,:) :: sh_flxprc       ! shallow interface gbm flux_convective_cloud_rain+snow (kg m^-2 s^-1) 
    real(r8), pointer, dimension(:,:) :: sh_flxsnw       ! shallow interface gbm flux_convective_cloud_snow (kg m^-2 s^-1)
    ! More pointers;  pbuf in stratiform.F90, getting from pbuf here
    ! a) added as output to pcond subroutine in cldwat.F90 and to nmicro_pcond subroutine in cldwat2m_micro.F90
    real(r8), pointer, dimension(:,:) :: ls_flxprc       ! stratiform interface gbm flux_cloud_rain+snow (kg m^-2 s^-1) 
    real(r8), pointer, dimension(:,:) :: ls_flxsnw       ! stratiform interface gbm flux_cloud_snow (kg m^-2 s^-1)
    
    ! cloud mixing ratio pointers (note: large-scale in state)
    ! More pointers;  pbuf in convect_shallow.F90 (cam4) or stratiform.F90 (cam5)
    ! calc in hk_conv.F90 (CAM4 should be 0!), uwshcu.F90 but then affected by micro so values from stratiform.F90 (CAM5)
    real(r8), pointer, dimension(:,:) :: sh_cldliq       ! shallow gbm cloud liquid water (kg/kg)
    real(r8), pointer, dimension(:,:) :: sh_cldice       ! shallow gbm cloud ice water (kg/kg)
    ! More pointers;  pbuf in zm_conv_intr.F90, calc in zm_conv.F90, 0 for CAM4 and CAM5 (same convection scheme)
    real(r8), pointer, dimension(:,:) :: dp_cldliq       ! deep gbm cloud liquid water (kg/kg)
    real(r8), pointer, dimension(:,:) :: dp_cldice       ! deep gmb cloud ice water (kg/kg)
    
    type(interp_type)  :: interp_wgts
    integer, parameter :: extrap_method = 1              ! sets extrapolation method to boundary value (1)

    type(size_distribution) :: sd_wk       ! Work size distribution used by radar simulator
                                           ! This is to avoid directly passing  sd_cs(lchnk) to
                                           ! subsample_and_optics which would
                                           ! fail runtime if compiled more strictly, like in debug mode 
    real(r8), dimension(npoints,pver) :: cld_s, cld_c
    real(r8), pointer, dimension(:,:,:,:) :: &
       crm_cld, crm_dtau, crm_emis, crm_qc, crm_qi, crm_rel, crm_rei

    ncol = state%ncol
    lchnk = state%lchnk

    cospIN%emsfc_lw      = emsfc_lw

    ! State variables
    psfc = 0
    do iy = 1,crm_ny_rad
      do ix = 1,crm_nx_rad
        do i = 1,ncol
          j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
          psfc(j) = state%ps(i)
        end do
      end do
    end do

    ! Cloud fractions
    if (use_MMF) then
       call pbuf_get_field(pbuf, pbuf_get_index('CRM_CLD_RAD'), crm_cld)
       cld_s = 0
       cld_c = 0
       do iz = 1,crm_nz
         do iy = 1,crm_ny_rad
           do ix = 1,crm_nx_rad
             do i = 1,ncol
               j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
               k = pver - iz + 1
               cld_s(j,k) = crm_cld(i,ix,iy,iz)
             end do
           end do
         end do
       end do
    else
       itim_old = pbuf_old_tim_idx()
       call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
       call pbuf_get_field(pbuf, concld_idx, concld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
       do k = 1,pver
         do i = 1,ncol
           cld_s(i,k) = cld(i,k)
           cld_c(i,k) = concld(i,k)
         end do
       end do
    end if  ! use_MMF   

    ! precipitation fluxes (use for both cam4 and cam5 for now....)
    if (use_MMF) then
       use_precipitation_fluxes = .false.
       snow_ls_interp = 0
       grpl_ls_interp = 0
       rain_cv_interp = 0
       snow_cv_interp = 0
       do iz = 1,crm_nz
         do iy = 1,crm_ny_rad
           do ix = 1,crm_nx_rad
             do i = 1,ncol
               j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
               k = pver - iz + 1
               rain_ls_interp(j,k) = 0._r8  ! TODO: fix this
               snow_ls_interp(j,k) = 0._r8  ! TODO: fix this
               grpl_ls_interp(j,k) = 0._r8  ! TODO: fix this
             end do
           end do
         end do
       end do
    else
       call pbuf_get_field(pbuf, dpflxprc_idx, dp_flxprc  )
       call pbuf_get_field(pbuf, dpflxsnw_idx, dp_flxsnw  )
       call pbuf_get_field(pbuf, shflxprc_idx, sh_flxprc  )
       call pbuf_get_field(pbuf, shflxsnw_idx, sh_flxsnw  )
       call pbuf_get_field(pbuf, lsflxprc_idx, ls_flxprc  )
       call pbuf_get_field(pbuf, lsflxsnw_idx, ls_flxsnw  )
      
       ! add together deep and shallow convection precipitation fluxes, recall *_flxprc variables are rain+snow
       rain_cv(1:ncol,1:pverp) = (sh_flxprc(1:ncol,1:pverp)-sh_flxsnw(1:ncol,1:pverp)) + &
            (dp_flxprc(1:ncol,1:pverp)-dp_flxsnw(1:ncol,1:pverp))
       snow_cv(1:ncol,1:pverp) = sh_flxsnw(1:ncol,1:pverp) + dp_flxsnw(1:ncol,1:pverp)
       
       ! All precip fluxes in COSP should be mid points, all values are grid-box mean ("gbm") (Yuying)
       ! Interpolate interface precip fluxes to mid points
       grpl_ls_interp(1:ncol,1:pver) = 0._r8  ! NOTE: graupel remains zero
       rain_ls_interp(1:ncol,1:pver) = 0._r8 
       snow_ls_interp(1:ncol,1:pver) = 0._r8
       do i = 1,ncol
          ! Find intepolation weights
          call lininterp_init(state%zi(i,1:pverp),pverp,state%zm(i,1:pver),pver,extrap_method,interp_wgts)
          ! Do interpolation
          call lininterp(rain_cv(i,1:pverp),pverp,rain_cv_interp(i,1:pver),pver,interp_wgts)
          call lininterp(snow_cv(i,1:pverp),pverp,snow_cv_interp(i,1:pver),pver,interp_wgts)
          call lininterp(ls_flxprc(i,1:pverp),pverp,rain_ls_interp(i,1:pver),pver,interp_wgts)
          call lininterp(ls_flxsnw(i,1:pverp),pverp,snow_ls_interp(i,1:pver),pver,interp_wgts)
          call lininterp_finish(interp_wgts)
          ! ls_flxprc is for rain+snow, find rain_ls_interp by subtracting off snow_ls_interp
          rain_ls_interp(i,1:pver)=rain_ls_interp(i,1:pver)-snow_ls_interp(i,1:pver)
       end do
       ! Make sure interpolated values are not less than 0 - COSP was complaining and resetting small 
       ! negative values to zero.
       ! ----- WARNING: COSP_CHECK_INPUT_2D: minimum value of rain_ls set to:      0.000000000000000 
       ! So I set negative values to zero here... 
       do k=1,pver
          do i=1,ncol
             if (rain_ls_interp(i,k) .lt. 0._r8) then
                rain_ls_interp(i,k)=0._r8
             end if
             if (snow_ls_interp(i,k) .lt. 0._r8) then
                snow_ls_interp(i,k)=0._r8
             end if
             if (rain_cv_interp(i,k) .lt. 0._r8) then
                rain_cv_interp(i,k)=0._r8
             end if
             if (snow_cv_interp(i,k) .lt. 0._r8) then
                snow_cv_interp(i,k)=0._r8
             end if
          end do
       end do
    end if  ! use_MMF

    ! Hydrometeor mixing ratios
    if (use_MMF) then
       ! TODO: convert mr_lsliq+mr_lsice to array like mrhydro with enum-like
       ! indices; handle preciptation mixing ratios here
       call pbuf_get_field(pbuf, pbuf_get_index('CRM_QC_RAD'), crm_qc)
       call pbuf_get_field(pbuf, pbuf_get_index('CRM_QI_RAD'), crm_qi)
       mr_lsliq = 0
       mr_lsice = 0
       mr_ccliq = 0
       mr_ccice = 0
       do iz = 1,crm_nz
         do iy = 1,crm_ny_rad
           do ix = 1,crm_nx_rad
             do i = 1,ncol
               j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
               k = pver - iz + 1
               ! Mixing ratios
               mr_lsliq(j,k) = crm_qc(i,ix,iy,iz)
               mr_lsice(j,k) = crm_qi(i,ix,iy,iz)
               mr_ccliq(j,k) = 0._r8
               mr_ccice(j,k) = 0._r8
             end do
           end do
         end do
       end do
    else
       ! CAM5 cloud mixing ratio calculations
       ! Note: Although CAM5 has non-zero convective cloud mixing ratios that affect the model state, 
       ! Convective cloud water is NOT part of radiation calculations.
       ! Get indices to radiative constituents
       call cnst_get_ind('CLDLIQ',ixcldliq)
       call cnst_get_ind('CLDICE',ixcldice)
       ! convective cloud mixing ratios (use for cam4 and cam5)
       call pbuf_get_field(pbuf, dpcldliq_idx, dp_cldliq  )
       call pbuf_get_field(pbuf, dpcldice_idx, dp_cldice  )
       ! get from pbuf in stratiform.F90
       call pbuf_get_field(pbuf, shcldliq1_idx, sh_cldliq  )
       call pbuf_get_field(pbuf, shcldice1_idx, sh_cldice  )
       mr_ccliq(1:ncol,1:pver)           = 0._r8
       mr_ccice(1:ncol,1:pver)           = 0._r8
       mr_lsliq(1:ncol,1:pver)           = 0._r8
       mr_lsice(1:ncol,1:pver)           = 0._r8
       do k=1,pver
          do i=1,ncol
             if (cld(i,k) .gt. 0._r8) then
                ! note: convective mixing ratio is the sum of shallow and deep convective clouds in CAM5
                mr_ccliq(i,k) = sh_cldliq(i,k) + dp_cldliq(i,k)
                mr_ccice(i,k) = sh_cldice(i,k) + dp_cldice(i,k)
                mr_lsliq(i,k) = state%q(i,k,ixcldliq)   ! state only includes stratiform (kg/kg)  
                mr_lsice(i,k) = state%q(i,k,ixcldice)   ! state only includes stratiform (kg/kg)
             else
                mr_ccliq(i,k) = 0._r8
                mr_ccice(i,k) = 0._r8
                mr_lsliq(i,k) = 0._r8
                mr_lsice(i,k) = 0._r8
             end if
          end do
       end do
    end if  ! use_MMF

    ! Effective radii
    if (use_MMF) then
       call pbuf_get_field(pbuf, pbuf_get_index('CRM_REL'), crm_rel)
       call pbuf_get_field(pbuf, pbuf_get_index('CRM_REI'), crm_rei)
       do iz = 1,crm_nz
         do iy = 1,crm_ny_rad
           do ix = 1,crm_nx_rad
             do i = 1,ncol
               j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
               k = pver - iz + 1
               reff_cosp(j,k,I_LSCLIQ) = crm_rel(i,ix,iy,iz) * 1e-6_r8  ! microns to meters
               reff_cosp(j,k,I_LSCICE) = crm_rei(i,ix,iy,iz) * 1e-6_r8  ! microns to meters
               reff_cosp(j,k,I_LSRAIN) = 0._r8
               reff_cosp(j,k,I_LSSNOW) = 0._r8
               reff_cosp(j,k,I_CVCLIQ) = 0._r8
               reff_cosp(j,k,I_CVCICE) = 0._r8
               reff_cosp(j,k,I_CVRAIN) = 0._r8
               reff_cosp(j,k,I_CVSNOW) = 0._r8
               reff_cosp(j,k,I_LSGRPL) = 0._r8
             end do
           end do
         end do
       end do
    else
       ! Previously, I had set use_reff=.false.
       ! use_reff = .false.  ! if you use this,all sizes use DEFAULT_LIDAR_REFF = 30.0e-6 meters
       ! The specification of reff_cosp now follows e-mail discussion with Yuying in January 2011. (see above)
       ! All of the values that I have assembled in the code are in microns...convert to meters since that is what COSP wants.
       call pbuf_get_field(pbuf, rel_idx, rel  )
       call pbuf_get_field(pbuf, rei_idx, rei)
       ! added some more sizes to physics buffer in stratiform.F90 for COSP inputs
       call pbuf_get_field(pbuf, lsreffrain_idx, ls_reffrain  )
       call pbuf_get_field(pbuf, lsreffsnow_idx, ls_reffsnow  )
       call pbuf_get_field(pbuf, cvreffliq_idx,  cv_reffliq   )
       call pbuf_get_field(pbuf, cvreffice_idx,  cv_reffice   )
       reff_cosp(1:ncol,1:pver,1:nhydro) = 0._r8
       reff_cosp(1:ncol,1:pver,I_LSCLIQ) = rel(1:ncol,1:pver)*1.e-6_r8          ! (same as effc and effliq in stratiform.F90)
       reff_cosp(1:ncol,1:pver,I_LSCICE) = rei(1:ncol,1:pver)*1.e-6_r8          ! (same as effi and effice in stratiform.F90)
       reff_cosp(1:ncol,1:pver,I_LSRAIN) = ls_reffrain(1:ncol,1:pver)*1.e-6_r8  ! (calculated in cldwat2m_micro.F90, passed to stratiform.F90)
       reff_cosp(1:ncol,1:pver,I_LSSNOW) = ls_reffsnow(1:ncol,1:pver)*1.e-6_r8  ! (calculated in cldwat2m_micro.F90, passed to stratiform.F90)
       reff_cosp(1:ncol,1:pver,I_CVCLIQ) = cv_reffliq(1:ncol,1:pver)*1.e-6_r8   ! (calculated in stratiform.F90, not actually used in radiation)
       reff_cosp(1:ncol,1:pver,I_CVCICE) = cv_reffice(1:ncol,1:pver)*1.e-6_r8   ! (calculated in stratiform.F90, not actually used in radiation)
       reff_cosp(1:ncol,1:pver,I_CVRAIN) = ls_reffrain(1:ncol,1:pver)*1.e-6_r8  ! (same as stratiform per Andrew)
       reff_cosp(1:ncol,1:pver,I_CVSNOW) = ls_reffsnow(1:ncol,1:pver)*1.e-6_r8  ! (same as stratiform per Andrew)
       reff_cosp(1:ncol,1:pver,I_LSGRPL) = 0._r8                                ! (using radar default reff)
    end if  ! use_MMF

    ! Set optical depth and emissivity
    if (use_MMF) then
       call pbuf_get_field(pbuf, pbuf_get_index('CRM_EMIS'), crm_emis)
       call pbuf_get_field(pbuf, pbuf_get_index('CRM_DTAU'), crm_dtau)
       dtau_s = 0
       dtau_c = 0
       dem_s = 0
       dem_c = 0
       do iz = 1,crm_nz
         do iy = 1,crm_ny_rad
           do ix = 1,crm_nx_rad
             do i = 1,ncol
               j = _IDX321(i, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
               k = pver - iz + 1
               dtau_s(j,k) = crm_dtau(i,ix,iy,iz)
               dtau_c(j,k) = crm_dtau(i,ix,iy,iz)
               dem_s (j,k) = crm_emis(i,ix,iy,iz)
               dem_c (j,k) = crm_emis(i,ix,iy,iz)
             end do
           end do
         end do
       end do
    else
       ! NOTES:
       ! 1) EAM assumes same radiative properties for stratiform and convective clouds, 
       ! 2) COSP wants in-cloud values. EAM values of cld_swtau are in-cloud means.
       ! 3) snow_tau and snow_emis are passed without modification to COSP
       dtau_s(1:ncol,1:pver)      = cld_swtau(1:ncol,1:pver)     ! 0.67 micron optical depth of stratiform (in-cloud)
       dtau_c(1:ncol,1:pver)      = cld_swtau(1:ncol,1:pver)     ! 0.67 micron optical depth of convective (in-cloud)
       dem_s(1:ncol,1:pver)       = emis(1:ncol,1:pver)          ! 10.5 micron longwave emissivity of stratiform (in-cloud)
       dem_c(1:ncol,1:pver)       = emis(1:ncol,1:pver)          ! 10.5 micron longwave emissivity of convective (in-cloud)
       if (present(snow_tau) .and. present(snow_emis)) then
          dem_s_snow(1:ncol,1:pver)  = snow_emis(1:ncol,1:pver)  ! 10.5 micron grid-box mean optical depth of stratiform snow
          dtau_s_snow(1:ncol,1:pver) = snow_tau(1:ncol,1:pver)   ! 0.67 micron grid-box mean optical depth of stratiform snow
       else
          dem_s_snow(1:ncol,1:pver) = 0._r8
          dtau_s_snow(1:ncol,1:pver) = 0._r8
       end if
    end if  ! use_MMF

    if (lradar_sim) then 
       cospIN%rcfg_cloudsat = rcfg_cs(lchnk)
       sd_wk = sd_cs(lchnk)
    end if

    call subsample_and_optics(npoints,pver,nscol_cosp,nhydro,overlap,           &
         use_precipitation_fluxes,lidar_ice_type,sd_wk,cld_s(1:npoints,1:pver), &
         cld_c(1:npoints,1:pver),rain_ls_interp(1:npoints,1:pver),              &
         snow_ls_interp(1:npoints,1:pver),grpl_ls_interp(1:npoints,1:pver),     &
         rain_cv_interp(1:npoints,1:pver),snow_cv_interp(1:npoints,1:pver),     &
         mr_lsliq(1:npoints,1:pver),mr_lsice(1:npoints,1:pver),                 &
         mr_ccliq(1:npoints,1:pver),mr_ccice(1:npoints,1:pver),                 &
         reff_cosp(1:npoints,1:pver,:),dtau_c(1:npoints,1:pver),                &
         dtau_s(1:npoints,1:pver),dem_c(1:npoints,1:pver),                      &
         dem_s(1:npoints,1:pver),dtau_s_snow(1:npoints,1:pver),                 &
         dem_s_snow(1:npoints,1:pver),psfc(1:npoints),cospstateIN,cospIN)
    if (lradar_sim) sd_cs(lchnk) = sd_wk

   end subroutine populate_cosp_subcol
#endif /* USE_COSP */
 

#ifdef USE_COSP
   ! Remask passive simulator output after call to COSP. This should NOT be necessary, and if it is
   ! we need to look into the version of COSP we are using and probably update. I've pushed some
   ! PRs to the COSP repo to fix this stuff in the past, so I think we are probably using a really
   ! old version of COSP here that does not include those fixes.
   ! TODO: revisit this
   subroutine cosp_remask_passive(cospstateIN, cospOUT)
      use mod_cosp, only: cosp_column_inputs, cosp_outputs
      type(cosp_column_inputs), intent(in) :: cospstateIN
      type(cosp_outputs), intent(in) :: cospOUT
      integer :: ncol, i, k

      ncol = size(cospstateIN%sunlit)

      ! ISCCP simulator
      if (lisccp_sim) then
         where(cospstateIN%sunlit(1:ncol) .eq. 0)
            cospOUT%isccp_totalcldarea(1:ncol)  = R_UNDEF
            cospOUT%isccp_meanptop(1:ncol)      = R_UNDEF
            cospOUT%isccp_meantaucld(1:ncol)    = R_UNDEF
            cospOUT%isccp_meanalbedocld(1:ncol) = R_UNDEF
            cospOUT%isccp_meantb(1:ncol)        = R_UNDEF
            cospOUT%isccp_meantbclr(1:ncol)     = R_UNDEF
         end where
         do i=1,nscol_cosp
            where (cospstateIN%sunlit(1:ncol) .eq. 0)
               cospOUT%isccp_boxtau(1:ncol,i)  = R_UNDEF
               cospOUT%isccp_boxptop(1:ncol,i) = R_UNDEF
            end where
         end do
         do i=1,nprs_cosp
            do k=1,ntau_cosp
               where(cospstateIN%sunlit(1:ncol) .eq. 0)
                  cospOUT%isccp_fq(1:ncol,k,i) = R_UNDEF
               end where
            end do
         end do
      end if

      ! MISR simulator
      if (lmisr_sim) then
         do i=1,nhtmisr_cosp
            do k=1,ntau_cosp
               where(cospstateIN%sunlit(1:ncol) .eq. 0)
                  cospOUT%misr_fq(1:ncol,k,i) = R_UNDEF
               end where
            end do
         end do
      end if

      ! MODIS simulator
      if (lmodis_sim) then
         where(cospstateIN%sunlit(1:ncol) .eq. 0)
            cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol)       = R_UNDEF
            cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol)       = R_UNDEF
            cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol)         = R_UNDEF
            cospOUT%modis_Cloud_Fraction_High_Mean(1:ncol)        = R_UNDEF
            cospOUT%modis_Cloud_Fraction_Mid_Mean(1:ncol)         = R_UNDEF
            cospOUT%modis_Cloud_Fraction_Low_Mean(1:ncol)         = R_UNDEF
            cospOUT%modis_Optical_Thickness_Total_Mean(1:ncol)    = R_UNDEF
            cospOUT%modis_Optical_Thickness_Water_Mean(1:ncol)    = R_UNDEF
            cospOUT%modis_Optical_Thickness_Ice_Mean(1:ncol)      = R_UNDEF
            cospOUT%modis_Optical_Thickness_Total_LogMean(1:ncol) = R_UNDEF
            cospOUT%modis_Optical_Thickness_Water_LogMean(1:ncol) = R_UNDEF
            cospOUT%modis_Optical_Thickness_Ice_LogMean(1:ncol)   = R_UNDEF
            cospOUT%modis_Cloud_Particle_Size_Water_Mean(1:ncol)  = R_UNDEF
            cospOUT%modis_Cloud_Particle_Size_Ice_Mean(1:ncol)    = R_UNDEF
            cospOUT%modis_Cloud_Top_Pressure_Total_Mean(1:ncol)   = R_UNDEF
            cospOUT%modis_Liquid_Water_Path_Mean(1:ncol)          = R_UNDEF
            cospOUT%modis_Ice_Water_Path_Mean(1:ncol)             = R_UNDEF
         end where
         do i=1,ntau_cosp_modis
            do k=1,nprs_cosp
               where(cospstateIN%sunlit(1:ncol) .eq. 0)
                  cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(1:ncol,i,k) = R_UNDEF 
               end where
            end do
            do k=1,numMODISReffIceBins
               where(cospstateIN%sunlit(1:ncol) .eq. 0)
                  cospOUT%modis_Optical_Thickness_vs_ReffICE(1:ncol,i,k) = R_UNDEF
               end where
            end do
            do k=1,numMODISReffLiqBins
               where(cospstateIN%sunlit(1:ncol) .eq. 0)
                  cospOUT%modis_Optical_Thickness_vs_ReffLIQ(1:ncol,i,k) = R_UNDEF
               end where
            end do
         end do
      end if
   end subroutine cosp_remask_passive
#endif

#ifdef USE_COSP
   subroutine cosp_write_outputs(state, cospIN, cospOUT)
      use physics_types, only: physics_state
      use mod_cosp, only: cosp_optical_inputs, cosp_outputs
      use cam_history, only: outfld

      type(physics_state)      , intent(in)    :: state
      type(cosp_optical_inputs), intent(in)    :: cospIN
      type(cosp_outputs)       , intent(inout) :: cospOUT  ! inout so we can fix/weight before writing
      integer :: ncol, lchnk
      real(r8) :: zeros(pcols)

      ! Number of columns and chunk ID
      ncol = state%ncol
      lchnk = state%lchnk

      ! ISCCP OUTPUTS
      if (lisccp_sim) then
         call outfld('FISCCP1_COSP', cospOUT%isccp_fq(1:ncol,1:ntau_cosp,1:nprs_cosp), ncol, lchnk)
         call outfld('CLDTOT_ISCCP', cospOUT%isccp_totalcldarea(1:ncol), ncol, lchnk)
         ! Weight outputs by cloud fraction
         where (cospOUT%isccp_totalcldarea(:ncol) .eq. R_UNDEF)
            cospOUT%isccp_meanalbedocld(:ncol) = R_UNDEF
            cospOUT%isccp_meanptop(:ncol)      = R_UNDEF
            cospOUT%isccp_meantaucld(:ncol)    = R_UNDEF
         elsewhere
            cospOUT%isccp_meanalbedocld(:ncol) = cospOUT%isccp_meanalbedocld(:ncol)*cospOUT%isccp_totalcldarea(:ncol)
            cospOUT%isccp_meanptop(:ncol)      = cospOUT%isccp_meanptop(:ncol)*cospOUT%isccp_totalcldarea(:ncol)
            cospOUT%isccp_meantaucld(:ncol)    = cospOUT%isccp_meantaucld(:ncol)*cospOUT%isccp_totalcldarea(:ncol)
         end where
         call outfld('MEANCLDALB_ISCCP', cospOUT%isccp_meanalbedocld(1:ncol), ncol, lchnk)
         call outfld('MEANPTOP_ISCCP'  , cospOUT%isccp_meanptop(1:ncol)     , ncol, lchnk)
         call outfld('MEANTAU_ISCCP'   , cospOUT%isccp_meantaucld(1:ncol)   , ncol, lchnk)
         call outfld('MEANTB_ISCCP'    , cospOUT%isccp_meantb(1:ncol)       , ncol, lchnk)
         call outfld('MEANTBCLR_ISCCP' , cospOUT%isccp_meantbclr(1:ncol)    , ncol, lchnk)
      end if
      
      ! CALIPSO SIMULATOR OUTPUTS
      if (llidar_sim) then
         ! Set missing values to 0 (clear air) to avoid issues with check_acccum (likely levels below sea level)
         call replace_values2d(cospOUT%calipso_lidarcld     (:ncol,:nht_cosp)  , R_UNDEF, 0._r8)
         call replace_values3d(cospOUT%calipso_cfad_sr      (:ncol,:,:)        , R_UNDEF, 0._r8)
         call replace_values2d(cospOUT%parasolGrid_refl     (:ncol,:nsza_cosp) , R_UNDEF, 0._r8)
         call replace_values3d(cospOUT%calipso_lidarcldphase(:ncol,:nht_cosp,:), R_UNDEF, 0._r8)
         call replace_values3d(cospOUT%calipso_lidarcldtmp  (:ncol,:nht_cosp,:), R_UNDEF, 0._r8)
         ! Send outputs to history buffer
         call outfld('CLDLOW_CAL'    , cospOUT%calipso_cldlayer(1:ncol,1)       , ncol, lchnk)
         call outfld('CLDMED_CAL'    , cospOUT%calipso_cldlayer(1:ncol,2)       , ncol, lchnk)
         call outfld('CLDHGH_CAL'    , cospOUT%calipso_cldlayer(1:ncol,3)       , ncol, lchnk)
         call outfld('CLDTOT_CAL'    , cospOUT%calipso_cldlayer(1:ncol,4)       , ncol, lchnk)
         call outfld('CLDTOT_CAL_ICE', cospOUT%calipso_cldlayerphase(1:ncol,4,1), ncol, lchnk)
         call outfld('CLDTOT_CAL_LIQ', cospOUT%calipso_cldlayerphase(1:ncol,4,2), ncol, lchnk)
         call outfld('CLDTOT_CAL_UN' , cospOUT%calipso_cldlayerphase(1:ncol,4,3), ncol, lchnk)
         call outfld('CLDHGH_CAL_ICE', cospOUT%calipso_cldlayerphase(1:ncol,3,1), ncol, lchnk)
         call outfld('CLDHGH_CAL_LIQ', cospOUT%calipso_cldlayerphase(1:ncol,3,2), ncol, lchnk)
         call outfld('CLDHGH_CAL_UN' , cospOUT%calipso_cldlayerphase(1:ncol,3,3), ncol, lchnk)
         call outfld('CLDMED_CAL_ICE', cospOUT%calipso_cldlayerphase(1:ncol,2,1), ncol, lchnk)
         call outfld('CLDMED_CAL_LIQ', cospOUT%calipso_cldlayerphase(1:ncol,2,2), ncol, lchnk)
         call outfld('CLDMED_CAL_UN' , cospOUT%calipso_cldlayerphase(1:ncol,2,3), ncol, lchnk)
         call outfld('CLDLOW_CAL_ICE', cospOUT%calipso_cldlayerphase(1:ncol,1,1), ncol, lchnk)
         call outfld('CLDLOW_CAL_LIQ', cospOUT%calipso_cldlayerphase(1:ncol,1,2), ncol, lchnk)
         call outfld('CLDLOW_CAL_UN' , cospOUT%calipso_cldlayerphase(1:ncol,1,3), ncol, lchnk)
         call outfld('CLD_CAL'       , cospOUT%calipso_lidarcld     (1:ncol,1:nht_cosp)  , ncol,lchnk)
         call outfld('MOL532_CAL'    , cospOUT%calipso_beta_mol     (1:ncol,1:nhtml_cosp), ncol, lchnk)
         call outfld('CFAD_SR532_CAL', cospOUT%calipso_cfad_sr      (1:ncol,1:nsr_cosp,1:nht_cosp),ncol,lchnk)
         call outfld('RFL_PARASOL'   , cospOUT%parasolGrid_refl     (1:ncol,1:nsza_cosp) , ncol, lchnk)
         call outfld('CLD_CAL_ICE'   , cospOUT%calipso_lidarcldphase(1:ncol,1:nht_cosp,1), ncol, lchnk)
         call outfld('CLD_CAL_LIQ'   , cospOUT%calipso_lidarcldphase(1:ncol,1:nht_cosp,2), ncol, lchnk)
         call outfld('CLD_CAL_UN'    , cospOUT%calipso_lidarcldphase(1:ncol,1:nht_cosp,3), ncol, lchnk)
         call outfld('CLD_CAL_TMP'   , cospOUT%calipso_lidarcldtmp  (1:ncol,1:nht_cosp,1), ncol, lchnk)
         call outfld('CLD_CAL_TMPUN' , cospOUT%calipso_lidarcldtmp  (1:ncol,1:nht_cosp,4), ncol, lchnk)
      end if
      
      ! RADAR SIMULATOR OUTPUTS
      if (lradar_sim) then
        
         ! Output the mixing-ratio for all hydrometeor types in Cloudsat near-surface precipitation diagnostics
         ! *NOTE* These fields are simply the native CAM mixing-ratios for each hydrometeor type used in the 
         !        CAM6 microphysics scheme, interpolated to the same vertical grid used by the Cloudsat
         !        simulator. These fields are not part of the radar simulator standard output, as these fields
         !        are entirely dependent on the host models microphysics, not the retrieval.

         ! fails check_accum if this is set... with ht_cosp set relative to sea level, mix of R_UNDEF and realvalue 
         call replace_values3d(cospOUT%cloudsat_cfad_ze(:ncol,:,:), R_UNDEF, 0._r8)
         call outfld('CFAD_DBZE94_CS',cospOUT%cloudsat_cfad_ze(:ncol,:,:), ncol, lchnk)
         call outfld('CLDTOT_CALCS',  cospOUT%radar_lidar_tcc(1:ncol), ncol, lchnk)
         zeros(1:ncol)  = 0._r8
         call outfld('CLDTOT_CS',     zeros(1:ncol), ncol, lchnk)  ! No longer in COSP; TODO: remove
         call outfld('CLDTOT_CS2',    zeros(1:ncol), ncol, lchnk)  ! No longer in COSP; TODO: remove
         call outfld('CLD_CAL_NOTCS', cospOUT%lidar_only_freq_cloud(1:ncol,1:nht_cosp), ncol, lchnk)

         ! Cloudsat near-surface precipitation diagnostics
         call outfld('CS_NOPRECIP', cospOUT%cloudsat_precip_cover(1:ncol,1 ),  ncol, lchnk)
         call outfld('CS_RAINPOSS', cospOUT%cloudsat_precip_cover(1:ncol,2 ),  ncol, lchnk)
         call outfld('CS_RAINPROB', cospOUT%cloudsat_precip_cover(1:ncol,3 ),  ncol, lchnk)
         call outfld('CS_RAINCERT', cospOUT%cloudsat_precip_cover(1:ncol,4 ),  ncol, lchnk)
         call outfld('CS_SNOWPOSS', cospOUT%cloudsat_precip_cover(1:ncol,5 ),  ncol, lchnk)
         call outfld('CS_SNOWCERT', cospOUT%cloudsat_precip_cover(1:ncol,6 ),  ncol, lchnk)
         call outfld('CS_MIXPOSS' , cospOUT%cloudsat_precip_cover(1:ncol,7 ),  ncol, lchnk)
         call outfld('CS_MIXCERT' , cospOUT%cloudsat_precip_cover(1:ncol,8 ),  ncol, lchnk)
         call outfld('CS_RAINHARD', cospOUT%cloudsat_precip_cover(1:ncol,9 ),  ncol, lchnk)
         call outfld('CS_UN'      , cospOUT%cloudsat_precip_cover(1:ncol,10),  ncol, lchnk)
         call outfld('CS_PIA'     , cospOUT%cloudsat_pia(1:ncol)            ,  ncol, lchnk)
      end if
      
      ! MISR SIMULATOR OUTPUTS
      if (lmisr_sim) then
         call outfld('CLD_MISR', cospOUT%misr_fq(1:ncol,1:ntau_cosp,1:nhtmisr_cosp), ncol, lchnk)
      end if
      
      ! MODIS SIMULATOR OUTPUTS
      if (lmodis_sim) then

         call outfld('CLTMODIS',cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol), ncol, lchnk)
         call outfld('CLWMODIS',cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol), ncol, lchnk)
         call outfld('CLIMODIS',cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol)  , ncol, lchnk)
         call outfld('CLHMODIS',cospOUT%modis_Cloud_Fraction_High_Mean(1:ncol) , ncol, lchnk)
         call outfld('CLMMODIS',cospOUT%modis_Cloud_Fraction_Mid_Mean(1:ncol)  , ncol, lchnk)
         call outfld('CLLMODIS',cospOUT%modis_Cloud_Fraction_Low_Mean(1:ncol)  , ncol, lchnk)
         
         ! where there is no cloud fraction or no retrieval, set to R_UNDEF, 
         ! otherwise weight retrievals by cloud fraction
         call weight_output(cospOUT%modis_Optical_Thickness_Total_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Optical_Thickness_Water_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Optical_Thickness_Ice_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Optical_Thickness_Total_LogMean(1:ncol), cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Optical_Thickness_Water_LogMean(1:ncol), cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Optical_Thickness_Ice_LogMean(1:ncol), cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Cloud_Particle_Size_Water_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Cloud_Particle_Size_Ice_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Cloud_Top_Pressure_Total_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Liquid_Water_Path_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol), R_UNDEF)
         call weight_output(cospOUT%modis_Ice_Water_Path_Mean(1:ncol), cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol), R_UNDEF)
         call outfld('TAUTMODIS'   , cospOUT%modis_Optical_Thickness_Total_Mean(1:ncol)   , ncol, lchnk)
         call outfld('TAUWMODIS'   , cospOUT%modis_Optical_Thickness_Water_Mean(1:ncol)   , ncol, lchnk)
         call outfld('TAUIMODIS'   , cospOUT%modis_Optical_Thickness_Ice_Mean(1:ncol)     , ncol, lchnk)
         call outfld('TAUTLOGMODIS', cospOUT%modis_Optical_Thickness_Total_LogMean(1:ncol), ncol, lchnk)
         call outfld('TAUWLOGMODIS', cospOUT%modis_Optical_Thickness_Water_LogMean(1:ncol), ncol, lchnk)
         call outfld('TAUILOGMODIS', cospOUT%modis_Optical_Thickness_Ice_LogMean(1:ncol)  , ncol, lchnk)
         call outfld('REFFCLWMODIS', cospOUT%modis_Cloud_Particle_Size_Water_Mean(1:ncol) , ncol, lchnk)
         call outfld('REFFCLIMODIS', cospOUT%modis_Cloud_Particle_Size_Ice_Mean(1:ncol)   , ncol, lchnk)
         call outfld('PCTMODIS'    , cospOUT%modis_Cloud_Top_Pressure_Total_Mean(1:ncol)  , ncol, lchnk)
         call outfld('LWPMODIS'    , cospOUT%modis_Liquid_Water_Path_Mean(1:ncol)         , ncol, lchnk)
         call outfld('IWPMODIS'    , cospOUT%modis_Ice_Water_Path_Mean(1:ncol)            , ncol, lchnk)
         call outfld('CLMODIS'     , cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(1:ncol,1:ntau_cosp_modis,1:nprs_cosp),ncol,lchnk) 
         call outfld('CLRIMODIS'   , cospOUT%modis_Optical_Thickness_vs_ReffICE(1:ncol,1:ntau_cosp_modis,1:numMODISReffIceBins), ncol, lchnk) 
         call outfld('CLRLMODIS'   , cospOUT%modis_Optical_Thickness_vs_ReffLIQ(1:ncol,1:ntau_cosp_modis,1:numMODISReffLiqBins), ncol, lchnk) 
      end if
      
      ! SUB-COLUMN OUTPUT
      if (lfrac_out) then
         call outfld('SCOPS_OUT',cospIN%frac_out(1:ncol,1:nscol_cosp,1:nhtml_cosp),ncol,lchnk)
         if (lisccp_sim) then
            call outfld('TAU_ISCCP', cospOUT%isccp_boxtau(1:ncol,1:nscol_cosp), ncol, lchnk)
            call outfld('CLDPTOP_ISCCP', cospOUT%isccp_boxptop(1:ncol,1:nscol_cosp), ncol, lchnk)
         end if
         if (llidar_sim) then
            call outfld('ATB532_CAL',cospOUT%calipso_beta_tot(1:ncol,1:nscol_cosp,1:nhtml_cosp),ncol,lchnk)
         end if
         if (lradar_sim) then
            call outfld('DBZE_CS',cospOUT%cloudsat_Ze_tot(1:ncol,:,:), ncol, lchnk)
         end if
      end if
   end subroutine cosp_write_outputs
#endif

#ifdef USE_COSP
    subroutine packed_average1d(nc, nx, ny, d_in, d_out)
      integer, intent(in) :: nc, nx, ny
      real(r8), intent(in) :: d_in(:)
      real(r8), intent(out) :: d_out(:)
      integer :: num_valid
      integer :: ip, ic, ix, iy
      call assert(size(d_in , 1) == nc * nx * ny, 'np /= nc * nx * ny')
      d_out = 0._r8
      do ic = 1,nc
        num_valid = 0
        do iy = 1,ny
          do ix = 1,nx
            ip = _IDX321(ic,ix,iy,nc,nx,ny)
            if (d_in(ip) /= R_UNDEF) then
              d_out(ic) = d_out(ic) + d_in(ip)
              num_valid = num_valid + 1
            end if
          end do
        end do
        ! If no valid values found, we need to set this column to fill value
        if (num_valid > 0) then
           d_out(ic) = d_out(ic) / num_valid
        else
           d_out(ic) = R_UNDEF
        end if
      end do
    end subroutine
    subroutine packed_average2d(nc, nx, ny, d_in, d_out)
      integer, intent(in) :: nc, nx, ny
      real(r8), intent(in) :: d_in(:,:)
      real(r8), intent(out) :: d_out(:,:)
      integer :: num_valid
      integer :: ip, ic, ix, iy, i2, i3
      call assert(size(d_in , 1) == nc * nx * ny, 'np /= nc * nx * ny')
      d_out = 0._r8
      do i2 = 1,size(d_in,2)
        do ic = 1,nc
          num_valid = 0
          do iy = 1,ny
            do ix = 1,nx
              ip = _IDX321(ic,ix,iy,nc,nx,ny)
              if (d_in(ip,i2) /= R_UNDEF) then
                d_out(ic,i2) = d_out(ic,i2) + d_in(ip,i2)
                num_valid = num_valid + 1
              end if
            end do
          end do
          ! If no valid values found, need to set to fill value
          if (num_valid > 0) then
            d_out(ic,i2) = d_out(ic,i2) / num_valid
          else
            d_out(ic,i2) = R_UNDEF
          end if
        end do
      end do
    end subroutine
    subroutine packed_average3d(nc, nx, ny, d_in, d_out)
      integer, intent(in) :: nc, nx, ny
      real(r8), intent(in) :: d_in(:,:,:)
      real(r8), intent(out) :: d_out(:,:,:)
      integer :: num_valid
      integer :: ip, ic, ix, iy, i2, i3
      call assert(size(d_in , 1) == nc * nx * ny, 'np /= nc * nx * ny')
      d_out = 0._r8
      do i3 = 1,size(d_in,3)
        do i2 = 1,size(d_in,2)
          do ic = 1,nc
            num_valid = 0
            do iy = 1,ny
              do ix = 1,nx
                ip = _IDX321(ic,ix,iy,nc,nx,ny)
                if (d_in(ip,i2,i3) /= R_UNDEF) then
                  d_out(ic,i2,i3) = d_out(ic,i2,i3) + d_in(ip,i2,i3)
                  num_valid = num_valid + 1
                end if
              end do
            end do
            ! If no valid values found, need to set to fill value
            if (num_valid > 0) then
              d_out(ic,i2,i3) = d_out(ic,i2,i3) / num_valid
            else
              d_out(ic,i2,i3) = R_UNDEF
            end if
          end do
        end do
      end do
    end subroutine
#endif

#ifdef USE_COSP
   subroutine cosp_unpack_outputs(nc, nx, ny, cospINpacked, cospOUTpacked, cospIN, cospOUT)
      use mod_cosp, only: cosp_optical_inputs, cosp_outputs
      integer, intent(in) :: nc, nx, ny  ! Packed dimensions
      type(cosp_optical_inputs), intent(in) :: cospINpacked
      type(cosp_outputs)       , intent(in) :: cospOUTpacked  ! inout so we can fix/weight before writing
      type(cosp_optical_inputs), intent(inout) :: cospIN
      type(cosp_outputs)       , intent(inout) :: cospOUT  ! inout so we can fix/weight before writing

      ! ISCCP OUTPUTS
      if (lisccp_sim) then
        call packed_average(nc, nx, ny, cospOUTpacked%isccp_fq           , cospOUT%isccp_fq           )
        call packed_average(nc, nx, ny, cospOUTpacked%isccp_totalcldarea , cospOUT%isccp_totalcldarea )
        call packed_average(nc, nx, ny, cospOUTpacked%isccp_meanalbedocld, cospOUT%isccp_meanalbedocld)
        call packed_average(nc, nx, ny, cospOUTpacked%isccp_meanptop     , cospOUT%isccp_meanptop     )
        call packed_average(nc, nx, ny, cospOUTpacked%isccp_meantaucld   , cospOUT%isccp_meantaucld   )
      end if

      ! CALIPSO SIMULATOR OUTPUTS
      if (llidar_sim) then
        call packed_average(nc, nx, ny, cospOUTpacked%calipso_cldlayer     , cospOUT%calipso_cldlayer     )
        call packed_average(nc, nx, ny, cospOUTpacked%calipso_cldlayerphase, cospOUT%calipso_cldlayerphase)
        call packed_average(nc, nx, ny, cospOUTpacked%calipso_lidarcld     , cospOUT%calipso_lidarcld     )
        call packed_average(nc, nx, ny, cospOUTpacked%calipso_beta_mol     , cospOUT%calipso_beta_mol     )
        call packed_average(nc, nx, ny, cospOUTpacked%calipso_cfad_sr      , cospOUT%calipso_cfad_sr      )
        call packed_average(nc, nx, ny, cospOUTpacked%parasolgrid_refl     , cospOUT%parasolgrid_refl     )
        call packed_average(nc, nx, ny, cospOUTpacked%calipso_lidarcldphase, cospOUT%calipso_lidarcldphase)
        call packed_average(nc, nx, ny, cospOUTpacked%calipso_lidarcldtmp  , cospOUT%calipso_lidarcldtmp  )
      end if
      
      ! RADAR SIMULATOR OUTPUTS
      if (lradar_sim) then
        call packed_average(nc, nx, ny, cospOUTpacked%cloudsat_cfad_ze     , cospOUT%cloudsat_cfad_ze     )
        call packed_average(nc, nx, ny, cospOUTpacked%radar_lidar_tcc      , cospOUT%radar_lidar_tcc      )
        call packed_average(nc, nx, ny, cospOUTpacked%lidar_only_freq_cloud, cospOUT%lidar_only_freq_cloud)
        call packed_average(nc, nx, ny, cospOUTpacked%cloudsat_precip_cover, cospOUT%cloudsat_precip_cover)
        call packed_average(nc, nx, ny, cospOUTpacked%cloudsat_pia         , cospOUT%cloudsat_pia         )
      end if
      
      ! MISR SIMULATOR OUTPUTS
      if (lmisr_sim) then
        call packed_average(nc, nx, ny, cospOUTpacked%misr_fq, cospOUT%misr_fq)
      end if
      
      ! MODIS SIMULATOR OUTPUTS
      if (lmodis_sim) then
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_fraction_total_mean  , cospOUT%modis_cloud_fraction_total_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_fraction_water_mean  , cospOUT%modis_cloud_fraction_water_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_fraction_ice_mean  , cospOUT%modis_cloud_fraction_ice_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_fraction_high_mean  , cospOUT%modis_cloud_fraction_high_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_fraction_mid_mean  , cospOUT%modis_cloud_fraction_mid_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_fraction_low_mean  , cospOUT%modis_cloud_fraction_low_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_total_mean  , cospOUT%modis_optical_thickness_total_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_water_mean  , cospOUT%modis_optical_thickness_water_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_ice_mean  , cospOUT%modis_optical_thickness_ice_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_total_logmean  , cospOUT%modis_optical_thickness_total_logmean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_water_logmean  , cospOUT%modis_optical_thickness_water_logmean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_ice_logmean  , cospOUT%modis_optical_thickness_ice_logmean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_particle_size_water_mean  , cospOUT%modis_cloud_particle_size_water_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_particle_size_ice_mean  , cospOUT%modis_cloud_particle_size_ice_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_cloud_top_pressure_total_mean  , cospOUT%modis_cloud_top_pressure_total_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_liquid_water_path_mean  , cospOUT%modis_liquid_water_path_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_ice_water_path_mean  , cospOUT%modis_ice_water_path_mean  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_vs_cloud_top_pressure, cospOUT%modis_optical_thickness_vs_cloud_top_pressure)
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_vs_reffice  , cospOUT%modis_optical_thickness_vs_reffice  )
        call packed_average(nc, nx, ny, cospOUTpacked%modis_optical_thickness_vs_reffliq  , cospOUT%modis_optical_thickness_vs_reffliq  )
      end if
      
      ! SUB-COLUMN OUTPUT
      ! NOTE: these are really not going to make much sense, but we'll go ahead
      ! and average them anyways for consistency. A better thing to do might be
      ! to repack to the subcol dimension, so ncol*nx*ny, nscol -> ncol, nx*ny*nscol
      if (lfrac_out) then
        call packed_average(nc, nx, ny, cospINpacked%frac_out, cospIN%frac_out)
        if (lisccp_sim) then
          call packed_average(nc, nx, ny, cospOUTpacked%isccp_boxtau, cospOUT%isccp_boxtau)
          call packed_average(nc, nx, ny, cospOUTpacked%isccp_boxptop, cospOUT%isccp_boxptop)
        end if
        if (llidar_sim) then
          call packed_average(nc, nx, ny, cospOUTpacked%calipso_beta_tot, cospOUT%calipso_beta_tot)
        end if
        if (lradar_sim) then
          call packed_average(nc, nx, ny, cospOUTpacked%cloudsat_ze_tot, cospOUT%cloudsat_ze_tot)
        end if
      end if
   end subroutine cosp_unpack_outputs
#endif


   ! Utility function to weight MODIS outputs by cloud fraction
   subroutine weight_output(arr, wgt, fillvalue)
      real(r8), intent(inout) :: arr(:)
      real(r8), intent(in)    :: wgt(:)
      real(r8), intent(in)    :: fillvalue
      where ((arr  .eq. fillvalue) .or. (wgt .eq. fillvalue))
         arr = fillvalue
      elsewhere
         arr = arr * wgt
      end where
   end subroutine weight_output

   ! Utility functions, used in cosp_write_outputs to replace fill values with
   ! zeros for 3d vars that may contain points below surface
   subroutine replace_values2d(arr, v1, v2)
      real(r8), intent(inout) :: arr(:,:)
      real(r8), intent(in) :: v1
      real(r8), intent(in) :: v2
      where (arr(:,:) == v1)
         arr(:,:) = v2
      end where
   end subroutine
   subroutine replace_values3d(arr, v1, v2)
      real(r8), intent(inout) :: arr(:,:,:)
      real(r8), intent(in) :: v1
      real(r8), intent(in) :: v2
      where (arr(:,:,:) == v1)
         arr(:,:,:) = v2
      end where
   end subroutine

#ifdef USE_COSP
   subroutine cosp_histfile_aux_out(state, cospstateIN, cospIN)
      use physics_types, only: physics_state
      use mod_cosp, only: cosp_column_inputs, cosp_optical_inputs
      use cam_history, only: outfld
      type(physics_state), intent(in) :: state
      type(cosp_column_inputs), intent(in) :: cospstateIN
      type(cosp_optical_inputs), intent(in) :: cospIN
      integer :: ncol, lchnk

      ! Number of columns and chunk ID
      ncol = state%ncol
      lchnk = state%lchnk
       
      ! 1D outputs
      call outfld('PS_COSP',        state%ps(1:ncol), ncol, lchnk)
      call outfld('TS_COSP',        cospstateIN%skt,  ncol, lchnk)
       
      ! 2D outputs
      call outfld('P_COSP',         cospstateIN%pfull,           ncol, lchnk)
      call outfld('PH_COSP',        cospstateIN%phalf,           ncol, lchnk)
      call outfld('ZLEV_COSP',      cospstateIN%hgt_matrix,      ncol, lchnk)
      call outfld('ZLEV_HALF_COSP', cospstateIN%hgt_matrix_half, ncol, lchnk)
      call outfld('T_COSP',         cospstateIN%at,              ncol, lchnk)
      call outfld('RH_COSP',        cospstateIN%qv,              ncol, lchnk)
      call outfld('Q_COSP',         cospstateIN%qv,              ncol, lchnk)

      ! 3D outputs
      call outfld('TAU_067',      cospIN%tau_067, ncol, lchnk)
      call outfld('EMISS_11',     cospIN%emiss_11, ncol, lchnk)
      call outfld('MODIS_asym',   cospIN%asym   , ncol, lchnk)
      call outfld('MODIS_ssa',    cospIN%ss_alb , ncol, lchnk)
      call outfld('MODIS_fracliq',cospIN%fracLiq, ncol, lchnk)
   end subroutine cosp_histfile_aux_out
#endif /* USE_COSP */

#ifdef USE_COSP
  ! ######################################################################################
  ! SUBROUTINE subsample_and_optics
  ! ######################################################################################
  subroutine subsample_and_optics(nPoints, nLevels, nColumns, nHydro,overlap,            &
                                  use_precipitation_fluxes, lidar_ice_type, sd, tca, cca,&
                                  fl_lsrainIN, fl_lssnowIN, fl_lsgrplIN, fl_ccrainIN,    &
                                  fl_ccsnowIN, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,   &
                                  reffIN, dtau_c, dtau_s, dem_c, dem_s, dtau_s_snow,     &
                                  dem_s_snow, sfcP, cospstateIN, cospIN)
    ! Dependencies
    use cosp_kinds,           only: wp
    use mod_rng,              only: rng_state, init_rng
    use mod_cosp_config,      only: R_UNDEF
    use mod_scops,            only: scops
    use mod_prec_scops,       only: prec_scops
    use mod_cosp_utils,       only: cosp_precip_mxratio
    use cosp_optics,          only: cosp_simulator_optics,lidar_optics,modis_optics,    &
                                    modis_optics_partition
    use mod_cosp_config,      only: Nlvgrid, vgrid_zl, vgrid_zu
    use mod_cosp_stats,       only: cosp_change_vertical_grid
    ! Inputs
    logical,intent(in) :: &
         use_precipitation_fluxes
    integer,intent(in) :: &
         nPoints,      & ! Number of gridpoints
         nLevels,      & ! Number of vertical levels
         nColumns,     & ! Number of subcolumns
         nHydro,       & ! Number pf hydrometeor types
         overlap,      & ! Overlap assumption (1/2/3)
         lidar_ice_type  ! Ice type assumption used by lidar optics
    real(wp),intent(in),dimension(nPoints,nLevels) :: &
         tca,          & ! Total cloud amount (0-1)
         cca,          & ! Convective cloud amount (0-1)
         mr_lsliq,     & ! Mixing ratio (kg/kg)
         mr_lsice,     & ! Mixing ratio (kg/kg)
         mr_ccliq,     & ! Mixing ratio (kg/kg)
         mr_ccice,     & ! Mixing ratio (kg/kg)
         dtau_c,       & ! 0.67-micron optical depth (convective)
         dtau_s,       & ! 0.67-micron optical depth (stratiform)
         dem_c,        & ! 11-micron emissivity (convective)
         dem_s,        & ! 11-micron emissivity (stratiform)
         fl_lsrainIN,  & ! Precipitation flux
         fl_lssnowIN,  & ! Precipitation flux
         fl_lsgrplIN,  & ! Precipitation flux
         fl_ccrainIN,  & ! Precipitation flux
         fl_ccsnowIN     ! Precipitation flux
    real(wp),intent(inout),dimension(nPoints,nLevels) :: &    
         dtau_s_snow,  & ! 0.67-micron optical depth (snow)
         dem_s_snow      ! 11-micron emissivity (snow)
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: &
         reffIN          !
    real(wp),intent(in),dimension(nPoints) :: &
         sfcP            ! Surface pressure 
    type(size_distribution), intent(inout) :: &
         sd
    
    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN
    
    ! Local variables
    integer :: i,j,k
    real(wp),dimension(nPoints,nLevels)      :: column_frac_out,column_prec_out,         &
                                                fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain, &
                                                fl_ccsnow
    real(wp),dimension(nPoints,nLevels,nHydro) :: ReffTemp                                                
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable         :: seed
    real(wp),dimension(:,:),allocatable      :: ls_p_rate,cv_p_rate,frac_ls,frac_cv,     &
                                                prec_ls,prec_cv,g_vol
    real(wp),dimension(:,:,:),  allocatable  :: frac_prec,&
                                                 MODIS_cloudWater,MODIS_cloudIce,        &
                                                 MODIS_watersize,MODIS_iceSize,          &
                                                 MODIS_snowSize,MODIS_cloudSnow,         &
                                                 MODIS_opticalThicknessLiq,              &
                                                 MODIS_opticalThicknessSnow,             &
                                                 MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable   :: mr_hydro,Reff,Np
             
    call t_startf("scops")
    if (Ncolumns .gt. 1) then
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! RNG used for subcolumn generation
       allocate(rngs(nPoints),seed(nPoints))
       seed = int(sfcP)
       if (Npoints .gt. 1) seed=(sfcP-int(sfcP))*1000000 
       call init_rng(rngs, seed)
   
       ! Call scops
       call scops(NPoints,Nlevels,Ncolumns,rngs,tca,cca,overlap,cospIN%frac_out,0)
       deallocate(seed,rngs)
       
       ! Sum up precipitation rates. If not using preciitation fluxes, mixing ratios are 
       ! stored in _rate variables.
       allocate(ls_p_rate(nPoints,nLevels),cv_p_rate(nPoints,Nlevels))
       if(use_precipitation_fluxes) then
          ls_p_rate(:,1:nLevels) = fl_lsrainIN + fl_lssnowIN + fl_lsgrplIN
          cv_p_rate(:,1:nLevels) = fl_ccrainIN + fl_ccsnowIN
       else
          ls_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow) + mixing_ratio (groupel)
          cv_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow)
       end if
       
       ! Call PREC_SCOPS
       allocate(frac_prec(nPoints,nColumns,nLevels))
       call prec_scops(nPoints,nLevels,nColumns,ls_p_rate,cv_p_rate,cospIN%frac_out,frac_prec)
       deallocate(ls_p_rate,cv_p_rate)
             
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute precipitation fraction in each gridbox
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
                if (frac_prec(j,i,k) .eq. 1)         prec_ls(j,k) = prec_ls(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 2)         prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)         prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)         prec_ls(j,k) = prec_ls(j,k)+1._wp
             enddo
             frac_ls(j,k)=frac_ls(j,k)/nColumns
             frac_cv(j,k)=frac_cv(j,k)/nColumns
             prec_ls(j,k)=prec_ls(j,k)/nColumns
             prec_cv(j,k)=prec_cv(j,k)/nColumns

             ! Adjust grid-box mean snow properties to local properties
             ! Convert longwave optical depth to longwave emissivity
             if (prec_ls(j,k) .ne. 0._r8 .and. dtau_s_snow(j,k) .gt. 0._r8) then
                dtau_s_snow(j,k) = dtau_s_snow(j,k)/prec_ls(j,k) 
             end if
             if (prec_ls(j,k) .ne. 0._r8 .and. dem_s_snow(j,k) .gt. 0._r8) then
                dem_s_snow(j,k) = dem_s_snow(j,k)/prec_ls(j,k)
                dem_s_snow(j,k) = 1._r8 - exp ( -1._r8*dem_s_snow(j,k))
             end if
          enddo
       enddo
             
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute mixing ratios, effective radii and precipitation fluxes for clouds
       ! and precipitation
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(mr_hydro(nPoints,nColumns,nLevels,nHydro),                               &
                Reff(nPoints,nColumns,nLevels,nHydro),                                   &
                Np(nPoints,nColumns,nLevels,nHydro))

       ! Initialize
       mr_hydro(:,:,:,:) = 0._wp
       Reff(:,:,:,:)     = 0._wp
       Np(:,:,:,:)       = 0._wp
       
       do k=1,nColumns
          ! Subcolumn clouds
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
       ! Convert the mixing ratio and precipitation fluxes from gridbox mean to
       ! the fraction-based values
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       do k=1,nLevels
          do j=1,nPoints
             ! Clouds
             if (frac_ls(j,k) .ne. 0._r8) then
                mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
             end if
             if (frac_cv(j,k) .ne. 0._r8) then
                mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
             end if
             
             ! Precipitation
             if (use_precipitation_fluxes) then
                if (prec_ls(j,k) .ne. 0._r8) then
                   fl_lsrain(j,k) = fl_lsrainIN(j,k)/prec_ls(j,k)
                   fl_lssnow(j,k) = fl_lssnowIN(j,k)/prec_ls(j,k)
                   fl_lsgrpl(j,k) = fl_lsgrplIN(j,k)/prec_ls(j,k)
                end if
                if (prec_cv(j,k) .ne. 0._r8) then
                   fl_ccrain(j,k) = fl_ccrainIN(j,k)/prec_cv(j,k)
                   fl_ccsnow(j,k) = fl_ccsnowIN(j,k)/prec_cv(j,k)
                end if
             else
                if (prec_ls(j,k) .ne. 0._r8) then
                   mr_hydro(j,:,k,I_LSRAIN) = mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                   mr_hydro(j,:,k,I_LSSNOW) = mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                   mr_hydro(j,:,k,I_LSGRPL) = mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
                end if
                if (prec_cv(j,k) .ne. 0._r8) then
                   mr_hydro(j,:,k,I_CVRAIN) = mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                   mr_hydro(j,:,k,I_CVSNOW) = mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
                end if
             end if
          enddo
       enddo
             
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert precipitation fluxes to mixing ratios
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (use_precipitation_fluxes) then
          ! LS rain
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSRAIN), n_bx(I_LSRAIN),         &
               alpha_x(I_LSRAIN), c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),       &
               a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN), gamma_2(I_LSRAIN),   &
               gamma_3(I_LSRAIN), gamma_4(I_LSRAIN), fl_lsrain,                          &
               mr_hydro(:,:,:,I_LSRAIN), Reff(:,:,:,I_LSRAIN))
          ! LS snow
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp,  n_ax(I_LSSNOW),  n_bx(I_LSSNOW),       &
               alpha_x(I_LSSNOW), c_x(I_LSSNOW),  d_x(I_LSSNOW),  g_x(I_LSSNOW),         &
               a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  gamma_2(I_LSSNOW),  &
               gamma_3(I_LSSNOW), gamma_4(I_LSSNOW), fl_lssnow,                          &
               mr_hydro(:,:,:,I_LSSNOW), Reff(:,:,:,I_LSSNOW))
          ! CV rain
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVRAIN),  n_bx(I_CVRAIN),        &
               alpha_x(I_CVRAIN), c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),       &
               a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN), gamma_2(I_CVRAIN),   &
               gamma_3(I_CVRAIN), gamma_4(I_CVRAIN), fl_ccrain,                          &
               mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN))
          ! CV snow
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVSNOW),  n_bx(I_CVSNOW),        &
               alpha_x(I_CVSNOW),  c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
               a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW), gamma_2(I_CVSNOW),   &
               gamma_3(I_CVSNOW), gamma_4(I_CVSNOW), fl_ccsnow,                          &
               mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW))
          ! LS groupel.
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSGRPL),  n_bx(I_LSGRPL),        &
               alpha_x(I_LSGRPL), c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),       &
               a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  gamma_2(I_LSGRPL),  &
               gamma_3(I_LSGRPL), gamma_4(I_LSGRPL), fl_lsgrpl,                          &
               mr_hydro(:,:,:,I_LSGRPL), Reff(:,:,:,I_LSGRPL))              
       end if

    else
       cospIN%frac_out(:,:,:) = 1  
       allocate(mr_hydro(nPoints, 1,nLevels,nHydro),Reff(nPoints,1,nLevels,nHydro),      &
                Np(nPoints,1,nLevels,nHydro))
       mr_hydro(:,1,:,I_LSCLIQ) = mr_lsliq
       mr_hydro(:,1,:,I_LSCICE) = mr_lsice
       mr_hydro(:,1,:,I_CVCLIQ) = mr_ccliq
       mr_hydro(:,1,:,I_CVCICE) = mr_ccice
       Reff(:,1,:,:)            = ReffIN
    end if
    call t_stopf("scops")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("cloudsat_optics")
    if (lradar_sim) then
      call cloudsat_optics( &
        Npoints, Ncolumns, Nlevels, Nlvgrid, Nhydro, &
        mr_hydro, Reff, Np, &
        sd, cospstateIN, cospIN &
      )
    end if
    call t_stopf("cloudsat_optics")
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CALIPSO Polarized optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("calipso_optics")
    if (Llidar_sim) then
       ReffTemp = ReffIN
       call lidar_optics(nPoints,nColumns,nLevels,5,lidar_ice_type,                      &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_LSCLIQ),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_LSCICE),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_CVCLIQ),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_CVCICE),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_LSSNOW),              &
                         ReffTemp(1:nPoints,1:nLevels,I_LSCLIQ),                         &
                         ReffTemp(1:nPoints,1:nLevels,I_LSCICE),                         &
                         ReffTemp(1:nPoints,1:nLevels,I_CVCLIQ),                         &
                         ReffTemp(1:nPoints,1:nLevels,I_CVCICE),                         & 
                         ReffTemp(1:nPoints,1:nLevels,I_LSSNOW),                         &
                         cospstateIN%pfull(1:nPoints,1:nLevels),                         &
                         cospstateIN%phalf(1:nPoints,1:nLevels+1),                       &
                         cospstateIN%at(1:nPoints,1:nLevels),                            &
                         cospIN%beta_mol_calipso(1:nPoints,1:nLevels),                   &
                         cospIN%betatot_calipso(1:nPoints,1:nColumns,1:nLevels),         &
                         cospIN%tau_mol_calipso(1:nPoints,1:nLevels),                    &
                         cospIN%tautot_calipso(1:nPoints,1:nColumns,1:nLevels),          &
                         cospIN%tautot_S_liq(1:nPoints,1:nColumns),                      &
                         cospIN%tautot_S_ice(1:nPoints,1:nColumns),                      &
                         cospIN%betatot_ice_calipso(1:nPoints,1:nColumns,1:nLevels),     &
                         cospIN%betatot_liq_calipso(1:nPoints,1:nColumns,1:nLevels),     &
                         cospIN%tautot_ice_calipso(1:nPoints,1:nColumns,1:nLevels),      &
                         cospIN%tautot_liq_calipso(1:nPoints,1:nColumns,1:nLevels)) 
    end if
    call t_stopf("calipso_optics")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Compute optical fields for passive simulators (i.e. only sunlit points)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 11 micron emissivity (needed by the ISCCP simulator)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("11micron_emissivity")
    if (Lisccp_sim) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dem_c,dem_s,  &
            cospIN%emiss_11)
       ! Add in contributions from radiative snow 
       do j=1,nColumns
          where(frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3)
             cospIN%emiss_11(:,j,:) = 1._wp - (1- cospIN%emiss_11(:,j,:))*(1-dem_s_snow)
          endwhere
       enddo
    end if
    call t_stopf("11micron_emissivity")
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 0.67 micron optical depth (needed by ISCCP, MISR and MODIS simulators)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("067tau")
    if (Lisccp_sim .or. Lmisr_sim .or. Lmodis_sim) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dtau_c,dtau_s,&
            cospIN%tau_067)
       
       ! Add in contributions from snow 
       do j=1,nColumns
          where((frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3) .and. &
               Reff(:,j,:,I_LSSNOW) .gt. 0._r8 .and. dtau_s_snow .gt. 0._r8)
             cospIN%tau_067(:,j,:)  = cospIN%tau_067(:,j,:)+dtau_s_snow
          endwhere
       enddo
    end if
    call t_stopf("067tau")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("modis_optics")
    if (lmodis_sim) then
       allocate(MODIS_cloudWater(nPoints,nColumns,nLevels),                              &
                MODIS_cloudIce(nPoints,nColumns,nLevels),                                &
                MODIS_cloudSnow(nPoints,nColumns,nLevels),                               &
                MODIS_waterSize(nPoints,nColumns,nLevels),                               &
                MODIS_iceSize(nPoints,nColumns,nLevels),                                 &
                MODIS_snowSize(nPoints,nColumns,nLevels),                                &
                MODIS_opticalThicknessLiq(nPoints,nColumns,nLevels),                     &
                MODIS_opticalThicknessIce(nPoints,nColumns,nLevels),                     &
                MODIS_opticalThicknessSnow(nPoints,nColumns,nLevels))

       ! Cloud water
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater)
       ! Cloud ice
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            mr_hydro(:,:,:,I_CVCICE),mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce)  
       ! Cloud water droplet size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            Reff(:,:,:,I_CVCLIQ),Reff(:,:,:,I_LSCLIQ),MODIS_waterSize)
       ! Cloud ice crystal size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            Reff(:,:,:,I_CVCICE),Reff(:,:,:,I_LSCICE),MODIS_iceSize)
       
       ! Cloud snow and size  
       MODIS_snowSize(:,:,:)  = Reff(:,:,:,I_LSSNOW)
       do j=1,nColumns
          where((frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3) .and. &
               Reff(:,j,:,I_LSSNOW) .gt. 0._r8 .and. dtau_s_snow .gt. 0._r8)
             MODIS_cloudSnow(:,j,:) = mr_hydro(:,j,:,I_LSSNOW)
             MODIS_snowSize(:,j,:)  = Reff(:,j,:,I_LSSNOW)
          elsewhere
             MODIS_snowSize(:,j,:)  = 0._wp
             MODIS_cloudSnow(:,j,:) = 0._wp
          endwhere
       enddo
       
       ! Partition optical thickness into liquid and ice parts
       call modis_optics_partition(nPoints, nLevels, nColumns, MODIS_cloudWater,     &
            MODIS_cloudIce, MODIS_cloudSnow, MODIS_waterSize, MODIS_iceSize,         &
            MODIS_snowSize, cospIN%tau_067, MODIS_opticalThicknessLiq,               &
            MODIS_opticalThicknessIce, MODIS_opticalThicknessSnow)                            
       
       ! Compute assymetry parameter and single scattering albedo 
       call modis_optics(nPoints, nLevels, nColumns, MODIS_opticalThicknessLiq,      &
            MODIS_waterSize*1.0e6_wp, MODIS_opticalThicknessIce,                     &
            MODIS_iceSize*1.0e6_wp, MODIS_opticalThicknessSnow,                      &
            MODIS_snowSize*1.0e6_wp, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)

    end if ! MODIS simulator optics
    call t_stopf("modis_optics")

  end subroutine subsample_and_optics


  subroutine cloudsat_optics(npoints, ncolumns, nlevels, nlvgrid, nhydro, mr_hydro, Reff, Np, sd, cospstateIN, cospIN)
    use cosp_kinds, only: wp
    use mod_cosp_config, only: R_UNDEF, vgrid_zl, vgrid_zu
    use mod_cosp, only: cosp_column_inputs, cosp_optical_inputs
    use mod_quickbeam_optics, only: size_distribution, gases, quickbeam_optics
    integer, intent(in) :: npoints, ncolumns, nlevels, nlvgrid, nhydro
    real(wp), dimension(npoints,ncolumns,nlevels,nhydro), intent(in) :: &
      mr_hydro, &  ! Mixing ratios
      Reff         ! Effective radii
    real(wp), dimension(npoints,ncolumns,nlevels,nhydro), intent(inout) :: &
      Np           ! Number concentrations
    type(size_distribution)  , intent(inout) :: sd
    type(cosp_optical_inputs), intent(inout) :: cospIN
    type(cosp_column_inputs) , intent(inout) :: cospstateIN
    real(wp), dimension(npoints,ncolumns,nlevels) :: fracPrecipIce
    real(wp), dimension(npoints,ncolumns,nlvgrid) :: fracPrecipIce_statGrid
    real(wp), dimension(npoints,nlevels)          :: g_vol
    integer :: i, j, k
                                               
       ! Compute gaseous absorption (assume identical for each subcolun)
       g_vol(:,:)=0._wp
       do i = 1, nPoints
          do j = 1, nLevels
             if (cospIN%rcfg_cloudsat%use_gas_abs == 1 .or. &
                (cospIN%rcfg_cloudsat%use_gas_abs == 2 .and. j == 1)) then
                g_vol(i,j) = gases(cospstateIN%pfull(i,j), cospstateIN%at(i,j),    &
                                   cospstateIN%qv(i,j), cospIN%rcfg_cloudsat%freq)
             end if
             cospIN%g_vol_cloudsat(i,:,j) = g_vol(i,j)
          end do
       end do

       ! Loop over all subcolumns
       fracPrecipIce(:,:,:) = 0._wp
       do k=1,nColumns
          call quickbeam_optics(sd, cospIN%rcfg_cloudsat, nPoints, nLevels, R_UNDEF, &
               mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,      &
               Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,                &
               cospstateIN%qv, cospIN%z_vol_cloudsat(1:nPoints,k,:),                 &
               cospIN%kr_vol_cloudsat(1:nPoints,k,:))

          ! At each model level, what fraction of the precipitation is frozen?
          where(mr_hydro(:,k,:,I_LSRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_LSSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_CVRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_CVSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_LSGRPL) .gt. 0)
             fracPrecipIce(:,k,:) = &
               (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + mr_hydro(:,k,:,I_LSGRPL)) / &
               (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + mr_hydro(:,k,:,I_LSGRPL) + &
                mr_hydro(:,k,:,I_LSRAIN) + mr_hydro(:,k,:,I_CVRAIN))
          elsewhere
             fracPrecipIce(:,k,:) = 0._wp
          endwhere
       enddo

       ! Regrid frozen fraction to Cloudsat/Calipso statistical grid
       fracPrecipIce_statGrid(:,:,:) = 0._wp
       call cosp_change_vertical_grid(npoints, ncolumns, nlevels, cospstateIN%hgt_matrix(:,nlevels:1:-1), &
            cospstateIN%hgt_matrix_half(:,nlevels:1:-1), fracPrecipIce(:,:,nlevels:1:-1), nlvgrid,  &
            vgrid_zl(nlvgrid:1:-1),  vgrid_zu(nlvgrid:1:-1), fracPrecipIce_statGrid(:,:,nlvgrid:1:-1))

       ! For near-surface diagnostics, we only need the frozen fraction at one layer.
       cospIN%fracPrecipIce(:,:) = fracPrecipIce_statGrid(:,:,cloudsat_preclvl)
       
  end subroutine cloudsat_optics
  
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
    
    allocate(y%tau_067(            npoints, ncolumns, nlevels),&
             y%emiss_11(           npoints, ncolumns, nlevels),&
             y%frac_out(           npoints, ncolumns, nlevels),&
             y%betatot_calipso(    npoints, ncolumns, nlevels),&
             y%betatot_ice_calipso(npoints, ncolumns, nlevels),&
             y%fracLiq(            npoints, ncolumns, nlevels),&
             y%betatot_liq_calipso(npoints, ncolumns, nlevels),&
             y%tautot_calipso(     npoints, ncolumns, nlevels),&
             y%tautot_ice_calipso( npoints, ncolumns, nlevels),&
             y%tautot_liq_calipso( npoints, ncolumns, nlevels),&
             y%z_vol_cloudsat(     npoints, ncolumns, nlevels),&
             y%kr_vol_cloudsat(    npoints, ncolumns, nlevels),&
             y%g_vol_cloudsat(     npoints, ncolumns, nlevels),&
             y%asym(               npoints, ncolumns, nlevels),&
             y%ss_alb(             npoints, ncolumns, nlevels),&
             y%beta_mol_calipso(   npoints,           nlevels),&
             y%tau_mol_calipso(    npoints,           nlevels),&
             y%tautot_S_ice(       npoints, ncolumns         ),&
             y%tautot_S_liq(       npoints, ncolumns)         ,&
             y%fracPrecipIce(npoints,   ncolumns))
    allocate(y%rcfg_cloudsat)

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
    
    y%npoints  = npoints
    y%ncolumns = ncolumns
    y%nlevels  = nlevels
    allocate(y%sunlit(npoints),y%skt(npoints),y%land(npoints),y%at(npoints,nlevels),     &
             y%pfull(npoints,nlevels),y%phalf(npoints,nlevels+1),y%qv(npoints,nlevels),  &
             y%o3(npoints,nlevels),y%hgt_matrix(npoints,nlevels),y%u_sfc(npoints),       &
             y%v_sfc(npoints),y%lat(npoints),y%lon(nPoints),y%emis_sfc(nchan),           &
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(nPoints),&
             y%fl_snow(nPoints,nLevels),y%fl_rain(nPoints,nLevels),y%seaice(npoints),    &
             y%tca(nPoints,nLevels),y%hgt_matrix_half(npoints,nlevels))

  end subroutine construct_cospstateIN
  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  ! ######################################################################################  
  subroutine construct_cosp_outputs(Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
    ! Inputs
    integer,intent(in) :: &
         Npoints,         & ! Number of sampled points
         Ncolumns,        & ! Number of subgrid columns
         Nlevels,         & ! Number of model levels
         Nlvgrid,         & ! Number of levels in L3 stats computation
         Nchan              ! Number of RTTOV channels  
    
    ! Outputs
    type(cosp_outputs),intent(out) :: &
         x           ! COSP output structure  
  
     ! ISCCP simulator outputs
    if (lisccp_sim) then
       allocate(x%isccp_boxtau(Npoints,Ncolumns)) 
       allocate(x%isccp_boxptop(Npoints,Ncolumns))
       allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
       allocate(x%isccp_totalcldarea(Npoints))
       allocate(x%isccp_meanptop(Npoints))
       allocate(x%isccp_meantaucld(Npoints))
       allocate(x%isccp_meantb(Npoints))
       allocate(x%isccp_meantbclr(Npoints))
       allocate(x%isccp_meanalbedocld(Npoints))
    end if

    ! MISR simulator
    if (lmisr_sim) then 
       allocate(x%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins))
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.
       allocate(x%misr_dist_model_layertops(Npoints,numMISRHgtBins))
       allocate(x%misr_meanztop(Npoints))
       allocate(x%misr_cldarea(Npoints))    
    end if
    
    ! MODIS simulator
    if (lmodis_sim) then
       allocate(x%modis_Cloud_Fraction_Total_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Water_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Ice_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_High_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Mid_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Low_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Total_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Water_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Ice_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Total_LogMean(Npoints))
       allocate(x%modis_Optical_Thickness_Water_LogMean(Npoints))
       allocate(x%modis_Optical_Thickness_Ice_LogMean(Npoints))
       allocate(x%modis_Cloud_Particle_Size_Water_Mean(Npoints))
       allocate(x%modis_Cloud_Particle_Size_Ice_Mean(Npoints))
       allocate(x%modis_Cloud_Top_Pressure_Total_Mean(Npoints))
       allocate(x%modis_Liquid_Water_Path_Mean(Npoints))
       allocate(x%modis_Ice_Water_Path_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins,numMODISPresBins))
       allocate(x%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins,numMODISReffLiqBins))   
       allocate(x%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins,numMODISReffIceBins))
    end if
    
    ! CALIPSO simulator
    if (llidar_sim) then
       allocate(x%calipso_beta_mol(Npoints,Nlevels))
       allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%calipso_srbval(SR_BINS+1))
       allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid))
       allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))  
       allocate(x%calipso_lidarcld(Npoints,Nlvgrid))
       allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))        
       allocate(x%calipso_lidarcldphase(Npoints,Nlvgrid,6))
       allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
       allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))     
       ! These 2 outputs are part of the calipso output type, but are not controlled by an 
       ! logical switch in the output namelist, so if all other fields are on, then allocate
       allocate(x%calipso_tau_tot(Npoints,Ncolumns,Nlevels))       
       allocate(x%calipso_temp_tot(Npoints,Nlevels))               
       ! Calipso opaque cloud diagnostics
!       allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
!       allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))  
!       allocate(x%calipso_cldtypemeanz(Npoints,2)) 
!       allocate(x%calipso_cldtypemeanzse(Npoints,3)) 
!       allocate(x%calipso_cldthinemis(Npoints))
!       allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))
    end if 
      
    ! PARASOL
    if (lparasol_sim) then
       allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
       allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    end if

    ! Cloudsat simulator
    if (lradar_sim) then
       allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%cloudsat_cfad_ze(Npoints,CLOUDSAT_DBZE_BINS,Nlvgrid))
       allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
       allocate(x%radar_lidar_tcc(Npoints))
       allocate(x%cloudsat_precip_cover(Npoints,nCloudsatPrecipClass))
       allocate(x%cloudsat_pia(Npoints))
    end if

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
    if (allocated(y%fracPrecipIce))       deallocate(y%fracPrecipIce)
    if (allocated(y%rcfg_cloudsat))       deallocate(y%rcfg_cloudsat)

  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%surfelev))        deallocate(y%surfelev)
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
     end if
     if (associated(y%calipso_temp_tot))          then
        deallocate(y%calipso_temp_tot)
        nullify(y%calipso_temp_tot)     
     end if
     if (associated(y%calipso_betaperp_tot))      then
        deallocate(y%calipso_betaperp_tot)
        nullify(y%calipso_betaperp_tot)     
     end if
     if (associated(y%calipso_beta_tot))          then
        deallocate(y%calipso_beta_tot)    
        nullify(y%calipso_beta_tot)     
     end if
     if (associated(y%calipso_tau_tot))           then
        deallocate(y%calipso_tau_tot) 
        nullify(y%calipso_tau_tot)     
     end if
     if (associated(y%calipso_lidarcldphase))     then
        deallocate(y%calipso_lidarcldphase)
        nullify(y%calipso_lidarcldphase)     
     end if
     if (associated(y%calipso_cldlayerphase))     then
        deallocate(y%calipso_cldlayerphase)
        nullify(y%calipso_cldlayerphase)     
     end if
     if (associated(y%calipso_lidarcldtmp))       then
        deallocate(y%calipso_lidarcldtmp)
        nullify(y%calipso_lidarcldtmp)     
     end if
     if (associated(y%calipso_cldlayer))          then
        deallocate(y%calipso_cldlayer)
        nullify(y%calipso_cldlayer)     
     end if
     if (associated(y%calipso_lidarcld))         then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)     
     end if
     if (associated(y%calipso_srbval))            then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)     
     end if
     if (associated(y%calipso_cfad_sr))          then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)     
     end if
     if (associated(y%parasolPix_refl))           then
        deallocate(y%parasolPix_refl)
        nullify(y%parasolPix_refl)     
     end if
     if (associated(y%parasolGrid_refl))          then
        deallocate(y%parasolGrid_refl) 
        nullify(y%parasolGrid_refl)     
     end if
     if (associated(y%cloudsat_Ze_tot))           then
        deallocate(y%cloudsat_Ze_tot) 
        nullify(y%cloudsat_Ze_tot)  
     end if
     if (associated(y%cloudsat_precip_cover)) then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     end if
     if (associated(y%cloudsat_pia)) then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     end if
     if (associated(y%cloudsat_cfad_ze))          then
        deallocate(y%cloudsat_cfad_ze)
        nullify(y%cloudsat_cfad_ze)     
     end if
     if (associated(y%radar_lidar_tcc))           then
        deallocate(y%radar_lidar_tcc) 
        nullify(y%radar_lidar_tcc)  
     end if
     if (associated(y%lidar_only_freq_cloud))     then
        deallocate(y%lidar_only_freq_cloud)
        nullify(y%lidar_only_freq_cloud)     
     end if
     if (associated(y%isccp_totalcldarea))        then
        deallocate(y%isccp_totalcldarea) 
        nullify(y%isccp_totalcldarea)  
     end if
     if (associated(y%isccp_meantb))              then
        deallocate(y%isccp_meantb) 
        nullify(y%isccp_meantb)     
     end if
     if (associated(y%isccp_meantbclr))           then
        deallocate(y%isccp_meantbclr)
        nullify(y%isccp_meantbclr)  
     end if
     if (associated(y%isccp_meanptop))            then
        deallocate(y%isccp_meanptop)
        nullify(y%isccp_meanptop)     
     end if
     if (associated(y%isccp_meantaucld))          then
        deallocate(y%isccp_meantaucld) 
        nullify(y%isccp_meantaucld)       
     end if
     if (associated(y%isccp_meanalbedocld))       then
        deallocate(y%isccp_meanalbedocld)
        nullify(y%isccp_meanalbedocld)     
     end if
     if (associated(y%isccp_boxtau))              then
        deallocate(y%isccp_boxtau)
        nullify(y%isccp_boxtau)       
     end if
     if (associated(y%isccp_boxptop))             then
        deallocate(y%isccp_boxptop)
        nullify(y%isccp_boxptop)     
     end if
     if (associated(y%isccp_fq))                  then
        deallocate(y%isccp_fq)
        nullify(y%isccp_fq)       
     end if
     if (associated(y%misr_fq))                   then
        deallocate(y%misr_fq) 
        nullify(y%misr_fq)     
     end if
     if (associated(y%misr_dist_model_layertops)) then
        deallocate(y%misr_dist_model_layertops)
        nullify(y%misr_dist_model_layertops)       
     end if
     if (associated(y%misr_meanztop))             then
        deallocate(y%misr_meanztop)
        nullify(y%misr_meanztop)     
     end if
     if (associated(y%misr_cldarea))              then
        deallocate(y%misr_cldarea)
        nullify(y%misr_cldarea)      
     end if
     if (associated(y%rttov_tbs))                 then
        deallocate(y%rttov_tbs)
        nullify(y%rttov_tbs)     
     end if
     if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Total_Mean)       
        nullify(y%modis_Cloud_Fraction_Total_Mean)       
     end if
     if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)     
        nullify(y%modis_Cloud_Fraction_Ice_Mean)     
     end if
     if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Water_Mean)           
        nullify(y%modis_Cloud_Fraction_Water_Mean)           
     end if
     if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
        deallocate(y%modis_Cloud_Fraction_High_Mean)     
        nullify(y%modis_Cloud_Fraction_High_Mean)     
     end if
     if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)       
        nullify(y%modis_Cloud_Fraction_Mid_Mean)       
     end if
     if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Low_Mean)     
        nullify(y%modis_Cloud_Fraction_Low_Mean)     
     end if
     if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Total_Mean)  
        nullify(y%modis_Optical_Thickness_Total_Mean)  
     end if
     if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Water_Mean)     
        nullify(y%modis_Optical_Thickness_Water_Mean)     
     end if
     if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
        deallocate(y%modis_Optical_Thickness_Ice_Mean)       
        nullify(y%modis_Optical_Thickness_Ice_Mean)       
     end if
     if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Total_LogMean)    
        nullify(y%modis_Optical_Thickness_Total_LogMean)    
     end if
     if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Water_LogMean)     
        nullify(y%modis_Optical_Thickness_Water_LogMean)     
     end if
     if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)     
        nullify(y%modis_Optical_Thickness_Ice_LogMean)     
     end if
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)       
        nullify(y%modis_Cloud_Particle_Size_Water_Mean)       
     end if
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)     
        nullify(y%modis_Cloud_Particle_Size_Ice_Mean)     
     end if
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)           
        nullify(y%modis_Cloud_Top_Pressure_Total_Mean)           
     end if
     if (associated(y%modis_Liquid_Water_Path_Mean))                         then
        deallocate(y%modis_Liquid_Water_Path_Mean)     
        nullify(y%modis_Liquid_Water_Path_Mean)     
     end if
     if (associated(y%modis_Ice_Water_Path_Mean))                            then
        deallocate(y%modis_Ice_Water_Path_Mean)       
        nullify(y%modis_Ice_Water_Path_Mean)       
     end if
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
        nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
     end if
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
        nullify(y%modis_Optical_thickness_vs_ReffLIQ)
     end if
     if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        nullify(y%modis_Optical_thickness_vs_ReffICE)
     end if
     if (associated(y%calipso_cldtype)) then
        deallocate(y%calipso_cldtype)
        nullify(y%calipso_cldtype)
     end if
     if (associated(y%calipso_cldtypetemp)) then
        deallocate(y%calipso_cldtypetemp) 
        nullify(y%calipso_cldtypetemp) 
     end if
     if (associated(y%calipso_cldtypemeanz)) then
        deallocate(y%calipso_cldtypemeanz) 
        nullify(y%calipso_cldtypemeanz) 
     end if
     if (associated(y%calipso_cldtypemeanzse)) then
        deallocate(y%calipso_cldtypemeanzse) 
        nullify(y%calipso_cldtypemeanzse) 
     end if
     if (associated(y%calipso_cldthinemis)) then
        deallocate(y%calipso_cldthinemis)
        nullify(y%calipso_cldthinemis)
     end if
     if (associated(y%calipso_lidarcldtype)) then
        deallocate(y%calipso_lidarcldtype)
        nullify(y%calipso_lidarcldtype)
     end if
        
   end subroutine destroy_cosp_outputs
#endif

   ! Assert just checks that condition is true, and aborts execution of the
   ! program if it is not.
   subroutine assert(condition, message)
      use cam_abortutils,       only: endrun
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         call endrun('Assertion failed: ' // message)
      end if
   end subroutine
   !-------------------------------------------------------------------------------

!#######################################################################
end module cospsimulator_intr
