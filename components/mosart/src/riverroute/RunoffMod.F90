module RunoffMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RunoffMod
!
! !DESCRIPTION:
! Module containing utilities for history file and coupler runoff data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mct_mod
  use RtmVar         , only : iulog, spval, heatflag, data_bgc_fluxes_to_ocean_flag
  use rof_cpl_indices, only : nt_rtm

! !PUBLIC TYPES:
  implicit none
  private

  type(mct_gsmap),public :: gsmap_r       ! gsmap for mosart decomposition

  ! sum of upstream data, dst gets sum of values from upstream
  type(mct_sMatP),public :: sMatP_upstrm  ! sparse matrix plus for downstream communication
  type(mct_avect),public :: avsrc_upstrm  ! src avect for SM mult downstream communication
  type(mct_avect),public :: avdst_upstrm  ! dst avect for SM mult downstream communication

  ! copy of downstream data, dst gets single value from downstream
  type(mct_sMatP),public :: sMatP_dnstrm  ! sparse matrix plus for upstream communication
  type(mct_avect),public :: avsrc_dnstrm  ! src avect for SM mult upstream communication
  type(mct_avect),public :: avdst_dnstrm  ! dst avect for SM mult upstream communication

  ! communication to ocean outlet, dst get direct values from upstream
  type(mct_sMatP),public :: sMatP_direct  ! sparse matrix plus for direct to outlet flow
  type(mct_avect),public :: avsrc_direct  ! src avect for SM mult direct to outlet flow
  type(mct_avect),public :: avdst_direct  ! dst avect for SM mult direct to outlet flow

  public :: runoff_flow
  type runoff_flow
     !    - local initialization
     real(r8), pointer :: lonc(:)          ! lon of cell
     real(r8), pointer :: latc(:)          ! lat of cell
     real(r8), pointer :: area(:)          ! area of cell
     integer , pointer :: mask(:)          ! general mask of cell 1=land, 2=ocean, 3=outlet
     integer , pointer :: gindex(:)        ! global index consistent with map file
     integer , pointer :: dsig(:)          ! downstream index, global index
     integer , pointer :: outletg(:)       ! outlet index, global index
     real(r8), pointer :: rmask(:)         ! general mask of cell 1=land, 2=ocean, 3=outlet
     real(r8), pointer :: rgindex(:)       ! global index consistent with map file
     real(r8), pointer :: rdsig(:)         ! downstream index, global index
     real(r8), pointer :: routletg(:)      ! outlet index, global index

     !    - global 
     real(r8), pointer :: rlon(:)          ! rtm longitude list, 1d
     real(r8), pointer :: rlat(:)          ! rtm latitude list, 1d
     real(r8)          :: totarea          ! global area
     integer           :: numr             ! rtm gdc global number of cells

     !    - local
     integer           :: begr,endr        ! local start/stop indices
     integer           :: lnumr            ! local number of cells
     integer , pointer :: iDown(:)         ! downstream index, local index
     integer , pointer :: nUp(:)           ! number of upstream units, maximum 8
     integer , pointer :: nUp_dstrm(:)      ! number of units flowing into the downstream units, maximum 8
     integer , pointer :: iUp(:,:)         ! indices of upstream units, local

     !    - local
     real(r8), pointer :: runofflnd(:,:)   ! runoff masked for land (m3 H2O/s)
     real(r8), pointer :: runoffocn(:,:)   ! runoff masked for ocn  (m3 H2O/s)
     real(r8), pointer :: runofftot(:,:)   ! total runoff masked for ocn  (m3 H2O/s)
     real(r8), pointer :: dvolrdt(:,:)     ! RTM change in storage (mm/s)
     real(r8), pointer :: dvolrdtlnd(:,:)  ! dvolrdt masked for land (mm/s)
     real(r8), pointer :: dvolrdtocn(:,:)  ! dvolrdt masked for ocn  (mm/s)
     real(r8), pointer :: volr(:,:)        ! RTM storage (m3)
     real(r8), pointer :: fthresh(:)       ! RTM water flood threshold
     real(r8), pointer :: inundffunit(:)   ! Inundation  water volume (m3) 
     real(r8), pointer :: inundwf(:)       ! Inundation floodplain water volume (m3)
     real(r8), pointer :: inundhf(:)       ! Inundation floodplain water depth (m)
     real(r8), pointer :: inundff(:)       ! Inundation floodplain water area fraction (no unit)
     real(r8), pointer :: concDIN(:)       ! Concentration of river outflow DIN (kg-N/kg-water)
     real(r8), pointer :: concDIP(:)       ! Concentration of river outflow DIP (kg-P/kg-water)
     real(r8), pointer :: concDON(:)       ! Concentration of river outflow DON (kg-N/kg-water)
     real(r8), pointer :: concDOP(:)       ! Concentration of river outflow DOP (kg-P/kg-water)
     real(r8), pointer :: concDOC(:)       ! Concentration of river outflow DOC (kg-C/kg-water)
     real(r8), pointer :: concPP(:)        ! Concentration of river outflow PP (kg-P/kg-water)
     real(r8), pointer :: concDSi(:)       ! Concentration of river outflow DSi (kg-Si/kg-water)
     real(r8), pointer :: concPOC(:)       ! Concentration of river outflow POC (kg-C/kg-water)
     real(r8), pointer :: concPN(:)        ! Concentration of river outflow PN (kg-N/kg-water)
     real(r8), pointer :: concDIC(:)       ! Concentration of river outflow DIC (kg-C/kg-water)
     real(r8), pointer :: concFe(:)        ! Concentration of river outflow Fe (kg-Fe/kg-water)

     !    - restarts
     real(r8), pointer :: wh(:,:)          ! MOSART hillslope surface water storage (m)
     real(r8), pointer :: wt(:,:)          ! MOSART sub-network water storage (m3)
     real(r8), pointer :: wr(:,:)          ! MOSART main channel water storage (m3)
     real(r8), pointer :: mr(:,:)          ! MOSART channel area
     real(r8), pointer :: yr(:,:)          ! MOSART channel water depth
     real(r8), pointer :: pr(:,:)          ! MOSART channel wetted p
     real(r8), pointer :: rr(:,:)          ! MOSART channel hydraulic r
     real(r8), pointer :: erout(:,:)       ! MOSART flow out of the main channel, instantaneous (m3/s) (negative is out)
     real(r8), pointer :: Tqsur(:)         ! MOSART hillslope surface runoff water temperature (K)
     real(r8), pointer :: Tqsub(:)         ! MOSART hillslope subsurface runoff water temperature (K)
     real(r8), pointer :: Tt(:)            ! MOSART sub-network water temperature (K)
     real(r8), pointer :: Tr(:)            ! MOSART main channel water temperature (K)
     real(r8), pointer :: Ha_rout(:)       ! MOSART heat flux out of the main channel, instantaneous (Watt)
     real(r8), pointer :: wt_al(:,:)       ! MOSART sub-network channel active layer storage (kg for mud and sand sediment)
     real(r8), pointer :: wr_al(:,:)       ! MOSART main channel active layer storage (kg for mud and sand sediment)
     real(r8), pointer :: wres(:,:)        ! MOSART reservoir storage (m3 for water, kg for mud and sand sediment)

     ! inputs
     real(r8), pointer :: qsur(:,:)        ! coupler surface forcing [m3/s]
     real(r8), pointer :: qsub(:,:)        ! coupler subsurface forcing [m3/s]
     real(r8), pointer :: qgwl(:,:)        ! coupler glacier/wetland/lake forcing [m3/s]
     real(r8), pointer :: qdto(:,:)        ! coupler diret-to-ocean forcing [m3/s]
     real(r8), pointer :: qdem(:,:)        ! coupler total demand diagnostic [m3/s]

     !    - outputs
     real(r8), pointer :: flood(:)         ! coupler return flood water sent back to clm [m3/s]
     real(r8), pointer :: runoff(:,:)      ! coupler return mosart basin derived flow [m3/s]
     real(r8), pointer :: direct(:,:)      ! coupler return direct flow [m3/s]
     real(r8), pointer :: inundinf(:)      ! coupler return drainage from floodplain inundation  [mm/s]

     !    - history (currently needed)
     real(r8), pointer :: runofflnd_nt1(:)
     real(r8), pointer :: runofflnd_nt2(:)
     real(r8), pointer :: runofflnd_nt3(:)
     real(r8), pointer :: runofflnd_nt4(:)
     real(r8), pointer :: runoffocn_nt1(:)
     real(r8), pointer :: runoffocn_nt2(:)
     real(r8), pointer :: runoffocn_nt3(:)
     real(r8), pointer :: runoffocn_nt4(:)
     real(r8), pointer :: runofftot_nt1(:)
     real(r8), pointer :: runofftot_nt2(:)
     real(r8), pointer :: runofftot_nt3(:)
     real(r8), pointer :: runofftot_nt4(:)
     real(r8), pointer :: runoffdir_nt1(:)
     real(r8), pointer :: runoffdir_nt2(:)
     real(r8), pointer :: runoffdir_nt3(:)
     real(r8), pointer :: runoffdir_nt4(:)
     real(r8), pointer :: dvolrdtlnd_nt1(:)
     real(r8), pointer :: dvolrdtlnd_nt2(:)
     real(r8), pointer :: dvolrdtlnd_nt3(:)
     real(r8), pointer :: dvolrdtlnd_nt4(:)
     real(r8), pointer :: dvolrdtocn_nt1(:)
     real(r8), pointer :: dvolrdtocn_nt2(:)
     real(r8), pointer :: wr_nt1(:)
     real(r8), pointer :: dvolrdtocn_nt3(:)
     real(r8), pointer :: dvolrdtocn_nt4(:)
     real(r8), pointer :: volr_nt1(:)
     real(r8), pointer :: volr_nt2(:)
     real(r8), pointer :: volr_nt3(:)
     real(r8), pointer :: volr_nt4(:)
     real(r8), pointer :: qsur_nt1(:)
     real(r8), pointer :: qsur_nt2(:)
     real(r8), pointer :: qsur_nt3(:)
     real(r8), pointer :: qsub_nt1(:)
     real(r8), pointer :: qsub_nt2(:)
     real(r8), pointer :: qgwl_nt1(:)
     real(r8), pointer :: qgwl_nt2(:)
     real(r8), pointer :: qdto_nt1(:)
     real(r8), pointer :: qdto_nt2(:)
     real(r8), pointer :: qdem_nt1(:)
     real(r8), pointer :: qdem_nt2(:)

     real(r8), pointer :: templand_Tqsur(:)
     real(r8), pointer :: templand_Tqsub(:)
     real(r8), pointer :: templand_Ttrib(:)
     real(r8), pointer :: templand_Tchanr(:)     
     real(r8), pointer :: templand_Tqsur_nt1(:)
     real(r8), pointer :: templand_Tqsub_nt1(:)
     real(r8), pointer :: templand_Ttrib_nt1(:)
     real(r8), pointer :: templand_Tchanr_nt1(:)
     real(r8), pointer :: templand_Tqsur_nt2(:)
     real(r8), pointer :: templand_Tqsub_nt2(:)
     real(r8), pointer :: templand_Ttrib_nt2(:)
     real(r8), pointer :: templand_Tchanr_nt2(:)

     real(r8), pointer :: ssh(:)
     real(r8), pointer :: yr_nt1(:)
     
  end type runoff_flow

  
  !== Hongyi
  ! constrol information 
  public :: Tcontrol
  type Tcontrol
     integer  :: NSTART           ! the # of the time step to start the routing. Previous NSTART - 1 steps will be passed over.
     integer  :: NSTEPS           ! number of time steps specified in the modeling
     integer  :: NWARMUP          ! time steps for model warming up
     real(r8) :: DATAH            ! time step of runoff generation in second provided by the user
     integer  :: Num_dt           ! number of sub-steps within the current step interval, 
                                  ! i.e., if the time step of the incoming runoff data is 3-hr, and num_dt is set to 10, 
                                  ! then deltaT = 3*3600/10 = 1080 seconds
     real(r8) :: DeltaT           ! Time step in seconds 
     real(r8) :: coupling_period  ! couping period of rof in seconds
     integer  :: DLevelH2R        ! The base number of channel routing sub-time-steps within one hillslope routing step. 
                                  ! Usually channel routing requires small time steps than hillslope routing.
     integer  :: DLevelR          ! The number of channel routing sub-time-steps at a higher level within one channel routing step at a lower level. 
     integer  :: Restart          ! flag, Restart=1 means starting from the state of last run, =0 means starting from model-inset initial state.
     integer  :: RoutingMethod    ! Flag for routing methods. 1 --> Kinematic wave routing method; 2 --> Diffusion wave method.
     integer  :: RoutingFlag      ! Flag for whether including hillslope and sub-network routing. 1--> include routing through hillslope, sub-network and main channel; 0--> main channel routing only.
 
     character(len=100) :: baseName    ! name of the case study, e.g., columbia
     character(len=200) :: ctlFile     ! the name of the control file
     character(len=100) :: ctlPath     ! the path of the control file
     character(len=100) :: runoffPath  ! the path of the runoff data
     character(len=100) :: outPath     ! the path of the output file(s)
     integer :: numStation             ! number of basins to be simulated
     character(len=200) :: staListFile ! name of the file containing station list
     integer, pointer :: out_ID(:)     ! the indices of the outlet subbasins whether the stations are located
     character(len=80), pointer :: out_name(:)  ! the name of the outlets  
     character(len=80) :: curOutlet    ! the name of the current outlet
   
     integer :: OPT_inund            ! Options for inundation, 0=inundation off, 1=inundation on
     integer :: OPT_trueDW           ! Options for diffusion wave channel routing method:
                                     !     1 -- True diffusion wave method for channel routing;
                                     !     2 -- False diffusion wave method: use riverbed slope as the surrogate for water surface slope. ( This is 
                                     !          a temporary treatment before the downstream-channel information can be retrieved. )
     integer :: OPT_calcNr           ! Options to calculate channel Manning roughness coefficients : 
                                     !     1 -- use channel depth (Luo et al. 2017 GMD); 
                                     !     2 -- use channel depth and exponent of 1/3 (Getirana et al. 2012 JHM); 
                                     !     3 -- use channel width (Decharme et al. 2010 JHM); 
                                     !     4 -- use one uniform value. 
                                     !     (Please see MOSARTinund_preProcs.F90 for references.)
     real(r8) :: nr_max              ! Max Manning coefficient for channels (when OPT_calcNr = 1, 2, 3) ( s*m^(-1/3) ).
     real(r8) :: nr_min              ! Min Manning coefficient for channels (when OPT_calcNr = 1, 2, 3) ( s*m^(-1/3) ).
     real(r8) :: nr_uniform          ! The uniform Manning coefficient for all channels (when OPT_calcNr = 4) ( s*m^(-1/3) ).
     real(r8) :: rdepth_max          ! Max channel depth (used when OPT_calcNr = 1, 2) (m).
     real(r8) :: rdepth_min          ! Min channel depth (used when OPT_calcNr = 1, 2) (m). 
     real(r8) :: rwidth_max          ! Max channel width (used when OPT_calcNr = 3) (m).
     real(r8) :: rwidth_min          ! Min channel width (used when OPT_calcNr = 3) (m). 
     real(r8) :: rslp_assume         ! Use this assumed riverbed slope when the input riverbed slope <= zero (dimensionless).
     real(r8) :: minL_tribRouting    ! Min tributary channel length for using tributary routing (m).  
  
     ! --------------------------------- 
     ! The following parameters are for the inundation scheme :
     ! --------------------------------- 
     integer :: OPT_elevProf         ! Options of elevation profile data: 1 -- Use real data; 2 -- Use hypothetical values.
     integer :: npt_elevProf         ! Number of dividing points in the elevation profile.
     real(r8) :: threshold_slpRatio  ! Threshold of the ratio of the lowest section's slope to the second lowest section's slope in 
                                     ! the elevation profile (used to alleviate the effect of DEM pits on elevation profiles).

     ! Elevations in the hypothetical elevation profile (m):
     ! (1) Gentle slope:
     real(r8) :: e_eprof_std(12) = (/ 0.0_r8, 1.0_r8, 2.0_r8, 3.0_r8, 4.0_r8, 5.0_r8, 7.0_r8, 11.0_r8, 19.0_r8, 35.0_r8, 75.0_r8, 10000.0_r8 /)

     ! (2) Steep slope:
     !real(r8) :: e_eprof_std(12) = (/ 0.0_r8, 15.0_r8, 35.0_r8, 60.0_r8, 90.0_r8, 125.0_r8, 165.0_r8, 205.0_r8, 245.0_r8, 285.0_r8, 325.0_r8, 10000.0_r8 /)

  end type Tcontrol
  
  ! --- Topographic and geometric properties, applicable for both grid- and subbasin-based representations
  public :: Tspatialunit
  type Tspatialunit
     ! grid properties
     integer , pointer :: mask(:)      ! mosart mask of mosart cell, 0=null, 1=land with dnID, 2=outlet
     integer , pointer :: ID0(:)         
     real(r8), pointer :: lat(:)       ! latitude of the centroid of the cell
     real(r8), pointer :: lon(:)       ! longitude of the centroid of the cell
     real(r8), pointer :: area(:)      ! area of local cell, [m2]
     real(r8), pointer :: areaTotal(:) ! total upstream drainage area, [m2]
     real(r8), pointer :: areaTotal2(:)! computed total upstream drainage area, [m2]
     real(r8), pointer :: rlenTotal(:) ! length of all reaches, [m]
     real(r8), pointer :: Gxr(:)       ! drainage density within the cell, [1/m]
     real(r8), pointer :: frac(:)      ! fraction of cell included in the study area, [-]
     real(r8), pointer :: domainfrac(:)! fraction of cell included in the study area from domain file, [-]
     logical , pointer :: euler_calc(:)! flag for calculating tracers in euler
     integer , pointer :: ocn_rof_coupling_ID(:)  ! ocn rof 2-way coupling ID, 0=off, 1=on
     real(r8), pointer :: vdatum_conversion(:)    ! ocn rof 2-way coupling vertical datum conversion

     ! hillslope properties
     real(r8), pointer :: nh(:)        ! manning's roughness of the hillslope (channel network excluded) 
     real(r8), pointer :: hslp(:)      ! slope of hillslope, [-]
     real(r8), pointer :: hslpsqrt(:)  ! sqrt of slope of hillslope, [-] 
     real(r8), pointer :: hlen(:)      ! length of hillslope within the cell, [m] 

     ! subnetwork channel properties
     real(r8), pointer :: tslp(:)      ! average slope of tributaries, [-]
     real(r8), pointer :: tslpsqrt(:)  ! sqrt of average slope of tributaries, [-] 
     real(r8), pointer :: tlen(:)      ! length of all sub-network reach within the cell, [m] 
     real(r8), pointer :: twidth(:)    ! bankfull width of the sub-reach, [m]
     real(r8), pointer :: nt(:)        ! manning's roughness of the subnetwork at hillslope  

     ! main channel properties
     real(r8), pointer :: rlen(:)      ! length of main river reach, [m]
     real(r8), pointer :: rslp(:)      ! slope of main river reach, [-]
     real(r8), pointer :: rslpsqrt(:)  ! sqrt of slope of main river reach, [-] 
     real(r8), pointer :: rwidth(:)    ! bankfull width of main reach, [m]
     real(r8), pointer :: rwidth0(:)   ! total width of the flood plain, [m]
     real(r8), pointer :: rdepth(:)    ! bankfull depth of river cross section, [m]
     real(r8), pointer :: nr(:)        ! manning's roughness of the main reach
     integer , pointer :: dnID(:)      ! IDs of the downstream units, corresponding to the subbasin ID in the input table
     integer , pointer :: nUp(:)       ! number of upstream units, maximum 8
     integer , pointer :: iUp(:,:)     ! IDs of upstream units, corresponding to the subbasin ID in the input table

     integer , pointer :: indexDown(:) ! indices of the downstream units in the ID array. sometimes subbasins IDs may not be continuous
     integer , pointer :: iDown(:)     ! indices of the downstream units, local
     integer , pointer :: nUp_dstrm(:)       ! number of upstream units, maximum 8
  
     integer , pointer :: numDT_r(:)   ! for a main reach, the number of sub-time-steps needed for numerical stability
     integer , pointer :: numDT_t(:)   ! for a subnetwork reach, the number of sub-time-steps needed for numerical stability
     real(r8), pointer :: phi_r(:)     ! the indicator used to define numDT_r
     real(r8), pointer :: phi_t(:)     ! the indicator used to define numDT_t
     real(r8), pointer :: rlen_dstrm(:)  ! Length of downstream channel (m).
     real(r8), pointer :: rslp_dstrm(:)  ! Bed slope of downstream channel (dimensionless).

!#ifdef INCLUDE_INUND   
     !real(r8), pointer :: rlen_dstrm(:)  ! Length of downstream channel (m).
     !real(r8), pointer :: rslp_dstrm(:)  ! Bed slope of downstream channel (dimensionless).
     real(r8), pointer :: wr_bf(:)       ! Water volume in the bankfull channel (i.e., channel storage capacity) (m^3).
  
     ! --------------------------------- 
     ! Parameters related to elevation profiles : 
     ! --------------------------------- 
     !real(r8) :: e_eprof_std(12)         ! Elevations in the hypothetical elevation profile (m).
     real(r8), pointer :: e_eprof_in2(:,:)  ! Absolute elevations in the input elevation profiles (m).
     real(r8), pointer :: a_eprof(:,:)   ! Area fractions of computation unit (grid cell or subbasin) in the elevation profiles (dimensionless).
     real(r8), pointer :: e_eprof(:,:)   ! Absolute elevations in the elevation profiles used in computation (m).
     real(r8), pointer :: a_chnl(:)      ! = channel area / computation unit area (dimensionless).
     real(r8), pointer :: e_chnl(:)      ! Channel banktop elevation (m).
     integer , pointer :: ipt_bl_bktp(:) ! The index of the point right below the banktop in the elevation profile.
    
     ! --------------------------------- 
     ! Parameters related to adjusted elevation profiles (where the section below banktop is replaced with level line) : 
     ! --------------------------------- 
     real(r8), pointer :: a_eprof3(:,:)  ! Area fractions of computation unit (grid cell or subbasin) in adjusted elevation profiles (dimensionless).
     real(r8), pointer :: e_eprof3(:,:)  ! Relative elevations in adjusted elevation profiles used in computation (m).
     integer , pointer :: npt_eprof3(:)  ! Number of points in the adjusted elevation profile.
     real(r8), pointer :: s_eprof3(:,:)  ! Total volume below the level through a point in the adjusted elevation profile (i.e., the
                                         ! volume between channel banktop and an elevation of the adjusted elevation profile) (m^3).

     ! --------------------------------- 
     ! In the inundation calculation, a quadratic equation is solved to derive the water depth based on water volume. 
     ! The following coefficients are for the quadratic equation. 
     ! (Please find more information in "MOSARTinund_Core_MOD.F90".)
     ! --------------------------------- 
     real(r8), pointer :: alfa3(:,:)     ! Coefficient (1/m).
     real(r8), pointer :: p3(:,:)        ! Coefficient (m).
     real(r8), pointer :: q3(:,:)        ! Coefficient (m^2).

  end type Tspatialunit

  ! status and flux variables
  public :: TstatusFlux
  type TstatusFlux
     ! hillsloope
     !! states
     real(r8), pointer :: wh(:,:)      ! storage of surface water, [m]
     real(r8), pointer :: dwh(:,:)     ! change of water storage, [m/s]              ( Unit is "m" when the inundation scheme is on. --Inund. )
     real(r8), pointer :: yh(:,:)      ! depth of surface water, [m]                  ( Not used when the inundation scheme is on; It is same as wh(:,:). --Inund. )
     real(r8), pointer :: wsat(:,:)    ! storage of surface water within saturated area at hillslope [m]
     real(r8), pointer :: wunsat(:,:)  ! storage of surface water within unsaturated area at hillslope [m]
     real(r8), pointer :: qhorton(:,:) ! Infiltration excess runoff generated from hillslope, [m/s]
     real(r8), pointer :: qdunne(:,:)  ! Saturation excess runoff generated from hillslope, [m/s]
     real(r8), pointer :: qsur(:,:)    ! Surface runoff generated from hillslope, [m/s]
     real(r8), pointer :: qsub(:,:)    ! Subsurface runoff generated from hillslope, [m/s]
     real(r8), pointer :: qdto(:,:)    ! Direct to Ocean runoff, [m/s]
     real(r8), pointer :: qgwl(:,:)    ! gwl runoff term from glacier, wetlands and lakes, [m/s]
     !! fluxes
     real(r8), pointer :: ehout(:,:)   ! overland flow from hillslope into the sub-channel, [m/s]          ( Note: outflow is negative. --Inund. )
     real(r8), pointer :: asat(:,:)    ! saturated area fraction from hillslope, [-]
     real(r8), pointer :: esat(:,:)    ! evaporation from saturated area fraction at hillslope, [m/s]
     real(r8), pointer :: ehexchange(:,:)    ! net influx from hillslope (soil) storage into overland flow, e.g., soil erosion [kg/s], or runoff re-infiltration [m/s]
     real(r8), pointer :: ehexch_avg(:,:)    ! net influx from hillslope (soil) storage into overland flow, e.g., soil erosion [kg/s], or runoff re-infiltration [m/s], average

     ! subnetwork channel
     !! states
     real(r8), pointer :: tarea(:,:)   ! area of channel water surface, [m2]
     real(r8), pointer :: wt(:,:)      ! storage of surface water, [m3]
     real(r8), pointer :: dwt(:,:)     ! change of water storage, [m3]
     real(r8), pointer :: yt(:,:)      ! water depth, [m]
     real(r8), pointer :: mt(:,:)      ! cross section area, [m2]
     real(r8), pointer :: rt(:,:)      ! hydraulic radii, [m]
     real(r8), pointer :: pt(:,:)      ! wetness perimeter, [m]
     real(r8), pointer :: vt(:,:)      ! flow velocity, [m/s]
     real(r8), pointer :: tt(:,:)      ! mean travel time of the water within the channel, [s]
     real(r8), pointer :: conc_t(:,:)  ! MOSART sub-network concentration of tracers such as sediment, C,N,P (kg/m**3)
     real(r8), pointer :: wt_al(:,:)   ! MOSART sub-network channel active layer storage (kg for mud and sand sediment)
     real(r8), pointer :: dwt_al(:,:)   ! change of MOSART sub-network channel active layer storage (kg for mud and sand sediment)
     !! fluxes
     real(r8), pointer :: tevap(:,:)   ! evaporation, [m/s]
     real(r8), pointer :: etin(:,:)    ! lateral inflow from hillslope, including surface and subsurface runoff generation components, [m3/s]
     real(r8), pointer :: etout(:,:)   ! discharge from sub-network into the main reach, [m3/s]          ( Note: outflow is negative. --Inund. )
     real(r8), pointer :: qdem(:,:)    ! irrigation demand [m/s] !added by Yuna 1/29/2018
     real(r8), pointer :: etexchange(:,:)    ! net influx from channel bank storage into channel, e.g., sediment erosion/deposition [kg/s], or groundwater/river water exchange [m/s]
     real(r8), pointer :: etexch_avg(:,:)    ! net influx from channel bank storage into channel, e.g., sediment erosion/deposition [kg/s], or groundwater/river water exchange [m/s], average

     ! main channel
     !! states
     real(r8), pointer :: rarea(:,:)   ! area of channel water surface, [m2] 
     real(r8), pointer :: wr(:,:)      ! storage of surface water, [m3]
     real(r8), pointer :: dwr(:,:)     ! change of water storage, [m3]
     real(r8), pointer :: yr(:,:)      ! water depth. [m]
     real(r8), pointer :: mr(:,:)      ! cross section area, [m2]
     real(r8), pointer :: rr(:,:)      ! hydraulic radius, [m]
     real(r8), pointer :: pr(:,:)      ! wetness perimeter, [m]
     real(r8), pointer :: vr(:,:)      ! flow velocity, [m/s]
     real(r8), pointer :: tr(:,:)      ! mean travel time of the water within the channel, [s]
     real(r8), pointer :: conc_r(:,:)  ! MOSART main channel concentration of tracers such as sediment, C,N,P (kg/m**3)
     real(r8), pointer :: wr_al(:,:)   ! MOSART main channel active layer storage (kg for mud and sand sediment)
     real(r8), pointer :: dwr_al(:,:)   ! change of MOSART main channel active layer storage (kg for mud and sand sediment)
     real(r8), pointer :: rslp_energy(:)! energy slope of channel water surface [-]
     real(r8), pointer :: wr_dstrm(:,:)  ! Downstream-channel water volume  (to constrain large upward flow from downstream channel to current channel ) (m^3 or kg).
     real(r8), pointer :: yr_dstrm(:)  ! Downstream-channel water depth (m).
     real(r8), pointer :: conc_r_dstrm(:,:) ! Downstream-channel concentration of BGC fluxes (kg/m**3).

     !! exchange fluxes
     real(r8), pointer :: erlg(:,:)    ! evaporation, [m/s]
     real(r8), pointer :: erlateral(:,:) ! lateral flow from hillslope, including surface and subsurface runoff generation components, [m3/s]
     real(r8), pointer :: erin(:,:)    ! inflow from upstream links, [m3/s]
     real(r8), pointer :: erout(:,:)   ! outflow into downstream links, [m3/s] (negative is out)
     real(r8), pointer :: eroup_lagi(:,:) ! outflow into downstream links from previous timestep, [m3/s]  (Note: average channel outflow in one MOSART sub-step. --Inund.)
     real(r8), pointer :: eroup_lagf(:,:) ! outflow into downstream links from current timestep, [m3/s]
     real(r8), pointer :: erowm_regi(:,:) ! initial outflow before dam regulation at current timestep, [m3/s]
     real(r8), pointer :: erowm_regf(:,:) ! final outflow after dam regulation at current timestep, [m3/s]
     real(r8), pointer :: eroutUp(:,:) ! outflow sum of upstream gridcells, instantaneous (m3/s)
     real(r8), pointer :: eroutUp_avg(:,:) ! outflow sum of upstream gridcells, average [m3/s]
     real(r8), pointer :: erlat_avg(:,:) ! erlateral average [m3/s]
     real(r8), pointer :: flow(:,:)    ! streamflow from the outlet of the reach, [m3/s] (positive is out)
     real(r8), pointer :: erin1(:,:)   ! inflow from upstream links during previous step, used for Muskingum method, [m3/s]
     real(r8), pointer :: erin2(:,:)   ! inflow from upstream links during current step, used for Muskingum method, [m3/s]
     real(r8), pointer :: ergwl(:,:)   ! flux item for the adjustment of water balance residual in glacie, wetlands and lakes dynamics [m3/s]
     real(r8), pointer :: erexchange(:,:)    ! net influx from channel bank storage into channel, e.g., sediment erosion/deposition [kg/s], or groundwater/river water exchange [m/s]
     real(r8), pointer :: erexch_avg(:,:)    ! net influx from channel bank storage into channel, e.g., sediment erosion/deposition [kg/s], or groundwater/river water exchange [m/s], average
     real(r8), pointer :: erin_dstrm(:,:)    ! total riverine inflow into the downstream link[m3/s]

     !! for Runge-Kutta algorithm
     real(r8), pointer :: wrtemp(:,:)  ! temporary storage item, for 4th order Runge-Kutta  algorithm;
     real(r8), pointer :: erintemp(:,:)
     real(r8), pointer :: erouttemp(:,:)
     real(r8), pointer :: k1(:,:)
     real(r8), pointer :: k2(:,:)
     real(r8), pointer :: k3(:,:)
     real(r8), pointer :: k4(:,:)
   
    !real(r8), pointer :: wr_ini(:)     ! Channel water volume at beginning of step (m^3).
    !real(r8), pointer :: yr_ini(:)     ! Channel water depth at beginning of step (m).
    real(r8), pointer :: wf_ini(:)      ! Floodplain water volume at beginning of step (m^3).
    real(r8), pointer :: hf_ini(:)      ! Floodplain max water depth (i.e., elevation difference between water level and channel banktop) at beginning of step (m).
    real(r8), pointer :: ff_ini(:)      ! Floodplain water area fraction at beginning of step (dimensionless).
    real(r8), pointer :: ffunit_ini(:)      ! Flooded water area fraction at beginning of step (dimensionless).

    real(r8), pointer :: netchange(:)   ! Amount of channel--floodplain exchange during one subcycle in MOSART timestep (positive: flow from channel to floodplain; vice versa ) (m^3).
    real(r8), pointer :: se_rf(:)       ! Amount of channel--floodplain exchange during one MOSART timestep(positive: flow from channel to floodplain; vice versa ) (m^3).
    real(r8), pointer :: ff_unit(:)       ! = area of inundated area (including channel area) divided by the computation-unit total area (dimensionless).

    real(r8), pointer :: ff_fp(:)       ! = area of inundated floodplain (not including channel area) divided by the computation-unit total area (dimensionless).
    real(r8), pointer :: fa_fp(:)       ! Area of inundated floodplain (not including channel area) (m^2).    
    real(r8), pointer :: wr_exchg(:)    ! Channel water volume after channel--floodplain exchange (m^3).
    real(r8), pointer :: wr_exchg_dstrm(:)  ! Downstream-channel water volume after channel--floodplain exchange (to 
                                            ! constrain large upward flow from downstream channel to current channel ) (m^3).
    real(r8), pointer :: yr_exchg(:)    ! Channel water depth after channel--floodplain exchange (m).
    real(r8), pointer :: yr_exchg_dstrm(:)  ! Downstream-channel water depth after channel--floodplain exchange (m).
    real(r8), pointer :: wf_exchg(:)    ! Floodplain water volume after channel--floodplain exchange (m^3).
    real(r8), pointer :: hf_exchg(:)    ! Floodplain max water depth after channel--floodplain exchange (m).     
    !real(r8), pointer :: delta_wr(:)   ! Change of channel water volume during channel routing (m^3).
    real(r8), pointer :: wr_rtg(:)      ! Channel water volume after channel routing (m^3).
    real(r8), pointer :: yr_rtg(:)      ! Channel water depth after channel routing (m).
   
  end type TstatusFlux
  !== Hongyi


  ! heat status and flux variables
  public :: TstatusFlux_heat
  type TstatusFlux_heat
      ! overall
      real(r8), pointer :: forc_t(:)      ! atmospheric temperature (Kelvin)
      real(r8), pointer :: forc_pbot(:)   ! atmospheric pressure (Pa)
      real(r8), pointer :: forc_vp(:)     ! atmospheric vapor pressure (Pa)
      real(r8), pointer :: forc_wind(:)   ! atmospheric wind speed (m/s)
      real(r8), pointer :: forc_lwrad(:)  ! downward infrared (longwave) radiation (W/m**2)
      real(r8), pointer :: forc_solar(:)  ! atmospheric incident solar (shortwave) radiation (W/m**2)
      real(r8), pointer :: coszen(:)      ! Cosine of Zenith angle (-)
      ! hillsloope
      !! states
      real(r8), pointer :: Tqsur(:)       ! temperature of surface runoff, [K]
      real(r8), pointer :: Tqsub(:)       ! temperature of subsurface runoff, [K]
      real(r8), pointer :: Tqice(:)       ! temperature of ice flow, [K]
      !! fluxes
      
      ! subnetwork channel
      !! states
      real(r8), pointer :: Tt(:)            ! temperature of subnetwork water, [K]
      !! fluxes
      !real(r8), pointer :: Ha_t(:)       ! advective heat flux through the subnetwork, [Watt]
      real(r8), pointer :: Ha_h2t(:)      ! advective heat flux from hillslope into the subnetwork, [Watt]
      real(r8), pointer :: Ha_t2r(:)      ! advective heat flux from subnetwork channel into the main channel, [Watt]
      real(r8), pointer :: Ha_lateral(:)  ! average advective heat flux from subnetwork channel into the main channel, [Watt], corresponding to TRunoff%erlateral
      real(r8), pointer :: Hs_t(:)        ! net solar short-wave radiation, [Watt]
      real(r8), pointer :: Hl_t(:)        ! net solar long-wave radiation, [Watt]
      real(r8), pointer :: He_t(:)        ! flux of latent heat, [Watt]
      real(r8), pointer :: Hh_t(:)        ! flux of sensible heat, [Watt]
      real(r8), pointer :: Hc_t(:)        ! conductive heat flux at the streambed, [Watt]
      real(r8), pointer :: deltaH_t(:)    ! net heat exchange with surroundings, [J]
      real(r8), pointer :: deltaM_t(:)    ! net heat change due to inflow, [J]

      ! main channel
      !! states
      real(r8), pointer :: Tr(:)            ! temperature of main channel water, [K]
      !! fluxes
      !real(r8), pointer :: Ha_r(:)       ! advective heat flux through the main channel, [Watt]
      real(r8), pointer :: Ha_rin(:)      ! advective heat flux from upstream into the main channel, [Watt]
      real(r8), pointer :: Ha_rout(:)     ! advective heat flux to downstream channel, [Watt]
      real(r8), pointer :: Ha_eroutUp(:) ! outflow sum of upstream gridcells, instantaneous (Watt)
      real(r8), pointer :: Ha_eroutUp_avg(:) ! outflow sum of upstream gridcells, average [Watt]
      real(r8), pointer :: Ha_erlat_avg(:) ! erlateral average [Watt]
      real(r8), pointer :: Hs_r(:)        ! net solar short-wave radiation, [Watt]
      real(r8), pointer :: Hl_r(:)        ! net solar long-wave radiation, [Watt]
      real(r8), pointer :: He_r(:)        ! flux of latent heat, [Watt]
      real(r8), pointer :: Hh_r(:)        ! flux of sensible heat, [Watt]
      real(r8), pointer :: Hc_r(:)        ! conductive heat flux at the streambed, [Watt]
      real(r8), pointer :: deltaH_r(:)    ! net heat exchange with surroundings, [J]
      real(r8), pointer :: deltaM_r(:)    ! net heat change due to inflow, [J]

      real(r8), pointer :: Tt_avg(:)      ! average temperature of subnetwork channel water, [K], for output purpose
      real(r8), pointer :: Tr_avg(:)      ! average temperature of main channel water, [K], for output purpose
      
  end type TstatusFlux_heat

 
  ! parameters to be calibrated. Ideally, these parameters are supposed to be uniform for one region
  public :: Tparameter
  type Tparameter
     real(r8), pointer :: c_nr(:)       ! coefficient to adjust the manning's roughness of channels
     real(r8), pointer :: c_nh(:)       ! coefficient to adjust the manning's roughness of overland flow across hillslopes
     real(r8), pointer :: c_twid(:)     ! coefficient to adjust the width of sub-reach channel
     real(r8), pointer :: t_alpha(:)    ! alpha parameter in air-water temperature relationship (S-curve)
     real(r8), pointer :: t_beta(:)     ! beta parameter in air-water temperature relationship (S-curve)
     real(r8), pointer :: t_gamma(:)    ! gamma parameter in air-water temperature relationship (S-curve)
     real(r8), pointer :: t_mu(:)       ! mu parameter in air-water temperature relationship (S-curve)
  end type Tparameter 

  !== Hongyi
  type (Tcontrol)    , public :: Tctl
  type (Tspatialunit), public :: TUnit
  type (TstatusFlux) , public :: TRunoff
  type (TstatusFlux_heat), public :: THeat
  type (Tparameter)  , public :: TPara
  !== Hongyi

  type (runoff_flow) , public :: rtmCTL

  public :: RunoffInit

contains

  subroutine RunoffInit(begr, endr, numr)

    integer, intent(in) :: begr, endr, numr

    integer :: ier

    allocate(rtmCTL%runoff(begr:endr,nt_rtm),     &
             rtmCTL%dvolrdt(begr:endr,nt_rtm),    &
             rtmCTL%runofflnd(begr:endr,nt_rtm),  &
             rtmCTL%dvolrdtlnd(begr:endr,nt_rtm), &
             rtmCTL%runoffocn(begr:endr,nt_rtm),  &
             rtmCTL%dvolrdtocn(begr:endr,nt_rtm), &
             rtmCTL%runofftot(begr:endr,nt_rtm),  &
             rtmCTL%area(begr:endr),              &
             rtmCTL%volr(begr:endr,nt_rtm),       &
             rtmCTL%lonc(begr:endr),              &
             rtmCTL%latc(begr:endr),              &
             rtmCTL%dsig(begr:endr),              &
             rtmCTL%rdsig(begr:endr),             &
             rtmCTL%iDown(begr:endr),             &
             rtmCTL%nUp(begr:endr),               &
             rtmCTL%nUp_dstrm(begr:endr),         &
             rtmCTL%iUp(begr:endr,8),             &
             rtmCTL%outletg(begr:endr),           &
             rtmCTL%routletg(begr:endr),          &
             rtmCTL%inundwf(begr:endr),           &
             rtmCTL%inundhf(begr:endr),           &
             rtmCTL%inundff(begr:endr),           &
             rtmCTL%inundffunit(begr:endr),       &
             rtmCTL%runofflnd_nt1(begr:endr),     &
             rtmCTL%runofflnd_nt2(begr:endr),     &
             rtmCTL%runofflnd_nt3(begr:endr),     &
             rtmCTL%runofflnd_nt4(begr:endr),     &
             rtmCTL%runoffocn_nt1(begr:endr),     &
             rtmCTL%runoffocn_nt2(begr:endr),     &
             rtmCTL%runoffocn_nt3(begr:endr),     &
             rtmCTL%runoffocn_nt4(begr:endr),     &
             rtmCTL%runofftot_nt1(begr:endr),     &
             rtmCTL%runofftot_nt2(begr:endr),     &
             rtmCTL%runoffdir_nt1(begr:endr),     &
             rtmCTL%runoffdir_nt2(begr:endr),     &
             rtmCTL%volr_nt1(begr:endr),          &
             rtmCTL%volr_nt2(begr:endr),          &
             rtmCTL%volr_nt3(begr:endr),          &
             rtmCTL%volr_nt4(begr:endr),          &
             rtmCTL%dvolrdtlnd_nt1(begr:endr),    &
             rtmCTL%dvolrdtlnd_nt2(begr:endr),    &
             rtmCTL%dvolrdtocn_nt1(begr:endr),    &
             rtmCTL%dvolrdtocn_nt2(begr:endr),    &
             rtmCTL%wr_nt1(begr:endr),            &
             rtmCTL%qsur_nt1(begr:endr),          &
             rtmCTL%qsur_nt2(begr:endr),          &
             rtmCTL%qsur_nt3(begr:endr),          &
             rtmCTL%qsub_nt1(begr:endr),          &
             rtmCTL%qsub_nt2(begr:endr),          &
             rtmCTL%qgwl_nt1(begr:endr),          &
             rtmCTL%qgwl_nt2(begr:endr),          &
             rtmCTL%qdto_nt1(begr:endr),          &
             rtmCTL%qdto_nt2(begr:endr),          &
             rtmCTL%qdem_nt1(begr:endr),          &
             rtmCTL%qdem_nt2(begr:endr),          &
             rtmCTL%mask(begr:endr),              &
             rtmCTL%rmask(begr:endr),             &
             rtmCTL%gindex(begr:endr),            &
             rtmCTL%rgindex(begr:endr),           &
             rtmCTL%fthresh(begr:endr),           &
             rtmCTL%flood(begr:endr),             &
             rtmCTL%direct(begr:endr,nt_rtm),     &
             rtmCTL%inundinf(begr:endr),          &
             rtmCTL%wh(begr:endr,nt_rtm),         &
             rtmCTL%wt(begr:endr,nt_rtm),         &
             rtmCTL%wr(begr:endr,nt_rtm),         &
             rtmCTL%mr(begr:endr,nt_rtm),         &
             rtmCTL%yr(begr:endr,nt_rtm),         &
             rtmCTL%pr(begr:endr,nt_rtm),         &
             rtmCTL%rr(begr:endr,nt_rtm),         &
             rtmCTL%erout(begr:endr,nt_rtm),      &
             rtmCTL%qsur(begr:endr,nt_rtm),       & 
             rtmCTL%qsub(begr:endr,nt_rtm),       &
             rtmCTL%qgwl(begr:endr,nt_rtm),       &
             rtmCTL%qdto(begr:endr,nt_rtm),       &
             rtmCTL%qdem(begr:endr,nt_rtm),       & 
             rtmCTL%yr_nt1(begr:endr),            &
             rtmCTL%ssh(begr:endr),               &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff local arrays'
       call shr_sys_abort
    end if

    if (data_bgc_fluxes_to_ocean_flag) then
      allocate(rtmCTL%concDIN(begr:endr),           &
               rtmCTL%concDIP(begr:endr),           &
               rtmCTL%concDON(begr:endr),           &
               rtmCTL%concDOP(begr:endr),           &
               rtmCTL%concDOC(begr:endr),           &
               rtmCTL%concPP(begr:endr),            &
               rtmCTL%concDSi(begr:endr),           &
               rtmCTL%concPOC(begr:endr),           &
               rtmCTL%concPN(begr:endr),            &
               rtmCTL%concDIC(begr:endr),           &
               rtmCTL%concFe(begr:endr),            &
               stat=ier)
      if (ier /= 0) then
         write(iulog,*)'Rtmini ERROR allocation of BGC local arrays'
         call shr_sys_abort
      end if
    end if

    rtmCTL%iDown(:) = 0
    rtmCTL%iUp(:,:) = 0
    rtmCTL%nUp(:) = 0
    rtmCTL%nUp_dstrm(:) = 0
    
    rtmCTL%runoff(:,:)     = 0._r8
    rtmCTL%runofflnd(:,:)  = spval
    rtmCTL%runoffocn(:,:)  = spval
    rtmCTL%runofftot(:,:)  = spval
    rtmCTL%dvolrdt(:,:)    = 0._r8
    rtmCTL%dvolrdtlnd(:,:) = spval
    rtmCTL%dvolrdtocn(:,:) = spval
    rtmCTL%volr(:,:)       = 0._r8
    rtmCTL%flood(:)        = 0._r8
    rtmCTL%direct(:,:)     = 0._r8
    rtmCTL%inundinf(:)     = 0._r8
    rtmCTL%inundwf(:)      = 0._r8
    rtmCTL%inundhf(:)      = 0._r8
    rtmCTL%inundff(:)      = 0._r8
    rtmCTL%inundffunit(:)  = 0._r8
    rtmCTL%wh(:,:)         = 0._r8
    rtmCTL%wt(:,:)         = 0._r8
    rtmCTL%wr(:,:)         = 0._r8
    rtmCTL%qsur(:,:)       = 0._r8
    rtmCTL%qsub(:,:)       = 0._r8
    rtmCTL%qgwl(:,:)       = 0._r8
    rtmCTL%qdto(:,:)       = 0._r8
    rtmCTL%qdem(:,:)       = 0._r8
    if (data_bgc_fluxes_to_ocean_flag) then
      rtmCTL%concDIN(:)      = 0._r8
      rtmCTL%concDIP(:)      = 0._r8
      rtmCTL%concDON(:)      = 0._r8
      rtmCTL%concDOP(:)      = 0._r8
      rtmCTL%concDOC(:)      = 0._r8
      rtmCTL%concPP(:)       = 0._r8
      rtmCTL%concDSi(:)      = 0._r8
      rtmCTL%concPOC(:)      = 0._r8
      rtmCTL%concPN(:)       = 0._r8
      rtmCTL%concDIC(:)      = 0._r8
      rtmCTL%concFe(:)       = 0._r8
    end if
    
    if (heatflag) then
      allocate(rtmCTL%Tqsur(begr:endr),                 &
               rtmCTL%Tqsub(begr:endr),                 &
               rtmCTL%Tt(begr:endr),                    &
               rtmCTL%Tr(begr:endr),                    &
               rtmCTL%Ha_rout(begr:endr),               &
               rtmCTL%templand_Tqsur(begr:endr),        &
               rtmCTL%templand_Tqsub(begr:endr),        &
               rtmCTL%templand_Ttrib(begr:endr),        &
               rtmCTL%templand_Tchanr(begr:endr),       &
               rtmCTL%templand_Tqsur_nt1(begr:endr),    &
               rtmCTL%templand_Tqsub_nt1(begr:endr),    &
               rtmCTL%templand_Ttrib_nt1(begr:endr),    &
               rtmCTL%templand_Tchanr_nt1(begr:endr),   &
               rtmCTL%templand_Tqsur_nt2(begr:endr),    &
               rtmCTL%templand_Tqsub_nt2(begr:endr),    &
               rtmCTL%templand_Ttrib_nt2(begr:endr),    &
               rtmCTL%templand_Tchanr_nt2(begr:endr),   &
               stat=ier)
      if (ier /= 0) then
         write(iulog,*)'Rtmini ERROR allocation of runoff local arrays'
         call shr_sys_abort
      end if
      
      rtmCTL%Tqsur(:)        = 273.15_r8
      rtmCTL%Tqsub(:)        = 273.15_r8
      rtmCTL%templand_Tqsur(:)  = spval
      rtmCTL%templand_Tqsub(:)  = spval
      rtmCTL%templand_Ttrib(:)  = spval
      rtmCTL%templand_Tchanr(:) = spval
      
    end if

  end subroutine RunoffInit

end module RunoffMod
