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
  use RtmVar         , only : iulog, spval
  use rof_cpl_indices, only : nt_rtm

! !PUBLIC TYPES:
  implicit none
  private

  type(mct_gsmap),public :: gsmap_r       ! gsmap for mosart decomposition

  type(mct_sMatP),public :: sMatP_dnstrm  ! sparse matrix plus for downstream advection
  type(mct_avect),public :: avsrc_dnstrm  ! src avect for SM mult downstream advection
  type(mct_avect),public :: avdst_dnstrm  ! dst avect for SM mult downstream advection

  type(mct_sMatP),public :: sMatP_direct  ! sparse matrix plus for direct to outlet flow
  type(mct_avect),public :: avsrc_direct  ! src avect for SM mult direct to outlet flow
  type(mct_avect),public :: avdst_direct  ! dst avect for SM mult direct to outlet flow

  type(mct_sMatP),public :: sMatP_eroutUp ! sparse matrix plus for eroutUp calc
  type(mct_avect),public :: avsrc_eroutUp ! src avect for SM mult eroutUp calc
  type(mct_avect),public :: avdst_eroutUp ! dst avect for SM mult eroutUp calc

  public :: runoff_flow
  type runoff_flow
     !    - local initialization
     real(r8), pointer :: lonc(:)          ! lon of cell
     real(r8), pointer :: latc(:)          ! lat of cell
     real(r8), pointer :: area(:)          ! area of cell
     integer , pointer :: gindex(:)        ! global index consistent with map file
     integer , pointer :: dsig(:)          ! downstream index, global index
     integer , pointer :: outletg(:)       ! outlet index, global index

     !    - global 
     integer , pointer :: mask(:)          ! general mask of cell 1=land, 2=ocean, 3=outlet
     real(r8), pointer :: rlon(:)          ! rtm longitude list, 1d
     real(r8), pointer :: rlat(:)          ! rtm latitude list, 1d
     real(r8)          :: totarea          ! global area
     integer           :: numr             ! rtm gdc global number of cells

     !    - local
     integer           :: begr,endr        ! local start/stop indices
     integer           :: lnumr            ! local number of cells

     !    - local
     real(r8), pointer :: runofflnd(:,:)   ! runoff masked for land (m3 H2O/s)
     real(r8), pointer :: runoffocn(:,:)   ! runoff masked for ocn  (m3 H2O/s)
     real(r8), pointer :: runofftot(:,:)   ! total runoff masked for ocn  (m3 H2O/s)
     real(r8), pointer :: dvolrdt(:,:)     ! RTM change in storage (mm/s)
     real(r8), pointer :: dvolrdtlnd(:,:)  ! dvolrdt masked for land (mm/s)
     real(r8), pointer :: dvolrdtocn(:,:)  ! dvolrdt masked for ocn  (mm/s)
     real(r8), pointer :: volr(:,:)        ! RTM storage (m3)
     real(r8), pointer :: fthresh(:)       ! RTM water flood threshold

     !    - restarts
     real(r8), pointer :: wh(:,:)          ! MOSART hillslope surface water storage (m)
     real(r8), pointer :: wt(:,:)          ! MOSART sub-network water storage (m3)
     real(r8), pointer :: wr(:,:)          ! MOSART main channel water storage (m3)
     real(r8), pointer :: erout(:,:)       ! MOSART flow out of the main channel, instantaneous (m3/s)

     ! inputs
     real(r8), pointer :: qsur(:,:)        ! coupler surface forcing [m3/s]
     real(r8), pointer :: qsub(:,:)        ! coupler subsurface forcing [m3/s]
     real(r8), pointer :: qgwl(:,:)        ! coupler glacier/wetland/lake forcing [m3/s]
     real(r8), pointer :: qdto(:,:)        ! coupler diret-to-ocean forcing [m3/s]

     !    - outputs
     real(r8), pointer :: flood(:)         ! coupler return flood water sent back to clm [m3/s]
     real(r8), pointer :: runoff(:,:)      ! coupler return mosart basin derived flow [m3/s]
     real(r8), pointer :: direct(:,:)      ! coupler return direct flow [m3/s]

     !    - history (currently needed)
     real(r8), pointer :: runofflnd_nt1(:)
     real(r8), pointer :: runofflnd_nt2(:)
     real(r8), pointer :: runoffocn_nt1(:)
     real(r8), pointer :: runoffocn_nt2(:)
     real(r8), pointer :: runofftot_nt1(:)
     real(r8), pointer :: runofftot_nt2(:)
     real(r8), pointer :: runoffdir_nt1(:)
     real(r8), pointer :: runoffdir_nt2(:)
     real(r8), pointer :: dvolrdtlnd_nt1(:)
     real(r8), pointer :: dvolrdtlnd_nt2(:)
     real(r8), pointer :: dvolrdtocn_nt1(:)
     real(r8), pointer :: dvolrdtocn_nt2(:)
     real(r8), pointer :: volr_nt1(:)
     real(r8), pointer :: volr_nt2(:)
     real(r8), pointer :: qsur_nt1(:)
     real(r8), pointer :: qsur_nt2(:)
     real(r8), pointer :: qsub_nt1(:)
     real(r8), pointer :: qsub_nt2(:)
     real(r8), pointer :: qgwl_nt1(:)
     real(r8), pointer :: qgwl_nt2(:)
     real(r8), pointer :: qdto_nt1(:)
     real(r8), pointer :: qdto_nt2(:)

  end type runoff_flow

  
  !== Hongyi
  ! constrol information 
  public :: Tcontrol
  type Tcontrol
     integer  :: NUnit            ! numer of Grides in the model domain, which is equal to the number of cells, nrows*ncols
     integer  :: NSTART           ! the # of the time step to start the routing. Previous NSTART - 1 steps will be passed over.
     integer  :: NSTEPS           ! number of time steps specified in the modeling
     integer  :: NWARMUP          ! time steps for model warming up
     real(r8) :: DATAH            ! time step of runoff generation in second provided by the user
     integer  :: Num_dt           ! number of sub-steps within the current step interval, 
                                  ! i.e., if the time step of the incoming runoff data is 3-hr, and num_dt is set to 10, 
                                  ! then deltaT = 3*3600/10 = 1080 seconds
     real(r8) :: DeltaT           ! Time step in seconds 
     integer  :: DLevelH2R        ! The base number of channel routing sub-time-steps within one hillslope routing step. 
                                  ! Usually channel routing requires small time steps than hillslope routing.
     integer  :: DLevelR          ! The number of channel routing sub-time-steps at a higher level within one channel routing step at a lower level. 
     integer  :: Restart          ! flag, Restart=1 means starting from the state of last run, =0 means starting from model-inset initial state.
     integer  :: RoutingMethod    ! Flag for routing methods. 1 --> variable storage method from SWAT model; 2 --> Muskingum method?
     integer  :: RoutingFlag      ! Flag for whether including hillslope and sub-network routing. 1--> include routing through hillslope, sub-network and main channel; 0--> main channel routing only.
 
     character(len=100) :: baseName    ! name of the case study, e.g., columbia
     character(len=200) :: ctlFile     ! the name of the control file
     character(len=100) :: ctlPath     ! the path of the control file
     character(len=200) :: paraFile    ! the path of the parameter files
     character(len=100) :: paraPath    ! the path of the parameter files
     character(len=100) :: runoffPath  ! the path of the runoff data
     character(len=100) :: outPath     ! the path of the output file(s)
     integer :: numStation             ! number of basins to be simulated
     character(len=200) :: staListFile ! name of the file containing station list
     integer, pointer :: out_ID(:)     ! the indices of the outlet subbasins whether the stations are located
     character(len=80), pointer :: out_name(:)  ! the name of the outlets  
     character(len=80) :: curOutlet    ! the name of the current outlet
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
     logical , pointer :: euler_calc(:) ! flag for calculating tracers in euler

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
  
     integer , pointer :: numDT_r(:)   ! for a main reach, the number of sub-time-steps needed for numerical stability
     integer , pointer :: numDT_t(:)   ! for a subnetwork reach, the number of sub-time-steps needed for numerical stability
     real(r8), pointer :: phi_r(:)     ! the indicator used to define numDT_r
     real(r8), pointer :: phi_t(:)     ! the indicator used to define numDT_t
  end type Tspatialunit

  ! status and flux variables
  public :: TstatusFlux
  type TstatusFlux
     ! hillsloope
     !! states
     real(r8), pointer :: wh(:,:)      ! storage of surface water, [m]
     real(r8), pointer :: dwh(:,:)     ! change of water storage, [m/s]
     real(r8), pointer :: yh(:,:)      ! depth of surface water, [m]
     real(r8), pointer :: wsat(:,:)    ! storage of surface water within saturated area at hillslope [m]
     real(r8), pointer :: wunsat(:,:)  ! storage of surface water within unsaturated area at hillslope [m]
     real(r8), pointer :: qhorton(:,:) ! Infiltration excess runoff generated from hillslope, [m/s]
     real(r8), pointer :: qdunne(:,:)  ! Saturation excess runoff generated from hillslope, [m/s]
     real(r8), pointer :: qsur(:,:)    ! Surface runoff generated from hillslope, [m/s]
     real(r8), pointer :: qsub(:,:)    ! Subsurface runoff generated from hillslope, [m/s]
     real(r8), pointer :: qdto(:,:)    ! Direct to Ocean runoff, [m/s]
     real(r8), pointer :: qgwl(:,:)    ! gwl runoff term from glacier, wetlands and lakes, [m/s]
     !! fluxes
     real(r8), pointer :: ehout(:,:)   ! overland flow from hillslope into the sub-channel, [m/s]
     real(r8), pointer :: asat(:,:)    ! saturated area fraction from hillslope, [-]
     real(r8), pointer :: esat(:,:)    ! evaporation from saturated area fraction at hillslope, [m/s]

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
     !! fluxes
     real(r8), pointer :: tevap(:,:)   ! evaporation, [m/s]
     real(r8), pointer :: etin(:,:)    ! lateral inflow from hillslope, including surface and subsurface runoff generation components, [m3/s]
     real(r8), pointer :: etout(:,:)   ! discharge from sub-network into the main reach, [m3/s]

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
     !! exchange fluxes
     real(r8), pointer :: erlg(:,:)    ! evaporation, [m/s]
     real(r8), pointer :: erlateral(:,:) ! lateral flow from hillslope, including surface and subsurface runoff generation components, [m3/s]
     real(r8), pointer :: erin(:,:)    ! inflow from upstream links, [m3/s]
     real(r8), pointer :: erout(:,:)   ! outflow into downstream links, [m3/s]
     real(r8), pointer :: erout_prev(:,:) ! outflow into downstream links from previous timestep, [m3/s]
     real(r8), pointer :: eroutUp(:,:) ! outflow sum of upstream gridcells, instantaneous (m3/s)
     real(r8), pointer :: eroutUp_avg(:,:) ! outflow sum of upstream gridcells, average [m3/s]
     real(r8), pointer :: erlat_avg(:,:) ! erlateral average [m3/s]
     real(r8), pointer :: flow(:,:)    ! streamflow from the outlet of the reach, [m3/s]
     real(r8), pointer :: erin1(:,:)   ! inflow from upstream links during previous step, used for Muskingum method, [m3/s]
     real(r8), pointer :: erin2(:,:)   ! inflow from upstream links during current step, used for Muskingum method, [m3/s]
     real(r8), pointer :: ergwl(:,:)   ! flux item for the adjustment of water balance residual in glacie, wetlands and lakes dynamics [m3/s]

     !! for Runge-Kutta algorithm
     real(r8), pointer :: wrtemp(:,:)  ! temporary storage item, for 4th order Runge-Kutta  algorithm;
     real(r8), pointer :: erintemp(:,:)
     real(r8), pointer :: erouttemp(:,:)
     real(r8), pointer :: k1(:,:)
     real(r8), pointer :: k2(:,:)
     real(r8), pointer :: k3(:,:)
     real(r8), pointer :: k4(:,:)
  end type TstatusFlux
  !== Hongyi
 
  ! parameters to be calibrated. Ideally, these parameters are supposed to be uniform for one region
  public :: Tparameter
  type Tparameter
     real(r8), pointer :: c_nr(:)       ! coefficient to adjust the manning's roughness of channels
     real(r8), pointer :: c_nh(:)       ! coefficient to adjust the manning's roughness of overland flow across hillslopes
     real(r8), pointer :: c_twid(:)     ! coefficient to adjust the width of sub-reach channel
  end type Tparameter 

  !== Hongyi
  type (Tcontrol)    , public :: Tctl
  type (Tspatialunit), public :: TUnit
  type (TstatusFlux) , public :: TRunoff
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
             rtmCTL%outletg(begr:endr),           &
             rtmCTL%runofflnd_nt1(begr:endr),     &
             rtmCTL%runofflnd_nt2(begr:endr),     &
             rtmCTL%runoffocn_nt1(begr:endr),     &
             rtmCTL%runoffocn_nt2(begr:endr),     &
             rtmCTL%runofftot_nt1(begr:endr),     &
             rtmCTL%runofftot_nt2(begr:endr),     &
             rtmCTL%runoffdir_nt1(begr:endr),     &
             rtmCTL%runoffdir_nt2(begr:endr),     &
             rtmCTL%volr_nt1(begr:endr),          &
             rtmCTL%volr_nt2(begr:endr),          &
             rtmCTL%dvolrdtlnd_nt1(begr:endr),    &
             rtmCTL%dvolrdtlnd_nt2(begr:endr),    &
             rtmCTL%dvolrdtocn_nt1(begr:endr),    &
             rtmCTL%dvolrdtocn_nt2(begr:endr),    &
             rtmCTL%qsur_nt1(begr:endr),          &
             rtmCTL%qsur_nt2(begr:endr),          &
             rtmCTL%qsub_nt1(begr:endr),          &
             rtmCTL%qsub_nt2(begr:endr),          &
             rtmCTL%qgwl_nt1(begr:endr),          &
             rtmCTL%qgwl_nt2(begr:endr),          &
             rtmCTL%qdto_nt1(begr:endr),          &
             rtmCTL%qdto_nt2(begr:endr),          &
             rtmCTL%mask(begr:endr),              &
             rtmCTL%gindex(begr:endr),            &
             rtmCTL%fthresh(begr:endr),           &
             rtmCTL%flood(begr:endr),             &
             rtmCTL%direct(begr:endr,nt_rtm),     &
             rtmCTL%wh(begr:endr,nt_rtm),         &
             rtmCTL%wt(begr:endr,nt_rtm),         &
             rtmCTL%wr(begr:endr,nt_rtm),         &
             rtmCTL%erout(begr:endr,nt_rtm),      &
             rtmCTL%qsur(begr:endr,nt_rtm),       & 
             rtmCTL%qsub(begr:endr,nt_rtm),       &
             rtmCTL%qgwl(begr:endr,nt_rtm),       &
             rtmCTL%qdto(begr:endr,nt_rtm),       &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff local arrays'
       call shr_sys_abort
    end if

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

    rtmCTL%qsur(:,:)        = 0._r8
    rtmCTL%qsub(:,:)        = 0._r8
    rtmCTL%qgwl(:,:)        = 0._r8
    rtmCTL%qdto(:,:)        = 0._r8

  end subroutine RunoffInit

end module RunoffMod
