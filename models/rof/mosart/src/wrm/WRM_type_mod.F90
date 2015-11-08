!
MODULE WRM_type_mod
! Description: module for public data structure
! 
! Developed by Nathalie Voisin 01/31/2012
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod  , only :  r8 => shr_kind_r8, SHR_KIND_CL

  implicit none
  private
  
! control information for subbasin-based representation
  type WRMcontrol_subw
     integer :: NDam             ! number of dams
     integer :: localNumDam      ! number of dams on decomposition

     integer :: month            ! month of the simulation
     integer :: year             ! year of the simulation
     integer :: WRMFlag          ! Flag for  using the water resources management model or not
     integer :: ExtractionFlag   ! Flag for whether including Water Demand extraction : 1--> water extraction from each subbasin , uses extraction module ; 0 --> no extraction
     integer :: ExtractionMainChannelFlag   ! Flag for whether including Water Demand extraction from the main channel, only if unregulated flow ( bells and whistle here based on RegulationFlag is 1, etc)
     integer :: RegulationFlag   ! Flag whether to use reseervoir regulation or not : 1 --> use reservoir module; 0--> natural flow
     integer :: ReturnFlowFlag   ! Flag for return flow, this flag indicates withdrwals and consumptive uses
     integer :: TotalDemandFlag  ! Flag to indicate if the demand includes irrigation and non irrigation demands
     integer :: GroundWaterFlag  ! Flag to know if demand needs to be separated with GW-SW

     character(len=350) :: paraFile         ! the path of the parameter files
     character(len=350) :: paraPath         ! the path of the parameter files
     character(len=350) :: demandPath       ! the path of the water demand data
     character(len=350) :: grndwaterPath    ! the path for the groundwater demand share
     character(len=350) :: outPath          ! the path of the output file(s)
     character(len=350) :: damListFile      ! name of the file containing Dam list
     integer, pointer :: out_ID(:)          ! the indices of the outlet subbasins whether the stations are located
     character(len=80), pointer :: out_name(:)  ! the name of the outlets  
     character(len=80) :: curOutlet         ! the name of the current outlet
  end type WRMcontrol_subw
  
     ! --- Topographic and geometric properties, applicable for both grid- and subbasin-based representations
  type WRMspatialunit
     ! grid properties
     integer , pointer :: mask(:)           ! (nd) mask of a cell, value = reservoir GRanD ID, 0=reservoir (later, can add run-on-the-river)
     integer , pointer :: icell(:)          ! (nd) subbasin ID the Ndam
     integer , pointer :: INVicell(:)       ! (b:e) give the reservoir ID 1->NDAM located in icell
     integer , pointer :: INVisubw(:)       ! (??) need the equivilent for subw for the extraction module
     integer           :: NUnitID           ! max ID number of the units - needed by WRM
     integer , pointer :: isDam(:)          ! (b:e) Dam ID
     character(len=100), pointer :: DamName(:) ! (nd) dam name
     integer , pointer :: dam_depend(:,:)   ! (nd,b:e) dependence of each dam to subbasins - array of IDs
     integer , pointer :: dam_Ndepend(:)    ! (b:e) give the number of dependent subbasins to each reservoir
     integer , pointer :: subw_depend(:,:)  ! (b:e,nd) dependence of each subbasin to a certain number of dams, map IDs
     integer , pointer :: subw_Ndepend(:)   ! (b:e) number of reservoir from which the subbasin depends
     integer , pointer :: subw_damlist(:)   ! (1:Ndams) global dams needed on this pe
     real(r8), pointer :: Surfarea(:)       ! (nd) surface area of the reservoir
     real(r8), pointer :: InstCap(:)        ! (nd) instance energy capacity (MW)
     real(r8), pointer :: StorCap(:)        ! (nd) maximum storage capacity of the reservoir
     real(r8), pointer :: Height(:)         ! (nd) height of the reservoir
     real(r8), pointer :: Length(:)         ! (nd) Length of the reservoir
     real(r8), pointer :: Depth(:)          ! (nd) depth of the reservoir
     real(r8), pointer :: MeanMthFlow(:,:)  ! (nd,13) long term mean monthly flow
     real(r8), pointer :: INVc(:)           ! (nd) inverse of c of Biemans 2011 ands Hanasaki 2006 RUnoff/Capacity 
     real(r8), pointer :: TotStorCapDepend(:) ! (b:e) sum of the reservoir capacities each subw depends on
     real(r8), pointer :: TotInflowDepend(:)  ! (b:e) sum of the total inflow to dependen reservoir each subw depends on
     real(r8), pointer :: StorTarget(:,:)     ! (nd,13) monthly storage target to be enforced if non zero
     integer , pointer :: StorageCalibFlag(:) ! (nd) Flag that indicates if the storage targets are enforced or generic. Default is 0, generic
     real(r8), pointer :: MinStorTarget(:)  ! (nd)
     real(r8), pointer :: MaxStorTarget(:)  ! (nd)

     integer , pointer :: YEAR(:)           ! (nd) year dam was constructed and operationnal
     integer , pointer :: use_Irrig(:)      ! (nd) reservoir purpose irrigation 
     integer , pointer :: use_Elec(:)       ! (nd) hydropower
     integer , pointer :: use_FCon(:)       ! (nd) flood control
     integer , pointer :: use_Supp(:)       ! (nd) water supply
     integer , pointer :: use_Fish(:)       ! (nd) fish protection
     integer , pointer :: use_Navi(:)       ! (nd) use navigation
     integer , pointer :: use_Rec(:)        ! (nd) use recreation
     real(r8), pointer :: Withdrawal(:)     ! (nd) ratio of consumptive use over withdrawal - used to calibrate releases
     real(r8), pointer :: Conveyance(:)     ! (nd) ratio loss / consumptive use - transport efficiency
     integer , pointer :: MthStOp(:)        ! (nd) month showing the start of the operationnal year
     real(r8), pointer :: StorMthStOp(:)    ! (nd) storage at beginning of the start of the operationnal year
     real(r8), pointer :: MeanMthDemand(:,:)! (nd,13) longterm mean monthly demand
     ! flood control
     integer , pointer :: MthStFC(:)        ! (nd) month showing the strat of the flood control
     integer , pointer :: MthNdFC(:)        ! (nd) month showing the stop of the flood control
     integer , pointer :: MthFCtrack(:)     ! (nd) track if within FC or not                
  end type WRMspatialunit

  ! status and flux variables for liquid water
  type WRMwater
     real(r8), pointer :: TmpStoRelease(:,:,:) ! temporary store released water for 5 days in order to make it
!available for "storage release extraction" . If extraction from main stem then need to worry about environmental
!constraints based on points of extractions along the maoin stem. Here ensure the release is the environmental flow.

     real(r8), pointer :: supply(:)         ! (b:e) supply of surface water [m3]
     real(r8), pointer :: deficit(:)        ! (b:e) unmet demand [m3]
     real(r8), pointer :: demand(:)         ! (b:e) total withdrawl water demand [m3]
     real(r8), pointer :: WithDemIrrig(:)   ! (b:e) withdrawal demand for irrigation [m3]
     real(r8), pointer :: WithDemNonIrrig(:)! (b:e) withdrawal demand non irrigation [m3]
     real(r8), pointer :: ConDemIrrig(:)    ! (b:e) consumptive demand irrigation [m3]
     real(r8), pointer :: ConDemNonIrrig(:) ! (b:e) consumptive demand non irrigation [m3]
     real(r8), pointer :: SuppIrrig(:)      ! (b:e) withdrawal supply for irrigation [m3]
     real(r8), pointer :: SuppNonIrrig(:)   ! (b:e) withdrawal supply non irrigation [m3]
     real(r8), pointer :: ReturnIrrig(:)    ! (b:e) consumptive supply irrigation [m3]
     real(r8), pointer :: ReturnNonIrrig(:) ! (b:e) consumptive supply non irrigation [m3]
     real(r8), pointer :: GWShareNonIrrig(:)! (b:e) share of non irrigation demand assigned to gw
     real(r8), pointer :: GWShareIrrig(:)   ! (b:e) share of irrigation demand assigned to gw

     real(r8), pointer :: storage(:)        ! (nd) storage in the reservoir uhits
     real(r8), pointer :: pre_release(:,:)  ! (nd,13) pre-release without the interannual fluctutation
     real(r8), pointer :: release(:)        ! (nd) pre-release with the interannual fluctutation
     real(r8), pointer :: FCrelease (:)     ! (nd) Flood control release to get to storage target
     real(r8), pointer :: pot_evap(:)       ! (b:e) potential evaporation in mm from the grid
     real(r8), pointer :: Conveyance (:)    ! (nd) Conveyance loss flux
  end type WRMwater

  ! parameters to be calibrated. Ideally, these parameters are supposed to be uniform for one region
  !type WRMparameter
  !end type WRMparameter 

  type (WRMcontrol_subw), public :: ctlSubwWRM
  type (WRMspatialunit) , public :: WRMUnit
  type (WRMwater)       , public :: StorWater
  !type (WRMparameter)  , public :: WRMpara

  
end MODULE WRM_type_mod
