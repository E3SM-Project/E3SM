module zm_microphysics_state
   !----------------------------------------------------------------------------
   ! Purpose: microphysics state structure definition and methods for ZM
   ! Original Author: Xialiang Song and Guang Zhang, June 2010
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_params, only: r8
#else
   use shr_kind_mod,     only: r8=>shr_kind_r8
#endif

   public :: zm_microp_st           ! structure to hold state and tendency of ZM microphysics
   public :: zm_microp_st_alloc     ! allocate zm_microp_st variables
   public :: zm_microp_st_dealloc   ! deallocate zm_microp_st variables
   public :: zm_microp_st_ini       ! intialize zm_microp_st variables
   public :: zm_microp_st_zero      ! zero out zm_microp_st variables for a single column
   public :: zm_microp_st_scatter   ! scatter the gathered microphysic array data

!===================================================================================================

type :: zm_microp_st
   real(r8), allocatable, dimension(:,:) :: wu           ! vertical velocity
   real(r8), allocatable, dimension(:,:) :: qliq         ! convective cloud liquid water
   real(r8), allocatable, dimension(:,:) :: qice         ! convective cloud ice
   real(r8), allocatable, dimension(:,:) :: qrain        ! convective rain water
   real(r8), allocatable, dimension(:,:) :: qsnow        ! convective snow
   real(r8), allocatable, dimension(:,:) :: qgraupel     ! convective graupel
   real(r8), allocatable, dimension(:,:) :: qnl          ! convective cloud liquid water num concen.
   real(r8), allocatable, dimension(:,:) :: qni          ! convective cloud ice num concen.
   real(r8), allocatable, dimension(:,:) :: qnr          ! convective rain water num concen.
   real(r8), allocatable, dimension(:,:) :: qns          ! convective snow num concen.
   real(r8), allocatable, dimension(:,:) :: qng          ! convective graupel num concen.
   real(r8), allocatable, dimension(:)   :: rice         ! reserved ice (not yet in cldce) for energy integrals
   real(r8), allocatable, dimension(:,:) :: sprd         ! snow production rate
   real(r8), allocatable, dimension(:,:) :: mudpcu       ! width parameter of droplet size distr
   real(r8), allocatable, dimension(:,:) :: lambdadpcu   ! slope of cloud liquid size distr
   real(r8), allocatable, dimension(:,:) :: qcde         ! tmp for detrainment - cld liq mixing ratio       [kg/kg]
   real(r8), allocatable, dimension(:,:) :: qide         ! tmp for detrainment - cld ice mixing ratio       [kg/kg]
   real(r8), allocatable, dimension(:,:) :: qsde         ! tmp for detrainment - snow mixing ratio          [kg/kg]
   real(r8), allocatable, dimension(:,:) :: ncde         ! tmp for detrainment - cld liq number conc        [1/kg]
   real(r8), allocatable, dimension(:,:) :: nide         ! tmp for detrainment - cld ice number conc        [1/kg]
   real(r8), allocatable, dimension(:,:) :: nsde         ! tmp for detrainment - snow number conc           [1/kg]
   real(r8), allocatable, dimension(:,:) :: dif          ! detrainment of conv cld ice water mixing ratio
   real(r8), allocatable, dimension(:,:) :: dsf          ! detrainment of conv snow mixing ratio
   real(r8), allocatable, dimension(:,:) :: dnlf         ! detrainment of conv cld liq water num concen
   real(r8), allocatable, dimension(:,:) :: dnif         ! detrainment of conv cld ice num concen
   real(r8), allocatable, dimension(:,:) :: dnsf         ! detrainment of snow num concen
   real(r8), allocatable, dimension(:,:) :: frz          ! heating rate due to freezing
   real(r8), allocatable, dimension(:,:) :: autolm       ! mass tendency due to autoconversion of droplets to rain
   real(r8), allocatable, dimension(:,:) :: accrlm       ! mass tendency due to accretion of droplets by rain
   real(r8), allocatable, dimension(:,:) :: bergnm       ! mass tendency due to Bergeron process
   real(r8), allocatable, dimension(:,:) :: fhtimm       ! mass tendency due to immersion freezing
   real(r8), allocatable, dimension(:,:) :: fhtctm       ! mass tendency due to contact freezing
   real(r8), allocatable, dimension(:,:) :: fhmlm        ! mass tendency due to homogeneous freezing
   real(r8), allocatable, dimension(:,:) :: hmpim        ! mass tendency due to HM process
   real(r8), allocatable, dimension(:,:) :: accslm       ! mass tendency due to accretion of droplets by snow
   real(r8), allocatable, dimension(:,:) :: dlfm         ! mass tendency due to detrainment of droplet
   real(r8), allocatable, dimension(:,:) :: dsfm         ! mass tendency due to detrainment of snow
   real(r8), allocatable, dimension(:,:) :: autoln       ! num tendency due to autoconversion of droplets to rain
   real(r8), allocatable, dimension(:,:) :: accrln       ! num tendency due to accretion of droplets by rain
   real(r8), allocatable, dimension(:,:) :: bergnn       ! num tendency due to Bergeron process
   real(r8), allocatable, dimension(:,:) :: fhtimn       ! num tendency due to immersion freezing
   real(r8), allocatable, dimension(:,:) :: fhtctn       ! num tendency due to contact freezing
   real(r8), allocatable, dimension(:,:) :: fhmln        ! num tendency due to homogeneous freezing
   real(r8), allocatable, dimension(:,:) :: accsln       ! num tendency due to accretion of droplets by snow
   real(r8), allocatable, dimension(:,:) :: activn       ! num tendency due to droplets activation
   real(r8), allocatable, dimension(:,:) :: dlfn         ! num tendency due to detrainment of droplet
   real(r8), allocatable, dimension(:,:) :: dsfn         ! num tendency due to detrainment of snow
   real(r8), allocatable, dimension(:,:) :: autoim       ! mass tendency due to autoconversion of cloud ice to snow
   real(r8), allocatable, dimension(:,:) :: accsim       ! mass tendency due to accretion of cloud ice by snow
   real(r8), allocatable, dimension(:,:) :: difm         ! mass tendency due to detrainment of cloud ice
   real(r8), allocatable, dimension(:,:) :: nuclin       ! num tendency due to ice nucleation
   real(r8), allocatable, dimension(:,:) :: autoin       ! num tendency due to autoconversion of cloud ice to snow
   real(r8), allocatable, dimension(:,:) :: accsin       ! num tendency due to accretion of cloud ice by snow
   real(r8), allocatable, dimension(:,:) :: hmpin        ! num tendency due to HM process
   real(r8), allocatable, dimension(:,:) :: difn         ! num tendency due to detrainment of cloud ice
   real(r8), allocatable, dimension(:,:) :: cmel         ! mass tendency due to condensation
   real(r8), allocatable, dimension(:,:) :: cmei         ! mass tendency due to deposition
   real(r8), allocatable, dimension(:,:) :: trspcm       ! LWC tendency due to convective transport
   real(r8), allocatable, dimension(:,:) :: trspcn       ! droplet num tendency due to convective transport
   real(r8), allocatable, dimension(:,:) :: trspim       ! IWC tendency due to convective transport
   real(r8), allocatable, dimension(:,:) :: trspin       ! ice crystal num tendency due to convective transport
   real(r8), allocatable, dimension(:,:) :: accgrm       ! mass tendency due to collection of rain by graupel
   real(r8), allocatable, dimension(:,:) :: accglm       ! mass tendency due to collection of droplets by graupel
   real(r8), allocatable, dimension(:,:) :: accgslm      ! mass tendency of graupel due to collection of droplets by snow
   real(r8), allocatable, dimension(:,:) :: accgsrm      ! mass tendency of graupel due to collection of rain by snow
   real(r8), allocatable, dimension(:,:) :: accgirm      ! mass tendency of graupel due to collection of rain by ice
   real(r8), allocatable, dimension(:,:) :: accgrim      ! mass tendency of graupel due to collection of ice by rain
   real(r8), allocatable, dimension(:,:) :: accgrsm      ! mass tendency due to collection of snow by rain
   real(r8), allocatable, dimension(:,:) :: accgsln      ! num tendency of graupel due to collection of droplets by snow
   real(r8), allocatable, dimension(:,:) :: accgsrn      ! num tendency of graupel due to collection of rain by snow
   real(r8), allocatable, dimension(:,:) :: accgirn      ! num tendency of graupel due to collection of rain by ice
   real(r8), allocatable, dimension(:,:) :: accsrim      ! mass tendency of snow due to collection of ice by rain
   real(r8), allocatable, dimension(:,:) :: acciglm      ! mass tendency of ice mult(splintering) due to acc droplets by graupel
   real(r8), allocatable, dimension(:,:) :: accigrm      ! mass tendency of ice mult(splintering) due to acc rain by graupel
   real(r8), allocatable, dimension(:,:) :: accsirm      ! mass tendency of snow due to collection of rain by ice
   real(r8), allocatable, dimension(:,:) :: accigln      ! num tendency of ice mult(splintering) due to acc droplets by graupel
   real(r8), allocatable, dimension(:,:) :: accigrn      ! num tendency of ice mult(splintering) due to acc rain by graupel
   real(r8), allocatable, dimension(:,:) :: accsirn      ! num tendency of snow due to collection of rain by ice
   real(r8), allocatable, dimension(:,:) :: accgln       ! num tendency due to collection of droplets by graupel
   real(r8), allocatable, dimension(:,:) :: accgrn       ! num tendency due to collection of rain by graupel
   real(r8), allocatable, dimension(:,:) :: accilm       ! mass tendency of cloud ice due to collection of droplet by cloud ice
   real(r8), allocatable, dimension(:,:) :: acciln       ! number conc tendency of cloud ice due to collection of droplet by cloud ice
   real(r8), allocatable, dimension(:,:) :: fallrm       ! mass tendency of rain fallout
   real(r8), allocatable, dimension(:,:) :: fallsm       ! mass tendency of snow fallout
   real(r8), allocatable, dimension(:,:) :: fallgm       ! mass tendency of graupel fallout
   real(r8), allocatable, dimension(:,:) :: fallrn       ! num tendency of rain fallout
   real(r8), allocatable, dimension(:,:) :: fallsn       ! num tendency of snow fallout
   real(r8), allocatable, dimension(:,:) :: fallgn       ! num tendency of graupel fallout
   real(r8), allocatable, dimension(:,:) :: fhmrm        ! mass tendency due to homogeneous freezing of rain
end type zm_microp_st

!===================================================================================================
contains
!===================================================================================================

subroutine zm_microp_st_alloc(microp_st_in,ncol_in,nlev_in)
   !----------------------------------------------------------------------------
   ! Purpose: allocate zm_microp_st variables
   !----------------------------------------------------------------------------
   ! Arguments
   type(zm_microp_st), intent(inout) :: microp_st_in  ! state and tendency of convective microphysics
   integer,            intent(in   ) :: ncol_in       ! number of atmospheric columns to allocate
   integer,            intent(in   ) :: nlev_in       ! number of atmospheric levels to allocate
   !----------------------------------------------------------------------------
   allocate( & 
      microp_st_in%wu         (ncol_in,nlev_in), &
      microp_st_in%qliq       (ncol_in,nlev_in), &
      microp_st_in%qice       (ncol_in,nlev_in), &
      microp_st_in%qrain      (ncol_in,nlev_in), &
      microp_st_in%qsnow      (ncol_in,nlev_in), &
      microp_st_in%qgraupel   (ncol_in,nlev_in), &
      microp_st_in%qnl        (ncol_in,nlev_in), &
      microp_st_in%qni        (ncol_in,nlev_in), &
      microp_st_in%qnr        (ncol_in,nlev_in), &
      microp_st_in%qns        (ncol_in,nlev_in), &
      microp_st_in%qng        (ncol_in,nlev_in), &
      microp_st_in%rice       (ncol_in), &
      microp_st_in%sprd       (ncol_in,nlev_in), &
      microp_st_in%mudpcu     (ncol_in,nlev_in), &
      microp_st_in%lambdadpcu (ncol_in,nlev_in), &
      microp_st_in%qcde       (ncol_in,nlev_in), &
      microp_st_in%qide       (ncol_in,nlev_in), &
      microp_st_in%qsde       (ncol_in,nlev_in), &
      microp_st_in%ncde       (ncol_in,nlev_in), &
      microp_st_in%nide       (ncol_in,nlev_in), &
      microp_st_in%nsde       (ncol_in,nlev_in), &
      microp_st_in%dif        (ncol_in,nlev_in), &
      microp_st_in%dsf        (ncol_in,nlev_in), &
      microp_st_in%dnlf       (ncol_in,nlev_in), &
      microp_st_in%dnif       (ncol_in,nlev_in), &
      microp_st_in%dnsf       (ncol_in,nlev_in), &
      microp_st_in%frz        (ncol_in,nlev_in), &
      microp_st_in%autolm     (ncol_in,nlev_in), &
      microp_st_in%accrlm     (ncol_in,nlev_in), &
      microp_st_in%bergnm     (ncol_in,nlev_in), &
      microp_st_in%fhtimm     (ncol_in,nlev_in), &
      microp_st_in%fhtctm     (ncol_in,nlev_in), &
      microp_st_in%fhmlm      (ncol_in,nlev_in), &
      microp_st_in%hmpim      (ncol_in,nlev_in), &
      microp_st_in%accslm     (ncol_in,nlev_in), &
      microp_st_in%dlfm       (ncol_in,nlev_in), &
      microp_st_in%dsfm       (ncol_in,nlev_in), &
      microp_st_in%autoln     (ncol_in,nlev_in), &
      microp_st_in%accrln     (ncol_in,nlev_in), &
      microp_st_in%bergnn     (ncol_in,nlev_in), &
      microp_st_in%fhtimn     (ncol_in,nlev_in), &
      microp_st_in%fhtctn     (ncol_in,nlev_in), &
      microp_st_in%fhmln      (ncol_in,nlev_in), &
      microp_st_in%accsln     (ncol_in,nlev_in), &
      microp_st_in%activn     (ncol_in,nlev_in), &
      microp_st_in%dlfn       (ncol_in,nlev_in), &
      microp_st_in%dsfn       (ncol_in,nlev_in), &
      microp_st_in%autoim     (ncol_in,nlev_in), &
      microp_st_in%accsim     (ncol_in,nlev_in), &
      microp_st_in%difm       (ncol_in,nlev_in), &
      microp_st_in%nuclin     (ncol_in,nlev_in), &
      microp_st_in%autoin     (ncol_in,nlev_in), &
      microp_st_in%accsin     (ncol_in,nlev_in), &
      microp_st_in%hmpin      (ncol_in,nlev_in), &
      microp_st_in%difn       (ncol_in,nlev_in), &
      microp_st_in%cmel       (ncol_in,nlev_in), &
      microp_st_in%cmei       (ncol_in,nlev_in), &
      microp_st_in%trspcm     (ncol_in,nlev_in), &
      microp_st_in%trspcn     (ncol_in,nlev_in), &
      microp_st_in%trspim     (ncol_in,nlev_in), &
      microp_st_in%trspin     (ncol_in,nlev_in), &
      microp_st_in%accgrm     (ncol_in,nlev_in), &
      microp_st_in%accglm     (ncol_in,nlev_in), &
      microp_st_in%accgslm    (ncol_in,nlev_in), &
      microp_st_in%accgsrm    (ncol_in,nlev_in), &
      microp_st_in%accgirm    (ncol_in,nlev_in), &
      microp_st_in%accgrim    (ncol_in,nlev_in), &
      microp_st_in%accgrsm    (ncol_in,nlev_in), &
      microp_st_in%accgsln    (ncol_in,nlev_in), &
      microp_st_in%accgsrn    (ncol_in,nlev_in), &
      microp_st_in%accgirn    (ncol_in,nlev_in), &
      microp_st_in%accsrim    (ncol_in,nlev_in), &
      microp_st_in%acciglm    (ncol_in,nlev_in), &
      microp_st_in%accigrm    (ncol_in,nlev_in), &
      microp_st_in%accsirm    (ncol_in,nlev_in), &
      microp_st_in%accigln    (ncol_in,nlev_in), &
      microp_st_in%accigrn    (ncol_in,nlev_in), &
      microp_st_in%accsirn    (ncol_in,nlev_in), &
      microp_st_in%accgln     (ncol_in,nlev_in), &
      microp_st_in%accgrn     (ncol_in,nlev_in), &
      microp_st_in%accilm     (ncol_in,nlev_in), &
      microp_st_in%acciln     (ncol_in,nlev_in), &
      microp_st_in%fallrm     (ncol_in,nlev_in), &
      microp_st_in%fallsm     (ncol_in,nlev_in), &
      microp_st_in%fallgm     (ncol_in,nlev_in), &
      microp_st_in%fallrn     (ncol_in,nlev_in), &
      microp_st_in%fallsn     (ncol_in,nlev_in), &
      microp_st_in%fallgn     (ncol_in,nlev_in), &
      microp_st_in%fhmrm      (ncol_in,nlev_in) )

end subroutine zm_microp_st_alloc

!===================================================================================================

subroutine zm_microp_st_dealloc(microp_st_in)
   !----------------------------------------------------------------------------
   ! Purpose: deallocate zm_microp_st variables
   !----------------------------------------------------------------------------
   ! Arguments
   type(zm_microp_st) :: microp_st_in ! state and tendency of convective microphysics
   !----------------------------------------------------------------------------
   deallocate( &
      microp_st_in%wu,        &
      microp_st_in%qliq,      &
      microp_st_in%qice,      &
      microp_st_in%qrain,     &
      microp_st_in%qsnow,     &
      microp_st_in%qgraupel,  &
      microp_st_in%qnl,       &
      microp_st_in%qni,       &
      microp_st_in%qnr,       &
      microp_st_in%qns,       &
      microp_st_in%qng,       &
      microp_st_in%rice,      &
      microp_st_in%sprd,      &
      microp_st_in%mudpcu,    &
      microp_st_in%lambdadpcu,&
      microp_st_in%qcde,      &
      microp_st_in%qide,      &
      microp_st_in%qsde,      &
      microp_st_in%ncde,      &
      microp_st_in%nide,      &
      microp_st_in%nsde,      &
      microp_st_in%dif,       &
      microp_st_in%dsf,       &
      microp_st_in%dnlf,      &
      microp_st_in%dnif,      &
      microp_st_in%dnsf,      &
      microp_st_in%frz,       &
      microp_st_in%autolm,    &
      microp_st_in%accrlm,    &
      microp_st_in%bergnm,    &
      microp_st_in%fhtimm,    &
      microp_st_in%fhtctm,    &
      microp_st_in%fhmlm ,    &
      microp_st_in%hmpim ,    &
      microp_st_in%accslm,    &
      microp_st_in%dlfm  ,    &
      microp_st_in%dsfm  ,    &
      microp_st_in%autoln,    &
      microp_st_in%accrln,    &
      microp_st_in%bergnn,    &
      microp_st_in%fhtimn,    &
      microp_st_in%fhtctn,    &
      microp_st_in%fhmln ,    &
      microp_st_in%accsln,    &
      microp_st_in%activn,    &
      microp_st_in%dlfn  ,    &
      microp_st_in%dsfn  ,    &
      microp_st_in%autoim,    &
      microp_st_in%accsim,    &
      microp_st_in%difm  ,    &
      microp_st_in%nuclin,    &
      microp_st_in%autoin,    &
      microp_st_in%accsin,    &
      microp_st_in%hmpin,     &
      microp_st_in%difn,      &
      microp_st_in%cmel,      &
      microp_st_in%cmei,      &
      microp_st_in%trspcm,    &
      microp_st_in%trspcn,    &
      microp_st_in%trspim,    &
      microp_st_in%trspin,    &
      microp_st_in%accgrm,    &
      microp_st_in%accglm,    &
      microp_st_in%accgslm,   &
      microp_st_in%accgsrm,   &
      microp_st_in%accgirm,   &
      microp_st_in%accgrim,   &
      microp_st_in%accgrsm,   &
      microp_st_in%accgsln,   &
      microp_st_in%accgsrn,   &
      microp_st_in%accgirn,   &
      microp_st_in%accsrim,   &
      microp_st_in%acciglm,   &
      microp_st_in%accigrm,   &
      microp_st_in%accsirm,   &
      microp_st_in%accigln,   &
      microp_st_in%accigrn,   &
      microp_st_in%accsirn,   &
      microp_st_in%accgln,    &
      microp_st_in%accgrn,    &
      microp_st_in%accilm,    &
      microp_st_in%acciln,    &
      microp_st_in%fallrm,    &
      microp_st_in%fallsm,    &
      microp_st_in%fallgm,    &
      microp_st_in%fallrn,    &
      microp_st_in%fallsn,    &
      microp_st_in%fallgn,    &
      microp_st_in%fhmrm )

end subroutine zm_microp_st_dealloc

!===================================================================================================

subroutine zm_microp_st_zero(microp_st_in,icol_in,nlev_in)
   !----------------------------------------------------------------------------
   ! Purpose: zero out zm_microp_st variables for a single column
   !----------------------------------------------------------------------------
   ! Arguments
   type(zm_microp_st), intent(inout) :: microp_st_in ! state and tendency of convective microphysics
   integer,            intent(in   ) :: icol_in      ! atmospheric column index
   integer,            intent(in   ) :: nlev_in      ! number of atmospheric levels to initialize
   !----------------------------------------------------------------------------
   microp_st_in%wu        (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qliq      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qice      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qrain     (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qsnow     (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qgraupel  (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qnl       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qni       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qnr       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qns       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qng       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%rice      (icol_in) = 0._r8
   microp_st_in%sprd      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%mudpcu    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%lambdadpcu(icol_in,1:nlev_in) = 0._r8
   microp_st_in%qcde      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qide      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%qsde      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%ncde      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%nide      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%nsde      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dif       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dsf       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dnlf      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dnif      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dnsf      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%frz       (icol_in,1:nlev_in) = 0._r8
   microp_st_in%autolm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accrlm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%bergnm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fhtimm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fhtctm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fhmlm     (icol_in,1:nlev_in) = 0._r8
   microp_st_in%hmpim     (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accslm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dlfm      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dsfm      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%autoln    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accrln    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%bergnn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fhtimn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fhtctn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fhmln     (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accsln    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%activn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dlfn      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%dsfn      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%cmel      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%autoim    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accsim    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%difm      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%cmei      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%nuclin    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%autoin    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accsin    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%hmpin     (icol_in,1:nlev_in) = 0._r8
   microp_st_in%difn      (icol_in,1:nlev_in) = 0._r8
   microp_st_in%trspcm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%trspcn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%trspim    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%trspin    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgrm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accglm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgslm   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgsrm   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgirm   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgrim   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgrsm   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgsln   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgsrn   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgirn   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accsrim   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%acciglm   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accigrm   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accsirm   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accigln   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accigrn   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accsirn   (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgln    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accgrn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%accilm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%acciln    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fallrm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fallsm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fallgm    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fallrn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fallsn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fallgn    (icol_in,1:nlev_in) = 0._r8
   microp_st_in%fhmrm     (icol_in,1:nlev_in) = 0._r8
end subroutine zm_microp_st_zero

!===================================================================================================

subroutine zm_microp_st_ini(microp_st_in,ncol_in,nlev_in)
   !----------------------------------------------------------------------------
   ! Purpose: initialize zm_microp_st variables
   !----------------------------------------------------------------------------
   ! Arguments
   type(zm_microp_st), intent(inout) :: microp_st_in  ! state and tendency of convective microphysics
   integer,            intent(in   ) :: ncol_in       ! number of atmospheric columns to initialize
   integer,            intent(in   ) :: nlev_in       ! number of atmospheric levels to initialize
   !----------------------------------------------------------------------------
   integer :: i
   !----------------------------------------------------------------------------
   do i = 1,ncol_in
      call zm_microp_st_zero(microp_st_in,i,nlev_in)
   end do
end subroutine zm_microp_st_ini

!===================================================================================================

subroutine zm_microp_st_scatter(microp_st_gth,microp_st_out,pcols,lengath,nlev_in,ideep)
   !----------------------------------------------------------------------------
   ! Purpose: gather microphysic arrays from microp_st to microp_st_in
   !----------------------------------------------------------------------------
   ! Arguments
   type(zm_microp_st), intent(in   ) :: microp_st_gth ! input gathered state and tendency of convective microphysics
   type(zm_microp_st), intent(inout) :: microp_st_out ! output scattered state and tendency of convective microphysics
   integer,            intent(in   ) :: lengath       ! number of gathered columns
   integer,            intent(in   ) :: pcols         ! number of columns for ideep
   integer,            intent(in   ) :: nlev_in       ! number of atmospheric levels to initialize
   integer,            intent(in   ) :: ideep(pcols)  ! flag for active columns
   !----------------------------------------------------------------------------
   integer :: i,k
   !----------------------------------------------------------------------------
   do i = 1,lengath
      microp_st_out%rice      (ideep(i)) = microp_st_gth%rice      (i)
      do k = 1,nlev_in   
         microp_st_out%wu        (ideep(i),k) = microp_st_gth%wu        (i,k)
         microp_st_out%qliq      (ideep(i),k) = microp_st_gth%qliq      (i,k)
         microp_st_out%qice      (ideep(i),k) = microp_st_gth%qice      (i,k)
         microp_st_out%qrain     (ideep(i),k) = microp_st_gth%qrain     (i,k)
         microp_st_out%qsnow     (ideep(i),k) = microp_st_gth%qsnow     (i,k)
         microp_st_out%qgraupel  (ideep(i),k) = microp_st_gth%qgraupel  (i,k)
         microp_st_out%qnl       (ideep(i),k) = microp_st_gth%qnl       (i,k)
         microp_st_out%qni       (ideep(i),k) = microp_st_gth%qni       (i,k)
         microp_st_out%qnr       (ideep(i),k) = microp_st_gth%qnr       (i,k)
         microp_st_out%qns       (ideep(i),k) = microp_st_gth%qns       (i,k)
         microp_st_out%qng       (ideep(i),k) = microp_st_gth%qng       (i,k)
         microp_st_out%sprd      (ideep(i),k) = microp_st_gth%sprd      (i,k)
         microp_st_out%mudpcu    (ideep(i),k) = microp_st_gth%mudpcu    (i,k)
         microp_st_out%lambdadpcu(ideep(i),k) = microp_st_gth%lambdadpcu(i,k)
         microp_st_out%qcde      (ideep(i),k) = microp_st_gth%qcde      (i,k)
         microp_st_out%qide      (ideep(i),k) = microp_st_gth%qide      (i,k)
         microp_st_out%qsde      (ideep(i),k) = microp_st_gth%qsde      (i,k)
         microp_st_out%ncde      (ideep(i),k) = microp_st_gth%ncde      (i,k)
         microp_st_out%nide      (ideep(i),k) = microp_st_gth%nide      (i,k)
         microp_st_out%nsde      (ideep(i),k) = microp_st_gth%nsde      (i,k)
         microp_st_out%dif       (ideep(i),k) = microp_st_gth%dif       (i,k)
         microp_st_out%dsf       (ideep(i),k) = microp_st_gth%dsf       (i,k)
         microp_st_out%dnlf      (ideep(i),k) = microp_st_gth%dnlf      (i,k)
         microp_st_out%dnif      (ideep(i),k) = microp_st_gth%dnif      (i,k)
         microp_st_out%dnsf      (ideep(i),k) = microp_st_gth%dnsf      (i,k)
         microp_st_out%frz       (ideep(i),k) = microp_st_gth%frz       (i,k)
         microp_st_out%autolm    (ideep(i),k) = microp_st_gth%autolm    (i,k)
         microp_st_out%accrlm    (ideep(i),k) = microp_st_gth%accrlm    (i,k)
         microp_st_out%bergnm    (ideep(i),k) = microp_st_gth%bergnm    (i,k)
         microp_st_out%fhtimm    (ideep(i),k) = microp_st_gth%fhtimm    (i,k)
         microp_st_out%fhtctm    (ideep(i),k) = microp_st_gth%fhtctm    (i,k)
         microp_st_out%fhmlm     (ideep(i),k) = microp_st_gth%fhmlm     (i,k)
         microp_st_out%hmpim     (ideep(i),k) = microp_st_gth%hmpim     (i,k)
         microp_st_out%accslm    (ideep(i),k) = microp_st_gth%accslm    (i,k)
         microp_st_out%dlfm      (ideep(i),k) = microp_st_gth%dlfm      (i,k)
         microp_st_out%dsfm      (ideep(i),k) = microp_st_gth%dsfm      (i,k)
         microp_st_out%autoln    (ideep(i),k) = microp_st_gth%autoln    (i,k)
         microp_st_out%accrln    (ideep(i),k) = microp_st_gth%accrln    (i,k)
         microp_st_out%bergnn    (ideep(i),k) = microp_st_gth%bergnn    (i,k)
         microp_st_out%fhtimn    (ideep(i),k) = microp_st_gth%fhtimn    (i,k)
         microp_st_out%fhtctn    (ideep(i),k) = microp_st_gth%fhtctn    (i,k)
         microp_st_out%fhmln     (ideep(i),k) = microp_st_gth%fhmln     (i,k)
         microp_st_out%accsln    (ideep(i),k) = microp_st_gth%accsln    (i,k)
         microp_st_out%activn    (ideep(i),k) = microp_st_gth%activn    (i,k)
         microp_st_out%dlfn      (ideep(i),k) = microp_st_gth%dlfn      (i,k)
         microp_st_out%dsfn      (ideep(i),k) = microp_st_gth%dsfn      (i,k)
         microp_st_out%cmel      (ideep(i),k) = microp_st_gth%cmel      (i,k)
         microp_st_out%autoim    (ideep(i),k) = microp_st_gth%autoim    (i,k)
         microp_st_out%accsim    (ideep(i),k) = microp_st_gth%accsim    (i,k)
         microp_st_out%difm      (ideep(i),k) = microp_st_gth%difm      (i,k)
         microp_st_out%cmei      (ideep(i),k) = microp_st_gth%cmei      (i,k)
         microp_st_out%nuclin    (ideep(i),k) = microp_st_gth%nuclin    (i,k)
         microp_st_out%autoin    (ideep(i),k) = microp_st_gth%autoin    (i,k)
         microp_st_out%accsin    (ideep(i),k) = microp_st_gth%accsin    (i,k)
         microp_st_out%hmpin     (ideep(i),k) = microp_st_gth%hmpin     (i,k)
         microp_st_out%difn      (ideep(i),k) = microp_st_gth%difn      (i,k)
         microp_st_out%trspcm    (ideep(i),k) = microp_st_gth%trspcm    (i,k)
         microp_st_out%trspcn    (ideep(i),k) = microp_st_gth%trspcn    (i,k)
         microp_st_out%trspim    (ideep(i),k) = microp_st_gth%trspim    (i,k)
         microp_st_out%trspin    (ideep(i),k) = microp_st_gth%trspin    (i,k)
         microp_st_out%accgrm    (ideep(i),k) = microp_st_gth%accgrm    (i,k)
         microp_st_out%accglm    (ideep(i),k) = microp_st_gth%accglm    (i,k)
         microp_st_out%accgslm   (ideep(i),k) = microp_st_gth%accgslm   (i,k)
         microp_st_out%accgsrm   (ideep(i),k) = microp_st_gth%accgsrm   (i,k)
         microp_st_out%accgirm   (ideep(i),k) = microp_st_gth%accgirm   (i,k)
         microp_st_out%accgrim   (ideep(i),k) = microp_st_gth%accgrim   (i,k)
         microp_st_out%accgrsm   (ideep(i),k) = microp_st_gth%accgrsm   (i,k)
         microp_st_out%accgsln   (ideep(i),k) = microp_st_gth%accgsln   (i,k)
         microp_st_out%accgsrn   (ideep(i),k) = microp_st_gth%accgsrn   (i,k)
         microp_st_out%accgirn   (ideep(i),k) = microp_st_gth%accgirn   (i,k)
         microp_st_out%accsrim   (ideep(i),k) = microp_st_gth%accsrim   (i,k)
         microp_st_out%acciglm   (ideep(i),k) = microp_st_gth%acciglm   (i,k)
         microp_st_out%accigrm   (ideep(i),k) = microp_st_gth%accigrm   (i,k)
         microp_st_out%accsirm   (ideep(i),k) = microp_st_gth%accsirm   (i,k)
         microp_st_out%accigln   (ideep(i),k) = microp_st_gth%accigln   (i,k)
         microp_st_out%accigrn   (ideep(i),k) = microp_st_gth%accigrn   (i,k)
         microp_st_out%accsirn   (ideep(i),k) = microp_st_gth%accsirn   (i,k)
         microp_st_out%accgln    (ideep(i),k) = microp_st_gth%accgln    (i,k)
         microp_st_out%accgrn    (ideep(i),k) = microp_st_gth%accgrn    (i,k)
         microp_st_out%accilm    (ideep(i),k) = microp_st_gth%accilm    (i,k)
         microp_st_out%acciln    (ideep(i),k) = microp_st_gth%acciln    (i,k)
         microp_st_out%fallrm    (ideep(i),k) = microp_st_gth%fallrm    (i,k)
         microp_st_out%fallsm    (ideep(i),k) = microp_st_gth%fallsm    (i,k)
         microp_st_out%fallgm    (ideep(i),k) = microp_st_gth%fallgm    (i,k)
         microp_st_out%fallrn    (ideep(i),k) = microp_st_gth%fallrn    (i,k)
         microp_st_out%fallsn    (ideep(i),k) = microp_st_gth%fallsn    (i,k)
         microp_st_out%fallgn    (ideep(i),k) = microp_st_gth%fallgn    (i,k)
         microp_st_out%fhmrm     (ideep(i),k) = microp_st_gth%fhmrm     (i,k)
      end do
   end do
end subroutine zm_microp_st_scatter

!===================================================================================================

end module zm_microphysics_state
