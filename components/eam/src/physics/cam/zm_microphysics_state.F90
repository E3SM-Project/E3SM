module  zm_microphysics_state
  !-----------------------------------------------------------------------------
  ! Purpose: microphysics state structure definition and methods for ZM
  ! Original Author: Xialiang Song and Guang Zhang, June 2010
  !-----------------------------------------------------------------------------
  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only: pcols, pver, pverp

  public :: zm_microp_st ! structure to hold state and tendency of ZM microphysics
  public :: zm_microp_st_alloc
  public :: zm_microp_st_dealloc
  public :: zm_microp_st_ini
  public :: zm_microp_st_gb

!===================================================================================================

type :: zm_microp_st
  real(r8), allocatable, dimension(:,:) :: wu       ! vertical velocity
  real(r8), allocatable, dimension(:,:) :: qliq     ! convective cloud liquid water.
  real(r8), allocatable, dimension(:,:) :: qice     ! convective cloud ice.
  real(r8), allocatable, dimension(:,:) :: qrain    ! convective rain water.
  real(r8), allocatable, dimension(:,:) :: qsnow    ! convective snow.
  real(r8), allocatable, dimension(:,:) :: qgraupel ! convective graupel.
  real(r8), allocatable, dimension(:,:) :: qnl      ! convective cloud liquid water num concen.
  real(r8), allocatable, dimension(:,:) :: qni      ! convective cloud ice num concen.
  real(r8), allocatable, dimension(:,:) :: qnr      ! convective rain water num concen.
  real(r8), allocatable, dimension(:,:) :: qns      ! convective snow num concen.
  real(r8), allocatable, dimension(:,:) :: qng      ! convective graupel num concen.
  real(r8), allocatable, dimension(:,:) :: autolm   ! mass tendency due to autoconversion of droplets to rain
  real(r8), allocatable, dimension(:,:) :: accrlm   ! mass tendency due to accretion of droplets by rain
  real(r8), allocatable, dimension(:,:) :: bergnm   ! mass tendency due to Bergeron process
  real(r8), allocatable, dimension(:,:) :: fhtimm   ! mass tendency due to immersion freezing
  real(r8), allocatable, dimension(:,:) :: fhtctm   ! mass tendency due to contact freezing
  real(r8), allocatable, dimension(:,:) :: fhmlm    ! mass tendency due to homogeneous freezing
  real(r8), allocatable, dimension(:,:) :: hmpim    ! mass tendency due to HM process
  real(r8), allocatable, dimension(:,:) :: accslm   ! mass tendency due to accretion of droplets by snow
  real(r8), allocatable, dimension(:,:) :: dlfm     ! mass tendency due to detrainment of droplet
  real(r8), allocatable, dimension(:,:) :: autoln   ! num tendency due to autoconversion of droplets to rain
  real(r8), allocatable, dimension(:,:) :: accrln   ! num tendency due to accretion of droplets by rain
  real(r8), allocatable, dimension(:,:) :: bergnn   ! num tendency due to Bergeron process
  real(r8), allocatable, dimension(:,:) :: fhtimn   ! num tendency due to immersion freezing
  real(r8), allocatable, dimension(:,:) :: fhtctn   ! num tendency due to contact freezing
  real(r8), allocatable, dimension(:,:) :: fhmln    ! num tendency due to homogeneous freezing
  real(r8), allocatable, dimension(:,:) :: accsln   ! num tendency due to accretion of droplets by snow
  real(r8), allocatable, dimension(:,:) :: activn   ! num tendency due to droplets activation
  real(r8), allocatable, dimension(:,:) :: dlfn     ! num tendency due to detrainment of droplet
  real(r8), allocatable, dimension(:,:) :: autoim   ! mass tendency due to autoconversion of cloud ice to snow
  real(r8), allocatable, dimension(:,:) :: accsim   ! mass tendency due to accretion of cloud ice by snow
  real(r8), allocatable, dimension(:,:) :: difm     ! mass tendency due to detrainment of cloud ice
  real(r8), allocatable, dimension(:,:) :: nuclin   ! num tendency due to ice nucleation
  real(r8), allocatable, dimension(:,:) :: autoin   ! num tendency due to autoconversion of cloud ice to snow
  real(r8), allocatable, dimension(:,:) :: accsin   ! num tendency due to accretion of cloud ice by snow
  real(r8), allocatable, dimension(:,:) :: hmpin    ! num tendency due to HM process
  real(r8), allocatable, dimension(:,:) :: difn     ! num tendency due to detrainment of cloud ice
  real(r8), allocatable, dimension(:,:) :: cmel     ! mass tendency due to condensation
  real(r8), allocatable, dimension(:,:) :: cmei     ! mass tendency due to deposition
  real(r8), allocatable, dimension(:,:) :: trspcm   ! LWC tendency due to convective transport
  real(r8), allocatable, dimension(:,:) :: trspcn   ! droplet num tendency due to convective transport
  real(r8), allocatable, dimension(:,:) :: trspim   ! IWC tendency due to convective transport
  real(r8), allocatable, dimension(:,:) :: trspin   ! ice crystal num tendency due to convective transport
  real(r8), allocatable, dimension(:,:) :: accgrm   ! mass tendency due to collection of rain by graupel
  real(r8), allocatable, dimension(:,:) :: accglm   ! mass tendency due to collection of droplets by graupel
  real(r8), allocatable, dimension(:,:) :: accgslm  ! mass tendency of graupel due to collection of droplets by snow
  real(r8), allocatable, dimension(:,:) :: accgsrm  ! mass tendency of graupel due to collection of rain by snow
  real(r8), allocatable, dimension(:,:) :: accgirm  ! mass tendency of graupel due to collection of rain by ice
  real(r8), allocatable, dimension(:,:) :: accgrim  ! mass tendency of graupel due to collection of ice by rain
  real(r8), allocatable, dimension(:,:) :: accgrsm  ! mass tendency due to collection of snow by rain
  real(r8), allocatable, dimension(:,:) :: accgsln  ! num tendency of graupel due to collection of droplets by snow
  real(r8), allocatable, dimension(:,:) :: accgsrn  ! num tendency of graupel due to collection of rain by snow
  real(r8), allocatable, dimension(:,:) :: accgirn  ! num tendency of graupel due to collection of rain by ice
  real(r8), allocatable, dimension(:,:) :: accsrim  ! mass tendency of snow due to collection of ice by rain
  real(r8), allocatable, dimension(:,:) :: acciglm  ! mass tendency of ice mult(splintering) due to acc droplets by graupel
  real(r8), allocatable, dimension(:,:) :: accigrm  ! mass tendency of ice mult(splintering) due to acc rain by graupel
  real(r8), allocatable, dimension(:,:) :: accsirm  ! mass tendency of snow due to collection of rain by ice
  real(r8), allocatable, dimension(:,:) :: accigln  ! num tendency of ice mult(splintering) due to acc droplets by graupel
  real(r8), allocatable, dimension(:,:) :: accigrn  ! num tendency of ice mult(splintering) due to acc rain by graupel
  real(r8), allocatable, dimension(:,:) :: accsirn  ! num tendency of snow due to collection of rain by ice
  real(r8), allocatable, dimension(:,:) :: accgln   ! num tendency due to collection of droplets by graupel
  real(r8), allocatable, dimension(:,:) :: accgrn   ! num tendency due to collection of rain by graupel
  real(r8), allocatable, dimension(:,:) :: accilm   ! mass tendency of cloud ice due to collection of droplet by cloud ice
  real(r8), allocatable, dimension(:,:) :: acciln   ! number conc tendency of cloud ice due to collection of droplet by cloud ice
  real(r8), allocatable, dimension(:,:) :: fallrm   ! mass tendency of rain fallout
  real(r8), allocatable, dimension(:,:) :: fallsm   ! mass tendency of snow fallout
  real(r8), allocatable, dimension(:,:) :: fallgm   ! mass tendency of graupel fallout
  real(r8), allocatable, dimension(:,:) :: fallrn   ! num tendency of rain fallout
  real(r8), allocatable, dimension(:,:) :: fallsn   ! num tendency of snow fallout
  real(r8), allocatable, dimension(:,:) :: fallgn   ! num tendency of graupel fallout
  real(r8), allocatable, dimension(:,:) :: fhmrm    ! mass tendency due to homogeneous freezing of rain
end type zm_microp_st

!===================================================================================================
contains
!===================================================================================================

subroutine zm_microp_st_alloc(loc_microp_st)

  type(zm_microp_st)  :: loc_microp_st ! state and tendency of convective microphysics

  allocate( &
    loc_microp_st%wu       (pcols,pver), &
    loc_microp_st%qliq     (pcols,pver), &
    loc_microp_st%qice     (pcols,pver), &
    loc_microp_st%qrain    (pcols,pver), &
    loc_microp_st%qsnow    (pcols,pver), &
    loc_microp_st%qgraupel (pcols,pver), &
    loc_microp_st%qnl      (pcols,pver), &
    loc_microp_st%qni      (pcols,pver), &
    loc_microp_st%qnr      (pcols,pver), &
    loc_microp_st%qns      (pcols,pver), &
    loc_microp_st%qng      (pcols,pver), &
    loc_microp_st%autolm   (pcols,pver), &
    loc_microp_st%accrlm   (pcols,pver), &
    loc_microp_st%bergnm   (pcols,pver), &
    loc_microp_st%fhtimm   (pcols,pver), &
    loc_microp_st%fhtctm   (pcols,pver), &
    loc_microp_st%fhmlm    (pcols,pver), &
    loc_microp_st%hmpim    (pcols,pver), &
    loc_microp_st%accslm   (pcols,pver), &
    loc_microp_st%dlfm     (pcols,pver), &
    loc_microp_st%autoln   (pcols,pver), &
    loc_microp_st%accrln   (pcols,pver), &
    loc_microp_st%bergnn   (pcols,pver), &
    loc_microp_st%fhtimn   (pcols,pver), &
    loc_microp_st%fhtctn   (pcols,pver), &
    loc_microp_st%fhmln    (pcols,pver), &
    loc_microp_st%accsln   (pcols,pver), &
    loc_microp_st%activn   (pcols,pver), &
    loc_microp_st%dlfn     (pcols,pver), &
    loc_microp_st%autoim   (pcols,pver), &
    loc_microp_st%accsim   (pcols,pver), &
    loc_microp_st%difm     (pcols,pver), &
    loc_microp_st%nuclin   (pcols,pver), &
    loc_microp_st%autoin   (pcols,pver), &
    loc_microp_st%accsin   (pcols,pver), &
    loc_microp_st%hmpin    (pcols,pver), &
    loc_microp_st%difn     (pcols,pver), &
    loc_microp_st%cmel     (pcols,pver), &
    loc_microp_st%cmei     (pcols,pver), &
    loc_microp_st%trspcm   (pcols,pver), &
    loc_microp_st%trspcn   (pcols,pver), &
    loc_microp_st%trspim   (pcols,pver), &
    loc_microp_st%trspin   (pcols,pver), &
    loc_microp_st%accgrm   (pcols,pver), &
    loc_microp_st%accglm   (pcols,pver), &
    loc_microp_st%accgslm  (pcols,pver), &
    loc_microp_st%accgsrm  (pcols,pver), &
    loc_microp_st%accgirm  (pcols,pver), &
    loc_microp_st%accgrim  (pcols,pver), &
    loc_microp_st%accgrsm  (pcols,pver), &
    loc_microp_st%accgsln  (pcols,pver), &
    loc_microp_st%accgsrn  (pcols,pver), &
    loc_microp_st%accgirn  (pcols,pver), &
    loc_microp_st%accsrim  (pcols,pver), &
    loc_microp_st%acciglm  (pcols,pver), &
    loc_microp_st%accigrm  (pcols,pver), &
    loc_microp_st%accsirm  (pcols,pver), &
    loc_microp_st%accigln  (pcols,pver), &
    loc_microp_st%accigrn  (pcols,pver), &
    loc_microp_st%accsirn  (pcols,pver), &
    loc_microp_st%accgln   (pcols,pver), &
    loc_microp_st%accgrn   (pcols,pver), &
    loc_microp_st%accilm   (pcols,pver), &
    loc_microp_st%acciln   (pcols,pver), &
    loc_microp_st%fallrm   (pcols,pver), &
    loc_microp_st%fallsm   (pcols,pver), &
    loc_microp_st%fallgm   (pcols,pver), &
    loc_microp_st%fallrn   (pcols,pver), &
    loc_microp_st%fallsn   (pcols,pver), &
    loc_microp_st%fallgn   (pcols,pver), &
    loc_microp_st%fhmrm    (pcols,pver) )

end subroutine zm_microp_st_alloc

!===================================================================================================

subroutine zm_microp_st_dealloc(loc_microp_st)

  type(zm_microp_st)  :: loc_microp_st ! state and tendency of convective microphysics

  deallocate( &
    loc_microp_st%wu,       &
    loc_microp_st%qliq,     &
    loc_microp_st%qice,     &
    loc_microp_st%qrain,    &
    loc_microp_st%qsnow,    &
    loc_microp_st%qgraupel, &
    loc_microp_st%qnl,      &
    loc_microp_st%qni,      &
    loc_microp_st%qnr,      &
    loc_microp_st%qns,      &
    loc_microp_st%qng,      &
    loc_microp_st%autolm,   &
    loc_microp_st%accrlm,   &
    loc_microp_st%bergnm,   &
    loc_microp_st%fhtimm,   &
    loc_microp_st%fhtctm,   &
    loc_microp_st%fhmlm ,   &
    loc_microp_st%hmpim ,   &
    loc_microp_st%accslm,   &
    loc_microp_st%dlfm  ,   &
    loc_microp_st%autoln,   &
    loc_microp_st%accrln,   &
    loc_microp_st%bergnn,   &
    loc_microp_st%fhtimn,   &
    loc_microp_st%fhtctn,   &
    loc_microp_st%fhmln ,   &
    loc_microp_st%accsln,   &
    loc_microp_st%activn,   &
    loc_microp_st%dlfn  ,   &
    loc_microp_st%autoim,   &
    loc_microp_st%accsim,   &
    loc_microp_st%difm  ,   &
    loc_microp_st%nuclin,   &
    loc_microp_st%autoin,   &
    loc_microp_st%accsin,   &
    loc_microp_st%hmpin,    &
    loc_microp_st%difn,     &
    loc_microp_st%cmel,     &
    loc_microp_st%cmei,     &
    loc_microp_st%trspcm,   &
    loc_microp_st%trspcn,   &
    loc_microp_st%trspim,   &
    loc_microp_st%trspin,   &
    loc_microp_st%accgrm,   &
    loc_microp_st%accglm,   &
    loc_microp_st%accgslm,  &
    loc_microp_st%accgsrm,  &
    loc_microp_st%accgirm,  &
    loc_microp_st%accgrim,  &
    loc_microp_st%accgrsm,  &
    loc_microp_st%accgsln,  &
    loc_microp_st%accgsrn,  &
    loc_microp_st%accgirn,  &
    loc_microp_st%accsrim,  &
    loc_microp_st%acciglm,  &
    loc_microp_st%accigrm,  &
    loc_microp_st%accsirm,  &
    loc_microp_st%accigln,  &
    loc_microp_st%accigrn,  &
    loc_microp_st%accsirn,  &
    loc_microp_st%accgln,   &
    loc_microp_st%accgrn,   &
    loc_microp_st%accilm,   &
    loc_microp_st%acciln,   &
    loc_microp_st%fallrm,   &
    loc_microp_st%fallsm,   &
    loc_microp_st%fallgm,   &
    loc_microp_st%fallrn,   &
    loc_microp_st%fallsn,   &
    loc_microp_st%fallgn,   &
    loc_microp_st%fhmrm )

end subroutine zm_microp_st_dealloc

!===================================================================================================

subroutine zm_microp_st_ini(microp_st,ncol)

  type(zm_microp_st)  :: microp_st     ! state and tendency of convective microphysics
  integer, intent(in) :: ncol          ! number of atmospheric columns

  microp_st%wu       (1:ncol,1:pver) = 0._r8
  microp_st%qliq     (1:ncol,1:pver) = 0._r8
  microp_st%qice     (1:ncol,1:pver) = 0._r8
  microp_st%qrain    (1:ncol,1:pver) = 0._r8
  microp_st%qsnow    (1:ncol,1:pver) = 0._r8
  microp_st%qgraupel (1:ncol,1:pver) = 0._r8
  microp_st%qnl      (1:ncol,1:pver) = 0._r8
  microp_st%qni      (1:ncol,1:pver) = 0._r8
  microp_st%qnr      (1:ncol,1:pver) = 0._r8
  microp_st%qns      (1:ncol,1:pver) = 0._r8
  microp_st%qng      (1:ncol,1:pver) = 0._r8
  microp_st%autolm   (1:ncol,1:pver) = 0._r8
  microp_st%accrlm   (1:ncol,1:pver) = 0._r8
  microp_st%bergnm   (1:ncol,1:pver) = 0._r8
  microp_st%fhtimm   (1:ncol,1:pver) = 0._r8
  microp_st%fhtctm   (1:ncol,1:pver) = 0._r8
  microp_st%fhmlm    (1:ncol,1:pver) = 0._r8
  microp_st%hmpim    (1:ncol,1:pver) = 0._r8
  microp_st%accslm   (1:ncol,1:pver) = 0._r8
  microp_st%dlfm     (1:ncol,1:pver) = 0._r8
  microp_st%autoln   (1:ncol,1:pver) = 0._r8
  microp_st%accrln   (1:ncol,1:pver) = 0._r8
  microp_st%bergnn   (1:ncol,1:pver) = 0._r8
  microp_st%fhtimn   (1:ncol,1:pver) = 0._r8
  microp_st%fhtctn   (1:ncol,1:pver) = 0._r8
  microp_st%fhmln    (1:ncol,1:pver) = 0._r8
  microp_st%accsln   (1:ncol,1:pver) = 0._r8
  microp_st%activn   (1:ncol,1:pver) = 0._r8
  microp_st%dlfn     (1:ncol,1:pver) = 0._r8
  microp_st%cmel     (1:ncol,1:pver) = 0._r8
  microp_st%autoim   (1:ncol,1:pver) = 0._r8
  microp_st%accsim   (1:ncol,1:pver) = 0._r8
  microp_st%difm     (1:ncol,1:pver) = 0._r8
  microp_st%cmei     (1:ncol,1:pver) = 0._r8
  microp_st%nuclin   (1:ncol,1:pver) = 0._r8
  microp_st%autoin   (1:ncol,1:pver) = 0._r8
  microp_st%accsin   (1:ncol,1:pver) = 0._r8
  microp_st%hmpin    (1:ncol,1:pver) = 0._r8
  microp_st%difn     (1:ncol,1:pver) = 0._r8
  microp_st%trspcm   (1:ncol,1:pver) = 0._r8
  microp_st%trspcn   (1:ncol,1:pver) = 0._r8
  microp_st%trspim   (1:ncol,1:pver) = 0._r8
  microp_st%trspin   (1:ncol,1:pver) = 0._r8
  microp_st%accgrm   (1:ncol,1:pver) = 0._r8
  microp_st%accglm   (1:ncol,1:pver) = 0._r8
  microp_st%accgslm  (1:ncol,1:pver) = 0._r8
  microp_st%accgsrm  (1:ncol,1:pver) = 0._r8
  microp_st%accgirm  (1:ncol,1:pver) = 0._r8
  microp_st%accgrim  (1:ncol,1:pver) = 0._r8
  microp_st%accgrsm  (1:ncol,1:pver) = 0._r8
  microp_st%accgsln  (1:ncol,1:pver) = 0._r8
  microp_st%accgsrn  (1:ncol,1:pver) = 0._r8
  microp_st%accgirn  (1:ncol,1:pver) = 0._r8
  microp_st%accsrim  (1:ncol,1:pver) = 0._r8
  microp_st%acciglm  (1:ncol,1:pver) = 0._r8
  microp_st%accigrm  (1:ncol,1:pver) = 0._r8
  microp_st%accsirm  (1:ncol,1:pver) = 0._r8
  microp_st%accigln  (1:ncol,1:pver) = 0._r8
  microp_st%accigrn  (1:ncol,1:pver) = 0._r8
  microp_st%accsirn  (1:ncol,1:pver) = 0._r8
  microp_st%accgln   (1:ncol,1:pver) = 0._r8
  microp_st%accgrn   (1:ncol,1:pver) = 0._r8 
  microp_st%accilm   (1:ncol,1:pver) = 0._r8
  microp_st%acciln   (1:ncol,1:pver) = 0._r8
  microp_st%fallrm   (1:ncol,1:pver) = 0._r8
  microp_st%fallsm   (1:ncol,1:pver) = 0._r8
  microp_st%fallgm   (1:ncol,1:pver) = 0._r8
  microp_st%fallrn   (1:ncol,1:pver) = 0._r8
  microp_st%fallsn   (1:ncol,1:pver) = 0._r8
  microp_st%fallgn   (1:ncol,1:pver) = 0._r8
  microp_st%fhmrm    (1:ncol,1:pver) = 0._r8

end subroutine zm_microp_st_ini

!===================================================================================================

subroutine zm_microp_st_gb(microp_st,loc_microp_st,ideep,lengath)
  !-----------------------------------------------------------------------------
  ! Purpose: Gather microphysic arrays from microp_st to loc_microp_st
  !-----------------------------------------------------------------------------
  type(zm_microp_st)  :: microp_st     ! state and tendency of convective microphysics
  type(zm_microp_st)  :: loc_microp_st ! state and tendency of convective microphysics
  integer ideep(pcols)                 ! holds position of gathered points vs longitude index.
  integer lengath
  integer i,k

  do k = 1,pver
    do i = 1,lengath
      microp_st%wu       (ideep(i),k) = loc_microp_st%wu       (i,k)
      microp_st%qliq     (ideep(i),k) = loc_microp_st%qliq     (i,k)
      microp_st%qice     (ideep(i),k) = loc_microp_st%qice     (i,k)
      microp_st%qrain    (ideep(i),k) = loc_microp_st%qrain    (i,k)
      microp_st%qsnow    (ideep(i),k) = loc_microp_st%qsnow    (i,k)
      microp_st%qgraupel (ideep(i),k) = loc_microp_st%qgraupel (i,k)
      microp_st%qnl      (ideep(i),k) = loc_microp_st%qnl      (i,k)
      microp_st%qni      (ideep(i),k) = loc_microp_st%qni      (i,k)
      microp_st%qnr      (ideep(i),k) = loc_microp_st%qnr      (i,k)
      microp_st%qns      (ideep(i),k) = loc_microp_st%qns      (i,k)
      microp_st%qng      (ideep(i),k) = loc_microp_st%qng      (i,k)
      microp_st%autolm   (ideep(i),k) = loc_microp_st%autolm   (i,k)
      microp_st%accrlm   (ideep(i),k) = loc_microp_st%accrlm   (i,k)
      microp_st%bergnm   (ideep(i),k) = loc_microp_st%bergnm   (i,k)
      microp_st%fhtimm   (ideep(i),k) = loc_microp_st%fhtimm   (i,k)
      microp_st%fhtctm   (ideep(i),k) = loc_microp_st%fhtctm   (i,k)
      microp_st%fhmlm    (ideep(i),k) = loc_microp_st%fhmlm    (i,k)
      microp_st%hmpim    (ideep(i),k) = loc_microp_st%hmpim    (i,k)
      microp_st%accslm   (ideep(i),k) = loc_microp_st%accslm   (i,k)
      microp_st%dlfm     (ideep(i),k) = loc_microp_st%dlfm     (i,k)
      microp_st%autoln   (ideep(i),k) = loc_microp_st%autoln   (i,k)
      microp_st%accrln   (ideep(i),k) = loc_microp_st%accrln   (i,k)
      microp_st%bergnn   (ideep(i),k) = loc_microp_st%bergnn   (i,k)
      microp_st%fhtimn   (ideep(i),k) = loc_microp_st%fhtimn   (i,k)
      microp_st%fhtctn   (ideep(i),k) = loc_microp_st%fhtctn   (i,k)
      microp_st%fhmln    (ideep(i),k) = loc_microp_st%fhmln    (i,k)
      microp_st%accsln   (ideep(i),k) = loc_microp_st%accsln   (i,k)
      microp_st%activn   (ideep(i),k) = loc_microp_st%activn   (i,k)
      microp_st%dlfn     (ideep(i),k) = loc_microp_st%dlfn     (i,k)
      microp_st%cmel     (ideep(i),k) = loc_microp_st%cmel     (i,k)
      microp_st%autoim   (ideep(i),k) = loc_microp_st%autoim   (i,k)
      microp_st%accsim   (ideep(i),k) = loc_microp_st%accsim   (i,k)
      microp_st%difm     (ideep(i),k) = loc_microp_st%difm     (i,k)
      microp_st%cmei     (ideep(i),k) = loc_microp_st%cmei     (i,k)
      microp_st%nuclin   (ideep(i),k) = loc_microp_st%nuclin   (i,k)
      microp_st%autoin   (ideep(i),k) = loc_microp_st%autoin   (i,k)
      microp_st%accsin   (ideep(i),k) = loc_microp_st%accsin   (i,k)
      microp_st%hmpin    (ideep(i),k) = loc_microp_st%hmpin    (i,k)
      microp_st%difn     (ideep(i),k) = loc_microp_st%difn     (i,k)
      microp_st%trspcm   (ideep(i),k) = loc_microp_st%trspcm   (i,k)
      microp_st%trspcn   (ideep(i),k) = loc_microp_st%trspcn   (i,k)
      microp_st%trspim   (ideep(i),k) = loc_microp_st%trspim   (i,k)
      microp_st%trspin   (ideep(i),k) = loc_microp_st%trspin   (i,k)
      microp_st%accgrm   (ideep(i),k) = loc_microp_st%accgrm   (i,k)
      microp_st%accglm   (ideep(i),k) = loc_microp_st%accglm   (i,k)
      microp_st%accgslm  (ideep(i),k) = loc_microp_st%accgslm  (i,k)
      microp_st%accgsrm  (ideep(i),k) = loc_microp_st%accgsrm  (i,k)
      microp_st%accgirm  (ideep(i),k) = loc_microp_st%accgirm  (i,k)
      microp_st%accgrim  (ideep(i),k) = loc_microp_st%accgrim  (i,k)
      microp_st%accgrsm  (ideep(i),k) = loc_microp_st%accgrsm  (i,k)
      microp_st%accgsln  (ideep(i),k) = loc_microp_st%accgsln  (i,k)
      microp_st%accgsrn  (ideep(i),k) = loc_microp_st%accgsrn  (i,k)
      microp_st%accgirn  (ideep(i),k) = loc_microp_st%accgirn  (i,k)
      microp_st%accsrim  (ideep(i),k) = loc_microp_st%accsrim  (i,k)
      microp_st%acciglm  (ideep(i),k) = loc_microp_st%acciglm  (i,k)
      microp_st%accigrm  (ideep(i),k) = loc_microp_st%accigrm  (i,k)
      microp_st%accsirm  (ideep(i),k) = loc_microp_st%accsirm  (i,k)
      microp_st%accigln  (ideep(i),k) = loc_microp_st%accigln  (i,k)
      microp_st%accigrn  (ideep(i),k) = loc_microp_st%accigrn  (i,k)
      microp_st%accsirn  (ideep(i),k) = loc_microp_st%accsirn  (i,k)
      microp_st%accgln   (ideep(i),k) = loc_microp_st%accgln   (i,k)
      microp_st%accgrn   (ideep(i),k) = loc_microp_st%accgrn   (i,k)
      microp_st%accilm   (ideep(i),k) = loc_microp_st%accilm   (i,k)
      microp_st%acciln   (ideep(i),k) = loc_microp_st%acciln   (i,k)
      microp_st%fallrm   (ideep(i),k) = loc_microp_st%fallrm   (i,k)
      microp_st%fallsm   (ideep(i),k) = loc_microp_st%fallsm   (i,k)
      microp_st%fallgm   (ideep(i),k) = loc_microp_st%fallgm   (i,k)
      microp_st%fallrn   (ideep(i),k) = loc_microp_st%fallrn   (i,k)
      microp_st%fallsn   (ideep(i),k) = loc_microp_st%fallsn   (i,k)
      microp_st%fallgn   (ideep(i),k) = loc_microp_st%fallgn   (i,k)
      microp_st%fhmrm    (ideep(i),k) = loc_microp_st%fhmrm    (i,k)
    end do
  end do
end subroutine zm_microp_st_gb

!===================================================================================================

end module zm_microphysics_state
