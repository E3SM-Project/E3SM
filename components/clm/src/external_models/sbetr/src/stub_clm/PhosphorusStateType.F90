module PhosphorusStateType
  use clm_varcon     , only : spval, ispval, c14ratio
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varpar     , only : ndecomp_pools, nlevdecomp_full
implicit none

  type, public :: phosphorusstate_type
    real(r8), pointer :: decomp_ppools_col            (:,:)  => null() ! col (gP/m2)  decomposing (litter, cwd, soil) P pools
    real(r8), pointer :: sminp_col                    (:)   => null()  ! col (gP/m2) soil mineral P
    real(r8), pointer :: solutionp_col                (:)  => null()       ! col (gP/m2) soil solution P
    real(r8), pointer :: labilep_col                  (:)   => null()      ! col (gP/m2) soil labile mineral P
    real(r8), pointer :: secondp_col                  (:)   => null()      ! col (gP/m2) soil secondary mineralP
    real(r8), pointer :: occlp_col                    (:)  => null()       ! col (gP/m2) soil occluded mineral P
    real(r8), pointer :: primp_col                    (:)  => null()       ! col (gP/m2) soil parimary mineral P
    real(r8), pointer :: totlitp_col                  (:)  => null()   ! col (gP/m2) total litter phosphorus
    real(r8), pointer :: totsomp_col                  (:) => null()    ! col (gP/m2) total soil organic matter phosphorus
    real(r8), pointer :: totecosysp_col               (:) => null()    ! col (gP/m2) total ecosystem phosphorus, incl veg
    real(r8), pointer :: totcolp_col                  (:) => null()    ! col (gP/m2) total column phosphorus, incl veg
    real(r8), pointer :: cwdp_col                     (:) => null()    ! col (gP/m2) Diagnostic: coarse woody debris P
    real(r8), pointer :: totlitp_1m_col               (:) => null()
    real(r8), pointer :: totsomp_1m_col               (:) => null()
    ! patch averaged to column variables
    real(r8), pointer :: totvegp_col                  (:)  => null()   ! col (gP/m2) total vegetation phosphorus (p2c)

    real(r8), pointer :: decomp_ppools_vr_col         (:,:,:) => null()    ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
    real(r8), pointer :: solutionp_vr_col             (:,:)  => null()     ! col (gP/m3) vertically-resolved soil solution P
    real(r8), pointer :: labilep_vr_col               (:,:) => null()      ! col (gP/m3) vertically-resolved soil labile mineral P
    real(r8), pointer :: secondp_vr_col               (:,:)  => null()     ! col (gP/m3) vertically-resolved soil secondary mineralP
    real(r8), pointer :: occlp_vr_col                 (:,:) => null()      ! col (gP/m3) vertically-resolved soil occluded mineral P
    real(r8), pointer :: primp_vr_col                 (:,:)  => null()     ! col (gP/m3) vertically-resolved soil parimary mineral P
    real(r8), pointer :: sminp_vr_col                 (:,:)  => null()     ! col (gP/m3) vertically-resolved soil total mineral P, diagnostic

  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type phosphorusstate_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(phosphorusstate_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate ( bounds )

    call this%InitCold ( bounds )

  end subroutine Init
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(phosphorusstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    allocate(this%decomp_ppools_col        (begc:endc,1:ndecomp_pools))   ; this%decomp_ppools_col        (:,:) = nan
    allocate(this%solutionp_col            (begc:endc))                   ; this%solutionp_col            (:)   = nan
    allocate(this%solutionp_col            (begc:endc))                   ; this%solutionp_col            (:)   = nan
    allocate(this%labilep_col              (begc:endc))                   ; this%labilep_col              (:)   = nan
    allocate(this%secondp_col              (begc:endc))                   ; this%secondp_col              (:)   = nan
    allocate(this%occlp_col                (begc:endc))                   ; this%occlp_col                (:)   = nan
    allocate(this%primp_col                (begc:endc))                   ; this%primp_col                (:)   = nan
    allocate(this%cwdp_col                 (begc:endc))                   ; this%cwdp_col                 (:)   = nan
    allocate(this%sminp_col                (begc:endc))                   ; this%sminp_col                (:)   = nan
    allocate(this%totecosysp_col           (begc:endc))                   ; this%totecosysp_col           (:)   = nan
    allocate(this%totcolp_col              (begc:endc))                   ; this%totcolp_col              (:)   = nan

    allocate(this%decomp_ppools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%decomp_ppools_vr_col(:,:,:)= nan
    allocate(this%solutionp_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr_col         (:,:) = nan
    allocate(this%labilep_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr_col           (:,:) = nan
    allocate(this%secondp_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr_col           (:,:) = nan
    allocate(this%occlp_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr_col             (:,:) = nan
    allocate(this%primp_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%primp_vr_col             (:,:) = nan
    allocate(this%sminp_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%sminp_vr_col             (:,:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use clm_varctl , only : nsrest, nsrStartup
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(phosphorusstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p,n,j,m            ! indices
    real(r8) ,pointer     :: gdp (:)                  ! global gdp data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: peatf (:)                ! global peatf data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: soilorder_rdin (:)       ! global soil order data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: abm (:)                  ! global abm data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: gti (:)                  ! read in - fmax (needs to be a pointer for use in ncdio)
    integer               :: dimid                    ! dimension id
    integer               :: ier                      ! error status
    type(file_desc_t)     :: ncid                     ! netcdf id
    logical               :: readvar
    character(len=256)    :: locfn                    ! local filename
    integer               :: begc, endc
    integer               :: begg, endg


  end subroutine initCold

end module PhosphorusStateType
