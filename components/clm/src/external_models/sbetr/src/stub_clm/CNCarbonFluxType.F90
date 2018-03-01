module CNCarbonFluxType
  use clm_varcon             , only : spval, ispval, c14ratio
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use clm_varcon             , only : spval
  use clm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: carbonflux_type
    real(r8), pointer :: rr_col                                    (:)   => null()  ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: rr_patch                                  (:)  => null()   ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
    real(r8), pointer :: annsum_npp_patch                          (:) => null() ! patch annual sum of NPP (gC/m2/yr)
    real(r8), pointer :: agnpp_patch                               (:)  => null()   ! (gC/m2/s) aboveground NPP
    real(r8), pointer :: bgnpp_patch                               (:)  => null()   ! (gC/m2/s) belowground NPP
    real(r8), pointer :: hr_col                                    (:)=> null()
    real(r8), pointer :: hr_vr_col                                 (:,:) => null()
    real(r8), pointer :: phenology_c_to_litr_met_c_col             (:,:)=> null()
    real(r8), pointer :: phenology_c_to_litr_cel_c_col             (:,:)=> null()
    real(r8), pointer :: phenology_c_to_litr_lig_c_col             (:,:)=> null()
    real(r8), pointer :: dwt_livecrootc_to_cwdc_col                (:,:)=> null()
    real(r8), pointer :: m_decomp_cpools_to_fire_vr_col            (:,:,:)=> null()
    real(r8), pointer :: dwt_deadcrootc_to_cwdc_col                (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_lig_c_col              (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_cel_c_col              (:,:)=> null()
    real(r8), pointer :: dwt_frootc_to_litr_met_c_col              (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_met_c_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_cel_c_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_litr_lig_c_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_c_to_cwdc_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_met_c_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_cel_c_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_litr_lig_c_col               (:,:)=> null()
    real(r8), pointer :: harvest_c_to_cwdc_col                     (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_met_fire_col                  (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_cel_fire_col                  (:,:)=> null()
    real(r8), pointer :: m_c_to_litr_lig_fire_col                  (:,:)=> null()
    real(r8), pointer :: fire_mortality_c_to_cwdc_col              (:,:)=> null()
    real(r8), pointer :: fire_decomp_closs_col                     (:) => null()
    real(r8), pointer :: som_c_leached_col                         (:) => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type carbonflux_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(carbonflux_type) :: this
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
    class(carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%rr_col                  (begc:endc))                  ; this%rr_col                    (:)  =nan
    allocate(this%rr_patch                (begp:endp))                  ; this%rr_patch                    (:)  =nan

    allocate(this%annsum_npp_patch      (begp:endp)) ; this%annsum_npp_patch      (:) = nan
    allocate(this%agnpp_patch                       (begp:endp)) ; this%agnpp_patch                               (:) = nan
    allocate(this%bgnpp_patch                       (begp:endp)) ; this%bgnpp_patch                               (:) = nan
    allocate(this%hr_col (begc:endc)); this%hr_col(:) = nan

    allocate(this%phenology_c_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_met_c_col(:,:) = nan
    allocate(this%phenology_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_cel_c_col(:,:) = nan
    allocate(this%dwt_livecrootc_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%dwt_livecrootc_to_cwdc_col(:,:) = nan
    allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:7)); this%m_decomp_cpools_to_fire_vr_col(:,:,:)=nan
    allocate(this%dwt_deadcrootc_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootc_to_cwdc_col(:,:) = nan
    allocate(this%dwt_frootc_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_lig_c_col(:,:) = nan
    allocate(this%dwt_frootc_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_cel_c_col(:,:) = nan
    allocate(this%dwt_frootc_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_met_c_col(:,:) = nan
    allocate(this%gap_mortality_c_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_met_c_col(:,:) = nan
    allocate(this%gap_mortality_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_cel_c_col(:,:) = nan
    allocate(this%gap_mortality_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_lig_c_col(:,:) = nan
    allocate(this%gap_mortality_c_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_cwdc_col(:,:) = nan
    allocate(this%harvest_c_to_litr_met_c_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_met_c_col(:,:) = nan
    allocate(this%harvest_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_cel_c_col(:,:) = nan
    allocate(this%harvest_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_lig_c_col(:,:) = nan
    allocate(this%harvest_c_to_cwdc_col(begc:endc,1:nlevdecomp_full)); this%harvest_c_to_cwdc_col(:,:) = nan
    allocate(this%m_c_to_litr_met_fire_col(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_met_fire_col(:,:) = nan
    allocate(this%m_c_to_litr_cel_fire_col(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_cel_fire_col(:,:) = nan
    allocate(this%m_c_to_litr_lig_fire_col(begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_lig_fire_col(:,:) = nan
    allocate(this%fire_mortality_c_to_cwdc_col (begc:endc,1:nlevdecomp_full)); this%fire_mortality_c_to_cwdc_col (:,:) = nan
    allocate(this%phenology_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_lig_c_col(:,:) = nan
    allocate(this%som_c_leached_col(begc:endc)); this%som_c_leached_col(:) = nan

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
    class(carbonflux_type) :: this
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


    this%phenology_c_to_litr_met_c_col(:,:) = 0._r8
    this%phenology_c_to_litr_cel_c_col(:,:) = 0._r8
    this%dwt_livecrootc_to_cwdc_col(:,:) = 0._r8
    this%m_decomp_cpools_to_fire_vr_col(:,:,:)= 0._r8
    this%dwt_deadcrootc_to_cwdc_col(:,:) = 0._r8
    this%dwt_frootc_to_litr_lig_c_col(:,:) = 0._r8
    this%dwt_frootc_to_litr_cel_c_col(:,:) = 0._r8
    this%dwt_frootc_to_litr_met_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_litr_met_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_litr_cel_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_litr_lig_c_col(:,:) = 0._r8
    this%gap_mortality_c_to_cwdc_col(:,:) = 0._r8
    this%harvest_c_to_litr_met_c_col(:,:) = 0._r8
    this%harvest_c_to_litr_cel_c_col(:,:) = 0._r8
    this%harvest_c_to_litr_lig_c_col(:,:) = 0._r8
    this%harvest_c_to_cwdc_col(:,:) = 0._r8
    this%m_c_to_litr_met_fire_col(:,:) = 0._r8
    this%m_c_to_litr_cel_fire_col(:,:) = 0._r8
    this%m_c_to_litr_lig_fire_col(:,:) = 0._r8
    this%fire_mortality_c_to_cwdc_col (:,:) = 0._r8
    this%phenology_c_to_litr_lig_c_col(:,:) = 0._r8
  end subroutine initCold

end module CNCarbonFluxType
