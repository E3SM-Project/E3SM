module PhosphorusFluxType
  use clm_varcon             , only : spval, ispval
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use clm_varpar             , only : nlevdecomp_full, ndecomp_pools
implicit none

  type, public :: phosphorusflux_type
    real(r8), pointer :: sminp_leached_col                         (:)    => null() !col inorganic P leaching loss, gP/m2/time step
    real(r8), pointer :: sminp_runoff_col                          (:)   => null()  !col inorganic P runoff loss, gP/m2/time step
    real(r8), pointer :: biochem_pmin_vr_col                       (:,:) => null()  ! col vertically-resolved total biochemical P mineralization (gP/m3/s)
    real(r8), pointer :: phenology_p_to_litr_met_p_col             (:,:)=> null()
    real(r8), pointer :: dwt_frootp_to_litr_met_p_col              (:,:)=> null()
    real(r8), pointer :: phenology_p_to_litr_cel_p_col             (:,:)=> null()
    real(r8), pointer :: dwt_livecrootp_to_cwdp_col                (:,:)=> null()
    real(r8), pointer :: phenology_p_to_litr_lig_p_col             (:,:)=> null()
    real(r8), pointer :: m_decomp_ppools_to_fire_vr_col            (:,:,:)=> null()
    real(r8), pointer :: dwt_deadcrootp_to_cwdp_col                (:,:)=> null()
    real(r8), pointer :: dwt_frootp_to_litr_lig_p_col              (:,:)=> null()
    real(r8), pointer :: dwt_frootp_to_litr_cel_p_col              (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_litr_met_p_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_cwdp_col               (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_litr_lig_p_col         (:,:)=> null()
    real(r8), pointer :: gap_mortality_p_to_litr_cel_p_col         (:,:)=> null()
    real(r8), pointer :: harvest_p_to_litr_met_p_col               (:,:)=> null()
    real(r8), pointer :: harvest_p_to_cwdp_col                     (:,:)=> null()
    real(r8), pointer :: harvest_p_to_litr_lig_p_col               (:,:)=> null()
    real(r8), pointer :: harvest_p_to_litr_cel_p_col               (:,:)=> null()
    real(r8), pointer :: m_p_to_litr_met_fire_col                  (:,:)=> null()
    real(r8), pointer :: m_p_to_litr_cel_fire_col                  (:,:)=> null()
    real(r8), pointer :: m_p_to_litr_lig_fire_col                  (:,:)=> null()
    real(r8), pointer :: fire_mortality_p_to_cwdp_col              (:,:)=> null()
    real(r8), pointer :: primp_to_labilep_vr_col                   (:,:)=> null()
    real(r8), pointer :: pdep_to_sminp_col                         (:)=> null()
    real(r8), pointer :: fert_p_to_sminp_col                       (:)=> null()
    real(r8), pointer :: sminp_to_plant_patch                      (:)=> null()
    real(r8), pointer :: supplement_to_sminp_col                   (:)=> null()
    real(r8), pointer :: secondp_to_occlp_col                      (:)=> null()
    real(r8), pointer :: fire_decomp_ploss_col                     (:)=> null()
    real(r8), pointer :: som_p_leached_col                         (:) => null()

  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
    procedure, private :: SetValues
  end type phosphorusflux_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(phosphorusflux_type) :: this
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
    class(phosphorusflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%sminp_leached_col     (begc:endc              ))     ;this%sminp_leached_col          (:)    = nan
    allocate(this%sminp_runoff_col      (begc:endc              ))      ;this%sminp_runoff_col         (:)    = nan
    allocate(this%biochem_pmin_vr_col   (begc:endc,1:nlevdecomp_full))   ;this%biochem_pmin_vr_col      (:,:) = nan

    allocate(this%phenology_p_to_litr_met_p_col(begc:endc, 1:nlevdecomp_full)); this%phenology_p_to_litr_met_p_col(:,:) = nan
    allocate(this%dwt_frootp_to_litr_met_p_col(begc:endc, 1:nlevdecomp_full)); this%dwt_frootp_to_litr_met_p_col(:,:) = nan
    allocate(this%phenology_p_to_litr_cel_p_col(begc:endc,1:nlevdecomp_full)); this%phenology_p_to_litr_cel_p_col(:,:) = nan
    allocate(this%dwt_livecrootp_to_cwdp_col(begc:endc,1:nlevdecomp_full)); this%dwt_livecrootp_to_cwdp_col(:,:) = nan
    allocate(this%phenology_p_to_litr_lig_p_col(begc:endc, 1:nlevdecomp_full)); this%phenology_p_to_litr_lig_p_col(:,:) = nan
    allocate(this%m_decomp_ppools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:7)); this%m_decomp_ppools_to_fire_vr_col(:,:,:) = nan
    allocate(this%dwt_deadcrootp_to_cwdp_col(begc:endc, 1:nlevdecomp_full)); this%dwt_deadcrootp_to_cwdp_col(:,:) = nan
    allocate(this%dwt_frootp_to_litr_lig_p_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootp_to_litr_lig_p_col(:,:) = nan
    allocate(this%dwt_frootp_to_litr_cel_p_col(begc:endc,1:nlevdecomp_full)); this%dwt_frootp_to_litr_cel_p_col(:,:) = nan
    allocate(this%gap_mortality_p_to_litr_met_p_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_litr_met_p_col(:,:) = nan
    allocate(this%gap_mortality_p_to_cwdp_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_cwdp_col(:,:) = nan
    allocate(this%gap_mortality_p_to_litr_lig_p_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_litr_lig_p_col(:,:) = nan
    allocate(this%gap_mortality_p_to_litr_cel_p_col(begc:endc,1:nlevdecomp_full)); this%gap_mortality_p_to_litr_cel_p_col(:,:) = nan
    allocate(this%harvest_p_to_litr_met_p_col(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_litr_met_p_col(:,:) = nan
    allocate(this%harvest_p_to_cwdp_col(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_cwdp_col(:,:) = nan
    allocate(this%harvest_p_to_litr_lig_p_col(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_litr_lig_p_col(:,:) = nan
    allocate(this%harvest_p_to_litr_cel_p_col(begc:endc,1:nlevdecomp_full)); this%harvest_p_to_litr_cel_p_col(:,:) = nan
    allocate(this%m_p_to_litr_met_fire_col(begc:endc,1:nlevdecomp_full)); this%m_p_to_litr_met_fire_col(:,:) = nan
    allocate(this%m_p_to_litr_cel_fire_col(begc:endc,1:nlevdecomp_full)); this%m_p_to_litr_cel_fire_col(:,:) = nan
    allocate(this%m_p_to_litr_lig_fire_col(begc:endc,1:nlevdecomp_full)); this%m_p_to_litr_lig_fire_col(:,:) = nan
    allocate(this%fire_mortality_p_to_cwdp_col(begc:endc,1:nlevdecomp_full)); this%fire_mortality_p_to_cwdp_col(:,:) = nan
    allocate(this%primp_to_labilep_vr_col(begc:endc,1:nlevdecomp_full)); this%primp_to_labilep_vr_col(:,:) = nan
    allocate(this%pdep_to_sminp_col(begc:endc)); this%pdep_to_sminp_col(begc:endc) = nan
    allocate(this%pdep_to_sminp_col(begc:endc)); this%pdep_to_sminp_col(begc:endc) = nan
    allocate(this%fert_p_to_sminp_col(begc:endc)); this%fert_p_to_sminp_col(begc:endc) = nan
    allocate(this%som_p_leached_col(begc:endc)); this%som_p_leached_col(:) = nan

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
    class(phosphorusflux_type) :: this
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


    this%phenology_p_to_litr_met_p_col(:,:) = 0._r8
    this%dwt_frootp_to_litr_met_p_col(:,:) = 0._r8
    this%phenology_p_to_litr_cel_p_col(:,:) = 0._r8
    this%dwt_livecrootp_to_cwdp_col(:,:) = 0._r8
    this%phenology_p_to_litr_lig_p_col(:,:) = 0._r8
    this%m_decomp_ppools_to_fire_vr_col(:,:,:) = 0._r8
    this%dwt_deadcrootp_to_cwdp_col(:,:) = 0._r8
    this%dwt_frootp_to_litr_lig_p_col(:,:) = 0._r8
    this%dwt_frootp_to_litr_cel_p_col(:,:) = 0._r8
    this%gap_mortality_p_to_litr_met_p_col(:,:) = 0._r8
    this%gap_mortality_p_to_cwdp_col(:,:) = 0._r8
    this%gap_mortality_p_to_litr_lig_p_col(:,:) = 0._r8
    this%gap_mortality_p_to_litr_cel_p_col(:,:) = 0._r8
    this%harvest_p_to_litr_met_p_col(:,:) = 0._r8
    this%harvest_p_to_cwdp_col(:,:) = 0._r8
    this%harvest_p_to_litr_lig_p_col(:,:) = 0._r8
    this%harvest_p_to_litr_cel_p_col(:,:) = 0._r8
    this%m_p_to_litr_met_fire_col(:,:) = 0._r8
    this%m_p_to_litr_cel_fire_col(:,:) = 0._r8
    this%m_p_to_litr_lig_fire_col(:,:) = 0._r8
    this%fire_mortality_p_to_cwdp_col(:,:) = 0._r8
    this%primp_to_labilep_vr_col(:,:) = 0._r8
    this%pdep_to_sminp_col(:) = 0._r8
    this%pdep_to_sminp_col(:) = 0._r8
    this%fert_p_to_sminp_col(:) = 0._r8

  end subroutine initCold


  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class (phosphorusflux_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column

    integer :: fi, i

    do fi = 1,num_column
       i = filter_column(fi)

    enddo
  end subroutine SetValues
end module PhosphorusFluxType
