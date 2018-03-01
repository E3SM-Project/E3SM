module CNStateType
  use clm_varcon     , only : spval, ispval, c14ratio
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varpar     , only : nlevdecomp_full
implicit none

  type, public :: cnstate_type
     real(r8), pointer :: rc14_atm_patch               (:)   => null()  ! patch C14O2/C12O2 in atmosphere
     real(r8), pointer :: froot_prof_patch             (:,:)=> null()
     real(r8), pointer :: nfixation_prof_col           (:,:) => null()  ! col (1/m) profile for N fixation additions
     integer  ,pointer :: isoilorder                   (:)  => null()   ! col global soil order data
     real(r8), pointer :: cn_scalar                    (:)  => null()   ! cn scaling factor for root n uptake kinetics (no units)
     real(r8), pointer :: cp_scalar                    (:)   => null()  ! cp scaling factor for root p uptake kinetics (no units)
     real(r8), pointer :: ndep_prof_col                (:,:)=> null()
     real(r8), pointer :: pdep_prof_col                (:,:)=> null()
     real(r8), pointer :: frac_loss_lit_to_fire_col    (:) => null()
     real(r8), pointer :: frac_loss_cwd_to_fire_col    (:) => null()
  contains

    procedure, public  :: Init
    procedure, private :: InitCold
    procedure, private :: InitAllocate
  end type cnstate_type


contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(cnstate_type) :: this
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
    class(cnstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    allocate(this%rc14_atm_patch              (begp:endp)) ;    this%rc14_atm_patch              (:) = nan
    allocate(this%froot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%froot_prof_patch    (:,:) = spval
    allocate(this%nfixation_prof_col  (begc:endc,1:nlevdecomp_full)) ; this%nfixation_prof_col  (:,:) = spval
    allocate(this%isoilorder            (begc:endc))                 ; this%isoilorder(:) = ispval
    allocate(this%ndep_prof_col (begc:endc, 1:nlevdecomp_full)); this%ndep_prof_col(:,:) = spval
    allocate(this%pdep_prof_col (begc:endc, 1:nlevdecomp_full)); this%pdep_prof_col(:,:) = spval

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
    class(cnstate_type) :: this
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


    do p = bounds%begp,bounds%endp
      this%rc14_atm_patch(p)              = c14ratio
    enddo
    this%ndep_prof_col(:,:) = 0._r8
    this%nfixation_prof_col(:,:) = 0._r8
  end subroutine initCold

end module CNStateType
