module lnd2iacMod

  !------------------------------------------------
  ! !DESCRIPTION:
  ! This module deals with arrays for exchanging data from land model
  ! to iac.  We need a separate field per pft, because we can only
  ! couple grid values; so our fields will be allocated (ngrid,npft)
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use domainMod       , only : ldomain
  use elm_varpar      , only : numpft
  use ColumnDataType  , only : col_cf   ! for hr
  use VegetationType  , only : veg_pp   ! pftwgt
  use VegetationDataType, only: veg_cf  ! for npp
  use VegetationDataType, only: veg_es  ! for t_ref2m
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! Base temperature for degree day accumulation: 18C in Kelvin
  real(r8), parameter :: t_base_K = 291.15_r8

  ! lnd -> iac variables structure
  ! Fields are dimensioned (ngrid,numpft+1)
  ! pftwgt is frac of actual grid cell (not frac of land)
  type, public :: lnd2iac_type
     real(r8), pointer :: hr(:,:) => null()
     real(r8), pointer :: npp(:,:) => null()
     real(r8), pointer :: pftwgt(:,:) => null()
     real(r8), pointer :: forc_hdm(:)    => null()
     ! Running cumulative degree day accumulators (K-days) since model start.
     ! Updated every ELM timestep; never reset in ELM.
     ! The IAC reads these each year and diffs consecutive snapshots to
     ! obtain annual increments, then averages over the GCAM period.
     !   HDD_accum += max(t_base_K - veg_es%t_ref2m, 0) * (dt/86400)  [heating]
     !   CDD_accum += max(veg_es%t_ref2m - t_base_K, 0) * (dt/86400)  [cooling]
     real(r8), pointer :: HDD_accum(:,:) => null()
     real(r8), pointer :: CDD_accum(:,:) => null()

   contains
     ! This object oriented stuff...
     procedure, public  :: Init
     procedure, public  :: update_lnd2iac
  end type lnd2iac_type

  ! !PUBLIC MEMBER FUNCTIONS:

contains

  !-------------------------------
  subroutine Init(this, bounds)

    ! !DESCRIPTION:
    ! Initialize land variables required by glc
    !
    ! !ARGUMENTS:
    class(lnd2iac_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    integer :: begg,endg

    begg = bounds%begg; endg = bounds%endg

    ! We add in bare ground as pft=0
    allocate(this%hr(begg:endg,0:numpft)) 
    allocate(this%npp(begg:endg,0:numpft))
    allocate(this%pftwgt(begg:endg,0:numpft)) 
    allocate(this%forc_hdm(begg:endg))
    allocate(this%HDD_accum(begg:endg,0:numpft))
    allocate(this%CDD_accum(begg:endg,0:numpft))

    this%hr(:,:)=0.0_r8
    this%npp(:,:)=0.0_r8
    this%pftwgt(:,:)=0.0_r8
    this%forc_hdm(:)   = 0.0_r8
    this%HDD_accum(:,:)= 0.0_r8
    this%CDD_accum(:,:)= 0.0_r8

  end subroutine Init

  !------------------------------------------------------
  subroutine update_lnd2iac(this, bounds)
    ! !DESCRIPTION:
    ! Stuff values into lnd2iac.
    ! Instantaneous fields (hr, npp, pftwgt, veg_es%t_ref2m) are overwritten each call.
    ! HDD_accum and CDD_accum are running cumulative sums from model start -
    ! they are NEVER reset in ELM.  The IAC diffs consecutive exported
    ! snapshots to compute annual increments and multi-year averages.
    !   HDD_accum += max(t_base_K - veg_es%t_ref2m, 0) * (dt/86400)  [K-days]
    !   CDD_accum += max(veg_es%t_ref2m - t_base_K, 0) * (dt/86400)  [K-days]
    !
    ! !ARGUMENTS:
    use elm_time_manager, only: get_step_size
    class(lnd2iac_type)       , intent(inout) :: this
    type(bounds_type)         , intent(in)   :: bounds

    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'update_lnd2iac'
    integer :: begg,endg
    integer :: begp,endp
    integer :: c,p,g,pft
    real(r8) :: dt_days  ! ELM timestep length in days

    begg = bounds%begg; endg = bounds%endg
    begp = bounds%begp; endp = bounds%endp

    dt_days = real(get_step_size(), r8) / 86400.0_r8

    ! Zero instantaneous fields; accumulators retain their running totals
    this%hr(begg:endg,:)=0.0_r8
    this%npp(begg:endg,:)=0.0_r8
    this%pftwgt(begg:endg,:)=0.0_r8

    ! Loop over patch index, extract fields by pft type and gridcell
    do p = begp, endp
       g=veg_pp%gridcell(p)
       pft=veg_pp%itype(p)
       c=veg_pp%column(p) ! for hr

       if (veg_pp%active(p)) then
         ! Instantaneous fields
         this%hr(g,pft)      = col_cf%hr(c)
         this%npp(g,pft)     = veg_cf%npp(p)
         this%pftwgt(g,pft)  = veg_pp%wtgcell(p) * ldomain%frac(g) * &
                               ldomain%mask(g)

         ! Accumulate degree day exceedances for this timestep (never reset)
         this%HDD_accum(g,pft) = this%HDD_accum(g,pft) + &
            max(t_base_K - veg_es%t_ref2m(p), 0.0_r8) * dt_days
         this%CDD_accum(g,pft) = this%CDD_accum(g,pft) + &
            max(veg_es%t_ref2m(p) - t_base_K, 0.0_r8) * dt_days
      end if
    end do
  end subroutine update_lnd2iac
end module lnd2iacMod
