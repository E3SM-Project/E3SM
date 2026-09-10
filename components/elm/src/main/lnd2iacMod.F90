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
  use lnd2atmType     , only : lnd2atm_type ! for  lnd2atm_vars
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! Base temperature for degree day calculation: 18C in Kelvin
  real(r8), parameter :: t_base_K = 291.15_r8

  ! lnd -> iac variables structure
  ! Fields are dimensioned (ngrid,numpft+1)
  ! pftwgt is frac of actual grid cell (not frac of land)
  type, public :: lnd2iac_type
     real(r8), pointer :: hr(:,:) => null()
     real(r8), pointer :: npp(:,:) => null()
     real(r8), pointer :: pftwgt(:,:) => null()
     real(r8), pointer :: forc_hdm(:) => null()
     real(r8), pointer :: hdd(:) => null()
     real(r8), pointer :: cdd(:) => null()

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
    allocate(this%hdd(begg:endg))
    allocate(this%cdd(begg:endg))

    this%hr(:,:)=0.0_r8
    this%npp(:,:)=0.0_r8
    this%pftwgt(:,:)=0.0_r8
    this%forc_hdm(:)= 0.0_r8
    this%hdd(:)= 0.0_r8
    this%cdd(:)= 0.0_r8

  end subroutine Init

  !------------------------------------------------------
  subroutine update_lnd2iac(this, bounds, lnd2atm_vars)
    ! !DESCRIPTION:
    ! Stuff values into lnd2iac.
    ! Instantaneous fields (hr, npp, pftwgt, hdd, cdd) are overwritten each call.
    !
    ! !ARGUMENTS:
    class(lnd2iac_type)       , intent(inout) :: this
    type(bounds_type)         , intent(in)   :: bounds
    type(lnd2atm_type)        , intent(in)   :: lnd2atm_vars ! elm land to atmosphere exchange data type

    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'update_lnd2iac'
    integer :: begg,endg
    integer :: begp,endp
    integer :: c,p,g,pft

    begg = bounds%begg; endg = bounds%endg
    begp = bounds%begp; endp = bounds%endp

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

      end if
    end do

    do g = begg, endg
      ! Estimate heating and cooling degree days
      this%hdd(g) = max(t_base_K - lnd2atm_vars%t_ref2m_grc(g), 0.0_r8)
      this%cdd(g) = max(lnd2atm_vars%t_ref2m_grc(g) - t_base_K, 0.0_r8)
    end do

  end subroutine update_lnd2iac
end module lnd2iacMod
