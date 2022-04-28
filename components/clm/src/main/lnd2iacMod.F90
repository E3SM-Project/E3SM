module lnd2iacMod

  !------------------------------------------------
  ! !DESCRIPTION:
  ! This module deals with arrays for exchanging data from land model
  ! to iac.  We need a separate field per pft, because we can only
  ! couple grid values; so our fields will be allocated (ngrid,npft)
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : get_proc_bounds, bounds_type
  use domainMod       , only : ldomain
  use clm_varctl      , only : iulog
  use clm_varpar      , only : numpft
  use abortutils      , only : endrun
  use ColumnType      , only : col_pp   ! 
  use ColumnDataType  , only : col_cf   ! for hr
  use VegetationType  , only : veg_pp   ! pftwgt
  use VegetationDataType, only: veg_cf  ! for npp
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! lnd -> iac variables structure
  ! Fields are dimensioned (ngrid,numpft+1)
  ! pftwgt is frac of actual grid cell (not frac of land)
  type, public :: lnd2iac_type
     real(r8), pointer :: hr(:,:) => null()
     real(r8), pointer :: npp(:,:) => null()
     real(r8), pointer :: pftwgt(:,:) => null()

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

    this%hr(:,:)=0.0_r8
    this%npp(:,:)=0.0_r8
    this%pftwgt(:,:)=0.0_r8

  end subroutine Init

  !------------------------------------------------------
  subroutine update_lnd2iac(this, bounds)
    ! !DESCRIPTION:
    ! Stuff values into lnd2iac
    !
    ! !ARGUMENTS:
    class(lnd2iac_type)       , intent(inout) :: this
    type(bounds_type)         , intent(in)   :: bounds

    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'update_lnd2iac'
    integer :: begg,endg
    integer :: begc,endc
    integer :: begp,endp

    integer :: c,p,g,pft

    begg = bounds%begg; endg = bounds%endg
    begp = bounds%begp; endp = bounds%endp

    ! Fill everything with zeros by default
    this%hr(begg:endg,:)=0.0_r8
    this%npp(begg:endg,:)=0.0_r8
    this%pftwgt(begg:endg,:)=0.0_r8

    ! Loop over patch index, extract fields by pft type and gridcell
    do p = begp, endp
       g=veg_pp%gridcell(p)
       pft=veg_pp%itype(p)

       c=veg_pp%column(p) ! for hr

       ! Assign values
       !write(iulog,*) 'TRS0: ', p, c, g, pft, begp, endp, begg, endg
       !write(iulog,*) 'TRS3: ', col_cf%hr(:)
       this%hr(g,pft) = col_cf%hr(c)   ! Every pft in this column gets this hr value
       this%npp(g,pft) = veg_cf%npp(p)
       ! this is the fraction of actual grid cell
       this%pftwgt(g,pft) = veg_pp%wtgcell(p) * ldomain%frac(g) * &
                             ldomain%mask(g)
    end do

    ! Ta da

  end subroutine update_lnd2iac
end module lnd2iacMod
