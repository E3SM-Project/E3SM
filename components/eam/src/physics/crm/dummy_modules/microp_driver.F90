module microp_driver
!-------------------------------------------------------------------------------
! Dummy module to avoid building the corresponding module in src/physics/cam
!-------------------------------------------------------------------------------
use shr_kind_mod,   only: r8 => shr_kind_r8
implicit none
private
save
public :: microp_driver_readnl
public :: microp_driver_init_cnst
public :: microp_driver_implements_cnst
!===============================================================================
contains
!===============================================================================
subroutine microp_driver_readnl(nlfile)
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
   ! Read in namelist for microphysics scheme
   !-----------------------------------------------------------------------
   return
end subroutine microp_driver_readnl
!===============================================================================
function microp_driver_implements_cnst(name)
   ! Return true if specified constituent is implemented by the
   ! microphysics package
   character(len=*), intent(in) :: name        ! constituent name
   logical :: microp_driver_implements_cnst    ! return value
   !-----------------------------------------------------------------------
   microp_driver_implements_cnst = .false.
end function microp_driver_implements_cnst
!===============================================================================
subroutine microp_driver_init_cnst(name, q, gcid)
   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.
   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------
   return
end subroutine microp_driver_init_cnst
!===============================================================================
end module microp_driver
