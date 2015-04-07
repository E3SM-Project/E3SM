module unittestSimpleSubgridSetupsMod

  ! This module provides wrappers to unittestSubgridMod, which give you a variety of
  ! simple subgrid setups.
  !
  ! Note that these routines do everything needed with the subgrid setup. So once you
  ! call these routines, you cannot add any more gridcells, landunits, etc.

  use unittestSubgridMod
  use shr_kind_mod , only : r8 => shr_kind_r8
  use landunit_varcon, only : istsoil

  implicit none
  private
  save

  ! Create a grid that has a single gridcell with a single vegetated patch
  public :: setup_single_veg_patch

contains

  !-----------------------------------------------------------------------
  subroutine setup_single_veg_patch(pft_type)
    !
    ! !DESCRIPTION:
    ! Create a grid that has a single gridcell with a single vegetated patch, with veg
    ! type given by the pft_type argument
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: pft_type  ! the type of the single vegetated patch
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'setup_single_veg_patch'
    !-----------------------------------------------------------------------
    
    call unittest_subgrid_setup_start()
    call unittest_add_gridcell()
    call unittest_add_landunit(my_gi=gi, ltype=istsoil, wtgcell=1.0_r8)
    call unittest_add_column(my_li=li, ctype=1, wtlunit=1.0_r8)
    call unittest_add_patch(my_ci=ci, ptype=pft_type, wtcol=1.0_r8)
    call unittest_subgrid_setup_end()

  end subroutine setup_single_veg_patch

end module unittestSimpleSubgridSetupsMod
