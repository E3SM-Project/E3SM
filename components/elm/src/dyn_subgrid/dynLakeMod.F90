module dynLakeMod

#include "shr_assert.h"

  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use shr_kind_mod , only : r8 => shr_kind_r8
  use landunit_varcon    , only : istdlak
  use LandunitType , only : lun_pp
  use ColumnType   , only : col_pp
  use decompMod    , only : bounds_type

  implicit none

  public :: dynlake_driver

contains

  subroutine dynlake_driver (bounds)
    !---------------------------------------------------------------------------
    ! Dynamic Lake Module
    !
    ! This module updates lake states and fluxes within the dynamic subgrid
    ! framework.
    !
    ! Input:
    !   bounds - bounds structure for the current processor
    !
    !---------------------------------------------------------------------------
    type(bounds_type), intent(in) :: bounds
    !
    integer :: c, l

    ! Loop over all columns assigned to this processor
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istdlak) then
          ! Update lake states and fluxes for this column
          if (col_pp%is_soil(c)) then
             col_pp%wtlunit(c) = col_pp%wtlunit(c) + 0.0001_r8
          endif
          if (col_pp%is_lake(c)) then
             col_pp%wtlunit(c) = col_pp%wtlunit(c) - 0.0001_r8
          endif
       end if
    end do

  end subroutine dynlake_driver

end module dynLakeMod
