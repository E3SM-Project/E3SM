! physics_utils: common elements for SCREAM physics models.
module physics_utils

#ifdef SCREAM_CONFIG_IS_CMAKE
  use iso_c_binding, only: c_double, c_float, c_bool
#else
  use shr_kind_mod,   only: rtype=>shr_kind_r8, itype=>shr_kind_i8
#endif

  implicit none
  private
  public:: calculate_drymmr_from_wetmmr, calculate_wetmmr_from_drymmr, &
       calculate_drymmr_from_wetmmr_1gridcell, calculate_wetmmr_from_drymmr_1gridcell
  save

#ifdef SCREAM_CONFIG_IS_CMAKE
#include "scream_config.f"

  integer,parameter,public :: rtype8 = c_double ! 8 byte real, compatible with c type double
  integer,parameter,public :: btype  = c_bool ! boolean type, compatible with c
  integer,parameter,public :: itype = selected_int_kind (13) ! 8 byte integer

#  ifdef SCREAM_DOUBLE_PRECISION
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
#  else
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type float
#  endif

#else
  integer,parameter,public :: btype = kind(.true.) ! native logical
  public :: rtype
  public :: itype
  integer,parameter,public :: rtype8 = selected_real_kind(15, 307) ! 8 byte real, compatible with c type double

#endif

contains

  pure function calculate_drymmr_from_wetmmr(ncols, pver, wetmmr, qv_wet) result (drymmr)
    !Computes drymmr (mass of a constituent divided by mass of dry air; commonly known as mixing ratio)
    !for any wetmmr constituent (mass of a constituent divided by mass of dry air plus water
    !vapor) using qv_wet (mass of water vapor divided by mass of dry air plus
    !water vapor; see specific humidity).

    implicit none

    !intent-ins
    integer,     intent(in) :: ncols, pver !number of columns and levels
    real(rtype), intent(in) :: wetmmr(:,:) !wet mmr of a constituent
    real(rtype), intent(in) :: qv_wet(:,:) !water vapor wet mass mixing ratio

    !return variable
    real(rtype) :: drymmr(ncols,pver) !dry mmr of a constituent

    integer :: klev, icol

    !Assign uninitialized columns using "huge"
    drymmr(:,:) = huge(1.0_rtype)

    !Compute drymmr
    do klev = 1, pver
       do icol = 1, ncols
          drymmr(icol,klev) = calculate_drymmr_from_wetmmr_1gridcell(wetmmr(icol,klev), qv_wet(icol,klev))
       enddo
    enddo

  end function calculate_drymmr_from_wetmmr

  !------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------

  pure function calculate_drymmr_from_wetmmr_1gridcell(wetmmr,qv_wet) result(drymmr)
    !Computes drymmr from wetmmr for a grid cell
    implicit none

    !intent-ins
    real(rtype), intent(in) :: wetmmr !wet mmr of a constituent
    real(rtype), intent(in) :: qv_wet !water vapor wet mass mixing ratio

    !return variable
    real(rtype) :: drymmr !dry mmr of a constituent

    drymmr = wetmmr/(1.0_rtype - qv_wet)

  end function calculate_drymmr_from_wetmmr_1gridcell

  !------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------

  pure function calculate_wetmmr_from_drymmr(ncols, pver, drymmr, qv_dry) result (wetmmr)

    !Computes wetmmr (mass of a constituent divided by mass of dry air plus water vapor)
    !for any drymmr constituent (mass of a constituent divided by mass of dry air;
    !commonly known as mixing ratio) using qv_dry (mass of water vapor divided by mass
    !of dry air)

    implicit none

    !intent-ins
    integer,     intent(in) :: ncols, pver !number of columns and levels
    real(rtype), intent(in) :: drymmr(:,:) !dry mmr of a constituent
    real(rtype), intent(in) :: qv_dry(:,:) !water vapor dry mass mixing ratio

    !return variable
    real(rtype) :: wetmmr(ncols,pver) !wet mmr of a constituent

    integer :: klev, icol

    !Assign uninitialized columns using "huge"
    wetmmr(:,:) = huge(1.0_rtype)

    !Compute drymmr
    do klev = 1, pver
       do icol = 1, ncols
          !Compute wetmmr
          wetmmr(icol,klev) = calculate_wetmmr_from_drymmr_1gridcell(drymmr(icol,klev), qv_dry(icol,klev))
       enddo
    enddo

  end function calculate_wetmmr_from_drymmr

  !------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------

  pure function calculate_wetmmr_from_drymmr_1gridcell(drymmr, qv_dry) result (wetmmr)
    !Computes wetmmr from drymmr for a grid cell
    implicit none

    !intent-ins
    real(rtype), intent(in) :: drymmr !dry mmr of a constituent
    real(rtype), intent(in) :: qv_dry !water vapor dry mass mixing ratio

    !return variable
    real(rtype) :: wetmmr !wet mmr of a constituent

    !Compute wetmmr
    wetmmr = drymmr/(1.0_rtype + qv_dry)

  end function calculate_wetmmr_from_drymmr_1gridcell

end module physics_utils
