! physics_utils: common elements for SCREAM physics models.
module physics_utils

#ifdef SCREAM_CONFIG_IS_CMAKE
  use iso_c_binding, only: c_double, c_float, c_bool
#else
  use shr_kind_mod,   only: rtype=>shr_kind_r8, itype=>shr_kind_i8
#endif

  implicit none
  private
  public:: calculate_drymmr_from_wetmmr, calculate_wetmmr_from_drymmr
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

  pure function calculate_drymmr_from_wetmmr(ncols, wetmmr, qv_wet) result (drymmr)
    !Compute drymmr for any wetmmr constituent using wet water vapor mixing ratio(qv_wet)

    use ppgrid, only: pcols, pver

    implicit none

    !intent-ins
    integer,     intent(in) :: ncols       !number of columns
    real(rtype), intent(in) :: wetmmr(:,:) !wet mmr of a constituent
    real(rtype), intent(in) :: qv_wet(:,:) !water vapor wet mass mixing ratio

    !return variable
    real(rtype) :: drymmr(pcols,pver) !dry mmr of a constituent

    !Compute drymmr
    drymmr(:ncols,:) = wetmmr(:ncols,:)/(1.0_rtype - qv_wet(:ncols,:))

    !Since pcols can be > ncols, assign uninitialized columns using "huge"
    drymmr(ncols+1:pcols,:) = huge(1.0_rtype)

  end function calculate_drymmr_from_wetmmr


  pure function calculate_wetmmr_from_drymmr(ncols, drymmr, qv_dry) result (wetmmr)
    !Compute wetmmr for any drymmr constituent using dry water vapor mixing ratio(qv_dry)
    
    use ppgrid, only: pcols, pver

    implicit none

    !intent-ins
    integer,     intent(in) :: ncols      !number of columns
    real(rtype), intent(in) :: drymmr(:,:)!dry mmr of a constituent
    real(rtype), intent(in) :: qv_dry(:,:)!water vapor dry mass mixing ratio

    !return variable
    real(rtype) :: wetmmr(pcols,pver) !wet mmr of a constituent

    !Compute wetmmr
    wetmmr(:ncols,:) = drymmr(:ncols,:)/(1.0_rtype + qv_dry(:ncols,:))

    !Since pcols can be > ncols, assign uninitialized columns using "huge"
    wetmmr(ncols+1:pcols,:) = huge(1.0_rtype)

  end function calculate_wetmmr_from_drymmr

end module physics_utils
