module column_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing landunit indices and associated variables and routines.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use landunit_varcon, only : isturb_MIN
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  !------------------------------------------------------------------
  ! Initialize column type constants
  !------------------------------------------------------------------

  ! urban column types

  integer, parameter, public :: icol_roof        = isturb_MIN*10 + 1
  integer, parameter, public :: icol_sunwall     = isturb_MIN*10 + 2
  integer, parameter, public :: icol_shadewall   = isturb_MIN*10 + 3
  integer, parameter, public :: icol_road_imperv = isturb_MIN*10 + 4
  integer, parameter, public :: icol_road_perv   = isturb_MIN*10 + 5

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: icemec_class_to_col_itype  ! convert an icemec class (1..maxpatch_glcmec) into col%itype
  public :: col_itype_to_icemec_class  ! convert col%itype into an icemec class (1..maxpatch_glcmec)

contains
  
  !-----------------------------------------------------------------------
  function icemec_class_to_col_itype(icemec_class) result(col_itype)
    !
    ! !DESCRIPTION:
    ! Convert an icemec class (1..maxpatch_glcmec) into col%itype
    !
    ! !USES:
    use clm_varpar, only : maxpatch_glcmec
    use landunit_varcon, only : istice_mec
    !
    ! !ARGUMENTS:
    integer :: col_itype                ! function result
    integer, intent(in) :: icemec_class ! icemec class, between 1 and maxpatch_glcmec
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'icemec_class_to_col_itype'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT((1 <= icemec_class .and. icemec_class <= maxpatch_glcmec), errMsg(__FILE__, __LINE__))

    col_itype = istice_mec*100 + icemec_class

  end function icemec_class_to_col_itype

  !-----------------------------------------------------------------------
  function col_itype_to_icemec_class(col_itype) result(icemec_class)
    !
    ! !DESCRIPTION:
    ! Convert a col%itype value (for an icemec landunit) into an icemec class (1..maxpatch_glcmec)
    !
    ! !USES:
    use clm_varpar, only : maxpatch_glcmec
    use landunit_varcon, only : istice_mec
    !
    ! !ARGUMENTS:
    integer :: icemec_class          ! function result
    integer, intent(in) :: col_itype ! col%itype value for an icemec landunit
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'col_itype_to_icemec_class'
    !-----------------------------------------------------------------------
    
    icemec_class = col_itype - istice_mec*100

    ! The following assertion is here to ensure that col_itype is really from an
    ! istice_mec landunit
    SHR_ASSERT((1 <= icemec_class .and. icemec_class <= maxpatch_glcmec), errMsg(__FILE__, __LINE__))

  end function col_itype_to_icemec_class

end module column_varcon
