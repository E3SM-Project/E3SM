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
  public :: is_hydrologically_active   ! returns true if the given column type is hydrologically active
  public :: icemec_class_to_col_itype  ! convert an icemec class (1..maxpatch_glcmec) into col_pp%itype
  public :: col_itype_to_icemec_class  ! convert col_pp%itype into an icemec class (1..maxpatch_glcmec)

contains
  
  !-----------------------------------------------------------------------
  function is_hydrologically_active(col_itype, lun_itype) &
       result(hydrologically_active)
    !
    ! !DESCRIPTION:
    ! Returns a logical value saying whether the given column type is hydrologically
    ! active
    !
    ! Note that calling this can be bad for performance, because it operates on a single
    ! point rather than a loop. So in performance-critical parts of the code (or just
    ! about anywhere, really), you should use the pre-set col%hydrologically_active(c).
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    logical :: hydrologically_active  ! function result
    integer, intent(in) :: col_itype  ! col%itype value
    integer, intent(in) :: lun_itype  ! lun%itype value for the landunit on which this column sits
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'is_hydrologically_active'
    !-----------------------------------------------------------------------

    ! If we had an easy way to figure out which landunit a column was on based on
    ! col_itype (which would be very helpful!), then we wouldn't need lun_itype.

    if (lun_itype == istsoil .or. lun_itype == istcrop) then
       hydrologically_active = .true.
    else if (col_itype == icol_road_perv) then
       hydrologically_active = .true.
    else
       hydrologically_active = .false.
    end if

  end function is_hydrologically_active

  !-----------------------------------------------------------------------
  function icemec_class_to_col_itype(icemec_class) result(col_itype)
    !
    ! !DESCRIPTION:
    ! Convert an icemec class (1..maxpatch_glcmec) into col_pp%itype
    !
    ! !USES:
    use elm_varpar, only : maxpatch_glcmec
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
    ! Convert a col_pp%itype value (for an icemec landunit) into an icemec class (1..maxpatch_glcmec)
    !
    ! !USES:
    use elm_varpar, only : maxpatch_glcmec
    use landunit_varcon, only : istice_mec
    !
    ! !ARGUMENTS:
    integer :: icemec_class          ! function result
    integer, intent(in) :: col_itype ! col_pp%itype value for an icemec landunit
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
