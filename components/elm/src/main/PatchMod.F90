module PatchMod

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varcon     , only : ispval
  use VegetationType , only : veg_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  public :: GetValuesForPatch
  public :: SetValuesForPatch
  public :: NumValuesPerPatch

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine GetValuesForPatch(p, values)
    !
    ! !DESCRIPTION:
    ! Returns all values for the p-th patch
    !
    ! NOTE: The order of the values retrieved/stored in
    !       GetValuesForPatch/SetValuesForPatch should be same.
    !
    ! !USES
    use abortutils       , only : endrun
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    integer  , intent(in)           :: p
    real(r8) , intent(out), pointer :: values(:)
    !
    ! !LOCAL VARIABLES:
    integer                         :: idx
    integer                         :: j
    integer                         :: nvalues

    ! Check the size of data
    call NumValuesPerPatch(nvalues)
    if (size(values) /= nvalues) then
       call endrun(msg="ERROR Size of data is incorrect "//errmsg(__FILE__, __LINE__))
    endif

    idx = 0

    idx = idx + 1;                                   values(idx) = real( veg_pp%gridcell (p))
    idx = idx + 1; if (.not. isnan(veg_pp%wtgcell(p))) values(idx) =       veg_pp%wtgcell  (p)
    idx = idx + 1;                                   values(idx) = real( veg_pp%landunit (p))
    idx = idx + 1; if (.not. isnan(veg_pp%wtlunit(p))) values(idx) =       veg_pp%wtlunit  (p)
    idx = idx + 1;                                   values(idx) = real( veg_pp%column   (p))
    idx = idx + 1; if (.not. isnan(veg_pp%wtcol(p)))   values(idx) =       veg_pp%wtcol    (p)
    idx = idx + 1;                                   values(idx) = real( veg_pp%itype    (p))
    idx = idx + 1;                                   values(idx) = real( veg_pp%mxy      (p))

    idx = idx + 1; if (veg_pp%active(p))               values(idx) = 1._r8

  end subroutine GetValuesForPatch

  !------------------------------------------------------------------------
  subroutine SetValuesForPatch(p, values)
    !
    ! !DESCRIPTION:
    ! Sets all values for the p-th patch.
    !
    ! NOTE: The order of the values retrieved/stored in
    !       GetValuesForPatch/SetValuesForPatch should be same.
    !
    ! !USES
    use abortutils       , only : endrun
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    integer  , intent(in)          :: p
    real(r8) , intent(in), pointer :: values(:)
    !
    integer                        :: idx
    integer                        :: j
    integer                        :: nvalues

    ! Check the size of data
    call NumValuesPerPatch(nvalues)
    if (size(values) /= nvalues) then
       call endrun(msg="ERROR Size of data is incorrect "//errmsg(__FILE__, __LINE__))
    endif

    idx = 0

    idx = idx + 1;                           veg_pp%gridcell (p) = int(values(idx))
    idx = idx + 1;                           veg_pp%wtgcell  (p) =     values(idx)
    idx = idx + 1;                           veg_pp%landunit (p) = int(values(idx))
    idx = idx + 1;                           veg_pp%wtlunit  (p) =     values(idx)
    idx = idx + 1;                           veg_pp%column   (p) = int(values(idx))
    idx = idx + 1;                           veg_pp%wtcol    (p) =     values(idx)
    idx = idx + 1;                           veg_pp%itype    (p) = int(values(idx))
    idx = idx + 1;                           veg_pp%mxy      (p) = int(values(idx))
    idx = idx + 1; if (values(idx) == 1._r8) veg_pp%active   (p) = .true.

  end subroutine SetValuesForPatch

  !------------------------------------------------------------------------
  subroutine NumValuesPerPatch(nvalues)
    !
    ! !DESCRIPTION:
    ! Returns the number of values for each patch within the patch_physical_properties_type
    !
    ! !ARGUMENTS:
    integer, intent(out)     :: nvalues

    nvalues = 9

  end subroutine NumValuesPerPatch

end module PatchMod
