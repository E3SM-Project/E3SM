module PatchMod

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval
  use PatchType      , only : patch_type, pft
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

    idx = idx + 1;                                   values(idx) = real( pft%gridcell (p))
    idx = idx + 1; if (.not. isnan(pft%wtgcell(p))) values(idx) =       pft%wtgcell  (p)
    idx = idx + 1;                                   values(idx) = real( pft%landunit (p))
    idx = idx + 1; if (.not. isnan(pft%wtlunit(p))) values(idx) =       pft%wtlunit  (p)
    idx = idx + 1;                                   values(idx) = real( pft%column   (p))
    idx = idx + 1; if (.not. isnan(pft%wtcol(p)))   values(idx) =       pft%wtcol    (p)
    idx = idx + 1;                                   values(idx) = real( pft%itype    (p))
    idx = idx + 1;                                   values(idx) = real( pft%mxy      (p))

    idx = idx + 1; if (pft%active(p))               values(idx) = 1._r8

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

    idx = idx + 1;                           pft%gridcell (p) = int(values(idx))
    idx = idx + 1;                           pft%wtgcell  (p) =     values(idx)
    idx = idx + 1;                           pft%landunit (p) = int(values(idx))
    idx = idx + 1;                           pft%wtlunit  (p) =     values(idx)
    idx = idx + 1;                           pft%column   (p) = int(values(idx))
    idx = idx + 1;                           pft%wtcol    (p) =     values(idx)
    idx = idx + 1;                           pft%itype    (p) = int(values(idx))
    idx = idx + 1;                           pft%mxy      (p) = int(values(idx))
    idx = idx + 1; if (values(idx) == 1._r8) pft%active   (p) = .true.

  end subroutine SetValuesForPatch

  !------------------------------------------------------------------------
  subroutine NumValuesPerPatch(nvalues)
    !
    ! !DESCRIPTION:
    ! Returns the number of values for each patch within the patch_type
    !
    ! !ARGUMENTS:
    integer, intent(out)     :: nvalues

    nvalues = 9

  end subroutine NumValuesPerPatch

end module PatchMod
