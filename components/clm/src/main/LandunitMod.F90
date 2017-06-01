module LandunitMod

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval
  use LandunitType   , only : landunit_type, lun
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  public :: GetValuesForLandunit
  public :: SetValuesForLandunit
  public :: NumValuesPerLandunit

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine GetValuesForLandunit(l, values)
    !
    ! !DESCRIPTION:
    ! Returns all values for the l-th landunit
    !
    ! NOTE: The order of the values retrieved/stored in
    !       GetValuesForLandunit/SetValuesForLandunit should be same.
    !
    ! !USES
    use abortutils       , only : endrun
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    integer  , intent(in)           :: l
    real(r8) , intent(out), pointer :: values(:)
    !
    ! !LOCAL VARIABLES:
    integer                         :: idx
    integer                         :: nvalues

    ! Check the size of data
    call NumValuesPerLandunit(nvalues)
    if (size(values) /= nvalues) then
       call endrun(msg="ERROR Size of data is incorrect "//errmsg(__FILE__, __LINE__))
    endif

    idx = 0

    idx = idx + 1;                                         values(idx) = real( lun%gridcell (l))
    idx = idx + 1;                                         values(idx) =       lun%wtgcell  (l)
    idx = idx + 1;                                         values(idx) = real( lun%coli     (l))
    idx = idx + 1;                                         values(idx) = real( lun%colf     (l))
    idx = idx + 1;                                         values(idx) = real( lun%ncolumns (l))
    idx = idx + 1;                                         values(idx) = real( lun%pfti     (l))
    idx = idx + 1;                                         values(idx) = real( lun%pftf     (l))
    idx = idx + 1;                                         values(idx) = real( lun%npfts    (l))
    idx = idx + 1;                                         values(idx) = real( lun%itype    (l))

    idx = idx + 1; if (lun%ifspecial (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun%lakpoi    (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun%urbpoi    (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun%glcmecpoi (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun%active    (l))                 values(idx) = 1._r8

    idx = idx + 1; if (.not. isnan(lun%canyon_hwr(l)   )) values(idx) = lun%canyon_hwr   (l)
    idx = idx + 1; if (.not. isnan(lun%wtroad_perv(l)  )) values(idx) = lun%wtroad_perv  (l)
    idx = idx + 1; if (.not. isnan(lun%ht_roof(l)      )) values(idx) = lun%ht_roof      (l)
    idx = idx + 1; if (.not. isnan(lun%wtlunit_roof(l) )) values(idx) = lun%wtlunit_roof (l)
    idx = idx + 1; if (.not. isnan(lun%z_0_town(l)     )) values(idx) = lun%z_0_town     (l)
    idx = idx + 1; if (.not. isnan(lun%z_d_town(l)     )) values(idx) = lun%z_d_town     (l)

  end subroutine GetValuesForLandunit

  !------------------------------------------------------------------------
  subroutine SetValuesForLandunit(l, values)
    !
    ! !DESCRIPTION:
    ! Sets all values for the l-th landunit.
    !
    ! NOTE: The order of the values retrieved/stored in
    !       GetValuesForLandunit/SetValuesForLandunit should be same.
    !
    ! !USES
    use abortutils       , only : endrun
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    integer , intent(in)          :: l
    real(r8), intent(in), pointer :: values(:)
    !
    ! !LOCAL VARIABLES:
    integer                       :: idx
    integer                       :: nvalues

    ! Check the size of data
    call NumValuesPerLandunit(nvalues)
    if (size(values) /= nvalues) then
       call endrun(msg="ERROR Size of data is incorrect "//errmsg(__FILE__, __LINE__))
    endif

    idx = 0

    idx = idx + 1;                           lun%gridcell     (l) = int(values(idx))
    idx = idx + 1;                           lun%wtgcell      (l) =     values(idx)
    idx = idx + 1;                           lun%coli         (l) = int(values(idx))
    idx = idx + 1;                           lun%colf         (l) = int(values(idx))
    idx = idx + 1;                           lun%ncolumns     (l) = int(values(idx))
    idx = idx + 1;                           lun%pfti         (l) = int(values(idx))
    idx = idx + 1;                           lun%pftf         (l) = int(values(idx))
    idx = idx + 1;                           lun%npfts        (l) = int(values(idx))
    idx = idx + 1;                           lun%itype        (l) = int(values(idx))

    idx = idx + 1; if (values(idx) == 1._r8) lun%ifspecial    (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun%lakpoi       (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun%urbpoi       (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun%glcmecpoi    (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun%active       (l) = .true.

    idx = idx + 1;                           lun%canyon_hwr   (l) = values(idx)
    idx = idx + 1;                           lun%wtroad_perv  (l) = values(idx)
    idx = idx + 1;                           lun%ht_roof      (l) = values(idx)
    idx = idx + 1;                           lun%wtlunit_roof (l) = values(idx)
    idx = idx + 1;                           lun%z_0_town     (l) = values(idx)
    idx = idx + 1;                           lun%z_d_town     (l) = values(idx)

  end subroutine SetValuesForLandunit

  !------------------------------------------------------------------------
  subroutine NumValuesPerLandunit(nvalues)
    !
    ! !DESCRIPTION:
    ! Returns the number of values for each landunit within the landunit_type
    !
    ! !ARGUMENTS:
    integer, intent(out)       :: nvalues

    nvalues = 20

  end subroutine NumValuesPerLandunit

end module LandunitMod
