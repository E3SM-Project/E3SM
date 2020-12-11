module LandunitMod

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varcon     , only : ispval
  use LandunitType   , only : lun_pp
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

    idx = idx + 1;                                         values(idx) = real( lun_pp%gridcell (l))
    idx = idx + 1;                                         values(idx) =       lun_pp%wtgcell  (l)
    idx = idx + 1;                                         values(idx) = real( lun_pp%coli     (l))
    idx = idx + 1;                                         values(idx) = real( lun_pp%colf     (l))
    idx = idx + 1;                                         values(idx) = real( lun_pp%ncolumns (l))
    idx = idx + 1;                                         values(idx) = real( lun_pp%pfti     (l))
    idx = idx + 1;                                         values(idx) = real( lun_pp%pftf     (l))
    idx = idx + 1;                                         values(idx) = real( lun_pp%npfts    (l))
    idx = idx + 1;                                         values(idx) = real( lun_pp%itype    (l))

    idx = idx + 1; if (lun_pp%ifspecial (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun_pp%lakpoi    (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun_pp%urbpoi    (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun_pp%glcmecpoi (l))                 values(idx) = 1._r8
    idx = idx + 1; if (lun_pp%active    (l))                 values(idx) = 1._r8

    idx = idx + 1; if (.not. isnan(lun_pp%canyon_hwr(l)   )) values(idx) = lun_pp%canyon_hwr   (l)
    idx = idx + 1; if (.not. isnan(lun_pp%wtroad_perv(l)  )) values(idx) = lun_pp%wtroad_perv  (l)
    idx = idx + 1; if (.not. isnan(lun_pp%ht_roof(l)      )) values(idx) = lun_pp%ht_roof      (l)
    idx = idx + 1; if (.not. isnan(lun_pp%wtlunit_roof(l) )) values(idx) = lun_pp%wtlunit_roof (l)
    idx = idx + 1; if (.not. isnan(lun_pp%z_0_town(l)     )) values(idx) = lun_pp%z_0_town     (l)
    idx = idx + 1; if (.not. isnan(lun_pp%z_d_town(l)     )) values(idx) = lun_pp%z_d_town     (l)

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

    idx = idx + 1;                           lun_pp%gridcell     (l) = int(values(idx))
    idx = idx + 1;                           lun_pp%wtgcell      (l) =     values(idx)
    idx = idx + 1;                           lun_pp%coli         (l) = int(values(idx))
    idx = idx + 1;                           lun_pp%colf         (l) = int(values(idx))
    idx = idx + 1;                           lun_pp%ncolumns     (l) = int(values(idx))
    idx = idx + 1;                           lun_pp%pfti         (l) = int(values(idx))
    idx = idx + 1;                           lun_pp%pftf         (l) = int(values(idx))
    idx = idx + 1;                           lun_pp%npfts        (l) = int(values(idx))
    idx = idx + 1;                           lun_pp%itype        (l) = int(values(idx))

    idx = idx + 1; if (values(idx) == 1._r8) lun_pp%ifspecial    (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun_pp%lakpoi       (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun_pp%urbpoi       (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun_pp%glcmecpoi    (l) = .true.
    idx = idx + 1; if (values(idx) == 1._r8) lun_pp%active       (l) = .true.

    idx = idx + 1;                           lun_pp%canyon_hwr   (l) = values(idx)
    idx = idx + 1;                           lun_pp%wtroad_perv  (l) = values(idx)
    idx = idx + 1;                           lun_pp%ht_roof      (l) = values(idx)
    idx = idx + 1;                           lun_pp%wtlunit_roof (l) = values(idx)
    idx = idx + 1;                           lun_pp%z_0_town     (l) = values(idx)
    idx = idx + 1;                           lun_pp%z_d_town     (l) = values(idx)

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
