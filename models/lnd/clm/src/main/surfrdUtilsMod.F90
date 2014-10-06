module surfrdUtilsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains utility methods that can be used when reading surface datasets or similar
  ! datasets (such as the landuse_timeseries dataset)
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use clm_varctl   , only : iulog
  use abortutils   , only : endrun
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: check_sums_equal_1  ! Confirm that sum(arr(n,:)) == 1 for all n
  
  !-----------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  subroutine check_sums_equal_1(arr, lb, name, caller)
    !
    ! !DESCRIPTION:
    ! Confirm that sum(arr(n,:)) == 1 for all n. If this isn't true for any n, abort with a message.
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: lb           ! lower bound of the first dimension of arr
    real(r8)        , intent(in) :: arr(lb:,:)   ! array to check
    character(len=*), intent(in) :: name         ! name of array
    character(len=*), intent(in) :: caller       ! identifier of caller, for more meaningful error messages
    !
    ! !LOCAL VARIABLES:
    logical :: found
    integer :: nl
    integer :: nindx
    real(r8), parameter :: eps = 1.e-14_r8
    !-----------------------------------------------------------------------

    found = .false.

    do nl = lbound(arr, 1), ubound(arr, 1)
       if (abs(sum(arr(nl,:)) - 1._r8) > eps) then
          found = .true.
          nindx = nl
          exit
       end if
    end do

    if (found) then
       write(iulog,*) trim(caller), ' ERROR: sum of ', trim(name), ' not 1.0 at nl=', nindx
       write(iulog,*) 'sum is: ', sum(arr(nindx,:))
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine check_sums_equal_1

end module surfrdUtilsMod
