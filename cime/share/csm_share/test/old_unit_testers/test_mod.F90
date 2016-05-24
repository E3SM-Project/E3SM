module test_mod

use shr_kind_mod, only : SHR_KIND_R8
use shr_sys_mod,  only : shr_sys_abort

implicit none

public test_init
public test_is
public test_close
public test_final

integer, save :: ntests       = 0
integer, save :: npass        = 0
integer, save :: num_expected = 0
logical, save :: num_expected_given = .false.
character(*), parameter :: formatTest  = '(A4, " ", i5.5, " - ", A)'
character(*), parameter :: formatArrayMatch = &
    '(" (all ", i5, " values match)")'
character(*), parameter :: formatArray2DMatch = &
    '(" (all ", i5, "x", i5, " values match)")'
character(*), parameter :: formatArrayMisMatch = &
    '(" (only ", i5, " values of ", i5, " values match)")'
character(*), parameter :: formatArray2DMisMatch = &
    '(" (only ", i5, " values of ", i5, "x", i5, " values match)")'
character(*), parameter :: formatRArrayClose   = &
    '(" (all ", i5, " values are within", 1pe9.1e2, " )")'
character(*), parameter :: formatRArrayNotClose = &
    '(" (only ", i5, " values of ", i5, " values are within", 1pe9.1e2, " max diff= ", 1pe9.1e2, ")")'
character(*), parameter :: formatRClose   = &
    '(" ( value within", 1pe9.1e2, " )")'
character(*), parameter :: formatRNotClose = &
    '(" ( value within", 1pe9.1e2, " diff= ", 1pe9.1e2, ")")'
 
interface test_is
  module procedure test_is_logical
  module procedure test_is_logical1D
  module procedure test_is_string
  module procedure test_is_integer
  module procedure test_is_integer1D
  module procedure test_is_real1D
  module procedure test_is_real2D
  module procedure test_is_realScalar
end interface test_is

interface test_close
  module procedure test_close_real1D
  module procedure test_close_realScalar
end interface test_close

private test_is_logical
private test_is_string
private test_is_integer
private test_is_integer1D
private test_is_real1D
private test_is_realScalar
private test_close_real1D

contains


subroutine test_init( num_expected_tests )
   integer, intent(IN), optional :: num_expected_tests

   if ( present(num_expected_tests) ) then
      num_expected = num_expected_tests
      num_expected_given = .true.
      write(*,formatTest) "1...", num_expected, "expected tests"
      write(*,*)
   end if

end subroutine test_init

subroutine test_is_logical( pass, description )

  implicit none

  logical,      intent(IN) :: pass        ! If matches or not
  character(*), intent(IN) :: description ! description of test

  character(4) :: status

  ntests = ntests + 1
  if ( pass )then
    npass = npass + 1
    status = "PASS"
  else
    status = "FAIL"
  end if
  write(*,formatTest) status, ntests, trim(description)

end subroutine test_is_logical

subroutine test_is_logical1D( value, expected, description )

  implicit none

  logical,      intent(IN) :: value(:)    ! test value
  logical,      intent(IN) :: expected(:) ! expected value
  character(*), intent(IN) :: description ! description of test

  logical :: pass
  integer :: nsize, nmatch
  character(256) :: descrip

  nsize = size(value)
  if ( all(value .eqv. expected) )then
     pass = .true.
     write(descrip,formatArrayMatch) nsize
  else
     nmatch = count(value .eqv. expected)
     write(descrip,formatArrayMisMatch) nmatch, nsize
     pass = .false.
  end if
  call test_is_logical( pass, trim(description)//trim(descrip) )

end subroutine test_is_logical1D
   

subroutine test_is_string( value, expected, description )

  implicit none

  character(len=*), intent(IN) :: value
  character(len=*), intent(IN) :: expected
  character(len=*), intent(IN) :: description ! description of test


  logical :: pass        ! If matches or not

  character(4) :: status

  if ( trim(value) == trim(expected) )then
     pass = .true.
  else
     pass = .false.
  end if
  ntests = ntests + 1
  if ( pass )then
    npass = npass + 1
    status = "PASS"
  else
    status = "FAIL"
  end if
  write(*,formatTest) status, ntests, trim(description)

end subroutine test_is_string

subroutine test_is_integer( value, expected, description )
  integer,      intent(IN) :: value       ! test value
  integer,      intent(IN) :: expected    ! expected value
  character(*), intent(IN) :: description ! description of test

  logical :: pass

  if ( value == expected )then
     pass = .true.
  else
     pass = .false.
  end if
  call test_is_logical( pass, description )

end subroutine test_is_integer

subroutine test_is_integer1D( value, expected, description )
  integer,      intent(IN) :: value(:)    ! test value
  integer,      intent(IN) :: expected(:) ! expected value
  character(*), intent(IN) :: description ! description of test

  logical :: pass
  integer :: nsize, nmatch
  character(256) :: descrip

  nsize = size(value)
  if ( all(value == expected) )then
     pass = .true.
     write(descrip,formatArrayMatch) nsize
  else
     nmatch = count(value == expected)
     write(descrip,formatArrayMisMatch) nmatch, nsize
     pass = .false.
  end if
  call test_is_logical( pass, trim(description)//trim(descrip) )

end subroutine test_is_integer1D

subroutine test_is_real1D( value, expected, description )
  real(SHR_KIND_R8), intent(IN) :: value(:)    ! test value
  real(SHR_KIND_R8), intent(IN) :: expected(:) ! expected value
  character(*),      intent(IN) :: description ! description of test

  logical :: pass
  integer :: nsize, nmatch
  character(256) :: descrip

  nsize = size(value)
  if ( all(value == expected) )then
     pass = .true.
     write(descrip,formatArrayMatch) nsize
  else
     nmatch = count(value == expected)
     write(descrip,formatArrayMisMatch) nmatch, nsize
     pass = .false.
  end if
  call test_is_logical( pass, trim(description)//trim(descrip) )

end subroutine test_is_real1D

subroutine test_is_real2D( value, expected, description )
  real(SHR_KIND_R8), intent(IN) :: value(:,:)    ! test value
  real(SHR_KIND_R8), intent(IN) :: expected(:,:) ! expected value
  character(*),      intent(IN) :: description ! description of test

  logical :: pass
  integer :: nsize1, nsize2, nmatch
  character(256) :: descrip

  nsize1 = size(value,1)
  nsize2 = size(value,2)
  if ( all(value == expected) )then
     pass = .true.
     write(descrip,formatArray2DMatch) nsize1, nsize2
  else
     nmatch = count(value == expected)
     write(descrip,formatArray2DMisMatch) nmatch, nsize1, nsize2
     pass = .false.
  end if
  call test_is_logical( pass, trim(description)//trim(descrip) )

end subroutine test_is_real2D

subroutine test_is_realScalar( value, expected, description )
  real(SHR_KIND_R8), intent(IN) :: value       ! test value
  real(SHR_KIND_R8), intent(IN) :: expected    ! expected value
  character(*),      intent(IN) :: description ! description of test

  logical :: pass

  if ( value == expected )then
     pass = .true.
  else
     pass = .false.
  end if
  call test_is_logical( pass, description )

end subroutine test_is_realScalar

subroutine test_close_real1D( value, expected, eps, description, rel_diff )
  real(SHR_KIND_R8), intent(IN) :: value(:)    ! test value
  real(SHR_KIND_R8), intent(IN) :: expected(:) ! expected value
  real(SHR_KIND_R8), intent(IN) :: eps         ! epsilon -- how close to be within
  character(*),      intent(IN) :: description ! description of test
  logical, optional, intent(IN) :: rel_diff    ! if should do relative difference or not

  logical :: pass, lreldiff
  integer :: nsize, nmatch, i, n0(1), nf(1)
  real(SHR_KIND_R8) :: within, diff
  character(256) :: descrip

  lreldiff = .false.
  if ( present(rel_diff) ) lreldiff = rel_diff
  nsize  = size(value)
  if ( nsize /= size(expected) )then
     call shr_sys_abort( "size of value and expected array is different" )
  end if
  if ( any(lbound(value) /= lbound(expected)) )then
     call shr_sys_abort( "lower bound of value and expected array is different" )
  end if
  nmatch = 0
  n0     = lbound(value)
  nf     = ubound(value)
  within = abs(value(n0(1)) - expected(n0(1)))
  if ( lreldiff .and. within > 0.0_SHR_KIND_R8 ) within = within / max( abs(value(n0(1))), abs(expected(n0(1))) )
  do i = n0(1), nf(1)
     diff   = abs(value(i) - expected(i))
     if ( lreldiff .and. diff > 0.0_SHR_KIND_R8 ) diff = diff / max(abs(value(i)),abs(expected(i)) )
     within = max( within, diff )
     if ( diff <= eps ) nmatch = nmatch + 1
  end do
  if( nmatch == nsize )then
     write(descrip,formatRArrayClose) nsize, eps
     pass = .true.
  else
     write(descrip,formatRArrayNotClose) nmatch, nsize, eps, within
     pass = .false.
  end if
  call test_is_logical( pass, trim(description)//trim(descrip) )

end subroutine test_close_real1D

subroutine test_close_realScalar( value, expected, eps, description )
  real(SHR_KIND_R8), intent(IN) :: value       ! test value
  real(SHR_KIND_R8), intent(IN) :: expected    ! expected value
  real(SHR_KIND_R8), intent(IN) :: eps         ! epsilon -- how close to be within
  character(*),      intent(IN) :: description ! description of test

  logical :: pass
  real(SHR_KIND_R8) :: diff
  character(256) :: descrip

  diff   = abs(value - expected)
  if ( diff <= eps ) then
     write(descrip,formatRClose) eps
     pass = .true.
  else
     write(descrip,formatRNotClose) eps, diff
     pass = .false.
  end if
  call test_is_logical( pass, trim(description)//trim(descrip) )

end subroutine test_close_realScalar

subroutine test_final( PassStatus )

  logical, intent(OUT), optional :: PassStatus

  character(4)  :: status
  character(50) :: desc

  write(*,*)
  status = "PASS"
  if ( present(PassStatus) ) PassStatus = .true.
  desc   = "All expected tests ran successfully"
  if ( num_expected_given .and. ntests /= num_expected )then
     status = "FAIL"
     desc   = "Different number of tests than expected"
     if ( present(PassStatus) ) PassStatus = .false.
  end if
  if ( npass  /= ntests       )then
     status = "FAIL"
     if ( present(PassStatus) ) PassStatus = .false.
     write(desc,'(A,i3,A)') "Not all tests passed (", &
                              ntests-npass, " tests failed)"
  end if
  write(*,formatTest) status, ntests, "tests run -- "//desc

end subroutine test_final

end module test_mod
