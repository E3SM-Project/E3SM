module test_mkutilsMod
! Module for testing mkutilsMod

  use mkutilsMod
  use test_mod
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: test_remove_small_cover
  public :: test_slightly_below
  public :: test_slightly_above

  character(len=*), parameter :: modname = 'test_mkutilsMod'

contains

!------------------------------------------------------------------------------
  subroutine test_remove_small_cover
    
    implicit none

    character(len=128) :: testname
    real(r8) :: pct_gcell
    real(r8) :: pct_gcell_t
    real(r8) :: pct_lunit(6:15)   ! length-10 array, but offset to make sure the routine doesn't break in the presence of non-1 lower bound
    real(r8) :: pct_lunit_t(6:15)
    real(r8), allocatable :: pct_lunit_empty(:)
    real(r8), allocatable :: pct_lunit_empty_t(:)
    integer  :: nsmall
    integer  :: nsmall_t

    real(r8), parameter :: eps = 1.e-16
    character(len=*), parameter :: subname = 'test_remove_small_cover'

    testname = 'zero cover'
    pct_gcell = 0._r8
    pct_lunit(:) = 0._r8
    pct_lunit(7) = 40._r8
    pct_lunit(9) = 60._r8
    call remove_small_cover(pct_gcell, pct_lunit, nsmall, too_small=1.e-12_r8)
    pct_gcell_t = 0._r8
    pct_lunit_t = pct_lunit
    nsmall_t = 0
    call check_results

    testname = 'only one pft > 0 and it is small'
    pct_gcell = 1.e-13_r8
    pct_lunit(:) = 0._r8
    pct_lunit(7) = 100._r8
    call remove_small_cover(pct_gcell, pct_lunit, nsmall, too_small=1.e-12_r8)
    pct_gcell_t = 0._r8
    pct_lunit_t = pct_lunit
    nsmall_t = 1
    call check_results

    testname = 'ten PFTs, 5 small, 5 zero'
    pct_gcell = 1.e-12_r8
    pct_lunit(:) = (/0._r8, 20._r8, 0._r8, 20._r8, 0._r8, 20._r8, 0._r8, 20._r8, 0._r8, 20._r8/)
    call remove_small_cover(pct_gcell, pct_lunit, nsmall, too_small=1.e-12_r8)
    pct_gcell_t = 0._r8
    pct_lunit_t = pct_lunit
    nsmall_t = 5
    call check_results

    testname = 'ten PFTs, all small'
    pct_gcell = 5.e-12_r8  ! by itself, this is greater than the too_small threshold, but it will be divided by 10
    pct_lunit(:) = (/10._r8, 10._r8, 10._r8, 10._r8, 10._r8, 10._r8, 10._r8, 10._r8, 10._r8, 10._r8/)
    call remove_small_cover(pct_gcell, pct_lunit, nsmall, too_small=1.e-12_r8)
    pct_gcell_t = 0._r8
    pct_lunit_t = pct_lunit
    nsmall_t = 10
    call check_results

    testname = 'ten PFTs, all > 0, 5 small, 5 big'
    pct_gcell = 1._r8
    pct_lunit(6:10)  = 20._r8 - 1.e-11_r8
    pct_lunit(11:15) = 1.e-11_r8  ! so 1.e-13 on the grid cell
    call remove_small_cover(pct_gcell, pct_lunit, nsmall, too_small=1.e-12_r8)
    pct_gcell_t = 1._r8 - 5*(1._r8 * 1.e-11_r8 / 100._r8)
    pct_lunit_t = (/20._r8, 20._r8, 20._r8, 20._r8, 20._r8, 0._r8, 0._r8, 0._r8, 0._r8, 0._r8/)
    nsmall_t = 5
    call check_results

    testname = 'ten PFTs, all either 0 or big'
    pct_gcell = 1._r8
    pct_lunit(6:10)  = 20._r8
    pct_lunit(11:15) = 0._r8
    call remove_small_cover(pct_gcell, pct_lunit, nsmall, too_small=1.e-12_r8)
    pct_gcell_t = pct_gcell
    pct_lunit_t = pct_lunit
    nsmall_t = 0
    call check_results

    testname = 'zero-length PFT array'
    pct_gcell = 0._r8
    allocate(pct_lunit_empty(0))
    allocate(pct_lunit_empty_t(0))
    call remove_small_cover(pct_gcell, pct_lunit_empty, nsmall, too_small=1.e-12_r8)
    pct_gcell_t = pct_gcell
    nsmall_t = 0
    call check_results

  contains
    subroutine check_results
      call test_close(pct_gcell, pct_gcell_t, eps, modname//' -- '//subname//' -- '//&
           trim(testname)//' -- pct_gcell')
      call test_close(pct_lunit, pct_lunit_t, eps, modname//' -- '//subname//' -- '//&
           trim(testname)//' -- pct_lunit')
      call test_is(nsmall, nsmall_t, modname//' -- '//subname//' -- '//&
           trim(testname)//' -- nsmall')
    end subroutine check_results

  end subroutine test_remove_small_cover
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
  subroutine test_slightly_below

    implicit none

    character(len=128) :: testname

    logical :: retval
    real(r8) :: a
    real(r8) :: b

    character(len=*), parameter :: subname = 'test_slightly_below'

    testname='basic-true'
    b = 3.0
    a = 3.0 - b*epsilon(b)
    retval = slightly_below(a,b)
    call test_is((retval .eqv. .true.), modname//' -- '//subname//' -- '//trim(testname))

    testname='far below'
    b = 3.0
    a = 2.0
    retval = slightly_below(a,b)
    call test_is((retval .eqv. .false.), modname//' -- '//subname//' -- '//trim(testname))
    
    testname='equal'
    b = 3.0
    a = 3.0
    retval = slightly_below(a,b)
    call test_is((retval .eqv. .false.), modname//' -- '//subname//' -- '//trim(testname))

    testname='above'
    b = 3.0
    a = 3.0 + epsilon(b)
    retval = slightly_below(a,b)
    call test_is((retval .eqv. .false.), modname//' -- '//subname//' -- '//trim(testname))

    testname='change epsilon to allow far below'
    b = 3.0
    a = 2.0
    retval = slightly_below(a,b,eps=0.75_r8)
    call test_is((retval .eqv. .true.), modname//' -- '//subname//' -- '//trim(testname))

  end subroutine test_slightly_below
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine test_slightly_above

    implicit none

    character(len=128) :: testname

    logical :: retval
    real(r8) :: a
    real(r8) :: b

    character(len=*), parameter :: subname = 'test_slightly_above'

    testname='basic-true'
    b = 3.0
    a = 3.0 + b*epsilon(b)
    retval = slightly_above(a,b)
    call test_is((retval .eqv. .true.), modname//' -- '//subname//' -- '//trim(testname))

    testname='far above'
    b = 3.0
    a = 4.0
    retval = slightly_above(a,b)
    call test_is((retval .eqv. .false.), modname//' -- '//subname//' -- '//trim(testname))
    
    testname='equal'
    b = 3.0
    a = 3.0
    retval = slightly_above(a,b)
    call test_is((retval .eqv. .false.), modname//' -- '//subname//' -- '//trim(testname))

    testname='below'
    b = 3.0
    a = 3.0 - epsilon(b)
    retval = slightly_above(a,b)
    call test_is((retval .eqv. .false.), modname//' -- '//subname//' -- '//trim(testname))

    testname='change epsilon to allow far above'
    b = 3.0
    a = 4.0
    retval = slightly_above(a,b,eps=0.75_r8)
    call test_is((retval .eqv. .true.), modname//' -- '//subname//' -- '//trim(testname))

  end subroutine test_slightly_above
!------------------------------------------------------------------------------

end module test_mkutilsMod


