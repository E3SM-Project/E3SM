module test_mkutilsMod
! Module for testing mkutilsMod

  use mkutilsMod
  use test_mod
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: test_slightly_below
  public :: test_slightly_above

  character(len=*), parameter :: modname = 'test_mkutilsMod'

contains
  
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


