module test_mkchecksMod
! Module for testing mkchecksMod

  use mkchecksMod
  use test_mod
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: test_min_bad
  public :: test_max_bad

  character(len=*), parameter :: modname = 'test_mkchecksMod'

contains

!------------------------------------------------------------------------------
  subroutine test_min_bad

    implicit none

    character(len=128) :: testname
    logical :: test_result

    character(len=*), parameter :: subname = 'test_min_bad'

    ! Tests for r8

    testname = 'r8 - pass'
    test_result = min_bad((/1._r8,2._r8,3._r8/), 0._r8, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    testname = 'r8 - pass on border'
    test_result = min_bad((/1._r8,2._r8,3._r8/), 1._r8, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    ! Note that we expect output to stdout from the following test that indicates an error
    testname = 'r8 - fail'
    test_result = min_bad((/1._r8,2._r8,3._r8/), 1.5_r8, 'testvar')
    call test_is(test_result==.true., modname//' -- '//subname//' -- '//trim(testname))

    ! Tests for int

    testname = 'int - pass'
    test_result = min_bad((/1,2,3/), 0, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    testname = 'int - pass on border'
    test_result = min_bad((/1,2,3/), 1, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    ! Note that we expect output to stdout from the following test that indicates an error
    testname = 'int - fail'
    test_result = min_bad((/1,2,3/), 2, 'testvar')
    call test_is(test_result==.true., modname//' -- '//subname//' -- '//trim(testname))

  end subroutine test_min_bad

!------------------------------------------------------------------------------
  subroutine test_max_bad

    implicit none

    character(len=128) :: testname
    logical :: test_result

    character(len=*), parameter :: subname = 'test_max_bad'

    ! Tests for r8

    testname = 'r8 - pass'
    test_result = max_bad((/1._r8,2._r8,3._r8/), 4._r8, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    testname = 'r8 - pass on border'
    test_result = max_bad((/1._r8,2._r8,3._r8/), 3._r8, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    ! Note that we expect output to stdout from the following test that indicates an error
    testname = 'r8 - fail'
    test_result = max_bad((/1._r8,2._r8,3._r8/), 2.5_r8, 'testvar')
    call test_is(test_result==.true., modname//' -- '//subname//' -- '//trim(testname))

    ! Tests for int

    testname = 'int - pass'
    test_result = max_bad((/1,2,3/), 4, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    testname = 'int - pass on border'
    test_result = max_bad((/1,2,3/), 3, 'testvar')
    call test_is(test_result==.false., modname//' -- '//subname//' -- '//trim(testname))

    ! Note that we expect output to stdout from the following test that indicates an error
    testname = 'int - fail'
    test_result = max_bad((/1,2,3/), 2, 'testvar')
    call test_is(test_result==.true., modname//' -- '//subname//' -- '//trim(testname))

  end subroutine test_max_bad
end module test_mkchecksMod
