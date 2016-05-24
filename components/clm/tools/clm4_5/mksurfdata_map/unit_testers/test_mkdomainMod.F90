module test_mkdomainMod
! Module for testing mkindexmapMod

  use mkdomainMod
  use test_mod
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: test_domain_read_dims

  character(len=*), parameter :: modname = 'test_mkdomainMod'
   
contains

!------------------------------------------------------------------------------
  subroutine test_domain_read_dims

    use mkncdio
    
    implicit none
    
    type(domain_type)  :: domain
    integer            :: ncid
    character(len=128) :: testname

    integer :: ni_t, nj_t, ns_t
    logical :: is_2d_t

    character(len=*), parameter :: subname = 'test_domain_read_dims'

    testname = 'lon'
    call check_ret(nf_open('unit_testers/inputs/test_domain_read_dims__lon.nc', 0, ncid), subname)
    ni_t = 2
    nj_t = 3
    ns_t = 6
    is_2d_t = .true.
    call domain_read_dims(domain, ncid)
    call check_results_2d

    testname = 'lsmlon'
    call check_ret(nf_open('unit_testers/inputs/test_domain_read_dims__lsmlon.nc', 0, ncid), subname)
    ni_t = 3
    nj_t = 4
    ns_t = 12
    is_2d_t = .true.
    call domain_read_dims(domain, ncid)
    call check_results_2d

    ! When we have both 'lon' and 'ni', should use 'ni'
    testname = 'lon_and_ni'
    call check_ret(nf_open('unit_testers/inputs/test_domain_read_dims__lon_and_ni.nc', 0, ncid), subname)
    ni_t = 4
    nj_t = 5
    ns_t = 20
    is_2d_t = .true.
    call domain_read_dims(domain, ncid)
    call check_results_2d

    ! test 1-d
    testname = 'num_pixels'
    call check_ret(nf_open('unit_testers/inputs/test_domain_read_dims__num_pixels.nc', 0, ncid), subname)
    ns_t = 17
    is_2d_t = .false.
    call domain_read_dims(domain, ncid)
    call check_results_1d

    ! When we have both 2-d and 1-d info, should use 2-d info
    testname = 'lon_and_num_pixels'
    call check_ret(nf_open('unit_testers/inputs/test_domain_read_dims__lon_and_num_pixels.nc', 0, ncid), subname)
    ni_t = 2
    nj_t = 3
    ns_t = 6
    is_2d_t = .true.
    call domain_read_dims(domain, ncid)
    call check_results_2d

  contains
    subroutine check_results_1d
      call test_is(domain%ns, ns_t, modname//' -- '//subname//' -- '//trim(testname)//' -- ns')
      call test_is((domain%is_2d == is_2d_t), modname//' -- '//subname//' -- '//trim(testname)//' -- is_2d')
    end subroutine check_results_1d

    subroutine check_results_2d
      call test_is(domain%ns, ns_t, modname//' -- '//subname//' -- '//trim(testname)//' -- ns')
      call test_is(domain%ni, ni_t, modname//' -- '//subname//' -- '//trim(testname)//' -- ni')
      call test_is(domain%nj, nj_t, modname//' -- '//subname//' -- '//trim(testname)//' -- nj')
      call test_is((domain%is_2d == is_2d_t), modname//' -- '//subname//' -- '//trim(testname)//' -- is_2d')
    end subroutine check_results_2d
  end subroutine test_domain_read_dims
end module test_mkdomainMod
    

    
