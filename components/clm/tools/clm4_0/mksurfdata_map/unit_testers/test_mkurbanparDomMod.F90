module test_mkurbanparDomMod
! Module for testing mkurbanparDomMod

   use mkurbanparDomMod
   use test_mod
   use shr_kind_mod, only : r8 => shr_kind_r8

   implicit none
   private

   public :: test_mkurban_dominant_density   

   character(len=*), parameter :: modname = 'test_mkurbanparDomMod'

contains
   
!------------------------------------------------------------------------------
   subroutine test_mkurban_dominant_density
      
      implicit none

      real(r8), allocatable :: urbn_by_dens_o(:,:)
      integer , allocatable :: dens_o(:)
      integer , allocatable :: dens_o_t(:)
      real(r8), allocatable :: urbn_o(:)
      real(r8), allocatable :: urbn_o_t(:)

      character(len=128) :: testname
      integer :: nodata

      character(len=*), parameter :: subname = 'test_mkurban_dominant_density'

      allocate(urbn_by_dens_o (4, 2), &
               dens_o   (4), &
               dens_o_t (4), &
               urbn_o   (4), &
               urbn_o_t (4))


      testname = 'basic test'
      nodata = -1
      ! This tests a few different things:
      ! (1) output 1 should end up as nodata
      ! (2) output 2 tests a tie
      ! (3) output 3 & 4 test "normal" cases

      urbn_by_dens_o = reshape(source=(/0, 10, 30, 10, &     ! column 1 (i.e., density class 1)
                                        0, 10, 10, 30/), &   ! column 2 (i.e., density class 2)
                               shape=(/4, 2/))
      
      dens_o_t = (/nodata, 1, 1, 2/)
      urbn_o_t = (/0, 20, 40, 40/)

      call mkurban_dominant_density(urbn_by_dens_o, nodata, dens_o, urbn_o)
      call check_results

      deallocate(urbn_by_dens_o, dens_o, dens_o_t, urbn_o, urbn_o_t)

   contains
      subroutine check_results
         call test_is(dens_o, dens_o_t, modname//' -- '//subname//' -- '//trim(testname)//&
                      ' -- dens_o')
         call test_is(urbn_o, urbn_o_t, modname//' -- '//subname//' -- '//trim(testname)//&
                      ' -- urbn_o')
      end subroutine check_results

   end subroutine test_mkurban_dominant_density
!------------------------------------------------------------------------------

end module test_mkurbanparDomMod
