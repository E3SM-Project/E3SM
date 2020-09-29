module test_mkurbanparMod
! Module for testing mkurbanparMod

   use mkurbanparMod
   use test_mod
   use shr_kind_mod, only : r8 => shr_kind_r8

   implicit none
   private

   public :: test_normalize_urbn_by_tot   

   character(len=*), parameter :: modname = 'test_mkurbanparMod'

contains
   
!------------------------------------------------------------------------------
   subroutine test_normalize_urbn_by_tot

      use mkutilsMod, only : normalize_classes_by_gcell
      
      implicit none

      character(len=128) :: testname

      real(r8), allocatable :: classes_pct_gcell_t(:,:)
      real(r8), allocatable :: classes_pct_gcell(:,:)
      real(r8), allocatable :: classes_pct_tot(:,:)
      real(r8), allocatable :: sums(:)

      integer :: n,nmax,nclass,totsize

      real(r8), parameter :: eps = 1.e-13_r8

      character(len=*), parameter :: subname = 'test_normalize_urbn_by_tot'


      ! This test does a basic check of both normalize_urbn_by_tot and
      ! normalize_classes_by_gcell, by ensuring that when the two are called in
      ! succession, the result is the same as the initial values
      ! (Note that it doesn't directly check the intermediate values -- i.e. the output
      ! produced by normalize_urbn_by_tot)
      testname = 'normalize_urbn_by_tot then normalize_classes_by_gcell'
      nmax = 7
      nclass = 3
      totsize = nmax*nclass
      allocate(classes_pct_gcell_t(nmax,nclass), &
               classes_pct_gcell  (nmax,nclass), &
               classes_pct_tot    (nmax,nclass), &
               sums               (nmax))

      ! The following values are designed to test a number of things, including summing
      ! to 100, summing to 0, some values 0 for a given n, and no values being 0 for a
      ! given n
      classes_pct_gcell_t(:,1) = (/  0.,  5., 0.,   0., 10.,  0., 10./)
      classes_pct_gcell_t(:,2) = (/  0.,  0., 0., 100., 30., 15., 50./)
      classes_pct_gcell_t(:,3) = (/100., 30., 0.,   0., 20.,  0., 40./)

      do n = 1, nmax
         sums(n) = sum(classes_pct_gcell_t(n,:))
      end do

      call normalize_urbn_by_tot(classes_pct_gcell_t, sums, classes_pct_tot)
      call normalize_classes_by_gcell(classes_pct_tot, sums, classes_pct_gcell)
      call test_close(reshape(classes_pct_gcell, (/totsize/)), &
                      reshape(classes_pct_gcell_t, (/totsize/)), &
                      eps, modname//' -- '//subname//' -- '//trim(testname), rel_diff=.true.)

      deallocate(classes_pct_gcell_t, classes_pct_gcell, classes_pct_tot, sums)


   end subroutine test_normalize_urbn_by_tot
!------------------------------------------------------------------------------

end module test_mkurbanparMod
