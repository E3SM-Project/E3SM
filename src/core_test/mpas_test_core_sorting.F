! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file distributed with this code, or at
! http://mpas-dev.github.com/license.html .
!
module test_core_sorting

   use mpas_derived_types
   use mpas_log

   private

   public :: test_core_test_sorting

   contains

   !***********************************************************************
   !
   !  routine test_core_test_sorting
   !
   !> \brief   Tests performance of sorting routines in framework
   !> \author  Michael Duda
   !> \date    1 September 2016
   !> \details 
   !>  This routine tests the performance of the sorting routines in
   !>  the MPAS framework on various input patterns. Timing information for
   !>  each pattern is written to stderrUnit.
   !
   !-----------------------------------------------------------------------
   subroutine test_core_test_sorting(domain, err)

      use mpas_sort, only : mpas_quicksort
      use mpas_timer, only : mpas_timer_start, mpas_timer_stop

      implicit none

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      integer :: i, l, r, pivot
      real (kind=RKIND) :: temp
      integer, parameter :: n = 122880
      real (kind=RKIND), dimension(n) :: vals
      integer :: count_start, count_stop, count_rate


      call mpas_log_write('Sorting tests:')
      
      !
      ! Random values
      !
      do i=1,n
          call random_number(vals(i))
      end do
      call system_clock(count=count_start)
      call mpas_timer_start('sorting: random')
      call mpas_quicksort(n, vals)
      call mpas_timer_stop('sorting: random')
      call system_clock(count=count_stop)
      call system_clock(count_rate=count_rate)
      call mpas_log_write('    random input timing         (quicksort): $r', realArgs=(/real(count_stop - count_start) / real(count_rate)/))
      
      !
      ! Values in sorted, ascending order
      !
      do i=1,n
          vals(i) = real(i)
      end do
      call system_clock(count=count_start)
      call mpas_timer_start('sorting: sorted')
      call mpas_quicksort(n, vals)
      call mpas_timer_stop('sorting: sorted')
      call system_clock(count=count_stop)
      call system_clock(count_rate=count_rate)
      call mpas_log_write('    sorted input timing         (quicksort): $r', realArgs=(/real(count_stop - count_start) / real(count_rate)/))
      
      !
      ! Values in sorted, descending order
      !
      do i=1,n
          vals(i) = real(n-i+1)
      end do
      call system_clock(count=count_start)
      call mpas_timer_start('sorting: reverse sorted')
      call mpas_quicksort(n, vals)
      call mpas_timer_stop('sorting: reverse sorted')
      call system_clock(count=count_stop)
      call system_clock(count_rate=count_rate)
      call mpas_log_write('    reverse sorted input timing (quicksort): $r', realArgs=(/real(count_stop - count_start) / real(count_rate)/))
      
      !
      ! Constant values
      !
      vals(:) = 42.0
      call system_clock(count=count_start)
      call mpas_timer_start('sorting: constant')
      call mpas_quicksort(n, vals)
      call mpas_timer_stop('sorting: constant')
      call system_clock(count=count_stop)
      call system_clock(count_rate=count_rate)
      call mpas_log_write('    constant input timing       (quicksort): $r', realArgs=(/real(count_stop - count_start) / real(count_rate)/))

      !
      ! Construct theoretically worst-case input for quicksort based on current
      ! method of chosing pivot element
      !
      do i=1,n
         vals(i) = real(i)
      end do
      r = n
      do l=n-1,1,-1
         ! Swap l and r
         temp = vals(l)
         vals(l) = vals(r)
         vals(r) = temp
      
         pivot = (l+r)/2
      
         ! Swap pivot and r
         temp = vals(pivot)
         vals(pivot) = vals(r)
         vals(r) = temp
      end do
      
      call system_clock(count=count_start)
      call mpas_timer_start('sorting: worst-case')
      call mpas_quicksort(n, vals)
      call mpas_timer_stop('sorting: worst-case')
      call system_clock(count=count_stop)
      call system_clock(count_rate=count_rate)
      call mpas_log_write('    ''worst-case'' input timing   (quicksort): $r', realArgs=(/real(count_stop - count_start) / real(count_rate)/))
      
      err = 0

   end subroutine test_core_test_sorting

end module test_core_sorting
