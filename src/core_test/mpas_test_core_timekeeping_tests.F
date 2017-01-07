! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module test_core_timekeeping_tests

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_field_routines
   use mpas_dmpar
   use mpas_threading
   use mpas_io_units
   use mpas_timer
   use mpas_timekeeping

   implicit none
   private

   public :: test_core_test_intervals

   contains

   !***********************************************************************
   !
   !  routine test_core_test_intervals
   !
   !> \brief   MPAS Test Core timekeeping interval tests
   !> \author  Doug Jacobsen
   !> \date    08/23/2016
   !> \details 
   !>  This routine performs tests related to timekeeping intervals
   !
   !-----------------------------------------------------------------------
   subroutine test_core_test_intervals(domain, threadErrs, err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, dimension(:), intent(out) :: threadErrs
      integer, intent(out) :: err

      character (len=StrKIND) :: int1_str, int2_str
      character (len=StrKIND) :: ref_str

      integer :: threadNum
      integer :: iErr, err_tmp
      integer :: divs

      call mpas_timer_start('timekeeping tests')

      threadNum = mpas_threading_get_thread_num()
      err = 0

      if ( threadNum == 0 ) then
         write(stderrUnit, *) ' Performing time interval tests'

         write(stderrUnit, *) '  Test 1:'
         call test_core_interval_test('0001-01-01_00:00:00', '0000-01-00_10:00:00', '0001_00:00:00', 31, '0000_10:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            write(stderrUnit, *) '   Result: PASSED'
         else
            write(stderrUnit, *) ' * Result: FAILED'
         end if

         write(stderrUnit, *) '  Test 2:'
         call test_core_interval_test('0001-01-01_00:00:00', '0000-01-00_00:00:00', '0001_00:00:00', 31, '0000_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            write(stderrUnit, *) '   Result: PASSED'
         else
            write(stderrUnit, *) ' * Result: FAILED'
         end if

         write(stderrUnit, *) '  Test 3:'
         call test_core_interval_test('0001-02-01_00:00:00', '0000-01-00_10:00:00', '0001_00:00:00', 28, '0000_10:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            write(stderrUnit, *) '   Result: PASSED'
         else
            write(stderrUnit, *) ' * Result: FAILED'
         end if

         write(stderrUnit, *) '  Test 4:'
         call test_core_interval_test('0001-02-01_00:00:00', '0000-01-00_00:00:00', '0001_00:00:00', 28, '0000_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            write(stderrUnit, *) '   Result: PASSED'
         else
            write(stderrUnit, *) ' * Result: FAILED'
         end if

         write(stderrUnit, *) '  Test 5:'
         call test_core_interval_test('0001-01-01_00:00:00', '0000-00-00_01:00:00', '0000_00:30:00', 2, '0000_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            write(stderrUnit, *) '   Result: PASSED'
         else
            write(stderrUnit, *) ' * Result: FAILED'
         end if

         write(stderrUnit, *) '  Test 6:'
         call test_core_interval_test('0001-01-01_00:00:00', '0001-01-00_00:00:00', '0001-00-00_00:00:00', 1, '0000-00-31_00:00:00', err_tmp)
         if ( err_tmp == 0 ) then
            write(stderrUnit, *) '   Result: PASSED'
         else
            write(stderrUnit, *) ' * Result: FAILED'
         end if



         write(stderrUnit, *) ' Completed time interval tests'
      end if

      call mpas_timer_stop('timekeeping tests')

   end subroutine test_core_test_intervals!}}}

   subroutine test_core_interval_test(ref_str, int1_str, int2_str, expected_divs, expected_remainder_str, ierr)!{{{
      character (len=*), intent(in) :: ref_str, int1_str, int2_str
      integer, intent(in) :: expected_divs
      character (len=*), intent(in) :: expected_remainder_str
      integer, intent(out) :: ierr

      integer :: divs

      character (len=StrKIND) :: remainder_str

      type (mpas_time_type) :: ref_time
      type (mpas_timeinterval_type) :: int1, int2, remainder
      type (mpas_timeinterval_type) :: expected_remainder

      integer :: err_tmp

      ierr = 0

      write(stderrUnit, *) ''
      write(stderrUnit, *) '  Testing time intervals:'
      write(stderrUnit, *) '      Reference time: ', trim(ref_str)
      write(stderrUnit, *) '      Interval 1: ', trim(int1_str)
      write(stderrUnit, *) '      Interval 2: ', trim(int2_str)

      call mpas_set_time(ref_time, dateTimeString=ref_str, ierr=err_tmp)
      call mpas_set_timeinterval(int1, timeString=int1_str, ierr=err_tmp)
      call mpas_set_timeinterval(int2, timeString=int2_str, ierr=err_tmp)
      call mpas_set_timeinterval(expected_remainder, timeString=expected_remainder_str, ierr=err_tmp)

      write(stderrUnit, *) '      -- Calling interval division'

      call mpas_interval_division(ref_time, int1, int2, divs, remainder)

      call mpas_get_timeinterval(remainder, startTimeIn=ref_time, timeString=remainder_str)

      write(stderrUnit, *) '      Interval Division summary'
      write(stderrUnit, *) '          Divisions: ', divs
      write(stderrUnit, *) '          Remainder: ', trim(remainder_str)
      write(stderrUnit, *) ''

      if ( divs == expected_divs ) then
         write(stderrUnit, *) '          Div Test: PASSED'
      else
         write(stderrUnit, *) ' **       Div Test: FAILED'
         ierr = 1
      end if

      if ( remainder == expected_remainder ) then
         write(stderrUnit, *) '          Remainder Test: PASSED'
      else
         write(stderrUnit, *) ' **       Remainder Test: FAILED'
         ierr = 1
      end if


   end subroutine test_core_interval_test!}}}

end module test_core_timekeeping_tests
