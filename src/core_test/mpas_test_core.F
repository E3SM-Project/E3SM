! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module test_core

   use mpas_framework
   use mpas_timekeeping

   type (MPAS_Clock_type), pointer :: clock

   contains

   !***********************************************************************
   !
   !  function test_core_init
   !
   !> \brief   MPAS Core Initialization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs the necessary initialization of a core. This can
   !>   involve things like setting up coeffiecients for remapping and advection,
   !>   but more generally can include anything the core needs to initialize.
   !
   !-----------------------------------------------------------------------
   function test_core_init(domain, startTimeStamp) result(iErr)!{{{
   
      use mpas_configure
      use mpas_derived_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      type (MPAS_Time_Type) :: startTime
      integer :: iErr


      iErr = 0

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, iErr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

   end function test_core_init!}}}


   !***********************************************************************
   !
   !  function test_core_run
   !
   !> \brief   MPAS Core Run Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs a run of the MPAS core. This can involve time
   !>  stepping if the core requires time stepping, but more generally includes
   !>  anything a core would typically do after it was initialized.
   !
   !-----------------------------------------------------------------------
   function test_core_run(domain) result(iErr)!{{{
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
      use test_core_sorting, only : test_core_test_sorting
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: iErr

      type (mpas_pool_type), pointer :: pool
      type (mpas_pool_iterator_type) :: itr

      !
      ! Test performance of framework sorting routines
      !
      call test_core_test_sorting(domain, iErr)
      if (iErr == 0) then
         write(stderrUnit,*) ' * Sorting tests: SUCCESS'
      else
         write(stderrUnit,*) ' * Sorting tests: FAILURE'
      end if

      iErr = 0

   end function test_core_run!}}}

   !***********************************************************************
   !
   !  function test_core_finalize
   !
   !> \brief   MPAS Core Finalization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function finalizes the MPAS core. It should at a minimum destroy the
   !>  simulation clock, but can perform any functions necessary to clean up the
   !>  MPAS core.
   !
   !-----------------------------------------------------------------------
   function test_core_finalize(domain) result(iErr)!{{{
   
      use mpas_derived_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: iErr
 

      iErr = 0

      call mpas_destroy_clock(clock, iErr)

   end function test_core_finalize!}}}

end module test_core
