! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_test_compute_interval
!
!> \brief MPAS ocean analysis core member: test_compute_interval
!> \author Mark Petersen
!> \date   May 2015
!> \details
!>  MPAS ocean analysis core member: test_compute_interval
!>  
!
!-----------------------------------------------------------------------

module ocn_test_compute_interval

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_setup_packages_test_compute_interval, &
             ocn_init_test_compute_interval, &
             ocn_compute_test_compute_interval, &
             ocn_restart_test_compute_interval, &
             ocn_finalize_test_compute_interval

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: am_test_compute_intervalTimer

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_setup_packages_test_compute_interval
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_test_compute_interval(configPool, packagePool, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      logical, pointer :: am_test_compute_intervalActive

      err = 0

      call mpas_pool_get_package(packagePool, 'am_test_compute_intervalActive', am_test_compute_intervalActive)

      ! turn on package for this analysis member
      am_test_compute_intervalActive = .true.

   end subroutine ocn_setup_packages_test_compute_interval!}}}


!***********************************************************************
!
!  routine ocn_init_test_compute_interval
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_test_compute_interval(domain, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: am_test_compute_intervalPool
      type (block_type), pointer :: block
      real (kind=RKIND), pointer :: test_compute_interval_counter

      err = 0

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'am_test_compute_interval', am_test_compute_intervalPool)
         call mpas_pool_get_array(am_test_compute_intervalPool, 'test_compute_interval_counter',test_compute_interval_counter)

         test_compute_interval_counter = 0

         block => block % next
      end do

   end subroutine ocn_init_test_compute_interval!}}}

!***********************************************************************
!
!  routine ocn_compute_test_compute_interval
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2015
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_test_compute_interval(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: am_test_compute_intervalPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: am_test_compute_interval
      character (len=StrKIND), pointer :: xtime

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, num_tracers
      integer :: iTracer, k, iCell
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real (kind=RKIND), pointer :: test_compute_interval_counter

      err = 0

      dminfo = domain % dminfo

      call mpas_timer_start("compute_test_compute_interval", .false., am_test_compute_intervalTimer)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'am_test_compute_interval', am_test_compute_intervalPool)
         call mpas_pool_get_array(am_test_compute_intervalPool, 'test_compute_interval_counter',test_compute_interval_counter)

         test_compute_interval_counter = test_compute_interval_counter + 1
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         block => block % next
      end do

      write(6,'(3A,F10.1)') 'in ocn_compute_test_compute_interval, time = ', trim(xtime), &
         ',   test_compute_interval_counter = ',test_compute_interval_counter

      call mpas_timer_stop("test_compute_interval", am_test_compute_intervalTimer)

   end subroutine ocn_compute_test_compute_interval!}}}

!***********************************************************************
!
!  routine ocn_restart_test_compute_interval
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2015
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_test_compute_interval(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_test_compute_interval!}}}

!***********************************************************************
!
!  routine ocn_finalize_test_compute_interval
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    May 2015
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_test_compute_interval(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_test_compute_interval!}}}

end module ocn_test_compute_interval

! vim: foldmethod=marker
