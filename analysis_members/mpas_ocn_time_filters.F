! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_filters
!
!> \brief MPAS ocean analysis mode member: time_filters
!> \author Phillip J. Wolfram
!> \date   07/17/2015
!> \details
!>  Performs time high and low pass filtering.
!>
!-----------------------------------------------------------------------

module ocn_time_filters

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines
#ifdef MPAS_DEBUG
   use mpas_constants
#endif

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

   public :: ocn_init_time_filters, &
             ocn_compute_time_filters, &
             ocn_restart_time_filters, &
             ocn_finalize_time_filters

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
#ifdef MPAS_DEBUG
   integer :: iEdgeOutput = 0, iBlockOutput = 0, iklevel = 1
   real (kind=RKIND) :: lonEdgePoint = (360.0_RKIND-7.5_RKIND)*pii/180.0_RKIND, latEdgePoint = 32.5_RKIND*pii/180.0_RKIND
#endif

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_time_filters
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    07/17/2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_time_filters(domain, err)!{{{

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
      type (block_type), pointer :: block
      logical, pointer :: initializeFilters
      type (mpas_pool_type), pointer :: timeFiltersAMPool, statePool
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalVelocityLowPass, normalVelocityHighPass
      logical, pointer :: config_AM_timeFilters_do_restart

#ifdef MPAS_DEBUG
      real (kind=RKIND), dimension(:), pointer :: lonEdge, latEdge
      real (kind=RKIND) :: dist, distmax = 1e9
      integer :: i, iBlock
      integer, pointer :: nEdgesSolve
#endif

      err = 0

      ! read in data on restart
      call mpas_pool_get_config(domain % configs, 'config_AM_timeFilters_do_restart', config_AM_timeFilters_do_restart)
      if ( config_AM_timeFilters_do_restart ) then
        call MPAS_stream_mgr_read(domain % streamManager, streamID='timeFiltersRestart', ierr=err)
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_AM_timeFilters_initialize_filters', initializeFilters)
      if (initializeFilters) then
#ifdef MPAS_DEBUG
        write(stderrUnit,*) 'initializing time filters'
#endif

        ! loop over all blocks and make assignments
        block => domain % blocklist
        do while (associated(block))

          ! get high and low pass velocity components
          call mpas_pool_get_subpool(block % structs, 'timeFiltersAM', timeFiltersAMPool)
          call mpas_pool_get_array(timeFiltersAMPool, 'normalVelocityLowPass', normalVelocityLowPass)
          call mpas_pool_get_array(timeFiltersAMPool, 'normalVelocityHighPass', normalVelocityHighPass)

          ! get normal velocity
          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=1)

          ! initialize normal velocities
          normalVelocityLowPass(:,:) = normalVelocity(:,:)
          normalVelocityHighPass(:,:) = normalVelocity(:,:)

          block => block % next
        end do

      end if

#ifdef MPAS_DEBUG
      ! get index for edge nearest to a location
      block => domain % blocklist
      iBlock = 0
      do while (associated(block))
        iBlock = iBlock + 1
        call mpas_pool_get_subpool(block % structs, 'mesh', statePool)
        call mpas_pool_get_array(statePool, 'latEdge', latEdge)
        call mpas_pool_get_array(statePool, 'lonEdge', lonEdge)
        call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)

        do i=1,nEdgesSolve
          dist = sqrt((latEdge(i) - latEdgePoint)**2 + (lonEdge(i) - lonEdgePoint)**2)
          if (dist < distmax) then
            distmax = dist
            iEdgeOutput = i
            iBlockOutput = iBlock
          end if
        end do

        block => block % next
      end do

      block => domain % blocklist
      ! get the right block number
      do i=1,iBlockOutput-1
        block => block % next
      end do
      call mpas_pool_get_subpool(block % structs, 'mesh', statePool)
      call mpas_pool_get_array(statePool, 'latEdge', latEdge)
      call mpas_pool_get_array(statePool, 'lonEdge', lonEdge)
      write(stderrUnit,*) 'lon = ', 180.0_RKIND/pii*lonEdge(iEdgeOutput), ' lat = ', 180.0_RKIND/pii*latEdge(iEdgeOutput), ' iklevel=',iklevel
#endif

#ifdef MPAS_DEBUG
        write(stderrUnit,*) 'finished initializing time filters'
#endif

   end subroutine ocn_init_time_filters!}}}

!***********************************************************************
!
!  routine ocn_compute_time_filters
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    07/17/2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_time_filters(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: timeFiltersAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: timeFiltersAM
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalVelocityLowPass, normalVelocityHighPass, normalVelocityTest
      integer, pointer :: nVertLevels, nEdgesSolve
      integer :: k, iEdge
      integer, dimension(:), pointer :: maxLevelEdgeBot

      type (MPAS_timeInterval_type) :: timeStepESMF
      character(len=StrKIND), pointer :: config_dt
      real (kind=RKIND) :: dt, tau
#ifdef MPAS_DEBUG
      integer :: iBlock
#endif

      err = 0

      dminfo = domain % dminfo

#ifdef MPAS_DEBUG
        write(stderrUnit,*) 'start computing time filters'
#endif

      ! get dt
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)
      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)
      ! get tau
      call mpas_pool_get_config(domain % configs, 'config_AM_timeFilters_tau', config_dt)
      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=tau)

#ifdef MPAS_DEBUG
          !write(stderrUnit,*) 'dt = ', dt, ' tau = ', tau
#endif

      block => domain % blocklist
#ifdef MPAS_DEBUG
      iBlock = 0
#endif
      do while (associated(block))
#ifdef MPAS_DEBUG
         iBlock = iBlock + 1
#endif
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'timeFiltersAM', timeFiltersAMPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)

         ! get high and low pass velocity components
         call mpas_pool_get_array(timeFiltersAMPool, 'normalVelocityLowPass', normalVelocityLowPass)
         call mpas_pool_get_array(timeFiltersAMPool, 'normalVelocityHighPass', normalVelocityHighPass)
         call mpas_pool_get_array(timeFiltersAMPool, 'normalVelocityFilterTest', normalVelocityTest)
         ! get normal velocity
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=1)

         ! perform filter computations (in place)
         do iEdge = 1,nEdgesSolve
            do k = 1, maxLevelEdgeBot(iEdge)
              normalVelocityLowPass(k,iEdge) = normalVelocityLowPass(k,iEdge)*(1.0_RKIND - dt/tau) + dt/tau*normalVelocity(k,iEdge)
              normalVelocityHighPass(k,iEdge) = normalVelocity(k,iEdge) - normalVelocityLowPass(k,iEdge)
              ! normalVelocityTest line can possibly be removed (needed for testing purposes)
              normalVelocityTest(k,iEdge) = normalVelocity(k,iEdge)
            end do
#ifdef MPAS_DEBUG
            if (iEdge == iEdgeOutput .and. iBlock == iBlockOutput) then
                write(stderrUnit,*) 'vl=', normalVelocityLowPass(iklevel, iEdge), ' v=', normalVelocity(iklevel, iEdge)
            end if
#endif
         end do

         block => block % next
      end do

#ifdef MPAS_DEBUG
        write(stderrUnit,*) 'finished computing time filters'
#endif

   end subroutine ocn_compute_time_filters!}}}

!***********************************************************************
!
!  routine ocn_restart_time_filters
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    07/17/2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_time_filters(domain, err)!{{{

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

   end subroutine ocn_restart_time_filters!}}}

!***********************************************************************
!
!  routine ocn_finalize_time_filters
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    07/17/2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_time_filters(domain, err)!{{{

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

   end subroutine ocn_finalize_time_filters!}}}

end module ocn_time_filters

! vim: foldmethod=marker
