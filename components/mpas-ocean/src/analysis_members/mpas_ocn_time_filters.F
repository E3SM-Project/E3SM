! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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

!#define TIME_FILTERS_DEBUG

module ocn_time_filters

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_vector_reconstruction

   use ocn_constants
   use ocn_config
#ifdef TIME_FILTERS_DEBUG
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
#ifdef TIME_FILTERS_DEBUG
   integer :: iEdgeOutput = 0, iBlockOutput = 0, iklevel = 1
   real (kind=RKIND) :: lonEdgePoint = 10.0_RKIND*pi/180.0_RKIND, latEdgePoint = 30.0_RKIND*pi/180.0_RKIND
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
      type (mpas_pool_type), pointer :: timeFiltersAMPool, statePool
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalVelocityLowPass, normalVelocityHighPass

#ifdef TIME_FILTERS_DEBUG
      real (kind=RKIND), dimension(:), pointer :: lonEdge, latEdge
      real (kind=RKIND) :: dist, distmax = 1e9
      integer :: i, iBlock
      integer, pointer :: nEdgesSolve
#endif

      err = 0

      if (config_AM_timeFilters_initialize_filters) then
#ifdef TIME_FILTERS_DEBUG
        call mpas_log_write( 'initializing time filters')
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

#ifdef TIME_FILTERS_DEBUG
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
      print *,  'lon = ', 180.0_RKIND/pi*lonEdge(iEdgeOutput), ' lat = ', 180.0_RKIND/pi*latEdge(iEdgeOutput), &
                          ' iklevel=',iklevel, ' iEdgeOutput=',iEdgeOutput, ' iBlockOutput = ', iBlockOutput
#endif

#ifdef TIME_FILTERS_DEBUG
        call mpas_log_write( 'finished initializing time filters')
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
      type (mpas_pool_type), pointer :: timeFiltersAM
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalVelocityLowPass, normalVelocityHighPass, &
                                                    normalVelocityTest

      type (field2DReal), pointer :: normalVelocityLowPassField, normalVelocityHighPassField

      real (kind=RKIND), dimension(:,:), pointer :: velocityZonalLowPass, velocityMeridionalLowPass, &
                                                    velocityXLowPass, velocityYLowPass, velocityZLowPass, &
                                                    velocityZonalHighPass, velocityMeridionalHighPass, &
                                                    velocityXHighPass, velocityYHighPass, velocityZHighPass
      integer, pointer :: nVertLevels, nEdgesSolve
      integer :: k, iEdge
      integer, dimension(:), pointer :: maxLevelEdgeBot

      type (MPAS_timeInterval_type) :: timeStepESMF
      real (kind=RKIND) :: dt, tau
#ifdef TIME_FILTERS_DEBUG
      integer :: iBlock
#endif

      err = 0

      dminfo = domain % dminfo

#ifdef TIME_FILTERS_DEBUG
        call mpas_log_write( 'start computing time filters')
#endif

      ! get dt
      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dt)
      ! get tau
      call mpas_set_timeInterval(timeStepESMF, timeString=config_AM_timeFilters_tau, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=tau)

#ifdef TIME_FILTERS_DEBUG
          !print *,  'dt = ', dt, ' tau = ', tau
#endif

      block => domain % blocklist
#ifdef TIME_FILTERS_DEBUG
      iBlock = 0
#endif
      do while (associated(block))
#ifdef TIME_FILTERS_DEBUG
         iBlock = iBlock + 1
#endif
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
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
#ifdef TIME_FILTERS_DEBUG
            if (iEdge == iEdgeOutput .and. iBlock == iBlockOutput) then
                print *,  'vl=', normalVelocityLowPass(iklevel, iEdge), ' v=', normalVelocity(iklevel, iEdge)
            end if
#endif
         end do
         ! exchange halo information in order to ensure that particles on halo are advected properly
         call mpas_pool_get_field(timeFiltersAMPool, 'normalVelocityLowPass', normalVelocityLowPassField)
         call mpas_pool_get_field(timeFiltersAMPool, 'normalVelocityHighPass', normalVelocityHighPassField)
         call mpas_dmpar_exch_halo_field(normalVelocityLowPassField)
         call mpas_dmpar_exch_halo_field(normalVelocityHighPassField)

         block => block % next
      end do

      ! do IO communications if this is an output time step
      if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='timeFiltersOutput', &
          direction=MPAS_STREAM_OUTPUT, ierr=err) .and. config_AM_timeFilters_compute_cell_centered_values) then
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'timeFiltersAM', timeFiltersAMPool)
           ! get variables for computations
           call mpas_pool_get_array(timeFiltersAMPool, 'normalVelocityLowPass', normalVelocityLowPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityZonalLowPass', velocityZonalLowPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityMeridionalLowPass', velocityMeridionalLowPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityXLowPass', velocityXLowPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityYLowPass', velocityYLowPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityZLowPass', velocityZLowPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'normalVelocityHighPass', normalVelocityHighPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityZonalHighPass', velocityZonalHighPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityMeridionalHighPass', velocityMeridionalHighPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityXHighPass', velocityXHighPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityYHighPass', velocityYHighPass)
           call mpas_pool_get_array(timeFiltersAMPool, 'velocityZHighPass', velocityZHighPass)
          ! must perform reconstruction for cell centered values
           call mpas_reconstruct(meshPool, normalVelocityLowPass,  &
                            velocityXLowPass, velocityYLowPass, velocityZLowPass,   &
                            velocityZonalLowPass, velocityMeridionalLowPass, &
                            includeHalos = .false.)
          ! must perform reconstruction for cell centered values
           call mpas_reconstruct(meshPool, normalVelocityHighPass,  &
                            velocityXHighPass, velocityYHighPass, velocityZHighPass,   &
                            velocityZonalHighPass, velocityMeridionalHighPass, &
                            includeHalos = .false.)
          block => block % next
        end do
      end if

#ifdef TIME_FILTERS_DEBUG
        call mpas_log_write( 'finished computing time filters')
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
