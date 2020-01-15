! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  ocn_lagrangian_particle_tracking_reset
!
!> \brief   LIGHT reset functionality
!> \author  Phillip J. Wolfram
!> \date    10/28/2015
!> \details
!> This module provides routines for performing particle resets in LIGHT.
!
!-----------------------------------------------------------------------
module ocn_lagrangian_particle_tracking_reset

  use mpas_derived_types
  use mpas_constants
  use mpas_timekeeping
  use mpas_stream_manager
  use mpas_pool_routines

  use ocn_constants

  implicit none
  private

  !-----------------------------------------------------------------
  ! public routines and interfaces
  !-----------------------------------------------------------------
  ! define publically accessible subroutines, functions, interfaces
  public :: ocn_setup_particle_reset_condition
  public :: ocn_evaluate_particle_reset_condition
  public :: ocn_finalize_particle_reset_condition

  contains

!***********************************************************************
!
!  routine ocn_setup_particle_reset_condition
!
!> \brief   Set up needed information for particle resets
!> \author  Phillip Wolfram
!> \date    10/28/2015
!> \details
!>  Purpose: Perform set up for particle resets.
!>  Input:   domain
!-----------------------------------------------------------------------
  subroutine ocn_setup_particle_reset_condition(domain, err) !{{{

     implicit none
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
     type (mpas_pool_type), pointer :: lagrPartTrackScalarPool, lagrPartTrackRegionsPool
     real (kind=RKIND), pointer :: globalResetTimeValue
     type (mpas_timeInterval_type) :: timeInterval
     character (len=StrKIND), pointer :: config_AM_lagrPartTrack_reset_global_timestamp
     character (len=StrKIND), pointer :: config_AM_lagrPartTrack_region_stream
     character (len=StrKIND), pointer :: config_AM_lagrPartTrack_reset_criteria
     type (field1DInteger), pointer :: resetInsideRegionMaskValue1Field, resetOutsideRegionMaskValue1Field
     integer, dimension(:), pointer :: resetInsideRegionMaskValue1, resetOutsideRegionMaskValue1

     err = 0

     ! get the configuration options
      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_reset_global_timestamp', &
        config_AM_lagrPartTrack_reset_global_timestamp)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_region_stream', &
        config_AM_lagrPartTrack_region_stream)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_reset_criteria', &
        config_AM_lagrPartTrack_reset_criteria)

     ! load in region masks streams (masks stored in pool)
     if (trim(config_AM_lagrPartTrack_reset_criteria) == 'region' .or. &
         trim(config_AM_lagrPartTrack_reset_criteria) == 'all' &
        ) then
       call MPAS_stream_mgr_read(domain % streamManager, streamID=trim(config_AM_lagrPartTrack_region_stream), ierr=err)
     end if

     ! convert input config_AM_lagrPartTrack_reset_global_timestamp into S for calculations
     block => domain % blocklist
     do while (associated(block))
       ! setup pointers / get block
       call mpas_pool_get_subpool(block % structs, 'lagrPartTrackScalars', lagrPartTrackScalarPool)
       call mpas_pool_get_array(lagrPartTrackScalarPool, 'globalResetTimeValue', globalResetTimeValue)

       ! convert config_AM_lagrPartTrack_reset_global_timestamp into seconds and store in globalResetTimeValue
       call mpas_set_timeInterval(timeInterval, timeString=trim(config_AM_lagrPartTrack_reset_global_timestamp))
       call mpas_get_timeInterval(timeInterval, dt=globalResetTimeValue)

       !print *,  'resetTimeValue = ', globalResetTimeValue

       block => block %  next
     end do

  end subroutine ocn_setup_particle_reset_condition!}}}

!***********************************************************************
!
!  routine ocn_evaluate_particle_reset_condition
!
!> \brief   Evaluate needed information for particle resets
!> \author  Phillip Wolfram
!> \date    10/30/2015
!> \details
!>  Purpose: Evaluate if particle resets should occur for a particle
!>  Input:   domain, particle
!>  Output:  boolean specifying whether the particles should be reset.
!-----------------------------------------------------------------------
  subroutine ocn_evaluate_particle_reset_condition(domain, block, particle, dt, iCell, resetParticle, err) !{{{

     implicit none
     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     real (kind=RKIND), intent(in) :: dt

     !-----------------------------------------------------------------
     !
     ! input/output variables
     !
     !-----------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     type (block_type), intent(inout), pointer :: block
     type (mpas_particle_type), pointer, intent(inout) :: particle
     integer, intent(inout) :: iCell

     !-----------------------------------------------------------------
     !
     ! output variables
     !
     !-----------------------------------------------------------------

     logical, intent(out) :: resetParticle
     integer, intent(out) :: err !< Output: error flag

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     type (mpas_pool_type), pointer :: lagrPartTrackScalarPool, lagrPartTrackRegionsPool
     integer, pointer :: transfered, numTimesReset
     integer, pointer :: currentBlock, currentBlockReset, currentCell, currentCellReset
     real (kind=RKIND), pointer :: xParticleReset, yParticleReset, zParticleReset, zLevelParticleReset
     real (kind=RKIND), pointer :: xParticle, yParticle, zParticle, zLevelParticle
     real (kind=RKIND), pointer :: timeSinceReset
!     real (kind=RKIND), pointer :: sumU, sumV, sumUU, sumUV, sumVV
     integer, pointer :: resetTime
     real (kind=RKIND), pointer :: globalResetTimeValue

     character (len=StrKIND), pointer :: config_AM_lagrPartTrack_reset_criteria
     logical, pointer :: config_AM_lagrPartTrack_reset_if_outside_region
     logical, pointer :: config_AM_lagrPartTrack_reset_if_inside_region
     integer, dimension(:), pointer :: resetInsideRegionMaskValue1, resetOutsideRegionMaskValue1

     ! initialize outputs
     err = 0
     resetParticle = .False.

     ! get config options
      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_reset_criteria', &
        config_AM_lagrPartTrack_reset_criteria)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_reset_if_outside_region', &
        config_AM_lagrPartTrack_reset_if_outside_region)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_reset_if_inside_region', &
        config_AM_lagrPartTrack_reset_if_inside_region)

     ! get variables
     call mpas_pool_get_array(particle % haloDataPool, 'timeSinceReset', timeSinceReset)
     call mpas_pool_get_array(particle % haloDataPool, 'resetTime', resetTime)

     call mpas_pool_get_subpool(block % structs, 'lagrPartTrackScalars', lagrPartTrackScalarPool)
     call mpas_pool_get_array(lagrPartTrackScalarPool, 'globalResetTimeValue', globalResetTimeValue)

     if (trim(config_AM_lagrPartTrack_reset_criteria) == 'region' .or. &
         trim(config_AM_lagrPartTrack_reset_criteria) == 'all' &
        ) then
      call mpas_pool_get_subpool(block % structs, 'lagrPartTrackRegions', lagrPartTrackRegionsPool)
      call mpas_pool_get_array(lagrPartTrackRegionsPool, 'resetInsideRegionMaskValue1', resetInsideRegionMaskValue1)
      call mpas_pool_get_array(lagrPartTrackRegionsPool, 'resetOutsideRegionMaskValue1', resetOutsideRegionMaskValue1)
    end if

     ! advance particle time
     timeSinceReset = timeSinceReset + dt

     ! determine whether reset should occur depending upon type of reset condition
     select case (trim(config_AM_lagrPartTrack_reset_criteria))

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! time based
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !{{{

       case ('particle_time')
         ! use particle's value for resetTime and timeSinceReset
         if (timeSinceReset > resetTime) then
           resetParticle = .True.
         end if

       case ('global_time')
         if (timeSinceReset > globalResetTimeValue) then
           resetParticle = .True.
         end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !}}}
       ! region based
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !{{{

       case ('region')
         ! outside region
         if (config_AM_lagrPartTrack_reset_if_outside_region .and. &
             resetOutsideRegionMaskValue1(iCell) == 0) then
           resetParticle = .True.
         end if
         ! inside region
         if (config_AM_lagrPartTrack_reset_if_inside_region .and. &
             resetInsideRegionMaskValue1(iCell) == 1) then
           resetParticle = .True.
         end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !}}}
       ! all conditions
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !{{{

       case ('all')
             ! particle time
         if ((timeSinceReset > resetTime) .or. &
             ! global time
             (timeSinceReset > globalResetTimeValue) .or. &
             ! outside region
             (config_AM_lagrPartTrack_reset_if_outside_region .and. &
              resetOutsideRegionMaskValue1(iCell) == 0) .or. &
             ! inside region
             (config_AM_lagrPartTrack_reset_if_inside_region .and. &
              resetInsideRegionMaskValue1(iCell) == 1)) then
           resetParticle = .True.
         end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !}}}
       ! default
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case default
         call mpas_log_write( 'WARNING: reset criteria in config_AM_lagrPartTrack_reset_criteria=' // &
           trim(config_AM_lagrPartTrack_reset_criteria) // ' unknown! Cannot restart.')

     end select

     ! reset particle block, cell, and position to reset values
     if (resetParticle) then

       !call mpas_log_write( 'reseting particle')

       ! get data
       call mpas_pool_get_array(particle % haloDataPool, 'currentBlock', currentBlock)
       call mpas_pool_get_array(particle % haloDataPool, 'currentBlockReset', currentBlockReset)
       call mpas_pool_get_array(particle % haloDataPool, 'currentCellReset', currentCellReset)
       call mpas_pool_get_array(particle % haloDataPool, 'xParticleReset', xParticleReset)
       call mpas_pool_get_array(particle % haloDataPool, 'yParticleReset', yParticleReset)
       call mpas_pool_get_array(particle % haloDataPool, 'zParticleReset', zParticleReset)
       call mpas_pool_get_array(particle % haloDataPool, 'zLevelParticleReset', zLevelParticleReset)
       call mpas_pool_get_array(particle % haloDataPool, 'xParticle', xParticle)
       call mpas_pool_get_array(particle % haloDataPool, 'yParticle', yParticle)
       call mpas_pool_get_array(particle % haloDataPool, 'zParticle', zParticle)
       call mpas_pool_get_array(particle % haloDataPool, 'zLevelParticle', zLevelParticle)
       call mpas_pool_get_array(particle % haloDataPool, 'numTimesReset', numTimesReset)
       call mpas_pool_get_array(particle % haloDataPool, 'transfered', transfered)
!       call mpas_pool_get_array(particle % haloDataPool, 'sumU', sumU)
!       call mpas_pool_get_array(particle % haloDataPool, 'sumV', sumV)
!       call mpas_pool_get_array(particle % haloDataPool, 'sumUU', sumUU)
!       call mpas_pool_get_array(particle % haloDataPool, 'sumUV', sumUV)
!       call mpas_pool_get_array(particle % haloDataPool, 'sumVV', sumVV)

       ! reset the time
       timeSinceReset = 0.0_RKIND

       ! increment counters
       if (currentBlock /= currentBlockReset) then
         transfered = transfered + 1
       end if
       numTimesReset = numTimesReset + 1

       ! reset the block and the current cell
       currentBlock = currentBlockReset
       ! this should be a -1 in general but could precache based on an initial decomposition for performance
       iCell = -1
       !iCell = currentCellReset

       ! reset positions
       xParticle = xParticleReset
       yParticle = yParticleReset
       zParticle = zParticleReset
       zLevelParticle = zLevelParticleReset

!       ! reset velocity sums
!       sumU = 0.0_RKIND
!       sumV = 0.0_RKIND
!       sumUU = 0.0_RKIND
!       sumUV = 0.0_RKIND
!       sumVV = 0.0_RKIND

       ! more variables may need to be reset in the future

     end if

   end subroutine ocn_evaluate_particle_reset_condition!}}}

!***********************************************************************
!
!  routine ocn_finalize_particle_reset_condition
!
!> \brief   Finalize information for particle resets
!> \author  Phillip Wolfram
!> \date    10/30/2015
!> \details
!>  Purpose: Finalize setup of particle resets
!>  Input:   domain
!-----------------------------------------------------------------------
  subroutine ocn_finalize_particle_reset_condition(domain, err) !{{{

     implicit none
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

     ! particle reset cleanup

  end subroutine ocn_finalize_particle_reset_condition!}}}

end module ocn_lagrangian_particle_tracking_reset

