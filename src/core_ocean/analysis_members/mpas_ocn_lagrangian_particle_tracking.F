! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_lagrangian_particle_tracking
!
!> \brief MPAS ocean analysis mode member: lagrangian_particle_tracking
!> \author Phillip J. Wolfram
!> \date   02/20/14 and 07/23/2015
!> \details
!>  MPAS ocean analysis core member: lagrangian particle tracking
!>  module computes Lagrangian particle trajectories and associated
!>  diagnostics
!-----------------------------------------------------------------------

module ocn_lagrangian_particle_tracking

   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants

   use ocn_particle_list
   use ocn_lagrangian_particle_tracking_interpolations

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

   public :: ocn_init_lagrangian_particle_tracking, &
             ocn_compute_lagrangian_particle_tracking, &
             ocn_restart_lagrangian_particle_tracking, &
             ocn_finalize_lagrangian_particle_tracking

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   ! timer
   type (timer_node), pointer :: timerTotalLPT, timerInitAlarms, timerInit, timerComputeStartup, timerCompute, timerWrite, timerRestart, timerFinalize, timerValidatedCell
   type (timer_node), pointer :: timerValidatedCell_out, timerValidatedCell_init, timerVerticalID, timerTransferParticles, timerVelTimeInterp
   type (timer_node), pointer :: timerVelTimeInterp_out, timerHorizMovement, timerUpdateIOHalo, timerConvertXYZLatLon, timerCellID, timerCellID_init, timerReconstFilter
   type (timer_node), pointer :: timerReconstFilter_init, timerVelocityPotDensity, timerTimeStepLPT, timerTransferParticles_init, timerTransferParticles_write, timerMemTasksLPT
   type (timer_node), pointer ::  timerSinglePartStats, timerParticleAssignment, timerHorizVelInterp
   ! neighboring processors and arrays for send numbers / recvs, total number of neighboring processors to a particular processor
   ! allocated in init, deallocated in finalize
   ! these globals could be moved to the framework component, as well
   ! as quite a few of the subroutines
   integer, dimension(:), pointer :: g_ProcNeighs => null(), g_ioProcNeighs=>null()

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_lagrangian_particle_tracking
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    02/20/14
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_lagrangian_particle_tracking(domain, err)!{{{

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
      logical, pointer :: config_do_restart

      err = 0

      write(stderrUnit,*) 'starting ocn_init_lagrangian_particle_tracking'
      call mpas_timer_start("totalLPT", .false., timerTotalLPT)
      call mpas_timer_start("initLPT", .false., timerInit)

      ! resets likely are unnecessary
      !call mpas_stream_mgr_reset_alarms(stream_manager, streamID='lagrPartTrackInput', ierr=err)
      !call mpas_stream_mgr_reset_alarms(stream_manager, streamID='lagrPartTrackRestart', ierr=err)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! build up the particles lists and fill them with their data
      call mpas_particle_list_build_and_assign_particle_list(domain, err)
      ! now we have built the particlelist for a given block
#ifdef MPAS_DEBUG
      write(stderrUnit,*) 'finished building and allocating particles in init'
#endif

      ! parallel code:

      ! get "MPI halos" for communication of particles in halo during computational step
      call mpas_particle_list_build_computation_halos(domain, err, g_ProcNeighs)
#ifdef MPAS_DEBUG
      write(stderrUnit,*) 'finished building and computational halos'
#endif

      ! get "MPI halos" for IO communication during write and restart steps (ioBlock to currentBlocks)
      call mpas_particle_list_build_io_halos(domain, err, 'currentBlock', g_ioProcNeighs)
#ifdef MPAS_DEBUG
      write(stderrUnit,*) 'g_ioProcNeighs=', g_ioProcNeighs
      write(stderrUnit,*) 'finished building io halos'
#endif

      ! transfer particles to their appropriate blocks (currentBlock) via MPI
      ! note, don't necessarily need to have g_ionSend and g_ionRecv comeout
#ifdef MPAS_DEBUG
      call mpas_timer_start("trans_from_block_to_blockLPT", .false., timerTransferParticles_init)
      call mpas_particle_list_test_numparticles_to_neighprocs(domain % dminfo % my_proc_id, g_ProcNeighs, g_ioProcNeighs)
#endif
      call mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, .True., .False., 'currentBlock', &
        g_ioProcNeighs)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("trans_from_block_to_blockLPT", timerTransferParticles_init)
      !call MPI_Barrier(domain % dminfo % comm, err)
#endif

      ! tests to make sure all the values are ok !{{{
#ifdef MPAS_DEBUG
      call mpas_particle_list_test_neighscalc(domain, err)
      call mpas_particle_list_test_numparticles_to_neighprocs(domain % dminfo % my_proc_id, g_ProcNeighs, g_ioProcNeighs)
      call mpas_particle_list_test_num_current_particlelist(domain)
#endif
      !}}}

      ! now set sums for autocorrelation calculation to be zero
      call zero_autocorrelation_sums(domain)

      ! previous compute startup calls
      call intialize_wachspress_coefficients(domain, err)
      call initalize_fields(domain, err)
      call compute_velocity_on_potentialdensity_surface(domain,err,1)
      call compute_velocity_on_potentialdensity_surface(domain,err,2)
      call initialize_particle_properties(domain,2,err)
      call write_lagrangian_particle_tracking(domain, err)

      write(stderrunit,*) 'finished ocn_init_lagrangian_particle_tracking'
      call mpas_timer_stop("initLPT", timerInit)

      call mpas_timer_stop("totalLPT", timerTotalLPT)

   end subroutine ocn_init_lagrangian_particle_tracking!}}}


!***********************************************************************
!
!  routine ocn_compute_lagrangian_particle_tracking
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    02/20/14
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_compute_lagrangian_particle_tracking(domain, timeLevel, err)!{{{

      implicit none
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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: lagrPartTrackFieldsPool, lagrPartTrackCellsPool, lagrPartTrackScratchPool
      integer, dimension(:), pointer :: cellOwnerBlock
      integer, pointer :: currentBlock, ioBlock, indexToParticleID, transfered
      type (mpas_particle_list_type), pointer :: particlelist
      type (mpas_particle_type), pointer :: particle

      real (kind=RKIND), dimension(3) :: particlePosition, particleVelocity
      real (kind=RKIND), pointer :: xParticle, yParticle, zParticle, lonVel, latVel, buoyancyParticle, sumU, sumV, sumUU, sumUV, sumVV
      real (kind=RKIND), dimension(3) :: xSubStep, diffSubStep, diffParticlePosition
      real (kind=RKIND), pointer :: zLevelParticle
      real (kind=RKIND), dimension(:,:), pointer :: zTop, vertVelocityTop, zMid, areaBArray
      real (kind=RKIND), dimension(:), pointer :: bottomDepth
      type (field2DReal), pointer :: normalVelocity, uVertexVelocity, vVertexVelocity, wVertexVelocity, layerThickness

      real (kind=RKIND), dimension(:,:), pointer :: uVertexVelocityArray, vVertexVelocityArray, wVertexVelocityArray, buoyancyTimeInterp, potentialDensity

      real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex

      integer, dimension(:,:), pointer :: verticesOnCell, boundaryVertex
      integer, dimension(:), pointer :: maxLevelCell

      integer theVertex, iLevel, iLevelBuoyancy, aVertex, &
        nSteps, timeStep, subStep, subStepOrder, timeInterpOrder, aTimeLevel, nCellVertices, blockProc, arrayIndex
      integer, pointer :: nCells, nVertLevels, iCell
      integer, dimension(:), pointer :: nCellVerticesArray
      integer, dimension(:,:), pointer :: cellsOnCell
      logical, dimension(:,:), pointer :: ioProcRecvList
      logical, dimension(:), pointer :: ioProcSendList
      logical, pointer :: onSphere

      real (kind=RKIND), dimension(4) :: kWeightK, kWeightKVert
      real (kind=RKIND), dimension(3,4) :: kWeightX
      real (kind=RKIND), dimension(4) :: kWeightXVert
      real (kind=RKIND), dimension(4) :: kWeightT, kWeightTVert
      real (kind=RKIND), dimension(3,5) :: kCoeff
      real (kind=RKIND), dimension(5) :: kCoeffVert
      real (kind=RKIND), dimension(2) :: timeCoeff
      real (kind=RKIND) :: dt, dtSim, tSubStep
      real (kind=RKIND), pointer :: dtParticle
      real (kind=RKIND) :: zSubStep
      real (kind=RKIND) :: diffSubStepVert, diffParticlePositionVert, particleVelocityVert, verticalVelocityInterp
      real (kind=RKIND) :: buoyancyInterp
      real (kind=RKIND), pointer :: sphereRadius
      integer, pointer :: verticalTreatment, vertexReconstMethod, timeIntegration, indexLevel, filterNum
      character(len=StrKIND), pointer :: config_dt
      type (MPAS_timeInterval_type) :: timeStepESMF
      integer :: err_tmp

      err = 0

      !! don't do compute for debugging purposes
      !write(stderrUnit,*) 'Computing Lagrangian Particle Tracking -- NO COMPUTATION!'
      !return

      dminfo = domain % dminfo

      write(stderrUnit,*) 'Computing Lagrangian Particle Tracking...'
      call mpas_timer_start("totalLPT", .false., timerTotalLPT)
      call mpas_timer_start("computeLPT", .false., timerCompute)

      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_filter_number', filterNum)

      allocate(ioProcRecvList(size(g_ioProcNeighs), domain % dminfo % nprocs))
      allocate(ioProcSendList(domain % dminfo % nprocs))
      ioProcRecvList = .False.
      ioProcSendList = .False.

      ! get the most recent velocities on the potential density surfaces
#ifdef MPAS_DEBUG
      call mpas_timer_start("velocity_pot_density_LPT", .false., timerVelocityPotDensity)
#endif
      call compute_velocity_on_potentialdensity_surface(domain,err,2)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("velocity_pot_density_LPT",timerVelocityPotDensity)
#endif

      block => domain % blocklist
      do while (associated(block)) !{{{
#ifdef MPAS_DEBUG
        call mpas_timer_start("memtasksLPT", .false., timerMemTasksLPT)
#endif
        ! allocate scratch memory / setup pointers / get block !{{{
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackFields', lagrPartTrackFieldsPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackScratch', lagrPartTrackScratchPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)

        ! particlelist should be stored in the structs pool probably (need to seriously rework the code!)
        particlelist => block % particlelist

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'zCell', zCell)
        call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
        call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
        call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'wachspressAreaB', areaBArray)
        call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
        call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
        call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
        call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', buoyancyTimeInterp)

        ! note, originally this was diagnostics % state % normalVelocity (without time level), but
        ! now there is a time level so selection of the correct time level appears to be tricky
        ! the issue is large, nearly NAN normalVelocities with timeLevel=2
        call mpas_pool_get_field(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
        ! potentially a costly function call, but ensures halos are ok
        !call mpas_dmpar_exch_halo_field(normalVelocity)
        call mpas_pool_get_field(statePool, 'layerThickness', layerThickness, timeLevel=timeLevel)

        call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'boundaryVertex', boundaryVertex)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nCellVerticesArray)
        call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

        call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
        call mpas_pool_get_config(meshPool, 'sphere_radius', sphereRadius)
        call mpas_pool_get_config(block % configs, 'config_dt', config_dt)
        call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
        !}}}

        !{{{
        ! flip previous time level to back (switching memory pointers)
        call mpas_pool_shift_time_levels(lagrPartTrackFieldsPool)
        call mpas_pool_get_field(lagrPartTrackFieldsPool, 'uVertexVelocity', uVertexVelocity, timeLevel=2)
        call mpas_pool_get_field(lagrPartTrackFieldsPool, 'vVertexVelocity', vVertexVelocity, timeLevel=2)
        call mpas_pool_get_field(lagrPartTrackFieldsPool, 'wVertexVelocity', wVertexVelocity, timeLevel=2)
#ifdef MPAS_DEBUG
        call mpas_timer_stop("memtasksLPT",timerMemTasksLPT)
#endif

#ifdef MPAS_DEBUG
        call mpas_timer_start("reconst_filter_LPT", .false., timerReconstFilter)
#endif
        ! need to handle periodicity within functions below for vertex reconstruction
        call ocn_vertex_reconstruction(filterNum, meshPool, lagrPartTrackScratchPool, lagrPartTrackCellsPool, &
          layerThickness % array, normalVelocity % array, &
          uVertexVelocity, vVertexVelocity, wVertexVelocity)
#ifdef MPAS_DEBUG
        call mpas_timer_stop("reconst_filter_LPT", timerReconstFilter)
        write(stderrUnit,*) 'uVertexVelocity=',uVertexVelocity % array
        write(stderrUnit,*) 'vVertexVelocity=',vVertexVelocity % array
        write(stderrUnit,*) 'wVertexVelocity=',wVertexVelocity % array
#endif
        !}}}

#ifdef MPAS_DEBUG
        write(stderrUnit,*) 'beginning particle loop with particlelist associated =  ', associated(particlelist)
        ! error is here
        call mpas_particle_list_test_num_current_particlelist(domain)
#endif

        !!!!!!!!!! LOOP OVER PARTICLES !!!!!!!!!!
        ! update the particle position (just from initialized value for now)
        ! this is a loop over particle list and its datastructures
        do while(associated(particlelist)) !{{{
          ! get pointers / option values
          particle => particlelist % particle

          ! get values {{{
#ifdef MPAS_DEBUG
          call mpas_timer_start("memtasksLPT", .false., timerMemTasksLPT)
#endif
          call mpas_pool_get_array(particle % haloDataPool, 'xParticle', xParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'yParticle', yParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'zParticle', zParticle)
          particlePosition(1) = xParticle
          particlePosition(2) = yParticle
          particlePosition(3) = zParticle

          call mpas_pool_get_array(particle % haloDataPool, 'zLevelParticle', zLevelParticle)

          call mpas_pool_get_array(particle % haloDataPool, 'verticalTreatment', verticalTreatment)
          call mpas_pool_get_array(particle % haloDataPool, 'vertexReconstMethod', vertexReconstMethod)
          call mpas_pool_get_array(particle % haloDataPool, 'indexLevel', indexLevel)
          call mpas_pool_get_array(particle % haloDataPool, 'timeIntegration', timeIntegration)
          call mpas_pool_get_array(particle % haloDataPool, 'dtParticle', dtParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'buoyancyParticle', buoyancyParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'indexToParticleID', indexToParticleID)
          call mpas_pool_get_array(particle % haloDataPool, 'currentCell', iCell)

          call mpas_pool_get_array(particle % haloDataPool, 'lonVel', lonVel)
          call mpas_pool_get_array(particle % haloDataPool, 'latVel', latVel)
          call mpas_pool_get_array(particle % haloDataPool, 'sumU', sumU)
          call mpas_pool_get_array(particle % haloDataPool, 'sumV', sumV)
          call mpas_pool_get_array(particle % haloDataPool, 'sumUU', sumUU)
          call mpas_pool_get_array(particle % haloDataPool, 'sumUV', sumUV)
          call mpas_pool_get_array(particle % haloDataPool, 'sumVV', sumVV)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("memtasksLPT",timerMemTasksLPT)
#endif
          ! process timers / particle reset functions here (need a timer and ability to reset particle to
          ! some set of initial values for the reset
          !}}}

          !!!!!!!!!! COMPUTE TIME STEP INFORMATION !!!!!!!!!!
#ifdef MPAS_DEBUG
          call mpas_timer_start("time_step_LPT", .false., timerTimeStepLPT)
#endif
          ! would eventually need to correspond to domain dt, but for now this is a
          ! global constant
          call mpas_get_timeInterval(timeStepESMF, dt=dtSim)
          ! adjust time step for consistency with integer number of steps
          nSteps = ceiling(dtSim/dtParticle)
          dt = dtSim/nSteps

          !!!!!!!!!! ASSIGN TEMPORAL INTEGRATION COEFFICIENTS !!!!!!!!!!
          ! kCoeff   is (3,subStepOrder+1)
          ! kWeightX is subStepOrder
          ! kWeightK is subStepOrder
          ! kWeightT is subStepOrder
          select case (timeIntegration) !{{{
            case(1) ! EE integration
              kWeightK(1) = 0.0_RKIND

              kWeightT(1) = 0.0_RKIND

              kWeightX(:,1) = 1.0_RKIND

              subStepOrder = 1
            case(2) ! RK2 integration
              kWeightK(1) = 0.0_RKIND
              kWeightK(2) = 0.5_RKIND

              kWeightT(1) = 0.0_RKIND
              kWeightT(2) = 0.5_RKIND

              kWeightX(:,1) = 0.0_RKIND
              kWeightX(:,2) = 1.0_RKIND

              subStepOrder = 2
            case(4) ! RK4 integration
              kWeightK(1) = 0.0_RKIND
              kWeightK(2) = 0.5_RKIND
              kWeightK(3) = 0.5_RKIND
              kWeightK(4) = 1.0_RKIND

              kWeightT(1) = 0.0_RKIND
              kWeightT(2) = 0.5_RKIND
              kWeightT(3) = 0.5_RKIND
              kWeightT(4) = 1.0_RKIND

              kWeightX(:,1) = 1.0_RKIND/6.0_RKIND
              kWeightX(:,2) = 1.0_RKIND/3.0_RKIND
              kWeightX(:,3) = 1.0_RKIND/3.0_RKIND
              kWeightX(:,4) = 1.0_RKIND/6.0_RKIND

              subStepOrder = 4
            case default ! RK2 integration
              kWeightK(1) = 0.0_RKIND
              kWeightK(2) = 0.5_RKIND

              kWeightT(1) = 0.0_RKIND
              kWeightT(2) = 0.5_RKIND

              kWeightX(:,1) = 0.0_RKIND
              kWeightX(:,2) = 1.0_RKIND

              subStepOrder = 2
          end select  !}}}

          ! use same integration coefficients for the vertical
          kWeightKVert = kWeightK
          kWeightXVert = kWeightX(1,:)
          kWeightTVert = kWeightT
          kCoeffVert   = kCoeff(1,:)

          !!!!!!!!!! LOOP OVER TIME STEPS !!!!!!!!!!
          do timeStep = 1, nSteps !{{{
            ! kCoeff   is (3,subStepOrder+1)
            kCoeff = 0.0_RKIND
            kCoeffVert = 0.0_RKIND
            ! compute first
            do subStep = 1, subStepOrder !{{{

              !!!!!!!!!! COMPUTE PARTICLE SUBSTEP POSITIONS USE FOR VELOCITY !!!!!!!!!!
              ! horizontal
              xSubStep = particlePosition
              diffSubStep = kWeightK(subStep) * kCoeff(:,subStep)

              if(kWeightK(subStep) /= 0.0_RKIND) then
                ! project substep to correct spherical shell because diffSubStep isn't 0 and particle moves
#ifdef MPAS_DEBUG
                call mpas_timer_start("particle_horizontal_movementLPT", .false., timerHorizMovement)
#endif
                call particle_horizontal_movement(meshPool, xSubStep, diffSubStep)
#ifdef MPAS_DEBUG
                call mpas_timer_stop("particle_horizontal_movementLPT", timerHorizMovement)
#endif
              end if

              ! vertical
              zSubStep = zLevelParticle
              diffSubStepVert = kWeightKVert(subStep) * kCoeffVert(subStep)

              if(kWeightKVert(subStep) /= 0.0_RKIND) then
                zSubStep = zSubStep + diffSubStepVert
              end if

              ! get new time step (tm = (timestep-1)*dt)
              tSubStep = (timeStep-1 + kWeightT(subStep)) * dt

              !!!!!!!!!! GET SPECIFIC CELL INDICES AND GEOMETRY !!!!!!!!!!

              ! determine cell location
#ifdef MPAS_DEBUG
              call mpas_timer_start("get_validated_cell_idLPT", .false., timerValidatedCell)
              write(stderrUnit, *) 'beginning of substeps'
#endif
              call get_validated_cell_id(nCells, xCell,yCell,zCell , xVertex,yVertex,zVertex, &
                xSubStep(1),xSubStep(2),xSubStep(3), meshPool, &
                nCellVerticesArray, verticesOnCell, iCell, nCellVertices, cellsOnCell)
#ifdef MPAS_DEBUG
              write(stderrUnit,*) 'iCell=',iCell
              call mpas_timer_stop("get_validated_cell_idLPT", timerValidatedCell)
#endif

              if(verticalTreatment /= 4) then
                ! other cases using zSubStep for iLevel and vertical interpolation
#ifdef MPAS_DEBUG
                call mpas_timer_start("mpas_get_vertical_idLPT", .false., timerVerticalID)
#endif
                iLevel = mpas_get_vertical_id(maxLevelCell(iCell), zSubStep, zMid(:,iCell))
#ifdef MPAS_DEBUG
                call mpas_timer_stop("mpas_get_vertical_idLPT", timerVerticalID)
                write(stderrUnit,*) 'iLevel=', iLevel
#endif
              end if

              !!!!!!!!!! TEMPORALLY INTERPOLATE TIME FIELD !!!!!!!!!!
              ! use these coefficients to just get the velocity at n
              !timeInterpOrder = 1
              !timeCoeff(1) =1.0_RKIND
              ! use these coefficients to just get the velocity at n+1
              !timeInterpOrder = 2
              !timeCoeff(1) =0.0_RKIND
              !timeCoeff(2) =1.0_RKIND
              ! get interpolation coefficients for linear interpolation in time
              timeInterpOrder = 2
              timeCoeff(1) = tSubStep / dtSim
              timeCoeff(2) = 1.0_RKIND - timeCoeff(1)

              ! ensure that buoyancy is fixed for each run
              buoyancyInterp = buoyancyParticle
              ! return interpolated horizontal velocity "particleVelocity" and vertical velocity "particleVelocityVert"

#ifdef MPAS_DEBUG
              call mpas_timer_start("velocity_time_interpolationLPT", .false., timerVelTimeInterp)
#endif
              call velocity_time_interpolation(particleVelocity, particleVelocityVert, &
                diagnosticsPool, lagrPartTrackFieldsPool, &
                timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
                verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
                xSubStep, zSubStep, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex, meshPool, areaBArray)
#ifdef MPAS_DEBUG
              call mpas_timer_stop("velocity_time_interpolationLPT", timerVelTimeInterp)
#endif

              !!!!!!!!!! FORM INTEGRATION WEIGHTS kj !!!!!!!!!!
              kCoeff(:,subStep+1) = dt * particleVelocity
              kCoeffVert(subStep+1) = dt * particleVelocityVert
            end do

            !!!!!!!!!! UPDATE PARTICLE POSITIONS !!!!!!!!!!
            !!!!!!!!!! HORIZONTAL AND VERTICAL CONSIDERED SEPARATELY !!!!!!!!!!
            diffParticlePosition = 0.0_RKIND
            diffParticlePositionVert = 0.0_RKIND
            do subStep = 1, subStepOrder
              ! first complete particle integration
              diffParticlePosition = diffParticlePosition + kWeightX(:,subStep) * kCoeff(:,subStep+1)
              diffParticlePositionVert = diffParticlePositionVert + kWeightXVert(subStep) * kCoeffVert(subStep+1)
            end do
            ! now, make sure particle position is still on same spherical shell as before
#ifdef MPAS_DEBUG
            call mpas_timer_start("particle_horizontal_movementLPT", .false., timerHorizMovement)
#endif
            call particle_horizontal_movement(meshPool, particlePosition, diffParticlePosition)
#ifdef MPAS_DEBUG
            call mpas_timer_stop("particle_horizontal_movementLPT", timerHorizMovement)
#endif
            ! now can do any vertical movements independent of the horizontal movement
            ! that was just calculated ( probably need to have more output from the vertical_treatment
            ! and aggregate here
            zLevelParticle = zLevelParticle + diffParticlePositionVert

          end do !}}}
#ifdef MPAS_DEBUG
          call mpas_timer_stop("time_step_LPT", timerTimeStepLPT)
#endif
          !}}}

          !!!!!!!!!! PERFORM SAMPLING (VELOCITY, TEMP, SALINITY, ETC) !!!!!!!!!!
          ! this could be done just before the output anyway

          ! need iCell computed for final position
          ! need scalar values interpolated in time to yield single value
          ! probably need to store zMid too, including flipping it
          ! get updated cell location
#ifdef MPAS_DEBUG
          call mpas_timer_start("get_validated_cell_idLPT", .false., timerValidatedCell_out)
          write(stderrUnit,*) 'do sampling'
#endif
          call get_validated_cell_id(nCells, xCell,yCell,zCell , xVertex,yVertex,zVertex, &
            particlePosition(1),particlePosition(2),particlePosition(3), meshPool, &
            nCellVerticesArray, verticesOnCell, iCell, nCellVertices, cellsOnCell)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("get_validated_cell_idLPT", timerValidatedCell_out)
#endif

          if(verticalTreatment == 4) then  !('buoyancySurface') !{{{
            !! determine index level (don't need validated version because that will "fix" values which we may not want
            !! however, if the particle's target buoyancy surface doesn't exist then we will need to
            !! ensure that vertical location is valid, placing particles outside of range of zMid back inside
            !! we don't validate this because we want the code to fail, at least initially, in the case that
            !! the particle is in a cell that does not have the proper buoyancy target because this implies
            !! that the buoyancy tracking mode has completely failed.
            !! need to make sure it is validated for buoyancy particles
            iLevelBuoyancy = mpas_get_vertical_id(maxLevelCell(iCell), buoyancyInterp, buoyancyTimeInterp(:,iCell))
            ! interpolate the scalars now (assumes that scalar value is constant within a particular cell)
            call interp_cell_scalars(iLevelBuoyancy, maxLevelCell(iCell), buoyancyInterp, buoyancyTimeInterp(:,iCell), &
              zMid(:,iCell), zLevelParticle)
            !}}}
          else
            ! make sure final zLevelParticle is ok so that it can't extent past zMid range
#ifdef MPAS_DEBUG
            call mpas_timer_start("mpas_get_vertical_idLPT", .false., timerVerticalID)
#endif
            iLevel = mpas_get_vertical_id(maxLevelCell(iCell), zLevelParticle, zMid(:,iCell))
#ifdef MPAS_DEBUG
            call mpas_timer_stop("mpas_get_vertical_idLPT", timerVerticalID)
#endif
          end if

          ! compute necessary information for autocorrelation !{{{
          ! get the updated velocity
          ! ensure that buoyancy is fixed for each run
          buoyancyInterp = buoyancyParticle
          ! we just need the last part of the velocity field interpolation because we are at the end of the timestep
          timeInterpOrder = 2
          timeCoeff(1) = 0.0_RKIND
          timeCoeff(2) = 1.0_RKIND
          ! return interpolated horizontal velocity "particleVelocity" and vertical velocity "particleVelocityVert"
          ! noting we use the final positions
#ifdef MPAS_DEBUG
          call mpas_timer_start("velocity_time_interpolationLPT", .false., timerVelTimeInterp_out)
#endif
          call velocity_time_interpolation(particleVelocity, particleVelocityVert, &
            diagnosticsPool, lagrPartTrackFieldsPool, &
            timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
            verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
            particlePosition, zLevelParticle, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex, &
            meshPool, areaBArray)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("velocity_time_interpolationLPT", timerVelTimeInterp_out)
#endif
          ! convert horizontal velocity to lat/lon velocity

          ! store velocity for use in computing normalized autocorrelation offline
#ifdef MPAS_DEBUG
          call mpas_timer_start("mpas_convert_xyz_velocity_to_latlonLPT", .false., timerConvertXYZLatLon)
#endif
          call mpas_convert_xyz_velocity_to_latlon(lonVel, latVel, particlePosition, particleVelocity)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("mpas_convert_xyz_velocity_to_latlonLPT", timerConvertXYZLatLon)
#endif

          ! now store components needed to compute integral timescale
#ifdef MPAS_DEBUG
          call mpas_timer_start("storeSingleParticleStats", .false., timerSinglePartStats)
#endif
          sumU = sumU + lonVel
          sumV = sumV + latVel
          sumUU = sumUU + lonVel*lonVel
          sumUV = sumUV + lonVel*latVel
          sumVV = sumVV + latVel*latVel
#ifdef MPAS_DEBUG
          call mpas_timer_stop("storeSingleParticleStats", timerSinglePartStats)
#endif
          !}}}

          ! properly store particle position (because we can't store arrays directly for particles and must
          ! work entirely in vectors
          xParticle = particlePosition(1)
          yParticle = particlePosition(2)
          zParticle = particlePosition(3)

          !!!!!!!!!! PASS PARTICLES FROM PROCESSOR TO PROCESSOR !!!!!!!!!!
          !{{{
          ! 1. determine if iCell is on halo (just set each particle's currentBlock to the correct currentBlock
          ! 2. determine owning block in halo, update particle's currentBlock

          ! determine currentBlock ownership of iCell
          ! and set cellOwnerBlock to be current block
#ifdef MPAS_DEBUG
          call mpas_timer_start("particleAssignments", .false., timerParticleAssignment)
#endif
          ! update halo fields
          call mpas_particle_list_update_computational_halos(domain, block, particle, 'lagrPartTrackCells', iCell, arrayIndex, ioProcRecvList, ioProcSendList, g_ioProcNeighs)
#ifdef MPAS_DEBUG
            call mpas_timer_stop("particleAssignments", timerParticleAssignment)
#endif
          !}}}

          ! get next particle to process on the list
          particlelist => particlelist % next
        end do !}}}

        ! get next block
        block => block % next
      end do !}}}

      !!!!!!!!!! PASS PARTICLES FROM PROCESSOR TO PROCESSOR !!!!!!!!!!
      ! MPI calls
      ! 3. place particle on temporary list to be sent to processor, removing particle from present list
      ! 4. pass list of particles to owning block (outside block loop so that all blocks can be processed)
      ! 5. delete all particles on the temporary lists (potentially if on different proc)
      !    because they have been permenantly moved to block owning the halo
      !
      ! Items 3-5 should be able to be described in terms of current code
      ! noting that the most important thing is to ensure that the particle's currentBlock is
      ! updated.  Then, a routine can be called to make sure particles are placed on their appropriate
      ! currentBlocks

      ! particle transfer can then occur from computational processor to computational processor
#ifdef MPAS_DEBUG
      call mpas_timer_start("trans_from_block_to_blockLPT", .false., timerTransferParticles)
#endif
      call mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, .True., .False., 'currentBlock', &
        g_ProcNeighs)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("trans_from_block_to_blockLPT", timerTransferParticles)
      call mpas_timer_start("update_io_haloLPT", .false., timerUpdateIOHalo)
#endif
      call mpas_particle_list_update_io_halos(domain, err, g_ioProcNeighs, ioProcSendList, ioProcRecvList)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("update_io_haloLPT",timerUpdateIOHalo)
#endif
      deallocate(ioProcSendList, ioProcRecvList)

      ! do IO communications if this is an output time step
      if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='lagrPartTrackOutput', direction=MPAS_STREAM_OUTPUT, ierr=err)) then
        call write_lagrangian_particle_tracking(domain, err)
      end if

      call mpas_timer_stop("computeLPT", timerCompute)
      call mpas_timer_stop("totalLPT", timerTotalLPT)

      write(stderrUnit,*) 'finished computing Lagrangian Particle Tracking'

   end subroutine ocn_compute_lagrangian_particle_tracking!}}}

!***********************************************************************
!
!  routine ocn_restart_lagrangian_particle_tracking
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    02/20/14 and 07/23/15
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_restart_lagrangian_particle_tracking(domain, err)!{{{

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

      type (domain_type), intent(in) :: domain

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

      call mpas_timer_start("totalLPT", .false., timerTotalLPT)
      call mpas_timer_start("restartLPT", .false., timerRestart)

      ! do restart if this is a restart step
      if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='lagrPartTrackRestart', direction=MPAS_STREAM_OUTPUT, ierr=err)) then
        call mpas_timer_start("restartLPT", .false., timerRestart)

      write(stderrUnit,*) 'start ocn_restart_lagrangian_particle_tracking'
      ! transfer particles to their appropriate blocks (ioBlock) via MPI
      ! note, don't necessarily need to have g_ionSend and g_ionRecv comeout
      call mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, .True., .True., 'ioBlock', &
        g_ioProcNeighs)

      ! write out all the data, sorting to make sure that shuffled particles
      ! are ouptut correctly (done separately in each function, could be
      ! pulled out as an optimization)
      ! write halo data out, but don't need nonhalo data because it is
      ! computed for output (diagnostic, not prognostic)
      call mpas_particle_list_write_halo_data(domain, err)
      !call mpas_particle_list_write_nonhalo_data(domain, err)

      ! need to now remove the io particles (remove particles that don't have the
      ! correct currentBlock)
      call mpas_particle_list_remove_particles_not_on_current_block(domain,err)

      write(stderrUnit,*) 'end ocn_restart_lagrangian_particle_tracking'
      call mpas_timer_stop("restartLPT", timerRestart)
      end if

      call mpas_timer_stop("totalLPT", timerTotalLPT)

   end subroutine ocn_restart_lagrangian_particle_tracking!}}}

!***********************************************************************
!
!  routine write_lagrangian_particle_tracking
!
!> \brief   Driver for MPAS-Ocean analysis output
!> \author  Phillip Wolfram
!> \date    02/20/14
!> \details
!>  This routine writes all output for this MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
   subroutine write_lagrangian_particle_tracking(domain, err)!{{{

      implicit none
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

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

      err = 0

      write(stderrUnit,*) 'start write_lagrangian_particle_tracking'
      call mpas_timer_start("totalLPT", .false., timerTotalLPT)
      call mpas_timer_start("writeLPT", .false., timerWrite)

      !call mpas_particle_list_test_num_current_particlelist(domain)

      ! transfer particles to their appropriate blocks (ioBlock) via MPI
      ! note, don't necessarily need to have g_ionSend and g_ionRecv comeout
      !write(stderrUnit,*) 'g_ioProcNeighs = ', g_ioProcNeighs
#ifdef MPAS_DEBUG
      call mpas_timer_start("trans_from_block_to_blockLPT", .false., timerTransferParticles_write)
#endif
      ! depreciated (can just use update_halo_io to keep g_ioProcNeighs up to date)
      !! get "MPI halos" for IO communication during write and restart steps (currentBlock to ioBlock)
      !call mpas_particle_list_build_io_halos(domain, err, 'ioBlock', g_ioProcNeighs)
      ! transfer the data
      call mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, .False., .True., 'ioBlock', &
        g_ioProcNeighs)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("trans_from_block_to_blockLPT", timerTransferParticles_write)
#endif

      ! write out all the data, sorting to make sure that shuffled particles
      ! are ouptut correctly (done separately in each function, could be
      ! pulled out as an optimization)
      call mpas_particle_list_write_halo_data(domain, err)
      call mpas_particle_list_write_nonhalo_data(domain, err)

#ifdef MPAS_DEBUG
      call mpas_particle_list_test_num_current_particlelist(domain)
#endif
      ! need to now remove the io particles (remove particles that don't have the
      ! correct currentBlock)
      call mpas_particle_list_remove_particles_not_on_current_block(domain,err)

#ifdef MPAS_DEBUG
      call mpas_particle_list_test_num_current_particlelist(domain)
#endif
      write(stderrUnit,*) 'end write_lagrangian_particle_tracking'
      call mpas_timer_stop("writeLPT", timerWrite)
      call mpas_timer_stop("totalLPT", timerTotalLPT)


   end subroutine write_lagrangian_particle_tracking!}}}

!***********************************************************************
!
!  routine ocn_finalize_lagrangian_particle_tracking
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Phillip J. Wolfram
!> \date    02/20/14
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_finalize_lagrangian_particle_tracking(domain, err)!{{{

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

      type (domain_type), intent(in) :: domain

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
      integer :: timeLev

      call mpas_timer_start("totalLPT", .false., timerTotalLPT)
      call mpas_timer_start("finalizeLPT", .false., timerFinalize)
      err = 0

      write(stderrUnit,*) 'start ocn_finalize_lagrangian_particle_tracking'
      block => domain % blocklist
      do while (associated(block))
        call mpas_particle_list_destroy_particle_list(block % particlelist)
        block => block %  next
      end do

      deallocate(g_ProcNeighs, g_ioProcNeighs)
      ! these following ones should be deallocated once rest of the code is sketched in
      !deallocate(g_nPartSend, g_nPartRecv, g_ionSend, g_ionRecv)

      write(stderrUnit,*) 'end ocn_finalize_lagrangian_particle_tracking'
      call mpas_timer_stop("finalizeLPT", timerFinalize)
      call mpas_timer_stop("totalLPT", timerTotalLPT)

   end subroutine ocn_finalize_lagrangian_particle_tracking!}}}

!-----------------------------------------------------------------------
!
! PRIVATE SUBROUTINES
!
!-----------------------------------------------------------------------
!{{{

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GET_VALIDATED_CELL_ID
   !
   ! Computes the validated cell ID for a particular location base on proximity to point.
   ! Phillip Wolfram 06/18/2014
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_validated_cell_id(nCells, xCell,yCell,zCell , xVertex,yVertex,zVertex, &
       xSubStep,ySubStep,zSubStep, meshPool, nCellVerticesArray, verticesOnCell, &
       iCell, nCellVertices, cellsOnCell)

     implicit none

     ! intent (in)
     integer, intent(in) :: nCells  !< number of cells
     integer, dimension(:,:), pointer, intent(in) :: verticesOnCell   !< vertex indices on cell
     integer, dimension(:), pointer, intent(in) :: nCellVerticesArray
     real (kind=RKIND), dimension(:), intent(in) :: xCell,yCell,zCell  !< spatial location of cell centers
     real (kind=RKIND), dimension(:), intent(in) :: xVertex,yVertex,zVertex  !< spatial location of cell vertices
     real (kind=RKIND), intent(in) :: xSubStep,ySubStep,zSubStep
     type (mpas_pool_type), intent(in), pointer :: meshPool   ! meshPool pointer
     integer, dimension(:,:), intent(in) :: cellsOnCell       ! cell connectivity

     !intent (out)
     integer, intent(inout) :: iCell
     integer, intent(out) :: nCellVertices
#ifdef MPAS_DEBUG
     logical, pointer :: is_periodic
     real(kind=RKIND), pointer :: x_period, y_period
     logical, pointer :: on_a_sphere
     real(kind=RKIND), dimension(:), pointer :: xtmp, ytmp
     integer :: i
     call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
     call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
     call mpas_pool_get_config(meshPool, 'x_period', x_period)
     call mpas_pool_get_config(meshPool, 'y_period', y_period)
#endif

     ! get cell index
!#ifdef MPAS_DEBUG
!     iCell = -1
!#endif
     call mpas_get_nearby_cell_index(nCells, xCell,yCell,zCell ,  &
       xSubStep,ySubStep,zSubStep, meshPool, iCell, cellsOnCell, nCellVerticesArray)

     nCellVertices = nCellVerticesArray(iCell)

#ifdef MPAS_DEBUG
     ! check to make sure the horizontal location is valid, otherwise report an error
     !write(stderrUnit,*) 'max verticesOnCell = ', maxval(verticesOnCell(:,iCell)), 'nVertices = ', size(xVertex)
     if (on_a_sphere .or. .not. is_periodic) then
       if(.not. point_in_cell(nCellVertices, &
         xVertex(verticesOnCell(1:nCellVertices,iCell)), &
         yVertex(verticesOnCell(1:nCellVertices,iCell)), &
         zVertex(verticesOnCell(1:nCellVertices,iCell)), &
         xSubStep,ySubStep,zSubStep , on_a_sphere)) then
           write(stderrUnit,*) 'Point (', xSubStep,ySubStep,zSubStep ,') is horizontally outside cell ', iCell
           write(stderrUnit,*) 'Cell  (',&
             xCell(iCell),yCell(iCell),zCell(iCell), ') with index ' , iCell
           write(stderrUnit,*) 'xVertex = ', xVertex(verticesOnCell(1:nCellVertices,iCell))
           write(stderrUnit,*) 'yVertex = ', yVertex(verticesOnCell(1:nCellVertices,iCell))
           write(stderrUnit,*) 'zVertex = ', zVertex(verticesOnCell(1:nCellVertices,iCell))
       end if
     else
       allocate(xtmp(nCellVertices), ytmp(nCellVertices))
       do i = 1, nCellVertices
         xtmp(i) = mpas_fix_periodicity(xVertex(verticesOnCell(i,iCell)), xSubStep, x_period)
         ytmp(i) = mpas_fix_periodicity(yVertex(verticesOnCell(i,iCell)), ySubStep, y_period)
       end do
       if(.not. point_in_cell(nCellVertices, xtmp, ytmp, &
         zVertex(verticesOnCell(1:nCellVertices,iCell)), &
         xSubStep,ySubStep,zSubStep , on_a_sphere)) then
           write(stderrUnit,*) 'Point (', xSubStep,ySubStep,zSubStep ,') is horizontally outside cell ', iCell
           write(stderrUnit,*) 'Cell  (',&
             xCell(iCell),yCell(iCell),zCell(iCell), ') with index ' , iCell
           write(stderrUnit,*) 'xVertex = ', xtmp
           write(stderrUnit,*) 'yVertex = ', ytmp
           write(stderrUnit,*) 'zVertex = ', zVertex(verticesOnCell(1:nCellVertices,iCell))
       end if
       deallocate(xtmp, ytmp)
     end if
#endif

   end subroutine get_validated_cell_id

#ifdef MPAS_DEBUG
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION POINT_IN_CELL
   !
   ! Check to make sure point (xp,yp,zp) is within cell iCell (implicit via xv,yv,zv)
   ! Phillip Wolfram 05/01/2014
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   logical function point_in_cell(nVertices, xv,yv,zv , xp,yp,zp, on_a_sphere) !{{{
     implicit none

     integer, intent(in) :: nVertices                         ! number of vertices for cell
     real (kind=RKIND), dimension(:), intent(in) :: xv,yv,zv  ! cell vertex locations
     real (kind=RKIND), intent(in) :: xp,yp,zp                ! point location
     logical, intent(in) :: on_a_sphere                       ! flag designating if we are on a sphere

     integer :: aPoint, v1, v0
     real (kind=RKIND) :: pointRadius   ! magnitude of a point radius
     real (kind=RKIND), dimension(3) :: pPoint                    ! point vector
     real (kind=RKIND), dimension(3, nVertices) :: pVertices      ! vertices vectors
     real (kind=RKIND), dimension(3) :: vec1, vec2, crossProd     ! temporary vectors
     integer, dimension(nVertices+1) :: cyclePlusOne

     ! normalize locations to same spherical shell (unit) for direct comparison
     if (on_a_sphere) then
       pointRadius = sqrt(xp*xp + yp*yp + zp*zp)
       pPoint = (/ xp/pointRadius , yp/pointRadius , zp/pointRadius /)
       do aPoint = 1, nVertices
         pointRadius = sqrt(xv(aPoint)*xv(aPoint) + yv(aPoint)*yv(aPoint) + zv(aPoint)*zv(aPoint))
         pVertices(:,aPoint) = (/ xv(aPoint), yv(aPoint), zv(aPoint) /) / pointRadius
       end do
     else
       pPoint = (/ xp,yp,zp /)
       do aPoint = 1, nVertices
         pVertices(:, aPoint) = (/ xv(aPoint), yv(aPoint), zv(aPoint) /)
       end do
     end if

     ! build up vertex cycle for the cell
     do aPoint = 1, nVertices-1
       cyclePlusOne(aPoint) = aPoint + 1
     end do
     cyclePlusOne(nVertices) = 1

     ! check, using cross-products, that point is within the cell, assuming it is to start
     point_in_cell = .true.
     do aPoint = 1, nVertices
       ! get indices of points
       v0 = aPoint
       v1 = cyclePlusOne(aPoint)

       ! compute the local vectors
       vec1 = pVertices(:,v1) - pVertices(:,v0)
       vec2 = pPoint - pVertices(:,v0)

       ! compute the cross product and dot with normal, if negative we are outside cell
       ! we only need to fail on a single test!
       call mpas_cross_product_in_r3(vec1,vec2,crossProd)
       if (on_a_sphere) then
         if(sum(crossProd*pVertices(:,v0)) < 0) then
           point_in_cell = .false.
         end if
       else
         if(crossProd(3) < 0) then
           point_in_cell = .false.
         end if
       end if

     end do

   end function point_in_cell!}}}
#endif

!***********************************************************************
!
!  routine initalize_fields
!
!> \brief   Initialize fields
!> \author  Phillip Wolfram
!> \date    05/22/2014
!> \details
!>  This routine inializes the fields necessary for particle tracking
!
!-----------------------------------------------------------------------
   subroutine initalize_fields(domain, err)!{{{

      implicit none

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: lagrPartTrackFieldsPool, lagrPartTrackScratchPool, lagrPartTrackCellsPool
      real(kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThickness
      integer :: timeLev
      integer, pointer :: filterNum
      type (field2DReal), pointer :: uVV, vVV, wVV
      real (kind=RKIND), dimension(:,:), pointer :: potentialDensity

      !write(stderrUnit,*) 'inialize_vertex_velocity start'

      call mpas_pool_get_config(ocnConfigs, 'config_AM_lagrPartTrack_filter_number', filterNum)

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! setup pointers / get block
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackFields', lagrPartTrackFieldsPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackScratch', lagrPartTrackScratchPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

        !! initialize field for RBF (needed depending upon calling pattern of analysis member)
        !call mpas_initialize_vectors(meshPool)
        !call mpas_init_reconstruct(meshPool)

        ! initialize fresh memory for both levels
        do timeLev = 1, 2
          ! connect variables to pointer array
          call mpas_pool_get_field(lagrPartTrackFieldsPool, 'uVertexVelocity', uVV, timeLevel=timeLev)
          call mpas_pool_get_field(lagrPartTrackFieldsPool, 'vVertexVelocity', vVV, timeLevel=timeLev)
          call mpas_pool_get_field(lagrPartTrackFieldsPool, 'wVertexVelocity', wVV, timeLevel=timeLev)

          ! initialize, but could potentially remove these lines
          uVV % array = 0.0_RKIND
          vVV % array = 0.0_RKIND
          wVV % array = 0.0_RKIND

         ! make sure memory has been allocated
#ifdef MPAS_DEBUG
         if(.not.associated(uVV) .or. &
            .not.associated(vVV) .or. &
            .not.associated(wVV)) then
            write(stderrUnit,*) '[u,v,w]VertexVelocity memory not allocated!'
          end if
#endif

          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLev)
          call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel=timeLev)

          ! get new time level velocity (linear RBF)
          !write(stderrUnit,*) 'uVV= ', uVV % array
#ifdef MPAS_DEBUG
          call mpas_timer_start("init_reconst_filter_LPT", .false., timerReconstFilter_init)
#endif
          call ocn_vertex_reconstruction(filterNum, meshPool, lagrPartTrackScratchPool, lagrPartTrackCellsPool, &
            layerThickness, normalVelocity, uVV, vVV, wVV)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("init_reconst_filter_LPT", timerReconstFilter_init)
#endif

        end do

        block => block %  next
      end do
      !write(stderrUnit,*) 'inialize_vertex_velocity end'

   end subroutine initalize_fields!}}}

!***********************************************************************
!
!  routine initalize_wachspress_coefficients
!
!> \brief   Initialize Wachspress coefficients
!> \author  Phillip Wolfram
!> \date    01/26/2015
!> \details
!>  This routine inializes the B_i Wachspress coefficients which are
!>  static in time
!
!-----------------------------------------------------------------------
   subroutine intialize_wachspress_coefficients(domain, err) !{{{

      implicit none

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: lagrPartTrackCellsPool
      type (mpas_pool_type), pointer :: meshPool
      integer  :: nVertices, iCell, i, im1, i0, ip1, iVertex
      integer, pointer :: nCells
      real (kind=RKIND), dimension(:), allocatable :: xv,yv,zv
      real (kind=RKIND), pointer :: radiusLocal
      integer, dimension(:,:), pointer :: verticesOnCell
      integer, dimension(:), pointer :: nCellVerticesArray
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), pointer :: areaBArray
      logical, pointer :: on_a_sphere, is_periodic
      real(kind=RKIND), pointer :: x_period, y_period

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! setup pointers / get block
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
        call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
        call mpas_pool_get_config(meshPool, 'x_period', x_period)
        call mpas_pool_get_config(meshPool, 'y_period', y_period)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nCellVerticesArray)
        call mpas_pool_get_config(meshPool, 'sphere_radius', radiusLocal)
        call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
        call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
        call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'zCell', zCell)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'wachspressAreaB', areaBArray)

        ! compute B_i coefficients
        areaBArray = 0.0_RKIND
        do iCell = 1, nCells
          nVertices = nCellVerticesArray(iCell)
          allocate(xv(nVertices), yv(nVertices), zv(nVertices))
          if (on_a_sphere .or. .not. is_periodic) then
            xv = xVertex(verticesOnCell(:,iCell))
            yv = yVertex(verticesOnCell(:,iCell))
            zv = zVertex(verticesOnCell(:,iCell))
          else
            do iVertex=1,nVertices
              xv(iVertex) = mpas_fix_periodicity(xVertex(verticesOnCell(iVertex,iCell)), &
                xCell(iCell), x_period)
              yv(iVertex) = mpas_fix_periodicity(yVertex(verticesOnCell(iVertex,iCell)), &
                yCell(iCell), y_period)
              zv(iVertex) = zVertex(verticesOnCell(iVertex,iCell))
            end do
          end if
          do i = 1, nVertices
            ! compute first area B_i
            ! get vertex indices
            im1 = mod(nVertices + i - 2, nVertices) + 1
            i0  = mod(nVertices + i - 1, nVertices) + 1
            ip1 = mod(nVertices + i    , nVertices) + 1

            ! precompute B_i areas
            ! always the same because B_i independent of xp,yp,zp
            areaBArray(iCell, i) = mpas_triangle_signed_area( (/ xv(im1),yv(im1),zv(im1) /) , &
              (/ xv(i0),yv(i0),zv(i0) /) , &
              (/ xv(ip1),yv(ip1),zv(ip1) /) , meshPool)
          end do
          deallocate(xv, yv, zv)

        end do


        block => block %  next
      end do

   end subroutine intialize_wachspress_coefficients !}}}

!***********************************************************************
!
!  routine compute_velocity_on_potentialdensity_surface
!
!> \brief   compute_velocity_on_potentialdensity_surface
!> \author  Phillip Wolfram
!> \date    09/15/2014
!> \details
!>  This routine interpolates the velocity field onto the potential
!>  density surface
!
!-----------------------------------------------------------------------
  subroutine compute_velocity_on_potentialdensity_surface(domain,err,aTimeLevel) !{{{

      implicit none

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain
      integer, intent(in) :: aTimeLevel

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool, diagnosticsPool, lagrPartTrackCellsPool
      real (kind=RKIND), dimension(:,:), pointer :: zonVel, merVel, depth, normalVelocityMer, normalVelocityZon
      real (kind=RKIND), dimension(:), pointer :: buoyancySurfaceValues
      integer, pointer :: nBuoyancySurfaces, nCells
      real (kind=RKIND), dimension(:,:), pointer :: buoyancy, zMid
      integer :: iLevel, aBuoyancySurface, iCell
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND) :: phiInterp                       !< location to interpolate
      real (kind=RKIND) :: alpha
      integer :: iHigh, iLow, aval
      real (kind=RKIND) :: eps=1e-14_RKIND

      !write(stderrUnit,*) 'compute_velocity_on_potentialdensity_surface start'

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! get pools
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        ! connect variables to pointer array
        call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', normalVelocityMer)
        call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', normalVelocityZon)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'buoyancySurfaceVelocityZonal', zonVel)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'buoyancySurfaceVelocityMeridional', merVel)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'buoyancySurfaceDepth', depth)
        call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', buoyancy)
        call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nBuoyancySurfaces', nBuoyancySurfaces)
        call mpas_pool_get_array(lagrPartTrackCellsPool,'buoyancySurfaceValues', buoyancySurfaceValues)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        !write(stderrUnit,*) 'nBuoyancySurfaces  =', nBuoyancySurfaces
        !write(stderrUnit,*) 'buoyancySurfaceValues=', buoyancySurfaceValues
        !write(stderrUnit,*) 'size(buoyancy)=',size(buoyancy)
        !write(stderrUnit,*) 'shape(buoyancy)=',shape(buoyancy)
        zonVel = -9999.0_RKIND
        merVel = -9999.0_RKIND
        ! for each buoyancy surface
        do aBuoyancySurface = 1, nBuoyancySurfaces
          ! for each cell
          do iCell = 1, nCells

            phiInterp = buoyancySurfaceValues(aBuoyancySurface)

            ! get correct vertical levels

#ifdef MPAS_DEBUG
            !call mpas_timer_start("mpas_get_vertical_idLPT", .false., timerVerticalID)
#endif
            iLevel = mpas_get_vertical_id(maxLevelCell(iCell), phiInterp, buoyancy(:,iCell))
#ifdef MPAS_DEBUG
            !call mpas_timer_stop("mpas_get_vertical_idLPT", timerVerticalID)
#endif

            if(iLevel < 1) then
              ! top level
              if (iLevel == 0) then
                aval = maxloc(buoyancy(1:maxLevelCell(iCell),iCell),1)
                ! bottom level
              else if (iLevel == -1) then
                aval = minloc(buoyancy(1:maxLevelCell(iCell),iCell),1)
              end if
              zonVel(aBuoyancySurface, iCell) = normalVelocityZon(aval,iCell)
              merVel(aBuoyancySurface, iCell) = normalVelocityMer(aval,iCell)
              depth(aBuoyancySurface, iCell)  = zMid(aval,iCell)
            else
              ! perform the interpolation
              call get_bounding_indices(iLow, iHigh, phiInterp, buoyancy(:,iCell), iLevel, maxLevelCell(iCell))
              ! get alpha between points
              if(abs(buoyancy(iHigh,iCell) - buoyancy(iLow,iCell)) < eps) then
                ! we really can't distinguish between each of these points numerically, just take the
                ! average of both
                alpha = 0.5_RKIND
              else
                alpha = (phiInterp - buoyancy(iLow,iCell))/(buoyancy(iHigh,iCell) - buoyancy(iLow,iCell))
              end if

              ! interpolate to the correct surface
              zonVel(aBuoyancySurface, iCell) = alpha * normalVelocityZon(iHigh, iCell) + &
                                  (1.0_RKIND - alpha) * normalVelocityZon(iLow, iCell)
              merVel(aBuoyancySurface, iCell) = alpha * normalVelocityMer(iHigh, iCell) + &
                                  (1.0_RKIND - alpha) * normalVelocityMer(iLow, iCell)
              depth(aBuoyancySurface, iCell)  = alpha * zMid(iHigh, iCell) + &
                                  (1.0_RKIND - alpha) * zMid(iLow, iCell)
            end if

          end do

        end do

        block => block %  next
      end do
      !write(stderrUnit,*) 'compute_velocity_on_potentialdensity_surface end'

   end subroutine compute_velocity_on_potentialdensity_surface !}}}

!***********************************************************************
!
!  routine initialize_particle_properties
!
!> \brief   Initialize particle properties
!> \author  Phillip Wolfram
!> \date    09/24/2014
!> \details
!>  This routine initializes particle data for
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
   subroutine initialize_particle_properties(domain, timeLevel, err)!{{{

      implicit none
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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: lagrPartTrackFieldsPool, lagrPartTrackCellsPool, lagrPartTrackScratchPool
      integer, dimension(:), pointer :: cellOwnerBlock
      integer, pointer :: currentBlock, ioBlock, indexToParticleID, transfered
      integer :: currentProc, ioProc
      type (mpas_particle_list_type), pointer :: particlelist
      type (mpas_particle_type), pointer :: particle

      real (kind=RKIND), dimension(3) :: particlePosition, particleVelocity
      real (kind=RKIND), pointer :: xParticle, yParticle, zParticle, lonVel, latVel, buoyancyParticle, sumU, sumV, sumUU, sumUV, sumVV
      real (kind=RKIND), dimension(3) :: xSubStep, diffSubStep, diffParticlePosition
      real (kind=RKIND), pointer :: zLevelParticle
      real (kind=RKIND), dimension(:,:), pointer :: zTop, vertVelocityTop, zMid, areaBArray
      real (kind=RKIND), dimension(:), pointer :: bottomDepth
      type (field2DReal), pointer :: normalVelocity, uVertexVelocity, vVertexVelocity, wVertexVelocity, layerThickness

      real (kind=RKIND), dimension(:,:), pointer :: uVertexVelocityArray, vVertexVelocityArray, wVertexVelocityArray, buoyancy, buoyancyTimeInterp, potentialDensity

      real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex

      integer, dimension(:,:), pointer :: verticesOnCell, boundaryVertex
      integer, dimension(:), pointer :: maxLevelCell

      integer theVertex, iLevel, iLevelBuoyancy, aVertex, &
        nSteps, timeStep, subStep, subStepOrder, timeInterpOrder, aTimeLevel, nCellVertices, blockProc, arrayIndex
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: nCellVerticesArray
      integer, dimension(:,:), pointer :: cellsOnCell
      logical, dimension(:,:), pointer :: ioProcRecvList
      logical, dimension(:), pointer :: ioProcSendList
      logical, pointer :: onSphere

      real (kind=RKIND), dimension(4) :: kWeightK, kWeightKVert
      real (kind=RKIND), dimension(3,4) :: kWeightX
      real (kind=RKIND), dimension(4) :: kWeightXVert
      real (kind=RKIND), dimension(4) :: kWeightT, kWeightTVert
      real (kind=RKIND), dimension(3,5) :: kCoeff
      real (kind=RKIND), dimension(5) :: kCoeffVert
      real (kind=RKIND), dimension(2) :: timeCoeff
      real (kind=RKIND) :: dt, dtSim, tSubStep
      real (kind=RKIND), pointer :: dtParticle
      real (kind=RKIND) :: zSubStep
      real (kind=RKIND) :: diffSubStepVert, diffParticlePositionVert, particleVelocityVert, verticalVelocityInterp
      real (kind=RKIND) :: buoyancyInterp
      real (kind=RKIND), pointer :: sphereRadius
      integer, pointer :: verticalTreatment, vertexReconstMethod, timeIntegration, indexLevel, filterNum, iCell

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block)) !{{{
        ! allocate scratch memory / setup pointers / get block
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackFields', lagrPartTrackFieldsPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackScratch', lagrPartTrackScratchPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)

        ! particlelist should be stored in the structs pool probably (need to seriously rework the code!)
        particlelist => block % particlelist

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'zCell', zCell)
        call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
        call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
        call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'wachspressAreaB', areaBArray)
        call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
        call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
        call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
        call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', buoyancyTimeInterp)

        call mpas_pool_get_field(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
        call mpas_dmpar_exch_halo_field(normalVelocity)
        call mpas_pool_get_field(statePool, 'layerThickness', layerThickness, timeLevel=timeLevel)

        call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'boundaryVertex', boundaryVertex)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nCellVerticesArray)

        call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
        call mpas_pool_get_config(meshPool, 'sphere_radius', sphereRadius)

        !!!!!!!!!! LOOP OVER PARTICLES !!!!!!!!!!
        ! update the particle position (just from initialized value for now)
        ! this is a loop over particle list and its datastructures
        do while(associated(particlelist)) !{{{
          ! get pointers / option values
          particle => particlelist % particle

          ! get values {{{
          call mpas_pool_get_array(particle % haloDataPool, 'xParticle', xParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'yParticle', yParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'zParticle', zParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'currentCell', iCell)
          particlePosition(1) = xParticle
          particlePosition(2) = yParticle
          particlePosition(3) = zParticle

          call mpas_pool_get_array(particle % haloDataPool, 'zLevelParticle', zLevelParticle)

          call mpas_pool_get_array(particle % haloDataPool, 'verticalTreatment', verticalTreatment)
          call mpas_pool_get_array(particle % haloDataPool, 'vertexReconstMethod', vertexReconstMethod)
          call mpas_pool_get_array(particle % haloDataPool, 'indexLevel', indexLevel)
          call mpas_pool_get_array(particle % haloDataPool, 'timeIntegration', timeIntegration)
          call mpas_pool_get_array(particle % haloDataPool, 'dtParticle', dtParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'buoyancyParticle', buoyancyParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'indexToParticleID', indexToParticleID)

          call mpas_pool_get_array(particle % haloDataPool, 'lonVel', lonVel)
          call mpas_pool_get_array(particle % haloDataPool, 'latVel', latVel)
          call mpas_pool_get_array(particle % haloDataPool, 'sumU', sumU)
          call mpas_pool_get_array(particle % haloDataPool, 'sumV', sumV)
          call mpas_pool_get_array(particle % haloDataPool, 'sumUU', sumUU)
          call mpas_pool_get_array(particle % haloDataPool, 'sumUV', sumUV)
          call mpas_pool_get_array(particle % haloDataPool, 'sumVV', sumVV)

          !}}}

          !!!!!!!!!! PERFORM SAMPLING (VELOCITY, TEMP, SALINITY, ETC) !!!!!!!!!!
          ! this could be done just before the output anyway

          ! need iCell computed for final position
          ! need scalar values interpolated in time to yield single value
          ! probably need to store zMid too, including flipping it
          ! get updated cell location
#ifdef MPAS_DEBUG
          call mpas_timer_start("get_validated_cell_idLPT", .false., timerValidatedCell_init)
          write(stderrUnit,*) 'sampling initialization'
#endif
          call get_validated_cell_id(nCells, xCell,yCell,zCell , xVertex,yVertex,zVertex, &
            particlePosition(1),particlePosition(2),particlePosition(3), meshPool, &
            nCellVerticesArray, verticesOnCell, iCell, nCellVertices, cellsOnCell)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("get_validated_cell_idLPT", timerValidatedCell_init)
#endif

          if(verticalTreatment == 4) then  !('buoyancySurface') !{{{
            iLevelBuoyancy = mpas_get_vertical_id(maxLevelCell(iCell), buoyancyParticle, buoyancyTimeInterp(:,iCell))

            ! interpolate the scalars now (assumes that scalar value is constant within a particular cell)
            call interp_cell_scalars(iLevelBuoyancy, maxLevelCell(iCell), buoyancyParticle, buoyancyTimeInterp(:,iCell), &
              zMid(:,iCell), zLevelParticle)
          else
            ! make sure final zLevelParticle is ok so that it can't extent past zMid range

#ifdef MPAS_DEBUG
            !call mpas_timer_start("mpas_get_vertical_idLPT", .false., timerVerticalID)
#endif
            iLevel = mpas_get_vertical_id(maxLevelCell(iCell), zLevelParticle, zMid(:,iCell))
#ifdef MPAS_DEBUG
            !call mpas_timer_stop("mpas_get_vertical_idLPT", timerVerticalID)
#endif
          end if

          ! compute necessary information for autocorrelation !{{{
          ! get the updated velocity
          ! ensure that buoyancy is fixed for each run
          buoyancyInterp = buoyancyParticle
          ! we just need the last part of the velocity field interpolation because we are at the end of the timestep
          timeInterpOrder = 1
          timeCoeff(1) = 1.0_RKIND
          timeCoeff(2) = 0.0_RKIND
          ! return interpolated horizontal velocity "particleVelocity" and vertical velocity "particleVelocityVert"
          ! noting we use the final positions
          call velocity_time_interpolation(particleVelocity, particleVelocityVert, &
            diagnosticsPool, lagrPartTrackFieldsPool, &
            timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
            verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
            particlePosition, zLevelParticle, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex,  &
            meshPool, areaBArray)
          ! convert horizontal velocity to lat/lon velocity

          !write(stderrUnit,*) iLevel, particleVelocity

          ! store velocity for use in computing normalized autocorrelation offline
#ifdef MPAS_DEBUG
          !call mpas_timer_start("mpas_convert_xyz_velocity_to_latlonLPT", .false., timerConvertXYZLatLon)
#endif

          call mpas_convert_xyz_velocity_to_latlon(lonVel, latVel, particlePosition, particleVelocity)

#ifdef MPAS_DEBUG
          !call mpas_timer_stop("mpas_convert_xyz_velocity_to_latlonLPT", timerConvertXYZLatLon)
#endif

          !}}}

          ! get next particle to process on the list
          particlelist => particlelist % next
        end do !}}}

        ! get next block
        block => block % next !}}}
      end do !}}}

   end subroutine initialize_particle_properties !}}}

!***********************************************************************
!
!  routine particle_vertical_treatment
!
!> \brief   Vertical treatment to obtain correct horizontal velocity field
!> \author  Phillip Wolfram
!> \date    03/31/2014
!> \details
!>  This routine returns the vertex values which will be used in the
!>  Wachspress interoplant (uvCell) based on
!>  vertex velocities uVertexVelocity, vVertexVelocity, wVertexVelocity
!>  for a given cell which has nCellVertices which are determined from
!>  the list verticesOnCell.
!>  The routine collapses the vertical to a scalar.
!
!-----------------------------------------------------------------------
   subroutine particle_vertical_treatment(verticalTreatment, indexLevel, nCellVertices, verticesOnCell, & !{{{
       uVertexVelocity, vVertexVelocity, wVertexVelocity, &
       uvCell, boundaryVertex, iLevel, nVertLevels, &
       zLoc, zMid, zTop, phiInterp, phiMid, vertVelocityTop, vertVelocityInterp)

      implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:,:), intent(in) :: &
       uVertexVelocity, vVertexVelocity, wVertexVelocity              !< vertex velocities
     integer, dimension(:), intent(in) :: verticesOnCell              !< list of vertex indices on cell
     integer, intent(in) :: nCellVertices                             !< current cell and the number of cell vertices
     integer, intent(in) :: iLevel                                    !< vertical level / cell of zLoc
     integer, intent(in) :: nVertLevels                               !< number of vertical levels
     integer, intent(in) :: verticalTreatment                         !< vertical treatment encoded as int
     integer, intent(in) :: indexLevel                                !< value of index for fixed index space
     integer, dimension(:), intent(in) :: boundaryVertex              !< boundary vertices for particular level
     real (kind=RKIND), intent(in) :: zLoc                            !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: zTop              !< elevation of cell top
     real (kind=RKIND), dimension(:), intent(in) :: zMid              !< elevation of cell middle
     real (kind=RKIND), intent(in) :: phiInterp                       !< buoyancy value to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: phiMid            !< buoyancy values at cell mid points
     real (kind=RKIND), dimension(:), intent(in) :: vertVelocityTop   !< velocity at top of cell

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:,:), intent(out) :: uvCell  !< components of vertex velocity (vertically selected)
     real (kind=RKIND), intent(out) :: vertVelocityInterp                         ! vertically interpolated velocity

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     integer :: aVertex, theVertex

     vertVelocityInterp = 0.0_RKIND

     verticalTreatmentCase: select case (verticalTreatment)

       case (1) verticalTreatmentCase !('indexLevel') !{{{
         ! get vertically interpolanted values for vertexes in cell
         ! and form polygon vertex values
         do aVertex = 1, nCellVertices
           theVertex = verticesOnCell(aVertex)
           ! assume that we only care about the top level for a surface drifter
           ! assumes that velocity is constant in top half of the cell
           uvCell(1,aVertex) = uVertexVelocity(indexLevel,theVertex)
           uvCell(2,aVertex) = vVertexVelocity(indexLevel,theVertex)
           uvCell(3,aVertex) = wVertexVelocity(indexLevel,theVertex)
         end do

         !! ensure that the boundary condition is enforced
         !call zero_boundary_nodal_values(nCellVertices, verticesOnCell, &
         !  boundaryVertex, uvUcell, uvVcell, uvWcell)

         ! no vertical motion, just horizontal motion
         return
         !}}}

       case (2) verticalTreatmentCase !('fixedZLevel') !{{{

         ! interpolate the horizontal velocity based on z-levels
         call interp_nodal_vectors(ncellvertices, verticesoncell, &
           ilevel, nVertLevels, zLoc, zmid,  &
           uvertexvelocity, vvertexvelocity, wvertexvelocity, uvCell)

         !! ensure that there is zero nodal velocity on the boundary
         !call zero_boundary_nodal_values(nCellVertices, verticesOnCell, &
         !  boundaryVertex, uvUcell, uvVcell, uvWcell)

         ! no vertical motion, just horizontal motion
         return
         !}}}

       case (3) verticalTreatmentCase !('passiveFloat') !{{{

         ! interpolate the vertical velocity
         vertVelocityInterp = interp_vert_velocity_to_zlevel( &
           iLevel, zLoc, zTop, vertVelocityTop)

         ! interpolate the horizontal velocity based on z-levels
         call interp_nodal_vectors(ncellvertices, verticesoncell, &
           ilevel, nVertLevels, zLoc, zmid, &
           uvertexvelocity, vvertexvelocity, wvertexvelocity, uvCell)

         !! ensure that there is zero nodal velocity on the boundary
         !call zero_boundary_nodal_values(nCellVertices, verticesOnCell, &
         !  boundaryVertex, uvUcell, uvVcell, uvWcell)

         return
         !}}}

       case (4) verticalTreatmentCase !('buoyancySurface') !{{{
         ! no vertical velocity required because there is not vertical integration for position

         ! interpolate the horizontal velocity
         call interp_nodal_vectors(ncellvertices, verticesoncell, &
           ilevel, nVertLevels, phiInterp, phiMid, &
           uvertexvelocity, vvertexvelocity, wvertexvelocity, uvCell)

         ! ensure that there is zero nodal velocity on the boundary
         !call zero_boundary_nodal_values(nCellVertices, verticesOnCell, &
         !  boundaryVertex, uvUcell, uvVcell, uvWcell)

         ! no vertical motion, just horizontal motion
         return
         !}}}

       case (5) verticalTreatmentCase !('argoFloat')  !{{{


         !}}}

       case default verticalTreatmentCase !{{{
         write(stderrUnit,*) 'Vertical treatment for particle integration unknwon (', verticalTreatment, ')!'
         return
         !}}}

     end select verticalTreatmentCase

   end subroutine particle_vertical_treatment!}}}

!***********************************************************************
!
!  routine velocity_time_interpolation
!
!> \brief   Compute velocity interpolations, including in time
!> \author  Phillip Wolfram
!> \date    09/12/2014
!> \details
!>  This routine interpolates velocity in time and space to a particular
!>  location xSubStep
!
!-----------------------------------------------------------------------
  subroutine velocity_time_interpolation(particleVelocity, particleVelocityVert, &
      diagnosticsPool,  lagrPartTrackFieldsPool, &
      timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
      verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
      xSubStep, zSubStep, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex, meshPool, areaBArray) !{{{
    implicit none

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (mpas_pool_type), pointer, intent(in) :: diagnosticsPool, lagrPartTrackFieldsPool
    integer, intent(in) :: timeInterpOrder
    real (kind=RKIND), dimension(2), intent(in) :: timeCoeff
    integer, intent(in) :: iCell
    integer, dimension(:), intent(in) :: maxLevelCell
    integer, intent(in) :: verticalTreatment
    integer, pointer, intent(in) :: indexLevel
    integer, intent(in) :: nCellVertices
    integer, dimension(:,:), intent(in) :: verticesOnCell
    integer, dimension(:,:), intent(in) :: boundaryVertex
    real (kind=RKIND), intent(in) :: zSubStep
    real (kind=RKIND), dimension(:,:), intent(in) :: zMid
    real (kind=RKIND), dimension(:,:), intent(in) :: zTop
    real (kind=RKIND), dimension(:,:), intent(in) :: vertVelocityTop
    real (kind=RKIND), dimension(:,:), intent(in) :: areaBArray
    real (kind=RKIND), dimension(:), intent(in) :: xVertex, yVertex, zVertex
    real (kind=RKIND), dimension(3), intent(in) :: xSubStep
    type (mpas_pool_type), pointer, intent(in) :: meshPool

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------
    real (kind=RKIND), intent(in) :: buoyancyInterp
    integer, intent(inout) :: iLevel

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------
    real (kind=RKIND), dimension(3), intent(out) :: particleVelocity
    real (kind=RKIND), intent(out) :: particleVelocityVert

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------
    integer :: aVertex, aTimeLevel
    real (kind=RKIND), dimension(:,:), pointer :: uVertexVelocityArray, vVertexVelocityArray, wVertexVelocityArray, buoyancy
    real (kind=RKIND) :: verticalVelocityInterp
    real(kind=RKIND), dimension(:), allocatable :: areaB
    real(kind=RKIND), dimension(:,:), allocatable :: vertCoords
    real(kind=RKIND), dimension(:,:), allocatable :: uvCell
    logical, pointer :: on_a_sphere, is_periodic
    real(kind=RKIND), pointer :: x_period, y_period
#ifdef MPAS_DEBUG
    call mpas_timer_start("velocity_time_interpolationLPT", .false., timerVelTimeInterp)
#endif

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
    call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
    call mpas_pool_get_config(meshPool, 'x_period', x_period)
    call mpas_pool_get_config(meshPool, 'y_period', y_period)

    ! allocations for particular cell !{{{
    allocate(vertCoords(3,nCellVertices), uvCell(3,nCellVertices), areaB(nCellVertices))
    !}}}

    ! get horizontal vertex locations (noting that there may be a
    ! bit of error because the particle could be at the top
    ! of the cell or at the bottom of the cell)
    do aVertex = 1, nCellVertices
      if (on_a_sphere .or. .not. is_periodic) then
        vertCoords(1,aVertex) = xVertex(verticesOnCell(aVertex,iCell))
        vertCoords(2,aVertex) = yVertex(verticesOnCell(aVertex,iCell))
        vertCoords(3,aVertex) = zVertex(verticesOnCell(aVertex,iCell))
      else
        vertCoords(1,aVertex) = mpas_fix_periodicity(xVertex(verticesOnCell(aVertex,iCell)), xSubStep(1), x_period)
        vertCoords(2,aVertex) = mpas_fix_periodicity(yVertex(verticesOnCell(aVertex,iCell)), xSubStep(2), y_period)
        vertCoords(3,aVertex) = zVertex(verticesOnCell(aVertex,iCell))
      end if
      areaB(aVertex) = areaBArray(iCell, aVertex)
    end do

    ! initialize velocities to 0
    particleVelocity = 0.0_RKIND
    particleVelocityVert = 0.0_RKIND

    ! general interpolation for the velocity field
    do aTimeLevel = 1, timeInterpOrder

      ! define arrays!{{{
      call mpas_pool_get_array(lagrPartTrackFieldsPool, 'uVertexVelocity', uVertexVelocityArray, timeLevel=aTimeLevel)
      call mpas_pool_get_array(lagrPartTrackFieldsPool, 'vVertexVelocity', vVertexVelocityArray, timeLevel=aTimeLevel)
      call mpas_pool_get_array(lagrPartTrackFieldsPool, 'wVertexVelocity', wVertexVelocityArray, timeLevel=aTimeLevel)
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', buoyancy)
      !}}}

      ! get final, interpolated particle velocity at this point (collapse to point)
      if (verticalTreatment == 4) then
        ! buoyancy case (not using zSubStep for interpolation / iLevel)
        ! use existing code noting we need to flip the order to get the right iLevel

#ifdef MPAS_DEBUG
        !call mpas_timer_start("mpas_get_vertical_idLPT", .false., timerVerticalID)
#endif
        iLevel = mpas_get_vertical_id(maxLevelCell(iCell), buoyancyInterp, buoyancy(:,iCell))
#ifdef MPAS_DEBUG
        !call mpas_timer_stop("mpas_get_vertical_idLPT", timerVerticalID)
        write(stderrUnit,*) 'iLevel=',iLevel
#endif
        ! note, if buoyancyInterp out of range this will try to reorient the particle to the top / bottom but there
        ! will definitely be some error with this type of computation because the buoyancy is not available at this location
        ! the time interpolation, as a consequency, can mix velocities from different buoyancy surfaces in order to advect
        ! the particle
        !write(stderrUnit,*) 'iLevel = ', iLevel
      end if

      call particle_vertical_treatment(verticalTreatment, indexLevel, nCellVertices, verticesOnCell(:,iCell), &
        uVertexVelocityArray, vVertexVelocityArray, wVertexVelocityArray, uvCell, boundaryVertex(iLevel,:), &
        iLevel, maxLevelCell(iCell), zSubStep, zMid(:,iCell), zTop(:,iCell), buoyancyInterp, buoyancy(:,iCell), &
        vertVelocityTop(:,iCell), verticalVelocityInterp)

      ! vertical
      particleVelocityVert = particleVelocityVert + &
        timeCoeff(aTimeLevel) * verticalVelocityInterp

      ! horizontal
      ! timer commented out because it is used in more than just compute...
#ifdef MPAS_DEBUG
      call mpas_timer_start("part_horiz_interpLPT", .false., timerHorizVelInterp)
      write(stderrUnit,*) 'particleVelocityVert=',particleVelocityVert
#endif
      particleVelocity = particleVelocity + &
        timeCoeff(aTimeLevel) * particle_horizontal_interpolation(nCellVertices, vertCoords, &
        xSubStep, uvCell, meshPool, areaB)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("part_horiz_interpLPT", timerHorizVelInterp)
      write(stderrUnit,*) 'particleVelocity=',particleVelocity
#endif


    end do

    ! deallocations of temp memory
    deallocate(vertCoords, uvCell, areaB)

#ifdef MPAS_DEBUG
    call mpas_timer_stop("velocity_time_interpolationLPT", timerVelTimeInterp)
#endif

  end subroutine velocity_time_interpolation !}}}

  subroutine zero_autocorrelation_sums(domain) !{{{
    implicit none

    ! input/output variables
    type (domain_type), intent(inout) :: domain
    ! local
    type (block_type), pointer :: block
    type (mpas_particle_list_type), pointer :: particlelist
    type (mpas_particle_type), pointer :: particle
    ! output variables (per particle)
    real (kind=RKIND), pointer :: sumU, sumV, sumUU, sumUV, sumVV
    integer, pointer :: currentCell

    ! get the appropriate pools
      block => domain % blocklist
      do while (associated(block)) !{{{
        particlelist => block % particlelist
        do while(associated(particlelist)) !{{{
          ! get pointers / option values
          particle => particlelist % particle

          ! get values (may want a flag for reinitialization in the future)
          !call mpas_pool_get_array(particle % haloDataPool, 'sumU', sumU)
          !call mpas_pool_get_array(particle % haloDataPool, 'sumV', sumV)
          !call mpas_pool_get_array(particle % haloDataPool, 'sumUU', sumUU)
          !call mpas_pool_get_array(particle % haloDataPool, 'sumUV', sumUV)
          !call mpas_pool_get_array(particle % haloDataPool, 'sumVV', sumVV)
          call mpas_pool_get_array(particle % haloDataPool, 'currentCell', currentCell)

          ! initialize the values
          !sumU = 0.0_RKIND
          !sumV = 0.0_RKIND
          !sumUU = 0.0_RKIND
          !sumUV = 0.0_RKIND
          !sumVV = 0.0_RKIND
          currentCell = -1

          ! get next particle to process on the list
          particlelist => particlelist % next
        end do !}}}

        ! get next block
        block => block % next
      end do !}}}

  end subroutine zero_autocorrelation_sums !}}}

!***********************************************************************
!
!  routine particle_horizontal_interpolation
!
!> \brief   Horizontal treatment to obtain correct velocity field at point
!> \author  Phillip Wolfram
!> \date    03/31/2014
!> \details
!>  This routine returns the point values which will be used in the
!>  particle interpolation time integration based on
!>  vertex velocities uVertexVelocity, vVertexVelocity, wVertexVelocity
!
!-----------------------------------------------------------------------
  function particle_horizontal_interpolation(nCellVertices, vertCoords, & !{{{
       pointInterp, uVertex, meshPool, areaB)

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     integer, intent(in) :: nCellVertices
     real (kind=RKIND), dimension(3, nCellVertices), intent(in) :: vertCoords
     real (kind=RKIND), dimension(3), intent(in) :: pointInterp
     real (kind=RKIND), dimension(3, nCellVertices), intent(in) :: uVertex
     real (kind=RKIND), dimension(nCellVertices), intent(in) :: areaB
     type (mpas_pool_type), pointer :: meshPool

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(nCellVertices) :: lambda

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(3) :: particle_horizontal_interpolation


     ! get lambda coordinate for particle
     lambda = mpas_wachspress_coordinates(nCellVertices, vertCoords , &
       pointInterp, meshPool, areaB)
!#ifdef MPAS_DEBUG
!write(stderrUnit,*) 'lambda=',lambda
!write(stderrUnit,*) 'uVertex=',uVertex(1,:)
!write(stderrUnit,*) 'vVertex=',uVertex(2,:)
!write(stderrUnit,*) 'wVertex=',uVertex(3,:)
!#endif

     ! update particle velocities via horizontal interpolation
     particle_horizontal_interpolation(1) = mpas_wachspress_interpolate(lambda, uVertex(1,:))
     particle_horizontal_interpolation(2) = mpas_wachspress_interpolate(lambda, uVertex(2,:))
     particle_horizontal_interpolation(3) = mpas_wachspress_interpolate(lambda, uVertex(3,:))

   end function particle_horizontal_interpolation !}}}

!***********************************************************************
!
!  routine particle_horizontal_movement
!
!> \brief   Compute horizontal movement for particle so particle stays
!>          in spherical shell
!> \author  Phillip Wolfram
!> \date    05/20/2014
!> \details
!>  This routine returns the particle position pParticle corresponding
!>  to an initial particle position pParticle for a Cartesian movemnt
!>  dpParticle.  If the calculation is onSphere, then the distance
!>  |dpParticle| must be along the great circle route of pParticle
!>  and the projection of pParticle + dpParticle on the spherical
!>  shell corresponding to pParticle.
!
!-----------------------------------------------------------------------
  subroutine particle_horizontal_movement(meshPool, pParticle, dpParticle) !{{{

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:), intent(in) :: dpParticle
     type (mpas_pool_type), intent(in), pointer :: meshPool

     !-----------------------------------------------------------------
     ! input / output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:), intent(inout) :: pParticle

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND) :: lenPath, arcLen
     real (kind=RKIND) :: radiusShell
     real (kind=RKIND), dimension(size(pParticle)) :: pParticleTemp
     real (kind=RKIND), dimension(size(pParticle)) :: pParticleInterp
     real (kind=RKIND) :: alpha
     real (kind=RKIND), parameter :: eps=1e-10_RKIND
     logical, pointer :: onSphere, is_periodic
     real(kind=RKIND), pointer :: x_period, y_period
     ! choosen based on the parameters, note that we loose about 6 - 7 units of precision because R is so large!
     ! therefore, eps = 1e-10 is conservative, if not too high!  this just helps with numerical stability
     !dpParticle    =   -4.2428037617887103E-011   4.3076544298828060E-011   5.0760704444480953E-011
     !pParticle     =    4444887.2990309987       -891565.00525021972        4476665.3916420965
     !pParticleTemp =    4444887.2990309987       -891565.00525021972        4476665.3916420965
     !mpas_arc_length =    0.0000000000000000      lenPath =    7.8945399869363434E-011

     call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
     call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
     call mpas_pool_get_config(meshPool, 'x_period', x_period)
     call mpas_pool_get_config(meshPool, 'y_period', y_period)

     ! may need a condition to determine if we need to project back to the sphere
     if(onSphere) then
       ! need to make sure new point is on the spherical shell

       ! get path length
       !write(stderrUnit, *) 'dpParticle = ', dpParticle
       lenPath = sqrt(sum(dpParticle*dpParticle))

       ! consider case of particle not moving (need to have this code here in general)
       !if (lenPath < eps) then
       ! this is ok because this is only the case if the points are the same.  If there is a
       ! numerical instability it probably should be handled differently.
       if (lenPath < eps) then
         return
       end if

       ! get radius of particle's horizontal shell
       radiusShell = sqrt(sum(pParticle*pParticle))

       ! project endpoint to spherical shell containing pParticle
       pParticleTemp = pParticle + dpParticle
       pParticleTemp = ( radiusShell / sqrt(sum(pParticleTemp*pParticleTemp)) ) * pParticleTemp

       ! compute alpha parameter for spherical interpolant / extrapolant
       !write(stderrUnit,*) 'dpParticle    = ', dpParticle
       !write(stderrUnit,*) 'pParticle     = ', pParticle
       !write(stderrUnit,*) 'pParticleTemp = ', pParticleTemp
       !write(stderrUnit,*) 'mpas_arc_length = ',  mpas_arc_length(pParticle(1),pParticle(2),pParticle(3) , pParticleTemp(1), pParticleTemp(2), pParticleTemp(3)), 'lenPath = ', lenPath
       arcLen = mpas_arc_length(pParticle(1),pParticle(2),pParticle(3) , pParticleTemp(1),pParticleTemp(2),pParticleTemp(3))
       if (arcLen > eps) then
         alpha = lenPath / arcLen
       else
         return
       endif

       ! compute final position based on spherical interpolant
       call mpas_spherical_linear_interp(pParticleInterp, pParticle, pParticleTemp, alpha)
       pParticle = pParticleInterp
     else
       ! we are just on a plane so there is no need for spherical interpolation to keep
       ! the new particle location on a spherical shell
       pParticle = pParticle + dpParticle

       ! periodic fix to make sure particle advection stays in domain
       if (is_periodic) then
         pParticle(1) = mpas_fix_periodicity(pParticle(1), x_period/2.0_RKIND, x_period)
         pParticle(2) = mpas_fix_periodicity(pParticle(2), y_period/2.0_RKIND, y_period)
         !pParticle(3) = pParticle(3)
       end if
     endif

   end subroutine particle_horizontal_movement!}}}

!***********************************************************************
!
!  routine interp_cell_scalars
!
!> \brief   Interpolate cell scalar vector based on a criteria (z-level, buoyancy, etc)
!> \author  Phillip Wolfram
!> \date    06/18/2014
!> \details
!>  This routine interpolates cell scalar vector to a particular scalar value
!>  depending upon a criteria such as z-level, buoyancy, etc.
!
!-----------------------------------------------------------------------
   subroutine interp_cell_scalars(iLevel, nVertLevels, zInterp, zVals, & !{{{
       phiVals, phiInterp)

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:), intent(in) :: zVals             !< scalar values (x) on cell for interpolant
     integer, intent(in) :: iLevel                                    !< vertical level / cell of phiInterp
     integer, intent(in) :: nVertLevels                               !< number of vertical levels
     real (kind=RKIND), intent(in) :: zInterp                       !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex velocities are defined)

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), intent(out) :: phiInterp        !< interpolated cell scalar

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND) :: alpha
     integer :: aVertex, theVertex, iHigh, iLow
     real (kind=RKIND) :: eps=1e-14

     if(iLevel < 1) then
       ! top level
       if (iLevel == 0) then
         phiInterp = phiVals(nVertLevels)
         ! bottom level
       else if (iLevel == -1) then
         phiInterp = phiVals(1)
       end if
     else
       call get_bounding_indices(iLow, iHigh, zInterp, zVals, iLevel, nVertLevels)

       ! interpolate to vertical level now
       if(abs(zVals(iHigh) - zVals(iLow)) < eps) then
         ! we really can't distinguish between each of these points numerically, just take the
         ! average of both
         alpha = 0.5_RKIND
       else
         ! interpolate to vertical level now
         alpha = (zInterp - zVals(iLow))/(zVals(iHigh) - zVals(iLow))
       end if

       ! interpolate to the vertical level
       phiInterp = alpha * phiVals(iHigh) + (1.0_RKIND - alpha) * phiVals(iLow)
     end if

   end subroutine interp_cell_scalars!}}}

!***********************************************************************
!
!  routine interp_nodal_scalars
!
!> \brief   Interpolate nodal scalar vector based on a criteria (z-level, buoyancy, etc)
!> \author  Phillip Wolfram
!> \date    05/27/2014
!> \details
!>  This routine interpolates nodal scalar vector to a particular scalar value
!>  depending upon a criteria such as z-level, buoyancy, etc.
!
!-----------------------------------------------------------------------
   subroutine interp_nodal_scalars(nCellVertices, verticesOnCell, & !{{{
       iLevel, nVertLevels, phiInterp, phiVals, &
       scalarVec, vertexScalar)

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:,:), intent(in) :: scalarVec       !< vertex scalar
     integer, dimension(:), intent(in) :: verticesOnCell              !< list of vertex indices on cell
     integer, intent(in) :: nCellVertices                             !< number of cell vertices
     integer, intent(in) :: iLevel                                    !< vertical level / cell of phiInterp
     integer, intent(in) :: nVertLevels                               !< number of vertical levels
     real (kind=RKIND), intent(in) :: phiInterp                       !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex velocities are defined)

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:), intent(out) :: vertexScalar      !< components of vertex scalar (interpolated)

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND) :: alpha
     integer :: aVertex, theVertex, iHigh, iLow

     call get_bounding_indices(iLow, iHigh, phiInterp, phiVals, iLevel, nVertLevels)

     ! interpolate to vertical level now
     alpha = (phiInterp - phiVals(iLow))/(phiVals(iHigh) - phiVals(iLow))

     ! interpolate to the vertical level
     do aVertex = 1, nCellVertices
       theVertex = verticesOnCell(aVertex)
       ! assume for now that we only care about the top level for a surface drifter
       vertexScalar(aVertex) = alpha * scalarVec(iHigh, theVertex) + (1.0_RKIND - alpha) * scalarVec(iLow, theVertex)
     end do

   end subroutine interp_nodal_scalars!}}}

!***********************************************************************
!
!  routine interp_nodal_vectors
!
!> \brief   Interpolate nodal vector to scalar based on a criteria (z-level, buoyancy, etc)
!> \author  Phillip Wolfram
!> \date    05/27/2014
!> \details
!>  This routine interpolates nodal vectors to a particular scalar value
!>  depending upon a criteria such as z-level, buoyancy, etc.
!
!-----------------------------------------------------------------------
   subroutine interp_nodal_vectors(nCellVertices, verticesOnCell, & !{{{
       iLevel, nVertLevels, phiInterp, phiVals, &
       uVertexVelocity, vVertexVelocity, wVertexVelocity, uvCell)

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:,:), intent(in) :: &
       uVertexVelocity, vVertexVelocity, wVertexVelocity              !< vertex velocities
     integer, dimension(:), intent(in) :: verticesOnCell              !< list of vertex indices on cell
     integer, intent(in) :: nCellVertices                             !< number of cell vertices
     integer, intent(in) :: iLevel                                    !< vertical level / cell of phiInterp
     integer, intent(in) :: nVertLevels                               !< number of vertical levels
     real (kind=RKIND), intent(in) :: phiInterp                       !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex velocities are defined)

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:,:), intent(out) :: uvCell !< components of vertex velocity (vertically selected)

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND) :: alpha
     integer :: aVertex, theVertex, iHigh, iLow
     real (kind=RKIND) :: eps=1e-14_RKIND

     uvCell = 0.0_RKIND

     !write(stderrUnit,*) 'interp'
     if(iLevel < 1) then
       if(iLevel == 0) then
         !write(stderrUnit,*) 'nCellVertices = ', nCellVertices, 'buoyancyInterp= ', phiInterp
         do aVertex = 1, nCellVertices
           theVertex = verticesOnCell(aVertex)
           !write(stderrUnit,*) maxloc(phiVals,1), phiVals(1:nVertLevels), uVertexVelocity(1:nVertLevels,theVertex)
           uvCell(1,aVertex) = uVertexVelocity(maxloc(phiVals(1:nVertLevels),1), theVertex)
           uvCell(2,aVertex) = vVertexVelocity(maxloc(phiVals(1:nVertLevels),1), theVertex)
           uvCell(3,aVertex) = wVertexVelocity(maxloc(phiVals(1:nVertLevels),1), theVertex)
         end do
       else if(iLevel == -1) then
         do aVertex = 1, nCellVertices
           theVertex = verticesOnCell(aVertex)
           uvCell(1,aVertex) = uVertexVelocity(minloc(phiVals(1:nVertLevels),1), theVertex)
           uvCell(2,aVertex) = vVertexVelocity(minloc(phiVals(1:nVertLevels),1), theVertex)
           uvCell(3,aVertex) = wVertexVelocity(minloc(phiVals(1:nVertLevels),1), theVertex)
         end do
       end if
     else
       call get_bounding_indices(iLow, iHigh, phiInterp, phiVals, iLevel, nVertLevels)

       ! interpolate to vertical level now
       if(abs(phiVals(iHigh) - phiVals(iLow)) < eps) then
         ! we really can't distinguish between each of these points numerically, just take the
         ! average of both
         alpha = 0.5_RKIND
       else
         alpha = (phiInterp - phiVals(iLow))/(phiVals(iHigh) - phiVals(iLow))
       end if

       ! interpolate to the vertical level
       do aVertex = 1, nCellVertices
         theVertex = verticesOnCell(aVertex)
         ! assume for now that we only care about the top level for a surface drifter
         uvCell(1,aVertex) = alpha * uVertexVelocity(iHigh, theVertex) + &
           (1.0_RKIND - alpha) * uVertexVelocity(iLow, theVertex)
         uvCell(2,aVertex) = alpha * vVertexVelocity(iHigh, theVertex) + &
           (1.0_RKIND - alpha) * vVertexVelocity(iLow, theVertex)
         uvCell(3,aVertex) = alpha * wVertexVelocity(iHigh, theVertex) + &
           (1.0_RKIND - alpha) * wVertexVelocity(iLow, theVertex)
       end do
     end if

   end subroutine interp_nodal_vectors!}}}

!***********************************************************************
!
!  routine zero_boundary_nodal_values
!
!> \brief   Enfore boundary condition for nodal value, setting to 0
!> \author  Phillip Wolfram
!> \date    05/27/2014
!> \details
!>  This routine ensures zero Dirchilet boundary conditions
!>  (commonly for the nodal velocity)
!
!-----------------------------------------------------------------------
   subroutine zero_boundary_nodal_values(nCellVertices, verticesOnCell, & !{{{
       boundaryVertex, uvCell)

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     integer, dimension(:), intent(in) :: verticesOnCell              !< list of vertex indices on cell
     integer, intent(in) :: nCellVertices                             !< number of cell vertices
     integer, dimension(:), intent(in) :: boundaryVertex              !< boundary vertices for particular level

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:,:), intent(out) :: uvCell  !< components of vertex velocity (vertically selected)

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     integer :: aVertex, theVertex

     ! make sure to mask all boundary vertexes to be zero to enforce boundary conditions
     do aVertex = 1, nCellVertices
       theVertex = verticesOnCell(aVertex)
       !write(stderrUnit,*) boundaryVertex(theVertex)
       ! make all the boundary values be zero to prevent particle from horizontally leaving cell
       uvCell(1,aVertex) = uvCell(1,aVertex) * (1-boundaryVertex(theVertex))
       uvCell(2,aVertex) = uvCell(2,aVertex) * (1-boundaryVertex(theVertex))
       uvCell(3,aVertex) = uvCell(3,aVertex) * (1-boundaryVertex(theVertex))
     end do

   end subroutine zero_boundary_nodal_values!}}}

!***********************************************************************
!
!  routine get_bounding_indices
!
!> \brief   Get indices for high and low values for interpolation
!> \author  Phillip Wolfram
!> \date    05/28/2014
!> \details
!>  This routine returns the indices (iLow, iHigh) on either side of phiInterp
!
!-----------------------------------------------------------------------
   subroutine get_bounding_indices(iLow, iHigh, phiInterp, phiVals, iLevel, nVertLevels) !{{{

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     integer, intent(in) :: iLevel                                    !< vertical level / cell of phiInterp
     integer, intent(in) :: nVertLevels                               !< number of vertical levels
     real (kind=RKIND), intent(in) :: phiInterp                       !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex velocities are defined)

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     integer, intent(out) :: iLow, iHigh                              !< interpolation indices

     ! assumes increasing index is in decreasing phi space, and phiInterp in range of phiVals
     !write(stderrUnit,*) 'phiInterp = ', phiInterp, 'phiVals = ', phiVals, 'nVertLevels = ', nVertLevels, 'iLevel = ', iLevel
     if(phiInterp > phiVals(iLevel)) then
       iHigh = iLevel + 1
       iLow  = iLevel
     else
       iHigh = iLevel
       iLow  = iLevel + 1
     end if

     ! check to make sure points are in range
     ! optimization point: smarter algorithm won't have to call this ever
     if(.not.((phiInterp <= phiVals(iHigh)) .and. (phiInterp >= phiVals(iLow)))) then
       !write(stderrUnit,*) 'fast interpolation failed, trying general, brute force search for interpolation bounds'
       !write(stderrUnit,*) 'iLow = ', iLow, ' iHigh = ', iHigh
       !write(stderrUnit,*) 'phiInterp = ', phiInterp , ' phiLow =',  phiVals(iLow), ' phiHigh = ', phiVals(iHigh)
       call get_bounding_indices_brute_force(nVertLevels, phiInterp, phiVals, iLow, iHigh)
     end if

   end subroutine get_bounding_indices !}}}

!***********************************************************************
!
!  routine get_bounding_indices_brute_force
!
!> \brief   Get the interpolation bounds via brute force
!> \author  Phillip Wolfram
!> \date    05/27/2014
!> \details
!>  This routine finds the interpolation bounds directly (brute force).
!
!-----------------------------------------------------------------------
   subroutine get_bounding_indices_brute_force(nVertLevels, phiInterp, phiVals, iLow, iHigh) !{{{

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     integer, intent(in) :: nVertLevels                               !< number of vertical levels
     real (kind=RKIND), intent(in) :: phiInterp                       !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex velocities are defined)

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     integer, intent(out) :: iLow, iHigh                          !< indexes for the high and low components for the interpolant

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     integer :: aLevel

     ! make no assumptions
     do aLevel = 1, nVertLevels-1
       if(phiVals(aLevel) <= phiInterp .and. phiInterp <= phiVals(aLevel+1)) then
         iLow  = aLevel
         iHigh = aLevel+1
         exit
       else if (phiVals(aLevel+1) <= phiInterp .and. phiInterp <= phiVals(aLevel)) then
         iLow  = aLevel+1
         iHigh = aLevel
         exit
       end if
     end do

#ifdef MPAS_DEBUG
     if(phiVals(iLow) <= phiInterp .and. phiInterp <= phiVals(iHigh)) then
       ! we are ok
       !write(stderrUnit,*) 'brute force interpolation successful'
       !write(stderrUnit,*) 'iLow = ', iLow, ' iHigh = ', iHigh
       !write(stderrUnit,*) 'phiInterp = ', phiInterp , ' phiLow =',  phiVals(iLow), ' phiHigh = ', phiVals(iHigh)
     else
       write(stderrUnit,*) 'brute force interpolation failed with phiInterp = ', phiInterp, ' phiLow = ', phiVals(iLow), ' phiHigh = ', phiVals(iHigh)
       !write(stderrUnit,*) ' phiVals = ', phiVals(1:nVertLevels)
     end if

     write(stderrUnit,*) 'Warning!: brute force interpolation used, boundary condition may be wrong!'
#endif

   end subroutine get_bounding_indices_brute_force!}}}

!***********************************************************************
!
!  routine interp_vert_velocity_to_zlevel
!
!> \brief   Interpolate the vertical velcity to a z level
!> \author  Phillip Wolfram
!> \date    05/08/2014
!> \details
!>  This routine interpolates the vertical velocity to a particular
!>  z-level.
!
!-----------------------------------------------------------------------
   real (kind=RKIND) function interp_vert_velocity_to_zlevel( & !{{{
       iLevel, zSubStep, zTop, vertVelocityTop)

      implicit none
     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     integer, intent(in) :: iLevel                                    !< vertical level / cell of zSubStep
     real (kind=RKIND), intent(in) :: zSubStep                        !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: zTop              !< elevation of cell top
     real (kind=RKIND), dimension(:), intent(in) :: vertVelocityTop   !< vertical velocity at top of cell

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     !real (kind=RKIND), intent(out) :: interp_vert_velocity_to_zlevel

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND) :: alpha

     if(iLevel < 1) then
       if(iLevel == 0) then
         interp_vert_velocity_to_zlevel = vertVelocityTop(maxloc(zTop,1))
       else if(iLevel == -1) then
         interp_vert_velocity_to_zlevel = vertVelocityTop(minloc(zTop,1))
       end if
     else
       ! interpolate the velocity [assumes zTop(iLevel+1) <= zSubStep <= zTop(iLevel)]
       if (zTop(iLevel+1) <= zSubStep .and. zSubStep <= zTop(iLevel)) then
         alpha = (zSubStep - zTop(iLevel+1))/(zTop(iLevel)- zTop(iLevel+1))
       else if (zTop(iLevel) <= zSubStep .and. zSubStep <= zTop(iLevel+1)) then
         alpha = (zSubStep - zTop(iLevel))/(zTop(iLevel+1)- zTop(iLevel))
#ifdef MPAS_DEBUG
       else
         write(stderrUnit,*) 'Error with vertical velocity interpolation!'
#endif
       end if
       interp_vert_velocity_to_zlevel = alpha * vertVelocityTop(iLevel)+ (1.0_RKIND - alpha) * vertVelocityTop(iLevel+1)
     end if

    end function interp_vert_velocity_to_zlevel!}}}

!***********************************************************************
!
!  routine time_interp_field
!
!> \brief   Interpolate a field in time
!> \author  Phillip Wolfram
!> \date    07/16/2014
!> \details
!>  This routine interpolates a field in time over multiple levels.
!
!-----------------------------------------------------------------------
   subroutine time_interp_field(basePool, timeInterpOrder, timeCoeff, field, fieldname) !{{{
     implicit none

     type (mpas_pool_type), pointer, intent(in) :: basePool
     integer, intent(in) :: timeInterpOrder
     real (kind=RKIND), dimension(:), intent(in) :: timeCoeff
     real (kind=RKIND), dimension(:,:), pointer, intent(out) :: field
     character(len=*), intent(in) :: fieldname

     real (kind=RKIND), dimension(:,:), pointer :: tempfield
     integer :: aTimeLevel

     field = 0.0_RKIND
     do aTimeLevel = 1, timeInterpOrder
       call mpas_pool_get_array(basePool, trim(fieldname), tempfield, timeLevel=aTimeLevel)
       field = field + timeCoeff(aTimeLevel) * tempfield
     end do

   end subroutine time_interp_field !}}}
!}}}

end module ocn_lagrangian_particle_tracking

! vim: foldmethod=marker
