! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
#define COMMA ,
#ifdef MPAS_DEBUG
#define LIGHT_DEBUG_WRITE(Z) !  print *, Z
#define LIGHT_DEBUG_ALL_WRITE(Z) ! print *, Z
#else
#define LIGHT_DEBUG_WRITE(Z)  ! print *, Z
#define LIGHT_DEBUG_ALL_WRITE(Z) ! print *, Z
#endif

#define LIGHT_ERROR_WRITE(Z)  print *, Z
#define LIGHT_WARNING_WRITE(Z)  print *, Z

module ocn_lagrangian_particle_tracking

   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables

   use ocn_particle_list
   use ocn_lagrangian_particle_tracking_interpolations
   use ocn_lagrangian_particle_tracking_reset

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

   ! neighboring processors and arrays for send numbers / recvs, total number of neighboring processors to a particular processor
   ! allocated in init, deallocated in finalize
   ! these globals could be moved to the framework component, as well
   ! as quite a few of the subroutines
   integer, dimension(:), pointer :: g_compProcNeighsNearby => null(), g_compProcNeighs => null(), g_ioProcNeighs=>null()

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
      logical :: resetParticles

      err = 0

      LIGHT_DEBUG_WRITE('starting ocn_init_lagrangian_particle_tracking')
      !call mpas_timer_start("totalLPT")
      call mpas_timer_start("initLPT")

      ! resets likely are unnecessary
      !call mpas_stream_mgr_reset_alarms(stream_manager, streamID='lagrPartTrackInput', ierr=err)
      !call mpas_stream_mgr_reset_alarms(stream_manager, streamID='lagrPartTrackRestart', ierr=err)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! build up the particles lists and fill them with their data
      call mpas_particle_list_build_and_assign_particle_list(domain, err)
      ! now we have built the particlelist for a given block
      LIGHT_DEBUG_WRITE('finished building and allocating particles in init')

      ! parallel code:

      ! get "MPI halos" for communication of particles in halo during computational step
      call mpas_particle_list_build_computation_halos(domain, err, g_compProcNeighsNearby)
      LIGHT_DEBUG_WRITE('finished building and computational halos')

      ! get "MPI halos" for IO communication during write and restart steps (ioBlock to currentBlocks)
      ! AllToAll Computation!
      call mpas_particle_list_build_halos(domain, err, 'currentBlock', g_ioProcNeighs)
      LIGHT_DEBUG_WRITE('g_ioProcNeighs=' COMMA g_ioProcNeighs)
      LIGHT_DEBUG_WRITE('finished building io halos')

      ! transfer particles to their appropriate blocks (currentBlock) via MPI
      ! note, don't necessarily need to have g_ionSend and g_ionRecv comeout
#ifdef MPAS_DEBUG
      call mpas_timer_start("trans_from_block_to_blockLPT")
      call mpas_particle_list_test_numparticles_to_neighprocs(domain % dminfo % my_proc_id, g_compProcNeighsNearby, g_ioProcNeighs)
#endif
      ! move particles to the 'currentBlock'
      call mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, .True., .False., 'currentBlock', &
        g_ioProcNeighs)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("trans_from_block_to_blockLPT")
      !call MPI_Barrier(domain % dminfo % comm, err)
#endif

      ! tests to make sure all the values are ok !{{{
#ifdef MPAS_DEBUG
      call mpas_particle_list_test_neighscalc(domain, err)
      call mpas_particle_list_test_numparticles_to_neighprocs(domain % dminfo % my_proc_id, g_compProcNeighsNearby, g_ioProcNeighs)
      call mpas_particle_list_test_num_current_particlelist(domain)
#endif
      LIGHT_DEBUG_WRITE('g_compProcNeighsNearby = ' COMMA g_compProcNeighsNearby)
      !}}}

!      ! now set sums for autocorrelation calculation to be zero
!      call zero_autocorrelation_sums(domain)

      ! previous compute startup calls
      call intialize_wachspress_coefficients(domain, err)
      call initalize_fields(domain, err)
      call compute_velocity_on_potentialdensity_surface(domain,err,1)
      call compute_velocity_on_potentialdensity_surface(domain,err,2)
      call initialize_particle_properties(domain,2,err)
      call write_lagrangian_particle_tracking(domain, err)

      ! set up particle reset condition
      call ocn_setup_particle_reset_condition(domain, err)

      LIGHT_DEBUG_WRITE('finished ocn_init_lagrangian_particle_tracking')
      call mpas_timer_stop("initLPT")

      !call mpas_timer_stop("totalLPT")

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
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: lagrPartTrackFieldsPool, lagrPartTrackCellsPool, lagrPartTrackScratchPool
      integer, dimension(:), pointer :: cellOwnerBlock
      integer, pointer :: currentBlock, ioBlock, indexToParticleID, transfered, currentCellGlobalID
      type (mpas_particle_list_type), pointer :: particlelist
      type (mpas_particle_type), pointer :: particle

      real (kind=RKIND), dimension(3) :: particlePosition, particleVelocity
      real (kind=RKIND), pointer :: xParticle, yParticle, zParticle, buoyancyParticle !, &
!        lonVel, latVel, sumU, sumV, sumUU, sumUV, sumVV
      real (kind=RKIND), pointer :: lonParticle, latParticle
      real (kind=RKIND), dimension(3) :: xSubStep, diffSubStep, diffParticlePosition
      real (kind=RKIND), pointer :: zLevelParticle
      real (kind=RKIND), dimension(:,:), pointer :: areaBArray
      real (kind=RKIND), dimension(:), pointer :: bottomDepth, ssh
      type (field2DReal), pointer :: normalVelocity, uVertexVelocity, vVertexVelocity, wVertexVelocity, layerThickness

      real (kind=RKIND), dimension(:,:), pointer :: uVertexVelocityArray, vVertexVelocityArray, wVertexVelocityArray, &
                                                    buoyancyTimeInterp

      real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex

      integer, dimension(:,:), pointer :: verticesOnCell, boundaryVertex
      integer, dimension(:), pointer :: maxLevelCell, indexToCellID

      integer theVertex, iLevel, iLevelBuoyancy, aVertex, &
        nSteps, timeStep, subStep, subStepOrder, timeInterpOrder, aTimeLevel, nCellVertices, blockProc, arrayIndex
      integer, pointer :: nCells, nVertLevels, iCell
      integer, dimension(:), pointer :: nCellVerticesArray
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:,:), pointer :: cellsOnVertex
      logical :: resetParticles

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
      integer, pointer :: verticalTreatment, indexLevel
      type (MPAS_timeInterval_type) :: timeStepESMF
      logical :: resetParticle, resetParticleAny
      integer :: err_tmp

      ! sampling variables
      ! active tracers
      integer, pointer :: index_temperature, &
                          index_salinity
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      real (kind=RKIND), pointer :: particleTemperature

      real (kind=RKIND), pointer :: particleSalinity

      ! BGC tracers
      integer, pointer :: index_DIC, &
                          index_ALK, &
                          index_PO4, &
                          index_NO3, &
                          index_SiO3, &
                          index_NH4, &
                          index_Fe, &
                          index_O2
      real (kind=RKIND), dimension(:,:,:), pointer :: ecosysTracers

      ! DIC
      real (kind=RKIND), pointer :: particleDIC

      ! ALK
      real (kind=RKIND), pointer :: particleALK

      ! PO4
      real (kind=RKIND), pointer :: particlePO4

      ! NO3
      real (kind=RKIND), pointer :: particleNO3

      ! SiO3
      real (kind=RKIND), pointer :: particleSiO3

      ! NH4
      real (kind=RKIND), pointer :: particleNH4

      ! Fe
      real (kind=RKIND), pointer :: particleFe

      ! O2
      real (kind=RKIND), pointer :: particleO2

      REAL (kind=RKIND), parameter :: eps=1e-12


      err = 0

      !! don't do compute for debugging purposes
      !call mpas_log_write( 'Computing Lagrangian Particle Tracking -- NO COMPUTATION!')
      !return

      dminfo = domain % dminfo

      LIGHT_DEBUG_WRITE('Computing Lagrangian Particle Tracking...')
      !call mpas_timer_start("totalLPT")
      call mpas_timer_start("computeLPT")

      if (trim(config_AM_lagrPartTrack_reset_criteria) == 'none') then
        resetParticles = .False.
      else
        resetParticles = .True.
      end if

      ! get the most recent velocities on the potential density surfaces
#ifdef MPAS_DEBUG
      call mpas_timer_start("velocity_pot_density_LPT")
#endif
      call compute_velocity_on_potentialdensity_surface(domain,err,2)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("velocity_pot_density_LPT")
#endif

      resetParticleAny = .False.

      block => domain % blocklist
      do while (associated(block)) !{{{
#ifdef MPAS_DEBUG
        call mpas_timer_start("memtasksLPT")
#endif
        ! allocate scratch memory / setup pointers / get block !{{{
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
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
        buoyancyTimeInterp => potentialDensity

        ! note, originally this was diagnostics % state % normalVelocity (without time level), but
        ! now there is a time level so selection of the correct time level appears to be tricky
        ! the issue is large, nearly NAN normalVelocities with timeLevel=2
        call mpas_pool_get_field(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
        call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel=timeLevel)
        ! potentially a costly function call, but ensures halos are ok
        !call mpas_dmpar_exch_halo_field(normalVelocity)
        call mpas_pool_get_field(statePool, 'layerThickness', layerThickness, timeLevel=timeLevel)

        ! sampling variables
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        if (config_use_ecosysTracers) then
          call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
          call mpas_pool_get_dimension(tracersPool, 'index_DIC', index_DIC)
          call mpas_pool_get_dimension(tracersPool, 'index_ALK', index_ALK)
          call mpas_pool_get_dimension(tracersPool, 'index_PO4', index_PO4)
          call mpas_pool_get_dimension(tracersPool, 'index_NO3', index_NO3)
          call mpas_pool_get_dimension(tracersPool, 'index_SiO3', index_SiO3)
          call mpas_pool_get_dimension(tracersPool, 'index_NH4', index_NH4)
          call mpas_pool_get_dimension(tracersPool, 'index_Fe', index_Fe)
          call mpas_pool_get_dimension(tracersPool, 'index_O2', index_O2)
        end if

        call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'boundaryVertex', boundaryVertex)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nCellVerticesArray)
        call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
        call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
        call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
        !}}}

        ! super-cycling:
        ! get the appropriate time scale for particle integration, e.g., dtSim
        if (trim(config_AM_lagrPartTrack_compute_interval) == 'dt') then
          ! from time step value
          ! would eventually need to correspond to domain dt, but for now this is a
          ! global constant
          call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
          ! using the output interval
        else if (trim(config_AM_lagrPartTrack_compute_interval) == 'output') then
          LIGHT_ERROR_WRITE('LIGHT output mode is not yet supported!')
          !! from output interval
          !call MPAS_stream_mgr_get_property(domain % streamManager, config_AM_lagrPartTrack_output_stream, &
          !                                  MPAS_STREAM_PROPERTY_REF_TIME, config_dt, err)
          !write(*,*) 'config_AM_lagrPartTrack_output_stream = ', config_AM_lagrPartTrack_output_stream
          !write(*,*) 'config_dt= ', config_dt
          !call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
        else
          ! using specified interval from timestring
          call mpas_set_timeInterval(timeStepESMF, timeString=trim(config_AM_lagrPartTrack_compute_interval), ierr=err)
        end if
        call mpas_get_timeInterval(timeStepESMF, dt=dtSim)
        LIGHT_DEBUG_WRITE('dtSim = ' COMMA dtSim)

        !{{{
        ! flip previous time level to back (switching memory pointers)
        call mpas_pool_shift_time_levels(lagrPartTrackFieldsPool)
        call mpas_pool_get_field(lagrPartTrackFieldsPool, 'uVertexVelocity', uVertexVelocity, timeLevel=2)
        call mpas_pool_get_field(lagrPartTrackFieldsPool, 'vVertexVelocity', vVertexVelocity, timeLevel=2)
        call mpas_pool_get_field(lagrPartTrackFieldsPool, 'wVertexVelocity', wVertexVelocity, timeLevel=2)
#ifdef MPAS_DEBUG
        call mpas_timer_stop("memtasksLPT")
#endif

#ifdef MPAS_DEBUG
        call mpas_timer_start("reconst_filter_LPT")
#endif
        ! need to handle periodicity within functions below for vertex reconstruction
        call ocn_vertex_reconstruction(config_AM_lagrPartTrack_filter_number, meshPool, lagrPartTrackScratchPool, lagrPartTrackCellsPool, &
          layerThickness % array, normalVelocity % array, &
          uVertexVelocity, vVertexVelocity, wVertexVelocity)
#ifdef MPAS_DEBUG
        call mpas_timer_stop("reconst_filter_LPT")
        LIGHT_DEBUG_ALL_WRITE('uVertexVelocity=' COMMA uVertexVelocity % array)
        LIGHT_DEBUG_ALL_WRITE('vVertexVelocity=' COMMA vVertexVelocity % array)
        LIGHT_DEBUG_ALL_WRITE('wVertexVelocity=' COMMA wVertexVelocity % array)
#endif
        !}}}

        LIGHT_DEBUG_WRITE('beginning particle loop with particlelist associated =  ' COMMA associated(particlelist))
#ifdef MPAS_DEBUG
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
          call mpas_timer_start("memtasksLPT")
#endif
          call mpas_pool_get_array(particle % haloDataPool, 'xParticle', xParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'yParticle', yParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'zParticle', zParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'lonParticle', lonParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'latParticle', latParticle)

          ! compute particle lat/lon position
          call convert_latlon_from_xyz(latParticle, lonParticle, xParticle, yParticle, zParticle)

          particlePosition(1) = xParticle
          particlePosition(2) = yParticle
          particlePosition(3) = zParticle

          call mpas_pool_get_array(particle % haloDataPool, 'zLevelParticle', zLevelParticle)

          call mpas_pool_get_array(particle % haloDataPool, 'verticalTreatment', verticalTreatment)
          call mpas_pool_get_array(particle % haloDataPool, 'indexLevel', indexLevel)
          call mpas_pool_get_array(particle % haloDataPool, 'dtParticle', dtParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'buoyancyParticle', buoyancyParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'indexToParticleID', indexToParticleID)
          call mpas_pool_get_array(particle % haloDataPool, 'currentCell', iCell)
          call mpas_pool_get_array(particle % haloDataPool, 'currentCellGlobalID', currentCellGlobalID)

          ! sampling variables
          call mpas_pool_get_array(particle % haloDataPool, 'particleTemperature', particleTemperature)
          call mpas_pool_get_array(particle % haloDataPool, 'particleSalinity', particleSalinity)

          if (config_use_ecosysTracers) then
            call mpas_pool_get_array(particle % haloDataPool, 'particleDIC', particleDIC)
            call mpas_pool_get_array(particle % haloDataPool, 'particleALK', particleALK)
            call mpas_pool_get_array(particle % haloDataPool, 'particlePO4', particlePO4)
            call mpas_pool_get_array(particle % haloDataPool, 'particleNO3', particleNO3)
            call mpas_pool_get_array(particle % haloDataPool, 'particleSiO3', particleSiO3)
            call mpas_pool_get_array(particle % haloDataPool, 'particleNH4', particleNH4)
            call mpas_pool_get_array(particle % haloDataPool, 'particleFe', particleFe)
            call mpas_pool_get_array(particle % haloDataPool, 'particleO2', particleO2)
          end if

!          call mpas_pool_get_array(particle % haloDataPool, 'lonVel', lonVel)
!          call mpas_pool_get_array(particle % haloDataPool, 'latVel', latVel)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumU', sumU)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumV', sumV)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumUU', sumUU)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumUV', sumUV)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumVV', sumVV)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("memtasksLPT")
#endif
          ! process timers / particle reset functions here (need a timer and ability to reset particle to
          ! some set of initial values for the reset
          !}}}

          !!!!!!!!!! COMPUTE TIME STEP INFORMATION !!!!!!!!!!
#ifdef MPAS_DEBUG
          call mpas_timer_start("time_step_LPT")
#endif
          ! adjust time step for consistency with integer number of steps
          nSteps = ceiling(dtSim/dtParticle)
          dt = dtSim/nSteps

          !!!!!!!!!! ASSIGN TEMPORAL INTEGRATION COEFFICIENTS !!!!!!!!!!
          ! kCoeff   is (3,subStepOrder+1)
          ! kWeightX is subStepOrder
          ! kWeightK is subStepOrder
          ! kWeightT is subStepOrder
          select case (config_AM_lagrPartTrack_timeIntegration) !{{{
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
                call mpas_timer_start("particle_horizontal_movementLPT")
#endif
                call particle_horizontal_movement(meshPool, xSubStep, diffSubStep)
#ifdef MPAS_DEBUG
                call mpas_timer_stop("particle_horizontal_movementLPT")
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
              call mpas_timer_start("get_validated_cell_idLPT")
#endif
              LIGHT_DEBUG_WRITE('beginning of substeps')
              call get_validated_cell_id(nCells, xCell,yCell,zCell , xVertex,yVertex,zVertex, &
                xSubStep(1),xSubStep(2),xSubStep(3), meshPool, &
                nCellVerticesArray, verticesOnCell, iCell, nCellVertices, cellsOnCell)
              LIGHT_DEBUG_WRITE('iCell=' COMMA iCell)
#ifdef MPAS_DEBUG
              call mpas_timer_stop("get_validated_cell_idLPT")
#endif

              if(verticalTreatment /= 4) then
                ! other cases using zSubStep for iLevel and vertical interpolation
#ifdef MPAS_DEBUG
                call mpas_timer_start("mpas_get_vertical_idLPT")
#endif
                ! must have inputs into mpas_get_vertical_id be positive to get right top/bottom
                iLevel = mpas_get_vertical_id(maxLevelCell(iCell), -zSubStep, -zMid(:,iCell))
#ifdef MPAS_DEBUG
                call mpas_timer_stop("mpas_get_vertical_idLPT")
#endif
                LIGHT_DEBUG_WRITE('iLevel=' COMMA iLevel)
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
              call mpas_timer_start("velocity_time_interpolationLPT")
#endif
              call velocity_time_interpolation(particleVelocity, particleVelocityVert, &
                lagrPartTrackFieldsPool, &
                timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
                verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
                xSubStep, zSubStep, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex, meshPool, areaBArray)
#ifdef MPAS_DEBUG
              call mpas_timer_stop("velocity_time_interpolationLPT")
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
            call mpas_timer_start("particle_horizontal_movementLPT")
#endif
            call particle_horizontal_movement(meshPool, particlePosition, diffParticlePosition)
#ifdef MPAS_DEBUG
            call mpas_timer_stop("particle_horizontal_movementLPT")
#endif
            ! now can do any vertical movements independent of the horizontal movement
            ! that was just calculated ( probably need to have more output from the vertical_treatment
            ! and aggregate here
            zLevelParticle = zLevelParticle + diffParticlePositionVert

          end do !}}}
#ifdef MPAS_DEBUG
          call mpas_timer_stop("time_step_LPT")
#endif
          !}}}

          !!!!!!!!!! PERFORM SAMPLING (VELOCITY, TEMP, SALINITY, ETC) AND ERROR CHECKING!!!!!!!!!!
          ! this could be done just before the output anyway

          ! need iCell computed for final position
          ! need scalar values interpolated in time to yield single value
          ! probably need to store zMid too, including flipping it
          ! get updated cell location
#ifdef MPAS_DEBUG
          call mpas_timer_start("get_validated_cell_idLPT")
#endif
          LIGHT_DEBUG_WRITE('do sampling')
          call get_validated_cell_id(nCells, xCell,yCell,zCell , xVertex,yVertex,zVertex, &
            particlePosition(1),particlePosition(2),particlePosition(3), meshPool, &
            nCellVerticesArray, verticesOnCell, iCell, nCellVertices, cellsOnCell)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("get_validated_cell_idLPT")
#endif

          call get_vertical_level_info(iLevel, zLevelParticle, &
            verticalTreatment, maxLevelCell, iCell, buoyancyInterp, buoyancyTimeInterp, &
            zMid, indexLevel, ssh, bottomDepth, indexToParticleID)

           !! sampling setup (if needed)
           !if (&
           !  config_AM_lagrPartTrack_sample_temperature .or. &
           !  config_AM_lagrPartTrack_sample_salinity) then

           !end if

           ! store global cell ID
           currentCellGlobalID = indexToCellID(iCell)

           if ( associated(activeTracers) ) then
             if (config_AM_lagrPartTrack_sample_temperature) then
              !call mpas_log_write(' testing output for temperature')
              if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                call interp_scalars(particleTemperature, activeTracers(index_temperature,:,:), &
                  verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                  particlePosition(1),particlePosition(2),particlePosition(3), &
                  xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                  yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                  zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                  xCell, yCell, zCell, &
                  meshPool, &
                  cellsOnVertex, &
                  verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                  zLevelParticle, zMid, indexLevel)
              else
                call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                  activeTracers(index_temperature,:,iCell), particleTemperature)
              end if
              !call mpas_log_write('done testing output for temperature')
             end if

             if (config_AM_lagrPartTrack_sample_salinity) then
              !call mpas_log_write(' testing output for salinity')
              if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                call interp_scalars(particleSalinity, activeTracers(index_salinity,:,:), &
                  verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                  particlePosition(1),particlePosition(2),particlePosition(3), &
                  xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                  yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                  zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                  xCell, yCell, zCell, &
                  meshPool, &
                  cellsOnVertex, &
                  verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                  zLevelParticle, zMid, indexLevel)
              else
                call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                  activeTracers(index_salinity,:,iCell), particleSalinity)
              end if
              !call mpas_log_write('done testing output for salinity')
             end if
           end if

           if (config_use_ecosysTracers) then
             call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
             if ( associated(ecosysTracers) ) then
               if (config_AM_lagrPartTrack_sample_DIC) then
                 !call mpas_log_write('testing output for DIC')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleDIC, ecosysTracers(index_DIC,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_DIC,:,iCell), particleDIC)
                 end if
                 !call mpas_log_write('done testing output for DIC')
               end if

               if (config_AM_lagrPartTrack_sample_ALK) then
                 !call mpas_log_write('testing output for ALK')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleALK, ecosysTracers(index_ALK,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_ALK,:,iCell), particleALK)
                 end if
                 !call mpas_log_write('done testing output for ALK')
               end if

               if (config_AM_lagrPartTrack_sample_PO4) then
                 !call mpas_log_write('testing output for PO4')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particlePO4, ecosysTracers(index_PO4,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_PO4,:,iCell), particlePO4)
                 end if
                 !call mpas_log_write('done testing output for PO4')
               end if

               if (config_AM_lagrPartTrack_sample_NO3) then
                 !call mpas_log_write('testing output for NO3')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleNO3, ecosysTracers(index_NO3,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_NO3,:,iCell), particleNO3)
                 end if
                 !call mpas_log_write('done testing output for NO3')
               end if

               if (config_AM_lagrPartTrack_sample_SiO3) then
                 !call mpas_log_write('testing output for SiO3')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleSiO3, ecosysTracers(index_SiO3,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_SiO3,:,iCell), particleSiO3)
                 end if
                 !call mpas_log_write('done testing output for SiO3')
               end if

               if (config_AM_lagrPartTrack_sample_NH4) then
                 !call mpas_log_write('testing output for NH4')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleNH4, ecosysTracers(index_NH4,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_NH4,:,iCell), particleNH4)
                 end if
                 !call mpas_log_write('done testing output for NH4')
               end if

               if (config_AM_lagrPartTrack_sample_Fe) then
                 !call mpas_log_write('testing output for Fe')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleFe, ecosysTracers(index_Fe,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_Fe,:,iCell), particleFe)
                 end if
                 !call mpas_log_write('done testing output for Fe')
               end if

               if (config_AM_lagrPartTrack_sample_O2) then
                 !call mpas_log_write('testing output for O2')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleO2, ecosysTracers(index_O2,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_O2,:,iCell), particleO2)
                 end if
                 !call mpas_log_write('done testing output for O2')
               end if
             end if
           end if

!          ! compute necessary information for autocorrelation !{{{
!          ! get the updated velocity
!          ! ensure that buoyancy is fixed for each run
!          buoyancyInterp = buoyancyParticle
!          ! we just need the last part of the velocity field interpolation because we are at the end of the timestep
!          timeInterpOrder = 2
!          timeCoeff(1) = 0.0_RKIND
!          timeCoeff(2) = 1.0_RKIND
!          ! return interpolated horizontal velocity "particleVelocity" and vertical velocity "particleVelocityVert"
!          ! noting we use the final positions
!#ifdef MPAS_DEBUG
!          call mpas_timer_start("velocity_time_interpolationLPT")
!#endif
!          call velocity_time_interpolation(particleVelocity, particleVelocityVert, &
!            lagrPartTrackFieldsPool, &
!            timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
!            verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
!            particlePosition, zLevelParticle, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex, &
!            meshPool, areaBArray)
!#ifdef MPAS_DEBUG
!          call mpas_timer_stop("velocity_time_interpolationLPT")
!#endif
          ! convert horizontal velocity to lat/lon velocity

!          ! store velocity for use in computing normalized autocorrelation offline
!#ifdef MPAS_DEBUG
!          call mpas_timer_start("mpas_convert_xyz_velocity_to_latlonLPT")
!#endif
!          call mpas_convert_xyz_velocity_to_latlon(lonVel, latVel, particlePosition, particleVelocity)
!#ifdef MPAS_DEBUG
!          call mpas_timer_stop("mpas_convert_xyz_velocity_to_latlonLPT")
!#endif

!          ! now store components needed to compute integral timescale
!#ifdef MPAS_DEBUG
!          call mpas_timer_start("storeSingleParticleStats")
!#endif
!          sumU = sumU + lonVel
!          sumV = sumV + latVel
!          sumUU = sumUU + lonVel*lonVel
!          sumUV = sumUV + lonVel*latVel
!          sumVV = sumVV + latVel*latVel
!#ifdef MPAS_DEBUG
!          call mpas_timer_stop("storeSingleParticleStats")
!#endif
          !}}}

          ! properly store particle position (because we can't store arrays directly for particles and must
          ! work entirely in vectors
          xParticle = particlePosition(1)
          yParticle = particlePosition(2)
          zParticle = particlePosition(3)

          ! compute particle lat/lon position
          call convert_latlon_from_xyz(latParticle, lonParticle, xParticle, yParticle, zParticle)

          !!!!!!!!!! PASS PARTICLES FROM PROCESSOR TO PROCESSOR !!!!!!!!!!
          ! note, assumes that integration above didn't move particle past halo
          !{{{
          ! 1. determine if iCell is on halo (just set each particle's currentBlock to the correct currentBlock
          ! 2. determine owning block in halo, update particle's currentBlock

          ! determine currentBlock ownership of iCell
          ! and set cellOwnerBlock to be current block
#ifdef MPAS_DEBUG
          call mpas_timer_start("particleAssignments")
#endif

          resetParticle = .False.
          if (resetParticles) then   ! need to link to currentBlockReset for communication
            ! determine if particles should be reset based on different criteria.  If so, reset them.
            call ocn_evaluate_particle_reset_condition(domain, block, particle, dtSim, iCell, resetParticle, err)
          end if
          resetParticleAny = resetParticleAny .or. resetParticle

          if (.not. resetParticle) then
            ! update halo fields for particles moving from adjacent computational halos
            call mpas_particle_list_update_particle_block(domain, block, particle, 'lagrPartTrackCells', iCell)
          end if

#ifdef MPAS_DEBUG
          call mpas_timer_stop("particleAssignments")
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
      call mpas_timer_start("trans_from_block_to_blockLPT")
#endif
      if (resetParticles .and. resetParticleAny) then
         ! need to link to currentBlockReset for communication
         call mpas_particle_list_build_halos(domain, err, 'currentBlockReset', g_compProcNeighs)
         ! take the union of this halo with the particle computation halos to build complete computational halo
         call mpas_particle_list_self_union_halo_lists(g_compProcNeighs, g_compProcNeighsNearby, domain % dminfo % nprocs, &
                                                       domain % dminfo % my_proc_id)
      else
        allocate(g_compProcNeighs(size(g_compProcNeighsNearby)))
        g_compProcNeighs = g_compProcNeighsNearby
      end if
      call mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, .True., .False., 'currentBlock', &
        g_compProcNeighs)
      deallocate(g_compProcNeighs)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("trans_from_block_to_blockLPT")
#endif

      ! do IO communications if this is an output time step
      if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='lagrPartTrackOutput', &
          direction=MPAS_STREAM_OUTPUT, ierr=err)) then
        call write_lagrangian_particle_tracking(domain, err)
      end if

      call mpas_timer_stop("computeLPT")
      !call mpas_timer_stop("totalLPT")

      LIGHT_DEBUG_WRITE('finished computing Lagrangian Particle Tracking')

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

      !call mpas_timer_start("totalLPT")

      ! do restart if this is a restart step
      if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='lagrPartTrackRestart', &
          direction=MPAS_STREAM_OUTPUT, ierr=err)) then

        call mpas_timer_start("restartLPT")
        LIGHT_DEBUG_WRITE('start ocn_restart_lagrangian_particle_tracking')

        ! note, difference here is that nonhalo data may not be needed for
        ! restart, but including here for now for simplity of the code
        call write_lagrangian_particle_tracking(domain, err)

        LIGHT_DEBUG_WRITE('end ocn_restart_lagrangian_particle_tracking')
        call mpas_timer_stop("restartLPT")
      end if

      !call mpas_timer_stop("totalLPT")

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

      LIGHT_DEBUG_WRITE('start write_lagrangian_particle_tracking')
      !call mpas_timer_start("totalLPT")
      call mpas_timer_start("writeLPT")

      !call mpas_particle_list_test_num_current_particlelist(domain)

      ! transfer particles to their appropriate blocks (ioBlock) via MPI
      ! note, don't necessarily need to have g_ionSend and g_ionRecv comeout
      !call mpas_log_write( 'g_ioProcNeighs = $i',MPAS_LOG_OUT, intArgs=(/  g_ioProcNeighs /) )
#ifdef MPAS_DEBUG
      call mpas_timer_start("trans_from_block_to_blockLPT")
#endif
      ! depreciated (can just use update_halo_io to keep g_ioProcNeighs up to date)
      !! get "MPI halos" for IO communication during write and restart steps (currentBlock to ioBlock)
      call mpas_particle_list_build_halos(domain, err, 'ioBlock', g_ioProcNeighs)
      ! transfer the data
      call mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, .False., .True., 'ioBlock', &
        g_ioProcNeighs)
      deallocate(g_ioProcNeighs)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("trans_from_block_to_blockLPT")
#endif

      ! write out all the data, sorting to make sure that shuffled particles
      ! are ouptut correctly (done separately in each function, could be
      ! pulled out as an optimization)
      call mpas_particle_list_write_halo_data(domain, err)
      ! there is not currently any halo data (don't need to communicate for now)
      !call mpas_particle_list_write_nonhalo_data(domain, err)

#ifdef MPAS_DEBUG
      call mpas_particle_list_test_num_current_particlelist(domain)
#endif
      ! need to now remove the io particles (remove particles that don't have the
      ! correct currentBlock)
      call mpas_particle_list_remove_particles_not_on_current_block(domain,err)

#ifdef MPAS_DEBUG
      call mpas_particle_list_test_num_current_particlelist(domain)
#endif
      LIGHT_DEBUG_WRITE('end write_lagrangian_particle_tracking')
      call mpas_timer_stop("writeLPT")
      !call mpas_timer_stop("totalLPT")


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

      !call mpas_timer_start("totalLPT")
      call mpas_timer_start("finalizeLPT")
      err = 0

      LIGHT_DEBUG_WRITE('start ocn_finalize_lagrangian_particle_tracking')
      block => domain % blocklist
      do while (associated(block))
        call mpas_particle_list_destroy_particle_list(block % particlelist)
        block => block %  next
      end do

      deallocate(g_compProcNeighsNearby)

      LIGHT_DEBUG_WRITE('end ocn_finalize_lagrangian_particle_tracking')
      call mpas_timer_stop("finalizeLPT")
      !call mpas_timer_stop("totalLPT")

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
     character (len=StrKIND) :: message
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
     ! call mpas_log_write( 'max verticesOnCell = $i nVertices = $i',MPAS_LOG_OUT, intArgs=(/ maxval(verticesOnCell(:,iCell)), size(xVertex) /) )
     if (on_a_sphere .or. .not. is_periodic) then
       if(.not. point_in_cell(nCellVertices, &
         xVertex(verticesOnCell(1:nCellVertices,iCell)), &
         yVertex(verticesOnCell(1:nCellVertices,iCell)), &
         zVertex(verticesOnCell(1:nCellVertices,iCell)), &
         xSubStep,ySubStep,zSubStep , on_a_sphere)) then
           write(message, *) 'Point (', xSubStep, ySubStep, zSubStep , ') is horizontally outside cell ', iCell
           LIGHT_DEBUG_WRITE(message)
           write(message, *) 'Cell  (', xCell(iCell), yCell(iCell), zCell(iCell), ') with index ', iCell
           LIGHT_DEBUG_WRITE(message)
           LIGHT_DEBUG_WRITE('xVertex = ' COMMA xVertex(verticesOnCell(1:nCellVertices,iCell)))
           LIGHT_DEBUG_WRITE('yVertex = ' COMMA yVertex(verticesOnCell(1:nCellVertices,iCell)))
           LIGHT_DEBUG_WRITE('zVertex = ' COMMA zVertex(verticesOnCell(1:nCellVertices,iCell)))
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
           write(message, *) 'Point (', xSubStep, ySubStep, zSubStep, ') is horizontally outside cell ', iCell
           LIGHT_DEBUG_WRITE(message)
           write(message, *) 'Cell  (', xCell(iCell), yCell(iCell), zCell(iCell), ') with index ', iCell
           LIGHT_DEBUG_WRITE(message)
           LIGHT_DEBUG_WRITE('xVertex = ' COMMA xtmp)
           LIGHT_DEBUG_WRITE('yVertex = ' COMMA ytmp)
           LIGHT_DEBUG_WRITE('zVertex = ' COMMA zVertex(verticesOnCell(1:nCellVertices,iCell)))
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
      type (mpas_pool_type), pointer :: lagrPartTrackFieldsPool, lagrPartTrackScratchPool, lagrPartTrackCellsPool
      real(kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThickness
      integer :: timeLev
      type (field2DReal), pointer :: uVV, vVV, wVV
      real (kind=RKIND), dimension(:,:), pointer :: potentialDensity

      !call mpas_log_write( 'inialize_vertex_velocity start')

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! setup pointers / get block
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackFields', lagrPartTrackFieldsPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackScratch', lagrPartTrackScratchPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)

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
            LIGHT_DEBUG_WRITE('[u,v,w]VertexVelocity memory not allocated!')
          end if
#endif

          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLev)
          call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel=timeLev)

          ! get new time level velocity (linear RBF)
          !call mpas_log_write( 'uVV= $f',MPAS_LOG_OUT, realArgs=(/ uVV % array /) ) ! This is an array, not sure if log will work.
#ifdef MPAS_DEBUG
          call mpas_timer_start("init_reconst_filter_LPT")
#endif
          call ocn_vertex_reconstruction(config_AM_lagrPartTrack_filter_number, meshPool, &
            lagrPartTrackScratchPool, lagrPartTrackCellsPool, &
            layerThickness, normalVelocity, uVV, vVV, wVV)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("init_reconst_filter_LPT")
#endif

        end do

        block => block %  next
      end do
      !call mpas_log_write( 'inialize_vertex_velocity end')

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
      real (kind=RKIND), dimension(3) :: v1, v2, v3
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
            v1(1) = xv(im1)
            v1(2) = yv(im1)
            v1(3) = zv(im1)
            v2(1) = xv(i0)
            v2(2) = yv(i0)
            v2(3) = zv(i0)
            v3(1) = xv(ip1)
            v3(2) = yv(ip1)
            v3(3) = zv(ip1)
            areaBArray(iCell, i) = mpas_triangle_signed_area(v1, v2, v3, meshPool)
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
      type (mpas_pool_type), pointer :: meshPool, lagrPartTrackCellsPool
      real (kind=RKIND), dimension(:,:), pointer :: zonVel, merVel, depth, normalVelocityMer, normalVelocityZon
      real (kind=RKIND), dimension(:), pointer :: buoyancySurfaceValues
      integer, pointer :: nBuoyancySurfaces, nCells
      real (kind=RKIND), dimension(:,:), pointer :: buoyancy
      integer :: iLevel, aBuoyancySurface, iCell
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND) :: phiInterp                       !< location to interpolate
      real (kind=RKIND) :: alpha
      integer :: iHigh, iLow, aval
      real (kind=RKIND) :: eps=1e-14_RKIND

      LIGHT_DEBUG_WRITE('compute_velocity_on_potentialdensity_surface start')

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! get pools
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)
        ! connect variables to pointer array
        normalVelocityMer => velocityMeridional
        normalVelocityZon => velocityZonal
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'buoyancySurfaceVelocityZonal', zonVel)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'buoyancySurfaceVelocityMeridional', merVel)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'buoyancySurfaceDepth', depth)
        buoyancy => potentialDensity
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nBuoyancySurfaces', nBuoyancySurfaces)
        call mpas_pool_get_array(lagrPartTrackCellsPool,'buoyancySurfaceValues', buoyancySurfaceValues)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        LIGHT_DEBUG_ALL_WRITE('nBuoyancySurfaces  =' COMMA nBuoyancySurfaces)
        LIGHT_DEBUG_ALL_WRITE('buoyancySurfaceValues=' COMMA buoyancySurfaceValues)
        LIGHT_DEBUG_ALL_WRITE('size(buoyancy)=' COMMA size(buoyancy))
        LIGHT_DEBUG_ALL_WRITE('shape(buoyancy)=' COMMA shape(buoyancy))
        zonVel = -9999.0_RKIND
        merVel = -9999.0_RKIND
        ! for each buoyancy surface
        do aBuoyancySurface = 1, nBuoyancySurfaces
          ! for each cell
          do iCell = 1, nCells

            phiInterp = buoyancySurfaceValues(aBuoyancySurface)

            ! get correct vertical levels

            iLevel = mpas_get_vertical_id(maxLevelCell(iCell), phiInterp, buoyancy(:,iCell))

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
      LIGHT_DEBUG_WRITE('compute_velocity_on_potentialdensity_surface end')

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
      type (mpas_pool_type), pointer :: lagrPartTrackFieldsPool, lagrPartTrackCellsPool, lagrPartTrackScratchPool
      integer, dimension(:), pointer :: cellOwnerBlock
      integer, pointer :: currentBlock, ioBlock, indexToParticleID, transfered
      integer :: currentProc, ioProc
      type (mpas_particle_list_type), pointer :: particlelist
      type (mpas_particle_type), pointer :: particle

      real (kind=RKIND), dimension(3) :: particlePosition, particleVelocity
      real (kind=RKIND), pointer :: xParticle, yParticle, zParticle, lonVel, latVel, buoyancyParticle !, &
!        lonVel, latVel, sumU, sumV, sumUU, sumUV, sumVV
      real (kind=RKIND), pointer :: lonParticle, latParticle
      real (kind=RKIND), dimension(3) :: xSubStep, diffSubStep, diffParticlePosition
      real (kind=RKIND), pointer :: zLevelParticle
      real (kind=RKIND), dimension(:,:), pointer :: areaBArray
      real (kind=RKIND), dimension(:), pointer :: bottomDepth
      type (field2DReal), pointer :: normalVelocity, uVertexVelocity, vVertexVelocity, wVertexVelocity, layerThickness

      real (kind=RKIND), dimension(:,:), pointer :: uVertexVelocityArray, vVertexVelocityArray, wVertexVelocityArray, buoyancy, &
                                                    buoyancyTimeInterp

      real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex

      integer, dimension(:,:), pointer :: verticesOnCell, boundaryVertex
      integer, dimension(:), pointer :: maxLevelCell, indexToCellID

      integer theVertex, iLevel, iLevelBuoyancy, aVertex, &
        nSteps, timeStep, subStep, subStepOrder, timeInterpOrder, aTimeLevel, nCellVertices, blockProc, arrayIndex
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: nCellVerticesArray
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:,:), pointer :: cellsOnVertex
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
      integer, pointer :: verticalTreatment, indexLevel, iCell, currentCellGlobalID

      type (mpas_pool_type), pointer :: tracersPool
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, pointer :: index_temperature, index_salinity
      real (kind=RKIND), pointer :: particleTemperature
      real (kind=RKIND), pointer :: particleSalinity

      real (kind=RKIND), dimension(:,:,:), pointer :: ecosysTracers
      integer, pointer :: index_DIC, &
                          index_ALK, &
                          index_PO4, &
                          index_NO3, &
                          index_SiO3, &
                          index_NH4, &
                          index_Fe, &
                          index_O2
      real (kind=RKIND), pointer :: particleDIC
      real (kind=RKIND), pointer :: particleALK
      real (kind=RKIND), pointer :: particlePO4
      real (kind=RKIND), pointer :: particleNO3
      real (kind=RKIND), pointer :: particleSiO3
      real (kind=RKIND), pointer :: particleNH4
      real (kind=RKIND), pointer :: particleFe
      real (kind=RKIND), pointer :: particleO2

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block)) !{{{
        ! allocate scratch memory / setup pointers / get block
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackFields', lagrPartTrackFieldsPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackScratch', lagrPartTrackScratchPool)
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackCellsPool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        ! particlelist should be stored in the structs pool probably (need to seriously rework the code!)
        particlelist => block % particlelist

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'zCell', zCell)
        call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
        call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
        call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
        call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
        call mpas_pool_get_array(lagrPartTrackCellsPool, 'wachspressAreaB', areaBArray)
        buoyancyTimeInterp => potentialDensity

        call mpas_pool_get_field(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
        call mpas_dmpar_exch_halo_field(normalVelocity)
        call mpas_pool_get_field(statePool, 'layerThickness', layerThickness, timeLevel=timeLevel)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel=timeLevel)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        if (config_use_ecosysTracers) then
          call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel=timeLevel)
          call mpas_pool_get_dimension(tracersPool, 'index_DIC', index_DIC)
          call mpas_pool_get_dimension(tracersPool, 'index_ALK', index_ALK)
          call mpas_pool_get_dimension(tracersPool, 'index_PO4', index_PO4)
          call mpas_pool_get_dimension(tracersPool, 'index_NO3', index_NO3)
          call mpas_pool_get_dimension(tracersPool, 'index_SiO3', index_SiO3)
          call mpas_pool_get_dimension(tracersPool, 'index_NH4', index_NH4)
          call mpas_pool_get_dimension(tracersPool, 'index_Fe', index_Fe)
          call mpas_pool_get_dimension(tracersPool, 'index_O2', index_O2)
        end if

        call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
        call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
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
          call mpas_pool_get_array(particle % haloDataPool, 'lonParticle', lonParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'latParticle', latParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'currentCell', iCell)
          call mpas_pool_get_array(particle % haloDataPool, 'currentCellGlobalID', currentCellGlobalID)
          iCell = -1

          particlePosition(1) = xParticle
          particlePosition(2) = yParticle
          particlePosition(3) = zParticle

          ! compute particle lat/lon position
          call convert_latlon_from_xyz(latParticle, lonParticle, xParticle, yParticle, zParticle)

          call mpas_pool_get_array(particle % haloDataPool, 'zLevelParticle', zLevelParticle)

          call mpas_pool_get_array(particle % haloDataPool, 'verticalTreatment', verticalTreatment)
          call mpas_pool_get_array(particle % haloDataPool, 'indexLevel', indexLevel)
          call mpas_pool_get_array(particle % haloDataPool, 'dtParticle', dtParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'buoyancyParticle', buoyancyParticle)
          call mpas_pool_get_array(particle % haloDataPool, 'indexToParticleID', indexToParticleID)

          call mpas_pool_get_array(particle % haloDataPool, 'particleTemperature', particleTemperature)
          call mpas_pool_get_array(particle % haloDataPool, 'particleSalinity', particleSalinity)
          if (config_use_ecosysTracers) then
            call mpas_pool_get_array(particle % haloDataPool, 'particleDIC', particleDIC)
            call mpas_pool_get_array(particle % haloDataPool, 'particleALK', particleALK)
            call mpas_pool_get_array(particle % haloDataPool, 'particlePO4', particlePO4)
            call mpas_pool_get_array(particle % haloDataPool, 'particleNO3', particleNO3)
            call mpas_pool_get_array(particle % haloDataPool, 'particleSiO3', particleSiO3)
            call mpas_pool_get_array(particle % haloDataPool, 'particleNH4', particleNH4)
            call mpas_pool_get_array(particle % haloDataPool, 'particleFe', particleFe)
            call mpas_pool_get_array(particle % haloDataPool, 'particleO2', particleO2)
          end if


!          call mpas_pool_get_array(particle % haloDataPool, 'lonVel', lonVel)
!          call mpas_pool_get_array(particle % haloDataPool, 'latVel', latVel)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumU', sumU)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumV', sumV)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumUU', sumUU)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumUV', sumUV)
!          call mpas_pool_get_array(particle % haloDataPool, 'sumVV', sumVV)

          !}}}

          !!!!!!!!!! PERFORM SAMPLING (VELOCITY, TEMP, SALINITY, ETC) !!!!!!!!!!
          ! this could be done just before the output anyway

          ! need iCell computed for final position
          ! need scalar values interpolated in time to yield single value
          ! probably need to store zMid too, including flipping it
          ! get updated cell location
#ifdef MPAS_DEBUG
          call mpas_timer_start("get_validated_cell_idLPT")
#endif
          LIGHT_DEBUG_WRITE('sampling initialization')
          call get_validated_cell_id(nCells, xCell,yCell,zCell , xVertex,yVertex,zVertex, &
            particlePosition(1),particlePosition(2),particlePosition(3), meshPool, &
            nCellVerticesArray, verticesOnCell, iCell, nCellVertices, cellsOnCell)
#ifdef MPAS_DEBUG
          call mpas_timer_stop("get_validated_cell_idLPT")
#endif

          call get_vertical_level_info(iLevel, zLevelParticle, &
            verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
            zMid, indexLevel)

          ! store global cell ID
          currentCellGlobalID = indexToCellID(iCell)

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
            lagrPartTrackFieldsPool, &
            timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
            verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
            particlePosition, zLevelParticle, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex,  &
            meshPool, areaBArray)
          ! convert horizontal velocity to lat/lon velocity

          LIGHT_DEBUG_WRITE(iLevel COMMA particleVelocity)

!          ! store velocity for use in computing normalized autocorrelation offline
!          call mpas_convert_xyz_velocity_to_latlon(lonVel, latVel, particlePosition, particleVelocity)
          !}}}
          if ( associated(activeTracers) ) then
            if (config_AM_lagrPartTrack_sample_temperature) then
             !call mpas_log_write(' testing output for temperature')
             if (config_AM_lagrPartTrack_sample_horizontal_interp) then
               call interp_scalars(particleTemperature, activeTracers(index_temperature,:,:), &
                   verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                   particlePosition(1),particlePosition(2),particlePosition(3), &
                   xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                   yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                   zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                   xCell, yCell, zCell, &
                   meshPool, &
                   cellsOnVertex, &
                   verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                   zLevelParticle, zMid, indexLevel)
             else
               call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                 activeTracers(index_temperature,:,iCell), particleTemperature)
             end if
             !call mpas_log_write('done testing output for temperature')
            end if

            if (config_AM_lagrPartTrack_sample_salinity) then
             !call mpas_log_write(' testing output for salinity')
             if (config_AM_lagrPartTrack_sample_horizontal_interp) then
               call interp_scalars(particleSalinity, activeTracers(index_salinity,:,:), &
                   verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                   particlePosition(1),particlePosition(2),particlePosition(3), &
                   xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                   yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                   zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                   xCell, yCell, zCell, &
                   meshPool, &
                   cellsOnVertex, &
                   verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                   zLevelParticle, zMid, indexLevel)
             else
               call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                 activeTracers(index_salinity,:,iCell), particleSalinity)
             end if
             !call mpas_log_write('done testing output for salinity')
            end if
          end if

          if (config_use_ecosysTracers) then
            call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
            if ( associated(ecosysTracers) ) then
              if (config_AM_lagrPartTrack_sample_DIC) then
                !call mpas_log_write('testing output for DIC')
                if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                  call interp_scalars(particleDIC, ecosysTracers(index_DIC,:,:), &
                    verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                    particlePosition(1),particlePosition(2),particlePosition(3), &
                    xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                    yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                    zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                    xCell, yCell, zCell, &
                    meshPool, &
                    cellsOnVertex, &
                    verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                    zLevelParticle, zMid, indexLevel)
                else
                  call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                    ecosysTracers(index_DIC,:,iCell), particleDIC)
                end if
                !call mpas_log_write('done testing output for DIC')
              end if

              if (config_AM_lagrPartTrack_sample_ALK) then
                !call mpas_log_write('testing output for ALK')
                if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                  call interp_scalars(particleALK, ecosysTracers(index_ALK,:,:), &
                    verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                    particlePosition(1),particlePosition(2),particlePosition(3), &
                    xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                    yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                    zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                    xCell, yCell, zCell, &
                    meshPool, &
                    cellsOnVertex, &
                    verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                    zLevelParticle, zMid, indexLevel)
                else
                  call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                    ecosysTracers(index_ALK,:,iCell), particleALK)
                end if
                !call mpas_log_write('done testing output for ALK')
              end if

               if (config_AM_lagrPartTrack_sample_PO4) then
                 !call mpas_log_write('testing output for PO4')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particlePO4, ecosysTracers(index_PO4,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_PO4,:,iCell), particlePO4)
                 end if
                 !call mpas_log_write('done testing output for PO4')
               end if

               if (config_AM_lagrPartTrack_sample_NO3) then
                 !call mpas_log_write('testing output for NO3')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleNO3, ecosysTracers(index_NO3,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_NO3,:,iCell), particleNO3)
                 end if
                 !call mpas_log_write('done testing output for NO3')
               end if

               if (config_AM_lagrPartTrack_sample_SiO3) then
                 !call mpas_log_write('testing output for SiO3')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleSiO3, ecosysTracers(index_SiO3,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_SiO3,:,iCell), particleSiO3)
                 end if
                 !call mpas_log_write('done testing output for SiO3')
               end if

               if (config_AM_lagrPartTrack_sample_NH4) then
                 !call mpas_log_write('testing output for NH4')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleNH4, ecosysTracers(index_NH4,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_NH4,:,iCell), particleNH4)
                 end if
                 !call mpas_log_write('done testing output for NH4')
               end if

               if (config_AM_lagrPartTrack_sample_Fe) then
                 !call mpas_log_write('testing output for Fe')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleFe, ecosysTracers(index_Fe,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_Fe,:,iCell), particleFe)
                 end if
                 !call mpas_log_write('done testing output for Fe')
               end if

               if (config_AM_lagrPartTrack_sample_O2) then
                 !call mpas_log_write('testing output for O2')
                 if (config_AM_lagrPartTrack_sample_horizontal_interp) then
                   call interp_scalars(particleO2, ecosysTracers(index_O2,:,:), &
                     verticesOnCell(:,iCell), nCellVerticesArray(iCell), &
                     particlePosition(1),particlePosition(2),particlePosition(3), &
                     xVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     yVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     zVertex(verticesOnCell(1:nCellVerticesArray(iCell),iCell)), &
                     xCell, yCell, zCell, &
                     meshPool, &
                     cellsOnVertex, &
                     verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
                     zLevelParticle, zMid, indexLevel)
                 else
                   call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
                     ecosysTracers(index_O2,:,iCell), particleO2)
                 end if
                 !call mpas_log_write('done testing output for O2')
               end if
            end if
          end if

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
     integer :: aVertex, theVertex, iLevelTmp

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
         !  boundaryVertex, uvCell)

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
         !  boundaryVertex, uvCell)

         ! no vertical motion, just horizontal motion
         return
         !}}}

       case (3) verticalTreatmentCase !('passiveFloat') !{{{

         ! interpolate the vertical velocity (can be different than zMid...,
         ! which is assumed input)
         iLevelTmp = mpas_get_vertical_id(nVertLevels, -zLoc, -zTop)
         vertVelocityInterp = interp_vert_velocity_to_zlevel( &
           iLevelTmp, -zLoc, -zTop, vertVelocityTop)

         ! interpolate the horizontal velocity based on z-levels
         call interp_nodal_vectors(ncellvertices, verticesoncell, &
           ilevel, nVertLevels, zLoc, zmid, &
           uVertexVelocity, vVertexVelocity, wVertexVelocity, uvCell)

         ! ensure that there is zero nodal velocity on the boundary
         !call zero_boundary_nodal_values(nCellVertices, verticesOnCell, &
         !  boundaryVertex, uvCell)

         return
         !}}}

       case (4) verticalTreatmentCase !('buoyancySurface') !{{{
         ! no vertical velocity required because there is not vertical integration for position

         ! interpolate the horizontal velocity
         call interp_nodal_vectors(ncellvertices, verticesoncell, &
           ilevel, nVertLevels, phiInterp, phiMid, &
           uVertexVelocity, vVertexVelocity, wVertexVelocity, uvCell)

         ! ensure that there is zero nodal velocity on the boundary
         !call zero_boundary_nodal_values(nCellVertices, verticesOnCell, &
         !  boundaryVertex, uvCell)

         ! no vertical motion, just horizontal motion
         return
         !}}}

       case (5) verticalTreatmentCase !('argoFloat')  !{{{


         !}}}

       case default verticalTreatmentCase !{{{
         LIGHT_ERROR_WRITE('Vertical treatment for particle integration unknown (' COMMA verticalTreatment COMMA ')!')
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
      lagrPartTrackFieldsPool, &
      timeInterpOrder, timeCoeff, iCell, iLevel, buoyancyInterp, maxLevelCell, &
      verticalTreatment, indexLevel, nCellVertices, verticesOnCell, boundaryVertex, &
      xSubStep, zSubStep, zMid, zTop, vertVelocityTop, xVertex, yVertex, zVertex, meshPool, areaBArray) !{{{
    implicit none

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (mpas_pool_type), pointer, intent(in) :: lagrPartTrackFieldsPool
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
    call mpas_timer_start("velocity_time_interpolationLPT")
#endif

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
    call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
    call mpas_pool_get_config(meshPool, 'x_period', x_period)
    call mpas_pool_get_config(meshPool, 'y_period', y_period)

    ! allocations for particular cell
    allocate(vertCoords(3,nCellVertices), uvCell(3,nCellVertices), areaB(nCellVertices))

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
      buoyancy => potentialDensity
      !}}}

      ! get final, interpolated particle velocity at this point (collapse to point)
      if (verticalTreatment == 4) then
        ! buoyancy case (not using zSubStep for interpolation / iLevel)
        ! use existing code noting we need to flip the order to get the right iLevel

        iLevel = mpas_get_vertical_id(maxLevelCell(iCell), buoyancyInterp, buoyancy(:,iCell))
        LIGHT_DEBUG_WRITE('iLevel=' COMMA iLevel)
        ! note, if buoyancyInterp out of range this will try to reorient the particle to the top / bottom but there
        ! will definitely be some error with this type of computation because the buoyancy is not available at this location
        ! the time interpolation, as a consequency, can mix velocities from different buoyancy surfaces in order to advect
        ! the particle
        !call mpas_log_write( 'iLevel = $i',MPAS_LOG_OUT, intArgs=(/ iLevel /) )
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
      call mpas_timer_start("part_horiz_interpLPT")
#endif
      LIGHT_DEBUG_WRITE('particleVelocityVert=' COMMA particleVelocityVert)
      particleVelocity = particleVelocity + &
        timeCoeff(aTimeLevel) * particle_horizontal_interpolation(nCellVertices, vertCoords, &
        xSubStep, uvCell, meshPool, areaB)
#ifdef MPAS_DEBUG
      call mpas_timer_stop("part_horiz_interpLPT")
#endif
      LIGHT_DEBUG_WRITE('particleVelocity=' COMMA particleVelocity)


    end do

    ! deallocations of temp memory
    deallocate(vertCoords, uvCell, areaB)

#ifdef MPAS_DEBUG
    call mpas_timer_stop("velocity_time_interpolationLPT")
#endif

  end subroutine velocity_time_interpolation !}}}

!  subroutine zero_autocorrelation_sums(domain) !{{{
!    implicit none
!
!    ! input/output variables
!    type (domain_type), intent(inout) :: domain
!    ! local
!    type (block_type), pointer :: block
!    type (mpas_particle_list_type), pointer :: particlelist
!    type (mpas_particle_type), pointer :: particle
!    ! output variables (per particle)
!    real (kind=RKIND), pointer :: sumU, sumV, sumUU, sumUV, sumVV
!    integer, pointer :: currentCell
!
!    ! get the appropriate pools
!      block => domain % blocklist
!      do while (associated(block)) !{{{
!        particlelist => block % particlelist
!        do while(associated(particlelist)) !{{{
!          ! get pointers / option values
!          particle => particlelist % particle
!
!          ! get values (may want a flag for reinitialization in the future)
!          !call mpas_pool_get_array(particle % haloDataPool, 'sumU', sumU)
!          !call mpas_pool_get_array(particle % haloDataPool, 'sumV', sumV)
!          !call mpas_pool_get_array(particle % haloDataPool, 'sumUU', sumUU)
!          !call mpas_pool_get_array(particle % haloDataPool, 'sumUV', sumUV)
!          !call mpas_pool_get_array(particle % haloDataPool, 'sumVV', sumVV)
!          call mpas_pool_get_array(particle % haloDataPool, 'currentCell', currentCell)
!
!          ! initialize the values
!          !sumU = 0.0_RKIND
!          !sumV = 0.0_RKIND
!          !sumUU = 0.0_RKIND
!          !sumUV = 0.0_RKIND
!          !sumVV = 0.0_RKIND
!          currentCell = -1
!
!          ! get next particle to process on the list
!          particlelist => particlelist % next
!        end do !}}}
!
!        ! get next block
!        block => block % next
!      end do !}}}
!
!  end subroutine zero_autocorrelation_sums !}}}

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
     LIGHT_DEBUG_ALL_WRITE('lambda=' COMMA lambda)
     LIGHT_DEBUG_ALL_WRITE('uVertex=' COMMA uVertex(1,:))
     LIGHT_DEBUG_ALL_WRITE('vVertex=' COMMA uVertex(2,:))
     LIGHT_DEBUG_ALL_WRITE('wVertex=' COMMA uVertex(3,:))

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
     character (len=StrKIND) :: message
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
       LIGHT_DEBUG_ALL_WRITE('dpParticle = ' COMMA dpParticle)
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
       LIGHT_DEBUG_ALL_WRITE('dpParticle    = ' COMMA dpParticle)
       LIGHT_DEBUG_ALL_WRITE('pParticle     = ' COMMA pParticle)
       LIGHT_DEBUG_ALL_WRITE('pParticleTemp = ' COMMA pParticleTemp)
       write(message, *) 'mpas_arc_length = ',  mpas_arc_length( pParticle(1), pParticle(2), pParticle(3), &
                         pParticleTemp(1), pParticleTemp(2), pParticleTemp(3)), 'lenPath = ', lenPath
       LIGHT_DEBUG_ALL_WRITE(message)
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
!  routine get_vertical_level_info
!
!> \brief   Returns iLevel and zLevelParticle information for vertical location
!> \author  Phillip Wolfram
!> \date    11/13/2018
!> \details
!>  This routine computes the iLevel and zLevelParticle associated with a particular
!>  verticalTreatment.
!
!-----------------------------------------------------------------------
   subroutine get_vertical_level_info(iLevel, zLevelParticle, &
       verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
       zMid, indexLevel, ssh, bottomDepth, indexToParticleID)  !{{{

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     integer, intent (in) :: verticalTreatment                        !< verticalTreatment used for particles
     integer, dimension(:), pointer :: maxLevelCell                   !< maximum level for each horizontal cell
     integer, intent(in) :: iCell                                     !< index for horizontal cell location
     real (kind=RKIND), intent(in) :: buoyancyParticle                !< buoyancy value for particle
     real (kind=RKIND), dimension(:,:), pointer :: buoyancyTimeInterp !< intepolated buoyancy value
     real (kind=RKIND), dimension(:,:), pointer :: zMid               !< mid depth z value for each cell
     integer, intent(in) :: indexLevel                                !< index level for particle vertical location
     real (kind=RKIND), dimension(:), pointer, optional :: bottomDepth, ssh  !< error checking variables
     integer, intent(in), pointer, optional :: indexToParticleID      !< index to particle ID


     !-----------------------------------------------------------------
     ! input/output variables
     !-----------------------------------------------------------------

     integer, intent(inout) :: iLevel                                 !< vertical level / cell of phiInterp
     real (kind=RKIND), intent(inout) :: zLevelParticle                     !< vertical level for particle

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     REAL (kind=RKIND), parameter :: eps=1e-12


     ! compute zLevelParticle and iLevel ! {{{
     if(verticalTreatment == 4) then  !('buoyancySurface')
       iLevel = mpas_get_vertical_id(maxLevelCell(iCell), buoyancyParticle, buoyancyTimeInterp(:,iCell))

       ! get zLevelParticle by interpolating the scalars now (assumes that
       ! scalar value is constant within a particular cell)
       call interp_cell_scalars(iLevel, maxLevelCell(iCell), buoyancyParticle, buoyancyTimeInterp(:,iCell), &
         zMid(:,iCell), zLevelParticle)
     else if (verticalTreatment == 1) then !('indexLevel')
       zLevelParticle = zMid(indexLevel, iCell)
       iLevel = indexLevel
     else if (verticalTreatment == 3 .or. verticalTreatment == 2) then !('passiveFloat' or 'fixedZLevel')
       ! make sure final zLevelParticle is ok so that it can't extent past zMid range
       ! note, zlevel can be modified near top / bottom of ocean to keep particle in bounds
       zLevelParticle = min(max(zMid(maxLevelCell(iCell),iCell), zLevelParticle), zMid(1,iCell))
       iLevel = mpas_get_vertical_id(maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell))
       !call mpas_log_write('iLevel = $i', intArgs=(/iLevel/))
     else
       call mpas_log_write('Vertical treatment is not defined!', MPAS_LOG_CRIT)
     end if
     ! }}}

     !verify particles are still in domain !{{{
     if ( present(ssh) .and. present(bottomDepth) .and. present(indexToParticleID) ) then
#ifdef MPAS_DEBUG
       call mpas_timer_start("verifyParticleDepthInDomainLPT")
#endif

      ! debug !{{{
!#ifdef MPAS_DEBUG
!       ! check and verify that free floating particle hasn't gone out of boundary
!       if (iLevel == -1) then
!         call mpas_log_write(' Warning: Particle $i location of $r is higher than zMid of $r', &
!           realArgs=(/ zLevelParticle, zMid(1, iCell) /), intArgs=(/iCell/))
!       end if
!       if (iLevel == 0) then
!         call mpas_log_write(' Warning: Particle $i location of $r is deeper than zMid at bottom $r', &
!           realArgs=(/ zLevelParticle, zMid(maxLevelCell(iCell),iCell) /), intArgs=(/iCell/))
!       end if
!#endif
      !}}}

       ! verify that it is less than SSH and greater than bottomDepth
       if (zLevelParticle + bottomDepth(iCell) < eps) then
         call mpas_log_write(' Error: $i ($i) Particle of type $i height $r is lower than bottom depth $r', &
           realArgs=(/ zLevelParticle, -bottomDepth(iCell) /), intArgs=(/iCell, indexToParticleID, verticalTreatment/),&
           messageType=MPAS_LOG_ERR)
       end if
       if (zLevelParticle - ssh(iCell) > eps) then
         call mpas_log_write(' Error: Particle $i height $r is higher than ssh of $r', &
           realArgs=(/ zLevelParticle, ssh(iCell) /), intArgs=(/iCell/), &
           messageType=MPAS_LOG_ERR)
       end if
#ifdef MPAS_DEBUG
       call mpas_timer_stop("verifyParticleDepthInDomainLPT")
#endif
     end if
     !}}}

   end subroutine get_vertical_level_info !}}}


!***********************************************************************
!
!  routine interp_scalars
!
!> \brief   Interpolate scalars horizontal with linear interpolation and
!>          based on a z level criteria (z-level, buoyancy, etc)
!> \author  Phillip Wolfram
!> \date    11/13/2018
!> \details
!>  This routine interpolates cell scalars linearly horizontally and veritcally
!>  depending upon a criteria such as z-level, buoyancy, etc.
!
!-----------------------------------------------------------------------
   subroutine interp_scalars(particleScalar, tracerarray, & !{{{
       verticesOnTheCell, nVertices, &
       xp,yp,zp, &
       xv,yv,zv, &
       xc,yc,zc, &
       meshPool, &
       cellsOnVertex, &
       verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
       zLevelParticleIn, zMid, indexLevel)

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     integer, intent (in) :: verticalTreatment                        !< verticalTreatment used for particles
     integer, dimension(:), pointer :: maxLevelCell                   !< maximum level for each horizontal cell
     real (kind=RKIND), intent(in) :: buoyancyParticle                !< buoyancy value for particle
     real (kind=RKIND), dimension(:,:), pointer :: buoyancyTimeInterp !< intepolated buoyancy value
     real (kind=RKIND), dimension(:,:), pointer :: zMid               !< mid depth z value for each cell
     integer, intent(in) :: indexLevel                                !< index level for particle vertical location
     integer, dimension(:), intent(in) :: verticesOnTheCell           !< vertex indices on cell
     integer, dimension(:,:), pointer, intent(in) :: cellsOnVertex    !< cells ajoined to vertex
     integer, intent(in) :: nVertices                                 !< number of vertices on the cell
     real (kind=RKIND), intent(in) :: xp,yp,zp                        !< particle location
     real (kind=RKIND), dimension(:), intent(in) :: xv,yv,zv          !< spatial location of cell vertices
     real (kind=RKIND), dimension(:), intent(in) :: xc,yc,zc          !< spatial location of cell centers
     type (mpas_pool_type), intent(in), pointer :: meshPool           !< meshPool pointer
     real (kind=RKIND), dimension(:,:), intent(in) :: tracerarray     !< tracer array to be interpolated
     real (kind=RKIND), pointer :: zLevelParticleIn                   !< vertical level z location

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), pointer :: particleScalar                     !< interpolated scalar value for the particle

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND) :: zLevelParticle                     !< vertical level z location
     integer :: iLevel                                       !< vertical level index
     integer :: iCell                                        !< local cell index for horizontal cell location
     integer :: triangleIndex                                !< index of triangle for interpolation
     real (kind=RKIND), dimension(3) :: scalarTriangle       !< interpolated values on triangle
     integer :: iVertex                                      !< vertex index variable
     real (kind=RKIND), dimension(3, nVertices) :: xVertex   !< coordinates of vertex
     logical, dimension(nVertices) :: interiorVertex         !< determines if vertex is on the interior
     real (kind=RKIND), dimension(3) :: xPoint               !< coordinates of transformed point for nearest neighbor calc
     real (kind=RKIND), dimension(3) :: lambda               !< barycentric weights for triangle

     real (kind=RKIND) :: pointRadius   ! magnitude of a point radius
     logical, pointer :: on_a_sphere    ! flag designating if we are on a sphere
     logical, pointer :: is_periodic    ! flag designating if periodicity is important
     real(kind=RKIND), pointer :: x_period, y_period
     integer, parameter :: vertexDegree = 3                  !< max number of cells on Vertex
     real(kind=RKIND) :: dx, dy, dz, r2, r2Min
     real (kind=RKIND), dimension(3,3) :: vertexTriangle
     integer :: aVertex, vertexID, aCell, aPoint, localCell
     integer, pointer :: nCells
     real(kind=RKIND), parameter :: eps = 1.0e-12_RKIND
     integer :: ierr = 0
     logical :: boundaryInterp

     ! assume that we cannot do spatial interpolation
     vertexID = -1

     !-----------------------------------------------------------------
     ! find the vertexID of the polygon containing the point aPoint
     !-----------------------------------------------------------------

     ! {{{
     ! for each vertex assess whether it is an interior vertex
     ! interior vertex until proven otherwise
     interiorVertex(:) = .true.
     do aPoint = 1, nVertices
       do aVertex = 1, vertexDegree
         !call mpas_log_write('CellsOnVertex($i, $i) = $i', intArgs = (/aPoint, aVertex, cellsOnVertex(aPoint, aVertex)/))
         if (cellsOnVertex(aVertex, verticesOnTheCell(aPoint)) == 0) then
           interiorVertex(aPoint) = .false.
         end if
       end do
     end do
     !! hard code as exterior vertex for error checking / debugging
     !interiorVertex(:) = .false.

     ! prepare vertex array xVertex needed to find closest vertex to point xPoint ! {{{
     call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
     call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
     call mpas_pool_get_config(meshPool, 'x_period', x_period)
     call mpas_pool_get_config(meshPool, 'y_period', y_period)
     call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

     ! check to see if point is in vertex polygon
     xPoint(1) = xp
     xPoint(2) = yp
     xPoint(3) = zp

     ! for each vertex assess whether it is valid to use
     do aVertex = 1, nVertices
       if (interiorVertex(aVertex)) then
         do aCell = 1, vertexDegree

           localCell = cellsOnVertex(aCell, verticesOnTheCell(aVertex))

           ! get coordinates of vertex for triangle for barycentric interpolation prepared
           if ( (.not. on_a_sphere) .and. is_periodic) then
             vertexTriangle(aCell, 1) = mpas_fix_periodicity(xc(localCell), xp, x_period)
             vertexTriangle(aCell, 2) = mpas_fix_periodicity(yc(localCell), yp, y_period)
             vertexTriangle(aCell, 3) = zc(localCell)
           else
             vertexTriangle(aCell, 1) = xc(localCell)
             vertexTriangle(aCell, 2) = yc(localCell)
             vertexTriangle(aCell, 3) = zc(localCell)
           end if
         end do

         ! check to make sure the point xPoint is in the vertexTriangle
         if ( mpas_point_in_polygon(xPoint, vertexTriangle, on_a_sphere) ) then
           vertexID = aVertex
         end if
         ! otherwise default to "boundary" or edge case
       end if
     end do

     !! normalize 3D points on sphere
     !if ( on_a_sphere .and. .not. is_periodic) then
     !  pointRadius = sqrt(xp*xp + yp*yp + zp*zp)
     !  xPoint = xPoint/pointRadius
     !  vertexTriangle = vertexTriangle/pointRadius
     !end if
     !}}}

     ! debugging output !{{{
!#ifdef MPAS_DEBUG
!     ! check to make sure nearest vertex found makes sense
!     call mpas_log_write('vertexID = $i, xPoint = $r $r $r', intArgs = (/ vertexID /), &
!       realArgs = (/ xPoint(1), xPoint(2), xPoint(3) /))
!     do aVertex = 1, nVertices
!       call mpas_log_write('$i xVertex = $r $r $r$', intArgs = (/ aVertex /), &
!         realArgs = (/ &
!         xVertex(1, aVertex), &
!         xVertex(2, aVertex), &
!         xVertex(3, aVertex)  &
!         /))
!     end do
!#endif
     !}}}
     !}}}
     ! }}}

     !-----------------------------------------------------------------
     ! for each vertex, determine how to interpolate from neighbooring
     ! cell values if an edge-case; if a boundary just do vertical
     ! interpolation
     !-----------------------------------------------------------------

     ! {{{
     ! cellsOnVertex assumes an "empty" neighbor (edge case) is 0
     do aCell = 1, vertexDegree

       boundaryInterp = .false.

       ! only handles the case where the vertex is an exterior vertex, otherwise treat as cell-based value
       ! (no horizontal interpolation)
       if (vertexID .eq. -1) then
         boundaryInterp = .true.
       else if (cellsOnVertex(aCell, verticesOnTheCell(vertexID)) == 0 .or. &
                cellsOnVertex(aCell, verticesOnTheCell(vertexID)) == nCells+1 .or. &
                (cellsOnVertex(1, verticesOnTheCell(vertexID)) == cellsOnVertex(2, verticesOnTheCell(vertexID))) .or. &
                (cellsOnVertex(2, verticesOnTheCell(vertexID)) == cellsOnVertex(3, verticesOnTheCell(vertexID))) .or. &
                (cellsOnVertex(3, verticesOnTheCell(vertexID)) == cellsOnVertex(1, verticesOnTheCell(vertexID)))) then
         boundaryInterp = .true.
       end if

       if (boundaryInterp) then

         !call mpas_log_write('Using column interpolation for $i because it is on an edge', intArgs = (/iCell/))

         ! note that zLevelParticle corresponding to adjacent cells will be a
         ! little different on boundary, hence the need to use this temporary variable
         ! so interpolation near  boundaries interpolates on the surface
         ! of the boundary (e.g., to make sure interpolated values are within interior
         ! away from dz/2 of column boundary)
         zLevelParticle = zLevelParticleIn

         ! in the case of any edge, just use cell value and return it
         call get_vertical_level_info(iLevel, zLevelParticle, &
           verticalTreatment, maxLevelCell, iCell, buoyancyParticle, buoyancyTimeInterp, &
           zMid, indexLevel)

         call interp_cell_scalars(iLevel, maxLevelCell(iCell), -zLevelParticle, -zMid(:,iCell), &
           tracerarray(:,iCell), particleScalar)

         return
       end if
     end do
     ! }}}

     !-----------------------------------------------------------------
     ! perform barycentric interpolation horizontal from
     ! vertical interpolants
     !-----------------------------------------------------------------

     ! {{{
     ! all cellsOnVertex(verticesOnTheCell(vertexID), 1:3) are real cells, for each of triangle vertices
     ! get cell values for the vertical interpolant for each point of triangle for interpolation
     vertexTriangle = 0.0_RKIND
     do aCell = 1, vertexDegree

       localCell = cellsOnVertex(aCell, verticesOnTheCell(vertexID))

       ! perform vertical interpolation for each cell
       ! return iLevel, zLevelParticle for vertical interpolation of the scalar
       ! note that zLevelParticle corresponding to adjacent cells will be a
       ! little different on boundary, hence the need to use this temporary variable
       ! so interpolation near  boundaries interpolates on the surface
       ! of the boundary (e.g., to make sure interpolated values are within interior
       ! away from dz/2 of column boundary)
       zLevelParticle = zLevelParticleIn
       call get_vertical_level_info(iLevel, zLevelParticle, &
         verticalTreatment, maxLevelCell, localCell, &
         buoyancyParticle, buoyancyTimeInterp, &
         zMid, indexLevel)

       ! perform vertical interpolation for the scalar
       call interp_cell_scalars(iLevel, maxLevelCell(localCell), &
         -zLevelParticle, -zMid(:,localCell), &
         tracerarray(:,localCell), scalarTriangle(aCell))

       ! get coordinates of vertex for triangle for barycentric interpolation prepared
       if ( (.not. on_a_sphere) .and. is_periodic) then
         vertexTriangle(1, aCell) = mpas_fix_periodicity(xc(localCell), xp, x_period)
         vertexTriangle(2, aCell) = mpas_fix_periodicity(yc(localCell), yp, y_period)
         vertexTriangle(3, aCell) = zc(localCell)
       else
         vertexTriangle(1, aCell) = xc(localCell)
         vertexTriangle(2, aCell) = yc(localCell)
         vertexTriangle(3, aCell) = zc(localCell)
       end if
     end do

     ! perform barycentric interpolation in the horizontal
     xPoint(1) = xp
     xPoint(2) = yp
     xPoint(3) = zp

     !! normalize 3D points on sphere
     !if ( on_a_sphere .and. .not. is_periodic) then
     !  pointRadius = sqrt(xp*xp + yp*yp + zp*zp)
     !  xPoint = xPoint/pointRadius
     !  vertexTriangle = vertexTriangle/pointRadius
     !end if

     !!! this wachspress routine is designed for
     !!! polygons with more than 3 sides (making it have higher cost too) ...
     lambda = mpas_wachspress_coordinates(vertexDegree, vertexTriangle, xPoint, meshPool)

     ! more robust (and less expensive) barycentric interpolation approach,
     ! assuming a PURE interpolation (which cannot be guaranteed here)
     !call mpas_calculate_barycentric_weights(xPoint, &
     !  vertexTriangle(:,1), vertexTriangle(:,2), vertexTriangle(:,3), &
     !  meshPool, &
     !  lambda(1), lambda(2), lambda(3), ierr)

! debugging output !{{{
#ifdef MPAS_DEBUG
     ! ensure interpolation is obtained
     if (minval(lambda) < -eps .or. maxval(lambda) > (1.0_RKIND + eps) .or. &
         abs(1.0_RKIND - sum(lambda)) > eps ) then
       call mpas_log_write('Not a valid barycentric interpolation: ' // &
         'lambda = $r $r $r vertexID = $i cellsOnVertex(:, verticesOnTheCell(vertexID)) = $i $i $i', MPAS_LOG_ERR, &
         realArgs=(/lambda(1), lambda(2), lambda(3)/), &
         intArgs=(/vertexID, &
         cellsOnVertex(1, verticesOnTheCell(vertexID)), &
         cellsOnVertex(2, verticesOnTheCell(vertexID)), &
         cellsOnVertex(3, verticesOnTheCell(vertexID))/))
     end if
#endif
!}}}

     particleScalar = scalarTriangle(1) * lambda(1) + &
                      scalarTriangle(2) * lambda(2) + &
                      scalarTriangle(3) * lambda(3)
     !}}}

     ! debugging output !{{{
!#ifdef MPAS_DEBUG
!     ! check to make sure nearest vertex found makes sense
!     call mpas_log_write('vertexID = $i, xPoint = $r $r $r', intArgs = (/ vertexID /), &
!       realArgs = (/ xPoint(1), xPoint(2), xPoint(3) /))
!     do aVertex = 1, nVertices
!       call mpas_log_write('$i xVertex = $r $r $r$', intArgs = (/ aVertex /), &
!         realArgs = (/ &
!         xVertex(1, aVertex), &
!         xVertex(2, aVertex), &
!         xVertex(3, aVertex) /))
!     end do
!     do aCell = 1, vertexDegree
!       localCell = cellsOnVertex(aCell, verticesOnTheCell(vertexID))
!       call mpas_log_write(' cell id $i $r $r $r', intArgs=(/ localCell/), &
!         realArgs = (/&
!         vertexTriangle(1, aCell), &
!         vertexTriangle(2, aCell), &
!         vertexTriangle(3, aCell)/))
!     end do
!     call mpas_log_write('scalar Triangle = $r $r $r', realArgs = &
!       (/scalarTriangle(1), &
!         scalarTriangle(2), &
!         scalarTriangle(3)/))
!     call mpas_log_write('lambda = $r $r $r', realArgs = &
!       (/lambda(1), &
!         lambda(2), &
!         lambda(3)/))
!     call mpas_log_write('ierr = $i', intArgs=(/ierr/))
!#endif
     !}}}

   end subroutine interp_scalars !}}}


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
     real (kind=RKIND), intent(in) :: zInterp                         !< location to interpolate
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex
                                                                      !< velocities are defined)

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
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex
                                                                      !< velocities are defined)

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
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex
                                                                      !< velocities are defined)

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

     !call mpas_log_write( 'interp' )
     if(iLevel < 1) then
       if(iLevel == 0) then
         !print *, 'nCellVertices = ', nCellVertices, 'buoyancyInterp= ', phiInterp
         do aVertex = 1, nCellVertices
           theVertex = verticesOnCell(aVertex)
           !print *, maxloc(phiVals,1), phiVals(1:nVertLevels), uVertexVelocity(1:nVertLevels,theVertex)
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
     real (kind=RKIND), dimension(:,:), intent(inout) :: uvCell  !< components of vertex velocity (vertically selected)

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     integer :: aVertex, theVertex

     ! make sure to mask all boundary vertexes to be zero to enforce boundary conditions
     do aVertex = 1, nCellVertices
       theVertex = verticesOnCell(aVertex)
       !print *, boundaryVertex(theVertex)
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
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where vertex
                                                                      !< velocities are defined)

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     integer, intent(out) :: iLow, iHigh                              !< interpolation indices

     ! assumes increasing index is in decreasing phi space, and phiInterp in range of phiVals
     !print *, 'phiInterp = ', phiInterp, 'phiVals = ', phiVals, 'nVertLevels = ', nVertLevels, 'iLevel = ', iLevel
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
       !call mpas_log_write( 'Warning: fast interpolation failed, trying general, brute force search for interpolation bounds')
       !print *,  'iLow = ', iLow, ' iHigh = ', iHigh
       !print *,  'phiInterp = ', phiInterp , ' phiLow =',  phiVals(iLow), ' phiHigh = ', phiVals(iHigh)
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
     real (kind=RKIND), dimension(:), intent(in) :: phiVals           !< values at elevation of cell middle (where
                                                                      !< vertex velocities are defined)

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     integer, intent(out) :: iLow, iHigh                          !< indexes for the high and low components for the interpolant

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     integer :: aLevel
     character (len=StrKIND) :: message

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
       else
         call mpas_log_write('Value for interpolation was not found!', MPAS_LOG_CRIT)
       end if
     end do

#ifdef MPAS_DEBUG
     if(phiVals(iLow) <= phiInterp .and. phiInterp <= phiVals(iHigh)) then
       ! we are ok
       LIGHT_DEBUG_ALL_WRITE('brute force interpolation successful')
       LIGHT_DEBUG_ALL_WRITE('iLow = ' COMMA iLow COMMA ' iHigh = ' COMMA iHigh)
       LIGHT_DEBUG_ALL_WRITE('phiInterp = ' COMMA phiInterp COMMA ' phiLow =' COMMA phiVals(iLow) COMMA ' phiHigh = ' COMMA phiVals(iHigh))
       LIGHT_DEBUG_ALL_WRITE(message)
     else
       LIGHT_DEBUG_ALL_WRITE('brute force interpolation failed with phiInterp = ' COMMA phiInterp COMMA ' phiLow = ' COMMA phiVals(iLow) COMMA ' phiHigh = ' COMMA phiVals(iHigh))
       LIGHT_DEBUG_ALL_WRITE(' phiVals = ' COMMA phiVals(1:nVertLevels))
     end if

     LIGHT_DEBUG_WRITE('Warning!: brute force interpolation used, boundary condition may be wrong!')
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
         call mpas_log_write('Error with vertical velocity interpolation! ' // &
           'iLevel = $i, zTop(iLevel) = $r, zTop(iLevel+1)=$r zTop(iLevel+2) = $r zSubstep=$r', &
           realArgs =(/zTop(iLevel), zTop(iLevel+1), zTop(iLevel+2), zSubStep/), intArgs=(/iLevel/), messageType=MPAS_LOG_CRIT)
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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  convert_latlon_from_xyz
!
!> \brief
!> \author Adrian K. Turner and Phillip J. Wolfram, LANL
!> \date 2013-2014, 10/08/2018
!> \details
!>  Adapted from seaice_latlon_from_xyz, with
!>          1. assume xyz points define r
!>          2. shift to be in range of 0 to 2*pi as in mesh tools
!>
!> Note, should be moved into framework.
!
!-----------------------------------------------------------------------
  subroutine convert_latlon_from_xyz(lat, lon, x, y, z) !{{{
    implicit none

    ! given xyz coordinates determine the latitude and longitude
    real(kind=RKIND), intent(out) :: &
         lat, & !< Output: latitude of point
         lon    !< Output: longitude of point

    real(kind=RKIND), intent(in) :: &
         x,   & !< Input: x position of point
         y,   & !< Input: y position of point
         z      !< Input: z position of point

    real(kind=RKIND) :: r      !< Computed: radius of earth

    r = sqrt(x*x + y*y + z*z)
    lon = atan2(y, x)
    lat = asin(z/r)

    ! ensure range in 0, 2*pi
    if (lon < 0.0_RKIND) then
      lon = 2*pi + lon
    end if

  end subroutine convert_latlon_from_xyz !}}}

!}}}

end module ocn_lagrangian_particle_tracking

! vim: foldmethod=marker
