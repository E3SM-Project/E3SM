! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for the split explicit
!>  time integration scheme
!
!-----------------------------------------------------------------------

module ocn_time_integration_split

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_threading
   use mpas_timekeeping
   use mpas_log

   use ocn_config
   use ocn_mesh
   use ocn_tendency
   use ocn_diagnostics_variables
   use ocn_diagnostics
   use ocn_gm
   use ocn_submesoscale_eddies

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_vertical_advection
   use ocn_vertical_remap
   use ocn_time_average_coupled

   use ocn_effective_density_in_land_ice
   use ocn_surface_land_ice_fluxes
   use ocn_transport_tests

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

   public :: ocn_time_integrator_split, &
             ocn_time_integration_split_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   character (len=*), parameter :: &
      subcycleGroupName = 'subcycleFields', &! name for cached halo exch
      finalBtrGroupName = 'finalBtrFields'   ! name for cached halo exch

   logical :: &
      unsplit           ! flag for unsplit_explicit time option

   integer, dimension(:), allocatable :: &
      numClinicIterations  ! number of baroclinic iterations for each
                           ! outer timestep iteration

   integer :: &
      neededHalos,     &! number of halo levels needed
      haloDecrement,   &! number of halos to invalidate each cycle
      numTSIterations, &! number of outer timestep iterations
      nBtrSubcycles     ! number of barotropic subcycles

   integer :: &
      thickEdgeFluxChoice ! choice of thickness flux type
   integer, parameter :: &
      thickEdgeFluxCenter   = 1, &! use mean thickness of cell neighbors
      thickEdgeFluxUpwind   = 2   ! use upwind cell thickness at edge
   real (kind=RKIND) :: pgf_sal_on

   real (kind=RKIND) :: &
      useVelocityCorrection,&! mask for velocity correction
      splitFact, &! mask for terms in split or unsplit cases
      self_attraction_and_loading_beta

!***********************************************************************

   contains

!***********************************************************************
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates the ocean model forward one master time
!>  step (dt) using a split-explicit method in which the fast
!>  barotropic mode is sub-cycled using a shorter explicit time step.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_split(domain, dt)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt              !< [in] time step (sec) to move forward

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: &
         block ! structure with subdomain data

      type (mpas_pool_type), pointer :: &
         statePool,         &! structure holding state variables
         tracersPool,       &! structure holding tracers
         meshPool,          &! structure holding mesh variables
         verticalMeshPool,  &! structure holding vertical mesh variables
         tendPool,          &! structure holding tendencies
         tracersTendPool,   &! structure holding tracer tendencies
         forcingPool,       &! structure holding forcing variables
         scratchPool,       &! structure holding temporary variables
         swForcingPool,     &! structure holding short-wave forcing vars
         tracersIdealAgeFieldsPool     ! structure holding idealAge-related fields

      logical :: &
         activeTracersOnly   ! only compute tendencies for active tracers

      logical, pointer :: &
         config_use_tracerGroup  ! flag for using each tracer group

      integer :: &
         iCell, iEdge, k, &! loop iterators for cell, edge, vert loops
         nCells, nEdges,  &! number of cells or edges (incl halos)
         i,j,             &! generic loop iterators
         cell1, cell2,    &! neighbor cell addresses
         eoe,             &! index for edge on edge
         err,             &! local error flag
         splitExplicitStep, &! loop index for outer timestep loop
         oldBtrSubcycleTime,&! time index for old barotropic value
         newBtrSubcycleTime,&! time index for new barotropic value
         uPerpTime,       &! time index to use for uPerp calculation
         BtrCorIter,      &! loop index for barotropic coriolis cycle
         stage1_tend_time,&! time index for stage 1 tendencies
         edgeHaloComputeCounter, &! halo counters to reduce
         cellHaloComputeCounter   ! halo updates during cycling

      real (kind=RKIND) :: &
         normalThicknessFluxSum, &! sum of thickness flux in column
         thicknessSum,     &! sum of thicknesses in column
         flux,             &! temp for computing flux for barotropic
         sshEdge,          &! sea surface height at edge
         CoriolisTerm,     &! temp for computing coriolis term (fuperp)
         normalVelocityCorrection, &! velocity correction
         temp,             &! temp for holding vars at new time
         temp_h,           &! temporary for phi
         temp_mask,        &! temporary mask (edgeMask)
         lat,              &! cell latitude
         sshCell1,         &! sea sfc height in neighboring cells
         sshCell2,         &
         sshGrad            ! temp for holding ssh gradient

      real (kind=RKIND), dimension(:,:), allocatable:: &
         uTemp

      real (kind=RKIND), dimension(:), allocatable :: &
         btrvel_temp,      &
         sshTemp,          &
         deltaSSH,         &
         wctEdge,          &! flux water column thickness at edge
         baroclinicThickness

      ! State Array Pointers
      real (kind=RKIND), dimension(:), pointer :: &
         sshSubcycleCur,              &! subcycl ssh    at current time
         sshSubcycleNew,              &! subcycl ssh    at new     time
         sshSubcycleCurWithTides,     &! subcycl ssh    at current time
         sshSubcycleNewWithTides,     &! subcycl ssh    at new     time
         normalBarotropicVelocitySubcycleCur,&! barotropic vel subcyc
         normalBarotropicVelocitySubcycleNew,&! at current, new times
         sshCur,                      &! sea sfc height at current time
         sshNew,                      &! sea sfc height at new     time
         normalBarotropicVelocityCur, &! barotropic vel at current time
         normalBarotropicVelocityNew   ! barotropic vel at new     time

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalBaroclinicVelocityCur, &! baroclinic vel at current time
         normalBaroclinicVelocityNew, &! baroclinic vel at new     time
         normalVelocityCur,           &! full velocity  at current time
         normalVelocityNew,           &! full velocity  at new     time
         layerThicknessCur,           &! layer thick    at current time
         layerThicknessNew,           &! layer thick    at new     time
         highFreqThicknessCur,        &! high frq thick at current time
         highFreqThicknessNew,        &! high frq thick at new     time
         lowFreqDivergenceCur,        &! low frq div    at current time
         lowFreqDivergenceNew,        &! low frq div    at new     time
         tracerGroupIdealAgeMask       ! mask for resetting surface ideal age

      type (field1DReal), pointer :: &
         effectiveDensityField         ! field pointer for halo update

      ! State/Tracer arrays, info
      real (kind=RKIND), dimension(:,:,:), pointer ::      &!
         tracersGroupCur,      &! old, new tracer arrays
         tracersGroupNew

      integer :: &
         startIndex, endIndex, &! start, end index for tracer groups
         indexSalinity          ! tracer index for salinity

      integer :: iTracer, numTracers, numActiveTracers

      type (mpas_pool_iterator_type) :: &
         groupItr               ! iterator for tracer groups

      character (len=StrKIND) :: &
         modifiedGroupName,    &! constructed tracer group names
         configName             ! constructed config names for tracers

      ! Tendency Array Pointers

      real (kind=RKIND), dimension(:), pointer :: &
         sshTend                ! sea-surface height tendency

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocityTend,    &! normal velocity tendency
         highFreqThicknessTend, &! thickness tendency for high-freq iter
         lowFreqDivergenceTend, &! thickness tendency for low freq
         layerThicknessTend      ! main thickness tendency

      real (kind=RKIND), dimension(:,:,:), pointer :: &
         tracersGroupTend,      &! tendencies for tracer groups
         activeTracersTend       ! active tracer tendencies

      ! Forcing pool
      real (kind=RKIND), dimension(:), pointer :: tidalPotentialEta

      real (kind=RKIND) :: totalDepth
      real (kind=RKIND), dimension(:), pointer :: tidalInputMask, tidalBCValue
      real (kind=RKIND), dimension(:,:), pointer :: restingThickness

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersNew

      ! Remap variables
      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThicknessLagNew

      ! These are public variables used from the diagnostics module
      ! indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      ! indexSSHGradientZonal, indexSSHGradientMeridional
      ! barotropicForcing, barotropicThicknessFlux
      ! layerThickEdgeFlux, layerThickEdgeMean,
      ! normalTransportVelocity, normalGMBolusVelocity
      ! vertAleTransportTop, normalMLEvelocity
      ! velocityX, velocityY, velocityZ
      ! velocityZonal, velocityMeridional
      ! gradSSH, gradSSHX, gradSSHY, gradSSHZ
      ! gradSSHZonal, gradSSHMeridional
      ! surfaceVelocity, SSHGradient
      ! temperatureShortWaveTendency
      ! activeTracerHorizontalAdvectionTendency
      ! activeTracerVerticalAdvectionTendency
      ! activeTracerSurfaceFluxTendency
      ! activeTracerNonLocalTendency
      ! activeTracerHorMixTendency
      ! activeTracerHorizontalAdvectionEdgeFlux
      ! vertVelocityTop

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Start main time step timer
      call mpas_timer_start("se timestep")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Prep variables before first iteration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se prep")

      !*** Retrieve model state, pools
      !*** No longer support sub-blocks so this retrieves the only block

      block => domain%blocklist
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block%structs, 'verticalMesh', &
                                                 verticalMeshPool)
      call mpas_pool_get_subpool(block%structs, 'state', &
                                                 statePool)
      call mpas_pool_get_subpool(block%structs, 'forcing', &
                                                 forcingPool)
      call mpas_pool_get_subpool(block%structs, 'shortwave', &
                                                 swForcingPool)
      call mpas_pool_get_subpool(block%structs, 'tend', &
                                                 tendPool)
      call mpas_pool_get_subpool(block%structs, 'scratch', &
                                                 scratchPool)
      call mpas_pool_get_subpool(statePool, 'tracers', &
                                             tracersPool)
      call mpas_pool_get_subpool(tendPool,  'tracersTend', &
                                             tracersTendPool)

      !*** Retrieve state variables at two time levels

      call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                           normalBaroclinicVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                           normalBarotropicVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityCur, 1)

      call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                           normalBaroclinicVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                           normalBarotropicVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityNew, 2)

      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessNew, 2)

      call mpas_pool_get_array(statePool, 'highFreqThickness', &
                                           highFreqThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'highFreqThickness', &
                                           highFreqThicknessNew, 2)

      call mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                           lowFreqDivergenceCur, 1)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                           lowFreqDivergenceNew, 2)

      call mpas_pool_get_dimension_scalar(tracersPool, 'index_salinity', &
                                                        indexSalinity)

      !*** Retrieve tendency variables

      call mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                          highFreqThicknessTend)
      call mpas_pool_get_array(tendPool, 'normalVelocity', &
                                          normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'ssh', &
                                          sshTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', &
                                          layerThicknessTend)
      call mpas_pool_get_array(tendPool, 'normalVelocity', &
                                          normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                          highFreqThicknessTend)
      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', &
                                          lowFreqDivergenceTend)

      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersNew, 2)

      allocate(baroclinicThickness(nEdgesAll+1))

      if (config_transport_tests_flow_id > 0) then
        ! This is a transport test. Write advection velocity from prescribed
        ! flow field.
        call ocn_transport_test_velocity(meshPool, daysSinceStartOfSim, &
          0.5*dt, normalVelocityCur)
      endif

#ifdef MPAS_OPENACC
      !$acc enter data copyin(normalVelocityCur, normalBarotropicVelocityCur, &
      !$acc   sshCur, layerThicknessCur)
      !$acc enter data create(normalBaroClinicVelocityCur, normalVelocityNew, &
      !$acc   normalBaroclinicVelocityNew, sshNew, layerThicknessNew)
      if (associated(highFreqThicknessNew)) then
         !$acc enter data copyin(highFreqThicknessCur)
         !$acc enter data create(highFreqThicknessNew)
      endif
      if (associated(lowFreqDivergenceNew)) then
         !$acc enter data copyin(lowFreqDivergenceCur)
         !$acc enter data create(lowFreqDivergenceNew)
      endif
#endif

      ! Initialize * variables that are used to compute baroclinic
      ! tendencies below.

      ! The baroclinic velocity needs be recomputed at the beginning
      ! of a timestep because the implicit vertical mixing is
      ! conducted on the total u.  We keep normalBarotropicVelocity
      ! from the previous timestep.
      ! Note that normalBaroclinicVelocity may now include a
      ! barotropic component, because the weights layerThickness
      ! have changed.  That is OK, because the barotropicForcing
      ! variable subtracts out the barotropic component from the
      ! baroclinic.

#ifdef MPAS_OPENACC
      !$acc parallel present(normalVelocityCur, normalBarotropicVelocityCur, &
      !$acc    sshCur, layerThicknessCur, normalBaroClinicVelocityCur, &
      !$acc    normalVelocityNew, normalBaroclinicVelocityNew, sshNew, &
      !$acc    layerThicknessNew, minLevelCell, maxLevelCell)
#endif

#ifdef MPAS_OPENACC
      !$acc loop collapse(2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1,nEdgesAll
      do k = 1,nVertLevels

         normalBaroclinicVelocityCur(k,iEdge) = &
                                 normalVelocityCur(k,iEdge) - &
                         normalBarotropicVelocityCur(iEdge)

         normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)

         normalBaroclinicVelocityNew(k,iEdge) = &
         normalBaroclinicVelocityCur(k,iEdge)
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc loop gang
#else
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCellsAll
         sshNew(iCell) = sshCur(iCell)
#ifdef MPAS_OPENACC
         !$acc loop vector
#endif
         do k = 1,nVertLevels
            layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
         end do
      end do
#ifdef MPAS_OPENACC
      !$acc end parallel
#else
      !$omp end do
      !$omp end parallel
#endif

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr))
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr%memberName,&
                                     tracersGroupCur, 1)
            call mpas_pool_get_array(tracersPool, groupItr%memberName,&
                                     tracersGroupNew, 2)

            if ( associated(tracersGroupCur) .and. &
                 associated(tracersGroupNew) ) then

#ifdef MPAS_OPENACC
               !$acc enter data create(tracersGroupNew)
               !$acc enter data copyin(tracersGroupCur)
               !$acc parallel loop gang present(minLevelCell, maxLevelCell, &
               !$acc    tracersGroupNew, tracersGroupCur)
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k)
#endif
               do iCell = 1, nCellsAll
               do k = minLevelCell(iCell), maxLevelCell(iCell)
#ifdef MPAS_OPENACC
               !$acc loop vector
#endif
                  do iTracer = 1,size(tracersGroupCur,1)
                     tracersGroupNew(iTracer,k,iCell) = &
                     tracersGroupCur(iTracer,k,iCell)
                  end do
               end do
               end do
#ifdef MPAS_OPENACC
               !$acc exit data copyout(tracersGroupNew)
               !$acc exit data delete(tracersGroupCur)
#else
               !$omp end do
               !$omp end parallel
#endif
            end if
         end if
      end do

      if (associated(highFreqThicknessNew)) then
#ifdef MPAS_OPENACC
         !$acc parallel loop gang &
         !$acc    present(highFreqThicknessCur, highFreqThicknessNew)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
         !$acc loop vector
#endif
         do k = 1,nVertLevels
            highFreqThicknessNew(k,iCell) = &
            highFreqThicknessCur(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      if (associated(lowFreqDivergenceNew)) then
#ifdef MPAS_OPENACC
         !$acc parallel loop gang &
         !$acc    present(lowFreqDivergenceCur, lowFreqDivergenceNew)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
         !$acc loop vector
#endif
         do k = 1,nVertLevels
            lowFreqDivergenceNew(k,iCell) = &
            lowFreqDivergenceCur(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      endif

      call mpas_timer_stop("se prep")

      call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

      if (landIcePressureOn) then 
         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
      endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(normalVelocityCur, normalBarotropicVelocityCur, &
      !$acc    sshCur, layerThicknessCur)
      !$acc exit data copyout(normalBaroClinicVelocityCur, normalVelocityNew, &
      !$acc    normalBaroclinicVelocityNew, sshNew, layerThicknessNew)
      if (associated(highFreqThicknessNew)) then
         !$acc exit data delete(highFreqThicknessCur)
         !$acc exit data copyout(highFreqThicknessNew)
      endif
      if (associated(lowFreqDivergenceNew)) then
         !$acc exit data delete(lowFreqDivergenceCur)
         !$acc exit data copyout(lowFreqDivergenceNew)
      endif
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large outer timestep iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do splitExplicitStep = 1, numTSIterations

         if (config_disable_thick_all_tend .and. &
             config_disable_vel_all_tend .and. &
             config_disable_tr_all_tend) then
            exit ! don't compute in loop meant to update velocity,
                 ! thickness, and tracers
         end if

         call mpas_timer_start('se loop')

         stage1_tend_time = min(splitExplicitStep,2)

         call mpas_pool_get_array(statePool, 'normalVelocity', &
                           normalVelocityCur, stage1_tend_time)

#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThicknessCur, normalVelocityCur, sshCur)
         !$acc update device(layerThickEdgeFlux)
         if (config_use_freq_filtered_thickness .or. associated(highFreqThicknessNew) ) then
            !$acc enter data create(highFreqThicknessNew)
            !$acc enter data copyin(highFreqThicknessCur, highFreqThicknessTend)
         endif
#endif

         ! ---  update halos for diagnostic ocean boundary layer depth
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("se halo diag obd")
            call mpas_dmpar_field_halo_exch(domain,'boundaryLayerDepth')
            call mpas_timer_stop("se halo diag obd")
         end if

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("se halo diag")

         call mpas_dmpar_field_halo_exch(domain, &
                                 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
            call mpas_dmpar_field_halo_exch(domain, 'divergence')
            call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("se halo diag")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  Stage 1: Baroclinic velocity (3D) prediction, explicit with
         !           long timestep
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config_use_freq_filtered_thickness) then

            call mpas_timer_start("se freq-filtered-thick computations")

            call ocn_tend_freq_filtered_thickness(tendPool, statePool, &
                                                  stage1_tend_time)

            call mpas_timer_stop("se freq-filtered-thick computations")

            call mpas_timer_start("se freq-filtered-thick halo update")

            call mpas_dmpar_field_halo_exch(domain,'tendHighFreqThickness')
            call mpas_dmpar_field_halo_exch(domain,'tendLowFreqDivergence')

            call mpas_timer_stop("se freq-filtered-thick halo update")

#ifdef MPAS_OPENACC
            !$acc parallel loop gang present(minLevelCell, maxLevelCell, &
            !$acc    highFreqThicknessNew, highFreqThicknessCur, highFreqThicknessTend)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! this is h^{hf}_{n+1}
               highFreqThicknessNew(k,iCell) = &
               highFreqThicknessCur(k,iCell) + dt* &
               highFreqThicknessTend(k,iCell)
            end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

         endif ! freq filtered thickness

         ! compute velocity tendencies, T(u*,w*,p*)
         call mpas_timer_start("se bcl vel")
         call mpas_timer_start('se bcl vel tend')

         ! compute vertAleTransportTop.  Use u (rather than &
         ! normalTransportVelocity) for momentum advection.
         ! Use the most recent time level available.

         if (associated(highFreqThicknessNew)) then
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalVelocityCur, sshCur, dt, &
                 vertAleTransportTop, err, highFreqThicknessNew)
         else
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalVelocityCur, sshCur, dt, &
                 vertAleTransportTop, err)
         endif

#ifdef MPAS_OPENACC
         !$acc exit data delete(layerThicknessCur, normalVelocityCur, sshCur)
         !$acc update host(vertAleTransportTop)
         if (config_use_freq_filtered_thickness .or. associated(highFreqThicknessNew) ) then
            !$acc exit data copyout(highFreqThicknessNew)
            !$acc exit data delete(highFreqThicknessCur, highFreqThicknessTend)
         endif
#endif

         call ocn_tend_vel(domain, tendPool, statePool, forcingPool, &
                           stage1_tend_time, domain % dminfo, dt)

         call mpas_timer_stop('se bcl vel tend')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do j=1,numClinicIterations(splitExplicitStep)

            call mpas_timer_start('bcl iters on linear Coriolis')

            ! Put f*normalBaroclinicVelocity^{perp} in
            ! normalVelocityNew as a work variable
            call ocn_fuperp(statePool, meshPool, 2)

            !TODO: look at loop optimizations here
            ! Only need to loop over owned cells, since there is a halo
            ! exchange immediately after this computation.

            allocate(uTemp(nVertLevels,nEdgesOwned))

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k, cell1, cell2, uTemp, sshGrad, &
            !$omp         normalThicknessFluxSum, thicknessSum)
            do iEdge = 1, nEdgesOwned
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
               uTemp(:,iEdge) = 0.0_RKIND
               sshGrad = splitFact * gravity * &
                         (sshNew(cell2) - sshNew(cell1)) &
                         /dcEdge(iEdge)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)

                  ! normalBaroclinicVelocityNew =
                  ! normalBaroclinicVelocityOld +
                  !                dt*(-f*normalBaroclinicVelocityPerp
                  !                    + T(u*,w*,p*) + g*grad(SSH*) )
                  ! Here uNew is a work variable containing
                  !  -fEdge(iEdge)*normalBaroclinicVelocityPerp(k,iEdge)
                  uTemp(k,iEdge) = normalBaroclinicVelocityCur(k,iEdge) &
                           + dt * (normalVelocityTend(k,iEdge) &
                           + normalVelocityNew(k,iEdge) &
                           + sshGrad )
               enddo ! vertical

               ! thicknessSum is initialized outside the loop because
               ! on land boundaries maxLevelEdgeTop=0, but we want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               normalThicknessFluxSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)* &
                                        uTemp(minLevelEdgeBot(iEdge),iEdge)
               thicknessSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)

               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                 layerThickEdgeFlux(k,iEdge)*uTemp(k,iEdge)
                  thicknessSum = thicknessSum + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo
               barotropicForcing(iEdge) = splitFact* &
                              normalThicknessFluxSum/thicknessSum/dt

               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  ! These two steps are together here:
                  !{\bf u}'_{k,n+1} =
                  !    {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                  !{\bf u}'_{k,n+1/2} = \frac{1}{2}\left(
                  !        {\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right)
                  ! so that normalBaroclinicVelocityNew is at time n+1/2
                  normalBaroclinicVelocityNew(k,iEdge) = 0.5_RKIND*( &
                  normalBaroclinicVelocityCur(k,iEdge) + uTemp(k,iEdge) - &
                         dt * barotropicForcing(iEdge))

               enddo

            enddo ! iEdge
            !$omp end do
            !$omp end parallel

            baroclinicThickness(:) = 0.0_RKIND
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(cell1, cell2, k)
            do iEdge = 1, nEdgesHalo(config_num_halos+1)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               baroclinicThickness(iEdge) = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)
               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  baroclinicThickness(iEdge) = baroclinicThickness(iEdge) + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo
            enddo ! iEdge
            !$omp end do
            !$omp end parallel

            deallocate(uTemp)

            call mpas_timer_start("se halo normalBaroclinicVelocity")
            call mpas_dmpar_field_halo_exch(domain, &
                              'normalBaroclinicVelocity', timeLevel=2)
            call mpas_timer_stop("se halo normalBaroclinicVelocity")

            call mpas_timer_stop('bcl iters on linear Coriolis')

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_start('se halo barotropicForcing')
         call mpas_dmpar_field_halo_exch(domain, 'barotropicForcing')
         call mpas_timer_stop('se halo barotropicForcing')

         call mpas_timer_stop("se bcl vel")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Stage 2: Barotropic velocity prediction, explicitly subcycled
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("se btr vel")

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2

         ! unsplit_explicit option
         if (unsplit) then

            call mpas_timer_start('btr vel ue')

            ! For Split_Explicit unsplit, simply set
            !    normalBarotropicVelocityNew=0,
            !    normalBarotropicVelocitySubcycle=0, and
            !    uNew=normalBaroclinicVelocityNew

            !$omp parallel
            !$omp do schedule(runtime) private(k)
            do iEdge = 1, nEdgesAll
               normalBarotropicVelocityNew(iEdge) = 0.0_RKIND
               do k = 1, nVertLevels
                  normalVelocityNew(k, iEdge) = &
                     normalBaroclinicVelocityNew(k, iEdge)

                  ! This is u used in advective terms for layerThickness and tracers
                  ! in tendency calls in stage 3.
                  normalTransportVelocity(k,iEdge) = normalBaroclinicVelocityNew(k,iEdge)


               enddo ! vertical
            end do  ! iEdge
            !$omp end do
            !$omp end parallel

            !add GM and MLE (submesoscale) velocities if requested
            call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdgesAll, &
                                             nVertLevels)
            call ocn_MLE_add_to_transport_vel(normalTransportvelocity, nEdgesAll)

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k)
            do iEdge=1, nEdgesAll
               do k = 1, nVertLevels
                  normalTransportVelocity(k,iEdge) = edgeMask(k,iEdge)* &
                                  normalTransportVelocity(k,iEdge)
               end do
            end do
            !$omp end do
            !$omp end parallel

            call mpas_timer_stop('btr vel ue')

         else ! split explicit option

            ! Initialize variables for barotropic subcycling
            call mpas_timer_start('btr vel se init')

            call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                     sshSubcycleCur, oldBtrSubcycleTime)
            call mpas_pool_get_array(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleCur, &
                                  oldBtrSubcycleTime)

            if (config_filter_btr_mode) then
               !$omp parallel
               !$omp do schedule(runtime)
               do iEdge = 1, nEdgesAll
                  barotropicForcing(iEdge) = 0.0_RKIND
               end do
               !$omp end do
               !$omp end parallel
            endif

            !$omp parallel
            !$omp do schedule(runtime)
            do iCell = 1, nCellsAll
               ! sshSubcycleOld = sshOld
               sshSubcycleCur(iCell) = sshCur(iCell)
            end do
            !$omp end do

            !$omp do schedule(runtime)
            do iEdge = 1, nEdgesAll

               ! normalBarotropicVelocitySubcycleOld =
               ! normalBarotropicVelocityOld
               normalBarotropicVelocitySubcycleCur(iEdge) = &
               normalBarotropicVelocityCur(iEdge)

               ! normalBarotropicVelocityNew = BtrOld
               !  This is the first for the summation
               normalBarotropicVelocityNew(iEdge) = &
               normalBarotropicVelocityCur(iEdge)

               ! barotropicThicknessFlux = 0
               barotropicThicknessFlux(iEdge) = 0.0_RKIND
            end do
            !$omp end do
            !$omp end parallel

            call mpas_timer_stop('btr vel se init')

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! BEGIN Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Allocate subcycled scratch fields before starting
            ! subcycle loop
            allocate(btrvel_temp(nEdgesAll+1))
            btrvel_temp(:) = 0.0_RKIND
            allocate(wctEdge(nEdgesAll+1))
            wctEdge(:) = 0.0_RKIND
            allocate(deltaSSH(nCellsAll))
            deltaSSH(:) = 0.0_RKIND

            cellHaloComputeCounter = 0
            edgeHaloComputeCounter = 0

            call mpas_timer_start('btr se subcycle loop')
            do j = 1, nBtrSubcycles * config_btr_subcycle_loop_factor

               ! Update halos if needed
               if (cellHaloComputeCounter < neededHalos) then

                  call mpas_timer_start('se halo subcycle')
                  call mpas_dmpar_exch_group_reuse_halo_exch(domain, &
                       subcycleGroupName, timeLevel=oldBtrSubcycleTime)
                  call mpas_timer_stop('se halo subcycle')

                  ! Reset halo counters after halo update
                  cellHaloComputeCounter = config_num_halos - &
                                           haloDecrement
                  edgeHaloComputeCounter = config_num_halos + 1 - &
                                           haloDecrement
               end if ! halo update

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: VELOCITY PREDICTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! only do this part if it is needed in next SSH solve
               if (config_btr_gam1_velWt1 > 1.0e-12_RKIND) then
                  uPerpTime = oldBtrSubcycleTime

                  call mpas_pool_get_array(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleCur, &
                                  uPerpTime)
                  call mpas_pool_get_array(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleNew, &
                                  newBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                  sshSubcycleCur, oldBtrSubcycleTime)

                  ! TODO: Modify this to eliminate pointer reassigning
                  ! Subtract tidal potential from ssh, if needed
                  ! Subtract the tidal potential from the current
                  !    subcycle ssh and store and a work array.
                  ! Then point sshSubcycleCur to the work array so the
                  ! tidal potential terms are included in the grad
                  ! operator inside the edge loop.
                  if (config_use_tidal_potential_forcing) then
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleCurWithTides', &
                                            sshSubcycleCurWithTides)
                     call mpas_pool_get_array(forcingPool, &
                                             'tidalPotentialEta', &
                                              tidalPotentialEta)

                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell = 1, nCellsAll
                        sshSubcycleCurWithTides(iCell) = &
                                 sshSubcycleCur(iCell) - &
                              tidalPotentialEta(iCell) - &
                           pgf_sal_on * pgf_sal(iCell) - &
                           (1.0_RKIND - pgf_sal_on) * self_attraction_and_loading_beta* &
                                 sshSubcycleCur(iCell)
                     end do
                     !$omp end do
                     !$omp end parallel

                     call mpas_pool_get_array(forcingPool, &
                                          'sshSubcycleCurWithTides', &
                                           sshSubcycleCur)

                  end if ! tidal potential forcing

                  ! direct application of tidal boundary condition
                  if (config_use_tidal_forcing .and. trim(config_tidal_forcing_type) == 'direct') then
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleCurWithTides', &
                                            sshSubcycleCurWithTides)
                     call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
                     call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)
                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell=1, nCellsAll
                        ! boolean mask for now, could generalize to tappered sponge layer
                        ! only modify layer thicknesses on tidal boundary
                        do k = minLevelCell(iCell), maxLevelCell(iCell)
                           sshSubcycleCurWithTides(iCell) = sshSubcycleCur(iCell) + &
                               (tidalInputMask(iCell)*(tidalBCValue(iCell) - sshSubcycleCur(iCell)))
                        end do
                     end do
                     !$omp end do
                     !$omp end parallel
                     call mpas_pool_get_array(forcingPool, &
                                          'sshSubcycleCurWithTides', &
                                           sshSubcycleCur)
                  end if
                  if (edgeHaloComputeCounter <= 1) then
                     nEdges = nEdgesOwned
                  else
                     nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
                  endif

                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(i, cell1, cell2, eoe, CoriolisTerm, &
                  !$omp         temp_mask)
                  do iEdge = 1, nEdges

                     temp_mask = edgeMask(1, iEdge)

                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0_RKIND
                     do i = 1, nEdgesOnEdge(iEdge)
                        eoe = edgesOnEdge(i,iEdge)
                        CoriolisTerm = CoriolisTerm + &
                                       weightsOnEdge(i,iEdge)* &
                             normalBarotropicVelocitySubcycleCur(eoe)* &
                                       fEdge(eoe)
                     end do

                     normalBarotropicVelocitySubcycleNew(iEdge) = &
                           temp_mask* &
                           (normalBarotropicVelocitySubcycleCur(iEdge) &
                            + dt/nBtrSubcycles*(CoriolisTerm - gravity &
                             *(sshSubcycleCur(cell2) - &
                               sshSubcycleCur(cell1))/dcEdge(iEdge) + &
                               barotropicForcing(iEdge)))

                  end do ! edges
                  !$omp end do
                  !$omp end parallel

               endif ! config_btr_gam1_velWt1>1.0e-12

               ! Halo from edges is corrupted here, so reduce the edge halo
               ! counter by 1
               edgeHaloComputeCounter = edgeHaloComputeCounter - 1

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: SSH PREDICTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                    sshSubcycleCur, &
                                                    oldBtrSubcycleTime)
               call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                    sshSubcycleNew, &
                                                    newBtrSubcycleTime)
               call mpas_pool_get_array(statePool, &
                             'normalBarotropicVelocitySubcycle', &
                              normalBarotropicVelocitySubcycleCur, &
                              oldBtrSubcycleTime)
               call mpas_pool_get_array(statePool, &
                             'normalBarotropicVelocitySubcycle', &
                              normalBarotropicVelocitySubcycleNew, &
                              newBtrSubcycleTime)

               if (cellHaloComputeCounter <= 1) then
                  nCells = nCellsOwned
               else
                  nCells = nCellsHalo(cellHaloComputeCounter-1)
               endif
               if (edgeHaloComputeCounter <= 1) then
                  nEdges = nEdgesOwned
               else
                  nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
               endif

               ! config_btr_gam1_velWt1 sets the forward weighting of
               !    velocity in the SSH computation
               ! config_btr_gam1_velWt1=  1
               !    flux = normalBarotropicVelocityNew*H
               ! config_btr_gam1_velWt1=0.5
               !    flux = 1/2*(normalBarotropicVelocityNew +
               !                normalBarotropicVelocityOld)*H
               ! config_btr_gam1_velWt1=  0
               !    flux = normalBarotropicVelocityOld*H

               !$omp parallel
               !$omp do schedule(runtime)
               do iEdge = 1, nEdgesAll
                  btrvel_temp(iEdge) = (1.0_RKIND - config_btr_gam1_velWt1) * &
                                       normalBarotropicVelocitySubcycleCur(iEdge) + &
                                       config_btr_gam1_velWt1 * &
                                       normalBarotropicVelocitySubcycleNew(iEdge)
               end do
               !$omp end do
               !$omp end parallel

               ! Compute barotropic thickness at the edge. Note this
               ! matches the sum of baroclinic thicknesses at this edge.
               !$omp parallel
               !$omp do schedule(runtime)
               do iCell = 1, nCellsAll
                  deltaSSH(iCell) = sshSubcycleCur(iCell) - sshNew(iCell)
               end do
               !$omp end do
               !$omp end parallel
               call ocn_diagnostic_solve_wctEdge(btrvel_temp, deltaSSH, &
                         baroclinicThickness, wctEdge)

               !$omp parallel
               !$omp do schedule(runtime) &
               !$omp private(i, iEdge, flux)
               do iCell = 1, nCells
                  sshTend(iCell) = 0.0_RKIND
                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
                     if (maxLevelEdgeTop(iEdge).eq.0) cycle

                     flux = btrvel_temp(iEdge) * wctEdge(iEdge)

                     sshTend(iCell) = sshTend(iCell) + &
                                      edgeSignOncell(i, iCell)*flux* &
                                      dvEdge(iEdge)

                  end do ! edges on cell

                  ! SSHnew = SSHold + dt/J*(-div(Flux))
                  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) &
                                        + dt/nBtrSubcycles* &
                                          sshTend(iCell)/areaCell(iCell)
               end do ! cells
               !$omp end do
               !$omp end parallel

               ! avoid redundant computations when
               ! config_btr_solve_SSH2 is true
               if (config_btr_solve_SSH2) then

                  ! If config_btr_solve_SSH2=.true.,
                  ! then do NOT accumulate barotropicThicknessFlux
                  ! in this SSH predictor section, because it will be
                  ! accumulated in the SSH corrector section.

               else

                  ! otherwise, DO accumulate barotropicThicknessFlux
                  ! in this SSH predictor section

                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(flux)
                  do iEdge = 1, nEdges
                     if (maxLevelEdgeTop(iEdge).eq.0) cycle

                     flux = btrvel_temp(iEdge) * wctEdge(iEdge)

                     barotropicThicknessFlux(iEdge) = &
                     barotropicThicknessFlux(iEdge) + flux
                  end do ! edges
                  !$omp end do
                  !$omp end parallel

               endif

               ! a cell halo layer is now corrupted - decrement counter
               cellHaloComputeCounter = cellHaloComputeCounter - 1

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: VELOCITY CORRECTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               do BtrCorIter = 1, config_n_btr_cor_iter
                  uPerpTime = newBtrSubcycleTime

                  call mpas_pool_get_array(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleCur, &
                                  oldBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleNew, &
                                  newBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                       sshSubcycleCur,&
                                                    oldBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                       sshSubcycleNew,&
                                                    newBtrSubcycleTime)

                  ! TODO: change to avoid pointer reassignment
                  ! Subtract tidal potential from ssh, if needed
                  ! Subtract the tidal potential from the current and
                  !    new subcycle ssh and store in work arrays.
                  ! Then point sshSubcycleCur and ssh SubcycleNew to
                  ! the work arrays so the tidal potential terms are
                  ! included in the grad operator inside the edge loop.
                  if (config_use_tidal_potential_forcing) then
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleCurWithTides',  &
                                            sshSubcycleCurWithTides)
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleNewWithTides',  &
                                            sshSubcycleNewWithTides)
                     call mpas_pool_get_array(forcingPool,  &
                                             'tidalPotentialEta',  &
                                              tidalPotentialEta)

                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell = 1, nCellsAll
                        sshSubcycleCurWithTides(iCell) = &
                                 sshSubcycleCur(iCell) - &
                              tidalPotentialEta(iCell) - &
                           pgf_sal_on * pgf_sal(iCell) - &
                           (1.0_RKIND - pgf_sal_on) * self_attraction_and_loading_beta* &
                                 sshSubcycleCur(iCell)

                        sshSubcycleNewWithTides(iCell) = &
                                 sshSubcycleNew(iCell) - &
                              tidalPotentialEta(iCell) - &
                           pgf_sal_on * pgf_sal(iCell) - &
                           (1.0_RKIND - pgf_sal_on) * self_attraction_and_loading_beta* &
                                 sshSubcycleNew(iCell)
                     end do
                     !$omp end do
                     !$omp end parallel

                     call mpas_pool_get_array(forcingPool, &
                                    'sshSubcycleCurWithTides',  &
                                     sshSubcycleCur)
                     call mpas_pool_get_array(forcingPool, &
                                    'sshSubcycleNewWithTides',  &
                                     sshSubcycleNew)
                  end if ! tidal potential forcing
                  ! direct application of tidal boundary condition
                  if (config_use_tidal_forcing .and. trim(config_tidal_forcing_type) == 'direct') then
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleCurWithTides',  &
                                            sshSubcycleCurWithTides)
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleNewWithTides',  &
                                            sshSubcycleNewWithTides)
                     call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
                     call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)
                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell=1, nCellsAll
                        ! boolean mask for now, could generalize to tappered sponge layer
                        ! only modify layer thicknesses on tidal boundary
                        do k = minLevelCell(iCell), maxLevelCell(iCell)
                           sshSubcycleCurWithTides(iCell) = sshSubcycleCur(iCell) + &
                               (tidalInputMask(iCell)*(tidalBCValue(iCell) - sshSubcycleCur(iCell)))
                           sshSubcycleNewWithTides(iCell) = sshSubcycleNew(iCell) + &
                               (tidalInputMask(iCell)*(tidalBCValue(iCell) - sshSubcycleNew(iCell)))
                        end do
                     end do
                     !$omp end do
                     !$omp end parallel
                     call mpas_pool_get_array(forcingPool, &
                                    'sshSubcycleCurWithTides',  &
                                     sshSubcycleCur)
                     call mpas_pool_get_array(forcingPool, &
                                    'sshSubcycleNewWithTides',  &
                                     sshSubcycleNew)
                  end if

                  ! Need to initialize btr_vel_temp over one more halo
                  ! than we are computing over
                  nEdges = nEdgesHalo(min(edgeHaloComputeCounter, &
                                          config_num_halos))

                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdges+1
                     btrvel_temp(iEdge) = &
                         normalBarotropicVelocitySubcycleNew(iEdge)
                  end do
                  !$omp end do
                  !$omp end parallel

                  if (edgeHaloComputeCounter <= 1) then
                     nEdges = nEdgesOwned
                  else
                     nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
                  endif

                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(i, cell1, cell2, eoe, coriolisTerm, &
                  !$omp         sshCell1, sshCell2, temp_mask)
                  do iEdge = 1, nEdges

                     temp_mask = edgeMask(1,iEdge)

                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0_RKIND
                     do i = 1, nEdgesOnEdge(iEdge)
                        eoe = edgesOnEdge(i,iEdge)
                        CoriolisTerm = CoriolisTerm &
                                     + weightsOnEdge(i,iEdge)* &
                                       btrvel_temp(eoe)*fEdge(eoe)
                     end do

                     ! In this final solve for velocity, SSH is a
                     ! linear combination of SSHold and SSHnew.
                     sshCell1 = (1 - config_btr_gam2_SSHWt1) * &
                                         sshSubcycleCur(cell1) &
                              +      config_btr_gam2_SSHWt1  * &
                                         sshSubcycleNew(cell1)
                     sshCell2 = (1 - config_btr_gam2_SSHWt1) * &
                                         sshSubcycleCur(cell2) &
                              +      config_btr_gam2_SSHWt1  * &
                                         sshSubcycleNew(cell2)

                     ! normalBarotropicVelocityNew =
                     ! normalBarotropicVelocityOld + dt/J*
                     !     (-f*normalBarotropicVelocityoldPerp
                     !      - g*grad(SSH) + G)
                     normalBarotropicVelocitySubcycleNew(iEdge) = &
                          temp_mask * &
                          (normalBarotropicVelocitySubcycleCur(iEdge) &
                           + dt/nBtrSubcycles &
                           *(CoriolisTerm - gravity* &
                             (sshCell2 - sshCell1)/dcEdge(iEdge) &
                           + barotropicForcing(iEdge)))

                  end do
                  !$omp end do
                  !$omp end parallel

                  edgeHaloComputeCounter = edgeHaloComputeCounter - 1
                  if (BtrCorIter >= 1 .or. &
                      .not. config_btr_solve_SSH2) &
                     cellHaloComputeCounter = cellHaloComputeCounter-1

               end do !do BtrCorIter=1,config_n_btr_cor_iter

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: SSH CORRECTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               if (config_btr_solve_SSH2) then

                  allocate(sshTemp(nCellsAll))

                  call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                       sshSubcycleCur,&
                                                    oldBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                       sshSubcycleNew,&
                                                    newBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, &
                                  'normalBarotropicVelocitySubcycle', &
                                   normalBarotropicVelocitySubcycleCur, &
                                   oldBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, &
                                  'normalBarotropicVelocitySubcycle', &
                                   normalBarotropicVelocitySubcycleNew, &
                                   newBtrSubcycleTime)


                  if (cellHaloComputeCounter <= 1) then
                     nCells = nCellsOwned
                  else
                     nCells = nCellsHalo(cellHaloComputeCounter-1)
                  endif
                  if (edgeHaloComputeCounter <= 1) then
                     nEdges = nEdgesOwned
                  else
                     nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
                  endif

                  ! config_btr_gam3_velWt2 sets the forward weighting
                  ! of velocity in the SSH computation
                  ! config_btr_gam3_velWt2=  1
                  !    flux = normalBarotropicVelocityNew*H
                  ! config_btr_gam3_velWt2=0.5
                  !    flux = 1/2*(normalBarotropicVelocityNew +
                  !                normalBarotropicVelocityOld)*H
                  ! config_btr_gam3_velWt2=  0
                  !    flux = normalBarotropicVelocityOld*H

                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iCell = 1, nCellsAll
                     ! SSH is linear combination of SSHold and SSHnew
                     sshTemp(iCell) = (1-config_btr_gam2_SSHWt1)* &
                                      sshSubcycleCur(iCell) + &
                                      config_btr_gam2_SSHWt1 * &
                                      sshSubcycleNew(iCell)
                     deltaSSH(iCell) = sshTemp(iCell) - sshNew(iCell)
                  end do ! cell loop for ssh
                  !$omp end do
                  !$omp end parallel

                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdgesAll
                     btrvel_temp(iEdge) = (1.0-config_btr_gam3_velWt2) * &
                                          normalBarotropicVelocitySubcycleCur(iEdge) + &
                                          config_btr_gam3_velWt2 * &
                                          normalBarotropicVelocitySubcycleNew(iEdge)
                  end do
                  !$omp end do
                  !$omp end parallel

                  call ocn_diagnostic_solve_wctEdge(btrvel_temp, deltaSSH, &
                         baroclinicThickness, wctEdge)

                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(i, iEdge, flux)
                  do iCell = 1, nCells
                     sshTend(iCell) = 0.0_RKIND
                     do i = 1, nEdgesOnCell(iCell)
                        iEdge = edgesOnCell(i, iCell)

                        if (maxLevelEdgeTop(iEdge).eq.0) cycle

                        flux = btrvel_temp(iEdge) * wctEdge(iEdge)

                        sshTend(iCell) = sshTend(iCell) + &
                                         edgeSignOnCell(i, iCell)* &
                                         flux*dvEdge(iEdge)

                     end do ! edges on cell

                     ! SSHnew = SSHold + dt/J*(-div(Flux))
                     sshSubcycleNew(iCell) = sshSubcycleCur(iCell) &
                                        + dt/nBtrSubcycles * &
                                          sshTend(iCell)/areaCell(iCell)
                  end do ! cell loop for ssh
                  !$omp end do
                  !$omp end parallel

                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iCell = 1, nCellsAll
                     ! SSH is linear combination of SSHold and SSHnew
                     sshTemp(iCell) = (1-config_btr_gam2_SSHWt1)* &
                                      sshSubcycleCur(iCell) + &
                                      config_btr_gam2_SSHWt1 * &
                                      sshSubcycleNew(iCell)
                     deltaSSH(iCell) = sshTemp(iCell) - sshNew(iCell)
                  end do ! cell loop for ssh
                  !$omp end do
                  !$omp end parallel

                  call ocn_diagnostic_solve_wctEdge(btrvel_temp, deltaSSH, &
                         baroclinicThickness, wctEdge)

                  ! compute barotropic thickness flux on edges
                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(flux)
                  do iEdge = 1, nEdges

                     flux = btrvel_temp(iEdge) * wctEdge(iEdge)

                     barotropicThicknessFlux(iEdge) = &
                     barotropicThicknessFlux(iEdge) + flux

                  end do ! edge barotropic thickness flux
                  !$omp end do
                  !$omp end parallel

                  edgeHaloComputeCounter = config_num_halos + 1

                  deallocate(sshTemp)

               endif ! config_btr_solve_SSH2

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: Accumulate running sums, advance
               !                      timestep pointers
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               call mpas_pool_get_array(statePool, &
                                  'normalBarotropicVelocitySubcycle', &
                                   normalBarotropicVelocitySubcycleNew, &
                                   newBtrSubcycleTime)

               ! normalBarotropicVelocityNew =
               ! normalBarotropicVelocityNew +
               ! normalBarotropicVelocitySubcycleNEW
               ! This accumulates the sum. If the Barotropic Coriolis
               ! iteration is limited to one, this could
               ! be merged with the above code.

               !$omp parallel
               !$omp do schedule(runtime)
               do iEdge = 1, nEdgesOwned
                  normalBarotropicVelocityNew(iEdge) = &
                  normalBarotropicVelocityNew(iEdge) + &
                  normalBarotropicVelocitySubcycleNew(iEdge)
               end do  ! iEdge
               !$omp end do
               !$omp end parallel

               ! advance time pointers
               oldBtrSubcycleTime = mod(oldBtrSubcycleTime,2)+1
               newBtrSubcycleTime = mod(newBtrSubcycleTime,2)+1

            end do ! j=1,nBtrSubcycles
            call mpas_timer_stop('btr se subcycle loop')

            deallocate(btrvel_temp)
            deallocate(wctEdge)
            deallocate(deltaSSH)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Normalize Barotropic subcycle sums: ssh,
            ! normalBarotropicVelocity, and F
            call mpas_timer_start('btr se norm')

            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdgesOwned
               barotropicThicknessFlux(iEdge) = &
               barotropicThicknessFlux(iEdge) &
                      /(nBtrSubcycles*config_btr_subcycle_loop_factor)

               normalBarotropicVelocityNew(iEdge) = &
               normalBarotropicVelocityNew(iEdge) &
                     / (nBtrSubcycles*config_btr_subcycle_loop_factor + 1)
            end do
            !$omp end do
            !$omp end parallel

            call mpas_timer_stop('btr se norm')

            ! boundary update on F
            call mpas_timer_start("se halo F and btr vel")
            call mpas_dmpar_exch_group_create(domain, finalBtrGroupName)

            call mpas_dmpar_exch_group_add_field(domain, &
                                              finalBtrGroupName, &
                                             'barotropicThicknessFlux')
            call mpas_dmpar_exch_group_add_field(domain, &
                                              finalBtrGroupName, &
                                             'normalBarotropicVelocity', &
                                              timeLevel=2)

            call mpas_dmpar_exch_group_full_halo_exch(domain, &
                                                      finalBtrGroupName)

            call mpas_dmpar_exch_group_destroy(domain,finalBtrGroupName)
            call mpas_timer_stop("se halo F and btr vel")

            ! Check that you can compute SSH using the total sum or the
            ! individual increments over the barotropic subcycles.
            ! efficiency: This next block of code is really a check for
            ! debugging, and can be removed later.
            call mpas_timer_start('btr se ssh verif')

            ! Correction velocity
            !   normalVelocityCorrection = (Flux - Sum(h u*))/H
            ! or, for the full latex version:
            !{\bf u}^{corr} = \left( {\overline {\bf F}}
            !  - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}
            ! {\bf u}_k^{avg} \right)
            ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}   \right.

            nEdges = nEdgesHalo(config_num_halos-1 )
            allocate(uTemp(nVertLevels,nEdges))

            !Compute uTemp
            ! velocity for normalVelocityCorrection is
            ! normalBarotropicVelocity +
            ! normalBaroclinicVelocity + uBolus
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               uTemp(:,iEdge) = normalBarotropicVelocityNew(iEdge) &
                            + normalBaroclinicVelocityNew(:,iEdge)
            end do
            !$omp end do
            !$omp end parallel

            call ocn_GM_add_to_transport_vel(uTemp, nEdges, &
                                             nVertLevels)
            call ocn_MLE_add_to_transport_vel(uTemp, nEdges)

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k)
            do iEdge = 1, nEdges
               do k = 1, nVertLevels
                  normalTransportVelocity(k,iEdge) =  &
                             normalBarotropicVelocityNew(iEdge)   + &
                             normalBaroclinicVelocityNew(k,iEdge)
               end do
            end do
            !$omp end do
            !$omp end parallel
            call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdges, &
                                             nVertLevels)
            call ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(cell1, cell2, k, normalThicknessFluxSum, &
            !$omp         thicknessSum, normalVelocityCorrection)
            do iEdge = 1, nEdges

               if (thickEdgeFluxChoice == thickEdgeFluxUpwind) then
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)
                  do k=1,nVertLevels
                     ! initialize layerThickEdgeFlux to avoid divide by
                     ! zero and NaN problems.
                     layerThickEdgeFlux(k,iEdge) = -1.0e34_RKIND
                  end do
                  do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                     if (uTemp(k,iEdge) > 0.0_RKIND) then
                        layerThickEdgeFlux(k,iEdge) = layerThicknessCur(k,cell1)
                     elseif (uTemp(k,iEdge) < 0.0_RKIND) then
                        layerThickEdgeFlux(k,iEdge) = layerThicknessCur(k,cell2)
                     else
                        layerThickEdgeFlux(k,iEdge) = &
                                            max(layerThicknessCur(k,cell1), &
                                                layerThicknessCur(k,cell2))
                     end if
                  end do
               endif
               ! thicknessSum is initialized outside the loop because
               ! on land boundaries maxLevelEdgeTop=0, but I want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               normalThicknessFluxSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)* &
                                        uTemp(minLevelEdgeBot(iEdge),iEdge)
               thicknessSum  = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)

               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                           layerThickEdgeFlux(k,iEdge)*&
                                           uTemp(k,iEdge)
                  thicknessSum = thicknessSum + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo

               normalVelocityCorrection = useVelocityCorrection* &
                             ((barotropicThicknessFlux(iEdge) - &
                               normalThicknessFluxSum)/thicknessSum)

               do k = 1, nVertLevels

                  ! normalTransportVelocity = normalBarotropicVelocity
                  !                         + normalBaroclinicVelocity
                  !                         + normalGMBolusVelocity
                  !                         + normalMLEvelocity 
                  !                         + normalVelocityCorrection
                  ! This is the velocity used in advective terms for layerThickness
                  ! and tracers in tendency calls in stage 3.
                  normalTransportVelocity(k,iEdge) = edgeMask(k,iEdge)*( &
                             normalTransportVelocity(k,iEdge) + &
                             normalVelocityCorrection)
               enddo

            end do ! iEdge
            !$omp end do
            !$omp end parallel

            deallocate(uTemp)

            call mpas_timer_stop('btr se ssh verif')

         endif ! split_explicit

         call mpas_timer_stop("se btr vel")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Stage 3: Tracer, density, pressure, vert velocity prediction
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! only compute tendencies for active tracers on last large iteration
         if (splitExplicitStep < numTSIterations) then
            activeTracersOnly = .true.
         else
            activeTracersOnly = .false.
         endif

         ! Thickness tendency computations and thickness halo updates
         ! are completed before tracer tendency computations to allow
         ! monotonic advection.

         call mpas_timer_start('se thick tend')

         ! compute vertAleTransportTop.  Use normalTransportVelocity
         ! for advection of layerThickness and tracers.
         ! Use time level 1 values of layerThickness and
         ! layerThickEdgeFlux because layerThickness has not yet
         ! been computed for time level 2.
         call mpas_timer_start('thick vert trans vel top')
#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThicknessCur, sshCur)
         !$acc update device(layerThickEdgeFlux, normalTransportVelocity)
#endif
         if (associated(highFreqThicknessNew)) then
#ifdef MPAS_OPENACC
            !$acc enter data copyin(highFreqThicknessNew)
#endif
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalTransportVelocity, sshCur, &
                 dt, vertAleTransportTop, err, highFreqThicknessNew)
#ifdef MPAS_OPENACC
            !$acc exit data delete(highFreqThicknessNew)
#endif
         else
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalTransportVelocity, sshCur, &
                 dt, vertAleTransportTop, err)
         endif
#ifdef MPAS_OPENACC
         !$acc exit data delete(layerThicknessCur, sshCur)
         !$acc update host(vertAleTransportTop, normalTransportVelocity)
#endif
         call mpas_timer_stop('thick vert trans vel top')

         call ocn_tend_thick(tendPool, forcingPool)

         call mpas_timer_stop('se thick tend')

         ! update halo for thickness tendencies
         call mpas_timer_start("se halo thickness")

         call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

         call mpas_timer_stop("se halo thickness")

         call mpas_timer_start('se tracer tend', .false.)
         call ocn_tend_tracer(tendPool, statePool, forcingPool, &
                              meshPool, swForcingPool, &
                              dt, activeTracersOnly, 2)

         call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                                tracersGroupCur, 1)
         call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                                tracersGroupNew, 2)
         call mpas_pool_get_array(statePool,   'normalVelocity', &
                                                normalVelocityCur, 1)
         call mpas_pool_get_array(tracersTendPool,'activeTracersTend', &
                                                   activeTracersTend)

#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
         !$acc enter data copyin(normalBarotropicVelocityNew, &
         !$acc                   normalBaroclinicVelocityNew)
         !$acc enter data copyin(activeTracersNew)
         !$acc update device(layerThickEdgeFlux)
         if (config_use_freq_filtered_thickness) then
            !$acc enter data copyin(lowFreqDivergenceNew)
            !$acc enter data copyin(lowFreqDivergenceCur)
            !$acc enter data copyin(lowFreqDivergenceTend)
         endif
         if (splitExplicitStep < numTSIterations) then
            !$acc update device (normalTransportVelocity)
            !$acc enter data copyin(atmosphericPressure, seaIcePressure)
            !$acc enter data copyin(sshNew)
            !$acc update device(tracersSurfaceValue)
            if ( associated(normalGMBolusVelocity) ) then
               !$acc update device (normalGMBolusVelocity)
            end if
            if ( associated(normalMLEVelocity) ) then
               !$acc update device (normalMLEvelocity)
            end if
            if ( associated(frazilSurfacePressure) ) then
               !$acc enter data copyin(frazilSurfacePressure)
            endif
            if (landIcePressureOn) then
               !$acc enter data copyin(landIcePressure)
               !$acc enter data copyin(landIceDraft)
            endif
            if (config_use_freq_filtered_thickness) then
               !$acc enter data copyin(highFreqThicknessNew)
               !$acc enter data copyin(highFreqThicknessCur)
            endif
         elseif (splitExplicitStep == numTSIterations) then
            !$acc enter data copyin(layerThicknessCur)
            !$acc enter data copyin(layerThicknessTend)
            !$acc enter data copyin(normalBaroclinicVelocityCur)
            if (config_compute_active_tracer_budgets) then
               !$acc update device(activeTracerHorizontalAdvectionEdgeFlux)
               !$acc update device(activeTracerHorizontalAdvectionTendency)
               !$acc update device(activeTracerVerticalAdvectionTendency)
               !$acc update device(activeTracerHorMixTendency)
               !$acc update device(activeTracerSurfaceFluxTendency)
               !$acc update device(temperatureShortWaveTendency)
               !$acc update device(activeTracerNonLocalTendency)
            endif
         endif
#endif

         call mpas_timer_stop('se tracer tend')

         ! update halo for tracer tendencies
         call mpas_timer_start("se halo tracers")

         call mpas_pool_begin_iteration(tracersTendPool)
         do while (mpas_pool_get_next_member(tracersTendPool, &
                                             groupItr) )
            if (groupItr%memberType == MPAS_POOL_FIELD ) then
               ! Only compute tendencies for active tracers if
               ! activeTracersOnly flag is true.
               if (.not. activeTracersOnly .or. &
                   trim(groupItr%memberName)=='activeTracersTend') then
                  call mpas_dmpar_field_halo_exch(domain, &
                                                  groupItr%memberName)
               end if
            end if
         end do
         call mpas_timer_stop("se halo tracers")

         call mpas_timer_start('se loop fini')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  If iterating, reset variables for next iteration
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (splitExplicitStep < numTSIterations) then

            ! Get indices for dynamic tracers (Includes T&S).
            call mpas_pool_get_dimension_scalar(tracersPool,'activeGRP_start',&
                                                             startIndex)
            call mpas_pool_get_dimension_scalar(tracersPool,'activeGRP_end', &
                                                             endIndex)

            ! Only need T & S for earlier iterations,
            ! then all the tracers needed the last time through.

!! Keeping this calc on the CPU for now
#ifdef MPAS_OPENACC
            !$acc update host(layerThicknessNew)
            !$acc update host(activeTracersNew)
            !$acc update host(tracersSurfaceValue)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(i, k, temp_h, temp)
#endif
            do iCell = 1, nCellsAll
            do k = minLevelCell(iCell), maxLevelCell(iCell)

               ! this is h_{n+1}
               temp_h = layerThicknessCur(k,iCell) + dt* &
               layerThicknessTend(k,iCell)

               ! this is h_{n+1/2}
               layerThicknessNew(k,iCell) = 0.5* &
               (layerThicknessCur(k,iCell) + temp_h)

               do i = startIndex, endIndex
                  ! This is Phi at n+1
                  temp = (tracersGroupCur(i,k,iCell)* &
                          layerThicknessCur(k,iCell) + dt* &
                          activeTracersTend(i,k,iCell))/temp_h

                  ! This is Phi at n+1/2
                  tracersGroupNew(i,k,iCell) = 0.5_RKIND* &
                            (tracersGroupCur(i,k,iCell) + temp)
               end do ! tracer index
            end do ! vertical
            end do ! iCell
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

#ifdef MPAS_OPENACC
            !$acc update device(layerThicknessNew)
            !$acc update device(activeTracersNew)
#endif

            if (config_use_freq_filtered_thickness) then

#ifdef MPAS_OPENACC
               !$acc parallel loop present(minLevelCell, maxLevelCell) &
               !$acc   present(highFreqThicknessNew) &
               !$acc   present(highFreqThicknessCur) &
               !$acc   present(lowFreqDivergenceNew) &
               !$acc   present(lowFreqDivergenceCur) &
               !$acc   present(lowFreqDivergenceTend)
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k, temp)
#endif
               do iCell = 1, nCellsAll
               do k = minLevelCell(iCell), maxLevelCell(iCell)

                  ! h^{hf}_{n+1} was computed in Stage 1

                  ! this is h^{hf}_{n+1/2}
                  highFreqThicknessNew(k,iCell) = 0.5_RKIND* &
                             (highFreqThicknessCur(k,iCell) + &
                              highFreqThicknessNew(k,iCell))

                  ! this is D^{lf}_{n+1}
                  temp = lowFreqDivergenceCur(k,iCell) + dt* &
                         lowFreqDivergenceTend(k,iCell)

                  ! this is D^{lf}_{n+1/2}
                  lowFreqDivergenceNew(k,iCell) = 0.5_RKIND* &
                           (lowFreqDivergenceCur(k,iCell) + temp)
               end do
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif
            end if

#ifdef MPAS_OPENACC
            !$acc parallel loop collapse(2) present(normalVelocityNew) &
            !$acc   present(normalBaroclinicVelocityNew, edgeMask) &
            !$acc   present(normalBarotropicVelocityNew)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iEdge = 1, nEdgesAll
            do k = 1, nVertLevels

               ! u = normalBarotropicVelocity + normalBaroclinicVelocity
               ! here normalBaroclinicVelocity is at time n+1/2
               ! This is u used in next iteration or step
               normalVelocityNew(k,iEdge) = edgeMask(k,iEdge)* &
                         ( normalBarotropicVelocityNew(iEdge) + &
                           normalBaroclinicVelocityNew(k,iEdge) )

            enddo
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

            ! Efficiency note: We really only need this to compute
            ! layerThickEdgeFlux, density, pressure, and SSH
            ! in this diagnostics solve.
            call ocn_diagnostic_solve(dt, statePool, forcingPool, &
                                      meshPool, verticalMeshPool, scratchPool, &
                                      tracersPool, 2, full=.false.)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If large iteration complete, compute all variables at
         ! time n+1
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         elseif (splitExplicitStep == numTSIterations) then

#ifdef MPAS_OPENACC
            !$acc parallel present(minLevelCell, maxLevelCell) &
            !$acc   present(layerThicknessTend) &
            !$acc   present(layerThicknessCur) &
            !$acc   present(minLevelEdgeBot, maxLevelEdgeTop) &
            !$acc   present(layerThickEdgeFlux) &
            !$acc   present(activeTracerHorizontalAdvectionEdgeFlux) &
            !$acc   present(layerThicknessNew) &
            !$acc   present(activeTracerHorizontalAdvectionTendency) &
            !$acc   present(activeTracerVerticalAdvectionTendency) &
            !$acc   present(activeTracerHorMixTendency) &
            !$acc   present(activeTracerSurfaceFluxTendency) &
            !$acc   present(temperatureShortWaveTendency) &
            !$acc   present(activeTracerNonLocalTendency)
#endif


#ifdef MPAS_OPENACC
            !$acc loop gang
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! this is h_{n+1}
               layerThicknessNew(k,iCell) = &
                                     layerThicknessCur(k,iCell) + &
                                  dt*layerThicknessTend(k,iCell)
            end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

            if (config_compute_active_tracer_budgets) then

#ifdef MPAS_OPENACC
               !$acc loop gang
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k)
#endif
               do iEdge = 1, nEdgesAll
#ifdef MPAS_OPENACC
               !$acc loop seq
#endif
               do k= minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
#ifdef MPAS_OPENACC
               !$acc loop vector
#endif
                  do iTracer = 1, size(activeTracerHorizontalAdvectionEdgeFlux,1)
                     activeTracerHorizontalAdvectionEdgeFlux(iTracer,k,iEdge) = &
                     activeTracerHorizontalAdvectionEdgeFlux(iTracer,k,iEdge) / &
                          layerThickEdgeFlux(k,iEdge)
                  enddo
               enddo
               enddo
#ifndef MPAS_OPENACC
               !$omp end do
#endif

#ifdef MPAS_OPENACC
               !$acc loop gang
#else
               !$omp do schedule(runtime) private(k)
#endif
               do iCell = 1, nCellsAll
               do k= minLevelCell(iCell), maxLevelCell(iCell)
#ifdef MPAS_OPENACC
               !$acc loop vector
                  do iTracer = 1, size(activeTracerHorizontalAdvectionTendency,1)
                     activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) = &
                     activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) / &
                               layerThicknessNew(k,iCell)

                     activeTracerVerticalAdvectionTendency(iTracer,k,iCell) = &
                     activeTracerVerticalAdvectionTendency(iTracer,k,iCell) / &
                               layerThicknessNew(k,iCell)

                     activeTracerHorMixTendency(iTracer,k,iCell) = &
                     activeTracerHorMixTendency(iTracer,k,iCell) / &
                                layerThicknessNew(k,iCell)

                     activeTracerSurfaceFluxTendency(iTracer,k,iCell) = &
                     activeTracerSurfaceFluxTendency(iTracer,k,iCell) / &
                                layerThicknessNew(k,iCell)

                     activeTracerNonLocalTendency(iTracer,k,iCell) = &
                     activeTracerNonLocalTendency(iTracer,k,iCell) / &
                                layerThicknessNew(k,iCell)
                  end do
#else
                  activeTracerHorizontalAdvectionTendency(:,k,iCell) = &
                  activeTracerHorizontalAdvectionTendency(:,k,iCell) / &
                            layerThicknessNew(k,iCell)

                  activeTracerVerticalAdvectionTendency(:,k,iCell) = &
                  activeTracerVerticalAdvectionTendency(:,k,iCell) / &
                            layerThicknessNew(k,iCell)

                  activeTracerHorMixTendency(:,k,iCell) = &
                  activeTracerHorMixTendency(:,k,iCell) / &
                             layerThicknessNew(k,iCell)

                  activeTracerSurfaceFluxTendency(:,k,iCell) = &
                  activeTracerSurfaceFluxTendency(:,k,iCell) / &
                             layerThicknessNew(k,iCell)

                  activeTracerNonLocalTendency(:,k,iCell) = &
                  activeTracerNonLocalTendency(:,k,iCell) / &
                             layerThicknessNew(k,iCell)
#endif

                  temperatureShortWaveTendency(k,iCell) = &
                  temperatureShortWaveTendency(k,iCell) / &
                             layerThicknessNew(k,iCell)
               end do
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif
            endif

#ifdef MPAS_OPENACC
            !$acc end parallel
#endif

! This tracer block is still computed on the CPU
#ifdef MPAS_OPENACC
            !$acc update host(layerThicknessNew)
            !$acc update host(layerThicknessCur, layerThicknessTend)
            !$acc update host(activeTracersNew)
            !$acc update host(tracersSurfaceValue)
#endif

            call mpas_pool_begin_iteration(tracersPool)
            do while (mpas_pool_get_next_member(tracersPool, &
                                                groupItr) )
               if (groupItr%memberType == MPAS_POOL_FIELD) then
                  configName = 'config_use_'//trim(groupItr%memberName)
                  call mpas_pool_get_config(domain%configs, &
                                configName, config_use_tracerGroup)

                  if ( config_use_tracerGroup ) then
                     call mpas_pool_get_array(tracersPool, &
                                              groupItr%memberName, &
                                              tracersGroupCur, 1)
                     call mpas_pool_get_array(tracersPool, &
                                              groupItr%memberName, &
                                              tracersGroupNew, 2)

                     modifiedGroupName = &
                             trim(groupItr % memberName) // 'Tend'
                     call mpas_pool_get_array(tracersTendPool, &
                                              modifiedGroupName, &
                                              tracersGroupTend)

                     !$omp parallel
                     !$omp do schedule(runtime) private(k)
                     do iCell = 1, nCellsAll
                     do k = minLevelCell(iCell), maxLevelCell(iCell)
                        tracersGroupNew(:,k,iCell) = &
                       (tracersGroupCur(:,k,iCell) * &
                        layerThicknessCur(k,iCell) + dt* &
                       tracersGroupTend(:,k,iCell))/ &
                        layerThicknessNew(k,iCell)
                     end do
                     end do
#ifndef MPAS_OPENACC
                     !$omp end do
                     !$omp end parallel
#endif

                     ! limit salinity in separate loop
                     if (trim(groupItr%memberName) == &
                         'activeTracers' ) then
                        !$omp parallel
                        !$omp do schedule(runtime) private(k)
                        do iCell = 1, nCellsAll
                        do k = minLevelCell(iCell), maxLevelCell(iCell)
                           tracersGroupNew(indexSalinity,k,iCell) = &
                           max(0.001_RKIND,  &
                           tracersGroupNew(indexSalinity,k,iCell))
                        end do
                        end do
#ifndef MPAS_OPENACC
                        !$omp end do
                        !$omp end parallel
#endif
                     end if

                     ! Reset debugTracers to fixed value at the surface
                     if (trim(groupItr % memberName) == &
                         'debugTracers' .and. &
                         config_reset_debugTracers_near_surface) then

                        !$omp parallel
                        !$omp do schedule(runtime) private(k, lat)
                        do iCell = 1, nCellsAll

                           ! Reset tracer1 to 2 in top n layers
                           do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                tracersGroupNew(1,k,iCell) = 2.0_RKIND
                           end do

                           ! Reset tracer2 to 2 in top n layers
                           ! in zonal bands, and 1 outside
                           lat = latCell(iCell)*180.0_RKIND/pi
                           if (     lat>-60.0.and.lat<-55.0 &
                                .or.lat>-40.0.and.lat<-35.0 &
                                .or.lat>- 2.5.and.lat<  2.5 &
                                .or.lat> 35.0.and.lat< 40.0 &
                                .or.lat> 55.0.and.lat< 60.0 ) then
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(2,k,iCell) = 2.0_RKIND
                              end do
                           else
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(2,k,iCell) = 1.0_RKIND
                              end do
                           end if

                           ! Reset tracer3 to 2 in top n layers
                           ! in zonal bands, and 1 outside
                           lat = latCell(iCell)*180.0_RKIND/pi
                           if (     lat>-55.0.and.lat<-50.0 &
                                .or.lat>-35.0.and.lat<-30.0 &
                                .or.lat>-15.0.and.lat<-10.0 &
                                .or.lat> 10.0.and.lat< 15.0 &
                                .or.lat> 30.0.and.lat< 35.0 &
                                .or.lat> 50.0.and.lat< 55.0 ) then
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(3,k,iCell) = 2.0_RKIND
                              end do
                           else
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(3,k,iCell) = 1.0_RKIND
                              end do
                           end if
                        end do ! cells
#ifndef MPAS_OPENACC
                        !$omp end do
                        !$omp end parallel
#endif
                     end if ! debug tracers
                  end if ! use tracer group
               end if ! tracer
            end do ! tracer group

#ifdef MPAS_OPENACC
            !$acc update device(layerThicknessNew)
            !$acc update device(activeTracersNew)
#endif

            if (config_use_freq_filtered_thickness) then
#ifdef MPAS_OPENACC
               !$acc parallel present(minLevelCell, maxLevelCell) &
               !$acc   present(lowFreqDivergenceNew) &
               !$acc   present(lowFreqDivergenceCur) &
               !$acc   present(lowFreqDivergenceTend)
#endif

#ifdef MPAS_OPENACC
               !$acc loop gang
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k)
#endif
               do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
               !$acc loop vector
#endif
               do k = minLevelCell(iCell), maxLevelCell(iCell)

                  ! h^{hf}_{n+1} was computed in Stage 1

                  ! this is D^{lf}_{n+1}
                  lowFreqDivergenceNew(k,iCell) = &
                  lowFreqDivergenceCur(k,iCell) + dt* &
                  lowFreqDivergenceTend(k,iCell)
               end do
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif

#ifdef MPAS_OPENACC
               !$acc end parallel
#endif

            end if

            ! Recompute final u to go on to next step.
            ! u_{n+1} = normalBarotropicVelocity_{n+1} +
            !           normalBaroclinicVelocity_{n+1}
            ! Right now normalBaroclinicVelocityNew is at time n+1/2,
            ! so back compute to get normalBaroclinicVelocity at
            ! time n+1 using normalBaroclinicVelocity_{n+1/2} =
            !            1/2*(normalBaroclinicVelocity_n + u_Bcl_{n+1})
            ! so the following lines are
            ! u_{n+1} = normalBarotropicVelocity_{n+1} +
            !         2*normalBaroclinicVelocity_{n+1/2} -
            !         normalBaroclinicVelocity_n
            ! note that normalBaroclinicVelocity is recomputed at the
            ! beginning of the next timestep due to Imp Vert mixing,
            ! so normalBaroclinicVelocity does not have to be
            ! recomputed here.

#ifdef MPAS_OPENACC
            !$acc parallel present(minLevelEdgeBot, maxLevelEdgeTop) &
            !$acc   present(normalVelocityNew) &
            !$acc   present(normalBarotropicVelocityNew) &
            !$acc   present(normalBaroclinicVelocityNew) &
            !$acc   present(normalBaroclinicVelocityCur)
#endif

#ifdef MPAS_OPENACC
            !$acc loop gang
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iEdge = 1, nEdgesAll
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               normalVelocityNew(k,iEdge) = &
                               normalBarotropicVelocityNew(iEdge) + &
                             2*normalBaroclinicVelocityNew(k,iEdge) - &
                               normalBaroclinicVelocityCur(k,iEdge)
            end do
            end do ! iEdges
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

#ifdef MPAS_OPENACC
            !$acc end parallel
#endif

         endif ! splitExplicitStep

#ifdef MPAS_OPENACC
         !$acc exit data copyout(layerThicknessNew)
         !$acc exit data copyout(normalVelocityNew)
         !$acc exit data delete(normalBarotropicVelocityNew, &
         !$acc                  normalBaroclinicVelocityNew)
         !$acc exit data copyout(activeTracersNew)
         !$acc update host(layerThickEdgeFlux)
         if (config_use_freq_filtered_thickness) then
            !$acc exit data copyout(lowFreqDivergenceNew)
            !$acc exit data delete(lowFreqDivergenceCur)
            !$acc exit data delete(lowFreqDivergenceTend)
         endif
         if (splitExplicitStep < numTSIterations) then
            if (config_use_gm) then
               !$acc update host(vertGMBolusVelocityTop)
            end if
            if (config_submesoscale_enable) then
               !$acc update host(vertMLEBolusVelocityTop)
            end if
            !$acc exit data delete (atmosphericPressure, seaIcePressure)
            !$acc exit data copyout(sshNew)
            !$acc update host(layerThickEdgeMean)
            !$acc update host(relativeVorticity, circulation)
            !$acc update host(vertTransportVelocityTop, &
            !$acc             relativeVorticityCell, &
            !$acc             divergence, &
            !$acc             kineticEnergyCell, &
            !$acc             tangentialVelocity, &
            !$acc             vertVelocityTop)
            !$acc update host(normRelVortEdge, normPlanetVortEdge, &
            !$acc             normalizedRelativeVorticityCell)
            !$acc update host (surfacePressure)
            !$acc update host(zMid, zTop)
            !$acc update host(tracersSurfaceValue)
            !$acc update host(normalVelocitySurfaceLayer)
            !$acc update host(density, potentialDensity, displacedDensity)
            !$acc update host(thermExpCoeff,  &
            !$acc&            salineContractCoeff)
            !$acc update host(montgomeryPotential, pressure)
            if (landIcePressureOn) then
               !$acc exit data delete(landIcePressure)
               !$acc exit data delete(landIceDraft)
            endif
            if (config_use_freq_filtered_thickness) then
               !$acc exit data copyout(highFreqThicknessNew)
               !$acc exit data delete(highFreqThicknessCur)
            endif
            if ( associated(frazilSurfacePressure) ) then
               !$acc exit data delete(frazilSurfacePressure)
            endif
         elseif (splitExplicitStep == numTSIterations) then
            !$acc exit data delete(layerThicknessCur, layerThicknessTend)
            !$acc exit data delete(normalBaroclinicVelocityCur)
            if (config_compute_active_tracer_budgets) then
               !$acc update host(activeTracerHorizontalAdvectionEdgeFlux)
               !$acc update host(activeTracerHorizontalAdvectionTendency)
               !$acc update host(activeTracerVerticalAdvectionTendency)
               !$acc update host(activeTracerHorMixTendency)
               !$acc update host(activeTracerSurfaceFluxTendency)
               !$acc update host(temperatureShortWaveTendency)
               !$acc update host(activeTracerNonLocalTendency)
            endif
         endif
#endif

         call mpas_timer_stop('se loop fini')
         call mpas_timer_stop('se loop')

      end do  ! outer timestep loop (splitExplicitStep)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se implicit vert mix")

      ! Call ocean diagnostic solve in preparation for vertical mixing.
      ! Note it is called again after vertical mixing, because u and
      ! tracers change. For Richardson vertical mixing, only density,
      ! layerThickEdgeFlux, and kineticEnergyCell need to be computed.
      ! For kpp, more variables may be needed.  Either way, this could
      ! be made more efficient by only computing what is needed for the
      ! implicit vmix routine that follows.
#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
      !$acc update device (normalTransportVelocity)
      if (config_use_gm) then
         !$acc update device (normalGMBolusVelocity)
      end if
      if (config_submesoscale_enable) then
         !$acc update device (normalMLEVelocity)
      end if
      !$acc enter data copyin(atmosphericPressure, seaIcePressure)
      !$acc enter data copyin(sshNew)
      !$acc enter data copyin(activeTracersNew)
      !$acc update device(tracersSurfaceValue)
      if ( associated(frazilSurfacePressure) ) then
         !$acc enter data copyin(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc enter data copyin(landIcePressure)
         !$acc enter data copyin(landIceDraft)
      endif
#endif
      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, &
                                scratchPool, tracersPool, 2)

      call mpas_dmpar_field_halo_exch(domain, 'surfaceFrictionVelocity')
      call mpas_dmpar_field_halo_exch(domain, 'tangentialVelocity')
#ifdef MPAS_OPENACC
      !$acc update host(layerThickEdgeFlux, layerThickEdgeMean)
      !$acc update host(relativeVorticity, circulation)
      if (config_use_gm) then
         !$acc update host (vertGMBolusVelocityTop)
      end if
      if (config_submesoscale_enable) then
         !$acc update host (vertMLEBolusVelocityTop)
      end if
      !$acc update host(vertTransportVelocityTop, &
      !$acc             relativeVorticityCell, &
      !$acc             divergence, &
      !$acc             kineticEnergyCell, &
      !$acc             tangentialVelocity, &
      !$acc             vertVelocityTop)
      !$acc update host(normRelVortEdge, normPlanetVortEdge, &
      !$acc             normalizedRelativeVorticityCell)
      !$acc update host (surfacePressure)
      !$acc update host(zMid, zTop)
      !$acc exit data copyout(sshNew)
      !$acc exit data delete(activeTracersNew)
      !$acc update host(tracersSurfaceValue)
      !$acc update host(normalVelocitySurfaceLayer)
      !$acc exit data delete (atmosphericPressure, seaIcePressure)
      !$acc update host(density, potentialDensity, displacedDensity)
      !$acc update host(thermExpCoeff,  &
      !$acc&            salineContractCoeff)
      !$acc update host(montgomeryPotential, pressure)
      !$acc update host(RiTopOfCell, &
      !$acc             BruntVaisalaFreqTop)
      !$acc update host(tracersSurfaceLayerValue, &
      !$acc             indexSurfaceLayerDepth, &
      !$acc             normalVelocitySurfaceLayer, &
      !$acc             sfcFlxAttCoeff, &
      !$acc             surfaceFluxAttenuationCoefficientRunoff)
      if ( associated(frazilSurfacePressure) ) then
         !$acc exit data delete(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc exit data delete(landIcePressure)
         !$acc exit data delete(landIceDraft)
      endif
      !$acc exit data delete(layerThicknessNew, normalVelocityNew)
#endif

      call ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, &
                             scratchPool, err, 2)

      ! Update halo on u and tracers, which were just updated for
      ! implicit vertical mixing.  If not done, this leads to lack of
      ! volume conservation.  It is required because halo updates in
      ! stage 3 are only conducted on tendencies, not on the velocity
      ! and tracer fields.  So this update is required to communicate
      ! the change due to implicit vertical mixing across the boundary.

      call mpas_timer_start('se vmix halos')

      call mpas_timer_start('se vmix halos normalVelFld')
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', &
                                      timeLevel=2)
      call mpas_timer_stop('se vmix halos normalVelFld')

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if (groupItr%memberType == MPAS_POOL_FIELD) then

            ! Reset iAge to zero where mask == 0
            if (config_use_idealAgeTracers.and.trim(groupItr%memberName) == 'idealAgeTracers') then
                call mpas_pool_get_array(tracersPool, &
                                         groupItr%memberName, &
                                         tracersGroupNew, 2)
                call mpas_pool_get_subpool(forcingPool, &
                                           'tracersIdealAgeFields', &
                                            tracersIdealAgeFieldsPool)
                call mpas_pool_get_array(tracersIdealAgeFieldsPool, &
                                         'idealAgeTracersIdealAgeMask', &
                                          tracerGroupIdealAgeMask)

                !$omp parallel
                !$omp do schedule(runtime)
                do iCell = 1, nCellsOwned
                   tracersGroupNew(1,1,iCell) = tracerGroupIdealAgeMask(1,iCell)*tracersGroupNew(1,1,iCell)
                end do ! cells
                !$omp end do
                !$omp end parallel
             endif

            ! Halo update on all tracer groups
            call mpas_dmpar_field_halo_exch(domain, groupItr%memberName, timeLevel=2)

         end if
      end do

      call mpas_timer_stop('se vmix halos')

      call mpas_timer_stop("se implicit vert mix")

      if ( configVertAdvMethod == vertAdvRemap ) then

         call mpas_timer_start("vertical remap")

         call mpas_pool_get_array(statePool, 'layerThicknessLag', layerThicknessLagNew, 2)

         ! Store the current, Lagrangian location for computation of the
         ! Eulerian vertical velocity induced by remapping
         layerThicknessLagNew = layerThicknessNew

         ! Perform the remapping
         call ocn_remap_vert_state(block, err)

         call mpas_timer_stop("vertical remap")

      end if

      call mpas_timer_start('se fini')

#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
      if (config_prescribe_velocity) then
         !$acc enter data copyin(normalVelocityCur)
      endif
      if (config_prescribe_thickness) then
         !$acc enter data copyin(layerThicknessCur)
      endif
      !$acc update device (normalTransportVelocity)
      if (config_use_gm) then
         !$acc update device (normalGMBolusVelocity)
      end if
      if (config_submesoscale_enable) then
         !$acc update device (normalMLEvelocity)
      end if
      !$acc enter data copyin(atmosphericPressure, seaIcePressure)
      !$acc enter data copyin(sshNew)
      !$acc enter data copyin(activeTracersNew)
      !$acc update device(tracersSurfaceValue)
      if ( associated(frazilSurfacePressure) ) then
         !$acc enter data copyin(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc enter data copyin(landIcePressure)
         !$acc enter data copyin(landIceDraft)
      endif
#endif
      if (config_prescribe_velocity) then
#ifdef MPAS_OPENACC
         !$acc parallel loop collapse(2) present(normalVelocityNew, normalVelocityCur)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iEdge = 1, nEdgesAll
         do k=1,nVertLevels
            normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      if (config_prescribe_thickness) then
#ifdef MPAS_OPENACC
         !$acc parallel loop collapse(2) present(layerThicknessNew, layerThicknessCur)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCellsAll
         do k=1,nVertLevels
            layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      ! direct application of tidal boundary condition
      if (config_use_tidal_forcing .and. trim(config_tidal_forcing_type) == 'direct') then
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
        call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
        call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)
        do iCell=1, nCells
          ! artificially assumes boolean mask for now, could generalize to tappered sponge layer
          if (tidalInputMask(iCell) == 1.0_RKIND) then
            ! compute total depth for relative thickness contribution
            totalDepth = 0.0_RKIND
            do k = minLevelCell(iCell), maxLevelCell(iCell)
              totalDepth = totalDepth + restingThickness(k,iCell)
            end do

            ! only modify layer thicknesses on tidal boundary
            do k = minLevelCell(iCell), maxLevelCell(iCell)
              layerThicknessNew(k, iCell) = tidalInputMask(iCell)* &
                                            (tidalBCValue(iCell) + bottomDepth(iCell))* &
                                            (restingThickness(k,iCell)/totalDepth)
              !(1.0_RKIND - tidalInputMask(iCell))*layerThicknessNew(k, iCell)  ! generalized tappered assumption code
            end do
          end if
        end do
      end if

      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, &
                                scratchPool, tracersPool, 2)

      ! Update the effective desnity in land ice if we're coupling to
      ! land ice
      call ocn_effective_density_in_land_ice_update(forcingPool, &
                                                    statePool, err)

      call mpas_timer_start('se final mpas reconstruct', .false.)

#ifdef MPAS_OPENACC
      !$acc enter data create(velocityX, velocityY, velocityZ)
      !$acc enter data create(velocityZonal, velocityMeridional)
      !$acc enter data create(gradSSHX, gradSSHY, gradSSHZ)
      !$acc enter data create(gradSSHZonal, gradSSHMeridional)
      !$acc enter data create(surfaceVelocity, SSHGradient)
#endif

      call mpas_reconstruct_gpu(meshPool, normalVelocityNew,       &
                                velocityX, velocityY, velocityZ,   &
                                velocityZonal, velocityMeridional, &
                                includeHalos = .true.)

      call mpas_reconstruct_gpu(meshPool, gradSSH,               &
                                gradSSHX, gradSSHY, gradSSHZ,    &
                                gradSSHZonal, gradSSHMeridional, &
                                includeHalos = .true.)

      call mpas_timer_stop('se final mpas reconstruct')

#ifdef MPAS_OPENACC
      !$acc parallel loop present(surfaceVelocity, velocityZonal, velocityMeridional, &
      !$acc   SSHGradient, gradSSHZonal, gradSSHMeridional)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         surfaceVelocity(indexSurfaceVelocityZonal,iCell) = &
                                   velocityZonal(minLevelCell(iCell),iCell)
         surfaceVelocity(indexSurfaceVelocityMeridional,iCell) = &
                                   velocityMeridional(minLevelCell(iCell),iCell)

         SSHGradient(indexSSHGradientZonal,iCell) = gradSSHZonal(iCell)
         SSHGradient(indexSSHGradientMeridional,iCell) = &
                                               gradSSHMeridional(iCell)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc update host(layerThickEdgeFlux, layerThickEdgeMean)
      !$acc update host(relativeVorticity, circulation)
      if (config_use_gm) then
         !$acc update host(vertGMBolusVelocityTop)
      end if
      if (config_submesoscale_enable) then
         !$acc update host(vertMLEBolusVelocityTop)
      end if
      !$acc update host(vertTransportVelocityTop, &
      !$acc             relativeVorticityCell, &
      !$acc             divergence, &
      !$acc             kineticEnergyCell, &
      !$acc             tangentialVelocity, &
      !$acc             vertVelocityTop)
      !$acc update host(normRelVortEdge, normPlanetVortEdge, &
      !$acc             normalizedRelativeVorticityCell)
      !$acc update host (surfacePressure)
      !$acc update host(zMid, zTop)
      !$acc exit data copyout(sshNew)
      !$acc exit data delete(activeTracersNew)
      !$acc update host(tracersSurfaceValue)
      !$acc update host(normalVelocitySurfaceLayer)
      !$acc exit data delete (atmosphericPressure, seaIcePressure)
      if ( associated(frazilSurfacePressure) ) then
         !$acc exit data delete(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc exit data delete(landIcePressure)
         !$acc exit data delete(landIceDraft)
      endif
      !$acc exit data copyout(layerThicknessNew, normalVelocityNew)
      !$acc update host(density, potentialDensity, displacedDensity)
      !$acc update host(thermExpCoeff,  &
      !$acc&            salineContractCoeff)
      !$acc update host(montgomeryPotential, pressure)
      !$acc update host(RiTopOfCell, &
      !$acc             BruntVaisalaFreqTop)
      !$acc update host(tracersSurfaceLayerValue, &
      !$acc             indexSurfaceLayerDepth, &
      !$acc             normalVelocitySurfaceLayer, &
      !$acc             sfcFlxAttCoeff, &
      !$acc             surfaceFluxAttenuationCoefficientRunoff)
      if (config_prescribe_velocity) then
         !$acc exit data delete(normalVelocityCur)
      endif
      if (config_prescribe_thickness) then
         !$acc exit data delete(layerThicknessCur)
      endif
      !$acc exit data copyout(velocityX, velocityY, velocityZ)
      !$acc exit data copyout(gradSSHX, gradSSHY, gradSSHZ)
      !$acc exit data copyout(velocityZonal, velocityMeridional)
      !$acc exit data copyout(gradSSHZonal, gradSSHMeridional)
      !$acc exit data copyout(surfaceVelocity, SSHGradient)
#endif

      call ocn_time_average_coupled_accumulate(statePool,forcingPool,2)

      if (config_use_GM .or. config_submesoscale_enable) then
         call ocn_reconstruct_eddy_vectors(meshPool)
      end if

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("se effective density halo")
         call mpas_pool_get_field(statePool, &
                                 'effectiveDensityInLandIce', &
                                  effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("se effective density halo")
      end if
      deallocate(baroclinicThickness)

      call mpas_timer_stop('se fini')
      call mpas_timer_stop("se timestep")

   end subroutine ocn_time_integrator_split!}}}

!***********************************************************************
!
!  routine ocn_time_integration_split_init
!
!> \brief   Initialize split-explicit time stepping within ocean
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes variables required for the split-explicit
!>  method of integrating the ocean model forward in time.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integration_split_init(domain)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain      !< [inout] data structure containing most variables

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: &
         statePool,      &! structure containing model state
         meshPool,       &! structure containing mesh fields
         tracersPool      ! structure containing tracer fields

      integer ::         &
         iCell, iEdge, k,&! loop indices for cell, edge and vertical
         kmax,           &! index of deepest active edge
         ierr,           &! local error flag
         cell1, cell2     ! neighbor cell indices across edge

      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) ::       &
         normalThicknessFluxSum, &! vertical sum of thick flux
         layerThicknessSum,      &! vertical sum of layer thickness
         layerThicknessEdge1      ! layer thickness on edge

      real (kind=RKIND), dimension(:), pointer :: &
         refBottomDepth,         &! reference bottom depth
         normalBarotropicVelocity ! normal barotropic velocity

      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThickness,          &! layer thickness cell center
         normalBaroclinicVelocity,&! normal baroclinic velocity
         normalVelocity            ! normal velocity (total)

      type (mpas_time_type) :: &
         nowTime          ! current model time

      type (mpas_timeInterval_type) :: &! various time step calculations
         fullTimeStep,       &! full model time step
         barotropicTimeStep, &! barotropic time step
         remainder,          &! remaining time after interval division
         zeroInterval         ! zero timestep for comparing remainder

      integer, pointer :: nVertLevels, nCells, nEdges

      integer (kind=I8KIND) :: nBtrSubcyclesI8

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Determine the number of barotropic subcycles based on the
      !*** ratio of time steps

      nowTime = mpas_get_clock_time(domain%clock, MPAS_NOW, ierr)

      call mpas_set_timeInterval( zeroInterval, S=0 )
      call mpas_set_timeInterval( fullTimeStep , timeString=config_dt )
      call mpas_set_timeInterval( barotropicTimeStep, &
                                             timeString=config_btr_dt )

      !*** transfer to I8 for division step
      nBtrSubcyclesI8 = nBtrSubcycles
      call mpas_interval_division( nowTime, fullTimeStep, &
                       barotropicTimeStep, nBtrSubcyclesI8, remainder )
      nBtrSubcycles = nBtrSubcyclesI8

      if ( remainder > zeroInterval ) then
         nBtrSubcycles = nBtrSubcycles + 1
      end if

      !*** Set mask for using velocity correction
      if (config_vel_correction) then
         useVelocityCorrection = 1.0_RKIND
      else
         useVelocityCorrection = 0.0_RKIND
      endif

      !*** Compute some halo needs to reduce the number
      !*** of halo updates during subcycling
      neededHalos = 1 + config_n_btr_cor_iter
      haloDecrement = mod(config_num_halos, neededHalos)

      !*** Determine the time integration type and set associated masks

      select case (trim(config_time_integrator))

      case ('unsplit_explicit')
         unsplit   = .true.
         splitFact = 0.0_RKIND

      case ('split_explicit')
         unsplit   = .false.
         splitFact = 1.0_RKIND
         call mpas_log_write( '*******************************************************************************')
         call mpas_log_write( 'The split explicit time integration is configured to use: $i barotropic subcycles', &
            intArgs=(/ nBtrSubcycles /) )
         call mpas_log_write( '*******************************************************************************')

      case default
         call mpas_log_write('Incorrect choice config_time_integrator',&
                             MPAS_LOG_CRIT)

      end select

      if (trim(config_thickness_flux_type) == 'upwind') then
         thickEdgeFluxChoice = thickEdgeFluxUpwind
      else
         thickEdgeFluxChoice = thickEdgeFluxCenter
      end if
      !*** set number of baroclinic iterations on each outer
      !*** time step iteration (number can be different on the
      !*** first and last time step iteration)

      numTSIterations     = config_n_ts_iter
      allocate(numClinicIterations(numTSIterations))

      numClinicIterations    = config_n_bcl_iter_mid ! most iterations
      numClinicIterations(1) = config_n_bcl_iter_beg ! first iteration
      numClinicIterations(numTSIterations)=config_n_bcl_iter_end !last

      !*** Initialize some variables if not a restart

      if (.not. config_do_restart) then

         !*** Initialize z-level mesh variables from h, read in from
         !*** input file.

         block => domain % blocklist
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

         call mpas_pool_get_array(statePool, 'layerThickness', &
                                              layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', &
                                              normalVelocity, 1)
         call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                              normalBarotropicVelocity, 1)
         call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                              normalBaroclinicVelocity, 1)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

         !*** Compute barotropic velocity at first timestep
         !*** This is only done upon start-up.
         if (unsplit) then

            do iEdge = 1, nEdgesAll
               normalBarotropicVelocity(  iEdge) = 0.0_RKIND
               normalBaroclinicVelocity(:,iEdge) = &
                         normalVelocity(:,iEdge)
            end do

         else ! split explicit

            if (config_filter_btr_mode) then
               do iCell = 1, nCells
                  layerThickness(minLevelCell(iCell),iCell) = refBottomDepth(minLevelCell(iCell))
               enddo
            endif

            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               kmax  = maxLevelEdgeTop(iEdge)

               ! normalBarotropicVelocity = sum(h*u)/sum(h) on each edge
               ! ocn_diagnostic_solve has not yet been called, so
               ! compute hEdge just for this edge.

               ! thicknessSum is initialized outside the loop because on
               ! land boundaries maxLevelEdgeTop=0, but we want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               layerThicknessEdge1 = 0.5_RKIND* &
                                     (layerThickness(minLevelCell(cell1),cell1) + &
                                      layerThickness(minLevelCell(cell2),cell2) )
               normalThicknessFluxSum = layerThicknessEdge1* &
                                        normalVelocity(minLevelEdgeBot(iEdge),iEdge)
               layerThicknessSum = layerThicknessEdge1

               do k=minLevelEdgeBot(iEdge)+1, kmax
                  layerThicknessEdge1 = 0.5_RKIND* &
                                       (layerThickness(k,cell1) + &
                                        layerThickness(k,cell2))

                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                           layerThicknessEdge1* &
                                           normalVelocity(k,iEdge)
                  layerThicknessSum = layerThicknessSum + &
                                      layerThicknessEdge1

               enddo
               normalBarotropicVelocity(iEdge) = &
                     normalThicknessFluxSum/layerThicknessSum

               ! normalBaroclinicVelocity = normalVelocity -
               !     normalBarotropicVelocity
               do k = minLevelEdgeBot(iEdge), kmax
                  normalBaroclinicVelocity(k,iEdge) = &
                            normalVelocity(k,iEdge) - &
                    normalBarotropicVelocity(iEdge)
               enddo

               ! normalBaroclinicVelocity=0,
               ! normalVelocity=0 on land cells
               do k = kmax+1, nVertLevels
                  normalBaroclinicVelocity(k,iEdge) = 0.0_RKIND
                  normalVelocity(k,iEdge) = 0.0_RKIND
               enddo
            enddo ! edge loop

            if (config_filter_btr_mode) then
               ! filter normalBarotropicVelocity out of initial condition

               normalVelocity(:,:) = normalBaroclinicVelocity(:,:)
               normalBarotropicVelocity(:) = 0.0_RKIND

            endif

         endif ! split explicit

      end if ! not restart

      ! Create reusable halo exchange for subcycling
      call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
      call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                           'sshSubcycle')
      call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                  'normalBarotropicVelocitySubcycle')
      call mpas_dmpar_exch_group_build_reusable_buffers(domain, &
                                                     subcycleGroupName)

      !-----------------------------------------------------------------
      if (config_use_self_attraction_loading) then
         ! set pgf_sal_on to 1
         pgf_sal_on = 1.0_RKIND
      else
         pgf_sal_on = 0.0_RKIND
      endif

      self_attraction_and_loading_beta = config_self_attraction_and_loading_beta

   end subroutine ocn_time_integration_split_init!}}}

!***********************************************************************

end module ocn_time_integration_split

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
