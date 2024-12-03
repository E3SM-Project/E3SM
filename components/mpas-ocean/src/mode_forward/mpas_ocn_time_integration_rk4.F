! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_rk4
!
!> \brief MPAS ocean RK4 Time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the RK4 time integration routine.
!
!-----------------------------------------------------------------------

module ocn_time_integration_rk4

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_gm
   use ocn_submesoscale_eddies

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average_coupled
   use ocn_wetting_drying

   use ocn_effective_density_in_land_ice
   use ocn_surface_land_ice_fluxes
   use ocn_transport_tests
   use ocn_time_varying_forcing

   use ocn_subgrid


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

   public :: ocn_time_integrator_rk4
   public :: ocn_time_integration_rk4_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_rk4
!
!> \brief MPAS ocean RK4 Time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates one timestep (dt) using an RK4 time integrator.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_rk4(domain, dt)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:))
   !                 plus mesh meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      real (kind=RKIND), intent(in) :: dt !< Input: timestep

      integer :: iCell, iEdge, iTracer, k, i, err, num_tracers
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: tracersTendPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: provisStatePool
      type (mpas_pool_type), pointer :: provisTracersPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: swForcingPool

      type (dm_info) :: dminfo

      integer :: rk_step

      type (mpas_pool_type), pointer :: nextProvisPool, prevProvisPool

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights, forcingTimeIncrementRK4

      real (kind=RKIND) :: coef
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell

      ! Dimensions
      integer, pointer :: nCells, nEdges, nVertLevels

      ! Config options
      logical, pointer :: config_prescribe_velocity, config_prescribe_thickness
      logical, pointer :: config_filter_btr_mode, config_use_freq_filtered_thickness
      logical, pointer :: config_use_GM
      logical, pointer :: config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup
      logical, pointer :: config_disable_thick_all_tend
      logical, pointer :: config_disable_vel_all_tend
      logical, pointer :: config_disable_tr_all_tend
      real (kind=RKIND), pointer :: config_mom_del4
      real (kind=RKIND), pointer :: config_drying_min_cell_height
      logical, pointer :: config_use_wetting_drying
      logical, pointer :: config_verify_not_dry
      logical, pointer :: config_prevent_drying
      logical, pointer :: config_zero_drying_velocity
      character (len=StrKIND), pointer :: config_land_ice_flux_mode

      ! Diagnostics Indices

      ! Mesh array pointers
      integer, dimension(:), pointer :: minLevelCell, maxLevelCell, maxLevelEdgeTop
      real (kind=RKIND), dimension(:), pointer :: bottomDepth

      ! Provis Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityProvis, layerThicknessProvis
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessProvis
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceProvis
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupProvis

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, &
                                                    layerThicknessTend
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend

      ! State Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew
      real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew

      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroup, tracersCur, tracersNew

      ! Diagnostics Field Pointers
      type (field1DReal), pointer :: boundaryLayerDepthField, effectiveDensityField
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField

      ! State/Tend Field Pointers
      type (field2DReal), pointer :: normalVelocityField, layerThicknessField
      type (field3DReal), pointer :: tracersGroupField

      ! Tracer Group Iteartion
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      ! Tidal boundary condition
      logical, pointer :: config_use_tidal_forcing
      character (len=StrKIND), pointer :: config_tidal_forcing_type
      real (kind=RKIND), pointer :: forcingTimeIncrement
      real (kind=RKIND), dimension(:), pointer :: tidalInputMask, tidalBCValue
      real (kind=RKIND), dimension(:,:), pointer :: restingThickness
      real (kind=RKIND) :: totalDepth

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersCur, activeTracersNew

      real(kind=RKIND), dimension(:,:), pointer :: subgridSshCellTableRange, subgridWetVolumeCellTable


      ! Get config options
      call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_GM', config_use_GM)
      call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      call mpas_pool_get_config(domain % configs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_use_wetting_drying', config_use_wetting_drying)
      call mpas_pool_get_config(domain % configs, 'config_prevent_drying', config_prevent_drying)
      call mpas_pool_get_config(domain % configs, 'config_verify_not_dry', config_verify_not_dry)
      call mpas_pool_get_config(domain % configs, 'config_drying_min_cell_height', config_drying_min_cell_height)
      call mpas_pool_get_config(domain % configs, 'config_zero_drying_velocity', config_zero_drying_velocity)
      call mpas_pool_get_config(domain % configs, 'config_use_tidal_forcing', config_use_tidal_forcing)
      call mpas_pool_get_config(domain % configs, 'config_tidal_forcing_type', config_tidal_forcing_type)

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with layerThickness in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

         call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         call mpas_pool_copy_pool(statePool, provisStatePool, 1)

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

         call mpas_pool_get_array(meshPool, 'subgridWetVolumeCellTable', &
                                subgridWetVolumeCellTable)

         call mpas_pool_get_array(meshPool, 'subgridSshCellTableRange', &
                                subgridSshCellTableRange)

         call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)

         forcingTimeIncrement = 0.0_RKIND

         ! Lower k-loop limit of 1 rather than minLevel* needed in *New = *Cur
         ! assignments below are needed to maintain bit-for-bit results

         !$omp parallel
         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdges
            do k = 1, maxLevelEdgeTop(iEdge)
               normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
            end do
         end do
         !$omp end do

         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
            end do
         end do
         !$omp end do
         !$omp end parallel

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )

            if ( groupItr % memberType == MPAS_POOL_FIELD ) then

               call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersCur, 1)
               call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersNew, 2)

               if ( associated(tracersCur) .and. associated(tracersNew) ) then
                  !$omp parallel
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells  ! couple tracers to thickness
                     do k = 1, maxLevelCell(iCell)
                        tracersNew(:, k, iCell) = tracersCur(:, k, iCell) * layerThicknessCur(k, iCell)
                     end do
                  end do
                  !$omp end do
                  !$omp end parallel
               end if
            end if
         end do

         if (associated(highFreqThicknessCur)) then
              !$omp parallel
              !$omp do schedule(runtime)
              do iCell = 1, nCells
                 highFreqThicknessNew(:, iCell) = highFreqThicknessCur(:, iCell)
              end do
              !$omp end do
              !$omp end parallel
         end if

         if (associated(lowFreqDivergenceCur)) then
              !$omp parallel
              !$omp do schedule(runtime)
              do iCell = 1, nCells
                 lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceCur(:, iCell)
              end do
              !$omp end do
              !$omp end parallel
         end if

         block => block % next
      end do

      ! Fourth-order Runge-Kutta, solving dy/dt = f(t,y) is typically written as follows
      ! where h = delta t is the large time step.  Here f(t,y) is the right hand side,
      ! called the tendencies in the code below.
      ! k_1 = h f(t_n        , y_n)
      ! k_2 = h f(t_n + 1/2 h, y_n + 1/2 k_1)
      ! k_3 = h f(t_n + 1/2 h, y_n + 1/2 k_2)
      ! k_4 = h f(t_n +     h, y_n +     k_3)
      ! y_{n+1} = y_n + 1/6 k_1 + 1/3 k_2 + 1/3 k_3 + 1/6 k_4

      ! in index notation:
      ! k_{j+1} = h f(t_n + a_j h, y_n + a_j k_j)
      ! y_{n+1} = y_n + sum ( b_j k_j )

      ! The coefficients of k_j are b_j = (1/6, 1/3, 1/3, 1/6) and are
      ! initialized here as delta t * b_j:

      rk_weights(1) = dt/6.0_RKIND
      rk_weights(2) = dt/3.0_RKIND
      rk_weights(3) = dt/3.0_RKIND
      rk_weights(4) = dt/6.0_RKIND

      ! The a_j coefficients of h in the computation of k_j are typically written (0, 1/2, 1/2, 1).
      ! However, in the algorithm below we pre-compute the state for the tendency one iteration early.
      ! That is, on j=1 (rk_step=1, below) we pre-compute y_n + 1/2 k_1 and save it in provis_state.
      ! Then we compute 1/6 k_1 and add it to state % time_levs(2).
      ! That is why the coefficients of h are one index early in the following, i.e.
      ! a = (1/2, 1/2, 1)

      rk_substep_weights(1) = dt/2.0_RKIND
      rk_substep_weights(2) = dt/2.0_RKIND
      rk_substep_weights(3) = dt
      rk_substep_weights(4) = dt ! a_4 only used for ALE step, otherwise it is skipped.

      ! these are time increments to evaluate the tidal forcing at the
      ! intermediate time-steps as required by RK4

      forcingTimeIncrementRK4(1) = 0.0_RKIND
      forcingTimeIncrementRK4(2) = dt/2.0_RKIND
      forcingTimeIncrementRK4(3) = dt/2.0_RKIND
      forcingTimeIncrementRK4(4) = dt

      call mpas_timer_start("RK4-main loop")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN RK loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do rk_step = 1, 4

        if (config_disable_thick_all_tend .and. config_disable_vel_all_tend .and. config_disable_tr_all_tend) then
          exit ! don't compute in loop meant to update velocity, thickness, and tracers
        end if

        if (config_transport_tests_flow_id > 0) then
          ! This is a transport test. Write advection velocity from prescribed
          ! flow field.
          call ocn_transport_test_velocity(meshPool, daysSinceStartOfSim, &
            rk_substep_weights(rk_step), normalVelocityNew)
        endif


        ! Update halos for diagnostic variables.
        if (config_use_cvmix_kpp) then
           call mpas_timer_start("RK4-boundary layer depth halo update")
           call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
           call mpas_timer_stop("RK4-boundary layer depth halo update")
        end if


        call mpas_timer_start("RK4-diagnostic halo update")


        call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
        if (config_mom_del4 > 0.0_RKIND) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
        end if
        call mpas_timer_stop("RK4-diagnostic halo update")

        ! Compute tendencies for high frequency thickness
        ! In RK4 notation, we are computing the right hand side f(t,y),
        ! which is the same as k_j / h.

        if (config_use_freq_filtered_thickness) then
           call mpas_timer_start("RK4-tendency computations")

           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call ocn_tend_freq_filtered_thickness(tendPool, provisStatePool, 1)
              block => block % next
           end do

           call mpas_timer_stop("RK4-tendency computations")

           call mpas_timer_start("RK4-prognostic halo update")

           call mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
           call mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

           call mpas_timer_stop("RK4-prognostic halo update")


           ! Compute next substep state for high frequency thickness.
           ! In RK4 notation, we are computing y_n + a_j k_j.
           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
              call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)
              call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)

              call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

              !$omp parallel
              !$omp do schedule(runtime)
              do iCell = 1, nCells
                 highFreqThicknessProvis(:, iCell) = highFreqThicknessCur(:, iCell) + rk_substep_weights(rk_step) &
                    * highFreqThicknessTend(:, iCell)
              end do
              !$omp end do
              !$omp end parallel
              block => block % next
           end do

        endif

        ! require that cells don't dry out
        if (config_use_wetting_drying) then
             call mpas_timer_start("RK4-prevent drying")

             ! compute wetting velocity to prevent drying of cell (sets up start of next iterate to not dry)
             if (config_prevent_drying) then
               block => domain % blocklist
                 do while (associated(block))
                   call ocn_prevent_drying_rk4(domain, block, dt, rk_substep_weights(rk_step), config_zero_drying_velocity, err)
                   block => block % next
                 end do
               ! exchange fields for parallelization
               call mpas_pool_get_field(statePool, 'normalVelocity', normalVelocityField, 1)
               call mpas_dmpar_exch_halo_field(normalVelocityField)
               call mpas_dmpar_exch_halo_field(wettingVelocityField)
             end if

             call mpas_timer_stop("RK4-prevent drying")
        end if

        ! Compute tendencies for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing the right hand side f(t,y),
        ! which is the same as k_j / h.
        call mpas_timer_start("RK4 vel/thick tendency computations")

        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
           call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
           forcingTimeIncrement = forcingTimeIncrementRK4(rk_step)
           block => block % next
        end do

        call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_rk4_compute_vel_tends(domain, block, dt, rk_substep_weights(rk_step), domain % dminfo, err )
           call ocn_time_integrator_rk4_compute_thick_tends( block, dt, rk_substep_weights(rk_step), err )
           block => block % next
        end do

        call mpas_timer_stop("RK4 vel/thick tendency computations")

        ! Update halos for prognostic variables.

        call mpas_timer_start("RK4 vel/thick prognostic halo update")

        call mpas_dmpar_field_halo_exch(domain, 'tendNormalVelocity')
        call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

        call mpas_timer_stop("RK4 vel/thick prognostic halo update")

        call mpas_timer_start("RK4 tracer tendency computations")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_rk4_compute_tracer_tends( block, dt, rk_substep_weights(rk_step), err )
           block => block % next
        end do

        call mpas_timer_stop("RK4 tracer tendency computations")

        call mpas_timer_start("RK4 tracer prognostic halo update")

        call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

        call mpas_pool_begin_iteration(tracersTendPool)
        do while ( mpas_pool_get_next_member(tracersTendPool, groupItr) )
           if ( groupItr % memberType == MPAS_POOL_FIELD ) then
              call mpas_dmpar_field_halo_exch(domain, trim(groupItr % memberName))
           end if
        end do

        call mpas_timer_stop("RK4 tracer prognostic halo update")

        ! Compute next substep state for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing y_n + a_j k_j.

        call mpas_timer_start("RK4-update diagnostic variables")

        if (rk_step < 4) then
           block => domain % blocklist
           do while (associated(block))
              call ocn_time_integrator_rk4_diagnostic_update(block, dt, rk_substep_weights(rk_step), err)
              block => block % next
           end do
        end if

        call mpas_timer_stop("RK4-update diagnostic variables")

        ! Accumulate update.
        ! In RK4 notation, we are computing b_j k_j and adding it to an accumulating sum so that we have
        !    y_{n+1} = y_n + sum ( b_j k_j )
        ! after the fourth iteration.

        call mpas_timer_start("RK4-accumulate update")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_rk4_accumulate_update(block, rk_weights(rk_step), err)

           block => block % next
        end do

        call mpas_timer_stop("RK4-accumulate update")

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END RK loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! verify that cells are not dry at conclusion of time step
      if (config_use_wetting_drying) then
        call mpas_timer_start("RK4- check wet cells")

        ! ensure existing layerThickness is valid
        if (config_verify_not_dry) then
          block => domain % blocklist
            do while (associated(block))
              call ocn_wetting_drying_verify(block, config_drying_min_cell_height, err)
              block => block % next
            end do
        end if

        call mpas_timer_stop("RK4- check wet cells")
      end if

      call mpas_timer_stop("RK4-main loop")

      !
      !  A little clean up at the end: rescale tracer fields and compute diagnostics for new state
      !
      call mpas_timer_start("RK4-cleanup phase")

      ! Rescale tracers
      block => domain % blocklist
      do while(associated(block))
        call ocn_time_integrator_rk4_cleanup(domain, block, dt, err)

        block => block % next
      end do

      call mpas_timer_start("RK4-implicit vert mix")
      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done,
      ! this leads to lack of volume conservation.  It is required because halo updates in RK4 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("RK4-implicit vert mix halos")

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         end if
      end do

      call mpas_timer_stop("RK4-implicit vert mix halos")

      call mpas_timer_stop("RK4-implicit vert mix")

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
         call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
         call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)

         call mpas_pool_get_array(tracersPool, 'activeTracers',activeTracersCur, 1)
         call mpas_pool_get_array(tracersPool, 'activeTracers',activeTracersNew, 2)

         if (config_prescribe_velocity) then
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge)
            end do
            !$omp end do
            !$omp end parallel
         end if

         if (config_disable_tr_all_tend) then
            num_tracers = size(activeTracersNew, dim=1)
            do iTracer = 1, num_tracers
               !$omp parallel
               !$omp do schedule(runtime)
               do iCell = 1, nCells
                  activeTracersNew(iTracer, :, iCell) = activeTracersCur(iTracer, :, iCell)
               end do
               !$omp end do
               !$omp end parallel
            end do
         end if

         if (config_prescribe_thickness) then
            !$omp parallel
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               layerThicknessNew(:, iCell) = layerThicknessCur(:, iCell)
            end do
            !$omp end do
            !$omp end parallel
         end if

         ! direct application of tidal boundary condition
         if (config_use_tidal_forcing .and. trim(config_tidal_forcing_type) == 'direct') then
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
                 layerThicknessNew(k, iCell) = tidalInputMask(iCell)*(tidalBCValue(iCell) + bottomDepth(iCell))*(restingThickness(k,iCell)/totalDepth)
                 !(1.0_RKIND - tidalInputMask(iCell))*layerThicknessNew(k, iCell)  ! generalized tappered assumption code
               end do

               if ( config_use_subgrid_wetting_drying ) then
                  call ocn_subgrid_layer_thickness_lookup( tidalBCValue(iCell), & 
                                         subgridWetVolumeCellTable(:,iCell), &
                                         subgridSshCellTableRange(:,iCell),&
                                         bottomDepth(iCell), &
                                         layerThicknessNew(1,iCell) )
               end if

             end if
           end do
         end if

         call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
         call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
         call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

         if (landIcePressureOn) then
            call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
            call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
         endif

#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
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
         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, scratchPool, tracersPool, 2)

         ! Update the effective density in land ice if we're coupling to land ice
         call ocn_effective_density_in_land_ice_update(forcingPool, &
                                                       statePool, err)

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
         !$acc exit data delete(layerThicknessNew, normalVelocityNew)
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
#endif

         ! ------------------------------------------------------------------
         ! Accumulating various parameterizations of the transport velocity
         ! ------------------------------------------------------------------
         !$omp parallel
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
         end do
         !$omp end do
         !$omp end parallel

         ! add submesoscale aend gm components if requested
         call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdges, nVertLevels)
         call ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

         ! ------------------------------------------------------------------
         ! End: Accumulating various parameterizations of the transport velocity
         ! ------------------------------------------------------------------

         call mpas_reconstruct(meshPool,  normalVelocityNew, &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

         call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)

         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(minLevelCell(iCell), iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(minLevelCell(iCell), iCell)

            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         end do
         !$omp end do
         !$omp end parallel

         call ocn_time_average_coupled_accumulate(statePool, forcingPool, 2)

         if (config_use_GM .or. config_submesoscale_enable) then
            call ocn_reconstruct_eddy_vectors(meshPool)
         end if

         block => block % next
      end do

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("RK4-effective density halo")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("RK4-effective density halo")
      end if

      call mpas_timer_stop("RK4-cleanup phase")

   end subroutine ocn_time_integrator_rk4!}}}

   subroutine ocn_time_integrator_rk4_compute_vel_tends(domain, block, dt, &
     rkSubstepWeight, dminfo, err)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkSubstepWeight
      integer, intent(out) :: err
      type (dm_info), intent(in) :: dminfo

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur,layerThickEdgeFlux, normalVelocityProvis, &
            sshCur, rkSubstepWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur,layerThickEdgeFlux, normalVelocityProvis, &
            sshCur, rkSubstepWeight, &
            vertAleTransportTop, err)
      endif

      call ocn_tend_vel(domain, tendPool, provisStatePool, forcingPool, 1, dminfo, dt)

   end subroutine ocn_time_integrator_rk4_compute_vel_tends!}}}

   subroutine ocn_time_integrator_rk4_compute_thick_tends(block, dt, rkSubstepWeight, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkSubstepWeight
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      integer :: iEdge
      integer, pointer :: nEdges, nVertLevels
      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

      if (config_use_GM .or. config_submesoscale_enable) then
         !$omp parallel
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityProvis(:, iEdge)
         end do
         !$omp end do
         !$omp end parallel
      end if

      ! add submesoscale and GM contributions if requested
      call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdges, nVertLevels)
      call ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur, layerThickEdgeFlux, normalTransportVelocity, &
            sshCur, rkSubstepWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur, layerThickEdgeFlux, normalTransportVelocity, &
            sshCur, rkSubstepWeight, &
            vertAleTransportTop, err)
      endif

      call ocn_tend_thick(tendPool, forcingPool)


   end subroutine ocn_time_integrator_rk4_compute_thick_tends!}}}

   subroutine ocn_time_integrator_rk4_compute_tracer_tends(block, dt, rkSubstepWeight, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkSubstepWeight
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: swForcingPool, tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      logical, pointer :: config_filter_btr_mode

      integer :: iEdge
      integer, pointer :: nEdges, nVertLevels
      err = 0

      if (config_disable_tr_all_tend) return

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

      if (config_use_GM .or. config_submesoscale_enable) then
         !$omp parallel
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityProvis(:, iEdge)
         end do
         !$omp end do
         !$omp end parallel
      end if

      call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdges, nVertLevels)
      call ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur, layerThickEdgeFlux, normalTransportVelocity, &
            sshCur, rkSubstepWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur, layerThickEdgeFlux, normalTransportVelocity, &
            sshCur, rkSubstepWeight, &
            vertAleTransportTop, err)
      endif

      if (config_filter_btr_mode) then
          call ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, meshPool, 1)
      endif

      call ocn_tend_tracer(tendPool, provisStatePool, forcingPool, meshPool, swForcingPool, &
                           dt, activeTracersOnlyIn=.false., timeLevelIn=1)

   end subroutine ocn_time_integrator_rk4_compute_tracer_tends!}}}

   subroutine ocn_time_integrator_rk4_compute_tends(domain, block, dt, rkWeight, dminfo, err)!{{{
      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkWeight
      integer, intent(out) :: err

      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool, provisStatePool
      type (mpas_pool_type), pointer :: swForcingPool, tracersPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityProvis, highFreqThicknessProvis

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur,layerThickEdgeFlux, normalVelocityProvis, &
            sshCur, rkWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur,layerThickEdgeFlux, normalVelocityProvis, &
            sshCur, rkWeight, &
            vertAleTransportTop, err)
      endif

      call ocn_tend_vel(domain, tendPool, provisStatePool, forcingPool, 1, dminfo, dt)

      if (associated(highFreqThicknessProvis)) then
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur, layerThickEdgeFlux, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err, highFreqThicknessProvis)
      else
         call ocn_vert_transport_velocity_top(verticalMeshPool, &
            layerThicknessCur, layerThickEdgeFlux, normalTransportVelocity, &
            sshCur, rkWeight, &
            vertAleTransportTop, err)
      endif

      call ocn_tend_thick(tendPool, forcingPool)

      if (config_filter_btr_mode) then
          call ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, meshPool, 1)
      endif

      call ocn_tend_tracer(tendPool, provisStatePool, forcingPool, meshPool, swForcingPool, &
                           dt, activeTracersOnlyIn=.false., timeLevelIn=1)

   end subroutine ocn_time_integrator_rk4_compute_tends!}}}

   subroutine ocn_time_integrator_rk4_diagnostic_update(block, dt, rkSubstepWeight, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkSubstepWeight
      integer, intent(out) :: err

      logical, pointer :: config_prescribe_velocity, config_prescribe_thickness, config_use_GM

      integer, pointer :: nCells, nEdges, nVertLevels
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, tendPool, meshPool, verticalMeshPool, scratchPool
      type (mpas_pool_type), pointer :: provisStatePool, forcingPool
      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool, provisTracersPool

      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityCur, normalVelocityProvis, normalVelocityTend
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, layerThicknessProvis, layerThicknessTend
      real (kind=RKIND), dimension(:, :), pointer :: lowFreqDivergenceCur, lowFreqDivergenceProvis, lowFreqDivergenceTend

      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupCur, tracersGroupProvis, tracersGroupTend

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell, maxLevelEdgeTop
      integer, dimension(:,:), pointer :: edgeMask

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft

      real (kind=RKIND), dimension(:), pointer :: sshCur

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersCur

      err = 0

      call mpas_pool_get_config(block % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(block % configs, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_get_config(block % configs, 'config_use_GM', config_use_GM)

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

      call mpas_pool_get_subpool(provisStatePool, 'tracers', provisTracersPool)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

      call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'layerThickness', layerThicknessProvis, 1)
      call mpas_pool_get_array(provisStatePool, 'lowFreqDivergence', lowFreqDivergenceProvis, 1)

      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersCur, 1)

      ! the loop below is on all vertlevels but we used edgemask to set
      ! velocities to zero below topo (esp. edges between topo and ocean)
      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            normalVelocityProvis(k, iEdge) = edgeMask(k, iEdge)*(normalVelocityCur(k, iEdge) + rkSubstepWeight &
                                           * normalVelocityTend(k, iEdge))
            normalVelocityProvis(k, iEdge) = normalVelocityProvis(k, iEdge) * (1.0_RKIND - wettingVelocityFactor(k, iEdge))
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = 1, maxLevelCell(iCell)
            layerThicknessProvis(k, iCell) = layerThicknessCur(k, iCell) + rkSubstepWeight &
                                           * layerThicknessTend(k, iCell)
         end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            configName = 'config_use_' // trim(groupItr % memberName)
            call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

            if ( config_use_tracerGroup ) then
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
               call mpas_pool_get_array(provisTracersPool, groupItr % memberName, tracersGroupProvis, 1)

               modifiedGroupName = trim(groupItr % memberName) // 'Tend'
               call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
               if ( associated(tracersGroupProvis) .and. associated(tracersGroupCur) .and. associated(tracersGroupTend) ) then
                  !$omp parallel
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupProvis(:, k, iCell) = ( layerThicknessCur(k, iCell) * tracersGroupCur(:, k, iCell)  &
                                                 + rkSubstepWeight * tracersGroupTend(:, k, iCell) &
                                                   ) / layerThicknessProvis(k, iCell)
                     end do

                  end do
                  !$omp end do
                  !$omp end parallel
               end if
            end if
         end if
      end do

      if (associated(lowFreqDivergenceCur)) then
         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            lowFreqDivergenceProvis(:, iCell) = lowFreqDivergenceCur(:, iCell) + rkSubstepWeight &
                                              * lowFreqDivergenceTend(:, iCell)
         end do
         !$omp end do
         !$omp end parallel
      end if

      if (config_prescribe_velocity) then
         !$omp parallel
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalVelocityProvis(:, iEdge) = normalVelocityCur(:, iEdge)
         end do
         !$omp end do
         !$omp end parallel
      end if

      if (config_prescribe_thickness) then
         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            layerThicknessProvis(:, iCell) = layerThicknessCur(:, iCell)
         end do
         !$omp end do
         !$omp end parallel
      end if

      call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

      if (landIcePressureOn) then
         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
      endif

#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThicknessCur, normalVelocityCur)
      !$acc update device (normalTransportVelocity)
      if (config_use_gm) then
         !$acc update device (normalGMBolusVelocity)
      end if
      if (config_submesoscale_enable) then
         !$acc update device (normalMLEvelocity)
      end if
      !$acc enter data copyin(atmosphericPressure, seaIcePressure)
      !$acc enter data copyin(sshCur)
      !$acc enter data copyin(activeTracersCur)
      !$acc update device(tracersSurfaceValue)
      if ( associated(frazilSurfacePressure) ) then
         !$acc enter data copyin(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc enter data copyin(landIcePressure)
         !$acc enter data copyin(landIceDraft)
      endif
#endif
      call ocn_diagnostic_solve(dt, provisStatePool, forcingPool, meshPool, verticalMeshPool, scratchPool, tracersPool, 1)
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
      !$acc exit data copyout(sshCur)
      !$acc exit data delete(activeTracersCur)
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
      !$acc exit data delete(layerThicknessCur, normalVelocityCur)
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
#endif

      ! ------------------------------------------------------------------
      ! Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------
      !$omp parallel
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         normalTransportVelocity(:, iEdge) = normalVelocityProvis(:, iEdge)
      end do
      !$omp end do
      !$omp end parallel

      ! add submesoscale and GM contributions if requested
      call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdges, nVertLevels)
      call ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

      ! ------------------------------------------------------------------
      ! End: Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------

   end subroutine ocn_time_integrator_rk4_diagnostic_update!}}}

   subroutine ocn_time_integrator_rk4_accumulate_update(block, rkWeight, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: rkWeight
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges, nVertLevels
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, tendPool, meshPool
      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool

      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityNew, normalVelocityTend
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessNew, layerThicknessTend
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThicknessNew, highFreqThicknessTend
      real (kind=RKIND), dimension(:, :), pointer :: lowFreqDivergenceNew, lowFreqDivergenceTend

      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupNew, tracersGroupTend

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell
      integer, dimension(:,:), pointer :: edgeMask

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      err = 0

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

      !call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      !call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      !call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
      !call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

      call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = minLevelCell(iCell), maxLevelCell(iCell)
            layerThicknessNew(k, iCell) = layerThicknessNew(k, iCell) + rkWeight * layerThicknessTend(k, iCell)
         end do
      end do
      !$omp end do
      
      ! the loop below is on all vertlevels but we used edgemask to set
      ! velocities to zero below topo (esp. edges between topo and ocean)
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            normalVelocityNew(k, iEdge) = edgeMask(k,iEdge) &
                                      *(normalVelocityNew(k, iEdge) + rkWeight * normalVelocityTend(k, iEdge))
            normalVelocityNew(k, iEdge) = normalVelocityNew(k, iEdge) * (1.0_RKIND - wettingVelocityFactor(k, iEdge))
         end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            configName = 'config_use_' // trim(groupItr % memberName)
            call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

            if ( config_use_tracerGroup ) then
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

               modifiedGroupName = trim(groupItr % memberName) // 'Tend'
               call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
               if ( associated(tracersGroupNew) .and. associated(tracersGroupTend) ) then
                  !$omp parallel
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = minLevelCell(iCell), maxLevelCell(iCell)
                        tracersGroupNew(:, k, iCell) = tracersGroupNew(:, k, iCell) + rkWeight &
                                                * tracersGroupTend(:, k, iCell)
                     end do
                  end do
                  !$omp end do
                  !$omp end parallel
               end if
            end if
         end if
      end do

      if (associated(highFreqThicknessNew)) then
         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            highFreqThicknessNew(:, iCell) = highFreqThicknessNew(:, iCell) + rkWeight * highFreqThicknessTend(:, iCell)
         end do
         !$omp end do
         !$omp end parallel
      end if

      if (associated(lowFreqDivergenceNew)) then
         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceNew(:, iCell) + rkWeight * lowFreqDivergenceTend(:, iCell)
         end do
         !$omp end do
         !$omp end parallel
      end if

   end subroutine ocn_time_integrator_rk4_accumulate_update!}}}

   subroutine ocn_time_integrator_rk4_cleanup(domain, block, dt, err)!{{{
      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, meshPool, verticalMeshPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessNew, normalVelocityNew
      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupNew

      integer, dimension(:), pointer :: minLevelCell, maxLevelCell

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      logical, pointer :: config_use_GM

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft

      real (kind=RKIND), dimension(:), pointer :: sshNew

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersNew

      err = 0

      call mpas_pool_get_config(block % configs, 'config_use_GM', config_use_GM)

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersNew, 2)

      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)

      if (landIcePressureOn) then
          call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
          call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
      endif

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
            if ( associated(tracersGroupNew) ) then
               !$omp parallel
               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                 do k = minLevelCell(iCell), maxLevelCell(iCell)
                   tracersGroupNew(:, k, iCell) = tracersGroupNew(:, k, iCell) / layerThicknessNew(k, iCell)
                 end do
               end do
               !$omp end do
               !$omp end parallel
            end if
         end if
      end do

#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
      !$acc update device (normalTransportVelocity, &
      !$acc                normalGMBolusVelocity,normalMLEvelocity)
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
      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, scratchPool, tracersPool, 2)
      call mpas_dmpar_field_halo_exch(domain, 'tangentialVelocity')

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
      !$acc exit data delete (layerThicknessNew, normalVelocityNew)
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
#endif

      call ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, scratchPool, err, 2)

   end subroutine ocn_time_integrator_rk4_cleanup!}}}

!***********************************************************************
!
!  routine ocn_time_integration_rk4_init
!
!> \brief   Initialize RK4 time stepping within ocean
!> \author  Carolyn Begeman
!> \date    April 2023
!> \details
!>  This routine checks config options for RK4 integrator
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integration_rk4_init(domain)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: &
         domain      !< [inout] data structure containing most variables
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool, statePool, provisStatePool
      type (mpas_pool_type), pointer :: nextProvisPool, prevProvisPool
      logical, pointer :: config_use_debugTracers
      integer, pointer :: nVertLevels
      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      call mpas_pool_get_config(domain % configs, 'config_use_debugTracers', config_use_debugTracers)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         if (config_use_debugTracers .and. nVertLevels == 1) then
            call mpas_log_write('Debug tracers may cause failures in a ' &
               // 'single layer case. Consider setting ' &
               // 'config_use_debugTracers to .false.', MPAS_LOG_WARN)
         endif

         call mpas_pool_create_pool(provisStatePool)
         call mpas_pool_clone_pool(statePool, provisStatePool, 1)
         call mpas_pool_add_subpool(block % structs, 'provis_state', provisStatePool)

         if (associated(block % prev)) then
            call mpas_pool_get_subpool(block % prev % structs, 'provis_state', prevProvisPool)
         else
            nullify(prevProvisPool)
         end if

         if (associated(block % next)) then
            call mpas_pool_get_subpool(block % next % structs, 'provis_state', nextProvisPool)
         else
            nullify(nextProvisPool)
         end if

         call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         if (associated(prevProvisPool) .and. associated(nextProvisPool)) then
            call mpas_pool_link_pools(provisStatePool, prevProvisPool, nextProvisPool)
         else if (associated(prevProvisPool)) then
            call mpas_pool_link_pools(provisStatePool, prevProvisPool)
         else if (associated(nextProvisPool)) then
            call mpas_pool_link_pools(provisStatePool, nextPool=nextProvisPool)
         else
            call mpas_pool_link_pools(provisStatePool)
         end if

         call mpas_pool_link_parinfo(block, provisStatePool)
         block => block % next
      end do

   end subroutine ocn_time_integration_rk4_init

end module ocn_time_integration_rk4

! vim: foldmethod=marker
