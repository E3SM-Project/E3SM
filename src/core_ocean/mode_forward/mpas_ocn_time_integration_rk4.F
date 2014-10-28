! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   use mpas_grid_types
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average
   use ocn_time_average_coupled
   use ocn_sea_ice

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

      integer :: iCell, k, i, err
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: provisStatePool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: averagePool

      integer :: rk_step

      type (mpas_pool_type), pointer :: nextProvisPool, prevProvisPool

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      real (kind=RKIND) :: coef
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:), allocatable:: A,C,uTemp
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      ! Dimensions
      integer, pointer :: nCells, nEdges, nVertLevels, num_tracers

      ! Config options
      logical, pointer :: config_prescribe_velocity, config_prescribe_thickness
      logical, pointer :: config_filter_btr_mode, config_use_freq_filtered_thickness
      logical, pointer :: config_use_standardGM
      logical, pointer :: config_use_cvmix_kpp
      real (kind=RKIND), pointer :: config_mom_del4

      ! State indices
      integer, pointer :: indexTemperature
      integer, pointer :: indexSalinity

      ! Diagnostics Indices
      integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

      ! Mesh array pointers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop

      ! Provis Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityProvis, layerThicknessProvis
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessProvis
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceProvis
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersProvis

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, layerThicknessTend
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersTend

      ! Diagnostics Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
      real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop
      real (kind=RKIND), dimension(:,:), pointer :: normalTransportVelocity, normalGMBolusVelocity
      real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:,:), pointer :: gradSSH
      real (kind=RKIND), dimension(:,:), pointer :: gradSSHX, gradSSHY, gradSSHZ
      real (kind=RKIND), dimension(:,:), pointer :: gradSSHZonal, gradSSHMeridional
      real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, sshGradient

      ! State Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew
      real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew

      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, tracersCur, tracersNew

      ! Forcing Array pointers
      real (kind=RKIND), dimension(:), pointer :: seaIceEnergy

      ! Diagnostics Field Pointers
      type (field1DReal), pointer :: boundaryLayerDepthField
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField

      ! State/Tend Field Pointers
      type (field2DReal), pointer :: highFreqThicknessField, lowFreqDivergenceField
      type (field2DReal), pointer :: normalVelocityField, layerThicknessField
      type (field3DReal), pointer :: tracersField

      ! Get config options
      call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_standardGM', config_use_standardGM)
      call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with layerThickness in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         
         allocate(provisStatePool)
         call mpas_pool_create_pool(provisStatePool)

         call mpas_pool_clone_pool(statePool, provisStatePool, 1)
         call mpas_pool_add_subpool(block % structs, 'provis_state', provisStatePool)

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
         call mpas_pool_get_array(statePool, 'tracers', tracersCur, 1)
         call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         normalVelocityNew(:,:) = normalVelocityCur(:,:)
         layerThicknessNew(:,:) = layerThicknessCur(:,:)

         do iCell = 1, nCells  ! couple tracers to thickness
            do k = 1, maxLevelCell(iCell)
               tracersNew(:,k,iCell) = tracersCur(:,k,iCell) * layerThicknessCur(k,iCell)
            end do
         end do

         if (associated(highFreqThicknessCur)) then
              highFreqThicknessNew(:,:) = highFreqThicknessCur(:,:)
         end if

         if (associated(lowFreqDivergenceCur)) then
              lowFreqDivergenceNew(:,:) = lowFreqDivergenceCur(:,:)
         end if

         block => block % next
      end do

      block => domain % blocklist
      do while(associated(block))
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

      rk_weights(1) = dt/6.
      rk_weights(2) = dt/3.
      rk_weights(3) = dt/3.
      rk_weights(4) = dt/6.

      ! The a_j coefficients of h in the computation of k_j are typically written (0, 1/2, 1/2, 1).
      ! However, in the algorithm below we pre-compute the state for the tendency one iteration early.
      ! That is, on j=1 (rk_step=1, below) we pre-compute y_n + 1/2 k_1 and save it in provis_state.
      ! Then we compute 1/6 k_1 and add it to state % time_levs(2).
      ! That is why the coefficients of h are one index early in the following, i.e.
      ! a = (1/2, 1/2, 1)

      rk_substep_weights(1) = dt/2.
      rk_substep_weights(2) = dt/2.
      rk_substep_weights(3) = dt
      rk_substep_weights(4) = dt ! a_4 only used for ALE step, otherwise it is skipped.

      call mpas_timer_start("RK4-main loop")
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do rk_step = 1, 4
        call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

        ! Update halos for diagnostic variables.
        call mpas_timer_start("RK4-boundary layer depth halo update")
        if (config_use_cvmix_kpp) then
           call mpas_pool_get_field(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepthField)
           call mpas_dmpar_exch_halo_field(boundaryLayerDepthField)
        end if
        call mpas_timer_stop("RK4-boundary layer depth halo update")

        call mpas_timer_start("RK4-diagnostic halo update")

        call mpas_pool_get_field(diagnosticsPool, 'normalizedRelativeVorticityEdge', normalizedRelativeVorticityEdgeField)
        call mpas_pool_get_field(diagnosticsPool, 'divergence', divergenceField)
        call mpas_pool_get_field(diagnosticsPool, 'relativeVorticity', relativeVorticityField)

        call mpas_dmpar_exch_halo_field(normalizedRelativeVorticityEdgeField)
        if (config_mom_del4 > 0.0) then
           call mpas_dmpar_exch_halo_field(divergenceField)
           call mpas_dmpar_exch_halo_field(relativeVorticityField)
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
              call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call ocn_tend_freq_filtered_thickness(tendPool, provisStatePool, diagnosticsPool, meshPool, 1)
              block => block % next
           end do
           call mpas_timer_stop("RK4-tendency computations")

           call mpas_timer_start("RK4-prognostic halo update")
           call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)

           call mpas_pool_get_field(tendPool, 'highFreqThickness', highFreqThicknessField)
           call mpas_pool_get_field(tendPool, 'lowFreqDivergence', lowFreqDivergenceField)

           call mpas_dmpar_exch_halo_field(highFreqThicknessField)
           call mpas_dmpar_exch_halo_field(lowFreqDivergenceField)
           call mpas_timer_stop("RK4-prognostic halo update")

           ! Compute next substep state for high frequency thickness.
           ! In RK4 notation, we are computing y_n + a_j k_j.

           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
              call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)
              call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)

              highFreqThicknessProvis(:,:) = highFreqThicknessCur(:,:) + rk_substep_weights(rk_step) * highFreqThicknessTend(:,:)
              block => block % next
           end do

        endif

        ! Compute tendencies for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing the right hand side f(t,y), 
        ! which is the same as k_j / h.
        call mpas_timer_start("RK4-tendency computations")
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
           call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
           call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

           call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
           call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
           call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

           call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
           call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
           call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

           call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
           call mpas_pool_get_array(provisStatePool, 'highFreqThickness', highFreqThicknessProvis, 1)

           ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
           call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, &
              layerThicknessCur,layerThicknessEdge, normalVelocityProvis, &
              sshCur, highFreqThicknessProvis, rk_substep_weights(rk_step), &
              vertAleTransportTop, err)

           call ocn_tend_vel(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1)

           call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, &
              layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
              sshCur, highFreqThicknessProvis, rk_substep_weights(rk_step), &
              vertAleTransportTop, err)

           call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

           if (config_filter_btr_mode) then
               call ocn_filter_btr_mode_tend_vel(tendPool, provisStatePool, diagnosticsPool, meshPool, 1)
           endif

           call ocn_tend_tracer(tendPool, provisStatePool, forcingPool, diagnosticsPool, meshPool, scratchPool, dt, 1)
           block => block % next
        end do
        call mpas_timer_stop("RK4-tendency computations")

        ! Update halos for prognostic variables.

        call mpas_timer_start("RK4-prognostic halo update")
        call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)

        call mpas_pool_get_field(tendPool, 'normalVelocity', normalVelocityField)
        call mpas_pool_get_field(tendPool, 'layerThickness', layerThicknessField)
        call mpas_pool_get_field(tendPool, 'tracers', tracersField)

        call mpas_dmpar_exch_halo_field(normalVelocityField)
        call mpas_dmpar_exch_halo_field(layerThicknessField)
        call mpas_dmpar_exch_halo_field(tracersField)
        call mpas_timer_stop("RK4-prognostic halo update")

        ! Compute next substep state for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing y_n + a_j k_j.

        call mpas_timer_start("RK4-update diagnostic variables")
        if (rk_step < 4) then
           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
              call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
              call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

              call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
              call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
              call mpas_pool_get_array(statePool, 'tracers', tracersCur, 1)
              call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

              call mpas_pool_get_array(provisStatePool, 'normalVelocity', normalVelocityProvis, 1)
              call mpas_pool_get_array(provisStatePool, 'layerThickness', layerThicknessProvis, 1)
              call mpas_pool_get_array(provisStatePool, 'tracers', tracersProvis, 1)
              call mpas_pool_get_array(provisStatePool, 'lowFreqDivergence', lowFreqDivergenceProvis, 1)

              call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
              call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
              call mpas_pool_get_array(tendPool, 'tracers', tracersTend)
              call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

              call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

              call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
              call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

              normalVelocityProvis(:,:) = normalVelocityCur(:,:) + rk_substep_weights(rk_step) * normalVelocityTend(:,:)

              layerThicknessProvis(:,:) = layerThicknessCur(:,:) + rk_substep_weights(rk_step) * layerThicknessTend(:,:)
              do iCell = 1, nCells
                 do k = 1, maxLevelCell(iCell)
                 tracersProvis(:,k,iCell) = ( layerThicknessCur(k,iCell) * tracersCur(:,k,iCell)  &
                                          + rk_substep_weights(rk_step) * tracersTend(:,k,iCell) &
                                            ) / layerThicknessProvis(k,iCell)
                 end do

              end do

              if (associated(lowFreqDivergenceCur)) then
                 lowFreqDivergenceProvis(:,:) = lowFreqDivergenceCur(:,:) + rk_substep_weights(rk_step) * lowFreqDivergenceTend(:,:)
              end if

              if (config_prescribe_velocity) then
                 normalVelocityProvis(:,:) = normalVelocityCur(:,:)
              end if

              if (config_prescribe_thickness) then
                 layerThicknessProvis(:,:) = layerThicknessCur(:,:)
              end if

              call ocn_diagnostic_solve(dt, provisStatePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 1)

              ! ------------------------------------------------------------------
              ! Accumulating various parametrizations of the transport velocity
              ! ------------------------------------------------------------------
              normalTransportVelocity(:,:) = normalVelocityProvis(:,:)

              ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
              if (config_use_standardGM) then
                 call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
              end if

              if (config_use_standardGM) then
                 normalTransportVelocity(:,:) = normalTransportVelocity(:,:) + normalGMBolusVelocity(:,:)
              end if
              ! ------------------------------------------------------------------
              ! End: Accumulating various parametrizations of the transport velocity
              ! ------------------------------------------------------------------

              block => block % next
           end do
        end if
        call mpas_timer_stop("RK4-update diagnostic variables")

        ! Accumulate update.
        ! In RK4 notation, we are computing b_j k_j and adding it to an accumulating sum so that we have
        !    y_{n+1} = y_n + sum ( b_j k_j ) 
        ! after the fourth iteration.

        call mpas_timer_start("RK4-RK4 accumulate update")
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

           call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
           call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersCur, 1)
           call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
           call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

           call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
           call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)
           call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
           call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

           call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
           call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
           call mpas_pool_get_array(tendPool, 'tracers', tracersTend)
           call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
           call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

           call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

           normalVelocityNew(:,:) = normalVelocityNew(:,:) + rk_weights(rk_step) * normalVelocityTend(:,:) 

           layerThicknessNew(:,:) = layerThicknessNew(:,:) + rk_weights(rk_step) * layerThicknessTend(:,:) 

           do iCell = 1, nCells
              do k = 1, maxLevelCell(iCell)
                 tracersNew(:,k,iCell) =  tracersNew(:,k,iCell) + rk_weights(rk_step) * tracersTend(:,k,iCell)
              end do
           end do

           if (associated(highFreqThicknessNew)) then
              highFreqThicknessNew(:,:) = highFreqThicknessNew(:,:) + rk_weights(rk_step) * highFreqThicknessTend(:,:) 
           end if

           if (associated(lowFreqDivergenceNew)) then
              lowFreqDivergenceNew(:,:) = lowFreqDivergenceNew(:,:) + rk_weights(rk_step) * lowFreqDivergenceTend(:,:) 
           end if

           block => block % next
        end do
        call mpas_timer_stop("RK4-RK4 accumulate update")

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_timer_stop("RK4-main loop")

      !
      !  A little clean up at the end: rescale tracer fields and compute diagnostics for new state
      !
      call mpas_timer_start("RK4-cleaup phase")

      ! Rescale tracers
      block => domain % blocklist
      do while(associated(block))
        call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

        call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)
        call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)

        do iCell = 1, nCells
          do k = 1, maxLevelCell(iCell)
            tracersNew(:, k, iCell) = tracersNew(:, k, iCell) / layerThicknessNew(k, iCell)
          end do
        end do

        call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 2)
        call ocn_sea_ice_formation(meshPool, indexTemperature, indexSalinity, layerThicknessNew, tracersNew, seaIceEnergy, err)
        block => block % next
      end do

      call mpas_timer_start("RK4-implicit vert mix")
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        ! Call ocean diagnostic solve in preparation for vertical mixing.  Note 
        ! it is called again after vertical mixing, because u and tracers change.
        ! For Richardson vertical mixing, only density, layerThicknessEdge, and kineticEnergyCell need to 
        ! be computed.  For kpp, more variables may be needed.  Either way, this
        ! could be made more efficient by only computing what is needed for the
        ! implicit vmix routine that follows. 
        call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 2)

        call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, err, 2)

        ! ------------------------------------------------------------------
        ! Accumulating various parametrizations of the transport velocity
        ! ------------------------------------------------------------------
        normalTransportVelocity(:,:) = normalVelocityNew(:,:)

        ! Compute normalGMBolusVelocity, slopeRelative and RediDiffVertCoef if respective flags are turned on
        ! QC Note: this routine is called here to get updated k33. normalTransportVelocity probably does not need to be updated at all here.
        if (config_use_standardGM) then
           call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
        end if

        if (config_use_standardGM) then
           normalTransportVelocity(:,:) = normalTransportVelocity(:,:) + normalGMBolusVelocity(:,:)
        end if
        ! ------------------------------------------------------------------
        ! End: Accumulating various parametrizations of the transport velocity
        ! ------------------------------------------------------------------

        block => block % next
      end do

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done, 
      ! this leads to lack of volume conservation.  It is required because halo updates in RK4 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to 
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("RK4-implicit vert mix halos")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

      call mpas_pool_get_field(statePool, 'normalVelocity', normalVelocityField, 2)
      call mpas_pool_get_field(statePool, 'tracers', tracersField, 2)

      call mpas_dmpar_exch_halo_field(normalVelocityField)
      call mpas_dmpar_exch_halo_field(tracersField)
      call mpas_timer_stop("RK4-implicit vert mix halos")

      call mpas_timer_stop("RK4-implicit vert mix")

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'average', averagePool)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityZonal', indexSurfaceVelocityZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityMeridional', indexSurfaceVelocityMeridional)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientZonal', indexSSHGradientZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientMeridional', indexSSHGradientMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
         call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)

         if (config_prescribe_velocity) then
            normalVelocityNew(:,:) = normalVelocityCur(:,:)
         end if

         if (config_prescribe_thickness) then
            layerThicknessNew(:,:) = layerThicknessCur(:,:)
         end if

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 2)

         ! ------------------------------------------------------------------
         ! Accumulating various parameterizations of the transport velocity
         ! ------------------------------------------------------------------
         normalTransportVelocity(:,:) = normalVelocityNew(:,:)

         ! Compute normalGMBolusVelocity and the tracer transport velocity
         if (config_use_standardGM) then
             call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
         end if

         if (config_use_standardGM) then
            normalTransportVelocity(:,:) = normalTransportVelocity(:,:) + normalGMBolusVelocity(:,:)
         end if
         ! ------------------------------------------------------------------
         ! End: Accumulating various parameterizations of the transport velocity
         ! ------------------------------------------------------------------

         call mpas_reconstruct(meshPool,  normalVelocityNew,                 &
                          velocityX, velocityY, velocityZ, &
                          velocityZonal, velocityMeridional      &
                         )

         call mpas_reconstruct(meshPool, gradSSH,         &
                          gradSSHX, gradSSHY, gradSSHZ,   &
                          gradSSHZonal, gradSSHMeridional &
                         )

         surfaceVelocity(indexSurfaceVelocityZonal, :) = velocityZonal(1, :)
         surfaceVelocity(indexSurfaceVelocityMeridional, :) = velocityMeridional(1, :)

         SSHGradient(indexSSHGradientZonal, :) = gradSSHZonal(1, :)
         SSHGradient(indexSSHGradientMeridional, :) = gradSSHMeridional(1, :)

         call ocn_time_average_accumulate(averagePool, statePool, diagnosticsPool, 2)
         call ocn_time_average_coupled_accumulate(diagnosticsPool, forcingPool)

         if (config_use_standardGM) then
            call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         end if

         block => block % next
      end do
      call mpas_timer_stop("RK4-cleaup phase")

      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         call mpas_pool_destroy_pool(provisStatePool)

         call mpas_pool_remove_subpool(block % structs, 'provis_state')
         block => block % next
      end do

   end subroutine ocn_time_integrator_rk4!}}}

end module ocn_time_integration_rk4

! vim: foldmethod=marker
