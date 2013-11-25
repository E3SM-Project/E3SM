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
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_tendency
   use ocn_diagnostics

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

      integer :: rk_step

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      integer :: nCells, nEdges, nVertLevels, num_tracers
      real (kind=RKIND) :: coef
      real (kind=RKIND), dimension(:,:), pointer :: &
        u, layerThickness, layerThicknessEdge, vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer, dimension(:), pointer :: & 
        maxLevelCell, maxLevelEdgeTop
      real (kind=RKIND), dimension(:), allocatable:: A,C,uTemp
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      call mpas_setup_provis_state(domain % blocklist)

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with layerThickness in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      do while (associated(block))

         block % state % time_levs(2) % state % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)
         block % state % time_levs(2) % state % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)

         do iCell=1,block % mesh % nCells  ! couple tracers to thickness
            do k=1,block % mesh % maxLevelCell % array(iCell)
               block % state % time_levs(2) % state % tracers % array(:,k,iCell) = block % state % time_levs(1) % state % tracers % array(:,k,iCell) &
                  * block % state % time_levs(1) % state % layerThickness % array(k,iCell)
            end do
         end do

         if (config_use_freq_filtered_thickness) then
              block % state % time_levs(2) % state % highFreqThickness % array(:,:) = block % state % time_levs(1) % state % highFreqThickness % array(:,:)
              block % state % time_levs(2) % state % lowFreqDivergence % array(:,:) = block % state % time_levs(1) % state % lowFreqDivergence % array(:,:)
         endif

         call mpas_copy_state(block % provis_state, block % state % time_levs(1) % state)

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

        ! Update halos for diagnostic variables.
        call mpas_timer_start("RK4-boundary layer depth halo update")
        if (config_use_cvmix_kpp) then
           call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % boundaryLayerDepth)
        end if
        call mpas_timer_stop("RK4-boundary layer depth halo update")

        call mpas_timer_start("RK4-diagnostic halo update")
        call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % normalizedRelativeVorticityEdge)
        if (config_mom_del4 > 0.0) then
           call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % divergence)
           call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % relativeVorticity)
        end if
        call mpas_timer_stop("RK4-diagnostic halo update")

        ! Compute tendencies for high frequency thickness
        ! In RK4 notation, we are computing the right hand side f(t,y), 
        ! which is the same as k_j / h.

        if (config_use_freq_filtered_thickness) then
           call mpas_timer_start("RK4-tendency computations")
           block => domain % blocklist
           do while (associated(block))
              call ocn_tend_freq_filtered_thickness(block % tend, block % provis_state, block % diagnostics, block % mesh)
              block => block % next
           end do
           call mpas_timer_stop("RK4-tendency computations")

           call mpas_timer_start("RK4-prognostic halo update")
           call mpas_dmpar_exch_halo_field(domain % blocklist % tend % highFreqThickness)
           call mpas_dmpar_exch_halo_field(domain % blocklist % tend % lowFreqDivergence)
           call mpas_timer_stop("RK4-prognostic halo update")

           ! Compute next substep state for high frequency thickness.
           ! In RK4 notation, we are computing y_n + a_j k_j.

           block => domain % blocklist
           do while (associated(block))
                 block % provis_state % highFreqThickness % array(:,:) = block % state % time_levs(1) % state % highFreqThickness % array(:,:)  &
                    + rk_substep_weights(rk_step) * block % tend % highFreqThickness % array(:,:)
              block => block % next
           end do

        endif

        ! Compute tendencies for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing the right hand side f(t,y), 
        ! which is the same as k_j / h.
        call mpas_timer_start("RK4-tendency computations")
        block => domain % blocklist
        do while (associated(block))

           ! advection of u uses u, while advection of layerThickness and tracers use uTransport.
           call ocn_vert_transport_velocity_top(block % mesh, block % verticalMesh, &
              block % state % time_levs(1) % state % layerThickness % array, &
              block % diagnostics % layerThicknessEdge % array, block % provis_state % normalVelocity % array, &
              block % state % time_levs(1) % state % ssh % array, block % provis_state % highFreqThickness % array, &
              rk_substep_weights(rk_step), block % diagnostics % vertTransportVelocityTop % array, err)

           call ocn_tend_vel(block % tend, block % provis_state, block % forcing, block % diagnostics, block % mesh, block % scratch)

           call ocn_vert_transport_velocity_top(block % mesh, block % verticalMesh, &
              block % state % time_levs(1) % state % layerThickness % array, &
              block % diagnostics % layerThicknessEdge % array, block % diagnostics % uTransport % array, &
              block % state % time_levs(1) % state % ssh % array, block % provis_state % highFreqThickness % array, &
              rk_substep_weights(rk_step), block % diagnostics % vertTransportVelocityTop % array, err)

           call ocn_tend_thick(block % tend, block % provis_state, block % forcing, block % diagnostics, block % mesh)

           if (config_filter_btr_mode) then
               call ocn_filter_btr_mode_tend_vel(block % tend, block % provis_state, block % diagnostics, block % mesh)
           endif

           call ocn_tend_tracer(block % tend, block % provis_state, block % forcing, block % diagnostics, block % mesh, dt)
           block => block % next
        end do
        call mpas_timer_stop("RK4-tendency computations")

        ! Update halos for prognostic variables.

        call mpas_timer_start("RK4-prognostic halo update")
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % normalVelocity)
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % layerThickness)
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % tracers)
        call mpas_timer_stop("RK4-prognostic halo update")

        ! Compute next substep state for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing y_n + a_j k_j.

        call mpas_timer_start("RK4-update diagnostic variables")
        if (rk_step < 4) then
           block => domain % blocklist
           do while (associated(block))

              block % provis_state % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)  &
                                                          + rk_substep_weights(rk_step) * block % tend % normalVelocity % array(:,:)

              block % provis_state % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)  &
                                                    + rk_substep_weights(rk_step) * block % tend % layerThickness % array(:,:)
              do iCell=1,block % mesh % nCells
                 do k=1,block % mesh % maxLevelCell % array(iCell)
                 block % provis_state % tracers % array(:,k,iCell) = ( block % state % time_levs(1) % state % layerThickness % array(k,iCell) * &
                                                                 block % state % time_levs(1) % state % tracers % array(:,k,iCell)  &
                                                             + rk_substep_weights(rk_step) * block % tend % tracers % array(:,k,iCell) &
                                                               ) / block % provis_state % layerThickness % array(k,iCell)
                 end do

              end do

              if (config_use_freq_filtered_thickness) then

                 block % provis_state % lowFreqDivergence % array(:,:) = block % state % time_levs(1) % state % lowFreqDivergence % array(:,:)  &
                    + rk_substep_weights(rk_step) * block % tend % lowFreqDivergence % array(:,:)

              end if

              if (config_prescribe_velocity) then
                 block % provis_state % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)
              end if

              if (config_prescribe_thickness) then
                 block % provis_state % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)
              end if

              call ocn_diagnostic_solve(dt, block % provis_state, block % forcing, block % mesh, block % diagnostics, block % scratch)

              ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
              block % diagnostics % uTransport % array(:,:) &
                    = block % provis_state % normalVelocity % array(:,:) &
                    + block % diagnostics % uBolusGM   % array(:,:)

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
           block % state % time_levs(2) % state % normalVelocity % array(:,:) = block % state % time_levs(2) % state % normalVelocity % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % normalVelocity % array(:,:) 

           block % state % time_levs(2) % state % layerThickness % array(:,:) = block % state % time_levs(2) % state % layerThickness % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % layerThickness % array(:,:) 

           do iCell=1,block % mesh % nCells
              do k=1,block % mesh % maxLevelCell % array(iCell)
                 block % state % time_levs(2) % state % tracers % array(:,k,iCell) =  &
                                                                        block % state % time_levs(2) % state % tracers % array(:,k,iCell) &
                                                                        + rk_weights(rk_step) * block % tend % tracers % array(:,k,iCell)
              end do
           end do

           if (config_use_freq_filtered_thickness) then
                block % state % time_levs(2) % state % highFreqThickness % array(:,:) &
              = block % state % time_levs(2) % state % highFreqThickness % array(:,:) &
                + rk_weights(rk_step) * block % tend % highFreqThickness % array(:,:) 

                block % state % time_levs(2) % state % lowFreqDivergence % array(:,:) &
              = block % state % time_levs(2) % state % lowFreqDivergence % array(:,:) &
                + rk_weights(rk_step) * block % tend % lowFreqDivergence % array(:,:) 
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
        do iCell = 1, block % mesh % nCells
          do k = 1, block % mesh % maxLevelCell % array(iCell)
            block % state % time_levs(2) % state % tracers % array(:, k, iCell) = block % state % time_levs(2) % state % tracers % array(:, k, iCell) &
                                                                                / block % state % time_levs(2) % state % layerThickness % array(k, iCell)
          end do
        end do

        call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % forcing, block % mesh, block % diagnostics, block % scratch)
        call ocn_sea_ice_formation(block % mesh, block % state % time_levs(2) % state % index_temperature, &
                                   block % state % time_levs(2) % state % index_salinity, block % state % time_levs(2) % state % layerThickness % array, &
                                   block % state % time_levs(2) % state % tracers % array, block % forcing % seaIceEnergy % array, err)
        block => block % next
      end do

      call mpas_timer_start("RK4-implicit vert mix")
      block => domain % blocklist
      do while(associated(block))

        ! Call ocean diagnostic solve in preparation for vertical mixing.  Note 
        ! it is called again after vertical mixing, because u and tracers change.
        ! For Richardson vertical mixing, only density, layerThicknessEdge, and kineticEnergyCell need to 
        ! be computed.  For kpp, more variables may be needed.  Either way, this
        ! could be made more efficient by only computing what is needed for the
        ! implicit vmix routine that follows. 
        call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % forcing, block % mesh, block % diagnostics, block % scratch)

        call ocn_vmix_implicit(dt, block % mesh, block % diagnostics, block % state % time_levs(2) % state, err)
        block => block % next
      end do

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done, 
      ! this leads to lack of volume conservation.  It is required because halo updates in RK4 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to 
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("RK4-implicit vert mix halos")
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % normalVelocity)
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % tracers)
      call mpas_timer_stop("RK4-implicit vert mix halos")

      call mpas_timer_stop("RK4-implicit vert mix")

      block => domain % blocklist
      do while (associated(block))
         if (config_prescribe_velocity) then
            block % state % time_levs(2) % state % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)
         end if

         if (config_prescribe_thickness) then
            block % state % time_levs(2) % state % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)
         end if

         call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % forcing, block % mesh, block % diagnostics, block % scratch)

         ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
         block % diagnostics % uTransport % array(:,:) &
               = block % state % time_levs(2) % state % normalVelocity % array(:,:) &
               + block % diagnostics % uBolusGM % array(:,:)

         call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % normalVelocity % array,          &
                          block % diagnostics % normalVelocityX % array,            &
                          block % diagnostics % normalVelocityY % array,            &
                          block % diagnostics % normalVelocityZ % array,            &
                          block % diagnostics % normalVelocityZonal % array,        &
                          block % diagnostics % normalVelocityMeridional % array    &
                         )

         call mpas_reconstruct(block % mesh, block % diagnostics % gradSSH % array,          &
                          block % diagnostics % gradSSHX % array,            &
                          block % diagnostics % gradSSHY % array,            &
                          block % diagnostics % gradSSHZ % array,            &
                          block % diagnostics % gradSSHZonal % array,        &
                          block % diagnostics % gradSSHMeridional % array    &
                         )

         block % diagnostics % surfaceVelocity % array(block % diagnostics % index_zonalSurfaceVelocity, :) = &
               block % diagnostics % normalVelocityZonal % array(1, :)
         block % diagnostics % surfaceVelocity % array(block % diagnostics % index_meridionalSurfaceVelocity, :) = &
               block % diagnostics % normalVelocityMeridional % array(1, :)

         block % diagnostics % SSHGradient % array(block % diagnostics % index_zonalSSHGradient, :) = &
               block % diagnostics % gradSSHZonal % array(1, :)
         block % diagnostics % SSHGradient % array(block % diagnostics % index_meridionalSSHGradient, :) = &
               block % diagnostics % gradSSHMeridional % array(1, :)

         call ocn_time_average_accumulate(block % average, block % state % time_levs(2) % state, block % diagnostics)
         call ocn_time_average_coupled_accumulate(block % diagnostics, block % forcing)

         block => block % next
      end do
      call mpas_timer_stop("RK4-cleaup phase")

      call mpas_deallocate_provis_state(domain % blocklist)

   end subroutine ocn_time_integrator_rk4!}}}

end module ocn_time_integration_rk4

! vim: foldmethod=marker
