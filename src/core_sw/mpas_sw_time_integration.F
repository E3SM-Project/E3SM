! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module sw_time_integration

   use mpas_vector_reconstruction
   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar


   contains


   subroutine sw_timestep(domain, dt, timeStamp)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (block_type), pointer :: block

      if (trim(config_time_integration) == 'RK4') then
         call sw_rk4(domain, dt)
      else
         write(0,*) 'Unknown time integration option '//trim(config_time_integration)
         write(0,*) 'Currently, only ''RK4'' is supported.'
         stop
      end if

      block => domain % blocklist
      do while (associated(block))
         block % state % time_levs(2) % state % xtime % scalar = timeStamp 
         block => block % next
      end do

   end subroutine sw_timestep


   subroutine sw_rk4(domain, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      integer :: iCell, k
      type (block_type), pointer :: block
      type (state_type), target :: provis_state
      type (state_type), pointer :: provis_state_ptr

      integer :: rk_step

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      call mpas_setup_provis_state(domain % blocklist)
   
     !
     ! Initialize time_levs(2) with state at current time
     ! Initialize first RK state
     ! Couple tracers time_levs(2) with h in time-levels
     ! Initialize RK weights
     !
     block => domain % blocklist
     do while (associated(block))

        block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
        block % state % time_levs(2) % state % h % array(:,:) = block % state % time_levs(1) % state % h % array(:,:)
        do iCell=1,block % mesh % nCells  ! couple tracers to h
          do k=1,block % mesh % nVertLevels
            block % state % time_levs(2) % state % tracers % array(:,k,iCell) = block % state % time_levs(1) % state % tracers % array(:,k,iCell) &
                                                                      * block % state % time_levs(1) % state % h % array(k,iCell)
           end do
        end do

        call mpas_copy_state(block % provis_state, block % state % time_levs(1) % state)

        block => block % next
     end do

     rk_weights(1) = dt/6.
     rk_weights(2) = dt/3.
     rk_weights(3) = dt/3.
     rk_weights(4) = dt/6.

     rk_substep_weights(1) = dt/2.
     rk_substep_weights(2) = dt/2.
     rk_substep_weights(3) = dt
     rk_substep_weights(4) = 0.


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     ! BEGIN RK loop 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     do rk_step = 1, 4

! --- update halos for diagnostic variables

        call mpas_dmpar_exch_halo_field(domain % blocklist % provis_state % pv_edge)

        if (config_h_mom_eddy_visc4 > 0.0) then
            call mpas_dmpar_exch_halo_field(domain % blocklist % provis_state % divergence)
            call mpas_dmpar_exch_halo_field(domain % blocklist % provis_state % vorticity)
       end if

! --- compute tendencies

       block => domain % blocklist
       do while (associated(block))
          call sw_compute_tend(block % tend, block % provis_state, block % mesh)
          call sw_compute_scalar_tend(block % tend, block % provis_state, block % mesh)
          call sw_enforce_boundary_edge(block % tend, block % mesh)
          block => block % next
       end do

! --- update halos for prognostic variables

       call mpas_dmpar_exch_halo_field(domain % blocklist % tend % u)
       call mpas_dmpar_exch_halo_field(domain % blocklist % tend % h)
       call mpas_dmpar_exch_halo_field(domain % blocklist % tend % tracers)

! --- compute next substep state

       if (rk_step < 4) then
          block => domain % blocklist
          do while (associated(block))
             block % provis_state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)  &
                                                   + rk_substep_weights(rk_step) * block % tend % u % array(:,:)
             block % provis_state % h % array(:,:) = block % state % time_levs(1) % state % h % array(:,:)  &
                                                   + rk_substep_weights(rk_step) * block % tend % h % array(:,:)
             do iCell=1,block % mesh % nCells
                do k=1,block % mesh % nVertLevels
                   block % provis_state % tracers % array(:,k,iCell) = ( block % state % time_levs(1) % state % h % array(k,iCell) * &
                                                                         block % state % time_levs(1) % state % tracers % array(:,k,iCell)  &
                                   + rk_substep_weights(rk_step) * block % tend % tracers % array(:,k,iCell) &
                                                                       ) / block % provis_state % h % array(k,iCell)
                end do
             end do
             if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                block % provis_state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
             end if
             call sw_compute_solve_diagnostics(dt, block % provis_state, block % mesh)
             block => block % next
          end do
       end if

!--- accumulate update (for RK4)

       block => domain % blocklist
       do while (associated(block))
          block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(2) % state % u % array(:,:) &
                                  + rk_weights(rk_step) * block % tend % u % array(:,:) 
          block % state % time_levs(2) % state % h % array(:,:) = block % state % time_levs(2) % state % h % array(:,:) &
                                  + rk_weights(rk_step) * block % tend % h % array(:,:) 
          do iCell=1,block % mesh % nCells
             do k=1,block % mesh % nVertLevels
                block % state % time_levs(2) % state % tracers % array(:,k,iCell) =  &
                                                                      block % state % time_levs(2) % state % tracers % array(:,k,iCell) &
                                              + rk_weights(rk_step) * block % tend % tracers % array(:,k,iCell)
             end do
          end do
          block => block % next
       end do

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! END RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      !
      !  A little clean up at the end: decouple new scalar fields and compute diagnostics for new state
      !
      block => domain % blocklist
      do while (associated(block))
         do iCell=1,block % mesh % nCells
            do k=1,block % mesh % nVertLevels
               block % state % time_levs(2) % state % tracers % array(:,k,iCell) = &
                                                                     block % state % time_levs(2) % state % tracers % array(:,k,iCell)  &
                                                                   / block % state % time_levs(2) % state % h % array(k,iCell)
            end do
         end do

         if (config_test_case == 1) then    ! For case 1, wind field should be fixed
            block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
         end if

         call sw_compute_solve_diagnostics(dt, block % state % time_levs(2) % state, block % mesh)

         call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % u % array,          &
                          block % state % time_levs(2) % state % uReconstructX % array,            &
                          block % state % time_levs(2) % state % uReconstructY % array,            &
                          block % state % time_levs(2) % state % uReconstructZ % array,            &
                          block % state % time_levs(2) % state % uReconstructZonal % array,        &
                          block % state % time_levs(2) % state % uReconstructMeridional % array    &
                         )

         block => block % next
      end do

      call mpas_deallocate_provis_state(domain % blocklist)

   end subroutine sw_rk4


   subroutine sw_compute_tend(tend, s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed tendencies for prognostic variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      real (kind=RKIND) :: flux, vorticity_abs, workpv, q, upstream_bias

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
                                                  meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, divergence, h_vertex
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND) :: r, u_diffusion

      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity

      real (kind=RKIND), dimension(:,:), pointer :: u_src
      real (kind=RKIND), parameter :: rho_ref = 1000.0
      real (kind=RKIND) :: ke_edge


      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      h_edge      => s % h_edge % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      pv_edge     => s % pv_edge % array
      vh          => s % vh % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      h_s               => grid % h_s % array
      fVertex           => grid % fVertex % array
      fEdge             => grid % fEdge % array

      tend_h      => tend % h % array
      tend_u      => tend % u % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      u_src => grid % u_src % array

      meshScalingDel2 => grid % meshScalingDel2 % array
      meshScalingDel4 => grid % meshScalingDel4 % array


      !
      ! Compute height tendency for each cell
      !
      tend_h(:,:) = 0.0
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
            tend_h(k,cell1) = tend_h(k,cell1) - flux
            tend_h(k,cell2) = tend_h(k,cell2) + flux
         end do
      end do 
      do iCell=1,grid % nCellsSolve
         do k=1,nVertLevels
            tend_h(k,iCell) = tend_h(k,iCell) / areaCell(iCell)
         end do
      end do

#ifdef LANL_FORMULATION
      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !
      tend_u(:,:) = 0.0
      do iEdge=1,grid % nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         
         do k=1,nVertLevels
            q = 0.0
            do j = 1,nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
               q = q + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * h_edge(k,eoe) 
            end do

            tend_u(k,iEdge) =       &
                              q     &
                              - (   ke(k,cell2) - ke(k,cell1) + &
                                    gravity * (h(k,cell2) + h_s(cell2) - h(k,cell1) - h_s(cell1)) &
                                  ) / dcEdge(iEdge)
         end do
      end do


#endif

#ifdef NCAR_FORMULATION
      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !
      tend_u(:,:) = 0.0
      do iEdge=1,grid % nEdgesSolve
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,nVertLevels
            vorticity_abs = fEdge(iEdge) + (circulation(k,vertex1) + circulation(k,vertex2)) / &
                                           (areaTriangle(vertex1) + areaTriangle(vertex2))

            workpv = 2.0 * vorticity_abs / (h(k,cell1) + h(k,cell2))

            tend_u(k,iEdge) = workpv * vh(k,iEdge) - &
                              (ke(k,cell2) - ke(k,cell1) + &
                                 gravity * (h(k,cell2) + h_s(cell2) - h(k,cell1) - h_s(cell1)) &
                              ) / &
                              dcEdge(iEdge)
         end do
      end do
#endif

     ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
     !                    only valid for visc == constant
     if (config_h_mom_eddy_visc2 > 0.0) then
        do iEdge=1,grid % nEdgesSolve
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)

           do k=1,nVertLevels
              u_diffusion =   ( divergence(k,cell2)  -  divergence(k,cell1) ) / dcEdge(iEdge) &
                   -(vorticity(k,vertex2)  - vorticity(k,vertex1) ) / dvEdge(iEdge)
              u_diffusion = meshScalingDel2(iEdge) * config_h_mom_eddy_visc2 * u_diffusion
              tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion
           end do
        end do
     end if

     !
     ! velocity tendency: del4 dissipation, -\nu_4 \nabla^4 u
     !   computed as \nabla^2 u = \nabla divergence + k \times \nabla vorticity
     !   applied recursively.
     !   strictly only valid for h_mom_eddy_visc4 == constant
     !
     if (config_h_mom_eddy_visc4 > 0.0) then
        allocate(delsq_divergence(nVertLevels, nCells+1))
        allocate(delsq_u(nVertLevels, nEdges+1))
        allocate(delsq_circulation(nVertLevels, nVertices+1))
        allocate(delsq_vorticity(nVertLevels, nVertices+1))

        delsq_u(:,:) = 0.0

        ! Compute \nabla^2 u = \nabla divergence + k \times \nabla vorticity
        do iEdge=1,grid % nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)

           do k=1,nVertLevels

              delsq_u(k,iEdge) = ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                   -( vorticity(k,vertex2) - vorticity(k,vertex1)) / dvEdge(iEdge)

           end do
        end do

        ! vorticity using \nabla^2 u
        delsq_circulation(:,:) = 0.0
        do iEdge=1,nEdges
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)
           do k=1,nVertLevels
              delsq_circulation(k,vertex1) = delsq_circulation(k,vertex1) &
                   - dcEdge(iEdge) * delsq_u(k,iEdge)
              delsq_circulation(k,vertex2) = delsq_circulation(k,vertex2) &
                   + dcEdge(iEdge) * delsq_u(k,iEdge)
           end do
        end do
        do iVertex=1,nVertices
           r = 1.0 / areaTriangle(iVertex)
           do k=1,nVertLevels
              delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
           end do
        end do

        ! Divergence using \nabla^2 u
        delsq_divergence(:,:) = 0.0
        do iEdge=1,nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           do k=1,nVertLevels
              delsq_divergence(k,cell1) = delsq_divergence(k,cell1) &
                   + delsq_u(k,iEdge)*dvEdge(iEdge)
              delsq_divergence(k,cell2) = delsq_divergence(k,cell2) &
                   - delsq_u(k,iEdge)*dvEdge(iEdge)
           end do
        end do
        do iCell = 1,nCells
           r = 1.0 / areaCell(iCell)
           do k = 1,nVertLevels
              delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
           end do
        end do

        ! Compute - \kappa \nabla^4 u 
        ! as  \nabla div(\nabla^2 u) + k \times \nabla ( k \cross curl(\nabla^2 u) )
        do iEdge=1,grid % nEdgesSolve
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)

           do k=1,nVertLevels

              u_diffusion = (  delsq_divergence(k,cell2) &
                   - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                   -(  delsq_vorticity(k,vertex2) &
                   - delsq_vorticity(k,vertex1) ) / dvEdge(iEdge)

              u_diffusion = meshScalingDel4(iEdge) * config_h_mom_eddy_visc4 * u_diffusion
              tend_u(k,iEdge) = tend_u(k,iEdge) - u_diffusion

           end do
        end do

        deallocate(delsq_divergence)
        deallocate(delsq_u)
        deallocate(delsq_circulation)
        deallocate(delsq_vorticity)

     end if

     ! Compute u (velocity) tendency from wind stress (u_src)
     if(config_wind_stress) then
         do iEdge=1,grid % nEdges
            tend_u(1,iEdge) =  tend_u(1,iEdge) &
                  + u_src(1,iEdge)/rho_ref/h_edge(1,iEdge)
         end do
     endif

     if (config_bottom_drag) then
         do iEdge=1,grid % nEdges
             ! bottom drag is the same as POP:
             ! -c |u| u  where c is unitless and 1.0e-3.
             ! see POP Reference guide, section 3.4.4.
             ke_edge = 0.5 * ( ke(1,cellsOnEdge(1,iEdge)) &
                   + ke(1,cellsOnEdge(2,iEdge)))

             tend_u(1,iEdge) = tend_u(1,iEdge)  &
                  - 1.0e-3*u(1,iEdge) &
                  *sqrt(2.0*ke_edge)/h_edge(1,iEdge)
         end do
     endif
 
   end subroutine sw_compute_tend


   subroutine sw_compute_scalar_tend(tend, s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iCell, iEdge, k, iTracer, cell1, cell2, i
      real (kind=RKIND) :: flux, tracer_edge, r
      real (kind=RKIND) :: invAreaCell1, invAreaCell2, tracer_turb_flux
      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:,:), allocatable :: boundaryMask
      real (kind=RKIND), dimension(:,:,:), allocatable:: delsq_tracer
      
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, tracer_tend
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryCell
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order
      real (kind=RKIND), dimension(:,:), pointer :: u, h_edge

      u           => s % u % array
      h_edge      => s % h_edge % array
      dcEdge      => grid % dcEdge % array
      deriv_two   => grid % deriv_two % array
      dvEdge      => grid % dvEdge % array
      tracers     => s % tracers % array
      cellsOnEdge => grid % cellsOnEdge % array
      boundaryCell=> grid % boundaryCell % array
      boundaryEdge=> grid % boundaryEdge % array
      areaCell    => grid % areaCell % array
      tracer_tend => tend % tracers % array

      coef_3rd_order = 0.
      if (config_tracer_adv_order == 3) coef_3rd_order = 1.0
      if (config_tracer_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25


      tracer_tend(:,:,:) = 0.0

      if (config_tracer_adv_order == 2) then

      do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then
               do k=1,grid % nVertLevels
                  do iTracer=1,grid % nTracers
                     tracer_edge = 0.5 * (tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))
                     flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge) * tracer_edge
                     tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux/areaCell(cell1)
                     tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux/areaCell(cell2)
                  end do 
               end do 
            end if
      end do 

      else if (config_tracer_adv_order == 3) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,grid % nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  do iTracer=1,grid % nTracers
 
                     !-- if not a boundary cell
                     if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                        d2fdx2_cell1 = deriv_two(1,1,iEdge) * tracers(iTracer,k,cell1)
                        d2fdx2_cell2 = deriv_two(1,2,iEdge) * tracers(iTracer,k,cell2)

                        !-- all edges of cell 1
                        do i=1, grid % nEdgesOnCell % array (cell1)
                                d2fdx2_cell1 = d2fdx2_cell1 + &
                                deriv_two(i+1,1,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell1))
                        end do

                        !-- all edges of cell 2
                        do i=1, grid % nEdgesOnCell % array (cell2)
                                d2fdx2_cell2 = d2fdx2_cell2 + &
                                deriv_two(i+1,2,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell2))
                        end do

                     endif

                     !-- if u > 0:
                     if (u(k,iEdge) > 0) then
                        flux = dvEdge(iEdge) * u(k,iEdge) * h_edge(k,iEdge) * (          &
                             0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                             -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                     !-- else u <= 0:
                     else
                        flux = dvEdge(iEdge) *  u(k,iEdge) * h_edge(k,iEdge) * (          &
                             0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                             +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                     end if

                     !-- update tendency
                     tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux/areaCell(cell1)
                     tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux/areaCell(cell2)
                  enddo
               end do
            end if
         end do

      else  if (config_tracer_adv_order == 4) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if an edge is not on the outer-most ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,grid % nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  do iTracer=1,grid % nTracers

                     !-- if not a boundary cell
                     if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                        d2fdx2_cell1 = deriv_two(1,1,iEdge) * tracers(iTracer,k,cell1)
                        d2fdx2_cell2 = deriv_two(1,2,iEdge) * tracers(iTracer,k,cell2)

                        !-- all edges of cell 1
                        do i=1, grid % nEdgesOnCell % array (cell1)
                                d2fdx2_cell1 = d2fdx2_cell1 + &
                                deriv_two(i+1,1,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell1))
                        end do

                        !-- all edges of cell 2
                        do i=1, grid % nEdgesOnCell % array (cell2)
                                d2fdx2_cell2 = d2fdx2_cell2 + &
                                deriv_two(i+1,2,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell2))
                        end do

                     endif

                     flux = dvEdge(iEdge) *  u(k,iEdge) * h_edge(k,iEdge) * (          &
                          0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )

                     !-- update tendency
                     tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux/areaCell(cell1)
                     tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux/areaCell(cell2)
                  enddo
               end do
            end if
         end do

      endif   ! if (config_tracer_adv_order == 2 )

      !
      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      if ( config_h_tracer_eddy_diff2 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(grid % nVertLevels, grid % nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invAreaCell1 = 1.0/areaCell(cell1)
            invAreaCell2 = 1.0/areaCell(cell2)

            do k=1,grid % nVertLevels
              do iTracer=1, grid % nTracers
                 ! \kappa_2 \nabla \phi on edge
                 tracer_turb_flux = config_h_tracer_eddy_diff2 &
                    *( tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge)

                 ! div(h \kappa_2 \nabla \phi) at cell center
                 flux = dvEdge(iEdge) * h_edge(k,iEdge) * tracer_turb_flux * boundaryMask(k, iEdge)
                 tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) + flux * invAreaCell1
                 tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) - flux * invAreaCell2
              end do
            end do

         end do

        deallocate(boundaryMask)

      end if

      !
      ! tracer tendency: del4 horizontal tracer diffusion, &
      !    div(h \kappa_4 \nabla [div(h \nabla \phi)])
      !
      if ( config_h_tracer_eddy_diff4 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(grid % nVertLevels, grid % nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         allocate(delsq_tracer(grid % nTracers, grid % nVertLevels, grid % nCells+1))

         delsq_tracer(:,:,:) = 0.

         ! first del2: div(h \nabla \phi) at cell center
         do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,grid % nVertLevels
              do iTracer=1, grid % nTracers
                 delsq_tracer(iTracer,k,cell1) = delsq_tracer(iTracer,k,cell1) &
                    + dvEdge(iEdge) * h_edge(k,iEdge) * (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge) * boundaryMask(k,iEdge)
                 delsq_tracer(iTracer,k,cell2) = delsq_tracer(iTracer,k,cell2) &
                    - dvEdge(iEdge) * h_edge(k,iEdge) * (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge) * boundaryMask(k,iEdge)
              end do
            end do

         end do

         do iCell = 1, grid % nCells
            r = 1.0 / grid % areaCell % array(iCell)
            do k=1,grid % nVertLevels
            do iTracer=1,grid % nTracers
               delsq_tracer(iTracer,k,iCell) = delsq_tracer(iTracer,k,iCell) * r
            end do
            end do
         end do

         ! second del2: div(h \nabla [delsq_tracer]) at cell center
         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            invAreaCell1 = 1.0 / grid % areaCell % array(cell1)
            invAreaCell2 = 1.0 / grid % areaCell % array(cell2)

            do k=1,grid % nVertLevels
            do iTracer=1,grid % nTracers
               tracer_turb_flux = config_h_tracer_eddy_diff4 * (delsq_tracer(iTracer,k,cell2) - delsq_tracer(iTracer,k,cell1)) / dcEdge(iEdge)
               flux = dvEdge(iEdge) * tracer_turb_flux
               tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux * invAreaCell1 * boundaryMask(k,iEdge)
               tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux * invAreaCell2 * boundaryMask(k,iEdge)
            end do
            enddo

         end do

         deallocate(delsq_tracer)
         deallocate(boundaryMask)

      end if

   end subroutine sw_compute_scalar_tend


   subroutine sw_compute_solve_diagnostics(dt, s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: grid - grid metadata
   !
   ! Output: s - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      real (kind=RKIND), intent(in) :: dt
      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, workpv

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, divergence, &
                                                    h_vertex, vorticity_cell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, boundaryEdge, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND) :: r, h1, h2
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order

      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      vh          => s % vh % array
      h_edge      => s % h_edge % array
      h_vertex    => s % h_vertex % array
      tend_h      => s % h % array
      tend_u      => s % u % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      pv_edge     => s % pv_edge % array
      pv_vertex   => s % pv_vertex % array
      pv_cell     => s % pv_cell % array
      vorticity_cell => s % vorticity_cell % array
      gradPVn     => s % gradPVn % array
      gradPVt     => s % gradPVt % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      h_s               => grid % h_s % array
      fVertex           => grid % fVertex % array
      fEdge             => grid % fEdge % array
      deriv_two         => grid % deriv_two % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      boundaryEdge => grid % boundaryEdge % array
      boundaryCell => grid % boundaryCell % array

      !
      ! Find those cells that have an edge on the boundary
      !
      boundaryCell(:,:) = 0
      do iEdge=1,nEdges
       do k=1,nVertLevels
         if(boundaryEdge(k,iEdge).eq.1) then
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           boundaryCell(k,cell1) = 1
           boundaryCell(k,cell2) = 1
         endif
       enddo
      enddo

      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed h_edge value
      !

      coef_3rd_order = 0.
      if (config_thickness_adv_order == 3) coef_3rd_order = 1.0
      if (config_thickness_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      if (config_thickness_adv_order == 2) then

         do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then
               do k=1,grid % nVertLevels
                  h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
               end do 
            end if
         end do 

      else if (config_thickness_adv_order == 3) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,grid % nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  !-- if not a boundary cell
                  if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                     d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                     d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                     !-- all edges of cell 1
                     do i=1, grid % nEdgesOnCell % array (cell1)
                             d2fdx2_cell1 = d2fdx2_cell1 + &
                             deriv_two(i+1,1,iEdge) * h(k,grid % CellsOnCell % array (i,cell1))
                     end do

                     !-- all edges of cell 2
                     do i=1, grid % nEdgesOnCell % array (cell2)
                             d2fdx2_cell2 = d2fdx2_cell2 + &
                             deriv_two(i+1,2,iEdge) * h(k,grid % CellsOnCell % array (i,cell2))
                     end do

                  endif

                  !-- if u > 0:
                  if (u(k,iEdge) > 0) then
                     h_edge(k,iEdge) =     &
                          0.5*(h(k,cell1) + h(k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                          -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
                  !-- else u <= 0:
                  else
                     h_edge(k,iEdge) =     &
                          0.5*(h(k,cell1) + h(k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                          +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
                  end if

               end do   ! do k
            end if      ! if (cell1 <=
         end do         ! do iEdge

      else  if (config_thickness_adv_order == 4) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,grid % nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  !-- if not a boundary cell
                  if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                     d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                     d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                     !-- all edges of cell 1
                     do i=1, grid % nEdgesOnCell % array (cell1)
                             d2fdx2_cell1 = d2fdx2_cell1 + &
                             deriv_two(i+1,1,iEdge) * h(k,grid % CellsOnCell % array (i,cell1))
                     end do

                     !-- all edges of cell 2
                     do i=1, grid % nEdgesOnCell % array (cell2)
                             d2fdx2_cell2 = d2fdx2_cell2 + &
                             deriv_two(i+1,2,iEdge) * h(k,grid % CellsOnCell % array (i,cell2))
                     end do

                  endif

                  h_edge(k,iEdge) =   &
                       0.5*(h(k,cell1) + h(k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

               end do   ! do k
            end if      ! if (cell1 <=
         end do         ! do iEdge

      endif   ! if(config_thickness_adv_order == 2)

      !
      ! set the velocity in the nEdges+1 slot to zero, this is a dummy address
      !    used to when reading for edges that do not exist
      !
      u(:,nEdges+1) = 0.0

      !
      ! Compute circulation and relative vorticity at each vertex
      !
      circulation(:,:) = 0.0
      do iEdge=1,nEdges
         do k=1,nVertLevels
            circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)
            circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
         end do
      end do
      do iVertex=1,nVertices
         do k=1,nVertLevels
            vorticity(k,iVertex) = circulation(k,iVertex) / areaTriangle(iVertex)
         end do
      end do


      !
      ! Compute the divergence at each cell center
      !
      divergence(:,:) = 0.0
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCells) then
            do k=1,nVertLevels
              divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
            enddo
         endif
         if(cell2 <= nCells) then
            do k=1,nVertLevels
              divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
            enddo
         end if
      end do
      do iCell = 1,nCells
        r = 1.0 / areaCell(iCell)
        do k = 1,nVertLevels
           divergence(k,iCell) = divergence(k,iCell) * r
        enddo
      enddo

      !
      ! Compute kinetic energy in each cell
      !
      ke(:,:) = 0.0
      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k=1,nVertLevels
               ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
            end do
         end do
         do k=1,nVertLevels
            ke(k,iCell) = ke(k,iCell) / areaCell(iCell)
         end do
      end do

      !
      ! Compute v (tangential) velocities
      !
      v(:,:) = 0.0
      do iEdge = 1,nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1,nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
         end do
      end do

#ifdef NCAR_FORMULATION
      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      vh(:,:) = 0.0
      do iEdge=1,grid % nEdgesSolve
         do j=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j,iEdge)
            do k=1,nVertLevels
               vh(k,iEdge) = vh(k,iEdge) + weightsOnEdge(j,iEdge) * u(k,eoe) * h_edge(k,eoe)
            end do
         end do
      end do
#endif


      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells and distance-1 ghost cells )
      !
      do iVertex = 1,nVertices
         do k=1,nVertLevels
            h_vertex(k,iVertex) = 0.0
            do i=1,grid % vertexDegree
               h_vertex(k,iVertex) = h_vertex(k,iVertex) + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex(k,iVertex) = h_vertex(k,iVertex) / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex(k,iVertex)
         end do
      end do


      !
      ! Compute gradient of PV in the tangent direction
      !   ( this computes gradPVt at all edges bounding real cells and distance-1 ghost cells )
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / &
                              dvEdge(iEdge)
         enddo
      enddo

      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells )
      !
      pv_edge(:,:) = 0.0
      do iVertex = 1,nVertices
        do i=1,grid % vertexDegree
           iEdge = edgesOnVertex(i,iVertex)
           do k=1,nVertLevels
              pv_edge(k,iEdge) =  pv_edge(k,iEdge)  + 0.5 * pv_vertex(k,iVertex)
           end do
        end do
      end do


      !
      ! Modify PV edge with upstream bias. 
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt(k,iEdge)
         enddo
      enddo


      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells and distance-1 ghost cells )
      !
      pv_cell(:,:) = 0.0
      vorticity_cell(:,:) = 0.0
      do iVertex = 1, nVertices
       do i=1,grid % vertexDegree
         iCell = cellsOnVertex(i,iVertex)
         if (iCell <= nCells) then
           do k = 1,nVertLevels
             pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) / areaCell(iCell)
             vorticity_cell(k,iCell) = vorticity_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * vorticity(k, iVertex) / areaCell(iCell)
           enddo
         endif
       enddo
      enddo


      !
      ! Compute gradient of PV in normal direction
      !   ( this computes gradPVn for all edges bounding real cells )
      !
      gradPVn(:,:) = 0.0
      do iEdge = 1,nEdges
        if( cellsOnEdge(1,iEdge) <= nCells .and. cellsOnEdge(2,iEdge) <= nCells) then
          do k = 1,nVertLevels
            gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / &
                                 dcEdge(iEdge)
          enddo
        endif
      enddo

      ! Modify PV edge with upstream bias.
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * u(k,iEdge) * dt * gradPVn(k,iEdge)
         enddo
      enddo

      !
      ! set pv_edge = fEdge / h_edge at boundary points
      !
   !  if (maxval(boundaryEdge).ge.0) then
   !  do iEdge = 1,nEdges
   !     cell1 = cellsOnEdge(1,iEdge)
   !     cell2 = cellsOnEdge(2,iEdge)
   !     do k = 1,nVertLevels
   !       if(boundaryEdge(k,iEdge).eq.1) then
   !         v(k,iEdge) = 0.0
   !         if(cell1.gt.0) then
   !            h1 = h(k,cell1)
   !            pv_edge(k,iEdge) = fEdge(iEdge) / h1
   !            h_edge(k,iEdge) = h1
   !         else
   !            h2 = h(k,cell2)
   !            pv_edge(k,iEdge) = fEdge(iEdge) / h2
   !            h_edge(k,iEdge) = h2
   !         endif
   !       endif
   !     enddo
   !  enddo
   !  endif


   end subroutine sw_compute_solve_diagnostics


   subroutine sw_enforce_boundary_edge(tend, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Enforce any boundary conditions on the normal velocity at each edge
   !
   ! Input: grid - grid metadata
   !
   ! Output: tend_u set to zero at boundaryEdge == 1 locations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      type (tend_type), intent(inout) :: tend
      type (mesh_type), intent(in) :: grid

      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:,:), pointer :: tend_u
      integer :: nCells, nEdges, nVertices, nVertLevels
      integer :: iEdge, k

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      boundaryEdge         => grid % boundaryEdge % array
      tend_u               => tend % u % array

      if(maxval(boundaryEdge).le.0) return

      do iEdge = 1,nEdges
        do k = 1,nVertLevels

          if(boundaryEdge(k,iEdge).eq.1) then
             tend_u(k,iEdge) = 0.0
          endif

        enddo
       enddo

   end subroutine sw_enforce_boundary_edge


end module sw_time_integration
