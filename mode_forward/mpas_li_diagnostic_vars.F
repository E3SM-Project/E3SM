! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_diagnostic_vars
!
!> \brief MPAS land ice module for calculating diagnostic variables
!> \author Matt Hoffman
!> \date   17 April 2011
!> \details
!>  This module contains various subroutines for
!>  calculating diagnostic (time-independent) variables
!>  for the land ice core.  
!>  These calculations should be general so as to be independent
!>  of time integration scheme.
!
!-----------------------------------------------------------------------

module li_diagnostic_vars

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use li_velocity
   use li_mask

   implicit none
   private

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
   public :: li_calculate_diagnostic_vars, li_calculate_apparent_diffusivity, &
             li_calculate_layerThickness

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_calculate_diagnostic_vars
!
!> \brief   Calculates diagnostic variables 
!> \author  Matthew Hoffman
!> \date    12 September 2013
!> \details
!>  This routine calculates diagnostic variables using the current prognostic
!>  variables.  This is done in 3 parts:
!>  1. diagnostic solve part 1;  2. solve velocity;  3. diagnostic solve part 2
!>  Note: If the velocity solver requires an initial guess, it will be taken 
!>  from the current value.  Therefore the normalVelocity in that time level
!>  should be updated with the guess prior to calling this subroutine, if necessary.
!
!-----------------------------------------------------------------------

   subroutine li_calculate_diagnostic_vars(domain, solveVelo, err)

      use li_thermal, only: li_heat_dissipation_sia, li_basal_friction_sia

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      logical, intent(in) :: solveVelo  !< Input: Whether or not to solve velocity

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
      ! Note: domain is passed in because halo updates are needed in this routine
      ! and halo updates have to happen outside block loops, which requires domain.

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
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: velocityPool

      character (len=StrKIND), pointer :: config_thermal_solver
      character (len=StrKIND), pointer :: config_velocity_solver

      integer :: err_tmp


      err = 0

      ! ===
      ! === Diagnostic solve of variables prior to velocity
      ! ===

      call mpas_timer_start("calc. diagnostic vars before vel")
      call diagnostic_solve_before_velocity(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calc. diagnostic vars before vel")


      ! ===
      ! === Diagnostic solve of velocity
      ! ===
      call li_velocity_solve(domain, solveVelo, err)    ! ****** Calculate Velocity ******
      err = ior(err, err_tmp)


      ! ===
      ! === Diagnostic solve of variables after velocity
      ! ===

      call mpas_timer_start("calc. diagnostic vars after vel")

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call diagnostic_solve_after_velocity(meshPool, geometryPool, velocityPool, err_tmp)  ! Some diagnostic variables require velocity to compute
         err = ior(err, err_tmp)

         block => block % next
      end do

      ! Calculate the heating terms needed by the thermal solver during the next time step.
      ! Note: The SIA heat dissipation could be computed before velocity, but for other approximations, 
      !       heat dissipation generally must be computed after velocity.
      !TODO - Make sure the heatDissipation and basalFriction fields are computed for the FO solver.

      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      if (trim(config_thermal_solver) == 'temperature' .or. trim(config_thermal_solver) == 'enthalpy') then

         if (trim(config_velocity_solver) == 'sia') then

            ! Compute interior heat dissipation
            call li_heat_dissipation_sia(domain, err_tmp)
            err = ior(err, err_tmp)

            ! Compute heat flux due to basal friction
            call li_basal_friction_sia(domain, err_tmp)
            err = ior(err, err_tmp)

         endif   ! sia

      endif   ! temperature or enthalpy

      ! Halo updates
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'layerThicknessEdge')
      call mpas_timer_stop("halo updates")

      call mpas_timer_stop("calc. diagnostic vars after vel")

      ! === error check and exit
      if (err == 1) then
          print *, "An error has occurred in li_calculate_diagnostic_vars. Aborting..."
          !call mpas_dmpar_global_abort(dminfo)
      endif

   !--------------------------------------------------------------------
   end subroutine li_calculate_diagnostic_vars



!***********************************************************************
!
!  subroutine li_calculate_apparent_diffusivity
!
!> \brief   Computes apparent diffusivity
!> \author  Matt Hoffman
!> \date    19 April 2012
!> \details 
!> This routine computes the apparent diffusivity.
!> Estimate diffusivity using the relation that the 2-d flux Q=-D grad h and Q=UH,
!> where h is surface elevation, D is diffusivity, U is 2-d velocity vector, and H is thickness
!> Solving for D = UH/-grad h
!> DCFL: dt = 0.5 * dx**2 / D = 0.5 * dx**2 * slopemag / flux_downslope
!
!-----------------------------------------------------------------------
   subroutine li_calculate_apparent_diffusivity(meshPool, velocityPool, scratchPool, geometryPool, allowableDiffDt)
      use mpas_vector_reconstruction

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(in) :: &
         velocityPool          !< Input: velocity information

      type (mpas_pool_type), intent(in) :: &
         scratchPool          !< Input: scratch information

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         geometryPool          !< Input: geometry information

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real(kind=RKIND), intent(out) :: allowableDiffDt !< Output: allowable timestep based on diffusive CFL

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, pointer :: config_print_thickness_advection_info
      real (kind=RKIND), dimension(:), pointer :: normalSlopeEdge
      type (field1dReal), pointer :: slopeReconstructXField, slopeReconstructYField, slopeReconstructZField  !< Only needed for calling mpas_reconstruct, but not actually used here
      type (field1dReal), pointer :: slopeCellAxis1Field
      type (field1dReal), pointer :: slopeCellAxis2Field
      real (kind=RKIND), dimension(:), pointer :: slopeCellAxis1, slopeCellAxis2
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructAxis1, uReconstructAxis2
      real (kind=RKIND), dimension(:), pointer :: apparentDiffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, pointer :: nCells, nVertLevels
      real (kind=RKIND) :: allowableDtHere
      real (kind=RKIND) :: fluxVeloAxis1, fluxVeloAxis2
      real (kind=RKIND) :: fluxDownslope
      real (kind=RKIND) :: slopeCellMagnitude
      real (kind=RKIND) :: dCell
      integer :: iCell, iEdge, iLevel
      real (kind=RKIND), parameter :: bigNumber = 1.0e16_RKIND  ! This is ~300 million years in seconds, but it is small enough not too overflow
      real (kind=RKIND), parameter :: smallNumber = 1.0e-36
      logical :: divideSingularityFound

      ! Note: This routine could be broken into 2: one to calculate diffusivity
      ! and another to get the diffusive CFL timestep.  In that case, the first (and possibly the second)
      ! could be moved to diagnostic_variable_solve_after_velocity.  However, since
      ! diffusivity is only used for this check, I don't think it makes sense to separate these
      ! calculations for now.

      ! get needed variables
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(geometryPool, 'normalSlopeEdge', normalSlopeEdge)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness, timeLevel=1)
      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel=1)
      call mpas_pool_get_array(geometryPool, 'apparentDiffusivity', apparentDiffusivity)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(velocityPool, 'uReconstructZonal', uReconstructAxis1)
      call mpas_pool_get_array(velocityPool, 'uReconstructMeridional', uReconstructAxis2)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_field(scratchPool, 'workCell', slopeReconstructXField)
      call mpas_allocate_scratch_field(slopeReconstructXField, .true.)
      call mpas_pool_get_field(scratchPool, 'workCell2', slopeReconstructYField)
      call mpas_allocate_scratch_field(slopeReconstructYField, .true.)
      call mpas_pool_get_field(scratchPool, 'workCell3', slopeReconstructZField)
      call mpas_allocate_scratch_field(slopeReconstructZField, .true.)
      call mpas_pool_get_field(scratchPool, 'slopeCellX', slopeCellAxis1Field)
      call mpas_allocate_scratch_field(slopeCellAxis1Field, .true.)
      slopeCellAxis1 => slopeCellAxis1Field % array
      call mpas_pool_get_field(scratchPool, 'slopeCellY', slopeCellAxis2Field)
      call mpas_allocate_scratch_field(slopeCellAxis2Field, .true.)
      slopeCellAxis2 => slopeCellAxis2Field % array

      ! given thickness, compute layerThickness
      call li_calculate_layerThickness(meshPool, thickness, layerThickness)

      ! Initialize output
      allowableDiffDt = bigNumber

      ! Approximate slope at cell centers
      ! reconstruct routines set uReconstructZonal = uReconstructX; uReconstructMeridional = uReconstructY
      ! for planar meshes, so those variables can be used as orthogonal components of the vector
      ! in either the plane or sphere.  This avoids needing to add logic for if we are on a sphere or not.
      call mpas_reconstruct(meshPool, normalSlopeEdge,        &
                            slopeReconstructXField % array, slopeReconstructYField % array, slopeReconstructZField % array, &
                            slopeCellAxis1, slopeCellAxis2)


      ! Approximate flux at cell centers
      divideSingularityFound = .false.
      do iCell = 1, nCells
         slopeCellMagnitude = sqrt(slopeCellAxis1(iCell)**2 + slopeCellAxis2(iCell)**2) + smallNumber

         if ( (slopeCellMagnitude < 1.0e-4_RKIND) .and. &
              (max(maxval(uReconstructAxis1(:,iCell)), maxval(uReconstructAxis2(:,iCell))) < 3.18e-8_RKIND) ) then  ! 3.18e-8=1 m/yr in m/s
            ! Ignore diffusivity near 'divide-singularities'
            apparentDiffusivity(iCell) = 0.0_RKIND
            divideSingularityFound = .true.
         else
            fluxDownslope = 0.0_RKIND
            do iLevel = 1, nVertLevels
               fluxVeloAxis1 = (uReconstructAxis1(iLevel, iCell) + uReconstructAxis1(iLevel+1, iCell)) * 0.5_RKIND
               fluxVeloAxis2 = (uReconstructAxis2(iLevel, iCell) + uReconstructAxis2(iLevel+1, iCell)) * 0.5_RKIND
               fluxDownslope = fluxDownslope + &
                               (-1.0_RKIND * slopeCellAxis1(iCell) * fluxVeloAxis1 - slopeCellAxis2(iCell) * fluxVeloAxis2) &
                               * layerThickness(iLevel, iCell) / slopeCellMagnitude
            enddo
            apparentDiffusivity(iCell) = abs(fluxDownslope) / slopeCellMagnitude
         endif

         ! Calculate allowable timestep based on DCFL
         if ( li_mask_is_grounded_ice(cellMask(iCell)) .and. li_mask_is_dynamic_ice(cellMask(iCell)) ) then
            ! Find shortest distance to a neighboring cell center, dCell
            dCell = minval(dcEdge(1:nEdgesOnCell(iCell)))
            allowableDtHere = 0.5_RKIND * dCell**2 / (apparentDiffusivity(iCell) + smallNumber)
         else
            allowableDtHere = bigNumber
         endif
         allowableDiffDt = min(allowableDiffDt, allowableDtHere)
      enddo

      if (divideSingularityFound .and. config_print_thickness_advection_info) then
         write (stderrUnit,*) 'Notice: In calculating apparentDiffusivity, one or more cells have been ignored due to flat slope and low velocity (assumed to be a divide where diffusivity is undefined).'
      endif

      call mpas_deallocate_scratch_field(slopeReconstructXField, .true.)
      call mpas_deallocate_scratch_field(slopeReconstructYField, .true.)
      call mpas_deallocate_scratch_field(slopeReconstructZField, .true.)
      call mpas_deallocate_scratch_field(slopeCellAxis1Field, .true.)
      call mpas_deallocate_scratch_field(slopeCellAxis2Field, .true.)
   !--------------------------------------------------------------------
   end subroutine li_calculate_apparent_diffusivity



!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************


!***********************************************************************
!
!  subroutine diagnostic_solve_before_velocity
!
!> \brief   Computes diagnostic variables prior to velocity
!> \author  Matt Hoffman
!> \date    19 April 2012
!> \details
!> This routine computes the diagnostic variables for land ice
!> that are needed before velocity is solved.
!
!-----------------------------------------------------------------------
   subroutine diagnostic_solve_before_velocity(domain, err)!{{{

      use mpas_geometry_utils, only: mpas_cells_to_points_using_baryweights
      use mpas_vector_operations, only: mpas_tangential_vector_1d
      use li_setup, only: li_cells_to_vertices_1dfield_using_kiteAreas
      use li_thermal

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

      type (domain_type), intent(inout) :: &
         domain         !< Input/Output: Domain object

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
      ! pointers to get from pools
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: velocityPool
      real (kind=RKIND), dimension(:), pointer :: thickness, upperSurface, &
        lowerSurface, bedTopography, upperSurfaceVertex, slopeEdge, &
        normalSlopeEdge, tangentSlopeEdge, dcEdge, dvEdge
      real (kind=RKIND), dimension(:), pointer :: layerCenterSigma
      real (kind=RKIND), dimension(:,:), pointer :: pmpTemperature
      real (kind=RKIND), dimension(:), pointer :: basalPmpTemperature
      integer, dimension(:), pointer :: cellMask, edgeMask, vertexMask
      integer, dimension(:), pointer :: vertexMaskOld, vertexMaskNew
      integer, dimension(:), pointer :: floatingEdges
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge
      integer, dimension(:,:), pointer :: baryCellsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: baryWeightsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      real (kind=RKIND), dimension(:,:), pointer :: temperature
      integer, pointer :: nCells, nVertices, nEdges
      integer, pointer :: anyDynamicVertexMaskChanged
      integer, pointer :: dirichletMaskChanged
      integer, dimension(:,:), pointer :: dirichletVelocityMaskOld, dirichletVelocityMaskNew
      real (kind=RKIND), dimension(:), pointer :: beta, betaSolve
      real (kind=RKIND), pointer :: config_sea_level, config_ice_density, config_ocean_density
      character (len=StrKIND), pointer :: config_velocity_solver, config_sia_tangent_slope_calculation
      character (len=StrKIND), pointer :: config_thickness_advection, config_tracer_advection
      logical, pointer :: config_adaptive_timestep_include_DCFL
      logical, pointer :: config_SGH
      ! truly local variables
      real (kind=RKIND) :: thisThk
      integer :: iCell, iLevel, iEdge, cell1, cell2
      integer :: blockDynamicVertexMaskChanged, procDynamicVertexMaskChanged
      integer :: blockDirichletMaskChanged, procDirichletMaskChanged

      integer :: err_tmp

      err = 0
      err_tmp = 0

      block => domain % blocklist
      do while (associated(block))
         ! Mesh information
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         ! Calculate masks - needs to happen before calculating lower surface so we know where the ice is floating
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)

         err = ior(err, err_tmp)

         block => block % next
      end do

      ! Update halos on masks - the outermost cells/edges/vertices may be wrong for mask components that need neighbor information
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'cellMask')
      call mpas_dmpar_field_halo_exch(domain, 'edgeMask')
      call mpas_dmpar_field_halo_exch(domain, 'vertexMask')
      call mpas_timer_stop("halo updates")

!!!      ! Update beta before the velocity solve occurs, now that we have the new state and its mask.
!!!      !    It's ok to overwrite the beta values with 0's here, because each time step
!!!      !    we get a fresh copy of the array from the annual_forcing subroutine.
!!!      !    Note: some velocity solvers may do this on their own, but we are doing it here for completeness.
!!!      where ( MASK_IS_FLOATING(cellMask) )
!!!         beta = 0.0_RKIND
!!!      end where

!!!      ! Give non-ice cells a temperature of 0
!!!      do iCell = 1, nCells
!!!         if ( MASK_IS_NOT_ICE(cellMask(iCell)) ) then
!!!            state % tracers % array( state%index_temperature, :, iCell) = 0.0_RKIND
!!!         end if
!!!      end do

      block => domain % blocklist
      do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
          call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)

          call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
          call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
          call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
          call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
          call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)
          call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
          call mpas_pool_get_config(liConfigs, 'config_tracer_advection', config_tracer_advection)
          call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)

          call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

          call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
          call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
          call mpas_pool_get_array(geometryPool, 'thickness', thickness)
          call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
          call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
          call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

          call mpas_pool_get_array(thermalPool, 'temperature', temperature)
          call mpas_pool_get_array(thermalPool, 'pmpTemperature', pmpTemperature)
          call mpas_pool_get_array(thermalPool, 'basalPmpTemperature', basalPmpTemperature)

          ! Lower surface is based on floatation for floating ice.  For grounded ice (and non-ice areas) it is the bed.
          where ( li_mask_is_floating_ice(cellMask) )
             lowerSurface = config_sea_level - thickness * (config_ice_density / config_ocean_density)
          elsewhere
             lowerSurface = bedTopography
          end where
          ! Make sure lowerSurface calculation is reasonable.  This check could be deleted once this has been throroughly tested.
          do iCell = 1, nCells
             if (lowerSurface(iCell) < bedTopography(iCell)) then
                write (stderrUnit,*) 'lowerSurface less than bedTopography at cell:', iCell
                err = 1
             endif
          end do

          ! Upper surface is the lower surface plus the thickness
          upperSurface(:) = lowerSurface(:) + thickness(:)


          ! Calculate diagnostic PMP temperature
          ! Note: this is only to allow these fields to be output;
          !       the thermal module performs these calculations internally as needed.
          ! Note: The PMP temperature calculated here will differ from the PMP values
          !       used in the thermal module, because the thermal solver occurs before
          !       advection occurs, and this calculation occurs after.
          call li_compute_pressure_melting_point_fields(nCells, thickness, layerCenterSigma, &
                pmpTemperature, basalPmpTemperature)


          ! Calculate SIA-related variables, if needed
          !   This first block for calculating normalSlopeEdge is also needed for the DCFL calculation that could occur with any velocity solver
          if ( (trim(config_velocity_solver) == 'sia') .or. (config_adaptive_timestep_include_DCFL) ) then
             call mpas_pool_get_array(geometryPool, 'normalSlopeEdge', normalSlopeEdge)
             call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
             call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
             call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
             call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

             ! Calculate normal slope
             do iEdge = 1, nEdges
                ! Only calculate slope for edges that have ice on at least one side.
                if ( li_mask_is_dynamic_ice(edgeMask(iEdge)) ) then
                   cell1 = cellsOnEdge(1,iEdge)
                   cell2 = cellsOnEdge(2,iEdge)
                   ! Calculate slope at edge
                   normalSlopeEdge(iEdge) = (upperSurface(cell1) - upperSurface(cell2) ) / dcEdge(iEdge)
                else
                   normalSlopeEdge(iEdge) = 0.0_RKIND
                endif
             end do  ! edges
          endif

          ! flowParamA only needed by the SIA solver or hydro model
          if (trim(config_velocity_solver) == 'sia' .or. config_SGH) then
             call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)

             ! Calculate flowA
             call calculate_flowParamA(meshPool, temperature(:,:), thickness, flowParamA, err_tmp)
             err = ior(err, err_tmp)
          endif

          ! The remaining calculations are only needed by the SIA solver
          if (trim(config_velocity_solver) == 'sia') then
             call mpas_pool_get_config(liConfigs, 'config_sia_tangent_slope_calculation', config_sia_tangent_slope_calculation)
             call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
             call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
             call mpas_pool_get_array(geometryPool, 'slopeEdge', slopeEdge)
             call mpas_pool_get_array(geometryPool, 'tangentSlopeEdge', tangentSlopeEdge)
             call mpas_pool_get_array(geometryPool, 'upperSurfaceVertex', upperSurfaceVertex)
             call mpas_pool_get_array(meshPool, 'baryCellsOnVertex', baryCellsOnVertex)
             call mpas_pool_get_array(meshPool, 'baryWeightsOnVertex', baryWeightsOnVertex)

             ! Calculate upperSurfaceVertex if needed
             select case (trim(config_sia_tangent_slope_calculation))
             case ('from_vertex_barycentric')
                call mpas_cells_to_points_using_baryweights(meshPool, baryCellsOnVertex(:, 1:nVertices), &
                   baryWeightsOnVertex(:, 1:nVertices), upperSurface, upperSurfaceVertex(1:nVertices), err_tmp)
                err = ior(err, err_tmp)
             case ('from_vertex_barycentric_kiteareas')
                call li_cells_to_vertices_1dfield_using_kiteAreas(meshPool, upperSurface, upperSurfaceVertex)
             end select

             ! Calculate tangent slope
             select case (trim(config_sia_tangent_slope_calculation))
             case ('from_vertex_barycentric', 'from_vertex_barycentric_kiteareas')
                do iEdge = 1, nEdges
                   ! Only calculate slope for edges that have ice on at least one side.
                   if ( li_mask_is_dynamic_ice(edgeMask(iEdge)) ) then
                      tangentSlopeEdge(iEdge) = ( upperSurfaceVertex(verticesOnEdge(1,iEdge)) -  &
                            upperSurfaceVertex(verticesOnEdge(2,iEdge)) ) / dvEdge(iEdge)
                   else
                      tangentSlopeEdge(iEdge) = 0.0_RKIND
                   endif
                end do  ! edges
             case ('from_normal_slope')
                call mpas_tangential_vector_1d(normalSlopeEdge, meshPool, &
                       includeHalo=.true., tangentialVector=tangentSlopeEdge)
             case default
                write (stdErrUnit,*) 'Error: Invalid value for config_sia_tangent_slope_calculation.'
                err = 1
             end select

             ! Now calculate the slope magnitude
             slopeEdge = sqrt(normalSlopeEdge**2 + tangentSlopeEdge**2)

             ! Note: the outer halo may be wrong, but that's ok as long as numhalos>1 because the velocity on the 0-halo will still be correct.

          end if ! SIA variables

          err = ior(err, err_tmp)

          ! This information is only needed by external dycores.
          if ( (trim(config_velocity_solver) == 'L1L2') .or. &
           (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then
             ! The interface expects an array where 1's are floating edges and 0's are non-floating edges.
             call mpas_pool_get_array(velocityPool, 'floatingEdges', floatingEdges)
             floatingEdges = li_mask_is_floating_ice_int(edgeMask)
             call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMask, timeLevel=1)
             call li_calculate_extrapolate_floating_edgemask(meshPool, vertexMask, floatingEdges)
          end if

          block => block % next
      end do

      ! This halo update is only needed by external dycores.
      if ( (trim(config_velocity_solver) == 'L1L2') .or. &
           (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then
         ! Update halos on masks - the outermost cells/edges/vertices may be wrong for mask components that need neighbor information
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'floatingEdges')
         call mpas_timer_stop("halo updates")
      endif

      ! This information is only needed for some external dycores.
      if ( (trim(config_velocity_solver) == 'L1L2') .or. &
           (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then
         procDynamicVertexMaskChanged = 0
         procDirichletMaskChanged = 0

         ! Note: External dycores don't support multiple blocks per proc., but checking across
         !  blocks anyway, in case some day they do.
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)
            ! Determine if the vertex mask changed during this time step for this block (needed for external dycores)
            ! TODO:  there may be some aspects of the mask that are ok change for external dycores, but for now just check the whole thing.
            ! TODO: if we ever have more than one time level, then this logic should be revisited.
            call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMaskOld, timeLevel=2)
            call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMaskNew, timeLevel=1)
            if ( sum(li_mask_is_dynamic_ice_int(vertexMaskNew) - li_mask_is_dynamic_ice_int(vertexMaskOld)) /= 0 ) then
                blockDynamicVertexMaskChanged = 1
            else
                blockDynamicVertexMaskChanged = 0
            endif
            !print *, 'blockVertexMaskChanged ', blockVertexMaskChanged
            ! Determine if any blocks on this processor had a change to the vertex mask
            procDynamicVertexMaskChanged = max(procDynamicVertexMaskChanged, blockDynamicVertexMaskChanged)
            !print *,'procVertexMaskChanged', procVertexMaskChanged

            ! Also check to see if the Dirichlet b.c. mask has changed
            call mpas_pool_get_array(velocityPool, 'dirichletVelocityMask', dirichletVelocityMaskOld, timeLevel=2)
            call mpas_pool_get_array(velocityPool, 'dirichletVelocityMask', dirichletVelocityMaskNew, timeLevel=1)
            if ( sum(dirichletVelocityMaskNew - dirichletVelocityMaskOld) /= 0 ) then
                blockDirichletMaskChanged = 1
            else
                blockDirichletMaskChanged = 0
            endif
            ! Determine if any blocks on this processor had a change to the vertex mask
            procDirichletMaskChanged = max(procDirichletMaskChanged, blockDirichletMaskChanged)

            ! -- Set beta to 0 under floating ice --
            call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
            call mpas_pool_get_array(velocityPool, 'beta', beta)
            call mpas_pool_get_array(velocityPool, 'betaSolve', betaSolve)
            where (li_mask_is_floating_ice(cellMask))
               betaSolve = 0.0_RKIND
            elsewhere
               betaSolve = beta
            end where

            block => block % next
         end do

         ! Determine if the vertex mask has changed on any processor and store the value for later use (need to exit the block loop to do so)
         call mpas_pool_get_array(velocityPool, 'anyDynamicVertexMaskChanged', anyDynamicVertexMaskChanged)
         call mpas_dmpar_max_int(domain % dminfo, procDynamicVertexMaskChanged, anyDynamicVertexMaskChanged)
         !print *,'anyDynamicVertexMaskChanged', anyDynamicVertexMaskChanged
         ! Do the same for the Dirichlet b.c. mask
         call mpas_pool_get_array(velocityPool, 'dirichletMaskChanged', dirichletMaskChanged)
         call mpas_dmpar_max_int(domain % dminfo, procDirichletMaskChanged, dirichletMaskChanged)
         !print *,'dirichletMaskChanged', dirichletMaskChanged
      end if   ! higher-order

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in diagnostic_solve_before_velocity."
      endif

   !--------------------------------------------------------------------
   end subroutine diagnostic_solve_before_velocity



!***********************************************************************
!
!  subroutine diagnostic_solve_after_velocity
!
!> \brief   Computes diagnostic variables that require knowing velocity
!> \author  Matt Hoffman
!> \date    19 April 2012
!> \details 
!> This routine computes the diagnostic variables that require knowing velocity for land ice
!
!-----------------------------------------------------------------------
   subroutine diagnostic_solve_after_velocity(meshPool, geometryPool, velocityPool, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(in) :: &
         velocityPool          !< Input: velocity information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: &
         geometryPool         !< Input/Output: geometry info

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
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, layerThicknessEdge, normalVelocity
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nEdges, nVertLevels
      character (len=StrKIND), pointer :: config_thickness_advection
      integer :: iEdge, cell1, cell2, k
      real (kind=RKIND) :: VelSign

      !WHL - debug
!      integer, pointer :: nCells
!      integer, dimension(:), pointer :: indexToCellID
!      logical, pointer :: config_print_thickness_advection_info
!      integer :: iCell

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)

      call mpas_pool_get_array(geometryPool, 'thickness', thickness, timeLevel = 1)
      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel = 1)
      call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)

      !WHL - debug
!      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
!      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
!      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

      ! Calculate h_edge.  This is used by both thickness and tracer advection on the following Forward Euler time step.  
      ! Note: FO-Upwind thickness advection does not explicitly use h_edge but a FO h_edge is implied.
      ! Note: SIA velocity solver uses its own local calculation of h_edge that is always 2nd order.
      ! Note: ocn_diagnostic_solve in mpas_ocn_tendency.F has 2, 3, & 4th order calculations for h_edge that can be used.  
      ! NOTE: This calculates FO upwind h edge

      ! given thickness, compute layerThickness
      call li_calculate_layerThickness(meshPool, thickness, layerThickness)

      if (config_thickness_advection == 'fo') then

         ! If using FO-Upwind then h_edge must be FO.
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1, nVertLevels
               ! Calculate h on edges using first order
               VelSign = sign(1.0_RKIND, normalVelocity(k, iEdge))
               layerThicknessEdge(k,iEdge) = max(VelSign * layerThickness(k, cell1), VelSign * (-1.0_RKIND) * layerThickness(k, cell2)) 
               ! + velocity goes from index 1 to 2 in the cellsOnEdge array.  
               !  Doug does the calculation as: h_edge = max(VelSign, 0.0) * h1 - min(VelSign, 0.0) * h2
               !!! ! Calculate h on edges using second order
               !!! layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k, cell1) + layerThickness(k, cell2))
            end do
            ! thickness_edge is not currently in registry and not currenly needed.  If it is, uncomment the next line
            !h_edge = max(thickness(cell1), thickness(cell2))
            !!!h_edge = (thickness(k) + thickness(k) ) / 2.0  ! 2nd order 
         end do
      else
          !write(stdoutUnit,*) 'layerThicknessEdge not calculated!'
      endif

      ! Note: the outmost layerThicknessEdge may be wrong if its upwind cell is off this block - halo update should be done if this variable will be used.
      !       This halo update is done near the end of subroutine li_calculate_diagnostic_vars.

      !WHL - debug - Commented out, but might be useful for future debugging
!      if (config_print_thickness_advection_info) then
!         write(stderrUnit,*) ' '
!         write(stderrUnit,*) 'End of timestep: iCell (global), new thickness:'
!         do iCell = 1, nCells
!            if (thickness(iCell) > 0.0_RKIND) then
!               write(stderrUnit,*) indexToCellID(iCell), thickness(iCell)
!            endif
!         enddo
!      endif

   end subroutine diagnostic_solve_after_velocity


!***********************************************************************
!
!  subroutine calculate_flowParamA
!
!> \brief   Calculates the flow law parameter A based on temperature
!> \author  Matt Hoffman
!> \date    23 Jan 2014
!> \details
!>  This routine calculates the flow law parameter A based on temperature
!>  depending on what option is chosen.
!>  The default option is a constant A assigned from config_default_flowParamA.
!>  The PB1982 option uses this equation from \emph{Paterson and Budd} [1982]
!>  and \emph{Paterson} [1994] (copied from CISM):
!>   \[
!>   A(T^{*})=A0 \exp \left(\frac{-Q}{RT^{*}}\right)
!>   \]
!>   This is equation 9 in {\em Payne and Dongelmans}. $A0$ is a constant of proportionality,
!>   $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
!>   The pressure-corrected temperature, $T^{*}$ is given by:
!>   \[
!>   T^{*} = T - T_{pmp} + T_0
!>   \]
!>   \[
!>   T_{pmp} = T_0 - \sigma \rho g H \Phi
!>   \]
!>   $T$ is the ice temperature, $T_0$ is the triple point of water,
!>   $\rho$ is the ice density, and $\Phi$ is the (constant) rate of change of
!>   melting point temperature with pressure.
!>
!>  The CP2010 option uses this equation from the 4th Edition of Physics of Glaciers (Eq. 3.35):
!>   \[
!>   A(T^{*})=A0 \exp \left(\frac{-Q}{R} ( \frac{1}{T^{*}} - \frac{1}{T_t})\right)
!>   \]
!>   where the variables are the same as above and $T_t$ is the pressure corrected
!>   transition temperature (-10 deg C at 0 pressure).
!>   Values for $A0, Q, \Phi$ differ from PB1982.
!>
!>  All options are adjusted by the enhancement factor (which defaults to 1.0).
!-----------------------------------------------------------------------
   subroutine calculate_flowParamA(meshPool, temperature, thickness, flowParamA, err)
      use mpas_constants, only: gravity
      use li_constants, only: idealGasConstant

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         temperature    !< Input: temperature
      real (kind=RKIND), dimension(:), intent(in) :: &
         thickness    !< Input: thickness

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      integer, intent(inout) :: err

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         flowParamA    !< Input: flowParamA

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer, pointer :: nCells, nVertLevels
      character (len=StrKIND), pointer :: config_flowParamA_calculation
      real (kind=RKIND), pointer :: config_default_flowParamA, &
                                    config_enhancementFactor,  &
                                    config_dynamic_thickness,  &
                                    config_ice_density
      integer :: iCell, iLevel, err_tmp
      real (kind=RKIND), dimension(:), pointer :: layerCenterSigma
      real (kind=RKIND) :: A0, Q, pressureMeltPointSlope
      real (kind=RKIND) :: temperatureCorrected, transitionTemperatureCorrected

      err_tmp = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

      call mpas_pool_get_config(liConfigs, 'config_flowParamA_calculation', config_flowParamA_calculation)
      call mpas_pool_get_config(liConfigs, 'config_enhancementFactor', config_enhancementFactor)
      call mpas_pool_get_config(liConfigs, 'config_default_flowParamA', config_default_flowParamA)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)


      select case(config_flowParamA_calculation)
      !-----------------------------------------------------------------
      case('constant')
          flowParamA = config_default_flowParamA
      !-----------------------------------------------------------------
      case('PB1982')
          pressureMeltPointSlope = 9.7456e-8_RKIND
          do iCell = 1, nCells
            if (thickness(iCell) > 0.0_RKIND) then
              do iLevel = 1, nVertLevels
                ! Calculate the pressure-corrected temperature
                temperatureCorrected = min(273.15_RKIND, temperature(iLevel,iCell) + pressureMeltPointSlope * &
                         thickness(iCell) * config_ice_density * gravity * layerCenterSigma(iLevel) )
                temperatureCorrected = max(223.15_RKIND, temperatureCorrected)
                ! Calculate flow A
                if (temperatureCorrected > 263.15_RKIND) then
                  A0 = 1.733e3_RKIND
                  Q  = 139.0e3_RKIND
                else
                  A0 = 3.613e-13_RKIND
                  Q  = 60.0e3_RKIND
                endif
                flowParamA(iLevel,iCell) = A0 * exp(-1.0_RKIND * Q / (idealGasConstant * temperatureCorrected))
              enddo ! levels
            else
              flowParamA(:,iCell) = 0.0_RKIND  ! non-ice cells get 0
            endif ! if dynamic ice
          enddo ! cells
      !-----------------------------------------------------------------
      case('CP2010')
          pressureMeltPointSlope = 7.0e-8_RKIND
          do iCell = 1, nCells
            if (thickness(iCell) > 0.0_RKIND) then  ! SIA solver could make use of A on thin ice if doing 2nd order averaging of flwa onto edges (otherwise this could be the dynamic thickness limit)
              do iLevel = 1, nVertLevels
                ! Calculate the pressure-corrected temperature
                temperatureCorrected = min(273.15_RKIND, temperature(iLevel,iCell) + pressureMeltPointSlope * &
                         thickness(iCell) * config_ice_density * gravity * layerCenterSigma(iLevel) )
                temperatureCorrected = max(223.15_RKIND, temperatureCorrected)
                transitionTemperatureCorrected = 263.15_RKIND + pressureMeltPointSlope * &
                         thickness(iCell) * config_ice_density * gravity * layerCenterSigma(iLevel)
                ! Calculate flow A
                A0 = 3.5e-25_RKIND
                if (temperatureCorrected > 263.15_RKIND) then
                  Q  = 115.0e3_RKIND
                else
                  Q  = 6.0e4_RKIND
                endif
                flowParamA(iLevel,iCell) = A0 * exp(-1.0_RKIND * Q / idealGasConstant * (1.0_RKIND/temperatureCorrected - 1.0_RKIND/transitionTemperatureCorrected))
              enddo ! levels
            else
              flowParamA(:,iCell) = 0.0_RKIND  ! non-ice cells get 0
            endif ! if dynamic ice
          enddo ! cells
      !-----------------------------------------------------------------
      end select

      !print *,'max flwa', maxval(flowParamA)
      !print *,'config_enhancementFactor', config_enhancementFactor

      ! Include enhancement factor
      flowParamA = flowParamA * config_enhancementFactor

      err = ior(err, err_tmp)

   end subroutine calculate_flowParamA

!***********************************************************************
!
!  subroutine li_calculate_layerThickness
!
!> \brief   Calculates the thickness of each layer, given the total thickness
!> \author  William Lipscomb
!> \date    16 Feb 2016
!> \details
!> This routine calculates layerThickness in each cell and column,
!> given thickness in each cell. The calculation is very simple, but is
!> packaged in a subroutine to ensure that it is done the same way 
!> in different parts of the code.
!-----------------------------------------------------------------------

   subroutine li_calculate_layerThickness(meshPool, thickness, layerThickness)

     !-----------------------------------------------------------------
     !
     ! input variables
     ! 
     !-----------------------------------------------------------------

     type (mpas_pool_type), intent(in) :: &
          meshPool                                     !< Input: mesh object

     real(kind=RKIND), dimension(:), intent(in) :: &
          thickness                                    !< Input: ice thickness

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

     real(kind=RKIND), dimension(:,:), intent(out) :: &
          layerThickness                               !< Output: thickness of each layer

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     integer, pointer :: &
          nCells,                    & ! number of cells
          nVertLevels                  ! number of vertical layers

     real (kind=RKIND), dimension(:), pointer :: &
          layerThicknessFractions      ! fractional thickness in each layer

     integer :: iCell, k

     call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
     call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

     call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)

     do iCell = 1, nCells
        do k = 1, nVertLevels
           layerThickness(k,iCell) = thickness(iCell) * layerThicknessFractions(k)
        enddo
     enddo

   end subroutine li_calculate_layerThickness


end module li_diagnostic_vars

