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
   public :: li_calculate_diagnostic_vars, li_calculate_apparent_diffusivity

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
      type (field2DReal), pointer :: layerThicknessEdgeField
      integer :: err_tmp


      err = 0

      ! ===
      ! === Diagnostic solve of variables prior to velocity 
      ! ===

      call mpas_timer_start("calc. diagnostic vars except vel")
      call diagnostic_solve_before_velocity(domain, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("calc. diagnostic vars except vel")


      ! ===
      ! === Diagnostic solve of velocity 
      ! ===
      if (solveVelo) then
          call li_velocity_solve(domain, err)    ! ****** Calculate Velocity ******
          err = ior(err, err_tmp)
      endif


      ! ===
      ! === Diagnostic solve of variables after velocity 
      ! ===

      call mpas_timer_start("calc. diagnostic vars except vel")

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call diagnostic_solve_after_velocity(meshPool, geometryPool, velocityPool, err_tmp)  ! Some diagnostic variables require velocity to compute
         err = ior(err, err_tmp)

         block => block % next
      end do

      call mpas_timer_start("halo updates")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_field(geometryPool, 'layerThicknessEdge', layerThicknessEdgeField)
      call mpas_dmpar_exch_halo_field(layerThicknessEdgeField)
      call mpas_timer_stop("halo updates")

      call mpas_timer_stop("calc. diagnostic vars except vel")

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
      real (kind=RKIND), dimension(:), pointer :: normalSlopeEdge
      type (field1dReal), pointer :: slopeReconstructXField, slopeReconstructYField, slopeReconstructZField  !< Only needed for calling mpas_reconstruct, but not actually used here
      type (field1dReal), pointer :: slopeCellAxis1Field
      type (field1dReal), pointer :: slopeCellAxis2Field
      real (kind=RKIND), dimension(:), pointer :: slopeCellAxis1, slopeCellAxis2
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
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(geometryPool, 'normalSlopeEdge', normalSlopeEdge)
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
               fluxDownslope = fluxDownslope + (-1.0_RKIND * slopeCellAxis1(iCell) * fluxVeloAxis1 - slopeCellAxis2(iCell) * fluxVeloAxis2) * layerThickness(iLevel, iCell) / slopeCellMagnitude
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

      if (divideSingularityFound) then
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
      
      !TODO - Move li_heat_dissipation_sia to a different module?
      use li_thermal, only: li_heat_dissipation_sia

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
      integer, dimension(:), pointer :: cellMask, edgeMask, vertexMask
      integer, dimension(:), pointer :: vertexMaskOld, vertexMaskNew
      integer, dimension(:), pointer :: floatingEdges
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge
      integer, dimension(:,:), pointer :: baryCellsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, baryWeightsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer, pointer :: index_temperature
      type (field1DInteger), pointer :: cellMaskField, edgeMaskField, vertexMaskField, floatingEdgesField
      integer, pointer :: nCells, nVertices, nEdges
      integer, pointer :: anyDynamicVertexMaskChanged
      integer, pointer :: dirichletMaskChanged
      integer, dimension(:,:), pointer :: dirichletVelocityMaskOld, dirichletVelocityMaskNew
      real (kind=RKIND), dimension(:), pointer :: beta, betaSolve
      real (kind=RKIND), pointer :: config_sea_level, config_ice_density, config_ocean_density
      character (len=StrKIND), pointer :: config_velocity_solver, config_sia_tangent_slope_calculation
      character (len=StrKIND), pointer :: config_thermal_solver
      logical, pointer :: config_adaptive_timestep_include_DCFL
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
      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_field(geometryPool, 'cellMask', cellMaskField)
      call mpas_pool_get_field(geometryPool, 'edgeMask', edgeMaskField)
      call mpas_pool_get_field(geometryPool, 'vertexMask', vertexMaskField, timeLevel=1)
      call mpas_dmpar_exch_halo_field(cellMaskField)
      call mpas_dmpar_exch_halo_field(edgeMaskField)
      call mpas_dmpar_exch_halo_field(vertexMaskField)
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
          call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
          call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)

          call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
          call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
          call mpas_pool_get_array(geometryPool, 'thickness', thickness)
          call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
          call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
          call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel = 1)
          call mpas_pool_get_array(thermalPool, 'tracers', tracers)
          call mpas_pool_get_dimension(thermalPool, 'index_temperature', index_temperature)
          call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)


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
          ! The remaining calculations are only needed by the SIA solver
          if (trim(config_velocity_solver) == 'sia') then
             call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
             call mpas_pool_get_config(liConfigs, 'config_sia_tangent_slope_calculation', config_sia_tangent_slope_calculation)
             call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
             call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
             call mpas_pool_get_array(geometryPool, 'slopeEdge', slopeEdge)
             call mpas_pool_get_array(geometryPool, 'tangentSlopeEdge', tangentSlopeEdge)
             call mpas_pool_get_array(geometryPool, 'upperSurfaceVertex', upperSurfaceVertex)
             call mpas_pool_get_array(meshPool, 'baryCellsOnVertex', baryCellsOnVertex)
             call mpas_pool_get_array(meshPool, 'baryWeightsOnVertex', baryWeightsOnVertex)
             call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)

             ! Calculate flowA
             call calculate_flowParamA(meshPool, tracers(index_temperature,:,:), thickness, flowParamA, err_tmp)
             err = ior(err, err_tmp)

             ! Calculate upperSurfaceVertex if needed
             select case (trim(config_sia_tangent_slope_calculation))
             case ('from_vertex_barycentric')
                call mpas_cells_to_points_using_baryweights(meshPool, baryCellsOnVertex(:, 1:nVertices), &
                   baryWeightsOnVertex(:, 1:nVertices), upperSurface, upperSurfaceVertex(1:nVertices), err_tmp)
                err = ior(err, err_tmp)
             case ('from_vertex_barycentric_kiteareas')
                call cells_to_vertices_1dfield_using_kiteAreas(meshPool, upperSurface, upperSurfaceVertex)
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


          ! Do vertical remapping of layerThickness and tracers
          !WHL - I think the last argument should be err_tmp.
          call vertical_remap(thickness, cellMask, meshPool, layerThickness, tracers, err)
          err = ior(err, err_tmp)

          ! This information is only needed by external dycores.
          if (config_velocity_solver /= 'sia') then
             ! The interface expects an array where 1's are floating edges and 0's are non-floating edges.
             call mpas_pool_get_array(velocityPool, 'floatingEdges', floatingEdges)
             floatingEdges = li_mask_is_floating_ice_int(edgeMask)
             call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMask, timeLevel=1)
             call li_calculate_extrapolate_floating_edgemask(meshPool, vertexMask, floatingEdges)
          end if

          block => block % next
      end do

      ! This halo update is only needed by external dycores.
      if (config_velocity_solver /= 'sia') then
         ! Update halos on masks - the outermost cells/edges/vertices may be wrong for mask components that need neighbor information
         call mpas_timer_start("halo updates")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)
         call mpas_pool_get_field(velocityPool, 'floatingEdges', floatingEdgesField)
         call mpas_dmpar_exch_halo_field(floatingEdgesField)
         call mpas_timer_stop("halo updates")
      endif

      ! This information is only needed for some external dycores.
      if (trim(config_velocity_solver) /= 'sia') then
         procDynamicVertexMaskChanged = 0

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
      end if

      ! Calculate heat dissipation, as needed by the thermal solver during the next time step.
      !TODO - Make sure heat dissipation is computed for the FO solver.

      if (trim(config_thermal_solver) == 'temperature' .or. trim(config_thermal_solver) == 'enthalpy') then

         if (trim(config_velocity_solver) == 'sia') then
            call li_heat_dissipation_sia(domain, err_tmp)
            err = ior(err, err_tmp)
         endif

      endif

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
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, layerThicknessEdge, normalVelocity
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nEdges, nVertLevels
      character (len=StrKIND), pointer :: config_thickness_advection
      integer :: iEdge, cell1, cell2, k
      real (kind=RKIND) :: VelSign

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)

      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel = 1)
      call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)

      ! Calculate h_edge.  This is used by both thickness and tracer advection on the following Forward Euler time step.  
      ! Note: FO-Upwind thickness advection does not explicitly use h_edge but a FO h_edge is implied.
      ! Note: SIA velocity solver uses its own local calculation of h_edge that is always 2nd order.
      ! Note: ocn_diagnostic_solve in mpas_ocn_tendency.F has 2, 3, & 4th order calculations for h_edge that can be used.  
      ! NOTE: This calculates FO upwind h edge

      ! Both thickness and layerThickness should be updated by this time.
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

   end subroutine diagnostic_solve_after_velocity



!***********************************************************************
!
!  subroutine vertical_remap_cism_loops
!
!> \brief   Vertical remapping of thickness and tracers
!> \author  Matt Hoffman
!> \date    06 February 2013
!> \details 
!>  This routine performs vertical remapping of thickness and tracers from one vertical
!>  coordinate system to another, as is required for our sigma coordinate system.  
!>  The remapping is first-order accurate.
!>  This uses code from the CISM glissade_transport.F90 module written by Bill Lipscomb.  
!>  I have altered the array structures to work with MPAS.  Indexing/looping order is a bit
!>  of a hodgepodge at the moment and should be optimized.
!>
!>  NOTE: This version (with _cism_loops suffix) uses the looping with minimal changes
!>  from the CISM code.  The other version has been refactored to be more efficient
!>  for MPAS.  This version is intended to be deleted once tracer advection is
!>  added to the model and the two versions can be compared.  (I can't test the new
!>  version until tracer advection exists!)
!
!-----------------------------------------------------------------------
   subroutine vertical_remap_cism_loops(layerThickness, thickness, tracers, meshPool, err)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         thickness    !< Input: 

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: LI mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThickness    !< Input:

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers    !< Input:

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
      ! pointers to mesh arrays
      real (kind=RKIND), dimension(:), pointer :: layerThicknessFractions, layerInterfaceSigma
      ! local arrays
      real (kind=RKIND), dimension(:), allocatable :: recipThickness
      real (kind=RKIND), dimension(:,:), allocatable :: layerInterfaceSigma_Input
      real (kind=RKIND), dimension(:,:,:), allocatable :: hTsum
      ! counters, mesh variables, index variables
      integer, pointer :: nCells, nVertLevels
      integer :: nTracers, iCell, k, k1, k2, nt
      ! stuff for making calculations
      real(kind=RKIND) :: thisThk, zhi, zlo, hOverlap

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tracers, 1)

      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)

      allocate(recipThickness(nCells+1))
      allocate(layerInterfaceSigma_Input(nVertLevels+1, nCells+1))
      allocate(hTsum(nCells+1, nTracers, nVertLevels))

      ! *** Calculate reciprocal thickness so we don't divide by 0
      where (thickness > 0.0_RKIND)
         recipThickness = 1.0_RKIND / thickness
      elsewhere
         recipThickness = 0.0_RKIND
      end where

      ! *** Calculate vertical sigma coordinates of each layer interface for the input non-sigma state and desired new sigma-based state (we already have that as mesh % layerInterfaceSigma)
      layerInterfaceSigma_Input(1,:) = 0.0_RKIND 
      do k = 2, nVertLevels
        layerInterfaceSigma_Input(k,:) = layerInterfaceSigma_Input(k-1,:) + layerThickness(k-1,:) * recipThickness(:)
      end do
      layerInterfaceSigma_Input(nVertLevels+1,:) = 1.0_RKIND 

      ! *** Compute new layer thicknesses (layerInterfaceSigma coordinates)
      do iCell = 1, nCells
         thisThk = thickness(iCell)
         do k = 1, nVertLevels
           layerThickness(k,iCell) = layerThicknessFractions(k) * thisThk
         end do
         ! Check for conservation of mass.  Put any residual in the top layer.
         layerThickness(1,iCell) = layerThickness(1,iCell) + (thisThk - sum(layerThickness(:,iCell)) )
      end do
      ! TODO This conservation check may make layerThicknesses inconsistent with the sigma levels (which are used below in tracer remapping.

       !-----------------------------------------------------------------
       ! Compute sum of h*T for each new layer (k2) by integrating
       ! over the regions of overlap with old layers (k1).
       ! Note: It might be worth trying a more efficient
       !       search algorithm if the number of layers is large.
       !       This algorithm scales as nlyr^2.
       !       Also, may want to rearrange loop order if there are many tracers.
       !-----------------------------------------------------------------

       do k2 = 1, nVertLevels
          hTsum(:,:,k2) = 0.d0 
          do k1 = 1, nVertLevels
             do nt = 1, nTracers
                do iCell = 1, nCells
                      zhi = min (layerInterfaceSigma_Input(k1+1,iCell), layerInterfaceSigma(k2+1)) 
                      zlo = max (layerInterfaceSigma_Input(k1,iCell), layerInterfaceSigma(k2))
                      hOverlap = max (zhi-zlo, 0.d0) * thickness(iCell)
                      hTsum(iCell,nt,k2) = htsum(iCell,nt,k2)    &
                                       +  tracers(nt,k1,iCell) * hOverlap
                enddo      ! iCell
             enddo         ! nt
          enddo            ! k1
       enddo               ! k2
 
       !-----------------------------------------------------------------
       ! Compute tracer values in new layers
       !-----------------------------------------------------------------
 
       do k = 1, nVertLevels
          do nt = 1, nTracers
                do iCell = 1, nCells
                   if (layerThickness(k, iCell) > 0.0_RKIND) then
                      tracers(nt,k,iCell) = hTsum(iCell,nt,k) / layerThickness(k, iCell)
                   else
                      tracers(nt,k,iCell) = 0.0_RKIND
                   endif
                enddo   ! iCell
          enddo         ! nt
       enddo            ! k

      deallocate(recipThickness)
      deallocate(layerInterfaceSigma_Input)
      deallocate(hTsum)

   end subroutine vertical_remap_cism_loops




!***********************************************************************
!
!  subroutine vertical_remap
!
!> \brief   Vertical remapping of thickness and tracers
!> \author  Matt Hoffman
!> \date    02 October 2013
!> \details 
!>  This routine performs vertical remapping of thickness and tracers from one vertical
!>  coordinate system to another, as is required for our sigma coordinate system.
!>  The remapping is first-order accurate.
!>  This uses code from the CISM glissade_transport.F90 module written by Bill Lipscomb.
!>  I have altered the array structures to work with MPAS and refactored it.
!>  It now does all calculations column-wise, so it can be vectorized using
!>  OpenMP over either blocks or cells.
!>  I have dealt with 0 thickness in non-ice cells using eps and the cellMask,
!>  rather than using if/where-statements.
!
!-----------------------------------------------------------------------
   subroutine vertical_remap(thickness, cellMask, meshPool, layerThickness, tracers, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         thickness    !< Input: ice thickness (already updated)

      integer, dimension(:), intent(in) :: &
         cellMask    !< Input: mask for cells (needed for determining presence/absence of ice)

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: LI mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThickness    !< Input/Output: thickness of layers (to be updated)

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers    !< Input/Output: tracer values (to be updated)

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
      ! pointers to mesh arrays
      real (kind=RKIND), dimension(:), pointer :: layerThicknessFractions, layerInterfaceSigma
      ! local arrays
      real (kind=RKIND), dimension(:), allocatable :: layerInterfaceSigma_Input
      real (kind=RKIND), dimension(:,:), allocatable :: hTsum
      ! counters, mesh variables, index variables
      integer, pointer :: nCells, nVertLevels
      integer :: nTracers, iCell, k, k1, k2, nt
      ! stuff for making calculations
      real(kind=RKIND) :: thisThk, zhi, zlo, hOverlap
      real(kind=RKIND), parameter :: eps = 1.0d-30  ! Make a tiny value so we don't corrupt any cells with ice

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tracers, 1)

      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)

      allocate(layerInterfaceSigma_Input(nVertLevels+1))
      allocate(hTsum(nTracers, nVertLevels))

      do iCell = 1, nCells  ! Do all calculations column-wise

         thisThk = thickness(iCell)

         ! *** Calculate vertical sigma coordinates of each layer interface for the input non-sigma state 
         !     (we already have the desired new sigma-based state as mesh % layerInterfaceSigma)
         layerInterfaceSigma_Input(1) = 0.0_RKIND
         do k = 2, nVertLevels
            layerInterfaceSigma_Input(k) = layerInterfaceSigma_Input(k-1) + layerThickness(k-1, iCell) / (thisThk + eps)
            ! Note: including eps to avoid divide by 0, but need to mask out cells with no ice below in calculating
            !       layerThickness and tracer values.
         end do
         layerInterfaceSigma_Input(nVertLevels+1) = 1.0_RKIND

         ! *** Compute new layer thicknesses (layerInterfaceSigma coordinates)

         do k = 1, nVertLevels
            layerThickness(k,iCell) = layerThicknessFractions(k) * (thisThk + eps)
            ! Note: layerThickness will have garbage values in non-ice cells due to eps
         end do
         ! Check for conservation of mass.  Put any residual in the top layer.
         layerThickness(1,iCell) = layerThickness(1,iCell) + (thisThk - sum(layerThickness(:,iCell)) )
         ! TODO This conservation check may make layerThicknesses inconsistent with the sigma levels (which are used below in tracer remapping.

         !-----------------------------------------------------------------
         ! Compute sum of h*T for each new layer (k2) by integrating
         ! over the regions of overlap with old layers (k1).
         ! Note: It might be worth trying a more efficient
         !       search algorithm if the number of layers is large.
         !       This algorithm scales as nlyr^2.
         !-----------------------------------------------------------------

         do k2 = 1, nVertLevels
            hTsum(:,k2) = 0.d0
            do k1 = 1, nVertLevels
               do nt = 1, nTracers
                  zhi = min (layerInterfaceSigma_Input(k1+1), layerInterfaceSigma(k2+1))
                  zlo = max (layerInterfaceSigma_Input(k1), layerInterfaceSigma(k2))
                  hOverlap = max (zhi-zlo, 0.d0) * thisThk
                  hTsum(nt,k2) = htsum(nt,k2) + tracers(nt,k1,iCell) * hOverlap
               enddo         ! nt
            enddo            ! k1
         enddo               ! k2

         !-----------------------------------------------------------------
         ! Compute tracer values in new layers
         !-----------------------------------------------------------------

         do k = 1, nVertLevels
            do nt = 1, nTracers
               tracers(nt,k,iCell) = hTsum(nt,k) / layerThickness(k, iCell) * real(li_mask_is_ice_int(cellMask(iCell)), RKIND)
               layerThickness(k,iCell) = layerThickness(k,iCell) * real(li_mask_is_ice_int(cellMask(iCell)), RKIND)
               ! Note: the last term is used to multiply by 0 where there is no ice in order to clear out garbage introduced by eps
            enddo         ! nt
         enddo            ! k

      enddo ! nCells


      deallocate(layerInterfaceSigma_Input)
      deallocate(hTsum)

   end subroutine vertical_remap


!***********************************************************************
!
!  subroutine cells_to_vertices_1dfield_using_kiteAreas
!
!> \brief   Converts a 1d scalar field from cells to vertices
!> \author  Matt Hoffman
!> \date    21 May 2012
!> \details
!>  This routine converts a 1d scalar field from cells to vertices.
!>  It will give garbage values on obtuse triangles!  But it does work
!>  on periodic meshes.
!>  TODO: It would be more efficient to calculate the weights once on init and then only
!>  perform the interp. in this routine.
!-----------------------------------------------------------------------
   subroutine cells_to_vertices_1dfield_using_kiteAreas(meshPool, fieldCells, fieldVertices)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      real (kind=RKIND), dimension(:), intent(in) :: &
         fieldCells    !< Input: field on cells

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), intent(out) :: &
         fieldVertices    !< Input: field on vertices

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      integer, dimension(:,:), pointer :: cellsOnVertex
      integer, pointer :: nVertices, vertexDegree
      integer :: iCell, icell2, iVertex, cellIndex
      real (kind=RKIND) :: fVertexAccum, baryweight, weightAccum

      ! Get needed items from mesh pool
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

      ! Calculate h on vertices using barycentric interpolation
      do iVertex = 1, nVertices  ! Loop over vertices
        fVertexAccum = 0.0_RKIND
        weightAccum = 0.0_RKIND
        ! Loop over cells on this vertex
        do iCell = 1, vertexDegree
          cellIndex = cellsOnVertex(iCell, iVertex)
          baryweight = 0.0_RKIND
          do iCell2 = 1, vertexDegree
            if (iCell2 /= icell) baryweight = baryweight + 0.5 * kiteAreasOnVertex(iCell2, iVertex)
          enddo
          fVertexAccum = fVertexAccum + baryweight * fieldCells(cellIndex)  ! add the contribution from this cell's kite
          weightAccum = weightAccum + kiteAreasOnVertex(iCell, iVertex)  ! This doesn't match areaTriangle for obtuse triangles!!!
        enddo
        fieldVertices(iVertex) = fVertexAccum / weightAccum  ! I assume this should never be 0...
      enddo

   end subroutine cells_to_vertices_1dfield_using_kiteAreas


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


end module li_diagnostic_vars

