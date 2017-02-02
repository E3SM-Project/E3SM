! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_soma
!
!> \brief MPAS ocean initialize case -- Simulating Ocean Mesoscale Activity (SOMA)
!> \author Todd Ringler
!> \date   10/08/2013
!> \details
!>  This module contains the routines for initializing the
!>  the idealized SOMA test case
!
!-----------------------------------------------------------------------

module ocn_init_soma

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_constants
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

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

   public :: ocn_init_setup_soma, &
             ocn_init_validate_soma

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_soma
!
!> \brief   Setup for soma test case
!> \author  Todd Ringler
!> \date    02/26/2014
!> \details
!>  This routine sets up the initial conditions for the
!>  SOMA configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_soma(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      ! local work variables
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool, forcingPool, tracersPool
      type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool, tracersInteriorRestoringFieldsPool

      integer :: iCell, iEdge, iVertex, k, idx
      real (kind=RKIND) :: distance, deltaLon, deltaLat, xDistance, yDistance, zMid, sphereRadius
      real (kind=RKIND) :: lonCurrent, latCurrent
      real (kind=RKIND) :: deltay, depth, factor, latCenter, lonCenter, windStress
      real (kind=RKIND) :: temperature, salinity
      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      ! Define config variable pointers
      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid

      ! SOMA test case run-time configuration parameters
      integer, pointer :: config_soma_vert_levels
      real (kind=RKIND), pointer :: config_eos_linear_alpha
      real (kind=RKIND), pointer :: config_soma_surface_salinity
      real (kind=RKIND), pointer :: config_soma_surface_temperature
      real (kind=RKIND), pointer :: config_soma_density_difference_linear
      real (kind=RKIND), pointer :: config_soma_thermocline_depth
      real (kind=RKIND), pointer :: config_soma_center_latitude
      real (kind=RKIND), pointer :: config_soma_center_longitude
      real (kind=RKIND), pointer :: config_soma_domain_width
      real (kind=RKIND), pointer :: config_soma_shelf_width
      real (kind=RKIND), pointer :: config_soma_shelf_depth
      real (kind=RKIND), pointer :: config_soma_bottom_depth
      real (kind=RKIND), pointer :: config_soma_phi
      real (kind=RKIND), pointer :: config_soma_ref_density
      real (kind=RKIND), pointer :: config_soma_density_difference
      real (kind=RKIND), pointer :: config_soma_surface_temp_restoring_at_center_latitude
      real (kind=RKIND), pointer :: config_soma_surface_temp_restoring_latitude_gradient
      real (kind=RKIND), pointer :: config_soma_restoring_temp_piston_vel

      ! Define dimension pointers
      integer, pointer :: nVertLevels, nCells, nVertLevelsP1, nCellsSolve, nEdgesSolve, nVerticesSolve
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      logical, pointer :: on_a_sphere, config_soma_use_surface_temp_restoring
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, bottomCell, refZMid, fCell, fEdge, fVertex
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND), dimension(:), pointer :: lonCell, latCell, latEdge, latVertex, bottomDepth
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
      real (kind=RKIND), dimension(:, :), pointer ::    activeTracersPistonVelocity, activeTracersSurfaceRestoringValue
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracersInteriorRestoringValue, activeTracersInteriorRestoringRate
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional

      ! assume no error
      iErr = 0

      ! test if SOMA is the desired configuration
      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('soma')) return

      ! get config variables
      call mpas_pool_get_config(domain % configs, 'config_eos_linear_alpha', config_eos_linear_alpha)
      call mpas_pool_get_config(domain % configs, 'config_soma_density_difference_linear', config_soma_density_difference_linear)
      call mpas_pool_get_config(domain % configs, 'config_soma_thermocline_depth', config_soma_thermocline_depth)
      call mpas_pool_get_config(domain % configs, 'config_soma_surface_temperature', config_soma_surface_temperature)
      call mpas_pool_get_config(domain % configs, 'config_soma_surface_salinity', config_soma_surface_salinity)
      call mpas_pool_get_config(domain % configs, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_get_config(domain % configs, 'config_soma_vert_levels', config_soma_vert_levels)
      call mpas_pool_get_config(domain % configs, 'config_soma_center_latitude', config_soma_center_latitude)
      call mpas_pool_get_config(domain % configs, 'config_soma_center_longitude', config_soma_center_longitude)
      call mpas_pool_get_config(domain % configs, 'config_soma_domain_width', config_soma_domain_width)
      call mpas_pool_get_config(domain % configs, 'config_soma_shelf_width', config_soma_shelf_width)
      call mpas_pool_get_config(domain % configs, 'config_soma_shelf_depth', config_soma_shelf_depth)
      call mpas_pool_get_config(domain % configs, 'config_soma_bottom_depth', config_soma_bottom_depth)
      call mpas_pool_get_config(domain % configs, 'config_soma_phi', config_soma_phi)
      call mpas_pool_get_config(domain % configs, 'config_soma_ref_density', config_soma_ref_density)
      call mpas_pool_get_config(domain % configs, 'config_soma_density_difference', config_soma_density_difference)
      call mpas_pool_get_config(domain % configs, 'config_soma_surface_temp_restoring_at_center_latitude', &
         config_soma_surface_temp_restoring_at_center_latitude)
      call mpas_pool_get_config(domain % configs, 'config_soma_surface_temp_restoring_latitude_gradient', &
         config_soma_surface_temp_restoring_latitude_gradient)
      call mpas_pool_get_config(domain % configs, 'config_soma_restoring_temp_piston_vel', config_soma_restoring_temp_piston_vel)
      call mpas_pool_get_config(domain % configs, 'config_soma_use_surface_temp_restoring', config_soma_use_surface_temp_restoring)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
      sphereRadius = sphere_radius

      ! error checking
      if(.not. on_a_sphere) then
        write(stderrUnit, *) 'ERROR: SOMA test case can only be defined on a spherical mesh.'
        iErr = 1
        return
      else
       write(stderrUnit, *) 'SOMA test case using spherical radius of size: ', sphereRadius
      end if

      ! assign config variables
      nVertLevels  = config_soma_vert_levels
      nVertLevelsP1 = nVertLevels + 1

      ! Define interface locations
      allocate( interfaceLocations( nVertLevelsP1 ) )
      call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

      ! set center of SOMA domain
      ! Convert center locations to radians from degrees
      latCenter = config_soma_center_latitude * pii / 180.0_RKIND
      lonCenter = config_soma_center_longitude * pii / 180.0_RKIND

      ! Setup the vertical grid and layerThickness initial condition
      write(stderrUnit,*) 'setting up vertical grid and layer thickness'
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
        call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)

        ! Set layerThickness and restingThickness
        ! Uniform layer thickness across lat/lon
        do k = 1, nVertLevels
           layerThickness(k, :) = config_soma_bottom_depth * ( interfaceLocations(k+1) - interfaceLocations(k) )
           restingThickness(k, :) = layerThickness(k, :)
        end do

         ! Set refBottomDepth
         do k = 1, nVertLevels
            refBottomDepth(k) = config_soma_bottom_depth * interfaceLocations(k+1)
            refZMid(k) = -config_soma_bottom_depth * (interfaceLocations(k)+interfaceLocations(k+1))/2.0_RKIND
         end do

        block_ptr => block_ptr % next

      end do

      ! Set bathymetry
      write(stderrUnit,*) 'setting up bathymetry'
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

        ! set bottomDepth
        bottomDepth(:) = 0.0_RKIND
        do iCell = 1, nCells
           lonCurrent = lonCell(iCell)
           latCurrent = latCell(iCell)

           deltaLon = abs(lonCurrent - lonCenter)
           if (deltaLon .gt. pii) deltaLon = deltaLon - 2.0_RKIND*pii
           deltaLat = latCurrent - latCenter
           xDistance = deltaLon * sphereRadius * cos(latCurrent)
           yDistance = deltaLat * sphereRadius
           distance = sqrt( xDistance**2 + yDistance**2 )
           factor = 1.0 - distance**2 / config_soma_domain_width**2

           if(factor > config_soma_shelf_width) then
              bottomDepth(iCell) = config_soma_shelf_depth + (config_soma_bottom_depth-config_soma_shelf_depth)/2.0_RKIND &
                                 * (1.0 + tanh(factor/config_soma_phi))
           else
              bottomDepth(iCell) = -1.0_RKIND
           endif

           ! Set maxLevelCell to -1 for cells to be culled
           if (bottomDepth(iCell) > 0.0_RKIND) then
             maxLevelCell(iCell) = 1
           else
             maxLevelCell(iCell) = -1
           endif

           ! Determine maxLevelCell based on bottomDepth and refBottomDepth
           ! Also set botomDepth based on refBottomDepth, since
           ! above bottomDepth was set with continuous analytical functions,
           ! and needs to be discrete
           if (maxLevelCell(iCell) > 0) then
             maxLevelCell(iCell) = nVertLevels
             if (nVertLevels .gt. 1) then
                do k = 1, nVertLevels
                  if (bottomDepth(iCell) < refBottomDepth(k) ) then
                      maxLevelCell(iCell) = k-1
                      bottomDepth(iCell) = refBottomDepth(k-1)
                      exit
                  end if
                end do
             end if
           end if

        enddo ! Looping through with iCell

        block_ptr => block_ptr % next

      enddo !  done setting bathymetry

      ! mark cells for culling
       block_ptr => domain % blocklist
       do while (associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call ocn_mark_maxlevelcell(meshPool, iErr)
          block_ptr => block_ptr % next
       end do

      ! Set forcing boundary conditions and initial conditions for temperature and salinity
      write(stderrUnit,*) 'setting up forcing and initial T/S conditions'
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
        call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
        call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional, 1)

        call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', activeTracersPistonVelocity, 1)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, &
          'activeTracersSurfaceRestoringValue', activeTracersSurfaceRestoringValue, 1)

        do iCell = 1, nCells
           lonCurrent = lonCell(iCell)
           latCurrent = latCell(iCell)

           ! Set initial temperature and salinity
           do k = 1, nVertLevels
              zMid = refZMid(k)

              distance = config_soma_ref_density &
                  - (1.0_RKIND - config_soma_density_difference_linear) * config_soma_density_difference &
                  * tanh(zMid / config_soma_thermocline_depth) &
                  - config_soma_density_difference_linear * config_soma_density_difference &
                  * zMid / config_soma_bottom_depth
              factor = (config_soma_ref_density - distance) / config_eos_linear_alpha
              temperature = config_soma_surface_temperature + factor

              factor = - zMid / config_soma_bottom_depth
              salinity = config_soma_surface_salinity + 2.0_RKIND * factor

              if ( associated(activeTracers) ) then
                 activeTracers(index_temperature, k, iCell) = temperature
                 activeTracers(index_salinity, k, iCell) = salinity
              end if
           enddo

           if (config_soma_use_surface_temp_restoring) then
             ! surface restoring
             idx = index_temperature
             activeTracersSurfaceRestoringValue(idx,iCell) = &
               config_soma_surface_temp_restoring_at_center_latitude &
               + config_soma_surface_temp_restoring_latitude_gradient &
               *(latCurrent*180.0_RKIND/pii - config_soma_center_latitude)
             activeTracersPistonVelocity(idx,iCell) = config_soma_restoring_temp_piston_vel

             idx = index_salinity
             activeTracersSurfaceRestoringValue(idx,iCell) = 34.0_RKIND
             activeTracersPistonVelocity(idx,iCell) = 0.0_RKIND
           end if


           ! Set up debugging tracers
           if ( associated(debugTracers) ) then
              debugTracers(index_tracer1, :, iCell) = 1.0_RKIND
           end if

        end do ! iCell = 1, nCells

        ! Set wind stress
        do iCell = 1, nCells
           lonCurrent = lonCell(iCell)
           latCurrent = latCell(iCell)

           deltay =  sphereRadius * ( latCurrent - latCenter )
           factor = 1.0_RKIND - 0.5_RKIND * deltay / config_soma_domain_width
           windstress = factor * 0.1_RKIND * exp( -(deltay / config_soma_domain_width)**2 ) &
                 * cos(pii * deltay / config_soma_domain_width)

           windStressZonal(iCell) = windStress
           windStressMeridional(iCell) = 0.0_RKIND

        end do

       ! Set Coriolis parameters
       do iCell = 1, nCellsSolve
          fCell(iCell) = 2.0_RKIND * omega * sin(latCell(iCell))
       end do
       do iEdge = 1, nEdgesSolve
          fEdge(iEdge) = 2.0_RKIND * omega * sin(latEdge(iEdge))
       end do
       do iVertex = 1, nVerticesSolve
          fVertex(iVertex) = 2.0_RKIND * omega * sin(latVertex(iVertex))
       end do

        block_ptr => block_ptr % next
      end do

      write(stderrUnit,*) 'exiting ocn_init_setup_soma'

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_soma!}}}

!***********************************************************************
!
!  routine ocn_init_validate_soma
!
!> \brief   Validation for SOMA test case
!> \author  Todd Ringler
!> \date    02/26/2014
!> \details
!>  This routine validates the configuration options for the SOMA test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_soma(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_soma_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('soma')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_soma_vert_levels', config_soma_vert_levels)

      if(config_vert_levels <= 0 .and. config_soma_vert_levels > 0) then
         config_vert_levels = config_soma_vert_levels
      else if (config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for SOMA. Not given a usable value for vertical levels.'
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_soma!}}}

end module ocn_init_soma

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
