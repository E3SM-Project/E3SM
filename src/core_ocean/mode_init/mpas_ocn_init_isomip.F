! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_isomip
!
!> \brief MPAS ocean initialize case -- isomip
!> \author Xylar Asay-Davis
!> \date   06/01/2015
!> \details
!>  This module contains the routines for initializing the
!>  the Ice Shelf/Ocean Model Intercomparision Project (ISOMIP) test cases
!
!-----------------------------------------------------------------------

module ocn_init_isomip

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

   use ocn_init_ssh_and_landIcePressure

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

   public :: ocn_init_setup_isomip, &
             ocn_init_validate_isomip

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_isomip
!
!> \brief   Setup for ISoMIP test cases
!> \author  Xylar Asay-Davis
!> \date    06/01/2015
!> \details
!>  This routine sets up the initial conditions for the ISOMIP test cases.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_isomip(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      ! Define config variable pointers
      character (len=StrKIND), pointer :: config_init_configuration, &
        config_isomip_vertical_level_distribution

      real (kind=RKIND), pointer :: config_isomip_bottom_depth, &
        config_isomip_temperature, &
        config_isomip_salinity, &
        config_isomip_restoring_temperature, &
        config_isomip_restoring_salinity, &
        config_isomip_temperature_piston_velocity, &
        config_isomip_salinity_piston_velocity, &
        config_isomip_coriolis_parameter, &
        config_isomip_southern_boundary, &
        config_isomip_northern_boundary, &
        config_isomip_western_boundary, &
        config_isomip_eastern_boundary, &
        config_isomip_angle, &
        config_isomip_y1, &
        config_isomip_z1, &
        config_isomip_ice_fraction1, &
        config_isomip_y2, &
        config_isomip_z2, &
        config_isomip_ice_fraction2, &
        config_isomip_y3, &
        config_isomip_z3, &
        config_isomip_ice_fraction3, &
        config_isomip_effective_density

      logical, pointer :: on_a_sphere

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool

      type (block_type), pointer :: block_ptr

      ! Define dimension pointers
      integer, pointer :: nCells, nEdgesSolve, nVertLevels
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      integer, dimension(:), pointer :: maxLevelCell, modifySSHMask, landIceMask
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, &
                                                  vertCoordMovementWeights, bottomDepth, &
                                                  fCell, fEdge, fVertex, dcEdge, &
                                                  landIceFraction, landIceSurfaceTemperature, &
                                                  refLayerThickness, refZMid, &
                                                  ssh
      !real (kind=RKIND), dimension(:), pointer :: temperatureRestore, salinityRestore, maskRestore
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, zMid
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
      real (kind=RKIND), dimension(:, :), pointer ::    activeTracersPistonVelocity, activeTracersSurfaceRestoringValue

      integer :: iCell, k, iFit

      real(kind=RKIND) :: x, y, ySouth, yNorth, xWest, xEast, &
        pressure, dcEdgeMinGlobal, dcEdgeMin

      real(kind=RKIND), parameter :: eps = 1e-3_RKIND

      real(kind=RKIND), dimension(5) :: yFit, zFit, fracFit

      real(kind=RKIND), dimension(:), pointer :: columnThicknessFraction

      iErr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

      if(trim(config_init_configuration) .ne. trim('isomip')) return

      ! Setup configuration
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_vertical_level_distribution', config_isomip_vertical_level_distribution)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_bottom_depth', config_isomip_bottom_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_temperature', config_isomip_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_salinity', config_isomip_salinity)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_restoring_temperature', config_isomip_restoring_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_restoring_salinity', config_isomip_restoring_salinity)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_temperature_piston_velocity', config_isomip_temperature_piston_velocity)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_salinity_piston_velocity', config_isomip_salinity_piston_velocity)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_coriolis_parameter', config_isomip_coriolis_parameter)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_southern_boundary', config_isomip_southern_boundary)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_northern_boundary', config_isomip_northern_boundary)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_western_boundary', config_isomip_western_boundary)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_eastern_boundary', config_isomip_eastern_boundary)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_y1', config_isomip_y1)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_z1', config_isomip_z1)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_ice_fraction1', config_isomip_ice_fraction1)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_y2', config_isomip_y2)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_z2', config_isomip_z2)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_ice_fraction2', config_isomip_ice_fraction2)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_y3', config_isomip_y3)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_z3', config_isomip_z3)
      call mpas_pool_get_config(ocnConfigs, 'config_isomip_ice_fraction3', config_isomip_ice_fraction3)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      if ( on_a_sphere ) call mpas_log_write(&
           'The ISOMIP configuration can only be applied to a planar mesh. Exiting...', MPAS_LOG_CRIT)

      ySouth = config_isomip_southern_boundary
      yNorth = config_isomip_northern_boundary
      xWest = config_isomip_western_boundary
      xEast = config_isomip_eastern_boundary

      dcEdgeMin = 1.0E10_RKIND
      ! Determine local min and max values.
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
        dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))
        block_ptr => block_ptr % next
      end do

      call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

      yFit(2) = config_isomip_y1
      yFit(3) = config_isomip_y2
      yFit(4) = config_isomip_y3
      zFit(2) = config_isomip_z1
      zFit(3) = config_isomip_z2
      zFit(4) = config_isomip_z3

      yFit(1) = min(yFit(2),ySouth)-eps
      zFit(1) = zFit(2)
      yFit(5) = max(yFit(4),yNorth)+eps
      zFit(5) = zFit(4)

      fracFit(2) = config_isomip_ice_fraction1
      fracFit(3) = config_isomip_ice_fraction2
      fracFit(4) = config_isomip_ice_fraction3
      fracFit(1) = fracFit(2)
      fracFit(5) = fracFit(4)

      allocate(columnThicknessFraction(nVertLevels))
      if(trim(config_isomip_vertical_level_distribution) == "constant") then
         columnThicknessFraction(:) = 1.0_RKIND/nVertLevels
      else if(trim(config_isomip_vertical_level_distribution) == "boundary_layer") then
         if(mod(nVertLevels,2) == 0) then
            columnThicknessFraction(nVertLevels/2) = 0.25_RKIND
            columnThicknessFraction(nVertLevels/2+1) = 0.25_RKIND
         else
            columnThicknessFraction(nVertLevels/2) = 0.125_RKIND
            columnThicknessFraction(nVertLevels/2+1) = 0.5_RKIND
            columnThicknessFraction(nVertLevels/2+2) = 0.125_RKIND
         end if
         do k = nVertLevels/2-1, 2, -1
            columnThicknessFraction(k) = 0.5_RKIND*columnThicknessFraction(k+1)
            columnThicknessFraction(nVertLevels-k+1) = 0.5_RKIND*columnThicknessFraction(nvertLevels-k)
         end do
         columnThicknessFraction(1) = columnThicknessFraction(2)
         columnThicknessFraction(nVertLevels) = columnThicknessFraction(nVertLevels-1)
      end if

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
        call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)
        call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)
        call mpas_pool_get_array(forcingPool, 'landIceSurfaceTemperature', landIceSurfaceTemperature)
        call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

        call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
        call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
        call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)
        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)


        call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', &
            activeTracersPistonVelocity, 1)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', &
            activeTracersSurfaceRestoringValue, 1)


        ! flat bottom
        maxLevelCell(:) = nVertLevels
        maxLevelCell(nCells+1) = -1
        bottomDepth(:) =  abs(config_isomip_bottom_depth)
        do iCell = 1, nCells
          do k = 1, nVertLevels
             restingThickness(k, iCell) = columnThicknessFraction(k)*bottomDepth(iCell)
          end do
        end do

        call ocn_mark_north_boundary(meshPool, yNorth, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, ySouth, dcEdgeMinGlobal, iErr)
        call ocn_mark_east_boundary(meshPool, xEast, dcEdgeMinGlobal, iErr)
        call ocn_mark_west_boundary(meshPool, xWest, dcEdgeMinGlobal, iErr)

        ! Set refBottomDepth
        refBottomDepth(1) =  columnThicknessFraction(1)*abs(config_isomip_bottom_depth)
        do k = 2, nVertLevels
            refBottomDepth(k) = refBottomDepth(k-1) + columnThicknessFraction(k)*abs(config_isomip_bottom_depth)
        end do

        ! Compute refLayerThickness and refZMid
        call ocn_compute_layerThickness_zMid_from_bottomDepth(refLayerThickness,refZMid, &
             refBottomDepth,refBottomDepth(nVertLevels), &
             nVertLevels,nVertLevels,iErr)

        fCell(:) = config_isomip_coriolis_parameter
        fEdge(:) = config_isomip_coriolis_parameter
        fVertex(:) = config_isomip_coriolis_parameter

        landIceFraction(:) = 0.0_RKIND
        landIceSurfaceTemperature(:) = -25.0_RKIND !doesn't matter because ice is insulating
        modifySSHMask(:) = 0
        landIceMask(:) = 0


        do iCell = 1, nCells
          ! tracers computed using restingThickness with no ice shelf
          x =  xCell(iCell)
          y =  yCell(iCell)

          ! Compute iceDraft by linear fit
          do iFit = 1, 4
            if((y >= yFit(iFit)) .and. (y <= yFit(iFit+1))) then
              ssh(iCell) = (zFit(iFit+1) - zFit(iFit))*(y - yFit(iFit)) &
                         / (yFit(iFit+1) - yFit(iFit)) + zFit(iFit)
              landIceFraction(iCell) = (fracFit(iFit+1) - fracFit(iFit))*(y - yFit(iFit)) &
                         / (yFit(iFit+1) - yFit(iFit)) + fracFit(iFit)
              exit
            end if
          end do
          if(landIceFraction(iCell) > 0.5_RKIND) then
            landIceMask(iCell) = 1
          else
            landIceFraction(iCell) = 0.0_RKIND
          end if

          if(ssh(iCell) < 0.0_RKIND) then
            modifySSHMask(iCell) = 1
          end if

          if(.not. associated(activeTracers)) then
            call mpas_log_write( 'isomip test case needs activeTracers package to be active.', MPAS_LOG_CRIT)
            return
          end if

          activeTracers(index_temperature, :, iCell) = config_isomip_temperature
          activeTracers(index_salinity, :, iCell) = config_isomip_salinity
          if(associated(debugTracers)) then
            debugTracers(index_tracer1, :, iCell) = 1.0_RKIND
          end if

           ! Set surface temperature restoring value and rate
           ! Value in units of C, piston velocity in units of m/s
           if ( associated(activeTracersSurfaceRestoringValue) ) then
              activeTracersSurfaceRestoringValue(index_temperature, iCell) = config_isomip_restoring_temperature
           end if
           if ( associated(activeTracersPistonVelocity) ) then
              ! only restore where there *isn't* land ice
              activeTracersPistonVelocity(index_temperature, iCell) = (1.0_RKIND - landIceFraction(iCell)) &
                 * config_isomip_temperature_piston_velocity
           end if

           ! Set surface salinity restoring value and rate
           ! Value in units of PSU, piston velocity in units of m/s
           if ( associated(activeTracersSurfaceRestoringValue) ) then
              activeTracersSurfaceRestoringValue(index_salinity, iCell) = config_isomip_restoring_salinity
           end if
           if ( associated(activeTracersPistonVelocity) ) then
              ! only restore where there *isn't* land ice
              activeTracersPistonVelocity(index_salinity, iCell) = (1.0_RKIND - landIceFraction(iCell)) &
                 * config_isomip_salinity_piston_velocity
           end if

        end do

        block_ptr => block_ptr % next
      end do

      deallocate(columnThicknessFraction)


      ! compute the vertical grid (layerThickness, restingThickness, maxLevelCell, zMid)
      !  based on ssh, bottomDepth and refBottomDepth
      call ocn_init_ssh_and_landIcePressure_vertical_grid(domain, iErr)

      if(iErr .ne. 0) then
        call mpas_log_write( 'ocn_init_ssh_and_landIcePressure_vertical_grid failed.', MPAS_LOG_CRIT)
        return
      end if

      ! compute or update the land-ice pressure (or possibly SSH), also computing density along the way
      ! If this is the initial guess, the vertical grid and activeTracers may also be recomputed based on SSH
      call ocn_init_ssh_and_landIcePressure_balance(domain, iErr)

      if(iErr .ne. 0) then
         call mpas_log_write( 'ocn_init_ssh_and_landIcePressure_balance failed.', MPAS_LOG_CRIT)
         return
      end if

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

        call ocn_compute_Haney_number(domain, iErr)

        if(iErr .ne. 0) then
           call mpas_log_write( 'ocn_compute_Haney_number failed.', MPAS_LOG_CRIT)
           return
        end if

        block_ptr => block_ptr % next
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_isomip!}}}

!***********************************************************************
!
!  routine ocn_init_validate_isomip
!
!> \brief   Validation for ISOMIP test cases
!> \author  Xylar Asay-Davis
!> \date    06/01/2015
!> \details
!>  This routine validates the configuration options for the ISOMIP test cases.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_isomip(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool, packagePool

      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration, &
        config_isomip_vertical_level_distribution

      integer, pointer :: config_vert_levels, config_isomip_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(trim(config_init_configuration) .ne. trim('isomip')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_isomip_vert_levels', config_isomip_vert_levels)

      if(config_vert_levels <= 0 .and. config_isomip_vert_levels > 0) then
         config_vert_levels = config_isomip_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for isomip. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

      call mpas_pool_get_config(configPool, 'config_isomip_vertical_level_distribution', config_isomip_vertical_level_distribution)

      if((trim(config_isomip_vertical_level_distribution) .ne. "constant") &
         .and. (trim(config_isomip_vertical_level_distribution) .ne. "boundary_layer"))  then
         call mpas_log_write( 'Validation failed for isomip. Bad vertical level distribution.', MPAS_LOG_CRIT)
         iErr = 1
         return
      end if



   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_isomip!}}}


!***********************************************************************

end module ocn_init_isomip

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
