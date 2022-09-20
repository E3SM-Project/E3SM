! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_tidal_boundary
!
!> \brief MPAS ocean initialize case -- Tests tidal boundary
!> \author  Phillip Wolfram
!> \date   04/05/2019
!> \details
!>  This module contains the routines for initializing the
!>  the tidal_boundary test case
!
!-----------------------------------------------------------------------

module ocn_init_tidal_boundary

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
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

   public :: ocn_init_setup_tidal_boundary, &
             ocn_init_validate_tidal_boundary

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_tidal_boundary
!
!> \brief   Setup for Tidal Boundary on drying slope
!> \author  Phillip Wolfram, Xylar Asay-Davis
!> \date    04/05/2019
!> \details
!>  This routine sets up the initial conditions for the tidal_boundary test case.
!
!-----------------------------------------------------------------------
  subroutine ocn_init_setup_tidal_boundary(domain, iErr)!{{{

    !--------------------------------------------------------------------

    implicit none

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr
    real (kind=RKIND) :: yMin, yMax, dcEdgeMin
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, dcEdgeMinGlobal
    real (kind=RKIND) :: plugWidth
    real (kind=RKIND) :: slopeCenter, slopeWidth, depCenter, depWidth

    type (block_type), pointer :: block_ptr

    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: forcingPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: tracersPool

    integer :: iCell, k, N

    ! Define dimensions
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity, index_tracer1

    ! Define arrays
    integer, dimension(:), pointer :: minLevelCell, maxLevelCell
    integer, dimension(:), pointer :: vegetationMask
    real (kind=RKIND), dimension(:), pointer :: vegetationHeight,vegetationDensity, vegetationDiameter
    real (kind=RKIND), dimension(:), pointer :: yCell, refBottomDepth, bottomDepth, vertCoordMovementWeights, dcEdge
    real (kind=RKIND), dimension(:), pointer :: tidalInputMask
    real (kind=RKIND), dimension(:), pointer :: bottomDrag
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations
    real (kind=RKIND),  parameter :: eps=1.0e-12
    real (kind=RKIND) :: cff1,cff2,cff3, dep_mark1, dep_mark2 ! intermediate variables

    iErr = 0

    call mpas_log_write('here')
    if(config_init_configuration .ne. trim('tidal_boundary')) return

    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid(config_vertical_grid, interfaceLocations)

    ! Initalize y values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max y value.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      yMin = min( yMin, minval(yCell(1:nCellssolve)))
      yMax = max( yMax, maxval(yCell(1:nCellssolve)))
      dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

      block_ptr => block_ptr % next
    end do

    ! Determine global min and max y value. This is so the domain
    ! can be split into right and left.
    call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

    plugWidth = (yMaxGlobal - yMinGlobal) * config_tidal_boundary_plug_width_frac

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
      call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)

      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
      call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
      call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)

      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
      call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)


      if (config_tidal_boundary_right_bottom_depth < config_tidal_boundary_left_bottom_depth) then
        call mpas_log_write('Right boundary must be deeper than left boundary!', MPAS_LOG_CRIT)
      end if

      ! Set refBottomDepth, bottomDepth, and maxLevelCell
      if (min(config_tidal_boundary_left_bottom_depth, config_tidal_boundary_right_bottom_depth) < 0.0_RKIND) then
        ! consider the case where there is wetting / drying and vertical mesh resolution is needed "on land"
        do k = 1, nVertLevels
          refBottomDepth(k) = config_tidal_boundary_left_bottom_depth + &
            (config_tidal_boundary_right_bottom_depth - config_tidal_boundary_left_bottom_depth)* interfaceLocations(k+1)
        end do
      else
        ! assumes we just need to build vertical mesh to deepest point (e.g., no "on land' vertical mesh)
        do k = 1, nVertLevels
          refBottomDepth(k) = &
          max(config_tidal_boundary_left_bottom_depth, config_tidal_boundary_right_bottom_depth) * interfaceLocations(k+1)
        end do
      end if

      if (config_tidal_boundary_use_distances) then
        yMin = config_tidal_boundary_left_value
        yMax = config_tidal_boundary_right_value
      end if

      do iCell = 1, nCellsSolve
        bottomDepth(iCell) = config_tidal_boundary_left_bottom_depth &
          + (yCell(iCell) - yMin) / (yMax - yMin)  * &
          (config_tidal_boundary_right_bottom_depth - config_tidal_boundary_left_bottom_depth)
      end do

      if (config_use_variable_drag) then
        call mpas_pool_get_array(forcingPool, 'bottomDrag', bottomDrag)
        do iCell = 1, nCellsSolve
          bottomDrag(iCell) = config_tidal_forcing_left_Cd_or_n &
            + (yCell(iCell) - yMin) / (yMax - yMin)  * &
            (config_tidal_forcing_right_Cd_or_n - config_tidal_forcing_left_Cd_or_n)
        end do
      end if

      if (config_use_vegetation_drag) then
        call mpas_pool_get_array(forcingPool, 'vegetationMask', vegetationMask)
        vegetationMask = 0
        call mpas_pool_get_array(forcingPool, 'vegetationDiameter', vegetationDiameter)
        vegetationDiameter = config_idealized_vegetation_diameter
        call mpas_pool_get_array(forcingPool, 'vegetationHeight', vegetationHeight)
        vegetationHeight = config_idealized_vegetation_height
        call mpas_pool_get_array(forcingPool, 'vegetationDensity', vegetationDensity)
        vegetationDensity = config_idealized_vegetation_density
      endif

      if (config_use_idealized_transect) then
        call mpas_pool_get_array(forcingPool, 'bottomDrag', bottomDrag)
        config_idealized_transect_Lmarsh = 1.0_RKIND - config_idealized_transect_Lcoast &
                                                    - config_idealized_transect_Lshore
        if (config_idealized_transect_Lmarsh .lt. 0.0_RKIND) then
          call mpas_log_write("Lshore+Lcoast cannot be bigger than 1.0")
          iErr = 1
        endif
        cff1 = yMax*config_idealized_transect_Lcoast
        cff2 = yMax*config_idealized_transect_Lmarsh
        cff3 = yMax*config_idealized_transect_Lshore
        ! by defining the slopes and left_bottom_depth, the pre-defined right_bottom_depth won't work.
        ! Need redefine it.
        config_tidal_boundary_right_bottom_depth = config_tidal_boundary_left_bottom_depth + &
                                                   cff1*config_idealized_transect_Scoast + &
                                                   cff2*config_idealized_transect_Smarsh + &
                                                   cff3*config_idealized_transect_Sshore
        do iCell = 1, nCellsSolve
          if (yCell(iCell) .lt. cff1) then
            bottomDepth(iCell) = config_tidal_boundary_left_bottom_depth &
              + (yCell(iCell)-yMin)*config_idealized_transect_Scoast
            dep_mark1 = bottomDepth(iCell)
            bottomDrag(iCell) = config_idealized_transect_roughness
          elseif (yCell(iCell) .lt. (cff1+cff2)) then
            bottomDepth(iCell) = dep_mark1 + (yCell(iCell)-cff1)*config_idealized_transect_Smarsh
            dep_mark2 = bottomDepth(iCell)
            bottomDrag(iCell) = config_idealized_transect_roughness_marsh
            if (config_use_vegetation_drag) vegetationMask(iCell) = 1
          else
            if (cff2 .eq. 0.0_RKIND) dep_mark2 = dep_mark1
            bottomDepth(iCell) = dep_mark2 + (yCell(iCell)-cff1-cff2)*config_idealized_transect_Sshore
            bottomDrag(iCell) = config_idealized_transect_roughness
          endif
        end do !! do iCell
      else
      ! if not config_idealized_transect, assign a constant slope
        do iCell = 1, nCellsSolve
          bottomDepth(iCell) = config_tidal_boundary_left_bottom_depth &
            + (yCell(iCell) - yMin) / (yMax - yMin)  * &
            (config_tidal_boundary_right_bottom_depth - config_tidal_boundary_left_bottom_depth)
        end do
      end if   !! if config_idealized_transect

      if (config_tidal_boundary_right_bottom_depth < config_tidal_boundary_left_bottom_depth) then
        call mpas_log_write('Right boundary must be deeper than left boundary!', MPAS_LOG_CRIT)
      end if

      ! Set refBottomDepth, bottomDepth, and maxLevelCell
      do k = 1, nVertLevels
        refBottomDepth(k) = config_tidal_boundary_right_bottom_depth * interfaceLocations(k+1)
      end do

      if (config_use_wetting_drying .and. config_tidal_start_dry .and. &
          trim(config_tidal_boundary_layer_type) == 'zstar') then
        do iCell = 1, nCellsSolve
          ! make sure bottomDepth not smaller than 3 layers)
          bottomDepth(iCell) = max(refBottomDepth(2), bottomDepth(iCell))
        end do
      end if

      if (.not. config_alter_ICs_for_pcs .and. trim(config_tidal_boundary_layer_type) == 'zstar') then
        do iCell = 1, nCellsSolve
          do k = 1,nVertLevels
            if (refBottomDepth(k) > bottomDepth(iCell)) then
              bottomDepth(iCell) = refBottomDepth(k)
              exit
            end if
          end do
        end do
      end if

      minLevelCell(:) = 1

      if(trim(config_tidal_boundary_layer_type) == 'zstar') then

        ! Determine maxLevelCell based on bottomDepth and refBottomDepth
        ! Also set botomDepth based on refBottomDepth, since
        ! above bottomDepth was set with continuous analytical functions,
        ! and needs to be discrete
        do iCell = 1, nCellsSolve
        maxLevelCell(iCell) = nVertLevels
          if (nVertLevels > 1) then
            do k = 1, nVertLevels
              if (bottomDepth(iCell) < refBottomDepth(k)) then
                maxLevelCell(iCell) = k-1
                bottomDepth(iCell) = refBottomDepth(k-1)
                exit
              end if
            end do
          end if

          ssh(iCell) = 0.0_RKIND
        end do

      else if (trim(config_tidal_boundary_layer_type) == 'sigma') then

        do iCell = 1, nCellsSolve
          maxLevelCell(iCell) = nVertLevels
        end do

      else if (trim(config_tidal_boundary_layer_type) == 'hybrid') then

        do iCell = 1, nCellsSolve
          maxLevelCell(iCell) = nVertLevels
          if (nVertLevels > 1) then
            do k = 1, nVertLevels
              if (bottomDepth(iCell) < refBottomDepth(k)) then
                if (k < config_tidal_boundary_min_vert_levels) then
                  ! Use sigma
                  maxLevelCell(iCell) = config_tidal_boundary_min_vert_levels
                  exit
                else
                  ! Use zstar
                  maxLevelCell(iCell) = k-1
                  bottomDepth(iCell) = refBottomDepth(k-1)
                  exit
                end if
              end if
            end do
          end if

          ssh(iCell) = 0.0_RKIND
        end do

      end if


      if (config_use_wetting_drying .and. config_tidal_start_dry) then
        do iCell = 1, nCellsSolve
          ssh(iCell) = -bottomDepth(iCell) + config_drying_min_cell_height*maxLevelCell(iCell)
          ssh(iCell) = MAX(ssh(iCell),-config_tidal_forcing_monochromatic_baseline)
          ! also computes zMid
          do k = 1, maxLevelCell(iCell)
            layerThickness(k,iCell) = (ssh(iCell) + bottomDepth(iCell))/maxLevelCell(iCell)
            restingThickness(k, iCell) = bottomDepth(iCell)/maxLevelCell(iCell)
          end do
        end do
      else
        do iCell = 1, nCellsSolve
          ! restingThickness is z-level, with ssh = 0
          call ocn_compute_z_level_layerThickness(restingThickness(:,iCell), refBottomDepth, 0.0_RKIND, &
                                                  bottomDepth(iCell), minLevelCell(iCell),              &
                                                  maxLevelCell(iCell), nVertLevels, iErr)

          ! stretch restingThickness to get the z-star layerThickness
          call ocn_compute_z_star_layerThickness(layerThickness(:,iCell), restingThickness(:,iCell),  &
                                                 ssh(iCell), bottomDepth(iCell), minLevelCell(iCell), &
                                                 maxLevelCell(iCell), nVertLevels, iErr)

          ! compute zMid
          call ocn_compute_zMid_from_layerThickness(zMid(:,iCell), layerThickness(:,iCell), ssh(iCell), &
                                                    minLevelCell(iCell), maxLevelCell(iCell), nVertLevels, iErr)
          do k = 1, maxLevelCell(iCell)
            restingThickness(k, iCell) = layerThickness(k, iCell)
          end do
        end do
      end if

      ! Set vertCoordMovementWeights
      vertCoordMovementWeights(:) = 1.0_RKIND

      ! Set tidal boundary mask
      do iCell = 1, nCellsSolve
        tidalInputMask(iCell) = 0.0_RKIND
        if (yCell(iCell) > (yMax - dcEdgeMinGlobal/2.0_RKIND)) then
          tidalInputMask(iCell) = 1.0_RKIND
        ! spread it over multiple cells
        !if (yCell(iCell) > (25.0e3 - 3*dcEdgeMinGlobal)) then
        !  tidalInputMask(iCell) = exp(-((yCell(iCell)-25.0e3)/dcEdgeMinGlobal)**2.0)
        end if
      end do

      ! check that there is some tidalInputMask 
      if (.not. sum(tidalInputMask) > 0) then
        call mpas_log_write('Input mask for tidal case is not set!', MPAS_LOG_CRIT)
      end if

      ! Set salinity
      if ( associated(activeTracers) ) then
        do iCell = 1, nCellsSolve
          activeTracers(index_salinity, :, iCell) = config_tidal_boundary_salinity
        end do
      end if

      ! Set temperature
      if ( associated(activeTracers) ) then
        do iCell = 1, nCellsSolve
          do k = 1, maxLevelCell(iCell)
            if(yCell(iCell) < yMinGlobal + plugWidth) then
              activeTracers(index_temperature, k, iCell) = config_tidal_boundary_plug_temperature
            else
              activeTracers(index_temperature, k, iCell) = config_tidal_boundary_domain_temperature
            end if
          end do
        end do
      end if

      ! Set debug tracer to match plug onshore
      if ( associated(debugTracers) ) then
        do iCell = 1, nCellsSolve
          do k = 1, maxLevelCell(iCell)
            if(yCell(iCell) < (yMinGlobal + plugWidth)) then
              debugTracers(index_tracer1, k, iCell) = 1.0_RKIND
            else
              debugTracers(index_tracer1, k, iCell) = 0.0_RKIND
            end if
          end do
        end do
      end if

      block_ptr => block_ptr % next
    end do

    deallocate(interfaceLocations)

    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_tidal_boundary!}}}


!***********************************************************************
!
!  routine ocn_init_validate_tidal_boundary
!
!> \brief   Validation for Wetting and drying by tidal boundary test case
!> \author  Phillip Wolfram
!> \date    04/05/2019
!> \details
!>  This routine validates the configuration options for the tidal_boundary test case.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_validate_tidal_boundary(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_tidal_boundary_vert_levels, config_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('tidal_boundary')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_vert_levels', config_tidal_boundary_vert_levels)

      if(config_vert_levels <= 0 .and. config_tidal_boundary_vert_levels > 0) then
         config_vert_levels = config_tidal_boundary_vert_levels
         !call mpas_log_write( 'Using value of $i', intArgs=(/ config_vert_levels /))
      else if(config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for tidal_boundary test case. ' &
           // 'Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_tidal_boundary!}}}


!***********************************************************************

end module ocn_init_tidal_boundary

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
