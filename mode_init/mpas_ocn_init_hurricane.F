! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_hurricane
!
!> \brief MPAS ocean initialize case -- hurricane
!> \author Steven Brus
!> \date   01/08/19
!> \details
!>  This module contains the routines for initializing a hurricane.
!>
!-----------------------------------------------------------------------

module ocn_init_hurricane

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

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

   public :: ocn_init_setup_hurricane, &
             ocn_init_validate_hurricane

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_hurricane
!
!> \brief   Setup for this initial condition
!> \author  Steven Brus, Phillip J. Wolfram
!> \date    09/04/2019
!> \details
!>  This routine sets up the initial conditions for this case,
!>  including Gaussian hump initial condition..
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_hurricane(domain, iErr)!{{{

    !--------------------------------------------------------------------
    implicit none

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: diagnosticsPool

    ! local variables
    integer :: iCell, k, idx, iEdge, iVertex
    real (kind=RKIND) :: dlon, dlat, dw
    real (kind=RKIND) :: maxBottomDepth, globalMaxBottomDepth, globalMinBottomDepth
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations
    real (kind=RKIND), parameter :: eps=1.0e-12
    real (kind=RKIND), parameter :: pi = 4.0_RKIND*ATAN(1.0_RKIND)

    ! Define config variable pointers
    character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
    integer, pointer :: config_hurricane_vert_levels
    real (kind=RKIND), pointer :: config_hurricane_min_depth
    real (kind=RKIND), pointer :: config_hurricane_max_depth
    logical, pointer :: config_hurricane_use_gaussian_hump
    real (kind=RKIND), pointer :: config_hurricane_gaussian_lon_center
    real (kind=RKIND), pointer :: config_hurricane_gaussian_lat_center
    real (kind=RKIND), pointer :: config_hurricane_gaussian_width
    real (kind=RKIND), pointer :: config_hurricane_gaussian_hump_amplitude
    real (kind=RKIND), pointer :: config_hurricane_gaussian_slr_amp
    logical, pointer :: config_use_wetting_drying
    real (kind=RKIND), pointer :: config_drying_min_cell_height

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, bottomDepthObserved, &
         fCell, fEdge, fVertex, &
         latCell, latEdge, latVertex, &
         lonCell
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

    iErr = 0

    call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

    if(config_init_configuration .ne. trim('hurricane')) return

    ! Get config flag settings
    call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_vert_levels', config_hurricane_vert_levels)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_min_depth', config_hurricane_min_depth)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_max_depth', config_hurricane_max_depth)

    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_use_gaussian_hump', config_hurricane_use_gaussian_hump)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_gaussian_lon_center', config_hurricane_gaussian_lon_center)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_gaussian_lat_center', config_hurricane_gaussian_lat_center)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_gaussian_width', config_hurricane_gaussian_width)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_gaussian_hump_amplitude', config_hurricane_gaussian_hump_amplitude)
    call mpas_pool_get_config(ocnConfigs, 'config_hurricane_gaussian_slr_amp', config_hurricane_gaussian_slr_amp)

    ! wetting and drying parameters
    call mpas_pool_get_config(ocnConfigs, 'config_drying_min_cell_height', config_drying_min_cell_height)
    call mpas_pool_get_config(ocnConfigs, 'config_use_wetting_drying', config_use_wetting_drying)

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! you may restrict your case geometry as follows:
    if ( .not. on_a_sphere ) call mpas_log_write('The hurricane configuration can only be applied ' &
      // 'to a spherical mesh. Exiting...', MPAS_LOG_CRIT)

    nVertLevels  = config_hurricane_vert_levels
    nVertLevelsP1 = nVertLevels + 1


    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations, ocnConfigs )

    ! Find max bottom depth
    maxBottomDepth = 9e10
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'bottomDepthObserved', bottomDepthObserved)

      if (MINVAL(bottomDepthObserved) < maxBottomDepth) then
        maxBottomDepth = MINVAL(bottomDepthObserved)
      end if

      block_ptr => block_ptr % next
    end do
    maxBottomDepth = -1.0_RKIND*maxBottomDepth

    ! Enforce max bottom depth from namelist option
    call mpas_dmpar_max_real(domain % dminfo, maxBottomDepth, globalMaxBottomDepth)
    if (config_hurricane_max_depth < globalMaxBottomDepth) then
       globalMaxBottomDepth = config_hurricane_max_depth
    end if

    ! Find min bottom depth
    maxBottomDepth = -9e10
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'bottomDepthObserved', bottomDepthObserved)

      if (MAXVAL(bottomDepthObserved) > maxBottomDepth) then
        maxBottomDepth = MAXVAL(bottomDepthObserved)
      end if

      block_ptr => block_ptr % next
    end do
    maxBottomDepth = -1.0_RKIND*maxBottomDepth

    call mpas_dmpar_max_real(domain % dminfo, maxBottomDepth, globalMinBottomDepth)
    !write(*,*) globalMaxBottomDepth, globalMinBottomDepth

    !--------------------------------------------------------------------
    ! Use this section to set initial values
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

       call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

       call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'bottomDepthObserved', bottomDepthObserved)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
       call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
       call mpas_pool_get_array(meshPool, 'fCell', fCell)
       call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

       call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
       call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
       call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
       call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

       do k = 1, nVertLevels
          refBottomDepth(k) = globalMaxBottomDepth * interfaceLocations(k+1)
          refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * globalMaxBottomDepth
       end do

       ! Set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND


       ! Set ssh
       if (config_hurricane_use_gaussian_hump) then
         do iCell = 1, nCellsSolve
           dlat = latCell(iCell) - config_hurricane_gaussian_lat_center*pii/180.0_RKIND
           dlon = abs(lonCell(iCell) - config_hurricane_gaussian_lon_center*pii/180.0_RKIND)
           dw = config_hurricane_gaussian_width*pii/180.0_RKIND
           if (dlon > pii) dlon = dlon - 2.0_RKIND**pii
           ssh(iCell) = config_hurricane_gaussian_slr_amp + &
             config_hurricane_gaussian_hump_amplitude * exp(-(dlat**2.0_RKIND + dlon**2.0_RKIND) / &
             (2.0_RKIND*dw**2.0_RKIND))
         end do
       else
         do iCell = 1, nCellsSolve
           ssh(iCell) = 0.0_RKIND
         end do
       end if


       do iCell = 1, nCellsSolve

          bottomDepth(iCell) = -1.0_RKIND*bottomDepthObserved(iCell)

          ! enforce minimum and maximum bottom depth
          bottomDepth(iCell) = max(bottomDepth(iCell), config_hurricane_min_depth)
          bottomDepth(iCell) = min(bottomDepth(iCell), globalMaxBottomDepth)

          !! Flatten out Bermuda
          !if ((lonCell(iCell) > (-66.0_RKIND+360.0_RKIND)*pi/180.0_RKIND .and. &
          !     lonCell(iCell) < (-62.0_RKIND+360.0_RKIND)*pi/180.0_RKIND)  .and. &
          !    (latCell(iCell) >  30.0_RKIND*pi/180.0_RKIND .and. &
          !     latCell(iCell) <  34.0_RKIND*pi/180.0_RKIND)) then
          !
          !  bottomDepth(iCell) = globalMaxBottomDepth
          !end if
        end do

        if (config_use_wetting_drying) then
          do iCell = 1, nCellsSolve
            ! Set up vertical grid
            maxLevelCell(iCell) = nVertLevels
            ! ensures only 3 layers can be obtained on mesh
            if (bottomDepth(iCell) < 0 .or. bottomDepth(iCell) < refBottomDepth(3)) then
              maxLevelCell(iCell) = 3
            else
              ! snap to zstar coordinate
              do k = 1, nVertLevels
                if (bottomDepth(iCell) < refBottomDepth(k)) then
                  maxLevelCell(iCell) = k-1
                  bottomDepth(iCell) = refBottomDepth(k-1)
                  exit
                end if
              end do
            end if
          end do

          ! shift coordinates to accomodate wetting and drying
          do k = 1, nVertLevels
            refBottomDepth(k) = refBottomDepth(k) - globalMinBottomDepth
          end do

          do iCell = 1, nCellsSolve

            bottomDepth(iCell) = bottomDepth(iCell) - globalMinBottomDepth
            ssh(iCell) = ssh(iCell) + globalMinBottomDepth

            ! make sure depth is thick enough via ssh = TOTAL_DEPTH - bottomDepth
            ssh(iCell) = - bottomDepth(iCell) + &
              max(ssh(iCell) + bottomDepth(iCell), &
              maxLevelCell(iCell)*(config_drying_min_cell_height + eps))

            do k = 1, maxLevelCell(iCell)
              layerThickness(k,iCell) = max(config_drying_min_cell_height + eps, &
                (ssh(iCell) + bottomDepth(iCell))/maxLevelCell(iCell))
              if (layerThickness(k,iCell) < config_drying_min_cell_height) then
                call mpas_log_write('layerThickness($i,$i)=$r', MPAS_LOG_CRIT, &
                  intArgs=(/k,iCell/), &
                  realArgs=(/layerThickness(k,iCell)/))
              end if
            end do

            do k = 1, maxLevelCell(iCell)
              restingThickness(k,iCell) = bottomDepth(iCell)/maxLevelCell(iCell)
            end do
          end do
        else
          ! case without wetting and drying
          do iCell = 1, nCellsSolve
            ! Set up vertical grid
            maxLevelCell(iCell) = nVertLevels
            if (nVertLevels > 10) then
              ! z-star
              do k = 1, nVertLevels
                if (bottomDepth(iCell) < refBottomDepth(k)) then
                  maxLevelCell(iCell) = k-1
                  bottomDepth(iCell) = refBottomDepth(k-1)
                  exit
                end if
              end do
              do k = 1, nVertLevels
                layerThickness(k,iCell) = globalMaxBottomDepth * (interfaceLocations(k+1) - interfaceLocations(k))
              end do
            else
              ! impermeable layers
              do k = 1, nVertLevels
                layerThickness(k,iCell) = bottomDepth(iCell)/real(nVertLevels,RKIND)
              end do
            end if

            ! set layer thicknesses
            do k = 1, nVertLevels
              restingThickness(k, iCell) = layerThickness(k, iCell)
            end do
          end do
        end if


        do iCell = 1, nCellsSolve
          ! Set temperature
          idx = index_temperature
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = 20.0_RKIND
          end do

          ! Set salinity
          idx = index_salinity
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = 33.0_RKIND
          end do

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

    deallocate(interfaceLocations)
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_hurricane!}}}

!***********************************************************************
!
!  routine ocn_init_validate_hurricane
!
!> \brief   Validation for this initial condition
!> \author  Steven Brus
!> \date    01/08/19
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_hurricane(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_hurricane_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('hurricane')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_hurricane_vert_levels', config_hurricane_vert_levels)

      if(config_vert_levels <= 0 .and. config_hurricane_vert_levels > 0) then
         config_vert_levels = config_hurricane_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for hurricane case. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_hurricane!}}}


!***********************************************************************

end module ocn_init_hurricane

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
