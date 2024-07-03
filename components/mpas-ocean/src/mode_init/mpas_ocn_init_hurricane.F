! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use ocn_config
   use ocn_init_vertical_grids
   use ocn_init_cell_markers
   use ocn_subgrid
   use ocn_init_subgrid
   use ocn_init_Buttermilk_bay

   use mpas_io
   use mpas_io_streams
   use mpas_stream_manager


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

   !! For netcdf topobathy input variables
   !integer :: nLatTopo, nLonTopo

   !type (field1DReal) :: topoLat, topoLon
   !type (field2DReal) :: topoIC

   !real(kind=RKIND), parameter:: eps = 1.0e-10_RKIND ;

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
    type (mpas_pool_type), pointer :: forcingPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool

    ! local variables
    integer :: iCell, k, idx, iEdge, iVertex
    real (kind=RKIND) :: dlon, dlat, dw
    real (kind=RKIND) :: maxBottomDepth, globalMaxBottomDepth, globalMinBottomDepth
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations
    real (kind=RKIND), parameter :: eps=1.0e-12
    real (kind=RKIND), parameter :: pi = 4.0_RKIND*ATAN(1.0_RKIND)

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:), pointer :: bottomDrag
    real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, bottomDepthObserved, &
         fCell, fEdge, fVertex, &
         latCell, latEdge, latVertex, &
         lonCell
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
    real (kind=RKIND), dimension(:), allocatable :: uVelocityAverage, vVelocityAverage

    iErr = 0

    if(config_init_configuration .ne. trim('hurricane')) return

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
    call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
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
    !write(*,*) globalMaxBottomDepth, globzalMinBottomDepth

    if (config_use_subgrid_wetting_drying) then

       call ocn_subgrid_init(domain,iErr)
       if (MPAS_stream_mgr_stream_exists(domain%streamManager, 'subgrid_init')) then

          allocate(uVelocityAverage(nEdgesSolve))
          allocate(vVelocityAverage(nEdgesSolve))

          if (config_Buttermilk_bay_topography_source == 'latlon_file' .or. &
                   config_Buttermilk_bay_topography_source == 'xy_file' ) then
             call mpas_log_write( 'Reading subgrid topography data from file.')
             call ocn_init_setup_Buttermilk_bay_read_topo(domain, iErr)
          endif


          call ocn_init_subgrid_calculations(domain, &
                                             ocn_init_Buttermilk_bay_bathymetry, &
                                             ocn_init_Buttermilk_bay_velocity, &
                                             ocn_init_Buttermilk_bay_ssh, &
                                             config_Buttermilk_bay_subgrid_refinement_level, &
                                             config_Buttermilk_bay_subgrid_edge_bathymetry_max_pixel, &
                                             config_Buttermilk_bay_subgrid_use_thin_layer, &
                                             uVelocityAverage, &
                                             vVelocityAverage, &
                                             iErr)
       else


         call mpas_log_write( 'Setting default subgrid table values')
         do iCell = 1, nCellsSolve
            subgridWetVolumeCellTable(nSubgridTableLevels,iCell) = -1e10
            subgridSshCellTableRange(2,iCell) = -1e10
         end do
         do iEdge = 1, nEdgesSolve
            subgridWetVolumeEdgeTable(nSubgridTableLevels,iEdge) = -1e10
            subgridSshEdgeTableRange(2,iEdge) = -1e10
         end do
         do iVertex = 1, nVerticesSolve
            subgridWetVolumeVertexTable(nSubgridTableLevels,iVertex) = -1e10
            subgridSshVertexTableRange(2,iVertex) = -1e10
         end do

       end if
    end if

    !--------------------------------------------------------------------
    ! Use this section to set initial values
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
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
         !do iCell = 1, nCellsSolve
         !  ssh(iCell) = 0.0_RKIND
         !end do
       end if


       do iCell = 1, nCellsSolve

          bottomDepth(iCell) = -1.0_RKIND*bottomDepthObserved(iCell)

          ! sea level rise adjustment
          bottomDepth(iCell) = bottomDepth(iCell) + config_hurricane_sea_level_rise_adjustment

          ! enforce minimum and maximum bottom depth
          bottomDepth(iCell) = max(bottomDepth(iCell), config_hurricane_min_depth)
          bottomDepth(iCell) = min(bottomDepth(iCell), globalMaxBottomDepth)

        end do

        if (config_use_wetting_drying) then
          do iCell = 1, nCellsSolve
            ! Set up vertical grid
            maxLevelCell(iCell) = nVertLevels
            !! ensures only 3 layers can be obtained on mesh
            !if (bottomDepth(iCell) < 0 .or. bottomDepth(iCell) < refBottomDepth(3)) then
            !  maxLevelCell(iCell) = 3
            !else
            !  ! snap to zstar coordinate
            !  do k = 1, nVertLevels
            !    if (bottomDepth(iCell) < refBottomDepth(k)) then
            !      maxLevelCell(iCell) = k-1
            !      bottomDepth(iCell) = refBottomDepth(k-1)
            !      exit
            !    end if
            !  end do
            !end if
          end do

          do iCell = 1, nCellsSolve

            if (config_use_subgrid_wetting_drying .and. MPAS_stream_mgr_stream_exists(domain%streamManager, 'subgrid_init')) then
              ! Initial contion for a subgrid run
              !call ocn_subgrid_ssh_lookup(layerThickness(1,iCell),&
              !                            subgridWetVolumeCellTable(:,iCell),&
              !                            subgridSshCellTableRange(:,iCell),&
              !                            bottomDepth(iCell),&
              !                            subgridCellBathymetryMin(iCell),&
              !                            ssh(iCell))
              ! if ( ssh(iCell) < subgridSshCellTableRange(3,iCell) ) then
              !    ssh(iCell) = subgridSshCellTableRange(3,iCell) ;   
              ! end if 

            else
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
            endif

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

       ! Set depth-variable drag
       if (config_use_variable_drag) then
         call mpas_pool_get_array(forcingPool, 'bottomDrag', bottomDrag)
         do iCell = 1, nCellsSolve
           if (bottomDepth(iCell) <= config_hurricane_land_z_limit) then
             bottomDrag(iCell) = config_hurricane_land_drag
           else if (config_hurricane_land_z_limit < bottomDepth(iCell) .and. bottomDepth(iCell) < config_hurricane_marsh_z_limit) then
             bottomDrag(iCell) = config_hurricane_marsh_drag
           else if (config_hurricane_marsh_z_limit <= bottomDepth(iCell)) then
             bottomDrag(iCell) = config_hurricane_channel_drag
           else
             call mpas_log_write('Default value for drag is not selected' &
               // ' properly for bottomDepth of $r in cell $i!', &
             realArgs=(/bottomDepth(iCell)/), intArgs=(/iCell/))
           end if
         end do
       end if


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
      type (mpas_io_context_type), intent(inout), target :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_hurricane_vert_levels
      integer, pointer :: config_subgrid_table_levels, config_Buttermilk_bay_subgrid_table_levels

      type (mpas_io_context_type), pointer :: iocontext_ptr
      type (MPAS_IO_Handle_type) :: inputFile
      character (len=StrKIND), pointer :: config_Buttermilk_bay_topography_source, &
                                          config_Buttermilk_bay_topography_file, &
                                          config_Buttermilk_bay_topography_nlat_dimname, &
                                          config_Buttermilk_bay_topography_nlon_dimname

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('hurricane')) return

      iocontext_ptr => iocontext

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_hurricane_vert_levels', config_hurricane_vert_levels)

      if(config_vert_levels <= 0 .and. config_hurricane_vert_levels > 0) then
         config_vert_levels = config_hurricane_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for hurricane case. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

      call mpas_pool_get_config(configPool, 'config_subgrid_table_levels', config_subgrid_table_levels)
      call mpas_pool_get_config(configPool, &
                                'config_Buttermilk_bay_subgrid_table_levels', &
                                config_Buttermilk_bay_subgrid_table_levels)
      if (config_subgrid_table_levels <=0 .and. config_Buttermilk_bay_subgrid_table_levels >0) then
         config_subgrid_table_levels = config_Buttermilk_bay_subgrid_table_levels
      else if (config_subgrid_table_levels <=0) then
         call mpas_log_write( 'Validation failed for Buttermilk bay.'// &
           'Not given a usable value for subgrid table levels.', MPAS_LOG_CRIT)
         iErr = 1 
      end if

      !---------- adapted from ocn_init_validate_global_ocean
      call mpas_pool_get_config(configPool, 'config_Buttermilk_bay_topography_source', &
                                config_Buttermilk_bay_topography_source)
      call mpas_pool_get_config(configPool, 'config_Buttermilk_bay_topography_file', &
                                config_Buttermilk_bay_topography_file)
      call mpas_pool_get_config(configPool, 'config_Buttermilk_bay_topography_nlat_dimname', &
                                config_Buttermilk_bay_topography_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_Buttermilk_bay_topography_nlon_dimname', &
                                config_Buttermilk_bay_topography_nlon_dimname)


      call mpas_log_write( config_Buttermilk_bay_topography_source  )
      call mpas_log_write( config_Buttermilk_bay_topography_file )
      call mpas_log_write( config_Buttermilk_bay_topography_nlat_dimname )
      call mpas_log_write( config_Buttermilk_bay_topography_nlon_dimname )

      if (config_Buttermilk_bay_topography_source /= 'latlon_file' .and. &
          config_Buttermilk_bay_topography_source /= 'xy_file') then
         call mpas_log_write( 'Unexpected value for &
                   config_Buttermilk_bay_topography_source: ' &
             // trim(config_Buttermilk_bay_topography_source), MPAS_LOG_CRIT)
         iErr = 1 
         return
      end if

      if (config_Buttermilk_bay_topography_file == 'none' .and. &
          (config_Buttermilk_bay_topography_source == 'latlon_file' .or. &
          config_Buttermilk_bay_topography_source == 'xy_file') ) then
         call mpas_log_write( 'Validation failed for Buttermilk bay test case. ' &
           // 'Invalid filename for config_Buttermilk_bay_topography_file ' &
           // config_Buttermilk_bay_topography_file // '  ' &
           // config_Buttermilk_bay_topography_source , MPAS_LOG_CRIT)
         iErr = 1 
         return
      end if

      call mpas_log_write( ' in ocn_init_validate_buttermilk_bay '//config_Buttermilk_bay_topography_source )

      if (config_Buttermilk_bay_topography_source == 'latlon_file' .or. &
           config_Buttermilk_bay_topography_source == 'xy_file'  ) then

         inputFile = MPAS_io_open(trim(config_Buttermilk_bay_topography_file), &
                        MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
         if (iErr /= 0) then
            call mpas_log_write( 'could not open file '// &
                   trim(config_Buttermilk_bay_topography_file), MPAS_LOG_CRIT)
            return
         end if

         call MPAS_io_inq_dim(inputFile, config_Buttermilk_bay_topography_nlat_dimname, nLatTopo, iErr)
         call MPAS_io_inq_dim(inputFile, config_Buttermilk_bay_topography_nlon_dimname, nLonTopo, iErr)

         call MPAS_io_close(inputFile, iErr)

      end if
      !----------


   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_hurricane!}}}


!***********************************************************************

end module ocn_init_hurricane

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
