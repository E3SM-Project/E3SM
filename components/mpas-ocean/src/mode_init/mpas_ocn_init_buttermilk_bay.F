! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_buttermilk_bay
!
!> \brief  initialize the Buttermilk Bay case
!> \author D. Wirasaet, S. Brus
!> \date   May-June 2022
!> \details
!>  This module contains the routines for initializing the Buttermilk
!>  Bay test case with or without subgrid corrections
!>
!-----------------------------------------------------------------------

module ocn_init_Buttermilk_bay

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

   use mpas_constants
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

   public :: ocn_init_setup_Buttermilk_bay, &
             ocn_init_validate_Buttermilk_bay, &
             ocn_init_Buttermilk_bay_bathymetry, &
             ocn_init_Buttermilk_bay_ssh, &
             ocn_init_Buttermilk_bay_velocity, &
             ocn_init_setup_Buttermilk_bay_read_topo

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   ! For netcdf topobathy input variables
   integer, public :: nLatTopo, nLonTopo

   type (field1DReal), public :: topoLat, topoLon
   type (field2DReal), public :: topoIC

   real(kind=RKIND), parameter:: eps = 1.0e-10_RKIND ;
!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_Buttermilk_bay
!
!> \brief   Setup for this initial condition
!> \author  D. Wirasaet and S. Brus
!> \date    May-June 2022
!> \details
!>  This routine sets up the initial conditions for this case.
!>   To be run in sigma vertical coordinates and single-layer
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_Buttermilk_bay(domain, iErr)!{{{
    use mpas_vector_operations  ! To calcutate edgeNormalVector

    implicit none
    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: forcingPool

    ! local variables
    integer :: iCell, iEdge, iVertex, k, idx
    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real (kind=RKIND) :: localVar1, localVar2
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity
    integer, pointer :: maxEdges

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: minLevelCell, maxLevelCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:,:), pointer :: verticesOnCell, verticesOnEdge
    integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, fCell, fEdge, fVertex, dcEdge
    real (kind=RKIND), dimension(:,:), pointer:: zMid

    real (kind=RKIND), dimension(:), pointer:: xEdge, yEdge, xVertex, yVertex
    real (kind=RKIND) :: minBottomDepth, maxBottomDepth, globalMaxBottomDepth, globalMinBottomDepth
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:), pointer :: areaCell
    real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
    ! Elevation Bcs
    real (kind=RKIND), dimension(:), pointer :: tidalInputMask


    real (kind=RKIND):: HH, uu, vv
    real (kind=RKIND):: RR, num, den
    real (kind=RKIND):: xshift = 0.0, yshift = 0.0
    real (kind=RKIND) :: layerThicknessEdgeAverage
    real (kind=RKIND), dimension(:,:), allocatable :: rSubgridPoints, sSubgridPoints
    real (kind=RKIND), dimension(:), allocatable :: subgridBathymetryValues, subgridAreas
    real (kind=RKIND), dimension(:), allocatable :: subgridSshValues
    real (kind=RKIND), dimension(:), allocatable :: subgridUValues, subgridVValues
    real (kind=RKIND), dimension(:), allocatable :: uVelocityAverage, vVelocityAverage
    integer :: nSubgridCell, nSubgridEdge, nSubgridVertex
    integer :: nSubgridTriPerSlice
    integer :: v1, v2
    integer :: c1, c2
    real (kind=RKIND) :: x(3), y(3)
    integer :: slice, nSlice
    real (kind=RKIND) :: deltaZ


    integer:: i, j, jj
    integer:: nsubgridCellEdge, iEdgeSegment
    real (kind=RKIND):: pi
    real (kind=RKIND), dimension(:,:), allocatable :: cellEdgeBathymetryValues
    real (kind=RKIND), dimension(:), allocatable:: dsEdge
    real (kind=RKIND), dimension(:), allocatable:: xSubgridCell, ySubgridCell
    real (kind=RKIND):: bathymetryMin, bathymetryMax
    iErr = 0

    if(config_init_configuration .ne. trim('buttermilk_bay')) return

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    nVertLevels  = config_Buttermilk_bay_vert_levels ;
    nVertLevelsP1 = nVertLevels + 1

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations, ocnConfigs ) ;
    !! Mental note: interfaceLocatons = (k-1)/N ;

    ! Initalize min/max values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    xMin = 1.0E10_RKIND
    xMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max values.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension( meshPool, 'nCellsSolve', nCellsSolve )
       call mpas_pool_get_dimension( meshPool, 'nEdgesSolve', nEdgesSolve )

       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

       yMin = min( yMin, minval(yCell(1:nCellsSolve)))
       yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
       xMin = min( xMin, minval(xCell(1:nCellsSolve)))
       xMax = max( xMax, maxval(xCell(1:nCellsSolve)))
       dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

       block_ptr => block_ptr % next
    end do

    ! Determine global min and max values.
    call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

    pi = acos(-1.0_RKIND)

    xshift = xMin
    yshift = (3.0_RKIND*yMin + dcEdgeMin*sin(pi/3.0_RKIND))/3.0_RKIND
    ! print*, "xMin, yMin = ", xMin, yMin, dcEdgeMin
    ! print*, "xshift, yshift = ", xshift, yshift

    !***********************************************************************
    !
    !  Topography
    !
    !***********************************************************************

    call mpas_log_write( 'Reading bathymetry from a NetCDF file')

    if (config_Buttermilk_bay_topography_source == 'latlon_file' .or. &
             config_Buttermilk_bay_topography_source == 'xy_file' ) then
       call mpas_log_write( 'Reading topography data from file.')
       call ocn_init_setup_Buttermilk_bay_read_topo(domain, iErr)
    endif

    !--------------------------------------------------------------------
    ! Use this section to set initial values
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
    call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
    call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
    call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

    call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
    call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)

    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels) ;
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve) ;
    call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve) ;
    call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve) ;
    call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)

    call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
    call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
    call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
    call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
    call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
    call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

    call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
    call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

    call mpas_pool_get_array(meshPool, 'fCell', fCell)
    call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
    call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

    call mpas_pool_get_array(meshPool, 'xEdge', xEdge )
    call mpas_pool_get_array(meshPool, 'yEdge', yEdge )
    call mpas_pool_get_array(meshPool, 'xVertex', xVertex )
    call mpas_pool_get_array(meshPool, 'yVertex', yVertex )

    call mpas_pool_get_array(statePool, 'zMid', zMid, 1) ;

    call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
    call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors ) ;
    call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity ) ;

    call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
    call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

    call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
    call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

    call mpas_pool_get_array(meshPool, 'subgridWetVolumeCellTable', &
                             subgridWetVolumeCellTable)
    call mpas_pool_get_array(meshPool, 'subgridWetVolumeEdgeTable', &
                             subgridWetVolumeEdgeTable)
    call mpas_pool_get_array(meshPool, 'subgridWetVolumeVertexTable', &
                             subgridWetVolumeVertexTable)
    call mpas_pool_get_array(meshPool, 'subgridWetFractionCellTable', &
                             subgridWetFractionCellTable)
    call mpas_pool_get_array(meshPool, 'subgridWetFractionEdgeTable', &
                             subgridWetFractionEdgeTable)
    call mpas_pool_get_array(meshPool,'subgridWetFractionVertexTable',&
                             subgridWetFractionVertexTable)
    call mpas_pool_get_array(meshPool, 'subgridSshCellTableRange', &
                             subgridSshCellTableRange)
    call mpas_pool_get_array(meshPool, 'subgridSshEdgeTableRange', &
                             subgridSshEdgeTableRange)
    call mpas_pool_get_array(meshPool, 'subgridSshVertexTableRange', &
                             subgridSshVertexTableRange)
    call mpas_pool_get_array(meshPool, 'subgridEdgeBathymetryMean', &
                             subgridEdgeBathymetryMean)
    call mpas_pool_get_array(meshPool, 'subgridVertexBathymetryMean', &
                             subgridVertexBathymetryMean)
    call mpas_pool_get_array(meshPool, 'subgridCellBathymetryMin', &
                             subgridCellBathymetryMin)
    call mpas_pool_get_array(meshPool, 'subgridEdgeBathymetryMin', &
                             subgridEdgeBathymetryMin)
    call mpas_pool_get_array(meshPool, 'subgridVertexBathymetryMin', &
                             subgridVertexBathymetryMin)
    call mpas_pool_get_dimension(meshPool, 'nSubgridTableLevels', &
                                 nSubgridTableLevels)
    call mpas_pool_get_array(meshPool, 'subgridLayerThicknessDebug', &
                             subgridLayerThicknessDebug)

    ! if config_buttermilk_bay_adjust_domain_center == .true.,
    ! Adjust center of the mesh so that its center is located at (0,0)
    if ( config_Buttermilk_bay_adjust_domain ) then
       xCell = xCell - xshift ;
       yCell = yCell - yshift ;

       xEdge = xEdge - xshift ;
       yEdge = yEdge - yshift ;

       xVertex = xVertex - xshift ;
       yVertex = yVertex - yshift ;

       ! get min, max coordinates of model domain   !
       ! after adjusting the coordinates             !
       yMin = min( yMin, minval(yCell(1:nCellsSolve)))
       yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
       xMin = min( xMin, minval(xCell(1:nCellsSolve)))
       xMax = max( xMax, maxval(xCell(1:nCellsSolve)))

       ! Determine global min and max values.
       call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
       call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
    end if

    ! Initlialze vector
    call mpas_initialize_vectors(meshPool) ;

    minLevelCell(:) = 1
    do iCell = 1, nCellsSolve
      ! Set up vertical grid
      maxLevelCell(iCell) = nVertLevels ; ! sigma coordinates
    end do


    do iCell = 1, nCellsSolve
       ! Set temperature
       activeTracers(index_temperature, :, iCell) = 10.0_RKIND

       ! Set salinity
       activeTracers(index_salinity, :, iCell) = 30.0_RKIND

       ! Set Coriolis parameters, if other than zero
       fCell(iCell) = 0.0_RKIND;
    end do

    do iEdge = 1, nEdgesSolve
       fEdge(iEdge) = 0.0_RKIND;
    end do

    do iVertex = 1, nVerticesSolve
       fVertex(iVertex) = 0.0_RKIND;
    end do

    allocate(uVelocityAverage(nEdgesSolve))
    allocate(vVelocityAverage(nEdgesSolve))

    if (config_use_subgrid_wetting_drying) then

       call ocn_subgrid_init(domain,iErr)
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
    end if

     ! Find max bottom depth
     maxBottomDepth = maxval( bottomDepth ) ;
     minBottomDepth = minval( bottomDepth ) ;
     call mpas_dmpar_max_real( domain % dminfo, maxBottomDepth, globalMaxBottomDepth ) ;
     call mpas_dmpar_min_real( domain % dminfo, minBottomDepth, globalMinBottomDepth ) ;

     ! Set refBottomDepth and refZMid
     do k = 1, nVertLevels
       refBottomDepth(k) = globalMaxBottomDepth*interfaceLocations(k+1) ;
       refZMid(k) = -0.5_RKIND*( interfaceLocations(k+1) + interfaceLocations(k))*globalMaxBottomDepth ;
     end do

     ! Set vertCoordMovementWeights
     vertCoordMovementWeights(:) = 1.0_RKIND

     ! Set velocity
     do iEdge = 1, nEdgesSolve

       if (config_use_subgrid_wetting_drying) then
          do k = 1, nVertLevels
             normalVelocity(k,iEdge) = uVelocityAverage(iEdge)*edgeNormalVectors(1,iEdge) &
                                     + vVelocityAverage(iEdge)*edgeNormalVectors(2,iEdge) ;
          end do
       else
          call ocn_init_Buttermilk_bay_velocity(xEdge(iEdge), yEdge(iEdge), uu, vv)
          do k = 1, nVertLevels
             normalVelocity(k,iEdge) = uu*edgeNormalVectors(1,iEdge) + vv*edgeNormalVectors(2,iEdge) ;
          end do
       end if
     end do

     ! Set layer thickness and ssh
     if (config_use_wetting_drying) then

       do iCell = 1, nCellsSolve
         ! Set up vertical grid
         maxLevelCell(iCell) = nVertLevels ; ! sigma coordinates
        end do

        do iCell = 1, nCellsSolve
          !
          ! make sure depth is thick enough via ssh = TOTAL_DEPTH - bottomDepth
          ! add a thin layer of nlayer*config_drying_min_cellhight
          !

          if (config_use_subgrid_wetting_drying) then
            ! Initial contion for a subgrid run
            call ocn_subgrid_ssh_lookup(layerThickness(1,iCell),&
                                        subgridWetVolumeCellTable(:,iCell),&
                                        subgridSshCellTableRange(:,iCell),&
                                        bottomDepth(iCell),&
                                        subgridCellBathymetryMin(iCell),&
                                        ssh(iCell))
            !call ocn_subgrid_layer_thickness_lookup(ssh(iCell), &
            !                               subgridWetVolumeCellTable(:,iCell), &
            !                               subgridSshCellTableRange(:,iCell),&
            !                               bottomDepth(iCell),&
            !                               LayerThickness(1,iCell))

          else
            ! Initial condition for a standard run
            call ocn_init_Buttermilk_bay_bathymetry(xCell(iCell), yCell(iCell), bottomDepth(iCell))

            ssh(iCell) = -bottomDepth(iCell) ;
            do k = 1, maxLevelCell(iCell)
              !
              layerThickness(k,iCell) = max( config_drying_min_cell_height, &
                  bottomDepth(iCell)/real(maxLevelCell(iCell),RKIND) )


              if (layerThickness(k,iCell) < config_drying_min_cell_height) then
                call mpas_log_write('layerThickness($i,$i)=$r', MPAS_LOG_CRIT, &
                  intArgs=(/k,iCell/), &
                  realArgs=(/layerThickness(k,iCell)/))
              end if

              ssh(iCell) = ssh(iCell) + layerThickness(k,iCell) ;
            end do
            !
          endif

          do k = 1, maxLevelCell(iCell)
            restingThickness(k,iCell) = bottomDepth(iCell)/maxLevelCell(iCell)
          end do

        end do

     end if


     if (config_use_subgrid_wetting_drying) then
        do iCell = 1,nCellsSolve
            call ocn_subgrid_layer_thickness_lookup(ssh(iCell), &
                                                 subgridWetVolumeCellTable(:,iCell), &
                                                 subgridSshCellTableRange(:,iCell),&
                                                 bottomDepth(iCell),&
                                                 subgridLayerThicknessDebug(iCell))
        enddo
     end if

     ! Set tidal boundary mask
     do iCell = 1, nCellsSolve
       tidalInputMask(iCell) = 0.0_RKIND
       if ( yCell(iCell) < (yMin+(dcEdgeMin*sin(pi/3.0_RKIND)/2.0_RKIND))  &
         .and.  yCell(iCell) > (yMin-(dcEDgeMin*sin(pi/3.0_RKIND)/2.0_RKIND)) ) then

          if ( (xCell(iCell) - dcEdgeMin/2.0_RKIND) > 2048.0_RKIND  &
           .and.  (xCell(iCell) + dcEdgeMin/2.0_RKIND) < 3072.0_RKIND ) then
              tidalInputMask(iCell) = 1.0_RKIND
          end if
        ! spread it over multiple cells
        ! if (yCell(iCell) > (25.0e3 - 3*dcEdgeMinGlobal)) then
        !  tidalInputMask(iCell) = exp(-((yCell(iCell)-25.0e3)/dcEdgeMinGlobal)**2.0)
       end if
     end do

    deallocate(interfaceLocations)
    if (config_global_ocean_topography_source == 'latlon_file') then
        call mpas_log_write( 'Cleaning up topography IC fields')
        call ocn_init_Buttermilk_bay_destroy_topo_fields()
    endif
    !--------------------------------------------------------------------

    print*, "****** End Butter milk bay init *****" ;

    return ;
  end subroutine ocn_init_setup_Buttermilk_bay!}}}

!***********************************************************************
!
!  routine ocn_init_Buttermilk_bay_bathymetry
!
!> \brief   Interpolate bathymetry
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Return the value of the bathymetry at a given x,y point
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_init_Buttermilk_bay_bathymetry(xin, yin, depth)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: xin, yin
      real (kind=RKIND), intent(out) :: depth

      integer:: ix(2)
      real (kind=RKIND) :: x, y
      real (kind=RKIND) :: xc(2)
      real (kind=RKIND) :: dx(2), x0(2), xN(2), bath(4), f(4), den
      real (kind=RKIND) :: xv(2), yv(2)

      !--------------------------------------------------------------------
      !
      ! DEM must be on a uniform grid

      y = yin 
      x = xin

      if (config_Buttermilk_bay_topography_latlon_degrees .and. &
          config_Buttermilk_bay_topography_source == 'latlon_file' ) then

          if (x > pii) then
              x = x - 2.0_RKIND*pii
          end if

      end if
     

      x0 = (/ topoLon%array(1), topoLat%array(1) /) ;
      xN = (/ topoLon%array(nLonTopo), topoLat%array(nLatTopo) /) ;
      dx = (/ topoLon%array(2) - topoLon%array(1), &
              topoLat%array(2) - topoLat%array(1) /) ;

      ! Bilienar interpolation
      xc = (/ x, y /) ;
      ix = floor( (xc - x0)/dx ) + 1 ;

      if ( ( ix(1) >= 1 .and. ix(1) < nLonTopo ) .and. &
              ( ix(2) >= 1 .and. ix(2) < nLatTopo ) ) then
        ! include the west and soutth bourdaries of
        ! a given uniform raster DEM
        bath(1) = TopoIC%array(ix(1),ix(2))  ;
        bath(2) = TopoIC%array(ix(1)+1,ix(2)) ;
        bath(3) = TopoIC%array(ix(1)+1,ix(2)+1) ;
        bath(4) = TopoIC%array(ix(1),ix(2)) ;

        xv(1) = topoLon%array(ix(1)) ;
        xv(2) = topoLon%array(ix(1)+1) ;

        yv(1) = topoLat%array(ix(2)) ;
        yv(2) = topoLat%array(ix(2)+1) ;

        den = dx(1)*dx(2) ;

        f(1) = ( xc(1) - xv(2) )*( xc(2) - yv(2) )/den ;
        f(2) = ( xc(1) - xv(1) )*( xc(2) - yv(2) )/(-den) ;
        f(3) = ( xc(1) - xv(1) )*( xc(2) - yv(1) )/den  ;
        f(4) = ( xc(1) - xv(2) )*( xc(2) - yv(1) )/(-den) ;

        depth = sum( bath*f ) ;
        !depth = -1.0_RKIND*depth

        if (minval(bath) < -99998) then
          depth = -1d10
        end if
        if (maxval(bath) > 99998) then
          depth = -1d10
        end if
      else
        !! nearest extrapolation !
        !ix(1) = merge( 1, ix(1), ix(1) < 1 ) ;
        !ix(1) = merge( nLonTopo, ix(1), ix(1) >= nLonTopo ) ;
        !ix(2) = merge( 1, ix(2), ix(2) < 1 )  ;
        !ix(2) = merge( nLatTopo, ix(2), ix(2) >= nLatTopo ) ;

        !depth = TopoIC%array(ix(1),ix(2)) ;
        depth = -1e10
      endif

      return ;
   end subroutine ocn_init_Buttermilk_bay_bathymetry!}}}

!***********************************************************************
!
!  routine ocn_init_Buttermilk_bay_ssh
!
!> \brief   Compute initial ssh field
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Use exact solution to compute ssh field for initial conditions
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_init_Buttermilk_bay_ssh(x, y, bottomDepth, ssh)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: x, y
      real (kind=RKIND), intent(in) :: bottomDepth
      real (kind=RKIND), intent(out) :: ssh
      real (kind=RKIND) :: RR

   !--------------------------------------------------------------------

      ssh = 0.0_RKIND ;
      ssh = - bottomDepth + max( ssh + bottomDepth, 0.0_RKIND ) ;

      return
   !--------------------------------------------------------------------

   end subroutine ocn_init_Buttermilk_bay_ssh!}}}

!***********************************************************************
!
!  routine ocn_init_Buttermilk_bay_velocity
!
!> \brief   Compute initial velocity field
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Use exact solution to comupte velocity field for initial conditions
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_init_Buttermilk_bay_velocity(x, y, u, v)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: x, y
      real (kind=RKIND), intent(out) :: u, v
      real (kind=RKIND) :: RR, HH

   !--------------------------------------------------------------------

      u = 0.0_RKIND ;
      v = 0.0_RKIND ;
   !--------------------------------------------------------------------

   end subroutine ocn_init_Buttermilk_bay_velocity!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_topo
!
!> \brief   Read the topography IC file
!> \author  D. Wirasaet, S. Brus
!> \date    Sep 2023
!> \details
!>  This routine reads the topography IC file, including latitude and longitude
!>   information for topography data.
!>
!>  Adapted from ocn_init_setup_global_ocean_read_topo
!-----------------------------------------------------------------------
    subroutine ocn_init_setup_Buttermilk_bay_read_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr
       integer :: i

       type (MPAS_Stream_type) :: topographyStream

       iErr = 0

       ! Define stream for depth levels
       call MPAS_createStream(topographyStream, domain % iocontext, &
                config_Buttermilk_bay_topography_file, MPAS_IO_NETCDF, &
                                                   MPAS_IO_READ, ierr=iErr)

       ! Setup topoLat, topoLon, and topoIC fields for stream to be read in
       topoLat % fieldName = trim(config_Buttermilk_bay_topography_lat_varname)
       topoLat % dimSizes(1) = nLatTopo
       topoLat % dimNames(1) = trim(config_Buttermilk_bay_topography_nlat_dimname)
       topoLat % isVarArray = .false.
       topoLat % isPersistent = .true.
       topoLat % isActive = .true.
       topoLat % hasTimeDimension = .false.
       topoLat % block => domain % blocklist
       allocate(topoLat % attLists(1))
       allocate(topoLat % array(nLatTopo))

       topoLon % fieldName = trim(config_Buttermilk_bay_topography_lon_varname)
       topoLon % dimSizes(1) = nLonTopo
       topoLon % dimNames(1) = trim(config_Buttermilk_bay_topography_nlon_dimname)
       topoLon % isVarArray = .false.
       topoLon % isPersistent = .true.
       topoLon % isActive = .true.
       topoLon % hasTimeDimension = .false.
       topoLon % block => domain % blocklist
       allocate(topoLon % attLists(1))
       allocate(topoLon % array(nLonTopo))

       topoIC % fieldName = trim(config_Buttermilk_bay_topography_varname)
       topoIC % dimSizes(1) = nLonTopo
       topoIC % dimSizes(2) = nLatTopo
       topoIC % dimNames(1) = trim(config_Buttermilk_bay_topography_nlon_dimname)
       topoIC % dimNames(2) = trim(config_Buttermilk_bay_topography_nlat_dimname)
       topoIC % isVarArray = .false.
       topoIC % isPersistent = .true.
       topoIC % isActive = .true.
       topoIC % hasTimeDimension = .false.
       topoIC % block => domain % blocklist
       allocate(topoIC % attLists(1))
       allocate(topoIC % array(nLonTopo, nLatTopo))


       ! Add topoLat, topoLon, and topoIC fields to stream
       call MPAS_streamAddField(topographyStream, topoLat, iErr)
       call MPAS_streamAddField(topographyStream, topoLon, iErr)
       call MPAS_streamAddField(topographyStream, topoIC, iErr)

       ! Read stream
       call MPAS_readStream(topographyStream, 1, iErr)
       topoIC%array = -topoIC%array ;


       ! Close stream
       call MPAS_closeStream(topographyStream)

       if ( config_Buttermilk_bay_topography_latlon_degrees .and. &
              config_Buttermilk_bay_topography_source == 'latlon_file' ) then
          !do i = 1,nLonTopo
          !   if (topoLon % array(i) < 0.0) then
          !      topoLon % array(i) = topoLon % array(i) + 360.0_RKIND
          !   end if
          !end do

          topoLat % array(:) = topoLat % array(:) * pii / 180.0_RKIND
          topoLon % array(:) = topoLon % array(:) * pii / 180.0_RKIND
       end if

    end subroutine ocn_init_setup_Buttermilk_bay_read_topo!}}}


!***********************************************************************
!
!  routine ocn_init_Buttermilk_bay_destroy_topo_fields
!
!> \brief   Topography field cleanup routine
!> \author  D. Wirasaet and S. Brus
!> \date    Sep 2023
!> \details
!>  This routine destroys the fields that were created to hold topography
!>  initial condition information
!>
!>  NOTE: adapteed from ocn_init_global_ocaen_destroy_topo_fileds
!-----------------------------------------------------------------------

    subroutine ocn_init_Buttermilk_bay_destroy_topo_fields()!{{{
       implicit none

       deallocate(topoIC % array)
       deallocate(topoLat % array)
       deallocate(topoLon % array)
    end subroutine ocn_init_Buttermilk_bay_destroy_topo_fields!}}}


!***********************************************************************
!
!  routine ocn_init_Buttermilk_bay
!
!> \brief   Validation for this initial condition
!> \author  D. Wirasaet and Steven Brus
!> \date   Sep 2022
!> \details
!>  This routine validates the configuration options for this case.
!>
!-----------------------------------------------------------------------
   subroutine ocn_init_validate_Buttermilk_bay(configPool, packagePool, iocontext, iErr)!{{{
      implicit none

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout), target :: iocontext

      integer, intent(out) :: iErr

      ! character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_Buttermilk_bay_vert_levels
      integer, pointer :: config_subgrid_table_levels, config_Buttermilk_bay_subgrid_table_levels


      type (mpas_io_context_type), pointer :: iocontext_ptr
      type (MPAS_IO_Handle_type) :: inputFile
      character (len=StrKIND), pointer :: config_init_configuration, &
                                          config_Buttermilk_bay_topography_source, &
                                          config_Buttermilk_bay_topography_file, &
                                          config_Buttermilk_bay_topography_nlat_dimname, &
                                          config_Buttermilk_bay_topography_nlon_dimname

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('buttermilk_bay')) return

      iocontext_ptr => iocontext

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool,  &
                                'config_Buttermilk_bay_vert_levels', &
                                config_Buttermilk_bay_vert_levels)

      if(config_vert_levels <= 0 .and. config_Buttermilk_bay_vert_levels > 0) then
         config_vert_levels = config_Buttermilk_bay_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for Buttermilk bay.'// &
           'Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
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

      call mpas_log_write( ' Done ocn_init_validate_buttermilk_bay ' )
   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_ButterMilk_bay!}}}

!***********************************************************************

end module ocn_init_Buttermilk_bay

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
