! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_subgrid
!
!> \brief  Calculate the subgrid information 
!> \author D. Wirasaet, S. Brus
!> \date   May-June 2022
!> \details
!>  This module contains the routines for calculating the information
!>  needed to use subgrid corrctions within a forward model run.
!>
!-----------------------------------------------------------------------

module ocn_init_subgrid

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_config
   use ocn_subgrid

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

   public :: ocn_init_subgrid_calculations, &
             ocn_init_subgrid_interpolate_topo, &
             ocn_init_subgrid_interpolate_drag, &
             ocn_init_read_subgrid_data

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real(kind=RKIND), parameter:: eps = 1.0e-10_RKIND ;   
   real(kind=RKIND) :: subgrid_thin_layer

   abstract interface
      subroutine bathymetry_function(x, y, b)
         use mpas_kind_types, only: rkind
         implicit none
         real (kind=RKIND), intent(in) :: x
         real (kind=RKIND), intent(in) :: y
         real (kind=RKIND), intent(out) :: b
      end subroutine
      subroutine velocity_function(x, y, u, v)
         use mpas_kind_types, only: rkind
         implicit none
         real (kind=RKIND), intent(in) :: x
         real (kind=RKIND), intent(in) :: y
         real (kind=RKIND), intent(out) :: u
         real (kind=RKIND), intent(out) :: v
      end subroutine
      subroutine ssh_function(x, y, b, z)
         use mpas_kind_types, only: rkind
         implicit none
         real (kind=RKIND), intent(in) :: x
         real (kind=RKIND), intent(in) :: y
         real (kind=RKIND), intent(in) :: b
         real (kind=RKIND), intent(out) :: z
      end subroutine
   end interface

   integer, parameter :: interp_bilinear = 1
   integer, parameter :: interp_nn = 2

   type (field1DReal), public :: topoLat, topoLon
   type (field2DReal), public :: topoData

   type (field1DReal), public :: lulcLat, lulcLon
   type (field2DReal), public :: lulcData

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_subgrid_calculations
!
!> \brief   Calculations for subgrid look-up tables 
!> \author  D. Wirasaet and S. Brus
!> \date    May-June 2022
!> \details
!>    This routine performs the vertical integration of the wet fraction
!>    over the cell, edge, and vertex control volumes to create
!>    lookup tables for the relationship between ssh and wet volume per
!>    unit area.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_subgrid_calculations(domain, problem_bathymetry, problem_drag, problem_velocity, problem_ssh, &
                                           subgrid_refinement_level, &
                                           subgrid_edge_bathymetry_max_pixel, &
                                           subgrid_use_thin_layer, &
                                           uVelocityAverage, vVelocityAverage, iErr)!{{{
    
    implicit none
    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr
    procedure(bathymetry_function) :: problem_bathymetry
    procedure(bathymetry_function) :: problem_drag
    procedure(velocity_function) :: problem_velocity
    procedure(ssh_function) :: problem_ssh
    integer, intent(in) :: subgrid_refinement_level
    logical, intent(in) :: subgrid_edge_bathymetry_max_pixel
    logical, intent(in) :: subgrid_use_thin_layer
    
    real (kind=RKIND), dimension(:), intent(out) :: uVelocityAverage, vVelocityAverage

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool

    ! local variables
    integer :: iCell, iEdge, iVertex, k, idx
    integer :: count
    real(kind=RKIND) :: lon, lat
    real(kind=RKIND) :: lon_min, lon_max, lat_min, lat_max
    real(kind=RKIND) :: bathyVertex

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve 
    integer, pointer :: maxEdges

    ! Define variable pointers
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:,:), pointer :: verticesOnCell, verticesOnEdge
    integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, bottomDepth 
    real (kind=RKIND), dimension(:,:), pointer:: edgeNormalVectors
    real (kind=RKIND), dimension(:), pointer:: xEdge, yEdge, xVertex, yVertex
    
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness
    real (kind=RKIND), dimension(:), pointer :: areaCell

    real (kind=RKIND) :: layerThicknessEdgeAverage
    real (kind=RKIND), dimension(:,:), allocatable :: rSubgridPoints, sSubgridPoints
    real (kind=RKIND), dimension(:), allocatable :: subgridBathymetryValues, subgridAreas
    real (kind=RKIND), dimension(:), allocatable :: subgridSshValues
    real (kind=RKIND), dimension(:), allocatable :: subgridUValues, subgridVValues
    real (kind=RKIND), dimension(:), allocatable :: subgridDragValues

    integer :: nSubgridCell, nSubgridEdge, nSubgridVertex
    integer :: nSubgridTriPerSlice
    integer :: v1, v2
    integer :: c1, c2
    real (kind=RKIND) :: x(3), y(3)
    integer :: slice, nSlice
    real (kind=RKIND) :: deltaZ
    logical, pointer :: on_a_sphere
    integer :: startIDx, endIdx



    integer:: i, j, jj
    integer:: nsubgridCellEdge, iEdgeSegment
    real (kind=RKIND):: pi
    real (kind=RKIND), dimension(:,:), allocatable :: cellEdgeBathymetryValues
    real (kind=RKIND), dimension(:), allocatable:: dsEdge
    real (kind=RKIND), dimension(:), allocatable:: xSubgridCell, ySubgridCell
    real (kind=RKIND):: bathymetryMin, bathymetryMax 
    iErr = 0


    block_ptr => domain % blocklist
    call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
    call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve) ;
    call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve) ; 
    call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve) ; 
    call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
    
    call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

    call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
    call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

    if (on_a_sphere) then
      call mpas_pool_get_array(meshPool, 'lonCell', xCell)
      call mpas_pool_get_array(meshPool, 'latCell', yCell)
      call mpas_pool_get_array(meshPool, 'lonEdge', xEdge ) 
      call mpas_pool_get_array(meshPool, 'latEdge', yEdge ) 
      call mpas_pool_get_array(meshPool, 'lonVertex', xVertex ) 
      call mpas_pool_get_array(meshPool, 'latVertex', yVertex ) 
    else
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge ) 
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge ) 
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex ) 
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex ) 
    endif

    call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
    call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors ) ; 
    call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity ) ; 
    call mpas_pool_get_array(statePool, 'layerThickness', layerThickness) ; 

    if (subgrid_use_thin_layer) then
       subgrid_thin_layer = config_drying_min_cell_height + eps
    else
       subgrid_thin_layer = 0.0_RKIND
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Cells
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, "Begin cells", nCellsSolve

    lon_min = -79.0
    lon_max = -70.0
    lat_min = 35.0
    lat_max = 42.0

    nSubgridTriPerSlice = subgrid_refinement_level**2
    allocate(rSubgridPoints(3,maxEdges*nSubgridTriPerSlice), sSubgridPoints(3,maxEdges*nSubgridTriPerSlice))
    call ocn_init_define_subgrid_points(subgrid_refinement_level, rSubgridPoints, sSubgridPoints)
    allocate(subgridBathymetryValues(maxEdges*nSubgridTriPerSlice), subgridAreas(maxEdges*nSubgridTriPerSlice))
    allocate(subgridDragValues(maxEdges*nSubgridTriPerSlice))
    allocate(subgridSshValues(maxEdges*nSubgridTriPerSlice))
    allocate(subgridUValues(maxEdges*nSubgridTriPerSlice))
    allocate(subgridVValues(maxEdges*nSubgridTriPerSlice))

    allocate(xSubgridCell(maxEdges*nSubgridTriPerSlice)) 
    allocate(ySubgridCell(maxEdges*nSubgridTriPerSlice)) 
    xSubgridCell = 0.0_RKIND  
    ySubgridCell = 0.0_RKIND 

cell:do iCell = 1,nCellsSolve

       do i = 1,nEdgesOnCell(iCell)
         v1 = verticesOnCell(i,iCell)
         
         lon = xVertex(v1)
         lat = yVertex(v1)

         call problem_bathymetry(lon, lat, bathyVertex)

         if (bathyVertex < -99998) then
            cycle cell 
         endif
       end do

      !
      ! Evaluate subgrid bathymetry at centers of sub-triangles for cell slices 
      ! (all subdivided triangles for each cell slice are gathered into 
      !  subgridBathymetryValues and subgridAreas)
      !---------------------------------------------------------------
      nSubgridCell = 0  ! Counter for all subgrid triangles over cell slices
      do slice = 1,nEdgesOnCell(iCell) ! Loop over cell slices

        v1 = verticesOnCell(slice,iCell)
        if (slice+1 <= nEdgesOnCell(iCell)) then
          v2 = verticesOnCell(slice+1,iCell)
        else
          v2 = verticesOnCell(1,iCell)
        endif

        ! Cell slice coordinates
        x(1) = xCell(iCell)
        y(1) = yCell(iCell)

        x(2) = xVertex(v1)
        y(2) = yVertex(v1)

        x(3) = xVertex(v2)
        y(3) = yVertex(v2)

        call ocn_init_evaluate_subgrid_data(x, y, nSubgridTriPerSlice, nSubgridCell, rSubgridPoints, sSubgridPoints, &
                                            problem_bathymetry, problem_velocity, problem_ssh, &  
                                            subgridBathymetryValues, subgridAreas, subgridSshValues, &
                                            subgridUValues, subgridVValues, xSubgridCell, ySubgridCell)
      enddo
    
      ! Evaluate bounds of look-up table range
      !---------------------------------------------------------------        
      bathymetryMin = maxval( subgridBathymetryValues(1:nSubgridCell) ) 
      bathymetryMax = minval( subgridBathymetryValues(1:nSubgridCell) )

      if (bathymetryMax < -99998) then
         cycle cell 
      endif
      if (bathymetryMin > 99998) then
         cycle cell
      endif

      call ocn_init_set_table_range(1, nSubgridCell, subgridBathymetryValues, subgridSshCellTableRange(:,iCell))

      ! Evaluate subgrid bathymetry
      !---------------------------------------------------------------
      bottomDepth(iCell) = sum(subgridBathymetryValues(1:nSubgridCell)*subgridAreas(1:nSubgridCell))/sum(subgridAreas(1:nSubgridCell))
      subgridCellBathymetryMin(iCell) = maxval(subgridBathymetryValues(1:nSubgridCell))

      ! Vertical integration of wet fraction
      !---------------------------------------------------------------

      call ocn_init_vertical_integration(iCell,subgridSshCellTableRange, nSubgridCell, subgridBathymetryValues, subgridAreas, &
                                         subgridWetVolumeCellTable, subgridWetFractionCellTable)

      ! Evaluate wet layerThickness average
      !---------------------------------------------------------------
        call ocn_init_grid_average(nSubgridCell, subgridBathymetryValues, subgridSshValues, subgridAreas, layerThickness(1,iCell))
!        call ocn_init_wet_average_ssh(nSubgridCell, subgridBathymetryValues, subgridSshValues, subgridAreas, ssh(iCell))!{{{

        call ocn_subgrid_ssh_lookup(layerThickness(1,iCell),&
                                    subgridWetVolumeCellTable(:,iCell),&
                                    subgridSshCellTableRange(:,iCell),&
                                    bottomDepth(iCell),&
                                    subgridCellBathymetryMin(iCell),&
                                    ssh(iCell))

!            call ocn_subgrid_layer_thickness_lookup(ssh(iCell), &
!                                           subgridWetVolumeCellTable(:,iCell), &
!                                           subgridSshCellTableRange(:,iCell),&
!                                           bottomDepth(iCell),&
!                                           LayerThickness(1,iCell))

       call ocn_subgrid_ssh_lookup( config_drying_min_cell_height, & 
                                     subgridWetVolumeCellTable(:,iCell), &
                                     subgridSshCellTableRange(:,iCell),&
                                     bottomDepth(iCell),&
                                     subgridCellBathymetryMin(iCell),&
                                     subgridSShCellTableRange(3,iCell) )

                             
       if ( ssh(iCell) < subgridSshCellTableRange(3,iCell) ) then
          ssh(iCell) = subgridSshCellTableRange(3,iCell) ;  

!          LayerThickness(1,iCell) = config_drying_min_cell_height ; 
       end if 

    end do cell

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Edges
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, "Begin edges", nEdgesSolve

    nSubgridCellEdge=subgrid_refinement_level 
    allocate( cellEdgeBathymetryValues(3,nSubgridCellEdge) )   
    allocate( dsEdge(nSubgridCellEdge) ) 
    cellEdgeBathymetryValues = -99999 ; 


edge:do iEdge = 1,nEdgesSolve

       do i = 1,2
         v1 = verticesOnEdge(i,iEdge)
         
         lon = xVertex(v1)
         lat = yVertex(v1)

         call problem_bathymetry(lon, lat, bathyVertex)

         if (bathyVertex < -99998) then
            cycle edge 
         endif
       end do
      !--------------------------------------------------------------
      ! Evaluate subgrid bathymetry at centers of sub-triangles for edge slices
      ! (all subdivided triangles for each edge slice are gathered into 
      !  subgridBathymetryValues and subgridAreas)
      !---------------------------------------------------------------

      nSlice = 0
      do slice = 1,2
        if (cellsOnEdge(slice,iEdge) <= nCellsSolve) then
          nSlice = nSlice + 1
        endif
      enddo

      nSubgridEdge = 0  ! Counter for all subgrid triangles over edge slices
      do slice = 1,nSlice ! Loop over edge slices

        ! Edge slice coordinates
        x(1) = xVertex(verticesOnEdge(1,iEdge))
        y(1) = yVertex(verticesOnEdge(1,iEdge))

        x(2) = xVertex(verticesOnEdge(2,iEdge))
        y(2) = yVertex(verticesOnEdge(2,iEdge))

        x(3) = xCell(cellsOnEdge(slice,iEdge))
        y(3) = yCell(cellsOnEdge(slice,iEdge))


        call ocn_init_evaluate_subgrid_data(x, y, nSubgridTriPerSlice, nSubgridEdge, rSubgridPoints, sSubgridPoints, &
                                            problem_bathymetry, problem_velocity, problem_ssh, &  
                                            subgridBathymetryValues, subgridAreas, subgridSshValues, subgridUValues, subgridVValues, &
                                            xSubgridCell, ySubgridCell, &
                                            problem_drag, subgridDragValues)
  
      enddo

      bathymetryMin = maxval( subgridBathymetryValues(1:nSubgridEdge) ) 
      bathymetryMax = minval( subgridBathymetryValues(1:nSubgridEdge) )

      if (bathymetryMax < -99998) then
         cycle edge
      endif
      if (bathymetryMin > 99998) then
         cycle edge
      endif

      ! Evaluate velocity average
      !---------------------------------------------------------------
      call ocn_init_grid_average(nSubgridEdge, subgridBathymetryValues, subgridSshValues, subgridAreas, layerThicknessEdgeAverage, &
                                subgridUValues, subgridVValues, uVelocityAverage(iEdge), vVelocityAverage(iEdge))

      if ( .NOT. subgrid_edge_bathymetry_max_pixel) then
         ! Evaluate bounds of look-up table range
         !---------------------------------------------------------------
         call ocn_init_set_table_range(1, nSubgridEdge, subgridBathymetryValues, subgridSshEdgeTableRange(:,iEdge))

         ! Evaluate subgrid bathymetry
         !---------------------------------------------------------------
         subgridEdgeBathymetryMean(iEdge) = sum(subgridBathymetryValues(1:nSubgridEdge))/real(nSubgridEdge,RKIND) 
         subgridEdgeBathymetryMin(iEdge)  = maxval(subgridBathymetryValues(1:nSubgridEdge))

         ! Vertical integration of wet fraction
         !---------------------------------------------------------------
         call ocn_init_vertical_integration(iEdge,subgridSshEdgeTableRange, nSubgridEdge, subgridBathymetryValues, subgridAreas, &
                                         subgridWetVolumeEdgeTable, subgridWetFractionEdgeTable)
      else
         ! DW: Use the higher values of the pair of subcells along the
         !     cell edge
         cellEdgeBathymetryValues(1,:) = subgridBathymetryValues(1:2*nSubgridCellEdge - 1:2) 
         if ( nslice > 1 ) then
           cellEdgeBathymetryValues(2,:) = &
              subgridBathymetryValues(nSubgridTriPerSlice+2*nsubgridCellEdge - 1:nSubgridTriPerSlice+1:-2)
         else
           cellEdgeBathymetryValues(2,:) = cellEdgeBathymetryValues(1,:)     
         endif

         do iEdgeSegment = 1, nSubgridCellEdge
            cellEdgeBathymetryValues(3,iEdgeSegment) = minval(cellEdgeBathymetryValues(1:2,iEdgeSegment) ) 
         end do 
         dsEdge(:) = sqrt( (x(2) - x(1))*(x(2) - x(1)) + (y(2) - y(1))*(y(2)- y(1)) )/nSubgridCellEdge 


         ! Evaluate bounds of look-up table range
         !---------------------------------------------------------------
         call ocn_init_set_table_range(1, nSubgridCellEdge, cellEdgeBathymetryValues(3,:), subgridSshEdgeTableRange(:,iEdge))

         ! Evaluate bounds of look-up table range
         !---------------------------------------------------------------
         subgridEdgeBathymetryMean(iEdge) = sum(cellEdgeBathymetryValues(3,1:nSubgridCellEdge))/real(nSubgridCellEdge,RKIND) 
         subgridEdgeBathymetryMin(iEdge)  = maxval(cellEdgeBathymetryValues(3,1:nSubgridCellEdge))

         ! Vertical integration of wet fraction
         !---------------------------------------------------------------
         call ocn_init_vertical_integration( iEdge, subgridSshEdgeTableRange, &
                    nSubgridCellEdge, cellEdgeBathymetryValues(3,:), dsEdge, subgridWetVolumeEdgeTable, subgridWetFractionEdgeTable )
      endif

      ! Vertical integration of drag
      !---------------------------------------------------------------
      call ocn_init_set_table_range(1, nSubgridEdge, subgridBathymetryValues, subgridSshDragEdgeTableRange(:,iEdge),lower_only=.true.)
      call ocn_init_vertical_integration_drag(subgridSshDragEdgeTableRange(:,iEdge), nSubgridEdge, & 
                                              subgridBathymetryValues, subgridDragValues, subgridAreas, &
                                              subgridDragEdgeTable(:,iEdge), subgridEdgeArea(iEdge))

    end do edge

    ! find an ssh value corresponding to drying_min_cell-height 
    ! of each edge 
    do iEdge = 1,nEdgesSolve
        call ocn_subgrid_ssh_lookup( config_drying_min_cell_height, & 
                                     subgridWetVolumeEdgeTable(:,iEdge),&
                                     subgridSshEdgeTableRange(:,iEdge),&
                                     subgridEdgeBathymetryMean(iEdge),&
                                     subgridEdgeBathymetryMin(iEdge),&
                                     subgridSshEdgeTableRange(3,iEdge) )
    end do

    return

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Vertex
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, "Begin vertex", nVerticesSolve

vertex: do iVertex = 1,nVerticesSolve

      ! Evaluate subgrid bathymetry at centers of sub-triangles for vertex triangle
      !---------------------------------------------------------------

      nSlice = 0
      do slice = 1,3
        if (cellsOnVertex(slice,iVertex) <= nCellsSolve) then
          nSlice = nSlice + 1
        endif
      enddo

      if (nSlice < 3) then
         cycle vertex
      endif

      nSubgridVertex = 0  ! Counter for all subgrid triangles over edge slices
      do slice = 1,nSlice

        c1 = cellsOnVertex(slice,iVertex)
        if (slice < 3) then
          c2 = cellsOnVertex(slice+1,iVertex)  
        else
          c2 = cellsOnVertex(1,iVertex)
        endif

        ! Vertex slice coordinates
        x(1) = xCell(c1)
        y(1) = yCell(c1)

        x(2) = xCell(c2)
        y(2) = yCell(c2)

        x(3) = xVertex(iVertex)
        y(3) = yVertex(iVertex)

        call ocn_init_evaluate_subgrid_data(x, y, nSubgridTriPerSlice, nSubgridVertex, rSubgridPoints, sSubgridPoints, &
                                            problem_bathymetry, problem_velocity, problem_ssh, &  
                                            subgridBathymetryValues, subgridAreas)
      enddo

      bathymetryMin = maxval( subgridBathymetryValues(1:nSubgridVertex) ) 
      bathymetryMax = minval( subgridBathymetryValues(1:nSubgridVertex) )

      if (bathymetryMax < -99998) then
         cycle 
      endif
      if (bathymetryMin > 99998) then
         cycle
      endif

      ! Evaluate bounds of look-up table range
      !---------------------------------------------------------------
      call ocn_init_set_table_range(1, nSubgridVertex, subgridBathymetryValues, subgridSshVertexTableRange(:,iVertex))

      ! Evaluate subgrid bathymetry
      !---------------------------------------------------------------
      subgridVertexBathymetryMean(iVertex) = sum(subgridBathymetryValues(1:nSubgridVertex)*subgridAreas(1:nSubgridVertex))/sum(subgridAreas(1:nSubgridVertex))
      subgridVertexBathymetryMin(iVertex)  = maxval(subgridBathymetryValues(1:nSubgridVertex))

      ! Vertical integration of wet fraction
      !---------------------------------------------------------------
      call ocn_init_vertical_integration(iVertex,subgridSshVertexTableRange, nSubgridVertex, subgridBathymetryValues, subgridAreas, &
                                         subgridWetVolumeVertexTable, subgridWetFractionVertexTable)

    enddo vertex

    ! find an ssh value corresponding to drying_min_cell-height 
    ! of each edge 
    do iVertex = 1,nVerticesSolve

        call ocn_subgrid_ssh_lookup( config_drying_min_cell_height, & 
                                     subgridWetVolumeVertexTable(:,iVertex),&
                                     subgridSshVertexTableRange(:,iVertex),&
                                     subgridVertexBathymetryMean(iVertex),&
                                     subgridVertexBathymetryMin(iVertex),&
                                     subgridSshVertexTableRange(3,iVertex) )
    end do

    return ;  
  end subroutine ocn_init_subgrid_calculations!}}}

!***********************************************************************
!
!  routine ocn_init_set_table_range
!
!> \brief   Set range for subgrid table
!> \author  Steven Brus 
!> \date    March 2025
!> \details Set the upper and lower limits of the subgrid lookup table
!>          range. Default to fully dry/fully wet limits, but restrict
!>          to expected range if default range is too large 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_set_table_range(subgridIdxStart, subgridIdxEnd, subgridBathymetryValues, subgridSshCellTableRange, lower_only) !{{{

      implicit none

      integer, intent(in) :: subgridIdxStart, subgridIdxEnd
      real (kind=RKIND), dimension(:), intent(in) :: subgridBathymetryValues
      real (kind=RKIND), dimension(:), intent(out) :: subgridSshCellTableRange
      logical, intent(in), optional :: lower_only

      real (kind=RKIND) :: SSHMin, SSHMax

   !--------------------------------------------------------------------

      SSHMin = -maxval( subgridBathymetryValues(subgridIdxStart:subgridIdxEnd) ) 
      SSHMax = -minval( subgridBathymetryValues(subgridIdxStart:subgridIdxEnd) )

      if ( abs(SSHMax - SSHMin) > 100.0_RKIND*eps ) then
          subgridSshCellTableRange(1) = SSHMin + subgrid_thin_layer
          subgridSshCellTableRange(2) = SSHMax
      else
          ! flat bathy !
          subgridSshCellTableRange(1) =  SSHMin + subgrid_thin_layer
          subgridSshCellTableRange(2) =  SSHMin + 2.0_RKIND*config_drying_min_cell_height  
      endif 

      if (subgridSshCellTableRange(1) < -8.0_RKIND) then
          subgridSshCellTableRange(1) = -8.0_RKIND
      endif

      if (subgridSshCellTableRange(2) > 8.0_RKIND) then
          subgridSshCellTableRange(2) = 8.0_RKIND
      endif

      if (present(lower_only)) then
         if (lower_only .eqv. .true.) then
            subgridSshCellTableRange(2) = 8.0_RKIND
         endif
      endif

   end subroutine ocn_init_set_table_range!}}}

!***********************************************************************
!
!  routine ocn_init_define_subgrid_points
!
!> \brief   Define subgrid points on reference triangle
!> \author  Steven Brus 
!> \date    November 2022 
!> \details Gives the r and s coordinates of the subgrid triangle vertices
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_define_subgrid_points(nSubgridLevel, rSubgridPoints, sSubgridPoints)!{{{
           
      implicit none 

      integer, intent(in) :: nSubgridLevel
      real (kind=RKIND), dimension(:,:), intent(inout) :: rSubgridPoints, sSubgridPoints

      integer :: i,j,k
      integer :: nSubgridTri
      real (kind=RKIND) :: dx
      real (kind=RKIND), dimension(:), allocatable :: xPoints

   !--------------------------------------------------------------------

      ! Equi-spaced nodes on -1,+1
      allocate(xPoints(nSubgridLevel+1))
      dx = 2.0_RKIND/real(nSubgridLevel,RKIND)
      xPoints(1) = -1.0_RKIND
      do i = 2,nSubgridLevel+1
        xPoints(i) = xPoints(i-1)+dx   
      enddo      

      !   * 
      !   |\
      !   | \
      !   |  \
      !   |   \
      ! s *----*
      !   |\ u |\
      !   | \  | \
      !   |  \ |  \
      !   | l \| l \
      !   *----*----*
      !        r

      ! Trianglulate tensor product of equi-spaced nodes on reference triangle
      ! r - horizontal coordinate (i index), s - vertical coordinate (j index)
      nSubgridTri = 0
      do j = 1,nSubgridLevel
        do i = 1,nSubgridLevel+1 - j 

          nSubgridTri = nSubgridTri+1 

          ! lower triangle in pair (see l in triangle above)
          rSubgridPoints(1,nSubgridTri) = xPoints(i)
          sSubgridPoints(1,nSubgridTri) = xPoints(j)
          
          rSubgridPoints(2,nSubgridTri) = xPoints(i+1)
          sSubgridPoints(2,nSubgridTri) = xPoints(j)

          rSubgridPoints(3,nSubgridTri) = xPoints(i)
          sSubgridPoints(3,nSubgridTri) = xPoints(j+1)

          ! upper triangle in pair (see u in triangle above, doesn't occur next to hypotenuse)          
          if (i < nSubgridLevel+1 - j) then
            
            nSubgridTri = nSubgridTri + 1

            rSubgridPoints(1,nSubgridTri) = xPoints(i+1)
            sSubgridPoints(1,nSubgridTri) = xPoints(j)
            
            rSubgridPoints(2,nSubgridTri) = xPoints(i+1)
            sSubgridPoints(2,nSubgridTri) = xPoints(j+1)

            rSubgridPoints(3,nSubgridTri) = xPoints(i)
            sSubgridPoints(3,nSubgridTri) = xPoints(j+1)

          endif

        enddo
      enddo

      print*, "exit ocn_init_define_subgrid_points"

   
   !--------------------------------------------------------------------

   end subroutine ocn_init_define_subgrid_points!}}}

!***********************************************************************
!
!  routine ocn_init_evaluate_subgrid_data
!
!> \brief   Evaluate subgrid infromation 
!> \author  Steven Brus
!> \date    November 2022
!> \details Evaluate data on subgrid triangles for a given triangluar grid cell region 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_evaluate_subgrid_data(xTri, yTri, nSubgridTri, nSubgridCV, rSubgridPoints, sSubgridPoints, & 
                                             problem_bathymetry, problem_velocity, problem_ssh, &
                                             subgridBathymetryValues, subgridAreas, &
                                             subgridSshValues, subgridUValues, subgridVValues, &
                                             xSubgridCell, ySubgridCell, &
                                             problem_drag, subgridDragValues) !{{{

           
      implicit none 
      real (kind=RKIND), intent(inout) :: xTri(3), yTri(3)
      integer, intent(in) :: nSubgridTri
      integer, intent(inout) :: nSubgridCV

      procedure(bathymetry_function) :: problem_bathymetry
      procedure(velocity_function) :: problem_velocity
      procedure(ssh_function) :: problem_ssh

      real (kind=RKIND), dimension(:,:), intent(in) :: rSubgridPoints, sSubgridPoints
      real (kind=RKIND), dimension(:), intent(inout) :: subgridBathymetryValues
      real (kind=RKIND), dimension(:), intent(inout) :: subgridAreas
      real (kind=RKIND), dimension(:), intent(inout),optional :: subgridSshValues
      real (kind=RKIND), dimension(:), intent(inout),optional :: subgridUValues, subgridVValues

      real (kind=RKIND), dimension(:), intent(inout), optional :: xSubgridCell, ySubgridCell
      procedure(bathymetry_function), optional :: problem_drag
      real (kind=RKIND), dimension(:), intent(inout), optional :: subgridDragValues

      real (kind=RKIND) :: rCenter, sCenter
      real (kind=RKIND) :: xCenter, yCenter
      real (kind=RKIND) :: xSubgridPoints(3), ySubgridPoints(3)
      real (kind=RKIND) :: x, y
      real (kind=RKIND) :: area


      integer :: iPt,i

   !--------------------------------------------------------------------

     ! Coordinates of physical triangle
     ! (ensure counter-clockwise numering)
     call ocn_init_tri_area(xTri(:), yTri(:), area)
     if (area < 0.0_RKIND) then
      x = xTri(1)
      y = yTri(1) 
      xTri(1) = xTri(2) 
      yTri(1) = yTri(2)
      xTri(2) = x ;
      yTri(2) = y ; 
      area = abs(area) ; 
     endif

      do iPt = 1,nSubgridTri

        ! Counter over all subcells within cell/edge/vertex control volume
        nSubgridCV = nSubgridCV + 1

        ! Center sub-triangle (on reference triangle)
        rCenter = sum(rSubgridPoints(:,iPt))/3.0_RKIND
        sCenter = sum(sSubgridPoints(:,iPt))/3.0_RKIND

        ! Transformation of sub-triangle center to physical coordinates
        call ocn_init_tri_coordinate_transform(rCenter, sCenter, xTri, yTri, xCenter, yCenter)

        ! Evaluate bathymetry
        call problem_bathymetry(xCenter, yCenter, subgridBathymetryValues(nSubgridCV))

        ! Evaluate drag
        if (present(problem_drag).and.present(subgridDragValues)) then
           call problem_drag(xCenter, yCenter, subgridDragValues(nSubgridCV))
        endif

        ! Transformation of sub-triangle vertices to physical coordinates
        do i = 1,3
           call ocn_init_tri_coordinate_transform(rSubgridPoints(i,iPt), sSubgridPoints(i,iPt), xTri, yTri, xSubgridPoints(i), ySubgridPoints(i))
        enddo
 
        ! Calculate area of sub-triangle
        call ocn_init_tri_area(xSubgridPoints(:), ySubgridPoints(:), subgridAreas(nSubgridCV))

        ! Optionally evalulate ssh !
        if (present(subgridsshValues)) then
           call problem_ssh(xCenter, yCenter, subgridBathymetryValues(nSubgridCV), subgridSshValues(nSubgridCV))
        endif

        if (present(subgridUValues).and.present(subgridVValues)) then
           call problem_velocity(xCenter, yCenter, subgridUValues(nSubgridCV), subgridVValues(nSubgridCV))
        endif   

        if ( present(xSubgridCell) ) then
            xSubgridCell(nSubgridCV) = xCenter ;
        endif

        if ( present(ySubgridCell) ) then
            ySubgridCell(nSubgridCV) = yCenter ;
        end if
      enddo

   !--------------------------------------------------------------------

   end subroutine ocn_init_evaluate_subgrid_data!}}}

!***********************************************************************
!
!  routine ocn_init_tri_area
!
!> \brief   Compute triangle area
!> \author  Steven Brus
!> \date    November 2022
!> \details Compute area of triangle given vertex coordinates
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_tri_area(x, y, area)!{{{
           
      implicit none 

      real (kind=RKIND), intent(in) :: x(3), y(3)
      real (kind=RKIND), intent(out) :: area

   !--------------------------------------------------------------------

      area = 0.5_RKIND*(x(2)*y(3) - y(2)*x(3) - x(1)*y(3) + y(1)*x(3) + x(1)*y(2) - y(1)*x(2))

      !xv = cosd( lat ).*cosd( lon ) ;
      !yv = cosd( lat ).*sind( lon ) ;
      !zv = sind( lat ) ;

      !avec = [ xv(:,1) yv(:,1) zv(:,1) ] ;
      !bvec = [ xv(:,2) yv(:,2) zv(:,2) ] ;
      !cvec = [ xv(:,3) yv(:,3) zv(:,3) ] ;

      !axb = cross_product( avec, bvec ) ;
      !cp = atan( sqrt(sum(axb.*axb ,2))./sum(avec.*bvec, 2)) ;

      !axb = cross_product( avec, cvec ) ;
      !bp = atan( sqrt(sum(axb.*axb,2))./sum(avec.*cvec,2 ) ) ;

      !axb = cross_product( bvec, cvec ) ;
      !ap = atan( sqrt(sum(axb.*axb,2))./sum(bvec.*cvec,2) ) ;


      !A = (cos(ap) - cos(bp).*cos(cp))./( sin(bp).*sin(cp) ) ;
      !B = (cos(bp) - cos(cp).*cos(ap))./( sin(cp).*sin(ap) ) ;
      !C = (cos(cp) - cos(ap).*cos(bp))./( sin(ap).*sin(bp) ) ;

      !A = acos( A ) ;
      !B = acos( B ) ;
      !C = acos( C ) ;

      !% 3. Return area
      !dA = A + B + C - pi ;

   !--------------------------------------------------------------------

   end subroutine ocn_init_tri_area!}}}

   function cross_product(a, b) result(c)
       implicit none
       real, intent(in) :: a(3), b(3)
       real :: c(3)
   
       c(1) = a(2) * b(3) - a(3) * b(2)
       c(2) = a(3) * b(1) - a(1) * b(3)
       c(3) = a(1) * b(2) - a(2) * b(1)
   
   end function cross_product

!***********************************************************************
!
!  routine ocn_init_tri_coordinate_transform
!
!> \brief   Transform reference triangle coordinates to mesh coordinates
!> \author  Steven Brus
!> \date    November 2022
!> \details Evaluate r,s reference coordinates in x,y mesh
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_tri_coordinate_transform(r, s, xTri, yTri, x, y)!{{{
           
      implicit none 

      real (kind=RKIND), intent(in) :: r, s
      real (kind=RKIND), intent(in) :: xTri(3), yTri(3)
      real (kind=RKIND), intent(out) :: x, y

   !--------------------------------------------------------------------

      ! Transformation of sub-triangle center to physical coordinates
      x = 0.5_RKIND*(-(r+s)*xTri(1) + (1.0_RKIND+r)*xTri(2) + (1.0_RKIND+s)*xTri(3))
      y = 0.5_RKIND*(-(r+s)*yTri(1) + (1.0_RKIND+r)*yTri(2) + (1.0_RKIND+s)*yTri(3))

   !--------------------------------------------------------------------

   end subroutine ocn_init_tri_coordinate_transform!}}}

!***********************************************************************
!
!  routine ocn_init_vertical_integration
!
!> \brief   Compute the wet volume per unit area lookup table
!> \author  Steven Brus
!> \date    November 2022
!> \details Integrate the wet fraction over discrete ssh values 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_vertical_integration(iCV, subgridSshTableRange, nSubgridCV, subgridBathymetryValues, subgridAreas, &
                                            subgridWetVolumeTable, subgridWetFractionTable)!{{{
           
      implicit none 

      integer, intent(in) :: iCV
      real (kind=RKIND), dimension(:,:), intent(in) :: subgridSshTableRange
      integer, intent(in) :: nSubgridCV
      real (kind=RKIND), dimension(:), intent(in) :: subgridBathymetryValues, subgridAreas
      real (kind=RKIND), dimension(:,:), intent(inout) :: subgridWetVolumeTable
      real (kind=RKIND), dimension(:,:), intent(inout) :: subgridWetFractionTable

      real (kind=RKIND) :: deltaZ, ssh, pVal
      integer :: lev, tri
      

   !--------------------------------------------------------------------
   
      deltaZ = (subgridSshTableRange(2,iCV)-subgridSshTableRange(1,iCV))/real(nSubgridTableLevels-1,RKIND);

      subgridWetVolumeTable(1,iCV) = subgrid_thin_layer
      ssh = subgridSshTableRange(1,iCV) + deltaZ
      subgridWetFractionTable(1,iCV) = 0.0_RKIND
      do lev = 2,nSubgridTableLevels

        subgridWetFractionTable(lev,iCV) = 0.0_RKIND
        subgridWetVolumeTable(lev,iCV) = 0.0_RKIND

        do tri = 1,nSubgridCV

          pVal = 0.0_RKIND
          if (subgridBathymetryValues(tri) + ssh > subgrid_thin_layer) then
            pVal = 1.0_RKIND
          endif

          subgridWetFractionTable(lev,iCV) = subgridWetFractionTable(lev,iCV) + pVal*subgridAreas(tri)
          subgridWetVolumeTable(lev,iCV) = subgridWetVolumeTable(lev,iCV) + max(ssh+subgridBathymetryValues(tri),subgrid_thin_layer)*subgridAreas(tri)
  
        enddo

        subgridWetFractionTable(lev,iCV) = subgridWetFractionTable(lev,iCV)/sum(subgridAreas(1:nSubgridCV))
        subgridWetVolumeTable(lev,iCV) = subgridWetVolumeTable(lev,iCV)/sum(subgridAreas(1:nSubgridCV))

        ssh = ssh + deltaZ
      enddo

   !--------------------------------------------------------------------

   end subroutine ocn_init_vertical_integration!}}}

   !***********************************************************************
!
!  routine ocn_init_vertical_integration_drag
!
!> \brief   Compute the wet volume per unit area lookup table
!> \author  Steven Brus
!> \date    November 2022
!> \details Integrate the wet fraction over discrete ssh values 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_vertical_integration_drag(subgridSshTableRange, nSubgridCV, subgridBathymetryValues, subgridManningsValues, subgridAreas, &
                                                 subgridDragTable, subgridEdgeArea)!{{{
           
      implicit none 

      real (kind=RKIND), dimension(:), intent(in) :: subgridSshTableRange
      integer, intent(in) :: nSubgridCV
      real (kind=RKIND), dimension(:), intent(in) :: subgridBathymetryValues, subgridManningsValues, subgridAreas
      real (kind=RKIND), dimension(:), intent(inout) :: subgridDragTable
      real (kind=RKIND), intent(inout) :: subgridEdgeArea

      real (kind=RKIND) :: deltaZ, ssh, H
      real (kind=RKIND) :: num, den
      integer :: lev, tri
      
   !--------------------------------------------------------------------
   
      deltaZ = (subgridSshTableRange(2)-subgridSshTableRange(1))/real(nSubgridTableLevels-1,RKIND);

      ssh = subgridSshTableRange(1)
      do lev = 1,nSubgridTableLevels

        num = 0.0_RKIND
        den = 0.0_RKIND
        do tri = 1,nSubgridCV

          H = max(ssh+subgridBathymetryValues(tri),subgrid_thin_layer) 
          if (H > 0.0_RKIND) then
             num = num + H*subgridAreas(tri)
             den = den + H**(5.0_RKIND/3.0_RKIND)/sqrt(gravity)/subgridManningsValues(tri)*subgridAreas(tri)
          endif

        enddo

        subgridEdgeArea = sum(subgridAreas(1:nSubgridCV))
        if (den > 0.0_RKIND) then
           subgridDragTable(lev) = (num**3/den**2)/subgridEdgeArea
        else
           subgridDragTable(lev) = 0.0_RKIND
        endif

        ssh = ssh + deltaZ
      enddo

   !--------------------------------------------------------------------

   end subroutine ocn_init_vertical_integration_drag!}}}

!***********************************************************************
!
!  routine ocn_init_grid_average
!
!> \brief   Compute thickness and velocity averages over wet area 
!> \author  Steven Brus
!> \date    November 2022
!> \details Compute the thickness and velocity averages over an area
!>          based on the subgrid wet area
!
!-----------------------------------------------------------------------

   subroutine ocn_init_grid_average(nSubgridCV, subgridBathymetryValues, subgridSshValues, subgridAreas, subgridThicknessAverage, &
                                   subgridUValues, subgridVValues, subgridUAverage, subgridVAverage)!{{{
           
      implicit none 

      integer, intent(in) :: nSubgridCV
      real (kind=RKIND), dimension(:), intent(in) :: subgridBathymetryValues, subgridAreas
      real (kind=RKIND), dimension(:), intent(in) :: subgridSshValues
      real (kind=RKIND), intent(inout) :: subgridThicknessAverage
      real (kind=RKIND), dimension(:), intent(in), optional :: subgridUValues, subgridVValues
      real (kind=RKIND), intent(inout), optional :: subgridUAverage, subgridVAverage

      real (kind=RKIND) :: deltaZ, ssh, pVal
      real (kind=RKIND) :: averageDepth, cellArea, layerThicknessValue
      real (kind=RKIND) :: H_int 
      real (kind=RKIND) :: HU_int, HV_int
      logical :: computeVelAverage
      integer :: lev, tri
      

   !--------------------------------------------------------------------
   
      averageDepth = sum(subgridBathymetryValues(1:nSubgridCV))/real(nSubgridCV,RKIND) 

      computeVelAverage = .false.
      if (present(subgridUValues) .and. present(subgridVValues) .and. &
          present(subgridUAverage) .and. present(subgridVAverage)) then
         computeVelAverage = .true.
      endif

      H_int = 0.0_RKIND
      cellArea = 0.0_RKIND
      HU_int = 0.0_RKIND
      HV_int = 0.0_RKIND
      do tri = 1,nSubgridCV

        layerThicknessValue = subgridBathymetryValues(tri) + subgridSshValues(tri)
        if (layerThicknessValue <= config_drying_min_cell_height + eps) then
          layerThicknessValue = config_drying_min_cell_height + eps
        endif
        H_int = H_int + layerThicknessValue*subgridAreas(tri)
        cellArea = cellArea + subgridAreas(tri)

        if (computeVelAverage) then
           HU_int = HU_int + layerThicknessValue*subgridUValues(tri)*subgridAreas(tri)
           HV_int = HV_int + layerThicknessValue*subgridVValues(tri)*subgridAreas(tri)
        endif
      enddo

      if (computeVelAverage) then
         subgridUAverage = HU_int/H_int
         subgridVAverage = HV_int/H_int
      endif
      subgridThicknessAverage = H_int/cellArea
      

   !--------------------------------------------------------------------

   end subroutine ocn_init_grid_average!}}}

!***********************************************************************
!
!  routine ocn_init_wet_average_ssh
!
!> \brief   
!> \author  Steven Brus
!> \date    November 2022
!> \details 
!>  
!
!-----------------------------------------------------------------------
   subroutine ocn_init_wet_average_ssh( nSubgridCV, subgridBathymetryValues, subgridSshValues, subgridAreas, sshWetAverage)!{{{
      implicit none

      integer, intent(in) :: nSubgridCV
      real (kind=RKIND), dimension(:), intent(in) :: subgridBathymetryValues, subgridAreas
      real (kind=RKIND), dimension(:), intent(in) :: subgridSshValues
      real (kind=RKIND), intent(inout) :: sshWetAverage

      real (kind=RKIND) :: deltaZ, ssh, pVal
      real (kind=RKIND) :: averageDepth, wetArea, layerThicknessValue
      integer :: lev, tri


      sshWetAverage = 0.0_RKIND 
      wetArea = 0.0_RKIND 

      do tri = 1, nsubgridCV
        layerThicknessValue = subgridBathymetryValues(tri) + subgridSshValues(tri)

        if ( layerThicknessValue > subgrid_thin_layer) then 
           sshWetAverage = sshWetAverage + subgridSshValues(tri)*subgridAreas(tri)
           wetArea = wetArea + subgridAreas(tri)
        endif
      end do

      if ( WetArea > 0.0_RKIND ) then
          sshWetAverage = sshWetAverage/WetArea ;
      else
         sshWetAverage = -maxval(subgridBathymetryValues(1:nsubgridCV)) + subgrid_thin_layer
      endif

      return ;
   end subroutine ocn_init_wet_average_ssh

!***********************************************************************
!
!  routine ocn_init_read_subgrid_data
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
    subroutine ocn_init_read_subgrid_data(domain, topoLon, topoLat, topoIC, &
                                          config_topography_file, config_topography_source, config_topography_latlon_degrees, config_topography_positive, &
                                          config_topography_lon_varname, config_topography_lat_varname, config_topography_nlon_dimname, config_topography_nlat_dimname, config_topography_varname, &
                                          iErr)!{{{
       type (domain_type), intent(inout) :: domain
       character (len=*), intent(in) :: config_topography_file
       character (len=*), intent(in) :: config_topography_source
       logical, intent(in) :: config_topography_latlon_degrees
       character (len=*), intent(in) :: config_topography_positive
       character (len=*), intent(in) :: config_topography_lon_varname
       character (len=*), intent(in) :: config_topography_lat_varname
       character (len=*), intent(in) :: config_topography_nlon_dimname
       character (len=*), intent(in) :: config_topography_nlat_dimname
       character (len=*), intent(in) :: config_topography_varname
       type (field1DReal), intent(inout) :: topoLon, topoLat
       type (field2DReal), intent(inout) :: topoIC
       integer, intent(out) :: iErr

       integer :: i
       integer :: nLonTopo, nLatTopo
       type (MPAS_Stream_type) :: topographyStream
       type (MPAS_IO_Handle_type) :: inputFile
       type (MPAS_IO_Context_type), pointer :: iocontext_ptr
       
       iErr = 0

       if (config_topography_source == 'latlon_file' .or. &
           config_topography_source == 'xy_file'  ) then

         inputFile = MPAS_io_open(trim(config_topography_file), &
                        MPAS_IO_READ, MPAS_IO_NETCDF, domain % iocontext, ierr=iErr)
         if (iErr /= 0) then
            call mpas_log_write( 'could not open file '// &
                   trim(config_topography_file), MPAS_LOG_CRIT)
            return
         end if

         call MPAS_io_inq_dim(inputFile, config_topography_nlat_dimname, nLatTopo, iErr)
         call MPAS_io_inq_dim(inputFile, config_topography_nlon_dimname, nLonTopo, iErr)
         call mpas_log_write('Topography file '//trim(config_topography_file))
         call mpas_log_write(' nLat = $i, nLon = $i', intArgs=(/nLatTopo, nLonTopo/))

         call MPAS_io_close(inputFile, iErr)

       end if

       ! Define stream for depth levels
       call MPAS_createStream(topographyStream, domain % iocontext, &
                              trim(config_topography_file), MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)
       call mpas_log_write('createStream Error = $i', intArgs=(/iErr/))

       ! Setup topoLat, topoLon, and topoIC fields for stream to be read in
       topoLat % fieldName = trim(config_topography_lat_varname)
       topoLat % dimSizes(1) = nLatTopo
       topoLat % dimNames(1) = trim(config_topography_nlat_dimname)
       topoLat % isVarArray = .false.
       topoLat % isPersistent = .true.
       topoLat % isActive = .true.
       topoLat % hasTimeDimension = .false.
       topoLat % block => domain % blocklist
       allocate(topoLat % attLists(1))
       allocate(topoLat % array(nLatTopo))

       topoLon % fieldName = trim(config_topography_lon_varname)
       topoLon % dimSizes(1) = nLonTopo
       topoLon % dimNames(1) = trim(config_topography_nlon_dimname)
       topoLon % isVarArray = .false.
       topoLon % isPersistent = .true.
       topoLon % isActive = .true.
       topoLon % hasTimeDimension = .false.
       topoLon % block => domain % blocklist
       allocate(topoLon % attLists(1))
       allocate(topoLon % array(nLonTopo))

       topoIC % fieldName = trim(config_topography_varname)
       topoIC % dimSizes(1) = nLonTopo
       topoIC % dimSizes(2) = nLatTopo
       topoIC % dimNames(1) = trim(config_topography_nlon_dimname)
       topoIC % dimNames(2) = trim(config_topography_nlat_dimname)
       topoIC % isVarArray = .false.
       topoIC % isPersistent = .true.
       topoIC % isActive = .true.
       topoIC % hasTimeDimension = .false.
       topoIC % block => domain % blocklist
       allocate(topoIC % attLists(1))
       allocate(topoIC % array(nLonTopo, nLatTopo))


       ! Add topoLat, topoLon, and topoIC fields to stream
       call MPAS_streamAddField(topographyStream, topoLat, iErr)
       call mpas_log_write('streamAddField lat Error = $i', intArgs=(/iErr/))
       call MPAS_streamAddField(topographyStream, topoLon, iErr)
       call mpas_log_write('streamAddField lon Error = $i', intArgs=(/iErr/))
       call MPAS_streamAddField(topographyStream, topoIC, iErr)
       call mpas_log_write('streamAddField IC Error = $i', intArgs=(/iErr/))

       ! Read stream
       call MPAS_readStream(topographyStream, 1, iErr)
       if ( adjustl(trim(config_topography_positive)) /= 'down' ) then
           topoIC%array = -topoIC%array ;
       end if


       ! Close stream
       call MPAS_closeStream(topographyStream)

       if ( config_topography_latlon_degrees .and. &
              config_topography_source == 'latlon_file' ) then
          !do i = 1,nLonTopo
          !   if (topoLon % array(i) < 0.0) then
          !      topoLon % array(i) = topoLon % array(i) + 360.0_RKIND
          !   end if
          !end do

          topoLat % array(:) = topoLat % array(:) * pii / 180.0_RKIND
          topoLon % array(:) = topoLon % array(:) * pii / 180.0_RKIND
       end if

    end subroutine ocn_init_read_subgrid_data!}}}

!***********************************************************************
!
!  routine ocn_init_subgrid_interpolate_topo
!
!> \brief   Interpolate bathymetry
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Return the value of the bathymetry at a given x,y point
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_init_subgrid_interpolate_topo(xin, yin, depth)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: xin, yin
      real (kind=RKIND), intent(out) :: depth

      integer :: ix(2)
      integer :: i, idx_lon, idx_lat
      integer :: nLonTopo, nLatTopo

      real (kind=RKIND) :: x, y
      real (kind=RKIND) :: xc(2)
      real (kind=RKIND) :: dx(2), bath(4), f(4), den, d(4)
      real (kind=RKIND) :: xv(2), yv(2)

      !--------------------------------------------------------------------
      !
      ! DEM must be on a uniform grid

      y = yin 
      x = xin

      if (config_subgrid_topography_latlon_degrees .and. &
          config_subgrid_topography_source == 'latlon_file' ) then

          if (x > pii) then
              x = x - 2.0_RKIND*pii
          end if

      end if

      nLonTopo = size(topoLon%array)
      nLatTopo = size(topoLat%array)

      ! Bilinear interpolation
      xc = (/ x, y /) ;

      ix(1) = bisection_index(x, topoLon%array, nLonTopo)
      ix(2) = bisection_index(y, topoLat%array, nLatTopo)

      !PRINT*, ix(1), ix(2)
      if ((ix(1) > 0) .and. (ix(2) > 0)) then

        dx = (/ topoLon%array(ix(1)+1) - topoLon%array(ix(1)), &
                topoLat%array(ix(2)+1) - topoLat%array(ix(2)) /) ;

        bath(1) = TopoData%array(ix(1),ix(2))  ;
        bath(2) = TopoData%array(ix(1)+1,ix(2)) ;
        bath(3) = TopoData%array(ix(1)+1,ix(2)+1) ;
        bath(4) = TopoData%array(ix(1),ix(2)+1) ;

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

        if ((bath(1) /= bath(1)) .or. (bath(2) /= bath(2)) .or. (bath(3) /= bath(3)) .or. (bath(4) /= bath(4))) then
            call extrapolate_topography( )
            return
        end if

        if ( minval(bath) < -(abs(config_subgrid_topography_fill_value) - 1)  ) then
          call extrapolate_topography( ) ;              
        end if    
        if ( minval(bath) >  (abs(config_subgrid_topography_fill_value) - 1)  ) then
           call extrapolate_topography( ) ;
        end if    

      else

         call extrapolate_topography( ) ;

      endif

      return ;
   contains
     
      subroutine extrapolate_topography(  )
        implicit none 

        ! depth = -1e10
        if ( adjustl(trim(config_subgrid_topography_extrapolation_method)) == 'constant' ) then
          depth = config_subgrid_topography_extrapolation_constant_value ;       
        else
          !! nearest extrapolation !
          ix(1) = merge( 1, ix(1), ix(1) < 1 ) ;
          ix(1) = merge( nLonTopo, ix(1), ix(1) >= nLonTopo ) ;
          ix(2) = merge( 1, ix(2), ix(2) < 1 )  ;
          ix(2) = merge( nLatTopo, ix(2), ix(2) >= nLatTopo ) ;
          depth = TopoData%array(ix(1),ix(2)) ;
        end if
           
      end subroutine extrapolate_topography        

   end subroutine ocn_init_subgrid_interpolate_topo!}}}

!***********************************************************************
!
!  routine ocn_init_subgrid_interpolate_drag
!
!> \brief   Interpolate bathymetry
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Return the value of the drag coefficient at a given x,y point
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_init_subgrid_interpolate_drag(xin, yin, drag)!{{{

      implicit none

      real (kind=RKIND), intent(in) :: xin, yin
      real (kind=RKIND), intent(out) :: drag

      integer :: ix(2)
      integer :: i, idx_lon, idx_lat
      integer :: min_idx
      integer :: nLonDrag, nLatDrag

      real (kind=RKIND) :: x, y
      real (kind=RKIND) :: xc(2)
      integer :: f(4)
      real (kind=RKIND) :: d(4)
      real (kind=RKIND) :: xv(2), yv(2)

      !--------------------------------------------------------------------
      !
      ! DEM must be on a uniform grid

      y = yin 
      x = xin

      if (config_subgrid_lulc_latlon_degrees .and. &
          config_subgrid_lulc_source == 'latlon_file' ) then

          if (x > pii) then
              x = x - 2.0_RKIND*pii
          end if

      end if

      nLonDrag = size(lulcLon%array)
      nLatDrag = size(lulcLat%array)

      xc = (/ x, y /) ;

      ix(1) = bisection_index(x, lulcLon%array, nLonDrag)
      ix(2) = bisection_index(y, lulcLat%array, nLatDrag)

      !PRINT*, ix(1), ix(2)
      if ((ix(1) > 0) .and. (ix(2) > 0)) then

        f(1) = lulcData%array(ix(1),ix(2))  ;
        f(2) = lulcData%array(ix(1)+1,ix(2)) ;
        f(3) = lulcData%array(ix(1)+1,ix(2)+1) ;
        f(4) = lulcData%array(ix(1),ix(2)+1) ;

        xv(1) = lulcLon%array(ix(1)) ;
        xv(2) = lulcLon%array(ix(1)+1) ;

        yv(1) = lulcLat%array(ix(2)) ;
        yv(2) = lulcLat%array(ix(2)+1) ;

        d(1) = ( xc(1) - xv(1) )**2 + ( xc(2) - yv(1) )**2
        d(2) = ( xc(1) - xv(2) )**2 + ( xc(2) - yv(1) )**2
        d(3) = ( xc(1) - xv(2) )**2 + ( xc(2) - yv(2) )**2
        d(4) = ( xc(1) - xv(1) )**2 + ( xc(2) - yv(2) )**2

        min_idx = minloc(d, dim=1)

        call ocn_init_lulc_mannings_mapping(f(min_idx), drag)
      else

        call ocn_init_lulc_mannings_mapping(1, drag)

      endif
      

   end subroutine ocn_init_subgrid_interpolate_drag!}}}

!***********************************************************************
!
!  routine ocn_init_lulc_mannings_mapping
!
!> \brief   Map land use/land cover to Manning n values
!> \author  D. Wirasaet, S. Brus
!> \date    Sep 2025
!> \details
!>  This routine maps land use/land cover types to
!>  appropriate Manning n values
!>
!-----------------------------------------------------------------------

subroutine ocn_init_lulc_mannings_mapping(lulc_type, manning_n)!{{{

   implicit none

   ! Local variables
   integer, intent(in) :: lulc_type
   real (kind=RKIND), intent(out) :: manning_n

      select case (lulc_type)
      case (2)
         manning_n = 0.120_RKIND   ! High Intensity Developed
      case (3)
         manning_n = 0.100_RKIND   ! Medium Intensity Developed
      case (4)
         manning_n = 0.070_RKIND   ! Low Intensity Developed
      case (5)
         manning_n = 0.035_RKIND   ! Developed Open Space
      case (6)
         manning_n = 0.100_RKIND   ! Cultivated Land
      case (7)
         manning_n = 0.055_RKIND   ! Pasture/Hay
      case (8)
         manning_n = 0.035_RKIND   ! Grassland
      case (9)
         manning_n = 0.160_RKIND   ! Deciduous Forest
      case (10)
         manning_n = 0.180_RKIND   ! Evergreen Forest
      case (11)
         manning_n = 0.170_RKIND   ! Mixed Forest
      case (12)
         manning_n = 0.080_RKIND   ! Scrub/Shrub
      case (13)
         manning_n = 0.150_RKIND   ! Palustrine Forested Wetland
      case (14)
         manning_n = 0.075_RKIND   ! Palustrine Scrub/Shrub Wetland
      case (15)
         manning_n = 0.070_RKIND   ! Palustrine Emergent Wetland
      case (16)
         manning_n = 0.150_RKIND   ! Estuarine Forested Wetland
      case (17)
         manning_n = 0.070_RKIND   ! Estuarine Scrub/Shrub Wetland
      case (18)
         manning_n = 0.050_RKIND   ! Estuarine Emergent Wetland
      case (19)
         manning_n = 0.030_RKIND   ! Unconsolidated Shore
      case (20)
         manning_n = 0.030_RKIND   ! Bare Land
      case (21)
         manning_n = 0.025_RKIND   ! Open Water
      case (22)
         manning_n = 0.035_RKIND   ! Palustrine Aquatic Bed
      case (23)
         manning_n = 0.030_RKIND   ! Estuarine Aquatic Bed
      case default
         manning_n = 0.025_RKIND   ! Default value
      end select

end subroutine ocn_init_lulc_mannings_mapping!}}}

end module ocn_init_subgrid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
