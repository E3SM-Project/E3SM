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

   public :: ocn_init_subgrid_calculations 

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

  subroutine ocn_init_subgrid_calculations(domain, problem_bathymetry, problem_velocity, problem_ssh, &
                                           subgrid_refinement_level, &
                                           subgrid_edge_bathymetry_max_pixel, &
                                           subgrid_use_thin_layer, &
                                           uVelocityAverage, vVelocityAverage, iErr)!{{{
    
    implicit none
    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr
    procedure(bathymetry_function) :: problem_bathymetry
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
    integer :: nSubgridCell, nSubgridEdge, nSubgridVertex
    integer :: nSubgridTriPerSlice
    integer :: v1, v2
    integer :: c1, c2
    real (kind=RKIND) :: x(3), y(3)
    integer :: slice, nSlice
    real (kind=RKIND) :: deltaZ
    logical, pointer :: on_a_sphere



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

    nSubgridTriPerSlice = subgrid_refinement_level**2
    allocate(rSubgridPoints(3,maxEdges*nSubgridTriPerSlice), sSubgridPoints(3,maxEdges*nSubgridTriPerSlice))
    call ocn_init_define_subgrid_points(subgrid_refinement_level, rSubgridPoints, sSubgridPoints)
    allocate(subgridBathymetryValues(maxEdges*nSubgridTriPerSlice), subgridAreas(maxEdges*nSubgridTriPerSlice))
    allocate(subgridSshValues(maxEdges*nSubgridTriPerSlice))
    allocate(subgridUValues(maxEdges*nSubgridTriPerSlice))
    allocate(subgridVValues(maxEdges*nSubgridTriPerSlice))

    allocate(xSubgridCell(maxEdges*nSubgridTriPerSlice)) 
    allocate(ySubgridCell(maxEdges*nSubgridTriPErSlice)) 
    xSubgridCell = 0.0_RKIND  
    ySubgridCell = 0.0_RKIND 

   do iCell = 1,nCellsSolve
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
                                            subgridUValues, subgridVValues, xSubgridCell, ysubgridCell) !{{{
      enddo
    
      ! Evaluate bounds of look-up table range
      !---------------------------------------------------------------        
      bathymetryMin = maxval( subgridBathymetryValues(1:nSubgridCell) ) 
      bathymetryMax = minval( subgridBathymetryValues(1:nSubgridCell) )

      if (bathymetryMax < -99998) then
         cycle 
      endif
      if (bathymetryMin > 99998) then
         cycle
      endif

      if ( abs(bathymetryMin - bathymetryMax) > 100.0*eps ) then
          subgridSshCellTableRange(1,iCell) = -maxval(subgridBathymetryValues(1:nSubgridCell)) + subgrid_thin_layer
          subgridSshCellTableRange(2,iCell) = -minval(subgridBathymetryValues(1:nSubgridCell))
      else
          ! flat bathy !
          subgridSshCellTableRange(1,iCell) =  -bathymetryMin + subgrid_thin_layer
          subgridSshCellTableRange(2,iCell) =  -bathymetryMin + 2.0*config_drying_min_cell_height  
      endif 

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

        call ocn_subgrid_ssh_lookup( config_drying_min_cell_height, & 
                                     subgridWetVolumeCellTable(:,iCell), &
                                     subgridSshCellTableRange(:,iCell),&
                                     bottomDepth(iCell),&
                                     subgridCellBathymetryMin(iCell),&
                                     subgridSShCellTableRange(3,iCell) )

       if ( ssh(iCell) < subgridSshCellTableRange(3,iCell) ) then
          ssh(iCell) = subgridSshCellTableRange(3,iCell) ;  
       end if 
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Edges
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, "Begin edges", nEdgesSolve

    nSubgridCellEdge=subgrid_refinement_level 
    allocate( cellEdgeBathymetryValues(3,nSubgridCellEdge) )   
    allocate( dsEdge(nSubgridCellEdge) ) 
    cellEdgeBathymetryValues = -99999 ; 


    do iEdge = 1,nEdgesSolve
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
                                            xSubgridCell, ySubgridCell)
  
      enddo

      bathymetryMin = maxval( subgridBathymetryValues(1:nSubgridEdge) ) 
      bathymetryMax = minval( subgridBathymetryValues(1:nSubgridEdge) )

      if (bathymetryMax < -99998) then
         cycle 
      endif
      if (bathymetryMin > 99998) then
         cycle
      endif

      ! Evaluate velocity average
      !---------------------------------------------------------------
      call ocn_init_grid_average(nSubgridCell, subgridBathymetryValues, subgridSshValues, subgridAreas, layerThicknessEdgeAverage, &
                                subgridUValues, subgridVValues, uVelocityAverage(iEdge), vVelocityAverage(iEdge))

      if ( .NOT. subgrid_edge_bathymetry_max_pixel) then
         ! Evaluate bounds of look-up table range
         !---------------------------------------------------------------
         subgridSshEdgeTableRange(1,iEdge) = -maxval(subgridBathymetryValues(1:nSubgridEdge)) + subgrid_thin_layer
         subgridSshEdgeTableRange(2,iEdge) = -minval(subgridBathymetryValues(1:nSubgridEdge))

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
         subgridSshEdgeTableRange(1,iEdge) = -maxval(CellEdgeBathymetryValues(3,1:nSubgridCellEdge)) + subgrid_thin_layer
         subgridSshEdgeTableRange(2,iEdge) = -minval(CellEdgeBathymetryValues(3,1:nSubgridCellEdge))

         ! Evaluate bounds of look-up table range
         !---------------------------------------------------------------
         subgridEdgeBathymetryMean(iEdge) = sum(cellEdgeBathymetryValues(3,1:nSubgridCellEdge))/real(nSubgridCellEdge,RKIND) 
         subgridEdgeBathymetryMin(iEdge)  = maxval(cellEdgeBathymetryValues(3,1:nSubgridCellEdge))

         ! Vertical integration of wet fraction
         !---------------------------------------------------------------
         call ocn_init_vertical_integration( iEdge, subgridSshEdgeTableRange, &
                    nSubgridCellEdge, cellEdgeBathymetryValues(3,:), dsEdge, subgridWetVolumeEdgeTable, subgridWetFractionEdgeTable )
      endif

    end do

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
      subgridSshVertexTableRange(1,iVertex) = -maxval(subgridBathymetryValues(1:nSubgridVertex)) + subgrid_thin_layer
      subgridSshVertexTableRange(2,iVertex) = -minval(subgridBathymetryValues(1:nSubgridVertex))

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
                                             xSubgridCell, ySubgridCell) !{{{

           
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

   !--------------------------------------------------------------------

   end subroutine ocn_init_tri_area!}}}

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
      subgridWetVolumeTable(2,iCV) = subgridWetVolumeTable(1,iCV)*sum(subgridAreas(1:nSubgridCV))
      subgridWetFractionTable(1,iCV) = 0.0_RKIND
      do lev = 2,nSubgridTableLevels
        do tri = 1,nSubgridCV

          pVal = 0.0_RKIND
          if (subgridBathymetryValues(tri) + ssh >= 0.0_RKIND) then
            pVal = 1.0_RKIND
          endif

          subgridWetVolumeTable(lev,iCV) = subgridWetVolumeTable(lev,iCV) + pVal*deltaZ*subgridAreas(tri)
          subgridWetFractionTable(lev,iCV) = subgridWetFractionTable(lev,iCV) + pVal*subgridAreas(tri)
  
        enddo

        if (lev < nSubgridTableLevels) then
          subgridWetVolumeTable(lev+1,iCV) = subgridWetVolumeTable(lev,iCV)
        endif

        subgridWetVolumeTable(lev,iCV) = subgridWetVolumeTable(lev,iCV)/sum(subgridAreas(1:nSubgridCV))
        subgridWetFractionTable(lev,iCV) = subgridWetFractionTable(lev,iCV)/sum(subgridAreas(1:nSubgridCV))

        ssh = ssh + deltaZ
      enddo

   !--------------------------------------------------------------------

   end subroutine ocn_init_vertical_integration!}}}


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


end module ocn_init_subgrid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
