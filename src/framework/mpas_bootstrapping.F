! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_bootstrapping

   use mpas_grid_types
   use mpas_dmpar
   use mpas_block_decomp
   use mpas_block_creator
   use mpas_sort
   use mpas_configure
   use mpas_timekeeping
   use mpas_io_streams
   use mpas_io_units


   integer :: readCellStart, readCellEnd, nReadCells
   integer :: readEdgeStart, readEdgeEnd, nReadEdges
   integer :: readVertexStart, readVertexEnd, nReadVertices


   contains


   !***********************************************************************
   !
   !  routine mpas_bootstrap_framework
   !
   !> \brief   Obtains mesh partition, builds halos, and allocates blocks.
   !> \author  Michael Duda
   !> \date    7 September 2014
   !> \details
   !>  This routine is responsible for reading basic decomposed field dimensions,
   !>  determining a block decomposition of the mesh, building halos for fields,
   !>  and, finally, allocating blocklist for a domain.
   !>
   !>  Dimensions required to be present in the grid file:
   !>      nCells
   !>      nEdges
   !>      nVertices
   !>      vertexDegree
   !>      maxEdges
   !>      ***** and, at present, all other dimensions so we can allocate all fields
   !>
   !>  Fields required to be present in the grid file:
   !>      indexToCellID
   !>      indexToEdgeID
   !>      indexToVertexID
   !>      {x,y,z}Cell
   !>      {x,y,z}Edge
   !>      {x,y,z}Vertex
   !>      nEdgesOnCell
   !>      cellsOnCell
   !>      edgesOnCell
   !>      verticesOnCell
   !>      cellsOnEdge
   !>      cellsOnVertex
   !>
   !>  Attributes required to be present in the grid file:
   !>      on_a_sphere
   !>      sphere_radius
   !>      ***** these are needed by mpas_block_creator_finalize_block_init()
   !>            so they can be set in the mesh pool and queried by, e.g., 
   !>            mpas_initialize_vectors()
   !
   !-----------------------------------------------------------------------
   subroutine mpas_bootstrap_framework(domain, mesh_filename) !{{{

      implicit none
   
      type (domain_type), pointer :: domain
      character(len=*), intent(in) :: mesh_filename

      type (block_type), pointer :: readingBlock
   
#include "dim_dummy_defines_noinput.inc"

      integer :: ierr
      type (MPAS_IO_Handle_type) :: inputHandle

      character (len=StrKIND) :: c_on_a_sphere
      real (kind=RKIND) :: r_sphere_radius
   
      type (field1dInteger), pointer :: indexToCellIDField
      type (field1dInteger), pointer :: indexToEdgeIDField
      type (field1dInteger), pointer :: indexToVertexIDField
      type (field1dInteger), pointer :: nEdgesOnCellField
      type (field2dInteger), pointer :: cellsOnCellField
      type (field2dInteger), pointer :: edgesOnCellField
      type (field2dInteger), pointer :: verticesOnCellField
      type (field2dInteger), pointer :: cellsOnEdgeField
      type (field2dInteger), pointer :: cellsOnVertexField

      type (field1dReal), pointer :: xCellField,   yCellField,   zCellField
      type (field1dReal), pointer :: xEdgeField,   yEdgeField,   zEdgeField
      type (field1dReal), pointer :: xVertexField, yVertexField, zVertexField

      type (field1dInteger), pointer :: nCellsSolveField
      type (field1dInteger), pointer :: nVerticesSolveField
      type (field1dInteger), pointer :: nEdgesSolveField

      type (field1DInteger), pointer :: indexToCellID_Block
      type (field1DInteger), pointer :: nEdgesOnCell_Block
      type (field2DInteger), pointer :: cellsOnCell_Block
      type (field2DInteger), pointer :: verticesOnCell_Block
      type (field2DInteger), pointer :: edgesOnCell_Block

      type (field1DInteger), pointer :: indexToVertexID_Block
      type (field2DInteger), pointer :: cellsOnVertex_Block

      type (field1DInteger), pointer :: indexToEdgeID_Block
      type (field2DInteger), pointer :: cellsOnEdge_Block

      integer, dimension(:), pointer :: local_cell_list
      integer, dimension(:), pointer :: block_id, block_start, block_count
      type (graph) :: partial_global_graph_info

      type (MPAS_Time_type) :: startTime
      character(len=StrKIND) :: timeStamp, restartTimeStamp
      character(len=StrKIND) :: filename

      integer, pointer :: config_num_halos, config_number_of_blocks
      logical, pointer :: config_explicit_proc_decomp
      character (len=StrKIND), pointer :: config_block_decomp_file_prefix, config_proc_decomp_file_prefix
      integer :: nHalos


      call mpas_pool_get_config(domain % configs, 'config_num_halos', config_num_halos)
      call mpas_pool_get_config(domain % configs, 'config_number_of_blocks', config_number_of_blocks)
      call mpas_pool_get_config(domain % configs, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
      call mpas_pool_get_config(domain % configs, 'config_block_decomp_file_prefix', config_block_decomp_file_prefix)
      call mpas_pool_get_config(domain % configs, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)

      nHalos = config_num_halos


      inputHandle = MPAS_io_open(trim(mesh_filename), MPAS_IO_READ, MPAS_IO_PNETCDF, ierr)
      if (ierr /= MPAS_IO_NOERR) then
         call mpas_dmpar_abort(domain % dminfo)
      end if


      !
      ! Read global number of cells/edges/vertices
      !
!MGD for bootstrapping, we really only need a small, fixed set of dimensions
!    HOWEVER: in order to allocate all fields, we need all non-derived dimensions...
#include "read_dimensions.inc"
   
      !
      ! Determine the range of cells/edges/vertices that a processor will initially read
      !   from the input file
      !
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nCells, readCellStart, readCellEnd)   
      nReadCells = readCellEnd - readCellStart + 1
   
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nEdges, readEdgeStart, readEdgeEnd)   
      nReadEdges = readEdgeEnd - readEdgeStart + 1
   
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nVertices, readVertexStart, readVertexEnd)   
      nReadVertices = readVertexEnd - readVertexStart + 1

      allocate(readingBlock)
      readingBlock % domain => domain
      readingBlock % blockID = domain % dminfo % my_proc_id
      readingBlock % localBlockID = 0

      !
      ! Allocate and read fields that we will need in order to ultimately work out
      !   which cells/edges/vertices are owned by each block, and which are ghost
      !

      call mpas_io_setup_cell_block_fields(inputHandle, nreadCells, readCellStart, readingBlock, maxEdges, &
                                           indexTocellIDField, xCellField, yCellField, zCellField, nEdgesOnCellField, &
                                           cellsOnCellField, edgesOnCellField, verticesOnCellField, nHalos)

      call mpas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeIDField, &
                                           xEdgeField, yEdgeField, zEdgeField, cellsOnEdgeField, nHalos)

      call mpas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, &
                                             indexToVertexIDField, xVertexField, yVertexField, zVertexField, cellsOnVertexField, &
                                             nHalos)
      !
      ! Set up a graph derived data type describing the connectivity for the cells 
      !   that were read by this process
      ! A partial description is passed to the block decomp module by each process,
      !   and the block decomp module returns with a list of global cell indices
      !   that belong to the block on this process
      !
      partial_global_graph_info % nVertices = nReadCells
      partial_global_graph_info % nVerticesTotal = nCells
      partial_global_graph_info % maxDegree = maxEdges
      partial_global_graph_info % ghostStart = nVertices+1
      allocate(partial_global_graph_info % vertexID(nReadCells))
      allocate(partial_global_graph_info % nAdjacent(nReadCells))
      allocate(partial_global_graph_info % adjacencyList(maxEdges, nReadCells))
   
      partial_global_graph_info % vertexID(:) = indexToCellIDField % array(:)
      partial_global_graph_info % nAdjacent(:) = nEdgesOnCellField % array(:)
      partial_global_graph_info % adjacencyList(:,:) = cellsOnCellField % array(:,:)
   
      ! TODO: Ensure (by renaming or exchanging) that initial cell range on each proc is contiguous
      !       This situation may occur when reading a restart file with cells/edges/vertices written
      !       in a scrambled order

      !
      ! Determine which cells are owned by this process
      !
      call mpas_block_decomp_cells_for_proc(domain % dminfo, partial_global_graph_info, local_cell_list, block_id, block_start, &
                                            block_count, config_number_of_blocks, config_explicit_proc_decomp, &
                                            config_block_decomp_file_prefix, config_proc_decomp_file_prefix)

      deallocate(partial_global_graph_info % vertexID)
      deallocate(partial_global_graph_info % nAdjacent)
      deallocate(partial_global_graph_info % adjacencyList)

      call mpas_block_creator_setup_blocks_and_0halo_cells(nHalos, domain, indexToCellID_Block, local_cell_list, block_id, &
                                                           block_start, block_count)
      call mpas_block_creator_build_0halo_cell_fields(nHalos, indexToCellIDField, nEdgesOnCellField, cellsOnCellField, &
                                                      verticesOnCellField, edgesOnCellField, indexToCellID_Block, &
                                                      nEdgesOnCell_Block, cellsOnCell_Block, verticesOnCell_Block, &
                                                      edgesOnCell_Block)

      call mpas_block_creator_build_0_and_1halo_edge_fields(nHalos, indexToEdgeIDField, cellsOnEdgeField, indexToCellID_Block, &
                                                            nEdgesOnCell_Block, edgesOnCell_Block, indexToEdgeID_Block, &
                                                            cellsOnEdge_Block, nEdgesSolveField)
      call mpas_block_creator_build_0_and_1halo_edge_fields(nHalos, indexToVertexIDField, cellsOnVertexField, indexToCellID_Block, &
                                                            nEdgesOnCell_Block, verticesOnCell_Block, indexToVertexID_Block, &
                                                            cellsOnVertex_Block, nVerticesSolveField)

      call mpas_block_creator_build_cell_halos(nHalos, indexToCellID_Block, nEdgesOnCell_Block, cellsOnCell_Block, &
                                               verticesOnCell_Block, edgesOnCell_Block, nCellsSolveField)

      call mpas_block_creator_build_edge_halos(nHalos, indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, &
                                               edgesOnCell_Block, indexToEdgeID_Block, cellsOnEdge_Block, nEdgesSolveField)
      call mpas_block_creator_build_edge_halos(nHalos, indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, &
                                               verticesOnCell_Block, indexToVertexID_Block, cellsOnVertex_Block, &
                                               nVerticesSolveField)

      !
      ! Before we can allocate blocks, we need the attributes on_a_sphere and sphere_radius so
      !    they can be propagated as configs to subpools
      !
      call MPAS_io_get_att(inputHandle, 'on_a_sphere', c_on_a_sphere, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
         write(stderrUnit,*) 'Warning: Attribute on_a_sphere not found in '//trim(mesh_filename)
         write(stderrUnit,*) '   Setting on_a_sphere to ''YES'''
         domain % on_a_sphere = .true.
      else
         if (index(c_on_a_sphere, 'YES') /= 0) then
            domain % on_a_sphere = .true.
         else
            domain % on_a_sphere = .false.
         end if
      end if

      call MPAS_io_get_att(inputHandle, 'sphere_radius', r_sphere_radius, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
         write(stderrUnit,*) 'Warning: Attribute sphere_radius not found in '//trim(mesh_filename)
         write(stderrUnit,*) '   Setting sphere_radius to 1.0'
         domain % sphere_radius = 1.0
      else
         domain % sphere_radius = r_sphere_radius
      end if


     ! Allocate blocks, and copy indexTo arrays into blocks
     call mpas_block_creator_finalize_block_init(nHalos, domain % blocklist, &
#include "dim_dummy_args.inc"         
                             , nCellsSolveField, nEdgesSolveField, nVerticesSolveField, indexToCellID_Block, indexToEdgeID_Block, &
                             indexToVertexID_Block)

      call mpas_link_fields(domain)

      call MPAS_io_close(inputHandle, ierr)


      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % sendList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % recvList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % copyList)

      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % sendList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % recvList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % copyList)

      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % sendList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % recvList)
      call mpas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % copyList)

      call mpas_deallocate_field(indexToCellIDField)
      call mpas_deallocate_field(indexToEdgeIDField)
      call mpas_deallocate_field(indexToVertexIDField)
      call mpas_deallocate_field(cellsOnCellField)

      call mpas_deallocate_field(edgesOnCellField)
      call mpas_deallocate_field(verticesOnCellField)
      call mpas_deallocate_field(cellsOnEdgeField)
      call mpas_deallocate_field(cellsOnVertexField)

      call mpas_deallocate_field(indexToCellID_Block)
      call mpas_deallocate_field(nEdgesOnCell_Block)
      call mpas_deallocate_field(cellsOnCell_Block)
      call mpas_deallocate_field(verticesOnCell_Block)
      call mpas_deallocate_field(edgesOnCell_Block)
      call mpas_deallocate_field(indexToVertexID_Block)
      call mpas_deallocate_field(cellsOnVertex_Block)
      call mpas_deallocate_field(indexToEdgeID_Block)
      call mpas_deallocate_field(cellsOnEdge_Block)

      call mpas_deallocate_field(nCellsSolveField)
      call mpas_deallocate_field(nVerticesSolveField)
      call mpas_deallocate_field(nEdgesSolveField)

      deallocate(local_cell_list)
      deallocate(block_id)
      deallocate(block_start)
      deallocate(block_count)
      deallocate(readingBlock)

   end subroutine mpas_bootstrap_framework !}}}


   subroutine mpas_io_setup_cell_block_fields(inputHandle, nReadCells, readCellStart, readingBlock, maxEdges, indexToCellID, &
                                              xCell, yCell, zCell, nEdgesOnCell, cellsOnCell, edgesOnCell, verticesOnCell, nHalos) !{{{

     implicit none

     type (MPAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadCells
     integer, intent(in) :: readCellStart
     integer, intent(in) :: maxEdges
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToCellID
     type (field1dReal), pointer :: xCell
     type (field1dReal), pointer :: yCell
     type (field1dReal), pointer :: zCell
     type (field1dInteger), pointer :: nEdgesOnCell
     type (field2dInteger), pointer :: cellsOnCell
     type (field2dInteger), pointer :: edgesOnCell
     type (field2dInteger), pointer :: verticesOnCell
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices

     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     ! Global cell indices
     allocate(indexToCellID)
     allocate(indexToCellID % ioinfo)
     indexToCellID % ioinfo % fieldName = 'indexToCellID'
     indexToCellID % ioinfo % start(1) = readCellStart
     indexToCellID % ioinfo % count(1) = nReadCells
     allocate(indexToCellID % array(nReadCells))
     allocate(readIndices(nReadCells))
     do i=1,nReadCells
        readIndices(i) = i + readCellStart - 1
     end do
     call MPAS_io_inq_var(inputHandle, 'indexToCellID', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'indexToCellID', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'indexToCellID', indexToCellID % array, ierr)
     indexToCellID % dimSizes(1) = nReadCells
     indexToCellID % block => readingBlock
     call mpas_dmpar_init_multihalo_exchange_list(indexToCellID % sendList, nHalos)
     call mpas_dmpar_init_multihalo_exchange_list(indexToCellID % recvList, nHalos)
     call mpas_dmpar_init_multihalo_exchange_list(indexToCellID % copyList, nHalos)
     nullify(indexToCellID % next)
   

     ! Number of cell/edges/vertices adjacent to each cell
     allocate(nEdgesOnCell)
     allocate(nEdgesOnCell % ioinfo)
     nEdgesOnCell % ioinfo % fieldName = 'nEdgesOnCell'
     nEdgesOnCell % ioinfo % start(1) = readCellStart
     nEdgesOnCell % ioinfo % count(1) = nReadCells
     allocate(nEdgesOnCell % array(nReadCells))
     call MPAS_io_inq_var(inputHandle, 'nEdgesOnCell', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'nEdgesOnCell', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'nEdgesOnCell', nEdgesOnCell % array, ierr)
     nEdgesOnCell % dimSizes(1) = nReadCells
     nEdgesOnCell % block => readingBlock
     nEdgesOnCell % sendList => indexToCellID % sendList
     nEdgesOnCell % recvList => indexToCellID % recvList
     nEdgesOnCell % copyList => indexToCellID % copyList
     nullify(nEdgesOnCell % next)
   
     ! Global indices of cells adjacent to each cell
     allocate(cellsOnCell)
     allocate(cellsOnCell % ioinfo)
     cellsOnCell % ioinfo % fieldName = 'cellsOnCell'
     cellsOnCell % ioinfo % start(1) = 1
     cellsOnCell % ioinfo % start(2) = readCellStart
     cellsOnCell % ioinfo % count(1) = maxEdges
     cellsOnCell % ioinfo % count(2) = nReadCells
     allocate(cellsOnCell % array(maxEdges,nReadCells))
     call MPAS_io_inq_var(inputHandle, 'cellsOnCell', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'cellsOnCell', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'cellsOnCell', cellsOnCell % array, ierr)
     cellsOnCell % dimSizes(1) = maxEdges
     cellsOnCell % dimSizes(2) = nReadCells
     cellsOnCell % block => readingBlock
     cellsOnCell % sendList => indexToCellID % sendList
     cellsOnCell % recvList => indexToCellID % recvList
     cellsOnCell % copyList => indexToCellID % copyList
     nullify(cellsOnCell % next)
   
     ! Global indices of edges adjacent to each cell
     allocate(edgesOnCell)
     allocate(edgesOnCell % ioinfo)
     edgesOnCell % ioinfo % fieldName = 'edgesOnCell'
     edgesOnCell % ioinfo % start(1) = 1
     edgesOnCell % ioinfo % start(2) = readCellStart
     edgesOnCell % ioinfo % count(1) = maxEdges
     edgesOnCell % ioinfo % count(2) = nReadCells
     allocate(edgesOnCell % array(maxEdges,nReadCells))
     call MPAS_io_inq_var(inputHandle, 'edgesOnCell', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'edgesOnCell', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'edgesOnCell', edgesOnCell % array, ierr)
     edgesOnCell % dimSizes(1) = maxEdges
     edgesOnCell % dimSizes(2) = nReadCells
     edgesOnCell % block => readingBlock
     edgesOnCell % sendList => indexToCellID % sendList
     edgesOnCell % recvList => indexToCellID % recvList
     edgesOnCell % copyList => indexToCellID % copyList
     nullify(edgesOnCell % next)
   
     ! Global indices of vertices adjacent to each cell
     allocate(verticesOnCell)
     allocate(verticesOnCell % ioinfo)
     verticesOnCell % ioinfo % fieldName = 'verticesOnCell'
     verticesOnCell % ioinfo % start(1) = 1
     verticesOnCell % ioinfo % start(2) = readCellStart
     verticesOnCell % ioinfo % count(1) = maxEdges
     verticesOnCell % ioinfo % count(2) = nReadCells
     allocate(verticesOnCell % array(maxEdges,nReadCells))
     call MPAS_io_inq_var(inputHandle, 'verticesOnCell', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'verticesOnCell', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'verticesOnCell', verticesOnCell % array, ierr)
     verticesOnCell % dimSizes(1) = maxEdges
     verticesOnCell % dimSizes(2) = nReadCells
     verticesOnCell % block => readingBlock
     verticesOnCell % sendList => indexToCellID % sendList
     verticesOnCell % recvList => indexToCellID % recvList
     verticesOnCell % copyList => indexToCellID % copyList
     nullify(verticesOnCell % next)

     deallocate(readIndices)
   
   end subroutine mpas_io_setup_cell_block_fields !}}}


   subroutine mpas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeID, &
                                              xEdge, yEdge, zEdge, cellsOnEdge, nHalos) !{{{

     implicit none

     type (MPAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadEdges
     integer, intent(in) :: readEdgeStart
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToEdgeID
     type (field1dReal), pointer :: xEdge
     type (field1dReal), pointer :: yEdge
     type (field1dReal), pointer :: zEdge
     type (field2dInteger), pointer :: cellsOnEdge
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices

     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     allocate(readIndices(nReadEdges))

     ! Global edge indices
     allocate(indexToEdgeID)
     allocate(indexToEdgeID % ioinfo)
     indexToEdgeID % ioinfo % fieldName = 'indexToEdgeID'
     indexToEdgeID % ioinfo % start(1) = readEdgeStart
     indexToEdgeID % ioinfo % count(1) = nReadEdges
     allocate(indexToEdgeID % array(nReadEdges))
     allocate(indexToEdgeID % array(nReadEdges))
     do i=1,nReadEdges
        readIndices(i) = i + readEdgeStart - 1
     end do
     call MPAS_io_inq_var(inputHandle, 'indexToEdgeID', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'indexToEdgeID', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'indexToEdgeID', indexToEdgeID % array, ierr)
     indexToEdgeID % dimSizes(1) = nREadEdges
     indexToEdgeID % block => readingBlock
     call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeID % sendList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeID % recvList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeID % copyList, nHalos+1)
     nullify(indexToEdgeID % next)
   
   
     ! Global indices of cells adjacent to each edge
     !    used for determining which edges are owned by a block, where 
     !    iEdge is owned iff cellsOnEdge(1,iEdge) is an owned cell
     allocate(cellsOnEdge)
     allocate(cellsOnEdge % ioinfo)
     cellsOnEdge % ioinfo % fieldName = 'cellsOnEdge'
     cellsOnEdge % ioinfo % start(1) = 1
     cellsOnEdge % ioinfo % start(2) = readEdgeStart
     cellsOnEdge % ioinfo % count(1) = 2
     cellsOnEdge % ioinfo % count(2) = nReadEdges
     allocate(cellsOnEdge % array(2,nReadEdges))
     call MPAS_io_inq_var(inputHandle, 'cellsOnEdge', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'cellsOnEdge', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'cellsOnEdge', cellsOnEdge % array, ierr)
     cellsOnEdge % dimSizes(1) = 2
     cellsOnEdge % dimSizes(2) = nReadEdges
     cellsOnEdge % block => readingBlock
     cellsOnEdge % sendList => indexToEdgeID % sendList
     cellsOnEdge % recvList => indexToEdgeID % recvList
     cellsOnEdge % copyList => indexToEdgeID % copyList
     nullify(cellsOnEdge % next)

     deallocate(readIndices)
   
   end subroutine mpas_io_setup_edge_block_fields !}}}


   subroutine mpas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, &
                                                indexToVertexID, xVertex, yVertex, zVertex, cellsOnVertex, nHalos) !{{{

     implicit none

     type (MPAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadVertices
     integer, intent(in) :: readVertexStart
     integer, intent(in) :: vertexDegree
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToVertexID
     type (field1dReal), pointer :: xVertex
     type (field1dReal), pointer :: yVertex
     type (field1dReal), pointer :: zVertex
     type (field2dInteger), pointer :: cellsOnVertex
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices

     ! Global vertex indices
     allocate(indexToVertexID)
     allocate(indexToVertexID % ioinfo)
     indexToVertexID % ioinfo % fieldName = 'indexToVertexID'
     indexToVertexID % ioinfo % start(1) = readVertexStart
     indexToVertexID % ioinfo % count(1) = nReadVertices
     allocate(indexToVertexID % array(nReadVertices))
     allocate(readIndices(nReadVertices))
     do i=1,nReadVertices
        readIndices(i) = i + readVertexStart - 1
     end do
     call MPAS_io_inq_var(inputHandle, 'indexToVertexID', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'indexToVertexID', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'indexToVertexID', indexToVertexID % array, ierr)
     indexToVertexID % dimSizes(1) = nReadVertices
     indexToVertexID % block => readingBlock
     call mpas_dmpar_init_multihalo_exchange_list(indexToVertexID % sendList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToVertexID % recvList, nHalos+1)
     call mpas_dmpar_init_multihalo_exchange_list(indexToVertexID % copyList, nHalos+1)
     nullify(indexToVertexID % next)
   
   
     ! Global indices of cells adjacent to each vertex
     !    used for determining which vertices are owned by a block, where 
     !    iVtx is owned iff cellsOnVertex(1,iVtx) is an owned cell
     allocate(cellsOnVertex)
     allocate(cellsOnVertex % ioinfo)
     cellsOnVertex % ioinfo % fieldName = 'cellsOnVertex'
     cellsOnVertex % ioinfo % start(1) = 1
     cellsOnVertex % ioinfo % start(2) = readVertexStart
     cellsOnVertex % ioinfo % count(1) = vertexDegree
     cellsOnVertex % ioinfo % count(2) = nReadVertices
     allocate(cellsOnVertex % array(vertexDegree,nReadVertices))
     call MPAS_io_inq_var(inputHandle, 'cellsOnVertex', ierr=ierr)
     call MPAS_io_set_var_indices(inputHandle, 'cellsOnVertex', readIndices, ierr=ierr)
     call mpas_io_get_var(inputHandle, 'cellsOnVertex', cellsOnVertex % array, ierr)
     cellsOnVertex % dimSizes(1) = vertexDegree
     cellsOnVertex % dimSizes(2) = nReadVertices
     cellsOnVertex % block => readingBlock
     cellsOnVertex % sendList => indexToVertexID % sendList
     cellsOnVertex % recvList => indexToVertexID % recvList
     cellsOnVertex % copyList => indexToVertexID % copyList
     nullify(cellsOnVertex % next)

     deallocate(readIndices)

   end subroutine mpas_io_setup_vertex_block_fields !}}}

#include "setup_immutable_streams.inc"

end module mpas_bootstrapping
