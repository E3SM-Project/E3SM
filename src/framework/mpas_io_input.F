! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_io_input

   use mpas_grid_types
   use mpas_dmpar
   use mpas_block_decomp
   use mpas_block_creator
   use mpas_sort
   use mpas_configure
   use mpas_timekeeping
   use mpas_io_streams
   use mpas_io_units

   integer, parameter :: STREAM_INPUT=1, STREAM_SFC=2, STREAM_RESTART=3
   integer, parameter :: idLength = 10

   type io_input_object
      character (len=StrKIND) :: filename
      integer :: rd_ncid
      integer :: stream

      integer :: time

      type (MPAS_Stream_type) :: io_stream

   end type io_input_object

   integer :: readCellStart, readCellEnd, nReadCells
   integer :: readEdgeStart, readEdgeEnd, nReadEdges
   integer :: readVertexStart, readVertexEnd, nReadVertices

   contains

   subroutine mpas_input_state_for_domain(domain)!{{{

      implicit none
   
      type (domain_type), pointer :: domain

      type (block_type), pointer :: block_ptr
      type (block_type), pointer :: readingBlock
   
      integer :: i, j, k
      type (io_input_object) :: input_obj
#include "dim_decls.inc"

      character (len=StrKIND) :: c_on_a_sphere
      real (kind=RKIND) :: r_sphere_radius

      integer :: ierr
      integer, dimension(:), pointer :: readIndices
      type (MPAS_IO_Handle_type) :: inputHandle
   
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

      type (field1DChar) :: xtime

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

      type (field1DReal), pointer :: xCell, yCell, zCell
      type (field1DReal), pointer :: xEdge, yEdge, zEdge
      type (field1DReal), pointer :: xVertex, yVertex, zVertex
   
      integer, dimension(:), pointer :: local_cell_list
      integer, dimension(:), pointer :: block_id, block_start, block_count
      type (graph) :: partial_global_graph_info

      type (MPAS_Time_type) :: startTime
      character(len=StrKIND) :: timeStamp, restartTimeStamp
      character(len=StrKIND) :: filename

      integer :: nHalos

      nHalos = config_num_halos

      if (config_do_restart) then
        ! this get followed by set is to ensure that the time is in standard format
        if(trim(config_start_time) == 'file') then
          open(22,file=trim(config_restart_timestamp_name),form='formatted',status='old')
          read(22,*) restartTimeStamp
          close(22)

        else
          restartTimeStamp = config_start_time
        end if

        write(stderrUnit,*) 'RestartTimeStamp ', trim(restartTimeStamp)
        call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp)
        call mpas_get_time(curr_time=startTime, dateTimeString=timeStamp)
        call mpas_insert_string_suffix(trim(config_restart_name), timeStamp, filename)

        input_obj % filename = trim(filename)
        input_obj % stream = STREAM_RESTART
      else
        input_obj % filename = trim(config_input_name)
        input_obj % stream = STREAM_INPUT
      end if
      inputHandle = MPAS_io_open(trim(input_obj % filename), MPAS_IO_READ, MPAS_IO_PNETCDF, ierr)
      if (ierr /= MPAS_IO_NOERR) then
        write(stderrUnit,*) ' '
        if (input_obj % stream == STREAM_RESTART) then
          write(stderrUnit,*) 'Error opening restart file ''', trim(input_obj % filename), ''''
        else if (input_obj % stream == STREAM_INPUT) then
          write(stderrUnit,*) 'Error opening input file ''', trim(input_obj % filename), ''''
        else if (input_obj % stream == STREAM_SFC) then
          write(stderrUnit,*) 'Error opening sfc file ''', trim(input_obj % filename), ''''
        end if
        write(stderrUnit,*) ' '
        call mpas_dmpar_abort(domain % dminfo)
      end if

      !
      ! Read global number of cells/edges/vertices
      !
#include "read_dims.inc"
   
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

      call mpas_io_setup_cell_block_fields(inputHandle, nreadCells, readCellStart, readingBlock, maxEdges, indexTocellIDField, xCellField, &
                                           yCellField, zCellField, nEdgesOnCellField, cellsOnCellField, edgesOnCellField, verticesOnCellField)

      call mpas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeIDField, xEdgeField, yEdgeField, zEdgeField, cellsOnEdgeField)

      call mpas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, indexToVertexIDField, &
                                             xVertexField, yVertexField, zVertexField, cellsOnVertexField)
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
   

      ! Determine which cells are owned by this process
      call mpas_block_decomp_cells_for_proc(domain % dminfo, partial_global_graph_info, local_cell_list, block_id, block_start, block_count)

      deallocate(partial_global_graph_info % vertexID)
      deallocate(partial_global_graph_info % nAdjacent)
      deallocate(partial_global_graph_info % adjacencyList)

      call mpas_block_creator_setup_blocks_and_0halo_cells(domain, indexToCellID_Block, local_cell_list, block_id, block_start, block_count)
      call mpas_block_creator_build_0halo_cell_fields(indexToCellIDField, nEdgesOnCellField, cellsOnCellField, verticesOnCellField, edgesOnCellField, indexToCellID_Block, nEdgesOnCell_Block, cellsOnCell_Block, verticesOnCell_Block, edgesOnCell_Block)

      call mpas_block_creator_build_0_and_1halo_edge_fields(indexToEdgeIDField, cellsOnEdgeField, indexToCellID_Block, nEdgesOnCell_Block, edgesOnCell_Block, indexToEdgeID_Block, cellsOnEdge_Block, nEdgesSolveField)
      call mpas_block_creator_build_0_and_1halo_edge_fields(indexToVertexIDField, cellsOnVertexField, indexToCellID_Block, nEdgesOnCell_Block, verticesOnCell_Block, indexToVertexID_Block, cellsOnVertex_Block, nVerticesSolveField)

      call mpas_block_creator_build_cell_halos(indexToCellID_Block, nEdgesOnCell_Block, cellsOnCell_Block, verticesOnCell_Block, edgesOnCell_Block, nCellsSolveField)

      call mpas_block_creator_build_edge_halos(indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, edgesOnCell_Block, indexToEdgeID_Block, cellsOnEdge_Block, nEdgesSolveField)
      call mpas_block_creator_build_edge_halos(indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, verticesOnCell_Block, indexToVertexID_Block, cellsOnVertex_Block, nVerticesSolveField)


     ! Allocate blocks, and copy indexTo arrays into blocks
     call mpas_block_creator_finalize_block_init(domain % blocklist, &
#include "dim_dummy_args.inc"         
                             , nCellsSolveField, nEdgesSolveField, nVerticesSolveField, indexToCellID_Block, indexToEdgeID_Block, indexToVertexID_Block)


      call mpas_io_input_init(input_obj, domain % blocklist, domain % dminfo)

      call MPAS_readStreamAtt(input_obj % io_stream, 'sphere_radius', r_sphere_radius, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
        write(stderrUnit,*) 'Warning: Attribute sphere_radius not found in '//trim(input_obj % filename)
        write(stderrUnit,*) '   Setting sphere_radius to 1.0'
        domain % blocklist % mesh % sphere_radius = 1.0
      else
        domain % blocklist % mesh % sphere_radius = r_sphere_radius
      end if

      call MPAS_readStreamAtt(input_obj % io_stream, 'on_a_sphere', c_on_a_sphere, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
        write(stderrUnit,*) 'Warning: Attribute on_a_sphere not found in '//trim(input_obj % filename)
        write(stderrUnit,*) '   Setting on_a_sphere to ''YES'''
        domain % blocklist % mesh % on_a_sphere = .true.
      else
        if (index(c_on_a_sphere, 'YES') /= 0) then
          domain % blocklist % mesh % on_a_sphere = .true.
        else
          domain % blocklist % mesh % on_a_sphere = .false.
        end if
      end if

#ifndef MPAS_CESM
      call MPAS_readStreamAtt(input_obj % io_stream, 'history', domain % history, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
        write(stderrUnit,*) 'Warning: Attribute History not found in '//trim(input_obj % filename)
        write(stderrUnit,*) '   Setting History to '''''
        domain % history = ""
      else
        ! Remove C String NULL characters, replace C String newlines with semicolons
        do i = 1, len(domain % history)
          if(iachar(domain % history(i:i)) == 0) then
            domain % history(i:i) = " "
          else if(iachar(domain % history(i:i)) == 10) then
            domain % history(i:i) = ";"
          end if
        end do
      end if
#else
      domain % history = "cesm_run"
#endif

      call MPAS_readStreamAtt(input_obj % io_stream, 'parent_id', domain % parent_id, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
         domain % parent_id = 'Unknown'
      end if

      call MPAS_readStreamAtt(input_obj % io_stream, 'mesh_id', domain % mesh_id, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
         domain % mesh_id = 'Unknown'
      else
         domain % parent_id = domain % mesh_id // "\n" // domain % parent_id
      end if

      call gen_random(idLength, domain % mesh_id)

      call MPAS_readStreamAtt(input_obj % io_stream, 'Conventions', domain % Conventions, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
        write(stderrUnit,*) 'Warning: Attribute Conventions not found in '//trim(input_obj % filename)
        write(stderrUnit,*) '   Setting Conventions to ''MPAS'''
        domain % Conventions = "MPAS"
      end if

      call MPAS_readStreamAtt(input_obj % io_stream, 'source', domain % source, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
        write(stderrUnit,*) 'Warning: Attribute source not found in '//trim(input_obj % filename)
        write(stderrUnit,*) '   Setting source to ''MPAS'''
        domain % source = "MPAS"
      end if

      call MPAS_readStreamAtt(input_obj % io_stream, 'mesh_spec', domain % mesh_spec, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
        write(stderrUnit,*) 'Warning: Attribute mesh_spec not found in '//trim(input_obj % filename)
        write(stderrUnit,*) '   Setting mesh_spec to ''Unknown'''
        domain % mesh_spec = "Unknown"
      end if

      block_ptr => domain % blocklist % next
      do while (associated(block_ptr))
        block_ptr % mesh % sphere_radius = domain % blocklist % mesh % sphere_radius
        block_ptr % mesh % on_a_sphere = domain % blocklist % mesh % on_a_sphere

        ! Link the sendList and recvList pointers in each field type to the appropriate lists 
        !   in parinfo, e.g., cellsToSend and cellsToRecv; in future, it can also be extended to 
        !   link blocks of fields to eachother
        call mpas_create_field_links(block_ptr)

        block_ptr => block_ptr % next
      end do

      if (.not. config_do_restart) then
        input_obj % time = 1
      else
        !
        ! If doing a restart, we need to decide which time slice to read from the 
        !   restart file
        !
        input_obj % time = MPAS_seekStream(input_obj % io_stream, restartTimeStamp, MPAS_STREAM_EXACT_TIME, timeStamp, ierr)
        if (ierr == MPAS_IO_ERR) then
          write(stderrUnit,*) 'Error: restart file '//trim(filename)//' did not contain time '//trim(restartTimeStamp)
          call mpas_dmpar_abort(domain % dminfo)
        end if

!       input_obj % time = MPAS_seekStream(input_obj % io_stream, config_start_time, MPAS_STREAM_EXACT_TIME, timeStamp, ierr)
!       if (ierr == MPAS_IO_ERR) then
!         write(stderrUnit,*) 'Error: restart file '//trim(filename)//' did not contain time '//trim(config_start_time)
!         call mpas_dmpar_abort(domain % dminfo)
!       end if
!write(stderrUnit,*) 'MGD DEBUGGING time = ', input_obj % time
        write(stderrUnit,*) 'Restarting model from time ', trim(timeStamp)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Do the actual work of reading all fields in from the input or restart file
      ! For each field:
      !   1) Each process reads a contiguous range of cell/edge/vertex indices, which
      !      may not correspond with the cells/edges/vertices that are owned by the
      !      process
      !   2) All processes then send the global indices that were read to the 
      !      processes that own those indices based on 
      !      {send,recv}{Cell,Edge,Vertex}List
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      call mpas_read_and_distribute_fields(input_obj)

      call mpas_io_input_finalize(input_obj, domain % dminfo)

      call MPAS_io_close(inputHandle, ierr)

      !
      ! Exchange halos for all of the fields that were read from the input file
      !
      call mpas_exch_input_field_halos(domain, input_obj)

      call mpas_block_creator_reindex_block_fields(domain % blocklist)

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

   end subroutine mpas_input_state_for_domain!}}}

   !CR:TODO: an identical subroutine is found in module_io_output - merge
   subroutine mpas_insert_string_suffix(stream, suffix, filename)!{{{

      implicit none

      character (len=*), intent(in) :: stream
      character (len=*), intent(in) :: suffix
      character (len=*), intent(out) :: filename
      integer :: length, i

      filename = trim(stream) // '.' // trim(suffix)

      length = len_trim(stream)
      do i=length-1,1,-1
         if(stream(i:i) == '.') then
            filename = trim(stream(:i)) // trim(suffix) // trim(stream(i:))
            exit
         end if
      end do

      do i=1,len_trim(filename)
         if (filename(i:i) == ':') filename(i:i) = '.'
      end do

   end subroutine mpas_insert_string_suffix!}}}

   subroutine mpas_read_and_distribute_fields(input_obj)!{{{
      
      implicit none

      type (io_input_object), intent(inout) :: input_obj

      integer :: ierr


      call MPAS_readStream(input_obj % io_stream, input_obj % time, ierr)


   end subroutine mpas_read_and_distribute_fields!}}}

   subroutine mpas_io_input_init(input_obj, blocklist, dminfo)!{{{
 
      implicit none

      type (io_input_object), intent(inout) :: input_obj
      type (block_type), intent(in) :: blocklist
      type (dm_info), intent(in) :: dminfo
 
      integer :: nferr
 
      call MPAS_createStream(input_obj % io_stream, trim(input_obj % filename), MPAS_IO_PNETCDF, MPAS_IO_READ, 1, nferr)
      if (nferr /= MPAS_STREAM_NOERR) then
         write(stderrUnit,*) ' '
         if (input_obj % stream == STREAM_RESTART) then
            write(stderrUnit,*) 'Error opening restart file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_INPUT) then
            write(stderrUnit,*) 'Error opening input file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_SFC) then
            write(stderrUnit,*) 'Error opening sfc file ''', trim(input_obj % filename), ''''
         end if
         write(stderrUnit,*) ' '
         call mpas_dmpar_abort(dminfo)
      end if

#include "add_input_fields.inc"

   end subroutine mpas_io_input_init!}}}

   subroutine mpas_exch_input_field_halos(domain, input_obj)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain
      type (io_input_object), intent(inout) :: input_obj

#include "exchange_input_field_halos.inc"

#include "non_decomp_copy_input_fields.inc"

   end subroutine mpas_exch_input_field_halos!}}}

   subroutine mpas_io_input_finalize(input_obj, dminfo)!{{{
 
      implicit none
 
      type (io_input_object), intent(inout) :: input_obj
      type (dm_info), intent(in) :: dminfo

      integer :: nferr
 
      call MPAS_closeStream(input_obj % io_stream, nferr)
 
   end subroutine mpas_io_input_finalize!}}}

   subroutine mpas_io_setup_cell_block_fields(inputHandle, nReadCells, readCellStart, readingBlock, maxEdges, indexToCellID, xCell, yCell, zCell, nEdgesOnCell, cellsOnCell, edgesOnCell, verticesOnCell)!{{{
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

     integer :: i, nHalos
     integer, dimension(:), pointer :: readIndices

     nHalos = config_num_halos
  
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
   
   end subroutine mpas_io_setup_cell_block_fields!}}}

   subroutine mpas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeID, xEdge, yEdge, zEdge, cellsOnEdge)!{{{
     type (MPAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadEdges
     integer, intent(in) :: readEdgeStart
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToEdgeID
     type (field1dReal), pointer :: xEdge
     type (field1dReal), pointer :: yEdge
     type (field1dReal), pointer :: zEdge
     type (field2dInteger), pointer :: cellsOnEdge

     integer :: i, nHalos
     integer, dimension(:), pointer :: readIndices

     nHalos = config_num_halos
  
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
   
   end subroutine mpas_io_setup_edge_block_fields!}}}

   subroutine mpas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, indexToVertexID, xVertex, yVertex, zVertex, cellsOnVertex)!{{{
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

     integer :: i, nHalos
     integer, dimension(:), pointer :: readIndices

     nHalos = config_num_halos
  
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

   end subroutine mpas_io_setup_vertex_block_fields!}}}

 
end module mpas_io_input
