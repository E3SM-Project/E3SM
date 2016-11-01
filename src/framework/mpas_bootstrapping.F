! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_bootstrapping

   use mpas_derived_types
   use mpas_field_routines
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_block_decomp
   use mpas_block_creator
   use mpas_sort
   use mpas_timekeeping
   use mpas_io_streams
   use mpas_io_units
   use mpas_abort, only : mpas_dmpar_global_abort


   integer :: readCellStart, readCellEnd, nReadCells
   integer :: readEdgeStart, readEdgeEnd, nReadEdges
   integer :: readVertexStart, readVertexEnd, nReadVertices


   contains


   !***********************************************************************
   !
   !  routine mpas_bootstrap_framework_phase1
   !
   !> \brief   Obtains mesh partition, builds halos, and allocates blocks.
   !> \author  Michael Duda, Doug Jacobsen
   !> \date    03/04/2015
   !> \details
   !>  This routine is responsible for reading basic decomposed field dimensions,
   !>  determining a block decomposition of the mesh, building halos for fields,
   !>  and, finally, allocating blocklist for a domain.
   !>
   !>  Phase1 does not currently allocate fields, or define any other
   !>  dimensions aside from the ones listed below.
   !>
   !>  Dimensions required to be present in the grid file:
   !>      nCells
   !>      nEdges
   !>      nVertices
   !>      vertexDegree
   !>      maxEdges
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
   !>      parent_id
   !>      mesh_spec
   !>      is_periodic
   !>      x_period
   !>      y_period
   !>      ***** these are needed by mpas_block_creator_finalize_block_init()
   !>            so they can be set in the mesh pool and queried by, e.g.,
   !>            mpas_initialize_vectors()
   !
   !-----------------------------------------------------------------------
   subroutine mpas_bootstrap_framework_phase1(domain, mesh_filename, mesh_iotype) !{{{

      implicit none

      type (domain_type), pointer :: domain
      character(len=*), intent(in) :: mesh_filename
      integer, intent(in) :: mesh_iotype

      type (block_type), pointer :: readingBlock

      integer :: nCells, nEdges, maxEdges, nVertices, vertexDegree

      integer :: ierr
      type (MPAS_IO_Handle_type) :: inputHandle

      character (len=StrKIND) :: c_on_a_sphere, c_is_periodic, c_parent_id, c_file_id, c_mesh_spec
      real (kind=RKIND) :: r_sphere_radius, r_x_period, r_y_period

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


      inputHandle = MPAS_io_open(trim(mesh_filename), MPAS_IO_READ, mesh_iotype, domain % ioContext, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
         write(stderrUnit,*) ' '
         write(stderrUnit,*) '************************************************************************************'
         write(stderrUnit,*) 'Error: Could not open input file '''//trim(mesh_filename)//''' to read mesh fields'
         write(stderrUnit,*) '************************************************************************************'
         call mpas_dmpar_abort(domain % dminfo)
      else
         write(stderrUnit,*) 'Bootstrapping framework with mesh fields from input file '''//trim(mesh_filename)//''''
      end if


      !
      ! Read global number of cells/edges/vertices
      !
      call mpas_io_inq_dim(inputHandle, 'nCells', nCells, ierr)
      call mpas_io_inq_dim(inputHandle, 'nEdges', nEdges, ierr)
      call mpas_io_inq_dim(inputHandle, 'maxEdges', maxEdges, ierr)
      call mpas_io_inq_dim(inputHandle, 'nVertices', nVertices, ierr)
      call mpas_io_inq_dim(inputHandle, 'vertexDegree', vertexDegree, ierr)

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

      !
      ! We need to read in the file_id, parent_id, and mesh spec attributes as well, to
      ! set them appropriately in the output files.
      !
      call MPAS_io_get_att(inputHandle, 'file_id', c_file_id, ierr=ierr)
      if ( ierr /= MPAS_IO_NOERR) then
         write(stderrUnit, *) 'Warning: Attribute file_id not found in '//trim(mesh_filename)
         write(stderrUnit, *) '  Using '''' for the previous file_id.'
         c_file_id = ''
      end if

      call MPAS_io_get_att(inputHandle, 'parent_id', c_parent_id, ierr=ierr)
      if ( ierr /= MPAS_IO_NOERR) then
         write(stderrUnit, *) 'Warning: Attribute parent_id not found in '//trim(mesh_filename)
         write(stderrUnit, *) '  Setting parent_id to '''''
         c_parent_id = ''
      end if

      if ( trim(c_file_id) == '' ) then
         domain % parent_id = c_parent_id
      else
         domain % parent_id = trim(c_file_id) // NEW_LINE('A') // trim(c_parent_id)
      end if

      call MPAS_io_get_att(inputHandle, 'mesh_spec', c_mesh_spec, ierr=ierr)
      if ( ierr /= MPAS_IO_NOERR) then
         write(stderrUnit, *) 'Warning: Attribute mesh_spec not found in '//trim(mesh_filename)
         write(stderrUnit, *) '  Setting mesh_spec to ''0.0'''
         domain % mesh_spec = '0.0'
      else
        domain % mesh_spec = c_mesh_spec
      end if

      !
      ! Get attributes based on planar mesh
      !
      if ( .not. domain % on_a_sphere ) then
         call MPAS_io_get_att(inputHandle, 'is_periodic', c_is_periodic, ierr=ierr)
         if ( ierr /= MPAS_IO_NOERR ) then
            write(stderrUnit, *) 'Warning: is_periodic attribute not found in '//trim(mesh_filename)
            write(stderrUnit, *) '   Setting is_periodic to .false.'
            domain % is_periodic = .false.
         else
            if (index(c_is_periodic, 'YES') /= 0) then
               domain % is_periodic = .true.
            else
               domain % is_periodic = .false.
            end if
         end if
      else
        domain % is_periodic = .false.
      end if

      if ( domain % is_periodic ) then
         call MPAS_io_get_att(inputHandle, 'x_period', r_x_period, ierr=ierr)
         if ( ierr /= MPAS_IO_NOERR ) then
            call mpas_dmpar_global_abort('ERROR: Required attribute ''x_period'' missing from periodic input mesh. Exiting...')
         else
            domain % x_period = r_x_period
         end if

         call MPAS_io_get_att(inputHandle, 'y_period', r_y_period, ierr=ierr)
         if ( ierr /= MPAS_IO_NOERR ) then
            call mpas_dmpar_global_abort('ERROR: Required attribute ''y_period'' missing from periodic input mesh. Exiting...')
         else
            domain % y_period = r_y_period
         end if
      else
         domain % x_period = 0
         domain % y_period = 0
      end if


     ! Allocate blocks, and copy indexTo arrays into blocks
     call mpas_block_creator_finalize_block_phase1(nHalos, domain % blocklist, &
             nCellsSolveField, nEdgesSolveField, nVerticesSolveField, vertexDegree, maxEdges, &
             indexToCellID_Block, indexToEdgeID_Block, indexToVertexID_Block)

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

   end subroutine mpas_bootstrap_framework_phase1 !}}}


   !***********************************************************************
   !
   !  routine mpas_bootstrap_framework_phase2
   !
   !> \brief   Defines dimensions and allocates fields
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    03/04/2015
   !> \details
   !>  This routine assumes blocks have already been setup correctly (some
   !>  fields have the correct indices, and exchange lists are created).
   !>  It will finish the setup of blocks by defining all remaining dimensions,
   !>  and allocating all fields and structs.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_bootstrap_framework_phase2(domain) !{{{

      use mpas_stream_manager
      use mpas_stream_list

      implicit none

      type (domain_type), pointer :: domain

      type (mpas_pool_type), pointer :: readableDimensions
      type (mpas_pool_type), pointer :: streamDimensions

      character (len=StrKIND) :: streamName, streamFilenameTemplate, streamFilenameInterval, streamReferenceTime, streamFilename
      character (len=StrKIND) :: fieldName
      integer :: streamDirection
      logical :: streamActive
      logical :: fieldActive
      integer :: ioType

      type (mpas_pool_iterator_type) :: poolItr

      type (MPAS_IO_Handle_type) :: inputHandle

      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, pointer :: dimValue
      integer :: tempDim
      integer :: i, err_level, err_local


      call mpas_pool_create_pool(readableDimensions)

      err_level = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)

      write(stderrUnit, *) ' '
      write(stderrUnit, *) ' '

      ! Reading dimensions from streams
      write(stderrUnit,'(a)') 'Reading dimensions from input streams ...'
      call mpas_stream_mgr_begin_iteration(domain % streamManager)
      do while ( mpas_stream_mgr_get_next_stream(domain % streamManager, streamID = streamName, directionProperty = streamDirection, &
                                                 activeProperty = streamActive) )

         if ( streamActive .and. ( streamDirection == MPAS_STREAM_INPUT .or. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) ) then

            call mpas_stream_mgr_begin_iteration(domain % streamManager, streamID=streamName)

            ! Build stream dimension pool from the list of fields
            call mpas_pool_create_pool(streamDimensions)

            do while ( mpas_stream_mgr_get_next_field(domain % streamManager, streamName, fieldName, isActive=fieldActive) )

               if (fieldActive) then
                  call get_dimlist_for_field(domain % blocklist % allFields, fieldName, dimNames)

                  do i=1,size(dimNames)
                     call mpas_pool_get_dimension(streamDimensions, dimNames(i), dimValue)
                     if ( .not. associated(dimValue) ) then
                        call mpas_pool_add_dimension(streamDimensions, dimNames(i), MPAS_MISSING_DIM)
                     end if
                  end do
                  deallocate(dimNames)
               end if

            end do

            ! Determine stream filename
            call mpas_get_stream_filename(domain % streamManager, streamID = streamName, filename = streamFilename, ierr = err_local)

            ! Determine stream io_type
            call MPAS_stream_mgr_get_property(domain % streamManager, streamName, &
                                              MPAS_STREAM_PROPERTY_IOTYPE, ioType, ierr = err_local)

            ! Try to open file
            inputHandle = MPAS_io_open(trim(streamFilename), MPAS_IO_READ, ioType, domain % ioContext, ierr = err_local)

            ! If to determine if file was opened or not.
            if ( err_local == MPAS_IO_NOERR ) then

               write(stderrUnit, *) ' '
               write(stderrUnit, *) '----- reading dimensions from stream '''//trim(streamName)//''' using file ' &
                                    //trim(streamFilename)

               ! Iterate over list of dimensions we determined we need from the above loop
               call mpas_pool_begin_iteration(streamDimensions)
               do while ( mpas_pool_get_next_member(streamDimensions, poolItr) )
                  if ( poolItr % memberType == MPAS_POOL_DIMENSION ) then
                     ! Try to read the dimension
                     call mpas_io_inq_dim(inputHandle, trim(poolItr % memberName), tempDim, ierr = err_local)

                     ! Check to see if the dimension has already been defined
                     call mpas_pool_get_dimension(readableDimensions, poolItr % memberName, dimValue)

                     ! If to see if dimension was read or not
                     if ( err_local == MPAS_IO_NOERR ) then
                        write(stderrUnit,'(a,a20,a,i8)') '       ', trim(poolItr % memberName), ' =', tempDim

                        if ( .not. associated(dimValue) ) then
                           call mpas_pool_add_dimension(readableDimensions, poolItr % memberName, tempDim)
                        else if ( dimValue /= tempDim .and. dimValue == MPAS_MISSING_DIM ) then
                           dimValue = tempDim
                        else if ( dimValue /= tempDim ) then
                           write(stderrUnit, *) '********************************************************************************'
                           write(stderrUnit, *) 'ERROR: Dimension ' // trim(poolItr % membername) &
                                                // ' was read with an inconsistent value.'
                           call mpas_dmpar_global_abort('********************************************************************************')
                        end if
                     else
                        write(stderrUnit,'(a,a20,a)') '       ', trim(poolItr % memberName), ' *** not found in stream ***'
                     end if

                  end if
               end do

               ! Close file
               call mpas_io_close(inputHandle)
            else
               write(stderrUnit, *) ' '
               write(stderrUnit, *) '  *** unable to open input file '//trim(streamFilename)//' for stream ''' &
                                    //trim(streamName)//''''
            end if

            ! Destroy pool that contains list of streams dimensions
            call mpas_pool_destroy_pool(streamDimensions)

         else if ( .not. streamActive .and. ( streamDirection == MPAS_STREAM_INPUT .or. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) ) then

            write(stderrUnit, *) ' '
            write(stderrUnit, *) '----- skipping inactive stream '''//trim(streamName)//''''
            
         end if

      end do

      write(stderrUnit, *) ' '
      write(stderrUnit, *) '----- done reading dimensions from input streams -----'
      write(stderrUnit, *) ' '
      write(stderrUnit, *) ' '

      call mpas_pool_set_error_level(err_level)

      ! Allocate blocks, and copy indexTo arrays into blocks
      call mpas_block_creator_finalize_block_phase2(domain % streamManager, domain % blocklist, readableDimensions)
      call mpas_link_fields(domain)

      call mpas_pool_destroy_pool(readableDimensions)

   end subroutine mpas_bootstrap_framework_phase2 !}}}


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


   !***********************************************************************
   !
   !  routine get_dimlist_for_field
   !
   !> \brief   Returns an array of dimension names for the specified field
   !> \author  Michael Duda
   !> \date    06 March 2015
   !> \details
   !>  Given a pool containing the field indicated by the 'fieldName' argument,
   !>  this routine allocates the output array, 'dimNames', to the
   !>  dimensionality of the field and fills the array with the dimension names.
   !>  
   !>  If the specified field is a var_array, the dimNames array will have one
   !>  fewer dimension than that of the var_array, and it will be filled with 
   !>  the outer-most dimension names.
   !>
   !>  If the specified field is is inactive (in a package sense) at the time 
   !>  of the call to this routine, the dimNames array will be returned with 
   !>  zero size.
   !
   !-----------------------------------------------------------------------
   subroutine get_dimlist_for_field(allFields, fieldName, dimNames)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: allFields
      character (len=*), intent(in) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: dimNames

      type (mpas_pool_field_info_type) :: fieldInfo
      type (field1DReal), pointer :: real1DField
      type (field2DReal), pointer :: real2DField
      type (field3DReal), pointer :: real3DField
      type (field4DReal), pointer :: real4DField
      type (field5DReal), pointer :: real5DField
      type (field1DInteger), pointer :: int1DField
      type (field2DInteger), pointer :: int2DField
      type (field3DInteger), pointer :: int3DField
      type (field1DChar), pointer :: char1DField


      call mpas_pool_get_field_info(allFields, fieldName, fieldInfo)

      if ( .not. fieldInfo % isActive ) then
         allocate(dimNames(0))
         return
      end if

      if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then

         if ( fieldInfo % nDims == 1 ) then
            call mpas_pool_get_field(allFields, fieldName, real1DField)
            if ( .not. real1DField % isVarArray ) then
               allocate(dimNames(1))
               dimNames(:) = real1DField % dimNames(:)
            else
               allocate(dimNames(0))
            end if

         else if ( fieldInfo % nDims == 2 ) then
            call mpas_pool_get_field(allFields, fieldName, real2DField)
            if ( .not. real2DField % isVarArray ) then
               allocate(dimNames(2))
               dimNames(:) = real2DField % dimNames(:)
            else
               allocate(dimNames(1))
               dimNames(1:1) = real2DField % dimNames(2:2)
            end if

         else if ( fieldInfo % nDims == 3 ) then
            call mpas_pool_get_field(allFields, fieldName, real3DField)
            if ( .not. real3DField % isVarArray ) then
               allocate(dimNames(3))
               dimNames(:) = real3DField % dimNames(:)
            else
               allocate(dimNames(2))
               dimNames(1:2) = real3DField % dimNames(2:3)
            end if

         else if ( fieldInfo % nDims == 4 ) then
            call mpas_pool_get_field(allFields, fieldName, real4DField)
            if ( .not. real4DField % isVarArray ) then
               allocate(dimNames(4))
               dimNames(:) = real4DField % dimNames(:)
            else
               allocate(dimNames(3))
               dimNames(1:3) = real4DField % dimNames(2:4)
            end if

         else if ( fieldInfo % nDims == 5 ) then
            call mpas_pool_get_field(allFields, fieldName, real5DField)
            if ( .not. real5DField % isVarArray ) then
               allocate(dimNames(5))
               dimNames(:) = real5DField % dimNames(:)
            else
               allocate(dimNames(4))
               dimNames(1:4) = real5DField % dimNames(2:5)
            end if

         else
            allocate(dimNames(0))
         end if

      else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then

         if ( fieldInfo % nDims == 1 ) then
            call mpas_pool_get_field(allFields, fieldName, int1DField)
            if ( .not. int1DField % isVarArray ) then
               allocate(dimNames(1))
               dimNames(:) = int1DField % dimNames(:)
            else
               allocate(dimNames(0))
            end if

         else if ( fieldInfo % nDims == 2 ) then
            call mpas_pool_get_field(allFields, fieldName, int2DField)
            if ( .not. int2DField % isVarArray ) then
               allocate(dimNames(2))
               dimNames(:) = int2DField % dimNames(:)
            else
               allocate(dimNames(1))
               dimNames(1:1) = int2DField % dimNames(2:2)
            end if

         else if ( fieldInfo % nDims == 3 ) then
            call mpas_pool_get_field(allFields, fieldName, int3DField)
            if ( .not. int3DField % isVarArray ) then
               allocate(dimNames(3))
               dimNames(:) = int3DField % dimNames(:)
            else
               allocate(dimNames(2))
               dimNames(1:2) = int3DField % dimNames(2:3)
            end if

         else
            allocate(dimNames(0))
         end if

      else if ( fieldInfo % fieldType == MPAS_POOL_CHARACTER ) then

         if ( fieldInfo % nDims == 1 ) then
            call mpas_pool_get_field(allFields, fieldName, char1DField)
            if ( .not. char1DField % isVarArray ) then
               allocate(dimNames(1))
               dimNames(:) = char1DField % dimNames(:)
            else
               allocate(dimNames(0))
            end if

         else
            allocate(dimNames(0))
         end if

      else

         allocate(dimNames(0))

      end if

   end subroutine get_dimlist_for_field!}}}

end module mpas_bootstrapping
