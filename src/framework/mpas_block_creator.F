! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_block_creator
!
!> \brief   This module is responsible for the intial creation and setup of the block data structures.
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!> This module provides routines for the creation of blocks, with both an
!> arbitrary number of blocks per processor and an arbitrary number of halos for
!> each block. The provided routines also setup the exchange lists for each
!> block.
!
!-----------------------------------------------------------------------

module mpas_block_creator

   use mpas_dmpar
   use mpas_block_decomp
   use mpas_hash
   use mpas_sort
   use mpas_derived_types
   use mpas_domain_routines
   use mpas_field_routines
   use mpas_pool_routines
   use mpas_abort, only : mpas_dmpar_global_abort

   contains

!***********************************************************************
!
!  routine mpas_block_creator_setup_blocks_and_0halo_cells
!
!> \brief   Initializes the list of blocks, and determines 0 halo cell indices.
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine sets up the linked list of blocks, and creates the
!>  indexToCellID field for the 0 halo. The information required to setup these
!>  structures is provided as input in cellList, blockID, blockStart, and
!>  blockCount.
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_setup_blocks_and_0halo_cells(nHalos, domain, indexToCellID, cellList, blockID, blockStart, blockCount)!{{{
     integer, intent(in) :: nHalos !< Input: Number of halos for cell fields
     type (domain_type), pointer :: domain !< Input: Domain information
     type (field1dInteger), pointer :: indexToCellID !< Input/Output: indexToCellID field
     integer, dimension(:), intent(in) :: cellList !< Input: List of cell indices owned by this processor
     integer, dimension(:), intent(in) :: blockID !< Input: List of block indices owned by this processor
     integer, dimension(:), intent(in) :: blockStart !< Input: Indices of starting cell id in cellList for each block
     integer, dimension(:), intent(in) :: blockCount !< Input: Number of cells from cellList owned by each block.
 
     type (block_type), pointer :: blockCursor
     type (field1dInteger), pointer :: fieldCursor
 
     integer :: i
     integer :: nBlocks
 
     nBlocks = size(blockID)

     ! Setup first block
     allocate(domain % blocklist)
     nullify(domain % blocklist % prev)
     nullify(domain % blocklist % next)
  
     ! Setup first block field
     allocate(indexToCellID)
     nullify(indexToCellID % next)
 
     ! Loop over blocks
     blockCursor => domain % blocklist
     fieldCursor => indexToCellID
     do i = 1, nBlocks
       ! Initialize block information
       blockCursor % blockID = blockID(i)
       blockCursor % localBlockID = i - 1
       blockCursor % domain => domain
  
       ! Link to block, and setup array size
       fieldCursor % block => blockCursor
       fieldCursor % dimSizes(1) = blockCount(i)
 
       ! Initialize exchange lists
       call mpas_dmpar_init_multihalo_exchange_list(fieldCursor % sendList, nHalos)
       call mpas_dmpar_init_multihalo_exchange_list(fieldCursor % recvList, nHalos)
       call mpas_dmpar_init_multihalo_exchange_list(fieldCursor % copyList, nHalos)
 
       ! Allocate array, and copy indices into array
       allocate(fieldCursor % array(fieldCursor % dimSizes(1)))
       fieldCursor % array(:) = cellList(blockStart(i)+1:blockStart(i)+blockCount(i))
       call mpas_quicksort(fieldCursor % dimSizes(1), fieldCursor % array)
  
       ! Advance cursors, and create new blocks as needed
       if(i < nBlocks) then
         allocate(blockCursor % next)
         allocate(fieldCursor % next)
 
         blockCursor % next % prev => blockCursor

         blockCursor => blockCursor % next
         fieldCursor => fieldCursor % next
       end if
 
       ! Nullify next pointers
       nullify(blockCursor % next)
       nullify(fieldCursor % next)
     end do
   end subroutine mpas_block_creator_setup_blocks_and_0halo_cells!}}}

!***********************************************************************
!
!  routine mpas_block_creator_build_0halo_cell_fields
!
!> \brief   Initializes 0 halo cell based fields requried to work out halos
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine uses the previously setup 0 halo cell field, and the blocks of
!>  data read in by other routhers to determine all of the connectivity for the 0
!>  halo cell fields on all blocks on a processor.
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_build_0halo_cell_fields(nHalos, indexToCellIDBlock, &
              nEdgesOnCellBlock, cellsOnCellBlock, verticesOnCellBlock, edgesOnCellBlock, &
              indexToCellID_0Halo, nEdgesOnCell_0Halo, cellsOnCell_0Halo, &
              verticesOnCell_0Halo, edgesOnCell_0Halo)!{{{
     integer, intent(in) :: nHalos !< Input: Number of halos for cell fields
     type(field1dInteger), pointer :: indexToCellIDBlock !< Input: Block of read in indexToCellID field
     type(field1dInteger), pointer :: nEdgesOnCellBlock !< Input: Block of read in nEdgesOnCell field
     type(field2dInteger), pointer :: cellsOnCellBlock !< Input: Block of read in cellsOnCell field
     type(field2dInteger), pointer :: verticesOnCellBlock !< Input: Block of read in verticesOnCell field
     type(field2dInteger), pointer :: edgesOnCellBlock !< Input: Block of read in edgesOnCellField

     type(field1dInteger), pointer :: indexToCellID_0Halo !< Input: 0-Halo indices for indexToCellID field
     type(field1dInteger), pointer :: nEdgesOnCell_0Halo !< Output: nEdgesOnCell field for 0-Halo cells
     type(field2dInteger), pointer :: cellsOnCell_0Halo !< Output: cellsOnCell field for 0-Halo cells
     type(field2dInteger), pointer :: verticesOnCell_0Halo !< Output: verticesOnCell field for 0-Halo cells
     type(field2dInteger), pointer :: edgesOnCell_0Halo !< Output: edgesOnCell field for 0-Halo cells

     type(field1dInteger), pointer :: indexCursor, nEdgesCursor
     type(field2dInteger), pointer :: cellsOnCellCursor, verticesOnCellCursor, edgesOnCellCursor

     integer, dimension(:), pointer :: sendingHaloLayers

     integer :: nCellsInBlock, maxEdges

     ! Only sending from halo layer 1 for setup
     allocate(sendingHaloLayers(1))
     sendingHaloLayers(1) = 1

     maxEdges = cellsOnCellBlock % dimSizes(1)

     ! Build exchange list from the block of read in data to each block's index fields.
     call mpas_dmpar_get_exch_list(1, indexToCellIDBlock, indexToCellID_0Halo)

     ! Setup header fields if at least 1 block exists
     allocate(nEdgesOnCell_0Halo)
     nullify(nEdgesOncell_0Halo % next)

     allocate(cellsOnCell_0Halo)
     nullify(cellsOnCell_0Halo % next)
  
     allocate(verticesOnCell_0Halo)
     nullify(verticesOnCell_0Halo % next)
  
     allocate(edgesOnCell_0Halo)
     nullify(edgesOnCell_0Halo % next)

     ! Loop over blocks
     indexCursor => indexToCellID_0Halo
     nEdgesCursor => nEdgesOnCell_0Halo
     cellsOnCellCursor => cellsOnCell_0Halo
     verticesOnCellCursor => verticesOnCell_0Halo
     edgesOnCellCursor => edgesOnCell_0Halo
     do while(associated(indexCursor))
       nCellsInBlock = indexCursor % dimSizes(1)

       ! Link to block structure
       nEdgesCursor % block => indexCursor % block
       cellsOnCellCursor % block => indexCursor % block
       verticesOnCellCursor % block => indexCursor % block
       edgesOnCellCursor % block => indexCursor % block

       ! Setup array sizes
       nEdgesCursor % dimSizes(1) = nCellsInBlock
       cellsOnCellCursor % dimSizes(1) = maxEdges
       cellsOnCellCursor % dimSizes(2) = nCellsInBlock
       verticesOnCellCursor % dimSizes(1) = maxEdges
       verticesOnCellCursor % dimSizes(2) = nCellsInBlock
       edgesOnCellCursor % dimSizes(1) = maxEdges
       edgesOnCellCursor % dimSizes(2) = nCellsInBlock

       ! Link exchange lists
       nEdgesCursor % sendList => indexCursor % sendList
       nEdgesCursor % recvList => indexCursor % recvList
       nEdgesCursor % copyList => indexCursor % copyList
       cellsOnCellCursor % sendList => indexCursor % sendList
       cellsOnCellCursor % recvList => indexCursor % recvList
       cellsOnCellCursor % copyList => indexCursor % copyList
       verticesOnCellCursor % sendList => indexCursor % sendList
       verticesOnCellCursor % recvList => indexCursor % recvList
       verticesOnCellCursor % copyList => indexCursor % copyList
       edgesOnCellCursor % sendList => indexCursor % sendList
       edgesOnCellCursor % recvList => indexCursor % recvList
       edgesOnCellCursor % copyList => indexCursor % copyList

       ! Allocate arrays
       allocate(nEdgesCursor % array(nEdgesCursor % dimSizes(1)))
       allocate(cellsOnCellCursor % array(cellsOnCellCursor % dimSizes(1), cellsOnCellCursor % dimSizes(2)))
       allocate(verticesOnCellCursor % array(verticesOnCellCursor % dimSizes(1), verticesOnCellCursor % dimSizes(2)))
       allocate(edgesOnCellCursor % array(edgesOnCellCursor % dimSizes(1), edgesOnCellCursor % dimSizes(2)))
       
       ! Create new blocks and advance cursors as needed
       indexCursor => indexCursor % next
       if(associated(indexCursor)) then
         allocate(nEdgesCursor % next)
         allocate(cellsOnCellCursor % next)
         allocate(verticesOnCellCursor % next)
         allocate(edgesOnCellCursor % next)

         nEdgesCursor => nEdgesCursor % next
         cellsOnCellCursor => cellsOnCellCursor % next
         verticesOnCellCursor => verticesOnCellCursor % next
         edgesOnCellCursor => edgesOnCellCursor % next

       end if

       ! Nullify next pointers
       nullify(nEdgesCursor % next)
       nullify(cellsOnCellCursor % next)
       nullify(verticesOnCellCursor % next)
       nullify(edgesOnCellCursor % next)
     end do ! indexCursor loop over blocks

     ! Communicate data from read in blocks to each block's fields
     call mpas_dmpar_alltoall_field(nEdgesOnCellBlock, nEdgesOnCell_0Halo, sendingHaloLayers)
     call mpas_dmpar_alltoall_field(cellsOnCellBlock, cellsOnCell_0Halo, sendingHaloLayers)
     call mpas_dmpar_alltoall_field(verticesOnCellBlock, verticesOnCell_0Halo, sendingHaloLayers)
     call mpas_dmpar_alltoall_field(edgesOnCellBlock, edgesOnCell_0Halo, sendingHaloLayers)
   end subroutine mpas_block_creator_build_0halo_cell_fields!}}}

!***********************************************************************
!
!  routine mpas_block_creator_build_0_and_1halo_edge_fields
!
!> \brief   Initializes 0 and 1 halo edge based fields requried to work out halos
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine uses the previously setup 0 halo cell fields, and the blocks of
!>  data read in by other routhers to determine which edges are in a blocks
!>  0 and 1 halo for all blocks on a processor.
!>  NOTE: This routine can be used on either edges or vertices
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_build_0_and_1halo_edge_fields(nHalos, indexToEdgeIDBlock, cellsOnEdgeBlock, indexToCellID_0Halo, nEdgesOnCell_0Halo, edgesOnCell_0Halo, indexToEdgeID_0Halo, cellsOnEdge_0Halo, nEdgesSolve)!{{{
     integer, intent(in) :: nHalos !< Input: Number of halos for cell fields
     type (field1dInteger), pointer :: indexToEdgeIDBlock !< Input: indexToEdgeID read in field
     type (field2dInteger), pointer :: cellsOnEdgeBlock !< Input: cellsOnEdge read in field
     type (field1dInteger), pointer :: indexToCellID_0Halo !< Input: indexToCellID field on 0 halo
     type (field1dInteger), pointer :: nEdgesOnCell_0Halo !< Input: nEdgesOnCell field on 0 halo
     type (field2dInteger), pointer :: edgesOnCell_0Halo !< Input: edgesOnCell field on 0 and 1 halos
     type (field1dInteger), pointer :: indexToEdgeID_0Halo !< Output: indexToEdgeID field on 0 and 1 halos
     type (field2dInteger), pointer :: cellsOnEdge_0Halo !< Output: CellsOnEdge field on 0 and 1 halos
     type (field1dInteger), pointer :: nEdgesSolve !< Output: Array with max index to edges in halos

     type (field0dInteger), pointer :: offSetField, edgeLimitField
     type (field1dInteger), pointer :: haloIndices

     type (field0dInteger), pointer :: offSetCursor, edgeLimitCursor
     type (field1dInteger), pointer :: indexToCellCursor, indexToEdgeCursor, nEdgesCursor, haloCursor, nEdgesSolveCursor
     type (field2dInteger), pointer :: edgesOnCellCursor, cellsOnEdgeCursor

     integer, dimension(:), pointer :: localEdgeList
     integer, dimension(:), pointer :: sendingHaloLayers
     integer :: nEdgesLocal, nCellsInBlock, maxEdges, edgeDegree
     integer :: haloStart

     ! Setup sendingHaloLayers
     allocate(sendingHaloLayers(1))
     sendingHaloLayers(1) = 1

     ! Get dimension information
     maxEdges = edgesOnCell_0Halo % dimSizes(1)
     edgeDegree = cellsOnEdgeBlock % dimSizes(1)

     ! Setup initial block for each field
     allocate(cellsOnEdge_0Halo)
     allocate(indexToEdgeID_0Halo)

     nullify(cellsOnEdge_0Halo % next)
     nullify(indexToEdgeID_0Halo % next)

     ! Loop over blocks
     indexToCellCursor => indexToCellID_0Halo
     edgesOnCellCursor => edgesOnCell_0Halo
     nEdgesCursor => nEdgesOnCell_0Halo
     indexToEdgeCursor => indexToEdgeID_0Halo
     cellsOnEdgeCursor => cellsOnEdge_0Halo
     do while(associated(indexToCellCursor))
       ! Determine number of cells in block
       nCellsInBlock = indexToCellCursor % dimSizes(1)

       ! Determine all edges in block
       call mpas_block_decomp_all_edges_in_block(maxEdges, nCellsInBlock, nEdgesCursor % array, edgesOnCellCursor % array, nEdgesLocal, localEdgeList)

       ! Setup indexToEdge block
       indexToEdgeCursor % block => indexToCellCursor % block
       indexToEdgeCursor % dimSizes(1) = nEdgesLocal
       allocate(indexToEdgeCursor % array(indexToEdgeCursor % dimSizes(1)))
       indexToEdgeCursor % array(:) = localEdgeList(:)

       ! Setup cellsOnEdge block
       cellsOnEdgeCursor % block => indexToCellCursor % block
       cellsOnEdgeCursor % dimSizes(1) = edgeDegree
       cellsOnEdgeCursor % dimSizes(2) = nEdgesLocal
       allocate(cellsOnEdgeCursor % array(cellsOnEdgeCursor % dimSizes(1), cellsOnEdgeCursor % dimSizes(2)))

       ! Setup exchange lists
       call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeCursor % sendList, nHalos+1)
       call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeCursor % recvList, nHalos+1)
       call mpas_dmpar_init_multihalo_exchange_list(indexToEdgeCursor % copyList, nHalos+1)

       ! Link exchange lists
       cellsOnEdgeCursor % sendList => indexToEdgeCursor % sendList
       cellsOnEdgeCursor % recvList => indexToEdgeCursor % recvList
       cellsOnEdgeCursor % copyList => indexToEdgeCursor % copyList
       
       ! Remove localEdgeList array
       deallocate(localEdgeList)

       ! Advance cursors, and create new blocks if needed
       indexToCellCursor => indexToCellCursor % next
       edgesOnCellCursor => edgesOnCellCursor % next
       nEdgescursor => nEdgesCursor % next
       if(associated(indexToCellCursor)) then
         allocate(indexToEdgeCursor % next)
         indexToEdgeCursor => indexToEdgeCursor % next

         allocate(cellsOnEdgeCursor % next)
         cellsOnEdgeCursor => cellsOnEdgeCursor % next
       end if

       ! Nullify next pointers
       nullify(indexToEdgeCursor % next)
       nullify(cellsOnEdgeCursor % next)
     end do ! indexToCursor loop over blocks

     ! Build exchangel ists from read in blocks to owned blocks.
     call mpas_dmpar_get_exch_list(1, indexToEdgeIDBlock, indexToEdgeID_0Halo)

     ! Perform all to all to get owned block data
     call mpas_dmpar_alltoall_field(cellsOnEdgeBlock, cellsOnEdge_0Halo, sendingHaloLayers)

     ! Setup first block's fields if there is at least 1 block.
     if(associated(indexToEdgeID_0Halo)) then
       allocate(haloIndices)
       allocate(offSetField)
       allocate(edgeLimitField)
       allocate(nEdgesSolve)
     else
       nullify(haloIndices)
       nullify(offSetField)
       nullify(edgeLimitField)
       nullify(nEdgesSolve)
     end if

     ! Loop over blocks
     indexToEdgeCursor => indexToEdgeID_0Halo
     cellsOnEdgeCursor => cellsOnEdge_0Halo
     indexToCellCursor => indexToCellID_0Halo
     haloCursor => haloIndices
     offSetCursor => offSetField
     edgeLimitCursor => edgeLimitField
     nEdgesSolveCursor => nEdgesSolve
     do while(associated(indexToEdgeCursor))
       ! Determine 0 and 1 halo edges
       call mpas_block_decomp_partitioned_edge_list(indexToCellCursor % dimSizes(1), indexToCellCursor % array, &
                                                    edgeDegree, indexToEdgeCursor % dimSizes(1), cellsOnEdgeCursor % array, &
                                                    indexToEdgeCursor % array, haloStart)

       ! Link blocks                                                
       haloCursor % block => indexToEdgeCursor % block
       offSetCursor % block => indexToEdgeCursor % block
       edgeLimitCursor % block => indexToEdgeCursor % block
       nEdgesSolveCursor % block => indexToEdgeCursor % block

       ! Setup haloIndices
       haloCursor % dimSizes(1) = indexToEdgeCursor % dimSizes(1) - (haloStart-1)
       allocate(haloCursor % array(haloCursor % dimSizes(1)))
       haloCursor % array(:) = indexToEdgeCursor % array(haloStart:indexToEdgeCursor % dimSizes(1))

       ! Link exchange lists
       haloCursor % sendList => indexToEdgeCursor % sendList
       haloCursor % recvList => indexToEdgeCursor % recvList
       haloCursor % copyList => indexToEdgeCursor % copyList

       ! Determine offSet and limit on 0 halo edges for exchange list creation
       offSetCursor % scalar = haloStart - 1
       edgeLimitCursor % scalar = haloStart - 1

       ! Setup nEdgesSolve
       nEdgesSolveCursor % dimSizes(1) = nHalos+2 
       allocate(nEdgesSolveCursor % array(nEdgesSolve % dimSizes(1)))
       nEdgesSolveCursor % array = -1
       nEdgesSolveCursor % array(1) = haloStart - 1
       nEdgesSolveCursor % array(2) = indexToEdgeCursor % dimSizes(1)

       ! Advance cursors, and create new blocks if needed
       indexToEdgeCursor => indexToEdgeCursor % next
       cellsOnEdgeCursor => cellsOnEdgeCursor % next
       indexToCellCursor => indexToCellCursor % next
       if(associateD(indexToEdgeCursor)) then
         allocate(haloCursor % next)
         haloCursor => haloCursor % next

         allocate(offSetcursor % next)
         offSetCursor => offSetCursor % next

         allocate(edgeLimitCursor % next)
         edgeLimitCursor => edgeLimitCursor % next

         allocate(nEdgesSolveCursor % next)
         nEdgesSolveCursor => nEdgesSolveCursor % next
       end if

       ! Nullify next pointers
       nullify(haloCursor % next)
       nullify(offSetCursor % next)
       nullify(edgeLimitCursor % next)
       nullify(nEdgesSolveCursor % next)
     end do

     ! Create exchange lists from 0 halo to 1 haloedges 
     call mpas_dmpar_get_exch_list(1, indexToEdgeID_0Halo, haloIndices, offSetField, edgeLimitField)

     ! Deallocate fields that are not needed anymore.
     call mpas_deallocate_field(haloIndices)
     call mpas_deallocate_field(offSetField)
     call mpas_deallocate_field(edgeLimitCursor)
     deallocate(sendingHaloLayers)

   end subroutine mpas_block_creator_build_0_and_1halo_edge_fields!}}}

!***********************************************************************
!
!  routine mpas_block_creator_build_cell_halos
!
!> \brief   Builds cell halos
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine uses the previously setup 0 halo cell fields to determine
!>  which cells fall in each halo layer for a block. During this process, each
!>  halo's exchange lists are created. This process is performed for all blocks on
!>  a processor.
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_build_cell_halos(nHalos, indexToCellID, nEdgesOnCell, cellsOnCell, verticesOnCell, edgesOnCell, nCellsSolve)!{{{
     integer, intent(in) :: nHalos !< Input: Number of halos for cell fields
     type (field1dInteger), pointer :: indexToCellID !< Input/Output: indexToCellID field for all halos
     type (field1dInteger), pointer :: nEdgesOnCell !< Input/Output: nEdgesOnCell field for all halos
     type (field2dInteger), pointer :: cellsOnCell !< Input/Output: cellsOnCell field for all halos
     type (field2dInteger), pointer :: verticesOnCell !< Input/Output: verticesOnCell field for all halos
     type (field2dInteger), pointer :: edgesOnCell !< Input/Output: edgesOnCell field for all halos
     type (field1dInteger), pointer :: nCellsSolve !< Output: Field with indices to end of each halo

     type (dm_info), pointer :: dminfo

     type (field1dInteger), pointer :: haloIndices

     type (field0dInteger), pointer :: offSetCursor, cellLimitCursor
     type (field1dInteger), pointer :: indexCursor, nEdgesCursor, haloCursor, nCellsSolveCursor
     type (field2dInteger), pointer :: cellsOnCellCursor, verticesOnCellCursor, edgesOnCellCursor

     type (field0dInteger), pointer :: offSetField
     type (field0dInteger), pointer :: cellLimitField

     integer, dimension(:), pointer :: sendingHaloLayers
     integer, dimension(:), pointer :: field1dArrayHolder
     integer, dimension(:,:), pointer :: field2dArrayHolder

     type (graph), pointer :: blockGraph, blockGraphWithHalo

     integer :: nCellsInBlock, nCellsInHalo, maxEdges
     integer :: iHalo

     dminfo => indexToCellID % block % domain % dminfo
     allocate(sendingHaloLayers(1))

     ! Setup header fields
     allocate(nCellsSolve)
     allocate(cellLimitField)
     allocate(offSetField)

     nullify(nCellsSolve % next)
     nullify(cellLimitField % next)
     nullify(offSetField % next)

     ! Loop over blocks
     offSetCursor => offsetField
     cellLimitCursor => cellLimitField
     indexCursor => indexToCellID
     nCellsSolveCursor => nCellsSolve
     do while (associated(indexCursor))
       ! Setup offset
       offSetCursor % scalar = indexCursor % dimSizes(1)
       offSetCursor % block => indexCursor % block

       ! Setup nCellsSolve
       nCellsSolveCursor % dimSizes(1) = nHalos+1
       allocate(nCellsSolveCursor % array(nCellsSolveCursor % dimSizes(1)))
       nCellsSolveCursor % array(1) = indexCursor % dimSizes(1)
       nCellsSolveCursor % block => indexCursor % block

       ! Setup owned cellLimit
       cellLimitCursor % scalar = indexCursor % dimSizes(1)
       cellLimitCursor % block => indexCursor % block

       ! Advance cursors and create new blocks if needed
       indexCursor => indexCursor % next
       if(associated(indexCursor)) then
         allocate(offSetCursor % next)
         offSetCursor => offSetCursor % next

         allocate(nCellsSolveCursor % next)
         nCellsSolveCursor => nCellsSolveCursor % next

         allocate(cellLimitCursor % next)
         cellLimitCursor => cellLimitCursor % next
       end if

       ! Nullify next pointers
       nullify(offSetCursor % next)
       nullify(nCellssolveCursor % next)
       nullify(cellLimitCursor % next)
     end do

     ! Loop over halos
     do iHalo = 1, nHalos
       ! Sending halo layer is the current halo
       sendingHaloLayers(1) = iHalo

       if(associated(indexToCellID)) then
         allocate(haloIndices)
         nullify(haloIndices % next)
       else
         nullify(haloIndices)
       end if

       ! Loop over blocks
       indexCursor => indexToCellID
       nEdgesCursor => nEdgesOnCell
       cellsOnCellCursor => cellsOnCell
       verticesOnCellCursor => verticesOnCell
       edgesOnCellCursor => edgesOnCell
       haloCursor => haloIndices
       offSetCursor => offSetField
       do while(associated(indexCursor))
         ! Determine block dimensions
         nCellsInBlock = indexCursor % dimSizes(1)
         maxEdges = cellsOnCellCursor % dimSizes(1)

         ! Setup offSet
         offSetCursor % scalar = nCellsInBlock 

         ! Setup block graphs
         allocate(blockGraphWithHalo)
         allocate(blockGraph)
         allocate(blockGraph % vertexID(nCellsInBlock))
         allocate(blockGraph % nAdjacent(nCellsInBlock))
         allocate(blockGraph % adjacencyList(maxEdges, nCellsInBlock))

         blockGraph % nVertices = nCellsInBlock
         blockGraph % nVerticesTotal = nCellsInBlock
         blockGraph % maxDegree = maxEdges
         blockGraph % ghostStart = nCellsInBlock + 1

         blockGraph % vertexID(:) = indexCursor % array(:)
         blockGraph % nAdjacent(:) = nEdgesCursor % array(:)
         blockGraph % adjacencyList(:,:) = cellsOnCellCursor % array(:,:)

         ! Determine all cell id's with the next halo added
         call mpas_block_decomp_add_halo(dminfo, blockGraph, blockGraphWithHalo)

         ! Setup haloIndices
         haloCursor % dimSizes(1) = blockGraphWithHalo % nVerticesTotal - blockGraphWithHalo % nVertices
         allocate(haloCursor % array(haloCursor % dimSizes(1)))
         haloCursor % array(:) = blockGraphWithHalo % vertexID(blockGraphWithHalo % nVertices+1:blockGraphWithHalo % nVerticesTotal)
         call mpas_quicksort(haloCursor % dimSizes(1), haloCursor % array)
         haloCursor % sendList => indexCursor % sendList
         haloCursor % recvList => indexCursor % recvList
         haloCursor % copyList => indexCursor % copyList
         haloCursor % block => indexCursor % block

         ! Deallocate block graphs
         deallocate(blockGraphWithHalo % vertexID)
         deallocate(blockGraphWithHalo % nAdjacent)
         deallocate(blockGraphWithHalo % adjacencyList)
         deallocate(blockGraphWithHalo)

         deallocate(blockGraph % vertexID)
         deallocate(blockGraph % nAdjacent)
         deallocate(blockGraph % adjacencyList)
         deallocate(blockGraph)

         ! Advance cursors and create new block if needed
         indexCursor => indexCursor % next
         nEdgesCursor => nEdgesCursor % next
         cellsOnCellCursor => cellsOnCellCursor % next
         verticesOnCellCursor => verticesOnCellCursor % next
         edgesOnCellCursor => edgesOnCellCursor % next
         offSetCursor => offSetCursor % next
         if(associated(indexCursor)) then
           allocate(haloCursor % next)
           haloCursor => haloCursor % next
         end if
         ! Nullify next pointer
         nullify(haloCursor % next)
       end do ! indexCursor loop over blocks

       ! Create exchange lists for current halo layer
       call mpas_dmpar_get_exch_list(iHalo, indexToCellID, haloIndices, offSetField, cellLimitField)

       ! Loop over blocks
       indexCursor => indexToCellID
       nEdgesCursor => nEdgesOnCell
       cellsOnCellCursor => cellsOnCell
       verticesOnCellCursor => verticesOnCell
       edgesOnCellCursor => edgesOnCell
       haloCursor => haloIndices
       nCellsSolveCursor => nCellsSolve
       do while(associated(indexCursor))
         ! Determine block dimensions
         nCellsInBlock = indexCursor % dimSizes(1)
         nCellsInHalo = haloCursor % dimSizes(1) 

         ! Setup new layer's nCellsSolve
         nCellsSolveCursor % array(iHalo+1) = nCellsInBlock + nCellsInHalo

         ! Copy cell indices into indexToCellID field
         field1dArrayHolder => indexCursor % array
         indexCursor % dimSizes(1) = nCellsSolveCursor % array(iHalo+1)
         allocate(indexCursor % array(indexCursor % dimSizes(1)))
         indexCursor % array(1:nCellsInBlock) = field1dArrayHolder(:)
         indexCursor % array(nCellsInBlock+1:nCellsSolveCursor % array(iHalo+1)) = haloCursor % array(1:nCellsInHalo)
         deallocate(field1dArrayHolder)

         ! Allocate space in nEdgesOnCell
         field1dArrayHolder => nEdgesCursor % array
         nEdgesCursor % dimSizes(1) = nCellsSolveCursor % array(iHalo+1)
         allocate(nEdgesCursor % array(nEdgesCursor % dimSizes(1)))
         nEdgesCursor % array = -1
         nEdgesCursor % array(1:nCellsInBlock) = field1dArrayHolder(:)
         deallocate(field1dArrayHolder)

         ! Allocate space in cellsOnCell
         field2dArrayHolder => cellsOnCellCursor % array
         cellsOnCellCursor  % dimSizes(2) = nCellsSolveCursor % array(iHalo+1)
         allocate(cellsOnCellCursor % array(cellsOnCellCursor % dimSizes(1), cellsOnCellCursor % dimSizes(2)))
         cellsOnCellCursor % array = -1
         cellsOnCellCursor % array(:,1:nCellsInBlock) = field2dArrayHolder(:,:)
         deallocate(field2dArrayHolder)

         ! Allocate space in verticesOnCell
         field2dArrayHolder => verticesOnCellCursor % array
         verticesOnCellCursor  % dimSizes(2) = nCellsSolveCursor % array(iHalo+1)
         allocate(verticesOnCellCursor % array(verticesOnCellCursor % dimSizes(1), verticesOnCellCursor % dimSizes(2)))
         verticesOnCellCursor % array = -1
         verticesOnCellCursor % array(:,1:nCellsInBlock) = field2dArrayHolder(:,:)
         deallocate(field2dArrayHolder)

         ! Allocate space in edgesOnCell
         field2dArrayHolder => edgesOnCellCursor % array
         edgesOnCellCursor  % dimSizes(2) = nCellsSolveCursor % array(iHalo+1)
         allocate(edgesOnCellCursor % array(edgesOnCellCursor % dimSizes(1), edgesOnCellCursor % dimSizes(2)))
         edgesOnCellCursor % array = -1
         edgesOnCellCursor % array(:,1:nCellsInBlock) = field2dArrayHolder(:,:)
         deallocate(field2dArrayHolder)
        
         indexCursor => indexCursor % next
         nEdgesCursor => nEdgesCursor % next
         cellsOnCellCursor => cellsOnCellCursor % next
         verticesOnCellCursor => verticesOnCellCursor % next
         edgesOnCellCursor => edgesOnCellCursor % next
         haloCursor => haloCursor % next
         nCellsSolveCursor => nCellsSolveCursor % next
       end do

       ! Perform allToAll communications
       call mpas_dmpar_alltoall_field(indexToCellID, indexToCellID, sendingHaloLayers)
       call mpas_dmpar_alltoall_field(nEdgesOnCell, nEdgesOncell, sendingHaloLayers)
       call mpas_dmpar_alltoall_field(cellsOnCell, cellsOnCell, sendingHaloLayers)
       call mpas_dmpar_alltoall_field(verticesOnCell, verticesOnCell, sendingHaloLayers)
       call mpas_dmpar_alltoall_field(edgesOnCell, edgesOnCell, sendingHaloLayers)

       ! Deallocate haloindices field
       call mpas_deallocate_field(haloIndices)
     end do ! iHalo loop over nHalos

     ! Deallocate array and field.
     deallocate(sendingHaloLayers)
     call mpas_deallocate_field(offSetField)

   end subroutine mpas_block_creator_build_cell_halos!}}}

!***********************************************************************
!
!  routine mpas_block_creator_build_edge_halos
!
!> \brief   Builds edge halos
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine uses the previously setup 0 and 1 edge fields and 0 halo cell fields to determine
!>  which edges fall in each halo layer for a block. During this process, each
!>  halo's exchange lists are created. This process is performed for all blocks on
!>  a processor. 
!>  NOTE: This routine can be used on either edges or edges
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_build_edge_halos(nHalos, indexToCellID, nEdgesOnCell, nCellsSolve, edgesOnCell, indexToEdgeID, cellsOnEdge, nEdgesSolve)!{{{
     integer, intent(in) :: nHalos !< Input: Number of halos for cell fields
     type (field1dInteger), pointer :: indexToCellID !< Input: indexToCellID field for all halos
     type (field1dInteger), pointer :: nEdgesOnCell !< Input: nEdgesOnCell field for all halos
     type (field1dInteger), pointer :: nCellsSolve !< Input: nCellsSolve field for all halos
     type (field2dInteger), pointer :: edgesOnCell !< Input/Output: edgesOnCell field for all halos
     type (field1dInteger), pointer :: indexToEdgeID !< Input/Output: indexToEdgeID field for halos 0 and 1, but output for all halos
     type (field2dInteger), pointer :: cellsOnEdge !< Output: cellsOnEdge field for all halos
     type (field1dInteger), pointer :: nEdgesSolve !< Input/Output: nEdgesSolve field for halos 0 and 1, but output for all halos

     type (field0dInteger), pointer :: offSetField, edgeLimitField
     type (field1dInteger), pointer :: haloIndices

     type (field0dInteger), pointer :: offSetCursor, edgeLimitCursor
     type (field1dInteger), pointer :: nEdgesCursor, nCellsSolveCursor, indexToEdgeCursor, nEdgesSolveCursor, haloCursor
     type (field2dInteger), pointer :: edgesOnCellCursor, cellsOnEdgeCursor

     integer, dimension(:), pointer :: sendingHaloLayers
     integer, dimension(:), pointer :: array1dHolder, localEdgeList
     integer, dimension(:,:), pointer :: array2dHolder

     integer :: iHalo, iBlock, i, j
     integer :: nBlocks, nCellsInBlock, nEdgesLocal, haloSize
     integer :: maxEdges, edgeDegree

     type (hashtable), dimension(:), pointer :: edgeList

     ! Determine dimensions
     maxEdges = edgesOnCell % dimSizes(1)
     edgeDegree = cellsOnEdge % dimSizes(1)

     ! Allocate some needed arrays and fields
     allocate(sendingHaloLayers(1))

     allocate(haloIndices)
     allocate(offSetField)
     allocate(edgeLimitField)

     nullify(haloIndices % next)
     nullify(offSetField % next)
     nullify(edgeLimitField % next)

     ! Determine number of blocks, and setup field lists
     ! Loop over blocks
     nBlocks = 0
     indexToEdgeCursor => indexToEdgeID
     haloCursor => haloIndices
     offSetCursor => offSetField
     edgeLimitCursor => edgeLimitField
     nEdgesSolveCursor => nEdgesSolve
     do while(associated(indexToEdgeCursor))
       nBlocks = nBlocks + 1

       ! Setup edgeLimit and offSet
       edgeLimitCursor % scalar = nEdgesSolveCursor % array(1)
       offSetCursor % scalar = nEdgesSolveCursor % array(2)

       ! Link blocks
       edgeLimitCursor % block => indexToEdgeCursor % block
       offSetCursor % block => indexToEdgeCursor % block
       haloCursor % block => indexToEdgeCursor % block

       ! Link exchange lists
       haloCursor % sendList => indexToEdgeCursor % sendList
       haloCursor % recvList => indexToEdgeCursor % recvList
       haloCursor % copyList => indexToEdgeCursor % copyList

       ! Advance cursors and create new blocks if needed
       indexToEdgeCursor => indexToEdgeCursor % next
       nEdgesSolveCursor => nEdgesSolveCursor % next
       if(associated(indexToEdgeCursor)) then
         allocate(haloCursor % next)
         haloCursor => haloCursor % next

         allocate(offSetCursor % next)
         offSetCursor => offSetCursor % next

         allocate(edgeLimitCursor % next)
         edgeLimitCursor =>edgeLimitCursor % next
       end if

       ! Nullify next pointers
       nullify(haloCursor % next)
       nullify(offSetCursor % next)
       nullify(edgeLimitCursor % next)
     end do

     ! Allocate and initialize hashtables
     allocate(edgeList(nBlocks))
     do iBlock = 1, nBlocks
       call mpas_hash_init(edgeList(iBlock))
     end do

     ! Build unique 0 and 1 halo list for each block
     indexToEdgeCursor => indexToEdgeID
     do while(associated(indexToEdgeCursor))
       iBlock = indexToEdgeCursor % block % localBlockID + 1

       do i = 1, indexToEdgeCursor % dimSizes(1)
         if(.not. mpas_hash_search(edgeList(iBlock), indexToEdgeCursor % array(i))) then
           call mpas_hash_insert(edgeList(iBlock), indexToEdgeCursor % array(i))
         end if
       end do

       indexToEdgeCursor => indexToEdgeCursor % next
     end do

     ! Append new unique edge id's to indexToEdgeID field.
     do iHalo = 3, nHalos+2
       sendingHaloLayers(1) = iHalo-1

       ! Loop over blocks
       indexToEdgeCursor => indexToEdgeID
       nEdgesCursor => nEdgesOnCell
       nCellsSolveCursor => nCellsSolve
       edgesOnCellCursor => edgesOnCell
       nEdgesSolveCursor => nEdgesSolve
       haloCursor => haloIndices
       offSetCursor => offSetField
       do while(associated(indexToEdgeCursor))
         iBlock = indexToEdgeCursor % block % localBlockID+1
         nCellsInBlock = nCellsSolveCursor % array(iHalo-1)
         offSetCursor % scalar = nEdgesSolveCursor % array(iHalo-1)
  
         ! Determine all edges in block
         call mpas_block_decomp_all_edges_in_block(maxEdges, nCellsInBlock, nEdgesCursor % array, edgesOnCellCursor % array, nEdgesLocal, localEdgeList)

         nEdgesSolveCursor % array(iHalo) = nEdgesLocal
         haloSize = nEdgesLocal - nEdgesSolveCursor % array(iHalo-1)
         haloCursor % dimSizes(1) = haloSize

         allocate(haloCursor % array(haloCursor % dimSizes(1)))

         ! Add all edges into block, and figure out which are new edges meaning they belong to the new halo layer
         j = 1
         do i = 1, nEdgesLocal
           if(.not. mpas_hash_search(edgeList(iBlock), localEdgeList(i))) then
             call mpas_hash_insert(edgeList(iBlock), localEdgeList(i))
             haloCursor % array(j) = localEdgeList(i)
             j = j + 1
           end if
         end do

         deallocate(localEdgeList)

         ! Advance Cursors
         indexToEdgeCursor => indexToEdgeCursor % next
         nEdgesCursor => nEdgesCursor % next
         nCellsSolveCursor => nCellsSolveCursor % next
         edgesOnCellCursor => edgesOnCellCursor % next
         nEdgesSolveCursor => nEdgesSolveCursor % next
         haloCursor => haloCursor % next
         offSetCursor => offSetCursor % next
       end do

       ! Build current layers exchange list
       call mpas_dmpar_get_exch_list(iHalo-1, indexToEdgeID, haloIndices, offSetField, edgeLimitField)

       ! Loop over blocks
       indexToEdgeCursor => indexToEdgeID
       cellsOnEdgeCursor => cellsOnEdge
       nEdgesSolveCursor => nEdgesSolve
       haloCursor => haloIndices
       do while(associated(indexToEdgeCursor))
         ! Copy in new halo indices
         array1dHolder => indexToEdgeCursor % array
         indexToEdgeCursor % dimSizes(1) = nEdgesSolveCursor % array(iHalo)
         allocate(indexToEdgeCursor % array(indexToEdgeCursor % dimSizes(1)))
         indexToEdgeCursor % array(1:nEdgesSolveCursor % array(iHalo-1)) = array1dHolder(:)
         indexToEdgeCursor % array(nEdgesSolveCursor % array(iHalo-1)+1:nEdgesSolveCursor % array(iHalo)) = haloCursor % array(:)
         deallocate(array1dHolder)

         ! Allocate space in cellsOnEdge
         array2dHolder => cellsOnEdgeCursor % array
         cellsOnEdgeCursor % dimSizes(2) = nEdgesSolveCursor % array(iHalo)
         allocate(cellsOnEdgeCursor % array(cellsOnEdgeCursor % dimSizes(1), cellsOnEdgeCursor % dimSizes(2)))
         cellsOnEdgeCursor % array(:,1:nEdgesSolveCursor % array(iHalo-1)) = array2dHolder(:,:)
         deallocate(array2dHolder)

         ! Deallocate haloCursor array
         deallocate(haloCursor % array)

         ! Advance cursors
         indexToEdgeCursor => indexToEdgeCursor % next
         cellsOnEdgeCursor => cellsOnEdgeCursor % next
         nEdgesSolveCursor => nEdgesSolveCursor % next
         haloCursor => haloCursor % next
       end do

       ! Performe allToAll communication
       call mpas_dmpar_alltoall_field(cellsOnEdge, cellsOnEdge, sendingHaloLayers)
     end do

     ! Deallocate fields, hashtables, and arrays
     call mpas_deallocate_field(haloIndices)
     call mpas_deallocate_field(edgeLimitField)
     call mpas_deallocate_field(offSetField)
     do iBlock=1,nBlocks
       call mpas_hash_destroy(edgeList(iBlock))
     end do
     deallocate(edgeList)
     deallocate(sendingHaloLayers)


   end subroutine mpas_block_creator_build_edge_halos!}}}

!***********************************************************************
!
!  routine mpas_block_creator_finalize_block_phase1
!
!> \brief   Phase 1 of block creation finalization
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine finalizes the block initialization processor. It calls
!>  mpas_block_allocate to allocate space for all fields in a block. Then the 0
!>  halo indices for each element and the exchange lists are copied into the
!>  appropriate block. A halo update is required after this routine is called
!>  to make sure all data in a block is valid.
!>
!>  Fields and non-framework required dimensions are not setup after this
!>  routine returns. To complete setup, a call to
!>  mpas_block_creator_finalize_block_phase2 must follow.
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_finalize_block_phase1(nHalos, blocklist, nCellsSolve, & !{{{
                                                       nEdgesSolve, nVerticesSolve, vertexDegree, maxEdges, &
                                                       indexToCellID, indexToEdgeID, &
                                                       indexToVertexID)
     integer, intent(in) :: nHalos !< Input: Number of halos for cell fields
     type (block_type), pointer :: blocklist !< Input/Output: Linked List of blocks
     type (field1dInteger), pointer :: nCellsSolve !< Input: nCellsSolve field information
     type (field1dInteger), pointer :: nEdgesSolve !< Input: nEdgesSolve field information
     type (field1dInteger), pointer :: nVerticesSolve !< Input: nVerticesSolve field information
     integer, intent(in) :: vertexDegree !< Input: vertexDegree dimension
     integer, intent(in) :: maxEdges !< Input: maxEdges dimension
     type (field1dInteger), pointer :: indexToCellID !< Input: indexToCellID field information
     type (field1dInteger), pointer :: indexToEdgeID !< Input: indexToEdgeID field information
     type (field1dInteger), pointer :: indexToVertexID !< Input: indexToVertexID field information

     type (domain_type), pointer :: domain

     type (field1dInteger), pointer :: indexToCellIDPoolField, indexToEdgeIDPoolField, indexToVertexIDPoolField
     integer, dimension(:), pointer :: indexToCellIDPool, indexToEdgeIDPool, indexToVertexIDPool

     type (block_type), pointer :: block_ptr
     type (field1dInteger), pointer :: nCellsCursor, nEdgesCursor, nVerticesCursor
     type (field1dInteger), pointer :: indexToCellCursor, indexToEdgeCursor, indexToVertexCursor

     integer :: nCellsSolve_0Halo, nVerticesSolve_0Halo, nEdgesSolve_0Halo
     integer :: blockID, localBlockID, err_level, iErr

     domain => blocklist % domain

     ! Loop over blocks
     block_ptr => blocklist
     nCellsCursor => nCellsSolve
     nEdgesCursor => nEdgesSolve
     nVerticesCursor => nVerticesSolve
     indexToCellCursor => indexToCellID
     indexToEdgeCursor => indexToEdgeID
     indexToVertexCursor => indexToVertexID
     do while(associated(block_ptr))
       ! Determine block dimensions
       nCells = nCellsCursor % array(nHalos+1)
       nEdges = nEdgesCursor % array(nHalos+2)
       nVertices = nVerticesCursor % array(nHalos+2)

       nCellsSolve_0Halo = nCellsCursor % array(1)
       nEdgesSolve_0Halo = nEdgesCursor % array(1)
       nVerticesSolve_0Halo = nVerticesCursor % array(1)

       ! Determine block IDs
       blockID = block_ptr % blockID
       localBlockID = block_ptr % localBlockID

       ! Allocate fields in block
       call mpas_allocate_block(nHalos, block_ptr, domain, blockID)

       ! Define block's static dimensions
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nCells', nCells)
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nEdges', nEdges)
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nVertices', nVertices)
       call mpas_pool_add_dimension(block_ptr % dimensions, 'vertexDegree', vertexDegree)
       call mpas_pool_add_dimension(block_ptr % dimensions, 'maxEdges', maxEdges)

       ! Set block's *Solve dimensions
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nCellsSolve', nCellsSolve_0Halo)
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nEdgesSolve', nEdgesSolve_0Halo)
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nVerticesSolve', nVerticesSolve_0Halo)

       call mpas_pool_add_dimension(block_ptr % dimensions, 'nCellsArray', nCellsCursor % array(:))
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nEdgesArray', nEdgesCursor % array(:))
       call mpas_pool_add_dimension(block_ptr % dimensions, 'nVerticesArray', nVerticesCursor % array(:))

       iErr = block_ptr % domain % core % setup_block(block_ptr)
       if ( iErr /= 0 ) then
          call mpas_dmpar_global_abort('ERROR: Block setup failed for core ' // trim(block_ptr % domain % core % coreName))
       end if

       ! Set block's local id
       block_ptr % localBlockID = localBlockID

       err_level = mpas_pool_get_error_level()
       call mpas_pool_set_error_level(MPAS_POOL_SILENT)

       ! Setup block's index fields
       allocate(indexToCellIDPoolField)
       indexToCellIDPoolField % block => block_ptr
       indexToCellIDPoolField % fieldName = 'indexToCellID'
       indexToCellIDPoolField % dimNames(1) = 'nCells'
       indexToCellIDPoolField % dimSizes(1) = nCells+1
       indexToCellIDPoolField % hasTimeDimension = .false.
       indexToCellIDPoolField % isActive = .true.
       indexToCellIDPoolField % isVarArray = .false.
       indexToCellIDPoolField % isPersistent = .true.
       nullify(indexToCellIDPoolField % next)
       nullify(indexToCellIDPoolField % prev)
       nullify(indexToCellIDPoolField % sendList)
       nullify(indexToCellIDPoolField % recvList)
       nullify(indexToCellIDPoolField % copyList)
       allocate(indexToCellIDPoolField % attLists(1))
       allocate(indexToCellIDPoolField % attLists(1) % attList)
       indexToCellIDPoolField % attLists(1) % attList % attName = ''
       indexToCellIDPoolField % attLists(1) % attList % attType = -1
       nullify(indexToCellIDPoolField % attLists(1) % attList % next)
       nullify(indexToCellIDPoolField % attLists(1) % attList % attValueIntA)
       nullify(indexToCellIDPoolField % attLists(1) % attList % attValueRealA)
       allocate(indexToCellIDPoolField % array(nCells+1))
       indexToCellIDPoolField % array(nCells+1) = 0

       call mpas_pool_add_field(block_ptr % allFields, 'indexToCellID_blk', indexToCellIDPoolField)

       allocate(indexToEdgeIDPoolField)
       indexToEdgeIDPoolField % block => block_ptr
       indexToEdgeIDPoolField % fieldName = 'indexToEdgeID'
       indexToEdgeIDPoolField % dimNames(1) = 'nEdges'
       indexToEdgeIDPoolField % dimSizes(1) = nEdges+1
       indexToEdgeIDPoolField % hasTimeDimension = .false.
       indexToEdgeIDPoolField % isActive = .true.
       indexToEdgeIDPoolField % isVarArray = .false.
       indexToEdgeIDPoolField % isPersistent = .true.
       nullify(indexToEdgeIDPoolField % next)
       nullify(indexToEdgeIDPoolField % prev)
       nullify(indexToEdgeIDPoolField % sendList)
       nullify(indexToEdgeIDPoolField % recvList)
       nullify(indexToEdgeIDPoolField % copyList)
       allocate(indexToEdgeIDPoolField % attLists(1))
       allocate(indexToEdgeIDPoolField % attLists(1) % attList)
       indexToEdgeIDPoolField % attLists(1) % attList % attName = ''
       indexToEdgeIDPoolField % attLists(1) % attList % attType = -1
       nullify(indexToEdgeIDPoolField % attLists(1) % attList % next)
       nullify(indexToEdgeIDPoolField % attLists(1) % attList % attValueIntA)
       nullify(indexToEdgeIDPoolField % attLists(1) % attList % attValueRealA)
       allocate(indexToEdgeIDPoolField % array(nEdges+1))
       indexToEdgeIDPoolField % array(nEdges+1) = 0

       call mpas_pool_add_field(block_ptr % allFields, 'indexToEdgeID_blk', indexToEdgeIDPoolField)

       allocate(indexToVertexIDPoolField)

       indexToVertexIDPoolField % block => block_ptr
       indexToVertexIDPoolField % fieldName = 'indexToVertexID'
       indexToVertexIDPoolField % dimNames(1) = 'nVertices'
       indexToVertexIDPoolField % dimSizes(1) = nVertices+1
       indexToVertexIDPoolField % hasTimeDimension = .false.
       indexToVertexIDPoolField % isActive = .true.
       indexToVertexIDPoolField % isVarArray = .false.
       indexToVertexIDPoolField % isPersistent = .true.
       nullify(indexToVertexIDPoolField % next)
       nullify(indexToVertexIDPoolField % prev)
       nullify(indexToVertexIDPoolField % sendList)
       nullify(indexToVertexIDPoolField % recvList)
       nullify(indexToVertexIDPoolField % copyList)
       allocate(indexToVertexIDPoolField % attLists(1))
       allocate(indexToVertexIDPoolField % attLists(1) % attList)
       indexToVertexIDPoolField % attLists(1) % attList % attName = ''
       indexToVertexIDPoolField % attLists(1) % attList % attType = -1
       nullify(indexToVertexIDPoolField % attLists(1) % attList % next)
       nullify(indexToVertexIDPoolField % attLists(1) % attList % attValueIntA)
       nullify(indexToVertexIDPoolField % attLists(1) % attList % attValueRealA)
       allocate(indexToVertexIDPoolField % array(nVertices+1))
       indexToVertexIDPoolField % array(nVertices+1) = 0

       call mpas_pool_add_field(block_ptr % allFields, 'indexToVertexID_blk', indexToVertexIDPoolField)

       ! Set block's 0 halo indices
       call mpas_pool_get_array(block_ptr % allFields, 'indexToCellID_blk', indexToCellIDPool)
       call mpas_pool_get_array(block_ptr % allFields, 'indexToEdgeID_blk', indexToEdgeIDPool)
       call mpas_pool_get_array(block_ptr % allFields, 'indexToVertexID_blk', indexToVertexIDPool)
       indexToCellIDPool(1:nCellsSolve_0Halo) = indexToCellCursor % array(1:nCellsSolve_0Halo)
       indexToEdgeIDPool(1:nEdgesSolve_0Halo) = indexToEdgeCursor % array(1:nEdgesSolve_0Halo)
       indexToVertexIDPool(1:nVerticesSolve_0Halo) = indexToVertexCursor % array(1:nVerticesSolve_0Halo)
       call mpas_pool_set_error_level(err_level)

       ! Set block's exchange lists and nullify unneeded exchange lists
       block_ptr % parinfo % cellsToSend => indexToCellCursor % sendList
       block_ptr % parinfo % cellsToRecv => indexToCellCursor % recvList
       block_ptr % parinfo % cellsToCopy => indexToCellCursor % copyList
       nullify(indexToCellCursor % sendList)
       nullify(indexToCellCursor % recvList)
       nullify(indexToCellCursor % copyList)

       block_ptr % parinfo % edgesToSend => indexToEdgeCursor % sendList
       block_ptr % parinfo % edgesToRecv => indexToEdgeCursor % recvList
       block_ptr % parinfo % edgesToCopy => indexToEdgeCursor % copyList
       nullify(indexToEdgeCursor % sendList)
       nullify(indexToEdgeCursor % recvList)
       nullify(indexToEdgeCursor % copyList)

       block_ptr % parinfo % verticesToSend => indexToVertexCursor % sendList
       block_ptr % parinfo % verticesToRecv => indexToVertexCursor % recvList
       block_ptr % parinfo % verticesToCopy => indexToVertexCursor % copyList
       nullify(indexToVertexCursor % sendList)
       nullify(indexToVertexCursor % recvList)
       nullify(indexToVertexCursor % copyList)

       ! Setup next/prev multihalo exchange list pointers
       !   (block 'next' pointers should be setup by now, but setting up 'next' pointers indirectly here)
       if ( associated(block_ptr % prev) ) then
          ! == Setup this block's 'prev' pointers ==
          ! 1. For Cell exchange lists
          block_ptr % parinfo % cellsToSend % prev => block_ptr % prev % parinfo % cellsToSend
          block_ptr % parinfo % cellsToRecv % prev => block_ptr % prev % parinfo % cellsToRecv
          block_ptr % parinfo % cellsToCopy % prev => block_ptr % prev % parinfo % cellsToCopy
          ! 2. For Edge exchange lists
          block_ptr % parinfo % edgesToSend % prev => block_ptr % prev % parinfo % edgesToSend
          block_ptr % parinfo % edgesToRecv % prev => block_ptr % prev % parinfo % edgesToRecv
          block_ptr % parinfo % edgesToCopy % prev => block_ptr % prev % parinfo % edgesToCopy
          ! 3. For Vertex exchange lists
          block_ptr % parinfo % verticesToSend % prev => block_ptr % prev % parinfo % verticesToSend
          block_ptr % parinfo % verticesToRecv % prev => block_ptr % prev % parinfo % verticesToRecv
          block_ptr % parinfo % verticesToCopy % prev => block_ptr % prev % parinfo % verticesToCopy
          ! == Setup the previous block's 'next' pointers ==
          ! 1. For Cell exchange lists
          block_ptr % prev % parinfo % cellsToSend % next => block_ptr % parinfo % cellsToSend
          block_ptr % prev % parinfo % cellsToRecv % next => block_ptr % parinfo % cellsToRecv
          block_ptr % prev % parinfo % cellsToCopy % next => block_ptr % parinfo % cellsToCopy
          ! 2. For Edge exchange lists
          block_ptr % prev % parinfo % edgesToSend % next => block_ptr % parinfo % edgesToSend
          block_ptr % prev % parinfo % edgesToRecv % next => block_ptr % parinfo % edgesToRecv
          block_ptr % prev % parinfo % edgesToCopy % next => block_ptr % parinfo % edgesToCopy
          ! 3. For Vertex exchange lists
          block_ptr % prev % parinfo % verticesToSend % next => block_ptr % parinfo % verticesToSend
          block_ptr % prev % parinfo % verticesToRecv % next => block_ptr % parinfo % verticesToRecv
          block_ptr % prev % parinfo % verticesToCopy % next => block_ptr % parinfo % verticesToCopy
          ! (the final block's 'next' pointer does not need to be dealt with because it was alredy nullified in mpas_dmpar_init_multihalo_exchange_list)
       end if

       ! Advance cursors
       block_ptr => block_ptr % next
       nCellsCursor => nCellsCursor % next
       nEdgesCursor => nEdgesCursor % next
       nVerticesCursor => nVerticesCursor % next
       indexToCellCursor => indexToCellCursor % next
       indexToEdgeCursor => indexToEdgeCursor % next
       indexToVertexCursor => indextoVertexcursor % next
     end do

   end subroutine mpas_block_creator_finalize_block_phase1!}}}

!***********************************************************************
!
!  routine mpas_block_creator_finalize_block_phase2
!
!> \brief   Phase 2 of block creation finalization
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details
!>  This routine finalizes the block setup process. It sets up all additional
!>  dimensions, and allocates space for fields. In addition, it groups all fields
!>  into their respective structures and puts them into the correct blocks.
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_finalize_block_phase2(stream_manager, blocklist, readableDimensions) !{{{

     use mpas_stream_manager
     use mpas_block_decomp

     type (mpas_streamManager_type), pointer :: stream_manager !< Input: Stream manager structure
     type (block_type), pointer :: blocklist !< Input/Output: Linked List of blocks
     type (mpas_pool_type), intent(inout) :: readableDimensions

     type (domain_type), pointer :: domain

     type (block_type), pointer :: block_ptr

     type (mpas_pool_iterator_type) :: itr

     type (field1dInteger), pointer :: indexField, indexFieldBlk, ownedIndices

     integer :: err_level, iErr
     integer, pointer :: dim0d
     integer, dimension(:), pointer :: dim1d

     domain => blocklist % domain

     ! Loop over blocks
     block_ptr => blocklist
     do while(associated(block_ptr))

       err_level = mpas_pool_get_error_level()
       call mpas_pool_set_error_level(MPAS_POOL_SILENT)
       call mpas_pool_begin_iteration(readableDimensions)
       do while (mpas_pool_get_next_member(readableDimensions, itr))
          if ( itr % memberType == MPAS_POOL_DIMENSION ) then
             if ( itr % nDims == 0 ) then
                call mpas_pool_get_dimension(block_ptr % dimensions, trim(itr % memberName), dim0d)
                if ( .not. associated(dim0d) ) then
                    call mpas_pool_get_dimension(readableDimensions, itr % memberName, dim0d)
                    call mpas_pool_add_dimension(block_ptr % dimensions, trim(itr % memberName), dim0d)
                end if
             else if ( itr % nDims == 1 ) then
                call mpas_pool_get_dimension(block_ptr % dimensions, trim(itr % memberName), dim1d)
                if ( .not. associated(dim1d) ) then
                    call mpas_pool_get_dimension(readableDimensions, itr % memberName, dim1d)
                    call mpas_pool_add_dimension(block_ptr % dimensions, trim(itr % memberName), dim1d)
                end if
             end if
          end if
       end do
       call mpas_pool_set_error_level(err_level)

       iErr = block_ptr % domain % core % setup_decomposed_dimensions(block_ptr, stream_manager, readableDimensions, block_ptr % dimensions, domain % dminfo % total_blocks)
       if ( iErr /= 0 ) then
          call mpas_dmpar_global_abort('ERROR: Decomposed dimension setup failed for core ' // trim(block_ptr % domain % core % coreName))
       end if

       iErr = block_ptr % domain % core % setup_derived_dimensions(readableDimensions, block_ptr % dimensions, block_ptr % configs)
       if ( iErr /= 0 ) then
          call mpas_dmpar_global_abort('ERROR: Derived dimension setup failed for core ' // trim(block_ptr % domain % core % coreName))
       end if

       call mpas_block_creator_allocate_pool_fields(block_ptr % structs, block_ptr % dimensions)

       err_level = mpas_pool_get_error_level()
       call mpas_pool_set_error_level(MPAS_POOL_SILENT)
       call mpas_pool_get_field(block_ptr % allFields, 'indexToCellID', indexField)
       call mpas_pool_get_field(block_ptr % allFields, 'indexToCellID_blk', indexFieldBlk)
       call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', dim0d)
       allocate(ownedIndices)
       ownedIndices % fieldName = 'nCellsOwnedIndices'
       ownedIndices % dimNames(1) = 'nCells'
       ownedIndices % dimSizes(1) = dim0d
       allocate(ownedIndices % array(dim0d))
       ownedIndices % isDecomposed = .false.
       ownedIndices % hasTimeDimension = .false.
       ownedIndices % isActive = .true.
       ownedIndices % isVarArray = .false.
       ownedIndices % isPersistent = .true.
       allocate(ownedIndices % attLists(1))
       allocate(ownedIndices % attLists(1) % attList)
       ownedIndices % attLists(1) % attList % attName = ''
       ownedIndices % attLists(1) % attList % attType = -1
       nullify(ownedIndices % attLists(1) % attList % next)
       nullify(ownedIndices % attLists(1) % attList % attValueIntA)
       nullify(ownedIndices % attLists(1) % attList % attValueRealA)
       call mpas_pool_add_field(block_ptr % allFields, 'nCellsOwnedIndices', ownedIndices)

       indexField % array(:) = indexFieldBlk % array(:)

       call mpas_pool_remove_field(block_ptr % allFields, 'indextoCellID_blk')
       call mpas_deallocate_field(indexFieldBlk)

       call mpas_pool_get_field(block_ptr % allFields, 'indexToEdgeID', indexField)
       call mpas_pool_get_field(block_ptr % allFields, 'indexToEdgeID_blk', indexFieldBlk)
       call mpas_pool_get_dimension(block_ptr % dimensions, 'nEdgesSolve', dim0d)
       allocate(ownedIndices)
       ownedIndices % fieldName = 'nEdgesOwnedIndices'
       ownedIndices % dimNames(1) = 'nEdges'
       ownedIndices % dimSizes(1) = dim0d
       allocate(ownedIndices % array(dim0d))
       ownedIndices % isDecomposed = .false.
       ownedIndices % hasTimeDimension = .false.
       ownedIndices % isActive = .true.
       ownedIndices % isVarArray = .false.
       ownedIndices % isPersistent = .true.
       allocate(ownedIndices % attLists(1))
       allocate(ownedIndices % attLists(1) % attList)
       ownedIndices % attLists(1) % attList % attName = ''
       ownedIndices % attLists(1) % attList % attType = -1
       nullify(ownedIndices % attLists(1) % attList % next)
       nullify(ownedIndices % attLists(1) % attList % attValueIntA)
       nullify(ownedIndices % attLists(1) % attList % attValueRealA)
       call mpas_pool_add_field(block_ptr % allFields, 'nEdgesOwnedIndices', ownedIndices)

       indexField % array(:) = indexFieldBlk % array(:)

       call mpas_pool_remove_field(block_ptr % allFields, 'indextoEdgeID_blk')
       call mpas_deallocate_field(indexFieldBlk)

       call mpas_pool_get_field(block_ptr % allFields, 'indexToVertexID', indexField)
       call mpas_pool_get_field(block_ptr % allFields, 'indexToVertexID_blk', indexFieldBlk)
       call mpas_pool_get_dimension(block_ptr % dimensions, 'nVerticesSolve', dim0d)
       allocate(ownedIndices)
       ownedIndices % fieldName = 'nVerticesOwnedIndices'
       ownedIndices % dimNames(1) = 'nVertices'
       ownedIndices % dimSizes(1) = dim0d
       allocate(ownedIndices % array(dim0d))
       ownedIndices % isDecomposed = .false.
       ownedIndices % hasTimeDimension = .false.
       ownedIndices % isActive = .true.
       ownedIndices % isVarArray = .false.
       ownedIndices % isPersistent = .true.
       allocate(ownedIndices % attLists(1))
       allocate(ownedIndices % attLists(1) % attList)
       ownedIndices % attLists(1) % attList % attName = ''
       ownedIndices % attLists(1) % attList % attType = -1
       nullify(ownedIndices % attLists(1) % attList % next)
       nullify(ownedIndices % attLists(1) % attList % attValueIntA)
       nullify(ownedIndices % attLists(1) % attList % attValueRealA)
       call mpas_pool_add_field(block_ptr % allFields, 'nVerticesOwnedIndices', ownedIndices)

       indexField % array(:) = indexFieldBlk % array(:)

       call mpas_pool_remove_field(block_ptr % allFields, 'indextoVertexID_blk')
       call mpas_deallocate_field(indexFieldBlk)
       call mpas_pool_set_error_level(err_level)

       block_ptr => block_ptr % next
     end do

   end subroutine mpas_block_creator_finalize_block_phase2!}}}


!***********************************************************************
!
!  routine mpas_block_creator_reindex_block_fields
!
!> \brief   Reindex mesh connectivity arrays
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine re-indexes the connectivity arrays for the mesh data
!>  structure. Prior to this routine, all indices are given as global index (which
!>  can later be found in the indexTo* arrays). After this routine is called,
!>  indices are provided as local indices now (1:nCells+1 ... etc).
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_reindex_block_fields(blocklist)!{{{
     type (block_type), pointer :: blocklist !< Input/Output: Linked list of blocks

     type (block_type), pointer :: block_ptr

     integer :: i, j, k
     integer, dimension(:,:), pointer :: cellIDSorted, edgeIDSorted, vertexIDSorted

     integer, pointer :: nCells, nEdges, nVertices, vertexDegree
     integer, dimension(:), pointer :: indexToCellID, indexToEdgeID, indexToVertexID
     integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
     integer, dimension(:,:), pointer :: cellsOnCell, cellsOnVertex, cellsOnEdge
     integer, dimension(:,:), pointer :: edgesOnCell, edgesOnVertex, edgesOnEdge
     integer, dimension(:,:), pointer :: verticesOnCell, verticesOnEdge

     ! Loop over blocks
     block_ptr => blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
       call mpas_pool_get_dimension(block_ptr % dimensions, 'nEdges', nEdges)
       call mpas_pool_get_dimension(block_ptr % dimensions, 'nVertices', nVertices)
       call mpas_pool_get_dimension(block_ptr % dimensions, 'vertexDegree', vertexDegree)

       call mpas_pool_get_array(block_ptr % allFields, 'indexToCellID', indexToCellID)
       call mpas_pool_get_array(block_ptr % allFields, 'indexToEdgeID', indexToEdgeID)
       call mpas_pool_get_array(block_ptr % allFields, 'indexToVertexID', indexToVertexID)
       call mpas_pool_get_array(block_ptr % allFields, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(block_ptr % allFields, 'nEdgesOnEdge', nEdgesOnEdge)
       call mpas_pool_get_array(block_ptr % allFields, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_array(block_ptr % allFields, 'cellsOnVertex', cellsOnVertex)
       call mpas_pool_get_array(block_ptr % allFields, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_array(block_ptr % allFields, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_array(block_ptr % allFields, 'edgesOnVertex', edgesOnVertex)
       call mpas_pool_get_array(block_ptr % allFields, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(block_ptr % allFields, 'verticesOnEdge', verticesOnEdge)
       call mpas_pool_get_array(block_ptr % allFields, 'edgesOnEdge', edgesOnEdge)

       !
       ! Rename vertices in cellsOnCell, edgesOnCell, etc. to local indices
       !
       allocate(cellIDSorted(2, nCells))
       allocate(edgeIDSorted(2, nEdges))
       allocate(vertexIDSorted(2, nVertices))
 
       do i = 1, nCells
         cellIDSorted(1,i) = indexToCellID(i)
         cellIDSorted(2,i) = i
       end do
       call mpas_quicksort(nCells, cellIDSorted)
 
       do i = 1, nEdges
         edgeIDSorted(1,i) = indexToEdgeID(i)
         edgeIDSorted(2,i) = i
       end do
       call mpas_quicksort(nEdges, edgeIDSorted)
 
       do i = 1, nVertices
         vertexIDSorted(1,i) = indexToVertexID(i)
         vertexIDSorted(2,i) = i
       end do
       call mpas_quicksort(nVertices, vertexIDSorted)
 
 
       do i = 1, nCells
         do j = 1, nEdgesOnCell(i)
           k = mpas_binary_search(cellIDSorted, 2, 1, nCells, cellsOnCell(j,i))
           if (k <= nCells) then
             cellsOnCell(j,i) = cellIDSorted(2,k)
           else
             cellsOnCell(j,i) = nCells + 1
           end if
 
           k = mpas_binary_search(edgeIDSorted, 2, 1, nEdges, edgesOnCell(j,i))
           if (k <= nEdges) then
             edgesOnCell(j,i) = edgeIDSorted(2,k)
           else
             edgesOnCell(j,i) = nEdges + 1
           end if
  
           k = mpas_binary_search(vertexIDSorted, 2, 1, nVertices, verticesOnCell(j,i))
           if (k <= nVertices) then
             verticesOnCell(j,i) = vertexIDSorted(2,k)
           else
             verticesOnCell(j,i) = nVertices + 1
           end if
         end do
       end do
  
       do i = 1, nEdges
         do j = 1, 2
  
           k = mpas_binary_search(cellIDSorted, 2, 1, nCells, cellsOnEdge(j,i))
           if (k <= nCells) then
             cellsOnEdge(j,i) = cellIDSorted(2,k)
           else
             cellsOnEdge(j,i) = nCells + 1
           end if
  
           k = mpas_binary_search(vertexIDSorted, 2, 1, nVertices, verticesOnEdge(j,i))
           if (k <= nVertices) then
             verticesOnEdge(j,i) = vertexIDSorted(2,k)
           else
             verticesOnEdge(j,i) = nVertices + 1
           end if
  
         end do
  
         do j = 1, nEdgesOnEdge(i)
  
           k = mpas_binary_search(edgeIDSorted, 2, 1, nEdges, edgesOnEdge(j,i))
           if (k <= nEdges) then
             edgesOnEdge(j,i) = edgeIDSorted(2,k)
           else
             edgesOnEdge(j,i) = nEdges + 1
           end if
         end do
       end do
  
       do i = 1, nVertices
         do j = 1, vertexDegree
  
           k = mpas_binary_search(cellIDSorted, 2, 1, nCells, cellsOnVertex(j,i))
           if (k <= nCells) then
             cellsOnVertex(j,i) = cellIDSorted(2,k)
           else
             cellsOnVertex(j,i) = nCells + 1
           end if
  
           k = mpas_binary_search(edgeIDSorted, 2, 1, nEdges, edgesOnVertex(j,i))
           if (k <= nEdges) then
             edgesOnVertex(j,i) = edgeIDSorted(2,k)
           else
             edgesOnVertex(j,i) = nEdges + 1
           end if
         end do
       end do
  
       deallocate(cellIDSorted)
       deallocate(edgeIDSorted)
       deallocate(vertexIDSorted)

       block_ptr => block_ptr % next
     end do

   end subroutine mpas_block_creator_reindex_block_fields!}}}

!***********************************************************************
!
!  routine mpas_block_creator_allocate_pool_fields
!
!> \brief   Pool field allocator
!> \author  Doug Jacobsen
!> \date    02/05/2015
!> \details
!>  This routine iterates over all pools, and allocates the fields in each pool
!>  based on their pre-defined dimension information.
!>  On each field, the dimNames array is used to query the correct dimension information
!>  and populate the dimSizes array.
!>  This routine also copies all dimensions from dimensionPool to currentPool
!
!-----------------------------------------------------------------------
   recursive subroutine mpas_block_creator_allocate_pool_fields(currentPool, dimensionPool)!{{{
      type (mpas_pool_type), pointer :: currentPool !< Input: Current pool to allocate and copy dimensions.
      type (mpas_pool_type), pointer :: dimensionPool !< Input: Pool of dimensions for the current block
      type (mpas_pool_type), pointer :: subPool
      type (mpas_pool_iterator_type) :: poolItr

      type (field1DReal), pointer :: real1DField
      type (field2DReal), pointer :: real2DField
      type (field3DReal), pointer :: real3DField
      type (field4DReal), pointer :: real4DField
      type (field5DReal), pointer :: real5DField
      type (field1DInteger), pointer :: int1DField
      type (field2DInteger), pointer :: int2DField
      type (field3DInteger), pointer :: int3DField
      type (field1DChar), pointer :: char1DField

      integer :: timeLev, iDim
      integer, pointer :: tempDim
      integer, dimension(:), pointer :: tempDim1D
      integer :: dimSize
      integer :: localErr

      call mpas_pool_begin_iteration(dimensionPool)
      do while( mpas_pool_get_next_member(dimensionPool, poolItr) )
         if (poolItr % memberType == MPAS_POOL_DIMENSION) then
            if (poolItr % nDims == 0) then
               call mpas_pool_get_dimension(dimensionPool, poolItr % memberName, tempDim)
               call mpas_pool_add_dimension(currentPool, poolItr % memberName, tempDim)
            else if (poolItr % nDims == 1) then
               call mpas_pool_get_dimension(dimensionPool, poolItr % memberName, tempDim1D)
               call mpas_pool_add_dimension(currentPool, poolItr % memberName, tempDim1D)
            end if
         end if
      end do

      call mpas_pool_begin_iteration(currentPool)
      do while ( mpas_pool_get_next_member(currentPool, poolItr) )
         if ( poolItr % memberType == MPAS_POOL_SUBPOOL ) then
            call mpas_pool_get_subpool(currentPool, poolItr % memberName, subPool)
            call mpas_block_creator_allocate_pool_fields(subPool, dimensionPool)
         else if ( poolItr % memberType == MPAS_POOL_FIELD ) then
            if ( poolItr % dataType == MPAS_POOL_REAL ) then
               if ( poolItr % nDims == 1 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, real1DField, timeLev)

                     if ( real1DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, real1DField % dimNames(iDim), tempDim)

                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(real1DField % dimNames(iDim), poolItr % memberName)
                           end if
                           real1DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                          real1DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( real1DField % isPersistent ) then
                           allocate(real1DField % array(real1DField % dimSizes(1)))
                           real1DField % array(:) = real1DField % defaultValue
                        end if
                     end if
                  end do
               else if ( poolItr % nDims == 2 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, real2DField, timeLev)

                     if ( real2DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, real2DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(real2DField % dimNames(iDim), poolItr % memberName)
                           end if
                           real2DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                          real2DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( real2DField % isPersistent ) then
                           allocate(real2DField % array(real2DField % dimSizes(1), real2DField % dimSizes(2)))
                           real2DField % array(:,:) = real2DField % defaultValue
                        end if
                     end if
                  end do
               else if ( poolItr % nDims == 3 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, real3DField, timeLev)

                     if ( real3DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, real3DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(real3DField % dimNames(iDim), poolItr % memberName)
                           end if
                           real3DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                          real3DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( real3DField % isPersistent ) then
                           allocate(real3DField % array(real3DField % dimSizes(1), real3DField % dimSizes(2), real3DField % dimSizes(3)))
                           real3DField % array(:,:,:) = real3DField % defaultValue
                        end if
                     end if
                  end do
               else if ( poolItr % nDims == 4 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, real4DField, timeLev)

                     if ( real4DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, real4DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(real4DField % dimNames(iDim), poolItr % memberName)
                           end if
                           real4DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                          real4DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( real4DField % isPersistent ) then
                           allocate(real4DField % array(real4DField % dimSizes(1), real4DField % dimSizes(2), &
                                                        real4DField % dimSizes(3), real4DField % dimSizes(4)))
                           real4DField % array(:,:,:,:) = real4DField % defaultValue
                        end if
                     end if
                  end do
               else if ( poolItr % nDims == 5 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, real5DField, timeLev)

                     if ( real5DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, real5DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(real5DField % dimNames(iDim), poolItr % memberName)
                           end if
                           real5DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                          real5DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( real5DField % isPersistent ) then
                           allocate(real5DField % array(real5DField % dimSizes(1), real5DField % dimSizes(2), &
                                                        real5DField % dimSizes(3), real5DField % dimSizes(4), &
                                                        real5DField % dimSizes(5)))
                           real5DField % array(:,:,:,:,:) = real5DField % defaultValue
                        end if
                     end if
                  end do
               end if
            else if ( poolItr % dataType == MPAS_POOL_INTEGER ) then
               if ( poolItr % nDims == 1 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, int1DField, timeLev)

                     if ( int1DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, int1DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(int1DField % dimNames(iDim), poolItr % memberName)
                           end if
                           int1DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                         int1DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( int1DField % isPersistent ) then
                           allocate(int1DField % array(int1DField % dimSizes(1)))
                           int1DField % array(:) = int1DField % defaultValue
                        end if
                     end if
                  end do
               else if ( poolItr % nDims == 2 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, int2DField, timeLev)

                     if ( int2DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, int2DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(int2DField % dimNames(iDim), poolItr % memberName)
                           end if
                           int2DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                         int2DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( int2DField % isPersistent ) then
                           allocate(int2DField % array(int2DField % dimSizes(1), int2DField % dimSizes(2)))
                           int2DField % array(:,:) = int2DField % defaultValue
                        end if
                     end if
                  end do
               else if ( poolItr % nDims == 3 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, int3DField, timeLev)

                     if ( int3DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, int3DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(int3DField % dimNames(iDim), poolItr % memberName)
                           end if
                           int3DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                         int3DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( int3DField % isPersistent ) then
                           allocate(int3DField % array(int3DField % dimSizes(1), int3DField % dimSizes(2), int3DField % dimSizes(3)))
                           int3DField % array(:,:,:) = int3DField % defaultValue
                        end if
                     end if
                  end do
               end if
            else if ( poolItr % dataType == MPAS_POOL_CHARACTER ) then
               if ( poolItr % nDims == 1 ) then
                  do timeLev = 1, poolItr % nTimeLevels
                     call mpas_pool_get_field(currentPool, poolItr % memberName, char1DField, timeLev)

                     if ( char1DField % isActive ) then
                        do iDim = 1, poolItr % nDims
                           call mpas_pool_get_dimension(currentPool, char1DField % dimNames(iDim), tempDim)
                           if ( .not. associated(tempDim) .or. tempDim == MPAS_MISSING_DIM ) then
                               call missing_dim_abort(char1DField % dimNames(iDim), poolItr % memberName)
                           end if
                           char1DField % dimSizes(iDim) = tempDim + mpas_dimension_num_garbage_elements( &
                                                          char1DField % dimNames(iDim), iErr = localErr)
                        end do

                        if ( char1DField % isPersistent ) then
                           allocate(char1DField % array(char1DField % dimSizes(1)))
                           char1DField % array(:) = char1DField % defaultValue
                        end if
                     end if
                  end do
               end if
            end if
         end if
      end do

   end subroutine mpas_block_creator_allocate_pool_fields!}}}


   subroutine missing_dim_abort(dimName, fieldName)

      implicit none

      character(len=*), intent(in) :: dimName
      character(len=*), intent(in) :: fieldName


      write(stderrUnit,*) '********************************************************************************'
      write(stderrUnit,*) 'ERROR: Dimension '''//trim(dimName)//''' required by field '''//trim(fieldName)//'''' &
                          //' was neither read nor defined.'
      call mpas_dmpar_global_abort('********************************************************************************')

   end subroutine missing_dim_abort 

end module mpas_block_creator
