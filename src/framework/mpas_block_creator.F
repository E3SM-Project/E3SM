! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
   use mpas_dmpar_types
   use mpas_block_decomp
   use mpas_hash
   use mpas_sort
   use mpas_grid_types
   use mpas_configure

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

   subroutine mpas_block_creator_setup_blocks_and_0halo_cells(domain, indexToCellID, cellList, blockID, blockStart, blockCount)!{{{
     type (domain_type), pointer :: domain !< Input: Domain information
     type (field1dInteger), pointer :: indexToCellID !< Input/Output: indexToCellID field
     integer, dimension(:), intent(in) :: cellList !< Input: List of cell indices owned by this processor
     integer, dimension(:), intent(in) :: blockID !< Input: List of block indices owned by this processor
     integer, dimension(:), intent(in) :: blockStart !< Input: Indices of starting cell id in cellList for each block
     integer, dimension(:), intent(in) :: blockCount !< Input: Number of cells from cellList owned by each block.
 
     integer :: nHalos
     type (block_type), pointer :: blockCursor
     type (field1dInteger), pointer :: fieldCursor
 
     integer :: i, iHalo
     integer :: nBlocks
 
     nBlocks = size(blockID)
     nHalos = config_num_halos

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
       nullify(fieldCursor % ioinfo)
 
       ! Initialize exchange lists
       call mpas_dmpar_init_mulithalo_exchange_list(fieldCursor % sendList, nHalos)
       call mpas_dmpar_init_mulithalo_exchange_list(fieldCursor % recvList, nHalos)
       call mpas_dmpar_init_mulithalo_exchange_list(fieldCursor % copyList, nHalos)
 
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

   subroutine mpas_block_creator_build_0halo_cell_fields(indexToCellIDBlock, nEdgesOnCellBlock, cellsOnCellBlock, verticesOnCellBlock, edgesOnCellBlock, indexToCellID_0Halo, nEdgesOnCell_0Halo, cellsOnCell_0Halo, verticesOnCell_0Halo, edgesOnCell_0Halo)!{{{
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

     type (mpas_exchange_list), pointer :: exchListPtr

     integer :: nCellsInBlock, maxEdges, nHalos
     integer :: i, iHalo

     nHalos = config_num_halos

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

       ! Nullify ioinfo, since this data is not read in
       nullify(nEdgesCursor % ioinfo)
       nullify(cellsOnCellCursor % ioinfo)
       nullify(verticesOnCellCursor % ioinfo)
       nullify(edgesOnCellCursor % ioinfo)

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
!>  NOTE: This routine can be used on either edges or edges
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_build_0_and_1halo_edge_fields(indexToEdgeIDBlock, cellsOnEdgeBlock, indexToCellID_0Halo, nEdgesOnCell_0Halo, edgesOnCell_0Halo, indexToEdgeID_0Halo, cellsOnEdge_0Halo, nEdgesSolve)!{{{
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
     type (field2dInteger), pointer :: edgesOnCellCursor, cellsOnEdgeCursor, cellsOnCellCursor

     type (mpas_exchange_list), pointer :: exchListPtr

     integer, dimension(:), pointer :: localEdgeList
     integer, dimension(:), pointer :: sendingHaloLayers
     integer :: nEdgesLocal, nCellsInBlock, maxEdges, edgeDegree, nHalos
     integer :: haloStart
     integer :: iBlock, i, j, k

     ! Setup sendingHaloLayers
     allocate(sendingHaloLayers(1))
     sendingHaloLayers(1) = 1

     ! Get dimension information
     maxEdges = edgesOnCell_0Halo % dimSizes(1)
     edgeDegree = cellsOnEdgeBlock % dimSizes(1)
     nHalos = config_num_halos

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
       nullify(indexToEdgeCursor % ioinfo)
       indexToEdgeCursor % dimSizes(1) = nEdgesLocal
       allocate(indexToEdgeCursor % array(indexToEdgeCursor % dimSizes(1)))
       indexToEdgeCursor % array(:) = localEdgeList(:)

       ! Setup cellsOnEdge block
       cellsOnEdgeCursor % block => indexToCellCursor % block
       nullify(cellsOnEdgeCursor % ioinfo)
       cellsOnEdgeCursor % dimSizes(1) = edgeDegree
       cellsOnEdgeCursor % dimSizes(2) = nEdgesLocal
       allocate(cellsOnEdgeCursor % array(cellsOnEdgeCursor % dimSizes(1), cellsOnEdgeCursor % dimSizes(2)))

       ! Setup exchange lists
       call mpas_dmpar_init_mulithalo_exchange_list(indexToEdgeCursor % sendList, nHalos+1)
       call mpas_dmpar_init_mulithalo_exchange_list(indexToEdgeCursor % recvList, nHalos+1)
       call mpas_dmpar_init_mulithalo_exchange_list(indexToEdgeCursor % copyList, nHalos+1)

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

       ! Nullify io info
       nullify(haloCursor % ioinfo)
       nullify(offSetCursor % ioinfo)
       nullify(edgeLimitCursor % ioinfo)
       nullify(nEdgesSolveCursor % ioinfo)

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

   subroutine mpas_block_creator_build_cell_halos(indexToCellID, nEdgesOnCell, cellsOnCell, verticesOnCell, edgesOnCell, nCellsSolve)!{{{
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

     type (mpas_exchange_list), pointer :: exchListPtr

     integer :: nHalos, nCellsInBlock, nCellsInHalo, maxEdges
     integer :: iHalo, iBlock, i

     nHalos = config_num_halos
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
       nullify(offSetCursor % ioinfo)

       ! Setup nCellsSolve
       nCellsSolveCursor % dimSizes(1) = nHalos+1
       allocate(nCellsSolveCursor % array(nCellsSolveCursor % dimSizes(1)))
       nCellsSolveCursor % array(1) = indexCursor % dimSizes(1)
       nCellsSolveCursor % block => indexCursor % block
       nullify(nCellsSolveCursor % ioinfo)

       ! Setup owned cellLimit
       cellLimitCursor % scalar = indexCursor % dimSizes(1)
       cellLimitCursor % block => indexCursor % block
       nullify(cellLimitCursor % ioinfo)

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
         nullify(haloCursor % ioinfo)

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

   subroutine mpas_block_creator_build_edge_halos(indexToCellID, nEdgesOnCell, nCellsSolve, edgesOnCell, indexToEdgeID, cellsOnEdge, nEdgesSolve)!{{{
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
     type (field1dInteger), pointer :: indexToCellCursor, nEdgesCursor, nCellsSolveCursor, indexToEdgeCursor, nEdgesSolveCursor, haloCursor
     type (field2dInteger), pointer :: edgesOnCellCursor, cellsOnEdgeCursor

     integer, dimension(:), pointer :: sendingHaloLayers
     integer, dimension(:), pointer :: array1dHolder, localEdgeList
     integer, dimension(:,:), pointer :: array2dHolder

     integer :: iHalo, iBlock, i, j, k
     integer :: nHalos, nBlocks, nCellsInBlock, nEdgesLocal, haloStart, haloEnd, haloSize
     integer :: maxEdges, edgeDegree

     type (hashtable), dimension(:), pointer :: edgeList

     ! Determine dimensions
     nHalos = config_num_halos
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

       ! Nullify ioinfo
       nullify(edgeLimitCursor % ioinfo)
       nullify(offSetCursor % ioinfo)
       nullify(haloCursor % ioinfo)

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
!  routine mpas_block_creator_finalize_block_init
!
!> \brief   Finalize block creation
!> \author  Doug Jacobsen
!> \date    05/31/12
!> \details 
!>  This routine finalizes the block initialization processor. It calls
!>  mpas_block_allocate to allocate space for all fields in a block. Then the 0
!>  halo indices for each element and the exchange lists are copied into the
!>  appropriate block. A halo update is required after this routien is called
!>  to make sure all data in a block is valid.
!
!-----------------------------------------------------------------------

   subroutine mpas_block_creator_finalize_block_init(blocklist, &  !{{{
#include "dim_dummy_args.inc"
                                                     , nCellsSolve, nEdgesSolve, nVerticesSolve, indexToCellID, indexToEdgeID, indexToVertexID)
     type (block_type), pointer :: blocklist !< Input/Output: Linked List of blocks
#include "dim_dummy_decls_inout.inc"
     type (field1dInteger), pointer :: nCellsSolve !< Input: nCellsSolve field information
     type (field1dInteger), pointer :: nEdgesSolve !< Input: nEdgesSolve field information
     type (field1dInteger), pointer :: nVerticesSolve !< Input: nVerticesSolve field information
     type (field1dInteger), pointer :: indexToCellID !< Input: indexToCellID field information
     type (field1dInteger), pointer :: indexToEdgeID !< Input: indexToEdgeID field information
     type (field1dInteger), pointer :: indexToVertexID !< Input: indexToVertexID field information

     type (domain_type), pointer :: domain

     type (block_type), pointer :: block_ptr
     type (field1dInteger), pointer :: nCellsCursor, nEdgesCursor, nVerticesCursor
     type (field1dInteger), pointer :: indexToCellCursor, indexToEdgeCursor, indexToVertexCursor

     integer :: nHalos
     integer :: nCellsSolve_0Halo, nVerticesSolve_0Halo, nEdgesSolve_0Halo
     integer :: blockID, localBlockID

     nHalos = config_num_halos
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
       call mpas_allocate_block(nHalos, block_ptr, domain, blockID, &
#include "dim_dummy_args.inc"
                               )

       allocate(block_ptr % mesh % nCellsArray(0:nHalos))
       allocate(block_ptr % mesh % nEdgesArray(0:nHalos+1))
       allocate(block_ptr % mesh % nVerticesArray(0:nHalos+1))

       block_ptr % mesh % nCellsArray(:) = nCellsCursor % array(:)
       block_ptr % mesh % nEdgesArray(:) = nEdgesCursor % array(:)
       block_ptr % mesh % nVerticesArray(:) = nVerticesCursor % array(:)

       ! Set block's local id
       block_ptr % localBlockID = localBlockID

       ! Set block's *Solve dimensions
       block_ptr % mesh % nCellsSolve = nCellsSolve_0Halo
       block_ptr % mesh % nEdgesSolve = nEdgesSolve_0Halo
       block_ptr % mesh % nVerticesSolve = nVerticesSolve_0Halo

       ! Set block's 0 halo indices
       block_ptr % mesh % indexToCellID % array(1:nCellsSolve_0Halo) = indexToCellCursor % array(1:nCellsSolve_0Halo)
       block_ptr % mesh % indexToEdgeID % array(1:nEdgesSolve_0Halo) = indexToEdgeCursor % array(1:nEdgesSolve_0Halo)
       block_ptr % mesh % indexToVertexID % array(1:nVerticesSolve_0Halo) = indexToVertexCursor % array(1:nVerticesSolve_0Halo)

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

       ! Advance cursors
       block_ptr => block_ptr % next
       nCellsCursor => nCellsCursor % next
       nEdgesCursor => nEdgesCursor % next
       nVerticesCursor => nVerticesCursor % next
       indexToCellCursor => indexToCellCursor % next
       indexToEdgeCursor => indexToEdgeCursor % next
       indexToVertexCursor => indextoVertexcursor % next
     end do

     ! Link fields between blocks
     block_ptr => blocklist
     do while(associated(block_ptr))
       call mpas_create_field_links(block_ptr)

       block_ptr => block_ptr % next
     end do
   end subroutine mpas_block_creator_finalize_block_init!}}}

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

     ! Loop over blocks
     block_ptr => blocklist
     do while(associated(block_ptr))
       !
       ! Rename vertices in cellsOnCell, edgesOnCell, etc. to local indices
       !
       allocate(cellIDSorted(2, block_ptr % mesh % nCells))
       allocate(edgeIDSorted(2, block_ptr % mesh % nEdges))
       allocate(vertexIDSorted(2, block_ptr % mesh % nVertices))
 
       do i=1,block_ptr % mesh % nCells
         cellIDSorted(1,i) = block_ptr % mesh % indexToCellID % array(i)
         cellIDSorted(2,i) = i
       end do
       call mpas_quicksort(block_ptr % mesh % nCells, cellIDSorted)
 
       do i=1,block_ptr % mesh % nEdges
         edgeIDSorted(1,i) = block_ptr % mesh % indexToEdgeID % array(i)
         edgeIDSorted(2,i) = i
       end do
       call mpas_quicksort(block_ptr % mesh % nEdges, edgeIDSorted)
 
       do i=1,block_ptr % mesh % nVertices
         vertexIDSorted(1,i) = block_ptr % mesh % indexToVertexID % array(i)
         vertexIDSorted(2,i) = i
       end do
       call mpas_quicksort(block_ptr % mesh % nVertices, vertexIDSorted)
 
 
       do i=1,block_ptr % mesh % nCells
         do j=1,block_ptr % mesh % nEdgesOnCell % array(i)
           k = mpas_binary_search(cellIDSorted, 2, 1, block_ptr % mesh % nCells, &
                                  block_ptr % mesh % cellsOnCell % array(j,i))
           if (k <= block_ptr % mesh % nCells) then
             block_ptr % mesh % cellsOnCell % array(j,i) = cellIDSorted(2,k)
           else
             block_ptr % mesh % cellsOnCell % array(j,i) = block_ptr % mesh % nCells + 1
           end if
 
           k = mpas_binary_search(edgeIDSorted, 2, 1, block_ptr % mesh % nEdges, &
                                  block_ptr % mesh % edgesOnCell % array(j,i))
           if (k <= block_ptr % mesh % nEdges) then
             block_ptr % mesh % edgesOnCell % array(j,i) = edgeIDSorted(2,k)
           else
             block_ptr % mesh % edgesOnCell % array(j,i) = block_ptr % mesh % nEdges + 1
           end if
  
           k = mpas_binary_search(vertexIDSorted, 2, 1, block_ptr % mesh % nVertices, &
                                  block_ptr % mesh % verticesOnCell % array(j,i))
           if (k <= block_ptr % mesh % nVertices) then
             block_ptr % mesh % verticesOnCell % array(j,i) = vertexIDSorted(2,k)
           else
             block_ptr % mesh % verticesOnCell % array(j,i) = block_ptr % mesh % nVertices + 1
           end if
         end do
       end do
  
       do i=1,block_ptr % mesh % nEdges
         do j=1,2
  
           k = mpas_binary_search(cellIDSorted, 2, 1, block_ptr % mesh % nCells, &
                                  block_ptr % mesh % cellsOnEdge % array(j,i))
           if (k <= block_ptr % mesh % nCells) then
             block_ptr % mesh % cellsOnEdge % array(j,i) = cellIDSorted(2,k)
           else
             block_ptr % mesh % cellsOnEdge % array(j,i) = block_ptr % mesh % nCells + 1
           end if
  
           k = mpas_binary_search(vertexIDSorted, 2, 1, block_ptr % mesh % nVertices, &
                                  block_ptr % mesh % verticesOnEdge % array(j,i))
           if (k <= block_ptr % mesh % nVertices) then
             block_ptr % mesh % verticesOnEdge % array(j,i) = vertexIDSorted(2,k)
           else
             block_ptr % mesh % verticesOnEdge % array(j,i) = block_ptr % mesh % nVertices + 1
           end if
  
         end do
  
         do j=1,block_ptr % mesh % nEdgesOnEdge % array(i)
  
           k = mpas_binary_search(edgeIDSorted, 2, 1, block_ptr % mesh % nEdges, &
                                  block_ptr % mesh % edgesOnEdge % array(j,i))
           if (k <= block_ptr % mesh % nEdges) then
             block_ptr % mesh % edgesOnEdge % array(j,i) = edgeIDSorted(2,k)
           else
             block_ptr % mesh % edgesOnEdge % array(j,i) = block_ptr % mesh % nEdges + 1
           end if
         end do
       end do
  
       do i=1,block_ptr % mesh % nVertices
         do j=1,block_ptr % mesh % vertexDegree
  
           k = mpas_binary_search(cellIDSorted, 2, 1, block_ptr % mesh % nCells, &
                                  block_ptr % mesh % cellsOnVertex % array(j,i))
           if (k <= block_ptr % mesh % nCells) then
             block_ptr % mesh % cellsOnVertex % array(j,i) = cellIDSorted(2,k)
           else
             block_ptr % mesh % cellsOnVertex % array(j,i) = block_ptr % mesh % nCells + 1
           end if
  
           k = mpas_binary_search(edgeIDSorted, 2, 1, block_ptr % mesh % nEdges, &
                             block_ptr % mesh % edgesOnVertex % array(j,i))
           if (k <= block_ptr % mesh % nEdges) then
             block_ptr % mesh % edgesOnVertex % array(j,i) = edgeIDSorted(2,k)
           else
             block_ptr % mesh % edgesOnVertex % array(j,i) = block_ptr % mesh % nEdges + 1
           end if
         end do
       end do
  
       deallocate(cellIDSorted)
       deallocate(edgeIDSorted)
       deallocate(vertexIDSorted)

       block_ptr => block_ptr % next
     end do

   end subroutine mpas_block_creator_reindex_block_fields!}}}

end module mpas_block_creator
