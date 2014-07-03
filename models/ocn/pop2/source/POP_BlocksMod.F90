!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_BlocksMod

!BOP
! !MODULE: POP_BlocksMod
!
! !DESCRIPTION: 
!  This module contains data types and tools for decomposing a global
!  horizontal domain into a set of 2d blocks.  It contains a data type 
!  for describing each block and contains routines for creating and 
!  querying the block decomposition for a global domain.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_BlocksMod.F90 69 2007-10-05 20:51:42Z pwjones $
!  2006-08-14: Phil Jones, John Dennis
!    new blocks module with new naming convention 
!    also generalized the Get and Set routines and added a
!      number of active points to block structure for support
!      of improvements by John Dennis (NCAR) elsewhere
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_DomainSizeMod

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: POP_block   ! block data type
      integer (POP_i4) ::   &
         blockID           ,&! global block number
         localID           ,&! local address of block in current distrib
         ib, ie, jb, je    ,&! begin,end indices for physical domain
         iBlock, jBlock    ,&! cartesian i,j position for bloc
         nxGlobal, nyGlobal,&! global domain extents
         numActivePoints     ! number of actual active points in block

      logical (POP_logical) :: &
         tripole             ! flag is true if block is at tripole bndy

      integer (POP_i4), dimension(:), pointer :: &
         iGlobal, jGlobal    ! global domain location for each point
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_BlocksCreate,       &
             POP_BlocksDestroy,      &
             POP_BlocksSet,          &
             POP_BlocksGetBlock,     &
             POP_BlocksGetNbrID,     & 
             POP_BlocksGetBlockInfo

! !DEFINED PARAMETERS:
  

   integer (POP_i4), parameter, public :: &
      POP_haloWidth = 2       ! number of ghost cells around each block

   ! size of block domain in i,j direction including ghost cells
   integer (POP_i4), parameter, public :: &
      POP_nxBlock = POP_blockSizeX + 2*POP_haloWidth,   &
      POP_nyBlock = POP_blockSizeY + 2*POP_haloWidth

   ! predefined directions for neighbor id routine
   integer (POP_i4), parameter, public :: &
      POP_blocksNorth          =  1,      & ! (i  ,j+1)
      POP_blocksSouth          =  2,      & ! (i  ,j-1)
      POP_blocksEast           =  3,      & ! (i+1,j  )
      POP_blocksWest           =  4,      & ! (i-1,j  )
      POP_blocksNorthEast      =  5,      & ! (i+1,j+1)
      POP_blocksNorthWest      =  6,      & ! (i-1,j+1)
      POP_blocksSouthEast      =  7,      & ! (i+1,j-1)
      POP_blocksSouthWest      =  8         ! (i-1,j-1)
   integer (POP_i4), parameter, public :: &
      POP_blocksNorth2         =  9,      & ! (i  ,j+2)
      POP_blocksSouth2         = 10,      & ! (i  ,j-2)
      POP_blocksEast2          = 11,      & ! (i+2,j  )
      POP_blocksWest2          = 12,      & ! (i-2,j  )
      POP_blocksNorthEast2     = 13,      & ! (i+2,j+2)
      POP_blocksNorthWest2     = 14,      & ! (i-2,j+2)
      POP_blocksSouthEast2     = 15,      & ! (i+2,j-2)
      POP_blocksSouthWest2     = 16         ! (i-2,j-2)
   integer (POP_i4), parameter, public :: &
      POP_blocksEastNorthEast  = 17,      & ! (i+2,j+1)
      POP_blocksEastSouthEast  = 18,      & ! (i+2,j-1)
      POP_blocksWestNorthWest  = 19,      & ! (i-2,j+1)
      POP_blocksWestSouthWest  = 20,      & ! (i-2,j-1)
      POP_blocksNorthNorthEast = 21,      & ! (i+1,j-2)
      POP_blocksSouthSouthEast = 22,      & ! (i+1,j-2)
      POP_blocksNorthNorthWest = 23,      & ! (i-1,j+2)
      POP_blocksSouthSouthWest = 24         ! (i-1,j-2)

! !PUBLIC DATA MEMBERS:

   integer (POP_i4), public :: &
      POP_numBlocks,    &! total number of blocks in decomposition
      POP_numBlocksX,   &! tot num blocks in i direction
      POP_numBlocksY     ! tot num blocks in j direction

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   type (POP_block), dimension(:), allocatable :: &
      allBlocks         ! block information for all blocks in domain

   integer (POP_i4), dimension(:,:),allocatable :: &
      allBlocksIJ       ! block index stored in Cartesian order
                        !   useful for determining block index
                        !   of neighbor blocks

   integer (POP_i4), dimension(:,:), allocatable, target :: &
      allIGlobal,         &! global i index for each point in each block
      allJGlobal           ! global j index for each point in each block

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_BlocksCreate
! !INTERFACE:

 subroutine POP_BlocksCreate(nxGlobal, nyGlobal, &
                             ewBoundaryType, nsBoundaryType, errorCode)

! !DESCRIPTION:
!  This subroutine decomposes the global domain into blocks and
!  fills the data structures with all the necessary block information.
!
! !REVISION HISTORY: 
!  same as module
!
! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      nxGlobal, nyGlobal           ! global domain size in x,y

   character (*), intent(in) :: &
      ewBoundaryType,  &! type of boundary in logical east-west dir
      nsBoundaryType    ! type of boundary in logical north-south dir

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) ::    &
      blockID,            &! block id for linear list of blocks
      i, j,               &! loop indices
      iBlock, jBlock,     &! block loop indices
      is, ie, js, je,     &! temp start, end indices
      istat                ! status flag for allocates

!----------------------------------------------------------------------
!
!  compute number of blocks and cartesian decomposition of blocks
!  if the requested block size does not divide the global domain
!  size evenly, add additional block space to accomodate padding
!
!----------------------------------------------------------------------

   errorCode = POP_success

   POP_numBlocksX   = (nxGlobal-1)/POP_blockSizeX + 1
   POP_numBlocksY   = (nyGlobal-1)/POP_blockSizeY + 1
   POP_numBlocks = POP_numBlocksX*POP_numBlocksY

!----------------------------------------------------------------------
!
!  allocate block arrays
!
!----------------------------------------------------------------------

   allocate(allBlocks(POP_numBlocks), stat=istat)
   if (istat /= 0) then
      call POP_ErrorSet(errorCode, &
                        'POP_BlocksCreate: error allocating allBlocks')
      return
   endif

   allocate(allIGlobal(POP_nxBlock,POP_numBlocks), &
            allJGlobal(POP_nyBlock,POP_numBlocks), stat=istat)
   if (istat /= 0) then
      call POP_ErrorSet(errorCode, &
                        'POP_BlocksCreate: error allocating ijGlobal')
      return
   endif

   allocate(allBlocksIJ(POP_numBlocksX,POP_numBlocksY), stat=istat)
   if (istat /= 0) then
      call POP_ErrorSet(errorCode, &
                        'POP_BlocksCreate: error allocating allBlocksIJ')
      return
   endif

!----------------------------------------------------------------------
!
!  fill block data structures for all blocks in domain
!
!----------------------------------------------------------------------

   blockID = 0  ! initialize block ID

   do jBlock=1,POP_numBlocksY

      !*** determine start, end of global physical domain in j 
      !*** direction for this row of blocks

      js = (jBlock-1)*POP_blockSizeY + 1
      je = js + POP_blockSizeY - 1
      if (je > nyGlobal) je = nyGlobal ! pad array
      if (js > nyGlobal) then
         call POP_ErrorSet(errorCode, &
          'POP_BlocksCreate: Bad block decomp: POP_nyBlock too large?')
         return
      endif

      do iBlock=1,POP_numBlocksX

         blockID = blockID + 1  ! increment global block id

         !*** determine start, end of global physical domain in i 
         !*** direction for this row of blocks

         is = (iBlock-1)*POP_blockSizeX + 1
         ie = is + POP_blockSizeX - 1
         if (ie > nxGlobal) ie = nxGlobal
         if (is > nxGlobal) then
            call POP_ErrorSet(errorCode, &
              'POP_BlocksCreate: Bad block decomp: POP_nxBlock too large?')
            return
         endif

         allBlocks(blockID)%blockID  = blockID
         allBlocks(blockID)%iBlock   = iBlock
         allBlocks(blockID)%jBlock   = jBlock
         allBlocks(blockID)%nxGlobal = nxGlobal
         allBlocks(blockID)%nyGlobal = nyGlobal

         if (jBlock == POP_numBlocksY .and. &
             nsBoundaryType == 'tripole') then
            allBlocks(blockID)%tripole = .true.
         else
            allBlocks(blockID)%tripole = .false.
         endif

         !*** set default values for start and end logical indices
         !*** for the physical domain.  ie,je may be reset later if
         !*** domains must be padded

         allBlocks(blockID)%ib       = POP_haloWidth + 1
         allBlocks(blockID)%jb       = POP_haloWidth + 1
         allBlocks(blockID)%ie       = POP_nxBlock - POP_haloWidth
         allBlocks(blockID)%je       = POP_nyBlock - POP_haloWidth

         allBlocksIJ(iBlock,jBlock) = blockID

         !*** for each point in this block, determine global indices

         do j=1,POP_nyBlock

            allJGlobal(j,blockID) = js - POP_haloWidth + j - 1

            !*** southern ghost cells

            if (allJGlobal(j,blockID) < 1) then
               select case (nsBoundaryType)
               case ('cyclic')
                  allJGlobal(j,blockID) = allJGlobal(j,blockID) + nyGlobal
               case ('closed')
                  allJGlobal(j,blockID) = 0
               case ('tripole')
                  allJGlobal(j,blockID) = 0
               case default
                  call POP_ErrorSet(errorCode, &
                     'POP_BlocksCreate: unknown n-s bndy type')
                  return
               end select
            endif

            !*** padding required

            if (allJGlobal(j,blockID) > nyGlobal + POP_haloWidth) then
               allJGlobal(j,blockID) = 0   ! padding

            !*** northern ghost cells

            else if (allJGlobal(j,blockID) > nyGlobal) then
               select case (nsBoundaryType)
               case ('cyclic')
                  allJGlobal(j,blockID) = allJGlobal(j,blockID) - nyGlobal
               case ('closed')
                  allJGlobal(j,blockID) = 0
               case ('tripole') ! use negative value to flag tripole
                  allJGlobal(j,blockID) = -allJGlobal(j,blockID)
               case default
                  call POP_ErrorSet(errorCode, &
                     'POP_BlocksCreate: unknown n-s bndy type')
                  return
               end select

            !*** set last physical point if padded domain

            else if (allJGlobal(j,blockID) == nyGlobal .and. &
                     j > allBlocks(blockID)%jb         .and. &
                     j < allBlocks(blockID)%je) then
               allBlocks(blockID)%je = j   ! last physical point in padded domain
            endif
         end do

         allBlocks(blockID)%jGlobal => allJGlobal(:,blockID)

         do i=1,POP_nxBlock
            allIGlobal(i,blockID) = is - POP_haloWidth + i - 1

            !*** western ghost cells

            if (allIGlobal(i,blockID) < 1) then
               select case (ewBoundaryType)
               case ('cyclic')
                  allIGlobal(i,blockID) = allIGlobal(i,blockID) + nxGlobal
               case ('closed')
                  allIGlobal(i,blockID) = 0
               case default
                  call POP_ErrorSet(errorCode, &
                     'POP_BlocksCreate: unknown e-w bndy type')
                  return
               end select
            endif

            !*** padded domain - fill padded region with zero

            if (allIGlobal(i,blockID) > nxGlobal + POP_haloWidth) then
               allIGlobal(i,blockID) = 0

            !*** eastern ghost cells

            else if (allIGlobal(i,blockID) > nxGlobal) then
               select case (ewBoundaryType)
               case ('cyclic')
                  allIGlobal(i,blockID) = allIGlobal(i,blockID) - nxGlobal
               case ('closed')
                  allIGlobal(i,blockID) = 0
               case default
                  call POP_ErrorSet(errorCode, &
                     'POP_BlocksCreate: unknown e-w bndy type')
                  return
               end select

            !*** last physical point in padded domain

            else if (allIGlobal(i,blockID) == nxGlobal .and. &
                     i > allBlocks(blockID)%ib         .and. &
                     i < allBlocks(blockID)%ie) then
               allBlocks(blockID)%ie = i
            endif
         end do

         allBlocks(blockID)%iGlobal => allIGlobal(:,blockID)

         !*** compute default number of active points
         !*** this will normally be reset later once the
         !*** physical domain is known and masked points are
         !*** eliminated

         allBlocks(blockID)%numActivePoints =                    &
            (allBlocks(blockID)%ie - allBlocks(blockID)%ib + 1)* &
            (allBlocks(blockID)%je - allBlocks(blockID)%jb + 1)
      end do
   end do

!EOC
!----------------------------------------------------------------------

end subroutine POP_BlocksCreate

!***********************************************************************
!BOP
! !IROUTINE: POP_BlocksGetBlock
! !INTERFACE:

 function POP_BlocksGetBlock(blockID, errorCode) &
                             result (outBlock)

! !DESCRIPTION:
!  This function returns the block data structure for the block
!  associated with the input block id.
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      blockID     ! global block id for requested block info

! !OUTPUT PARAMETERS:

   type (POP_block) :: &
      outBlock    ! block information returned for requested block

   integer (POP_i4), intent(out) :: &
      errorCode   ! returned error code

!EOP
!BOC
!----------------------------------------------------------------------
!
!  check for valid id.  if valid, return block info for requested block
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   if (blockID < 1 .or. blockID > POP_numBlocks) then
      call POP_ErrorSet(errorCode, &
              'POP_BlocksGetBlock: invalid blockID')
      return
   endif

   outBlock = allBlocks(blockID)

!----------------------------------------------------------------------
!EOC

 end function POP_BlocksGetBlock

!***********************************************************************
!BOP
! !IROUTINE: POP_BlocksGetNbrID
! !INTERFACE:

 function POP_BlocksGetNbrID(blockID, direction, iBoundary, jBoundary, &
                             errorCode) &
                             result (nbrID)

! !DESCRIPTION:
!  This function returns the block id of a neighboring block in a
!  requested direction.  Supported directions currently include:
!      POP\_blocksNorth             (i  ,j+1)
!      POP\_blocksSouth             (i  ,j-1)
!      POP\_blocksEast              (i+1,j  )
!      POP\_blocksWest              (i-1,j  )
!      POP\_blocksNorthEast         (i+1,j+1)
!      POP\_blocksNorthWest         (i-1,j+1)
!      POP\_blocksSouthEast         (i  ,j-1)
!      POP\_blocksSouthWest         (i-1,j-1)
!      POP\_blocksNorth2            (i  ,j+2)
!      POP\_blocksSouth2            (i  ,j-2)
!      POP\_blocksEast2             (i+2,j  )
!      POP\_blocksWest2             (i-2,j  )
!      POP\_blocksNorthEast2        (i+2,j+2)
!      POP\_blocksNorthWest2        (i-2,j+2)
!      POP\_blocksSouthEast2        (i+2,j-2)
!      POP\_blocksSouthWest2        (i-2,j-2)
!      POP\_blocksEastNorthEast     (i+2,j+1)
!      POP\_blocksEastSouthEast     (i+2,j-1)
!      POP\_blocksWestNorthWest     (i-2,j+1)
!      POP\_blocksWestSouthWest     (i-2,j-1)
!      POP\_blocksNorthNorthEast    (i+1,j-2)
!      POP\_blocksSouthSouthEast    (i+1,j-2)
!      POP\_blocksNorthNorthWest    (i-1,j+2)
!      POP\_blocksSouthSouthWest    (i-1,j-2)
!

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in)  :: &
      blockID,       &! id of block for which neighbor id requested
      direction       ! direction for which to look for neighbor -
                      !   must be one of the predefined module
                      !   variables for block direction

   character (*), intent(in) :: &
      iBoundary,     &! determines what to do at edges of domain
      jBoundary       !  options are - closed, cyclic, tripole

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

   integer (POP_i4) :: &
      nbrID           ! block ID of neighbor in requested dir

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
    
   integer (POP_i4) :: &
      iBlock, jBlock,  &! i,j block location of current block
      inbr,   jnbr      ! i,j block location of neighboring block

!----------------------------------------------------------------------
!
!  retrieve info for current block
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_BlocksGetBlockInfo(blockID, errorCode, &
                               iBlock=iBlock, jBlock=jBlock)

   if (errorCode /= POP_success) then
      call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: error getting block info')
      return
   endif

!----------------------------------------------------------------------
!
!  compute i,j block location of neighbor
!
!----------------------------------------------------------------------

   select case(direction)

   case (POP_blocksNorth)

      inbr = iBlock
      jnbr = jBlock + 1
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock + 1 
            jnbr = -jBlock
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksSouth)

      inbr = iBlock
      jnbr = jBlock - 1
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksEast )

      inbr = iBlock + 1
      jnbr = jBlock
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif

   case (POP_blocksWest )
      jnbr = jBlock
      inbr = iBlock - 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif

   case (POP_blocksNorthEast)

      inbr = iBlock + 1
      jnbr = jBlock + 1
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock 
            if (inbr == 0) inbr = POP_numBlocksX
            jnbr = -jBlock
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksNorthWest)

      inbr = iBlock - 1
      jnbr = jBlock + 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock + 2 
            if (inbr > POP_numBlocksX) inbr = 1
            jnbr = -jBlock
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksSouthEast )

      inbr = iBlock + 1
      jnbr = jBlock - 1
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksSouthWest )
      inbr = iBlock - 1
      jnbr = jBlock - 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksNorth2)

      inbr = iBlock
      jnbr = jBlock + 2
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr - POP_numBlocksY
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock + 1 
            jnbr = -(POP_numBlocksY - (jnbr - POP_numBlocksY - 1))
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksSouth2)

      inbr = iBlock
      jnbr = jBlock - 2
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY + jnbr
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksEast2)

      inbr = iBlock + 2
      jnbr = jBlock
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif

   case (POP_blocksWest2)
      jnbr = jBlock
      inbr = iBlock - 2
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX + inbr
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif

   case (POP_blocksNorthEast2)

      inbr = iBlock + 2
      jnbr = jBlock + 2
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr - POP_numBlocksY
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock - 1 
            if (inbr == 0) inbr = POP_numBlocksX
            jnbr = -(POP_numBlocksY - (jnbr - POP_numBlocksY - 1))
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksNorthWest2)

      inbr = iBlock - 2
      jnbr = jBlock + 2
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX + inbr
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock + 3
            if (inbr > POP_numBlocksX) inbr = 0
            jnbr = -(POP_numBlocksY - (jnbr - POP_numBlocksY - 1))
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksSouthEast2)

      inbr = iBlock + 2
      jnbr = jBlock - 2
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY + jnbr
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksSouthWest2)
      inbr = iBlock - 2
      jnbr = jBlock - 2
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX + inbr
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY + jnbr
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksEastNorthEast)

      inbr = iBlock + 2
      jnbr = jBlock + 1
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr - POP_numBlocksY
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock - 1 
            if (inbr == 0) inbr = POP_numBlocksX
            jnbr = -(POP_numBlocksY - (jnbr - POP_numBlocksY - 1))
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksWestNorthWest)

      inbr = iBlock - 2
      jnbr = jBlock + 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX + inbr
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr + POP_numBlocksY
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock + 3
            if (inbr > POP_numBlocksX) inbr = 0
            jnbr = -(POP_numBlocksY - (jnbr - POP_numBlocksY - 1))
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksEastSouthEast)

      inbr = iBlock + 2
      jnbr = jBlock - 1
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY + jnbr
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksWestSouthWest)
      inbr = iBlock - 2
      jnbr = jBlock - 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX + inbr
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY + jnbr
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksNorthNorthEast)

      inbr = iBlock + 1
      jnbr = jBlock + 2
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr - POP_numBlocksY
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock
            if (inbr == 0) inbr = POP_numBlocksX
            jnbr = -(POP_numBlocksY - (jnbr - POP_numBlocksY - 1))
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksNorthNorthWest)

      inbr = iBlock - 1
      jnbr = jBlock + 2
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX + inbr
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr > POP_numBlocksY) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr - POP_numBlocksY
         case ('tripole')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  POP_numBlocksX - iBlock + 2
            if (inbr > POP_numBlocksX) inbr = 0
            jnbr = -(POP_numBlocksY - (jnbr - POP_numBlocksY - 1))
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown north boundary')
         end select
      endif

   case (POP_blocksSouthSouthEast)

      inbr = iBlock + 1
      jnbr = jBlock - 2
      if (inbr > POP_numBlocksX) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - POP_numBlocksX
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown east boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY + jnbr
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case (POP_blocksSouthSouthWest)
      inbr = iBlock - 1
      jnbr = jBlock - 2
      if (inbr < 1) then
         select case(iBoundary)
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = POP_numBlocksX + inbr
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown west boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = POP_numBlocksY + jnbr
         case ('tripole')
            jnbr = 0
         case default
            call POP_ErrorSet(errorCode, &
               'POP_BlocksGetNbrID: unknown south boundary')
         end select
      endif

   case default

      call POP_ErrorSet(errorCode, &
          'POP_BlocksGetNbrID: unknown direction')
      return

   end select

!----------------------------------------------------------------------
!
!  now get block id for this neighbor block
!
!----------------------------------------------------------------------

   if (inbr > 0 .and. jnbr > 0) then
      nbrID = allBlocksIJ(inbr,jnbr)
   else if (inbr > 0 .and. jnbr < 0) then  ! tripole upper boundary
      !*** return negative value to flag tripole
      nbrID = -allBlocksIJ(inbr,abs(jnbr))
   else
      nbrID = 0   ! neighbor outside domain
   endif

!----------------------------------------------------------------------
!EOC

 end function POP_BlocksGetNbrID

!**********************************************************************
!BOP
! !IROUTINE: POP_BlocksGetBlockInfo
! !INTERFACE:

 subroutine POP_BlocksGetBlockInfo(blockID, errorCode,      &
                                   localID, ib, ie, jb, je, &
                                   numActivePoints,         &
                                   nxGlobal, nyGlobal,      &
                                   tripole,                 &
                                   iBlock, jBlock, iGlobal, jGlobal)

! !DESCRIPTION:
!  This routine returns requested parts of the block data type
!  for the block associated with the input block id
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      blockID   ! global block id for which parameters are requested

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode ! returned error code

   !(optional) parts of block data type to extract if requested

   integer (POP_i4), intent(out), optional :: &
      localID,          &! local id assigned to block in current distrb
      ib, ie, jb, je,   &! begin,end indices for physical domain
      iBlock, jBlock,   &! cartesian i,j position for block
      nxGlobal, nyGlobal, &! global domain size info
      numActivePoints    ! number of actual active points in block

   logical (POP_logical), intent(out), optional :: &
      tripole              ! flag is true if block on tripole bndy

   integer (POP_i4), dimension(:), pointer, optional :: &
      iGlobal, jGlobal     ! global domain location for each point

!EOP
!BOC
!----------------------------------------------------------------------
!
!  check for valid blockID
!
!----------------------------------------------------------------------

   errorCode = POP_Success
   if (blockID < 1 .or. blockID > POP_numBlocks) then
      call POP_ErrorSet(errorCode, &
               'POP_BlocksGetBlockInfo: invalid blockID')
      return
   endif

!----------------------------------------------------------------------
!
!  extract each component of data type if requested
!
!----------------------------------------------------------------------

   if (present(numActivePoints)) numActivePoints = &
              allBlocks(blockID)%numActivePoints
   if (present(localID)) localID =  allBlocks(blockID)%localID
   if (present(ib     )) ib      =  allBlocks(blockID)%ib
   if (present(ie     )) ie      =  allBlocks(blockID)%ie
   if (present(jb     )) jb      =  allBlocks(blockID)%jb
   if (present(je     )) je      =  allBlocks(blockID)%je
   if (present(iBlock )) iBlock  =  allBlocks(blockID)%iBlock
   if (present(jBlock )) jBlock  =  allBlocks(blockID)%jBlock
   if (present(tripole)) tripole =  allBlocks(blockID)%tripole
   if (present(nxGlobal)) nxGlobal=  allBlocks(blockID)%nxGlobal
   if (present(nyGlobal)) nyGlobal=  allBlocks(blockID)%nyGlobal
   if (present(iGlobal)) iGlobal => allBlocks(blockID)%iGlobal
   if (present(jGlobal)) jGlobal => allBlocks(blockID)%jGlobal

!----------------------------------------------------------------------
!EOC

 end subroutine POP_BlocksGetBlockInfo

!***********************************************************************
!BOP
! !IROUTINE: POP_BlocksSet
! !INTERFACE:

   subroutine POP_BlocksSet(blockID, errorCode,                 &
                            localID, numActivePoints,           &
                            ib, ie, jb, je, iBlock, jBlock,     &
                            nxGlobal, nyGlobal, tripole,        &
                            iGlobal, jGlobal)

! !DESCRIPTION:
!  This function sets or resets any parameter of a particular block
!  given the block id.  Most parameters are set correctly by the
!  BlocksCreate routine and should not be routinely reset.  The two 
!  exceptions are the localID (changes with each block distribution) 
!  and numActivePoints (should be set after details of the physical 
!  domain and ocean mask are known).
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   !*** required argument
   integer (POP_i4), intent(in) :: &
      blockID ! global block id for block to change

   !*** optional arguments for parameters to set
   integer (POP_i4), intent(in), optional ::   &
      localID,           &! local address of block in current distrib
      numActivePoints,   &! number of actual active points in block
      ib, ie, jb, je,    &! begin,end indices for physical domain
      iBlock, jBlock,    &! cartesian i,j position for bloc
      nxGlobal, nyGlobal  ! global domain information

   logical (POP_logical), intent(in), optional :: &
      tripole             ! flag is true if block on tripole bndy

   integer (POP_i4), dimension(:), intent(in), optional :: &
      iGlobal, jGlobal    ! global domain location for each point

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode               ! output error code

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: i,j  ! dummy loop indices

!----------------------------------------------------------------------
!
!  check for valid block id
!
!----------------------------------------------------------------------

   errorCode = POP_success

   if (blockID < 0 .or. blockID > POP_numBlocks) then
      call POP_ErrorSet(errorCode, &
         'POP_BlocksSet: invalid block id')
      return
   endif

!----------------------------------------------------------------------
!
!  for each parameter present, reset the value in the main allBlocks
!  array containing the block info
!
!----------------------------------------------------------------------

   if (present(numActivePoints)) &
                 allBlocks(blockID)%numActivePoints = numActivePoints
   if (present(localID)) allBlocks(blockID)%localID = localID
   if (present(ib))      allBlocks(blockID)%ib      = ib
   if (present(ie))      allBlocks(blockID)%ie      = ie
   if (present(jb))      allBlocks(blockID)%jb      = jb
   if (present(je))      allBlocks(blockID)%je      = je
   if (present(iBlock))  allBlocks(blockID)%iBlock  = iBlock
   if (present(jBlock))  allBlocks(blockID)%jBlock  = jBlock
   if (present(tripole)) allBlocks(blockID)%tripole = tripole
   if (present(nxGlobal)) allBlocks(blockID)%nxGlobal = nxGlobal
   if (present(nyGlobal)) allBlocks(blockID)%nyGlobal = nyGlobal

   if (present(iGlobal)) then
      if (size(iGlobal) == POP_nxBlock) then
         do i=1,POP_nxBlock
            allIGlobal(i,blockID) = iGlobal(i) 
         end do
      else
         call POP_ErrorSet(errorCode, &
             'POP_BlocksSet: iGlobal not correct size')
         return
      endif
   endif

   if (present(jGlobal)) then
      if (size(jGlobal) == POP_nyBlock) then
         do j=1,POP_nyBlock
            allIGlobal(j,blockID) = iGlobal(j) 
         end do
      else
         call POP_ErrorSet(errorCode, &
             'POP_BlocksSet: jGlobal not correct size')
         return
      endif
   endif

!EOC
!----------------------------------------------------------------------

  end subroutine POP_BlocksSet

!**********************************************************************
!BOP
! !IROUTINE: POP_BlocksDestroy
! !INTERFACE:

 subroutine POP_BlocksDestroy(errorCode)

! !DESCRIPTION:
!  This subroutine deallocates all arrays allocated by this module.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: istat  ! status flag for deallocate

!----------------------------------------------------------------------
!
!  deallocate arrays
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   deallocate(allBlocks, stat=istat)
   if (istat /= 0) then
      call POP_ErrorSet(errorCode, &
         'POP_BlocksDestroy: error deallocating allBlocks')
      return
   endif

   deallocate(allBlocksIJ, stat=istat)
   if (istat /= 0) then
      call POP_ErrorSet(errorCode, &
         'POP_BlocksDestroy: error deallocating allBlocksIJ')
      return
   endif


   deallocate(allIGlobal, allJGlobal, stat=istat)
   if (istat /= 0) then
      call POP_ErrorSet(errorCode, &
         'POP_BlocksDestroy: error deallocating ijGlobal')
      return
   endif

!EOC
!----------------------------------------------------------------------

 end subroutine POP_BlocksDestroy

!***********************************************************************

 end module POP_BlocksMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
