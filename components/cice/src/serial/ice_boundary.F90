!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_boundary

 module ice_boundary

! !DESCRIPTION:
!  This module contains data types and routines for updating halo
!  regions (ghost cells) 
!
! !REVISION HISTORY:
!  SVN:$Id$
!  2007-07-19: Phil Jones, Yoshi Yoshida, John Dennis
!              new naming conventions, optimizations during
!              initialization, true multi-dimensional updates 
!              (rather than serial call to two-dimensional updates), 
!              fixes for non-existent blocks
!  2008-01-28: Elizabeth Hunke replaced old routines with new POP
!              infrastructure

! !USES:

   use ice_kinds_mod
   use ice_communicate, only: my_task
   use ice_constants, only: field_type_scalar, &
           field_type_vector, field_type_angle, &
           field_loc_center,  field_loc_NEcorner, &
           field_loc_Nface, field_loc_Eface
   use ice_global_reductions, only: global_maxval
   use ice_exit

   use ice_blocks, only: nx_block, ny_block, nghost, &
           nblocks_tot, ice_blocksNorth, &
           ice_blocksSouth, ice_blocksEast, ice_blocksWest, &
           ice_blocksEast2, ice_blocksWest2, &
           ice_blocksNorthEast, ice_blocksNorthWest, &
           ice_blocksEastNorthEast, ice_blocksWestNorthWest, &
           ice_blocksSouthEast, ice_blocksSouthWest, &
           ice_blocksGetNbrID, get_block_parameter
   use ice_distribution, only: distrb, &
           ice_distributionGetBlockLoc, ice_distributionGet

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: ice_halo
      integer (int_kind) ::  &
         communicator,     &! communicator to use for update messages
         numLocalCopies,   &! num local copies for halo update
         tripoleRows        ! number of rows in tripole buffer

      logical (log_kind) ::  &
         tripoleTFlag       ! NS boundary is a tripole T-fold

      integer (int_kind), dimension(:,:), pointer :: &
         srcLocalAddr,     &! src addresses for each local copy
         dstLocalAddr       ! dst addresses for each local copy

   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: ice_HaloCreate,  &
             ice_HaloUpdate,  &
             ice_HaloExtrapolate

   interface ice_HaloUpdate  ! generic interface
      module procedure ice_HaloUpdate2DR8, &
                       ice_HaloUpdate2DR4, &
                       ice_HaloUpdate2DI4, &
                       ice_HaloUpdate3DR8, &
                       ice_HaloUpdate3DR4, &
                       ice_HaloUpdate3DI4, &
                       ice_HaloUpdate4DR8, &
                       ice_HaloUpdate4DR4, &
                       ice_HaloUpdate4DI4
   end interface

   interface ice_HaloExtrapolate  ! generic interface
      module procedure ice_HaloExtrapolate2DR8 !, &
!                       ice_HaloExtrapolate2DR4, &  ! not yet
!                       ice_HaloExtrapolate2DI4, &  ! implemented
   end interface

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  global buffers for tripole boundary
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), allocatable :: &
      bufTripoleI4

   real (real_kind), dimension(:,:), allocatable :: &
      bufTripoleR4

   real (dbl_kind), dimension(:,:), allocatable :: &
      bufTripoleR8

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloCreate
! !INTERFACE:

 function ice_HaloCreate(dist, nsBoundaryType, ewBoundaryType, &
                         nxGlobal)  result(halo)

! !DESCRIPTION:
!  This routine creates a halo type with info necessary for
!  performing a halo (ghost cell) update. This info is computed
!  based on the input block distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      dist             ! distribution of blocks across procs

   character (*), intent(in) :: &
      nsBoundaryType,   &! type of boundary to use in logical ns dir
      ewBoundaryType     ! type of boundary to use in logical ew dir

   integer (int_kind), intent(in) :: &
      nxGlobal           ! global grid extent for tripole grids

! !OUTPUT PARAMETERS:

   type (ice_halo) :: &
      halo               ! a new halo type with info for halo updates

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::             &
      i,j,k,l,n,m,&
      istat,                       &! allocate status flag
      numProcs,                    &! num of processors involved
      communicator,                &! communicator for message passing
      iblock,                      &! block counter
      eastBlock, westBlock,        &! block id  east,  west neighbors
      northBlock, southBlock,      &! block id north, south neighbors
      neBlock, nwBlock,            &! block id northeast, northwest nbrs
      seBlock, swBlock,            &! block id southeast, southwest nbrs
      srcProc, dstProc,            &! source, dest processor locations
      srcLocalID, dstLocalID,      &! local block index of src,dst blocks
      blockSizeX,                  &! size of default physical domain in X 
      blockSizeY,                  &! size of default physical domain in Y 
      eastMsgSize, westMsgSize,    &! nominal sizes for e-w msgs
      northMsgSize, southMsgSize,  &! nominal sizes for n-s msgs
      tripoleRows,                 &! number of rows in tripole buffer
      cornerMsgSize, msgSize        ! nominal size for corner msg

   integer (int_kind), dimension(:), allocatable :: &
      sendCount, recvCount          ! count number of words to each proc

   logical (log_kind) :: &
      tripoleFlag,          &! flag for allocating tripole buffers
      tripoleBlock,         &! flag for identifying north tripole blocks
      tripoleTFlag           ! flag for processing tripole buffer as T-fold

!-----------------------------------------------------------------------
!
!  Initialize some useful variables and return if this task not
!  in the current distribution.
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist,          &
                            nprocs = numProcs,       &
                            communicator = communicator)

   if (my_task >= numProcs) return

   halo%communicator = communicator

   blockSizeX = nx_block - 2*nghost
   blockSizeY = ny_block - 2*nghost
   eastMsgSize  = nghost*blockSizeY
   westMsgSize  = nghost*blockSizeY
   southMsgSize = nghost*blockSizeX
   cornerMsgSize = nghost*nghost
   tripoleRows = nghost+1

   if (nsBoundaryType == 'tripole' .or. nsBoundaryType == 'tripoleT') then
      tripoleFlag = .true.
      tripoleTFlag = (nsBoundaryType == 'tripoleT')
      if (tripoleTflag) tripoleRows = tripoleRows+1
      northMsgSize = tripoleRows*blockSizeX

      !*** allocate tripole message buffers if not already done

      if (.not. allocated(bufTripoleR8)) then
         allocate (bufTripoleI4(nxGlobal, tripoleRows), &
                   bufTripoleR4(nxGlobal, tripoleRows), &
                   bufTripoleR8(nxGlobal, tripoleRows), &
                   stat=istat)

         if (istat > 0) then
            call abort_ice( &
               'ice_HaloCreate: error allocating tripole buffers')
            return
         endif
      endif

   else
      tripoleFlag = .false.
      tripoleTFlag = .false.
      northMsgSize = nghost*blockSizeX
   endif
   halo%tripoleTFlag = tripoleTFlag
   halo%tripoleRows = tripoleRows

!-----------------------------------------------------------------------
!
!  Count the number of messages to send/recv from each processor
!  and number of words in each message.  These quantities are
!  necessary for allocating future arrays.
!
!-----------------------------------------------------------------------

   allocate (sendCount(numProcs), recvCount(numProcs), stat=istat)

   if (istat > 0) then
      call abort_ice('ice_HaloCreate: error allocating count arrays')
      return
   endif

   sendCount  = 0
   recvCount  = 0

   msgCountLoop: do iblock=1,nblocks_tot

      call ice_distributionGetBlockLoc(dist, iblock, srcProc, &
                                       srcLocalID)

      !*** find north neighbor block and add to message count
      !***  also set tripole block flag for later special cases

      northBlock = ice_blocksGetNbrID(iblock, ice_blocksNorth,        &
                                      ewBoundaryType, nsBoundaryType)
      if (northBlock > 0) then
         tripoleBlock = .false.
         call ice_distributionGetBlockLoc(dist, northBlock, dstProc, &
                                          dstLocalID)
      else if (northBlock < 0) then ! tripole north row, count block
         tripoleBlock = .true.
         call ice_distributionGetBlockLoc(dist, abs(northBlock), &
                                 dstProc, dstLocalID)
      else
         tripoleBlock = .false.
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,           &
                                     srcProc, dstProc, northMsgSize)

      !*** if a tripole boundary block, also create a local
      !*** message into and out of tripole buffer 

      if (tripoleBlock) then
         !*** copy in
         call ice_HaloIncrementMsgCount(sendCount, recvCount,        &
                                        srcProc, srcProc,            &
                                        northMsgSize)

         !*** copy out of tripole buffer - includes halo
         call ice_HaloIncrementMsgCount(sendCount, recvCount,     &
                                           srcProc, srcProc,      & 
                                           (nghost+1)*nx_block)
      endif

      !*** find south neighbor block and add to message count

      southBlock = ice_blocksGetNbrID(iblock, ice_blocksSouth,        &
                                      ewBoundaryType, nsBoundaryType)

      if (southBlock > 0) then
         call ice_distributionGetBlockLoc(dist, southBlock, dstProc, &
                                          dstLocalID)
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,           &
                                     srcProc, dstProc, southMsgSize)

      !*** find east neighbor block and add to message count

      eastBlock = ice_blocksGetNbrID(iblock, ice_blocksEast,         &
                                     ewBoundaryType, nsBoundaryType)

      if (eastBlock > 0) then
         call ice_distributionGetBlockLoc(dist, eastBlock, dstProc, &
                                          dstLocalID)
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, eastMsgSize)

      !*** if a tripole boundary block, non-local east neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

!echmod      if (tripoleBlock .and. dstProc /= srcProc) then
      if (tripoleBlock) then
         call ice_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, northMsgSize)
      endif

      !*** find west neighbor block and add to message count

      westBlock = ice_blocksGetNbrID(iblock, ice_blocksWest,         &
                                     ewBoundaryType, nsBoundaryType)

      if (westBlock > 0) then
         call ice_distributionGetBlockLoc(dist, westBlock, dstProc, &
                                          dstLocalID)
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, westMsgSize)

      !*** if a tripole boundary block, non-local west neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

!echmod      if (tripoleBlock .and. dstProc /= srcProc) then
      if (tripoleBlock) then
         call ice_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, northMsgSize)
      endif

      !*** find northeast neighbor block and add to message count

      neBlock = ice_blocksGetNbrID(iblock, ice_blocksNorthEast,    &
                                   ewBoundaryType, nsBoundaryType)

      if (neBlock > 0) then
         msgSize = cornerMsgSize  ! normal corner message 

         call ice_distributionGetBlockLoc(dist, neBlock, dstProc, &
                                          dstLocalID)

      else if (neBlock < 0) then ! tripole north row
         msgSize = northMsgSize  ! tripole needs whole top row of block

         call ice_distributionGetBlockLoc(dist, abs(neBlock), dstProc, &
                                          dstLocalID)
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,      &
                                     srcProc, dstProc, msgSize)

      !*** find northwest neighbor block and add to message count

      nwBlock = ice_blocksGetNbrID(iblock, ice_blocksNorthWest,    &
                                   ewBoundaryType, nsBoundaryType)

      if (nwBlock > 0) then
         msgSize = cornerMsgSize ! normal NE corner update

         call ice_distributionGetBlockLoc(dist, nwBlock, dstProc, &
                                          dstLocalID)

      else if (nwBlock < 0) then ! tripole north row, count block
         msgSize = northMsgSize ! tripole NE corner update - entire row needed

         call ice_distributionGetBlockLoc(dist, abs(nwBlock), dstProc, &
                                          dstLocalID)

      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,      &
                                     srcProc, dstProc, msgSize)

      !*** find southeast neighbor block and add to message count

      seBlock = ice_blocksGetNbrID(iblock, ice_blocksSouthEast,    &
                                   ewBoundaryType, nsBoundaryType)

      if (seBlock > 0) then
         call ice_distributionGetBlockLoc(dist, seBlock, dstProc, &
                                          dstLocalID)
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,            &
                                     srcProc, dstProc, cornerMsgSize)

      !*** find southwest neighbor block and add to message count

      swBlock = ice_blocksGetNbrID(iblock, ice_blocksSouthWest,    &
                                   ewBoundaryType, nsBoundaryType)

      if (swBlock > 0) then
         call ice_distributionGetBlockLoc(dist, swBlock, dstProc, &
                                          dstLocalID)
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,            &
                                     srcProc, dstProc, cornerMsgSize)

   end do msgCountLoop

!-----------------------------------------------------------------------
!
!  if messages are received from the same processor, the message is 
!  actually a local copy - count them and reset to zero
!
!-----------------------------------------------------------------------

   halo%numLocalCopies = recvCount(my_task+1)

   sendCount(my_task+1) = 0
   recvCount(my_task+1) = 0

!-----------------------------------------------------------------------
!
!  allocate arrays for message information and initialize
!
!-----------------------------------------------------------------------

   allocate(halo%srcLocalAddr(3,halo%numLocalCopies), &
            halo%dstLocalAddr(3,halo%numLocalCopies), &
            stat = istat)

   if (istat > 0) then
      call abort_ice( &
         'ice_HaloCreate: error allocating halo message info arrays')
      return
   endif

   halo%srcLocalAddr = 0
   halo%dstLocalAddr = 0

   deallocate(sendCount, recvCount, stat=istat)

   if (istat > 0) then
      call abort_ice( &
         'ice_HaloCreate: error deallocating count arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  repeat loop through blocks but this time, determine all the
!  required message information for each message or local copy
!
!-----------------------------------------------------------------------

   !*** reset halo scalars to use as counters

   halo%numLocalCopies = 0

   msgConfigLoop: do iblock=1,nblocks_tot

      call ice_distributionGetBlockLoc(dist, iblock, srcProc, &
                                       srcLocalID)

      !*** find north neighbor block

      northBlock = ice_blocksGetNbrID(iblock, ice_blocksNorth,        &
                                      ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, northBlock, 'north')

      !*** set tripole flag and add two copies for inserting
      !*** and extracting info from the tripole buffer

      if (northBlock < 0) then
         tripoleBlock = .true.
         call ice_HaloMsgCreate(halo, dist, iblock, -iblock, 'north')
         call ice_HaloMsgCreate(halo, dist, -iblock, iblock, 'north')
      else
         tripoleBlock = .false.
      endif

      !*** find south neighbor block

      southBlock = ice_blocksGetNbrID(iblock, ice_blocksSouth,        &
                                      ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, southBlock, 'south')

      !*** find east neighbor block

      eastBlock = ice_blocksGetNbrID(iblock, ice_blocksEast,         &
                                     ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, eastBlock, 'east')

      !*** for tripole grids, send a north tripole message to
      !*** the east block to make sure enough information is
      !*** available for tripole manipulations

      if (tripoleBlock) then
         call ice_HaloMsgCreate(halo, dist, iblock, -eastBlock, 'north')
      endif

      !*** find west neighbor block

      westBlock = ice_blocksGetNbrID(iblock, ice_blocksWest,         &
                                     ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, westBlock, 'west')

      !*** for tripole grids, send a north tripole message to
      !*** the west block to make sure enough information is
      !*** available for tripole manipulations
 
      if (tripoleBlock) then
         call ice_HaloMsgCreate(halo, dist, iblock, -westBlock, 'north')
      endif

      !*** find northeast neighbor block

      neBlock = ice_blocksGetNbrID(iblock, ice_blocksNorthEast,    &
                                   ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, neBlock, 'northeast')

      !*** find northwest neighbor block

      nwBlock = ice_blocksGetNbrID(iblock, ice_blocksNorthWest,    &
                                   ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, nwBlock, 'northwest')

      !*** find southeast neighbor block

      seBlock = ice_blocksGetNbrID(iblock, ice_blocksSouthEast,    &
                                   ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, seBlock, 'southeast')

      !*** find southwest neighbor block

      swBlock = ice_blocksGetNbrID(iblock, ice_blocksSouthWest,    &
                                   ewBoundaryType, nsBoundaryType)

      call ice_HaloMsgCreate(halo, dist, iblock, swBlock, 'southwest')

   end do msgConfigLoop

!-----------------------------------------------------------------------
!EOC

 end function ice_HaloCreate

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate2DR8
! !INTERFACE:

 subroutine ice_HaloUpdate2DR8(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 2d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (dbl_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,n,nmsg,                &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (dbl_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_dbl_kind
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR8)) then
      nxGlobal = size(bufTripoleR8,dim=1)
      bufTripoleR8 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            array(iDst,jDst,dstBlock) = &
            array(iSrc,jSrc,srcBlock)
         else if (dstBlock < 0) then ! tripole copy into buffer
            bufTripoleR8(iDst,jDst) = &
            array(iSrc,jSrc,srcBlock)
         endif
      else if (srcBlock == 0) then
         array(iDst,jDst,dstBlock) = fill
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate2DR8: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripoleR8(i   ,halo%tripoleRows)
              x2 = bufTripoleR8(iDst,halo%tripoleRows)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripoleR8(i   ,halo%tripoleRows) = xavg
              bufTripoleR8(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripoleR8(i   ,halo%tripoleRows)
              x2 = bufTripoleR8(iDst,halo%tripoleRows)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripoleR8(i   ,halo%tripoleRows) = xavg
              bufTripoleR8(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate2DR8: Unknown field location')
        end select

      else ! tripole u-fold

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripoleR8(i   ,halo%tripoleRows)
              x2 = bufTripoleR8(iDst,halo%tripoleRows)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripoleR8(i   ,halo%tripoleRows) = xavg
              bufTripoleR8(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripoleR8(i   ,halo%tripoleRows)
              x2 = bufTripoleR8(iDst,halo%tripoleRows)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripoleR8(i   ,halo%tripoleRows) = xavg
              bufTripoleR8(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate2DR8: Unknown field location')
        end select

      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               array(iDst,jDst,dstBlock) = isign*bufTripoleR8(iSrc,jSrc)
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate2DR8

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate2DR4
! !INTERFACE:

 subroutine ice_HaloUpdate2DR4(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 2d horizontal arrays of single precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (real_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (real_kind), dimension(:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,n,nmsg,                &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (real_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_real_kind
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR4)) then
      nxGlobal = size(bufTripoleR4,dim=1)
      bufTripoleR4 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            array(iDst,jDst,dstBlock) = &
            array(iSrc,jSrc,srcBlock)
         else if (dstBlock < 0) then ! tripole copy into buffer
            bufTripoleR4(iDst,jDst) = &
            array(iSrc,jSrc,srcBlock)
         endif
      else if (srcBlock == 0) then
         array(iDst,jDst,dstBlock) = fill
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate2DR4: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripoleR4(i   ,halo%tripoleRows)
              x2 = bufTripoleR4(iDst,halo%tripoleRows)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripoleR4(i   ,halo%tripoleRows) = xavg
              bufTripoleR4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripoleR4(i   ,halo%tripoleRows)
              x2 = bufTripoleR4(iDst,halo%tripoleRows)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripoleR4(i   ,halo%tripoleRows) = xavg
              bufTripoleR4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate2DR4: Unknown field location')
        end select        

      else ! tripole u-fold
  
        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripoleR4(i   ,halo%tripoleRows)
              x2 = bufTripoleR4(iDst,halo%tripoleRows)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripoleR4(i   ,halo%tripoleRows) = xavg
              bufTripoleR4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripoleR4(i   ,halo%tripoleRows)
              x2 = bufTripoleR4(iDst,halo%tripoleRows)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripoleR4(i   ,halo%tripoleRows) = xavg
              bufTripoleR4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate2DR4: Unknown field location')
        end select

      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               array(iDst,jDst,dstBlock) = isign*bufTripoleR4(iSrc,jSrc)
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate2DR4

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate2DI4
! !INTERFACE:

 subroutine ice_HaloUpdate2DI4(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 2d horizontal integer arrays.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   integer (int_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,n,nmsg,                &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (int_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0_int_kind
   endif

   nxGlobal = 0
   if (allocated(bufTripoleI4)) then
      nxGlobal = size(bufTripoleI4,dim=1)
      bufTripoleI4 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            array(iDst,jDst,dstBlock) = &
            array(iSrc,jSrc,srcBlock)
         else if (dstBlock < 0) then ! tripole copy into buffer
            bufTripoleI4(iDst,jDst) = &
            array(iSrc,jSrc,srcBlock)
         endif
      else if (srcBlock == 0) then
         array(iDst,jDst,dstBlock) = fill
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate2DI4: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripoleI4(i   ,halo%tripoleRows)
              x2 = bufTripoleI4(iDst,halo%tripoleRows)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripoleI4(i   ,halo%tripoleRows) = xavg
              bufTripoleI4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripoleI4(i   ,halo%tripoleRows)
              x2 = bufTripoleI4(iDst,halo%tripoleRows)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripoleI4(i   ,halo%tripoleRows) = xavg
              bufTripoleI4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate2DI4: Unknown field location')
        end select

      else ! tripole u-fold  
  
        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripoleI4(i   ,halo%tripoleRows)
              x2 = bufTripoleI4(iDst,halo%tripoleRows)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripoleI4(i   ,halo%tripoleRows) = xavg
              bufTripoleI4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripoleI4(i   ,halo%tripoleRows)
              x2 = bufTripoleI4(iDst,halo%tripoleRows)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripoleI4(i   ,halo%tripoleRows) = xavg
              bufTripoleI4(iDst,halo%tripoleRows) = isign*xavg
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate2DI4: Unknown field location')
        end select

      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               array(iDst,jDst,dstBlock) = isign*bufTripoleI4(iSrc,jSrc)
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate2DI4

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate3DR8
! !INTERFACE:

 subroutine ice_HaloUpdate3DR8(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 3d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (dbl_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,n,nmsg,              &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of array in 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (dbl_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (dbl_kind), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_dbl_kind
   endif

   nz = size(array, dim=3)

   nxGlobal = 0
   if (allocated(bufTripoleR8)) then
      nxGlobal = size(bufTripoleR8,dim=1)
      allocate(bufTripole(nxGlobal,halo%tripoleRows,nz))
      bufTripole = fill
   endif 

!-----------------------------------------------------------------------
!
!  do local copies
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            do k=1,nz
               array(iDst,jDst,k,dstBlock) = &
               array(iSrc,jSrc,k,srcBlock)
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do k=1,nz
               bufTripole(iDst,jDst,k) = &
               array(iSrc,jSrc,k,srcBlock)
            end do
         endif
      else if (srcBlock == 0) then
         do k=1,nz
            array(iDst,jDst,k,dstBlock) = fill
         end do
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate3DR8: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do k=1,nz
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate3DR8: Unknown field location')
        end select 
  
      else ! tripole u-fold
  
        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate3DR8: Unknown field location')
        end select

      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               do k=1,nz
                  array(iDst,jDst,k,dstBlock) = isign*    &
                                  bufTripole(iSrc,jSrc,k)
               end do
            endif

         endif
      end do

   endif

   if (allocated(bufTripole)) deallocate(bufTripole)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate3DR8

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate3DR4
! !INTERFACE:

 subroutine ice_HaloUpdate3DR4(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 3d horizontal arrays of single precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (real_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (real_kind), dimension(:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,n,nmsg,              &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of array in 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (real_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (real_kind), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_real_kind
   endif

   nz = size(array, dim=3)

   nxGlobal = 0
   if (allocated(bufTripoleR4)) then
      nxGlobal = size(bufTripoleR4,dim=1)
      allocate(bufTripole(nxGlobal,halo%tripoleRows,nz))
      bufTripole = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            do k=1,nz
               array(iDst,jDst,k,dstBlock) = &
               array(iSrc,jSrc,k,srcBlock)
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do k=1,nz
               bufTripole(iDst,jDst,k) = &
               array(iSrc,jSrc,k,srcBlock)
            end do
         endif
      else if (srcBlock == 0) then
         do k=1,nz
            array(iDst,jDst,k,dstBlock) = fill
         end do
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate3DR4: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do k=1,nz
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate3DR4: Unknown field location')
        end select  
  
      else ! tripole u-fold  
  
        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate3DR4: Unknown field location')
        end select

      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               do k=1,nz
                  array(iDst,jDst,k,dstBlock) = isign*    &
                                  bufTripole(iSrc,jSrc,k)
               end do
            endif

         endif
      end do

   endif

   if (allocated(bufTripole)) deallocate(bufTripole)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate3DR4

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate3DI4
! !INTERFACE:

 subroutine ice_HaloUpdate3DI4(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 3d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   integer (int_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,n,nmsg,              &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of array in 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (int_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (int_kind), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0_int_kind
   endif

   nz = size(array, dim=3)

   nxGlobal = 0
   if (allocated(bufTripoleI4)) then
      nxGlobal = size(bufTripoleI4,dim=1)
      allocate(bufTripole(nxGlobal,halo%tripoleRows,nz))
      bufTripole = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            do k=1,nz
               array(iDst,jDst,k,dstBlock) = &
               array(iSrc,jSrc,k,srcBlock)
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do k=1,nz
               bufTripole(iDst,jDst,k) = &
               array(iSrc,jSrc,k,srcBlock)
            end do
         endif
      else if (srcBlock == 0) then
         do k=1,nz
            array(iDst,jDst,k,dstBlock) = fill
         end do
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then
 
      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate3DI4: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do k=1,nz
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate3DI4: Unknown field location')
        end select
  
      else ! tripole u-fold  

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k)
              x2 = bufTripole(iDst,halo%tripoleRows,k)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k) = xavg
              bufTripole(iDst,halo%tripoleRows,k) = isign*xavg
           end do
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate3DI4: Unknown field location')
        end select
 
      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               do k=1,nz
                  array(iDst,jDst,k,dstBlock) = isign*    &
                                  bufTripole(iSrc,jSrc,k)
               end do
            endif

         endif
      end do

   endif

   if (allocated(bufTripole)) deallocate(bufTripole)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate3DI4

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate4DR8
! !INTERFACE:

 subroutine ice_HaloUpdate4DR8(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 4d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (dbl_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,l,n,nmsg,            &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of array in 3rd,4th dimensions
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (dbl_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (dbl_kind), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_dbl_kind
   endif

   nz = size(array, dim=3)
   nt = size(array, dim=4)

   nxGlobal = 0
   if (allocated(bufTripoleR8)) then
      nxGlobal = size(bufTripoleR8,dim=1)
      allocate(bufTripole(nxGlobal,halo%tripoleRows,nz,nt))
      bufTripole = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            do l=1,nt
            do k=1,nz
               array(iDst,jDst,k,l,dstBlock) = &
               array(iSrc,jSrc,k,l,srcBlock)
            end do
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do l=1,nt
            do k=1,nz
               bufTripole(iDst,jDst,k,l) = &
               array(iSrc,jSrc,k,l,srcBlock)
            end do
            end do
         endif
      else if (srcBlock == 0) then
         do l=1,nt
         do k=1,nz
            array(iDst,jDst,k,l,dstBlock) = fill
         end do
         end do
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then
 
      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate4DR8: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do l=1,nt
           do k=1,nz
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate4DR8: Unknown field location')
        end select  
  
      else ! tripole u-fold  

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_dbl_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate4DR8: Unknown field location')
        end select
 
      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               do l=1,nt
               do k=1,nz
                  array(iDst,jDst,k,l,dstBlock) = isign*    &
                                  bufTripole(iSrc,jSrc,k,l)
               end do
               end do
            endif

         endif
      end do

   endif

   if (allocated(bufTripole)) deallocate(bufTripole)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate4DR8

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate4DR4
! !INTERFACE:

 subroutine ice_HaloUpdate4DR4(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 4d horizontal arrays of single precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (real_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (real_kind), dimension(:,:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,l,n,nmsg,            &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of array in 3rd,4th dimensions
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (real_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (real_kind), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_real_kind
   endif

   nz = size(array, dim=3)
   nt = size(array, dim=4)

   nxGlobal = 0
   if (allocated(bufTripoleR4)) then
      nxGlobal = size(bufTripoleR4,dim=1)
      allocate(bufTripole(nxGlobal,halo%tripoleRows,nz,nt))
      bufTripole = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            do l=1,nt
            do k=1,nz
               array(iDst,jDst,k,l,dstBlock) = &
               array(iSrc,jSrc,k,l,srcBlock)
            end do
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do l=1,nt
            do k=1,nz
               bufTripole(iDst,jDst,k,l) = &
               array(iSrc,jSrc,k,l,srcBlock)
            end do
            end do
         endif
      else if (srcBlock == 0) then
         do l=1,nt
         do k=1,nz
            array(iDst,jDst,k,l,dstBlock) = fill
         end do
         end do
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then
 
      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate4DR4: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do l=1,nt
           do k=1,nz
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate4DR4: Unknown field location')
        end select  
  
      else ! tripole u-fold  

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = 0.5_real_kind*(x1 + isign*x2)
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate4DR4: Unknown field location')
        end select
 
      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               do l=1,nt
               do k=1,nz
                  array(iDst,jDst,k,l,dstBlock) = isign*    &
                                  bufTripole(iSrc,jSrc,k,l)
               end do
               end do
            endif

         endif
      end do

   endif

   if (allocated(bufTripole)) deallocate(bufTripole)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate4DR4

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate4DI4
! !INTERFACE:

 subroutine ice_HaloUpdate4DI4(array, halo,                    &
                               fieldLoc, fieldKind, &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  POP\_HaloUpdate.  This routine is the specific interface
!  for 4d horizontal integer arrays.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT PARAMETERS:

   type (ice_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   integer (int_kind), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   integer (int_kind), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,l,n,nmsg,            &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of array in 3rd,4th dimensions
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (int_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (int_kind), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0_int_kind
   endif

   nz = size(array, dim=3)
   nt = size(array, dim=4)

   nxGlobal = 0
   if (allocated(bufTripoleI4)) then
      nxGlobal = size(bufTripoleI4,dim=1)
      allocate(bufTripole(nxGlobal,halo%tripoleRows,nz,nt))
      bufTripole = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies
!  if srcBlock is zero, that denotes an eliminated land block or a 
!    closed boundary where ghost cell values are undefined
!  if srcBlock is less than zero, the message is a copy out of the
!    tripole buffer and will be treated later
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numLocalCopies
      iSrc     = halo%srcLocalAddr(1,nmsg)
      jSrc     = halo%srcLocalAddr(2,nmsg)
      srcBlock = halo%srcLocalAddr(3,nmsg)
      iDst     = halo%dstLocalAddr(1,nmsg)
      jDst     = halo%dstLocalAddr(2,nmsg)
      dstBlock = halo%dstLocalAddr(3,nmsg)

      if (srcBlock > 0) then
         if (dstBlock > 0) then
            do l=1,nt
            do k=1,nz
               array(iDst,jDst,k,l,dstBlock) = &
               array(iSrc,jSrc,k,l,srcBlock)
            end do
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do l=1,nt
            do k=1,nz
               bufTripole(iDst,jDst,k,l) = &
               array(iSrc,jSrc,k,l,srcBlock)
            end do
            end do
         endif
      else if (srcBlock == 0) then
         do l=1,nt
         do k=1,nz
            array(iDst,jDst,k,l,dstBlock) = fill
         end do
         end do
      endif
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top nghost+1 rows (u-fold) or nghost+2 rows
!  (T-fold) of physical domain for entire (global) top row
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call abort_ice( &
            'ice_HaloUpdate4DI4: Unknown field kind')
      end select

      if (halo%tripoleTFlag) then

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = -1
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do l=1,nt
           do k=1,nz
           do i = 2,nxGlobal/2
              iDst = nxGlobal - i + 2
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 0
           joffset = 1
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value

           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = -1
           joffset = 1
  
        case default
           call abort_ice( &
              'ice_HaloUpdate4DI4: Unknown field location')
        end select 
  
      else ! tripole u-fold  

        select case (fieldLoc)
        case (field_loc_center)   ! cell center location
  
           ioffset = 0
           joffset = 0
  
        case (field_loc_NEcorner)   ! cell corner location
  
           ioffset = 1
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2 - 1
              iDst = nxGlobal - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case (field_loc_Eface)   ! cell center location
  
           ioffset = 1
           joffset = 0
  
        case (field_loc_Nface)   ! cell corner (velocity) location
  
           ioffset = 0
           joffset = 1
  
           !*** top row is degenerate, so must enforce symmetry
           !***   use average of two degenerate points for value
  
           do l=1,nt
           do k=1,nz
           do i = 1,nxGlobal/2
              iDst = nxGlobal + 1 - i
              x1 = bufTripole(i   ,halo%tripoleRows,k,l)
              x2 = bufTripole(iDst,halo%tripoleRows,k,l)
              xavg = nint(0.5_dbl_kind*(x1 + isign*x2))
              bufTripole(i   ,halo%tripoleRows,k,l) = xavg
              bufTripole(iDst,halo%tripoleRows,k,l) = isign*xavg
           end do
           end do
           end do
  
        case default
           call abort_ice( &
              'ice_HaloUpdate4DI4: Unknown field location')
        end select
 
      endif

      !*** copy out of global tripole buffer into local
      !*** ghost cells

      !*** look through local copies to find the copy out
      !*** messages (srcBlock < 0)

      do nmsg=1,halo%numLocalCopies
         srcBlock = halo%srcLocalAddr(3,nmsg)

         if (srcBlock < 0) then

            iSrc     = halo%srcLocalAddr(1,nmsg) ! tripole buffer addr
            jSrc     = halo%srcLocalAddr(2,nmsg)

            iDst     = halo%dstLocalAddr(1,nmsg) ! local block addr
            jDst     = halo%dstLocalAddr(2,nmsg)
            dstBlock = halo%dstLocalAddr(3,nmsg)

            !*** correct for offsets
            iSrc = iSrc - ioffset
            jSrc = jSrc - joffset
            if (iSrc == 0) iSrc = nxGlobal
            if (iSrc > nxGlobal) iSrc = iSrc - nxGlobal

            !*** for center and Eface on u-fold, and NE corner and Nface
            !*** on T-fold, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= halo%tripoleRows .and. jSrc>0 .and. jDst>0) then
               do l=1,nt
               do k=1,nz
                  array(iDst,jDst,k,l,dstBlock) = isign*    &
                                  bufTripole(iSrc,jSrc,k,l)
               end do
               end do
            endif

         endif
      end do

   endif

   if (allocated(bufTripole)) deallocate(bufTripole)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate4DI4

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloIncrementMsgCount
! !INTERFACE:

   subroutine ice_HaloIncrementMsgCount(sndCounter, rcvCounter,    &
                                        srcProc, dstProc, msgSize)

! !DESCRIPTION:
!  This is a utility routine to increment the arrays for counting
!  whether messages are required.  It checks the source and destination
!  task to see whether the current task needs to send, receive or
!  copy messages to fill halo regions (ghost cells).

! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      srcProc,               &! source processor for communication
      dstProc,               &! destination processor for communication
      msgSize                 ! number of words for this message

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      sndCounter,       &! array for counting messages to be sent
      rcvCounter         ! array for counting messages to be received

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  error check
!
!-----------------------------------------------------------------------

   if (srcProc < 0 .or. dstProc < 0 .or. &
       srcProc > size(sndCounter)   .or. &
       dstProc > size(rcvCounter)) then
      call abort_ice( &
         'ice_HaloIncrementMsgCount: invalid processor number')
      return
   endif

!-----------------------------------------------------------------------
!
!  if destination all land or outside closed boundary (dstProc = 0), 
!  then no send is necessary, so do the rest only for dstProc /= 0
!
!-----------------------------------------------------------------------

   if (dstProc == 0) return

!-----------------------------------------------------------------------
!
!  if the current processor is the source, must send data 
!  local copy if dstProc = srcProc
!
!-----------------------------------------------------------------------

   if (srcProc == my_task + 1) sndCounter(dstProc) = &
                               sndCounter(dstProc) + msgSize

!-----------------------------------------------------------------------
!
!  if the current processor is the destination, must receive data 
!  local copy if dstProc = srcProc
!
!-----------------------------------------------------------------------

   if (dstProc == my_task + 1) then

      if (srcProc > 0) then  
         !*** the source block has ocean points
         !*** count as a receive from srcProc

         rcvCounter(srcProc) = rcvCounter(srcProc) + msgSize

      else
         !*** if the source block has been dropped, create
         !*** a local copy to fill halo with a fill value

         rcvCounter(dstProc) = rcvCounter(dstProc) + msgSize

      endif
   endif
!-----------------------------------------------------------------------
!EOC

   end subroutine ice_HaloIncrementMsgCount

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloMsgCreate
! !INTERFACE:

   subroutine ice_HaloMsgCreate(halo, dist, srcBlock, dstBlock, direction)

! !DESCRIPTION:
!  This is a utility routine to determine the required address and
!  message information for a particular pair of blocks.

! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      dist             ! distribution of blocks across procs

   integer (int_kind), intent(in) :: &
      srcBlock,   dstBlock   ! source,destination block id

   character (*), intent(in) :: &
      direction              ! direction of neighbor block
                             !  (north,south,east,west,
                             !   and NE, NW, SE, SW)

! !INPUT/OUTPUT PARAMETERS:

   type (ice_halo), intent(inout) :: &
      halo                   ! data structure containing halo info

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      srcProc, srcLocalID,   &! source block location in distribution
      dstProc, dstLocalID,   &! source block location in distribution
      msgIndx,               &! message counter and index into msg array
      ibSrc, ieSrc, jbSrc, jeSrc, &! phys domain info for source block
      ibDst, ieDst, jbDst, jeDst, &! phys domain info for dest   block
      nxGlobal,              &! size of global domain in e-w direction
      i,j,n                   ! dummy loop index

   integer (int_kind), dimension(:), pointer :: &
      iGlobal                 ! global i index for location in tripole

!-----------------------------------------------------------------------
!
!  initialize
!
!-----------------------------------------------------------------------

   if (allocated(bufTripoleR8)) nxGlobal = size(bufTripoleR8,dim=1)

!-----------------------------------------------------------------------
!
!  find source and destination block locations
!
!-----------------------------------------------------------------------

   if (srcBlock /= 0) then
      call ice_DistributionGetBlockLoc(dist, abs(srcBlock), srcProc, &
                                       srcLocalID)
   else
      srcProc    = 0
      srcLocalID = 0
   endif

   if (dstBlock /= 0) then
      call ice_DistributionGetBlockLoc(dist, abs(dstBlock), dstProc, &
                                       dstLocalID)
   else
      dstProc    = 0
      dstLocalID = 0
   endif

!-----------------------------------------------------------------------
!
!  if destination all land or outside closed boundary (dstProc = 0), 
!  then no send is necessary, so do the rest only for dstProc /= 0
!
!-----------------------------------------------------------------------

   if (dstProc == 0) return

!-----------------------------------------------------------------------
!
!  get block information if either block is local
!
!-----------------------------------------------------------------------

   if (srcProc == my_task+1 .or. dstProc == my_task+1) then

      if (srcBlock >= 0 .and. dstBlock >= 0) then
         call get_block_parameter(srcBlock, &
                                     ilo=ibSrc, ihi=ieSrc,   &
                                     jlo=jbSrc, jhi=jeSrc)
      else ! tripole - need iGlobal info
         call get_block_parameter(abs(srcBlock), &
                                     ilo=ibSrc, ihi=ieSrc,        &
                                     jlo=jbSrc, jhi=jeSrc,        &
                                     i_glob=iGlobal)

      endif

      if (dstBlock /= 0) then
         call get_block_parameter(abs(dstBlock), &
                                     ilo=ibDst, ihi=ieDst,   &
                                     jlo=jbDst, jhi=jeDst)
      endif

   endif

!-----------------------------------------------------------------------
!
!  if both blocks are local, create a local copy to fill halo
!
!-----------------------------------------------------------------------

   if (srcProc == my_task+1 .and. &
       dstProc == my_task+1) then   

      !*** compute addresses based on direction

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = ieSrc - nghost + i
            halo%srcLocalAddr(2,msgIndx) = jbSrc + j - 1
            halo%srcLocalAddr(3,msgIndx) = srcLocalID

            halo%dstLocalAddr(1,msgIndx) = i
            halo%dstLocalAddr(2,msgIndx) = jbDst + j - 1
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
            halo%srcLocalAddr(2,msgIndx) = jbSrc + j - 1
            halo%srcLocalAddr(3,msgIndx) = srcLocalID

            halo%dstLocalAddr(1,msgIndx) = ieDst + i
            halo%dstLocalAddr(2,msgIndx) = jbDst + j - 1
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('north')

         !*** copy northern physical domain of src
         !*** into southern halo of dst

         if (srcBlock > 0 .and. dstBlock > 0) then  ! normal north boundary

            msgIndx = halo%numLocalCopies

            do j=1,nghost
            do i=1,ieSrc-ibSrc+1

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc - nghost + j
               halo%srcLocalAddr(3,msgIndx) = srcLocalID

               halo%dstLocalAddr(1,msgIndx) = ibDst + i - 1
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         else if (srcBlock > 0 .and. dstBlock < 0) then

            !*** tripole grid - copy info into tripole buffer
            !*** copy physical domain of top halo+1 rows
            !*** into global buffer at src location

            !*** perform an error check to make sure the
            !*** block has enough points to perform a tripole
            !*** update

            if (jeSrc - jbSrc + 1 < halo%tripoleRows) then
               call abort_ice( &
               'ice_HaloMsgCreate: not enough points in block for tripole')
               return
            endif 

            msgIndx = halo%numLocalCopies

            do j=1,halo%tripoleRows
            do i=1,ieSrc-ibSrc+1

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc-halo%tripoleRows+j
               halo%srcLocalAddr(3,msgIndx) = srcLocalID

               halo%dstLocalAddr(1,msgIndx) = iGlobal(ibSrc + i - 1)
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = -dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         else if (srcBlock < 0 .and. dstBlock > 0) then

            !*** tripole grid - set up for copying out of 
            !*** tripole buffer into ghost cell domains
            !*** include e-w ghost cells

            msgIndx = halo%numLocalCopies

            do j=1,halo%tripoleRows
            do i=1,ieSrc+nghost

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = nxGlobal - iGlobal(i) + 1
               halo%srcLocalAddr(2,msgIndx) = nghost + 3 - j
               halo%srcLocalAddr(3,msgIndx) = -srcLocalID

               halo%dstLocalAddr(1,msgIndx) = i
               if (j.gt.nghost+1) then
                 halo%dstLocalAddr(2,msgIndx) = -1 ! never used
               else
                 halo%dstLocalAddr(2,msgIndx) = jeSrc + j - 1
               endif
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,nghost
         do i=1,ieSrc-ibSrc+1

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
            halo%srcLocalAddr(2,msgIndx) = jbSrc + j - 1
            halo%srcLocalAddr(3,msgIndx) = srcLocalID

            halo%dstLocalAddr(1,msgIndx) = ibDst + i - 1
            halo%dstLocalAddr(2,msgIndx) = jeDst + j
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('northeast')

         !*** normal northeast boundary - just copy NE corner
         !*** of physical domain into SW halo of NE nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

            do j=1,nghost
            do i=1,nghost

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ieSrc - nghost + i
               halo%srcLocalAddr(2,msgIndx) = jeSrc - nghost + j
               halo%srcLocalAddr(3,msgIndx) = srcLocalID

               halo%dstLocalAddr(1,msgIndx) = i
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         else

            !*** tripole grid - copy entire top halo+1 
            !*** rows into global buffer at src location

            msgIndx = halo%numLocalCopies

            do j=1,nghost+1
            do i=1,ieSrc-ibSrc+1

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc-1-nghost+j
               halo%srcLocalAddr(3,msgIndx) = srcLocalID

               halo%dstLocalAddr(1,msgIndx) = iGlobal(ibSrc + i - 1)
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = -dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         endif

      case ('northwest')

         !*** normal northeast boundary - just copy NW corner
         !*** of physical domain into SE halo of NW nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

            do j=1,nghost
            do i=1,nghost

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc - nghost + j
               halo%srcLocalAddr(3,msgIndx) = srcLocalID

               halo%dstLocalAddr(1,msgIndx) = ieDst + i
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         else

            !*** tripole grid - copy entire top halo+1 
            !*** rows into global buffer at src location

            msgIndx = halo%numLocalCopies

            do j=1,nghost+1
            do i=1,ieSrc-ibSrc+1

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc-1-nghost+j
               halo%srcLocalAddr(3,msgIndx) = srcLocalID

               halo%dstLocalAddr(1,msgIndx) = iGlobal(ibSrc + i - 1)
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = -dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,nghost
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = ieSrc - nghost + i
            halo%srcLocalAddr(2,msgIndx) = jbSrc + j - 1
            halo%srcLocalAddr(3,msgIndx) = srcLocalID

            halo%dstLocalAddr(1,msgIndx) = i
            halo%dstLocalAddr(2,msgIndx) = jeDst + j
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,nghost
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
            halo%srcLocalAddr(2,msgIndx) = jbSrc + j - 1
            halo%srcLocalAddr(3,msgIndx) = srcLocalID

            halo%dstLocalAddr(1,msgIndx) = ieDst + i
            halo%dstLocalAddr(2,msgIndx) = jeDst + j
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case default

         call abort_ice( &
            'ice_HaloMsgCreate: unknown direction local copy')
         return

      end select

!-----------------------------------------------------------------------
!
!  if dest block is local and source block does not exist, create a 
!  local copy to fill halo with a fill value
!
!-----------------------------------------------------------------------

   else if (srcProc == 0 .and. dstProc == my_task+1) then   

      !*** compute addresses based on direction

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = 0
            halo%srcLocalAddr(2,msgIndx) = 0
            halo%srcLocalAddr(3,msgIndx) = 0

            halo%dstLocalAddr(1,msgIndx) = i
            halo%dstLocalAddr(2,msgIndx) = jbDst + j - 1
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = 0
            halo%srcLocalAddr(2,msgIndx) = 0
            halo%srcLocalAddr(3,msgIndx) = 0

            halo%dstLocalAddr(1,msgIndx) = ieDst + i
            halo%dstLocalAddr(2,msgIndx) = jbDst + j - 1
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('north')

         !*** copy northern physical domain of src
         !*** into southern halo of dst

         if (dstBlock > 0) then  ! normal north boundary

            msgIndx = halo%numLocalCopies

            do j=1,nghost
            do i=1,ieSrc-ibSrc+1

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = 0
               halo%srcLocalAddr(2,msgIndx) = 0
               halo%srcLocalAddr(3,msgIndx) = 0

               halo%dstLocalAddr(1,msgIndx) = ibDst + i - 1
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,nghost
         do i=1,ieSrc-ibSrc+1

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = 0
            halo%srcLocalAddr(2,msgIndx) = 0
            halo%srcLocalAddr(3,msgIndx) = 0

            halo%dstLocalAddr(1,msgIndx) = ibDst + i - 1
            halo%dstLocalAddr(2,msgIndx) = jeDst + j
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('northeast')

         !*** normal northeast boundary - just copy NE corner
         !*** of physical domain into SW halo of NE nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

            do j=1,nghost
            do i=1,nghost

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = 0
               halo%srcLocalAddr(2,msgIndx) = 0
               halo%srcLocalAddr(3,msgIndx) = 0

               halo%dstLocalAddr(1,msgIndx) = i
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         endif

      case ('northwest')

         !*** normal northeast boundary - just copy NW corner
         !*** of physical domain into SE halo of NW nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

            do j=1,nghost
            do i=1,nghost

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = 0
               halo%srcLocalAddr(2,msgIndx) = 0
               halo%srcLocalAddr(3,msgIndx) = 0

               halo%dstLocalAddr(1,msgIndx) = ieDst + i
               halo%dstLocalAddr(2,msgIndx) = j
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

            halo%numLocalCopies = msgIndx

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,nghost
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = 0
            halo%srcLocalAddr(2,msgIndx) = 0
            halo%srcLocalAddr(3,msgIndx) = 0

            halo%dstLocalAddr(1,msgIndx) = i
            halo%dstLocalAddr(2,msgIndx) = jeDst + j
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

         msgIndx = halo%numLocalCopies

         do j=1,nghost
         do i=1,nghost

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = 0
            halo%srcLocalAddr(2,msgIndx) = 0
            halo%srcLocalAddr(3,msgIndx) = 0

            halo%dstLocalAddr(1,msgIndx) = ieDst + i
            halo%dstLocalAddr(2,msgIndx) = jeDst + j
            halo%dstLocalAddr(3,msgIndx) = dstLocalID

         end do
         end do

         halo%numLocalCopies = msgIndx

      case default

         call abort_ice( &
            'ice_HaloMsgCreate: unknown direction local copy')
         return

      end select

!-----------------------------------------------------------------------
!
!  if none of the cases above, no message info required for this
!  block pair
!
!-----------------------------------------------------------------------

   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine ice_HaloMsgCreate

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloExtrapolate
! !INTERFACE:

 subroutine ice_HaloExtrapolate2DR8(ARRAY,dist,ew_bndy_type,ns_bndy_type)

! !DESCRIPTION:
!  This subroutine extrapolates ARRAY values into the first row or column 
!  of ghost cells, and is intended for grid variables whose ghost cells 
!  would otherwise be set using the default boundary conditions (Dirichlet 
!  or Neumann).
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface ice_HaloExtrapolate

! !USES:

   use ice_blocks
   use ice_constants
   use ice_distribution

! !INPUT PARAMETERS:

    character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

! !OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc,                &! ghost cell column or row
     npad                 ! padding column/row counter

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,iblk) = c2*ARRAY(2,j,iblk) - ARRAY(3,j,iblk)
            enddo
         endif

      elseif (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, - 1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,iblk) = c2*ARRAY(ibc-1,j,iblk) - ARRAY(ibc-2,j,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,iblk) = c2*ARRAY(i,2,iblk) - ARRAY(i,3,iblk)
            enddo
         endif

      elseif (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, - 1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,iblk) = c2*ARRAY(i,ibc-1,iblk) - ARRAY(i,ibc-2,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloExtrapolate2DR8

!***********************************************************************

end module ice_boundary

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
