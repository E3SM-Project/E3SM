!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_boundary

 module ice_boundary

! !DESCRIPTION:
!  This module contains data types and routines for updating halo
!  regions (ghost cells) using MPI calls
!
! !REVISION HISTORY:
!  SVN:$Id: ice_boundary.F90 231 2009-09-01 23:16:05Z eclare $
!  2007-07-19: Phil Jones, Yoshi Yoshida, John Dennis
!              new naming conventions, optimizations during
!              initialization, true multi-dimensional updates 
!              (rather than serial call to two-dimensional updates), 
!              fixes for non-existent blocks
!  2008-01-28: Elizabeth Hunke replaced old routines with new POP
!              infrastructure

! !USES:

   use ice_kinds_mod
   use ice_communicate, only: my_task, mpiR4, mpiR8, mpitagHalo
   use ice_domain_size, only: max_blocks
   use ice_constants, only: field_type_scalar, &
          field_type_vector, field_type_angle, &
         field_loc_center,  field_loc_NEcorner, &
         field_loc_Nface, field_loc_Eface
   use ice_global_reductions, only: global_maxval
   use ice_fileunits, only: nu_diag, flush_fileunit
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
   include 'mpif.h'
   save

! !PUBLIC TYPES:

   type, public :: ice_halo
      integer (int_kind) ::  &
         communicator,     &! communicator to use for update messages
         numMsgSend,       &! number of messages to send halo update
         numMsgRecv,       &! number of messages to recv halo update
         numLocalCopies,   &! num local copies for halo update
         tripoleRows        ! number of rows in tripole buffer

      logical (log_kind) ::  &
         tripoleTFlag       ! NS boundary is a tripole T-fold

      integer (int_kind), dimension(:), pointer :: &
         recvTask,         &! task from which to recv each msg
         sendTask,         &! task to   which to send each msg
         sizeSend,         &! size of each sent message
         sizeRecv,         &! size of each recvd message
         tripSend,         &! send msg tripole flag, 0=non-zipper block
         tripRecv           ! recv msg tripole flag, for masked halos

      integer (int_kind), dimension(:,:), pointer :: &
         srcLocalAddr,     &! src addresses for each local copy
         dstLocalAddr       ! dst addresses for each local copy

      integer (int_kind), dimension(:,:,:), pointer :: &
         sendAddr,         &! src addresses for each sent message
         recvAddr           ! dst addresses for each recvd message

   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: ice_HaloCreate,  &
             ice_HaloMask,  &
             ice_HaloDestroy, &
             ice_HaloUpdate,  &
             ice_HaloUpdate_stress,  &
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
!  to prevent frequent allocate-deallocate for 2d halo updates, create
!  a static 2d buffer to be allocated once at creation.  if future 
!  creation needs larger buffer, resize during the creation.
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufSizeSend,    &! max buffer size for send messages
      bufSizeRecv      ! max buffer size for recv messages

   integer (int_kind), dimension(:,:), allocatable :: &
      bufSendI4,     &! buffer for use to send in 2D i4 halo updates
      bufRecvI4       ! buffer for use to recv in 2D i4 halo updates

   real (real_kind), dimension(:,:), allocatable :: &
      bufSendR4,     &! buffer for use to send in 2D r4 halo updates
      bufRecvR4       ! buffer for use to recv in 2D r4 halo updates

   real (dbl_kind), dimension(:,:), allocatable :: &
      bufSendR8,     &! buffer for use to send in 2D r8 halo updates
      bufRecvR8       ! buffer for use to recv in 2D r8 halo updates

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
      maxTmp,                      &! temp for global maxval      
      blockSizeX,                  &! size of default physical domain in X 
      blockSizeY,                  &! size of default physical domain in Y 
      maxSizeSend, maxSizeRecv,    &! max buffer sizes
      numMsgSend, numMsgRecv,      &! number of messages for this halo
      eastMsgSize, westMsgSize,    &! nominal sizes for e-w msgs
      northMsgSize, southMsgSize,  &! nominal sizes for n-s msgs
      tripoleMsgSize,              &! size for tripole messages
      tripoleMsgSizeOut,           &! size for tripole messages
      tripoleRows,                 &! number of rows in tripole buffer
      cornerMsgSize, msgSize        ! nominal size for corner msg

   integer (int_kind), dimension(:), allocatable :: &
      sendCount, recvCount          ! count number of words to each proc

   logical (log_kind) :: &
      resize,               &! flag for resizing buffers
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
   northMsgSize = nghost*blockSizeX
   cornerMsgSize = nghost*nghost
   tripoleRows = nghost+1

   if (nsBoundaryType == 'tripole' .or. nsBoundaryType == 'tripoleT') then
      tripoleFlag = .true.
      tripoleTFlag = (nsBoundaryType == 'tripoleT')
      if (tripoleTflag) tripoleRows = tripoleRows+1

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
   endif
   halo%tripoleTFlag = tripoleTFlag
   halo%tripoleRows = tripoleRows
   tripoleMsgSize = tripoleRows*blockSizeX
   tripoleMsgSizeOut = tripoleRows*nx_block

!-----------------------------------------------------------------------
!
!  Count the number of messages to send/recv from each processor
!  and number of words in each message.  These quantities are
!  necessary for allocating future arrays.
!
!-----------------------------------------------------------------------

   allocate (sendCount(numProcs), recvCount(numProcs), stat=istat)

   if (istat > 0) then
      call abort_ice( &
                        'ice_HaloCreate: error allocating count arrays')
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
         msgSize = northMsgSize
         call ice_distributionGetBlockLoc(dist, northBlock, dstProc, &
                                          dstLocalID)
      else if (northBlock < 0) then ! tripole north row, count block
         tripoleBlock = .true.
         msgSize = tripoleMsgSize
         call ice_distributionGetBlockLoc(dist, abs(northBlock), &
                                 dstProc, dstLocalID)
      else
         tripoleBlock = .false.
         msgSize = northMsgSize
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloIncrementMsgCount(sendCount, recvCount,           &
                                     srcProc, dstProc, msgSize)

      !*** if a tripole boundary block, also create a local
      !*** message into and out of tripole buffer 

      if (tripoleBlock) then
         !*** copy out of tripole buffer - includes halo
         call ice_HaloIncrementMsgCount(sendCount, recvCount,        &
                                        srcProc, srcProc,            &
                                        tripoleMsgSizeOut)

         !*** copy in only required if dstProc not same as srcProc
         if (dstProc /= srcProc) then
            call ice_HaloIncrementMsgCount(sendCount, recvCount,  &
                                           srcProc, srcProc,      & 
                                           msgSize)
         endif
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

      if (tripoleBlock .and. dstProc /= srcProc) then
         call ice_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, tripoleMsgSize)
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

      if (tripoleBlock .and. dstProc /= srcProc) then
         call ice_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, tripoleMsgSize)
      endif

      !*** find northeast neighbor block and add to message count

      neBlock = ice_blocksGetNbrID(iblock, ice_blocksNorthEast,    &
                                   ewBoundaryType, nsBoundaryType)

      if (neBlock > 0) then
         msgSize = cornerMsgSize  ! normal corner message 

         call ice_distributionGetBlockLoc(dist, neBlock, dstProc, &
                                          dstLocalID)

      else if (neBlock < 0) then ! tripole north row
         msgSize = tripoleMsgSize  ! tripole needs whole top row of block

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
         msgSize = tripoleMsgSize ! tripole NE corner update - entire row needed

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

      !*** for tripole grids with padded domain, padding will
      !*** prevent tripole buffer from getting all the info
      !*** it needs - must extend footprint at top boundary

      if (tripoleBlock                  .and. & !tripole
          mod(nxGlobal,blockSizeX) /= 0) then   !padding

         !*** find east2 neighbor block and add to message count

         eastBlock = ice_blocksGetNbrID(iBlock, ice_blocksEast2,     &
                                     ewBoundaryType, nsBoundaryType)

         if (eastBlock > 0) then
            call ice_distributionGetBlockLoc(dist, eastBlock, dstProc, &
                                             dstLocalID)
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloIncrementMsgCount(sendCount, recvCount,       &
                                     srcProc, dstProc, tripoleMsgSize)
         endif

         !*** find EastNorthEast neighbor block and add to message count

         neBlock = ice_blocksGetNbrID(iBlock, ice_blocksEastNorthEast, &
                                     ewBoundaryType, nsBoundaryType)

         if (neBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call ice_distributionGetBlockLoc(dist, abs(neBlock), dstProc, &
                                             dstLocalID)
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloIncrementMsgCount(sendCount, recvCount,   &
                                        srcProc, dstProc, msgSize)
         endif

         !*** find west2 neighbor block and add to message count

         westBlock = ice_blocksGetNbrID(iBlock, ice_blocksWest2,     &
                                     ewBoundaryType, nsBoundaryType)

         if (westBlock > 0) then
            call ice_distributionGetBlockLoc(dist, westBlock, dstProc, &
                                             dstLocalID)
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloIncrementMsgCount(sendCount, recvCount,       &
                                     srcProc, dstProc, tripoleMsgSize)
         endif

         !*** find WestNorthWest neighbor block and add to message count

         nwBlock = ice_blocksGetNbrID(iBlock, ice_blocksWestNorthWest, &
                                     ewBoundaryType, nsBoundaryType)

         if (nwBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call ice_distributionGetBlockLoc(dist, abs(nwBlock), dstProc, &
                                             dstLocalID)
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloIncrementMsgCount(sendCount, recvCount,   &
                                        srcProc, dstProc, msgSize)
         endif

      endif

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
!  now count the number of actual messages to be sent and received
!
!-----------------------------------------------------------------------

   numMsgSend = count(sendCount /= 0)
   numMsgRecv = count(recvCount /= 0)
   halo%numMsgSend = numMsgSend
   halo%numMsgRecv = numMsgRecv

!-----------------------------------------------------------------------
!
!  allocate buffers for 2-d halo updates to save time later
!  if the buffers are already allocated by previous create call,
!   check to see if they need to be re-sized
!
!-----------------------------------------------------------------------
   
   maxTmp = maxval(sendCount)
   maxSizeSend = global_maxval(maxTmp, dist)
   maxTmp = maxval(recvCount)
   maxSizeRecv = global_maxval(maxTmp, dist)

   if (.not. allocated(bufSendR8)) then

      bufSizeSend = maxSizeSend
      bufSizeRecv = maxSizeRecv

      allocate(bufSendI4(bufSizeSend, numMsgSend), &
               bufRecvI4(bufSizeRecv, numMsgRecv), &
               bufSendR4(bufSizeSend, numMsgSend), &
               bufRecvR4(bufSizeRecv, numMsgRecv), &
               bufSendR8(bufSizeSend, numMsgSend), &
               bufRecvR8(bufSizeRecv, numMsgRecv), stat=istat)

      if (istat > 0) then
         call abort_ice( &
            'ice_HaloCreate: error allocating 2d buffers')
         return
      endif

   else

      resize = .false.

      if (maxSizeSend > bufSizeSend) then
         resize = .true.
         bufSizeSend = maxSizeSend
      endif
      if (maxSizeRecv > bufSizeRecv) then
         resize = .true.
         bufSizeRecv = maxSizeRecv
      endif

      if (numMsgSend > size(bufSendR8,dim=2)) resize = .true.
      if (numMsgRecv > size(bufRecvR8,dim=2)) resize = .true.

      if (resize) then
         deallocate(bufSendI4, bufRecvI4, bufSendR4, &
                    bufRecvR4, bufSendR8, bufRecvR8, stat=istat)

         if (istat > 0) then
            call abort_ice( &
               'ice_HaloCreate: error deallocating 2d buffers')
            return
         endif

         allocate(bufSendI4(bufSizeSend, numMsgSend), &
                  bufRecvI4(bufSizeRecv, numMsgRecv), &
                  bufSendR4(bufSizeSend, numMsgSend), &
                  bufRecvR4(bufSizeRecv, numMsgRecv), &
                  bufSendR8(bufSizeSend, numMsgSend), &
                  bufRecvR8(bufSizeRecv, numMsgRecv), stat=istat)

         if (istat > 0) then
            call abort_ice( &
               'ice_HaloCreate: error reallocating 2d buffers')
            return
         endif

      endif

   endif

!-----------------------------------------------------------------------
!
!  allocate arrays for message information and initialize
!
!-----------------------------------------------------------------------

   allocate(halo%sendTask(numMsgSend), &
            halo%recvTask(numMsgRecv), &
            halo%sizeSend(numMsgSend), &
            halo%sizeRecv(numMsgRecv), &
            halo%tripSend(numMsgSend), &
            halo%tripRecv(numMsgRecv), &
            halo%sendAddr(3,bufSizeSend,numMsgSend), &
            halo%recvAddr(3,bufSizeRecv,numMsgRecv), &
            halo%srcLocalAddr(3,halo%numLocalCopies), &
            halo%dstLocalAddr(3,halo%numLocalCopies), &
            stat = istat)

   if (istat > 0) then
      call abort_ice( &
         'ice_HaloCreate: error allocating halo message info arrays')
      return
   endif

   halo%sendTask = 0
   halo%recvTask = 0
   halo%sizeSend = 0
   halo%sizeRecv = 0
   halo%tripSend = 0
   halo%tripRecv = 0
   halo%sendAddr = 0
   halo%recvAddr = 0
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

   halo%numMsgSend     = 0
   halo%numMsgRecv     = 0
   halo%numLocalCopies = 0

   msgConfigLoop: do iblock=1,nblocks_tot

      call ice_distributionGetBlockLoc(dist, iblock, srcProc, &
                                       srcLocalID)

      !*** find north neighbor block and set msg info
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

      call ice_HaloMsgCreate(halo, iblock,     srcProc, srcLocalID, &
                                   northBlock, dstProc, dstLocalID, &
                                   'north')

      !*** if a tripole boundary block, also create a local
      !*** message into and out of tripole buffer 

      if (tripoleBlock) then
         !*** copy out of tripole buffer - includes halo
         call ice_HaloMsgCreate(halo,-iblock, srcProc, srcLocalID, &
                                      iblock, srcProc, srcLocalID, &
                                      'north')

         !*** copy in only required if dstProc not same as srcProc
         if (dstProc /= srcProc) then
            call ice_HaloMsgCreate(halo, iblock, srcProc, srcLocalID, &
                                        -iblock, srcProc, srcLocalID, &
                                         'north')

         endif
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

      call ice_HaloMsgCreate(halo, iblock,     srcProc, srcLocalID, &
                                   southBlock, dstProc, dstLocalID, &
                                   'south')

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

      call ice_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                   eastBlock, dstProc, dstLocalID, &
                                   'east')

      !*** if a tripole boundary block, non-local east neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

      if (tripoleBlock .and. dstProc /= srcProc) then
         call ice_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                     -eastBlock, dstProc, dstLocalID, &
                                      'north')

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

      call ice_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                   westBlock, dstProc, dstLocalID, &
                                   'west')


      !*** if a tripole boundary block, non-local west neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

      if (tripoleBlock .and. dstProc /= srcProc) then
         call ice_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                     -westBlock, dstProc, dstLocalID, &
                                      'north')

      endif

      !*** find northeast neighbor block and add to message count

      neBlock = ice_blocksGetNbrID(iblock, ice_blocksNorthEast,    &
                                   ewBoundaryType, nsBoundaryType)

      if (neBlock /= 0) then
         call ice_distributionGetBlockLoc(dist, abs(neBlock), dstProc, &
                                          dstLocalID)

      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   neBlock, dstProc, dstLocalID, &
                                   'northeast')

      !*** find northwest neighbor block and add to message count

      nwBlock = ice_blocksGetNbrID(iblock, ice_blocksNorthWest,    &
                                   ewBoundaryType, nsBoundaryType)

      if (nwBlock /= 0) then
         call ice_distributionGetBlockLoc(dist, abs(nwBlock), dstProc, &
                                          dstLocalID)

      else
         dstProc = 0
         dstLocalID = 0
      endif

      call ice_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   nwBlock, dstProc, dstLocalID, &
                                   'northwest')

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

      call ice_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   seBlock, dstProc, dstLocalID, &
                                   'southeast')

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

      call ice_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   swBlock, dstProc, dstLocalID, &
                                   'southwest')

      !*** for tripole grids with padded domain, padding will
      !*** prevent tripole buffer from getting all the info
      !*** it needs - must extend footprint at top boundary

      if (tripoleBlock                  .and. & !tripole
          mod(nxGlobal,blockSizeX) /= 0) then   !padding

         !*** find east2 neighbor block and add to message count

         eastBlock = ice_blocksGetNbrID(iBlock, ice_blocksEast2,     &
                                     ewBoundaryType, nsBoundaryType)

         if (eastBlock > 0) then
            call ice_distributionGetBlockLoc(dist, eastBlock, dstProc, &
                                             dstLocalID)

         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                        -eastBlock, dstProc, dstLocalID, &
                                         'north')

         endif

         !*** find EastNorthEast neighbor block and add to message count

         neBlock = ice_blocksGetNbrID(iBlock, ice_blocksEastNorthEast, &
                                     ewBoundaryType, nsBoundaryType)

         if (neBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call ice_distributionGetBlockLoc(dist, abs(neBlock), dstProc, &
                                             dstLocalID)

         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                         neBlock, dstProc, dstLocalID, &
                                         'north')
         endif

         !*** find west2 neighbor block and add to message count

         westBlock = ice_blocksGetNbrID(iBlock, ice_blocksWest2,     &
                                     ewBoundaryType, nsBoundaryType)

         if (westBlock > 0) then
            call ice_distributionGetBlockLoc(dist, westBlock, dstProc, &
                                             dstLocalID)

         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                        -westBlock, dstProc, dstLocalID, &
                                         'north')

         endif

         !*** find WestNorthWest neighbor block and add to message count

         nwBlock = ice_blocksGetNbrID(iBlock, ice_blocksWestNorthWest, &
                                     ewBoundaryType, nsBoundaryType)

         if (nwBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call ice_distributionGetBlockLoc(dist, abs(nwBlock), dstProc, &
                                             dstLocalID)

         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call ice_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                         nwBlock, dstProc, dstLocalID, &
                                         'north')

         endif

      endif

   end do msgConfigLoop

!   write(nu_diag,'(a,5i8)') 'ice_HaloCreate ',my_task,halo%numMsgSend,halo%numMsgRecv, &
!      sum(halo%sizeSend(1:halo%numMsgSend)),sum(halo%sizeRecv(1:halo%numMsgRecv))

!-----------------------------------------------------------------------
!EOC

 end function ice_HaloCreate

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloMask
! !INTERFACE:

 subroutine ice_HaloMask(halo, basehalo, mask)

! !DESCRIPTION:
!  This routine creates a halo type with info necessary for
!  performing a halo (ghost cell) update. This info is computed
!  based on a base halo already initialized and a mask
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (ice_halo) :: &
      basehalo            ! basehalo to mask
   integer (int_kind), intent(in) ::  &
      mask(nx_block,ny_block,max_blocks)   ! mask of live points

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

   integer (int_kind) ::           &
      n,nmsg,scnt,                 &! counters
      icel,jcel,nblock,            &! gridcell index
      istat,                       &! allocate status flag
      communicator,                &! communicator for message passing
      numMsgSend, numMsgRecv,      &! number of messages for this halo
      numLocalCopies,              &! num local copies for halo update
      tripoleRows,                 &! number of rows in tripole buffer
      lbufSizeSend,                &! buffer size for send messages
      lbufSizeRecv                  ! buffer size for recv messages
   logical (log_kind) :: &
      tripoleTFlag,                &! flag for processing tripole buffer as T-fold
      tmpflag                       ! temporary flag for if tests

!-----------------------------------------------------------------------
!
!  allocate and initialize halo
!  always keep tripole zipper msgs
!
!-----------------------------------------------------------------------

      communicator   = basehalo%communicator
      tripoleRows    = basehalo%tripoleRows
      tripoleTFlag   = basehalo%tripoleTFlag
      numMsgSend     = basehalo%numMsgSend
      numMsgRecv     = basehalo%numMsgRecv
      numLocalCopies = basehalo%numLocalCopies
      lbufSizeSend   = size(basehalo%sendAddr,dim=2)
      lbufSizeRecv   = size(basehalo%recvAddr,dim=2)

      allocate(halo%sendTask(numMsgSend), &
               halo%recvTask(numMsgRecv), &
               halo%sizeSend(numMsgSend), &
               halo%sizeRecv(numMsgRecv), &
               halo%tripSend(numMsgSend), &
               halo%tripRecv(numMsgRecv), &
               halo%sendAddr(3,lbufSizeSend,numMsgSend), &
               halo%recvAddr(3,lbufSizeRecv,numMsgRecv), &
               halo%srcLocalAddr(3,numLocalCopies), &
               halo%dstLocalAddr(3,numLocalCopies), &
               stat = istat)

      if (istat > 0) then
         call abort_ice( &
            'ice_HaloMask: error allocating halo message info arrays')
         return
      endif

      halo%communicator   = communicator
      halo%tripoleRows    = tripoleRows
      halo%tripoleTFlag   = tripoleTFlag
!      halo%numMsgSend     = numMsgSend
!      halo%numMsgRecv     = numMsgRecv
      halo%numLocalCopies = numLocalCopies

!      halo%recvTask       = basehalo%recvTask
!      halo%sendTask       = basehalo%sendTask
!      halo%sizeSend       = basehalo%sizeSend
!      halo%sizeRecv       = basehalo%sizeRecv
!      halo%tripSend       = basehalo%tripSend
!      halo%tripRecv       = basehalo%tripRecv
      halo%srcLocalAddr   = basehalo%srcLocalAddr
      halo%dstLocalAddr   = basehalo%dstLocalAddr
!      halo%sendAddr       = basehalo%sendAddr
!      halo%recvAddr       = basehalo%recvAddr

   numMsgSend = 0
   do nmsg=1,basehalo%numMsgSend
      scnt = 0
      do n=1,basehalo%sizeSend(nmsg)
         icel     = basehalo%sendAddr(1,n,nmsg)
         jcel     = basehalo%sendAddr(2,n,nmsg)
         nblock   = basehalo%sendAddr(3,n,nmsg)

! the following line fails bounds check for mask when tripSend /= 0
!         if (mask(icel,jcel,abs(nblock)) /= 0 .or. basehalo%tripSend(nmsg) /= 0) then
         tmpflag = .false.
         if (basehalo%tripSend(nmsg) /= 0) then
            tmpflag = .true.
         elseif (mask(icel,jcel,abs(nblock)) /= 0) then
            tmpflag = .true.
         endif
         
         if (tmpflag) then
            scnt = scnt + 1
            if (scnt == 1) then
               numMsgSend = numMsgSend + 1
               halo%sendTask(numMsgSend) = basehalo%sendTask(nmsg)
               halo%tripSend(numMsgSend) = basehalo%tripSend(nmsg)
            endif
            halo%sendAddr(1,scnt,numMsgSend) = icel
            halo%sendAddr(2,scnt,numMsgSend) = jcel
            halo%sendAddr(3,scnt,numMsgSend) = nblock
            halo%sizeSend(numMsgSend) = scnt
         endif
      enddo
   enddo
   halo%numMsgSend = numMsgSend     

   numMsgRecv = 0
   do nmsg=1,basehalo%numMsgRecv
      scnt = 0
      do n=1,basehalo%sizeRecv(nmsg)
         icel     = basehalo%recvAddr(1,n,nmsg)
         jcel     = basehalo%recvAddr(2,n,nmsg)
         nblock   = basehalo%recvAddr(3,n,nmsg)

! the following line fails bounds check for mask when tripSend /= 0
!         if (mask(icel,jcel,abs(nblock)) /= 0 .or. basehalo%tripRecv(nmsg) /= 0) then
         tmpflag = .false.
         if (basehalo%tripRecv(nmsg) /= 0) then
            tmpflag = .true.
         elseif (mask(icel,jcel,abs(nblock)) /= 0) then
            tmpflag = .true.
         endif

         if (tmpflag) then
            scnt = scnt + 1
            if (scnt == 1) then
               numMsgRecv = numMsgRecv + 1
               halo%recvTask(numMsgRecv) = basehalo%recvTask(nmsg)
               halo%tripRecv(numMsgRecv) = basehalo%tripRecv(nmsg)
            endif
            halo%recvAddr(1,scnt,numMsgRecv) = icel
            halo%recvAddr(2,scnt,numMsgRecv) = jcel
            halo%recvAddr(3,scnt,numMsgRecv) = nblock
            halo%sizeRecv(numMsgRecv) = scnt
         endif
      enddo
   enddo
   halo%numMsgRecv = numMsgRecv     

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloMask

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloDestroy
! !INTERFACE:

 subroutine ice_HaloDestroy(halo)

! !DESCRIPTION:
!  This routine creates a halo type with info necessary for
!  performing a halo (ghost cell) update. This info is computed
!  based on the input block distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   type (ice_halo) :: &
      halo               ! a new halo type with info for halo updates

!EOP
!BOC
   integer (int_kind) ::           &
      istat                      ! error or status flag for MPI,alloc
!-----------------------------------------------------------------------

   deallocate(halo%sendTask, stat=istat)
   deallocate(halo%recvTask, stat=istat)
   deallocate(halo%sizeSend, stat=istat)
   deallocate(halo%sizeRecv, stat=istat)
   deallocate(halo%tripSend, stat=istat)
   deallocate(halo%tripRecv, stat=istat)
   deallocate(halo%srcLocalAddr, stat=istat)
   deallocate(halo%dstLocalAddr, stat=istat)
   deallocate(halo%sendAddr, stat=istat)
   deallocate(halo%recvAddr, stat=istat)

end subroutine ice_HaloDestroy
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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   real (dbl_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (int_kind) ::  len  ! length of messages

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
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate2DR8: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)
      call MPI_IRECV(bufRecvR8(1:len,nmsg), len, mpiR8, &
                     halo%recvTask(nmsg),               &
                     mpitagHalo + halo%recvTask(nmsg),  &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         bufSendR8(n,nmsg) = array(iSrc,jSrc,srcBlock)
      end do
      do n=halo%sizeSend(nmsg)+1,bufSizeSend
         bufSendR8(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)
      call MPI_ISEND(bufSendR8(1:len,nmsg), len, mpiR8, &
                     halo%sendTask(nmsg),               &
                     mpitagHalo + my_task,              &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do


!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:) = fill
      array(1:nx_block,ny_block-j+1,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:) = fill
      array(nx_block-i+1,1:ny_block,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            array(iDst,jDst,dstBlock) = bufRecvR8(n,nmsg)
         else if (dstBlock < 0) then !tripole
            bufTripoleR8(iDst,jDst) = bufRecvR8(n,nmsg)
         endif
      end do
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
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate2DR8: error deallocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_HaloUpdate2DR8

!***********************************************************************
!BOP
! !IROUTINE: ice_HaloUpdate_stress
! !INTERFACE:

 subroutine ice_HaloUpdate_stress(array1, array2, halo, &
                               fieldLoc, fieldKind,     &
                               fillValue)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array using
!  a second array as needed by the stress fields.
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
      array1           ,&  ! array containing field for which halo
                           ! needs to be updated
      array2               ! array containing field for which halo
                           ! in array1 needs to be updated

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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   real (dbl_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (int_kind) ::  len  ! length of messages

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
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate_stress: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)
      call MPI_IRECV(bufRecvR8(1:len,nmsg), len, mpiR8, &
                     halo%recvTask(nmsg),               &
                     mpitagHalo + halo%recvTask(nmsg),  &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         bufSendR8(n,nmsg) = array2(iSrc,jSrc,srcBlock)
      end do
      do n=halo%sizeSend(nmsg)+1,bufSizeSend
         bufSendR8(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)
      call MPI_ISEND(bufSendR8(1:len,nmsg), len, mpiR8, &
                     halo%sendTask(nmsg),               &
                     mpitagHalo + my_task,              &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  do NOT zero the halo out, this halo update just updates
!  the tripole zipper as needed for stresses.  if you zero
!  it out, all halo values will be wiped out.
!-----------------------------------------------------------------------
!   do j = 1,nghost
!      array1(1:nx_block,           j,:) = fill
!      array1(1:nx_block,ny_block-j+1,:) = fill
!   enddo
!   do i = 1,nghost
!      array1(i,           1:ny_block,:) = fill
!      array1(nx_block-i+1,1:ny_block,:) = fill
!   enddo

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
         if (dstBlock < 0) then ! tripole copy into buffer
            bufTripoleR8(iDst,jDst) = &
            array2(iSrc,jSrc,srcBlock)
         endif
      else if (srcBlock == 0) then
         array1(iDst,jDst,dstBlock) = fill
     endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock < 0) then !tripole
            bufTripoleR8(iDst,jDst) = bufRecvR8(n,nmsg)
         endif
      end do
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
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
            'ice_HaloUpdate_stress: Unknown field kind')
      end select

      select case (fieldLoc)
      case (field_loc_center)   ! cell center location

         ioffset = 0
         joffset = 0

      case (field_loc_NEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

      case (field_loc_Eface) 

         ioffset = 1
         joffset = 0

      case (field_loc_Nface) 

         ioffset = 0
         joffset = 1

      case default
         call abort_ice( &
               'ice_HaloUpdate_stress: Unknown field location')
      end select

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
 
            !*** for center and Eface, do not need to replace
            !*** top row of physical domain, so jSrc should be
            !*** out of range and skipped
            !*** otherwise do the copy

            if (jSrc <= nghost+1) then
               array1(iDst,jDst,dstBlock) = isign*bufTripoleR8(iSrc,jSrc)
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate_stress: error deallocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!EOC
 end subroutine ice_HaloUpdate_stress

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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   real (real_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

    integer (int_kind) :: len  ! length of messages

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
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate2DR4: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)
      call MPI_IRECV(bufRecvR4(1:len,nmsg), len, mpiR4, &
                     halo%recvTask(nmsg),               &
                     mpitagHalo + halo%recvTask(nmsg),  &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         bufSendR4(n,nmsg) = array(iSrc,jSrc,srcBlock)
      end do
      do n=halo%sizeSend(nmsg)+1,bufSizeSend
         bufSendR4(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)
      call MPI_ISEND(bufSendR4(1:len,nmsg), len, mpiR4, &
                     halo%sendTask(nmsg),               &
                     mpitagHalo + my_task,              &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:) = fill
      array(1:nx_block,ny_block-j+1,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:) = fill
      array(nx_block-i+1,1:ny_block,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            array(iDst,jDst,dstBlock) = bufRecvR4(n,nmsg)
         else if (dstBlock < 0) then !tripole
            bufTripoleR4(iDst,jDst) = bufRecvR4(n,nmsg)
         endif
      end do
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
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate2DR4: error deallocating req,status arrays')
      return
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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   integer (int_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (int_kind) :: len ! length of messages

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
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate2DI4: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)
      call MPI_IRECV(bufRecvI4(1:len,nmsg), len, MPI_INTEGER, &
                     halo%recvTask(nmsg),                     &
                     mpitagHalo + halo%recvTask(nmsg),        &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         bufSendI4(n,nmsg) = array(iSrc,jSrc,srcBlock)
      end do
      do n=halo%sizeSend(nmsg)+1,bufSizeSend
         bufSendI4(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)
      call MPI_ISEND(bufSendI4(1:len,nmsg), len, MPI_INTEGER, &
                     halo%sendTask(nmsg),                     &
                     mpitagHalo + my_task,                    &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:) = fill
      array(1:nx_block,ny_block-j+1,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:) = fill
      array(nx_block-i+1,1:ny_block,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            array(iDst,jDst,dstBlock) = bufRecvI4(n,nmsg)
         else if (dstBlock < 0) then !tripole
            bufTripoleI4(iDst,jDst) = bufRecvI4(n,nmsg)
         endif
      end do
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
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate2DI4: error deallocating req,status arrays')
      return
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
                               fillValue, mode)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  ice\_HaloUpdate.  This routine is the specific interface
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

   character(*), intent(in), optional :: &
      mode                 ! optional value to specify mode of comm

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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   real (dbl_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (dbl_kind), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 3d send,recv buffers

   real (dbl_kind), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

   integer (int_kind) :: len ! length of message 
   logical :: do_send, do_recv
   save

!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   do_send = .true.
   do_recv = .true.
   if (present(mode)) then
      if (trim(mode) == 'send') then
         do_send=.true.
         do_recv=.false.
      elseif (trim(mode) == 'recv') then
         do_send=.false.
         do_recv=.true.
      else
         call abort_ice( &
         'ice_HaloUpdate3DR8: illegal mode : '//trim(mode))
      endif
   endif

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_dbl_kind
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR8)) nxGlobal = size(bufTripoleR8,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

 if (do_send) then

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR8: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  allocate 3D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)

   allocate(bufSend(bufSizeSend*nz, halo%numMsgSend), &
            bufRecv(bufSizeRecv*nz, halo%numMsgRecv), &
            bufTripole(nxGlobal, halo%tripoleRows, nz), &
            stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR8: error allocating buffers')
      return
   endif

   bufTripole = fill

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)*nz
      call MPI_IRECV(bufRecv(1:len,nmsg), len, mpiR8,   &
                     halo%recvTask(nmsg),               &
                     mpitagHalo + halo%recvTask(nmsg),  &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      i=0
      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         do k=1,nz
            i = i + 1
            bufSend(i,nmsg) = array(iSrc,jSrc,k,srcBlock)
         end do
      end do
      do n=i+1,bufSizeSend*nz
         bufSend(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)*nz
      call MPI_ISEND(bufSend(1:len,nmsg), len, mpiR8, &
                     halo%sendTask(nmsg),             &
                     mpitagHalo + my_task,            &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:,:) = fill
      array(1:nx_block,ny_block-j+1,:,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:,:) = fill
      array(nx_block-i+1,1:ny_block,:,:) = fill
   enddo

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
 endif   ! do_send

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

 if (do_recv) then
   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      i = 0
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            do k=1,nz
               i = i + 1
               array(iDst,jDst,k,dstBlock) = bufRecv(i,nmsg)
            end do
         else if (dstBlock < 0) then !tripole
            do k=1,nz
               i = i + 1
               bufTripole(iDst,jDst,k) = bufRecv(i,nmsg)
            end do
         endif
      end do
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
!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR8: error deallocating req,status arrays')
      return
   endif

   deallocate(bufSend, bufRecv, bufTripole, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR8: error deallocating 3d buffers')
      return
   endif
 endif    ! do_recv

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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   real (real_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (real_kind), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 3d send,recv buffers

   real (real_kind), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

   integer (int_kind) :: len ! length of message 

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
   if (allocated(bufTripoleR4)) nxGlobal = size(bufTripoleR4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR4: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  allocate 3D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)

   allocate(bufSend(bufSizeSend*nz, halo%numMsgSend),  &
            bufRecv(bufSizeRecv*nz, halo%numMsgRecv),  &
            bufTripole(nxGlobal, halo%tripoleRows, nz), &
            stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR4: error allocating buffers')
      return
   endif

   bufTripole = fill

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)*nz
      call MPI_IRECV(bufRecv(1:len,nmsg), len, mpiR4,   &
                     halo%recvTask(nmsg),               &
                     mpitagHalo + halo%recvTask(nmsg),  &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      i=0
      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         do k=1,nz
            i = i + 1
            bufSend(i,nmsg) = array(iSrc,jSrc,k,srcBlock)
         end do
      end do
      do n=i+1,bufSizeSend*nz
         bufSend(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)*nz
      call MPI_ISEND(bufSend(1:len,nmsg), len, mpiR4, &
                     halo%sendTask(nmsg),             &
                     mpitagHalo + my_task,            &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:,:) = fill
      array(1:nx_block,ny_block-j+1,:,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:,:) = fill
      array(nx_block-i+1,1:ny_block,:,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      i = 0
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            do k=1,nz
               i = i + 1
               array(iDst,jDst,k,dstBlock) = bufRecv(i,nmsg)
            end do
         else if (dstBlock < 0) then !tripole
            do k=1,nz
               i = i + 1
               bufTripole(iDst,jDst,k) = bufRecv(i,nmsg)
            end do
         endif
      end do
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

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR4: error deallocating req,status arrays')
      return
   endif

   deallocate(bufSend, bufRecv, bufTripole, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DR4: error deallocating 3d buffers')
      return
   endif

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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   integer (int_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (int_kind), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 3d send,recv buffers

   integer (int_kind), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

   integer (int_kind) :: len ! length of message

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
   if (allocated(bufTripoleI4)) nxGlobal = size(bufTripoleI4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DI4: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  allocate 3D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)

   allocate(bufSend(bufSizeSend*nz, halo%numMsgSend),  &
            bufRecv(bufSizeRecv*nz, halo%numMsgRecv),  &
            bufTripole(nxGlobal, halo%tripoleRows, nz), &
            stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DI4: error allocating buffers')
      return
   endif

   bufTripole = fill

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)*nz
      call MPI_IRECV(bufRecv(1:len,nmsg), len, MPI_INTEGER, &
                     halo%recvTask(nmsg),                   &
                     mpitagHalo + halo%recvTask(nmsg),      &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      i=0
      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         do k=1,nz
            i = i + 1
            bufSend(i,nmsg) = array(iSrc,jSrc,k,srcBlock)
         end do
      end do
      do n=i+1,bufSizeSend*nz
         bufSend(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)*nz
      call MPI_ISEND(bufSend(1:len,nmsg), len, MPI_INTEGER, &
                     halo%sendTask(nmsg),                   &
                     mpitagHalo + my_task,                  &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:,:) = fill
      array(1:nx_block,ny_block-j+1,:,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:,:) = fill
      array(nx_block-i+1,1:ny_block,:,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      i = 0
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            do k=1,nz
               i = i + 1
               array(iDst,jDst,k,dstBlock) = bufRecv(i,nmsg)
            end do
         else if (dstBlock < 0) then !tripole
            do k=1,nz
               i = i + 1
               bufTripole(iDst,jDst,k) = bufRecv(i,nmsg)
            end do
         endif
      end do
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

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DI4: error deallocating req,status arrays')
      return
   endif

   deallocate(bufSend, bufRecv, bufTripole, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate3DI4: error deallocating 3d buffers')
      return
   endif

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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   real (dbl_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (dbl_kind), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 4d send,recv buffers

   real (dbl_kind), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

   integer (int_kind) :: len ! length of message

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
   if (allocated(bufTripoleR8)) nxGlobal = size(bufTripoleR8,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR8: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  allocate 4D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)
   nt = size(array, dim=4)

   allocate(bufSend(bufSizeSend*nz*nt, halo%numMsgSend),   &
            bufRecv(bufSizeRecv*nz*nt, halo%numMsgRecv),   &
            bufTripole(nxGlobal, halo%tripoleRows, nz, nt), &
            stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR8: error allocating buffers')
      return
   endif

   bufTripole = fill

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)*nz*nt
      call MPI_IRECV(bufRecv(1:len,nmsg), len, mpiR8,  &
                     halo%recvTask(nmsg),              &
                     mpitagHalo + halo%recvTask(nmsg), &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      i=0
      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         do l=1,nt
         do k=1,nz
            i = i + 1
            bufSend(i,nmsg) = array(iSrc,jSrc,k,l,srcBlock)
         end do
         end do
      end do

      do n=i+1,bufSizeSend*nz*nt
         bufSend(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)*nz*nt
      call MPI_ISEND(bufSend(1:len,nmsg), len, mpiR8, &
                     halo%sendTask(nmsg),             &
                     mpitagHalo + my_task,            &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:,:,:) = fill
      array(1:nx_block,ny_block-j+1,:,:,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:,:,:) = fill
      array(nx_block-i+1,1:ny_block,:,:,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      i = 0
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            do l=1,nt
            do k=1,nz
               i = i + 1
               array(iDst,jDst,k,l,dstBlock) = bufRecv(i,nmsg)
            end do
            end do
         else if (dstBlock < 0) then !tripole
            do l=1,nt
            do k=1,nz
               i = i + 1
               bufTripole(iDst,jDst,k,l) = bufRecv(i,nmsg)
            end do
            end do
         endif
      end do
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

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR8: error deallocating req,status arrays')
      return
   endif

   deallocate(bufSend, bufRecv, bufTripole, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR8: error deallocating 4d buffers')
      return
   endif

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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   real (real_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (real_kind), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 4d send,recv buffers

   real (real_kind), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

   integer (int_kind) :: len ! length of message

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
   if (allocated(bufTripoleR4)) nxGlobal = size(bufTripoleR4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR4: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  allocate 4D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)
   nt = size(array, dim=4)

   allocate(bufSend(bufSizeSend*nz*nt, halo%numMsgSend),   &
            bufRecv(bufSizeRecv*nz*nt, halo%numMsgRecv),   &
            bufTripole(nxGlobal, halo%tripoleRows, nz, nt), &
            stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR4: error allocating buffers')
      return
   endif

   bufTripole = fill

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)*nz*nt
      call MPI_IRECV(bufRecv(1:len,nmsg), len, mpiR4,  &
                     halo%recvTask(nmsg),              &
                     mpitagHalo + halo%recvTask(nmsg), &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      i=0
      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         do l=1,nt
         do k=1,nz
            i = i + 1
            bufSend(i,nmsg) = array(iSrc,jSrc,k,l,srcBlock)
         end do
         end do
      end do

      do n=i+1,bufSizeSend*nz*nt
         bufSend(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)*nz*nt
      call MPI_ISEND(bufSend(1:len,nmsg), len, mpiR4, &
                     halo%sendTask(nmsg),             &
                     mpitagHalo + my_task,            &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:,:,:) = fill
      array(1:nx_block,ny_block-j+1,:,:,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:,:,:) = fill
      array(nx_block-i+1,1:ny_block,:,:,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      i = 0
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            do l=1,nt
            do k=1,nz
               i = i + 1
               array(iDst,jDst,k,l,dstBlock) = bufRecv(i,nmsg)
            end do
            end do
         else if (dstBlock < 0) then !tripole
            do l=1,nt
            do k=1,nz
               i = i + 1
               bufTripole(iDst,jDst,k,l) = bufRecv(i,nmsg)
            end do
            end do
         endif
      end do
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

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR4: error deallocating req,status arrays')
      return
   endif

   deallocate(bufSend, bufRecv, bufTripole, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DR4: error deallocating 4d buffers')
      return
   endif

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
!  ice\_HaloUpdate.  This routine is the specific interface
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

   integer (int_kind), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      sndStatus,       &! MPI status flags
      rcvStatus         ! MPI status flags

   integer (int_kind) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (int_kind), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 4d send,recv buffers

   integer (int_kind), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

   integer (int_kind) :: len  ! length of messages

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
   if (allocated(bufTripoleI4)) nxGlobal = size(bufTripoleI4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request and status arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            sndStatus(MPI_STATUS_SIZE,halo%numMsgSend), &
            rcvStatus(MPI_STATUS_SIZE,halo%numMsgRecv), stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DI4: error allocating req,status arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  allocate 4D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)
   nt = size(array, dim=4)

   allocate(bufSend(bufSizeSend*nz*nt, halo%numMsgSend),   &
            bufRecv(bufSizeRecv*nz*nt, halo%numMsgRecv),   &
            bufTripole(nxGlobal, halo%tripoleRows, nz, nt), &
            stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DI4: error allocating buffers')
      return
   endif

   bufTripole = fill

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      len = halo%SizeRecv(nmsg)*nz*nt
      call MPI_IRECV(bufRecv(1:len,nmsg), len, MPI_INTEGER, &
                     halo%recvTask(nmsg),                   &
                     mpitagHalo + halo%recvTask(nmsg),      &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgSend

      i=0
      do n=1,halo%sizeSend(nmsg)
         iSrc     = halo%sendAddr(1,n,nmsg)
         jSrc     = halo%sendAddr(2,n,nmsg)
         srcBlock = halo%sendAddr(3,n,nmsg)

         do l=1,nt
         do k=1,nz
            i = i + 1
            bufSend(i,nmsg) = array(iSrc,jSrc,k,l,srcBlock)
         end do
         end do
      end do

      do n=i+1,bufSizeSend*nz*nt
         bufSend(n,nmsg) = fill  ! fill remainder of buffer
      end do

      len = halo%SizeSend(nmsg)*nz*nt
      call MPI_ISEND(bufSend(1:len,nmsg), len, MPI_INTEGER, &
                     halo%sendTask(nmsg),                   &
                     mpitagHalo + my_task,                  &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  while messages are being communicated,
!  fill out halo region, needed for masked halos to ensure
!  halo values are fill for halo gridcells that are not updated
!
!-----------------------------------------------------------------------

   do j = 1,nghost
      array(1:nx_block,           j,:,:,:) = fill
      array(1:nx_block,ny_block-j+1,:,:,:) = fill
   enddo
   do i = 1,nghost
      array(i,           1:ny_block,:,:,:) = fill
      array(nx_block-i+1,1:ny_block,:,:,:) = fill
   enddo

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
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, rcvStatus, ierr)

   do nmsg=1,halo%numMsgRecv
      i = 0
      do n=1,halo%sizeRecv(nmsg)
         iDst     = halo%recvAddr(1,n,nmsg)
         jDst     = halo%recvAddr(2,n,nmsg)
         dstBlock = halo%recvAddr(3,n,nmsg)

         if (dstBlock > 0) then
            do l=1,nt
            do k=1,nz
               i = i + 1
               array(iDst,jDst,k,l,dstBlock) = bufRecv(i,nmsg)
            end do
            end do
         else if (dstBlock < 0) then !tripole
            do l=1,nt
            do k=1,nz
               i = i + 1
               bufTripole(iDst,jDst,k,l) = bufRecv(i,nmsg)
            end do
            end do
         endif
      end do
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

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(halo%numMsgSend, sndRequest, sndStatus, ierr)

   deallocate(sndRequest, rcvRequest, sndStatus, rcvStatus, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DI4: error deallocating req,status arrays')
      return
   endif

   deallocate(bufSend, bufRecv, bufTripole, stat=ierr)

   if (ierr > 0) then
      call abort_ice( &
         'ice_HaloUpdate4DI4: error deallocating 4d buffers')
      return
   endif

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

   subroutine ice_HaloMsgCreate(halo, srcBlock, srcProc, srcLocalID, &
                                      dstBlock, dstProc, dstLocalID, &
                                      direction)

! !DESCRIPTION:
!  This is a utility routine to determine the required address and
!  message information for a particular pair of blocks.

! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      srcBlock,   dstBlock,   & ! source,destination block id
      srcProc,    dstProc,    & ! source,destination processor location
      srcLocalID, dstLocalID    ! source,destination local index

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
      msgIndx,               &! message counter and index into msg array
      blockIndx,             &! block counter and index into msg array
      bufSize,               &! size of message buffer
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

      msgIndx = halo%numLocalCopies

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call abort_ice( &
            'ice_HaloMsgCreate: msg count > array size')
         return
      endif

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

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

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

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

      case ('north')

         !*** copy northern physical domain of src
         !*** into southern halo of dst

         if (srcBlock > 0 .and. dstBlock > 0) then  ! normal north boundary

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

         else if (srcBlock < 0 .and. dstBlock > 0) then

            !*** tripole grid - set up for copying out of 
            !*** tripole buffer into ghost cell domains
            !*** include e-w ghost cells

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

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

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

      case ('northeast')

         !*** normal northeast boundary - just copy NE corner
         !*** of physical domain into SW halo of NE nbr block

         if (dstBlock > 0) then

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

         else

            !*** tripole grid - this local copy should already
            !*** have taken place for the north boundary

         endif

      case ('northwest')

         !*** normal northeast boundary - just copy NW corner
         !*** of physical domain into SE halo of NW nbr block

         if (dstBlock > 0) then

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

         else

            !*** tripole grid - this local copy should already
            !*** have taken place for the north boundary

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

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

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

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

      case default

         call abort_ice( &
            'ice_HaloMsgCreate: unknown direction local copy')
         return

      end select

      halo%numLocalCopies = msgIndx

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call abort_ice( &
            'ice_HaloMsgCreate: msg count > array size')
         return
      endif

!-----------------------------------------------------------------------
!
!  if dest block is local and source block does not exist, create a 
!  local copy to fill halo with a fill value
!
!-----------------------------------------------------------------------

   else if (srcProc == 0 .and. dstProc == my_task+1) then   

      msgIndx = halo%numLocalCopies

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call abort_ice( &
            'ice_HaloMsgCreate: msg count > array size')
         return
      endif

      !*** compute addresses based on direction

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

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

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

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

      case ('north')

         !*** copy northern physical domain of src
         !*** into southern halo of dst

         if (dstBlock > 0) then  ! normal north boundary

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

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

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

      case ('northeast')

         !*** normal northeast boundary - just copy NE corner
         !*** of physical domain into SW halo of NE nbr block

         if (dstBlock > 0) then

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

         endif

      case ('northwest')

         !*** normal northeast boundary - just copy NW corner
         !*** of physical domain into SE halo of NW nbr block

         if (dstBlock > 0) then

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

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

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

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

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

      case default

         call abort_ice( &
            'ice_HaloMsgCreate: unknown direction local copy')
         return

      end select

      halo%numLocalCopies = msgIndx

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call abort_ice( &
            'ice_HaloMsgCreate: msg count > array size')
         return
      endif

!-----------------------------------------------------------------------
!
!  if source block local and dest block remote, send a message
!
!-----------------------------------------------------------------------

   else if (srcProc == my_task+1 .and. &
            dstProc /= my_task+1 .and. dstProc > 0) then

      !*** first check to see if a message to this processor has
      !*** already been defined
      !*** if not, update counters and indices

      msgIndx = 0

      srchSend: do n=1,halo%numMsgSend
         if (halo%sendTask(n) == dstProc - 1) then
            msgIndx = n
            bufSize = halo%sizeSend(n)
            exit srchSend
         endif
      end do srchSend 

      if (msgIndx == 0) then
         msgIndx = halo%numMsgSend + 1
         halo%numMsgSend = msgIndx
         halo%sendTask(msgIndx) = dstProc - 1
         bufSize = 0
      endif

      !*** now compute message info based on msg direction

      select case (direction)
      case ('east')

         !*** send easternmost physical domain of src
         !*** into westernmost halo of dst

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            bufSize = bufSize + 1

            halo%sendAddr(1,bufSize,msgIndx) = ieSrc - nghost + i
            halo%sendAddr(2,bufSize,msgIndx) = jbSrc + j - 1
            halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

         end do
         end do

         halo%sizeSend(msgIndx) = bufSize

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            bufSize = bufSize + 1

            halo%sendAddr(1,bufSize,msgIndx) = ibSrc + i - 1
            halo%sendAddr(2,bufSize,msgIndx) = jbSrc + j - 1
            halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

         end do
         end do

         halo%sizeSend(msgIndx) = bufSize

      case ('north')

         if (dstBlock > 0) then

            !*** copy northern physical domain of src
            !*** into southern halo of dst

            do j=1,nghost
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx) = ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx) = jeSrc-nghost+j
               halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         else 

            !*** tripole block - send top halo%tripoleRows rows of phys domain

            halo%tripSend(msgIndx) = 1
            do j=1,halo%tripoleRows
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx)=ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx)=jeSrc-halo%tripoleRows+j
               halo%sendAddr(3,bufSize,msgIndx)=srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         do j=1,nghost
         do i=1,ieSrc-ibSrc+1

            bufSize = bufSize + 1

            halo%sendAddr(1,bufSize,msgIndx) = ibSrc + i - 1
            halo%sendAddr(2,bufSize,msgIndx) = jbSrc + j - 1
            halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

         end do
         end do

         halo%sizeSend(msgIndx) = bufSize

      case ('northeast')


         if (dstBlock > 0) then

            !*** normal northeast corner
            !*** copy northeast corner of src physical domain
            !*** into southwestern halo of dst

            do j=1,nghost
            do i=1,nghost

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx) = ieSrc-nghost+i
               halo%sendAddr(2,bufSize,msgIndx) = jeSrc-nghost+j
               halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         else 

            !*** tripole block - send top halo%tripoleRows rows of phys domain

            halo%tripSend(msgIndx) = 1
            do j=1,halo%tripoleRows
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx)=ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx)=jeSrc-halo%tripoleRows+j
               halo%sendAddr(3,bufSize,msgIndx)=srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         endif

      case ('northwest')

         if (dstBlock > 0) then

            !*** normal northwest corner
            !*** copy northwest corner of src physical domain
            !*** into southeastern halo of dst

            do j=1,nghost
            do i=1,nghost

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx) = ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx) = jeSrc-nghost+j
               halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         else 

            !*** tripole block - send top halo%tripoleRows rows of phys domain

            halo%tripSend(msgIndx) = 1
            do j=1,halo%tripoleRows
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx)=ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx)=jeSrc-halo%tripoleRows+j
               halo%sendAddr(3,bufSize,msgIndx)=srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

         do j=1,nghost
         do i=1,nghost

            bufSize = bufSize + 1

            halo%sendAddr(1,bufSize,msgIndx) = ieSrc - nghost + i
            halo%sendAddr(2,bufSize,msgIndx) = jbSrc + j - 1
            halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

         end do
         end do

         halo%sizeSend(msgIndx) = bufSize

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

         do j=1,nghost
         do i=1,nghost

            bufSize = bufSize + 1

            halo%sendAddr(1,bufSize,msgIndx) = ibSrc + i - 1
            halo%sendAddr(2,bufSize,msgIndx) = jbSrc + j - 1
            halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

         end do
         end do

         halo%sizeSend(msgIndx) = bufSize

      case default

         !*** already checked in previous case construct

      end select

!-----------------------------------------------------------------------
!
!  if source block remote and dest block local, recv a message
!
!-----------------------------------------------------------------------

   else if (dstProc == my_task+1 .and. &
            srcProc /= my_task+1 .and. srcProc > 0) then

      !*** first check to see if a message from this processor has
      !*** already been defined
      !*** if not, update counters and indices

      msgIndx = 0

      srchRecv: do n=1,halo%numMsgRecv
         if (halo%recvTask(n) == srcProc - 1) then
            msgIndx = n
            bufSize = halo%sizeRecv(n)
            exit srchRecv
         endif
      end do srchRecv 

      if (msgIndx == 0) then
         msgIndx = halo%numMsgRecv + 1
         halo%numMsgRecv = msgIndx
         halo%recvTask(msgIndx) = srcProc - 1
         bufSize = 0
      endif

      !*** now compute message info based on msg direction

      select case (direction)
      case ('east')

         !*** send easternmost physical domain of src
         !*** into westernmost halo of dst

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            bufSize = bufSize + 1

            halo%recvAddr(1,bufSize,msgIndx) = i
            halo%recvAddr(2,bufSize,msgIndx) = jbDst + j - 1
            halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

         end do
         end do

         halo%sizeRecv(msgIndx) = bufSize

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

         do j=1,jeSrc-jbSrc+1
         do i=1,nghost

            bufSize = bufSize + 1

            halo%recvAddr(1,bufSize,msgIndx) = ieDst + i
            halo%recvAddr(2,bufSize,msgIndx) = jbDst + j - 1
            halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

         end do
         end do

         halo%sizeRecv(msgIndx) = bufSize

      case ('north')

         if (dstBlock > 0) then

            !*** copy northern physical domain of src
            !*** into southern halo of dst

            do j=1,nghost
            do i=1,ieDst-ibDst+1

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = ibDst + i - 1
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         else

            !*** tripole block - receive into tripole buffer

            halo%tripRecv(msgIndx) = 1
            do j=1,halo%tripoleRows
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = iGlobal(ibSrc + i - 1)
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = -dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         do j=1,nghost
         do i=1,ieSrc-ibSrc+1

            bufSize = bufSize + 1

            halo%recvAddr(1,bufSize,msgIndx) = ibDst + i - 1
            halo%recvAddr(2,bufSize,msgIndx) = jeDst + j
            halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

         end do
         end do

         halo%sizeRecv(msgIndx) = bufSize

      case ('northeast')

         if (dstBlock > 0) then

            !*** normal northeast neighbor
            !*** copy northeast physical domain into
            !*** into southwest halo of dst

            do j=1,nghost
            do i=1,nghost

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = i
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         else

            !*** tripole block - receive into tripole buffer

            halo%tripRecv(msgIndx) = 1
            do j=1,halo%tripoleRows
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = iGlobal(ibSrc + i - 1)
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = -dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         endif

      case ('northwest')

         if (dstBlock > 0) then

            !*** normal northwest neighbor
            !*** copy northwest physical domain into
            !*** into southeast halo of dst

            do j=1,nghost
            do i=1,nghost

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = ieDst + i
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         else

            !*** tripole block - receive into tripole buffer

            halo%tripRecv(msgIndx) = 1
            do j=1,halo%tripoleRows
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = iGlobal(ibSrc + i - 1)
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = -dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

         do j=1,nghost
         do i=1,nghost

            bufSize = bufSize + 1

            halo%recvAddr(1,bufSize,msgIndx) = i
            halo%recvAddr(2,bufSize,msgIndx) = jeDst + j
            halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

         end do
         end do

         halo%sizeRecv(msgIndx) = bufSize

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

         do j=1,nghost
         do i=1,nghost

            bufSize = bufSize + 1

            halo%recvAddr(1,bufSize,msgIndx) = ieDst + i
            halo%recvAddr(2,bufSize,msgIndx) = jeDst + j
            halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

         end do
         end do

         halo%sizeRecv(msgIndx) = bufSize

      case default

         !*** already checked in previous case construct

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
