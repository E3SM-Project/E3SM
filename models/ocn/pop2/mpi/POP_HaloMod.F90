!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: POP_HaloMod

!pw++
#define _SAVE_3D_BUFFERS 1
#define _SAVE_4D_BUFFERS 1
!pw
!pw #define _ALTREQ 1
!pw #define _WAITANY_2D 1
!pw #define _WAITANY_3D 1
!pw #define _WAITANY_4D 1
!pw--

 module POP_HaloMod

! !DESCRIPTION:
!  This module contains data types and routines for updating halo
!  regions (ghost cells) using MPI calls
!
! !REVISION HISTORY:
!  SVN:$Id$
!  2007-07-19: Phil Jones, Yoshi Yoshida, John Dennis
!              new naming conventions, optimizations during
!              initialization, true multi-dimensional updates 
!              (rather than serial call to two-dimensional updates), 
!              fixes for non-existent blocks
!  2008-01-30: Phil Jones, Elizabeth Hunke
!              fixed some bugs Elizabeth found with one assumption
!                of halo width 2, a typo on sw/se nbr, and
!                tripole buffers uninitialized

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod
   use POP_CommMod
   use POP_BlocksMod
   use POP_ReductionsMod
   use POP_DistributionMod
   use POP_FieldMod
   use POP_GridHorzMod
!pw++
   use perf_mod
!pw--

   implicit none
   private
   save

   include 'mpif.h'

! !PUBLIC TYPES:

   type, public :: POP_halo
      integer (POP_i4) ::  &
         communicator,     &! communicator to use for update messages
         numMsgSend,       &! number of messages to send halo update
         numMsgRecv,       &! number of messages to recv halo update
         numLocalCopies     ! num local copies for halo update

      integer (POP_i4), dimension(:), pointer :: &
         recvTask,         &! task from which to recv each msg
         sendTask,         &! task to   which to send each msg
         sizeSend,         &! size of each sent message
         sizeRecv           ! size of each recvd message

      integer (POP_i4), dimension(:,:), pointer :: &
         srcLocalAddr,     &! src addresses for each local copy
         dstLocalAddr       ! dst addresses for each local copy

      integer (POP_i4), dimension(:,:,:), pointer :: &
         sendAddr,         &! src addresses for each sent message
         recvAddr           ! dst addresses for each recvd message

   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_HaloCreate,  &
             POP_HaloDestroy, &
             POP_HaloUpdate,  &
             POP_HaloPrintStats

   interface POP_HaloUpdate  ! generic interface
      module procedure POP_HaloUpdate2DR8, &
                       POP_HaloUpdate2DR4, &
                       POP_HaloUpdate2DI4, &
                       POP_HaloUpdate3DR8, &
                       POP_HaloUpdate3DR4, &
                       POP_HaloUpdate3DI4, &
                       POP_HaloUpdate4DR8, &
                       POP_HaloUpdate4DR4, &
                       POP_HaloUpdate4DI4
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

   integer (POP_i4) :: &
      bufSizeSend,    &! max buffer size for send messages
      bufSizeRecv      ! max buffer size for recv messages

   integer (POP_i4), dimension(:,:), allocatable :: &
      bufSendI4,     &! buffer for use to send in 2D i4 halo updates
      bufRecvI4       ! buffer for use to recv in 2D i4 halo updates

   real (POP_r4), dimension(:,:), allocatable :: &
      bufSendR4,     &! buffer for use to send in 2D r4 halo updates
      bufRecvR4       ! buffer for use to recv in 2D r4 halo updates

   real (POP_r8), dimension(:,:), allocatable :: &
      bufSendR8,     &! buffer for use to send in 2D r8 halo updates
      bufRecvR8       ! buffer for use to recv in 2D r8 halo updates

!-----------------------------------------------------------------------
!
!  global buffers for tripole boundary
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(:,:), allocatable :: &
      bufTripoleI4

   real (POP_r4), dimension(:,:), allocatable :: &
      bufTripoleR4

   real (POP_r8), dimension(:,:), allocatable :: &
      bufTripoleR8

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloCreate
! !INTERFACE:

 function POP_HaloCreate(distrb, nsBoundaryType, ewBoundaryType, &
                         nxGlobal, errorCode)  result(halo)

! !DESCRIPTION:
!  This routine creates a halo type with info necessary for
!  performing a halo (ghost cell) update. This info is computed
!  based on the input block distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      distrb             ! distribution of blocks across procs

   character (*), intent(in) :: &
      nsBoundaryType,   &! type of boundary to use in logical ns dir
      ewBoundaryType     ! type of boundary to use in logical ew dir

   integer (POP_i4), intent(in) :: &
      nxGlobal           ! global grid extent for tripole grids

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

   type (POP_halo) :: &
      halo               ! a new halo type with info for halo updates

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::             &
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
      cornerMsgSize, msgSize        ! nominal size for corner msg

   integer (POP_i4), dimension(:), allocatable :: &
      sendCount, recvCount          ! count number of words to each proc

   logical (POP_logical) :: &
      resize,               &! flag for resizing buffers
      tripoleFlag,          &! flag for allocating tripole buffers
      tripoleBlock           ! flag for identifying north tripole blocks

!-----------------------------------------------------------------------
!
!  Initialize some useful variables and return if this task not
!  in the current distribution.
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(distrb, errorCode,         &
                            numProcs = numProcs,       &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
                        'POP_HaloCreate: error getting distrb info')
      return
   endif

   if (POP_myTask >= numProcs) return

   halo%communicator = communicator

   blockSizeX = POP_nxBlock - 2*POP_haloWidth
   blockSizeY = POP_nyBlock - 2*POP_haloWidth
   eastMsgSize  = POP_haloWidth*blockSizeY
   westMsgSize  = POP_haloWidth*blockSizeY
   southMsgSize = POP_haloWidth*blockSizeX
   northMsgSize = POP_haloWidth*blockSizeX
   cornerMsgSize = POP_haloWidth*POP_haloWidth
   tripoleMsgSize = (POP_haloWidth+1)*blockSizeX
   tripoleMsgSizeOut = (POP_haloWidth+1)*POP_nxBlock

   if (nsBoundaryType == 'tripole') then
      tripoleFlag = .true.

      !*** allocate tripole message buffers if not already done

      if (.not. allocated(bufTripoleR8)) then
         allocate (bufTripoleI4(nxGlobal, POP_haloWidth+1), &
                   bufTripoleR4(nxGlobal, POP_haloWidth+1), &
                   bufTripoleR8(nxGlobal, POP_haloWidth+1), &
                   stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error allocating tripole buffers')
            return
         endif
      endif

   else
      tripoleFlag = .false.
   endif

!-----------------------------------------------------------------------
!
!  Count the number of messages to send/recv from each processor
!  and number of words in each message.  These quantities are
!  necessary for allocating future arrays.
!
!-----------------------------------------------------------------------

   allocate (sendCount(numProcs), recvCount(numProcs), stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
                        'POP_HaloCreate: error allocating count arrays')
      return
   endif

   sendCount  = 0
   recvCount  = 0

   msgCountLoop: do iblock=1,POP_numBlocks

      call POP_DistributionGetBlockLoc(distrb, iblock, srcProc, &
                                       srcLocalID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding source block location')
         return
      endif

      !*** find north neighbor block and add to message count
      !***  also set tripole block flag for later special cases

      northBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorth,        &
                                      ewBoundaryType, nsBoundaryType, &
                                      errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding north neighbor block')
         return
      endif

      if (northBlock > 0) then
         tripoleBlock = .false.
         msgSize = northMsgSize
         call POP_DistributionGetBlockLoc(distrb, northBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding north block location')
            return
         endif
      else if (northBlock < 0) then ! tripole north row, count block
         tripoleBlock = .true.
         msgSize = tripoleMsgSize
         call POP_DistributionGetBlockLoc(distrb, abs(northBlock), &
                                 dstProc, dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding north block location')
            return
         endif
      else
         tripoleBlock = .false.
         msgSize = northMsgSize
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,           &
                                     srcProc, dstProc, msgSize,      &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing north message')
         return
      endif

      !*** if a tripole boundary block, also create a local
      !*** message into and out of tripole buffer 

      if (tripoleBlock) then
         !*** copy out of tripole buffer - includes halo
         call POP_HaloIncrementMsgCount(sendCount, recvCount,        &
                                        srcProc, srcProc,            &
                                        tripoleMsgSizeOut, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing tripole copy out count')
            return
         endif

         !*** copy in only required if dstProc not same as srcProc
         if (dstProc /= srcProc) then
            call POP_HaloIncrementMsgCount(sendCount, recvCount,  &
                                           srcProc, srcProc,      & 
                                           msgSize, errorCode)
            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error incrementing tripole copy in')
               return
            endif

         endif
      endif

      !*** find south neighbor block and add to message count

      southBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSouth,        &
                                      ewBoundaryType, nsBoundaryType, &
                                      errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding south neighbor block')
         return
      endif

      if (southBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, southBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding south block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,           &
                                     srcProc, dstProc, southMsgSize, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing south message count')
         return
      endif

      !*** find east neighbor block and add to message count

      eastBlock = POP_BlocksGetNbrID(iblock, POP_BlocksEast,         &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding east neighbor block')
         return
      endif

      if (eastBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, eastBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding east block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, eastMsgSize, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing east message count')
         return
      endif

      !*** if a tripole boundary block, non-local east neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

      if (tripoleBlock .and. dstProc /= srcProc) then
         call POP_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, tripoleMsgSize, &
                                     errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing tripole east msg count')
            return
         endif
      endif

      !*** find west neighbor block and add to message count

      westBlock = POP_BlocksGetNbrID(iblock, POP_BlocksWest,         &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding west neighbor block')
         return
      endif

      if (westBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, westBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding west block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, westMsgSize, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing west message count')
         return
      endif

      !*** if a tripole boundary block, non-local west neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

      if (tripoleBlock .and. dstProc /= srcProc) then
         call POP_HaloIncrementMsgCount(sendCount, recvCount,          &
                                     srcProc, dstProc, tripoleMsgSize, &
                                     errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing tripole west msg count')
            return
         endif
      endif

      !*** find northeast neighbor block and add to message count

      neBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorthEast,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding northeast neighbor block')
         return
      endif

      if (neBlock > 0) then
         msgSize = cornerMsgSize  ! normal corner message 

         call POP_DistributionGetBlockLoc(distrb, neBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding northeast block location')
            return
         endif

      else if (neBlock < 0) then ! tripole north row
         msgSize = tripoleMsgSize  ! tripole needs whole top row of block

         call POP_DistributionGetBlockLoc(distrb, abs(neBlock), dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding northeast block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,      &
                                     srcProc, dstProc, msgSize, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing northeast message count')
         return
      endif

      !*** find northwest neighbor block and add to message count

      nwBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorthWest,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding northwest neighbor block')
         return
      endif

      if (nwBlock > 0) then
         msgSize = cornerMsgSize ! normal NE corner update

         call POP_DistributionGetBlockLoc(distrb, nwBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding northwest block location')
            return
         endif

      else if (nwBlock < 0) then ! tripole north row, count block
         msgSize = tripoleMsgSize ! tripole NE corner update - entire row needed

         call POP_DistributionGetBlockLoc(distrb, abs(nwBlock), dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding northwest block location')
            return
         endif

      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,      &
                                     srcProc, dstProc, msgSize, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing northwest message count')
         return
      endif

      !*** find southeast neighbor block and add to message count

      seBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSoutheast,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding southeast neighbor block')
         return
      endif

      if (seBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, seBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding southeast block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,            &
                                     srcProc, dstProc, cornerMsgSize, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing southeast message count')
         return
      endif

      !*** find southwest neighbor block and add to message count

      swBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSouthWest,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding southwest neighbor block')
         return
      endif

      if (swBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, swBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding southwest block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloIncrementMsgCount(sendCount, recvCount,            &
                                     srcProc, dstProc, cornerMsgSize, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing southwest message count')
         return
      endif

      !*** for tripole grids with padded domain, padding will
      !*** prevent tripole buffer from getting all the info
      !*** it needs - must extend footprint at top boundary

      if (tripoleBlock                  .and. & !tripole
          mod(nxGlobal,blockSizeX) /= 0) then   !padding

         !*** find east2 neighbor block and add to message count

         eastBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksEast2,     &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding east2 neighbor block')
            return
         endif

         if (eastBlock > 0) then
            call POP_DistributionGetBlockLoc(distrb, eastBlock, dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding east2 block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloIncrementMsgCount(sendCount, recvCount,       &
                                     srcProc, dstProc, tripoleMsgSize, &
                                     errorCode)
            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error incrementing tripole east2 msg count')
               return
            endif
         endif

         !*** find EastNorthEast neighbor block and add to message count

         neBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksEastNorthEast, &
                                     ewBoundaryType, nsBoundaryType,     &
                                     errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding EastNorthEast neighbor block')
            return
         endif

         if (neBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call POP_DistributionGetBlockLoc(distrb, abs(neBlock), dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding EastNorthEast block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloIncrementMsgCount(sendCount, recvCount,   &
                                        srcProc, dstProc, msgSize, &
                                        errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error incrementing northeast2 message count')
               return
            endif
         endif

         !*** find west2 neighbor block and add to message count

         westBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksWest2,     &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding west2 neighbor block')
            return
         endif

         if (westBlock > 0) then
            call POP_DistributionGetBlockLoc(distrb, westBlock, dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding west2 block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloIncrementMsgCount(sendCount, recvCount,       &
                                     srcProc, dstProc, tripoleMsgSize, &
                                     errorCode)
            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error incrementing tripole west2 msg count')
               return
            endif
         endif

         !*** find WestNorthWest neighbor block and add to message count

         nwBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksWestNorthWest, &
                                     ewBoundaryType, nsBoundaryType,   &
                                     errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding WestNorthWest neighbor block')
            return
         endif

         if (nwBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call POP_DistributionGetBlockLoc(distrb, abs(nwBlock), dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding northwest2 block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloIncrementMsgCount(sendCount, recvCount,   &
                                        srcProc, dstProc, msgSize, &
                                        errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error incrementing northwest2 message count')
               return
            endif
         endif

      endif

   end do msgCountLoop

!-----------------------------------------------------------------------
!
!  if messages are received from the same processor, the message is 
!  actually a local copy - count them and reset to zero
!
!-----------------------------------------------------------------------

   halo%numLocalCopies = recvCount(POP_myTask+1)

   sendCount(POP_myTask+1) = 0
   recvCount(POP_myTask+1) = 0

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
   maxSizeSend = POP_GlobalMaxval(maxTmp, distrb, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
           'POP_HaloCreate: error computing size of send buffer')
      return
   endif

   maxTmp = maxval(recvCount)
   maxSizeRecv = POP_GlobalMaxval(maxTmp, distrb, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
           'POP_HaloCreate: error computing size of recv buffer')
      return
   endif

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
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error allocating 2d buffers')
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
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error deallocating 2d buffers')
            return
         endif

         allocate(bufSendI4(bufSizeSend, numMsgSend), &
                  bufRecvI4(bufSizeRecv, numMsgRecv), &
                  bufSendR4(bufSizeSend, numMsgSend), &
                  bufRecvR4(bufSizeRecv, numMsgRecv), &
                  bufSendR8(bufSizeSend, numMsgSend), &
                  bufRecvR8(bufSizeRecv, numMsgRecv), stat=istat)

         if (istat > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error reallocating 2d buffers')
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
            halo%sendAddr(3,bufSizeSend,numMsgSend), &
            halo%recvAddr(3,bufSizeRecv,numMsgRecv), &
            halo%srcLocalAddr(3,halo%numLocalCopies), &
            halo%dstLocalAddr(3,halo%numLocalCopies), &
            stat = istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloCreate: error allocating halo message info arrays')
      return
   endif

   halo%sendTask = 0
   halo%recvTask = 0
   halo%sizeSend = 0
   halo%sizeRecv = 0
   halo%sendAddr = 0
   halo%recvAddr = 0
   halo%srcLocalAddr = 0
   halo%dstLocalAddr = 0

   deallocate(sendCount, recvCount, stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloCreate: error deallocating count arrays')
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

   msgConfigLoop: do iblock=1,POP_numBlocks

      call POP_DistributionGetBlockLoc(distrb, iblock, srcProc, &
                                       srcLocalID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding source block location')
         return
      endif

      !*** find north neighbor block and set msg info
      !***  also set tripole block flag for later special cases

      northBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorth,        &
                                      ewBoundaryType, nsBoundaryType, &
                                      errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding north neighbor block')
         return
      endif

      if (northBlock > 0) then
         tripoleBlock = .false.
         call POP_DistributionGetBlockLoc(distrb, northBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding north block location')
            return
         endif
      else if (northBlock < 0) then ! tripole north row, count block
         tripoleBlock = .true.
         call POP_DistributionGetBlockLoc(distrb, abs(northBlock), &
                                 dstProc, dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding north block location')
            return
         endif
      else
         tripoleBlock = .false.
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,     srcProc, srcLocalID, &
                                   northBlock, dstProc, dstLocalID, &
                                   'north', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating north message')
         return
      endif

      !*** if a tripole boundary block, also create a local
      !*** message into and out of tripole buffer 

      if (tripoleBlock) then
         !*** copy out of tripole buffer - includes halo
         call POP_HaloMsgCreate(halo,-iblock, srcProc, srcLocalID, &
                                      iblock, srcProc, srcLocalID, &
                                      'north', errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating tripole copy out of buf')
            return
         endif

         !*** copy in only required if dstProc not same as srcProc
         if (dstProc /= srcProc) then
            call POP_HaloMsgCreate(halo, iblock, srcProc, srcLocalID, &
                                        -iblock, srcProc, srcLocalID, &
                                         'north', errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error creating tripole copy into buf')
               return
            endif

         endif
      endif

      !*** find south neighbor block and add to message count

      southBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSouth,        &
                                      ewBoundaryType, nsBoundaryType, &
                                      errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding south neighbor block')
         return
      endif

      if (southBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, southBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding south block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,     srcProc, srcLocalID, &
                                   southBlock, dstProc, dstLocalID, &
                                   'south', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating south message')
         return
      endif

      !*** find east neighbor block and add to message count

      eastBlock = POP_BlocksGetNbrID(iblock, POP_BlocksEast,         &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding east neighbor block')
         return
      endif

      if (eastBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, eastBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding east block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                   eastBlock, dstProc, dstLocalID, &
                                   'east', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating east message')
         return
      endif

      !*** if a tripole boundary block, non-local east neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

      if (tripoleBlock .and. dstProc /= srcProc) then
         call POP_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                     -eastBlock, dstProc, dstLocalID, &
                                      'north', errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating tripole east msg')
            return
         endif
      endif

      !*** find west neighbor block and add to message count

      westBlock = POP_BlocksGetNbrID(iblock, POP_BlocksWest,         &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding west neighbor block')
         return
      endif

      if (westBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, westBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding west block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                   westBlock, dstProc, dstLocalID, &
                                   'west', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating west message')
         return
      endif

      !*** if a tripole boundary block, non-local west neighbor
      !*** needs a chunk of the north boundary, so add a message
      !*** for that

      if (tripoleBlock .and. dstProc /= srcProc) then
         call POP_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                     -westBlock, dstProc, dstLocalID, &
                                      'north', errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating tripole west msg')
            return
         endif
      endif

      !*** find northeast neighbor block and add to message count

      neBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorthEast,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding northeast neighbor block')
         return
      endif

      if (neBlock /= 0) then
         call POP_DistributionGetBlockLoc(distrb, abs(neBlock), dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding northeast block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   neBlock, dstProc, dstLocalID, &
                                   'northeast', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating northeast message')
         return
      endif

      !*** find northwest neighbor block and add to message count

      nwBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorthWest,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding northwest neighbor block')
         return
      endif

      if (nwBlock /= 0) then
         call POP_DistributionGetBlockLoc(distrb, abs(nwBlock), dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding northwest block location')
            return
         endif

      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   nwBlock, dstProc, dstLocalID, &
                                   'northwest', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating northwest message')
         return
      endif

      !*** find southeast neighbor block and add to message count

      seBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSoutheast,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding southeast neighbor block')
         return
      endif

      if (seBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, seBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding southeast block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   seBlock, dstProc, dstLocalID, &
                                   'southeast', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating southeast message')
         return
      endif

      !*** find southwest neighbor block and add to message count

      swBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSouthWest,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding southwest neighbor block')
         return
      endif

      if (swBlock > 0) then
         call POP_DistributionGetBlockLoc(distrb, swBlock, dstProc, &
                                          dstLocalID, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding southwest block location')
            return
         endif
      else
         dstProc = 0
         dstLocalID = 0
      endif

      call POP_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                   swBlock, dstProc, dstLocalID, &
                                   'southwest', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating southwest message')
         return
      endif

      !*** for tripole grids with padded domain, padding will
      !*** prevent tripole buffer from getting all the info
      !*** it needs - must extend footprint at top boundary

      if (tripoleBlock                  .and. & !tripole
          mod(nxGlobal,blockSizeX) /= 0) then   !padding

         !*** find east2 neighbor block and add to message count

         eastBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksEast2,     &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding east2 neighbor block')
            return
         endif

         if (eastBlock > 0) then
            call POP_DistributionGetBlockLoc(distrb, eastBlock, dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding east2 block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                        -eastBlock, dstProc, dstLocalID, &
                                         'north', errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error creating tripole east2 msg')
               return
            endif
         endif

         !*** find EastNorthEast neighbor block and add to message count

         neBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksEastNorthEast, &
                                     ewBoundaryType, nsBoundaryType,     &
                                     errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding EastNorthEast neighbor block')
            return
         endif

         if (neBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call POP_DistributionGetBlockLoc(distrb, abs(neBlock), dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding EastNorthEast block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                         neBlock, dstProc, dstLocalID, &
                                         'north', errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error creating EastNorthEast message')
               return
            endif

         endif

         !*** find west2 neighbor block and add to message count

         westBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksWest2,     &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding west2 neighbor block')
            return
         endif

         if (westBlock > 0) then
            call POP_DistributionGetBlockLoc(distrb, westBlock, dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding west2 block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloMsgCreate(halo, iblock,    srcProc, srcLocalID, &
                                        -westBlock, dstProc, dstLocalID, &
                                         'north', errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error creating tripole west msg')
               return
            endif
         endif

         !*** find WestNorthWest neighbor block and add to message count

         nwBlock = POP_BlocksGetNbrID(iBlock, POP_BlocksWestNorthWest, &
                                     ewBoundaryType, nsBoundaryType,   &
                                     errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error finding WestNorthWest neighbor block')
            return
         endif

         if (nwBlock < 0) then ! tripole north row
            msgSize = tripoleMsgSize  ! tripole needs whole top row of block

            call POP_DistributionGetBlockLoc(distrb, abs(nwBlock), dstProc, &
                                             dstLocalID, errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error finding WestNorthWest block location')
               return
            endif
         else
            dstProc = 0
            dstLocalID = 0
         endif

         if (dstProc /= srcProc) then
            call POP_HaloMsgCreate(halo, iblock,  srcProc, srcLocalID, &
                                         nwBlock, dstProc, dstLocalID, &
                                         'north', errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'POP_HaloCreate: error creating WestNorthWest message')
               return
            endif
         endif

      endif

   end do msgConfigLoop

!-----------------------------------------------------------------------
!EOC

 end function POP_HaloCreate

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloDestroy
! !INTERFACE:

 subroutine POP_HaloDestroy(halo, errorCode)

! !DESCRIPTION:
!  This routine destroys a halo structure by deallocating all memory
!  associated with the halo and nullifying pointers.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (POP_halo), intent(inout) :: &
      halo          ! boundary structure to be destroyed

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode     ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local status flag for deallocate
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: istat

!-----------------------------------------------------------------------
!
!  reset all scalars
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   halo%communicator   = 0
   halo%numMsgSend     = 0
   halo%numMsgRecv     = 0
   halo%numLocalCopies = 0

!-----------------------------------------------------------------------
!
!  deallocate all pointers
!
!-----------------------------------------------------------------------

   deallocate(halo%recvTask, halo%sendTask, &
              halo%sizeSend, halo%sizeRecv, &
              halo%srcLocalAddr, halo%dstLocalAddr, &
              halo%sendAddr, halo%recvAddr,         &
              stat = istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloDestroy: error deallocating halo')
      return
   endif

!-----------------------------------------------------------------------
!
!  nullify all pointers
!
!-----------------------------------------------------------------------

   nullify(halo%recvTask)
   nullify(halo%sendTask)
   nullify(halo%sizeSend)
   nullify(halo%sizeRecv)
   nullify(halo%srcLocalAddr)
   nullify(halo%dstLocalAddr)
   nullify(halo%sendAddr)
   nullify(halo%recvAddr)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_HaloDestroy

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate2DR8
! !INTERFACE:

 subroutine POP_HaloUpdate2DR8(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (POP_r8), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,n,nmsg,step,           &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

#ifdef _ALTREQ
   integer (POP_i4) :: &
      sndRequest(halo%numMsgSend),  &! MPI request ids
      rcvRequest(halo%numMsgRecv)    ! MPI request ids
#else
   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids
#endif

   real (POP_r8) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

!pw++
   call t_startf("HaloUpdate2DR8")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_POP_r8
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR8)) then
      nxGlobal = size(bufTripoleR8,dim=1)
      bufTripoleR8 = fill
   endif

#ifndef _ALTREQ
!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate2DR8: error allocating req arrays')
      return
   endif
#endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecvR8(1:msgSize,nmsg), msgSize, POP_mpiR8, &
                     halo%recvTask(nmsg),                           &
                     POP_mpitagHalo + halo%recvTask(nmsg),          &
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

      msgSize = halo%sizeSend(nmsg)
      call MPI_ISEND(bufSendR8(1:msgSize,nmsg), msgSize, POP_mpiR8, &
                     halo%sendTask(nmsg),                           &
                     POP_mpitagHalo + POP_myTask,                   &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

#ifndef _WAITANY_2D
   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)
#endif

#ifdef _WAITANY_2D
   do step=1,halo%numMsgRecv
      nmsg = -1
      call MPI_WAITANY(halo%numMsgRecv, rcvRequest, nmsg, &
                       MPI_STATUS_IGNORE, ierr)

      if ((nmsg < 1) .or. (nmsg > halo%numMsgRecv)) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate2DR8: error waiting for messages')
         return
      endif
#else
   do nmsg=1,halo%numMsgRecv
#endif
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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate2DR8: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripoleR8(i   ,POP_haloWidth+1)
            x2 = bufTripoleR8(iDst,POP_haloWidth+1)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripoleR8(i   ,POP_haloWidth+1) = isign*sign(xavg, x2)
            bufTripoleR8(iDst,POP_haloWidth+1) = isign*sign(xavg, x1)
         end do
         bufTripoleR8(nxGlobal,POP_haloWidth+1) = isign* &
         bufTripoleR8(nxGlobal,POP_haloWidth+1)

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripoleR8(i   ,POP_haloWidth+1)
            x2 = bufTripoleR8(iDst,POP_haloWidth+1)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripoleR8(i   ,POP_haloWidth+1) = isign*sign(xavg, x2)
            bufTripoleR8(iDst,POP_haloWidth+1) = isign*sign(xavg, x1)
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate2DR8: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

#ifndef _ALTREQ
   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate2DR8: error deallocating req arrays')
      return
   endif
#endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate2DR8")
!pw--

 end subroutine POP_HaloUpdate2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate2DR4
! !INTERFACE:

 subroutine POP_HaloUpdate2DR4(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (POP_r4), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,n,nmsg,step,           &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   real (POP_r4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

!pw++
   call t_startf("HaloUpdate2DR4")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_POP_r4
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR4)) then
      nxGlobal = size(bufTripoleR4,dim=1)
      bufTripoleR4 = fill
   endif

!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate2DR4: error allocating req arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecvR4(1:msgSize,nmsg), msgSize, POP_mpiR4, &
                     halo%recvTask(nmsg),                           &
                     POP_mpitagHalo + halo%recvTask(nmsg),          &
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

      msgSize = halo%sizeSend(nmsg)
      call MPI_ISEND(bufSendR4(1:msgSize,nmsg), msgSize, POP_mpiR4, &
                     halo%sendTask(nmsg),                           &
                     POP_mpitagHalo + POP_myTask,                   &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)

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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate2DR4: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripoleR4(i   ,POP_haloWidth+1)
            x2 = bufTripoleR4(iDst,POP_haloWidth+1)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripoleR4(i   ,POP_haloWidth+1) = isign*sign(xavg, x2)
            bufTripoleR4(iDst,POP_haloWidth+1) = isign*sign(xavg, x1)
         end do
         bufTripoleR4(nxGlobal,POP_haloWidth+1) = isign* &
         bufTripoleR4(nxGlobal,POP_haloWidth+1)

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripoleR4(i   ,POP_haloWidth+1)
            x2 = bufTripoleR4(iDst,POP_haloWidth+1)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripoleR4(i   ,POP_haloWidth+1) = isign*sign(xavg, x2)
            bufTripoleR4(iDst,POP_haloWidth+1) = isign*sign(xavg, x1)
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate2DR4: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate2DR4: error deallocating req arrays')
      return
   endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate2DR4")
!pw--

 end subroutine POP_HaloUpdate2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate2DI4
! !INTERFACE:

 subroutine POP_HaloUpdate2DI4(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   integer (POP_i4), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,n,nmsg,step,           &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (POP_i4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

!pw++
   call t_startf("HaloUpdate2DI4")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0_POP_i4
   endif

   nxGlobal = 0
   if (allocated(bufTripoleI4)) then
      nxGlobal = size(bufTripoleI4,dim=1)
      bufTripoleI4 = fill
   endif

!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate2DI4: error allocating req arrays')
      return
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecvI4(1:msgSize,nmsg), msgSize, MPI_INTEGER, &
                     halo%recvTask(nmsg),                             &
                     POP_mpitagHalo + halo%recvTask(nmsg),            &
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

      msgSize = halo%sizeSend(nmsg)
      call MPI_ISEND(bufSendI4(1:msgSize,nmsg), msgSize, MPI_INTEGER, &
                     halo%sendTask(nmsg),                             &
                     POP_mpitagHalo + POP_myTask,                     &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)

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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate2DI4: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripoleI4(i   ,POP_haloWidth+1)
            x2 = bufTripoleI4(iDst,POP_haloWidth+1)
            xavg = nint(0.5_POP_r8*(abs(x1) + abs(x2)))
            bufTripoleI4(i   ,POP_haloWidth+1) = isign*sign(xavg, x2)
            bufTripoleI4(iDst,POP_haloWidth+1) = isign*sign(xavg, x1)
         end do
         bufTripoleI4(nxGlobal,POP_haloWidth+1) = isign* &
         bufTripoleI4(nxGlobal,POP_haloWidth+1)

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripoleI4(i   ,POP_haloWidth+1)
            x2 = bufTripoleI4(iDst,POP_haloWidth+1)
            xavg = nint(0.5_POP_r8*(abs(x1) + abs(x2)))
            bufTripoleI4(i   ,POP_haloWidth+1) = isign*sign(xavg, x2)
            bufTripoleI4(iDst,POP_haloWidth+1) = isign*sign(xavg, x1)
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate2DI4: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate2DI4: error deallocating req arrays')
      return
   endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate2DI4")
!pw--

 end subroutine POP_HaloUpdate2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate3DR8
! !INTERFACE:

 subroutine POP_HaloUpdate3DR8(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (POP_r8), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,k,n,nmsg,step,         &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of array in 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

#ifdef _ALTREQ
   integer (POP_i4) :: &
      sndRequest(halo%numMsgSend),  &! MPI request ids
      rcvRequest(halo%numMsgRecv)    ! MPI request ids
#else
   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids
#endif

   real (POP_r8) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

#ifdef _SAVE_3D_BUFFERS
   real (POP_r8), dimension(:,:), allocatable, save :: &
      bufSend, bufRecv            ! 3d send,recv buffers

   real (POP_r8), dimension(:,:,:), allocatable, save :: &
      bufTripole                  ! 3d tripole buffer

   integer (POP_i4), save ::  &
      save_nxGlobal,          &
      save_nz_a,              &
      save_nz_b,              &
      save_nz_c,              &
      save_bufSizeSend,       &
      save_bufSizeRecv,       &
      save_numMsgSend,        &
      save_numMsgRecv

   logical (POP_logical) ::   &
      do_allocate       ! flag used to control alloc of 3D buffers
#else
   real (POP_r8), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 3d send,recv buffers

   real (POP_r8), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer
#endif

!pw++
   call t_startf("HaloUpdate3DR8")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_POP_r8
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR8)) nxGlobal = size(bufTripoleR8,dim=1)

#ifndef _ALTREQ
!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR8: error allocating req arrays')
      return
   endif
#endif

!-----------------------------------------------------------------------
!
!  allocate 3D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)

#ifdef _SAVE_3D_BUFFERS
   do_allocate = .false.
   if (.not. allocated(bufSend)) then
      do_allocate = .true.
   else
      if ((save_nz_a        /= nz              ) .or. &
          (save_bufSizeSend /= bufSizeSend     ) .or. &
          (save_numMsgSend  /= halo%numMsgSend )) then

         deallocate(bufSend, stat=ierr)

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate3DR8: error deallocating bufSend')
            return
         endif

         do_allocate = .true.
      endif
   endif

   if (do_allocate) then
      allocate(bufSend(bufSizeSend*nz, halo%numMsgSend), &
               stat=ierr)
      save_nz_a        = nz
      save_bufSizeSend = bufSizeSend
      save_numMsgSend  = halo%numMsgSend

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: error allocating bufSend')
         return
      endif
   endif

   do_allocate = .false.
   if (.not. allocated(bufRecv)) then
      do_allocate = .true.
   else
      if ((save_nz_b        /= nz              ) .or. &
          (save_bufSizeRecv /= bufSizeRecv     ) .or. &
          (save_numMsgRecv  /= halo%numMsgRecv )) then

         deallocate(bufRecv, stat=ierr)

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate3DR8: error deallocating bufRecv')
            return
         endif

         do_allocate = .true.
      endif
   endif

   if (do_allocate) then
      allocate(bufRecv(bufSizeRecv*nz, halo%numMsgRecv), &
               stat=ierr)
      save_nz_b        = nz
      save_bufSizeRecv = bufSizeRecv
      save_numMsgRecv  = halo%numMsgRecv

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: error allocating bufRecv')
         return
      endif
   endif

   do_allocate = .false.
   if (.not. allocated(bufTripole)) then
      do_allocate = .true.
   else
      if ((save_nz_c     /= nz       ) .or. &
          (save_nxGlobal /= nxGlobal )) then
         deallocate(bufTripole, stat=ierr)

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate3DR8: error deallocating bufTripole')
            return
         endif

         do_allocate = .true.
      endif
   endif

   if (do_allocate) then
      if (nxGlobal > 0) then
         allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz), &
                  stat=ierr)

         save_nz_c     = nz
         save_nxGlobal = nxGlobal

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate3DR8: error allocating bufTripole')
            return
         endif
      endif
   endif

   if (nxGlobal > 0) then
      bufTripole = fill
   endif
#else
   allocate(bufSend(bufSizeSend*nz, halo%numMsgSend), &
            bufRecv(bufSizeRecv*nz, halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR8: error allocating buffers')
      return
   endif

   if (nxGlobal > 0) then
      allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz), &
               stat=ierr)
      bufTripole = fill

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: error allocating buffers')
         return
      endif
   endif
#endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = nz*halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecv(1:msgSize,nmsg), msgSize, POP_mpiR8,   &
                     halo%recvTask(nmsg),                           &
                     POP_mpitagHalo + halo%recvTask(nmsg),          &
                     halo%communicator, rcvRequest(nmsg), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(nmsg,i,n,iSrc,jSrc,srcBlock,k)
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
   end do
   !$OMP END PARALLEL DO

   do nmsg=1,halo%numMsgSend

      msgSize = nz*halo%sizeSend(nmsg)
      call MPI_ISEND(bufSend(1:msgSize,nmsg), msgSize, POP_mpiR8, &
                     halo%sendTask(nmsg),                         &
                     POP_mpitagHalo + POP_myTask,                 &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

#ifndef _WAITANY_3D
   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)
#endif

#ifdef _WAITANY_3D
   do step=1,halo%numMsgRecv
      nmsg = -1
      call MPI_WAITANY(halo%numMsgRecv, rcvRequest, nmsg, &
                       MPI_STATUS_IGNORE, ierr)

      if ((nmsg < 1) .or. (nmsg > halo%numMsgRecv)) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: error waiting for messages')
         return
      endif
#else
   !$OMP PARALLEL DO PRIVATE(nmsg,i,n,iDst,jDst,dstBlock,k)
   do nmsg=1,halo%numMsgRecv
#endif
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
#ifndef _WAITANY_3D
   !$OMP END PARALLEL DO
#endif

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k)
            x2 = bufTripole(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         bufTripole(nxGlobal,POP_haloWidth+1,k) = isign* &
         bufTripole(nxGlobal,POP_haloWidth+1,k)
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k)
            x2 = bufTripole(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

#ifndef _ALTREQ
   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR8: error deallocating req arrays')
      return
   endif
#endif

#ifndef _SAVE_3D_BUFFERS
   deallocate(bufSend, bufRecv, stat=ierr)
   if ((ierr == 0) .and. (allocated(bufTripole))) &
      deallocate(bufTripole, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR8: error deallocating 3d buffers')
      return
   endif
#endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate3DR8")
!pw--

 end subroutine POP_HaloUpdate3DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate3DR4
! !INTERFACE:

 subroutine POP_HaloUpdate3DR4(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (POP_r4), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,k,n,nmsg,step,         &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of array in 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   real (POP_r4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (POP_r4), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 3d send,recv buffers

   real (POP_r4), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

!pw++
   call t_startf("HaloUpdate3DR4")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_POP_r4
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR4)) nxGlobal = size(bufTripoleR4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR4: error allocating req arrays')
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
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR4: error allocating buffers')
      return
   endif

   if (nxGlobal > 0) then
      allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz), &
               stat=ierr)
      bufTripole = fill

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: error allocating buffers')
         return
      endif
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = nz*halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecv(1:msgSize,nmsg), msgSize, POP_mpiR4,   &
                     halo%recvTask(nmsg),                           &
                     POP_mpitagHalo + halo%recvTask(nmsg),          &
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

      msgSize = nz*halo%sizeSend(nmsg)
      call MPI_ISEND(bufSend(1:msgSize,nmsg), msgSize, POP_mpiR4, &
                     halo%sendTask(nmsg),                         &
                     POP_mpitagHalo + POP_myTask,                 &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)

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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k)
            x2 = bufTripole(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         bufTripole(nxGlobal,POP_haloWidth+1,k) = isign* &
         bufTripole(nxGlobal,POP_haloWidth+1,k)
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k)
            x2 = bufTripole(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR4: error deallocating req arrays')
      return
   endif

   deallocate(bufSend, bufRecv, stat=ierr)
   if (allocated(bufTripole)) deallocate(bufTripole, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DR4: error deallocating 3d buffers')
      return
   endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate3DR4")
!pw--

 end subroutine POP_HaloUpdate3DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate3DI4
! !INTERFACE:

 subroutine POP_HaloUpdate3DI4(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   integer (POP_i4), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,k,n,nmsg,step,         &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of array in 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (POP_i4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (POP_i4), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 3d send,recv buffers

   integer (POP_i4), dimension(:,:,:), allocatable :: &
      bufTripole                  ! 3d tripole buffer

!pw++
   call t_startf("HaloUpdate3DI4")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0_POP_i4
   endif

   nxGlobal = 0
   if (allocated(bufTripoleI4)) nxGlobal = size(bufTripoleI4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DI4: error allocating req arrays')
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
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DI4: error allocating buffers')
      return
   endif

   if (nxGlobal > 0) then
      allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz), &
               stat=ierr)
      bufTripole = fill

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: error allocating buffers')
         return
      endif
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = nz*halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecv(1:msgSize,nmsg), msgSize, MPI_INTEGER, &
                     halo%recvTask(nmsg),                           &
                     POP_mpitagHalo + halo%recvTask(nmsg),          &
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

      msgSize = nz*halo%sizeSend(nmsg)
      call MPI_ISEND(bufSend(1:msgSize,nmsg), msgSize, MPI_INTEGER, &
                     halo%sendTask(nmsg),                           &
                     POP_mpitagHalo + POP_myTask,                   &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)

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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DI4: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k)
            x2 = bufTripole(iDst,POP_haloWidth+1,k)
            xavg = nint(0.5_POP_r8*(abs(x1) + abs(x2)))
            bufTripole(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         bufTripole(nxGlobal,POP_haloWidth+1,k) = isign* &
         bufTripole(nxGlobal,POP_haloWidth+1,k)
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k)
            x2 = bufTripole(iDst,POP_haloWidth+1,k)
            xavg = nint(0.5_POP_r8*(abs(x1) + abs(x2)))
            bufTripole(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DI4: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DI4: error deallocating req arrays')
      return
   endif

   deallocate(bufSend, bufRecv, stat=ierr)
   if (allocated(bufTripole)) deallocate(bufTripole, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate3DI4: error deallocating 3d buffers')
      return
   endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate3DI4")
!pw--

 end subroutine POP_HaloUpdate3DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate4DR8
! !INTERFACE:

 subroutine POP_HaloUpdate4DR8(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (POP_r8), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,k,l,n,nmsg,step,       &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of array in 3rd,4th dimensions
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

#ifdef _ALTREQ
   integer (POP_i4) :: &
      sndRequest(halo%numMsgSend),  &! MPI request ids
      rcvRequest(halo%numMsgRecv)    ! MPI request ids
#else
   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids
#endif

   real (POP_r8) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

#ifdef _SAVE_4D_BUFFERS
   real (POP_r8), dimension(:,:), allocatable, save :: &
      bufSend, bufRecv            ! 4d send,recv buffers

   real (POP_r8), dimension(:,:,:,:), allocatable, save :: &
      bufTripole                  ! 4d tripole buffer

   integer (POP_i4), save :: &
      save_nxGlobal,         &
      save_nz_a,             &
      save_nz_b,             &
      save_nz_c,             &
      save_nt_a,             &
      save_nt_b,             &
      save_nt_c,             &
      save_bufSizeSend,      &
      save_bufSizeRecv,      &
      save_numMsgSend,       &
      save_numMsgRecv

   logical (POP_logical) :: &
      do_allocate       ! flag used to control alloc of 4D buffers
#else
   real (POP_r8), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 4d send,recv buffers

   real (POP_r8), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer
#endif

!pw++
   call t_startf("HaloUpdate4DR8")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_POP_r8
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR8)) nxGlobal = size(bufTripoleR8,dim=1)

#ifndef _ALTREQ
!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR8: error allocating req arrays')
      return
   endif
#endif

!-----------------------------------------------------------------------
!
!  allocate 4D buffers
!
!-----------------------------------------------------------------------

   nz = size(array, dim=3)
   nt = size(array, dim=4)

#ifdef _SAVE_4D_BUFFERS
   do_allocate = .false.
   if (.not. allocated(bufSend)) then
      do_allocate = .true.
   else
      if ((save_nz_a        /= nz              ) .or. &
          (save_nt_a        /= nt              ) .or. &
          (save_bufSizeSend /= bufSizeSend     ) .or. &
          (save_numMsgSend  /= halo%numMsgSend )) then

         deallocate(bufSend, stat=ierr)

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate4DR8: error deallocating bufSend')
            return
         endif

         do_allocate = .true.
      endif
   endif

   if (do_allocate) then
      allocate(bufSend(bufSizeSend*nz*nt, halo%numMsgSend),   &
               stat=ierr)
      save_nz_a        = nz
      save_nt_a        = nt
      save_bufSizeSend = bufSizeSend
      save_numMsgSend  = halo%numMsgSend

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR8: error allocating bufSend')
         return
      endif
   endif

   do_allocate = .false.
   if (.not. allocated(bufRecv)) then
      do_allocate = .true.
   else
      if ((save_nz_b        /= nz              ) .or. &
          (save_nt_b        /= nt              ) .or. &
          (save_bufSizeRecv /= bufSizeRecv     ) .or. &
          (save_numMsgRecv  /= halo%numMsgRecv )) then

         deallocate(bufRecv, stat=ierr)

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate4DR8: error deallocating bufRecv')
            return
         endif

         do_allocate = .true.
      endif
   endif

   if (do_allocate) then
      allocate(bufRecv(bufSizeRecv*nz*nt, halo%numMsgRecv),   &
               stat=ierr)
      save_nz_b        = nz
      save_nt_b        = nt
      save_bufSizeRecv = bufSizeRecv
      save_numMsgRecv  = halo%numMsgRecv

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR8: error allocating 4d buffers')
         return
      endif
   endif

   do_allocate = .false.
   if (.not. allocated(bufTripole)) then
      do_allocate = .true.
   else
      if ((save_nz_c     /= nz       ) .or. &
          (save_nt_c     /= nt       ) .or. &
          (save_nxGlobal /= nxGlobal )) then
         deallocate(bufTripole, stat=ierr)

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate4DR8: error deallocating bufTripole')
            return
         endif

         do_allocate = .true.
      endif
   endif

   if (do_allocate) then
      if (nxGlobal > 0) then
         allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz, nt), &
                  stat=ierr)

         save_nz_c     = nz
         save_nt_c     = nt
         save_nxGlobal = nxGlobal

         if (ierr > 0) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloUpdate4DR8: error allocating bufTripole')
            return
         endif
      endif
   endif

   if (nxGlobal > 0) then
      bufTripole = fill
   endif
#else
   allocate(bufSend(bufSizeSend*nz*nt, halo%numMsgSend),   &
            bufRecv(bufSizeRecv*nz*nt, halo%numMsgRecv),   &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR8: error allocating buffers')
      return
   endif

   if (nxGlobal > 0) then
      allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz, nt), &
               stat=ierr)
      bufTripole = fill

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: error allocating buffers')
         return
      endif
   endif
#endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = nz*nt*halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecv(1:msgSize,nmsg), msgSize, POP_mpiR8, &
                     halo%recvTask(nmsg),                         &
                     POP_mpitagHalo + halo%recvTask(nmsg),        &
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

      msgSize = nz*nt*halo%sizeSend(nmsg)
      call MPI_ISEND(bufSend(1:msgSize,nmsg), msgSize, POP_mpiR8, &
                     halo%sendTask(nmsg),                         &
                     POP_mpitagHalo + POP_myTask,                 &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

#ifndef _WAITANY_4D
   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)
#endif

#ifdef _WAITANY_4D
   do step=1,halo%numMsgRecv
      nmsg = -1
      call MPI_WAITANY(halo%numMsgRecv, rcvRequest, nmsg, &
                       MPI_STATUS_IGNORE, ierr)

      if ((nmsg < 1) .or. (nmsg > halo%numMsgRecv)) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR8: error waiting for messages')
         return
      endif
#else
   do nmsg=1,halo%numMsgRecv
#endif
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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR8: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do l=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k,l)
            x2 = bufTripole(iDst,POP_haloWidth+1,k,l)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k,l) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k,l) = isign*sign(xavg, x1)
         end do
         bufTripole(nxGlobal,POP_haloWidth+1,k,l) = isign* &
         bufTripole(nxGlobal,POP_haloWidth+1,k,l)
         end do
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do l=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k,l)
            x2 = bufTripole(iDst,POP_haloWidth+1,k,l)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k,l) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k,l) = isign*sign(xavg, x1)
         end do
         end do
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR8: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

#ifndef _ALTREQ
   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR8: error deallocating req arrays')
      return
   endif
#endif

#ifndef _SAVE_4D_BUFFERS
   deallocate(bufSend, bufRecv, stat=ierr)
   if ((ierr == 0) .and. (allocated(bufTripole))) &
      deallocate(bufTripole, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR8: error deallocating 4d buffers')
      return
   endif
#endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate4DR8")
!pw--

 end subroutine POP_HaloUpdate4DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate4DR4
! !INTERFACE:

 subroutine POP_HaloUpdate4DR4(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   real (POP_r4), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,k,l,n,nmsg,step,       &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of array in 3rd,4th dimensions
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   real (POP_r4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (POP_r4), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 4d send,recv buffers

   real (POP_r4), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

!pw++
   call t_startf("HaloUpdate4DR4")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0.0_POP_r4
   endif

   nxGlobal = 0
   if (allocated(bufTripoleR4)) nxGlobal = size(bufTripoleR4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR4: error allocating req arrays')
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
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR4: error allocating buffers')
      return
   endif

   if (nxGlobal > 0) then
      allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz, nt), &
               stat=ierr)
      bufTripole = fill

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: error allocating buffers')
         return
      endif
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = nz*nt*halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecv(1:msgSize,nmsg), msgSize, POP_mpiR4, &
                     halo%recvTask(nmsg),                         &
                     POP_mpitagHalo + halo%recvTask(nmsg),        &
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

      msgSize = nz*nt*halo%sizeSend(nmsg)
      call MPI_ISEND(bufSend(1:msgSize,nmsg), msgSize, POP_mpiR4, &
                     halo%sendTask(nmsg),                         &
                     POP_mpitagHalo + POP_myTask,                 &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)

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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR4: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do l=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k,l)
            x2 = bufTripole(iDst,POP_haloWidth+1,k,l)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k,l) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k,l) = isign*sign(xavg, x1)
         end do
         bufTripole(nxGlobal,POP_haloWidth+1,k,l) = isign* &
         bufTripole(nxGlobal,POP_haloWidth+1,k,l)
         end do
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do l=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k,l)
            x2 = bufTripole(iDst,POP_haloWidth+1,k,l)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole(i   ,POP_haloWidth+1,k,l) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k,l) = isign*sign(xavg, x1)
         end do
         end do
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR4: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR4: error deallocating req arrays')
      return
   endif

   deallocate(bufSend, bufRecv, stat=ierr)
   if (allocated(bufTripole)) deallocate(bufTripole, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DR4: error deallocating 4d buffers')
      return
   endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate4DR4")
!pw--

 end subroutine POP_HaloUpdate4DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloUpdate4DI4
! !INTERFACE:

 subroutine POP_HaloUpdate4DI4(array, halo,                    &
                               fieldLoc, fieldKind, errorCode, &
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

   type (POP_halo), intent(in) :: &
      halo                 ! precomputed halo structure containing all
                           !  information needed for halo update

   character (*), intent(in) :: &
      fieldKind,          &! id for type of field (scalar, vector, angle)
      fieldLoc             ! id for location on horizontal grid
                           !  (center, NEcorner, Nface, Eface)

   integer (POP_i4), intent(in), optional :: &
      fillValue            ! optional value to put in ghost cells
                           !  where neighbor points are unknown
                           !  (e.g. eliminated land blocks or
                           !   closed boundaries)

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:,:,:), intent(inout) :: &
      array                ! array containing field for which halo
                           ! needs to be updated

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::           &
      i,j,k,l,n,nmsg,step,       &! dummy loop indices
      ierr,                      &! error or status flag for MPI,alloc
      msgSize,                   &! size of an individual message
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of array in 3rd,4th dimensions
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest,      &! MPI request ids
      rcvRequest        ! MPI request ids

   integer (POP_i4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (POP_i4), dimension(:,:), allocatable :: &
      bufSend, bufRecv            ! 4d send,recv buffers

   integer (POP_i4), dimension(:,:,:,:), allocatable :: &
      bufTripole                  ! 4d tripole buffer

!pw++
   call t_startf("HaloUpdate4DI4")
!pw--
!-----------------------------------------------------------------------
!
!  initialize error code and fill value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(fillValue)) then
      fill = fillValue
   else
      fill = 0_POP_i4
   endif

   nxGlobal = 0
   if (allocated(bufTripoleI4)) nxGlobal = size(bufTripoleI4,dim=1)

!-----------------------------------------------------------------------
!
!  allocate request arrays for messages
!
!-----------------------------------------------------------------------

   allocate(sndRequest(halo%numMsgSend), &
            rcvRequest(halo%numMsgRecv), &
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DI4: error allocating req arrays')
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
            stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DI4: error allocating buffers')
      return
   endif

   if (nxGlobal > 0) then
      allocate(bufTripole(nxGlobal, POP_haloWidth+1, nz, nt), &
               stat=ierr)
      bufTripole = fill

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: error allocating buffers')
         return
      endif
   endif

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do nmsg=1,halo%numMsgRecv

      msgSize = nz*nt*halo%sizeRecv(nmsg)
      call MPI_IRECV(bufRecv(1:msgSize,nmsg), msgSize, MPI_INTEGER, &
                     halo%recvTask(nmsg),                           &
                     POP_mpitagHalo + halo%recvTask(nmsg),          &
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

      msgSize = nz*nt*halo%sizeSend(nmsg)
      call MPI_ISEND(bufSend(1:msgSize,nmsg), msgSize, MPI_INTEGER, &
                     halo%sendTask(nmsg),                           &
                     POP_mpitagHalo + POP_myTask,                   &
                     halo%communicator, sndRequest(nmsg), ierr)
   end do

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

   call MPI_WAITALL(halo%numMsgRecv, rcvRequest, MPI_STATUSES_IGNORE, &
                    ierr)

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
!  bufTripole array contains the top haloWidth+1 rows of physical
!    domain for entire (global) top row 
!
!-----------------------------------------------------------------------

   if (nxGlobal > 0) then

      select case (fieldKind)
      case (POP_fieldKindScalar)
         isign =  1
      case (POP_fieldKindVector)
         isign = -1
      case (POP_fieldKindAngle)
         isign = -1
      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DI4: Unknown field kind')
      end select

      select case (fieldLoc)
      case (POP_gridHorzLocCenter)   ! cell center location

         ioffset = 0
         joffset = 0

      case (POP_gridHorzLocNEcorner)   ! cell corner location

         ioffset = 1
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do l=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k,l)
            x2 = bufTripole(iDst,POP_haloWidth+1,k,l)
            xavg = nint(0.5_POP_r8*(abs(x1) + abs(x2)))
            bufTripole(i   ,POP_haloWidth+1,k,l) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k,l) = isign*sign(xavg, x1)
         end do
         bufTripole(nxGlobal,POP_haloWidth+1,k,l) = isign* &
         bufTripole(nxGlobal,POP_haloWidth+1,k,l)
         end do
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value

         do l=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole(i   ,POP_haloWidth+1,k,l)
            x2 = bufTripole(iDst,POP_haloWidth+1,k,l)
            xavg = nint(0.5_POP_r8*(abs(x1) + abs(x2)))
            bufTripole(i   ,POP_haloWidth+1,k,l) = isign*sign(xavg, x2)
            bufTripole(iDst,POP_haloWidth+1,k,l) = isign*sign(xavg, x1)
         end do
         end do
         end do

      case default
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DI4: Unknown field location')
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

            if (jSrc <= POP_haloWidth+1) then
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

   call MPI_WAITALL(halo%numMsgSend, sndRequest, MPI_STATUSES_IGNORE, &
                    ierr)

   deallocate(sndRequest, rcvRequest, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DI4: error deallocating req arrays')
      return
   endif

   deallocate(bufSend, bufRecv, stat=ierr)
   if (allocated(bufTripole)) deallocate(bufTripole, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloUpdate4DI4: error deallocating 4d buffers')
      return
   endif

!-----------------------------------------------------------------------
!EOC
!pw++
   call t_stopf("HaloUpdate4DI4")
!pw--

 end subroutine POP_HaloUpdate4DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloIncrementMsgCount
! !INTERFACE:

   subroutine POP_HaloIncrementMsgCount(sndCounter, rcvCounter,    &
                                        srcProc, dstProc, msgSize, &
                                        errorCode)

! !DESCRIPTION:
!  This is a utility routine to increment the arrays for counting
!  whether messages are required.  It checks the source and destination
!  task to see whether the current task needs to send, receive or
!  copy messages to fill halo regions (ghost cells).

! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcProc,               &! source processor for communication
      dstProc,               &! destination processor for communication
      msgSize                 ! number of words for this message

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:), intent(inout) :: &
      sndCounter,       &! array for counting messages to be sent
      rcvCounter         ! array for counting messages to be received

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  error check
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (srcProc < 0 .or. dstProc < 0 .or. &
       srcProc > size(sndCounter)   .or. &
       dstProc > size(rcvCounter)) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloIncrementMsgCount: invalid processor number')
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

   if (srcProc == POP_myTask + 1) sndCounter(dstProc) = &
                                  sndCounter(dstProc) + msgSize

!-----------------------------------------------------------------------
!
!  if the current processor is the destination, must receive data 
!  local copy if dstProc = srcProc
!
!-----------------------------------------------------------------------

   if (dstProc == POP_myTask + 1) then

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

   end subroutine POP_HaloIncrementMsgCount

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloMsgCreate
! !INTERFACE:

   subroutine POP_HaloMsgCreate(halo, srcBlock, srcProc, srcLocalID, &
                                      dstBlock, dstProc, dstLocalID, &
                                      direction, errorCode)

! !DESCRIPTION:
!  This is a utility routine to determine the required address and
!  message information for a particular pair of blocks.

! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcBlock,   dstBlock,   & ! source,destination block id
      srcProc,    dstProc,    & ! source,destination processor location
      srcLocalID, dstLocalID    ! source,destination local index

   character (*), intent(in) :: &
      direction              ! direction of neighbor block
                             !  (north,south,east,west,
                             !   and NE, NW, SE, SW)

! !INPUT/OUTPUT PARAMETERS:

   type (POP_halo), intent(inout) :: &
      halo                   ! data structure containing halo info

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      msgIndx,               &! message counter and index into msg array
      blockIndx,             &! block counter and index into msg array
      bufSize,               &! size of message buffer
      ibSrc, ieSrc, jbSrc, jeSrc, &! phys domain info for source block
      ibDst, ieDst, jbDst, jeDst, &! phys domain info for dest   block
      nxGlobal,              &! size of global domain in e-w direction
      i,j,n                   ! dummy loop index

   integer (POP_i4), dimension(:), pointer :: &
      iGlobal                 ! global i index for location in tripole

!-----------------------------------------------------------------------
!
!  initialize
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

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

   if (srcProc == POP_myTask+1 .or. dstProc == POP_myTask+1) then

      if (srcBlock >= 0 .and. dstBlock >= 0) then
         call POP_BlocksGetBlockInfo(srcBlock, errorCode,  &
                                     ib=ibSrc, ie=ieSrc,   &
                                     jb=jbSrc, je=jeSrc)
      else ! tripole - need iGlobal info
         call POP_BlocksGetBlockInfo(abs(srcBlock), errorCode,  &
                                     ib=ibSrc, ie=ieSrc,        &
                                     jb=jbSrc, je=jeSrc,        &
                                     iGlobal=iGlobal)

      endif

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: error getting source block info')
         return
      endif

      if (dstBlock /= 0) then
         call POP_BlocksGetBlockInfo(abs(dstBlock), errorCode,  &
                                     ib=ibDst, ie=ieDst,   &
                                     jb=jbDst, je=jeDst)
      endif

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: error getting dest block info')
         return
      endif

   endif

!-----------------------------------------------------------------------
!
!  if both blocks are local, create a local copy to fill halo
!
!-----------------------------------------------------------------------

   if (srcProc == POP_myTask+1 .and. &
       dstProc == POP_myTask+1) then   

      !*** compute addresses based on direction

      msgIndx = halo%numLocalCopies

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: msg count > array size')
         return
      endif

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

         do j=1,jeSrc-jbSrc+1
         do i=1,POP_haloWidth

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = ieSrc - POP_haloWidth + i
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
         do i=1,POP_haloWidth

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

            do j=1,POP_haloWidth
            do i=1,ieSrc-ibSrc+1

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc - POP_haloWidth + j
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

            if (jeSrc - jbSrc + 1 < POP_haloWidth + 1) then
               call POP_ErrorSet(errorCode, &
               'POP_HaloMsgCreate: not enough points in block for tripole')
               return
            endif 

            do j=1,POP_haloWidth+1
            do i=1,ieSrc-ibSrc+1

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc-1-POP_haloWidth+j
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

            do j=1,POP_haloWidth+1
            do i=1,ieSrc+POP_haloWidth

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = nxGlobal - iGlobal(i) + 1
               halo%srcLocalAddr(2,msgIndx) = POP_haloWidth + 3 - j
               halo%srcLocalAddr(3,msgIndx) = -srcLocalID

               halo%dstLocalAddr(1,msgIndx) = i
               halo%dstLocalAddr(2,msgIndx) = jeSrc + j - 1
               halo%dstLocalAddr(3,msgIndx) = dstLocalID

            end do
            end do

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         do j=1,POP_haloWidth
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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ieSrc - POP_haloWidth + i
               halo%srcLocalAddr(2,msgIndx) = jeSrc - POP_haloWidth + j
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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

               msgIndx = msgIndx + 1

               halo%srcLocalAddr(1,msgIndx) = ibSrc + i - 1
               halo%srcLocalAddr(2,msgIndx) = jeSrc - POP_haloWidth + j
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

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

            msgIndx = msgIndx + 1

            halo%srcLocalAddr(1,msgIndx) = ieSrc - POP_haloWidth + i
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

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

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

         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: unknown direction local copy')
         return

      end select

      halo%numLocalCopies = msgIndx

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: msg count > array size')
         return
      endif

!-----------------------------------------------------------------------
!
!  if dest block is local and source block does not exist, create a 
!  local copy to fill halo with a fill value
!
!-----------------------------------------------------------------------

   else if (srcProc == 0 .and. dstProc == POP_myTask+1) then   

      msgIndx = halo%numLocalCopies

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: msg count > array size')
         return
      endif

      !*** compute addresses based on direction

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

         do j=1,jeSrc-jbSrc+1
         do i=1,POP_haloWidth

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
         do i=1,POP_haloWidth

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

            do j=1,POP_haloWidth
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

         do j=1,POP_haloWidth
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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

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

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

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

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

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

         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: unknown direction local copy')
         return

      end select

      halo%numLocalCopies = msgIndx

      if (msgIndx > size(halo%srcLocalAddr,dim=2) .or. &
          msgIndx > size(halo%dstLocalAddr,dim=2)) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: msg count > array size')
         return
      endif

!-----------------------------------------------------------------------
!
!  if source block local and dest block remote, send a message
!
!-----------------------------------------------------------------------

   else if (srcProc == POP_myTask+1 .and. &
            dstProc /= POP_myTask+1 .and. dstProc > 0) then

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
         do i=1,POP_haloWidth

            bufSize = bufSize + 1

            halo%sendAddr(1,bufSize,msgIndx) = ieSrc - POP_haloWidth + i
            halo%sendAddr(2,bufSize,msgIndx) = jbSrc + j - 1
            halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

         end do
         end do

         halo%sizeSend(msgIndx) = bufSize

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

         do j=1,jeSrc-jbSrc+1
         do i=1,POP_haloWidth

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

            do j=1,POP_haloWidth
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx) = ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx) = jeSrc-POP_haloWidth+j
               halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         else 

            !*** tripole block - send top three rows of phys domain

            do j=1,POP_haloWidth+1
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx)=ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx)=jeSrc-POP_haloWidth+j-1
               halo%sendAddr(3,bufSize,msgIndx)=srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         do j=1,POP_haloWidth
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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx) = ieSrc-POP_haloWidth+i
               halo%sendAddr(2,bufSize,msgIndx) = jeSrc-POP_haloWidth+j
               halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         else 

            !*** tripole block - send top three rows of phys domain

            do j=1,POP_haloWidth+1
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx)=ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx)=jeSrc-POP_haloWidth+j-1
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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx) = ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx) = jeSrc-POP_haloWidth+j
               halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         else 

            !*** tripole block - send top three rows of phys domain

            do j=1,POP_haloWidth+1
            do i=1,ieSrc-ibSrc+1

               bufSize = bufSize + 1

               halo%sendAddr(1,bufSize,msgIndx)=ibSrc + i - 1
               halo%sendAddr(2,bufSize,msgIndx)=jeSrc-POP_haloWidth+j-1
               halo%sendAddr(3,bufSize,msgIndx)=srcLocalID

            end do
            end do

            halo%sizeSend(msgIndx) = bufSize

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

            bufSize = bufSize + 1

            halo%sendAddr(1,bufSize,msgIndx) = ieSrc - POP_haloWidth + i
            halo%sendAddr(2,bufSize,msgIndx) = jbSrc + j - 1
            halo%sendAddr(3,bufSize,msgIndx) = srcLocalID

         end do
         end do

         halo%sizeSend(msgIndx) = bufSize

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

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

   else if (dstProc == POP_myTask+1 .and. &
            srcProc /= POP_myTask+1 .and. srcProc > 0) then

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
         do i=1,POP_haloWidth

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
         do i=1,POP_haloWidth

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

            do j=1,POP_haloWidth
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

            do j=1,POP_haloWidth+1
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

         do j=1,POP_haloWidth
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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = i
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         else

            !*** tripole block - receive into tripole buffer

            do j=1,POP_haloWidth+1
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

            do j=1,POP_haloWidth
            do i=1,POP_haloWidth

               bufSize = bufSize + 1

               halo%recvAddr(1,bufSize,msgIndx) = ieDst + i
               halo%recvAddr(2,bufSize,msgIndx) = j
               halo%recvAddr(3,bufSize,msgIndx) = dstLocalID

            end do
            end do

            halo%sizeRecv(msgIndx) = bufSize

         else

            !*** tripole block - receive into tripole buffer

            do j=1,POP_haloWidth+1
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

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

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

         do j=1,POP_haloWidth
         do i=1,POP_haloWidth

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

   end subroutine POP_HaloMsgCreate

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloPrintStats
! !INTERFACE:

   subroutine POP_HaloPrintStats(halo, distrb, errorCode)

! !DESCRIPTION:
!  This routine compiles some message statistics for a given halo
!  updates and writes them to stdout.  The routine only outputs
!  information for 2D halos of type r8.  Other statistics can
!  be obtained by scaling 2D r8 values accordingly.
!
! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   type (POP_halo), intent(in)   :: &
      halo               ! defined halo for which stats requested

   type (POP_distrb), intent(in) :: &
      distrb             ! associated block distribution for halo

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
 
   integer (POP_i4) ::              &
      bytesSend,     bytesRecv,     &
      maxBytesSend,  maxBytesRecv,  &
      minBytesSend,  minBytesRecv,  & 
      minNumMsgRecv, maxNumMsgRecv, &
      numProcs, n

   real (POP_r8) ::  &
      avgBytesSend, avgBytesRecv

!-----------------------------------------------------------------------
!
!  initialize num procs
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(distrb,errorCode,numProcs = numProcs)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting num procs')
      return
   endif

!-----------------------------------------------------------------------
!
!  determine number of messages received
!
!-----------------------------------------------------------------------

   maxNumMsgRecv = POP_GlobalMaxval(halo%numMsgRecv, distrb, errorCode) 
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting max number of messages')
      return
   endif

   minNumMsgRecv = POP_GlobalMinval(halo%numMsgRecv, distrb, errorCode) 
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting min number of messages')
      return
   endif

!-----------------------------------------------------------------------
!
!  compute local number of bytes sent or received, then determine
!  global statistics
!
!-----------------------------------------------------------------------

   bytesSend = 0
   do n=1,halo%numMsgSend
      bytesSend = bytesSend + 8*halo%sizeSend(n)
   end do

   bytesRecv = 0
   do n=1,halo%numMsgRecv
      bytesRecv = bytesRecv + 8*halo%sizeRecv(n)
   end do

   maxBytesSend = POP_GlobalMaxval(bytesSend, distrb, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting max bytes sent')
      return
   endif

   minBytesSend = POP_GlobalMinval(bytesSend, distrb, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting min bytes sent')
      return
   endif

   avgBytesSend = POP_GlobalSum(real(bytesSend), distrb, errorCode)/ &
	          real(numProcs)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting avg bytes sent')
      return
   endif

   maxBytesRecv = POP_GlobalMaxval(bytesRecv, distrb, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting max bytes received')
      return
   endif

   minBytesRecv = POP_GlobalMinval(bytesRecv, distrb, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting min bytes received')
      return
   endif

   avgBytesRecv = POP_GlobalSum(real(bytesRecv), distrb, errorCode)/ &
	          real(numProcs)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloPrintStats: error getting avg bytes received')
      return
   endif

!-----------------------------------------------------------------------
!
!  output data to stdout
!
!-----------------------------------------------------------------------

   if (POP_mastertask == POP_mytask) then 

      write(POP_stdout,'(a30,2(i13))') 'bufSize{Recv,Send} [words] : ',&
                                        bufSizeRecv, bufSizeSend
      write(POP_stdout,'(a30,2(i13))') 'num messages: {min,max}: ',    &
                                    minNumMsgRecv, maxNumMsgRecv
      write(POP_stdout,'(a45,2(i13),2x,e10.3)')                        &
                     'Bytes RECV for 2D bndy exch {min,max,avg}: ',    &
                      minBytesRecv, maxBytesRecv, avgBytesRecv
      write(POP_stdout,'(a45,2(i13),2x,e10.3)')                        &
                     'Bytes SEND for 2D bndy exch {min,max,avg}: ',    &
                      minBytesSend, maxBytesSend, avgBytesSend

   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine  POP_HaloPrintStats

!***********************************************************************

end module POP_HaloMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
