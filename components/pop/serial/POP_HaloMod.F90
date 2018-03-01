!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: POP_HaloMod

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
!              fixed some bugs including a hard-wired assumption of
!                 halo width 2 and some redundant messaging

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_FieldMod
   use POP_GridHorzMod

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: POP_halo
      integer (POP_i4) ::  &
         communicator,     &! communicator to use for update messages
         numLocalCopies     ! num local copies for halo update

      integer (POP_i4), dimension(:,:), pointer :: &
         srcLocalAddr,     &! src addresses for each local copy
         dstLocalAddr       ! dst addresses for each local copy

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
      istat,                       &! allocate status flag
      numProcs,                    &! num of processors involved
      communicator,                &! communicator for message passing
      iblock,                      &! block counter
      eastBlock, westBlock,        &! block id  east,  west neighbors
      northBlock, southBlock,      &! block id north, south neighbors
      neBlock, nwBlock,            &! block id northeast, northwest nbrs
      seBlock, swBlock,            &! block id southeast, southwest nbrs
      blockSizeX,                  &! size of default phys domain in X
      blockSizeY,                  &! size of default phys domain in Y
      srcProc, dstProc,            &! source, dest processor locations
      srcLocalID, dstLocalID,      &! local block index of src,dst blocks
      eastMsgSize, westMsgSize,    &! nominal sizes for e-w msgs
      northMsgSize, southMsgSize,  &! nominal sizes for n-s msgs
      tripoleMsgSize,              &! tripole message size
      tripoleMsgSizeOut,           &! tripole message size for copy out
      cornerMsgSize, msgSize        ! nominal size for corner msg

   integer (POP_i4), dimension(:), allocatable :: &
      sendCount, recvCount          ! count number of words to each proc

   logical (POP_logical) :: &
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

   blockSizeX  = POP_nxBlock - 2*POP_haloWidth
   blockSizeY  = POP_nyBlock - 2*POP_haloWidth
   eastMsgSize  = POP_haloWidth*blockSizeY
   westMsgSize  = POP_haloWidth*blockSizeY
   southMsgSize = POP_haloWidth*blockSizeX
   northMsgSize = POP_haloWidth*blockSizeX
   cornerMsgSize = POP_haloWidth*POP_haloWidth
   tripoleMsgSize    = (POP_haloWidth + 1)*blockSizeX
   tripoleMsgSizeOut = (POP_haloWidth + 1)*POP_nxBlock

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
               'POP_HaloCreate: error allocating count arrays')
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

      call POP_HaloIncrementMsgCount(sendCount, recvCount,      &
                                     srcProc, dstProc, msgSize, &
                                     errorCode)

      !*** if a tripole boundary block, also create a local
      !*** message into and out of tripole buffer 

      if (tripoleBlock) then
         call POP_HaloIncrementMsgCount(sendCount, recvCount, &
                                     srcProc, srcProc,        &
                                     tripoleMsgSize, errorCode)
      endif

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing tripole copy in count')
         return
      endif

      if (tripoleBlock) then
         call POP_HaloIncrementMsgCount(sendCount, recvCount,  &
                                     srcProc, srcProc,         &
                                     tripoleMsgSizeOut,        &
                                     errorCode)
      endif

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error incrementing tripole copy out count')
         return
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
!  allocate arrays for message information and initialize
!
!-----------------------------------------------------------------------

   allocate(halo%srcLocalAddr(3,halo%numLocalCopies), &
            halo%dstLocalAddr(3,halo%numLocalCopies), &
            stat = istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_HaloCreate: error allocating halo message info arrays')
      return
   endif

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

   halo%numLocalCopies = 0

   msgConfigLoop: do iblock=1,POP_numBlocks

      !***
      !*** find north neighbor block
      !***

      northBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorth,        &
                                      ewBoundaryType, nsBoundaryType, &
                                      errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding north neighbor block')
         return
      endif

      call POP_HaloMsgCreate(halo, distrb, iblock, northBlock, &
                             'north', errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating north message')
         return
      endif

      !*** set tripole flag and add a two copies for inserting
      !*** and extracting info from the tripole buffer

      if (northBlock < 0) then
         tripoleBlock = .true.

         call POP_HaloMsgCreate(halo, distrb, iblock, -iblock, &
                                'north', errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating north tripole copy in')
            return
         endif

         call POP_HaloMsgCreate(halo, distrb, -iblock, iblock, &
                                'north', errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating north tripole copy out')
            return
         endif

      else
         tripoleBlock = .false.
      endif

      !***
      !*** find south neighbor block
      !***

      southBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSouth,        &
                                      ewBoundaryType, nsBoundaryType, &
                                      errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding south neighbor block')
         return
      endif

      call POP_HaloMsgCreate(halo, distrb, iblock, southBlock, &
                             'south', errorCode)
 
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating south message')
         return
      endif

      !***
      !*** find east neighbor block
      !***

      eastBlock = POP_BlocksGetNbrID(iblock, POP_BlocksEast,         &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding east neighbor block')
         return
      endif

      call POP_HaloMsgCreate(halo, distrb, iblock, eastBlock, &
                             'east', errorCode)
 
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating east message')
         return
      endif

      !***
      !*** find west neighbor block
      !***

      westBlock = POP_BlocksGetNbrID(iblock, POP_BlocksWest,         &
                                     ewBoundaryType, nsBoundaryType, &
                                     errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding west neighbor block')
         return
      endif

      call POP_HaloMsgCreate(halo, distrb, iblock, westBlock, &
                             'west', errorCode)
 
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating west message')
         return
      endif

      !***
      !*** find northeast neighbor block
      !***

      neBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorthEast,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding northeast neighbor block')
         return
      endif

      if (neBlock > 0) then
         call POP_HaloMsgCreate(halo, distrb, iblock, neBlock, &
                                'northeast', errorCode)
 
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error creating northeast message')
            return
         endif
      endif

      !***
      !*** find northwest neighbor block
      !***

      nwBlock = POP_BlocksGetNbrID(iblock, POP_BlocksNorthWest,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding northwest neighbor block')
         return
      endif

      if (nwBlock > 0) then
         call POP_HaloMsgCreate(halo, distrb, iblock, nwBlock, &
                                'northwest', errorCode)
 
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_HaloCreate: error creating northwest message')
            return
         endif
      endif

      !***
      !*** find southeast neighbor block
      !***

      seBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSouthEast,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding southeast neighbor block')
         return
      endif

      call POP_HaloMsgCreate(halo, distrb, iblock, seBlock, &
                             'southeast', errorCode)
 
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating southeast message')
         return
      endif

      !***
      !*** find southwest neighbor block
      !***

      swBlock = POP_BlocksGetNbrID(iblock, POP_BlocksSouthWest,    &
                                   ewBoundaryType, nsBoundaryType, &
                                   errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error finding southwest neighbor block')
         return
      endif

      call POP_HaloMsgCreate(halo, distrb, iblock, swBlock, &
                             'southwest', errorCode)
 
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloCreate: error creating southwest message')
         return
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
   halo%numLocalCopies = 0

!-----------------------------------------------------------------------
!
!  deallocate all pointers
!
!-----------------------------------------------------------------------

   deallocate(halo%srcLocalAddr, halo%dstLocalAddr, &
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

   nullify(halo%srcLocalAddr)
   nullify(halo%dstLocalAddr)

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
      i,nmsg,                    &! dummy loop indices
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (POP_r8) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

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

!-----------------------------------------------------------------------
!
!  do local copies to fill halo region
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

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
!EOC

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
      i,nmsg,                    &! dummy loop indices
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (POP_r4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

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
!  do local copies to fill halo region
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

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
!EOC

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
      i,nmsg,                    &! dummy loop indices
      nxGlobal,                  &! global domain size in x (tripole)
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

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
!  do local copies to fill halo region
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripoleI4(i   ,POP_haloWidth+1)
            x2 = bufTripoleI4(iDst,POP_haloWidth+1)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripoleI4(i   ,POP_haloWidth+1)
            x2 = bufTripoleI4(iDst,POP_haloWidth+1)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
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
!EOC

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
      i,k,nmsg,                  &! dummy loop indices
      istat,                     &! status flag for allocate
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (POP_r8) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (POP_r8), dimension(:,:,:), allocatable :: &
      bufTripole3DR8    ! tripole buffer

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

   nz = size(array,dim=3)

   nxGlobal = 0
   if (allocated(bufTripoleR8)) then
      nxGlobal = size(bufTripoleR8,dim=1)
      allocate(bufTripole3DR8(nxGlobal,POP_haloWidth+1,nz), stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: error allocating tripole buffer')
         return 
      endif
      bufTripole3DR8 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies to fill halo region
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
               bufTripole3DR8(iDst,jDst,k) = &
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole3DR8(i   ,POP_haloWidth+1,k)
            x2 = bufTripole3DR8(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole3DR8(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole3DR8(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         bufTripole3DR8(nxGlobal,POP_haloWidth+1,k) = isign* &
         bufTripole3DR8(nxGlobal,POP_haloWidth+1,k)
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole3DR8(i   ,POP_haloWidth+1,k)
            x2 = bufTripole3DR8(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole3DR8(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole3DR8(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
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
                  array(iDst,jDst,k,dstBlock) = isign* &
                                  bufTripole3DR8(iSrc,jSrc,k)
               end do
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!
!  clean up buffers
!
!-----------------------------------------------------------------------

   if (allocated(bufTripole3DR8)) then
      deallocate(bufTripole3DR8, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR8: error deallocating tripole buffer')
         return 
      endif
   endif

!-----------------------------------------------------------------------
!EOC

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
      i,k,nmsg,                  &! dummy loop indices
      istat,                     &! status flag for allocate
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (POP_r4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (POP_r4), dimension(:,:,:), allocatable :: &
      bufTripole3DR4    ! tripole buffer

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

   nz = size(array,dim=3)

   nxGlobal = 0
   if (allocated(bufTripoleR4)) then
      nxGlobal = size(bufTripoleR4,dim=1)
      allocate(bufTripole3DR4(nxGlobal,POP_haloWidth+1,nz), stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: error allocating tripole buffer')
         return 
      endif
      bufTripole3DR4 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies to fill halo region
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
               bufTripole3DR4(iDst,jDst,k) = &
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole3DR4(i   ,POP_haloWidth+1,k)
            x2 = bufTripole3DR4(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole3DR4(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole3DR4(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         bufTripole3DR4(nxGlobal,POP_haloWidth+1,k) = isign* &
         bufTripole3DR4(nxGlobal,POP_haloWidth+1,k)
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole3DR4(i   ,POP_haloWidth+1,k)
            x2 = bufTripole3DR4(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole3DR4(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole3DR4(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
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
                  array(iDst,jDst,k,dstBlock) = isign* &
                                  bufTripole3DR4(iSrc,jSrc,k)
               end do
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!
!  clean up buffers
!
!-----------------------------------------------------------------------

   if (allocated(bufTripole3DR4)) then
      deallocate(bufTripole3DR4, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DR4: error deallocating tripole buffer')
         return 
      endif
   endif

!-----------------------------------------------------------------------
!EOC

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
      i,k,nmsg,                  &! dummy loop indices
      istat,                     &! status flag for allocate
      nxGlobal,                  &! global domain size in x (tripole)
      nz,                        &! size of 3rd dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (POP_i4), dimension(:,:,:), allocatable :: &
      bufTripole3DI4    ! tripole buffer

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

   nz = size(array,dim=3)

   nxGlobal = 0
   if (allocated(bufTripoleI4)) then
      nxGlobal = size(bufTripoleI4,dim=1)
      allocate(bufTripole3DI4(nxGlobal,POP_haloWidth+1,nz), stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DI4: error allocating tripole buffer')
         return 
      endif
      bufTripole3DI4 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies to fill halo region
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
               bufTripole3DI4(iDst,jDst,k) = &
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole3DI4(i   ,POP_haloWidth+1,k)
            x2 = bufTripole3DI4(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole3DI4(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole3DI4(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
         end do
         bufTripole3DI4(nxGlobal,POP_haloWidth+1,k) = isign* &
         bufTripole3DI4(nxGlobal,POP_haloWidth+1,k)
         end do

      case (POP_gridHorzLocEface)   ! cell center location

         ioffset = 1
         joffset = 0

      case (POP_gridHorzLocNface)   ! cell corner (velocity) location

         ioffset = 0
         joffset = 1

         !*** top row is degenerate, so must enforce symmetry
         !***   use average of two degenerate points for value
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole3DI4(i   ,POP_haloWidth+1,k)
            x2 = bufTripole3DI4(iDst,POP_haloWidth+1,k)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole3DI4(i   ,POP_haloWidth+1,k) = isign*sign(xavg, x2)
            bufTripole3DI4(iDst,POP_haloWidth+1,k) = isign*sign(xavg, x1)
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
                  array(iDst,jDst,k,dstBlock) = isign* &
                                  bufTripole3DI4(iSrc,jSrc,k)
               end do
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!
!  clean up buffers
!
!-----------------------------------------------------------------------

   if (allocated(bufTripole3DI4)) then
      deallocate(bufTripole3DI4, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate3DI4: error deallocating tripole buffer')
         return 
      endif
   endif

!-----------------------------------------------------------------------
!EOC

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
      i,k,n,nmsg,                &! dummy loop indices
      istat,                     &! status flag for allocate
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of 3rd, 4th dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (POP_r8) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (POP_r8), dimension(:,:,:,:), allocatable :: &
      bufTripole4DR8    ! tripole buffer

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

   nz = size(array,dim=3)
   nt = size(array,dim=4)

   nxGlobal = 0
   if (allocated(bufTripoleR8)) then
      nxGlobal = size(bufTripoleR8,dim=1)
      allocate(bufTripole4DR8(nxGlobal,POP_haloWidth+1,nz,nt), stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR8: error allocating tripole buffer')
         return 
      endif
      bufTripole4DR8 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies to fill halo region
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
            do n=1,nt
            do k=1,nz
               array(iDst,jDst,k,n,dstBlock) = &
               array(iSrc,jSrc,k,n,srcBlock)
            end do
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do n=1,nt
            do k=1,nz
               bufTripole4DR8(iDst,jDst,k,n) = &
               array(iSrc,jSrc,k,n,srcBlock)
            end do
            end do
         endif
      else if (srcBlock == 0) then
         do n=1,nt
         do k=1,nz
            array(iDst,jDst,k,n,dstBlock) = fill
         end do
         end do
      endif
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do n=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole4DR8(i   ,POP_haloWidth+1,k,n)
            x2 = bufTripole4DR8(iDst,POP_haloWidth+1,k,n)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole4DR8(i   ,POP_haloWidth+1,k,n) = isign*sign(xavg, x2)
            bufTripole4DR8(iDst,POP_haloWidth+1,k,n) = isign*sign(xavg, x1)
         end do
         bufTripole4DR8(nxGlobal,POP_haloWidth+1,k,n) = isign* &
         bufTripole4DR8(nxGlobal,POP_haloWidth+1,k,n)
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do n=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole4DR8(i   ,POP_haloWidth+1,k,n)
            x2 = bufTripole4DR8(iDst,POP_haloWidth+1,k,n)
            xavg = 0.5_POP_r8*(abs(x1) + abs(x2))
            bufTripole4DR8(i   ,POP_haloWidth+1,k,n) = isign*sign(xavg, x2)
            bufTripole4DR8(iDst,POP_haloWidth+1,k,n) = isign*sign(xavg, x1)
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
               do n=1,nt
               do k=1,nz
                  array(iDst,jDst,k,n,dstBlock) = isign* &
                                  bufTripole4DR8(iSrc,jSrc,k,n)
               end do
               end do
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!
!  clean up buffers
!
!-----------------------------------------------------------------------

   if (allocated(bufTripole4DR8)) then
      deallocate(bufTripole4DR8, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR8: error deallocating tripole buffer')
         return 
      endif
   endif

!-----------------------------------------------------------------------
!EOC

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
      i,k,n,nmsg,                &! dummy loop indices
      istat,                     &! status flag for allocate
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of 3rd, 4th dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   real (POP_r4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   real (POP_r4), dimension(:,:,:,:), allocatable :: &
      bufTripole4DR4    ! tripole buffer

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

   nz = size(array,dim=3)
   nt = size(array,dim=4)

   nxGlobal = 0
   if (allocated(bufTripoleR4)) then
      nxGlobal = size(bufTripoleR4,dim=1)
      allocate(bufTripole4DR4(nxGlobal,POP_haloWidth+1,nz,nt), stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR4: error allocating tripole buffer')
         return 
      endif
      bufTripole4DR4 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies to fill halo region
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
            do n=1,nt
            do k=1,nz
               array(iDst,jDst,k,n,dstBlock) = &
               array(iSrc,jSrc,k,n,srcBlock)
            end do
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do n=1,nt
            do k=1,nz
               bufTripole4DR4(iDst,jDst,k,n) = &
               array(iSrc,jSrc,k,n,srcBlock)
            end do
            end do
         endif
      else if (srcBlock == 0) then
         do n=1,nt
         do k=1,nz
            array(iDst,jDst,k,n,dstBlock) = fill
         end do
         end do
      endif
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do n=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole4DR4(i   ,POP_haloWidth+1,k,n)
            x2 = bufTripole4DR4(iDst,POP_haloWidth+1,k,n)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole4DR4(i   ,POP_haloWidth+1,k,n) = isign*sign(xavg, x2)
            bufTripole4DR4(iDst,POP_haloWidth+1,k,n) = isign*sign(xavg, x1)
         end do
         bufTripole4DR4(nxGlobal,POP_haloWidth+1,k,n) = isign* &
         bufTripole4DR4(nxGlobal,POP_haloWidth+1,k,n)
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do n=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole4DR4(i   ,POP_haloWidth+1,k,n)
            x2 = bufTripole4DR4(iDst,POP_haloWidth+1,k,n)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole4DR4(i   ,POP_haloWidth+1,k,n) = isign*sign(xavg, x2)
            bufTripole4DR4(iDst,POP_haloWidth+1,k,n) = isign*sign(xavg, x1)
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
               do n=1,nt
               do k=1,nz
                  array(iDst,jDst,k,n,dstBlock) = isign* &
                                  bufTripole4DR4(iSrc,jSrc,k,n)
               end do
               end do
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!
!  clean up buffers
!
!-----------------------------------------------------------------------

   if (allocated(bufTripole4DR4)) then
      deallocate(bufTripole4DR4, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DR4: error deallocating tripole buffer')
         return 
      endif
   endif

!-----------------------------------------------------------------------
!EOC

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
      i,k,n,nmsg,                &! dummy loop indices
      istat,                     &! status flag for allocate
      nxGlobal,                  &! global domain size in x (tripole)
      nz, nt,                    &! size of 3rd, 4th dimension
      iSrc,jSrc,                 &! source addresses for message
      iDst,jDst,                 &! dest   addresses for message
      srcBlock,                  &! local block number for source
      dstBlock,                  &! local block number for destination
      ioffset, joffset,          &! address shifts for tripole
      isign                       ! sign factor for tripole grids

   integer (POP_i4) :: &
      fill,            &! value to use for unknown points
      x1,x2,xavg        ! scalars for enforcing symmetry at U pts

   integer (POP_i4), dimension(:,:,:,:), allocatable :: &
      bufTripole4DI4    ! tripole buffer

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

   nz = size(array,dim=3)
   nt = size(array,dim=4)

   nxGlobal = 0
   if (allocated(bufTripoleI4)) then
      nxGlobal = size(bufTripoleI4,dim=1)
      allocate(bufTripole4DI4(nxGlobal,POP_haloWidth+1,nz,nt), stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DI4: error allocating tripole buffer')
         return 
      endif
      bufTripole4DI4 = fill
   endif

!-----------------------------------------------------------------------
!
!  do local copies to fill halo region
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
            do n=1,nt
            do k=1,nz
               array(iDst,jDst,k,n,dstBlock) = &
               array(iSrc,jSrc,k,n,srcBlock)
            end do
            end do
         else if (dstBlock < 0) then ! tripole copy into buffer
            do n=1,nt
            do k=1,nz
               bufTripole4DI4(iDst,jDst,k,n) = &
               array(iSrc,jSrc,k,n,srcBlock)
            end do
            end do
         endif
      else if (srcBlock == 0) then
         do n=1,nt
         do k=1,nz
            array(iDst,jDst,k,n,dstBlock) = fill
         end do
         end do
      endif
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do n=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal - i
            x1 = bufTripole4DI4(i   ,POP_haloWidth+1,k,n)
            x2 = bufTripole4DI4(iDst,POP_haloWidth+1,k,n)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole4DI4(i   ,POP_haloWidth+1,k,n) = isign*sign(xavg, x2)
            bufTripole4DI4(iDst,POP_haloWidth+1,k,n) = isign*sign(xavg, x1)
         end do
         bufTripole4DI4(nxGlobal,POP_haloWidth+1,k,n) = isign* &
         bufTripole4DI4(nxGlobal,POP_haloWidth+1,k,n)
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
         !*** swap locations with symmetric points so buffer has
         !***   correct values during the copy out

         do n=1,nt
         do k=1,nz
         do i = 1,nxGlobal/2
            iDst = nxGlobal + 1 - i
            x1 = bufTripole4DI4(i   ,POP_haloWidth+1,k,n)
            x2 = bufTripole4DI4(iDst,POP_haloWidth+1,k,n)
            xavg = 0.5_POP_r4*(abs(x1) + abs(x2))
            bufTripole4DI4(i   ,POP_haloWidth+1,k,n) = isign*sign(xavg, x2)
            bufTripole4DI4(iDst,POP_haloWidth+1,k,n) = isign*sign(xavg, x1)
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
               do n=1,nt
               do k=1,nz
                  array(iDst,jDst,k,n,dstBlock) = isign* &
                                  bufTripole4DI4(iSrc,jSrc,k,n)
               end do
               end do
            endif

         endif
      end do

   endif

!-----------------------------------------------------------------------
!
!  clean up buffers
!
!-----------------------------------------------------------------------

   if (allocated(bufTripole4DI4)) then
      deallocate(bufTripole4DI4, stat=istat)
      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloUpdate4DI4: error deallocating tripole buffer')
         return 
      endif
   endif

!-----------------------------------------------------------------------
!EOC

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

   subroutine POP_HaloMsgCreate(halo, distrb, srcBlock, dstBlock, &
                                direction, errorCode)

! !DESCRIPTION:
!  This is a utility routine to determine the required address and
!  message information for a particular pair of blocks.

! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      distrb                 ! block distribution

   integer (POP_i4), intent(in) :: &
      srcBlock,            & ! source      block id
      dstBlock               ! destination block id

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
      srcProc, srcLocalID,   &! source block location in distribution
      dstProc, dstLocalID,   &! source block location in distribution
      msgIndx,               &! message counter and index into msg array
      ibSrc, ieSrc, jbSrc, jeSrc, &! phys domain info for source block
      ibDst, ieDst, jbDst, jeDst, &! phys domain info for dest   block
      nxGlobal,              &! size of global domain in e-w direction
      i,j                     ! dummy loop index

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
!  find source and destination block locations
!
!-----------------------------------------------------------------------

   if (srcBlock /= 0) then
      call POP_DistributionGetBlockLoc(distrb, abs(srcBlock), srcProc, &
                                       srcLocalID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: error finding source block location')
         return
      endif
   else
      srcProc    = 0
      srcLocalID = 0
   endif

   if (dstBlock /= 0) then
      call POP_DistributionGetBlockLoc(distrb, abs(dstBlock), dstProc, &
                                       dstLocalID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: error finding dest block location')
         return
      endif
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

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('north')

         !*** copy northern physical domain of src
         !*** into southern halo of dst

         if (srcBlock > 0 .and. dstBlock > 0) then  ! normal north boundary

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

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

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         else if (srcBlock < 0 .and. dstBlock > 0) then

            !*** tripole grid - set up for copying out of 
            !*** tripole buffer into ghost cell domains
            !*** include e-w ghost cells

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('northeast')

         !*** normal northeast boundary - just copy NE corner
         !*** of physical domain into SW halo of NE nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         else

            !*** tripole grid - copy entire top halo+1 
            !*** rows into global buffer at src location

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         endif

      case ('northwest')

         !*** normal northeast boundary - just copy NW corner
         !*** of physical domain into SE halo of NW nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         else

            !*** tripole grid - copy entire top halo+1 
            !*** rows into global buffer at src location

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case default

         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: unknown direction local copy')
         return

      end select

!-----------------------------------------------------------------------
!
!  if dest block is local and source block does not exist, create a 
!  local copy to fill halo with a fill value
!
!-----------------------------------------------------------------------

   else if (srcProc == 0 .and. dstProc == POP_myTask+1) then   

      !*** compute addresses based on direction

      select case (direction)
      case ('east')

         !*** copy easternmost physical domain of src
         !*** into westernmost halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('west')

         !*** copy westernmost physical domain of src
         !*** into easternmost halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('north')

         !*** copy northern physical domain of src
         !*** into southern halo of dst

         if (dstBlock > 0) then  ! normal north boundary

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         endif

      case ('south')

         !*** copy southern physical domain of src
         !*** into northern halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('northeast')

         !*** normal northeast boundary - just copy NE corner
         !*** of physical domain into SW halo of NE nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         endif

      case ('northwest')

         !*** normal northeast boundary - just copy NW corner
         !*** of physical domain into SE halo of NW nbr block

         if (dstBlock > 0) then

            msgIndx = halo%numLocalCopies

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

            halo%numLocalCopies = msgIndx

         endif

      case ('southeast')

         !*** copy southeastern corner of src physical domain
         !*** into northwestern halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case ('southwest')

         !*** copy southwestern corner of src physical domain
         !*** into northeastern halo of dst

         msgIndx = halo%numLocalCopies

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

         halo%numLocalCopies = msgIndx

      case default

         call POP_ErrorSet(errorCode, &
            'POP_HaloMsgCreate: unknown direction local copy')
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

   end subroutine POP_HaloMsgCreate

!***********************************************************************
!BOP
! !IROUTINE: POP_HaloPrintStats
! !INTERFACE:

   subroutine POP_HaloPrintStats(halo, distrb, errorCode)

! !DESCRIPTION:
!  This routine compiles some message statistics for a given halo
!  updates and writes them to stdout.  For serial execution, this
!  routine does nothing.
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
!  do nothing
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!EOC

   end subroutine  POP_HaloPrintStats

!***********************************************************************

end module POP_HaloMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
