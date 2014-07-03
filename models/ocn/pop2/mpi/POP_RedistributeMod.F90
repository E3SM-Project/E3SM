!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: POP_RedistributeMod

 module POP_RedistributeMod

! !DESCRIPTION:
!  This module contains several routines for changing the distribution
!  of data on processors.  Two routines are supplied for gathering
!  or scattering distributed data into/out of a local global-sized
!  array.  Another routine is supplied for moving data from one block
!  distribution to another (e.g. from between baroclinic and barotropic
!  distributions).
!
! !REVISION HISTORY:
!  SVN: $Id $
!  2006-10-19: Phil Jones
!              new module for redistributing data that includes
!                the old gather/scatter routines for gathering
!                global arrays and the redistribution routine for
!                redistributing blocks.  more general routines may
!                follow later

! !USES:

   use POP_KindsMod
   use POP_CommMod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_ErrorMod

   implicit none
   private
   save

   include 'mpif.h'

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_RedistributeGather,      &
             POP_RedistributeScatter,     &
             POP_RedistributeBlocks

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  overload module functions
!
!-----------------------------------------------------------------------

   interface POP_RedistributeGather
     module procedure POP_RedistributeGatherR8, &
                      POP_RedistributeGatherR4, &
                      POP_RedistributeGatherI4
   end interface

   interface POP_RedistributeScatter
     module procedure POP_RedistributeScatterR8, &
                      POP_RedistributeScatterR4, &
                      POP_RedistributeScatterI4
   end interface

   interface POP_RedistributeBlocks
     module procedure POP_RedistributeBlocksR8, &
                      POP_RedistributeBlocksR4, &
                      POP_RedistributeBlocksI4
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeGatherR8
! !INTERFACE:

 subroutine POP_RedistributeGatherR8(arrayGlobal, arrayDistrb, &
                                     dstTask, distribution, errorCode, &
                                     fillValue)

! !DESCRIPTION:
!  This subroutine gathers a distributed 2d array to a global-sized
!  array on the processor dstTask.  An optional fillValue can be supplied
!  to fill data for land blocks that have been eliminated.  Otherwise,
!  zero is used to fill data in land blocks.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface POP\_RedistributeGather.  

! !USES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      dstTask        ! task/node to which array should be gathered

   type (POP_distrb), intent(in) :: &
      distribution   ! distribution of blocks for the source array

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      arrayDistrb    ! distributed array of data to be gathered

   real (POP_r8), intent(in), optional :: &
      fillValue      ! optional value for filling land blocks

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:), intent(out) :: &
      arrayGlobal    ! array containing global horizontal field on dstTask

   integer (POP_i4), intent(out) :: &
      errorCode      ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,n,           &! dummy loop counters
      srcTask,         &! task or processor location of a block
      blockIndex,      &! local index of block in distribution
      nsends,          &! number of messages sent
      nrecvs,          &! number of messages to recv
      ierr              ! MPI error flag

   real (POP_r8) ::    &
      fill              ! fill value to use for missing blocks

   type (POP_block) :: &
      thisBlock         ! block structure for current block


   integer (POP_i4), dimension(MPI_STATUS_SIZE) :: &
      mpiStatus         ! MPI status array for async sends

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest        ! MPI request array for async sends

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvInfo,         &! list of src tasks, block ids for recv
      sndStatus         ! MPI status array for async sends

   real (POP_r8), dimension(:,:), allocatable :: &
      msgBuffer         ! receive buffer

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

!-----------------------------------------------------------------------
!
!  if this task is the dstTask, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (POP_myTask == dstTask) then

      !***
      !*** do local copies to give time for messages to arrive
      !*** and save info for the expected receives
      !*** 

      allocate (rcvInfo(2,POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR8: error allocating rcvInfo')
         return
      endif

      nrecvs = 0

      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       srcTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR8: error getting block location')
            return
         endif

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(n, errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR8: error getting block info')
            return
         endif

         !*** copy local blocks

         if (srcTask == POP_myTask + 1) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayGlobal(thisBlock%iGlobal(i),   &
                           thisBlock%jGlobal(j)) = &
               arrayDistrb(i,j,blockIndex)
            end do
            end do

         !*** fill land blocks with fill, Phil

         else if (srcTask == 0) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayGlobal(thisBlock%iGlobal(i), &
                           thisBlock%jGlobal(j)) = fill
            end do
            end do

         !*** otherwise must recv a message - save info so we
         !*** can do the receives later

         else

            nrecvs = nrecvs + 1
            rcvInfo(1,nrecvs) = srcTask
            rcvInfo(2,nrecvs) = n        ! block id

         endif

      end do

      !*** 
      !*** now receive blocks to fill up the rest
      !*** 

      allocate (msgBuffer(POP_nxBlock,POP_nyBlock), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR8: error allocating buffer')
         return
      endif

      do n=1,nrecvs

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(rcvInfo(2,n), errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR8: error getting block info')
            return
         endif

         !*** if this is remote, receive a message

         call MPI_RECV(msgBuffer, size(msgBuffer),        &
                       POP_mpiR8, rcvInfo(1,n)-1,         & 
                       3*POP_mpitagRedist+rcvInfo(2,n),   &
                       POP_communicator, mpiStatus, ierr)

         if (ierr /= 0) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR8: error receiving msg')
            return
         endif

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i), &
                        thisBlock%jGlobal(j)) = msgBuffer(i,j)
         end do
         end do
      end do

      deallocate (msgBuffer, rcvInfo, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR8: error deallocating buffer')
         return
      endif

!-----------------------------------------------------------------------
!
!  otherwise send data to dstTask
!
!-----------------------------------------------------------------------

   else

      allocate(sndRequest(POP_numBlocks), &
               sndStatus (MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR8: error allocating MPI send status')
         return
      endif

      nsends = 0
      do n=1,POP_numBlocks

         !*** find block location

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       srcTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR8: error getting block location')
            return
         endif

         !*** if location is remote, must send data

         if (srcTask == POP_myTask+1) then

            nsends = nsends + 1

            call MPI_ISEND(arrayDistrb(1,1,blockIndex),        &
                     POP_nxBlock*POP_nyBlock,                  &
                     POP_mpiR8, dstTask, 3*POP_mpitagRedist+n, &
                     POP_communicator, sndRequest(nsends), ierr)
         endif
      end do

      if (nsends > 0) &
         call MPI_WAITALL(nsends, sndRequest, sndStatus, ierr)
      deallocate(sndRequest, sndStatus, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR8: error deallocating MPI status')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeGatherR8

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeGatherR4
! !INTERFACE:

 subroutine POP_RedistributeGatherR4(arrayGlobal, arrayDistrb, &
                                     dstTask, distribution, errorCode, &
                                     fillValue)

! !DESCRIPTION:
!  This subroutine gathers a distributed 2d array to a global-sized
!  array on the processor dstTask.  An optional fillValue can be supplied
!  to fill data for land blocks that have been eliminated.  Otherwise,
!  zero is used to fill data in land blocks.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for single precision arrays
!  corresponding to the generic interface POP\_RedistributeGather.  

! !USES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      dstTask        ! task/node to which array should be gathered

   type (POP_distrb), intent(in) :: &
      distribution   ! distribution of blocks for the source array

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      arrayDistrb    ! distributed array of data to be gathered

   real (POP_r4), intent(in), optional :: &
      fillValue      ! optional value for filling land blocks

! !OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:), intent(out) :: &
      arrayGlobal    ! array containing global horizontal field on dstTask

   integer (POP_i4), intent(out) :: &
      errorCode      ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,n,           &! dummy loop counters
      srcTask,         &! task or processor location of a block
      blockIndex,      &! local index of block in distribution
      nsends,          &! number of messages sent
      nrecvs,          &! number of messages to recv
      ierr              ! MPI error flag

   real (POP_r4) ::    &
      fill              ! fill value to use for missing blocks

   type (POP_block) :: &
      thisBlock         ! block structure for current block


   integer (POP_i4), dimension(MPI_STATUS_SIZE) :: &
      mpiStatus         ! MPI status array for async sends

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest        ! MPI request array for async sends

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvInfo,         &! list of src tasks, block ids for recv
      sndStatus         ! MPI status array for async sends

   real (POP_r4), dimension(:,:), allocatable :: &
      msgBuffer         ! receive buffer

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

!-----------------------------------------------------------------------
!
!  if this task is the dstTask, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (POP_myTask == dstTask) then

      !***
      !*** do local copies to give time for messages to arrive
      !*** and save info for the expected receives
      !*** 

      allocate (rcvInfo(2,POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR4: error allocating rcvInfo')
         return
      endif

      nrecvs = 0

      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       srcTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR4: error getting block location')
            return
         endif

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(n, errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR4: error getting block info')
            return
         endif

         !*** copy local blocks

         if (srcTask == POP_myTask + 1) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayGlobal(thisBlock%iGlobal(i),   &
                           thisBlock%jGlobal(j)) = &
               arrayDistrb(i,j,blockIndex)
            end do
            end do

         !*** fill land blocks with fill, Phil

         else if (srcTask == 0) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayGlobal(thisBlock%iGlobal(i), &
                           thisBlock%jGlobal(j)) = fill
            end do
            end do

         !*** otherwise must recv a message - save info so we
         !*** can do the receives later

         else

            nrecvs = nrecvs + 1
            rcvInfo(1,nrecvs) = srcTask
            rcvInfo(2,nrecvs) = n        ! block id

         endif

      end do

      !*** 
      !*** now receive blocks to fill up the rest
      !*** 

      allocate (msgBuffer(POP_nxBlock,POP_nyBlock), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR4: error allocating buffer')
         return
      endif

      do n=1,nrecvs

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(rcvInfo(2,n), errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR4: error getting block info')
            return
         endif

         !*** if this is remote, receive a message

         call MPI_RECV(msgBuffer, size(msgBuffer),        &
                       POP_mpiR4, rcvInfo(1,n)-1,         & 
                       3*POP_mpitagRedist+rcvInfo(2,n),   &
                       POP_communicator, mpiStatus, ierr)

         if (ierr /= 0) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR4: error receiving msg')
            return
         endif

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i), &
                        thisBlock%jGlobal(j)) = msgBuffer(i,j)
         end do
         end do
      end do

      deallocate (msgBuffer, rcvInfo, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR4: error deallocating buffer')
         return
      endif

!-----------------------------------------------------------------------
!
!  otherwise send data to dstTask
!
!-----------------------------------------------------------------------

   else

      allocate(sndRequest(POP_numBlocks), &
               sndStatus (MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR4: error allocating MPI send status')
         return
      endif

      nsends = 0
      do n=1,POP_numBlocks

         !*** find block location

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       srcTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherR4: error getting block location')
            return
         endif

         !*** if location is remote, must send data

         if (srcTask == POP_myTask+1) then

            nsends = nsends + 1

            call MPI_ISEND(arrayDistrb(1,1,blockIndex),        &
                     POP_nxBlock*POP_nyBlock,                  &
                     POP_mpiR4, dstTask, 3*POP_mpitagRedist+n, &
                     POP_communicator, sndRequest(nsends), ierr)
         endif
      end do

      if (nsends > 0) &
         call MPI_WAITALL(nsends, sndRequest, sndStatus, ierr)
      deallocate(sndRequest, sndStatus, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherR4: error deallocating MPI status')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeGatherR4

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeGatherI4
! !INTERFACE:

 subroutine POP_RedistributeGatherI4(arrayGlobal, arrayDistrb, &
                                     dstTask, distribution, errorCode, &
                                     fillValue)

! !DESCRIPTION:
!  This subroutine gathers a distributed 2d array to a global-sized
!  array on the processor dstTask.  An optional fillValue can be supplied
!  to fill data for land blocks that have been eliminated.  Otherwise,
!  zero is used to fill data in land blocks.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for integer arrays
!  corresponding to the generic interface POP\_RedistributeGather.  

! !USES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      dstTask        ! task/node to which array should be gathered

   type (POP_distrb), intent(in) :: &
      distribution   ! distribution of blocks for the source array

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      arrayDistrb    ! distributed array of data to be gathered

   integer (POP_i4), intent(in), optional :: &
      fillValue      ! optional value for filling land blocks

! !OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:), intent(out) :: &
      arrayGlobal    ! array containing global horizontal field on dstTask

   integer (POP_i4), intent(out) :: &
      errorCode      ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,n,           &! dummy loop counters
      srcTask,         &! task or processor location of a block
      blockIndex,      &! local index of block in distribution
      nsends,          &! number of messages sent
      nrecvs,          &! number of messages to recv
      ierr              ! MPI error flag

   integer (POP_i4) ::    &
      fill              ! fill value to use for missing blocks

   type (POP_block) :: &
      thisBlock         ! block structure for current block


   integer (POP_i4), dimension(MPI_STATUS_SIZE) :: &
      mpiStatus         ! MPI status array for async sends

   integer (POP_i4), dimension(:), allocatable :: &
      sndRequest        ! MPI request array for async sends

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvInfo,         &! list of src tasks, block ids for recv
      sndStatus         ! MPI status array for async sends

   integer (POP_i4), dimension(:,:), allocatable :: &
      msgBuffer         ! receive buffer

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

!-----------------------------------------------------------------------
!
!  if this task is the dstTask, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (POP_myTask == dstTask) then

      !***
      !*** do local copies to give time for messages to arrive
      !*** and save info for the expected receives
      !*** 

      allocate (rcvInfo(2,POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherI4: error allocating rcvInfo')
         return
      endif

      nrecvs = 0

      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       srcTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherI4: error getting block location')
            return
         endif

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(n, errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherI4: error getting block info')
            return
         endif

         !*** copy local blocks

         if (srcTask == POP_myTask + 1) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayGlobal(thisBlock%iGlobal(i),   &
                           thisBlock%jGlobal(j)) = &
               arrayDistrb(i,j,blockIndex)
            end do
            end do

         !*** fill land blocks with fill, Phil

         else if (srcTask == 0) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayGlobal(thisBlock%iGlobal(i), &
                           thisBlock%jGlobal(j)) = fill
            end do
            end do

         !*** otherwise must recv a message - save info so we
         !*** can do the receives later

         else

            nrecvs = nrecvs + 1
            rcvInfo(1,nrecvs) = srcTask
            rcvInfo(2,nrecvs) = n        ! block id

         endif

      end do

      !*** 
      !*** now receive blocks to fill up the rest
      !*** 

      allocate (msgBuffer(POP_nxBlock,POP_nyBlock), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherI4: error allocating buffer')
         return
      endif

      do n=1,nrecvs

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(rcvInfo(2,n), errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherI4: error getting block info')
            return
         endif

         !*** if this is remote, receive a message

         call MPI_RECV(msgBuffer, size(msgBuffer),        &
                       MPI_INTEGER, rcvInfo(1,n)-1,       & 
                       3*POP_mpitagRedist+rcvInfo(2,n),   &
                       POP_communicator, mpiStatus, ierr)

         if (ierr /= 0) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherI4: error receiving msg')
            return
         endif

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i), &
                        thisBlock%jGlobal(j)) = msgBuffer(i,j)
         end do
         end do
      end do

      deallocate (msgBuffer, rcvInfo, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherI4: error deallocating buffer')
         return
      endif

!-----------------------------------------------------------------------
!
!  otherwise send data to dstTask
!
!-----------------------------------------------------------------------

   else

      allocate(sndRequest(POP_numBlocks), &
               sndStatus (MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherI4: error allocating MPI send status')
         return
      endif

      nsends = 0
      do n=1,POP_numBlocks

         !*** find block location

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       srcTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeGatherI4: error getting block location')
            return
         endif

         !*** if location is remote, must send data

         if (srcTask == POP_myTask+1) then

            nsends = nsends + 1

            call MPI_ISEND(arrayDistrb(1,1,blockIndex),          &
                     POP_nxBlock*POP_nyBlock,                    &
                     MPI_INTEGER, dstTask, 3*POP_mpitagRedist+n, &
                     POP_communicator, sndRequest(nsends), ierr)
         endif
      end do

      if (nsends > 0) &
         call MPI_WAITALL(nsends, sndRequest, sndStatus, ierr)
      deallocate(sndRequest, sndStatus, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeGatherI4: error deallocating MPI status')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeGatherI4

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeScatterR8
! !INTERFACE:

 subroutine POP_RedistributeScatterR8(arrayDistrb, arrayGlobal, &
                                      srcTask, distribution, errorCode) 

! !DESCRIPTION:
!  This subroutine scatters data from a global-sized array on the 
!  processor srcTask to a distribution of blocks given by distribution.
!  {\bf NOTE: Only the physical domain of each block receives data.
!  If ghost cells/halo points need to be updated, a call to the
!  halo update routine is required.}
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface POP\_RedistributeScatter.

! !USES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask       ! task from which array should be scattered

   type (POP_distrb), intent(in) :: &
      distribution  ! distribution of blocks for distributed array

   real (POP_r8), dimension(:,:), intent(in) :: &
      arrayGlobal   ! array containing global field on src_task

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(out) :: &
      arrayDistrb   ! distributed array to hold result

   integer (POP_i4), intent(out) :: &
      errorCode    ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,n,           &! dummy loop indices
      nrecvs,          &! actual number of messages received
      dstTask,         &! proc/task location of block in distribution
      blockIndex,      &! local block index in distribution
      ierr              ! MPI and allocate error flag

   type (POP_block) :: &
      thisBlock         ! block info for current block

   integer (POP_i4), dimension(MPI_STATUS_SIZE) :: &
      mpiStatus         ! mpi array for async messages

   integer (POP_i4), dimension(:), allocatable :: &
      rcvRequest        ! request array for receives

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvStatus         ! status array for receives

   real (POP_r8), dimension(:,:), allocatable :: &
      msgBuffer         ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   arrayDistrb = 0.0_POP_r8

!-----------------------------------------------------------------------
!
!  if this task is the srcTask, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (POP_myTask == srcTask) then

      !*** send non-local blocks away

      allocate (msgBuffer(POP_nxBlock,POP_nyBlock), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR8: error allocating buffer')
         return
      endif

      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       dstTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterR8: error getting block location')
            return
         endif

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(n, errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterR8: error getting block info')
            return
         endif

         !*** if a non-local block, send the data to the
         !*** proper processor

         if (dstTask > 0 .and. dstTask - 1 /= POP_myTask) then

            !*** copy data into the send buffer

            msgBuffer = 0.0_POP_r8

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               msgBuffer(i,j) = arrayGlobal(thisBlock%iGlobal(i),&
                                            thisBlock%jGlobal(j))
            end do
            end do

            call MPI_SEND(msgBuffer, POP_nxBlock*POP_nyBlock,         &
                          POP_mpiR8, dstTask-1, 3*POP_mpitagRedist+n, &
                          POP_communicator, mpiStatus, ierr)

         !*** if a local block, copy the data directly

         else if (dstTask > 0 .and. dstTask - 1 == POP_myTask) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayDistrb(i,j,blockIndex) = &
               arrayGlobal(thisBlock%iGlobal(i),&
                           thisBlock%jGlobal(j))
            end do
            end do

         endif

      end do

      deallocate(msgBuffer, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR8: error deallocating buffer')
         return
      endif

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

      allocate (rcvRequest(POP_numBlocks), &
                rcvStatus(MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR8: error allocating MPI status')
         return
      endif

      rcvRequest = 0
      rcvStatus  = 0

      nrecvs = 0
      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                      dstTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterR8: error getting block location')
            return
         endif

         if (dstTask == POP_myTask+1) then
            nrecvs = nrecvs + 1
            call MPI_IRECV(arrayDistrb(1,1,blockIndex),          &
                       POP_nxBlock*POP_nyBlock,                  &
                       POP_mpiR8, srcTask, 3*POP_mpitagRedist+n, &
                       POP_communicator, rcvRequest(nrecvs), ierr)
         endif
      end do

      if (nrecvs > 0) &
        call MPI_WAITALL(nrecvs, rcvRequest, rcvStatus, ierr)

      deallocate(rcvRequest, rcvStatus, stat=ierr)
      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR8: error deallocating MPI status')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeScatterR8

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeScatterR4
! !INTERFACE:

 subroutine POP_RedistributeScatterR4(arrayDistrb, arrayGlobal, &
                                      srcTask, distribution, errorCode) 

! !DESCRIPTION:
!  This subroutine scatters data from a global-sized array on the 
!  processor srcTask to a distribution of blocks given by distribution.
!  {\bf NOTE: Only the physical domain of each block receives data.
!  If ghost cells/halo points need to be updated, a call to the
!  halo update routine is required.}
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for single precision arrays
!  corresponding to the generic interface POP\_RedistributeScatter.

! !USES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask       ! task from which array should be scattered

   type (POP_distrb), intent(in) :: &
      distribution  ! distribution of blocks for distributed array

   real (POP_r4), dimension(:,:), intent(in) :: &
      arrayGlobal   ! array containing global field on src_task

! !OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(out) :: &
      arrayDistrb   ! distributed array to hold result

   integer (POP_i4), intent(out) :: &
      errorCode    ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,n,           &! dummy loop indices
      nrecvs,          &! actual number of messages received
      dstTask,         &! proc/task location of block in distribution
      blockIndex,      &! local block index in distribution
      ierr              ! MPI and allocate error flag

   type (POP_block) :: &
      thisBlock         ! block info for current block

   integer (POP_i4), dimension(MPI_STATUS_SIZE) :: &
      mpiStatus         ! mpi array for async messages

   integer (POP_i4), dimension(:), allocatable :: &
      rcvRequest        ! request array for receives

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvStatus         ! status array for receives

   real (POP_r4), dimension(:,:), allocatable :: &
      msgBuffer         ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   arrayDistrb = 0.0_POP_r4

!-----------------------------------------------------------------------
!
!  if this task is the srcTask, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (POP_myTask == srcTask) then

      !*** send non-local blocks away

      allocate (msgBuffer(POP_nxBlock,POP_nyBlock), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR4: error allocating buffer')
         return
      endif

      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       dstTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterR4: error getting block location')
            return
         endif

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(n, errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterR4: error getting block info')
            return
         endif

         !*** if a non-local block, send the data to the
         !*** proper processor

         if (dstTask > 0 .and. dstTask - 1 /= POP_myTask) then

            !*** copy data into the send buffer

            msgBuffer = 0.0_POP_r4

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               msgBuffer(i,j) = arrayGlobal(thisBlock%iGlobal(i),&
                                            thisBlock%jGlobal(j))
            end do
            end do

            call MPI_SEND(msgBuffer, POP_nxBlock*POP_nyBlock,         &
                          POP_mpiR4, dstTask-1, 3*POP_mpitagRedist+n, &
                          POP_communicator, mpiStatus, ierr)

         !*** if a local block, copy the data directly

         else if (dstTask > 0 .and. dstTask - 1 == POP_myTask) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayDistrb(i,j,blockIndex) = &
               arrayGlobal(thisBlock%iGlobal(i),&
                           thisBlock%jGlobal(j))
            end do
            end do

         endif

      end do

      deallocate(msgBuffer, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR4: error deallocating buffer')
         return
      endif

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

      allocate (rcvRequest(POP_numBlocks), &
                rcvStatus(MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR4: error allocating MPI status')
         return
      endif

      rcvRequest = 0
      rcvStatus  = 0

      nrecvs = 0
      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                      dstTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterR4: error getting block location')
            return
         endif

         if (dstTask == POP_myTask+1) then
            nrecvs = nrecvs + 1
            call MPI_IRECV(arrayDistrb(1,1,blockIndex),          &
                       POP_nxBlock*POP_nyBlock,                  &
                       POP_mpiR4, srcTask, 3*POP_mpitagRedist+n, &
                       POP_communicator, rcvRequest(nrecvs), ierr)
         endif
      end do

      if (nrecvs > 0) &
        call MPI_WAITALL(nrecvs, rcvRequest, rcvStatus, ierr)

      deallocate(rcvRequest, rcvStatus, stat=ierr)
      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR4: error deallocating MPI status')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeScatterR4

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeScatterI4
! !INTERFACE:

 subroutine POP_RedistributeScatterI4(arrayDistrb, arrayGlobal, &
                                      srcTask, distribution, errorCode) 

! !DESCRIPTION:
!  This subroutine scatters data from a global-sized array on the 
!  processor srcTask to a distribution of blocks given by distribution.
!  {\bf NOTE: Only the physical domain of each block receives data.
!  If ghost cells/halo points need to be updated, a call to the
!  halo update routine is required.}
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for integer arrays
!  corresponding to the generic interface POP\_RedistributeScatter.

! !USES:

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      srcTask       ! task from which array should be scattered

   type (POP_distrb), intent(in) :: &
      distribution  ! distribution of blocks for distributed array

   integer (POP_i4), dimension(:,:), intent(in) :: &
      arrayGlobal   ! array containing global field on src_task

! !OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(out) :: &
      arrayDistrb   ! distributed array to hold result

   integer (POP_i4), intent(out) :: &
      errorCode    ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,n,           &! dummy loop indices
      nrecvs,          &! actual number of messages received
      dstTask,         &! proc/task location of block in distribution
      blockIndex,      &! local block index in distribution
      ierr              ! MPI and allocate error flag

   type (POP_block) :: &
      thisBlock         ! block info for current block

   integer (POP_i4), dimension(MPI_STATUS_SIZE) :: &
      mpiStatus         ! mpi array for async messages

   integer (POP_i4), dimension(:), allocatable :: &
      rcvRequest        ! request array for receives

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvStatus         ! status array for receives

   integer (POP_i4), dimension(:,:), allocatable :: &
      msgBuffer         ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   arrayDistrb = 0_POP_i4

!-----------------------------------------------------------------------
!
!  if this task is the srcTask, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (POP_myTask == srcTask) then

      !*** send non-local blocks away

      allocate (msgBuffer(POP_nxBlock,POP_nyBlock), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterI4: error allocating buffer')
         return
      endif

      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                       dstTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterI4: error getting block location')
            return
         endif

         !*** get block information for this block

         thisBlock = POP_BlocksGetBlock(n, errorCode)
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterI4: error getting block info')
            return
         endif

         !*** if a non-local block, send the data to the
         !*** proper processor

         if (dstTask > 0 .and. dstTask - 1 /= POP_myTask) then

            !*** copy data into the send buffer

            msgBuffer = 0_POP_i4

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               msgBuffer(i,j) = arrayGlobal(thisBlock%iGlobal(i),&
                                            thisBlock%jGlobal(j))
            end do
            end do

            call MPI_SEND(msgBuffer, POP_nxBlock*POP_nyBlock,           &
                          MPI_INTEGER, dstTask-1, 3*POP_mpitagRedist+n, &
                          POP_communicator, mpiStatus, ierr)

         !*** if a local block, copy the data directly

         else if (dstTask > 0 .and. dstTask - 1 == POP_myTask) then

            do j=thisBlock%jb,thisBlock%je
            do i=thisBlock%ib,thisBlock%ie
               arrayDistrb(i,j,blockIndex) = &
               arrayGlobal(thisBlock%iGlobal(i),&
                           thisBlock%jGlobal(j))
            end do
            end do

         endif

      end do

      deallocate(msgBuffer, stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterI4: error deallocating buffer')
         return
      endif

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

      allocate (rcvRequest(POP_numBlocks), &
                rcvStatus(MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterI4: error allocating MPI status')
         return
      endif

      rcvRequest = 0
      rcvStatus  = 0

      nrecvs = 0
      do n=1,POP_numBlocks

         !*** find location of this block in the distribution

         call POP_DistributionGetBlockLoc(distribution, n, &
                                      dstTask, blockIndex, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_RedistributeScatterI4: error getting block location')
            return
         endif

         if (dstTask == POP_myTask+1) then
            nrecvs = nrecvs + 1
            call MPI_IRECV(arrayDistrb(1,1,blockIndex),            &
                       POP_nxBlock*POP_nyBlock,                    &
                       MPI_INTEGER, srcTask, 3*POP_mpitagRedist+n, &
                       POP_communicator, rcvRequest(nrecvs), ierr)
         endif
      end do

      if (nrecvs > 0) &
        call MPI_WAITALL(nrecvs, rcvRequest, rcvStatus, ierr)

      deallocate(rcvRequest, rcvStatus, stat=ierr)
      if (ierr > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterI4: error deallocating MPI status')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeScatterI4

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeBlocksR8
! !INTERFACE:

 subroutine POP_RedistributeBlocksR8(dstArray, dstDistribution, &
                                     srcArray, srcDistribution, &
                                     errorCode)

! !DESCRIPTION:
!  This subroutine redistributes data from an array in which the
!  blocks are distributed in one decomposition to an array in which the
!  blocks are distributed differently. 
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface POP\_RedistributeBlocks.

! !USES:

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      srcDistribution, &! distribution of blocks for source      array
      dstDistribution   ! distribution of blocks for destination array

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      srcArray     ! array containing data in source distribution

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(out) :: &
      dstArray     ! array containing data in dest distribution

   integer (POP_i4), intent(out) :: &
      errorCode    ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! block loop index
      srcIndex,        &! local index for source      distribution
      dstIndex,        &! local index for destination distribution
      srcTask,         &! processor loc for block in source distribution
      dstTask,         &! processor loc for block in dest   distribution
      numSends,        &! number of messages sent from this task
      numRecvs,        &! number of messages received by this task
      ierr              ! MPI error flag

   integer (POP_i4), dimension(:), allocatable :: &
      rcvRequest,      &! request array for receives
      sndRequest        ! request array for sends

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvStatus,       &! status array for receives
      sndStatus         ! status array for sends

!BOC
!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcvRequest(POP_numBlocks), &
             sndRequest(POP_numBlocks), &
             rcvStatus(MPI_STATUS_SIZE, POP_numBlocks), &
             sndStatus(MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_RedistributeBlocksR8: error allocating status arrays')
      return
   endif

   rcvRequest = 0
   sndRequest = 0
   rcvStatus  = 0
   sndStatus  = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   numRecvs = 0
   numSends = 0

   do n=1,POP_numBlocks

      !*** find location of this block in each distribution

      call POP_DistributionGetBlockLoc(srcDistribution, n, &
                                       srcTask, srcIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeBlocksR8: error getting source location')
         return
      endif

      call POP_DistributionGetBlockLoc(dstDistribution, n, &
                                       dstTask, dstIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeBlocksR8: error getting dest location')
         return
      endif

      !*** if this destination is local and source is not, post a 
      !***  receive for this blocks

      if (dstTask == POP_myTask+1 .and. srcTask /= POP_myTask+1) then

         numRecvs = numRecvs + 1

         call MPI_IRECV(dstArray(1,1,dstIndex),                     &
                        POP_nxBlock*POP_nyBlock,                    &
                        POP_mpiR8, srcTask-1, 3*POP_mpitagRedist+n, &
                        POP_communicator, rcvRequest(numRecvs), ierr)
      endif

      !*** if this source is local and destination is not, post a 
      !***  send for this block

      if (srcTask == POP_myTask+1 .and. dstTask /= POP_myTask+1) then

         numSends = numSends + 1

         call MPI_ISEND(srcArray(1,1,srcIndex),                     &
                        POP_nxBlock*POP_nyBlock,                    &
                        POP_mpiR8, dstTask-1, 3*POP_mpitagRedist+n, &
                        POP_communicator, sndRequest(numSends), ierr)
      endif

      !*** if both blocks are local, simply copy the blocks

      if (srcTask == POP_myTask+1 .and. dstTask == POP_myTask+1) then

         dstArray(:,:,dstIndex) = srcArray(:,:,srcIndex)

      endif

   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (numSends /= 0) &
      call MPI_WAITALL(numSends, sndRequest, sndStatus, ierr)
   if (numRecvs /= 0) &
      call MPI_WAITALL(numRecvs, rcvRequest, rcvStatus, ierr)

   deallocate (rcvRequest, sndRequest, rcvStatus, sndStatus, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_RedistributeBlocksR8: error deallocating status arrays')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeBlocksR8

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeBlocksR4
! !INTERFACE:

 subroutine POP_RedistributeBlocksR4(dstArray, dstDistribution, &
                                     srcArray, srcDistribution, &
                                     errorCode)

! !DESCRIPTION:
!  This subroutine redistributes data from an array in which the
!  blocks are distributed in one decomposition to an array in which the
!  blocks are distributed differently. 
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for single precision arrays
!  corresponding to the generic interface POP\_RedistributeBlocks.

! !USES:

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      srcDistribution, &! distribution of blocks for source      array
      dstDistribution   ! distribution of blocks for destination array

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      srcArray     ! array containing data in source distribution

! !OUTPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(out) :: &
      dstArray     ! array containing data in dest distribution

   integer (POP_i4), intent(out) :: &
      errorCode    ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! block loop index
      srcIndex,        &! local index for source      distribution
      dstIndex,        &! local index for destination distribution
      srcTask,         &! processor loc for block in source distribution
      dstTask,         &! processor loc for block in dest   distribution
      numSends,        &! number of messages sent from this task
      numRecvs,        &! number of messages received by this task
      ierr              ! MPI error flag

   integer (POP_i4), dimension(:), allocatable :: &
      rcvRequest,      &! request array for receives
      sndRequest        ! request array for sends

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvStatus,       &! status array for receives
      sndStatus         ! status array for sends

!BOC
!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcvRequest(POP_numBlocks), &
             sndRequest(POP_numBlocks), &
             rcvStatus(MPI_STATUS_SIZE, POP_numBlocks), &
             sndStatus(MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_RedistributeBlocksR4: error allocating status arrays')
      return
   endif

   rcvRequest = 0
   sndRequest = 0
   rcvStatus  = 0
   sndStatus  = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   numRecvs = 0
   numSends = 0

   do n=1,POP_numBlocks

      !*** find location of this block in each distribution

      call POP_DistributionGetBlockLoc(srcDistribution, n, &
                                       srcTask, srcIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeBlocksR4: error getting source location')
         return
      endif

      call POP_DistributionGetBlockLoc(dstDistribution, n, &
                                       dstTask, dstIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeBlocksR4: error getting dest location')
         return
      endif

      !*** if this destination is local and source is not, post a 
      !***  receive for this blocks

      if (dstTask == POP_myTask+1 .and. srcTask /= POP_myTask+1) then

         numRecvs = numRecvs + 1

         call MPI_IRECV(dstArray(1,1,dstIndex),                     &
                        POP_nxBlock*POP_nyBlock,                    &
                        POP_mpiR4, srcTask-1, 3*POP_mpitagRedist+n, &
                        POP_communicator, rcvRequest(numRecvs), ierr)
      endif

      !*** if this source is local and destination is not, post a 
      !***  send for this block

      if (srcTask == POP_myTask+1 .and. dstTask /= POP_myTask+1) then

         numSends = numSends + 1

         call MPI_ISEND(srcArray(1,1,srcIndex),                     &
                        POP_nxBlock*POP_nyBlock,                    &
                        POP_mpiR4, dstTask-1, 3*POP_mpitagRedist+n, &
                        POP_communicator, sndRequest(numSends), ierr)
      endif

      !*** if both blocks are local, simply copy the blocks

      if (srcTask == POP_myTask+1 .and. dstTask == POP_myTask+1) then

         dstArray(:,:,dstIndex) = srcArray(:,:,srcIndex)

      endif

   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (numSends /= 0) &
      call MPI_WAITALL(numSends, sndRequest, sndStatus, ierr)
   if (numRecvs /= 0) &
      call MPI_WAITALL(numRecvs, rcvRequest, rcvStatus, ierr)

   deallocate (rcvRequest, sndRequest, rcvStatus, sndStatus, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_RedistributeBlocksR4: error deallocating status arrays')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeBlocksR4

!***********************************************************************
!BOP
! !IROUTINE: POP_RedistributeBlocksI4
! !INTERFACE:

 subroutine POP_RedistributeBlocksI4(dstArray, dstDistribution, &
                                     srcArray, srcDistribution, &
                                     errorCode)

! !DESCRIPTION:
!  This subroutine redistributes data from an array in which the
!  blocks are distributed in one decomposition to an array in which the
!  blocks are distributed differently. 
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for integer arrays
!  corresponding to the generic interface POP\_RedistributeBlocks.

! !USES:

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      srcDistribution, &! distribution of blocks for source      array
      dstDistribution   ! distribution of blocks for destination array

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      srcArray     ! array containing data in source distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(out) :: &
      dstArray     ! array containing data in dest distribution

   integer (POP_i4), intent(out) :: &
      errorCode    ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,               &! block loop index
      srcIndex,        &! local index for source      distribution
      dstIndex,        &! local index for destination distribution
      srcTask,         &! processor loc for block in source distribution
      dstTask,         &! processor loc for block in dest   distribution
      numSends,        &! number of messages sent from this task
      numRecvs,        &! number of messages received by this task
      ierr              ! MPI error flag

   integer (POP_i4), dimension(:), allocatable :: &
      rcvRequest,      &! request array for receives
      sndRequest        ! request array for sends

   integer (POP_i4), dimension(:,:), allocatable :: &
      rcvStatus,       &! status array for receives
      sndStatus         ! status array for sends

!BOC
!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcvRequest(POP_numBlocks), &
             sndRequest(POP_numBlocks), &
             rcvStatus(MPI_STATUS_SIZE, POP_numBlocks), &
             sndStatus(MPI_STATUS_SIZE, POP_numBlocks), stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_RedistributeBlocksI4: error allocating status arrays')
      return
   endif

   rcvRequest = 0
   sndRequest = 0
   rcvStatus  = 0
   sndStatus  = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   numRecvs = 0
   numSends = 0

   do n=1,POP_numBlocks

      !*** find location of this block in each distribution

      call POP_DistributionGetBlockLoc(srcDistribution, n, &
                                       srcTask, srcIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeBlocksI4: error getting source location')
         return
      endif

      call POP_DistributionGetBlockLoc(dstDistribution, n, &
                                       dstTask, dstIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeBlocksI4: error getting dest location')
         return
      endif

      !*** if this destination is local and source is not, post a 
      !***  receive for this blocks

      if (dstTask == POP_myTask+1 .and. srcTask /= POP_myTask+1) then

         numRecvs = numRecvs + 1

         call MPI_IRECV(dstArray(1,1,dstIndex),                       &
                        POP_nxBlock*POP_nyBlock,                      &
                        MPI_INTEGER, srcTask-1, 3*POP_mpitagRedist+n, &
                        POP_communicator, rcvRequest(numRecvs), ierr)
      endif

      !*** if this source is local and destination is not, post a 
      !***  send for this block

      if (srcTask == POP_myTask+1 .and. dstTask /= POP_myTask+1) then

         numSends = numSends + 1

         call MPI_ISEND(srcArray(1,1,srcIndex),                       &
                        POP_nxBlock*POP_nyBlock,                      &
                        MPI_INTEGER, dstTask-1, 3*POP_mpitagRedist+n, &
                        POP_communicator, sndRequest(numSends), ierr)
      endif

      !*** if both blocks are local, simply copy the blocks

      if (srcTask == POP_myTask+1 .and. dstTask == POP_myTask+1) then

         dstArray(:,:,dstIndex) = srcArray(:,:,srcIndex)

      endif

   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (numSends /= 0) &
      call MPI_WAITALL(numSends, sndRequest, sndStatus, ierr)
   if (numRecvs /= 0) &
      call MPI_WAITALL(numRecvs, rcvRequest, rcvStatus, ierr)

   deallocate (rcvRequest, sndRequest, rcvStatus, sndStatus, stat=ierr)

   if (ierr > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_RedistributeBlocksI4: error deallocating status arrays')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeBlocksI4

!***********************************************************************

 end module POP_RedistributeMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
