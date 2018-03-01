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
      blockIndex        ! local index of block in distribution

   real (POP_r8) :: &
      fill              ! fill value to use for missing blocks

   type (POP_block) :: &
      thisBlock         ! block structure for current block

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
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

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

      !*** since this is serial implementation, ignore the
      !*** processor id and just do local copies

      if (blockIndex /= 0) then  ! active block

         
         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i),   &
                        thisBlock%jGlobal(j)) = &
            arrayDistrb(i,j,blockIndex)
         end do
         end do

      else !*** eliminated land block - fill with fill, Phil

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i), &
                        thisBlock%jGlobal(j)) = fill
         end do
         end do

      endif

   end do

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
      blockIndex        ! local index of block in distribution

   real (POP_r4) :: &
      fill              ! fill value to use for missing blocks

   type (POP_block) :: &
      thisBlock         ! block structure for current block

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
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

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

      !*** since this is serial implementation, ignore the
      !*** processor id and just do local copies

      if (blockIndex /= 0) then  ! active block

         
         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i),   &
                        thisBlock%jGlobal(j)) = &
            arrayDistrb(i,j,blockIndex)
         end do
         end do

      else !*** eliminated land block - fill with fill, Phil

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i), &
                        thisBlock%jGlobal(j)) = fill
         end do
         end do

      endif

   end do

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
      blockIndex        ! local index of block in distribution

   integer (POP_i4) :: &
      fill              ! fill value to use for missing blocks

   type (POP_block) :: &
      thisBlock         ! block structure for current block

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
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

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

      !*** since this is serial implementation, ignore the
      !*** processor id and just do local copies

      if (blockIndex /= 0) then  ! active block

         
         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i),   &
                        thisBlock%jGlobal(j)) = &
            arrayDistrb(i,j,blockIndex)
         end do
         end do

      else !*** eliminated land block - fill with fill, Phil

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayGlobal(thisBlock%iGlobal(i), &
                        thisBlock%jGlobal(j)) = fill
         end do
         end do

      endif

   end do

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
      dstTask,         &! proc/task location of block in distribution
      blockIndex        ! local block index in distribution

   type (POP_block) :: &
      thisBlock         ! block info for current block

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   arrayDistrb = 0.0_POP_r8

!-----------------------------------------------------------------------
!
!  distribute the data based on distribution information
!
!-----------------------------------------------------------------------

   do n=1,POP_numBlocks

      !*** find location of this block in the distribution

      call POP_DistributionGetBlockLoc(distribution, n, &
                                       dstTask, blockIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR8: error getting block location')
         return
      endif

      if (blockIndex > size(arrayDistrb,dim=3)) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR8: distributed array not big enough')
         return
      endif

      !*** get block information for this block

      thisBlock = POP_BlocksGetBlock(n, errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR8: error getting block info')
         return
      endif

      !*** since this is serial implementation, ignore the
      !*** processor location and just do local copies

      if (blockIndex /= 0) then  !*** active block

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayDistrb(i,j,blockIndex) = &
            arrayGlobal(thisBlock%iGlobal(i),thisBlock%jGlobal(j))
         end do
         end do

      endif
   end do

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
      dstTask,         &! proc/task location of block in distribution
      blockIndex        ! local block index in distribution

   type (POP_block) :: &
      thisBlock         ! block info for current block

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   arrayDistrb = 0.0_POP_r4

!-----------------------------------------------------------------------
!
!  distribute the data based on distribution information
!
!-----------------------------------------------------------------------

   do n=1,POP_numBlocks

      !*** find location of this block in the distribution

      call POP_DistributionGetBlockLoc(distribution, n, &
                                       dstTask, blockIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR4: error getting block location')
         return
      endif

      if (blockIndex > size(arrayDistrb,dim=3)) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR4: distributed array not big enough')
         return
      endif

      !*** get block information for this block

      thisBlock = POP_BlocksGetBlock(n, errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterR4: error getting block info')
         return
      endif

      !*** since this is serial implementation, ignore the
      !*** processor location and just do local copies

      if (blockIndex /= 0) then  !*** active block

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayDistrb(i,j,blockIndex) = &
            arrayGlobal(thisBlock%iGlobal(i),thisBlock%jGlobal(j))
         end do
         end do

      endif
   end do

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
      dstTask,         &! proc/task location of block in distribution
      blockIndex        ! local block index in distribution

   type (POP_block) :: &
      thisBlock         ! block info for current block

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   arrayDistrb = 0_POP_i4

!-----------------------------------------------------------------------
!
!  distribute the data based on distribution information
!
!-----------------------------------------------------------------------

   do n=1,POP_numBlocks

      !*** find location of this block in the distribution

      call POP_DistributionGetBlockLoc(distribution, n, &
                                       dstTask, blockIndex, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterI4: error getting block location')
         return
      endif

      if (blockIndex > size(arrayDistrb,dim=3)) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterI4: distributed array not big enough')
         return
      endif

      !*** get block information for this block

      thisBlock = POP_BlocksGetBlock(n, errorCode)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_RedistributeScatterI4: error getting block info')
         return
      endif

      !*** since this is serial implementation, ignore the
      !*** processor location and just do local copies

      if (blockIndex /= 0) then  !*** active block

         do j=thisBlock%jb,thisBlock%je
         do i=thisBlock%ib,thisBlock%ie
            arrayDistrb(i,j,blockIndex) = &
            arrayGlobal(thisBlock%iGlobal(i),thisBlock%jGlobal(j))
         end do
         end do

      endif
   end do

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
      dstTask           ! processor loc for block in dest   distribution

!-----------------------------------------------------------------------
!
!  copy blocks from one distribution to another
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   dstArray  = 0.0_POP_r8

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

      if (srcIndex /= 0) then  !*** active block

         dstArray(:,:,dstIndex) = srcArray(:,:,srcIndex)

      endif
   end do

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
      dstTask           ! processor loc for block in dest   distribution

!-----------------------------------------------------------------------
!
!  copy blocks from one distribution to another
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   dstArray  = 0.0_POP_r4

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

      if (srcIndex /= 0) then  !*** active block

         dstArray(:,:,dstIndex) = srcArray(:,:,srcIndex)

      endif
   end do

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
      dstTask           ! processor loc for block in dest   distribution

!-----------------------------------------------------------------------
!
!  copy blocks from one distribution to another
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   dstArray  = 0_POP_i4

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

      if (srcIndex /= 0) then  !*** active block

         dstArray(:,:,dstIndex) = srcArray(:,:,srcIndex)

      endif
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_RedistributeBlocksI4

!***********************************************************************

 end module POP_RedistributeMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
