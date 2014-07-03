!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_ReductionsMod

!BOP
! !MODULE: POP_ReductionsMod
! !DESCRIPTION:
!  This module contains all the routines for performing global
!  reductions like global sums, minvals, maxvals, etc.
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:

   use POP_KindsMod
   use POP_CommMod
   use POP_ErrorMod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_GridHorzMod
   use registry

   implicit none
   private
   save

   include 'mpif.h'

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_GlobalSum,      &
             POP_GlobalSumProd,  &
             POP_GlobalCount,    &
             POP_GlobalMaxval,   &
             POP_GlobalMinval,   &
             POP_GlobalMaxloc,   &
             POP_GlobalMinloc,   &
             POP_initReductions

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface POP_GlobalSum
     module procedure POP_GlobalSum2DR8,     &
                      POP_GlobalSum2DR4,     &
                      POP_GlobalSum2DI4,     &
                      POP_GlobalSumScalarR8, &
                      POP_GlobalSumScalarR4, &
                      POP_GlobalSumScalarI4, &
                      POP_GlobalSumNfields2DR8
   end interface

   interface POP_GlobalSumProd
     module procedure POP_GlobalSumProd2DR8, &
                      POP_GlobalSumProd2DR4, &
                      POP_GlobalSumProd2DI4
   end interface

   interface POP_GlobalCount
     module procedure POP_GlobalCount2DR8,   &
                      POP_GlobalCount2DR4,   &
                      POP_GlobalCount2DI4,   &
                      POP_GlobalCount2DLogical
   end interface

   interface POP_GlobalMaxval
     module procedure POP_GlobalMaxval2DR8,     &
                      POP_GlobalMaxval2DR4,     &
                      POP_GlobalMaxval2DI4,     &
                      POP_GlobalMaxvalScalarR8, &
                      POP_GlobalMaxvalScalarR4, &
                      POP_GlobalMaxvalScalarI4
   end interface

   interface POP_GlobalMinval
     module procedure POP_GlobalMinval2DR8,     &
                      POP_GlobalMinval2DR4,     &
                      POP_GlobalMinval2DI4,     &
                      POP_GlobalMinvalScalarR8, &
                      POP_GlobalMinvalScalarR4, &
                      POP_GlobalMinvalScalarI4
   end interface

   interface POP_GlobalMaxloc
     module procedure POP_GlobalMaxloc2DR8,     &
                      POP_GlobalMaxloc2DR4,     &
                      POP_GlobalMaxloc2DI4
   end interface

   interface POP_GlobalMinloc
     module procedure POP_GlobalMinloc2DR8,     &
                      POP_GlobalMinloc2DR4,     &
                      POP_GlobalMinloc2DI4
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

    logical (POP_Logical) :: b4b

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_initReductions
! !INTERFACE:

 subroutine POP_initReductions

! !DESCRIPTION:
!  Initializes flags for global reductions.
!
! !REVISION HISTORY:
!  same as module
!EOP
!BOC
!-----------------------------------------------------------------------

   b4b = registry_match('b4b_flag')

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_initReductions

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSum
! !INTERFACE:

 function POP_GlobalSum2DR8(array, dist, fieldLoc, errorCode, &
                            mMask, lMask)                      &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to double precision arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

   real (POP_r8), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   real (POP_r16) :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision global sum
#else
   real (POP_r8) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains
   real (POP_r8), dimension(:), allocatable :: &
      blockSum_array_loc, &! sum of local blocks
      blockSum_array_glo   ! sum of all blocks
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numProcs,        &! number of processor participating
      numBlocks,       &! number of local blocks
      communicator      ! communicator for this distribution

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   localSum  = 0.0_POP_r16
#else
   localSum  = 0.0_POP_r8
   if (b4b) then
      allocate(blockSum_array_loc(POP_numBlocks), &
               blockSum_array_glo(POP_numBlocks))
      blockSum_array_loc = 0.0_POP_r8
   endif
#endif
   globalSum = 0.0_POP_r8

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

#ifdef REPRODUCIBLE
      blockSum = 0.0_POP_r16
#else
      blockSum = 0.0_POP_r8
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = &
                     blockSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

#ifndef REPRODUCIBLE
      if (b4b) then
         blockSum_array_loc(thisBlock%blockID) = blockSum
      endif
#endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSumTmp, 1, &
                         POP_mpiR16, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (.not. b4b) then
      if (POP_myTask < numProcs) then
         call MPI_ALLREDUCE(localSum, globalSum, 1, &
                            POP_mpiR8, MPI_SUM, communicator, ierr)
      endif
   else
      if (POP_myTask < numProcs) then
         call MPI_ALLREDUCE(blockSum_array_loc, blockSum_array_glo, POP_numBlocks, &
                            POP_mpiR8, MPI_SUM, communicator, ierr)
         globalSum = 0.0_POP_r8
         do iblock=1,POP_numBlocks
            globalSum = globalSum + blockSum_array_glo(iblock)
         end do
      endif
      deallocate(blockSum_array_loc, blockSum_array_glo)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSum2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSum
! !INTERFACE:

 function POP_GlobalSum2DR4(array, dist, fieldLoc, errorCode, &
                            mMask, lMask)                      &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to single precision arrays.

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

   real (POP_r4), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r4) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   real (POP_r8) :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision global sum
#else
   real (POP_r4) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   localSum  = 0.0_POP_r8
#else
   localSum  = 0.0_POP_r4
#endif
   globalSum = 0.0_POP_r4

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DR4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DR4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DR4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je


#ifdef REPRODUCIBLE
      blockSum = 0.0_POP_r8
#else
      blockSum = 0.0_POP_r4
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = &
                     blockSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSumTmp, 1, &
                         POP_mpiR8, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         POP_mpiR4, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSum2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSum
! !INTERFACE:

 function POP_GlobalSum2DI4(array, dist, fieldLoc, errorCode, &
                            mMask, lMask)                      &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to integer arrays.

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

   integer (POP_i4), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   integer (POP_i4) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localSum  = 0_POP_i4
   globalSum = 0_POP_i4

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DI4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DI4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DI4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je


      blockSum = 0_POP_i4

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = &
                     blockSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = blockSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSum2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSum
! !INTERFACE:

 function POP_GlobalSumNfields2DR8(array, dist, fieldLoc, errorCode, &
                                   mMask, lMask)                     &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a set of
!  2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to a stack of double precision arrays.  The
!  generic interface is identical but will handle real and integer 
!  2-d slabs and real, integer, and double precision scalars.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:), intent(in) :: &
      array                ! array to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

   real (POP_r8), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8), dimension(size(array,dim=3)) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   real (POP_r16), dimension(size(array,dim=3)) :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision global sum
#else
   real (POP_r8), dimension(size(array,dim=3)) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains
   real (POP_r8), dimension(:,:), allocatable :: &
      blockSum_array_loc, &! sum of local blocks
      blockSum_array_glo   ! sum of all blocks
#endif

   integer (POP_i4) :: &
      i,j,n,iblock,    &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numFields,       &! number of 2d arrays to sum
      numProcs,        &! number of processor participating
      numBlocks,       &! number of local blocks
      communicator      ! communicator for this distribution

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   localSum  = 0.0_POP_r16
#else
   localSum  = 0.0_POP_r8
   if (b4b) then
      allocate(blockSum_array_loc(size(array,dim=3),POP_numBlocks), &
               blockSum_array_glo(size(array,dim=3),POP_numBlocks))
      blockSum_array_loc = 0.0_POP_r8
   endif
#endif
   globalSum = 0.0_POP_r8

   numFields = size(array,dim=3)

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSum2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

#ifdef REPRODUCIBLE
      blockSum = 0.0_POP_r16
#else
      blockSum = 0.0_POP_r8
#endif

      if (present(mMask)) then
         do n=1,numFields
         do j=jb,je
         do i=ib,ie
            blockSum(n) = &
            blockSum(n) + array(i,j,n,iblock)*mMask(i,j,iblock)
         end do
         end do
         end do
      else if (present(lMask)) then
         do n=1,numFields
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum(n) = &
               blockSum(n) + array(i,j,n,iblock)
            endif
         end do
         end do
         end do
      else
         do n=1,numFields
         do j=jb,je
         do i=ib,ie
            blockSum(n) = blockSum(n) + array(i,j,n,iblock)
         end do
         end do
         end do
      endif

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            if (present(mMask)) then
               do n=1,numFields
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum(n) = &
                     blockSum(n) - array(i,j,n,iblock)*mMask(i,j,iblock)
                  endif
               end do
               end do
            else if (present(lMask)) then
               do n=1,numFields
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     blockSum(n) = blockSum(n) - array(i,j,n,iblock)
                  endif
               end do
               end do
            else
               do n=1,numFields
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum(n) = blockSum(n) - array(i,j,n,iblock)
                  endif
               end do
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum(:) = localSum(:) + blockSum(:)

#ifndef REPRODUCIBLE
      if (b4b) then
         blockSum_array_loc(:,thisBlock%blockID) = blockSum(:)
      endif
#endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSumTmp, numFields, &
                         POP_mpiR16, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp 
   endif
#else
   if (.not. b4b) then
      if (POP_myTask < numProcs) then
         call MPI_ALLREDUCE(localSum, globalSum, numFields, &
                            POP_mpiR8, MPI_SUM, communicator, ierr)
      endif
   else
      if (POP_myTask < numProcs) then
         call MPI_ALLREDUCE(blockSum_array_loc, blockSum_array_glo, numFields*POP_numBlocks, &
                            POP_mpiR8, MPI_SUM, communicator, ierr)
         globalSum = 0.0_POP_r8
         do iblock=1,POP_numBlocks
            globalSum(:) = globalSum(:) + blockSum_array_glo(:,iblock)
         end do
      endif
      deallocate(blockSum_array_loc, blockSum_array_glo)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSumNfields2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSum
! !INTERFACE:

 function POP_GlobalSumScalarR8(scalar, dist, errorCode) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to double precision scalars.

! !INPUT PARAMETERS:

   real (POP_r8), intent(in) :: &
      scalar               ! scalar to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution 

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

#ifdef REPRODUCIBLE
   real (POP_r16) :: &
      scalarTmp, globalSumTmp  ! higher precision for reproducibility
#endif
!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumScalarR8: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (POP_myTask < numProcs) then
      scalarTmp = scalar
      call MPI_ALLREDUCE(scalarTmp, globalSumTmp, 1, &
                         POP_mpiR16, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalSum, 1, &
                         POP_mpiR8, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSumScalarR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSum
! !INTERFACE:

 function POP_GlobalSumScalarR4(scalar, dist, errorCode) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to single precision scalars.

! !INPUT PARAMETERS:

   real (POP_r4), intent(in) :: &
      scalar               ! scalar to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution 

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r4) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

#ifdef REPRODUCIBLE
   real (POP_r8) :: &
      scalarTmp, globalSumTmp  ! higher precision for reproducibility
#endif
!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumScalarR4: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (POP_myTask < numProcs) then
      scalarTmp = scalar
      call MPI_ALLREDUCE(scalarTmp, globalSumTmp, 1, &
                         POP_mpiR8, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalSum, 1, &
                         POP_mpiR4, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSumScalarR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSum
! !INTERFACE:

 function POP_GlobalSumScalarI4(scalar, dist, errorCode) &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to integer scalars.

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      scalar               ! scalar to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution 

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   integer (POP_i4) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumScalarI4: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSumScalarI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSumProd
! !INTERFACE:

 function POP_GlobalSumProd2DR8(array1, array2, dist, fieldLoc, &
                                errorCode, mMask, lMask)        &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalSumProd function corresponding to double precision 
!  arrays.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

   real (POP_r8), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   real (POP_r16) :: &
      blockSum,      &! sum of local block domain
      localSum,      &! sum of all local block domains
      globalSumTmp    ! higher precision for reproducibility
#else
   real (POP_r8) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   localSum  = 0.0_POP_r16
#else
   localSum  = 0.0_POP_r8
#endif
   globalSum = 0.0_POP_r8

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumProd2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSumProd2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSumProd2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

#ifdef REPRODUCIBLE
      blockSum = 0.0_POP_r16
#else
      blockSum = 0.0_POP_r8
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = &
                     blockSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                        blockSum = blockSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = blockSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSumTmp, 1, &
                         POP_mpiR16, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         POP_mpiR8, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSumProd2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSumProd
! !INTERFACE:

 function POP_GlobalSumProd2DR4(array1, array2, dist, fieldLoc, &
                                errorCode, mMask, lMask)        &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalSumProd function corresponding to single precision 
!  arrays.

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

   real (POP_r4), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r4) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   real (POP_r8) :: &
      blockSum,     &! sum of local block domain
      localSum,     &! sum of all local block domains
      globalSumTmp   ! higher precision for reproducibility
#else
   real (POP_r4) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   localSum  = 0.0_POP_r8
#else
   localSum  = 0.0_POP_r4
#endif
   globalSum = 0.0_POP_r4

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumProd2DR4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSumProd2DR4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSumProd2DR4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

#ifdef REPRODUCIBLE
      blockSum = 0.0_POP_r8
#else
      blockSum = 0.0_POP_r4
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = &
                     blockSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                        blockSum = blockSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = blockSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSumTmp, 1, &
                         POP_mpiR8, MPI_SUM, communicator, ierr)
      globalSum = globalSumTmp
   endif
#else
   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         POP_mpiR4, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSumProd2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalSumProd
! !INTERFACE:

 function POP_GlobalSumProd2DI4(array1, array2, dist, fieldLoc, &
                                errorCode, mMask, lMask)        &
          result(globalSum)

! !DESCRIPTION:
!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalSumProd function corresponding to integer arrays.

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

   integer (POP_i4), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   integer (POP_i4) :: &
      globalSum            ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localSum  = 0_POP_i4
   globalSum = 0_POP_i4

   call POP_DistributionGet(dist, errorCode,     &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumProd2DI4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSumProd2DI4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalSumProd2DI4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je


      blockSum = 0_POP_i4

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            blockSum = &
            blockSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockSum = &
               blockSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = &
                     blockSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                        blockSum = blockSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     blockSum = blockSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalSumProd2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalCount
! !INTERFACE:

 function POP_GlobalCount2DR8 (mask, dist, fieldLoc, errorCode) &
          result (globalCount)

! !DESCRIPTION:
!  This function returns the number of true or non-zero elements
!  in the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic
!  POP_GlobalCount function corresponding to double precision arrays.

!INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      mask                 ! array for which non-zero elements
                           !   are to be counted

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for MASK

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

!OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code 

   integer (POP_i4) :: &
      globalCount          ! resulting global count

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::   &
      ib, ie, jb, je,    &! start,end of physical domain
      ierr,              &! mpi error flag
      blockCount,        &! count of local block
      localCount,        &! count of all local blocks
      blockID,           &! block location
      numBlocks,       &! number of local blocks
      numProcs,          &! number of processor participating
      communicator,      &! communicator for this distribution
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localCount  = 0
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock) /= 0.0_POP_r8) then
            blockCount = blockCount + 1
         endif
      end do
      end do

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            do i=ib,ie
               if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                  if (mask(i,j,iblock) /= 0.0_POP_r8) &
                     blockCount = blockCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      localCount = localCount + blockCount

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localCount, globalCount, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalCount2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalCount
! !INTERFACE:

 function POP_GlobalCount2DR4 (mask, dist, fieldLoc, errorCode) &
          result (globalCount)

! !DESCRIPTION:
!  This function returns the number of true or non-zero elements
!  in the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic
!  POP_GlobalCount function corresponding to single precision arrays.

!INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      mask                 ! array for which non-zero elements
                           !   are to be counted

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for MASK

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

!OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code 

   integer (POP_i4) :: &
      globalCount          ! resulting global count

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::   &
      ib, ie, jb, je,    &! start,end of physical domain
      ierr,              &! mpi error flag
      blockCount,        &! count of local block
      localCount,        &! count of all local blocks
      blockID,           &! block id
      numBlocks,       &! number of local blocks
      numProcs,          &! number of processor participating
      communicator,      &! communicator for this distribution
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localCount  = 0
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,     &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DR4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DR4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DR4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock) /= 0.0_POP_r4) then
            blockCount = blockCount + 1
         endif
      end do
      end do

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            do i=ib,ie
               if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                  if (mask(i,j,iblock) /= 0.0_POP_r4) &
                     blockCount = blockCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      localCount = localCount + blockCount

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localCount, globalCount, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalCount2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalCount
! !INTERFACE:

 function POP_GlobalCount2DI4 (mask, dist, fieldLoc, errorCode) &
          result (globalCount)

! !DESCRIPTION:
!  This function returns the number of true or non-zero elements
!  in the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic
!  POP_GlobalCount function corresponding to integer arrays.

!INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      mask                 ! array for which non-zero elements
                           !   are to be counted

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for MASK

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

!OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code 

   integer (POP_i4) :: &
      globalCount          ! resulting global count

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::   &
      ib, ie, jb, je,    &! start,end of physical domain
      ierr,              &! mpi error flag
      blockCount,        &! count of local block
      localCount,        &! count of all local blocks
      blockID,           &! block id
      numBlocks,       &! number of local blocks
      numProcs,          &! number of processor participating
      communicator,      &! communicator for this distribution
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localCount  = 0
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DI4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DI4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DI4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock) /= 0_POP_i4) then
            blockCount = blockCount + 1
         endif
      end do
      end do

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            do i=ib,ie
               if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                  if (mask(i,j,iblock) /= 0_POP_i4) &
                     blockCount = blockCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      localCount = localCount + blockCount

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localCount, globalCount, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalCount2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalCount
! !INTERFACE:

 function POP_GlobalCount2DLogical (mask, dist, fieldLoc, errorCode) &
          result (globalCount)

! !DESCRIPTION:
!  This function returns the number of true elements
!  in the physical domain of a 2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic
!  POP_GlobalCount function corresponding to logical arrays.

!INPUT PARAMETERS:

   logical (POP_logical), dimension(:,:,:), intent(in) :: &
      mask                 ! array for which non-zero elements
                           !   are to be counted

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for MASK

   character (*), intent(in) :: &
      fieldLoc             ! grid stagger location for this field

!OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code 

   integer (POP_i4) :: &
      globalCount          ! resulting global count

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::   &
      ib, ie, jb, je,    &! start,end of physical domain
      ierr,              &! mpi error flag
      blockCount,        &! count of local block
      localCount,        &! count of all local blocks
      blockID,           &! block id
      numBlocks,       &! number of local blocks
      numProcs,          &! number of processor participating
      communicator,      &! communicator for this distribution
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localCount  = 0
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DLogical: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DLogical: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalCount2DLogical: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock)) then
            blockCount = blockCount + 1
         endif
      end do
      end do

      !*** if this block along tripole boundary and field 
      !*** located on north face and northeast corner points
      !*** must eliminate redundant points from global sum

      if (thisBlock%tripole) then
         if (fieldLoc == POP_gridHorzLocNface .or. &
             fieldLoc == POP_gridHorzLocNEcorner) then

            j = je

            do i=ib,ie
               if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                  if (mask(i,j,iblock)) &
                     blockCount = blockCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      localCount = localCount + blockCount

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localCount, globalCount, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalCount2DLogical

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxval
! !INTERFACE:

 function POP_GlobalMaxval2DR8(array, dist, errorCode, lMask) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of the physical domain of a 
!  2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxval function corresponding to double precision arrays.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8) :: &
      globalMaxval         ! resulting maximum valu of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (POP_r8) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMaxval  = -HUGE(0.0_POP_r8)
   globalMaxval = -HUGE(0.0_POP_r8)

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxval2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxval2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxval2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockMaxval = -HUGE(0.0_POP_r8)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         POP_mpiR8, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMaxval2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxval
! !INTERFACE:

 function POP_GlobalMaxval2DR4(array, dist, errorCode, lMask) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of the physical domain of a 
!  2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxval function corresponding to single precision arrays.

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r4) :: &
      globalMaxval         ! resulting maximum valu of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (POP_r4) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMaxval  = -HUGE(0.0_POP_r4)
   globalMaxval = -HUGE(0.0_POP_r4)

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxval2DR4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxval2DR4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxval2DR4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockMaxval = -HUGE(0.0_POP_r4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         POP_mpiR4, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMaxval2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxval
! !INTERFACE:

 function POP_GlobalMaxval2DI4(array, dist, errorCode, lMask) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of the physical domain of a 
!  2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxval function corresponding to integer arrays.

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   integer (POP_i4) :: &
      globalMaxval         ! resulting maximum valu of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMaxval  = -HUGE(0_POP_i4)
   globalMaxval = -HUGE(0_POP_i4)

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxval2DI4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxval2DI4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxval2DI4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockMaxval = -HUGE(0_POP_i4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMaxval2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinval
! !INTERFACE:

 function POP_GlobalMinval2DR8(array, dist, errorCode, lMask) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of the physical domain of a 
!  2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinval function corresponding to double precision arrays.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8) :: &
      globalMinval         ! resulting minimum valu of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (POP_r8) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMinval  = HUGE(0.0_POP_r8)
   globalMinval = HUGE(0.0_POP_r8)

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinval2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinval2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinval2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockMinval = HUGE(0.0_POP_r8)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         POP_mpiR8, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMinval2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinval
! !INTERFACE:

 function POP_GlobalMinval2DR4(array, dist, errorCode, lMask) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of the physical domain of a 
!  2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinval function corresponding to single precision arrays.

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r4) :: &
      globalMinval         ! resulting minimum valu of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (POP_r4) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMinval  = HUGE(0.0_POP_r4)
   globalMinval = HUGE(0.0_POP_r4)

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinval2DR4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinval2DR4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinval2DR4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockMinval = HUGE(0.0_POP_r4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         POP_mpiR4, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMinval2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinval
! !INTERFACE:

 function POP_GlobalMinval2DI4(array, dist, errorCode, lMask) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of the physical domain of a 
!  2-d array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinval function corresponding to integer arrays.

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   integer (POP_i4) :: &
      globalMinval         ! resulting minimum valu of array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of all local block domains

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMinval  = HUGE(0_POP_i4)
   globalMinval = HUGE(0_POP_i4)

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinval2DI4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinval2DI4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinval2DI4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      blockMinval = HUGE(0_POP_i4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         MPI_INTEGER, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMinval2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxval
! !INTERFACE:

 function POP_GlobalMaxvalScalarR8 (scalar, dist, errorCode) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxval function corresponding to double precision scalars.

! !INPUT PARAMETERS:

   real (POP_r8), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! current block distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8) :: &
      globalMaxval         ! resulting maximum value of scalar

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxvalScalarR8: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         POP_mpiR8, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMaxvalScalarR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxval
! !INTERFACE:

 function POP_GlobalMaxvalScalarR4 (scalar, dist, errorCode) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxval function corresponding to single precision scalars.

! !INPUT PARAMETERS:

   real (POP_r4), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! current block distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r4) :: &
      globalMaxval         ! resulting maximum value of scalar

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxvalScalarR4: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         POP_mpiR4, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMaxvalScalarR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxval
! !INTERFACE:

 function POP_GlobalMaxvalScalarI4 (scalar, dist, errorCode) &
          result(globalMaxval)

! !DESCRIPTION:
!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxval function corresponding to integer scalars.

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! current block distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   integer (POP_i4) :: &
      globalMaxval         ! resulting maximum value of scalar

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxvalScalarI4: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMaxvalScalarI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinval
! !INTERFACE:

 function POP_GlobalMinvalScalarR8 (scalar, dist, errorCode) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinval function corresponding to double precision scalars.

! !INPUT PARAMETERS:

   real (POP_r8), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! current block distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r8) :: &
      globalMinval         ! resulting minimum value of scalar

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinvalScalarR8: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         POP_mpiR8, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMinvalScalarR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinval
! !INTERFACE:

 function POP_GlobalMinvalScalarR4 (scalar, dist, errorCode) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinval function corresponding to single precision scalars.

! !INPUT PARAMETERS:

   real (POP_r4), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! current block distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   real (POP_r4) :: &
      globalMinval         ! resulting minimum value of scalar

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinvalScalarR4: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         POP_mpiR4, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMinvalScalarR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinval
! !INTERFACE:

 function POP_GlobalMinvalScalarI4 (scalar, dist, errorCode) &
          result(globalMinval)

! !DESCRIPTION:
!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinval function corresponding to integer scalars.

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (POP_distrb), intent(in) :: &
      dist                 ! current block distribution

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error flag

   integer (POP_i4) :: &
      globalMinval         ! resulting minimum value of scalar

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

!-----------------------------------------------------------------------
!
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(dist, errorCode,     &
                            numProcs = numProcs, &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinvalScalarI4: error getting communicator')
      return
   endif

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         MPI_INTEGER, MPI_MIN, communicator, ierr)
   endif

!-----------------------------------------------------------------------
!EOC

 end function POP_GlobalMinvalScalarI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxloc
! !INTERFACE:

 subroutine POP_GlobalMaxloc2DR8(array, dist, &
                                 iLoc, jLoc, maxValue, errorCode, lMask)

! !DESCRIPTION:
!  This routine finds the location of the global maximum for the
!  physical domain of a 2-d field and returns the global domain
!  index and maximum value of that location.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxloc function corresponding to double precision arrays.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      array                ! array for which maxloc required

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (POP_r8), intent(out) :: &
       maxValue            ! maximum value of the field

   integer (POP_i4), intent(out) :: &
       errorCode,         &! returned errorCode
       iLoc, jLoc          ! global i,j location of maximum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

   real (POP_r8) ::    &
      localMaxval       ! max value from all local blocks

   integer (POP_i4), dimension(2) :: &
      localMaxAddr,    &! global address of local maxval
      globalMaxAddr     ! global address of global maxval

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMaxval = -HUGE(0.0_POP_r8)
   maxValue = -HUGE(0.0_POP_r8)
   iLoc     = 0
   jLoc     = 0
   localMaxAddr(:) = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxloc2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxloc2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxloc2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               if (array(i,j,iblock) > localMaxval) then
                  localMaxval = array(i,j,iblock)
                  localMaxAddr(1) = thisBlock%iGlobal(i)
                  localMaxAddr(2) = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) > localMaxval) then
               localMaxval = array(i,j,iblock)
               localMaxAddr(1) = thisBlock%iGlobal(i)
               localMaxAddr(2) = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to find global maxval from local maxval
!  then find address of that value
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMaxval, maxValue, 1, &
                         POP_mpiR8, MPI_MAX, communicator, ierr)

      if (localMaxval /= maxValue) then
         localMaxAddr(:) = 0
      endif

      call MPI_ALLREDUCE(localMaxAddr, globalMaxAddr, 2, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)

      iLoc = globalMaxAddr(1)
      jLoc = globalMaxAddr(2)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_GlobalMaxloc2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxloc
! !INTERFACE:

 subroutine POP_GlobalMaxloc2DR4(array, dist, &
                                 iLoc, jLoc, maxValue, errorCode, lMask)

! !DESCRIPTION:
!  This routine finds the location of the global maximum for the
!  physical domain of a 2-d field and returns the global domain
!  index and maximum value of that location.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxloc function corresponding to single precision arrays.

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      array                ! array for which maxloc required

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (POP_r4), intent(out) :: &
       maxValue            ! maximum value of the field

   integer (POP_i4), intent(out) :: &
       errorCode,         &! returned errorCode
       iLoc, jLoc          ! global i,j location of maximum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

   real (POP_r4) ::    &
      localMaxval       ! max value from all local blocks

   integer (POP_i4), dimension(2) :: &
      localMaxAddr,    &! global address of local maxval
      globalMaxAddr     ! global address of global maxval

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMaxval = -HUGE(0.0_POP_r4)
   maxValue    = -HUGE(0.0_POP_r4)
   iLoc     = 0
   jLoc     = 0
   localMaxAddr(:) = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxloc2DR4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxloc2DR4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxloc2DR4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               if (array(i,j,iblock) > localMaxval) then
                  localMaxval = array(i,j,iblock)
                  localMaxAddr(1) = thisBlock%iGlobal(i)
                  localMaxAddr(2) = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) > localMaxval) then
               localMaxval = array(i,j,iblock)
               localMaxAddr(1) = thisBlock%iGlobal(i)
               localMaxAddr(2) = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to find global maxval from local maxval
!  then find address of that value
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMaxval, maxValue, 1, &
                         POP_mpiR4, MPI_MAX, communicator, ierr)

      if (localMaxval /= maxValue) then
         localMaxAddr(:) = 0
      endif

      call MPI_ALLREDUCE(localMaxAddr, globalMaxAddr, 2, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)

      iLoc = globalMaxAddr(1)
      jLoc = globalMaxAddr(2)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_GlobalMaxloc2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMaxloc
! !INTERFACE:

 subroutine POP_GlobalMaxloc2DI4(array, dist, &
                                 iLoc, jLoc, maxValue, errorCode, lMask)

! !DESCRIPTION:
!  This routine finds the location of the global maximum for the
!  physical domain of a 2-d field and returns the global domain
!  index and maximum value of that location.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMaxloc function corresponding to integer arrays.

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      array                ! array for which maxloc required

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
       maxValue            ! maximum value of the field

   integer (POP_i4), intent(out) :: &
       errorCode,         &! returned errorCode
       iLoc, jLoc          ! global i,j location of maximum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

   integer (POP_i4) ::    &
      localMaxval       ! max value from all local blocks

   integer (POP_i4), dimension(2) :: &
      localMaxAddr,    &! global address of local maxval
      globalMaxAddr     ! global address of global maxval

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMaxval = -HUGE(0_POP_i4)
   maxValue    = -HUGE(0_POP_i4)
   iLoc     = 0
   jLoc     = 0
   localMaxAddr(:) = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxloc2DI4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxloc2DI4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMaxloc2DI4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               if (array(i,j,iblock) > localMaxval) then
                  localMaxval = array(i,j,iblock)
                  localMaxAddr(1) = thisBlock%iGlobal(i)
                  localMaxAddr(2) = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) > localMaxval) then
               localMaxval = array(i,j,iblock)
               localMaxAddr(1) = thisBlock%iGlobal(i)
               localMaxAddr(2) = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to find global maxval from local maxval
!  then find address of that value
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMaxval, maxValue, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)

      if (localMaxval /= maxValue) then
         localMaxAddr(:) = 0
      endif

      call MPI_ALLREDUCE(localMaxAddr, globalMaxAddr, 2, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)

      iLoc = globalMaxAddr(1)
      jLoc = globalMaxAddr(2)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_GlobalMaxloc2DI4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinloc
! !INTERFACE:

 subroutine POP_GlobalMinloc2DR8(array, dist, &
                                 iLoc, jLoc, minValue, errorCode, lMask)

! !DESCRIPTION:
!  This routine finds the location of the global minimum for the
!  physical domain of a 2-d field and returns the global domain
!  index and minimum value of that location.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinloc function corresponding to double precision arrays.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      array                ! array for which minloc required

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (POP_r8), intent(out) :: &
       minValue            ! minimum value of the field

   integer (POP_i4), intent(out) :: &
       errorCode,         &! returned errorCode
       iLoc, jLoc          ! global i,j location of minimum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

   real (POP_r8) ::    &
      localMinval       ! max value from all local blocks

   integer (POP_i4), dimension(2) :: &
      localMinAddr,    &! global address of local maxval
      globalMinAddr     ! global address of global maxval

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMinval = HUGE(0.0_POP_r8)
   minValue    = HUGE(0.0_POP_r8)
   iLoc     = 0
   jLoc     = 0
   localMinAddr(:) = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinloc2DR8: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinloc2DR8: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinloc2DR8: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               if (array(i,j,iblock) < localMinval) then
                  localMinval = array(i,j,iblock)
                  localMinAddr(1) = thisBlock%iGlobal(i)
                  localMinAddr(2) = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) < localMinval) then
               localMinval = array(i,j,iblock)
               localMinAddr(1) = thisBlock%iGlobal(i)
               localMinAddr(2) = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to find global minval from local minval
!  then find address of that value
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMinval, minValue, 1, &
                         POP_mpiR8, MPI_MIN, communicator, ierr)

      if (localMinval /= minValue) then
         localMinAddr(:) = 0
      endif

      call MPI_ALLREDUCE(localMinAddr, globalMinAddr, 2, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)

      iLoc = globalMinAddr(1)
      jLoc = globalMinAddr(2)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_GlobalMinloc2DR8

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinloc
! !INTERFACE:

 subroutine POP_GlobalMinloc2DR4(array, dist, &
                                 iLoc, jLoc, minValue, errorCode, lMask)

! !DESCRIPTION:
!  This routine finds the location of the global minimum for the
!  physical domain of a 2-d field and returns the global domain
!  index and minimum value of that location.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinloc function corresponding to single precision arrays.

! !INPUT PARAMETERS:

   real (POP_r4), dimension(:,:,:), intent(in) :: &
      array                ! array for which minloc required

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   real (POP_r4), intent(out) :: &
       minValue            ! minimum value of the field

   integer (POP_i4), intent(out) :: &
       errorCode,         &! returned errorCode
       iLoc, jLoc          ! global i,j location of minimum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

   real (POP_r4) ::    &
      localMinval       ! max value from all local blocks

   integer (POP_i4), dimension(2) :: &
      localMinAddr,    &! global address of local maxval
      globalMinAddr     ! global address of global maxval

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMinval = HUGE(0.0_POP_r4)
   minValue    = HUGE(0.0_POP_r4)
   iLoc     = 0
   jLoc     = 0
   localMinAddr(:) = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinloc2DR4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinloc2DR4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinloc2DR4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               if (array(i,j,iblock) < localMinval) then
                  localMinval = array(i,j,iblock)
                  localMinAddr(1) = thisBlock%iGlobal(i)
                  localMinAddr(2) = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) < localMinval) then
               localMinval = array(i,j,iblock)
               localMinAddr(1) = thisBlock%iGlobal(i)
               localMinAddr(2) = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to find global maxval from local maxval
!  then find address of that value
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMinval, minValue, 1, &
                         POP_mpiR4, MPI_MIN, communicator, ierr)

      if (localMinval /= minValue) then
         localMinAddr(:) = 0
      endif

      call MPI_ALLREDUCE(localMinAddr, globalMinAddr, 2, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)

      iLoc = globalMinAddr(1)
      jLoc = globalMinAddr(2)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_GlobalMinloc2DR4

!***********************************************************************
!BOP
! !IROUTINE: POP_GlobalMinloc
! !INTERFACE:

 subroutine POP_GlobalMinloc2DI4(array, dist, &
                                 iLoc, jLoc, minValue, errorCode, lMask)

! !DESCRIPTION:
!  This routine finds the location of the global minimum for the
!  physical domain of a 2-d field and returns the global domain
!  index and minimum value of that location.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  POP_GlobalMinloc function corresponding to integer arrays.

! !INPUT PARAMETERS:

   integer (POP_i4), dimension(:,:,:), intent(in) :: &
      array                ! array for which minloc required

   type (POP_distrb), intent(in) :: &
      dist                 ! block distribution for array

   logical (POP_logical), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
       minValue            ! minimum value of the field

   integer (POP_i4), intent(out) :: &
       errorCode,         &! returned errorCode
       iLoc, jLoc          ! global i,j location of minimum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

   integer (POP_i4) ::    &
      localMinval       ! max value from all local blocks

   integer (POP_i4), dimension(2) :: &
      localMinAddr,    &! global address of local maxval
      globalMinAddr     ! global address of global maxval

!-----------------------------------------------------------------------

   errorCode = POP_Success
   localMinval = HUGE(0_POP_i4)
   minValue    = HUGE(0_POP_i4)
   iLoc     = 0
   jLoc     = 0
   localMinAddr(:) = 0

   call POP_DistributionGet(dist, errorCode,            &
                            numLocalBlocks = numBlocks, &
                            numProcs = numProcs,        &
                            communicator = communicator)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinloc2DI4: error getting communicator')
      return
   endif

   do iblock=1,numBlocks
      call POP_DistributionGetBlockID(dist, iblock, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinloc2DI4: error getting block id')
         return
      endif

      thisBlock = POP_BlocksGetBlock(blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_GlobalMinloc2DI4: error getting block')
         return
      endif

      ib = thisBlock%ib
      ie = thisBlock%ie
      jb = thisBlock%jb
      je = thisBlock%je

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               if (array(i,j,iblock) < localMinval) then
                  localMinval = array(i,j,iblock)
                  localMinAddr(1) = thisBlock%iGlobal(i)
                  localMinAddr(2) = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) < localMinval) then
               localMinval = array(i,j,iblock)
               localMinAddr(1) = thisBlock%iGlobal(i)
               localMinAddr(2) = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to find global minval from local minval
!  then find address of that value
!
!-----------------------------------------------------------------------

   if (POP_myTask < numProcs) then
      call MPI_ALLREDUCE(localMinval, minValue, 1, &
                         MPI_INTEGER, MPI_MIN, communicator, ierr)

      if (localMinval /= minValue) then
         localMinAddr(:) = 0
      endif

      call MPI_ALLREDUCE(localMinAddr, globalMinAddr, 2, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)

      iLoc = globalMinAddr(1)
      jLoc = globalMinAddr(2)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_GlobalMinloc2DI4

!***********************************************************************

 end module POP_ReductionsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
