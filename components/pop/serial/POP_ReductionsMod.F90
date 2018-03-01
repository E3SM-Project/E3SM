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

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_GlobalSum,        &
             POP_GlobalSumProd,    &
             POP_GlobalCount,      &
             POP_GlobalMaxval,     &
             POP_GlobalMinval,     &
             POP_GlobalMaxloc,     &
             POP_GlobalMinloc

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface POP_GlobalSum
     module procedure POP_GlobalSum2DR8,        &
                      POP_GlobalSum2DR4,        &
                      POP_GlobalSum2DI4,        &
                      POP_GlobalSumScalarR8,    &
                      POP_GlobalSumScalarR4,    &
                      POP_GlobalSumScalarI4,    &
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

!EOC
!***********************************************************************

 contains

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
      localSum,      &! sum of local block domain
      globalSumTmp    ! quad version of global sum
#else
   real (POP_r8) :: &
      localSum       ! sum of local block domain
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   globalSumTmp = 0.0_POP_r16
#else
   globalSum = 0.0_POP_r8
#endif

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DR8: error getting distribution info')
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
      localSum = 0.0_POP_r16
#else
      localSum = 0.0_POP_r8
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array(i,j,iblock)
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
                     localSum = &
                     localSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
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
      localSum,     &! sum of local block domain
      globalSumTmp   ! hold higher precision global sum
#else
   real (POP_r4) :: &
      localSum       ! sum of local block domain
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   globalSumTmp = 0.0_POP_r8
#else
   globalSum = 0.0_POP_r4
#endif

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DR4: error getting distribution info')
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
      localSum = 0.0_POP_r8
#else
      localSum = 0.0_POP_r4
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array(i,j,iblock)
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
                     localSum = &
                     localSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
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
      localSum       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalSum = 0_POP_i4

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DI4: error getting distribution info')
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


      localSum = 0_POP_i4

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array(i,j,iblock)
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
                     localSum = &
                     localSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      globalSum = globalSum + localSum

   end do

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
!  Computes the global sum of the physical domain of a several
!  2-d arrays simultaneously.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic POP_GlobalSum
!  function corresponding to a stack of double precision arrays.  The
!  generic interface is identical but will handle real and integer 2-d
!  slabs and real, integer, and double precision scalars.

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:,:), intent(in) :: &
      array                ! set of arrays to be summed

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
      localSum,     &! sum of local block domain
      globalSumTmp   ! high precision form of global sum
#else
   real (POP_r8), dimension(size(array,dim=3)) :: &
      localSum       ! sum of local block domain
#endif

   integer (POP_i4) :: &
      i,j,n,iblock,    &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      nFields,         &! number of 2d fields to be summed
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   nFields   = size(array,dim=3)
#ifdef REPRODUCIBLE
   globalSumTmp = 0.0_POP_r16
#else
   globalSum = 0.0_POP_r8
#endif

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSum2DR8: error getting distribution info')
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
      localSum = 0.0_POP_r16
#else
      localSum = 0.0_POP_r8
#endif

      if (present(mMask)) then
         do n=1,nFields
         do j=jb,je
         do i=ib,ie
            localSum(n) = &
            localSum(n) + array(i,j,n,iblock)*mMask(i,j,iblock)
         end do
         end do
         end do
      else if (present(lMask)) then
         do n=1,nFields
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum(n) = &
               localSum(n) + array(i,j,n,iblock)
            endif
         end do
         end do
         end do
      else
         do n=1,nFields
         do j=jb,je
         do i=ib,ie
            localSum(n) = localSum(n) + array(i,j,n,iblock)
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
               do n=1,nFields
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum(n) = &
                     localSum(n) - array(i,j,n,iblock)*mMask(i,j,iblock)
                  endif
               end do
               end do
            else if (present(lMask)) then
               do n=1,nFields
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                     localSum(n) = localSum(n) - array(i,j,n,iblock)
                  endif
               end do
               end do
            else
               do n=1,nFields
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum(n) = localSum(n) - array(i,j,n,iblock)
                  endif
               end do
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
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
!  no operation needed for serial execution
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalSum = scalar

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
!  no operation needed for serial execution
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalSum = scalar

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
!  no operation needed for serial execution
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalSum = scalar

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
      localSum,      &! sum of local block domain
      globalSumTmp    ! higher precision global sum
#else
   real (POP_r8) :: &
      localSum       ! sum of local block domain
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   globalSumTmp = 0.0_POP_r16
#else
   globalSum = 0.0_POP_r8
#endif

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumProd2DR8: error getting distribution info')
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
      localSum = 0.0_POP_r16
#else
      localSum = 0.0_POP_r8
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array1(i,j,iblock)*array2(i,j,iblock)
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
                     localSum = &
                     localSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                        localSum = localSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum = localSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
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
      localSum,     &! sum of local block domain
      globalSumTmp   ! higher precision form of global sum
#else
   real (POP_r4) :: &
      localSum       ! sum of local block domain
#endif

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
#ifdef REPRODUCIBLE
   globalSumTmp = 0.0_POP_r8
#else
   globalSum = 0.0_POP_r4
#endif

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumProd2DR4: error getting distribution info')
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
      localSum = 0.0_POP_r8
#else
      localSum = 0.0_POP_r4
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array1(i,j,iblock)*array2(i,j,iblock)
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
                     localSum = &
                     localSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                        localSum = localSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum = localSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
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
      localSum       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalSum = 0_POP_i4

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalSumProd2DI4: error getting distribution info')
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


      localSum = 0_POP_i4

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array1(i,j,iblock)*array2(i,j,iblock)
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
                     localSum = &
                     localSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     if (lMask(i,j,iblock)) &
                        localSum = localSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (thisBlock%iGlobal(i) > thisBlock%nxGlobal/2) then
                     localSum = localSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      globalSum = globalSum + localSum

   end do

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
      localCount,        &! count of local block
      numBlocks,         &! number of local blocks in distribution
      blockID,           &! block id
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DR8: error getting distribution info')
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

      localCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock) /= 0.0_POP_r8) then
            localCount = localCount + 1
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
                     localCount = localCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      globalCount = globalCount + localCount

   end do

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
      localCount,        &! count of local block
      numBlocks,         &! number of local blocks in distribution
      blockID,           &! block id
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DR4: error getting distribution info')
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

      localCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock) /= 0.0_POP_r4) then
            localCount = localCount + 1
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
                     localCount = localCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      globalCount = globalCount + localCount

   end do

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
      localCount,        &! count of local block
      numBlocks,         &! number of local blocks in distribution
      blockID,           &! block id
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DI4: error getting distribution info')
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

      localCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock) /= 0_POP_i4) then
            localCount = localCount + 1
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
                     localCount = localCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      globalCount = globalCount + localCount

   end do

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
      localCount,        &! count of local block
      numBlocks,         &! number of local blocks in distribution
      blockID,           &! block id
      i,j,iblock          ! dummy counters

   type (POP_block) :: &
      thisBlock           ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalCount = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalCount2DLogical: error getting distribution info')
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

      localCount = 0

      do j=jb,je
      do i=ib,ie
         if (mask(i,j,iblock)) then
            localCount = localCount + 1
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
                     localCount = localCount - 1
               endif
            end do
         endif
      endif

      !*** now add block count to global count

      globalCount = globalCount + localCount

   end do

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

   real (POP_r8) :: &
      localMaxval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMaxval = -HUGE(0.0_POP_r8)

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxval2DR8: error getting distribution info')
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

      localMaxval = -HUGE(0.0_POP_r8)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMaxval = max(localMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMaxval = max(localMaxval,array(i,j,iblock))
         end do
         end do
      endif

      globalMaxval = max(globalMaxval,localMaxval)

   end do

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

   real (POP_r4) :: &
      localMaxval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMaxval = -HUGE(0.0_POP_r4)

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxval2DR4: error getting distribution info')
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

      localMaxval = -HUGE(0.0_POP_r4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMaxval = max(localMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMaxval = max(localMaxval,array(i,j,iblock))
         end do
         end do
      endif

      globalMaxval = max(globalMaxval,localMaxval)

   end do

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
      localMaxval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMaxval = -HUGE(0_POP_i4)

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxval2DI4: error getting distribution info')
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

      localMaxval = -HUGE(0_POP_i4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMaxval = max(localMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMaxval = max(localMaxval,array(i,j,iblock))
         end do
         end do
      endif

      globalMaxval = max(globalMaxval,localMaxval)

   end do

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

   real (POP_r8) :: &
      localMinval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMinval = HUGE(0.0_POP_r8)

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinval2DR8: error getting distribution info')
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

      localMinval = HUGE(0.0_POP_r8)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMinval = min(localMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMinval = min(localMinval,array(i,j,iblock))
         end do
         end do
      endif

      globalMinval = min(globalMinval,localMinval)

   end do

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

   real (POP_r4) :: &
      localMinval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMinval = HUGE(0.0_POP_r4)

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinval2DR4: error getting distribution info')
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

      localMinval = HUGE(0.0_POP_r4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMinval = min(localMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMinval = min(localMinval,array(i,j,iblock))
         end do
         end do
      endif

      globalMinval = min(globalMinval,localMinval)

   end do

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
      localMinval       ! sum of local block domain

   integer (POP_i4) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMinval = HUGE(0_POP_i4)

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinval2DI4: error getting distribution info')
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

      localMinval = HUGE(0_POP_i4)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMinval = min(localMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMinval = min(localMinval,array(i,j,iblock))
         end do
         end do
      endif

      globalMinval = min(globalMinval,localMinval)

   end do

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
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMaxval = scalar

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
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMaxval = scalar

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
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMaxval = scalar

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
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMinval = scalar

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
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMinval = scalar

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
!  no operations required for serial execution - return input value
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   globalMinval = scalar

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
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   maxValue = -HUGE(0.0_POP_r8)
   iLoc     = 0
   jLoc     = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxloc2DR8: error getting distribution info')
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
               if (array(i,j,iblock) > maxValue) then
                  maxValue = array(i,j,iblock)
                  iLoc = thisBlock%iGlobal(i)
                  jLoc = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) > maxValue) then
               maxValue = array(i,j,iblock)
               iLoc = thisBlock%iGlobal(i)
               jLoc = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

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
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   maxValue = -HUGE(0.0_POP_r4)
   iLoc     = 0
   jLoc     = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxloc2DR4: error getting distribution info')
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
               if (array(i,j,iblock) > maxValue) then
                  maxValue = array(i,j,iblock)
                  iLoc = thisBlock%iGlobal(i)
                  jLoc = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) > maxValue) then
               maxValue = array(i,j,iblock)
               iLoc = thisBlock%iGlobal(i)
               jLoc = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

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
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   maxValue = -HUGE(0_POP_i4)
   iLoc     = 0
   jLoc     = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMaxloc2DI4: error getting distribution info')
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
               if (array(i,j,iblock) > maxValue) then
                  maxValue = array(i,j,iblock)
                  iLoc = thisBlock%iGlobal(i)
                  jLoc = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) > maxValue) then
               maxValue = array(i,j,iblock)
               iLoc = thisBlock%iGlobal(i)
               jLoc = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

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
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   minValue = HUGE(0.0_POP_r8)
   iLoc     = 0
   jLoc     = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinloc2DR8: error getting distribution info')
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
               if (array(i,j,iblock) < minValue) then
                  minValue = array(i,j,iblock)
                  iLoc = thisBlock%iGlobal(i)
                  jLoc = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) < minValue) then
               minValue = array(i,j,iblock)
               iLoc = thisBlock%iGlobal(i)
               jLoc = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

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
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   minValue = HUGE(0.0_POP_r4)
   iLoc     = 0
   jLoc     = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinloc2DR4: error getting distribution info')
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
               if (array(i,j,iblock) < minValue) then
                  minValue = array(i,j,iblock)
                  iLoc = thisBlock%iGlobal(i)
                  jLoc = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) < minValue) then
               minValue = array(i,j,iblock)
               iLoc = thisBlock%iGlobal(i)
               jLoc = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

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
      numBlocks,       &! number of local blocks in distribution
      blockID           ! block location

   type (POP_Block) :: &
      thisBlock         ! block information for local block

!-----------------------------------------------------------------------

   errorCode = POP_Success
   minValue = HUGE(0_POP_i4)
   iLoc     = 0
   jLoc     = 0

   call POP_DistributionGet(dist, errorCode,          &
                            numLocalBlocks = numBlocks)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_GlobalMinloc2DI4: error getting distribution info')
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
               if (array(i,j,iblock) < minValue) then
                  minValue = array(i,j,iblock)
                  iLoc = thisBlock%iGlobal(i)
                  jLoc = thisBlock%jGlobal(j)
               endif
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            if (array(i,j,iblock) < minValue) then
               minValue = array(i,j,iblock)
               iLoc = thisBlock%iGlobal(i)
               jLoc = thisBlock%jGlobal(j)
            endif
         end do
         end do
      endif

   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_GlobalMinloc2DI4

!***********************************************************************

 end module POP_ReductionsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
