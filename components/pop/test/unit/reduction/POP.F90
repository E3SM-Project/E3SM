!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 program POP_ReductionTest

!----------------------------------------------------------------------
!
!  this program tests POP global reduction operations, including
!  global sums, minval, maxval, etc.  This test code relies on the
!  redistribution mod working correctly for gather/scatter operations.
!
!----------------------------------------------------------------------

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_DomainSizeMod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_RedistributeMod
   use POP_ReductionsMod
   use POP_GridHorzMod

   implicit none

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: &
      errorCode,       &! error flag
      istat,           &! allocate error flag
      nProcsClinic,    &! number of processors in clinic distribution
      blockID,         &! global block id
      ib,ie,jb,je,     &! beg, end of physical domain
      i,j,n             ! loop indices

   type (POP_distrb) :: &
      distrbClinic       ! baroclinic distribution

   integer (POP_i4), dimension(:), allocatable :: &
      workPerBlock  ! test values

   integer (POP_i4), dimension(2) :: &
      iLoc          ! location index for minloc, maxloc

   integer (POP_i4), dimension(:), pointer :: &
      iGlobal       ! global index for this block

   integer (POP_i4), dimension(:), pointer :: &
      jGlobal       ! global index for this block

   integer (POP_i4) :: & 
      i4Test, i4Expect

   real (POP_r4) :: &
      r4Test, r4Expect  ! test values

   real (POP_r8) :: &
      r8Test, r8Expect  ! test values

   integer (POP_i4), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksClinic) :: &
      i4Array2D, i4Mask

   real (POP_r4), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksClinic) :: &
      r4Array2D, r4Mask

   real (POP_r8), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksClinic) :: &
      r8Array2D, r8Mask

   logical (POP_logical), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksClinic) :: &
      mask

   integer (POP_i4), &
      dimension(POP_nxGlobal, POP_nyGlobal) :: &
      i4ArrayG

   real (POP_r4), &
      dimension(POP_nxGlobal, POP_nyGlobal) :: &
      r4ArrayG

   real (POP_r8), &
      dimension(POP_nxGlobal, POP_nyGlobal) :: &
      r8ArrayG

!----------------------------------------------------------------------
!
!  initialize communcation environment
!
!----------------------------------------------------------------------

   call POP_CommInitMessageEnvironment
   call POP_CommInit

   errorCode = POP_Success

!----------------------------------------------------------------------
!
!  create block decomposition and distributions
!
!----------------------------------------------------------------------

   !*** create block decomposition
 
   call POP_BlocksCreate(POP_nxGlobal, POP_nyGlobal, &
                         'cyclic', 'tripole', errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error creating block decomposition')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   !*** initialize global arrays with land in the middle

   do j=1,POP_nyGlobal
   do i=1,POP_nxGlobal
      if (i > POP_nxGlobal/2 + POP_nxBlock .or. &
          i < POP_nxGlobal/2 - POP_nxBlock .or. & 
          j > POP_nyGlobal/2 + POP_nyBlock .or. &
          j < POP_nyGlobal/2 - POP_nyBlock) then
         i4ArrayG(i,j) = (i+j)*10_POP_i4
         r4ArrayG(i,j) = (i+j)*100.0_POP_r4
         r8ArrayG(i,j) = (i+j)*1000.0_POP_r8
      else
         i4ArrayG(i,j) = 0_POP_i4
         r4ArrayG(i,j) = 0.0_POP_r4
         r8ArrayG(i,j) = 0.0_POP_r8
      endif
   end do
   end do

   !*** create artificial work per block with varying work and
   !*** at least one empty block 


   allocate(workPerBlock(POP_numBlocks), stat=istat)
   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error allocating work per block')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   do blockID=1,POP_numBlocks
      call POP_BlocksGetBlockInfo(blockID, errorCode,         &
                                  ib=ib, ie=ie, jb=jb, je=je, &
                                  iGlobal=iGlobal, jGlobal=jGlobal)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ReductionTest: error getting block info')
         call POP_ErrorPrint(errorCode)
         stop
      endif

      workPerBlock(blockID) = 0
      do j=jb,je
      do i=ib,ie
         if (i4ArrayG(iGlobal(i),jGlobal(j)) /= 0) then
            workPerBlock(blockID) = workPerBlock(blockID) + 1
         endif
      end do
      end do
   end do

   !*** create baroclinic distribution

   nprocsClinic = POP_CommGetNumProcs(POP_Communicator)

   distrbClinic = POP_DistributionCreate(POP_distributionMethodRake, &
                                         nProcsClinic, &
                                         workPerBlock, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error creating clinic distribution')
      call POP_ErrorPrint(errorCode)
      stop
   endif

!----------------------------------------------------------------------
!
!  initialize arrays by scattering global arrays to distribution
!
!----------------------------------------------------------------------

   call POP_RedistributeScatter(i4Array2D, i4ArrayG, &
                              POP_masterTask, distrbClinic, errorCode) 

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 scatter')
   endif

   call POP_RedistributeScatter(r4Array2D, r4ArrayG, &
                              POP_masterTask, distrbClinic, errorCode) 

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 scatter')
   endif

   call POP_RedistributeScatter(r8Array2D, r8ArrayG, &
                              POP_masterTask, distrbClinic, errorCode) 

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 scatter')
   endif

   mask = (i4Array2D /= 0)

!----------------------------------------------------------------------
!
!  first test scalar sum operations - use myTask as scalar
!
!----------------------------------------------------------------------

   i4Expect = 0_POP_i4
   do n=1,nprocsClinic
      i4Expect = i4Expect + n
   end do
   r4Expect = i4Expect*10.0_POP_r4
   r8Expect = i4Expect*100.0_POP_r8

   i4Test = POP_GlobalSum((POP_myTask + 1), distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 scalar sum')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 scalar sum')
   endif

   r4Test = POP_GlobalSum((POP_myTask + 1)*10.0_POP_r4, &
                                            distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 scalar sum')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 scalar sum')
   endif

   r8Test = POP_GlobalSum((POP_myTask + 1)*100.0_POP_r8, &
                                            distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 scalar sum')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 scalar sum')
   endif

!----------------------------------------------------------------------
!
!  now test scalar maxval operations - use myTask as scalar
!
!----------------------------------------------------------------------

   i4Expect = nprocsClinic
   r4Expect = i4Expect*10.0_POP_r4
   r8Expect = i4Expect*100.0_POP_r8

   i4Test = POP_GlobalMaxval((POP_myTask + 1), distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 scalar maxval')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 scalar maxval')
   endif

   r4Test = POP_GlobalMaxval((POP_myTask + 1)*10.0_POP_r4, &
                                            distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 scalar maxval')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 scalar maxval')
   endif

   r8Test = POP_GlobalMaxval((POP_myTask + 1)*100.0_POP_r8, &
                                            distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 scalar maxval')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 scalar maxval')
   endif

!----------------------------------------------------------------------
!
!  now test scalar minval operations - use myTask as scalar
!
!----------------------------------------------------------------------

   i4Expect = 1
   r4Expect = i4Expect*10.0_POP_r4
   r8Expect = i4Expect*100.0_POP_r8

   i4Test = POP_GlobalMinval((POP_myTask + 1), distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 scalar minval')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 scalar minval')
   endif

   r4Test = POP_GlobalMinval((POP_myTask + 1)*10.0_POP_r4, &
                                            distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 scalar minval')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 scalar minval')
   endif

   r8Test = POP_GlobalMinval((POP_myTask + 1)*100.0_POP_r8, &
                                            distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 scalar minval')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 scalar minval')
   endif

!----------------------------------------------------------------------
!
!  now test global maxval functions
!
!----------------------------------------------------------------------

   i4Expect = maxval(i4ArrayG)
   r4Expect = maxval(r4ArrayG)
   r8Expect = maxval(r8ArrayG)

   i4Test = POP_GlobalMaxval(i4Array2D, distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 array maxval')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 array maxval')
   endif

   r4Test = POP_GlobalMaxval(r4Array2D, distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 array maxval')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 array maxval')
   endif

   r8Test = POP_GlobalMaxval(r8Array2D, distrbClinic, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 array maxval')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 array maxval')
   endif

!----------------------------------------------------------------------
!
!  test maxloc
!
!----------------------------------------------------------------------

   i4Expect = maxval(i4ArrayG)
   r4Expect = maxval(r4ArrayG)
   r8Expect = maxval(r8ArrayG)
   iLoc = maxloc(i4ArrayG)

   call POP_GlobalMaxloc(i4Array2D, distrbClinic, &
                         i, j, i4Test, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 array maxloc')
   endif
   if (i4Test /= i4Expect .or. &
       iLoc(1) /= i .or. iLoc(2) /= j) then
      print *,POP_myTask,i,j,iLoc(1),iLoc(2)
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 array maxloc')
   endif

   call POP_GlobalMaxloc(r4Array2D, distrbClinic, &
                         i, j, r4Test, errorCode)
   if (errorCode /= POP_Success) then
      print *,POP_myTask,i,j,iLoc(1),iLoc(2)
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 array maxloc')
   endif
   if (r4Test /= r4Expect .or. &
       iLoc(1) /= i .or. iLoc(2) /= j) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 array maxloc')
   endif

   call POP_GlobalMaxloc(r8Array2D, distrbClinic, &
                         i, j, r8Test, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 array maxloc')
   endif
   if (r8Test /= r8Expect .or. &
       iLoc(1) /= i .or. iLoc(2) /= j) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 array maxloc')
   endif

!----------------------------------------------------------------------
!
!  now test global minval functions
!
!----------------------------------------------------------------------

   i4Expect = minval(i4ArrayG, mask=(i4ArrayG /= 0_POP_i4))
   r4Expect = minval(r4ArrayG, mask=(r4ArrayG /= 0.0_POP_r4))
   r8Expect = minval(r8ArrayG, mask=(r8ArrayG /= 0.0_POP_r8))

   i4Test = POP_GlobalMinval(i4Array2D, distrbClinic, errorCode, &
                             lMask = mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 array minval')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 array minval')
   endif

   r4Test = POP_GlobalMinval(r4Array2D, distrbClinic, errorCode, &
                             lMask = mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 array minval')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 array minval')
   endif

   r8Test = POP_GlobalMinval(r8Array2D, distrbClinic, errorCode, &
                             lMask = mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 array minval')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 array minval')
   endif

!----------------------------------------------------------------------
!
!  test minloc
!
!----------------------------------------------------------------------

   i4Expect = minval(i4ArrayG, mask=(i4ArrayG /= 0_POP_i4))
   r4Expect = minval(r4ArrayG, mask=(r4ArrayG /= 0.0_POP_r4))
   r8Expect = minval(r8ArrayG, mask=(r8ArrayG /= 0.0_POP_r8))
   iLoc = minloc(i4ArrayG, mask=(i4ArrayG /= 0_POP_i4))

   call POP_GlobalMinloc(i4Array2D, distrbClinic, &
                         i, j, i4Test, errorCode, lMask = mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 array minloc')
   endif
   if (i4Test /= i4Expect .or. &
       iLoc(1) /= i .or. iLoc(2) /= j) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 array minloc')
   endif

   call POP_GlobalMinloc(r4Array2D, distrbClinic, &
                         i, j, r4Test, errorCode, lMask = mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 array minloc')
   endif
   if (r4Test /= r4Expect .or. &
       iLoc(1) /= i .or. iLoc(2) /= j) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 array minloc')
   endif

   call POP_GlobalMinloc(r8Array2D, distrbClinic, &
                         i, j, r8Test, errorCode, lMask = mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 array minloc')
   endif
   if (r8Test /= r8Expect .or. &
       iLoc(1) /= i .or. iLoc(2) /= j) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 array minloc')
   endif

!----------------------------------------------------------------------
!
!  now test global count functions
!
!----------------------------------------------------------------------

   !*** test all kinds at cell center

   i4Expect = count(i4ArrayG /= 0_POP_i4)

   i4Test = POP_GlobalCount(mask, distrbClinic, &
                            POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in logical count')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in logical count')
   endif

   i4Test = POP_GlobalCount(i4Array2D, distrbClinic, &
                            POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 count')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 count')
   endif

   i4Test = POP_GlobalCount(r4Array2D, distrbClinic, &
                            POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 count')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 count')
   endif

   i4Test = POP_GlobalCount(r8Array2D, distrbClinic, &
                            POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 count')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 count')
   endif

   !*** test r8 function at NFace and NEcorner for tripole

   i4Expect = count(r8ArrayG /= 0.0_POP_r8)

   j=size(r8ArrayG, dim=2)
   ie = size(r8ArrayG, dim=1)
   ib = ie/2 + 1
   do i=ib,ie
      if (r8ArrayG(i,j) /= 0.0_POP_r8) i4Expect = i4Expect - 1
   end do

   i4Test = POP_GlobalCount(r8Array2D, distrbClinic, &
                            POP_gridHorzLocNFace, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 NFace count')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 NFace count')
   endif

   i4Test = POP_GlobalCount(r8Array2D, distrbClinic, &
                            POP_gridHorzLocNECorner, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 NECorner count')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 NECorner count')
   endif

!----------------------------------------------------------------------
!
!  now test global sum operations
!  first test without masks
!
!----------------------------------------------------------------------

   i4Expect = sum(i4ArrayG)
   r4Expect = sum(r4ArrayG)
   r8Expect = sum(r8ArrayG)

   i4Test = POP_GlobalSum(i4Array2D, distrbClinic,        &
                          POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 sum')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 sum')
   endif

   r4Test = POP_GlobalSum(r4Array2D, distrbClinic,        &
                          POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 sum')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 sum')
   endif

   r8Test = POP_GlobalSum(r8Array2D, distrbClinic,        &
                          POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 sum')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 sum')
   endif

!----------------------------------------------------------------------
!
!  now test with masks and test the sum with product - all should give
!  same answer
!
!----------------------------------------------------------------------

   where (mask)
      i4Mask = 1_POP_i4
      r4Mask = 1.0_POP_r4
      r8Mask = 1.0_POP_r8
   elsewhere
      i4Mask = 0_POP_i4
      r4Mask = 0.0_POP_r4
      r8Mask = 0.0_POP_r8
   end where

   !*** test with logical mask

   i4Test = POP_GlobalSum(i4Array2D, distrbClinic,          &
                          POP_gridHorzLocCenter, errorCode, &
                          lMask=mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 sum with log mask')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 sum with log mask')
   endif

   r4Test = POP_GlobalSum(r4Array2D, distrbClinic,          &
                          POP_gridHorzLocCenter, errorCode, &
                          lMask=mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 sum with log mask')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 sum with log mask')
   endif

   r8Test = POP_GlobalSum(r8Array2D, distrbClinic,          &
                          POP_gridHorzLocCenter, errorCode, &
                          lMask=mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 sum with log mask')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 sum with log mask')
   endif

   !*** test with multiplicative mask

   i4Test = POP_GlobalSum(i4Array2D, distrbClinic,          &
                          POP_gridHorzLocCenter, errorCode, &
                          mMask = i4Mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 sum with mult mask')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 sum with mult mask')
   endif

   r4Test = POP_GlobalSum(r4Array2D, distrbClinic,          &
                          POP_gridHorzLocCenter, errorCode, &
                          mMask = r4Mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 sum with mult mask')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 sum with mult mask')
   endif

   r8Test = POP_GlobalSum(r8Array2D, distrbClinic,          &
                          POP_gridHorzLocCenter, errorCode, &
                          mMask = r8Mask)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 sum with mult mask')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 sum with mult mask')
   endif

   !*** now sum with product

   i4Test = POP_GlobalSumProd(i4Array2D, i4Mask, distrbClinic, &
                              POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in i4 sum prod')
   endif
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in i4 sum prod')
   endif

   r4Test = POP_GlobalSumProd(r4Array2D, r4Mask, distrbClinic, &
                              POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r4 sum prod')
   endif
   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r4 sum prod')
   endif

   r8Test = POP_GlobalSumProd(r8Array2D, r8Mask, distrbClinic, &
                              POP_gridHorzLocCenter, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 sum prod')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 sum prod')
   endif

!----------------------------------------------------------------------
!
!  finally, test sum function on tripole NFace and NECorner
!
!----------------------------------------------------------------------

   j=size(r8ArrayG, dim=2)
   ie = size(r8ArrayG, dim=1)
   ib = ie/2 + 1
   do i=ib,ie
      r8Expect = r8Expect - r8ArrayG(i,j)
   end do

   r8Test = POP_GlobalSum(r8Array2D, distrbClinic, &
                          POP_gridHorzLocNFace, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 NFace sum')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 NFace sum')
   endif

   r8Test = POP_GlobalSum(r8Array2D, distrbClinic, &
                          POP_gridHorzLocNECorner, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: error in r8 NECorner sum')
   endif
   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ReductionTest: bad answer in r8 NECorner sum')
   endif

!----------------------------------------------------------------------
!
!  clean up
!
!----------------------------------------------------------------------

   call POP_ErrorPrint(errorCode)
   call POP_CommExitMessageEnvironment

!----------------------------------------------------------------------

 end program POP_ReductionTest

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
