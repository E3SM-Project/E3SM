!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 program POP_RedistributeTest

!----------------------------------------------------------------------
!
!  this program tests POP redistribution communications, including
!  gather/scatter operations and moving blocks between distributions
!
!----------------------------------------------------------------------

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_DomainSizeMod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_RedistributeMod

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
      nProcsTropic,    &! number of processors in tropic distribution
      nBlocksClinic,   &! number of local blocks in clinic distribution
      nBlocksTropic,   &! number of local blocks in tropic distribution
      blockID,         &! global block id
      ib,ie,jb,je,     &! beg, end of physical domain
      i,j,n             ! loop indices

   type (POP_distrb) :: &
      distrbClinic,     &! baroclinic distribution
      distrbTropic       ! barotropic distribution

   integer (POP_i4), dimension(:), allocatable :: &
      workPerBlock  ! test values

   integer (POP_i4), dimension(:), pointer :: &
      iGlobal       ! global index for this block

   integer (POP_i4), dimension(:), pointer :: &
      jGlobal       ! global index for this block

   integer (POP_i4), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksClinic) :: &
      i4Test, i4Expect  ! test values

   real (POP_r4), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksClinic) :: &
      r4Test, r4Expect  ! test values

   real (POP_r8), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksClinic) :: &
      r8Test, r8Expect  ! test values

   integer (POP_i4), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksTropic) :: &
      i4TestTropic, i4ExpectTropic  ! test values

   real (POP_r4), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksTropic) :: &
      r4TestTropic, r4ExpectTropic  ! test values

   real (POP_r8), &
      dimension(POP_nxBlock, POP_nyBlock, POP_maxBlocksTropic) :: &
      r8TestTropic, r8ExpectTropic  ! test values

   integer (POP_i4), &
      dimension(POP_nxGlobal, POP_nyGlobal) :: &
      i4TestGlobal, i4ExpectGlobal  ! test values

   real (POP_r4), &
      dimension(POP_nxGlobal, POP_nyGlobal) :: &
      r4TestGlobal, r4ExpectGlobal  ! test values

   real (POP_r8), &
      dimension(POP_nxGlobal, POP_nyGlobal) :: &
      r8TestGlobal, r8ExpectGlobal  ! test values

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
                         'cyclic', 'closed', errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error creating block decomposition')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   !*** initialize global arrays with land in the middle

   do j=1,POP_nyGlobal
   do i=1,POP_nxGlobal
      i4TestGlobal(i,j) = 0_POP_i4
      r4TestGlobal(i,j) = 0.0_POP_r4
      r8TestGlobal(i,j) = 0.0_POP_r8
      if (i > POP_nxGlobal/2 + POP_nxBlock .or. &
          i < POP_nxGlobal/2 - POP_nxBlock .or. & 
          j > POP_nyGlobal/2 + POP_nyBlock .or. &
          j < POP_nyGlobal/2 - POP_nyBlock) then
         i4ExpectGlobal(i,j) = (i+j)*10_POP_i4
         r4ExpectGlobal(i,j) = (i+j)*100.0_POP_r4
         r8ExpectGlobal(i,j) = (i+j)*1000.0_POP_r8
      else
         i4ExpectGlobal(i,j) = 0_POP_i4
         r4ExpectGlobal(i,j) = 0.0_POP_r4
         r8ExpectGlobal(i,j) = 0.0_POP_r8
      endif
   end do
   end do

   !*** create artificial work per block with varying work and
   !*** at least one empty block 

   allocate(workPerBlock(POP_numBlocks), stat=istat)
   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error allocating work per block')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   do blockID=1,POP_numBlocks
      call POP_BlocksGetBlockInfo(blockID, errorCode,         &
                                  ib=ib, ie=ie, jb=jb, je=je, &
                                  iGlobal=iGlobal, jGlobal=jGlobal)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: error getting block info')
         call POP_ErrorPrint(errorCode)
         stop
      endif

      workPerBlock(blockID) = 0
      do j=jb,je
      do i=ib,ie
         if (i4ExpectGlobal(iGlobal(i),jGlobal(j)) /= 0) then
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
         'RedistTest: error creating clinic distribution')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   !*** create barotropic distribution

   nprocsTropic = max(1,nprocsClinic/2)

   distrbTropic = POP_DistributionCreate(POP_distributionMethodCartesian, &
                                         nProcsTropic, &
                                         workPerBlock, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error creating tropic distribution')
      call POP_ErrorPrint(errorCode)
      stop
   endif

!----------------------------------------------------------------------
!
!  initialize test arrays
!
!----------------------------------------------------------------------

   call POP_DistributionGet(distrbClinic, errorCode,            &
                            numLocalBlocks = nBlocksClinic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error getting nBlocksClinic')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   call POP_DistributionGet(distrbTropic, errorCode,            &
                            numLocalBlocks = nBlocksTropic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error getting nBlocksTropic')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   !*** initialize clinic arrays

   do n=1,nBlocksClinic

      call POP_DistributionGetBlockID(distrbClinic, n, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: error getting blockID for clinic')
         call POP_ErrorPrint(errorCode)
         stop
      endif

      call POP_BlocksGetBlockInfo(blockID, errorCode,         &
                                  ib=ib, ie=ie, jb=jb, je=je, &
                                  iGlobal=iGlobal, jGlobal=jGlobal)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: error getting block info clinic')
         call POP_ErrorPrint(errorCode)
         stop
      endif

      i4Test(:,:,n) = 0
      r4Test(:,:,n) = 0.0_POP_r4
      r8Test(:,:,n) = 0.0_POP_r8
      i4Expect(:,:,n) = 0
      r4Expect(:,:,n) = 0.0_POP_r4
      r8Expect(:,:,n) = 0.0_POP_r8

      do j=jb,je
      do i=ib,ie
         if (i4ExpectGlobal(iGlobal(i),jGlobal(j)) /= 0) then
            i4Expect(i,j,n) = (iGlobal(i)+jGlobal(j))*10_POP_i4
            r4Expect(i,j,n) = (iGlobal(i)+jGlobal(j))*100.0_POP_r4
            r8Expect(i,j,n) = (iGlobal(i)+jGlobal(j))*1000.0_POP_r8
         endif
      end do
      end do
   end do

   !*** initialize tropic arrays

   do n=1,nBlocksTropic

      call POP_DistributionGetBlockID(distrbTropic, n, &
                                      blockID, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: error getting blockID for tropic')
         call POP_ErrorPrint(errorCode)
         stop
      endif

      call POP_BlocksGetBlockInfo(blockID, errorCode,         &
                                  ib=ib, ie=ie, jb=jb, je=je, &
                                  iGlobal=iGlobal, jGlobal=jGlobal)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: error getting block info tropic')
         call POP_ErrorPrint(errorCode)
         stop
      endif

      do j=jb,je
      do i=ib,ie
         i4TestTropic(i,j,n) = 0
         r4TestTropic(i,j,n) = 0.0_POP_r4
         r8TestTropic(i,j,n) = 0.0_POP_r8
         i4ExpectTropic(i,j,n) = 0
         r4ExpectTropic(i,j,n) = 0.0_POP_r4
         r8ExpectTropic(i,j,n) = 0.0_POP_r8

         if (i4ExpectGlobal(iGlobal(i),jGlobal(j)) /= 0) then
            i4ExpectTropic(i,j,n) = (iGlobal(i)+jGlobal(j))*10_POP_i4
            r4ExpectTropic(i,j,n) = (iGlobal(i)+jGlobal(j))*100.0_POP_r4
            r8ExpectTropic(i,j,n) = (iGlobal(i)+jGlobal(j))*1000.0_POP_r8
         endif
      end do
      end do
   end do

!----------------------------------------------------------------------
!
!  test scatter functions
!
!----------------------------------------------------------------------

   !*** scatter global array to local arrays

   call POP_RedistributeScatter(i4Test, i4ExpectGlobal, &
                              POP_masterTask, distrbClinic, errorCode) 

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error in i4 scatter')
   endif

   call POP_RedistributeScatter(r4Test, r4ExpectGlobal, &
                              POP_masterTask, distrbClinic, errorCode) 

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error in r4 scatter')
   endif

   call POP_RedistributeScatter(r8Test, r8ExpectGlobal, &
                              POP_masterTask, distrbClinic, errorCode) 

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error in r8 scatter')
   endif

   !***  check for proper values

   if (count(i4Test /= i4Expect) > 0) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: bad values from i4 scatter')
   endif

   if (count(r4Test /= r4Expect) > 0) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: bad values from r4 scatter')
   endif

   if (count(r8Test /= r8Expect) > 0) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: bad values from r8 scatter')
   endif

!----------------------------------------------------------------------
!
!  now test gather operations
!
!----------------------------------------------------------------------

   call POP_RedistributeGather(i4TestGlobal, i4Expect, &
                               POP_masterTask, distrbClinic, &
                               errorCode, fillValue = 0_POP_i4)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error in i4 gather')
   endif

   call POP_RedistributeGather(r4TestGlobal, r4Expect, &
                               POP_masterTask, distrbClinic, &
                               errorCode, fillValue = 0.0_POP_r4)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error in r4 gather')
   endif

   call POP_RedistributeGather(r8TestGlobal, r8Expect, &
                               POP_masterTask, distrbClinic, &
                               errorCode, fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error in r8 gather')
   endif

   !*** check values

   if (POP_myTask == POP_masterTask) then
      if (count(i4TestGlobal /= i4ExpectGlobal) > 0) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: bad values from i4 gather')
      endif

      if (count(r4TestGlobal /= r4ExpectGlobal) > 0) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: bad values from r4 gather')
      endif

      if (count(r8TestGlobal /= r8ExpectGlobal) > 0) then
         call POP_ErrorSet(errorCode, &
            'RedistTest: bad values from r8 gather')
      endif
   endif

!----------------------------------------------------------------------
!
!  now test block redistribution
!
!----------------------------------------------------------------------

   call POP_RedistributeBlocks(i4TestTropic, distrbTropic, &
                               i4Expect,     distrbClinic, &
                               errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error redistributing i4 blocks')
   endif

   call POP_RedistributeBlocks(r4TestTropic, distrbTropic, &
                               r4Expect,     distrbClinic, &
                               errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error redistributing r4 blocks')
   endif

   call POP_RedistributeBlocks(r8TestTropic, distrbTropic, &
                               r8Expect,     distrbClinic, &
                               errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: error redistributing r8 blocks')
   endif

   !*** check values

   if (count(i4TestTropic /= i4ExpectTropic) > 0) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: bad values from i4 block redist')
   endif

   if (count(r4TestTropic /= r4ExpectTropic) > 0) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: bad values from r4 block redist')
   endif

   if (count(r8TestTropic /= r8ExpectTropic) > 0) then
      call POP_ErrorSet(errorCode, &
         'RedistTest: bad values from r8 block redist')
   endif

!----------------------------------------------------------------------
!
!  clean up
!
!----------------------------------------------------------------------

   call POP_ErrorPrint(errorCode)
   call POP_CommExitMessageEnvironment

!----------------------------------------------------------------------

 end program POP_RedistributeTest

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
