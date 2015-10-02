!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 program POP_BlockDistributionTest

!----------------------------------------------------------------------
!
!  this program tests how POP distributes blocks based on an input
!  topography and processor count
!
!----------------------------------------------------------------------

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_IOUnitsMod
   use POP_DomainSizeMod
   use POP_BlocksMod
   use POP_DistributionMod

   use netcdf

   implicit none

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   !integer :: recordLength

   integer (POP_i4) :: &
      errorCode,       &! error flag
      istat,           &! allocate error flag
      nProcsClinic,    &! number of processors in clinic distribution
      nProcsTropic,    &! number of processors in tropic distribution
      blockID,         &! global block id
      ib,ie,jb,je,     &! beg, end of physical domain
      iunit,           &! IO unit
      i,j,n             ! loop indices

   type (POP_distrb) :: &
      distrbClinic,     &! baroclinic distribution
      distrbTropic       ! barotropic distribution

   character (POP_charLength) :: &
      topographyFile,     &! filename (and path) for input topo file
      distrbMethodClinic, &! choice for distributing clinic blocks
      distrbMethodTropic   ! choice for distributing tropic blocks

   integer (POP_i4), dimension(:), allocatable :: &
      workPerBlock     ! work per block

   integer (POP_i4), dimension(:), pointer :: &
      blockLocClinic, &! location of each block in clinic distrb
      blockLocTropic   ! location of each block in tropic distrb

   integer (POP_i4), dimension(:), pointer :: &
      iGlobal       ! global index for this block

   integer (POP_i4), dimension(:), pointer :: &
      jGlobal       ! global index for this block

   integer (POP_i4), dimension(:,:), allocatable :: &
      kmt          ! global topography

   namelist /domainNML/ nprocsClinic, nprocsTropic, topographyFile, &
                        distrbMethodClinic, distrbMethodTropic

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
!  read input namelist for distribution options and processor counts
!
!----------------------------------------------------------------------

   call POP_IOUnitsGet(iunit)
   open(iunit, file='pop_in', form='formatted', status='old')
   read(iunit, nml=domainNML)
   close(iunit)
   call POP_IOUnitsRelease(iunit)

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

!----------------------------------------------------------------------
!
!  read in topography and compute work per block
!
!----------------------------------------------------------------------

   allocate(kmt(POP_nxGlobal,POP_nyGlobal), stat=istat)
   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'distTest: error allocating kmt')
   endif

   kmt = 0

   !inquire (IOLENGTH=recordLength) kmt

   open(15, file=trim(topographyFile), form='unformatted', &
        convert='big_endian', access='direct', status='old', &
        recl=POP_nxGlobal*POP_nyGlobal*4)
        !recl=istat)
   read(15,rec=1) kmt
   close(15)

   !*** create work per block based on input topography

   allocate(workPerBlock(POP_numBlocks), stat=istat)
   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'BlockDistributionTest: error allocating work per block')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   do blockID=1,POP_numBlocks
      call POP_BlocksGetBlockInfo(blockID, errorCode,         &
                                  ib=ib, ie=ie, jb=jb, je=je, &
                                  iGlobal=iGlobal, jGlobal=jGlobal)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'BlockDistributionTest: error getting block info')
         call POP_ErrorPrint(errorCode)
         stop
      endif

      workPerBlock(blockID) = 0
      do j=jb,je
      do i=ib,ie
         if (kmt(iGlobal(i),jGlobal(j)) /= 0) then
            workPerBlock(blockID) = workPerBlock(blockID) + 1
         endif
      end do
      end do

      !*** rescale so that all ocean blocks have same
      !*** amount of work

      where (workPerBlock /= 0) workPerBlock = 10
   end do

   deallocate(kmt)

!----------------------------------------------------------------------
!
!  create block distributions
!
!----------------------------------------------------------------------

   !*** create baroclinic distribution

   select case (distrbMethodClinic)
   case ('rake')
      distrbClinic = POP_DistributionCreate(POP_distributionMethodRake, &
                                            nProcsClinic, &
                                            workPerBlock, errorCode)
   case ('cartesian')
      distrbClinic = POP_DistributionCreate(POP_distributionMethodCartesian, &
                                            nProcsClinic, &
                                            workPerBlock, errorCode)
   case default
      call POP_ErrorSet(errorCode, &
         'BlockDistributionTest: unknown clinic distribution method')
      call POP_ErrorPrint(errorCode)
   end select

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BlockDistributionTest: error creating clinic distribution')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   !*** create barotropic distribution

   select case (distrbMethodTropic)
   case ('rake')
      distrbTropic = POP_DistributionCreate(POP_distributionMethodRake, &
                                            nProcsTropic, &
                                            workPerBlock, errorCode)
   case ('cartesian')
      distrbTropic = POP_DistributionCreate(POP_distributionMethodCartesian, &
                                            nProcsTropic, &
                                            workPerBlock, errorCode)
   case default
      call POP_ErrorSet(errorCode, &
         'BlockDistributionTest: unknown tropic distribution method')
      call POP_ErrorPrint(errorCode)
   end select


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BlockDistributionTest: error creating tropic distribution')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   deallocate(workPerBlock)

!----------------------------------------------------------------------
!
!  find block locations
!
!----------------------------------------------------------------------

   call POP_DistributionGet(distrbClinic, errorCode,            &
                            blockLocation = blockLocClinic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BlockDistributionTest: error retrieving clinic locations')
      call POP_ErrorPrint(errorCode)
      stop
   endif

   call POP_DistributionGet(distrbTropic, errorCode,            &
                            blockLocation = blockLocTropic)
                                                                                
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BlockDistributionTest: error retrieving tropic locations')
      call POP_ErrorPrint(errorCode)
      stop
   endif

!----------------------------------------------------------------------
!
!  write out distribution information
!
!----------------------------------------------------------------------

   write(POP_stdout,'(a23)') 'Baroclinic Distribution'
   write(POP_stdout,"('blockID',2x,'blockLoc')")
   write(POP_stdout,"('-------',2x,'--------')")
   do n=1,POP_numBlocks
      write(POP_stdout,'(i7,2x,i8)') n,blockLocClinic(n)
   end do

   write(POP_stdout,'(a23)') 'Barotropic Distribution'
   write(POP_stdout,"('blockID',2x,'blockLoc')")
   write(POP_stdout,"('-------',2x,'--------')")
   do n=1,POP_numBlocks
      write(POP_stdout,'(i7,2x,i8)') n,blockLocTropic(n)
   end do

!----------------------------------------------------------------------
!
!  clean up and exit
!
!----------------------------------------------------------------------

   call POP_ErrorPrint(errorCode)
   call POP_CommExitMessageEnvironment

!----------------------------------------------------------------------

 end program POP_BlockDistributionTest

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
