!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_SolversMod

!BOP
! !MODULE: POP_SolversMod
!
! !DESCRIPTION:
!  This module contains routines and operators for solving the elliptic
!  system for surface pressure in the barotropic mode.
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_ConfigMod
   use POP_IOUnitsMod
   use POP_BlocksMod
   use POP_DistributionMod
   use POP_GridHorzMod
   use POP_ReductionsMod
   use POP_RedistributeMod
   use POP_HaloMod
   use POP_GridHorzMod
   use POP_FieldMod
   use POP_DomainSizeMod
   use domain
   use grid
   use perf_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_SolversInit,     &
             POP_SolversRun,      &
             POP_SolversDiagonal, &
             POP_SolversGetDiagnostics

! !PUBLIC DATA MEMBERS:


!-----------------------------------------------------------------------
!
!  other operator and preconditioner weights for barotropic operator
!
!-----------------------------------------------------------------------

   real (POP_r8), dimension (:,:,:), allocatable, public :: & 
      mMaskTropic,     &! land mask in barotropic distribution 
      btropWgtCenter,  &! barotropic operater center coefficient
      btropWgtNorth,   &! barotropic operater north  coefficient
      btropWgtEast,    &! barotropic operater east   coefficient
      btropWgtNE        ! barotropic operater northeast coefficient

   real (POP_r8), dimension (:,:,:), allocatable, public :: & 
      centerWgtClinicIndep, &! time indep  center wgt on clinic distrb
      centerWgtClinic        ! time depend center wgt on clinic distrb  

!EOP
!BOC
   !*** preconditioner operator coefficients (ninept operator)

   real (POP_r8), dimension (:,:,:), allocatable :: & 
      precondCenter,              &
      precondNorth, precondSouth, &
      precondEast,  precondWest,  &
      precondNE,    precondSE,    &
      precondNW,    precondSW

!-----------------------------------------------------------------------
!
!  supported solvers and preconditioners
!
!-----------------------------------------------------------------------

   character (POP_charLength) :: &
      solverChoice

   character (3), parameter :: &
      solverChoicePCG = 'pcg'
   character (9), parameter :: &
      solverChoiceChronGear = 'ChronGear'

   character (POP_charLength) :: &
      preconditionerChoice

   character (8), parameter :: &
      precondChoiceDiag = 'diagonal'
   character (4), parameter :: &
      precondChoiceFile = 'file'

   logical (POP_logical) :: &
      usePreconditioner

!-----------------------------------------------------------------------
!
!  scalar convergence-related variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::      &
      maxIterations,        &! max number of solver iterations
      convergenceCheckFreq   ! check convergence every freq steps

   real (POP_r8) ::         &
      convergenceCriterion, &! convergence error criterion
      residualNorm           ! residual normalization

   !*** convergence diagnostics

   integer (POP_i4), public ::      &
      numIterations          ! accumulated no of iterations (diagnostic)

   real (POP_r8) ::         &
      rmsResidual            ! residual (also a diagnostic)

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversRun
! !INTERFACE:

 subroutine POP_SolversRun(sfcPressure, rhsClinic, errorCode)

! !DESCRIPTION:
!  Solves the elliptic equation for surface pressure by calling
!  the requested solver routine.  Also redistributes necessary
!  array to the barotropic distribution of blocks for better performance
!  of the solver.
!  The elliptic equation is
!  \begin{equation}
!     AF = B
!  \end{equation}
!  where $F$ is a field (eg surface pressure), $B$ is the right hand side
!  and $A$ is the operator defined as
!  \begin{equation}
!     AF = a \nabla\cdot(H \nabla F)
!  \end{equation}
!  where $a$ is the cell area.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      rhsClinic         ! right-hand-side of linear system
                        !  for blocks in baroclinic distribution

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
     sfcPressure              ! on input,  initial guess in baroclinic distrb
                        ! on output, final solution for sfc pressure

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      numProcs           ! number of processors in barotropic distrib

   real (POP_r8), dimension(size(sfcPressure,dim=1), &
                            size(sfcPressure,dim=2), &
                            size(sfcPressure,dim=3)) :: &
      pressTropic,     &! surface pressure in barotropic distribution
      rhsTropic         ! right hand side  in barotropic distribution

!-----------------------------------------------------------------------
!
!  switch to the barotropic distribution for iterative solvers
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_RedistributeBlocks(btropWgtCenter,  POP_distrbTropic, &
                               centerWgtClinic, POP_distrbClinic, errorCode)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing center operator weight')
      return
   endif

   call POP_RedistributeBlocks(pressTropic, POP_distrbTropic, &
                               sfcPressure, POP_distrbClinic, errorCode)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing pressure')
      return
   endif

   call POP_RedistributeBlocks(rhsTropic, POP_distrbTropic, &
                               rhsClinic, POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing right hand side')
      return
   endif

!-----------------------------------------------------------------------
!
!  call proper routine based on user choice of solver
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numProcs = numProcs)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error getting num procs')
      return
   endif

   if (POP_myTask < numProcs) then
      select case(trim(solverChoice))
      case (solverChoicePCG)

         call pcg(pressTropic, rhsTropic, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolverRun: error in PCG')
            return
         endif

      case (solverChoiceChronGear)

         call ChronGear(pressTropic, rhsTropic, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolverRun: error in ChronGear')
            return
         endif

      end select
   endif

!-----------------------------------------------------------------------
!
!  switch solution back to the baroclinic distribution
!
!-----------------------------------------------------------------------

   call POP_RedistributeBlocks(sfcPressure, POP_distrbClinic, &
                               pressTropic, POP_distrbTropic, errorCode)



   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversRun: error redistributing pressure back')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversRun

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversInit
! !INTERFACE:

 subroutine POP_SolversInit(errorCode)

! !DESCRIPTION:
!  This routine initializes choice of solver, calculates the 
!  coefficients of the 9-point stencils for the barotropic operator and
!  reads in a preconditioner if requested.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode       ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!       {X,Y}{NE,SE,NW,SW} = contribution to {ne,se,nw,sw} coefficients 
!         from {x,y} components of divergence
!       HU = depth at U points
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,             &! dummy loop counters
      configUnit,      &! unit for configuration file
      numBlocksClinic, &! num local blocks in baroclinic distribution
      numBlocksTropic, &! num local blocks in barotropic distribution
      iblock,          &! block counter
      istat             ! status flag for allocates

   character (POP_charLength) :: &
      preconditionerFile  ! file containing preconditioner

   real (POP_r8) ::    &
      xne,xse,xnw,xsw, &! contribution to coefficients from x,y
      yne,yse,ynw,ysw, &!   components of divergence
      ase,anw,asw

   real (POP_r8), dimension(:,:,:), allocatable :: &
      work0,           &! temp space for computing barotropic
      workNorth,       &!
      workEast,        &!
      workNE,          &!
      mMaskTmp          ! land mask in barotropic distribution


!-----------------------------------------------------------------------
!
!  read solver choice and solver constants from configuration file
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (POP_myTask == POP_masterTask) then
      write(POP_stdout,POP_blankFormat)
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,'(a35)') ' Solver options (barotropic solver)'
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,POP_blankFormat)
   endif

   call POP_ConfigOpen(configUnit, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error opening config file')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'solverChoice',    &
   call POP_ConfigRead(configUnit, 'solvers', 'solverchoice',    &
                       solverChoice, solverChoicePCG, errorCode, &
                       outStringBefore = 'Solver choice: ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading solver choice')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'convergenceCriterion', &
   call POP_ConfigRead(configUnit, 'solvers', 'convergencecriterion', &
                   convergenceCriterion, 1.e-12_POP_r8, errorCode,    &
                   outStringBefore = 'Solver converged for err < ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading solver convergence criterion')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'maxIterations',    &
   call POP_ConfigRead(configUnit, 'solvers', 'maxiterations',    &
                       maxIterations, 1000, errorCode, &
                       outStringBefore = 'Solver maximum iterations: ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading solver max iterations')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'convergenceCheckFreq', &
   call POP_ConfigRead(configUnit, 'solvers', 'convergencecheckfreq', &
                       convergenceCheckFreq, 10, errorCode,           &
                       outStringBefore = 'Check convergence every ',  &
                       outStringAfter  = ' iterations')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading convergence check frequency')
      return
   endif

!  call POP_ConfigRead(configUnit, 'solvers', 'preconditionerChoice', &
   call POP_ConfigRead(configUnit, 'solvers', 'preconditionerchoice', &
                       preconditionerChoice, precondChoiceDiag,       &
                       errorCode,                                     &
                       outStringBefore = 'Preconditioner choice: ')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error reading preconditioner choice')
      return
   endif

   if (trim(preconditionerChoice) == precondChoiceFile) then
!     call POP_ConfigRead(configUnit, 'solvers', 'preconditionerFile', &
      call POP_ConfigRead(configUnit, 'solvers', 'preconditionerfile', &
                 preconditionerFile, 'UnknownPrecondFile', errorCode,  &
                 outStringBefore = 'Reading preconditioner from file: ')

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversInit: error reading preconditioner file name')
         return
      endif
   endif

   call POP_ConfigClose(configUnit, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error closing config file')
      return
   endif

!-----------------------------------------------------------------------
!
!  check inputs
!
!-----------------------------------------------------------------------

   select case(trim(solverChoice))
   case(solverChoicePCG)
   case(solverChoiceChronGear)
   case default
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: unknown solver - must be pcg, ChronGear')
      return
   end select

   select case (trim(preconditionerChoice))
   case(precondChoiceDiag)
      usePreconditioner = .false.   ! default is diagonal
   case(precondChoiceFile)
      usePreconditioner = .true.
   case default
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: unknown preconditioner choice')
      return
   end select

!-----------------------------------------------------------------------
!
!  compute nine point operator coefficients: compute on baroclinic
!  decomposition first where grid info defined and redistribute
!  to barotropic distribution
!  leave center coefficients in baroclinic distribution to facilitate 
!  easy time-dependent changes in barotropic routine
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbClinic, errorCode, &
                            numLocalBlocks = numBlocksClinic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error retrieving clinic local block count')
      return
   endif

   allocate(work0     (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workNorth (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workEast  (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            workNE    (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            mMaskTmp  (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
       centerWgtClinicIndep(POP_nxBlock,POP_nyBlock,numBlocksClinic), &
       centerWgtClinic     (POP_nxBlock,POP_nyBlock,numBlocksClinic), &
            stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error allocating temporary arrays')
      return
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,i,j,xne,xse,xnw,xsw,yne,yse,ynw,ysw, &
   !$OMP                     ase,anw,asw)

   do iblock = 1,numBlocksClinic

      work0                (:,:,iblock) = 0.0_POP_r8
      workNorth            (:,:,iblock) = 0.0_POP_r8
      workEast             (:,:,iblock) = 0.0_POP_r8
      workNE               (:,:,iblock) = 0.0_POP_r8
      mMaskTmp             (:,:,iblock) = 0.0_POP_r8
      centerWgtClinicIndep (:,:,iblock) = 0.0_POP_r8

      do j=2,POP_nyBlock
      do i=2,POP_nxBlock

         xne = 0.25_POP_r8*HU(i  ,j  ,iblock)*DXUR(i  ,j  ,iblock)* &
                                              DYU (i  ,j  ,iblock)
         xse = 0.25_POP_r8*HU(i  ,j-1,iblock)*DXUR(i  ,j-1,iblock)* &
                                              DYU (i  ,j-1,iblock)
         xnw = 0.25_POP_r8*HU(i-1,j  ,iblock)*DXUR(i-1,j  ,iblock)* &
                                              DYU (i-1,j  ,iblock)
         xsw = 0.25_POP_r8*HU(i-1,j-1,iblock)*DXUR(i-1,j-1,iblock)* &
                                              DYU (i-1,j-1,iblock)

         yne = 0.25_POP_r8*HU(i  ,j  ,iblock)*DYUR(i  ,j  ,iblock)* &
                                              DXU (i  ,j  ,iblock)
         yse = 0.25_POP_r8*HU(i  ,j-1,iblock)*DYUR(i  ,j-1,iblock)* &
                                              DXU (i  ,j-1,iblock)
         ynw = 0.25_POP_r8*HU(i-1,j  ,iblock)*DYUR(i-1,j  ,iblock)* &
                                              DXU (i-1,j  ,iblock)
         ysw = 0.25_POP_r8*HU(i-1,j-1,iblock)*DYUR(i-1,j-1,iblock)* &
                                              DXU (i-1,j-1,iblock)

         workNE(i,j,iblock) = xne + yne
         ase                = xse + yse
         anw                = xnw + ynw
         asw                = xsw + ysw
 
         workEast (i,j,iblock)  = xne + xse - yne - yse
         workNorth(i,j,iblock)  = yne + ynw - xne - xnw

         centerWgtClinicIndep(i,j,iblock) = &
                        -(workNE(i,j,iblock) + ase + anw + asw)

         work0   (i,j,iblock) = TAREA (i,j,iblock)**2
         mMaskTmp(i,j,iblock) = RCALCT(i,j,iblock)

      end do
      end do
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  redistribute operator weights and mask to barotropic distribution
!
!-----------------------------------------------------------------------

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocksTropic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error retrieving tropic local block count')
      return
   endif

   allocate(btropWgtCenter (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            btropWgtNorth  (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            btropWgtEast   (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            btropWgtNE     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            mMaskTropic    (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
            stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error allocating operator weights')
      return
   endif

   call POP_RedistributeBlocks(btropWgtNorth, POP_distrbTropic, &
                               workNorth,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error redistributing north operator weight')
      return
   endif

   call POP_RedistributeBlocks(btropWgtEast, POP_distrbTropic, &
                               workEast,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error redistributing east operator weight')
      return
   endif

   call POP_RedistributeBlocks(btropWgtNE, POP_distrbTropic, &
                               workNE,     POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error redistributing NE operator weight')
      return
   endif

   call POP_RedistributeBlocks(mMaskTropic, POP_distrbTropic, &
                               mMaskTmp,    POP_distrbClinic, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error redistributing mask')
      return
   endif

!-----------------------------------------------------------------------
!
!  calculate normalization constant (darea,darea) for rmsResidual
!  in cgr routine.
!
!-----------------------------------------------------------------------

   residualNorm = 1.0_POP_r8/POP_GlobalSum(work0, POP_distrbClinic, &
                                           POP_gridHorzLocCenter,   &
                                           errorCode,               &
                                           mMask = mMaskTmp)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error computing normalization factor')
      return
   endif

   convergenceCriterion = convergenceCriterion**2/residualNorm

   deallocate(mMaskTmp, stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error deallocating temp mask')
      return
   endif

!-----------------------------------------------------------------------
!
!  setup preconditioner if required
!
!-----------------------------------------------------------------------

   if (usePreconditioner) then

      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: preconditioner not supported')
      return

      allocate(                                                   &
         precondCenter (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondNorth  (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondSouth  (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondEast   (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondWest   (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondNE     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondSE     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondNW     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         precondSW     (POP_nxBlock,POP_nyBlock,numBlocksTropic), &
         stat = istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversInit: error allocating preconditioner')
         return
      endif

!****
!**** OLD CODE FOR READING PRECONDITIONER - NEEDS TO BE REVISED
!****
!     if (POP_myTask == POP_masterTask) then
!       write(stdout,*) ' Preconditioner read from file: ', &
!                         trim(precond_file)
!     endif
!
!     allocate(WORKS (POP_nxBlock,POP_nyBlock,nblocks_clinic), &
!              WORKW (POP_nxBlock,POP_nyBlock,nblocks_clinic), &
!              WORKSE(POP_nxBlock,POP_nyBlock,nblocks_clinic), &
!              WORKSW(POP_nxBlock,POP_nyBlock,nblocks_clinic))
!
!     allocate(icheck(nblocks_clinic))
!
!-----------------------------------------------------------------------
!
!    read preconditioner and check that it is consistent with
!    KMU field
!
!-----------------------------------------------------------------------
!
!     call open_parallel_file(nu,precond_file,recl_dbl)
!     call read_array(nu,workC)
!     call read_array(nu,workN)
!     call read_array(nu,WORKS)
!     call read_array(nu,workE)
!     call read_array(nu,WORKW)
!     call read_array(nu,workNE)
!     call read_array(nu,workNW)
!     call read_array(nu,WORKSE)
!     call read_array(nu,WORKSW)
!     call close_parallel_file(nu)
!
!     if (POP_myTask == POP_masterTask) then
!       write(stdout,blank_fmt)
!       write(stdout,*) ' file read: ', trim(precond_file)
!     endif
!
!-----------------------------------------------------------------------
!
!    check that PC is consistent with KMU field
!
!-----------------------------------------------------------------------
!
!     do iblock = 1,nblocks_clinic
!
!       icheck(iblock) = 0
!
!       do j=1,POP_nyBlock
!       do i=1,POP_nxBlock
!
!         mlandne = .false.
!         mlandnw = .false.
!         mlandse = .false.
!         mlandsw = .false.
!         if (KMU(i  ,j  ,iblock) == 0) mlandne = .true.
!         if (KMU(i-1,j  ,iblock) == 0) mlandnw = .true.
!         if (KMU(i  ,j-1,iblock) == 0) mlandse = .true.
!         if (KMU(i-1,j-1,iblock) == 0) mlandsw = .true.
!
!         if (mlandne .and. workNE(i,j,iblock) /= 0.0_POP_r8)  &
!                           icheck(iblock) = icheck(iblock) + 1
!
!         if (mlandnw .and. workNW(i,j,iblock) /= 0.0_POP_r8)  &
!                           icheck(iblock) = icheck(iblock) + 1
!
!         if (mlandse .and. WORKSE(i,j,iblock) /= 0.0_POP_r8)  &
!                           icheck(iblock) = icheck(iblock) + 1
!
!         if (mlandsw .and. WORKSW(i,j,iblock) /= 0.0_POP_r8)  &
!                           icheck(iblock) = icheck(iblock) + 1
!      
!         if (mlandne .and. mlandnw .and. (workN(i,j,iblock) /= 0.0_POP_r8)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandne .and. mlandse .and. (workE(i,j,iblock) /= 0.0_POP_r8)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandnw .and. mlandsw .and. (WORKW(i,j,iblock) /= 0.0_POP_r8)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandse .and. mlandsw .and. (WORKS(i,j,iblock) /= 0.0_POP_r8)) &
!                           icheck(iblock) = icheck(iblock) + 1
!         if (mlandne .and. mlandse .and.                            &
!             mlandnw .and. mlandsw .and. (workC(i,j,iblock) /= 0.0_POP_r8)) &
!                           icheck(iblock) = icheck(iblock) + 1
!       end do
!       end do
!     end do
!
!     ncheck = sum(icheck)
!     if (POP_GlobalSum(ncheck, POP_distrbClinic, errorCode) /= 0) then
!        call POP_ErrorSet(errorCode, &
!           'POP_SolversInit: PC and KMU are incompatible')
!        return
!     endif
!
!     deallocate(icheck)
!
!     call redistribute_blocks(precondCenter ,distrb_tropic,workC ,distrb_clinic)
!     call redistribute_blocks(precondNorth ,distrb_tropic,workN ,distrb_clinic)
!     call redistribute_blocks(precondEast ,distrb_tropic,workE ,distrb_clinic)
!     call redistribute_blocks(precondSouth ,distrb_tropic,WORKS ,distrb_clinic)
!     call redistribute_blocks(precondWest ,distrb_tropic,WORKW ,distrb_clinic)
!     call redistribute_blocks(precondNE,distrb_tropic,workNE,distrb_clinic)
!     call redistribute_blocks(precondNW,distrb_tropic,workNW,distrb_clinic)
!     call redistribute_blocks(precondSE,distrb_tropic,WORKSE,distrb_clinic)
!     call redistribute_blocks(precondSW,distrb_tropic,WORKSW,distrb_clinic)
!
!     deallocate(WORKS, WORKW, workNW, WORKSE, WORKSW)
!
   endif

!-----------------------------------------------------------------------
!
!  clean up temporary arrays
!
!-----------------------------------------------------------------------

   deallocate(work0, workNorth, workEast, workNE, stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversInit: error deallocating temp mask')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversInit

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversDiagonal
! !INTERFACE:

 subroutine POP_SolversDiagonal(diagonalCorrection, blockIndx, errorCode)

! !DESCRIPTION:
!  This routine corrects the center barotropic operator by 
!  subtracting off the time dependent diagnonal term computed in
!  the barotropic driver.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:), intent(in) :: &
      diagonalCorrection   ! time dependent diagonal term

   integer (POP_i4), intent(in) :: &
      blockIndx            ! local block index

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  subtract the time dependent diagonal term from the center operator
!  coefficient
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   centerWgtClinic(:,:,blockIndx) =                               &
                            centerWgtClinicIndep(:,:,blockIndx) - &
                              diagonalCorrection(:,:)
                     
!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversDiagonal

!***********************************************************************
!BOP
! !IROUTINE: POP_SolversGetDiagnostics
! !INTERFACE:

 subroutine POP_SolversGetDiagnostics(iterationCount, residual, &
                                      errorCode)

! !DESCRIPTION:
!  This routine returns the latest iteration count and residual
!  for diagnostic purposes.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   real (POP_r8), intent(out) :: &
      residual          ! latest residual from solver

   integer (POP_i4), intent(out) :: &
      iterationCount,  &! latest iteration count from solver
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  return the diagnostic quantities
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   residual = rmsResidual
   iterationCount = numIterations

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_SolversGetDiagnostics

!***********************************************************************
!BOP
! !IROUTINE: pcg
! !INTERFACE:

 subroutine pcg(X,B,errorCode)

! !DESCRIPTION:
!  This routine uses a preconditioned conjugate-gradient solver to
!  solve the equation $Ax=b$.  Both the operator $A$ and preconditioner
!  are nine-point stencils. If no preconditioner has been supplied,
!  a diagonal preconditioner is applied.  Convergence is checked
!  every {\em ncheck} steps.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
      X                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system

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
      i,j,m,           &! local iteration counter
      nx, ny,          &! horizontal extents
      numBlocks,       &! number of local blocks
      iblock            ! local block     counter

   real (POP_r8) ::          &
      eta0,eta1,rr        ! scalar inner product results

   real (POP_r8), dimension(size(X,dim=1),size(X,dim=2), &
                                          size(X,dim=3)) :: &
      R,                 &! residual (b-Ax)
      S,                 &! conjugate direction vector
      Q,work0,work1       ! various cg intermediate results

!-----------------------------------------------------------------------
!
!  compute initial residual and initialize S
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocks)

   nx = size(X,dim=1)
   ny = size(X,dim=2)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversPCG: error retrieving local block count')
      return
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,i,j)

   do iblock=1,numBlocks
      call btropOperator(S,X,iblock)
      do j=1,ny
      do i=1,nx
         R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock)
         S(i,j,iblock) = 0.0_POP_r8
      end do
      end do
   end do ! block loop

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  initialize fields and scalars
!
!-----------------------------------------------------------------------

   call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                       POP_fieldKindScalar, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversPCG: error updating initial residual halo')
      return
   endif

   eta0 =1.0_POP_r8 
   numIterations = maxIterations
 
!-----------------------------------------------------------------------
!
!  iterate
!
!-----------------------------------------------------------------------

   iterationLoop: do m = 1, maxIterations

!-----------------------------------------------------------------------
!
!     calculate (PC)r 
!     diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)

      do iblock=1,numBlocks

         if (usePreconditioner) then
            call preconditioner(work1,R,iblock)
         else
            do j=1,ny
            do i=1,nx
               if (btropWgtCenter(i,j,iblock) /= 0.0_POP_r8) then
                  work1(i,j,iblock) = R(i,j,iblock)/ &
                                      btropWgtCenter(i,j,iblock)
               else
                  work1(i,j,iblock) = 0.0_POP_r8
               endif
            end do
            end do
         endif

         do j=1,ny
         do i=1,nx
            work0(i,j,iblock) = R(i,j,iblock)*work1(i,j,iblock)
         end do
         end do
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     update conjugate direction vector s
!
!-----------------------------------------------------------------------

      if (usePreconditioner) then
         call POP_HaloUpdate(work1, POP_haloTropic, &
                             POP_gridHorzLocCenter, &
                             POP_fieldKindScalar, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversPCG: error updating halo for preconditioner')
            return
         endif
      endif

      !*** (r,(PC)r)
      eta1 = POP_GlobalSum(work0, POP_distrbTropic, &
                           POP_gridHorzLocCenter,   &
                           errorCode, mMask = mMaskTropic)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: error in initial dot product')
         return
      endif

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)

      do iblock=1,numBlocks

         do j=1,ny
         do i=1,nx
            S(i,j,iblock) = work1(i,j,iblock) + &
                                S(i,j,iblock)*(eta1/eta0) 
         end do
         end do

!-----------------------------------------------------------------------
!
!        compute As
!
!-----------------------------------------------------------------------

         call btropOperator(Q,S,iblock)
         do j=1,ny
         do i=1,nx
            work0(i,j,iblock) = Q(i,j,iblock)*S(i,j,iblock)
         end do
         end do

      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     compute next solution and residual
!
!-----------------------------------------------------------------------

      call POP_HaloUpdate(Q, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: error updating Q halo')
         return
      endif

      eta0 = eta1
      eta1 = eta0/POP_GlobalSum(work0, POP_distrbTropic,          &
                                POP_gridHorzLocCenter, errorCode, &
                                mMask = mMaskTropic)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: error in dot product')
         return
      endif

      !$OMP PARALLEL DO PRIVATE(iblock,i,j)

      do iblock=1,numBlocks

         do j=1,ny
         do i=1,nx
            X(i,j,iblock) = X(i,j,iblock) + eta1*S(i,j,iblock)
            R(i,j,iblock) = R(i,j,iblock) - eta1*Q(i,j,iblock)
         end do
         end do

         if (mod(m,convergenceCheckFreq) == 0) then

            call btropOperator(R,X,iblock)
            do j=1,ny
            do i=1,nx
               R(i,j,iblock) = B(i,j,iblock) - R(i,j,iblock)
               work0(i,j,iblock) = R(i,j,iblock)*R(i,j,iblock)
            end do
            end do
         endif
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!     test for convergence
!
!-----------------------------------------------------------------------

      if (mod(m,convergenceCheckFreq) == 0) then

         call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                             POP_fieldKindScalar, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversPCG: error updating residual halo in convrg')
            return
         endif

         rr = POP_GlobalSum(work0, POP_distrbTropic, &
                            POP_gridHorzLocCenter,   &
                            errorCode, mMask = mMaskTropic)   ! (r,r)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversPCG: error computing convergence dot prod')
            return
         endif

         if (rr < convergenceCriterion) then
            numIterations = m
            exit iterationLoop
         endif

      endif

   enddo iterationLoop

   rmsResidual = sqrt(rr*residualNorm)

   if (numIterations == maxIterations) then
      if (convergenceCriterion /= 0.0_POP_r8) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversPCG: solver not converged')
         return
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg

!***********************************************************************
!BOP
! !IROUTINE: ChronGear
! !INTERFACE:

 subroutine ChronGear(X,B,errorCode)

! !DESCRIPTION:
!  This routine implements the Chronopoulos-Gear conjugate-gradient 
!  solver with preconditioner for solving the linear system $Ax=b$.
!  It is a rearranged conjugate gradient solver that reduces the 
!  number of inner products per iteration from two to one (not 
!  counting convergence check). Both the operator $A$ and 
!  preconditioner are nine-point stencils. If no preconditioner has 
!  been supplied, a diagonal preconditioner is applied.  Convergence 
!  is checked every {\em ncheck} steps.
!
!
!  References:
!     Dongarra, J. and V. Eijkhout. LAPACK Working Note 159.
!        Finite-choice algorithm optimization in conjugate gradients.
!        Tech. Rep. ut-cs-03-502. Computer Science Department.
!        University of Tennessee, Knoxville. 2003.
!
!     D Azevedo, E.F., V.L. Eijkhout, and C.H. Romine. LAPACK Working
!        Note 56. Conjugate gradient algorithms with reduced
!        synchronization overhead on distributed memory multiprocessors.
!        Tech. Rep. CS-93-185.  Computer Science Department.
!        University of Tennessee, Knoxville. 1993.
!
! !REVISION HISTORY:
!  this routine implemented by Frank Bryan et al., NCAR

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(inout) :: &
      X                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,m,           &! local iteration counter
      nx, ny,          &! horizontal extents
      numBlocks,       &! number of local blocks
      iblock            ! local block counter

   real (POP_r8) :: & ! scalar results
      cgAlpha, cgBeta, cgSigma, cgDelta, cgRhoOld, cgRho, rr

   real (POP_r8), dimension(size(X,dim=1),size(X,dim=2),    &
                                          size(X,dim=3)) :: &
      R,                  &! residual (b-Ax)
      S,                  &! conjugate direction vector
      Q,Z,AZ,WORK0,       &! various cg intermediate results
      A0R                  ! diagonal preconditioner

   real (POP_r8), dimension(size(X,dim=1),size(X,dim=2),    &
                                        2,size(X,dim=3)) :: & 
      WORKN              ! WORK array 

   real (POP_r8), dimension(2) :: &
      sumN               ! global sum results for multiple arrays

!-----------------------------------------------------------------------
!
!  initialize some scalars
!
!-----------------------------------------------------------------------

   errorCode = POP_Success


   call POP_DistributionGet(POP_distrbTropic, errorCode, &
                            numLocalBlocks = numBlocks)


   nx = size(X,dim=1)
   ny = size(X,dim=2)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error retrieving local block count')
      return
   endif

   cgRho = 1.0_POP_r8
   numIterations = maxIterations


!-----------------------------------------------------------------------
!
!  compute initial residual and initialize other arrays
!
!-----------------------------------------------------------------------

   call t_startf("CG_loop1")
   !$OMP PARALLEL DO PRIVATE(iblock,i,j)

   do iblock=1,numBlocks

      do j=1,ny
      do i=1,nx
         R    (i,j,iblock)   = 0.0_POP_r8
         S    (i,j,iblock)   = 0.0_POP_r8
         Z    (i,j,iblock)   = 0.0_POP_r8
         Q    (i,j,iblock)   = 0.0_POP_r8
         AZ   (i,j,iblock)   = 0.0_POP_r8
         WORK0(i,j,iblock)   = 0.0_POP_r8
         WORKN(i,j,1,iblock) = 0.0_POP_r8
         WORKN(i,j,2,iblock) = 0.0_POP_r8
      end do
      end do

      if (.not. usePreconditioner) then
     
         !--- diagonal preconditioner if preconditioner not specified
         do j=1,ny
         do i=1,nx
            if (btropWgtCenter(i,j,iblock) /= 0.0_POP_r8) then
               A0R(i,j,iblock) = 1.0_POP_r8/btropWgtCenter(i,j,iblock)
            else
               A0R(i,j,iblock) = 0.0_POP_r8
            endif
         end do
         end do
      endif

     

      ! use S as a temp here for Ax

      call btropOperator(S,X,iblock)
      do j=1,ny
      do i=1,nx
         R(i,j,iblock) = B(i,j,iblock) - S(i,j,iblock) ! b-Ax
      end do
      end do
   end do ! block loop

   !$OMP END PARALLEL DO
   call t_stopf("CG_loop1")

   call t_startf("CG_Halo_Upd1")
   call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error updating initial residual halo')
      return
   endif
   call t_stopf("CG_Halo_Upd1")

!-----------------------------------------------------------------------
!
!    take one pass of standard algorithm
!
!-----------------------------------------------------------------------

   call t_startf("CG_loop2")
   !$OMP PARALLEL DO PRIVATE(iblock,i,j)

   do iblock=1,numBlocks

      !---- calculate (PC)r store in Z
      if (usePreconditioner) then
         call preconditioner(Z,R,iblock)
      else   ! use diagonal preconditioner
         do j=1,ny
         do i=1,nx
            Z(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
         end do
         end do
      endif


      !---- Compute intermediate result for dot product
      !---- update conjugate direction vector S
      do j=1,ny
      do i=1,nx
         WORKN(i,j,1,iblock) = R(i,j,iblock)*Z(i,j,iblock)
         S(i,j,iblock) =  Z(i,j,iblock)
      end do
      end do

      !---- compute Q = A * S
      call btropOperator(Q,S,iblock)

      !---- compute intermediate result for dot product
      do j=1,ny
      do i=1,nx
         WORKN(i,j,2,iblock) = S(i,j,iblock)*Q(i,j,iblock)
      end do
      end do

   end do
   !$OMP END PARALLEL DO
   call t_stopf("CG_loop2")

   call t_startf("CG_Halo_Upd2")
   call POP_HaloUpdate(Q, POP_haloTropic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode)


   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error updating Q halo')
      return
   endif
   call t_stopf("CG_Halo_Upd2")


   !---- Form dot products
   call t_startf("CG_globalsum1")
   sumN = POP_GlobalSum(WORKN, POP_distrbTropic,        &
                               POP_gridHorzLocCenter,   &
                               errorCode, mMask = mMaskTropic)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_SolversChronGear: error in initial dot products')
      return
   endif
   call t_stopf("CG_globalsum1")

   call t_startf("CG_loop3")
   cgRhoOld = sumN(1) !(r,PCr)
   cgSigma  = sumN(2) !(s,As)
   cgAlpha  = cgRhoOld/cgSigma

   !---- compute first solution and residual
   !$OMP PARALLEL DO PRIVATE(iblock,i,j)
   do iblock=1,numBlocks

      do j=1,ny
      do i=1,nx
         X(i,j,iblock) = X(i,j,iblock) + cgAlpha*S(i,j,iblock)
         R(i,j,iblock) = R(i,j,iblock) - cgAlpha*Q(i,j,iblock)
      end do
      end do

   end do
   !$OMP END PARALLEL DO
   call t_stopf("CG_loop3")

!-----------------------------------------------------------------------
!
!     iterate
!
!-----------------------------------------------------------------------

   iterationLoop: do m = 1, maxIterations

!-----------------------------------------------------------------------
!
!     calculate (PC)r and A*(Pc)r
!
!-----------------------------------------------------------------------

      call t_startf("CG_loop4")
      !$OMP PARALLEL DO PRIVATE(iblock,i,j)
      do iblock=1,numBlocks

         if (usePreconditioner) then
            call preconditioner(Z,R,iblock)
         else
            do j=1,ny
            do i=1,nx
               Z(i,j,iblock) = R(i,j,iblock)*A0R(i,j,iblock)
            end do
            end do
         endif

         call btropOperator(AZ,Z,iblock)

         !--- intermediate results for inner products

         do j=1,ny
         do i=1,nx
            WORKN(i,j,1,iblock) =  R(i,j,iblock)*Z(i,j,iblock)
            WORKN(i,j,2,iblock) = AZ(i,j,iblock)*Z(i,j,iblock)
         end do
         end do

      end do
      !$OMP END PARALLEL DO
      call t_stopf("CG_loop4")

      call t_startf("CG_Halo_Upd3")
      call POP_HaloUpdate(AZ, POP_haloTropic, POP_gridHorzLocCenter, &
                              POP_fieldKindScalar, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversChronGear: error updating AZ halo')
         return
      endif
      call t_stopf("CG_Halo_Upd3")

      call t_startf("CG_globalsum2")
      sumN = POP_GlobalSum(WORKN, POP_distrbTropic,        &
                                  POP_gridHorzLocCenter,   &
                                  errorCode, mMask = mMaskTropic)


      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_SolversChronGear: error in dot products')
         return
      endif
      call t_stopf("CG_globalsum2")

      call t_startf("CG_loop5")
      cgRho    = sumN(1)     ! (r,(PC)r)
      cgDelta  = sumN(2)     ! (A (PC)r,(PC)r)
      cgBeta   = cgRho/cgRhoOld
      cgSigma  = cgDelta - (cgBeta**2)*cgSigma
      cgAlpha  = cgRho/cgSigma
      cgRhoOld = cgRho

!-----------------------------------------------------------------------
!
!     compute S and Q
!     compute next solution and residual
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock, i,j)
      do iblock=1,numBlocks

         do j=1,ny
         do i=1,nx
            S(i,j,iblock) =  Z(i,j,iblock) + cgBeta *S(i,j,iblock)
            Q(i,j,iblock) = AZ(i,j,iblock) + cgBeta *Q(i,j,iblock)
            X(i,j,iblock) =  X(i,j,iblock) + cgAlpha*S(i,j,iblock)
            R(i,j,iblock) =  R(i,j,iblock) - cgAlpha*Q(i,j,iblock)
         end do
         end do

         !--- recompute residual as b-Ax for convergence check
         if (mod(m,convergenceCheckFreq) == 0) then

            !--- Reset residual using r = b - Ax
            !--- (r,r) for norm of residual
            call btropOperator(Z,X,iblock)
            do j=1,ny
            do i=1,nx
               R(i,j,iblock) = B(i,j,iblock) - Z(i,j,iblock)
               WORK0(i,j,iblock) = R(i,j,iblock)*R(i,j,iblock)
            end do
            end do

         endif
      end do
      !$OMP END PARALLEL DO
      call t_stopf("CG_loop5")

!-----------------------------------------------------------------------
!
!     test for convergence if it is time
!
!-----------------------------------------------------------------------

      if (mod(m,convergenceCheckFreq) == 0) then

         call t_startf("CG_Halo_Upd4")
         !--- update ghost cells for next iteration
         call POP_HaloUpdate(R, POP_haloTropic, POP_gridHorzLocCenter, &
                                POP_fieldKindScalar, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversChronGear: error updating residual halo in convrg')
            return
         endif
         call t_stopf("CG_Halo_Upd4")

         !--- residual norm for convergence 
         call t_startf("CG_globalsum3")
         rr = POP_GlobalSum(work0, POP_distrbTropic,        &! (r,r)
                                   POP_gridHorzLocCenter,   &
                                   errorCode, mMask = mMaskTropic)


         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_SolversChronGear: error computing convergence dot prod')
            return
         endif
         call t_stopf("CG_globalsum3")

         if (rr < convergenceCriterion) then
            numIterations = m
            exit iterationLoop
         endif

      endif

   end do iterationLoop

   rmsResidual = sqrt(rr*residualNorm)

   if (numIterations == maxIterations) then
      if (convergenceCriterion /= 0.0_POP_r8) then

         call POP_ErrorSet(errorCode, &
            'POP_SolversChronGear: solver not converged')
         return
      endif

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ChronGear

!***********************************************************************
!BOP
! !IROUTINE: preconditioner
! !INTERFACE:

 subroutine preconditioner(PX,X,bid)

! !DESCRIPTION:
!  This function applies a precomputed preconditioner as a nine-point
!  stencil operator.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: & 
      X                     ! array to be operated on 

   integer (POP_i4), intent(in) :: &
      bid                    ! local block address for this block

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(out) :: &
      PX                  ! nine point operator result

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j                ! dummy counters

!-----------------------------------------------------------------------

   PX(:,:,bid) = 0.0_POP_r8

   do j=2,size(X,dim=2)-1
   do i=2,size(X,dim=1)-1
      PX(i,j,bid) = precondNE    (i,j,bid)*X(i+1,j+1,bid) + &
                    precondNW    (i,j,bid)*X(i-1,j+1,bid) + &
                    precondSE    (i,j,bid)*X(i+1,j-1,bid) + &
                    precondSW    (i,j,bid)*X(i-1,j-1,bid) + &
                    precondNorth (i,j,bid)*X(i  ,j+1,bid) + &
                    precondSouth (i,j,bid)*X(i  ,j-1,bid) + &
                    precondEast  (i,j,bid)*X(i+1,j  ,bid) + &
                    precondWest  (i,j,bid)*X(i-1,j  ,bid) + &
                    precondCenter(i,j,bid)*X(i  ,j  ,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine preconditioner

!***********************************************************************
!BOP
! !IROUTINE: btropOperator
! !INTERFACE:

 subroutine btropOperator(AX,X,bid)

! !DESCRIPTION:
!  This routine applies the nine-point stencil operator for the
!  barotropic solver.  It takes advantage of some 9pt weights being 
!  shifted versions of others.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(in) :: & 
      X                  ! array to be operated on 

   integer (POP_i4), intent(in) :: &
      bid                    ! local block address for this block

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(:,:,:), intent(out) :: &
      AX                     ! nine point operator result (Ax)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j                ! dummy counters

!-----------------------------------------------------------------------

   AX(:,:,bid) = 0.0_POP_r8

   do j=2,size(X,dim=2)-1
   do i=2,size(X,dim=1)-1
      AX(i,j,bid) = btropWgtCenter(i  ,j  ,bid)*X(i  ,j  ,bid) + &
                    btropWgtNorth (i  ,j  ,bid)*X(i  ,j+1,bid) + &
                    btropWgtNorth (i  ,j-1,bid)*X(i  ,j-1,bid) + &
                    btropWgtEast  (i  ,j  ,bid)*X(i+1,j  ,bid) + &
                    btropWgtEast  (i-1,j  ,bid)*X(i-1,j  ,bid) + &
                    btropWgtNE    (i  ,j  ,bid)*X(i+1,j+1,bid) + &
                    btropWgtNE    (i  ,j-1,bid)*X(i+1,j-1,bid) + &
                    btropWgtNE    (i-1,j  ,bid)*X(i-1,j+1,bid) + &
                    btropWgtNE    (i-1,j-1,bid)*X(i-1,j-1,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine btropOperator

!***********************************************************************

 end module POP_SolversMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
