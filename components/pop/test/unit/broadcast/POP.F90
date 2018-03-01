!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 program POP_BroadcastTest

!----------------------------------------------------------------------
!
!  this program tests POP broadcast communications
!
!----------------------------------------------------------------------

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_BroadcastMod

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4), parameter :: &
      nx = 3, ny = 3, nz = 5, nl = 2, nm = 2, nn = 2

   integer (POP_i4) :: &
      errorCode,       &
      i,j,k,l,m,n,     &! loop indices
      i4Test, i4Expect  ! scalar test values

   real (POP_r4) :: &
      r4Test, r4Expect  ! scalar test values

   real (POP_r8) :: &
      r8Test, r8Expect  ! scalar test values

   character (POP_charLength) :: &
      charTest, charExpect  ! scalar test values

   logical (POP_logical) :: &
      logTest, logExpect  ! scalar test values

   integer (POP_i4), dimension(nx) :: &
      i4Test1D, i4Expect1D  ! 1D test values

   real (POP_r4), dimension(nx) :: &
      r4Test1D, r4Expect1D  ! 1D test values

   real (POP_r8), dimension(nx) :: &
      r8Test1D, r8Expect1D  ! 1D test values

   character (POP_charLength), dimension(nx) :: &
      charTest1D, charExpect1D  ! 1D test values

   logical (POP_logical), dimension(nx) :: &
      logTest1D, logExpect1D  ! 1D test values

   integer (POP_i4), dimension(nx,ny) :: &
      i4Test2D, i4Expect2D  ! 2D test values

   real (POP_r4), dimension(nx,ny) :: &
      r4Test2D, r4Expect2D  ! 2D test values

   real (POP_r8), dimension(nx,ny) :: &
      r8Test2D, r8Expect2D  ! 2D test values

   logical (POP_logical), dimension(nx,ny) :: &
      logTest2D, logExpect2D  ! 2D test values

   integer (POP_i4), dimension(nx,ny,nz) :: &
      i4Test3D, i4Expect3D  ! 3D test values

   real (POP_r4), dimension(nx,ny,nz) :: &
      r4Test3D, r4Expect3D  ! 3D test values

   real (POP_r8), dimension(nx,ny,nz) :: &
      r8Test3D, r8Expect3D  ! 3D test values

   logical (POP_logical), dimension(nx,ny,nz) :: &
      logTest3D, logExpect3D  ! 3D test values

   integer (POP_i4), dimension(nx,ny,nz,nl) :: &
      i4Test4D, i4Expect4D  ! 4D test values

   real (POP_r4), dimension(nx,ny,nz,nl) :: &
      r4Test4D, r4Expect4D  ! 4D test values

   real (POP_r8), dimension(nx,ny,nz,nl) :: &
      r8Test4D, r8Expect4D  ! 4D test values

   logical (POP_logical), dimension(nx,ny,nz,nl) :: &
      logTest4D, logExpect4D  ! 4D test values

   integer (POP_i4), dimension(nx,ny,nz,nl,nm) :: &
      i4Test5D, i4Expect5D  ! 5D test values

   real (POP_r4), dimension(nx,ny,nz,nl,nm) :: &
      r4Test5D, r4Expect5D  ! 5D test values

   real (POP_r8), dimension(nx,ny,nz,nl,nm) :: &
      r8Test5D, r8Expect5D  ! 5D test values

   logical (POP_logical), dimension(nx,ny,nz,nl,nm) :: &
      logTest5D, logExpect5D  ! 5D test values

   integer (POP_i4), dimension(nx,ny,nz,nl,nm,nn) :: &
      i4Test6D, i4Expect6D  ! 6D test values

   real (POP_r4), dimension(nx,ny,nz,nl,nm,nn) :: &
      r4Test6D, r4Expect6D  ! 6D test values

   real (POP_r8), dimension(nx,ny,nz,nl,nm,nn) :: &
      r8Test6D, r8Expect6D  ! 6D test values

   logical (POP_logical), dimension(nx,ny,nz,nl,nm,nn) :: &
      logTest6D, logExpect6D  ! 6D test values

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
!  initialize test arrays
!
!----------------------------------------------------------------------

   i4Expect   =   10
   r4Expect   =  100.0_POP_r4
   r8Expect   = 1000.0_POP_r8
   charExpect = 'defined'
   logExpect  = .true.

   do i=1,nx
      i4Expect1D(i) = i*i4Expect
      r4Expect1D(i) = i*r4Expect
      r8Expect1D(i) = i*r8Expect
      charExpect1D(i) = charExpect
      logExpect1D(i) = logExpect
   end do

   do j=1,ny
   do i=1,nx
      i4Expect2D(i,j) = (i+j)*i4Expect
      r4Expect2D(i,j) = (i+j)*r4Expect
      r8Expect2D(i,j) = (i+j)*r8Expect
      if (mod(i+j,2) == 0) then
         logExpect2D(i,j) = logExpect
      else
         logExpect2D(i,j) = .false.
      endif
   end do
   end do

   do k=1,nz
   do j=1,ny
   do i=1,nx
      i4Expect3D(i,j,k) = (i+j+k)*i4Expect
      r4Expect3D(i,j,k) = (i+j+k)*r4Expect
      r8Expect3D(i,j,k) = (i+j+k)*r8Expect
      if (mod(i+j+k,2) == 0) then
         logExpect3D(i,j,k) = logExpect
      else
         logExpect3D(i,j,k) = .false.
      endif
   end do
   end do
   end do

   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      i4Expect4D(i,j,k,l) = (i+j+k+l)*i4Expect
      r4Expect4D(i,j,k,l) = (i+j+k+l)*r4Expect
      r8Expect4D(i,j,k,l) = (i+j+k+l)*r8Expect
      if (mod(i+j+k+l,2) == 0) then
         logExpect4D(i,j,k,l) = logExpect
      else
         logExpect4D(i,j,k,l) = .false.
      endif
   end do
   end do
   end do
   end do

   do m=1,nm
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      i4Expect5D(i,j,k,l,m) = (i+j+k+l+m)*i4Expect
      r4Expect5D(i,j,k,l,m) = (i+j+k+l+m)*r4Expect
      r8Expect5D(i,j,k,l,m) = (i+j+k+l+m)*r8Expect
      if (mod(i+j+k+l+m,2) == 0) then
         logExpect5D(i,j,k,l,m) = logExpect
      else
         logExpect5D(i,j,k,l,m) = .false.
      endif
   end do
   end do
   end do
   end do
   end do

   do n=1,nn
   do m=1,nm
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      i4Expect6D(i,j,k,l,m,n) = (i+j+k+l+m+n)*i4Expect
      r4Expect6D(i,j,k,l,m,n) = (i+j+k+l+m+n)*r4Expect
      r8Expect6D(i,j,k,l,m,n) = (i+j+k+l+m+n)*r8Expect
      if (mod(i+j+k+l+m+n,2) == 0) then
         logExpect6D(i,j,k,l,m,n) = logExpect
      else
         logExpect6D(i,j,k,l,m,n) = .false.
      endif
   end do
   end do
   end do
   end do
   end do
   end do

   if (POP_myTask == POP_masterTask) then
      i4Test    = i4Expect
      r4Test    = r4Expect
      r8Test    = r8Expect
      charTest  = charExpect
      logTest   = logExpect
      i4Test1D  = i4Expect1D
      r4Test1D  = r4Expect1D
      r8Test1D  = r8Expect1D
      charTest1D = charExpect1D
      logTest1D = logExpect1D
      i4Test2D  = i4Expect2D
      r4Test2D  = r4Expect2D
      r8Test2D  = r8Expect2D
      logTest2D = logExpect2D
      i4Test3D  = i4Expect3D
      r4Test3D  = r4Expect3D
      r8Test3D  = r8Expect3D
      logTest3D = logExpect3D
      i4Test4D  = i4Expect4D
      r4Test4D  = r4Expect4D
      r8Test4D  = r8Expect4D
      logTest4D = logExpect4D
      i4Test5D  = i4Expect5D
      r4Test5D  = r4Expect5D
      r8Test5D  = r8Expect5D
      logTest5D = logExpect5D
      i4Test6D  = i4Expect6D
      r4Test6D  = r4Expect6D
      r8Test6D  = r8Expect6D
      logTest6D = logExpect6D
   else
      i4Test    = 0
      r4Test    = 0.0_POP_r4
      r8Test    = 0.0_POP_r8
      charTest  = ' '
      logTest   = .false.
      i4Test1D  = 0
      r4Test1D  = 0.0_POP_r4
      r8Test1D  = 0.0_POP_r8
      charTest1D = ' '
      logTest1D = .false.
      i4Test2D  = 0
      r4Test2D  = 0.0_POP_r4
      r8Test2D  = 0.0_POP_r8
      logTest2D = .false.
      i4Test3D  = 0
      r4Test3D  = 0.0_POP_r4
      r8Test3D  = 0.0_POP_r8
      logTest3D = .false.
      i4Test4D  = 0
      r4Test4D  = 0.0_POP_r4
      r8Test4D  = 0.0_POP_r8
      logTest4D = .false.
      i4Test5D  = 0
      r4Test5D  = 0.0_POP_r4
      r8Test5D  = 0.0_POP_r8
      logTest5D = .false.
      i4Test6D  = 0
      r4Test6D  = 0.0_POP_r4
      r8Test6D  = 0.0_POP_r8
      logTest6D = .false.
   endif

!----------------------------------------------------------------------
!
!  test broadcast of each size, type
!
!----------------------------------------------------------------------

   call POP_Broadcast(  r8Test, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r8 scalar broadcast')
   endif
   call POP_Broadcast(  r4Test, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r4 scalar broadcast')
   endif
   call POP_Broadcast(  i4Test, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in i4 scalar broadcast')
   endif
   call POP_Broadcast( logTest, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in log scalar broadcast')
   endif
   call POP_Broadcast(charTest, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in char scalar broadcast')
   endif
   call POP_Broadcast(  r8Test1D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r8 1D broadcast')
   endif
   call POP_Broadcast(  r4Test1D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r4 1D broadcast')
   endif
   call POP_Broadcast(  i4Test1D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in i4 1D broadcast')
   endif
   call POP_Broadcast( logTest1D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in log 1D broadcast')
   endif
   call POP_Broadcast(charTest1D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in char 1D broadcast')
   endif
   call POP_Broadcast( r8Test2D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r8 2D broadcast')
   endif
   call POP_Broadcast( r4Test2D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r4 2D broadcast')
   endif
   call POP_Broadcast( i4Test2D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in i4 2D broadcast')
   endif
   call POP_Broadcast(logTest2D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in log 2D broadcast')
   endif
   call POP_Broadcast( r8Test3D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r8 3D broadcast')
   endif
   call POP_Broadcast( r4Test3D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r4 3D broadcast')
   endif
   call POP_Broadcast( i4Test3D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in i4 3D broadcast')
   endif
   call POP_Broadcast(logTest3D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in log 3D broadcast')
   endif
   call POP_Broadcast( r8Test4D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r8 4D broadcast')
   endif
   call POP_Broadcast( r4Test4D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r4 4D broadcast')
   endif
   call POP_Broadcast( i4Test4D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in i4 4D broadcast')
   endif
   call POP_Broadcast(logTest4D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in log 4D broadcast')
   endif
   call POP_Broadcast( r8Test5D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r8 5D broadcast')
   endif
   call POP_Broadcast( r4Test5D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r4 5D broadcast')
   endif
   call POP_Broadcast( i4Test5D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in i4 5D broadcast')
   endif
   call POP_Broadcast(logTest5D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in log 5D broadcast')
   endif
   call POP_Broadcast( r8Test6D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r8 6D broadcast')
   endif
   call POP_Broadcast( r4Test6D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in r4 6D broadcast')
   endif
   call POP_Broadcast( i4Test6D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in i4 6D broadcast')
   endif
   call POP_Broadcast(logTest6D, POP_masterTask, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: error in log 6D broadcast')
   endif

!----------------------------------------------------------------------
!
!  check for proper values
!
!----------------------------------------------------------------------

   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r8 scalar broadcast')
   endif

   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r4 scalar broadcast')
   endif

   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from i4 scalar broadcast')
   endif

   if (logTest /= logExpect) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from log scalar broadcast')
   endif

   if (charTest /= charExpect) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from char scalar broadcast')
   endif

   if (count(r8Test1D /= r8Expect1D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r8 1D broadcast')
   endif

   if (count(r4Test1D /= r4Expect1D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r4 1D broadcast')
   endif

   if (count(i4Test1D /= i4Expect1D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from i4 1D broadcast')
   endif

   if (count(logTest1D /= logExpect1D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from log 1D broadcast')
   endif

   if (count(charTest1D /= charExpect1D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from char 1D broadcast')
   endif

   if (count(r8Test2D /= r8Expect2D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r8 2D broadcast')
   endif

   if (count(r4Test2D /= r4Expect2D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r4 2D broadcast')
   endif

   if (count(i4Test2D /= i4Expect2D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from i4 2D broadcast')
   endif

   if (count(logTest2D /= logExpect2D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from log 2D broadcast')
   endif

   if (count(r8Test3D /= r8Expect3D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r8 3D broadcast')
   endif

   if (count(r4Test3D /= r4Expect3D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r4 3D broadcast')
   endif

   if (count(i4Test3D /= i4Expect3D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from i4 3D broadcast')
   endif

   if (count(logTest3D /= logExpect3D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from log 3D broadcast')
   endif

   if (count(r8Test4D /= r8Expect4D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r8 4D broadcast')
   endif

   if (count(r4Test4D /= r4Expect4D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r4 4D broadcast')
   endif

   if (count(i4Test4D /= i4Expect4D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from i4 4D broadcast')
   endif

   if (count(logTest4D /= logExpect4D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from log 4D broadcast')
   endif

   if (count(r8Test5D /= r8Expect5D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r8 5D broadcast')
   endif

   if (count(r4Test5D /= r4Expect5D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r4 5D broadcast')
   endif

   if (count(i4Test5D /= i4Expect5D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from i4 5D broadcast')
   endif

   if (count(logTest5D /= logExpect5D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from log 5D broadcast')
   endif

   if (count(r8Test6D /= r8Expect6D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r8 6D broadcast')
   endif

   if (count(r4Test6D /= r4Expect6D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from r4 6D broadcast')
   endif

   if (count(i4Test6D /= i4Expect6D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from i4 6D broadcast')
   endif

   if (count(logTest6D /= logExpect6D) > 0) then
      call POP_ErrorSet(errorCode, &
         'BcastTest: bad values from log 6D broadcast')
   endif

!----------------------------------------------------------------------
!
!  clean up
!
!----------------------------------------------------------------------

   call POP_ErrorPrint(errorCode)
   call POP_CommExitMessageEnvironment

!----------------------------------------------------------------------

 end program POP_BroadcastTest

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
