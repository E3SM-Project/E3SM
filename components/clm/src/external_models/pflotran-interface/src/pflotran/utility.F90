module Utility_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private


  interface DotProduct
    module procedure DotProduct1
    module procedure DotProduct2
    module procedure DotProduct3
  end interface
  
  interface CrossProduct
    module procedure CrossProduct1
  end interface
  
  interface reallocateRealArray
    module procedure reallocateRealArray1D
    module procedure reallocateRealArray2D
  end interface
  
  interface UtilityReadArray
    module procedure UtilityReadIntArray
    module procedure UtilityReadRealArray
  end interface

  interface DeallocateArray
    ! TODO(geh) replace deallocations with the below
    module procedure DeallocateArray1DInteger
    module procedure DeallocateArray2DInteger
    module procedure DeallocateArray3DInteger
    module procedure DeallocateArray1DReal
    module procedure DeallocateArray2DReal
    module procedure DeallocateArray3DReal
    module procedure DeallocateArray1DLogical
    module procedure DeallocateArray2DLogical
    module procedure DeallocateArray3DLogical
    module procedure DeallocateArray1DString
    module procedure DeallocateArray2DString
  end interface
  
  interface InterfaceApprox
    module procedure InterfaceApproxWithDeriv
    module procedure InterfaceApproxWithoutDeriv
  end interface
  
  interface CalcParallelSUM
    module procedure CalcParallelSUM1
    module procedure CalcParallelSUM2
  end interface

  public :: GetRndNumFromNormalDist, &
            DotProduct, &
            CrossProduct, &
            reallocateRealArray, &
            reallocateIntArray, &
            UtilityReadArray, &
            DeallocateArray, &
            InterfaceApprox, &
            Interpolate, &
            GradientLinear, &
            InterpolateBilinear, &
            SearchOrderedArray, &
            ludcmp, &
            lubksb, &
            FileExists, &
            Equal, &
            BestFloat, &
            QuadraticPolynomialSetup, &
            QuadraticPolynomialEvaluate, &
            CubicPolynomialSetup, &
            CubicPolynomialEvaluate, &
            ConvertMatrixToVector, &
            Kron, &
            Transposer, &
            Determinant, &
            InterfaceApproxWithDeriv, &
            InterfaceApproxWithoutDeriv, &
            PrintProgressBarInt, &
            InverseNorm, &
            Erf_, &
            DigitsOfAccuracy, &
            CalcParallelSum, &
            MatCompare

  public :: HFunctionSmooth, &             ! F.-M. Yuan (2017-03-09)
            where_checkerr                 ! F.-M. Yuan (2018-04-11)
contains

! ************************************************************************** !

function rnd()

  implicit none

  integer*8, save :: iseed = 1
  PetscReal :: rnd

  iseed = iseed*125
  iseed = iseed - (iseed/2796203) * 2796203
  rnd   = iseed/2796203.0
  return
end function rnd

! ************************************************************************** !

function ran1(idum)
      
  implicit none
      
  save

!-----returns a random number in the range (0,1). Set idum to neg.
!     value to initialize

  PetscReal :: ran1
  PetscReal :: r(97),rm1,rm2
  PetscInt :: idum,iff,ix1,ix2,ix3,j,m1,ia1,ic1,m2,ia2,ic2,m3,ia3,ic3

  parameter (M1  = 259200)
  parameter (IA1 = 7141)
  parameter (IC1 = 54773)
  parameter (RM1 = 1.0/M1)
  parameter (M2  = 134456)
  parameter (IA2 = 8121)
  parameter (IC2 = 28411)
  parameter (RM2 = 1.0/M2)
  parameter (M3  = 243000)
  parameter (IA3 = 4561)
  parameter (IC3 = 51349)

  data iff/0/

  if (idum.lt.0 .or. iff.eq.0) then
    iff=1
    ix1=mod(IC1-idum,M1)
    ix1=mod(IA1*ix1+IC1,M1)
    ix2=mod(ix1,M2)
    ix1=mod(IA1*ix1+IC1,M1)
    ix3=mod(ix1,M3)
    do j=1,97
      ix1=mod(IA1*ix1+IC1,M1)
      ix2=mod(IA2*ix2+IC2,M2)
      r(j)=(float(ix1)+float(ix2)*RM2)*RM1
    enddo
    idum=1
  endif
  ix1=mod(IA1*ix1+IC1,M1)
  ix2=mod(IA2*ix2+IC2,M2)
  ix3=mod(IA3*ix3+IC3,M3)
  j=1+(97*ix3)/M3
!  if (j.gt.97 .or. j.lt.1) pause

  ran1=r(j)
  r(j)=(float(ix1)+float(ix2)*RM2)*RM1
      
  return
end function ran1

! ************************************************************************** !

function ran2(idum)
      
  implicit none

!-----Minimal random number generator of Park and Miller 
!     in the range (0,1)

  PetscReal :: ran2, AM, EPS, RNMX, temp
  PetscInt :: IA, IM, IQ, IR, NTAB, idum, iy, j, k, iv(32), NDIV

  parameter (IA = 16807)
  parameter (IM = 2147483647)
  parameter (AM = 1.0/IM)
  parameter (IQ = 127773)
  parameter (IR = 2836)
  parameter (NTAB = 32)
  parameter (NDIV = 1+(IM-1)/NTAB)
  parameter (EPS  = 1.2e-7)
  parameter (RNMX = 1.0-EPS)

  !dimension iv(NTAB)

  iy = 0
  if (idum.le.0 .or. iy.eq.0) then
    if (-idum .lt. 1) then
      idum = 1
    else 
      idum = -idum
    endif
    do j = NTAB+7,0,-1
      k = idum/IQ
      idum = IA*(idum-k*IQ)-IR*k
      if (idum .lt. 0) idum = idum+IM
      if (j .lt. NTAB) iv(j) = idum
    enddo
    iy = iv(1)
  endif
  k = idum/IQ
  idum = IA*(idum-k*IQ)-IR*k
  if (idum .lt. 0) idum = idum+IM
  j= iy/NDIV
  iy = iv(j)
  iv(j) = idum
  temp = AM*iy
  if (temp .gt. RNMX) then
    ran2 = RNMX
  else 
    ran2 = temp
  endif

  return
end function ran2

! ************************************************************************** !

subroutine GetRndNumFromNormalDist(mean,st_dev,number)
  ! 
  ! Generates a random number that is normally distributed, as defined by the
  ! mean and standard deviation given. This subroutine uses the Box-Muller
  ! transform.
  ! G. E. P. Box and M. E. Muller (1958), A note on the generation of random
  ! normal deviates, The Annals of Mathematical Statistics, Vol. 29, No. 2,
  ! pp. 610-611.
  ! 
  ! Author: Jenn Frederick
  ! Date: 2/12/2016

  implicit none
  
  PetscReal :: mean, st_dev, number

  PetscBool, save :: switch
  PetscReal, save :: z0, z1
  PetscReal :: u1, u2
  PetscReal :: TWO_PI

  switch = .not.switch
  TWO_PI = 2*3.14159265358979323846264338327950288419716939937510582

  if (.not.switch) then

    ! Generate two random numbers between (0,1)
    u1 = rnd()
    u2 = rnd()
    
    z0 = sqrt(-2.0*log(u1)) * cos(TWO_PI*u2)
    z1 = sqrt(-2.0*log(u1)) * sin(TWO_PI*u2) 
    
    number = z0*st_dev + mean

  else
    
    number = z1*st_dev + mean

  endif

end subroutine GetRndNumFromNormalDist
! ************************************************************************** !

subroutine Natural2LocalIndex(ir, nl, llist, llength)
  implicit none
  PetscInt :: nl, ir,na, l_search, itt, llength
  PetscInt :: llist(*)
  
  PetscInt ::  nori0, nori1, nori
  
  
  nl=-1
  l_search = llength
  
  na = ir!-1 
  itt=0
  nori0 =1
  nori1 = llength
  if (na>=llist(1) .and. na <= llist(llength))then
  do while(l_search > 1 .and.itt<=50)
  
    itt=itt+1
    if (na == llist(nori0))then
      nl = nori0
      exit
    elseif (na == llist(nori1))then
       nl = nori1
      exit
    endif   
     
     ! nori = int((real(nori0 + nori1))/ 2.) + mod ( nori0 + nori1,2 )
    nori =  int(floor(real(nori0+nori1)/2D0 + .75D0))
    if ( na > llist(nori)) then
      nori0 = nori
    elseif (na < llist(nori))then
      nori1 = nori
    else
      if (na == llist(nori))then
        nl = nori
        exit
      else
        print *, 'wrong index', na, nori, llist(nori); stop
      endif  
    endif
    l_search = nori1-nori0
    if (itt>=40)then
      print *, na, nori0,nori1,nori, llist(nori0), llist(nori1)
      if (itt>=50) stop
    endif
  enddo  
 endif         
          
end subroutine Natural2LocalIndex

! ************************************************************************** !

subroutine reallocateIntArray(array,size)
  ! 
  ! Reallocates an integer array to a larger size and copies
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  ! 

  implicit none

  PetscInt, pointer :: array(:)
  PetscInt :: size
  
  PetscInt, allocatable :: array2(:)
  
  allocate(array2(size))
  array2(1:size) = array(1:size)
  deallocate(array)
  allocate(array(2*size))
  array = 0
  array(1:size) = array2(1:size)
  size = 2*size
  deallocate(array2)

end subroutine reallocateIntArray

! ************************************************************************** !

subroutine reallocateRealArray1D(array,size)
  ! 
  ! reallocateRealArray2D: Reallocates a 2D real array to a larger size and
  ! copies values over.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/29/07, 10/03/13
  ! 

  implicit none

  PetscReal, pointer :: array(:)
  PetscInt :: size
  
  PetscReal, allocatable :: array2(:)
  
  allocate(array2(size))
  array2(1:size) = array(1:size)
  deallocate(array)
  allocate(array(2*size))
  array = 0.d0
  array(1:size) = array2(1:size)
  size = 2*size
  deallocate(array2)

end subroutine reallocateRealArray1D

! ************************************************************************** !

subroutine reallocateRealArray2D(array,rank2_size)
  ! 
  ! Reallocates a 2D real array to a larger size in last
  ! dimension and copies values over.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/03/13
  ! 

  implicit none

  PetscReal, pointer :: array(:,:)
  PetscInt :: rank1_size, rank2_size
  
  PetscReal, allocatable :: array2(:,:)
  
  rank1_size = size(array,1)
  allocate(array2(rank1_size,rank2_size))
  array2(:,1:rank2_size) = array(:,1:rank2_size)
  deallocate(array)
  allocate(array(rank1_size,2*rank2_size))
  array = 0.d0
  array(:,1:rank2_size) = array2(:,1:rank2_size)
  rank2_size = 2*rank2_size
  deallocate(array2)

end subroutine reallocateRealArray2D

! ************************************************************************** !

subroutine ludcmp(A,N,INDX,D)
  ! 
  ! Given an NxN matrix A, with physical dimension NP, this routine replaces it
  ! by the LU decomposition of a rowwise permutation of itself.
  ! A and N are input. A is output; INDX is output vector which records the
  ! row permutation effected by the partial pivoting; D id output as +1 or -1
  ! depending on whether the number of row interchanges was odd or even,
  ! respectively. This routine is used in combination with lubksb to solve
  ! linear equations or invert a matrix.
  ! 

  implicit none

  PetscInt :: N
  PetscReal, parameter :: tiny=1.0d-20
  PetscReal :: A(N,N),VV(N)
  PetscInt :: INDX(N)
  PetscInt :: D

  PetscInt :: i, j, k, imax
  PetscReal :: aamax, sum, dum
  PetscMPIInt ::  rank
  PetscErrorCode :: ierr

  D=1
  do i=1,N
    aamax=0.d0
    do j=1,N
      if (abs(A(i,j)).gt.aamax) aamax=abs(A(i,j))
    enddo
    if (aamax <= 0.d0) then
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      print *, "ERROR: Singular value encountered in ludcmp() on processor: ", rank, ' aamax = ',aamax,' row = ',i
      do k = 1, N
        print *, "Jacobian: ",k,(j,A(k,j),j=1,N)
      enddo
      call MPI_Abort(MPI_COMM_WORLD,ONE_INTEGER_MPI,ierr)
      call MPI_Finalize(ierr)
      stop
    endif
    VV(i)=1./aamax
  enddo
  do j=1,N
    do i=1,j-1
      sum=A(i,j)
      do k=1,i-1
        sum=sum-A(i,k)*A(k,j)
      enddo
      A(i,j)=sum
    enddo
    aamax=0
    do i=j,N
      sum=A(i,j)
      do k=1,j-1
        sum=sum-A(i,k)*A(k,j)
      enddo
      A(i,j)=sum
      dum=VV(i)*abs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
    enddo
    if (j.ne.imax) then
      do k=1,N
        dum=A(imax,k)
        A(imax,k)=A(j,k)
        A(j,k)=dum
      enddo
      D=-D
      VV(imax)=VV(j)
    endif
    INDX(j)=imax
    if (A(j,j).eq.0.d0) A(j,j)=tiny
    if (j.ne.N) then
      dum=1.d0/A(j,j)
      do i=j+1,N
        A(i,j)=A(i,j)*dum
      enddo
    endif
  enddo
  return

end subroutine ludcmp

! ************************************************************************** !

subroutine lubksb(A,N,INDX,B)
  ! 
  ! Solves the set of N linear equations A.X=D. Here A is input, not as a matrix
  ! A but rather as its LU decomposition. INDX is the input as the permutation
  ! vector returned by ludcmp. B is input as the right-hand side vector B, and
  ! returns with the solution vector X.
  ! 

  implicit none

  PetscInt :: N
  PetscReal :: A(N,N),B(N)
  PetscInt :: INDX(N)

  PetscInt :: i, j, ii, ll
  PetscReal :: sum


  ii=0
  do i=1,N
    ll=INDX(i)
    sum=B(ll)
    B(ll)=B(i)
    if (ii.ne.0) then
      do j=ii,i-1
        sum=sum-A(i,j)*B(j)
      enddo
    else if (sum.ne.0.d0) then
      ii=i
    endif
    B(i)=sum
  enddo
  do i=N,1,-1
    sum=B(i)
    if (i.lt.N) then
      do j=i+1,N
        sum=sum-A(i,j)*B(j)
      enddo
    endif
    B(i)=sum/A(i,i)
  enddo
  return

end subroutine lubksb

! ************************************************************************** !

subroutine ludcmp_chunk(A,N,INDX,D,chunk_size,ithread,num_threads)
  ! 
  ! Given an NxN matrix A, with physical dimension NP, this routine replaces it
  ! by the LU decomposition of a rowwise permutation of itself.
  ! A and N are input. A is output; INDX is output vector which records the
  ! row permutation effected by the partial pivoting; D id output as +1 or -1
  ! depending on whether the number of row interchanges was odd or even,
  ! respectively. This routine is used in combination with lubksb to solve
  ! linear equations or invert a matrix.
  ! 

  implicit none

  PetscInt :: N
  PetscInt :: chunk_size
  PetscInt :: num_threads
  PetscReal, parameter :: tiny=1.0d-20
  PetscReal :: A(chunk_size,num_threads,N,N),VV(chunk_size,num_threads,N)
  PetscInt :: INDX(chunk_size,num_threads,N)
  PetscInt :: D(chunk_size,num_threads)
  PetscInt :: ithread

  PetscInt :: i, j, k, imax
  PetscReal :: aamax, sum, dum
  PetscMPIInt ::  rank
  PetscErrorCode :: ierr
  
  PetscInt :: ichunk

  do ichunk = 1, chunk_size

  D(ichunk,ithread)=1
  do i=1,N
    aamax=0
    do j=1,N
      if (abs(A(ichunk,ithread,i,j)).gt.aamax) aamax=abs(A(ichunk,ithread,i,j))
    enddo
    if (aamax.eq.0.d0) then
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      print *, "ERROR: Singular value encountered in ludcmp() on processor", rank, ichunk,ithread
      call MPI_Abort(MPI_COMM_WORLD,ONE_INTEGER_MPI,ierr)
      call MPI_Finalize(ierr)
      stop
    endif
    VV(ichunk,ithread,i)=1./aamax
  enddo
  do j=1,N
    do i=1,j-1
      sum=A(ichunk,ithread,i,j)
      do k=1,i-1
        sum=sum-A(ichunk,ithread,i,k)*A(ichunk,ithread,k,j)
      enddo
      A(ichunk,ithread,i,j)=sum
    enddo
    aamax=0
    do i=j,N
      sum=A(ichunk,ithread,i,j)
      do k=1,j-1
        sum=sum-A(ichunk,ithread,i,k)*A(ichunk,ithread,k,j)
      enddo
      A(ichunk,ithread,i,j)=sum
      dum=VV(ichunk,ithread,i)*abs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
    enddo
    if (j.ne.imax) then
      do k=1,N
        dum=A(ichunk,ithread,imax,k)
        A(ichunk,ithread,imax,k)=A(ichunk,ithread,j,k)
        A(ichunk,ithread,j,k)=dum
      enddo
      D(ichunk,ithread)=-D(ichunk,ithread)
      VV(ichunk,ithread,imax)=VV(ichunk,ithread,j)
    endif
    INDX(ichunk,ithread,j)=imax
    if (A(ichunk,ithread,j,j).eq.0.d0) A(ichunk,ithread,j,j)=tiny
    if (j.ne.N) then
      dum=1./A(ichunk,ithread,j,j)
      do i=j+1,N
        A(ichunk,ithread,i,j)=A(ichunk,ithread,i,j)*dum
      enddo
    endif
  enddo
  
  enddo ! chunk loop
  
  return

end subroutine ludcmp_chunk

! ************************************************************************** !

subroutine lubksb_chunk(A,N,INDX,B,chunk_size,ithread,num_threads)
  ! 
  ! Solves the set of N linear equations A.X=D. Here A is input, not as a matrix
  ! A but rather as its LU decomposition. INDX is the input as the permutation
  ! vector returned bu ludcmp. B is input as the right-hand side vector B, and
  ! returns with the solution vector X.
  ! 

  implicit none

  PetscInt :: N
  PetscInt :: chunk_size
  PetscInt :: num_threads
  PetscReal :: A(chunk_size,num_threads,N,N),B(chunk_size,num_threads,N)
  PetscInt :: INDX(chunk_size,num_threads,N)
  PetscInt :: ithread

  PetscInt :: i, j, ii, ll
  PetscReal :: sum

  PetscInt :: ichunk

  do ichunk = 1, chunk_size
  
  ii=0
  do i=1,N
    ll=INDX(ichunk,ithread,i)
    sum=B(ichunk,ithread,ll)
    B(ichunk,ithread,ll)=B(ichunk,ithread,i)
    if (ii.ne.0) then
      do j=ii,i-1
        sum=sum-A(ichunk,ithread,i,j)*B(ichunk,ithread,j)
      enddo
    else if (sum.ne.0.d0) then
      ii=i
    endif
    B(ichunk,ithread,i)=sum
  enddo
  do i=N,1,-1
    sum=B(ichunk,ithread,i)
    if (i.lt.N) then
      do j=i+1,N
        sum=sum-A(ichunk,ithread,i,j)*B(ichunk,ithread,j)
      enddo
    endif
    B(ichunk,ithread,i)=sum/A(ichunk,ithread,i,i)
  enddo
  
  enddo ! chunk loop
  
  return

end subroutine lubksb_chunk

! ************************************************************************** !

subroutine Interpolate(x_high,x_low,x,y_high,y_low,y)
  ! 
  ! Interpolates values between two reference values
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/09/09
  ! 

  implicit none

  PetscReal :: x_high, x_low, x
  PetscReal :: y_high, y_low, y
  
  PetscReal :: weight
  PetscReal :: x_diff
  
  x_diff = x_high-x_low
  if (dabs(x_diff) < 1.d-10) then
    y = y_low
  else
    weight = (x-x_low)/x_diff
    y = y_low + weight*(y_high-y_low)
  endif

end subroutine Interpolate

! ************************************************************************** !

subroutine GradientLinear(x_high,x_low,y_high,y_low,dy_dx)
  ! 
  ! Computes linear gradient given two reference values
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/12/18
  ! 

  implicit none

  PetscReal, intent(in) :: x_high, x_low
  PetscReal, intent(in) :: y_high, y_low
  PetscReal, intent(out) :: dy_dx
  
  PetscReal :: x_diff
  
  x_diff = x_high-x_low
  if (dabs(x_diff) < 1.d-10) then
    dy_dx = 0.0
  else
    dy_dx = (y_high - y_low) / x_diff
  endif

end subroutine GradientLinear

! ************************************************************************** !

function InterpolateBilinear(x,y,x1,x2,y1,y2,z1,z2,z3,z4)
  ! 
  ! Interpolates values between four reference values in 2D
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  implicit none

  PetscReal :: x,y,x1,x2,y1,y2,z1,z2,z3,z4
  PetscReal :: InterpolateBilinear
  
  
  !  x1,y2,z3 ------ x2,y2,z4
  !     |               |
  !     |               |
  !     |   x,y         |
  !     |               |
  !  x1,y1,z1 ------ x2,y1,z2
  
  
  InterpolateBilinear = (z1*(x2-x)*(y2-y)+z2*(x-x1)*(y2-y)+ &
                         z3*(x2-x)*(y-y1)+z4*(x-x1)*(y-y1))/ &
                        ((x2-x1)*(y2-y1))

end function InterpolateBilinear

! ************************************************************************** !

function DotProduct1(v1,v2)
  ! 
  ! Computes the dot product between two 3d vectors
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/28/07
  ! 

  implicit none
  
  PetscReal :: v1(3), v2(3)
  
  PetscReal :: DotProduct1
  
  DotProduct1 = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

end function DotProduct1

! ************************************************************************** !

function DotProduct2(v1,v2x,v2y,v2z)
  ! 
  ! Computes the dot product between two 3d vectors
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/28/07
  ! 

  implicit none
  
  PetscReal :: v1(3), v2x, v2y, v2z
  
  PetscReal :: DotProduct2
  
  DotProduct2 = v1(1)*v2x+v1(2)*v2y+v1(3)*v2z

end function DotProduct2

! ************************************************************************** !

function DotProduct3(v1x,v1y,v1z,v2x,v2y,v2z)
  ! 
  ! Computes the dot product between components of two 3d
  ! vectors
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/28/07
  ! 

  implicit none
  
  PetscReal :: v1x, v1y, v1z, v2x, v2y, v2z
  
  PetscReal :: DotProduct3
  
  DotProduct3 = v1x*v2x+v1y*v2y+v1z*v2z

end function DotProduct3

! ************************************************************************** !

function CrossProduct1(v1,v2)
  ! 
  ! Computes the cross product between two 3d vectors
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/09
  ! 

  implicit none
  
  PetscReal :: v1(3), v2(3)
  
  PetscReal :: CrossProduct1(3)
  
  CrossProduct1(1) = v1(2)*v2(3)-v1(3)*v2(2)
  CrossProduct1(2) = v1(3)*v2(1)-v1(1)*v2(3)
  CrossProduct1(3) = v1(1)*v2(2)-v1(2)*v2(1)

end function CrossProduct1

! ************************************************************************** !

function Erf_(x)
  ! 
  ! Computes an approximate to erf(x)
  ! from: http://jin.ece.uiuc.edu/routines/merror.for
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/20/09
  ! 

  implicit none
  
  PetscReal :: x
  
  PetscReal :: Erf_

  PetscReal, parameter :: EPS = 1.d-15
  PetscReal, parameter :: PI=3.141592653589793d0
  PetscReal :: x2, er, r, co
  PetscInt :: k
  
  x2=x*x
  if (dabs(x) < 3.5d0) then
    er=1.d0
    r=1.d0
    do k = 1, 50
      r=r*x2/(dble(k)+0.5d0)
      er=er+r
      if (dabs(r) < dabs(er)*EPS) exit
    enddo
    co=2.d0/sqrt(PI)*x*exp(-x2)
    Erf_=co*er
  else
    er=1.d0
    r=1.d0
    do k = 1, 12
      r=-r*(dble(k)-0.5d0)/x2
      er=er+r
    enddo
    co=exp(-x2)/(dabs(x)*sqrt(PI))
    Erf_=1.d0-co*er
    if (x < 0.d0) Erf_=-Erf_
  endif

end function Erf_

! ************************************************************************** !

subroutine InverseNorm(p,invnormdist,calculate_derivative,dinvnormdist_dp)
  ! This function returns the scaled inverse normal distribution
  ! which can be related to the inverse complementary error function.
  ! input range: 0 < x < 2
  ! erfc^{-1}(x) = -InverseNorm(x/2)/sqrt(2.0)
  !
  ! adapted from
  ! #
  ! # Lower tail quantile for standard normal distribution function.
  ! #
  ! # This function returns an approximation of the inverse cumulative
  ! # standard normal distribution function.  I.e., given P, it returns
  ! # an approximation to the X satisfying P = Pr{Z <= X} where Z is a
  ! # random variable from the standard normal distribution.
  ! #
  ! # The algorithm uses a minimax approximation by rational functions
  ! # and the result has a relative error whose absolute value is less
  ! # than 1.15e-9.
  ! #
  ! # Author:      Peter J. Acklam
  ! # Time-stamp:  2000-07-19 18:26:14
  ! # E-mail:      pjacklam@online.no
  ! # WWW URL:     http://home.online.no/~pjacklam
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/20/09
  ! 

  implicit none
  
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: invnormdist
  PetscBool, intent(in) :: calculate_derivative
  PetscReal, intent(out) :: dinvnormdist_dp

  PetscReal :: X, Z
  PetscReal :: dX_dq, dZ_dq
  PetscReal :: dX_dr, dZ_dr
  PetscReal :: dr_dq
  PetscReal :: dq_dp
  
 ! Coefficients in rational approximations.
  PetscReal, parameter :: A(6) = (/-3.969683028665376d+1,2.209460984245205d+2, &
                                  -2.759285104469687d+2,1.383577518672690d+2, &
                                  -3.066479806614716d+1,2.506628277459239d+0/)
  PetscReal, parameter :: B(5) = (/-5.447609879822406d+1,1.615858368580409d+2, &
                                  -1.556989798598866d+2,6.680131188771972d+1, &
                                  -1.328068155288572d+1/)
  PetscReal, parameter :: C(6) = (/-7.784894002430293d-3, &
                                  -3.223964580411365d-1, &
                                  -2.400758277161838d+0,-2.549732539343734d+0, &
                                  4.374664141464968d+0,2.938163982698783d+0/)
  PetscReal, parameter :: D(4) = (/7.784695709041462d-03, &
                                  3.224671290700398d-01, &
                                  2.445134137142996d+00,  3.754408661907416d+0/)

  ! Define break-points.
  PetscReal, parameter :: PLOW  = 0.02425d0;
  PetscReal, parameter :: PHIGH = 0.97575d0 ! 1 - PLOW;
  PetscReal :: q, r

  ! Rational approximation for lower region:
  if (p < PLOW) then
    q = sqrt(-2.d0*log(p))
    X = (((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6))
    Z = ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1.d0)
    invnormdist = X/Z
    if (calculate_derivative) then
      dq_dp = -1.d0/(q*p)
      dX_dq = (((5.d0*C(1)*q+4.d0*C(2))*q+3.d0*C(3))*q+2.d0*C(4))*q+C(5)
      dZ_dq = ((4.d0*D(1)*q+3.d0*D(2))*q+2.d0*D(3))*q+D(4)
      dinvnormdist_dp = (dX_dq/Z-invnormdist/Z*dZ_dq)*dq_dp
    endif
  ! Rational approximation for upper region:
  elseif (PHIGH < p) then
    q = sqrt(-2.d0*log(1.d0-p))
    X = (((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6))
    Z = ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1.d0)
    invnormdist = -1.d0*X/Z
    if (calculate_derivative) then
      dq_dp = 1.d0/(q*(1.d0-p))
      dX_dq = (((5.d0*C(1)*q+4.d0*C(2))*q+3.d0*C(3))*q+2.d0*C(4))*q+C(5)
      dZ_dq = ((4.d0*D(1)*q+3.d0*D(2))*q+2.d0*D(3))*q+D(4)
      dinvnormdist_dp = -1.d0*(dX_dq/Z+invnormdist/Z*dZ_dq)*dq_dp
    endif
  ! Rational approximation for central region:
  else
    q = p - 0.5d0;
    r = q*q;
    X = (((((A(1)*r+A(2))*r+A(3))*r+A(4))*r+A(5))*r+A(6))
    Z = (((((B(1)*r+B(2))*r+B(3))*r+B(4))*r+B(5))*r+1.d0)
    invnormdist = X*q/Z
    if (calculate_derivative) then
      dq_dp = 1.d0
      dr_dq = 2.d0*q
      dX_dr = (((5.d0*A(1)*r+4.d0*A(2))*r+3.d0*A(3))*r+2.d0*A(4))*r+A(5)
      dZ_dr = (((5.d0*B(1)*r+4.d0*B(2))*r+3.d0*B(3))*r+2.d0*B(4))*r+B(5)
      dinvnormdist_dp = (invnormdist/q+ &
                         (dX_dr*q/Z-invnormdist/Z*dZ_dr)*dr_dq)* &
                        dq_dp
    endif
  endif

end subroutine InverseNorm

! ************************************************************************** !

subroutine UtilityReadIntArray(array,array_size,comment,input,option)
  ! 
  ! Reads an array of integers from an input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/30/11
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: comment
  PetscInt :: array_size
  PetscInt, pointer :: array(:)
  
  PetscInt :: i, num_values, icount
  type(input_type), pointer :: input2
  character(len=MAXSTRINGLENGTH) :: string2
  character(len=MAXWORDLENGTH) :: word, word2, word3
  character(len=1) :: backslash
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscBool :: continuation_flag
  PetscInt :: value
  PetscInt, pointer :: temp_array(:)
  PetscInt :: temp_array_size

  err_string = trim(comment) // ',UtilityReadIntArray'
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  temp_array_size = 1000
  if (array_size > 0) then
    temp_array_size = array_size
  endif
  allocate(temp_array(temp_array_size))
  temp_array = 0
  
  input%ierr = 0
  if (len_trim(input%buf) > 0) then
    string2 = trim(input%buf)
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'file or value','UtilityReadIntArray')
    call StringToLower(word)
    if (StringCompare(word,'file',FOUR_INTEGER)) then
      call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
      input%err_buf = 'filename'
      input%err_buf2 = comment
      call InputErrorMsg(input,option)
      input2 => InputCreate(input,string2,option)
    else
      input2 => input
      input%buf = string2
    endif
  else
    input2 => input
  endif
  
  if (.not. len_trim(input2%buf) > 1) then
    call InputReadPflotranString(input2,option)
    call InputReadStringErrorMsg(input2,option,comment)
  endif
  
  icount = 0
  do

    continuation_flag = PETSC_FALSE
    if (index(input2%buf,backslash) > 0) &
      continuation_flag = PETSC_TRUE

    do 
      call InputReadWord(input2,option,word,PETSC_TRUE)
      if (InputError(input2) .or. &
          StringCompare(word,backslash,ONE_INTEGER)) exit
      i = index(word,'*')
      if (i == 0) i = index(word,'@')
      if (i /= 0) then
        word2 = word(1:i-1)
        word3 = word(i+1:len_trim(word))
        string2 = word2
        call InputReadInt(string2,option,num_values,input2%ierr)
        call InputErrorMsg(input2,option,'# values',err_string)
        string2 = word3
        call InputReadInt(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value',err_string)
        do while (icount+num_values > temp_array_size)
          ! careful.  reallocateRealArray double temp_array_size every time.
          call reallocateIntArray(temp_array,temp_array_size) 
        enddo
        do i=1, num_values
          icount = icount + 1
          temp_array(icount) = value
        enddo
      else
        string2 = word
        call InputReadInt(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value',err_string)
        icount = icount + 1
        if (icount > temp_array_size) then
          ! careful.  reallocateRealArray double temp_array_size every time.
          call reallocateIntArray(temp_array,temp_array_size) 
        endif
        temp_array(icount) = value
      endif
    enddo

    if (continuation_flag) then
      call InputReadPflotranString(input2,option)
      call InputReadStringErrorMsg(input2,option,comment)
    else
      if (array_size > 0) then
        if (icount == 1) then
          temp_array = temp_array(icount)
        else if (icount > array_size .or. icount < array_size) then
          write(word,*) icount
          write(word2,*) array_size
          if (len_trim(comment) > 0) then
            option%io_buffer = 'Incorrect number of values read in &
              &UtilityReadIntArray() for ' // trim(comment) // '.'
          else
            option%io_buffer = 'Incorrect number of values read in &
              &UtilityReadIntArray().'
          endif
          option%io_buffer = trim(option%io_buffer) // &
            '  Expected ' // trim(adjustl(word2)) // &
            ' but read ' // trim(adjustl(word)) // '.'
          call printErrMsg(option)
        endif
        exit
      else if (icount == 0) then
        if (len_trim(comment) > 0) then
          option%io_buffer = 'No values read in UtilityReadIntArray() &
            &for ' // trim(comment) // '.'
        else
          option%io_buffer = 'No values read in UtilityReadIntArray().'
        endif
        call printErrMsg(option)
      else
        exit
      endif
    endif
  enddo
  
  if (array_size > 0) icount = array_size
  
  if (.not.associated(input2,input)) call InputDestroy(input2)
  nullify(input2)
  
  if (associated(array)) deallocate(array)
  allocate(array(icount))
  array(1:icount) = temp_array(1:icount)
  deallocate(temp_array)
  nullify(temp_array)

end subroutine UtilityReadIntArray

! ************************************************************************** !

subroutine UtilityReadRealArray(array,array_size,comment,input,option)
  ! 
  ! Reads an array of double precision numbers from the
  ! input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/21/09
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: comment
  PetscInt :: array_size
  PetscReal, pointer :: array(:)
  
  PetscInt :: i, num_values, icount
  type(input_type), pointer :: input2
  character(len=MAXSTRINGLENGTH) :: string2
  character(len=MAXWORDLENGTH) :: word, word2, word3
  character(len=1) :: backslash
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscBool :: continuation_flag
  PetscReal :: value
  PetscReal, pointer :: temp_array(:)
  PetscInt :: temp_array_size

  err_string = trim(comment) // ',UtilityReadRealArray'
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  temp_array_size = 1000
  if (array_size > 0) then
    temp_array_size = array_size
  endif
  allocate(temp_array(temp_array_size))
  temp_array = 0.d0
  
  input%ierr = 0
  if (len_trim(input%buf) > 0) then
    string2 = trim(input%buf)
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'file or value','UtilityReadRealArray')
    call StringToLower(word)
    if (StringCompare(word,'file',FOUR_INTEGER)) then
      call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
      input%err_buf = 'filename'
      input%err_buf2 = comment
      call InputErrorMsg(input,option)
      input2 => InputCreate(input,string2,option)
    else
      input2 => input
      input%buf = string2
    endif
  else
    input2 => input
  endif
  
  if (.not. len_trim(input2%buf) > 1) then
    call InputReadPflotranString(input2,option)
    call InputReadStringErrorMsg(input2,option,comment)
  endif

  icount = 0
  do

    continuation_flag = PETSC_FALSE
    if (index(input2%buf,backslash) > 0) &
      continuation_flag = PETSC_TRUE

    do 
      call InputReadWord(input2,option,word,PETSC_TRUE)
      if (InputError(input2) .or. &
          StringCompare(word,backslash,ONE_INTEGER)) exit
      i = index(word,'*')
      if (i == 0) i = index(word,'@')
      if (i /= 0) then
        word2 = word(1:i-1)
        word3 = word(i+1:len_trim(word))
        string2 = word2
        call InputReadInt(string2,option,num_values,input2%ierr)
        call InputErrorMsg(input2,option,'# values',err_string)
        string2 = word3
        call InputReadDouble(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value',err_string)
        do while (icount+num_values > temp_array_size)
          ! careful.  reallocateRealArray double temp_array_size every time.
          call reallocateRealArray(temp_array,temp_array_size) 
        enddo
        do i=1, num_values
          icount = icount + 1
          temp_array(icount) = value
        enddo
      else
        string2 = word
        call InputReadDouble(string2,option,value,input2%ierr)
        call InputErrorMsg(input2,option,'value',err_string)
        icount = icount + 1
        if (icount > temp_array_size) then
          ! careful.  reallocateRealArray double temp_array_size every time.
          call reallocateRealArray(temp_array,temp_array_size) 
        endif
        temp_array(icount) = value
      endif
    enddo

    if (continuation_flag) then
      call InputReadPflotranString(input2,option)
      call InputReadStringErrorMsg(input2,option,comment)
    else
      if (array_size > 0) then
        if (icount == 1) then
          temp_array = temp_array(icount)
        else if (icount > array_size .or. icount < array_size) then
          write(word,*) icount
          write(word2,*) array_size
          if (len_trim(comment) > 0) then
            option%io_buffer = 'Incorrect number of values read in &
              &UtilityReadRealArray() for ' // trim(comment) // '.'
          else
            option%io_buffer = 'Incorrect number of values read in &
              &UtilityReadRealArray().'
          endif
          option%io_buffer = trim(option%io_buffer) // &
            '  Expected ' // trim(adjustl(word2)) // &
            ' but read ' // trim(adjustl(word)) // '.'
          call printErrMsg(option)
        endif
        exit
      else if (icount == 0) then
        if (len_trim(comment) > 0) then
          option%io_buffer = 'No values read in UtilityReadRealArray() &
            &for ' // trim(comment) // '.'
        else
          option%io_buffer = 'No values read in UtilityReadRealArray().'
        endif
        call printErrMsg(option)
      else
        exit
      endif
    endif
  enddo
  
  if (array_size > 0) icount = array_size
  
  if (.not.associated(input2,input)) call InputDestroy(input2)
  nullify(input2)
  
  if (associated(array)) deallocate(array)
  allocate(array(icount))
  array(1:icount) = temp_array(1:icount)
  deallocate(temp_array)
  nullify(temp_array)

end subroutine UtilityReadRealArray

! ************************************************************************** !

function SearchOrderedArray(array,array_length,int_value)
  ! 
  ! Locates an integer value in an ordered array and
  ! returned the index
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/09
  ! 

  implicit none

  PetscInt :: array_length 
  PetscInt :: array(array_length)
  PetscInt :: int_value

  PetscInt :: SearchOrderedArray
  PetscInt :: i
  PetscInt :: array_value
  PetscInt :: upper_bound, lower_bound

  SearchOrderedArray = -1

  upper_bound = array_length
  lower_bound = 1

  i = array_length/2
  if (i == 0) i = 1

  do 
    array_value = array(i)
    if (array_value == int_value) then
      SearchOrderedArray = i
      return
    endif
    if (array_value > int_value) then
      upper_bound = i 
    else
      lower_bound = i 
    endif
    i = lower_bound + (upper_bound-lower_bound) / 2
    if (i == lower_bound) then
      if (array(lower_bound) == int_value) SearchOrderedArray = lower_bound
      if (array(upper_bound) == int_value) SearchOrderedArray = upper_bound
      return
    endif
  enddo

end function SearchOrderedArray

! ************************************************************************** !

function FileExists(filename)
  ! 
  ! Returns PETSC_TRUE if file exists
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/27/11
  ! 

  implicit none
  
  PetscBool :: FileExists

  character(len=*) :: filename
  
  inquire(file=filename,exist=FileExists)

end function FileExists

! ************************************************************************** !

function Equal(value1, value2)
  ! 
  ! Returns PETSC_TRUE if values are equal
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/27/11
  ! 

  implicit none
  
  PetscBool :: Equal

  PetscReal :: value1, value2

  Equal = PETSC_FALSE
  ! using "abs(x) < spacing(y)/2.0" consistently gives same response as "x == y" for reals
  ! using both gfortran and intel compilers for y around 0.0 and 1.0
  ! this is setup assuming the "correct value" is on the RHS (second arg)
  if (dabs(value1 - value2) < spacing(value2)/2.0)  Equal = PETSC_TRUE
  
end function Equal

! ************************************************************************** !

function BestFloat(float,upper_bound,lower_bound)
  ! 
  ! Returns the best format for a floating point number
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/21/11
  ! 

  implicit none
  
  PetscReal :: float
  PetscReal :: upper_bound
  PetscReal :: lower_bound

  character(len=MAXWORDLENGTH) :: BestFloat
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i
  
100 format(f12.3)
101 format(es12.2)
102 format(es12.4)

  if (dabs(float) <= upper_bound .and. dabs(float) >= lower_bound) then
    write(word,100) float
    word = adjustl(word)
    do i = len_trim(word), 1, -1
      if (word(i:i) == '0') then
        word(i:i) = ' '
      else
        exit
      endif
    enddo
  else if (dabs(float) < lower_bound) then
    write(word,101) float
  else
    write(word,102) float
  endif
  
  BestFloat = adjustl(word)
  
end function BestFloat

! ************************************************************************** !

subroutine QuadraticPolynomialSetup(value_1,value_2,coefficients, &
                                    derivative_at_1)
  ! 
  ! Sets up a quadratic polynomial for smoothing discontinuous functions
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/25/14
  ! 

  implicit none

  PetscReal :: value_1
  PetscReal :: value_2
  PetscReal :: coefficients(3)
  PetscBool :: derivative_at_1
  
  PetscReal :: A(3,3)
  PetscInt :: indx(3)
  PetscInt :: d

  A(1,1) = 1.d0
  A(2,1) = 1.d0
  A(3,1) = 0.d0
  
  A(1,2) = value_1
  A(2,2) = value_2
  A(3,2) = 1.d0
  
  A(1,3) = value_1**2.d0
  A(2,3) = value_2**2.d0
  if (derivative_at_1) then
    A(3,3) = 2.d0*value_1
  else
    A(3,3) = 2.d0*value_2
  endif
  
  ! coefficients(1): value at 1
  ! coefficients(2): value at 2
  ! coefficients(3): derivative at 1 or 2
  
  call ludcmp(A,THREE_INTEGER,indx,d)
  call lubksb(A,THREE_INTEGER,indx,coefficients)

end subroutine QuadraticPolynomialSetup

! ************************************************************************** !

subroutine QuadraticPolynomialEvaluate(coefficients,x,f,df_dx)
  ! 
  ! Evaluates value in quadratic polynomial
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/12/12
  ! 

  implicit none

  PetscReal :: coefficients(3)
  PetscReal :: x
  PetscReal :: f
  PetscReal :: df_dx

  f = coefficients(1) + &
      coefficients(2)*x + &
      coefficients(3)*x*x
  
  df_dx = coefficients(2) + &
          coefficients(3)*2.d0*x
  
end subroutine QuadraticPolynomialEvaluate

! ************************************************************************** !

subroutine CubicPolynomialSetup(value_1,value_2,coefficients)
  ! 
  ! Sets up a cubic polynomial for smoothing discontinuous functions
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/12/12

  implicit none

  PetscReal :: value_1
  PetscReal :: value_2
  PetscReal :: coefficients(4)
  
  PetscReal :: A(4,4)
  PetscInt :: indx(4)
  PetscInt :: d

  A(1,1) = 1.d0
  A(2,1) = 1.d0
  A(3,1) = 0.d0
  A(4,1) = 0.d0
  
  A(1,2) = value_1
  A(2,2) = value_2
  A(3,2) = 1.d0
  A(4,2) = 1.d0
  
  A(1,3) = value_1**2.d0
  A(2,3) = value_2**2.d0
  A(3,3) = 2.d0*value_1
  A(4,3) = 2.d0*value_2
  
  A(1,4) = value_1**3.d0
  A(2,4) = value_2**3.d0
  A(3,4) = 3.d0*value_1**2.d0
  A(4,4) = 3.d0*value_2**2.d0
  
  ! coefficients(1): value at 1
  ! coefficients(2): value at 2
  ! coefficients(3): derivative at 1
  ! coefficients(4): derivative at 2
  
  call ludcmp(A,FOUR_INTEGER,indx,d)
  call lubksb(A,FOUR_INTEGER,indx,coefficients)

end subroutine CubicPolynomialSetup

! ************************************************************************** !

subroutine CubicPolynomialEvaluate(coefficients,x,f,df_dx)
  ! 
  ! Evaluates value in cubic polynomial
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/12/12
  ! 

  implicit none

  PetscReal :: coefficients(4)
  PetscReal :: x
  PetscReal :: f
  PetscReal :: df_dx

  PetscReal :: x_squared
  
  x_squared = x*x
  
  f = coefficients(1) + &
      coefficients(2)*x + &
      coefficients(3)*x_squared + &
      coefficients(4)*x_squared*x
  
  df_dx = coefficients(2) + &
          coefficients(3)*2.d0*x + &
          coefficients(4)*3.d0*x_squared
  
end subroutine CubicPolynomialEvaluate

! ************************************************************************** !

subroutine DeallocateArray1DInteger(array)
  ! 
  ! Deallocates a 1D integer array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscInt, pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DInteger

! ************************************************************************** !

subroutine DeallocateArray2DInteger(array)
  ! 
  ! Deallocates a 2D integer array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscInt, pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DInteger

! ************************************************************************** !

subroutine DeallocateArray3DInteger(array)
  ! 
  ! Deallocates a 3D integer array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscInt, pointer :: array(:,:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray3DInteger

! ************************************************************************** !

subroutine DeallocateArray1DReal(array)
  ! 
  ! Deallocates a 1D real array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscReal, pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DReal

! ************************************************************************** !

subroutine DeallocateArray2DReal(array)
  ! 
  ! Deallocates a 2D real array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscReal, pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DReal

! ************************************************************************** !

subroutine DeallocateArray3DReal(array)
  ! 
  ! Deallocates a 3D real array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscReal, pointer :: array(:,:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray3DReal

! ************************************************************************** !

subroutine DeallocateArray1DLogical(array)
  ! 
  ! Deallocates a 1D logical array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscBool, pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DLogical

! ************************************************************************** !

subroutine DeallocateArray2DLogical(array)
  ! 
  ! Deallocates a 2D logical array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscBool, pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DLogical

! ************************************************************************** !

subroutine DeallocateArray3DLogical(array)
  ! 
  ! Deallocates a 3D logical array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  PetscBool, pointer :: array(:,:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray3DLogical

! ************************************************************************** !

subroutine DeallocateArray1DString(array)
  ! 
  ! Deallocates a 1D array of character strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: array(:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray1DString

! ************************************************************************** !

subroutine DeallocateArray2DString(array)
  ! 
  ! Deallocates a 2D array of character strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 

  implicit none
  
  character(len=MAXWORDLENGTH), pointer :: array(:,:)
  
  if (associated(array)) deallocate(array)
  nullify(array)

end subroutine DeallocateArray2DString

! ************************************************************************** !

subroutine ConvertMatrixToVector(A,vecA)
  ! 
  ! Converts a given matrix A to a vec
  ! This vec is different from PETSc Vec
  ! A = [a1 a2 a3 .... am], where ai, i = 1, m are the columns
  ! then vec(A) = [a1
  ! a2
  ! .
  ! .
  ! .
  ! am]
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 6/19/2013
  ! 

  PetscReal :: A(:,:)
  PetscReal, allocatable :: vecA(:,:)
  PetscInt :: m, n
  
  m = size(A,1)
  n = size(A,2)
  
  allocate(vecA(m*n,ONE_INTEGER))
  
  vecA = reshape(A,(/m*n,ONE_INTEGER/))

end subroutine ConvertMatrixToVector

! ************************************************************************** !

subroutine Kron(A,B,K)
  ! 
  ! Returns the Kronecker product of two matrices A, B
  ! Reference: The ubiquitous Kronecker product, by Charles F.Van Loan
  ! for basics of Kronecker product
  ! Also see wikipedia page: http://en.wikipedia.org/wiki/Kronecker_product
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 6/19/2013
  ! 

  PetscReal :: A(:,:),B(:,:)
  PetscReal, allocatable :: K(:,:)
  PetscInt :: mA,nA,mB,nB
  PetscInt :: iA,jA,iB,jB,iK,jK
  
  mA = size(A,1)
  nA = size(A,2)
  mB = size(B,1)
  nB = size(B,2)
  
  allocate(K(mA*mB,nA*nB))
  
  do iB = 1, mB
    do jB = 1, nB
      do iA = 1, mA
        do jA = 1, nA
          iK = iB + (iA-1)*mB
          jK = jB + (jA-1)*nB
          K(iK,jK) = A(iA,jA)*B(iB,jB)
        enddo
      enddo
    enddo
  enddo
  
end subroutine Kron

! ************************************************************************** !

subroutine Transposer(m,n,T)
  ! 
  ! Transposer Converts vec of a matrix to vec of its transpose
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 6/19/2013
  ! 

  PetscReal, allocatable :: T(:,:)
  PetscInt :: m,n
  PetscInt :: i,j
  PetscReal :: A(m,n)
  PetscReal, allocatable :: vecA(:,:)
  
  allocate(T(m*n,m*n))
  T = 0.d0
  
  do i = 1,m
    do j = 1,n
      A = 0.d0
      A(i,j) = 1.d0
      call ConvertMatrixToVector(transpose(A),vecA)
      T(:,i+m*(j-1)) = vecA(:,1)
      deallocate(vecA)
    enddo
  enddo
  
end subroutine Transposer

! ************************************************************************** !

subroutine Determinant(A,detA)
  ! 
  ! Determinant of a 3x3 matrix
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 6/24/2013
  ! 

  PetscReal :: A(3,3)
  PetscReal :: detA  
 
  detA = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) &
       + A(1,2)*(A(3,1)*A(2,3) - A(2,1)*A(3,3))  &
       + A(1,3)*(A(2,1)*A(3,2) - A(3,1)*A(2,2))
  

end subroutine Determinant

! ************************************************************************** !
subroutine InterfaceApproxWithDeriv(v_up, v_dn, dv_up, dv_dn, dv_up2dn, &
                                    approx_type, v_interf, &
                                    dv_interf_dv_up, dv_interf_dv_dn)
  ! 
  ! Approximates interface value and it's derivative from values specified
  ! up and down of a face based on the approximation type
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 05/05/2014
  ! 

  implicit none
  
  PetscReal, intent(in) :: v_up, v_dn
  PetscReal, intent(in) :: dv_up, dv_dn
  PetscReal, intent(in) :: dv_up2dn
  PetscInt, intent(in) :: approx_type
  PetscReal, intent(out) :: v_interf, dv_interf_dv_up, dv_interf_dv_dn

  PetscReal :: denom
  PetscReal :: eps = 1.d-15


  if (dv_up2dn > 0.d0) then
    v_interf = v_up
    dv_interf_dv_up = dv_up
    dv_interf_dv_dn = 0.d0
  else
    v_interf = v_dn
    dv_interf_dv_up = 0.d0
    dv_interf_dv_up = dv_dn
  endif

  select case (approx_type)

    case (UPWIND)
      if (dv_up2dn > 0.d0) then
        v_interf = v_up
        dv_interf_dv_up = dv_up
        dv_interf_dv_dn = 0.d0
      else
        v_interf = v_dn
        dv_interf_dv_up = 0.d0
        dv_interf_dv_dn = dv_dn
      endif

    case (HARMONIC)
      if (v_up < eps .or. v_dn < eps) then
        v_interf = 0.d0
        dv_interf_dv_up = 0.d0
        dv_interf_dv_dn = 0.d0
      else
        denom = (v_up + v_dn)
        v_interf = 2.d0*v_up*v_dn/denom
        dv_interf_dv_up = 2.d0*(denom*dv_up*v_dn - v_up*v_dn*dv_up)/(denom**2.d0)
        dv_interf_dv_dn = 2.d0*(denom*v_up*dv_dn - v_up*v_dn*dv_dn)/(denom**2.d0)
      endif

  end select

end subroutine InterfaceApproxWithDeriv

! ************************************************************************** !

subroutine InterfaceApproxWithoutDeriv(v_up, v_dn, dv_up2dn, &
                                       approx_type, v_interf)
  ! 
  ! Approximates interface value from values specified
  ! up and down of a face based on the approximation type
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 05/05/2014
  ! 

  implicit none
  
  PetscReal, intent(in) :: v_up, v_dn
  PetscReal, intent(in) :: dv_up2dn
  PetscInt, intent(in) :: approx_type
  PetscReal, intent(out) :: v_interf

  PetscReal :: dummy_in
  PetscReal :: dummy_out

  dummy_in = 1.d0

  call InterfaceApproxWithDeriv(v_up, v_dn, dummy_in, dummy_in, dv_up2dn, &
                                approx_type, v_interf, dummy_out, dummy_out)

end subroutine InterfaceApproxWithoutDeriv

! ************************************************************************** !

subroutine PrintProgressBarInt(max_value,increment,current)
  ! 
  ! Prints a piece of a progress bar to the screen based on the maximum
  ! value, the increment of progress (must be given in percent), and the
  ! current value.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/16/2016
  ! 

  implicit none
  
  PetscReal :: max_value
  PetscInt :: increment
  PetscInt :: current

  PetscInt :: max_value_int
  PetscInt :: g, j, chunk
  character(len=MAXWORDLENGTH) :: percent_num

  max_value_int = floor(max_value)
  if (max_value_int < increment) then
    max_value_int = max_value_int*(increment/max_value_int)
    current = current*(increment/max_value_int)
  endif

  chunk = floor(max_value_int*(increment/100.0))

  if (mod(current,chunk) == 0) then
    j = current/chunk
    g = 0
    do while(g < (100))
      g = g + increment
      if (g/increment == j) then
        write(percent_num,*) g
        write(*,'(a4)',advance='no') trim(adjustl(percent_num)) // '%-'
        if (g == 100) write(*,*) '  Done.'
      endif
    enddo   
  endif

end subroutine PrintProgressBarInt

! ************************************************************************** !

function DigitsOfAccuracy(num1,num2)

  implicit none
  
  PetscReal :: num1
  PetscReal :: num2
  
  character(len=2) :: DigitsOfAccuracy
  
  PetscReal :: tempreal
  PetscReal :: relative_difference
  PetscInt :: tempint
  
  DigitsOfAccuracy = ' 0'
  if (dabs(num1) > 0.d0 .and. dabs(num2) > 0.d0) then
    relative_difference = dabs((num1-num2)/num2)
    if (relative_difference < 1.d-17) then
      ! accuracy is beyond double precision
      DigitsOfAccuracy = '99'
    else
      tempreal = 1.d0 / relative_difference
      tempint = 0
      do
        if (tempreal < 10.d0) exit
        tempreal = tempreal / 10.d0
        tempint = tempint + 1
      enddo
      write(DigitsOfAccuracy,'(i2)') tempint
    endif
  else if (dabs(num1) > 0.d0 .or. dabs(num2) > 0.d0) then
    ! change this value if you want to report something difference for
    ! either one being zero.
  else
    ! change this value if you want to report something difference for
    ! double zeros.
    DigitsOfAccuracy = '  '
  endif
    
end function DigitsOfAccuracy

! ************************************************************************** !

subroutine CalcParallelSUM1(option,rank_list,local_val,global_sum)
  ! 
  ! Calculates global sum for a MPI_DOUBLE_PRECISION number (local_val).
  ! This function uses only MPI_Send and MPI_Recv functions and does not need 
  ! a communicator object other than option%mycomm. It reduces communication 
  ! to the processes that are included in the rank_list array rather than using
  ! a call to MPI_Allreduce.
  ! 
  ! Author: Jenn Frederick
  ! Date: 07/05/2017
  
  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option
  PetscInt :: rank_list(:)
  PetscReal :: local_val
  PetscReal :: global_sum
  
  PetscReal :: passed_local_val(1)
  PetscReal :: passed_global_val(1)
  
  passed_local_val(1) = local_val
  passed_global_val(1) = 0.d0
  
  call CalcParallelSUM2(option,rank_list,passed_local_val,passed_global_val)
  
  global_sum = passed_global_val(1)

end subroutine CalcParallelSUM1

! ************************************************************************** !

subroutine CalcParallelSUM2(option,rank_list,local_val,global_sum)
  ! 
  ! Calculates global sum for a MPI_DOUBLE_PRECISION number (local_val).
  ! This function uses only MPI_Send and MPI_Recv functions and does not need 
  ! a communicator object other than option%mycomm. It reduces communication 
  ! to the processes that are included in the rank_list array rather than using
  ! a call to MPI_Allreduce.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/23/17, 07/05/2017
  
  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option
  PetscInt :: rank_list(:)
  PetscReal :: local_val(:)
  PetscReal :: global_sum(:)

  PetscReal, pointer :: temp_array(:)
  PetscInt :: num_ranks, val_size
  PetscInt :: m, j
  PetscInt :: TAG
  PetscErrorCode :: ierr
  
  num_ranks = size(rank_list)
  val_size = size(local_val)
  allocate(temp_array(num_ranks))
  TAG = 0
  
  if (num_ranks > 1) then
  !------------------------------------------
    temp_array = 0.d0
    do j = 1,val_size
  
      if (option%myrank .ne. rank_list(1)) then
        call MPI_Send(local_val(j),ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      rank_list(1),TAG,option%mycomm,ierr)
      else
        temp_array(1) = local_val(j)
        do m = 2,num_ranks
          call MPI_Recv(local_val(j),ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                        rank_list(m),TAG,option%mycomm,MPI_STATUS_IGNORE,ierr)
          temp_array(m) = local_val(j)
        enddo
        global_sum(j) = sum(temp_array)
      endif
      if (option%myrank == rank_list(1)) then
        do m = 2,num_ranks
          call MPI_Send(global_sum(j),ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                        rank_list(m),TAG,option%mycomm,ierr)
        enddo
      else
        call MPI_Recv(global_sum(j),ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                      rank_list(1),TAG,option%mycomm,MPI_STATUS_IGNORE,ierr)
      endif             
    
    enddo
  !------------------------------------------        
  else 
    global_sum = local_val
  endif
  
  deallocate(temp_array)

end subroutine CalcParallelSUM2

! ************************************************************************** !
! something like Nathan Collier's H Function or Guoping Tang's tail cut-off approach
! F.-M. Yuan: 2017-03-09. It's useful for truncation or interpolation, with monotonic smoothing curve
Subroutine HfunctionSmooth(x, x_1, x_0, H, dH)
  ! Usage: H(x) = 1, dH(x)=0 @ x_1
  !        H(x) = 0, dH(x)=0 @ x_0
  !        and, dH(x) is smoothed and either positive or negative Heaveside function, with peak at ~ 3-quater point
  !    So H(x) is monotonic.
  !
  !  How to tail-smooth a curve (e.g. f(x)) monotonically:
  !       F(x) = f(x)*H(f), and dF(x) = f(x)*dH(x)+ df(x)*H(x),
  !      From: x_cutoff1 - F(x)=f(x), dF(x)=df(x), i.e. exactly matching with f(x) @ x_cutoff1
  !      To:   x_cutoff0 - F(x)=0, dF(x)=0.
  !
  implicit none

  PetscReal, intent(in) :: x
  PetscReal, intent(in) :: x_1    ! starting point with H = 1.0
  PetscReal, intent(in) :: x_0    ! ending point with H = 0.0
  PetscReal, intent(out):: H
  PetscReal, intent(out):: dH
  PetscReal :: x_star

  !----------------------------------------------------------

  ! real Heaveside function (step function)
  if (abs(x_1-x_0)<1.d-50) then
    H  = sign(0.5d0, (x-x_1))+0.5d0  ! if x<x_1, H=0; otherwise H=1
    dH = 0.d0

    return

  else
  ! smoothed H function

    if (((x-x_0)/(x_1-x_0))<0.d0) then      ! beyond 'x_0'
      H  = 0.0d0
      dH = 0.0d0

    elseif (((x-x_0)/(x_1-x_0))>1.d0) then  ! beyond 'x_1'
      H  = 1.0d0
      dH = 0.0d0

    else

      x_star  =  1.0d0 - (x-x_0)*(x-x_0)                  &
                        /(x_1-x_0)/(x_1-x_0)

      H  = 1.0d0 - x_star*x_star                                      ! so it's a special quadratic polynomal function
      ! so, dH = -2*x_star*d(x_star), with d(x_star)=-2*(x-x_0)/(x_1-x_0)^2
      !    and, due to 'x_star' ranging 0~1 and (x_1-x_0)^2 positive,
      !         sign(dH) is upon 'x-x_0' only guaranted monotonic H curve
      dH = 4.0d0 * x_star * (x-x_0)/(x_1-x_0)/(x_1-x_0)

    endif

  endif

end subroutine HfunctionSmooth

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !SUBROUTINE: where_checkerr(ierr)
  !
  ! !INTERFACE:
  subroutine where_checkerr(ierr, subname, filename, line)
  !
  ! !DESCRIPTION:
  ! When using PETSc functions, it usually throws an error code for checking.
  ! BUT it won't show where the error occurs in the first place, therefore it's hardly useful.
  !
  ! !USES:

    implicit none

  ! !ARGUMENTS:
    character(len=*), intent(IN) :: subname  ! subroutine name called this
    character(len=*), intent(IN) :: filename ! filename called this
    integer, intent(IN) :: line              ! line number triggered this
    PetscErrorCode, intent(IN) :: ierr       ! petsc error code

  !EOP
  !-----------------------------------------------------------------------

    if (ierr /= 0) then
       print *, ' PETSc ERROR: @Subroutine - ' // trim(subname)
       print *, ' PETSc ERROR: @File - ' // trim(filename)
       print *, ' PETSc ERROR: @Line -', line
    end if
    CHKERRQ(ierr)

  end subroutine where_checkerr

!--------------------------------------------------------------------------------------

! ************************************************************************** !

subroutine MatCompare(a1, a2, n, m, tol, do_rel_err)

  !! Daniel Stone, March 2018
  !! Just output warnings and provide place
  !! for breakpoints.
  !! Used in testing analytical derivatives
  !! and comparing with numerical.

  implicit none
  PetscInt :: n, m
  PetscReal, dimension(1:n, 1:m) :: a1, a2
  PetscReal :: tol
  PetscBool :: do_rel_err 

  PetscInt :: i, j
  PetscReal :: dff

  do i = 1,n
    do j = 1,m
      dff = abs(a1(i,j) - a2(i,j)) 
      if (do_rel_err) then
        dff = dff/abs(a1(i,j))
      endif
      if (dff > tol) then
        print *, "difference in matrices at ", i, ", ", j, ", value ", dff
        print *, a1(i,j), " compare to ", a2(i,j)
        print *, "..."
      endif
    end do
  end do 

end subroutine MatCompare

! ************************************************************************** !

end module Utility_module
