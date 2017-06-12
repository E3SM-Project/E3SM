!-----------------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! Copyright 2017; all rights reserved
!=======================================================================

module FortranVector
  !-----------------------------------------------------------------------
  ! Description: simple Fortran user-defined type for example interface
  !-----------------------------------------------------------------------

  use kinds,                  only: real_kind, long_kind, int_kind
  use dimensions_mod,         only: np, nlev, nelem

  implicit none
  private
  integer, public, parameter :: timelevels = 3

! =========== PRIMITIVE-EQUATION DATA-STRUCTURES =====================
  type, public :: Fvec
  

    real (kind=real_kind) :: u              ! horizontal velocity x-component
    real (kind=real_kind) :: v              ! horizontal velocity y-component
    real (kind=real_kind) :: w              ! vertical velocity
    real (kind=real_kind) :: theta_dp_cp    ! potential temperature
    real (kind=real_kind) :: phi            ! geopotential
    real (kind=real_kind) :: dp3d           ! delta p on levels

  end type Fvec
  !=======================================================================
end module FortranVector


#ifdef UNUSED_EXAMPLE_CODE

subroutine FNVExtPrint(x)
  !-----------------------------------------------------------------------
  ! Print routine for EXT vector
  ! 
  ! Note: this function is not required by ARKode (or any of SUNDIALS) -- 
  !       it is merely here for convenience when debugging
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(Fvec), intent(inout), target :: x(:)

  !=======Internals ============

  ! print vector data
  print *, "wassup", size(x)

  return
end subroutine FNVExtPrint
!=======================================================================



subroutine FNVExtClone(x, y_C)
  !-----------------------------------------------------------------------
  ! Clone routine for EXT vector -- allocates memory for a new FVec y that 
  ! that matches structure for x.
  !-----------------------------------------------------------------------
  use FortranVector
  use iso_c_binding
  implicit none
  type(FVec), intent(in)   :: x(:)
  type(C_PTR), intent(out) :: y_C
  type(Fvec), pointer      :: y(:)

  !=======Internals ============

  integer n,m
  ! allocate new vector (in this demonstration code, there's no need to 
  ! actually look at x since we know the structure of these vectors)
  n = size(x)
  ! initialize values to zero
  do m=1,n
   y(m)%u           = 0.d0
   y(m)%v           = 0.d0
   y(m)%w           = 0.d0
   y(m)%phi         = 0.d0
   y(m)%theta_dp_cp = 0.d0
   y(m)%dp3d        = 0.d0
  end do 
  ! associate y_C with yvec
  y_C = c_loc(y)

  return
end subroutine FNVExtClone
!=======================================================================



subroutine FNVExtDestroy(x)
  !-----------------------------------------------------------------------
  ! Destroy routine for EXT vector -- deallocates memory for a given FVec.
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), pointer :: x(:)
  integer :: ierr

  !=======Internals ============

  ! deallocate memory for x
  deallocate(x, stat=ierr)
  if (ierr /= 0)  write(0,*) 'FNVExtDestroy: deallocation request denied'
  
  return
end subroutine FNVExtDestroy
!=======================================================================



subroutine FNVExtLinearSum(aval, x, bval, y, z)
  !-----------------------------------------------------------------------
  ! z = a*x + b*y
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  real*8, intent(in)              :: aval
  type(FVec), intent(in), pointer :: x(:)
  real*8, intent(in)              :: bval
  type(FVec), intent(in)          :: y(:)
  type(FVec), intent(inout)       :: z(:)

  !=======Internals ============
  integer :: n,m
  ! perform vector operation
  n = size(x)
  do m=1,n
    z(m)%u           = aval*x(m)%u           + bval*y(m)%u
    z(m)%v           = aval*x(m)%v           + bval*y(m)%v
    z(m)%w           = aval*x(m)%w           + bval*y(m)%w
    z(m)%phi         = aval*x(m)%phi         + bval*y(m)%phi
    z(m)%theta_dp_cp = aval*x(m)%theta_dp_cp + bval*y(m)%theta_dp_cp
    z(m)%dp3d        = aval*x(m)%dp3d        + bval*y(m)%dp3d
  end do
  return
end subroutine FNVExtLinearSum
!=======================================================================



subroutine FNVExtConst(cval, z)
  !-----------------------------------------------------------------------
  ! z = c 
  !-----------------------------------------------------------------------
  use FortranVector
  use iso_c_binding
  use dimensions_mod,         only: np, nlev
  implicit none
  real*8, intent(in)        :: cval
  type(FVec), intent(inout) :: z(:)

  !=======Internals ============
  integer :: n,m
  n=size(z)
  ! perform vector operation 
  do m=1,n
    z(m)%u           = cval
    z(m)%v           = cval
    z(m)%w           = cval
    z(m)%phi         = cval
    z(m)%theta_dp_cp = cval
    z(m)%dp3d        = cval
  end do 
  return
end subroutine FNVExtConst
!=======================================================================



subroutine FNVExtProd(x, y, z)
  !-----------------------------------------------------------------------
  ! z = x.*y (Matlab notation)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in)    :: x(:)
  type(FVec), intent(in)    :: y(:)
  type(FVec), intent(inout) :: z(:)

  !=======Internals ============
  integer :: n,m
  ! perform vector operation
  n = size(x)
  do m =1,n
    z(m)%u           = x(m)%u * y(m)%u
    z(m)%v           = x(m)%v * y(m)%v
    z(m)%w           = x(m)%w * y(m)%w
    z(m)%phi         = x(m)%phi * y(m)%phi
    z(m)%theta_dp_cp = x(m)%theta_dp_cp * y(m)%theta_dp_cp
    z(m)%dp3d        = x(m)%dp3d * y(m)%dp3d
  end do

  return
end subroutine FNVExtProd
!=======================================================================



subroutine FNVExtDiv(x, y, z)
  !-----------------------------------------------------------------------
  ! z = x./y (Matlab notation; doesn't check for legal denominator)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in)    :: x(:)
  type(FVec), intent(in)    :: y(:)
  type(FVec), intent(inout) :: z(:)

  !=======Internals ============
  integer :: n,m
  
  ! perform vector operation
  n = size(x)
  do m=1,n
    z(m)%u           = x(m)%u / y(m)%u
    z(m)%v           = x(m)%v / y(m)%v
    z(m)%w           = x(m)%w / y(m)%w
    z(m)%phi         = x(m)%phi / y(m)%phi
    z(m)%theta_dp_cp = x(m)%theta_dp_cp / y(m)%theta_dp_cp
    z(m)%dp3d        = x(m)%dp3d / y(m)%dp3d
  end do

  return
end subroutine FNVExtDiv
!=======================================================================



subroutine FNVExtScale(cval, x, z)
  !-----------------------------------------------------------------------
  ! z = c*x
  !-----------------------------------------------------------------------
  use FortranVector
  use iso_c_binding
  implicit none
  real*8, intent(in)        :: cval
  type(FVec), intent(in)    :: x(:)
  type(FVec), intent(inout) :: z(:)

  !=======Internals ============
  integer :: n,m
  ! perform vector operation
  n = size(x)
  do m=1,n
    z(m)%u           = cval * x(m)%u
    z(m)%v           = cval * x(m)%v
    z(m)%w           = cval * x(m)%w
    z(m)%phi         = cval * x(m)%phi
    z(m)%theta_dp_cp = cval * x(m)%theta_dp_cp
    z(m)%dp3d        = cval * x(m)%dp3d
  end do 
  return
end subroutine FNVExtScale
!=======================================================================



subroutine FNVExtAbs(x, z)
  !-----------------------------------------------------------------------
  ! z = |x|  (componentwise)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in)    :: x(:)
  type(FVec), intent(inout) :: z(:)

  !=======Internals ============
  integer :: n,m
  ! perform vector operation
  n = size(x)
  do m=1,n
    z(m)%u           = dabs( x(m)%u)
    z(m)%v           = dabs( x(m)%v)
    z(m)%w           = dabs( x(m)%w)
    z(m)%phi         = dabs( x(m)%phi)
    z(m)%theta_dp_cp = dabs( x(m)%theta_dp_cp)
    z(m)%dp3d        = dabs( x(m)%dp3d)
  end do
  return
end subroutine FNVExtAbs
!=======================================================================



subroutine FNVExtInv(x, z)
  !-----------------------------------------------------------------------
  ! z = 1./x (Matlab notation: doesn't check for division by zero)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in)    :: x(:)
  type(FVec), intent(inout) :: z(:)

  !=======Internals ============
  integer :: n,m
  ! perform vector operation
  n = size(x)
  do m=1,n
    z(m)%u           = 1.d0 / x(m)%u
    z(m)%v           = 1.d0 / x(m)%v
    z(m)%w           = 1.d0 / x(m)%w
    z(m)%phi         = 1.d0 / x(m)%phi
    z(m)%theta_dp_cp = 1.d0 / x(m)%theta_dp_cp
    z(m)%dp3d        = 1.d0 / x(m)%dp3d
  end do 

  return
end subroutine FNVExtInv
!=======================================================================



subroutine FNVExtAddConst(cval, x, z)
  !-----------------------------------------------------------------------
  ! z = x+c (add c to each component of the vector x)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  real*8, intent(in)        :: cval
  type(FVec), intent(in)    :: x(:)
  type(FVec), intent(inout) :: z(:)

  !=======Internals ============
  integer :: n,m
  ! perform vector operation
  do m=1,n
    z(m)%u           = cval + x(m)%u
    z(m)%v           = cval + x(m)%v
    z(m)%w           = cval + x(m)%w
    z(m)%phi         = cval + x(m)%phi
    z(m)%theta_dp_cp = cval + x(m)%theta_dp_cp
    z(m)%dp3d        = cval + x(m)%dp3d
  end do

  return
end subroutine FNVExtAddConst
!=======================================================================



subroutine FNVExtDotProd(x, y, cval)
  !-----------------------------------------------------------------------
  ! c = <x,y>  (only include 'active' data; no ghost cells, etc.)
  !-----------------------------------------------------------------------
  use FortranVector
 
  implicit none
  type(FVec), intent(in) :: x(:)
  type(FVec), intent(in) :: y(:)
  real*8, intent(out)    :: cval

  !=======Internals ============
  integer :: n,m
  ! perform vector operation
  cval = 0.d0
  n =size(x)
  do m=1,n
    cval  = cval+x(m)%u*y(m)%u                  & 
            + x(m)%v*y(m)%v                     &
            + x(m)%w*y(m)%w                     &     
            + x(m)%phi*y(m)%phi                 &
            + x(m)%theta_dp_cp*y(m)%theta_dp_cp &
            + x(m)%dp3d*y(m)%dp3d 
  end do
  return
end subroutine FNVExtDotProd
!=======================================================================



subroutine FNVExtMaxNorm(x, cval)
  !-----------------------------------------------------------------------
  ! c = max(|x|)
  !-----------------------------------------------------------------------
  use FortranVector

  implicit none
  type(FVec), intent(in) :: x(:)
  real*8, intent(out)    :: cval

  !=======Internals ============
  ! perform vector operation
  integer :: n,m
  real*8 :: cvaltemp
  n=size(x)
  cval=0.d0
  do m=1,n
    cvaltemp = max( dabs(x(m)%u), max( dabs(x(m)%v),  &
      max(dabs(x(m)%w),max(dabs(x(m)%phi),            &
      max(dabs(x(m)%theta_dp_cp),dabs(x(m)%dp3d))))))
      if (cvaltemp.gt.cval) then
        cval=cvaltemp
      end if
  end do

  return
end subroutine FNVExtMaxNorm
!=======================================================================



subroutine FNVExtWrmsNorm(x, w, cval)
  !-----------------------------------------------------------------------
  ! cval = sqrt(||x.*w||^2_2 / Ntotal)
  !-----------------------------------------------------------------------
  use FortranVector

  implicit none
  type(FVec), intent(in), target :: x(:)
  type(FVec), intent(in), target :: w(:)
  real*8, intent(out)            :: cval

  !=======Internals ============
  integer :: i,j,k,n,m
  n = size(x)
  cval=0.d0
  do m=1,n
    ! perform vector operation
    cval = cval + (x(m)%u * w(m)%u)**2         + &
      (x(m)%v * w(m)%v)**2                     + &
      (x(m)%w * w(m)%w)**2                     + &
      (x(m)%phi * w(m)%phi)**2                 + &
      (x(m)%theta_dp_cp * w(m)%theta_dp_cp)**2 + &
      (x(m)%dp3d * w(m)%dp3d)**2
  end do
  cval =dsqrt(cval/(6*n))
  return
end subroutine FNVExtWrmsNorm
!=======================================================================



subroutine FNVExtWrmsNormMask(x, w, l, cval)
  !-----------------------------------------------------------------------
  ! cval = ||x.*w||_2 / Ntotal (but only where l > 0)
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL')
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in), target :: x(:)
  type(FVec), intent(in), target :: w(:)
  type(FVec), intent(in), target :: l(:)
  real*8, intent(out) :: cval
  integer :: ntot

  !=======Internals ============

  ! perform vector operation
!  cval = 0.d0
!  ntot = 0
!  if (l(:)%v(:,:,1,:) > 0.d0) then
!     cval = cval + (x(:)%v(:,:,1,:) * w(:)%v(:,:,1,:))**2
!     ntot = ntot + 1
!  endif
!  if (l(:)%v(:,:,2,:) > 0.d0) then
!     cval = cval + (x(:)%v(:,:,2,:) * w(:)%v(:,:,2,:))**2
!     ntot = ntot + 1
!  endif
!  if (l(:)%w(:,:,:) > 0.d0) then
!     cval = cval + (x(:)%w(:,:,:) * w(:)%w(:,:,:))**2
!     ntot = ntot + 1
!  endif
!  if (l(:)%phi(:,:,:) > 0.d0) then
!     cval = cval + (x(:)%phi(:,:,:) * w(:)%phi(:,:,:))**2
!     ntot = ntot + 1
!  endif
!  if (l(:)%theta_dp_cp(:,:,:) > 0.d0) then
!     cval = cval + (x(:)%theta_dp_cp(:,:,:) * w(:)%theta_dp_cp(:,:,:))**2
!     ntot = ntot + 1
!  endif
!  if (l(:)%dp3d(:,:,:) > 0.d0) then
!     cval = cval + (x(:)%dp3d(:,:,:) * w(:)%dp3d(:,:,:))**2
!     ntot = ntot + 1
!  endif


  cval = dsqrt(cval / ntot)

  return
end subroutine FNVExtWrmsNormMask
!=======================================================================



subroutine FNVExtMin(x, cval)
  !-----------------------------------------------------------------------
  ! cval = min(xvec)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in), target :: x(:)
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
!  cval = min( dabs(x(:)%v(:,:,1,:)), min( dabs(x(:)%v(:,:,1,:)),&
!    min(dabs(x(:)%w(:,:,:)),min(dabs(x(:)%phi(:,:,:)),         &
!    min(dabs(x(:)%theta_dp_cp(:,:,:)),dabs(x(:)%dp3d(:,:,:)))))))
  return
end subroutine FNVExtMin
!=======================================================================



subroutine FNVExtWl2Norm(x, w, cval)
  !-----------------------------------------------------------------------
  ! c = ||x.*w||_2
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL')
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in), target :: x(:)
  type(FVec), intent(in), target :: w(:)
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
!  cval = dsqrt( (x%u * w%u)**2 + (x%v * w%v)**2 + (x%w * w%w)**2 &
!   + (x%phi + w%phi)**2 + (x%theta_dp_cp * w%theta_dp_cp)**2 +   &
!       (x%dp3d * w%dp3d)**2)

  return
end subroutine FNVExtWl2Norm
!=======================================================================



subroutine FNVExtl1Norm(x, cval)
  !-----------------------------------------------------------------------
  ! c = ||x||_1
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL')
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in), target :: x(:)
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
!  cval = dabs(x%u) + dabs(x%v) + dabs(x%w) + dabs(x%phi) + dabs(x%theta_dp_cp) &
!    + dabs(x%dp3d)

  return
end subroutine FNVExtl1Norm
!=======================================================================



subroutine FNVExtCompare(cval, x, z)
  !-----------------------------------------------------------------------
  ! z(i) = 1 if x(i) > c, otherwise z(i) = 0
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL')
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  real*8, intent(in) :: cval
  type(FVec), intent(in), target  :: x(:)
  type(FVec), intent(out), target :: z(:)

  !=======Internals ============

  ! perform vector operation
!  if (x%u > cval)           z%u = 1.d0
!  if (x%v > cval)           z%v = 1.d0
!  if (x%w > cval)           z%w = 1.d0
!  if (x%phi > cval)         z%phi = 1.d0
!  if (x%theta_dp_cp > cval) z%theta_dp_cp = 1.d0
!  if (x%dp3d > cval)        z%dp3d = 1.d0

  return
end subroutine FNVExtCompare
!=======================================================================



subroutine FNVExtInvTest(x, z, cval)
  !-----------------------------------------------------------------------
  ! z = 1./x, if all x nonzero, cval = 0, else cval = 1
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL')
  !-----------------------------------------------------------------------
  use FortranVector
  use iso_c_binding
  implicit none
  type(FVec), intent(in), target  :: x(:)
  type(FVec), intent(out), target :: z(:)
  integer(C_INT), intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
!  cval = 0
!  if (x%u /= 0.d0) then
!     z%u = 1.d0 / x%u
!  else
!     z%u = 0.d0
!     cval = 1
!  endif
!  if (x%v /= 0.d0) then
!     z%v = 1.d0 / x%v
!  else
!     z%v = 0.d0
!     cval = 1
!  endif
!  if (x%w /= 0.d0) then
!     z%w = 1.d0 / x%w
!  else
!     z%w = 0.d0
!     cval = 1
!  endif
 ! if (x%phi /= 0.d0) then
 !    z%phi = 1.d0 / x%phi
 ! else
 !    z%phi = 0.d0
 !    cval = 1
 ! endif
 ! if (x%theta_dp_cp /= 0.d0) then
 !    z%theta_dp_cp = 1.d0 / x%theta_dp_cp
!  else
 !    z%theta_dp_cp = 0.d0
 !    cval = 1
!  endif
!  if (x%dp3d /= 0.d0) then
!     z%dp3d = 1.d0 / x%dp3d
!  else
!     z%dp3d = 0.d0
!     cval = 1
!  endif


  return
end subroutine FNVExtInvTest
!=======================================================================



subroutine FNVExtConstrMask(c, x, m, testval)
  !-----------------------------------------------------------------------
  ! Returns testval = 1 if any element fails constraint test,
  ! otherwise testval = 0.  Test is as follows:
  ! if c(i) =  2, then x(i) >  0
  ! if c(i) =  1, then x(i) >= 0
  ! if c(i) = -1, then x(i) <= 0
  ! if c(i) = -2, then x(i) <  0
  ! fills m(i) = 1 where test fails, m(i) = 0 elsewhere.
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL')
  !-----------------------------------------------------------------------
  use FortranVector
  use iso_c_binding
  implicit none
  type(FVec), intent(in), target  :: c(:)
  type(FVec), intent(in), target  :: x(:)
  type(FVec), intent(out), target :: m(:)
  integer(C_INT), intent(out) :: testval

  !=======Internals ============

  ! perform vector operation
  testval = 0

!  m%u = 0.d0
!  if (dabs(c%u) > 1.5d0) then
!     if (x%u*c%u <= 0.d0) then
!        testval = 1
!        m%u = 1.d0
!     endif
!  else if (dabs(c%u) > 0.5d0) then
!     if (x%u*c%u < 0.d0) then
 !       testval = 1
  !!      m%u = 1.d0
!     endif
!  endif

!  m%v = 0.d0
!  if (dabs(c%v) > 1.5d0) then
!     if (x%v*c%v <= 0.d0) then
!        testval = 1
!        m%v = 1.d0
!     endif
!  else if (dabs(c%v) > 0.5d0) then
!     if (x%v*c%v < 0.d0) then
!        testval = 1
!        m%v = 1.d0
!     endif
!  endif

!  m%w = 0.d0
!  if (dabs(c%w) > 1.5d0) then
!     if (x%w*c%w <= 0.d0) then
!        testval = 1
!        m%w = 1.d0
!     endif
!  else if (dabs(c%w) > 0.5d0) then
!     if (x%w*c%w < 0.d0) then
!        testval = 1
!        m%w = 1.d0
!     endif
!  endif

!  m%phi = 0.d0
!  if (dabs(c%phi) > 1.5d0) then
!     if (x%phi*c%phi <= 0.d0) then
!        testval = 1
!	m%phi = 1.d0
!     endif
!  else if (dabs(c%phi) > 0.5d0) then
!     if (x%u*c%phi < 0.d0) then
!        testval = 1
!        m%phi = 1.d0
!     endif
!  endif

!  m%theta_dp_cp = 0.d0
!  if (dabs(c%theta_dp_cp) > 1.5d0) then
 !    if (x%theta_dp_cp*c%theta_dp_cp <= 0.d0) then
!        testval = 1
!        m%theta_dp_cp = 1.d0
!     endif
!  else if (dabs(c%theta_dp_cp) > 0.5d0) then
!     if (x%theta_dp_cp*c%theta_dp_cp < 0.d0) then
!        testval = 1
!        m%theta_dp_cp = 1.d0
!     endif
!  endif

!  m%dp3d = 0.d0
!  if (dabs(c%dp3d) > 1.5d0) then
!     if (x%dp3d*c%dp3d <= 0.d0) then
!     testval = 1
!        m%dp3d = 1.d0
!     endif
!  else if (dabs(c%dp3d) > 0.5d0) then
!     if (x%dp3d*c%dp3d < 0.d0) then
!     testval = 1
!        m%dp3d = 1.d0
!     endif
 ! endif


  return
end subroutine FNVExtConstrMask
!=======================================================================



subroutine FNVExtMinQuotient(x, y, cval)
  !-----------------------------------------------------------------------
  ! c = min(x./y), over all y /= 0
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL')
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in), target :: x(:)
  type(FVec), intent(in), target :: y(:)
  real*8, intent(out) :: cval
  logical :: notyet

  !=======Internals ============

  ! perform vector operation
  cval = 0.d0
!  notyet = .true.
!  if (y%u /= 0.d0) then
!     if (notyet) then
!        cval = x%u / y%u
!        notyet = .false.
!     else
!        cval = min( cval, x%u / y%u )
!     endif
!  endif
!  if (y%v /= 0.d0) then
!     if (notyet) then
!        cval = x%v / y%v
!        notyet = .false.
!     else
!        cval = min( cval, x%v / y%v )
!     endif
!  endif
!  if (y%w /= 0.d0) then
!     if (notyet) then
!       cval = x%w / y%w
!        notyet = .false.
!     else
!        cval = min( cval, x%w / y%w )
!     endif
!  endif
!  if (y%phi /= 0.d0) then
!     if (notyet) then
!        cval = x%phi / y%phi
!        notyet = .false.
!     else
!        cval = min( cval, x%phi / y%phi )
!     endif
!  endif
!  if (y%theta_dp_cp /= 0.d0) then
!     if (notyet) then
!        cval = x%theta_dp_cp / y%theta_dp_cp
!        notyet = .false.
!     else
!        cval = min( cval, x%theta_dp_cp / y%theta_dp_cp )
!     endif
!  endif
!!  if (y%dp3d /= 0.d0) then
!    if (notyet) then
!        cval = x%dp3d / y%dp3d
!        notyet = .false.
!     else
!        cval = min( cval, x%dp3d / y%dp3d )
!     endif
!  endif


  return
end subroutine FNVExtMinQuotient
!=======================================================================

#endif
