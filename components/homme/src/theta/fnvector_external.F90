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

  type FVec
     real*8 :: u
     real*8 :: v
     real*8 :: w
  end type FVec
  !=======================================================================
end module FortranVector



subroutine FNVExtPrint(x)
  !-----------------------------------------------------------------------
  ! Print routine for EXT vector
  ! 
  ! Note: this function is not required by ARKode (or any of SUNDIALS) -- 
  !       it is merely here for convenience when debugging
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x

  !=======Internals ============

  ! print vector data
  print '(4x,3(es17.10,1x))', x%u, x%v, x%w

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
  type(FVec), intent(in) :: x
  type(C_PTR), intent(out) :: y_C
  type(FVec), pointer :: y

  !=======Internals ============

  ! allocate new vector (in this demonstration code, there's no need to 
  ! actually look at x since we know the structure of these vectors)
  allocate(y)

  ! initialize values to zero
  y%u = 0.d0
  y%v = 0.d0
  y%w = 0.d0

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
  type(FVec), pointer :: x
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
  real*8, intent(in) :: aval
  type(FVec), intent(in) :: x
  real*8, intent(in) :: bval
  type(FVec), intent(in) :: y
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = aval*x%u + bval*y%u
  z%v = aval*x%v + bval*y%v
  z%w = aval*x%w + bval*y%w

  return
end subroutine FNVExtLinearSum
!=======================================================================



subroutine FNVExtConst(cval, z)
  !-----------------------------------------------------------------------
  ! z = c 
  !-----------------------------------------------------------------------
  use FortranVector
  use iso_c_binding
  implicit none
  real*8, intent(in) :: cval
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = cval
  z%v = cval
  z%w = cval

  return
end subroutine FNVExtConst
!=======================================================================



subroutine FNVExtProd(x, y, z)
  !-----------------------------------------------------------------------
  ! z = x.*y (Matlab notation)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x
  type(FVec), intent(in) :: y
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = x%u * y%u
  z%v = x%v * y%v
  z%w = x%w * y%w

  return
end subroutine FNVExtProd
!=======================================================================



subroutine FNVExtDiv(x, y, z)
  !-----------------------------------------------------------------------
  ! z = x./y (Matlab notation; doesn't check for legal denominator)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x
  type(FVec), intent(in) :: y
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = x%u / y%u
  z%v = x%v / y%v
  z%w = x%w / y%w

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
  real*8, intent(in)  :: cval
  type(FVec), intent(in) :: x
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = cval * x%u
  z%v = cval * x%v
  z%w = cval * x%w

  return
end subroutine FNVExtScale
!=======================================================================



subroutine FNVExtAbs(x, z)
  !-----------------------------------------------------------------------
  ! z = |x|  (componentwise)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = dabs(x%u)
  z%v = dabs(x%v)
  z%w = dabs(x%w)

  return
end subroutine FNVExtAbs
!=======================================================================



subroutine FNVExtInv(x, z)
  !-----------------------------------------------------------------------
  ! z = 1./x (Matlab notation: doesn't check for division by zero)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = 1.d0 / x%u
  z%v = 1.d0 / x%v
  z%w = 1.d0 / x%w

  return
end subroutine FNVExtInv
!=======================================================================



subroutine FNVExtAddConst(cval, x, z)
  !-----------------------------------------------------------------------
  ! z = x+c (add c to each component of the vector x)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  real*8, intent(in) :: cval
  type(FVec), intent(in) :: x
  type(FVec), intent(inout) :: z

  !=======Internals ============

  ! perform vector operation
  z%u = x%u + cval
  z%v = x%v + cval
  z%w = x%w + cval

  return
end subroutine FNVExtAddConst
!=======================================================================



subroutine FNVExtDotProd(x, y, cval)
  !-----------------------------------------------------------------------
  ! c = <x,y>  (only include 'active' data; no ghost cells, etc.)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x
  type(FVec), intent(in) :: y
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
  cval  =  x%u * y%u  +  x%v * y%v  +  x%w * y%w

  return
end subroutine FNVExtDotProd
!=======================================================================



subroutine FNVExtMaxNorm(x, cval)
  !-----------------------------------------------------------------------
  ! c = max(|x|)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
  cval = max( dabs(x%u), max( dabs(x%v), dabs(x%w)) )

  return
end subroutine FNVExtMaxNorm
!=======================================================================



subroutine FNVExtWrmsNorm(x, w, cval)
  !-----------------------------------------------------------------------
  ! cval = sqrt(||x.*w||^2_2 / Ntotal)
  !-----------------------------------------------------------------------
  use FortranVector
  implicit none
  type(FVec), intent(in) :: x
  type(FVec), intent(in) :: w
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
  cval = dsqrt( ((x%u * w%u)**2 + (x%v * w%v)**2 + (x%w * w%w)**2) / 3.d0 )
  
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
  type(FVec), intent(in) :: x
  type(FVec), intent(in) :: w
  type(FVec), intent(in) :: l
  real*8, intent(out) :: cval
  integer :: ntot

  !=======Internals ============

  ! perform vector operation
  cval = 0.d0
  ntot = 0
  if (l%u > 0.d0) then
     cval = cval + (x%u * w%u)**2
     ntot = ntot + 1
  endif
  if (l%v > 0.d0) then
     cval = cval + (x%v * w%v)**2
     ntot = ntot + 1
  endif
  if (l%w > 0.d0) then
     cval = cval + (x%w * w%w)**2
     ntot = ntot + 1
  endif
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
  type(FVec), intent(in) :: x
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
  cval = min( x%u, min( x%v, x%w) )

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
  type(FVec), intent(in) :: x
  type(FVec), intent(in) :: w
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
  cval = dsqrt( (x%u * w%u)**2 + (x%v * w%v)**2 + (x%w * w%w)**2 )

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
  type(FVec), intent(in) :: x
  real*8, intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
  cval = dabs(x%u) + dabs(x%v) + dabs(x%w)

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
  type(FVec), intent(in)  :: x
  type(FVec), intent(out) :: z

  !=======Internals ============

  ! perform vector operation
  if (x%u > cval)  z%u = 1.d0
  if (x%v > cval)  z%v = 1.d0
  if (x%w > cval)  z%w = 1.d0

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
  type(FVec), intent(in)  :: x
  type(FVec), intent(out) :: z
  integer(C_INT), intent(out) :: cval

  !=======Internals ============

  ! perform vector operation
  cval = 0
  if (x%u /= 0.d0) then
     z%u = 1.d0 / x%u
  else
     z%u = 0.d0
     cval = 1
  endif
  if (x%v /= 0.d0) then
     z%v = 1.d0 / x%v
  else
     z%v = 0.d0
     cval = 1
  endif
  if (x%w /= 0.d0) then
     z%w = 1.d0 / x%w
  else
     z%w = 0.d0
     cval = 1
  endif

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
  type(FVec), intent(in)  :: c
  type(FVec), intent(in)  :: x
  type(FVec), intent(out) :: m
  integer(C_INT), intent(out) :: testval

  !=======Internals ============

  ! perform vector operation
  testval = 0

  m%u = 0.d0
  if (dabs(c%u) > 1.5d0) then
     if (x%u*c%u <= 0.d0) then
        testval = 1
        m%u = 1.d0
     endif
  else if (dabs(c%u) > 0.5d0) then
     if (x%u*c%u < 0.d0) then
        testval = 1
        m%u = 1.d0
     endif
  endif

  m%v = 0.d0
  if (dabs(c%v) > 1.5d0) then
     if (x%v*c%v <= 0.d0) then
        testval = 1
        m%v = 1.d0
     endif
  else if (dabs(c%v) > 0.5d0) then
     if (x%v*c%v < 0.d0) then
        testval = 1
        m%v = 1.d0
     endif
  endif

  m%w = 0.d0
  if (dabs(c%w) > 1.5d0) then
     if (x%w*c%w <= 0.d0) then
        testval = 1
        m%w = 1.d0
     endif
  else if (dabs(c%w) > 0.5d0) then
     if (x%w*c%w < 0.d0) then
        testval = 1
        m%w = 1.d0
     endif
  endif

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
  type(FVec), intent(in) :: x
  type(FVec), intent(in) :: y
  real*8, intent(out) :: cval
  logical :: notyet

  !=======Internals ============

  ! perform vector operation
  cval = 0.d0
  notyet = .true.
  if (y%u /= 0.d0) then
     if (notyet) then
        cval = x%u / y%u
        notyet = .false.
     else
        cval = min( cval, x%u / y%u )
     endif
  endif
  if (y%v /= 0.d0) then
     if (notyet) then
        cval = x%v / y%v
        notyet = .false.
     else
        cval = min( cval, x%v / y%v )
     endif
  endif
  if (y%w /= 0.d0) then
     if (notyet) then
        cval = x%w / y%w
        notyet = .false.
     else
        cval = min( cval, x%w / y%w )
     endif
  endif

  return
end subroutine FNVExtMinQuotient
!=======================================================================
