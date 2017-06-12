!-----------------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! Copyright 2017; all rights reserved
!=======================================================================


module HommeNVector
  !-----------------------------------------------------------------------
  ! Description: simple Fortran user-defined type for example interface
  !-----------------------------------------------------------------------
  use element_mod,    only: element_t
  use hybvcoord_mod,  only: hvcoord_t
  use hybrid_mod,     only: hybrid_t
  use derivative_mod, only: derivative_t
  use element_state,  only: timelevels
  use, intrinsic :: iso_c_binding
  public

  type :: NVec_t
     type(element_t),    pointer :: elem(:)
     type(hvcoord_t),    pointer :: hvcoord
     type(hybrid_t),     pointer :: hybrid
     type(derivative_t), pointer :: deriv
     integer                     :: nets
     integer                     :: nete
     integer                     :: qn0
     integer                     :: tl_idx
  end type NVec_t

  save

  integer, parameter :: RegistryLength=timelevels
  logical :: HommeNVectorRegistry(RegistryLength)=.false.

  !-------------------------------------

contains
  
  integer function ReserveHommeNVectorRegistryIdx()
    ! Marks first available vector from registry as used (if any are available)
    integer :: i
    do i=1,RegistryLength
       if (.not. HommeNVectorRegistry(i)) then
          HommeNVectorRegistry(i) = .true.
          ReserveHommeNVectorRegistryIdx = i
          return
       end if
    end do
    ReserveHommeNVectorRegistryIdx = 0
    return
  end function ReserveHommeNVectorRegistryIdx

  subroutine MakeHommeNVector(elem, hvcoord, hybrid, deriv, nets, nete, qn0, tl_idx, vec, ier)
    ! Function to create an NVec_t wrapper (vec) to the Homme solution 
    !   structure, and reserve its index in the vector registry.
    ! The return value is 0 if succssful, 1 if failure (if the registry 
    !   index is already taken)
    use element_mod,    only: element_t
    use hybvcoord_mod,  only: hvcoord_t
    use hybrid_mod,     only: hybrid_t
    use derivative_mod, only: derivative_t
    type (element_t),    target    :: elem(:)
    type (hvcoord_t),    target    :: hvcoord
    type (hybrid_t),     target    :: hybrid
    type (derivative_t), target    :: deriv
    integer, intent(in)            :: nets
    integer, intent(in)            :: nete
    integer, intent(in)            :: qn0
    integer, intent(in)            :: tl_idx
    type(NVec_t), intent(out)      :: vec
    integer, intent(out)           :: ier

    ! if registry index is already taken, return failure
    ier = 0
    if (HommeNVectorRegistry(tl_idx)) then
       ier = 1
       return
    end if
    
    ! otherwise, set up wrapper and reserve registry entry
    vec%elem => elem
    vec%hvcoord => hvcoord
    vec%hybrid => hybrid
    vec%deriv => deriv
    vec%nets = nets
    vec%nete = nete
    vec%qn0 = qn0
    vec%tl_idx = tl_idx
    HommeNVectorRegistry(tl_idx) = .true.
    return
  end subroutine MakeHommeNVector

  !=======================================================================
end module HommeNVector



subroutine FNVExtPrint(x_C)
  !-----------------------------------------------------------------------
  ! Print routine for EXT vector
  ! 
  ! Note: this function is not required by ARKode (or any of SUNDIALS) -- 
  !       it is merely here for convenience when debugging
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr) :: x_C
  type(NVec_t), pointer :: x => NULL()

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! print vector data (not implemented)

  return
end subroutine FNVExtPrint
!=======================================================================



subroutine FNVExtClone(x_C, y_C)
  !-----------------------------------------------------------------------
  ! Clone routine for EXT vector -- allocates memory for a new NVec_t y that 
  ! that matches structure for x.
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t, ReserveHommeNVectorRegistryIdx
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr) :: x_C
  type(c_ptr) :: y_C
  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  integer :: tl_idx

  !=======Internals ============

  ! grab next available registry index; 
  ! if none are available, return with a NULL value for y_C 
  tl_idx = ReserveHommeNVectorRegistryIdx()
  if (tl_idx == 0) then  
     y_C = C_NULL_PTR
     return
  end if

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! allocate new vector
  allocate(y)

  ! setup y to replicate x, but with unique registry index
  y%elem => x%elem
  y%hvcoord => x%hvcoord
  y%hybrid => x%hybrid
  y%deriv => x%deriv
  y%nets = x%nets
  y%nete = x%nete
  y%qn0 = x%qn0
  y%tl_idx = tl_idx
  
  ! associate y_C pointer with y
  y_C = c_loc(y)

  return
end subroutine FNVExtClone
!=======================================================================



subroutine FNVExtDestroy(x_C)
  !-----------------------------------------------------------------------
  ! Destroy routine for EXT vector -- nullifies data for a given NVec_t
  ! and returns its vector index to the registry
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t, HommeNVectorRegistry
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr) :: x_C
  type(NVec_t), pointer :: x => NULL()
  integer :: ierr

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! return index to registry
  HommeNVectorRegistry(x%tl_idx) = .false.
  
  ! detach pointers, deallocate vector and nullify x_C
  if (associated(x%elem))     nullify(x%elem)
  if (associated(x%hvcoord))  nullify(x%hvcoord)
  if (associated(x%hybrid))   nullify(x%hybrid)
  if (associated(x%deriv))    nullify(x%deriv)
  deallocate(x)
  x_C = C_NULL_PTR

  return
end subroutine FNVExtDestroy
!=======================================================================



subroutine FNVExtLinearSum(aval, x_C, bval, y_C, z_C)
  !-----------------------------------------------------------------------
  ! z = a*x + b*y
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: aval
  type(c_ptr),    intent(in) :: x_C
  real(c_double), intent(in) :: bval
  type(c_ptr),    intent(in) :: y_C
  type(c_ptr),    intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  type(NVec_t), pointer :: z => NULL()
  integer :: ie

  !=======Internals ============

  ! dereference pointer for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=x%nets,x%nete
     z%elem(ie)%state%v(:,:,:,:,x%tl_idx) = &
          aval * z%elem(ie)%state%v(:,:,:,:,z%tl_idx) + &
          bval * y%elem(ie)%state%v(:,:,:,:,y%tl_idx)
     z%elem(ie)%state%w(:,:,:,x%tl_idx) = &
          aval * z%elem(ie)%state%w(:,:,:,z%tl_idx) + &
          bval * y%elem(ie)%state%w(:,:,:,y%tl_idx)
     z%elem(ie)%state%phi(:,:,:,x%tl_idx) = &
          aval * z%elem(ie)%state%phi(:,:,:,z%tl_idx) + &
          bval * y%elem(ie)%state%phi(:,:,:,y%tl_idx)
     z%elem(ie)%state%theta_dp_cp(:,:,:,x%tl_idx) = &
          aval * z%elem(ie)%state%theta_dp_cp(:,:,:,z%tl_idx) + &
          bval * y%elem(ie)%state%theta_dp_cp(:,:,:,y%tl_idx)
     z%elem(ie)%state%dp3d(:,:,:,x%tl_idx) = &
          aval * z%elem(ie)%state%dp3d(:,:,:,z%tl_idx) + &
          bval * y%elem(ie)%state%dp3d(:,:,:,y%tl_idx)
     !z%elem(ie)%state%psv(:,:,x%tl_idx) = &
     !     aval * z%elem(ie)%state%psv(:,:,z%tl_idx) + &
     !     bval * y%elem(ie)%state%psv(:,:,y%tl_idx)
  end do

  return
end subroutine FNVExtLinearSum
!=======================================================================



subroutine FNVExtConst(cval, z_C)
  !-----------------------------------------------------------------------
  ! z = c 
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: cval
  type(c_ptr),    intent(in) :: z_C

  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference z_C pointer for NVec_t object
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtConst
!=======================================================================



subroutine FNVExtProd(x_C, y_C, z_C)
  !-----------------------------------------------------------------------
  ! z = x.*y (Matlab notation)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: y_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtProd
!=======================================================================



subroutine FNVExtDiv(x_C, y_C, z_C)
  !-----------------------------------------------------------------------
  ! z = x./y (Matlab notation; doesn't check for legal denominator)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: y_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtDiv
!=======================================================================



subroutine FNVExtScale(cval, x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = c*x
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: cval
  type(c_ptr),    intent(in) :: x_C
  type(c_ptr),    intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtScale
!=======================================================================



subroutine FNVExtAbs(x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = |x|  (componentwise)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtAbs
!=======================================================================



subroutine FNVExtInv(x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = 1./x (Matlab notation: doesn't check for division by zero)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtInv
!=======================================================================



subroutine FNVExtAddConst(cval, x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = x+c (add c to each component of the vector x)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: cval
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtAddConst
!=======================================================================



subroutine FNVExtDotProd(x_C, y_C, cval)
  !-----------------------------------------------------------------------
  ! c = <x,y>  (only include 'active' data; no ghost cells, etc.)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: y_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtDotProd
!=======================================================================



subroutine FNVExtMaxNorm(x_C, cval)
  !-----------------------------------------------------------------------
  ! c = max(|x|)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtMaxNorm
!=======================================================================



subroutine FNVExtWrmsNorm(x_C, w_C, cval)
  !-----------------------------------------------------------------------
  ! cval = sqrt(||x.*w||^2_2 / Ntotal)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: w_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: w => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(w_C, w)

  ! perform vector operation
  cval = 0.d0
  
  return
end subroutine FNVExtWrmsNorm
!=======================================================================



subroutine FNVExtWrmsNormMask(x_C, w_C, l_C, cval)
  !-----------------------------------------------------------------------
  ! cval = ||x.*w||_2 / Ntotal (but only where l > 0)
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL()')
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: w_C
  type(c_ptr),    intent(in)  :: l_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: w => NULL()
  type(NVec_t), pointer :: l => NULL()
  integer :: ntot

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(w_C, w)
  call c_f_pointer(l_C, l)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtWrmsNormMask
!=======================================================================



subroutine FNVExtMin(x_C, cval)
  !-----------------------------------------------------------------------
  ! cval = min(xvec)
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtMin
!=======================================================================



subroutine FNVExtWl2Norm(x_C, w_C, cval)
  !-----------------------------------------------------------------------
  ! c = ||x.*w||_2
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL()')
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: w_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: w => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(w_C, w)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtWl2Norm
!=======================================================================



subroutine FNVExtl1Norm(x_C, cval)
  !-----------------------------------------------------------------------
  ! c = ||x||_1
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL()')
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtl1Norm
!=======================================================================



subroutine FNVExtCompare(cval, x_C, z_C)
  !-----------------------------------------------------------------------
  ! z(i) = 1 if x(i) > c, otherwise z(i) = 0
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL()')
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: cval
  type(c_ptr),    intent(in) :: x_C
  type(c_ptr),    intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation

  return
end subroutine FNVExtCompare
!=======================================================================



subroutine FNVExtInvTest(x_C, z_C, cval)
  !-----------------------------------------------------------------------
  ! z = 1./x, if all x nonzero, cval = 0, else cval = 1
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL()')
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: z_C
  integer(c_int), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtInvTest
!=======================================================================



subroutine FNVExtConstrMask(c_C, x_C, m_C, testval)
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
  !       (set the corresponding function pointer in the 'ops' to 'NULL()')
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: c_C
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: m_C
  integer(c_int), intent(out) :: testval

  type(NVec_t), pointer :: c => NULL()
  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: m => NULL()

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(c_C, c)
  call c_f_pointer(x_C, x)
  call c_f_pointer(m_C, m)

  ! perform vector operation
  testval = 0.d0

  return
end subroutine FNVExtConstrMask
!=======================================================================



subroutine FNVExtMinQuotient(x_C, y_C, cval)
  !-----------------------------------------------------------------------
  ! c = min(x./y), over all y /= 0
  ! 
  ! Note: this function is not required by ARKode (it is used for other
  !       SUNDIALS solvers).  As such it need not be implemented, and 
  !       it may be removed from the associated nvector_external files
  !       (set the corresponding function pointer in the 'ops' to 'NULL()')
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: y_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  logical :: notyet

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtMinQuotient
!=======================================================================
