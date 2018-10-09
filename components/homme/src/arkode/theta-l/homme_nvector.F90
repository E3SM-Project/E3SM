!-----------------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! Copyright 2017; all rights reserved
!
! Modified by Christopher J. Vogl (LLNL) 2018
!=======================================================================

module HommeNVector
  !-----------------------------------------------------------------------
  ! Description: simple Fortran user-defined type for example interface
  !-----------------------------------------------------------------------
  use element_mod,    only: element_t
  use element_state,  only: timelevels
  use parallel_mod,   only: parallel_t
  use, intrinsic :: iso_c_binding
  public

  type :: NVec_t
     type(element_t), pointer :: elem(:)
     integer                  :: nets
     integer                  :: nete
     integer                  :: tl_idx
  end type NVec_t

  save

  integer, parameter :: RegistryLength=timelevels
  logical :: HommeNVectorRegistry(RegistryLength)=.false.
  type(parallel_t), pointer :: par_ptr

  !-------------------------------------

contains

  subroutine SetHommeNVectorPar(par)
    ! sets current par object pointer for MPI_Allreduce and wrap_repro_sum calls
    implicit none
    type(parallel_t), target, intent(in) :: par
    par_ptr => par
  end subroutine

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

  subroutine MakeHommeNVector(elem, nets, nete, tl_idx, vec, ierr)
    ! Function to create an NVec_t wrapper (vec) to the Homme solution
    !   structure, and reserve its index in the vector registry.
    ! The return value is 0 if succssful, 1 if failure (if the registry
    !   index is already taken)
    use element_mod,    only: element_t
    type (element_t),    target    :: elem(:)
    integer, intent(in)            :: nets
    integer, intent(in)            :: nete
    integer, intent(in)            :: tl_idx
    type(NVec_t), intent(out)      :: vec
    integer(C_INT), intent(out)    :: ierr

    ! if registry index is already taken, return failure
    ierr = 0
    if (HommeNVectorRegistry(tl_idx)) then
       ierr = 1
       return
    end if

    ! otherwise, set up wrapper and reserve registry entry
    vec%elem => elem
    vec%nets = nets
    vec%nete = nete
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
  use HommeNVector,     only: NVec_t, par_ptr
  use dimensions_mod,   only: np, nlev, nlevp
  use MPI,              only: MPI_comm_rank

  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr) :: x_C
  type(NVec_t), pointer :: x => NULL()

  integer :: ie, inpx, inpy, inlev, rank, ierr

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! get rank
  call MPI_comm_rank(par_ptr%comm,rank,ierr)

  ! print vector data
  do inlev=1,nlev
    do ie=x%nets,x%nete
      do inpx=1,np
        do inpy=1,np
          print '(/,"proc ",i2,",", " elem ",i4,",", " u(",i1,",",i1,",",i2,") = ",f15.5)', &
            rank, ie, inpx, inpy, inlev, x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)
          print '("proc ",i2,",", " elem ",i4,",", " v(",i1,",",i1,",",i2,") = ",f15.5)', &
            rank, ie, inpx, inpy, inlev, x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)
          print '("proc ",i2,",", " elem ",i4,",", " vtheta_dp(",i1,",",i1,",",i2,") = ",f15.5)', &
            rank, ie, inpx, inpy, inlev, x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)
          print '("proc ",i2,",", " elem ",i4,",", " dp3d(",i1,",",i1,",",i2,") = ",f15.5,/)', &
            rank, ie, inpx, inpy, inlev, x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)
        end do
      end do
    end do
  end do
  do inlev=1,nlevp
    do ie=x%nets,x%nete
      do inpx=1,np
        do inpy=1,np
          print '("proc ",i2,",", " elem ",i4,",", " w(",i1,",",i1,",",i2,") = ",f15.5)', &
            rank, ie, inpx, inpy, inlev, x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)
          print '("proc ",i2,",", " elem ",i4,",", " phinh(",i1,",",i1,",",i2,") = ",f15.5)', &
            rank, ie, inpx, inpy, inlev, x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)
        end do
      end do
    end do
  end do

  return
end subroutine FNVExtPrint
!=======================================================================



subroutine FNVExtClone(x_C, y_C, ierr)
  !-----------------------------------------------------------------------
  ! Clone routine for EXT vector -- allocates memory for a new NVec_t y that
  ! that matches structure for x.
  !-----------------------------------------------------------------------
  use HommeNVector, only: NVec_t, ReserveHommeNVectorRegistryIdx
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(out) :: y_C
  integer(c_int), intent(out) :: ierr
  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  integer :: tl_idx

  !=======Internals ============

  ! grab next available registry index;
  ! if none are available, set error flag and return
  ierr = 0
  tl_idx = ReserveHommeNVectorRegistryIdx()
  if (tl_idx == 0) then
     ierr = 1
     return
  end if

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! allocate new vector
  allocate(y)

  ! setup y to replicate x, but with unique registry index
  y%elem => x%elem
  y%nets = x%nets
  y%nete = x%nete
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

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! return index to registry
  HommeNVectorRegistry(x%tl_idx) = .false.

  ! detach pointers, deallocate vector and nullify x_C
  if (associated(x%elem))     nullify(x%elem)
  deallocate(x)
  x_C = C_NULL_PTR

  return
end subroutine FNVExtDestroy
!=======================================================================



subroutine FNVExtLinearSum(aval, x_C, bval, y_C, z_C)
  !-----------------------------------------------------------------------
  ! z = a*x + b*y
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev, nlevp
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

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointer for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%v(inpx,inpy,1,inlev,z%tl_idx) = &
            aval * x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx) + &
            bval * y%elem(ie)%state%v(inpx,inpy,1,inlev,y%tl_idx)
          z%elem(ie)%state%v(inpx,inpy,2,inlev,z%tl_idx) = &
            aval * x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx) + &
            bval * y%elem(ie)%state%v(inpx,inpy,2,inlev,y%tl_idx)
          z%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,z%tl_idx) = &
            aval * x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx) + &
            bval * y%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,y%tl_idx)
          z%elem(ie)%state%dp3d(inpx,inpy,inlev,z%tl_idx) = &
            aval * x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx) + &
            bval * y%elem(ie)%state%dp3d(inpx,inpy,inlev,y%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%w_i(inpx,inpy,inlev,z%tl_idx) = &
            aval * x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx) + &
            bval * y%elem(ie)%state%w_i(inpx,inpy,inlev,y%tl_idx)
          z%elem(ie)%state%phinh_i(inpx,inpy,inlev,z%tl_idx) = &
            aval * x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx) + &
            bval * y%elem(ie)%state%phinh_i(inpx,inpy,inlev,y%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  return
end subroutine FNVExtLinearSum
!=======================================================================



subroutine FNVExtConst(cval, z_C)
  !-----------------------------------------------------------------------
  ! z = c
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: cval
  type(c_ptr),    intent(in) :: z_C

  type(NVec_t), pointer :: z => NULL()

  integer :: ie

  !=======Internals ============

  ! dereference z_C pointer for NVec_t object
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
     z%elem(ie)%state%v(:,:,:,:,z%tl_idx) = cval
     z%elem(ie)%state%w_i(:,:,:,z%tl_idx) = cval
     z%elem(ie)%state%phinh_i(:,:,:,z%tl_idx) = cval
     z%elem(ie)%state%vtheta_dp(:,:,:,z%tl_idx) = cval
     z%elem(ie)%state%dp3d(:,:,:,z%tl_idx) = cval
  end do

  return
end subroutine FNVExtConst
!=======================================================================



subroutine FNVExtProd(x_C, y_C, z_C)
  !-----------------------------------------------------------------------
  ! z = x.*y (Matlab notation)
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev, nlevp
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: y_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  type(NVec_t), pointer :: z => NULL()

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%v(inpx,inpy,1,inlev,z%tl_idx) = &
            x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)* &
            y%elem(ie)%state%v(inpx,inpy,1,inlev,y%tl_idx)
          z%elem(ie)%state%v(inpx,inpy,2,inlev,z%tl_idx) = &
            x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)* &
            y%elem(ie)%state%v(inpx,inpy,2,inlev,y%tl_idx)
          z%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)* &
            y%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,y%tl_idx)
          z%elem(ie)%state%dp3d(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)* &
            y%elem(ie)%state%dp3d(inpx,inpy,inlev,y%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%w_i(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)* &
            y%elem(ie)%state%w_i(inpx,inpy,inlev,y%tl_idx)
          z%elem(ie)%state%phinh_i(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)* &
            y%elem(ie)%state%phinh_i(inpx,inpy,inlev,y%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  return
end subroutine FNVExtProd
!=======================================================================



subroutine FNVExtDiv(x_C, y_C, z_C)
  !-----------------------------------------------------------------------
  ! z = x./y (Matlab notation; doesn't check for legal denominator)
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev, nlevp
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: y_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()
  type(NVec_t), pointer :: z => NULL()

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%v(inpx,inpy,1,inlev,z%tl_idx) = &
            x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)/ &
            y%elem(ie)%state%v(inpx,inpy,1,inlev,y%tl_idx)
          z%elem(ie)%state%v(inpx,inpy,2,inlev,z%tl_idx) = &
            x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)/ &
            y%elem(ie)%state%v(inpx,inpy,2,inlev,y%tl_idx)
          z%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)/ &
            y%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,y%tl_idx)
          z%elem(ie)%state%dp3d(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)/ &
            y%elem(ie)%state%dp3d(inpx,inpy,inlev,y%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%w_i(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)/ &
            y%elem(ie)%state%w_i(inpx,inpy,inlev,y%tl_idx)
          z%elem(ie)%state%phinh_i(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)/ &
            y%elem(ie)%state%phinh_i(inpx,inpy,inlev,y%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  return
end subroutine FNVExtDiv
!=======================================================================



subroutine FNVExtScale(cval, x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = c*x
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev, nlevp
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: cval
  type(c_ptr),    intent(in) :: x_C
  type(c_ptr),    intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  integer :: ie, inlev, inpy, inpx

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%v(inpx,inpy,1,inlev,z%tl_idx) = &
            cval * x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)
          z%elem(ie)%state%v(inpx,inpy,2,inlev,z%tl_idx) = &
            cval * x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)
          z%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,z%tl_idx) = &
            cval * x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)
          z%elem(ie)%state%dp3d(inpx,inpy,inlev,z%tl_idx) = &
            cval * x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%w_i(inpx,inpy,inlev,z%tl_idx) = &
            cval * x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)
          z%elem(ie)%state%phinh_i(inpx,inpy,inlev,z%tl_idx) = &
            cval * x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  return
end subroutine FNVExtScale
!=======================================================================



subroutine FNVExtAbs(x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = |x|  (componentwise)
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev, nlevp
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%v(inpx,inpy,1,inlev,z%tl_idx) = &
            abs(x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx))
          z%elem(ie)%state%v(inpx,inpy,2,inlev,z%tl_idx) = &
            abs(x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx))
          z%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,z%tl_idx) = &
            abs(x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx))
          z%elem(ie)%state%dp3d(inpx,inpy,inlev,z%tl_idx) = &
            abs(x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx))
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%w_i(inpx,inpy,inlev,z%tl_idx) = &
            abs(x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx))
          z%elem(ie)%state%phinh_i(inpx,inpy,inlev,z%tl_idx) = &
            abs(x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx))
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  return
end subroutine FNVExtAbs
!=======================================================================



subroutine FNVExtInv(x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = 1./x (Matlab notation: doesn't check for division by zero)
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev, nlevp
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%v(inpx,inpy,1,inlev,z%tl_idx) = &
            1.d0 / x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)
          z%elem(ie)%state%v(inpx,inpy,2,inlev,z%tl_idx) = &
            1.d0 / x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)
          z%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,z%tl_idx) = &
            1.d0 / x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)
          z%elem(ie)%state%dp3d(inpx,inpy,inlev,z%tl_idx) = &
            1.d0 / x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%w_i(inpx,inpy,inlev,z%tl_idx) = &
            1.d0 / x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)
          z%elem(ie)%state%phinh_i(inpx,inpy,inlev,z%tl_idx) = &
            1.d0 / x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  return
end subroutine FNVExtInv
!=======================================================================


subroutine FNVExtAddConst(cval, x_C, z_C)
  !-----------------------------------------------------------------------
  ! z = x+c (add c to each component of the vector x)
  !-----------------------------------------------------------------------
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev, nlevp
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double), intent(in) :: cval
  type(c_ptr), intent(in) :: x_C
  type(c_ptr), intent(in) :: z_C

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: z => NULL()

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(z_C, z)

  ! perform vector operation
  do ie=z%nets,z%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%v(inpx,inpy,1,inlev,z%tl_idx) = &
            x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx) + cval
          z%elem(ie)%state%v(inpx,inpy,2,inlev,z%tl_idx) = &
            x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx) + cval
          z%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx) + cval
          z%elem(ie)%state%dp3d(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx) + cval
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          z%elem(ie)%state%w_i(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx) + cval
          z%elem(ie)%state%phinh_i(inpx,inpy,inlev,z%tl_idx) = &
            x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx) + cval
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  return
end subroutine FNVExtAddConst
!=======================================================================



subroutine FNVExtDotProd(x_C, y_C, cval)
  !-----------------------------------------------------------------------
  ! c = <x,y>  (only include 'active' data; no ghost cells, etc.)
  !-----------------------------------------------------------------------
  use HommeNVector,     only: NVec_t, par_ptr
  use dimensions_mod,   only: np, nlev, nlevp
  use parallel_mod,     only: global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: y_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: y => NULL()

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)

  ! NOTE: this use of spheremp acknowledges the fact that unknowns are
  ! duplicated, in that the sum including spheremp performs the
  ! integral over the domain
  do ie=x%nets,x%nete
    global_shared_buf(ie,1) = 0.d0
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          global_shared_buf(ie,1) = global_shared_buf(ie,1) + &
            (x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)* &
              y%elem(ie)%state%v(inpx,inpy,1,inlev,y%tl_idx) + &
            x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)* &
              y%elem(ie)%state%v(inpx,inpy,2,inlev,y%tl_idx) + &
            x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)* &
              y%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,y%tl_idx) + &
            x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)* &
              y%elem(ie)%state%dp3d(inpx,inpy,inlev,y%tl_idx)) &
            * x%elem(ie)%spheremp(inpx,inpy)
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          global_shared_buf(ie,1) = global_shared_buf(ie,1) + &
            (x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)* &
              y%elem(ie)%state%w_i(inpx,inpy,inlev,y%tl_idx) + &
            x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)* &
              y%elem(ie)%state%phinh_i(inpx,inpy,inlev,y%tl_idx)) &
            * x%elem(ie)%spheremp(inpx,inpy)
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  ! accumulate sum using wrap_repro_sum and then copy to cval
  ! Q: should we divide by 4*pi to account for the integral?
  call wrap_repro_sum(nvars=1, comm=par_ptr%comm)
  cval = global_shared_sum(1)

  return
end subroutine FNVExtDotProd
!=======================================================================



subroutine FNVExtMaxNorm(x_C, cval)
  !-----------------------------------------------------------------------
  ! c = max(|x|)
  !-----------------------------------------------------------------------
  use HommeNVector,     only: NVec_t, par_ptr
  use dimensions_mod,   only: np, nlev, nlevp
  use MPI,              only: MPI_MAX
  use parallel_mod,     only: MPIreal_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()

  integer :: ie, inlev, inpx, inpy
  double precision :: cval_loc, cval_max

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! NOTE: we do not use spheremp here since we want the 'max'
  ! and all 'inactive' data are just duplicates of 'active' data
  cval_loc = 0.d0
  do ie=x%nets,x%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          cval_loc = max(cval_loc, &
            abs(x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)))
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          cval_loc = max(cval_loc, &
            abs(x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)))
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  ! accumulate in a local "double precision" variable
  ! just to be safe, and then copy to cval
  call MPI_Allreduce(cval_loc, cval_max, 1, MPIreal_t, &
       MPI_MAX, par_ptr%comm, ie)
  cval = cval_max

  return
end subroutine FNVExtMaxNorm
!=======================================================================



subroutine FNVExtWrmsNorm(x_C, w_C, cval)
  !-----------------------------------------------------------------------
  ! cval = sqrt( sum_i smp(i)*[x(i)*w(i)]^2 / [4*pi*6*nlev] )
  !-----------------------------------------------------------------------
  use HommeNVector,       only: NVec_t, par_ptr
  use dimensions_mod,     only: np, nlev, nlevp
  use physical_constants, only: dd_pi
  use parallel_mod,       only: abortmp, global_shared_buf, global_shared_sum
  use global_norms_mod,   only: wrap_repro_sum
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  type(c_ptr),    intent(in)  :: w_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()
  type(NVec_t), pointer :: w => NULL()

  integer :: ie, inlev, inpx, inpy

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(w_C, w)

  ! NOTE: this use of spheremp acknowledges the fact that unknowns are
  ! duplicated, in that the sum including spheremp performs the
  ! integral over the domain.  We also use the fact that the overall
  ! spherical domain has area 4*pi
  do ie=x%nets,x%nete
    global_shared_buf(ie,1) = 0.d0
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          global_shared_buf(ie,1) = global_shared_buf(ie,1) + &
            ( &
              (x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)* &
                w%elem(ie)%state%v(inpx,inpy,1,inlev,w%tl_idx))**2 + &
              (x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)* &
                w%elem(ie)%state%v(inpx,inpy,2,inlev,w%tl_idx))**2 + &
              (x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)* &
                w%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,w%tl_idx))**2 + &
              (x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)* &
                w%elem(ie)%state%dp3d(inpx,inpy,inlev,w%tl_idx))**2 &
            ) * x%elem(ie)%spheremp(inpx,inpy)
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          global_shared_buf(ie,1) = global_shared_buf(ie,1) + &
            ( &
              (x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)* &
                w%elem(ie)%state%w_i(inpx,inpy,inlev,w%tl_idx))**2 + &
              (x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)* &
                w%elem(ie)%state%phinh_i(inpx,inpy,inlev,w%tl_idx))**2 &
            ) * x%elem(ie)%spheremp(inpx,inpy)
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  ! accumulate sum using wrap_repro_sum and then copy to cval and divide
  ! note that if all x_C entries are alpha and all w_C entries are 1/beta,
  ! then ||x_C||_wrms < 1 implies that alpha < beta, because
  ! sum_i smp(i)*[x(i)*w(i)]^2 = 4*pi*(nlev*6+2)*alpha^2/beta^2)
  call wrap_repro_sum(nvars=1, comm=par_ptr%comm)
  cval = sqrt(global_shared_sum(1)/(4.d0*dd_pi*(6.d0*nlev+2.d0)))

  return
end subroutine FNVExtWrmsNorm
!=======================================================================



subroutine FNVExtMin(x_C, cval)
  !-----------------------------------------------------------------------
  ! cval = min(|xvec|)
  !-----------------------------------------------------------------------
  use HommeNVector,     only: NVec_t, par_ptr
  use dimensions_mod,   only: np, nlev, nlevp
  use MPI,              only: MPI_MIN
  use parallel_mod,     only: MPIreal_t
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr),    intent(in)  :: x_C
  real(c_double), intent(out) :: cval

  type(NVec_t), pointer :: x => NULL()

  integer :: ie, inlev, inpx, inpy
  double precision :: cval_loc, cval_min

  !=======Internals ============

  ! dereference x_C pointer for NVec_t object
  call c_f_pointer(x_C, x)

  ! NOTE: we do not use spheremp here since we want the 'min'
  ! and all 'inactive' data are just duplicates of 'active' data
  cval_loc = abs(x%elem(x%nets)%state%v(1,1,1,1,x%tl_idx))
  do ie=x%nets,x%nete
    do inlev=1,nlev
      do inpy=1,np
        do inpx=1,np
          cval_loc = min(cval_loc, &
            abs(x%elem(ie)%state%v(inpx,inpy,1,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%v(inpx,inpy,2,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%vtheta_dp(inpx,inpy,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%dp3d(inpx,inpy,inlev,x%tl_idx)))
        end do ! inpx
      end do ! inpy
    end do ! inlev
    do inlev=1,nlevp
      do inpy=1,np
        do inpx=1,np
          cval_loc = min(cval_loc, &
            abs(x%elem(ie)%state%w_i(inpx,inpy,inlev,x%tl_idx)), &
            abs(x%elem(ie)%state%phinh_i(inpx,inpy,inlev,x%tl_idx)))
        end do ! inpx
      end do ! inpy
    end do ! inlev
  end do ! ie

  ! accumulate in a local "double precision" variable
  ! just to be safe, and then copy to cval
  call MPI_Allreduce(cval_loc, cval_min, 1, MPIreal_t, &
       MPI_MIN, par_ptr%comm, ie)
  cval = cval_min

  return
end subroutine FNVExtMin
!=======================================================================





!-----------------------------------------------------------------------
! NOTE: the remaining subroutines are not required when interfacing
! with ARKode, and are hence not currently implemented.  If this will
! utilize alternate SUNDIALS solvers in the future, some additional
! subset of these may need to be implemented as well.
!-----------------------------------------------------------------------




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
  cval = -1

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
  testval = -1

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

  !=======Internals ============

  ! dereference pointers for NVec_t objects
  call c_f_pointer(x_C, x)
  call c_f_pointer(y_C, y)

  ! perform vector operation
  cval = 0.d0

  return
end subroutine FNVExtMinQuotient
!=======================================================================
