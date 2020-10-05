!-----------------------------------------------------------------------------
! Daniel R. Reynolds, SMU Mathematics
! Christopher J. Vogl, LLNL
! David J. Gardner, LLNL
!
! Copyright 2017; all rights reserved
!-----------------------------------------------------------------------------
! NOTE: Only the vector operations required by ARKode have been
! implemented. If this module is utilized with other SUNDIALS solvers
! in the future, some additional vector operations may need to be
! implemented as well.
!-----------------------------------------------------------------------------

module HommeNVector
  !-----------------------------------------------------------------------------
  ! Fortran user-defined NVector interface
  !-----------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod

  use element_mod,   only: element_t
  use element_state, only: timelevels
  use parallel_mod,  only: parallel_t

  implicit none
  public

  ! Fortran object stored as N_Vector content
  type :: NVec_t
     type(element_t), pointer :: elem(:)  ! element data
     integer                  :: nets     ! starting index
     integer                  :: nete     ! ending index
     integer                  :: tl_idx   ! time level index
     logical                  :: own_data ! vec owns data
  end type NVec_t

  save

  ! number of vectors available in the registry (time levels)
  integer, parameter :: RegistryLength=timelevels

  ! indicates if a slot in the registry is in use
  logical :: HommeNVectorRegistry(RegistryLength)=.false.

  ! parallelization data
  type(parallel_t), pointer :: par_ptr

contains

  !=============================================================================

  subroutine SetHommeNVectorPar(par)
    !---------------------------------------------------------------------------
    ! sets current par object pointer for MPI_Allreduce and wrap_repro_sum calls
    !---------------------------------------------------------------------------
    implicit none
    type(parallel_t), target, intent(in) :: par
    par_ptr => par
  end subroutine SetHommeNVectorPar

  !=============================================================================

  integer function ReserveHommeNVectorRegistryIdx()
    !---------------------------------------------------------------------------
    ! Marks first available vector from registry as used (if any are available)
    !---------------------------------------------------------------------------
    implicit none
    integer :: i

    !======= Internals ============

    ! search for first free slot in the registry
    do i=1,RegistryLength
       if (.not. HommeNVectorRegistry(i)) then
          HommeNVectorRegistry(i) = .true.
          ReserveHommeNVectorRegistryIdx = i
          return
       end if
    end do

    ! no open slot found
    ReserveHommeNVectorRegistryIdx = 0

    return

  end function ReserveHommeNVectorRegistryIdx

  !=============================================================================

  function MakeHommeNVector(elem, nets, nete, tl_idx, ierr) result(sunvec_y)
    !---------------------------------------------------------------------------
    ! Function to create an N_Vector wrapper (sunvec_y) to the HOMME solution
    !   structure, and reserve its index in the vector registry.
    ! The value of ierr is 0 if succssful, 1 if failure (if the registry
    !   index is already taken)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use element_mod, only: element_t

    implicit none
    type (element_t), target    :: elem(:)
    integer, intent(in)         :: nets
    integer, intent(in)         :: nete
    integer, intent(in)         :: tl_idx
    integer(c_int), intent(out) :: ierr
    type(N_Vector), pointer     :: sunvec_y

    type(N_Vector_Ops), pointer :: ops
    type(NVec_t), pointer       :: content

    !======= Internals ============

    ! initialize error flag
    ierr = 0

    ! if registry index is already taken, return failure
    if (HommeNVectorRegistry(tl_idx)) then
       ierr = 1
       sunvec_y => null()
       return
    end if

    ! allocate output N_Vector structure
    sunvec_y => FN_VNewEmpty()

    ! allocate and fill content structure
    allocate(content)
    content%elem     => elem
    content%nets     =  nets
    content%nete     =  nete
    content%tl_idx   =  tl_idx
    content%own_data =  .false.

    ! mark the index in the registry as used
    HommeNVectorRegistry(tl_idx) = .true.

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(content)

    ! access the N_Vector ops structure
    call c_f_pointer(sunvec_y%ops, ops)

    ! set function pointers to vector operations
    ops%nvgetvectorid = c_funloc(FN_VGetVectorID_HOMME)
    ops%nvclone       = c_funloc(FN_VClone_HOMME)
    ops%nvdestroy     = c_funloc(FN_VDestroy_HOMME)
    ops%nvconst       = c_funloc(FN_VConst_HOMME)
    ops%nvdotprod     = c_funloc(FN_VDotProd_HOMME)
    ops%nvlinearsum   = c_funloc(FN_VLinearSum_HOMME)
    ops%nvprod        = c_funloc(FN_VProd_HOMME)
    ops%nvdiv         = c_funloc(FN_VDiv_HOMME)
    ops%nvscale       = c_funloc(FN_VScale_HOMME)
    ops%nvabs         = c_funloc(FN_VAbs_HOMME)
    ops%nvinv         = c_funloc(FN_VInv_HOMME)
    ops%nvaddconst    = c_funloc(FN_VAddConst_HOMME)
    ops%nvmaxnorm     = c_funloc(FN_VMaxNorm_HOMME)
    ops%nvwrmsnorm    = c_funloc(FN_VWRMSNorm_HOMME)
    ops%nvmin         = c_funloc(FN_VMin_HOMME)
    ops%nvprint       = c_funloc(FN_VPrint_HOMME)

    return

  end function MakeHommeNVector

  !=============================================================================

  function FN_VGetContent(sunvec_y) result(y)
    !---------------------------------------------------------------------------
    ! Get the Fortran object from the N_Vector content
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none
    type(N_Vector)          :: sunvec_y
    type(NVec_t), pointer   :: y

    !======= Internals ============

    ! extract Fortran data from N_Vector
    call c_f_pointer(sunvec_y%content, y)

    return

  end function FN_VGetContent

  !=============================================================================

  subroutine FN_VPrint_HOMME(sunvec_x)
    !---------------------------------------------------------------------------
    ! Print routine
    !
    ! Note: this function is not required by ARKode (or any of SUNDIALS)
    !       it is merely here for convenience when debugging
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp
    use MPI,            only: MPI_comm_rank

    implicit none
    type(N_Vector)        :: sunvec_x
    type(NVec_t), pointer :: x => NULL()

    integer :: ie, inpx, inpy, inlev, rank, ierr

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)

    ! get rank
    call MPI_comm_rank(par_ptr%comm, rank, ierr)

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

  end subroutine FN_VPrint_HOMME

  !=============================================================================

  integer(N_Vector_ID) function FN_VGetVectorID_HOMME(sunvec_y) &
       result(id) bind(C)
    !---------------------------------------------------------------------------
    ! Return the vector ID
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none
    type(N_Vector) :: sunvec_y

    !======= Internals ============

    id = SUNDIALS_NVEC_CUSTOM

    return

  end function FN_VGetVectorID_HOMME

  !=============================================================================

  function FN_VClone_HOMME(sunvec_x) result(y_ptr) bind(C)
    !---------------------------------------------------------------------------
    ! Clone routine for EXT vector -- allocates memory for a new NVec_t y
    ! that that matches structure for x.
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    use parallel_mod, only: abortmp

    implicit none
    type(N_Vector) :: sunvec_x
    type(c_ptr)    :: y_ptr

    type(N_Vector), pointer :: sunvec_y
    integer(c_int)          :: retval
    type(NVec_t), pointer   :: x => NULL()
    type(NVec_t), pointer   :: y => NULL()
    integer                 :: tl_idx

    !======= Internals ============

    ! grab next available registry index
    tl_idx = ReserveHommeNVectorRegistryIdx()
    if (tl_idx == 0) then
      call abortmp('FN_VClone_HOMME failed')
    end if

    ! allocate output N_Vector structure
    sunvec_y => FN_VNewEmpty()
    if (.not.associated(sunvec_y)) then
      call abortmp('FN_VClone_HOMME failed (unassociated vector)')
    end if

    ! copy operations from x into y
    retval = FN_VCopyOps(sunvec_x, sunvec_y)
    if (retval /= 0) then
      call abortmp('FN_VCopyOps failed')
    end if

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)

    ! allocate and fill content structure
    allocate(y)
    y%elem     => x%elem
    y%nets     =  x%nets
    y%nete     =  x%nete
    y%tl_idx   =  tl_idx
    y%own_data = .true.

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(y)

    ! set the c_ptr output
    y_ptr = c_loc(sunvec_y)

    return

  end function FN_VClone_HOMME

  !=============================================================================

  subroutine FN_VDestroy_HOMME(sunvec_x) bind(C)
    !---------------------------------------------------------------------------
    ! Destroy routine for EXT vector -- nullifies data for a given NVec_t
    ! and returns its vector index to the registry
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none
    type(N_Vector), target :: sunvec_x
    type(NVec_t), pointer  :: x => NULL()

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)

    ! mark the index in the registry as unused
    if (x%own_data) HommeNVectorRegistry(x%tl_idx) = .false.

    ! detach pointer
    if (associated(x%elem)) nullify(x%elem)

    ! deallocate the underlying Fortran object (the content)
    deallocate(x)

    ! set N_Vector structure members to NULL
    sunvec_x%content = c_null_ptr

    ! deallocate overall N_Vector structure
    call FN_VFreeEmpty(sunvec_x)

    return

  end subroutine FN_VDestroy_HOMME

  !=============================================================================

  subroutine FN_VConst_HOMME(cval, sunvec_z) bind(C)
    !---------------------------------------------------------------------------
    ! z = c
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none
    type(N_Vector)             :: sunvec_z
    real(c_double), value      :: cval

    type(NVec_t), pointer :: z => NULL()

    integer :: ie

    !======= Internals ============

    ! access NVec_t object
    z => FN_VGetContent(sunvec_z)

    ! perform vector operation
    do ie=z%nets,z%nete
       z%elem(ie)%state%v(:,:,:,:,z%tl_idx) = cval
       z%elem(ie)%state%w_i(:,:,:,z%tl_idx) = cval
       z%elem(ie)%state%phinh_i(:,:,:,z%tl_idx) = cval
       z%elem(ie)%state%vtheta_dp(:,:,:,z%tl_idx) = cval
       z%elem(ie)%state%dp3d(:,:,:,z%tl_idx) = cval
    end do

    return

  end subroutine FN_VConst_HOMME

  !=============================================================================

  real(c_double) function FN_VDotProd_HOMME(sunvec_x, sunvec_y) &
       result(cval) bind(C)
    !---------------------------------------------------------------------------
    ! c = <x,y>  (only include 'active' data; no ghost cells, etc.)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod,   only: np, nlev, nlevp
    use parallel_mod,     only: global_shared_buf, global_shared_sum
    use global_norms_mod, only: wrap_repro_sum

    implicit none
    type(N_Vector)             :: sunvec_x
    type(N_Vector)             :: sunvec_y

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: y => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    y => FN_VGetContent(sunvec_y)

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

  end function FN_VDotProd_HOMME

  !=============================================================================

  subroutine FN_VLinearSum_HOMME(aval, sunvec_x, bval, sunvec_y, sunvec_z) &
       bind(C)
    !---------------------------------------------------------------------------
    ! z = a*x + b*y
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp

    implicit none
    real(c_double), value :: aval
    type(N_Vector)        :: sunvec_x
    real(c_double), value :: bval
    type(N_Vector)        :: sunvec_y
    type(N_Vector)        :: sunvec_z

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: y => NULL()
    type(NVec_t), pointer :: z => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    y => FN_VGetContent(sunvec_y)
    z => FN_VGetContent(sunvec_z)

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

  end subroutine FN_VLinearSum_HOMME

  !=============================================================================

  subroutine FN_VProd_HOMME(sunvec_x, sunvec_y, sunvec_z) bind(C)
    !---------------------------------------------------------------------------
    ! z = x.*y (Matlab notation)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp

    implicit none
    type(N_Vector)             :: sunvec_x
    type(N_Vector)             :: sunvec_y
    type(N_Vector)             :: sunvec_z

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: y => NULL()
    type(NVec_t), pointer :: z => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    y => FN_VGetContent(sunvec_y)
    z => FN_VGetContent(sunvec_z)

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

  end subroutine FN_VProd_HOMME

  !=============================================================================

  subroutine FN_VDiv_HOMME(sunvec_x, sunvec_y, sunvec_z) bind(C)
    !---------------------------------------------------------------------------
    ! z = x./y (Matlab notation; doesn't check for legal denominator)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp

    implicit none
    type(N_Vector)             :: sunvec_x
    type(N_Vector)             :: sunvec_y
    type(N_Vector)             :: sunvec_z

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: y => NULL()
    type(NVec_t), pointer :: z => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    y => FN_VGetContent(sunvec_y)
    z => FN_VGetContent(sunvec_z)

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

  end subroutine FN_VDiv_HOMME

  !=============================================================================

  subroutine FN_VScale_HOMME(cval, sunvec_x, sunvec_z) bind(C)
    !---------------------------------------------------------------------------
    ! z = c*x
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp

    implicit none
    real(c_double), value      :: cval
    type(N_Vector)             :: sunvec_x
    type(N_Vector)             :: sunvec_z

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: z => NULL()

    integer :: ie, inlev, inpy, inpx

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    z => FN_VGetContent(sunvec_z)

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

  end subroutine FN_VScale_HOMME

  !=============================================================================

  subroutine FN_VAbs_HOMME(sunvec_x, sunvec_z) bind(C)
    !---------------------------------------------------------------------------
    ! z = |x|  (componentwise)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp

    implicit none
    type(N_Vector)             :: sunvec_x
    type(N_Vector)             :: sunvec_z

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: z => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    z => FN_VGetContent(sunvec_z)

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

  end subroutine FN_VAbs_HOMME

  !=============================================================================

  subroutine FN_VInv_HOMME(sunvec_x, sunvec_z) bind(C)
    !---------------------------------------------------------------------------
    ! z = 1./x (Matlab notation: doesn't check for division by zero)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp

    implicit none
    type(N_Vector)             :: sunvec_x
    type(N_Vector)             :: sunvec_z

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: z => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    z => FN_VGetContent(sunvec_z)

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

  end subroutine FN_VInv_HOMME

  !=============================================================================

  subroutine FN_VAddConst_HOMME(cval, sunvec_x, sunvec_z) bind(C)
    !---------------------------------------------------------------------------
    ! z = x+c (add c to each component of the vector x)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod, only: np, nlev, nlevp

    implicit none
    real(c_double), intent(in) :: cval
    type(N_Vector)             :: sunvec_x
    type(N_Vector)             :: sunvec_z

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: z => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    z => FN_VGetContent(sunvec_z)

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

  end subroutine FN_VAddConst_HOMME

  !=============================================================================

  real(C_double) function FN_VMaxNorm_HOMME(sunvec_x) result(cval) bind(C)
    !---------------------------------------------------------------------------
    ! c = max(|x|)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod,   only: np, nlev, nlevp
    use MPI,              only: MPI_MAX
    use parallel_mod,     only: MPIreal_t

    implicit none
    type(N_Vector)             :: sunvec_x

    type(NVec_t), pointer :: x => NULL()

    integer :: ie, inlev, inpx, inpy
    double precision :: cval_loc, cval_max

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)

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

  end function FN_VMaxNorm_HOMME

  !=============================================================================

  real(c_double) function FN_VWrmsNorm_HOMME(sunvec_x, sunvec_w) &
       result(cval) bind(C)
    !---------------------------------------------------------------------------
    ! cval = sqrt( sum_i smp(i)*[x(i)*w(i)]^2 / [4*pi*6*nlev] )
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod,     only: np, nlev, nlevp
    use physical_constants, only: dd_pi
    use parallel_mod,       only: abortmp, global_shared_buf, global_shared_sum
    use global_norms_mod,   only: wrap_repro_sum

    implicit none
    type(N_Vector), intent(in)  :: sunvec_x
    type(N_Vector), intent(in)  :: sunvec_w

    type(NVec_t), pointer :: x => NULL()
    type(NVec_t), pointer :: w => NULL()

    integer :: ie, inlev, inpx, inpy

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)
    w => FN_VGetContent(sunvec_w)

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

  end function FN_VWrmsNorm_HOMME

  !=============================================================================

  real(c_double) function FN_VMin_HOMME(sunvec_x) result(cval) bind(C)
    !---------------------------------------------------------------------------
    ! cval = min(|xvec|)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use dimensions_mod,   only: np, nlev, nlevp
    use MPI,              only: MPI_MIN
    use parallel_mod,     only: MPIreal_t

    implicit none
    type(N_Vector)             :: sunvec_x

    type(NVec_t), pointer :: x => NULL()

    integer :: ie, inlev, inpx, inpy
    double precision :: cval_loc, cval_min

    !======= Internals ============

    ! access NVec_t object
    x => FN_VGetContent(sunvec_x)

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

  end function FN_VMin_HOMME

  !=============================================================================

end module HommeNVector
