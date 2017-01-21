#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef TRILINOS

!#define DEBUG_PRINT_ON   

module prim_preconditioner_mod
  
  use iso_c_binding, only : c_int

  use parallel_mod, only : haltmp, abortmp, iam
  use perf_mod, only     : t_startf, t_stopf

  use prim_jacobian_operator_mod, only : prim_jacobian_operator_init, &
                                         prim_jacobian_operator_finalize, &
                                         prim_jacobian_operator_update

  use prim_jacobian_dense_mod, only : prim_jacobian_dense_init, &
                                      prim_jacobian_dense_finalize, &
                                      compute_jacobian

  use prim_jacobian_sparse_mod, only : prim_jacobian_sparse_init, &
                                       prim_jacobian_sparse_finalize

  use prim_jacobian_check_mod, only : prim_jacobian_check_init, &
                                      prim_jacobian_check_finalize


  implicit none
  private
  save

  integer(c_int) :: ptype

  public :: prim_preconditioner_init, prim_preconditioner_finalize
  public :: haltrun, abortrun
 
contains

  subroutine prim_preconditioner_init(par, elem)
    ! ---------------------------------------------------------------------------
    ! initialize preconditioner
    ! ---------------------------------------------------------------------------
    use parallel_mod, only : parallel_t
    use element_mod, only  : element_t

    implicit none

    ! --------------------------------Interfaces---------------------------------
    interface
       subroutine get_ptype(pflag) bind(C,name='get_ptype')
         use, intrinsic :: iso_c_binding
         integer(c_int) :: pflag
       end subroutine get_ptype
    end interface
    ! ---------------------------------------------------------------------------

    ! --------------------------------Arguments----------------------------------
    type(parallel_t), intent(in) :: par
    type(element_t),  intent(in), target :: elem(:)
    ! ---------------------------------------------------------------------------

    if (iam == 1) write(*,*)  "HOMME: prim_preconditioner_init"

    ! get preconditioner type from trilinos
    call get_ptype(ptype)

    if (iam == 1) write(*,*)  "preconditioner type",ptype

    if (ptype == 1) then
       call prim_jacobian_operator_init(par, elem)

    else if (ptype == 2) then       
       call prim_jacobian_dense_init(par, elem)

    else if (ptype == 3) then
       call prim_jacobian_sparse_init(par, elem)

    else if (ptype == -1) then
       call prim_jacobian_dense_init(par, elem)
       call prim_jacobian_sparse_init(par, elem)
       call prim_jacobian_check_init(par, elem)

    end if

#ifdef DEBUG_PRINT_ON   
    if (iam == 1) write(*,*)  "<<< HOMME: prim_preconditioner_init"
#endif

  end subroutine prim_preconditioner_init



  subroutine prim_preconditioner_finalize()
    ! ---------------------------------------------------------------------------
    ! Free preconditioner data
    ! ---------------------------------------------------------------------------

    implicit none

    if (ptype == 1) then 
       call prim_jacobian_operator_finalize()
       
    else if (ptype == 2) then
       call prim_jacobian_dense_finalize()

    else if (ptype == 3) then 
       call prim_jacobian_sparse_finalize()

    else if (ptype == -1) then 
       call prim_jacobian_dense_finalize()
       call prim_jacobian_sparse_finalize()
       call prim_jacobian_check_finalize()

    end if
    
  end subroutine prim_preconditioner_finalize



  subroutine update_prec_state(xstate, nxstate, c_ptr_to_object) &
       bind(C,name='update_prec_state')
    ! ---------------------------------------------------------------------------
    ! update current state with most recent iteration value and compute the blocks
    ! of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelemd
    use prim_derived_type_mod ,only : derived_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    integer :: nets, nete, np1
    integer :: i, j, k, ie, lx
    ! ---------------------------------------------------------------------------

    call t_startf('Update_Prec_State')

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

    if (ptype == 1) then 
       call prim_jacobian_operator_update(xstate, nxstate, c_ptr_to_object)
    end if

    if (ptype == 2 .or. ptype == 3) then 
       
       nets = fptr%nets
       nete = fptr%nete
       np1  = fptr%tl%np1
       
       ! store current Newton iterate in np1 timelevel

       ! initialize 1d array index
       lx  = 1

       ! copy V1
       do ie=nets,nete
          do k=1,nlev
             do j=1,np
                do i=1,np
                   fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                   lx = lx+1
                end do
             end do
          end do
       end do
       
       ! copy V2
       do ie=nets,nete
          do k=1,nlev
             do j=1,np
                do i=1,np
                   fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
                   lx = lx+1
                end do
             end do
          end do
       end do
             
       ! copy Ps
       do ie=nets,nete
          do j=1,np
             do i=1,np
                fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do

       ! copy T
       do ie=nets,nete
          do k=1,nlev
             do j=1,np
                do i=1,np
                   fptr%base(ie)%state%T(i,j,k,np1) = xstate(lx)
                   lx = lx+1
                end do
             end do
          end do
       end do
    end if

    if (ptype == 2) then
       ! compute local dense Jacobian matrix at current Newton iterate 
       call compute_jacobian(fptr)
    end if

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "<<< HOMME: update_prec_state"
#endif
    
    call t_stopf('Update_Prec_State')
    
  end subroutine update_prec_state


  ! =============================================================================
  ! preconditioner testing subroutines
  ! =============================================================================


  subroutine prim_ident(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_ident')
    ! ---------------------------------------------------------------------------
    ! identity function
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding 
    
    implicit none 
    
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(in)        :: xstate(nxstate)
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    fxstate = xstate
    
  end subroutine prim_ident



  subroutine prim_dti_ident(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_dti_ident')
    ! ---------------------------------------------------------------------------
    ! scaled identity function
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding 

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use control_mod, only           : tstep_type
    
    implicit none 
    
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(in)        :: xstate(nxstate)
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind) :: dti ! 1/dt

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    dti = 1.0d0/fptr%dt 
        
    if (tstep_type == 12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    end if

    fxstate = dti * xstate
    
  end subroutine prim_dti_ident



  subroutine prim_dt_ident(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_dt_ident')
    ! ---------------------------------------------------------------------------
    ! scaled identity function
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding 

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use control_mod, only           : tstep_type
    
    implicit none 
    
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(in)        :: xstate(nxstate)
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind) :: dt

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    dt = fptr%dt 
        
    if (tstep_type == 12) then
       dt = (2.0d0/3.0d0)*dt ! BDF2 has a factor of 3/2
    end if

    fxstate = dt * xstate
    
  end subroutine prim_dt_ident


  ! =============================================================================
  ! subroutines to stop run
  ! =============================================================================


  subroutine haltrun() bind(C,name='haltrun')
    ! ---------------------------------------------------------------------------
    ! wrapper to MPI finalize
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    call haltmp("halt called from c++")

  end subroutine haltrun


  subroutine abortrun() bind(C,name='abortrun')
    ! ---------------------------------------------------------------------------
    ! wrapper to MPI about
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    call abortmp("abort called from c++")

  end subroutine abortrun
   
end module prim_preconditioner_mod

#endif
! TRILINOS
