#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef TRILINOS

! ==============================================================================
! This module contains subroutines for comparing the analytic Jacobian
! implemented in prim_jacobian_dense_mod.f90 against a second order centered 
! finite difference approximation. 
! ==============================================================================

! output using local index (otherwise uses global index)
!#define LOCAL_IDX

! perform DSS exchange when computing finite difference
! #define DSS_ON

! print analytic Jacobian
#define PRINT_AJ

! print finite difference Jacobian
#define PRINT_FDJ

module prim_jacobian_check_mod

  use kinds, only        : real_kind
  use edgetype_mod, only : EdgeBuffer_t
  use edge_mod, only     : edgevpack, edgevunpack
  use bndry_mod, only    : bndry_exchangev
  use parallel_mod, only : abortmp, haltmp, iam
  use perf_mod, only     : t_startf, t_stopf

  use prim_jacobian_dense_mod

  use prim_jacobian_sparse_mod, only : Uoffset, Voffset, Toffset, Psoffset, &
       print_global_IDs
  
  implicit none 
  private

  ! exchange buffer
  type(EdgeBuffer_t) :: edge3p1

  ! power on epsilon 10^eps_pow_start -> 10^eps_pow_stop
  integer, parameter :: eps_pow_start = -1
  integer, parameter :: eps_pow_stop  = -1
    
  public :: check_jacobian
  public :: prim_jacobian_check_init, prim_jacobian_check_finalize

  ! define interfaces for pointers to residual functions
  abstract interface

     ! residual components on model levels
     subroutine resid(xstate, nxstate, fptr, utens, vtens, ttens, pstens)
       use, intrinsic :: iso_c_binding

       use kinds, only                 : real_kind
       use dimensions_mod, only        : np, nlev, nvar, nelemd
       use prim_derived_type_mod, only : derived_type

       implicit none
       
       real(c_double), intent(in)        :: xstate(nxstate)
       integer(c_int), intent(in), value :: nxstate
       type(derived_type), pointer       :: fptr
       
       real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: utens, vtens, ttens
       real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: pstens
     end subroutine resid

     ! auxiliary residual components on model levels
     subroutine resid_aux(xstate, nxstate, fptr, ftens)
       use, intrinsic :: iso_c_binding

       use kinds, only                 : real_kind
       use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
       use prim_derived_type_mod, only : derived_type

       implicit none
       
       real(c_double), intent(in)        :: xstate(nxstate)
       integer(c_int), intent(in), value :: nxstate
       type(derived_type), pointer       :: fptr
       
       real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: ftens
     end subroutine resid_aux

     ! auxiliary residual components on model interfaces
     subroutine resid_auxi(xstate, nxstate, fptr, ftens)
       use, intrinsic :: iso_c_binding

       use kinds, only                 : real_kind
       use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
       use prim_derived_type_mod, only : derived_type

       implicit none
       
       real(c_double), intent(in)        :: xstate(nxstate)
       integer(c_int), intent(in), value :: nxstate
       type(derived_type), pointer       :: fptr
       
       real(kind=real_kind), dimension(np,np,nlev+1,nelemd), intent(out) :: ftens
     end subroutine resid_auxi

     ! Jacobian blocks
     subroutine jac(fptr)
       use prim_derived_type_mod, only : derived_type

       implicit none

       type(derived_type), pointer, intent(in) :: fptr
     end subroutine jac

     ! auxiliary Jacobian calculations at model levels
     subroutine jac_aux(du_vals, dv_vals, dps_vals, fptr, ie)
       use kinds, only                 : real_kind
       use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
       use prim_derived_type_mod, only : derived_type

       implicit none

       real(kind=real_kind), intent(out) :: du_vals(np,np,nlev,np,np,nlev)
       real(kind=real_kind), intent(out) :: dv_vals(np,np,nlev,np,np,nlev)
       real(kind=real_kind), intent(out) :: dps_vals(np,np,np,np,nlev)    

       type(derived_type), pointer, intent(in) :: fptr
       integer, intent(in) :: ie 
     end subroutine jac_aux

     ! auxiliary Jacobian calculations at model interfaces
     subroutine jac_auxi(du_vals, dv_vals, dps_vals, fptr, ie)
       use kinds, only                 : real_kind
       use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
       use prim_derived_type_mod, only : derived_type

       implicit none

       real(kind=real_kind), intent(out) :: du_vals(np,np,nlev,np,np,nlev+1)
       real(kind=real_kind), intent(out) :: dv_vals(np,np,nlev,np,np,nlev+1)
       real(kind=real_kind), intent(out) :: dps_vals(np,np,np,np,nlev+1)    

       type(derived_type), pointer, intent(in) :: fptr
       integer, intent(in) :: ie 
     end subroutine jac_auxi

  end interface
 
contains

  subroutine prim_jacobian_check_init(par, elem)

    use parallel_mod, only   : parallel_t
    use element_mod, only    : element_t
    use edge_mod, only       : initEdgeBuffer
    use dimensions_mod, only : np, nlev, nelemd
    use control_mod, only    : rsplit

    implicit none

    type(parallel_t), intent(in) :: par
    type(element_t),  intent(in), target :: elem(:)

    ! allocate exchange buffe
    if (rsplit==0) then
       call initEdgeBuffer(par, edge3p1, elem, 3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par, edge3p1, elem, 4*nlev+1)
    endif

  end subroutine prim_jacobian_check_init



  subroutine prim_jacobian_check_finalize

    use edge_mod, only : FreeEdgeBuffer

    implicit none

    call FreeEdgeBuffer(edge3p1)

  end subroutine prim_jacobian_check_finalize


  
  subroutine set_state_vector(xstate, nxstate, c_ptr_to_object) &
              bind(C,name='set_state_vector')
    ! ---------------------------------------------------------------------------
    ! Overwrite state vector with random values for testing Jacobian assembly
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
   
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelem, nelemd
    use prim_derived_type_mod, only : derived_type

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(inout)     :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()
    
    real(kind=real_kind) :: rval(nxstate) ! array of random numbers
   
    integer :: nets, nete, np1, n0, nm1
    integer :: i, j, k, ie, lx
    integer :: a, b, c, h, l, ierr
    ! ---------------------------------------------------------------------------

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

    if (fptr%hybrid%masterthread) &
         write(*,*) ">>> Setting Random State Vector <<<"
    
    ! element indices and time level
    nets = fptr%nets
    nete = fptr%nete
    np1  = fptr%tl%np1   
    n0   = fptr%tl%n0
    nm1  = fptr%tl%nm1
      
    ! array of uniformly distributed random values in [-1,1)
    call random_number(rval)
    rval = 2.0d0 * rval - 1.0d0

    ! u velocity values
    lx = 1
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                ! xstate(lx) = 50.0d0 + 5.0d0 * rval(lx)
                xstate(lx) = 20.0d0 + 10.0d0 * rval(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! v velocity values
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                ! xstate(lx) = 50.0d0 + 5.0d0 * rval(lx)
                xstate(lx) = 20.0d0 + 10.0d0 * rval(lx)
                lx = lx+1
             end do
          end do
       end do
    end do

    ! surface pressure values
    do ie=nets,nete
       do j=1,np
          do i=1,np
             ! xstate(lx) = 50.0d0 + 5.0d0 * rval(lx)
             xstate(lx) = 1.0d5 + 500.0d0 * rval(lx)
             lx = lx+1
          end do
       end do
    end do

    ! temperature values
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                ! xstate(lx) = 50.0d0 + 5.0d0 * rval(lx)
                xstate(lx) = 260.0d0 + 40.0d0 * rval(lx)
                lx = lx+1
             end do
          end do
       end do
    end do

    ! set Ps0 value
    ! fptr%hvcoord%ps0 = 100.0d0

    ! set Phi_s value
    ! fptr%base(ie)%state%phis

    ! copy state vector into pointer for HOMME derived type
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)
    call copy_state_to_HOMME(xstate, nxstate, fptr, n0)
    call copy_state_to_HOMME(xstate, nxstate, fptr, nm1)

  end subroutine set_state_vector


  
  subroutine check_jacobian(xstate, nxstate, c_ptr_to_object) &
       bind(C,name='check_jacobian')
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix against a finite 
    ! difference approximation. The default is to use a second order centered 
    ! difference, uncomment '#define FD1' at the top of this file to use a first
    ! order forward difference approximation
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
   
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelem, nelemd
    use prim_derived_type_mod, only : derived_type

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    ! jacobian and residual subroutines   
    procedure(jac),   pointer :: jacobian_fn 
    procedure(resid), pointer :: residual_fn 

    ! auxiliary jacobian and residual subroutines
    procedure(jac_aux),   pointer :: jacobian_fn_aux 
    procedure(resid_aux), pointer :: residual_fn_aux 

    ! auxiliary interface jacobian and residual subroutines
    procedure(jac_auxi),   pointer :: jacobian_fn_auxi 
    procedure(resid_auxi), pointer :: residual_fn_auxi
    
    integer :: nets, nete, np1, n0, nm1
    integer :: i, j, k, ie, lx
    integer :: a, b, c, h, l, ierr

    character(len=1024) :: Jname ! name of Jacobian subroutine
    ! ---------------------------------------------------------------------------

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

    if (fptr%hybrid%masterthread) &
         write(*,*) ">>> Outputting Jacobian to File <<<"

    ! print element data
    call print_elem_data(fptr)

    ! ------------------------------------------------------------------------
    ! check aux terms for div_dpdn_vel
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'div_dpdn_vel'
    write(*,'(a)') 'div_dpdn_vel'
    
    jacobian_fn_aux => Jac_div_dpdn_vel
    residual_fn_aux => residual_aux_div_dpdn_vel
    call print_jacobian_aux(xstate, nxstate, fptr, residual_fn_aux, jacobian_fn_aux, Jname)

    ! ------------------------------------------------------------------------
    ! check aux interface terms for etadot_dpdn
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'etadot_dpdn'
    write(*,'(a)') 'etadot_dpdn'

    jacobian_fn_auxi => Jac_etadot_dpdn
    residual_fn_auxi => residual_aux_etadot_dpdn
    call print_jacobian_auxi(xstate, nxstate, fptr, residual_fn_auxi, jacobian_fn_auxi, Jname)

    ! ------------------------------------------------------------------------
    ! check temporal terms
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'time'
    write(*,'(a)') 'time'

    residual_fn => residual_time
    jacobian_fn => Jac_TimeDeriv
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check vorticity terms
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'vort'
    write(*,'(a)') 'vort'

    residual_fn => residual_vorticity
    jacobian_fn => Jac_Vorticity
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check coriolis terms
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'cor'
    write(*,'(a)') 'cor'

    residual_fn => residual_coriolis
    jacobian_fn => Jac_Coriolis
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check kinetic energy
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'KE'
    write(*,'(a)') 'KE'

    residual_fn => residual_kinetic_energy
    jacobian_fn => Jac_Kinetic_Energy
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check geopotential
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'GP'
    write(*,'(a)') 'GP'

    residual_fn => residual_geopotential
    jacobian_fn => Jac_Geopotential
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check pressure gradient
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'gradp'
    write(*,'(a)') 'gradp'

    residual_fn => residual_pressure_grad
    jacobian_fn => Jac_press_grad
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check vertical advection
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'vertadv'
    write(*,'(a)') 'vertadv'

    residual_fn => residual_vertadv
    jacobian_fn => Jac_VertAdv
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check horizontal temperature advection
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'horizadv'
    write(*,'(a)') 'horizadv'

    residual_fn => residual_HorizAdv_temp
    jacobian_fn => Jac_HorizAdv_Temp
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check pressure vertical velocity
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'pvertvel'
    write(*,'(a)') 'pvertvel'
    
    residual_fn => residual_pressure_vert_vel
    jacobian_fn => Jac_pressure_vert_vel
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check surface pressure equation
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'intps'
    write(*,'(a)') 'intps'
    
    residual_fn => residual_int_div_dpdn_vel
    jacobian_fn => Jac_int_div_dpdn_vel
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)

    ! ------------------------------------------------------------------------
    ! check full Jacobian matrix
    ! ------------------------------------------------------------------------
    write(Jname,'(a)') 'Jac'
    write(*,'(a)') 'Jac'

    residual_fn => residual_full
    jacobian_fn => compute_jacobian
    call print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)
        
  end subroutine check_jacobian



  subroutine print_jacobian(xstate, nxstate, fptr, residual_fn, jacobian_fn, Jname)
    ! ---------------------------------------------------------------------------
    ! Outputs finite difference Jacobian to file
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    procedure(resid), pointer :: residual_fn ! pointer to Residual subroutine
    procedure(jac),   pointer :: jacobian_fn ! pointer to Jacobian subroutine

    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(c_double) :: ystate(nxstate)

    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: ttens1, ttens2
    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: diff_ttens

    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: utens1,utens2
    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: diff_utens

    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: vtens1, vtens2
    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: diff_vtens

    real(kind=real_kind), dimension(np,np,nelemd) :: pstens1, pstens2
    real(kind=real_kind), dimension(np,np,nelemd) :: diff_pstens

    real(kind=real_kind) :: eps   ! scaling for finite difference perturbation
    real(kind=real_kind) :: sigma ! finite difference perturbation

    integer :: nets, nete, np1    ! element indices, time level
    integer :: i, a, b, c, lx, ie ! loop indices
    ! ---------------------------------------------------------------------------
   
#ifdef PRINT_AJ
    ! compute and output analytic Jacobian 
    np1 = fptr%tl%np1

    call zero_jacobian()
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)
    call jacobian_fn(fptr)
    call print_analytic_jacobian(fptr, Jname)
#endif

#ifdef PRINT_FDJ 
    ! prints 1 column at a time

    ! starting and ending element index
    nets = fptr%nets
    nete = fptr%nete

    ! make copy the original unpreturbed state vector
    ystate = xstate

    do i = eps_pow_start, eps_pow_stop
       
       ! perturbation size for finite difference 
       eps = 10.0d0**i

       ! initialize index in 1D state vector
       lx = 0

       ! Derivatives with respect to U
       do ie = nets,nete
          do c = 1,nlev
             do b = 1,np
                do a = 1,np

                   ! update position in 1D state vector
                   lx = lx+1

                   ! compute perturbed residual values
                   sigma = abs(xstate(lx)) * eps

                   ystate(lx) = xstate(lx) + sigma
                   call residual_fn(ystate, nxstate, fptr, &
                        utens1, vtens1, ttens1, pstens1)

                   ystate(lx) = xstate(lx) - sigma
                   call residual_fn(ystate, nxstate, fptr, &
                        utens2, vtens2, ttens2, pstens2)

                   sigma = 2.0d0 * sigma

                   ! finite difference Jacobian
                   diff_utens  = (utens1 - utens2)/sigma
                   diff_vtens  = (vtens1 - vtens2)/sigma
                   diff_ttens  = (ttens1 - ttens2)/sigma
                   diff_pstens = (pstens1 - pstens2)/sigma

                   ! print finite difference Jacobian
                   call print_du(Jname, fptr, a, b, c, ie, &
                        diff_utens, diff_vtens, diff_ttens, diff_pstens)

                   ! restore perturbed value
                   ystate(lx) = xstate(lx)
                
                end do
             end do
          end do
       end do
  
       ! ------------------------------------------------------------------------
       ! Derivatives with respect to V
       ! ------------------------------------------------------------------------
       do ie = nets,nete
          do c = 1,nlev
             do b = 1,np
                do a = 1,np

                   ! update position in 1D state vector
                   lx = lx+1

                   ! compute perturbed residual values
                   sigma = abs(xstate(lx)) * eps

                   ystate(lx) = xstate(lx) + sigma
                   call residual_fn(ystate, nxstate, fptr, &
                        utens1, vtens1, ttens1, pstens1)

                   ystate(lx) = xstate(lx) - sigma
                   call residual_fn(ystate, nxstate, fptr, &
                        utens2, vtens2, ttens2, pstens2)

                   sigma = 2.0d0 * sigma

                   ! finite difference Jacobian
                   diff_utens  = (utens1 - utens2)/sigma
                   diff_vtens  = (vtens1 - vtens2)/sigma
                   diff_ttens  = (ttens1 - ttens2)/sigma
                   diff_pstens = (pstens1 - pstens2)/sigma

                   ! print finite difference Jacobian
                   call print_dv(Jname, fptr, a, b, c, ie, & 
                        diff_utens, diff_vtens, diff_ttens, diff_pstens)

                   ! restore perturbed value
                   ystate(lx) = xstate(lx)

                end do
             end do
          end do
       end do

       ! ------------------------------------------------------------------------
       ! Derivatives with respect to Ps
       ! ------------------------------------------------------------------------
       do ie = nets,nete
          do b = 1,np
             do a = 1,np

                ! update position in 1D state vector
                lx = lx+1

                ! compute perturbed residual values
                sigma = abs(xstate(lx)) * eps

                ystate(lx) = xstate(lx) + sigma
                call residual_fn(ystate, nxstate, fptr, &
                     utens1, vtens1, ttens1, pstens1)

                ystate(lx) = xstate(lx) - sigma
                call residual_fn(ystate, nxstate, fptr, &
                     utens2, vtens2, ttens2, pstens2)

                sigma = 2.0d0 * sigma

                ! finite difference Jacobian
                diff_utens  = (utens1 - utens2)/sigma
                diff_vtens  = (vtens1 - vtens2)/sigma
                diff_ttens  = (ttens1 - ttens2)/sigma
                diff_pstens = (pstens1 - pstens2)/sigma

                ! print finite difference Jacobian
                call print_dPs(Jname, fptr, a, b, ie, &
                     diff_utens, diff_vtens, diff_ttens, diff_pstens)

                ! restore perturbed value
                ystate(lx) = xstate(lx)

             end do
          end do
       end do

       ! ------------------------------------------------------------------------
       ! Derivatives with respect to T
       ! ------------------------------------------------------------------------
       do ie = nets,nete
          do c = 1,nlev
             do b = 1,np
                do a = 1,np

                   ! update position in 1D state vector
                   lx = lx+1

                   ! compute perturbed residual values
                   sigma = abs(xstate(lx)) * eps

                   ystate(lx) = xstate(lx) + sigma
                   call residual_fn(ystate, nxstate, fptr, &
                        utens1, vtens1, ttens1, pstens1)

                   ystate(lx) = xstate(lx) - sigma
                   call residual_fn(ystate, nxstate, fptr, &
                        utens2, vtens2, ttens2, pstens2)

                   sigma = 2.0d0 * sigma

                   ! finite difference Jacobian
                   diff_utens  = (utens1 - utens2)/sigma
                   diff_vtens  = (vtens1 - vtens2)/sigma
                   diff_ttens  = (ttens1 - ttens2)/sigma
                   diff_pstens = (pstens1 - pstens2)/sigma

                   ! print finite difference Jacobian
                   call print_dT(Jname, fptr, a, b, c, ie, &
                        diff_utens, diff_vtens, diff_ttens, diff_pstens)

                   ! restore perturbed value
                   ystate(lx) = xstate(lx)

                end do
             end do
          end do
       end do

    end do ! loop over epsilon

#endif
    
  end subroutine print_jacobian



  subroutine print_jacobian_aux(xstate, nxstate, fptr, residual_fn_aux, &
       jacobian_fn_aux, Jname)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    procedure(resid_aux), pointer :: residual_fn_aux ! Residual subroutine
    procedure(jac_aux),   pointer :: jacobian_fn_aux ! Jacobian subroutine

    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    character(len=1024) :: fname ! name of output file

    logical :: fexist      ! if output file exists
    integer :: fid         ! file ID for output

    real(c_double) :: ystate(nxstate)

    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: ftens1, ftens2
    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: diff_ftens

    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev) :: du_vals, dv_vals
    real(kind=real_kind), dimension(np,np,np,np,nlev)      :: dps_vals

    real(kind=real_kind) :: eps   ! scaling for finite difference perturbation
    real(kind=real_kind) :: sigma ! finite difference perturbation

    integer :: ie, nets, nete, np1    ! element indices, time level
    integer :: i, j, k, ridx
    integer :: a, b, c, cidx
    integer :: ii, lx, e
    ! ---------------------------------------------------------------------------

#ifdef PRINT_AJ
    ! compute and output analytic Jacobian 
    np1 = fptr%tl%np1

    do ie = 1,nelemd
       
       call copy_state_to_HOMME(xstate, nxstate, fptr, np1)
       call jacobian_fn_aux(du_vals, dv_vals, dps_vals, fptr, ie)

       ! ------------------------------------------------------------------------
       ! print du to file 
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'du_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np
                         
                         ! local column index
                         cidx = a + (b-1)*np + (c-1)*np*np
                         write(fid,'(i23,4x,i23,4x,e23.16)') &
                              ridx, cidx, du_vals(a,b,c,i,j,k) 
                         
                      end do
                   end do
                end do
                
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! print dv to file 
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'dv_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np
                         
                         ! local column index
                         cidx = a + (b-1)*np + (c-1)*np*np
                         write(fid,'(i23,4x,i23,4x,e23.16)') &
                              ridx, cidx, dv_vals(a,b,c,i,j,k) 
                         
                      end do
                   end do
                end do
                
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! print dps to file 
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'dps_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
                
                do b = 1,np
                   do a = 1,np
                      
                      ! local column index
                      cidx = a + (b-1)*np
                      write(fid,'(i23,4x,i23,4x,e23.16)') &
                           ridx, cidx, dps_vals(a,b,i,j,k) 
                      
                   end do
                end do
                
             end do
          end do
       end do

       close(fid)
      
    end do ! end element loop
#endif

#ifdef PRINT_FDJ 
    ! prints 1 column at a time

    ! starting and ending element index
    nets = fptr%nets
    nete = fptr%nete

    ! copy original state
    ystate = xstate
    
    do ii = eps_pow_start, eps_pow_stop

       ! perturbation for finite difference 
       eps = 10.0d0**ii

       ! initialize position in 1D state vector
       lx = 0
   
       ! ------------------------------------------------------------------------
       ! Compute Jacobian entries with respect to U velocity
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'du_FD_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if

       do ie = nets,nete
          do c = 1,nlev
             do b = 1,np
                do a = 1,np
                   
                   ! update position in 1D state vector
                   lx = lx+1
                   
                   ! compute finite difference approximation 
                   sigma = abs(xstate(lx)) * eps
                   
                   ystate(lx) = xstate(lx) + sigma
                   call residual_fn_aux(ystate, nxstate, fptr, ftens1)

                   ystate(lx) = xstate(lx) - sigma
                   call residual_fn_aux(ystate, nxstate, fptr, ftens2)

                   sigma = 2.0d0 * sigma

                   diff_ftens = (ftens1 - ftens2)/sigma

                   ! restore perturbed value
                   ystate(lx) = xstate(lx)

                   ! local column index
                   cidx = a + (b-1)*np + (c-1)*np*np
                   
                   do e = nets,nete
                      do k = 1,nlev
                         do j = 1,np
                            do i = 1,np
                               
                               ! local row index
                               ridx = i + (j-1)*np + (k-1)*np*np

                               write(fid,'(i23,4x,i23,4x,e23.16)') &
                                    ridx, cidx, diff_ftens(i,j,k,e) 
                               
                            end do
                         end do
                      end do
                   end do
                      
                end do
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! Compute Jacobian entries with respect to V velocity
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'dv_FD_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if

       do ie = nets,nete
          do c = 1,nlev
             do b = 1,np
                do a = 1,np

                   ! update position in 1D state vector
                   lx = lx+1

                   ! compute finite difference approximation 
                   sigma = abs(xstate(lx)) * eps

                   ystate(lx) = xstate(lx) + sigma
                   call residual_fn_aux(ystate, nxstate, fptr, ftens1)

                   ystate(lx) = xstate(lx) - sigma
                   call residual_fn_aux(ystate, nxstate, fptr, ftens2)

                   sigma = 2.0d0 * sigma

                   diff_ftens = (ftens1 - ftens2)/sigma

                   ! restore perturbed value
                   ystate(lx) = xstate(lx)

                   ! local column index
                   cidx = a + (b-1)*np + (c-1)*np*np
                
                   do e = nets,nete
                      do k = 1,nlev
                         do j = 1,np
                            do i = 1,np
                               
                               ! local row index
                               ridx = i + (j-1)*np + (k-1)*np*np

                               write(fid,'(i23,4x,i23,4x,e23.16)') &
                                    ridx, cidx, diff_ftens(i,j,k,e) 
                               
                            end do
                         end do
                      end do
                   end do

                end do
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! Compute Jacobian entries with respect to Ps
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
            'dps_FD_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if

       do ie = nets,nete
          do b = 1,np
             do a = 1,np
                
                ! update position in 1D state vector
                lx = lx+1
                
                ! compute finite difference approximation 
                sigma = abs(xstate(lx)) * eps
                
                ystate(lx) = xstate(lx) + sigma
                call residual_fn_aux(ystate, nxstate, fptr, ftens1)
                
                ystate(lx) = xstate(lx) - sigma
                call residual_fn_aux(ystate, nxstate, fptr, ftens2)
                
                sigma = 2.0d0 * sigma
                
                diff_ftens = (ftens1 - ftens2)/sigma
                
                ! restore perturbed value
                ystate(lx) = xstate(lx)

                ! local column index
                cidx = a + (b-1)*np

                do e = nets,nete
                   do k = 1,nlev
                      do j = 1,np
                         do i = 1,np
                            
                            ! local row index
                            ridx = i + (j-1)*np + (k-1)*np*np

                            write(fid,'(i23,4x,i23,4x,e23.16)') &
                                 ridx, cidx, diff_ftens(i,j,k,ie) 
                            
                         end do
                      end do
                   end do
                end do
                   
             end do
          end do
       end do

       close(fid)

    end do ! loop over epsilon
#endif
    
  end subroutine print_jacobian_aux



  subroutine print_jacobian_auxi(xstate, nxstate, fptr, residual_fn_auxi, &
       jacobian_fn_auxi, Jname)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    procedure(resid_auxi), pointer :: residual_fn_auxi ! Residual subroutine
    procedure(jac_auxi),   pointer :: jacobian_fn_auxi ! Jacobian subroutine

    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    character(len=1024) :: fname ! name of output file

    logical :: fexist      ! if output file exists
    integer :: fid         ! file ID for output

    real(c_double) :: ystate(nxstate)

    real(kind=real_kind), dimension(np,np,nlev+1,nelemd) :: ftens1, ftens2
    real(kind=real_kind), dimension(np,np,nlev+1,nelemd) :: diff_ftens

    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: du_vals, dv_vals
    real(kind=real_kind), dimension(np,np,np,np,nlev+1)      :: dps_vals

    real(kind=real_kind) :: eps   ! scaling for finite difference perturbation
    real(kind=real_kind) :: sigma ! finite difference perturbation

    integer :: ie, nets, nete, np1    ! element indices, time level
    integer :: i, j, k, ridx
    integer :: a, b, c, cidx
    integer :: ii, lx, e
    ! ---------------------------------------------------------------------------

#ifdef PRINT_AJ
    ! compute and output analytic Jacobian 
    np1 = fptr%tl%np1

    do ie = 1,nelemd
       
       call copy_state_to_HOMME(xstate, nxstate, fptr, np1)
       call jacobian_fn_auxi(du_vals, dv_vals, dps_vals, fptr, ie)

       ! ------------------------------------------------------------------------
       ! print to file 
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'du_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if
       
       do k = 1,nlev+1
          do j = 1,np
             do i = 1,np
                
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np
                         
                         ! local column index
                         cidx = a + (b-1)*np + (c-1)*np*np
                         write(fid,'(i23,4x,i23,4x,e23.16)') &
                              ridx, cidx, du_vals(a,b,c,i,j,k) 
                         
                      end do
                   end do
                end do
                
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! print to file 
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'dv_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if
       
       do k = 1,nlev+1
          do j = 1,np
             do i = 1,np
                
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np
                         
                         ! local column index
                         cidx = a + (b-1)*np + (c-1)*np*np
                         write(fid,'(i23,4x,i23,4x,e23.16)') &
                              ridx, cidx, dv_vals(a,b,c,i,j,k) 
                         
                      end do
                   end do
                end do
                
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! print to file 
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'dps_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if
       
       do k = 1,nlev+1
          do j = 1,np
             do i = 1,np
                
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
                
                do b = 1,np
                   do a = 1,np
                      
                      ! local column index
                      cidx = a + (b-1)*np
                      write(fid,'(i23,4x,i23,4x,e23.16)') &
                           ridx, cidx, dps_vals(a,b,i,j,k) 
                      
                   end do
                end do
                
             end do
          end do
       end do

       close(fid)
      
    end do ! end element loop
#endif

#ifdef PRINT_FDJ 
    ! prints 1 column at a time

    ! starting and ending element index
    nets = fptr%nets
    nete = fptr%nete

    ! copy original state
    ystate = xstate
    
    do ii = eps_pow_start, eps_pow_stop

       ! perturbation for finite difference 
       eps = 10.0d0**ii

       ! initialize position in 1D state vector
       lx = 0
   
       ! ------------------------------------------------------------------------
       ! Compute Jacobian entries with respect to U velocity
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'du_FD_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if

       do ie = nets,nete
          do c = 1,nlev
             do b = 1,np
                do a = 1,np
                   
                   ! update position in 1D state vector
                   lx = lx+1
                   
                   ! compute finite difference approximation 
                   sigma = abs(xstate(lx)) * eps
                   
                   ystate(lx) = xstate(lx) + sigma
                   call residual_fn_auxi(ystate, nxstate, fptr, ftens1)

                   ystate(lx) = xstate(lx) - sigma
                   call residual_fn_auxi(ystate, nxstate, fptr, ftens2)

                   sigma = 2.0d0 * sigma

                   diff_ftens = (ftens1 - ftens2)/sigma

                   ! restore perturbed value
                   ystate(lx) = xstate(lx)

                   ! local column index
                   cidx = a + (b-1)*np + (c-1)*np*np
                   
                   do e = nets,nete
                      do k = 1,nlev+1
                         do j = 1,np
                            do i = 1,np
                               
                               ! local row index
                               ridx = i + (j-1)*np + (k-1)*np*np
                               write(fid,'(i23,4x,i23,4x,e23.16)') &
                                    ridx, cidx, diff_ftens(i,j,k,e) 
                               
                            end do
                         end do
                      end do
                   end do
                      
                end do
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! Compute Jacobian entries with respect to V velocity
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'dv_FD_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if

       do ie = nets,nete
          do c = 1,nlev
             do b = 1,np
                do a = 1,np

                   ! update position in 1D state vector
                   lx = lx+1

                   ! compute finite difference approximation 
                   sigma = abs(xstate(lx)) * eps

                   ystate(lx) = xstate(lx) + sigma
                   call residual_fn_auxi(ystate, nxstate, fptr, ftens1)

                   ystate(lx) = xstate(lx) - sigma
                   call residual_fn_auxi(ystate, nxstate, fptr, ftens2)

                   sigma = 2.0d0 * sigma

                   diff_ftens = (ftens1 - ftens2)/sigma

                   ! restore perturbed value
                   ystate(lx) = xstate(lx)

                   ! local column index
                   cidx = a + (b-1)*np + (c-1)*np*np
                
                   do e = nets,nete
                      do k = 1,nlev+1
                         do j = 1,np
                            do i = 1,np
                               
                               ! local row index
                               ridx = i + (j-1)*np + (k-1)*np*np
                               write(fid,'(i23,4x,i23,4x,e23.16)') &
                                    ridx, cidx, diff_ftens(i,j,k,e) 
                               
                            end do
                         end do
                      end do
                   end do

                end do
             end do
          end do
       end do

       close(fid)

       ! ------------------------------------------------------------------------
       ! Compute Jacobian entries with respect to Ps
       ! ------------------------------------------------------------------------
       fid = 1000 + iam
       write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
            'dps_FD_',trim(Jname),'_proc',iam,'.txt'
       
       inquire(file=fname, exist=fexist)
       if (fexist) then
          open(fid, file=fname, form='formatted', status="old", &
               position="append", action="write")
       else
          open(fid, file=fname, form='formatted', status="new", action="write")
          write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
       end if

       do ie = nets,nete
          do b = 1,np
             do a = 1,np
                
                ! update position in 1D state vector
                lx = lx+1
                
                ! compute finite difference approximation 
                sigma = abs(xstate(lx)) * eps
                
                ystate(lx) = xstate(lx) + sigma
                call residual_fn_auxi(ystate, nxstate, fptr, ftens1)
                
                ystate(lx) = xstate(lx) - sigma
                call residual_fn_auxi(ystate, nxstate, fptr, ftens2)
                
                sigma = 2.0d0 * sigma
                
                diff_ftens = (ftens1 - ftens2)/sigma
                
                ! restore perturbed value
                ystate(lx) = xstate(lx)

                ! local column index
                cidx = a + (b-1)*np

                do e = nets,nete
                   do k = 1,nlev+1
                      do j = 1,np
                         do i = 1,np
                            
                            ! local row index
                            ridx = i + (j-1)*np + (k-1)*np*np
                            write(fid,'(i23,4x,i23,4x,e23.16)') &
                                 ridx, cidx, diff_ftens(i,j,k,ie) 
                            
                         end do
                      end do
                   end do
                end do
                   
             end do
          end do
       end do

       close(fid)

    end do ! loop over epsilon
#endif
    
  end subroutine print_jacobian_auxi



  ! =============================================================================
  ! Finite difference output subroutines called by print_jacobian
  ! =============================================================================



  subroutine print_du(Jname, fptr, a, b, c, ie, diff_utens, diff_vtens, &
       diff_ttens, diff_pstens)
    ! ---------------------------------------------------------------------------
    ! Outputs finite difference Jacobian to file
    ! ---------------------------------------------------------------------------
    use kinds, only                 : real_kind, long_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine

    type(derived_type), pointer :: fptr

    integer, intent(in) :: a, b, c, ie ! local derivative (column) index

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_utens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_vtens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_ttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(in) :: diff_pstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    character(len=1024) :: fname ! name of output file

    logical :: fexist      ! if output file exists
    integer :: fid         ! file ID for output
    
    integer :: i, j, k, e  ! equation loop
    integer :: ridx, cidx  ! row and column index
    ! ---------------------------------------------------------------------------

#ifdef LOCAL_IDX
    ! local cloumn index
    cidx = a + (b-1)*np + (c-1)*np*np
#else
    ! global column index
    cidx = fptr%base(ie)%gdofP(a,b) + Uoffset(c)
#endif

    ! ---------------------------------------------------------------------------
    ! A block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Ablock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dU_dU
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(e)%gdofP(i,j) + Uoffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_utens(i,j,k,e)

             end do
          end do
       end do
    end do

    ! dV_dU
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
      
#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index         
                ridx = fptr%base(e)%gdofP(i,j) + Voffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_vtens(i,j,k,e)
                
             end do
          end do
       end do
    end do
    
    close(fid)

    ! ---------------------------------------------------------------------------
    ! D Block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Dblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dPs_dU
    do e = 1,nelemd
       do j = 1,np
          do i = 1,np

#ifdef LOCAL_IDX
             ! local row index
             ridx = 2*np*np*nlev + i + (j-1)*np
#else
             ! global row index                                                    
             ridx = fptr%base(e)%gdofP(i,j) + Psoffset
#endif
             write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_pstens(i,j,e) 
             
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! F Block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Fblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dT_dU
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index                        
                ridx = fptr%base(e)%gdofP(i,j) + Toffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_ttens(i,j,k,e)
                
             end do
          end do
       end do
    end do

    close(fid)

  end subroutine print_du



  subroutine print_dv(Jname, fptr, a, b, c, ie, diff_utens, diff_vtens, &
       diff_ttens, diff_pstens)
    ! ---------------------------------------------------------------------------
    ! Outputs finite difference Jacobian to file
    ! ---------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine

    type(derived_type), pointer :: fptr

    integer, intent(in) :: a, b, c, ie

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_utens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_vtens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_ttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(in) :: diff_pstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    character(len=1024) :: fname ! name of output file

    logical :: fexist      ! if output file exists
    integer :: fid         ! file ID for output

    integer :: i, j, k, e  ! local equation (row) index 
    integer :: ridx, cidx  ! global row and column index
    ! ---------------------------------------------------------------------------

#ifdef LOCAL_IDX
    ! local cloumn index
    cidx = np*np*nlev + a + (b-1)*np + (c-1)*np*np
#else
    ! global column index
    cidx = fptr%base(ie)%gdofP(a,b) + Voffset(c)
#endif 

    ! ---------------------------------------------------------------------------
    ! A Block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Ablock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dU_dV
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(e)%gdofP(i,j) + Uoffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_utens(i,j,k,e)

             end do
          end do
       end do
    end do

    ! dV_dV
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index                        
                ridx = fptr%base(e)%gdofP(i,j) + Voffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_vtens(i,j,k,e)
                
             end do
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! D Block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Dblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dPs_dV
    do e = 1,nelemd
       do j = 1,np
          do i = 1,np

#ifdef LOCAL_IDX
             ! local row index
             ridx = 2*np*np*nlev + i + (j-1)*np
#else
             ! global row index
             ridx = fptr%base(e)%gdofP(i,j) + Psoffset
#endif
             write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_pstens(i,j,e) 
             
          end do
       end do
    end do
    
    close(fid)

    ! ---------------------------------------------------------------------------
    ! F Block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Fblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dT_dV
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index                                        
                ridx = fptr%base(e)%gdofP(i,j) + Toffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_ttens(i,j,k,e)
                
             end do
          end do
       end do
    end do

    close(fid)

  end subroutine print_dv



  subroutine print_dps(Jname, fptr, a, b, ie, diff_utens, diff_vtens, &
       diff_ttens, diff_pstens)
    ! ---------------------------------------------------------------------------
    ! Outputs finite difference Jacobian to file
    ! ---------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine

    type(derived_type), pointer :: fptr

    integer, intent(in) :: a, b, ie

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_utens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_vtens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_ttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(in) :: diff_pstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    character(len=1024) :: fname ! name of output file

    logical :: fexist      ! if output file exists
    integer :: fid         ! file ID for output

    integer :: i, j, k, e  ! equation loop
    integer :: ridx, cidx  ! row and column index
    ! ---------------------------------------------------------------------------

#ifdef LOCAL_IDX
    ! local cloumn index
    cidx = 2*np*np*nlev + a + (b-1)*np
#else
    ! global column index
    cidx = fptr%base(ie)%gdofP(a,b) + Psoffset
#endif

    ! ---------------------------------------------------------------------------
    ! B block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Bblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dU_dPs
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
           
#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index    
                ridx = fptr%base(e)%gdofP(i,j) + Uoffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_utens(i,j,k,e)

             end do
          end do
       end do
    end do

    ! dV_dPs
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
      
#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index                  
                ridx = fptr%base(e)%gdofP(i,j) + Voffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_vtens(i,j,k,e)
                
             end do
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! E block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Eblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dPs_dPs
    do e = 1,nelemd
       do j = 1,np
          do i = 1,np

#ifdef LOCAL_IDX
             ! local row index
             ridx = 2*np*np*nlev + i + (j-1)*np
#else
             ! global row index            
             ridx = fptr%base(e)%gdofP(i,j) + Psoffset
#endif
             write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_pstens(i,j,e) 
             
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! G block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Gblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dT_dPs
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(e)%gdofP(i,j) + Toffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_ttens(i,j,k,e)
                
             end do
          end do
       end do
    end do

    close(fid)


  end subroutine print_dps
  


  subroutine print_dT(Jname, fptr, a, b, c, ie, diff_utens, diff_vtens, &
       diff_ttens, diff_pstens)
    ! ---------------------------------------------------------------------------
    ! Outputs finite difference Jacobian to file
    ! ---------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine

    type(derived_type), pointer :: fptr

    integer, intent(in) :: a, b, c, ie

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_utens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_vtens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(in) :: diff_ttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(in) :: diff_pstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    character(len=1024) :: fname ! name of output file

    logical :: fexist      ! if output file exists
    integer :: fid         ! file ID for output

    integer :: i, j, k, e  ! equation loop
    integer :: ridx, cidx  ! row and column index
    ! ---------------------------------------------------------------------------

#ifdef LOCAL_IDX
    ! local cloumn index
    cidx = 2*np*np*nlev + np*np + a + (b-1)*np + (c-1)*np*np
#else
    ! global column index
    cidx = fptr%base(ie)%gdofP(a,b) + Toffset(c)
#endif

    ! ---------------------------------------------------------------------------
    ! C block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Cblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if

    ! dU_dT
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(e)%gdofP(i,j) + Uoffset(k)
#endif                
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_utens(i,j,k,e)

             end do
          end do
       end do
    end do

    ! dV_dT
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index                        
                ridx = fptr%base(e)%gdofP(i,j) + Voffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_vtens(i,j,k,e)
                
             end do
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! H block
    ! ---------------------------------------------------------------------------
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') &
         'Hblock_FD_',trim(Jname),'_proc',iam,'.txt'
    
    inquire(file=fname, exist=fexist)
    if (fexist) then
       open(fid, file=fname, form='formatted', status="old", &
            position="append", action="write")
    else
       open(fid, file=fname, form='formatted', status="new", action="write")
       write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    end if
    
    ! dT_dT
    do e = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index                                        
                ridx = fptr%base(e)%gdofP(i,j) + Toffset(k)
#endif
                write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_ttens(i,j,k,e)
                
             end do
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! write(fid,*)
    ! write(fid,'(a)') 'dPs_dT'

    ! do e = 1,nelemd
    !    do j = 1,np
    !       do i = 1,np
             
    !          ridx = fptr%base(e)%gdofP(i,j) + Psoffset
    !          write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, diff_pstens(i,j,e) 
             
    !       end do
    !    end do
    ! end do

  end subroutine print_dT

  
  
  ! =============================================================================
  ! Utility subroutines called by print_jacobian
  ! =============================================================================



  subroutine print_elem_data(fptr)
    ! ---------------------------------------------------------------------------
    ! Print values for data in HOMME element type
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
   
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    type(derived_type), pointer :: fptr
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: fid

    character(len=1024) :: fname

    integer :: a, b, c
    integer :: h, i, j, k, l
    integer :: ie, nets, nete
    integer :: np1

    real(kind=real_kind) :: temp1, temp2, temp3, val
    ! ---------------------------------------------------------------------------

    ! print global IDs
    call print_global_IDs(fptr)

    ! get start and end element indices, time level
    nets = fptr%nets
    nete = fptr%nete
    np1  = fptr%tl%np1

    ! assume each proc has one element in testing
    ie = 1 

    ! open output file
    fid = 1000 + iam
    write(fname,'(A,I2.2,A)') 'proc',iam,'.txt'      
    open(fid, file=fname, form='formatted')
    
    ! write out dimensions, parameter values
    write(fid,*) "nelem   = ", nelem, "nelemd = ",nelemd
    write(fid,*) "nets    = ", nete,  "nete = ", nets
    write(fid,*) "np      = ", np
    write(fid,*) "nlev    = ", nlev
    write(fid,*) "rrearth = ", rrearth

    ! write out mapping values for D
    write(fid,*)
    write(fid,*)
    do j = 1,np
       do i = 1,np
          write(fid,'(A2,I1,A1,I1,A4,F20.16,F20.16,A,9X,F20.16,F20.16,A)') &
               "D(",i,",",j,") = ", &
               fptr%base(ie)%D(i,j,1,1), fptr%base(ie)%D(i,j,1,2), &
               new_line('a'), &
               fptr%base(ie)%D(i,j,2,1), fptr%base(ie)%D(i,j,2,2), &
               new_line('a')
       end do
    end do

    ! write out mapping values for D^{-1}
    write(fid,*)
    write(fid,*)
    do j = 1,np
       do i = 1,np
          write(fid,'(A5,I1,A1,I1,A4,F20.16,F20.16,A,12X,F20.16,F20.16,A)') &
               "Dinv(",i,",",j,") = ", &
               fptr%base(ie)%Dinv(i,j,1,1), fptr%base(ie)%Dinv(i,j,1,2), &
               new_line('a'), &
               fptr%base(ie)%Dinv(i,j,2,1), fptr%base(ie)%Dinv(i,j,2,2), &
               new_line('a')
       end do
    end do

    ! write SEM derivative matrix values
    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "Dvv(i,j) = "
    do i = 1,np
       do j = 1,np
          write(fid,'(F20.16)',advance='no') fptr%deriv%Dvv(i,j)
       end do
       write(fid,*) 
    end do

    ! write out metric tensor
    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "metdet(i,j) = "
    do i = 1,np
       do j = 1,np
          write(fid,'(F20.16)',advance='no') fptr%base(ie)%metdet(i,j)
       end do
       write(fid,*) 
    end do

    ! write out 1/metric tensor
    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "rmetdet(i,j) = "
    do i = 1,np
       do j = 1,np
          write(fid,'(F20.16)',advance='no') fptr%base(ie)%rmetdet(i,j)
       end do
       write(fid,*) 
    end do

    ! write out layer interface weights
    write(fid,*)
    write(fid,*)   
    write(fid,'(A)') "Ai(k) and Bi(k) = "
    do k = 1,nlev+1
       write(fid,'(I4,2F20.16)') k, fptr%hvcoord%hyai(k), fptr%hvcoord%hybi(k)
    end do

    ! write out layer midpoint weights
    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "Am(k) and Bm(k) = "
    do k = 1,nlev
       write(fid,'(I4,2F20.16)') k, fptr%hvcoord%hyam(k), fptr%hvcoord%hybm(k)
    end do

    ! write out weights for projection operator (DSS)
    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "spheremp(i,j) = "
    do i = 1,np
       do j = 1,np
          write(fid,'(F20.16)',advance='no') fptr%base(ie)%spheremp(i,j)
       end do
       write(fid,*) 
    end do

    ! write out 1/weights for projection operator (DSS)
    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "rspheremp(i,j) = "
    do i = 1,np
       do j = 1,np
          write(fid,'(F20.16)',advance='no') fptr%base(ie)%rspheremp(i,j)
       end do
       write(fid,*) 
    end do

    ! write out global dof index
    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "gdofP(i,j) = "
    do i = 1,np
       do j = 1,np
          write(fid,'(I0)',advance='no') fptr%base(ie)%gdofp(i,j)
       end do
       write(fid,*) 
    end do

    write(fid,*)
    write(fid,*)
    do k = 1,nlev
       write(fid,*)
       write(fid,'(A,I0,A)') "U(i,j,",k,") = "   
       do i = 1,np
          do j = 1,np
             write(fid, '(E23.16)',advance='no') fptr%base(ie)%state%v(i,j,1,k,np1)
          end do
          write(fid,*)
       end do
       write(fid,*)
    end do

    write(fid,*)
    write(fid,*)
    do k = 1,nlev
       write(fid,*)
       write(fid,'(A,I0,A)') "V(i,j,",k,") = "   
       do i = 1,np
          do j = 1,np
             write(fid, '(E23.16)',advance='no') fptr%base(ie)%state%v(i,j,2,k,np1)
          end do
          write(fid,*)
       end do
       write(fid,*)
    end do

    write(fid,*)
    write(fid,*)
    do k = 1,nlev
       write(fid,*)
       write(fid,'(A,I0,A)') "T(i,j,",k,") = "   
       do i = 1,np
          do j = 1,np
             write(fid, '(E23.16)',advance='no') fptr%base(ie)%state%T(i,j,k,np1)
          end do
          write(fid,*)
       end do
       write(fid,*)
    end do

    write(fid,*)
    write(fid,*)
    write(fid,'(A)') "Ps(i,j) = "   
    do i = 1,np
       do j = 1,np
          write(fid, '(E23.16)',advance='no') fptr%base(ie)%state%ps_v(i,j,np1)
       end do
       write(fid,*)
    end do

    ! write reference pressure
    write(fid,*)
    write(fid,*) "ps0 = ",fptr%hvcoord%ps0   
    write(fid,*)

    close(fid)

    ! open output file
    write(fname,'(A,I2.2,A)') 'Element_',iam,'.txt'      
    open(fid, file=fname, form='formatted')

    ! write out global dof index
    write(fid,*)
    write(fid,'(A)') "global node index:"
    write(fid,*)
    do j = np,1,-1
       write(fid,'(I4.1,A5)',advance='no') j,'  |  '
       do i = 1,np
          write(fid,'(I5)',advance='no') fptr%base(ie)%gdofp(i,j)
       end do
       write(fid,*) 
    end do
    do i = 1,30
       write(fid,'(A1)',advance='no') '-'
    end do
    write(fid,*)
    write(fid,'(4X,A5)',advance='no') '  |  '
    do i = 1,np
       write(fid,'(I5)',advance='no') i
    end do
    write(fid,*)

    ! write out weights for projection operator (DSS)
    write(fid,*)
    write(fid,'(A)') "spheremp: "
    write(fid,*)
    do j = np,1,-1
       write(fid,'(I4.1,A5)',advance='no') j,'  |  '
       do i = 1,np
          write(fid,'(F21.16)',advance='no') fptr%base(ie)%spheremp(i,j)
       end do
       write(fid,*) 
    end do
    do i = 1,100
       write(fid,'(A1)',advance='no') '-'
    end do
    write(fid,*)
    write(fid,'(4X,A5)',advance='no') '  |  '
    do i = 1,np
       write(fid,'(10X,I1,10X)',advance='no') i
    end do
    write(fid,*)

    ! write out 1/weights for projection operator (DSS)
    write(fid,*)
    write(fid,'(A)') "rspheremp: "
    write(fid,*)
    do j = np,1,-1
       write(fid,'(I4.1,A5)',advance='no') j,'  |  '
       do i = 1,np
          write(fid,'(F21.16)',advance='no') fptr%base(ie)%rspheremp(i,j)
       end do
       write(fid,*) 
    end do
    do i = 1,100
       write(fid,'(A1)',advance='no') '-'
    end do
    write(fid,*)
    write(fid,'(4X,A5)',advance='no') '  |  '
    do i = 1,np
       write(fid,'(10X,I1,10X)',advance='no') i
    end do
    write(fid,*)

    ! write out weights for projection operator (DSS)
    write(fid,*)
    write(fid,'(A)') "spheremp * rspheremp: "
    write(fid,*)
    do j = np,1,-1
       write(fid,'(I4.1,A5)',advance='no') j,'  |  '
       do i = 1,np
          write(fid,'(F21.16)',advance='no') &
               fptr%base(ie)%spheremp(i,j) * fptr%base(ie)%rspheremp(i,j)
       end do
       write(fid,*) 
    end do
    do i = 1,100
       write(fid,'(A1)',advance='no') '-'
    end do
    write(fid,*)
    write(fid,'(4X,A5)',advance='no') '  |  '
    do i = 1,np
       write(fid,'(10X,I1,10X)',advance='no') i
    end do
    write(fid,*)

    ! write out weights for projection operator (DSS)
    write(fid,*)
    write(fid,'(A)') "1 / (spheremp * rspheremp): "
    write(fid,*)
    do j = np,1,-1
       write(fid,'(I4.1,A5)',advance='no') j,'  |  '
       do i = 1,np
          write(fid,'(F21.16)',advance='no') &
               1.0d0 / (fptr%base(ie)%spheremp(i,j) * fptr%base(ie)%rspheremp(i,j))
       end do
       write(fid,*) 
    end do
    do i = 1,100
       write(fid,'(A1)',advance='no') '-'
    end do
    write(fid,*)
    write(fid,'(4X,A5)',advance='no') '  |  '
    do i = 1,np
       write(fid,'(10X,I1,10X)',advance='no') i
    end do
    write(fid,*)

    close(fid)
    
  end subroutine print_elem_data



  ! =============================================================================
  ! The following subroutines implement each term in the residual equation to 
  ! to allow for verifying the corresponding Jacobian subroutines as well as the
  ! full residual without the DSS boundary exchange inorder to check local 
  ! Jacobian entries
  ! =============================================================================



  subroutine residual_full(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! compute the RHS, accumulate into a residual for each dependent variable 
    ! and apply DSS
    !
    ! For the Backward Euler first order option 
    !
    !           F(u)  = (u(np1) - u(n0))/dt  -  DSS[ RHS(u(np1)) ]
    !
    ! For the BDF2 second order option 
    !
    !           F(u)  = (u(np1) - u(n0))/dt  -  DSS[ RHS(u(np1)) ]
    !
    ! This subroutine is normally called to compute a leapfrog timestep
    ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
    ! accomodated.  For example, setting nm1=np1=n0 this routine will
    ! take a forward euler step, overwriting the input with the output.
    !
    !    qn0 = timelevel used to access Qdp() in order to compute 
    !          virtual Temperature
    !    qn0 = -1 for the dry case
    !
    ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
    !
    ! Combining the RHS and DSS pack operation in one routine
    ! allows us to fuse these two loops for more cache reuse
    !
    ! Combining the dt advance and DSS unpack operation in one routine
    ! allows us to fuse these two loops for more cache reuse
    !
    ! note: for prescribed velocity case, velocity will be computed at
    ! "real_time", which should be the time of timelevel n0.
    !
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use hybvcoord_mod, only         : hvcoord_t
    use hybrid_mod, only            : hybrid_t
    use element_mod, only           : element_t
    use physical_constants, only    : cp, cpwater_vapor, Rgas, kappa
    use prim_derived_type_mod, only : derived_type

    use control_mod, only           : moisture, qsplit, use_cpstar, rsplit
    use physics_mod, only           : virtual_specific_heat, virtual_temperature
    use prim_si_mod, only           : preq_vertadv, preq_omega_ps, preq_hydrostatic
    use derivative_mod, only        : derivative_t, gradient_sphere, divergence_sphere, vorticity_sphere

    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use time_mod, only              : TimeLevel_t

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: fvtens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real (kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real (kind=real_kind), pointer, dimension(:,:)   :: ps  ! surface pressure
    real (kind=real_kind), pointer, dimension(:,:,:) :: phi ! geopotential 

    real (kind=real_kind), dimension(np,np,nlev)   :: omega_p
    real (kind=real_kind), dimension(np,np,nlev)   :: T_v
    real (kind=real_kind), dimension(np,np,nlev)   :: divdp
    real (kind=real_kind), dimension(np,np,nlev+1) :: eta_dot_dpdn ! half level vertical velocity on p-grid
    real (kind=real_kind), dimension(np,np)        :: sdot_sum     ! temporary field
    real (kind=real_kind), dimension(np,np,2)      :: vtemp        ! generic gradient storage
    real (kind=real_kind), dimension(np,np)        :: vgrad_T      ! v.grad(T)
    real (kind=real_kind), dimension(np,np)        :: Ephi         ! kinetic energy + PHI term
    real (kind=real_kind), dimension(np,np,2)      :: grad_ps      ! lat-lon coord version
    real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p       
    real (kind=real_kind), dimension(np,np,nlev)   :: vort         ! vorticity
    real (kind=real_kind), dimension(np,np,nlev)   :: p            ! pressure
    real (kind=real_kind), dimension(np,np,nlev)   :: dp           ! delta pressure
    real (kind=real_kind), dimension(np,np,nlev)   :: rdp          ! inverse of delta pressure
    real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv       ! temperature vertical advection
    real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p      ! v.grad(p)
    real (kind=real_kind), dimension(np,np,nlev+1) :: ph           ! half level pressures on p-grid
    real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv       ! velocity vertical advection

    real (kind=real_kind) ::  kappa_star(np,np,nlev)
    real (kind=real_kind) ::  vtens1(np,np,nlev)
    real (kind=real_kind) ::  vtens2(np,np,nlev)
    real (kind=real_kind) ::  ttens(np,np,nlev)

    ! real (kind=real_kind) ::  fvtens(np,np,2,nlev,nelemd)
    ! real (kind=real_kind) ::  fttens(np,np,nlev,nelemd)
    ! real (kind=real_kind) ::  fdptens(np,np,nlev,nelemd)
    ! real (kind=real_kind) ::  fpstens(np,np,nelemd)

    real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
    real (kind=real_kind) ::  glnps1, glnps2, gpterm
    real (kind=real_kind) ::  gam
    integer :: i,j,k,kptr,ie, n_Q

    type (hvcoord_t)    :: hvcoord
    type (hybrid_t)     :: hybrid
    type (derivative_t) :: deriv
    type (TimeLevel_t)  :: tl

    integer :: method, nstep

    integer :: nets, nete, n0, np1, nm1, nm, qn0, lx, n
    real (kind=real_kind) :: dt2, dti
    logical  :: compute_diagnostics

    ! type(c_ptr)                        :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! set these to match the compute_and_apply_rhs for evenutal reincorporation, so
    ! BE    0 = x(np1)-x(nm1) - (RSS*x(n0))
    ! BDF2  0 = (1.5) (x(np1)-x(nm1)) - (0.5) (x(nm1)-x(nm)) - (RSS*x(n0))
    !       0 = (1.5) x(np1) - 2 x(nm1) + 0.5 x(nm) - (RSS*x(n0))

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1    

    nstep = fptr%tl%nstep 

    method     = fptr%method
    dt2        = fptr%dt
    dti        = 1.0d0/dt2
    hvcoord    = fptr%hvcoord
    hybrid     = fptr%hybrid
    deriv      = fptr%deriv
    compute_diagnostics = fptr%compute_diagnostics
    qn0        = fptr%n_Q

    ! unpack state vector
    lx = 1
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

    do ie=nets,nete
       do j=1,np
          do i=1,np
             fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
             lx = lx+1
          end do
       end do
    end do

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


    ! initialize residual
    futens  = 0.0d0
    fvtens  = 0.0d0
    fttens  = 0.0d0
    fpstens = 0.0d0

    ! ---------------------------------------------------------------------------
    ! compute residual
    ! ---------------------------------------------------------------------------

    if (method==11) then
       gam = 0.0d0  ! use BE for BDF as well for the first time step
    else if (method==12) then
       if (nstep==0) then
          gam = 0.0d0  ! BE
       else 
          gam = 0.5d0  ! BDF2 
       end if
    else 
       call abortmp('ERROR: method > 12 not yet coded for implicit time integration')
    end if

    do ie=nets,nete

       ! ==================================================
       ! compute pressure (p) on half levels from ps
       ! using the hybrid coordinates relationship, i.e.
       ! e.g. equation (3.a.92) of the CCM-2 description,
       ! (NCAR/TN-382+STR), June 1993, p. 24.
       ! ==================================================
       ! vertically eulerian only needs grad(ps)
       if (rsplit==0) &
            grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,n0),deriv,fptr%base(ie)%Dinv)

       ! ============================
       ! compute p and delta p
       ! ============================
       do k=1,nlev
          if (rsplit==0) then
             dp(:,:,k) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*fptr%base(ie)%state%ps_v(:,:,n0)) &
                  - (hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*fptr%base(ie)%state%ps_v(:,:,n0))
             p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*fptr%base(ie)%state%ps_v(:,:,n0)
             grad_p(:,:,:,k) = hvcoord%hybm(k)*grad_ps(:,:,:)

          else
             ! vertically lagrangian code: we advect dp3d instead of ps_v
             ! we also need grad(p) at all levels (not just grad(ps))
             !p(k)= hyam(k)*ps0 + hybm(k)*ps
             !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
             !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
             !
             ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
             !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
             dp(:,:,k) = fptr%base(ie)%state%dp3d(:,:,k,n0)
             if (k==1) then
                p(:,:,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp(:,:,k)/2
             else
                p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
             endif
             grad_p(:,:,:,k) = gradient_sphere(p(:,:,k),deriv,fptr%base(ie)%Dinv)
          endif

          rdp(:,:,k) = 1.0D0/dp(:,:,k)

          ! ============================
          ! compute vgrad_lnps
          ! ============================
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)

                vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
                vtemp(i,j,1) = v1*dp(i,j,k)
                vtemp(i,j,2) = v2*dp(i,j,k)
             end do
          end do

          ! =========================================
          !
          ! Compute relative vorticity and divergence
          !
          ! =========================================
          divdp(:,:,k)=divergence_sphere(vtemp,deriv,fptr%base(ie))
          vort(:,:,k)=vorticity_sphere(fptr%base(ie)%state%v(:,:,:,k,n0),deriv,fptr%base(ie))

       enddo

       ! compute T_v for timelevel n0
       !if ( moisture /= "dry") then
       if (qn0 == -1 ) then

          do k=1,nlev
             do j=1,np
                do i=1,np
                   T_v(i,j,k) = fptr%base(ie)%state%T(i,j,k,n0)
                   kappa_star(i,j,k) = kappa
                end do
             end do
          end do
       else

          do k=1,nlev
             do j=1,np
                do i=1,np
                   ! Qt = fptr%base(ie)%state%Q(i,j,k,1)
                   Qt = fptr%base(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
                   T_v(i,j,k) = Virtual_Temperature(fptr%base(ie)%state%T(i,j,k,n0),Qt)
                   if (use_cpstar==1) then
                      kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                   else
                      kappa_star(i,j,k) = kappa
                   endif
                end do
             end do
          end do
       end if

       ! ====================================================
       ! Compute Hydrostatic equation, modeld after CCM-3
       ! ====================================================
       !call geopotential_t(p,dp,T_v,Rgas,fptr%base(ie)%derived%phi)
       call preq_hydrostatic(fptr%base(ie)%derived%phi,fptr%base(ie)%state%phis,T_v,p,dp)

       ! ====================================================
       ! Compute omega_p according to CCM-3
       ! ====================================================
       call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)

       ! ==================================================
       ! zero partial sum for accumulating sum
       !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
       ! used by eta_dot_dpdn and lnps tendency
       ! ==================================================
       sdot_sum=0.0d0

       ! ==================================================
       ! Compute eta_dot_dpdn
       ! save sdot_sum as this is the -RHS of ps_v equation
       ! ==================================================
       if (rsplit>0) then
          ! VERTICALLY LAGRANGIAN:   no vertical motion
          eta_dot_dpdn=0.0d0
          T_vadv=0.0d0
          v_vadv=0.0d0
       else
          do k=1,nlev
             ! ==================================================
             ! add this term to PS equation so we exactly conserve dry mass
             ! ==================================================
             sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k)
             eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
          end do

          ! ===========================================================
          ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
          ! compute at interfaces:
          !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
          ! for reference: at mid layers we have:
          !    omega = v grad p  - integral_etatop^eta[ divdp ]
          ! ===========================================================
          do k=1,nlev-1
             eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
          end do

          eta_dot_dpdn(:,:,1     ) = 0.0D0
          eta_dot_dpdn(:,:,nlev+1) = 0.0D0

          ! ===========================================================
          ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
          ! ==============================================
          call preq_vertadv(fptr%base(ie)%state%T(:,:,:,n0),fptr%base(ie)%state%v(:,:,:,:,n0), &
               eta_dot_dpdn,rdp,T_vadv,v_vadv)
       endif

       ! ==============================================
       ! Compute phi + kinetic energy term: 10*np*np Flops
       ! ==============================================
       do k=1,nlev
          do j=1,np
             do i=1,np
                v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                E = 0.5D0*( v1*v1 + v2*v2 )
                Ephi(i,j)=E+fptr%base(ie)%derived%phi(i,j,k)+fptr%base(ie)%derived%pecnd(i,j,k)
             end do
          end do
          ! ================================================
          ! compute gradp term (ps/p)*(dp/dps)*T
          ! ================================================
          vtemp(:,:,:)   = gradient_sphere(fptr%base(ie)%state%T(:,:,k,n0),deriv,fptr%base(ie)%Dinv)
          do j=1,np
             do i=1,np
                v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
             end do
          end do


          ! vtemp = grad ( E + PHI )
          vtemp = gradient_sphere(Ephi(:,:),deriv,fptr%base(ie)%Dinv)

          do j=1,np
             do i=1,np

                gpterm = T_v(i,j,k)/p(i,j,k)
                glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
                glnps2 = Rgas*gpterm*grad_p(i,j,2,k)

                v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2     = fptr%base(ie)%state%v(i,j,2,k,n0)

                vtens1(i,j,k) =   - v_vadv(i,j,1,k) &
                     + v2*(fptr%base(ie)%fcor(i,j) + vort(i,j,k)) &
                     - vtemp(i,j,1) - glnps1

                vtens2(i,j,k) =   - v_vadv(i,j,2,k) &
                     - v1*(fptr%base(ie)%fcor(i,j) + vort(i,j,k)) &
                     - vtemp(i,j,2) - glnps2

                ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) &
                     + kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)

             end do
          end do

       end do

       ! Implicit Euler gam = 0.0
       ! BDF2           gam = 0.5
       do k=1,nlev

          fttens(:,:,k,ie) = (1.0d0 + gam) * (fptr%base(ie)%state%T(:,:,k,np1) &
               - fptr%base(ie)%state%T(:,:,k,nm1)) * dti - &
               gam * (fptr%base(ie)%state%T(:,:,k,nm1) &
               - fptr%base(ie)%state%T(:,:,k,nm))*dti - ttens(:,:,k)

          futens(:,:,k,ie) = (1.0d0 + gam) * (fptr%base(ie)%state%v(:,:,1,k,np1) &
               - fptr%base(ie)%state%v(:,:,1,k,nm1)) * dti - &
               gam * (fptr%base(ie)%state%v(:,:,1,k,nm1) &
               - fptr%base(ie)%state%v(:,:,1,k,nm)) * dti - vtens1(:,:,k)

          fvtens(:,:,k,ie) = (1.0d0 + gam) * (fptr%base(ie)%state%v(:,:,2,k,np1) &
               - fptr%base(ie)%state%v(:,:,2,k,nm1)) * dti - &
               gam * (fptr%base(ie)%state%v(:,:,2,k,nm1) - &
               fptr%base(ie)%state%v(:,:,2,k,nm))*dti - vtens2(:,:,k)
       enddo

       fpstens(:,:,ie) = (1.0d0 + gam) * (fptr%base(ie)%state%ps_v(:,:,np1) &
            - fptr%base(ie)%state%ps_v(:,:,nm1)) * dti - &
            gam * (fptr%base(ie)%state%ps_v(:,:,nm1) &
            - fptr%base(ie)%state%ps_v(:,:,nm)) * dti + sdot_sum(:,:)
    end do


#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_full


  
  subroutine residual_time(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: fvtens
    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind) :: dti, gam

    integer :: nets, nete, n0, np1, nm1, nm, nstep, method
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1    

    nstep  = fptr%tl%nstep 
    method = fptr%method

    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute time derivative contribution to residual
    ! ---------------------------------------------------------------------------

    dti = 1.0d0/fptr%dt

    if (method==11) then
       gam=0.0d0  ! use BE for BDF as well for the first time step
    else if (method==12) then
       if (nstep==0) then
          gam=0.0d0  ! BE
       else 
          gam=0.5d0  ! BDF2 
       end if
    else 
       call abortmp('ERROR: method > 12 not yet coded for implicit time integration')
    end if
    
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = ((1+gam) * (fptr%base(ie)%state%T(:,:,k,np1) &
               - fptr%base(ie)%state%T(:,:,k,nm1) )*dti - &
               gam*(fptr%base(ie)%state%T(:,:,k,nm1) &
               - fptr%base(ie)%state%T(:,:,k,nm))*dti )

          futens(:,:,k,ie) = ((1+gam) * (fptr%base(ie)%state%v(:,:,1,k,np1) &
               - fptr%base(ie)%state%v(:,:,1,k,nm1) )*dti - &
               gam*(fptr%base(ie)%state%v(:,:,1,k,nm1) &
               - fptr%base(ie)%state%v(:,:,1,k,nm))*dti )

          fvtens(:,:,k,ie) = ((1+gam) * (fptr%base(ie)%state%v(:,:,2,k,np1) &
               - fptr%base(ie)%state%v(:,:,2,k,nm1) )*dti - &
               gam*(fptr%base(ie)%state%v(:,:,2,k,nm1) &
               - fptr%base(ie)%state%v(:,:,2,k,nm))*dti )
       enddo
       fpstens(:,:,ie) = ((1+gam) * (fptr%base(ie)%state%ps_v(:,:,np1) &
            - fptr%base(ie)%state%ps_v(:,:,nm1) )*dti - &
            gam*(fptr%base(ie)%state%ps_v(:,:,nm1) &
            - fptr%base(ie)%state%ps_v(:,:,nm))*dti )
    end do  

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_time




  subroutine residual_vorticity(xstate, nxstate, fptr, futens, fvtens, fttens, &
       fpstens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : vorticity_sphere
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind), dimension(np,np,nlev) :: vort
    
    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from vorticity
    ! ---------------------------------------------------------------------------
    
    do ie = nets,nete
       do k = 1,nlev

          vort(:,:,k) = vorticity_sphere(fptr%base(ie)%state%v(:,:,:,k,n0), &
               fptr%deriv, fptr%base(ie))
          
          do j = 1,np
             do i = 1,np
                
                futens(i,j,k,ie) = -1.0d0 * fptr%base(ie)%state%v(i,j,2,k,n0)*vort(i,j,k)

                fvtens(i,j,k,ie) = fptr%base(ie)%state%v(i,j,1,k,n0)*vort(i,j,k)

             end do
          end do
       end do
    end do  

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_vorticity

  

  subroutine residual_coriolis(xstate, nxstate, fptr, futens, fvtens, fttens, &
       fpstens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod ,only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from coriolis
    ! ---------------------------------------------------------------------------
    
    do ie = nets,nete
       
       do k = 1,nlev          
          do j = 1,np
             do i = 1,np
                
                futens(i,j,k,ie) = -1.0d0 * fptr%base(ie)%state%v(i,j,2,k,n0)*fptr%base(ie)%fcor(i,j)

                fvtens(i,j,k,ie) = fptr%base(ie)%state%v(i,j,1,k,n0)*fptr%base(ie)%fcor(i,j)
                
             end do
          end do
       end do
    end do  

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_coriolis



  subroutine residual_Kinetic_Energy(xstate, nxstate, fptr, futens, fvtens, &
       fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : gradient_sphere
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------  
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind) :: v1, v2, E, Ephi(np,np), vtemp(np,np,2)
 
    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from kinetic energy
    ! ---------------------------------------------------------------------------
    do ie = nets,nete
       do k = 1,nlev          
          do j = 1,np
             do i = 1,np
                
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)
                E  = 0.5D0*( v1*v1 + v2*v2 )
                Ephi(i,j) = E 
                
             end do
          end do

          vtemp = gradient_sphere(Ephi(:,:), fptr%deriv, fptr%base(ie)%Dinv)

          futens(:,:,k,ie) = vtemp(:,:,1)

          fvtens(:,:,k,ie) = vtemp(:,:,2)

       end do
    end do  

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_Kinetic_Energy



  subroutine residual_Geopotential(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : gradient_sphere
    use prim_si_mod, only           : preq_hydrostatic
    use physical_constants, only    : kappa, Rgas
    use physics_mod, only           : virtual_specific_heat, virtual_temperature
    use control_mod, only           : use_cpstar
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------  
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind), dimension(np,np,nlev) :: p, dp, T_v, kappa_star
    real(kind=real_kind), dimension(np,np,2)    :: vtemp
    real(kind=real_kind), dimension(np,np)      :: Ephi
    real(kind=real_kind) :: Qt
 
    integer :: nets, nete, n0, np1, nm1, nm, qn0
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1

    qn0 = fptr%n_Q
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from kinetic energy
    ! ---------------------------------------------------------------------------
    do ie = nets,nete
       do k = 1,nlev          
          
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,n0)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,n0))
          
          p(:,:,k) = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,n0)

       end do
       
       ! compute T_v for timelevel n0
       !if ( moisture /= "dry") then
       if ( qn0 == -1 ) then
          do k=1,nlev
             do j=1,np
                do i=1,np
                   T_v(i,j,k) = fptr%base(ie)%state%T(i,j,k,n0)
                   kappa_star(i,j,k) = kappa
                end do
             end do
          end do
       else
          do k=1,nlev
             do j=1,np
                do i=1,np
                   ! Qt = fptr%base(ie)%state%Q(i,j,k,1)
                   Qt = fptr%base(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
                   T_v(i,j,k) = Virtual_Temperature(fptr%base(ie)%state%T(i,j,k,n0),Qt)
                   if (use_cpstar==1) then
                      kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                   else
                      kappa_star(i,j,k) = kappa
                   endif
                end do
             end do
          end do
       end if
       
       call preq_hydrostatic(fptr%base(ie)%derived%phi,fptr%base(ie)%state%phis,T_v,p,dp)
       
       do k = 1,nlev          
          do j = 1,np
             do i = 1,np
                
                Ephi(i,j) = fptr%base(ie)%derived%phi(i,j,k)
                
             end do
          end do

          vtemp = gradient_sphere(Ephi(:,:), fptr%deriv, fptr%base(ie)%Dinv)

          futens(:,:,k,ie) = vtemp(:,:,1)

          fvtens(:,:,k,ie) = vtemp(:,:,2)

       end do
    end do  

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_Geopotential



  subroutine residual_pressure_grad(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Compute the pressure gradient part of the residual 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : gradient_sphere
    use prim_si_mod, only           : preq_hydrostatic
    use physical_constants, only    : kappa, Rgas
    use physics_mod, only           : virtual_specific_heat, virtual_temperature
    use control_mod, only           : use_cpstar
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------  
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p
    real(kind=real_kind), dimension(np,np,nlev)   :: p, dp, T_v, kappa_star
    real(kind=real_kind), dimension(np,np,2)      :: grad_ps
    real(kind=real_kind) :: Qt, gpterm, glnps1, glnps2
 
    integer :: nets, nete, n0, np1, nm1, nm, qn0
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    qn0 = fptr%n_Q
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from pressure gradient
    ! ---------------------------------------------------------------------------
    
    do ie = nets,nete

       grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,n0), fptr%deriv, fptr%base(ie)%Dinv)

       do k = 1,nlev          
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,n0)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,n0))
                    
          p(:,:,k) = fptr%hvcoord%hyam(k)*fptr%hvcoord%ps0 + fptr%hvcoord%hybm(k)*fptr%base(ie)%state%ps_v(:,:,n0)

          grad_p(:,:,:,k) = fptr%hvcoord%hybm(k)*grad_ps(:,:,:)

       end do
       
       ! compute T_v for timelevel n0
       !if ( moisture /= "dry") then
       if (qn0 == -1 ) then
          do k=1,nlev
             do j=1,np
                do i=1,np
                   T_v(i,j,k) = fptr%base(ie)%state%T(i,j,k,n0)
                   kappa_star(i,j,k) = kappa
                end do
             end do
          end do
       else
          do k=1,nlev
             do j=1,np
                do i=1,np
                   ! Qt = fptr%base(ie)%state%Q(i,j,k,1)
                   Qt = fptr%base(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
                   T_v(i,j,k) = Virtual_Temperature(fptr%base(ie)%state%T(i,j,k,n0),Qt)
                   if (use_cpstar==1) then
                      kappa_star(i,j,k) = Rgas / Virtual_Specific_Heat(Qt)
                   else
                      kappa_star(i,j,k) = kappa
                   endif
                end do
             end do
          end do
       end if
              
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
       
                gpterm = T_v(i,j,k) / p(i,j,k)
                glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
                glnps2 = Rgas*gpterm*grad_p(i,j,2,k)
        
                futens(i,j,k,ie) = glnps1

                fvtens(i,j,k,ie) = glnps2
                
             end do
          end do
       end do

    end do

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_pressure_grad




  subroutine residual_VertAdv(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Computes the vertical advection parts of the residual 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : divergence_sphere
    use prim_si_mod, only           : preq_vertadv
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate) 
    integer(c_int), intent(in), value :: nxstate         
    type(derived_type), pointer       :: fptr            

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------  
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind), dimension(np,np,nlev+1) :: eta_dot_dpdn
    real(kind=real_kind), dimension(np,np,nlev)   :: dp, rdp, divdp
    real(kind=real_kind), dimension(np,np,2)      :: vtemp
    real(kind=real_kind), dimension(np,np)        :: sdot_sum
    real(kind=real_kind) :: v1, v2

    real(kind=real_kind), dimension(np,np,2,nlev) :: v_vadv
    real(kind=real_kind), dimension(np,np,nlev)   :: T_vadv

    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,ie,lx
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from vertical advection
    ! ---------------------------------------------------------------------------
    do ie = nets,nete
    
       do k=1,nlev

          ! compute dp/dn
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,n0)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,n0))

          ! compute 1/(dp/dn)
          rdp(:,:,k) = 1.0d0 / dp(:,:,k)

          ! compute (1/(dp/dn)) * v
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)

                vtemp(i,j,1) = v1*dp(i,j,k)
                vtemp(i,j,2) = v2*dp(i,j,k)
             end do
          end do

          ! compute divergence( (1/(dp/dn)) * v )
          divdp(:,:,k) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))
       end do

       ! compute vertical integrals necessary to form \dot{n} * dp/dn
       sdot_sum = 0.0d0
       do k=1,nlev
          sdot_sum(:,:)         = sdot_sum(:,:) + divdp(:,:,k)
          eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
       end do

       ! compute \dot{n} * dp/dn
       do k=1,nlev-1
          eta_dot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1) * sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
       end do

       eta_dot_dpdn(:,:,1     ) = 0.0d0
       eta_dot_dpdn(:,:,nlev+1) = 0.0d0

       ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
       call preq_vertadv(fptr%base(ie)%state%T(:,:,:,n0), fptr%base(ie)%state%v(:,:,:,:,n0), &
            eta_dot_dpdn, rdp, T_vadv, v_vadv)

       do k=1,nlev
          do j=1,np
             do i=1,np
                futens(i,j,k,ie) = v_vadv(i,j,1,k)

                fvtens(i,j,k,ie) = v_vadv(i,j,2,k)

                fttens(i,j,k,ie) = T_vadv(i,j,k)
             end do
          end do
       end do

    end do ! element loop 

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_VertAdv
  
  

  subroutine residual_HorizAdv_Temp(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Computes the horizontal advection of temperature part of the residual
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod ,only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : gradient_sphere
        
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace------------------------------- 
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind), dimension(np,np,2) :: vtemp
    real(kind=real_kind), dimension(np,np)   :: vgrad_T
    real(kind=real_kind) :: v1, v2

    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from horizontal temperature advection
    ! ---------------------------------------------------------------------------
    do ie = nets,nets
       do k = 1,nlev
          
          ! gradient of temperature
          vtemp = gradient_sphere(fptr%base(ie)%state%T(:,:,k,n0), fptr%deriv, fptr%base(ie)%Dinv)

          ! v dot grad(T)
          do j = 1,np
             do i = 1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)

                vgrad_T(i,j) = v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
             end do
          end do

          ! temperature tendency due to horizontal temperature advection
          do j = 1,np
             do i = 1,np
                fttens(i,j,k,ie) = vgrad_T(i,j)
             end do
          end do

       end do
    end do

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_HorizAdv_Temp



  subroutine residual_pressure_vert_vel(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : gradient_sphere, divergence_sphere
    use physical_constants, only    : kappa, Rgas
    use physics_mod, only           : virtual_specific_heat, virtual_temperature
    use control_mod, only           : use_cpstar
    use prim_si_mod, only           : preq_omega_ps
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------  
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p
    real(kind=real_kind), dimension(np,np,nlev)   :: p, dp, vgrad_p, divdp
    real(kind=real_kind), dimension(np,np,nlev)   :: T_v, kappa_star
    real(kind=real_kind), dimension(np,np,nlev)   :: omega_p
    real(kind=real_kind), dimension(np,np,2)      :: grad_ps, vtemp
    real(kind=real_kind) :: v1, v2, Qt
 
    integer :: nets, nete, n0, np1, nm1, nm, qn0
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1

    qn0 = fptr%n_Q
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from pressure vertical velocity
    ! ---------------------------------------------------------------------------  
    do ie = nets,nete

       grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,n0), fptr%deriv, fptr%base(ie)%Dinv)

       do k = 1,nlev          

          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,n0)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,n0))
                    
          p(:,:,k) = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,n0)

          grad_p(:,:,:,k) = fptr%hvcoord%hybm(k)*grad_ps(:,:,:)

          ! compute v dot grad(p) and (1/(dp/dn)) * v
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)

                vgrad_p(i,j,k) = v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k)

                vtemp(i,j,1) = v1*dp(i,j,k)
                vtemp(i,j,2) = v2*dp(i,j,k)
             end do
          end do

          ! compute divergence( (1/(dp/dn)) * v )
          divdp(:,:,k) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))

       end do

       ! compute T_v for timelevel n0
       !if ( moisture /= "dry") then
       if (qn0 == -1 ) then

          do k=1,nlev
             do j=1,np
                do i=1,np

                   T_v(i,j,k)        = fptr%base(ie)%state%T(i,j,k,n0)
                   kappa_star(i,j,k) = kappa

                end do
             end do
          end do

       else

          do k=1,nlev
             do j=1,np
                do i=1,np

                   Qt         = fptr%base(ie)%state%Qdp(i,j,k,1,qn0) / dp(i,j,k)
                   T_v(i,j,k) = Virtual_Temperature(fptr%base(ie)%state%T(i,j,k,n0), Qt)

                   if (use_cpstar==1) then
                      kappa_star(i,j,k) = Rgas / Virtual_Specific_Heat(Qt)
                   else
                      kappa_star(i,j,k) = kappa

                   endif
                end do
             end do
          end do

       end if

       ! compute pressure vertical velocity (omega_p)
       call preq_omega_ps(omega_p, fptr%hvcoord, p, vgrad_p, divdp)

       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fttens(i,j,k,ie) = -1.0d0 * kappa_star(i,j,k) * T_v(i,j,k) * omega_p(i,j,k)
             end do
          end do
       end do

    end do

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_pressure_vert_vel
  
  
  
  subroutine residual_int_div_dpdn_vel(xstate, nxstate, fptr, futens, fvtens, fttens, fpstens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : divergence_sphere
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: futens, fvtens, fttens
    real(kind=real_kind), dimension(np,np,nelemd),      intent(out) :: fpstens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(kind=real_kind) :: fvtens_vec(np,np,2,nlev,nelemd)

    real(kind=real_kind), dimension(np,np,nlev) :: dp, divdp
    real(kind=real_kind), dimension(np,np,2)    :: vtemp
    real(kind=real_kind), dimension(np,np)      :: sdot_sum
    real(kind=real_kind) :: v1, v2
 
    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    integer :: kptr
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize residual
    futens  = 0.0d0    
    fvtens  = 0.0d0
    fttens  = 0.0d0    
    fpstens = 0.0d0   

    ! ---------------------------------------------------------------------------
    ! compute residual contribution from pressure vertical velocity
    ! ---------------------------------------------------------------------------
    do ie = nets,nete
       
       do k=1,nlev

          ! compute dp/dn
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,n0)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,n0))

          ! compute (1/(dp/dn)) * v
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)
                
                vtemp(i,j,1) = v1*dp(i,j,k)
                vtemp(i,j,2) = v2*dp(i,j,k)
             end do
          end do
          
          ! compute divergence( (1/(dp/dn)) * v )
          divdp(:,:,k) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))
       end do
       
       ! compute vertical integral
       sdot_sum = 0.0d0
       do k=1,nlev
          sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k)
       end do
       
       do j=1,np
          do i=1,np
             fpstens(i,j,ie) = sdot_sum(i,j)
          end do
       end do

    end do

#ifdef DSS_ON
    do ie = nets,nete
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:) * fttens(:,:,k,ie) 
          fvtens_vec(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:) * futens(:,:,k,ie)
          fvtens_vec(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:) * fvtens(:,:,k,ie)
       enddo

       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:) * fpstens(:,:,ie)
       
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)
       
       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)
       
       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens_vec(:,:,:,:,ie),2*nlev,kptr,ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    do ie=nets,nete

       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)
       
       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)
       
       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens_vec(:,:,:,:,ie), 2*nlev, kptr, ie)

       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          futens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,1,k,ie)
          fvtens(:,:,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens_vec(:,:,2,k,ie)
       end do
       
       fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
    end do
#endif

  end subroutine residual_int_div_dpdn_vel



  subroutine residual_aux_vel_grad_p(xstate, nxstate, fptr, ftens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : gradient_sphere
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: ftens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p
    real(kind=real_kind), dimension(np,np,2)      :: grad_ps
    real(kind=real_kind), dimension(np,np,nlev)   :: vgrad_p
    real(kind=real_kind), dimension(np,np,2)      :: vtemp
    real(kind=real_kind) :: v1, v2
 
    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize output
    ftens = 0.0d0    

    do ie = nets,nete
       
       grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,n0), fptr%deriv, fptr%base(ie)%Dinv)
       
       do k=1,nlev
          grad_p(:,:,:,k) = fptr%hvcoord%hybm(k)*grad_ps(:,:,:)
       
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)
                
                vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
             end do
          end do
       
       end do

       do k = 1,nlev            
          do j=1,np
             do i=1,np
                ftens(i,j,k,ie) = vgrad_p(i,j,k)
             end do
          end do
       end do

    end do

  end subroutine residual_aux_vel_grad_p



  subroutine residual_aux_div_dpdn_vel(xstate, nxstate, fptr, ftens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : divergence_sphere
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev,nelemd), intent(out) :: ftens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev) :: dp, divdp
    real(kind=real_kind), dimension(np,np,2)    :: vtemp
    real(kind=real_kind) :: v1, v2
 
    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize output
    ftens = 0.0d0    

    do ie = nets,nete
       
       do k=1,nlev

          ! compute dp/dn
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,n0)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,n0))

          ! compute dp/dn * v
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)
                
                vtemp(i,j,1) = v1*dp(i,j,k)
                vtemp(i,j,2) = v2*dp(i,j,k)
             end do
          end do
          
          ! compute divergence(dp/dn * v)
          divdp(:,:,k) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))
              
          do j=1,np
             do i=1,np
                ftens(i,j,k,ie) = divdp(i,j,k)
             end do
          end do
          
       end do

    end do

  end subroutine residual_aux_div_dpdn_vel



  subroutine residual_aux_etadot_dpdn(xstate, nxstate, fptr, ftens)
    ! ---------------------------------------------------------------------------
    ! Checks the elements of the analytic Jacobian matrix 
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    use perf_mod, only              : t_startf, t_stopf
    use derivative_mod, only        : divergence_sphere
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr

    real(kind=real_kind), dimension(np,np,nlev+1,nelemd), intent(out) :: ftens
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev+1) :: eta_dot_dpdn
    real(kind=real_kind), dimension(np,np,nlev)   :: dp, divdp
    real(kind=real_kind), dimension(np,np,2)      :: vtemp
    real(kind=real_kind), dimension(np,np)        :: sdot_sum
    real(kind=real_kind) :: v1, v2
 
    integer :: nets, nete, n0, np1, nm1, nm
    integer :: i,j,k,lx,ie
    ! ---------------------------------------------------------------------------

    ! unpack element indices
    nets = fptr%nets
    nete = fptr%nete

    ! unpack time levels
    np1 = fptr%tl%np1  
    nm1 = fptr%tl%n0    
    n0  = fptr%tl%np1    
    nm  = fptr%tl%nm1
    
    call copy_state_to_HOMME(xstate, nxstate, fptr, np1)

    ! initialize output
    ftens = 0.0d0    

    do ie = nets,nete
       
       do k=1,nlev

          ! compute dp/dn
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,n0)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,n0))

          ! compute (1/(dp/dn)) * v
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)
                
                vtemp(i,j,1) = v1*dp(i,j,k)
                vtemp(i,j,2) = v2*dp(i,j,k)
             end do
          end do
          
          ! compute divergence( dp/dn * v )
          divdp(:,:,k) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))
       end do

       ! compute vertical integral
       sdot_sum = 0.0d0

       do k=1,nlev
          sdot_sum(:,:)         = sdot_sum(:,:) + divdp(:,:,k)
          eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
       end do

       do k=1,nlev-1
          eta_dot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1)*sdot_sum(:,:) &
               - eta_dot_dpdn(:,:,k+1)
       end do

       eta_dot_dpdn(:,:,1     ) = 0.0D0
       eta_dot_dpdn(:,:,nlev+1) = 0.0D0

       do k = 1,nlev+1
          do j = 1,np
             do i = 1,np
                ftens(i,j,k,ie) = eta_dot_dpdn(i,j,k)
             end do
          end do
       end do

    end do

  end subroutine residual_aux_etadot_dpdn



  subroutine copy_state_from_HOMME(xstate, nxstate, fptr, idx)
    ! ---------------------------------------------------------------------------
    ! Copys array xstate into Fortran pointer for HOMME derived type
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(out)       :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr
    integer, intent(in)               :: idx
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: nets, nete
    integer :: i, j, k, ie, lx
    ! ---------------------------------------------------------------------------
    
    nets = fptr%nets
    nete = fptr%nete

    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                xstate(lx) = fptr%base(ie)%state%v(i,j,1,k,idx)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                xstate(lx) = fptr%base(ie)%state%v(i,j,2,k,idx)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             xstate(lx) = fptr%base(ie)%state%ps_v(i,j,idx)
             lx = lx+1
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                xstate(lx) = fptr%base(ie)%state%T(i,j,k,idx)
                lx = lx+1
             end do
          end do
       end do
    end do

  end subroutine copy_state_from_HOMME



  subroutine copy_state_to_HOMME(xstate, nxstate, fptr, idx)
    ! ---------------------------------------------------------------------------
    ! Copys array xstate into Fortran pointer for HOMME derived type
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelem, nelemd
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(derived_type), pointer       :: fptr
    integer, intent(in)               :: idx
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: nets, nete
    integer :: i, j, k, ie, lx
    ! ---------------------------------------------------------------------------
    
    nets = fptr%nets
    nete = fptr%nete

    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,1,k,idx) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,2,k,idx) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fptr%base(ie)%state%ps_v(i,j,idx) = xstate(lx)
             lx = lx+1
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%T(i,j,k,idx) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do

  end subroutine copy_state_to_HOMME

  
end module prim_jacobian_check_mod

#endif
