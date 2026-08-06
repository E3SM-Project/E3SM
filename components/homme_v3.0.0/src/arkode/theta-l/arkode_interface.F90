!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
!
! This file contains all required subroutines for interfacing to
! the ARKode Fortran interface, as well as a small number of
! auxiliary routines to assist.  This uses a custom NVector
! interface, implemented in Fortran.
!
! The following subroutines are 'users' of the ARKode Fortran
! interface:
!    arkode_init -- initializes NVectors and ARKode for subsequent
!                   solves (only called once)
!
! The following user-supplied subroutines are required by the
! ARKode Fortran interface:
!    farkifun -- implements all implicit portions of the ODE
!                right-hand side function
!    farkefun -- implements all explicit portions of the ODE
!                right-hand side function
!
! The following user-supplied subroutines are optional within
! the ARKode Fortran interface:
!    FColumnSolSolve -- routine to perform a Fortran-supplied,
!                       columnwise linear solver
!    farkjtsetup -- prepares for Jacobian-vector products, J*v,
!                   where J(u) is the Jacobian of farkifun
!                   with respect to u.
!    farkjtimes -- implements the Jacobian-vector product, J*v,
!                  where J(u) is the Jacobian of farkifun
!                  with respect to u.
!    farkpset -- performs setup for the preconditioner matrix
!                (called infrequently)
!    farkpsol -- performs the preconditioner solve (called every
!                linear iteration)
!
! The following are 'helper' subroutines that I often use when
! interfacing to SUNDIALS solvers from Fortran:
!    farkdiags -- writes ARKode solver diagnostics to the screen
!
! Copyright 2017; all rights reserved
!
! Modified by Christopher J. Vogl (LLNL) 2020
!=================================================================

function farkifun(t, y, ydot, hvcoord, hybrid, deriv, imex, dt, eta_ave_w) &
  result(ierr)
  !-----------------------------------------------------------------
  ! Description: farkifun provides the implicit portion of the right
  !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
  !
  ! Arguments:
  !       t - (dbl, input) current time
  !       y - (NVec_t) NVec_t containing current solution
  !    ydot - (NVec_t) NVec_t to hold right-hand side function
  !-----------------------------------------------------------------

  !======= Inclusions ===========
  use fsundials_nvector_mod, only: N_Vector
  use kinds,                 only: real_kind
  use HommeNVector,          only: NVec_t
  use hybrid_mod,            only: hybrid_t
  use derivative_mod,        only: derivative_t
  use dimensions_mod,        only: nlevp
  use hybvcoord_mod,         only: hvcoord_t
  use prim_advance_mod,      only: compute_andor_apply_rhs

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,             intent(in)    :: t
  type(NVec_t),       intent(in)    :: y
  type(NVec_t),       intent(inout) :: ydot
  type(hvcoord_t),    intent(in)    :: hvcoord
  type(hybrid_t),     intent(in)    :: hybrid
  type(derivative_t), intent(in)    :: deriv
  integer,            intent(in)    :: imex
  real(real_kind),    intent(in)    :: dt
  real(real_kind),    intent(in)    :: eta_ave_w
  integer                           :: ierr

  ! local variables
  real(real_kind) :: bval, cval, scale1, scale2
  integer         :: ie

  !======= Internals ============

  ! set scale factors depending on whether using implicit, explicit, or IMEX
  if (imex == 0) then
    scale1 = 1.d0
    scale2 = 1.d0
  else if (imex == 1) then
    scale1 = 0.d0
    scale2 = 0.d0
  else if (imex == 2) then
    scale1 = 0.d0
    scale2 = 1.d0
  end if

  ! set return value to success
  ierr = 0

  ! TODO: obtain b value for current 'stage number' (only affect tracers)
  bval = 1.d0

  ! The function call to compute_andor_apply_rhs is as follows:
  !  compute_andor_apply_rhs(np1, nm1, n0, dt2, elem, hvcoord, hybrid, &
  !     deriv, nets, nete, compute_diagnostics, eta_ave_w, scale1, scale2, scale3)
  !
  !  This call returns the following:
  !
  !   u(np1) = scale3*u(nm1) + dt2*DSS[ nonstiffRHS(u(n0))*scale1 + stiffRHS(un0)*scale2 ]
  !
  !   nonstiffRHS and the stiffRHS are determined within the function and can be change by
  !   multiplying different terms by scale1 and scale2
  !
  !  Setting scale1=scale2=1.0, scale3=0.0, and dt2=1.0 returns the full rhs
  !
  !  DSS is the averaging procedure for the active and inactive nodes
  !

  call compute_andor_apply_rhs(ydot%tl_idx, ydot%tl_idx, y%tl_idx, &
        1.d0, y%elem, hvcoord, hybrid, deriv, y%nets, y%nete, &
        .false., bval*eta_ave_w, scale1, scale2, 0.d0)

  ! Zero out RHS vector elements to maintain (d/dt)phinh_i(nlevp) = 0
  do ie=ydot%nets,ydot%nete
    ydot%elem(ie)%state%phinh_i(:,:,nlevp,ydot%tl_idx) = 0.d0
  end do

  return
end function farkifun

!=================================================================

function farkefun(t, y, ydot, hvcoord, hybrid, deriv, imex, dt, eta_ave_w) &
  result(ierr)
  !-----------------------------------------------------------------
  ! Description: farkefun provides the explicit portion of the right
  !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
  !
  ! Arguments:
  !       t - (dbl, input) current time
  !       y - (NVec_t) NVec_t containing current solution
  !    ydot - (NVec_t) NVec_t to hold right-hand side function
  !-----------------------------------------------------------------

  !======= Inclusions ===========
  use fsundials_nvector_mod, only: N_Vector
  use kinds,                 only: real_kind
  use HommeNVector,          only: NVec_t
  use hybrid_mod,            only: hybrid_t
  use derivative_mod,        only: derivative_t
  use dimensions_mod,        only: nlevp
  use hybvcoord_mod,         only: hvcoord_t
  use prim_advance_mod,      only: compute_andor_apply_rhs

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,             intent(in)    :: t
  type(NVec_t),       intent(in)    :: y
  type(NVec_t),       intent(inout) :: ydot
  type(hvcoord_t),    intent(in)    :: hvcoord
  type(hybrid_t),     intent(in)    :: hybrid
  type(derivative_t), intent(in)    :: deriv
  integer,            intent(in)    :: imex
  real(real_kind),    intent(in)    :: dt
  real(real_kind),    intent(in)    :: eta_ave_w
  integer                           :: ierr

  ! local variables
  real(real_kind) :: bval, cval, scale1, scale2
  integer         :: ie

  !======= Internals ============

  ! set scale factors depending on whether using implicit, explicit, or IMEX
  if (imex == 0) then
    scale1 = 0.d0
    scale2 = 0.d0
  else if (imex == 1) then
    scale1 = 1.d0
    scale2 = 1.d0
  else if (imex == 2) then
    scale1 = 1.d0
    scale2 = 0.d0
  end if

  ! set return value to success
  ierr = 0

  ! TODO: obtain b value for current 'stage number' (only affects tracers)
  bval = 1.d0

  ! The function call to compute_andor_apply_rhs is as follows:
  !  compute_andor_apply_rhs(np1, nm1, n0, dt2, elem, hvcoord, hybrid, &
  !     deriv, nets, nete, compute_diagnostics, eta_ave_w, scale1, scale2, scale3)
  !
  !  This call returns the following:
  !
  !   u(np1) = scale3*u(nm1) + dt2*DSS[ nonstiffRHS(u(n0))*scale1 + stiffRHS(un0)*scale2 ]
  !
  !   nonstiffRHS and the stiffRHS are determined within the function and can be change by
  !   multiplying different terms by scale1 and scale2
  !
  !  Setting scale1=scale2=1.0, scale3=0.0, and dt2=1.0 returns the full rhs
  !
  !  DSS is the averaging procedure for the active and inactive nodes
  !

  call compute_andor_apply_rhs(ydot%tl_idx, ydot%tl_idx, y%tl_idx, &
        1.d0, y%elem, hvcoord, hybrid, deriv, y%nets, y%nete, &
        .false., bval*eta_ave_w, scale1, scale2, 0.d0)

  ! Zero out RHS vector elements to maintain (d/dt)phinh_i(nlevp) = 0
  do ie=ydot%nets,ydot%nete
    ydot%elem(ie)%state%phinh_i(:,:,nlevp,ydot%tl_idx) = 0.d0
  end do

  return
end function farkefun


function FColumnSolSolve(b, t, y, gamma, x) result(ierr)
  !-----------------------------------------------------------------
  ! Description: FColumnSolSolve is the routine called by ARKode to
  !     perform the linear solve at the state defined by (t,y),
  !     where b stores the right-hand side vector on input, and
  !     the solution vector on output.
  !
  ! Arguments:
  !       b - (NVec_t, input) NVec_t containing linear system RHS
  !       t - (dbl, input) current time
  !       y - (NVec_t) NVec_t containing current state
  !   gamma - (dbl, input) scaling factor for Jacobian in system
  !            matrix, A = I-gamma*J
  !       x - (NVec_t, output) NVec_t containing solution
  !-----------------------------------------------------------------

  !======= Inclusions ===========
  use arkode_mod,         only: get_hvcoord_ptr
  use control_mod,        only: theta_hydrostatic_mode
  use eos,                only: pnh_and_exner_from_eos
  use imex_mod,           only: get_dirk_jacobian
  use kinds,              only: real_kind
  use HommeNVector,       only: NVec_t
  use hybvcoord_mod,      only: hvcoord_t
  use physical_constants, only: g
  use dimensions_mod,     only: np, nlev, nlevp
  use parallel_mod,       only: abortmp
  use iso_c_binding

  !======= Declarations =========
  implicit none

  ! calling variables
  type(NVec_t),      intent(in) :: b
  real*8,            intent(in)    :: t
  type(NVec_t),      intent(in)    :: y
  real*8,            intent(in)    :: gamma
  type(NVec_t),      intent(out)   :: x
  integer                          :: ierr

  ! local variables
  type(hvcoord_t)       :: hvcoord
  real (kind=real_kind), pointer, dimension(:,:,:) :: phi_np1
  real (kind=real_kind), pointer, dimension(:,:,:) :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:) :: vtheta_dp
  real (kind=real_kind), pointer, dimension(:,:)   :: phis
  real (kind=real_kind) :: JacD(nlev,np,np)
  real (kind=real_kind) :: JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np)
  real (kind=real_kind) :: JacU2(nlev-2,np,np)
  real (kind=real_kind) :: pnh(np,np,nlev)
  real (kind=real_kind) :: pnh_i(np,np,nlevp)
  real (kind=real_kind) :: dphi(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: b1(nlev,np,np)
  integer :: i, j, k, ie, info(np,np), Ipiv(nlev,np,np)

  !======= Internals ============
  call get_hvcoord_ptr(hvcoord)

  ! set return value to success
  ierr = 0

  !--------------------------------
  ! perform solve
  !
  ! Notes: prim_advance_mod::compute_stage_value_dirk() performs a full Newton
  ! iteration on the implicit system
  !    w = w0 + g*dt2*(1-dpdpi(phi))
  !    phi = phi0 + g*dt2*w
  ! where w0 and phi0 contain the explicit portions of each time-dependent equation,
  ! g is the gravitational constant, dt2 is the scaled time step size (includes the
  ! diagonal coefficient for the DIRK method: dt*Ai(i,i)), and dpdpi = (dp/dpi).  Due to its
  ! structure, this is solved via substitution:
  !
  !    phi = phi0 + g*dt2*(w0 + g*dt2*(1-dpdpi(phi)))
  ! <=>
  !    phi = phi0 + g*dt2*w0 + (g*dt2)^2 - (g*dt2)^2*dpdpi(phi)
  ! <=>
  !    phi + (g*dt2)^2*dpdpi(phi) - phi0 - g*dt2*w0 - (g*dt2)^2 = 0
  !
  ! This nonlinear equation has tridiagonal Jacobian matrix, M = I + (g*dt2)^2*J,
  ! where J is the Jacobian matrix of dpdpi(phi).
  ! The diagonals of M are constructed in the call get_dirk_jacobian().
  !
  ! For our problem, we consider a 'state vector'
  !    U = [v1; v2; w; phi; vtheta_dp; dp3d]
  ! and we solve an IMEX problem of the form  U' = fE(U) + fI(U), where the implicit
  ! portion has structure
  !    fI(U) = [0; 0; fI_w(phi); fI_phi(w); 0; 0], where
  !    fI_W(phi) = g*(1 - dpdpi(phi))
  !    fI_phi(w) = g*w
  ! The resulting Jacobian for our nonlinear implicit solve has structure
  !    A = I-gamma*[0  0   0    0   0  0]
  !                [0  0   0    0   0  0]
  !                [0  0   0  -g*J  0  0]
  !                [0  0  g*I   0   0  0]
  !                [0  0   0    0   0  0]
  !                [0  0   0    0   0  0]
  !    gamma = dt*Ai(i,i) = dt2
  ! or equivalently, the linear system A*x = b corresponds to
  !      x_v1 = b_v1,  x_v2 = b_v2,  x_theta = b_theta,  x_dp3d = b_dp3d,
  ! and
  !       [       I     gamma*g*J ] [ x_w   ] = [ b_w   ]
  !       [ -gamma*g*I       I    ] [ x_phi ] = [ b_phi ]
  ! This 2x2 block linear system is equivalent to
  !       x_w + gamma*g*J*x_phi = b_w
  !       -gamma*g*x_w + x_phi  = b_phi
  !   <=>
  !       x_w = b_w - gamma*g*J*x_phi
  !       x_phi - gamma*g*(b_w - gamma*g*J*x_phi) = b_phi
  !   <=>
  !       x_w = b_w - gamma*g*J*x_phi
  !       (I + (gamma*g)^2*J) x_phi = b_phi + gamma*g*b_w
  !   <=>
  !       M*x_phi = (b_phi + gamma*g*b_w)
  !       x_w =  b_w + [x_phi - M*x_phi]/(gamma*g)
  ! So to solve our full linear system, we use the following steps:
  !    (a) construct M = I + (gamma*g)^2*J via a call to get_dirk_jacobian with dt2=gamma
  !    (b) construct right-hand side vector:  b1 = (b_phi + gamma*g*b_w)
  !    (c) solve M*x_phi = b1 for x_phi via calls to DGTTRF and DGTTRS
  !    (d) compute x_w = [x_phi - b_phi]/(gamma*g)
  !    (e) set x_v1 = b_v1, x_v2 = b_v2, x_theta = b_theta, x_dp3d = b_dp3d

  ! outer loop
  do ie=y%nets,y%nete

     !-------------
     ! step (a) construct M = I + (gamma*g)^2*J via a call to get_dirk_jacobian with dt2=gamma
     !-------------
     ! set pointers for this ie
     dp3d  => y%elem(ie)%state%dp3d(:,:,:,y%tl_idx)
     vtheta_dp  => y%elem(ie)%state%vtheta_dp(:,:,:,y%tl_idx)
     phi_np1 => y%elem(ie)%state%phinh_i(:,:,:,y%tl_idx)
     phis => y%elem(ie)%state%phis(:,:)

     ! compute intermediate variables for Jacobian calculation
     call pnh_and_exner_from_eos(hvcoord, vtheta_dp, dp3d, phi_np1, &
             pnh, exner, dpnh_dp_i, pnh_i_out=pnh_i)

     ! compute dp3d_i -- note that compute_stage_value_dirk
     ! computes this twice (with potential memory reference issues in the second pass)
     do k=1,nlev
       dphi(:,:,k)=phi_np1(:,:,k+1)-phi_np1(:,:,k)
     enddo


     ! get tridigonal matrix M
     info(:,:) = 0

     call get_dirk_jacobian(JacL,JacD,JacU,gamma,dp3d,dphi,pnh,1)

     ! loop over components of this ie
#if (defined COLUMN_OPENMP)
  !$omp parallel do private(i,j,k) collapse(2)
#endif
     do i=1,np
        do j=1,np
           !-------------
           ! step (b) construct right-hand side vector:  b1 = (b_phi + gamma*g*b_w)
           !-------------
           b1(:,i,j) = b%elem(ie)%state%phinh_i(i,j,1:nlev,b%tl_idx) &
                    + gamma*g*b%elem(ie)%state%w_i(i,j,1:nlev,b%tl_idx)

           !-------------
           ! step (c) solve M*x_phi = b1 for x_phi via calls to DGTTRF and DGTTRS
           !          note solution will be in b1 on output of DGTTRS
           !-------------
           call DGTTRF( nlev, JacL(:,i,j), JacD(:,i,j), JacU(:,i,j), JacU2(:,i,j), &
                        Ipiv(:,i,j), info(i,j) )
           call DGTTRS( 'N', nlev, 1, JacL(:,i,j), JacD(:,i,j), JacU(:,i,j), &
                        JacU2(:,i,j), Ipiv(:,i,j), b1(:,i,j), nlev, info(i,j) )
           x%elem(ie)%state%phinh_i(i,j,1:nlev,x%tl_idx) = b1(:,i,j)

           !-------------
           ! step (d) compute x_w = [x_phi - b_phi]/(gamma*g), then copy x_phi
           !          into b_phi (for output purposes)
           !-------------
           x%elem(ie)%state%w_i(i,j,1:nlev,x%tl_idx) = &
             (x%elem(ie)%state%phinh_i(i,j,1:nlev,x%tl_idx) - &
              b%elem(ie)%state%phinh_i(i,j,1:nlev,b%tl_idx)) / (gamma*g)

           !-------------
           ! step (e) set x_v1 = b_v1, x_v2 = b_v2, x_theta = b_theta, x_dp3d = b_dp3d
           !-------------
           x%elem(ie)%state%v(i,j,1,1:nlev,x%tl_idx) = &
             b%elem(ie)%state%v(i,j,1,1:nlev,b%tl_idx)
           x%elem(ie)%state%v(i,j,2,1:nlev,x%tl_idx) = &
             b%elem(ie)%state%v(i,j,2,1:nlev,b%tl_idx)
           x%elem(ie)%state%vtheta_dp(i,j,1:nlev,x%tl_idx) = &
             b%elem(ie)%state%vtheta_dp(i,j,1:nlev,b%tl_idx)
           x%elem(ie)%state%dp3d(i,j,1:nlev,x%tl_idx) = &
             b%elem(ie)%state%dp3d(i,j,1:nlev,b%tl_idx)

        end do  ! end j loop
     end do  ! end i loop
  end do   ! end ie loop

  return
end subroutine FColumnSolSolve

!=================================================================


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! NOTE: All of the remaining subroutines are not required when
! interfacing with ARKode, and so they are not currently
! implemented.  These may be provided to improve ARKode
! performance on this application; if so they should be 'enabled'
! by calling the relevant ARKode interface routine.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


subroutine farkjtimes(v_C, Jv_C, t, y_C, fy_C, h, ipar, rpar, v1_C, ierr)
  !-----------------------------------------------------------------
  ! Description: farkjtimes provides the Jacobian-vector product
  !    routine for the linearized Newton system.
  !
  !  Arguments:
  !       v_C - (ptr) C pointer to NVec_t vector to multiply
  !      Jv_C - (ptr) C pointer to NVec_t for result of Jacobian-vector product
  !         t - (dbl, input) current time
  !       y_C - (ptr) C pointer to NVec_t containing current solution
  !      fy_C - (ptr) C pointer to NVec_t containing current implicit ODE rhs
  !         h - (dbl, input) time step size for last internal step
  !      ipar - (long int(*), input) integer user parameter data
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here,
  !             unused)
  !      v1_C - (ptr) C pointer to NVec_t scratch vector
  !      ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  type(c_ptr),     intent(in), target :: v_C
  type(c_ptr),     intent(in), target :: Jv_C
  real*8,          intent(in)         :: t
  type(c_ptr),     intent(in), target :: y_C
  type(c_ptr),     intent(in), target :: fy_C
  real*8,          intent(in)         :: h
  integer(C_LONG), intent(in)         :: ipar(1)
  real*8,          intent(in)         :: rpar(1)
  type(c_ptr),     intent(in), target :: v1_C
  integer(C_INT),  intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: v  => NULL()
  type(NVec_t), pointer :: Jv => NULL()
  type(NVec_t), pointer :: y  => NULL()
  type(NVec_t), pointer :: fy => NULL()
  type(NVec_t), pointer :: v1 => NULL()


  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(v_C, v)
  call c_f_pointer(Jv_C, Jv)
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)
  call c_f_pointer(v1_C, v1)

  ! perform matrix-vector product, inserting result into Jv

  return
end subroutine farkjtimes
!=================================================================






subroutine farkjtsetup(t, y_C, fy_C, h, ipar, rpar, ierr)
  !-----------------------------------------------------------------
  ! Description: farkjtsetup performs any preparations for
  !    subsequent calls to farkjtimes.
  !
  !  Arguments:
  !         t - (dbl, input) current time
  !       y_C - (ptr) C pointer to NVec_t containing current solution
  !      fy_C - (ptr) C pointer to NVec_t containing current implicit ODE rhs
  !         h - (dbl, input) time step size for last internal step
  !      ipar - (long int(*), input) integer user parameter data
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here,
  !             unused)
  !      ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)         :: t
  type(c_ptr),     intent(in), target :: y_C
  type(c_ptr),     intent(in), target :: fy_C
  real*8,          intent(in)         :: h
  integer(C_LONG), intent(in)         :: ipar(1)
  real*8,          intent(in)         :: rpar(1)
  integer(C_INT),  intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: y  => NULL()
  type(NVec_t), pointer :: fy => NULL()


  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)

  ! perform setup for matrix-vector products

  return
end subroutine farkjtsetup
!=================================================================




subroutine farkpset(t, y_C, fy_C, jok, jcur, gamma, h, ipar, rpar, ierr)
  !-----------------------------------------------------------------
  ! Description: farkpset provides the preconditioner setup
  !    routine for the linearized Newton system.
  !
  !  Arguments:
  !         t - (dbl, input) current time
  !       y_C - (ptr) C pointer to NVec_t containing current solution
  !      fy_C - (ptr) C pointer to NVec_t containing current implicit ODE rhs
  !       jok - (int, input) flag denoting whether to recompute
  !              Jacobian-related data: 0=>recompute, 1=>unnecessary
  !      jcur - (int, output) output flag to say if Jacobian data
  !              was recomputed: 1=>was recomputed, 0=>was not
  !     gamma - (dbl, input) the scalar appearing in the Newton matrix
  !              A = M-gamma*J
  !         h - (dbl, input) time step size for last internal step
  !      ipar - (long int(*), input) integer user parameter data
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here,
  !             unused)
  !      ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)         :: t
  type(c_ptr),     intent(in), target :: y_C
  type(c_ptr),     intent(in), target :: fy_C
  integer(C_INT),  intent(in)         :: jok
  integer(C_INT),  intent(out)        :: jcur
  real*8,          intent(in)         :: gamma
  real*8,          intent(in)         :: h
  integer(C_LONG), intent(in)         :: ipar(1)
  real*8,          intent(in)         :: rpar(1)
  integer(C_INT),  intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: y  => NULL()
  type(NVec_t), pointer :: fy => NULL()

  !======= Internals ============

  ! initialize return value to success, jcur to not-recomputed
  ierr = 0
  jcur = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)

  ! return if no preconditioner update is required
  if (jok == 1)  return

  ! update the preconditioner

  ! set Jacobian recomputation flag
  jcur = 1

  return
end subroutine farkpset
!=================================================================




subroutine farkpsol(t, y_C, fy_C, r_C, z_C, gamma, delta, lr, ipar, &
                    rpar, ierr)
  !-----------------------------------------------------------------
  ! Description: farkpsol provides the preconditioner solve routine
  !    for the preconditioning of the linearized Newton system.
  !         i.e. solves P*z = r
  !
  ! Arguments:
  !         t - (dbl, input) current time
  !       y_C - (ptr) C pointer to NVec_t containing current solution
  !      fy_C - (ptr) C pointer to NVec_t containing current implicit ODE rhs
  !       r_C - (ptr) C pointer to NVec_t rhs vector of prec. system
  !       z_C - (ptr) C pointer to NVec_t solution vector of prec. system
  !     gamma - (dbl, input) scalar appearing in the Newton Matrix
  !             A = M-gamma*J
  !     delta - (dbl, input) desired tolerance if using an iterative
  !             method.  In that case, solve until
  !                  Sqrt[Sum((r-Pz).*ewt)^2] < delta
  !             where the ewt vector is obtainable by calling
  !             FARKGETERRWEIGHTS()
  !        lr - (int, input) flag indicating preconditioning type to
  !             apply: 1 => left,  2 => right
  !      ipar - (long int(*), input) integer user parameter data
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here,
  !             unused)
  !      ierr - (int, output) return flag: 0=>success,
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding
  use HommeNVector,   only: NVec_t
  use dimensions_mod, only: np, nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)         :: t
  type(c_ptr),     intent(in), target :: y_C
  type(c_ptr),     intent(in), target :: fy_C
  type(c_ptr),     intent(in), target :: r_C
  type(c_ptr),     intent(in), target :: z_C
  real*8,          intent(in)         :: gamma
  real*8,          intent(in)         :: delta
  integer(C_INT),  intent(in)         :: lr
  integer(C_LONG), intent(in)         :: ipar(1)
  real*8,          intent(in)         :: rpar(1)
  integer(C_INT),  intent(out)        :: ierr

  ! local variables
  type(NVec_t), pointer :: y  => NULL()
  type(NVec_t), pointer :: fy => NULL()
  type(NVec_t), pointer :: r  => NULL()
  type(NVec_t), pointer :: z  => NULL()

  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! dereference pointer for NVec_t objects
  call c_f_pointer(y_C, y)
  call c_f_pointer(fy_C, fy)
  call c_f_pointer(r_C, r)
  call c_f_pointer(z_C, z)

  ! perform preconditioner solve to fill z

  return
end subroutine farkpsol
!=================================================================
