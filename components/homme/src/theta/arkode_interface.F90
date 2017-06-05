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
!    farkjtimes -- implements a Jacobian-vector product, J*v, 
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
!=================================================================



subroutine arkode_init(t0, dt, y, rtol, atol, iout, rout, ierr)
  !-----------------------------------------------------------------
  ! Description: arkode_init initializes the ARKode solver.
  !   Arguments:
  !       t0 - (dbl, input) initial time
  !       dt - (dbl, input) time step size to use (first step)
  !        y - (vec, input) template solution vector
  !     rtol - (dbl, input) relative tolerance (for iterative solves)
  !     atol - (vec, input) absolute tolerance (for iterative solves),
  !            here I use a vector of the same shape/type as y, but
  !            this could instead be scalar-valued for all of y
  !     iout - (int*, input) integer solver parameter storage
  !     rout - (dbl*, input) real solver parameter storage
  !     ierr - (int, output) return flag: 0=>success, 
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use FortranVector
  use iso_c_binding

  !======= Declarations =========
  implicit none
  
  ! calling variables
  real*8,          intent(in)  :: t0
  real*8,          intent(in)  :: dt
  type(FVec),      intent(in)  :: y
  real*8,          intent(in)  :: rtol
!  type(FVec),      intent(in)  :: atol
  real*8,          intent(in)  :: atol
  integer(C_LONG), intent(in)  :: iout(40)
  real*8,          intent(in)  :: rout(40)
  integer(C_INT),  intent(out) :: ierr

  ! local variables
  integer(C_INT)  :: idef, iatol, imex, precLR, gstype, maxl
  real*8          :: rpar(1), lintol
  integer(C_LONG) :: lidef, ipar(1)
       

  !======= Internals ============

  ! initialize error flag
  ierr = 0

  ! initialize ARKode data & operators
  !    Nvector specs
  idef = 4    ! flag specifying which SUNDIALS solver will be used (4=ARKode)
  call fnvextinit(idef, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: fnvextinit failed'
  endif

  !    ARKode dataspace
  iatol = 1    ! specify type for atol: 1=scalar, 2=array
  imex = 1     ! specify problem type: 0=implicit, 1=explicit, 2=imex
  call farkmalloc(t0, y, imex, iatol, rtol, atol, &
                  iout, rout, ipar, rpar, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farkmalloc failed'
  endif

  !    ARKode options
  print *, 'Setting ARKode options'

  !      Indicate that we will set time step sizes ourselves, and the step 
  !      size to use on the first time step (disable adaptivity)
  call farksetrin('FIXED_STEP', dt, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farksetrin failed'
  endif

  !      Requested order of accuracy for ARK method (3,4,5 are supported).
  !      Alternately, a user can supply a custom Butcher table pair to 
  !      define their ARK method, by calling FARKSETARKTABLES()
  lidef = 4
  call farksetiin('ORDER', lidef, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farksetiin failed'
  endif

  !      To indicate that the implicit problem is linear, make the following 
  !      call.  The argument specifies whether the linearly implicit problem
  !      changes as the problem evolves (1) or not (0)
  lidef = 0
  call farksetiin('LINEAR', lidef, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farksetiin failed'
  endif

  !      Indicate use of the GMRES linear solver, the arguments indicate:
  !      precLR -- type of preconditioning: 0=none, 1=left, 2=right, 3=left+right
  !      gstype -- type of Gram-Schmidt orthogonalization: 1=modified, 2=classical
  !      maxl -- maximum size of Krylov subspace (# of iterations/vectors)
  !      lintol -- linear convergence tolerance factor (0 indicates default); this 
  !                example is very stiff so it requires tight linear solves
  precLR = 0
  gstype = 1
  maxl = 50
  lintol = 1.d-3
  call farkspgmr(precLR, gstype, maxl, lintol, ierr)
  if (ierr /= 0) then
     write(0,*) ' arkode_init: farkspgmr failed'
  endif

  !      Indicate to use our own Jacobian-vector product routine (otherwise it 
  !      uses a finite-difference approximation)
  !idef = 1
  !call farkspilssetjac(idef, ierr)
  !if (ierr /= 0) then
  !   write(0,*) ' arkode_init: farkspilssetjac failed'
  !endif

  !      Indicate to use our own preconditioner setup/solve routines (otherwise 
  !      preconditioning is disabled)
  !idef = 1
  !call farkspilssetprec(idef, ierr)
  !if (ierr /= 0) then
  !   write(0,*) ' arkode_init: farkspilssetprec failed'
  !endif

  ! output ARKode solver parameters to screen
  print *, '  '
  call farkwriteparameters(ierr)

  print *, 'Finished ARKode initialization'

  return
end subroutine arkode_init
!=================================================================



subroutine farkifun(t, y, fy, ipar, rpar, ierr)
  !-----------------------------------------------------------------
  ! Description: farkifun provides the implicit portion of the right
  !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
  !
  ! Arguments:
  !       t - (dbl, input) current time
  !       y - (vec, input) current solution
  !      fy - (vec, output) right-hand side function
  !    ipar - (long int(*), input) integer user parameter data
  !           (passed back here, unused)
  !    rpar - (dbl(*), input) real user parameter data (passed here, 
  !           unused)
  !    ierr - (int, output) return flag: 0=>success, 
  !            1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use FortranVector
  use iso_c_binding
  use prim_advance_mod, only: arkode_pars,compute_andor_apply_rhs
  use element_mod,      only: element_t
  use hybrid_mod,       only: hybrid_t
  use derivative_mod,   only: derivative_t
  use hybvcoord_mod,    only: hvcoord_t
  use dimensions_mod,   only: np,nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)            :: t
  type(FVec),      intent(in),    target :: y(:)
  type(FVec),      intent(inout), target :: fy(:)
  integer(C_LONG), intent(in)            :: ipar(1)
  type(arkode_pars), intent(in)          :: rpar
  integer(C_INT),  intent(out)           :: ierr
 
 
  ! local variables
  type (element_t), allocatable :: elem(:)
  integer :: ie,nets,nete,qn0
  !======= Internals ============
  real*8, allocatable :: Fvectemp(:,:,:,:,:)
  type (hvcoord_t)                   :: hvcoord
  type (hybrid_t)                    :: hybrid
  type (derivative_t)                :: deriv

  allocate(Fvectemp(nete-nets+1,np,np,nlev,6))

  ! set constants, extract solution components
  ! fill implicit portion of RHS
   nets=rpar%nets
   nete=rpar%nete
   allocate(elem(nete-nets+1))
   qn0=rpar%qn0

   Fvectemp(:,:,:,:,1)=reshape(y(:)%u,(/nete-nets+1,np,np,nlev/))
   Fvectemp(:,:,:,:,2)=reshape(y(:)%v,(/nete-nets+1,np,np,nlev/))
   Fvectemp(:,:,:,:,3)=reshape(y(:)%w,(/nete-nets+1,np,np,nlev/))
   Fvectemp(:,:,:,:,4)=reshape(y(:)%phi,(/nete-nets+1,np,np,nlev/))
   Fvectemp(:,:,:,:,5)=reshape(y(:)%theta_dp_cp,(/nete-nets+1,np,np,nlev/))
   Fvectemp(:,:,:,:,6)=reshape(y(:)%dp3d,(/nete-nets+1,np,np,nlev/))

   do ie=nets,nete
     elem(ie)%state%v(:,:,1,:,1)         = Fvectemp(nete-nets+ie,:,:,:,1)
     elem(ie)%state%v(:,:,2,:,1)         = Fvectemp(nete-nets+ie,:,:,:,2)
     elem(ie)%state%w(:,:,:,1)           = Fvectemp(nete-nets+ie,:,:,:,3)
     elem(ie)%state%phi(:,:,:,1)         = Fvectemp(nete-nets+ie,:,:,:,4)
     elem(ie)%state%theta_dp_cp(:,:,:,1) = Fvectemp(nete-nets+ie,:,:,:,5)
     elem(ie)%state%dp3d(:,:,:,1)        = Fvectemp(nete-nets+ie,:,:,:,6)
   end do 

   call compute_andor_apply_rhs(1,1,1,qn0,1.d0,elem,hvcoord,hybrid,&
       deriv,nets,nete,.false.,1.d0,0.d0,1.d0,0.d0)

   do ie=nets,nete
     Fvectemp(nete-nets+ie,:,:,:,1) = elem(ie)%state%v(:,:,1,:,1)
     Fvectemp(nete-nets+ie,:,:,:,2) = elem(ie)%state%v(:,:,2,:,1)
     Fvectemp(nete-nets+ie,:,:,:,3) = elem(ie)%state%w(:,:,:,1)
     Fvectemp(nete-nets+ie,:,:,:,4) = elem(ie)%state%phi(:,:,:,1)
     Fvectemp(nete-nets+ie,:,:,:,5) = elem(ie)%state%theta_dp_cp(:,:,:,1)
     Fvectemp(nete-nets+ie,:,:,:,6) = elem(ie)%state%dp3d(:,:,:,1)
   end do
 
   fy(:)%u           = reshape(Fvectemp(:,:,:,:,1),(/(nete-nets+1)*np*np*nlev/))
   fy(:)%v           = reshape(Fvectemp(:,:,:,:,2),(/(nete-nets+1)*np*np*nlev/))
   fy(:)%w           = reshape(Fvectemp(:,:,:,:,3),(/(nete-nets+1)*np*np*nlev/))
   fy(:)%phi         = reshape(Fvectemp(:,:,:,:,4),(/(nete-nets+1)*np*np*nlev/))
   fy(:)%theta_dp_cp = reshape(Fvectemp(:,:,:,:,5),(/(nete-nets+1)*np*np*nlev/))
   fy(:)%dp3d        = reshape(Fvectemp(:,:,:,:,6),(/(nete-nets+1)*np*np*nlev/))
    deallocate(Fvectemp,elem)
   ierr = 0

  return
end subroutine farkifun
!=================================================================




subroutine farkefun(t, y, fy, ipar, rpar, ierr)
  !-----------------------------------------------------------------
  ! Description: farkefun provides the explicit portion of the right
  !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
  !
  ! Arguments:
  !       t - (dbl, input) current time
  !       y - (vec, input) current solution
  !      fy - (vec, output) right-hand side function
  !    ipar - (long int(*), input) integer user parameter data
  !           (passed back here, unused)
  !    rpar - (dbl(*), input) real user parameter data (passed here, 
  !           unused)
  !    ierr - (int, output) return flag: 0=>success, 
  !            1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  use FortranVector
  use iso_c_binding
  use prim_advance_mod, only: arkode_pars,compute_andor_apply_rhs
  use element_mod,      only: element_t
  use hybrid_mod,       only: hybrid_t
  use derivative_mod,   only: derivative_t
  use hybvcoord_mod,    only: hvcoord_t
  use dimensions_mod,   only: np,nlev

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)            :: t
  type(FVec),      intent(in),    target :: y(:)
  type(FVec),      intent(inout), target :: fy(:)
  integer(C_LONG), intent(in)            :: ipar(1)
  type(arkode_pars), intent(in)          :: rpar
  integer(C_INT),  intent(out)           :: ierr
 
 

  ! local variables
  type (element_t), allocatable :: elem(:)
  integer :: ie,nets,nete,qn0
  !======= Internals ============
  real*8, allocatable         :: Fvectemp(:,:,:,:,:)
  type (hvcoord_t)    :: hvcoord
  type (hybrid_t)     :: hybrid
  type (derivative_t) :: deriv

  ! set constants, extract solution components
  ! fill implicit portion of RHS
   nets=rpar%nets
   nete=rpar%nete
   allocate(elem(nete-nets+1),Fvectemp(nete-nets+1,np,np,nlev,6))
   qn0=rpar%qn0

   do ie=nets,nete
     elem(ie)%state%v(:,:,1,:,1)         = Fvectemp(nete-nets+ie,:,:,:,1)
     elem(ie)%state%v(:,:,2,:,1)         = Fvectemp(nete-nets+ie,:,:,:,2)
     elem(ie)%state%w(:,:,:,1)           = Fvectemp(nete-nets+ie,:,:,:,3)
     elem(ie)%state%phi(:,:,:,1)         = Fvectemp(nete-nets+ie,:,:,:,4)
     elem(ie)%state%theta_dp_cp(:,:,:,1) = Fvectemp(nete-nets+ie,:,:,:,5)
     elem(ie)%state%dp3d(:,:,:,1)        = Fvectemp(nete-nets+ie,:,:,:,6)
   end do 

   call compute_andor_apply_rhs(1,1,1,qn0,1.d0,elem,hvcoord,hybrid,&
       deriv,nets,nete,.false.,1.d0,1.d0,0d0,0.d0)

   do ie=nets,nete
     Fvectemp(nete-nets+ie,:,:,:,1) = elem(ie)%state%v(:,:,1,:,1)
     Fvectemp(nete-nets+ie,:,:,:,2) = elem(ie)%state%v(:,:,2,:,1)
     Fvectemp(nete-nets+ie,:,:,:,3) = elem(ie)%state%w(:,:,:,1)
     Fvectemp(nete-nets+ie,:,:,:,4) = elem(ie)%state%phi(:,:,:,1)
     Fvectemp(nete-nets+ie,:,:,:,5) = elem(ie)%state%theta_dp_cp(:,:,:,1)
     Fvectemp(nete-nets+ie,:,:,:,6) = elem(ie)%state%dp3d(:,:,:,1)
   end do

  fy(:)%u           = reshape(Fvectemp(:,:,:,:,1),(/(nete-nets+1)*np*np*nlev/))
  fy(:)%v           = reshape(Fvectemp(:,:,:,:,2),(/(nete-nets+1)*np*np*nlev/))
  fy(:)%w           = reshape(Fvectemp(:,:,:,:,3),(/(nete-nets+1)*np*np*nlev/))
  fy(:)%phi         = reshape(Fvectemp(:,:,:,:,4),(/(nete-nets+1)*np*np*nlev/))
  fy(:)%theta_dp_cp = reshape(Fvectemp(:,:,:,:,5),(/(nete-nets+1)*np*np*nlev/))
  fy(:)%dp3d        = reshape(Fvectemp(:,:,:,:,6),(/(nete-nets+1)*np*np*nlev/))
  deallocate(elem,Fvectemp)
  ierr = 0

  return
end subroutine farkefun
!=================================================================


! ============= should this include hybrid,deriv, and hvcoord as inputs?

subroutine farkjtimes(v, Jv, t, y, fy, h, ipar, rpar, v1, ierr)
  !-----------------------------------------------------------------
  ! Description: farkjtimes provides the Jacobian-vector product
  !    routine for the linearized Newton system.
  ! 
  !  Arguments:
  !         v - (vec, input) vector to multiply
  !        Jv - (vec, output) result of Jacobian-vector product
  !         t - (dbl, input) current time
  !         y - (vec, input) current solution
  !        fy - (vec, input) current implicit ODE rhs of ODE, fi(t,y)
  !         h - (dbl, input) time step size for last internal step
  !      ipar - (long int(*), input) integer user parameter data 
  !             (passed back here, unused)
  !      rpar - (dbl(*), input) real user parameter data (passed here, 
  !             unused)
  !        v1 - (vec) scratch vector with same size as y
  !      ierr - (int, output) return flag: 0=>success, 
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use FortranVector
  use iso_c_binding

  !======= Declarations =========
  implicit none

  ! calling variables
  type(FVec),      intent(in)    :: v
  type(FVec),      intent(out)   :: Jv
  real*8,          intent(in)    :: t
  type(FVec),      intent(in)    :: y
  type(FVec),      intent(in)    :: fy
  real*8,          intent(in)    :: h
  integer(C_LONG), intent(in)    :: ipar(1)
  real*8,          intent(in)    :: rpar(1)
  type(FVec),      intent(inout) :: v1
  integer(C_INT),  intent(out)   :: ierr

  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! perform matrix-vector product, inserting result into Jv

  return
end subroutine farkjtimes
!=================================================================




subroutine farkpset(t, y, fy, jok, jcur, gamma, h, ipar, rpar, &
                    v1, v2, v3, ierr)
  !-----------------------------------------------------------------
  ! Description: farkpset provides the preconditioner setup 
  !    routine for the linearized Newton system.
  ! 
  !  Arguments:
  !         t - (dbl, input) current time
  !         y - (vec, input) current solution
  !        fy - (vec, input) current implicit ODE rhs of ODE, fi(t,y)
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
  !        v1 - (vec) scratch vector with same size as y
  !        v2 - (vec) scratch vector with same size as y
  !        v3 - (vec) scratch vector with same size as y
  !      ierr - (int, output) return flag: 0=>success, 
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use FortranVector
  use iso_c_binding

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)    :: t
  type(FVec),      intent(in)    :: y
  type(FVec),      intent(in)    :: fy
  integer(C_INT),  intent(in)    :: jok
  integer(C_INT),  intent(out)   :: jcur
  real*8,          intent(in)    :: gamma
  real*8,          intent(in)    :: h
  integer(C_LONG), intent(in)    :: ipar(1)
  real*8,          intent(in)    :: rpar(1)
  type(FVec),      intent(inout) :: v1
  type(FVec),      intent(inout) :: v2
  type(FVec),      intent(inout) :: v3
  integer(C_INT),  intent(out)   :: ierr
  
  !======= Internals ============

  ! initialize return value to success, jcur to not-recomputed
  ierr = 0
  jcur = 0

  ! return if no preconditioner update is required
  if (jok == 1)  return
  
  ! update the preconditioner 

  ! set Jacobian recomputation flag
  jcur = 1

  return
end subroutine farkpset
!=================================================================




subroutine farkpsol(t, y, fy, r, z, gamma, delta, lr, ipar, rpar, vt, ierr)
  !-----------------------------------------------------------------
  ! Description: farkpsol provides the preconditioner solve routine 
  !    for the preconditioning of the linearized Newton system.
  !         i.e. solves P*z = r
  !
  ! Arguments:
  !         t - (dbl, input) current time
  !         y - (vec, input) current solution
  !        fy - (vec, input) current implicit ODE rhs of ODE, fi(t,y)
  !         r - (vec, input) rhs vector of prec. system 
  !         z - (vec, output) solution vector of prec. system
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
  !        vt - (vec) scratch vector with same size as y
  !      ierr - (int, output) return flag: 0=>success, 
  !             1=>recoverable error, -1=>non-recoverable error
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use FortranVector
  use iso_c_binding

  !======= Declarations =========
  implicit none

  ! calling variables
  real*8,          intent(in)    :: t
  type(FVec),      intent(in)    :: y
  type(FVec),      intent(in)    :: fy
  type(FVec),      intent(in)    :: r
  type(FVec),      intent(out)   :: z
  real*8,          intent(in)    :: gamma
  real*8,          intent(in)    :: delta
  integer(C_INT),  intent(in)    :: lr
  integer(C_LONG), intent(in)    :: ipar(1)
  real*8,          intent(in)    :: rpar(1)
  type(FVec),      intent(inout) :: vt
  integer(C_INT),  intent(out)   :: ierr

  !======= Internals ============

  ! set return value to success
  ierr = 0

  ! perform preconditioner solve to fill z

  return
end subroutine farkpsol
!=================================================================




subroutine farkdiags(iout, rout)
  !-----------------------------------------------------------------
  ! Description: subroutine to output arkode diagnostics
  !
  ! Arguments:
  !        iout - (long int*, input) integer optional outputs 
  !        rout - (dbl*, input) real optional outputs
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use iso_c_binding

  !======= Declarations =========
  implicit none

  integer(C_LONG), intent(in) :: iout(40)
  real*8,          intent(in) :: rout(40)

  !======= Internals ============

  ! general solver statistics
  print *, '   '
  print *,  ' ARKode output values:'
  print '(4x,A,i9)','Total internal steps taken =',iout(3)
  print '(4x,A,i9)','   stability-limited steps =',iout(4)
  print '(4x,A,i9)','    accuracy-limited steps =',iout(5)
  print '(4x,A,i9)','  internal steps attempted =',iout(6)
  print '(4x,A,i9)','Total explicit rhs calls   =',iout(7)
  print '(4x,A,i9)','Total implicit rhs calls   =',iout(8)

  print '(4x,A,i9)','Total nonlinear iterations =',iout(11)
  if (iout(13) > 0)  print '(4x,A,i9)','Total root function calls  =',iout(13)

  ! linear solver statistics
  print '(4x,A,i9)','Total linear iterations    =',iout(21)
  if (iout(9) > 0)  print '(4x,A,i9)','Num lin solver setup calls =',iout(9)
  if (iout(17) > 0) print '(4x,A,i9)','Num lin solver rhs calls   =',iout(17)
  if (iout(18) > 0) print '(4x,A,i9)','Num lin solver Jac calls   =',iout(18)
  if (iout(19) > 0) print '(4x,A,i9)','Num PSet routine calls     =',iout(19)
  if (iout(20) > 0) print '(4x,A,i9)','Num PSolve routine calls   =',iout(20)

  ! error statistics
  if (iout(10) > 0) print '(4x,A,i9)','Num error test failures    =',iout(10)
  if (iout(12) > 0) print '(4x,A,i9)','Num nonlin conv failures   =',iout(12)
  if (iout(22) > 0) print '(4x,A,i9)','Num linear conv failures   =',iout(22)


  ! general time-stepping information
  print '(4x,A,es12.5)','First internal step size   =',rout(1)
  print '(4x,A,es12.5)','Last internal step size    =',rout(2)
  print '(4x,A,es12.5)','Next internal step size    =',rout(3)
  print '(4x,A,es12.5)','Current internal time      =',rout(4)
  print '(4x,A,es12.5)','Suggested tol scale factor =',rout(5)
  print *, '   '

  return
end subroutine farkdiags
!=================================================================
