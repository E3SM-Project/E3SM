module sphere_op_interface_mod

  use dimensions_mod,       only : np
  use kinds,                only : real_kind
  use derivative_mod_base,  only : derivative_t
  use element_mod,          only : element_t

  implicit none
  private

! ======================================
! Public Interfaces
! ======================================

  public  :: gradient_sphere_c_callable
  public  :: curl_sphere_wk_testcov_c_callable
  public  :: gradient_sphere_wk_testcov_c_callable
  public  :: divergence_sphere_wk_c_callable
  public  :: vorticity_sphere_c_callable
  public  :: divergence_sphere_c_callable
  public  :: laplace_sphere_wk_c_callable
  public  :: laplace_simple_c_callable
  public  :: vlaplace_sphere_wk_cartesian_c_callable
  public  :: vlaplace_sphere_wk_contra_c_callable

contains

  subroutine gradient_sphere_c_callable(s, dvv, dinv, grad) bind(c)
    use derivative_mod_base, only : gradient_sphere
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in)  :: s(np, np)
    real(kind=real_kind), intent(in)  :: dvv(np, np)
    real(kind=real_kind), intent(in)  :: dinv(np, np, 2, 2)
    real(kind=real_kind), intent(out) :: grad(np, np, 2)
    !
    ! Locals
    !
    type(derivative_t) :: deriv

    deriv%dvv = dvv

    grad = gradient_sphere(s, deriv, dinv)

  end subroutine gradient_sphere_c_callable


  subroutine curl_sphere_wk_testcov_c_callable(s,dvv,D,mp,ds) bind(c)
    use derivative_mod_base, only : curl_sphere_wk_testcov
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in) :: s(np,np)
    real(kind=real_kind), intent(in) :: dvv(np, np)
    real(kind=real_kind), intent(in) :: D(np, np, 2, 2)
    real(kind=real_kind), intent(in) :: mp(np, np)
    !
    ! Locals
    !
    type(derivative_t)   :: deriv
    type (element_t)     :: elem
    real(kind=real_kind) :: ds(np,np,2)

    deriv%dvv = dvv
    elem%D = D
    elem%mp = mp

    ds = curl_sphere_wk_testcov(s,deriv,elem)
          
  end subroutine curl_sphere_wk_testcov_c_callable

  subroutine gradient_sphere_wk_testcov_c_callable(s,dvv,metinv,metdet,D,mp,ds) bind(c)
    use derivative_mod_base, only : gradient_sphere_wk_testcov
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in)  :: s(np,np)
    real(kind=real_kind), intent(in)  :: dvv(np, np)
    real(kind=real_kind), intent(in)  :: D(np, np, 2, 2)
    real(kind=real_kind), intent(in)  :: metdet(np, np)
    real(kind=real_kind), intent(in)  :: metinv(np, np, 2, 2)
    real(kind=real_kind), intent(in)  :: mp(np, np)
    real(kind=real_kind), intent(out) :: ds(np,np,2)
    !
    ! Locals
    !
    type (derivative_t) :: deriv
    type (element_t) :: elem

    deriv%dvv = dvv
    elem%D = D
    elem%mp = mp
    elem%metinv = metinv
    elem%metdet = metdet

    ds = gradient_sphere_wk_testcov(s,deriv,elem)

  end subroutine gradient_sphere_wk_testcov_c_callable

  subroutine divergence_sphere_wk_c_callable(v, dvv, spheremp, Dinv, div) bind(c)
    use derivative_mod_base, only : divergence_sphere_wk
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in) :: v(np, np, 2)
    real(kind=real_kind), intent(in) :: dvv(np, np)
    real(kind=real_kind), intent(in) :: spheremp(np, np)
    real(kind=real_kind), intent(in) :: Dinv(np, np, 2, 2)
    real(kind=real_kind), intent(out) :: div(np, np)
    !
    ! Locals
    !
    type(derivative_t) :: deriv
    type(element_t) :: elem

    deriv%dvv = dvv
    elem%Dinv = Dinv
    elem%spheremp = spheremp

    div = divergence_sphere_wk(v, deriv, elem)

  end subroutine divergence_sphere_wk_c_callable

  subroutine vorticity_sphere_c_callable(v, dvv, rmetdet, d, vort) bind(c)
    use derivative_mod_base, only : vorticity_sphere
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in)  :: v(np, np, 2)
    real(kind=real_kind), intent(in)  :: dvv(np, np)
    real(kind=real_kind), intent(in)  :: rmetdet(np, np)
    real(kind=real_kind), intent(in)  :: d(np, np, 2, 2)
    real(kind=real_kind), intent(out) :: vort(np, np)
    !
    ! Locals
    !
    type(derivative_t) :: deriv
    type(element_t) :: elem

    deriv%dvv = dvv
    elem%D = d
    elem%rmetdet = rmetdet

    vort = vorticity_sphere(v, deriv, elem)

  end subroutine vorticity_sphere_c_callable

  subroutine divergence_sphere_c_callable(v, dvv, metdet, dinv, div) bind(c)
    use derivative_mod_base, only : divergence_sphere
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in) :: v(np, np, 2)
    real(kind=real_kind), intent(in) :: dvv(np, np)
    real(kind=real_kind), intent(in) :: metdet(np, np)
    real(kind=real_kind), intent(in) :: dinv(np, np, 2, 2)
    real(kind=real_kind), intent(out) :: div(np, np)
    !
    ! Locals
    !
    type(derivative_t) :: deriv
    type(element_t) :: elem

    deriv%dvv = dvv
    elem%Dinv = dinv
    elem%metdet = metdet
    elem%rmetdet = 1.0D0/metdet

    div = divergence_sphere(v, deriv, elem)

  end subroutine divergence_sphere_c_callable

! I have to check for this:
! For logical types, please note that the Fortran standard only guarantees
! interoperability between C99’s _Bool and Fortran’s C_Bool-kind logicals and
! C99 defines that true has the value 1 and false the value 0. Using any other
! integer value with GNU Fortran’s LOGICAL (with any kind parameter) gives an
! undefined result. (Passing other integer values than 0 and 1 to GCC’s _Bool is
! also undefined, unless the integer is explicitly or implicitly casted to
! _Bool.)

  subroutine laplace_sphere_wk_c_callable(s,dvv,dinv,spheremp,tensorVisc,&
             hvpower, hvscaling, var_coef,laplace) bind(c)
    use derivative_mod_base, only : laplace_sphere_wk
    !
    ! Inputs
    !
    use control_mod, only: hypervis_power, hypervis_scaling
    real(kind=real_kind), intent(in) :: s(np,np)
    real(kind=real_kind), intent(in) :: dvv(np, np)
    real(kind=real_kind), intent(in) :: dinv(np, np, 2, 2)
    real(kind=real_kind), intent(in) :: spheremp(np, np)
    real(kind=real_kind), intent(in) :: tensorVisc(np, np, 2, 2)
    logical, value, intent(in) :: var_coef
    real(kind=real_kind), intent(in) :: hvpower, hvscaling 
    real(kind=real_kind),intent(out)     :: laplace(np,np)
    !
    ! Locals
    !
    type (derivative_t) :: deriv
    type (element_t) :: elem

!redefining params from control_mod, not the usual homme practice, but...
    hypervis_power = hvpower
    hypervis_scaling = hvscaling
    deriv%dvv = dvv
    elem%Dinv = Dinv
    elem%spheremp = spheremp
    elem%tensorVisc = tensorVisc

    laplace=laplace_sphere_wk(s,deriv,elem,var_coef)

  end subroutine laplace_sphere_wk_c_callable


! not a homme function, for debugging cxx
  subroutine laplace_simple_c_callable(s,dvv,dinv,spheremp,laplace) bind(c)
    use derivative_mod_base, only : gradient_sphere, divergence_sphere_wk
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in)  :: s(np,np)
    real(kind=real_kind), intent(in)  :: dvv(np,np)
    real(kind=real_kind), intent(in)  :: dinv(np,np,2,2)
    real(kind=real_kind), intent(in)  :: spheremp(np, np)
    real(kind=real_kind), intent(out) :: laplace(np,np)
    !
    ! Locals
    !
    type (element_t) :: elem
    type (derivative_t) :: deriv
    real(kind=real_kind) :: grads(np,np,2)

    deriv%dvv = dvv
    elem%Dinv = dinv
    elem%spheremp = spheremp

    grads=gradient_sphere(s,deriv,dinv)
    laplace=divergence_sphere_wk(grads,deriv,elem)

  end subroutine laplace_simple_c_callable

!OG logics around hvpower, ... var_coef is not clear, but cleaning it
!would mean a different *nl for F and C, so, keeping these vars for now.
  subroutine vlaplace_sphere_wk_cartesian_c_callable(v, dvv, dinv, spheremp, &
             tensorVisc, vec_sph2cart, hvpower, hvscaling, var_coef, laplace) bind(c)
    use derivative_mod_base, only : vlaplace_sphere_wk_cartesian
    use control_mod,         only : hypervis_power, hypervis_scaling
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in) :: v(np,np,2)
    real(kind=real_kind), intent(in) :: dvv(np, np)
    real(kind=real_kind), intent(in) :: dinv(np, np, 2, 2)
    real(kind=real_kind), intent(in) :: spheremp(np, np)
    real(kind=real_kind), intent(in) :: tensorVisc(np, np, 2, 2)
    real(kind=real_kind), intent(in) :: vec_sph2cart(np, np, 3, 2)
    logical, value, intent(in) :: var_coef
    real(kind=real_kind), intent(in) :: hvpower, hvscaling
    real(kind=real_kind), intent(out)     :: laplace(np,np,2)
    !
    ! Locals
    !
    type (derivative_t) :: deriv
    type (element_t) :: elem

!redefining params from control_mod, not the usual homme practice, but...
    hypervis_power = hvpower
    hypervis_scaling = hvscaling
    deriv%dvv = dvv
    elem%Dinv = Dinv
    elem%spheremp = spheremp
    elem%tensorVisc = tensorVisc
    elem%vec_sphere2cart = vec_sph2cart

    laplace=vlaplace_sphere_wk_cartesian(v,deriv,elem,var_coef)

  end subroutine vlaplace_sphere_wk_cartesian_c_callable


!needs whatever is in strong div, vorticity, grad_testcov, curl_testcov  
!and nu_ratio
!and dvv
!Also, make this interface only for the case of const HV, so, no logic
!for var_coef hyperviscosity.
  subroutine vlaplace_sphere_wk_contra_c_callable(v, dvv, d, dinv, mp, spheremp, metinv,&
                                                  metdet, rmetdet, nu_ratio, laplace) bind(c)
    use derivative_mod_base, only : vlaplace_sphere_wk_contra
    !
    ! Inputs
    !
    real(kind=real_kind), intent(in) :: v(np,np,2)
    real(kind=real_kind), intent(in) :: nu_ratio
    real(kind=real_kind), intent(in) :: dvv(np, np)
    real(kind=real_kind), intent(in) :: D(np, np, 2, 2)
    real(kind=real_kind), intent(in) :: mp(np, np)
    real(kind=real_kind), intent(in) :: spheremp(np, np)
    real(kind=real_kind), intent(in) :: metdet(np, np)
    real(kind=real_kind), intent(in) :: metinv(np, np, 2, 2)
    real(kind=real_kind), intent(in) :: rmetdet(np, np)
    real(kind=real_kind), intent(in) :: dinv(np, np, 2, 2)
    real(kind=real_kind), intent(out) :: laplace(np,np,2)
    !
    ! Locals
    !
    type (derivative_t) :: deriv
    type (element_t) :: elem

    deriv%dvv = dvv
    elem%D = D
    elem%mp = mp
    elem%spheremp = spheremp
    elem%metinv = metinv
    elem%metdet = metdet
    elem%Dinv = dinv
    elem%rmetdet = rmetdet

    laplace = vlaplace_sphere_wk_contra(v,deriv,elem,.false.,nu_ratio)

  end subroutine vlaplace_sphere_wk_contra_c_callable

end module sphere_op_interface_mod
