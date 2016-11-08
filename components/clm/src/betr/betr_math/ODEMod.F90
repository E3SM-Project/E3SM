module ODEMod
  !
  ! !DESCRIPTION:
  ! ode integrators for the biogeochemistry model
  ! Jinyun Tang, 2013
  !
  ! !USES:
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use clm_varctl            , only : iulog
  implicit none
  save
  private

  public :: ode_mbbks1, ode_adapt_mbbks1
  public :: ode_ebbks1
  public :: ode_rk4
  public :: ode_rk2
  real(r8), parameter :: tiny = 1.e-23_r8

  type, private:: mbkks_type
     real(r8), pointer :: aj(:)
     real(r8) :: iJ
     integer :: nJ
  end type mbkks_type
  logical,public :: ldebug_ode=.false.
  type(mbkks_type), private :: mbkks_data
  interface get_rerr
     module procedure get_rerr_v, get_rerr_s
  end interface get_rerr

contains

  !-------------------------------------------------------------------------------
  subroutine ode_ebbks1(odefun, y0, nprimeq, neq, t, dt, y, pscal)
    ! !DESCRIPTION:
    !first order accurate explicit BBKS fixed time step positive preserving ode integrator
    !reference: Broekhuizen et al., 2008
    !!
    !
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)            :: nprimeq !number of primary equations that are subject to positive constraint
    integer,  intent(in)            :: neq     !total number of equations
    real(r8), intent(in)            :: y0(neq) !initial values
    real(r8), intent(in)            :: t       !current time
    real(r8), intent(in)            :: dt      !time step
    real(r8), intent(out)           :: y(neq)  !return values
    real(r8), optional, intent(out) :: pscal   !scaling factor
    external :: odefun

    ! !LOCAL VARIABLES:
    real(r8) :: f(neq)
    real(r8) :: pscal_loc
    call odefun(y0, dt, t, nprimeq, neq, f)

    call ebbks(y0, f, nprimeq, neq, dt, y,pscal_loc)
    if(present(pscal))pscal=pscal_loc
  end subroutine ode_ebbks1
  !-------------------------------------------------------------------------------
  subroutine ode_ebbks2(odefun, y0, nprimeq, neq, t, dt, y)
    ! !DESCRIPTION:
    !second order accurate explicit BBKS fixed time step positive preserving ode integrator
    !reference: Broekhuizen et al., 2008
    !!
    !
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)  :: nprimeq    !number of primary equations that are subject to positive constraint
    integer,  intent(in)  :: neq        !total number of equations
    real(r8), intent(in)  :: y0(neq)    !initial values
    real(r8), intent(in)  :: t          !current time
    real(r8), intent(in)  :: dt         !time step
    real(r8), intent(out) :: y(neq)     !return value

    external :: odefun

    ! !LOCAL VARIABLES:
    real(r8) :: f(neq)
    real(r8) :: f1(neq)
    real(r8) :: y1(neq)
    real(r8) :: ti
    integer  :: n

    call odefun(y0, dt, t, nprimeq, neq, f)
    call ebbks(y0, f, nprimeq, neq, dt, y1)
    ti=t+dt
    call odefun(y1, dt, ti, nprimeq, neq, f1)
    do n = 1, neq
       f(n) = (f(n)+f1(n))*0.5_r8
    enddo
    call ebbks(y0, f, nprimeq, neq, dt, y)
  end subroutine ode_ebbks2
  !-------------------------------------------------------------------------------
  subroutine ode_mbbks1(odefun, y0, nprimeq, neq, t, dt, y, pscal)
    ! !DESCRPTION:
    !first order accurate implicit BBKS fixed time step positive preserving ode integrator
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)            :: nprimeq
    integer,  intent(in)            :: neq
    real(r8), intent(in)            :: y0(neq)
    real(r8), intent(in)            :: t
    real(r8), intent(in)            :: dt
    real(r8), intent(out)           :: y(neq)
    real(r8), optional, intent(out) :: pscal
    external :: odefun

    ! !LOCAL VARIABLES:
    real(r8) :: f(neq)
    real(r8) :: pscal1

    call odefun(y0, dt, t, nprimeq, neq, f)

    call mbbks(y0, f, nprimeq, neq, dt, y, pscal1)

    if(present(pscal))pscal=pscal1
  end subroutine ode_mbbks1
  !-------------------------------------------------------------------------------
  subroutine get_tscal(rerr,dt_scal,acc)
    !
    ! !DESCRIPTION:
    !obtain the time step scalar for adaptive ode
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: rerr        !input relative err
    real(r8), intent(out):: dt_scal     !output dt_scal, 2, 1, 0.5
    logical,  intent(out):: acc         !true or false

    ! !LOCAL VARIABLES:
    real(r8), parameter  :: rerr_thr=1.e-4_r8   !relative error threshold

    if(rerr<0.5*rerr_thr)then
       dt_scal = 2._r8
       acc = .true.
    elseif(rerr<rerr_thr)then
       dt_scal = 1._r8
       acc = .true.
    elseif(rerr<2._r8*rerr_thr)then
       dt_scal = 0.5_r8
       acc = .true.
    else
       dt_scal = 0.5_r8
       acc = .false.
    endif

  end subroutine get_tscal
  !-------------------------------------------------------------------------------
  subroutine ode_mbbks2(odefun, y0, nprimeq, neq, t, dt, y)
    !
    ! !DESCRIPTION:
    !second order implicit bkks ode integration with the adaptive time stepping
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)  :: nprimeq
    integer,  intent(in)  :: neq
    real(r8), intent(in)  :: y0(neq)
    real(r8), intent(in)  :: t
    real(r8), intent(in)  :: dt
    real(r8), intent(out) :: y(neq)
    external :: odefun


    ! !LOCAL VARIABLES:
    real(r8) :: f(neq)
    real(r8) :: f1(neq)
    real(r8) :: y1(neq)
    real(r8) :: ti
    integer  :: n
    real(r8) :: nJ, pp
    real(r8) :: pscal

    call odefun(y0, dt, t, nprimeq, neq, f)

    call mbbks(y0, f, nprimeq, neq, dt, y1, pscal)
    ti = t + dt
    call odefun(y1, dt, ti,nprimeq, neq, f1)

    pp = 1._r8
    nJ = 0._r8

    do n = 1, neq
       f(n) = (f(n) + f1(n))*0.5_r8
       if(f(n)<0._r8 .and. n<=nprimeq)then
          pp = pp * y0(n) / y1(n)
          nJ = nJ + 1._r8
       endif
    enddo
    if(nJ>0._r8)then
       do n = 1, neq
          f(n) = f(n) * pp**(1._r8/nJ)
       enddo
    endif

    call mbbks(y0, f, nprimeq, neq, dt, y, pscal)

  end subroutine ode_mbbks2
  !-------------------------------------------------------------------------------
  subroutine mbbks(y0, f, nprimeq, neq, dt, y, pscal)
    ! !DESCRIPTION:
    ! mbbks update
    !
    ! !USES:
    use MathfuncMod   , only : safe_div
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: f(neq)   ! derivative
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(out) :: y(neq)   ! updated state variable
    real(r8), intent(out) :: pscal

    ! !LOCAL VARIABLES:
    real(r8), pointer :: aj(:)
    real(r8) :: pmax
    real(r8) :: pm
    real(r8) :: a
    integer  :: n, nJ

    allocate(mbkks_data%aj(neq))
    aj => mbkks_data%aj
    nJ = 0
    pmax = 0._r8
    do n = 1, nprimeq
       if(f(n)<0._r8)then
          nJ = nJ  + 1
          pm = -y0(n)/(f(n)*dt)

          pm = min(pm,1.e30_r8)
          aj(nJ) = -safe_div(1._r8,pm)
          if(nJ==1)then
             pmax= pm
          else
             pmax = min(pm, pmax)
          endif
       endif
    enddo
    if(nJ>0)then
       pmax=min(1._r8,pmax**(nJ))


       !solve the gradient modifier function
       mbkks_data%nJ=nJ
       mbkks_data%iJ=1._r8/nJ
       if(pmax<1.e-8_r8)then
          pscal=pmax
       else
          pscal=GetGdtScalar(aj,nJ,pmax)
          pscal=pscal**(1._r8/nJ)
       endif
       !reduce the chance of negative y(n) from roundoff error
       pscal=pscal*0.9999_r8
    else
       pscal=1._r8
    endif

    y(:)=y0(:)
    a=pscal*dt
    !daxpy(N,DA,DX,INCX,DY,INCY)
    call daxpy(neq, a, f, 1, y, 1)
    deallocate(mbkks_data%aj)


  end subroutine mbbks

  !-------------------------------------------------------------------------------
  subroutine ode_adapt_mbbks1(odefun, y0, nprimeq, neq, t, dt, y)
    ! !DESCRIPTION:
    !first order implicit bkks ode integration with the adaptive time stepping
    !This could be used as an example for the implementation of time-adaptive
    !mbbks1.
    ! !NOTE:
    ! this code should only be used for mass positive ODE integration
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: t        ! time stamp
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(out) :: y(neq)   ! updated state variable
    external :: odefun

    ! !LOCAL VARIABLES:
    real(r8) :: yc(neq)    !coarse time stepping solution
    real(r8) :: yf(neq)    !fine time stepping solution
    real(r8) :: ycp(neq)   !temporary variable
    real(r8) :: f(neq)   ! derivative
    real(r8) :: dt2
    real(r8) :: dtr
    real(r8) :: dt05
    real(r8) :: dtmin
    real(r8) :: tt,tt2     !temporary variables
    logical  :: acc
    real(r8) :: rerr, dt_scal, pscal
    integer  :: n, nJ

    dt2=dt
    dtmin=dt/64._r8
    dtr=dt
    tt=0._r8
    !make a copy of the solution at the current time step
    y=y0
    do
       if(dt2<=dtmin)then
          call odefun(y, dt2, tt, nprimeq, neq, f)
          call mbbks(y, f, nprimeq, neq, dt2, yc, pscal)
          dtr=dtr-dt2
          tt=tt+dt2
          y=yc
       else
          !get coarse grid solution
          call odefun(y, dt2, tt, nprimeq, neq, f)
          call mbbks(y, f, nprimeq, neq, dt2, yc, pscal)

          !get fine grid solution
          dt05=dt2*0.5_r8
          call mbbks(y,f,nprimeq, neq,dt05, yf, pscal)
          tt2=tt+dt05
          ycp=yf
          call odefun(ycp, dt05, tt, nprimeq, neq, f)
          call mbbks(ycp,f,nprimeq, neq,dt05,yf,pscal)

          !determine the relative error
          rerr=get_rerr_v(yc,yf, neq)

          !determine time scalar factor
          call get_tscal(rerr,dt_scal,acc)

          if(acc)then
             dtr=dtr-dt2
             tt=tt+dt2
             y=yf
          endif
          dt2=dt2*dt_scal
          dt2=min(dt2,dtr)
       endif
       if(abs(dtr/dt)<1.e-4_r8)exit
    enddo

  end subroutine ode_adapt_mbbks1

  !-------------------------------------------------------------------------------
  function get_rerr_v(yc,yf,neq)result(rerr)
    !
    ! !DESCRIPTION:
    ! obtain the relative error
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: yc(neq) !coarse solution
    real(r8), intent(in) :: yf(neq) !fine solution
    integer, intent(in)  :: neq     !number of equations

    ! !LOCAL VARIABLES:
    real(r8) :: rerr
    real(r8) :: rtmp
    integer :: n
    rerr=abs(yc(1)-yf(1))/(abs(yf(1))+1.e-20_r8)
    do n = 2, neq
       rtmp=abs(yc(n)-yf(n))/(abs(yf(n))+1.e-20_r8)
       rerr=max(rerr,rtmp)
    enddo

  end function get_rerr_v
  !-------------------------------------------------------------------------------
  function get_rerr_s(yc,yf)result(rerr)
    !
    ! DESCRIPTION:
    ! obtain the relative error
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: yc   !coarse solution
    real(r8), intent(in) :: yf   !fine solution
    ! !LOCAL VARIABLES
    real(r8) :: rerr
    real(r8) :: rtmp
    integer :: n

    rerr=abs(yc-yf)/(abs(yf)+1.e-20_r8)

  end function get_rerr_s

  !-------------------------------------------------------------------------------
  function GetGdtScalar(aj,nJ,pmax)result(pp)
    ! !DESCRIPTION:
    !get the gradient scaling factor for bkks integrator
    !
    ! !USES:
    use FindRootMod, only : brent
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: aj(nJ)
    real(r8), intent(in) :: pmax
    integer,  intent(in) :: nJ
    ! !LOCAL VARIABLES:
    real(r8) :: iJ
    real(r8) :: f1, f2
    real(r8), parameter :: macheps = 1.e-8_r8
    real(r8), parameter :: tol = 1.e-8_r8

    real(r8) :: pp

    call gfunc_mbkks(0._r8, f1)
    call gfunc_mbkks(pmax, f2)
    call brent(pp, 0._r8, pmax, f1, f2, macheps, tol, gfunc_mbkks)
  end function GetGdtScalar
  !-------------------------------------------------------------------------------

  subroutine gfunc_mbkks(p, value)
    ! !DESCRIPTION:
    !the bkks function
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: p
    real(r8), intent(out):: value

    ! !LOCAL VARIABLES:
    integer :: jj
    real(r8), pointer :: aj(:)
    integer :: nJ
    real(r8) :: iJ

    aj => mbkks_data%aj
    nJ = mbkks_data%nJ
    iJ = mbkks_data%iJ
    value = 1._r8
    do jj = 1, nJ
       value = value * (1._r8 + aj(jj) * p**(iJ))
    enddo
    value = value - p
    if(abs(value)<1.e-20_r8)value=0._r8
  end subroutine gfunc_mbkks

  !-------------------------------------------------------------------------------
  subroutine ebbks(y0, f, nprimeq, neq, dt, y,ps)
    ! !DESCRIPTION:
    !ebbks update
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: y0(neq)
    real(r8), intent(in) :: f(neq)
    real(r8), intent(in) :: dt
    integer,  intent(in) :: nprimeq
    integer,  intent(in) :: neq
    real(r8), intent(out):: y(neq)
    real(r8), optional, intent(out):: ps
    ! !LOCAL VARIABLES:
    real(r8), parameter :: beta=0.999_r8  !scaling parameter
    real(r8) :: js, jsmin
    real(r8) :: p
    integer :: n, nJ

    nJ=0
    do n = 1, nprimeq
       if(f(n)<0._r8)then
          js = y0(n)/(-f(n)*dt)
          nJ=nJ+1
          if(nJ==1)then
             jsmin=js
          else
             jsmin=min(jsmin,js)
          endif
          if(ldebug_ode)then
             write(*,'(A,X,I3,3(X,E20.10))')'debbkb',n,f(n),js,y0(n)
          endif
       endif
    enddo
    p=1._r8
    if(nJ>0)then
       p = min(jsmin*beta,1._r8)
    endif

    y(:) = y0(:)
    if(present(ps))ps=p
    p = p * dt
    call daxpy(neq, p, f, 1, y, 1)

  end subroutine ebbks



  !-------------------------------------------------------------------------------
  subroutine ode_rk4(odefun, y0, neq, t, dt, y )
    !
    ! !DESCRIPTION:
    !  4-th order runge-kutta method for ode integration
    !  Solve differential equations with a non-adaptive method of order 4.
    !   call rk4(y, ODEFUN,t, dt,Y0, neq) integrates
    !   the system of differential equations y' = f(t,y) by stepping from T to
    !   t+dt. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
    !   The vector Y0 is the initial conditions at T0. Each row in the solution
    !   array Y corresponds to a time specified in TSPAN.
    !
    !   This is a non-adaptive solver. The step sequence is determined by TSPAN
    !   but the derivative function ODEFUN is evaluated multiple times per step.
    !
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)  :: neq
    real(r8), intent(in)  :: y0(neq)
    real(r8), intent(in)  :: t
    real(r8), intent(in)  :: dt
    real(r8), intent(out) :: y(neq)
    ! !LOCAL VARIABLES:
    real(r8) :: k1(neq)
    real(r8) :: k2(neq)
    real(r8) :: k3(neq)
    real(r8) :: k4(neq)
    real(r8) :: kt(neq)
    real(r8) :: ti, dt05, a
    integer :: n
    external :: odefun

    ti = t
    dt05 = dt * 0.5_r8

    call odefun(y0, dt05, ti, neq, k1)

    y(:) = y0(:)
    call daxpy(neq, dt05, k1, 1, y, 1)

    ti = t + dt05
    call odefun(y, dt05, ti, neq, k2)

    y(:) = y0(:)
    call daxpy(neq, dt05, k2, 1, y, 1)

    ti = t + dt05
    call odefun( y, dt05, ti, neq, k3)

    y(:) = y0(:)
    call daxpy(neq, dt, k3, 1, y, 1)

    ti = t + dt
    call odefun(y, dt, ti, neq, k4)

    do n = 1, neq
       kt(n) = k1(n)+2._r8*K2(n)+2._r8*k3(n)+k4(n)
    enddo
    a = dt / 6._r8

    y(:) = y0(:)
    call daxpy(neq, a, kt, 1, y, 1)

  end subroutine ode_rk4


  !-------------------------------------------------------------------------------
  subroutine ode_rk2(odefun, y0, neq, t, dt, y )
    !
    ! !DESCRIPTION:
    !  2-th order runge-kutta method for ode integration
    !   Solve differential equations with a non-adaptive method of order 2.
    !   call rk2(y, ODEFUN,t, dt,Y0, neq) integrates
    !   the system of differential equations y' = f(t,y) by stepping from T to
    !   t+dt. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
    !   The vector Y0 is the initial conditions at T0. Each row in the solution
    !   array Y corresponds to a time specified in TSPAN.

    !
    !   This is a non-adaptive solver. The step sequence is determined by TSPAN
    !   but the derivative function ODEFUN is evaluated multiple times per step.
    !
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)  :: neq
    real(r8), intent(in)  :: y0(neq)
    real(r8), intent(in)  :: t
    real(r8), intent(in)  :: dt
    real(r8), intent(out) :: y(neq)
    ! !LOCAL VARIABLES:
    real(r8) :: k1(neq)
    real(r8) :: k2(neq)
    real(r8) :: ti, dt05
    integer :: n
    external :: odefun

    ti = t
    dt05 = dt * 0.5_r8

    call odefun(y0, dt, ti, neq, k1)

    y(:) = y0(:)
    call daxpy(neq, dt05, k1, 1, y, 1)

    ti = t + dt05
    call odefun(y, dt05, ti, neq, k2)

    y(:) = y0(:)
    call daxpy(neq, dt, k2, 1, y, 1)
  end subroutine ode_rk2

end module ODEMod
