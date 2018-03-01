module ODEMod
  !
  ! !DESCRIPTION:
  ! ode integrators for the biogeochemistry model
  ! Jinyun Tang, 2013
  !
  ! !USES:
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use betr_ctrl     , only : iulog => biulog
  use gbetrType     , only : gbetr_type
  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  real(r8),        parameter :: tiny = 1.e-23_r8
  logical, public, parameter :: ldebug_ode=.false.

  public :: ode_mbbks1, ode_adapt_mbbks1
  public :: ode_ebbks1, ode_adapt_ebbks1
  public :: ode_rk4
  public :: ode_rk2
  public :: get_tscal
  public :: ebbks
  public :: get_rerr


  interface get_rerr
     module procedure get_rerr_v
     module procedure get_rerr_s
  end interface get_rerr

contains

  !-------------------------------------------------------------------------------
  subroutine ode_ebbks1(odefun, extra, y0, nprimeq, neq, t, dt, y, pscal)
    ! !DESCRIPTION:
    !first order accurate explicit BBKS fixed time step positive preserving ode integrator
    !reference: Broekhuizen et al., 2008
    !!
    !
    implicit none
    ! !ARGUMENTS:
    class(gbetr_type),  intent(inout)  :: extra
    integer,            intent(in)  :: nprimeq !number of primary equations that are subject to positive constraint
    integer,            intent(in)  :: neq     !total number of equations
    real(r8),           intent(in)  :: y0(neq) !initial values
    real(r8),           intent(in)  :: t       !current time
    real(r8),           intent(in)  :: dt      !time step
    real(r8),           intent(out) :: y(neq)  !return values
    real(r8), optional, intent(out) :: pscal   !scaling factor
    external :: odefun

    ! !LOCAL VARIABLES:
    real(r8) :: f(neq)
    real(r8) :: pscal_loc
    call odefun(extra, y0, dt, t, nprimeq, neq, f)

    call ebbks(y0, f, nprimeq, neq, dt, y,pscal_loc)
    if(present(pscal))pscal=pscal_loc
  end subroutine ode_ebbks1
  !-------------------------------------------------------------------------------
  subroutine ode_ebbks2(odefun, extra, y0, nprimeq, neq, t, dt, y)
    ! !DESCRIPTION:
    !second order accurate explicit BBKS fixed time step positive preserving ode integrator
    !reference: Broekhuizen et al., 2008
    !!
    !
    implicit none
    ! !ARGUMENTS:
    class(gbetr_type),  intent(inout)  :: extra
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

    call odefun(extra, y0, dt, t, nprimeq, neq, f)
    call ebbks(y0, f, nprimeq, neq, dt, y1)
    ti=t+dt
    call odefun(extra, y1, dt, ti, nprimeq, neq, f1)
    do n = 1, neq
       f(n) = (f(n)+f1(n))*0.5_r8
    enddo
    call ebbks(y0, f, nprimeq, neq, dt, y)
  end subroutine ode_ebbks2
  !-------------------------------------------------------------------------------
  subroutine ode_mbbks1(odefun, extra, y0, nprimeq, neq, t, dt, y, pscal, bstatus)
    ! !DESCRPTION:
    !first order accurate implicit BBKS fixed time step positive preserving ode integrator
    use BetrStatusType    , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(gbetr_type),  intent(inout)  :: extra
    integer,            intent(in)  :: nprimeq
    integer,            intent(in)  :: neq
    real(r8),           intent(in)  :: y0(neq)
    real(r8),           intent(in)  :: t
    real(r8),           intent(in)  :: dt
    real(r8),           intent(out) :: y(neq)
    real(r8), optional, intent(out) :: pscal
    type(betr_status_type), intent(out) :: bstatus

    external :: odefun

    ! !LOCAL VARIABLES:
    real(r8) :: f(neq)
    real(r8) :: pscal1

    call bstatus%reset()
    call odefun(extra, y0, dt, t, nprimeq, neq, f)

    call mbbks(y0, f, nprimeq, neq, dt, y, pscal1, bstatus)

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
  subroutine ode_mbbks2(odefun, extra, y0, nprimeq, neq, t, dt, y, bstatus)
    !
    ! !DESCRIPTION:
    !second order implicit bkks ode integration with the adaptive time stepping
    use BetrStatusType    , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(gbetr_type),  intent(inout)  :: extra
    integer,  intent(in)  :: nprimeq
    integer,  intent(in)  :: neq
    real(r8), intent(in)  :: y0(neq)
    real(r8), intent(in)  :: t
    real(r8), intent(in)  :: dt
    real(r8), intent(out) :: y(neq)
    type(betr_status_type), intent(out) :: bstatus
    external :: odefun

    ! !LOCAL VARIABLES:
    real(r8) :: f(neq)
    real(r8) :: f1(neq)
    real(r8) :: y1(neq)
    real(r8) :: ti
    integer  :: n
    real(r8) :: nJ, pp
    real(r8) :: pscal
    call bstatus%reset()
    call odefun(extra, y0, dt, t, nprimeq, neq, f)

    call mbbks(y0, f, nprimeq, neq, dt, y1, pscal, bstatus)
    if(bstatus%check_status())return
    ti = t + dt
    call odefun(extra, y1, dt, ti,nprimeq, neq, f1)

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

    call mbbks(y0, f, nprimeq, neq, dt, y, pscal, bstatus)

  end subroutine ode_mbbks2
  !-------------------------------------------------------------------------------
  subroutine mbbks(y0, f, nprimeq, neq, dt, y, pscal, bstatus)
    ! !DESCRIPTION:
    ! mbbks update
    !
    ! !USES:
    use MathfuncMod       , only : safe_div
    use func_data_type_mod, only : func_data_type
    use BetrStatusType    , only : betr_status_type
    implicit none

    ! !ARGUMENTS:
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: f(neq)   ! derivative
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
    real(r8), intent(out) :: y(neq)   ! updated state variable
    real(r8), intent(out) :: pscal
    type(betr_status_type), intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    real(r8), pointer :: aj(:)
    real(r8) :: pmax
    real(r8) :: pm
    real(r8) :: a
    integer  :: n, nJ
    type(func_data_type) :: mbkks_data

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
          mbkks_data%nJ = nJ
          call GetGdtScalar(mbkks_data, pmax, pscal, bstatus)
          if(bstatus%check_status())return
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
  subroutine ode_adapt_mbbks1(odefun, extra, y0, nprimeq, neq, t, dt, y, bstatus)
    ! !DESCRIPTION:
    !first order implicit bkks ode integration with the adaptive time stepping
    !This could be used as an example for the implementation of time-adaptive
    !mbbks1.
    ! !NOTE:
    ! this code should only be used for mass positive ODE integration
    use BetrStatusType    , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(gbetr_type),  intent(inout)  :: extra
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: t        ! time stamp
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
    real(r8), intent(out) :: y(neq)   ! updated state variable
    type(betr_status_type), intent(out) :: bstatus
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

    call bstatus%reset()
    ! remove compiler warning about unused dummy arg
    if (t > 0.0_r8) continue

    dt2=dt
    dtmin=dt/64._r8
    dtr=dt
    tt=0._r8
    !make a copy of the solution at the current time step
    y=y0
    do
       if(dt2<=dtmin)then
          call odefun(extra, y, dt2, tt, nprimeq, neq, f)
          call mbbks(y, f, nprimeq, neq, dt2, yc, pscal, bstatus)
          if(bstatus%check_status())return
          dtr=dtr-dt2
          tt=tt+dt2
          y=yc
       else
          !get coarse grid solution
          call odefun(extra, y, dt2, tt, nprimeq, neq, f)
          call mbbks(y, f, nprimeq, neq, dt2, yc, pscal, bstatus)
          if(bstatus%check_status())return
          !get fine grid solution
          dt05=dt2*0.5_r8
          call mbbks(y,f,nprimeq, neq,dt05, yf, pscal, bstatus)
          if(bstatus%check_status())return
          tt2=tt+dt05
          ycp=yf
          call odefun(extra, ycp, dt05, tt, nprimeq, neq, f)
          call mbbks(ycp,f,nprimeq, neq,dt05,yf,pscal, bstatus)
          if(bstatus%check_status())return

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
    integer,  intent(in) :: neq     !number of equations
    real(r8), intent(in) :: yc(neq) !coarse solution
    real(r8), intent(in) :: yf(neq) !fine solution

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
  subroutine GetGdtScalar(mbkks_data, pmax, pp, bstatus)
    ! !DESCRIPTION:
    !get the gradient scaling factor for bkks integrator
    !
    ! !USES:
    use FindRootMod        , only : brent
    use func_data_type_mod , only : func_data_type
    use BetrstatusType     , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    type(func_data_type) , intent(in) :: mbkks_data
    real(r8)             , intent(in) :: pmax
    real(r8)             , intent(out):: pp
    type(betr_status_type), intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: iJ
    real(r8) :: f1, f2
    real(r8), parameter :: macheps = 1.e-8_r8
    real(r8), parameter :: tol = 1.e-8_r8

    call bstatus%reset()

    call gfunc_mbkks(0._r8, mbkks_data, f1)
    call gfunc_mbkks(pmax, mbkks_data, f2)
    call brent(pp, 0._r8, pmax, f1, f2, macheps, tol, mbkks_data, gfunc_mbkks, bstatus)

  end subroutine GetGdtScalar

  !-------------------------------------------------------------------------------
  subroutine gfunc_mbkks(p, mbkks_data, gf_value)
    ! !DESCRIPTION:
    !the bkks function
    !
    !USES
    use func_data_type_mod, only : func_data_type
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: p
    type(func_data_type), intent(in) :: mbkks_data
    real(r8), intent(out):: gf_value

    ! !LOCAL VARIABLES:
    integer :: jj
    real(r8), pointer :: aj(:)
    integer :: nJ
    real(r8) :: iJ

    aj => mbkks_data%aj
    nJ = mbkks_data%nJ
    iJ = mbkks_data%iJ
    gf_value = 1._r8
    do jj = 1, nJ
       gf_value = gf_value * (1._r8 + aj(jj) * p**(iJ))
    enddo
    gf_value = gf_value - p
    if(abs(gf_value) < 1.e-20_r8) then
       gf_value = 0._r8
    end if
  end subroutine gfunc_mbkks

  !-------------------------------------------------------------------------------
  subroutine ebbks(y0, f, nprimeq, neq, dt, y,ps)
    ! !DESCRIPTION:
    !ebbks update
    !
    !USES
    implicit none
    ! !ARGUMENTS:
    integer,            intent(in)  :: neq
    real(r8),           intent(in)  :: y0(neq)
    real(r8),           intent(in)  :: f(neq)
    real(r8),           intent(in)  :: dt
    integer,            intent(in)  :: nprimeq
    real(r8),           intent(out) :: y(neq)
    real(r8), optional, intent(out) :: ps
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
             write(*,'(A,5X,I3,3(5X,E20.10))')'debbkb',n,f(n),js,y0(n)
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
  subroutine ode_rk4(odefun,extra, y0, neq, t, dt, y )
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
    class(gbetr_type),  intent(inout)  :: extra
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

    call odefun(extra, y0, dt05, ti, neq, k1)

    y(:) = y0(:)
    call daxpy(neq, dt05, k1, 1, y, 1)

    ti = t + dt05
    call odefun(extra, y, dt05, ti, neq, k2)

    y(:) = y0(:)
    call daxpy(neq, dt05, k2, 1, y, 1)

    ti = t + dt05
    call odefun(extra, y, dt05, ti, neq, k3)

    y(:) = y0(:)
    call daxpy(neq, dt, k3, 1, y, 1)

    ti = t + dt
    call odefun(extra, y, dt, ti, neq, k4)

    do n = 1, neq
       kt(n) = k1(n)+2._r8*K2(n)+2._r8*k3(n)+k4(n)
    enddo
    a = dt / 6._r8

    y(:) = y0(:)
    call daxpy(neq, a, kt, 1, y, 1)

  end subroutine ode_rk4


  !-------------------------------------------------------------------------------
  subroutine ode_rk2(odefun, extra, y0, neq, t, dt, y )
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
    class(gbetr_type),  pointer  :: extra
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

    call odefun(extra, y0, dt, ti, neq, k1)

    y(:) = y0(:)
    call daxpy(neq, dt05, k1, 1, y, 1)

    ti = t + dt05
    call odefun(extra, y, dt05, ti, neq, k2)

    y(:) = y0(:)
    call daxpy(neq, dt, k2, 1, y, 1)
  end subroutine ode_rk2

  !-------------------------------------------------------------------------------
  subroutine ode_adapt_ebbks1(odefun, extra, y0, nprimeq, neq, t, dt, y)
    ! !DESCRIPTION:
    !first order implicit bkks ode integration with the adaptive time stepping
    !This could be used as an example for the implementation of time-adaptive
    !mbbks1.
    ! !NOTE:
    ! this code should only be used for mass positive ODE integration
    implicit none
    ! !ARGUMENTS:
    class(gbetr_type), target :: extra
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: t        ! time stamp
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
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
       print*,'ode',dt2,dtmin
       if(dt2<=dtmin)then
          call odefun(extra, y, dt2, tt, nprimeq, neq, f)
          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)
          dtr=dtr-dt2
          tt=tt+dt2
          y=yc
       else
          print*,'get coarse grid solution'
          call odefun(extra, y, dt2, tt, nprimeq, neq, f)
          print*,'ebbks'
          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)

          !get fine grid solution
          dt05=dt2*0.5_r8
          call ebbks(y,f,nprimeq, neq,dt05, yf, pscal)
          tt2=tt+dt05
          ycp=yf
          call odefun(extra,ycp, dt05, tt, nprimeq, neq, f)
          call ebbks(ycp,f,nprimeq, neq,dt05,yf,pscal)

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
  end subroutine ode_adapt_ebbks1
end module ODEMod
