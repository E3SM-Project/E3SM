#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module imex_mod
  use kinds,              only: iulog, real_kind
  use dimensions_mod,     only: nlev, nlevp, np
  use hybvcoord_mod,      only: hvcoord_t
  use hybrid_mod,         only: hybrid_t
  use element_mod,        only: element_t
  use derivative_mod,     only: derivative_t
  use time_mod,           only: timelevel_t, timelevel_qdp
  use physical_constants, only: g, kappa
  use eos,                only: pnh_and_exner_from_eos, pnh_and_exner_from_eos2, phi_from_eos
  use element_state,      only: max_itercnt, max_deltaerr, max_reserr
  use control_mod,        only: theta_hydrostatic_mode, qsplit
  use perf_mod,           only: t_startf, t_stopf
#ifdef HOMMEXX_BFB_TESTING
  use iso_c_binding,      only: c_loc
#endif

  implicit none

contains

  subroutine compute_gwphis(gwh_i,dp3d,v,gradphis,hvcoord)
    !
    !  compute a vertical velocity induced by surface topography
    !  wh_i =  ubar grad phis 
    !
    real (kind=real_kind) :: gwh_i(np,np,nlevp)
    real (kind=real_kind) :: dp3d(np,np,nlev)
    real (kind=real_kind) :: v(np,np,2,nlev)
    real (kind=real_kind) :: gradphis(np,np,2)
    type (hvcoord_t)     , intent(in) :: hvcoord  

    ! local
    integer :: k
    real (kind=real_kind) :: v_i(np,np,2,nlevp)

    ! add wh_i(n0) term to RHS
    do k=2,nlev
       v_i(:,:,1,k) = (dp3d(:,:,k)*v(:,:,1,k) + &
            dp3d(:,:,k-1)*v(:,:,1,k-1) ) / (dp3d(:,:,k)+dp3d(:,:,k-1))
       v_i(:,:,2,k) = (dp3d(:,:,k)*v(:,:,2,k) + &
            dp3d(:,:,k-1)*v(:,:,2,k-1) ) / (dp3d(:,:,k)+dp3d(:,:,k-1))
    end do
    gwh_i(:,:,1)=0
    gwh_i(:,:,nlevp)=0  ! should be elem(ie)%state%w_i(:,:,nlevp,n0), but not used
    do k=2,nlev
       gwh_i(:,:,k) = (v_i(:,:,1,k)*gradphis(:,:,1) + v_i(:,:,2,k)*gradphis(:,:,2))&
            *hvcoord%hybi(k)
    enddo
  end subroutine compute_gwphis





  subroutine compute_stage_value_dirk(nm1,alphadt_nm1,n0,alphadt_n0,np1,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,itercount,itererr,verbosity_in)
    !===================================================================================
    ! this subroutine solves a stage value equation for a DIRK method which takes the form
    !
    ! gi = un0 + dt* sum(1:i-1)(aij n(gj)+a2ij s(gj)) + dt *a2ii s(gi) := y + dt a2ii s(gi)
    !
    ! It is assumed that un0 has the value of y and the computed value of gi is stored at
    ! unp1
    ! 
    ! w_0 = w(np1)
    ! phi_0 = phi(np1)
    ! Then solve (ovewriting w(np1),phi(np1): 
    ! w(np1) = w_0 + alphadt_n0 * SW(n0)  + alphadt_nm1*SW(nm1) +  dt2*SW(np1)
    ! phi(np1) = phi_0 + alphadt_n0 * SPHI(n0)  + alphadt_nm1*SPHI(nm1) +  dt2*SPHI(np1)
    !
    ! SW(nt) = g*dpnh_dp_i(nt)-1
    ! SPHI(nt) = g*w(nt) -g*a(k) u(nt) dot grad_phis
    !
    ! pnh_and_exner_from_eos()   used to compute dpnh_dp_i(nt) in SW(nt)
    ! compute_gwphis()            used to compute g*a(k) u(nt) dot grad_phis  in SPHI(nt)   
    !
    ! We then precompute:
    !  w_rhs = w_0 + alphadt_n0 * SW(n0)  + alphadt_nm1*SW(nm1) 
    !  phi_rhs = w_0 + alphad_n0 * SPHI(n0)  + alphadt1_nm1*SPHI(nm1) -  dt2*compute_gwphi(np1)
    ! and solve, via Newton iteration:
    !   w(np1) = w_rhs + dt2*SW(np1)
    !   phi(np1) = phi_rhs + dt2*g*w(np1)
    !
    !===================================================================================

    integer, intent(in) :: nm1,n0,np1,nets,nete
    real (kind=real_kind), intent(in) :: dt2
    integer :: itercount
    real (kind=real_kind) :: itererr
    real (kind=real_kind), intent(in) :: alphadt_n0,alphadt_nm1
    integer, intent(in), optional :: verbosity_in

    type (hvcoord_t)     , intent(in) :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (derivative_t)  , intent(in) :: deriv


    ! local
    real (kind=real_kind), pointer, dimension(:,:,:)   :: phi_np1
    real (kind=real_kind), pointer, dimension(:,:,:)   :: w_np1
    real (kind=real_kind) :: JacD(np,np,nlev)  , JacL(np,np,nlev-1)
    real (kind=real_kind) :: JacU(np,np,nlev-1), JacU2(nlev-2,np,np)
    real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
    real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
    real (kind=real_kind) :: w_n0(np,np,nlevp)    
    real (kind=real_kind) :: dphi(np,np,nlev)    
    real (kind=real_kind) :: dphi_n0(np,np,nlev)    
    real (kind=real_kind) :: phi_n0(np,np,nlevp)    
    real (kind=real_kind) :: Fn(np,np,nlev),x(np,np,nlev)
    real (kind=real_kind) :: gwh_i(np,np,nlevp)  ! w hydrostatic
    real (kind=real_kind) :: v_i(np,np,2,nlevp)  ! w hydrostatic

    real (kind=real_kind) :: Jac2D(np,np,nlev)  , Jac2L(np,np,nlev-1)
    real (kind=real_kind) :: Jac2U(np,np,nlev-1)

    real (kind=real_kind) :: wmax
    integer :: maxiter
    real (kind=real_kind) :: deltatol,restol,deltaerr,reserr,rcond,min_rcond,anorm,dt3,alpha
    real (kind=real_kind) :: dw,dx,alpha_k,alphas(np,np)
    integer :: i,j,k,l,ie,nt
    integer :: nsafe,verbosity

#undef NEWTONCOND
#ifdef NEWTONCOND
    real*8 :: DLANGT  ! external
    real (kind=real_kind) :: work(2*nlev)
    integer :: iwork(nlev),info2
#endif

    verbosity = 1
    if (present(verbosity_in)) verbosity = verbosity_in

    call t_startf('compute_stage_value_dirk')

    if (theta_hydrostatic_mode) then
       itercount=0
       itererr=0
       do ie=nets,nete
          call phi_from_eos(hvcoord,elem(ie)%state%phis,elem(ie)%state%vtheta_dp(:,:,:,np1),&
               elem(ie)%state%dp3d(:,:,:,np1),elem(ie)%state%phinh_i(:,:,:,np1))
          elem(ie)%state%w_i(:,:,:,np1)=0
       enddo
       return
    endif

       
    ! dirk settings
    maxiter=20
#ifdef HOMMEXX_BFB_TESTING
    deltatol=1.0e-6_real_kind ! In BFB testing we can't converge due to calls of zeroulp
#else
    deltatol=1.0e-11_real_kind  ! exit if newton increment < deltatol
#endif

    !restol=1.0e-13_real_kind    ! exit if residual < restol  
    ! condition number and thus residual depends strongly on dt and min(dz)
    ! more work needed to exit iteration early based on residual error
    min_rcond=1.0e20_real_kind

    do ie=nets,nete
       phi_n0 = elem(ie)%state%phinh_i(:,:,:,np1)
       w_n0 = elem(ie)%state%w_i(:,:,:,np1)
       wmax=max(1d0,maxval(abs(w_n0)))

       phi_np1 => elem(ie)%state%phinh_i(:,:,:,np1)
       w_np1 => elem(ie)%state%w_i(:,:,:,np1)


       if (alphadt_n0.ne.0d0) then ! add dt*alpha*S(un0) to the rhs
          dt3=alphadt_n0
          nt=n0
          call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,nt), &
               elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%phinh_i(:,:,:,nt),pnh,   &
               exner,dpnh_dp_i,caller='dirkn0')
          w_n0(:,:,1:nlev)     = w_n0(:,:,1:nlev) + &
               dt3*g*(dpnh_dp_i(:,:,1:nlev)-1)

          call compute_gwphis(gwh_i,elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%v(:,:,:,:,nt),&
               elem(ie)%derived%gradphis,hvcoord)
          phi_n0(:,:,1:nlev) = phi_n0(:,:,1:nlev) + &
               dt3*g*elem(ie)%state%w_i(:,:,1:nlev,nt) -  dt3*gwh_i(:,:,1:nlev)
       end if


       if (alphadt_nm1.ne.0d0) then ! add dt*alpha*S(unm1) to the rhs
          dt3=alphadt_nm1
          nt=nm1
          call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,nt), &
               elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%phinh_i(:,:,:,nt),pnh,   &
               exner,dpnh_dp_i,caller='dirknm1')
          w_n0(:,:,1:nlev)     = w_n0(:,:,1:nlev) + &
               dt3*g*(dpnh_dp_i(:,:,1:nlev)-1)

          call compute_gwphis(gwh_i,elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%v(:,:,:,:,nt),&
               elem(ie)%derived%gradphis,hvcoord)
          phi_n0(:,:,1:nlev) = phi_n0(:,:,1:nlev) + &
               dt3*g*elem(ie)%state%w_i(:,:,1:nlev,nt) -  dt3*gwh_i(:,:,1:nlev)
       end if

       ! add just the wphis(np1) term to the RHS
       nt=np1
       call compute_gwphis(gwh_i,elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%v(:,:,:,:,nt),&
            elem(ie)%derived%gradphis,hvcoord)
       phi_n0(:,:,1:nlev) = phi_n0(:,:,1:nlev) -  dt2*gwh_i(:,:,1:nlev)

#if 0
       ! w(np1) as initial guess:
       phi_np1(:,:,1:nlev) =  phi_n0(:,:,1:nlev) +  dt2*g*w_np1(:,:,1:nlev)
#endif
#if 0
       ! phi_np1 as initial guess:
       w_np1(:,:,1:nlev) = (phi_np1(:,:,1:nlev) -  phi_n0(:,:,1:nlev) )/(dt2*g)
#endif
#if 0
       ! wh_i as initial guess:
       w_np1(:,:,1:nlev)=gwh_i(:,:,1:nlev)/g  
       phi_np1(:,:,1:nlev) =  phi_n0(:,:,1:nlev) +  dt2*g*w_np1(:,:,1:nlev)
#endif
#if 1
       ! use hydrostatic for initial guess
       call phi_from_eos(hvcoord,elem(ie)%state%phis,elem(ie)%state%vtheta_dp(:,:,:,np1),&
            elem(ie)%state%dp3d(:,:,:,np1),phi_np1)
       w_np1(:,:,1:nlev) = (phi_np1(:,:,1:nlev) -  phi_n0(:,:,1:nlev) )/(dt2*g)
#endif

       ! initial residual
       do k=1,nlev
          dphi(:,:,k)=phi_np1(:,:,k+1)-phi_np1(:,:,k)
          dphi_n0(:,:,k)=phi_n0(:,:,k+1)-phi_n0(:,:,k)
       enddo

       nsafe=0
       do k=1,nlev
          do j=1,np
             do i=1,np
                if ( dphi(i,j,k)  > -g) then
                   if (verbosity > 0) then
                   write(iulog,*) 'WARNING:IMEX limiting initial guess. ie,i,j,k=',ie,i,j,k
                   write(iulog,*) 'dphi(i,j,k)=  ',dphi(i,j,k)
                   endif
                   dphi(i,j,k)=-g
                   nsafe=1
                endif
             enddo
          enddo
       enddo
       if (nsafe==1) then
          ! in rare cases when limter was triggered, just recompute:
          do k=nlev,1,-1  ! scan
             phi_np1(:,:,k) = phi_np1(:,:,k+1)-dphi(:,:,k)
          enddo
          w_np1(:,:,1:nlev) = (phi_np1(:,:,1:nlev) -  phi_n0(:,:,1:nlev) )/(dt2*g)
       endif
       call pnh_and_exner_from_eos2(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,np1),elem(ie)%state%dp3d(:,:,:,np1),&
            dphi,pnh,exner,dpnh_dp_i,'dirk1')
       Fn(:,:,1:nlev) = w_np1(:,:,1:nlev) - &
            (w_n0(:,:,1:nlev) + g*dt2 * (dpnh_dp_i(:,:,1:nlev)-1))


       itercount=0
       do while (itercount < maxiter) 
          ! numerical J:
          !call get_dirk_jacobian(JacL,JacD,JacU,dt2,elem(ie)%state%dp3d(:,:,:,np1),dphi,pnh,0,1d-4,hvcoord,dpnh_dp_i,vtheta_dp) 
          ! analytic J:
          call get_dirk_jacobian(JacL,JacD,JacU,dt2,elem(ie)%state%dp3d(:,:,:,np1),dphi,pnh,1) 

          x(:,:,1:nlev) = -Fn(:,:,1:nlev)

          call solve_strict_diag_dominant_tridiag(JacL, JacD, JacU, x)

          do k = 1,nlev-1
             dphi(:,:,k) = dphi_n0(:,:,k) + &
                  dt2*g*((w_np1(:,:,k+1) - w_np1(:,:,k)) + &
                          (x(:,:,k+1) - x(:,:,k)))
          end do
          dphi(:,:,nlev) = dphi_n0(:,:,nlev) - dt2*g*(w_np1(:,:,nlev) + x(:,:,nlev))

          alphas = 1
          if (any(dphi(:,:,1:nlev) >= 0)) then
             do j=1,np
                do i=1,np
                   if (.not. any(dphi(i,j,1:nlev) >= 0)) cycle

                   ! Step halfway to the distance at which at least one dphi is 0.
                   alpha = 1
                   do k = 1,nlev
                      if (k < nlev) then
                         dx = x(i,j,k+1) - x(i,j,k)
                         dw = w_np1(i,j,k+1) - w_np1(i,j,k)
                      else
                         dx = -x(i,j,k)
                         dw = -w_np1(i,j,k)
                      end if
                      if (dx /= 0) then
                         alpha_k = -(dphi_n0(i,j,k) + dt2*g*dw)/(dt2*g*dx)
                         if (alpha_k >= 0) alpha = min(alpha, alpha_k)
                      end if
                   end do
                   alpha = alpha/2
                   alphas(i,j) = alpha

                   write(iulog,*) 'WARNING:IMEX is reducing step length from 1 to',alpha

                   do k = 1,nlev-1
                      dphi(i,j,k) = dphi_n0(i,j,k) + &
                           dt2*g*((w_np1(i,j,k+1) - w_np1(i,j,k)) + &
                           alpha*(x(i,j,k+1) - x(i,j,k)))
                   end do
                   dphi(i,j,nlev) = dphi_n0(i,j,nlev) - dt2*g*(w_np1(i,j,nlev) + alpha*x(i,j,nlev))
                end do
             end do
          end if
          do k = 1,nlev
             w_np1(:,:,k) = w_np1(:,:,k) + alphas*x(:,:,k)
          end do

          call pnh_and_exner_from_eos2(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,np1),&
               elem(ie)%state%dp3d(:,:,:,np1),dphi,pnh,exner,dpnh_dp_i,'dirk2')
          Fn(:,:,1:nlev) = w_np1(:,:,1:nlev) - (w_n0(:,:,1:nlev) + g*dt2 * (dpnh_dp_i(:,:,1:nlev)-1))

          ! this is not used in this loop, so move out of loop
          !reserr=maxval(abs(Fn))/(wmax*abs(dt2)) 
          deltaerr=maxval(abs(x))/wmax

          ! update iteration count and error measure
          itercount=itercount+1
          !if (reserr < restol) exit
          if (deltaerr<deltatol) exit
       end do ! end do for the do while loop

       ! update phi:
       phi_np1(:,:,1:nlev) =  phi_n0(:,:,1:nlev) +  dt2*g*w_np1(:,:,1:nlev)


       ! keep track of running  max iteraitons and max error (reset after each diagnostics output)
       max_itercnt=max(itercount,max_itercnt)
       max_deltaerr=max(deltaerr,max_deltaerr)
       reserr=maxval(abs(Fn))/(wmax*abs(dt2))
       max_reserr=max(reserr,max_reserr)
#ifdef NEWTONCOND
       min_rcond=min(rcond,min_rcond)
#endif
       itererr=min(max_reserr,max_deltaerr) ! passed back to ARKODE


       if (itercount >= maxiter) then
          if (verbosity > 0) then
          write(iulog,*) 'WARNING:IMEX solver failed b/c max iteration count was met',deltaerr,reserr
          end if
       end if
    end do ! end do for the ie=nets,nete loop
#ifdef NEWTONCOND
    if (hybrid%masterthread) print *,'max J condition number (mpi task0): ',1/min_rcond
#endif

    call t_stopf('compute_stage_value_dirk')
  end subroutine compute_stage_value_dirk


  subroutine get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,dphi,pnh,exact,&
     epsie,hvcoord,dpnh_dp_i,vtheta_dp)
  !================================================================================
  ! compute Jacobian of F(phi) = sum(dphi) +const + (dt*g)^2 *(1-dp/dpi) column wise
  ! with respect to phi
  !
  ! This subroutine forms the tridiagonal analytic Jacobian (we actually form the diagonal, sub-, and super-diagonal)
  ! J for use in a LApack tridiagonal LU factorization and solver to solve  J * x = -f either exactly or
  ! approximately
  !
  !  input:  
  !  exact jacobian:        dphi, dp3d, pnh
  !  matrix free jacobian:  dphi, dp3d, vtheta_dp, hvcoord, dpnh_dp_i
  !
  ! epsie == 1 means exact Jacobian, epsie ~= 1 means finite difference approximate jacobian
  ! exact,epsie,hvcoord,dpnh_dp,vtheta_dp,pnh,exner,exner_i are only needed as inputs
  ! if epsie ~=1
  !  
  ! The rule-of-thumb optimal epsie  is epsie = norm(elem)*sqrt(macheps)
  !
  !===================================================================================
    real (kind=real_kind), intent(out) :: JacD(np,np,nlev)
    real (kind=real_kind), intent(out) :: JacL(np,np,nlev-1),JacU(np,np,nlev-1)
    real (kind=real_kind), intent(in)  :: dp3d(np,np,nlev)
    real (kind=real_kind), intent(inout) :: pnh(np,np,nlev)
    real (kind=real_kind), intent(in)  :: dphi(np,np,nlev)
    real (kind=real_kind), intent(in)  :: dt2

    real (kind=real_kind), intent(in), optional :: epsie ! epsie is the differencing size in the approx. Jacobian
    real (kind=real_kind), intent(in), optional :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind), intent(in), optional :: vtheta_dp(np,np,nlev)
    type (hvcoord_t)     , intent(in),  optional    :: hvcoord

    integer, intent(in) :: exact

    ! local
    real (kind=real_kind) :: alpha1(np,np),alpha2(np,np)
    real (kind=real_kind) :: e(np,np,nlev),dphi_temp(np,np,nlev),exner(np,np,nlev)
    real (kind=real_kind) :: dpnh2(np,np,nlev),dpnh_dp_i_epsie(np,np,nlevp)
    real (kind=real_kind) :: ds(np,np,nlev),delta_mu(np,np,nlevp)
    real (kind=real_kind) :: a,b(np,np),ck(np,np),ckm1(np,np)
    !
    integer :: k,l,k2
    if (exact.eq.1) then ! use exact Jacobian
       ! this code will need to change when the equation of state is changed.
       ! add special cases for k==1 and k==nlev+1

       a  = (dt2*g)**2/(1-kappa)
       k  = 1 ! Jacobian row 1
       b  = a/dp3d(:,:,k)
       ck = pnh(:,:,k)/dphi(:,:,k)
       JacU(:,:,k) = 2*b*ck
       JacD(:,:,k) = 1 - JacU(:,:,k)
       ckm1 = ck
       do k = 2,nlev-1 ! Jacobian row k
          b  = 2*a/(dp3d(:,:,k-1) + dp3d(:,:,k))
          ck = pnh(:,:,k)/dphi(:,:,k)
          JacL(:,:,k-1) = b*ckm1
          JacU(:,:,k  ) = b*ck
          JacD(:,:,k  ) = 1 - JacL(:,:,k-1) - JacU(:,:,k)
          ckm1 = ck
       end do
       k  = nlev ! Jacobian row nlev
       b  = 2*a/(dp3d(:,:,k) + dp3d(:,:,k-1))
       ck = pnh(:,:,k)/dphi(:,:,k)
       JacL(:,:,k-1) = b*ckm1
       JacD(:,:,k  ) = 1 - JacL(:,:,k-1) - b*ck

    else ! use finite difference approximation to Jacobian with differencing size espie
      ! compute Jacobian of F(dphi) = phi +const + (dt*g)^2 *(1-dp/dpi) column wise
      ! if NEWTON_DPHI is defined, compute Jacobian of:
      !                   dF(dphi) = dphi +const + (dt*g)^2 *d(1-dp/dpi) 
      ! we only form the tridagonal entries and this code can easily be modified to
      ! accomodate sparse non-tridigonal and dense Jacobians, however, testing only
      ! the tridiagonal of a Jacobian is probably sufficient for testing purpose
      do k=1,nlev
        e=0
        e(:,:,k)=1
        dphi_temp(:,:,:) = dphi(:,:,:)
        ! dphi_tmp(k)= ( phi(k+1)-(phi(k)+eps ) = dphi(k)-eps
        ! dphi_tmp(k-1)= ( phi(k)+eps-(phi(k-1) ) = dphi(k-1)+eps
        dphi_temp(:,:,k) = dphi(:,:,k) - epsie*e(:,:,k)
        if (k>1) dphi_temp(:,:,k-1) = dphi(:,:,k-1) + epsie*e(:,:,k)

        if (theta_hydrostatic_mode) then
           !dpnh_dp_i_epsie(:,:,:)=1.d0
           delta_mu=0
        else
           call pnh_and_exner_from_eos2(hvcoord,vtheta_dp,dp3d,dphi_temp,pnh,exner,dpnh_dp_i_epsie,'get_dirk_jacobian')
           delta_mu(:,:,:)=(g*dt2)**2*(dpnh_dp_i(:,:,:)-dpnh_dp_i_epsie(:,:,:))/epsie
        end if

        JacD(:,:,k) = 1 +  delta_mu(:,:,k)
        if (k.eq.1) then
           JacL(:,:,k) =   delta_mu(:,:,k+1)
        elseif (k.eq.nlev) then
           JacU(:,:,k-1) = delta_mu(:,:,k-1)
        else
           JacL(:,:,k)   = delta_mu(:,:,k+1)
           JacU(:,:,k-1) = delta_mu(:,:,k-1)
        end if

      end do
    end if
  end subroutine get_dirk_jacobian

  subroutine test_imex_jacobian(elem,hybrid,hvcoord,tl,nets,nete)
    ! the following code compares the analytic vs exact imex Jacobian
    ! can test over more elements if desired
    type(element_t)   , intent(in) :: elem(:)
    type(hybrid_t)    , intent(in) :: hybrid
    type (hvcoord_t)  , intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl  
    integer                        :: nets,nete

    real (kind=real_kind) :: JacD(np,np,nlev)  , JacL(np,np,nlev-1)
    real (kind=real_kind) :: JacU(np,np,nlev-1)
    real (kind=real_kind) :: Jac2D(np,np,nlev)  , Jac2L(np,np,nlev-1)
    real (kind=real_kind) :: Jac2U(np,np,nlev-1)

    real (kind=real_kind) :: dp3d(np,np,nlev), phis(np,np)
    real (kind=real_kind) :: phi_i(np,np,nlevp)
    real (kind=real_kind) :: dphi(np,np,nlev)
    real (kind=real_kind) :: vtheta_dp(np,np,nlev)
    real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind) :: exner(np,np,nlev)
    real (kind=real_kind) :: pnh(np,np,nlev), pnh_i(np,np,nlevp)
    real (kind=real_kind) :: norminfJ0(np,np)

    real (kind=real_kind) :: dt,epsie,jacerrorvec(6),minjacerr
    integer :: k,ie,i,j
    minjacerr=0
    if (hybrid%masterthread) write(iulog,*)'Running IMEX Jacobian unit test...'
    do ie=nets,nete
       do k=1,nlev
          dp3d(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
       enddo
       vtheta_dp(:,:,:) = elem(ie)%state%vtheta_dp(:,:,:,tl%n0)
       phi_i(:,:,:)         = elem(ie)%state%phinh_i(:,:,:,tl%n0)
       phis(:,:)          = elem(ie)%state%phis(:,:)
       call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_i,&
               pnh,exner,dpnh_dp_i,pnh_i_out=pnh_i)

       dt=100.0

       do k=1,nlev
          dphi(:,:,k)=phi_i(:,:,k+1)-phi_i(:,:,k)
       enddo
       call get_dirk_jacobian(JacL,JacD,JacU,dt,dp3d,dphi,pnh,1)

      ! compute infinity norm of the initial Jacobian 
       norminfJ0=0.d0
       do k=1,nlev
       do i=1,np
       do j=1,np
          if (k.eq.1) then
            norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacD(i,j,k))+abs(JacU(i,j,k))))
          elseif (k.eq.nlev) then
            norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(i,j,k-1))+abs(JacD(i,j,k))))
          else
            norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(i,j,k-1))+abs(JacD(i,j,k))+ &
              abs(JacU(i,j,k))))
          end if
      end do
      end do
      end do

       jacerrorvec=0
       do j=1,6
          ! compute numerical jacobian with 5 different epsilons and take the smallest error
          ! the function that we are trying to estimate the numerical jacobian of
          ! phi + const + (dt*g)^2 (1-dp/dpi)is dt*g^2 dpnh/dpi = O(10,000)
          ! =================================================================
          ! PLEASE NOTE:  We take the minimum rather than the maximum error
          ! since the error of finite difference approximations to the 
          ! Jacobian in finite precision arithmetic first decrease and then
          ! increase as the perturbation size decreases due to round-off error.
          ! So to test the "correctness" of the exact Jacobian, we need to find
          ! that the sweetspot where the finite difference error is minimized is
          ! =================================================================
          epsie=10.d0/(10.d0)**(j+1)
          do k=1,nlev
             dphi(:,:,k)=phi_i(:,:,k+1)-phi_i(:,:,k)
          enddo
          call get_dirk_jacobian(Jac2L,Jac2D,Jac2U,dt,dp3d,dphi,pnh,0,&
             epsie,hvcoord,dpnh_dp_i,vtheta_dp)

          if (maxval(abs(JacD(:,:,:)-Jac2D(:,:,:))) > jacerrorvec(j)) then 
             jacerrorvec(j) = maxval(abs(JacD(:,:,:)-Jac2D(:,:,:)))
          end if
          if (maxval(abs(JacL(:,:,:)-Jac2L(:,:,:))) > jacerrorvec(j)) then
             jacerrorvec(j) = maxval(abs(JacL(:,:,:)-Jac2L(:,:,:)))
          end if
          if (maxval(abs(JacU(:,:,:)-Jac2U(:,:,:))) > jacerrorvec(j)) then
             jacerrorvec(j) = maxval(abs(JacU(:,:,:)-Jac2U(:,:,:)))
          end if
       end do
       minjacerr = max( minval(jacerrorvec(:))/maxval(norminfJ0)   ,minjacerr)
       !     minjacerr = minval(jacerrorvec(:))
    end do
    if (minjacerr > 1.0e-3_real_kind) then 
       write(iulog,*)'WARNING:  Analytic and exact Jacobian differ by ', minjacerr
       write(iulog,*)'Please check that the IMEX exact Jacobian in eos.F90 is actually exact'
    else
       if (hybrid%masterthread) write(iulog,*)&
            'PASS. max error of analytic and exact Jacobian: ',minjacerr
    end if

    call test_tridiag_solver(hybrid, JacL, JacD, JacU)
  end subroutine test_imex_jacobian

  subroutine solve_strict_diag_dominant_tridiag(dl, d, du, x)
    ! Solve a strictly diagonally dominant tridiagonal system. On
    ! input, x is the RHS; on output, the solution. dl contains rows
    ! 2:nlev of the lower diagonal; d, the diagonal; du, rows 1:nlev-1
    ! of the upper diagonal. d is overwritten.

    use dimensions_mod, only: npsq

    real (kind=real_kind), intent(in) :: dl(npsq,nlev-1), du(npsq,nlev-1)
    real (kind=real_kind), intent(inout) :: d(npsq,nlev), x(npsq,nlev)

    real (kind=real_kind) :: dlk(npsq)
    integer :: k

    do k = 1,nlev-1
       dlk = dl(:,k) / d(:,k)
       d(:,k+1) = d(:,k+1) - dlk * du(:,k)
       x(:,k+1) = x(:,k+1) - dlk * x(:,k)
    end do
    x(:,nlev) = x(:,nlev) / d(:,nlev)
    do k = nlev-1,1,-1
       x(:,k) = (x(:,k) - du(:,k) * x(:,k+1)) / d(:,k)
    end do
  end subroutine solve_strict_diag_dominant_tridiag

  subroutine test_tridiag_solver(hybrid, dlc, dc, duc)
    type(hybrid_t), intent(in) :: hybrid
    real (kind=real_kind), intent(in) :: dlc(np,np,nlev-1), dc(np,np,nlev), duc(np,np,nlev-1)

    real (kind=real_kind) :: &
         dl(np,np,nlev-1), d(np,np,nlev), du(np,np,nlev-1), &
         b(np,np,nlev), x(np,np,nlev), y(np,np,nlev), &
         relerr(np,np)

    integer :: i, j, k
    
    ! Make a RHS.
    dl = dlc; d = dc; du = duc
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             b(i,j,k) = k + real(i, real_kind)/2 - real(j,real_kind)/3
          end do
       end do
    end do
    x = b

    call solve_strict_diag_dominant_tridiag(dl, d, du, x)

    ! y = A*x
    y = 0
    do k = 1,nlev-1
       y(:,:,k+1) = dlc(:,:,k) * x(:,:,k)
    end do
    do k = 1,nlev
       y(:,:,k) = y(:,:,k) + dc(:,:,k) * x(:,:,k)
    end do
    do k = 1,nlev-1
       y(:,:,k) = y(:,:,k) + duc(:,:,k) * x(:,:,k+1)
    end do

    ! Compare y and b
    relerr = maxval(abs(y - b), 3)/maxval(abs(b), 3)
    if (hybrid%masterthread .and. maxval(relerr) > 1e5*epsilon(1.0_real_kind)) then
       write(iulog,*) 'FAIL test_tridiag_solver', maxval(relerr)
    end if
  end subroutine test_tridiag_solver

end module imex_mod
