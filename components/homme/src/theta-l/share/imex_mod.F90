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
#ifdef XX_BFB_TESTING
  use bfb_mod,            only: tridiag_diagdom_bfb_a1x1
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

  subroutine compute_stage_value_dirk(n0,np1,alphadt,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,itercount,itererr,nm1)
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
    ! w(np1) = w_0 + alphadt1 * SW(n0)  + alphadt1*SW(nm1) +  dt2*SW(np1)
    ! phi(np1) = phi_0 + alphadt1 * SPHI(n0)  + alphadt1*SPHI(nm1) +  dt2*SPHI(np1)
    !
    ! SW(nt) = g*dpnh_dp_i(nt)-1
    ! SPHI(nt) = g*w(nt) -g*a(k) u(nt) dot grad_phis
    !
    ! pnh_and_exner_from_eos()   used to compute dpnh_dp_i(nt) in SW(nt)
    ! compute_gwphis()            used to compute g*a(k) u(nt) dot grad_phis  in SPHI(nt)   
    !
    ! We then precompute:
    !  w_rhs = w_0 + alphadt1 * SW(n0)  + alphadt1*SW(nm1) 
    !  phi_rhs = w_0 + alphadt1 * SPHI(n0)  + alphadt1*SPHI(nm1) -  dt2*compute_gwphi(np1)
    ! and solve, via Newton iteraiton:
    !   w(np1) = w_rhs + dt2*SW(np1)
    !   phi(np1) = phi_rhs + dt2*g*w(np1)
    !
    !===================================================================================

    integer, intent(in) :: n0,np1,qn0,nets,nete
    real (kind=real_kind), intent(in) :: dt2
    integer :: itercount
    real (kind=real_kind) :: itererr
    real (kind=real_kind), intent(in) :: alphadt

    type (hvcoord_t)     , intent(in) :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer, optional :: nm1


    ! local
    real (kind=real_kind), pointer, dimension(:,:,:)   :: phi_np1
    real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
    real (kind=real_kind), pointer, dimension(:,:,:)   :: vtheta_dp
    real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
    real (kind=real_kind) :: JacU(nlev-1,np,np), JacU2(nlev-2,np,np)
    real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
    real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
    real (kind=real_kind) :: w_n0(np,np,nlevp)    
    real (kind=real_kind) :: dphi(np,np,nlev)    
    real (kind=real_kind) :: dphi_n0(np,np,nlev)    
    real (kind=real_kind) :: phi_n0(np,np,nlevp)    
    real (kind=real_kind) :: delta_phi(np,np,nlevp)    ! phi_np1-phi_n0
    real (kind=real_kind) :: Ipiv(nlev,np,np)
    real (kind=real_kind) :: Fn(np,np,nlev),x(nlev,np,np)
    real (kind=real_kind) :: gwh_i(np,np,nlevp)  ! w hydrostatic
    real (kind=real_kind) :: v_i(np,np,2,nlevp)  ! w hydrostatic

    real (kind=real_kind) :: Jac2D(nlev,np,np)  , Jac2L(nlev-1,np,np)
    real (kind=real_kind) :: Jac2U(nlev-1,np,np)


    real (kind=real_kind) :: wgdtmax
    integer :: maxiter
    real*8 :: deltatol,restol,deltaerr,reserr,rcond,min_rcond,anorm,dt3,alpha

    integer :: i,j,k,l,ie,info(np,np),nt
    integer :: nsafe
#undef NEWTONCOND
#ifdef NEWTONCOND
    real*8 :: DLANGT  ! external
    real (kind=real_kind) :: work(2*nlev)
    integer :: iwork(nlev),info2
#endif

    call t_startf('compute_stage_value_dirk')

    ! dirk settings
    maxiter=20
    deltatol=1.0e-13_real_kind  ! exit if newton increment < deltatol
    !restol=1.0e-13_real_kind    ! exit if residual < restol  
    ! condition number and thus residual depends strongly on dt and min(dz)
    ! more work needed to exit iteration early based on residual error
    delta_phi(:,:,nlevp)=0
    min_rcond=1.0e20_real_kind

    do ie=nets,nete
       w_n0 = elem(ie)%state%w_i(:,:,:,np1)
       wgdtmax=max(1d0,maxval(abs(w_n0)))*abs(dt2)*g
       ! approximate the initial error of f(x) \approx 0
       vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,np1)
       phi_np1 => elem(ie)%state%phinh_i(:,:,:,np1)

       phi_n0 = phi_np1

       if (alphadt.ne.0d0) then ! add dt*alpha*S(un0) to the rhs
          dt3=alphadt
          nt=n0
          call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,nt), &
               elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%phinh_i(:,:,:,nt),pnh,   &
               exner,dpnh_dp_i,caller='dirk0')
          w_n0(:,:,1:nlev)     = w_n0(:,:,1:nlev) + &
               dt3*g*(dpnh_dp_i(:,:,1:nlev)-1d0)

          call compute_gwphis(gwh_i,elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%v(:,:,:,:,nt),&
               elem(ie)%derived%gradphis,hvcoord)
          phi_n0(:,:,1:nlev) = phi_n0(:,:,1:nlev) + &
               dt3*g*elem(ie)%state%w_i(:,:,1:nlev,nt) -  dt3*gwh_i(:,:,1:nlev)
       end if


       if (present(nm1)) then ! add dt*alpha*S(unm1) to the rhs
          dt3=alphadt
          nt=nm1
          call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,nt), &
               elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%phinh_i(:,:,:,nt),pnh,   &
               exner,dpnh_dp_i,caller='dirk0')
          w_n0(:,:,1:nlev)     = w_n0(:,:,1:nlev) + &
               dt3*g*(dpnh_dp_i(:,:,1:nlev)-1d0)

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

       dp3d  => elem(ie)%state%dp3d(:,:,:,np1)
#if 0
       ! use hydrostatic for initial guess
       call phi_from_eos(hvcoord,elem(ie)%state%phis,vtheta_dp,dp3d,phi_np1)
#endif

       ! newton iteration will iterate of d(phi)/deta.  initialize:
       do k=1,nlev
          dphi_n0(:,:,k)=phi_n0(:,:,k+1)-phi_n0(:,:,k)
          dphi(:,:,k)=phi_np1(:,:,k+1)-phi_np1(:,:,k)
       enddo

       do k=1,nlev
          do j=1,np
             do i=1,np
                if ( dphi(i,j,k)  > -g) then
                   write(iulog,*) 'WARNING:IMEX inital PHI is bad, delta z < 1m. ie,i,j,k=',ie,i,j,k
                   write(iulog,*) 'dphi(i,j,k)=  ',dphi(i,j,k)
                   dphi(i,j,k)=-g
                endif
             enddo
          enddo
       enddo

       ! initial residual
       call pnh_and_exner_from_eos2(hvcoord,vtheta_dp,dp3d,dphi,pnh,exner,dpnh_dp_i,'dirk1')
       elem(ie)%state%w_i(:,:,1:nlev,np1) = w_n0(:,:,1:nlev) - g*dt2 * &
            (1.0-dpnh_dp_i(:,:,1:nlev))
       do k=1,nlev
          Fn(:,:,k) = (phi_np1(:,:,k) - phi_n0(:,:,k)) - dt2*g*(elem(ie)%state%w_i(:,:,k,np1))
       enddo


       itercount=0
       do while (itercount < maxiter) 

          info(:,:) = 0
          ! numerical J:
          !call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,dphi,pnh,0,1d-4,hvcoord,dpnh_dp_i,vtheta_dp) 
          ! analytic J:
          call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,dphi,pnh,1) 

          do i=1,np
             do j=1,np
                x(1:nlev,i,j) = -Fn(i,j,1:nlev)  !+Fn(i,j,nlev+1:2*nlev,1)/(g*dt2))
#ifdef XX_BFB_TESTING
                ! Note: the C function is designed to accept both single and double precision,
                !       so we need to pass also the size of a real (last argument)
                call tridiag_diagdom_bfb_a1x1(nlev, JacL(:,i,j), jacD(:,i,j), jacU(:,i,j), x(:,i,j),INT(SIZEOF(JacL)/SIZE(JacL),4))
#else
#ifdef NEWTONCOND
                ! nlev condition number: 500e3 with phi, 850e3 with dphi
                anorm=DLANGT('1-norm', nlev, JacL(:,i,j),jacD(:,i,j),jacU(:,i,j))
                call DGTTRF(nlev, JacL(:,i,j), JacD(:,i,j),JacU(:,i,j),JacU2(:,i,j), Ipiv(:,i,j), info )
                call DGTCON('1',nlev,JacL(:,i,j),JacD(:,i,j),JacU(:,i,j),jacU2(:,i,j),Ipiv(:,i,j),&
                     ANORM, RCOND, WORK, IWORK, info2 )
#else
                call DGTTRF(nlev, JacL(:,i,j), JacD(:,i,j),JacU(:,i,j),JacU2(:,i,j), Ipiv(:,i,j), info(i,j) )
#endif
                ! Tridiagonal solve
                call DGTTRS( 'N', nlev,1, JacL(:,i,j), JacD(:,i,j), JacU(:,i,j), JacU2(:,i,j), Ipiv(:,i,j),x(:,i,j), nlev, info(i,j) )
#endif
                ! update approximate solution of phi
                do k=1,nlev-1
                   dphi(i,j,k)=dphi(i,j,k) + x(k+1,i,j)-x(k,i,j)
                enddo
                dphi(i,j,nlev)=dphi(i,j,nlev) + (0 - x(nlev,i,j) )

                alpha = 0
                do nsafe=1,8
                   if (all(dphi(i,j,1:nlev) < 0 ))  exit
                   ! remove the last netwon increment, try reduced increment
                   alpha = 1.0_real_kind/(2**nsafe)
                   do k=1,nlev-1
                      dphi(i,j,k)=dphi(i,j,k) - (x(k+1,i,j)-x(k,i,j))*alpha
                   enddo
                   dphi(i,j,nlev)=dphi(i,j,nlev) - (0-x(nlev,i,j))*alpha
                enddo
                if (nsafe>1) write(iulog,*) 'WARNING:IMEX reducing newton increment, nsafe=',nsafe
                ! if nsafe>8, code will crash in next call to pnh_and_exner_from_eos

                phi_np1(i,j,:) = phi_np1(i,j,:) + (1 - alpha)*x(:,i,j)
             end do
          end do
          do k=1,nlev
             dphi(:,:,k) = phi_np1(:,:,k+1) - phi_np1(:,:,k)
          end do
          call pnh_and_exner_from_eos2(hvcoord,vtheta_dp,dp3d,dphi,pnh,exner,dpnh_dp_i,'dirk2')

          ! update approximate solution of w
          elem(ie)%state%w_i(:,:,1:nlev,np1) = w_n0(:,:,1:nlev) - g*dt2 * &
               (1.0-dpnh_dp_i(:,:,1:nlev))
          do k=1,nlev
             Fn(:,:,k) = (phi_np1(:,:,k) - phi_n0(:,:,k)) - dt2*g*(elem(ie)%state%w_i(:,:,k,np1))
          enddo
          reserr=maxval(abs(Fn))/wgdtmax   ! residual error in phi tendency, relative to source term gw

          deltaerr=0
          do k=1,nlev
             ! delta residual:
             do i=1,np
                do j=1,np
                   deltaerr=max(deltaerr, abs(x(k,i,j))/max(g,abs(dphi_n0(i,j,k))) )
                enddo
             enddo
          enddo

          ! update iteration count and error measure
          itercount=itercount+1
          !if (reserr < restol) exit
          if (deltaerr<deltatol) exit
       end do ! end do for the do while loop

       ! keep track of running  max iteraitons and max error (reset after each diagnostics output)
       max_itercnt=max(itercount,max_itercnt)
       max_deltaerr=max(deltaerr,max_deltaerr)
       max_reserr=max(reserr,max_reserr)
#ifdef NEWTONCOND
       min_rcond=min(rcond,min_rcond)
#endif
       itererr=min(max_reserr,max_deltaerr) ! passed back to ARKODE


! For BFB testing on GPU, we zero out some bits in calls to pow functions,
! to preserve bfb GPU vs F90, so we *expect* issues in convergence.
! In any other build, this should be concerning though.
#if !(defined(XX_BFB_TESTING) && defined(CUDA_BUILD))
       if (itercount >= maxiter) then
          write(iulog,*) 'WARNING:IMEX solver failed b/c max iteration count was met',deltaerr,reserr
          do k=1,nlev
             i=1 ; j=1
             !print *,k,( abs(Fn(i,j,k))/wgdtmax
          enddo
       end if
#endif
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
    real (kind=real_kind), intent(out) :: JacD(nlev,np,np)
    real (kind=real_kind), intent(out) :: JacL(nlev-1,np,np),JacU(nlev-1,np,np)
    real (kind=real_kind), intent(in)    :: dp3d(np,np,nlev)
    real (kind=real_kind), intent(inout) :: pnh(np,np,nlev)
    real (kind=real_kind), intent(in) :: dphi(np,np,nlev)
    real (kind=real_kind), intent(in)    :: dt2

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
       JacU(k,:,:) = 2*b*ck
       JacD(k,:,:) = 1 - JacU(k,:,:)
       ckm1 = ck
       do k = 2,nlev-1 ! Jacobian row k
          b  = 2*a/(dp3d(:,:,k-1) + dp3d(:,:,k))
          ck = pnh(:,:,k)/dphi(:,:,k)
          JacL(k-1,:,:) = b*ckm1
          JacU(k,  :,:) = b*ck
          JacD(k,  :,:) = 1 - JacL(k-1,:,:) - JacU(k,:,:)
          ckm1 = ck
       end do
       k  = nlev ! Jacobian row nlev
       b  = 2*a/(dp3d(:,:,k) + dp3d(:,:,k-1))
       ck = pnh(:,:,k)/dphi(:,:,k)
       JacL(k-1,:,:) = b*ckm1
       JacD(k  ,:,:) = 1 - JacL(k-1,:,:) - b*ck

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

        JacD(k,:,:) = 1 +  delta_mu(:,:,k)
        if (k.eq.1) then
           JacL(k,:,:) =   delta_mu(:,:,k+1)
        elseif (k.eq.nlev) then
           JacU(k-1,:,:) = delta_mu(:,:,k-1)
        else
           JacL(k,:,:)   = delta_mu(:,:,k+1)
           JacU(k-1,:,:) = delta_mu(:,:,k-1)
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

    real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
    real (kind=real_kind) :: JacU(nlev-1,np,np)
    real (kind=real_kind) :: Jac2D(nlev,np,np)  , Jac2L(nlev-1,np,np)
    real (kind=real_kind) :: Jac2U(nlev-1,np,np)

    real (kind=real_kind) :: dp3d(np,np,nlev), phis(np,np)
    real (kind=real_kind) :: phi_i(np,np,nlevp)
    real (kind=real_kind) :: dphi(np,np,nlev)
    real (kind=real_kind) :: vtheta_dp(np,np,nlev)
    real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind) :: exner(np,np,nlev)
    real (kind=real_kind) :: pnh(np,np,nlev),	pnh_i(np,np,nlevp)
    real (kind=real_kind) :: norminfJ0(np,np)

    real (kind=real_kind) :: dt,epsie,jacerrorvec(6),minjacerr
    integer :: k,ie,qn0,i,j
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
       call TimeLevel_Qdp(tl, qsplit, qn0)
       call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_i,&
               pnh,exner,dpnh_dp_i,pnh_i_out=pnh_i)

       dt=100.0

       do k=1,nlev
          dphi(:,:,k)=phi_i(:,:,k+1)-phi_i(:,:,k)
       enddo
       call get_dirk_jacobian(JacL,JacD,JacU,dt,dp3d,dphi,pnh,1)

      ! compute infinity norm of the initial Jacobian 
       norminfJ0=0.d0
       do i=1,np
       do j=1,np
         do k=1,nlev
          if (k.eq.1) then
            norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacD(k,i,j))+abs(JacU(k,i,j))))
          elseif (k.eq.nlev) then
            norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k-1,i,j))+abs(JacD(k,i,j))))
          else
            norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k-1,i,j))+abs(JacD(k,i,j))+ &
              abs(JacU(k,i,j))))
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
  end subroutine test_imex_jacobian

end module imex_mod
