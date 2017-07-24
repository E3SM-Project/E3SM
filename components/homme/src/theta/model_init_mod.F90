#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  model specific initialization code called from prim_init2
!  (threaded initialization code)
!
!  most models do nothing.  introduced for preqx_acc to initialize
!  GPU related data
!
module model_init_mod
  implicit none

contains

  subroutine model_init2( elem , deriv,hvcoord,tl )
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t
    use	hybvcoord_mod, 	only: hvcoord_t
    use dimensions_mod, only: np,nlev,nlevp,nelemd
    use eos          ,  only: get_pnh_and_exner,get_dry_phinh,get_dirk_jacobian
    use element_ops,   only: get_kappa_star
    use kinds,         only: real_kind
    use control_mod,   only: qsplit,theta_hydrostatic_mode
    use time_mod,      only: timelevel_qdp, timelevel_t

    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    type (hvcoord_t)  , intent(in), optional :: hvcoord
    type (TimeLevel_t)   , intent(in)            :: tl

    real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
    real (kind=real_kind) :: JacU(nlev-1,np,np)
    real (kind=real_kind) :: Jac2D(nlev,np,np)  , Jac2L(nlev-1,np,np)
    real (kind=real_kind) :: Jac2U(nlev-1,np,np)

    real (kind=real_kind) :: kappa_star(np,np,nlev),kappa_star_i(np,np,nlevp)
    real (kind=real_kind) :: dp3d(np,np,nlev), phi(np,np,nlev), phis(np,np)
    real (kind=real_kind) :: theta_dp_cp(np,np,nlev), dpnh_dp(np,np,nlev),dpnh(np,np,nlev)
    real (kind=real_kind) :: exner(np,np,nlev), exner_i(np,np,nlevp)
    real (kind=real_kind) :: pnh(np,np,nlev),	pnh_i(np,np,nlevp)

    real (kind=real_kind) :: dt,epsie,maxjacerrorvec(5),maxjacerr
    integer :: k,ie,qn0,j
    ! the following code compares the analytic vs exact imex Jacobian
    ! can test over more elements if desired
    maxjacerr=1d-1
    do ie=1,1
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
      do k=1,nlev
        dp3d(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
      enddo
      theta_dp_cp(:,:,:) = elem(ie)%state%theta_dp_cp(:,:,:,tl%n0)
      phi(:,:,:)         = elem(ie)%state%phi(:,:,:,tl%n0)
      phis(:,:)          = elem(ie)%state%phis(:,:)
      call TimeLevel_Qdp(tl, qsplit, qn0)
      call get_kappa_star(kappa_star,elem(ie)%state%Qdp(:,:,:,1,qn0),dp3d)
      if (theta_hydrostatic_mode) then
        dpnh_dp(:,:,:)=1.d0
      else
        call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
          kappa_star,pnh,dpnh,exner,exner_i,pnh_i)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
      end if
      do k=1,nlev
        kappa_star_i(:,:,k+1) = 0.5D0* (kappa_star(:,:,k+1)+kappa_star(:,:,k))
      end do

      dt=100.d0
      maxjacerrorvec=0.d0
      do j=1,5
        JacL=0.d0
        JacD=0.d0
        JacU=0.d0
        call get_dirk_jacobian(JacL,JacD,JacU,dt,dp3d,phi,phis,kappa_star_i,pnh_i,1)
        epsie=1.d0/(10.d0)**(j+1)
        call get_dirk_jacobian(Jac2L,Jac2D,Jac2U,dt,dp3d,phi,phis,kappa_star_i,pnh_i,0,&
          epsie,hvcoord,dpnh_dp,theta_dp_cp,kappa_star,pnh,exner,exner_i)
        if (maxval(abs(JacD(:,:,:)-Jac2D(:,:,:))) > maxjacerrorvec(j)) then 
          maxjacerrorvec(j) = maxval(abs(JacD(:,:,:)-Jac2D(:,:,:)))
        end if
       	if (maxval(abs(JacL(:,:,:)-Jac2L(:,:,:))) > maxjacerrorvec(j)) then
       	  maxjacerrorvec(j) = maxval(abs(JacL(:,:,:)-Jac2L(:,:,:)))
    	end if
       	if (maxval(abs(JacU(:,:,:)-Jac2U(:,:,:))) > maxjacerrorvec(j)) then
       	  maxjacerrorvec(j) = maxval(abs(JacU(:,:,:)-Jac2U(:,:,:)))
        end if
      end do

      maxjacerr = min(minval(maxjacerrorvec(:)),maxjacerr)
    end do 
    ! do nothing
    if (maxjacerr > 1.d-1) then 
      print *, 'WARNING:  Analytic and exact Jacobian differ by ', maxjacerr
      print *, 'Please check that the IMEX exact Jacobian in eos.F90 is actually exact'
    else
      print *, 'CONGRATULATIONS: Tridiagonal of analytic and exact Jacobian differ by ', maxjacerr
    end if
  
  end subroutine 

end module 
