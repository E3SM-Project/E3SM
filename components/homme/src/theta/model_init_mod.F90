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

  use element_mod, only: element_t
  use derivative_mod, only: derivative_t
  use hybvcoord_mod, 	only: hvcoord_t
  use hybrid_mod,         only: hybrid_t
  use dimensions_mod, only: np,nlev,nlevp,nelemd
  use eos          ,  only: get_pnh_and_exner,get_dry_phinh,get_dirk_jacobian
  use element_ops,   only: get_kappa_star
  use kinds,         only: real_kind,iulog
  use control_mod,   only: qsplit,theta_hydrostatic_mode
  use time_mod,      only: timelevel_qdp, timelevel_t

  implicit none
  
contains

  subroutine vertical_mesh_init2(elem, nets, nete, hybrid, hvcoord)

    ! additional solver specific initializations (called from prim_init2)

    type (element_t),			intent(inout), target :: elem(:)! array of element_t structures
    integer,				intent(in) :: nets,nete		! start and end element indices
    type (hybrid_t),			intent(in) :: hybrid		! mpi/omp data struct
    type (hvcoord_t),			intent(inout)	:: hvcoord	! hybrid vertical coord data struct

  end subroutine vertical_mesh_init2



  subroutine model_init2(elem,hybrid,deriv,hvcoord,tl,nets,nete )

    type(element_t)   , intent(in) :: elem(:)
    type(hybrid_t)    , intent(in) :: hybrid
    type(derivative_t), intent(in) :: deriv
    type (hvcoord_t)  , intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl
    integer                        :: nets,nete

    ! unit test for analytic jacobian used by IMEX methods
    if (.not. theta_hydrostatic_mode) &
         call test_imex_jacobian(elem,hybrid,hvcoord,tl,nets,nete)

    ! other theta specific model initialization should go here
  
  end subroutine 




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
  
  real (kind=real_kind) :: kappa_star(np,np,nlev),kappa_star_i(np,np,nlevp)
  real (kind=real_kind) :: dp3d(np,np,nlev), phi(np,np,nlev), phis(np,np)
  real (kind=real_kind) :: theta_dp_cp(np,np,nlev), dpnh_dp(np,np,nlev),dpnh(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev),pnh_i(np,np,nlevp)
  real (kind=real_kind) :: norminfJ0(np,np)
  
  real (kind=real_kind) :: dt,epsie,maxjacerrorvec(5),maxjacerr
  integer :: k,ie,qn0,i,j
  maxjacerr=0
  if (hybrid%masterthread) write(iulog,*)'Running IMEX Jacobian unit test...'
  do ie=nets,nete
     do k=1,nlev
        dp3d(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
     enddo
     theta_dp_cp(:,:,:) = elem(ie)%state%theta_dp_cp(:,:,:,tl%n0)
     phi(:,:,:)         = elem(ie)%state%phinh(:,:,:,tl%n0)
     phis(:,:)          = elem(ie)%state%phis(:,:)
     call TimeLevel_Qdp(tl, qsplit, qn0)
     call get_kappa_star(kappa_star,elem(ie)%state%Q(:,:,:,1))
     if (theta_hydrostatic_mode) then
        dpnh_dp(:,:,:)=1.d0
     else
        call get_pnh_and_exner(hvcoord,theta_dp_cp,dp3d,phi,elem(ie)%state%phis,&
             kappa_star,pnh,dpnh,exner,pnh_i_out=pnh_i)
        dpnh_dp(:,:,:) = dpnh(:,:,:)/dp3d(:,:,:)
     end if
     do k=1,nlev
        kappa_star_i(:,:,k+1) = 0.5D0* (kappa_star(:,:,k+1)+kappa_star(:,:,k))
     end do

     dt=100.d0
     call get_dirk_jacobian(JacL,JacD,JacU,dt,dp3d,phi,phis,kappa_star_i,pnh_i,1)

    ! compute infinity norm of the initial Jacobian 
     norminfJ0=0.d0
     do i=1,np
     do j=1,np
       do k=1,nlev
        if (k.eq.1) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacD(k,i,j))+abs(JacU(k,i,j))))
        elseif (k.eq.nlev) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k,i,j))+abs(JacD(k,i,j))))
        else
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k,i,j))+abs(JacD(k,i,j))+ &
            abs(JacU(k,i,j))))
        end if
      end do
    end do
    end do

     
     
     maxjacerrorvec=0.d0
     do j=1,5
        ! compute numerical jacobian with 5 different epsilons and take the smallest error
        ! the function that we are trying to estimate the numerical jacobian of
        ! phi + const + (dt*g)^2 (1-dp/dpi)is dt*g^2 dpnh/dpi = O(10,000)
        epsie=10.d0/(10.d0)**(j+1)
        call get_dirk_jacobian(Jac2L,Jac2D,Jac2U,dt,dp3d,phi,phis,kappa_star_i,pnh_i,0,&
             epsie,hvcoord,dpnh_dp,theta_dp_cp,kappa_star,pnh,exner)
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
     maxjacerr = max( minval(maxjacerrorvec(:))/maxval(norminfJ0)   ,maxjacerr)
  end do
  if (maxjacerr > 1e-3) then 
     write(iulog,*)'WARNING:  Analytic and exact Jacobian differ by ', maxjacerr
     write(iulog,*)'Please check that the IMEX exact Jacobian in eos.F90 is actually exact'
  else
     if (hybrid%masterthread) write(iulog,*)&
          'PASS. max error of analytic and exact Jacobian: ',maxjacerr
  end if
end subroutine test_imex_jacobian



end module 
