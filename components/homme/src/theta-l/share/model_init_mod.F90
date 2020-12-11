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
!  2018/10 MT: adding TOM pressure-based sponge layer dissipation from P. Lauritzen
!
module model_init_mod

  use element_mod,        only: element_t
  use derivative_mod,     only: derivative_t,gradient_sphere
  use hybvcoord_mod, 	  only: hvcoord_t
  use hybrid_mod,         only: hybrid_t
  use dimensions_mod,     only: np,nlev,nlevp
  use element_ops,        only: set_theta_ref
  use element_state,      only: timelevels, nu_scale_top, nlev_tom
  use viscosity_mod,      only: make_c0_vector
  use kinds,              only: real_kind,iulog
  use control_mod,        only: qsplit,theta_hydrostatic_mode
  use time_mod,           only: timelevel_qdp, timelevel_t
  use physical_constants, only: g, TREF, Rgas, kappa
  use imex_mod,           only: test_imex_jacobian
  use eos,                only: phi_from_eos
 
  implicit none
  
contains

  subroutine model_init2(elem,hybrid,deriv,hvcoord,tl,nets,nete )

    type(element_t)   , intent(inout) :: elem(:)
    type(hybrid_t)    , intent(in)    :: hybrid
    type(derivative_t), intent(in)    :: deriv
    type (hvcoord_t)  , intent(in)    :: hvcoord
    type (TimeLevel_t), intent(in)    :: tl
    integer,            intent(in)    :: nets,nete

    ! local variables
    integer :: ie,t,k
    real (kind=real_kind) :: gradtemp(np,np,2,nets:nete)
    real (kind=real_kind) :: temp(np,np,nlev),ps_ref(np,np)
    real (kind=real_kind) :: ptop_over_press


    ! other theta specific model initialization should go here    
    do ie=nets,nete
       gradtemp(:,:,:,ie) = gradient_sphere( elem(ie)%state%phis(:,:), deriv, elem(ie)%Dinv)
    enddo
    call make_C0_vector(gradtemp,elem,hybrid,nets,nete)
    
    do ie=nets,nete
      elem(ie)%derived%gradphis(:,:,:) = gradtemp(:,:,:,ie)
      ! compute w_i(nlevp)
      if (theta_hydrostatic_mode) then
        elem(ie)%state%w_i = 0.0D0
      else
        elem(ie)%state%w_i(:,:,nlevp,tl%n0) = (&
           elem(ie)%state%v(:,:,1,nlev,tl%n0)*elem(ie)%derived%gradphis(:,:,1) + &
           elem(ie)%state%v(:,:,2,nlev,tl%n0)*elem(ie)%derived%gradphis(:,:,2))/g
      endif

      ! assign phinh_i(nlevp) to be phis at all timelevels
      do t=1,timelevels
         elem(ie)%state%phinh_i(:,:,nlevp,t) = elem(ie)%state%phis(:,:)
      enddo

      ! initialize reference states used by hyberviscosity
#define HV_REFSTATES_V2
#ifdef HV_REFSTATES_V0
      elem(ie)%derived%dp_ref=0
      elem(ie)%derived%phi_ref=0
      elem(ie)%derived%theta_ref=0
#endif
#ifdef HV_REFSTATES_V1
      ps_ref(:,:) = hvcoord%ps0 * exp ( -elem(ie)%state%phis(:,:)/(Rgas*TREF)) 
      do k=1,nlev
         elem(ie)%derived%dp_ref(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
              (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:)
      enddo
      call set_theta_ref(hvcoord,elem(ie)%derived%dp_ref,elem(ie)%derived%theta_ref)
      temp=elem(ie)%derived%theta_ref*elem(ie)%derived%dp_ref
      call phi_from_eos(hvcoord,elem(ie)%state%phis,&
           temp,elem(ie)%derived%dp_ref,elem(ie)%derived%phi_ref)
      elem(ie)%derived%theta_ref=0
#endif
#ifdef HV_REFSTATES_V2
      ps_ref(:,:) = hvcoord%ps0 * exp ( -elem(ie)%state%phis(:,:)/(Rgas*TREF)) 
      do k=1,nlev
         elem(ie)%derived%dp_ref(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
              (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:)
      enddo
      call set_theta_ref(hvcoord,elem(ie)%derived%dp_ref,elem(ie)%derived%theta_ref)
      temp=elem(ie)%derived%theta_ref*elem(ie)%derived%dp_ref
      call phi_from_eos(hvcoord,elem(ie)%state%phis,&
           temp,elem(ie)%derived%dp_ref,elem(ie)%derived%phi_ref)
      elem(ie)%derived%theta_ref=0
      elem(ie)%derived%dp_ref=0
#endif
    enddo 


    ! unit test for analytic jacobian used by IMEX methods
#if 0
    if (.not. theta_hydrostatic_mode) &
         call test_imex_jacobian(elem,hybrid,hvcoord,tl,nets,nete)
#endif

    !$omp master
    ! 
    ! compute scaling of sponge layer damping 
    !
    nlev_tom=0
    if (hybrid%masterthread) write(iulog,*) "sponge layer nu_top viscosity scaling factor"
    do k=1,nlev
       !press = (hvcoord%hyam(k)+hvcoord%hybm(k))*hvcoord%ps0
       !ptop  = hvcoord%hyai(1)*hvcoord%ps0
       ! sponge layer starts at p=4*ptop 
       ! 
       ! some test cases have ptop=200mb
       if (hvcoord%etai(1)==0) then
          ! pure sigma coordinates could have etai(1)=0
          ptop_over_press = hvcoord%etam(1) / hvcoord%etam(k)  
       else
          ptop_over_press = hvcoord%etai(1) / hvcoord%etam(k)  
       endif

       ! active for p<10*ptop (following cd_core.F90 in CAM-FV)
       ! CAM 26L and 30L:  top 3 levels 
       ! E3SM 72L:  top 6 levels
       !original cam formula
       !nu_scale_top(k) = 8*(1+tanh(log(ptop_over_press))) ! active for p<4*ptop
       nu_scale_top(k) = 16*ptop_over_press**2 / (ptop_over_press**2 + 1)

       if (nu_scale_top(k)<0.15d0) nu_scale_top(k)=0

       !nu_scale_top(k) = 8*(1+.911*tanh(log(ptop_over_press))) ! active for p<6.5*ptop
       !if (nu_scale_top(k)<1d0) nu_scale_top(k)=0

       ! original CAM3/preqx formula
       !if (k==1) nu_scale_top(k)=4
       !if (k==2) nu_scale_top(k)=2
       !if (k==3) nu_scale_top(k)=1
       !if (k>3) nu_scale_top(k)=0

       if (nu_scale_top(k)>0) nlev_tom=k

       if (hybrid%masterthread) then
          if (nu_scale_top(k)>0) write(iulog,*) "  nu_scale_top ",k,nu_scale_top(k)
       end if
    end do
    if (hybrid%masterthread) then
       write(iulog,*) "  nlev_tom ",nlev_tom
    end if
    !$omp end master
    !$omp barrier

  end subroutine 

end module 
