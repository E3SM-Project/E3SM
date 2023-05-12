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
  use derivative_mod,     only: derivative_t,gradient_sphere, laplace_sphere_wk
  use hybvcoord_mod, 	  only: hvcoord_t
  use hybrid_mod,         only: hybrid_t
  use dimensions_mod,     only: np,nlev,nlevp
  use element_ops,        only: initialize_reference_states
  use element_state,      only: timelevels, nu_scale_top, nlev_tom
  use viscosity_mod,      only: make_c0_vector
  use kinds,              only: real_kind,iulog
  use control_mod,        only: qsplit,theta_hydrostatic_mode, hv_ref_profiles, &
       hv_theta_correction, tom_sponge_start
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
      call initialize_reference_states(hvcoord, elem(ie)%state%phis, &
                                       elem(ie)%derived%dp_ref,      &
                                       elem(ie)%derived%theta_ref,   &
                                       elem(ie)%derived%phi_ref)

      if (hv_theta_correction/=0) then
         ! compute weak laplace of p_ref
         ps_ref(:,:) = hvcoord%ps0 * exp ( -elem(ie)%state%phis(:,:)/(Rgas*TREF))
         do k=1,nlev
            ! compute lap(p) for z surface correction
            temp(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps_ref(:,:)
            elem(ie)%derived%lap_p_wk(:,:,k)=laplace_sphere_wk(temp(:,:,k),deriv,elem(ie),&
                 var_coef=.false.)
         enddo
      endif


    enddo


    ! unit test for analytic jacobian and tri-diag solve used by IMEX methods
    if (.not. theta_hydrostatic_mode) &
         call test_imex_jacobian(elem,hybrid,hvcoord,tl,nets,nete)

    !$omp master
    !
    ! compute scaling of sponge layer damping
    !
    nlev_tom=0
    if (hybrid%masterthread) write(iulog,*) "sponge layer nu_top viscosity scaling factor"
    do k=1,nlev
       if (tom_sponge_start==0) then
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
       else
          ptop_over_press = (tom_sponge_start/1d3) / hvcoord%etam(k)
          nu_scale_top(k)=0.15d0 * ptop_over_press**2
       endif

       if (nu_scale_top(k)<0.15d0) nu_scale_top(k)=0
       if (nu_scale_top(k)>8d0) nu_scale_top(k)=8d0


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
