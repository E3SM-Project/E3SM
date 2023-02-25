#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!
!  Man dynamics routines for "theta" nonhydrostatic model
!  Original version: Mark Taylor 2017/1
!  
!  2018/8 TOM sponge layer scaling from P. Lauritzen
!  09/2018: O. Guba  code for new ftypes
!  2018/12: M. Taylor apply forcing assuming nearly constant p 
!  2019/5:  M. Taylor time-split TOM dissipation and hyperviscsity
!  2019/7:  M. Taylor add dp3d limiter to prevent zero thickness layers
!
module prim_advance_mod

  use bndry_mod,          only: bndry_exchangev
  use control_mod,        only: dcmip16_mu, dcmip16_mu_s, hypervis_order, hypervis_subcycle,&
    integration, nu, nu_div, nu_p, nu_s, nu_top, prescribed_wind, qsplit, rsplit, test_case,&
    theta_hydrostatic_mode, tstep_type, theta_advect_form, hypervis_subcycle_tom, pgrad_correction,&
    vtheta_thresh, dp3d_thresh
  use derivative_mod,     only: derivative_t, divergence_sphere, gradient_sphere, laplace_sphere_wk,&
    laplace_z, vorticity_sphere, vlaplace_sphere_wk 
  use derivative_mod,     only: subcell_div_fluxes, subcell_dss_fluxes
  use dimensions_mod,     only: max_corner_elem, nlev, nlevp, np, qsize
  use edge_mod,           only: edge_g, edgevpack_nlyr, edgevunpack_nlyr
  use edgetype_mod,       only: EdgeBuffer_t,  EdgeDescriptor_t, edgedescriptor_t
  use element_mod,        only: element_t
  use element_state,      only: nu_scale_top, nlev_tom, max_itercnt, max_deltaerr,max_reserr
  use element_ops,        only: state0, get_R_star, tref_lapse_rate
  use eos,                only: pnh_and_exner_from_eos,pnh_and_exner_from_eos2,phi_from_eos
  use hybrid_mod,         only: hybrid_t
  use hybvcoord_mod,      only: hvcoord_t
  use kinds,              only: iulog, real_kind
  use perf_mod,           only: t_adj_detailf, t_barrierf, t_startf, t_stopf ! _EXTERNAL
  use parallel_mod,       only: abortmp, global_shared_buf, global_shared_sum, iam, parallel_t
  use physical_constants, only: Cp, cp, cpwater_vapor, g, kappa, Rgas, Rwater_vapor, p0, TREF
  use physics_mod,        only: virtual_specific_heat, virtual_temperature
  use prim_si_mod,        only: preq_vertadv_v1
  use reduction_mod,      only: parallelmax, reductionbuffer_ordered_1d_t
  use time_mod,           only: timelevel_qdp, timelevel_t
  use prim_state_mod,     only: prim_diag_scalars, prim_energy_halftimes
#if !defined(CAM) && !defined(SCREAM)
  use test_mod,           only: set_prescribed_wind
#endif
  use viscosity_theta,    only: biharmonic_wk_theta

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif
 
#ifdef HOMMEXX_BFB_TESTING
  use bfb_mod,        only: cxx_log
#endif

  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init1, advance_hypervis, &
            applycamforcing_dynamics, compute_andor_apply_rhs, limiter_dp3d_k

contains





  subroutine prim_advance_init1(par, elem,integration)
        
    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i
    integer :: ie


  end subroutine prim_advance_init1



  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)
    use imex_mod, only: compute_stage_value_dirk
#ifdef ARKODE
    use arkode_mod,     only: parameter_list, evolve_solution, &
                              calc_nonlinear_stats, update_nonlinear_stats, &
                              rel_tol, abs_tol
    use arkode_tables,  only: table_list, butcher_table_set, set_Butcher_tables
    use iso_c_binding
#endif

    type (element_t),      intent(inout), target :: elem(:)
    type (derivative_t),   intent(in)            :: deriv
    type (hvcoord_t)                             :: hvcoord
    type (hybrid_t),       intent(in)            :: hybrid
    real (kind=real_kind), intent(in)            :: dt
    type (TimeLevel_t)   , intent(in)            :: tl
    integer              , intent(in)            :: nets
    integer              , intent(in)            :: nete
    logical,               intent(in)            :: compute_diagnostics

    real (kind=real_kind) :: dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) :: itertol,a1,a2,a3,a4,a5,a6,ahat1,ahat2
    real (kind=real_kind) :: ahat3,ahat4,ahat5,ahat6,dhat1,dhat2,dhat3,dhat4
    real (kind=real_kind) ::  gamma,delta,ap,aphat,dhat5,offcenter

    integer :: ie,nm1,n0,np1,nstep,qsplit_stage,k
    integer :: n,i,j,maxiter

#ifdef ARKODE 
    type(parameter_list)    :: arkode_parameters
    type(table_list)        :: arkode_table_list
    type(butcher_table_set) :: arkode_table_set
    integer(C_INT)          :: ierr
#endif

    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep


! integration = "explicit"
!
!   tstep_type=1  RK2
!   tstep_type=4  Kinnmark&Gray RK 5 stage 2nd order            CFL=4.00
!   tstep_type=5  Kinnmark&Gray RK 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                 From Paul Ullrich.  3rd order for nonlinear terms also
!                 K&G method is only 3rd order for linear
!   tstep_type=7  KG5+BE      KG5(2nd order, 4.0CFL) + BE.  1st order max stability IMEX
!   tstep_type=8  KG3+BE/CN   KG3 2nd order explicit, 1st order off-centering implicit
!   tstep_type=9  KGU53+BE/CN KGU53 3rd order explicit, 2st order implicit
!   tstep_type=10 KGU42+BE/optimized, from O. Guba
!

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

!   this should not be needed, but in case physics update u without updating w b.c.:
    do ie=nets,nete
       elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
            elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
    enddo
 
#if !defined(CAM) && !defined(SCREAM)
    ! if "prescribed wind" set dynamics explicitly and skip time-integration
    if (prescribed_wind ==1 ) then
       call set_prescribed_wind(elem,deriv,hybrid,hvcoord,dt,tl,nets,nete,eta_ave_w)
       call t_stopf('prim_advance_exp')
       return
    endif
#endif

    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt
    if (tstep_type==1) then 
       ! RK2                                                                                                              
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))                                               
       call compute_andor_apply_rhs(np1,n0,n0,dt/2,elem,hvcoord,hybrid,&                                              
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)                                                      
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))                                                     
       call compute_andor_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&                                               
            deriv,nets,nete,.false.,eta_ave_w,1.d0,1.d0,1.d0)                                                             


    else if (tstep_type==4) then ! explicit table from IMEX-KG254  method                                                              
      call compute_andor_apply_rhs(np1,n0,n0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,dt/6,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,3*dt/8,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w*1d0,1.d0,1.d0,1.d0)



    else if (tstep_type==5) then
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call compute_andor_apply_rhs(nm1,n0,n0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,1.d0,1.d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_andor_apply_rhs(np1,n0,nm1,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_andor_apply_rhs(np1,n0,np1,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_andor_apply_rhs(np1,n0,np1,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%vtheta_dp(:,:,:,nm1)= (5*elem(ie)%state%vtheta_dp(:,:,:,nm1) &
               - elem(ie)%state%vtheta_dp(:,:,:,n0) )/4
          elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          elem(ie)%state%w_i(:,:,1:nlevp,nm1)= (5*elem(ie)%state%w_i(:,:,1:nlevp,nm1) &
                  - elem(ie)%state%w_i(:,:,1:nlevp,n0) )/4
          elem(ie)%state%phinh_i(:,:,1:nlev,nm1)= (5*elem(ie)%state%phinh_i(:,:,1:nlev,nm1) &
                  - elem(ie)%state%phinh_i(:,:,1:nlev,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_andor_apply_rhs(np1,nm1,np1,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,1.d0,1.d0)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
!=========================================================================================
    else if (tstep_type==7) then ! KG5(2nd order CFL=4) + BE  MAX STABILITY
      a1=0d0
      a2=1-a1
      dt2=dt/4
      call compute_andor_apply_rhs(np1,n0,n0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,a1*dt2,np1,a2*dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      dt2=dt/6
      call compute_andor_apply_rhs(nm1,n0,np1,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,a1*dt2,nm1,a2*dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      dt2=3*dt/8
      call compute_andor_apply_rhs(np1,n0,nm1,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,a1*dt2,np1,a2*dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)


      dt2=dt/2
      call compute_andor_apply_rhs(np1,n0,np1,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,a1*dt2,np1,a2*dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      call compute_andor_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w*1d0,1.d0,0d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,a1*dt,np1,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      !  u0 saved in elem(n0)
      !  u2 saved in elem(nm1)
      !  u4 saved in elem(np1)
      !  u5 = u0 + dt*N(u4) + dt*6/22*S(u0) + dt*6/22 S(u1) + dt*10/22* S(u5)
!===================================================================================
   elseif (tstep_type == 8 ) then ! KG3 + CN + offcentering

      ! introduce 1st order offcentering
      offcenter = 0.5d0
      aphat = 0.5d0-offcenter
      dhat3 = 0.5d0+offcenter

      call compute_andor_apply_rhs(np1,n0,n0,dt/2,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0) !   aphat/ap,1d0)               

      call compute_stage_value_dirk(nm1,0d0,n0,aphat*dt/2,np1,dhat3*dt/2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      call compute_andor_apply_rhs(np1,n0,np1,dt/2,elem,hvcoord,hybrid,&
        deriv,nets,nete,.false.,0d0,1d0,0d0,1d0)

      call compute_stage_value_dirk(nm1,0d0,n0,aphat*dt/2,np1,dhat3*dt/2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      ! introduce 1st order offcentering
      offcenter = 0.01d0
      aphat = 0.5d0-offcenter
      dhat3 = 0.5d0+offcenter

      call compute_andor_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
           deriv,nets,nete,.false.,eta_ave_w,1d0,0d0,1d0)
      call compute_stage_value_dirk(nm1,0d0,n0,aphat*dt,np1,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

    else if (tstep_type==9) then 
       ! KGU5-3 (3rd order) with IMEX backward euler (2nd order)
       ! 
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call compute_andor_apply_rhs(nm1,n0,n0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,0.d0,1.d0)
       call compute_stage_value_dirk(nm1,0d0,n0,0d0,nm1,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,maxiter,itertol)

       ! u2 = u0 + dt/5 RHS(u1)
       call compute_andor_apply_rhs(np1,n0,nm1,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       call compute_stage_value_dirk(nm1,0d0,n0,0d0,np1,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,maxiter,itertol)

       ! u3 = u0 + dt/3 RHS(u2)
       call compute_andor_apply_rhs(np1,n0,np1,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       call compute_stage_value_dirk(nm1,0d0,n0,0d0,np1,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,maxiter,itertol)

       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_andor_apply_rhs(np1,n0,np1,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
       call compute_stage_value_dirk(nm1,0d0,n0,0d0,np1,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,maxiter,itertol)


       ! u5 = u1 + dt 3/4 RHS(u4)
       call compute_andor_apply_rhs(np1,nm1,np1,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,0.d0,1.d0)
       ! u(np1) = [u1 + 3dt/4 RHS(u4)] +  1/4 (u1 - u0)    STABLE
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,np1)=elem(ie)%state%v(:,:,:,:,np1)+&
               (elem(ie)%state%v(:,:,:,:,nm1)-elem(ie)%state%v(:,:,:,:,n0))/4
          elem(ie)%state%vtheta_dp(:,:,:,np1)=elem(ie)%state%vtheta_dp(:,:,:,np1)+&
               (elem(ie)%state%vtheta_dp(:,:,:,nm1)-elem(ie)%state%vtheta_dp(:,:,:,n0))/4
          elem(ie)%state%dp3d(:,:,:,np1)=elem(ie)%state%dp3d(:,:,:,np1)+&
               (elem(ie)%state%dp3d(:,:,:,nm1)-elem(ie)%state%dp3d(:,:,:,n0))/4
          elem(ie)%state%w_i(:,:,1:nlevp,np1)=elem(ie)%state%w_i(:,:,1:nlevp,np1)+&
              (elem(ie)%state%w_i(:,:,1:nlevp,nm1)-elem(ie)%state%w_i(:,:,1:nlevp,n0))/4
          elem(ie)%state%phinh_i(:,:,1:nlev,np1)=elem(ie)%state%phinh_i(:,:,1:nlev,np1)+&
               (elem(ie)%state%phinh_i(:,:,1:nlev,nm1)-elem(ie)%state%phinh_i(:,:,1:nlev,n0))/4
          call limiter_dp3d_k(elem(ie)%state%dp3d(:,:,:,np1),elem(ie)%state%vtheta_dp(:,:,:,np1),&
               elem(ie)%spheremp,hvcoord%dp0)
       enddo

       !  n0          nm1       np1 
       ! u0*5/18  + u1*5/18  + u5*8/18
       a1=5*dt/18
       a2=dt/36    ! 5/18 - 1/4 (due to the 1/4*u1 added above)
       a3=8*dt/18
       call compute_stage_value_dirk(nm1,a2,n0,a1,np1,a3,elem,hvcoord,hybrid,&
            deriv,nets,nete,maxiter,itertol)

    else if (tstep_type==10) then ! KG5(2nd order CFL=4) + optimized
      dt2=dt/4
      call compute_andor_apply_rhs(nm1,n0,n0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,0d0,nm1,dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      dt2=dt/6
      call compute_andor_apply_rhs(np1,n0,nm1,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,0d0,np1,dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      dt2=3*dt/8
      call compute_andor_apply_rhs(np1,n0,np1,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,0d0,np1,dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)


      dt2=dt/2
      call compute_andor_apply_rhs(np1,n0,np1,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,0.d0,1.d0)
      call compute_stage_value_dirk(nm1,0d0,n0,0d0,np1,dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)

      call compute_andor_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w*1d0,1.d0,0d0,1.d0)


      a1=.24362d0   
      a2=.34184d0 
      a3=1-(a1+a2)
      call compute_stage_value_dirk(nm1,a2*dt,n0,a1*dt,np1,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      !  u0 saved in elem(n0)
      !  u1 saved in elem(nm1)
      !  u4 saved in elem(np1)



#ifdef ARKODE
    else if (tstep_type==20) then ! ARKode RK2
      call set_Butcher_tables(arkode_table_set, arkode_table_list%RK2)

    else if (tstep_type==21) then ! ARKode Kinnmark, Gray, Ullrich 3rd-order, 5-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%KGU35)

    else if (tstep_type==22) then ! ARKode Ascher 2nd/2nd/2nd-order, 3-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARS232)

    else if (tstep_type==23) then ! ARKode Candidate ARK453 Method
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARK453)

    else if (tstep_type==24) then ! ARKode Ascher 2nd/2nd/2nd-order, 3-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARS222)

    else if (tstep_type==25) then ! ARKode Ascher 3rd/4th/3rd-order, 3-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARS233)

    else if (tstep_type==26) then ! ARKode Ascher 3rd/3rd/3rd-order, 4-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARS343)

    else if (tstep_type==27) then ! ARKode Ascher 3rd/3rd/3rd-order, 5-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARS443)

    else if (tstep_type==28) then ! ARKode Kennedy 3rd/3rd/3rd-order, 4-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARK324)

    else if (tstep_type==29) then ! ARKode Kennedy 4th/4th/4th-order, 6-stage
      call set_Butcher_tables(arkode_table_set, arkode_table_list%ARK436)

    else if (tstep_type==30) then ! ARKode Conde et al ssp3(3,3,3)a (renamed here)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%SSP3333B)

    else if (tstep_type==31) then ! ARKode Conde et al ssp3(3,3,3)b (renamed here)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%SSP3333C)

    else if (tstep_type==32) then ! ARKode IMKG 2nd-order, 4 stage (2 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG232)

    else if (tstep_type==33) then ! ARKode IMKG 2nd-order, 5 stage (2 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG242)

    else if (tstep_type==34) then ! ARKode IMKG 2nd-order, 5 stage (3 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG243)

    else if (tstep_type==35) then ! ARKode IMKG 2nd-order, 6 stage (2 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG252)

    else if (tstep_type==36) then ! ARKode IMKG 2nd-order, 6 stage (3 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG253)

    else if (tstep_type==37) then ! ARKode IMKG 2nd-order, 6 stage (4 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG254)

    else if (tstep_type==38) then ! ARKode IMKG 3rd-order, 5 stage (2 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG342)

    else if (tstep_type==39) then ! ARKode IMKG 3rd-order, 5 stage (3 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG343)

    else if (tstep_type==40) then ! ARKode IMKG 3rd-order, 6 stage (3 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG353)

    else if (tstep_type==41) then ! ARKode IMKG 3rd-order, 6 stage (4 implicit)
      call set_Butcher_tables(arkode_table_set, arkode_table_list%IMKG354)

    else 
       call abortmp('ERROR: bad choice of tstep_type')
    endif

    ! Use ARKode to advance solution
    if (tstep_type >= 20) then

      ! If implicit solves are involved, set corresponding parameters
      if (arkode_table_set%imex /= 1) then
        ! Iteration tolerances (appear in WRMS array as rtol*|u_i| + atol_i)
        arkode_parameters%rtol = rel_tol
        if (abs_tol < 0.d0) then
          arkode_parameters%atol(1) = 1.d1*arkode_parameters%rtol ! assumes u ~ 1e1
          arkode_parameters%atol(2) = 1.d1*arkode_parameters%rtol ! assumes v ~ 1e1
          arkode_parameters%atol(3) = 1.d1*arkode_parameters%rtol ! assumes w_i ~ 1e1
          arkode_parameters%atol(4) = 1.d5*arkode_parameters%rtol ! assumes phinh_i ~ 1e5
          arkode_parameters%atol(5) = 1.d6*arkode_parameters%rtol ! assumes vtheta_dp ~ 1e6
          arkode_parameters%atol(6) = 1.d0*arkode_parameters%rtol ! assumes dp3d ~ 1e0
        else
          arkode_parameters%atol(:) = abs_tol
        end if
      end if

      ! use ARKode solver to evolve solution
      ierr = evolve_solution(elem, nets, nete, deriv, hvcoord, hybrid, &
                             dt, eta_ave_w, n0, np1,  arkode_parameters, &
                             arkode_table_set)
      if (ierr /= 0) then
        call abortmp('ARKode evolve failed')
      endif
      if (calc_nonlinear_stats) then
        call update_nonlinear_stats()
      end if
    end if
#else
    else 
       call abortmp('ERROR: bad choice of tstep_type')
    endif
#endif


    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    ! forward-in-time, hypervis applied to dp3d
    if (compute_diagnostics) then
       call t_startf("prim_diag")
       call prim_energy_halftimes(elem,hvcoord,tl,5,.false.,nets,nete)
       call prim_diag_scalars(elem,hvcoord,tl,5,.false.,nets,nete)
       call t_stopf("prim_diag")
    endif

    if (hypervis_order == 2 .and. nu>0) &
         call advance_hypervis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)


    ! warning: advance_physical_vis currently requires levels that are equally spaced in z
    if (dcmip16_mu>0) call advance_physical_vis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt,dcmip16_mu_s,dcmip16_mu)

    if (compute_diagnostics) then
       call t_startf("prim_diag")
       call prim_energy_halftimes(elem,hvcoord,tl,6,.false.,nets,nete)
       call prim_diag_scalars(elem,hvcoord,tl,6,.false.,nets,nete)
       call t_stopf("prim_diag")
    endif
    call t_stopf('prim_advance_exp')
  end subroutine prim_advance_exp

!----------------------------- APPLYCAMFORCING-DYNAMICS ----------------------------

  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt,nets,nete)

  type (element_t)     ,  intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete

  integer :: k,ie
  do ie=nets,nete

     elem(ie)%state%vtheta_dp(:,:,:,np1) = elem(ie)%state%vtheta_dp(:,:,:,np1) + dt*elem(ie)%derived%FVTheta(:,:,:)
     elem(ie)%state%phinh_i(:,:,1:nlev,np1) = elem(ie)%state%phinh_i(:,:,1:nlev,np1) + dt*elem(ie)%derived%FPHI(:,:,1:nlev)

     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt*elem(ie)%derived%FM(:,:,1:2,:)

#ifndef CAM
     elem(ie)%state%w_i(:,:,1:nlev,np1) = elem(ie)%state%w_i(:,:,1:nlev,np1) + dt*elem(ie)%derived%FM(:,:,3,:)
#endif

     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,np1) = (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo
  
  end subroutine applyCAMforcing_dynamics


!----------------------------- ADVANCE-HYPERVIS ----------------------------

  subroutine advance_hypervis(elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  integer :: k2,k,kptr,i,j,ie,ic,nt,nlyr_tot,nlyr_tom,ssize
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,4,nets:nete)      :: stens  ! dp3d,theta,w,phi


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np,4) :: lap_s  ! dp3d,theta,w,phi
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: exner0(nlev)
  real (kind=real_kind) :: heating(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlevp)    
  real (kind=real_kind) :: temp(np,np,nlev)    
  real (kind=real_kind) :: temp_i(np,np,nlevp)    
  real (kind=real_kind) :: dt,xfac

  integer :: l1p,l2p,l1n,l2n,l
  call t_startf('advance_hypervis')

#ifdef HOMMEXX_BFB_TESTING
  ! Exchange all vars even in hydro mode, for the sake of bfb comparison with xx code
  nlyr_tot=6*nlev  ! total amount of data for DSS
  nlyr_tom=6*nlev_tom
  ssize=4*nlev
#else
  if (theta_hydrostatic_mode) then
     nlyr_tot=4*nlev        ! dont bother to dss w_i and phinh_i
     nlyr_tom=4*nlev_tom
     ssize=2*nlev
  else
     nlyr_tot=6*nlev  ! total amount of data for DSS
     nlyr_tom=6*nlev_tom
     ssize=4*nlev
  endif
#endif
  
  do k=1,nlev
     exner0(k) = (hvcoord%etam(k)*hvcoord%ps0/p0 )**kappa
  enddo


  do ie=nets,nete
     ! convert vtheta_dp -> theta
     do k=1,nlev
        elem(ie)%state%vtheta_dp(:,:,k,nt)=&
             elem(ie)%state%vtheta_dp(:,:,k,nt)/elem(ie)%state%dp3d(:,:,k,nt)
     enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt=dt2/hypervis_subcycle
  do ic=1,hypervis_subcycle
     do ie=nets,nete
        ! remove ref state
        elem(ie)%state%vtheta_dp(:,:,:,nt)=elem(ie)%state%vtheta_dp(:,:,:,nt)-&
             elem(ie)%derived%theta_ref(:,:,:)
        elem(ie)%state%phinh_i(:,:,:,nt)=elem(ie)%state%phinh_i(:,:,:,nt)-&
             elem(ie)%derived%phi_ref(:,:,:)
        elem(ie)%state%dp3d(:,:,:,nt)=elem(ie)%state%dp3d(:,:,:,nt)-&
             elem(ie)%derived%dp_ref(:,:,:)
     enddo
     
     call biharmonic_wk_theta(elem,stens,vtens,deriv,edge_g,hybrid,nt,nets,nete)
     
     do ie=nets,nete
        !add ref state back
        elem(ie)%state%vtheta_dp(:,:,:,nt)=elem(ie)%state%vtheta_dp(:,:,:,nt)+&
             elem(ie)%derived%theta_ref(:,:,:)
        elem(ie)%state%phinh_i(:,:,:,nt)=elem(ie)%state%phinh_i(:,:,:,nt)+&
             elem(ie)%derived%phi_ref(:,:,:)
        elem(ie)%state%dp3d(:,:,:,nt)=elem(ie)%state%dp3d(:,:,:,nt)+&
             elem(ie)%derived%dp_ref(:,:,:)
        
        
        ! comptue mean flux
        if (nu_p>0) then
           elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)/hypervis_subcycle
           elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                eta_ave_w*stens(:,:,:,1,ie)/hypervis_subcycle
        endif
        do k=1,nlev
           vtens(:,:,:,k,ie)=-nu  *vtens(:,:,:,k,ie) ! u,v
           stens(:,:,k,1,ie)=-nu_p*stens(:,:,k,1,ie) ! dp3d
           stens(:,:,k,2,ie)=-nu  *stens(:,:,k,2,ie) ! theta
           stens(:,:,k,3,ie)=-nu  *stens(:,:,k,3,ie) ! w
           stens(:,:,k,4,ie)=-nu_s*stens(:,:,k,4,ie) ! phi
        enddo
        
        kptr=0;      call edgeVpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,nlyr_tot)
        kptr=2*nlev; call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),ssize,kptr,nlyr_tot)
        
     enddo

     call t_startf('ahdp_bexchV2')
     call bndry_exchangeV(hybrid,edge_g)
     call t_stopf('ahdp_bexchV2')
     
     do ie=nets,nete
        
        kptr=0
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,nlyr_tot)
        kptr=2*nlev
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),ssize,kptr,nlyr_tot)
        
        
        ! apply inverse mass matrix, accumulate tendencies
        do k=1,nlev
#ifdef HOMMEXX_BFB_TESTING
           vtens(:,:,1,k,ie)=vtens(:,:,1,k,ie)*(elem(ie)%rspheremp(:,:)*dt)  ! u
           vtens(:,:,2,k,ie)=vtens(:,:,2,k,ie)*(elem(ie)%rspheremp(:,:)*dt)  ! v
           stens(:,:,k,1,ie)=stens(:,:,k,1,ie)*(elem(ie)%rspheremp(:,:)*dt)  ! dp3d
           stens(:,:,k,2,ie)=stens(:,:,k,2,ie)*(elem(ie)%rspheremp(:,:)*dt)  ! theta
           stens(:,:,k,3,ie)=stens(:,:,k,3,ie)*(elem(ie)%rspheremp(:,:)*dt)  ! w
           stens(:,:,k,4,ie)=stens(:,:,k,4,ie)*(elem(ie)%rspheremp(:,:)*dt)  ! phi
#else
           vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)  ! u
           vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)  ! v
           stens(:,:,k,1,ie)=dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)  ! dp3d
           stens(:,:,k,2,ie)=dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)  ! theta
           stens(:,:,k,3,ie)=dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)  ! w
           stens(:,:,k,4,ie)=dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)  ! phi
#endif
        enddo
        
        
        ! apply heating after updating state.  using updated v gives better results in PREQX model
        !
        ! d(IE)/dt =  cp*exner*d(Theta)/dt + phi d(dp3d)/dt   (Theta = dp3d*theta)
        !   Our eqation:  d(theta)/dt = diss(theta) - heating
        !   Assuming no diffusion on dp3d, we can approximate by:
        !   d(IE)/dt = exner*cp*dp3d * diss(theta)  - exner*cp*dp3d*heating               
        !
        ! KE dissipaiton will be given by:
        !   d(KE)/dt = dp3d*U dot diss(U)
        ! we want exner*cp*dp3d*heating = dp3d*U dot diss(U)
        ! and thus heating =  U dot diss(U) / exner*cp
        ! 
        ! compute exner needed for heating term and IE scaling
        ! this is using a mixture of data before viscosity and after viscosity 
#if 0
        temp(:,:,:)=elem(ie)%state%vtheta_dp(:,:,:,nt)*elem(ie)%state%dp3d(:,:,:,nt)
        call pnh_and_exner_from_eos(hvcoord,temp,&
             elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%phinh_i(:,:,:,nt),&
             pnh,exner,temp_i,caller='advance_hypervis')
        
        do k=1,nlev
           k2=min(k+1,nlev)
           if (theta_hydrostatic_mode) then
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie) ) / &
                   (exner(:,:,k)*Cp)

           else
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie)  +&
                   (elem(ie)%state%w_i(:,:,k,nt)*stens(:,:,k,3,ie)  +&
                     elem(ie)%state%w_i(:,:,k2,nt)*stens(:,:,k2,3,ie))/2 ) /  &
                   (exner(:,:,k)*Cp)  
           endif
           !elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
           !     +stens(:,:,k,2,ie)*hvcoord%dp0(k)*exner0(k)/(exner(:,:,k)*elem(ie)%state%dp3d(:,:,k,nt)&
           !     )  -heating(:,:,k)
           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
                  -heating(:,:,k)
        enddo
#endif

        do k=1,nlev
           elem(ie)%state%v(:,:,:,k,nt)=elem(ie)%state%v(:,:,:,k,nt) + &
                vtens(:,:,:,k,ie)
           elem(ie)%state%w_i(:,:,k,nt)=elem(ie)%state%w_i(:,:,k,nt) &
                +stens(:,:,k,3,ie)
           
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
                +stens(:,:,k,1,ie)
           
           elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt) &
                +stens(:,:,k,4,ie)

           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
                +stens(:,:,k,2,ie)
        enddo




     enddo ! ie
  enddo  ! subcycle

! convert vtheta_dp -> theta
  do ie=nets,nete            
     elem(ie)%state%vtheta_dp(:,:,:,nt)=&
          elem(ie)%state%vtheta_dp(:,:,:,nt)*elem(ie)%state%dp3d(:,:,:,nt)
    
     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,nt) = (elem(ie)%state%v(:,:,1,nlev,nt)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,nt)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  sponge layer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (nu_top>0 .and. hypervis_subcycle_tom>0 ) then
  dt=dt2/hypervis_subcycle_tom
  do ic=1,hypervis_subcycle_tom

     do ie=nets,nete
        do k=1,nlev_tom
           ! add regular diffusion near top
           lap_s(:,:,1)=laplace_sphere_wk(elem(ie)%state%dp3d     (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
           lap_s(:,:,2)=laplace_sphere_wk(elem(ie)%state%vtheta_dp(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
           lap_s(:,:,3)=laplace_sphere_wk(elem(ie)%state%w_i      (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
           lap_s(:,:,4)=laplace_sphere_wk(elem(ie)%state%phinh_i  (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
           lap_v=vlaplace_sphere_wk(elem(ie)%state%v            (:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
           
           xfac=dt*nu_scale_top(k)*nu_top

           vtens(:,:,:,k,ie)=xfac*lap_v(:,:,:)
           stens(:,:,k,1,ie)=xfac*lap_s(:,:,1)  ! dp3d
           stens(:,:,k,2,ie)=xfac*lap_s(:,:,2)  ! vtheta_dp
           stens(:,:,k,3,ie)=xfac*lap_s(:,:,3)  ! w_i
           stens(:,:,k,4,ie)=xfac*lap_s(:,:,4)  ! phi_i
        enddo
        
        kptr=0;      
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev_tom,kptr,nlyr_tom)
        kptr=2*nlev_tom; 
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,1,ie),nlev_tom,kptr,nlyr_tom)
        kptr=kptr+nlev_tom
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,2,ie),nlev_tom,kptr,nlyr_tom)
        if (.not.theta_hydrostatic_mode) then
           kptr=kptr+nlev_tom
           call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,3,ie),nlev_tom,kptr,nlyr_tom)
           kptr=kptr+nlev_tom
           call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,4,ie),nlev_tom,kptr,nlyr_tom)
        endif
     enddo
     
     call t_startf('ahdp_bexchV2')
     call bndry_exchangeV(hybrid,edge_g)
     call t_stopf('ahdp_bexchV2')
     
     do ie=nets,nete
        
        kptr=0
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev_tom,kptr,nlyr_tom)
        kptr=2*nlev_tom
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,1,ie),nlev_tom,kptr,nlyr_tom)
        kptr=kptr+nlev_tom
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,2,ie),nlev_tom,kptr,nlyr_tom)
        if (.not.theta_hydrostatic_mode) then
           kptr=kptr+nlev_tom
           call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,3,ie),nlev_tom,kptr,nlyr_tom)
           kptr=kptr+nlev_tom
           call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,4,ie),nlev_tom,kptr,nlyr_tom)
        endif
        
        
        ! apply inverse mass matrix, add tendency
        do k=1,nlev_tom
           elem(ie)%state%v(:,:,1,k,nt)=elem(ie)%state%v(:,:,1,k,nt) + &
                vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
           elem(ie)%state%v(:,:,2,k,nt)=elem(ie)%state%v(:,:,2,k,nt) + &
                vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
           elem(ie)%state%w_i(:,:,k,nt)=elem(ie)%state%w_i(:,:,k,nt) &
                +stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)
           
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
                +stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)
           elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt) &
                +stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)

           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
                +stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)
        enddo
     enddo ! ie
  enddo  ! subcycle

  endif



  call t_stopf('advance_hypervis')

  end subroutine advance_hypervis





  subroutine advance_physical_vis(elem,hvcoord,hybrid,deriv,nt,nets,nete,dt,mu_s,mu)
  !
  !  take one timestep of of physical viscosity (single laplace operator) for
  !  all state variables in both horizontal and vertical
  !  
  !  as of 2017/5, used only for the supercell test case
  !  so for now:
  !     dont bother to optimize
  !     apply only to perturbation from background state (supercell initial condition)
  !     uniform spacing in z with delz = 20km/nlev
  !
  !

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt, mu_s, mu
  integer :: nt,nets,nete

  ! local
  integer :: k,kptr,ie
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,4,nets:nete)      :: stens  ! dp3d,theta,w,phi
  real (kind=real_kind), dimension(np,np,nlevp,2,nets:nete)     :: stens_i 


  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: delz,delz_i

  real (kind=real_kind) :: theta_ref(np,np,nlev)

  real (kind=real_kind) :: theta_prime(np,np,nlev)
  real (kind=real_kind) :: phi_prime(np,np,nlevp)
  real (kind=real_kind) :: dp_prime(np,np,nlev)
  real (kind=real_kind) :: w_prime(np,np,nlevp)
  real (kind=real_kind) :: u_prime(np,np,2,nlev)

  !if(test_case .ne. 'dcmip2016_test3') call abortmp("dcmip16_mu is currently limited to dcmip16 test 3")

  call t_startf('advance_physical_vis')
  delz = 20d3/nlev
  delz_i = 20d3/nlevp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute reference states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete

     theta_ref(:,:,:) = state0(ie)%vtheta_dp(:,:,:,1)/state0(ie)%dp3d(:,:,:,1)
     elem(ie)%state%vtheta_dp(:,:,:,nt)=&
             elem(ie)%state%vtheta_dp(:,:,:,nt)/elem(ie)%state%dp3d(:,:,:,nt)

     ! perturbation variables
     u_prime(:,:,:,:)  = elem(ie)%state%v(:,:,:,:,nt)          -state0(ie)%v(:,:,:,:,1)
     w_prime(:,:,:)    = elem(ie)%state%w_i(:,:,:,nt)     -state0(ie)%w_i(:,:,:,1)
     dp_prime(:,:,:)   = elem(ie)%state%dp3d(:,:,:,nt)         -state0(ie)%dp3d(:,:,:,1)
     phi_prime(:,:,:)  = elem(ie)%state%phinh_i(:,:,:,nt) -state0(ie)%phinh_i(:,:,:,1)
     theta_prime(:,:,:)= elem(ie)%state%vtheta_dp(:,:,:,nt)-theta_ref(:,:,:)

     ! vertical viscosity
     call laplace_z(u_prime,    vtens(:,:,:,:,ie),2,nlev,delz)
     call laplace_z(dp_prime,   stens(:,:,:,1,ie),1,nlev,delz)
     call laplace_z(theta_prime,stens(:,:,:,2,ie),1,nlev,delz)
     call laplace_z(w_prime,    stens_i(:,:,:,1,ie),1,nlevp,delz_i)
     call laplace_z(phi_prime,  stens_i(:,:,:,2,ie),1,nlevp,delz_i)

     ! add in horizontal viscosity
     ! multiply by mass matrix for DSS
     ! horiz viscosity already has mass matrix built in
     ! for interface quantities, only use 1:nlev (dont apply at surface)
     do k=1,nlev
        lap_v = vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)

        vtens(:,:,1,k,ie) = (vtens(:,:,1,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,1))
        vtens(:,:,2,k,ie) = (vtens(:,:,2,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,2))

        stens(:,:,k,1,ie) = (stens(:,:,k,1,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,2,ie) = (stens(:,:,k,2,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%vtheta_dp(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,3,ie) = (stens_i(:,:,k,1,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%w_i(:,:,k,nt),deriv,elem(ie),var_coef=.false.) )

        stens(:,:,k,4,ie) = (stens_i(:,:,k,2,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%phinh_i(:,:,k,nt),deriv,elem(ie),var_coef=.false.) ) 

     enddo

     kptr=0
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,6*nlev)
     kptr=2*nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),4*nlev,kptr,6*nlev)
     
  enddo

  call bndry_exchangeV(hybrid,edge_g)
  
  do ie=nets,nete
     
     kptr=0
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,6*nlev)
     kptr=2*nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),4*nlev,kptr,6*nlev)
     
     ! apply inverse mass matrix, accumulate tendencies
     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,nt)=elem(ie)%state%v(:,:,1,k,nt) + &
             mu*dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%v(:,:,2,k,nt)=elem(ie)%state%v(:,:,2,k,nt) + &
             mu*dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%w_i(:,:,k,nt)=elem(ie)%state%w_i(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)
        
     enddo
  enddo


  ! convert vtheta_dp -> theta
  do ie=nets,nete            
     do k=1,nlev
        elem(ie)%state%vtheta_dp(:,:,k,nt)=&
             elem(ie)%state%vtheta_dp(:,:,k,nt)*elem(ie)%state%dp3d(:,:,k,nt)
     enddo

     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,nt) = (elem(ie)%state%v(:,:,1,nlev,nt)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,nt)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo


  call t_stopf('advance_physical_vis')

  end subroutine advance_physical_vis





!============================ stiff and or non-stiff ============================================

 subroutine compute_andor_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,scale1,scale2,scale3)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !   u(np1) = scale3*u(nm1) + dt2*DSS[ nonstiffRHS(u(n0))*scale1 + stiffRHS(un0)*scale2 ]
  !
  ! This subroutine was orgininally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  ! ===================================

  integer,              intent(in) :: np1,nm1,n0,nets,nete
  real*8,               intent(in) :: dt2
  logical,              intent(in) :: compute_diagnostics
  type (hvcoord_t),     intent(in) :: hvcoord
  type (hybrid_t),      intent(in) :: hybrid
  type (element_t),     intent(inout), target :: elem(:)
  type (derivative_t),  intent(in) :: deriv

  real (kind=real_kind) :: eta_ave_w,scale1,scale2,scale3  ! weighting for eta_dot_dpdn mean flux, scale of unm1

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:) :: phi_i
  real (kind=real_kind), pointer, dimension(:,:,:) :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:) :: vtheta_dp
   
  real (kind=real_kind) :: vtheta(np,np,nlev)
  real (kind=real_kind) :: vtheta_i(np,np,nlevp)
  real (kind=real_kind) :: omega_i(np,np,nlevp)
  real (kind=real_kind) :: omega(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)           ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)     
  real (kind=real_kind) :: phi(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)               ! nh (nonydro) pressure
  real (kind=real_kind) :: dp3d_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)         ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)    ! dpnh / dp3d at interfaces
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces
  real (kind=real_kind) :: KE(np,np,nlev)             ! Kinetic energy
  real (kind=real_kind) :: gradexner(np,np,2,nlev)    ! grad(p^kappa)
  real (kind=real_kind) :: gradphinh_i(np,np,2,nlevp) ! gradphi at interfaces
  real (kind=real_kind) :: mgrad(np,np,2,nlev)        ! gradphi metric term at cell centers
  real (kind=real_kind) :: gradKE(np,np,2,nlev)       ! grad(0.5 u^T u )
  real (kind=real_kind) :: wvor(np,np,2,nlev)         ! w vorticity term

  real (kind=real_kind) :: gradw_i(np,np,2,nlevp)
  real (kind=real_kind) :: v_gradw_i(np,np,nlevp)     
  real (kind=real_kind) :: v_gradtheta(np,np,nlev)     
  real (kind=real_kind) :: v_theta(np,np,2,nlev)
  real (kind=real_kind) :: div_v_theta(np,np,nlev)
  real (kind=real_kind) :: v_gradphinh_i(np,np,nlevp) ! v*gradphi at interfaces
  real (kind=real_kind) :: v_i(np,np,2,nlevp)

  real (kind=real_kind) :: v_vadv(np,np,2,nlev)     ! velocity vertical advection
  real (kind=real_kind) :: theta_vadv(np,np,nlev)   ! w,phi, theta  vertical advection term
  real (kind=real_kind) :: w_vadv_i(np,np,nlevp)      ! w,phi, theta  vertical advection term
  real (kind=real_kind) :: phi_vadv_i(np,np,nlevp)    ! w,phi, theta  vertical advection term

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: stens(np,np,nlev,3) ! tendencies w,phi,theta
                                               ! w,phi tendencies not computed at nlevp
  real (kind=real_kind) :: w_tens(np,np,nlevp)  ! need to update w at surface as well
  real (kind=real_kind) :: theta_tens(np,np,nlev)
  real (kind=real_kind) :: phi_tens(np,np,nlevp)
                                               

  real (kind=real_kind) :: pi(np,np,nlev)                ! hydrostatic pressure
  real (kind=real_kind) :: pi_i(np,np,nlevp)             ! hydrostatic pressure interfaces
  real (kind=real_kind), dimension(np,np,nlev) :: vgrad_p

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind) ::  vtemp(np,np,2,nlev)       ! generic gradient storage
  real (kind=real_kind), dimension(np,np) :: sdot_sum ! temporary field
  real (kind=real_kind) ::  v1,v2,w,d_eta_dot_dpdn_dn, T0
  integer :: i,j,k,kptr,ie, nlyr_tot

  call t_startf('compute_andor_apply_rhs')

  if (theta_hydrostatic_mode) then
     nlyr_tot=4*nlev        ! dont bother to dss w_i and phinh_i
  else
     nlyr_tot=5*nlev+nlevp  ! total amount of data for DSS
  endif
     
  do ie=nets,nete
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
     vtheta(:,:,:) = vtheta_dp(:,:,:)/dp3d(:,:,:)
     phi_i => elem(ie)%state%phinh_i(:,:,:,n0)

#ifdef ENERGY_DIAGNOSTICS
     if (.not. theta_hydrostatic_mode) then
        ! check w b.c.
        temp(:,:,1) =  (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
        do j=1,np
        do i=1,np
           if ( abs(temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,n0)) >1e-10) then
              write(iulog,*) 'WARNING:CAAR w(n0) does not satisfy b.c.',ie,i,j,k
              write(iulog,*) 'val1 = ',temp(i,j,1)
              write(iulog,*) 'val2 = ',elem(ie)%state%w_i(i,j,nlevp,n0)
              write(iulog,*) 'diff: ',temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,n0)
           endif
        enddo
        enddo
        ! w boundary condition. just in case:
        !elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
        !     elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g

        ! check for layer spacing <= 1m
        do k=1,nlev
        do j=1,np
        do i=1,np
           if ((phi_i(i,j,k)-phi_i(i,j,k+1)) < g) then
              write(iulog,*) 'WARNING:CAAR before ADV, delta z < 1m. ie,i,j,k=',ie,i,j,k
              write(iulog,*) 'phi(i,j,k)=  ',phi_i(i,j,k)
              write(iulog,*) 'phi(i,j,k+1)=',phi_i(i,j,k+1)
           endif
        enddo
        enddo
        enddo
     endif
#endif
     ! this routine will set dpnh_dp_i(nlevp)=1 - a very good approximation, that will
     ! then be corrected below, after the DSS.  
     call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_i,pnh,exner,dpnh_dp_i,caller='CAAR')

     dp3d_i(:,:,1) = dp3d(:,:,1)
     dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
     do k=2,nlev
        dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
     end do

     ! special averaging for velocity for energy conservation
     v_i(:,:,1:2,1) = elem(ie)%state%v(:,:,1:2,1,n0)  
     v_i(:,:,1:2,nlevp) = elem(ie)%state%v(:,:,1:2,nlev,n0)
     do k=2,nlev
        v_i(:,:,1,k) = (dp3d(:,:,k)*elem(ie)%state%v(:,:,1,k,n0) + &
             dp3d(:,:,k-1)*elem(ie)%state%v(:,:,1,k-1,n0) ) / (2*dp3d_i(:,:,k))
        v_i(:,:,2,k) = (dp3d(:,:,k)*elem(ie)%state%v(:,:,2,k,n0) + &
             dp3d(:,:,k-1)*elem(ie)%state%v(:,:,2,k-1,n0) ) / (2*dp3d_i(:,:,k))
     end do
     
     if (theta_hydrostatic_mode) then
        do k=nlev,1,-1          ! traditional Hydrostatic integral
           phi_i(:,:,k)=phi_i(:,:,k+1)+&
                Rgas*vtheta_dp(:,:,k)*exner(:,:,k)/pnh(:,:,k)
        enddo
        ! in H mode, ignore w contibutions to KE term
        ! set to zero so H and NH can share code and reduce if statements
        elem(ie)%state%w_i(:,:,:,n0)=0   
     endif

     do k=1,nlev
        phi(:,:,k) = (phi_i(:,:,k)+phi_i(:,:,k+1))/2  ! for diagnostics

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        vtemp(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*dp3d(:,:,k)
        vtemp(:,:,2,k) = elem(ie)%state%v(:,:,2,k,n0)*dp3d(:,:,k)
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:,k)

        divdp(:,:,k)=divergence_sphere(vtemp(:,:,:,k),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))
     enddo

     ! Compute omega =  Dpi/Dt   Used only as a DIAGNOSTIC
     pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
     omega_i(:,:,1)=0
     do k=1,nlev
        pi_i(:,:,k+1)=pi_i(:,:,k) + dp3d(:,:,k)
        omega_i(:,:,k+1)=omega_i(:,:,k)+divdp(:,:,k)
     enddo
     do k=1,nlev
#ifdef HOMMEXX_BFB_TESTING
        pi(:,:,k)=(pi_i(:,:,k) + pi_i(:,:,k+1))/2
#else
        pi(:,:,k)=pi_i(:,:,k) + dp3d(:,:,k)/2
#endif
        vtemp(:,:,:,k) = gradient_sphere( pi(:,:,k), deriv, elem(ie)%Dinv);
        vgrad_p(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1,k)+&
             elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2,k)
        omega(:,:,k) = (vgrad_p(:,:,k) - ( omega_i(:,:,k)+omega_i(:,:,k+1))/2) 
     enddo        

     ! ==================================================
     ! Compute eta_dot_dpdn
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0
        w_vadv_i=0
        phi_vadv_i=0
        theta_vadv=0
        v_vadv=0
     else
        sdot_sum=0
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k)
           eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
        end do


        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
        do k=1,nlev-1
           eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
        end do

        eta_dot_dpdn(:,:,1     ) = 0
        eta_dot_dpdn(:,:,nlevp)  = 0
        vtheta_i(:,:,1) =       0
        vtheta_i(:,:,nlevp) =   0


        ! ===========================================================
        ! Compute vertical advection of v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv_v1(elem(ie)%state%v(:,:,:,:,n0),eta_dot_dpdn,dp3d,v_vadv)

        ! compute (cp*theta) at interfaces
        ! for energy conservation, use averaging consistent with EOS
        ! dont bother to compute at surface and top since it will be multiplied by eta-dot
#if 0           
           do k=2,nlev  ! simple averaging
              vtheta_i(:,:,k) = (vtheta(:,:,k)+vtheta(:,:,k-1))/2
           enddo
#else
           ! E conserving average, but much more dissipative
           do k=2,nlev
#ifdef HOMMEXX_BFB_TESTING
              vtheta_i(:,:,k) = -dpnh_dp_i(:,:,k)* ((phi(:,:,k)-phi(:,:,k-1))/&
                                                    (exner(:,:,k)-exner(:,:,k-1)) / Cp)
#else
              vtheta_i(:,:,k) = -dpnh_dp_i(:,:,k)*(phi(:,:,k)-phi(:,:,k-1))/&
                   (exner(:,:,k)-exner(:,:,k-1)) / Cp
#endif
           enddo
#endif           



        do k=1,nlev
           ! average interface quantity to midpoints:
           temp(:,:,k) = (( eta_dot_dpdn(:,:,k)+eta_dot_dpdn(:,:,k+1))/2)*&
                (elem(ie)%state%w_i(:,:,k+1,n0)-elem(ie)%state%w_i(:,:,k,n0))
           
           ! theta vadv term at midoints
           theta_vadv(:,:,k)= eta_dot_dpdn(:,:,k+1)*vtheta_i(:,:,k+1) - &
                eta_dot_dpdn(:,:,k)*vtheta_i(:,:,k)
        enddo
        ! compute ave( ave(etadot) d/dx )
        do k=2,nlev
           w_vadv_i(:,:,k)  =(temp(:,:,k-1)+temp(:,:,k))/2
           phi_vadv_i(:,:,k)=eta_dot_dpdn(:,:,k)*(phi(:,:,k)-phi(:,:,k-1))
        end do
        w_vadv_i(:,:,1) = temp(:,:,1)
        w_vadv_i(:,:,nlevp) = temp(:,:,nlev)
        phi_vadv_i(:,:,1) = 0
        phi_vadv_i(:,:,nlevp) = 0

        ! final form of SB81 vertical advection operator:
        w_vadv_i=w_vadv_i/dp3d_i
        phi_vadv_i=phi_vadv_i/dp3d_i
     endif


     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
     do k=1,nlev  !  Loop index added (AAM)
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,k) + eta_ave_w*eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega_p(:,:,k) = &
             elem(ie)%derived%omega_p(:,:,k) + eta_ave_w*omega(:,:,k)
     enddo
     elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)

     ! ================================================
     ! w,phi tendencies including surface
     ! ================================================  
     do k=1,nlev
        ! compute gradphi at interfaces and then average to levels
        gradphinh_i(:,:,:,k)   = gradient_sphere(phi_i(:,:,k),deriv,elem(ie)%Dinv)   
           
        gradw_i(:,:,:,k)   = gradient_sphere(elem(ie)%state%w_i(:,:,k,n0),deriv,elem(ie)%Dinv)
        v_gradw_i(:,:,k) = v_i(:,:,1,k)*gradw_i(:,:,1,k) + v_i(:,:,2,k)*gradw_i(:,:,2,k)
        ! w - tendency on interfaces 
        w_tens(:,:,k) = (-w_vadv_i(:,:,k) - v_gradw_i(:,:,k))*scale1 - scale2*g*(1-dpnh_dp_i(:,:,k) )

        ! phi - tendency on interfaces
        ! vtemp(:,:,:,k) = gradphinh_i(:,:,:,k) + &
        !    (scale2-1)*hvcoord%hybi(k)*elem(ie)%derived%gradphis(:,:,:)
        v_gradphinh_i(:,:,k) = v_i(:,:,1,k)*gradphinh_i(:,:,1,k) &
             +v_i(:,:,2,k)*gradphinh_i(:,:,2,k) 
        phi_tens(:,:,k) =  (-phi_vadv_i(:,:,k) - v_gradphinh_i(:,:,k))*scale1 &
          + scale2*g*elem(ie)%state%w_i(:,:,k,n0)
        if (scale1/=scale2) then
           ! add imex phi_h splitting 
           ! use approximate phi_h = hybi*phis 
           ! could also use true hydrostatic pressure, but this requires extra DSS in dirk()
           phi_tens(:,:,k) =  phi_tens(:,:,k)+(scale1-scale2)*(&
                v_i(:,:,1,k)*elem(ie)%derived%gradphis(:,:,1) + &
                v_i(:,:,2,k)*elem(ie)%derived%gradphis(:,:,2) )*hvcoord%hybi(k)
        endif
     end do


     ! k =nlevp case, all terms in the imex methods are treated explicitly at the boundary
     k =nlevp 
    ! compute gradphi at interfaces and then average to levels
    gradphinh_i(:,:,:,k)   = gradient_sphere(phi_i(:,:,k),deriv,elem(ie)%Dinv)
    gradw_i(:,:,:,k)   = gradient_sphere(elem(ie)%state%w_i(:,:,k,n0),deriv,elem(ie)%Dinv)
    v_gradw_i(:,:,k) = v_i(:,:,1,k)*gradw_i(:,:,1,k) + v_i(:,:,2,k)*gradw_i(:,:,2,k)
    ! w - tendency on interfaces
    w_tens(:,:,k) = (-w_vadv_i(:,:,k) - v_gradw_i(:,:,k))*scale1 - scale1*g*(1-dpnh_dp_i(:,:,k) )

    ! phi - tendency on interfaces
    v_gradphinh_i(:,:,k) = v_i(:,:,1,k)*gradphinh_i(:,:,1,k) &
     +v_i(:,:,2,k)*gradphinh_i(:,:,2,k)
    phi_tens(:,:,k) =  (-phi_vadv_i(:,:,k) - v_gradphinh_i(:,:,k))*scale1 &
    + scale1*g*elem(ie)%state%w_i(:,:,k,n0)
    




     ! ================================================                                                                 
     ! v1,v2 tendencies:                                                                                          
     ! ================================================           
     do k=1,nlev
        ! theta - tendency on levels
        if (theta_advect_form==0) then
           v_theta(:,:,1,k)=elem(ie)%state%v(:,:,1,k,n0)*vtheta_dp(:,:,k)
           v_theta(:,:,2,k)=elem(ie)%state%v(:,:,2,k,n0)*vtheta_dp(:,:,k)
           div_v_theta(:,:,k)=divergence_sphere(v_theta(:,:,:,k),deriv,elem(ie))
        else
           ! alternate form, non-conservative, better HS topography results
           v_theta(:,:,:,k) = gradient_sphere(vtheta(:,:,k),deriv,elem(ie)%Dinv)
           div_v_theta(:,:,k)=vtheta(:,:,k)*divdp(:,:,k) + &
                dp3d(:,:,k)*elem(ie)%state%v(:,:,1,k,n0)*v_theta(:,:,1,k) + &
                dp3d(:,:,k)*elem(ie)%state%v(:,:,2,k,n0)*v_theta(:,:,2,k) 
        endif
#ifdef HOMMEXX_BFB_TESTING
        theta_tens(:,:,k)=(-theta_vadv(:,:,k)-div_v_theta(:,:,k))
#else
        theta_tens(:,:,k)=(-theta_vadv(:,:,k)-div_v_theta(:,:,k))*scale1
#endif

        ! w vorticity correction term
        temp(:,:,k) = (elem(ie)%state%w_i(:,:,k,n0)**2 + &
             elem(ie)%state%w_i(:,:,k+1,n0)**2)/4
        wvor(:,:,:,k) = gradient_sphere(temp(:,:,k),deriv,elem(ie)%Dinv)
        wvor(:,:,1,k) = wvor(:,:,1,k) - (elem(ie)%state%w_i(:,:,k,n0)*gradw_i(:,:,1,k) +&
             elem(ie)%state%w_i(:,:,k+1,n0)*gradw_i(:,:,1,k+1))/2
        wvor(:,:,2,k) = wvor(:,:,2,k) - (elem(ie)%state%w_i(:,:,k,n0)*gradw_i(:,:,2,k) +&
             elem(ie)%state%w_i(:,:,k+1,n0)*gradw_i(:,:,2,k+1))/2

        KE(:,:,k) = ( elem(ie)%state%v(:,:,1,k,n0)**2 + elem(ie)%state%v(:,:,2,k,n0)**2)/2
        gradKE(:,:,:,k) = gradient_sphere(KE(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)
#if 0
        ! another form: (good results in dcmip2012 test2.0)  max=0.195
        ! but bad results with HS topo
        !  grad(exner) =( grad(theta*exner) - exner*grad(theta))/theta
        vtemp(:,:,:,k) = gradient_sphere(vtheta(:,:,k)*exner(:,:,k),deriv,elem(ie)%Dinv)
        v_theta(:,:,:,k) = gradient_sphere(vtheta(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,1,k) = (vtemp(:,:,1,k)-exner(:,:,k)*v_theta(:,:,1,k))/&
             vtheta(:,:,k)
        gradexner(:,:,2,k) = (vtemp(:,:,2,k)-exner(:,:,k)*v_theta(:,:,2,k))/&
             vtheta(:,:,k)
#endif
#if 0
        ! entropy form: dcmip2012 test2.0 best: max=0.130  (0.124 with conservation form theta)
        vtemp(:,:,:,k) = gradient_sphere(vtheta(:,:,k)*exner(:,:,k),deriv,elem(ie)%Dinv)
        v_theta(:,:,:,k) = gradient_sphere(log(vtheta(:,:,k)),deriv,elem(ie)%Dinv)
        gradexner(:,:,1,k) = (vtemp(:,:,1,k)-exner(:,:,k)*vtheta(:,:,k)*v_theta(:,:,1,k))/&
             vtheta(:,:,k)
        gradexner(:,:,2,k) = (vtemp(:,:,2,k)-exner(:,:,k)*vtheta(:,:,k)*v_theta(:,:,2,k))/&
             vtheta(:,:,k)
#endif
#if 0
        ! another form:  terrible results in dcmip2012 test2.0
        ! grad(exner) = grad(p) * kappa * exner / p
        gradexner(:,:,:,k) = gradient_sphere(pnh(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,1,k) = gradexner(:,:,1,k)*(Rgas/Cp)*exner(:,:,k)/pnh(:,:,k)
        gradexner(:,:,2,k) = gradexner(:,:,2,k)*(Rgas/Cp)*exner(:,:,k)/pnh(:,:,k)
#endif

        ! special averaging of dpnh/dpi grad(phi) for E conservation
        mgrad(:,:,1,k) = (dpnh_dp_i(:,:,k)*gradphinh_i(:,:,1,k)+ &
              dpnh_dp_i(:,:,k+1)*gradphinh_i(:,:,1,k+1))/2
        mgrad(:,:,2,k) = (dpnh_dp_i(:,:,k)*gradphinh_i(:,:,2,k)+ &
              dpnh_dp_i(:,:,k+1)*gradphinh_i(:,:,2,k+1))/2

        if (pgrad_correction==1) then
           T0 = TREF-tref_lapse_rate*TREF*Cp/g     ! = 97  
#ifdef HOMMEXX_BFB_TESTING
           ! For BFB testing, calculate log(exner) using cxx_log()
           ! and then call gradient sphere.
           do j=1,np
             do i=1,np
               temp(i,j,k) = cxx_log(exner(i,j,k))
             end do
           end do

           vtemp(:,:,:,k)=gradient_sphere(temp(:,:,k),deriv,elem(ie)%Dinv)
#else
           vtemp(:,:,:,k)=gradient_sphere(log(exner(:,:,k)),deriv,elem(ie)%Dinv)
#endif
           mgrad(:,:,1,k)=mgrad(:,:,1,k) + Cp*T0*(vtemp(:,:,1,k)-gradexner(:,:,1,k)/exner(:,:,k))
           mgrad(:,:,2,k)=mgrad(:,:,2,k) + Cp*T0*(vtemp(:,:,2,k)-gradexner(:,:,2,k)/exner(:,:,k))
        endif


        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)

#ifdef HOMMEXX_BFB_TESTING
              vtens1(i,j,k) = ( - Cp*vtheta(i,j,k)*gradexner(i,j,1,k) &
                                - (v_vadv(i,j,1,k) + gradKE(i,j,1,k)) &
                                - (mgrad(i,j,1,k) + wvor(i,j,1,k))    &
                                + v2*(elem(ie)%fcor(i,j) + vort(i,j,k)) )

              vtens2(i,j,k) = ( - Cp*vtheta(i,j,k)*gradexner(i,j,2,k) &
                                - (v_vadv(i,j,2,k) + gradKE(i,j,2,k)) &
                                - (mgrad(i,j,2,k) + wvor(i,j,2,k))    &
                                - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) )
#else
              vtens1(i,j,k) = (-v_vadv(i,j,1,k) &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - gradKE(i,j,1,k) - mgrad(i,j,1,k) &
                  -Cp*vtheta(i,j,k)*gradexner(i,j,1,k)&
                  -wvor(i,j,1,k) )*scale1


              vtens2(i,j,k) = (-v_vadv(i,j,2,k) &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) &
                   - gradKE(i,j,2,k) - mgrad(i,j,2,k) &
                  -Cp*vtheta(i,j,k)*gradexner(i,j,2,k) &
                  -wvor(i,j,2,k) )*scale1
#endif
           end do
        end do     
     end do 



     
#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     ! diagnostics. not performance critical, dont thread
     ! =========================================================
     if (compute_diagnostics) then
        elem(ie)%accum%KEu_horiz1=0
        elem(ie)%accum%KEu_horiz2=0
        elem(ie)%accum%KEu_vert1=0
        elem(ie)%accum%KEu_vert2=0
        elem(ie)%accum%KEw_horiz1=0
        elem(ie)%accum%KEw_horiz2=0
        elem(ie)%accum%KEw_horiz3=0
        elem(ie)%accum%KEw_vert1=0
        elem(ie)%accum%KEw_vert2=0

        elem(ie)%accum%PEhoriz1=0
        elem(ie)%accum%PEhoriz2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%PEvert1=0
        elem(ie)%accum%PEvert2=0
        elem(ie)%accum%T01=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S2=0
        elem(ie)%accum%P1=0
        elem(ie)%accum%P2=0

        do k =1,nlev
          do j=1,np
            do i=1,np                
               d_eta_dot_dpdn_dn=(eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k))
               !  Form horiz advection of KE-u
               elem(ie)%accum%KEu_horiz1(i,j)=elem(ie)%accum%KEu_horiz1(i,j) &
                    -dp3d(i,j,k)*( &
                    elem(ie)%state%v(i,j,1,k,n0)*gradKE(i,j,1,k) + &
                    elem(ie)%state%v(i,j,2,k,n0)*gradKE(i,j,2,k) )
               elem(ie)%accum%KEu_horiz2(i,j)=elem(ie)%accum%KEu_horiz2(i,j)              &
                    -KE(i,j,k)*divdp(i,j,k)
               !  Form horiz advection of KE-w
               elem(ie)%accum%KEw_horiz1(i,j)=elem(ie)%accum%KEw_horiz1(i,j)-   &
                    dp3d(i,j,k) * (&
                    elem(ie)%state%w_i(i,j,k,n0) * v_gradw_i(i,j,k)    + &
                    elem(ie)%state%w_i(i,j,k+1,n0) * v_gradw_i(i,j,k+1) )/2
               elem(ie)%accum%KEw_horiz2(i,j)=elem(ie)%accum%KEw_horiz2(i,j)-   &
                    divdp(i,j,k)*(elem(ie)%state%w_i(i,j,k,n0)**2 + &
                    elem(ie)%state%w_i(i,j,k+1,n0)**2 ) /4
               elem(ie)%accum%KEw_horiz3(i,j)=elem(ie)%accum%KEw_horiz3(i,j)   &
                    -dp3d(i,j,k) * (elem(ie)%state%v(i,j,1,k,n0) * wvor(i,j,1,k) +  &
                    elem(ie)%state%v(i,j,2,k,n0) * wvor(i,j,2,k))
               !  Form vertical advection of KE-u 
               elem(ie)%accum%KEu_vert1(i,j)=elem(ie)%accum%KEu_vert1(i,j)- &
                    (elem(ie)%state%v(i,j,1,k,n0) * v_vadv(i,j,1,k) +            &
                    elem(ie)%state%v(i,j,2,k,n0) *v_vadv(i,j,2,k))*dp3d(i,j,k)
               elem(ie)%accum%KEu_vert2(i,j)=elem(ie)%accum%KEu_vert2(i,j)- &
                    0.5*((elem(ie)%state%v(i,j,1,k,n0))**2 +                     &
                    (elem(ie)%state%v(i,j,2,k,n0))**2)*d_eta_dot_dpdn_dn
               !  Form vertical advection of KE-w
               elem(ie)%accum%KEw_vert1(i,j)=elem(ie)%accum%KEw_vert1(i,j) - &
                    dp3d(i,j,k) * &
                    (w_vadv_i(i,j,k)*elem(ie)%state%w_i(i,j,k,n0)+ &
                    w_vadv_i(i,j,k+1)*elem(ie)%state%w_i(i,j,k+1,n0))/2
               
               elem(ie)%accum%KEw_vert2(i,j)=elem(ie)%accum%KEw_vert2(i,j)      &
                    -d_eta_dot_dpdn_dn* &
                    (.5*elem(ie)%state%w_i(i,j,k,n0)**2 +&
                    .5*elem(ie)%state%w_i(i,j,k+1,n0)**2)/2
               
               !  Form IEvert1
               elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j)      &
                    -Cp*exner(i,j,k)*theta_vadv(i,j,k)                        
               ! Form IEvert2 
               ! here use of dpnh_dp_i on boundry (with incorrect data)
               ! is harmess becuase eta_dot_dpdn=0
               elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j)      &
                    + ( dpnh_dp_i(i,j,k)*eta_dot_dpdn(i,j,k)+ &
                        dpnh_dp_i(i,j,k+1)*eta_dot_dpdn(i,j,k+1)) &
                    *(phi_i(i,j,k+1)-phi_i(i,j,k))/2
               
               !  Form PEhoriz1
               elem(ie)%accum%PEhoriz1(i,j)=(elem(ie)%accum%PEhoriz1(i,j))  &
                    -phi(i,j,k)*divdp(i,j,k) 
               !  Form PEhoriz2
               elem(ie)%accum%PEhoriz2(i,j)=elem(ie)%accum%PEhoriz2(i,j)    &
                    -dp3d(i,j,k)* &
                    (elem(ie)%state%v(i,j,1,k,n0)*                          &
                    (gradphinh_i(i,j,1,k)+gradphinh_i(i,j,1,k+1))/2  +      &
                    elem(ie)%state%v(i,j,2,k,n0)*                           &
                    (gradphinh_i(i,j,2,k)+gradphinh_i(i,j,2,k+1))/2  )
               
               !  Form PEvert1
               elem(ie)%accum%PEvert1(i,j) = elem(ie)%accum%PEvert1(i,j)    &
                    -phi(i,j,k)*d_eta_dot_dpdn_dn                                 
               elem(ie)%accum%PEvert2(i,j) = elem(ie)%accum%PEvert2(i,j)     &
                    -dp3d(i,j,k)*(phi_vadv_i(i,j,k)+phi_vadv_i(i,j,k+1))/2
               
               !  Form T01
               elem(ie)%accum%T01(i,j)=elem(ie)%accum%T01(i,j)               &
                    -(Cp*elem(ie)%state%vtheta_dp(i,j,k,n0))                       &
                    *(gradexner(i,j,1,k)*elem(ie)%state%v(i,j,1,k,n0) +           &
                    gradexner(i,j,2,k)*elem(ie)%state%v(i,j,2,k,n0))              
               !  Form S1 
               elem(ie)%accum%S1(i,j)=elem(ie)%accum%S1(i,j)                 &
                    -Cp*exner(i,j,k)*div_v_theta(i,j,k)

               !  Form P1  = -P2  (no reason to compute P2?)
               elem(ie)%accum%P1(i,j)=elem(ie)%accum%P1(i,j) -g*dp3d(i,j,k)* &
                    ( elem(ie)%state%w_i(i,j,k,n0) + &
                    elem(ie)%state%w_i(i,j,k+1,n0) )/2
               !  Form P2
               elem(ie)%accum%P2(i,j)=elem(ie)%accum%P2(i,j) + g*dp3d(i,j,k)*&
                    ( elem(ie)%state%w_i(i,j,k,n0) + &
                    elem(ie)%state%w_i(i,j,k+1,n0) )/2
            enddo
         enddo
      enddo

      ! these terms are better easier to compute by summing interfaces
      do k=2,nlev
         elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
              (g*elem(ie)%state%w_i(:,:,k,n0)-v_gradphinh_i(:,:,k)) &
               * dpnh_dp_i(:,:,k)*dp3d_i(:,:,k)
      enddo
      ! boundary terms
      do k=1,nlevp,nlev
         elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
           (g*elem(ie)%state%w_i(:,:,k,n0)-v_gradphinh_i(:,:,k)) &
           * dpnh_dp_i(:,:,k)*dp3d_i(:,:,k)/2
      enddo
      ! boundary term is incorrect.  save the term so we can correct it
      ! once we have coorect value of dpnh_dp_i:
      elem(ie)%accum%T2_nlevp_term(:,:)=&
           (g*elem(ie)%state%w_i(:,:,nlevp,n0)-v_gradphinh_i(:,:,nlevp)) &
           * dp3d_i(:,:,nlevp)/2

   endif
#endif


     do k=1,nlev
#ifdef HOMMEXX_BFB_TESTING
        elem(ie)%state%v(:,:,1,k,np1) = (elem(ie)%spheremp(:,:)*scale3) * elem(ie)%state%v(:,:,1,k,nm1) &
          + (scale1*dt2*elem(ie)%spheremp(:,:))*vtens1(:,:,k)
        elem(ie)%state%v(:,:,2,k,np1) = (elem(ie)%spheremp(:,:)*scale3) * elem(ie)%state%v(:,:,2,k,nm1) &
          +  (scale1*dt2*elem(ie)%spheremp(:,:))*vtens2(:,:,k)
        elem(ie)%state%vtheta_dp(:,:,k,np1) = (elem(ie)%spheremp(:,:)*scale3) * elem(ie)%state%vtheta_dp(:,:,k,nm1) &
          + (scale1*dt2*elem(ie)%spheremp(:,:))*theta_tens(:,:,k)
        if ( .not. theta_hydrostatic_mode ) then
           elem(ie)%state%w_i(:,:,k,np1)    = (elem(ie)%spheremp(:,:)*scale3) * elem(ie)%state%w_i(:,:,k,nm1)   &
                + (elem(ie)%spheremp(:,:)*dt2)*w_tens(:,:,k)
           elem(ie)%state%phinh_i(:,:,k,np1)   = (elem(ie)%spheremp(:,:)*scale3) * elem(ie)%state%phinh_i(:,:,k,nm1) &
                + (elem(ie)%spheremp(:,:)*dt2)*phi_tens(:,:,k)
        endif

        elem(ie)%state%dp3d(:,:,k,np1) = &
             (elem(ie)%spheremp(:,:) * scale3) * elem(ie)%state%dp3d(:,:,k,nm1) - &
             (scale1*dt2*elem(ie)%spheremp(:,:)) * (divdp(:,:,k) + (eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))

#else
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%v(:,:,1,k,nm1) &
          + dt2*vtens1(:,:,k) )
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%v(:,:,2,k,nm1) &
          +  dt2*vtens2(:,:,k) )
        elem(ie)%state%vtheta_dp(:,:,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%vtheta_dp(:,:,k,nm1) &
          + dt2*theta_tens(:,:,k))

        if ( .not. theta_hydrostatic_mode ) then
           elem(ie)%state%w_i(:,:,k,np1)    = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%w_i(:,:,k,nm1)   &
                + dt2*w_tens(:,:,k))
           elem(ie)%state%phinh_i(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%phinh_i(:,:,k,nm1) & 
                + dt2*phi_tens(:,:,k))
        endif

        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (scale3 * elem(ie)%state%dp3d(:,:,k,nm1) - &
             scale1*dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
#endif
     enddo
     if ( .not. theta_hydrostatic_mode ) then
        k=nlevp
#ifdef HOMMEXX_BFB_TESTING
        elem(ie)%state%w_i(:,:,k,np1)= (elem(ie)%spheremp(:,:)*scale3) * elem(ie)%state%w_i(:,:,k,nm1)   &
        + (elem(ie)%spheremp(:,:)*dt2)*w_tens(:,:,k)
#else
        elem(ie)%state%w_i(:,:,k,np1)=elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%w_i(:,:,k,nm1)   &
        + dt2*w_tens(:,:,k))
#endif
     endif


     kptr=0
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,nlyr_tot)
     kptr=kptr+2*nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, nlyr_tot)
     kptr=kptr+nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%vtheta_dp(:,:,:,np1),nlev,kptr,nlyr_tot)
     if (.not. theta_hydrostatic_mode) then
        kptr=kptr+nlev
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%w_i(:,:,:,np1),nlevp,kptr,nlyr_tot)
        kptr=kptr+nlevp
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%phinh_i(:,:,:,np1),nlev,kptr,nlyr_tot)
     endif

   end do ! end do for the ie=nets,nete loop

  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge_g)
  call t_stopf('caar_bexchV')

  do ie=nets,nete
     kptr=0
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,nlyr_tot)
     kptr=kptr+2*nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,nlyr_tot)
     kptr=kptr+nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%vtheta_dp(:,:,:,np1),nlev,kptr,nlyr_tot)
     if (.not. theta_hydrostatic_mode) then
        kptr=kptr+nlev
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%w_i(:,:,:,np1),nlevp,kptr,nlyr_tot)
        kptr=kptr+nlevp
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%phinh_i(:,:,:,np1),nlev,kptr,nlyr_tot)
     endif
      
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%vtheta_dp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%vtheta_dp(:,:,k,np1)
        if ( .not. theta_hydrostatic_mode ) then
           elem(ie)%state%w_i(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w_i(:,:,k,np1)
           elem(ie)%state%phinh_i(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%phinh_i(:,:,k,np1)
        endif
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
     k=nlevp
     if ( .not. theta_hydrostatic_mode ) &
          elem(ie)%state%w_i(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w_i(:,:,k,np1)


     ! now we can compute the correct dphn_dp_i() at the surface:
     if (.not. theta_hydrostatic_mode) then
        ! solve for (dpnh_dp_i-1)
        dpnh_dp_i(:,:,nlevp) = 1 + (  &
             ((elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g - &
             elem(ie)%state%w_i(:,:,nlevp,np1)) / &
             (g + ( elem(ie)%derived%gradphis(:,:,1)**2 + &
             elem(ie)%derived%gradphis(:,:,2)**2)/(2*g))   )  / dt2

        ! update solution with new dpnh_dp_i value:
        elem(ie)%state%w_i(:,:,nlevp,np1) = elem(ie)%state%w_i(:,:,nlevp,np1) +&
             scale1*dt2*g*(dpnh_dp_i(:,:,nlevp)-1)
        elem(ie)%state%v(:,:,1,nlev,np1) =  elem(ie)%state%v(:,:,1,nlev,np1) -&
             scale1*dt2*(dpnh_dp_i(:,:,nlevp)-1)*elem(ie)%derived%gradphis(:,:,1)/2
        elem(ie)%state%v(:,:,2,nlev,np1) =  elem(ie)%state%v(:,:,2,nlev,np1) -&
             scale1*dt2*(dpnh_dp_i(:,:,nlevp)-1)*elem(ie)%derived%gradphis(:,:,2)/2

#ifdef ENERGY_DIAGNOSTICS
        ! add in boundary term to T2 and S2 diagnostics:
        if (compute_diagnostics) then
           elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
                elem(ie)%accum%T2_nlevp_term(:,:)*(dpnh_dp_i(:,:,nlevp)-1)
           elem(ie)%accum%S2(:,:)=-elem(ie)%accum%T2(:,:)      
        endif

        ! check w b.c.
        temp(:,:,1) =  (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
        do j=1,np
        do i=1,np
           if ( abs(temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,np1)) >1e-10) then
              write(iulog,*) 'WARNING:CAAR w(np1) does not satisfy b.c.',ie,i,j,k
              write(iulog,*) 'val1 = ',temp(i,j,1)
              write(iulog,*) 'val2 = ',elem(ie)%state%w_i(i,j,nlevp,np1)
              write(iulog,*) 'diff: ',temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,np1)
           endif
        enddo
        enddo

        ! check for layer spacing <= 1m
        if (scale3 /= 0) then
        do k=1,nlev
        do j=1,np
        do i=1,np
           if ((elem(ie)%state%phinh_i(i,j,k,np1)-elem(ie)%state%phinh_i(i,j,k+1,np1)) < g) then
              write(iulog,*) 'WARNING:CAAR after ADV, delta z < 1m. ie,i,j,k=',ie,i,j,k
              write(iulog,*) 'phi(i,j,k)=  ',elem(ie)%state%phinh_i(i,j,k,np1)
              write(iulog,*) 'phi(i,j,k+1)=',elem(ie)%state%phinh_i(i,j,k+1,np1)
           endif
        enddo
        enddo
        enddo
        endif
#endif
     endif
     if (scale3 /= 0) then
       call limiter_dp3d_k(elem(ie)%state%dp3d(:,:,:,np1),elem(ie)%state%vtheta_dp(:,:,:,np1),&
            elem(ie)%spheremp,hvcoord%dp0)
     endif
  end do
  call t_stopf('compute_andor_apply_rhs')

  end subroutine compute_andor_apply_rhs


  subroutine limiter_dp3d_k(dp3d,vtheta_dp,spheremp,dp0)
  ! mass conserving column limiter (1D only)
  !
  ! if dp3d < dp3d_thresh*hvcoord%dp0 then apply vertical mixing 
  ! to prevent layer from getting too thin
  !
  ! This is rarely triggered and is mostly for safety when using 
  ! long remap timesteps
  !
  implicit none
  real (kind=real_kind), intent(inout) :: dp3d(np,np,nlev)
  real (kind=real_kind), intent(inout) :: vtheta_dp(np,np,nlev)
  real (kind=real_kind), intent(in) :: dp0(nlev)
  real (kind=real_kind), intent(in) :: spheremp(np,np)  !  density

  ! local
  real (kind=real_kind) :: Qcol(nlev)
  real (kind=real_kind) :: mass,mass_new
  logical :: warn
  integer i,j,k

  ! first check if limter is needed, and print warning
  warn=.false. 
  do k=1,nlev
     if ( minval(dp3d(:,:,k)) < dp3d_thresh*dp0(k)) then
#ifndef HOMMEXX_BFB_TESTING
        ! In bfb unit tests, we use (semi-)random inputs, so we expect to hit this.
        ! Still, we don't want to fill up the console output
        write(iulog,*) 'WARNING:CAAR: dp3d too small. dt_remap may be too large'
        write(iulog,*) 'k,dp3d(k), dp0: ',k,minval(dp3d(:,:,k)),dp0(k)
#endif
        warn=.true.
     endif
  enddo

  if (warn) then
#ifndef HOMMEXX_BFB_TESTING
    vtheta_dp(:,:,:)=vtheta_dp(:,:,:)/dp3d(:,:,:)
#endif
    do j = 1 , np
       do i = 1 , np
          if ( minval(dp3d(i,j,:) - dp3d_thresh*dp0(:)) < 0 ) then
#ifdef HOMMEXX_BFB_TESTING
             vtheta_dp(i,j,:)=vtheta_dp(i,j,:)/dp3d(i,j,:)
#endif
             ! subtract min, multiply in by weights
             Qcol(:) = (dp3d(i,j,:) - dp3d_thresh*dp0(:))*spheremp(i,j)
             mass = 0
             do k = 1,nlev
                mass = mass + Qcol(k)
             enddo

             ! negative mass.  so reduce all postive values to zero
             ! then increase negative values as much as possible
             if ( mass < 0 ) Qcol = -Qcol
             mass_new = 0
             do k=1,nlev
                if ( Qcol(k) < 0 ) then
                   Qcol(k) = 0
                else
                   mass_new = mass_new + Qcol(k)
                endif
             enddo
             ! now scale the all positive values to restore mass
             if ( mass_new > 0 ) Qcol(:) = Qcol(:) * (abs(mass) / mass_new)
             if ( mass     < 0 ) Qcol(:) = -Qcol(:)
             !
             dp3d(i,j,:) = Qcol(:)/spheremp(i,j) + dp3d_thresh*dp0(:)
#ifdef HOMMEXX_BFB_TESTING
             vtheta_dp(i,j,:)=vtheta_dp(i,j,:)*dp3d(i,j,:)
#endif
          endif
       enddo
    enddo
#ifndef HOMMEXX_BFB_TESTING
    vtheta_dp(:,:,:)=vtheta_dp(:,:,:)*dp3d(:,:,:)
#endif
  endif

#if 1
  ! check for theta < 10K                                                                                                       
  warn=.false.
  do k=1,nlev
     if ( minval(vtheta_dp(:,:,k)-vtheta_thresh*dp3d(:,:,k))   <  0) then
#ifndef HOMMEXX_BFB_TESTING
        ! In bfb unit tests, we use (semi-)random inputs, so we expect to hit this.
        ! Still, we don't want to fill up the console output
        write(iulog,*) 'WARNING:CAAR: theta<',vtheta_thresh,' applying limiter'
        write(iulog,*) 'k,vtheta(k): ',k,minval(vtheta_dp(:,:,k)/dp3d(:,:,k))
#endif
        warn=.true.
     endif
  enddo
  if (warn) then
    do k=1,nlev
      do j = 1 , np
         do i = 1 , np
            if ( (vtheta_dp(i,j,k) - vtheta_thresh*dp3d(i,j,k)) < 0 ) then
               vtheta_dp(i,j,k)=vtheta_thresh*dp3d(i,j,k)
            endif
         enddo
      enddo
    enddo
  endif
#endif




  end subroutine limiter_dp3d_k

end module prim_advance_mod
