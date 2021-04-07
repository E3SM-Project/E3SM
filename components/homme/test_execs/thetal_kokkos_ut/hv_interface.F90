module hv_interface

  use iso_c_binding,  only: c_int, c_bool, c_double, c_ptr, c_f_pointer
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: nlev, nlevp, np
  use element_mod,    only: element_t
  use kinds,          only: real_kind 
  use parallel_mod,   only: abortmp

  implicit none

  public :: init_hv_f90
  public :: biharmonic_wk_theta_f90
  public :: advance_hypervis_f90

contains

  subroutine init_hv_f90 (ne, hyai, hybi, hyam, hybm, dvv, mp, ps0, hv_subcycle, &
                          hv_nu, hv_nu_div, hv_nu_top, hv_nu_p, hv_nu_s) bind(c)
    use control_mod,    only: hypervis_subcycle, nu, nu_div, nu_top, nu_p, nu_s
    use thetal_test_interface, only: init_f90
    use edge_mod_base, only: initEdgeBuffer, edge_g
    use element_state, only: nlev_tom, nu_scale_top
    use geometry_interface_mod, only: par, elem
    use physical_constants,     only: scale_factor, scale_factor_inv, laplacian_rigid_factor, rearth, rrearth

    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: ne, hv_subcycle
    real (kind=real_kind), intent(in) :: hyai(nlevp), hybi(nlevp), hyam(nlev), hybm(nlev)
    real (kind=real_kind), intent(in) :: ps0, hv_nu, hv_nu_div, hv_nu_top, hv_nu_p, hv_nu_s
    real (kind=real_kind), intent(out) :: dvv(np,np), mp(np,np)
    !
    ! Locals
    !
    integer :: ie

    scale_factor = rearth
    scale_factor_inv = rrearth
    laplacian_rigid_factor = rrearth
    
    call init_f90(ne, hyai, hybi, hyam, hybm, dvv, mp, ps0)

    ! There are 4 scalar fields (dp,theta,phi,w) and 1 vector field (v)
    ! Recall that the interface quantities are NOT exchanged at the surface
    call initEdgeBuffer(par,edge_g,elem,6*nlev)

    hypervis_subcycle = hv_subcycle
    nu = hv_nu
    nu_div = hv_nu_div
    nu_top = hv_nu_top
    nu_p = hv_nu_p
    nu_s = hv_nu_s

    nlev_tom = 3
    nu_scale_top(1) = 4D0
    nu_scale_top(2) = 2D0
    nu_scale_top(3) = 1D0
  end subroutine init_hv_f90

  subroutine biharmonic_wk_theta_f90(np1, hv_scaling, hydrostatic, dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr, &
                                     dptens_ptr, ttens_ptr, wtens_ptr, phitens_ptr, vtens_ptr) bind(c)
    use control_mod,            only: hypervis_scaling, theta_hydrostatic_mode
    use dimensions_mod,         only: nelemd
    use edge_mod_base,          only: edge_g
    use element_state,          only: timelevels
    use geometry_interface_mod, only: par,elem,hybrid
    use hybrid_mod,             only: hybrid_t, hybrid_create
    use thetal_test_interface,  only: deriv
    use viscosity_theta,        only: biharmonic_wk_theta
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: np1
    real (kind=c_double), intent(in) :: hv_scaling
    type (c_ptr), intent(in) :: dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr
    type (c_ptr), intent(in) :: dptens_ptr, ttens_ptr, wtens_ptr, phitens_ptr, vtens_ptr
    logical (kind=c_bool), intent(in) :: hydrostatic
    !
    ! Locals
    !
    real (kind=real_kind), pointer :: dp        (:,:,:,:,:)
    real (kind=real_kind), pointer :: vtheta_dp (:,:,:,:,:)
    real (kind=real_kind), pointer :: w_i       (:,:,:,:,:)
    real (kind=real_kind), pointer :: phi_i     (:,:,:,:,:)
    real (kind=real_kind), pointer :: v         (:,:,:,:,:,:)

    real (kind=real_kind), pointer :: dptens    (:,:,:,:)
    real (kind=real_kind), pointer :: ttens     (:,:,:,:)
    real (kind=real_kind), pointer :: wtens     (:,:,:,:)
    real (kind=real_kind), pointer :: phitens   (:,:,:,:)
    real (kind=real_kind), pointer :: vtens     (:,:,:,:,:)
    real (kind=real_kind) :: stens (np,np,nlev,4,nelemd)

    integer :: ie

    hypervis_scaling = hv_scaling
    theta_hydrostatic_mode = hydrostatic

    call c_f_pointer(dp_ptr,        dp,        [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(vtheta_dp_ptr, vtheta_dp, [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(w_i_ptr,       w_i,       [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(phi_i_ptr,     phi_i,     [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(v_ptr,         v,         [np,np,2,nlev,  timelevels, nelemd])

    call c_f_pointer(dptens_ptr,    dptens,    [np,np,  nlev,nelemd])
    call c_f_pointer(ttens_ptr,     ttens,     [np,np,  nlev,nelemd])
    call c_f_pointer(wtens_ptr,     wtens,     [np,np,  nlev,nelemd])
    call c_f_pointer(phitens_ptr,   phitens,   [np,np,  nlev,nelemd])
    call c_f_pointer(vtens_ptr,     vtens,     [np,np,2,nlev,nelemd])

    do ie=1,nelemd
      ! Copy inputs in the elem%state
      elem(ie)%state%dp3d(:,:,:,:) = dp(:,:,:,:,ie)
      elem(ie)%state%vtheta_dp(:,:,:,:) = vtheta_dp(:,:,:,:,ie)
      elem(ie)%state%w_i(:,:,:,:) = w_i(:,:,:,:,ie)
      elem(ie)%state%phinh_i(:,:,:,:) = phi_i(:,:,:,:,ie)
      elem(ie)%state%v(:,:,:,:,:) = v(:,:,:,:,:,ie)
    enddo

    ! Call the biharmonic routine
    call biharmonic_wk_theta(elem,stens,vtens,deriv,edge_g,hybrid,np1,1,nelemd)

    do ie=1,nelemd
      ! Copy outputs from stens into the individual *tens variables
      dptens(:,:,:,ie)  = stens(:,:,:,1,ie)
      ttens(:,:,:,ie)   = stens(:,:,:,2,ie)
      wtens(:,:,:,ie)   = stens(:,:,:,3,ie)
      phitens(:,:,:,ie) = stens(:,:,:,4,ie)
    enddo

  end subroutine biharmonic_wk_theta_f90

  subroutine advance_hypervis_f90(np1, dt, eta_ave_w, hv_scaling, hydrostatic, &
                                  dp_ref_ptr, theta_ref_ptr, phi_ref_ptr,      &
                                  v_ptr,w_ptr,vtheta_ptr,dp_ptr,phinh_ptr) bind(c)
    use control_mod,            only: hypervis_scaling, theta_hydrostatic_mode
    use prim_advance_mod,       only: advance_hypervis
    use geometry_interface_mod, only: elem, hybrid
    use dimensions_mod,         only: nelemd
    use element_state,          only: timelevels
    use thetal_test_interface,  only: deriv, hvcoord
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: np1
    type (c_ptr), intent(in) :: dp_ptr, vtheta_ptr, w_ptr, phinh_ptr, v_ptr
    type (c_ptr), intent(in) :: dp_ref_ptr, theta_ref_ptr, phi_ref_ptr
    real (kind=c_double), intent(in) :: dt, eta_ave_w, hv_scaling
    logical (kind=c_bool), intent(in) :: hydrostatic
    !
    ! Locals
    !
    integer :: ie
    real (kind=real_kind), pointer :: dp     (:,:,:,:,:)
    real (kind=real_kind), pointer :: vtheta (:,:,:,:,:)
    real (kind=real_kind), pointer :: w      (:,:,:,:,:)
    real (kind=real_kind), pointer :: phinh  (:,:,:,:,:)
    real (kind=real_kind), pointer :: v      (:,:,:,:,:,:)

    real (kind=real_kind), pointer :: dp_ref    (:,:,:,:)
    real (kind=real_kind), pointer :: theta_ref (:,:,:,:)
    real (kind=real_kind), pointer :: phi_ref   (:,:,:,:)

    hypervis_scaling = hv_scaling
    theta_hydrostatic_mode = hydrostatic

    call c_f_pointer(v_ptr,      v,      [np,np,2,nlev,  timelevels, nelemd])
    call c_f_pointer(w_ptr,      w,      [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(vtheta_ptr, vtheta, [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(dp_ptr,     dp,     [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(phinh_ptr,  phinh,  [np,np,  nlevp, timelevels, nelemd])

    call c_f_pointer(dp_ref_ptr,    dp_ref,    [np,np, nlev,  nelemd])
    call c_f_pointer(theta_ref_ptr, theta_ref, [np,np, nlev,  nelemd])
    call c_f_pointer(phi_ref_ptr,   phi_ref,   [np,np, nlevp, nelemd])

    do ie=1,nelemd
      ! Copy input-outputs in the elem%state
      elem(ie)%state%v(:,:,:,:,:)       = v(:,:,:,:,:,ie)
      elem(ie)%state%w_i(:,:,:,:)       = w(:,:,:,:,ie)
      elem(ie)%state%vtheta_dp(:,:,:,:) = vtheta(:,:,:,:,ie)
      elem(ie)%state%dp3d(:,:,:,:)      = dp(:,:,:,:,ie)
      elem(ie)%state%phinh_i(:,:,:,:)   = phinh(:,:,:,:,ie)

      elem(ie)%derived%dp_ref(:,:,:)    = dp_ref(:,:,:,ie)
      elem(ie)%derived%theta_ref(:,:,:) = theta_ref(:,:,:,ie)
      elem(ie)%derived%phi_ref(:,:,:)   = phi_ref(:,:,:,ie)
    enddo

    call advance_hypervis(elem,hvcoord,hybrid,deriv,np1,1,nelemd,dt,eta_ave_w)

    do ie=1,nelemd
      ! Copy elem%state into input-outputs
      v(:,:,:,:,:,ie)    = elem(ie)%state%v(:,:,:,:,:)       
      w(:,:,:,:,ie)      = elem(ie)%state%w_i(:,:,:,:)       
      vtheta(:,:,:,:,ie) = elem(ie)%state%vtheta_dp(:,:,:,:) 
      dp(:,:,:,:,ie)     = elem(ie)%state%dp3d(:,:,:,:)      
      phinh(:,:,:,:,ie)  = elem(ie)%state%phinh_i(:,:,:,:)   
    enddo

  end subroutine advance_hypervis_f90

end module hv_interface
