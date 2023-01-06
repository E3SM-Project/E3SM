module caar_interface

  use kinds,          only: real_kind 

  implicit none

contains

  subroutine init_caar_f90 (ne, hyai, hybi, hyam, hybm, dvv, mp, ps0) bind(c)
    use iso_c_binding,          only: c_int
    use thetal_test_interface,  only: init_f90
    use dimensions_mod,         only: nlev, nlevp, np
    use edge_mod_base,          only: initEdgeBuffer, edge_g
    use geometry_interface_mod, only: par, elem
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: ne
    real (kind=real_kind), intent(in) :: hyai(nlevp), hybi(nlevp), hyam(nlev), hybm(nlev)
    real (kind=real_kind), intent(in) :: ps0
    real (kind=real_kind), intent(out) :: dvv(np,np), mp(np,np)
    !
    ! Locals
    !

    call init_f90(ne, hyai, hybi, hyam, hybm, dvv, mp, ps0)

    ! single global edge buffer for all models:
    ! hydrostatic 4*nlev      NH:  6*nlev+1  
    ! if this is too small, code will abort with an error message
    call initEdgeBuffer(par,edge_g,elem,6*nlev+1)
  end subroutine init_caar_f90

  subroutine run_caar_f90 (nm1, n0, np1, dt, eta_ave_w, scale1, scale2, scale3, &
                           hydrostatic, adv_conservative, rsplit_in,            &
                           dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr,    &
                           vn0_ptr, etadot_dpdn_ptr, omega_p_ptr) bind(c)
    use iso_c_binding,          only: c_ptr, c_f_pointer, c_int, c_bool
    use control_mod,            only: theta_hydrostatic_mode, theta_advect_form, rsplit
    use dimensions_mod,         only: nelemd, nlev, nlevp, np
    use prim_advance_mod,       only: compute_andor_apply_rhs
    use thetal_test_interface,  only: deriv, hvcoord
    use geometry_interface_mod, only: hybrid, elem
    use element_state,          only: timelevels
    !
    ! Input(s)
    !
    real (kind=real_kind), intent(in) :: dt, eta_ave_w, scale1, scale2, scale3
    integer (kind=c_int),  intent(in) :: nm1, n0, np1, rsplit_in
    logical (kind=c_bool), intent(in) :: hydrostatic, adv_conservative
    type (c_ptr),          intent(in) :: dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr
    type (c_ptr),          intent(in) :: vn0_ptr, etadot_dpdn_ptr, omega_p_ptr
    !
    ! Locals
    !
    real (kind=real_kind), pointer :: dp           (:,:,:,:,:)
    real (kind=real_kind), pointer :: vtheta_dp    (:,:,:,:,:)
    real (kind=real_kind), pointer :: w_i          (:,:,:,:,:)
    real (kind=real_kind), pointer :: phi_i        (:,:,:,:,:)
    real (kind=real_kind), pointer :: v            (:,:,:,:,:,:)
    real (kind=real_kind), pointer :: vn0          (:,:,:,:,:)
    real (kind=real_kind), pointer :: eta_dot_dpdn (:,:,:,:)
    real (kind=real_kind), pointer :: omega_p      (:,:,:,:)
    integer :: ie

    call c_f_pointer(dp_ptr,        dp,        [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(vtheta_dp_ptr, vtheta_dp, [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(w_i_ptr,       w_i,       [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(phi_i_ptr,     phi_i,     [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(v_ptr,         v,         [np,np,2,nlev,  timelevels, nelemd])

    call c_f_pointer(vn0_ptr,         vn0,          [np,np,2,nlev, nelemd])
    call c_f_pointer(etadot_dpdn_ptr, eta_dot_dpdn, [np,np,  nlevp,nelemd])
    call c_f_pointer(omega_p_ptr,     omega_p,      [np,np,  nlev, nelemd])

    do ie=1,nelemd
      ! Copy inputs in the elem%state
      elem(ie)%state%dp3d(:,:,:,:) = dp(:,:,:,:,ie)
      elem(ie)%state%vtheta_dp(:,:,:,:) = vtheta_dp(:,:,:,:,ie)
      elem(ie)%state%w_i(:,:,:,:) = w_i(:,:,:,:,ie)
      elem(ie)%state%phinh_i(:,:,:,:) = phi_i(:,:,:,:,ie)
      elem(ie)%state%v(:,:,:,:,:) = v(:,:,:,:,:,ie)

      elem(ie)%derived%vn0(:,:,:,:)        = vn0(:,:,:,:,ie)
      elem(ie)%derived%eta_dot_dpdn(:,:,:) = eta_dot_dpdn(:,:,:,ie)
      elem(ie)%derived%omega_p(:,:,:)      = omega_p(:,:,:,ie)
    enddo

    ! set control variables
    rsplit = rsplit_in
    theta_hydrostatic_mode = hydrostatic
    if (adv_conservative) then
      theta_advect_form = 0
    else
      theta_advect_form = 1
    endif

    call compute_andor_apply_rhs(np1,nm1,n0,dt,            &
                                 elem,hvcoord,hybrid,deriv,  &
                                 1,nelemd,.false.,eta_ave_w, &
                                 scale1,scale2,scale3)

    do ie=1,nelemd
      ! Copy elem%state into outputs
      dp(:,:,:,:,ie)        = elem(ie)%state%dp3d(:,:,:,:)
      vtheta_dp(:,:,:,:,ie) = elem(ie)%state%vtheta_dp(:,:,:,:)
      w_i(:,:,:,:,ie)       = elem(ie)%state%w_i(:,:,:,:)
      phi_i(:,:,:,:,ie)     = elem(ie)%state%phinh_i(:,:,:,:)
      v(:,:,:,:,:,ie)       = elem(ie)%state%v(:,:,:,:,:)

      vn0(:,:,:,:,ie)         =  elem(ie)%derived%vn0(:,:,:,:)
      eta_dot_dpdn(:,:,:,ie)  =  elem(ie)%derived%eta_dot_dpdn(:,:,:)
      omega_p(:,:,:,ie)       =  elem(ie)%derived%omega_p(:,:,:)
    enddo

  end subroutine run_caar_f90

  subroutine run_limiter_f90 (np1, dp_ptr, vtheta_dp_ptr) bind(c)
    use iso_c_binding,          only: c_ptr, c_f_pointer, c_int, c_bool
    use control_mod,            only: theta_hydrostatic_mode, theta_advect_form, rsplit
    use dimensions_mod,         only: nelemd, nlev, np
    use prim_advance_mod,       only: limiter_dp3d_k
    use element_state,          only: timelevels
    use geometry_interface_mod, only: elem
    use thetal_test_interface,  only: hvcoord
    !
    ! Input(s)
    !
    integer (kind=c_int),  intent(in) :: np1
    type (c_ptr),          intent(in) :: dp_ptr, vtheta_dp_ptr
    !
    ! Locals
    !
    real (kind=real_kind), pointer :: dp           (:,:,:,:,:)
    real (kind=real_kind), pointer :: vtheta_dp    (:,:,:,:,:)
    integer :: ie

    call c_f_pointer(dp_ptr,        dp,        [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(vtheta_dp_ptr, vtheta_dp, [np,np,  nlev,  timelevels, nelemd])

    do ie=1,nelemd
      call limiter_dp3d_k(dp(:,:,:,np1,ie),vtheta_dp(:,:,:,np1,ie),elem(ie)%spheremp,hvcoord%dp0)
    enddo

  end subroutine run_limiter_f90

end module caar_interface
