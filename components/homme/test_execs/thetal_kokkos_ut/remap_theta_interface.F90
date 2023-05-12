module remap_interface

  use kinds,                  only: real_kind 
  use geometry_interface_mod, only: elem

  implicit none

  public :: init_phis_f90
  public :: run_remap_f90
contains

  subroutine init_remap_f90 (ne, hyai, hybi, hyam, hybm, dvv, mp, ps0) bind(c)
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

    ! I don't reall need it for this test, but cleanup_f90 deallocates the edge buffer,
    ! so not init-ing it would give a runtime error.
    call initEdgeBuffer(par,edge_g,elem,1)
  end subroutine init_remap_f90

  subroutine init_phis_f90 (phis_ptr, gradphis_ptr) bind(c)
    use iso_c_binding,  only: c_ptr, c_f_pointer
    use dimensions_mod, only: nelemd, np
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: phis_ptr, gradphis_ptr
    !
    ! Locals
    !
    integer :: ie
    real (kind=real_kind), pointer :: phis (:,:,:)
    real (kind=real_kind), pointer :: grad_phis (:,:,:,:)

    call c_f_pointer(phis_ptr,     phis,      [np,np,  nelemd])
    call c_f_pointer(gradphis_ptr, grad_phis, [np,np,2,nelemd])
    do ie=1,nelemd
      elem(ie)%derived%gradphis = grad_phis(:,:,:,ie)
      elem(ie)%state%phis       = phis(:,:,ie)
    enddo

  end subroutine init_phis_f90

  subroutine run_remap_f90 (np1, np1_qdp, dt, rsplit_in, qsize_in, vr_alg,    &
                            dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr, &
                            ps_ptr, eta_dot_dpdn_ptr, qdp_ptr) bind(c)
    use iso_c_binding,          only: c_ptr, c_f_pointer, c_int, c_bool
    use control_mod,            only: vert_remap_q_alg, vert_remap_u_alg, rsplit
    use dimensions_mod,         only: nelemd, nlev, nlevp, np, qsize, qsize_d
    use thetal_test_interface,  only: hvcoord
    use geometry_interface_mod, only: hybrid
    use element_state,          only: timelevels
    use vertremap_mod,          only: vertical_remap
    !
    ! Input(s)
    !
    real (kind=real_kind), intent(in) :: dt
    integer (kind=c_int),  intent(in) :: np1, np1_qdp, rsplit_in, vr_alg, qsize_in
    type (c_ptr),          intent(in) :: dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr
    type (c_ptr),          intent(in) :: ps_ptr, eta_dot_dpdn_ptr, qdp_ptr
    !
    ! Locals
    !
    real (kind=real_kind), pointer :: dp           (:,:,:,:,:)
    real (kind=real_kind), pointer :: vtheta_dp    (:,:,:,:,:)
    real (kind=real_kind), pointer :: w_i          (:,:,:,:,:)
    real (kind=real_kind), pointer :: phi_i        (:,:,:,:,:)
    real (kind=real_kind), pointer :: v            (:,:,:,:,:,:)
    real (kind=real_kind), pointer :: ps           (:,:,:,:)
    real (kind=real_kind), pointer :: eta_dot_dpdn (:,:,:,:)
    real (kind=real_kind), pointer :: qdp          (:,:,:,:,:,:)
    integer :: ie

    call c_f_pointer(dp_ptr,           dp,           [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(vtheta_dp_ptr,    vtheta_dp,    [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(w_i_ptr,          w_i,          [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(phi_i_ptr,        phi_i,        [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(v_ptr,            v,            [np,np,2,nlev,  timelevels, nelemd])
    call c_f_pointer(ps_ptr,           ps,           [np,np,         timelevels, nelemd])
    call c_f_pointer(eta_dot_dpdn_ptr, eta_dot_dpdn, [np,np,  nlevp,             nelemd])
    call c_f_pointer(qdp_ptr,          qdp,          [np,np,  nlev,  qsize_d, 2, nelemd])

    do ie=1,nelemd
      ! Copy inputs in the elem%state
      elem(ie)%state%dp3d(:,:,:,:) = dp(:,:,:,:,ie)
      elem(ie)%state%vtheta_dp(:,:,:,:) = vtheta_dp(:,:,:,:,ie)
      elem(ie)%state%w_i(:,:,:,:) = w_i(:,:,:,:,ie)
      elem(ie)%state%phinh_i(:,:,:,:) = phi_i(:,:,:,:,ie)
      elem(ie)%state%v(:,:,:,:,:) = v(:,:,:,:,:,ie)
      elem(ie)%state%ps_v(:,:,:) = ps(:,:,:,ie)
      elem(ie)%state%qdp(:,:,:,:,:) = qdp(:,:,:,:,:,ie)

      elem(ie)%derived%eta_dot_dpdn(:,:,:) = eta_dot_dpdn(:,:,:,ie)
    enddo

    ! set control variables
    rsplit = rsplit_in
    vert_remap_q_alg = vr_alg
    vert_remap_u_alg = vr_alg
    qsize = qsize_in

    ! Call f90 vertical remap
    call vertical_remap(hybrid,elem,hvcoord,dt,np1,np1_qdp,1,nelemd)

    do ie=1,nelemd
      ! Copy elem%state into outputs
      dp(:,:,:,:,ie)        = elem(ie)%state%dp3d(:,:,:,:)
      vtheta_dp(:,:,:,:,ie) = elem(ie)%state%vtheta_dp(:,:,:,:)
      w_i(:,:,:,:,ie)       = elem(ie)%state%w_i(:,:,:,:)
      phi_i(:,:,:,:,ie)     = elem(ie)%state%phinh_i(:,:,:,:)
      v(:,:,:,:,:,ie)       = elem(ie)%state%v(:,:,:,:,:)
      qdp(:,:,:,:,:,ie)     = elem(ie)%state%qdp(:,:,:,:,:)
    enddo

  end subroutine run_remap_f90

end module remap_interface
