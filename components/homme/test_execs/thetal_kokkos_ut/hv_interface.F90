module hv_interface

  use iso_c_binding,  only: c_int, c_bool, c_double, c_ptr, c_f_pointer
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: nlev, nlevp, np
  use edge_mod_base,  only: edge_g
  use element_mod,    only: element_t
  use kinds,          only: real_kind 
  use hybvcoord_mod,  only: hvcoord_t 
  use parallel_mod,   only: abortmp

  implicit none

  type(hvcoord_t) :: hvcoord
  type (derivative_t) :: deriv

  public :: init_f90
  public :: biharmonic_wk_theta_f90

contains

  subroutine init_f90 (ne, hyai, hybi, hyam, hybm, dvv, mp, ps0, hv_subcycle, &
                       hv_nu, hv_nu_div, hv_nu_top, hv_nu_p, hv_nu_s) bind(c)
    use control_mod,    only: hypervis_subcycle, nu, nu_div, nu_top, nu_p, nu_s, cubed_sphere_map
    use cube_mod,       only: cube_init_atomic, set_corner_coordinates, set_area_correction_map0
    use derivative_mod, only: derivinit
    use dimensions_mod, only: nelemd
    use geometry_interface_mod, only: initmp_f90, init_cube_geometry_f90, init_connectivity_f90
    use geometry_interface_mod, only: par, elem
    use hybvcoord_mod, only: set_layer_locations
    use edge_mod_base, only: initEdgeBuffer
    use mass_matrix_mod, only: mass_matrix
    use quadrature_mod, only: gausslobatto, quadrature_t
    use element_state, only: allocate_element_arrays, setup_element_pointers_ie, nlev_tom, nu_scale_top
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
    type (quadrature_t) :: gp

    call derivinit(deriv)

    call initmp_f90()
    call init_cube_geometry_f90(ne)
    call init_connectivity_f90()

    cubed_sphere_map = 0
    gp=gausslobatto(np)  ! GLL points
    do ie=1,nelemd
      call set_corner_coordinates(elem(ie))
    end do
    do ie=1,nelemd
      call cube_init_atomic(elem(ie),gp%points)
    enddo
    call allocate_element_arrays(nelemd)

    call mass_matrix(par,elem)
    call set_area_correction_map0(elem, nelemd, par, gp)
    call mass_matrix(par,elem)

    ! Copy refFE matrices back to C
    mp = elem(1)%mp
    dvv = deriv%dvv

    do ie=1,nelemd
      call setup_element_pointers_ie(ie,elem(ie)%state, elem(ie)%derived, elem(ie)%accum)
    enddo

    hvcoord%hyai = hyai
    hvcoord%hybi = hybi
    hvcoord%hyam = hyam
    hvcoord%hybm = hybm
    hvcoord%ps0 = ps0
    call set_layer_locations (hvcoord,.false.,.false.)

    deriv%dvv = dvv
    hypervis_subcycle = hv_subcycle
    nu = hv_nu
    nu_div = hv_nu_div
    nu_top = hv_nu_top
    nu_p = hv_nu_p
    nu_s = hv_nu_s

    nlev_tom = 3
    nu_scale_top(1) = 4
    nu_scale_top(2) = 2
    nu_scale_top(3) = 1

    ! There are 4 scalar fields (dp,theta,phi,w) and 1 vector field (v)
    ! Recall that the interface quantities are NOT exchanged at the surface
    call initEdgeBuffer(par,edge_g,elem,6*nlev)

  end subroutine init_f90

  subroutine init_geo_views_f90 (d_ptr, dinv_ptr,       &
                       phis_ptr, gradphis_ptr,          &
                       spmp_ptr, rspmp_ptr, tVisc_ptr,  &
                       sph2c_ptr,mdet_ptr,minv_ptr) bind(c)
    use dimensions_mod, only: nelemd
    use geometry_interface_mod, only: par, elem
    use element_state, only: allocate_element_arrays, setup_element_pointers_ie
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: d_ptr, dinv_ptr, spmp_ptr, rspmp_ptr, tVisc_ptr
    type (c_ptr), intent(in) :: sph2c_ptr, mdet_ptr, minv_ptr, phis_ptr, gradphis_ptr
    !
    ! Locals
    !
    integer :: ie
    real (kind=real_kind), pointer :: scalar2d (:,:,:)
    real (kind=real_kind), pointer :: vector2d (:,:,:,:)
    real (kind=real_kind), pointer :: tensor2d (:,:,:,:,:)

    ! Set all geometric views (we use 1 tensor, 1 vector, and 1 scalar temps)
    call c_f_pointer(d_ptr, tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(mdet_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%d      
      scalar2d(:,:,ie)     = elem(ie)%metdet 
    enddo

    call c_f_pointer(dinv_ptr, tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(spmp_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%dinv
      scalar2d(:,:,ie)     = elem(ie)%spheremp
    enddo

    call c_f_pointer(minv_ptr, tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(rspmp_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%metinv
      scalar2d(:,:,ie)     = elem(ie)%rspheremp
    enddo

    call c_f_pointer(sph2c_ptr, tensor2d, [np,np,3,2,nelemd])
    call c_f_pointer(mdet_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%vec_sphere2cart
      scalar2d(:,:,ie)     = elem(ie)%metdet
    enddo

    call c_f_pointer(phis_ptr, scalar2d, [np,np,nelemd])
    call c_f_pointer(gradphis_ptr, vector2d, [np,np,2,nelemd])
    call c_f_pointer(tVisc_ptr, tensor2d, [np,np,2,2,nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%tensorVisc
      elem(ie)%derived%gradphis = vector2d(:,:,:,ie)
      elem(ie)%state%phis       = scalar2d(:,:,ie)
    enddo

  end subroutine init_geo_views_f90

  subroutine biharmonic_wk_theta_f90(np1, hv_scaling, hydrostatic, dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr, &
                                     dptens_ptr, ttens_ptr, wtens_ptr, phitens_ptr, vtens_ptr) bind(c)
    use control_mod,            only: hypervis_scaling, theta_hydrostatic_mode
    use dimensions_mod,         only: nelemd
    use element_state,          only: timelevels
    use geometry_interface_mod, only: par,elem,hybrid
    use hybrid_mod,             only: hybrid_t, hybrid_create
    use viscosity_theta,        only: biharmonic_wk_theta
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: np1
    real (kind=c_double), intent(in) :: hv_scaling
    type (c_ptr), intent(in) :: dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr
    type (c_ptr), intent(in) :: dptens_ptr, ttens_ptr, wtens_ptr, phitens_ptr, vtens_ptr
    logical (kind=c_bool), intent(in) :: hydrostatic

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
                                  v_ptr,w_ptr,vtheta_ptr,dp_ptr,phinh_ptr) bind(c)
    use control_mod,            only: hypervis_scaling, theta_hydrostatic_mode
    use prim_advance_mod,       only: advance_hypervis
    use geometry_interface_mod, only: elem, hybrid
    use dimensions_mod,         only: nelemd
    use element_state,          only: timelevels
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: np1
    type (c_ptr), intent(in) :: dp_ptr, vtheta_ptr, w_ptr, phinh_ptr, v_ptr
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

    hypervis_scaling = hv_scaling
    theta_hydrostatic_mode = hydrostatic

    call c_f_pointer(v_ptr,      v,      [np,np,2,nlev,  timelevels, nelemd])
    call c_f_pointer(w_ptr,      w,      [np,np,  nlevp, timelevels, nelemd])
    call c_f_pointer(vtheta_ptr, vtheta, [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(dp_ptr,     dp,     [np,np,  nlev,  timelevels, nelemd])
    call c_f_pointer(phinh_ptr,  phinh,  [np,np,  nlevp, timelevels, nelemd])

    do ie=1,nelemd
      ! Copy input-outputs in the elem%state
      elem(ie)%state%v(:,:,:,:,:)       = v(:,:,:,:,:,ie)
      elem(ie)%state%w_i(:,:,:,:)       = w(:,:,:,:,ie)
      elem(ie)%state%vtheta_dp(:,:,:,:) = vtheta(:,:,:,:,ie)
      elem(ie)%state%dp3d(:,:,:,:)      = dp(:,:,:,:,ie)
      elem(ie)%state%phinh_i(:,:,:,:)   = phinh(:,:,:,:,ie)
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

  subroutine cleanup_f90 () bind(c)
    use edge_mod_base, only : FreeEdgeBuffer
    use element_state, only : deallocate_element_arrays
    use geometry_interface_mod, only: cleanup_geometry_f90
    use parallel_mod, only: rrequest, srequest, status

    call FreeEdgeBuffer(edge_g)
    call cleanup_geometry_f90()

    ! Deallocate all module-scope variables, so the next iteration of catch2
    ! can successfully call init_f90 (which will try to allocate these again)
    deallocate(rrequest)
    deallocate(srequest)
    deallocate(status)

    call deallocate_element_arrays()

  end subroutine cleanup_f90

end module hv_interface
