module hv_interface

  use iso_c_binding,  only: c_int, c_bool, c_ptr, c_f_pointer
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: nlev, nlevp, np
  use edgetype_mod,   only: EdgeBuffer_t
  use element_mod,    only: element_t
  use kinds,          only: real_kind 
  use hybvcoord_mod,  only: hvcoord_t 
  use parallel_mod,   only: abortmp

  implicit none

  type(hvcoord_t) :: hvcoord
  type (EdgeBuffer_t), public :: edge
  type (derivative_t) :: deriv

  public :: init_f90
  public :: biharmonic_wk_theta_f90

contains

  subroutine init_f90 (ne, hyai, dvv, mp, ps0, hv_subcycle, &
                       hv_nu, hv_nu_div, hydrostatic) bind(c)
    use control_mod,    only: hypervis_subcycle, nu, nu_div, theta_hydrostatic_mode
    use dimensions_mod, only: nelemd
    use geometry_interface_mod, only: initmp_f90, init_cube_geometry_f90, init_connectivity_f90
    use geometry_interface_mod, only: par, elem
    use edge_mod_base, only: initEdgeBuffer
    use element_state, only: allocate_element_arrays, setup_element_pointers_ie
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: ne, hv_subcycle
    real (kind=real_kind), intent(in) :: hyai(nlevp)
    real (kind=real_kind), intent(in) :: dvv(np,np), mp(np,np)
    real (kind=real_kind), intent(in) :: ps0, hv_nu, hv_nu_div
    logical (kind=c_bool), intent(in) :: hydrostatic
    !
    ! Locals
    !
    integer :: ie

    call initmp_f90()
    call init_cube_geometry_f90(ne)
    call init_connectivity_f90()
    call allocate_element_arrays(nelemd)

    do ie=1,nelemd
      call setup_element_pointers_ie(ie,elem(ie)%state, elem(ie)%derived, elem(ie)%accum)
      elem(ie)%mp = mp
    enddo

    hvcoord%hyai = hyai
    hvcoord%ps0 = ps0
    deriv%dvv = dvv
    hypervis_subcycle = hv_subcycle
    nu = hv_nu
    nu_div = hv_nu_div
    theta_hydrostatic_mode = hydrostatic

    ! There are 4 scalar fields (dp,theta,phi,w) and 1 vector field (v)
    ! Recall that the interface quantities are NOT exchanged at the surface
    call initEdgeBuffer(par,edge,elem,6*nlev)

  end subroutine init_f90

  subroutine init_geo_views_f90 (d_ptr, dinv_ptr, &
                       spmp_ptr, rspmp_ptr, tVisc_ptr,     &
                       sph2c_ptr,mdet_ptr,minv_ptr) bind(c)
    use dimensions_mod, only: nelemd
    use geometry_interface_mod, only: par, elem
    use edge_mod_base, only: initEdgeBuffer
    use element_state, only: allocate_element_arrays, setup_element_pointers_ie
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: d_ptr, dinv_ptr, spmp_ptr, rspmp_ptr, tVisc_ptr
    type (c_ptr), intent(in) :: sph2c_ptr, mdet_ptr, minv_ptr
    !
    ! Locals
    !
    integer :: ie
    real (kind=real_kind), pointer :: scalar2d (:,:,:)
    real (kind=real_kind), pointer :: tensor2d (:,:,:,:,:)

    ! Set all geometric views (we use 1 tensor and 1 scalar temps)
    call c_f_pointer(d_ptr, tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(mdet_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      elem(ie)%d = tensor2d(:,:,:,:,ie)
      elem(ie)%metdet = scalar2d(:,:,ie)
    enddo

    call c_f_pointer(dinv_ptr, tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(spmp_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      elem(ie)%dinv = tensor2d(:,:,:,:,ie)
      elem(ie)%spheremp = scalar2d(:,:,ie)
    enddo

    call c_f_pointer(minv_ptr, tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(rspmp_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      elem(ie)%metinv = tensor2d(:,:,:,:,ie)
      elem(ie)%rspheremp = scalar2d(:,:,ie)
    enddo

    call c_f_pointer(sph2c_ptr, tensor2d, [np,np,3,2,nelemd])
    call c_f_pointer(mdet_ptr, scalar2d, [np,np,nelemd])
    do ie=1,nelemd
      elem(ie)%vec_sphere2cart = tensor2d(:,:,:,:,ie)
      elem(ie)%metdet = scalar2d(:,:,ie)
      elem(ie)%rmetdet = 1.0 / elem(ie)%metdet
    enddo

    call c_f_pointer(tVisc_ptr, tensor2d, [np,np,2,2,nelemd])
    do ie=1,nelemd
      elem(ie)%tensorVisc = tensor2d(:,:,:,:,ie)
    enddo

  end subroutine init_geo_views_f90

  subroutine biharmonic_wk_theta_f90(np1, hv_scaling, dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr, &
                                     dptens_ptr, ttens_ptr, wtens_ptr, phitens_ptr, vtens_ptr) bind(c)
    use control_mod,            only: hypervis_scaling
    use dimensions_mod,         only: nelemd
    use element_state,          only: timelevels
    use geometry_interface_mod, only: par,elem,hybrid
    use hybrid_mod,             only: hybrid_t, hybrid_create
    use viscosity_theta,        only: biharmonic_wk_theta
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: np1
    real (kind=real_kind), intent(in) :: hv_scaling
    type (c_ptr), intent(in) :: dp_ptr, vtheta_dp_ptr, w_i_ptr, phi_i_ptr, v_ptr
    type (c_ptr), intent(in) :: dptens_ptr, ttens_ptr, wtens_ptr, phitens_ptr, vtens_ptr

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
    call biharmonic_wk_theta(elem,stens,vtens,deriv,edge,hybrid,np1,1,nelemd)

    do ie=1,nelemd
      ! Copy outputs from stens into the individual *tens variables
      dptens(:,:,:,ie)  = stens(:,:,:,1,ie)
      ttens(:,:,:,ie)   = stens(:,:,:,2,ie)
      wtens(:,:,:,ie)   = stens(:,:,:,3,ie)
      phitens(:,:,:,ie) = stens(:,:,:,4,ie)
    enddo

  end subroutine biharmonic_wk_theta_f90

  subroutine save_state_f90 () bind(c)
  end subroutine save_state_f90

  subroutine cleanup_f90 () bind(c)
    use edge_mod_base, only : FreeEdgeBuffer
    use element_state, only : deallocate_element_arrays
    use geometry_interface_mod, only: cleanup_geometry_f90
    use parallel_mod, only: rrequest, srequest, status

    call FreeEdgeBuffer(edge)
    call cleanup_geometry_f90()

    ! Deallocate all module-scope variables, so the next iteration of catch2
    ! can successfully call init_f90 (which will try to allocate these again)
    deallocate(rrequest)
    deallocate(srequest)
    deallocate(status)

    call deallocate_element_arrays()

  end subroutine cleanup_f90

end module hv_interface
