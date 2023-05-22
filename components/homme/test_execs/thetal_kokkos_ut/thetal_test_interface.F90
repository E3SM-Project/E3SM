module thetal_test_interface

  use iso_c_binding,  only: c_int, c_bool, c_double, c_ptr, c_f_pointer
  use derivative_mod, only: derivative_t
  use element_mod,    only: element_t
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use parallel_mod,   only: abortmp
  use geometry_mod,   only: set_area_correction_map0

  implicit none

  type(hvcoord_t) :: hvcoord
  type (derivative_t) :: deriv

  public :: init_f90
  public :: cleanup_f90
  public :: init_geo_views_f90
  public :: initialize_reference_states_f90

contains

  subroutine init_f90 (ne, hyai, hybi, hyam, hybm, dvv, mp, ps0) bind(c)
    use control_mod,            only: cubed_sphere_map
    use cube_mod,               only: cube_init_atomic, set_corner_coordinates
    use derivative_mod,         only: derivinit
    use dimensions_mod,         only: nelemd, nlev, nlevp, np
    use geometry_interface_mod, only: initmp_f90, init_cube_geometry_f90, init_connectivity_f90
    use geometry_interface_mod, only: par, elem
    use hybvcoord_mod,          only: set_layer_locations
    use element_state,          only: allocate_element_arrays, setup_element_pointers_ie, &
                                      nlev_tom, nu_scale_top
    use mass_matrix_mod,        only: mass_matrix
    use quadrature_mod,         only: gausslobatto, quadrature_t
    use physical_constants,     only: scale_factor, scale_factor_inv, laplacian_rigid_factor, rearth, rrearth
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
    integer :: ie, k
    type (quadrature_t) :: gp

    scale_factor = rearth
    scale_factor_inv = rrearth
    laplacian_rigid_factor = rrearth

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
    do k=1,nlev
      hvcoord%dp0(k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*ps0 + &
                       (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps0
    enddo

    call set_layer_locations (hvcoord,.false.,.false.)

    deriv%dvv = dvv

  end subroutine init_f90

  subroutine init_geo_views_f90 (d_ptr, dinv_ptr,        &
                       phis_ptr, gradphis_ptr, fcor_ptr, &
                       spmp_ptr, rspmp_ptr, tVisc_ptr,   &
                       sph2c_ptr,mdet_ptr,minv_ptr) bind(c)
    use dimensions_mod, only: nelemd, np
    use geometry_interface_mod, only: elem
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: d_ptr, dinv_ptr, spmp_ptr, rspmp_ptr, tVisc_ptr, fcor_ptr
    type (c_ptr), intent(in) :: sph2c_ptr, mdet_ptr, minv_ptr, phis_ptr, gradphis_ptr
    !
    ! Locals
    !
    integer :: ie
    real (kind=real_kind), pointer :: scalar2d (:,:,:)
    real (kind=real_kind), pointer :: vector2d (:,:,:,:)
    real (kind=real_kind), pointer :: tensor2d (:,:,:,:,:)

    ! Set all geometric views (we use 1 tensor, 1 vector, and 1 scalar temps)
    call c_f_pointer(d_ptr,    tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(mdet_ptr, scalar2d, [np,np,    nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%d
      scalar2d(:,:,ie)     = elem(ie)%metdet
    enddo

    call c_f_pointer(dinv_ptr, tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(spmp_ptr, scalar2d, [np,np,    nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%dinv
      scalar2d(:,:,ie)     = elem(ie)%spheremp
    enddo

    call c_f_pointer(minv_ptr,  tensor2d, [np,np,2,2,nelemd])
    call c_f_pointer(rspmp_ptr, scalar2d, [np,np,    nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%metinv
      scalar2d(:,:,ie)     = elem(ie)%rspheremp
    enddo

    call c_f_pointer(sph2c_ptr, tensor2d, [np,np,3,2,nelemd])
    call c_f_pointer(mdet_ptr,  scalar2d, [np,np,    nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie) = elem(ie)%vec_sphere2cart
      scalar2d(:,:,ie)     = elem(ie)%metdet
    enddo

    call c_f_pointer(phis_ptr,     scalar2d, [np,np,    nelemd])
    call c_f_pointer(gradphis_ptr, vector2d, [np,np,2,  nelemd])
    call c_f_pointer(tVisc_ptr,    tensor2d, [np,np,2,2,nelemd])
    do ie=1,nelemd
      tensor2d(:,:,:,:,ie)      = elem(ie)%tensorVisc
      elem(ie)%derived%gradphis = vector2d(:,:,:,ie)
      elem(ie)%state%phis       = scalar2d(:,:,ie)
    enddo

    call c_f_pointer(fcor_ptr, scalar2d, [np, np, nelemd])
    do ie=1,nelemd
      scalar2d(:,:,ie) = elem(ie)%fcor
    enddo

  end subroutine init_geo_views_f90

  subroutine cleanup_f90 () bind(c)
    use edge_mod_base, only : FreeEdgeBuffer, edge_g
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


  subroutine initialize_reference_states_f90(phis_ptr, dp_ref_ptr, theta_ref_ptr, phi_ref_ptr) bind(c)
    use element_ops,    only: initialize_reference_states
    use dimensions_mod, only: nelemd, nlev, nlevp, np
    use theta_f2c_mod,  only : init_reference_states_c
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: phis_ptr, dp_ref_ptr, theta_ref_ptr, phi_ref_ptr
    !
    ! Locals
    !
    real (kind=real_kind), dimension(:,:,:),  pointer :: phis
    real (kind=real_kind), dimension(:,:,:,:),  pointer :: dp_ref, theta_ref, phi_ref
    integer :: ie

    call c_f_pointer(phis_ptr,      phis,      [np,np,      nelemd])
    call c_f_pointer(dp_ref_ptr,    dp_ref,    [np,np,nlev, nelemd])
    call c_f_pointer(theta_ref_ptr, theta_ref, [np,np,nlev, nelemd])
    call c_f_pointer(phi_ref_ptr,   phi_ref,   [np,np,nlevp,nelemd])

    ! Compute reference states
    do ie=1,nelemd
      call initialize_reference_states(hvcoord,             &
                                       phis(:,:,ie),        &
                                       dp_ref(:,:,:,ie),    &
                                       theta_ref(:,:,:,ie), &
                                       phi_ref(:,:,:,ie))
    end do

    ! Initialize reference states in C++
    call init_reference_states_c(dp_ref_ptr,theta_ref_ptr,phi_ref_ptr)

  end subroutine initialize_reference_states_f90

end module thetal_test_interface
