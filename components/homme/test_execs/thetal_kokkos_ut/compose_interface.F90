module compose_interface
  use kinds, only: real_kind
  use iso_c_binding, only: c_bool, c_int, c_double, c_ptr, c_loc
  use dimensions_mod, only: nlev, nlevp, np, nelemd, ne, qsize, qsize_d
  use geometry_interface_mod, only: par, elem
  implicit none

contains

  subroutine init_compose_f90(ne, hyai, hybi, hyam, hybm, ps0, dvv, mp, qsize_in, hv_q, &
       lim, cdr_check, is_sphere, nearest_point, halo, traj_nsubstep) bind(c)
    use hybvcoord_mod, only: set_layer_locations
    use thetal_test_interface, only: init_f90
    use theta_f2c_mod, only: init_elements_c
    use edge_mod_base, only: initEdgeBuffer, edge_g
    use control_mod, only: transport_alg, semi_lagrange_cdr_alg, semi_lagrange_cdr_check, &
         semi_lagrange_hv_q, limiter_option, nu_q, hypervis_subcycle_q, hypervis_order, &
         vert_remap_q_alg, qsplit, rsplit, dt_remap_factor, dt_tracer_factor, &
         theta_hydrostatic_mode, semi_lagrange_nearest_point_lev, semi_lagrange_halo, &
         semi_lagrange_trajectory_nsubstep
    use geometry_interface_mod, only: GridVertex
    use bndry_mod, only: sort_neighbor_buffer_mapping
    use reduction_mod, only: initreductionbuffer, red_sum, red_min, red_max
    use parallel_mod, only: global_shared_buf, nrepro_vars
    use compose_mod, only: compose_init, cedr_set_ie2gci, compose_set_null_bufs
    use sl_advection, only: sl_init1

    real (real_kind), intent(in) :: hyai(nlevp), hybi(nlevp), hyam(nlev), hybm(nlev)
    integer (c_int), value, intent(in) :: ne, qsize_in, hv_q, lim, halo, traj_nsubstep
    real (real_kind), value, intent(in) :: ps0
    real (real_kind), intent(out) :: dvv(np,np), mp(np,np)
    logical (c_bool), value, intent(in) :: cdr_check, is_sphere, nearest_point

    integer :: ie, edgesz

    if (.not. is_sphere) print *, "NOT IMPL'ED YET"

    transport_alg = 12
    semi_lagrange_cdr_alg = 3
    semi_lagrange_cdr_check = cdr_check
    qsize = qsize_in
    limiter_option = lim
    vert_remap_q_alg = 10
    qsplit = 1
    rsplit = 1
    dt_tracer_factor = -1
    dt_remap_factor = -1
    theta_hydrostatic_mode = .true.
    semi_lagrange_nearest_point_lev = -1
    if (nearest_point) semi_lagrange_nearest_point_lev = 100000
    semi_lagrange_halo = halo
    semi_lagrange_trajectory_nsubstep = traj_nsubstep

    hypervis_order = 2
    semi_lagrange_hv_q = hv_q
    nu_q = min(1.0e13_real_kind, 1.0e15_real_kind*(30.0_real_kind/ne)**3.2_real_kind)
    hypervis_subcycle_q = 6

    call init_f90(ne, hyai, hybi, hyam, hybm, dvv, mp, ps0)
    call init_elements_c(nelemd)

    edgesz = max((qsize+3)*nlev+2,6*nlev+1)
    call initEdgeBuffer(par, edge_g, elem, edgesz)

    call initReductionBuffer(red_sum,5)
    call initReductionBuffer(red_min,1)
    call initReductionBuffer(red_max,1)
    allocate(global_shared_buf(nelemd, nrepro_vars))
    
    call sort_neighbor_buffer_mapping(par, elem, 1, nelemd)
    call compose_init(par, elem, GridVertex, init_kokkos=.false.)
    do ie = 1, nelemd
       call cedr_set_ie2gci(ie, elem(ie)%vertex%number)
    end do
    call sl_init1(par, elem)
    call compose_set_null_bufs()
  end subroutine init_compose_f90

  subroutine init_geometry_f90() bind(c)
    use coordinate_systems_mod, only: cartesian3D_t
    use sl_advection, only: sphere2cart
    use theta_f2c_mod, only: init_elements_2d_c, init_geopotential_c

    real (real_kind), target, dimension(np,np)     :: elem_mp, elem_fcor, elem_spheremp, &
         elem_rspheremp, elem_metdet, elem_state_phis
    real (real_kind), target, dimension(np,np,2)   :: elem_gradphis
    real (real_kind), target, dimension(np,np,2,2) :: elem_D, elem_Dinv, elem_metinv, elem_tensorvisc
    real (real_kind), target, dimension(np,np,3,2) :: elem_vec_sph2cart
    type (c_ptr) :: elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr, elem_spheremp_ptr, &
         elem_rspheremp_ptr, elem_metdet_ptr, elem_metinv_ptr, elem_tensorvisc_ptr, &
         elem_vec_sph2cart_ptr, elem_state_phis_ptr, elem_gradphis_ptr

    type (cartesian3D_t) :: sphere_cart
    real (kind=real_kind) :: sphere_cart_vec(3,np,np), sphere_latlon_vec(2,np,np)

    integer :: ie, i, j

    elem_D_ptr            = c_loc(elem_D)
    elem_Dinv_ptr         = c_loc(elem_Dinv)
    elem_fcor_ptr         = c_loc(elem_fcor)
    elem_spheremp_ptr     = c_loc(elem_spheremp)
    elem_rspheremp_ptr    = c_loc(elem_rspheremp)
    elem_metdet_ptr       = c_loc(elem_metdet)
    elem_metinv_ptr       = c_loc(elem_metinv)
    elem_tensorvisc_ptr   = c_loc(elem_tensorvisc)
    elem_vec_sph2cart_ptr = c_loc(elem_vec_sph2cart)
    elem_state_phis_ptr   = c_loc(elem_state_phis)
    elem_gradphis_ptr     = c_loc(elem_gradphis)
    do ie = 1,nelemd
      elem_D            = elem(ie)%D
      elem_Dinv         = elem(ie)%Dinv
      elem_fcor         = elem(ie)%fcor
      elem_spheremp     = elem(ie)%spheremp
      elem_rspheremp    = elem(ie)%rspheremp
      elem_metdet       = elem(ie)%metdet
      elem_metinv       = elem(ie)%metinv
      elem_state_phis   = elem(ie)%state%phis
      elem_gradphis     = elem(ie)%derived%gradphis
      elem_tensorvisc   = elem(ie)%tensorVisc
      elem_vec_sph2cart = elem(ie)%vec_sphere2cart
      do j = 1,np
         do i = 1,np
            call sphere2cart(elem(ie)%spherep(i,j), sphere_cart)
            sphere_cart_vec(1,i,j) = sphere_cart%x
            sphere_cart_vec(2,i,j) = sphere_cart%y
            sphere_cart_vec(3,i,j) = sphere_cart%z
            sphere_latlon_vec(1,i,j) = elem(ie)%spherep(i,j)%lat
            sphere_latlon_vec(2,i,j) = elem(ie)%spherep(i,j)%lon
         end do
      end do
      call init_elements_2d_c(ie-1, elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr, &
           elem_spheremp_ptr, elem_rspheremp_ptr, elem_metdet_ptr, elem_metinv_ptr, &
           elem_tensorvisc_ptr, elem_vec_sph2cart_ptr, sphere_cart_vec, &
           sphere_latlon_vec)
      call init_geopotential_c(ie-1, elem_state_phis_ptr, elem_gradphis_ptr)
    enddo
  end subroutine init_geometry_f90

  subroutine cleanup_compose_f90() bind(c)
    use compose_mod, only: compose_finalize
    use thetal_test_interface, only: cleanup_f90

    call compose_finalize(finalize_kokkos=.false.)
    call cleanup_f90()
  end subroutine cleanup_compose_f90

  subroutine run_compose_standalone_test_f90(nmax_out, eval) bind(c)
    use thetal_test_interface, only: deriv, hvcoord
    use compose_test_mod, only: compose_test
    use domain_mod, only: domain1d_t
    use control_mod, only: transport_alg, statefreq
    use time_mod, only: nmax
    use thread_mod, only: hthreads, vthreads
    use dimensions_mod, only: nlev, qsize

    integer(c_int), intent(inout) :: nmax_out
    real(c_double), intent(out) :: eval((nlev+1)*qsize)

    type (domain1d_t), pointer :: dom_mt(:)
    real(real_kind) :: buf((nlev+1)*qsize)

    integer :: i

    hthreads = 1
    vthreads = 1
    allocate(dom_mt(0:0))
    dom_mt(0)%start = 1
    dom_mt(0)%end = nelemd
    transport_alg = 19
    if (nmax_out <= 1) then
       nmax = 7*ne
       nmax_out = nmax
    else
       nmax = nmax_out
    end if
    statefreq = 2*ne
    call compose_test(par, hvcoord, dom_mt, elem, buf)
    do i = 1,size(buf)
       eval(i) = buf(i)
    end do
    transport_alg = 12
    deallocate(dom_mt)
  end subroutine run_compose_standalone_test_f90

  subroutine run_trajectory_f90(t0, t1, independent_time_steps, dep, dprecon) bind(c)
    use time_mod, only: timelevel_t, timelevel_init_default
    use control_mod, only: qsplit, geometry
    use physical_constants, only: Sx, Sy, Lx, Ly
    use hybrid_mod, only: hybrid_t, hybrid_create
    use thetal_test_interface, only: deriv, hvcoord
    use compose_test_mod, only: compose_stt_init, compose_stt_fill_v, compose_stt_clear, &
         compose_stt_fill_uniform_density
    use sl_advection, only: calc_trajectory, dep_points_all

    real(c_double), value, intent(in) :: t0, t1
    logical(c_bool), value, intent(in) :: independent_time_steps
    real(c_double), intent(out) :: dep(3,np,np,nlev,nelemd), dprecon(np,np,nlev,nelemd)

    type (timelevel_t) :: tl
    type (hybrid_t) :: hybrid
    real(real_kind) :: dt
    integer :: ie, i, j, k, d, testno, geometry_type
    logical :: its

    call timelevel_init_default(tl)
    testno = 0
    geometry_type = 0 ! sphere
    if (trim(geometry) == "plane") geometry_type = 1
    call compose_stt_init(np, nlev, qsize, qsize_d, nelemd, testno, &
         geometry_type, Sx, Sy, Lx, Ly)
    do ie = 1, nelemd
       call compose_stt_fill_uniform_density(ie, tl%np1, elem(ie)%state%dp3d, &
            elem(ie)%derived%dp)
       call compose_stt_fill_v(ie, elem(ie)%spherep, t0, &
            elem(ie)%derived%vstar)
       call compose_stt_fill_v(ie, elem(ie)%spherep, t1, &
            elem(ie)%state%v(:,:,:,:,tl%np1))
    end do
    hybrid = hybrid_create(par, 0, 1)
    dt = t1 - t0
    its = independent_time_steps
    call calc_trajectory(elem, deriv, hvcoord, hybrid, dt, tl, its, 1, nelemd)

    do ie = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                do d = 1, 3
                   dep(d,i,j,k,ie) = dep_points_all(d,i,j,k,ie)
                end do
                dprecon(i,j,k,ie) = elem(ie)%derived%divdp(i,j,k)
             end do
          end do
       end do
    end do
  end subroutine run_trajectory_f90

  subroutine run_sl_vertical_remap_bfb_f90(diagnostic) bind(c)
    use time_mod, only: timelevel_t, timelevel_init_default, timelevel_qdp
    use sl_advection, only: sl_vertically_remap_tracers
    use hybrid_mod, only: hybrid_t, hybrid_create
    use control_mod, only: qsplit, dt_tracer_factor
    
    real (c_double), intent(out) :: diagnostic

    real (real_kind), parameter :: c1 = 0.491827d0, c2 = 0.2432109d0, c3 = 0.1234567d0, c4 = 0.0832d0
    
    real (real_kind) :: unused
    type (timelevel_t) :: tl
    type (hybrid_t) :: hybrid
    integer :: ie, i, j, k, krev, q, nq, n0_qdp, np1_qdp

    hybrid = hybrid_create(par, 0, 1)
    call timelevel_init_default(tl)
    tl%np1 = 1
    tl%nstep = 42
    call timelevel_qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)
    
    nq = size(elem(1)%state%q,4)
    unused = 0

    ! manufactured values
    do ie = 1,nelemd
       do q = 1,nq
          do k = 1,nlev
             krev = nlev - k + 1
             do j = 1,np
                do i = 1,np
                   if (q == 1) then
                      elem(ie)%state%dp3d(i,j,k,tl%np1) = k + c1*i + c2*j*i + c3*(j*j + j*k)
                      elem(ie)%derived%divdp(i,j,k) = krev + c1*i + c2*j*i + c3*(j*j + j*krev)
                   end if
                   elem(ie)%state%qdp(i,j,k,q,np1_qdp) = ie + c1*i + c2*j*i + c3*(j*j + j*k) + c4*q
                end do
             end do
          end do
       end do
    end do
    
    call sl_vertically_remap_tracers(hybrid, elem, 1, nelemd, tl, unused)

    diagnostic = 0
    do ie = 1,nelemd
       do q = 1,nq
          do k = 1,nlev
             krev = nlev - k + 1
             do j = 1,np
                do i = 1,np
                   diagnostic = diagnostic + elem(ie)%state%q(i,j,k,q)
                end do
             end do
          end do
       end do
    end do
  end subroutine run_sl_vertical_remap_bfb_f90
  
end module compose_interface
