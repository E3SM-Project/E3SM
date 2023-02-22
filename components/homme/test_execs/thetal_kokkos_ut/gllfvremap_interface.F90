module physgrid_interface
  use kinds, only: real_kind
  use iso_c_binding, only: c_bool, c_int, c_double, c_ptr, c_loc
  use dimensions_mod, only: nlev, nlevp, np, nelemd, ne, qsize, qsize_d
  use geometry_interface_mod, only: par, elem
  implicit none
  
contains

  subroutine init_gllfvremap_f90(ne, hyai, hybi, hyam, hybm, ps0, dvv, mp, qsize_in, &
       is_sphere) bind(c)
    use iso_c_binding, only: c_double
    use hybvcoord_mod, only: set_layer_locations
    use thetal_test_interface, only: init_f90
    use theta_f2c_mod, only: init_elements_c
    use edge_mod_base, only: initEdgeBuffer, edge_g, initEdgeSBuffer
    use prim_advection_base, only: edgeAdvQminmax
    use bndry_mod, only: sort_neighbor_buffer_mapping
    use reduction_mod, only: initreductionbuffer, red_sum, red_min, red_max
    use parallel_mod, only: global_shared_buf, nrepro_vars
    use control_mod, only: use_moisture, theta_hydrostatic_mode

    real (c_double), intent(in) :: hyai(nlevp), hybi(nlevp), hyam(nlev), hybm(nlev)
    integer (c_int), value, intent(in) :: ne, qsize_in
    real (c_double), value, intent(in) :: ps0
    real (c_double), intent(out) :: dvv(np,np), mp(np,np)
    logical (c_bool), value, intent(in) :: is_sphere

    integer :: edgesz

    if (.not. is_sphere) print *, "NOT IMPL'ED YET"

    qsize = qsize_in
    use_moisture = .true.
    theta_hydrostatic_mode = .false.

    call init_f90(ne, hyai, hybi, hyam, hybm, dvv, mp, ps0)
    call init_elements_c(nelemd)

    edgesz = max((qsize+3)*nlev+2,6*nlev+1)
    call initEdgeBuffer(par, edge_g, elem, edgesz)
    call initEdgeSBuffer(par, edgeAdvQminmax, elem, qsize*nlev*2)

    call initReductionBuffer(red_sum,5)
    call initReductionBuffer(red_min,1)
    call initReductionBuffer(red_max,1)
    allocate(global_shared_buf(nelemd, nrepro_vars))
  end subroutine init_gllfvremap_f90

  subroutine run_gfr_test(nerr) bind(c)
    use thetal_test_interface, only: deriv, hvcoord
    use domain_mod, only: domain1d_t
    use hybrid_mod, only: hybrid_t, hybrid_create
    use gllfvremap_mod, only: gfr_test

    integer (c_int), intent(out) :: nerr

    type (hybrid_t) :: hybrid
    type (domain1d_t) :: dom_mt(0:0)

    dom_mt(0)%start = 1
    dom_mt(0)%end = nelemd
    hybrid = hybrid_create(par, 0, 1)

    nerr = gfr_test(hybrid, dom_mt, hvcoord, deriv, elem)
  end subroutine run_gfr_test

  subroutine gfr_init_f90(nf, thm) bind(c)
    use gllfvremap_mod, only: gfr_init
    use control_mod, only: theta_hydrostatic_mode
    
    integer (c_int), value, intent(in) :: nf
    logical (c_bool), value, intent(in) :: thm

    theta_hydrostatic_mode = thm
    call gfr_init(par, elem, nf, 2, .false.)
  end subroutine gfr_init_f90
  
  subroutine gfr_finish_f90() bind(c)
    use gllfvremap_mod, only: gfr_finish

    call gfr_finish()
  end subroutine gfr_finish_f90

  subroutine run_gfr_check_api(thm, nerr) bind(c)
    use thetal_test_interface, only: hvcoord
    use hybrid_mod, only: hybrid_t, hybrid_create
    use control_mod, only: theta_hydrostatic_mode
    use gllfvremap_util_mod, only: gfr_check_api

    logical (c_bool), value, intent(in) :: thm
    integer (c_int), intent(out) :: nerr

    type (hybrid_t) :: hybrid
    logical :: thm_save
    
    hybrid = hybrid_create(par, 0, 1)
    thm_save = theta_hydrostatic_mode
    theta_hydrostatic_mode = thm
    nerr = gfr_check_api(hybrid, 1, nelemd, hvcoord, elem)
    theta_hydrostatic_mode = thm_save
  end subroutine run_gfr_check_api

  subroutine limiter1_clip_and_sum_f90(n, spheremp, qmin, qmax, dp, q) bind(c)
    use gllfvremap_mod, only: limiter1_clip_and_sum

    integer (c_int), value, intent(in) :: n
    real (c_double), intent(in) :: spheremp(n*n), dp(n*n)
    real (c_double), intent(inout) :: qmin, qmax, q(n*n)

    call limiter1_clip_and_sum(n*n, spheremp, qmin, qmax, dp, q)
  end subroutine limiter1_clip_and_sum_f90

  subroutine calc_dp_fv_f90(nf, ps, dp_fv) bind(c)
    use thetal_test_interface, only: hvcoord
    use gllfvremap_mod, only: calc_dp_fv

    integer (c_int), value, intent(in) :: nf
    real (c_double), intent(in) :: ps(nf*nf)
    real (c_double), intent(out) :: dp_fv(nf*nf,nlev)

    call calc_dp_fv(nf, hvcoord, ps, dp_fv)
  end subroutine calc_dp_fv_f90

  subroutine get_temperature_f90(ie, nt, thm, T) bind(c)
    use thetal_test_interface, only: hvcoord
    use element_ops, only: get_temperature
    use control_mod, only: theta_hydrostatic_mode
    
    integer (c_int), value, intent(in) :: ie, nt
    real (c_double), intent(out) :: T(np,np,nlev)
    logical (c_bool), value, intent(in) :: thm

    logical :: thm_save

    thm_save = theta_hydrostatic_mode
    theta_hydrostatic_mode = thm
    call get_temperature(elem(ie), T(:,:,:), hvcoord, nt)
    theta_hydrostatic_mode = thm_save
  end subroutine get_temperature_f90

  subroutine init_dyn_data_f90(nk, nkp, nq, ps, phis, dp3d, vthdp, uv, omega, q, phinh_i) bind(c)
    use element_state, only: timelevels

    integer (c_int), value, intent(in) :: nk, nkp, nq
    real (c_double), intent(in) :: ps(np,np,timelevels,nelemd), phis(np,np,nelemd), &
         dp3d(nk,np,np,timelevels,nelemd), vthdp(nk,np,np,timelevels,nelemd), &
         uv(nk,np,np,2,timelevels,nelemd), omega(nk,np,np,nelemd), q(nk,np,np,nq,nelemd), &
         phinh_i(nkp,np,np,timelevels,nelemd)
    
    integer :: ie, k, tl, iq

    do ie = 1,nelemd
       elem(ie)%state%phis = phis(:,:,ie)
       elem(ie)%state%ps_v = ps(:,:,:,ie)
       do k = 1,nlev
          do tl = 1,timelevels
             elem(ie)%state%dp3d(:,:,k,tl) = dp3d(k,:,:,tl,ie)
             elem(ie)%state%vtheta_dp(:,:,k,tl) = vthdp(k,:,:,tl,ie)
             elem(ie)%state%v(:,:,:,k,tl) = uv(k,:,:,:,tl,ie)
             elem(ie)%state%phinh_i(:,:,k,tl) = phinh_i(k,:,:,tl,ie)
             if (k == nlev) elem(ie)%state%phinh_i(:,:,k+1,tl) = phinh_i(k+1,:,:,tl,ie)
          end do
          elem(ie)%derived%omega_p(:,:,k) = omega(k,:,:,ie)
          do iq = 1,nq
             elem(ie)%state%q(:,:,k,iq) = q(k,:,:,iq,ie)
          end do
       end do
    end do
  end subroutine init_dyn_data_f90

  subroutine gfr_dyn_to_fv_phys_f90(nf, nt, ps, phis, T, uv, omega_p, q) bind(c)
    use thetal_test_interface, only: hvcoord
    use hybrid_mod, only: hybrid_t, hybrid_create
    use gllfvremap_mod, only: gfr_dyn_to_fv_phys, gfr_get_nphys
    
    integer (c_int), value, intent(in) :: nf, nt
    real (c_double), intent(out) :: ps(nf*nf,nelemd), phis(nf*nf,nelemd), &
         T(nf*nf,nlev,nelemd), uv(nf*nf,2,nlev,nelemd), omega_p(nf*nf,nlev,nelemd), &
         q(nf*nf,nlev,qsize,nelemd)

    type (hybrid_t) :: hybrid

    if (nf /= gfr_get_nphys()) print *, 'ERROR: nf vs gfr%nphys:', nf, gfr_get_nphys()

    hybrid = hybrid_create(par, 0, 1)
    call gfr_dyn_to_fv_phys(hybrid, nt, hvcoord, elem, 1, nelemd, &
         ps, phis, T, uv, omega_p, q)
  end subroutine gfr_dyn_to_fv_phys_f90

  subroutine gfr_fv_phys_to_dyn_f90(nf, nt, T, uv, q) bind(c)
    use thetal_test_interface, only: hvcoord
    use hybrid_mod, only: hybrid_t, hybrid_create
    use gllfvremap_mod, only: gfr_fv_phys_to_dyn, gfr_get_nphys, gfr_f2g_dss

    integer (c_int), value, intent(in) :: nf, nt
    real (c_double), intent(in) :: T(nf*nf,nlev,nelemd), uv(nf*nf,2,nlev,nelemd), &
         q(nf*nf,nlev,qsize,nelemd)

    type (hybrid_t) :: hybrid

    if (nf /= gfr_get_nphys()) print *, 'ERROR: nf vs gfr%nphys:', nf, gfr_get_nphys()

    hybrid = hybrid_create(par, 0, 1)
    call gfr_fv_phys_to_dyn(hybrid, nt, hvcoord, elem, 1, nelemd, T, uv, q)
    call gfr_f2g_dss(hybrid, elem, 1, nelemd)
  end subroutine gfr_fv_phys_to_dyn_f90

  function equal(a, b) result(e)
    real (c_double), intent(in) :: a, b
    logical :: e
    real (c_double), parameter :: eps = epsilon(1.0_c_double)
#ifdef HOMMEXX_BFB_TESTING
    e = a == b
#else
    e = abs(a - b) <= 1d10*eps*abs(a)
#endif
  end function equal

  subroutine cmp_dyn_data_f90(nk, nq, ft, fm, q, fq, nerr) bind(c)
    integer (c_int), value, intent(in) :: nk, nq
    real (c_double), intent(in) :: ft(nk,np,np,nelemd), &
         fm(nk,np,np,3,nelemd), q(nk,np,np,nq,nelemd), fq(nk,np,np,nq,nelemd)
    integer (c_int), intent(out) :: nerr

    integer, parameter :: outmax = 20

    integer :: ie, i, j, k, iq, d

    nerr = 0
    do ie = 1,nelemd
       do i = 1,np
          do j = 1,np
             do k = 1,nlev
                if (.not. equal(elem(ie)%derived%FT(j,i,k), ft(k,j,i,ie))) then
                   nerr = nerr+1
                   if (nerr < outmax) then
                      print '(a,i4,i2,i2,i3,es23.15,es23.15)', &
                           'ft: ie,i,j,k',ie,i,j,k, &
                           elem(ie)%derived%FT(j,i,k), ft(k,j,i,ie)
                   end if
                end if
                do d = 1,2
                   if (.not. equal(elem(ie)%derived%FM(j,i,d,k), fm(k,j,i,d,ie))) then
                      nerr = nerr+1
                      if (nerr < outmax) then
                         print '(a,i4,i2,i2,i2,i3,es23.15,es23.15)', &
                              'fm: ie,d,i,j,k',ie,d,i,j,k, &
                              elem(ie)%derived%FM(j,i,d,k), fm(k,j,i,d,ie)
                      end if
                   end if
                end do
                do iq = 1,qsize
                   if (.not. equal(elem(ie)%state%Q(j,i,k,iq), q(k,j,i,iq,ie))) then
                      nerr = nerr+1
                      if (nerr < outmax) then
                         print '(a,i4,i3,i2,i2,i3,es23.15,es23.15)', &
                              'q: ie,iq,i,j,k',ie,iq,i,j,k, &
                              elem(ie)%state%Q(j,i,k,iq), q(k,j,i,iq,ie)
                      end if
                   end if
                   if (.not. equal(elem(ie)%derived%FQ(j,i,k,iq), fq(k,j,i,iq,ie))) then
                      nerr = nerr+1
                      if (nerr < outmax) then
                         print '(a,i4,i3,i2,i2,i3,es23.15,es23.15)', &
                              'fq: ie,iq,i,j,k',ie,iq,i,j,k, &
                              elem(ie)%derived%FQ(j,i,k,iq), fq(k,j,i,iq,ie)
                      end if
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine cmp_dyn_data_f90
end module physgrid_interface
