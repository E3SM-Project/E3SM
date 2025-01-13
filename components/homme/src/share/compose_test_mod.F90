#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module compose_test_mod

  implicit none

interface

#ifdef HOMME_ENABLE_COMPOSE
   subroutine compose_unittest() bind(c)
   end subroutine compose_unittest

   subroutine compose_stt_init(np, nlev, qsize, qsize_d, nelemd, testno, &
        geometry, Sx, Sy, Lx, Ly) bind(c)
     use iso_c_binding, only: c_int, c_double
     integer (kind=c_int), value, intent(in) :: np, nlev, qsize, qsize_d, nelemd, &
          testno, geometry
     real (kind=c_double), value, intent(in) :: Sx, Sy, Lx, Ly
   end subroutine compose_stt_init

   subroutine compose_stt_fill_uniform_density(ie, np1, dp3d, dp) bind(c)
     use iso_c_binding, only: c_int, c_double
     use dimensions_mod, only: nlev, np
     use element_state, only: timelevels
     integer (kind=c_int), value, intent(in) :: ie, np1
     real (kind=c_double), intent(in) :: dp3d(np,np,nlev,timelevels), &
          dp(np,np,nlev)
   end subroutine compose_stt_fill_uniform_density

   subroutine compose_stt_fill_ics(ie, p_elem, dp, n0_qdp, qdp) bind(c)
     use iso_c_binding, only: c_int, c_double
     use dimensions_mod, only: nlev, np, qsize_d
     use coordinate_systems_mod, only: spherical_polar_t
     integer (kind=c_int), value, intent(in) :: ie, n0_qdp
     type(spherical_polar_t), intent(in) :: p_elem(np,np)
     real (kind=c_double), intent(in) :: dp(np,np,nlev), &
          qdp(np,np,nlev,qsize_d,2)
   end subroutine compose_stt_fill_ics

   subroutine compose_stt_fill_v(ie, p_elem, t, v) bind(c)
     use iso_c_binding, only: c_int, c_double
     use dimensions_mod, only: nlev, np
     use coordinate_systems_mod, only: spherical_polar_t
     integer (kind=c_int), value, intent(in) :: ie
     type(spherical_polar_t), intent(in) :: p_elem(np,np)
     real (kind=c_double), intent(in) :: t, v(np,np,2,nlev)
   end subroutine compose_stt_fill_v

   subroutine compose_stt_begin_record() bind(c)
   end subroutine compose_stt_begin_record

   subroutine compose_stt_record_q(ie, p_elem, spheremp, np1, dp3d, n0_qdp, qdp) bind(c)
     use iso_c_binding, only: c_int, c_double
     use dimensions_mod, only: nlev, np, qsize_d
     use element_state, only: timelevels
     use coordinate_systems_mod, only: spherical_polar_t
     integer (kind=c_int), value, intent(in) :: ie, np1, n0_qdp
     type(spherical_polar_t), intent(in) :: p_elem(np,np)
     real (kind=c_double), intent(in) :: spheremp(np,np), &
          dp3d(np,np,nlev,timelevels), qdp(np,np,nlev,qsize_d,2)
   end subroutine compose_stt_record_q

   subroutine compose_stt_clear() bind(c)
   end subroutine compose_stt_clear

   subroutine compose_stt_finish(comm, root, rank, eval) bind(c)
     use iso_c_binding, only: c_int, c_double
     use dimensions_mod, only: nlev, qsize
     integer (kind=c_int), value, intent(in) :: comm, root, rank
     real (kind=c_double), intent(out) :: eval((nlev+1)*qsize)
   end subroutine compose_stt_finish
#endif

end interface

contains

  ! For comprehensive testing.
  subroutine compose_test(par, hvcoord, dom_mt, elem, eval)
    use kinds, only: real_kind
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use element_mod, only: element_t
    use thread_mod, only: hthreads, vthreads, omp_set_num_threads, omp_get_thread_num
    use hybrid_mod, only: hybrid_t, hybrid_create
    use derivative_mod, only: derivative_t, derivinit
    use hybvcoord_mod, only: hvcoord_t
    use time_mod, only: nEndStep
    use control_mod, only: transport_alg
    use sl_advection, only: sl_unittest
#ifdef HOMME_ENABLE_COMPOSE
    use compose_mod, only: cedr_unittest
#endif

    type (parallel_t), intent(in) :: par
    type (hvcoord_t) , intent(in) :: hvcoord
    type (domain1d_t), pointer, intent(in) :: dom_mt(:)
    type (element_t), intent(inout) :: elem(:)
    real (real_kind), optional, intent(out) :: eval(:)

    type (hybrid_t) :: hybrid
    type (derivative_t) :: deriv
    integer :: ithr, nets, nete, nerr

#ifdef HOMME_ENABLE_COMPOSE
    if (transport_alg == 19) then
       nEndStep = -1
    else
       return
    end if

    call derivinit(deriv)

    if (par%masterproc) print *, '~*~ Comprehensively test COMPOSE ~*~'

    ! 1. Unit tests.
    call compose_unittest()
    call sl_unittest(par, hvcoord)
    nerr = 0
    call cedr_unittest(par%comm, nerr)
    if (nerr /= 0) print *, 'cedr_unittest returned', nerr

#if (defined HORIZ_OPENMP)
    !$omp parallel num_threads(hthreads), default(SHARED), private(ithr,nets,nete,hybrid)
    call omp_set_num_threads(vthreads)
#endif
    ithr = omp_get_thread_num()
    hybrid = hybrid_create(par, ithr, hthreads)
    nets = dom_mt(ithr)%start
    nete = dom_mt(ithr)%end

    ! 2. Standalone tracer advection test, useful as a basic but comprehensive
    ! correctness test and also as part of a convergence test.
    call compose_stt(hybrid, dom_mt, nets, nete, hvcoord, deriv, elem, eval)
#if (defined HORIZ_OPENMP)
    !$omp end parallel
#endif
#endif
  end subroutine compose_test

  subroutine print_software_statistics(hybrid, nets, nete)
    use hybrid_mod, only: hybrid_t
    use kinds, only: real_kind, iulog
    use parallel_mod, only: global_shared_buf, global_shared_sum
    use global_norms_mod, only: wrap_repro_sum
    use reduction_mod, only: parallelmax, parallelmin

    integer :: GPTLget_memusage

    type(hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nets, nete
    integer :: ok, size, rss_int, share, text, datastack, ie
    real(kind=real_kind) :: rss, rss_min, rss_max, rss_mean

    ok = GPTLget_memusage(size, rss_int, share, text, datastack)
    rss = rss_int
    rss_min = parallelmin(rss, hybrid)
    rss_max = parallelmax(rss, hybrid)
    do ie = nets, nete
       global_shared_buf(ie,1) = 0
    end do
    if (nets == 1) then
       global_shared_buf(1,1) = rss
    end if
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
    rss_mean = global_shared_sum(1) / hybrid%par%nprocs
    if (hybrid%par%masterproc) then
       write(iulog,'(a10,3(f14.2))') "rss   = ", rss_min, rss_max, rss_mean
    end if
  end subroutine print_software_statistics

  subroutine compose_stt(hybrid, dom_mt, nets, nete, hvcoord, deriv, elem, eval)
    use iso_c_binding, only: c_loc
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybrid_mod, only: hybrid_t
    use element_mod, only: element_t
    use time_mod, only: timelevel_t, timelevel_init_default, timelevel_qdp
    use kinds, only: real_kind
    use thread_mod, only: hthreads, vthreads, omp_set_num_threads, omp_get_thread_num
    use derivative_mod, only: derivative_t, derivinit
    use dimensions_mod, only: ne, np, nlev, qsize, qsize_d, nelemd
    use coordinate_systems_mod, only: spherical_polar_t
    use control_mod, only: qsplit, statefreq, se_fv_phys_remap_alg, geometry
    use physical_constants, only: Sx, Sy, Lx, Ly
    use time_mod, only: nmax
    use hybvcoord_mod, only: hvcoord_t
    use perf_mod
    use sl_advection
    use gllfvremap_mod
    use gllfvremap_util_mod
#ifdef HOMME_ENABLE_COMPOSE
    use compose_mod, only: compose_h2d, compose_d2h
#endif

    type (hybrid_t), intent(in) :: hybrid
    type (domain1d_t), pointer, intent(in) :: dom_mt(:)
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(inout) :: elem(:)
    type (hvcoord_t) , intent(in) :: hvcoord
    integer, intent(in) :: nets, nete
    real (real_kind), optional, intent(out) :: eval(:)

    real (kind=real_kind), parameter :: twelve_days = 3600.d0 * 24 * 12

    type (timelevel_t) :: tl
    integer :: nsteps, n0_qdp, np1_qdp, ie, i, j, geometry_type, nerr
    real (kind=real_kind) :: dt, tprev, t, unused((nlev+1)*qsize)

    nerr = 0
    
    if (se_fv_phys_remap_alg == -1) then
       nerr = nerr + gfr_test(hybrid, dom_mt, hvcoord, deriv, elem)
       nerr = nerr + gfr_check_api(hybrid, nets, nete, hvcoord, elem)
       return
    end if

#ifdef HOMME_ENABLE_COMPOSE  
    call t_startf('compose_stt')
    geometry_type = 0 ! sphere
    if (trim(geometry) == "plane") geometry_type = 1
    ! Set up time stepping and initialize q and density.
    call timelevel_init_default(tl)
    call timelevel_qdp(tl, qsplit, np1_qdp, n0_qdp)
    call compose_stt_init(np, nlev, qsize, qsize_d, nelemd, 0, &
         geometry_type, Sx, Sy, Lx, Ly)
    do ie = nets, nete
       call compose_stt_fill_uniform_density(ie, tl%np1, elem(ie)%state%dp3d, &
            elem(ie)%derived%dp)
       call compose_stt_fill_ics(ie, elem(ie)%spherep, elem(ie)%derived%dp, &
            n0_qdp, elem(ie)%state%qdp)
    end do
    ! Time step.
    !   For now, support only the nondivergent flow field. Supporting the
    ! divergent flow field will require providing a useful density.
    !   6*ne gives the large time step, which is said to be CN ~5.5. The
    ! corresponding qsplit is 15.
    !   nsteps = nint(6*ne*(15.d0/qsplit))
    nsteps = nmax
    if (hybrid%par%masterproc .and. hybrid%masterthread) then
       print '(a,i6,a,i5)', 'COMPOSE> nsteps ', nsteps, ' ne ', ne
    end if
    dt = twelve_days / nsteps
    call t_barrierf('compose_stt_step_start_barrier', hybrid%par%comm)
    call t_startf('compose_stt_step')
    do i = 1, nsteps
       compose_h2d = i == 1
       compose_d2h = i == 1 .or. i == nsteps
       tprev = dt*(i-1)
       t = dt*i
       do ie = nets, nete
          call compose_stt_fill_v(ie, elem(ie)%spherep, tprev, &
               elem(ie)%derived%vstar)
          call compose_stt_fill_v(ie, elem(ie)%spherep, t, &
               elem(ie)%state%v(:,:,:,:,tl%np1))
       end do
       tl%nstep = tl%nstep + qsplit
       call prim_advec_tracers_remap_ale(elem, deriv, hvcoord, hybrid, dt, tl, nets, nete)
       if (mod(i,statefreq) == 0) then
          call print_software_statistics(hybrid, nets, nete)
       end if
    end do
    call t_barrierf('compose_stt_step_stop_barrier', hybrid%par%comm)
    call t_stopf('compose_stt_step')
    ! Record final q values.
    call compose_stt_begin_record()
    call timelevel_qdp(tl, qsplit, n0_qdp, np1_qdp)
    do ie = nets, nete
       call compose_stt_record_q(ie, elem(ie)%spherep, elem(ie)%spheremp, &
            tl%np1, elem(ie)%state%dp3d, np1_qdp, elem(ie)%state%qdp)
    end do
    ! Do the global reductions, print diagnostic information, and clean up.
    if (present(eval)) then
       call compose_stt_finish(hybrid%par%comm, hybrid%par%root, hybrid%par%rank, eval)
    else
       call compose_stt_finish(hybrid%par%comm, hybrid%par%root, hybrid%par%rank, unused)
    end if
    call t_stopf('compose_stt')
#endif
  end subroutine compose_stt

end module compose_test_mod
