module compose_test_mod

interface
   subroutine compose_unittest()
   end subroutine compose_unittest
end interface

contains

  ! For comprehensive testing.
  subroutine compose_test(par, dom_mt, elem, tl)
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use element_mod, only: element_t
    use time_mod, only: timelevel_t

    use thread_mod, only: hthreads, vthreads, omp_set_num_threads, omp_get_thread_num
    use hybrid_mod, only: hybrid_t, hybrid_create
    use sl_advection

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (element_t), intent(in) :: elem(:)
    type (timelevel_t), intent(in) :: tl

    type (hybrid_t) :: hybrid
    integer :: ithr, nets, nete

#if (defined HORIZ_OPENMP)
    !$omp parallel num_threads(hthreads), default(SHARED), private(ithr,nets,nete,hybrid)
    call omp_set_num_threads(vthreads)
#endif
    ithr = omp_get_thread_num()
    hybrid = hybrid_create(par, ithr, hthreads)
    nets = dom_mt(ithr)%start
    nete = dom_mt(ithr)%end

    if (par%masterproc) print *, '~*~ Comprehensively test Compose ~*~'
    call compose_unittest()

#if (defined HORIZ_OPENMP)
    !$omp end parallel
#endif
  end subroutine compose_test

end module compose_test_mod
