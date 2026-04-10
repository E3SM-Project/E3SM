module eamxx_repro_sum_mod

  implicit none

contains

  subroutine eamxx_repro_sum(send, recv, nlocal, nfld, comm) bind(c)
    use iso_c_binding, only: c_int, c_double

    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc

    integer(kind=c_int), value, intent(in) :: nlocal, nfld, comm
    real(kind=c_double), intent(in) :: send(nlocal,nfld)
    real(kind=c_double), intent(out) :: recv(nfld)

    call repro_sum(send, recv, nlocal, nlocal, nfld, commid=comm)
  end subroutine eamxx_repro_sum

end module eamxx_repro_sum_mod
