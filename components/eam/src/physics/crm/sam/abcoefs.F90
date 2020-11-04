
module abcoefs_mod
  implicit none

contains

  subroutine abcoefs(ncrms)
    !      coefficients for the Adams-Bashforth scheme
    use grid
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) alpha, beta

    if(nstep.ge.3.and.nadams.eq.3.or.nrestart.eq.2) then
      alpha = dt3(nb) / dt3(na)
      beta = dt3(nc) / dt3(na)
      ct = (2.D0+3.D0* alpha) / (6.D0* (alpha + beta) * beta)
      bt = -(1.D0+2.D0*(alpha + beta) * ct)/(2.D0 * alpha)
      at = 1.D0 - bt - ct
    else if(nstep.ge.2) then
      at = 3.D0/2.D0
      bt = -1.D0/2.D0
      ct = 0.D0
    else
      at = 1.D0
      bt = 0.D0
      ct = 0.D0
    end if

  end subroutine abcoefs

end module abcoefs_mod
