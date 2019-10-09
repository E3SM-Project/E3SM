module mo_rxt_rates_conv
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  private
  public :: set_rates
contains
   subroutine set_rates( rxt_rates, sol, ncol )
      real(r8), intent(inout) :: rxt_rates(:,:,:)
      real(r8), intent(in) :: sol(:,:,:)
      integer, intent(in) :: ncol
      rxt_rates(:ncol,:,     1) = rxt_rates(:ncol,:,     1)*sol(:ncol,:,     1)                                                ! rate_const*H2O2
                                                                                                                               ! rate_const
      rxt_rates(:ncol,:,     3) = rxt_rates(:ncol,:,     3)*sol(:ncol,:,     1)                                                ! rate_const*OH*H2O2
      rxt_rates(:ncol,:,     4) = rxt_rates(:ncol,:,     4)*sol(:ncol,:,     3)                                                ! rate_const*OH*SO2
      rxt_rates(:ncol,:,     5) = rxt_rates(:ncol,:,     5)*sol(:ncol,:,     4)                                                ! rate_const*OH*DMS
      rxt_rates(:ncol,:,     6) = rxt_rates(:ncol,:,     6)*sol(:ncol,:,     4)                                                ! rate_const*OH*DMS
      rxt_rates(:ncol,:,     7) = rxt_rates(:ncol,:,     7)*sol(:ncol,:,     4)                                                ! rate_const*NO3*DMS
      rxt_rates(:ncol,:,     8) = rxt_rates(:ncol,:,     8)*sol(:ncol,:,     5)                                                ! rate_const*OH*NH3
  end subroutine set_rates
end module mo_rxt_rates_conv
