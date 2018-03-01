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
      rxt_rates(:ncol,:,     1) = rxt_rates(:ncol,:,     1)*sol(:ncol,:,     1)                                                ! rate_const*CH4
      rxt_rates(:ncol,:,     2) = rxt_rates(:ncol,:,     2)*sol(:ncol,:,     2)                                                ! rate_const*N2O
      rxt_rates(:ncol,:,     3) = rxt_rates(:ncol,:,     3)*sol(:ncol,:,     3)                                                ! rate_const*CFC11
      rxt_rates(:ncol,:,     4) = rxt_rates(:ncol,:,     4)*sol(:ncol,:,     4)                                                ! rate_const*CFC12
      rxt_rates(:ncol,:,     5) = rxt_rates(:ncol,:,     5)*sol(:ncol,:,     5)                                                ! rate_const*H2O
  end subroutine set_rates
end module mo_rxt_rates_conv
