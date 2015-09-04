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
      rxt_rates(:ncol,:,     1) = rxt_rates(:ncol,:,     1)*sol(:ncol,:,     1)                                                ! rate_const*O3
      rxt_rates(:ncol,:,     2) = rxt_rates(:ncol,:,     2)*sol(:ncol,:,     4)                                                ! rate_const*H2O2
      rxt_rates(:ncol,:,     3) = rxt_rates(:ncol,:,     3)*sol(:ncol,:,     6)                                                ! rate_const*NO2
      rxt_rates(:ncol,:,     4) = rxt_rates(:ncol,:,     4)*sol(:ncol,:,     9)                                                ! rate_const*CH2O
      rxt_rates(:ncol,:,     5) = rxt_rates(:ncol,:,     5)*sol(:ncol,:,     9)                                                ! rate_const*CH2O
      rxt_rates(:ncol,:,     6) = rxt_rates(:ncol,:,     6)*sol(:ncol,:,    11)                                                ! rate_const*CH3OOH
      rxt_rates(:ncol,:,     7) = rxt_rates(:ncol,:,     7)*sol(:ncol,:,     1)*sol(:ncol,:,     2)                            ! rate_const*O3*OH
      rxt_rates(:ncol,:,     8) = rxt_rates(:ncol,:,     8)*sol(:ncol,:,     3)*sol(:ncol,:,     1)                            ! rate_const*HO2*O3
      rxt_rates(:ncol,:,     9) = rxt_rates(:ncol,:,     9)*sol(:ncol,:,     3)*sol(:ncol,:,     2)                            ! rate_const*HO2*OH
      rxt_rates(:ncol,:,    10) = rxt_rates(:ncol,:,    10)*sol(:ncol,:,     3)*sol(:ncol,:,     3)                            ! rate_const*HO2*HO2
      rxt_rates(:ncol,:,    11) = rxt_rates(:ncol,:,    11)*sol(:ncol,:,     4)*sol(:ncol,:,     2)                            ! rate_const*H2O2*OH
      rxt_rates(:ncol,:,    12) = rxt_rates(:ncol,:,    12)*sol(:ncol,:,     5)*sol(:ncol,:,     1)                            ! rate_const*NO*O3
      rxt_rates(:ncol,:,    13) = rxt_rates(:ncol,:,    13)*sol(:ncol,:,     3)*sol(:ncol,:,     5)                            ! rate_const*HO2*NO
      rxt_rates(:ncol,:,    14) = rxt_rates(:ncol,:,    14)*sol(:ncol,:,     6)*sol(:ncol,:,     2)                            ! rate_const*M*NO2*OH
      rxt_rates(:ncol,:,    15) = rxt_rates(:ncol,:,    15)*sol(:ncol,:,     2)                                                ! rate_const*CH4*OH
      rxt_rates(:ncol,:,    16) = rxt_rates(:ncol,:,    16)*sol(:ncol,:,     8)*sol(:ncol,:,     2)                            ! rate_const*CO*OH
      rxt_rates(:ncol,:,    17) = rxt_rates(:ncol,:,    17)*sol(:ncol,:,     9)*sol(:ncol,:,     2)                            ! rate_const*CH2O*OH
      rxt_rates(:ncol,:,    18) = rxt_rates(:ncol,:,    18)*sol(:ncol,:,    10)*sol(:ncol,:,     3)                            ! rate_const*CH3O2*HO2
      rxt_rates(:ncol,:,    19) = rxt_rates(:ncol,:,    19)*sol(:ncol,:,    11)*sol(:ncol,:,     2)                            ! rate_const*CH3OOH*OH
      rxt_rates(:ncol,:,    20) = rxt_rates(:ncol,:,    20)*sol(:ncol,:,    11)*sol(:ncol,:,     2)                            ! rate_const*CH3OOH*OH
      rxt_rates(:ncol,:,    21) = rxt_rates(:ncol,:,    21)*sol(:ncol,:,    10)*sol(:ncol,:,     5)                            ! rate_const*CH3O2*NO
      rxt_rates(:ncol,:,    22) = rxt_rates(:ncol,:,    22)*sol(:ncol,:,    10)*sol(:ncol,:,    10)                            ! rate_const*CH3O2*CH3O2
      rxt_rates(:ncol,:,    23) = rxt_rates(:ncol,:,    23)*sol(:ncol,:,     6)                                                ! rate_const*H2O*NO2
      rxt_rates(:ncol,:,    24) = rxt_rates(:ncol,:,    24)*sol(:ncol,:,    12)*sol(:ncol,:,     2)                            ! rate_const*DMS*OH
      rxt_rates(:ncol,:,    25) = rxt_rates(:ncol,:,    25)*sol(:ncol,:,    12)*sol(:ncol,:,     2)                            ! rate_const*DMS*OH
      rxt_rates(:ncol,:,    26) = rxt_rates(:ncol,:,    26)*sol(:ncol,:,    13)*sol(:ncol,:,     2)                            ! rate_const*SO2*OH
      rxt_rates(:ncol,:,    27) = rxt_rates(:ncol,:,    27)*sol(:ncol,:,    14)*sol(:ncol,:,     2)                            ! rate_const*ISOP*OH
      rxt_rates(:ncol,:,    28) = rxt_rates(:ncol,:,    28)*sol(:ncol,:,    14)*sol(:ncol,:,     1)                            ! rate_const*ISOP*O3
  end subroutine set_rates
end module mo_rxt_rates_conv
