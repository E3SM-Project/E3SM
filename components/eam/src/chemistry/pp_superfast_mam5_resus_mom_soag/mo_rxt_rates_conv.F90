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
                                                                                                                               ! 
                                                                                                                               ! rate_const*O2
      rxt_rates(:ncol,:,     3) = rxt_rates(:ncol,:,     3)*sol(:ncol,:,     4)                                                ! rate_const*H2O2
      rxt_rates(:ncol,:,     4) = rxt_rates(:ncol,:,     4)*sol(:ncol,:,     5)                                                ! rate_const*CH2O
      rxt_rates(:ncol,:,     5) = rxt_rates(:ncol,:,     5)*sol(:ncol,:,     5)                                                ! rate_const*CH2O
      rxt_rates(:ncol,:,     6) = rxt_rates(:ncol,:,     6)*sol(:ncol,:,     7)                                                ! rate_const*CH3OOH
      rxt_rates(:ncol,:,     7) = rxt_rates(:ncol,:,     7)*sol(:ncol,:,    22)                                                ! rate_const*C2H5OOH
      rxt_rates(:ncol,:,     8) = rxt_rates(:ncol,:,     8)*sol(:ncol,:,     9)                                                ! rate_const*NO2
      rxt_rates(:ncol,:,     9) = rxt_rates(:ncol,:,     9)*sol(:ncol,:,    10)                                                ! rate_const*NO3
      rxt_rates(:ncol,:,    10) = rxt_rates(:ncol,:,    10)*sol(:ncol,:,    10)                                                ! rate_const*NO3
      rxt_rates(:ncol,:,    11) = rxt_rates(:ncol,:,    11)*sol(:ncol,:,    11)                                                ! rate_const*N2O5
      rxt_rates(:ncol,:,    12) = rxt_rates(:ncol,:,    12)*sol(:ncol,:,    11)                                                ! rate_const*N2O5
      rxt_rates(:ncol,:,    13) = rxt_rates(:ncol,:,    13)*sol(:ncol,:,    12)                                                ! rate_const*HNO3
      rxt_rates(:ncol,:,    14) = rxt_rates(:ncol,:,    14)*sol(:ncol,:,    13)                                                ! rate_const*HO2NO2
      rxt_rates(:ncol,:,    15) = rxt_rates(:ncol,:,    15)*sol(:ncol,:,    13)                                                ! rate_const*HO2NO2
      rxt_rates(:ncol,:,    16) = rxt_rates(:ncol,:,    16)*sol(:ncol,:,    23)                                                ! rate_const*CH3CHO
      rxt_rates(:ncol,:,    17) = rxt_rates(:ncol,:,    17)*sol(:ncol,:,    14)                                                ! rate_const*PAN
      rxt_rates(:ncol,:,    18) = rxt_rates(:ncol,:,    18)*sol(:ncol,:,    20)                                                ! rate_const*CH3COCH3
      rxt_rates(:ncol,:,    19) = rxt_rates(:ncol,:,    19)*sol(:ncol,:,    29)                                                ! rate_const*MVKMACR
      rxt_rates(:ncol,:,    20) = rxt_rates(:ncol,:,    20)*sol(:ncol,:,    54)                                                ! rate_const*soa_a1
      rxt_rates(:ncol,:,    21) = rxt_rates(:ncol,:,    21)*sol(:ncol,:,    55)                                                ! rate_const*soa_a2
      rxt_rates(:ncol,:,    22) = rxt_rates(:ncol,:,    22)*sol(:ncol,:,    56)                                                ! rate_const*soa_a3
      rxt_rates(:ncol,:,    23) = rxt_rates(:ncol,:,    23)*sol(:ncol,:,     1)                                                ! rate_const*H2O*O3
      rxt_rates(:ncol,:,    24) = rxt_rates(:ncol,:,    24)*sol(:ncol,:,     1)                                                ! rate_const*H2*O3
      rxt_rates(:ncol,:,    25) = rxt_rates(:ncol,:,    25)*sol(:ncol,:,     1)*sol(:ncol,:,    34)                            ! rate_const*O3*CH4LNZ
      rxt_rates(:ncol,:,    26) = rxt_rates(:ncol,:,    26)*sol(:ncol,:,    15)*sol(:ncol,:,     2)                            ! rate_const*CO*OH
      rxt_rates(:ncol,:,    27) = rxt_rates(:ncol,:,    27)*sol(:ncol,:,    15)*sol(:ncol,:,     2)                            ! rate_const*M*CO*OH
      rxt_rates(:ncol,:,    28) = rxt_rates(:ncol,:,    28)*sol(:ncol,:,     2)                                                ! rate_const*H2*OH
      rxt_rates(:ncol,:,    29) = rxt_rates(:ncol,:,    29)*sol(:ncol,:,    34)*sol(:ncol,:,     2)                            ! rate_const*CH4LNZ*OH
      rxt_rates(:ncol,:,    30) = rxt_rates(:ncol,:,    30)*sol(:ncol,:,    16)*sol(:ncol,:,     2)                            ! rate_const*C2H6*OH
      rxt_rates(:ncol,:,    31) = rxt_rates(:ncol,:,    31)*sol(:ncol,:,    17)*sol(:ncol,:,     2)                            ! rate_const*C3H8*OH
      rxt_rates(:ncol,:,    32) = rxt_rates(:ncol,:,    32)*sol(:ncol,:,    18)*sol(:ncol,:,     2)                            ! rate_const*M*C2H4*OH
      rxt_rates(:ncol,:,    33) = rxt_rates(:ncol,:,    33)*sol(:ncol,:,    18)*sol(:ncol,:,     1)                            ! rate_const*C2H4*O3
      rxt_rates(:ncol,:,    34) = rxt_rates(:ncol,:,    34)*sol(:ncol,:,    25)*sol(:ncol,:,     1)                            ! rate_const*ISOP*O3
      rxt_rates(:ncol,:,    35) = rxt_rates(:ncol,:,    35)*sol(:ncol,:,    25)*sol(:ncol,:,     2)                            ! rate_const*ISOP*OH
      rxt_rates(:ncol,:,    36) = rxt_rates(:ncol,:,    36)*sol(:ncol,:,     5)*sol(:ncol,:,     2)                            ! rate_const*CH2O*OH
      rxt_rates(:ncol,:,    37) = rxt_rates(:ncol,:,    37)*sol(:ncol,:,     2)*sol(:ncol,:,     1)                            ! rate_const*OH*O3
      rxt_rates(:ncol,:,    38) = rxt_rates(:ncol,:,    38)*sol(:ncol,:,     2)*sol(:ncol,:,     2)                            ! rate_const*OH*OH
      rxt_rates(:ncol,:,    39) = rxt_rates(:ncol,:,    39)*sol(:ncol,:,     3)*sol(:ncol,:,     1)                            ! rate_const*HO2*O3
      rxt_rates(:ncol,:,    40) = rxt_rates(:ncol,:,    40)*sol(:ncol,:,     3)*sol(:ncol,:,     2)                            ! rate_const*HO2*OH
      rxt_rates(:ncol,:,    41) = rxt_rates(:ncol,:,    41)*sol(:ncol,:,     3)*sol(:ncol,:,     3)                            ! rate_const*M*HO2*HO2
      rxt_rates(:ncol,:,    42) = rxt_rates(:ncol,:,    42)*sol(:ncol,:,     3)*sol(:ncol,:,     3)                            ! rate_const*HO2*HO2
      rxt_rates(:ncol,:,    43) = rxt_rates(:ncol,:,    43)*sol(:ncol,:,     2)*sol(:ncol,:,     2)                            ! rate_const*M*OH*OH
      rxt_rates(:ncol,:,    44) = rxt_rates(:ncol,:,    44)*sol(:ncol,:,     4)*sol(:ncol,:,     2)                            ! rate_const*H2O2*OH
      rxt_rates(:ncol,:,    45) = rxt_rates(:ncol,:,    45)*sol(:ncol,:,     8)*sol(:ncol,:,     1)                            ! rate_const*NO*O3
      rxt_rates(:ncol,:,    46) = rxt_rates(:ncol,:,    46)*sol(:ncol,:,     8)*sol(:ncol,:,     3)                            ! rate_const*NO*HO2
      rxt_rates(:ncol,:,    47) = rxt_rates(:ncol,:,    47)*sol(:ncol,:,     9)*sol(:ncol,:,     1)                            ! rate_const*NO2*O3
      rxt_rates(:ncol,:,    48) = rxt_rates(:ncol,:,    48)*sol(:ncol,:,     2)*sol(:ncol,:,    10)                            ! rate_const*OH*NO3
      rxt_rates(:ncol,:,    49) = rxt_rates(:ncol,:,    49)*sol(:ncol,:,     8)*sol(:ncol,:,    10)                            ! rate_const*NO*NO3
      rxt_rates(:ncol,:,    50) = rxt_rates(:ncol,:,    50)*sol(:ncol,:,    13)*sol(:ncol,:,     2)                            ! rate_const*HO2NO2*OH
      rxt_rates(:ncol,:,    51) = rxt_rates(:ncol,:,    51)*sol(:ncol,:,    12)*sol(:ncol,:,     2)                            ! rate_const*HNO3*OH
      rxt_rates(:ncol,:,    52) = rxt_rates(:ncol,:,    52)*sol(:ncol,:,    12)*sol(:ncol,:,     2)                            ! rate_const*HNO3*OH
      rxt_rates(:ncol,:,    53) = rxt_rates(:ncol,:,    53)*sol(:ncol,:,     9)*sol(:ncol,:,     2)                            ! rate_const*M*NO2*OH
      rxt_rates(:ncol,:,    54) = rxt_rates(:ncol,:,    54)*sol(:ncol,:,     9)*sol(:ncol,:,     3)                            ! rate_const*M*NO2*HO2
      rxt_rates(:ncol,:,    55) = rxt_rates(:ncol,:,    55)*sol(:ncol,:,     9)*sol(:ncol,:,    10)                            ! rate_const*M*NO2*NO3
      rxt_rates(:ncol,:,    56) = rxt_rates(:ncol,:,    56)*sol(:ncol,:,    24)*sol(:ncol,:,     9)                            ! rate_const*M*CH3CO3*NO2
      rxt_rates(:ncol,:,    57) = rxt_rates(:ncol,:,    57)*sol(:ncol,:,    13)                                                ! rate_const*M*HO2NO2
      rxt_rates(:ncol,:,    58) = rxt_rates(:ncol,:,    58)*sol(:ncol,:,    11)                                                ! rate_const*M*N2O5
      rxt_rates(:ncol,:,    59) = rxt_rates(:ncol,:,    59)*sol(:ncol,:,    14)                                                ! rate_const*M*PAN
      rxt_rates(:ncol,:,    60) = rxt_rates(:ncol,:,    60)*sol(:ncol,:,     6)*sol(:ncol,:,     3)                            ! rate_const*CH3O2*HO2
      rxt_rates(:ncol,:,    61) = rxt_rates(:ncol,:,    61)*sol(:ncol,:,     6)*sol(:ncol,:,     8)                            ! rate_const*CH3O2*NO
      rxt_rates(:ncol,:,    62) = rxt_rates(:ncol,:,    62)*sol(:ncol,:,     6)*sol(:ncol,:,     6)                            ! rate_const*CH3O2*CH3O2
      rxt_rates(:ncol,:,    63) = rxt_rates(:ncol,:,    63)*sol(:ncol,:,     7)*sol(:ncol,:,     2)                            ! rate_const*CH3OOH*OH
      rxt_rates(:ncol,:,    64) = rxt_rates(:ncol,:,    64)*sol(:ncol,:,    21)*sol(:ncol,:,     8)                            ! rate_const*C2H5O2*NO
      rxt_rates(:ncol,:,    65) = rxt_rates(:ncol,:,    65)*sol(:ncol,:,    21)*sol(:ncol,:,    21)                            ! rate_const*C2H5O2*C2H5O2
      rxt_rates(:ncol,:,    66) = rxt_rates(:ncol,:,    66)*sol(:ncol,:,    21)*sol(:ncol,:,     6)                            ! rate_const*C2H5O2*CH3O2
      rxt_rates(:ncol,:,    67) = rxt_rates(:ncol,:,    67)*sol(:ncol,:,    21)*sol(:ncol,:,     3)                            ! rate_const*C2H5O2*HO2
      rxt_rates(:ncol,:,    68) = rxt_rates(:ncol,:,    68)*sol(:ncol,:,    22)*sol(:ncol,:,     2)                            ! rate_const*C2H5OOH*OH
      rxt_rates(:ncol,:,    69) = rxt_rates(:ncol,:,    69)*sol(:ncol,:,    22)*sol(:ncol,:,     2)                            ! rate_const*C2H5OOH*OH
      rxt_rates(:ncol,:,    70) = rxt_rates(:ncol,:,    70)*sol(:ncol,:,    23)*sol(:ncol,:,     2)                            ! rate_const*CH3CHO*OH
      rxt_rates(:ncol,:,    71) = rxt_rates(:ncol,:,    71)*sol(:ncol,:,    23)*sol(:ncol,:,    10)                            ! rate_const*CH3CHO*NO3
      rxt_rates(:ncol,:,    72) = rxt_rates(:ncol,:,    72)*sol(:ncol,:,    24)*sol(:ncol,:,     8)                            ! rate_const*CH3CO3*NO
      rxt_rates(:ncol,:,    73) = rxt_rates(:ncol,:,    73)*sol(:ncol,:,    24)*sol(:ncol,:,     6)                            ! rate_const*CH3CO3*CH3O2
      rxt_rates(:ncol,:,    74) = rxt_rates(:ncol,:,    74)*sol(:ncol,:,    24)*sol(:ncol,:,    24)                            ! rate_const*CH3CO3*CH3CO3
      rxt_rates(:ncol,:,    75) = rxt_rates(:ncol,:,    75)*sol(:ncol,:,    20)*sol(:ncol,:,     2)                            ! rate_const*CH3COCH3*OH
      rxt_rates(:ncol,:,    76) = rxt_rates(:ncol,:,    76)*sol(:ncol,:,    20)*sol(:ncol,:,     2)                            ! rate_const*CH3COCH3*OH
      rxt_rates(:ncol,:,    77) = rxt_rates(:ncol,:,    77)*sol(:ncol,:,    19)*sol(:ncol,:,     8)                            ! rate_const*ROHO2*NO
      rxt_rates(:ncol,:,    78) = rxt_rates(:ncol,:,    78)*sol(:ncol,:,    19)*sol(:ncol,:,     3)                            ! rate_const*ROHO2*HO2
      rxt_rates(:ncol,:,    79) = rxt_rates(:ncol,:,    79)*sol(:ncol,:,    19)*sol(:ncol,:,     6)                            ! rate_const*ROHO2*CH3O2
      rxt_rates(:ncol,:,    80) = rxt_rates(:ncol,:,    80)*sol(:ncol,:,    27)*sol(:ncol,:,     8)                            ! rate_const*ISOPO2*NO
      rxt_rates(:ncol,:,    81) = rxt_rates(:ncol,:,    81)*sol(:ncol,:,    27)*sol(:ncol,:,     3)                            ! rate_const*ISOPO2*HO2
      rxt_rates(:ncol,:,    82) = rxt_rates(:ncol,:,    82)*sol(:ncol,:,    27)*sol(:ncol,:,     6)                            ! rate_const*ISOPO2*CH3O2
      rxt_rates(:ncol,:,    83) = rxt_rates(:ncol,:,    83)*sol(:ncol,:,    29)*sol(:ncol,:,     1)                            ! rate_const*MVKMACR*O3
      rxt_rates(:ncol,:,    84) = rxt_rates(:ncol,:,    84)*sol(:ncol,:,    29)*sol(:ncol,:,     2)                            ! rate_const*MVKMACR*OH
      rxt_rates(:ncol,:,    85) = rxt_rates(:ncol,:,    85)*sol(:ncol,:,    30)*sol(:ncol,:,     8)                            ! rate_const*MVKO2*NO
      rxt_rates(:ncol,:,    86) = rxt_rates(:ncol,:,    86)*sol(:ncol,:,    30)*sol(:ncol,:,     3)                            ! rate_const*MVKO2*HO2
      rxt_rates(:ncol,:,    87) = rxt_rates(:ncol,:,    87)*sol(:ncol,:,    31)                                                ! rate_const*E90
      rxt_rates(:ncol,:,    88) = rxt_rates(:ncol,:,    88)*sol(:ncol,:,    36)*sol(:ncol,:,     2)                            ! rate_const*DMS*OH
      rxt_rates(:ncol,:,    89) = rxt_rates(:ncol,:,    89)*sol(:ncol,:,    36)*sol(:ncol,:,     2)                            ! rate_const*DMS*OH
      rxt_rates(:ncol,:,    90) = rxt_rates(:ncol,:,    90)*sol(:ncol,:,    37)*sol(:ncol,:,     2)                            ! rate_const*SO2*OH
      rxt_rates(:ncol,:,    91) = rxt_rates(:ncol,:,    91)*sol(:ncol,:,    36)*sol(:ncol,:,    10)                            ! rate_const*DMS*NO3
      rxt_rates(:ncol,:,    92) = rxt_rates(:ncol,:,    92)*sol(:ncol,:,    39)                                                ! rate_const*prsd_OH*SOAG0
      rxt_rates(:ncol,:,    93) = rxt_rates(:ncol,:,    93)*sol(:ncol,:,    40)                                                ! rate_const*prsd_OH*SOAG15
      rxt_rates(:ncol,:,    94) = rxt_rates(:ncol,:,    94)*sol(:ncol,:,    41)                                                ! rate_const*prsd_OH*SOAG24
      rxt_rates(:ncol,:,    95) = rxt_rates(:ncol,:,    95)*sol(:ncol,:,    26)                                                ! rate_const*prsd_OH*ISOP_VBS
      rxt_rates(:ncol,:,    96) = rxt_rates(:ncol,:,    96)*sol(:ncol,:,    28)                                                ! rate_const*prsd_OH*C10H16
      rxt_rates(:ncol,:,    97) = rxt_rates(:ncol,:,    97)*sol(:ncol,:,    26)                                                ! rate_const*prsd_O3*ISOP_VBS
      rxt_rates(:ncol,:,    98) = rxt_rates(:ncol,:,    98)*sol(:ncol,:,    28)                                                ! rate_const*prsd_O3*C10H16
      rxt_rates(:ncol,:,    99) = rxt_rates(:ncol,:,    99)*sol(:ncol,:,    26)                                                ! rate_const*prsd_NO3*ISOP_VBS
      rxt_rates(:ncol,:,   100) = rxt_rates(:ncol,:,   100)*sol(:ncol,:,    28)                                                ! rate_const*prsd_NO3*C10H16
      rxt_rates(:ncol,:,   101) = rxt_rates(:ncol,:,   101)*sol(:ncol,:,    46)                                                ! rate_const*prsd_OH*SOAG35
      rxt_rates(:ncol,:,   102) = rxt_rates(:ncol,:,   102)*sol(:ncol,:,    45)                                                ! rate_const*prsd_OH*SOAG34
      rxt_rates(:ncol,:,   103) = rxt_rates(:ncol,:,   103)*sol(:ncol,:,    44)                                                ! rate_const*prsd_OH*SOAG33
      rxt_rates(:ncol,:,   104) = rxt_rates(:ncol,:,   104)*sol(:ncol,:,    43)                                                ! rate_const*prsd_OH*SOAG32
  end subroutine set_rates
end module mo_rxt_rates_conv
