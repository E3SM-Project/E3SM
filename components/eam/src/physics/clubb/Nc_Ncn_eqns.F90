!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module Nc_Ncn_eqns

  ! Description:
  ! Equations are provided to perform calculations back-and-forth between Nc and
  ! Ncn, where Nc is cloud droplet concentration and Ncn is simplified cloud
  ! nuclei concentration.  The equation that relates the two is:
  !
  ! Nc = Ncn * H(chi);
  !
  ! where chi is extended liquid water mixing ratio, which is equal to cloud
  ! water mixing ratio, rc, when both are positive.  However, chi is negative in
  ! subsaturated air.
  !
  ! Equation are provided relating mean cloud droplet concentration (overall),
  ! Ncm, and/or mean cloud droplet concentration (in-cloud), Nc_in_cloud, to
  ! mean simplified cloud nuclei concentration, Ncnm.

  ! Notes:
  !
  ! Meaning of Nc flag combinations:
  !
  ! l_const_Nc_in_cloud:
  ! When this flag is enabled, cloud droplet concentration (in-cloud) is
  ! constant (spatially) at a grid level (it is constant over the subgrid
  ! domain, but could vary over time depending on the value of l_predict_Nc).
  ! The value of in-cloud Nc does not vary at a grid level.  This also means
  ! that Ncn is constant across the entire grid level.  When this flag is turned
  ! off, both in-cloud Nc and Ncn vary at a grid level.
  !
  ! l_predict_Nc:
  ! When this flag is enabled, Nc_in_cloud (or alternatively Ncm) is predicted.
  ! It is advanced every time step by a predictive equation, and can change
  ! at every time step at a grid level.  When this flag is turned off,
  ! Nc_in_cloud does not change at a grid level over the course of a model run.
  !
  ! 1) l_predict_Nc turned on and l_const_Nc_in_cloud turned on:
  !    The value of Nc_in_cloud (mean in-cloud Nc) is predicted and can change
  !    at every timestep at a grid level.  However, the value of in-cloud Nc is
  !    constant (spatially) at a grid level (no subgrid variability).
  !
  ! 2) l_predict_Nc turned on and l_const_Nc_in_cloud turned off:
  !    The value of Nc_in_cloud (mean in-cloud Nc) is predicted and can change
  !    at every timestep at a grid level.  The value of in-cloud Nc also varies
  !    (spatially) at a grid level (subgrid variability around mean
  !    in-cloud Nc).
  !
  ! 3) l_predict_Nc turned off and l_const_Nc_in_cloud turned on:
  !    The value of Nc_in_cloud (mean in-cloud Nc) is constant over time at a
  !    grid level.  It retains its initial value.  Additionally, the value of
  !    in-cloud Nc is constant (spatially) at a grid level (no subgrid
  !    variability).  This configuration is used most often in idealized cases.
  !
  ! 4) l_predict_Nc turned off and l_const_Nc_in_cloud turned off:
  !    The value of Nc_in_cloud (mean in-cloud Nc) is constant over time at a
  !    grid level.  It retains its initial value.  However, the value of
  !    in-cloud Nc varies (spatially) at a grid level (subgrid variability
  !    around mean in-cloud Nc).
  !
  !
  !
  ! Nc_in_cloud/Nc - Ncn flow chart of CLUBB code:
  !
  ! (Please update when warranted).
  !
  !
  ! Ncm/Nc-in-cloud                                         Ncnm/Ncn PDF params.
  ! --->
  ! |  |                    Start of CLUBB main time step loop
  ! |  |
  ! |  |                    advance_clubb_core
  ! |  |
  ! |  |
  ! |  |\
  ! |  | \
  ! |  |  (intent in)-------setup_pdf_parameters-------->calc. Ncnm (local)
  ! |  |                                                       |
  ! |  |                                                      \ /
  ! |  |                                             mu_Ncn_i, sigma_Ncn_i,
  ! |  |                                                  corr_xNcn_i
  ! |  |                                                       |
  ! |  |                                                      \ /
  ! |  |                                               PDF param. arrays:
  ! |  |                                             mu_x_i_n, sigma_x_i_n,
  ! |  |                                                 corr_array_i_n
  ! |  |                                                  (intent out)
  ! |  |                                                       |
  ! |  |                                                       |
  ! |  |                                                       |
  ! |  |                                                       |
  ! |  |                                                       |
  ! |  |                                                       |
  ! |  |--(intent in)---calc_microphys_scheme_tendcies----(intent in)
  ! |  |                            |
  ! |  |                            |
  ! |  |                call a microphysics scheme
  ! |  |                            | 
  ! |  |  Local micro. scheme-----------Latin Hypercube-----------Upscaled KK
  ! |  |         |                            |                        |
  ! |  |  Ncm/Nc-in-cloud:          Populate sample points      Use PDF params.
  ! |  |  used to find micro.       using PDF params (Ncn).         of Ncn
  ! |  |    tendencies.             At every sample point:      (mu_Ncn_i, etc.)
  ! |  |         |                    Nc = Ncn * H(chi).         to find micro.
  ! |  |         |                  Use sample-point Nc to         tendencies.
  ! |  |         |                  find micro. tendencies             |
  ! |  |         |                when calling micro. scheme.          |
  ! |  |         |                            |                        |
  ! |  |    hydromet_mc/-----------------hydromet_mc/-------------hydromet_mc
  ! |  |    Ncm_mc (intent out)     |    Ncm_mc (intent out)      (intent out)
  ! |  |                            |
  ! |  |                            |
  ! |  |                            |
  ! |  |                            |
  ! |  |                            |
  ! |  |                            |
  ! |  |                            |
  ! |  |                        (intent in)
  ! |  |                            |
  ! |  |--(intent inout)----advance_microphys
  ! |  |
  ! |  |
  ! |  | advance microphysics variables (hydromet, Nc_in_cloud/Ncm) one timestep
  ! |  |                
  ! |  |                   l_predict_Nc = true:
  ! |  |            Nc_in_cloud/Ncm necessary for starting
  ! |  |            value of Nc_in_cloud/Ncm when advancing
  ! |  |            one timestep using predictive equation.
  ! |  |
  ! |  |
  ! |  |                    End of CLUBB main time step loop
  ! <---

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  private ! default scope

  public :: Ncnm_to_Nc_in_cloud, &
            Nc_in_cloud_to_Ncnm, &
            Ncnm_to_Ncm, &
            Ncm_to_Ncnm

  private :: bivar_NL_chi_Ncn_mean, &
             bivar_Ncnm_eqn_comp

  contains

  !=============================================================================
  elemental function Ncnm_to_Nc_in_cloud( mu_chi_1, mu_chi_2, mu_Ncn_1, &
                                          mu_Ncn_2, sigma_chi_1, sigma_chi_2, &
                                          sigma_Ncn_1, sigma_Ncn_2, &
                                          sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                          corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                          mixt_frac, cloud_frac_1, &
                                          cloud_frac_2 ) &
  result( Nc_in_cloud )

    ! Description:
    ! The in-cloud mean of cloud droplet concentration is calculated from the
    ! PDF parameters involving simplified cloud nuclei concentration, Ncn, and
    ! cloud fraction.  At any point, cloud droplet concentration, Nc, is given
    ! by:
    !
    ! Nc = Ncn * H(chi);
    !
    ! where extended liquid water mixing ratio, chi, is equal to cloud water
    ! ratio, rc, when positive.  When the atmosphere is saturated at this point,
    ! cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    ! and Nc = 0.
    !
    ! The overall mean of cloud droplet concentration, <Nc>, is calculated from
    ! the PDF parameters involving Ncn.  The in-cloud mean of cloud droplet
    ! concentration is calculated from <Nc> and cloud fraction.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        cloud_frac_min

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_1,         & ! Mean of chi (old s) (1st PDF component)      [kg/kg]
      mu_chi_2,         & ! Mean of chi (old s) (2nd PDF component)      [kg/kg]
      mu_Ncn_1,         & ! Mean of Ncn (1st PDF component)             [num/kg]
      mu_Ncn_2,         & ! Mean of Ncn (2nd PDF component)             [num/kg]
      sigma_chi_1,      & ! Standard deviation of chi (1st PDF comp.)    [kg/kg]
      sigma_chi_2,      & ! Standard deviation of chi (2nd PDF comp.)    [kg/kg]
      sigma_Ncn_1,      & ! Standard deviation of Ncn (1st PDF comp.)   [num/kg]
      sigma_Ncn_2,      & ! Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
      sigma_Ncn_1_n,    & ! Standard deviation of ln Ncn (1st PDF component) [-]
      sigma_Ncn_2_n,    & ! Standard deviation of ln Ncn (2nd PDF component) [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)    [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
      mixt_frac,        & ! Mixture fraction                                 [-]
      cloud_frac_1,     & ! Cloud fraction (1st PDF component)               [-]
      cloud_frac_2        ! Cloud fraction (2nd PDF component)               [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      Nc_in_cloud    ! Mean cloud droplet concentration (in-cloud)      [num/kg]

    ! Local Variable
    real( kind = core_rknd ) :: &
      Ncm,        & ! Mean cloud droplet concentration (overall)    [num/kg]
      cloud_frac    ! Cloud fraction                                [-]
 

    ! Calculate overall cloud fraction as calculated by the PDF.
    ! The variable cloud_frac is not used here because it is altered by factors
    ! such as the trapezoidal rule calculation.
    ! Cloud fraction can be recalculated here from cloud_frac_1 and cloud_frac_2
    ! as long neither of these variables are altered by any factor.  They can
    ! only be calculated from PDF.
    cloud_frac = mixt_frac * cloud_frac_1 + ( one - mixt_frac ) * cloud_frac_2

    if ( cloud_frac > cloud_frac_min ) then

       ! There is cloud found at this grid level.  Calculate Nc_in_cloud.
       Ncm = Ncnm_to_Ncm( mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, &
                          sigma_chi_1, sigma_chi_2, sigma_Ncn_1, &
                          sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                          corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, mixt_frac )

       Nc_in_cloud = Ncm / cloud_frac

    else ! cloud_frac <= cloud_frac_min

       ! This level is entirely clear.  Set Nc_in_cloud to <Ncn>.
       ! Since <Ncn> = mu_Ncn_1 = mu_Ncn_2, use mu_Ncn_1 here.
       Nc_in_cloud = mu_Ncn_1

    endif


    return

  end function Ncnm_to_Nc_in_cloud

  !=============================================================================
  elemental function Nc_in_cloud_to_Ncnm( mu_chi_1, mu_chi_2, sigma_chi_1, &
                                          sigma_chi_2, mixt_frac, Nc_in_cloud, &
                                          cloud_frac_1, cloud_frac_2, &
                                          const_Ncnp2_on_Ncnm2, &
                                          const_corr_chi_Ncn ) &
  result( Ncnm )

    ! Description:
    ! The overall mean of simplified cloud nuclei concentration, <Ncn>, is
    ! calculated from the in-cloud mean of cloud droplet concentration, <Nc>,
    ! cloud fraction, and some of the PDF parameters.
    !
    ! At any point, cloud droplet concentration, Nc, is given by:
    !
    ! Nc = Ncn * H(chi);
    !
    ! where extended liquid water mixing ratio, chi, is equal to cloud water
    ! ratio, rc, when positive.  When the atmosphere is saturated at this point,
    ! cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    ! and Nc = 0.
    !
    ! The overall mean of cloud droplet concentration, <Nc>, is calculated from
    ! Nc_in_cloud and cloud fraction.  The value of <Ncn> is calculated from
    ! <Nc> and PDF parameters.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        cloud_frac_min, &
        eps

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_1,    & ! Mean of chi (old s) (1st PDF component)           [kg/kg]
      mu_chi_2,    & ! Mean of chi (old s) (2nd PDF component)           [kg/kg]
      sigma_chi_1, & ! Standard deviation of chi (1st PDF component)     [kg/kg]
      sigma_chi_2, & ! Standard deviation of chi (2nd PDF component)     [kg/kg]
      mixt_frac      ! Mixture fraction                                      [-]

    real( kind = core_rknd ), intent(in) :: &
      Nc_in_cloud,          & ! Mean cloud droplet conc. (in-cloud)     [num/kg]
      cloud_frac_1,         & ! Cloud fraction (1st PDF component)           [-]
      cloud_frac_2,         & ! Cloud fraction (2nd PDF component)           [-]
      const_Ncnp2_on_Ncnm2, & ! Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
      const_corr_chi_Ncn      ! Prescribed correlation of chi and Ncn        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      Ncnm    ! Mean simplified cloud nuclei concentration (overall)    [num/kg]

    ! Local Variable
    real( kind = core_rknd ) :: &
      Ncm,        & ! Mean cloud droplet concentration (overall)      [num/kg]
      cloud_frac    ! Cloud fraction                                  [-]


    ! Calculate overall cloud fraction as calculated by the PDF.
    ! The variable cloud_frac is not used here because it is altered by factors
    ! such as the trapezoidal rule calculation.
    ! Cloud fraction can be recalculated here from cloud_frac_1 and cloud_frac_2
    ! as long neither of these variables are altered by any factor.  They can
    ! only be calculated from the PDF.
    cloud_frac = mixt_frac * cloud_frac_1 + ( one - mixt_frac ) * cloud_frac_2

    if ( cloud_frac > cloud_frac_min &
         .and. abs(const_corr_chi_Ncn * const_Ncnp2_on_Ncnm2) > eps ) then

       ! There is cloud found at this grid level.  Additionally, Ncn varies.
       ! Calculate Nc_in_cloud.
       Ncm = Nc_in_cloud * cloud_frac

       Ncnm = Ncm_to_Ncnm( mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, &
                           mixt_frac, Ncm, const_Ncnp2_on_Ncnm2, &
                           const_corr_chi_Ncn, Nc_in_cloud )

    else ! cloud_frac <= cloud_frac_min .or. const_Ncnp2_on_Ncnm2 = 0

       ! When Ncn is constant a a grid level, it is equal to Nc_in_cloud.
       ! Additionally, when a level is entirely clear, <Ncn>, which is based on
       ! Nc_in_cloud, here, must be set to something.  Set <Ncn> to Nc_in_cloud.
       Ncnm = Nc_in_cloud

    endif


    return

  end function Nc_in_cloud_to_Ncnm

  !=============================================================================
  elemental function Ncnm_to_Ncm( mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, &
                                  sigma_chi_1, sigma_chi_2, sigma_Ncn_1, &
                                  sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                  corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, &
                                  mixt_frac ) &
  result( Ncm )

    ! Description:
    ! The overall mean of cloud droplet concentration, <Nc>, is calculated from
    ! the PDF parameters involving the simplified cloud nuclei concentration,
    ! Ncn.  At any point, cloud droplet concentration, Nc, is given by:
    !
    ! Nc = Ncn * H(chi);
    !
    ! where extended liquid water mixing ratio, chi, is equal to cloud water
    ! ratio, rc, when positive.  When the atmosphere is saturated at this point,
    ! cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    ! and Nc = 0.
    !
    ! The overall mean of cloud droplet concentration, <Nc>, is found by
    ! integrating over the PDF of chi and Ncn, such that:
    !
    ! <Nc> = INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P(chi,Ncn) dNcn dchi;
    !
    ! which can also be written as:
    !
    ! <Nc> = SUM(i=1,n) mixt_frac_i
    !        * INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;
    !
    ! where n is the number of multivariate joint PDF components, mixt_frac_i is
    ! the weight of the ith PDF component, and P_i is the functional form of the
    ! multivariate joint PDF in the ith PDF component.
    !
    ! This equation is rewritten as:
    !
    ! <Nc> = SUM(i=1,n) mixt_frac_i
    !        * INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi.
    !
    ! When both chi and Ncn vary in the ith PDF component, the integral is
    ! evaluated and the result is:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = (1/2) * exp{ mu_Ncn_i_n + (1/2) * sigma_Ncn_i_n^2 }
    !   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
    !                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) );
    !
    ! which can be reduced to:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = (1/2) * mu_Ncn_i
    !   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
    !                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) ).
    !
    ! When chi is constant, but Ncn varies, in the ith PDF component, the
    ! integral is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
    !
    ! when mu_chi_i > 0; and
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
    !
    ! when mu_chi_i <= 0.
    !
    ! When chi varies, but Ncn is constant, in the ith PDF component, the
    ! integral is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = mu_Ncn_i * (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! When both chi and Ncn are constant in the ith PDF component, the integral
    ! is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
    !
    ! when mu_chi_i > 0; and
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
    !
    ! when mu_chi_i <= 0.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_1,         & ! Mean of chi (old s) (1st PDF component)      [kg/kg]
      mu_chi_2,         & ! Mean of chi (old s) (2nd PDF component)      [kg/kg]
      mu_Ncn_1,         & ! Mean of Ncn (1st PDF component)             [num/kg]
      mu_Ncn_2,         & ! Mean of Ncn (2nd PDF component)             [num/kg]
      sigma_chi_1,      & ! Standard deviation of chi (1st PDF comp.)    [kg/kg]
      sigma_chi_2,      & ! Standard deviation of chi (2nd PDF comp.)    [kg/kg]
      sigma_Ncn_1,      & ! Standard deviation of Ncn (1st PDF comp.)   [num/kg]
      sigma_Ncn_2,      & ! Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
      sigma_Ncn_1_n,    & ! Standard deviation of ln Ncn (1st PDF component) [-]
      sigma_Ncn_2_n,    & ! Standard deviation of ln Ncn (2nd PDF component) [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi and ln Ncn (1st PDF comp.)    [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi and ln Ncn (2nd PDF comp.)    [-]
      mixt_frac           ! Mixture fraction                                 [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      Ncm    ! Mean cloud droplet concentration (overall)    [num/kg] 


    ! Calculate mean cloud droplet concentration (overall), <Nc>.
    Ncm &
    = mixt_frac &
      * bivar_NL_chi_Ncn_mean( mu_chi_1, mu_Ncn_1, sigma_chi_1, &
                               sigma_Ncn_1, sigma_Ncn_1_n, corr_chi_Ncn_1_n ) &
      + ( one - mixt_frac ) &
        * bivar_NL_chi_Ncn_mean( mu_chi_2, mu_Ncn_2, sigma_chi_2, &
                                 sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n )


    return

  end function Ncnm_to_Ncm

  !=============================================================================
  elemental function Ncm_to_Ncnm( mu_chi_1, mu_chi_2, sigma_chi_1, &
                                  sigma_chi_2, mixt_frac, Ncm, &
                                  const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn, &
                                  Ncnm_val_denom_0 ) &
  result( Ncnm )

    ! Description:
    ! The overall mean of simplified cloud nuclei concentration, <Ncn>, is
    ! calculated from the overall mean of cloud droplet concentration, <Nc>, and
    ! some of the PDF parameters.
    !
    ! At any point, cloud droplet concentration, Nc, is given by:
    !
    ! Nc = Ncn * H(chi);
    !
    ! where extended liquid water mixing ratio, chi, is equal to cloud water
    ! ratio, rc, when positive.  When the atmosphere is saturated at this point,
    ! cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    ! and Nc = 0.
    !
    ! The overall mean of cloud droplet concentration, <Nc>, is found by
    ! integrating over the PDF of chi and Ncn, such that:
    !
    ! <Nc> = INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P(chi,Ncn) dNcn dchi;
    !
    ! which can also be written as:
    !
    ! <Nc> = SUM(i=1,n) mixt_frac_i
    !        * INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;
    !
    ! where n is the number of multivariate joint PDF components, mixt_frac_i is
    ! the weight of the ith PDF component, and P_i is the functional form of the
    ! multivariate joint PDF in the ith PDF component.
    !
    ! This equation is rewritten as:
    !
    ! <Nc> = SUM(i=1,n) mixt_frac_i
    !        * INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi.
    !
    ! When both chi and Ncn vary in the ith PDF component, the integral is
    ! evaluated and the result is:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = (1/2) * exp{ mu_Ncn_i_n + (1/2) * sigma_Ncn_i_n^2 }
    !   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
    !                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) );
    !
    ! which can be reduced to:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = (1/2) * mu_Ncn_i
    !   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
    !                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) ).
    !
    ! When chi is constant, but Ncn varies, in the ith PDF component, the
    ! integral is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
    !
    ! when mu_chi_i > 0; and
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
    !
    ! when mu_chi_i <= 0.
    !
    ! When chi varies, but Ncn is constant, in the ith PDF component, the
    ! integral is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = mu_Ncn_i * (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! When both chi and Ncn are constant in the ith PDF component, the integral
    ! is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
    !
    ! when mu_chi_i > 0; and
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
    !
    ! when mu_chi_i <= 0.
    !
    !
    ! Solving for <Ncn>
    ! =================
    !
    ! The individual marginal for simplified cloud nuclei concentration, Ncn, is
    ! a single lognormal distribution over the entire horizontal domain.  In
    ! order to accomplish this in a two-component PDF structure, the PDF
    ! parameters involving Ncn are set equal between the two components.  This
    ! results in:
    !
    ! mu_Ncn_1 = mu_Ncn_2 = mu_Ncn_i = <Ncn>;
    ! mu_Ncn_1_n = mu_Ncn_2_n = mu_Ncn_i_n;
    ! sigma_Ncn_1 = sigma_Ncn_2 = sigma_Ncn_i = sqrt( <Ncn'^2> );
    ! sigma_Ncn_1_n = sigma_Ncn_2_n = sigma_Ncn_i_n;
    ! rho_chi_Ncn_1 = rho_chi_Ncn_2 = rho_chi_Ncn_i = rho_chi_Ncn; and
    ! rho_chi_Ncn_1_n = rho_chi_Ncn_2_n = rho_chi_Ncn_i_n.
    !
    ! Additionally, the equation for sigma_Ncn_i_n is:
    !
    ! sigma_Ncn_i_n = sqrt( ln( 1 + ( sigma_Ncn_i^2 / mu_Ncn_i^2 ) ) );
    !
    ! and the equation for rho_chi_Ncn_i_n is:
    !
    ! rho_chi_Ncn_i_n
    ! = rho_chi_Ncn_i * sqrt( exp{ sigma_Ncn_i_n^2 } - 1 ) / sigma_Ncn_i_n.
    !
    ! The product of rho_chi_Ncn_i_n and sigma_Ncn_i_n is:
    !
    ! rho_chi_Ncn_i_n * sigma_Ncn_i_n
    ! = rho_chi_Ncn_i * sqrt( exp{ sigma_Ncn_i_n^2 } - 1 ).
    !
    ! After substituting for sigma_Ncn_i_n^2, the equation for the product of
    ! rho_chi_Ncn_i_n and sigma_Ncn_i_n is:
    !
    ! rho_chi_Ncn_i_n * sigma_Ncn_i_n
    ! = rho_chi_Ncn_i * sqrt( sigma_Ncn_i^2 / mu_Ncn_i^2 );
    !
    ! which can be rewritten as:
    !
    ! rho_chi_Ncn_i_n * sigma_Ncn_i_n
    ! = rho_chi_Ncn * sqrt( <Ncn'^2> / <Ncn>^2 ).
    !
    ! Substituting all of this into the equation for <Nc>, the equation for <Nc>
    ! becomes:
    !
    ! <Nc> = <Ncn>
    !        * SUM(i=1,n) mixt_frac_i
    !          ---
    !          | (1/2) * erfc( - ( 1 / sqrt(2) )
    !          |                 * ( ( mu_chi_i / sigma_chi_i )
    !          |                     + rho_chi_Ncn * sqrt(<Ncn'^2>/<Ncn>^2) ) );
    !          | where sigma_chi_i > 0 and <Ncn'^2> > 0;
    !          |
    !        * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !          | where sigma_chi_i > 0 and <Ncn'^2> = 0;
    !          |
    !          | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
    !          |
    !          | 0; where sigma_chi_i = 0 and mu_chi_i <= 0.
    !          ---
    !
    ! In order to isolate <Ncn>, the value of <Ncn'^2>/<Ncn>^2 is set to a
    ! constant value, const_Ncn.  The value of this constant does not depend on
    ! <Ncn>.  Likewise, the value of rho_chi_Ncn does not depend on <Ncn>.
    ! Solving for <Ncn>, the equation becomes:
    !
    ! <Ncn>
    ! = <Nc> / ( SUM(i=1,n) mixt_frac_i
    !              ---
    !              | (1/2) * erfc( - ( 1 / sqrt(2) )
    !              |                 * ( ( mu_chi_i / sigma_chi_i )
    !              |                     + rho_chi_Ncn * sqrt( const_Ncn ) ) );
    !              | where sigma_chi_i > 0 and const_Ncn > 0;
    !              |
    !            * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !              | where sigma_chi_i > 0 and const_Ncn = 0;
    !              |
    !              | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
    !              |
    !              | 0; where sigma_chi_i = 0 and mu_chi_i <= 0               ).
    !              ---
    !
    ! When the denominator term is 0, there is only clear air.  Both the
    ! numerator (<Nc>) and the denominator have a value of 0, and <Ncn> is set
    ! to an appropriate value.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_1,    & ! Mean of chi (old s) (1st PDF component)           [kg/kg]
      mu_chi_2,    & ! Mean of chi (old s) (2nd PDF component)           [kg/kg]
      sigma_chi_1, & ! Standard deviation of chi (1st PDF component)     [kg/kg]
      sigma_chi_2, & ! Standard deviation of chi (2nd PDF component)     [kg/kg]
      mixt_frac      ! Mixture fraction                                      [-]

    real( kind = core_rknd ), intent(in) :: &
      Ncm,                  & ! Mean cloud droplet conc. (overall)      [num/kg]
      const_Ncnp2_on_Ncnm2, & ! Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
      const_corr_chi_Ncn,   & ! Prescribed correlation of chi and Ncn        [-]
      Ncnm_val_denom_0        ! Ncnm value -- denominator in eqn. is 0  [num/kg]

    ! Return Variable
    real( kind = core_rknd ) :: &
      Ncnm    ! Mean simplified cloud nuclei concentration (overall)    [num/kg]

    ! Local Variable
    real( kind = core_rknd ) :: &
      denominator_term    ! Denominator in the equation for <Ncn>            [-]


    denominator_term &
    = mixt_frac &
      * bivar_Ncnm_eqn_comp( mu_chi_1, sigma_chi_1, &
                             const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn ) &
      + ( one - mixt_frac ) &
        * bivar_Ncnm_eqn_comp( mu_chi_2, sigma_chi_2, &
                               const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn )


    if ( denominator_term > zero ) then

       Ncnm = Ncm / denominator_term

    else ! denominator_term = 0

       ! When the denominator is 0, it is usually because there is only clear
       ! air.  In that scenario, Ncm should also be 0.  Set Ncnm to a value that
       ! is usual or typical
       Ncnm = Ncnm_val_denom_0

    endif ! denominator_term > 0


    return

  end function Ncm_to_Ncnm

  !=============================================================================
  elemental function bivar_NL_chi_Ncn_mean( mu_chi_i, mu_Ncn_i, sigma_chi_i, &
                                            sigma_Ncn_i, sigma_Ncn_i_n, &
                                            corr_chi_Ncn_i_n )

    ! Description:
    ! The double integral over Ncn * H(chi) multiplied by the
    ! bivariate normal-lognormal joint PDF of chi and Ncn is evaluated.  The
    ! integral is given by:
    !
    ! INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;
    !
    ! which reduces to:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi;
    !
    ! where the individual marginal distribution of chi is normal in the ith PDF
    ! component and the individual marginal distribution of Ncn is lognormal in
    ! the ith PDF component.
    !
    ! When both chi and Ncn vary in the ith PDF component, the integral is
    ! evaluated and the result is:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = (1/2) * exp{ mu_Ncn_i_n + (1/2) * sigma_Ncn_i_n^2 }
    !   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
    !                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) );
    !
    ! which can be reduced to:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = (1/2) * mu_Ncn_i
    !   * erfc( - ( 1 / sqrt(2) ) * ( ( mu_chi_i / sigma_chi_i )
    !                                 + rho_chi_Ncn_i_n * sigma_Ncn_i_n ) ).
    !
    ! When chi is constant, but Ncn varies, in the ith PDF component, the
    ! integral is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
    !
    ! when mu_chi_i > 0; and
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
    !
    ! when mu_chi_i <= 0.
    !
    ! When chi varies, but Ncn is constant, in the ith PDF component, the
    ! integral is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi
    ! = mu_Ncn_i * (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! When both chi and Ncn are constant in the ith PDF component, the integral
    ! is evaluated and results in:
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = mu_Ncn_i;
    !
    ! when mu_chi_i > 0; and
    !
    ! INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi = 0;
    !
    ! when mu_chi_i <= 0.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        sqrt_2,   & ! Constant(s)
        one,      &
        one_half, &
        zero,     &
        chi_tol,  &
        Ncn_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_i,         & ! Mean of chi (old s) (ith PDF component)      [kg/kg]
      mu_Ncn_i,         & ! Mean of Ncn (ith PDF component)             [num/kg]
      sigma_chi_i,      & ! Standard deviation of chi (ith PDF comp.)    [kg/kg]
      sigma_Ncn_i,      & ! Standard deviation of Ncn (ith PDF comp.)   [num/kg]
      sigma_Ncn_i_n,    & ! Standard deviation of ln Ncn (ith PDF component) [-]
      corr_chi_Ncn_i_n    ! Correlation of chi and ln Ncn (ith PDF comp.)    [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      bivar_NL_chi_Ncn_mean


    if ( sigma_chi_i <= chi_tol .and. sigma_Ncn_i <= Ncn_tol ) then

       ! The ith PDF component variances of both chi and Ncn are 0.

       if ( mu_chi_i > zero ) then

          bivar_NL_chi_Ncn_mean = mu_Ncn_i

       else ! mu_chi_i <= 0

          bivar_NL_chi_Ncn_mean = zero

       endif


    elseif ( sigma_chi_i <= chi_tol ) then

       ! The ith PDF component variance of chi is 0.

       if ( mu_chi_i > zero ) then

          bivar_NL_chi_Ncn_mean = mu_Ncn_i

       else ! mu_chi_i <= 0

          bivar_NL_chi_Ncn_mean = zero

       endif


    elseif ( sigma_Ncn_i <= Ncn_tol ) then

       ! The ith PDF component variance of Ncn is 0.

       bivar_NL_chi_Ncn_mean &
       = mu_Ncn_i * one_half * erfc( - ( mu_chi_i / ( sqrt_2 * sigma_chi_i ) ) )


    else

       ! Both chi and Ncn vary in the ith PDF component. 

       bivar_NL_chi_Ncn_mean &
       = one_half * mu_Ncn_i &
         * erfc( - ( one / sqrt_2 ) &
                   * ( ( mu_chi_i / sigma_chi_i ) &
                       + corr_chi_Ncn_i_n * sigma_Ncn_i_n ) )


    endif


    return

  end function bivar_NL_chi_Ncn_mean

  !=============================================================================
  elemental function bivar_Ncnm_eqn_comp( mu_chi_i, sigma_chi_i, &
                                          const_Ncnp2_on_Ncnm2, &
                                          const_corr_chi_Ncn )

    ! Description:
    ! When <Ncn> is found based on the value of <Nc>, the following equation is
    ! used:
    !
    ! <Ncn>
    ! = <Nc> / ( SUM(i=1,n) mixt_frac_i
    !              ---
    !              | (1/2) * erfc( - ( 1 / sqrt(2) )
    !              |                 * ( ( mu_chi_i / sigma_chi_i )
    !              |                     + rho_chi_Ncn * sqrt( const_Ncn ) ) );
    !              | where sigma_chi_i > 0 and const_Ncn > 0;
    !              |
    !            * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !              | where sigma_chi_i > 0 and const_Ncn = 0;
    !              |
    !              | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
    !              |
    !              | 0; where sigma_chi_i = 0 and mu_chi_i <= 0               ).
    !              ---
    !
    ! In the above equation, const_Ncn = <Ncn'^2> / <Ncn>^2.  It is a constant,
    ! prescribed parameter.  Likewise, rho_chi_Ncn is a parameter that is not
    ! based on the value of <Ncn>.
    !
    ! When the denominator term is 0, there is only clear air.  Both the
    ! numerator (<Nc>) and the denominator have a value of 0, and <Ncn> is set
    ! to an appropriate value.
    !
    ! The contribution of the ith PDF component to the denominator term in the
    ! equation is calculated here.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        sqrt_2,   & ! Constant(s)
        one,      &
        one_half, &
        zero,     &
        chi_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_i,    & ! Mean of chi (old s) (ith PDF component)           [kg/kg]
      sigma_chi_i    ! Standard deviation of chi (ith PDF component)     [kg/kg]

    real( kind = core_rknd ), intent(in) :: &
      const_Ncnp2_on_Ncnm2, & ! Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
      const_corr_chi_Ncn      ! Prescribed correlation of chi and Ncn        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      bivar_Ncnm_eqn_comp


    if ( sigma_chi_i <= chi_tol ) then

       ! The ith PDF component variances of chi is 0.  The value of the ith PDF
       ! component variance of Ncn does not matter in this scenario.

       if ( mu_chi_i > zero ) then

          bivar_Ncnm_eqn_comp = one

       else ! mu_chi_i <= 0

          bivar_Ncnm_eqn_comp = zero

       endif

    else

       ! Both chi and Ncn vary in the ith PDF component. 

       bivar_Ncnm_eqn_comp &
       = one_half * erfc( - ( one / sqrt_2 ) * ( ( mu_chi_i / sigma_chi_i ) &
                       + const_corr_chi_Ncn * sqrt( const_Ncnp2_on_Ncnm2 ) ) )


    endif


    return

  end function bivar_Ncnm_eqn_comp

!===============================================================================

end module Nc_Ncn_eqns
