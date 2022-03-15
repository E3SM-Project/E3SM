!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module diagnose_correlations_module 

  use clubb_precision, only: &
      core_rknd

  implicit none 

  public :: calc_mean, calc_varnce, calc_w_corr, &
            calc_cholesky_corr_mtx_approx, &
            cholesky_to_corr_mtx_approx, setup_corr_cholesky_mtx, &
            diagnose_correlations
            

  private :: diagnose_corr, rearrange_corr_array, &
             corr_array_assertion_checks

  private ! Default scope
  contains 

!-----------------------------------------------------------------------
  subroutine diagnose_correlations( pdf_dim, corr_array_pre, & ! Intent(in)
                                    l_calc_w_corr, & ! Intent(in)
                                    corr_array )                   ! Intent(out)
    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it
    !   into SILHS microphysics.

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

!    use array_index, only: &
!        iiPDF_w ! Variable(s)

    use constants_clubb, only: &
        zero

    implicit none

    intrinsic :: max, sqrt, transpose

    ! Input Variables
    integer, intent(in) :: &
      pdf_dim  ! number of diagnosed correlations

    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim), intent(in) :: &
      corr_array_pre   ! Prescribed correlations

    logical, intent(in) :: &
      l_calc_w_corr ! Calculate the correlations between w and the hydrometeors

    ! Output variables
    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim), intent(out) :: &
      corr_array

    ! Local Variables
    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim) :: &
      corr_array_pre_swapped, &
      corr_array_swapped

    ! We actually don't need this right now
    real( kind = core_rknd ), dimension(pdf_dim) :: &
      sigma2_on_mu2_ip_array  ! Ratios: sigma_x^2/mu_x^2 (ith PDF comp.) ip [-]

    integer :: i ! Loop iterator

    !-------------------- Begin code --------------------

    ! Initialize sigma2_on_mu2_ip_array
    do i = 1, pdf_dim
       sigma2_on_mu2_ip_array(i) = zero
    end do 

    ! Swap the w-correlations to the first row for the prescribed correlations
    call rearrange_corr_array( pdf_dim, corr_array_pre, & ! Intent(in)
                               corr_array_pre_swapped)        ! Intent(inout)

    ! diagnose correlations
    
    if ( .not. l_calc_w_corr ) then
       corr_array_swapped = corr_array_pre_swapped
    endif

    call diagnose_corr( pdf_dim, sqrt(sigma2_on_mu2_ip_array), & ! intent(in)
                        corr_array_pre_swapped, & ! intent(in)
                        corr_array_swapped ) ! intent(inout)

    ! Swap rows back
    call rearrange_corr_array( pdf_dim, corr_array_swapped, & ! Intent(in)
                               corr_array)        ! Intent(out)

  end subroutine diagnose_correlations


  !-----------------------------------------------------------------------
  subroutine diagnose_corr( n_variables, sqrt_sigma2_on_mu2_ip, & ! intent(in)
                            corr_matrix_prescribed, & !intent(in)
                            corr_matrix_approx ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix for each timestep.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        max_mag_correlation

    implicit none

    intrinsic :: &
      sqrt, abs, sign

    ! Input Variables
    integer, intent(in) :: &
      n_variables  ! number of variables in the correlation matrix [-]
    
    real( kind = core_rknd ), dimension(n_variables), intent(in) :: & 
      sqrt_sigma2_on_mu2_ip  ! sqrt of sigma_x^2/mu_x^2 (ith PDF comp.) ip [-]

    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(in) :: &
      corr_matrix_prescribed ! correlation matrix [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(inout) :: &
      corr_matrix_approx ! correlation matrix [-]


    ! Local Variables
    integer :: i, j ! Loop iterator

    real( kind = core_rknd ) :: &
      f_ij
!      f_ij_o

    real( kind = core_rknd ), dimension(n_variables) :: &
      s_1j ! s_1j = sqrt(1-c_1j^2)


    !-------------------- Begin code --------------------

    ! Remove compiler warnings about unused variables.
    if ( .false. ) then
       print *, "sqrt_sigma2_on_mu2_ip = ", sqrt_sigma2_on_mu2_ip
    endif

    ! calculate all square roots
    do i = 1, n_variables

       s_1j(i) = sqrt(1._core_rknd-corr_matrix_approx(i,1)**2)

    end do


    ! Diagnose the missing correlations (upper triangle)
    do j = 2, (n_variables-1)
      do i = (j+1), n_variables

        ! formula (16) in the ref. paper (Larson et al. (2011))
        !f_ij = alpha_corr * sqrt_sigma2_on_mu2_ip(i) * sqrt_sigma2_on_mu2_ip(j) &
        !        * sign(1.0_core_rknd,corr_matrix_approx(1,i)*corr_matrix_approx(1,j))

        ! If the predicting c1i's are small then cij will be closer to the prescribed value. If
        ! the c1i's are bigger, then cij will be closer to formular (15) from the ref. paper. See
        ! clubb:ticket:514:comment:61 for details.
      !f_ij = (1-abs(corr_matrix_approx(1,i)*corr_matrix_approx(1,j)))*corr_matrix_prescribed(i,j) &
      !       + abs(corr_matrix_approx(1,i)*corr_matrix_approx(1,j))*f_ij_o

        f_ij = corr_matrix_prescribed(i,j)

        ! make sure -1 < f_ij < 1
        if ( f_ij < -max_mag_correlation ) then

           f_ij = -max_mag_correlation

        else if ( f_ij > max_mag_correlation ) then

           f_ij = max_mag_correlation

        end if


        ! formula (15) in the ref. paper (Larson et al. (2011))
        corr_matrix_approx(i,j) = corr_matrix_approx(i,1) * corr_matrix_approx(j,1) &
        + f_ij * s_1j(i) * s_1j(j)

      end do ! do j
    end do ! do i
    
  end subroutine diagnose_corr 


  !-----------------------------------------------------------------------
!  subroutine approx_w_corr( nz, pdf_dim, pdf_params, & ! Intent(in)
!                            rrm, Nrm, Ncnm, &
!                            stdev_w, sigma_rr_1, &
!                            sigma_Nr_1, sigma_Ncn_1, &
!                            corr_array) ! Intent(out)
!    ! Description:
!    ! Approximate the correlations of w with the hydrometeors.
!
!    ! References:
!    ! clubb:ticket:514
!    !-----------------------------------------------------------------------
!
!    use clubb_precision, only: &
!        core_rknd ! Variable(s)
!
!    use pdf_parameter_module, only:  &
!        pdf_parameter  ! Type
!
!    use constants_clubb, only:  &
!        one,          & ! Constant(s)
!        rr_tol,       &
!        Nr_tol,       &
!        Ncn_tol,      &
!        w_tol,        & ! [m/s]
!        chi_tol    ! [kg/kg]
!
!    implicit none
!
!    ! Input Variables
!    integer, intent(in) :: &
!      pdf_dim, & ! Number of diagnosed correlations
!      nz             ! Number of model vertical grid levels
!
!    type(pdf_parameter), dimension(nz), intent(in) :: &
!      pdf_params    ! PDF parameters                         [units vary]
!
!    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
!      rrm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
!      Nrm,             & ! Mean rain drop concentration, < N_r >    [num/kg]
!      Ncnm,            & ! Mean cloud nuclei conc., < N_cn >        [num/kg]
!      stdev_w            ! Standard deviation of w                              [m/s]
!
!    real( kind = core_rknd ), intent(in) :: &
!      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
!      sigma_Nr_1,    & ! Standard deviation of Nr (2nd PDF component)   [num/kg]
!      sigma_rr_1       ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
!
!    ! Output Variables
!    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim, nz), intent(out) :: &
!      corr_array
!
!    ! Local Variables
!    real( kind = core_rknd ), dimension(nz) :: &
!      corr_chi_w,       & ! Correlation between w and chi(s_mellor) (both components)       [-]
!      corr_wrr,      & ! Correlation between w and rr (both components)      [-]
!      corr_wNr,      & ! Correlation between w and Nr (both components)      [-]
!      corr_wNcn        ! Correlation between w and Ncn (both components)     [-]
!
!    real( kind = core_rknd ), dimension(nz) ::  &
!      wpchip_zt,   & ! Covariance of chi and w on the zt-grid    [(m/s)(kg/kg)]
!      wprrp_zt,  & ! Covariance of r_r and w on the zt-grid  [(m/s)(kg/kg)]
!      wpNrp_zt,  & ! Covariance of N_r and w on the zt-grid  [(m/s)(#/kg)]
!      wpNcnp_zt    ! Covariance of N_cn and w on the zt-grid  [(m/s)(#/kg)]
!
!    real( kind = core_rknd ) :: &
!      chi_m,      & ! Mean of chi (s_mellor)                             [kg/kg]
!      stdev_chi     ! Standard deviation of chi (s_mellor)               [kg/kg]
!
!    integer :: k ! vertical loop iterator
!
!    ! ----- Begin Code -----
!
!    call approx_w_covar( nz, pdf_params, rrm, Nrm, Ncnm, & ! Intent(in)
!                         wpchip_zt, wprrp_zt, wpNrp_zt, wpNcnp_zt ) ! Intent(out)
!
!    do k = 1, nz
!
!       chi_m &
!       = calc_mean( pdf_params(k)%mixt_frac, pdf_params(k)%chi_1, &
!                    pdf_params(k)%chi_2 )
!
!       stdev_chi &
!       = sqrt( pdf_params(k)%mixt_frac &
!               * ( ( pdf_params(k)%chi_1 - chi_m )**2 &
!                     + pdf_params(k)%stdev_chi_1**2 ) &
!             + ( one - pdf_params(k)%mixt_frac ) &
!               * ( ( pdf_params(k)%chi_2 - chi_m )**2 &
!                     + pdf_params(k)%stdev_chi_2**2 ) &
!             )
!
!       corr_chi_w(k) &
!       = calc_w_corr( wpchip_zt(k), stdev_w(k), stdev_chi, &
!                      w_tol, chi_tol )
!
!       corr_wrr(k) &
!       = calc_w_corr( wprrp_zt(k), stdev_w(k), sigma_rr_1, w_tol, rr_tol )
!
!       corr_wNr(k) &
!       = calc_w_corr( wpNrp_zt(k), stdev_w(k), sigma_Nr_1, w_tol, Nr_tol )
!
!       corr_wNcn(k) &
!       = calc_w_corr( wpNcnp_zt(k), stdev_w(k), sigma_Ncn_1, w_tol, Ncn_tol )
!
!    enddo
!
!    call set_w_corr( nz, pdf_dim, & ! Intent(in)
!                         corr_chi_w, corr_wrr, corr_wNr, corr_wNcn, &
!                         corr_array ) ! Intent(inout)
!
!  end subroutine approx_w_corr


  !-----------------------------------------------------------------------
!  subroutine approx_w_covar( nz, pdf_params, rrm, Nrm, Ncnm, Kh_zm, &   ! Intent(in)
!                             wpchip_zt, wprrp_zt, wpNrp_zt, wpNcnp_zt ) ! Intent(out)
!    ! Description:
!    ! Approximate the covariances of w with the hydrometeors using Eddy
!    ! diffusivity.
!
!    ! References:
!    ! clubb:ticket:514
!    !-----------------------------------------------------------------------
!
!    use clubb_precision, only: &
!      core_rknd ! Variable(s)
!
!    use grid_class, only: &
!        gr,  & ! Variable(s)
!        zm2zt,  & ! Procedure(s)
!        zt2zm
!
!    use pdf_parameter_module, only:  &
!        pdf_parameter  ! Type
!
!    use constants_clubb, only: &
!        one ! Constant(s)
!
!    use advance_windm_edsclrm_module, only: &
!        xpwp_fnc ! Procedure(s)
!
!    implicit none
!
!    ! Input Variables
!    integer, intent(in) :: &
!      nz          ! Number of model vertical grid levels
!
!    type(pdf_parameter), dimension(nz), intent(in) :: &
!      pdf_params    ! PDF parameters                         [units vary]
!
!    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
!      rrm,   & ! Mean rain water mixing ratio, < r_r >      [kg/kg]
!      Nrm,   & ! Mean rain drop concentration, < N_r >      [num/kg]
!      Ncnm,  & ! Mean cloud nuclei concentration, < N_cn >  [num/kg]
!      Kh_zm    ! Eddy diffusivity coef. on momentum levels  [m^2/s]
!
!    ! Output Variables
!    real( kind = core_rknd ), dimension(nz), intent(out) ::  &
!      wpchip_zt,   & ! Covariance of chi(s) and w on the zt-grid     [(m/s)(kg/kg)]
!      wprrp_zt,  & ! Covariance of r_r and w on the zt-grid   [(m/s)(kg/kg)]
!      wpNrp_zt,  & ! Covariance of N_r and w on the zt-grid   [(m/s)(#/kg)]
!      wpNcnp_zt    ! Covariance of N_cn and w on the zt-grid  [(m/s)(#/kg)]
!
!    ! Local Variables
!    real( kind = core_rknd ), dimension(nz) ::  &
!      wpchip_zm,   & ! Covariance of chi(s) and w on the zm-grid     [(m/s)(kg/kg)]
!      wprrp_zm,  & ! Covariance of r_r and w on the zm-grid   [(m/s)(kg/kg)]
!      wpNrp_zm,  & ! Covariance of N_r and w on the zm-grid   [(m/s)(#/kg)]
!      wpNcnp_zm    ! Covariance of N_cn and w on the zm-grid  [(m/s)(#/kg)]
!
!    integer :: k ! vertical loop iterator
!
!    ! ----- Begin Code -----
!
!    ! calculate the covariances of w with the hydrometeors
!    do k = 1, nz
!      wpchip_zm(k) = pdf_params(k)%mixt_frac &
!                   * ( one - pdf_params(k)%mixt_frac ) &
!                   * ( pdf_params(k)%chi_1 - pdf_params(k)%chi_2 ) &
!                   * ( pdf_params(k)%w_1 - pdf_params(k)%w_2 )
!    enddo
!
!! same for wpNrp
!!    wprrp_zm(1:nz-1) &
!!    = xpwp_fnc( -c_K_hm * Kh_zm(1:nz-1), &
!!                rrm(1:nz-1) / max( precip_frac(1:nz-1), eps ), &
!!                rrm(2:nz) / max( precip_frac(2:nz), eps ), &
!!                gr%invrs_dzm(1:nz-1) )
!
!    wprrp_zm(1:nz-1) &
!    = xpwp_fnc( -c_K_hm * Kh_zm(1:nz-1), &
!                rrm(1:nz-1), rrm(2:nz), &
!                gr%invrs_dzm(1:nz-1) )
!
!    wpNrp_zm(1:nz-1) &
!    = xpwp_fnc( -c_K_hm * Kh_zm(1:nz-1), &
!                Nrm(1:nz-1), Nrm(2:nz), &
!                gr%invrs_dzm(1:nz-1) )
!
!    wpNcnp_zm(1:nz-1) = xpwp_fnc( -c_K_hm * Kh_zm(1:nz-1), Ncnm(1:nz-1), &
!                                  Ncnm(2:nz), gr%invrs_dzm(1:nz-1) )
!
!    ! Boundary conditions; We are assuming constant flux at the top.
!    wprrp_zm(nz)  = wprrp_zm(nz-1)
!    wpNrp_zm(nz)  = wpNrp_zm(nz-1)
!    wpNcnp_zm(nz) = wpNcnp_zm(nz-1)
!
!    ! interpolate back to zt-grid
!    wpchip_zt   = zm2zt(wpchip_zm)
!    wprrp_zt  = zm2zt(wprrp_zm)
!    wpNrp_zt  = zm2zt(wpNrp_zm)
!    wpNcnp_zt = zm2zt(wpNcnp_zm)
!
!  end subroutine approx_w_covar

  !-----------------------------------------------------------------------
  function calc_w_corr( wpxp, stdev_w, stdev_x, w_tol, x_tol )
    ! Description:
    ! Compute the correlations of w with the hydrometeors.

    ! References:
    ! clubb:ticket:514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        max_mag_correlation

    implicit none

    intrinsic :: max

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      stdev_w,  & ! standard deviation of w [m/s]
      stdev_x,  & ! standard deviation of x [units vary]
      wpxp,     & ! Covariances of w with the hydrometeors [units vary]
      w_tol,    & ! tolerance for w [m/s]
      x_tol       ! tolerance for x [units vary]

    real( kind = core_rknd ) :: &
      calc_w_corr

    ! --- Begin Code ---

    calc_w_corr = wpxp / ( max(stdev_x, x_tol) * max(stdev_w, w_tol) )

    ! Make sure the correlation is in [-1,1]
    if ( calc_w_corr < -max_mag_correlation ) then

      calc_w_corr = -max_mag_correlation

    else if ( calc_w_corr > max_mag_correlation ) then

      calc_w_corr = max_mag_correlation

    end if
   
  end function calc_w_corr


  !-----------------------------------------------------------------------
  function calc_varnce( mixt_frac, x1, x2, xm, x1p2, x2p2 )

    ! Description:
    ! Calculate the variance xp2 from the components x1, x2.

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
    !   page 3535
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mixt_frac, &  ! mixing ratio [-]
      x1, &         ! first component of the double gaussian [units vary]
      x2, &         ! second component of the double gaussian [units vary]
      xm, &         ! mean of x [units vary]
      x1p2, &       ! variance of the first component [units vary]
      x2p2          ! variance of the second component [units vary]

    ! Return Variable
    real( kind = core_rknd ) :: &
      calc_varnce ! variance of x (both components) [units vary]

    ! --- Begin Code ---

    calc_varnce &
    = mixt_frac * ( ( x1 - xm )**2 + x1p2 ) &
      + ( 1.0_core_rknd - mixt_frac ) * ( ( x2 - xm )**2 + x2p2 )

    return
  end function calc_varnce

  !-----------------------------------------------------------------------
  function calc_mean( mixt_frac, x1, x2 )

    ! Description:
    ! Calculate the mean xm from the components x1, x2.

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
    !   page 3535
    !-----------------------------------------------------------------------
    
    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mixt_frac, &  ! mixing ratio [-]
      x1, &         ! first component of the double gaussian [units vary]
      x2            ! second component of the double gaussian [units vary]

    ! Return Variable
    real( kind = core_rknd ) :: &
      calc_mean  ! mean of x (both components) [units vary]

    ! --- Begin Code ---

    calc_mean = mixt_frac * x1 + (1.0_core_rknd - mixt_frac) * x2

    return
  end function calc_mean


  !-----------------------------------------------------------------------
  subroutine calc_cholesky_corr_mtx_approx &
                     ( n_variables, corr_matrix, & ! intent(in)
                       corr_cholesky_mtx, corr_mtx_approx ) ! intent(out)

    ! Description:
    !   This subroutine calculates the transposed correlation cholesky matrix
    !   from the correlation matrix
    !
    ! References:
    !   1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   2 CLUBB Trac ticket#514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        zero ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_variables  ! number of variables in the correlation matrix [-]

    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(in) :: &
      corr_matrix ! correlation matrix [-]

    ! Output Variables

    ! correlation cholesky matrix transposed L', C = LL'; see reference 1 formula 10
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(out) :: &
      corr_cholesky_mtx, & ! Transposed correlation cholesky matrix    [-]
      corr_mtx_approx      ! Approximated correlation matrix (C = LL') [-]

    ! Local Variables
    integer :: i, j ! Loop iterators

    ! Swapped means that the w-correlations are swapped to the first row
    real( kind = core_rknd ), dimension(n_variables,n_variables) :: &
      corr_cholesky_mtx_swap, & ! Swapped correlation cholesky matrix    [-]
      corr_mtx_approx_swap,   & ! Swapped correlation matrix (approx.)   [-]
      corr_mtx_swap             ! Swapped correlation matrix             [-]

    !-------------------- Begin code --------------------

    call rearrange_corr_array( n_variables, corr_matrix, & ! Intent(in)
                               corr_mtx_swap )             ! Intent(inout)

    call setup_corr_cholesky_mtx( n_variables, corr_mtx_swap, & ! intent(in)
                                  corr_cholesky_mtx_swap )      ! intent(out)

    call rearrange_corr_array( n_variables, corr_cholesky_mtx_swap, & ! Intent(in)
                               corr_cholesky_mtx )                    ! Intent(inout)

    call cholesky_to_corr_mtx_approx( n_variables, corr_cholesky_mtx_swap, & ! intent(in)
                                      corr_mtx_approx_swap )                 ! intent(out)

    call rearrange_corr_array( n_variables, corr_mtx_approx_swap, &  ! Intent(in)
                               corr_mtx_approx )                     ! Intent(inout)

    call corr_array_assertion_checks( n_variables, corr_mtx_approx ) ! intent(in)

    ! Set lower triangle to zero for conformity
    do i = 2, n_variables
       do j = 1, i-1
          corr_mtx_approx(j,i) = zero
       end do
    end do

    return

  end subroutine calc_cholesky_corr_mtx_approx
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine setup_corr_cholesky_mtx( n_variables, corr_matrix, & ! intent(in)
                                      corr_cholesky_mtx_t )       ! intent(out)

    ! Description:
    !   This subroutine calculates the transposed correlation cholesky matrix
    !   from the correlation matrix
    !
    ! References:
    !   1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   2 CLUBB Trac ticket#514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        zero, & ! Variable(s)
        one

    implicit none

    intrinsic :: sqrt

    ! Input Variables
    integer, intent(in) :: &
      n_variables  ! number of variables in the correlation matrix [-]

    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(in) :: &
      corr_matrix ! correlation matrix [-]

    ! Output Variables

    ! correlation cholesky matrix transposed L', C = LL'; see reference 1 formula 10
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(out) :: &
      corr_cholesky_mtx_t ! transposed correlation cholesky matrix [-]

    ! Local Variables
    integer :: i, j, k ! Loop iterators

    real( kind = core_rknd ), dimension(n_variables, n_variables) :: &
      s ! s(i,j) = sqrt(1-c(i,j)^2); see ref 1

    !-------------------- Begin code --------------------

    ! calculate all necessary square roots
    do i = 1, n_variables-1
       do j = i+1, n_variables

          s(j,i) = sqrt(1._core_rknd - corr_matrix(j,i)**2)

       end do
    end do

    !!! calculate transposed correlation cholesky matrix; ref 1 formula 10

    ! initialize matrix to zero
    do i = 1, n_variables
       do j = 1, n_variables

          corr_cholesky_mtx_t(j,i) = zero

       end do
    end do

    ! initialize upper triangle and diagonal to one
    do i = 1, n_variables
       do j = i, n_variables

          corr_cholesky_mtx_t(j,i) = one

       end do
    end do

    ! set diagonal elements
    do j = 2, n_variables
       do i = 1, j-1

          corr_cholesky_mtx_t(j,j) = corr_cholesky_mtx_t(j,j)*s(j,i)
         ! print *, "s(", j, ",", i, ") = ", s(j,i)

       end do
    end do

    ! set first row
    do j = 2, n_variables

       corr_cholesky_mtx_t(j,1) = corr_matrix(j,1)

    end do

    ! set upper triangle
    do i = 2, n_variables-1
       do j = i+1, n_variables
          do k = 1, i-1

             corr_cholesky_mtx_t(j,i) = corr_cholesky_mtx_t(j,i)*s(j,k)

          end do

          corr_cholesky_mtx_t(j,i) = corr_cholesky_mtx_t(j,i)*corr_matrix(j,i)

       end do
    end do

    return

  end subroutine setup_corr_cholesky_mtx
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  subroutine cholesky_to_corr_mtx_approx( n_variables, corr_cholesky_mtx_t, & ! intent(in)
                                          corr_matrix_approx )                ! intent(out)

    ! Description:
    !   This subroutine approximates the correlation matrix from the correlation
    !   cholesky matrix
    !
    ! References:
    !   1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   2 CLUBB Trac ticket#514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    intrinsic :: matmul, transpose

    ! Input Variables
    integer, intent(in) :: &
      n_variables  ! number of variables in the correlation matrix [-]

    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(in) :: &
      corr_cholesky_mtx_t ! transposed correlation cholesky matrix [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(out) :: &
      corr_matrix_approx ! correlation matrix [-]

    !-------------------- Begin code --------------------

    ! approximate the correlation matrix; see ref 1 formula (8)
    corr_matrix_approx = matmul(corr_cholesky_mtx_t, transpose(corr_cholesky_mtx_t))

    return

  end subroutine cholesky_to_corr_mtx_approx
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  subroutine corr_array_assertion_checks( n_variables, corr_array )

    ! Description:
    !   This subroutine does the assertion checks for the corr_array.

    ! References:
    !
    !
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        max_mag_correlation ! Variable(s)

    use constants_clubb, only: &
        one ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level  ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_variables  ! number of variables in the correlation matrix [-]

    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(in) :: &
      corr_array ! correlation matrix [-]

    ! Local Variables
    integer :: i, j ! Loop iterator

    real( kind = core_rknd ), parameter :: &
    tol = 1.e-6_core_rknd ! Maximum acceptable tolerance for the difference of the diagonal
                          ! elements of corr_array to one

    !-------------------- Begin code --------------------

    if ( clubb_at_least_debug_level( 1 ) ) then

       do i = 1, n_variables - 1
          do j = i+1, n_variables

             ! Check if upper and lower triangle values are within the correlation boundaries
             if ( ( corr_array(i,j) < -max_mag_correlation ) &
                  .or. ( corr_array(i,j) > max_mag_correlation ) &
                  .or. ( corr_array(j,i) < -max_mag_correlation ) &
                  .or. ( corr_array(j,i) > max_mag_correlation ) ) &
             then

                error stop "Error: A value in the correlation matrix is out of range."

             endif

          enddo
       enddo

    endif

    if ( clubb_at_least_debug_level( 2 ) ) then

       do i = 1, n_variables
          ! Check if the diagonal elements are one (up to a tolerance)
          if ( ( corr_array(i,i) > one + tol ) .or. (corr_array(i,i) < one - tol ) ) then

             error stop "Error: Diagonal element(s) of the correlation matrix are unequal to one."

          endif
       enddo

    endif

    return

  end subroutine corr_array_assertion_checks


!-----------------------------------------------------------------------
  subroutine rearrange_corr_array( pdf_dim, corr_array, & ! Intent(in)
                                   corr_array_swapped)        ! Intent(out)
    ! Description:
    !   This subroutine swaps the w-correlations to the first row if the input
    !   matrix is in the same order as the *_corr_array_cloud.in files. It swaps
    !   the rows back to the order of the *_corr_array_cloud.in files if the
    !   input matrix is already swapped (first row w-correlations).
    !
    ! References:
    !
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use array_index, only: &
        iiPDF_w ! Variable(s)

    implicit none

    intrinsic :: max, sqrt, transpose

    ! Input Variables
    integer, intent(in) :: &
      pdf_dim  ! number of diagnosed correlations

    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim), intent(in) :: &
      corr_array   ! Correlation matrix

    ! Output variables
    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim), intent(out) :: &
      corr_array_swapped   ! Swapped correlation matrix

    ! Local Variables
    real( kind = core_rknd ), dimension(pdf_dim) :: &
      swap_array

    !-------------------- Begin code --------------------


    ! Swap the w-correlations to the first row for the prescribed correlations
    corr_array_swapped = corr_array
    swap_array = corr_array_swapped (:,1)
    corr_array_swapped(1:iiPDF_w, 1) = corr_array_swapped(iiPDF_w, iiPDF_w:1:-1)
    corr_array_swapped((iiPDF_w+1):pdf_dim, 1) = corr_array_swapped( &
                                                    (iiPDF_w+1):pdf_dim, iiPDF_w)
    corr_array_swapped(iiPDF_w, 1:iiPDF_w) = swap_array(iiPDF_w:1:-1)
    corr_array_swapped((iiPDF_w+1):pdf_dim, iiPDF_w) = swap_array((iiPDF_w+1):pdf_dim)

    return

  end subroutine rearrange_corr_array
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
!  subroutine set_w_corr( nz, pdf_dim, & ! Intent(in)
!                         corr_chi_w, corr_wrr, corr_wNr, corr_wNcn, &
!                         corr_array ) ! Intent(inout)
!
!    ! Description:
!    ! Set the first row of corr_array to the according w-correlations.
!
!    ! References:
!    ! clubb:ticket:514
!    !-----------------------------------------------------------------------
!
!    use clubb_precision, only: &
!      core_rknd ! Variable(s)
!
!    use array_index, only: &
!      iiPDF_w,           & ! Variable(s)
!      iiPDF_chi,         &
!      iiPDF_rr,          &
!      iiPDF_Nr,          &
!      iiPDF_Ncn
!
!    implicit none
!
!    ! Input Variables
!    integer, intent(in) :: &
!      nz,          & ! Number of model vertical grid levels
!      pdf_dim   ! Number of Variables to be diagnosed
!
!    real( kind = core_rknd ), dimension(nz), intent(in) :: &
!      corr_chi_w,       & ! Correlation between chi (s) & w (both components)         [-]
!      corr_wrr,      & ! Correlation between rr & w (both components)        [-]
!      corr_wNr,      & ! Correlation between Nr & w (both components)        [-]
!      corr_wNcn        ! Correlation between Ncn & w (both components)       [-]
!
!    ! Input/Output Variables
!    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim, nz), &
!    intent(inout) :: &
!      corr_array
!
!    ! ----- Begin Code -----
!
!      corr_array(iiPDF_w, iiPDF_chi, :) = corr_chi_w
!      corr_array(iiPDF_w, iiPDF_rr, :) = corr_wrr
!      corr_array(iiPDF_w, iiPDF_Nr, :) = corr_wNr
!      corr_array(iiPDF_w, iiPDF_Ncn, :) = corr_wNcn
!
!  end subroutine set_w_corr

  !=============================================================================
!  subroutine unpack_correlations( pdf_dim, corr_array, & ! Intent(in)
!                                  corr_w_chi, corr_wrr, corr_wNr, corr_wNcn, &
!                                  corr_chi_eta, corr_chi_rr, corr_chi_Nr, corr_chi_Ncn, &
!                                  corr_eta_rr, corr_eta_Nr, corr_eta_Ncn, corr_rrNr )  
!
!    ! Description:
!
!    ! References:
!    !-----------------------------------------------------------------------

!    use clubb_precision, only: &
!        core_rknd ! Variable(s)

!    use array_index, only: &
!        iiPDF_w,        & ! Variable(s)
!        iiPDF_chi,      &
!        iiPDF_eta,      &
!        iiPDF_rr,       &
!        iiPDF_Nr,       &
!        iiPDF_Ncn

!    implicit none

!    intrinsic :: max, sqrt, transpose

!    ! Input Variables
!    integer, intent(in) :: &
!      pdf_dim  ! number of diagnosed correlations

!    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim), intent(in) :: &
!      corr_array   ! Prescribed correlations

!    ! Output variables
!    real( kind = core_rknd ), intent(out) :: &
!      corr_w_chi,     & ! Correlation between w and chi(s) (1st PDF component)     [-]
!      corr_wrr,    & ! Correlation between w and rr (1st PDF component) ip [-]
!      corr_wNr,    & ! Correlation between w and Nr (1st PDF component) ip [-]
!      corr_wNcn,   & ! Correlation between w and Ncn (1st PDF component)   [-]
!      corr_chi_eta,     & ! Correlation between chi(s) and eta(t) (1st PDF component)     [-]
!      corr_chi_rr,    & ! Correlation between chi(s) and rr (1st PDF component) ip [-]
!      corr_chi_Nr,    & ! Correlation between chi(s) and Nr (1st PDF component) ip [-]
!      corr_chi_Ncn,   & ! Correlation between chi(s) and Ncn (1st PDF component)   [-]
!      corr_eta_rr,    & ! Correlation between eta(t) and rr (1st PDF component) ip [-]
!      corr_eta_Nr,    & ! Correlation between eta(t) and Nr (1st PDF component) ip [-]
!      corr_eta_Ncn,   & ! Correlation between (t) and Ncn (1st PDF component)   [-]
!      corr_rrNr      ! Correlation between rr & Nr (1st PDF component) ip  [-]

    ! ---- Begin Code ----

!    corr_w_chi   = corr_array(iiPDF_w, iiPDF_chi)
!    corr_wrr  = corr_array(iiPDF_w, iiPDF_rr)
!    corr_wNr  = corr_array(iiPDF_w, iiPDF_Nr)
!    corr_wNcn = corr_array(iiPDF_w, iiPDF_Ncn)
!    corr_chi_eta   = corr_array(iiPDF_chi, iiPDF_eta)
!    corr_chi_rr  = corr_array(iiPDF_chi, iiPDF_rr)
!    corr_chi_Nr  = corr_array(iiPDF_chi, iiPDF_Nr)
!    corr_chi_Ncn = corr_array(iiPDF_chi, iiPDF_Ncn)
!    corr_eta_rr  = corr_array(iiPDF_eta, iiPDF_rr)
!    corr_eta_Nr  = corr_array(iiPDF_eta, iiPDF_Nr)
!    corr_eta_Ncn = corr_array(iiPDF_eta, iiPDF_Ncn)
!    corr_rrNr = corr_array(iiPDF_rr, iiPDF_Nr)

!    corr_w_chi   = corr_array(iiPDF_chi, iiPDF_w)
!    corr_wrr  = corr_array(iiPDF_rr, iiPDF_w)
!    corr_wNr  = corr_array(iiPDF_Nr, iiPDF_w)
!    corr_wNcn = corr_array(iiPDF_Ncn, iiPDF_w)
!    corr_chi_eta   = corr_array(iiPDF_eta, iiPDF_chi)
!    corr_chi_rr  = corr_array(iiPDF_rr, iiPDF_chi)
!    corr_chi_Nr  = corr_array(iiPDF_Nr, iiPDF_chi)
!    corr_chi_Ncn = corr_array(iiPDF_Ncn, iiPDF_chi)
!    corr_eta_rr  = corr_array(iiPDF_rr, iiPDF_eta)
!    corr_eta_Nr  = corr_array(iiPDF_Nr, iiPDF_eta)
!    corr_eta_Ncn = corr_array(iiPDF_Ncn, iiPDF_eta)
!    corr_rrNr = corr_array(iiPDF_rr, iiPDF_Nr)

!  end subroutine unpack_correlations

!===============================================================================

end module diagnose_correlations_module
