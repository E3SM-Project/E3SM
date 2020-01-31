! $Id$
module diagnose_correlations_module 

  use clubb_precision, only: &
      core_rknd

  implicit none 

  public :: diagnose_KK_corr, diagnose_LH_corr, &
            calc_mean, calc_varnce, calc_w_corr
            

  private :: diagnose_corr 


  contains 

!-----------------------------------------------------------------------
  subroutine diagnose_KK_corr( Ncm, rrainm, Nrm, & ! intent(in)
                               Ncp2_on_Ncm2, rrp2_on_rrm2, Nrp2_on_Nrm2, &
                               corr_ws, corr_wrr, corr_wNr, corr_wNc, &
                               pdf_params, &
                               corr_rrNr_p, corr_srr_p, corr_sNr_p, corr_sNc_p, &
                               corr_rrNr, corr_srr, corr_sNr, corr_sNc ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it 
    !   into KK microphysics.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB-Trac:ticket:514)
    !-----------------------------------------------------------------------


    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Type

    use constants_clubb, only: &
        w_tol,         & ! [m/s]
        s_mellor_tol,  & ! [kg/kg]
        Nc_tol,        & ! [num/kg]
        rr_tol,        & ! [kg/kg] 
        Nr_tol           ! [num/kg]

    use stats_type, only: & 
        stat_update_var_pt  ! Procedure(s)

    implicit none

    intrinsic :: sqrt

    ! Local Constants
    integer, parameter :: &
      n_variables = 5

    ! Input Variables

    real( kind = core_rknd ), intent(in) :: &
      Ncm,            &  ! Cloud droplet number conc.            [num/kg]
      rrainm,         &  ! rain water mixing ratio               [kg/kg]
      Nrm,            &  ! Mean rain drop concentration          [num/kg]
      Ncp2_on_Ncm2,   &  ! Variance of Nc divided by Ncm^2       [-]
      rrp2_on_rrm2,   &  ! Variance of rrain divided by rrainm^2 [-]
      Nrp2_on_Nrm2,   &  ! Variance of Nr divided by Nrm^2       [-]
      corr_ws,        &  ! Correlation between s_mellor and w    [-]
      corr_wrr,       &  ! Correlation between rrain and w       [-]
      corr_wNr,       &  ! Correlation between Nr and w          [-]
      corr_wNc,       &  ! Correlation between Nc and w          [-]
      corr_rrNr_p,    &  ! Prescribed correlation between rrain and Nr [-]
      corr_srr_p,     &  ! Prescribed correlation between s and rrain  [-]
      corr_sNr_p,     &  ! Prescribed correlation between s and Nr     [-]
      corr_sNc_p         ! Prescribed correlation between s and Nc     [-]
      
    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters  [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout) :: &
      corr_rrNr,   &  ! Correlation between rrain and Nr [-]
      corr_srr,    &  ! Correlation between s and rrain  [-]
      corr_sNr,    &  ! Correlation between s and Nr     [-]
      corr_sNc        ! Correlation between s and Nc     [-]



    ! Local Variables
    real( kind = core_rknd ), dimension(n_variables, n_variables) :: &
      corr_matrix_approx, & ! [-]
      corr_matrix_prescribed ! [-]

    real( kind = core_rknd ), dimension(n_variables) :: &
      sqrt_xp2_on_xm2, & ! sqrt of x_variance / x_mean^2          [units vary]
      xm                 ! means of the hydrometeors              [units vary]

    ! Indices of the hydrometeors
    integer :: &
      ii_w = 1, &
      ii_s = 2, &
      ii_rrain = 3, &
      ii_Nr = 4, &
      ii_Nc = 5

    integer :: i, j ! Loop Iterators


    !-------------------- Begin code --------------------

    ! set up xp2_on_xm2

    ! TODO Why is wp2_on_wm2=1
    ! S_i is set to 1 for s_mellor and w, because s_mellorm could be 0
    sqrt_xp2_on_xm2(ii_w) = 1._core_rknd
    sqrt_xp2_on_xm2(ii_s) = 1._core_rknd

    sqrt_xp2_on_xm2(ii_rrain) = sqrt(rrp2_on_rrm2)
    sqrt_xp2_on_xm2(ii_Nr) = sqrt(Nrp2_on_Nrm2)
    sqrt_xp2_on_xm2(ii_Nc) = sqrt(Ncp2_on_Ncm2)

    ! initialize the correlation matrix with 0
    do i=1, n_variables
       do j=1, n_variables
          corr_matrix_approx(i,j) = 0._core_rknd
          corr_matrix_prescribed(i,j) = 0._core_rknd
       end do
    end do

    ! set diagonal of the correlation matrix to 1
    do i = 1, n_variables
       corr_matrix_approx(i,i) = 1._core_rknd
       corr_matrix_prescribed(i,i) = 1._core_rknd
    end do


    ! set the first row to the corresponding prescribed correlations
    corr_matrix_approx(ii_s,1) = corr_ws
    corr_matrix_approx(ii_rrain,1) = corr_wrr
    corr_matrix_approx(ii_Nr,1) = corr_wNr
    corr_matrix_approx(ii_Nc,1) = corr_wNc

    !corr_matrix_prescribed = corr_matrix_approx

    ! set up the prescribed correlation matrix
    if( ii_rrain > ii_Nr ) then
      corr_matrix_prescribed(ii_rrain, ii_Nr) = corr_rrNr_p
    else
      corr_matrix_prescribed(ii_Nr, ii_rrain) = corr_rrNr_p
    end if

    if ( ii_s > ii_rrain ) then
      corr_matrix_prescribed(ii_s, ii_rrain) = corr_srr_p
    else
      corr_matrix_prescribed(ii_rrain, ii_s) = corr_srr_p
    end if

    if ( ii_s > ii_Nr ) then
      corr_matrix_prescribed(ii_s, ii_Nr) = corr_sNr_p
    else
      corr_matrix_prescribed(ii_Nr, ii_s) = corr_sNr_p
    end if

    if ( ii_s > ii_Nc ) then
      corr_matrix_prescribed(ii_s, ii_Nc) = corr_sNc_p
    else
      corr_matrix_prescribed(ii_Nc, ii_s) = corr_sNc_p
    end if
    
    call diagnose_corr( n_variables, sqrt_xp2_on_xm2, corr_matrix_prescribed, & !intent(in)
                        corr_matrix_approx ) ! intent(inout)    

    if( ii_rrain > ii_Nr ) then
      corr_rrNr = corr_matrix_approx(ii_rrain, ii_Nr)
    else
      corr_rrNr = corr_matrix_approx(ii_Nr, ii_rrain)
    end if

    if ( ii_s > ii_rrain ) then
      corr_srr = corr_matrix_approx(ii_s, ii_rrain)
    else
      corr_srr = corr_matrix_approx(ii_rrain, ii_s)
    end if

    if ( ii_s > ii_Nr ) then
      corr_sNr = corr_matrix_approx(ii_s, ii_Nr)
    else
      corr_sNr = corr_matrix_approx(ii_Nr, ii_s)
    end if

    if ( ii_s > ii_Nc ) then
      corr_sNc = corr_matrix_approx(ii_s, ii_Nc)
    else
      corr_sNc = corr_matrix_approx(ii_Nc, ii_s)
    end if

  end subroutine diagnose_KK_corr

!-----------------------------------------------------------------------
  subroutine diagnose_LH_corr( xp2_on_xm2, d_variables, corr_matrix_prescribed, & !intent(in)
                               corr_array ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it 
    !   into SILHS microphysics.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use corr_matrix_module, only: &
      iiLH_w ! Variable(s)

    implicit none

    intrinsic :: max, sqrt, transpose

    ! Input Variables
    integer, intent(in) :: d_variables

    real( kind = core_rknd ), dimension(d_variables, d_variables), intent(in) :: &
      corr_matrix_prescribed

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      xp2_on_xm2 ! ratios of x_variance over x_mean^2

    ! Input/Output variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), intent(inout) :: &
      corr_array

    ! Local Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables) :: &
      corr_matrix_pre_swapped

    real( kind = core_rknd ), dimension(d_variables) :: &
      swap_array

    !-------------------- Begin code --------------------

    ! Swap the w-correlations to the first row
    swap_array = corr_array(:, 1)
    corr_array(1:iiLH_w, 1) = corr_array(iiLH_w, iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, 1) = corr_array((iiLH_w+1):d_variables, iiLH_w)
    corr_array(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

    corr_matrix_pre_swapped = corr_matrix_prescribed
    swap_array = corr_matrix_pre_swapped (:,1)
    corr_matrix_pre_swapped(1:iiLH_w, 1) = corr_matrix_pre_swapped(iiLH_w, iiLH_w:1:-1)
    corr_matrix_pre_swapped((iiLH_w+1):d_variables, 1) = corr_matrix_pre_swapped( &
                                                         (iiLH_w+1):d_variables, iiLH_w)
    corr_matrix_pre_swapped(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_matrix_pre_swapped((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

    ! diagnose correlations
    call diagnose_corr( d_variables, sqrt(xp2_on_xm2), corr_matrix_pre_swapped, &
                        corr_array)

    ! Swap rows back
    swap_array = corr_array(:, 1)
    corr_array(1:iiLH_w, 1) = corr_array(iiLH_w, iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, 1) = corr_array((iiLH_w+1):d_variables, iiLH_w)
    corr_array(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

  end subroutine diagnose_LH_corr

!-----------------------------------------------------------------------
  subroutine diagnose_corr( n_variables, sqrt_xp2_on_xm2, corr_matrix_prescribed, & !intent(in)
                            corr_matrix_approx ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix for each timestep.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use parameters_tunable, only:  & 
        alpha_corr ! Constant(s)

    use constants_clubb, only: &
      max_mag_correlation

    implicit none

    intrinsic :: &
      sqrt, abs, sign

    ! Input Variables
    integer, intent(in) :: &
      n_variables  ! number of variables in the correlation matrix [-]
    
    real( kind = core_rknd ), dimension(n_variables), intent(in) :: & 
      sqrt_xp2_on_xm2    ! sqrt of x_variance / x_mean^2 [units vary]

    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(in) :: &
      corr_matrix_prescribed ! correlation matrix [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(inout) :: &
      corr_matrix_approx ! correlation matrix [-]


    ! Local Variables
    integer :: i, j ! Loop iterator

    real( kind = core_rknd ) :: &
      f_ij, &
      f_ij_o

    real( kind = core_rknd ), dimension(n_variables) :: &
      s_1j ! s_1j = sqrt(1-c_1j^2)


    !-------------------- Begin code --------------------

    ! calculate all square roots
    do i = 1, n_variables

       s_1j(i) = sqrt(1._core_rknd-corr_matrix_approx(i,1)**2)

    end do


    ! Diagnose the missing correlations (upper triangle)
    do j = 2, (n_variables-1)
      do i = (j+1), n_variables

        ! formula (16) in the ref. paper (Larson et al. (2011))
        !f_ij = alpha_corr * sqrt_xp2_on_xm2(i) * sqrt_xp2_on_xm2(j) &
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

    calc_varnce = mixt_frac * ((x1 - xm)**2 + x1p2) + (1.0_core_rknd - mixt_frac) * ((x2 - xm)**2 + x2p2)

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

end module diagnose_correlations_module 
