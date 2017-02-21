!-------------------------------------------------------------------------------
! $Id: latin_hypercube_driver_module.F90 8082 2016-05-05 22:48:48Z raut@uwm.edu $
!===============================================================================
module latin_hypercube_driver_module

  use clubb_precision, only: &
    core_rknd           ! Constant

  implicit none

  ! Constant Parameters
  logical, parameter, private :: &
    l_output_2D_lognormal_dist   = .false., & ! Output a 2D netCDF file of the lognormal variates
    l_output_2D_uniform_dist     = .false.    ! Output a 2D netCDF file of the uniform distribution

  private ! Default scope

  type lh_clipped_variables_type

    real( kind = core_rknd ) :: &
      rt,      & ! Total water mixing ratio            [kg/kg]
      thl,     & ! Liquid potential temperature        [K]
      rc,      & ! Cloud water mixing ratio            [kg/kg]
      rv,      & ! Vapor water mixing ratio            [kg/kg]
      Nc         ! Cloud droplet number concentration  [#/kg]

  end type lh_clipped_variables_type

  public :: lh_clipped_variables_type

#ifdef SILHS
  public :: latin_hypercube_2D_output, &
    latin_hypercube_2D_close, stats_accumulate_lh, lh_subcolumn_generator, &
    copy_X_nl_into_hydromet_all_pts, clip_transform_silhs_output

  private :: stats_accumulate_uniform_lh

  contains

!-------------------------------------------------------------------------------
  subroutine lh_subcolumn_generator &
             ( iter, d_variables, num_samples, sequence_length, nz, & ! In
               pdf_params, delta_zm, rcm, Lscale, & ! In
               rho_ds_zt, mu1, mu2, sigma1, sigma2, & ! In
               corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
               hydromet_pdf_params, & ! In
               X_nl_all_levs, X_mixt_comp_all_levs, & ! Out
               lh_sample_point_weights ) ! Out

! Description:
!   Generate sample points of moisture, temperature, et cetera for the purpose
!   of computing tendencies with a microphysics or radiation scheme.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use corr_varnce_module, only: &
      iiPDF_chi    ! Variables

    use transform_to_pdf_module, only: &
      transform_uniform_sample_to_pdf      ! Procedure

    use output_2D_samples_module, only: &
      output_2D_lognormal_dist_file, & ! Procedure(s)
      output_2D_uniform_dist_file

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use constants_clubb, only: &
      fstderr, & ! Constant(s)
      zero_threshold, &
      zero, &
      one, &
      cloud_frac_min

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use clubb_precision, only: &
      core_rknd, &
      stat_rknd

    use parameters_silhs, only: &
      l_lh_importance_sampling

    implicit none

    ! External
    intrinsic :: allocated, mod, maxloc, epsilon, transpose

    ! Parameter Constants

    integer, parameter :: &
      d_uniform_extra = 2   ! Number of variables that are included in the uniform sample but not in
                            ! the lognormal sample. Currently:
                            !
                            ! d_variables+1: Mixture component, for choosing PDF component
                            ! d_variables+2: Precipitation fraction, for determining precipitation

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nz                 ! Number of vertical model levels

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      delta_zm, &  ! Difference in moment. altitudes    [m]
      rcm          ! Liquid water mixing ratio          [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Lscale       ! Turbulent mixing length            [m]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    
    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(num_samples) :: &
      lh_sample_point_weights

    ! More Input Variables!
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(d_variables,nz), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    ! Local variables
    real( kind = core_rknd ), dimension(nz,num_samples,(d_variables+d_uniform_extra)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    integer :: &
      k_lh_start, & ! Height for preferentially sampling within cloud
      k, sample, i  ! Loop iterators

    logical, dimension(nz,num_samples) :: &
      l_in_precip   ! Whether sample is in precipitation

    logical :: l_error, l_error_in_sub

    ! Precipitation fraction in a component of the PDF, for each sample
    real( kind = core_rknd ), dimension(num_samples) :: precip_frac_i

    ! ---- Begin Code ----

    l_error = .false.

    ! Sanity checks for l_lh_importance_sampling
    if ( l_lh_importance_sampling .and. sequence_length /= 1 ) then
      write(fstderr,*) "Cloud weighted sampling requires sequence length be equal to 1."
      stop "Fatal error."
    end if

    !--------------------------------------------------------------
    ! Latin hypercube sampling
    !--------------------------------------------------------------

    ! Compute k_lh_start, the starting level for SILHS sampling
    k_lh_start = compute_k_lh_start( nz, rcm, pdf_params )

    ! Generate a uniform sample at k_lh_start
    call generate_uniform_k_lh_start &
         ( iter, d_variables, d_uniform_extra, num_samples, sequence_length, & ! Intent(in)
           pdf_params(k_lh_start), hydromet_pdf_params(k_lh_start), &          ! Intent(in)
           X_u_all_levs(k_lh_start,:,:), lh_sample_point_weights )             ! Intent(out)

    ! Generate uniform sample at other height levels by vertically correlating them
    call vertical_overlap_driver &
         ( nz, d_variables, d_uniform_extra, num_samples, &     ! Intent(in)
           k_lh_start, delta_zm, rcm, Lscale, rho_ds_zt, &      ! Intent(in)
           X_u_all_levs )                                       ! Intent(inout)

    do k = 1, nz
      ! Determine mixture component for all levels
      where ( in_mixt_comp_1( X_u_all_levs(k,:,d_variables+1), pdf_params(k)%mixt_frac ) )
        X_mixt_comp_all_levs(k,:) = 1
      else where
        X_mixt_comp_all_levs(k,:) = 2
      end where

      ! Determine precipitation fraction
      where ( X_mixt_comp_all_levs(k,:) == 1 )
        precip_frac_i(:) = hydromet_pdf_params(k)%precip_frac_1
      else where
        precip_frac_i(:) = hydromet_pdf_params(k)%precip_frac_2
      end where

      ! Determine precipitation for all levels
      where ( in_precipitation( X_u_all_levs(k,:,d_variables+2), &
                  precip_frac_i(:) ) )
        l_in_precip(k,:) = .true.
      else where
        l_in_precip(k,:) = .false.
      end where

    end do ! k = 1 .. nz

    call stats_accumulate_uniform_lh( nz, num_samples, l_in_precip, X_mixt_comp_all_levs, &
                                      X_u_all_levs(:,:,iiPDF_chi), pdf_params, &
                                      lh_sample_point_weights, k_lh_start )

    ! Check to ensure uniform variates are in the appropriate range
    do sample=1, num_samples
      do k=1, nz
        do i=1, d_variables+d_uniform_extra
          if ( X_u_all_levs(k,sample,i) >= one ) then
            X_u_all_levs(k,sample,i) = one - epsilon( X_u_all_levs(k,sample,i) )
          else if ( X_u_all_levs(k,sample,i) <= zero ) then
            X_u_all_levs(k,sample,i) = epsilon( X_u_all_levs(k,sample,i) )
          end if
        end do
      end do
    end do

    ! Sample loop
    do k = 1, nz
      ! Generate LH sample, represented by X_u and X_nl, for level k
      do sample = 1, num_samples, 1
        call transform_uniform_sample_to_pdf &
             ( d_variables, d_uniform_extra, & ! In
               mu1(:,k), mu2(:,k), sigma1(:,k), sigma2(:,k), & ! In
               corr_cholesky_mtx_1(:,:,k), & ! In
               corr_cholesky_mtx_2(:,:,k), & ! In
               X_u_all_levs(k,sample,:), X_mixt_comp_all_levs(k,sample), & ! In
               l_in_precip(k,sample), & ! In
               X_nl_all_levs(k,sample,:) ) ! Out
      end do ! sample = 1, num_samples, 1
    end do ! k = 1, nz

    if ( l_output_2D_lognormal_dist ) then
      ! Eric Raut removed lh_rt and lh_thl from call to output_2D_lognormal_dist_file
      ! because they are no longer generated in lh_subcolumn_generator.
      call output_2D_lognormal_dist_file( nz, num_samples, d_variables, &
                                          real(X_nl_all_levs, kind = stat_rknd) )
    end if
    if ( l_output_2D_uniform_dist ) then
      call output_2D_uniform_dist_file( nz, num_samples, d_variables+2, &
                                        X_u_all_levs, &
                                        X_mixt_comp_all_levs, &
                                        lh_sample_point_weights )
    end if

    ! Various nefarious assertion checks
    if ( clubb_at_least_debug_level( 2 ) ) then

      ! Simple assertion check to ensure uniform variates are in the appropriate
      ! range
      if ( any( X_u_all_levs <= zero .or. X_u_all_levs >= one ) ) then
        write(fstderr,*) "A uniform variate was not in the correct range."
        l_error = .true.
      end if

      do k=2, nz

        call assert_consistent_cloud_frac( pdf_params(k), l_error_in_sub )
        l_error = l_error .or. l_error_in_sub

        ! Check for correct transformation in normal space
        call assert_correct_cloud_normal( num_samples, X_u_all_levs(k,:,iiPDF_chi), & ! In
                                          X_nl_all_levs(k,:,iiPDF_chi), & ! In
                                          X_mixt_comp_all_levs(k,:), & ! In
                                          pdf_params(k)%cloud_frac_1, & ! In
                                          pdf_params(k)%cloud_frac_2, & ! In
                                          l_error_in_sub ) ! Out
        l_error = l_error .or. l_error_in_sub

      end do ! k=2, nz

    end if ! clubb_at_least_debug_level( 2 )

    ! Stop the run if an error occurred
    if ( l_error ) then
      stop "Fatal error in lh_subcolumn_generator"
    end if

    return
  end subroutine lh_subcolumn_generator
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine generate_uniform_k_lh_start &
             ( iter, d_variables, d_uniform_extra, num_samples, sequence_length, &
               pdf_params, hydromet_pdf_params, &
               X_u_k_lh_start, lh_sample_point_weights )

  ! Description:
  !   Generates a uniform sample, X_u, at a single height level, applying Latin
  !   hypercube and importance sampling where configured.

  ! References:
  !   V. E. Larson and D. P. Schanen, 2013. The Subgrid Importance Latin
  !   Hypercube Sampler (SILHS): a multivariate subcolumn generator
  !----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd                   ! Precision

    use constants_clubb, only: &
      one, fstderr                ! Constant(s)

    use parameters_silhs, only: &
      l_lh_straight_mc, &         ! Variable(s)
      l_lh_importance_sampling

    use pdf_parameter_module, only: &
      pdf_parameter               ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter      ! Type

    use generate_uniform_sample_module, only: &
      rand_uniform_real, &        ! Procedure(s)
      generate_uniform_lh_sample

    use silhs_importance_sample_module, only: &
      importance_sampling_driver, & ! Procedure(s)
      cloud_weighted_sampling_driver

    use latin_hypercube_arrays, only: &
      one_height_time_matrix      ! Variable

    use corr_varnce_module, only: &
      iiPDF_chi                   ! Variable

    implicit none

    ! Local Constants
    logical, parameter :: &
      l_lh_old_cloud_weighted  = .false. ! Use the old method of importance sampling that
                                         ! places one point in cloud and one point out of
                                         ! cloud

    ! Input Variables
    integer, intent(in) :: &
      iter,              &        ! Model iteration number
      d_variables,       &        ! Number of variates in CLUBB's PDF
      d_uniform_extra,   &        ! Uniform variates included in uniform sample but not
                                  !  in normal/lognormal sample
      num_samples,       &        ! Number of SILHS sample points
      sequence_length             ! Number of timesteps before new sample points are picked

    type(pdf_parameter), intent(in) :: &
      pdf_params                  ! The PDF parameters at k_lh_start

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params

    ! Output Variables
    real( kind = core_rknd ), dimension(num_samples,d_variables+d_uniform_extra), intent(out) :: &
      X_u_k_lh_start              ! Uniform sample at k_lh_start

    real( kind = core_rknd ), dimension(num_samples), intent(out) :: &
      lh_sample_point_weights     ! Weight of each sample point (all equal to one if importance
                                  ! sampling is not used)

    ! Local Variables
    integer :: &
      i, sample

  !----------------------------------------------------------------------
    !----- Begin Code -----

    ! Sanity check
    if ( l_lh_old_cloud_weighted .and. mod( num_samples, 2 ) /= 0 ) then
      write(fstderr,*) "Old cloud weighted sampling requires num_samples to be divisible by 2."
      stop "Fatal error."
    end if

    if ( l_lh_straight_mc ) then

      ! Do a straight Monte Carlo sample without LH or importance sampling.
      do i=1, d_variables+d_uniform_extra
        do sample=1, num_samples
          X_u_k_lh_start(sample,i) = rand_uniform_real( )
        end do
      end do

      ! Importance sampling is not performed, so all sample points have the same weight!!
      lh_sample_point_weights(1:num_samples)  =  one

    else ! .not. l_lh_straight_mc

      ! Generate a uniformly distributed Latin hypercube sample
      call generate_uniform_lh_sample( iter, num_samples, sequence_length, &  ! Intent(in)
                                       d_variables+d_uniform_extra, &         ! Intent(in)
                                       X_u_k_lh_start(:,:) )                  ! Intent(out)

      if ( l_lh_importance_sampling ) then

        if ( l_lh_old_cloud_weighted ) then

          call cloud_weighted_sampling_driver &
               ( num_samples, one_height_time_matrix(:,iiPDF_chi), & ! In
                 one_height_time_matrix(:,d_variables+1), & ! In
                 pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, & ! In
                 pdf_params%mixt_frac, & ! In
                 X_u_k_lh_start(:,iiPDF_chi), & ! In/Out
                 X_u_k_lh_start(:,d_variables+1), & ! In/Out
                 lh_sample_point_weights ) ! Out

        else ! .not. l_lh_old_cloud_weighted

          call importance_sampling_driver &
               ( num_samples, pdf_params, hydromet_pdf_params, & ! In
                 X_u_k_lh_start(:,iiPDF_chi), & ! In/Out
                 X_u_k_lh_start(:,d_variables+1), & ! In/Out
                 X_u_k_lh_start(:,d_variables+2), & ! In/Out
                 lh_sample_point_weights ) ! Out

        end if ! l_lh_old_cloud_weighted

      else

        ! No importance sampling is performed, so all sample points have the same weight.
        lh_sample_point_weights(1:num_samples) = one

      end if ! l_lh_importance_sampling

    end if ! l_lh_straight_mc

    return
  end subroutine generate_uniform_k_lh_start
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine vertical_overlap_driver &
             ( nz, d_variables, d_uniform_extra, num_samples, &
               k_lh_start, delta_zm, rcm, Lscale, rho_ds_zt, &
               X_u_all_levs )

  ! Description:
  !   Takes a uniform sample at k_lh_start and correlates it vertically
  !   to other height levels

  ! References:
  !   none
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd      ! Precision

    use constants_clubb, only: &
      fstderr, &     ! Constant(s)
      one,     &
      zero

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use grid_class, only: &
      gr             ! Variable

    use corr_varnce_module, only: &
      iiPDF_chi      ! Variable

    use parameters_silhs, only: &
      l_Lscale_vert_avg  ! Variable

    use fill_holes, only: &
      vertical_avg  ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,              &
      d_variables,     &
      d_uniform_extra, &
      num_samples,     &
      k_lh_start

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      delta_zm,        &  ! Difference in altitude between momentum levels    [m]
      rcm,             &  ! Liquid water mixing ratio                         [kg/kg]
      Lscale,          &  ! Turbulent mixing length                           [m]
      rho_ds_zt           ! Dry, static density on thermodynamic levels       [kg/m^3]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz,num_samples,d_variables+d_uniform_extra), &
        intent(inout) :: &
      X_u_all_levs        ! A full uniform sample

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      Lscale_vert_avg, &  ! 3pt vertical average of Lscale                    [m]
      X_vert_corr         ! Vertical correlations between height levels       [-]

    integer :: k, km1, kp1, sample, ivar

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    if ( l_Lscale_vert_avg ) then
      ! Determine 3pt vertically averaged Lscale
      do k = 1, nz, 1
        kp1 = min( k+1, nz )
        km1 = max( k-1, 1 )
        Lscale_vert_avg(k) = vertical_avg &
                             ( (kp1-km1+1), rho_ds_zt(km1:kp1), &
                               Lscale(km1:kp1), gr%invrs_dzt(km1:kp1) )
      end do
    else
        Lscale_vert_avg = Lscale 
    end if

    X_vert_corr(1:nz) = compute_vert_corr( nz, delta_zm, Lscale_vert_avg, rcm )

    ! Assertion check for the vertical correlation
    if ( clubb_at_least_debug_level( 2 ) ) then
      if ( any( X_vert_corr > one ) .or. any( X_vert_corr < zero ) ) then
        write(fstderr,*) "The vertical correlation in latin_hypercube_driver"// &
          "is not in the correct range"
        do k = 1, nz
          write(fstderr,*) "k = ", k,  "Vert. correlation = ", X_vert_corr(k)
        end do
        stop "Fatal error in vertical_overlap_driver"
      end if ! Some correlation isn't between [0,1]
    end if ! clubb_at_least_debug_level 1

    do sample = 1, num_samples
      ! Correlate chi vertically
      call compute_arb_overlap &
           ( nz, k_lh_start, &  ! In
             X_vert_corr, & ! In
             X_u_all_levs(:,sample,iiPDF_chi) ) ! Inout
      ! Correlate the d+1 variate vertically (used to compute the mixture
      ! component later)
      call compute_arb_overlap &
           ( nz, k_lh_start, &  ! In
             X_vert_corr, & ! In
             X_u_all_levs(:,sample,d_variables+1) ) ! Inout

      ! Correlate the d+2 variate vertically (used to determine precipitation
      ! later)
      call compute_arb_overlap &
           ( nz, k_lh_start, &  ! In
             X_vert_corr, & ! In
             X_u_all_levs(:,sample,d_variables+2) ) ! Inout

      ! Use these lines to make all variates vertically correlated, using the
      ! same correlation we used above for chi and the d+1 variate
      do ivar = 1, d_variables
        if ( ivar /= iiPDF_chi ) then
          call compute_arb_overlap &
               ( nz, k_lh_start, &  ! In
                 X_vert_corr, & ! In
                 X_u_all_levs(:,sample,ivar) ) ! Inout
        end if
      end do ! 1..d_variables
    end do ! 1..num_samples

    return
  end subroutine vertical_overlap_driver
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_k_lh_start( nz, rcm, pdf_params ) result( k_lh_start )

  ! Description:
  !   Determines the starting SILHS sample level

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd      ! Constant

    use constants_clubb, only: &
      cloud_frac_min ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use pdf_utilities, only: &
      compute_mean_binormal  ! Procedure

    use math_utilities, only: &
      rand_integer_in_range  ! Procedure

    use parameters_silhs, only: &
      l_rcm_in_cloud_k_lh_start,  &  ! Variable(s)
      l_random_k_lh_start

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm         ! Liquid water mixing ratio               [kg/kg]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params  ! PDF parameters       [units vary]

    ! Output Variable
    integer :: &
      k_lh_start  ! Starting SILHS sample level

    ! Local Variables
    integer :: &
      k_lh_start_rcm_in_cloud, &
      k_lh_start_rcm

    real( kind = core_rknd ), dimension(nz) :: &
      rcm_pdf, cloud_frac_pdf

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_rcm_in_cloud_k_lh_start .or. l_random_k_lh_start ) then
      rcm_pdf = compute_mean_binormal( pdf_params%rc_1, pdf_params%rc_2, pdf_params%mixt_frac )
      cloud_frac_pdf = compute_mean_binormal( pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, &
                                              pdf_params%mixt_frac )
      k_lh_start_rcm_in_cloud = maxloc( rcm_pdf / max( cloud_frac_pdf, cloud_frac_min ), 1 )
    end if

    if ( .not. l_rcm_in_cloud_k_lh_start .or. l_random_k_lh_start ) then
      k_lh_start_rcm    = maxloc( rcm, 1 )
    end if

    if ( l_random_k_lh_start ) then
      if ( k_lh_start_rcm_in_cloud == k_lh_start_rcm ) then
        k_lh_start = k_lh_start_rcm
      else
        ! Pick a random height level between k_lh_start_rcm and
        ! k_lh_start_rcm_in_cloud
        if ( k_lh_start_rcm_in_cloud > k_lh_start_rcm ) then
          k_lh_start = rand_integer_in_range( k_lh_start_rcm, k_lh_start_rcm_in_cloud )
        else if ( k_lh_start_rcm > k_lh_start_rcm_in_cloud ) then
          k_lh_start = rand_integer_in_range( k_lh_start_rcm_in_cloud, k_lh_start_rcm )
        end if
      end if
    else if ( l_rcm_in_cloud_k_lh_start ) then
      k_lh_start = k_lh_start_rcm_in_cloud
    else ! .not. l_random_k_lh_start .and. .not. l_rcm_in_cloud_k_lh_start
      k_lh_start = k_lh_start_rcm
    end if

    ! If there's no cloud k_lh_start appears to end up being 1. Check if
    ! k_lh_start is 1 or nz and set it to the middle of the domain in that
    ! case.
    if ( k_lh_start == nz .or. k_lh_start == 1 ) then
      k_lh_start = nz / 2
    end if

    return
  end function compute_k_lh_start
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine clip_transform_silhs_output( nz, num_samples, d_variables, &         ! In
                                          X_mixt_comp_all_levs, X_nl_all_levs, &  ! In
                                          pdf_params, l_use_Ncn_to_Nc, &          ! In
                                          lh_clipped_vars )                       ! Out

  ! Description:
  !   Derives from the SILHS sampling structure X_nl_all_levs the variables
  !   rt, thl, rc, rv, and Nc, for all sample points and height levels.

  ! References:
  !   ticket:751
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      zero, &   ! Constant(s)
      rt_tol

    use pdf_parameter_module, only: &
      pdf_parameter ! Type

    use corr_varnce_module, only: &
      iiPDF_chi, &  ! Variable(s)
      iiPDF_eta, &
      iiPDF_Ncn

    use transform_to_pdf_module, only: &
      chi_eta_2_rtthl ! Awesome procedure

    implicit none

    ! Input Variables
    logical, intent(in) :: &
      l_use_Ncn_to_Nc        ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                             ! Ncn_to_Nc might cause problems with the MG microphysics 
                             ! since the changes made here (Nc-tendency) are not fed into 
                             ! the microphysics

    integer, intent(in) :: &
      nz,          &         ! Number of vertical levels
      num_samples, &         ! Number of SILHS sample points
      d_variables            ! Number of variates in X_nl

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs   ! Which component this sample is in (1 or 2)

    real( kind = core_rknd ), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs          ! A SILHS sample

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params             ! **The** PDF parameters!

    ! Output Variables
    type(lh_clipped_variables_type), dimension(nz,num_samples), intent(out) :: &
      lh_clipped_vars        ! SILHS clipped and transformed variables

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      lh_rt,   &             ! Total water mixing ratio            [kg/kg]
      lh_thl,  &             ! Liquid potential temperature        [K]
      lh_rc,   &             ! Cloud water mixing ratio            [kg/kg]
      lh_rv,   &             ! Vapor water mixing ratio            [kg/kg]
      lh_Nc                  ! Cloud droplet number concentration  [#/kg]

    real( kind = core_rknd ) :: &
      rt_1,    &
      rt_2,    &
      thl_1,   &
      thl_2,   &
      crt_1,   &
      crt_2,   &
      cthl_1,  &
      cthl_2,  &
      chi_1,   &
      chi_2,   &
      lh_chi,  &
      lh_eta

    integer :: &
      isample, k

    logical :: &
      l_rt_clipped, &
      l_rv_clipped

  !-----------------------------------------------------------------------

    l_rt_clipped = .false.
    l_rv_clipped = .false.

    !----- Begin Code -----

    ! Loop over all thermodynamic levels above the model lower boundary.
    do k=2, nz

      ! Enter the PDF parameters!!
      rt_1   = pdf_params(k)%rt_1
      rt_2   = pdf_params(k)%rt_2
      thl_1  = pdf_params(k)%thl_1
      thl_2  = pdf_params(k)%thl_2
      crt_1  = pdf_params(k)%crt_1
      crt_2  = pdf_params(k)%crt_2
      cthl_1 = pdf_params(k)%cthl_1
      cthl_2 = pdf_params(k)%cthl_2
      chi_1  = pdf_params(k)%chi_1
      chi_2  = pdf_params(k)%chi_2

      do isample=1, num_samples

        lh_chi = X_nl_all_levs(k,isample,iiPDF_chi)
        lh_eta = X_nl_all_levs(k,isample,iiPDF_eta)

        ! Compute lh_rt and lh_thl
        call chi_eta_2_rtthl( rt_1, thl_1, rt_2, thl_2, &                        ! Intent(in)
                              crt_1, cthl_1, crt_2, cthl_2, &                    ! Intent(in)
                              chi_1, chi_2, &                                    ! Intent(in)
                              lh_chi, lh_eta, X_mixt_comp_all_levs(k,isample), & ! Intent(in)
                              lh_rt(k,isample), lh_thl(k,isample) )              ! Intent(out)

        ! If necessary, clip rt
        if ( lh_rt(k,isample) < rt_tol ) then
          lh_rt(k,isample) = rt_tol
          l_rt_clipped = .true.
        end if

        ! Compute lh_rc
        lh_rc(k,isample) = chi_to_rc( X_nl_all_levs(k,isample,iiPDF_chi) )
        ! Clip lh_rc.
        if ( lh_rc(k,isample) > lh_rt(k,isample) - rt_tol ) then
          lh_rc(k,isample) = lh_rt(k,isample) - rt_tol
        end if

        ! Compute lh_rv
        lh_rv(k,isample) = lh_rt(k,isample) - lh_rc(k,isample)

        ! Compute lh_Nc
        if ( l_use_Ncn_to_Nc ) then
           lh_Nc(k,isample) = Ncn_to_Nc( X_nl_all_levs(k,isample,iiPDF_Ncn), &
                                         X_nl_all_levs(k,isample,iiPDF_chi) )
        else
           lh_Nc(k,isample) = X_nl_all_levs(k,isample,iiPDF_Ncn)
        endif ! l_use_Ncn_to_Nc

      end do ! isample=1, num_samples

    end do ! k=2, nz

    ! These parameters are not computed at the model lower level.
    lh_rt(1,:) = zero
    lh_thl(1,:) = zero
    lh_rc(1,:) = zero
    lh_rv(1,:) = zero
    lh_Nc(1,:) = zero

    ! Pack into output structure
    lh_clipped_vars(:,:)%rt  = lh_rt(:,:)
    lh_clipped_vars(:,:)%thl = lh_thl(:,:)
    lh_clipped_vars(:,:)%rc  = lh_rc(:,:)
    lh_clipped_vars(:,:)%rv  = lh_rv(:,:)
    lh_clipped_vars(:,:)%Nc  = lh_Nc(:,:)

    return
  end subroutine clip_transform_silhs_output
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine assert_consistent_cloud_frac( pdf_params, l_error )

  ! Description:
  !   Performs an assertion check that cloud_frac_i is consistent with chi_i and
  !   stdev_chi_i in pdf_params for each PDF component.

  ! References:
  !   Eric Raut
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr          ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter    ! Type

    implicit none

    ! Input Variables
    type(pdf_parameter), intent(in) :: &
      pdf_params       ! PDF parameters, containing distribution of chi     [units vary]
                       ! and cloud fraction

    ! Output Variables
    logical, intent(out) :: &
      l_error          ! True if the assertion check fails

    ! Local Variables
    logical :: &
      l_error_in_sub

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    ! Perform assertion check for PDF component 1
    call assert_consistent_cf_component &
         ( pdf_params%chi_1, pdf_params%stdev_chi_1, pdf_params%cloud_frac_1, & ! Intent(in)
           l_error_in_sub )                                                    ! Intent(out)

    l_error = l_error .or. l_error_in_sub
    if ( l_error_in_sub ) then
      write(fstderr,*) "Cloud fraction is inconsistent in PDF component 1"
    end if

    ! Perform assertion check for PDF component 2
    call assert_consistent_cf_component &
         ( pdf_params%chi_2, pdf_params%stdev_chi_2, pdf_params%cloud_frac_2, & ! Intent(in)
           l_error_in_sub )                                                    ! Intent(out)

    l_error = l_error .or. l_error_in_sub
    if ( l_error_in_sub ) then
      write(fstderr,*) "Cloud fraction is inconsistent in PDF component 2"
    end if

    return
  end subroutine assert_consistent_cloud_frac
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  subroutine assert_consistent_cf_component( mu_chi_i, sigma_chi_i, cloud_frac_i, &
                                             l_error )

  ! Description:
  !   Performs an assertion check that cloud_frac_i is consistent with chi_i and
  !   stdev_chi_i for a PDF component.
  !
  !   The SILHS sample generation process relies on precisely a cloud_frac
  !   amount of mass in the cloudy portion of the PDF of chi, that is, where
  !   chi > 0. In other words, the probability that chi > 0 should be exactly
  !   cloud_frac.
  !
  !   Stated even more mathematically, CDF_chi(0) = 1 - cloud_frac, where
  !   CDF_chi is the cumulative distribution function of chi. This can be
  !   expressed as invCDF_chi(1 - cloud_frac) = zero.
  !
  !   This subroutine uses ltqnorm, which is apparently a fancy name for the
  !   inverse cumulative distribution function of the standard normal
  !   distribution.

  ! References:
  !   Eric Raut
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd        ! Constant(s)

    use constants_clubb, only: &
      fstderr, &       ! Constant(s)
      zero, &
      one

    use pdf_parameter_module, only: &
      pdf_parameter    ! Type

    use transform_to_pdf_module, only: &
      ltqnorm          ! Procedure

    implicit none

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      ! Values below ltqnorm_min_arg (or above 1-ltqnorm_min_arg) will not be
      ! supplied as arguments to the ltqnorm function.
      ltqnorm_min_arg = 1.0e-5_core_rknd, &
      ! It will be verified that all values of the chi uniform variate will not
      ! trigger the in/out of cloud assertion error, except for those values
      ! within a box with the following half-width, centered around the cloud
      ! fraction in this component.
      box_half_width = 5.0e-6_core_rknd

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_i,      & ! Mean of chi in a PDF component
      sigma_chi_i,   & ! Standard deviation of chi in a PDF component
      cloud_frac_i     ! Cloud fraction in a PDF component

    ! Output Variables
    logical, intent(out) :: &
      l_error          ! True if the assertion check fails

    ! Local Variables
    real( kind = core_rknd ) :: chi, chi_std_normal, one_minus_ltqnorm_arg

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    ! Check left end of box.
    one_minus_ltqnorm_arg = cloud_frac_i - box_half_width
    ! Do not bother to check this end of the box if it dips below the minimum
    ! ltqnorm argument.
    if ( one_minus_ltqnorm_arg >= ltqnorm_min_arg ) then

      if ( one_minus_ltqnorm_arg > (one - ltqnorm_min_arg) ) then
        one_minus_ltqnorm_arg = one - ltqnorm_min_arg
      end if

      chi_std_normal = ltqnorm( one - one_minus_ltqnorm_arg )
      chi = chi_std_normal * sigma_chi_i + mu_chi_i
      if ( chi <= zero ) then
        l_error = .true.
        write(fstderr,*) 'chi (left side of box) = ', chi
      end if
    end if ! one_minus_ltqnorm_arg >= ltqnorm_min_arg

    ! Check right end of box.
    one_minus_ltqnorm_arg = cloud_frac_i + box_half_width
    ! Do not bother to check this end of the box if it exceeds the maximum
    ! ltqnorm argument
    if ( one_minus_ltqnorm_arg <= one-ltqnorm_min_arg ) then

      if ( one_minus_ltqnorm_arg < ltqnorm_min_arg ) then
        one_minus_ltqnorm_arg = ltqnorm_min_arg
      end if

      chi_std_normal = ltqnorm( one - one_minus_ltqnorm_arg )
      chi = chi_std_normal * sigma_chi_i + mu_chi_i
      if ( chi > zero ) then
        l_error = .true.
        write(fstderr,*) 'chi (right side of box) = ', chi
      end if
    end if ! one_minus_ltqnorm_arg >= ltqnorm_min_arg

    if ( l_error ) then
      write(fstderr,*) "In assert_consistent_cf_component, cloud_frac_i is inconsistent with &
                       &mu_chi_i and stdev_chi_i."
      write(fstderr,*) "mu_chi_i = ", mu_chi_i
      write(fstderr,*) "sigma_chi_i = ", sigma_chi_i
      write(fstderr,*) "cloud_frac_i = ", cloud_frac_i
    end if

    return
  end subroutine assert_consistent_cf_component
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine assert_correct_cloud_normal( num_samples, X_u_chi, X_nl_chi, X_mixt_comp, &
                                          cloud_frac_1, cloud_frac_2, &
                                          l_error )

  ! Description:
  !   Asserts that all SILHS sample points that are in cloud in uniform space
  !   are in cloud in normal space, and that all SILHS sample points that are
  !   in clear air in uniform space are in clear air in normal space.
  
  ! References:
  !   None
  !-----------------------------------------------------------------------
  
    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one,     &
      zero, &      ! Constant(s)
      fstderr

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples            ! Number of SILHS sample points

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      X_u_chi,  &            ! Samples of chi in uniform space
      X_nl_chi               ! Samples of chi in normal space

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp            ! PDF component of each sample

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac_1,   &      ! Cloud fraction in PDF component 1
      cloud_frac_2           ! Cloud fraction in PDF component 2

    ! Output Variables
    logical, intent(out) :: &
      l_error                ! True if the assertion check fails

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      error_threshold = 1000._core_rknd * epsilon( X_nl_chi ) ! A threshold to determine whether a
                                                              ! rogue value triggers the assertion
                                                              ! check.

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac_i

    integer :: sample

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    do sample = 1, num_samples, 1

      ! Determine the appropriate cloud fraction
      if ( X_mixt_comp(sample) == 1 ) then
        cloud_frac_i = cloud_frac_1
      else if ( X_mixt_comp(sample) == 2 ) then
        cloud_frac_i = cloud_frac_2
      end if

      if ( X_u_chi(sample) < (one - cloud_frac_i) ) then

        ! The uniform sample is in clear air
        if ( X_nl_chi(sample) > error_threshold ) then
          l_error = .true.
        end if

      else if ( X_u_chi(sample) >= (one - cloud_frac_i) .and. &
                X_u_chi(sample) < one ) then

        ! The uniform sample is in cloud
        if ( X_nl_chi(sample) <= -error_threshold ) then
          l_error = .true.
        end if

      else
        stop "X_u_chi not in correct range in assert_correct_cloud_normal"
      end if

    end do ! 1..num_samples

    if ( l_error ) then
      write(fstderr,*) "In assert_correct_cloud_normal:"
      write(fstderr,*) "The 'cloudiness' of points in uniform and normal space is not consistent"
      write(fstderr,'(4X,A,A)')  "X_u_chi         ", "X_nl_chi "
      do sample = 1, num_samples, 1
        write(fstderr,'(I4,2G20.4)') &
          sample, X_u_chi(sample), X_nl_chi(sample)
      end do
      ! This will hopefully stop the run at some unknown point in the future
      l_error = .true.
    end if  ! in_cloud_points /= out_of_cloud_points

    return
  end subroutine assert_correct_cloud_normal
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_output &
             ( fname_prefix, fdir, stats_tout, nz, &
               stats_zt, time_initial, num_samples )
!-------------------------------------------------------------------------------

    use corr_varnce_module, only: &
      iiPDF_chi, & ! Variables
      iiPDF_eta, &
      iiPDF_w, &
      iiPDF_rr, & 
      iiPDF_ri, &
      iiPDF_rs, &
      iiPDF_rg, &
      iiPDF_Nr, &
      iiPDF_Ni, &
      iiPDF_Ns, &
      iiPDF_Ng, &
      iiPDF_Ncn

    use clubb_precision, only: &
      time_precision, & ! Constant
      core_rknd

    use output_2D_samples_module, only: &
      open_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Instance of a type
      uniform_sample_file

    use corr_varnce_module, only: &
      d_variables ! Variable


    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      fname_prefix, & ! Prefix for file name
      fdir            ! Directory for output

    real(kind=core_rknd), intent(in) :: &
      stats_tout    ! Frequency to write to disk        [s]

    real(kind=time_precision), intent(in) :: &
      time_initial  ! Initial time                      [s]

    integer, intent(in) :: &
      nz ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      stats_zt ! Altitudes [m]

    integer, intent(in) :: num_samples

    ! Local Variables
    character(len=100), allocatable, dimension(:) :: &
      variable_names, variable_descriptions, variable_units

    integer :: i

    ! ---- Begin Code ----

    if ( l_output_2D_lognormal_dist ) then

      allocate( variable_names(d_variables), variable_descriptions(d_variables), &
                variable_units(d_variables) )

      variable_names(iiPDF_chi)        = "chi"
      variable_descriptions(iiPDF_chi) = "The variable 's' from Mellor 1977"
      variable_units(iiPDF_chi)        = "kg/kg"

      variable_names(iiPDF_eta)        = "eta"
      variable_descriptions(iiPDF_eta) = "The variable 't' from Mellor 1977"
      variable_units(iiPDF_eta)        = "kg/kg"

      variable_names(iiPDF_w)        = "w"
      variable_descriptions(iiPDF_w) = "Vertical velocity"
      variable_units(iiPDF_w)        = "m/s"

      if ( iiPDF_rr > 0 ) then
        variable_names(iiPDF_rr)        = "rr"
        variable_descriptions(iiPDF_rr) = "Rain water mixing ratio"
        variable_units(iiPDF_rr)        = "kg/kg"
      end if
      if ( iiPDF_ri > 0 ) then
        variable_names(iiPDF_ri)        = "ri"
        variable_descriptions(iiPDF_ri) = "Ice water mixing ratio"
        variable_units(iiPDF_ri)        = "kg/kg"
      end if
      if ( iiPDF_rs > 0 ) then
        variable_names(iiPDF_rs)        = "rs"
        variable_descriptions(iiPDF_rs) = "Snow water mixing ratio"
        variable_units(iiPDF_rs)        = "kg/kg"
      end if
      if ( iiPDF_rg > 0 ) then
        variable_names(iiPDF_rg)        = "rg"
        variable_descriptions(iiPDF_rg) = "Graupel water mixing ratio"
        variable_units(iiPDF_rg)        = "kg/kg"
      end if

      if ( iiPDF_Nr > 0 ) then
        variable_names(iiPDF_Nr)        = "Nr"
        variable_descriptions(iiPDF_Nr) = "Rain droplet number concentration"
        variable_units(iiPDF_Nr)        = "count/kg"
      end if
      if ( iiPDF_Ncn > 0 ) then
        variable_names(iiPDF_Ncn)        = "Ncn"
        variable_descriptions(iiPDF_Ncn) = "Cloud nuclei concentration (simplified)"
        variable_units(iiPDF_Ncn)        = "count/kg"
      end if
      if ( iiPDF_Ni > 0 ) then
        variable_names(iiPDF_Ni)        = "Ni"
        variable_descriptions(iiPDF_Ni) = "Ice number concentration"
        variable_units(iiPDF_Ni)        = "count/kg"
      end if
      if ( iiPDF_Ns > 0 ) then
        variable_names(iiPDF_Ns)        = "Ns"
        variable_descriptions(iiPDF_Ns) = "Snow number concentration"
        variable_units(iiPDF_Ns)        = "count/kg"
      end if
      if ( iiPDF_Ng > 0 ) then
        variable_names(iiPDF_Ng)        = "Ng"
        variable_descriptions(iiPDF_Ng) = "Graupel number concentration"
        variable_units(iiPDF_Ng)        = "count/kg"
      end if

      call open_2D_samples_file( nz, num_samples, d_variables, & ! In
                                 trim( fname_prefix )//"_nl", fdir, & ! In
                                 time_initial, stats_tout, stats_zt, variable_names, & ! In
                                 variable_descriptions, variable_units, & ! In
                                 lognormal_sample_file ) ! In/Out

      deallocate( variable_names, variable_descriptions, variable_units )

    end if

    if ( l_output_2D_uniform_dist ) then

      allocate( variable_names(d_variables+4), variable_descriptions(d_variables+4), &
                variable_units(d_variables+4) )

      ! The uniform distribution corresponds to all the same variables as X_nl,
      ! except the d+1 component is the mixture component.

      variable_names(iiPDF_chi)        = "chi"
      variable_descriptions(iiPDF_chi) = "Uniform dist of the variable 's' from Mellor 1977"

      variable_names(iiPDF_eta)        = "eta"
      variable_descriptions(iiPDF_eta) = "Uniform dist of the variable 't' from Mellor 1977"

      variable_names(iiPDF_w)        = "w"
      variable_descriptions(iiPDF_w) = "Uniform dist of the vertical velocity"


      if ( iiPDF_rr > 0 ) then
        variable_names(iiPDF_rr)        = "rr"
        variable_descriptions(iiPDF_rr) = "Rain water mixing ratio"
        variable_units(iiPDF_rr)        = "kg/kg"
      end if
      if ( iiPDF_ri > 0 ) then
        variable_names(iiPDF_ri)        = "ri"
        variable_descriptions(iiPDF_ri) = "Ice water mixing ratio"
        variable_units(iiPDF_ri)        = "kg/kg"
      end if
      if ( iiPDF_rs > 0 ) then
        variable_names(iiPDF_rs)        = "rs"
        variable_descriptions(iiPDF_rs) = "Snow water mixing ratio"
        variable_units(iiPDF_rs)        = "kg/kg"
      end if
      if ( iiPDF_rg > 0 ) then
        variable_names(iiPDF_rg)        = "rg"
        variable_descriptions(iiPDF_rg) = "Graupel water mixing ratio"
        variable_units(iiPDF_rg)        = "kg/kg"
      end if

      if ( iiPDF_Nr > 0 ) then
        variable_names(iiPDF_Nr)        = "Nr"
        variable_descriptions(iiPDF_Nr) = "Rain droplet number concentration"
        variable_units(iiPDF_Nr)        = "count/kg"
      end if
      if ( iiPDF_Ncn > 0 ) then
        variable_names(iiPDF_Ncn)        = "Ncn"
        variable_descriptions(iiPDF_Ncn) = "Cloud nuclei concentration (simplified)"
        variable_units(iiPDF_Ncn)        = "count/kg"
      end if
      if ( iiPDF_Ni > 0 ) then
        variable_names(iiPDF_Ni)        = "Ni"
        variable_descriptions(iiPDF_Ni) = "Ice number concentration"
        variable_units(iiPDF_Ni)        = "count/kg"
      end if
      if ( iiPDF_Ns > 0 ) then
        variable_names(iiPDF_Ns)        = "Ns"
        variable_descriptions(iiPDF_Ns) = "Snow number concentration"
        variable_units(iiPDF_Ns)        = "count/kg"
      end if
      if ( iiPDF_Ng > 0 ) then
        variable_names(iiPDF_Ng)        = "Ng"
        variable_descriptions(iiPDF_Ng) = "Graupel number concentration"
        variable_units(iiPDF_Ng)        = "count/kg"
      end if

      i = d_variables + 1
      variable_names(i) = "dp1"
      variable_descriptions(i) = "Uniform distribution for the mixture component"

      i = d_variables + 2
      variable_names(i) = "dp2"
      variable_descriptions(i) = "Uniform variate used to determine precipitation!"

      i = d_variables + 3
      variable_names(i) = "X_mixt_comp"
      variable_descriptions(i) = "Mixture component (should be 1 or 2)"

      i = d_variables + 4
      variable_names(i) = "lh_sample_point_weights"
      variable_descriptions(i) = "Weight of each sample point"

      ! Set all the units
      variable_units(:) = "count" ! Unidata units format for a dimensionless quantity

      call open_2D_samples_file( nz, num_samples, i, & ! In
                                 trim( fname_prefix )//"_u", fdir, & ! In
                                 time_initial, stats_tout, stats_zt, & ! In
                                 variable_names(1:i), variable_descriptions(1:i), & ! In
                                 variable_units(1:i), & ! In
                                 uniform_sample_file ) ! In/Out

      deallocate( variable_names, variable_descriptions, variable_units )

    end if

    return
  end subroutine latin_hypercube_2D_output

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_close
! Description:
!   Close a 2D sample file

! References:
!   None
!-------------------------------------------------------------------------------
    use output_2D_samples_module, only: &
      close_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Variable(s)
      uniform_sample_file

    implicit none

    ! ---- Begin Code ----

    if ( l_output_2D_lognormal_dist ) then
      call close_2D_samples_file( lognormal_sample_file )
    end if
    if ( l_output_2D_uniform_dist ) then
      call close_2D_samples_file( uniform_sample_file )
    end if

    return
  end subroutine latin_hypercube_2D_close
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
  elemental function in_mixt_comp_1( X_u_dp1_element, frac )

! Description:
!   Determine if we're in mixture component 1

! References:
!   None
!----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    real(kind=core_rknd), intent(in) :: &
      X_u_dp1_element, & ! Element of X_u telling us which mixture component we're in
      frac               ! The mixture fraction

    logical :: in_mixt_comp_1

    ! ---- Begin Code ----

    if ( X_u_dp1_element < frac ) then
      in_mixt_comp_1 = .true.
    else
      in_mixt_comp_1 = .false.
    end if

    return
  end function in_mixt_comp_1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  elemental function in_precipitation( rnd, precip_frac ) result( l_in_precip )

  ! Description:
  !   Determines if a sample is in precipitation

  ! References:
  !   None
  !-----------------------------------------------------------------------------

    use clubb_precision, only: core_rknd

    implicit none

    ! Input Variables
    real( kind=core_rknd ), intent(in) :: &
      rnd, &         ! Random number between 0 and 1
      precip_frac    ! Precipitation fraction

    ! Output Variable
    logical :: &
      l_in_precip    ! Whether the sample is in precipitation

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( rnd < precip_frac ) then
      l_in_precip = .true.
    else
      l_in_precip = .false.
    end if

    return
  end function in_precipitation
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine compute_arb_overlap( nz, k_lh_start, &
                                  vert_corr, &
                                  X_u_one_var_all_levs )
! Description:
!   Re-computes X_u (uniform sample) for a single variate (e.g. chi) using
!   an arbitrary correlation specified by the input vert_corr variable (which
!   can vary with height).
!   This is an improved algorithm that doesn't require us to convert from a
!   unifrom distribution to a Gaussian distribution and back again.

! References:
!   None
!-------------------------------------------------------------------------------

    use generate_uniform_sample_module, only: &
      rand_uniform_real ! Procedure

    use clubb_precision, only: &
      core_rknd ! Precision

    use constants_clubb, only: &
      zero, &   ! Constants
      one, &
      two, &
      fstderr

    implicit none

    ! Parameter Constants

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels [-]
      k_lh_start   ! Starting k level          [-]

    real(kind=core_rknd), dimension(nz), intent(in) :: &
      vert_corr ! Vertical correlation between k points in range [0,1]   [-]

    ! Output Variables
    real(kind=core_rknd), dimension(nz), intent(inout) :: &
      X_u_one_var_all_levs ! Uniform distribution of 1 variate at all levels [-]
                           ! The value of this variate at k_lh_start should already be populated
                           ! in this array and will be used to fill in the other levels.

    ! Local Variables
    real(kind=core_rknd) :: rand, min_val, half_width, offset, unbounded_point

    integer :: k, kp1, km1 ! Loop iterators

    ! ---- Begin Code ----

    ! Upwards loop
    do k = k_lh_start, nz-1

      kp1 = k+1 ! This is the level we're computing

      if ( vert_corr(kp1) < zero .or. vert_corr(kp1) > one ) then
        write(fstderr,*) "vert_corr(kp1) not between 0 and 1"
        write(fstderr,*) "vert_corr(",kp1,") = ", vert_corr(kp1)
        stop "Fatal error in compute_arb_overlap (SILHS)"
      end if

      half_width = one - vert_corr(kp1)
      min_val = X_u_one_var_all_levs(k) - half_width

      rand = rand_uniform_real( )
      offset = two * half_width * rand

      unbounded_point = min_val + offset

      ! If unbounded_point lies outside the range [0,1],
      ! fold it back so that it is between [0,1]
      if ( unbounded_point > one ) then
        X_u_one_var_all_levs(kp1) = two - unbounded_point
      else if ( unbounded_point < zero ) then
        X_u_one_var_all_levs(kp1) = - unbounded_point
      else
        X_u_one_var_all_levs(kp1) = unbounded_point
      end if

    end do ! k_lh_start..nz-1

    ! Downwards loop
    do k = k_lh_start, 2, -1

      km1 = k-1 ! This is the level we're computing

      if ( vert_corr(km1) < zero .or. vert_corr(km1) > one ) then
        stop "vert_corr(km1) not between 0 and 1"
      end if

      half_width = one - vert_corr(km1)
      min_val = X_u_one_var_all_levs(k) - half_width

      rand = rand_uniform_real( )
      offset = two * half_width * rand

      unbounded_point = min_val + offset

      ! If unbounded_point lies outside the range [0,1],
      ! fold it back so that it is between [0,1]
      if ( unbounded_point > one ) then
        X_u_one_var_all_levs(km1) = two - unbounded_point
      else if ( unbounded_point < zero ) then
        X_u_one_var_all_levs(km1) = - unbounded_point
      else
        X_u_one_var_all_levs(km1) = unbounded_point
      end if

    end do ! k_lh_start..2 decrementing

    return
  end subroutine compute_arb_overlap
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  function compute_vert_corr( nz, delta_zm, Lscale_vert_avg, rcm ) result( vert_corr )
! Description:
!   This function computes the vertical correlation for arbitrary overlap, using
!   density weighted 3pt averaged Lscale and the difference in height levels
!   (delta_zm).
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Constant

    use constants_clubb, only: &
      rc_tol, &
      one

    use parameters_silhs, only: &
      l_max_overlap_in_cloud, & ! Variable(s)
      vert_decorr_coef

    implicit none

    ! External
    intrinsic :: exp

    ! Input Variables
    integer, intent(in) :: &
      nz ! Number of vertical levels  [-]

    real( kind = core_rknd ), intent(in), dimension(nz) :: &
      delta_zm, &        ! Difference between altitudes    [m]
      Lscale_vert_avg, & ! Vertically averaged Lscale      [m]
      rcm                ! Cloud water mixing ratio        [kg/kg]

    ! Output Variable
    real( kind = core_rknd ), dimension(nz) :: &
      vert_corr ! The vertical correlation      [-]

    ! ---- Begin Code ----
    vert_corr(1:nz) = exp( -vert_decorr_coef &
                            * ( delta_zm(1:nz) / Lscale_vert_avg(1:nz) ) )

    if ( l_max_overlap_in_cloud ) then
      where ( rcm > rc_tol )
        vert_corr = one
      end where
    end if

    return
  end function compute_vert_corr

!-------------------------------------------------------------------------------
  subroutine stats_accumulate_lh &
             ( nz, num_samples, d_variables, rho_ds_zt, &
               lh_sample_point_weights, X_nl_all_levs, &
               lh_clipped_vars )

! Description:
!   Clip subcolumns from latin hypercube and create stats for diagnostic
!   purposes.

! References:
!   None
!-------------------------------------------------------------------------------

    use parameters_model, only: hydromet_dim ! Variable

    use grid_class, only: gr

    use stats_variables, only: &
      l_stats_samp, & ! Variable(s)
      ilh_rrm, &
      ilh_Nrm, &
      ilh_rim, &
      ilh_Nim, &
      ilh_rsm, &
      ilh_Nsm, &
      ilh_rgm, &
      ilh_Ngm, &
      ilh_thlm, &
      ilh_rcm, &
      ilh_Ncm, &
      ilh_Ncnm, &
      ilh_rvm, &
      ilh_wm, &
      ilh_cloud_frac, &
      ilh_cloud_frac_unweighted, &
      ilh_chi, &
      ilh_chip2, &
      ilh_eta

    use stats_variables, only: &
      ilh_wp2_zt, &  ! Variable(s)
      ilh_Nrp2_zt, &
      ilh_Ncp2_zt, &
      ilh_Ncnp2_zt, &
      ilh_rcp2_zt, &
      ilh_rtp2_zt, &
      ilh_thlp2_zt, &
      ilh_rrp2_zt, &
      ilh_vwp, &
      ilh_lwp, &
      stats_lh_zt, &
      stats_lh_sfc

    use math_utilities, only: &
      compute_sample_mean, & ! Procedure(s)
      compute_sample_variance

    use stats_type_utilities, only: &
      stat_update_var, & ! Procedure(s)
      stat_update_var_pt

    use array_index, only: &
      iirrm, & ! Variables
      iirsm, & 
      iirim, & 
      iirgm, & 
      iiNrm, &
      iiNsm, &
      iiNim, &
      iiNgm

    use corr_varnce_module, only: &
      iiPDF_chi, & ! Variable(s)
      iiPDF_eta, &
      iiPDF_w, &
      iiPDF_Ncn

    use constants_clubb, only: &
      zero_threshold, &    ! Constant(s)
      zero, &
      one

    use clubb_precision, only: & 
      core_rknd    ! Constant

   use fill_holes, only: &
     vertical_integral ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz                 ! Number of vertical model levels

    real( kind = core_rknd ), intent(in), dimension(nz) :: &
      rho_ds_zt  ! Dry, static density (thermo. levs.) [kg/m^3]

    real( kind = core_rknd ), intent(in), dimension(num_samples) :: &
      lh_sample_point_weights

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    type(lh_clipped_variables_type), dimension(nz,num_samples), intent(in) :: &
      lh_clipped_vars   ! SILHS variables

    ! Local variables
    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      rc_all_points,  & ! Cloud water mixing ratio for all levels        [kg/kg]
      Nc_all_points,  & ! Cloud droplet conc. for all levels              [#/kg]
      Ncn_all_points, & ! Cloud nuclei conc. for all levs.; Nc=Ncn*H(chi) [#/kg]
      rv_all_points     ! Vapor mixing ratio for all levels              [kg/kg]

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real( kind = core_rknd ), dimension(nz) :: &
      lh_thlm,       & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,        & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_Ncm,        & ! Average value of the latin hypercube est. of Nc                [num/kg]
      lh_Ncnm,       & ! Average value of the latin hypercube est. of Ncn               [num/kg]
      lh_rvm,        & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,         & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt,     & ! Average value of the variance of the LH est. of vert. vel.     [m^2/s^2]
      lh_rrp2_zt, & ! Average value of the variance of the LH est. of rr.         [(kg/kg)^2]
      lh_rcp2_zt,    & ! Average value of the variance of the LH est. of rc.            [(kg/kg)^2]
      lh_rtp2_zt,    & ! Average value of the variance of the LH est. of rt             [kg^2/kg^2]
      lh_thlp2_zt,   & ! Average value of the variance of the LH est. of thetal         [K^2]
      lh_Nrp2_zt,    & ! Average value of the variance of the LH est. of Nr.            [#^2/kg^2]
      lh_Ncp2_zt,    & ! Average value of the variance of the LH est. of Nc.            [#^2/kg^2]
      lh_Ncnp2_zt,   & ! Average value of the variance of the LH est. of Ncn.           [#^2/kg^2]
      lh_cloud_frac, & ! Average value of the latin hypercube est. of cloud fraction    [-]
      lh_chi,   & ! Average value of the latin hypercube est. of Mellor's s        [kg/kg]
      lh_eta,   & ! Average value of the latin hypercube est. of Mellor's t        [kg/kg]
      lh_chip2           ! Average value of the variance of the LH est. of chi       [kg/kg]


    real(kind=core_rknd) :: xtmp

    integer :: sample, ivar

    ! ---- Begin Code ----

    rc_all_points = lh_clipped_vars%rc

    if ( l_stats_samp ) then

      ! For all cases where l_lh_importance_sampling is false, the weights
      ! will be 1 (all points equally weighted)

      if ( ilh_rcm + ilh_rcp2_zt + ilh_lwp > 0 ) then
        lh_rcm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      rc_all_points )
        call stat_update_var( ilh_rcm, lh_rcm, stats_lh_zt )

        if ( ilh_lwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 lh_rcm(2:gr%nz), gr%invrs_dzt(2:gr%nz) )

          call stat_update_var_pt( ilh_lwp, 1, xtmp, stats_lh_sfc )
        end if
      end if

      if ( ilh_thlm + ilh_thlp2_zt > 0 ) then
        lh_thlm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                       real( lh_clipped_vars%thl, kind = core_rknd ) )
        call stat_update_var( ilh_thlm, lh_thlm, stats_lh_zt )
      end if

      if ( ilh_rvm + ilh_rtp2_zt > 0 ) then
        rv_all_points = lh_clipped_vars(:,:)%rv
        lh_rvm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      rv_all_points )
        call stat_update_var( ilh_rvm, lh_rvm, stats_lh_zt )
        if ( ilh_vwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 lh_rvm(2:gr%nz), gr%invrs_dzt(2:gr%nz) )

          call stat_update_var_pt( ilh_vwp, 1, xtmp, stats_lh_sfc )
        end if
      end if

      if ( ilh_wm + ilh_wp2_zt > 0 ) then
        lh_wm  = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      real( X_nl_all_levs(:,:,iiPDF_w), kind = core_rknd) )
        call stat_update_var( ilh_wm, lh_wm, stats_lh_zt )
      end if

      if ( ilh_rrm + ilh_Nrm + ilh_rim + ilh_Nim + ilh_rsm + ilh_Nsm + &
           ilh_rgm + ilh_Ngm + ilh_Ncnm + ilh_Ncm > 0 ) then

        lh_hydromet = 0._core_rknd
        call copy_X_nl_into_hydromet_all_pts( nz, d_variables, num_samples, & ! In
                                      X_nl_all_levs, &  ! In
                                      lh_hydromet, & ! In
                                      hydromet_all_points, &  ! Out
                                      Ncn_all_points ) ! Out

        Nc_all_points = lh_clipped_vars%Nc

        ! Get rid of an annoying compiler warning.
        ivar = 1
        ivar = ivar

        forall ( ivar = 1:hydromet_dim )
          lh_hydromet(:,ivar) = compute_sample_mean( nz, num_samples, lh_sample_point_weights,&
                                                     hydromet_all_points(:,:,ivar) )
        end forall ! 1..hydromet_dim

      end if

      if ( ilh_Ncnm > 0 ) then
        lh_Ncnm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                       Ncn_all_points(:,:) )
        call stat_update_var( ilh_Ncnm, lh_Ncnm, stats_lh_zt )
      end if

      if ( ilh_Ncm > 0 ) then
        lh_Ncm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      Nc_all_points(:,:) )
        call stat_update_var( ilh_Ncm, lh_Ncm, stats_lh_zt )
      end if

      ! Latin hypercube estimate of cloud fraction
      if ( ilh_cloud_frac > 0 ) then
        lh_cloud_frac(:) = zero
        do sample = 1, num_samples
          where ( X_nl_all_levs(:,sample,iiPDF_chi) > zero )
            lh_cloud_frac(:) = lh_cloud_frac(:) + one * lh_sample_point_weights(sample)
          end where
        end do
        lh_cloud_frac(:) = lh_cloud_frac(:) / real( num_samples, kind = core_rknd )

        call stat_update_var( ilh_cloud_frac, lh_cloud_frac, stats_lh_zt )
      end if

      ! Sample of lh_cloud_frac that is not weighted
      if ( ilh_cloud_frac_unweighted > 0 ) then
        lh_cloud_frac(:) = zero
        do sample = 1, num_samples
          where ( X_nl_all_levs(:,sample,iiPDF_chi) > zero )
            lh_cloud_frac(:) = lh_cloud_frac(:) + one
          end where
        end do
        lh_cloud_frac(:) = lh_cloud_frac(:) / real( num_samples, kind = core_rknd )

        call stat_update_var( ilh_cloud_frac_unweighted, lh_cloud_frac, stats_lh_zt )
      end if

      ! Latin hypercube estimate of chi
      if ( ilh_chi > 0 ) then
        lh_chi(1:nz) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               X_nl_all_levs(1:nz, 1:num_samples, iiPDF_chi) )
        call stat_update_var( ilh_chi, lh_chi, stats_lh_zt )
      end if

      ! Latin hypercube estimate of variance of chi
      if ( ilh_chip2 > 0 ) then
        lh_chip2(1:nz) &
        = compute_sample_variance( nz, num_samples, &
                                   X_nl_all_levs(:,:,iiPDF_chi), &
                                   lh_sample_point_weights, lh_chi(1:nz) )
        call stat_update_var( ilh_chip2, lh_chip2, stats_lh_zt )
      end if

      ! Latin hypercube estimate of eta
      if ( ilh_eta > 0 ) then
        lh_eta(1:nz) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               X_nl_all_levs(1:nz, 1:num_samples, iiPDF_eta) )

        call stat_update_var( ilh_eta, lh_eta, stats_lh_zt )
      end if

      if ( ilh_wp2_zt > 0 ) then
        ! Compute the variance of vertical velocity
        lh_wp2_zt = compute_sample_variance( nz, num_samples, &
                                             X_nl_all_levs(:,:,iiPDF_w), &
                                             lh_sample_point_weights, lh_wm )
        call stat_update_var( ilh_wp2_zt, lh_wp2_zt, stats_lh_zt )
      end if

      if ( ilh_rcp2_zt  > 0 ) then
        ! Compute the variance of cloud water mixing ratio
        lh_rcp2_zt = compute_sample_variance &
                     ( nz, num_samples, rc_all_points, &
                       lh_sample_point_weights, lh_rcm )
        call stat_update_var( ilh_rcp2_zt, lh_rcp2_zt, stats_lh_zt )
      end if

      if ( ilh_rtp2_zt > 0 ) then
        ! Compute the variance of total water
        lh_rtp2_zt = compute_sample_variance &
                     ( nz, num_samples, &
                       lh_clipped_vars%rt, lh_sample_point_weights, &
                       lh_rvm+lh_rcm )
        call stat_update_var( ilh_rtp2_zt, lh_rtp2_zt, stats_lh_zt )
      end if

      if ( ilh_thlp2_zt > 0 ) then
        ! Compute the variance of liquid potential temperature
        lh_thlp2_zt = compute_sample_variance( nz, num_samples, &
                        lh_clipped_vars%thl, lh_sample_point_weights, &
                        lh_thlm )
        call stat_update_var( ilh_thlp2_zt, lh_thlp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of rain water mixing ratio
      if ( iirrm > 0 .and. ilh_rrp2_zt > 0 ) then
        lh_rrp2_zt = compute_sample_variance &
                        ( nz, num_samples, hydromet_all_points(:,:,iirrm), &
                          lh_sample_point_weights, lh_hydromet(:,iirrm) )
        call stat_update_var( ilh_rrp2_zt, lh_rrp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of cloud nuclei concentration (simplifed)
      if ( iiPDF_Ncn > 0 .and. ilh_Ncnp2_zt > 0 ) then
        lh_Ncnp2_zt = compute_sample_variance &
                      ( nz, num_samples, Ncn_all_points(:,:), &
                        lh_sample_point_weights, lh_Ncnm(:) )
        call stat_update_var( ilh_Ncnp2_zt, lh_Ncnp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of cloud droplet concentration
      if ( ilh_Ncp2_zt > 0 ) then
        lh_Ncp2_zt = compute_sample_variance &
                     ( nz, num_samples, Nc_all_points(:,:), &
                       lh_sample_point_weights, lh_Ncm(:) )
        call stat_update_var( ilh_Ncp2_zt, lh_Ncp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of rain droplet number concentration
      if ( iiNrm > 0 .and. ilh_Nrp2_zt > 0 ) then
        lh_Nrp2_zt = compute_sample_variance( nz, num_samples, hydromet_all_points(:,:,iiNrm),&
                                              lh_sample_point_weights, lh_hydromet(:,iiNrm) )
        call stat_update_var( ilh_Nrp2_zt, lh_Nrp2_zt, stats_lh_zt )
      end if

      ! Averages of points being fed into the microphysics
      ! These are for diagnostic purposes, and are not needed for anything
      if ( iirrm > 0 ) then
        call stat_update_var( ilh_rrm, lh_hydromet(:,iirrm), stats_lh_zt )
      end if
      if ( iiNrm > 0 ) then
        call stat_update_var( ilh_Nrm, lh_hydromet(:,iiNrm), stats_lh_zt )
      end if
      if ( iirim > 0 ) then
        call stat_update_var( ilh_rim, lh_hydromet(:,iirim), stats_lh_zt )
      end if
      if ( iiNim > 0 ) then
        call stat_update_var( ilh_Nim, lh_hydromet(:,iiNim), stats_lh_zt )
      end if
      if ( iirsm > 0 ) then
        call stat_update_var( ilh_rsm, lh_hydromet(:,iirsm), stats_lh_zt )
      end if
      if ( iiNsm > 0 ) then
        call stat_update_var( ilh_Nsm, lh_hydromet(:,iiNsm), stats_lh_zt )
      end if
      if ( iirgm > 0 ) then
        call stat_update_var( ilh_rgm, lh_hydromet(:,iirgm), stats_lh_zt )
      end if
      if ( iiNgm > 0 ) then
        call stat_update_var( ilh_Ngm, lh_hydromet(:,iiNgm), stats_lh_zt )
      end if

    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_lh

  !-----------------------------------------------------------------------
  subroutine stats_accumulate_uniform_lh( nz, num_samples, l_in_precip_all_levs, &
                                          X_mixt_comp_all_levs, X_u_chi_all_levs, pdf_params, &
                                          lh_sample_point_weights, k_lh_start )

  ! Description:
  !   Samples statistics that cannot be deduced from the normal-lognormal
  !   SILHS sample (X_nl_all_levs)

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd            ! Constant

    use stats_type_utilities, only: &
      stat_update_var,   & ! Procedure(s)
      stat_update_var_pt

    use stats_variables, only: &
      l_stats_samp, &      ! Variable(s)
      ilh_precip_frac, &
      ilh_mixt_frac, &
      ilh_precip_frac_unweighted, &
      ilh_mixt_frac_unweighted, &
      ik_lh_start, &
      ilh_samp_frac_category, &
      stats_lh_zt, &
      stats_lh_sfc

    use math_utilities, only: &
      compute_sample_mean ! Procedure

    use constants_clubb, only: &
      one, &              ! Constant(s)
      zero

    use pdf_parameter_module, only: &
      pdf_parameter       ! Type

    use silhs_importance_sample_module, only: &
      importance_category_type, &  ! Type
      num_importance_categories, & ! Constant
      define_importance_categories

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &         ! Number of vertical levels
      num_samples ! Number of SILHS sample points

    logical, dimension(nz,num_samples), intent(in) :: &
      l_in_precip_all_levs ! Boolean variables indicating whether a sample is in
                           ! precipitation at a given height level

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Integers indicating which mixture component a
                           ! sample is in at a given height level

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      X_u_chi_all_levs     ! Uniform value of chi

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params           ! The official PDF parameters!

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! The weight of each sample

    integer, intent(in) :: &
      k_lh_start           ! Vertical level for sampling preferentially within       [-]
                           ! cloud

    ! Local Variables
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories

    real( kind = core_rknd ), dimension(nz) :: &
      lh_precip_frac, &
      lh_mixt_frac

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      int_in_precip, & ! '1' for samples in precipitation, '0' otherwise
      int_mixt_comp    ! '1' for samples in the first PDF component, '0' otherwise

    real( kind = core_rknd ), dimension(num_samples) :: &
      one_weights

    real( kind = core_rknd ), dimension(nz,num_importance_categories) :: &
      lh_samp_frac

    real( kind = core_rknd ) :: &
      cloud_frac_i

    logical :: &
      l_in_cloud, &
      l_in_comp_1

    integer, dimension(num_importance_categories) :: &
      category_counts  ! Count of number of samples in each importance category

    integer :: k, isample, icategory
  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_stats_samp ) then
      ! Estimate of lh_precip_frac
      if ( ilh_precip_frac > 0 ) then
        where ( l_in_precip_all_levs )
          int_in_precip = 1.0_core_rknd
        else where
          int_in_precip = 0.0_core_rknd
        end where
        lh_precip_frac(:) = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                                 int_in_precip )
        call stat_update_var( ilh_precip_frac, lh_precip_frac, stats_lh_zt )
      end if

      ! Unweighted estimate of lh_precip_frac
      if ( ilh_precip_frac_unweighted > 0 ) then
        where ( l_in_precip_all_levs )
          int_in_precip = 1.0_core_rknd
        else where
          int_in_precip = 0.0_core_rknd
        end where
        one_weights = one
        lh_precip_frac(:) = compute_sample_mean( nz, num_samples, one_weights, &
                                                 int_in_precip )
        call stat_update_var( ilh_precip_frac_unweighted, lh_precip_frac, stats_lh_zt )
      end if

      ! Estimate of lh_mixt_frac
      if ( ilh_mixt_frac > 0 ) then
        where ( X_mixt_comp_all_levs == 1 )
          int_mixt_comp = 1.0_core_rknd
        else where
          int_mixt_comp = 0.0_core_rknd
        end where
        lh_mixt_frac(:) = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                               int_mixt_comp )
        call stat_update_var( ilh_mixt_frac, lh_mixt_frac, stats_lh_zt )
      end if

      ! Unweighted estimate of lh_mixt_frac
      if ( ilh_mixt_frac_unweighted > 0 ) then
        where ( X_mixt_comp_all_levs == 1 )
          int_mixt_comp = 1.0_core_rknd
        else where
          int_mixt_comp = 0.0_core_rknd
        end where
        one_weights = one
        lh_mixt_frac(:) = compute_sample_mean( nz, num_samples, one_weights, &
                                               int_mixt_comp )
        call stat_update_var( ilh_mixt_frac_unweighted, lh_mixt_frac, stats_lh_zt )
      end if

      ! k_lh_start is an integer, so it would be more appropriate to sample it
      ! as an integer, but as far as I can tell our current sampling
      ! infrastructure mainly supports sampling real numbers.
      call stat_update_var_pt( ik_lh_start, 1, real( k_lh_start, kind=core_rknd ), stats_lh_sfc )

      if ( allocated( ilh_samp_frac_category ) ) then
        if ( ilh_samp_frac_category(1) > 0 ) then

          importance_categories = define_importance_categories( )

          do k=1, nz
            category_counts(:) = 0

            do isample=1, num_samples

              if ( X_mixt_comp_all_levs(k,isample) == 1 ) then
                l_in_comp_1 = .true.
                cloud_frac_i = pdf_params(k)%cloud_frac_1
              else
                l_in_comp_1 = .false.
                cloud_frac_i = pdf_params(k)%cloud_frac_2
              end if

              l_in_cloud = X_u_chi_all_levs(k,isample) > (one - cloud_frac_i)

              do icategory=1, num_importance_categories
                if ( (l_in_cloud .eqv. importance_categories(icategory)%l_in_cloud) .and. &
                     (l_in_precip_all_levs(k,isample) .eqv. importance_categories(icategory)%&
                                                           l_in_precip) .and. &
                     (l_in_comp_1 .eqv. importance_categories(icategory)%l_in_component_1) ) then

                  category_counts(icategory) = category_counts(icategory) + 1
                  exit

                end if
              end do

            end do ! isample=1, num_samples

            do icategory=1, num_importance_categories
              lh_samp_frac(k,icategory) = real( category_counts(icategory), kind=core_rknd ) / &
                                          real( num_samples, kind=core_rknd )
            end do

          end do ! k=2, nz

          ! Microphysics is not run at lower level
          lh_samp_frac(1,:) = zero

          do icategory=1, num_importance_categories
            call stat_update_var( ilh_samp_frac_category(icategory), lh_samp_frac(:,icategory), &
                                  stats_lh_zt )
          end do ! icategory=1, num_importance_categories

        end if ! ilh_samp_frac_category(1) > 0
      end if ! allocated( ilh_samp_frac_category )

    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_uniform_lh

  !-----------------------------------------------------------------------------
  subroutine copy_X_nl_into_hydromet_all_pts( nz, d_variables, num_samples, &
                                      X_nl_all_levs, &
                                      hydromet, &
                                      hydromet_all_points, &
                                      Ncn_all_points )

  ! Description:
  !   Copy the points from the latin hypercube sample to an array with just the
  !   hydrometeors
  ! References:
  !   None
  !-----------------------------------------------------------------------------
    use parameters_model, only: &
      hydromet_dim ! Variable

    use array_index, only: &
      iirrm, & ! Variables
      iirsm, & 
      iirim, & 
      iirgm, & 
      iiNrm, &
      iiNsm, &
      iiNim, &
      iiNgm

    use corr_varnce_module, only: &
      iiPDF_rr, &
      iiPDF_rs, &
      iiPDF_ri, &
      iiPDF_rg, &
      iiPDF_Nr, &
      iiPDF_Ns, &
      iiPDF_Ng, &
      iiPDF_Ncn, &
      iiPDF_Ni

    use clubb_precision, only: &
      core_rknd

    implicit none

    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      d_variables,   & ! Number of variates
      num_samples    ! Number of calls to microphysics

    real( kind = core_rknd ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim), intent(out) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      Ncn_all_points    ! Cloud nuclei conc. (simplified); Nc=Ncn*H(chi)  [#/kg]

    integer :: sample, ivar

    do sample = 1, num_samples
      ! Copy the sample points into the temporary arrays
      do ivar = 1, hydromet_dim, 1
        if ( ivar == iirrm .and. iiPDF_rr > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rr), kind = core_rknd )

        else if ( ivar == iirsm .and. iiPDF_rs > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rs), kind = core_rknd )

        else if ( ivar == iirim .and. iiPDF_ri > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_ri), kind = core_rknd )

        else if ( ivar == iirgm .and. iiPDF_rg > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rg), kind = core_rknd )

        else if ( ivar == iiNrm .and. iiPDF_Nr > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Nr), kind = core_rknd )

        else if ( ivar == iiNsm .and. iiPDF_Ns > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Ns), kind = core_rknd )

        else if ( ivar == iiNgm .and. iiPDF_Ng > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Ng), kind = core_rknd )

        else if ( ivar == iiNim .and. iiPDF_Ni > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Ni), kind = core_rknd )

        else ! Use the mean field, rather than a sample point
          ! This is the case for hail and graupel in the Morrison microphysics
          ! currently -dschanen 23 March 2010
          hydromet_all_points(:,sample,ivar) = hydromet(:,ivar)

        end if
      end do ! 1..hydromet_dim
      ! Copy Ncn into Ncn all points
      if ( iiPDF_Ncn > 0 ) then
        Ncn_all_points(:,sample) = &
          real( X_nl_all_levs(:,sample,iiPDF_Ncn), kind=core_rknd )
      end if
    end do ! 1..num_samples

    return
  end subroutine copy_X_nl_into_hydromet_all_pts
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  elemental function Ncn_to_Nc( Ncn, chi ) result ( Nc )

  ! Description:
  !   Converts a sample of Ncn to a sample of Nc, where
  !   Nc = Ncn * H(chi)
  !   and H(x) is the Heaviside step function.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd    ! Our awesome generalized precision (constant)

    use constants_clubb, only: &
      zero         ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Ncn,  &  ! Simplified cloud nuclei concentration N_cn  [num/kg]
      chi      ! Extended cloud water mixing ratio           [kg/kg]

    ! Output Variable
    real( kind = core_rknd ) :: &
      Nc       ! Cloud droplet concentration                 [num/kg]

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( chi > zero ) then
      Nc = Ncn
    else
      Nc = zero
    end if

    return
  end function Ncn_to_Nc
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  elemental function chi_to_rc( chi ) result ( rc )

  ! Description:
  !   Converts a sample of chi to a sample of rc, where
  !   rc = chi * H(chi)
  !   and H(x) is the Heaviside step function.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd    ! Our awesome generalized precision (constant)

    use constants_clubb, only: &
      zero         ! Constant

    implicit none

    ! Input Variable
    real( kind = core_rknd ), intent(in) :: &
      chi      ! Extended cloud water mixing ratio           [kg/kg]

    ! Output Variable
    real( kind = core_rknd ) :: &
      rc       ! Cloud water mixing ratio                    [kg/kg]

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( chi > zero ) then
      rc = chi
    else
      rc = zero
    end if

    return
  end function chi_to_rc
  !-----------------------------------------------------------------------


#endif /* SILHS */

end module latin_hypercube_driver_module
