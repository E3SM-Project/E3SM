!-------------------------------------------------------------------------------
! $Id: silhs_importance_sample_module.F90 8171 2016-07-06 21:38:34Z raut@uwm.edu $
!===============================================================================

module silhs_importance_sample_module

  implicit none

  integer, parameter, public :: &
    num_importance_categories = 8 ! Number of importance sampling categories
                                  ! ( e.g., (cloud,precip,comp1) )

  private ! Default scope

  public :: importance_category_type, importance_sampling_driver, define_importance_categories, &
            compute_category_real_probs, generate_strat_uniform_variate, pick_sample_categories, &
            cloud_weighted_sampling_driver, determine_sample_categories

  type importance_category_type

    logical :: &
      l_in_cloud,  &
      l_in_precip, &
      l_in_component_1

  end type importance_category_type

  contains

!-----------------------------------------------------------------------
  subroutine importance_sampling_driver &
             ( num_samples, pdf_params, hydromet_pdf_params,        &
               X_u_chi_one_lev, X_u_dp1_one_lev, X_u_dp2_one_lev,   &
               lh_sample_point_weights )

  ! Description:
  !   Applies importance sampling to a single vertical level !

  ! References:
  !   clubb:ticket:736 !
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd         ! Constant(s)

    use constants_clubb, only: &
      fstderr           ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter     ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use error_code, only: &
      clubb_at_least_debug_level  ! Procedure

    use parameters_silhs, only: &
      eight_cluster_allocation_opt, & ! Constant(s)
      four_cluster_allocation_opt, &
      two_cluster_cp_nocp_opt, &
      l_lh_clustered_sampling, & ! Variable(s)
      l_lh_limit_weights, &
      cluster_allocation_strategy, &
      l_lh_normalize_weights

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples       ! Number of SILHS sample points

    type(pdf_parameter), intent(in) :: &
      pdf_params

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(num_samples), intent(inout) :: &
      X_u_chi_one_lev,  & ! These uniform variates are scaled by the importance
      X_u_dp1_one_lev,  & ! sampling process.
      X_u_dp2_one_lev

    ! Output Variables
    real( kind = core_rknd ), dimension(num_samples), intent(out) :: &
      lh_sample_point_weights

    ! Local Variables
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories ! A vector containing the different importance categories

    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_real_probs,       & ! The real PDF probabilities for each category
      category_prescribed_probs, & ! Prescribed probability for each category
      category_sample_weights      ! Sample weight for each category

    real( kind = core_rknd ), dimension(num_samples) :: &
      rand_vect

    real( kind = core_rknd ) :: weights_sum

    integer, dimension(num_samples) :: &
      int_sample_category  ! An integer for each sample corresponding to the
                           ! category picked for the sample

    integer :: sample
    logical :: l_error

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    importance_categories = define_importance_categories( )

    category_real_probs = compute_category_real_probs &
                          ( importance_categories, pdf_params, hydromet_pdf_params )

    if ( l_lh_clustered_sampling ) then

      ! A cluster allocation strategy is selected based on the value of the
      ! integer parameter.
      select case ( cluster_allocation_strategy )
      case ( eight_cluster_allocation_opt )
        category_prescribed_probs = eight_cluster_allocation &
                                    ( importance_categories, category_real_probs )
      case ( four_cluster_allocation_opt )
        category_prescribed_probs = four_cluster_no_precip &
                                    ( importance_categories, category_real_probs )
      case ( two_cluster_cp_nocp_opt )
        category_prescribed_probs = two_cluster_cp_nocp &
                                    ( importance_categories, category_real_probs )
      case default
        write(fstderr,*) "Unsupported allocation strategy:", cluster_allocation_strategy
        stop "Fatal error in importance_sampling_driver"
      end select

      if ( l_lh_limit_weights ) then
        call limit_category_weights( category_real_probs, & ! In
                                     category_prescribed_probs ) ! In/Out
      end if

    else ! .not. l_lh_clustered_sampling

      ! Importance sampling strategy that places half of all sample points in cloud!
      category_prescribed_probs = cloud_importance_sampling &
                                  ( importance_categories, category_real_probs, &
                                    pdf_params )

    end if ! l_lh_clustered_sampling

    ! Compute weight of each sample category
    category_sample_weights = compute_category_sample_weights &
                              ( category_real_probs, category_prescribed_probs )

    ! Generate a stratified sample to be used to pick categories for the samples
    rand_vect = generate_strat_uniform_variate( num_samples )

    ! Pick the sample points!
    int_sample_category = pick_sample_categories( num_samples, category_prescribed_probs, &
                                                  rand_vect )

    ! Loop over sample points
    do sample=1, num_samples

      ! Scale and translate point to reside in its respective cateogry
      call scale_sample_to_category &
           ( importance_categories(int_sample_category(sample)), & ! In
             pdf_params, hydromet_pdf_params, & ! In
             X_u_chi_one_lev(sample), X_u_dp1_one_lev(sample), X_u_dp2_one_lev(sample) ) ! In/Out

      ! Pick a weight for the sample point
      lh_sample_point_weights(sample) = &
        category_sample_weights(int_sample_category(sample))

    end do ! sample=1, num_samples

    ! Normalize weights (if enabled)
    if ( l_lh_normalize_weights ) then
      weights_sum = sum( lh_sample_point_weights(1:num_samples) )

      lh_sample_point_weights(1:num_samples) = &
        lh_sample_point_weights(1:num_samples) * real( num_samples, kind = core_rknd ) / weights_sum
    end if

    if ( clubb_at_least_debug_level( 2 ) ) then

      call importance_sampling_assertions &
           ( num_samples, importance_categories, category_real_probs, & ! In
             category_prescribed_probs, category_sample_weights, X_u_chi_one_lev, & ! In
             X_u_dp1_one_lev, X_u_dp2_one_lev, lh_sample_point_weights, int_sample_category, & ! In
             pdf_params, hydromet_pdf_params, & ! In
             l_error ) ! Out

      if ( l_error ) then
        stop "Fatal error in importance_sampling_driver"
      end if ! l_error

    end if

    return
  end subroutine importance_sampling_driver
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function define_importance_categories( ) &

  result( importance_categories )

  ! Description:
  !   Creates a vector of size num_importance_categories that defines the eight
  !   importance sampling categories.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Output Variable
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories ! A vector containing the different importance categories

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    importance_categories(1)%l_in_cloud       = .true.
    importance_categories(1)%l_in_precip      = .true.
    importance_categories(1)%l_in_component_1 = .true.

    importance_categories(2)%l_in_cloud       = .true.
    importance_categories(2)%l_in_precip      = .true.
    importance_categories(2)%l_in_component_1 = .false.

    importance_categories(3)%l_in_cloud       = .false.
    importance_categories(3)%l_in_precip      = .true.
    importance_categories(3)%l_in_component_1 = .true.

    importance_categories(4)%l_in_cloud       = .false.
    importance_categories(4)%l_in_precip      = .true.
    importance_categories(4)%l_in_component_1 = .false.

    importance_categories(5)%l_in_cloud       = .true.
    importance_categories(5)%l_in_precip      = .false.
    importance_categories(5)%l_in_component_1 = .true.

    importance_categories(6)%l_in_cloud       = .true.
    importance_categories(6)%l_in_precip      = .false.
    importance_categories(6)%l_in_component_1 = .false.

    importance_categories(7)%l_in_cloud       = .false.
    importance_categories(7)%l_in_precip      = .false.
    importance_categories(7)%l_in_component_1 = .true.

    importance_categories(8)%l_in_cloud       = .false.
    importance_categories(8)%l_in_precip      = .false.
    importance_categories(8)%l_in_component_1 = .false.

    return
  end function define_importance_categories
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_category_real_probs( importance_categories, pdf_params, &
                                        hydromet_pdf_params ) &

  result( category_real_probs )

  ! Description:
  !   Computes the real PDF probability associated with each importance sampling
  !   category.
  !
  !   For example, if a category is in cloud, out of precipitation, and in mixture
  !   component two, then without importance sampling, the probability that a
  !   point will appear in that category is:
  !   P(cloud,noprecip,comp2) = cloud_frac_2 * (1-precip_frac_2) * (1-mixt_frac)

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one    ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter           ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter  ! Type

    implicit none

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories  ! A list of importance categories

    type(pdf_parameter), intent(in) :: &
      pdf_params             ! The PDF parameters!

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params    ! The hydrometeor PDF parameters!

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_real_probs ! The real PDF probabilities for each category

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac_i,        &
      precip_frac_i,       &
      cloud_factor,        &
      precip_factor,       &
      component_factor

    integer :: icategory

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do icategory = 1, num_importance_categories

      ! Determine component of category
      if ( importance_categories(icategory)%l_in_component_1 ) then
        cloud_frac_i     = pdf_params%cloud_frac_1
        precip_frac_i    = hydromet_pdf_params%precip_frac_1
        component_factor = pdf_params%mixt_frac
      else
        cloud_frac_i     = pdf_params%cloud_frac_2
        precip_frac_i    = hydromet_pdf_params%precip_frac_2
        component_factor = (one-pdf_params%mixt_frac)
      end if

      ! Determine cloud factor
      if ( importance_categories(icategory)%l_in_cloud ) then
        cloud_factor = cloud_frac_i
      else
        cloud_factor = (one-cloud_frac_i)
      end if

      ! Determine precip factor
      if ( importance_categories(icategory)%l_in_precip ) then
        precip_factor = precip_frac_i
      else
        precip_factor = (one-precip_frac_i)
      end if

      ! Compute the category probability
      category_real_probs(icategory) = component_factor * cloud_factor * precip_factor

    end do ! icategory = 1, num_importance_categories

    return
  end function compute_category_real_probs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_category_sample_weights( category_real_probs, category_prescribed_probs ) &

  result( category_sample_weights )

  ! Description:
  !   Compute the sample point weights for a sample point in each category based
  !   on the PDF probability and the modified probability from importance
  !   sampling

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd     ! Constant

    use constants_clubb, only: &
      zero, &       ! Constant
      unused_var

    implicit none

    ! Input Variables

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs,  &   ! The actual PDF probability of each category
      category_prescribed_probs ! The modified probability of each category due to
                                ! importance sampling

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_sample_weights   ! Sample weight for each category

    ! Local Variable
    integer :: icategory

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do icategory=1, num_importance_categories

      if ( category_prescribed_probs(icategory) == zero ) then
        ! If a category has no probability of being sampled, then its weight is irrevelant.
        category_sample_weights(icategory) = unused_var
      else
        category_sample_weights(icategory) = &
          category_real_probs(icategory) / category_prescribed_probs(icategory)
      end if

    end do

    return
  end function compute_category_sample_weights
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine limit_category_weights( category_real_probs, category_prescribed_probs )

  ! Description:
  !   Modifies the category prescribed probabilities such that an
  !   invariant maximum weight is achieved. The weight for each category
  !   j is equal to category_real_probs(j) / category_prescribed_probs(j).
  !   In this subroutine, category_prescribed_probs(j) is increased in
  !   each necessary category such that no category has a weight that
  !   exceeds the maximum.

  ! References:
  !   none
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd      ! Constant

    use constants_clubb, only: &
      zero, &        ! Constant(s)
      fstderr

    use error_code, only: &
      clubb_at_least_debug_level  ! Procedure

    implicit none

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      max_weight                    = 2.0_core_rknd, &
      real_prob_thresh_for_transfer = 1.0e-8_core_rknd

    ! Input Variables
    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs    ! The PDF probability for each importance category

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(num_importance_categories), intent(inout) :: &
      category_prescribed_probs  ! Prescribed probability for each category; these will
                                 ! be modified to ensure the minimum weight.

    ! Local Variables
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      min_presc_probs, &         ! Minimum prescribed probabilities necessary for each category
      min_presc_prob_diff        ! category_prescribed_probs(:) - min_presc_probs(:)

    real( kind = core_rknd ) :: &
      total_diff_under, &
      total_diff_over, &
      weight

    integer :: icategory
  !-----------------------------------------------------------------------

    !----- Begin Code -----
    total_diff_under = zero

    do icategory=1, num_importance_categories
      if ( category_real_probs(icategory) >= real_prob_thresh_for_transfer ) then
        min_presc_probs(icategory) = category_real_probs(icategory) / max_weight
      else
        min_presc_probs(icategory) = category_real_probs(icategory)
      end if

      min_presc_prob_diff(icategory) = category_prescribed_probs(icategory) - &
                                       min_presc_probs(icategory)
      if ( min_presc_prob_diff(icategory) < zero ) then
        total_diff_under = total_diff_under + abs( min_presc_prob_diff(icategory) )
      end if
    end do ! icategory=1, num_importance_categories

    if ( total_diff_under > zero ) then

      ! Find the sum of the positive differences in categories
      total_diff_over = zero
      do icategory=1, num_importance_categories
        if ( min_presc_prob_diff(icategory) > zero ) then
          total_diff_over = total_diff_over + min_presc_prob_diff(icategory)
        end if
      end do ! icategory=1, num_importance_categories

      if ( total_diff_under > total_diff_over ) then
        ! The prescribed probabilities cannot be adjusted to achieve the maximum
        ! weight. This could happen if the maximum weight is too low.
        write(fstderr,*) "The sample point weights could not be limited to the &
                         &maximum value."
        stop "Fatal error in limit_category_weights"
      end if

      ! Adjust the prescribed probabilities to achieve the minimum.
      do icategory=1, num_importance_categories
        if ( min_presc_prob_diff(icategory) > zero ) then
          category_prescribed_probs(icategory) = category_prescribed_probs(icategory) - &
              ( min_presc_prob_diff(icategory) / total_diff_over ) * total_diff_under
        else if ( min_presc_prob_diff(icategory) < zero ) then
          category_prescribed_probs(icategory) = min_presc_probs(icategory)
        end if ! min_presc_prob_diff(icategory) > zero
      end do ! icategory=1, num_importance_categories

    end if ! total_diff_under > zero

    ! As an assertion check, make sure this subroutine actually
    ! performed its task successfully.
    if ( clubb_at_least_debug_level( 2 ) ) then
      do icategory=1, num_importance_categories
        if ( category_prescribed_probs(icategory) > zero .and. &
             category_real_probs(icategory) >= real_prob_thresh_for_transfer ) then
          weight = category_real_probs(icategory) / category_prescribed_probs(icategory)
          if ( weight > max_weight ) then
            write(fstderr,*) "In limit_category_weights, a weight was not limited."
            write(fstderr,*) "category_real_prob = ", category_real_probs(icategory)
            write(fstderr,*) "category_prescribed_prob = ", category_prescribed_probs(icategory)
            write(fstderr,*) "weight = ", weight
            write(fstderr,*) "min_presc_probs = ", min_presc_probs(icategory)
            write(fstderr,*) "min_presc_prob_diff = ", min_presc_prob_diff(icategory)
            stop "Fatal error in limit_category_weights"
          end if
        end if ! category_prescribed_probs(icategory) > zero
      end do ! icategory=1, num_importance_categories
    end if ! clubb_at_least_debug_level( 2 )

    return
  end subroutine limit_category_weights
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function pick_sample_categories( num_samples, category_prescribed_probs, rand_vect ) &

  result( int_sample_category )

  ! Description:
  !   Picks a category for each sample point, based on the given probabilities,
  !   such that the distribution of categories of the sample points
  !   approximates the probabilities that are given.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd      ! Konstant

    use constants_clubb, only: &
      fstderr, &     ! Constant(s)
      one, &
      zero

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples          ! Number of sample points to be picked

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_prescribed_probs ! Prescribed probability for each category

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      rand_vect            ! A sample of num_samples values from the uniform distribution
                           ! in the range (0,1). This will be used to pick the
                           ! categories.

    ! Output Variable
    integer, dimension(num_samples) :: &
      int_sample_category  ! An integer for each sample corresponding to the
                           ! category picked for the sample

    ! Local Variables
    integer :: sample,category      ! Looping variable(s)

    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_cumulative_probs

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    !--------------------------------------------------------------------------
    ! In order to facilitate picking categories for the sample points, a new
    ! array, category_cumulative_probs, is created.
    !
    ! Each element of category_cumulative_probs is simply the sum of the
    ! probabilities of all the previous categories. For example,
    !
    ! category_cumulative_probs(1) = 0
    ! category_cumulative_probs(2) = category_prescribed_probs(1)
    ! category_cumulative_probs(3) = category_prescribed_probs(1) + category_prescribed_probs(2)
    ! ...
    ! category_cumulative_probs(num_importance_categories) =
    !             category_prescribed_probs(1) + category_prescribed_probs(2) + ... +
    !             category_prescribed_probs(num_importance_categories-1)
    !--------------------------------------------------------------------------
    category_cumulative_probs(1) = zero
    do category=2, num_importance_categories
      category_cumulative_probs(category) = category_cumulative_probs(category-1) + &
                                            category_prescribed_probs(category-1)
    end do

    !--------------------------------------------------------------------------
    ! Pick categories based on the values of rand_vect.
    !--------------------------------------------------------------------------
    do sample=1, num_samples
      ! Initialize int_sample_category(sample) for error checking purposes.
      int_sample_category(sample) = 0
      do category=1, num_importance_categories

        if ( category < num_importance_categories ) then

          if ( rand_vect(sample) >= category_cumulative_probs(category) .and. &
               rand_vect(sample) <  category_cumulative_probs(category+1) ) then

            int_sample_category(sample) = category
            exit   ! Break out of the loop over categories, since we have found the category

          end if

        else if ( category == num_importance_categories ) then

          ! If the random number is greater than category_cumulative_probs(num_imp_categories-1)
          ! and less than 1, then it belongs in the last category
          if ( rand_vect(sample) >= category_cumulative_probs(category) .and. &
               rand_vect(sample) < one ) then
            int_sample_category(sample) = category
          end if

        end if ! category < num_importance_categories

      end do ! category=1, num_importance_categories-1

      ! We should have picked a category by now.
      if ( int_sample_category(sample) == 0 ) then
        write(fstderr,*) "Invalid rand_vect number in pick_sample_categories"
        write(fstderr,*) "rand_vect(sample) = ", rand_vect(sample)
        stop "Fatal error"
      end if

    end do ! sample=1, num_samples

    return
  end function pick_sample_categories
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine scale_sample_to_category( category, pdf_params, hydromet_pdf_params, &
                                       X_u_chi, X_u_dp1, X_u_dp2 )

  ! Description:
  !   Scale and transpose a sample point to reside in the specified category

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd          ! Konstant

    use constants_clubb, only: &
      one                ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter

    implicit none

    ! Input Variables
    type(importance_category_type), intent(in) :: &
      category             ! Scale the sample point to reside in this category

    type(pdf_parameter), intent(in) :: &
      pdf_params

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params

    ! Input/Output Variable
    ! These uniform samples, upon input, are uniformly distributed in the
    ! range (0,1).
    real( kind = core_rknd ), intent(inout) :: &
      X_u_chi, & ! Uniform sample of extend cloud water mixing ratio
      X_u_dp1, & ! Uniform sample of the d+1 variate
      X_u_dp2    ! Uniform sample of the d+2 variate

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac_i, & 
      precip_frac_i, &
      mixt_frac

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    mixt_frac = pdf_params%mixt_frac

    !--------------------------------------------------------
    ! Scale dp1 variate to be in component 1 or 2
    !--------------------------------------------------------
    if ( category%l_in_component_1 ) then

      ! Samples in component 1 have a dp1 variate that satisfies
      ! 0 < X_u_dp1 < mixt_frac.
      ! Scale X_u_dp1 to lie in (0,mixt_frac)
      X_u_dp1 = X_u_dp1 * mixt_frac

      ! Choose appropriate component cloud and precipitation fractions
      cloud_frac_i  = pdf_params%cloud_frac_1
      precip_frac_i = hydromet_pdf_params%precip_frac_1

    else  ! in component 2

      ! Scale and translate X_u_dp1 to lie in (mixt_frac, 1)
      X_u_dp1 = X_u_dp1 * (one - mixt_frac) + mixt_frac

      ! Choose appropriate component cloud and precipitation fractions
      cloud_frac_i  = pdf_params%cloud_frac_2
      precip_frac_i = hydromet_pdf_params%precip_frac_2

    end if ! category%l_in_component_1

    !--------------------------------------------------------
    ! Scale dp2 variate to be in or out of precipitation
    !--------------------------------------------------------
    if ( category%l_in_precip ) then
      ! Samples in precipitation have a dp2 variate that satisfies
      ! 0 < X_u_dp2 < precip_frac_i
      ! Scale X_u_dp2 to lie in (0,precip_frac_i)
      X_u_dp2 = X_u_dp2 * precip_frac_i
    else
      ! Scale and translate X_u_dp2 to lie in (precip_frac_i,1)
      X_u_dp2 = X_u_dp2 * (one - precip_frac_i) + precip_frac_i
    end if

    !--------------------------------------------------------
    ! Scale chi variate to be in or out of cloud
    !--------------------------------------------------------
    if ( category%l_in_cloud ) then
      ! Samples in cloud have a chi variate that satisfies
      ! (1.0 - cloud_frac_i) < X_u_chi < 1
      ! Scale and translate X_u_chi to lie in ( (1 - cloud_frac_i) , 1 )
      X_u_chi = X_u_chi * cloud_frac_i + (one - cloud_frac_i)
    else
      ! Scale X_u_chi to lie in ( 0, (1 - cloud_frac_i) )
      X_u_chi = X_u_chi * (one - cloud_frac_i)
    end if

    return
  end subroutine scale_sample_to_category
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function two_cluster_cp_nocp( importance_categories, category_real_probs ) &

  result( category_prescribed_probs )

  ! Description:
  !   Clusters importance categories into two clusters: categories that contain
  !   either cloud or precipitation (or both), and clusters that contain
  !   neither.

  ! References:
  !   clubb:ticket:740
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd      ! Constant

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    implicit none

    integer, parameter :: &
      num_clusters = 2, &
      max_num_categories_in_cluster = 6

    integer, parameter :: &
      iicld_or_precip = 1, &
      iinocld_precip = 2

    real( kind = core_rknd ), parameter :: &
      prob_cld_or_precip = 1.0_core_rknd, &
      prob_nocld_precip  = 0.0_core_rknd

    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories   ! A list of importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs     ! The real probability for each category

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs ! The prescribed probability for each category

    ! Local Variables
    integer, dimension(num_clusters,max_num_categories_in_cluster) :: &
      cluster_categories

    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_prescribed_probs

    integer, dimension(num_clusters) :: &
      num_categories_in_cluster

    integer :: icategory

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    num_categories_in_cluster(:) = 0
    cluster_categories(:,:) = 0
    do icategory=1, num_importance_categories
      if ( importance_categories(icategory)%l_in_cloud .or. &
           importance_categories(icategory)%l_in_precip ) then
        num_categories_in_cluster(iicld_or_precip) = &
          num_categories_in_cluster(iicld_or_precip) + 1
        cluster_categories(iicld_or_precip,num_categories_in_cluster(iicld_or_precip)) = icategory
      else
        num_categories_in_cluster(iinocld_precip) = num_categories_in_cluster(iinocld_precip) + 1
        cluster_categories(iinocld_precip,num_categories_in_cluster(iinocld_precip)) = icategory
      end if
    end do ! icategory=1, num_importance_categories

    if ( clubb_at_least_debug_level( 2 ) ) then
      if ( num_categories_in_cluster(iinocld_precip) /= 2 .or. &
           num_categories_in_cluster(iicld_or_precip) /= 6 ) then
        stop "Invalid categories in two_cluster_cp_nocp"
      end if
    end if

    cluster_prescribed_probs(iicld_or_precip) = prob_cld_or_precip
    cluster_prescribed_probs(iinocld_precip)  = prob_nocld_precip

    category_prescribed_probs = compute_clust_category_probs &
      ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
        num_categories_in_cluster, cluster_categories, cluster_prescribed_probs )

    return
  end function two_cluster_cp_nocp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function eight_cluster_allocation( importance_categories, category_real_probs ) &

  result( category_prescribed_probs )

  ! Description:
  !   Clusters importance categories such that each of the eight importance
  !   categories has its own cluster. Effectively, there are no clusters.

  ! References:
  !   clubb:ticket:752
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd     ! Constant

    use parameters_silhs, only: &
      eight_cluster_presc_probs   ! Variable(s)

    implicit none

    ! Local Constants
    integer, parameter :: &
      num_clusters = 8, &
      max_num_categories_in_cluster = 1

    integer, dimension(num_clusters), parameter :: &
      num_categories_in_cluster = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories   ! A list of importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs     ! The real probability for each category

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs ! The prescribed probability for each category

    ! Local Variables
    integer, dimension(num_clusters,max_num_categories_in_cluster) :: &
      cluster_categories

    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_prescribed_probs

    logical :: l_in_cloud, l_in_component_1, l_in_precip

    integer :: icategory

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    do icategory=1, num_importance_categories

      cluster_categories(icategory,1) = icategory

      l_in_cloud       = importance_categories(icategory)%l_in_cloud
      l_in_component_1 = importance_categories(icategory)%l_in_component_1
      l_in_precip      = importance_categories(icategory)%l_in_precip

      if ( l_in_cloud .and. l_in_precip .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%cloud_precip_comp1

      else if ( l_in_cloud .and. l_in_precip .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%cloud_precip_comp2

      else if ( (.not. l_in_cloud) .and. l_in_precip .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%nocloud_precip_comp1

      else if ( (.not. l_in_cloud) .and. l_in_precip .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%nocloud_precip_comp2

      else if ( l_in_cloud .and. (.not. l_in_precip) .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%cloud_noprecip_comp1

      else if ( l_in_cloud .and. (.not. l_in_precip) .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%cloud_noprecip_comp2

      else if ( (.not. l_in_cloud) .and. (.not. l_in_precip) .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%nocloud_noprecip_comp1

      else if ( (.not. l_in_cloud) .and. (.not. l_in_precip) .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = eight_cluster_presc_probs%nocloud_noprecip_comp2

      else
        stop "Invalid category in eight_cluster_allocation"
      end if

    end do ! icategory=1, num_importance_categories

    category_prescribed_probs = compute_clust_category_probs &
      ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
        num_categories_in_cluster, cluster_categories, cluster_prescribed_probs )

    return
  end function eight_cluster_allocation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function four_cluster_no_precip( importance_categories, category_real_probs ) &

  result( category_prescribed_probs )

  ! Description:
  !   Clusters categories into four clusters for the four combinations of
  !   cloud/no cloud and comp 1/comp 2. Precip fraction is effectively
  !   ignored

  ! References:
  !   clubb:ticket:752
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd     ! Constant

    use constants_clubb, only: &
      fstderr       ! Constant

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    implicit none

    ! Local Constants
    integer, parameter :: &
      num_clusters = 4, &
      max_num_categories_in_cluster = 2

    integer, parameter :: &
      iicld_comp1  = 1, &
      iicld_comp2  = 2, &
      iincld_comp1 = 3, &
      iincld_comp2 = 4

    !!! Prescribed probability definitions
    real( kind = core_rknd ), parameter :: &
      cloud_comp1      = 0.30_core_rknd, &
      cloud_comp2      = 0.30_core_rknd, &
      nocloud_comp1    = 0.20_core_rknd, &
      nocloud_comp2    = 0.20_core_rknd

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories   ! A list of importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs     ! The real probability for each category

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs ! The prescribed probability for each category

    ! Local Variables
    integer, dimension(num_clusters,max_num_categories_in_cluster) :: &
      cluster_categories

    integer, dimension(num_clusters) :: &
      num_categories_in_cluster

    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_prescribed_probs

    logical :: l_in_cloud, l_in_component_1

    integer :: icategory

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    num_categories_in_cluster(:) = 0

    do icategory=1, num_importance_categories

      l_in_cloud       = importance_categories(icategory)%l_in_cloud
      l_in_component_1 = importance_categories(icategory)%l_in_component_1

      if ( l_in_cloud .and. l_in_component_1 ) then
        cluster_prescribed_probs(iicld_comp1) = cloud_comp1
        num_categories_in_cluster(iicld_comp1) = num_categories_in_cluster(iicld_comp1) + 1
        cluster_categories(iicld_comp1,num_categories_in_cluster(iicld_comp1)) = icategory
        
      else if ( l_in_cloud .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(iicld_comp2) = cloud_comp2
        num_categories_in_cluster(iicld_comp2) = num_categories_in_cluster(iicld_comp2) + 1
        cluster_categories(iicld_comp2,num_categories_in_cluster(iicld_comp2)) = icategory

      else if ( (.not. l_in_cloud) .and. l_in_component_1 ) then
        cluster_prescribed_probs(iincld_comp1) = nocloud_comp1
        num_categories_in_cluster(iincld_comp1) = num_categories_in_cluster(iincld_comp1) + 1
        cluster_categories(iincld_comp1,num_categories_in_cluster(iincld_comp1)) = icategory

      else if ( (.not. l_in_cloud) .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(iincld_comp2) = nocloud_comp2
        num_categories_in_cluster(iincld_comp2) = num_categories_in_cluster(iincld_comp2) + 1
        cluster_categories(iincld_comp2,num_categories_in_cluster(iincld_comp2)) = icategory

      else
        stop "Invalid category in four_cluster_no_precip"
      end if

    end do ! icategory=1, num_importance_categories

    if ( clubb_at_least_debug_level( 2 ) ) then
      if ( any( num_categories_in_cluster /= 2 ) ) then
        write(fstderr,*) "Not all clusters have two categories"
        stop "Fatal error in four_cluster_no_precip"
      end if
    end if

    category_prescribed_probs = compute_clust_category_probs &
      ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
        num_categories_in_cluster, cluster_categories, cluster_prescribed_probs )

    return
  end function four_cluster_no_precip
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_clust_category_probs &
           ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
             num_categories_in_cluster, cluster_categories, cluster_fractions ) &

  result( category_prescribed_probs )

  ! Description:
  !   Calls clust_cat_probs_frm_presc_prb or clust_cat_probs_frm_var_fracs!

  ! References:
  !   clubb:ticket:740
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd       ! Constant

    use parameters_silhs, only: &
      l_lh_var_frac

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs           ! The real probability for each category

    integer, intent(in) :: &
      num_clusters, &                   ! The number of clusters to sample from
      max_num_categories_in_cluster    ! The max number of categories in each cluster

    integer, dimension(num_clusters), intent(in) :: &
      num_categories_in_cluster        ! The number of categories in each cluster

    integer, dimension(num_clusters,max_num_categories_in_cluster), intent(in) :: &
      cluster_categories            ! An integer matrix containing indices corresponding
                                    ! to the members of the clusters

    real( kind = core_rknd ), dimension(num_clusters), intent(in) :: &
      cluster_fractions             ! Prescribed fraction of some sort for each cluster

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs     ! Resulting prescribed probability for each individual category

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    if ( l_lh_var_frac ) then
      category_prescribed_probs = clust_cat_probs_frm_var_fracs &
           ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
             num_categories_in_cluster, cluster_categories, cluster_fractions )
    else
      category_prescribed_probs = clust_cat_probs_frm_presc_prb &
           ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
             num_categories_in_cluster, cluster_categories, cluster_fractions )
    end if

    return
  end function compute_clust_category_probs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function clust_cat_probs_frm_var_fracs &
           ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
             num_categories_in_cluster, cluster_categories, cluster_variance_fractions ) &

  result( category_prescribed_probs )

  ! Description:
  !   This is a generalized algorithm that takes as input a set of "clusters"
  !   of the importance categories and a variance fraction      for each
  !   cluster, and computes the prescribed probabilities for each category.

  ! References:
  !   clubb:ticket:740
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd       ! Constant

    use constants_clubb, only: &
      zero            ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs           ! The real probability for each category

    integer, intent(in) :: &
      num_clusters, &                   ! The number of clusters to sample from
      max_num_categories_in_cluster    ! The max number of categories in each cluster

    integer, dimension(num_clusters), intent(in) :: &
      num_categories_in_cluster        ! The number of categories in each cluster

    integer, dimension(num_clusters,max_num_categories_in_cluster), intent(in) :: &
      cluster_categories            ! An integer matrix containing indices corresponding
                                    ! to the members of the clusters

    real( kind = core_rknd ), dimension(num_clusters), intent(in) :: &
      cluster_variance_fractions    ! Prescribed variance fraction for each cluster

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs     ! Resulting prescribed probability for each individual category

    ! Local Variables
    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_real_probs, &         ! Total PDF probability for each cluster
      cluster_prescribed_probs      ! Total prescribed probability for each cluster

    real( kind = core_rknd ) :: &
      pdf_prob_var_frac_prod_sum

    integer :: icluster, icategory, cat_idx

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    ! Compute the total PDF probability for each cluster.
    cluster_real_probs(:) = zero
    do icluster=1, num_clusters
      do icategory=1, num_categories_in_cluster(icluster)
        cluster_real_probs(icluster) = cluster_real_probs(icluster) + &
            category_real_probs(cluster_categories(icluster,icategory))
      end do
    end do

    ! Compute the sum of p_j * f_j across clusters
    pdf_prob_var_frac_prod_sum = sum( cluster_real_probs(:) * cluster_variance_fractions(:) )

    ! Compute the prescribed probability for each cluster!
    if ( pdf_prob_var_frac_prod_sum == zero ) then
      ! No variance prescribed in clusters with non-zero PDF probability!
      ! Fall back to no importance sampling!
      cluster_prescribed_probs(:) = cluster_real_probs(:)
    else
      cluster_prescribed_probs(:) = cluster_real_probs(:) * cluster_variance_fractions(:) / &
                                    pdf_prob_var_frac_prod_sum
    end if

    ! Now, split into categories
    do icluster=1, num_clusters
      do icategory=1, num_categories_in_cluster(icluster)
        cat_idx = cluster_categories(icluster,icategory)
        ! Scale category probability based on the cluster probability
        if ( cluster_real_probs(icluster) == zero ) then
          category_prescribed_probs(cat_idx) = zero
        else
          category_prescribed_probs(cat_idx) = ( category_real_probs(cat_idx) / &
            cluster_real_probs(icluster) ) * cluster_prescribed_probs(icluster)
        end if
      end do ! icategory=1, num_categories_in_cluster(icluster)
    end do ! icluster=1, num_clusters


    return
  end function clust_cat_probs_frm_var_fracs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function clust_cat_probs_frm_presc_prb &
           ( category_real_probs, num_clusters, max_num_categories_in_cluster, &
             num_categories_in_cluster, cluster_categories, cluster_prescribed_probs ) &

  result( category_prescribed_probs )

  ! Description:
  !   This is a generalized algorithm that takes as input a set of "clusters"
  !   of the importance categories and a prescribed probability for each
  !   cluster, and computes the prescribed probabilities for each category
  !   such that the sum of the prescribed probabilities of every category
  !   within a cluster is equal to the prescribed probability of the cluster.

  ! References:
  !   clubb:ticket:752
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd       ! Constant

    use constants_clubb, only: &
      zero            ! Constant

    use parameters_silhs, only: &
      importance_prob_thresh  ! Variable

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs           ! The real probability for each category

    integer, intent(in) :: &
      num_clusters, &                   ! The number of clusters to sample from
      max_num_categories_in_cluster    ! The max number of categories in each cluster

    integer, dimension(num_clusters), intent(in) :: &
      num_categories_in_cluster        ! The number of categories in each cluster

    integer, dimension(num_clusters,max_num_categories_in_cluster), intent(in) :: &
      cluster_categories            ! An integer matrix containing indices corresponding
                                    ! to the members of the clusters

    real( kind = core_rknd ), dimension(num_clusters), intent(in) :: &
      cluster_prescribed_probs      ! Prescribed probability sum for each cluster

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs     ! Resulting prescribed probability for each individual category

    ! Local Variables
    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_real_probs            ! Total PDF probability for each cluster

    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_prescribed_probs_mod  ! Prescribed probability sum for each cluster, modified to
                                    ! take real probability thresholding into account

    logical, dimension(num_clusters) :: &
      l_cluster_presc_prob_modified ! Whether a given cluster prescribed probability was modified
                                    ! due to thresholding

    real( kind = core_rknd ) :: &
      nonzero_real_clust_sum, &     ! Sum of PDF probabilities for each non-modified cluster
      presc_prob_difference         ! Extra prescribed mass for given cluster to be distributed

    integer :: icluster, jcluster, icategory, cat_idx

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    ! Compute the total PDF probability for each cluster.
    cluster_real_probs(:) = zero
    do icluster=1, num_clusters
      do icategory=1, num_categories_in_cluster(icluster)
        cluster_real_probs(icluster) = cluster_real_probs(icluster) + &
            category_real_probs(cluster_categories(icluster,icategory))
      end do
    end do

    ! Apply thresholding to ensure that clusters with extremely small PDF
    ! probability are not importance sampled.
    do icluster=1, num_clusters
      if ( cluster_real_probs(icluster) < importance_prob_thresh ) then
        ! Thresholding is necessary for this cluster. The prescribed probability for
        ! this cluster will be set equal to the PDF probability of the cluster (that
        ! is, no importance sampling).
        cluster_prescribed_probs_mod(icluster) = cluster_real_probs(icluster)
        l_cluster_presc_prob_modified(icluster) = .true.
      else
        ! Thresholding is not necessary
        cluster_prescribed_probs_mod(icluster) = cluster_prescribed_probs(icluster)
        l_cluster_presc_prob_modified(icluster) = .false.
      end if ! cluster_real_probs(icluster) < importance_prob_thresh .and. ...
    end do ! icluster=1, num_clusters

    ! Distribute any "extra" prescribed probability weight from thresholding to
    ! other clusters.
    if ( any( l_cluster_presc_prob_modified ) ) then

      ! Compute the sum of prescribed probabilities for all non-modified clusters
      nonzero_real_clust_sum = zero
      do icluster=1, num_clusters
        if ( .not. l_cluster_presc_prob_modified(icluster) ) then
          nonzero_real_clust_sum = nonzero_real_clust_sum + cluster_real_probs(icluster)
        end if
      end do ! icluster=1, num_clusters

      ! Transfer extra prescribed probability mass to other clusters.
      do icluster=1, num_clusters
        if ( l_cluster_presc_prob_modified(icluster) ) then

          ! Note that presc_prob_difference may be negative.
          presc_prob_difference = cluster_prescribed_probs(icluster) - &
                                  cluster_prescribed_probs_mod(icluster)

          do jcluster=1, num_clusters
            if ( .not. l_cluster_presc_prob_modified(jcluster) ) then
              cluster_prescribed_probs_mod(jcluster) = cluster_prescribed_probs_mod(jcluster) + &
                ( presc_prob_difference * cluster_real_probs(jcluster) / nonzero_real_clust_sum )
            end if
          end do

        end if ! l_cluster_presc_prob_modified(icluster)
      end do ! icluster=1, num_clusters

    end if ! any( l_cluster_presc_prob_modified )

    ! Finally, compute the prescribed probabilities for each category based on the cluster
    ! probabilities.
    do icluster=1, num_clusters
      do icategory=1, num_categories_in_cluster(icluster)
        cat_idx = cluster_categories(icluster,icategory)
        if ( l_cluster_presc_prob_modified(icluster) ) then
          ! No scaling needs to be done, since the cluster's prescribed probability
          ! equals its PDF probability.
          category_prescribed_probs(cat_idx) = category_real_probs(cat_idx)
        else

          ! Scale category probability based on the cluster probability
          category_prescribed_probs(cat_idx) = ( category_real_probs(cat_idx) / &
            cluster_real_probs(icluster) ) * cluster_prescribed_probs_mod(icluster)

        end if ! l_cluster_presc_prob_modified(icluster)
      end do ! icategory=1, num_categories_in_cluster(icluster)
    end do ! icluster=1, num_clusters

    return
  end function clust_cat_probs_frm_presc_prb
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function cloud_importance_sampling( importance_categories, category_real_probs, &
                                      pdf_params ) &

  result( category_prescribed_probs )

  ! Description:
  !   Applies cloud weighted sampling such that approximately half of all
  !   sample points land in cloud and half land out of cloud!

  ! References:
  !   None :(
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd  ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use constants_clubb, only: &
      one,  &    ! Constant(s)
      two,  &
      zero, &
      fstderr

    use pdf_utilities, only: &
      compute_mean_binormal ! Procedure

    implicit none

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      cloud_frac_min_samp = 0.001_core_rknd, &  ! Minimum cloud fraction for sampling
                                                ! preferentially within cloud
      cloud_frac_max_samp = 0.5_core_rknd       ! Maximum cloud fraction (exclusive) for
                                                ! sampling preferentially within cloud

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories   ! A list of importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs     ! The actual PDF probability for each category

    type(pdf_parameter), intent(in) :: &
      pdf_params

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs ! Probability of each category, scaled such that approximately half
                                ! of all sample points will appear in cloud

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac
    integer :: icategory

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    cloud_frac = compute_mean_binormal( pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, &
                                        pdf_params%mixt_frac )

    if ( cloud_frac >= cloud_frac_min_samp .and. cloud_frac < cloud_frac_max_samp ) then

      ! In-cloud categories ought to be divided by 2*cloud_frac
      ! Out-of-cloud categories should be divided by 2*(1-cloud_frac)

      do icategory=1, num_importance_categories
        if ( importance_categories(icategory)%l_in_cloud ) then
          category_prescribed_probs(icategory) = &
            category_real_probs(icategory) / (two*cloud_frac)
        else
          category_prescribed_probs(icategory) = &
            category_real_probs(icategory) / (two*(one-cloud_frac))
        end if ! importance_categories(icategory)%l_in_cloud
      end do

    else ! cloud_frac < cloud_frac_min_samp .or. cloud_frac >= cloud_frac_max_samp

      ! Do not perform cloud weighted sampling. Let the probabilities remain unmodified.
      category_prescribed_probs = category_real_probs

    end if ! cloud_frac < cloud_frac_min_samp .or. cloud_frac >= cloud_frac_max_samp

    return
  end function cloud_importance_sampling
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine importance_sampling_assertions &
             ( num_samples, importance_categories, category_real_probs, &
               category_prescribed_probs, category_sample_weights, X_u_chi_one_lev, &
               X_u_dp1_one_lev, X_u_dp2_one_lev, lh_sample_point_weights, int_sample_category, &
               pdf_params, hydromet_pdf_params, &
               l_error )

  ! Description:
  !   Various assertion checks for importance sampling are performed here.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd             ! Precision

    use constants_clubb, only: &
      zero, &               ! Constant(s)
      one, &
      fstderr

    use pdf_parameter_module, only: &
      pdf_parameter

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter

    use parameters_silhs, only: &
      l_lh_normalize_weights

    implicit none

    ! Input Variables
    integer :: &
      num_samples                         ! Number of SILHS sample points

    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories               ! The defined importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs,        &       ! The real PDF probabilities for each category
      category_prescribed_probs,  &       ! Prescribed probability for each category
      category_sample_weights             ! Sample weight for each category

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      X_u_chi_one_lev, &                  ! Samples of chi in uniform space
      X_u_dp1_one_lev, &                  ! Samples of the dp1 variate
      X_u_dp2_one_lev, &                  ! Samples of the dp2 variate
      lh_sample_point_weights             ! Weights of samples

    integer, dimension(num_samples), intent(in) :: &
      int_sample_category                 ! An integer for each sample corresponding to the
                                          ! category picked for the sample

    type(pdf_parameter), intent(in) :: &
      pdf_params

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params

    ! Output Variables
    logical, intent(out) :: &
      l_error                             ! True if the assertion check fails.

    ! Local Variables
    real( kind = core_rknd ) :: &
      category_sum, tolerance, weights_avg, num_samples_real

    real( kind = core_rknd ) :: &
      cloud_frac_i, precip_frac_i

    integer :: isample

    type(importance_category_type) :: category

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    !----------------------------------------------------------
    ! Assert that each of the probability groups sum to 1.0
    !----------------------------------------------------------
    tolerance = real( num_importance_categories, kind=core_rknd ) * epsilon( category_sum )

    category_sum = sum( category_real_probs )
    if ( abs( category_sum - one ) > tolerance ) then
      write(fstderr,*) "The real category probabilities do not sum to one."
      write(fstderr,*) "sum( category_real_probs ) = ", category_sum
      l_error = .true.
    end if

    category_sum = sum( category_prescribed_probs )
    if ( abs( category_sum - one ) > tolerance ) then
      write(fstderr,*) "The prescribed category probabilities do not sum to one."
      write(fstderr,*) "sum( category_prescribed_probs ) = ", category_sum
      write(fstderr,*) 'category_real_probs = ', category_real_probs
      write(fstderr,*) 'category_prescribed_probs = ', category_prescribed_probs
      l_error = .true.
    end if

    category_sum = sum( category_sample_weights * category_prescribed_probs )
    if ( abs( category_sum - one ) > tolerance ) then
      write(fstderr,*) "The weighted prescribed category probabilities do not sum to one."
      write(fstderr,*) "sum( category_sample_weights * category_prescribed_probs ) = ", category_sum
      l_error = .true.
    end if

    !------------------------------------------------------------------
    ! Assert that sample point weights average to 1.0 (if enabled)
    !------------------------------------------------------------------
    if ( l_lh_normalize_weights ) then
      num_samples_real = real( num_samples, kind=core_rknd )
      weights_avg = sum( lh_sample_point_weights(:) ) / num_samples_real
      if ( abs( weights_avg - one ) > num_samples_real*epsilon( weights_avg ) ) then
        write(fstderr,*) "The sample point weights do not average to one."
        write(fstderr,*) "num_samples = ", num_samples
        write(fstderr,*) "sum( lh_sample_point_weights(:) ) = ", sum( lh_sample_point_weights(:) )
        write(fstderr,*) "avg( lh_sample_point_weights(:) ) = ", weights_avg
        l_error = .true.
      end if
    end if

    !---------------------------------------------------------------------
    ! Verify that samples have been correctly scaled to reside in the
    ! appropriate category
    !---------------------------------------------------------------------
    do isample = 1, num_samples

      category = importance_categories(int_sample_category(isample))

      ! Verification of component
      if ( category%l_in_component_1 ) then
        if ( X_u_dp1_one_lev(isample) > pdf_params%mixt_frac ) then
          write(fstderr,*) "The component of a sample is incorrect."
          l_error = .true.
        end if
        cloud_frac_i = pdf_params%cloud_frac_1
        precip_frac_i = hydromet_pdf_params%precip_frac_1
      else ! .not. category%l_in_component_1
        if ( X_u_dp1_one_lev(isample) < pdf_params%mixt_frac ) then
          write(fstderr,*) "The component of a sample is incorrect."
          l_error = .true.
        end if
        cloud_frac_i = pdf_params%cloud_frac_2
        precip_frac_i = hydromet_pdf_params%precip_frac_2
      end if ! category%l_in_component_1

      ! Verification of cloud
      if ( category%l_in_cloud ) then
        if ( X_u_chi_one_lev(isample) < (one - cloud_frac_i) ) then
          write(fstderr,*) "The chi element of a sample is incorrect."
          l_error = .true.
        end if
      else
        if ( X_u_chi_one_lev(isample) > (one - cloud_frac_i) ) then
          write(fstderr,*) "The chi element of a sample is incorrect."
          l_error = .true.
        end if
      end if

      ! Verification of precipitation
      if ( category%l_in_precip ) then
        if ( X_u_dp2_one_lev(isample) > precip_frac_i ) then
          write(fstderr,*) "The in-precipitation status of a sample is incorrect."
          l_error = .true.
        end if
      else
        if ( X_u_dp2_one_lev(isample) < precip_frac_i ) then
          write(fstderr,*) "The in-precipitation status of a sample is incorrect."
          l_error = .true.
        end if
      end if

    end do ! isample = 1, num_samples

    return
  end subroutine importance_sampling_assertions
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine cloud_weighted_sampling_driver &
             ( num_samples, p_matrix_chi, p_matrix_dp1, &
               cloud_frac_1, cloud_frac_2, mixt_frac, &
               X_u_chi, X_u_dp1, &
               lh_sample_point_weights )

  ! Description:
  !   Performs importance sampling such that half of sample points are in cloud

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd          ! Precision(s)

    use generate_uniform_sample_module, only: &
      choose_permuted_random ! Procedure

    use pdf_utilities, only: &
      compute_mean_binormal

    use constants_clubb, only: &
      one, &
      two

    implicit none

    ! Parameter Constants
    real( kind = core_rknd ), parameter :: &
      cloud_frac_min_samp = 0.001_core_rknd, &  ! Minimum cloud fraction for sampling
                                                ! preferentially within cloud
      cloud_frac_max_samp = 0.5_core_rknd       ! Maximum cloud fraction (exclusive) for
                                                ! sampling preferentiallywithin cloud

    ! Input Variables
    integer :: &
      num_samples                         ! Number of SILHS sample points

    integer, dimension(num_samples), intent(in) :: &
      p_matrix_chi, &                     ! Permutation of the integers 0..num_samples
                                          ! from p_matrix for chi
      p_matrix_dp1                        ! Elements from p_matrix for dp1 element

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac_1, &                     ! Cloud fraction in PDF component 1
      cloud_frac_2, &                     ! Cloud fraction in PDF component 2
      mixt_frac                           ! Weight of first gaussian component

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(num_samples), intent(inout) :: &
      X_u_chi, &                          ! Samples of chi in uniform space
      X_u_dp1                             ! Samples of the dp1 variate for determining mixture
                                          ! component

    ! Output Variables
    real( kind = core_rknd ), dimension(num_samples), intent(out) :: &
      lh_sample_point_weights             ! Weight of SILHS sample points (these must be applied
                                          ! when averaging results from, e.g., calling microphysics

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac                          ! Cloud fraction at k_lh_start

    real( kind = core_rknd ), dimension(num_samples/2,1) :: &
      mixt_rand_cloud, &                  ! Stratified random variable for determining mixture
                                          ! component in cloud
      mixt_rand_clear                     ! Stratified random variable for determining mixture
                                          ! component in clear air

    real( kind = core_rknd ) :: &
      mixt_rand_element

    real( kind = core_rknd ) :: &
      lh_sample_cloud_weight,     &       ! Weight of a cloudy sample point
      lh_sample_clear_weight              ! Weight of a clear  sample point

    integer, dimension(num_samples/2,1) :: &
      mixt_permuted_cloud, &              ! Permuted random numbers for mixt_rand_cloud
      mixt_permuted_clear                 ! Permuted random numbers for mixt_rand_clear

    integer :: sample, n_cloudy_samples, n_clear_samples
    logical :: l_cloudy_sample

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    cloud_frac = compute_mean_binormal( cloud_frac_1, cloud_frac_2, mixt_frac )

    if ( cloud_frac >= cloud_frac_min_samp .and. cloud_frac < cloud_frac_max_samp ) then

      n_cloudy_samples = 0
      n_clear_samples  = 0
      ! Pick two stratified random numbers for determining a mixture fraction.
      do sample=1, num_samples
        ! We arbitrarily choose that p_matrix elements less than num_samples/2
        ! will be used for mixt_rand_cloud
        if ( p_matrix_dp1(sample) < ( num_samples / 2 ) ) then
          n_cloudy_samples = n_cloudy_samples + 1
          mixt_permuted_cloud(n_cloudy_samples,1) = p_matrix_dp1(sample)
        else if ( p_matrix_dp1(sample) >= ( num_samples / 2 ) ) then
          n_clear_samples = n_clear_samples + 1
          mixt_permuted_clear(n_clear_samples,1) = p_matrix_dp1(sample) - (num_samples / 2)
        end if ! p_matrix_dp1(sample) < ( num_samples / 2 )
      end do

      ! Generate the stratified uniform numbers!
      do sample = 1, num_samples/2
        mixt_rand_cloud(sample,1) = &
          choose_permuted_random( num_samples/2, mixt_permuted_cloud(sample,1) )
        mixt_rand_clear(sample,1) = &
          choose_permuted_random( num_samples/2, mixt_permuted_clear(sample,1) )
      end do

      n_cloudy_samples = 0
      n_clear_samples  = 0

      ! Compute weights for each type of point
      lh_sample_cloud_weight = two * cloud_frac
      lh_sample_clear_weight = two - lh_sample_cloud_weight

      do sample=1, num_samples

        ! Detect which half of the sample points are in clear air and which half are in
        ! the cloudy air
        if ( p_matrix_chi(sample) < ( num_samples / 2 ) ) then

          l_cloudy_sample = .false.
          lh_sample_point_weights(sample) = lh_sample_clear_weight
          n_clear_samples = n_clear_samples + 1
          mixt_rand_element = mixt_rand_clear(n_clear_samples,1)
        else

          l_cloudy_sample = .true.
          lh_sample_point_weights(sample) = lh_sample_cloud_weight
          n_cloudy_samples = n_cloudy_samples + 1
          mixt_rand_element = mixt_rand_cloud(n_cloudy_samples,1)
        end if

        ! Transpose and scale the points to be in or out of cloud
        call choose_X_u_scaled &
             ( l_cloudy_sample, & ! In
               p_matrix_chi(sample), num_samples, & ! In
               cloud_frac_1, cloud_frac_2, & ! In
               mixt_frac, mixt_rand_element, & !In
               X_u_dp1(sample), X_u_chi(sample) ) ! Out

      end do ! sample=1, num_samples

    else ! cloud_frac < cloud_frac_min_samp .or. cloud_frac >= cloud_frac_max_samp

      ! Do not perform cloud weighted sampling.
      lh_sample_point_weights(:) = one

    end if ! cloud_frac >= cloud_frac_min_samp .and. cloud_frac < cloud_frac_max_samp

    return
  end subroutine cloud_weighted_sampling_driver
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  function generate_strat_uniform_variate( num_samples ) &

  result( stratified_variate )

  ! Description:
  !   Generates a stratified uniform sample for a single variable

  ! References:
  !   None
  !-----------------------------------

    use clubb_precision, only: &
      core_rknd    ! Constant

    use generate_uniform_sample_module, only: &
      rand_permute ! Procedure

    use generate_uniform_sample_module, only: &
      choose_permuted_random

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      num_samples ! Number of SILHS sample points

    ! Output Variable
    real( kind = core_rknd ), dimension(num_samples) :: &
      stratified_variate      ! Uniform samples stratified in (0,1)

    ! Local Variables

    ! Vector of the integers 0,1,2,...,num_samples in random order
    integer, dimension(num_samples) :: pvect

    integer :: sample

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    !-------------------------------------------------------------------
    ! Draw the integers 0,1,2,...,num_samples in random order
    !-------------------------------------------------------------------
    call rand_permute( num_samples, pvect )

    !----------------------------------------------------------------------
    ! For each permuted integer (each box), determine a random real number
    !----------------------------------------------------------------------
    do sample=1, num_samples
      stratified_variate(sample) = &
        real( choose_permuted_random( num_samples, pvect(sample) ), kind = core_rknd )
    end do

    return
  end function generate_strat_uniform_variate
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine choose_X_u_scaled &
             ( l_cloudy_sample, &
               p_matrix_element, num_samples, &
               cloud_frac_1, cloud_frac_2, &
               mixt_frac, mixt_rand_element, &
               X_u_dp1_element, X_u_chi_element )

! Description:
!   Find a clear or cloudy point for sampling.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use generate_uniform_sample_module, only: &
      choose_permuted_random    ! Procedure

    use clubb_precision, only: &
      core_rknd ! Konstant

    use constants_clubb, only: &
      one

    implicit none

    ! Input Variables
    logical, intent(in) :: &
      l_cloudy_sample ! Whether his is a cloudy or clear air sample point

    integer, intent(in) :: &
      p_matrix_element, & ! Integer from 0..num_samples for this sample
      num_samples       ! Total number of calls to the microphysics

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac_1, &    ! Cloud fraction associated with mixture component 1     [-]
      cloud_frac_2, &    ! Cloud fraction associated with mixture component 2     [-]
      mixt_frac         ! Mixture fraction                                       [-]

    real( kind = core_rknd ), intent(in) :: &
      mixt_rand_element ! Random number (0,1) for determining mixture component

    ! Output Variables
    real(kind=core_rknd), intent(out) :: &
      X_u_dp1_element, X_u_chi_element ! Elements from X_u (uniform dist.)

    ! Local Variables
    real(kind=core_rknd) :: cloud_frac_i, conditional_mixt_frac

    real( kind = core_rknd ) :: &
      cld_comp1_frac,    &          ! Fraction of points in component 1 and cloud
      cld_comp2_frac,    &          ! Fraction of points in component 2 and cloud
      nocld_comp1_frac,  &          ! Fraction of points in component 1 and clear air
      nocld_comp2_frac              ! Fraction of points in component 2 and clear air

    integer :: X_mixt_comp_one_lev, p_matrix_element_ranged

    real( kind = core_rknd ) :: mixt_rand_element_scaled, chi_rand_element

    ! ---- Begin code ----

    cld_comp1_frac = mixt_frac*cloud_frac_1
    cld_comp2_frac = (one-mixt_frac)*cloud_frac_2

    !---------------------------------------------------------------------
    ! Determine the conditional mixture fraction, given whether we are in
    ! cloud
    !---------------------------------------------------------------------
    if ( l_cloudy_sample ) then

      conditional_mixt_frac = cld_comp1_frac / ( cld_comp1_frac + cld_comp2_frac )

    else ! .not. l_cloudy_sample

      nocld_comp1_frac = mixt_frac - cld_comp1_frac
      nocld_comp2_frac = (one - mixt_frac) - cld_comp2_frac

      conditional_mixt_frac = nocld_comp1_frac / (nocld_comp1_frac + nocld_comp2_frac)

    end if ! l_cloudy_sample

    !---------------------------------------------------------------------
    ! Determine mixture component given the conditional mixture fraction
    !---------------------------------------------------------------------
!    if ( in_mixt_comp_1( mixt_rand_element, conditional_mixt_frac ) ) then
    if ( mixt_rand_element < conditional_mixt_frac ) then
      X_mixt_comp_one_lev = 1
    else
      X_mixt_comp_one_lev = 2
    end if

    !---------------------------------------------------------------------
    ! Determine dp1 element given mixture component
    !---------------------------------------------------------------------
    if ( X_mixt_comp_one_lev == 1 ) then
      ! mixt_rand_element is scaled to give real number stratified in (0,1)
      mixt_rand_element_scaled = mixt_rand_element / conditional_mixt_frac
      X_u_dp1_element = mixt_rand_element_scaled * mixt_frac

    else if ( X_mixt_comp_one_lev == 2 ) then
      ! mixt_rand_element is scaled to give real number stratified in (0,1)
      mixt_rand_element_scaled = (mixt_rand_element - conditional_mixt_frac) / &
                                   (one - conditional_mixt_frac) 
      X_u_dp1_element = mixt_rand_element_scaled * (one - mixt_frac) + mixt_frac

    else
      stop "Should not be here"
    end if ! X_mixt_comp_one_lev == 1

    !---------------------------------------------------------------------
    ! Determine chi element given mixture component and l_cloudy_sample
    !---------------------------------------------------------------------
    ! Get p_matrix_element in the proper range
    if ( p_matrix_element >= ( num_samples / 2 ) ) then
      p_matrix_element_ranged = p_matrix_element - ( num_samples / 2 )
    else
      p_matrix_element_ranged = p_matrix_element
    end if
    
    ! Get a stratified number in (0,1) using the element of p_matrix!
    chi_rand_element = choose_permuted_random( num_samples/2, p_matrix_element_ranged )

    ! Determine cloud fraction
    if ( X_mixt_comp_one_lev == 1 ) then
      cloud_frac_i = cloud_frac_1
    else
      cloud_frac_i = cloud_frac_2
    end if

    if ( l_cloudy_sample ) then
      ! Scale and translate sample point to reside in cloud
      X_u_chi_element = cloud_frac_i * chi_rand_element + &
                          (one - cloud_frac_i )
    else
      ! Scale and translate sample point to reside in clear air (no cloud)
      X_u_chi_element = (one-cloud_frac_i) * chi_rand_element
    end if ! l_cloudy_sample

    return
  end subroutine choose_X_u_scaled
!----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function determine_sample_categories( num_samples, d_variables, X_nl_one_lev, &
                                        X_mixt_comp_one_lev, importance_categories ) &
  result( int_sample_category )

  ! Description:
  !   Determines the importance category of each sample.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd       ! Constant

    use constants_clubb, only: &
      zero     ! Constant(s)

    use corr_varnce_module, only: &
      iiPDF_chi, &
      iiPDF_rr

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples, &        ! Number of SILHS sample points
      d_variables           ! Number of variates in X_nl

    real( kind = core_rknd ), dimension(num_samples,d_variables), intent(in) :: &
      X_nl_one_lev          ! SILHS sample vector at one height level

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp_one_lev   ! Category of each sample point

    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories

    ! Output Variable
    integer, dimension(num_samples) :: &
      int_sample_category   ! Category of each sample

    type(importance_category_type) :: sample_category

    integer :: isample, icategory, found_category_index

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do isample=1, num_samples

      if ( X_nl_one_lev(isample,iiPDF_chi) < zero ) then
        sample_category%l_in_cloud = .false.
      else
        sample_category%l_in_cloud = .true.
      end if

      if ( iiPDF_rr == -1 ) then
        stop "iiPDF_rr must be greater than zero for the category sampler to work."
      end if

      if ( X_nl_one_lev(isample,iiPDF_rr) > zero ) then
        sample_category%l_in_precip = .true.
      else
        sample_category%l_in_precip = .false.
      end if

      if ( X_mixt_comp_one_lev(isample) == 1 ) then
        sample_category%l_in_component_1 = .true.
      else
        sample_category%l_in_component_1 = .false.
      end if

      found_category_index = -1

      do icategory=1, num_importance_categories

        if ( (importance_categories(icategory)%l_in_cloud .eqv. sample_category%l_in_cloud) .and. &
             (importance_categories(icategory)%l_in_precip .eqv. sample_category%l_in_precip) .and.&
             (importance_categories(icategory)%l_in_component_1 .eqv. &
              sample_category%l_in_component_1) ) then

          found_category_index = icategory
          exit

        end if

      end do ! icategory=1, num_importance_categories

      if ( found_category_index == -1 ) then
        stop "Fatal error determining category in determine_sample_categories"
      end if

      int_sample_category(isample) = found_category_index

    end do ! isample=1, num_samples

    return
  end function determine_sample_categories
  !-----------------------------------------------------------------------

end module silhs_importance_sample_module
