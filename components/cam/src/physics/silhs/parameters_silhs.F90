!-----------------------------------------------------------------------
! $Id: parameters_silhs.F90 7864 2015-07-25 20:40:58Z raut@uwm.edu $
!===============================================================================
module parameters_silhs

! Description:
!   Parameters for SILHS!

! References:
!   None
!-------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd     ! Constant

  implicit none

  ! Cluster allocation strategies!!!
  integer, parameter, public :: &
    ! All eight categories, effectively no clustering
    eight_cluster_allocation_opt = 1, &
    ! Four clusters for the combinations of cloud/no cloud and component 1/2.
    ! Precipitation fraction is ignored.
    four_cluster_allocation_opt  = 2, &
    ! Two clusters, one containing all categories with either cloud or precip,
    ! and the other containing categories with neither
    two_cluster_cp_nocp_opt      = 3

  integer, public :: &
    cluster_allocation_strategy = eight_cluster_allocation_opt

  !$omp threadprivate( cluster_allocation_strategy )

  ! The following type defines parameters that control the sample point
  ! allocation for the clustered sampling scheme
  ! (l_lh_clustered_sampling = .true.).
  type eight_cluster_presc_probs_type

    real( kind = core_rknd ) :: &
      cloud_precip_comp1      = 0.15_core_rknd, &
      cloud_precip_comp2      = 0.15_core_rknd, &
      nocloud_precip_comp1    = 0.15_core_rknd, &
      nocloud_precip_comp2    = 0.15_core_rknd, &
      cloud_noprecip_comp1    = 0.15_core_rknd, &
      cloud_noprecip_comp2    = 0.15_core_rknd, &
      nocloud_noprecip_comp1  = 0.05_core_rknd, &
      nocloud_noprecip_comp2  = 0.05_core_rknd

  end type eight_cluster_presc_probs_type

  ! Flags for the SILHS sampling code 
  logical, public :: &
    l_lh_importance_sampling  = .true., &     ! Limit noise by performing importance sampling
    l_Lscale_vert_avg         = .true., &     ! Calculate Lscale_vert_avg in lh_subcolumn_generator
    l_lh_straight_mc          = .false.,&     ! Use true Monte Carlo sampling with no Latin
                                              !  hypercube sampling and no importance sampling
    l_lh_clustered_sampling   = .true., &     ! Use the "new" SILHS importance sampling
                                              !  scheme with prescribed probabilities
    l_rcm_in_cloud_k_lh_start = .false.,&     ! Determine k_lh_start based on maximum within-cloud
                                              !  rcm
    l_random_k_lh_start       = .false.       ! Place k_lh_start at a random grid level between
                                              !  maximum rcm and maximum rcm_in_cloud

  !$omp threadprivate( l_lh_importance_sampling, l_Lscale_vert_avg, l_lh_straight_mc, &
  !$omp                l_lh_clustered_sampling, l_rcm_in_cloud_k_lh_start, l_random_k_lh_start )

  type(eight_cluster_presc_probs_type), public, save :: &
    eight_cluster_presc_probs                 ! Prescribed probabilities for
                                              ! l_lh_clustered_sampling = .true.

  !$omp threadprivate( eight_cluster_presc_probs )

  logical, public :: &
    l_lh_limit_weights = .false., &           ! Limit SILHS sample point weights for stability
    l_lh_var_frac      = .false.              ! Prescribe variance fractions

  !$omp threadprivate( l_lh_limit_weights, l_lh_var_frac )

  real( kind = core_rknd ), public :: &
    importance_prob_thresh = 5.0e-3_core_rknd ! Minimum PDF probability of category for importance
                                              ! sampling

  !$omp threadprivate( importance_prob_thresh )

  private ! Default Scope

  public :: eight_cluster_presc_probs_type

end module parameters_silhs
