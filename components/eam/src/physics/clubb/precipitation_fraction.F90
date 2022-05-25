!-------------------------------------------------------------------------
! $Id$
!===============================================================================
module precipitation_fraction

  ! Description:
  ! Sets overall precipitation fraction as well as the precipitation fraction
  ! in each PDF component.

  implicit none

  private

  public :: precip_fraction

  private :: component_precip_frac_weighted, &
             component_precip_frac_specify,  &
             precip_frac_assert_check

  integer, parameter, public :: &
    precip_frac_calc_type = 2  ! Option used to calculate component precip_frac

  contains

  !=============================================================================
  subroutine precip_fraction( nz, ngrdcol, &
                              hydromet, cloud_frac, cloud_frac_1, &
                              cloud_frac_2, ice_supersat_frac, &
                              ice_supersat_frac_1, ice_supersat_frac_2, &
                              mixt_frac, clubb_params, l_stats_samp, &
                              stats_sfc, & 
                              precip_frac, &
                              precip_frac_1, &
                              precip_frac_2, &
                              precip_frac_tol )

    ! Description:
    ! Determines (overall) precipitation fraction over the horizontal domain, as
    ! well as the precipitation fraction within each PDF component, at every
    ! vertical grid level.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        cloud_frac_min, &
        fstderr

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        iupsilon_precip_frac_rat

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        l_mix_rat_hm, & ! Variable(s)
        l_frozen_hm,  &
        hydromet_tol

    use stats_variables, only: &
        iprecip_frac_tol

    use stats_type_utilities, only: &
        stat_update_var_pt  ! Procedure(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level, &   ! Procedure
        err_code, &                     ! Error Indicator
        clubb_fatal_error               ! Constant

    use stats_type, only: stats ! Type

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,       & ! Number of model vertical grid levels
      ngrdcol     ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      cloud_frac,          & ! Cloud fraction (overall)                      [-]
      cloud_frac_1,        & ! Cloud fraction (1st PDF component)            [-]
      cloud_frac_2,        & ! Cloud fraction (2nd PDF component)            [-]
      ice_supersat_frac,   & ! Ice supersaturation fraction (overall)        [-]
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2, & ! Ice supersaturation fraction (2nd PDF comp.)  [-]
      mixt_frac              ! Mixture fraction                              [-]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.
      
    ! Inout Variables
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_sfc

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)               [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    ! Local Variables

    ! "Maximum allowable" hydrometeor mixing ratio in-precip component mean.
    real( kind = core_rknd ), parameter :: &
      max_hm_ip_comp_mean = 0.0025_core_rknd  ! [kg/kg]

    real( kind = core_rknd ), parameter :: &
      precip_frac_tol_coef = 0.1_core_rknd  ! Coefficient for precip_frac_tol

    integer :: &
      j, k, ivar   ! Loop indices


    ! Initialize the precipitation fraction variables (precip_frac,
    ! precip_frac_1, and precip_frac_2) to 0.
    precip_frac(:,:)   = zero
    precip_frac_1(:,:) = zero
    precip_frac_2(:,:) = zero

    ! Set the minimum allowable precipitation fraction when hydrometeors are
    ! found at a grid level.
    if ( any( l_frozen_hm ) ) then
      ! Ice microphysics included.
      do j = 1, ngrdcol
        precip_frac_tol(j) &
        = max( precip_frac_tol_coef &
               * max( maxval( cloud_frac(j,:) ), maxval( ice_supersat_frac(j,:) ) ), &
               cloud_frac_min )
      end do
    else
      ! Warm microphysics.
      do j = 1, ngrdcol
        precip_frac_tol(j) = max( precip_frac_tol_coef * maxval( cloud_frac(j,:) ), &
                                cloud_frac_min )
      end do
    endif
    
    !!! Find overall precipitation fraction.
    ! The precipitation fraction is the greatest cloud fraction at or above a
    ! vertical level.
    if ( any( l_frozen_hm ) ) then
      
      ! Ice microphysics included.
      precip_frac(:,nz) = max( cloud_frac(:,nz), ice_supersat_frac(:,nz) )
      
      do k = nz-1, 1, -1
        precip_frac(:,k) = max( precip_frac(:,k+1), cloud_frac(:,k), &
                                ice_supersat_frac(:,k) )
      end do
      
    else
      
      ! Warm microphysics.
      precip_frac(:,nz) = cloud_frac(:,nz)
      
      do k = nz-1, 1, -1
        precip_frac(:,k) = max( precip_frac(:,k+1), cloud_frac(:,k) )
      end do
      
    end if

    !!! Special checks for overall precipitation fraction
    do k = 1, nz, 1
      do j = 1, ngrdcol

        if ( any( hydromet(j,k,:) >= hydromet_tol(:) ) &
            .and. precip_frac(j,k) < precip_frac_tol(j) ) then

          ! In a scenario where we find any hydrometeor at this grid level, but
          ! no cloud at or above this grid level, set precipitation fraction to
          ! a minimum threshold value.
          precip_frac(j,k) = precip_frac_tol(j)

        elseif ( all( hydromet(j,k,:) < hydromet_tol(:) ) ) then

          ! The means (overall) of every precipitating hydrometeor are all less
          ! than their respective tolerance amounts.  They are all considered to
          ! have values of 0.  There are not any hydrometeor species found at
          ! this grid level.  There is also no cloud at or above this grid
          ! level, so set precipitation fraction to 0.
          precip_frac(j,k) = zero

        endif
        
      end do
    enddo ! Special checks for overall precipitation fraction loop: k = 1, nz, 1


    !!! Find precipitation fraction within each PDF component.
    !
    ! The overall precipitation fraction, f_p, is given by the equation:
    !
    ! f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
    !
    ! where "a" is the mixture fraction (weight of PDF component 1), f_p(1) is
    ! the precipitation fraction within PDF component 1, and f_p(2) is the
    ! precipitation fraction within PDF component 2.  Overall precipitation
    ! fraction is found according the method above, and mixture fraction is
    ! already determined, leaving f_p(1) and f_p(2) to be solved for.  The
    ! values for f_p(1) and f_p(2) must satisfy the above equation.
    if ( precip_frac_calc_type == 1 ) then

      ! Calculatate precip_frac_1 and precip_frac_2 based on the greatest
      ! weighted cloud_frac_1 at or above a grid level.
      call component_precip_frac_weighted( nz, ngrdcol, & ! intent(in)
                                           hydromet(:,:,:), precip_frac(:,:), & ! intent(in)
                                           cloud_frac_1(:,:), cloud_frac_2(:,:), & ! intent(in)
                                           ice_supersat_frac_1(:,:), & ! intent(in)
                                           ice_supersat_frac_2(:,:), mixt_frac(:,:), & !intent(in)
                                           precip_frac_tol(:), & ! intent(in)
                                           precip_frac_1(:,:), precip_frac_2(:,:) ) ! intent(out)
                                            
    elseif ( precip_frac_calc_type == 2 ) then

      ! Specified method.
      call component_precip_frac_specify( nz, ngrdcol, & ! intent(in)
                                          clubb_params(iupsilon_precip_frac_rat), & ! intent(in)
                                          hydromet(:,:,:), precip_frac(:,:), & ! intent(in)
                                          mixt_frac(:,:), precip_frac_tol(:), & ! intent(in)
                                          precip_frac_1(:,:), precip_frac_2(:,:) ) ! intent(out)

    else ! Invalid option selected.

       write(fstderr,*) "Invalid option to calculate precip_frac_1 " &
                        // "and precip_frac_2."
       err_code = clubb_fatal_error
       return

    endif ! precip_frac_calc_type


    ! Increase Precipiation Fraction under special conditions.
    !
    ! There are scenarios that sometimes occur that require precipitation
    ! fraction to be boosted.  Precipitation fraction is calculated from cloud
    ! fraction and ice supersaturation fraction.  For numerical reasons, CLUBB's
    ! PDF may become entirely subsaturated with respect to liquid and ice,
    ! resulting in both a cloud fraction of 0 and an ice supersaturation
    ! fraction of 0.  When this happens, precipitation fraction drops to 0 when
    ! there aren't any hydrometeors present at that grid level, or to
    ! precip_frac_tol when there is at least one hydrometeor present at that
    ! grid level.  However, sometimes there are large values of hydrometeors
    ! found at that grid level.  When this occurs, the PDF component in-precip
    ! mean of a hydrometeor can become ridiculously large.  This is because the
    ! ith PDF component in-precip mean of a hydrometeor, mu_hm_i,  is given by
    ! the equation:
    !
    ! mu_hm_i = hm_i / precip_frac_i;
    !
    ! where hm_i is the overall ith PDF component mean of the hydrometeor, and
    ! precip_frac_i is the ith PDF component precipitation fraction.  When
    ! precip_frac_i has a value of precip_frac_tol and hm_i is large, mu_hm_i
    ! can be huge.  This can cause enormous microphysical process rates and
    ! result in numerical instability.  It is also very inaccurate.
    !
    ! In order to limit this problem, the ith PDF component precipitation
    ! fraction is increased in order to decrease mu_hm_i.  First, an "upper
    ! limit" is set for mu_hm_i when the hydrometeor is a mixing ratio.  This is
    ! called max_hm_ip_comp_mean.  At every vertical level and for every
    ! hydrometeor mixing ratio, a check is made to try to prevent mu_hm_i from
    ! exceeding the "upper limit".  The check is:
    !
    ! hm_i / precip_frac_i ( which = mu_hm_i )  >  max_hm_ip_comp_mean,
    !
    ! which can be rewritten:
    !
    ! hm_i > precip_frac_i * max_hm_ip_comp_mean.
    !
    ! Since hm_i has not been calculated yet, the assumption for this check is
    ! that all of the hydrometeor is found in one PDF component, which is the
    ! worst-case scenario in violating this limit.  The check becomes:
    !
    ! <hm> / ( mixt_frac * precip_frac_1 )  >  max_hm_ip_comp_mean;
    !    in PDF comp. 1; and
    ! <hm> / ( ( 1 - mixt_frac ) * precip_frac_2 )  >  max_hm_ip_comp_mean;
    !    in PDF comp. 2.
    !
    ! These limits can be rewritten as:
    !
    ! <hm>  >  mixt_frac * precip_frac_1 * max_hm_ip_comp_mean;
    !    in PDF comp. 1; and
    ! <hm>  >  ( 1 - mixt_frac ) * precip_frac_2 * max_hm_ip_comp_mean;
    !    in PDF comp. 2.
    !
    ! When component precipitation fraction is found to be in excess of the
    ! limit, precip_frac_i is increased to:
    !
    ! <hm> / ( mixt_frac * max_hm_ip_comp_mean );
    !    when the limit is exceeded in PDF comp. 1; and
    ! <hm> / ( ( 1 - mixt_frac ) * max_hm_ip_comp_mean );
    !    when the limit is exceeded in PDF comp. 2.
    !
    ! Of course, precip_frac_i is not allowed to exceed 1, so when
    ! <hm> / mixt_frac (or <hm> / ( 1 - mixt_frac )) is already greater than
    ! max_hm_ip_comp_mean, mu_hm_i will also have to be greater than
    ! max_hm_ip_comp_mean.  However, the value of mu_hm_i is still reduced when
    ! compared to what it would have been using precip_frac_tol.  In the event
    ! that multiple hydrometeor mixing ratios violate the check, the code is set
    ! up so that precip_frac_i is increased based on the highest hm_i.
    

   do ivar = 1, hydromet_dim
      if ( l_mix_rat_hm(ivar) ) then
        do k = 1, nz
          do j = 1, ngrdcol

            ! The hydrometeor is a mixing ratio.

            if ( hydromet(j,k,ivar) >= hydromet_tol(ivar) .and. &
                 hydromet(j,k,ivar) > mixt_frac(j,k) * precip_frac_1(j,k) &
                                     * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 1st PDF component.
                precip_frac_1(j,k) &
                = min( hydromet(j,k,ivar) &
                       / ( mixt_frac(j,k) * max_hm_ip_comp_mean ), one )

                ! The value of precip_frac_1 must be at least precip_frac_tol
                ! when precipitation is found in the 1st PDF component.
                precip_frac_1(j,k) = max( precip_frac_1(j,k), precip_frac_tol(j) )

            endif ! <hm>/(mixt_frac*precip_frac_1) > max_hm_ip_comp_mean

            if ( hydromet(j,k,ivar) >= hydromet_tol(ivar) .and. &
                  hydromet(j,k,ivar) > ( one - mixt_frac(j,k) ) * precip_frac_2(j,k) &
                                     * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 2nd PDF component.
                precip_frac_2(j,k) &
                = min( hydromet(j,k,ivar) &
                       / ( ( one - mixt_frac(j,k) ) * max_hm_ip_comp_mean ), one )

                ! The value of precip_frac_2 must be at least precip_frac_tol
                ! when precipitation is found in the 2nd PDF component.
                precip_frac_2(j,k) = max( precip_frac_2(j,k), precip_frac_tol(j) )

            endif ! <hm>/((1-mixt_frac)*precip_frac_2) > max_hm_ip_comp_mean

          end do
        end do ! k = 1, nz
      end if ! l_mix_rat_hm(ivar)
    end do ! ivar = 1, hydromet_dim

    ! Recalculate overall precipitation fraction for consistency.
    precip_frac(:,:) = mixt_frac(:,:) * precip_frac_1(:,:) &
                       + ( one - mixt_frac(:,:) ) * precip_frac_2(:,:)

    ! Double check that precip_frac_tol <= precip_frac <= 1 when hydrometeors
    ! are found at a grid level.
    ! PLEASE DO NOT ALTER precip_frac, precip_frac_1, or precip_frac_2 anymore
    ! after this point in the code.
    do k = 1, nz, 1
      do j = 1, ngrdcol
        if ( any( hydromet(j,k,:) >= hydromet_tol(:) ) ) then
          precip_frac(j,k) = min( max( precip_frac(j,k), precip_frac_tol(j) ), one )
        end if ! any( hydromet(k,:) >= hydromet_tol(:) )
      end do
    end do ! k = 1, nz, 1


    ! Statistics
    if ( l_stats_samp ) then
      if ( iprecip_frac_tol > 0 ) then
        do j = 1, ngrdcol
          call stat_update_var_pt( iprecip_frac_tol, 1, precip_frac_tol(j), & ! intent(in)
                                   stats_sfc(j) ) ! intent(inout)
        end do
      end if ! iprecip_frac_tol
    end if ! l_stats_samp


    ! Assertion check for precip_frac, precip_frac_1, and precip_frac_2.
    if ( clubb_at_least_debug_level( 2 ) ) then
      do j = 1, ngrdcol
        call precip_frac_assert_check( nz, hydromet(j,:,:), mixt_frac(j,:), precip_frac(j,:), & !in
                                       precip_frac_1(j,:), precip_frac_2(j,:), & ! intent(in)
                                       precip_frac_tol(j) ) ! intent(in)
      end do
    endif

    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine component_precip_frac_weighted( nz, ngrdcol, &
                                             hydromet, precip_frac, &
                                             cloud_frac_1, cloud_frac_2, &
                                             ice_supersat_frac_1, &
                                             ice_supersat_frac_2, mixt_frac, &
                                             precip_frac_tol, &
                                             precip_frac_1, precip_frac_2 )

    ! Description:
    ! Set precipitation fraction in each component of the PDF.  The weighted 1st
    ! PDF component precipitation fraction (weighted_pfrac_1) at a grid level is
    ! calculated by the greatest value of mixt_frac * cloud_frac_1 at or above
    ! the relevant grid level.  Likewise, the weighted 2nd PDF component
    ! precipitation fraction (weighted_pfrac_2) at a grid level is calculated by
    ! the greatest value of ( 1 - mixt_frac ) * cloud_frac_2 at or above the
    ! relevant grid level.
    !
    ! The fraction weighted_pfrac_1 / ( weighted_pfrac_1 + weighted_pfrac_2 ) is
    ! the weighted_pfrac_1 fraction.  Multiplying this fraction by overall
    ! precipitation fraction and then dividing by mixt_frac produces the 1st PDF
    ! component precipitation fraction (precip_frac_1).  Then, calculate the 2nd
    ! PDF component precipitation fraction (precip_frac_2) accordingly.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        l_frozen_hm,  & ! Variable(s)
        hydromet_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of model vertical grid levels
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      precip_frac,         & ! Precipitation fraction (overall)              [-]
      cloud_frac_1,        & ! Cloud fraction (1st PDF component)            [-]
      cloud_frac_2,        & ! Cloud fraction (2nd PDF component)            [-]
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2, & ! Ice supersaturation fraction (2nd PDF comp.)  [-]
      mixt_frac              ! Mixture fraction                              [-]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      weighted_pfrac_1, & ! Product of mixt_frac and cloud_frac_1           [-]
      weighted_pfrac_2    ! Product of ( 1 - mixt_frac ) and cloud_frac_2   [-]

    integer :: k, j  ! Loop index


    !!! Find precipitation fraction within PDF component 1.
    ! The method used to find overall precipitation fraction will also be to
    ! find precipitation fraction within PDF component 1 and PDF component 2.
    ! In order to do so, it is assumed (poorly) that PDF component 1 overlaps
    ! PDF component 1 and that PDF component 2 overlaps PDF component 2 at every
    ! vertical level in the vertical profile.  
    
    ! The weighted precipitation fraction in PDF component 1 is the greatest
    ! value of the product of mixture fraction (mixt_frac) and 1st PDF
    ! component cloud fraction at or above a vertical level.  Likewise, the
    ! weighted precipitation fraction in PDF component 2 is the greatest
    ! value of the product of ( 1 - mixt_frac ) and 2nd PDF component cloud
    ! fraction at or above a vertical level.
    if ( any( l_frozen_hm ) ) then
      
      ! Ice microphysics included.

      ! Weighted precipitation fraction in PDF component 1.
      weighted_pfrac_1(:,nz) &
      = max( mixt_frac(:,nz) * cloud_frac_1(:,nz), &
             mixt_frac(:,nz) * ice_supersat_frac_1(:,nz) )

      ! Weighted precipitation fraction in PDF component 2.
      weighted_pfrac_2(:,nz) &
      = max( ( one - mixt_frac(:,nz) ) * cloud_frac_2(:,nz), &
             ( one - mixt_frac(:,nz) ) * ice_supersat_frac_2(:,nz) )
             
      do k = nz, 1, -1
        do j = 1, ngrdcol
          ! Weighted precipitation fraction in PDF component 1.
          weighted_pfrac_1(j,k) &
          = max( weighted_pfrac_1(j,k+1), &
                 mixt_frac(j,k) * cloud_frac_1(j,k), &
                 mixt_frac(j,k) * ice_supersat_frac_1(j,k) )

          ! Weighted precipitation fraction in PDF component 2.
          weighted_pfrac_2(j,k) &
          = max( weighted_pfrac_2(j,k+1), &
                 ( one - mixt_frac(j,k) ) * cloud_frac_2(j,k), &
                 ( one - mixt_frac(j,k) ) * ice_supersat_frac_2(j,k) )
        end do
      end do
      
      
    else
      
       ! Warm microphysics.

       ! Weighted precipitation fraction in PDF component 1.
       weighted_pfrac_1(:,nz) = mixt_frac(:,nz) * cloud_frac_1(:,nz)

       ! Weighted precipitation fraction in PDF component 2.
       weighted_pfrac_2(:,nz) = ( one - mixt_frac(:,nz) ) * cloud_frac_2(:,nz)
       
      do k = nz, 1, -1
        do j = 1, ngrdcol
          ! Weighted precipitation fraction in PDF component 1.
          weighted_pfrac_1(j,k) &
          = max( weighted_pfrac_1(j,k+1), &
                 mixt_frac(j,k) * cloud_frac_1(j,k) )

          ! Weighted precipitation fraction in PDF component 2.
          weighted_pfrac_2(j,k) &
          = max( weighted_pfrac_2(j,k+1), &
                 ( one - mixt_frac(j,k) ) * cloud_frac_2(j,k) )
        end do
      end do
      
    end if
    

    ! Calculate precip_frac_1 and special cases for precip_frac_1.
    do k = 1, nz, 1
      do j = 1, ngrdcol
        ! Calculate precipitation fraction in the 1st PDF component.
        if ( weighted_pfrac_1(j,k) + weighted_pfrac_2(j,k) > zero ) then

          ! Adjust weighted 1st PDF component precipitation fraction by
          ! multiplying it by a factor.  That factor is overall precipitation
          ! fraction divided by the sum of the weighted 1st PDF component
          ! precipitation fraction and the weighted 2nd PDF component
          ! precipitation fraction.  The 1st PDF component precipitation
          ! fraction is then found by dividing the adjusted weighted 1st PDF
          ! component precipitation fraction by mixture fraction.
          precip_frac_1(j,k) &
          = weighted_pfrac_1(j,k) &
            * ( precip_frac(j,k) &
                / ( weighted_pfrac_1(j,k) + weighted_pfrac_2(j,k) ) ) &
            / mixt_frac(j,k)
        else

          ! Usually, the sum of the weighted 1st PDF component precipitation
          ! fraction and the 2nd PDF component precipitation fraction go to 0
          ! when overall precipitation fraction goes to 0.  Since 1st PDF
          ! component weighted precipitation fraction is 0, 1st PDF component
          ! precipitation fraction also 0.
          precip_frac_1(j,k) = zero

        end if
      end do
    end do
      
    do k = 1, nz, 1
      do j = 1, ngrdcol
        ! Special cases for precip_frac_1.
        if ( any( hydromet(j,k,:) >= hydromet_tol(:) ) &
            .and. precip_frac_1(j,k) &
                  > min( one, precip_frac(j,k) / mixt_frac(j,k) ) ) then

          ! Using the above method, it is possible for precip_frac_1 to be
          ! greater than 1.  For example, the mixture fraction at level k+1 is
          ! 0.10 and the cloud_frac_1 at level k+1 is 1, resulting in a
          ! weighted_pfrac_1 of 0.10.  This product is greater than the product
          ! of mixt_frac and cloud_frac_1 at level k.  The mixture fraction at
          ! level k is 0.05, resulting in a precip_frac_1 of 2.  The value of
          ! precip_frac_1 is limited at 1.  The leftover precipitation fraction
          ! (a result of the decreasing weight of PDF component 1 between the
          ! levels) is applied to PDF component 2.
          ! Additionally, when weighted_pfrac_1 at level k is greater than
          ! overall precipitation fraction at level k, the resulting calculation
          ! of precip_frac_2 at level k will be negative.
          precip_frac_1(j,k) = min( one, precip_frac(j,k) / mixt_frac(j,k) )

        elseif ( any( hydromet(j,k,:) >= hydromet_tol(:) ) &
                .and. precip_frac_1(j,k) > zero &
                .and. precip_frac_1(j,k) < precip_frac_tol(j) ) then

          ! In a scenario where we find precipitation in the 1st PDF component
          ! (it is allowed to have a value of 0 when all precipitation is found
          ! in the 2nd PDF component) but it is tiny (less than tolerance
          ! level), boost 1st PDF component precipitation fraction to tolerance
          ! level.
          precip_frac_1(j,k) = precip_frac_tol(j)

        elseif ( all( hydromet(j,k,:) < hydromet_tol(:) ) ) then

          ! The means (overall) of every precipitating hydrometeor are all less
          ! than their respective tolerance amounts.  They are all considered to
          ! have values of 0.  There are not any hydrometeor species found at
          ! this grid level.  There is also no cloud at or above this grid
          ! level, so set 1st component precipitation fraction to 0.
          precip_frac_1(j,k) = zero

        end if
      end do
    end do


    !!! Find precipitation fraction within PDF component 2.
    ! The equation for precipitation fraction within PDF component 2 is:
    !
    ! f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
    !
    ! given the overall precipitation fraction, f_p (calculated above), the
    ! precipitation fraction within PDF component 1, f_p(1) (calculated above),
    ! and mixture fraction, a.  Any leftover precipitation fraction from
    ! precip_frac_1 will be included in this calculation of precip_frac_2.
    do k = 1, nz, 1
      do j = 1, ngrdcol

        if ( any( hydromet(j,k,:) >= hydromet_tol(:) ) ) then

          ! Calculate precipitation fraction in the 2nd PDF component.
          precip_frac_2(j,k) &
          = max( ( precip_frac(j,k) - mixt_frac(j,k) * precip_frac_1(j,k) ) &
                 / ( one - mixt_frac(j,k) ), &
                 zero )

          ! Special cases for precip_frac_2.
          if ( precip_frac_2(j,k) > one ) then

             ! Again, it is possible for precip_frac_2 to be greater than 1.
             ! For example, the mixture fraction at level k+1 is 0.10 and the
             ! cloud_frac_1 at level k+1 is 1, resulting in a weighted_pfrac_1
             ! of 0.10.  This product is greater than the product of mixt_frac
             ! and cloud_frac_1 at level k.  Additionally, precip_frac (overall)
             ! is 1 for level k.  The mixture fraction at level k is 0.5,
             ! resulting in a precip_frac_1 of 0.2.  Using the above equation,
             ! precip_frac_2 is calculated to be 1.8.  The value of
             ! precip_frac_2 is limited at 1.  The leftover precipitation
             ! fraction (as a result of the increasing weight of component 1
             ! between the levels) is applied to PDF component 1.
             precip_frac_2(j,k) = one

             ! Recalculate precipitation fraction in the 1st PDF component.
             precip_frac_1(j,k) &
             = ( precip_frac(j,k) - ( one - mixt_frac(j,k) ) ) / mixt_frac(j,k)

             ! Double check precip_frac_1
             if ( precip_frac_1(j,k) > one ) then

                precip_frac_1(j,k) = one

                precip_frac_2(j,k) = ( precip_frac(j,k) - mixt_frac(j,k) ) &
                                   / ( one - mixt_frac(j,k) )

             elseif ( precip_frac_1(j,k) > zero .and. precip_frac_1(j,k) < precip_frac_tol(j) ) then

                precip_frac_1(j,k) = precip_frac_tol(j)

                ! fp = a*fp1+(1-a)*fp2 solving for fp2
                precip_frac_2(j,k) = precip_frac_1(j,k) &
                                     * ( ( ( precip_frac(j,k) / precip_frac_1(j,k)) &
                                         - mixt_frac(j,k) ) / ( one - mixt_frac(j,k) ) )
                
             end if

          elseif ( precip_frac_2(j,k) > zero &
                   .and. precip_frac_2(j,k) < precip_frac_tol(j) ) then

            ! In a scenario where we find precipitation in the 2nd PDF
            ! component (it is allowed to have a value of 0 when all
            ! precipitation is found in the 1st PDF component) but it is tiny
            ! (less than tolerance level), boost 2nd PDF component
            ! precipitation fraction to tolerance level.
            precip_frac_2(j,k) = precip_frac_tol(j)

            ! Recalculate precipitation fraction in the 1st PDF component.
            precip_frac_1(j,k) &
             = ( precip_frac(j,k) - ( one - mixt_frac(j,k) ) * precip_frac_2(j,k) ) &
               / mixt_frac(j,k)

            ! Double check precip_frac_1
            if ( precip_frac_1(j,k) > one ) then

              precip_frac_1(j,k) = one

              precip_frac_2(j,k) = ( precip_frac(j,k) - mixt_frac(j,k) ) &
                                   / ( one - mixt_frac(j,k) )

            elseif ( precip_frac_1(j,k) > zero .and. precip_frac_1(j,k) < precip_frac_tol(j) ) then

              precip_frac_1(j,k) = precip_frac_tol(j) 

              ! fp = a*fp1+(1-a)*fp2 solving for fp2
              precip_frac_2(j,k) = precip_frac_1(j,k) &
                                   * ( ( ( precip_frac(j,k) / precip_frac_1(j,k)) &
                                       - mixt_frac(j,k) ) / ( one - mixt_frac(j,k) ) )

            end if

          end if ! Special cases for precip_frac_2

        else ! all( hydromet(k,:) < hydromet_tol(:) )

          ! The means (overall) of every precipitating hydrometeor are all less
          ! than their respective tolerance amounts.  They are all considered to
          ! have values of 0.  There are not any hydrometeor species found at
          ! this grid level.  There is also no cloud at or above this grid
          ! level, so set 2nd component precipitation fraction to 0.
          precip_frac_2(j,k) = zero

        end if ! any( hydromet(k,:) > hydromet_tol(:) )
      end do
    end do

    return

  end subroutine component_precip_frac_weighted

  !=============================================================================
  subroutine component_precip_frac_specify( nz, ngrdcol, &
                                            upsilon_precip_frac_rat, &
                                            hydromet, precip_frac, &
                                            mixt_frac, precip_frac_tol, &
                                            precip_frac_1, precip_frac_2 )

    ! Description:
    ! Calculates the precipitation fraction in each PDF component.
    !
    ! The equation for precipitation fraction is:
    !
    ! f_p = mixt_frac * f_p(1) + ( 1 - mixt_frac ) * f_p(2);
    !
    ! where f_p is overall precipitation fraction, f_p(1) is precipitation
    ! fraction in the 1st PDF component, f_p(2) is precipitation fraction in the
    ! 2nd PDF component, and mixt_frac is the mixture fraction.  Using this
    ! method, a new specified parameter is introduced, upsilon, where:
    !
    ! upsilon = mixt_frac * f_p(1) / f_p; and where 0 <= upsilon <= 1.
    !
    ! In other words, upsilon is the ratio of mixt_frac * f_p(1) to f_p.  Since
    ! f_p and mixt_frac are calculated previously, and upsilon is specified,
    ! f_p(1) can be calculated by:
    !
    ! f_p(1) = upsilon * f_p / mixt_frac;
    !
    ! and has an upper limit of 1.  The value of f_p(2) can then be calculated
    ! by:
    !
    ! f_p(2) = ( f_p - mixt_frac * f_p(1) ) / ( 1 - mixt_frac );
    !
    ! and also has an upper limit of 1.  When upsilon = 1, all of the
    ! precipitation is found in the 1st PDF component (as long as
    ! f_p <= mixt_frac, otherwise it would cause f_p(1) to be greater than 1).
    ! When upsilon = 0, all of the precipitation is found in the 2nd PDF
    ! component (as long as f_p <= 1 - mixt_frac, otherwise it would cause
    ! f_p(2) to be greater than 1).  When upsilon is between 0 and 1,
    ! precipitation is split between the two PDF components accordingly.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero, &
        eps

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        hydromet_tol  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of model vertical grid levels
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), intent(in) :: &
      upsilon_precip_frac_rat    ! ratio mixt_frac*precip_frac_1/precip_frac [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      precip_frac, & ! Precipitation fraction (overall)                      [-]
      mixt_frac      ! Mixture fraction                                      [-]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    integer :: k, j  ! Loop index.


    ! There are hydrometeors found at this grid level.
    if ( abs(upsilon_precip_frac_rat - one) < &
          abs(upsilon_precip_frac_rat + one) / 2 * eps ) then


      ! Loop over all vertical levels.
      do k = 1, nz, 1
        do j = 1, ngrdcol

          if ( any( hydromet(j,k,:) >= hydromet_tol(:) ) ) then
             
            if ( precip_frac(j,k) <= mixt_frac(j,k) ) then

              ! All the precipitation is found in the 1st PDF component.
              precip_frac_1(j,k) = precip_frac(j,k) / mixt_frac(j,k)
              precip_frac_2(j,k) = zero

            else ! precip_frac(k) > mixt_frac(k)

              ! Some precipitation is found in the 2nd PDF component.
              precip_frac_1(j,k) = one
              precip_frac_2(j,k) = ( precip_frac(j,k) - mixt_frac(j,k) ) &
                                 / ( one - mixt_frac(j,k) )

              if ( precip_frac_2(j,k) > one &
                   .and. abs(precip_frac(j,k) - one) < abs(precip_frac(j,k) + one) / 2 * eps ) then

                 ! Set precip_frac_2 = 1.
                 precip_frac_2(j,k) = one

              elseif ( precip_frac_2(j,k) < precip_frac_tol(j) ) then

                ! Since precipitation is found in the 2nd PDF component, it
                ! must have a value of at least precip_frac_tol.
                precip_frac_2(j,k) = precip_frac_tol(j)

                ! Recalculate precip_frac_1
                precip_frac_1(j,k) &
                 = ( precip_frac(j,k) &
                     - ( one - mixt_frac(j,k) ) * precip_frac_2(j,k) ) &
                   / mixt_frac(j,k)

                ! Double check precip_frac_1
                if ( precip_frac_1(j,k) > one ) then

                  precip_frac_1(j,k) = one

                  precip_frac_2(j,k) = ( precip_frac(j,k) - mixt_frac(j,k) ) &
                                       / ( one - mixt_frac(j,k) )

                elseif ( precip_frac_1(j,k) < precip_frac_tol(j) ) then

                  precip_frac_1(j,k) = precip_frac_tol(j)

                  ! fp = a*fp1+(1-a)*fp2 solving for fp2
                  precip_frac_2(j,k) = precip_frac_1(j,k) * &
                                       ( ( ( precip_frac(j,k) / precip_frac_1(j,k)) &
                                       - mixt_frac(j,k) ) / ( one - mixt_frac(j,k) ) )

                end if

              end if ! precip_frac_2(k) < precip_frac_tol

            endif ! precip_frac(k) <= mixt_frac(k)
            
          else ! all( hydromet(k,:) < hydromet_tol(:) )
            ! There aren't any hydrometeors found at the grid level.
            precip_frac_1(j,k) = zero
            precip_frac_2(j,k) = zero
          end if
          
        end do
      end do

    elseif ( abs(upsilon_precip_frac_rat - zero) < &
          abs(upsilon_precip_frac_rat + zero) / 2 * eps ) then

      do k = 1, nz, 1
        do j = 1, ngrdcol

          if ( any( hydromet(j,k,:) >= hydromet_tol(:) ) ) then
            
            if ( precip_frac(j,k) <= ( one - mixt_frac(j,k) ) ) then

              ! All the precipitation is found in the 2nd PDF component.
              precip_frac_1(j,k) = zero
              precip_frac_2(j,k) = precip_frac(j,k) / ( one - mixt_frac(j,k) )

            else ! precip_frac(k) > ( 1 - mixt_frac(k) )

              ! Some precipitation is found in the 1st PDF component.
              precip_frac_1(j,k) = ( precip_frac(j,k) - ( one - mixt_frac(j,k) ) ) &
                                    / mixt_frac(j,k)
              precip_frac_2(j,k) = one

              if ( precip_frac_1(j,k) > one &
                   .and. abs(precip_frac(j,k) - one) < abs(precip_frac(j,k) + one) / 2 * eps ) then

                ! Set precip_frac_1 = 1.
                precip_frac_1(j,k) = one

              elseif ( precip_frac_1(j,k) < precip_frac_tol(j) ) then

                ! Since precipitation is found in the 1st PDF component, it
                ! must have a value of at least precip_frac_tol.
                precip_frac_1(j,k) = precip_frac_tol(j)

                ! Recalculate precip_frac_2
                precip_frac_2(j,k) = ( precip_frac(j,k) &
                                       - mixt_frac(j,k) * precip_frac_1(j,k) ) &
                                      / ( one - mixt_frac(j,k) )

                ! Double check precip_frac_2
                if ( precip_frac_2(j,k) > one ) then

                  precip_frac_2(j,k) = one

                  precip_frac_1(j,k) = ( ( precip_frac(j,k) - one ) + mixt_frac(j,k) ) &
                                       / mixt_frac(j,k)

                elseif ( precip_frac_2(j,k) < precip_frac_tol(j) ) then

                  precip_frac_2(j,k) = precip_frac_tol(j)

                  ! fp = a*fp1+(1-a)*fp2 solving for fp1
                  precip_frac_1(j,k) = ( precip_frac(j,k) - precip_frac_2(j,k) ) / mixt_frac(j,k) &
                                       + precip_frac_2(j,k)

                endif

              end if ! precip_frac_1(k) < precip_frac_tol

            endif ! precip_frac(k) <= ( 1 - mixt_frac(k) )
            
          else ! all( hydromet(k,:) < hydromet_tol(:) )
            ! There aren't any hydrometeors found at the grid level.
            precip_frac_1(j,k) = zero
            precip_frac_2(j,k) = zero
          end if
          
        end do
      end do

    else  ! 0 < upsilon_precip_frac_rat < 1

      do k = 1, nz, 1
        do j = 1, ngrdcol

          if ( any( hydromet(j,k,:) >= hydromet_tol(:) ) ) then
            ! Precipitation is found in both PDF components.  Each component
            ! must have a precipitation fraction that is at least
            ! precip_frac_tol and that does not exceed 1.

            ! Calculate precipitation fraction in the 1st PDF component.
            precip_frac_1(j,k) &
            = upsilon_precip_frac_rat * precip_frac(j,k) / mixt_frac(j,k)

            ! Special cases for precip_frac_1
            if ( precip_frac_1(j,k) > one ) then
              precip_frac_1(j,k) = one
            elseif ( precip_frac_1(j,k) < precip_frac_tol(j) ) then
              precip_frac_1(j,k) = precip_frac_tol(j)
            endif

            ! Calculate precipitation fraction in the 2nd PDF component.
            precip_frac_2(j,k) = ( precip_frac(j,k) &
                                   - mixt_frac(j,k) * precip_frac_1(j,k) ) &
                                 / ( one - mixt_frac(j,k) )

            ! Special case for precip_frac_2
            if ( precip_frac_2(j,k) > one ) then

              ! Set precip_frac_2 to 1.
              precip_frac_2(j,k) = one

              ! Recalculate precipitation fraction in the 1st PDF component.
              precip_frac_1(j,k) &
              = ( precip_frac(j,k) - ( one - mixt_frac(j,k) ) ) / mixt_frac(j,k)

              ! Double check precip_frac_1
              if ( precip_frac_1(j,k) > one ) then

                precip_frac_1(j,k) = one

                precip_frac_2(j,k) = ( precip_frac(j,k) - mixt_frac(j,k) ) &
                                     / ( one - mixt_frac(j,k) )

              elseif ( precip_frac_1(j,k) < precip_frac_tol(j) ) then

                precip_frac_1(j,k) = precip_frac_tol(j)

                ! fp = a*fp1+(1-a)*fp2 solving for fp2
                precip_frac_2(j,k) = precip_frac_1(j,k) &
                                     * ( ( ( precip_frac(j,k) / precip_frac_1(j,k)) &
                                     - mixt_frac(j,k) ) / ( one - mixt_frac(j,k) ) )

              endif

            elseif ( precip_frac_2(j,k) < precip_frac_tol(j) ) then

              ! Set precip_frac_2 to precip_frac_tol.
              precip_frac_2(j,k) = precip_frac_tol(j)

              ! Recalculate precipitation fraction in the 1st PDF component.
              precip_frac_1(j,k) = ( precip_frac(j,k) &
                                     - ( one - mixt_frac(j,k) ) * precip_frac_2(j,k) ) &
                                   / mixt_frac(j,k)

              ! Double check precip_frac_1
              if ( precip_frac_1(j,k) > one ) then

                precip_frac_1(j,k) = one

                precip_frac_2(j,k) = ( precip_frac(j,k) - mixt_frac(j,k) ) &
                                     / ( one - mixt_frac(j,k) )

              elseif ( precip_frac_1(j,k) < precip_frac_tol(j) ) then

                precip_frac_1(j,k) = precip_frac_tol(j)

                ! fp = a*fp1+(1-a)*fp2 solving for fp2
                precip_frac_2(j,k) = precip_frac_1(j,k) &
                                     * ( ( ( precip_frac(j,k) / precip_frac_1(j,k)) &
                                     - mixt_frac(j,k) ) / ( one - mixt_frac(j,k) ) )

              end if

            end if ! Special cases for precip_frac_2
         
          else ! all( hydromet(k,:) < hydromet_tol(:) )
            ! There aren't any hydrometeors found at the grid level.
            precip_frac_1(j,k) = zero
            precip_frac_2(j,k) = zero
          end if
          
        end do
      end do

    end if  ! upsilon_precip_frac_rat

    return

  end subroutine component_precip_frac_specify

  !=============================================================================
  subroutine precip_frac_assert_check( nz, hydromet, mixt_frac, precip_frac, &
                                       precip_frac_1, precip_frac_2, &
                                       precip_frac_tol )

    ! Description:
    ! Assertion check for the precipitation fraction code.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,     & ! Constant(s)
        zero,    &
        fstderr, &
        eps

    use array_index, only: &
        hydromet_tol  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use error_code, only: &
        err_code, &                     ! Error Indicator
        clubb_fatal_error               ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      mixt_frac,     & ! Mixture fraction                               [-]
      precip_frac,   & ! Precipitation fraction (overall)               [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    ! Local Variables
    integer :: k  ! Loop index


    ! Loop over all vertical levels.
    do k = 1, nz, 1

       if ( any( hydromet(k,:) >= hydromet_tol(:) ) ) then

          ! Overall precipitation fraction cannot be less than precip_frac_tol
          ! when a hydrometeor is present at a grid level.
          if ( precip_frac(k) < precip_frac_tol ) then
             write(fstderr,*) "precip_frac < precip_frac_tol when " &
                              // "a hydrometeor is present"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac = ", precip_frac(k), &
                              "precip_frac_tol = ", precip_frac_tol
             err_code = clubb_fatal_error
             return
          endif

          ! Overall precipitation fraction cannot exceed 1.
          if ( precip_frac(k) > one ) then
             write(fstderr,*) "precip_frac > 1"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac = ", precip_frac(k)
             err_code = clubb_fatal_error
             return
          endif

          ! Precipitation fraction in the 1st PDF component is allowed to be 0
          ! when all the precipitation is found in the 2nd PDF component.
          ! Otherwise, it cannot be less than precip_frac_tol when a hydrometeor
          ! is present at a grid level.  In other words, it cannot have a value
          ! that is greater than 0 but less than precip_frac_tol
          if ( precip_frac_1(k) > zero &
               .and. precip_frac_1(k) < precip_frac_tol ) then
             write(fstderr,*) "0 < precip_frac_1 < precip_frac_tol"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_1 = ", precip_frac_1(k), &
                              "precip_frac_tol = ", precip_frac_tol
             err_code = clubb_fatal_error
             return
          endif

          ! Precipitation fraction in the 1st PDF component cannot exceed 1.
          if ( precip_frac_1(k) > one ) then
             write(fstderr,*) "precip_frac_1 > 1"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_1 = ", precip_frac_1(k)
             err_code = clubb_fatal_error
             return
          endif

          ! Precipiation fraction in the 1st PDF component cannot be negative.
          if ( precip_frac_1(k) < zero ) then
             write(fstderr,*) "precip_frac_1 < 0"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_1 = ", precip_frac_1(k)
             err_code = clubb_fatal_error
             return
          endif

          ! Precipitation fraction in the 2nd PDF component is allowed to be 0
          ! when all the precipitation is found in the 1st PDF component.
          ! Otherwise, it cannot be less than precip_frac_tol when a hydrometeor
          ! is present at a grid level.  In other words, it cannot have a value
          ! that is greater than 0 but less than precip_frac_tol
          if ( precip_frac_2(k) > zero &
               .and. precip_frac_2(k) < precip_frac_tol ) then
             write(fstderr,*) "0 < precip_frac_2 < precip_frac_tol"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_2 = ", precip_frac_2(k), &
                              "precip_frac_tol = ", precip_frac_tol
             err_code = clubb_fatal_error
             return
          endif

          ! Precipitation fraction in the 2nd PDF component cannot exceed 1.
          if ( precip_frac_2(k) > one ) then
             write(fstderr,*) "precip_frac_2 > 1"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_2 = ", precip_frac_2(k)
             err_code = clubb_fatal_error
             return
          endif

          ! Precipiation fraction in the 2nd PDF component cannot be negative.
          if ( precip_frac_2(k) < zero ) then
             write(fstderr,*) "precip_frac_2 < 0"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_2 = ", precip_frac_2(k)
             err_code = clubb_fatal_error
             return
          endif

       else  ! all( hydromet(k,:) < hydromet_tol(:) )

          ! Overall precipitation fraction must be 0 when no hydrometeors are
          ! found at a grid level.
          if ( abs(precip_frac(k)) > eps) then
             write(fstderr,*) "precip_frac /= 0 when no hydrometeors are found"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac = ", precip_frac(k)
             err_code = clubb_fatal_error
             return
          endif

          ! Precipitation fraction in the 1st PDF component must be 0 when no
          ! hydrometeors are found at a grid level.
          if ( abs(precip_frac_1(k)) > eps) then
             write(fstderr,*) "precip_frac_1 /= 0 when no hydrometeors " &
                              // "are found"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_1 = ", precip_frac_1(k)
             err_code = clubb_fatal_error
             return
          endif

          ! Precipitation fraction in the 2nd PDF component must be 0 when no
          ! hydrometeors are found at a grid level.
          if ( abs(precip_frac_2(k)) > eps) then
             write(fstderr,*) "precip_frac_2 /= 0 when no hydrometeors " &
                              // "are found"
             write(fstderr,*) "level = ", k
             write(fstderr,*) "precip_frac_2 = ", precip_frac_2(k)
             err_code = clubb_fatal_error
             return
          endif

       endif  ! any( hydromet(k,:) >= hydromet_tol(:) )

       ! The precipitation fraction equation is:
       !
       ! precip_frac
       !    = mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2;
       !
       ! which means that:
       !
       ! precip_frac
       ! - ( mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2 )
       ! = 0.
       !
       ! Check that this is true with numerical round off.
       if ( ( precip_frac(k) &
              - ( mixt_frac(k) * precip_frac_1(k) &
                  + ( one - mixt_frac(k) ) * precip_frac_2(k) ) ) &
            > ( epsilon( precip_frac(k) ) * precip_frac(k) ) ) then
          write(fstderr,*) "mixt_frac * precip_frac_1 " &
                           // "+ ( 1 - mixt_frac ) * precip_frac_2 " &
                           // "/= precip_frac within numerical roundoff"
          write(fstderr,*) "level = ", k
          write(fstderr,*) "mixt_frac * precip_frac_1 " &
                           // "+ ( 1 - mixt_frac ) * precip_frac_2 = ", &
                           mixt_frac(k) * precip_frac_1(k) &
                           + ( one - mixt_frac(k) ) * precip_frac_2(k)
          write(fstderr,*) "precip_frac = ", precip_frac(k)
          err_code = clubb_fatal_error
          return
       endif

    enddo  ! k = 1, nz, 1


    return

  end subroutine precip_frac_assert_check

!===============================================================================

end module precipitation_fraction
