!-------------------------------------------------------------------------------
! $Id: clip_explicit.F90 7315 2014-09-30 20:49:54Z schemena@uwm.edu $
!===============================================================================
module clip_explicit

  implicit none

  private

  public :: clip_covars_denom, &
            clip_covar, & 
            clip_covar_level, & 
            clip_variance, & 
            clip_skewness, &
            clip_skewness_core

  ! Named constants to avoid string comparisons
  integer, parameter, public :: &
    clip_rtp2 = 1, &         ! Named constant for rtp2 clipping
    clip_thlp2 = 2, &        ! Named constant for thlp2 clipping
    clip_rtpthlp = 3, &      ! Named constant for rtpthlp clipping
    clip_up2 = 5, &          ! Named constant for up2 clipping
    clip_vp2 = 6, &          ! Named constant for vp2 clipping
!    clip_scalar = 7, &       ! Named constant for scalar clipping
    clip_wprtp = 8, &        ! Named constant for wprtp clipping
    clip_wpthlp = 9, &       ! Named constant for wpthlp clipping
    clip_upwp = 10, &        ! Named constant for upwp clipping
    clip_vpwp = 11, &        ! Named constant for vpwp clipping
    clip_wp2 = 12, &         ! Named constant for wp2 clipping
    clip_wpsclrp = 13, &     ! Named constant for wp scalar clipping
    clip_sclrp2 = 14, &      ! Named constant for sclrp2 clipping
    clip_sclrprtp = 15, &    ! Named constant for sclrprtp clipping
    clip_sclrpthlp = 16, &   ! Named constant for sclrpthlp clipping
    clip_wphydrometp = 17    ! Named constant for wphydrometp clipping

  contains

  !=============================================================================
  subroutine clip_covars_denom( dt, rtp2, thlp2, up2, vp2, wp2, &
                                sclrp2, wprtp_cl_num, wpthlp_cl_num, &
                                wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num, &
                                wprtp, wpthlp, upwp, vpwp, wpsclrp )

    ! Description:
    ! Some of the covariances found in the CLUBB model code need to be clipped
    ! multiple times during each timestep to ensure that the correlation between
    ! the two relevant variables stays between -1 and 1 at all times during the
    ! model run.  The covariances that need to be clipped multiple times are
    ! w'r_t', w'th_l', w'sclr', u'w', and v'w'.  One of the times that each one
    ! of these covariances is clipped is immediately after each one is set.
    ! However, each covariance still needs to be clipped two more times during
    ! each timestep (once after advance_xp2_xpyp is called and once after
    ! advance_wp2_wp3 is called).  This subroutine handles the times that the
    ! covariances are clipped away from the time that they are set.  In other
    ! words, this subroutine clips the covariances after the denominator terms
    ! in the relevant correlation equation have been altered, ensuring that
    ! all correlations will remain between -1 and 1 at all times.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr ! Variable(s)

    use parameters_model, only: &
        sclr_dim ! Variable(s)

    use model_flags, only: &
        l_tke_aniso ! Logical

    use clubb_precision, only: & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: &
       stat_modify ! Procedure(s)

    use stats_variables, only: & 
        iwprtp_bt, &  ! Variable(s)
        iwpthlp_bt, &
        stats_zm, &
        l_stats_samp

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      dt ! Timestep [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rtp2,  & ! r_t'^2         [(kg/kg)^2]
      thlp2, & ! theta_l'^2     [K^2]
      up2,   & ! u'^2           [m^2/s^2]
      vp2,   & ! v'^2           [m^2/s^2]
      wp2      ! w'^2           [m^2/s^2]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(in) :: &
      sclrp2 ! sclr'^2  [{units vary}^2]

    integer, intent(in) :: &
      wprtp_cl_num,   &
      wpthlp_cl_num,  &
      wpsclrp_cl_num, &
      upwp_cl_num,    &
      vpwp_cl_num

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      wprtp,  & ! w'r_t'        [(kg/kg) m/s]
      wpthlp, & ! w'theta_l'    [K m/s]
      upwp,   & ! u'w'          [m^2/s^2]
      vpwp      ! v'w'          [m^2/s^2]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(inout) :: &
      wpsclrp ! w'sclr'         [units m/s]

    ! Local Variables
    logical :: & 
      l_first_clip_ts, & ! First instance of clipping in a timestep.
      l_last_clip_ts     ! Last instance of clipping in a timestep.

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wprtp_chnge,  & ! Net change in w'r_t' due to clipping  [(kg/kg) m/s]
      wpthlp_chnge, & ! Net change in w'th_l' due to clipping [K m/s]
      upwp_chnge,   & ! Net change in u'w' due to clipping    [m^2/s^2]
      vpwp_chnge      ! Net change in v'w' due to clipping    [m^2/s^2]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      wpsclrp_chnge   ! Net change in w'sclr' due to clipping [{units vary}]

    integer :: i  ! scalar array index.

    ! ---- Begin Code ----

    !!! Clipping for w'r_t'
    !
    ! Clipping w'r_t' at each vertical level, based on the
    ! correlation of w and r_t at each vertical level, such that:
    ! corr_(w,r_t) = w'r_t' / [ sqrt(w'^2) * sqrt(r_t'^2) ];
    ! -1 <= corr_(w,r_t) <= 1.
    !
    ! Since w'^2, r_t'^2, and w'r_t' are each advanced in different
    ! subroutines from each other in advance_clubb_core, clipping for w'r_t'
    ! is done three times during each timestep (once after each variable has
    ! been updated).
    !
    ! This subroutine handles the first and third instances of
    ! w'r_t' clipping.
    ! The first instance of w'r_t' clipping takes place after
    ! r_t'^2 is updated in advance_xp2_xpyp.
    ! The third instance of w'r_t' clipping takes place after
    ! w'^2 is updated in advance_wp2_wp3.

    ! Include effect of clipping in wprtp time tendency budget term.
    if ( l_stats_samp ) then
    
      ! if wprtp_cl_num == 1 do nothing since
      ! iwprtp_bt stat_begin_update is called outside of this method
      
      if ( wprtp_cl_num == 2 ) then
        ! wprtp total time tendency (effect of clipping)
        call stat_modify( iwprtp_bt,  -wprtp / dt,  & ! intent(in)
                          stats_zm )                               ! intent(inout)
      elseif ( wprtp_cl_num == 3 ) then
        ! wprtp total time tendency (effect of clipping)
        call stat_modify( iwprtp_bt, -wprtp / dt,  & ! intent(in)
                          stats_zm )                               ! intent(inout)
      endif
    endif

    ! Used within subroutine clip_covar.
    if ( wprtp_cl_num == 1 ) then
      l_first_clip_ts = .true.
      l_last_clip_ts  = .false.
    elseif ( wprtp_cl_num == 2 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .false.
    elseif ( wprtp_cl_num == 3 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .true.
    endif

    ! Clip w'r_t'
    call clip_covar( clip_wprtp, l_first_clip_ts,   & ! intent(in) 
                     l_last_clip_ts, dt, wp2, rtp2, & ! intent(in)
                     wprtp, wprtp_chnge )             ! intent(inout)

    if ( l_stats_samp ) then
      if ( wprtp_cl_num == 1 ) then
        ! wprtp total time tendency (effect of clipping)
        call stat_modify( iwprtp_bt,  wprtp / dt,  & ! intent(in)
                          stats_zm )                              ! intent(inout)
      elseif ( wprtp_cl_num == 2 ) then
        ! wprtp total time tendency (effect of clipping)
        call stat_modify( iwprtp_bt, wprtp / dt,  & ! intent(in)
                          stats_zm )                              ! intent(inout)
      ! if wprtp_cl_num == 3 do nothing since
      ! iwprtp_bt stat_end_update is called outside of this method
      
      endif
    endif


    !!! Clipping for w'th_l'
    !
    ! Clipping w'th_l' at each vertical level, based on the
    ! correlation of w and th_l at each vertical level, such that:
    ! corr_(w,th_l) = w'th_l' / [ sqrt(w'^2) * sqrt(th_l'^2) ];
    ! -1 <= corr_(w,th_l) <= 1.
    !
    ! Since w'^2, th_l'^2, and w'th_l' are each advanced in different
    ! subroutines from each other in advance_clubb_core, clipping for w'th_l'
    ! is done three times during each timestep (once after each variable has
    ! been updated).
    !
    ! This subroutine handles the first and third instances of
    ! w'th_l' clipping.
    ! The first instance of w'th_l' clipping takes place after
    ! th_l'^2 is updated in advance_xp2_xpyp.
    ! The third instance of w'th_l' clipping takes place after
    ! w'^2 is updated in advance_wp2_wp3.

    ! Include effect of clipping in wpthlp time tendency budget term.
    if ( l_stats_samp ) then
    
      ! if wpthlp_cl_num == 1 do nothing since
      ! iwpthlp_bt stat_begin_update is called outside of this method
      
      if ( wpthlp_cl_num == 2 ) then
        ! wpthlp total time tendency (effect of clipping)
        call stat_modify( iwpthlp_bt, -wpthlp / dt,  & ! intent(in)
                          stats_zm )                                 ! intent(inout)
      elseif ( wpthlp_cl_num == 3 ) then
        ! wpthlp total time tendency (effect of clipping)
        call stat_modify( iwpthlp_bt, -wpthlp / dt,  & ! intent(in)
                          stats_zm )                                 ! intent(inout)
      endif
    endif

    ! Used within subroutine clip_covar.
    if ( wpthlp_cl_num == 1 ) then
      l_first_clip_ts = .true.
      l_last_clip_ts  = .false.
    elseif ( wpthlp_cl_num == 2 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .false.
    elseif ( wpthlp_cl_num == 3 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .true.
    endif

    ! Clip w'th_l'
    call clip_covar( clip_wpthlp, l_first_clip_ts,   & ! intent(in)
                     l_last_clip_ts, dt, wp2, thlp2, & ! intent(in)
                     wpthlp, wpthlp_chnge )            ! intent(inout)


    if ( l_stats_samp ) then
      if ( wpthlp_cl_num == 1 ) then
        ! wpthlp total time tendency (effect of clipping)
        call stat_modify( iwpthlp_bt, wpthlp / dt,  & ! intent(in)
                          stats_zm )                                ! intent(inout)
      elseif ( wpthlp_cl_num == 2 ) then
        ! wpthlp total time tendency (effect of clipping)
        call stat_modify( iwpthlp_bt, wpthlp / dt,  & ! intent(in)
                          stats_zm )                                ! intent(inout)
                          
      ! if wpthlp_cl_num == 3 do nothing since
      ! iwpthlp_bt stat_end_update is called outside of this method
      
      endif
    endif


    !!! Clipping for w'sclr'
    !
    ! Clipping w'sclr' at each vertical level, based on the
    ! correlation of w and sclr at each vertical level, such that:
    ! corr_(w,sclr) = w'sclr' / [ sqrt(w'^2) * sqrt(sclr'^2) ];
    ! -1 <= corr_(w,sclr) <= 1.
    !
    ! Since w'^2, sclr'^2, and w'sclr' are each advanced in different
    ! subroutines from each other in advance_clubb_core, clipping for w'sclr'
    ! is done three times during each timestep (once after each variable has
    ! been updated).
    !
    ! This subroutine handles the first and third instances of
    ! w'sclr' clipping.
    ! The first instance of w'sclr' clipping takes place after
    ! sclr'^2 is updated in advance_xp2_xpyp.
    ! The third instance of w'sclr' clipping takes place after
    ! w'^2 is updated in advance_wp2_wp3.

    ! Used within subroutine clip_covar.
    if ( wpsclrp_cl_num == 1 ) then
      l_first_clip_ts = .true.
      l_last_clip_ts  = .false.
    elseif ( wpsclrp_cl_num == 2 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .false.
    elseif ( wpsclrp_cl_num == 3 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .true.
    endif

    ! Clip w'sclr'
    do i = 1, sclr_dim, 1
      call clip_covar( clip_wpsclrp, l_first_clip_ts,           & ! intent(in)
                       l_last_clip_ts, dt, wp2(:), sclrp2(:,i), & ! intent(in)
                       wpsclrp(:,i), wpsclrp_chnge(:,i) )         ! intent(inout)
    enddo


    !!! Clipping for u'w'
    !
    ! Clipping u'w' at each vertical level, based on the
    ! correlation of u and w at each vertical level, such that:
    ! corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
    ! -1 <= corr_(u,w) <= 1.
    !
    ! Since w'^2, u'^2, and u'w' are each advanced in different
    ! subroutines from each other in advance_clubb_core, clipping for u'w'
    ! is done three times during each timestep (once after each variable has
    ! been updated).
    !
    ! This subroutine handles the first and second instances of
    ! u'w' clipping.
    ! The first instance of u'w' clipping takes place after
    ! u'^2 is updated in advance_xp2_xpyp.
    ! The second instance of u'w' clipping takes place after
    ! w'^2 is updated in advance_wp2_wp3.

    ! Used within subroutine clip_covar.
    if ( upwp_cl_num == 1 ) then
      l_first_clip_ts = .true.
      l_last_clip_ts  = .false.
    elseif ( upwp_cl_num == 2 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .false.
    elseif ( upwp_cl_num == 3 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .true.
    endif

    ! Clip u'w'
    if ( l_tke_aniso ) then
      call clip_covar( clip_upwp, l_first_clip_ts,   & ! intent(in)
                       l_last_clip_ts, dt, wp2, up2, & ! intent(in)
                       upwp, upwp_chnge )              ! intent(inout)
    else
      ! In this case, up2 = wp2, and the variable `up2' does not interact
      call clip_covar( clip_upwp, l_first_clip_ts,   & ! intent(in)
                       l_last_clip_ts, dt, wp2, wp2, & ! intent(in)
                       upwp, upwp_chnge )              ! intent(inout)
    end if



    !!! Clipping for v'w'
    !
    ! Clipping v'w' at each vertical level, based on the
    ! correlation of v and w at each vertical level, such that:
    ! corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
    ! -1 <= corr_(v,w) <= 1.
    !
    ! Since w'^2, v'^2, and v'w' are each advanced in different
    ! subroutines from each other in advance_clubb_core, clipping for v'w'
    ! is done three times during each timestep (once after each variable has
    ! been updated).
    !
    ! This subroutine handles the first and second instances of
    ! v'w' clipping.
    ! The first instance of v'w' clipping takes place after
    ! v'^2 is updated in advance_xp2_xpyp.
    ! The second instance of v'w' clipping takes place after
    ! w'^2 is updated in advance_wp2_wp3.

    ! Used within subroutine clip_covar.
    if ( vpwp_cl_num == 1 ) then
      l_first_clip_ts = .true.
      l_last_clip_ts  = .false.
    elseif ( vpwp_cl_num == 2 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .false.
    elseif ( vpwp_cl_num == 3 ) then
      l_first_clip_ts = .false.
      l_last_clip_ts  = .true.
    endif

    if ( l_tke_aniso ) then
      call clip_covar( clip_vpwp, l_first_clip_ts,   & ! intent(in)
                       l_last_clip_ts, dt, wp2, vp2, & ! intent(in)
                       vpwp, vpwp_chnge )              ! intent(inout)
    else
      ! In this case, vp2 = wp2, and the variable `vp2' does not interact
      call clip_covar( clip_vpwp, l_first_clip_ts,   & ! intent(in)
                       l_last_clip_ts, dt, wp2, wp2, & ! intent(in)
                       vpwp, vpwp_chnge )              ! intent(inout)
    end if


    return
  end subroutine clip_covars_denom

  !=============================================================================
  subroutine clip_covar( solve_type, l_first_clip_ts,  & 
                         l_last_clip_ts, dt, xp2, yp2,  & 
                         xpyp, xpyp_chnge )

    ! Description:
    ! Clipping the value of covariance x'y' based on the correlation between x
    ! and y.
    !
    ! The correlation between variables x and y is:
    !
    ! corr_(x,y) = x'y' / [ sqrt(x'^2) * sqrt(y'^2) ];
    !
    ! where x'^2 is the variance of x, y'^2 is the variance of y, and x'y' is
    ! the covariance of x and y.
    !
    ! The correlation of two variables must always have a value between -1
    ! and 1, such that:
    !
    ! -1 <= corr_(x,y) <= 1.
    !
    ! Therefore, there is an upper limit on x'y', such that:
    !
    ! x'y' <=  [ sqrt(x'^2) * sqrt(y'^2) ];
    !
    ! and a lower limit on x'y', such that:
    !
    ! x'y' >= -[ sqrt(x'^2) * sqrt(y'^2) ].
    !
    ! The values of x'y', x'^2, and y'^2 are all found on momentum levels.
    !
    ! The value of x'y' may need to be clipped whenever x'y', x'^2, or y'^2 is
    ! updated.
    !
    ! The following covariances are found in the code:
    !
    ! w'r_t', w'th_l', w'sclr', (computed in advance_xm_wpxp);
    ! r_t'th_l', sclr'r_t', sclr'th_l', (computed in advance_xp2_xpyp);
    ! u'w', v'w', w'edsclr' (computed in advance_windm_edsclrm);
    ! and w'hm' (computed in setup_pdf_parameters).

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only: &
        max_mag_correlation ! Constant(s)

    use clubb_precision, only: & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_begin_update,  & ! Procedure(s)
        stat_modify, & 
        stat_end_update

    use stats_variables, only: & 
        stats_zm,  & ! Variable(s)
        iwprtp_cl, & 
        iwpthlp_cl, & 
        irtpthlp_cl, & 
        l_stats_samp

    implicit none

    ! Input Variables
    integer, intent(in) :: & 
      solve_type       ! Variable being solved; used for STATS.

    logical, intent(in) :: & 
      l_first_clip_ts, & ! First instance of clipping in a timestep.
      l_last_clip_ts     ! Last instance of clipping in a timestep.

    real( kind = core_rknd ), intent(in) ::  & 
      dt     ! Model timestep; used here for STATS           [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      xp2, & ! Variance of x, x'^2 (momentum levels)         [{x units}^2]
      yp2    ! Variance of y, y'^2 (momentum levels)         [{y units}^2]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      xpyp   ! Covariance of x and y, x'y' (momentum levels) [{x units}*{y units}]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      xpyp_chnge  ! Net change in x'y' due to clipping [{x units}*{y units}]


    ! Local Variable
    integer :: k  ! Array index

    integer :: & 
      ixpyp_cl

    ! ---- Begin Code ----

    select case ( solve_type )
    case ( clip_wprtp )   ! wprtp clipping budget term
      ixpyp_cl = iwprtp_cl
    case ( clip_wpthlp )   ! wpthlp clipping budget term
      ixpyp_cl = iwpthlp_cl
    case ( clip_rtpthlp )   ! rtpthlp clipping budget term
      ixpyp_cl = irtpthlp_cl
    case default   ! scalars (or upwp/vpwp) are involved
      ixpyp_cl = 0
    end select


    if ( l_stats_samp ) then
      if ( l_first_clip_ts ) then
        call stat_begin_update( ixpyp_cl, xpyp / dt, stats_zm )
      else
        call stat_modify( ixpyp_cl, -xpyp / dt, stats_zm )
      endif
    endif

    ! The value of x'y' at the surface (or lower boundary) is a set value that
    ! is either specified or determined elsewhere in a surface subroutine.  It
    ! is ensured elsewhere that the correlation between x and y at the surface
    ! (or lower boundary) is between -1 and 1.  Thus, the covariance clipping
    ! code does not need to be invoked at the lower boundary.  Likewise, the
    ! value of x'y' is set at the upper boundary, so the covariance clipping
    ! code does not need to be invoked at the upper boundary.
    ! Note that if clipping were applied at the lower boundary, momentum will
    ! not be conserved, therefore it should never be added.
    do k = 2, gr%nz-1, 1

      ! Clipping for xpyp at an upper limit corresponding with a correlation
      ! between x and y of max_mag_correlation.
      if ( xpyp(k) >  max_mag_correlation * sqrt( xp2(k) * yp2(k) ) ) then

        xpyp_chnge(k) =  max_mag_correlation * sqrt( xp2(k) * yp2(k) ) - xpyp(k)

        xpyp(k) =  max_mag_correlation * sqrt( xp2(k) * yp2(k) )

        ! Clipping for xpyp at a lower limit corresponding with a correlation
        ! between x and y of -max_mag_correlation.
      elseif ( xpyp(k) < -max_mag_correlation * sqrt( xp2(k) * yp2(k) ) ) then

        xpyp_chnge(k) = -max_mag_correlation * sqrt( xp2(k) * yp2(k) ) - xpyp(k)

        xpyp(k) = -max_mag_correlation * sqrt( xp2(k) * yp2(k) )

      else

        xpyp_chnge(k) = 0.0_core_rknd

      endif

    enddo ! k = 2..gr%nz

    ! Since there is no covariance clipping at the upper or lower boundaries,
    ! the change in x'y' due to covariance clipping at those levels is 0.
    xpyp_chnge(1)       = 0.0_core_rknd
    xpyp_chnge(gr%nz) = 0.0_core_rknd

    if ( l_stats_samp ) then
      if ( l_last_clip_ts ) then
        call stat_end_update( ixpyp_cl, xpyp / dt, stats_zm )
      else
        call stat_modify( ixpyp_cl, xpyp / dt, stats_zm )
      endif
    endif


    return
  end subroutine clip_covar

  !=============================================================================
  subroutine clip_covar_level( solve_type, level, l_first_clip_ts,  & 
                               l_last_clip_ts, dt, xp2, yp2,  & 
                               xpyp, xpyp_chnge )

    ! Description:
    ! Clipping the value of covariance x'y' based on the correlation between x
    ! and y.  This is all done at a single vertical level.
    !
    ! The correlation between variables x and y is:
    !
    ! corr_(x,y) = x'y' / [ sqrt(x'^2) * sqrt(y'^2) ];
    !
    ! where x'^2 is the variance of x, y'^2 is the variance of y, and x'y' is
    ! the covariance of x and y.
    !
    ! The correlation of two variables must always have a value between -1
    ! and 1, such that:
    !
    ! -1 <= corr_(x,y) <= 1.
    !
    ! Therefore, there is an upper limit on x'y', such that:
    !
    ! x'y' <=  [ sqrt(x'^2) * sqrt(y'^2) ];
    !
    ! and a lower limit on x'y', such that:
    !
    ! x'y' >= -[ sqrt(x'^2) * sqrt(y'^2) ].
    !
    ! The values of x'y', x'^2, and y'^2 are all found on momentum levels.
    !
    ! The value of x'y' may need to be clipped whenever x'y', x'^2, or y'^2 is
    ! updated.
    !
    ! The following covariances are found in the code:
    !
    ! w'r_t', w'th_l', w'sclr', (computed in advance_xm_wpxp);
    ! r_t'th_l', sclr'r_t', sclr'th_l', (computed in advance_xp2_xpyp);
    ! u'w', v'w', w'edsclr' (computed in advance_windm_edsclrm);
    ! and w'hm' (computed in setup_pdf_parameters).

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        max_mag_correlation, & ! Constant(s)
        zero

    use clubb_precision, only: & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_begin_update_pt, & ! Procedure(s)
        stat_modify_pt,       & 
        stat_end_update_pt

    use stats_variables, only: & 
        stats_zm,  & ! Variable(s)
        iwprtp_cl, & 
        iwpthlp_cl, & 
        irtpthlp_cl, & 
        l_stats_samp

    implicit none

    ! Input Variables
    integer, intent(in) :: & 
      solve_type, & ! Variable being solved; used for STATS
      level         ! Vertical level index

    logical, intent(in) :: & 
      l_first_clip_ts, & ! First instance of clipping in a timestep.
      l_last_clip_ts     ! Last instance of clipping in a timestep.

    real( kind = core_rknd ), intent(in) ::  & 
      dt     ! Model timestep; used here for STATS        [s]

    real( kind = core_rknd ), intent(in) :: & 
      xp2, & ! Variance of x, <x'^2>                      [{x units}^2]
      yp2    ! Variance of y, <y'^2>                      [{y units}^2]

    ! Output Variable
    real( kind = core_rknd ), intent(inout) :: & 
      xpyp   ! Covariance of x and y, <x'y'>              [{x units}*{y units}]

    real( kind = core_rknd ), intent(out) :: &
      xpyp_chnge  ! Net change in <x'y'> due to clipping  [{x units}*{y units}]


    ! Local Variable
    integer :: & 
      ixpyp_cl    ! Statistics index


    select case ( solve_type )
    case ( clip_wprtp )   ! wprtp clipping budget term
      ixpyp_cl = iwprtp_cl
    case ( clip_wpthlp )   ! wpthlp clipping budget term
      ixpyp_cl = iwpthlp_cl
    case ( clip_rtpthlp )   ! rtpthlp clipping budget term
      ixpyp_cl = irtpthlp_cl
    case default   ! scalars (or upwp/vpwp) are involved
      ixpyp_cl = 0
    end select


    if ( l_stats_samp ) then
       if ( l_first_clip_ts ) then
          call stat_begin_update_pt( ixpyp_cl, level, &
                                     xpyp / dt, stats_zm )
       else
          call stat_modify_pt( ixpyp_cl, level, &
                               -xpyp / dt, stats_zm )
       endif
    endif

    ! The value of x'y' at the surface (or lower boundary) is a set value that
    ! is either specified or determined elsewhere in a surface subroutine.  It
    ! is ensured elsewhere that the correlation between x and y at the surface
    ! (or lower boundary) is between -1 and 1.  Thus, the covariance clipping
    ! code does not need to be invoked at the lower boundary.  Likewise, the
    ! value of x'y' is set at the upper boundary, so the covariance clipping
    ! code does not need to be invoked at the upper boundary.
    ! Note that if clipping were applied at the lower boundary, momentum will
    ! not be conserved, therefore it should never be added.

    ! Clipping for xpyp at an upper limit corresponding with a correlation
    ! between x and y of max_mag_correlation.
    if ( xpyp >  max_mag_correlation * sqrt( xp2 * yp2 ) ) then

        xpyp_chnge =  max_mag_correlation * sqrt( xp2 * yp2 ) - xpyp

        xpyp =  max_mag_correlation * sqrt( xp2 * yp2 )

    ! Clipping for xpyp at a lower limit corresponding with a correlation
    ! between x and y of -max_mag_correlation.
    elseif ( xpyp < -max_mag_correlation * sqrt( xp2 * yp2 ) ) then

        xpyp_chnge = -max_mag_correlation * sqrt( xp2 * yp2 ) - xpyp

        xpyp = -max_mag_correlation * sqrt( xp2 * yp2 )

    else

        xpyp_chnge = zero

    endif

    if ( l_stats_samp ) then
       if ( l_last_clip_ts ) then
          call stat_end_update_pt( ixpyp_cl, level, &
                                   xpyp / dt, stats_zm )
       else
          call stat_modify_pt( ixpyp_cl, level, &
                               xpyp / dt, stats_zm )
       endif
    endif


    return
  end subroutine clip_covar_level

  !=============================================================================
  subroutine clip_variance( solve_type, dt, threshold, &
                            xp2 )

    ! Description:
    ! Clipping the value of variance x'^2 based on a minimum threshold value.
    ! The threshold value must be greater than or equal to 0.
    !
    ! The values of x'^2 are found on the momentum levels.
    !
    ! The following variances are found in the code:
    !
    ! r_t'^2, th_l'^2, u'^2, v'^2, sclr'^2, (computed in advance_xp2_xpyp);
    ! w'^2 (computed in advance_wp2_wp3).

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use clubb_precision, only: & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_begin_update,  & ! Procedure(s)
        stat_end_update

    use stats_variables, only: & 
        stats_zm,  & ! Variable(s)
        iwp2_cl, & 
        irtp2_cl, & 
        ithlp2_cl, & 
        iup2_cl, & 
        ivp2_cl, & 
        l_stats_samp

    implicit none

    ! Input Variables
    integer, intent(in) :: & 
      solve_type  ! Variable being solved; used for STATS.

    real( kind = core_rknd ), intent(in) :: & 
      dt          ! Model timestep; used here for STATS     [s]

    real( kind = core_rknd ), intent(in) :: & 
      threshold   ! Minimum value of x'^2                   [{x units}^2]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      xp2         ! Variance of x, x'^2 (momentum levels)   [{x units}^2]

    ! Local Variables
    integer :: k   ! Array index


    integer :: & 
      ixp2_cl

    ! ---- Begin Code ----

    select case ( solve_type )
    case ( clip_wp2 )   ! wp2 clipping budget term
      ixp2_cl = iwp2_cl
    case ( clip_rtp2 )   ! rtp2 clipping budget term
      ixp2_cl = irtp2_cl
    case ( clip_thlp2 )   ! thlp2 clipping budget term
      ixp2_cl = ithlp2_cl
    case ( clip_up2 )   ! up2 clipping budget term
      ixp2_cl = iup2_cl
    case ( clip_vp2 )   ! vp2 clipping budget term
      ixp2_cl = ivp2_cl
    case default   ! scalars are involved
      ixp2_cl = 0
    end select


    if ( l_stats_samp ) then
      call stat_begin_update( ixp2_cl, xp2 / dt, stats_zm )
    endif

    ! Limit the value of x'^2 at threshold.
    ! The value of x'^2 at the surface (or lower boundary) is a set value that
    ! is determined elsewhere in a surface subroutine.  Thus, the variance
    ! clipping code does not need to be invoked at the lower boundary.
    ! Likewise, the value of x'^2 is set at the upper boundary, so the variance
    ! clipping code does not need to be invoked at the upper boundary.
    !
    ! charlass on 09/11/2013: I changed the clipping so that also the surface
    ! level is clipped. I did this because we discovered that there are slightly
    ! negative values in thlp2(1) and rtp2(1) when running quarter_ss case with
    ! WRF-CLUBB (see wrf:ticket:51#comment:33) 
    do k = 1, gr%nz-1, 1
      if ( xp2(k) < threshold ) then
        xp2(k) = threshold
      endif
    enddo

    if ( l_stats_samp ) then
      call stat_end_update( ixp2_cl, xp2 / dt, stats_zm )
    endif


    return
  end subroutine clip_variance

  !=============================================================================
  subroutine clip_skewness( dt, sfc_elevation, wp2_zt, wp3 )

    ! Description:
    ! Clipping the value of w'^3 based on the skewness of w, Sk_w.
    !
    ! Aditionally, to prevent possible crashes due to wp3 growing too large, 
    ! abs(wp3) will be clipped to 100.
    !
    ! The skewness of w is:
    !
    ! Sk_w = w'^3 / (w'^2)^(3/2).
    !
    ! The value of Sk_w is limited to a range between an upper limit and a lower
    ! limit.  The values of the limits depend on whether the level altitude is
    ! within 100 meters of the surface.
    !
    ! For altitudes less than or equal to 100 meters above ground level (AGL):
    !
    ! -0.2_core_rknd*sqrt(2) <= Sk_w <= 0.2_core_rknd*sqrt(2);
    !
    ! while for all altitudes greater than 100 meters AGL:
    !
    ! -4.5_core_rknd <= Sk_w <= 4.5_core_rknd.
    !
    ! Therefore, there is an upper limit on w'^3, such that:
    !
    ! w'^3  <=  threshold_magnitude * (w'^2)^(3/2);
    !
    ! and a lower limit on w'^3, such that:
    !
    ! w'^3  >= -threshold_magnitude * (w'^2)^(3/2).
    !
    ! The values of w'^3 are found on the thermodynamic levels, while the values
    ! of w'^2 are found on the momentum levels.  Therefore, the values of w'^2
    ! are interpolated to the thermodynamic levels before being used to
    ! calculate the upper and lower limits for w'^3.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
      gr ! Variable(s)

    use clubb_precision, only: & 
      core_rknd ! Variable(s)

    use stats_type_utilities, only: &
      stat_begin_update,  & ! Procedure(s)
      stat_end_update

    use stats_variables, only: & 
      stats_zt,  & ! Variable(s)
      iwp3_cl, & 
      l_stats_samp     

    implicit none

    ! External
    intrinsic :: sign, sqrt, real

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      dt               ! Model timestep; used here for STATS        [s]

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation    ! Elevation of ground level                  [m AMSL]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wp2_zt           ! w'^2 interpolated to thermodyamic levels   [m^2/s^2]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      wp3              ! w'^3 (thermodynamic levels)                [m^3/s^3]

    ! ---- Begin Code ----

    if ( l_stats_samp ) then
      call stat_begin_update( iwp3_cl, wp3 / dt, stats_zt )
    endif

    call clip_skewness_core( sfc_elevation, wp2_zt, wp3 )

    if ( l_stats_samp ) then
      call stat_end_update( iwp3_cl, wp3 / dt, stats_zt )
    endif

    return
  end subroutine clip_skewness

!=============================================================================
  subroutine clip_skewness_core( sfc_elevation, wp2_zt, wp3 )
!
    use grid_class, only: & 
      gr ! Variable(s)

    use constants_clubb, only: &
      Skw_max_mag_sqd ! [-]

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sign, sqrt, real

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation    ! Elevation of ground level                  [m AMSL]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wp2_zt           ! w'^2 interpolated to thermodyamic levels   [m^2/s^2]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      wp3              ! w'^3 (thermodynamic levels)                [m^3/s^3]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      wp2_zt_cubed, & ! Variance of vertical velocity cubed (w^2_{zt}^3)   [m^6/s^6]
      wp3_lim_sqd     ! Keeps absolute value of Sk_w from becoming > limit [m^6/s^6]

    integer :: k       ! Vertical array index.

    real( kind = core_rknd ), parameter :: &  
      wp3_max = 100._core_rknd ! Threshold for wp3 [m^3/s^3]      

    ! ---- Begin Code ----

    ! Compute the upper and lower limits of w'^3 at every level,
    ! based on the skewness of w, Sk_w, such that:
    ! Sk_w = w'^3 / (w'^2)^(3/2);
    ! -4.5 <= Sk_w <= 4.5;
    ! or, if the level altitude is within 100 meters of the surface,
    ! -0.2*sqrt(2) <= Sk_w <= 0.2*sqrt(2).

    ! The normal magnitude limit of skewness of w in the CLUBB code is 4.5.
    ! However, according to Andre et al. (1976b & 1978), wp3 should not exceed
    ! [2*(wp2^3)]^(1/2) at any level.  However, this term should be multiplied
    ! by 0.2 close to the surface to include surface effects.  There already is
    ! a wp3 clipping term in place for all other altitudes, but this term will
    ! be included for the surface layer only.  Therefore, the lowest level wp3
    ! should not exceed 0.2 * sqrt(2) * wp2^(3/2).  Brian Griffin.  12/18/05.

    ! To lower compute time, we squared both sides of the equation and compute
    ! wp2^3 only once. -dschanen 9 Oct 2008

    wp2_zt_cubed(1:gr%nz) = wp2_zt(1:gr%nz)**3

    do k = 1, gr%nz, 1
      if ( gr%zt(k) - sfc_elevation <= 100.0_core_rknd ) then ! Clip for 100 m. AGL.
       !wp3_upper_lim(k) =  0.2_core_rknd * sqrt_2 * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
       !wp3_lower_lim(k) = -0.2_core_rknd * sqrt_2 * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
        wp3_lim_sqd(k) = 0.08_core_rknd * wp2_zt_cubed(k) ! Where 0.08_core_rknd
                              ! == (sqrt(2)*0.2_core_rknd)**2 known magic number
      else                          ! Clip skewness consistently with a.
       !wp3_upper_lim(k) =  4.5_core_rknd * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
       !wp3_lower_lim(k) = -4.5_core_rknd * wp2_zt(k)**(3.0_core_rknd/2.0_core_rknd)
        wp3_lim_sqd(k) = Skw_max_mag_sqd * wp2_zt_cubed(k) ! Skw_max_mag = 4.5_core_rknd^2
      endif
    enddo

    ! Clipping for w'^3 at an upper and lower limit corresponding with
    ! the appropriate value of Sk_w.
    where ( wp3**2 > wp3_lim_sqd ) &
      ! Set the magnitude to the wp3 limit and apply the sign of the current wp3
      wp3 = sign( sqrt( wp3_lim_sqd ), wp3 )

    ! Clipping abs(wp3) to 100. This keeps wp3 from growing too large in some 
    ! deep convective cases, which helps prevent these cases from blowing up.
    where ( abs(wp3) > wp3_max ) &
      wp3 = sign( wp3_max , wp3 ) ! Known magic number

  end subroutine clip_skewness_core

!===============================================================================

end module clip_explicit
