!------------------------------------------------------------------------
! $Id: numerical_check.F90 7309 2014-09-20 17:06:28Z betlej@uwm.edu $
!===============================================================================
module numerical_check

  implicit none

!       Made is_nan_2d public so it may be used
!       for finding code that cause NaNs
!       Joshua Fasching November 2007

!       *_check subroutines were added to ensure that the
!       subroutines they are checking perform correctly
!       Joshua Fasching February 2008

!       rad_clipping has been replaced by rad_check as the new
!       subroutine only reports if there are invalid values.
!       Joshua Fasching March 2008

  private ! Default scope

  public :: invalid_model_arrays, is_nan_2d,  & 
            rad_check, parameterization_check, & 
            surface_varnce_check, pdf_closure_check, & 
            length_check, is_nan_sclr, calculate_spurious_source

  private :: check_negative, check_nan


  ! Abstraction of check_nan
  interface check_nan
    module procedure check_nan_sclr, check_nan_2d
  end interface

  ! Abstraction of check_negative
  interface check_negative
    module procedure check_negative_total, check_negative_index
  end interface


  contains
!---------------------------------------------------------------------------------
  subroutine length_check( Lscale, Lscale_up, Lscale_down, err_code )
!
!        Description: This subroutine determines if any of the output
!        variables for the length_new subroutine carry values that
!        are NaNs.
!
!        Joshua Fasching February 2008
!---------------------------------------------------------------------------------
    use grid_class, only: & 
        gr ! Variable

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters
    character(*), parameter :: proc_name = "compute_length"

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      Lscale,     & ! Mixing length                 [m]
      Lscale_up,  & ! Upward mixing length          [m]
      Lscale_down   ! Downward mixing length        [m]

    ! Output Variable
    integer, intent(inout) :: & 
      err_code

!-----------------------------------------------------------------------------

    call check_nan( Lscale, "Lscale", proc_name, err_code )
    call check_nan( Lscale_up, "Lscale_up", proc_name, err_code )
    call check_nan( Lscale_down, "Lscale_down", proc_name, err_code )

    return
  end subroutine length_check

!---------------------------------------------------------------------------
  subroutine pdf_closure_check( wp4, wprtp2, wp2rtp, wpthlp2, & 
                          wp2thlp, cloud_frac, rcm, wpthvp, wp2thvp, & 
                          rtpthvp, thlpthvp, wprcp, wp2rcp, & 
                          rtprcp, thlprcp, rcp2, wprtpthlp, & 
                          crt_1, crt_2, cthl_1, cthl_2, pdf_params, &
                          sclrpthvp, sclrprcp, wpsclrp2, & 
                          wpsclrprtp, wpsclrpthlp, wp2sclrp, &
                          err_code )

! Description: This subroutine determines if any of the output
!   variables for the pdf_closure subroutine carry values that
!   are NaNs.
!
! Joshua Fasching February 2008
!---------------------------------------------------------------------------

    use parameters_model, only: & 
      sclr_dim ! Variable

    use pdf_parameter_module, only:  &
      pdf_parameter  ! type

    use stats_variables, only: &
      iwp4,       & ! Variables
      ircp2,      &
      iwprtp2,    &
      iwprtpthlp, &
      iwpthlp2

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Parameter Constants
    character(len=*), parameter :: proc_name = &
      "pdf_closure"

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      wp4,             & ! w'^4                  [m^4/s^4]
      wprtp2,          & ! w' r_t'               [(m kg)/(s kg)]
      wp2rtp,          & ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
      wpthlp2,         & ! w' th_l'^2            [(m K^2)/s]
      wp2thlp,         & ! w'^2 th_l'            [(m^2 K)/s^2]
      cloud_frac,      & ! Cloud fraction        [-]
      rcm,             & ! Mean liquid water     [kg/kg]
      wpthvp,          & ! Buoyancy flux         [(K m)/s] 
      wp2thvp,         & ! w'^2 th_v'            [(m^2 K)/s^2]
      rtpthvp,         & ! r_t' th_v'            [(kg K)/kg]
      thlpthvp,        & ! th_l' th_v'           [K^2]
      wprcp,           & ! w' r_c'               [(m kg)/(s kg)]
      wp2rcp,          & ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
      rtprcp,          & ! r_t' r_c'             [(kg^2)/(kg^2)]
      thlprcp,         & ! th_l' r_c'            [(K kg)/kg]
      rcp2,            & ! r_c'^2                [(kg^2)/(kg^2)]
      wprtpthlp,       & ! w' r_t' th_l'         [(m kg K)/(s kg)]
      crt_1, crt_2,  & 
      cthl_1, cthl_2

    type(pdf_parameter), intent(in) ::  & 
      pdf_params        ! PDF parameters          [units vary]

    ! Input (Optional passive scalar variables)
    real( kind = core_rknd ), dimension(sclr_dim), intent(in) ::  & 
      sclrpthvp,  & 
      sclrprcp,  & 
      wpsclrp2, & 
      wpsclrprtp, & 
      wpsclrpthlp, & 
      wp2sclrp

    ! Output Variable
    integer, intent(inout) ::  & 
      err_code          ! Returns appropriate error code

!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    if ( iwp4 > 0 ) call check_nan( wp4,"wp4", proc_name, err_code )
    if ( iwprtp2 > 0 ) call check_nan( wprtp2,"wprtp2", proc_name, err_code )
    call check_nan( wp2rtp,"wp2rtp", proc_name, err_code )
    if ( iwpthlp2 > 0 ) call check_nan( wpthlp2,"wpthlp2", proc_name, err_code )
    call check_nan( wp2thlp,"wp2thlp", proc_name, err_code )
    call check_nan( cloud_frac,"cloud_frac", proc_name, err_code )
    call check_nan( rcm,"rcm", proc_name, err_code )
    call check_nan( wpthvp, "wpthvp", proc_name, err_code )
    call check_nan( wp2thvp, "wp2thvp", proc_name, err_code )
    call check_nan( rtpthvp, "rtpthvp", proc_name, err_code )
    call check_nan( thlpthvp, "thlpthvp", proc_name, err_code )
    call check_nan( wprcp, "wprcp", proc_name, err_code )
    call check_nan( wp2rcp, "wp2rcp", proc_name, err_code )
    call check_nan( rtprcp, "rtprcp", proc_name, err_code )
    call check_nan( thlprcp, "thlprcp", proc_name, err_code )
    if ( ircp2 >  0 ) call check_nan( rcp2, "rcp2", proc_name, err_code)
    if ( iwprtpthlp > 0 ) call check_nan( wprtpthlp, "wprtpthlp", proc_name, err_code )
    call check_nan( crt_1, "crt_1", proc_name, err_code )
    call check_nan( crt_2, "crt_2", proc_name, err_code )
    call check_nan( cthl_1, "cthl_1", proc_name, err_code )
    call check_nan( cthl_2, "cthl_2", proc_name, err_code )
    ! Check each PDF parameter at the grid level sent in.
    call check_nan( pdf_params%w_1, "pdf_params%w_1", proc_name, err_code )
    call check_nan( pdf_params%w_2, "pdf_params%w_2", proc_name, err_code )
    call check_nan( pdf_params%varnce_w_1, "pdf_params%varnce_w_1", proc_name, err_code )
    call check_nan( pdf_params%varnce_w_2, "pdf_params%varnce_w_2", proc_name, err_code )
    call check_nan( pdf_params%rt_1, "pdf_params%rt_1", proc_name, err_code )
    call check_nan( pdf_params%rt_2, "pdf_params%rt_2", proc_name, err_code )
    call check_nan( pdf_params%varnce_rt_1, "pdf_params%varnce_rt_1", proc_name, err_code )
    call check_nan( pdf_params%varnce_rt_2, "pdf_params%varnce_rt_2", proc_name, err_code )
    call check_nan( pdf_params%thl_1, "pdf_params%thl_1", proc_name, err_code )
    call check_nan( pdf_params%thl_2, "pdf_params%thl_2", proc_name, err_code )
    call check_nan( pdf_params%varnce_thl_1, "pdf_params%varnce_thl_1", proc_name, err_code )
    call check_nan( pdf_params%varnce_thl_2, "pdf_params%varnce_thl_2", proc_name, err_code )
    call check_nan( pdf_params%mixt_frac, "pdf_params%mixt_frac", proc_name, err_code )
    call check_nan( pdf_params%rrtthl, "pdf_params%rrtthl", proc_name, err_code )
    call check_nan( pdf_params%rc_1, "pdf_params%rc_1", proc_name, err_code )
    call check_nan( pdf_params%rc_2, "pdf_params%rc_2", proc_name, err_code )
    call check_nan( pdf_params%rsatl_1, "pdf_params%rsatl_1", proc_name, err_code )
    call check_nan( pdf_params%rsatl_2, "pdf_params%rsatl_2", proc_name, err_code )
    call check_nan( pdf_params%cloud_frac_1, "pdf_params%cloud_frac_1", proc_name, err_code )
    call check_nan( pdf_params%cloud_frac_2, "pdf_params%cloud_frac_2", proc_name, err_code )
    call check_nan( pdf_params%chi_1, "pdf_params%chi_1", proc_name, err_code )
    call check_nan( pdf_params%chi_2, "pdf_params%chi_2", proc_name, err_code )
    call check_nan( pdf_params%stdev_chi_1, "pdf_params%stdev_chi_1", proc_name, err_code )
    call check_nan( pdf_params%stdev_chi_2, "pdf_params%stdev_chi_2", proc_name, err_code )
    call check_nan( pdf_params%alpha_thl, "pdf_params%alpha_thl", proc_name, err_code )
    call check_nan( pdf_params%alpha_rt, "pdf_params%alpha_rt", proc_name, err_code )

    if ( sclr_dim > 0 ) then
      call check_nan( sclrpthvp,"sclrpthvp", & 
                      proc_name, err_code)
      call check_nan( sclrprcp, "sclrprcp", & 
                      proc_name, err_code )
      call check_nan( wpsclrprtp, "wpsclrprtp",  & 
                      proc_name, err_code )
      call check_nan( wpsclrp2, "wpsclrp2",  & 
                      proc_name, err_code )
      call check_nan( wpsclrpthlp, "wpsclrtlp",  & 
                      proc_name, err_code )
      call check_nan( wp2sclrp, "wp2sclrp",  & 
                      proc_name, err_code )
    end if

    return
  end subroutine pdf_closure_check

!-------------------------------------------------------------------------------
  subroutine parameterization_check & 
             ( thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
               wm_zm, wm_zt, p_in_Pa, rho_zm, rho, exner, &
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
               invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
               um, upwp, vm, vpwp, up2, vp2, &
               rtm, wprtp, thlm, wpthlp, &
               wp2, wp3, rtp2, thlp2, rtpthlp, &
               prefix, &
               wpsclrp_sfc, wpedsclrp_sfc, & 
               sclrm, wpsclrp, sclrp2, sclrprtp, sclrpthlp, &
               sclrm_forcing, edsclrm, edsclrm_forcing, err_code )
!
! Description:
!   This subroutine determines what input variables may have NaN values.
!   In addition it checks to see if rho_zm, rho, exner, up2, vp2, rtm, thlm,
!   wp2, rtp2, thlp2, or tau_zm have negative values.
!-------------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable

    use parameters_model, only: & 
        sclr_dim,  & ! Variable
        edsclr_dim

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters
    ! Name of the procedure using parameterization_check
    character(len=25), parameter ::  & 
      proc_name = "parameterization_timestep"

    ! Input variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      p_in_Pa,         & ! Air pressure (thermodynamic levels)       [Pa]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      exner,           & ! Exner function (thermodynamic levels)     [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo. levs.  [K]

    real( kind = core_rknd ), intent(in) ::  & 
      wpthlp_sfc,   & ! w' theta_l' at surface.   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface.       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface.          [m^2/s^2]
      vpwp_sfc        ! v'w' at surface.          [m^2/s^2]

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    character(len=*), intent(in) :: prefix ! Location where subroutine is called

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: & 
      wpsclrp_sfc    ! Scalar flux at surface [units m/s]

    real( kind = core_rknd ), intent(in), dimension(edsclr_dim) :: & 
      wpedsclrp_sfc ! Eddy-Scalar flux at surface      [units m/s]

    real( kind = core_rknd ), intent(in),dimension(gr%nz,sclr_dim) :: & 
      sclrm,         & ! Passive scalar mean      [units vary]
      wpsclrp,       & ! w'sclr'                  [units vary]
      sclrp2,        & ! sclr'^2                  [units vary]
      sclrprtp,      & ! sclr'rt'                 [units vary]
      sclrpthlp,     & ! sclr'thl'                [units vary]
      sclrm_forcing    ! Passive scalar forcing   [units / s]

    real( kind = core_rknd ), intent(in),dimension(gr%nz,edsclr_dim) :: & 
      edsclrm,         & ! Eddy passive scalar mean    [units vary]
      edsclrm_forcing    ! Eddy passive scalar forcing [units / s]

    ! In / Out Variables
    integer, intent(inout) :: &
      err_code ! Error code

    ! Local Variables
    integer :: i ! Loop iterator for the scalars

!-------- Input Nan Check ----------------------------------------------

    call check_nan( thlm_forcing, "thlm_forcing", prefix//proc_name, err_code)
    call check_nan( rtm_forcing,"rtm_forcing", prefix//proc_name, err_code )
    call check_nan( um_forcing,"um_forcing", prefix//proc_name, err_code )
    call check_nan( vm_forcing,"vm_forcing", prefix//proc_name, err_code )

    call check_nan( wm_zm, "wm_zm", prefix//proc_name, err_code )
    call check_nan( wm_zt, "wm_zt", prefix//proc_name, err_code )
    call check_nan( p_in_Pa, "p_in_Pa", prefix//proc_name, err_code )
    call check_nan( rho_zm, "rho_zm", prefix//proc_name, err_code )
    call check_nan( rho, "rho", prefix//proc_name, err_code )
    call check_nan( exner, "exner", prefix//proc_name, err_code )
    call check_nan( rho_ds_zm, "rho_ds_zm", prefix//proc_name, err_code )
    call check_nan( rho_ds_zt, "rho_ds_zt", prefix//proc_name, err_code )
    call check_nan( invrs_rho_ds_zm, "invrs_rho_ds_zm", prefix//proc_name, err_code )
    call check_nan( invrs_rho_ds_zt, "invrs_rho_ds_zt", prefix//proc_name, err_code )
    call check_nan( thv_ds_zm, "thv_ds_zm", prefix//proc_name, err_code )
    call check_nan( thv_ds_zt, "thv_ds_zt", prefix//proc_name, err_code )

    call check_nan( um, "um", prefix//proc_name, err_code )
    call check_nan( upwp, "upwp", prefix//proc_name, err_code )
    call check_nan( vm, "vm", prefix//proc_name, err_code )
    call check_nan( vpwp, "vpwp", prefix//proc_name, err_code )
    call check_nan( up2, "up2", prefix//proc_name, err_code )
    call check_nan( vp2, "vp2", prefix//proc_name, err_code )
    call check_nan( rtm, "rtm", prefix//proc_name, err_code )
    call check_nan( wprtp, "wprtp", prefix//proc_name, err_code )
    call check_nan( thlm, "thlm", prefix//proc_name, err_code )
    call check_nan( wpthlp, "wpthlp", prefix//proc_name, err_code )
    call check_nan( wp2, "wp2", prefix//proc_name, err_code )
    call check_nan( wp3, "wp3", prefix//proc_name, err_code )
    call check_nan( rtp2, "rtp2", prefix//proc_name, err_code )
    call check_nan( thlp2, "thlp2", prefix//proc_name, err_code )
    call check_nan( rtpthlp, "rtpthlp", prefix//proc_name, err_code )

    call check_nan( wpthlp_sfc, "wpthlp_sfc", prefix//proc_name, err_code )
    call check_nan( wprtp_sfc, "wprtp_sfc", prefix//proc_name, err_code )
    call check_nan( upwp_sfc, "upwp_sfc", prefix//proc_name, err_code )
    call check_nan( vpwp_sfc, "vpwp_sfc", prefix//proc_name, err_code )

    do i = 1, sclr_dim

      call check_nan( sclrm_forcing(:,i),"sclrm_forcing",  & 
                      prefix//proc_name, err_code )

      call check_nan( wpsclrp_sfc(i),"wpsclrp_sfc",  & 
                      prefix//proc_name, err_code )

      call check_nan( sclrm(:,i),"sclrm", prefix//proc_name, err_code )
      call check_nan( wpsclrp(:,i),"wpsclrp", prefix//proc_name, err_code )
      call check_nan( sclrp2(:,i),"sclrp2", prefix//proc_name, err_code )
      call check_nan( sclrprtp(:,i),"sclrprtp", prefix//proc_name, err_code )
      call check_nan( sclrpthlp(:,i),"sclrpthlp", prefix//proc_name, err_code )

    end do


    do i = 1, edsclr_dim

      call check_nan( edsclrm_forcing(:,i),"edsclrm_forcing", prefix//proc_name, err_code )

      call check_nan( wpedsclrp_sfc(i),"wpedsclrp_sfc",  & 
                      prefix//proc_name, err_code )

      call check_nan( edsclrm(:,i),"edsclrm", prefix//proc_name, err_code )

    enddo

!---------------------------------------------------------------------


    call check_negative( rtm, gr%nz ,"rtm", prefix//proc_name, err_code )
    call check_negative( p_in_Pa, gr%nz ,"p_in_Pa", prefix//proc_name, err_code )
    call check_negative( rho, gr%nz ,"rho", prefix//proc_name, err_code )
    call check_negative( rho_zm, gr%nz ,"rho_zm", prefix//proc_name, err_code )
    call check_negative( exner, gr%nz ,"exner", prefix//proc_name, err_code )
    call check_negative( rho_ds_zm, gr%nz ,"rho_ds_zm", prefix//proc_name, err_code )
    call check_negative( rho_ds_zt, gr%nz ,"rho_ds_zt", prefix//proc_name, err_code )
    call check_negative( invrs_rho_ds_zm, gr%nz ,"invrs_rho_ds_zm", &
                         prefix//proc_name, err_code )
    call check_negative( invrs_rho_ds_zt, gr%nz ,"invrs_rho_ds_zt", &
                         prefix//proc_name, err_code )
    call check_negative( thv_ds_zm, gr%nz ,"thv_ds_zm", prefix//proc_name, err_code )
    call check_negative( thv_ds_zt, gr%nz ,"thv_ds_zt", prefix//proc_name, err_code )
    call check_negative( up2, gr%nz ,"up2", prefix//proc_name, err_code )
    call check_negative( vp2, gr%nz ,"vp2", prefix//proc_name, err_code )
    call check_negative( wp2, gr%nz ,"wp2", prefix//proc_name, err_code )
    call check_negative( rtm, gr%nz ,"rtm", prefix//proc_name, err_code )
    call check_negative( thlm, gr%nz ,"thlm", prefix//proc_name, err_code )
    call check_negative( rtp2, gr%nz ,"rtp2", prefix//proc_name, err_code )
    call check_negative( thlp2, gr%nz ,"thlp2", prefix//proc_name, err_code )

    return
  end subroutine parameterization_check

!-----------------------------------------------------------------------
  subroutine surface_varnce_check( wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, & 
           rtp2_sfc, rtpthlp_sfc, &
           sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc, &
           err_code )
!
!       Description:This subroutine determines if any of the output
!       variables for the surface_varnce subroutine carry values that
!       are nans.
!
!       Joshua Fasching February 2008
!
!
!-----------------------------------------------------------------------
    use parameters_model, only: & 
        sclr_dim ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters
    ! Name of the subroutine calling the check
    character(len=*), parameter :: &
      proc_name = "surface_varnce"

    ! Input Variables
    real( kind = core_rknd ),intent(in) ::  & 
      wp2_sfc,     & ! Vertical velocity variance        [m^2/s^2]
      up2_sfc,     & ! u'^2                              [m^2/s^2]
      vp2_sfc,     & ! u'^2                              [m^2/s^2]
      thlp2_sfc,   & ! thetal variance                   [K^2]
      rtp2_sfc,    & ! rt variance                       [(kg/kg)^2]
      rtpthlp_sfc    ! thetal rt covariance              [kg K/kg]


    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: & 
      sclrp2_sfc,    & ! Passive scalar variance                 [units^2]
      sclrprtp_sfc,  & ! Passive scalar r_t covariance           [units kg/kg]
      sclrpthlp_sfc ! Passive scalar theta_l covariance       [units K]

    ! Input/Output Variable
    integer, intent(inout) :: err_code    ! Are these outputs valid?

!-----------------------------------------------------------------------

    ! ---- Begin Code ----

    call check_nan( wp2_sfc, "wp2_sfc", proc_name, err_code)
    call check_nan( up2_sfc, "up2_sfc", proc_name, err_code)
    call check_nan( vp2_sfc, "vp2_sfc", proc_name, err_code)
    call check_nan( thlp2_sfc, "thlp2_sfc", proc_name, err_code)
    call check_nan( rtp2_sfc, "rtp2_sfc", proc_name, err_code)
    call check_nan( rtpthlp_sfc, "rtpthlp_sfc",  & 
                    proc_name, err_code)

    if ( sclr_dim > 0 ) then
      call check_nan( sclrp2_sfc, "sclrp2_sfc", & 
                      proc_name, err_code )

      call check_nan( sclrprtp_sfc, "sclrprtp_sfc", & 
                      proc_name, err_code )

      call check_nan( sclrpthlp_sfc, "sclrpthlp_sfc",  & 
                      proc_name, err_code )
    end if

    return
  end subroutine surface_varnce_check

!-----------------------------------------------------------------------
  subroutine rad_check( thlm, rcm, rtm, rim,  & 
                        cloud_frac, p_in_Pa, exner, rho_zm )
! Description:
!   Checks radiation input variables. If they are < 0 it reports
!   to the console.
!------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters
    character(len=*), parameter ::  & 
      proc_name = "Before BUGSrad."

    ! Input/Output variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      thlm,           & ! Liquid Water Potential Temperature   [K/s]
      rcm,            & ! Liquid Water Mixing Ratio            [kg/kg]
      rtm,            & ! Total Water Mixing Ratio             [kg/kg]
      rim,          & ! Ice Water Mixing Ratio               [kg/kg]
      cloud_frac,     & ! Cloud Fraction                       [-]
      p_in_Pa,        & ! Pressure                             [Pa]
      exner,          & ! Exner Function                       [-]
      rho_zm            ! Air Density                          [kg/m^3]

    ! Local variables
    real( kind = core_rknd ),dimension(gr%nz) :: rvm

!-------------------------------------------------------------------------

    rvm = rtm - rcm

    call check_negative( thlm, gr%nz ,"thlm", proc_name )
    call check_negative( rcm, gr%nz ,"rcm", proc_name )
    call check_negative( rtm, gr%nz ,"rtm", proc_name )
    call check_negative( rvm, gr%nz ,"rvm", proc_name )
    call check_negative( rim, gr%nz ,"rim", proc_name )
    call check_negative( cloud_frac, gr%nz ,"cloud_frac", proc_name )
    call check_negative( p_in_Pa, gr%nz ,"p_in_Pa", proc_name )
    call check_negative( exner, gr%nz ,"exner", proc_name )
    call check_negative( rho_zm, gr%nz ,"rho_zm", proc_name )

    return

  end subroutine rad_check

!-----------------------------------------------------------------------
  logical function invalid_model_arrays( )

!       Description:
!       Checks for invalid floating point values in select model arrays.

!       References:
!       None
!------------------------------------------------------------------------

    use variables_diagnostic_module, only: & 
      hydromet,  & ! Variable(s)
      wp2thvp, & 
      rtpthvp, & 
      thlpthvp

    use variables_prognostic_module, only: & 
      um,  & ! Variable(s)
      vm, & 
      wp2, & 
      wp3, & 
      rtm, & 
      thlm, & 
      rtp2, & 
      thlp2, & 
      wprtp, & 
      wpthlp, & 
      rtpthlp, & 
      sclrm, & 
      edsclrm

    use constants_clubb, only: & 
      fstderr   ! Constant(s)

    use parameters_model, only: & 
      sclr_dim,  & ! Variable(s)
      edsclr_dim, &
      hydromet_dim

    use array_index, only: &
      hydromet_list ! Variable(s)

    implicit none

    ! Local Variables
    integer :: i

    invalid_model_arrays = .false.

    ! Check whether any variable array contains a NaN for
    ! um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp,
    ! wp2, & wp3.
    if ( is_nan_2d( um ) ) then
      write(fstderr,*) "NaN in um model array"
!         write(fstderr,*) "um= ", um
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( vm ) ) then
      write(fstderr,*) "NaN in vm model array"
!         write(fstderr,*) "vm= ", vm
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wp2 ) ) then
      write(fstderr,*) "NaN in wp2 model array"
!         write(fstderr,*) "wp2= ", wp2
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wp3 ) ) then
      write(fstderr,*) "NaN in wp3 model array"
!         write(fstderr,*) "wp3= ", wp3
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtm ) ) then
      write(fstderr,*) "NaN in rtm model array"
!         write(fstderr,*) "rtm= ", rtm
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( thlm ) ) then
      write(fstderr,*) "NaN in thlm model array"
!         write(fstderr,*) "thlm= ", thlm
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtp2 ) ) then
      write(fstderr,*) "NaN in rtp2 model array"
!         write(fstderr,*) "rtp2= ", rtp2
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( thlp2 ) ) then
      write(fstderr,*) "NaN in thlp2 model array"
!         write(fstderr,*) "thlp2= ", thlp2
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wprtp ) ) then
      write(fstderr,*) "NaN in wprtp model array"
!         write(fstderr,*) "wprtp= ", wprtp
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( wpthlp ) ) then
      write(fstderr,*) "NaN in wpthlp model array"
!         write(fstderr,*) "wpthlp= ", wpthlp
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtpthlp ) ) then
      write(fstderr,*) "NaN in rtpthlp model array"
!         write(fstderr,*) "rtpthlp= ", rtpthlp
      invalid_model_arrays = .true.
!         return
    end if

    if ( hydromet_dim > 0 ) then
      do i = 1, hydromet_dim, 1
        if ( is_nan_2d( hydromet(:,i) ) ) then
          write(fstderr,*) "NaN in a hydrometeor model array "// &
            trim( hydromet_list(i) )
!             write(fstderr,*) "hydromet= ", hydromet
          invalid_model_arrays = .true.
!             return
        end if
      end do
    end if

!       if ( is_nan_2d( wm_zt ) ) then
!         write(fstderr,*) "NaN in wm_zt model array"
!         write(fstderr,*) "wm_zt= ", wm_zt
!         invalid_model_arrays = .true.
!         return
!       end if

    if ( is_nan_2d( wp2thvp ) ) then
      write(fstderr,*) "NaN in wp2thvp model array"
!         write(fstderr,*) "wp2thvp = ", wp2thvp
      invalid_model_arrays = .true.
!         return
    end if

    if ( is_nan_2d( rtpthvp ) ) then
      write(fstderr,*) "NaN in rtpthvp model array"
!         write(fstderr,*) "rtpthvp = ", rtpthvp
      invalid_model_arrays = .true.
    end if

    if ( is_nan_2d( thlpthvp ) ) then
      write(fstderr,*) "NaN in thlpthvp model array"
!         write(fstderr,*) "thlpthvp = ", thlpthvp
      invalid_model_arrays = .true.
    end if

    do i = 1, sclr_dim, 1
      if ( is_nan_2d( sclrm(:,i) ) ) then
        write(fstderr,*) "NaN in sclrm", i, "model array"
!           write(fstderr,'(a6,i2,a1)') "sclrm(", i, ")"
!           write(fstderr,*) sclrm(:,i)
        invalid_model_arrays = .true.
      end if
    end do

    do i = 1, edsclr_dim, 1
      if ( is_nan_2d( edsclrm(:,i) ) ) then
        write(fstderr,*) "NaN in edsclrm", i, "model array"
!           write(fstderr,'(a8,i2,a1)') "edsclrm(", i, ")"
!           write(fstderr,*) edsclrm(:,i)
        invalid_model_arrays = .true.
      end if
    end do

    return
  end function invalid_model_arrays

!------------------------------------------------------------------------
  logical function is_nan_sclr( xarg )

! Description:
!   Checks if a given scalar real is a NaN, +inf or -inf.

! Notes:
!   I was advised by Andy Vaught to use a data statement and the transfer( )
!   intrinsic rather than using a hex number in a parameter for portability.

!   Certain compiler optimizations may cause variables with invalid
!   results to flush to zero.  Avoid these!
!  -dschanen 16 Dec 2010

!------------------------------------------------------------------------

#ifndef __GFORTRAN__
    use parameters_model, only: &
      PosInf ! Variable(s)
#endif

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: xarg

#ifdef __GFORTRAN__  /* if the isnan extension is available, we use it here */
    is_nan_sclr = isnan( xarg )
#else
    ! ---- Begin Code ---

    ! This works on compilers with standardized floating point,
    ! because the IEEE 754 spec defines that subnormals and nans
    ! should not equal themselves.
    ! However, all compilers do not seem to follow this.
    if (xarg /= xarg ) then
      is_nan_sclr = .true.

      ! This a second check, assuming the above does not work as
      ! expected.
    else if ( xarg == PosInf ) then
      is_nan_sclr = .true.

    else
      is_nan_sclr = .false. ! Our result should be a standard float

    end if
#endif

    return
  end function is_nan_sclr
!------------------------------------------------------------------------

!------------------------------------------------------------------------
  logical function is_nan_2d( x2d )

! Description:
!   Checks if a given real vector is a NaN, +inf or -inf.

!------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: any

    ! Input Variables
    real( kind = core_rknd ), dimension(:), intent(in) :: x2d

    ! Local Variables
    integer :: k

    ! ---- Begin Code ----

    is_nan_2d = .false.

    do k = 1, size( x2d )
      if ( is_nan_sclr( x2d(k) ) ) then
        is_nan_2d = .true.
        exit
      end if
    end do

    return

  end function is_nan_2d

!------------------------------------------------------------------------
  subroutine check_negative_total & 
            ( var, varname, operation, err_code )
!
! Description:
!   Checks for negative values in the var array and reports them.
!
!-----------------------------------------------------------------------
    use constants_clubb, only: & 
        fstderr ! Variable(s)

    use error_code, only:  & 
        clubb_var_less_than_zero ! Variable(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: any, present

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(:) :: var

    character(len=*), intent(in)::  & 
    varname,     & ! Varible being examined
    operation   ! Procedure calling check_zero

    ! Optional In/Out Variable
    integer, optional, intent(inout) :: err_code

    if ( any( var < 0.0_core_rknd ) ) then

      write(fstderr,*) varname, " < 0 in ", operation
      if ( present( err_code ) ) then
        if (err_code < clubb_var_less_than_zero ) then
          err_code = clubb_var_less_than_zero
        end if
      end if

    end if ! any ( var < 0 )

    return

  end subroutine check_negative_total


!------------------------------------------------------------------------
  subroutine check_negative_index & 
            ( var, ndim, varname, operation, err_code )
!
! Description:
!   Checks for negative values in the var array and reports
!   the index in which the negative values occur.
!
!-----------------------------------------------------------------------
    use constants_clubb, only: & 
        fstderr ! Variable

    use error_code, only:  & 
        clubb_var_less_than_zero ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: any, present

    ! Input Variables
    integer, intent(in) :: ndim

    real( kind = core_rknd ), intent(in), dimension(ndim) :: var

    character(len=*), intent(in)::  & 
    varname,     & ! Varible being examined
    operation   ! Procedure calling check_zero

    ! Optional In/Out Variable
    integer, optional, intent(inout) :: err_code

    ! Local Variable
    integer :: k ! Loop iterator

    do k=1,ndim,1

      if ( var(k) < 0.0_core_rknd )  then

        write(fstderr,*) varname, " < 0 in ", operation,  & 
                         " at k = ", k

        if ( present( err_code ) ) then
          if (err_code < clubb_var_less_than_zero ) then
            err_code = clubb_var_less_than_zero
          end if
        end if

      end if

    end do ! 1..n

    return

  end subroutine check_negative_index


!------------------------------------------------------------------------
  subroutine check_nan_2d( var, varname, operation, err_code )
!
!  Description:
!    Checks for a NaN in the var array and reports it.
!
!
!------------------------------------------------------------------------
    use constants_clubb, only:  & 
        fstderr ! Variable(s)
    use error_code, only:  & 
        clubb_var_equals_NaN ! Variable(s)
    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: present

    ! Input variables
    real( kind = core_rknd ), intent(in), dimension(:) :: var ! Variable being examined

    character(len=*), intent(in)::  & 
      varname,  & ! Name of variable
      operation   ! Procedure calling check_nan

    ! Optional In/Out Variable
    integer, optional, intent(inout) :: err_code

    if ( is_nan_2d( var ) ) then
      write(fstderr,*) varname, " is NaN in ",operation
      if ( present( err_code ) ) then
        if( err_code < clubb_var_equals_NaN ) then
          err_code = clubb_var_equals_NaN
        end if
      end if
    end if

    return
  end subroutine check_nan_2d

!-----------------------------------------------------------------------
  subroutine check_nan_sclr( var, varname, operation, err_code )
!
! Description:
!   Checks for a NaN in the scalar var then reports it.
!
!-----------------------------------------------------------------------
    use constants_clubb, only:  & 
        fstderr ! Variable
    use error_code, only:  & 
        clubb_var_equals_NaN ! Variable
    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: present

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: var        ! Variable being examined

    character(len=*), intent(in)::  & 
      varname,    & ! Name of variable being examined
      operation  ! Procedure calling check_nan

    ! Optional In/Out variable
    integer, optional, intent(inout) :: err_code
!--------------------------------------------------------------------
    if ( is_nan_sclr( var ) ) then
      write(fstderr,*) varname, " is NaN in ",operation
      if ( present( err_code ) ) then
        if( err_code < clubb_var_equals_NaN ) then
          err_code = clubb_var_equals_NAN
        end if
      end if
    end if

    return

  end subroutine check_nan_sclr
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
  pure function calculate_spurious_source( integral_after, integral_before, &
                                           flux_top, flux_sfc, & 
                                           integral_forcing, dt ) &
  result( spurious_source )
!
! Description:
!   Checks whether there is conservation within the column and returns any
!   imbalance as spurious_source where spurious_source is defined negative
!   for a spurious sink.
!
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      integral_after, &   ! Vertically-integrated quantity after dt time  [units vary]
      integral_before, &  ! Vertically-integrated quantity before dt time [units vary]
      flux_top, &         ! Total flux at the top of the domain           [units vary]
      flux_sfc, &         ! Total flux at the bottom of the domain        [units vary]
      integral_forcing, & ! Vertically-integrated forcing                 [units vary] 
      dt                  ! Timestep size                                 [s]

    ! Return Variable
    real( kind = core_rknd ) :: spurious_source ! [units vary]

!--------------------------------------------------------------------

    ! ---- Begin Code ----

    spurious_source = (integral_after - integral_before) / dt & 
                        + flux_top - flux_sfc - integral_forcing

    return

  end function calculate_spurious_source
!-------------------------------------------------------------------------
end module numerical_check
