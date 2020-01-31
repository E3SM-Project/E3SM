!-----------------------------------------------------------------------
! $Id: stats_LH_zt.F90 5997 2012-12-18 20:47:09Z raut@uwm.edu $

module stats_LH_zt

  implicit none

  private ! Default Scope

  public :: stats_init_LH_zt

! Constant parameters
  integer, parameter, public :: nvarmax_LH_zt = 100 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_LH_zt( vars_LH_zt, l_error )

! Description:
!   Initializes array indices for zt

! Note:
!   All code that is within subroutine stats_init_zt, including variable
!   allocation code, is not called if l_stats is false.  This subroutine is
!   called only when l_stats is true.

!-----------------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr ! Constant(s)

    use stats_variables, only: & 
      LH_zt    ! Variable

    use stats_variables, only: & 
      iAKm, &  ! Variable(s)
      iLH_AKm, & 
      iAKstd, & 
      iAKstd_cld, & 
      iAKm_rcm, & 
      iAKm_rcc

    use stats_variables, only: &
      iLH_thlm_mc, &  ! Variable(s)
      iLH_rvm_mc, & 
      iLH_rcm_mc, & 
      iLH_Ncm_mc, & 
      iLH_rrainm_mc, & 
      iLH_Nrm_mc, & 
      iLH_rsnowm_mc, & 
      iLH_Nsnowm_mc, & 
      iLH_rgraupelm_mc, & 
      iLH_Ngraupelm_mc, & 
      iLH_ricem_mc, & 
      iLH_Nim_mc, & 
      iLH_Vrr, &
      iLH_VNr, &
      iLH_rcm_avg

    use stats_variables, only: &
      iLH_rrainm, & ! Variable(s)
      iLH_Nrm, &
      iLH_ricem, &
      iLH_Nim, &
      iLH_rsnowm, &
      iLH_Nsnowm, &
      iLH_rgraupelm, &
      iLH_Ngraupelm, &
      iLH_thlm, &
      iLH_rcm, &
      iLH_Ncm, &
      iLH_rvm, &
      iLH_wm, &
      iLH_wp2_zt, &
      iLH_rcp2_zt, &
      iLH_rtp2_zt, &
      iLH_thlp2_zt, &
      iLH_rrainp2_zt, &
      iLH_Nrp2_zt, &
      iLH_Ncp2_zt, &
      iLH_cloud_frac, &
      iLH_rrainm_auto, &
      iLH_rrainm_accr


    use stats_type, only: & 
      stat_assign ! Procedure

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_LH_zt), intent(in) :: vars_LH_zt

    ! Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for LH_zt

    iAKm            = 0  ! analytic Kessler.  Vince Larson 22 May 2005
    iLH_AKm         = 0  ! LH Kessler.  Vince Larson 22 May 2005
    iAKstd          = 0
    iAKstd_cld      = 0
    iAKm_rcm        = 0
    iAKm_rcc        = 0

    iLH_thlm_mc      = 0
    iLH_rvm_mc       = 0
    iLH_rcm_mc       = 0
    iLH_Ncm_mc       = 0
    iLH_rrainm_mc    = 0
    iLH_Nrm_mc       = 0
    iLH_rsnowm_mc    = 0
    iLH_Nsnowm_mc    = 0   
    iLH_rgraupelm_mc = 0 
    iLH_Ngraupelm_mc = 0 
    iLH_ricem_mc     = 0
    iLH_Nim_mc       = 0

    iLH_rcm_avg = 0

    iLH_Vrr = 0
    iLH_VNr = 0

    iLH_rrainm = 0
    iLH_ricem = 0
    iLH_rsnowm = 0
    iLH_rgraupelm = 0

    iLH_Nrm = 0
    iLH_Nim = 0
    iLH_Nsnowm = 0
    iLH_Ngraupelm = 0

    iLH_thlm = 0
    iLH_rcm = 0
    iLH_rvm = 0
    iLH_wm = 0
    iLH_cloud_frac = 0

    iLH_wp2_zt = 0
    iLH_rcp2_zt = 0
    iLH_rtp2_zt = 0
    iLH_thlp2_zt = 0
    iLH_rrainp2_zt = 0
    iLH_Nrp2_zt = 0
    iLH_Ncp2_zt = 0

    iLH_rrainm_auto = 0
    iLH_rrainm_accr = 0

    ! Assign pointers for statistics variables zt

    k = 1
    do i=1,LH_zt%nn

      select case ( trim(vars_LH_zt(i)) )
      case ( 'AKm' )           ! Vince Larson 22 May 2005
        iAKm = k
        call stat_assign( iAKm, "AKm", & 
             "Analytic Kessler ac [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_AKm' )       ! Vince Larson 22 May 2005
        iLH_AKm = k

        call stat_assign( iLH_AKm, "LH_AKm", & 
             "LH Kessler estimate  [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'AKstd' )
        iAKstd = k

        call stat_assign( iAKstd, "AKstd", & 
             "Exact standard deviation of gba Kessler [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'AKstd_cld' )
        iAKstd_cld = k

        call stat_assign( iAKstd_cld, "AKstd_cld", & 
             "Exact w/in cloud std of gba Kessler [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'AKm_rcm' )
        iAKm_rcm = k

        call stat_assign( iAKm_rcm, "AKm_rcm", & 
             "Exact local gba auto based on rcm [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'AKm_rcc' )
        iAKm_rcc = k

        call stat_assign( iAKm_rcc, "AKm_rcc", & 
             "Exact local gba based on w/in cloud rc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_rvm_mc' )
        iLH_rvm_mc = k

        call stat_assign( iLH_rvm_mc, "LH_rvm_mc", & 
             "Latin hypercube estimate of rvm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_thlm_mc' )
        iLH_thlm_mc = k

        call stat_assign( iLH_thlm_mc, "LH_thlm_mc", & 
             "Latin hypercube estimate of thlm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_rcm_mc' )
        iLH_rcm_mc = k

        call stat_assign( iLH_rcm_mc, "LH_rcm_mc", & 
             "Latin hypercube estimate of rcm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_Ncm_mc' )
        iLH_Ncm_mc = k

        call stat_assign( iLH_Ncm_mc, "LH_Ncm_mc", & 
             "Latin hypercube estimate of Ncm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_rrainm_mc' )
        iLH_rrainm_mc = k

        call stat_assign( iLH_rrainm_mc, "LH_rrainm_mc", & 
             "Latin hypercube estimate of rrainm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_Nrm_mc' )
        iLH_Nrm_mc = k

        call stat_assign( iLH_Nrm_mc, "LH_Nrm_mc", & 
             "Latin hypercube estimate of Nrm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case('LH_rsnowm_mc')
        iLH_rsnowm_mc = k

        call stat_assign( iLH_rsnowm_mc, "LH_rsnowm_mc", & 
             "Latin hypercube estimate of rsnowm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_Nsnowm_mc' )
        iLH_Nsnowm_mc = k

        call stat_assign( iLH_Nsnowm_mc, "LH_Nsnowm_mc", & 
             "Latin hypercube estimate of Nsnowm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_rgraupelm_mc' )
        iLH_rgraupelm_mc = k

        call stat_assign( iLH_rgraupelm_mc, "LH_rgraupelm_mc", & 
             "Latin hypercube estimate of rgraupelm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_Ngraupelm_mc' )
        iLH_Ngraupelm_mc = k

        call stat_assign( iLH_Ngraupelm_mc, "LH_Ngraupelm_mc", & 
             "Latin hypercube estimate of Ngraupelm_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_ricem_mc' )
        iLH_ricem_mc = k

        call stat_assign( iLH_ricem_mc, "LH_ricem_mc", & 
             "Latin hypercube estimate of ricem_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_Nim_mc' )
        iLH_Nim_mc = k

        call stat_assign( iLH_Nim_mc, "LH_Nim_mc", & 
             "Latin hypercube estimate of Nim_mc [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_Vrr' )
        iLH_Vrr = k

        call stat_assign( iLH_Vrr, "LH_Vrr", & 
             "Latin hypercube estimate of rrainm sedimentation velocity [m/s]", "m/s", LH_zt )
        k = k + 1

      case ( 'LH_VNr' )
        iLH_VNr = k

        call stat_assign( iLH_VNr, "LH_VNr", & 
             "Latin hypercube estimate of Nrm sedimentation velocity [m/s]", "m/s", LH_zt )
        k = k + 1

      case ( 'LH_rcm_avg' )
        iLH_rcm_avg = k

        call stat_assign( iLH_rcm_avg, "LH_rcm_avg", & 
             "Latin hypercube average estimate of rcm [kg/kg]", "kg/kg", LH_zt )

        k = k + 1

      case ( 'LH_rrainm' )
        iLH_rrainm = k

        call stat_assign( iLH_rrainm, "LH_rrainm", & 
             "Latin hypercube estimate of rrainm [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_Nrm' )
        iLH_Nrm = k

        call stat_assign( iLH_Nrm, "LH_Nrm", & 
             "Latin hypercube estimate of Nrm [count/kg]", "count/kg", LH_zt )
        k = k + 1

      case ( 'LH_ricem' )
        iLH_ricem = k

        call stat_assign( iLH_ricem, "LH_ricem", & 
             "Latin hypercube estimate of ricem [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_Nim' )
        iLH_Nim = k

        call stat_assign( iLH_Nim, "LH_Nim", & 
             "Latin hypercube estimate of Nim [count/kg]", "count/kg", LH_zt )
        k = k + 1

      case ( 'LH_rsnowm' )
        iLH_rsnowm = k

        call stat_assign( iLH_rsnowm, "LH_rsnowm", & 
             "Latin hypercube estimate of rsnowm [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_Nsnowm' )
        iLH_Nsnowm = k

        call stat_assign( iLH_Nsnowm, "LH_Nsnowm", & 
             "Latin hypercube estimate of Nsnowm [count/kg]", "count/kg", LH_zt )
        k = k + 1


      case ( 'LH_rgraupelm' )
        iLH_rgraupelm = k

        call stat_assign( iLH_rgraupelm, "LH_rgraupelm", & 
             "Latin hypercube estimate of rgraupelm [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_Ngraupelm' )
        iLH_Ngraupelm = k

        call stat_assign( iLH_Ngraupelm, "LH_Ngraupelm", & 
             "Latin hypercube estimate of Ngraupelm [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_thlm' )
        iLH_thlm = k

        call stat_assign( iLH_thlm, "LH_thlm", & 
             "Latin hypercube estimate of thlm [K]", "K", LH_zt )
        k = k + 1

      case ( 'LH_rcm' )
        iLH_rcm = k

        call stat_assign( iLH_rcm, "LH_rcm", & 
             "Latin hypercube estimate of rcm [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_Ncm' )
        iLH_Ncm = k

        call stat_assign( iLH_Ncm, "LH_Ncm", & 
             "Latin hypercube estimate of Ncm [count/kg]", "count/kg", LH_zt )
        k = k + 1


      case ( 'LH_rvm' )
        iLH_rvm = k

        call stat_assign( iLH_rvm, "LH_rvm", & 
             "Latin hypercube estimate of rvm [kg/kg]", "kg/kg", LH_zt )
        k = k + 1

      case ( 'LH_wm' )
        iLH_wm = k

        call stat_assign( iLH_wm, "LH_wm", & 
             "Latin hypercube estimate of vertical velocity [m/s]", "m/s", LH_zt )
        k = k + 1

      case ( 'LH_cloud_frac' )
        iLH_cloud_frac = k

        ! Note: count is the udunits compatible unit
        call stat_assign( iLH_cloud_frac, "LH_cloud_frac", & 
             "Latin hypercube estimate of cloud fraction [count]", "count", LH_zt )
        k = k + 1

      case ( 'LH_wp2_zt' )
        iLH_wp2_zt = k
        call stat_assign( iLH_wp2_zt, "LH_wp2_zt", & 
             "Variance of the latin hypercube estimate of w [m^2/s^2]", "m^2/s^2", LH_zt )
        k = k + 1

      case ( 'LH_Ncp2_zt' )
        iLH_Ncp2_zt = k
        call stat_assign( iLH_Ncp2_zt, "LH_Ncp2_zt", & 
          "Variance of the latin hypercube estimate of Nc [count^2/kg^2]", "count^2/kg^2", LH_zt )
        k = k + 1

      case ( 'LH_Nrp2_zt' )
        iLH_Nrp2_zt = k
        call stat_assign( iLH_Nrp2_zt, "LH_Nrp2_zt", & 
          "Variance of the latin hypercube estimate of Nr [count^2/kg^2]", "count^2/kg^2", LH_zt )
        k = k + 1

      case ( 'LH_rcp2_zt' )
        iLH_rcp2_zt = k
        call stat_assign( iLH_rcp2_zt, "LH_rcp2_zt", & 
             "Variance of the latin hypercube estimate of rc [kg^2/kg^2]", "kg^2/kg^2", LH_zt )
        k = k + 1

      case ( 'LH_rtp2_zt' )
        iLH_rtp2_zt = k
        call stat_assign( iLH_rtp2_zt, "LH_rtp2_zt", & 
             "Variance of the latin hypercube estimate of rt [kg^2/kg^2]", "kg^2/kg^2", LH_zt )
        k = k + 1

      case ( 'LH_thlp2_zt' )
        iLH_thlp2_zt = k
        call stat_assign( iLH_thlp2_zt, "LH_thlp2_zt", & 
             "Variance of the latin hypercube estimate of thl [K^2]", "K^2", LH_zt )
        k = k + 1

      case ( 'LH_rrainp2_zt' )
        iLH_rrainp2_zt = k
        call stat_assign( iLH_rrainp2_zt, "LH_rrainp2_zt", & 
             "Variance of the latin hypercube estimate of rrain [kg^2/kg^2]", "kg^2/kg^2", LH_zt )
        k = k + 1

      case ( 'LH_rrainm_auto' )
        iLH_rrainm_auto = k
        call stat_assign( iLH_rrainm_auto, "LH_rrainm_auto", & 
             "Latin hypercube estimate of autoconversion [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case ( 'LH_rrainm_accr' )
        iLH_rrainm_accr = k
        call stat_assign( iLH_rrainm_accr, "LH_rrainm_accr", & 
             "Latin hypercube estimate of accretion [kg/kg/s]", "kg/kg/s", LH_zt )
        k = k + 1

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_LH_zt:  ', trim( vars_LH_zt(i) )

        l_error = .true.  ! This will stop the run.

      end select

    end do

    return
  end subroutine stats_init_LH_zt

end module stats_LH_zt
