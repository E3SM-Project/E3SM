!-------------------------------------------------------------------------------
! $Id: parameters_model.F90 7226 2014-08-19 15:52:41Z betlej@uwm.edu $
!===============================================================================
module parameters_model

! Description:
!   Contains model parameters that are determined at run time rather than
!   compile time.
!
! References:
!   None
!-------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Default scope

  ! Maximum magnitude of PDF parameter 'mixt_frac'. 
  real( kind = core_rknd ), public :: mixt_frac_max_mag

!$omp threadprivate(mixt_frac_max_mag)

  ! Model parameters and constraints setup in the namelists
  real( kind = core_rknd ), public ::  & 
    T0 = 300._core_rknd,       & ! Reference temperature (usually 300)  [K]
    ts_nudge = 0._core_rknd      ! Timescale of u/v nudging             [s]

#ifdef GFDL
 real( kind = core_rknd ), public ::  &   ! h1g, 2010-06-15
    cloud_frac_min    ! minimum cloud fraction for droplet #
!$omp threadprivate(cloud_frac_min)
#endif


!$omp threadprivate(T0, ts_nudge)

  real( kind = core_rknd), public :: &
    rtm_min = epsilon( rtm_min ), &             ! Value below which rtm will be nudged [kg/kg]
    rtm_nudge_max_altitude = 10000._core_rknd ! Highest altitude at which to nudge rtm [m]
!$omp threadprivate(rtm_min, rtm_nudge_max_altitude)

  integer, public :: & 
    sclr_dim = 0,        & ! Number of passive scalars
    edsclr_dim = 0,      & ! Number of eddy-diff. passive scalars
    hydromet_dim = 0       ! Number of hydrometeor species

!$omp threadprivate(sclr_dim, edsclr_dim, hydromet_dim)

  real( kind = core_rknd ), dimension(:), allocatable, public :: & 
    sclr_tol ! Threshold(s) on the passive scalars  [units vary]

!$omp threadprivate(sclr_tol)

  real( kind = selected_real_kind(6) ), public :: PosInf

!$omp threadprivate(PosInf)

  public :: setup_parameters_model 

  contains

!-------------------------------------------------------------------------------
  subroutine setup_parameters_model &
             ( T0_in, ts_nudge_in, &
               hydromet_dim_in, & 
               sclr_dim_in, sclr_tol_in, edsclr_dim_in &
#ifdef GFDL
              , cloud_frac_min_in &    ! hlg, 2010-6-15
#endif
    
              )

! Description:
!   Sets parameters to their initial values
!
! References:
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only: Skw_max_mag, Skw_max_mag_sqd

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt, allocated, transfer

    ! Constants
    integer(kind=4), parameter :: nanbits = 2139095040

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      T0_in,        & ! Ref. temperature             [K]
      ts_nudge_in     ! Timescale for u/v nudging    [s]

#ifdef GFDL
    real( kind = core_rknd ), intent(in) ::  cloud_frac_min_in  ! h1g, 2010-06-15
#endif


    integer, intent(in) :: & 
      hydromet_dim_in,  & ! Number of hydrometeor species
      sclr_dim_in,      & ! Number of passive scalars
      edsclr_dim_in       ! Number of eddy-diff. passive scalars

    real( kind = core_rknd ), intent(in), dimension(sclr_dim_in) :: & 
      sclr_tol_in     ! Threshold on passive scalars

    ! --- Begin Code --- 
     
    ! Formula from subroutine pdf_closure, where sigma_sqd_w = 0.4 and Skw =
    ! Skw_max_mag in this formula.  Note that this is constant, but can't appear
    ! with a Fortran parameter attribute, so we define it here. 
    mixt_frac_max_mag = 1.0_core_rknd &
      - ( 0.5_core_rknd * ( 1.0_core_rknd - Skw_max_mag / &
      sqrt( 4.0_core_rknd * ( 1.0_core_rknd - 0.4_core_rknd )**3 &
      + Skw_max_mag_sqd ) ) ) ! Known magic number

    T0       = T0_in
    ts_nudge = ts_nudge_in

    hydromet_dim = hydromet_dim_in
    sclr_dim     = sclr_dim_in
    edsclr_dim   = edsclr_dim_in

    ! In a tuning run, this array has the potential to be allocated already
    if ( .not. allocated( sclr_tol ) ) then
      allocate( sclr_tol(1:sclr_dim) )
    else
      deallocate( sclr_tol )
      allocate( sclr_tol(1:sclr_dim) )
    end if

    sclr_tol(1:sclr_dim) = sclr_tol_in(1:sclr_dim)

    PosInf = transfer( nanbits, PosInf )

#ifdef GFDL
     cloud_frac_min = cloud_frac_min_in  ! h1g, 2010-06-15
#endif

    return
  end subroutine setup_parameters_model
!-------------------------------------------------------------------------------

end module parameters_model
