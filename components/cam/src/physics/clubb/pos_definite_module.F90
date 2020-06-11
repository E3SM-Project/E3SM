!-------------------------------------------------------------------------
!$Id: pos_definite_module.F90 7140 2014-07-31 19:14:05Z betlej@uwm.edu $
!===============================================================================
module pos_definite_module

  implicit none

  public :: pos_definite_adj

  private ! Default Scope

  contains
!-----------------------------------------------------------------------
  subroutine pos_definite_adj & 
            ( dt, field_grid, field_np1,  & 
              flux_np1, field_n, field_pd, flux_pd )
! Description:
!   Applies a  flux conservative positive definite scheme to a variable

!   There are two possible grids:
!   (1) flux on zm  field on zt
!   then
!   flux_zt(k) = ( flux_zm(k) + flux_zm(k-1) ) / 2

!         CLUBB grid                  Smolarkiewicz grid
!   m +-- flux  zm(k)  --+               flux        k + 1/2
!   t +-- field zt(k)  --+               field, fout k
!   m +-- flux  zm(k-1) --+              flux        k - 1/2
!   t +-- field zt(k-1) --+

!   (1) flux on zt field on zm
!   then
!   flux_zm(k) = ( flux_zt(k) + flux_zt(k+1) ) / 2

!         CLUBB grid                  Smolarkiewicz grid
!   m +-- field  (k+1)  --+
!   t +-- flux   (k+1)  --+               flux        k + 1/2
!   m +-- field  (k)    --+               field, fout k
!   t +-- flux   (k)    --+               flux        k - 1/2


! References:
!   ``A Positive Definite Advection Scheme Obtained by
!     Nonlinear Renormalization of the Advective Fluxes'' Smolarkiewicz (1989)
!     Monthly Weather Review, Vol. 117, pp. 2626--2632
!-----------------------------------------------------------------------

    use grid_class, only: & 
      gr,   & ! Variable(s)
      ddzt, & ! Function
      ddzm    ! Function

    use constants_clubb, only :  & 
      eps, & ! Variable(s)
      zero_threshold

    use clubb_precision, only:  & 
      core_rknd ! Variable(s)

    use error_code, only: &
      clubb_at_least_debug_level

    implicit none

    ! External
    intrinsic :: eoshift, kind, any, min, max

    ! Input variables
    real( kind = core_rknd ), intent(in) :: & 
      dt ! Timestep    [s]

    character(len=2), intent(in) :: & 
      field_grid ! The grid of the field, either zt or zm

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      field_n ! The field (e.g. rtm) at n, prior to n+1

    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  & 
      flux_pd,  & ! Budget of the change in the flux term due to the scheme
      field_pd    ! Budget of the change in the mean term due to the scheme

    ! Output Variables

    real( kind = core_rknd ), intent(inout), dimension(gr%nz) :: & 
      field_np1,   & ! Field at n+1 (e.g. rtm in [kg/kg])
      flux_np1    ! Flux applied to field

    ! Local Variables
    integer :: & 
      kabove,   & ! # of vertical levels the flux higher point resides
      kbelow   ! # of vertical levels the flux lower point resides

    integer ::  & 
      k, kmhalf, kp1, kphalf ! Loop indices

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      flux_plus, flux_minus, & ! [F_i+1/2]^+ [F_i+1/2]^- in Smolarkiewicz 
      fout,                  & ! (A4) F_i{}^OUT, or the sum flux_plus+flux_minus
      flux_lim,              & ! Correction applied to flux at n+1
      field_nonlim          ! Temporary variable for calculation

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      dz_over_dt ! Conversion factor  [m/s]


!-----------------------------------------------------------------------

    ! If all the values are positive or the values at the previous
    ! timestep were negative, then just return
    if ( .not. any( field_np1 < 0._core_rknd ) .or. any( field_n < 0._core_rknd ) ) then
      flux_pd  = 0._core_rknd
      field_pd = 0._core_rknd
      return
    end if

    if ( field_grid == "zm" ) then
      kabove = 0
      kbelow = 1
    else if ( field_grid == "zt" ) then
      kabove = 1
      kbelow = 0
    else
      ! This is only necessary to avoid a compiler warning in g95
      kabove = -1
      kbelow = -1
      ! Joshua Fasching June 2008

      stop "Error in pos_def_adj"
    end if

    if ( clubb_at_least_debug_level( 1 ) ) then
      print *, "Correcting flux"
    end if

    do k = 1, gr%nz, 1

      ! Def. of F+ and F- from eqn 2 Smolarkowicz
      flux_plus(k)  =  max( zero_threshold, flux_np1(k) ) ! defined on flux levels
      flux_minus(k) = -min( zero_threshold, flux_np1(k) ) ! defined on flux levels

      if ( field_grid == "zm" ) then
        dz_over_dt(k) = ( 1._core_rknd/gr%invrs_dzm(k) ) / dt

      else if ( field_grid == "zt" ) then
        dz_over_dt(k) = ( 1._core_rknd/gr%invrs_dzt(k) ) / dt

      end if

    end do

    do k = 1, gr%nz, 1
      ! If the scalar variable is on the kth t-level, then
      ! Smolarkowicz's k+1/2 flux level is the kth m-level in CLUBB.

      ! If the scalar variable is on the kth m-level, then
      ! Smolarkowicz's k+1/2 flux level is the k+1 t-level in CLUBB.

      kphalf = min( k+kabove, gr%nz ) ! k+1/2 flux level
      kmhalf = max( k-kbelow, 1 )       ! k-1/2 flux level

      ! Eqn A4 from Smolarkowicz
      ! We place a limiter of eps to prevent a divide by zero, and
      !   after this calculation fout is on the scalar level, and
      !   fout is the total outward flux for the scalar level k.

      fout(k) = max( flux_plus(kphalf) + flux_minus(kmhalf), eps )

    end do


    do k = 1, gr%nz, 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! FIXME:
      ! We haven't tested this for negative values at the gr%nz level
      ! -dschanen 13 June 2008
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      kphalf = min( k+kabove, gr%nz ) ! k+1/2 flux level
      kp1    = min( k+1, gr%nz )      ! k+1 scalar level

      ! Eqn 10 from Smolarkowicz (1989)

      flux_lim(kphalf) & 
      = max( min( flux_np1(kphalf), & 
                  ( flux_plus(kphalf)/fout(k) ) * field_n(k) & 
                    * dz_over_dt(k) & 
                ), & 
             -( ( flux_minus(kphalf)/fout(kp1) ) * field_n(kp1) & 
                  * dz_over_dt(k) ) & 
           )
    end do

    ! Boundary conditions
    flux_lim(1) = flux_np1(1)
    flux_lim(gr%nz) = flux_np1(gr%nz)

    flux_pd = ( flux_lim - flux_np1 ) / dt

    field_nonlim = field_np1

    ! Apply change to field at n+1
    if ( field_grid == "zt" ) then

      field_np1 = -dt * ddzm( flux_lim - flux_np1 ) + field_np1

    else if ( field_grid == "zm" ) then

      field_np1 = -dt * ddzt( flux_lim - flux_np1 ) + field_np1

    end if

    ! Determine the total time tendency in field due to this calculation
    ! (for diagnostic purposes)
    field_pd = ( field_np1 - field_nonlim ) / dt

    ! Replace the non-limited flux with the limited flux
    flux_np1 = flux_lim

    return
  end subroutine pos_definite_adj

end module pos_definite_module
