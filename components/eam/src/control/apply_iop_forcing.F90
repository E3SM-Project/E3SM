module apply_iop_forcing_mod

use shr_kind_mod,   only: r8 => shr_kind_r8, i8 => shr_kind_i8
use pmgrid
use constituents,   only: pcnst, cnst_get_ind
use pspect
use physconst,      only: rair,cpair
use cam_logfile,    only: iulog
use iop_data_mod

implicit none

public advance_iop_forcing
public advance_iop_nudging
public advance_iop_subsidence

!=========================================================================
contains
!=========================================================================

subroutine advance_iop_forcing(scm_dt, ps_in, &             ! In
                    u_in, v_in, t_in, q_in, t_phys_frc,&    ! In
                    u_update, v_update, t_update, q_update) ! Out

!-----------------------------------------------------------------------
!
! Purpose:
! Apply large scale forcing for t, q, u, and v as provided by the
!   case IOP forcing file.
!
! Author:
! Original version: Adopted from CAM3.5/CAM5
! Updated version for E3SM: Peter Bogenschutz (bogenschutz1@llnl.gov)
!  and replaces the forecast.F90 routine in CAM3.5/CAM5/CAM6/E3SMv1/E3SMv2
!
!-----------------------------------------------------------------------

  ! Input arguments
  real(r8), intent(in) :: ps_in             ! surface pressure [Pa]
  real(r8), intent(in) :: u_in(plev)        ! zonal wind [m/s]
  real(r8), intent(in) :: v_in(plev)        ! meridional wind [m/s]
  real(r8), intent(in) :: t_in(plev)        ! temperature [K]
  real(r8), intent(in) :: q_in(plev,pcnst)  ! q tracer array [units vary]
  real(r8), intent(in) :: t_phys_frc(plev)  ! temperature forcing from physics [K/s]
  real(r8), intent(in) :: scm_dt            ! model time step [s]

  ! Output arguments
  real(r8), intent(out) :: t_update(plev)      ! updated temperature [K]
  real(r8), intent(out) :: q_update(plev,pcnst)! updated q tracer array [units vary]
  real(r8), intent(out) :: u_update(plev)      ! updated zonal wind [m/s]
  real(r8), intent(out) :: v_update(plev)      ! updated meridional wind [m/s]

  ! Local variables
  real(r8) pmidm1(plev)  ! pressure at model levels
  real(r8) pintm1(plevp) ! pressure at model interfaces
  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
  real(r8) t_lsf(plev)       ! storage for temperature large scale forcing
  real(r8) q_lsf(plev,pcnst) ! storage for moisture large scale forcing
  real(r8) fac, t_expan

  integer i,k,m           ! longitude, level, constituent indices
  integer nlon

  !! Get vertical level profiles

  nlon = 1 ! number of columns for plevs0 routine
  call plevs0(nlon, plon, plev, ps_in, pintm1, pmidm1, pdelm1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Advance T and Q due to large scale forcing

  if (use_3dfrc) then
    t_lsf(:plev) = divt3d(:plev)
    q_lsf(:plev,:pcnst) = divq3d(:plev,:pcnst)
  else
    t_lsf(:plev) = divt(:plev)
    q_lsf(:plev,:pcnst) = divq(:plev,:pcnst)
  endif

  do k=1,plev
    ! Initialize thermal expansion term to zero.  This term is only
    !  considered if using the preq-x dycore and if three dimensional
    !  forcing is not provided by IOP forcing file.
    t_expan = 0._r8
#ifndef MODEL_THETA_L
    ! this term is already taken into account through
    !  LS vertical advection in theta-l dycore
    if (.not. use_3dfrc) then
      t_expan = scm_dt*wfld(k)*t_in(k)*rair/(cpair*pmidm1(k))
    endif
#endif
    t_update(k) = t_in(k) + t_expan + scm_dt*(t_phys_frc(k) + t_lsf(k))
    do m=1,pcnst
      q_update(k,m) = q_in(k,m) + scm_dt*q_lsf(k,m)
    end do
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Set U and V fields

  if (have_v .and. have_u .and. .not. dp_crm) then
    do k=1,plev
      u_update(k) = uobs(k)
      v_update(k) = vobs(k)
    enddo
  else
    do k=1,plev
      u_update(k) = u_in(k)
      v_update(k) = v_in(k)
    enddo
  endif

end subroutine advance_iop_forcing

!=========================================================================

subroutine advance_iop_nudging(scm_dt, ps_in, t_in, q_in, &       ! In
                             t_update, q_update, relaxt, relaxq ) ! Out

!-----------------------------------------------------------------------
!
! Purpose:
! Option to nudge t and q to observations as specified by the IOP file
!-----------------------------------------------------------------------

  ! Input arguments
  real(r8), intent(in) :: scm_dt           ! model time step [s]
  real(r8), intent(in) :: ps_in            ! surface pressure [Pa]
  real(r8), intent(in) :: t_in(plev)       ! temperature [K]
  real(r8), intent(in) :: q_in(plev)       ! water vapor mixing ratio [kg/kg]

  ! Output arguments
  real(r8), intent(out) :: t_update(plev)  ! updated temperature [K]
  real(r8), intent(out) :: q_update(plev)  ! updated water vapor [kg/kg]
  real(r8), intent(out) :: relaxt(plev)    ! relaxation of temperature [K/s]
  real(r8), intent(out) :: relaxq(plev)    ! relaxation of vapor [kg/kg/s]

  ! Local variables
  integer :: k, nlon
  real(r8) rtau(plev)
  real(r8) pmidm1(plev)  ! pressure at model levels
  real(r8) pintm1(plevp) ! pressure at model interfaces
  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)

  nlon = 1 ! number of columns for plevs0 routine
  call plevs0(nlon, plon, plev, ps_in, pintm1, pmidm1, pdelm1)

  ! Set relaxation arrays to zero
  do k=1,plev
    relaxt(k) = 0.0_r8
    relaxq(k) = 0.0_r8
  end do

  do k=1,plev

    if (pmidm1(k) .le. iop_nudge_tq_low*100._r8 .and. &
      pmidm1(k) .ge. iop_nudge_tq_high*100._r8) then

      ! Set the relaxation time scale
      rtau(k)   = iop_nudge_tscale
      rtau(k)   = max(scm_dt,rtau(k))
      relaxt(k) = -(t_update(k) - tobs(k))/rtau(k)
      relaxq(k) = -(q_update(k) - qobs(k))/rtau(k)

      t_update(k) = t_update(k) + relaxt(k)*scm_dt
      q_update(k) = q_update(k) + relaxq(k)*scm_dt

    endif

  end do

end subroutine advance_iop_nudging

!=========================================================================

subroutine advance_iop_subsidence(scm_dt, ps_in, &                    ! In
                             u_in, v_in, t_in, q_in, &                ! In
                             u_update, v_update, t_update, q_update ) ! Out

!-----------------------------------------------------------------------
!
! Purpose:
! Option to compute effects of large scale subsidence on T, q, u, and v.
! Code originated from CAM3.5/CAM5 Eulerian subsidence computation for SCM
! in the old forecast.f90 routine.
!-----------------------------------------------------------------------

  ! Input arguments
  real(r8), intent(in) :: scm_dt           ! model time step [s]
  real(r8), intent(in) :: ps_in            ! surface pressure [Pa]
  real(r8), intent(in) :: u_in(plev)       ! zonal wind [m/s]
  real(r8), intent(in) :: v_in(plev)       ! meridional wind [m/s]
  real(r8), intent(in) :: t_in(plev)       ! temperature [K]
  real(r8), intent(in) :: q_in(plev,pcnst) ! tracer [vary]

  ! Output arguments
  real(r8), intent(out) :: u_update(plev)      ! zonal wind [m/s]
  real(r8), intent(out) :: v_update(plev)      ! meridional wind [m/s]
  real(r8), intent(out) :: t_update(plev)      ! temperature [m/s]
  real(r8), intent(out) :: q_update(plev,pcnst)! tracer [vary]

  ! Local variables
  integer :: k, m, nlon
  real(r8) :: fac, weight
  real(r8) :: wfldint(plevp)
  real(r8) pmidm1(plev)  ! pressure at model levels
  real(r8) pintm1(plevp) ! pressure at model interfaces
  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Perform weight averaged in pressure of the omega profile to get from
  !   the midpoint grid to the interface grid

  nlon = 1 ! number of columns for plevs0 routine
  call plevs0(nlon, plon, plev, ps_in, pintm1, pmidm1, pdelm1)

  wfldint(1) = 0.0_r8

  do k=2,plev
     weight = (pintm1(k) - pmidm1(k-1))/(pmidm1(k) - pmidm1(k-1))
     wfldint(k) = (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
  end do

  wfldint(plevp) = 0.0_r8

  do k=2,plev-1
    fac = scm_dt/(2.0_r8*pdelm1(k))
    u_update(k) = u_in(k) &
      - fac*(wfldint(k+1)*(u_in(k+1) - u_in(k)) &
      + wfldint(k)*(u_in(k) - u_in(k-1)))
    v_update(k) = v_in(k) &
      - fac*(wfldint(k+1)*(v_in(k+1) - v_in(k)) &
      + wfldint(k)*(v_in(k) - v_in(k-1)))
    t_update(k) = t_in(k) &
      - fac*(wfldint(k+1)*(t_in(k+1) - t_in(k)) &
      + wfldint(k)*(t_in(k) - t_in(k-1)))

    do m=1,pcnst
      q_update(k,m) = q_in(k,m) &
        - fac*(wfldint(k+1) * (q_in(k+1,m) - q_in(k,m)) &
        + wfldint(k)*(q_in(k,m) - q_in(k-1,m)))
    enddo
  enddo

  ! Top and bottom levels next
  k = 1
  fac = scm_dt/(2.0_r8*pdelm1(k))
  u_update(k) = u_in(k) - fac*(wfldint(k+1)*(u_in(k+1) - u_in(k)))
  v_update(k) = v_in(k) - fac*(wfldint(k+1)*(v_in(k+1) - v_in(k)))
  t_update(k) = t_in(k) - fac*(wfldint(k+1)*(t_in(k+1) - t_in(k)))
  do m=1,pcnst
    q_update(k,m) = q_in(k,m) - fac*(wfldint(k+1)*(q_in(k+1,m) - q_in(k,m)))
  enddo

  k = plev
  fac = scm_dt/(2.0_r8*pdelm1(plev))
  u_update(k) = u_in(k) - fac*(wfldint(k)*(u_in(k) - u_in(k-1)))
  v_update(k) = v_in(k) - fac*(wfldint(k)*(v_in(k) - v_in(k-1)))
  t_update(k) = t_in(k) - fac*(wfldint(k)*(t_in(k) - t_in(k-1)))
  do m=1,pcnst
    q_update(k,m) = q_in(k,m) - fac*(wfldint(k)*(q_in(k,m) - q_in(k-1,m)))
  enddo

  ! thermal expansion term due to LS vertical advection
  do k=1,plev
    t_update(k) = t_update(k) + scm_dt*wfld(k)*t_in(k)*rair/(cpair*pmidm1(k))
  enddo

end subroutine advance_iop_subsidence

!
!-----------------------------------------------------------------------
!

end module apply_iop_forcing_mod


