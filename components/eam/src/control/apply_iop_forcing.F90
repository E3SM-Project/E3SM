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

subroutine advance_iop_subsidence(scm_dt, ps_in, &
                             u_in, v_in, t_in, q_in, &
                             u_update, v_update, t_update, q_update)

  implicit none

  ! Input arguments
  real(r8), intent(in) :: scm_dt, ps_in
  real(r8), intent(in) :: u_in(plev), v_in(plev), t_in(plev)
  real(r8), intent(in) :: q_in(plev,pcnst)

  ! Output arguments
  real(r8), intent(out) :: u_update(plev), v_update(plev), t_update(plev), q_update(plev,pcnst)

  ! Local variables
  integer :: k, m, nlon
  real(r8) :: pintm1(plevp), pmidm1(plev), pdelm1(plev)
  real(r8) :: p_dep, phi_dep
  real(r8) :: omega_k

  ! Compute pressure levels
  nlon = 1
  call plevs0(nlon, plon, plev, ps_in, pintm1, pmidm1, pdelm1)

  ! Loop over vertical levels
  do k = 1, plev
    omega_k = wfld(k)    ! omega at level k [Pa/s]
    p_dep = pmidm1(k) - omega_k * scm_dt

    ! === U velocity ===
    call linear_interp_pressure(pmidm1, u_in, plev, p_dep, u_update(k))
    call linear_interp_pressure(pmidm1, v_in, plev, p_dep, v_update(k))
    call linear_interp_pressure(pmidm1, t_in, plev, p_dep, t_update(k))

    do m = 1, pcnst
      call linear_interp_pressure(pmidm1, q_in(:,m), plev, p_dep, q_update(k,m))
    end do

    ! Apply explicit thermal expansion (same as before)
    t_update(k) = t_update(k) + scm_dt * omega_k * t_in(k) * rair / (cpair * pmidm1(k))
  end do

end subroutine advance_iop_subsidence

subroutine linear_interp_pressure(p, phi, n, p_target, phi_interp)
  implicit none
  integer, intent(in) :: n
  real(r8), intent(in) :: p(n), phi(n), p_target
  real(r8), intent(out) :: phi_interp

  integer :: i
  real(r8) :: w

  ! Clamp if outside bounds
  if (p_target <= p(1)) then
    phi_interp = phi(1)
  else if (p_target >= p(n)) then
    phi_interp = phi(n)
  else
    do i = 1, n-1
      if (p_target >= p(i) .and. p_target < p(i+1)) then
        w = (p_target - p(i)) / (p(i+1) - p(i))
        phi_interp = (1.0_r8 - w) * phi(i) + w * phi(i+1)
        return
      end if
    end do
  end if
end subroutine linear_interp_pressure


!
!-----------------------------------------------------------------------
!

end module apply_iop_forcing_mod


