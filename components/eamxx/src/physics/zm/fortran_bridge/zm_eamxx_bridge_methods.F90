module zm_eamxx_bridge_methods
  !-----------------------------------------------------------------------------
  use zm_eamxx_bridge_params, only: r8, pcols, pver, pverp, top_lev
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods copied from EAM
  public :: cldfrc_fice
  public :: zm_tend_init
  public :: zm_physics_update
  public :: vertinterp

!===================================================================================================
contains
!===================================================================================================

! Copied from cloud_fraction.F90 (and adjusted indentation)
subroutine cldfrc_fice(ncol, t, fice, fsnow)
  !
  ! Compute the fraction of the total cloud water which is in ice phase.
  ! The fraction depends on temperature only. 
  ! This is the form that was used for radiation, the code came from cldefr originally
  ! 
  ! Author: B. A. Boville Sept 10, 2002
  !  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
  !-----------------------------------------------------------------------
  use zm_eamxx_bridge_physconst, only: tmelt

  ! Arguments
  integer,  intent(in)  :: ncol                 ! number of active columns
  real(r8), intent(in)  :: t(pcols,pver)        ! temperature

  real(r8), intent(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
  real(r8), intent(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

  ! Local variables
  real(r8) :: tmax_fice                         ! max temperature for cloud ice formation
  real(r8) :: tmin_fice                         ! min temperature for cloud ice formation
  real(r8) :: tmax_fsnow                        ! max temperature for transition to convective snow
  real(r8) :: tmin_fsnow                        ! min temperature for transition to convective snow

  integer :: i,k                                ! loop indexes

  !-----------------------------------------------------------------------

  tmax_fice   = tmelt     - 10._r8   ! max temperature for cloud ice formation
  tmin_fice   = tmax_fice - 30._r8   ! min temperature for cloud ice formation
  tmax_fsnow  = tmelt                ! max temperature for transition to convective snow
  tmin_fsnow  = tmelt - 5._r8        ! min temperature for transition to convective snow
  fice(:,:top_lev-1)  = 0._r8
  fsnow(:,:top_lev-1) = 0._r8

  ! Define fractional amount of cloud that is ice
  do k=top_lev,pver
    do i=1,ncol
      ! If warmer than tmax then water phase
      if (t(i,k) > tmax_fice) then
        fice(i,k) = 0.0_r8
      ! If colder than tmin then ice phase
      else if (t(i,k) < tmin_fice) then
        fice(i,k) = 1.0_r8
      ! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
      else 
        fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
      end if
      ! snow fraction partitioning
      ! If warmer than tmax then water phase
      if (t(i,k) > tmax_fsnow) then
        fsnow(i,k) = 0.0_r8
      ! If colder than tmin then ice phase
      else if (t(i,k) < tmin_fsnow) then
        fsnow(i,k) = 1.0_r8
      ! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
      else 
        fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
      end if
    end do
  end do
end subroutine cldfrc_fice

!===================================================================================================

! This mimics the functionality of physics_ptend_init()
subroutine zm_tend_init( ncol, pcols, pver, tend_s, tend_q, tend_u, tend_v )
  !-----------------------------------------------------------------------------
  ! Arguments
  integer,                         intent(in   ) :: ncol    ! number of local columns
  integer,                         intent(in   ) :: pcols   ! max number of columns for variable dimensions
  integer,                         intent(in   ) :: pver    ! number of local columns
  real(r8), dimension(pcols,pver), intent(inout) :: tend_s  ! tendency of dry static energy
  real(r8), dimension(pcols,pver), intent(inout) :: tend_q  ! tendency of water vapor
  real(r8), dimension(pcols,pver), intent(inout) :: tend_u  ! tendency of zonal wind
  real(r8), dimension(pcols,pver), intent(inout) :: tend_v  ! tendency of meridional wind
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: i,k
  !-----------------------------------------------------------------------------
  ! initialize tendency variables
  do i = 1,ncol
    do k = 1,pver
      tend_s(i,k) = 0
      tend_q(i,k) = 0
      tend_u(i,k) = 0
      tend_v(i,k) = 0
    end do
  end do
end subroutine

!===================================================================================================

! This combines functionality of:
! - physics_update()      [see physics_update_mod.F90]
! - physics_update_main() [see physics_types.F90]
subroutine zm_physics_update( ncol, pcols, dt, state_phis, state_zm, state_zi, &
                              state_p_mid, state_p_int, state_p_del, &
                              state_t, state_qv, ptend_s, ptend_q)
  use zm_eamxx_bridge_physconst, only: cpair
  !-----------------------------------------------------------------------------
  ! Arguments
  integer,                         intent(in   ) :: ncol            ! number of local columns
  integer,                         intent(in   ) :: pcols           ! max number of columns for variable dimensions
  real(r8),                        intent(in   ) :: dt              ! time step
  real(r8), dimension(pcols),      intent(in   ) :: state_phis      ! input state surface geopotential height
  real(r8), dimension(pcols,pver), intent(inout) :: state_zm        ! input state altitude at mid-levels
  real(r8), dimension(pcols,pverp),intent(inout) :: state_zi        ! input state altitude at interfaces
  real(r8), dimension(pcols,pver), intent(in   ) :: state_p_mid     ! input state mid-point pressure
  real(r8), dimension(pcols,pverp),intent(in   ) :: state_p_int     ! input state interface pressure
  real(r8), dimension(pcols,pver), intent(in   ) :: state_p_del     ! input state pressure thickness
  ! real(r8), dimension(pcols,pver), intent(inout) :: state_dse       ! input state dry static energy
  real(r8), dimension(pcols,pver), intent(inout) :: state_t         ! input state temperature
  real(r8), dimension(pcols,pver), intent(inout) :: state_qv        ! input state water vapor
  real(r8), dimension(pcols,pver), intent(in   ) :: ptend_s         ! tendency of dry static energy
  real(r8), dimension(pcols,pver), intent(in   ) :: ptend_q         ! tendency of water vapor
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: i,k
  !-----------------------------------------------------------------------------
  do i = 1,ncol
    do k = 1, pver
      ! update water vapor
      state_qv(i,k) = state_qv(i,k) + ptend_q(i,k) * dt
      ! update temperature - assume that dS is really dEn, En=enthalpy=c_p*T, then  dT = dEn/c_p, so, state%t += ds/c_p.
      state_t(i,k) = state_t(i,k) + ptend_s(i,k)/cpair * dt
    end do
  end do
    
  call zm_geopotential_t( ncol, pcols, state_p_int, state_p_mid, state_p_del, state_t, state_qv, state_zi, state_zm )

  ! skip DSE update for EAMxx
  ! do i = 1,ncol
  !   do k = 1, pver
  !     ! update dry static energy using updated temperature
  !     state_dse(i,k) = state_t(i,k)*cpair + gravit*state_zm(i,k) + state_phis(i)
  !   end do
  ! end do
  
end subroutine zm_physics_update

!===================================================================================================

! copied and modified from geopotential.F90
subroutine zm_geopotential_t( ncol, pcols, pint, pmid, pdel, t, q, zi, zm )
  use zm_eamxx_bridge_physconst, only: zvir, rair, gravit
  !----------------------------------------------------------------------- 
  ! Purpose: Compute the geopotential height (above the surface) at the
  ! midpoints and interfaces using the input temperatures and pressures
  !-----------------------------------------------------------------------------
  ! Arguments
  integer,                         intent(in   ) :: ncol    ! Number of columns
  integer,                         intent(in   ) :: pcols
  real(r8), dimension(pcols,pverp),intent(in   ) :: pint    ! Interface pressures
  real(r8), dimension(pcols,pver), intent(in   ) :: pmid    ! Midpoint pressures
  real(r8), dimension(pcols,pver), intent(in   ) :: pdel    ! layer thickness
  real(r8), dimension(pcols,pver), intent(in   ) :: t       ! temperature
  real(r8), dimension(pcols,pver), intent(in   ) :: q       ! specific humidity 
  real(r8), dimension(pcols,pverp),intent(inout) :: zi      ! Height above surface at interfaces
  real(r8), dimension(pcols,pver), intent(inout) :: zm      ! Geopotential height at mid level
  !-----------------------------------------------------------------------------
  ! Local variables
  integer  :: i,k                     ! Lon, level indices
  real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
  real(r8) :: hkl(ncol)               ! off-diagonal element
  real(r8) :: tv                      ! virtual temperature
  real(r8) :: tvfac                   ! Tv/T
  !-----------------------------------------------------------------------------
  ! The lowest height interface is at the surface, and thus zero by definition
  do i = 1,ncol
    zi(i,pverp) = 0.0_r8
  end do
  ! Compute zi, zm from bottom up
  do k = pver, 1, -1
    ! First set hydrostatic elements consistent with dynamics   
    do i = 1,ncol
      hkl(i) = pdel(i,k) / pmid(i,k)
      hkk(i) = 0.5_r8 * hkl(i)
    end do
    ! Now compute tv, zm, zi
    do i = 1,ncol
      tvfac   = 1._r8 + zvir * q(i,k)
      tv      = t(i,k) * tvfac
      zm(i,k) = zi(i,k+1) + (rair/gravit) * tv * hkk(i)
      zi(i,k) = zi(i,k+1) + (rair/gravit) * tv * hkl(i)
    end do
  end do
  return
end subroutine zm_geopotential_t

!===================================================================================================

! copied from components/eam/src/control/interpolate_data.F90
subroutine vertinterp(ncol, ncold, nlev, pmid, pout, arrin, arrout)
    !-----------------------------------------------------------------------
    ! Purpose: Vertically interpolate input array to output pressure level
    ! Copy values at boundaries.
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    integer , intent(in)  :: ncol              ! column dimension
    integer , intent(in)  :: ncold             ! declared column dimension
    integer , intent(in)  :: nlev              ! vertical dimension
    real(r8), intent(in)  :: pmid(ncold,nlev)  ! input level pressure levels
    real(r8), intent(in)  :: pout              ! output pressure level
    real(r8), intent(in)  :: arrin(ncold,nlev) ! input  array
    real(r8), intent(out) :: arrout(ncold)     ! output array (interpolated)
    !---------------------------Local variables-----------------------------
    integer i,k               ! indices
    integer kupper(ncold)     ! Level indices for interpolation
    real(r8) dpu              ! upper level pressure difference
    real(r8) dpl              ! lower level pressure difference
    logical found(ncold)      ! true if input levels found
    logical error             ! error flag
    !-----------------------------------------------------------------
    !
    ! Initialize index array and logical flags
    !
    do i=1,ncol
       found(i)  = .false.
       kupper(i) = 1
    end do
    error = .false.
    !
    ! Store level indices for interpolation.
    ! If all indices for this level have been found,
    ! do the interpolation
    !
    do k=1,nlev-1
       do i=1,ncol
          if ((.not. found(i)) .and. pmid(i,k)<pout .and. pout<=pmid(i,k+1)) then
             found(i) = .true.
             kupper(i) = k
          end if
       end do
    end do
    !
    ! If we've fallen through the k=1,nlev-1 loop, we cannot interpolate and
    ! must extrapolate from the bottom or top data level for at least some
    ! of the longitude points.
    !
    do i=1,ncol
       if (pout <= pmid(i,1)) then
          arrout(i) = arrin(i,1)
       else if (pout >= pmid(i,nlev)) then
          arrout(i) = arrin(i,nlev)
       else if (found(i)) then
          dpu = pout - pmid(i,kupper(i))
          dpl = pmid(i,kupper(i)+1) - pout
          arrout(i) = (arrin(i,kupper(i)  )*dpl + arrin(i,kupper(i)+1)*dpu)/(dpl + dpu)
       else
          error = .true.
       end if
    end do
    !
    ! Error check
    !
    ! if (error) then
    !    call endrun ('VERTINTERP: ERROR FLAG')
    ! end if

    return
  end subroutine vertinterp

!===================================================================================================

end module zm_eamxx_bridge_methods
