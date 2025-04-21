module zm_eamxx_bridge_methods
  !-----------------------------------------------------------------------------
  use zm_eamxx_bridge_params, only: r8, pcols, pver, pverp, top_lev
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods copied from EAM
  public :: cldfrc_fice

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

  tmax_fice = tmelt - 10._r8        ! max temperature for cloud ice formation
  tmin_fice = tmax_fice - 30._r8    ! min temperature for cloud ice formation
  tmax_fsnow = tmelt                ! max temperature for transition to convective snow
  tmin_fsnow = tmelt - 5._r8        ! min temperature for transition to convective snow
  fice(:,:top_lev-1) = 0._r8
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

end module zm_eamxx_bridge_methods
