module qbo

! Stub version of qbo module

implicit none
private
save

!---------------------------------------------------------------------
! Public methods
!---------------------------------------------------------------------
public               :: qbo_readnl             ! read namelist
public               :: qbo_init               ! initialize qbo package
public               :: qbo_timestep_init      ! interpolate to current time
public               :: qbo_relax              ! relax zonal mean wind

logical, public, parameter :: qbo_use_forcing  = .FALSE.

contains

subroutine qbo_readnl(nlfile)

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Stub; do nothing.

end subroutine qbo_readnl

subroutine qbo_init

  ! Stub; do nothing.

end subroutine qbo_init

subroutine qbo_timestep_init

  ! Stub; do nothing.

end subroutine qbo_timestep_init

subroutine qbo_relax( state, pbuf, ptend )

  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
!--------------------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------------------
  type(physics_state), intent(in)    :: state                ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)              ! Physics buffer
  type(physics_ptend), intent(out)   :: ptend                ! individual parameterization tendencies

  ! Stub; do nothing except init unused ptend.
  call physics_ptend_init(ptend, state%psetcols, 'qbo (stub)')

end subroutine qbo_relax

end module qbo
