!================================================================================
! stub unit driver
!================================================================================
module unit_driver

  use shr_kind_mod, only: r8=>SHR_KIND_R8

  implicit none
  private
  save

  public :: unit_driver_run
  public :: unit_driver_init

contains

  subroutine unit_driver_run(phys_state, pbuf2d, cam_out, cam_in, dtime, recno)
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk
    use camsrfexch,       only: cam_out_t, cam_in_t     
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    real(r8),            intent(in)    :: dtime
    integer,             intent(in)    :: recno

  end subroutine unit_driver_run

  subroutine unit_driver_init
  end subroutine unit_driver_init

end module unit_driver
