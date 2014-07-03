!================================================================================
! radiation unit driver
!================================================================================
module unit_driver

  use shr_kind_mod, only: r8=>SHR_KIND_R8
  use abortutils,   only: endrun
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver, pverp, begchunk, endchunk
  use cam_logfile,  only: iulog

  implicit none
  private
  save

  public :: unit_driver_run
  public :: unit_driver_init

  real(r8), pointer :: landm(:,:)       ! land fraction ramp

contains

!================================================================================
!================================================================================
  subroutine unit_driver_run( phys_state, pbuf2d, cam_out, cam_in, dtime, recno)
    use physics_types,    only: physics_state
    use physics_types,    only: physics_ptend
    use ppgrid,           only: begchunk, endchunk
    use time_manager,     only: get_step_size, timemgr_set_date_time
    use drv_input_data,   only: dates, secs, ntimes
    use camsrfexch,       only: cam_out_t, cam_in_t
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_chunk
    use radiation,        only: radiation_tend
    use rad_data_input,   only: get_rad_data_input

    implicit none

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    real(r8),            intent(in)    :: dtime
    integer,             intent(in)    :: recno

    integer :: c                                 ! chunk index

    type(physics_ptend) :: ptend
    real(r8) :: fsns(pcols)      ! Surface solar absorbed flux
    real(r8) :: fsnt(pcols)      ! Net column abs solar flux at model top
    real(r8) :: flns(pcols)      ! Srf longwave cooling (up-down) flux
    real(r8) :: flnt(pcols)      ! Net outgoing lw flux at model top
    real(r8) :: fsds(pcols)      ! Surface solar down flux
    real(r8) :: net_flx(pcols)
    type(physics_buffer_desc), pointer :: pbuf(:)

    call timemgr_set_date_time(dates(recno),secs(recno))

    call get_rad_data_input( phys_state, pbuf2d, cam_in, landm, recno=recno )

!$OMP PARALLEL DO PRIVATE (c, pbuf ,fsns,fsnt,flns,flnt,fsds )
    do c=begchunk,endchunk

       pbuf => pbuf_get_chunk(pbuf2d, c)

       call radiation_tend( phys_state(c), ptend, pbuf, cam_out(c), cam_in(c), &
                            cam_in(c)%landfrac,landm(:,c),cam_in(c)%icefrac, cam_in(c)%snowhland, &
                            fsns, fsnt, flns, flnt, fsds, net_flx)

    end do

  end subroutine unit_driver_run

!================================================================================
!================================================================================
  subroutine unit_driver_init
    use rad_data_input,   only: init_rad_data_input

    implicit none

    allocate ( landm( pcols, begchunk:endchunk ))

    call init_rad_data_input()

  end subroutine unit_driver_init

end module unit_driver
