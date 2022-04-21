module test_helper_mod

  implicit none

  public :: init_test_params_f90
  public :: cleanup_test_f90

contains

  ! This routine marks the params as inited in the F90 interface to homme,
  ! as well as hard codes some parameters that prescribe a cubed-sphere mesh
  subroutine init_test_params_f90 () bind(c)
    use control_mod,       only: topology, cubed_sphere_map, partmethod
    use params_mod,        only: SFCURVE
    use homme_context_mod, only: is_params_inited
    use physical_constants, only : domain_size, DD_PI

    ! Hard coded choices for cubed sphere
    topology = 'cube'
    cubed_sphere_map = 0
    partmethod = SFCURVE
    domain_size = 4.0D0*DD_PI

    is_params_inited = .true.
  end subroutine init_test_params_f90

  subroutine cleanup_test_f90 () bind(c)
    use schedtype_mod,     only: schedule
    use parallel_mod,      only: rrequest, srequest, global_shared_buf, status
    use homme_context_mod, only: is_parallel_inited
    use prim_driver_base,  only: prim_init1_cleanup

    ! Cleanup the schedule structure
    deallocate(Schedule(1)%SendCycle)
    deallocate(Schedule(1)%RecvCycle)
    deallocate(Schedule(1)%MoveCycle)
    deallocate(Schedule(1)%pIndx)
    deallocate(Schedule(1)%gIndx)
    deallocate(Schedule(1)%Local2Global)
    deallocate(Schedule)

    deallocate(rrequest)
    deallocate(srequest)
    deallocate(status)
    deallocate(global_shared_buf)

    call prim_init1_cleanup()

    is_parallel_inited = .false.
  end subroutine cleanup_test_f90

end module test_helper_mod
