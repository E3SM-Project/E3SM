module geometry_interface_mod

  implicit none

  public :: init_cube_geometry_f90
  public :: cleanup_geometry_f90

contains

  subroutine init_cube_geometry_f90 (ne_in) bind(c)
    use iso_c_binding,  only: c_int
    use dimensions_mod, only: ne, npart
    use control_mod,    only: topology, cubed_sphere_map, partmethod
    use params_mod,     only: SFCURVE
    use homme_grid_mod, only: init_geometry_f90
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: ne_in

    ! Hard coded choices for cubed sphere
    topology = 'cube'
    cubed_sphere_map = 0
    partmethod = SFCURVE

    ! Set desired resolution
    ne = ne_in

    ! Rely on homme_grid_mod to init geometry
    ! Note: this does a bit more than we need, but it saves lines of code
    call init_geometry_f90 ()

  end subroutine init_cube_geometry_f90

  subroutine cleanup_geometry_f90 () bind(c)
    use schedtype_mod,     only: schedule
    use homme_context_mod, only: is_geometry_inited, elem, dom_mt
    use prim_driver_base,  only: prim_init1_cleanup

    ! Clean up homme_context stuff
    deallocate(elem)
    deallocate(dom_mt)

    ! Cleanup the schedule structure
    deallocate(Schedule(1)%SendCycle)
    deallocate(Schedule(1)%RecvCycle)
    deallocate(Schedule(1)%MoveCycle)
    deallocate(Schedule(1)%pIndx)
    deallocate(Schedule(1)%gIndx)
    deallocate(Schedule(1)%Local2Global)
    deallocate(Schedule)

    is_geometry_inited = .false.
  end subroutine cleanup_geometry_f90

end module geometry_interface_mod
