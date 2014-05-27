module cice_dynamics

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics, &
            run_dynamics

contains

  !--------------------------------------------------------------------------
  
  subroutine init_dynamics(block, dt)

    use cice_dynamics_hex, only: init_dynamics_hex

    type (block_type), intent(inout) :: block
    real (kind=RKIND), intent(in) :: dt

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         hexweak, &
         hexvar, &
         hexwach, &
         hexpwl, &
         configs

    character(len=strKIND), pointer :: &
         config_dynamics_grid_type

    configs => block % configs

    call MPAS_pool_get_config(configs, "config_dynamics_grid_type", config_dynamics_grid_type)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "hexweak", hexweak)
    call MPAS_pool_get_subpool(block % structs, "hexvar", hexvar)
    call MPAS_pool_get_subpool(block % structs, "hexwach", hexwach)
    call MPAS_pool_get_subpool(block % structs, "hexpwl", hexpwl)

    if (trim(config_dynamics_grid_type) == "hex") then

       call init_dynamics_hex(&
            mesh,                  &
            dt,                    &
            boundary,              &
            hexweak,               &
            hexvar,                &
            hexwach,               &
            hexpwl,                &
            configs)

    else

       write(*,*) "init_dynamics: Unsupported grid type"
       stop

    end if

  end subroutine init_dynamics

  !--------------------------------------------------------------------------

  subroutine run_dynamics(block, dt, nstep)

    use cice_dynamics_hex, only: run_dynamics_hex

    type (block_type), intent(inout) :: block
    real (kind=RKIND), intent(in) :: dt
    integer, intent(in) :: nstep

    type (MPAS_pool_type), pointer :: &
         mesh,          &
         configs,       &
         tracers,       &
         icestate,      &
         boundary,      &
         atmos_forcing, &
         ocean_forcing, &
         hexdyn,        &
         hexweak,       &
         hexvar,        &
         hexwach,       &
         hexpwl

    character(len=strKIND), pointer :: &
         config_dynamics_grid_type

    configs => block % configs

    call MPAS_pool_get_config(configs, "config_dynamics_grid_type", config_dynamics_grid_type)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
    call MPAS_pool_get_subpool(block % structs, "ocean_forcing", ocean_forcing)
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "hexweak", hexweak)
    call MPAS_pool_get_subpool(block % structs, "hexvar", hexvar)
    call MPAS_pool_get_subpool(block % structs, "hexwach", hexwach)
    call MPAS_pool_get_subpool(block % structs, "hexpwl", hexpwl)

    if (trim(config_dynamics_grid_type) == "hex") then

       call run_dynamics_hex(&
            mesh, &
            tracers,       icestate,      &
            boundary,                     &
            atmos_forcing, ocean_forcing, &
            hexdyn,                       &
            hexweak,       hexvar,        &
            hexwach,       hexpwl,        &
            dt,            configs,       &
            nstep)
  
    end if

  end subroutine run_dynamics

  !--------------------------------------------------------------------------

end module cice_dynamics

