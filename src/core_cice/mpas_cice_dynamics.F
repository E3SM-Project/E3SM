module cice_dynamics

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics, &
            run_dynamics

  !--------------------------------------------------------------------------
  ! Dynamics options
  !--------------------------------------------------------------------------
  
  ! "hex" or ??
  character(len=200), parameter :: dynamicsGridType = "hex";

  ! "weak", "variational" or "pwl" 
  character(len=200), parameter :: divergenceFormulation = "weak"
  !character(len=200), parameter :: divergenceFormulation = "variational"
  !character(len=200), parameter :: divergenceFormulation = "pwl"

  ! EVP Damping - ".false." or ".true."
  logical, parameter :: evpDamping = .true.
  !logical, parameter :: evpDamping = .false.

  ! Grid curvature - "planar" or "spherical"
  !character(len=200), parameter :: gridCurvature = "planar"
  character(len=200), parameter :: gridCurvature = "spherical"

contains

  !--------------------------------------------------------------------------
  
  subroutine init_dynamics(block, dt)

    use cice_dynamics_hex, only: init_dynamics_hex

    type (block_type), intent(inout) :: block
    real (kind=RKIND), intent(in) :: dt

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         normal, &
         hexvar, &
         hexwach, &
         hexpwl

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "normal", normal)
    call MPAS_pool_get_subpool(block % structs, "hexvar", hexvar)
    call MPAS_pool_get_subpool(block % structs, "hexwach", hexwach)
    call MPAS_pool_get_subpool(block % structs, "hexpwl", hexpwl)

    if (trim(dynamicsGridType) == "hex") then

       call init_dynamics_hex(&
            mesh,                  &
            dt,                    &
            boundary,              &
            normal,                &
            hexvar,                &
            hexwach,               &
            hexpwl,                &
            divergenceFormulation, &
            gridCurvature)

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
         tracers,       &
         icestate,      &
         boundary,      &
         normal,        &
         atmos_forcing, &
         hexdyn,        &
         hexweak,       &
         hexvar,        &
         hexwach,       &
         hexpwl,        &
         hexfor

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "normal", normal)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "hexweak", hexweak)
    call MPAS_pool_get_subpool(block % structs, "hexvar", hexvar)
    call MPAS_pool_get_subpool(block % structs, "hexwach", hexwach)
    call MPAS_pool_get_subpool(block % structs, "hexpwl", hexpwl)
    call MPAS_pool_get_subpool(block % structs, "hexfor", hexfor)

    if (trim(dynamicsGridType) == "hex") then

       call run_dynamics_hex(&
            mesh,                      &
            tracers,                   &
            icestate,                  &
            boundary,                  &
            normal,                    &
            atmos_forcing,             &
            hexdyn,                    &
            hexweak,                   &
            hexvar,                    &
            hexwach,                   &
            hexpwl,                    &
            hexfor,                    &
            dt, divergenceFormulation, &
            evpDamping, &
            nstep)
  
    end if

  end subroutine run_dynamics

  !--------------------------------------------------------------------------

end module cice_dynamics

