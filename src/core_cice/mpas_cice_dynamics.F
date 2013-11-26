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

    if (trim(dynamicsGridType) == "hex") then

       call init_dynamics_hex(&
            block % mesh,          &
            dt,                    &
            block % boundary,      &
            block % normal,        &
            block % hexvar,        &
            block % hexpwl,        &
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

    if (trim(dynamicsGridType) == "hex") then

       call run_dynamics_hex(&
            block % mesh,          &
            block % icestate,      &
            block % boundary,      &
            block % normal,        &
            block % atmos_forcing, &
            block % hexdyn,        &
            block % hexweak,       &
            block % hexvar,        &
            block % hexpwl,        &
            block % hexfor,        &
            dt, divergenceFormulation, &
            evpDamping, &
            nstep)
  
    end if

  end subroutine run_dynamics

  !--------------------------------------------------------------------------

end module cice_dynamics

