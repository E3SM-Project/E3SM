module BeTRSimulationFactory
  !
  ! !DESCRIPTION:
  !  factory to load the specific betr simulator
  !
  ! !USES:
  !
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use shr_log_mod , only : errMsg => shr_log_errMsg

  implicit none
  save
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  public :: create_betr_simulation
contains

!-------------------------------------------------------------------------------
  function create_betr_simulation(simulator_name) result(simulator)
    !DESCRIPTION
    ! create a betr simulation object
    !
    !USES
    use BeTRSimulation           , only : betr_simulation_type
    use BeTRSimulationStandalone , only : create_betr_simulation_standalone
    use BeTRSimulationCLM        , only : create_betr_simulation_clm
    use BeTRSimulationALM        , only : create_betr_simulation_alm
    implicit none
    !ARGUMENTS
    character(len=*), intent(in)         :: simulator_name
    class(betr_simulation_type), pointer :: simulator

    select case(trim(simulator_name))
       case ("standalone")
          allocate(simulator, source=create_betr_simulation_standalone())
       case ("clm")
          allocate(simulator, source=create_betr_simulation_clm())
       case ("alm")
          allocate(simulator, source=create_betr_simulation_alm())
       case default
          write(*, *) "ERROR: unknown simulator type '", &
               trim(simulator_name), "'."
          call endrun(msg=errMsg(mod_filename, __LINE__))
    end select
  end function create_betr_simulation


end module BeTRSimulationFactory
