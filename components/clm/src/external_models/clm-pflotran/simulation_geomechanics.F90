module Simulation_Geomechanics_class

  use Option_module
  use Simulation_Subsurface_class
  use Geomechanics_Regression_module
  use PMC_Base_class
  use PMC_Subsurface_class
  use PMC_Geomechanics_class
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use PFLOTRAN_Constants_module
  use Waypoint_module
  use Simulation_Aux_module
  use Output_Aux_module 
  
  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  type, public, extends(simulation_subsurface_type) :: &
    simulation_geomechanics_type
    ! pointer to geomechanics coupler
    class(pmc_geomechanics_type), pointer :: geomech_process_model_coupler
    class(realization_geomech_type), pointer :: geomech_realization
    type(waypoint_list_type), pointer :: waypoint_list_geomechanics
    type(geomechanics_regression_type), pointer :: geomech_regression
  contains
    procedure, public :: Init => GeomechanicsSimulationInit
    procedure, public :: InitializeRun => GeomechanicsSimulationInitializeRun
    procedure, public :: InputRecord => GeomechanicsSimInputRecord
    procedure, public :: ExecuteRun => GeomechanicsSimulationExecuteRun
    procedure, public :: FinalizeRun => GeomechanicsSimulationFinalizeRun
    procedure, public :: Strip => GeomechanicsSimulationStrip
  end type simulation_geomechanics_type
  
  public :: GeomechanicsSimulationCreate, &
            GeomechanicsSimulationDestroy
  
contains

! ************************************************************************** !

function GeomechanicsSimulationCreate(option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Option_module

  implicit none

  type(option_type), pointer :: option

  class(simulation_geomechanics_type), pointer :: GeomechanicsSimulationCreate

  print *,'GeomechanicsSimulationCreate'

  allocate(GeomechanicsSimulationCreate)
  call GeomechanicsSimulationCreate%Init(option)

end function GeomechanicsSimulationCreate

! ************************************************************************** !

subroutine GeomechanicsSimulationInit(this, option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! Modified: Satish Karra, 06/01/2016
  ! 
  use Waypoint_module
  use Option_module

  implicit none

  class(simulation_geomechanics_type) :: this
  type(option_type), pointer :: option

  call SubsurfaceSimulationInit(this, option)
  nullify(this%geomech_realization)
  nullify(this%geomech_regression)
  this%waypoint_list_geomechanics => WaypointListCreate()

end subroutine GeomechanicsSimulationInit

! ************************************************************************** !

subroutine GeomechanicsSimulationInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Output_module
  use PMC_Geomechanics_class

  implicit none

  class(simulation_geomechanics_type) :: this

  call printMsg(this%option,'Simulation%InitializeRun()')
  call this%process_model_coupler_list%InitializeRun()

  if (this%option%restart_flag) then
    call printErrMsg(this%option,'add code for restart of GeomechanicsSimulation')
  endif

end subroutine GeomechanicsSimulationInitializeRun

! ************************************************************************** !

subroutine GeomechanicsSimInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  use Output_module
  
  implicit none
  
  class(simulation_geomechanics_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
 
  write(id,'(a29)',advance='no') 'simulation type: '
  write(id,'(a)') 'geomechanics'

  ! print output file information
  call OutputInputRecord(this%output_option,this%waypoint_list_geomechanics)

end subroutine GeomechanicsSimInputRecord

! ************************************************************************** !

subroutine GeomechanicsSimulationExecuteRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Waypoint_module
  use Simulation_Base_class
  use Timestepper_Base_class, only : TS_CONTINUE

  implicit none
  
  class(simulation_geomechanics_type) :: this
  
  PetscReal :: time
  PetscReal :: final_time
  PetscReal :: dt
  PetscViewer :: viewer

  time = this%option%time

  final_time = SimulationGetFinalWaypointTime(this)

  call printMsg(this%option,'GeomechanicsSimulationExecuteRun()')

  if (.not.associated(this%geomech_realization)) then
    call this%RunToTime(final_time)

  else

    ! If simulation is decoupled subsurfac-geomech simulation, set
    ! dt_coupling to be dt_max
    if (this%geomech_realization%dt_coupling == 0.d0) then
      this%option%io_buffer = 'Set non-zero COUPLING_TIME_SIZE in GEOMECHANICS_TIME.'
      call printErrMsg(this%option)
    else
      do
        if (time + this%geomech_realization%dt_coupling > final_time) then
          dt = final_time-time
        else
          dt = this%geomech_realization%dt_coupling
        endif

        time = time + dt
        this%geomech_process_model_coupler%timestepper%dt = dt
        call this%RunToTime(time)

        if (this%stop_flag /= TS_CONTINUE) exit ! end simulation
 
        if (time >= final_time) exit
      enddo
    endif
  endif

end subroutine GeomechanicsSimulationExecuteRun


! ************************************************************************** !

subroutine GeomechanicsSimulationFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! Modified by Satish Karra, 06/22/16

  use Timestepper_Steady_class

  implicit none

  class(simulation_geomechanics_type) :: this
  class(timestepper_steady_type), pointer :: geomech_timestepper

  call printMsg(this%option,'GeomechanicsSimulationFinalizeRun')

  call SubsurfaceFinalizeRun(this)
  !call GeomechanicsFinalizeRun(this)
  nullify(geomech_timestepper)
  if (associated(this%geomech_process_model_coupler)) then
    select type(ts => this%geomech_process_model_coupler%timestepper)
      class is(timestepper_steady_type)
        geomech_timestepper => ts
    end select
  endif
  call GeomechanicsRegressionOutput(this%geomech_regression, &
                                    this%geomech_realization, &
                                    geomech_timestepper)  


end subroutine GeomechanicsSimulationFinalizeRun

! ************************************************************************** !

subroutine GeomechanicsSimulationStrip(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! Modified by Satish Karra, 06/01/16
  ! 

  implicit none
  
  class(simulation_geomechanics_type) :: this
  
  call printMsg(this%option,'GeomechanicsSimulationStrip()')
  
  call SubsurfaceSimulationStrip(this)
  call GeomechanicsRegressionDestroy(this%geomech_regression)
  call WaypointListDestroy(this%waypoint_list_subsurface)
  call WaypointListDestroy(this%waypoint_list_geomechanics)
   
end subroutine GeomechanicsSimulationStrip

! ************************************************************************** !

subroutine GeomechanicsSimulationDestroy(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  implicit none
  
  class(simulation_geomechanics_type), pointer :: simulation
  
  call printMsg(simulation%option,'GeomehanicsSimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine GeomechanicsSimulationDestroy

end module Simulation_Geomechanics_class
