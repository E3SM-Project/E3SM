module Factory_Surface_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Surface_class

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  public :: SurfaceInitialize, &
            SurfaceInitializePostPETSc, &
            SurfaceReadInput, &
            SurfaceJumpStart

contains

! ************************************************************************** !

subroutine SurfaceInitialize(simulation_base,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  use Input_Aux_module
  use Timestepper_Base_class
  use Simulation_Base_class
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(simulation_surface_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => SurfaceSimulationCreate(option)
  call SurfaceInitializePostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine SurfaceInitialize

! ************************************************************************** !

subroutine SurfaceInitializePostPETSc(simulation, option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  use Init_Common_module
  
  implicit none
  
  class(simulation_surface_type) :: simulation
  type(option_type), pointer :: option
  
  call SurfaceJumpStart(simulation)

end subroutine SurfaceInitializePostPETSc

! ************************************************************************** !

subroutine SurfaceJumpStart(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Realization_Surface_class
  use Option_module
  use Timestepper_Surface_class
  use Output_Aux_module
  use Output_module
  use Logging_module  
  use Condition_Control_module
  use Checkpoint_Surface_module
  use Output_Surface_module, only : OutputSurface, OutputSurfaceInit

  implicit none

  type(simulation_surface_type) :: simulation
  
  class(realization_surface_type), pointer :: surf_realization
  class(timestepper_surface_type), pointer :: master_timestepper
  class(timestepper_surface_type), pointer :: surf_flow_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: snapshot_plot_flag,observation_plot_flag,massbal_plot_flag
  PetscBool :: surf_flow_read
  PetscBool :: failure
  PetscErrorCode :: ierr
  PetscReal :: surf_flow_prev_dt

  surf_realization => simulation%surf_realization
  output_option => simulation%output_option
  
  select type(ts => simulation%surf_flow_process_model_coupler%timestepper)
    class is(timestepper_surface_type)
      surf_flow_timestepper => ts
  end select
  nullify(master_timestepper)
  
  option => surf_realization%option

  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr);CHKERRQ(ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported for &
                       &surface-flow.'
    call printErrMsg(option)
    return
  endif
  
  master_timestepper => surf_flow_timestepper

  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  surf_flow_read = PETSC_FALSE
  failure = PETSC_FALSE
  
#if 0
  if (option%restart_flag) then
    call SurfaceRestart(surf_realization,surf_flow_prev_dt,surf_flow_read)

    if (option%time /= option%surf_flow_time) then
      option%io_buffer = 'option%time does not match option%surf_flow_time' // &
        ' while restarting simulation. Check the restart files.'
      call printErrMsg(option)
    endif

    if (surf_flow_read) then
      surf_flow_timestepper%prev_dt = surf_flow_prev_dt
      surf_flow_timestepper%target_time = option%surf_flow_time
      call TSSetTime(surf_flow_timestepper%solver%ts,option%surf_flow_time, &
                     ierr);CHKERRQ(ierr)
    endif

  endif
#endif

  ! pushed in Init()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

  ! popped in TimestepperFinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)

  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_timestepper%max_time_step < 0) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputSurfaceInit(master_timestepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_timestepper%max_time_step >= 0) then
    if (output_option%print_initial_snap) snapshot_plot_flag = PETSC_TRUE
    if (output_option%print_initial_obs) observation_plot_flag = PETSC_TRUE
    if (output_option%print_initial_massbal) massbal_plot_flag = PETSC_FALSE
    !call OutputSurface(surf_realization,snapshot_plot_flag, &
    !                   observation_plot_flag,massbal_plot_flag)
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_timestepper%max_time_step < 1) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'') 
    return
  endif

  ! increment plot number so that 000 is always the initial condition, 
  ! and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (.not.associated(surf_flow_timestepper%cur_waypoint)) then
    option%io_buffer = &
      'Null flow waypoint list; final time likely equal to start time.'
    call printMsg(option)
    return
  else
    surf_flow_timestepper%dt_max = surf_flow_timestepper%cur_waypoint%dt_max
  endif
          
  surf_flow_timestepper%start_time_step = surf_flow_timestepper%steps + 1
  
  if (surf_realization%debug%print_couplers) then
  !  call OutputPrintSurfaceCouplers(surf_realization,ZERO_INTEGER)
  endif

end subroutine SurfaceJumpStart

! ************************************************************************** !

subroutine SurfaceReadInput(surf_realization,surf_flow_solver,waypoint_list, &
                            input)
  ! 
  ! This routine reads surface flow data from the input file
  ! grids.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/09/12
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Surface_Material_module
  use Realization_Surface_class
  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Dataset_Base_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Grid_Unstructured_Aux_module
  use Discretization_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Checkpoint_module
  use Strata_module
  use Debug_module
  use Units_module
  use Waypoint_module
  use Patch_module
  use Solver_module
  use Output_Aux_module
  use Output_Surface_module
  use Utility_module, only : DeallocateArray, UtilityReadArray

  implicit none

  class(realization_surface_type), pointer :: surf_realization
  type(solver_type) :: surf_flow_solver
  type(waypoint_list_type) :: waypoint_list
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  type(discretization_type),pointer :: discretization
  type(grid_type), pointer :: grid
  type(grid_unstructured_type), pointer :: un_str_sfgrid
  type(surface_material_property_type),pointer :: surf_material_property
  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  class(dataset_base_type), pointer :: dataset

  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  PetscReal :: units_conversion
  PetscReal, pointer :: temp_real_array(:)
  PetscInt :: i

  character(len=MAXWORDLENGTH) :: word 
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: temp_string
  character(len=MAXWORDLENGTH) :: card
  character(len=1) :: backslash

  PetscBool :: velocities
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate

  type(waypoint_type), pointer :: waypoint
  PetscReal :: temp_real, temp_real2

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++

  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  option => surf_realization%option
  
  discretization => surf_realization%discretization
  output_option => surf_realization%output_option

  patch => surf_realization%patch

  if (associated(patch)) grid => patch%grid

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SURFACE_FLOW')
    call StringToUpper(word)
    write(*,*) 'word :: ',trim(word)

    select case(trim(word))
      !.........................................................................
      ! Read surface grid information
      case ('SURF_GRID')
        call InputSkipToEND(input,option,trim(word))
      !.........................................................................
      case ('SURF_FLOW_FORMULATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call StringToUpper(word)
        select case(trim(word))
          !case ('KINEMATIC')
          !  option%surface_flow_formulation = KINEMATIC_WAVE
          case ('DIFFUSIVE')
            option%surface_flow_formulation = DIFFUSION_WAVE
          case default
            call InputKeywordUnrecognized(word, &
                  'SURFACE_FLOW,SURF_FLOW_FORMULATION',option)
        end select
      !.........................................................................
      case ('SURF_MAX_MANNING_VELOCITY')
        call InputReadDouble(input,option,temp_real)
        option%max_manning_velocity = temp_real
      !.........................................................................
      case ('SURF_MAX_INFILTRATION_VELOCITY')
        call InputReadDouble(input,option,temp_real)
        option%max_infiltration_velocity = temp_real
      !.........................................................................
      ! Read surface material information
      case ('SURF_MATERIAL_PROPERTY')
        surf_material_property => SurfaceMaterialPropertyCreate()

        call InputReadWord(input,option,surf_material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
        call SurfaceMaterialPropertyRead(surf_material_property,input,option)
        call SurfaceMaterialPropertyAddToList(surf_material_property, &
                                      surf_realization%surf_material_properties)
        nullify(surf_material_property)

      !.........................................................................
      case ('SURF_REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SURF_REGION')
        call printMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,surf_realization%surf_regions)
        nullify(region)

      !.........................................................................
      case ('SURF_FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'SURF_FLOW_CONDITION','name')
        call printMsg(option,flow_condition%name)
        call FlowConditionRead(flow_condition,input,option)
        call FlowConditionAddToList(flow_condition, &
                                    surf_realization%surf_flow_conditions)
        nullify(flow_condition)

      !.........................................................................
      case ('SURF_BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name')
        call CouplerRead(coupler,input,option)
        call RealizSurfAddCoupler(surf_realization,coupler)
        nullify(coupler)

      !.........................................................................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call RealizSurfAddStrata(surf_realization,strata)
        nullify(strata)
        
      !.........................................................................
      case ('SURF_INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name') 
        call CouplerRead(coupler,input,option)
        call RealizSurfAddCoupler(surf_realization,coupler)
        nullify(coupler)        

      !.........................................................................
      case ('SURF_SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call CouplerRead(coupler,input,option)
        call RealizSurfAddCoupler(surf_realization,coupler)
        nullify(coupler)        

      !.........................................................................
      case ('SURF_SUBSURFACE_COUPLING')
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call StringToUpper(word)
        select case(trim(word))
          case('DECOUPLED')
            option%subsurf_surf_coupling = DECOUPLED
          case('FULLY_COUPLED')
            option%subsurf_surf_coupling = FULLY_COUPLED
          case('SEQ_COUPLED')
            option%subsurf_surf_coupling = SEQ_COUPLED
          case default
            option%io_buffer = 'Invalid value for SURF_SUBSURFACE_COUPLING'
            call printErrMsg(option)
        end select
        call InputSkipToEND(input,option,trim(word))

      !.........................................................................
      case ('SURF_DEBUG')
        call DebugRead(surf_realization%debug,input,option)

      !.........................................................................
      case ('SURF_OUTPUT')
        velocities = PETSC_FALSE
        mass_flowrate = PETSC_FALSE
        energy_flowrate = PETSC_FALSE
        aveg_mass_flowrate = PETSC_FALSE
        aveg_energy_flowrate = PETSC_FALSE
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','SURF_OUTPUT')
          call StringToUpper(word)
          select case(trim(word))
            case('NO_FINAL','NO_PRINT_FINAL')
              output_option%print_final_snap = PETSC_FALSE
              output_option%print_final_obs = PETSC_FALSE
              output_option%print_final_massbal = PETSC_FALSE
            case('NO_INITIAL','NO_PRINT_INITIAL')
              output_option%print_initial_snap = PETSC_FALSE
              output_option%print_initial_obs = PETSC_FALSE
              output_option%print_initial_massbal = PETSC_FALSE
            case('PERMEABILITY')
              option%io_buffer = 'PERMEABILITY output must now be entered &
                                 &under OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_permeability = PETSC_TRUE
            case('POROSITY')
              option%io_buffer = 'POROSITY output must now be entered under &
                                 &OUTPUT/VARIABLES card.'
              call printErrMsg(option)
!              output_option%print_porosity = PETSC_TRUE
            case('PRINT_COLUMN_IDS')
              output_option%print_column_ids = PETSC_TRUE
            case('TIMES')
              internal_units = 'sec'
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'units','SURF_OUTPUT')
              units_conversion = UnitsConvertToInternal(word, &
                                 internal_units,option)
              temp_string = 'SURF_OUTPUT,TIMES'
              nullify(temp_real_array)
              call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
                                    temp_string,input,option)
              do i = 1, size(temp_real_array)
                waypoint => WaypointCreate()
                waypoint%time = temp_real_array(i)*units_conversion
                waypoint%print_snap_output = PETSC_TRUE
                call WaypointInsertInList(waypoint,waypoint_list)
              enddo
              call DeallocateArray(temp_real_array)
            case('OUTPUT_FILE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'SURF_OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_file = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'SURF_OUTPUT,PERIODIC,OUTPUT_FILE')
                case default
                  call InputKeywordUnrecognized(word, &
                    'SURF_OUTPUT,PERIODIC,OUTPUT_FILE',option)
              end select
            case('SCREEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment','OUTPUT,SCREEN')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_screen = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%screen_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'SURF_OUTPUT,PERIODIC,SCREEN')
                case default
                  call InputKeywordUnrecognized(word, &
                    'SURF_OUTPUT,PERIODIC,SCREEN',option)
              end select
            case('PERIODIC')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'SURF_OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  internal_units = 'sec'
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'SURF_OUTPUT,PERIODIC,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'SURF_OUTPUT,PERIODIC,TIME')
                  units_conversion = UnitsConvertToInternal(word, &
                                     internal_units,option)
                  output_option%periodic_snap_output_time_incr = temp_real* &
                                                            units_conversion
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then
                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'start time units', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word, &
                                         internal_units,option)
                      temp_real = temp_real * units_conversion
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'and', &
                                          'SURF_OUTPUT,PERIODIC,TIME"')
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'end time units', &
                                         'SURF_OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word, &
                                         internal_units,option)
                      temp_real2 = temp_real2 * units_conversion
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_snap_output = PETSC_TRUE
                        call WaypointInsertInList(waypoint,waypoint_list)
                        temp_real = temp_real + &
                                    output_option%periodic_snap_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_snap_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'between', &
                                          'SURF_OUTPUT,PERIODIC,TIME')
                    endif
                  endif
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_snap_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'SURF_OUTPUT,PERIODIC,TIMESTEP')
                case default
                  call InputKeywordUnrecognized(word, &
                    'SURF_OUTPUT,PERIODIC,TIMESTEP',option)
              end select
            case('PERIODIC_OBSERVATION')
              output_option%print_observation = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                'SURF_OUTPUT, PERIODIC_OBSERVATION')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  internal_units = 'sec'
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'SURF_OUTPUT,PERIODIC_OBSERVATION,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'SURF_OUTPUT,PERIODIC_OBSERVATION,TIME')
                  units_conversion = UnitsConvertToInternal(word, &
                                     internal_units,option) 
                  output_option%periodic_obs_output_time_incr = temp_real* &
                                                               units_conversion
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_obs_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'SURF_OUTPUT,PERIODIC_OBSERVATION,&
                                     &TIMESTEP')
                case default
                  call InputKeywordUnrecognized(word, &
                    'SURF_OUTPUT,PERIODIC_OBSERVATION,TIMESTEP',option)
              end select
            case('FORMAT')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword','OUTPUT,FORMAT') 
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputDefaultMsg(input,option, &
                                       'OUTPUT,FORMAT,HDF5,# FILES')
                  if (len_trim(word) > 1) then 
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                        output_option%times_per_h5_file = 1
                        call InputReadWord(input,option,word,PETSC_TRUE)
                        if (len_trim(word)>0) then
                          select case(trim(word))
                            case('TIMES_PER_FILE')
                              call InputReadInt(input,option, &
                                              output_option%times_per_h5_file)
                              call InputErrorMsg(input,option,'timestep &
                                     &increment','OUTPUT,FORMAT,MULTIPLE_FILES&
                                     &,TIMES_PER_FILE')
                            case default
                              call InputKeywordUnrecognized(word, &
                                'SURF_OUTPUT,FORMAT,HDF5,MULTIPLE_FILES',option)
                          end select
                        endif
                      case default
                        call InputKeywordUnrecognized(word, &
                          'SURF_OUTPUT,FORMAT,HDF5',option)
                    end select
                  endif
                case ('MAD')
                  output_option%print_mad = PETSC_TRUE
                case default
                  call InputKeywordUnrecognized(word, &
                         'SURF_OUTPUT,FORMAT',option)
              end select

            case('VELOCITY_AT_CENTER')
              velocities = PETSC_TRUE
            case ('HDF5_WRITE_GROUP_SIZE')
              call InputReadInt(input,option,option%hdf5_write_group_size)
              call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE', &
                                 'Group size')
            case('HYDROGRAPH')
              output_option%print_hydrograph = PETSC_TRUE
            case('PROCESSOR_ID')
              option%io_buffer = 'PROCESSOR_ID output must now be entered &
                                 &under OUTPUT/VARIABLES card as PROCESS_ID.'
              call printErrMsg(option)
!              output_option%print_iproc = PETSC_TRUE
            case('FLOWRATES','FLOWRATE')
              mass_flowrate = PETSC_TRUE
              energy_flowrate = PETSC_TRUE
            case('MASS_FLOWRATE')
              mass_flowrate = PETSC_TRUE
            case('ENERGY_FLOWRATE')
              energy_flowrate = PETSC_TRUE
            case('AVERAGE_FLOWRATES','AVERAGE_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
              aveg_energy_flowrate = PETSC_TRUE
            case('AVERAGE_MASS_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
            case('AVERAGE_ENERGY_FLOWRATE')
              aveg_energy_flowrate = PETSC_TRUE
            case('VARIABLES')
              call OutputSurfaceVariableRead(input,option, &
                     output_option%output_variable_list)
            case('AVERAGE_VARIABLES')
              call OutputSurfaceVariableRead(input,option, &
                     output_option%aveg_output_variable_list)
            case default
              call InputKeywordUnrecognized(word,'SURF_OUTPUT',option)
          end select
        enddo

        if (velocities) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_vel_cent = PETSC_TRUE
          if (output_option%print_hdf5) &
            output_option%print_hdf5_vel_cent = PETSC_TRUE
          if (output_option%print_vtk) &
            output_option%print_vtk_vel_cent = PETSC_TRUE
        endif
        if (mass_flowrate.or.energy_flowrate.or.aveg_mass_flowrate.or. &
            aveg_energy_flowrate) then
          if (output_option%print_hdf5) then
            output_option%print_hdf5_mass_flowrate = mass_flowrate
            output_option%print_hdf5_energy_flowrate = energy_flowrate
            output_option%print_hdf5_aveg_mass_flowrate = aveg_mass_flowrate
            output_option%print_hdf5_aveg_energy_flowrate = aveg_energy_flowrate
            if (aveg_mass_flowrate.or.aveg_energy_flowrate) then
              if (Equal(output_option%periodic_snap_output_time_incr,0.d0)) then
                option%io_buffer = 'Keyword: AVEGRAGE_FLOWRATES/ ' // &
                  'AVEGRAGE_MASS_FLOWRATE/ENERGY_FLOWRATE defined without' // &
                  ' PERIODIC TIME being set.'
                call printErrMsg(option)
              endif
            endif
            option%flow%store_fluxes = PETSC_TRUE
          else
            option%io_buffer='Output FLOWRATES/MASS_FLOWRATE/ENERGY_FLOWRATE ' // &
              'only available in HDF5 format'
            call printErrMsg(option)
          endif
        endif

      !.........................................................................
      case('NEWTON_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            call SolverReadNewton(surf_flow_solver,input,option)
        end select

      !.........................................................................
      case ('SURF_TIME')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word','SURF_TIME')
          internal_units = 'sec'
          select case(trim(word))
            case ('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Initial Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'Initial Timestep Size Time Units','TIME')
              surf_realization%dt_init = &
                temp_real*UnitsConvertToInternal(word,internal_units,option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Maximum Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'Maximum Timestep Size Time Units','TIME')
              surf_realization%dt_max = & 
                temp_real*UnitsConvertToInternal(word,internal_units,option)
            case('COUPLING_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Coupling Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'Coupling Timestep Size Time Units','TIME')
              surf_realization%dt_coupling = &
                temp_real*UnitsConvertToInternal(word,internal_units,option)
            case default
              call InputKeywordUnrecognized(word,'TIME',option)
            end select
        enddo
      !.........................................................................
      case ('SURF_DATASET')
      nullify(dataset)
      call DatasetRead(input,dataset,option)
      call DatasetBaseAddToList(dataset,surf_realization%datasets)
      nullify(dataset)
      
      !.........................................................................
      case ('SURF_RESTART')
        option%io_buffer = 'The SURF_RESTART card within SURFACE_FLOW &
                           &block has been deprecated.'      
        !option%surf_restart_flag = PETSC_TRUE
        !call InputReadNChars(input,option,option%surf_restart_filename, &
        !                     MAXSTRINGLENGTH,PETSC_TRUE)
        !call InputErrorMsg(input,option,'SURF_RESTART','Surface restart &
        !                                               &file name') 
        !call InputReadDouble(input,option,option%restart_time)
        !if (input%ierr == 0) then
        !  call printErrMsg(option,'Setting time to value not supported in &
        !                          &surface-flow')
        !endif
        !option%first_step_after_restart = PETSC_TRUE

!......................

      case ('SURF_CHECKPOINT')
        option%io_buffer = 'The SURF_CHECKPOINT card within SURFACE_FLOW &
                           &block has been deprecated.'      
!        call CheckpointRead(input,option,surf_realization%checkpoint_option, &
!                            surf_realization%waypoint_list)
        
!......................

      case('END_SURFACE_FLOW')
        exit
        
      case default
        call InputKeywordUnrecognized(word,'SURFACE_FLOW',option)
    end select
  enddo

  if (option%restart_flag .neqv. option%surf_restart_flag) then
    option%io_buffer='option%restart_flag /= option%surf_restart_flag'
    call printErrMsg(option)
  endif

end subroutine SurfaceReadInput

end module Factory_Surface_module
