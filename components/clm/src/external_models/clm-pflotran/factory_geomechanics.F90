
module Factory_Geomechanics_module

  use Simulation_Geomechanics_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: GeomechanicsInitialize

contains

! ************************************************************************** !

subroutine GeomechanicsInitialize(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  implicit none
  
  class(simulation_geomechanics_type) :: simulation

  ! NOTE: PETSc must already have been initialized here!
  call GeomechanicsInitializePostPETSc(simulation)
  
end subroutine GeomechanicsInitialize

! ************************************************************************** !

subroutine GeomechanicsInitializePostPETSc(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 02/10/15
  ! 

  use Simulation_Geomechanics_class
  use Simulation_Subsurface_class
  use Factory_Subsurface_module
  use Init_Common_module
  use Option_module
  use PM_Base_class
  use PM_Base_Pointer_module
  use PM_Geomechanics_Force_class
  use PMC_Base_class
  use PMC_Geomechanics_class
  use PFLOTRAN_Constants_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module
  use Geomechanics_Realization_class
  use Geomechanics_Regression_module
  use Simulation_Aux_module
  use Realization_Subsurface_class
  use Timestepper_Steady_class
  use Input_Aux_module
  use Logging_module
  use Output_Aux_module

  implicit none
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(simulation_geomechanics_type) :: simulation

  type(option_type), pointer :: option
  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization
  class(pmc_base_type), pointer :: cur_process_model_coupler
  type(gmdm_ptr_type), pointer :: dm_ptr
  class(pm_base_type), pointer :: cur_pm, prev_pm
  class(pm_geomech_force_type), pointer :: pm_geomech
  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(timestepper_steady_type), pointer :: timestepper
  character(len=MAXSTRINGLENGTH) :: string
  type(waypoint_type), pointer :: waypoint
  type(input_type), pointer :: input
  PetscErrorCode :: ierr

  option => simulation%option
  
  nullify(prev_pm)
  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_geomech_force_type)
        pm_geomech => cur_pm
        if (associated(prev_pm)) then
          prev_pm%next => cur_pm%next
        else
          simulation%process_model_list => cur_pm%next
        endif
        exit
      class default
    end select
    prev_pm => cur_pm
    cur_pm => cur_pm%next
  enddo
  
  call SubsurfaceInitializePostPetsc(simulation)  
  simulation%process_model_coupler_list%is_master = PETSC_TRUE
    
  if (option%geomech_on) then
    simulation%geomech_realization => GeomechRealizCreate(option)
    geomech_realization => simulation%geomech_realization
    subsurf_realization => simulation%realization
    subsurf_realization%output_option => OutputOptionDuplicate(simulation%output_option)
    geomech_realization%input => InputCreate(IN_UNIT,option%input_filename,option)
    call GeomechicsInitReadRequiredCards(geomech_realization)
    pmc_geomech => PMCGeomechanicsCreate()
    pmc_geomech%name = 'PMCGeomech'
    simulation%geomech_process_model_coupler => pmc_geomech
    pmc_geomech%option => option
    pmc_geomech%checkpoint_option => simulation%checkpoint_option
    pmc_geomech%waypoint_list => simulation%waypoint_list_subsurface
    pmc_geomech%pm_list => pm_geomech
    pmc_geomech%pm_ptr%pm => pm_geomech
    pmc_geomech%geomech_realization => simulation%geomech_realization
    pmc_geomech%subsurf_realization => simulation%realization
    timestepper => TimestepperSteadyCreate()
    pmc_geomech%timestepper => timestepper
    ! set up logging stage
    string = trim(pmc_geomech%name) // 'Geomechanics'
    call LoggingCreateStage(string,pmc_geomech%stage)

    input => InputCreate(IN_UNIT,option%input_filename,option)    
    string = 'GEOMECHANICS'
    call InputFindStringInFile(input,option,string)
    call InputFindStringErrorMsg(input,option,string)  
    geomech_realization%output_option => OutputOptionDuplicate(simulation%output_option)
    nullify(geomech_realization%output_option%output_snap_variable_list)
    nullify(geomech_realization%output_option%output_obs_variable_list)    
    geomech_realization%output_option%output_snap_variable_list => OutputVariableListCreate()
    geomech_realization%output_option%output_obs_variable_list => OutputVariableListCreate()
    call GeomechanicsInitReadInput(simulation,timestepper%solver,input)
    pm_geomech%output_option => geomech_realization%output_option



    ! Hijack subsurface waypoint to geomechanics waypoint
    ! Subsurface controls the output now
    ! Always have snapshot on at t=0
    pmc_geomech%waypoint_list%first%print_snap_output = PETSC_TRUE
    

    ! link geomech and flow timestepper waypoints to geomech way point list
    if (associated(simulation%geomech_process_model_coupler)) then
      if (associated(simulation%geomech_process_model_coupler% &
                     timestepper)) then
        simulation%geomech_process_model_coupler%timestepper%cur_waypoint => &
          pmc_geomech%waypoint_list%first
      endif

      if (associated(simulation%flow_process_model_coupler%timestepper)) then
        simulation%flow_process_model_coupler%timestepper%cur_waypoint => &
          pmc_geomech%waypoint_list%first
      endif
    endif
    
    ! print the waypoints when debug flag is on
    if (geomech_realization%geomech_debug%print_waypoints) then
      call WaypointListPrint(pmc_geomech%waypoint_list,option, &
                             geomech_realization%output_option)      
    endif

    ! initialize geomech realization
    call GeomechInitSetupRealization(simulation)

    call pm_geomech%PMGeomechForceSetRealization(geomech_realization)
    call pm_geomech%Setup()

    call pmc_geomech%SetupSolvers()

    ! Here I first calculate the linear part of the jacobian and store it
    ! since the jacobian is always linear with geomech (even when coupled with
    ! flow since we are performing sequential coupling). Although
    ! SNESSetJacobian is called, nothing is done there and PETSc just 
    ! re-uses the linear Jacobian at all iterations and times
    call GeomechForceJacobianLinearPart(timestepper%solver%J,geomech_realization)
    nullify(simulation%process_model_coupler_list)                             
  endif
  ! sim_aux: Create PETSc Vectors and VectorScatters
  if (option%ngeomechdof > 0) then

    call GeomechCreateGeomechSubsurfVec(subsurf_realization, &
                                        geomech_realization)
    call SimAuxCopySubsurfVec(simulation%sim_aux,subsurf_realization%field%work)

    call GeomechCreateSubsurfStressStrainVec(subsurf_realization, &
                                             geomech_realization)
    call SimAuxCopySubsurfGeomechVec(simulation%sim_aux, &
          geomech_realization%geomech_field%strain_subsurf)

    call GeomechRealizMapSubsurfGeomechGrid(subsurf_realization, &
                                            geomech_realization, &
                                            option)

    dm_ptr => GeomechDiscretizationGetDMPtrFromIndex( &
                geomech_realization%geomech_discretization, ONEDOF)

    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                              SUBSURF_TO_GEOMECHANICS)
    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                              GEOMECHANICS_TO_SUBSURF)
  endif

  call GeomechanicsRegressionCreateMapping(simulation%geomech_regression, &
                                           geomech_realization)

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if (associated(simulation%rt_process_model_coupler)) &
    simulation%rt_process_model_coupler%sim_aux => simulation%sim_aux
  if (option%ngeomechdof>0 .and. &
     associated(simulation%geomech_process_model_coupler)) &
    simulation%geomech_process_model_coupler%sim_aux => simulation%sim_aux
 
  ! set geomech as not master
  simulation%geomech_process_model_coupler%is_master = PETSC_FALSE
  ! link geomech and master
  simulation%process_model_coupler_list => &
    simulation%geomech_process_model_coupler
  ! link subsurface flow as peer
  simulation%process_model_coupler_list%peer => &
    simulation%flow_process_model_coupler  
    

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  if (associated(cur_process_model_coupler%peer)) then
    cur_process_model_coupler => cur_process_model_coupler%peer
    call cur_process_model_coupler%GetAuxData()
    call cur_process_model_coupler%SetAuxData()
    select type(pmc => cur_process_model_coupler)
      class is(pmc_geomechanics_type)
        call GeomechStoreInitialPressTemp(pmc%geomech_realization)
    end select
  endif    

  call GeomechanicsJumpStart(simulation)
  call InputDestroy(geomech_realization%input)
  
end subroutine GeomechanicsInitializePostPETSc

! ************************************************************************** !

subroutine GeomechanicsJumpStart(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Geomechanics_Realization_class
  use Option_module
  use Timestepper_Steady_class
  use Output_Aux_module
  use Output_module, only : Output, OutputPrintCouplers
  use Output_Geomechanics_module
  use Logging_module
  use Condition_Control_module

  implicit none

  type(simulation_geomechanics_type) :: simulation

  class(realization_geomech_type), pointer :: geomch_realization
  class(timestepper_steady_type), pointer :: master_timestepper
  class(timestepper_steady_type), pointer :: geomech_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: snapshot_plot_flag,observation_plot_flag,massbal_plot_flag
  PetscBool :: geomech_read
  PetscBool :: failure
  PetscErrorCode :: ierr

  geomch_realization => simulation%geomech_realization

  select type(ts => simulation%geomech_process_model_coupler%timestepper)
    class is(timestepper_steady_type)
      geomech_timestepper => ts
  end select
  nullify(master_timestepper)

  option => geomch_realization%option

  call PetscOptionsHasName(PETSC_NULL_OBJECT, &
                           PETSC_NULL_CHARACTER, "-vecload_block_size", &
                           failure, ierr);CHKERRQ(ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported for &
                       &surface-flow.'
    call printErrMsg(option)
    return
  endif
  
  geomech_timestepper%name = 'GEOMECHANICS'
 
  master_timestepper => geomech_timestepper

  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  geomech_read = PETSC_FALSE
  failure = PETSC_FALSE
  
  call OutputGeomechInit(master_timestepper%steps)

  ! pushed in INIT_STAGE()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

  ! popped in TS_STAGE()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)

end subroutine GeomechanicsJumpStart

! ************************************************************************** !

subroutine GeomechicsInitReadRequiredCards(geomech_realization)
  ! 
  ! Reads the required input file cards
  ! related to geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  ! 

  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Option_module
  
  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  
  option         => geomech_realization%option  
  input          => geomech_realization%input
  
! Read in select required cards
!.........................................................................

  ! GEOMECHANICS information
  string = "GEOMECHANICS"
  call InputFindStringInFile(input,option,string)
  if (InputError(input)) return
  option%ngeomechdof = 3  ! displacements in x, y, z directions
  option%n_stress_strain_dof = 6
  
  string = "GEOMECHANICS_GRID"
  call InputFindStringInFile(input,option,string)
  call GeomechanicsInit(geomech_realization,input,option)  


end subroutine GeomechicsInitReadRequiredCards

! ************************************************************************** !

subroutine GeomechanicsInit(geomech_realization,input,option)
  ! 
  ! Reads the required geomechanics data from input file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_module
  
  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  type(grid_unstructured_type), pointer :: ugrid
  character(len=MAXWORDLENGTH) :: card
  
  geomech_discretization       => geomech_realization%geomech_discretization
       
  input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)
    
    select case(trim(word))
      case ('TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','TYPE')
        call StringToUpper(word)

        select case(trim(word))
          case ('UNSTRUCTURED')
            geomech_discretization%itype = UNSTRUCTURED_GRID
            call InputReadNChars(input,option, &
                                 geomech_discretization%filename, &
                                 MAXSTRINGLENGTH, &
                                 PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','filename')

            geomech_discretization%grid  => GMGridCreate()
            ugrid => UGridCreate()
            call UGridRead(ugrid,geomech_discretization%filename,option)
            call UGridDecompose(ugrid,option)
            call CopySubsurfaceGridtoGeomechGrid(ugrid, &
                                                 geomech_discretization%grid, &
                                                 option)
            patch => GeomechanicsPatchCreate()
            patch%geomech_grid => geomech_discretization%grid
            geomech_realization%geomech_patch => patch
          case default
            option%io_buffer = 'Geomechanics supports only unstructured grid'
            call printErrMsg(option)
        end select
      case ('GRAVITY')
        call InputReadDouble(input,option,option%geomech_gravity(X_DIRECTION))
        call InputErrorMsg(input,option,'x-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomech_gravity(Y_DIRECTION))
        call InputErrorMsg(input,option,'y-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomech_gravity(Z_DIRECTION))
        call InputErrorMsg(input,option,'z-direction','GEOMECH GRAVITY')
        if (option%myrank == option%io_rank .and. &
            option%print_to_screen) &
            write(option%fid_out,'(/," *GEOMECH_GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1p3e12.4 &
            & )') option%geomech_gravity(1:3)
      case default    
        call InputKeywordUnrecognized(word,'GEOMECHANICS_GRID',option)
    end select
  enddo
    
end subroutine GeomechanicsInit

! ************************************************************************** !

subroutine GeomechanicsInitReadInput(simulation,geomech_solver, &
                                     input)
  ! 
  ! Reads the geomechanics input data
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  ! 
  use Simulation_Geomechanics_class
  use Option_module
  use Input_Aux_module
  use String_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Material_module
  use Geomechanics_Region_module
  use Geomechanics_Debug_module
  use Geomechanics_Strata_module
  use Geomechanics_Condition_module
  use Geomechanics_Coupler_module
  use Geomechanics_Regression_module
  use Output_Aux_module
  use Output_Tecplot_module
  use Solver_module
  use Units_module
  use Waypoint_module
  use Utility_module, only : DeallocateArray, UtilityReadArray

  ! Still need to add other geomech modules for output, etc once created
  
  implicit none
  
  class(simulation_geomechanics_type) :: simulation
  type(solver_type) :: geomech_solver
  type(input_type), pointer :: input
  
  class(realization_geomech_type), pointer :: geomech_realization
  type(option_type), pointer :: option
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_material_property_type),pointer :: geomech_material_property
  type(waypoint_type), pointer :: waypoint
  type(geomech_grid_type), pointer :: grid
  type(gm_region_type), pointer :: region
  type(geomech_debug_type), pointer :: debug
  type(geomech_strata_type), pointer :: strata
  type(geomech_condition_type), pointer :: condition
  type(geomech_coupler_type), pointer :: coupler
  type(output_option_type), pointer :: output_option
  type(waypoint_list_type), pointer :: waypoint_list
  PetscReal :: units_conversion

  character(len=MAXWORDLENGTH) :: word, internal_units
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: string
  character(len=1) :: backslash
  
  PetscReal :: temp_real, temp_real2
  PetscReal, pointer :: temp_real_array(:)
  PetscInt :: i
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''
  
  waypoint_list => simulation%waypoint_list_geomechanics
  geomech_realization => simulation%geomech_realization
  option => simulation%option
  geomech_discretization => geomech_realization%geomech_discretization
  output_option => simulation%output_option
    
  if (associated(geomech_realization%geomech_patch)) grid => &
    geomech_realization%geomech_patch%geomech_grid
    
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)
    option%io_buffer = 'word :: ' // trim(word)
    call printMsg(option)   

    select case(trim(word))
    
      !.........................................................................
      ! Read geomechanics grid information
      case ('GEOMECHANICS_GRID')
        call InputSkipToEND(input,option,trim(word))
        
      !.........................................................................
      ! Read geomechanics material information
      case ('GEOMECHANICS_MATERIAL_PROPERTY')
        geomech_material_property => GeomechanicsMaterialPropertyCreate()

        call InputReadWord(input,option,geomech_material_property%name, &
                           PETSC_TRUE)
                           
        call InputErrorMsg(input,option,'name','GEOMECHANICS_MATERIAL_PROPERTY')
        call GeomechanicsMaterialPropertyRead(geomech_material_property,input, &
                                              option)
        call GeomechanicsMaterialPropertyAddToList(geomech_material_property, &
                                geomech_realization%geomech_material_properties)
        nullify(geomech_material_property)

      !.........................................................................
      case ('GEOMECHANICS_REGION')
        region => GeomechRegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','GEOMECHANICS_REGION')
        call printMsg(option,region%name)
        call GeomechRegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call GeomechRegionAddToList(region,geomech_realization%geomech_region_list)
        nullify(region)
 
      !.........................................................................
      case ('GEOMECHANICS_CONDITION')
        condition => GeomechConditionCreate(option)
        call InputReadWord(input,option,condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'GEOMECHANICS_CONDITION','name')
        call printMsg(option,condition%name)
        call GeomechConditionRead(condition,input,option)
        call GeomechConditionAddToList(condition,geomech_realization%geomech_conditions)
        nullify(condition)
        
     !.........................................................................
      case ('GEOMECHANICS_BOUNDARY_CONDITION')
        coupler =>  GeomechCouplerCreate(GM_BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Geomech Boundary Condition name')
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)
        
      !.........................................................................
      case ('GEOMECHANICS_SRC_SINK')
        coupler => GeomechCouplerCreate(GM_SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)
                 
      !.........................................................................
      case('NEWTON_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('GEOMECHANICS')
            call SolverReadNewton(geomech_solver,input,option)
        end select     
        
     !....................
      case ('LINEAR_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('GEOMECHANICS')
            call SolverReadLinear(geomech_solver,input,option)
        end select

      !.....................
      case ('GEOMECHANICS_REGRESSION')
        call GeomechanicsRegressionRead(simulation%geomech_regression,input,option)

      !.........................................................................
      case ('GEOMECHANICS_TIME')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word','GEOMECHANICS_TIME')
          select case(trim(word))
            case('COUPLING_TIMESTEP_SIZE')
              internal_units = 'sec'
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option, &
                                 'Coupling Timestep Size','GEOMECHANICS_TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                        'Coupling Timestep Size Time Units','GEOMECHANICS_TIME')
              geomech_realization%dt_coupling = &
                            temp_real*UnitsConvertToInternal(word, &
                            internal_units,option)
            case default
              call InputKeywordUnrecognized(word,'GEOMECHANICS_TIME',option)
            end select
        enddo
                
      !.........................................................................
      case ('GEOMECHANICS_DEBUG')
        call GeomechDebugRead(geomech_realization%geomech_debug,input,option)

      !.........................................................................
      case ('GEOMECHANICS_SUBSURFACE_COUPLING')
        option%geomech_subsurf_coupling = -1
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (word)
          case ('ONE_WAY_COUPLED')      
            option%geomech_subsurf_coupling = GEOMECH_ONE_WAY_COUPLED 
          case ('TWO_WAY_COUPLED')      
            option%geomech_subsurf_coupling = GEOMECH_TWO_WAY_COUPLED 
          case default
            call InputKeywordUnrecognized(word, &
                   'GEOMECHANICS_SUBSURFACE_COUPLING',option)
        end select
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','MAPPING_FILE')
          call InputReadNChars(input,option, &
                               grid%mapping_filename, &
                               MAXSTRINGLENGTH, &
                               PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','mapping_file')
          call GeomechSubsurfMapFromFilename(grid,grid%mapping_filename, &
                                             option)
        enddo
      !.........................................................................
      case ('GEOMECHANICS_OUTPUT')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT')
          call StringToUpper(word)
          select case(trim(word))
            case('TIMES')
              option%io_buffer = 'Subsurface times are now used for ' // &
              'geomechanics as well. No need for TIMES keyword under ' // &
              'GEOMECHANICS_OUTPUT.'
              call printWrnMsg(option)
            case('FORMAT')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT,&
                                                         &FORMAT') 
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputDefaultMsg(input,option, &
                                       'GEOMECHANICS_OUTPUT,FORMAT,HDF5,&
                                        &# FILES')
                  if (len_trim(word) > 1) then 
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                      case default
                        option%io_buffer = 'HDF5 keyword (' // trim(word) // &
                          ') not recongnized.  Use "SINGLE_FILE" or ' // &
                          '"MULTIPLE_FILES".'
                        call printErrMsg(option)
                    end select
                  endif
                case ('MAD')
                  output_option%print_mad = PETSC_TRUE
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'TECPLOT','GEOMECHANICS_OUTPUT,FORMAT') 
                  call StringToUpper(word)
                  output_option%tecplot_format = TECPLOT_FEQUADRILATERAL_FORMAT ! By default it is unstructured
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  call InputKeywordUnrecognized(word, &
                                 'GEOMECHANICS_OUTPUT,FORMAT',option)
              end select
            case default
              call InputKeywordUnrecognized(word, &
                             'GEOMECHANICS_OUTPUT',option)
          end select
        enddo
                        
      !.........................................................................
      case ('GEOMECHANICS_STRATIGRAPHY','GEOMECHANICS_STRATA')
        strata => GeomechStrataCreate()
        call GeomechStrataRead(strata,input,option)
        call GeomechRealizAddStrata(geomech_realization,strata)
        nullify(strata)       
      !.........................................................................
      case ('END_GEOMECHANICS')
        exit       
        
      !.........................................................................
      case default
        call InputKeywordUnrecognized(word, &
                                 'GeomechanicsInitReadInput',option)
    end select
  enddo
  
end subroutine GeomechanicsInitReadInput

! ************************************************************************** !

subroutine GeomechInitMatPropToGeomechRegions(geomech_realization)
  ! 
  ! This routine assigns geomech material
  ! properties to associated regions
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Strata_module
  use Geomechanics_Region_module
  use Geomechanics_Material_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module
  use Option_module

  implicit none
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  class(realization_geomech_type) :: geomech_realization
  
  PetscReal, pointer :: vec_p(:)
  
  PetscInt :: ivertex, local_id, ghosted_id, natural_id, geomech_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(geomech_grid_type), pointer :: grid
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: field
  type(geomech_strata_type), pointer :: strata
  type(geomech_patch_type), pointer :: patch  

  type(geomech_material_property_type), pointer :: geomech_material_property
  type(geomech_material_property_type), pointer :: null_geomech_material_property
  type(gm_region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  PetscReal, pointer :: imech_loc_p(:)
  
  option => geomech_realization%option
  geomech_discretization => geomech_realization%geomech_discretization
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch

  ! loop over all patches and allocation material id arrays
  if (.not.associated(patch%imat)) then
    allocate(patch%imat(patch%geomech_grid%ngmax_node))
    ! initialize to "unset"
    patch%imat = UNINITIALIZED_INTEGER
  endif

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  grid => patch%geomech_grid
  strata => patch%geomech_strata_list%first
  do
    if (.not.associated(strata)) exit
    ! Read in cell by cell material ids if they exist
    if (.not.associated(strata%region) .and. strata%active) then
      option%io_buffer = 'Reading of material prop from file for' // &
        ' geomech is not implemented.'
      call printErrMsgByRank(option)
    ! Otherwise, set based on region
    else if (strata%active) then
      update_ghosted_material_ids = PETSC_TRUE
      region => strata%region
      geomech_material_property => strata%material_property
      if (associated(region)) then
        istart = 1
        iend = region%num_verts
      else
        istart = 1
        iend = grid%nlmax_node
      endif
      do ivertex = istart, iend
        if (associated(region)) then
          local_id = region%vertex_ids(ivertex)
        else
          local_id = ivertex
        endif
        ghosted_id = grid%nL2G(local_id)
        patch%imat(ghosted_id) = geomech_material_property%id
      enddo
    endif
    strata => strata%next
  enddo
    
  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call GeomechRealizLocalToLocalWithArray(geomech_realization, &
                                            MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_geomech_material_property => GeomechanicsMaterialPropertyCreate()
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax_node
    ghosted_id = grid%nL2G(local_id)
    geomech_material_id = patch%imat(ghosted_id)
    if (geomech_material_id == 0) then ! accomodate inactive cells
      geomech_material_property = null_geomech_material_property
    else if ( geomech_material_id > 0 .and. &
              geomech_material_id <= &
              size(geomech_realization%geomech_material_property_array)) then
      geomech_material_property => &
         geomech_realization% &
           geomech_material_property_array(geomech_material_id)%ptr
      if (.not.associated(geomech_material_property)) then
        write(dataset_name,*) geomech_material_id
        option%io_buffer = 'No material property for geomech material id ' // &
                            trim(adjustl(dataset_name)) &
                            //  ' defined in input file.'
        call printErrMsgByRank(option)
      endif
    else if (Uninitialized(geomech_material_id)) then 
      write(dataset_name,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'Uninitialized geomech material id in patch at cell ' // &
                         trim(adjustl(dataset_name))
      call printErrMsgByRank(option)
    else if (geomech_material_id > size(geomech_realization% &
      geomech_material_property_array)) then
      write(option%io_buffer,*) geomech_material_id
      option%io_buffer = 'Unmatched geomech material id in patch:' // &
        adjustl(trim(option%io_buffer))
      call printErrMsgByRank(option)
    else
      option%io_buffer = 'Something messed up with geomech material ids. ' // &
        ' Possibly material ids not assigned to all grid cells. ' // &
        ' Contact Glenn/Satish!'
      call printErrMsgByRank(option)
    endif
    imech_loc_p(ghosted_id) = geomech_material_property%id
  enddo ! local_id - loop
  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  
  call GeomechanicsMaterialPropertyDestroy(null_geomech_material_property)
  nullify(null_geomech_material_property)
  
  call GeomechDiscretizationLocalToLocal(geomech_discretization,field%imech_loc, &
                                         field%imech_loc,ONEDOF)
  
end subroutine GeomechInitMatPropToGeomechRegions

! ************************************************************************** !

subroutine GeomechInitSetupRealization(simulation)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Simulation_Geomechanics_class
  use Geomechanics_Realization_class
  use Geomechanics_Global_module
  use Geomechanics_Force_module
  use Realization_Subsurface_class
  
  use Option_module
  use Waypoint_module
  
  implicit none
  
  class(simulation_geomechanics_type) :: simulation
  
  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization
  type(option_type), pointer :: option
  
  subsurf_realization => simulation%realization
  geomech_realization => simulation%geomech_realization
  option => subsurf_realization%option
  
  call GeomechRealizCreateDiscretization(geomech_realization)

  if (option%geomech_subsurf_coupling /= 0) then
    call GeomechCreateGeomechSubsurfVec(subsurf_realization, &
                                        geomech_realization)
    call GeomechCreateSubsurfStressStrainVec(subsurf_realization, &
                                              geomech_realization)

    call GeomechRealizMapSubsurfGeomechGrid(subsurf_realization, &
                                            geomech_realization, &
                                            option)
  endif
  call GeomechRealizLocalizeRegions(geomech_realization)
  call GeomechRealizPassFieldPtrToPatch(geomech_realization)
  call GeomechRealizProcessMatProp(geomech_realization)
  call GeomechRealizProcessGeomechCouplers(geomech_realization)
  call GeomechRealizProcessGeomechConditions(geomech_realization)
  call GeomechInitMatPropToGeomechRegions(geomech_realization)
  call GeomechRealizInitAllCouplerAuxVars(geomech_realization)  
  call GeomechRealizPrintCouplers(geomech_realization)  
  call GeomechGridElemSharedByNodes(geomech_realization,option)
  call GeomechForceSetup(geomech_realization)
  call GeomechGlobalSetup(geomech_realization)
    
  ! SK: We are solving quasi-steady state solution for geomechanics.
  ! Initial condition is not needed, hence CondControlAssignFlowInitCondGeomech
  ! is not needed, at this point.
  call GeomechForceUpdateAuxVars(geomech_realization)
  
end subroutine GeomechInitSetupRealization

! ************************************************************************** !

subroutine GeomechInitSetupSolvers(geomech_realization,realization, &
                                   convergence_context,solver)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use Option_module
  
  use Solver_module
  use Convergence_module
  use Discretization_module
  use Geomechanics_Force_module
  use Geomechanics_Discretization_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscpc.h"
  
  class(realization_geomech_type), pointer :: geomech_realization
  class(realization_subsurface_type), pointer :: realization
  type(convergence_context_type), pointer :: convergence_context
  type(solver_type), pointer :: solver

  type(option_type), pointer :: option
  SNESLineSearch :: linesearch
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  option => realization%option
  
  call printMsg(option,"  Beginning setup of GEOMECH SNES ")
    
  if (solver%J_mat_type == MATAIJ) then
    option%io_buffer = 'AIJ matrix not supported for geomechanics.'
    call printErrMsg(option)
  endif

  call SolverCreateSNES(solver,option%mycomm)  
  call SNESSetOptionsPrefix(solver%snes, "geomech_", &
                            ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)
        
  if (solver%Jpre_mat_type == '') then
    solver%Jpre_mat_type = solver%J_mat_type
  endif
  call GeomechDiscretizationCreateJacobian(geomech_realization% &
                                            geomech_discretization,NGEODOF, &
                                            solver%Jpre_mat_type, &
                                            solver%Jpre,option)

  solver%J = solver%Jpre
  call MatSetOptionsPrefix(solver%Jpre,"geomech_", &
                            ierr);CHKERRQ(ierr)
    
#if 1
  call SNESSetFunction(solver%snes,geomech_realization%geomech_field%disp_r, &
                       GeomechForceResidual, &
                       geomech_realization,ierr);CHKERRQ(ierr)

  call SNESSetJacobian(solver%snes,solver%J, &
                       solver%Jpre,GeomechForceJacobian, &
                       geomech_realization,ierr);CHKERRQ(ierr)
#endif

  ! by default turn off line search
  call SNESGetLineSearch(solver%snes,linesearch, ierr);CHKERRQ(ierr)
  call SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC, &
                              ierr);CHKERRQ(ierr)


  ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
  if (option%verbosity >= 1) then
    string = '-geomech_snes_view'
    call PetscOptionsInsertString(PETSC_NULL_OBJECT, &
                                   string, ierr);CHKERRQ(ierr)
  endif

  call SolverSetSNESOptions(solver,option)

  option%io_buffer = 'Solver: ' // trim(solver%ksp_type)
  call printMsg(option)
  option%io_buffer = 'Preconditioner: ' // trim(solver%pc_type)
  call printMsg(option)

  ! shell for custom convergence test.  The default SNES convergence test
  ! is call within this function.
  !TODO(geh): free this convergence context somewhere!
  option%io_buffer = 'DEALLOCATE GEOMECH CONVERGENCE CONTEXT somewhere!!!'
  convergence_context => ConvergenceContextCreate(solver,option, &
                                                  realization%patch%grid)
  call SNESSetConvergenceTest(solver%snes,ConvergenceTest, &
                              convergence_context, &
                              PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)                                                  

  call printMsg(option,"  Finished setting up GEOMECH SNES ")
    
end subroutine GeomechInitSetupSolvers

end module Factory_Geomechanics_module

