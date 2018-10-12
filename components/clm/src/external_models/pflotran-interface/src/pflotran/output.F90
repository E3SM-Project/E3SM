module Output_module

#include "petsc/finclude/petscdm.h"
  use petscdm
  use Logging_module 
  use Output_Aux_module

  use Output_HDF5_module
  use Output_Observation_module
  
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

#if defined(SCORPIO_WRITE)
  include "scorpiof.h"
#endif

  PetscInt, parameter ::  HDF5_FILE = 1

  
  PetscBool :: observation_first
  PetscBool :: hdf5_first
  PetscBool :: mass_balance_first

  public :: OutputInit, &
            Output, &
            OutputPrintCouplers, &
            OutputPrintRegions, &
            OutputVariableRead, &
            OutputFileRead, &
            OutputInputRecord, &
            OutputEnsureVariablesExist, &
            OutputFindNaNOrInfInVec

contains

! ************************************************************************** !

subroutine OutputInit(option,num_steps)
  ! 
  ! Initializes variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/22/09
  ! 
  use Option_module
  use Output_Common_module
  use Output_EKG_module

  implicit none
  
  type(option_type) :: option
  PetscInt :: num_steps
  
  call OutputCommonInit()
  call OutputObservationInit(num_steps)
  call OutputHDF5Init(num_steps)
  call OutputEKGInit(option,num_steps)

end subroutine OutputInit

! ************************************************************************** !

subroutine OutputFileRead(input,realization,output_option, &
                          waypoint_list,block_name)
  ! 
  ! Reads the *_FILE block within the OUTPUT block.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 02/23/2016
  ! 

  use Option_module
  use Input_Aux_module
  use Output_Aux_module
  use String_module
  use Realization_Subsurface_class
  use Waypoint_module
  use Units_module
  use Utility_module
  use Grid_module
  use Patch_module
  use Region_module

  implicit none

  type(input_type), pointer :: input
  class(realization_subsurface_type), pointer :: realization
  type(output_option_type), pointer :: output_option
  type(waypoint_list_type), pointer :: waypoint_list
  character(len=*) :: block_name
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(waypoint_type), pointer :: waypoint
  type(mass_balance_region_type), pointer :: new_massbal_region
  type(mass_balance_region_type), pointer :: cur_mbr
  PetscReal, pointer :: temp_real_array(:)

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units, internal_units
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: temp_real,temp_real2
  PetscReal :: units_conversion
  PetscInt :: k
  PetscBool :: added
  PetscBool :: vel_cent, vel_face
  PetscBool :: fluxes
  PetscBool :: mass_flowrate, energy_flowrate
  PetscBool :: aveg_mass_flowrate, aveg_energy_flowrate

  option => realization%option
  patch => realization%patch
  if (associated(patch)) grid => patch%grid

  vel_cent = PETSC_FALSE
  vel_face = PETSC_FALSE
  fluxes = PETSC_FALSE
  mass_flowrate = PETSC_FALSE
  energy_flowrate = PETSC_FALSE
  aveg_mass_flowrate = PETSC_FALSE
  aveg_energy_flowrate = PETSC_FALSE
  k = 0
  nullify(temp_real_array)

  select case(trim(block_name))
    case('SNAPSHOT_FILE')
    case('OBSERVATION_FILE')
      output_option%print_observation = PETSC_TRUE
    case('MASS_BALANCE_FILE')
      option%compute_mass_balance_new = PETSC_TRUE
  end select

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    string = 'OUTPUT,' // trim(block_name)
    call InputErrorMsg(input,option,'keyword',string)
    call StringToUpper(word)
    
    select case(trim(word))
!......................................
      case('NO_FINAL','NO_PRINT_FINAL')
        select case(trim(block_name))
          case('OBSERVATION_FILE')
            output_option%print_final_obs = PETSC_FALSE
          case('SNAPSHOT_FILE')
            output_option%print_final_snap = PETSC_FALSE
          case('MASS_BALANCE_FILE')
            output_option%print_final_massbal = PETSC_FALSE
        end select

!..........................................
      case('NO_INITIAL','NO_PRINT_INITIAL')
        select case(trim(block_name))
          case('OBSERVATION_FILE')
            output_option%print_initial_obs = PETSC_FALSE
          case('SNAPSHOT_FILE')
            output_option%print_initial_snap = PETSC_FALSE
          case('MASS_BALANCE_FILE')
            output_option%print_initial_massbal = PETSC_FALSE
        end select
        
!...............................
      case('TOTAL_MASS_REGIONS')
        select case(trim(block_name))
          case('OBSERVATION_FILE')
            option%io_buffer = 'TOTAL_MASS_REGIONS cannot be specified for &
                               &OUTPUT,OBSERVATION_FILE block.'
            call printErrMsg(option)
          case('SNAPSHOT_FILE')
            option%io_buffer = 'TOTAL_MASS_REGIONS cannot be specified for &
                               &OUTPUT,SNAPSHOT_FILE block.'
            call printErrMsg(option)
          case('MASS_BALANCE_FILE')
            string = 'OUTPUT,' // trim(block_name) // ',TOTAL_MASS_REGIONS'
            output_option%mass_balance_region_flag = PETSC_TRUE
            do
              ! Read region name:
              call InputReadPflotranString(input,option)
              call InputReadStringErrorMsg(input,option,string)
              if (InputCheckExit(input,option)) exit
              ! Region name found; read the region name
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword',string) 
              ! Create a new mass balance region
              new_massbal_region => OutputMassBalRegionCreate()
              new_massbal_region%region_name = trim(word)
              ! Add the new mass balance region to the list
              added = PETSC_FALSE
              if (.not.associated(output_option%mass_balance_region_list)) then
                output_option%mass_balance_region_list => new_massbal_region
              else
                cur_mbr => output_option%mass_balance_region_list
                do
                  if (.not.associated(cur_mbr)) exit
                  if (.not.associated(cur_mbr%next)) then
                    cur_mbr%next => new_massbal_region
                    added = PETSC_TRUE
                  endif
                  if (added) exit
                  cur_mbr => cur_mbr%next
                enddo
              endif
              nullify(new_massbal_region)
            enddo ! Read loop
        end select

!..................
      case('TIMES')
        string = 'OUTPUT,' // trim(block_name) // ',TIMES'
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'units',string)
        internal_units = 'sec'
        units_conversion = &
             UnitsConvertToInternal(word,internal_units,option)
        call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
             string,input,option)
        do k = 1, size(temp_real_array)
          waypoint => WaypointCreate()
          waypoint%time = temp_real_array(k)*units_conversion
          select case(trim(block_name))
            case('SNAPSHOT_FILE')
              waypoint%print_snap_output = PETSC_TRUE
            case('OBSERVATION_FILE')
              waypoint%print_obs_output = PETSC_TRUE
            case('MASS_BALANCE_FILE')
              waypoint%print_msbl_output = PETSC_TRUE
          end select    
          call WaypointInsertInList(waypoint,waypoint_list)
        enddo
        call DeallocateArray(temp_real_array)

!.....................
      case('PERIODIC')
        string = 'OUTPUT,' // trim(block_name) // ',PERIODIC'
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'periodic time increment type',string)
        call StringToUpper(word)
        select case(trim(word))
        !.............
          case('TIME')
            string = 'OUTPUT,' // trim(block_name) // ',PERIODIC,TIME'
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'time increment',string)
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'time increment units',string)
            internal_units = 'sec'
            units_conversion = UnitsConvertToInternal(word, &
                 internal_units,option) 
            select case(trim(block_name))
              case('SNAPSHOT_FILE')
                output_option%periodic_snap_output_time_incr = temp_real* &
                     units_conversion
              case('OBSERVATION_FILE')
                output_option%periodic_obs_output_time_incr = temp_real* &
                     units_conversion
              case('MASS_BALANCE_FILE')
                output_option%periodic_msbl_output_time_incr = temp_real* &
                     units_conversion
            end select
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr == 0) then
              if (StringCompareIgnoreCase(word,'between')) then
                call InputReadDouble(input,option,temp_real)
                call InputErrorMsg(input,option,'start time',string)
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'start time units',string)
                units_conversion = UnitsConvertToInternal(word, &
                     internal_units,option) 
                temp_real = temp_real * units_conversion
                call InputReadWord(input,option,word,PETSC_TRUE)
                if (.not.StringCompareIgnoreCase(word,'and')) then
                  input%ierr = 1
                endif
                call InputErrorMsg(input,option,'and',string)
                call InputReadDouble(input,option,temp_real2)
                call InputErrorMsg(input,option,'end time',string)
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'end time units',string)
                units_conversion = UnitsConvertToInternal(word, &
                     internal_units,option) 
                temp_real2 = temp_real2 * units_conversion
                select case(trim(block_name))
                  case('SNAPSHOT_FILE')
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
                  case('OBSERVATION_FILE')
                    do
                      waypoint => WaypointCreate()
                      waypoint%time = temp_real
                      waypoint%print_obs_output = PETSC_TRUE
                      call WaypointInsertInList(waypoint,waypoint_list)
                      temp_real = temp_real + &
                           output_option%periodic_obs_output_time_incr
                      if (temp_real > temp_real2) exit
                    enddo
                    output_option%periodic_obs_output_time_incr = 0.d0
                  case('MASS_BALANCE_FILE')
                    do
                      waypoint => WaypointCreate()
                      waypoint%time = temp_real
                      waypoint%print_msbl_output = PETSC_TRUE
                      call WaypointInsertInList(waypoint,waypoint_list)
                      temp_real = temp_real + &
                           output_option%periodic_msbl_output_time_incr
                      if (temp_real > temp_real2) exit
                    enddo
                    output_option%periodic_msbl_output_time_incr = 0.d0
                end select
              else
                input%ierr = 1
                call InputErrorMsg(input,option,'between',string)
              endif
            endif
        !.................
          case('TIMESTEP')
            string = 'OUTPUT,' // trim(block_name) // ',TIMESTEP'
            select case(trim(block_name))
              case('SNAPSHOT_FILE')
                call InputReadInt(input,option, &
                     output_option%periodic_snap_output_ts_imod)
              case('OBSERVATION_FILE')
                call InputReadInt(input,option, &
                     output_option%periodic_obs_output_ts_imod)
              case('MASS_BALANCE_FILE')
                call InputReadInt(input,option, &
                     output_option%periodic_msbl_output_ts_imod)
            end select
            call InputErrorMsg(input,option,'timestep increment',string)
        !.............
          case default
            call InputKeywordUnrecognized(word,'OUTPUT,PERIODIC',option)
        end select

!...................
      case('SCREEN')
        string = 'OUTPUT,' // trim(block_name) // ',SCREEN'
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'time increment',string)
        call StringToUpper(word)
        select case(trim(word))
          case('OFF')
            option%print_to_screen = PETSC_FALSE
          case('PERIODIC')
            string = trim(string) // ',PERIODIC'
            call InputReadInt(input,option,output_option%screen_imod)
            call InputErrorMsg(input,option,'timestep increment',string)
          case default
            call InputKeywordUnrecognized(word,string,option)
        end select

!...................
      case('FORMAT')
        string = 'OUTPUT,' // trim(block_name) // ',FORMAT'
        select case(trim(block_name))
          case('OBSERVATION_FILE')
            option%io_buffer = 'FORMAT cannot be specified within &
                 &the OUTPUT,OBSERVATION_FILE block. Observation output is &
                 &written in TECPLOT format only.'
            call printErrMsg(option)
          case('MASS_BALANCE_FILE')
            option%io_buffer = 'FORMAT cannot be specified within &
                 &the OUTPUT,MASS_BALANCE_FILE block. Mass balance output is &
                 &written in TECPLOT format only.'
            call printErrMsg(option)
        end select
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',string) 
        call StringToUpper(word)
        select case(trim(word))
        !..............
          case ('HDF5')
            string = trim(string) // ',HDF5'
            output_option%print_hdf5 = PETSC_TRUE
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr /= 0) then
              call InputDefaultMsg(input,option,string)
              output_option%print_single_h5_file = PETSC_TRUE
            else
              call StringToUpper(word)
              select case(trim(word))
              !....................
                case('SINGLE_FILE')
                  output_option%print_single_h5_file = PETSC_TRUE
              !.......................
                case('MULTIPLE_FILES')
                  string = trim(string) // ',MULTIPLE_FILES'
                  output_option%print_single_h5_file = PETSC_FALSE
                  output_option%times_per_h5_file = 1
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (input%ierr == 0) then
                    select case(trim(word))
                      case('TIMES_PER_FILE')
                        string = trim(string) // ',TIMES_PER_FILE'
                        call InputReadInt(input,option, &
                             output_option%times_per_h5_file)
                        call InputErrorMsg(input,option,'timestep increment', &
                                           string)
                      case default
                        call InputKeywordUnrecognized(word,string,option)
                    end select
                  endif
              !.............
                case default
                  call InputKeywordUnrecognized(word,string,option)
              end select
            endif
        !.............
          case ('MAD')
            output_option%print_mad = PETSC_TRUE
        !.............
          case default
            call InputKeywordUnrecognized(word,string,option)
        end select

!...................................
      case ('HDF5_WRITE_GROUP_SIZE')
        string = 'OUTPUT,' // trim(block_name) // ',HDF5_WRITE_GROUP_SIZE'
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'group size',string)

!......................
      case('VARIABLES')
        select case(trim(block_name))
          case('SNAPSHOT_FILE')           
            call OutputVariableRead(input,option, &
                 output_option%output_snap_variable_list)
          case('OBSERVATION_FILE')           
            call OutputVariableRead(input,option, &
                 output_option%output_obs_variable_list)
          case('MASS_BALANCE_FILE')
            option%io_buffer = 'A variable list cannot be specified within &
                 &the MASS_BALANCE_FILE block. Mass balance variables are &
                 &determined internally.'
            call printErrMsg(option)
        end select
        
!.............................
      case('PRINT_COLUMN_IDS')
        output_option%print_column_ids = PETSC_TRUE
        
!.............................
      case('DETAILED')
        select case(trim(block_name))
          case('MASS_BALANCE_FILE') 
            option%mass_bal_detailed = PETSC_TRUE
        end select

!...............................
      case('VELOCITY_AT_CENTER')
        vel_cent = PETSC_TRUE
      case('VELOCITY_AT_FACE')
        vel_face = PETSC_TRUE

!...................
      case('FLUXES')
        fluxes = PETSC_TRUE
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

!.................
      case default
        string = 'OUTPUT,' // trim(block_name)
        call InputKeywordUnrecognized(word,string,option)
    end select
  enddo
  

  if (vel_cent) then
    if (output_option%print_tecplot) &
         output_option%print_tecplot_vel_cent = PETSC_TRUE
    if (output_option%print_hdf5) &
         output_option%print_hdf5_vel_cent = PETSC_TRUE
    if (output_option%print_vtk) &
         output_option%print_vtk_vel_cent = PETSC_TRUE
  endif

  if (vel_face) then
    if (output_option%print_tecplot) &
         output_option%print_tecplot_vel_face = PETSC_TRUE
    if (output_option%print_hdf5) &
         output_option%print_hdf5_vel_face = PETSC_TRUE
  endif

  if (fluxes) then
    output_option%print_fluxes = PETSC_TRUE
  endif

  if(output_option%aveg_output_variable_list%nvars>0) then
    if(Equal(output_option%periodic_snap_output_time_incr,0.d0)) then
      option%io_buffer = 'Keyword: AVERAGE_VARIABLES defined without &
                         &PERIODIC TIME being set.'
      call printErrMsg(option)
    endif
    if(.not.output_option%print_hdf5) then
      option%io_buffer = 'Keyword: AVERAGE_VARIABLES only defined for &
                         &FORMAT HDF5'
      call printErrMsg(option)
    endif
  endif

  if (mass_flowrate.or.energy_flowrate.or.aveg_mass_flowrate &
       .or.aveg_energy_flowrate) then
    if (output_option%print_hdf5) then
      output_option%print_hdf5_mass_flowrate = mass_flowrate
      output_option%print_hdf5_energy_flowrate = energy_flowrate
      output_option%print_hdf5_aveg_mass_flowrate = aveg_mass_flowrate
      output_option%print_hdf5_aveg_energy_flowrate = aveg_energy_flowrate
      if(aveg_mass_flowrate.or.aveg_energy_flowrate) then
        if(Equal(output_option%periodic_snap_output_time_incr,0.d0)) then
          option%io_buffer = 'Keyword: AVEGRAGE_FLOWRATES/&
                             &AVEGRAGE_MASS_FLOWRATE/ENERGY_FLOWRATE &
                             &defined without PERIODIC TIME being set.'
          call printErrMsg(option)
        endif
      endif
      option%flow%store_fluxes = PETSC_TRUE
    endif
    if (associated(grid%unstructured_grid%explicit_grid)) then
      option%flow%store_fluxes = PETSC_TRUE
      output_option%print_explicit_flowrate = mass_flowrate
    endif
  endif
  
end subroutine OutputFileRead

! ************************************************************************** !

subroutine OutputVariableRead(input,option,output_variable_list)
  ! 
  ! This routine reads a variable from the input file.
  ! 
  ! Author: Gautam Bisht, LBNL; Glenn Hammond PNNL/SNL
  ! Date: 12/21/12
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Variables_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(output_variable_list_type), pointer :: output_variable_list
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_type), pointer :: output_variable
  PetscInt :: temp_int

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','VARIABLES')
    call StringToUpper(word)
    
    select case(trim(word))
      case ('MAXIMUM_PRESSURE')
        name = 'Maximum Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     MAXIMUM_PRESSURE)
      case ('LIQUID_PRESSURE')
        name = 'Liquid Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     LIQUID_PRESSURE)
      case ('LIQUID_SATURATION')
        name = 'Liquid Saturation'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_SATURATION,units, &
                                     LIQUID_SATURATION)
      case ('LIQUID_HEAD')
        name = 'Liquid Head'
        units = 'm'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_HEAD)
        
      case ('LIQUID_DENSITY')
        name = 'Liquid Density'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'MOLAR')) then
            units = 'kmol/m^3'
            temp_int = LIQUID_DENSITY_MOL
          else
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,LIQUID_DENSITY')
          endif
        else
          units = 'kg/m^3'
          temp_int = LIQUID_DENSITY
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     temp_int)
      case ('LIQUID_MOBILITY')
        name = 'Liquid Mobility'
        units = '1/Pa-s'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MOBILITY)
      case ('LIQUID_VISCOSITY')
        name = 'Liquid Viscosity'
        units = 'Pa-s'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_VISCOSITY)
      case ('LIQUID_ENERGY')
        name = 'Liquid Energy'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'PER_VOLUME')) then
            units = 'MJ/m^3'
            temp_int = ONE_INTEGER
          else
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,LIQUID_ENERGY')
          endif
        else
          units = 'MJ/kmol'
          temp_int = ZERO_INTEGER
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_ENERGY,temp_int)
    
      case ('GAS_PRESSURE')
        name = 'Gas Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     GAS_PRESSURE)
      case ('GAS_SATURATION')
        name = 'Gas Saturation'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_SATURATION,units, &
                                     GAS_SATURATION)
      case ('GAS_DENSITY')
        name = 'Gas Density'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'MOLAR')) then
            units = 'kmol/m^3'
            temp_int = GAS_DENSITY_MOL
          else
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,GAS_DENSITY')
          endif
        else
          units = 'kg/m^3'
          temp_int = GAS_DENSITY
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     temp_int)
      case ('GAS_MOBILITY')
        name = 'Gas Mobility'
        units = '1/Pa-s'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MOBILITY)
      case ('GAS_ENERGY')
        name = 'Gas Energy'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'PER_VOLUME')) then
            units = 'MJ/m^3'
            temp_int = ONE_INTEGER
          else
            input%ierr = 1
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,GAS_ENERGY')
          endif
        else
          units = 'MJ/kmol'
          temp_int = ZERO_INTEGER
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_ENERGY,temp_int)
      case ('OIL_PRESSURE')
        name = 'Oil Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     OIL_PRESSURE)
      case ('OIL_SATURATION')
        name = 'Oil Saturation'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_SATURATION,units, &
                                     OIL_SATURATION)
      case ('OIL_DENSITY')
        name = 'Oil Density'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'MOLAR')) then
            units = 'kmol/m^3'
            temp_int = OIL_DENSITY_MOL
          else
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,OIL_DENSITY')
          endif
        else
          units = 'kg/m^3'
          temp_int = OIL_DENSITY
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     temp_int)
      case ('OIL_MOBILITY')
        name = 'Oil Mobility'
        units = '1/Pa-s'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     OIL_MOBILITY)
      case ('OIL_VISCOSITY')
        name = 'Oil Viscosity'
        units = 'Pa-s'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     OIL_VISCOSITY)
      case ('OIL_ENERGY')
        name = 'Oil Energy'
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (StringCompareIgnoreCase(word,'PER_VOLUME')) then
            units = 'MJ/m^3'
            temp_int = ONE_INTEGER
          else
            input%ierr = 1
            call InputErrorMsg(input,option,'optional keyword', &
                               'VARIABLES,OIL_ENERGY')
          endif
        else
          units = 'MJ/kmol'
          temp_int = ZERO_INTEGER
        endif
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     OIL_ENERGY,temp_int)

      case ('BUBBLE_POINT')
        name = 'Bubble Point'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     BUBBLE_POINT)

      case ('ICE_SATURATION')
        name = 'Ice Saturation'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_SATURATION,units, &
                                     ICE_SATURATION)

      case ('LIQUID_MOLE_FRACTIONS')
        name = 'X_g^l'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MOLE_FRACTION, &
                                     option%air_id)
        name = 'X_l^l'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MOLE_FRACTION, &
                                     option%water_id)
      case ('GAS_MOLE_FRACTIONS')
        name = 'X_g^g'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MOLE_FRACTION, &
                                     option%air_id)
        name = 'X_l^g'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MOLE_FRACTION, &
                                     option%water_id)
      case ('LIQUID_MASS_FRACTIONS')
        name = 'w_g^l'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MASS_FRACTION, &
                                     option%air_id)
        name = 'w_l^l'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     LIQUID_MASS_FRACTION, &
                                     option%water_id)
      case ('GAS_MASS_FRACTIONS')
        name = 'w_g^g'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MASS_FRACTION, &
                                     option%air_id)
        name = 'w_l^g'
        units = ''
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_MASS_FRACTION, &
                                     option%water_id)
      case ('AIR_PRESSURE')
        name = 'Air Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     AIR_PRESSURE)
      case ('CAPILLARY_PRESSURE')
        name = 'Capillary Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     CAPILLARY_PRESSURE)
      case ('VAPOR_PRESSURE')
        name = 'Vapor Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     VAPOR_PRESSURE)
      case ('SATURATION_PRESSURE')
        name = 'Saturation Pressure'
        units = 'Pa'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_PRESSURE,units, &
                                     SATURATION_PRESSURE)
      case('THERMODYNAMIC_STATE')
        name = 'Thermodynamic State'
         units = ''
         output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                                 units,STATE)
         ! toggle output off for observation
!geh: nope, this can change over time.
!geh         output_variable%plot_only = PETSC_TRUE 

         output_variable%iformat = 1 ! integer
         call OutputVariableAddToList(output_variable_list,output_variable)
         nullify(output_variable)
      case ('TEMPERATURE')
        name = 'Temperature'
        units = 'C'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     TEMPERATURE)
      case ('RESIDUAL')
        units = ''
        do temp_int = 1, option%nflowdof
          write(word,*) temp_int
          name = 'Residual_' // trim(adjustl(word))
          call OutputVariableAddToList(output_variable_list,name, &
                                       OUTPUT_GENERIC,units, &
                                       RESIDUAL,temp_int)
        enddo
      case ('POROSITY')
        units = ''
        name = 'Porosity'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     POROSITY)
      case ('MINERAL_POROSITY')
        units = ''
        name = 'Mineral Porosity'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     MINERAL_POROSITY)
      case ('EFFECTIVE_POROSITY')
        units = ''
        name = 'Effective Porosity'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     EFFECTIVE_POROSITY)
      case ('TORTUOSITY')
        units = ''
        name = 'Tortuosity'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     TORTUOSITY)
      case ('PERMEABILITY','PERMEABILITY_X')
        units = 'm^2'
        name = 'Permeability X'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     PERMEABILITY)
      case ('PERMEABILITY_Y')
        units = 'm^2'
        name = 'Permeability Y'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     PERMEABILITY_Y)
      case ('PERMEABILITY_Z')
        units = 'm^2'
        name = 'Permeability Z'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     PERMEABILITY_Z)
      case ('GAS_PERMEABILITY','GAS_PERMEABILITY_X')
        units = 'm^2'
        name = 'Gas Permeability X'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_PERMEABILITY)
      case ('GAS_PERMEABILITY_Y')
        units = 'm^2'
        name = 'Gas Permeability Y'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_PERMEABILITY_Y)
      case ('GAS_PERMEABILITY_Z')
        units = 'm^2'
        name = 'Gas Permeability Z'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     GAS_PERMEABILITY_Z)
      case ('SOIL_COMPRESSIBILITY')
        units = ''
        name = 'Compressibility'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     SOIL_COMPRESSIBILITY)
      case ('SOIL_REFERENCE_PRESSURE')
        units = 'Pa'
        name = 'Soil Reference Pressure'
        call OutputVariableAddToList(output_variable_list,name, &
                                     OUTPUT_GENERIC,units, &
                                     SOIL_REFERENCE_PRESSURE)
      case ('NATURAL_ID')
        units = ''
        name = 'Natural ID'
        output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                                units,NATURAL_ID)
        output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
        output_variable%iformat = 1 ! integer
        call OutputVariableAddToList(output_variable_list,output_variable)
      case ('PROCESS_ID')
        units = ''
        name = 'Process ID'
        output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                                units,PROCESS_ID)
        output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
        output_variable%iformat = 1 ! integer
        call OutputVariableAddToList(output_variable_list,output_variable)
      case ('VOLUME')
        units = 'm^3'
        name = 'Volume'
        output_variable => OutputVariableCreate(name,OUTPUT_GENERIC, &
                                                units,VOLUME)
        !output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
        output_variable%iformat = 0 ! double
        call OutputVariableAddToList(output_variable_list,output_variable)
      case ('MATERIAL_ID')
        units = ''
        name = 'Material ID'
        output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                                units,MATERIAL_ID)
        output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
        output_variable%iformat = 1 ! integer
        call OutputVariableAddToList(output_variable_list,output_variable)
      case ('FRACTURE')
        units = ''
        name = 'Fracture Flag'
        output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                                units,FRACTURE)
        output_variable%iformat = 1 ! integer
        call OutputVariableAddToList(output_variable_list,output_variable)
      case ('MATERIAL_ID_KLUDGE_FOR_VISIT')
        units = ''
        name = 'Kludged material ids for VisIt'
        output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                                units,MATERIAL_ID)
        output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
        output_variable%iformat = 1 ! integer
        call OutputVariableAddToList(output_variable_list,output_variable)
      case('NO_FLOW_VARIABLES')
        output_variable_list%flow_vars = PETSC_FALSE
      case('NO_ENERGY_VARIABLES')
        output_variable_list%energy_vars = PETSC_FALSE
      case default
        call InputKeywordUnrecognized(word,'VARIABLES',option)
    end select

  enddo

end subroutine OutputVariableRead

! ************************************************************************** !

subroutine Output(realization_base,snapshot_plot_flag,observation_plot_flag, &
                  massbal_plot_flag)
  ! 
  ! Main driver for all output subroutines
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! Notes: Modified by Jenn Frederick, 2/23/2016
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module, only : OptionCheckTouch, option_type, printMsg, printErrMsg
  
  implicit none
  
  class(realization_base_type) :: realization_base
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend
  type(option_type), pointer :: option

  option => realization_base%option
  
  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)

  ! check for plot request from active directory
  if (.not.snapshot_plot_flag) then

    if (option%use_touch_options) then
      string = 'plot'
      if (OptionCheckTouch(option,string)) then
        realization_base%output_option%plot_name = 'plot'
        snapshot_plot_flag = PETSC_TRUE
      endif
    endif

  endif

!.................................
  if (snapshot_plot_flag) then

    if (realization_base%output_option%print_hdf5) then
      call PetscTime(tstart,ierr);CHKERRQ(ierr)
      call PetscLogEventBegin(logging%event_output_hdf5,ierr);CHKERRQ(ierr)
      if (realization_base%discretization%itype == UNSTRUCTURED_GRID) then
        select case (realization_base%discretization%grid%itype)
          case (EXPLICIT_UNSTRUCTURED_GRID)
             call OutputHDF5UGridXDMFExplicit(realization_base, &
                  INSTANTANEOUS_VARS)
          case (IMPLICIT_UNSTRUCTURED_GRID)
            call OutputHDF5UGridXDMF(realization_base,INSTANTANEOUS_VARS)
          case (POLYHEDRA_UNSTRUCTURED_GRID)
            call printErrMsg(option,'Add code for HDF5 output for &
                                    &Polyhedra mesh')
        end select
      else
        call OutputHDF5(realization_base,INSTANTANEOUS_VARS)
      endif      
      call PetscLogEventEnd(logging%event_output_hdf5,ierr);CHKERRQ(ierr)
      call PetscTime(tend,ierr);CHKERRQ(ierr)
#ifdef SCORPIO_WRITE
      if (option%myrank == 0) write (*,'(" Parallel IO Write method is used in &
        &writing the output, HDF5_WRITE_GROUP_SIZE = ",i5)') &
        option%hdf5_write_group_size
#endif
      write(option%io_buffer,'(f10.2," Seconds to write HDF5 file.")') &
            tend-tstart
#ifndef CLM_PFLOTRAN
      call printMsg(option)
#endif
    endif

    if (realization_base%output_option%print_explicit_flowrate) then
      call PetscTime(tstart,ierr);CHKERRQ(ierr)
      call PetscLogEventBegin(logging%event_output_tecplot,ierr);CHKERRQ(ierr)
      !call OutputPrintExplicitFlowrates(realization_base)
      call PetscLogEventEnd(logging%event_output_tecplot,ierr);CHKERRQ(ierr)
      call PetscTime(tend,ierr);CHKERRQ(ierr)
      write(option%io_buffer,'(f10.2," Seconds to write to Rates file.")') &
            tend-tstart
      call printMsg(option)        
    endif
      
    if (realization_base%output_option%print_mad) then
      call PetscTime(tstart,ierr);CHKERRQ(ierr)
      call PetscLogEventBegin(logging%event_output_mad,ierr);CHKERRQ(ierr)
      call OutputMAD(realization_base)

      call PetscLogEventEnd(logging%event_output_mad,ierr);CHKERRQ(ierr)
      call PetscTime(tend,ierr);CHKERRQ(ierr)
      write(option%io_buffer,'(f10.2," Seconds to write to MAD HDF5 &
                             &file(s)")') tend-tstart
      call printMsg(option) 
    endif

    if (option%compute_statistics) then
      call ComputeFlowCellVelocityStats(realization_base)
      call ComputeFlowFluxVelocityStats(realization_base)
    endif

  endif
  
!.................................
  if (observation_plot_flag) then
    call OutputObservation(realization_base)
  endif

!.................................
  if (massbal_plot_flag) then
    call OutputMassBalance(realization_base)
  endif

  ! Output temporally average variables 
  call OutputAvegVars(realization_base)

  if (snapshot_plot_flag) then
    realization_base%output_option%plot_number = &
      realization_base%output_option%plot_number + 1
  endif

  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  realization_base%output_option%plot_name = ''

  call PetscLogStagePop(ierr);CHKERRQ(ierr)

end subroutine Output

! ************************************************************************** !

subroutine OutputInputRecord(output_option,waypoint_list)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  !  
  use Output_Aux_module
  use Waypoint_module

  implicit none

  type(output_option_type), pointer :: output_option
  type(waypoint_list_type), pointer :: waypoint_list
  
  type(waypoint_type), pointer :: cur_waypoint
  type(output_variable_type), pointer :: cur_variable
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: snap_string,obs_string,msbl_string
  PetscBool :: snap_output_found,obs_output_found,msbl_output_found
  PetscInt :: id = INPUT_RECORD_UNIT  
  character(len=10) :: Format
  
  Format = '(ES14.7)'

  write(id,'(a)') ' '
    write(id,'(a)') '---------------------------------------------------------&
                    &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'OUTPUT FILES'

  write(id,'(a29)',advance='no') 'periodic screen: '
  if (output_option%screen_imod /= 0) then
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'screen increment: '
    write(word,*) output_option%screen_imod
    write(id,'(a)') adjustl(trim(word))
  else
    write(id,'(a)') 'OFF'
  endif

  write(id,'(a29)',advance='no') 'output time unit: '
  write(id,'(a)') trim(output_option%tunit)

  snap_string = ''
  obs_string = ''
  msbl_string = ''
  snap_output_found = PETSC_FALSE
  obs_output_found = PETSC_FALSE
  msbl_output_found = PETSC_FALSE
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%print_snap_output) then
      snap_output_found = PETSC_TRUE
      write(word,Format) cur_waypoint%time / output_option%tconv
      snap_string = trim(snap_string) // adjustl(trim(word)) // ','
    endif
    if (cur_waypoint%print_obs_output) then
      obs_output_found = PETSC_TRUE
      write(word,Format) cur_waypoint%time / output_option%tconv
      obs_string = trim(obs_string) // adjustl(trim(word)) // ','
    endif
    if (cur_waypoint%print_msbl_output) then
      msbl_output_found = PETSC_TRUE
      write(word,Format) cur_waypoint%time / output_option%tconv
      msbl_string = trim(msbl_string) // adjustl(trim(word)) // ','
    endif
    cur_waypoint => cur_waypoint%next
  enddo

  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'snapshot file output'
  if (output_option%print_hdf5) then
    write(id,'(a29)',advance='no') 'format: '
    if (output_option%print_single_h5_file) then
      write(id,'(a)') 'hd5f, single file'
    endif
    if (output_option%times_per_h5_file /= 1) then
      write(word,*) output_option%times_per_h5_file
      write(id,'(a)') 'hdf5, ' // trim(word) // ' times per file'
    endif
    if (output_option%print_hdf5_vel_cent) then
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') 'velocity on cell centers'
    endif
    if (output_option%print_hdf5_vel_face) then
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') 'velocity on cell faces'
    endif
    if (output_option%print_hdf5_mass_flowrate) then
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') 'mass flow rate'
    endif
    if (output_option%print_hdf5_energy_flowrate) then
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') 'energy flow rate'
    endif
    if (output_option%print_hdf5_aveg_mass_flowrate) then
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') 'average mass flow rate'
    endif
    if (output_option%print_hdf5_aveg_energy_flowrate) then
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') 'average energy flow rate'
    endif
    if (output_option%print_explicit_flowrate) then
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') 'explicit flow rate'
    endif
  endif
  if (output_option%print_mad) then
    write(id,'(a29)',advance='no') 'format: '
    write(id,'(a)') 'mad'
  endif
  if (output_option%print_vtk) then
    write(id,'(a29)',advance='no') 'format: '
    write(id,'(a)') 'vtk'
  endif
  write(id,'(a29)',advance='no') 'periodic timestep: '
  if (output_option%periodic_snap_output_ts_imod == 100000000) then
    write(id,'(a)') 'OFF'
  else
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'timestep increment: '
    write(word,'(i9)') output_option%periodic_snap_output_ts_imod
    write(id,'(a)') adjustl(trim(word))
  endif
  write(id,'(a29)',advance='no') 'periodic time: '
  if (output_option%periodic_snap_output_time_incr <= 0) then
    write(id,'(a)') 'OFF'
  else
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'time increment: '
    write(word,Format) output_option%periodic_snap_output_time_incr / &
                  output_option%tconv
    write(id,'(a)') adjustl(trim(word)) // &
                    adjustl(trim(output_option%tunit))
  endif
  write(id,'(a29)',advance='no') 'specific times: '
  if (snap_output_found) then
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'times (' // &
                                    trim(output_option%tunit) // '): '
    write(id,'(a)') trim(snap_string)
  else
    write(id,'(a)') 'OFF'
  endif
  if (associated(output_option%output_snap_variable_list%first)) then
    write(id,'(a29)',advance='no') 'variable list: '
    cur_variable => output_option%output_snap_variable_list%first
    write(id,'(a)') trim(cur_variable%name) // ' [' // &
                    trim(cur_variable%units) // ']'
    cur_variable => cur_variable%next
    do
      if (.not.associated(cur_variable)) exit
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') trim(cur_variable%name) // ' [' // &
           trim(cur_variable%units) // ']'
      cur_variable => cur_variable%next
    enddo
  endif
  write(id,'(a29)',advance='no') 'print initial time: '
  if (output_option%print_initial_snap) then
    write(id,'(a)') 'ON'
  else
    write(id,'(a)') 'OFF'
  endif
  write(id,'(a29)',advance='no') 'print final time: '
  if (output_option%print_final_snap) then
    write(id,'(a)') 'ON'
  else
    write(id,'(a)') 'OFF'
  endif

  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'observation file output'
  write(id,'(a29)',advance='no') 'format: '
  write(id,'(a)') 'tecplot'
  write(id,'(a29)',advance='no') 'periodic timestep: '
  if (output_option%periodic_obs_output_ts_imod == 100000000) then
    write(id,'(a)') 'OFF'
  else
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'timestep increment: '
    write(word,'(i9)') output_option%periodic_obs_output_ts_imod
    write(id,'(a)') adjustl(trim(word))
  endif
  write(id,'(a29)',advance='no') 'periodic time: '
  if (output_option%periodic_obs_output_time_incr <= 0) then
    write(id,'(a)') 'OFF'
  else
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'time increment: '
    write(word,Format) output_option%periodic_obs_output_time_incr / &
                  output_option%tconv
    write(id,'(a)') adjustl(trim(word)) // &
                    adjustl(trim(output_option%tunit))
  endif
  write(id,'(a29)',advance='no') 'specific times: '
  if (obs_output_found) then
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'times (' // &
                                    trim(output_option%tunit) // '): '
    write(id,'(a)') trim(obs_string)
  else
    write(id,'(a)') 'OFF'
  endif
  if (associated(output_option%output_obs_variable_list%first)) then
    write(id,'(a29)',advance='no') 'variable list: '
    cur_variable => output_option%output_obs_variable_list%first
    write(id,'(a)') trim(cur_variable%name)
    cur_variable => cur_variable%next
    do
      if (.not.associated(cur_variable)) exit
      write(id,'(a29)',advance='no') ' '
      write(id,'(a)') trim(cur_variable%name) // ' [' // &
           trim(cur_variable%units) // ']'
      cur_variable => cur_variable%next
    enddo
  endif
  write(id,'(a29)',advance='no') 'print initial time: '
  if (output_option%print_initial_obs) then
    write(id,'(a)') 'ON'
  else
    write(id,'(a)') 'OFF'
  endif
  write(id,'(a29)',advance='no') 'print final time: '
  if (output_option%print_final_obs) then
    write(id,'(a)') 'ON'
  else
    write(id,'(a)') 'OFF'
  endif

  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'mass balance file output'
  write(id,'(a29)',advance='no') 'format: '
  write(id,'(a)') 'tecplot'
  write(id,'(a29)',advance='no') 'periodic timestep: '
  if (output_option%periodic_msbl_output_ts_imod == 100000000) then
    write(id,'(a)') 'OFF'
  else
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'timestep increment: '
    write(word,'(i7)') output_option%periodic_msbl_output_ts_imod
    write(id,'(a)') adjustl(trim(word))
  endif
  write(id,'(a29)',advance='no') 'periodic time: '
  if (output_option%periodic_msbl_output_time_incr <= 0) then
    write(id,'(a)') 'OFF'
  else
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'time increment: '
    write(word,Format) output_option%periodic_msbl_output_time_incr / &
                  output_option%tconv
    write(id,'(a)') adjustl(trim(word)) // &
                    adjustl(trim(output_option%tunit))
  endif
  write(id,'(a29)',advance='no') 'specific times: '
  if (msbl_output_found) then
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'times (' // &
                                    trim(output_option%tunit) // '): '
    write(id,'(a)') trim(msbl_string)
  else
    write(id,'(a)') 'OFF'
  endif
  write(id,'(a29)',advance='no') 'print initial time: '
  if (output_option%print_initial_massbal) then
    write(id,'(a)') 'ON'
  else
    write(id,'(a)') 'OFF'
  endif
  write(id,'(a29)',advance='no') 'print final time: '
  if (output_option%print_final_massbal) then
    write(id,'(a)') 'ON'
  else
    write(id,'(a)') 'OFF'
  endif
  

end subroutine OutputInputRecord

! ************************************************************************** !

subroutine OutputMAD(realization_base)
  ! 
  ! Print to HDF5 file for MAD final output
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use Realization_Base_class, only : realization_base_type, &
                                     RealizationGetVariable
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Reaction_Aux_module
  use Variables_module
  use Output_Common_module, only : OutputGetVariableArray
 
#if !defined(PETSC_HAVE_HDF5)
  implicit none
  
  class(realization_base_type) :: realization_base

  call printMsg(realization_base%option,'')
  write(realization_base%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to ", &
        &"write HDF5 formatted structured grids.")')
  call printErrMsg(realization_base%option)
#else

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use hdf5
  use HDF5_module
  
  implicit none

  class(realization_base_type) :: realization_base

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer(HSIZE_T) :: dims(3)
  
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  
  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: array(:)
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_flag 
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr

  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  reaction => realization_base%reaction
  output_option => realization_base%output_option

#define ALL
#ifdef ALL
  write(string,'(i6)') option%mygroup_id
  filename = trim(option%global_prefix) // '-MAD-G' // trim(adjustl(string)) // '.h5'
!  filename = trim(option%global_prefix) // '-MAD.h5'

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  ! turn off error reporting
  call h5eset_auto_f(OFF,hdf5_err)
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                     prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)
#else
  filename = trim(option%global_prefix) // trim(option%group_prefix) // '.h5'

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err,H5P_DEFAULT_F, &
                   prop_id)
  call h5pclose_f(prop_id,hdf5_err)
#endif

  ! write out data sets 
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)   

  ! pressure
  call RealizationGetVariable(realization_base,global_vec,LIQUID_PRESSURE, &
                              ZERO_INTEGER)
#ifdef ALL
  string = 'Pressure' // trim(option%group_prefix)
#else
  string = 'Pressure'
#endif
  call HDF5WriteStructDataSetFromVec(string,realization_base,global_vec, &
                                     file_id,H5T_NATIVE_DOUBLE)

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif
end subroutine OutputMAD

! ************************************************************************** !

subroutine ComputeFlowCellVelocityStats(realization_base)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/11/08
  ! 

  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Patch_module
  use Discretization_module

  implicit none
  
  class(realization_base_type) :: realization_base
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(discretization_type), pointer :: discretization
  type(output_option_type), pointer :: output_option
  PetscInt :: iconn, i, direction, iphase, sum_connection
  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscReal :: flux
  Vec :: global_vec, global_vec2

  PetscReal :: average, sum, max, min, std_dev
  PetscInt :: max_loc, min_loc
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:), den_loc_p(:)
  PetscReal, allocatable :: sum_area(:)
  PetscErrorCode :: ierr
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  discretization => realization_base%discretization
    
  allocate(sum_area(grid%nlmax))
  call DiscretizationDuplicateVector(discretization,field%work,global_vec)
  call DiscretizationDuplicateVector(discretization,field%work,global_vec2)

  do iphase = 1,option%nphase

    do direction = 1,3
    
      sum_area(1:grid%nlmax) = 0.d0
      call VecSet(global_vec,0.d0,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)

      ! interior velocities  
      connection_set_list => grid%internal_connection_set_list
      cur_connection_set => connection_set_list%first
      sum_connection = 0
      do 
        if (.not.associated(cur_connection_set)) exit
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          ghosted_id_up = cur_connection_set%id_up(iconn)
          ghosted_id_dn = cur_connection_set%id_dn(iconn)
          local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
          local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
          ! velocities are stored as the downwind face of the upwind cell
          flux = patch%internal_velocities(iphase,sum_connection)* &
                   cur_connection_set%area(iconn)* &
                   cur_connection_set%dist(direction,iconn)
          if (local_id_up > 0) then
            vec_ptr(local_id_up) = vec_ptr(local_id_up) - flux
          endif
          if (local_id_dn > 0) then
            vec_ptr(local_id_dn) = vec_ptr(local_id_dn) + flux
          endif
        enddo
        cur_connection_set => cur_connection_set%next
      enddo

      ! boundary velocities
      boundary_condition => patch%boundary_condition_list%first
      sum_connection = 0
      do
        if (.not.associated(boundary_condition)) exit
        cur_connection_set => boundary_condition%connection_set
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          local_id = cur_connection_set%id_dn(iconn)
          vec_ptr(local_id) = vec_ptr(local_id)+ &
                              cur_connection_set%dist(direction,iconn)* &
                              patch%boundary_velocities(iphase,sum_connection)* &
                              cur_connection_set%area(iconn)
        enddo
        boundary_condition => boundary_condition%next
      enddo

      call VecRestoreArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)

      call VecSum(global_vec,sum,ierr);CHKERRQ(ierr)
      average = sum/real(grid%nmax)
      call VecSet(global_vec2,average,ierr);CHKERRQ(ierr)
      call VecMax(global_vec,max_loc,max,ierr);CHKERRQ(ierr)
      call VecMin(global_vec,min_loc,min,ierr);CHKERRQ(ierr)
      call VecAYPX(global_vec2,-1.d0,global_vec,ierr);CHKERRQ(ierr)
      call VecNorm(global_vec2,NORM_2,std_dev,ierr);CHKERRQ(ierr)
      select case(direction)
        case(X_DIRECTION)
          string = 'X-Direction,'
        case(Y_DIRECTION)
          string = 'Y-Direction,'
        case(Z_DIRECTION)
          string = 'Z-Direction,'
      end select
      select case(iphase)
        case(LIQUID_PHASE)
          string = trim(string) // ' Liquid Phase'
        case(GAS_PHASE)
          string = trim(string) // ' Gas Phase'
      end select
      string = trim(string) // ' Velocity Statistics [m/' // &
               trim(output_option%tunit) // ']:'

      if (option%myrank == option%io_rank) then
        write(*,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
        write(option%fid_out,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
      endif

    enddo
  enddo
  
  if (allocated(sum_area)) deallocate(sum_area)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec2,ierr);CHKERRQ(ierr)

end subroutine ComputeFlowCellVelocityStats

! ************************************************************************** !

subroutine ComputeFlowFluxVelocityStats(realization_base)
  ! 
  ! Print flux statistics
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/11/08
  ! 
!geh - specifically, the flow velocities at the interfaces between cells
 
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Grid_module
  use Option_module
  use Field_module
  use Connection_module
  use Patch_module
  
  implicit none

  class(realization_base_type) :: realization_base
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: discretization  
  type(output_option_type), pointer :: output_option
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscInt :: iphase
  PetscInt :: direction
  PetscInt :: local_id, ghosted_id
  PetscInt :: iconn, sum_connection
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, global_vec2
  PetscReal :: sum, average, max, min , std_dev
  PetscInt :: max_loc, min_loc
  PetscErrorCode :: ierr

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
    
  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  field => realization_base%field
  output_option => realization_base%output_option
  
  call DiscretizationDuplicateVector(discretization,field%work,global_vec) 
  call DiscretizationDuplicateVector(discretization,field%work,global_vec2) 

  do iphase = 1,option%nphase
    do direction = 1,3
    
      call VecZeroEntries(global_vec,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
      
      ! place interior velocities in a vector
      connection_set_list => grid%internal_connection_set_list
      cur_connection_set => connection_set_list%first
      sum_connection = 0
      do 
        if (.not.associated(cur_connection_set)) exit
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          ghosted_id = cur_connection_set%id_up(iconn)
          local_id = grid%nG2L(ghosted_id) ! = zero for ghost nodes
          ! velocities are stored as the downwind face of the upwind cell
          if (local_id <= 0 .or. &
              dabs(cur_connection_set%dist(direction,iconn)) < 0.99d0) cycle
          vec_ptr(local_id) = patch%internal_velocities(iphase,sum_connection)
        enddo
        cur_connection_set => cur_connection_set%next
      enddo

      call VecRestoreArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
      
      ! compute stats
      call VecSum(global_vec,sum,ierr);CHKERRQ(ierr)
      average = sum/real(grid%nmax)
      call VecSet(global_vec2,average,ierr);CHKERRQ(ierr)
      call VecMax(global_vec,max_loc,max,ierr);CHKERRQ(ierr)
      call VecMin(global_vec,min_loc,min,ierr);CHKERRQ(ierr)
      call VecAYPX(global_vec2,-1.d0,global_vec,ierr);CHKERRQ(ierr)
      call VecNorm(global_vec2,NORM_2,std_dev,ierr);CHKERRQ(ierr)
      select case(direction)
        case(X_DIRECTION)
          string = 'X-Direction,'
        case(Y_DIRECTION)
          string = 'Y-Direction,'
        case(Z_DIRECTION)
          string = 'Z-Direction,'
      end select
      select case(iphase)
        case(LIQUID_PHASE)
          string = trim(string) // ' Liquid Phase'
        case(GAS_PHASE)
          string = trim(string) // ' Gas Phase'
      end select
      string = trim(string) // ' Flux Velocity Statistics [m/' // &
               trim(output_option%tunit) // ']:'
      if (option%myrank == option%io_rank) then
        write(*,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
        write(option%fid_out,'(/,a,/, &
                     &"Average:",1es12.4,/, &
                     &"Max:    ",1es12.4,"  Location:",i11,/, &
                     &"Min:    ",1es12.4,"  Location:",i11,/, &
                     &"Std Dev:",1es12.4,/)') trim(string), &
                                              average,max,max_loc+1, &
                                              min,min_loc+1,std_dev
      endif
    enddo
  enddo
  
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec2,ierr);CHKERRQ(ierr)
  
end subroutine ComputeFlowFluxVelocityStats

! ************************************************************************** !

subroutine OutputPrintCouplers(realization_base,istep)
  ! 
  ! Prints values of auxiliary variables associated with
  ! couplers (boundary and initial conditions, source
  ! sinks).  Note that since multiple connections for
  ! couplers can exist for a single cell, the latter will
  ! overwrite the former.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/11
  ! 

  use Realization_Base_class, only : realization_base_type
  use Coupler_module
  use Connection_module
  use Option_module
  use Debug_module
  use Field_module
  use Patch_module
  use Grid_module
  use Input_Aux_module

  class(realization_base_type) :: realization_base
  PetscInt :: istep
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: cur_patch
  type(field_type), pointer :: field
  type(coupler_type), pointer :: coupler
  type(debug_type), pointer :: flow_debug
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string, coupler_string
  type(connection_set_type), pointer :: cur_connection_set
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id, iconn, iaux
  PetscInt, allocatable :: iauxvars(:)
  character(len=MAXWORDLENGTH), allocatable :: auxvar_names(:)
  PetscErrorCode :: ierr
  
  
  option => realization_base%option
  flow_debug => realization_base%debug
  field => realization_base%field

  if (len_trim(flow_debug%coupler_string) == 0) then
    option%io_buffer = &
      'Coupler debugging requested, but no string of coupler names was included.'
    call printErrMsg(option)
  endif

  select case(option%iflowmode)
    case(TH_MODE)
      allocate(iauxvars(2),auxvar_names(2))
      iauxvars(1) = TH_PRESSURE_DOF
      auxvar_names(1) = 'liquid_pressure'
      iauxvars(2) = TH_TEMPERATURE_DOF
      auxvar_names(2) = 'temperature'
    case default
      option%io_buffer = &
        'OutputPrintCouplers() not yet supported for this flow mode'
      call printErrMsg(option)
  end select
  
  coupler_string = flow_debug%coupler_string
  ierr = 0
  do
    call InputReadWord(coupler_string,word,PETSC_TRUE,ierr)
    if (ierr /= 0) exit
    
    do iaux = 1, size(iauxvars)
      cur_patch => realization_base%patch_list%first
      do
        if (.not.associated(cur_patch)) exit
        grid => cur_patch%grid
        coupler => CouplerGetPtrFromList(word,cur_patch%boundary_condition_list)
        call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
        call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        if (associated(coupler)) then
          cur_connection_set => coupler%connection_set
          do iconn = 1, cur_connection_set%num_connections
            local_id = cur_connection_set%id_dn(iconn)
            if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
            vec_ptr(local_id) = coupler%flow_aux_real_var(iauxvars(iaux),iconn)
          enddo
        endif
        call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        cur_patch => cur_patch%next
      enddo

      if (istep > 0) then
        write(string,*) istep
        string = adjustl(string)
        string = trim(word) // '_' // trim(auxvar_names(iaux)) // '_' // &
                 trim(string)
      else 
        string = trim(word) // '_' // trim(auxvar_names(iaux))
      endif
      if (len_trim(option%group_prefix) > 1) then
        string = trim(string) // trim(option%group_prefix)
      endif
      string = trim(string) // '.tec'
      !call OutputVectorTecplot(string,word,realization_base,field%work)
    enddo
      
  enddo

  deallocate(iauxvars)
  deallocate(auxvar_names)

end subroutine OutputPrintCouplers

! ************************************************************************** !

subroutine OutputPrintRegions(realization_base)
  ! 
  ! Prints out the number of connections to each cell in a region.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/03/16
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module
  use Debug_module
  use Field_module
  use Patch_module
  use Region_module

  class(realization_base_type) :: realization_base
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(debug_type), pointer :: flow_debug
  type(region_type), pointer :: cur_region
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr
  
  
  option => realization_base%option
  flow_debug => realization_base%debug
  field => realization_base%field

  cur_region => realization_base%patch%region_list%first
  do
    if (.not.associated(cur_region)) exit
    call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    do i = 1, cur_region%num_cells
      vec_ptr(cur_region%cell_ids(i)) = vec_ptr(cur_region%cell_ids(i)) + 1.d0
    enddo
    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    string = 'region_' // trim(cur_region%name) // '.tec'
    word = 'region'
    !call OutputVectorTecplot(string,word,realization_base,field%work)
    cur_region => cur_region%next
  enddo
  
end subroutine OutputPrintRegions

! ************************************************************************** !

subroutine OutputAvegVars(realization_base)
  ! 
  ! This routine temporally averages variables and outputs thems
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/10/13
  ! 

  use Realization_Base_class, only : realization_base_type
  use Option_module, only : OptionCheckTouch, option_type, printMsg
  use Output_Aux_module
  use Output_Common_module, only : OutputGetVariableArray  
  use Field_module

  implicit none
  
  class(realization_base_type) :: realization_base

  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  type(field_type), pointer :: field  

  PetscReal :: dtime
  PetscBool :: aveg_plot_flag
  PetscInt :: ivar
  PetscReal,pointer :: aval_p(:),ival_p(:)
  PetscErrorCode :: ierr  
  PetscLogDouble :: tstart, tend

  option => realization_base%option
  output_option => realization_base%output_option
  field => realization_base%field

  ! 
  if (option%time<1.d-10) return
  
  dtime = option%time-output_option%aveg_var_time
  output_option%aveg_var_dtime = output_option%aveg_var_dtime + dtime
  output_option%aveg_var_time = output_option%aveg_var_time + dtime
  
  if (abs(output_option%aveg_var_dtime - &
          output_option%periodic_snap_output_time_incr)<1.d0) then
    aveg_plot_flag=PETSC_TRUE
  else
    aveg_plot_flag=PETSC_FALSE
  endif

  if (.not.associated(output_option%aveg_output_variable_list%first)) then
    if (output_option%print_hdf5_aveg_mass_flowrate.or. &
       output_option%print_hdf5_aveg_energy_flowrate) then
      ! There is a possibility to output average-flowrates, thus
      ! call output subroutine depending on mesh type
      if (realization_base%discretization%itype == UNSTRUCTURED_GRID) then
        call OutputHDF5UGridXDMF(realization_base,AVERAGED_VARS)
      else
      !  call OutputHDF5(realization_base,AVERAGED_VARS)
      endif
    endif
    return
  endif
  
  ivar = 0
  cur_variable => output_option%aveg_output_variable_list%first
  do
    if (.not.associated(cur_variable)) exit

    ! Get the variable
    call OutputGetVariableArray(realization_base,field%work,cur_variable)

    ! Cumulatively add the variable*dtime
    ivar = ivar + 1
    call VecGetArrayF90(field%work,ival_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%avg_vars_vec(ivar),aval_p,ierr);CHKERRQ(ierr)
    aval_p = aval_p + ival_p*dtime
    call VecRestoreArrayF90(field%work,ival_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%avg_vars_vec(ivar),aval_p, &
                            ierr);CHKERRQ(ierr)

    ! Check if it is time to output the temporally average variable
    if (aveg_plot_flag) then

      ! Divide vector values by 'time'
      call VecGetArrayF90(field%avg_vars_vec(ivar),aval_p,ierr);CHKERRQ(ierr)
      aval_p = aval_p/output_option%periodic_snap_output_time_incr
      call VecRestoreArrayF90(field%avg_vars_vec(ivar),aval_p, &
                              ierr);CHKERRQ(ierr)

    endif
    
    cur_variable => cur_variable%next
  enddo

  if (aveg_plot_flag) then

    if (realization_base%output_option%print_hdf5) then
      call PetscTime(tstart,ierr);CHKERRQ(ierr)
      call PetscLogEventBegin(logging%event_output_hdf5,ierr);CHKERRQ(ierr)
      if (realization_base%discretization%itype == UNSTRUCTURED_GRID) then
        call OutputHDF5UGridXDMF(realization_base,AVERAGED_VARS)
      else
        call OutputHDF5(realization_base,AVERAGED_VARS)
      endif      
      call PetscLogEventEnd(logging%event_output_hdf5,ierr);CHKERRQ(ierr)
      call PetscTime(tend,ierr);CHKERRQ(ierr)
      write(option%io_buffer,'(f10.2," Seconds to write HDF5 file.")') tend-tstart
#ifndef CLM_PFLOTRAN
      call printMsg(option)
#endif
    endif

    ! Reset the vectors to zero
    do ivar=1,output_option%aveg_output_variable_list%nvars
      call VecSet(field%avg_vars_vec(ivar),0.d0,ierr);CHKERRQ(ierr)
    enddo

    output_option%aveg_var_dtime=0.d0

  endif


end subroutine OutputAvegVars

! ************************************************************************** !

subroutine OutputEnsureVariablesExist(output_option,option)
  ! 
  ! Loop over output variables to ensure that they exist in the simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/17
  ! 
  use Option_module

  implicit none

  type(output_option_type) :: output_option
  type(option_type) :: option

  call OutputListEnsureVariablesExist(output_option%output_variable_list, &
                                      option)
  call OutputListEnsureVariablesExist(output_option%output_snap_variable_list, &
                                      option)
  call OutputListEnsureVariablesExist(output_option%output_obs_variable_list, &
                                      option)
  call OutputListEnsureVariablesExist(output_option%aveg_output_variable_list, &
                                      option)

end subroutine OutputEnsureVariablesExist

! ************************************************************************** !

subroutine OutputListEnsureVariablesExist(output_variable_list,option)
  ! 
  ! Loop over output variables to ensure that they exist in the simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/17
  ! 
  use Option_module
  use Material_Aux_class, only : soil_compressibility_index, &
                                 soil_reference_pressure_index
  use Variables_module

  implicit none

  type(output_variable_list_type), pointer :: output_variable_list
  type(option_type) :: option

  type(output_variable_type), pointer :: cur_variable
  PetscBool :: error_flag
  PetscInt :: error_count

  cur_variable => output_variable_list%first
  error_count =  0
  do
    if (.not.associated(cur_variable)) exit
    error_flag = PETSC_FALSE
    select case(cur_variable%ivar)
      case(SOIL_COMPRESSIBILITY)
        if (soil_compressibility_index == 0) error_flag = PETSC_TRUE
      case(SOIL_REFERENCE_PRESSURE)
        if (soil_reference_pressure_index == 0) error_flag = PETSC_TRUE
    end select
    if (error_flag) then
      error_count = error_count + 1
      if (error_count == 1) then
        if (OptionPrintToScreen(option)) then
          print *
          print *, 'The following OUTPUT VARIABLES are undefined in this &
            &simulation:'
          print *
        endif
      endif
      if (OptionPrintToScreen(option)) then
        print *, '  ' // trim(cur_variable%name)
      endif
    endif
    cur_variable => cur_variable%next
  enddo
  if (error_count > 0) then
    option%io_buffer = 'Simulation was stopped due to undefined output &
                       &variables.'
    call printErrMsg(option)
  endif

end subroutine OutputListEnsureVariablesExist

! ************************************************************************** !

subroutine OutputFindNaNOrInfInVec(vec,grid,option)
  ! 
  ! Reports Infs or NaNs in a vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/08/18
  ! 
  use Grid_module
  use Option_module
!geh: ieee_arithmetic is not yet supported by gfortran 4.x or lower
!  use ieee_arithmetic

  implicit none

  Vec :: vec
  type(grid_type), pointer :: grid
  type(option_type) :: option

  PetscReal, pointer :: vec_p(:)
  PetscInt :: i, idof, icell, block_size, local_size, local_count, exscan_count
  PetscInt, parameter :: max_number_to_print = 10
  PetscInt :: iarray(2,max_number_to_print)
  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr

  iarray = 0
  call VecGetLocalSize(vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetBlockSize(vec,block_size,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(vec,vec_p,ierr);CHKERRQ(ierr)
  local_count = 0
  do i = 1, local_size
     if (PetscIsInfOrNanReal(vec_p(i))) then
!    if (ieee_is_nan(vec_p(i)) .or. .not.ieee_is_finite(vec_p(i))) then
      local_count = local_count + 1
      icell = int(float(i-1)/float(block_size))+1
      iarray(1,local_count) = grid%nG2A(grid%nL2G(icell))
      idof = i-(icell-1)*block_size
!      if (ieee_is_nan(vec_p(i))) idof = -idof
      iarray(2,local_count) = idof
    endif
  enddo
  call VecRestoreArrayReadF90(vec,vec_p,ierr);CHKERRQ(ierr)

  exscan_count = 0
  call MPI_Exscan(local_count,exscan_count,ONE_INTEGER_MPI, &
                MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  do i = 1, min(max_number_to_print-exscan_count,local_count)
    idof = iarray(2,i)
    if (idof > 0) then
      option%io_buffer = 'NaN'
    else
      option%io_buffer = 'Inf'
    endif
    write(word,*) iarray(1,i)
    option%io_buffer = trim(option%io_buffer) // ' at cell ' // &
      trim(adjustl(word)) // ' and dof'
    write(word,*) iabs(idof)
    option%io_buffer = trim(option%io_buffer) // ' ' // &
      trim(adjustl(word)) //  '.'
    call printMsgByRank(option)
  enddo

end subroutine OutputFindNaNOrInfInVec

end module Output_module
