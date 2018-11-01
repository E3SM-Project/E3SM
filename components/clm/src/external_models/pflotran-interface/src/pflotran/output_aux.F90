module Output_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: INSTANTANEOUS_VARS = 1
  PetscInt, parameter, public :: AVERAGED_VARS = 2
  
  PetscInt, parameter, public :: CHECKPOINT_BINARY = 1
  PetscInt, parameter, public :: CHECKPOINT_HDF5 = 2
  PetscInt, parameter, public :: CHECKPOINT_BOTH = 3

  type, public :: checkpoint_option_type
    character(len=MAXWORDLENGTH) :: tunit
    PetscReal :: tconv
    PetscReal :: periodic_time_incr
    PetscInt :: periodic_ts_incr
    PetscInt :: format
  end type checkpoint_option_type
  
  type, public :: output_option_type

    character(len=MAXWORDLENGTH) :: tunit
    PetscReal :: tconv

    PetscBool :: print_initial_obs
    PetscBool :: print_final_obs
    PetscBool :: print_initial_snap
    PetscBool :: print_final_snap
    PetscBool :: print_initial_massbal
    PetscBool :: print_final_massbal
  
    PetscBool :: print_hdf5
    PetscBool :: extend_hdf5_time_format
    PetscBool :: print_hdf5_vel_cent
    PetscBool :: print_hdf5_vel_face
    PetscBool :: print_single_h5_file
    PetscInt :: times_per_h5_file
    PetscBool :: print_hdf5_mass_flowrate
    PetscBool :: print_hdf5_energy_flowrate
    PetscBool :: print_hdf5_aveg_mass_flowrate
    PetscBool :: print_hdf5_aveg_energy_flowrate
    PetscBool :: print_explicit_flowrate

    PetscBool :: print_tecplot 
    PetscInt :: tecplot_format
    PetscBool :: print_tecplot_vel_cent
    PetscBool :: print_tecplot_vel_face
    PetscBool :: print_fluxes
    
    PetscBool :: print_vtk 
    PetscBool :: print_vtk_vel_cent

    PetscBool :: print_observation 
    PetscBool :: print_column_ids

    PetscBool :: print_mad 

    PetscBool :: print_explicit_primal_grid    ! prints primal grid if true
    PetscBool :: print_explicit_dual_grid      ! prints voronoi (dual) grid if true

    PetscInt :: screen_imod
    PetscInt :: output_file_imod
    
    PetscInt :: periodic_snap_output_ts_imod
    PetscInt :: periodic_obs_output_ts_imod
    PetscInt :: periodic_msbl_output_ts_imod
    
    PetscReal :: periodic_snap_output_time_incr
    PetscReal :: periodic_obs_output_time_incr
    PetscReal :: periodic_msbl_output_time_incr
    
    PetscBool :: filter_non_state_variables

    PetscInt :: xmf_vert_len
    
    type(output_variable_list_type), pointer :: output_variable_list ! (master)
    type(output_variable_list_type), pointer :: output_snap_variable_list
    type(output_variable_list_type), pointer :: output_obs_variable_list
    type(output_variable_list_type), pointer :: aveg_output_variable_list
    
    type(mass_balance_region_type), pointer :: mass_balance_region_list
    PetscBool :: mass_balance_region_flag

    PetscReal :: aveg_var_time
    PetscReal :: aveg_var_dtime
    
    PetscInt :: plot_number
    character(len=MAXWORDLENGTH) :: plot_name

    PetscBool :: print_hydrograph
    PetscInt :: surf_xmf_vert_len

  end type output_option_type
  
  type, public :: output_variable_list_type
    type(output_variable_type), pointer :: first
    type(output_variable_type), pointer :: last
    PetscInt :: nvars
    PetscBool :: flow_vars
    PetscBool :: energy_vars
  end type output_variable_list_type
  
  type, public :: output_variable_type
    character(len=MAXWORDLENGTH) :: name   ! string that appears in hdf5 file
    character(len=MAXWORDLENGTH) :: units
    ! jmf: change to snapshot_plot_only?
    PetscBool :: plot_only
    PetscInt :: iformat   ! 0 = for REAL values; 1 = for INTEGER values
    PetscInt :: icategory ! category for variable-specific regression testing
    PetscInt :: ivar
    PetscInt :: isubvar
    PetscInt :: isubsubvar
    type(output_variable_type), pointer :: next
  end type output_variable_type
  
  type, public :: mass_balance_region_type
    character(len=MAXWORDLENGTH) :: region_name
    PetscInt :: num_cells
    PetscInt, pointer :: region_cell_ids(:)
    PetscReal :: total_mass
    type(mass_balance_region_type), pointer :: next
  end type mass_balance_region_type

!  type, public, EXTENDS (output_variable_type) :: aveg_output_variable_type
!    PetscReal :: time_interval
!  end type aveg_output_variable_type
  
  interface OutputVariableCreate
    module procedure OutputVariableCreate1
    module procedure OutputVariableCreate2
    module procedure OutputVariableCreate3
  end interface OutputVariableCreate
  
  interface OutputVariableAddToList
    module procedure OutputVariableAddToList1
    module procedure OutputVariableAddToList2
  end interface OutputVariableAddToList
  
  ! Output categories
  PetscInt, parameter, public :: OUTPUT_GENERIC = 0
  PetscInt, parameter, public :: OUTPUT_PRESSURE = 1
  PetscInt, parameter, public :: OUTPUT_SATURATION = 2
  PetscInt, parameter, public :: OUTPUT_CONCENTRATION = 3
  PetscInt, parameter, public :: OUTPUT_RATE = 4
  PetscInt, parameter, public :: OUTPUT_VOLUME_FRACTION = 5
  PetscInt, parameter, public :: OUTPUT_DISCRETE = 6
  PetscInt, parameter, public :: OUTPUT_DISPLACEMENT = 7
  PetscInt, parameter, public :: OUTPUT_STRESS = 8
  PetscInt, parameter, public :: OUTPUT_STRAIN = 9

  public :: OutputOptionCreate, &
            OutputOptionDuplicate, &
            OutputVariableCreate, &
            OutputVariableInit, &
            OutputMassBalRegionCreate, &
            OutputVariableListCreate, &
            OutputVariableListDuplicate, &
            OutputMassBalRegListDuplicate, &
            OutputVariableAddToList, &
            OutputWriteToHeader, &
            OutputWriteVariableListToHeader, &
            OutputVariableToCategoryString, &
            OutputVariableAppendDefaults, &
            OpenAndWriteInputRecord, &
            OutputOptionDestroy, &
            OutputVariableListDestroy, &
            CheckpointOptionCreate, &
            CheckpointOptionDestroy

contains

! ************************************************************************** !

function OutputOptionCreate()
  ! 
  ! Creates output options object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(output_option_type), pointer :: OutputOptionCreate

  type(output_option_type), pointer :: output_option
  
  allocate(output_option)
  output_option%print_hdf5 = PETSC_FALSE
  output_option%extend_hdf5_time_format = PETSC_FALSE
  output_option%print_hdf5_vel_cent = PETSC_FALSE
  output_option%print_hdf5_vel_face = PETSC_FALSE
  output_option%print_single_h5_file = PETSC_TRUE
  output_option%times_per_h5_file = 0
  output_option%print_hdf5_mass_flowrate = PETSC_FALSE
  output_option%print_hdf5_energy_flowrate = PETSC_FALSE
  output_option%print_hdf5_aveg_mass_flowrate = PETSC_FALSE
  output_option%print_hdf5_aveg_energy_flowrate = PETSC_FALSE
  output_option%print_explicit_flowrate = PETSC_FALSE
  output_option%print_tecplot = PETSC_FALSE
  output_option%tecplot_format = 0
  output_option%print_tecplot_vel_cent = PETSC_FALSE
  output_option%print_fluxes = PETSC_FALSE
  output_option%print_tecplot_vel_face = PETSC_FALSE
  output_option%print_vtk = PETSC_FALSE
  output_option%print_vtk_vel_cent = PETSC_FALSE
  output_option%print_observation = PETSC_FALSE
  output_option%print_column_ids = PETSC_FALSE
  output_option%print_mad = PETSC_FALSE
  output_option%print_explicit_primal_grid = PETSC_FALSE
  output_option%print_explicit_dual_grid = PETSC_FALSE
  output_option%print_initial_obs = PETSC_TRUE
  output_option%print_final_obs = PETSC_TRUE
  output_option%print_initial_snap = PETSC_TRUE
  output_option%print_final_snap = PETSC_TRUE
  output_option%print_initial_massbal = PETSC_FALSE
  output_option%print_final_massbal = PETSC_TRUE
  output_option%plot_number = 0
  output_option%screen_imod = 1
  output_option%output_file_imod = 1
  output_option%periodic_snap_output_ts_imod  = 100000000
  output_option%periodic_obs_output_ts_imod  = 100000000
  output_option%periodic_msbl_output_ts_imod  = 100000000
  output_option%periodic_snap_output_time_incr = 0
  output_option%periodic_obs_output_time_incr = 0
  output_option%periodic_msbl_output_time_incr = 0
  output_option%plot_name = ""
  output_option%aveg_var_time = 0.d0
  output_option%aveg_var_dtime = 0.d0
  output_option%xmf_vert_len = UNINITIALIZED_INTEGER
  output_option%filter_non_state_variables = PETSC_TRUE

  nullify(output_option%output_variable_list) ! master
  output_option%output_variable_list => OutputVariableListCreate() ! master
  nullify(output_option%output_snap_variable_list)
  output_option%output_snap_variable_list => OutputVariableListCreate()
  nullify(output_option%output_obs_variable_list)
  output_option%output_obs_variable_list => OutputVariableListCreate()
  nullify(output_option%aveg_output_variable_list)
  output_option%aveg_output_variable_list => OutputVariableListCreate()
  
  nullify(output_option%mass_balance_region_list)
  output_option%mass_balance_region_flag = PETSC_FALSE
  
  output_option%tconv = 1.d0
  output_option%tunit = ''
  
  output_option%print_hydrograph = PETSC_FALSE

  OutputOptionCreate => output_option
  
end function OutputOptionCreate

! ************************************************************************** !

function OutputOptionDuplicate(output_option)
  ! 
  ! Creates a copy of output options object
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/2016
  ! 

  implicit none
  
  type(output_option_type), pointer :: output_option

  type(output_option_type), pointer :: OutputOptionDuplicate

  type(output_option_type), pointer :: output_option2
  
  allocate(output_option2)

  output_option2%print_hdf5 = output_option%print_hdf5
  output_option2%extend_hdf5_time_format = &
    output_option%extend_hdf5_time_format
  output_option2%print_hdf5_vel_cent = output_option%print_hdf5_vel_cent
  output_option2%print_hdf5_vel_face = output_option%print_hdf5_vel_face
  output_option2%print_single_h5_file = output_option%print_single_h5_file
  output_option2%times_per_h5_file = output_option%times_per_h5_file
  output_option2%print_hdf5_mass_flowrate = &
    output_option%print_hdf5_mass_flowrate
  output_option2%print_hdf5_energy_flowrate = &
    output_option%print_hdf5_energy_flowrate
  output_option2%print_hdf5_aveg_mass_flowrate = &
    output_option%print_hdf5_aveg_mass_flowrate
  output_option2%print_hdf5_aveg_energy_flowrate = &
    output_option%print_hdf5_aveg_energy_flowrate
  output_option2%print_explicit_flowrate = &
    output_option%print_explicit_flowrate
  output_option2%print_tecplot = output_option%print_tecplot
  output_option2%tecplot_format = output_option%tecplot_format
  output_option2%print_tecplot_vel_cent = output_option%print_tecplot_vel_cent
  output_option2%print_fluxes = output_option%print_fluxes
  output_option2%print_tecplot_vel_face = output_option%print_tecplot_vel_face
  output_option2%print_vtk = output_option%print_vtk
  output_option2%print_vtk_vel_cent = output_option%print_vtk_vel_cent
  output_option2%print_observation = output_option%print_observation
  output_option2%print_column_ids = output_option%print_column_ids
  output_option2%print_mad = output_option%print_mad
  output_option2%print_initial_obs = output_option%print_initial_obs
  output_option2%print_final_obs = output_option%print_final_obs
  output_option2%print_initial_snap = output_option%print_initial_snap
  output_option2%print_final_snap = output_option%print_final_snap
  output_option2%print_initial_massbal = output_option%print_initial_massbal
  output_option2%print_final_massbal = output_option%print_final_massbal
  output_option2%plot_number = output_option%plot_number
  output_option2%screen_imod = output_option%screen_imod
  output_option2%output_file_imod = output_option%output_file_imod
  output_option2%periodic_snap_output_ts_imod = &
    output_option%periodic_snap_output_ts_imod
  output_option2%periodic_obs_output_ts_imod = &
    output_option%periodic_obs_output_ts_imod
  output_option2%periodic_msbl_output_ts_imod = &
    output_option%periodic_msbl_output_ts_imod
  output_option2%periodic_snap_output_time_incr = &
    output_option%periodic_snap_output_time_incr
  output_option2%periodic_obs_output_time_incr = &
    output_option%periodic_obs_output_time_incr
  output_option2%periodic_msbl_output_time_incr = &
    output_option%periodic_msbl_output_time_incr
  output_option2%plot_name = output_option%plot_name
  output_option2%aveg_var_time = output_option%aveg_var_time
  output_option2%aveg_var_dtime = output_option%aveg_var_dtime
  output_option2%xmf_vert_len = output_option%xmf_vert_len
  output_option2%filter_non_state_variables = &
    output_option%filter_non_state_variables

  nullify(output_option2%output_variable_list)
  nullify(output_option2%output_snap_variable_list)
  nullify(output_option2%output_obs_variable_list)
  nullify(output_option2%aveg_output_variable_list)
  
  output_option2%output_variable_list => &
       OutputVariableListDuplicate(output_option%output_variable_list)
  output_option2%output_snap_variable_list => &
       OutputVariableListDuplicate(output_option%output_snap_variable_list)
  output_option2%output_obs_variable_list => &
       OutputVariableListDuplicate(output_option%output_obs_variable_list)
  output_option2%aveg_output_variable_list => &
       OutputVariableListDuplicate(output_option%aveg_output_variable_list)
       
  nullify(output_option2%mass_balance_region_list)
  if (associated(output_option%mass_balance_region_list)) then
    output_option2%mass_balance_region_list => &
       OutputMassBalRegListDuplicate(output_option%mass_balance_region_list)
  endif
  output_option2%mass_balance_region_flag = &
    output_option%mass_balance_region_flag
  
  output_option2%tconv = output_option%tconv
  output_option2%tunit = output_option%tunit
  
  output_option2%print_hydrograph = output_option%print_hydrograph

  OutputOptionDuplicate => output_option2
  
end function OutputOptionDuplicate

! ************************************************************************** !

function CheckpointOptionCreate()
  ! 
  ! Creates output options object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(checkpoint_option_type), pointer :: CheckpointOptionCreate

  type(checkpoint_option_type), pointer :: checkpoint_option
  
  allocate(checkpoint_option)
  checkpoint_option%tunit = ''
  checkpoint_option%tconv = 0.d0
  checkpoint_option%periodic_time_incr = UNINITIALIZED_DOUBLE
  checkpoint_option%periodic_ts_incr = 0
  !checkpoint_option%periodic_ts_incr = huge(checkpoint_option%periodic_ts_incr)
  checkpoint_option%format = CHECKPOINT_BINARY

  CheckpointOptionCreate => checkpoint_option
  
end function CheckpointOptionCreate 
  
! ************************************************************************** !

function OutputVariableCreate1()
  ! 
  ! initializes output variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_type), pointer :: OutputVariableCreate1
  
  type(output_variable_type), pointer :: output_variable
  
  allocate(output_variable)
  call OutputVariableInit(output_variable)
  
  OutputVariableCreate1 => output_variable
  
end function OutputVariableCreate1

! ************************************************************************** !

function OutputVariableCreate2(name,icategory,units,ivar,isubvar,isubsubvar)
  ! 
  ! initializes output variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  character(len=*) :: name
  PetscInt :: icategory ! note that I tuck it inbetween the strings to avoid
                        ! errors
  character(len=*) :: units
  PetscInt :: ivar
  PetscInt, intent(in), optional :: isubvar
  PetscInt, intent(in), optional :: isubsubvar

  type(output_variable_type), pointer :: OutputVariableCreate2
  
  type(output_variable_type), pointer :: output_variable
  
  output_variable => OutputVariableCreate()
  output_variable%name = trim(adjustl(name))
  output_variable%icategory = icategory
  output_variable%units = trim(adjustl(units))
  output_variable%ivar = ivar
  if (present(isubvar)) then
    output_variable%isubvar = isubvar
  endif
  if (present(isubsubvar)) then
    output_variable%isubsubvar = isubsubvar
  endif
  nullify(output_variable%next)
  
  OutputVariableCreate2 => output_variable
  
end function OutputVariableCreate2

! ************************************************************************** !

function OutputVariableCreate3(output_variable)
  ! 
  ! initializes output variable object from an existing
  ! output variabl object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_type), pointer :: output_variable

  type(output_variable_type), pointer :: OutputVariableCreate3
  
  type(output_variable_type), pointer :: new_output_variable
  
  new_output_variable => OutputVariableCreate()
  new_output_variable%name = output_variable%name
  new_output_variable%units = output_variable%units
  new_output_variable%plot_only = output_variable%plot_only
  new_output_variable%iformat = output_variable%iformat
  new_output_variable%icategory = output_variable%icategory
  new_output_variable%ivar = output_variable%ivar
  new_output_variable%isubvar = output_variable%isubvar
  new_output_variable%isubsubvar = output_variable%isubsubvar
  nullify(new_output_variable%next)
  
  OutputVariableCreate3 => new_output_variable
  
end function OutputVariableCreate3

! ************************************************************************** !

subroutine OutputVariableInit(output_variable)
  ! 
  ! initializes output variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/17
  ! 
  implicit none
  
  type(output_variable_type) :: output_variable
  
  output_variable%name = ''
  output_variable%units = ''
  output_variable%plot_only = PETSC_FALSE
  output_variable%iformat = 0
  output_variable%icategory = OUTPUT_GENERIC
  output_variable%ivar = 0
  output_variable%isubvar = 0
  output_variable%isubsubvar = 0
  nullify(output_variable%next)
  
end subroutine OutputVariableInit

! ************************************************************************** !

function OutputVariableListCreate()
  ! 
  ! initializes output variable list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type), pointer :: OutputVariableListCreate
  
  type(output_variable_list_type), pointer :: output_variable_list
  
  allocate(output_variable_list)
  nullify(output_variable_list%first)
  nullify(output_variable_list%last)
  output_variable_list%nvars = 0
  output_variable_list%flow_vars = PETSC_TRUE
  output_variable_list%energy_vars = PETSC_TRUE
  
  OutputVariableListCreate => output_variable_list
  
end function OutputVariableListCreate

! ************************************************************************** !

function OutputMassBalRegionCreate()
  ! 
  ! Creates and initializes a mass balance region list object
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/26/2016
  ! 

  implicit none
  
  type(mass_balance_region_type), pointer :: OutputMassBalRegionCreate
   
  allocate(OutputMassBalRegionCreate)
  OutputMassBalRegionCreate%region_name =''
  nullify(OutputMassBalRegionCreate%region_cell_ids)
  OutputMassBalRegionCreate%num_cells = 0
  OutputMassBalRegionCreate%total_mass = 0.d0
  nullify(OutputMassBalRegionCreate%next)
  
end function OutputMassBalRegionCreate

! ************************************************************************** !

function OutputVariableListDuplicate(old_list)
  ! 
  ! initializes output variable list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type) :: old_list
  
  type(output_variable_list_type), pointer :: OutputVariableListDuplicate
  
  type(output_variable_list_type), pointer :: new_list
  type(output_variable_type), pointer :: cur_variable
  
  allocate(new_list)
  nullify(new_list%first)
  nullify(new_list%last)
  new_list%nvars = old_list%nvars
  
  cur_variable => old_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputVariableAddToList(new_list,OutputVariableCreate(cur_variable))
    cur_variable => cur_variable%next
  enddo

  OutputVariableListDuplicate => new_list
  
end function OutputVariableListDuplicate

! ************************************************************************** !

function OutputMassBalRegListDuplicate(old_list)
  ! 
  ! Duplicates a mass balance region list object
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/27/2016
  ! 

  implicit none
  
  type(mass_balance_region_type), pointer :: old_list
  
  type(mass_balance_region_type), pointer :: new_list
  type(mass_balance_region_type), pointer :: new_mbr
  type(mass_balance_region_type), pointer :: cur_mbr
  type(mass_balance_region_type), pointer :: OutputMassBalRegListDuplicate
  PetscBool :: added
  
  nullify(new_list)

  do
    if (.not.associated(old_list)) exit
    new_mbr => OutputMassBalRegionCreate()
    new_mbr%region_name = old_list%region_name
    new_mbr%num_cells = old_list%num_cells
    new_mbr%region_cell_ids => old_list%region_cell_ids
    new_mbr%total_mass = old_list%total_mass
    ! Add new mass balance region to new list
    if (.not.associated(new_list)) then
      new_list => new_mbr
    else
      cur_mbr => new_list
      do
        if (.not.associated(cur_mbr)) exit
        if (.not.associated(cur_mbr%next)) then
          cur_mbr%next => new_mbr
          added = PETSC_TRUE
        endif
        if (added) exit
        cur_mbr => cur_mbr%next
      enddo
    endif
    old_list => old_list%next
    nullify(new_mbr)
  enddo

  OutputMassBalRegListDuplicate => new_list
  
end function OutputMassBalRegListDuplicate

! ************************************************************************** !

subroutine OutputVariableAddToList1(list,variable)
  ! 
  ! adds variable to list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type) :: list
  type(output_variable_type), pointer :: variable
  
  if (.not. associated(list%first)) then
    list%first => variable
  else
    list%last%next => variable
  endif
  list%last => variable
  
  list%nvars = list%nvars+1
  
end subroutine OutputVariableAddToList1

! ************************************************************************** !

subroutine OutputVariableAddToList2(list,name,icategory,units,ivar, &
                                    isubvar,isubsubvar)
  ! 
  ! creates variable and adds to list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type) :: list
  character(len=*) :: name
  character(len=*) :: units
  PetscInt :: icategory
  PetscInt :: ivar
  PetscInt, intent(in), optional :: isubvar
  PetscInt, intent(in), optional :: isubsubvar
  
  type(output_variable_type), pointer :: variable
  
  if (present(isubvar)) then
    if (present(isubsubvar)) then
      variable => OutputVariableCreate(name,icategory,units, &
                                       ivar,isubvar,isubsubvar)
    else
      variable => OutputVariableCreate(name,icategory,units, &
                                       ivar,isubvar)
    endif
  else
    variable => OutputVariableCreate(name,icategory,units,ivar)
  endif
  call OutputVariableAddToList1(list,variable)
  
end subroutine OutputVariableAddToList2

! ************************************************************************** !

subroutine OutputWriteVariableListToHeader(fid,variable_list,cell_string, &
                                           icolumn,plot_file,variable_count)
  ! 
  ! Converts a variable list to a header string
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  use Option_module
  
  implicit none
  
  PetscInt :: fid
  type(output_variable_list_type) :: variable_list
  character(len=*) :: cell_string
  PetscInt :: icolumn
  PetscBool :: plot_file
  PetscInt :: variable_count
  
  type(output_variable_type), pointer :: cur_variable
  character(len=MAXWORDLENGTH) :: variable_name, units
  
  variable_count = 0
  cur_variable => variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    if (.not. plot_file .and. cur_variable%plot_only) then
      cur_variable => cur_variable%next
      cycle
    endif
    variable_name = cur_variable%name
    units = cur_variable%units
    call OutputWriteToHeader(fid,variable_name,units,cell_string,icolumn)
    variable_count = variable_count + 1
    cur_variable => cur_variable%next
  enddo

end subroutine OutputWriteVariableListToHeader

! ************************************************************************** !

subroutine OutputWriteToHeader(fid,variable_string,units_string, &
                               cell_string, icolumn)
  ! 
  ! Appends formatted strings to header string
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/27/11
  ! 

  implicit none

  PetscInt :: fid
  character(len=*) :: variable_string, units_string, cell_string
  character(len=MAXWORDLENGTH) :: column_string
  character(len=MAXWORDLENGTH) :: variable_string_adj, units_string_adj
  character(len=MAXSTRINGLENGTH) :: cell_string_adj
  PetscInt :: icolumn, len_cell_string, len_units

  character(len=MAXSTRINGLENGTH) :: string

  variable_string_adj = variable_string
  units_string_adj = units_string
  cell_string_adj = cell_string

  !geh: Shift to left.  Cannot perform on same string since len=*
  variable_string_adj = adjustl(variable_string_adj)
  units_string_adj = adjustl(units_string_adj)
  cell_string_adj = adjustl(cell_string_adj)

  if (icolumn > 0) then
    icolumn = icolumn + 1
    write(column_string,'(i4,''-'')') icolumn
    column_string = trim(adjustl(column_string))
  else
    column_string = ''
  endif

  !geh: this is all to remove the lousy spaces
  len_units = len_trim(units_string)
  len_cell_string = len_trim(cell_string)
  if (len_units > 0 .and. len_cell_string > 0) then
    write(string,'('',"'',a,a,'' ['',a,''] '',a,''"'')') trim(column_string), &
          trim(variable_string_adj), trim(units_string_adj), &
          trim(cell_string_adj)
  else if (len_units > 0 .or. len_cell_string > 0) then
    if (len_units > 0) then
      write(string,'('',"'',a,a,'' ['',a,'']"'')') trim(column_string), &
            trim(variable_string_adj), trim(units_string_adj)
    else
      write(string,'('',"'',a,a,'' '',a,''"'')') trim(column_string), &
            trim(variable_string_adj), trim(cell_string_adj)
    endif
  else
    write(string,'('',"'',a,a,''"'')') trim(column_string), &
          trim(variable_string_adj)
  endif
  write(fid,'(a)',advance="no") trim(string)

end subroutine OutputWriteToHeader

! ************************************************************************** !

function OutputVariableToCategoryString(icategory)
  ! 
  ! returns a string associated with an
  ! output variable category
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  PetscInt :: icategory
  
  character(len=MAXWORDLENGTH) :: OutputVariableToCategoryString
  
  character(len=MAXWORDLENGTH) :: string
  
  select case(icategory)
    case(OUTPUT_GENERIC)
      string = 'GENERIC'
    case(OUTPUT_PRESSURE)
      string = 'PRESSURE'
    case(OUTPUT_SATURATION)
      string = 'SATURATION'
    case(OUTPUT_CONCENTRATION)
      string = 'CONCENTRATION'
    case(OUTPUT_RATE)
      string = 'RATE'
    case(OUTPUT_VOLUME_FRACTION)
      string = 'VOLUME_FRACTION'
    case(OUTPUT_DISCRETE)
      string = 'DISCRETE'
    case(OUTPUT_DISPLACEMENT)
      string = 'DISPLACEMENT'
    case(OUTPUT_STRESS)
      string = 'STRESS'
    case(OUTPUT_STRAIN)
      string = 'STRAIN'
    case default
      string = 'GENERIC'
  end select

  OutputVariableToCategoryString = string

end function OutputVariableToCategoryString

! ************************************************************************** !

subroutine OutputVariableAppendDefaults(output_variable_list,option)
  ! 
  ! Adds default output variables to list
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/21/12
  ! 

  use Option_module
  use Variables_module

  implicit none

  type(output_variable_list_type), pointer :: output_variable_list
  type(option_type), pointer :: option
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_type), pointer :: output_variable

  ! Material IDs
  units = ''
  name = 'Material ID'
  output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                          units,MATERIAL_ID)
  output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList(output_variable_list,output_variable)
  
end subroutine OutputVariableAppendDefaults

! ************************************************************************** !

subroutine OpenAndWriteInputRecord(option)
  ! 
  ! Opens the input record file and begins to write to it.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 

  use Option_module

  implicit none
  
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=8) :: date_word
  character(len=10) :: time_word
  character(len=5) :: zone_word
  PetscInt :: id

  id = option%fid_inputrecord
  ! the input record file has a .rec extension:
  filename = trim(option%global_prefix) // trim(option%group_prefix) // '.rec'
  open(unit=id,file=filename,action="write",status="replace")
!geh: this call does not work with IBM
!  call fdate(word)
  call date_and_time(date_word,time_word,zone_word)
  if (OptionPrintToFile(option)) then
    write(id,'(a)') '---------------------------------------------------------&
                    &-----------------------'
    write(id,'(a)') '---------------------------------------------------------&
                    &-----------------------'
    write(id,'(a)') ' PFLOTRAN INPUT RECORD    ' // date_word(5:6) // '/' //  &
                    date_word(7:8) // '/' // date_word(1:4) // ' ' //         &
                    time_word(1:2) // ':' // time_word(3:4) // ' (' //        &
                    zone_word(1:3) // ':' // zone_word(4:5) // ' UTC)'
    write(id,'(a)') '---------------------------------------------------------&
                    &-----------------------'
    write(id,'(a)') '---------------------------------------------------------&
                    &-----------------------'
  
    write(id,'(a18)',advance='no') 'input file: '  
    write(id,*) trim(option%global_prefix) // '.in' 
    
    write(id,'(a18)',advance='no') 'group: ' 
    write(id,*) trim(option%group_prefix)
  
    write(word,*) option%global_commsize
    write(id,'(a18)',advance='no') 'n processors: ' 
    write(id,*) trim(adjustl(word))
  endif

end subroutine OpenAndWriteInputRecord

! ************************************************************************** !

subroutine OutputVariableListDestroy(output_variable_list)
  ! 
  ! Deallocates an output variable list object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_list_type), pointer :: output_variable_list
  
  if (.not.associated(output_variable_list)) return

  nullify(output_variable_list%last)
  call OutputVariableDestroy(output_variable_list%first)
  
  deallocate(output_variable_list)
  nullify(output_variable_list)
  
end subroutine OutputVariableListDestroy

! ************************************************************************** !

recursive subroutine OutputVariableDestroy(output_variable)
  ! 
  ! Deallocates an output variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 

  implicit none
  
  type(output_variable_type), pointer :: output_variable
  
  if (.not.associated(output_variable)) return
  
  call OutputVariableDestroy(output_variable%next)
  
  deallocate(output_variable)
  nullify(output_variable)
  
end subroutine OutputVariableDestroy

! ************************************************************************** !

subroutine CheckpointOptionDestroy(checkpoint_option)
  ! 
  ! Deallocates an output option
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(checkpoint_option_type), pointer :: checkpoint_option
  
  if (.not.associated(checkpoint_option)) return
  
  deallocate(checkpoint_option)
  nullify(checkpoint_option)
  
end subroutine CheckpointOptionDestroy

! ************************************************************************** !

recursive subroutine OutputMassBalRegDestroy(mass_balance_region)
  ! 
  ! Nullifies and deallocates a mass balance region object
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/27/2016
  ! 

  implicit none
  
  type(mass_balance_region_type), pointer :: mass_balance_region
  
  if (associated(mass_balance_region)) then
    ! do not deallocate because the region owns the cell_ids array,
    ! not the mass_balance_region, so just nullify it
    nullify(mass_balance_region%region_cell_ids)
    if (associated(mass_balance_region%next)) then
      call OutputMassBalRegDestroy(mass_balance_region%next)
    endif
    deallocate(mass_balance_region)
  endif
  
end subroutine OutputMassBalRegDestroy

! ************************************************************************** !

subroutine OutputOptionDestroy(output_option)
  ! 
  ! Deallocates an output option
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(output_option_type), pointer :: output_option
  
  if (.not.associated(output_option)) return

  if (associated(output_option%output_variable_list, &
                 output_option%output_snap_variable_list)) then
    nullify(output_option%output_snap_variable_list)
  endif

  if (associated(output_option%output_variable_list, &
                 output_option%output_obs_variable_list)) then
    nullify(output_option%output_obs_variable_list)
  endif
  
  call OutputVariableListDestroy(output_option%output_variable_list)
  call OutputVariableListDestroy(output_option%output_snap_variable_list)
  call OutputVariableListDestroy(output_option%output_obs_variable_list)
  call OutputVariableListDestroy(output_option%aveg_output_variable_list)
  
  call OutputMassBalRegDestroy(output_option%mass_balance_region_list)
    
  deallocate(output_option)
  nullify(output_option)
  
end subroutine OutputOptionDestroy

end module Output_Aux_module
