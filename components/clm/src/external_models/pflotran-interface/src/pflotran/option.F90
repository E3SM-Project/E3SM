module Option_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!


#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Option_Flow_module
  use Option_Transport_module

  implicit none

  private

  type, public :: option_type

    type(flow_option_type), pointer :: flow
    type(transport_option_type), pointer :: transport

    PetscInt :: id                         ! id of realization
    PetscInt :: successful_exit_code       ! code passed out of PFLOTRAN
                                           ! indicating successful completion
                                           ! of simulation
    PetscMPIInt :: global_comm             ! MPI_COMM_WORLD
    PetscMPIInt :: global_rank             ! rank in MPI_COMM_WORLD
    PetscMPIInt :: global_commsize         ! size of MPI_COMM_WORLD
    PetscMPIInt :: global_group            ! id of group for MPI_COMM_WORLD

    PetscMPIInt :: mycomm                  ! PETSC_COMM_WORLD
    PetscMPIInt :: myrank                  ! rank in PETSC_COMM_WORLD
    PetscMPIInt :: mycommsize              ! size of PETSC_COMM_WORLD
    PetscMPIInt :: mygroup                 ! id of group for PETSC_COMM_WORLD
    PetscMPIInt :: mygroup_id

! don't place a character string near here.  It causes the Windows Intel compiler
! to crash.  Don't know why....

    PetscMPIInt :: io_rank
    PetscMPIInt :: hdf5_read_group_size, hdf5_write_group_size
    PetscBool :: broadcast_read

#if defined(SCORPIO)
    PetscMPIInt :: ioread_group_id, iowrite_group_id
#endif

    character(len=MAXSTRINGLENGTH) :: io_buffer

    PetscInt :: fid_out
    PetscInt :: fid_inputrecord

    ! defines the mode (e.g. mph, richards, vadose, etc.
    character(len=MAXWORDLENGTH) :: flowmode
    PetscInt :: iflowmode
    PetscInt :: iflow_sub_mode
    character(len=MAXWORDLENGTH) :: tranmode
    PetscInt :: itranmode

    PetscInt :: nphase
    PetscInt :: liquid_phase
    PetscInt :: gas_phase
    PetscInt :: oil_phase
    PetscInt, pointer :: phase_map(:)
    PetscInt :: nflowdof
    PetscInt :: nflowspec
    PetscInt :: nmechdof
    PetscInt :: nsec_cells
    PetscInt :: nwells
    PetscInt :: neos_table_indices
    PetscBool :: use_th_freezing

    PetscBool :: surf_flow_on
    PetscInt :: nsurfflowdof
    PetscInt :: subsurf_surf_coupling
    PetscInt :: surface_flow_formulation
    PetscReal :: surf_flow_time, surf_flow_dt
    PetscReal :: surf_subsurf_coupling_time
    PetscReal :: surf_subsurf_coupling_flow_dt
    PetscReal :: surf_restart_time
    PetscBool :: surf_restart_flag
    character(len=MAXSTRINGLENGTH) :: surf_initialize_flow_filename
    character(len=MAXSTRINGLENGTH) :: surf_restart_filename

    PetscBool :: geomech_on
    PetscBool :: geomech_initial
    PetscInt :: ngeomechdof
    PetscInt :: n_stress_strain_dof
    PetscReal :: geomech_time
    PetscInt :: geomech_subsurf_coupling
    PetscReal :: geomech_gravity(3)
    PetscBool :: sec_vars_update
    PetscInt :: air_pressure_id
    PetscInt :: capillary_pressure_id
    PetscInt :: vapor_pressure_id
    PetscInt :: saturation_pressure_id
    PetscInt :: water_id  ! index of water component dof
    PetscInt :: air_id  ! index of air component dof
    PetscInt :: oil_id  ! index of oil component dof
    PetscInt :: energy_id  ! index of energy dof

    PetscInt :: ntrandof

    PetscInt :: iflag
    PetscInt :: status
    PetscBool :: input_record
    !geh: remove once legacy code is gone.
!    PetscBool :: init_stage
    ! these flags are for printing outside of time step loop
    PetscBool :: print_to_screen
    PetscBool :: print_to_file
    ! these flags are for printing within time step loop where printing may
    ! need to be temporarily turned off to accommodate periodic screen outout.
    PetscBool :: print_screen_flag
    PetscBool :: print_file_flag
    PetscInt :: verbosity  ! Values >0 indicate additional console output.

    PetscReal :: uniform_velocity(3)

    ! Program options
    PetscBool :: use_matrix_free  ! If true, do not form the Jacobian.

    PetscBool :: use_isothermal
    PetscBool :: use_mc           ! If true, multiple continuum formulation is used.
    PetscBool :: set_secondary_init_temp  ! If true, then secondary init temp is different from prim. init temp
    PetscBool :: set_secondary_init_conc

    PetscBool :: update_flow_perm ! If true, permeability changes due to pressure

    PetscInt :: ice_model         ! specify water/ice/vapor phase partitioning model
    PetscReal:: frzthw_halfwidth  ! freezing-thawing smoothing half-width (oC)
      
    PetscReal :: flow_time, tran_time, time  ! The time elapsed in the simulation.
    PetscReal :: flow_dt ! The size of the time step.
    PetscReal :: tran_dt
    PetscReal :: dt
    PetscReal :: dt_min
    PetscBool :: match_waypoint
    PetscReal :: refactor_dt

    PetscReal :: gravity(3)

    PetscReal :: scale

    PetscReal :: m_nacl

    PetscInt :: ideriv
    PetscInt :: idt_switch
    PetscReal :: reference_temperature
    PetscReal :: reference_pressure
    PetscReal :: reference_density(2)
    PetscReal :: reference_porosity
    PetscReal :: reference_saturation

    PetscBool :: converged
    PetscInt :: convergence

    PetscReal :: infnorm_res_sec  ! inf. norm of secondary continuum rt residual

    PetscReal :: minimum_hydrostatic_pressure

!   table lookup
    PetscInt :: itable
    PetscInt :: co2eos
    character(len=MAXSTRINGLENGTH) :: co2_database_filename

    PetscBool :: restart_flag
    PetscReal :: restart_time
    character(len=MAXSTRINGLENGTH) :: restart_filename
    character(len=MAXSTRINGLENGTH) :: input_filename

    PetscLogDouble :: start_time
    PetscBool :: wallclock_stop_flag
    PetscLogDouble :: wallclock_stop_time

    PetscInt :: log_stage(10)

    PetscBool :: numerical_derivatives_multi_coupling
    PetscBool :: compute_statistics
    PetscBool :: compute_mass_balance_new
    PetscBool :: mass_bal_detailed
    PetscBool :: use_touch_options
    PetscBool :: overwrite_restart_transport
    PetscBool :: overwrite_restart_flow
    PetscInt :: io_handshake_buffer_size

    character(len=MAXSTRINGLENGTH) :: initialize_flow_filename
    character(len=MAXSTRINGLENGTH) :: initialize_transport_filename

    character(len=MAXSTRINGLENGTH) :: input_prefix
    character(len=MAXSTRINGLENGTH) :: global_prefix
    character(len=MAXWORDLENGTH) :: group_prefix
    !PO
    character(len=MAXSTRINGLENGTH) :: output_file_name_prefix
    character(len=MAXSTRINGLENGTH) :: output_dir
    !PO end

    PetscBool :: steady_state
    PetscBool :: use_matrix_buffer
    PetscBool :: force_newton_iteration
    PetscBool :: use_upwinding
    PetscBool :: out_of_table

    ! Specify secondary continuum solver
    PetscInt :: secondary_continuum_solver     ! Specify secondary continuum solver

    PetscInt :: subsurface_simulation_type

    ! For WIPP_type pc-sat characteristic curves that use Pct
    PetscBool :: pct_updated

    ! Type of averaging scheme for relative permeability
    PetscInt :: rel_perm_aveg
    PetscBool :: first_step_after_restart

    ! value of a cutoff for Manning's/Infiltration velocity
    PetscReal :: max_manning_velocity
    PetscReal :: max_infiltration_velocity

    ! when the scaling factor is too small, stop in reactive transport
    PetscReal :: min_allowable_scale

    PetscBool :: print_ekg

    ! flag to use inline surface flow in Richards mode
    PetscBool :: inline_surface_flow
    PetscReal :: inline_surface_Mannings_coeff
    character(len=MAXSTRINGLENGTH) :: inline_surface_region_name
    
#ifdef CLM_PFLOTRAN
    character(len=MAXSTRINGLENGTH) :: input_dir
    PetscBool :: mapping_files
#endif

    PetscReal :: debug_tol
    PetscBool :: matcompare_reldiff
    PetscBool :: use_GP

  end type option_type

  PetscInt, parameter, public :: SUBSURFACE_SIM_TYPE = 1
  PetscInt, parameter, public :: MULTISIMULATION_SIM_TYPE = 2
  PetscInt, parameter, public :: STOCHASTIC_SIM_TYPE = 3

  interface printMsg
    module procedure printMsg1
    module procedure printMsg2
  end interface

  interface printMsgAnyRank
    module procedure printMsgAnyRank1
    module procedure printMsgAnyRank2
  end interface

  interface printMsgByRank
    module procedure printMsgByRank1
    module procedure printMsgByRank2
  end interface

  interface printErrMsgByRank
    module procedure printErrMsgByRank1
    module procedure printErrMsgByRank2
  end interface

  interface printErrMsgNoStopByRank
    module procedure printErrMsgNoStopByRank1
    module procedure printErrMsgNoStopByRank2
  end interface

  interface printErrMsg
    module procedure printErrMsg1
    module procedure printErrMsg2
  end interface

  interface printWrnMsg
    module procedure printWrnMsg1
    module procedure printWrnMsg2
  end interface

  interface OptionInitMPI
    module procedure OptionInitMPI1
    module procedure OptionInitMPI2
  end interface

  public :: OptionCreate, &
            OptionCheckCommandLine, &
            printErrMsg, &
            printErrMsgByRank, &
            printWrnMsg, &
            printMsg, &
            printMsgAnyRank, &
            printMsgByRank, &
            printMsgByCell, &
            printErrMsgNoStopByRank, &
            printVerboseMsg, &
            OptionCheckTouch, &
            OptionPrintToScreen, &
            OptionPrintToFile, &
            OptionInitRealization, &
            OptionMeanVariance, &
            OptionMaxMinMeanVariance, &
            OptionInitMPI, &
            OptionInitPetsc, &
            OptionDivvyUpSimulations, &
            OptionCreateProcessorGroups, &
            OptionBeginTiming, &
            OptionEndTiming, &
            OptionFinalize, &
            OptionDestroy

contains

! ************************************************************************** !

function OptionCreate()
  !
  ! Allocates and initializes a new Option object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  type(option_type), pointer :: OptionCreate

  type(option_type), pointer :: option

  allocate(option)
  option%flow => OptionFlowCreate()
  option%transport => OptionTransportCreate()

  ! DO NOT initialize members of the option type here.  One must decide
  ! whether the member needs initialization once for all stochastic
  ! simulations or initialization for every realization (e.g. within multiple
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionInitAll(option)
  OptionCreate => option

end function OptionCreate

! ************************************************************************** !

subroutine OptionInitAll(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  type(option_type) :: option

  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)

  call OptionFlowInitAll(option%flow)
  call OptionTransportInitAll(option%transport)

  option%id = 0
  option%successful_exit_code = 0

  option%global_comm = 0
  option%global_rank = 0
  option%global_commsize = 0
  option%global_group = 0

  option%mycomm = 0
  option%myrank = 0
  option%mycommsize = 0
  option%mygroup = 0
  option%mygroup_id = 0

  option%input_prefix = 'pflotran'
  option%group_prefix = ''
  option%global_prefix = ''
  option%output_file_name_prefix = ''
  option%output_dir = ''

  option%broadcast_read = PETSC_FALSE
  option%io_rank = 0
  option%hdf5_read_group_size = 0
  option%hdf5_write_group_size = 0

  option%input_record = PETSC_FALSE
  option%print_screen_flag = PETSC_FALSE
  option%print_file_flag = PETSC_FALSE
  option%print_to_screen = PETSC_TRUE
  option%print_to_file = PETSC_TRUE
  option%verbosity = 0

  option%input_filename = ''

  option%use_upwinding = PETSC_TRUE

  option%out_of_table = PETSC_FALSE

  option%subsurface_simulation_type = SUBSURFACE_SIM_TYPE

  option%rel_perm_aveg = UPWIND
  option%first_step_after_restart = PETSC_FALSE

  call OptionInitRealization(option)

end subroutine OptionInitAll

! ************************************************************************** !

subroutine OptionInitRealization(option)
  !
  ! Initializes option variables specific to a single
  ! realization
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  type(option_type) :: option

  ! These variables should be initialized once at the beginning of every
  ! PFLOTRAN realization or simulation of a single realization
  call OptionFlowInitRealization(option%flow)
  call OptionTransportInitRealization(option%transport)


  option%fid_out = OUT_UNIT
  option%fid_inputrecord = INPUT_RECORD_UNIT

  option%iflag = 0
  option%io_buffer = ''

  option%use_isothermal = PETSC_FALSE
  option%use_matrix_free = PETSC_FALSE
  option%use_mc = PETSC_FALSE
  option%set_secondary_init_temp = PETSC_FALSE
  option%ice_model = PAINTER_EXPLICIT
  option%frzthw_halfwidth = UNINITIALIZED_DOUBLE
  option%set_secondary_init_conc = PETSC_FALSE

  option%update_flow_perm = PETSC_FALSE

  option%flowmode = ""
  option%iflowmode = NULL_MODE
  option%iflow_sub_mode = NULL_MODE
  option%nflowdof = 0
  option%nmechdof = 0
  option%nsec_cells = 0
  option%nwells = 0
  option%neos_table_indices = 0
  option%use_th_freezing = PETSC_FALSE

  option%nsurfflowdof = 0
  option%surf_flow_on = PETSC_FALSE
  option%subsurf_surf_coupling = DECOUPLED
  option%surface_flow_formulation = DIFFUSION_WAVE
  option%surf_flow_dt = 0.d0
  option%surf_flow_time =0.d0
  option%surf_subsurf_coupling_time = 0.d0
  option%surf_subsurf_coupling_flow_dt = 0.d0
  option%surf_initialize_flow_filename = ""
  option%surf_restart_filename = ""
  option%surf_restart_flag = PETSC_FALSE
  option%surf_restart_time = UNINITIALIZED_DOUBLE

  option%geomech_on = PETSC_FALSE
  option%geomech_initial = PETSC_FALSE
  option%ngeomechdof = 0
  option%n_stress_strain_dof = 0
  option%geomech_time = 0.d0
  option%geomech_subsurf_coupling = 0
  option%geomech_gravity(:) = 0.d0
  option%geomech_gravity(3) = -1.d0*EARTH_GRAVITY    ! m/s^2

  option%tranmode = ""
  option%itranmode = NULL_MODE
  option%ntrandof = 0

  nullify(option%phase_map)
  option%nphase = 0
  option%liquid_phase = 0
  option%oil_phase = 0
  option%gas_phase = 0

  option%air_pressure_id = 0
  option%capillary_pressure_id = 0
  option%vapor_pressure_id = 0
  option%saturation_pressure_id = 0

  option%water_id = 0
  option%air_id = 0
  option%energy_id = 0

  option%uniform_velocity = 0.d0

!-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
  !TODO(geh): move to option%flow.F90
  option%reference_pressure = 101325.d0
  option%reference_temperature = 25.d0
  option%reference_density = 0.d0
  option%reference_porosity = 0.25d0
  option%reference_saturation = 1.d0

  option%converged = PETSC_FALSE
  option%convergence = CONVERGENCE_OFF

  option%infnorm_res_sec = 0.d0

  option%minimum_hydrostatic_pressure = -1.d20

  !set scale factor for heat equation, i.e. use units of MJ for energy
  option%scale = 1.d-6

  option%ideriv = 1

  option%gravity(:) = 0.d0
  option%gravity(3) = -1.d0*EARTH_GRAVITY ! m/s^2

  !physical constants and defult variables
!  option%difaq = 1.d-9 ! m^2/s read from input file
!  option%difaq = 0.d0
!  option%delhaq = 12.6d0 ! kJ/mol read from input file
!  option%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  ! default brine concentrations
  option%m_nacl = 0.d0

!  option%disp = 0.d0

  option%restart_flag = PETSC_FALSE
  option%restart_filename = ""
  option%restart_time = UNINITIALIZED_DOUBLE

  option%start_time = 0.d0
  option%wallclock_stop_flag = PETSC_FALSE
  option%wallclock_stop_time = 0.d0

  option%log_stage = 0

  option%numerical_derivatives_multi_coupling = PETSC_FALSE
  option%compute_statistics = PETSC_FALSE
  option%compute_mass_balance_new = PETSC_FALSE

!fmy: mass_balance for bc/ss IS needed by default if coupled with CLM
#ifdef CLM_PFLOTRAN
  option%compute_mass_balance_new = PETSC_TRUE
  option%input_dir = ""
  option%mapping_files = PETSC_FALSE
  ! user-defined CLM-PFLOTRAN mesh maps NOT provided (default)
#endif
!fmy: mass_balance for bc/ss IS needed by default if coupled with CLM

  option%mass_bal_detailed = PETSC_FALSE

  option%use_touch_options = PETSC_FALSE
  option%overwrite_restart_transport = PETSC_FALSE
  option%overwrite_restart_flow = PETSC_FALSE

  option%time = 0.d0
  option%flow_dt = 0.d0
  option%tran_dt = 0.d0
  option%dt = 0.d0
  option%dt_min = 1.d-20   ! Ten zeptoseconds
  option%refactor_dt = 0.d0
  option%match_waypoint = PETSC_FALSE

  option%io_handshake_buffer_size = 0

  option%initialize_flow_filename = ''
  option%initialize_transport_filename = ''

  option%steady_state = PETSC_FALSE

  option%itable = 0

  option%idt_switch = -1

  option%use_matrix_buffer = PETSC_FALSE
  option%status = PROCEED
  option%force_newton_iteration = PETSC_FALSE
  option%secondary_continuum_solver = 1

  ! initially set to a large value to effectively disable
  option%max_manning_velocity = 1.d20
  option%max_infiltration_velocity = 1.d20

  ! when the scaling factor is too small, stop in reactive transport
  option%min_allowable_scale = 1.0d-10

  option%print_ekg = PETSC_FALSE

  option%pct_updated = PETSC_FALSE

  option%inline_surface_flow           = PETSC_FALSE
  option%inline_surface_Mannings_coeff = 0.02d0
  option%inline_surface_region_name    = ""

  option%debug_tol = 1.d0
  option%matcompare_reldiff = PETSC_FALSE
  option%use_GP = PETSC_FALSE

end subroutine OptionInitRealization

! ************************************************************************** !

subroutine OptionCheckCommandLine(option)
  !
  ! Checks all PETSc options on input
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option

  PetscBool :: option_found
  PetscInt :: temp_int
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string

  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-buffer_matrix", &
                           option%use_matrix_buffer, ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-snes_mf", &
                           option%use_matrix_free, ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-use_isothermal", &
                           option%use_isothermal, ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-use_mc", &
                           option%use_mc, ierr);CHKERRQ(ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                             '-restart', option%restart_filename, &
                             option%restart_flag, ierr);CHKERRQ(ierr)
  ! check on possible modes
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-use_richards", &
                           option_found, ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "richards"
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-use_thc", &
                           option_found, ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "thc"
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-use_mph", &
                           option_found, ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "mph"
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                           PETSC_NULL_CHARACTER, "-use_flash2", &
                           option_found, ierr);CHKERRQ(ierr)
  if (option_found) option%flowmode = "flash2"

end subroutine OptionCheckCommandLine

! ************************************************************************** !

subroutine printErrMsg1(option)
  !
  ! Prints the error message from p0 and stops
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option

  call printErrMsg2(option,option%io_buffer)

end subroutine printErrMsg1

! ************************************************************************** !

subroutine printErrMsg2(option,string)
  !
  ! Prints the error message from p0 and stops
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  PetscBool :: petsc_initialized
  PetscErrorCode :: ierr

  if (OptionPrintToScreen(option)) then
    print *
    print *, 'ERROR: ' // trim(string)
    print *
    print *, 'Stopping!'
  endif
  call MPI_Barrier(option%mycomm,ierr)
  call PetscInitialized(petsc_initialized, ierr);CHKERRQ(ierr)
  if (petsc_initialized) then
    call PetscFinalize(ierr);CHKERRQ(ierr)
  endif
  stop

end subroutine printErrMsg2

! ************************************************************************** !

subroutine printErrMsgByRank1(option)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option

  call printErrMsgByRank2(option,option%io_buffer)

end subroutine printErrMsgByRank1

! ************************************************************************** !

subroutine printErrMsgByRank2(option,string)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  character(len=MAXWORDLENGTH) :: word

  if (option%print_to_screen) then
    write(word,*) option%myrank
    print *
    print *, 'ERROR(' // trim(adjustl(word)) // '): ' // trim(string)
    print *
    print *, 'Stopping!'
  endif
  stop

end subroutine printErrMsgByRank2

! ************************************************************************** !

! ************************************************************************** !

subroutine printErrMsgNoStopByRank1(option)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option

  call printErrMsgNoStopByRank2(option,option%io_buffer)

end subroutine printErrMsgNoStopByRank1

! ************************************************************************** !

subroutine printErrMsgNoStopByRank2(option,string)
  !
  ! Prints the error message from processor with error along
  ! with rank
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/11
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  character(len=MAXWORDLENGTH) :: word

  if (option%print_to_screen) then
    write(word,*) option%myrank
    print *
    print *, 'ERROR(' // trim(adjustl(word)) // '): ' // trim(string)
    print *
  endif

end subroutine printErrMsgNoStopByRank2

! ************************************************************************** !

subroutine printWrnMsg1(option)
  !
  ! Prints the warning message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option

  call printWrnMsg2(option,option%io_buffer)

end subroutine printWrnMsg1

! ************************************************************************** !

subroutine printWrnMsg2(option,string)
  !
  ! Prints the warning message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  if (OptionPrintToScreen(option)) print *, 'WARNING: ' // trim(string)

end subroutine printWrnMsg2

! ************************************************************************** !

subroutine printMsg1(option)
  !
  ! Prints the message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !

  implicit none

  type(option_type) :: option

  call printMsg2(option,option%io_buffer)

end subroutine printMsg1

! ************************************************************************** !

subroutine printMsg2(option,string)
  !
  ! Prints the message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  if (OptionPrintToScreen(option)) print *, trim(string)

end subroutine printMsg2

! ************************************************************************** !

subroutine printMsgAnyRank1(option)
  !
  ! Prints the message from any processor core
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  !

  implicit none

  type(option_type) :: option

  if (option%print_to_screen) call printMsgAnyRank2(option%io_buffer)

end subroutine printMsgAnyRank1

! ************************************************************************** !

subroutine printMsgAnyRank2(string)
  !
  ! Prints the message from any processor core
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  !

  implicit none

  character(len=*) :: string
  
  print *, trim(string)

end subroutine printMsgAnyRank2

! ************************************************************************** !

subroutine printMsgByRank1(option)
  !
  ! Prints a message from processor along with rank
  !
  ! Author: Glenn Hammond
  ! Date: 03/27/12
  !

  implicit none

  type(option_type) :: option

  call printMsgByRank2(option,option%io_buffer)

end subroutine printMsgByRank1

! ************************************************************************** !

subroutine printMsgByRank2(option,string)
  !
  ! Prints a message from processor along with rank
  !
  ! Author: Glenn Hammond
  ! Date: 03/27/12
  !

  implicit none

  type(option_type) :: option
  character(len=*) :: string

  character(len=MAXWORDLENGTH) :: word

  if (option%print_to_screen) then
    write(word,*) option%myrank
    print *, '(' // trim(adjustl(word)) // '): ' // trim(string)
  endif

end subroutine printMsgByRank2

! ************************************************************************** !

subroutine printMsgByCell(option,cell_id,string)
  !
  ! Prints the message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !

  implicit none

  type(option_type) :: option
  PetscInt :: cell_id
  character(len=*) :: string

  character(len=MAXWORDLENGTH) :: word

  write(word,*) cell_id
  word = adjustl(word)
  option%io_buffer = trim(string) // ' for cell ' // trim(word) // '.'
  call printMsgByRank(option)

end subroutine printMsgByCell

! ************************************************************************** !

subroutine printVerboseMsg(option)
  !
  ! Prints the message from p0
  !
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  !

  implicit none

  type(option_type) :: option

  if (option%verbosity > 0) then
    call printMsg(option,option%io_buffer)
  endif

end subroutine printVerboseMsg

! ************************************************************************** !

function OptionCheckTouch(option,filename)
  !
  ! Users can steer the code by touching files.
  !
  ! Author: Glenn Hammond
  ! Date: 03/04/08
  !

  implicit none

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: filename

  PetscInt :: ios
  PetscInt :: fid = 86
  PetscBool :: OptionCheckTouch
  PetscErrorCode :: ierr

  OptionCheckTouch = PETSC_FALSE

  if (option%myrank == option%io_rank) &
    open(unit=fid,file=trim(filename),status='old',iostat=ios)
  call MPI_Bcast(ios,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                 option%mycomm,ierr)

  if (ios == 0) then
    if (option%myrank == option%io_rank) close(fid,status='delete')
    OptionCheckTouch = PETSC_TRUE
  endif

end function OptionCheckTouch

! ************************************************************************** !

function OptionPrintToScreen(option)
  !
  ! Determines whether printing should occur
  !
  ! Author: Glenn Hammond
  ! Date: 12/09/08
  !

  implicit none

  type(option_type) :: option

  PetscBool :: OptionPrintToScreen

  if (option%myrank == option%io_rank .and. option%print_to_screen) then
    OptionPrintToScreen = PETSC_TRUE
  else
    OptionPrintToScreen = PETSC_FALSE
  endif

end function OptionPrintToScreen

! ************************************************************************** !

function OptionPrintToFile(option)
  !
  ! Determines whether printing to file should occur
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/09
  !

  implicit none

  type(option_type) :: option

  PetscBool :: OptionPrintToFile

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    OptionPrintToFile = PETSC_TRUE
  else
    OptionPrintToFile = PETSC_FALSE
  endif

end function OptionPrintToFile

! ************************************************************************** !

subroutine OptionMaxMinMeanVariance(value,max,min,mean,variance, &
                                    calculate_variance,option)
  !
  ! Calculates the maximum, minumum, mean and
  ! optionally variance of a number across processor
  ! cores
  !
  ! Author: Glenn Hammond
  ! Date: 06/01/10
  !

  implicit none

  type(option_type) :: option
  PetscReal :: value
  PetscReal :: max
  PetscReal :: min
  PetscReal :: mean
  PetscReal :: variance
  PetscBool :: calculate_variance

  PetscReal :: temp_real_in(2), temp_real_out(2)
  PetscErrorCode :: ierr

  temp_real_in(1) = value
  temp_real_in(2) = -1.d0*value
  call MPI_Allreduce(temp_real_in,temp_real_out,TWO_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MAX,option%mycomm,ierr)
  max = temp_real_out(1)
  min = -1.d0*temp_real_out(2)

  call OptionMeanVariance(value,mean,variance,calculate_variance,option)

end subroutine OptionMaxMinMeanVariance

! ************************************************************************** !

subroutine OptionMeanVariance(value,mean,variance,calculate_variance,option)
  !
  ! Calculates the mean and optionally variance of a number
  ! across processor cores
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/10
  !

  implicit none

  type(option_type) :: option
  PetscReal :: value
  PetscReal :: mean
  PetscReal :: variance
  PetscBool :: calculate_variance

  PetscReal :: temp_real
  PetscErrorCode :: ierr

  call MPI_Allreduce(value,temp_real,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr)
  mean = temp_real / dble(option%mycommsize)

  if (calculate_variance) then
    temp_real = value-mean
    temp_real = temp_real*temp_real
    call MPI_Allreduce(temp_real,variance,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_SUM,option%mycomm,ierr)
    variance = variance / dble(option%mycommsize)
  endif

end subroutine OptionMeanVariance

! ************************************************************************** !

subroutine OptionInitMPI1(option)
  !
  ! Initializes base MPI communicator
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  !

  implicit none

  type(option_type) :: option

  PetscErrorCode :: ierr

  ! when coupling with CLM, MPI has already initialized, i.e. using MPI_COMM_WORLD directly.
  ! (actually this subroutine is not called in pflotran_clm_main.F90)
#ifndef CLM_PFLOTRAN
  call MPI_Init(ierr)
#endif

  call OptionInitMPI2(option,MPI_COMM_WORLD)

end subroutine OptionInitMPI1

! ************************************************************************** !

subroutine OptionInitMPI2(option,communicator)
  !
  ! Initializes base MPI communicator
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  !

  implicit none

  type(option_type) :: option

  PetscMPIInt :: communicator
  PetscErrorCode :: ierr

  option%global_comm = communicator
  call MPI_Comm_rank(communicator,option%global_rank, ierr)
  call MPI_Comm_size(communicator,option%global_commsize,ierr)
  call MPI_Comm_group(communicator,option%global_group,ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group

end subroutine OptionInitMPI2

! ************************************************************************** !

subroutine OptionInitPetsc(option)
  !
  ! Initialization of PETSc.
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  !

  use Logging_module

  implicit none

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  PETSC_COMM_WORLD = option%mycomm
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)    !fmy: tiny memory leak here (don't know why)

  if (option%verbosity > 0) then
    call PetscLogDefaultBegin(ierr);CHKERRQ(ierr)
    string = '-log_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                  string, ierr);CHKERRQ(ierr)
  endif

  call LoggingCreate()

end subroutine OptionInitPetsc

! ************************************************************************** !

subroutine OptionBeginTiming(option)
  !
  ! Start outer timing.
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  !

  use Logging_module

  implicit none

#include "petsc/finclude/petsclog.h"

  type(option_type) :: option

  PetscLogDouble :: timex_wall
  PetscErrorCode :: ierr

  call PetscTime(timex_wall, ierr);CHKERRQ(ierr)
  option%start_time = timex_wall

end subroutine OptionBeginTiming

! ************************************************************************** !

subroutine OptionEndTiming(option)
  !
  ! End timing.
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  !

  use Logging_module

  implicit none

#include "petsc/finclude/petsclog.h"

  type(option_type) :: option

  PetscLogDouble :: timex_wall
  PetscErrorCode :: ierr

  ! Final Time
  call PetscTime(timex_wall, ierr);CHKERRQ(ierr)

  if (option%myrank == option%io_rank) then

    if (option%print_to_screen) then
      write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall-option%start_time, &
        (timex_wall-option%start_time)/60.d0, &
        (timex_wall-option%start_time)/3600.d0
    endif
    if (option%print_to_file) then
      write(option%fid_out,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall-option%start_time, &
        (timex_wall-option%start_time)/60.d0, &
        (timex_wall-option%start_time)/3600.d0
    endif
  endif

end subroutine OptionEndTiming

! ************************************************************************** !

subroutine OptionDivvyUpSimulations(option,filenames)
  !
  ! Divides simulation in to multple simulations with
  ! multiple input decks
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  !

  implicit none

  type(option_type) :: option

  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  i = size(filenames)
  call OptionCreateProcessorGroups(option,i)
  option%input_filename = filenames(option%mygroup_id)
  i = index(option%input_filename,'.',PETSC_TRUE)
  if (i > 1) then
    i = i-1
  else
    ! for some reason len_trim doesn't work on MS Visual Studio in
    ! this location
    i = len(trim(option%input_filename))
  endif
  option%global_prefix = option%input_filename(1:i)
  write(string,*) option%mygroup_id
  option%group_prefix = 'G' // trim(adjustl(string))

end subroutine OptionDivvyUpSimulations

! ************************************************************************** !

subroutine OptionCreateProcessorGroups(option,num_groups)
  !
  ! Splits MPI_COMM_WORLD into N separate
  ! processor groups
  !
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  !

  implicit none

  type(option_type) :: option
  PetscInt :: num_groups

  PetscInt :: local_commsize
  PetscInt :: offset, delta, remainder
  PetscInt :: igroup
  PetscMPIInt :: mycolor_mpi, mykey_mpi
  PetscErrorCode :: ierr

  local_commsize = option%global_commsize / num_groups
  remainder = option%global_commsize - num_groups * local_commsize
  offset = 0
  do igroup = 1, num_groups
    delta = local_commsize
    if (igroup < remainder) delta = delta + 1
    if (option%global_rank >= offset .and. &
        option%global_rank < offset + delta) exit
    offset = offset + delta
  enddo
  mycolor_mpi = igroup
  option%mygroup_id = igroup
  mykey_mpi = option%global_rank - offset
  call MPI_Comm_split(MPI_COMM_WORLD,mycolor_mpi,mykey_mpi,option%mycomm,ierr)
  call MPI_Comm_group(option%mycomm,option%mygroup,ierr)

  call MPI_Comm_rank(option%mycomm,option%myrank, ierr)
  call MPI_Comm_size(option%mycomm,option%mycommsize,ierr)

end subroutine OptionCreateProcessorGroups

! ************************************************************************** !

subroutine OptionFinalize(option)
  !
  ! End the simulation.
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  !

  use Logging_module

  implicit none

  type(option_type), pointer :: option

  PetscInt :: iflag
  PetscErrorCode :: ierr

  call LoggingDestroy()
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-options_left','no',ierr);CHKERRQ(ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-objects_left','yes',ierr);CHKERRQ(ierr)
  call MPI_Barrier(option%global_comm,ierr)
  iflag = option%successful_exit_code
  call OptionDestroy(option)

  ! NOTE: MPI communitcator from CLM should have already initialized beyond PFLOTRAN.
  !       SO it must also be finalized beyond PFLOTRAN. Otherwise it causes MPI_Finalize failure.
#ifndef CLM_PFLOTRAN
  call PetscFinalize(ierr);CHKERRQ(ierr)
  call MPI_Finalize(ierr)
  call exit(iflag)
#endif
  
end subroutine OptionFinalize

! ************************************************************************** !

subroutine OptionDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  !

  implicit none

  type(option_type), pointer :: option

  call OptionFlowDestroy(option%flow)
  call OptionTransportDestroy(option%transport)
  ! all kinds of stuff needs to be added here.
  if (associated(option%phase_map) ) then
    deallocate(option%phase_map)
    nullify(option%phase_map)
  end if
  ! all the below should be placed somewhere other than option.F90

  deallocate(option)
  nullify(option)

end subroutine OptionDestroy

end module Option_module
