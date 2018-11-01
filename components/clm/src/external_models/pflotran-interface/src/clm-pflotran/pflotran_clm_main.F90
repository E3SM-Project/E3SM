#ifdef CLM_PFLOTRAN

module pflotran_clm_main_module

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
  use petscsys
  use petscvec

  use PFLOTRAN_Constants_module
  use Option_module, only : option_type
  use Simulation_Base_class, only : simulation_base_type
  use Multi_Simulation_module, only : multi_simulation_type
  use Realization_Base_class, only : realization_base_type
  use Utility_module, only : where_checkerr

  use Mapping_module
  use clm_pflotran_interface_data


  implicit none

  private

  type, public :: pflotran_model_type
    class(simulation_base_type),  pointer :: simulation
    type(multi_simulation_type), pointer :: multisimulation
    type(option_type),      pointer :: option
    PetscReal :: pause_time_1
    PetscReal :: pause_time_2
    type(mapping_type),                pointer :: map_clm_sub_to_pf_sub
    type(mapping_type),                pointer :: map_clm_2dtop_to_pf_2dtop
    type(mapping_type),                pointer :: map_pf_sub_to_clm_sub
    type(mapping_type),                pointer :: map_pf_2dtop_to_clm_2dtop

    type(mapping_type),                pointer :: map_clm_2dbot_to_pf_2dbot
    type(mapping_type),                pointer :: map_pf_2dbot_to_clm_2dbot
     
    PetscInt :: nlclm
    PetscInt :: ngclm

  end type pflotran_model_type
  !
  public::pflotranModelCreate,               &
       ! PF running
       pflotranModelStepperRunInit,          &
       pflotranModelUpdateFinalWaypoint,     &
       pflotranModelStepperRunTillPauseTime, &
       pflotranModelSetupRestart,            &
       pflotranModelStepperRunFinalize,      &
       pflotranModelStepperCheckpoint,       &
       pflotranModelDestroy,                 &
       ! soil domain
       pflotranModelSetSoilDimension,         &
       ! Soil physical properties
       pflotranModelSetSoilProp,              &
       pflotranModelResetSoilPorosityFromCLM, &
       pflotranModelGetSoilPropFromPF,        &
       pflotranModelUpdateTHfromCLM,          &      ! from CLM to PF's global vars to drive BGC
       pflotranModelGetTemperatureFromPF,     &
       pflotranModelGetSaturationFromPF,      &
       ! BGC subroutines
       pflotranModelGetRTspecies,               &
       pflotranModelSetSOMKfromCLM,             &
       pflotranModelSetBGCRatesFromCLM,         &
       pflotranModelUpdateAqConcFromCLM,        &
       pflotranModelUpdateAqGasesFromCLM,       &
       pflotranModelSetBgcConcFromCLM,          &
       pflotranModelGetBgcVariablesFromPF,      &
       ! TH subroutines
       pflotranModelSetInternalTHStatesfromCLM, &    ! T/H states from CLM to PFLOTRAN flow mode's field%**
       pflotranModelUpdateHSourceSink,          &    ! water src/sink (e.g., ET)
       pflotranModelUpdateSubsurfTCond,         &    ! thermal BC
       pflotranModelSetSoilHbcsFromCLM,         &    ! water BC
       !
       pflotranModelGetBCMassBalanceDeltaFromPF      ! Mass-Balance at BCs

  private :: &
       pflotranModelInsertWaypoint,          &
       pflotranModelDeleteWaypoint

!------------------------------------------------------------

  PetscReal, parameter :: xeps0_c = 1.0d-50
  PetscReal, parameter :: xeps0_n = 1.0d-51

  character(len=MAXWORDLENGTH) :: subname
!------------------------------------------------------------

contains

! ************************************************************************************ !

  subroutine pflotranModelCreate(mpicomm, pflotran_inputdir, pflotran_prefix, model)
  ! 
  ! Allocates and initializes the pflotranModel object.
  ! It performs the same sequence of commands as done in pflotran.F90
  ! before model integration is performed by the call to StepperRun()
  ! routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    use Option_module
    use Simulation_Base_class
    use Multi_Simulation_module
    use Factory_PFLOTRAN_module
    use Factory_Subsurface_module
    use PFLOTRAN_Constants_module
    use Output_Aux_module, only : INSTANTANEOUS_VARS
    use PFLOTRAN_Provenance_module, only : PrintProvenanceToScreen
  
    implicit none

    PetscInt, intent(in) :: mpicomm
    character(len=256), intent(in) :: pflotran_inputdir
    character(len=256), intent(in) :: pflotran_prefix

    type(pflotran_model_type),      pointer :: model

    allocate(model)

    nullify(model%simulation)
    nullify(model%multisimulation)
    nullify(model%option)

    model%option => OptionCreate()
    call OptionInitMPI(model%option, mpicomm)
    call PFLOTRANInitializePrePetsc(model%multisimulation, model%option)

    ! NOTE(bja) 2013-06-25 : external driver must provide an input
    ! prefix string. If the driver wants to use pflotran.in, then it
    ! should explicitly request that with 'pflotran'.
    if (len(trim(pflotran_prefix)) > 1) then
      model%option%input_prefix = trim(pflotran_prefix)

      if (len(trim(pflotran_inputdir)) > 1) then
        model%option%input_dir = trim(pflotran_inputdir)
        model%option%input_prefix = trim(pflotran_inputdir) // '/' // trim(pflotran_prefix)

      else
        model%option%input_dir = "."
        model%option%input_prefix = trim(pflotran_prefix)
      endif

      model%option%input_filename = trim(model%option%input_prefix) // '.in'
      model%option%global_prefix = trim(pflotran_prefix)
    else
      model%option%io_buffer = 'The external driver must provide the ' // &
           'pflotran input file prefix.'
      call printErrMsg(model%option)
    end if

    call OptionInitPetsc(model%option)
    if (model%option%myrank == model%option%io_rank .and. &
        model%option%print_to_screen) then
      call PrintProvenanceToScreen()
    end if

    ! NOTE(bja, 2013-07-19) GB's Hack to get communicator correctly
    ! setup on mpich/mac. should be generally ok, but may need an
    ! apple/mpich ifdef if it cause problems elsewhere.
    PETSC_COMM_SELF = MPI_COMM_SELF
    PETSC_COMM_WORLD = MPI_COMM_WORLD

    call PFLOTRANInitializePostPetsc(model%simulation, model%multisimulation, model%option)

    model%pause_time_1 = -1.0d0
    model%pause_time_2 = -1.0d0

  end subroutine pflotranModelCreate

! ************************************************************************** !

  subroutine pflotranModelStepperRunInit(model)
  ! 
  ! It performs the same execution of commands
  ! that are carried out in StepperRun() before the model integration
  ! begins over the entire simulation time interval
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    type(pflotran_model_type), pointer, intent(inout) :: model

    call model%simulation%InitializeRun()

  end subroutine pflotranModelStepperRunInit

! ************************************************************************** !
  !
  ! pflotranModelUpdateFinalWaypoint:
  !  Get CLM final timestep and converts it to PFLOTRAN final way point.
  !  And also set an option for turning on/off PF printing

  subroutine pflotranModelUpdateFinalWaypoint(model, waypoint_time, waypoint_dtmax, isprintout)

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Subsurface_class, only : realization_subsurface_type, &
      RealizationAddWaypointsToList

    use Waypoint_module, only : waypoint_type, waypoint_list_type, &
      WaypointCreate, &
      WaypointDeleteFromList, WaypointInsertInList, &
      WaypointListFillIn, WaypointListRemoveExtraWaypnts
    use Units_module, only : UnitsConvertToInternal
    use Option_module, only : printErrMsg

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in)              :: waypoint_time, waypoint_dtmax
    PetscBool, intent(in)              :: isprintout

    type(waypoint_list_type), pointer  :: waypoint_list
    type(waypoint_type), pointer       :: waypoint, waypoint1
    character(len=MAXWORDLENGTH)       :: word, internal_units

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

!-------------------------------------------------------------------------
    option => model%option
    select type (modsim => model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modsim
        realization => simulation%realization
        waypoint_list => simulation%waypoint_list_subsurface

      class default
        option%io_buffer = "pflotranModelUPdateFinalWaypoint is " // &
              "Not support in this mode."
        call printErrMsg(option)
    end select

    ! new final waypoint
    word = 's'
    internal_units = 'sec'
    waypoint1 => WaypointCreate()
    waypoint1%time          = waypoint_time * UnitsConvertToInternal(word, internal_units, option)
    waypoint1%print_snap_output = PETSC_TRUE
    waypoint1%print_checkpoint  = PETSC_FALSE
    waypoint1%final             = PETSC_TRUE
    waypoint1%dt_max        = waypoint_dtmax * UnitsConvertToInternal(word, internal_units, option)

    ! update subsurface-realization final waypoint
    if (associated(realization) .and. associated(simulation)) then
      ! remove original final waypoint
      waypoint => waypoint_list%first
      do
        if (.not.associated(waypoint)) exit
        if (waypoint%final) then
          waypoint%final = PETSC_FALSE
          exit

        else
           waypoint => waypoint%next
        endif
      enddo

      ! insert new final waypoint
      call WaypointInsertInList(waypoint1, waypoint_list)
      call RealizationAddWaypointsToList(realization, waypoint_list)
      call WaypointListFillIn(waypoint_list, option)
      call WaypointListRemoveExtraWaypnts(waypoint_list, option)

      ! turn off the 'print out' if required from CLM
      if(.not.isprintout) then
        if (model%option%io_rank == model%option%myrank) then
          write(model%option%fid_out, *) 'NOTE: h5 output at input-defined interval ' // &
            'for subsurface flow from PFLOTRAN IS OFF! '
        endif
        waypoint => waypoint_list%first
        do
          if (.not.associated(waypoint)) exit
          waypoint%print_snap_output = PETSC_FALSE
          waypoint => waypoint%next
        enddo
      endif

    endif

  end subroutine pflotranModelUpdateFinalWaypoint

  !-------------------------------------------------------------------!

  subroutine pflotranModelInsertWaypoint(model, waypoint_time, waypoint_dtmax, waypoint_final, isprintout)
  !
  ! Inserts a waypoint within the waypoint list
  ! so that the model integration can be paused when that waypoint is
  ! reached
  ! NOTE: It is assumed the 'waypoint_time' is in seconds
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! Revised by Fengming YUAN

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Base_class, only : realization_base_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointInsertInList
    use Units_module, only : UnitsConvertToInternal
    use Option_module, only : printErrMsg

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in)              :: waypoint_time
    PetscReal, intent(in)              :: waypoint_dtmax
    PetscBool, intent(in)              :: waypoint_final
    PetscBool, intent(in)              :: isprintout

    type(waypoint_type), pointer       :: waypoint
    character(len=MAXWORDLENGTH) :: word, internal_units

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    !------------------------
    option => model%option

    word = 's'
    internal_units = 'sec'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word, internal_units, option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%print_snap_output = isprintout
    waypoint%print_obs_output  = PETSC_FALSE
    waypoint%print_checkpoint  = PETSC_FALSE
    waypoint%print_msbl_output = PETSC_FALSE
    waypoint%final             = waypoint_final
    waypoint%dt_max            = waypoint_dtmax * UnitsConvertToInternal(word, internal_units, option)

    select type (modsim => model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modsim
        realization => simulation%realization

        call WaypointInsertInList(waypoint, simulation%waypoint_list_subsurface)

      class default
         model%option%io_buffer = "pflotranModelInsertWaypoint only " // &
              "works on subsurface simulations."
         call printErrMsg(model%option)

    end select

  end subroutine pflotranModelInsertWaypoint

  !-------------------------------------------------------------------!

  subroutine pflotranModelDeleteWaypoint(model, waypoint_time)

    use Option_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Base_class, only : realization_base_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointDeleteFromList
    use Units_module, only : UnitsConvertToInternal

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    type(waypoint_type), pointer       :: waypoint
    PetscReal                          :: waypoint_time
    character(len=MAXWORDLENGTH) :: word, internal_units

    ! ------------------------------------------ !
    option => model%option
    select type (modsim => model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modsim
        realization => simulation%realization

      class default
         option%io_buffer = "pflotranModelDeleteWaypoint only " // &
              "works on subsurface simulations."
         call printErrMsg(option)
    end select
    !
    word = 's'
    internal_units = 'sec'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time * UnitsConvertToInternal(word, internal_units, option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%dt_max            = 86400.0d0

    if (associated(realization) .and. associated(simulation)) then
       call WaypointDeleteFromList(waypoint, simulation%waypoint_list_subsurface)
    end if

    ! when call 'WaypointCreate()', 'waypoint' will be allocated memory,
    ! and the 'WaypointDeleteFromList''s destroy appears not work( TODO checking)
    ! which causes memory leak continuously - it's very dangerous to system if runs long
    if(associated(waypoint)) deallocate(waypoint)
  end subroutine pflotranModelDeleteWaypoint

! ************************************************************************** !

  subroutine pflotranModelSetupRestart(model, restart_stamp)
  !
  ! pflotranModelSetupRestart()
  ! This checks to see if a restart file stamp was provided by the
  ! driver. If so, we set the restart flag and reconstruct the
  ! restart file name. The actual restart is handled by the standard
  ! pflotran mechanism in TimeStepperInitializeRun()
  ! NOTE: this must be called between pflotranModelCreate() and
  ! pflotranModelStepperRunInit()
  !

    use Option_module
    use String_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    character(len=MAXWORDLENGTH) :: restart_stamp

    model%option%io_buffer = 'restart is not implemented in clm-pflotran.' // &
       'AND, pflotran will be initialized from CLM'

    if (.not. StringNull(restart_stamp)) then
       model%option%restart_flag = PETSC_TRUE
       model%option%restart_filename = &
            trim(model%option%global_prefix) // &
            trim(model%option%group_prefix) // &
            '-' // trim(restart_stamp) // '.chk'

       model%option%io_buffer = 'restart file is: ' // &
            trim(model%option%restart_filename)

    end if

    call printWrnMsg(model%option)

  end subroutine pflotranModelSetupRestart

! ************************************************************************** !

  subroutine pflotranModelStepperRunTillPauseTime(model, pause_time, dtime, isprintout)
  !
  ! It performs the model integration
  ! till the specified pause_time.
  ! NOTE: It is assumed 'pause_time' is in seconds.
  ! NOTE(bja, 2013-07) the strange waypoint insertion of t+30min /
  ! deletion of t is to ensure that we always have a valid waypoint in
  ! front of us, but pflotran does not delete them, so we don't want
  ! to accumulate too large of a list.
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! Revised by Fengming YUAN

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in) :: pause_time
    PetscReal, intent(in) :: dtime
    PetscBool, intent(in) :: isprintout

    PetscReal :: pause_time1

    if(isprintout) then
      if (model%option%io_rank == model%option%myrank) then
        write(model%option%fid_out, *) '>>>> Inserting waypoint at pause_time (s) = ', pause_time
        write(model%option%fid_out, *) '>>>> for CLM timestep: ', pause_time/dtime
      endif
    endif

    pause_time1 = pause_time + dtime
    call pflotranModelInsertWaypoint(model, pause_time1, dtime, PETSC_FALSE, isprintout)

    call model%simulation%RunToTime(pause_time)

    call pflotranModelDeleteWaypoint(model, pause_time)

  end subroutine pflotranModelStepperRunTillPauseTime

! ************************************************************************** !

  subroutine pflotranModelStepperCheckpoint(model, id_stamp)
  !
  ! wrapper around StepperCheckpoint
  ! NOTE(bja, 2013-06-27) : the date stamp from clm is 32 characters
  !

    use Option_module

    implicit none

    type(pflotran_model_type), pointer :: model
    character(len=MAXSTRINGLENGTH), intent(in) :: id_stamp

    if (associated(model%simulation%process_model_coupler_list%checkpoint_option)) then
      call model%simulation%process_model_coupler_list%Checkpoint(id_stamp)
    endif

  end subroutine pflotranModelStepperCheckpoint

! ************************************************************************** !

  subroutine pflotranModelStepperRunFinalize(model)
  !
  ! It performs the same execution of commands
  ! that are carried out in StepperRun() once the model integration is
  ! finished
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    implicit none

    type(pflotran_model_type), pointer :: model

    call model%simulation%FinalizeRun()

  end subroutine pflotranModelStepperRunFinalize

! ************************************************************************** !

  subroutine pflotranModelDestroy(model)
  !
  ! Deallocates the pflotranModel object
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    use Factory_PFLOTRAN_module, only : PFLOTRANFinalize
    use Option_module, only : OptionFinalize
    use clm_pflotran_interface_data

    implicit none

    type(pflotran_model_type), pointer :: model

    call model%simulation%FinalizeRun()
    !call model%simulation%Strip()    ! this causes petsc error of seq. fault issue, although doesn't matter.
    if(associated(model%simulation)) deallocate(model%simulation)
    nullify(model%simulation)

    call PFLOTRANFinalize(model%option)
    call OptionFinalize(model%option)

    call CLMPFLOTRANIDataDestroy()

    if (associated(model%map_clm_sub_to_pf_sub)) &
      call MappingDestroy(model%map_clm_sub_to_pf_sub)
    if (associated(model%map_clm_2dtop_to_pf_2dtop)) &
      call MappingDestroy(model%map_clm_2dtop_to_pf_2dtop)
    if (associated(model%map_clm_2dbot_to_pf_2dbot)) &
      call MappingDestroy(model%map_clm_2dbot_to_pf_2dbot)

    if (associated(model%map_pf_sub_to_clm_sub)) &
      call MappingDestroy(model%map_pf_sub_to_clm_sub)
    if (associated(model%map_pf_2dtop_to_clm_2dtop)) &
      call MappingDestroy(model%map_pf_2dtop_to_clm_2dtop)
    if (associated(model%map_pf_2dbot_to_clm_2dbot)) &
      call MappingDestroy(model%map_pf_2dbot_to_clm_2dbot)

    if (associated(model)) deallocate(model)
    nullify(model)

  end subroutine pflotranModelDestroy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
! THE FOLLOWING BLOCKS OF CODES ARE FOR CLM-PFLOTRAN TH & BGC COUPLING
!
!  Soil Domain Dimension (mesh in PFLOTRAN ~ CLM grid+soil-layers)
!  (normally NOT variable, unless specified)
!
! ************************************************************************** !

  subroutine pflotranModelSetSoilDimension(pflotran_model)
  !
  ! Force soil dimension from CLM to PFLOTRAN subsurface domain.

    use Discretization_module
    use Patch_module
    use Field_module
    use Grid_module
    use Grid_Structured_module
    use Option_module
    use Material_module
    use Material_Aux_class

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    use Variables_module, only : VOLUME

    implicit none

#include "geodesic.inc"

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(discretization_type), pointer        :: discretization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization


    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, i, j, k

    PetscScalar, pointer :: dlon_pf_loc(:),dlat_pf_loc(:), dzsoil_pf_loc(:)           ! soil cell length/width/thickness (deg/deg/m) for 3-D PF cells
    PetscScalar, pointer :: lonc_pf_loc(:),latc_pf_loc(:), zisoil_pf_loc(:)           ! soil cell coordinates (deg/deg/m) for 3-D PF cells
    PetscScalar, pointer :: topface_pf_loc(:)
    PetscScalar, pointer :: cellid_clm_loc(:), cellid_pf_loc(:)
    PetscReal            :: lon_c, lat_c, lon_e, lon_w, lat_s, lat_n
    PetscReal            :: x_global, y_global
    PetscReal            :: tempreal

    ! for calling functions in 'geodesic.for'
    double precision a, f, lat1, lon1, azi1, lat2, lon2, azi2, s12,   &
        dummy1, dummy2, dummy3, dummy4, dummy5
    double precision lats(4), lons(4)
    integer omask
    a = 6378137.0d0           ! major-axis length of Earth Ellipsoid in metres in WGS-84
    f = 1.d0/298.257223563d0  ! flatening of Earth Ellipsoid in WGS-84
    omask = 0

    subname = 'pflotranModelSetSoilDimension'

!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select

    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    discretization  => realization%discretization

    ! 2D - top layer grid id
    if(associated(pflotran_model%map_clm_2dtop_to_pf_2dtop)) then
      call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%cellid_2dtop_clmp, &
                                    clm_pf_idata%cellid_2dtop_pfs)
    endif

    ! 3D - cellid
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%cellid_clmp, &
                                    clm_pf_idata%cellid_pfs)
    !
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%dzsoil_clmp, &
                                    clm_pf_idata%dzsoil_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%zisoil_clmp, &
                                    clm_pf_idata%zisoil_pfs)
    call VecGetArrayReadF90(clm_pf_idata%dzsoil_pfs, dzsoil_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayReadF90(clm_pf_idata%zisoil_pfs, zisoil_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
         (.not.option%mapping_files) ) then

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%dxsoil_clmp, &
                                    clm_pf_idata%dxsoil_pfs)

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%dysoil_clmp, &
                                    clm_pf_idata%dysoil_pfs)

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%xsoil_clmp, &
                                    clm_pf_idata%xsoil_pfs)

      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%ysoil_clmp, &
                                    clm_pf_idata%ysoil_pfs)

    !
      call VecGetArrayReadF90(clm_pf_idata%dxsoil_pfs, dlon_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayReadF90(clm_pf_idata%dysoil_pfs, dlat_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayReadF90(clm_pf_idata%xsoil_pfs, lonc_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayReadF90(clm_pf_idata%ysoil_pfs, latc_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    end if

    !
    select case(grid%itype)
      case(STRUCTURED_GRID)
         option%io_buffer = "INFO: CLM column dimension will over-ride PF structured grid mesh."
         call printMsg(option)

         select case(grid%structured_grid%itype)
           case(CARTESIAN_GRID)
             option%io_buffer = "INFO: CLM column dimension will over-ride PF structured CARTESIAN_GRID."
             call printMsg(option)

             if (option%io_rank == option%myrank) then
               write(option%fid_out, &
              '(/," Requested processors and decomposition = ", i5,", npx,y,z= ",3i4)') &
               option%mycommsize, &
               grid%structured_grid%npx, &
               grid%structured_grid%npy, &
               grid%structured_grid%npz
               write(option%fid_out,'(" Actual decomposition: npx,y,z= ",3i4,/)') &
               grid%structured_grid%npx_final, &
               grid%structured_grid%npy_final, &
               grid%structured_grid%npz_final
             endif

           case default
             option%io_buffer = "ERROR: Currently only works on structured  CARTESIAN_GRID mesh."
             call printErrMsg(option)
         end select

      case default
         option%io_buffer = "ERROR: Currently only works on structured grids."
         call printErrMsg(option)
    end select

    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)

      ! hack cell's 3-D dimensions

      ! adjusting (x,y) if runs with 2D CLM grid (usually in lat/lon paired grid)
      if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
           (.not.option%mapping_files) ) then

        call StructGridGetIJKFromGhostedID(grid%structured_grid,ghosted_id,i,j,k)

        !        4--(lonc, latc+1/2*dy)-3              ^
        !       /            |           \             |
        !      /             |            \            |
        !     +--------(lonc, latc)--------+     dyg_local(j)
        !    <--------dxg_local(i)---------->          |
        !   /                |               \         |
        !  /                 |                \        |
        ! 1---------(lonc, latc-1/2*dy)--------2       V

        lon_c = lonc_pf_loc(ghosted_id)
        lat_c = latc_pf_loc(ghosted_id)
        if(dlon_pf_loc(ghosted_id) >0.d0 .and. dlat_pf_loc(ghosted_id) >0.d0) then
          ! the following assumes an isoceles trapezoid grid from CLM unstructured or 1-D gridcells
          ! may have distortion for area estimation, but should be good for distance calculation
          lon_e = lon_c + dlon_pf_loc(ghosted_id)/2.0d0  ! East
          lon_w = lon_c - dlon_pf_loc(ghosted_id)/2.0d0  ! West
          lat_s = lat_c - dlat_pf_loc(ghosted_id)/2.0d0  ! South
          lat_n = lat_c + dlat_pf_loc(ghosted_id)/2.0d0  ! North

        else
          lon_e = lon_c + grid%structured_grid%dxg_local(i)/2.0d0  ! East
          lon_w = lon_c - grid%structured_grid%dxg_local(i)/2.0d0  ! West
          lat_s = lat_c - grid%structured_grid%dyg_local(j)/2.0d0  ! South
          lat_n = lat_c + grid%structured_grid%dyg_local(j)/2.0d0  ! North
        endif

        lats(1) = lat_s
        lons(1) = lon_w
        lats(2) = lat_s
        lons(2) = lon_e
        lats(3) = lat_n
        lons(3) = lon_e
        lats(4) = lat_n
        lons(4) = lon_w

        ! mid-longitudal meteric length of trapezoid (x-axis) -
        lat1 = lat_c
        lon1 = lon_w
        lat2 = lat_c
        lon2 = lon_e
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        grid%structured_grid%dx(ghosted_id) = s12

        ! mid-latitudal meteric height of trapezoid (y-axis) -
        lat1 = lat_s
        lon1 = lon_c
        lat2 = lat_n
        lon2 = lon_c
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        grid%structured_grid%dy(ghosted_id) = s12

        ! some checking
        ! areas of grid (x,y) in meters
        call area(a, f, lats, lons, 4, dummy1, dummy2)
        tempreal = grid%structured_grid%dx(ghosted_id)*grid%structured_grid%dy(ghosted_id)/dummy1
        if (abs(tempreal-1.d0)>1.e-5 .and. option%mapping_files) then
          option%io_buffer = "Warning: remarkably large gaps in grid areas btw two approaches FOR cell: "
          call printMsg(option)
        end if

        ! bottom/top segment line length
        lat1 = lats(1)
        lon1 = lons(1)
        lat2 = lats(2)
        lon2 = lons(2)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = s12

        lat1 = lats(3)
        lon1 = lons(3)
        lat2 = lats(4)
        lon2 = lons(4)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = tempreal+s12
        tempreal = 0.5d0*tempreal/grid%structured_grid%dx(ghosted_id)
        if (abs(tempreal-1.d0)>1.e-5) then   ! mathematically, dx = 0.5*(a+b)
          option%io_buffer = "Warning: remarkably large gaps in longitudal-length FOR a cell: "
          call printMsg(option)
        end if

        ! isoscele side line length
        lat1 = lats(2)
        lon1 = lons(2)
        lat2 = lats(3)
        lon2 = lons(3)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = s12

        lat1 = lats(1)
        lon1 = lons(1)
        lat2 = lats(4)
        lon2 = lons(4)
        call invers(a, f, lat1, lon1, lat2, lon2, s12, azi1, azi2, omask,     &
          dummy1, dummy2, dummy3, dummy4 , dummy5)
        tempreal = tempreal/s12
        if (abs(tempreal-1.d0)>1.e-5) then   ! mathematically, c=d
          option%io_buffer = "Warning: remarkably large gaps in isoscele latitudal-length FOR a cell: "
          call printMsg(option)
        end if

      end if ! if (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1 .and. .not.mapping_files)

      ! vertical (z-axis) (always from CLM to PF)
      grid%structured_grid%dz(ghosted_id) = dzsoil_pf_loc(ghosted_id)
      grid%z(ghosted_id)   = zisoil_pf_loc(ghosted_id)    ! directly over-ride PF 'z' coordinate from CLM soil layer 'zi'
    enddo

    if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
         (.not.option%mapping_files) ) then
      call VecRestoreArrayReadF90(clm_pf_idata%dxsoil_pfs, dlon_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayReadF90(clm_pf_idata%dysoil_pfs, dlat_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayReadF90(clm_pf_idata%xsoil_pfs, lonc_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayReadF90(clm_pf_idata%ysoil_pfs, latc_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    end if

    call VecRestoreArrayReadF90(clm_pf_idata%dzsoil_pfs, dzsoil_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayReadF90(clm_pf_idata%zisoil_pfs, zisoil_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! re-do some dimension calculations after changes above
    call MPI_Barrier(option%mycomm,ierr)

    call GridComputeVolumes(grid, field%volume0, option)      ! cell volumes
    call GridComputeInternalConnect(grid, option)             ! cell internal connection distances
    call PatchProcessCouplers(patch,realization%flow_conditions,             &  ! BC/IC/SrcSink connection (face) areas
                              realization%transport_conditions,              &
                              realization%option)

    ! re-assign updated field%volume0 to material_auxvar%volume
    call DiscretizationGlobalToLocal(discretization,field%volume0, &
                                   field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material, &
                               field%work_loc,VOLUME,ZERO_INTEGER)



    ! --------------
    ! the following variable is directly used in 'sandbox_somdec'
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%zsoil_clmp, &
                                    clm_pf_idata%zsoil_pfs)

    ! the following are 'wtgcell' and 'landfrac' adjusted 'TOP' face area (not yet figured out how to use it for multiple columns in ONE grid)
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%area_top_face_clmp, &
                                    clm_pf_idata%area_top_face_pfs)

    if ( (clm_pf_idata%nxclm_mapped >= 1 .and. clm_pf_idata%nyclm_mapped >= 1) .and. &
         (.not.option%mapping_files) ) then
      ! inactive cells with weighted top-surface area of 0 (i.e. CLM grid of inactive or zero coverage of land)
      ! by setting their 'material' id to -999
      call VecGetArrayReadF90(clm_pf_idata%area_top_face_pfs, topface_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      do ghosted_id = 1, grid%ngmax
        local_id = grid%nG2L(ghosted_id)
        if (ghosted_id <= 0 .or. local_id <= 0) cycle

        if (topface_pf_loc(ghosted_id)<=1.d-20 .and. associated(patch%imat)) then
          patch%imat(ghosted_id) = -999
        endif
      end do

      call VecRestoreArrayReadF90(clm_pf_idata%area_top_face_pfs, topface_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    end if

#ifdef CLM_PF_DEBUG
    ! the following is for checking
    call VecGetArrayF90(clm_pf_idata%cellid_pfp, cellid_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%cellid_pfs, cellid_clm_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      cellid_pf_loc(local_id) = grid%nG2A(ghosted_id)

      if(ghosted_id>0) &
        write(option%myrank+200,*) 'PF: natural id -', option%myrank, ghosted_id, local_id, &
        cellid_pf_loc(local_id),grid%nG2A(ghosted_id), cellid_clm_loc(ghosted_id)

    end do
    call VecRestoreArrayF90(clm_pf_idata%cellid_pfp, cellid_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%cellid_pfs, cellid_clm_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%cellid_pfp, &
                                    clm_pf_idata%cellid_clms)
#endif

  end subroutine pflotranModelSetSoilDimension
! ************************************************************************** !
!
!   (BLANK AS INTENDED)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
! THE FOLLOWING BLOCKS OF CODES ARE NEEDED FOR BOTH CLM-PFLOTRAN TH & BGC COUPLING
!
!  Soil Properties: (1) thermal-hydraulic properties (invariable)
!                   (2) effective porosity (dynamical)
!                   (3) essential states: saturation, pressure, matric potential & temperature (dynamical)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  subroutine pflotranModelSetSoilProp(pflotran_model)
  !
  ! Converts hydraulic properties from CLM units
  ! into PFLOTRAN units.
  !
  ! Author: Gautam Bisht
  ! Date: 10/22/2010
  !
  ! (fmy) 4/28/2014: modifications after updating to pflotran-dev

    use Option_module
    use Realization_Base_class
    use Simulation_Base_class
    use Discretization_module
    use Patch_module
    use Grid_module
    use Field_module
    use Material_module
    use Material_Aux_class
    use Coupler_module
    use Connection_module
    use Init_Subsurface_module

    use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_YZ, PERMEABILITY_XZ, &
                               TORTUOSITY, POROSITY

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use TH_Aux_module
    use Characteristic_Curves_module   ! this is used by TH_module
    use Characteristic_Curves_Base_module
    use Characteristic_Curves_Common_module

    use clm_pflotran_interface_data

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(discretization_type), pointer        :: discretization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field

    class(characteristic_curves_type), pointer:: characteristic_curves
    type(th_auxvar_type), pointer             :: th_auxvars(:), th_auxvars_bc(:), th_auxvars_ss(:)
    type(th_auxvar_type), pointer             :: th_auxvar
    type(TH_parameter_type), pointer          :: TH_parameter
    class(material_auxvar_type), pointer      :: material_auxvars(:)
    class(material_property_type), pointer    :: material_property

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(coupler_type), pointer :: boundary_condition, source_sink
    type(connection_set_type), pointer :: cur_connection_set

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, iconn, sum_connection, material_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
    PetscInt           :: sf_func_type, rpf_func_type, ithrm, iphas, icap
    ! Saturation/Permeability functions currently supported by coupled CLM-PFLOTRAN
    PetscInt, parameter :: VAN_GENUCHTEN = 1   ! not yet (TODO)
    PetscInt, parameter :: BROOKS_COREY = 2
    PetscInt, parameter :: BURDINE = 1
    PetscInt, parameter :: MUALEM = 2          ! not yet (TODO)

    PetscScalar, pointer :: tkwet_pf_loc(:)   ! thermal conductivity at saturated (W/m/K)
    PetscScalar, pointer :: tkdry_pf_loc(:)   ! thermal conductivity - dry (W/m/K)
    PetscScalar, pointer :: tkfrz_pf_loc(:)   ! thermal conductivity - frozen (W/m/K)
    PetscScalar, pointer :: hcapvs_pf_loc(:)   ! volume specific heat capacity for soil particle only (J/m3/K)

    PetscScalar, pointer :: hksat_x_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_pf_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: sucsat_pf_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: bsw_pf_loc(:)     ! Clapp and Hornberger "b"

    PetscScalar, pointer :: bd_dry_pf_loc(:)

    PetscScalar, pointer :: vec_p(:)

    ! -------------------------------------------

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = EARTH_GRAVITY      ! [m/S^2]


    subname = 'pflotranModelSetSoilProp'
!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select

    patch           => realization%patch
    discretization  => realization%discretization
    grid            => patch%grid
    field           => realization%field
    material_auxvars=> patch%aux%Material%auxvars

    select case(option%iflowmode)
      case(TH_MODE)
        th_auxvars      => patch%aux%TH%auxvars
        th_auxvars_bc   => patch%aux%TH%auxvars_bc
        th_auxvars_ss   => patch%aux%TH%auxvars_ss
        th_parameter    => patch%aux%TH%TH_parameter
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if(option%ntrandof.le.0) then
            option%io_buffer =  &
               'Current PFLOTRAN mode not supported by pflotranModelSetSoilProp'
            call printErrMsg(option)
        endif
    end select

    ! ---------------------
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%tkwet_clmp, &
                                    clm_pf_idata%tkwet_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%tkdry_clmp, &
                                    clm_pf_idata%tkdry_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%tkfrz_clmp, &
                                    clm_pf_idata%tkfrz_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%hcvsol_clmp, &
                                    clm_pf_idata%hcvsol_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%hksat_x_clmp, &
                                    clm_pf_idata%hksat_x_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%hksat_y_clmp, &
                                    clm_pf_idata%hksat_y_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%hksat_z_clmp, &
                                    clm_pf_idata%hksat_z_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%sucsat_clmp, &
                                    clm_pf_idata%sucsat_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%bsw_clmp, &
                                    clm_pf_idata%bsw_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%watsat_clmp, &
                                    clm_pf_idata%watsat_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%watfc_clmp, &
                                    clm_pf_idata%watfc_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%bulkdensity_dry_clmp, &
                                    clm_pf_idata%bulkdensity_dry_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%effporosity_clmp, &
                                    clm_pf_idata%effporosity_pfs)

    !
    !---------------------
    !
    call VecGetArrayF90(clm_pf_idata%tkwet_pfs, tkwet_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%tkdry_pfs, tkdry_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%tkfrz_pfs, tkfrz_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%hcvsol_pfs, hcapvs_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%sucsat_pfs,  sucsat_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%bsw_pfs,     bsw_pf_loc,     ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(clm_pf_idata%bulkdensity_dry_pfs,  bd_dry_pf_loc,     ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(field%porosity0, porosity_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    if(option%iflowmode==TH_MODE) then
      ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
      call VecGetArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayF90(field%icap_loc,  icap_loc_p, ierr);
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr);
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! --------------------------------------------------------------------------------------
    ! for all internal grid cells
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (ghosted_id <= 0 .or. local_id <= 0) cycle
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) < 0) cycle    ! imat maybe 0, which causes issue
      endif

      !(TODO) need a better way to generate MVG parameters from CLM inputs

      !F.-M. Yuan: (1) the following IS to pass CLM soil hydraulic data into 'saturation_function';
      !            (2) data-passing IS by from 'ghosted_id' to PF's 'local_id'.
      if(option%iflowmode == TH_MODE) then

        ! TH_MODE now are using 'charateristic_curves' module
        characteristic_curves => patch%  &
            characteristic_curves_array(patch%sat_func_id(ghosted_id))%ptr  ! MUST be in 'ghosted_id' for 'sat_func_id(:)'.

        select type(sf => characteristic_curves%saturation_function)
          !class is(sat_func_VG_type)
             ! not-yet (TODO)

          class is(sat_func_BC_type)
            sf_func_type = BROOKS_COREY

            ! currently BC-Burdine saturation/permisivity function type in PFLOTRAN,
            ! with specified values to match with Clapp-Hornberger Eq. used in CLM biogeophysics

            ! Clapp-Hornberger: soilpsi = sucsat * (-9.81) * (fsattmp)**(-bsw)  ! mm H2O Head --> -pa
            !                   K = Ks*fsattmp**(3+2*bsw)
            !         vs.
            ! BC-Burdine: pc =  (Se**(-1.d0/lambda))/alpha, with Se=(lsat-Sr)/(1-Sr)
            !             relative_perm = Se**power, with power = 3+2/lamda


            sf%alpha  = 1.d0/(9.81d0*sucsat_pf_loc(ghosted_id))
            sf%lambda = 1.d0/bsw_pf_loc(ghosted_id)
            ! A NOTE here:
            ! 'lambda' of < 0.16 (or 'bsw'>6) causes large residual saturation(Sr/pcmax) in SF_BC function.
            ! one unreasonable result of this may be large liq. water saturation under frozen condition
            ! (TODO - it's from high SOM soil layers, implying further work on Pedo-Transfer function for peat
            !  e.g. Letts et al. 2000. Fibric b=2.7, Sr=0.04/0.93;
            !                          Hemic  b=6.1, Sr=0.15/0.88;
            !                          Sapric b=12., Sr=0.22/0.83.
            sf%Sr     = 0.0d0

            if(associated(sf%sat_poly) .and. associated(sf%pres_poly)) then
              call sf%SetupPolynomials(option, 'Error for setup smoothing saturation function')
            endif

          class default
            option%io_buffer = 'Currently ONLY support Brooks_COREY saturation function type' // &
              ' when coupled with CLM.'
              call printErrMsg(option)

        end select

        select type(rpf => characteristic_curves%liq_rel_perm_function)
          !class is(rpf_Mualem_VG_liq_type)
              ! (TODO)

          class is(rpf_Burdine_BC_liq_type)
            rpf_func_type = BURDINE

            rpf%lambda = 1.d0/bsw_pf_loc(ghosted_id)
            ! A NOTE here:
            ! 'lambda' of < 0.16 (or 'bsw'>6) causes large residual saturation(Sr/pcmax) in SF_BC function.
            ! one unreasonable result of this may be large liq. water saturation under frozen condition
            ! (TODO - it's from high SOM soil layers, implying further work on Pedo-Transfer function for peat
            !  e.g. Letts et al. 2000. Fibric b=2.7, Sr=0.04/0.93;
            !                          Hemic  b=6.1, Sr=0.15/0.88;
            !                          Sapric b=12., Sr=0.22/0.83.
            rpf%Sr     = 0.0d0

            if(associated(rpf%poly)) then
              call rpf%SetupPolynomials(option, 'Error for setup smoothing liq. permeability function')
            endif

          class default
            option%io_buffer = 'Currently ONLY support Brooks_COREY-Burdine liq. ' // &
             ' permissivity function type when coupled with CLM.'
            call printErrMsg(option)

        end select

        !
        select case(option%iflowmode)
          case(TH_MODE)
            th_auxvar   => th_auxvars(ghosted_id)

            ithrm = int(ithrm_loc_p(ghosted_id))
            iphas = int(iphase_loc_p(ghosted_id))
            icap  = int(icap_loc_p(ghosted_id))

            th_parameter%alpha(ithrm)    = 1.d0/(9.81d0*sucsat_pf_loc(ghosted_id))
            th_parameter%ckwet(ithrm)    = tkwet_pf_loc(ghosted_id)*option%scale   ! W/m/K --> MW/m/K
            !(note: option%scale multiplier is done in TH.F90: setuppatch(), so it's needed here too)
            th_parameter%ckdry(ithrm)    = tkdry_pf_loc(ghosted_id)*option%scale   ! W/m/K --> MW/m/K
            th_parameter%ckfrozen(ithrm) = tkfrz_pf_loc(ghosted_id)*option%scale   ! W/m/K --> MW/m/K

            th_parameter%dencpr(ithrm)   = hcapvs_pf_loc(ghosted_id)*option%scale  ! J/m3-particle/K --> MJ/m3-particle/K

        end select

      endif

      ! hydraulic conductivity => permissivity IS going to 'field%'
      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      if(option%iflowmode==TH_MODE) then
           ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
           perm_xx_loc_p(local_id) = hksat_x_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
           perm_yy_loc_p(local_id) = hksat_y_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
           perm_zz_loc_p(local_id) = hksat_z_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
      endif

      porosity_loc_p(local_id) = watsat_pf_loc(ghosted_id)

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (watsat):
      !  (turn it on with similar output in clm_pflotran_interfaceMod.F90 and reaction_sandbox_denitrification.F90)
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to PF-internal (field%)-vec 'local_id';
      write(option%myrank+200,*) 'checking pflotran-model prior to set soil properties: ', &
        'rank=',option%myrank, 'ngmax=',grid%ngmax, 'nlmax=',grid%nlmax, &
        'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'pfp_porosity(local_id)=',porosity_loc_p(local_id), &
        'clms_watsat(ghosted_id)=',watsat_pf_loc(ghosted_id)
#endif

      ! material_property updates from CLM
      ! TIP: unlike 'characteristic_curves' above, 'material_property' ID is not directly assigned to patch's grid ID.
      !      So, the above data-passing did NOT modify material_property.
      material_id = patch%imat(ghosted_id)
      material_property => patch%material_property_array(material_id)%ptr

      ! have to do one by one as needed
      material_property%rock_density = &
        bd_dry_pf_loc(ghosted_id)/(1.d0-watsat_pf_loc(ghosted_id))      ! kg soil particle /m3 bulk soils

      material_property%permeability(1,1) = hksat_x_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
      material_property%permeability(2,2) = hksat_y_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
      material_property%permeability(3,3) = hksat_z_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0

      material_property%porosity =  watsat_pf_loc(ghosted_id)

      material_property%specific_heat = hcapvs_pf_loc(ghosted_id)/ &
                                        material_property%rock_density  ! J/m^3-K ==> J/kg rock-K

      material_property%thermal_conductivity_dry = tkdry_pf_loc(ghosted_id)
      material_property%thermal_conductivity_wet = tkwet_pf_loc(ghosted_id)
      material_property%alpha = 1.d0/(9.81d0*sucsat_pf_loc(ghosted_id))
      material_property%thermal_conductivity_frozen = tkfrz_pf_loc(ghosted_id)
      !material_property%alpha_fr = material_property_default%alpha_fr

    enddo


    ! -------------
    call VecRestoreArrayF90(clm_pf_idata%tkwet_pfs, tkwet_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%tkdry_pfs, tkdry_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%tkfrz_pfs, tkfrz_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%hcvsol_pfs, hcapvs_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_pfs,  sucsat_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%bsw_pfs,     bsw_pf_loc,     ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecRestoreArrayF90(clm_pf_idata%bulkdensity_dry_pfs,  bd_dry_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecRestoreArrayF90(field%porosity0, porosity_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    if(option%iflowmode==TH_MODE) then
      ! F.-M. Yuan: without flowmode, the folllowing will throw out segementation fault error
      call VecRestoreArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayF90(field%icap_loc,  icap_loc_p, ierr);
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr);
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    call MPI_Barrier(option%mycomm,ierr)

    ! re-do subsurface assigning material properties due to modification above
    call InitSubsurfAssignMatProperties(realization)


    ! --------------------------------------------------------------------------------------
    ! for all boundary cells already defined and updated above
    ! NOTE: here assumed that boundary cells are ALL or Partial entire grid cells in PF mesh.
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
       if (.not.associated(boundary_condition)) exit
       cur_connection_set => boundary_condition%connection_set

       do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          if (ghosted_id <= 0 .or. local_id <= 0) cycle
          if (patch%imat(ghosted_id) < 0) cycle

          select case(option%iflowmode)
            case(TH_MODE)
              call THAuxVarCopy(th_auxvars(ghosted_id),               &   ! 'th_auxvars' have already updated above
                                th_auxvars_bc(sum_connection), option)


          end select

       enddo
       boundary_condition => boundary_condition%next
    enddo

    ! --------------------------------------------------------------------------------------
    ! for all src/sink cells already defined
    ! NOTE: here assumed that src/sink cells are ALL or Partial entire grid cells in PF mesh.
    source_sink => patch%source_sink_list%first

    sum_connection = 0
    do
      if (.not.associated(source_sink)) exit
      cur_connection_set => source_sink%connection_set

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1

        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (ghosted_id <= 0 .or. local_id <= 0) cycle
        if (patch%imat(ghosted_id) < 0) cycle

        select case(option%iflowmode)
          case(TH_MODE)
            call THAuxVarCopy(th_auxvars(ghosted_id),               &   ! 'th_auxvars' have already updated above
                              th_auxvars_ss(sum_connection), option)
        end select

      enddo
      source_sink => source_sink%next
    enddo

  end subroutine pflotranModelSetSoilProp

! ************************************************************************** !

  subroutine pflotranModelResetSoilPorosityFromCLM(pflotran_model)
  !
  ! Resetting soil porosity in pflotran's internal vecs due to changes from CLM
  ! Note: this is used to adjust porosity of ice from total, when Thermal mode is NOT used in PFLOTRAN
  ! F.-M. YUAN:  4/28/2014

    use Realization_Base_class
    use Discretization_module
    use Patch_module
    use Grid_module
    use Field_module
    use Option_module
    use Material_module
    use Material_Aux_class
    use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, POROSITY

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type


    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(discretization_type), pointer        :: discretization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: porosity0_loc_p(:)    ! this is from 'field%porosity0'
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: unitconv, perm_adj, tempreal

    PetscScalar, pointer :: porosity_pfs_loc(:), porosity_pfp_loc(:)  ! these are from 'clm-pf-idata%'
    PetscScalar, pointer :: hksat_x_pf_loc(:), hksat_y_pf_loc(:), hksat_z_pf_loc(:)
    PetscScalar, pointer :: watsat_pf_loc(:), bsw_pf_loc(:)

    !---------------------------------------------------------------------------------

    subname = 'ModelResetSoilPorosityFromCLM'

    !-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select

    patch           => realization%patch
    discretization  => realization%discretization
    grid            => patch%grid
    field           => realization%field

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%effporosity_clmp, &
                                    clm_pf_idata%effporosity_pfs)
    ! for adjusting porosity
    call VecGetArrayF90(clm_pf_idata%effporosity_pfs,  porosity_pfs_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(field%porosity0, porosity0_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! for adjusting permissivity
    if(option%iflowmode==TH_MODE) then

        unitconv  = 0.001002d0/(998.2d0*EARTH_GRAVITY)/1000.d0    ! from hydraulic conductivity (mmH2O/sec) to permissivity (kg/sec)
        perm_adj  = 1.0d0

        call VecGetArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(clm_pf_idata%bsw_pfs,  bsw_pf_loc,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        call VecGetArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (ghosted_id <= 0 .or. local_id <= 0) cycle
      if (patch%imat(ghosted_id) < 0) cycle    ! imat maybe 0, which causes issue

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (ice-adjusted porosity):
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to 'local_id'
      write(option%myrank+200,*) 'checking pflotran-model prior to resetting porosity:', &
        'rank=',option%myrank, &
        'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'porosity0(local_id)=',porosity0_loc_p(local_id),'adjporo(ghosted_id)=',porosity_pfs_loc(ghosted_id)
#endif

      porosity0_loc_p(local_id) = porosity_pfs_loc(ghosted_id)

      if(option%iflowmode==TH_MODE) then
           ! Ksat is based on actaul porosity, so when porosity is using the effective one, Ksat should be effective as well
           ! This will prevent large hydraulic conductivity in PFLOTRAN when shrinking pore size
           ! because PFLOTRAN uses pressure (saturation) in its rel. perm calculation.
           tempreal = porosity_pfs_loc(ghosted_id)/watsat_pf_loc(ghosted_id)
           perm_adj = tempreal**(2.0d0*bsw_pf_loc(ghosted_id)+3.0d0)        ! assuming shrunk pore as VWC to estimate K, by Clapp-Hornberger Eq.
           perm_adj = max(0.d0, min(perm_adj*perm_adj, 1.0d0))
           perm_xx_loc_p(local_id) = perm_adj*hksat_x_pf_loc(ghosted_id)*unitconv
           perm_yy_loc_p(local_id) = perm_adj*hksat_y_pf_loc(ghosted_id)*unitconv
           perm_zz_loc_p(local_id) = perm_adj*hksat_z_pf_loc(ghosted_id)*unitconv

      endif

    enddo

    call VecRestoreArrayF90(clm_pf_idata%effporosity_pfs,  porosity_pfs_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(field%porosity0, porosity0_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    if(option%iflowmode==TH_MODE) then
        call VecRestoreArrayF90(clm_pf_idata%hksat_x_pfs, hksat_x_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(clm_pf_idata%hksat_y_pfs, hksat_y_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(clm_pf_idata%hksat_z_pfs, hksat_z_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(clm_pf_idata%watsat_pfs,  watsat_pf_loc,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(clm_pf_idata%bsw_pfs,  bsw_pf_loc,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        call VecRestoreArrayF90(field%perm0_xx,  perm_xx_loc_p,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(field%perm0_yy,  perm_yy_loc_p,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(field%perm0_zz,  perm_zz_loc_p,  ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    !
    call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                               field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,ZERO_INTEGER)

    if(option%iflowmode==TH_MODE) then
        call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
        call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
        call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)
        call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)
    endif

  end subroutine pflotranModelResetSoilPorosityFromCLM


! ************************************************************************** !

  subroutine pflotranModelGetSoilPropFromPF(pflotran_model)
  !
  ! Pass soil physical properties from PFLOTRAN to CLM, if needed
  !
  ! Author: Fengming Yuan
  ! Date: 1/30/2014
  !

    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Field_module
    use Option_module

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Characteristic_Curves_module
    use Characteristic_Curves_Base_module
    use Characteristic_Curves_Common_module

    use TH_Aux_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    class(characteristic_curves_type), pointer :: characteristic_curves

    type(th_auxvar_type), pointer :: th_auxvar

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: tempreal, tempreal2

    ! pf internal variables
    PetscReal, pointer :: porosity_loc_p(:)

    ! clm-pf-interface Vecs for PF thermal-hydroloical parameters used in CLM-PFLOTRAN interface
    PetscScalar, pointer :: porosity_loc_pfp(:)  ! soil porosity
    PetscScalar, pointer :: sr_pcwmax_loc_pfp(:) ! soil vwc at max. capillary pressure (note: not 'Sr')
    PetscScalar, pointer :: pcwmax_loc_pfp(:)    ! max. capillary pressure

    character(len=MAXSTRINGLENGTH) :: error_string
    PetscInt :: cur_sat_func_id

    subname = 'pflotranModelGetSoilPropFromPF'

    !-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select

    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field

    call VecGetArrayF90(field%porosity_t,porosity_loc_p,ierr)     ! current porosity (checking ?? 'porosity_t' or 'porosity_tpdt')
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%sr_pcwmax_pfp, sr_pcwmax_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%pcwmax_pfp, pcwmax_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (ghosted_id <= 0 .or. local_id <= 0) cycle
      if (associated(patch%imat)) then
         if (patch%imat(ghosted_id) < 0) cycle    ! imat maybe 0, which causes issue
      endif

      ! PF's porosity
      porosity_loc_pfp(local_id) = porosity_loc_p(local_id)

      ! soil hydraulic properties ID for current cell
      cur_sat_func_id = patch%sat_func_id(ghosted_id)

      !
      if(option%iflowmode==TH_MODE) then
        th_auxvar => patch%aux%TH%auxvars(ghosted_id)
      endif

      ! TH_MODE now are using 'charateristic_curves' module
      characteristic_curves => patch% &
          characteristic_curves_array(cur_sat_func_id)%ptr

      select type(sf => characteristic_curves%saturation_function)
        !class is(sat_func_VG_type)
             ! not-yet (TODO)
        class is(sat_func_BC_type)
          ! PF's limits on soil matrix potential (Capillary pressure)
          pcwmax_loc_pfp(local_id) = sf%pcmax

          ! PF's limits on soil water at pcwmax (NOT: not 'Sr', at which PC is nearly 'inf')
          call sf%Saturation(sf%pcmax, tempreal, tempreal2, option)
          sr_pcwmax_loc_pfp(local_id) = tempreal

        class default
            option%io_buffer = 'Currently ONLY support Brooks_COREY saturation function type' // &
              ' when coupled with CLM.'
         call printErrMsg(option)
      end select

    enddo

    call VecRestoreArrayF90(field%porosity_t,porosity_loc_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%sr_pcwmax_pfp, sr_pcwmax_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%pcwmax_pfp, pcwmax_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%effporosity_pfp, &
                                    clm_pf_idata%effporosity_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%sr_pcwmax_pfp, &
                                    clm_pf_idata%sr_pcwmax_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%pcwmax_pfp, &
                                    clm_pf_idata%pcwmax_clms)

    ! reference pressure
    clm_pf_idata%pressure_reference = option%reference_pressure

  end subroutine pflotranModelGetSoilPropFromPF

! ************************************************************************************* !
  ! This routine Updates T/H drivers (PF global vars) for PFLOTRAN BGC/Flow (TH)-mode,
  ! that is to say, CLM passes soil temperature or moisture (liq. pressure) or both to PF global auxvars,
  ! if either T ('pf_tmode') or H ('pf_hmode') or both NOT invoked in PFLOTRAN.
  !
  subroutine pflotranModelUpdateTHfromCLM(pflotran_model, pf_hmode, pf_tmode)

    use Option_module
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    logical, intent(in)                       :: pf_hmode, pf_tmode

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_auxvars(:)

    PetscErrorCode     :: ierr

    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: soillsat_pf_loc(:), soilisat_pf_loc(:)
    PetscReal, pointer :: soilt_pf_loc(:)
    PetscReal, pointer :: soilpress_pf_loc(:)

    subname = 'pflotranModelUpdateTHfromCLM'

!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select

    patch           => realization%patch
    grid            => patch%grid
    global_auxvars  => patch%aux%Global%auxvars

    ! Save the liq saturation values from CLM to PFLOTRAN, if needed
    if (.not.pf_hmode) then
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soillsat_clmp, &
                                    clm_pf_idata%soillsat_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilpsi_clmp, &
                                    clm_pf_idata%soilpsi_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%press_clmp, &
                                    clm_pf_idata%press_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilliq_clmp, &
                                    clm_pf_idata%soilliq_pfs)

        call VecGetArrayF90(clm_pf_idata%soillsat_pfs, soillsat_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        do ghosted_id=1, grid%ngmax
          local_id=grid%nG2L(ghosted_id)
          if (ghosted_id<=0 .or. local_id <=0) cycle
          !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

#ifdef CLM_PF_DEBUG
       ! F.-M. Yuan: the following check proves DATA-passing from CLM to PF MUST BE done by ghosted_id --> ghosted_id
       ! if passing to 'global_auxvars'
          write(option%myrank+200,*) 'checking pflotran-model 1 (CLM->PF lsat): ', &
              'local_id=',local_id, 'ghosted_id=',ghosted_id, &
              'sat_globalvars(ghosted_id)=',global_auxvars(ghosted_id)%sat(1), &
              'sat_pfs(ghosted_id)=',soillsat_pf_loc(ghosted_id)
#endif
          global_auxvars(ghosted_id)%sat(1)=soillsat_pf_loc(ghosted_id)
          global_auxvars(ghosted_id)%pres(1)=soilpress_pf_loc(ghosted_id)
        enddo

        call VecRestoreArrayF90(clm_pf_idata%soillsat_pfs, soillsat_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        !
        ! for exactly using moisture and other response functions of decomposition from CLM-CN
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%w_scalar_clmp, &
                                    clm_pf_idata%w_scalar_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%o_scalar_clmp, &
                                    clm_pf_idata%o_scalar_pfs)

    endif

    ! Save soil temperature values from CLM to PFLOTRAN, if needed
    if (.not.pf_tmode) then
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilt_clmp, &
                                    clm_pf_idata%soilt_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilice_clmp, &
                                    clm_pf_idata%soilice_pfs)

        call VecGetArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        do ghosted_id=1, grid%ngmax
            local_id=grid%nG2L(ghosted_id)
            if (ghosted_id<=0 .or. local_id<=0) cycle
            !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

            global_auxvars(ghosted_id)%temp=soilt_pf_loc(ghosted_id)
        enddo
        call VecRestoreArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! for exactly using temperature response function of decomposition from CLM-CN
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%t_scalar_clmp, &
                                    clm_pf_idata%t_scalar_pfs)

    endif

  end subroutine pflotranModelUpdateTHfromCLM


! ************************************************************************** !
  subroutine pflotranModelGetSaturationFromPF(pflotran_model)
  !
  ! Extract soil saturation values simulated by
  ! PFLOTRAN in a PETSc vector.
  !
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  !
  ! 4/28/2014: fmy - updates

    use Option_module
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Field_module
    use Global_Aux_module
    use TH_Aux_module

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use TH_module, only : THUpdateAuxVars

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(global_auxvar_type), pointer         :: global_auxvars(:)

    type(TH_auxvar_type), pointer             :: th_auxvars(:)

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: soillsat_pf_p(:)       ! 0 - 1 of porosity
    PetscReal, pointer :: soilisat_pf_p(:)       !
    PetscReal, pointer :: soilliq_pf_p(:)        ! kg/m^3 bulk soil
    PetscReal, pointer :: soilice_pf_p(:)
    PetscReal, pointer :: press_pf_p(:)
    PetscReal, pointer :: soilpsi_pf_p(:)
    PetscReal, pointer :: porosity0_loc_p(:)     ! soil porosity in field%porosity0
    PetscScalar, pointer :: porosity_loc_pfp(:)  ! soil porosity saved in clm-pf-idata
    PetscReal          :: liq_kgm3, ice_kgm3

    subname = 'pflotranModelGetSaturationFromPF'
!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
!-------------------------------------------------------------------------

    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    global_auxvars  => patch%aux%Global%auxvars

    select case(option%iflowmode)
      case (TH_MODE)
         call THUpdateAuxVars(realization)
         th_auxvars => patch%aux%TH%auxvars
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if (option%ntrandof .le. 0) then
            option%io_buffer='pflotranModelGetUpdatedStates ' // &
             'implmentation in this mode is not supported!'

            call printErrMsg(option)
        endif

    end select

    ! Save the saturation/pc/pressure values
    call VecGetArrayF90(clm_pf_idata%soillsat_pfp, soillsat_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%soilliq_pfp, soilliq_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%press_pfp, press_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%soilpsi_pfp, soilpsi_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! save porosity for estimating actual water content from saturation, when needed
    call VecGetArrayF90(field%porosity0,porosity0_loc_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    do local_id=1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      if (ghosted_id <=0 ) cycle

      soillsat_pf_p(local_id)=global_auxvars(ghosted_id)%sat(1)
      press_pf_p(local_id)   =global_auxvars(ghosted_id)%pres(1)

      ! PF's field porosity pass to clm-pf-idata and saved
      porosity_loc_pfp(local_id) = porosity0_loc_p(local_id)

      ! calculate water mass and pass to clm-pf-idata
      liq_kgm3 = global_auxvars(ghosted_id)%den_kg(1) ! water den = kg/m^3
      soilliq_pf_p(local_id) = global_auxvars(ghosted_id)%sat(1)* &
                               porosity0_loc_p(local_id)*liq_kgm3

#ifdef CLM_PF_DEBUG
! F.-M. Yuan: the following check proves DATA-passing from PF to CLM MUST BE done by ghosted_id --> local_id
! if passing from 'global_auxvars'
write(option%myrank+200,*) 'checking pflotran-model 2 (PF->CLM lsat):  ', &
        'local_id=',local_id, 'ghosted_id=',ghosted_id,  &
        'sat_globalvar(ghosted_id)=',global_auxvars(ghosted_id)%sat(1), &
        'idata%sat_pfp(local_id)=',soillsat_pf_p(local_id)
#endif

      if (option%iflowmode == TH_MODE) then
        soilpsi_pf_p(local_id) = -th_auxvars(ghosted_id)%pc
      endif
    enddo
    call VecRestoreArrayF90(clm_pf_idata%soillsat_pfp, soillsat_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%soilliq_pfp, soilliq_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%press_pfp, press_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%soilpsi_pfp, soilpsi_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(field%porosity0,porosity0_loc_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%effporosity_pfp, porosity_loc_pfp, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! mapping to CLM vecs (seq)
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%soillsat_pfp, &
                                    clm_pf_idata%soillsat_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%soilliq_pfp,  &
                                    clm_pf_idata%soilliq_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%press_pfp, &
                                    clm_pf_idata%press_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%soilpsi_pfp, &
                                    clm_pf_idata%soilpsi_clms)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%effporosity_pfp, &
                                    clm_pf_idata%effporosity_clms)


    if (option%iflowmode == TH_MODE .and. &
        option%use_th_freezing) then

      TH_auxvars => patch%aux%TH%auxvars

      call VecGetArrayF90(clm_pf_idata%soilisat_pfp, soilisat_pf_p, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayF90(clm_pf_idata%soilice_pfp, soilice_pf_p, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecGetArrayF90(field%porosity0,porosity0_loc_p,ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (ghosted_id <=0 ) cycle

        soilisat_pf_p(local_id) = TH_auxvars(ghosted_id)%ice%sat_ice

        ice_kgm3 = TH_auxvars(ghosted_id)%ice%den_ice        ! ice den = kg/m^3
        soilice_pf_p(local_id) = TH_auxvars(ghosted_id)%ice%sat_ice* &
                                 porosity0_loc_p(local_id)*ice_kgm3

      enddo
      call VecRestoreArrayF90(clm_pf_idata%soilisat_pfp, soilisat_pf_p, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayF90(clm_pf_idata%soilice_pfp, soilice_pf_p, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecRestoreArrayF90(field%porosity0,porosity0_loc_p,ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                      option, &
                                      clm_pf_idata%soilisat_pfp, &
                                      clm_pf_idata%soilisat_clms)

      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                      option, &
                                      clm_pf_idata%soilice_pfp, &
                                      clm_pf_idata%soilice_clms)
    endif

  end subroutine pflotranModelGetSaturationFromPF

! ************************************************************************** !

  subroutine pflotranModelGetTemperatureFromPF(pflotran_model)
  !
  ! This routine get updated states evoloved by PFLOTRAN.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 5/14/2013
  !

    use Option_module
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use TH_Aux_module

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use TH_module, only : THUpdateAuxVars

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    class(simulation_subsurface_type), pointer :: simulation
    class(realization_subsurface_type), pointer:: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid

    type(global_auxvar_type), pointer         :: global_auxvars(:)
    type(th_auxvar_type), pointer             :: th_auxvars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: soilt_pf_p(:)

    subname = 'ModelGetTemperatureFromPF'
!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
!-------------------------------------------------------------------------

    patch           => realization%patch
    grid            => patch%grid
    global_auxvars  => patch%aux%Global%auxvars

    select case(option%iflowmode)
      case (TH_MODE)
         call THUpdateAuxVars(realization)
         th_auxvars => patch%aux%TH%auxvars
      case default
        !F.-M. Yuan: if no flowmode AND no reactive-transport, then let crash the run
        ! because 'uniform_velocity' with [0, 0, 0] velocity transport doesn't need specific flowmode,
        ! which is used for BGC coupling only (i.e., no TH coupling)
        if (option%ntrandof .le. 0) then
            option%io_buffer='pflotranModelGetUpdatedStates ' // &
             'implmentation in this mode is not supported!'
            call printErrMsg(option)
        endif
    end select

    call VecGetArrayF90(clm_pf_idata%soilt_pfp, soilt_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    do local_id=1,grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (ghosted_id>0) then
        soilt_pf_p(local_id) = global_auxvars(ghosted_id)%temp
      endif
    enddo
    call VecRestoreArrayF90(clm_pf_idata%soilt_pfp, soilt_pf_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%soilt_pfp, &
                                    clm_pf_idata%soilt_clms)

  end subroutine pflotranModelGetTemperatureFromPF
! ************************************************************************** !

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
! THE FOLLOWING BLOCKS OF CODES ARE FOR CLM-PFLOTRAN BGC COUPLING
!
!

  ! ************************************************************************** !
  ! pflotranModelGetRTspecies:
  !  PF RT bgc species Name and index (idof)
  !  Then, all indices is saved in 'clm_pf_idata%' for using in this interface module
  !  so if modification needed, only this subroutine
  !     and the 'ispec_*' and 'name_*' put in 'clm_pflotran_interface_data.F90' are modified.
  !
  subroutine pflotranModelGetRTspecies(pflotran_model)

    use Option_module
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer          :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(patch_type), pointer :: patch
    type(grid_type), pointer  :: grid

    PetscErrorCode  :: ierr
    PetscInt        :: k

    character(len=MAXWORDLENGTH) :: word

    !-------------------------------------------------------------------------
    subname = 'pflotranModelGetRTSpecies'
    !-------------------------------------------------------------------------
    option=> pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    patch => realization%patch
    grid  => patch%grid

    !-------------------------------------------------------------------------

    if (option%ntrandof <= 0) return

    !
    !immobile species for liter and SOM (decomposing pools)

    ! for total HR, Nmin, Nimm, and Nimmp (NOTE: do this 'total' before 'individual')
    word = clm_pf_idata%name_hrim
    clm_pf_idata%ispec_hrimm  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = clm_pf_idata%name_Nmin
    clm_pf_idata%ispec_nmin  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = clm_pf_idata%name_nimm
    clm_pf_idata%ispec_nimm  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

    word = clm_pf_idata%name_nimp
    clm_pf_idata%ispec_nimp  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)


    do k=1, clm_pf_idata%ndecomp_pools

      ! NOTE: the PF soil bgc sandbox 'SomDec' has a naming protocol as following
      ! (1) fixed-CN ratio decomposing pool: only C pool name is defined, while N is not needed;
      ! (2) varying-CN ratio decomposing pool: 2 pool names are defined with ending letter 'C' or 'N'

      if (clm_pf_idata%floating_cn_ratio(k)) then
        word = trim(clm_pf_idata%decomp_pool_name(k)) // "C"
        clm_pf_idata%ispec_decomp_c(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
        !
        if (clm_pf_idata%ispec_decomp_c(k) <= 0) then
          option%io_buffer = 'CLM decomposing pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
          call printErrMsg(option)
        endif

        word = trim(clm_pf_idata%decomp_pool_name(k)) // "N"
        clm_pf_idata%ispec_decomp_n(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

        if (clm_pf_idata%ispec_decomp_n(k) <= 0) then
          option%io_buffer = 'CLM decomposing pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
          call printErrMsg(option)
        endif
      !

      else

        word = trim(clm_pf_idata%decomp_pool_name(k))
        clm_pf_idata%ispec_decomp_c(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

        if (clm_pf_idata%ispec_decomp_c(k) <= 0) then
          option%io_buffer = 'CLM decomposing pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
          call printErrMsg(option)
        endif

        !
        clm_pf_idata%ispec_decomp_n(k) = UNINITIALIZED_INTEGER

      endif

      ! 'hr'
      word = trim(clm_pf_idata%decomp_pool_name(k)) // "CHR"
      clm_pf_idata%ispec_decomp_hr(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
      !
      if (clm_pf_idata%ispec_decomp_hr(k) <= 0 .and. clm_pf_idata%ispec_hrimm <= 0 ) then
        option%io_buffer = 'CLM decomposing pool HR: ' // &
          trim(word) // &
          ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
        call printErrMsg(option)
      endif

      ! 'nmin'
      word = trim(clm_pf_idata%decomp_pool_name(k)) // "NMIN"
      clm_pf_idata%ispec_decomp_nmin(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
      !
      if (clm_pf_idata%ispec_decomp_nmin(k) <= 0 .and. clm_pf_idata%ispec_nmin <= 0 ) then
        option%io_buffer = 'CLM decomposing pool NMIN: ' // &
          trim(word) // &
          ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
        call printErrMsg(option)
      endif

      ! 'nimm'
      word = trim(clm_pf_idata%decomp_pool_name(k)) // "NIMM"
      clm_pf_idata%ispec_decomp_nimm(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
      !
      if (clm_pf_idata%ispec_decomp_nimm(k) <= 0 .and. clm_pf_idata%ispec_nimm <= 0 ) then
        option%io_buffer = 'CLM decomposing pool NIMM: ' // &
          trim(word) // &
          ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
        call printErrMsg(option)
      endif

      ! 'nimp'
      word = trim(clm_pf_idata%decomp_pool_name(k)) // "NIMP"
      clm_pf_idata%ispec_decomp_nimp(k) = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
      !
      if (clm_pf_idata%ispec_decomp_nimp(k) <= 0 .and. clm_pf_idata%ispec_nimp <= 0 ) then
        option%io_buffer = 'CLM decomposing pool NIMP: ' // &
          trim(word) // &
          ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
        call printErrMsg(option)
      endif

    end do

    ! aq. species in soil solution
    word = clm_pf_idata%name_co2aq
    clm_pf_idata%ispec_co2aq  = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_co2aq<=0) then
      option%io_buffer = 'CLM co2 (aq) species ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_n2aq
    clm_pf_idata%ispec_co2aq  = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_co2aq<=0) then
      option%io_buffer = 'CLM n2 (aq) species ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_n2oaq
    clm_pf_idata%ispec_co2aq  = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_co2aq<=0) then
      option%io_buffer = 'CLM n2o (aq) species ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_nh4
    clm_pf_idata%ispec_nh4  = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_nh4 > 0) then
      word = clm_pf_idata%name_nh4sorb
      clm_pf_idata%ispec_nh4sorb = GetImmobileSpeciesIDFromName(word, &                        ! for using sandbox of absorption
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    else
      word = clm_pf_idata%name_nh4s
      clm_pf_idata%ispec_nh4s = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

      if (clm_pf_idata%ispec_nh4s <= 0) then
        option%io_buffer = 'CLM N pool ' // &
          ' either, ' // trim(clm_pf_idata%name_nh4) // &
          ' or,' // trim(clm_pf_idata%name_nh4s) // &
          ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
        call printErrMsg(option)
      endif
    endif

    word = clm_pf_idata%name_no3
    clm_pf_idata%ispec_no3    = GetPrimarySpeciesIDFromName(word, &
                  realization%reaction,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_no3 <=0 ) then
      word = clm_pf_idata%name_no3s
      clm_pf_idata%ispec_no3s = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)

      if (clm_pf_idata%ispec_no3s <= 0) then
        option%io_buffer = 'CLM N pool ' // &
          ' either, ' // trim(clm_pf_idata%name_no3) // &
          ' or,' // trim(clm_pf_idata%name_no3s) // &
          ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
        call printErrMsg(option)
      endif
    endif
    !
    ! species for gases (now as immobile species)
    word = clm_pf_idata%name_co2
    clm_pf_idata%ispec_co2  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_co2<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_n2
    clm_pf_idata%ispec_n2  = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_n2<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_n2o
    clm_pf_idata%ispec_n2o = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_n2o<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    ! other tracking variables
    word = clm_pf_idata%name_plantndemand
    clm_pf_idata%ispec_plantndemand = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_plantndemand<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printMsg(option)         ! ONLY warning
    endif

    word = clm_pf_idata%name_plantno3uptake
    clm_pf_idata%ispec_plantno3uptake = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_plantno3uptake<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_plantnh4uptake
    clm_pf_idata%ispec_plantnh4uptake = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_plantnh4uptake<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    !
    word = clm_pf_idata%name_ngasmin
    clm_pf_idata%ispec_ngasmin = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_ngasmin<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_ngasnitr
    clm_pf_idata%ispec_ngasnitr = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_ngasnitr<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

    word = clm_pf_idata%name_ngasdeni
    clm_pf_idata%ispec_ngasdeni = GetImmobileSpeciesIDFromName(word, &
                  realization%reaction%immobile,PETSC_FALSE,realization%option)
    if (clm_pf_idata%ispec_ngasdeni<=0) then
      option%io_buffer = 'CLM-PF bgc pool ' // &
            trim(word) // &
            ' in PFLOTRAN_CLM_MAIN interface not found in list of PF chemical species pools.'
      call printErrMsg(option)
    endif

  end subroutine pflotranModelGetRTspecies

  ! *********************************************************************************** !
  ! This routine Pass CLM SOM pools and decomposition rate constants for PFLOTRAN bgc
  ! So that both are consistent
  !
  ! @author
  ! F.-M. Yuan
  !
  ! Date: 12/5/2014
  !
  subroutine pflotranModelSetSOMKfromCLM(pflotran_model)

    use Option_module
    use Realization_Base_class
    use Reaction_module
    use Reaction_Sandbox_module
    use Reaction_Sandbox_Base_class
    use Reaction_Sandbox_SomDec_class

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer          :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    class(reaction_sandbox_base_type), pointer  :: cur_rtsandbox
    class(reaction_sandbox_somdec_type), pointer:: rtsandbox_somdec

    !
    PetscInt :: irxn, jdown, ki, kj
    PetscReal:: sum_cfrac, sum_nfrac

    !-------------------------------------------------------------------------
    subname = 'pflotranModelSetSOMKformCLM'
    !-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    !-------------------------------------------------------------------------

    if (associated(rxn_sandbox_list)) then   ! note: 'rxn_sandbox_list' IS PUBLIC in 'Reaction_module'
      cur_rtsandbox => rxn_sandbox_list
      do
        if (.not.associated(cur_rtsandbox)) exit

        select type (cur_rtsandbox)
          class is (reaction_sandbox_somdec_type)
            rtsandbox_somdec => cur_rtsandbox

            !
            do irxn = 1, rtsandbox_somdec%nrxn

              do ki = 1, clm_pf_idata%ndecomp_pools
                if (clm_pf_idata%ispec_decomp_c(ki) == &
                     rtsandbox_somdec%upstream_c_id(irxn)) then
                   exit
                else
                   if (ki==clm_pf_idata%ndecomp_pools) then   ! didn't find the decomposing pool
                     option%io_buffer &
                       = "ERROR: CLM decomposing pool NOT found in sandbox of SOMDECOMP"
                     call printErrMsg(option)
                   end if
                end if
              end do

              rtsandbox_somdec%rate_constant(irxn)      = clm_pf_idata%ck_decomp_c(ki)
              rtsandbox_somdec%rate_decomposition(irxn) = UNINITIALIZED_DOUBLE
              rtsandbox_somdec%rate_ad_factor(irxn)     = clm_pf_idata%adfactor_ck_c(ki)
              rtsandbox_somdec%upstream_is_varycn(irxn) = clm_pf_idata%floating_cn_ratio(ki)
              rtsandbox_somdec%upstream_nc(irxn)        = clm_pf_idata%decomp_element_ratios(ki,2) &
                                                         /clm_pf_idata%decomp_element_ratios(ki,1)
              rtsandbox_somdec%mineral_c_stoich(irxn)   = clm_pf_idata%fr_decomp_c(ki,ki)  ! this is the decomposed C fraction as respirted CO2)

              sum_cfrac = 0.d0
              sum_nfrac = 0.d0
              do jdown = 1, rtsandbox_somdec%n_downstream_pools(irxn)
                do kj = 1, clm_pf_idata%ndecomp_pools
                  if (clm_pf_idata%ispec_decomp_c(kj) == &
                    rtsandbox_somdec%downstream_c_id(irxn,jdown)) then
                    exit
                  else
                    if (kj==clm_pf_idata%ndecomp_pools) then   ! didn't find the decomposing downstream pool
                      option%io_buffer &
                        = "ERROR: CLM decomposing downstream pool NOT found in sandbox of SOMDECOMP, for: " // &
                          trim(clm_pf_idata%decomp_pool_name(ki))
                      call printErrMsg(option)
                    end if
                  end if
                end do

                rtsandbox_somdec%downstream_is_varycn(irxn, jdown) = clm_pf_idata%floating_cn_ratio(kj)
                ! note: the following is the initial NC ratios, which for 'varying_cn' species will be modified
                rtsandbox_somdec%downstream_nc(irxn, jdown)        = clm_pf_idata%decomp_element_ratios(kj,2) &
                                                                    /clm_pf_idata%decomp_element_ratios(kj,1)
                rtsandbox_somdec%downstream_stoich(irxn, jdown)    = clm_pf_idata%fr_decomp_c(ki,kj)
                sum_cfrac = sum_cfrac + rtsandbox_somdec%downstream_stoich(irxn, jdown)
                sum_nfrac = sum_nfrac + rtsandbox_somdec%downstream_stoich(irxn, jdown) * &
                                        rtsandbox_somdec%downstream_nc(irxn, jdown)

              end do

              if (abs(sum_cfrac+rtsandbox_somdec%mineral_c_stoich(irxn)-1.d0)>1.0d-30) then
                option%io_buffer &
                  = "ERROR: fraction of CLM decomposing downstream pools NOT summed to 1.0, for: " // &
                          trim(clm_pf_idata%decomp_pool_name(ki))
                call printErrMsg(option)
              else
                rtsandbox_somdec%mineral_c_stoich(irxn) = 1.0d0 - sum_cfrac    ! just in case (may not be needed)

                rtsandbox_somdec%mineral_n_stoich(irxn) = rtsandbox_somdec%upstream_nc(irxn) - sum_nfrac

              endif

            end do

            ! finally, exit the do-loop expcilitly so that NO duplicate sandbox checking
            exit

          class default
            if (.not.associated(cur_rtsandbox%next)) then
              option%io_buffer &
                = "ERROR: SetSOMK from CLM currently only supported sandbox of SOMDECOMP"
              call printErrMsg(option)
            end if

        end select


        cur_rtsandbox => cur_rtsandbox%next
      enddo

    else
      option%io_buffer &
          = "ERROR: SetSOMK from CLM currently only supported sandbox approach"
      call printErrMsg(option)

    end if

  end subroutine pflotranModelSetSOMKfromCLM

  ! ************************************************************************** !
  !
  ! pflotranModelSetBgcConc:
  ! Get CLM concentrations (C, N) and set them into PFLOTRAN bgc species
  !
  subroutine pflotranModelSetBgcConcFromCLM(pflotran_model)

    use Global_Aux_module
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer          :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(field_type), pointer                   :: field
    type(patch_type), pointer                   :: patch
    type(grid_type), pointer                    :: grid

    type(global_auxvar_type), pointer :: global_auxvars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id
    PetscInt           :: ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: decomp_cpools_vr_pf_loc(:)      ! (molesC/m3)
    PetscScalar, pointer :: decomp_npools_vr_pf_loc(:)      ! (molesN/m3)
    PetscScalar, pointer :: smin_no3_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4sorb_vr_pf_loc(:)       ! (molesN/m3)

    PetscReal, pointer :: porosity_loc_p(:)

    PetscInt  :: offset, offsetim
    PetscReal :: porosity, saturation, theta ! for concentration conversion from mol/m3 to mol/L
    PetscReal :: xmass, den_kg_per_L         ! for from mol/L to mol/kg

    Vec                  :: vec_clmp
    Vec                  :: vec_pfs
    PetscScalar, pointer :: array_clmp(:), array_pfs(:), array_temp(:)
    PetscInt             :: j, k, vec_offset

    !-------------------------------------------------------------------------
    subname = 'pflotranModelGetBgcConcFromCLM'
    !-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    patch => realization%patch
    grid  => patch%grid
    field => realization%field

    global_auxvars  => patch%aux%Global%auxvars

    !-----------------------------------------------------------------

    ! create temporary vecs/arrays for each 'decomp_pool' data-mapping
    call VecDuplicate(clm_pf_idata%zsoil_clmp, vec_clmp,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDuplicate(clm_pf_idata%zsoil_pfs, vec_pfs,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! decomp'C'
    if (associated(clm_pf_idata%ispec_decomp_c)) then
      ! assembly the 'vec_clmp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%decomp_cpools_vr_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_cpools_vr_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_clmp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlclm_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'

        call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_clmp, array_temp, ierr)
        call VecGetArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_clmp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlclm_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_clmp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec for one species
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option,                    &
                                    vec_clmp,                                 &
                                    vec_pfs)

        ! insert 'vec_pfs' into the whole '_pfs' vec
        vec_offset = (k-1)*clm_pf_idata%ngpf_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_pfs, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngpf_sub) = array_pfs

        call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_pfs, array_temp, ierr)
        call VecRestoreArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_pfs' vec
      call VecAssemblyBegin(clm_pf_idata%decomp_cpools_vr_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_cpools_vr_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! decomp_'N'
    if (associated(clm_pf_idata%ispec_decomp_n)) then
      call VecAssemblyBegin(clm_pf_idata%decomp_npools_vr_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_npools_vr_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a segment of data from the whole '_clmp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlclm_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_clmp, array_temp, ierr)
        call VecGetArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_clmp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlclm_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_clmp, array_temp, ierr)
        call VecRestoreArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. Vec
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option,                    &
                                    vec_clmp,                                 &
                                    vec_pfs)

        ! insert 'vec_pfs' into the whole '_pfs' vec
        vec_offset = (k-1)*clm_pf_idata%ngpf_sub       ! Seq. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_pfs, array_temp, ierr)
        call VecGetArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngpf_sub) = array_pfs

        call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_pfs, array_temp, ierr)
        call VecRestoreArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_pfs' vec, which then pass to PF's field%
      call VecAssemblyBegin(clm_pf_idata%decomp_npools_vr_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_npools_vr_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! clear-up of temporary vecs/arrarys
    call VecDestroy(vec_clmp,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDestroy(vec_pfs,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)


    !-----------------------------------------------------------------

    if(clm_pf_idata%ispec_no3 > 0 .or. clm_pf_idata%ispec_no3s > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%smin_no3_vr_clmp, &
                                    clm_pf_idata%smin_no3_vr_pfs)
    endif

    if(clm_pf_idata%ispec_nh4 > 0 .or. clm_pf_idata%ispec_nh4s > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%smin_nh4_vr_clmp, &
                                    clm_pf_idata%smin_nh4_vr_pfs)

       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%smin_nh4sorb_vr_clmp, &
                                    clm_pf_idata%smin_nh4sorb_vr_pfs)
    endif

    !----------------------------------------------------------------------------

    if (associated(clm_pf_idata%ispec_decomp_c)) then
      call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_pfs, decomp_cpools_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif
    if (associated(clm_pf_idata%ispec_decomp_n)) then
      call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_pfs, decomp_npools_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    !
    if(clm_pf_idata%ispec_no3 > 0 .or. clm_pf_idata%ispec_nh4s > 0) then
      call VecGetArrayReadF90(clm_pf_idata%smin_no3_vr_pfs, smin_no3_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if(clm_pf_idata%ispec_nh4 > 0 .or. clm_pf_idata%ispec_nh4s > 0) then
      call VecGetArrayReadF90(clm_pf_idata%smin_nh4_vr_pfs, smin_nh4_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      call VecGetArrayReadF90(clm_pf_idata%smin_nh4sorb_vr_pfs, smin_nh4sorb_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayReadF90(field%porosity0, porosity_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)


    !----------------------------------------------------------------------------------------------------
    do ghosted_id=1, grid%ngmax
      local_id=grid%nG2L(ghosted_id)
      if (ghosted_id<=0 .or. local_id <= 0) cycle ! bypass ghosted corner cells
      !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

      offset = (local_id - 1)*realization%reaction%ncomp

      saturation = global_auxvars(ghosted_id)%sat(1)          ! using 'ghosted_id' if from 'global_auxvars'
      porosity = porosity_loc_p(local_id)                     ! using 'local_id' if from 'field%???'
      theta = saturation * porosity

      ! needs the following (water density) to do convertion (esp. if it's upon estimating by EOSwater module in PF)
      ! 'xx_p' is in molality (moles/kgH2o) for aq. species, while 'clm_pf_idata%???' is in moles/m3.
      ! but they're same for immobile species
      xmass = 1.d0
      if (associated(global_auxvars(ghosted_id)%xmass)) then
        xmass = global_auxvars(ghosted_id)%xmass(1)
      endif
      den_kg_per_L = global_auxvars(ghosted_id)%den_kg(1)*xmass*1.d-3

      if(clm_pf_idata%ispec_no3 > 0) then
         xx_p(offset + clm_pf_idata%ispec_no3) = max(xeps0_n,                     &      ! from 'ghosted_id' to field%xx_p's local
                                                smin_no3_vr_pf_loc(ghosted_id)    &
                                                / (theta*1000.d0*den_kg_per_L) )    ! moles/m3 /(m3/m3 * L/m3 * kg/L) = moles/kgh2o

      elseif (clm_pf_idata%ispec_no3s > 0) then  ! in immobile form
         xx_p(offsetim + clm_pf_idata%ispec_no3s) = max(xeps0_n,                  &
                                        smin_no3_vr_pf_loc(ghosted_id) )

      endif

      if(clm_pf_idata%ispec_nh4 > 0) then
         xx_p(offset + clm_pf_idata%ispec_nh4) = max(xeps0_n,                     &
                                                smin_nh4_vr_pf_loc(ghosted_id)    &
                                                / (theta*1000.d0*den_kg_per_L) )

      elseif (clm_pf_idata%ispec_nh4s > 0) then  ! in immobile form
         xx_p(offsetim + clm_pf_idata%ispec_nh4s) = max(xeps0_n,                  &
                                        smin_nh4_vr_pf_loc(ghosted_id) )

      endif

      !
      offsetim = offset + realization%reaction%offset_immobile

      if(clm_pf_idata%ispec_nh4sorb > 0) then   ! for absorbed NH4 as immobile species used in sandbox of absorption
         xx_p(offsetim + clm_pf_idata%ispec_nh4sorb) = max(xeps0_n,               &
                                        smin_nh4sorb_vr_pf_loc(ghosted_id) )
      endif

      !
      do k=1,clm_pf_idata%ndecomp_pools
        vec_offset = (k-1)*clm_pf_idata%ngpf_sub       ! Seq. decomp_pfs vec: 'cell' first, then 'species'

        if(clm_pf_idata%ispec_decomp_c(k) > 0) then
          xx_p(offsetim + clm_pf_idata%ispec_decomp_c(k)) = max( xeps0_c,    &               ! field%tran_xx vec IS arranged 'species' first and then 'cell'
                       decomp_cpools_vr_pf_loc(vec_offset+ghosted_id) )      ! Seq. decomp_pfs vec: 'cell' first, then 'species'
        endif

        if(clm_pf_idata%ispec_decomp_n(k) > 0) then
          xx_p(offsetim + clm_pf_idata%ispec_decomp_n(k)) = max( xeps0_n,    &               ! field%tran_xx vec IS arranged 'species' first and then 'cell'
                       decomp_npools_vr_pf_loc(vec_offset+ghosted_id) )      ! Seq. decomp_pfs vec: 'cell' first, then 'species'
        endif


#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (som4c pool):
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to 'local_id'
        if (k==8) then
          write(option%myrank+200,*) 'checking bgc - pflotran-model setting init. conc.: '
          write(option%myrank+200,*) 'rank=',option%myrank, &
          'local_id=',local_id, 'ghosted_id=',ghosted_id, &
          'xxp_som4_id', (offsetim+clm_pf_idata%ispec_decomp_c(k)), &
          'som4_pfs(ghosted_id)=', decomp_cpools_vr_pf_loc(vec_offset+ghosted_id), &
          'xx_p(xxp_som4_id)=',xx_p(offsetim + clm_pf_idata%ispec_decomp_c(k))
        endif
#endif

      enddo  ! k=1,clm_pf_idata%ndecomp_pools

    enddo

    !----------------------------------------------------------------------------------------------------

    if (associated(clm_pf_idata%ispec_decomp_c)) then
      call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_pfs, decomp_cpools_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif
    if (associated(clm_pf_idata%ispec_decomp_n)) then
      call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_pfs, decomp_npools_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    !
    if(clm_pf_idata%ispec_no3 > 0 .or. clm_pf_idata%ispec_no3s > 0) then
      call VecRestoreArrayReadF90(clm_pf_idata%smin_no3_vr_pfs, smin_no3_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if(clm_pf_idata%ispec_nh4 > 0 .or. clm_pf_idata%ispec_nh4s > 0) then
      call VecRestoreArrayReadF90(clm_pf_idata%smin_nh4_vr_pfs, smin_nh4_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      call VecRestoreArrayReadF90(clm_pf_idata%smin_nh4sorb_vr_pfs, smin_nh4sorb_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(field%porosity0, porosity_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  end subroutine pflotranModelSetBgcConcFromCLM

  ! ************************************************************************** !
  ! pflotranModelSetBGCRatesFromCLM:
  !  Get CLM litter, som, mineral N production and plant demand rates
  !  Convert from CLM units into PFLOTRAN units.
  !  Set values in PFLOTRAN
  ! ************************************************************************** !
  subroutine pflotranModelSetBGCRatesFromCLM(pflotran_model)

    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module
    use Data_Mediator_Base_class, only : data_mediator_base_type
    use Data_Mediator_Dataset_class, only : data_mediator_dataset_type
    use Data_Mediator_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer  :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(field_type), pointer           :: field
    type(patch_type), pointer           :: patch
    type(grid_type), pointer            :: grid

    class(data_mediator_base_type), pointer :: cur_data_mediator

    PetscErrorCode      :: ierr
    PetscInt            :: ghosted_id,local_id
    PetscInt            :: offset, offsetim

    Vec                  :: vec_clmp
    Vec                  :: vec_pfs
    PetscScalar, pointer :: array_clmp(:), array_pfs(:), array_temp(:)
    PetscInt             :: j, k, vec_offset
    PetscBool            :: found_rtmasstr

    PetscScalar, pointer :: rate_pf_loc(:)   !
    PetscReal,   pointer :: volume_p(:)

    !-------------------------------------------------------------------------
    subname = 'pflotranModelSetBgcRatesFromCLM'
    !-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    patch => realization%patch
    grid  => patch%grid
    field => realization%field

    !-----------------------------------------------------------------

    ! create temporary vecs/arrays for rate of each 'decomp_pool''s data-mapping
    call VecDuplicate(clm_pf_idata%zsoil_clmp, vec_clmp,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDuplicate(clm_pf_idata%zsoil_pfs, vec_pfs,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! rate of decomp_'C' pools from CLM-CN
    if (associated(clm_pf_idata%ispec_decomp_c)) then
      ! assembly the 'vec_clmp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%rate_decomp_c_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%rate_decomp_c_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_clmp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlclm_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%rate_decomp_c_clmp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_clmp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlclm_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%rate_decomp_c_clmp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option,                    &
                                    vec_clmp,                                 &
                                    vec_pfs)

        ! insert 'vec_pfs' into the whole '_pfs' vec
        vec_offset = (k-1)*clm_pf_idata%ngpf_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%rate_decomp_c_pfs, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngpf_sub) = array_pfs

        call VecRestoreArrayF90(clm_pf_idata%rate_decomp_c_pfs, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_pfs' vec
      call VecAssemblyBegin(clm_pf_idata%rate_decomp_c_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%rate_decomp_c_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! rate of decomp_'N' pools from CLM-CN
    if (associated(clm_pf_idata%ispec_decomp_n)) then

      call VecAssemblyBegin(clm_pf_idata%rate_decomp_n_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%rate_decomp_n_clmp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

      do k=1,clm_pf_idata%ndecomp_pools
        ! get a segment of data from the whole '_clmp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlclm_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%rate_decomp_n_clmp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_clmp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlclm_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%rate_decomp_n_clmp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_clmp, array_clmp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. Vec
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option,                    &
                                    vec_clmp,                                 &
                                    vec_pfs)

        ! insert 'vec_pfs' into the whole '_pfs' vec
        vec_offset = (k-1)*clm_pf_idata%ngpf_sub       ! Seq. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%rate_decomp_n_pfs, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngpf_sub) = array_pfs

        call VecRestoreArrayF90(clm_pf_idata%rate_decomp_n_pfs, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_pfs, array_pfs, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_pfs' vec, which then pass to PF's field%
      call VecAssemblyBegin(clm_pf_idata%rate_decomp_n_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%rate_decomp_n_pfs, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! clear-up of temporary vecs/arrarys
    call VecDestroy(vec_clmp,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDestroy(vec_pfs,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !-----------------------------------------------------------------
    ! the following is a site-scalar to adjust decomposition rate constants
    ! due to its dependency upon 'time', it must be called each CLM time-step.
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%kscalar_decomp_c_clmp, &
                                    clm_pf_idata%kscalar_decomp_c_pfs)


    !-----------------------------------------------------------------

    ! NOTE: direct data passing from interface to PF for N demand
    if(clm_pf_idata%ispec_plantndemand >0) then
      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%rate_plantndemand_clmp, &
                                    clm_pf_idata%rate_plantndemand_pfs)
    endif

    if(clm_pf_idata%ispec_no3 >0 .or. clm_pf_idata%ispec_no3s > 0) then
      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%rate_smin_no3_clmp, &
                                    clm_pf_idata%rate_smin_no3_pfs)
    endif

    if(clm_pf_idata%ispec_nh4 >0 .or. clm_pf_idata%ispec_no3s > 0) then
      call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%rate_smin_nh4_clmp, &
                                    clm_pf_idata%rate_smin_nh4_pfs)
    endif

    !----------------------------------------------------------------------

    ! get cell volume to convert mass transfer rate unit from moles/m3/s to moles/s
    call VecGetArrayReadF90(field%volume0,volume_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    offsetim = realization%reaction%offset_immobile

    cur_data_mediator => realization%tran_data_mediator_list
    do
      if (.not.associated(cur_data_mediator)) exit


      select type (cur_data_mediator)
        class is (data_mediator_dataset_type)
          !
          found_rtmasstr = PETSC_FALSE
          vec_offset = 0
          !
          if(cur_data_mediator%idof == clm_pf_idata%ispec_nh4 .or. &
             cur_data_mediator%idof == clm_pf_idata%ispec_nh4s) then
            call VecGetArrayReadF90(clm_pf_idata%rate_smin_nh4_pfs, rate_pf_loc, ierr)
            call where_checkerr(ierr, subname, __FILE__, __LINE__)
            found_rtmasstr = PETSC_TRUE

          elseif(cur_data_mediator%idof == clm_pf_idata%ispec_no3 .or. &
             cur_data_mediator%idof == clm_pf_idata%ispec_no3s) then
            call VecGetArrayReadF90(clm_pf_idata%rate_smin_no3_pfs, rate_pf_loc, ierr)
            call where_checkerr(ierr, subname, __FILE__, __LINE__)
            found_rtmasstr = PETSC_TRUE

          else
          !--------------------------------------------------------------------------
            if (associated(clm_pf_idata%ispec_decomp_c) .or. &
               associated(clm_pf_idata%ispec_decomp_n)) then
             do k=1,clm_pf_idata%ndecomp_pools
               if( cur_data_mediator%idof == (offsetim + clm_pf_idata%ispec_decomp_c(k)) ) then
                 vec_offset = (k-1)*clm_pf_idata%ngpf_sub       ! Seq. decomp_pfs vec: 'cell' first, then 'species'
                 call VecGetArrayReadF90(clm_pf_idata%rate_decomp_c_pfs, rate_pf_loc, ierr)
                 call where_checkerr(ierr, subname, __FILE__, __LINE__)
                 found_rtmasstr = PETSC_TRUE

                 exit   ! exit the 'do k=1, ...' loop

               elseif( cur_data_mediator%idof == (offsetim + clm_pf_idata%ispec_decomp_n(k)) ) then
                 vec_offset = (k-1)*clm_pf_idata%ngpf_sub       ! Seq. decomp_pfs vec: 'cell' first, then 'species'
                 call VecGetArrayReadF90(clm_pf_idata%rate_decomp_n_pfs, rate_pf_loc, ierr)
                 call where_checkerr(ierr, subname, __FILE__, __LINE__)
                 found_rtmasstr = PETSC_TRUE

                 exit   ! exit the 'do k=1, ...' loop
               endif
             enddo
           endif
         !--------------------------------------------------------------------------
         endif

         cur_data_mediator%dataset%rarray(:) = 0.d0

         ! found a RT mass-transfer variable in the list
         if(found_rtmasstr) then

           do local_id = 1, grid%nlmax
             ghosted_id = grid%nL2G(local_id)
             if (local_id<=0 .or. ghosted_id<=0) cycle
             !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

#ifdef CLM_PF_DEBUG
      !F.-M. Yuan: the following IS a checking, comparing CLM passed data (mass transfer rate):
      ! Conclusions: (1) local_id runs from 1 ~ grid%nlmax; and ghosted_id is obtained by 'nL2G' as corrected above;
      !              OR, ghosted_id runs from 1 ~ grid%ngmax; and local_id is obtained by 'nG2L'.
      !              (2) data-passing IS by from 'ghosted_id' to 'local_id'
      if (cur_data_mediator%idof == clm_pf_idata%ispec_nh4) &
      write(option%myrank+200,*) 'checking bgc-mass-rate - pflotran_model: ', &
        'rank=',option%myrank, 'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'rate_nh4_pfs(ghosted_id)=',rate_pf_loc(ghosted_id), &
        'masstransfer_nh4_predataset(local_id)=',cur_data_mediator%dataset%rarray(local_id)/volume_p(local_id)
#endif

               cur_data_mediator%dataset%rarray(local_id) = &
                         rate_pf_loc(vec_offset+ghosted_id) &    ! for 'decomp' species, must be offset the vecs' starting location
                       * volume_p(local_id)                      ! from moles/m3/s --> moles/s

           enddo  !do local_id = 1, grid%nlmax

           ! close the open vecs
           if(cur_data_mediator%idof == clm_pf_idata%ispec_nh4 .or. &
              cur_data_mediator%idof == clm_pf_idata%ispec_nh4s) then
             call VecRestoreArrayReadF90(clm_pf_idata%rate_smin_nh4_pfs, rate_pf_loc, ierr)
             call where_checkerr(ierr, subname, __FILE__, __LINE__)
           elseif(cur_data_mediator%idof == clm_pf_idata%ispec_no3 .or. &
                  cur_data_mediator%idof == clm_pf_idata%ispec_no3s) then
             call VecRestoreArrayReadF90(clm_pf_idata%rate_smin_no3_pfs, rate_pf_loc, ierr)
             call where_checkerr(ierr, subname, __FILE__, __LINE__)

           else
           !--------------------------------------------------------------------------
             if (associated(clm_pf_idata%ispec_decomp_c) .or. &
                 associated(clm_pf_idata%ispec_decomp_n)) then
               do k=1,clm_pf_idata%ndecomp_pools
                 if( cur_data_mediator%idof == (offsetim + clm_pf_idata%ispec_decomp_c(k)) ) then
                   call VecRestoreArrayReadF90(clm_pf_idata%rate_decomp_c_pfs, rate_pf_loc, ierr)
                   call where_checkerr(ierr, subname, __FILE__, __LINE__)
                   exit   ! exit the 'do k=1, ...' loop

                 elseif( cur_data_mediator%idof == (offsetim + clm_pf_idata%ispec_decomp_n(k)) ) then
                   call VecRestoreArrayReadF90(clm_pf_idata%rate_decomp_n_pfs, rate_pf_loc, ierr)
                   call where_checkerr(ierr, subname, __FILE__, __LINE__)
                   exit   ! exit the 'do k=1, ...' loop
                 endif
               enddo
             endif
         !--------------------------------------------------------------------------

           endif

           ! MUST call the following subroutine, OTHERWISE there is one time-step delay of data-passing
           call DataMediatorUpdate(realization%tran_data_mediator_list,  &
                            realization%field%tran_mass_transfer, &
                            realization%option)

          endif ! if(found_rtmasstr) block

        ! end of class is (data_mediator_dataset_type)

        class default
          option%io_buffer = ' tran_mass_DATASET is not of type of ' // &
            'data_mediator_dataset_type.'
          call printErrMsg(option)

      end select

      cur_data_mediator => cur_data_mediator%next
    end do

    call VecRestoreArrayReadF90(field%volume0,volume_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

  end subroutine pflotranModelSetBGCRatesFromCLM

  ! ************************************************************************** !
  !
  ! pflotranModelUpdateAqConcFromCLM:
  !  Get CLM aqueous nutrient concentrations (NH4, NO3 at this momment),
  !  converts from CLM units into PFLOTRAN units, and reset concentrations in PFLOTRAN
  !
  !   notes: when NOT coupled with PF Hydrology, forcing CLM water saturation to reset
  !          PF's global saturation status WOULD cause aq. phase element mass balance issue
  !
  !
  subroutine pflotranModelUpdateAqConcFromCLM(pflotran_model)

    use Global_Aux_module
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Discretization_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer          :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(field_type), pointer                 :: field
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid

    type(global_auxvar_type), pointer         :: global_auxvars(:)
    type(aq_species_type), pointer            :: cur_aq_spec

    PetscErrorCode     :: ierr
    PetscInt           :: local_id
    PetscInt           :: ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscReal, pointer   :: porosity0_loc_p(:)              ! current CLM-updated porosity

    PetscReal, pointer :: porosity_pre_pf_loc(:)            !previous time-step porosity (m3/m3 bulk soil)
    PetscReal, pointer :: soillsat_pre_pf_loc(:)            !previous time-step soil liq. water saturation (0 - 1)

    PetscInt :: offset, offset_aq

    character(len=MAXWORDLENGTH) :: word
    PetscReal :: porosity, saturation, theta ! for concentration conversion from mol/m3 to mol/L
    PetscReal :: theta_pre

    !-------------------------------------------------------------------------
    subname = 'pflotranModelUpdateAqConcFromCLM'
    !-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    patch => realization%patch
    grid  => patch%grid
    field => realization%field

    global_auxvars  => patch%aux%Global%auxvars

    !-----------------------------------------------------------------

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayReadF90(field%porosity0, porosity0_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! the previous time-step porosity and lsat from PFLOTRAN, saved in clm_pf_idata
    ! NOTE: make sure NOT modified by CLM
    call VecGetArrayReadF90(clm_pf_idata%effporosity_pfp, porosity_pre_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayReadF90(clm_pf_idata%soillsat_pfp, soillsat_pre_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    do local_id = 1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      if (ghosted_id<=0 .or. local_id <= 0) cycle ! bypass ghosted corner cells
      !if (patch%imat(local_id) <= 0) cycle  !(TODO) imat IS 0 for some cells when decomposing domain in X and Y directions.

      offset = (local_id - 1)*realization%reaction%ncomp
      offset_aq = offset + realization%reaction%offset_aqueous

      ! at this moment, both 'saturtion' and 'porosity' have been updated from CLM and pass to PF
      ! otherwise, it's no use.
      saturation = global_auxvars(ghosted_id)%sat(1)          ! using 'ghosted_id' if from 'global_auxvars'
      porosity = porosity0_loc_p(local_id)                    ! using 'local_id' if from 'field%???'
      theta = saturation * porosity

      ! adjusting aq. species conc. due to CLM pass-in theta (porosity X saturation)

      ! previous timestep saved PF's 'porosity' and 'saturation' in clm-pf-idata%???_pfp
      theta_pre = porosity_pre_pf_loc(local_id)*  &
                  soillsat_pre_pf_loc(local_id)

      cur_aq_spec => realization%reaction%primary_species_list
      do
        if (.not.associated(cur_aq_spec)) exit

        if (theta_pre> 1.d-20 .and. theta > 1.0d-20) then
          xx_p(offset_aq + cur_aq_spec%id) = xx_p(offset_aq + cur_aq_spec%id) &
                                            * theta_pre / theta
        end if
        cur_aq_spec => cur_aq_spec%next
      enddo

    enddo

    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(field%porosity0, porosity0_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecRestoreArrayReadF90(clm_pf_idata%effporosity_pfp, porosity_pre_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayReadF90(clm_pf_idata%soillsat_pfp, soillsat_pre_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  end subroutine pflotranModelUpdateAqConcFromCLM

  ! **************************************************************************************** !
  !> This routine reset gas concenctration (as immobile species) after emission adjusting
  !> in CLM-PFLOTRAN interface
  !! Note: this is a temporary work-around, because PF doesn't have BGC gas
  !!       processes at this moment
  !> @author
  !! F.-M. Yuan
  !!
  !! date: 3/19/2014
  !
  subroutine pflotranModelUpdateAqGasesFromCLM(pflotran_model)

    use Global_Aux_module
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Discretization_module
    use Reaction_Aux_module
    use Reaction_Immobile_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer          :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(field_type), pointer                 :: field
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid

    type(global_auxvar_type), pointer :: global_auxvars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: gco2_vr_pf_loc(:)              ! (molC/m3 bulk soil)
    PetscScalar, pointer :: gn2_vr_pf_loc(:)               ! (molN2-N/m3 bulk soil)
    PetscScalar, pointer :: gn2o_vr_pf_loc(:)              ! (molN2O-N/m3 bulk soil)

    PetscInt :: offset

    character(len=MAXWORDLENGTH) :: word

    !-------------------------------------------------------------------------
    subname = 'ModelUpdateAqGasesFromCLM'
    !-------------------------------------------------------------------------
    option=> pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    patch => realization%patch
    grid  => patch%grid
    field => realization%field

    global_auxvars  => patch%aux%Global%auxvars

    !-----------------------------------------------------------------

    ! mapping CLM vecs to PF vecs
    if(clm_pf_idata%ispec_co2 > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%gco2_vr_clmp, &
                                    clm_pf_idata%gco2_vr_pfs)
    endif

    if(clm_pf_idata%ispec_n2o > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%gn2o_vr_clmp, &
                                    clm_pf_idata%gn2o_vr_pfs)
    endif

    if(clm_pf_idata%ispec_n2 > 0) then
       call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%gn2_vr_clmp, &
                                    clm_pf_idata%gn2_vr_pfs)
    endif

    ! (iii) get the 'PF' vecs for resetting data
    call VecGetArrayF90(clm_pf_idata%gco2_vr_pfs, gco2_vr_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%gn2_vr_pfs, gn2_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%gn2o_vr_pfs, gn2o_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)  ! extract data from pflotran internal portion

    do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (ghosted_id <= 0 .or. local_id <= 0) cycle
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) < 0) cycle
        endif

        offset = (local_id-1) * realization%reaction%ncomp &
                + realization%reaction%offset_immobile

        if(clm_pf_idata%ispec_co2 > 0) then
            xx_p(offset + clm_pf_idata%ispec_co2) = max(gco2_vr_pf_loc(ghosted_id), xeps0_c)
        endif

        if(clm_pf_idata%ispec_n2 > 0) then
            xx_p(offset + clm_pf_idata%ispec_n2) = max(gn2_vr_pf_loc(ghosted_id), xeps0_n)
        endif

        if(clm_pf_idata%ispec_n2o > 0) then
            xx_p(offset + clm_pf_idata%ispec_n2o) = max(gn2o_vr_pf_loc(ghosted_id), xeps0_n)
        endif

    enddo

    call VecRestoreArrayF90(clm_pf_idata%gco2_vr_pfs, gco2_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%gn2_vr_pfs, gn2_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%gn2o_vr_pfs, gn2o_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  end subroutine pflotranModelUpdateAqGasesFromCLM

  ! **********************************************************************************!
  !
  ! This routine get updated bgc states/fluxes evoloved by PFLOTRAN.
  !
  !!
  subroutine pflotranModelGetBgcVariablesFromPF(pflotran_model)

    use Global_Aux_module
    use Material_Aux_class
    use Realization_Base_class
    use Patch_module
    use Grid_module
    use Option_module
    use Field_module
    use Reaction_module
    use Reaction_Aux_module
    use Reactive_Transport_module, only : RTUpdateAuxVars
    use Reactive_Transport_Aux_module
    use Reaction_Immobile_Aux_module
    use Discretization_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer          :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(field_type), pointer                 :: field
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid

    type(reaction_type), pointer              :: reaction
    type(global_auxvar_type), pointer         :: global_auxvar
    type(material_auxvar_type), pointer       :: material_auxvar
    type(reactive_transport_auxvar_type), pointer :: rt_auxvar

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: xx_p(:)

    PetscScalar, pointer :: decomp_cpools_vr_pf_loc(:)      ! (molesC/m3)
    PetscScalar, pointer :: decomp_npools_vr_pf_loc(:)      ! (molesN/m3)
    PetscScalar, pointer :: accextrnh4_vr_pf_loc(:)         ! (molesN/m3)
    PetscScalar, pointer :: accextrno3_vr_pf_loc(:)         ! (molesN/m3)
    PetscScalar, pointer :: smin_no3_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: smin_nh4sorb_vr_pf_loc(:)       ! (molesN/m3)
    PetscScalar, pointer :: gco2_vr_pf_loc(:)               ! (molC/m3)
    PetscScalar, pointer :: gn2_vr_pf_loc(:)                ! (molN2/m3)
    PetscScalar, pointer :: gn2o_vr_pf_loc(:)               ! (molN2O/m3)
    PetscScalar, pointer :: acchr_vr_pf_loc(:)              ! (molesC/m3)
    PetscScalar, pointer :: accnmin_vr_pf_loc(:)            ! (molesN/m3)
    PetscScalar, pointer :: accnimmp_vr_pf_loc(:)           ! (molesN/m3)
    PetscScalar, pointer :: accnimm_vr_pf_loc(:)            ! (molesN/m3)
    PetscScalar, pointer :: accngasmin_vr_pf_loc(:)         ! (molesN/m3)
    PetscScalar, pointer :: accngasnitr_vr_pf_loc(:)        ! (molesN/m3)
    PetscScalar, pointer :: accngasdeni_vr_pf_loc(:)        ! (molesN/m3)

    PetscReal, pointer :: porosity_loc_p(:)

    PetscReal :: porosity, saturation, theta ! for concentration conversion from mol/m3 to mol/L
    PetscReal :: den_kg_per_L, xmass

    PetscReal :: conc
    PetscInt  :: offset, offsetim

    Vec                  :: vec_pfp
    Vec                  :: vec_clms
    PetscScalar, pointer :: array_pfp(:), array_clms(:), array_temp(:)
    PetscInt             :: j, k, vec_offset

    PetscReal :: zeroing_conc = 1.0d-50

    !-------------------------------------------------------------------------
    subname = 'ModelSetBgcVariablesFromPF'
    !-------------------------------------------------------------------------
    option=> pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    patch => realization%patch
    grid  => patch%grid
    field => realization%field

    reaction => realization%reaction

    ! using user-input zeroing concentration, if available
    if (Initialized(reaction%truncated_concentration)) then
      zeroing_conc = reaction%truncated_concentration
    endif

    !-----------------------------------------------------------------

    !
    ! (i) open the original 'pf' vecs for updating
    !

    call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_pfp, decomp_cpools_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_pfp, decomp_npools_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecGetArrayF90(clm_pf_idata%smin_no3_vr_pfp, smin_no3_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%smin_nh4_vr_pfp, smin_nh4_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%smin_nh4sorb_vr_pfp, smin_nh4sorb_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecGetArrayF90(clm_pf_idata%gco2_vr_pfp, gco2_vr_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%gn2_vr_pfp, gn2_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%gn2o_vr_pfp, gn2o_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecGetArrayF90(clm_pf_idata%accextrnh4_vr_pfp, accextrnh4_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%accextrno3_vr_pfp, accextrno3_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    if (clm_pf_idata%ispec_hrimm>0) then
      call VecGetArrayF90(clm_pf_idata%acctothr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecGetArrayF90(clm_pf_idata%acchr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nmin>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else   ! can have both, but now NOT
      call VecGetArrayF90(clm_pf_idata%accnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimp>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else   ! can have both, but now NOT
      call VecGetArrayF90(clm_pf_idata%accnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimm>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else   ! can have both, but now NOT
      call VecGetArrayF90(clm_pf_idata%accnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    call VecGetArrayF90(clm_pf_idata%accngasmin_vr_pfp, accngasmin_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%accngasnitr_vr_pfp, accngasnitr_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%accngasdeni_vr_pfp, accngasdeni_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    ! (ii) pass the data from internal to CLM-PFLOTRAN interface vecs
    !

    call VecGetArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)  ! extract data from pflotran internal vecs
    call VecGetArrayReadF90(field%porosity0, porosity_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! loop over cells
    do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (ghosted_id <= 0 .or. local_id <= 0) cycle
        if (associated(patch%imat)) then
           if (patch%imat(ghosted_id) < 0) cycle
        endif

        global_auxvar    => patch%aux%Global%auxvars(ghosted_id)
        rt_auxvar        => patch%aux%RT%auxvars(ghosted_id)
        material_auxvar  => patch%aux%Material%auxvars(ghosted_id)

        ! for convertion btw liq. water mass and volume
        xmass = 1.d0
        if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(1)
        den_kg_per_L = global_auxvar%den_kg(1)*xmass*1.d-3      ! kg/L: kg/m3 *scaler* m3/L

        saturation = global_auxvar%sat(1)
        !porosity = porosity_loc_p(local_id)
        porosity = material_auxvar%porosity
        theta = saturation * porosity

        offset = (local_id - 1)*realization%reaction%ncomp

        offsetim = offset + realization%reaction%offset_immobile

        !-----------------------------------------------------------------
        ! loop over 'decomp_pools casecade'
        do k=1,clm_pf_idata%ndecomp_pools
          vec_offset = (k-1)*clm_pf_idata%nlpf_sub              ! MPI decomp_pfs vec: 'cell' first, then 'species'

          if(clm_pf_idata%ispec_decomp_c(k) > 0) then
             ! field%tran_xx vec IS arranged 'species' first and then 'cell'
             !conc = xx_p(offsetim + ispec_decomp_c(k))                      ! unit: M (mol/m3)
             conc = rt_auxvar%immobile(clm_pf_idata%ispec_decomp_c(k))
             ! MPI decomp_pfp vec: 'cell' first, then 'species'
             decomp_cpools_vr_pf_loc(vec_offset+local_id)= max(conc, 0.d0)

          endif

          if(clm_pf_idata%ispec_decomp_n(k) > 0) then
             ! field%tran_xx vec IS arranged 'species' first and then 'cell'
!             conc = xx_p(offsetim + ispec_decomp_n(k))                      ! unit: M (mol/m3)
             conc = rt_auxvar%immobile(clm_pf_idata%ispec_decomp_n(k))
             ! MPI decomp_pfp vec: 'cell' first, then 'species'
             decomp_npools_vr_pf_loc(vec_offset+local_id)= max(conc, 0.d0)
          endif

          if(clm_pf_idata%ispec_decomp_hr(k) > 0) then
             if(clm_pf_idata%ispec_hrimm <= 0) then   ! only when NO total HR tracked in PF, do the data-passing
               ! field%tran_xx vec IS arranged 'species' first and then 'cell'
!               conc = xx_p(offsetim + ispec_decomp_hr(k))                      ! unit: M (mol/m3)
               conc = rt_auxvar%immobile(clm_pf_idata%ispec_decomp_hr(k))
               ! MPI acchr_vr_pfp vec: 'cell' first, then 'species'
               acchr_vr_pf_loc(vec_offset+local_id)= max(conc-zeroing_conc, 0.d0)
             endif

             ! resetting the tracking variable state so that cumulative IS for the time-step
             xx_p(offsetim + clm_pf_idata%ispec_decomp_hr(k)) = zeroing_conc
          endif

          if(clm_pf_idata%ispec_decomp_nmin(k) > 0) then
            if (clm_pf_idata%ispec_nmin <= 0) then   ! only when NO total NMIN tracked in PF, then do data-passing
              ! field%tran_xx vec IS arranged 'species' first and then 'cell'
!             conc = xx_p(offsetim + ispec_decomp_nmin(k))                      ! unit: M (mol/m3)
              conc = rt_auxvar%immobile(clm_pf_idata%ispec_decomp_nmin(k))
              ! MPI accnmin_vr_pfp vec: 'cell' first, then 'species'
              accnmin_vr_pf_loc(vec_offset+local_id)= max(conc-zeroing_conc, 0.d0)
            endif

            ! resetting the tracking variable state so that cumulative IS for the time-step
            xx_p(offsetim + clm_pf_idata%ispec_decomp_nmin(k)) = zeroing_conc
          endif

          if(clm_pf_idata%ispec_decomp_nimp(k) > 0) then
            if (clm_pf_idata%ispec_nimp <= 0) then   ! only when NO total NIMP tracked in PF, then do data-passing
              ! field%tran_xx vec IS arranged 'species' first and then 'cell'
!             conc = xx_p(offsetim + ispec_decomp_nimp(k))                      ! unit: M (mol/m3)
              conc = rt_auxvar%immobile(clm_pf_idata%ispec_decomp_nimp(k))
              ! MPI accnmin_vr_pfp vec: 'cell' first, then 'species'
              accnimmp_vr_pf_loc(vec_offset+local_id)= max(conc-zeroing_conc, 0.d0)
            endif

            ! resetting the tracking variable state so that cumulative IS for the time-step
            xx_p(offsetim + clm_pf_idata%ispec_decomp_nimp(k)) = zeroing_conc
          endif

          if(clm_pf_idata%ispec_decomp_nimm(k) > 0) then
            if(clm_pf_idata%ispec_nimm <= 0) then   ! only when NO total NIMM tracked in PF, then do data-passing
              ! field%tran_xx vec IS arranged 'species' first and then 'cell'
!             conc = xx_p(offsetim + ispec_decomp_nimm(k))                      ! unit: M (mol/m3)
              conc = rt_auxvar%immobile(clm_pf_idata%ispec_decomp_nimm(k))
              ! MPI accnmin_vr_pfp vec: 'cell' first, then 'species'
              accnimm_vr_pf_loc(vec_offset+local_id)= max(conc-zeroing_conc, 0.d0)
            endif

            ! resetting the tracking variable state so that cumulative IS for the time-step
            xx_p(offsetim + clm_pf_idata%ispec_decomp_nimm(k)) = zeroing_conc
          endif

        enddo  ! end of loop of 'decomp_pools'

        ! total HR from all decomposition
        if (clm_pf_idata%ispec_hrimm > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_hrimm)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_hrimm)
           acchr_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + clm_pf_idata%ispec_hrimm) = zeroing_conc
        endif

        !
        if(clm_pf_idata%ispec_nh4 > 0) then
           !conc = xx_p(offset + ispec_nh4) * theta * 1000.0d0      ! 7-8-2015: this is NOT right.

            ! the following approach appears more like what output module does in pflotran
            ! but needs further checking if it efficient as directly read from 'xx_p' as above
            ! 7-8-2015: checked that the directly reading from 'xx_p' is incorrect, which may
            !           cause mass-balance errors. (TODO - for immobile species appears OK)
            ! 7-14-2015: if taking data from 'xx_p', the unit always is in moles/m3 bulk volume
            !            if taking data from 'rt_auxvar%total', for aq. species, is in moles/L water
            !        but, still has errors - best solution: fixed water density to 1.d3 kg/m3 in input cards
            !                                      which will eliminate the difference btw two approaches

            conc = rt_auxvar%total(clm_pf_idata%ispec_nh4, 1) * (theta * 1000.0d0 * den_kg_per_L)

            smin_nh4_vr_pf_loc(local_id) = max(conc, 0.d0)

            if (clm_pf_idata%ispec_nh4sorb>0) then    ! kinetic-languir adsorption reaction sandbox used for soil NH4+ absorption
                !conc = xx_p(offsetim + clm_pf_idata%ispec_nh4sorb)                 ! unit: M (mol/m3)
                conc = rt_auxvar%immobile(clm_pf_idata%ispec_nh4sorb)
                smin_nh4sorb_vr_pf_loc(local_id) = max(conc, 0.d0)
            elseif (reaction%neqsorb > 0) then  ! equilibrium-sorption reactions used
                conc = rt_auxvar%total_sorb_eq(clm_pf_idata%ispec_nh4)
                smin_nh4sorb_vr_pf_loc(local_id) = max(conc, 0.d0)
            endif

        elseif(clm_pf_idata%ispec_nh4s > 0) then ! as immobile form
           !conc = xx_p(offsetim + clm_pf_idata%ispec_nh4s)                   ! unit: M (molN/m3)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_nh4s)
           smin_nh4_vr_pf_loc(local_id)   = max(conc, 0.d0)

        endif

        if(clm_pf_idata%ispec_no3 > 0) then
           !conc = xx_p(offset + clm_pf_idata%ispec_no3) * theta * 1000.0d0                        ! 7-14-2015: converting from /m3 bulk.
           conc = rt_auxvar%total(clm_pf_idata%ispec_no3, 1) * (theta * 1000.0d0 * den_kg_per_L)     !
           smin_no3_vr_pf_loc(local_id)   = max(conc, 0.d0)

        elseif(clm_pf_idata%ispec_no3s > 0) then ! as immobile form
           !conc = xx_p(offsetim + clm_pf_idata%ispec_no3s)                   ! unit: M (molN/m3)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_no3s)
           smin_no3_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        ! immobile gas conc in mol/m3 bulk soil to aovid 'theta' inconsistence (due to porosity) during unit conversion
        if(clm_pf_idata%ispec_co2 > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_co2)                    ! unit: M (molC/m3)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_co2)
           gco2_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        if(clm_pf_idata%ispec_n2 > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_n2)                     ! unit: M (molN2/m3)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_n2)
           gn2_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        if(clm_pf_idata%ispec_n2o > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_n2o)                    ! unit: M (molN2O/m3)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_n2o)
           gn2o_vr_pf_loc(local_id)   = max(conc, 0.d0)
        endif

        ! tracking N bgc reaction fluxes

        ! total NMIN from decomposition
        if (clm_pf_idata%ispec_nmin > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_nmin)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_nmin)
           accnmin_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + clm_pf_idata%ispec_nmin) = zeroing_conc
        endif
        ! total NIMP from decomposition
        if (clm_pf_idata%ispec_nimp > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_nimp)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_nimp)
           accnimmp_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + clm_pf_idata%ispec_nimp) = zeroing_conc
        endif
        ! total NIMM from decomposition
        if (clm_pf_idata%ispec_nmin > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_nimm)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_nimm)
           accnimm_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + clm_pf_idata%ispec_nimm) = zeroing_conc
        endif

        if (clm_pf_idata%ispec_plantndemand > 0) then
           ! NO need to pass back to CLM

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + clm_pf_idata%ispec_plantndemand) = zeroing_conc
        endif

        if (clm_pf_idata%ispec_plantnh4uptake > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_plantnh4uptake)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_plantnh4uptake)
           accextrnh4_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + clm_pf_idata%ispec_plantnh4uptake) = zeroing_conc
        endif

        if (clm_pf_idata%ispec_plantno3uptake > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_plantno3uptake)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_plantno3uptake)
           accextrno3_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step only
           xx_p(offsetim + clm_pf_idata%ispec_plantno3uptake) = zeroing_conc
        endif

        if(clm_pf_idata%ispec_ngasmin > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_ngasmin)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_ngasmin)
           accngasmin_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + clm_pf_idata%ispec_ngasmin) = zeroing_conc
        endif

        if(clm_pf_idata%ispec_ngasnitr > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_ngasnitr)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_ngasnitr)
           accngasnitr_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + clm_pf_idata%ispec_ngasnitr) = zeroing_conc
        endif

        if(clm_pf_idata%ispec_ngasdeni > 0) then
           !conc = xx_p(offsetim + clm_pf_idata%ispec_ngasdeni)
           conc = rt_auxvar%immobile(clm_pf_idata%ispec_ngasdeni)
           accngasdeni_vr_pf_loc(local_id) = max(conc-zeroing_conc, 0.d0)

           ! resetting the tracking variable state so that cumulative IS for the time-step
           xx_p(offsetim + clm_pf_idata%ispec_ngasdeni) = zeroing_conc
        endif

    enddo

    !
    call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_pfp, decomp_cpools_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_pfp, decomp_npools_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecRestoreArrayF90(clm_pf_idata%smin_no3_vr_pfp, smin_no3_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4_vr_pfp, smin_nh4_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%smin_nh4sorb_vr_pfp, smin_nh4sorb_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecRestoreArrayF90(clm_pf_idata%gco2_vr_pfp, gco2_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%gn2_vr_pfp, gn2_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%gn2o_vr_pfp, gn2o_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecRestoreArrayF90(clm_pf_idata%accextrnh4_vr_pfp, accextrnh4_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%accextrno3_vr_pfp, accextrno3_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    if (clm_pf_idata%ispec_hrimm>0) then
      call VecGetArrayF90(clm_pf_idata%acctothr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%acchr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nmin>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%accnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimp>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%accnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimm>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%accnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    call VecRestoreArrayF90(clm_pf_idata%accngasmin_vr_pfp, accngasmin_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%accngasnitr_vr_pfp, accngasnitr_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%accngasdeni_vr_pfp, accngasdeni_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !
    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayReadF90(field%porosity0, porosity_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

     ! resetting the tracked variable states, which zeroed above
    call DiscretizationGlobalToLocal(realization%discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)
    call VecCopy(field%tran_xx,field%tran_yy,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

    !
    ! (iii) pass (mapping) the '_pfp' vecs to '_clms' vecs, which then can be passed to CLMCN
    ! (implemented in 'clm_pflotran_interfaceMod'
    !

    !-----------------------------------------------------------------
    ! for 'decomp_pools', NEED to do looping over each pool

    ! create temporary vecs/arrays for each 'decomp_pool' data-mapping
    call VecDuplicate(clm_pf_idata%zsoil_pfp, vec_pfp,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDuplicate(clm_pf_idata%zsoil_clms, vec_clms,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! decomp'C'
    if (associated(clm_pf_idata%ispec_decomp_c)) then
      ! assembly the 'vec_pfp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%decomp_cpools_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_cpools_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_pfp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%decomp_cpools_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_pfp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlpf_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%decomp_cpools_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec
        call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                        option,                &
                                        vec_pfp,                              &
                                        vec_clms)

        ! insert 'vec_clms' into the whole '_clms' vec
        vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%decomp_cpools_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngclm_sub) = array_clms

        call VecRestoreArrayF90(clm_pf_idata%decomp_cpools_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_clms' vec, after k compositing pools updated
      call VecAssemblyBegin(clm_pf_idata%decomp_cpools_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_cpools_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! decomp'N'
    if (associated(clm_pf_idata%ispec_decomp_n)) then
      ! assembly the 'vec_pfp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%decomp_npools_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_npools_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_pfp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%decomp_npools_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_pfp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlpf_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%decomp_npools_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec
        call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                        option,                &
                                        vec_pfp,                              &
                                        vec_clms)

        ! insert 'vec_clms' into the whole '_clms' vec
        vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%decomp_npools_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngclm_sub) = array_clms

        call VecRestoreArrayF90(clm_pf_idata%decomp_npools_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_clms' vec, after k compositing pools updated
      call VecAssemblyBegin(clm_pf_idata%decomp_npools_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%decomp_npools_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! HR from decomp'C'
    if (associated(clm_pf_idata%ispec_decomp_c) .and. clm_pf_idata%ispec_hrimm<=0) then
      ! assembly the 'vec_pfp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%acchr_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%acchr_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_pfp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%acchr_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_pfp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlpf_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%acchr_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec
        call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                        option,                &
                                        vec_pfp,                              &
                                        vec_clms)

        ! insert 'vec_clms' into the whole '_clms' vec
        vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%acchr_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngclm_sub) = array_clms

        call VecRestoreArrayF90(clm_pf_idata%acchr_vr_clms, array_temp, ierr)
         call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_clms' vec, after k compositing pools updated
      call VecAssemblyBegin(clm_pf_idata%acchr_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%acchr_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! NMIN from decomp'N'
    if (associated(clm_pf_idata%ispec_decomp_n) .and. clm_pf_idata%ispec_nmin<=0) then
      ! assembly the 'vec_pfp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%accnmin_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%accnmin_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_pfp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%accnmin_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_pfp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlpf_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%accnmin_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec
        call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                        option,                &
                                        vec_pfp,                              &
                                        vec_clms)

        ! insert 'vec_clms' into the whole '_clms' vec
        vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%accnmin_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngclm_sub) = array_clms

        call VecRestoreArrayF90(clm_pf_idata%accnmin_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_clms' vec, after k compositing pools updated
      call VecAssemblyBegin(clm_pf_idata%accnmin_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%accnmin_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! NIMP from decomp'N'
    if (associated(clm_pf_idata%ispec_decomp_n) .and. clm_pf_idata%ispec_nimp<=0) then
      ! assembly the 'vec_pfp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%accnimmp_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%accnimmp_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_pfp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%accnimmp_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_pfp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlpf_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%accnimmp_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec
        call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                        option,                &
                                        vec_pfp,                              &
                                        vec_clms)

        ! insert 'vec_clms' into the whole '_clms' vec
        vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%accnimmp_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngclm_sub) = array_clms

        call VecRestoreArrayF90(clm_pf_idata%accnimmp_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_clms' vec, after k compositing pools updated
      call VecAssemblyBegin(clm_pf_idata%accnimmp_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%accnimmp_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! NIMM from decomp'N'
    if (associated(clm_pf_idata%ispec_decomp_n) .and. clm_pf_idata%ispec_nimm<=0) then
      ! assembly the 'vec_pfp' (?? not sure if needed, though 'PETSC' manual said so)
      call VecAssemblyBegin(clm_pf_idata%accnimm_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%accnimm_vr_pfp, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      do k=1,clm_pf_idata%ndecomp_pools
        ! get a seg. of data from the whole '_pfp' vec for the 'k'th pool
        vec_offset = (k-1)*clm_pf_idata%nlpf_sub       ! MPI decomp_clmp vec: 'cell' first, then 'species'
        call VecGetArrayReadF90(clm_pf_idata%accnimm_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_pfp = array_temp(vec_offset+1:vec_offset+clm_pf_idata%nlpf_sub)

        call VecRestoreArrayReadF90(clm_pf_idata%accnimm_vr_pfp, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayF90(vec_pfp, array_pfp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        ! mapping from MPI vec to Seq. vec
        call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                        option,                &
                                        vec_pfp,                              &
                                        vec_clms)

        ! insert 'vec_clms' into the whole '_clms' vec
        vec_offset = (k-1)*clm_pf_idata%ngclm_sub       ! SEQ. decomp_pfs vec: 'cell' first, then 'species'
        call VecGetArrayF90(clm_pf_idata%accnimm_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecGetArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

        array_temp(vec_offset+1:vec_offset+clm_pf_idata%ngclm_sub) = array_clms

        call VecRestoreArrayF90(clm_pf_idata%accnimm_vr_clms, array_temp, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)
        call VecRestoreArrayReadF90(vec_clms, array_clms, ierr)
        call where_checkerr(ierr, subname, __FILE__, __LINE__)

      enddo

      ! assembly the whole '_clms' vec, after k compositing pools updated
      call VecAssemblyBegin(clm_pf_idata%accnimm_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
      call VecAssemblyEnd(clm_pf_idata%accnimm_vr_clms, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)

    endif

    ! clear-up of temporary vecs/arrarys
    call VecDestroy(vec_pfp,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDestroy(vec_clms,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !-----------------------------------------------------------------
    ! for single element, mapping approach directly

    !
    if(clm_pf_idata%ispec_no3 > 0 .or. clm_pf_idata%ispec_no3s > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%smin_no3_vr_pfp, &
                                    clm_pf_idata%smin_no3_vr_clms)
    endif

    if(clm_pf_idata%ispec_nh4 > 0 .or. clm_pf_idata%ispec_nh4s > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%smin_nh4_vr_pfp, &
                                    clm_pf_idata%smin_nh4_vr_clms)

      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%smin_nh4sorb_vr_pfp, &
                                    clm_pf_idata%smin_nh4sorb_vr_clms)
    endif

    !
    if (clm_pf_idata%ispec_co2 > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%gco2_vr_pfp, &
                                    clm_pf_idata%gco2_vr_clms)
    endif

    !
    if(clm_pf_idata%ispec_plantnh4uptake > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%accextrnh4_vr_pfp, &
                                    clm_pf_idata%accextrnh4_vr_clms)
    endif

    if(clm_pf_idata%ispec_plantno3uptake > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%accextrno3_vr_pfp, &
                                    clm_pf_idata%accextrno3_vr_clms)
    endif

    if(clm_pf_idata%ispec_ngasmin > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%accngasmin_vr_pfp, &
                                    clm_pf_idata%accngasmin_vr_clms)
    endif

    if(clm_pf_idata%ispec_ngasnitr > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%accngasnitr_vr_pfp, &
                                    clm_pf_idata%accngasnitr_vr_clms)
    endif

    if(clm_pf_idata%ispec_ngasdeni > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%accngasdeni_vr_pfp, &
                                    clm_pf_idata%accngasdeni_vr_clms)
    endif

    ! if total HR/NMIN/NIMM/NIMP tracked, otherwise have already done above
    if(clm_pf_idata%ispec_hrimm > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%acctothr_vr_pfp, &
                                    clm_pf_idata%acctothr_vr_clms)
    endif

    if(clm_pf_idata%ispec_nmin > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%acctotnmin_vr_pfp, &
                                    clm_pf_idata%acctotnmin_vr_clms)
    endif

    if(clm_pf_idata%ispec_nimm > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%acctotnimm_vr_pfp, &
                                    clm_pf_idata%acctotnimm_vr_clms)
    endif

    if(clm_pf_idata%ispec_nimp > 0) then
      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    option, &
                                    clm_pf_idata%acctotnimmp_vr_pfp, &
                                    clm_pf_idata%acctotnimmp_vr_clms)
    endif
    !
    ! -----------------------------------------------------------------------------------------
    !
    ! (iv) after data-passing of "tracking variables", may need to ZERO them for next clm timestep
    !     (may not be needed due to update at end of every timestep, but just in case)
    !

    if (clm_pf_idata%ispec_hrimm>0) then
      call VecGetArrayF90(clm_pf_idata%acctothr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecGetArrayF90(clm_pf_idata%acchr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nmin>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecGetArrayF90(clm_pf_idata%accnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimp>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecGetArrayF90(clm_pf_idata%accnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimm>0) then
      call VecGetArrayF90(clm_pf_idata%acctotnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecGetArrayF90(clm_pf_idata%accnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    call VecGetArrayF90(clm_pf_idata%accextrnh4_vr_pfp, accextrnh4_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%accextrno3_vr_pfp, accextrno3_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%accngasmin_vr_pfp, accngasmin_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%accngasnitr_vr_pfp, accngasnitr_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%accngasdeni_vr_pfp, accngasdeni_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    acchr_vr_pf_loc(:)      = 0.d0
    accnmin_vr_pf_loc(:)    = 0.d0
    accnimmp_vr_pf_loc(:)   = 0.d0
    accnimm_vr_pf_loc(:)    = 0.d0
    accextrnh4_vr_pf_loc(:) = 0.d0
    accextrno3_vr_pf_loc(:) = 0.d0
    accngasmin_vr_pf_loc(:) = 0.d0
    accngasnitr_vr_pf_loc(:)= 0.d0
    accngasdeni_vr_pf_loc(:)= 0.d0

    if (clm_pf_idata%ispec_hrimm>0) then
      call VecRestoreArrayF90(clm_pf_idata%acctothr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%acchr_vr_pfp, acchr_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nmin>0) then
      call VecRestoreArrayF90(clm_pf_idata%acctotnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%accnmin_vr_pfp, accnmin_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimp>0) then
      call VecRestoreArrayF90(clm_pf_idata%acctotnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%accnimmp_vr_pfp, accnimmp_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    if (clm_pf_idata%ispec_nimm>0) then
      call VecRestoreArrayF90(clm_pf_idata%acctotnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    else
      call VecRestoreArrayF90(clm_pf_idata%accnimm_vr_pfp, accnimm_vr_pf_loc, ierr)
      call where_checkerr(ierr, subname, __FILE__, __LINE__)
    endif

    call VecRestoreArrayF90(clm_pf_idata%accextrnh4_vr_pfp, accextrnh4_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%accextrno3_vr_pfp, accextrno3_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%accngasmin_vr_pfp, accngasmin_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%accngasnitr_vr_pfp, accngasnitr_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%accngasdeni_vr_pfp, accngasdeni_vr_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    !

  end subroutine pflotranModelGetBgcVariablesFromPF
! ************************************************************************** !

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!   (BLANK AS INTENDED)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
! THE FOLLOWING BLOCKS OF CODES ARE FOR CLM-PFLOTRAN Thermal-Hydrology (TH) COUPLING
! (TODO)
!
  ! ************************************************************************** !
  !
  ! pflotranModelSetInternalTHStatesfromCLM: Set initial TH States from CLM
  !
  ! Note: This subroutine directly set initial soil temperature and saturation from CLM
  !       It's needed because of uniform initialization of TH states in PFLOTRAN, which
  !       are from the input card.
  ! (This is different from the 'pflotranModelUpdateTHfromCLM', which pass TH from CLM to
  !   pflotran's global variables and will not affect the internal vec of TH mode).

  ! author: Fengming YUAN
  ! date: 9/23/2013
  ! ************************************************************************** !
subroutine pflotranModelSetInternalTHStatesfromCLM(pflotran_model, PRESSURE_DATAPASSING)

    use Option_module
    use Patch_module
    use Grid_module
    use Field_module
    use Discretization_module
    use TH_Aux_module
    use Material_Aux_class
    use Global_Aux_module

    use Realization_Base_class
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use Characteristic_Curves_module
    use Characteristic_Curves_Base_module
    use Characteristic_Curves_Common_module
    use TH_module, only : THUpdateAuxVars

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    logical, intent(in)                       :: PRESSURE_DATAPASSING
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(global_auxvar_type), pointer         :: global_auxvars(:)
    class(material_auxvar_type), pointer      :: material_auxvars(:)

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    class(characteristic_curves_type), pointer  :: characteristic_curves

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id, istart, iend
    PetscInt           :: cur_sat_func_id
    PetscReal          :: liquid_saturation, capillary_pressure, dx, porosity
    PetscReal          :: liq_kgm3
    PetscReal, pointer :: xx_loc_p(:)

    PetscScalar, pointer :: soilt_pf_loc(:)      ! temperature [oC]
    PetscScalar, pointer :: soilpress_pf_loc(:)  ! water pressure (Pa)
    PetscScalar, pointer :: soilliq_pf_loc(:)    ! liq. water mass (kg/m3)
    PetscScalar, pointer :: soilice_pf_loc(:)    ! ice water mass (kg/m3)

    subname = 'ModelSetInternalTHStatesFromCLM'

!-------------------------------------------------------------------------
    option => pflotran_model%option

    select case(option%iflowmode)
      case (TH_MODE)
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soillsat_clmp, &
                                    clm_pf_idata%soillsat_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilisat_clmp, &
                                    clm_pf_idata%soilisat_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilt_clmp, &
                                    clm_pf_idata%soilt_pfs)
        !
        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilliq_clmp, &
                                    clm_pf_idata%soilliq_pfs)

        call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%soilice_clmp, &
                                    clm_pf_idata%soilice_pfs)
      case default
        if(option%ntrandof.le.0) then
            option%io_buffer='pflotranModelSetInitialTHStatesfromCLM ' // &
              'not implmented for this mode.'
            call printErrMsg(option)
        else
            ! reactive-transport without flow-mode on
            return
        endif
    end select

    !
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    global_auxvars  => patch%aux%Global%auxvars
    material_auxvars=> patch%aux%Material%auxvars


    call VecGetArrayF90(field%flow_xx, xx_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%soilliq_pfs, soilliq_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%soilice_pfs, soilice_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (ghosted_id <= 0 .or. local_id <= 0) cycle
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) < 0) cycle
       endif

       iend = local_id*option%nflowdof
       istart = iend-option%nflowdof+1

       porosity = material_auxvars(ghosted_id)%porosity
       liq_kgm3 = global_auxvars(ghosted_id)%den_kg(1) ! water den = kg/m^3

       ! soil hydraulic properties ID for current cell
       cur_sat_func_id = patch%sat_func_id(ghosted_id)
       characteristic_curves => patch% &
         characteristic_curves_array(cur_sat_func_id)%ptr

       if (PRESSURE_DATAPASSING) then
         xx_loc_p(istart)  = soilpress_pf_loc(ghosted_id)

         ! may need to recalculate 'saturation' from pressure
         capillary_pressure = option%reference_pressure - xx_loc_p(istart)
         select type(sf => characteristic_curves%saturation_function)
           !class is(sat_func_VG_type)
             ! not-yet (TODO)
           class is(sat_func_BC_type)
             call sf%Saturation(capillary_pressure, liquid_saturation, dx, option)

           class default
             option%io_buffer = 'Currently ONLY support Brooks_COREY saturation function type' // &
               ' when coupled with CLM.'
             call printErrMsg(option)
         end select

       else
         ! need to recalculate 'pressure' from saturation/water-mass
         liquid_saturation = soilliq_pf_loc(ghosted_id)/liq_kgm3/porosity
         select type(sf => characteristic_curves%saturation_function)
           !class is(sat_func_VG_type)
             ! not-yet (TODO)
           class is(sat_func_BC_type)
             call sf%CapillaryPressure(liquid_saturation, capillary_pressure, dx, option)

             xx_loc_p(istart) = option%reference_pressure - capillary_pressure

           class default
             option%io_buffer = 'Currently ONLY support Brooks_COREY saturation function type' // &
               ' when coupled with CLM.'
             call printErrMsg(option)
         end select

       end if

       if (option%iflowmode .eq. TH_MODE)  then
         xx_loc_p(istart+1)= soilt_pf_loc(ghosted_id)
       end if
    enddo

    call VecRestoreArrayF90(field%flow_xx, xx_loc_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%soilt_pfs, soilt_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%press_pfs, soilpress_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%soilliq_pfs, soilliq_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%soilice_pfs, soilice_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call DiscretizationGlobalToLocal(realization%discretization, field%flow_xx, &
         field%flow_xx_loc, NFLOWDOF)
    call VecCopy(field%flow_xx, field%flow_yy, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    select case(option%iflowmode)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
      case default
        if(option%ntrandof.le.0) then
           option%io_buffer='pflotranModelSetInitialTHStatesfromCLM ' // &
                 'not implmented for this mode.'
           call printErrMsg(option)
        endif
    end select

end subroutine pflotranModelSetInternalTHStatesfromCLM

  ! ************************************************************************** !
  ! pflotranModelSetSoilHbcs()
  ! refresh Hydrological BC variables from CLM to PF
  !
  ! by 1-18-2013: only water pressure-head type (dirichlet) available
  ! by 4-11-2013: dirichlet/neumman both available
  ! ************************************************************************** !
  subroutine pflotranModelSetSoilHbcsFromCLM(pflotran_model)

    use Realization_Base_class
    use Option_module
    use Patch_module
    use Grid_module
    use Coupler_module
    use Connection_module

    use TH_Aux_module

    use String_module

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use TH_module, only : THUpdateAuxVars

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    type(coupler_type), pointer :: boundary_condition
    type(connection_set_type), pointer :: cur_connection_set
    PetscInt :: ghosted_id, local_id, press_dof, iconn, sum_connection

    PetscBool:: HAVE_QFLUX_TOPBC, HAVE_PRESS_TOPBC, HAVE_EXFIL_TOPBC

    PetscErrorCode     :: ierr

    PetscScalar, pointer :: press_maxponding_pf_loc(:)  ! subsurface top boundary max. ponding pressure (Pa) (seepage BC)
    PetscScalar, pointer :: press_subsurf_pf_loc(:)     ! subsurface top boundary pressure-head (Pa) (dirichlet BC)
    PetscScalar, pointer :: qfluxw_subsurf_pf_loc(:)    ! subsurface top boundary infiltration rate (m/s) (neumann BC)
    PetscScalar, pointer :: qfluxv_subsurf_pf_loc(:)    ! subsurface top boundary evaporation rate (m/s) (neumann BC)
    PetscScalar, pointer :: press_subbase_pf_loc(:)     ! bottom boundary pressure-head (Pa) (dirichlet BC)
    PetscScalar, pointer :: qfluxw_subbase_pf_loc(:)     ! botoom boundary drainage flow rate (m/s) (neumann BC)

    PetscScalar, pointer :: toparea_p(:)                ! subsurface top area saved

    !------------------------------------------------------------------------------------

    subname = 'pflotranModelSetSoilHbcsFromCLM'

!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    patch           => realization%patch
    grid            => patch%grid

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%press_subsurf_clmp, &
                                    clm_pf_idata%press_subsurf_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%qfluxw_subsurf_clmp, &
                                    clm_pf_idata%qfluxw_subsurf_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%qfluxev_subsurf_clmp, &
                                    clm_pf_idata%qfluxev_subsurf_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%press_maxponding_clmp, &
                                    clm_pf_idata%press_maxponding_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dbot_to_pf_2dbot, &
                                    option, &
                                    clm_pf_idata%press_subbase_clmp, &
                                    clm_pf_idata%press_subbase_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dbot_to_pf_2dbot, &
                                    option, &
                                    clm_pf_idata%qfluxw_subbase_clmp, &
                                    clm_pf_idata%qfluxw_subbase_pfs)

    ! interface vecs of PF
    call VecGetArrayF90(clm_pf_idata%press_subsurf_pfs,  press_subsurf_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%qfluxw_subsurf_pfs,  qfluxw_subsurf_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%qfluxev_subsurf_pfs,  qfluxv_subsurf_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%press_subbase_pfs,  press_subbase_pf_loc,  ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%qfluxw_subbase_pfs, qfluxw_subbase_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%press_maxponding_pfs, press_maxponding_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(clm_pf_idata%area_top_face_pfp, toparea_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! passing from interface to internal
    select case(option%iflowmode)
      case (TH_MODE)
        press_dof = TH_PRESSURE_DOF
      case default
        option%io_buffer='pflotranModelSetTHbcs ' // &
          'not implmented for this mode.'
        call printErrMsg(option)
    end select

    ! need to check the BC list first, so that we have necessary BCs in a consistent way
    HAVE_QFLUX_TOPBC = PETSC_FALSE  ! top BC: (water) flux type (NEUMANN)
    HAVE_PRESS_TOPBC = PETSC_FALSE  ! top BC: (water) pressure-head type (DIRICHLET)
    HAVE_EXFIL_TOPBC = PETSC_FALSE  ! top BC: (water) pressure-type one-way (SEEPAGE: hydrostatic, but one-way) for upward/outlet water flux

    boundary_condition => patch%boundary_condition_list%first
    do
      if (.not.associated(boundary_condition)) exit

      if(StringCompare(boundary_condition%name,'clm_gwflux_bc')) then
        if (boundary_condition%flow_condition%itype(press_dof) == NEUMANN_BC) then
          HAVE_QFLUX_TOPBC = PETSC_TRUE
        else
          option%io_buffer='pflotranModelSetTHbcs -  ' // &
              ' for CLM-PFLOTRAN coupling - flow condition MUST be named as following: ' // &
              ' "clm_gwflux_bc/NEUMANN " for subsurface-top TYPE I  '
          call printErrMsg(option)
        endif

      elseif(StringCompare(boundary_condition%name,'clm_gpress_bc')) then
        if (boundary_condition%flow_condition%itype(press_dof) == DIRICHLET_BC) then
          HAVE_PRESS_TOPBC = PETSC_TRUE
        else
          option%io_buffer='pflotranModelSetTHbcs -  ' // &
               ' for CLM-PFLOTRAN coupling - flow condition MUST be named as following: ' // &
               ' "clm_gpress_bc/DIRICHLET " for subsurface-top TYPE II  '
          call printErrMsg(option)
        endif

      elseif(StringCompare(boundary_condition%name,'exfiltration')) then
        if (boundary_condition%flow_condition%itype(press_dof) == SEEPAGE_BC) then
          HAVE_EXFIL_TOPBC = PETSC_TRUE
        else
          option%io_buffer='pflotranModelSetTHbcs -  ' // &
               ' for CLM-PFLOTRAN coupling - flow condition MUST be named as following: ' // &
               ' "exfiltration/SEEPAGE " for subsurface-top TYPE II - upward/outlet flow  '
          call printErrMsg(option)
        endif

      endif
      boundary_condition => boundary_condition%next
    end do
    if(.not.HAVE_QFLUX_TOPBC) then
      option%io_buffer='pflotranModelSetTHbcs -  ' // &
               ' for CLM-PFLOTRAN coupling - BC flow conditions DO NOT have : ' // &
               ' "clm_gwflux_bc/NEUMANN " for subsurface-top TYPE I  '
      !call printMsg(option)
    endif
    if(HAVE_PRESS_TOPBC .and. HAVE_EXFIL_TOPBC) then
      option%io_buffer='pflotranModelSetTHbcs -  ' // &
               ' for CLM-PFLOTRAN coupling - BC flow conditions are having both : ' // &
               ' "exfiltration/SEEPAGE " and "clm_gpress_bc/DIRICHLET" for subsurface-top '
      !call printMsg(option)
    endif

    ! assign data to BCs from CLM
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1

        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (ghosted_id <= 0 .or. local_id <= 0) cycle
        if (patch%imat(ghosted_id) < 0) cycle

        if(StringCompare(boundary_condition%name,'clm_gwflux_bc')) then
          if (boundary_condition%flow_condition%itype(press_dof) == NEUMANN_BC) then
            boundary_condition%flow_aux_real_var(press_dof,iconn)= &
              qfluxw_subsurf_pf_loc(iconn)

            if (HAVE_PRESS_TOPBC .and. HAVE_QFLUX_TOPBC) then
              cur_connection_set%area(iconn) = toparea_p(local_id)     ! normally it's ON (MPI vec, it's from 'local_id')
              if(press_subsurf_pf_loc(iconn) > clm_pf_idata%pressure_reference) then         ! shut-off the BC by resetting the BC 'area' to a tiny value
                cur_connection_set%area(iconn) = 0.d0
              endif
            endif

          endif
        endif

        if(StringCompare(boundary_condition%name,'clm_gpress_bc')) then
          if (boundary_condition%flow_condition%itype(press_dof) == DIRICHLET_BC) then
            boundary_condition%flow_aux_real_var(press_dof,iconn)= &
              press_subsurf_pf_loc(iconn)

            if (HAVE_PRESS_TOPBC .and. HAVE_QFLUX_TOPBC) then
              cur_connection_set%area(iconn) = 0.d0               ! normally shut-off this BC
              if(press_subsurf_pf_loc(iconn) > clm_pf_idata%pressure_reference) then         ! turn on the BC by resetting the BC 'area' to real value
                cur_connection_set%area(iconn) = toparea_p(local_id)
              endif
            endif

          endif

        endif

        ! for soil evaporation
        if(StringCompare(boundary_condition%name,'clm_gevflux_bc')) then
          if (boundary_condition%flow_condition%itype(press_dof) == NEUMANN_BC) then
            boundary_condition%flow_aux_real_var(press_dof,iconn)= &
              qfluxv_subsurf_pf_loc(iconn)
          endif
        endif

        ! for exfiltration
        if(StringCompare(boundary_condition%name,'clm_exfiltration_bc')) then
          if (boundary_condition%flow_condition%itype(press_dof) == SEEPAGE_BC) then
            boundary_condition%flow_aux_real_var(press_dof,iconn) = &
              max(option%reference_pressure, press_maxponding_pf_loc(iconn))
          endif
        endif

        ! bottom water flux
        if(StringCompare(boundary_condition%name,'clm_bflux_bc')) then
          if (boundary_condition%flow_condition%itype(press_dof) == DIRICHLET_BC) then
            boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       press_subbase_pf_loc(iconn)
          else if (boundary_condition%flow_condition%itype(press_dof) == NEUMANN_BC) then
            boundary_condition%flow_aux_real_var(press_dof,iconn)= &
                       qfluxw_subbase_pf_loc(iconn)

          else if (boundary_condition%flow_condition%itype(press_dof) /= ZERO_GRADIENT_BC) then
            option%io_buffer='pflotranModelSetTHbcs -  ' // &
                  ' for CLM-PFLOTRAN coupling - flow condition MUST be named as following: ' // &
                  ' "clm_bflux_bc/NEUMANN or ZERO_GRADIENT" for subsurface-base TYPE I;  ' // &
                  ' "clm_bflux_bc/DIRICHLET " for subsurface-base TYPE II;  '
            call printErrMsg(option)

          end if
        endif

      enddo

      boundary_condition => boundary_condition%next

    enddo

    call VecRestoreArrayF90(clm_pf_idata%press_subsurf_pfs, press_subsurf_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%qfluxw_subsurf_pfs, qfluxw_subsurf_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%qfluxev_subsurf_pfs, qfluxv_subsurf_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%press_subbase_pfs, press_subbase_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%qfluxw_subbase_pfs, qfluxw_subbase_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%press_maxponding_pfs, press_maxponding_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecRestoreArrayF90(clm_pf_idata%area_top_face_pfp, toparea_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    select case(option%iflowmode)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
      case default
        option%io_buffer='pflotranModelSetTHbcs ' // &
          'not implmented for this mode.'
        call printErrMsg(option)
    end select

  end subroutine pflotranModelSetSoilHbcsFromCLM

  ! ************************************************************************** !

  subroutine pflotranModelUpdateHSourceSink(pflotran_model)
  !
  ! Update the source/sink term of hydrology
  !
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  ! Revised by Fengming YUAN

    use Connection_module
    use Coupler_module
    use Patch_module
    use Grid_module
    use Material_Aux_class
    use Option_module
    use String_module

    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(coupler_type), pointer               :: source_sink

    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    class(material_auxvar_type), pointer      :: material_auxvars(:)

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    type(connection_set_type), pointer        :: cur_connection_set

    PetscScalar, pointer                      :: qflx_pf_loc(:), qflxt_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn, local_id, ghosted_id, sum_connection
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof, temperature_dof

    subname = 'pflotranModelUpdateHSourceSink'
!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    patch            => realization%patch
    grid             => patch%grid
    material_auxvars => patch%aux%Material%auxvars

!-------------------------------------------------------------------------

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%qflow_clmp, &
                                    clm_pf_idata%qflow_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%qflowt_clmp, &
                                    clm_pf_idata%qflowt_pfs)

    ! Find value of pressure-dof depending on flow mode
    select case (option%iflowmode)
      case (TH_MODE)
        press_dof       = TH_PRESSURE_DOF
        temperature_dof = TH_TEMPERATURE_DOF
      case default
        option%io_buffer = 'Unsupported Flow mode'
        call printErrMsg(option)
    end select

    ! Update the 'clm_et_ss' source/sink term
    call VecGetArrayF90(clm_pf_idata%qflow_pfs,qflx_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%qflowt_pfs,qflxt_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    found = PETSC_FALSE

    source_sink => realization%patch%source_sink_list%first

    sum_connection = 0
    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate Source/Sink from the list of Source/Sinks
      if(StringCompare(source_sink%name,'clm_et_ss')) then

        found = PETSC_TRUE
        if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS) then
          call printErrMsg(option,'clm_et_ss is not of ' // &
                           'HET_MASS_RATE_SS for water flow (RATE) ')
        endif

        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          source_sink%flow_aux_real_var(press_dof,iconn) = qflx_pf_loc(ghosted_id) &
                               *material_auxvars(ghosted_id)%volume                    ! kg/m3/sec --> kg/sec

          if(option%iflowmode == TH_MODE) then

            if(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF) == DIRICHLET_BC) then
              source_sink%flow_condition%temperature%dataset%rarray(1) = qflxt_pf_loc(ghosted_id)

            elseif(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF) == HET_DIRICHLET_BC) then
              source_sink%flow_aux_real_var(TWO_INTEGER,iconn) = qflxt_pf_loc(ghosted_id)

            elseif(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF) /= ZERO_GRADIENT_BC) then
              call printErrMsg(option,'clm_et_ss is not of ' // &
                           'DIRCHLET_BC or HET_DIRICHLET_BC or ZERO_GRADIENT_BC for temperature')
            endif

          endif

#ifdef CLM_PF_DEBUG
      ! the following checking shows data passing IS from 'ghosted_id' to 'iconn (local_id)' (multiple processors)
      write(option%myrank+200,*) 'checking H-et ss. -pf_model-UpdateSrcSink:', &
        'rank=',option%myrank, 'local_id=',local_id, 'ghosted_id=',ghosted_id, &
        'iconn=',iconn, 'qflx_pfs_loc(iconn)=',qflx_pf_loc(iconn), &
        'qflx_pfs_loc(ghosted_id)=',qflx_pf_loc(ghosted_id)
#endif

        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(clm_pf_idata%qflow_pfs,qflx_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%qflowt_pfs,qflxt_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !if(.not.found) &
      !call printMsg(option,'clm_et_ss not found in ' // &
      !                 'source-sink list of subsurface model.')

  end subroutine pflotranModelUpdateHSourceSink


! ************************************************************************** !

  subroutine pflotranModelUpdateSubsurfTCond(pflotran_model)
  !
  ! This routine updates subsurface boundary condtions of PFLOTRAN related to
  ! energy equation.
  !
  ! Author: Fengming YUAN, CCSI/ESD-ORNL
  ! Date: 08/10/2016
  !

    use Patch_module
    use Grid_module
    use Material_Aux_class
    use Option_module
    use String_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type
    use Connection_module
    use Coupler_module

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer        :: pflotran_model
    type(coupler_type), pointer               :: boundary_condition
    type(coupler_type), pointer               :: source_sink
    type(connection_set_type), pointer        :: cur_connection_set

    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    class(material_auxvar_type), pointer      :: material_auxvars(:)

    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization

    PetscScalar, pointer                      :: geflx_subsurf_pf_loc(:)   ! unit: MJ/m2/sec (all forms)
    PetscScalar, pointer                      :: geflxr_subsurf_pf_loc(:)  ! unit: MJ/m2/sec (net radiation)
    PetscScalar, pointer                      :: geflxl_subsurf_pf_loc(:)  ! unit: MJ/m2/sec (soil evaporation LE)
    PetscScalar, pointer                      :: gtemp_subsurf_pf_loc(:)
    PetscScalar, pointer                      :: geflx_subbase_pf_loc(:)   ! unit: MJ/m2/sec (geo-thermal flux)
    PetscScalar, pointer                      :: gtemp_subbase_pf_loc(:)

    PetscScalar, pointer                      :: geflow_sub_pf_loc(:)      ! unit: MJ/m3/sec (energy flow rate)

    PetscBool                                 :: HAVE_GTEMP_TOPBC, HAVE_GEV_TOPBC

    PetscInt                                  :: iconn, sum_connection
    PetscInt                                  :: local_id, ghosted_id
    PetscErrorCode                            :: ierr

    subname = 'pflotranModelUpdateSubsurfTCond'
!-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    patch            => realization%patch
    grid             => patch%grid
    material_auxvars => patch%aux%Material%auxvars

!-------------------------------------------------------------------------

    if (clm_pf_idata%nlpf_2dtop <= 0 .and. clm_pf_idata%ngpf_2dtop <= 0 ) return

    ! Map ground-heat flux from CLM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%eflux_subsurf_clmp, &
                                    clm_pf_idata%eflux_subsurf_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%efluxr_subsurf_clmp, &
                                    clm_pf_idata%efluxr_subsurf_pfs)

    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%efluxl_subsurf_clmp, &
                                    clm_pf_idata%efluxl_subsurf_pfs)

    ! Map ground-soil interface temperature from CLM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_clm_2dtop_to_pf_2dtop, &
                                    option, &
                                    clm_pf_idata%gtemp_subsurf_clmp, &
                                    clm_pf_idata%gtemp_subsurf_pfs)

    ! Map base-heat flux from CLM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_clm_2dbot_to_pf_2dbot, &
                                    option, &
                                    clm_pf_idata%eflux_subbase_clmp, &
                                    clm_pf_idata%eflux_subbase_pfs)

    ! Map base temperature from CLM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_clm_2dbot_to_pf_2dbot, &
                                    option, &
                                    clm_pf_idata%gtemp_subbase_clmp, &
                                    clm_pf_idata%gtemp_subbase_pfs)

    ! Map energy flow rate from CLM--to--PF cell
    call MappingSourceToDestination(pflotran_model%map_clm_sub_to_pf_sub, &
                                    option, &
                                    clm_pf_idata%eflow_clmp, &
                                    clm_pf_idata%eflow_pfs)

    ! Update the heat flux/ground-temperature BC term
    call VecGetArrayF90(clm_pf_idata%eflux_subsurf_pfs,geflx_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%efluxr_subsurf_pfs,geflxr_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%efluxl_subsurf_pfs,geflxl_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%gtemp_subsurf_pfs,gtemp_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! Update the 'clm_bflux_bc' base heat/energy flux/base-temperature BC term
    call VecGetArrayF90(clm_pf_idata%eflux_subbase_pfs,geflx_subbase_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%gtemp_subbase_pfs,gtemp_subbase_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! energy flow rate (source/sink)
    call VecGetArrayF90(clm_pf_idata%eflow_pfs,geflow_sub_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! need to check the BC list first, so that we have necessary BCs in a consistent way
    HAVE_GEV_TOPBC   = PETSC_FALSE  ! top BC: (water and heat) having soil evaporation flux type (NEUMANN)
    HAVE_GTEMP_TOPBC = PETSC_FALSE  ! top BC: (heat) having thermal state (temperature) type (DIRICHLET)

    boundary_condition => patch%boundary_condition_list%first
    do
      if (.not.associated(boundary_condition)) exit

      if(StringCompare(boundary_condition%name,'clm_gevflux_bc')) then
        if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) == NEUMANN_BC) then
          HAVE_GEV_TOPBC = PETSC_TRUE
        else
          option%io_buffer='pflotranModelSetTcond -  ' // &
               ' for CLM-PFLOTRAN coupling - flow condition MUST be named as following: ' // &
               ' "clm_gevflux_bc: ENERGY_FLUX neumann " for subsurface-top TYPE I  '
          call printErrMsg(option)
        endif

      elseif(StringCompare(boundary_condition%name,'clm_gtemp_bc')) then
        if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) == DIRICHLET_BC) then
          HAVE_GTEMP_TOPBC = PETSC_TRUE
        else
          option%io_buffer='pflotranModelSetTcond -  ' // &
               ' for CLM-PFLOTRAN coupling - flow condition MUST be named as following: ' // &
               ' "clm_gtemp_bc: TEMPERATURE dirichlet " for subsurface-top TYPE II  '
          call printErrMsg(option)
        endif

      endif
      boundary_condition => boundary_condition%next
    end do

    ! BC conditions update from CLM
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      ! Find appropriate BC from the list of boundary conditions
      ! TOP of subsurface
      ! the following thermal BC is accompanying for BC water flow
      if(StringCompare(boundary_condition%name,'clm_gpress_bc')) then

        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == DIRICHLET_BC) then
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gtemp_subsurf_pf_loc(iconn)                          ! so, this MUST be the temperature above first soil

          elseif(boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == ZERO_GRADIENT_BC) then
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gtemp_subsurf_pf_loc(iconn)                          ! so, this shall be the infiltration/water body temperature
                                                                               !    (TODO: currently just taking that same as above)
          else
            option%io_buffer='pflotranModelSetTcond -  ' // &
               ' for CLM-PFLOTRAN coupling - BC flow/temperature condition IS NOT: ' // &
               ' "TEMPERATURE dirichlet or zero_gradient" for subsurface-top '
            call printErrMsg(option)
          end if

        enddo
      endif
      !
      if(StringCompare(boundary_condition%name,'clm_gwflux_bc')) then

        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == NEUMANN_BC) then                    ! for liq. water flow BC, thermal conduction OFF, but with its own thermal states
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gtemp_subsurf_pf_loc(iconn)                          ! so, this shall be the infiltration/water body temperature
                                                                               !    (TODO: currently just taking that same as above)

          elseif (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == DIRICHLET_BC) then                    ! for liq. water flow BC, thermal conduction ON
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                gtemp_subsurf_pf_loc(iconn)

          elseif(boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == ZERO_GRADIENT_BC) then
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gtemp_subsurf_pf_loc(iconn)                          ! so, this shall be the infiltration/water body temperature
                                                                               !    (TODO: currently just taking that same as above)
          else
            option%io_buffer='pflotranModelSetTcond -  ' // &
               ' for CLM-PFLOTRAN coupling - BC flow/temperature condition "clm_gwflux_bc" MUST be: ' // &
               ' TEMPERATURE neumann/dirichlet or TEMPERATURE zero_gradient '
            call printErrMsg(option)

          end if

        enddo

      endif
      !
      if(StringCompare(boundary_condition%name,'clm_gevflux_bc')) then
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == NEUMANN_BC) then                    ! for liq. water flow BC, thermal conduction OFF
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                geflxl_subsurf_pf_loc(iconn)

          elseif (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == DIRICHLET_BC) then                  ! for liq. water flow BC, water boday temperature external
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                gtemp_subsurf_pf_loc(iconn)
          else

            option%io_buffer='pflotranModelSetTcond -  ' // &
               ' for CLM-PFLOTRAN coupling - BC flow/temperature condition "clm_gevflux_bc" MUST be: ' // &
               ' TEMPERATURE nenumann or dirichlet'
            call printErrMsg(option)

          end if

        enddo
      endif

      ! the following thermal BC is usually for energy only, i.e. no water flux
      if(StringCompare(boundary_condition%name,'clm_gtemp_bc')) then

        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == DIRICHLET_BC) then
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gtemp_subsurf_pf_loc(iconn)                          ! so, this MUST be the temperature above first soil

          else
            option%io_buffer='pflotranModelSetTHbcs -  ' // &
               ' for CLM-PFLOTRAN coupling - BC flow/temperature condition MUST be like : ' // &
               ' "TEMPERATURE dirichlet" for subsurface-top  TYPE II '
            call printErrMsg(option)
          end if

        enddo
      endif
      !
      if(StringCompare(boundary_condition%name,'clm_geflux_bc')) then

        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == NEUMANN_BC) then

            ! the net energy flux of all forms, G = Rnet + (-SH) + (-LE)
            ! (+ into soil, - out soil)

            if (HAVE_GEV_TOPBC .and. HAVE_GTEMP_TOPBC) then
              ! 'clm_gevflux_bc' is for LE, and 'clm_gtemp_bc' is for SH
              ! (i.e. gtemp MUST be temperature above first soil)
              boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          geflxr_subsurf_pf_loc(iconn)                         ! net radiation flux into first soil, i.e. excluding SH + LE

            elseif (HAVE_GEV_TOPBC) then
              ! LE shall not be included: G-(-LE)
              boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          geflx_subsurf_pf_loc(iconn) - &                      ! net heat flux into first soil, i.e. excluding SH + LE
                          geflxl_subsurf_pf_loc(iconn)                         ! excluding (negative) soil evap LE

            elseif (HAVE_GTEMP_TOPBC) then
              ! SH shall not be included: G-(-SH) = Rnet+(-LE)
              boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          geflxr_subsurf_pf_loc(iconn) + &                     ! net radiation flux into first soil
                          geflxl_subsurf_pf_loc(iconn)                         ! including (negative) soil evap LE

            else
              ! the net energy flux of all forms, i.e. ground heat flux G at soil interface, if ONLY known energy flux as BC
              boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          geflx_subsurf_pf_loc(iconn)                          ! so, this MUST be the all energy flux into first soil
            endif

          else
            option%io_buffer='pflotranModelSetTHbcs -  ' // &
               ' for CLM-PFLOTRAN coupling - BC flow/temperature condition MUST be like : ' // &
               ' "ENERGY_FLUX neumann" for subsurface-top  TYPE I '
            call printErrMsg(option)
          end if

        enddo
      endif

      ! BOTTOM (BASE) of subsurface
      ! the following thermal BC may or may not be accompanying for BC water flow
      if(StringCompare(boundary_condition%name,'clm_bflux_bc')) then
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
                == DIRICHLET_BC) then                  ! for water flow BC, thermal boundary MUST be 'dirichlet' type
            boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          gtemp_subbase_pf_loc(iconn)

          end if

        enddo
      endif

      boundary_condition => boundary_condition%next
    enddo

    ! Source/Sink terms
    source_sink => realization%patch%source_sink_list%first
    sum_connection = 0
    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate SrcSink from the list of flow conditions
      if(StringCompare(source_sink%name,'clm_ghf_ss')) then     ! heat flow rate as term of source/sink
        do iconn = 1, cur_connection_set%num_connections
          sum_connection = sum_connection + 1

          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          if (patch%imat(ghosted_id) < 0) cycle

          if (source_sink%flow_condition%itype(TH_TEMPERATURE_DOF) == ENERGY_RATE_SS) then

            source_sink%flow_condition%energy_rate%dataset%rarray(1) = &
                          geflow_sub_pf_loc(iconn)*material_auxvars(ghosted_id)%volume

          else if (source_sink%flow_condition%itype(TH_TEMPERATURE_DOF) == HET_ENERGY_RATE_SS) then

            source_sink%flow_aux_real_var(TH_PRESSURE_DOF,iconn) = 0.d0

            source_sink%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                          geflow_sub_pf_loc(iconn)*material_auxvars(ghosted_id)%volume

          end if

        enddo

      endif

      source_sink => source_sink%next

    enddo

    call VecRestoreArrayF90(clm_pf_idata%eflux_subsurf_pfs,geflx_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%efluxr_subsurf_pfs,geflxr_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%efluxl_subsurf_pfs,geflxl_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%gtemp_subsurf_pfs,gtemp_subsurf_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%eflux_subbase_pfs,geflx_subbase_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%gtemp_subbase_pfs,gtemp_subbase_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%eflow_pfs,geflow_sub_pf_loc,ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

  end subroutine pflotranModelUpdateSubsurfTCond

! ************************************************************************** !
  subroutine pflotranModelGetBCMassBalanceDeltaFromPF(pflotran_model)
  !
  ! Calculate mass balance at BC for passing flow rates to CLM
  !
  ! Author: Fengming Yuan
  ! Date: 03/14/2014
  !
    use Patch_module
    use Grid_module
    use Option_module
    use Connection_module
    use Coupler_module
    use Utility_module
    use String_module

    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use Global_Aux_module
    use Reactive_Transport_Aux_module
    use Reaction_Aux_module

    use clm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(Option_type), pointer :: option
    type(pflotran_model_type), pointer          :: pflotran_model
    class(simulation_subsurface_type), pointer  :: simulation
    class(realization_subsurface_type), pointer :: realization
    type(patch_type), pointer :: patch
    type(grid_type), pointer  :: grid

    type(coupler_type), pointer :: boundary_condition
    type(connection_set_type), pointer :: cur_connection_set
    type(global_auxvar_type), pointer :: global_auxvars_bc(:)
    type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)

    PetscReal, pointer :: qinfl_subsurf_pf_loc(:)
    PetscReal, pointer :: qsurf_subsurf_pf_loc(:)
    PetscReal, pointer :: qflux_subbase_pf_loc(:)
    PetscReal, pointer :: f_nh4_subsurf_pf_loc(:)
    PetscReal, pointer :: f_no3_subsurf_pf_loc(:)
    PetscReal, pointer :: f_nh4_subbase_pf_loc(:)
    PetscReal, pointer :: f_no3_subbase_pf_loc(:)

    PetscInt :: local_id, ghosted_id, iconn
    PetscInt :: offset
    PetscErrorCode :: ierr

    subname = 'ModelGetBCMassBalanceDeltaFromPF'
    !-------------------------------------------------------------------------
    option => pflotran_model%option
    select type (modelsim => pflotran_model%simulation)
      class is (simulation_subsurface_type)
        simulation  => modelsim
        realization => simulation%realization

      class default
        option%io_buffer = " subroutine is " // trim(subname) // &
              "currently is Not support in this simulation."
        call printErrMsg(option)
    end select
    !
    patch => realization%patch
    grid  => patch%grid
    option=> realization%option
    !-------------------------------------------------------------------------

    if (clm_pf_idata%nlpf_2dtop <= 0 .and. clm_pf_idata%ngpf_2dtop <= 0    &
        .and. clm_pf_idata%nlpf_2dbot <= 0 .and. clm_pf_idata%ngpf_2dbot <= 0) then
        return
    endif

    !
    call VecGetArrayF90(clm_pf_idata%qinfl_subsurf_pfp, qinfl_subsurf_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%qsurf_subsurf_pfp, qsurf_subsurf_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecGetArrayF90(clm_pf_idata%qflux_subbase_pfp, qflux_subbase_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    qinfl_subsurf_pf_loc(:) = 0.d0
    qsurf_subsurf_pf_loc(:) = 0.d0
    qflux_subbase_pf_loc(:) = 0.d0

    !
    boundary_condition => patch%boundary_condition_list%first
    global_auxvars_bc => patch%aux%Global%auxvars_bc
    if (option%ntrandof > 0) then
       rt_auxvars_bc => patch%aux%RT%auxvars_bc
    endif

    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      offset = cur_connection_set%offset

      if (option%nflowdof > 0) then

          ! retrieving H2O flux at top BC
          do iconn = 1, cur_connection_set%num_connections
             local_id = cur_connection_set%id_dn(iconn)
             ghosted_id = grid%nL2G(local_id)
             if (ghosted_id <= 0 .or. local_id <= 0) cycle
             if (patch%imat(ghosted_id) < 0) cycle

             if(StringCompare(boundary_condition%name,'clm_gpress_bc')) then          ! infilitration (+)
                qinfl_subsurf_pf_loc(iconn) = &
                               -global_auxvars_bc(offset+iconn)%mass_balance(1,1)

                ! 'mass_balance' IS accumulative, so need to reset to Zero for next desired time-step
                global_auxvars_bc(offset+iconn)%mass_balance(1,1) = 0.d0

             endif

             if(StringCompare(boundary_condition%name,'exfiltration')) then          ! exfiltration (-)
                qsurf_subsurf_pf_loc(iconn) = &
                               -global_auxvars_bc(offset+iconn)%mass_balance(1,1)

                 ! 'mass_balance' IS accumulative, so need to reset to Zero for next desired time-step
                global_auxvars_bc(offset+iconn)%mass_balance(1,1) = 0.d0

             endif

             ! retrieving H2O flux at bottom BC
             if(StringCompare(boundary_condition%name,'clm_bflux_bc')) then          ! bottom water flux
                qflux_subbase_pf_loc(iconn) = &
                               -global_auxvars_bc(offset+iconn)%mass_balance(1,1)

                ! 'mass_balance' IS accumulative, so need to reset to Zero for next desired time-step
                global_auxvars_bc(offset+iconn)%mass_balance(1,1) = 0.d0

             endif

          enddo


      endif

      boundary_condition => boundary_condition%next

    enddo

    call VecRestoreArrayF90(clm_pf_idata%qinfl_subsurf_pfp, qinfl_subsurf_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%qsurf_subsurf_pfp, qsurf_subsurf_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecRestoreArrayF90(clm_pf_idata%qflux_subbase_pfp, qflux_subbase_pf_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! pass vecs to CLM
    if (clm_pf_idata%nlpf_2dtop > 0 .and. clm_pf_idata%ngpf_2dtop > 0 ) then
      call MappingSourceToDestination(pflotran_model%map_pf_2dtop_to_clm_2dtop, &
                                    option, &
                                    clm_pf_idata%qinfl_subsurf_pfp, &
                                    clm_pf_idata%qinfl_subsurf_clms)

      call MappingSourceToDestination(pflotran_model%map_pf_2dtop_to_clm_2dtop, &
                                    option, &
                                    clm_pf_idata%qsurf_subsurf_pfp, &
                                    clm_pf_idata%qsurf_subsurf_clms)
    endif

    if (clm_pf_idata%nlpf_2dbot > 0 .and. clm_pf_idata%ngpf_2dbot > 0 ) then
      call MappingSourceToDestination(pflotran_model%map_pf_2dbot_to_clm_2dbot, &
                                    option, &
                                    clm_pf_idata%qflux_subbase_pfp, &
                                    clm_pf_idata%qflux_subbase_clms)
    endif
  end subroutine pflotranModelGetBCMassBalanceDeltaFromPF

  ! ************************************************************************** !
  
end module pflotran_clm_main_module

#endif
