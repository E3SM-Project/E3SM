module PM_Waste_Form_class

  use PM_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Geometry_module
  use Data_Mediator_Vec_class
  use Dataset_Base_class
  use Region_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscBool, public :: bypass_warning_message = PETSC_FALSE

! --------------- waste form species packages ---------------------------------
  type, public :: rad_species_type
   PetscReal :: formula_weight
   PetscReal :: decay_constant
   PetscReal :: mass_fraction
   PetscReal :: inst_release_fraction
   PetscInt :: daugh_id
   character(len=MAXWORDLENGTH) :: daughter
   PetscInt :: column_id
   PetscInt :: ispecies
   character(len=MAXWORDLENGTH) :: name
  end type rad_species_type

! --------------- waste form mechanism types ----------------------------------
  type :: wf_mechanism_base_type
    type(rad_species_type), pointer :: rad_species_list(:)
    PetscInt :: num_species
    PetscBool :: canister_degradation_model
    PetscReal :: vitality_rate_mean
    PetscReal :: vitality_rate_stdev
    PetscReal :: vitality_rate_trunc
    PetscReal :: canister_material_constant
    PetscReal :: matrix_density                 ! kg/m^3
    PetscReal :: specific_surface_area          ! m^2/kg
    character(len=MAXWORDLENGTH) :: name
    class(wf_mechanism_base_type), pointer :: next
  contains
    procedure, public :: Dissolution => WFMechBaseDissolution
  end type wf_mechanism_base_type

  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_glass_type
    PetscReal :: dissolution_rate    ! kg-glass/m^2/sec
    PetscReal :: k0                  ! k-glass/m^2/day
    PetscReal :: k_long              ! k-glass/m^2/day
    PetscReal :: nu                  ! [-]
    PetscReal :: Ea                  ! [J/mol]
    PetscReal :: Q
    PetscReal :: K
    PetscReal :: v
    PetscReal :: pH
    PetscBool :: use_pH
    PetscBool :: use_Q
    PetscInt :: h_ion_id
    PetscInt :: SiO2_id
  contains
    procedure, public :: Dissolution => WFMechGlassDissolution
  end type wf_mechanism_glass_type

  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_dsnf_type
    PetscReal :: frac_dissolution_rate    ! 1/sec
  contains
    procedure, public :: Dissolution => WFMechDSNFDissolution
  end type wf_mechanism_dsnf_type
  
  ! the WIPP mechanism is the same as an instantaneous dissolution type
  ! when selecting for DSNF and WIPP together, class is() can be used
  ! when selecting for either DSNF or WIPP, type is() should be used
  type, public, extends(wf_mechanism_dsnf_type) :: wf_mechanism_wipp_type
  end type wf_mechanism_wipp_type

  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_fmdm_type
    PetscReal :: dissolution_rate         ! kg-matrix/m^2/sec
    PetscReal :: frac_dissolution_rate    ! 1/sec
    PetscReal :: burnup                   ! GWd/MTHM (kg-matrix/m^2/sec)
    PetscInt :: num_grid_cells_in_waste_form  ! hardwired to 40
    ! mapping of fmdm species into fmdm concentration array:
    PetscInt, pointer :: mapping_fmdm(:)
    ! mapping of species in fmdm concentration array to pflotran:
    PetscInt, pointer :: mapping_fmdm_to_pflotran(:)
    PetscReal, pointer :: concentration(:,:)
    PetscInt :: num_concentrations            ! hardwired to 11
    PetscInt :: iUO2_2p
    PetscInt :: iUCO3_2n
    PetscInt :: iUO2
    PetscInt :: iCO3_2n
    PetscInt :: iO2
    PetscInt :: iH2O2
    PetscInt :: iFe_2p
    PetscInt :: iH2
    PetscInt :: iUO2_sld
    PetscInt :: iUO3_sld
    PetscInt :: iUO4_sld
  contains
    procedure, public :: Dissolution => WFMechFMDMDissolution
  end type wf_mechanism_fmdm_type
  
  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_custom_type
    PetscReal :: dissolution_rate         ! kg-matrix/m^2/sec
    PetscReal :: frac_dissolution_rate    ! 1/sec
  contains
    procedure, public :: Dissolution => WFMechCustomDissolution
  end type wf_mechanism_custom_type

! --------------- waste form types --------------------------------------------
  type :: waste_form_base_type
    PetscInt :: id
    PetscMPIInt :: myMPIgroup_id
    PetscMPIInt :: myMPIcomm
    type(point3d_type) :: coordinate
    character(len=MAXWORDLENGTH) :: region_name
    type(region_type), pointer :: region
    PetscReal, pointer :: scaling_factor(:)             ! [-]
    PetscReal :: init_volume                            ! m^3
    PetscReal :: volume                                 ! m^3
    PetscReal :: exposure_factor                        ! unitless 
    PetscReal :: eff_dissolution_rate                   ! kg-matrix/sec
    PetscReal, pointer :: instantaneous_mass_rate(:)    ! mol/sec
    PetscReal, pointer :: cumulative_mass(:)            ! mol
    PetscReal, pointer :: rad_mass_fraction(:)          ! g-rad/g-matrix
    PetscReal, pointer :: rad_concentration(:)          ! mol-rad/g-matrix
    PetscReal, pointer :: inst_release_amount(:)        ! of rad
    PetscBool :: canister_degradation_flag
    PetscReal :: canister_vitality                      ! %
    PetscReal :: canister_vitality_rate
    PetscReal :: eff_canister_vit_rate
    PetscReal :: breach_time                            ! sec
    PetscBool :: breached
    character(len=MAXWORDLENGTH) :: mech_name
    class(wf_mechanism_base_type), pointer :: mechanism
    class(waste_form_base_type), pointer :: next
  end type waste_form_base_type

! --------------- waste form process model ------------------------------------
  type, public, extends(pm_base_type) :: pm_waste_form_type
    class(realization_subsurface_type), pointer :: realization
    character(len=MAXWORDLENGTH) :: data_mediator_species
    class(data_mediator_vec_type), pointer :: data_mediator
    class(waste_form_base_type), pointer :: waste_form_list
    class(wf_mechanism_base_type), pointer :: mechanism_list
    PetscBool :: print_mass_balance
  contains
    procedure, public :: PMWFSetRealization
    procedure, public :: Setup => PMWFSetup
    procedure, public :: Read => PMWFRead
    procedure, public :: InitializeRun => PMWFInitializeRun
    procedure, public :: InitializeTimestep => PMWFInitializeTimestep
    procedure, public :: FinalizeTimestep => PMWFFinalizeTimestep
    procedure, public :: UpdateSolution => PMWFUpdateSolution
    procedure, public :: Solve => PMWFSolve
    procedure, public :: Checkpoint => PMWFCheckpoint    
    procedure, public :: Restart => PMWFRestart  
    procedure, public :: InputRecord => PMWFInputRecord
    procedure, public :: Destroy => PMWFDestroy
  end type pm_waste_form_type
  
  public :: PMWFCreate, &
            PMWFSetup, &
            MechanismGlassCreate, &
            MechanismDSNFCreate, &
            MechanismWIPPCreate, &
            MechanismCustomCreate, &
            MechanismFMDMCreate, &
            RadSpeciesCreate
  
contains

! ************************************************************************** !

subroutine MechanismInit(this)
  ! 
  ! Initializes the base waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

  class(wf_mechanism_base_type) :: this

  nullify(this%next)
  nullify(this%rad_species_list)
  this%num_species = 0
  this%matrix_density = UNINITIALIZED_DOUBLE
  this%specific_surface_area = UNINITIALIZED_DOUBLE
  this%name = ''
 !---- canister degradation model ----------------------
  this%canister_degradation_model = PETSC_FALSE
  this%vitality_rate_mean = UNINITIALIZED_DOUBLE
  this%vitality_rate_stdev = UNINITIALIZED_DOUBLE
  this%vitality_rate_trunc = UNINITIALIZED_DOUBLE
  this%canister_material_constant = UNINITIALIZED_DOUBLE
 !------------------------------------------------------

end subroutine MechanismInit

! ************************************************************************** !

function MechanismGlassCreate()
  ! 
  ! Creates the glass waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_glass_type), pointer :: MechanismGlassCreate
  
  allocate(MechanismGlassCreate)
  call MechanismInit(MechanismGlassCreate)
  MechanismGlassCreate%dissolution_rate = 0.d0        ! [kg/m^2/sec]
  MechanismGlassCreate%k0 = UNINITIALIZED_DOUBLE      ! [kg/m^2/sec]
  MechanismGlassCreate%k_long = UNINITIALIZED_DOUBLE  ! [kg/m^2/sec]
  MechanismGlassCreate%nu = UNINITIALIZED_DOUBLE      ! [-]
  MechanismGlassCreate%Ea = UNINITIALIZED_DOUBLE      ! [J/mol]
  MechanismGlassCreate%Q = UNINITIALIZED_DOUBLE      
  MechanismGlassCreate%K = UNINITIALIZED_DOUBLE
  MechanismGlassCreate%v = UNINITIALIZED_DOUBLE    
  MechanismGlassCreate%pH = UNINITIALIZED_DOUBLE   
  MechanismGlassCreate%use_pH = PETSC_FALSE  
  MechanismGlassCreate%use_Q = PETSC_FALSE
  MechanismGlassCreate%h_ion_id = 0
  MechanismGlassCreate%SiO2_id = 0

end function MechanismGlassCreate

! ************************************************************************** !

function MechanismDSNFCreate()
  ! 
  ! Creates the DSNF (DOE Spent Nuclear Fuel) waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_dsnf_type), pointer :: MechanismDSNFCreate
  
  allocate(MechanismDSNFCreate)
  call MechanismInit(MechanismDSNFCreate)
  MechanismDSNFCreate%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/sec

end function MechanismDSNFCreate

! ************************************************************************** !

function MechanismWIPPCreate()
  ! 
  ! Creates the WIPP (Waste Isolation Pilot Plant) waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 012/7/2016

  implicit none
  
  class(wf_mechanism_wipp_type), pointer :: MechanismWIPPCreate
  
  allocate(MechanismWIPPCreate)
  call MechanismInit(MechanismWIPPCreate)
  MechanismWIPPCreate%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/sec

end function MechanismWIPPCreate

! ************************************************************************** !

function MechanismFMDMCreate()
  ! 
  ! Creates the FMDM waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_fmdm_type), pointer :: MechanismFMDMCreate
  
  allocate(MechanismFMDMCreate)
  call MechanismInit(MechanismFMDMCreate)
  
  MechanismFMDMCreate%dissolution_rate = UNINITIALIZED_DOUBLE       ! kg/m^2/sec
  MechanismFMDMCreate%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/day
  MechanismFMDMCreate%burnup = UNINITIALIZED_DOUBLE     ! GWd/MTHM or (kg/m^2/sec)
  
  MechanismFMDMCreate%num_grid_cells_in_waste_form = 40  ! hardwired
  
  nullify(MechanismFMDMCreate%concentration)
  MechanismFMDMCreate%num_concentrations = 11         
  MechanismFMDMCreate%iUO2_2p = 1
  MechanismFMDMCreate%iUCO3_2n = 2
  MechanismFMDMCreate%iUO2 = 3
  MechanismFMDMCreate%iCO3_2n = 4
  MechanismFMDMCreate%iO2 = 5
  MechanismFMDMCreate%iH2O2 = 6
  MechanismFMDMCreate%iFe_2p = 7
  MechanismFMDMCreate%iH2 = 8
  MechanismFMDMCreate%iUO2_sld = 9
  MechanismFMDMCreate%iUO3_sld = 10
  MechanismFMDMCreate%iUO4_sld = 11
  
  allocate(MechanismFMDMCreate%mapping_fmdm_to_pflotran( &
           MechanismFMDMCreate%num_concentrations))
  MechanismFMDMCreate%mapping_fmdm_to_pflotran = UNINITIALIZED_INTEGER
  
  ! concentration can be allocated here because we hardwired
  ! the num_grid_cells_in_waste_form value, but if it becomes
  ! user defined, then allocation must be delayed until PMWFSetup
  allocate(MechanismFMDMCreate%concentration( &
             MechanismFMDMCreate%num_concentrations, &
             MechanismFMDMCreate%num_grid_cells_in_waste_form))
  MechanismFMDMCreate%concentration = 1.d-13
  
  allocate(MechanismFMDMCreate%mapping_fmdm(4))
  MechanismFMDMCreate%mapping_fmdm = [MechanismFMDMCreate%iO2, &
                                      MechanismFMDMCreate%iCO3_2n, &
                                      MechanismFMDMCreate%iH2, &
                                      MechanismFMDMCreate%iFe_2p]

end function MechanismFMDMCreate

! ************************************************************************** !

function MechanismCustomCreate()
  ! 
  ! Creates the 'custom' waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_custom_type), pointer :: MechanismCustomCreate
  
  allocate(MechanismCustomCreate)
  call MechanismInit(MechanismCustomCreate)
  MechanismCustomCreate%dissolution_rate = UNINITIALIZED_DOUBLE    ! kg/m^2/sec
  MechanismCustomCreate%frac_dissolution_rate = UNINITIALIZED_DOUBLE    ! 1/sec

end function MechanismCustomCreate


! ************************************************************************** !

function RadSpeciesCreate()
  ! 
  ! Creates a radioactive species in the waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/09/16

  implicit none
  
  type(rad_species_type) :: RadSpeciesCreate

  RadSpeciesCreate%name = ''
  RadSpeciesCreate%daughter = ''
  RadSpeciesCreate%daugh_id = UNINITIALIZED_INTEGER
  RadSpeciesCreate%formula_weight = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%decay_constant = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%mass_fraction = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%inst_release_fraction = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%column_id = UNINITIALIZED_INTEGER
  RadSpeciesCreate%ispecies = UNINITIALIZED_INTEGER

end function RadSpeciesCreate

! ************************************************************************** !

function WasteFormCreate()
  ! 
  ! Creates a waste form and initializes all parameters
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

  type(waste_form_base_type), pointer :: WasteFormCreate

  allocate(WasteFormCreate)
  WasteFormCreate%id = UNINITIALIZED_INTEGER
  WasteFormCreate%myMPIgroup_id = 0
  WasteFormCreate%myMPIcomm = 0
  WasteFormCreate%coordinate%x = UNINITIALIZED_DOUBLE
  WasteFormCreate%coordinate%y = UNINITIALIZED_DOUBLE
  WasteFormCreate%coordinate%z = UNINITIALIZED_DOUBLE
  nullify(WasteFormCreate%region)
  WasteFormCreate%region_name = ''
  nullify(WasteFormCreate%scaling_factor) ! [-]
  WasteFormCreate%init_volume = UNINITIALIZED_DOUBLE
  WasteFormCreate%volume = UNINITIALIZED_DOUBLE
  WasteFormCreate%exposure_factor = 1.0d0
  WasteFormCreate%eff_dissolution_rate = UNINITIALIZED_DOUBLE
  WasteFormCreate%mech_name = ''
  nullify(WasteFormCreate%instantaneous_mass_rate) ! mol-rad/sec
  nullify(WasteFormCreate%cumulative_mass) ! mol-rad
  nullify(WasteFormCreate%rad_mass_fraction) ! g-rad/g-matrix
  nullify(WasteFormCreate%rad_concentration) ! mol-rad/g-matrix
  nullify(WasteFormCreate%inst_release_amount) ! of rad
  nullify(WasteFormCreate%mechanism)
  nullify(WasteFormCreate%next)
 !------- canister degradation model -----------------
  WasteFormCreate%canister_degradation_flag = PETSC_FALSE
  WasteFormCreate%breached = PETSC_FALSE
  WasteFormCreate%breach_time = UNINITIALIZED_DOUBLE
  WasteFormCreate%canister_vitality = 0.d0
  WasteFormCreate%canister_vitality_rate = UNINITIALIZED_DOUBLE
  WasteFormCreate%eff_canister_vit_rate = UNINITIALIZED_DOUBLE
 !----------------------------------------------------

end function WasteFormCreate

! ************************************************************************** !

function PMWFCreate()
  ! 
  ! Creates and initializes the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15, 07/20/15
  ! Notes: Modified by Jenn Frederick 03/24/2016

  implicit none
  
  class(pm_waste_form_type), pointer :: PMWFCreate
  
  allocate(PMWFCreate)
  nullify(PMWFCreate%realization)
  nullify(PMWFCreate%data_mediator)
  nullify(PMWFCreate%waste_form_list)
  nullify(PMWFCreate%mechanism_list)  
  PMWFCreate%print_mass_balance = PETSC_FALSE
  PMWFCreate%name = 'waste form general'

  call PMBaseInit(PMWFCreate)

end function PMWFCreate

! ************************************************************************** !

subroutine PMWFRead(this,input)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick, 03/24/2016

  use Input_Aux_module
  use Option_module
  use String_module
  use Region_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  
  class(waste_form_base_type), pointer :: cur_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: matched

  option => this%option
  input%ierr = 0
  error_string = 'WASTE_FORM_GENERAL'

  option%io_buffer = 'pflotran card:: ' // trim(error_string)
  call printMsg(option)

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadMechanism(this,input,option,word,error_string,found)
    if (found) cycle
    
    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadWasteForm(this,input,option,word,error_string,found)
    if (found) cycle
   
    select case(trim(word))
    !-------------------------------------
      case('PRINT_MASS_BALANCE')
        this%print_mass_balance = PETSC_TRUE
    !-------------------------------------
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !-------------------------------------
    end select
  enddo

  ! Assign chosen mechanism to each waste form
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    matched = PETSC_FALSE
    cur_mechanism => this%mechanism_list
    do
      if (.not.associated(cur_mechanism)) exit
      if (StringCompare(trim(cur_waste_form%mech_name), &
                        trim(cur_mechanism%name))) then
        cur_waste_form%mechanism => cur_mechanism
        matched = PETSC_TRUE
      endif
      if (matched) exit
      cur_mechanism => cur_mechanism%next
    enddo
    ! error messaging: ----------------------------------------------
    if (.not.associated(cur_waste_form%mechanism)) then
      option%io_buffer = 'WASTE_FORM MECHANISM ' // &
                         trim(cur_waste_form%mech_name) // &
                         ' not found amoung given mechanism names.'
      call printErrMsg(option)
    endif
    
    if (.not.cur_waste_form%mechanism%canister_degradation_model) then
      ! canister vitality specified, but can.deg. model is off:
      if (initialized(cur_waste_form%canister_vitality_rate)) then
        option%io_buffer = 'WASTE_FORM MECHANISM ' // &
          trim(cur_waste_form%mech_name) // ' does not have the canister &
          &degradation model turned on, but at least one of the waste forms &
          &assigned to this mechanism specifies a canister vitality rate.'
        call printErrMsg(option)
      endif
      ! canister breach time specified, but can.deg. model is off:
      if (initialized(cur_waste_form%breach_time)) then
        option%io_buffer = 'WASTE_FORM MECHANISM ' // &
          trim(cur_waste_form%mech_name) // ' does not have the canister &
          &degradation model turned on, but at least one of the waste forms &
          &assigned to this mechanism specifies a canister breach time.'
        call printErrMsg(option)
      endif
    endif

    ! both waste form and mechanism canister vitality rate parameters 
    ! are specified:
    if (initialized(cur_waste_form%canister_vitality_rate) .and. &
        ( initialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
          initialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
          initialized(cur_waste_form%mechanism%vitality_rate_trunc) )) then
      option%io_buffer = 'Either CANISTER_VITALITY_RATE within the &
        &WASTE_FORM blocks -or- the VITALITY_LOG10_MEAN, &
        &VITALITY_LOG10_STDEV, and VITALITY_UPPER_TRUNCATION within &
        &the WASTE_FORM MECHANISM ' // trim(cur_waste_form%mechanism%name) &
        // ' block should be specified, but not both.'
      call printErrMsg(option)
    endif
    
    ! the canister degradation model is on, but there are problems with
    ! the parameters provided:
    if (cur_waste_form%mechanism%canister_degradation_model) then 
      ! all parameters are missing:
      if ( (uninitialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
            uninitialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
            uninitialized(cur_waste_form%mechanism%vitality_rate_trunc) ) .and. &
          uninitialized(cur_waste_form%canister_vitality_rate) .and. &
          uninitialized(cur_waste_form%breach_time)                 )  then 
        option%io_buffer = 'CANISTER_VITALITY_RATE within the WASTE_FORM &
          &blocks -or- CANISTER_BREACH_TIME within the WASTE_FORM blocks &
          &-or- the VITALITY_LOG10_MEAN, VITALITY_LOG10_STDEV, and &
          &VITALITY_UPPER_TRUNCATION within the WASTE_FORM MECHANISM ' // &
          trim(cur_waste_form%mechanism%name) // ' is missing.'
        call printErrMsg(option)
      endif
      ! all parameters are given:
      if ( (initialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
            initialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
            initialized(cur_waste_form%mechanism%vitality_rate_trunc) ) .and. &
          initialized(cur_waste_form%canister_vitality_rate) .and. &
          initialized(cur_waste_form%breach_time)                 )  then 
        option%io_buffer = 'CANISTER_VITALITY_RATE within the WASTE_FORM &
          &blocks -or- CANISTER_BREACH_TIME within the WASTE_FORM blocks &
          &-or- the VITALITY_LOG10_MEAN, VITALITY_LOG10_STDEV, and &
          &VITALITY_UPPER_TRUNCATION within the WASTE_FORM MECHANISM ' // &
          trim(cur_waste_form%mechanism%name) // ' should be specified, &
          &but not all.'
        call printErrMsg(option)
      endif
      ! both breach time and can. deg. rate were given
      if (initialized(cur_waste_form%canister_vitality_rate) .and. &
          initialized(cur_waste_form%breach_time)) then
        option%io_buffer = 'Either CANISTER_VITALITY_RATE -or- &
          &CANISTER_BREACH_TIME within the WASTE_FORM block with &
          &WASTE_FORM MECHANISM ' // trim(cur_waste_form%mechanism%name) &
          // ' should be specified, but not both.'
        call printErrMsg(option)
      endif
      ! both breach time and can. deg. distribution were given
      if ((initialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
           initialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
           initialized(cur_waste_form%mechanism%vitality_rate_trunc)) .and. &
          initialized(cur_waste_form%breach_time)) then
        option%io_buffer = 'Either CANISTER_BREACH_TIME within the &
          &WASTE_FORM block with WASTE_FORM MECHANISM ' &
          // trim(cur_waste_form%mechanism%name) // ' -or- the &
          &VITALITY_LOG10_MEAN, VITALITY_LOG10_STDEV, and &
          &VITALITY_UPPER_TRUNCATION within the WASTE_FORM MECHANISM ' // &
          trim(cur_waste_form%mechanism%name) // ' should be specified, &
          &but not both.'
        call printErrMsg(option)
      endif
    endif
    
    cur_waste_form => cur_waste_form%next
  enddo
    
end subroutine PMWFRead

! ************************************************************************** !

subroutine PMWFReadMechanism(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the waste form mechanism
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016
  !
  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use String_module
  use Units_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscBool :: added
  character(len=MAXWORDLENGTH) :: word, units, internal_units
  character(len=MAXSTRINGLENGTH) :: temp_buf, string
  type(rad_species_type), pointer :: temp_species_array(:)
  class(wf_mechanism_base_type), pointer :: new_mechanism, cur_mechanism
  PetscInt :: k, j, icol
  PetscReal :: double

  error_string = trim(error_string) // ',MECHANISM'
  found = PETSC_TRUE
  added = PETSC_FALSE
  input%ierr = 0
  allocate(temp_species_array(50))
  k = 0

  select case(trim(keyword))
  !-------------------------------------
    case('MECHANISM')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'mechanism type',error_string)
      call StringToUpper(word)
      select case(trim(word))
      !---------------------------------
        case('GLASS')
          error_string = trim(error_string) // ' GLASS'
          allocate(new_mechanism)
          new_mechanism => MechanismGlassCreate()
      !---------------------------------
        case('DSNF')
          error_string = trim(error_string) // ' DSNF'
          allocate(new_mechanism)
          new_mechanism => MechanismDSNFCreate()
      !---------------------------------
        case('WIPP')
          error_string = trim(error_string) // ' WIPP'
          allocate(new_mechanism)
          new_mechanism => MechanismWIPPCreate()
      !---------------------------------
        case('FMDM')
          ! for now, set bypass_warning_message = TRUE so we can run 
          ! the fmdm model even though its not included/linked 
          bypass_warning_message = PETSC_TRUE
#ifndef FMDM_MODEL
          this%option%io_buffer = 'Preprocessing statement FMDM_MODEL must &
            &be defined and the ANL FMDM library must be linked to PFLOTRAN &
            &to employ the fuel matrix degradation model.'
          if (.not.bypass_warning_message) then
            call printErrMsg(this%option)
          endif
#endif
          error_string = trim(error_string) // ' FMDM'
          allocate(new_mechanism)
          new_mechanism => MechanismFMDMCreate()
      !---------------------------------
        case('CUSTOM')
          error_string = trim(error_string) // ' CUSTOM'
          allocate(new_mechanism)
          new_mechanism => MechanismCustomCreate()
      !---------------------------------
        case default
          option%io_buffer = 'Unrecognized mechanism type &
                             &in the ' // trim(error_string) // ' block.'
          call printErrMsg(option)
      !---------------------------------
      end select
      
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !--------------------------
          case('NAME')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'mechanism name',error_string)
            call StringToUpper(word)
            new_mechanism%name = trim(word)
        !--------------------------
          case('SPECIFIC_SURFACE_AREA')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'specific surface area', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'m^2/kg', &
                            trim(error_string)//',specific surface area',option)
            select type(new_mechanism)
              class is(wf_mechanism_dsnf_type)
                ! applies to dsnf & wipp types
                option%io_buffer = 'SPECIFIC_SURFACE_AREA cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
              class default
                new_mechanism%specific_surface_area = double
            end select
        !--------------------------
          case('MATRIX_DENSITY')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'matrix density',error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^3', &
                                   trim(error_string)//',matrix density',option)
            select type(new_mechanism)
              class default
                new_mechanism%matrix_density = double
            end select
        !--------------------------
          case('FRACTIONAL_DISSOLUTION_RATE')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'fractional dissolution rate', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'unitless/sec', &
                      trim(error_string)//',fractional dissolution rate',option)
            select type(new_mechanism)
              type is(wf_mechanism_custom_type)
                new_mechanism%frac_dissolution_rate = double
              class default
                option%io_buffer = 'FRACTIONAL_DISSOLUTION_RATE cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('DISSOLUTION_RATE')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'dissolution rate',error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^2-sec', &
                                 trim(error_string)//',dissolution_rate',option)
            select type(new_mechanism)
              type is(wf_mechanism_custom_type)
                new_mechanism%dissolution_rate = double
              class default
                option%io_buffer = 'DISSOLUTION_RATE cannot be specified for ' &
                                   // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('K0')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'K0 (intrinsic dissolution rate)', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^2-sec', &
                  trim(error_string)//',K0 (intrinsic dissolution rate)',option)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%k0 = double
              class default
                option%io_buffer = 'K0 (intrinsic dissolution rate) cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('K_LONG')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'K_LONG (dissolution rate)', &
                               error_string)
            call InputReadAndConvertUnits(input,double,'kg/m^2-sec', &
                    trim(error_string)//',K_LONG (dissolution rate)',option)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%k_long = double
              class default
                option%io_buffer = 'K_LONG (dissolution rate) cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('NU')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'NU (pH dependence parameter)', &
                               error_string)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%nu = double
              class default
                option%io_buffer = 'NU (pH dependence parameter) cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('EA')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'EA (effective activation energy)',&
                               error_string)
            call InputReadAndConvertUnits(input,double,'J/mol', &
                 trim(error_string)//',EA (effective activation energy)',option)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%Ea = double
              class default
                option%io_buffer = 'EA (effective activation energy) cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('Q')
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                temp_buf = input%buf
                call InputReadDouble(input,option,double)
                if (InputError(input)) then
                  word = adjustl(trim(temp_buf))
                  call StringToUpper(word)
                  if (trim(word) == 'AS_CALCULATED') then 
                    new_mechanism%use_Q = PETSC_TRUE
                  else
                    option%io_buffer = 'Q value (ion activity product) was &
                                       &not provided, or Q instructions not &
                                       &understood for ' // trim(error_string)
                    call printErrMsg(option)
                  endif
                endif
                if (new_mechanism%use_Q) then
                  new_mechanism%Q = 0.d0  ! initializes to value
                else
                  new_mechanism%Q = double
                endif
              class default
                option%io_buffer = 'Q (ion activity product) cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('K')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'K (equilibrium constant)',&
                               error_string)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%K = double
              class default
                option%io_buffer = 'K (equilibrium constant) cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('V')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'V (exponent parameter)',&
                               error_string)
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%v = double
              class default
                option%io_buffer = 'V (exponent parameter) cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('PH')
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                temp_buf = input%buf
                call InputReadDouble(input,option,double)
                if (InputError(input)) then
                  word = adjustl(trim(temp_buf))
                  call StringToUpper(word)
                  if (trim(word) == 'AS_CALCULATED') then 
                    new_mechanism%use_pH = PETSC_TRUE
                  else
                    option%io_buffer = 'PH value was not provided, or PH &
                                       &instructions not understood for ' &
                                       // trim(error_string) // '.'
                    call printErrMsg(option)
                  endif
                endif
                if (new_mechanism%use_pH) then
                  new_mechanism%pH = 7.d0  ! initializes to neutral value
                else
                  new_mechanism%pH = double
                endif
              class default
                option%io_buffer = 'PH cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('KIENZLER_DISSOLUTION')
            select type(new_mechanism)
              type is(wf_mechanism_glass_type)
                new_mechanism%k0 = 560.d0/(24.d0*3600.d0)  ! kg/m^2-sec
                new_mechanism%k_long = 0.d0
                new_mechanism%nu = 0.d0
                new_mechanism%Ea = 7397.d0*8.314d0
                new_mechanism%Q = 0.d0
                new_mechanism%K = 1.d0     ! This value doesn't matter since Q=0
                new_mechanism%v = 1.d0
                new_mechanism%pH = 0.d0
              class default
                option%io_buffer = 'KIENZLER_DISSOLUTION cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('BURNUP')
            select type(new_mechanism)
              type is(wf_mechanism_fmdm_type)
                call InputReadDouble(input,option,new_mechanism%burnup)
                call InputErrorMsg(input,option,'burnup',error_string)
#ifndef FMDM_MODEL
                ! if fmdm model is not on, then burnup is dissolution rate
                call InputReadAndConvertUnits(input,new_mechanism%burnup, &
                                'kg/m^2-sec',trim(error_string)//',burnup', &
                                option)
                option%io_buffer = 'Warning: FMDM is not linked, but an &
                                   &FMDM mechanism was defined. BURNUP &
                                   &will be used for fuel dissolution rate.'
                call printMsg(option)
#else
                option%io_buffer = 'FMDM is linked.'
                call printMsg(option)
#endif
              class default
                option%io_buffer = 'BURNUP cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('SPECIES')
            do
              call InputReadPflotranString(input,option)
              if (InputCheckExit(input,option)) exit
              k = k + 1
              if (k > 50) then
                option%io_buffer = 'More than 50 radionuclide species are &
                                 &provided in the ' // trim(error_string) // &
                                 ', SPECIES block. Reduce to less than 50 &
                                 &species, or e-mail pflotran-dev at &
                                 &googlegroups dot com.'
                call printErrMsg(option)
              endif
              temp_species_array(k) = RadSpeciesCreate() 
              ! read species name
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'SPECIES name',error_string)
              temp_species_array(k)%name = trim(word)
              ! read species formula weight [g-species/mol-species]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES formula weight', &
                                 error_string)
              temp_species_array(k)%formula_weight = double
              ! read species decay constant [1/sec]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES decay rate constant', &
                                 error_string)
              temp_species_array(k)%decay_constant = double
              ! read species initial mass fraction [g-species/g-bulk]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES initial mass fraction', &
                                 error_string)
              temp_species_array(k)%mass_fraction = double
              ! read species instant release fraction [-]
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES instant release &
                                 &fraction',error_string)
              temp_species_array(k)%inst_release_fraction = double
              ! read species daughter
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                temp_species_array(k)%daughter = trim(word)
              else
                temp_species_array(k)%daughter = 'no_daughter'
              endif
              new_mechanism%num_species = k
            enddo
            if (k == 0) then
              option%io_buffer = 'At least one radionuclide species must be &
                                 &provided in the ' // trim(error_string) // &
                                 ', SPECIES block.'
              call printErrMsg(option)
            endif
            allocate(new_mechanism%rad_species_list(k))
            new_mechanism%rad_species_list(1:k) = temp_species_array(1:k)
            deallocate(temp_species_array)
            k = 0
            do while (k < new_mechanism%num_species)
              k = k + 1
              if (trim(new_mechanism%rad_species_list(k)%daughter) == &
                  'no_daughter') then
                new_mechanism%rad_species_list(k)%daugh_id = 0
              else
                j = 0
                do while (j < new_mechanism%num_species)
                  j = j + 1
                  if (trim(new_mechanism%rad_species_list(k)%daughter) == &
                       trim(new_mechanism%rad_species_list(j)%name)) then
                    new_mechanism%rad_species_list(k)%daugh_id = j
                    exit
                  endif
                enddo
              endif
            enddo
        !--------------------------
          case('CANISTER_DEGRADATION_MODEL')
            new_mechanism%canister_degradation_model = PETSC_TRUE
            do
              call InputReadPflotranString(input,option)
              if (InputCheckExit(input,option)) exit
              call InputReadWord(input,option,word,PETSC_TRUE)
              call StringToUpper(word)
              select case(trim(word))
              case('VITALITY_LOG10_MEAN')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_mean)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &mean value',error_string)
              case('VITALITY_LOG10_STDEV')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_stdev)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &st. dev. value',error_string)
              case('VITALITY_UPPER_TRUNCATION')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_trunc)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &upper truncation value',error_string)
              case('CANISTER_MATERIAL_CONSTANT')
                call InputReadDouble(input,option, &
                                     new_mechanism%canister_material_constant)
                call InputErrorMsg(input,option,'canister material constant', &
                                   error_string)
              case default
                option%io_buffer = 'Keyword ' // trim(word) // &
                                   ' not recognized in the ' // &
                                   trim(error_string) // &
                                   ' CANISTER_DEGRADATION_MODEL block.'
                call printErrMsg(option)
              end select
            enddo
        !--------------------------
          case default
            call InputKeywordUnrecognized(word,error_string,option)
        !--------------------------
        end select
      enddo

     !----------- error messaging ----------------------------------------------
      if (new_mechanism%name == '') then
        option%io_buffer = 'NAME must be specified in ' // trim(error_string) &
                           // ' block.'
        call printErrMsg(option)
      endif
      select type(new_mechanism)
        type is(wf_mechanism_glass_type)
          if (uninitialized(new_mechanism%specific_surface_area)) then
            option%io_buffer = 'SPECIFIC_SURFACE_AREA must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%k0)) then
            option%io_buffer = 'K0 must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%k_long)) then
            option%io_buffer = 'K_LONG must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%nu)) then
            option%io_buffer = 'NU must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%Ea)) then
            option%io_buffer = 'EA must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%Q)) then
            option%io_buffer = 'Q must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%K)) then
            option%io_buffer = 'K must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%pH)) then
            option%io_buffer = 'PH must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%v)) then
            option%io_buffer = 'V must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block, or choose &
                               &the KIENZLER_DISSOLUTION option.'
            call printErrMsg(option)
          endif
        type is(wf_mechanism_custom_type)
          if (uninitialized(new_mechanism%specific_surface_area) .and. &
              uninitialized(new_mechanism%dissolution_rate) .and. &
              uninitialized(new_mechanism%frac_dissolution_rate)) then
            option%io_buffer = 'FRACTIONAL_DISSOLUTION_RATE or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
          if ( (initialized(new_mechanism%frac_dissolution_rate) .and. &
                initialized(new_mechanism%dissolution_rate)    ) .or. &
               (uninitialized(new_mechanism%frac_dissolution_rate) .and. &
                uninitialized(new_mechanism%dissolution_rate)    ) ) then
            option%io_buffer = 'Either FRACTIONAL_DISSOLUTION_RATE or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block. &
                               &Both types of dissolution rates cannot be &
                               &specified.'
            call printErrMsg(option)
          endif
          if ( (initialized(new_mechanism%specific_surface_area) .and. &
                uninitialized(new_mechanism%dissolution_rate)  ) .or. &
               (uninitialized(new_mechanism%specific_surface_area) .and. &
                initialized(new_mechanism%dissolution_rate)      ) ) then
            option%io_buffer = 'FRACTIONAL_DISSOLUTION_RATE or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
        type is(wf_mechanism_fmdm_type)
          if (uninitialized(new_mechanism%burnup)) then
            option%io_buffer = 'BURNUP must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
          if (uninitialized(new_mechanism%specific_surface_area)) then
            option%io_buffer = 'SPECIFIC_SURFACE_AREA must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
      end select
      if (uninitialized(new_mechanism%matrix_density)) then
        option%io_buffer = 'MATRIX_DENSITY must be specified in ' // &
                           trim(error_string) // ' ' // &
                           trim(new_mechanism%name) // ' block.'
        call printErrMsg(option)
      endif

      if (new_mechanism%canister_degradation_model .and. &
          uninitialized(new_mechanism%canister_material_constant)) then
        option%io_buffer = 'CANISTER_MATERIAL_CONSTANT must be given in the '&
                           // trim(error_string) // ' ' // &
                           trim(new_mechanism%name) // &
                           ', CANISTER_DEGRADATION_MODEL block.'
        call printErrMsg(option)
      endif

      if (.not.associated(new_mechanism%rad_species_list)) then
        option%io_buffer = 'At least one SPECIES must be specified in the ' // &
          trim(error_string) // ' ' // trim(new_mechanism%name) // ' block.'
        call printErrMsg(option)
      endif

      if (.not.associated(this%mechanism_list)) then
        this%mechanism_list => new_mechanism
      else
        cur_mechanism => this%mechanism_list
        do
          if (.not.associated(cur_mechanism)) exit
          if (.not.associated(cur_mechanism%next)) then
            cur_mechanism%next => new_mechanism
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_mechanism => cur_mechanism%next
        enddo
      endif
      nullify(new_mechanism)
  !-------------------------------------    
    case default !(MECHANISM keyword not found)
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMWFReadMechanism

! ************************************************************************** !

subroutine PMWFReadWasteForm(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the waste form 
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016
  !
  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use Dataset_Ascii_class 
  use String_module
  use Units_module
  use Region_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscBool :: added
  character(len=MAXWORDLENGTH) :: word, internal_units
  class(waste_form_base_type), pointer :: new_waste_form, cur_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism

  error_string = trim(error_string) // ',WASTE_FORM'
  found = PETSC_TRUE
  added = PETSC_FALSE

  select case(trim(keyword))
  !-------------------------------------
    case('WASTE_FORM')
      allocate(new_waste_form)
      new_waste_form => WasteFormCreate()
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('EXPOSURE_FACTOR')
            call InputReadDouble(input,option,new_waste_form%exposure_factor)
            call InputErrorMsg(input,option,'exposure factor',error_string)
        !-----------------------------
          case('VOLUME')
            call InputReadDouble(input,option,new_waste_form%volume)
            call InputErrorMsg(input,option,'volume',error_string)
            call InputReadAndConvertUnits(input,new_waste_form%volume, &
                                          'm^3',trim(error_string)//',volume', &
                                          option)
            new_waste_form%init_volume = new_waste_form%volume
        !-----------------------------
          case('COORDINATE')
            call GeometryReadCoordinate(input,option, &
                                        new_waste_form%coordinate,error_string)
        !-----------------------------
          case('REGION')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'region assignment',error_string)
            new_waste_form%region_name = trim(word)
        !-----------------------------
          case('MECHANISM_NAME')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'mechanism assignment',error_string)
            call StringToUpper(word)
            new_waste_form%mech_name = trim(word)
        !-----------------------------
          case('CANISTER_VITALITY_RATE')
            call InputReadDouble(input,option, &
                                 new_waste_form%canister_vitality_rate)
            call InputErrorMsg(input,option,'canister vitality rate', &
                               error_string)
            call InputReadAndConvertUnits(input, &
                        new_waste_form%canister_vitality_rate,'unitless/sec', &
                        trim(error_string)//'canister vitality rate',option)
        !-----------------------------
          case('CANISTER_BREACH_TIME')
            call InputReadDouble(input,option, &
                                 new_waste_form%breach_time)
            call InputErrorMsg(input,option,'CANISTER_BREACH_TIME',error_string)
            call InputReadAndConvertUnits(input,new_waste_form%breach_time, &
                           'sec',trim(error_string)//',CANISTER_BREACH_TIME', &
                           option)
        !-----------------------------    
          case default
            call InputKeywordUnrecognized(word,error_string,option)
        !-----------------------------
        end select
      enddo
      
     ! ----------------- error messaging -------------------------------------
      if (Uninitialized(new_waste_form%volume)) then
        option%io_buffer = 'VOLUME must be specified for all waste forms.'
        call printErrMsg(option)
      endif
      if (Uninitialized(new_waste_form%coordinate%z) .and. &
          (len(trim(new_waste_form%region_name)) == 0)) then
        option%io_buffer = 'Either COORDINATE or REGION must be specified &
                           &for all waste forms.'
        call printErrMsg(option)
      endif
      if (Initialized(new_waste_form%coordinate%z) .and. &
          (len(trim(new_waste_form%region_name)) > 0)) then
        option%io_buffer = 'Either COORDINATE or REGION must be specified &
                           &for all waste forms, but not both.'
        call printErrMsg(option)
      endif
      if (new_waste_form%mech_name == '') then
        option%io_buffer = 'MECHANISM_NAME must be specified for &
                           &all waste forms.'
        call printErrMsg(option)
      endif
      !note: do not throw error if EXPOSURE_FACTOR isn't specified (default = 1)
      
      if (.not.associated(this%waste_form_list)) then
        this%waste_form_list => new_waste_form
      else
        cur_waste_form => this%waste_form_list
        do
          if (.not.associated(cur_waste_form)) exit
          if (.not.associated(cur_waste_form%next)) then
            cur_waste_form%next => new_waste_form
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_waste_form => cur_waste_form%next
        enddo
      endif
      nullify(new_waste_form)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (.not.associated(this%waste_form_list)) then
    option%io_buffer = 'At least one WASTE_FORM must be specified in the &
                       &WASTE_FORM_GENERAL block.'
    call printErrMsg(option)
  endif

end subroutine PMWFReadWasteForm

! ************************************************************************** !

subroutine PMWFAssociateRegion(this,region_list)
  ! 
  ! Associates the waste form to its assigned region via the REGION keyword
  ! or the COORDINATE keyword.
  ! 
  ! Author: Jenn Frederick
  ! Date: 10/24/2016
  !

  use Region_module
  use Option_module
  use String_module
  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module

  implicit none
  
  class(pm_waste_form_type) :: this
  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: cur_region
  type(region_type), pointer :: new_region
  class(waste_form_base_type), pointer :: cur_waste_form
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: word1, word2
  PetscBool :: matched
  PetscReal :: x, y, z
  PetscInt :: i, j, k
  PetscInt :: local_id(1)
  PetscInt :: coordinate_counter
  
  option => this%option
  grid => this%realization%patch%grid
  coordinate_counter = 0
  
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    ! if COORDINATE was given, auto-create a region for it
    if (Initialized(cur_waste_form%coordinate%z)) then
      local_id(1) = -1
      coordinate_counter = coordinate_counter + 1
      x = cur_waste_form%coordinate%x
      y = cur_waste_form%coordinate%y
      z = cur_waste_form%coordinate%z
      select case(grid%itype)
        case(STRUCTURED_GRID)
          call StructGridGetIJKFromCoordinate(grid%structured_grid,x,y,z, &
                                              i,j,k)
          if (i > 0 .and. j > 0 .and. k > 0) then
            local_id(1) = i + (j-1)*grid%structured_grid%nlx + &
                        (k-1)*grid%structured_grid%nlxy
          endif
        case(IMPLICIT_UNSTRUCTURED_GRID)
          call UGridGetCellFromPoint(x,y,z, &
                                     grid%unstructured_grid,option,local_id(1))
        case default
          option%io_buffer = 'Only STRUCTURED_GRID and &
                    &IMPLICIT_UNSTRUCTURED_GRID types supported in PMWasteForm.'
          call printErrMsg(option)
      end select
      ! create the region only if the current process owns the waste form
      if (local_id(1) > 0) then
        new_region => RegionCreate(local_id)
        write(word1,'(i6)') coordinate_counter
        write(word2,'(i6)') option%myrank
        new_region%name = 'WF_COORDINATE_' // trim(adjustl(word1)) // '_p' //  &
                          trim(adjustl(word2))
        cur_waste_form%region => new_region
        allocate(cur_waste_form%scaling_factor(1))
        cur_waste_form%scaling_factor(1) = 1.d0
      endif
    else
      cur_region => region_list%first
      do
        if (.not.associated(cur_region)) exit
        matched = PETSC_FALSE
        if (StringCompare(trim(cur_region%name), &
                          trim(cur_waste_form%region_name))) then
          cur_waste_form%region => cur_region
          matched = PETSC_TRUE
        endif
        if (matched) exit
        cur_region => cur_region%next
      enddo
      if (.not.associated(cur_waste_form%region)) then
        option%io_buffer = 'WASTE_FORM REGION ' // &
                           trim(cur_waste_form%region_name) // ' not found.'
        call printErrMsg(option)
      endif
      !call PMWFSetRegionScaling(this,cur_waste_form)
    endif
    !
    cur_waste_form => cur_waste_form%next
  enddo
  
end subroutine PMWFAssociateRegion
  
! ************************************************************************** !

subroutine PMWFSetRegionScaling(this,waste_form)
  ! 
  ! Calculates and sets the scaling factor vector for each of the waste forms
  ! that have assigned regions. This function is called only if a region was
  ! just associated with it. It assumes the volume of the cells that make up
  ! the region do not change over the course of the simulation.
  ! 
  ! Author: Jenn Frederick
  ! Date: 10/21/2016
  !

  use Material_Aux_class
  use Grid_module

  implicit none
  
  class(pm_waste_form_type) :: this
  class(waste_form_base_type), pointer :: waste_form
  
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  PetscInt :: k, cell_id
  PetscReal :: total_volume_local, total_volume_global
  PetscErrorCode :: ierr
  
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid
  allocate(waste_form%scaling_factor(waste_form%region%num_cells))
  total_volume_local = 0.d0
  total_volume_global = 0.d0
  
  ! scale by cell volume
  do k = 1,waste_form%region%num_cells
    cell_id = grid%nL2G(waste_form%region%cell_ids(k))
    waste_form%scaling_factor(k) = material_auxvars(cell_id)%volume ! [m^3]
    total_volume_local = total_volume_local &
                         + material_auxvars(cell_id)%volume  ! [m^3]
  enddo
  call MPI_Allreduce(total_volume_local,total_volume_global,ONE_INTEGER_MPI, &
              MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
  waste_form%scaling_factor = waste_form%scaling_factor/total_volume_global 
  
end subroutine PMWFSetRegionScaling

! ************************************************************************** !

subroutine PMWFSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Realization_Subsurface_class

  implicit none
  
  class(pm_waste_form_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWFSetRealization

! ************************************************************************** !

subroutine PMWFSetup(this)
  ! 
  ! Associates the waste forms to their regions and sets the waste form id.
  ! Creates an MPI group/communicator for processes that own a waste form.
  ! Throws out waste forms on processes that do not own the waste form region.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
  ! Notes: Updated/modified by Jennifer Frederick, 10/24/2016.

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Option_module
  use Reaction_Aux_module
  use Utility_module, only : GetRndNumFromNormalDist
  use String_module

  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: species_name
  character(len=MAXWORDLENGTH), pointer :: names(:)
  class(waste_form_base_type), pointer :: cur_waste_form, prev_waste_form
  class(waste_form_base_type), pointer :: next_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism
  PetscInt :: waste_form_id
  PetscInt :: i
  PetscBool :: local, found
  PetscErrorCode :: ierr
  PetscMPIInt :: newcomm
  
  option => this%realization%option
  reaction => this%realization%reaction
  
  ! point the waste form region to the desired region 
  call PMWFAssociateRegion(this,this%realization%patch%region_list)
  
  waste_form_id = 0
  nullify(prev_waste_form)
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    waste_form_id = waste_form_id + 1
    local = PETSC_FALSE
   !--------- canister degradation model --------------------
    if (cur_waste_form%mechanism%canister_degradation_model) then
      cur_waste_form%canister_degradation_flag = PETSC_TRUE
      cur_waste_form%canister_vitality = 1.d0
      ! waste form breach time specified:
      if (initialized(cur_waste_form%breach_time) .and. &
          uninitialized(cur_waste_form%canister_vitality_rate)) then
        cur_waste_form%eff_canister_vit_rate = &
          (1.d0/cur_waste_form%breach_time)
      ! distribution for canister degradation rate specified:
      elseif (uninitialized(cur_waste_form%canister_vitality_rate) .and. &
              uninitialized(cur_waste_form%breach_time)) then
        ! call to random number generator must be done while each processor
        ! knows about every other processor's waste forms, otherwise the
        ! memory of the random number generator will not be global
        call GetRndNumFromNormalDist( &
             cur_waste_form%mechanism%vitality_rate_mean, &
             cur_waste_form%mechanism%vitality_rate_stdev,&
             cur_waste_form%canister_vitality_rate)
        if (cur_waste_form%canister_vitality_rate > &
            cur_waste_form%mechanism%vitality_rate_trunc) then
          cur_waste_form%canister_vitality_rate = &
            cur_waste_form%mechanism%vitality_rate_trunc
        endif
        ! Given rates are in units of log-10/yr, so convert to 1/yr:
        cur_waste_form%canister_vitality_rate = &
          10.0**(cur_waste_form%canister_vitality_rate)
        ! Convert rates from 1/yr to internal units of 1/sec
        cur_waste_form%canister_vitality_rate = &
          cur_waste_form%canister_vitality_rate * &
          (1.0/365.0/24.0/3600.0)
      endif
    endif
   !----------------------------------------------------------       
    if (associated(cur_waste_form%region)) then
      if (cur_waste_form%region%num_cells > 0) then
          local = PETSC_TRUE
      endif
    endif
    if (local) then
      ! assign the waste form ID number and MPI communicator
      cur_waste_form%id = waste_form_id
      cur_waste_form%myMPIgroup_id = waste_form_id
      call MPI_Comm_split(option%mycomm,cur_waste_form%myMPIgroup_id, &
                          option%myrank,newcomm,ierr)
    else
      cur_waste_form%id = 0
      cur_waste_form%myMPIgroup_id = 0
      call MPI_Comm_split(option%mycomm,MPI_UNDEFINED,option%myrank, &
                          newcomm,ierr)
    endif
    cur_waste_form%myMPIcomm = newcomm
    if (local) then
      call PMWFSetRegionScaling(this,cur_waste_form)
      prev_waste_form => cur_waste_form
      cur_waste_form => cur_waste_form%next
    else 
      ! remove waste form because it is not local
      next_waste_form => cur_waste_form%next
      if (associated(prev_waste_form)) then
        prev_waste_form%next => next_waste_form
      else
        this%waste_form_list => next_waste_form
      endif
      deallocate(cur_waste_form)
      cur_waste_form => next_waste_form
    endif
  enddo
  
  ! check if the mechanism list includes fmdm or glass mechanisms:
  cur_mechanism => this%mechanism_list
  do
    if (.not.associated(cur_mechanism)) exit
    select type(cur_mechanism)
      ! set up indexing of solute concentrations for fmdm model:
      type is(wf_mechanism_fmdm_type)
        species_name = 'O2(aq)'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iO2) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'HCO3-'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iCO3_2n) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'H2(aq)'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iH2) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
        species_name = 'Fe++'
        cur_mechanism%mapping_fmdm_to_pflotran(cur_mechanism%iFe_2p) = &
          GetPrimarySpeciesIDFromName(species_name,reaction,option)
      type is(wf_mechanism_glass_type)
        if (cur_mechanism%use_pH) then
          if ((associated(this%realization%reaction%species_idx) .and. &
              (this%realization%reaction%species_idx%h_ion_id == 0)) .or. &
              (.not.associated(this%realization%reaction%species_idx))) then
            option%io_buffer = 'pH may not be calculated when H+ is not &
                               &defined as a species and/or the card&
                               &USE_FULL_GEOCHEMISTRY is not specified in the &
                               &CHEMISTRY block - (MECHANISM GLASS).'
            call printErrMsg(option)
          else
            cur_mechanism%h_ion_id = &
                                 this%realization%reaction%species_idx%h_ion_id                                
          endif  
        endif
        if (cur_mechanism%use_Q) then   
          species_name = 'SiO2(aq)'
          if (associated(this%realization%reaction)) then
            ! search through the species names so that the generic error
            ! message from GetPrimarySpeciesIDFromName is not thrown first
            ! when SiO2 is missing
            allocate(names(GetPrimarySpeciesCount(this%realization%reaction)))
            names => GetPrimarySpeciesNames(this%realization%reaction)
            i = 0
            found = PETSC_FALSE
            do while (i < len(names))
              i = i + 1
              if (adjustl(trim(species_name)) == adjustl(trim(names(i)))) then
                cur_mechanism%SiO2_id = &
                                     GetPrimarySpeciesIDFromName(species_name, &
                                     this%realization%reaction,option)
                found = PETSC_TRUE
              endif
              if (found) exit
              if ((.not.found) .and. (i == len(names))) then
                deallocate(names)
                allocate(names(GetSecondarySpeciesCount(this%realization%reaction)))
                names => GetSecondarySpeciesNames(this%realization%reaction)
                i = 0
                do while (i < len(names))
                  i = i + 1
                  if (adjustl(trim(species_name)) == &
                      adjustl(trim(names(i)))) then
                    cur_mechanism%SiO2_id = &
                                   GetSecondarySpeciesIDFromName(species_name, &
                                   this%realization%reaction,option)
                    found = PETSC_TRUE
                  endif
                  if (found) exit
                  if ((.not.found) .and. (i == len(names))) then
                    option%io_buffer = 'Q may not be calculated when SiO2(aq) &
                              &is not defined as a primary or secondary &
                              &species - (MECHANISM GLASS).'
                    call printErrMsg(option)
                  endif
                enddo
              endif
            enddo
            deallocate(names)
          else
            option%io_buffer = 'Q may not be calculated when SiO2(aq) is not &
                               &defined as a species and/or the card &
                               &USE_FULL_GEOCHEMISTRY is not specified in the &
                               &CHEMISTRY block - (MECHANISM GLASS).'
            call printErrMsg(option)
          endif
        endif
    end select
    cur_mechanism => cur_mechanism%next
  enddo
  
  
end subroutine PMWFSetup

! ************************************************************************** !

 subroutine PMWFInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/25/15
  use Reaction_Aux_module
  use Realization_Base_class
  
  implicit none

#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_type) :: this
  
  IS :: is
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: num_waste_form_cells
  PetscInt :: num_species
  PetscInt :: size_of_vec
  PetscInt :: i, j, k
  PetscInt :: data_mediator_species_id
  PetscInt, allocatable :: species_indices_in_residual(:)
  PetscErrorCode :: ierr

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    num_species = cur_waste_form%mechanism%num_species
    allocate(cur_waste_form%instantaneous_mass_rate(num_species))
    allocate(cur_waste_form%cumulative_mass(num_species))
    cur_waste_form%instantaneous_mass_rate = 0.d0
    cur_waste_form%cumulative_mass = 0.d0
    allocate(cur_waste_form%rad_mass_fraction(num_species))
    allocate(cur_waste_form%rad_concentration(num_species))
    allocate(cur_waste_form%inst_release_amount(num_species))
    cur_waste_form%rad_mass_fraction = &
      cur_waste_form%mechanism%rad_species_list%mass_fraction
    cur_waste_form%rad_concentration = 0.d0
    cur_waste_form%inst_release_amount = 0.d0
    do j = 1, num_species
      cur_waste_form%mechanism%rad_species_list(j)%ispecies = &
        GetPrimarySpeciesIDFromName( &
        cur_waste_form%mechanism%rad_species_list(j)%name, &
        this%realization%reaction,this%option)
    enddo
    cur_waste_form => cur_waste_form%next
  enddo
  
  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
  endif
  
  if (.not.this%option%restart_flag .and. this%print_mass_balance) then
    call PMWFOutputHeader(this)
    call PMWFOutput(this)
  endif

  ! set up mass transfer
  call RealizCreateTranMassTransferVec(this%realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(this%realization%tran_data_mediator_list)
  ! create a Vec sized by # waste packages * # waste package cells in region *
  ! # primary dofs influenced by waste package
  ! count of waste form cells
  cur_waste_form => this%waste_form_list
  num_waste_form_cells = 0
  size_of_vec = 0
  do
    if (.not.associated(cur_waste_form)) exit
    size_of_vec = size_of_vec + (cur_waste_form%mechanism%num_species * &
                                 cur_waste_form%region%num_cells)
    num_waste_form_cells = num_waste_form_cells + 1
    cur_waste_form => cur_waste_form%next
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,size_of_vec, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  if (num_waste_form_cells > 0) then
    allocate(species_indices_in_residual(size_of_vec))
    species_indices_in_residual = 0
    cur_waste_form => this%waste_form_list
    i = 0
    do
      if (.not.associated(cur_waste_form)) exit
      do k = 1,cur_waste_form%region%num_cells
        do j = 1,cur_waste_form%mechanism%num_species
          i = i + 1
          species_indices_in_residual(i) = &
              (cur_waste_form%region%cell_ids(k)-1)*this%option%ntrandof + &
              cur_waste_form%mechanism%rad_species_list(j)%ispecies
        enddo
      enddo
      cur_waste_form => cur_waste_form%next
    enddo                             ! zero-based indexing
    species_indices_in_residual(:) = species_indices_in_residual(:) - 1
    ! set to global petsc index
    species_indices_in_residual(:) = species_indices_in_residual(:) + &
      this%realization%patch%grid%global_offset*this%option%ntrandof
  endif
  call ISCreateGeneral(this%option%mycomm,size_of_vec, &
                       species_indices_in_residual, &
                       PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
  if (allocated(species_indices_in_residual)) &
    deallocate(species_indices_in_residual)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_OBJECT, &
                        this%realization%field%tran_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)

  call PMWFSolve(this,0.d0,ierr)
  
end subroutine PMWFInitializeRun

! ************************************************************************** !

subroutine PMWFInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick 03/28/2016

  use Global_Aux_module
  use Material_Aux_class
  use Field_module
  use Option_module
  use Grid_module
  use Patch_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  class(pm_waste_form_type) :: this
  
  class(waste_form_base_type), pointer :: cur_waste_form
  class(wf_mechanism_base_type), pointer :: cwfm
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  PetscReal :: rate
  PetscReal :: dV
  PetscReal :: dt
  PetscReal :: avg_temp
  PetscInt :: i, k, p, g, d, f
  PetscInt :: num_species
  PetscErrorCode :: ierr
  PetscInt :: cell_id, idof
  PetscReal, allocatable :: Coeff(:)
  PetscReal, allocatable :: concentration_old(:)
  PetscReal :: inst_release_molality
  PetscReal, parameter :: conversion = 1.d0/(24.d0*3600.d0)
  PetscReal, pointer :: xx_p(:)

  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  field => this%realization%field
  option => this%option
  grid => this%realization%patch%grid
  dt = option%tran_dt
  
  if (option%print_screen_flag) then
    write(*,'(/,2("=")," WASTE FORM MODEL ",60("="))')
  endif

  ! zero entries from previous time step
  call VecZeroEntries(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  cur_waste_form => this%waste_form_list
  do 
    if (.not.associated(cur_waste_form)) exit
    cwfm => cur_waste_form%mechanism
    num_species = cwfm%num_species
    allocate(Coeff(num_species))
    allocate(concentration_old(num_species))
    ! ------ update mass balances after transport step ---------------------
    select type(cwfm => cur_waste_form%mechanism)
      class is(wf_mechanism_dsnf_type)
        ! note: do nothing here because the cumulative mass update for dsnf
        ! and/or wipp mechanisms has already occured (if breached)
      class default
        cur_waste_form%cumulative_mass = cur_waste_form%cumulative_mass + &
          cur_waste_form%instantaneous_mass_rate*dt
    end select
    ! ------ update matrix volume ------------------------------------------
    select type(cwfm => cur_waste_form%mechanism)
      class is(wf_mechanism_dsnf_type)
        ! note: do nothing here because the volume update for dsnf and/or
        ! wipp mechanisms has already occured (if breached)
      class default
         dV = cur_waste_form%eff_dissolution_rate / &   ! kg-matrix/sec
           cwfm%matrix_density * &                      ! kg-matrix/m^3-matrix
           dt                                           ! sec
         cur_waste_form%volume = cur_waste_form%volume - dV
    end select
    if (cur_waste_form%volume <= 1.d-8) then
      cur_waste_form%volume = 0.d0
    endif
    
    ! ------ get species concentrations from mass fractions ----------------
    do k = 1,num_species
      if (cur_waste_form%volume <= 0.d0) then
        cur_waste_form%rad_concentration(k) = 0.d0
        cur_waste_form%rad_mass_fraction(k) = 0.d0
      else
        cur_waste_form%rad_concentration(k) = &
          cur_waste_form%rad_mass_fraction(k) / &
          cwfm%rad_species_list(k)%formula_weight
      endif
    enddo

    !---------------- vitality degradation function ------------------------
    if (cur_waste_form%canister_degradation_flag .and. &
        (cur_waste_form%canister_vitality > 1.d-3)) then
      if (.not.cur_waste_form%breached .and. &
          initialized(cur_waste_form%breach_time)) then
        ! do not modify eff_canister_vit_rate from what it was set to
        cur_waste_form%eff_canister_vit_rate = &
          cur_waste_form%eff_canister_vit_rate   
      else
        i = 0
        avg_temp = 0.d0
        do while (i < cur_waste_form%region%num_cells)
          i = i + 1
          avg_temp = avg_temp + &  ! Celcius
            global_auxvars(grid%nL2G(cur_waste_form%region%cell_ids(i)))%temp* &
            cur_waste_form%scaling_factor(i)
        enddo
        call MPI_Allreduce(MPI_IN_PLACE,avg_temp,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,cur_waste_form%myMPIcomm,ierr)
        avg_temp = avg_temp+273.15d0   ! Kelvin
        cur_waste_form%eff_canister_vit_rate = &
          cur_waste_form%canister_vitality_rate * &
          exp( cwfm%canister_material_constant * ( (1.d0/333.15d0) - &
          (1.d0/(avg_temp))) )
      endif
      cur_waste_form%canister_vitality = cur_waste_form%canister_vitality &
                                 - (cur_waste_form%eff_canister_vit_rate*dt)
      if (cur_waste_form%canister_vitality <= 1.d-3) then
        cur_waste_form%canister_vitality = 0.d0
        cur_waste_form%eff_canister_vit_rate = 0.d0
        cur_waste_form%canister_vitality_rate = 0.d0
      endif
    endif

    !------- instantaneous release ----------------------------------------- 
    if ((.not.cur_waste_form%breached .and. &
         cur_waste_form%canister_vitality < 1.d-3) .or. &
        (.not.cur_waste_form%breached .and. &
         initialized(cur_waste_form%breach_time) .and. &
         option%time > cur_waste_form%breach_time)) then
      call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
      do k = 1,num_species
        cur_waste_form%inst_release_amount(k) = &
           (cwfm%rad_species_list(k)%inst_release_fraction * &
            cur_waste_form%rad_concentration(k))
        cur_waste_form%rad_concentration(k) = &
           cur_waste_form%rad_concentration(k) - &
           cur_waste_form%inst_release_amount(k)
        ! update mass fractions after instantaneous release
        cur_waste_form%rad_mass_fraction(k) = &
           cur_waste_form%rad_concentration(k) * &
           cwfm%rad_species_list(k)%formula_weight
        ! update transport solution vector with mass injection molality
        ! as an alternative to a source term (issue with tran_dt changing)
        do f = 1, cur_waste_form%region%num_cells
          cell_id = grid%nL2G(cur_waste_form%region%cell_ids(f))
          inst_release_molality = &                    ! [mol-rad/kg-water]
            ! [mol-rad]
            (cur_waste_form%inst_release_amount(k) * & ! [mol-rad/g-matrix]
             cur_waste_form%volume * &                 ! [m^3-matrix]
             cwfm%matrix_density * &                   ! [kg-matrix/m^3-matrix] 
             1.d3) / &                               ! [kg-matrix] -> [g-matrix]
             ! [kg-water]
            (material_auxvars(cell_id)%porosity * &         ! [-]
             global_auxvars(cell_id)%sat(LIQUID_PHASE) * &  ! [-]
             material_auxvars(cell_id)%volume * &           ! [m^3]
             global_auxvars(cell_id)%den_kg(LIQUID_PHASE))  ! [kg/m^3-water]
          idof = cwfm%rad_species_list(k)%ispecies + &
                 ((cur_waste_form%region%cell_ids(f) - 1) * option%ntrandof) 
          xx_p(idof) = xx_p(idof) + & 
                       (inst_release_molality*cur_waste_form%scaling_factor(f))
        enddo
      enddo
      cur_waste_form%breached = PETSC_TRUE 
      cur_waste_form%breach_time = option%time
      call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    endif
    
    ! Save the concentration after inst. release for the decay step
    concentration_old = cur_waste_form%rad_concentration

    if (cur_waste_form%volume >= 0.d0) then
      !------- decay the radionuclide species --------------------------------
      ! FIRST PASS =====================
      do d = 1,num_species
        ! Update the initial value of the species coefficient
        Coeff(d) = cur_waste_form%rad_concentration(d)
        do p = 1,num_species
          ! If the daughter has a parent(s):
          if (d == cwfm%rad_species_list(p)%daugh_id) then
            Coeff(d) = Coeff(d) - &
              (cwfm%rad_species_list(p)%decay_constant * &
               concentration_old(p)) / &
              (cwfm%rad_species_list(d)%decay_constant - &
               cwfm%rad_species_list(p)%decay_constant)
            do g = 1,num_species
              ! If the daughter has a grandparent(s):
              if (p == cwfm%rad_species_list(g)%daugh_id) then
                Coeff(d) = Coeff(d) - &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant))) + &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(p)%decay_constant)))
              endif
            enddo ! grandparent loop
          endif
        enddo ! parent loop
      enddo
      ! SECOND PASS ====================
      do d = 1,num_species
        ! Decay the species
        cur_waste_form%rad_concentration(d) = Coeff(d) * exp(-1.d0 * &
          cwfm%rad_species_list(d)%decay_constant * dt)
        do p = 1,num_species
          ! If the daughter has a parent(s):
          if (d == cwfm%rad_species_list(p)%daugh_id) then
            cur_waste_form%rad_concentration(d) = &
              cur_waste_form%rad_concentration(d) + &
              (((cwfm%rad_species_list(p)%decay_constant* &
                 concentration_old(p)) / &
                (cwfm%rad_species_list(d)%decay_constant - &
                 cwfm%rad_species_list(p)%decay_constant)) * &
               exp(-1.d0 * cwfm%rad_species_list(p)%decay_constant * dt)) 
            do g = 1,num_species
              ! If the daughter has a grandparent(s):
              if (p == cwfm%rad_species_list(g)%daugh_id) then
                cur_waste_form%rad_concentration(d) = &
                  cur_waste_form%rad_concentration(d) - &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)*exp(-1.d0* &
                    cwfm%rad_species_list(p)%decay_constant*dt)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(p)%decay_constant))) + &
                  ((cwfm%rad_species_list(p)%decay_constant* &
                    cwfm%rad_species_list(g)%decay_constant* &
                    concentration_old(g)*exp(-1.d0* &
                    cwfm%rad_species_list(g)%decay_constant*dt)) / &
                  ((cwfm%rad_species_list(p)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)* &
                   (cwfm%rad_species_list(d)%decay_constant - &
                    cwfm%rad_species_list(g)%decay_constant)))
              endif
            enddo ! grandparent loop
          endif
        enddo ! parent loop
      enddo     

      ! ------ update species mass fractions ---------------------------------
      do k = 1,num_species
        cur_waste_form%rad_mass_fraction(k) = &       ! [g-rad/g-wf]
        cur_waste_form%rad_concentration(k) * &       ! [mol-rad/g-wf]
          cur_waste_form%mechanism%rad_species_list(k)%formula_weight
        ! to avoid errors in plotting data when conc is very very low:  
        if (cur_waste_form%rad_mass_fraction(k) <= 1d-40) then
          cur_waste_form%rad_mass_fraction(k) = 0.d0
        endif
      enddo
    endif
    deallocate(concentration_old)
    deallocate(Coeff)
    cur_waste_form => cur_waste_form%next
  enddo

  if (this%print_mass_balance) then
    call PMWFOutput(this)
  endif

end subroutine PMWFInitializeTimestep

! ************************************************************************** !

subroutine PMWFSolve(this,time,ierr)
  !
  ! Updates the source term based on the dissolution model chosen
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Updated/modified by Jenn Frederick 04/2016
  !
  ! Notes: The species loop must be the inner loop, while the grid cell loop
  ! must be the outer loop, in order for the vec_p(i) indexing to work.
  
  use Global_Aux_module
  use Material_Aux_class
  use Grid_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: i, j, k
  PetscInt :: num_species
  PetscInt :: cell_id
  PetscInt :: idof
  PetscReal :: inst_diss_molality
  PetscReal, pointer :: vec_p(:)  
  PetscReal, pointer :: xx_p(:)
  PetscInt :: fmdm_count_global, fmdm_count_local
  character(len=MAXWORDLENGTH) :: word
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type), pointer :: grid
  
  fmdm_count_global = 0
  fmdm_count_local = 0
  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  grid => this%realization%patch%grid

  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%realization%field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  
  cur_waste_form => this%waste_form_list
  i = 0
  do 
    if (.not.associated(cur_waste_form)) exit
    num_species = cur_waste_form%mechanism%num_species  
    !
    if ((cur_waste_form%volume > 0.d0) .and. &
        (cur_waste_form%canister_vitality <= 1.d-40)) then
    !---------------------------------------------------------------------------
      ! calculate the mechanism-specific eff_dissolution_rate [kg-matrix/sec]:
      call cur_waste_form%mechanism%Dissolution(cur_waste_form,this,ierr)
      select type(cwfm => cur_waste_form%mechanism)
      !-----------------------------------------------------------------------
        ! ignore source term if dsnf/wipp type, and directly update the
        ! solution vector instead (see note in WFMech[DSNF/WIPP]Dissolution):
        class is(wf_mechanism_dsnf_type)
          do k = 1,cur_waste_form%region%num_cells
            cell_id = grid%nL2G(cur_waste_form%region%cell_ids(k))
            do j = 1,num_species
              i = i + 1
              cur_waste_form%instantaneous_mass_rate(j) = &
                (cur_waste_form%eff_dissolution_rate / &      ! kg-matrix/sec
                 cur_waste_form%mechanism%rad_species_list(j)%formula_weight * &! kg-rad/kmol-rad
                 cur_waste_form%rad_mass_fraction(j) * &      ! kg-rad/kg-matrix
                 1.d3) 
              inst_diss_molality = &                          ! mol-rad/kg-water
                cur_waste_form%instantaneous_mass_rate(j) * & ! mol-rad/sec
                this%realization%option%tran_dt / &           ! sec
                ! [kg-water]
                (material_auxvars(cell_id)%porosity * &        ! [-]
                 global_auxvars(cell_id)%sat(LIQUID_PHASE) * & ! [-]
                 material_auxvars(cell_id)%volume * &          ! [m^3]
                 global_auxvars(cell_id)%den_kg(LIQUID_PHASE)) ! [kg/m^3-water]
              idof = cwfm%rad_species_list(j)%ispecies + &
                     ((cur_waste_form%region%cell_ids(k) - 1) * &
                      this%option%ntrandof)
              xx_p(idof) = xx_p(idof) + &                     ! mol-rad/kg-water
                           (inst_diss_molality*cur_waste_form%scaling_factor(k))  
              vec_p(i) = 0.d0
              if (k == 1) then
                ! update the cumulative mass now, not at next timestep:
                cur_waste_form%cumulative_mass(j) = &
                  cur_waste_form%cumulative_mass(j) + &          ! mol-rad
                  cur_waste_form%instantaneous_mass_rate(j) * &  ! mol-rad/sec
                  this%realization%option%tran_dt                ! sec
                ! update the volume now, not at next timestep:
                cur_waste_form%volume = 0.d0 
              endif              
            enddo 
          enddo
      !-----------------------------------------------------------------------
        ! for all other waste form types, load the source term, and update
        ! the cumulative mass and volume at next timestep:
        class default
          do k = 1,cur_waste_form%region%num_cells
            do j = 1,num_species
              i = i + 1
              cur_waste_form%instantaneous_mass_rate(j) = &
                (cur_waste_form%eff_dissolution_rate / &      ! kg-matrix/sec
                 cur_waste_form%mechanism%rad_species_list(j)%formula_weight * &! kg-rad/kmol-rad
                 cur_waste_form%rad_mass_fraction(j) * &      ! kg-rad/kg-matrix
                 1.d3)
              vec_p(i) = cur_waste_form%instantaneous_mass_rate(j) * &
                         cur_waste_form%scaling_factor(k)  ! mol/sec * [-]
            enddo 
          enddo
      !-------------------------------------------------------------------------
      end select
      ! count the number of times FMDM was called:
      select type(cwfm => cur_waste_form%mechanism)
        type is(wf_mechanism_fmdm_type)
          fmdm_count_local = fmdm_count_local + 1
      end select
    !---------------------------------------------------------------------------
    else ! (canister not breached, or all waste form has dissolved already)
      vec_p((i+1):(i+num_species*cur_waste_form%region%num_cells)) = 0.d0
      i = i + num_species*cur_waste_form%region%num_cells
      cur_waste_form%eff_dissolution_rate = 0.d0
      cur_waste_form%instantaneous_mass_rate = 0.d0
    endif
    !
    cur_waste_form => cur_waste_form%next
  enddo
 
  ! ideally, this print statement would go inside the dissolution subroutine
  call MPI_Allreduce(fmdm_count_local,fmdm_count_global,ONE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,this%realization%option%mycomm,ierr)
  if ((fmdm_count_global > 0) .and. &
      this%realization%option%print_screen_flag) then
    write(word,'(i5)') fmdm_count_global
    write(*,'(/,2("=")," FMDM ",72("="))')
  ! ** START (this can be removed after FMDM profiling is finished) **
    write(*,'(a)') '== ' // adjustl(trim(word)) // ' calls to FMDM.'
  ! ** END (this can be removed after FMDM profiling is finished) **
  endif
  
  call VecRestoreArrayF90(this%realization%field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMWFSolve

! ************************************************************************** !

subroutine WFMechBaseDissolution(this,waste_form,pm,ierr) 
  ! 
  ! Calculates the waste form dissolution rate; must be extended
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  implicit none
  
  class(wf_mechanism_base_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr

  ! This routine must be extended.
  print *, 'subroutine WFMechBaseDissolution must be extended!'
  stop

end subroutine WFMechBaseDissolution

! ************************************************************************** !

subroutine WFMechGlassDissolution(this,waste_form,pm,ierr) 
  ! 
  ! Calculates the glass waste form dissolution rate
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module  
  use Global_Aux_module
  use String_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(wf_mechanism_glass_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
                                            ! 1/day -> 1/sec
  PetscReal, parameter :: time_conversion = 1.d0/(24.d0*3600.d0)
  PetscReal :: avg_temp
  PetscReal :: avg_pri_molal, avg_sec_molal
  PetscReal :: avg_pri_act_coef, avg_sec_act_coef
  PetscInt :: i

  grid => pm%realization%patch%grid
  global_auxvars => pm%realization%patch%aux%Global%auxvars
  rt_auxvars => pm%realization%patch%aux%RT%auxvars
  
  ierr = 0

  ! Glass dissolution equation: Kienzler et al. (2012) Eq. 6 pg. 17
  ! Kienzler, B., M. Altmaier, et al. (2012) Radionuclide Source Term form
  ! HLW Glass, Spent Nuclear Fuel, and Compacted Hulls and End Pieces
  ! (CSD-C Waste). KIT Scientific Reports 7624. Karlsruhe Institute of
  ! Technology, Baden-Wurttemberg, Germany.
  ! Generalized glass dissolution equation comes from Eq. 2.3.7-6 in
  ! Yucca Mountain Repository SAR, Section 2.3.7, DOE/RW-0573 Rev.0
  
  i = 0
  avg_temp = 0.d0
  do while (i < waste_form%region%num_cells)
    i = i + 1
    avg_temp = avg_temp + &  ! Celcius
               global_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))%temp * &
               waste_form%scaling_factor(i)
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,avg_temp,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
  avg_temp = avg_temp+273.15d0   ! Kelvin
              
  if (this%use_pH) then
    if (this%h_ion_id > 0) then   ! primary species
      i = 0
      avg_pri_molal = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_pri_molal = avg_pri_molal + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        pri_molal(this%h_ion_id)*waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_pri_molal,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      i = 0
      avg_pri_act_coef = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_pri_act_coef = avg_pri_act_coef + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        pri_act_coef(this%h_ion_id)*waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_pri_act_coef,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      this%pH = -log10(avg_pri_molal*avg_pri_act_coef)
    elseif (this%h_ion_id < 0) then   ! secondary species
      i = 0
      avg_sec_molal = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_sec_molal = avg_sec_molal + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        sec_molal(this%h_ion_id)*waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_sec_molal,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      i = 0
      avg_sec_act_coef = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_sec_act_coef = avg_sec_act_coef + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        sec_act_coef(this%h_ion_id)*waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_sec_act_coef,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      this%pH = -log10(avg_sec_molal*avg_sec_act_coef)
    endif
  endif
  
  if (this%use_Q) then
    if (this%SiO2_id > 0) then   ! primary species
      i = 0
      avg_pri_molal = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_pri_molal = avg_pri_molal + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        pri_molal(this%SiO2_id)*waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_pri_molal,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      i = 0
      avg_pri_act_coef = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_pri_act_coef = avg_pri_act_coef + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        pri_act_coef(this%SiO2_id)*waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_pri_act_coef,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      this%Q = avg_pri_molal*avg_pri_act_coef
    elseif (this%SiO2_id < 0) then   ! secondary species
      i = 0
      avg_sec_molal = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_sec_molal = avg_sec_molal + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        sec_molal(abs(this%SiO2_id))* &
                        waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_sec_molal,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      i = 0
      avg_sec_act_coef = 0.d0
      do while (i < waste_form%region%num_cells)
        i = i + 1
        avg_sec_act_coef = avg_sec_act_coef + &
                        rt_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))% &
                        sec_act_coef(abs(this%SiO2_id))* &
                        waste_form%scaling_factor(i)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,avg_sec_act_coef,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
      this%Q = avg_sec_molal*avg_sec_act_coef
    endif
  endif
  
  ! kg-glass/m^2/sec
  this%dissolution_rate = this%k0 * (10.d0**(this%nu*this%pH)) * &
                          exp(-this%Ea/(8.314d0*avg_temp)) * &
                          (1.d0 - (this%Q/this%K)**(1/this%v)) + this%k_long

  ! kg-glass/sec
  waste_form%eff_dissolution_rate = &
    this%dissolution_rate * &          ! kg-glass/m^2/sec
    this%specific_surface_area * &     ! m^2/kg glass
    this%matrix_density * &            ! kg-glass/m^3-glass
    waste_form%volume * &              ! m^3-glass
    waste_form%exposure_factor         ! [-]

end subroutine WFMechGlassDissolution

! ************************************************************************** !

subroutine WFMechDSNFDissolution(this,waste_form,pm,ierr) 
  ! 
  ! Calculates the DSNF waste form dissolution rate
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  implicit none

  class(wf_mechanism_dsnf_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm  
  PetscErrorCode :: ierr
  
  ierr = 0
  
  ! Because the DSNF dissolution rate is instantaneous, the amount of
  ! released isotopes gets updated directly in the solution vector after
  ! this routine is called, within PMWFSolve.
  ! Doing the direct update to the solution vector resolves the potential
  ! error that may occur if the next timestep size is different from the
  ! current timestep size, when the dissolution rate would have been
  ! calculated. This potential error is greatly reduced in magnitude for
  ! the other dissolution models, so we only do the direct update for DSNF.
  
  ! the entire waste form dissolves in the current timestep:
  this%frac_dissolution_rate = 1.d0 / pm%realization%option%tran_dt

  ! kg-matrix/sec
  waste_form%eff_dissolution_rate = &
    this%frac_dissolution_rate * &           ! 1/sec
    this%matrix_density * &                  ! kg matrix/m^3 matrix
    waste_form%volume * &                    ! m^3 matrix
    waste_form%exposure_factor               ! [-]

end subroutine WFMechDSNFDissolution

! ************************************************************************** !

subroutine WFMechWIPPDissolution(this,waste_form,pm,ierr) 
  ! 
  ! Calculates the WIPP waste form dissolution rate
  ! 
  ! Author: Jenn Frederick
  ! Date: 012/8/2016

  implicit none

  class(wf_mechanism_wipp_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm  
  PetscErrorCode :: ierr
  
  ! This subroutine is only a placeholder.
  ! The WIPP waste form mechanism is an extension of the DSNF waste form
  ! mechanism and does not have its own dissolution routine, however, this
  ! placeholder exists in case a different one should be implemented.

end subroutine WFMechWIPPDissolution

! ************************************************************************** !

subroutine WFMechFMDMDissolution(this,waste_form,pm,ierr)
  !
  ! Calculates the FMDM waste form dissolution rate using the FMDM model
  !
  ! Author: Jenn Frederick (with old code by Glenn Hammond)
  ! Date: 05/05/2016

  use Grid_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Option_module

  implicit none
  
 ! FMDM model:
 !=================================================================== 
  interface
    subroutine AMP_step ( burnup, sTme, temperature_C, conc, &
                          initialRun, fuelDisRate, Usource, success )
      real ( kind = 8), intent( in ) :: burnup   
      real ( kind = 8), intent( in ) :: sTme   
      real ( kind = 8), intent( in ) :: temperature_C   
      real ( kind = 8), intent( inout ),  dimension (:,:) :: conc
      logical ( kind = 4), intent( in ) :: initialRun
      ! sum of fluxes of 3 uranium compounds (UO2,2+;UCO3,2+;UO2)
      ! units: g/m^2/yr where g = sum of uranium compound mass
      real ( kind = 8), intent(out) :: fuelDisRate 
      ! flux of just the uranium from the 3 uranium compounds 
      ! units: g/m^2/yr where g = uranium mass
      real ( kind = 8), intent(out) :: Usource
      integer ( kind = 4), intent(out) :: success
    end subroutine
  end interface  
 !===================================================================

  class(wf_mechanism_fmdm_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: i, k
  PetscInt :: icomp_fmdm
  PetscInt :: icomp_pflotran
  PetscInt :: ghosted_id
  
 ! FMDM model: 
 !=======================================================
  integer ( kind = 4) :: success
  logical ( kind = 4) :: initialRun
  PetscReal :: time
  PetscReal :: Usource
  PetscReal :: avg_temp
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(option_type), pointer :: option
 !========================================================
  
  grid => pm%realization%patch%grid
  rt_auxvars => pm%realization%patch%aux%RT%auxvars
  global_auxvars => pm%realization%patch%aux%Global%auxvars
  option => pm%realization%option

  ierr = 0
  
  do k = 1,waste_form%region%num_cells
    ghosted_id = grid%nL2G(waste_form%region%cell_ids(k))
    ! overwrite the components in mapping_pflotran array
    do i = 1, size(this%mapping_fmdm)
      icomp_fmdm = this%mapping_fmdm(i)
      icomp_pflotran = this%mapping_fmdm_to_pflotran(icomp_fmdm)
      this%concentration(icomp_fmdm,:) = &
        rt_auxvars(ghosted_id)%total(icomp_pflotran,LIQUID_PHASE)
    enddo
  enddo
  
  ! convert total component concentration from mol/L to mol/m3 (*1.d3)
  this%concentration = this%concentration*1.d3
  
  if (waste_form%volume /= waste_form%init_volume) then
    initialRun = PETSC_FALSE
  else
    initialRun = PETSC_TRUE
  endif 
  
#ifdef FMDM_MODEL  
 ! FMDM model calculates this%dissolution_rate and Usource [g/m^2/yr]:
 !====================================================================
  time = option%time
  
  i = 0
  avg_temp = 0.d0
  do while (i < waste_form%region%num_cells)
    i = i + 1
    avg_temp = avg_temp + &  ! Celcius
               global_auxvars(grid%nL2G(waste_form%region%cell_ids(i)))%temp * &
               waste_form%scaling_factor(i)
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,avg_temp,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,waste_form%myMPIcomm,ierr)
  call AMP_step(this%burnup, time, avg_temp, this%concentration, &
                initialRun, this%dissolution_rate, Usource, success) 
  write(*,*) this%dissolution_rate
  ! convert total component concentration from mol/m3 back to mol/L (/1.d3)
  this%concentration = this%concentration/1.d3
  ! convert this%dissolution_rate from fmdm to pflotran units:
  ! g/m^2/yr => kg/m^2/sec
  this%dissolution_rate = this%dissolution_rate / (1000.0*24.0*3600.0*365)
  Usource = Usource / (1000.0*24.0*3600.0*365)
 !====================================================================
#else
  ! if no FMDM model, use the burnup as this%dissolution_rate:
  ! if no FMDM model, the units of burnup should already be kg-matrix/m^2/sec:
  success = 1
  this%dissolution_rate = this%burnup
  Usource = this%burnup
#endif

  if (success == 0) then
    ierr = 1
    return
  endif
  
  !==================
  this%frac_dissolution_rate = &    ! 1/sec
    this%dissolution_rate * &       ! kg-matrix/m^2/sec
    this%specific_surface_area      ! m^2/kg-matrix
  !==================
  
  ! kg-matrix / sec
  waste_form%eff_dissolution_rate = &
     this%dissolution_rate * &         ! kg-matrix/m^2/sec
     this%specific_surface_area * &    ! m^2/kg-matrix
     this%matrix_density * &           ! kg-matrix/m^3-matrix
     waste_form%volume * &             ! m^3-matrix
     waste_form%exposure_factor        ! [-]
  
end subroutine WFMechFMDMDissolution

! ************************************************************************** !

subroutine WFMechCustomDissolution(this,waste_form,pm,ierr) 
  ! 
  ! Calculates the "custom" waste form dissolution rate
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module  
  use Global_Aux_module

  implicit none

  class(wf_mechanism_custom_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm
  PetscErrorCode :: ierr

  ! Note: Units for dissolution rates have already been converted to
  ! internal units within the PMWFRead routine.
  
  ierr = 0

  if (uninitialized(this%frac_dissolution_rate)) then
    ! kg-matrix / sec
    waste_form%eff_dissolution_rate = &
       this%dissolution_rate * &         ! kg-matrix/m^2/sec
       this%specific_surface_area * &    ! m^2/kg-matrix
       this%matrix_density * &           ! kg-matrix/m^3-matrix
       waste_form%volume * &             ! m^3-matrix
       waste_form%exposure_factor        ! [-]
  else
    ! kg-matrix / sec
    waste_form%eff_dissolution_rate = &
       this%frac_dissolution_rate * &     ! [-]/sec
       this%matrix_density * &            ! kg-matrix/m^3-matrix
       waste_form%volume * &              ! m^3 matrix
       waste_form%exposure_factor         ! [-]
  endif

end subroutine WFMechCustomDissolution

! ************************************************************************** !

subroutine PMWFFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_type) :: this
  
end subroutine PMWFFinalizeTimestep

! ************************************************************************** !

subroutine PMWFUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_type) :: this
  
  PetscErrorCode :: ierr
  
  ! update glass mass here?

end subroutine PMWFUpdateSolution

! ************************************************************************** !

recursive subroutine PMWFFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  
  implicit none
  
  class(pm_waste_form_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMWFFinalizeRun

! ************************************************************************** !

subroutine PMWFOutput(this)
  ! 
  ! Sets up output for a waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module
  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(waste_form_base_type), pointer :: cur_waste_form
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: i
  
  if (.not.associated(this%waste_form_list)) return
  
100 format(100es18.8)
101 format(1I6.1)

  option => this%realization%option
  output_option => this%realization%output_option
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  
  fid = 86
  filename = PMWFOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the reactive transport step
  write(fid,100,advance="no") option%time / output_option%tconv
  
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit  
    write(fid,101,advance="no") cur_waste_form%id
    do i = 1, cur_waste_form%mechanism%num_species
      write(fid,100,advance="no") cur_waste_form%cumulative_mass(i), &
                                  cur_waste_form%instantaneous_mass_rate(i), &
                                  cur_waste_form%rad_mass_fraction(i)
    enddo
    write(fid,100,advance="no") cur_waste_form%eff_dissolution_rate, &
                                cur_waste_form%volume, &
                                cur_waste_form%eff_canister_vit_rate, &
                                cur_waste_form%canister_vitality*100.0
    cur_waste_form => cur_waste_form%next
  enddo
  close(fid)
  
end subroutine PMWFOutput

! ************************************************************************** !

function PMWFOutputFilename(option)
  ! 
  ! Generates filename for waste form output
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: PMWFOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  PMWFOutputFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       '-wf_mass-' // trim(adjustl(word)) // '.dat'
  
end function PMWFOutputFilename  

! ************************************************************************** !

subroutine PMWFOutputHeader(this)
  ! 
  ! Writes header for waste form output file
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Output_Aux_module
  use Grid_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(output_option_type), pointer :: output_option
  type(grid_type), pointer :: grid
  class(waste_form_base_type), pointer :: cur_waste_form
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: icolumn, i
  
  if (.not.associated(this%waste_form_list)) return
  
  output_option => this%realization%output_option
  grid => this%realization%patch%grid
  
  fid = 86
  filename = PMWFOutputFilename(this%option)
  open(unit=fid,file=filename,action="write",status="replace")  
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif 
  
  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    if (initialized(cur_waste_form%coordinate%z)) then
      ! cell natural id
      write(cell_string,*) grid%nG2A(grid%nL2G(cur_waste_form%region%cell_ids(1)))
      cell_string = ' (' // trim(adjustl(cell_string)) // ')'
      ! coordinate of waste form
      x_string = BestFloat(cur_waste_form%coordinate%x,1.d4,1.d-2)
      y_string = BestFloat(cur_waste_form%coordinate%y,1.d4,1.d-2)
      z_string = BestFloat(cur_waste_form%coordinate%z,1.d4,1.d-2)
      cell_string = trim(cell_string) // &
               ' (' // trim(adjustl(x_string)) // &
               ' ' // trim(adjustl(y_string)) // &
               ' ' // trim(adjustl(z_string)) // ')'
    else
      cell_string = trim(cur_waste_form%region_name)
    endif
    variable_string = 'WF ID#'
    units_string = ''
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    do i = 1, cur_waste_form%mechanism%num_species
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Cum. Mass Flux'
      ! cumulative
      units_string = 'mol'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Inst. Mass Flux'
      ! instantaneous
      units_string = 'mol/s' !// trim(adjustl(output_option%tunit))
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)       
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Mass Frac.'
      units_string = 'g-rad/g-matrix' 
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
    enddo
    variable_string = 'WF Dissolution Rate'
    units_string = 'kg/s'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'WF Volume'
    units_string = 'm^3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'WF Vitality Degradation Rate'
    units_string = '1/yr'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'WF Canister Vitality'
    units_string = '%' 
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)

    cur_waste_form => cur_waste_form%next
  enddo
  
  close(fid)
  
end subroutine PMWFOutputHeader

! ***************************************************************************** !

subroutine PMWFCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
  use Option_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_type) :: this
  PetscViewer :: viewer
  
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: maximum_waste_form_id
  PetscInt :: local_waste_form_count
  PetscInt :: temp_int
  
  Vec :: local, global
  PetscErrorCode :: ierr
  
  this%option%io_buffer = 'PMWFCheckpoint not implemented.'
  call printErrMsg(this%option)
  
  ! calculate maximum waste form id
  maximum_waste_form_id = 0
  local_waste_form_count = 0
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    local_waste_form_count = local_waste_form_count + 1
    maximum_waste_form_id = max(maximum_waste_form_id,cur_waste_form%id)
    cur_waste_form => cur_waste_form%next
  enddo
  call MPI_Allreduce(maximum_waste_form_id,temp_int,ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,this%option%mycomm,ierr)
!  call VecCreateMPI(this%option%mycomm,local_waste_form_count,PETSC_DETERMINE,ierr)
  
                     
end subroutine PMWFCheckpoint

! ***************************************************************************** !

subroutine PMWFRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_type) :: this
  PetscViewer :: viewer
  
  this%option%io_buffer = 'PMWFRestart not implemented.'
  call printErrMsg(this%option)
!  call RestartFlowProcessModel(viewer,this%realization)
!  call this%UpdateAuxVars()
!  call this%UpdateSolution()
  
end subroutine PMWFRestart

! ************************************************************************** !

subroutine PMWFInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_waste_form_type) :: this

  character(len=MAXWORDLENGTH) :: word
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: id
  PetscInt :: k

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMWFInputRecord

! ************************************************************************** !

subroutine PMWFStrip(this)
  ! 
  ! Strips the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick, 03/28/2016

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_waste_form_type) :: this
  
  class(waste_form_base_type), pointer :: cur_waste_form, prev_waste_form

  nullify(this%realization)
  nullify(this%data_mediator)

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    prev_waste_form => cur_waste_form
    cur_waste_form => cur_waste_form%next
    call DeallocateArray(prev_waste_form%rad_mass_fraction)
    call DeallocateArray(prev_waste_form%rad_concentration)
    call DeallocateArray(prev_waste_form%inst_release_amount)
    call DeallocateArray(prev_waste_form%instantaneous_mass_rate)
    call DeallocateArray(prev_waste_form%cumulative_mass)
    call DeallocateArray(prev_waste_form%scaling_factor)
    nullify(prev_waste_form%mechanism)
    nullify(prev_waste_form%region)
    deallocate(prev_waste_form)
    nullify(prev_waste_form)
  enddo
  nullify(this%waste_form_list)
  call PMWFMechanismStrip(this)

end subroutine PMWFStrip

! ************************************************************************** !

subroutine PMWFMechanismStrip(this)
  ! 
  ! Strips the waste form mechanisms in the waste form process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  class(pm_waste_form_type) :: this
  
  class(wf_mechanism_base_type), pointer :: cur_mechanism, prev_mechanism

  cur_mechanism => this%mechanism_list
  do
    if (.not.associated(cur_mechanism)) exit
    prev_mechanism => cur_mechanism
    cur_mechanism => cur_mechanism%next
    deallocate(prev_mechanism%rad_species_list)
    nullify(prev_mechanism%rad_species_list)
    select type(prev_mechanism)
      type is(wf_mechanism_fmdm_type)
        call DeallocateArray(prev_mechanism%concentration)
        call DeallocateArray(prev_mechanism%mapping_fmdm)
        call DeallocateArray(prev_mechanism%mapping_fmdm_to_pflotran)
    end select
    deallocate(prev_mechanism)
    nullify(prev_mechanism)
  enddo
  nullify(this%mechanism_list)

end subroutine PMWFMechanismStrip
  
! ************************************************************************** !

subroutine PMWFDestroy(this)
  ! 
  ! Destroys the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_type) :: this
  
  call PMWFStrip(this)
  
end subroutine PMWFDestroy

! ************************************************************************** !
  
end module PM_Waste_Form_class
