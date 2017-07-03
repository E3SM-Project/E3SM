module PM_UFD_Decay_class

  use PM_Base_class
  use Realization_Subsurface_class
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  
  type, public, extends(pm_base_type) :: pm_ufd_decay_type
    class(realization_subsurface_type), pointer :: realization
    PetscInt, pointer :: element_isotopes(:,:)
    PetscInt, pointer :: isotope_to_primary_species(:)
    PetscInt, pointer :: isotope_to_mineral(:)
    PetscReal, pointer :: isotope_decay_rate(:)
    PetscInt, pointer :: isotope_daughters(:,:)
    PetscReal, pointer :: isotope_daughter_stoich(:,:)
    PetscInt, pointer :: isotope_parents(:,:)
    PetscReal, pointer :: element_solubility(:)
    PetscReal, pointer :: element_Kd(:,:)
    PetscInt :: num_elements
    PetscInt :: num_isotopes
    PetscBool :: implicit_solution
    character(len=MAXWORDLENGTH), pointer :: element_name(:)
    character(len=MAXWORDLENGTH), pointer :: isotope_name(:)
    type(isotope_type), pointer :: isotope_list
    type(element_type), pointer :: element_list
  contains
!geh: commented out subroutines can only be called externally
    procedure, public :: Setup => PMUFDDecayInit
    procedure, public :: Read => PMUFDDecayRead
    procedure, public :: PMUFDDecaySetRealization
    procedure, public :: InitializeRun => PMUFDDecayInitializeRun
!!    procedure, public :: FinalizeRun => PMUFDDecayFinalizeRun
    procedure, public :: InitializeTimestep => PMUFDDecayInitializeTimestep
    procedure, public :: FinalizeTimestep => PMUFDDecayFinalizeTimestep
!    procedure, public :: PreSolve => PMUFDDecayPreSolve
    procedure, public :: Solve => PMUFDDecaySolve
!    procedure, public :: PostSolve => PMUFDDecayPostSolve
!    procedure, public :: AcceptSolution => PMUFDDecayAcceptSolution
!    procedure, public :: TimeCut => PMUFDDecayTimeCut
!    procedure, public :: UpdateSolution => PMUFDDecayUpdateSolution
!    procedure, public :: UpdateAuxVars => PMUFDDecayUpdateAuxVars
!    procedure, public :: Checkpoint => PMUFDDecayCheckpoint    
!    procedure, public :: Restart => PMUFDDecayRestart  
    procedure, public :: InputRecord => PMUFDDecayInputRecord
    procedure, public :: Destroy => PMUFDDecayDestroy
  end type pm_ufd_decay_type
  
  type :: isotope_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: element
    PetscInt :: iisotope
    PetscInt :: ielement
    PetscReal :: decay_rate
    type(daughter_type), pointer :: daughter_list
    type(isotope_type), pointer :: next
  end type isotope_type
  
  type :: daughter_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: stoichiometry
    type(daughter_type), pointer :: next
  end type
  
  type :: element_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: ielement
    PetscReal :: solubility
    PetscReal, pointer :: Kd(:)
    character(len=MAXWORDLENGTH), pointer :: Kd_material_name(:)
    type(element_type), pointer :: next
  end type
  
  public :: PMUFDDecayCreate, &
            PMUFDDecayInit !, &
!            PMUFDDecayInitializeTimestepA, &
!            PMUFDDecayInitializeTimestepB, &
!            PMUFDDecayInitializeRun, &
!            PMUFDDecayUpdateSolution, &
!            PMUFDDecayUpdatePropertiesNI, &
!            PMUFDDecayTimeCut, &
!            PMUFDDecayDestroy
  
contains

! ************************************************************************** !

function PMUFDDecayCreate()
  ! 
  ! Creates the UFD decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type), pointer :: PMUFDDecayCreate
  
  allocate(PMUFDDecayCreate)
  call PMBaseInit(PMUFDDecayCreate)

  PMUFDDecayCreate%num_isotopes = 0
  PMUFDDecayCreate%num_elements = 0
  PMUFDDecayCreate%implicit_solution = PETSC_FALSE
  nullify(PMUFDDecayCreate%realization)
  nullify(PMUFDDecayCreate%element_isotopes)
  nullify(PMUFDDecayCreate%isotope_to_primary_species)
  nullify(PMUFDDecayCreate%isotope_to_mineral)
  nullify(PMUFDDecayCreate%isotope_decay_rate)
  nullify(PMUFDDecayCreate%isotope_daughters)
  nullify(PMUFDDecayCreate%isotope_daughter_stoich)
  nullify(PMUFDDecayCreate%isotope_parents)
  nullify(PMUFDDecayCreate%element_solubility)
  nullify(PMUFDDecayCreate%element_Kd)
  nullify(PMUFDDecayCreate%element_name)
  nullify(PMUFDDecayCreate%isotope_name)
  nullify(PMUFDDecayCreate%isotope_list)
  nullify(PMUFDDecayCreate%element_list)

  call PMBaseInit(PMUFDDecayCreate)

end function PMUFDDecayCreate

! ************************************************************************** !

subroutine PMUFDDecayRead(this,input)
  ! 
  ! Reads input file parameters associated with the ufd decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(isotope_type), pointer :: isotope, prev_isotope
  type(daughter_type), pointer :: daughter, prev_daughter
  type(element_type), pointer :: element, prev_element
  PetscInt :: i
  PetscInt, parameter :: MAX_KD_SIZE = 100
  character(len=MAXWORDLENGTH) :: Kd_material_name(MAX_KD_SIZE)
  PetscReal :: Kd(MAX_KD_SIZE)
  PetscReal :: tempreal

  option => this%option
  
  option%io_buffer = 'pflotran card:: UFD_Decay'
  call printMsg(option)
  
  input%ierr = 0
  nullify(prev_isotope)
  nullify(prev_element)
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    select case(trim(word))
      case('ELEMENT')
        error_string = 'UFD Decay, Element'
        element => ElementCreate()
        call InputReadWord(input,option,element%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
            case('SOLUBILITY')
              call InputReadDouble(input,option,element%solubility)
              call InputErrorMsg(input,option,'solubility',error_string)
            case('KD')
              i = 0
              Kd(:) = UNINITIALIZED_DOUBLE
              Kd_material_name(:) = ''
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                i = i + 1
                if (i > MAX_KD_SIZE) then
                  write(word,*) i-1
                  option%io_buffer = 'Kd array in PMUFDDecayRead() must be &
                    &allocated larger than ' // trim(adjustl(word)) // '.'
                  call printErrMsg(option)
                endif
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'Kd material name', &
                                   error_string)
                Kd_material_name(i) = word
                call InputReadDouble(input,option,Kd(i))
                call InputErrorMsg(input,option,'Kd',error_string)
              enddo
              if (i == 0) then
                option%io_buffer = 'No KD/material combinations specified &
                  &under ' // trim(error_string) // '.'
                call printErrMsg(option)
              endif
              allocate(element%Kd(i))
              element%Kd = Kd(1:i)
              allocate(element%Kd_material_name(i))
              element%Kd_material_name = Kd_material_name(1:i)
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          end select
        enddo
        if (associated(prev_element)) then
          prev_element%next => element
        else
          this%element_list => element
        endif
        prev_element => element
        nullify(element)
      case('ISOTOPE')
        nullify(prev_daughter)
        error_string = 'UFD Decay, Isotope'
        isotope => IsotopeCreate()
        call InputReadWord(input,option,isotope%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name',error_string)
        error_string = 'UFD Decay, Isotope, ' // trim(isotope%name)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword',error_string)
          call StringToUpper(word)
          select case(trim(word))
            case('ELEMENT')
              call InputReadWord(input,option,isotope%element,PETSC_TRUE)
              call InputErrorMsg(input,option,'element name',error_string)
            case('DECAY_RATE')
              call InputReadDouble(input,option,isotope%decay_rate)
              call InputErrorMsg(input,option,'decay rate',error_string)
              call InputReadAndConvertUnits(input,isotope%decay_rate,'1/sec', &
                                            trim(error_string)//',decay rate', &
                                            option)
            case('HALF_LIFE')
              call InputReadDouble(input,option,tempreal)
              call InputErrorMsg(input,option,'half life',error_string)
              call InputReadAndConvertUnits(input,tempreal,'sec', &
                                            trim(error_string)//',half life', &
                                            option)
              ! convert half life to rate constant
              isotope%decay_rate = -1.d0*log(0.5d0)/tempreal
            case('DAUGHTER')
              daughter => IsotopeDaughterCreate()
              call InputReadWord(input,option,daughter%name,PETSC_TRUE)
              call InputErrorMsg(input,option,'daughter name',error_string)
              call InputReadDouble(input,option,daughter%stoichiometry)
              call InputErrorMsg(input,option,'daughter stoichiometry', &
                                 error_string)
              if (associated(prev_daughter)) then
                prev_daughter%next => daughter
              else
                isotope%daughter_list => daughter
              endif
              prev_daughter => daughter
              nullify(daughter)
            case default
              call InputKeywordUnrecognized(word,error_string,option)
          end select
        enddo
        if (associated(prev_isotope)) then
          prev_isotope%next => isotope
        else
          this%isotope_list => isotope
        endif
        prev_isotope => isotope
        nullify(isotope)
      case('IMPLICIT_SOLUTION')
        this%implicit_solution = PETSC_TRUE
      case default
        error_string = 'UFD Decay'
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMUFDDecayRead


! ************************************************************************** !

function ElementCreate()
  ! 
  ! Creates isotope object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/20/15
  ! 
  implicit none

  type(element_type), pointer :: ElementCreate
  
  allocate(ElementCreate)
  ElementCreate%name = ''
  ElementCreate%ielement = UNINITIALIZED_INTEGER
  ElementCreate%solubility = UNINITIALIZED_DOUBLE
  nullify(ElementCreate%Kd)
  nullify(ElementCreate%Kd_material_name)
  nullify(ElementCreate%next)
  
end function ElementCreate

! ************************************************************************** !

function IsotopeCreate()
  ! 
  ! Creates isotope object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/20/15
  ! 
  implicit none

  type(isotope_type), pointer :: IsotopeCreate
  
  allocate(IsotopeCreate)
  IsotopeCreate%name = ''
  IsotopeCreate%element = ''
  IsotopeCreate%iisotope = UNINITIALIZED_INTEGER
  IsotopeCreate%ielement = UNINITIALIZED_INTEGER
  IsotopeCreate%decay_rate = UNINITIALIZED_DOUBLE
  nullify(IsotopeCreate%daughter_list)
  nullify(IsotopeCreate%next)
  
end function IsotopeCreate

! ************************************************************************** !

function IsotopeDaughterCreate()
  ! 
  ! Creates daughter object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/20/15
  ! 
  implicit none

  type(daughter_type), pointer :: IsotopeDaughterCreate
  
  allocate(IsotopeDaughterCreate)
  IsotopeDaughterCreate%name = ''
  IsotopeDaughterCreate%stoichiometry = UNINITIALIZED_DOUBLE
  nullify(IsotopeDaughterCreate%next)

end function IsotopeDaughterCreate

! ************************************************************************** !

subroutine PMUFDDecayInit(this)
  ! 
  ! Initializes variables associated with the UFD decay process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Option_module
  use String_module
  use Grid_module
  use Reaction_Aux_module
  use Reaction_Mineral_Aux_module
  use Reactive_Transport_Aux_module
  use Material_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(grid_type), pointer :: grid
  type(isotope_type), pointer :: isotope, isotope2
  type(daughter_type), pointer :: daughter
  type(element_type), pointer :: element
  PetscInt, allocatable :: num_isotopes_per_element(:)
  character(len=MAXWORDLENGTH) :: word
  type(material_property_ptr_type), pointer :: material_property_array(:)
  type(material_property_type), pointer :: material_property
  
  PetscInt :: icount
  PetscInt :: ghosted_id
  PetscInt :: max_daughters_per_isotope
  PetscInt :: max_parents_per_isotope
  PetscBool :: found
  PetscInt :: iisotope, ielement
  PetscInt :: g, ig, p, ip, d, id
  
  option => this%realization%option
  grid => this%realization%patch%grid
  reaction => this%realization%reaction
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  material_property_array => this%realization%patch%material_property_array
  
  do ghosted_id = 1, grid%ngmax
    allocate(rt_auxvars(ghosted_id)%total_sorb_eq(reaction%naqcomp))
    rt_auxvars(ghosted_id)%total_sorb_eq = 0.d0
  enddo
  
  max_daughters_per_isotope = 0
  max_parents_per_isotope = 0

  ! sum the number of isotopes, elements, max number of isotopes per element
  
  element => this%element_list
  icount = 0
  do
    if (.not.associated(element)) exit
    icount = icount + 1
    element%ielement = icount
    element => element%next
  enddo
  this%num_elements = icount
  allocate(num_isotopes_per_element(this%num_elements))
  num_isotopes_per_element = 0
  allocate(this%element_solubility(this%num_elements))
  this%element_solubility = 0.d0
  allocate(this%element_name(this%num_elements))
  this%element_name = ''
  allocate(this%element_Kd(this%num_elements,size(material_property_array)))
  this%element_Kd = UNINITIALIZED_DOUBLE
  element => this%element_list
  do
    if (.not.associated(element)) exit
    if (element%solubility > 0.d0) then
      this%element_solubility(element%ielement) = element%solubility
    else
      option%io_buffer = 'Element "' // trim(element%name) // '" in &
        &UFD_DECAY block must have a solubility greater than zero.'
      call printErrMsg(option)
    endif
    this%element_name(element%ielement) = element%name
    if (.not.associated(element%Kd)) then
      write(word,*) size(material_property_array)
      option%io_buffer = trim(adjustl(word)) // ' Kds must be defined for &
        &element "' // trim(element%name) // '", one for each &
        &MATERIAL_PROPERTY in the format "<string> <double>".'
      call printErrMsg(option)
    endif
    if (size(element%Kd) /= size(material_property_array)) then
      write(word,*) size(element%Kd)
      option%io_buffer = 'Incorrect number of Kds (' // &
        trim(adjustl(word)) // ') specified for number of materials ('
      write(word,*) size(material_property_array)
      option%io_buffer = trim(option%io_buffer) // &
        trim(adjustl(word)) // ') for UFD Decay element "' // &
        trim(element%name) // '".'
      call printErrMsg(option)
    endif
    do icount = 1, size(element%Kd_material_name)
      material_property => &
        MaterialPropGetPtrFromArray(element%Kd_material_name(icount), &
                                    material_property_array)
      this%element_Kd(element%ielement,material_property%internal_id) = &
        element%Kd(icount)
    enddo
    do icount = 1, size(material_property_array)
      if (UnInitialized(this%element_Kd(element%ielement,icount))) then
        option%io_buffer = 'Uninitialized KD in UFD Decay element "' // &
          trim(element%name) // '" for material "' // &
          trim(material_property_array(icount)%ptr%name) // '".'
        call printErrMsg(option)
      endif
    enddo
    element => element%next
  enddo  
  
  this%num_isotopes = 0
  isotope => this%isotope_list
  do
    if (.not.associated(isotope)) exit
    this%num_isotopes = this%num_isotopes + 1
    isotope%iisotope = this%num_isotopes
    found = PETSC_FALSE
    element => this%element_list
    do
      if (.not.associated(element)) exit
      if (StringCompare(isotope%element,element%name)) then
        found = PETSC_TRUE
        isotope%ielement = element%ielement
        num_isotopes_per_element(element%ielement) = &
          num_isotopes_per_element(element%ielement) + 1
        exit
      endif
      element => element%next
    enddo
    if (.not.found) then
      option%io_buffer = 'Element "' // trim(isotope%element) // &
        '" of isotope "' // trim(isotope%name) // &
        '" not found among list of elements.'
      call printErrMsg(option)
    endif
    daughter => isotope%daughter_list
    icount = 0
    do
      if (.not.associated(daughter)) exit
      icount = icount + 1
      daughter => daughter%next
    enddo
    max_daughters_per_isotope = max(max_daughters_per_isotope,icount)
    isotope => isotope%next
  enddo
  
  allocate(this%isotope_name(this%num_isotopes))
  this%isotope_name = ''
  allocate(this%isotope_to_primary_species(this%num_isotopes))
  this%isotope_to_primary_species = UNINITIALIZED_INTEGER
  allocate(this%isotope_to_mineral(this%num_isotopes))
  this%isotope_to_mineral = UNINITIALIZED_INTEGER
  allocate(this%element_isotopes(0:maxval(num_isotopes_per_element), &
           this%num_elements))
  this%element_isotopes = UNINITIALIZED_INTEGER
  this%element_isotopes(0,:) = 0
  allocate(this%isotope_decay_rate(this%num_isotopes))
  this%isotope_decay_rate = UNINITIALIZED_DOUBLE
  allocate(this%isotope_daughters(0:max_daughters_per_isotope, &
                                  this%num_isotopes))
  this%isotope_daughters = UNINITIALIZED_INTEGER
  this%isotope_daughters(0,:) = 0
  allocate(this%isotope_daughter_stoich(max_daughters_per_isotope, &
                                        this%num_isotopes))
  this%isotope_daughter_stoich(:,:) = UNINITIALIZED_DOUBLE
  
  isotope => this%isotope_list
  do
    if (.not.associated(isotope)) exit
    found = PETSC_FALSE
    this%isotope_name(isotope%iisotope) = isotope%name
    this%isotope_to_primary_species(isotope%iisotope) = &
      GetPrimarySpeciesIDFromName(isotope%name,reaction,option)
    word = isotope%name
    word = trim(word) // '(s)'
    this%isotope_to_mineral(isotope%iisotope) = &
      GetKineticMineralIDFromName(word,reaction%mineral,option)
    this%element_isotopes(0,isotope%ielement) = &
      this%element_isotopes(0,isotope%ielement) + 1
    this%element_isotopes(this%element_isotopes(0,isotope%ielement), &
                          isotope%ielement) = isotope%iisotope
    this%isotope_decay_rate(isotope%iisotope) = isotope%decay_rate
    daughter => isotope%daughter_list
    icount = 0
    do
      if (.not.associated(daughter)) exit
      icount = icount + 1
      isotope2 => this%isotope_list
      found = PETSC_FALSE
      do
        if (.not.associated(isotope2)) exit
        if (StringCompare(daughter%name,isotope2%name)) then
          found = PETSC_TRUE
          this%isotope_daughters(icount,isotope%iisotope) = isotope2%iisotope
          this%isotope_daughters(0,isotope%iisotope) = icount
          this%isotope_daughter_stoich(icount,isotope%iisotope) = &
            daughter%stoichiometry
          exit
        endif
        isotope2 => isotope2%next
      enddo
      if (.not.found) then
        option%io_buffer = 'Daughter "' // trim(daughter%name) // &
                           '" not found among isotope list.'
        call printErrMsg(option)
      endif
      daughter => daughter%next
    enddo
    isotope => isotope%next
  enddo
  
  allocate(this%isotope_parents(1,this%num_isotopes))
  this%isotope_parents = 0
  do iisotope = 1, this%num_isotopes
    do icount = 1, this%isotope_daughters(0,iisotope)
      this%isotope_parents(1,this%isotope_daughters(icount,iisotope)) = &
        this%isotope_parents(1,this%isotope_daughters(icount,iisotope)) + 1
    enddo
  enddo
  max_parents_per_isotope = maxval(this%isotope_parents)
  deallocate(this%isotope_parents)
  allocate(this%isotope_parents(0:max_parents_per_isotope,this%num_isotopes))
  this%isotope_parents = 0
  do iisotope = 1, this%num_isotopes
    do icount = 1, this%isotope_daughters(0,iisotope)
      this%isotope_parents(0,this%isotope_daughters(icount,iisotope)) = &
        this%isotope_parents(0,this%isotope_daughters(icount,iisotope)) + 1
      this%isotope_parents(this%isotope_parents(0, &
                                    this%isotope_daughters(icount,iisotope)), &
                           this%isotope_daughters(icount,iisotope)) = iisotope
    enddo
  enddo
  ! error checking
  do ielement = 1, this%num_elements
    ! nothing yet
  enddo
  do iisotope = 1, this%num_isotopes
    ! ensure that a decay rate is defined
    if (Uninitialized(this%isotope_decay_rate(iisotope))) then
      option%io_buffer = 'A decay rate must be defined for isotope "' // &
        trim(this%isotope_name(iisotope)) // '".'
      call printErrMsg(option)
    endif
    ! ensure that a stoichiometry is defined for all daughters
    do d = 1, this%isotope_daughters(0,iisotope)
      id = this%isotope_daughters(d,iisotope)
      if (Uninitialized(this%isotope_daughter_stoich(d,iisotope))) then
        option%io_buffer = 'A stoichiomtry must be defined for isotope ' // &
          trim(this%isotope_name(iisotope)) // "'s daughter " // '"' // &
          trim(this%isotope_name(id)) // '".'
        call printErrMsg(option)
      endif
    enddo
    ! ensure that a daughter is not the same as a parent or grandparent. this 
    ! will produce NaNs through a divide by zero.
    do p = 1, this%isotope_parents(0,iisotope)
      ip = this%isotope_parents(p,iisotope)
      if (ip == iisotope) then
        option%io_buffer = 'PM UFD_DECAY isotope "' // &
          trim(this%isotope_name(iisotope)) // &
          '" is the same as its parent.'
        call printErrMsg(option)
      endif
      do g = 1, this%isotope_parents(0,ip)
        ig = this%isotope_parents(g,ip)
        if (ig == iisotope) then
          option%io_buffer = 'PM UFD_DECAY isotope "' // &
            trim(this%isotope_name(iisotope)) // &
            '" is the same as its grandparent.'
          call printErrMsg(option)
        endif
      enddo
    enddo
  enddo
  
#if 0  
  this%num_elements = 3
  max_num_isotopes_per_element = 3
  this%num_isotopes = 6
  this%isotope_to_primary_species = [1, 2, 3, 4, 5, 6]
  this%element_isotopes = 0
  this%element_isotopes(0:1,1) = [1,1]
  this%element_isotopes(0:3,2) = [3,2,3,4]
  this%element_isotopes(0:2,3) = [2,5,6]
  this%isotope_decay_rate = 0.d0
  this%isotope_decay_rate(:) = [1.29d-15, 5.08d-11, 1.03d-14, 1.38d-13, 2.78d-12,2.78d-13]
  allocate(this%isotope_daughters(0:max_daughters,this%num_isotopes))
  this%isotope_daughters = 0
  this%isotope_daughters(:,:) = 0
  this%isotope_daughters(0,1) = 1  ! 1 -> 3
  this%isotope_daughters(1,1) = 3
  this%isotope_daughters(0,3) = 1   ! 3 -> 5
  this%isotope_daughters(1,3) = 5
  allocate(this%element_solubility(this%num_elements))
  this%element_solubility = 0.d0
#endif

end subroutine PMUFDDecayInit

! ************************************************************************** !

subroutine PMUFDDecaySetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Realization_Subsurface_class

  implicit none
  
  class(pm_ufd_decay_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMUFDDecaySetRealization

! ************************************************************************** !

recursive subroutine PMUFDDecayInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Patch_module
  use Grid_module
  use Reactive_Transport_Aux_module
  
  implicit none

  class(pm_ufd_decay_type) :: this
  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  
  PetscReal :: kd_kgw_m3b
  PetscInt :: local_id, ghosted_id
  PetscInt :: iele, iiso, ipri, i, imat
  
  patch => this%realization%patch
  grid => patch%grid
  rt_auxvars => patch%aux%RT%auxvars
  
  ! set initial sorbed concentration in equilibrium with aqueous phase
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    imat = patch%imat(ghosted_id) 
    if (imat <= 0) cycle
    do iele = 1, this%num_elements
      kd_kgw_m3b = this%element_Kd(iele,imat)
      do i = 1, this%element_isotopes(0,iele)
        iiso = this%element_isotopes(i,iele)
        ipri = this%isotope_to_primary_species(iiso)
        rt_auxvars(ghosted_id)%total_sorb_eq(ipri) = &
          rt_auxvars(ghosted_id)%pri_molal(ipri) * kd_kgw_m3b
      enddo
    enddo      
  enddo

  if (maxval(this%isotope_daughter_stoich) > 1.d0) then
    this%option%io_buffer = 'Daughter stoichiometries have not been set up &
      &in pm_ufd_decay.F90.'
    call printErrMsg(this%option)
  endif

end subroutine PMUFDDecayInitializeRun

! ************************************************************************** !

subroutine PMUFDDecayInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Global_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," USED FUEL DISPOSITION DECAY MODEL ",43("="))')
  endif

end subroutine PMUFDDecayInitializeTimestep

! ************************************************************************** !

subroutine PMUFDDecayPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  use Grid_module
  use Global_Aux_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  rt_auxvars => this%realization%patch%aux%RT%auxvars
  
end subroutine PMUFDDecayPreSolve

! ************************************************************************** !

subroutine PMUFDDecaySolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  !
  use Option_module
  use Reaction_Aux_module
  use Patch_module
  use Grid_module
  use Field_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Utility_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_ufd_decay_type) :: this
  
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iele, i, p, g, ip, ig, iiso, ipri, imnrl, imat
  PetscReal :: dt
  PetscReal :: vol, por, sat, den_w_kg, vps
  PetscReal :: conc_iso_aq0, conc_iso_sorb0, conc_iso_ppt0
  PetscReal :: conc_ele_aq1, conc_ele_sorb1, conc_ele_ppt1
  PetscReal :: mass_iso_aq0, mass_iso_sorb0, mass_iso_ppt0
  PetscReal :: mass_ele_aq1, mass_ele_sorb1, mass_ele_ppt1, mass_cmp_tot1
  PetscReal :: mass_iso_tot0(this%num_isotopes), mass_iso_tot_star(this%num_isotopes)
  PetscReal :: mass_iso_tot1(this%num_isotopes), delta_mass_iso_tot(this%num_isotopes)
  PetscReal :: mass_ele_tot1
  PetscReal :: coeff(this%num_isotopes)
  PetscReal :: mass_old(this%num_isotopes)
  PetscReal :: mol_fraction_iso(this%num_isotopes)
  PetscReal :: kd_kgw_m3b
  PetscBool :: above_solubility
  PetscReal, pointer :: xx_p(:)
  PetscReal :: norm

  PetscReal :: residual(this%num_isotopes)
  PetscReal :: solution(this%num_isotopes)
  PetscReal :: rhs(this%num_isotopes)
  PetscInt :: indices(this%num_isotopes)
  PetscReal :: Jacobian(this%num_isotopes,this%num_isotopes)
  PetscReal :: rate, rate_constant, stoich, one_over_dt
  PetscReal, parameter :: tolerance = 1.d-12
  PetscInt :: idaughter
  PetscInt :: it

  ierr = 0
  
  option => this%realization%option
  reaction => this%realization%reaction
  patch => this%realization%patch
  field => this%realization%field
  grid => patch%grid
  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  dt = option%tran_dt
  one_over_dt = 1.d0 / dt
  call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    vol = material_auxvars(ghosted_id)%volume
    den_w_kg = global_auxvars(ghosted_id)%den_kg(1)
    por = material_auxvars(ghosted_id)%porosity
    sat = global_auxvars(ghosted_id)%sat(1)
    vps = vol * por * sat ! m^3 water
    
    ! sum up mass of each isotope across phases and decay
    do iele = 1, this%num_elements
      do i = 1, this%element_isotopes(0,iele)
        iiso = this%element_isotopes(i,iele)
        ipri = this%isotope_to_primary_species(iiso)
        imnrl = this%isotope_to_mineral(iiso)
        ! # indicated time level (0 = prev time level, 1 = new time level) 
        conc_iso_aq0 = xx_p((local_id-1)*reaction%ncomp+ipri) * &
                       den_w_kg / 1000.d0  ! mol/L
        !conc_iso_aq0 = rt_auxvars(ghosted_id)%total(ipri,1) ! mol/L
        conc_iso_sorb0 = rt_auxvars(ghosted_id)%total_sorb_eq(ipri) ! mol/m^3 bulk
        conc_iso_ppt0 = rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) ! m^3 mnrl/m^3 bulk
        mass_iso_aq0 = conc_iso_aq0*vps*1.d3 ! mol/L * m^3 water * 1000 L /m^3 = mol
        mass_iso_sorb0 = conc_iso_sorb0 * vol ! mol/m^3 bulk * m^3 bulk = mol
        mass_iso_ppt0 = conc_iso_ppt0 * vol / &  ! m^3 mnrl/m^3 bulk * m^3 bulk / (m^3 mnrl/mol mnrl) = mol
                        reaction%mineral%kinmnrl_molar_vol(imnrl)
        mass_iso_tot0(iiso) = mass_iso_aq0 + mass_iso_sorb0 + mass_iso_ppt0
      enddo
    enddo 
    
    ! save the mass from the previous time step:
    mass_old(:) = mass_iso_tot0(:)

    if (.not.this%implicit_solution) then

      ! FIRST PASS decay ==============================================
      do i = 1,this%num_isotopes
        ! update the initial value of the isotope coefficient:
        coeff(i) = mass_old(i)
        ! loop through the isotope's parents:
        do p = 1,this%isotope_parents(0,i)
          ip = this%isotope_parents(p,i)
          coeff(i) = coeff(i) - (this%isotope_decay_rate(ip) * mass_old(ip)) / &
            (this%isotope_decay_rate(i) - this%isotope_decay_rate(ip))
          ! loop through the isotope's parent's parents:
          do g = 1,this%isotope_parents(0,ip)
            ig = this%isotope_parents(g,ip)
            coeff(i) = coeff(i) - &
              ((this%isotope_decay_rate(ip) * this%isotope_decay_rate(ig) * &        
              mass_old(ig)) / ((this%isotope_decay_rate(ip) - &
              this%isotope_decay_rate(ig)) * (this%isotope_decay_rate(i) - &
              this%isotope_decay_rate(ig)))) + ((this%isotope_decay_rate(ip) * &
              this%isotope_decay_rate(ig) * mass_old(ig)) / &
              ((this%isotope_decay_rate(ip) - this%isotope_decay_rate(ig)) * &
              (this%isotope_decay_rate(i) - this%isotope_decay_rate(ip))))
          enddo ! grandparent loop
        enddo ! parent loop
      enddo ! isotope loop
      ! SECOND PASS decay =============================================
      do i = 1,this%num_isotopes
        ! decay the isotope species:
        mass_iso_tot1(i) = coeff(i)*exp(-1.d0*this%isotope_decay_rate(i)*dt)
        ! loop through the isotope's parents:
        do p = 1,this%isotope_parents(0,i)
          ip = this%isotope_parents(p,i)
          mass_iso_tot1(i) = mass_iso_tot1(i) + &
                (((this%isotope_decay_rate(ip) * mass_old(ip)) / &
                (this%isotope_decay_rate(i) - this%isotope_decay_rate(ip))) * &
                 exp(-1.d0 * this%isotope_decay_rate(ip) * dt))
          ! loop through the isotope's parent's parents:
          do g = 1,this%isotope_parents(0,ip)
            ig = this%isotope_parents(g,ip)
            mass_iso_tot1(i) = mass_iso_tot1(i) - &
              ((this%isotope_decay_rate(ip) * this%isotope_decay_rate(ig) * &
              mass_old(ig) * exp(-1.d0 * this%isotope_decay_rate(ip) * dt)) / &
              ((this%isotope_decay_rate(ip) - this%isotope_decay_rate(ig)) * &
              (this%isotope_decay_rate(i) - this%isotope_decay_rate(ip)))) + &
              ((this%isotope_decay_rate(ip) * this%isotope_decay_rate(ig) * &
              mass_old(ig) * exp(-1.d0 * this%isotope_decay_rate(ig) * dt)) / &
            ((this%isotope_decay_rate(ip) - this%isotope_decay_rate(ig)) * &
            (this%isotope_decay_rate(i) - this%isotope_decay_rate(ig))))
          enddo ! grandparent loop
        enddo ! parent loop
      enddo ! isotope loop

    else
      ! implicit solution approach
      residual = 1.d0
      solution = mass_iso_tot0
      it = 0
      do ! nonlinear loop
        if (dot_product(residual,residual) < tolerance) exit
        it = it + 1
        residual = 0.d0
        Jacobian = 0.d0
        do iiso = 1, this%num_isotopes
          rate_constant = this%isotope_decay_rate(iiso)
          ! accumulation term
          residual(iiso) = residual(iiso) + &
                           (solution(iiso) - mass_iso_tot0(iiso)) * one_over_dt
          Jacobian(iiso,iiso) = Jacobian(iiso,iiso) + &
                                one_over_dt
          rate = rate_constant * solution(iiso)
          residual(iiso) = residual(iiso) + rate
          Jacobian(iiso,iiso) = Jacobian(iiso,iiso) + rate_constant
          do i = 1, this%isotope_daughters(0,iiso)
            idaughter = this%isotope_daughters(i,iiso)
            stoich = this%isotope_daughter_stoich(i,iiso)
            residual(idaughter) = residual(idaughter) - rate * stoich
            Jacobian(idaughter,iiso) = Jacobian(idaughter,iiso) - &
                                       rate_constant * stoich
          enddo
        enddo
        ! scale Jacobian
        do iiso = 1, this%num_isotopes
          norm = max(1.d0,maxval(abs(Jacobian(iiso,:))))
          norm = 1.d0/norm
          rhs(iiso) = residual(iiso)*norm
          Jacobian(iiso,:) = Jacobian(iiso,:)*norm
        enddo 
        ! log formulation for derivatives
        do iiso = 1, this%num_isotopes
          Jacobian(:,iiso) = Jacobian(:,iiso)*solution(iiso)
        enddo
        call ludcmp(Jacobian,this%num_isotopes,indices,i)
        call lubksb(Jacobian,this%num_isotopes,indices,rhs)
        rhs = dsign(1.d0,rhs)*min(dabs(rhs),10.d0)
        solution = solution*exp(-rhs)
      enddo
      mass_iso_tot1 = solution 
    endif

    mass_iso_tot1 = max(mass_iso_tot1,1.d-90)

    do iele = 1, this%num_elements
      ! calculate mole fractions
      mass_ele_tot1 = 0.d0
      mol_fraction_iso = 0.d0
      do i = 1, this%element_isotopes(0,iele)
        iiso = this%element_isotopes(i,iele)
        mass_ele_tot1 = mass_ele_tot1 + mass_iso_tot1(iiso)
      enddo
      do i = 1, this%element_isotopes(0,iele)
        iiso = this%element_isotopes(i,iele)
        mol_fraction_iso(i) = mass_iso_tot1(iiso) / mass_ele_tot1
      enddo

      ! split mass between phases
      kd_kgw_m3b = this%element_Kd(iele,imat)
      conc_ele_aq1 = mass_ele_tot1 / (1.d0+kd_kgw_m3b/(den_w_kg*por*sat)) / &
                         (vps*1.d3)
      above_solubility = conc_ele_aq1 > this%element_solubility(iele)
      if (above_solubility) then
        conc_ele_aq1 = this%element_solubility(iele)
      endif
      ! assume identical sorption for all isotopes in element
      conc_ele_sorb1 = conc_ele_aq1 / den_w_kg * 1.d3 * kd_kgw_m3b
      mass_ele_aq1 = conc_ele_aq1*vps*1.d3
      mass_ele_sorb1 = conc_ele_sorb1 * vol
      ! roundoff can erroneously result in precipitate. this conditional avoids
      ! such an issue
      if (above_solubility) then
        mass_ele_ppt1 = max(mass_ele_tot1 - mass_ele_aq1 - mass_ele_sorb1,0.d0)
      else
        mass_ele_ppt1 = 0.d0
      endif
      conc_ele_ppt1 = mass_ele_ppt1 * &
                      reaction%mineral%kinmnrl_molar_vol(imnrl) / vol
      ! store mass in data structures
      do i = 1, this%element_isotopes(0,iele)
        iiso = this%element_isotopes(i,iele)
        ipri = this%isotope_to_primary_species(iiso)
        imnrl = this%isotope_to_mineral(iiso)
        rt_auxvars(ghosted_id)%total(ipri,1) = &
          conc_ele_aq1 * mol_fraction_iso(i)
        rt_auxvars(ghosted_id)%pri_molal(ipri) = &
          conc_ele_aq1 / den_w_kg * 1.d3 * mol_fraction_iso(i)
        rt_auxvars(ghosted_id)%total_sorb_eq(ipri) = &
          conc_ele_sorb1 * mol_fraction_iso(i)
        rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) = &
          conc_ele_ppt1 * mol_fraction_iso(i)
        ! need to copy primary molalities back into transport solution Vec
        xx_p((local_id-1)*reaction%ncomp+ipri) = &
          rt_auxvars(ghosted_id)%pri_molal(ipri)
      enddo
    enddo      
  enddo

  call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  if (reaction%use_log_formulation) then
    call VecCopy(field%tran_xx,field%tran_log_xx,ierr);CHKERRQ(ierr)
    call VecLog(field%tran_log_xx,ierr);CHKERRQ(ierr)
  endif  
!  call DiscretizationGlobalToLocal(this%realization%discretization, &
!                                   field%tran_xx,field%tran_xx_loc,NTRANDOF)
  
end subroutine PMUFDDecaySolve

! ************************************************************************** !

subroutine PMUFDDecayPostSolve(this)
  ! 
  ! PMUFDDecayUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  ! 
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayPostSolve

! ************************************************************************** !

function PMUFDDecayAcceptSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscBool :: PMUFDDecayAcceptSolution
  
  ! do nothing
  PMUFDDecayAcceptSolution = PETSC_TRUE
  
end function PMUFDDecayAcceptSolution

! ************************************************************************** !

subroutine PMUFDDecayUpdatePropertiesTS(this)
  ! 
  ! Updates parameters/properties at each Newton iteration
  !
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayUpdatePropertiesTS

! ************************************************************************** !

subroutine PMUFDDecayTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine PMUFDDecayTimeCut

! ************************************************************************** !

subroutine PMUFDDecayFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
end subroutine PMUFDDecayFinalizeTimestep

! ************************************************************************** !

subroutine PMUFDDecayUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  PetscErrorCode :: ierr

end subroutine PMUFDDecayUpdateSolution  

! ************************************************************************** !

subroutine PMUFDDecayUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
  
  class(pm_ufd_decay_type) :: this

  this%option%io_buffer = 'PMUFDDecayUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMUFDDecayUpdateAuxVars   

! ************************************************************************** !

subroutine PMUFDDecayCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_ufd_decay_type) :: this
  PetscViewer :: viewer
  
end subroutine PMUFDDecayCheckpoint

! ************************************************************************** !

subroutine PMUFDDecayRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_ufd_decay_type) :: this
  PetscViewer :: viewer
  
end subroutine PMUFDDecayRestart

! ************************************************************************** !

recursive subroutine PMUFDDecayFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  
  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMUFDDecayFinalizeRun

! ************************************************************************** !

subroutine PMUFDDecayInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  use Material_module
  
  implicit none
  
  class(pm_ufd_decay_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id
  PetscInt :: iele
  PetscInt :: iiso
  PetscInt :: i
  PetscInt :: iparent, idaughter
  type(material_property_ptr_type), pointer :: material_property_array(:)

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  material_property_array => this%realization%patch%material_property_array

  do iele = 1, this%num_elements
    write(id,'(2x,"Element: ",a)') this%element_name(iele)
    write(id,'(4x,"Solubility:",es13.5)') this%element_solubility(iele)
    write(id,'(4x,"KDs")')
    do i = 1, size(this%element_Kd,2)
      write(id,'(6x,a32,es13.5)') material_property_array(i)%ptr%name, &
        this%element_Kd(iele,i)
    enddo 
    write(id,'(4x,"Isotopes")')
    do i = 1, this%element_isotopes(0,iele)
      iiso = this%element_isotopes(i,iele)
      write(id,'(6x,a)') this%isotope_name(iiso)
    enddo
  enddo

  do iiso = 1, this%num_isotopes
    write(id,'(2x,"Isotope: ",a)') this%isotope_name(iiso)
    write(id,'(4x,"Primary Species: ",a)') &
      this%realization%reaction%primary_species_names( &
        this%isotope_to_primary_species(iiso))
    write(id,'(4x,"Decay Rate:",es13.5)') this%isotope_decay_rate(iiso)
    write(id,'(4x,"Parent(s)")')
    if (this%isotope_parents(0,iiso) > 0) then
      do i = 1, this%isotope_parents(0,iiso)
        iparent = this%isotope_parents(i,iiso)
        write(id,'(6x,a)') this%isotope_name(iparent)
      enddo
    else
        write(id,'(6x,"None")')
    endif
    write(id,'(4x,"Daughter(s), stoichiometry")')
    if (this%isotope_daughters(0,iiso) > 0) then
      do i = 1, this%isotope_daughters(0,iiso)
        idaughter = this%isotope_daughters(i,iiso)
        write(id,'(6x,a32,es13.5)') this%isotope_name(idaughter), &
          this%isotope_daughter_stoich(i,iiso)
      enddo
    else
        write(id,'(6x,"None")')
    endif
  enddo

end subroutine PMUFDDecayInputRecord

! ************************************************************************** !

subroutine PMUFDDecayDestroy(this)
  ! 
  ! Destroys Subsurface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/24/15
  use Utility_module, only : DeallocateArray
  use Option_module

  implicit none
  
  class(pm_ufd_decay_type) :: this
  
  type(element_type), pointer :: cur_element, prev_element
  type(isotope_type), pointer :: cur_isotope, prev_isotope
  type(daughter_type), pointer :: cur_daughter, prev_daughter
    
  call DeallocateArray(this%element_isotopes)
  call DeallocateArray(this%isotope_to_primary_species)
  call DeallocateArray(this%isotope_to_mineral)
  call DeallocateArray(this%isotope_decay_rate)
  call DeallocateArray(this%isotope_daughters)
  call DeallocateArray(this%isotope_daughter_stoich)
  call DeallocateArray(this%isotope_parents)
  call DeallocateArray(this%isotope_name)
  call DeallocateArray(this%element_solubility)
  call DeallocateArray(this%element_Kd)
  call DeallocateArray(this%element_name)
  
  cur_isotope => this%isotope_list
  do
    if (.not.associated(cur_isotope)) exit
    cur_daughter => cur_isotope%daughter_list
    do
      if (.not.associated(cur_daughter)) exit
      prev_daughter => cur_daughter
      cur_daughter => cur_daughter%next
      nullify(prev_daughter%next)
      deallocate(prev_daughter)
      nullify(prev_daughter)
    enddo
    prev_isotope => cur_isotope
    cur_isotope => cur_isotope%next
    nullify(prev_isotope%next)
    deallocate(prev_isotope)
    nullify(prev_isotope)
  enddo
  
  cur_element => this%element_list
  do
    if (.not.associated(cur_element)) exit
    prev_element => cur_element
    cur_element => cur_element%next
    call DeallocateArray(prev_element%Kd)
    call DeallocateArray(prev_element%Kd_material_name)
    nullify(prev_element%next)
    deallocate(prev_element)
    nullify(prev_element)
  enddo
  
end subroutine PMUFDDecayDestroy
  
end module PM_UFD_Decay_class
