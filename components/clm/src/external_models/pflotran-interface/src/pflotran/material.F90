module Material_module
 
#include "petsc/finclude/petscsys.h"
   use petscsys
  use Dataset_Base_class

  use PFLOTRAN_Constants_module
  use Material_Aux_class
  use Utility_module, only : Equal
  
  implicit none

  private


  type, public :: material_property_type
    PetscInt :: external_id
    PetscInt :: internal_id
    PetscBool :: active
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: permeability(3,3)
    PetscBool :: isotropic_permeability
    PetscReal :: vertical_anisotropy_ratio ! (vertical / horizontal)
    PetscReal :: permeability_scaling_factor
!    character(len=MAXWORDLENGTH) :: permeability_dataset_name
    class(dataset_base_type), pointer :: permeability_dataset
    class(dataset_base_type), pointer :: permeability_dataset_y
    class(dataset_base_type), pointer :: permeability_dataset_z
    PetscReal :: porosity
!    character(len=MAXWORDLENGTH) :: porosity_dataset_name
    class(dataset_base_type), pointer :: porosity_dataset
    class(dataset_base_type), pointer :: tortuosity_dataset
    PetscReal :: tortuosity
    PetscBool :: tortuosity_function_of_porosity
    PetscInt :: saturation_function_id
    character(len=MAXWORDLENGTH) :: saturation_function_name
    PetscReal :: rock_density ! kg/m^3
    PetscReal :: specific_heat ! J/kg-K
    PetscReal :: thermal_conductivity_dry
    PetscReal :: thermal_conductivity_wet
    PetscReal :: alpha    ! conductivity saturation relation exponent

    character(len=MAXWORDLENGTH) :: soil_compressibility_function
    PetscReal :: soil_compressibility
    PetscReal :: soil_reference_pressure
    PetscBool :: soil_reference_pressure_initial
    class(dataset_base_type), pointer :: soil_reference_pressure_dataset
!    character(len=MAXWORDLENGTH) :: compressibility_dataset_name
    class(dataset_base_type), pointer :: compressibility_dataset

    ! ice properties
    PetscReal :: thermal_conductivity_frozen
    PetscReal :: alpha_fr

    PetscReal :: thermal_expansitivity   
    PetscReal :: dispersivity(3)
    PetscReal :: tortuosity_pwr
    PetscReal :: tortuosity_func_porosity_pwr
    PetscReal :: min_pressure
    PetscReal :: max_pressure
    PetscReal :: max_permfactor
    !geh: minral surface area power functions must be defined on a per
    !     mineral basis, look in reaction_aux.F90
    !PetscReal :: mnrl_surf_area_volfrac_pwr
    !PetscReal :: mnrl_surf_area_porosity_pwr
    PetscReal :: permeability_pwr
    PetscReal :: permeability_crit_por
    PetscReal :: permeability_min_scale_fac
    character(len=MAXWORDLENGTH) :: secondary_continuum_name
    PetscReal :: secondary_continuum_length
    PetscReal :: secondary_continuum_matrix_block_size
    PetscReal :: secondary_continuum_fracture_spacing
    PetscReal :: secondary_continuum_radius
    PetscReal :: secondary_continuum_area
    PetscInt :: secondary_continuum_ncells
    PetscReal :: secondary_continuum_epsilon
    PetscReal :: secondary_continuum_aperture
    PetscReal :: secondary_continuum_init_temp
    PetscReal :: secondary_continuum_init_conc
    PetscReal :: secondary_continuum_porosity
    PetscReal :: secondary_continuum_diff_coeff
    PetscReal :: secondary_continuum_mnrl_volfrac
    PetscReal :: secondary_continuum_mnrl_area 
    PetscBool :: secondary_continuum_log_spacing
    PetscReal :: secondary_continuum_outer_spacing
    PetscReal :: secondary_continuum_area_scaling
    type(material_property_type), pointer :: next
  end type material_property_type
  
  type, public :: material_property_ptr_type
    type(material_property_type), pointer :: ptr
  end type material_property_ptr_type
  
  public :: MaterialPropertyCreate, &
            MaterialPropertyDestroy, &
            MaterialPropertyAddToList, &
            MaterialPropGetPtrFromList, &
            MaterialPropGetPtrFromArray, &
            MaterialPropConvertListToArray, &
            MaterialAnisotropyExists, &
            MaterialSetAuxVarScalar, &
            MaterialSetAuxVarVecLoc, &
            MaterialGetAuxVarVecLoc, &
            MaterialAuxVarCommunicate, &
            MaterialPropertyRead, &
            MaterialInitAuxIndices, &
            MaterialAssignPropertyToAux, &
            MaterialSetup, &
            MaterialUpdateAuxVars, &
            MaterialStoreAuxVars, &
            MaterialWeightAuxVars, &
            MaterialGetMaxExternalID, &
            MaterialCreateIntToExtMapping, &
            MaterialCreateExtToIntMapping, &
            MaterialApplyMapping, &
            MaterialPropInputRecord
  
contains

! ************************************************************************** !

function MaterialPropertyCreate()
  ! 
  ! Creates a material property
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 
  
  implicit none

  type(material_property_type), pointer :: MaterialPropertyCreate
  
  type(material_property_type), pointer :: material_property
  
  allocate(material_property)
  material_property%external_id = 0
  material_property%internal_id = 0
  material_property%active = PETSC_TRUE
  material_property%name = ''
  ! initialize to UNINITIALIZED_DOUBLE to catch bugs
  material_property%permeability = UNINITIALIZED_DOUBLE
  material_property%isotropic_permeability = PETSC_TRUE
  material_property%vertical_anisotropy_ratio = UNINITIALIZED_DOUBLE
  material_property%permeability_scaling_factor = 0.d0
  material_property%permeability_pwr = 1.d0
  material_property%permeability_crit_por = 0.d0
  material_property%permeability_min_scale_fac = 1.d0
!  material_property%permeability_dataset_name = ''
  nullify(material_property%permeability_dataset)
  nullify(material_property%permeability_dataset_y)
  nullify(material_property%permeability_dataset_z)
  ! initialize to UNINITIALIZED_DOUBLE to catch bugs
  material_property%porosity = UNINITIALIZED_DOUBLE
!  material_property%porosity_dataset_name = ''
  nullify(material_property%porosity_dataset)
  nullify(material_property%tortuosity_dataset)
  material_property%tortuosity_function_of_porosity = PETSC_FALSE
  material_property%tortuosity = 1.d0
  material_property%tortuosity_pwr = 0.d0
  material_property%tortuosity_func_porosity_pwr = UNINITIALIZED_DOUBLE
  material_property%saturation_function_id = 0
  material_property%saturation_function_name = ''
  material_property%rock_density = UNINITIALIZED_DOUBLE
  material_property%specific_heat = UNINITIALIZED_DOUBLE
  material_property%thermal_conductivity_dry = UNINITIALIZED_DOUBLE
  material_property%thermal_conductivity_wet = UNINITIALIZED_DOUBLE
  material_property%alpha = 0.45d0

  material_property%soil_compressibility_function = ''
  material_property%soil_compressibility = UNINITIALIZED_DOUBLE
  material_property%soil_reference_pressure = UNINITIALIZED_DOUBLE
  material_property%soil_reference_pressure_initial = PETSC_FALSE
  nullify(material_property%soil_reference_pressure_dataset)
!  material_property%compressibility_dataset_name = ''
  nullify(material_property%compressibility_dataset)

  material_property%thermal_conductivity_frozen = 0.d0
  material_property%alpha_fr = 0.95d0

  material_property%thermal_expansitivity = 0.d0  
  material_property%dispersivity = 0.d0
  material_property%min_pressure = 0.d0
  material_property%max_pressure = 1.d6
  material_property%max_permfactor = 1.d0
  material_property%secondary_continuum_name = ''
  material_property%secondary_continuum_length = 0.d0
  material_property%secondary_continuum_matrix_block_size = 0.d0
  material_property%secondary_continuum_fracture_spacing = 0.d0
  material_property%secondary_continuum_radius = 0.d0
  material_property%secondary_continuum_area = 0.d0
  material_property%secondary_continuum_epsilon = 1.d0
  material_property%secondary_continuum_aperture = 0.d0
  material_property%secondary_continuum_init_temp = 100.d0
  material_property%secondary_continuum_init_conc = 0.d0
  material_property%secondary_continuum_porosity = 0.5d0
  material_property%secondary_continuum_diff_coeff = 1.d-9
  material_property%secondary_continuum_mnrl_volfrac = 0.d0
  material_property%secondary_continuum_mnrl_area = 0.d0
  material_property%secondary_continuum_ncells = 0
  material_property%secondary_continuum_log_spacing = PETSC_FALSE
  material_property%secondary_continuum_outer_spacing = 1.d-3
  material_property%secondary_continuum_area_scaling = 1.d0
  nullify(material_property%next)
  MaterialPropertyCreate => material_property

end function MaterialPropertyCreate

! ************************************************************************** !

subroutine MaterialPropertyRead(material_property,input,option)
  ! 
  ! Reads in contents of a material_property card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Dataset_module
  use Units_module
  
  implicit none
  
  type(material_property_type) :: material_property
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, internal_units
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: buffer_save

  PetscInt :: length
  PetscBool :: therm_k_frz
  PetscBool :: therm_k_exp_frz
  PetscReal :: tempreal
  PetscInt, parameter :: TMP_SOIL_COMPRESSIBILITY = 1
  PetscInt, parameter :: TMP_BULK_COMPRESSIBILITY = 2
  PetscInt, parameter :: TMP_POROSITY_COMPRESSIBILITY = 3
  PetscInt :: soil_or_bulk_compressibility

  therm_k_frz = PETSC_FALSE
  therm_k_exp_frz = PETSC_FALSE
  soil_or_bulk_compressibility = UNINITIALIZED_INTEGER

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','MATERIAL_PROPERTY')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('NAME') 
        call InputReadWord(input,option,material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
      case('ID') 
        call InputReadInt(input,option,material_property%external_id)
        call InputErrorMsg(input,option,'id','MATERIAL_PROPERTY')
        if (material_property%external_id == UNINITIALIZED_INTEGER) then
          write(string,*) UNINITIALIZED_INTEGER
          option%io_buffer = 'Material ID "' // trim(adjustl(string)) // &
            '" is reserved for uninitialized materials.  Please choose a &
            &different value.'
        endif
      case('ACTIVE')
        material_property%active = PETSC_TRUE
      case('INACTIVE')
        material_property%active = PETSC_FALSE
      case('SATURATION_FUNCTION','CHARACTERISTIC_CURVES') 
        call InputReadWordDbaseCompatible(input,option, &
                           material_property%saturation_function_name, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'saturation function name', &
                           'MATERIAL_PROPERTY')
      case('ROCK_DENSITY') 
        call InputReadDouble(input,option,material_property%rock_density)
        call InputErrorMsg(input,option,'rock density','MATERIAL_PROPERTY')
        call InputReadAndConvertUnits(input,material_property%rock_density, &
                          'kg/m^3','MATERIAL_PROPERTY,rock density',option)
      case('SPECIFIC_HEAT','HEAT_CAPACITY') 
        call InputReadDouble(input,option,material_property%specific_heat)
        call InputErrorMsg(input,option,'specific heat','MATERIAL_PROPERTY')
        call InputReadAndConvertUnits(input,material_property%specific_heat, &
                          'J/kg-C','MATERIAL_PROPERTY,specific heat',option)
      case('LONGITUDINAL_DISPERSIVITY') 
        call InputReadDouble(input,option,material_property%dispersivity(1))
        call InputErrorMsg(input,option,'longitudinal_dispersivity', &
                           'MATERIAL_PROPERTY')
      case('TRANSVERSE_DISPERSIVITY_H') 
        call InputReadDouble(input,option,material_property%dispersivity(2))
        call InputErrorMsg(input,option,'transverse_dispersivity_h', &
                           'MATERIAL_PROPERTY')
      case('TRANSVERSE_DISPERSIVITY_V') 
        call InputReadDouble(input,option,material_property%dispersivity(3))
        call InputErrorMsg(input,option,'transverse_dispersivity_v', &
                           'MATERIAL_PROPERTY')
      case('THERMAL_CONDUCTIVITY_DRY') 
        call InputReadDouble(input,option, &
                             material_property%thermal_conductivity_dry)
        call InputErrorMsg(input,option,'dry thermal conductivity', &
                           'MATERIAL_PROPERTY')
        call InputReadAndConvertUnits(input, &
                   material_property%thermal_conductivity_dry, &
                   'W/m-C','MATERIAL_PROPERTY,dry thermal conductivity',option)
      case('THERMAL_CONDUCTIVITY_WET') 
        call InputReadDouble(input,option, &
                             material_property%thermal_conductivity_wet)
        call InputErrorMsg(input,option,'wet thermal conductivity', &
                           'MATERIAL_PROPERTY')
        call InputReadAndConvertUnits(input, &
                   material_property%thermal_conductivity_wet, &
                   'W/m-C','MATERIAL_PROPERTY,wet thermal conductivity',option)
      case('THERMAL_COND_EXPONENT') 
        call InputReadDouble(input,option, &
                             material_property%alpha)
        call InputErrorMsg(input,option,'thermal conductivity exponent', &
                           'MATERIAL_PROPERTY')
      case('THERMAL_CONDUCTIVITY_FROZEN') 
        therm_k_frz = PETSC_TRUE
        call InputReadDouble(input,option, &
                             material_property%thermal_conductivity_frozen)
        call InputErrorMsg(input,option,'frozen thermal conductivity', &
                           'MATERIAL_PROPERTY')
        call InputReadAndConvertUnits(input, &
                 material_property%thermal_conductivity_frozen, &
                 'W/m-C','MATERIAL_PROPERTY,frozen thermal conductivity',option)
      case('THERMAL_COND_EXPONENT_FROZEN') 
        therm_k_exp_frz = PETSC_TRUE
        call InputReadDouble(input,option, &
                             material_property%alpha_fr)
        call InputErrorMsg(input,option, &
                           'thermal conductivity frozen exponent', &
                           'MATERIAL_PROPERTY')
      case('SOIL_COMPRESSIBILITY_FUNCTION')
        call InputReadWord(input,option, &
                           material_property%soil_compressibility_function, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'soil compressibility function', &
                           'MATERIAL_PROPERTY')
      case('SOIL_COMPRESSIBILITY','BULK_COMPRESSIBILITY','POROSITY_COMPRESSIBILITY') 
        select case(keyword)
          case('SOIL_COMPRESSIBILITY') 
            soil_or_bulk_compressibility = TMP_SOIL_COMPRESSIBILITY
          case('BULK_COMPRESSIBILITY') 
            soil_or_bulk_compressibility = TMP_BULK_COMPRESSIBILITY
          case('POROSITY_COMPRESSIBILITY') 
            soil_or_bulk_compressibility = TMP_POROSITY_COMPRESSIBILITY
        end select
        call DatasetReadDoubleOrDataset(input,material_property% &
                                          soil_compressibility, &
                                   material_property%compressibility_dataset, &
                                        'soil compressibility', &
                                        'MATERIAL_PROPERTY',option)
      case('SOIL_REFERENCE_PRESSURE') 
        string = trim(input%buf)
        ! first read the word to determine if it is the keyword 
        ! INITIAL_CELL_PRESSURE.
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'soil reference pressure', &
                           'MATERIAL_PROPERTY')
        length = 16
        if (StringCompare(word,'INITIAL_PRESSURE',length)) then
          material_property%soil_reference_pressure_initial = PETSC_TRUE
        else
          ! if not the keyword above, copy back into buffer to be read as a
          ! double precision or dataset.
          input%buf = string
          call DatasetReadDoubleOrDataset(input, &
                  material_property%soil_reference_pressure, &
                  material_property%soil_reference_pressure_dataset, &
                  'soil reference pressure','MATERIAL_PROPERTY',option)
        endif
      case('THERMAL_EXPANSITIVITY') 
        call InputReadDouble(input,option, &
                             material_property%thermal_expansitivity)
        call InputErrorMsg(input,option,'thermal expansitivity', &
                           'MATERIAL_PROPERTY')
      case('POROSITY')
        call DatasetReadDoubleOrDataset(input,material_property%porosity, &
                                        material_property%porosity_dataset, &
                                        'porosity','MATERIAL_PROPERTY',option)
      case('TORTUOSITY')
        call DatasetReadDoubleOrDataset(input,material_property%tortuosity, &
                                        material_property%tortuosity_dataset, &
                                        'tortuosity','MATERIAL_PROPERTY',option)
      case('TORTUOSITY_FUNCTION_OF_POROSITY')
        material_property%tortuosity_function_of_porosity = PETSC_TRUE
        material_property%tortuosity = UNINITIALIZED_DOUBLE
        call InputReadDouble(input,option, &
                             material_property%tortuosity_func_porosity_pwr)
        call InputErrorMsg(input,option,'tortuosity power as a function of &
                           &porosity','MATERIAL_PROPERTY')
      case('PERMEABILITY')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'MATERIAL_PROPERTY,PERMEABILITY')
          
          if (InputCheckExit(input,option)) exit
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'MATERIAL_PROPERTY,PERMEABILITY')   
          select case(trim(word))
            case('ANISOTROPIC')
              material_property%isotropic_permeability = PETSC_FALSE
            case('VERTICAL_ANISOTROPY_RATIO')
              material_property%isotropic_permeability = PETSC_FALSE
              call InputReadDouble(input,option, &
                                   material_property%vertical_anisotropy_ratio)
              call InputErrorMsg(input,option,'vertical anisotropy ratio', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('ISOTROPIC')
              material_property%isotropic_permeability = PETSC_TRUE
            case('PERMEABILITY_SCALING_FACTOR')
              call InputReadDouble(input,option, &
                                  material_property%permeability_scaling_factor)
              call InputErrorMsg(input,option,'permeability scaling factor', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_X')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,1))
              call InputErrorMsg(input,option,'x permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_Y')
              call InputReadDouble(input,option, &
                                   material_property%permeability(2,2))
              call InputErrorMsg(input,option,'y permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_Z')
              call InputReadDouble(input,option, &
                                   material_property%permeability(3,3))
              call InputErrorMsg(input,option,'z permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_X_LOG10')
              call InputReadDouble(input,option, tempreal)
              call InputErrorMsg(input,option,'log10 x permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(1,1) = 10.d0**tempreal
            case('PERM_Y_LOG10')
              call InputReadDouble(input,option, tempreal)
              call InputErrorMsg(input,option,'log10 y permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(2,2) = 10.d0**tempreal
            case('PERM_Z_LOG10')
              call InputReadDouble(input,option, tempreal)
              call InputErrorMsg(input,option,'log10 z permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(3,3) = 10.d0**tempreal
            case('PERM_XZ')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,3))
              call InputErrorMsg(input,option,'xz permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_XY')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,2))
              call InputErrorMsg(input,option,'xy permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_YZ')
              call InputReadDouble(input,option, &
                                   material_property%permeability(2,3))
              call InputErrorMsg(input,option,'yz permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_XZ_LOG10')
              call InputReadDouble(input,option, tempreal)
              call InputErrorMsg(input,option,'log10 xz permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(1,3) = 10.d0**tempreal
            case('PERM_XY_LOG10')
              call InputReadDouble(input,option, tempreal)
              call InputErrorMsg(input,option,'log10 xy permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(1,2) = 10.d0**tempreal
            case('PERM_YZ_LOG10')
              call InputReadDouble(input,option, tempreal)
              call InputErrorMsg(input,option,'log10 yz permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(2,3) = 10.d0**tempreal
            case('PERM_ISO_LOG10')
              call InputReadDouble(input,option, tempreal)
              call InputErrorMsg(input,option,'log10 isotropic permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(1,1) = 10.d0**tempreal
              material_property%permeability(2,2) = 10.d0**tempreal
              material_property%permeability(3,3) = 10.d0**tempreal
            case('PERM_ISO')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,1))
              call InputErrorMsg(input,option,'isotropic permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(2,2) = &
                material_property%permeability(1,1)
              material_property%permeability(3,3) = &
                material_property%permeability(1,1)
            case('RANDOM_DATASET')
              option%io_buffer = 'RANDOM_DATASET is no longer supported.  ' // &
                'Please use the new DATASET object in the input file and ' // &
                'reference that dataset through "DATASET name" within ' // &
                'the PERMEABILITY card.'
              call printErrMsg(option)
            case('DATASET')
              material_property%permeability_dataset => DatasetBaseCreate()
              call InputReadNChars(input,option, &
                                   material_property% &
                                     permeability_dataset%name, &
                                   MAXWORDLENGTH,PETSC_TRUE)
              call InputErrorMsg(input,option,'DATASET,NAME', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')   
            case default
              call InputKeywordUnrecognized(word, &
                     'MATERIAL_PROPERTY,PERMEABILITY',option)
          end select
        enddo
        if (dabs(material_property%permeability(1,1) - &
                 material_property%permeability(2,2)) > 1.d-40 .or. &
            dabs(material_property%permeability(1,1) - &
                 material_property%permeability(3,3)) > 1.d-40) then
          material_property%isotropic_permeability = PETSC_FALSE
        endif
      case('PERM_FACTOR') 
      ! Permfactor is the multiplier to permeability to increase perm
      ! The perm increase could be due to pressure or other variable
      ! Added by Satish Karra, LANL, 1/8/12
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'MATERIAL_PROPERTY,PERM_FACTOR')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'MATERIAL_PROPERTY,PERM_FACTOR')   
          select case(trim(word))
          ! Assuming only ramp function for now
          ! The permfactor ramps from 1 to max_permfactor at max_pressure
          ! and remains same
            case('MIN_PRESSURE')       
              call InputReadDouble(input,option,material_property%min_pressure)
              call InputErrorMsg(input,option,'min pressure','PERM_FACTOR')  
            case('MAX_PRESSURE')       
              call InputReadDouble(input,option,material_property%max_pressure)
              call InputErrorMsg(input,option,'max pressure','PERM_FACTOR')
            case('MAX_PERMFACTOR')       
              call InputReadDouble(input,option,material_property%max_permfactor)
              call InputErrorMsg(input,option,'max permfactor','PERM_FACTOR')
            case default
              call InputKeywordUnrecognized(word, &
                     'MATERIAL_PROPERTY,PERM_FACTOR',option)
          end select
        enddo
      case('PERMEABILITY_POWER')
        call InputReadDouble(input,option, &
                             material_property%permeability_pwr)
        call InputErrorMsg(input,option,'permeability power', &
                           'MATERIAL_PROPERTY')
      case('PERMEABILITY_CRITICAL_POROSITY')
        call InputReadDouble(input,option, &
                             material_property%permeability_crit_por)
        call InputErrorMsg(input,option,'permeability critical porosity', &
                           'MATERIAL_PROPERTY')
      case('PERMEABILITY_MIN_SCALE_FACTOR')
        call InputReadDouble(input,option, &
                             material_property%permeability_min_scale_fac)
        call InputErrorMsg(input,option,'permeability min scale factor', &
                           'MATERIAL_PROPERTY')
      case('TORTUOSITY_POWER')
        call InputReadDouble(input,option, &
                             material_property%tortuosity_pwr)
        call InputErrorMsg(input,option,'tortuosity power','MATERIAL_PROPERTY')
      case('MINERAL_SURFACE_AREA_POWER')
        option%io_buffer = 'Adjustment of mineral surface area based on ' // &
          'mineral volume fraction or porosity must be performed on a ' // &
          'per mineral basis under the MINERAL_KINETICS card.  See ' // &
          'reaction_aux.F90.'
          call printErrMsg(option)
      case('SECONDARY_CONTINUUM')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'MATERIAL_PROPERTY,SECONDARY_CONTINUUM')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'MATERIAL_PROPERTY,SECONDARY_CONTINUUM')   
          select case(trim(word))
            case('TYPE')
              call InputReadNChars(input,option, &
                                   material_property%secondary_continuum_name,&
                                   MAXWORDLENGTH,PETSC_TRUE)
              call InputErrorMsg(input,option,'type', &
                                'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('MATRIX_BLOCK_SIZE')
              call InputReadDouble(input,option, &
                        material_property%secondary_continuum_matrix_block_size)
              call InputErrorMsg(input,option,'matrix_block_size', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('FRACTURE_SPACING')
              call InputReadDouble(input,option, &
                        material_property%secondary_continuum_fracture_spacing)
              call InputErrorMsg(input,option,'fracture_spacing', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('RADIUS')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_radius)
              call InputErrorMsg(input,option,'radius', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('LENGTH')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_length)
              call InputErrorMsg(input,option,'length', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('AREA')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_area)
              call InputErrorMsg(input,option,'area', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('NUM_CELLS')
              call InputReadInt(input,option, &
                                   material_property%secondary_continuum_ncells)
              call InputErrorMsg(input,option,'number of cells', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('EPSILON')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_epsilon)
              call InputErrorMsg(input,option,'epsilon', &
                           'MATERIAL_PROPERTY')
            case('APERTURE')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_aperture)
              call InputErrorMsg(input,option,'aperture', &
                           'MATERIAL_PROPERTY')
            case('TEMPERATURE')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_init_temp)
              call InputErrorMsg(input,option,'secondary continuum init temp', &
                           'MATERIAL_PROPERTY')
              option%set_secondary_init_temp = PETSC_TRUE
            case('CONCENTRATION')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_init_conc)
              call InputErrorMsg(input,option,'secondary continuum init conc', &
                           'MATERIAL_PROPERTY')
              option%set_secondary_init_conc = PETSC_TRUE
            case('POROSITY')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_porosity)
              call InputErrorMsg(input,option,'secondary continuum porosity', &
                           'MATERIAL_PROPERTY')
            case('DIFFUSION_COEFFICIENT')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_diff_coeff)
              call InputErrorMsg(input,option, &
                                 'secondary continuum diff coeff', &
                                 'MATERIAL_PROPERTY')
            case('MINERAL_VOLFRAC')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_mnrl_volfrac)
              call InputErrorMsg(input,option,'secondary cont. mnrl volfrac.', &
                           'MATERIAL_PROPERTY')  
            case('MINERAL_AREA')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_mnrl_area)
              call InputErrorMsg(input,option,'secondary cont. mnrl area', &
                           'MATERIAL_PROPERTY')
            case('LOG_GRID_SPACING')
              material_property%secondary_continuum_log_spacing = PETSC_TRUE
            case('OUTER_SPACING')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_outer_spacing)
              call InputErrorMsg(input,option,'secondary cont. outer spacing', &
                           'MATERIAL_PROPERTY')
            case('AREA_SCALING_FACTOR')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_area_scaling)
              call InputErrorMsg(input,option,'secondary area scaling factor', &
                           'MATERIAL_PROPERTY')
            case default
              call InputKeywordUnrecognized(word, &
                     'MATERIAL_PROPERTY,SECONDARY_CONTINUUM',option)
          end select
        enddo

      case default
        call InputKeywordUnrecognized(keyword,'MATERIAL_PROPERTY',option)
    end select 
  enddo
  
  if (material_property%tortuosity_function_of_porosity) then
    if (associated(material_property%tortuosity_dataset)) then
      option%io_buffer = 'A TORTUOSITY dataset may not be assigned in &
        &combination with TORTUOSITY_FUNCTION_OF_POROSITY.'
      call printErrMsg(option)
    endif
    if (.not.associated(material_property%porosity_dataset)) then
      material_property%tortuosity = material_property%porosity** &
        material_property%tortuosity_func_porosity_pwr
    endif
  endif

  if (associated(material_property%permeability_dataset) .and. &
      .not.material_property%isotropic_permeability .and. &
      Uninitialized(material_property%vertical_anisotropy_ratio)) then
    material_property%permeability_dataset_y => DatasetBaseCreate()
    material_property%permeability_dataset_z => DatasetBaseCreate()
    material_property%permeability_dataset_y%name = &
      trim(material_property%permeability_dataset%name) // 'Y'
    material_property%permeability_dataset_z%name = &
      trim(material_property%permeability_dataset%name) // 'Z'
    material_property%permeability_dataset%name = &
      trim(material_property%permeability_dataset%name) // 'X'
  endif

  if (option%iflowmode == TH_MODE) then
     if (option%use_th_freezing .eqv. PETSC_TRUE) then
        if (.not. therm_k_frz) then
           option%io_buffer = 'THERMAL_CONDUCTIVITY_FROZEN must be set &
             &in inputdeck for MODE TH(C) ICE'
           call printErrMsg(option)
        endif
        if (.not. therm_k_exp_frz) then
           option%io_buffer = 'THERMAL_COND_EXPONENT_FROZEN must be set &
             &in inputdeck for MODE TH(C) ICE'
           call printErrMsg(option)
        endif
     endif
  endif

  if (len_trim(material_property%soil_compressibility_function) > 0) then
    word = material_property%soil_compressibility_function
    select case(word)
      case('BRAGFLO','BULK_EXPONENTIAL')
        if (soil_or_bulk_compressibility /= TMP_BULK_COMPRESSIBILITY) then
          option%io_buffer = 'A BULK_COMPRESSIBILITY should be entered &
            &instead of a SOIL_COMPRESSIBILITY in MATERIAL_PROPERTY "' // &
            trim(material_property%name) // '" since a BRAGFLO or WIPP &
            &SOIL_COMPRESSIBILITY function is defined.'
          call printErrMsg(option)
        endif
        word = 'BULK_COMPRESSIBILITY'
      case('POROSITY_EXPONENTIAL')
        if (soil_or_bulk_compressibility /= TMP_POROSITY_COMPRESSIBILITY) then
          option%io_buffer = 'A POROSITY_COMPRESSIBILITY should be entered &
            &in MATERIAL_PROPERTY "' // &
            trim(material_property%name) // '" since a POROSITY_EXPONENTIAL &
            &not POROSITY_COMPRESSIBILITY function is defined.'
          call printErrMsg(option)
        endif
        word = 'POROSITY_COMPRESSIBILITY'
      case('LEIJNSE','DEFAULT')
        if (soil_or_bulk_compressibility /= TMP_SOIL_COMPRESSIBILITY) then
          option%io_buffer = 'A SOIL_COMPRESSIBILITY should be entered &
            &instead of a BULK_COMPRESSIBILITY in MATERIAL_PROPERTY "' // &
            trim(material_property%name) // '" since a LEIJNSE or DEFAULT &
            &SOIL_COMPRESSIBILITY function is defined.'
          call printErrMsg(option)
        endif
        word = 'SOIL_COMPRESSIBILITY'
      case default
    end select
    option%flow%transient_porosity = PETSC_TRUE
    if (Uninitialized(material_property%soil_compressibility) .and. &
        .not.associated(material_property%compressibility_dataset)) then
      option%io_buffer = 'SOIL_COMPRESSIBILITY_FUNCTION is specified in &
        &inputdeck for MATERIAL_PROPERTY "' // &
        trim(material_property%name) // &
        '", but a ' // trim(word) // ' is not defined.'
      call printErrMsg(option)
    endif
    if (Uninitialized(material_property%soil_reference_pressure) .and. &
        .not.associated(material_property% &
                          soil_reference_pressure_dataset) .and. &
        .not.material_property%soil_reference_pressure_initial) then
      option%io_buffer = 'SOIL_COMPRESSIBILITY_FUNCTION is specified in &
        &inputdeck for MATERIAL_PROPERTY "' // &
        trim(material_property%name) // &
        '", but a SOIL_REFERENCE_PRESSURE is not defined.'
      call printErrMsg(option)
    endif
    if ((Initialized(material_property%soil_reference_pressure) .or. &
         associated(material_property%soil_reference_pressure_dataset)) .and. &
        material_property%soil_reference_pressure_initial) then
      option%io_buffer = 'SOIL_REFERENCE_PRESSURE may not be defined by the &
        &initial pressure and a specified pressure in material "' // &
        trim(material_property%name) // '".'
      call printErrMsg(option)
    endif
  endif

  ! material id must be > 0
  if (material_property%external_id <= 0) then
    write(word,*) material_property%external_id
    option%io_buffer = 'Material ID in MATERIAL_PROPERTY "' // &
      trim(material_property%name) // '" must be > 0 (' // &
      trim(adjustl(word)) // '). If you would like to inactivate a &
      &material, please do so by adding INACTIVE to the STRATA to which &
      &the MATERIAL_PROPERTY is coupled.'
    call printErrMsg(option)
  endif

end subroutine MaterialPropertyRead

! ************************************************************************** !

subroutine MaterialPropertyAddToList(material_property,list)
  ! 
  ! Adds a material property to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  implicit none
  
  type(material_property_type), pointer :: material_property
  type(material_property_type), pointer :: list

  type(material_property_type), pointer :: cur_material_property
  
  if (associated(list)) then
    cur_material_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_material_property%next)) exit
      cur_material_property => cur_material_property%next
    enddo
    cur_material_property%next => material_property
    material_property%internal_id = abs(cur_material_property%internal_id) + 1
  else
    list => material_property
    material_property%internal_id = 1
  endif
  if (.not.material_property%active) then
    material_property%internal_id = -1*material_property%internal_id
  endif
  
end subroutine MaterialPropertyAddToList

! ************************************************************************** !

subroutine MaterialPropConvertListToArray(list,array,option)
  ! 
  ! Creates an array of pointers to the
  ! material_properties in the list
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/07
  ! 

  use Option_module
  use String_module

  implicit none
  
  type(material_property_type), pointer :: list
  type(material_property_ptr_type), pointer :: array(:)
  type(option_type) :: option
    
  type(material_property_type), pointer :: cur_material_property
  type(material_property_type), pointer :: prev_material_property
  type(material_property_type), pointer :: next_material_property
  PetscInt :: i, j, length1,length2, max_internal_id, max_external_id
  PetscInt, allocatable :: id_count(:)
  PetscBool :: error_flag
  character(len=MAXSTRINGLENGTH) :: string

#if 0
! don't necessary need right now, but maybe in future
  ! reorder into ascending order
  swapped = PETSC_FALSE
  do
    if (.not.swapped) exit
    cur_material_property => list
    do 
      if (.not.associated(cur_material_property)) exit
      next_material_property => cur_material_property%next
      if (associated(next_material_property)) then
        if (cur_material_property%id > next_material_property%id) then
          ! swap
          if (associated(prev_material_property)) then
            prev_material_property%next => next_material_property
          else
            list => next_material_property
          endif
          cur_material_property%next => next_material_property%next
          next_material_property%next => cur_material_property
          swapped = PETSC_TRUE
        endif
      endif
      prev_material_property => cur_material_property
      cur_material_property => next_material_property
    enddo
  enddo
#endif

  ! check to ensure that max internal id is equal to the number of 
  ! material properties and that internal ids are contiguous
  max_internal_id = 0
  max_external_id = 0
  cur_material_property => list
  do 
    if (.not.associated(cur_material_property)) exit
    max_internal_id = max_internal_id + 1
    max_external_id = max(max_external_id,cur_material_property%external_id)
    if (max_internal_id /= abs(cur_material_property%internal_id)) then
      write(string,*) cur_material_property%external_id
      option%io_buffer = 'Non-contiguous internal material id for ' // &
        'material named "' // trim(cur_material_property%name) // &
        '" with external id "' // trim(adjustl(string)) // '" '
      write(string,*) cur_material_property%internal_id
      option%io_buffer = trim(option%io_buffer) // &
        'and internal id "' // trim(adjustl(string)) // '".'
      call printErrMsg(option)
    endif
    cur_material_property => cur_material_property%next
  enddo
  
  if (associated(array)) deallocate(array)
  allocate(array(max_internal_id))
  do i = 1, max_internal_id
    nullify(array(i)%ptr)
  enddo
  
  ! use id_count to ensure that an id is not duplicated
  allocate(id_count(max_external_id))
  id_count = 0
  
  cur_material_property => list
  do 
    if (.not.associated(cur_material_property)) exit
    id_count(cur_material_property%external_id) = &
      id_count(cur_material_property%external_id) + 1
    array(abs(cur_material_property%internal_id))%ptr => cur_material_property
    cur_material_property => cur_material_property%next
  enddo
  
  ! check to ensure that an id is not duplicated
  error_flag = PETSC_FALSE
  do i = 1, max_external_id
    if (id_count(i) > 1) then
      write(string,*) i
      option%io_buffer = 'Material ID ' // trim(adjustl(string)) // &
        ' is duplicated in input file.'
      call printMsg(option)
      error_flag = PETSC_TRUE
    endif
  enddo

  deallocate(id_count)

  if (error_flag) then
    option%io_buffer = 'Duplicate Material IDs.'
    call printErrMsg(option)
  endif
  
  ! ensure unique material names
  error_flag = PETSC_FALSE
  do i = 1, size(array)
    if (associated(array(i)%ptr)) then
      length1 = len_trim(array(i)%ptr%name)
      do j = 1, i-1
        if (associated(array(j)%ptr)) then
          length2 = len_trim(array(j)%ptr%name)
          if (length1 /= length2) cycle
          if (StringCompare(array(i)%ptr%name,array(j)%ptr%name,length1)) then
            option%io_buffer = 'Material name "' // &
              trim(adjustl(array(i)%ptr%name)) // &
              '" is duplicated in input file.'
            call printMsg(option)
            error_flag = PETSC_TRUE
          endif
        endif
      enddo
    endif
  enddo

  if (error_flag) then
    option%io_buffer = 'Duplicate Material names.'
    call printErrMsg(option)
  endif
  
end subroutine MaterialPropConvertListToArray

! ************************************************************************** !

function MaterialGetMaxExternalID(material_property_array)
  ! 
  ! Maps internal material ids to external for I/O, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/14
  ! 
  implicit none
  
  type(material_property_ptr_type) :: material_property_array(:)
  
  PetscInt :: MaterialGetMaxExternalID
  
  PetscInt :: i

  MaterialGetMaxExternalID = UNINITIALIZED_INTEGER
  do i = 1, size(material_property_array)
    MaterialGetMaxExternalID = max(MaterialGetMaxExternalID, &
                                  (material_property_array(i)%ptr%external_id))
  enddo

end function MaterialGetMaxExternalID

! ************************************************************************** !

subroutine MaterialCreateIntToExtMapping(material_property_array,mapping)
  ! 
  ! Maps internal material ids to external for I/O, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/14
  ! 
  implicit none
  
  type(material_property_ptr_type) :: material_property_array(:)
  PetscInt, pointer :: mapping(:)

  PetscInt :: i
  
  allocate(mapping(0:size(material_property_array)))
  mapping = UNINITIALIZED_INTEGER
  mapping(0) = 0
  
  do i = 1, size(material_property_array)
    mapping(abs(material_property_array(i)%ptr%internal_id)) = &
      material_property_array(i)%ptr%external_id
  enddo

end subroutine MaterialCreateIntToExtMapping

! ************************************************************************** !

subroutine MaterialCreateExtToIntMapping(material_property_array,mapping)
  ! 
  ! Maps external material ids to internal for setup. This array should be 
  ! temporary and never stored for the duration of the simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/14
  ! 
  implicit none
  
  type(material_property_ptr_type) :: material_property_array(:)
  PetscInt, pointer :: mapping(:)
  
  PetscInt :: i
  
  allocate(mapping(0:MaterialGetMaxExternalID(material_property_array)))
  mapping = -888
  mapping(0) = 0
  
  do i = 1, size(material_property_array)
    mapping(material_property_array(i)%ptr%external_id) = &
      material_property_array(i)%ptr%internal_id
  enddo  

end subroutine MaterialCreateExtToIntMapping

! ************************************************************************** !

subroutine MaterialApplyMapping(mapping,array)
  ! 
  ! Maps internal material ids to external for I/O, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/14
  ! 
  implicit none
  
  PetscInt :: mapping(0:)
  PetscInt :: array(:)

  PetscInt :: i
  PetscInt :: mapping_size
  PetscInt :: mapped_id

  mapping_size = size(mapping)-1 ! subtract 1 for 0 index
  do i = 1, size(array)
    if (array(i) <= mapping_size) then
      mapped_id = mapping(array(i))
    else
      mapped_id = -888 ! indicates corresponding mapped value does not exist.
    endif
    array(i) = mapped_id
  enddo

end subroutine MaterialApplyMapping

! ************************************************************************** !

subroutine MaterialSetup(material_parameter, material_property_array, &
                         characteristic_curves_array, option)
  ! 
  ! Creates arrays for material parameter object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/14
  !
  use Option_module
  use Characteristic_Curves_module
  
  implicit none
  
  type(material_parameter_type) :: material_parameter
  type(material_property_ptr_type) :: material_property_array(:)
  type(characteristic_curves_ptr_type) :: characteristic_curves_array(:)
  type(option_type), pointer :: option
  
  PetscInt :: num_characteristic_curves
  PetscInt :: num_mat_prop
  PetscInt :: i
  
  num_mat_prop = size(material_property_array)
  num_characteristic_curves = size(characteristic_curves_array)
  
  allocate(material_parameter%soil_residual_saturation(option%nphase, &
                                                   num_characteristic_curves))
  material_parameter%soil_residual_saturation = UNINITIALIZED_DOUBLE
  do i = 1, num_characteristic_curves
    if (associated(characteristic_curves_array(i)%ptr)) then
      material_parameter%soil_residual_saturation(:,i) = &
        CharCurvesGetGetResidualSats(characteristic_curves_array(i)%ptr,option)
    endif
  enddo

  if (option%iflowmode == TH_MODE) then
    allocate(material_parameter%soil_heat_capacity(num_mat_prop))
    allocate(material_parameter%soil_thermal_conductivity(2,num_mat_prop))
    material_parameter%soil_heat_capacity = UNINITIALIZED_DOUBLE
    material_parameter%soil_thermal_conductivity = UNINITIALIZED_DOUBLE
    do i = 1, num_mat_prop
      if (associated(material_property_array(i)%ptr)) then
        ! kg rock/m^3 rock * J/kg rock-K * 1.e-6 MJ/J
        material_parameter%soil_heat_capacity(i) = &
          material_property_array(i)%ptr%specific_heat * option%scale ! J -> MJ
        material_parameter%soil_thermal_conductivity(1,i) = &
          material_property_array(i)%ptr%thermal_conductivity_dry
        material_parameter%soil_thermal_conductivity(2,i) = &
          material_property_array(i)%ptr%thermal_conductivity_wet
      endif
    enddo
  endif
  
end subroutine MaterialSetup
  
! ************************************************************************** !

function MaterialPropGetPtrFromList(material_property_name, &
                                    material_property_list)
  ! 
  ! Returns a pointer to the material property
  ! matching material_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  use String_module
  
  implicit none
  
  type(material_property_type), pointer :: MaterialPropGetPtrFromList
  character(len=MAXWORDLENGTH) :: material_property_name
  type(material_property_type), pointer :: material_property_list
  PetscInt :: length
  type(material_property_type), pointer :: material_property
    
  nullify(MaterialPropGetPtrFromList)
  material_property => material_property_list
  
  do 
    if (.not.associated(material_property)) exit
    length = len_trim(material_property_name)
    if (length == len_trim(material_property%name) .and. &
        StringCompare(material_property%name,material_property_name,length)) then
      MaterialPropGetPtrFromList => material_property
      return
    endif
    material_property => material_property%next
  enddo
  
end function MaterialPropGetPtrFromList

! ************************************************************************** !

function MaterialPropGetPtrFromArray(material_property_name, &
                                     material_property_array)
  ! 
  ! Returns a pointer to the material property
  ! matching material_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  use String_module

  implicit none
  
  type(material_property_type), pointer :: MaterialPropGetPtrFromArray
  character(len=MAXWORDLENGTH) :: material_property_name
  type(material_property_ptr_type), pointer :: material_property_array(:)
  PetscInt :: length
  PetscInt :: imaterial_property
    
  nullify(MaterialPropGetPtrFromArray)
  
  do imaterial_property = 1, size(material_property_array)
    length = len_trim(material_property_name)
    if (.not.associated(material_property_array(imaterial_property)%ptr)) cycle
    if (length == &
        len_trim(material_property_array(imaterial_property)%ptr%name) .and. &
        StringCompare(material_property_array(imaterial_property)%ptr%name, &
                        material_property_name,length)) then
      MaterialPropGetPtrFromArray => &
        material_property_array(imaterial_property)%ptr
      return
    endif
  enddo
  
end function MaterialPropGetPtrFromArray

! ************************************************************************** !

function MaterialAnisotropyExists(material_property_list)
  ! 
  ! Determines whether any of the material
  ! properties are anisotropic
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/13
  ! 

  implicit none
  
  type(material_property_type), pointer :: material_property_list

  PetscBool :: MaterialAnisotropyExists
  
  type(material_property_type), pointer :: cur_material_property
    
  MaterialAnisotropyExists = PETSC_FALSE
  
  cur_material_property => material_property_list
  do 
    if (.not.associated(cur_material_property)) exit
    if (.not. cur_material_property%isotropic_permeability) then
      MaterialAnisotropyExists = PETSC_TRUE
      return
    endif
    cur_material_property => cur_material_property%next
  enddo
  
end function MaterialAnisotropyExists


! ************************************************************************** !

subroutine MaterialInitAuxIndices(material_property_ptrs,option)
  !
  ! Initializes the pointer used to index material property arrays
  !
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  !
  use Material_Aux_class
  use String_module
  use Option_module
  
  implicit none
  
  type(material_property_ptr_type), pointer :: material_property_ptrs(:)
  type(option_type) :: option

  PetscInt :: i
  PetscInt :: icount = 0
  PetscInt :: num_soil_compress_func = 0
  PetscInt :: num_soil_compress = 0
  PetscInt :: num_soil_ref_press = 0
  PetscInt :: num_material_properties
  
  procedure(MaterialCompressSoilDummy), pointer :: &
    MaterialCompressSoilPtrTmp 
  
  soil_compressibility_index = 0
  soil_reference_pressure_index = 0
  soil_heat_capacity_index = 0
  soil_thermal_conductivity_index = 0
  max_material_index = 0
  
  num_material_properties = size(material_property_ptrs)
  ! must be nullified here to avoid an error message on subsequent calls
  ! on stochastic simulations
  MaterialCompressSoilPtr => null()
  
  do i = 1, num_material_properties
    MaterialCompressSoilPtrTmp => null()
    if (len_trim(material_property_ptrs(i)%ptr% &
                   soil_compressibility_function) > 1) then
      call StringToUpper(material_property_ptrs(i)%ptr% &
                           soil_compressibility_function)
      select case(material_property_ptrs(i)%ptr%soil_compressibility_function)
        case('BRAGFLO','BULK_EXPONENTIAL')
          MaterialCompressSoilPtrTmp => MaterialCompressSoilBRAGFLO
        case('POROSITY_EXPONENTIAL')
          MaterialCompressSoilPtrTmp => MaterialCompressSoilPoroExp
        case('QUADRATIC')
          MaterialCompressSoilPtrTmp => MaterialCompressSoilQuadratic
        case('LEIJNSE','DEFAULT')
          MaterialCompressSoilPtrTmp => MaterialCompressSoilLeijnse
        case('LINEAR')
          MaterialCompressSoilPtrTmp => MaterialCompressSoilLinear
        case default
          option%io_buffer = 'Soil compressibility function "' // &
            trim(material_property_ptrs(i)%ptr% &
                   soil_compressibility_function) // &
            '" not recognized.'
          call printErrMsg(option)
      end select
      num_soil_compress_func = num_soil_compress_func + 1
    endif  
    if (.not.associated(MaterialCompressSoilPtr)) then
      MaterialCompressSoilPtr => MaterialCompressSoilPtrTmp
    else if (.not.associated(MaterialCompressSoilPtr, &
                             MaterialCompressSoilPtrTmp)) then
      option%io_buffer = 'All MATERIAL_PROPERTIES must specify the ' // &
        'same soil compressibility function.'
      call printErrMsg(option)
    endif
    if (Initialized(material_property_ptrs(i)%ptr%soil_compressibility) .or. &
        associated(material_property_ptrs(i)%ptr%compressibility_dataset)) then
      if (soil_compressibility_index == 0) then
        icount = icount + 1
        soil_compressibility_index = icount
      endif
      num_soil_compress = num_soil_compress + 1
    endif
    if (Initialized(material_property_ptrs(i)%ptr%&
                      soil_reference_pressure) .or. &
        associated(material_property_ptrs(i)%ptr%&
                      soil_reference_pressure_dataset) .or. &
        material_property_ptrs(i)%ptr%soil_reference_pressure_initial) then
      if (soil_reference_pressure_index == 0) then
        icount = icount + 1
        soil_reference_pressure_index = icount
      endif
      num_soil_ref_press = num_soil_ref_press + 1
    endif
    if (material_property_ptrs(i)%ptr%specific_heat > 0.d0 .and. &
        soil_heat_capacity_index == 0) then
      icount = icount + 1
      soil_heat_capacity_index = icount
    endif
    if (material_property_ptrs(i)%ptr%thermal_conductivity_wet > 0.d0 .and. &
        soil_thermal_conductivity_index == 0) then
      icount = icount + 1
      soil_thermal_conductivity_index = icount
    endif
  enddo
  max_material_index = icount
  
  if (.not.associated(MaterialCompressSoilPtr)) then
    MaterialCompressSoilPtr => MaterialCompressSoilLeijnse
  endif

  ! check of uninitialized values
  if (num_soil_compress_func > 0 .and. &
      num_soil_compress_func /= num_material_properties) then
    option%io_buffer = 'SOIL_COMPRESSIBILITY_FUNCTION must be defined for all &
      &materials.'
    call printErrMsg(option)
  endif
  if (soil_compressibility_index > 0 .and. &
      num_soil_compress /= num_material_properties) then
    option%io_buffer = 'SOIL_COMPRESSIBILITY must be defined for all &
      &materials.'
    call printErrMsg(option)
  endif
  if (soil_reference_pressure_index > 0 .and. &
      num_soil_ref_press /= num_material_properties) then
    option%io_buffer = 'SOIL_REFERENCE_PRESSURE must be defined for all &
      &materials.'
    call printErrMsg(option)
  endif
  if (soil_compressibility_index > 0 .and. &
      soil_reference_pressure_index == 0) then
    option%io_buffer = 'SOIL_REFERENCE_PRESSURE must be defined to model &
      &soil compressibility.'
    call printErrMsg(option)
  endif
  
end subroutine MaterialInitAuxIndices

! ************************************************************************** !

subroutine MaterialAssignPropertyToAux(material_auxvar,material_property, &
                                       option)
  !
  ! Initializes the pointer used to index material property arrays
  !
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  !
  use Material_Aux_class
  use Option_module
  
  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  type(material_property_type) :: material_property
  type(option_type) :: option

  if (Initialized(material_property%rock_density)) then
    material_auxvar%soil_particle_density = &
      material_property%rock_density
  endif

  if (soil_compressibility_index > 0) then
    material_auxvar%soil_properties(soil_compressibility_index) = &
      material_property%soil_compressibility
  endif
  if (soil_reference_pressure_index > 0) then
    ! soil reference pressure may be assigned as the initial cell pressure, and
    ! in that case, it will be assigned elsewhere
    if (Initialized(material_property%soil_reference_pressure)) then
      material_auxvar%soil_properties(soil_reference_pressure_index) = &
        material_property%soil_reference_pressure
    endif
  endif
  if (soil_heat_capacity_index > 0) then
    material_auxvar%soil_properties(soil_heat_capacity_index) = &
      material_property%specific_heat
  endif
  if (soil_thermal_conductivity_index > 0) then
    material_auxvar%soil_properties(soil_thermal_conductivity_index) = &
      material_property%thermal_conductivity_wet
  endif
  
end subroutine MaterialAssignPropertyToAux

! ************************************************************************** !

subroutine MaterialSetAuxVarScalar(Material,value,ivar)
  ! 
  ! Sets values of a material auxvar data using a scalar value.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Variables_module
  
  implicit none

  type(material_type) :: Material ! from realization%patch%aux%Material
  PetscReal :: value
  PetscInt :: ivar

  PetscInt :: i
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
!  material_auxvars => Material%auxvars
!geh: can't use this pointer as gfortran does not like it.  Must use
!     Material%auxvars%....

  select case(ivar)
    case(VOLUME)
      do i=1, Material%num_aux
        Material%auxvars(i)%volume = value
      enddo
    case(POROSITY)
      do i=1, Material%num_aux
        Material%auxvars(i)%porosity = value
      enddo
    case(TORTUOSITY)
      do i=1, Material%num_aux
        Material%auxvars(i)%tortuosity = value
      enddo
    case(PERMEABILITY_X)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_xx_index) = value
      enddo
    case(PERMEABILITY_Y)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_yy_index) = value
      enddo
    case(PERMEABILITY_Z)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_zz_index) = value
      enddo
    case(PERMEABILITY_XY)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_xy_index) = value
      enddo
    case(PERMEABILITY_YZ)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_yz_index) = value
      enddo
    case(PERMEABILITY_XZ)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_xz_index) = value
      enddo
  end select
  
end subroutine MaterialSetAuxVarScalar

! ************************************************************************** !

subroutine MaterialSetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)
  ! 
  ! Sets values of material auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Variables_module
  
  implicit none

  type(material_type) :: Material ! from realization%patch%aux%Material
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscErrorCode :: ierr
  
!  material_auxvars => Material%auxvars
!geh: can't use this pointer as gfortran does not like it.  Must use
!     Material%auxvars%....
  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(SOIL_COMPRESSIBILITY)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)% &
          soil_properties(soil_compressibility_index) = vec_loc_p(ghosted_id)
      enddo
    case(SOIL_REFERENCE_PRESSURE)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)% &
          soil_properties(soil_reference_pressure_index) = vec_loc_p(ghosted_id)
      enddo
    case(VOLUME)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%volume = vec_loc_p(ghosted_id)
      enddo
    case(POROSITY)
      select case(isubvar)
        case(POROSITY_CURRENT)
          do ghosted_id=1, Material%num_aux
            Material%auxvars(ghosted_id)%porosity = vec_loc_p(ghosted_id)
          enddo
        case(POROSITY_MINERAL)
          do ghosted_id=1, Material%num_aux
            Material%auxvars(ghosted_id)%porosity_base = vec_loc_p(ghosted_id)
          enddo
      end select
    case(TORTUOSITY)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%tortuosity = vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_X)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_xx_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_Y)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_yy_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_Z)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_zz_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_XY)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_xy_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_YZ)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_yz_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_XZ)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_xz_index) = &
          vec_loc_p(ghosted_id)
      enddo
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine MaterialSetAuxVarVecLoc

! ************************************************************************** !

subroutine MaterialGetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)
  ! 
  ! Gets values of material auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Variables_module
  
  implicit none

  type(material_type) :: Material ! from realization%patch%aux%Material
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscErrorCode :: ierr
  
!  material_auxvars => Material%auxvars
!geh: can't use this pointer as gfortran does not like it.  Must use
!     Material%auxvars%....
  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(SOIL_COMPRESSIBILITY)
      if (soil_compressibility_index > 0) then
        do ghosted_id=1, Material%num_aux
          vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)% &
                                     soil_properties(soil_compressibility_index)
        enddo
      else
        vec_loc_p(:) = UNINITIALIZED_DOUBLE
      endif
    case(SOIL_REFERENCE_PRESSURE)
      if (soil_reference_pressure_index > 0) then
        do ghosted_id=1, Material%num_aux
          vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)% &
                                  soil_properties(soil_reference_pressure_index)
        enddo
      else
        vec_loc_p(:) = UNINITIALIZED_DOUBLE
      endif
    case(VOLUME)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)%volume
      enddo
    case(POROSITY)
      select case(isubvar)
        case(POROSITY_CURRENT)
          do ghosted_id=1, Material%num_aux
            vec_loc_p(ghosted_id) = &
              Material%auxvars(ghosted_id)%porosity
          enddo
        case(POROSITY_MINERAL)
          do ghosted_id=1, Material%num_aux
            vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)%porosity_base
          enddo
      end select
    case(TORTUOSITY)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)%tortuosity
      enddo
    case(PERMEABILITY_X)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_xx_index)
      enddo
    case(PERMEABILITY_Y)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_yy_index)
      enddo
    case(PERMEABILITY_Z)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_zz_index)
      enddo
    case(PERMEABILITY_XY)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_xy_index)
      enddo
    case(PERMEABILITY_YZ)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_yz_index)
      enddo
    case(PERMEABILITY_XZ)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_xz_index)
      enddo
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine MaterialGetAuxVarVecLoc

! ************************************************************************** !

subroutine MaterialWeightAuxVars(Material,weight,field,comm1)
  ! 
  ! Updates the porosities in auxiliary variables associated with 
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/17/14
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Field_module
  use Communicator_Base_module
  use Variables_module, only : POROSITY
  
  implicit none

  type(material_type) :: Material
  type(field_type) :: field
  PetscReal :: weight
  class(communicator_type) :: comm1
  
  PetscErrorCode :: ierr
  
!  material_auxvars => Material%auxvars
!geh: can't use this pointer as gfortran does not like it.  Must use
!     Material%auxvars%....
  call VecCopy(field%porosity_t,field%work,ierr)
  call VecAXPBY(field%work,weight,1.d0-weight, &
                field%porosity_tpdt,ierr);CHKERRQ(ierr)
  call comm1%GlobalToLocal(field%work,field%work_loc)
  call MaterialSetAuxVarVecLoc(Material,field%work_loc,POROSITY, &
                               POROSITY_CURRENT)
  
end subroutine MaterialWeightAuxVars
 
! ************************************************************************** !

subroutine MaterialStoreAuxVars(Material,time)
  ! 
  ! Moves material properties from TIME_TpDT -> TIME_T in storage arrays
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/14
  ! 

  use Option_module

  implicit none

  type(material_type) :: Material
  PetscReal :: time
  
  PetscInt :: ghosted_id
  
  Material%time_t = time
  
  do ghosted_id=1, Material%num_aux
!    Material%auxvars(ghosted_id)%porosity_store(TIME_T) = &
!      Material%auxvars(ghosted_id)%porosity_store(TIME_TpDT)
  enddo

end subroutine MaterialStoreAuxVars

! ************************************************************************** !

subroutine MaterialUpdateAuxVars(Material,comm1,vec_loc,time_level,time)
  ! 
  ! Updates material aux var variables for use in reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/09
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Communicator_Base_module
  use Variables_module, only : POROSITY

  implicit none
  
  type(material_type) :: Material
  class(communicator_type) :: comm1
  Vec :: vec_loc
  PetscReal :: time
  PetscInt :: time_level
  
  select case(time_level)
    case(TIME_T)
      Material%time_t = time
    case(TIME_TpDT)
      Material%time_tpdt = time
  end select  
  
  print *, 'MaterialUpdateAuxVars not implemented.'
  stop
  ! porosity
!  call MaterialGetAuxVarVecLoc(Material,vec_loc,POROSITY,ZERO_INTEGER)
!  call comm1%LocalToLocal(vec_loc,vec_loc)
  ! note that 'time_level' is not ZERO_INTEGER.  thus, this differs
  ! from MaterialAuxVarCommunicate.
!  call MaterialSetAuxVarVecLoc(Material,vec_loc,POROSITY,time_level)

end subroutine MaterialUpdateAuxVars

! ************************************************************************** !

subroutine MaterialAuxVarCommunicate(comm,Material,vec_loc,ivar,isubvar)
  ! 
  ! Sets values of material auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Communicator_Base_module
  
  implicit none

  class(communicator_type), pointer :: comm
  type(material_type) :: Material ! from realization%patch%aux%Material
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  call MaterialGetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)
  call comm%LocalToLocal(vec_loc,vec_loc)
  call MaterialSetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)

end subroutine MaterialAuxVarCommunicate

! ************************************************************************** !

subroutine MaterialUpdatePorosity(Material,global_auxvars,porosity_loc)
  ! 
  ! Gets values of material auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Variables_module
  use Global_Aux_module
  
  implicit none

  type(material_type) :: Material ! from realization%patch%aux%Material
  type(global_auxvar_type) :: global_auxvars(:)
  Vec :: porosity_loc
  
  PetscReal, pointer :: porosity_loc_p(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ghosted_id
  PetscReal :: compressed_porosity
  PetscReal :: dcompressed_porosity_dp
  PetscErrorCode :: ierr
  
  if (soil_compressibility_index > 0) then
    material_auxvars => Material%auxvars
    call VecGetArrayReadF90(porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
    do ghosted_id = 1, Material%num_aux
      material_auxvars(ghosted_id)%porosity = porosity_loc_p(ghosted_id)
      call MaterialCompressSoil(material_auxvars(ghosted_id), &
                                maxval(global_auxvars(ghosted_id)%pres), &
                                compressed_porosity,dcompressed_porosity_dp)
      material_auxvars(ghosted_id)%porosity = compressed_porosity
      material_auxvars(ghosted_id)%dporosity_dp = dcompressed_porosity_dp
    enddo
    call VecRestoreArrayReadF90(porosity_loc,porosity_loc_p, &
                                ierr);CHKERRQ(ierr)
  endif
  
end subroutine MaterialUpdatePorosity

! **************************************************************************** !

subroutine MaterialPropInputRecord(material_property_list)
  ! 
  ! Prints ingested material property information to the input record file
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/08/2016
  ! 

  implicit none

  type(material_property_type), pointer :: material_property_list
  
  type(material_property_type), pointer :: cur_matprop
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'MATERIAL PROPERTIES'
  
  cur_matprop => material_property_list
  do
    if (.not.associated(cur_matprop)) exit
    
    write(id,'(a29)',advance='no') 'material property name: '
    write(id,'(a)') adjustl(trim(cur_matprop%name))
    
    if (Initialized(cur_matprop%external_id)) then
      write(id,'(a29)',advance='no') 'material id: '
      write(word1,*) cur_matprop%external_id
      write(id,'(a)') adjustl(trim(word1))
    endif
    
    write(id,'(a29)',advance='no') 'material property is: '
    if (cur_matprop%active) then
      write(id,'(a)') 'active'
    else
      write(id,'(a)') 'inactive'
    endif
    
    write(id,'(a29)',advance='no') 'permeability: '
    if (associated(cur_matprop%permeability_dataset)) then
      write(id,'(a)') cur_matprop%permeability_dataset%name
      write(id,'(a29)',advance='no') 'from file: '
      write(id,'(a)') cur_matprop%permeability_dataset%filename
    else
      if (cur_matprop%isotropic_permeability) then
        write(id,'(a)') 'isotropic'
      else      
        write(id,'(a)') 'anisotropic'
        if (Initialized(cur_matprop%vertical_anisotropy_ratio)) then
          write(id,'(a29)',advance='no') 'vertical anisotropy ratio: '
          write(word1,*) cur_matprop%vertical_anisotropy_ratio
          write(id,'(a)') adjustl(trim(word1)) 
        endif
      endif
      write(id,'(a29)',advance='no') 'k_xx: '
      write(word1,*) cur_matprop%permeability(1,1)
      write(id,'(a)') adjustl(trim(word1)) // ' m^2'
      write(id,'(a29)',advance='no') 'k_yy: '
      write(word1,*) cur_matprop%permeability(2,2)
      write(id,'(a)') adjustl(trim(word1)) // ' m^2'
      write(id,'(a29)',advance='no') 'k_zz: '
      write(word1,*) cur_matprop%permeability(3,3)
      write(id,'(a)') adjustl(trim(word1)) // ' m^2'
    endif
    if (cur_matprop%permeability_scaling_factor > 0.d0) then
      write(id,'(a29)',advance='no') 'permeability scaling factor: '
      write(word1,*) cur_matprop%permeability_scaling_factor
      write(id,'(a)') adjustl(trim(word1)) 
    endif
    if (.not. Equal(cur_matprop%permeability_pwr,1.d0)) then
      write(id,'(a29)',advance='no') 'permeability power: '
      write(word1,*) cur_matprop%permeability_pwr
      write(id,'(a)') adjustl(trim(word1))
    endif
    if (cur_matprop%permeability_crit_por > 0.d0) then
      write(id,'(a29)',advance='no') 'permeability critical por.: '
      write(word1,*) cur_matprop%permeability_crit_por
      write(id,'(a)') adjustl(trim(word1))
    endif
    
    write(id,'(a29)',advance='no') 'tortuosity: '
    write(word1,*) cur_matprop%tortuosity
    write(id,'(a)') adjustl(trim(word1))
    
    if (Initialized(cur_matprop%rock_density)) then
      write(id,'(a29)',advance='no') 'rock density: '
      write(word1,*) cur_matprop%rock_density
      write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    endif
    
    write(id,'(a29)',advance='no') 'porosity: '
    if (associated(cur_matprop%porosity_dataset)) then
      write(id,'(a)') adjustl(trim(cur_matprop%porosity_dataset%name))
      write(id,'(a29)',advance='no') 'from file: '
      write(id,'(a)') adjustl(trim(cur_matprop%porosity_dataset%filename))
    else
      write(word1,*) cur_matprop%porosity
      write(id,'(a)') adjustl(trim(word1))
    endif
    
    write(id,'(a29)',advance='no') 'tortuosity: '
    if (associated(cur_matprop%tortuosity_dataset)) then
      write(id,'(a)') adjustl(trim(cur_matprop%tortuosity_dataset%name))
      write(id,'(a29)',advance='no') 'from file: '
      write(id,'(a)') adjustl(trim(cur_matprop%tortuosity_dataset%filename))
    else
      write(word1,*) cur_matprop%tortuosity
      write(id,'(a)') adjustl(trim(word1))
    endif

    if (Initialized(cur_matprop%specific_heat)) then
      write(id,'(a29)',advance='no') 'specific heat capacity: '
      write(word1,*) cur_matprop%specific_heat
      write(id,'(a)') adjustl(trim(word1)) // ' J/kg-C'
    endif
    
    if (Initialized(cur_matprop%thermal_conductivity_dry)) then
      write(id,'(a29)',advance='no') 'dry th. conductivity: '
      write(word1,*) cur_matprop%thermal_conductivity_dry
      write(id,'(a)') adjustl(trim(word1)) // ' W/m-C'
    endif
    if (Initialized(cur_matprop%thermal_conductivity_wet)) then
      write(id,'(a29)',advance='no') 'wet th. conductivity: '
      write(word1,*) cur_matprop%thermal_conductivity_wet
      write(id,'(a)') adjustl(trim(word1)) // ' W/m-C'
    endif
    if (cur_matprop%thermal_conductivity_frozen > 0.d0) then
      write(id,'(a29)',advance='no') 'frozen th. conductivity: '
      write(word1,*) cur_matprop%thermal_conductivity_frozen
      write(id,'(a)') adjustl(trim(word1)) // ' W/m-C'
    endif
    
    if (len_trim(cur_matprop%soil_compressibility_function) > 0) then
      write(id,'(a29)',advance='no') 'soil compressibility func.: '
      write(id,'(a)') adjustl(trim(cur_matprop%soil_compressibility_function)) 
    endif
    if (Initialized(cur_matprop%soil_compressibility)) then
      write(id,'(a29)',advance='no') 'soil compressibility: '
      write(word1,*) cur_matprop%soil_compressibility
      write(id,'(a)') adjustl(trim(word1)) 
    endif
    if (Initialized(cur_matprop%soil_reference_pressure)) then
      write(id,'(a29)',advance='no') 'soil reference pressure: '
      write(word1,*) cur_matprop%soil_reference_pressure
      write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    endif
    if (cur_matprop%soil_reference_pressure_initial) then
      write(id,'(a29)',advance='no') 'soil reference pressure: '
      write(id,'(a)') 'initial cell pressure'
    endif
    
    if (cur_matprop%dispersivity(1) > 0.d0 .or. &
        cur_matprop%dispersivity(2) > 0.d0 .or. &
        cur_matprop%dispersivity(3) > 0.d0) then
      write(id,'(a29)',advance='no') 'longitudinal dispersivity: '
      write(word1,*) cur_matprop%dispersivity(1)
      write(id,'(a)') adjustl(trim(word1)) // ' m'
      write(id,'(a29)',advance='no') 'transverse h dispersivity: '
      write(word1,*) cur_matprop%dispersivity(2)
      write(id,'(a)') adjustl(trim(word1)) // ' m'
      write(id,'(a29)',advance='no') 'transverse v dispersivity: '
      write(word1,*) cur_matprop%dispersivity(2)
      write(id,'(a)') adjustl(trim(word1)) // ' m'
    endif
    
    write(id,'(a29)',advance='no') 'cc / saturation function: '
    write(id,'(a)') adjustl(trim(cur_matprop%saturation_function_name))
    
    write(id,'(a29)') '---------------------------: '
    cur_matprop => cur_matprop%next
  enddo
  
end subroutine MaterialPropInputRecord
  
! ************************************************************************** !

recursive subroutine MaterialPropertyDestroy(material_property)
  ! 
  ! Destroys a material_property
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  implicit none
  
  type(material_property_type), pointer :: material_property
  
  if (.not.associated(material_property)) return
  
  call MaterialPropertyDestroy(material_property%next)
  
  ! simply nullify since the datasets reside in a list within realization
  nullify(material_property%permeability_dataset)
  nullify(material_property%permeability_dataset_y)
  nullify(material_property%permeability_dataset_z)
  nullify(material_property%porosity_dataset)
  nullify(material_property%tortuosity_dataset)
  nullify(material_property%compressibility_dataset)
  nullify(material_property%soil_reference_pressure_dataset)
    
  deallocate(material_property)
  nullify(material_property)
  
end subroutine MaterialPropertyDestroy

end module Material_module
