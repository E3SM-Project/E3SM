module Reaction_module

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module  
  use Global_Aux_module
  use Material_Aux_class
  
  use Reaction_Mineral_module
  use Reaction_Immobile_module
  use Reaction_Gas_module

  use Reaction_Mineral_Aux_module
  use Reaction_Immobile_Aux_module
  use Reaction_Gas_Aux_module

  use Reaction_Sandbox_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
 
  private

#include "petsc/finclude/petscsys.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: ReactionInit, &
            ReactionReadPass1, &
            ReactionReadPass2, &
            ReactionReadOutput, &
            ReactionReadRedoxSpecies, &
            RTotal, &
            RTotalAqueous, &
            RActivityCoefficients, &
            RReaction, &
            RReactionDerivative, &
            ReactionProcessConstraint, &
            ReactionEquilibrateConstraint, &
            ReactionPrintConstraint, &
            ReactionFitLogKCoef, &
            ReactionInitializeLogK, &
            ReactionComputeKd, &
            RAccumulationSorb, &
            RAccumulationSorbDerivative, &
            RJumpStartKineticSorption, &
            RAge, &
            RReact, &
            RTAuxVarCompute, &
            RTAccumulation, &
            RTAccumulationDerivative, &
            RTPrintAuxVar, &
            RTSetPlotVariables, &
            ReactionInterpolateLogK_hpt, &
            ReactionInitializeLogK_hpt, &
            RUpdateKineticState, &
            RUpdateTempDependentCoefs, &
            RZeroSorb, &
            RCO2MoleFraction

contains

! ************************************************************************** !

subroutine ReactionInit(reaction,input,option)
  ! 
  ! ReactionReadPass1: Initializes the reaction object, creating object and
  ! reading first pass of CHEMISTRY input file block
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/03/13
  ! 

  use Option_module
  use Input_Aux_module

  implicit none
  
  type(reaction_type), pointer :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option
  
  reaction => ReactionCreate()
  
  ! must be called prior to the first pass
  call RSandboxInit(option)
 
  call ReactionReadPass1(reaction,input,option)
  reaction%primary_species_names => GetPrimarySpeciesNames(reaction)
  ! PCL add in colloid dofs
  option%ntrandof = GetPrimarySpeciesCount(reaction)
  option%ntrandof = option%ntrandof + GetColloidCount(reaction)
  option%ntrandof = option%ntrandof + GetImmobileCount(reaction)
  reaction%ncomp = option%ntrandof  
  reaction%nphase = option%transport%nphase

end subroutine ReactionInit

! ************************************************************************** !

subroutine ReactionReadPass1(reaction,input,option)
  ! 
  ! Reads chemistry (first pass)
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/02/08
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module
  use Variables_module, only : PRIMARY_MOLALITY, PRIMARY_MOLARITY, &
                               TOTAL_MOLALITY, TOTAL_MOLARITY, &
                               SECONDARY_MOLALITY, SECONDARY_MOLARITY
  use Generic_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXWORDLENGTH) :: kd_units
  type(aq_species_type), pointer :: species, prev_species
  type(gas_species_type), pointer :: gas, prev_gas
  type(immobile_species_type), pointer :: immobile_species
  type(immobile_species_type), pointer :: prev_immobile_species
  type(colloid_type), pointer :: colloid, prev_colloid
  type(ion_exchange_rxn_type), pointer :: ionx_rxn, prev_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cation, prev_cation
  type(general_rxn_type), pointer :: general_rxn, prev_general_rxn
  type(radioactive_decay_rxn_type), pointer :: radioactive_decay_rxn
  type(radioactive_decay_rxn_type), pointer :: prev_radioactive_decay_rxn
  type(kd_rxn_type), pointer :: kd_rxn, prev_kd_rxn
  type(kd_rxn_type), pointer :: sec_cont_kd_rxn, sec_cont_prev_kd_rxn
  type(generic_parameter_type), pointer :: generic_list
  PetscInt :: i, temp_int
  PetscReal :: temp_real
  PetscBool :: found
  PetscBool :: reaction_sandbox_read
  PetscBool :: reaction_clm_read

  nullify(prev_species)
  nullify(prev_gas)
  nullify(prev_immobile_species)
  nullify(prev_colloid)
  nullify(prev_cation)
  nullify(prev_general_rxn)
  nullify(prev_radioactive_decay_rxn)
  nullify(prev_kd_rxn)
  nullify(prev_ionx_rxn)
  
  if (option%use_mc) then
    nullify(sec_cont_prev_kd_rxn)
  endif
  
  reaction_sandbox_read = PETSC_FALSE

  kd_units = ''
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY')
    call StringToUpper(word)
    
    select case(trim(word))
    
      case('PRIMARY_SPECIES')
        nullify(prev_species)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          reaction%naqcomp = reaction%naqcomp + 1
          
          species => AqueousSpeciesCreate()
          call InputReadWord(input,option,species%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,&
                             &PRIMARY_SPECIES')    
          if (.not.associated(reaction%primary_species_list)) then
            reaction%primary_species_list => species
            species%id = 1
          endif
          if (associated(prev_species)) then
            prev_species%next => species
            species%id = prev_species%id + 1
          endif
          prev_species => species
          nullify(species)
        enddo
      case('SECONDARY_SPECIES')
        nullify(prev_species)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          reaction%neqcplx = reaction%neqcplx + 1
          
          species => AqueousSpeciesCreate()
          call InputReadWord(input,option,species%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,&
                             &SECONDARY_SPECIES')
          if (.not.associated(reaction%secondary_species_list)) then
            reaction%secondary_species_list => species
            species%id = 1
          endif
          if (associated(prev_species)) then
            prev_species%next => species
            species%id = prev_species%id + 1            
          endif
          prev_species => species
          nullify(species)
        enddo
      case('AQUEOUS_DIFFUSION_COEFFICIENTS','GAS_DIFFUSION_COEFFICIENTS')
        card = word
        select case(card)
          case('AQUEOUS_DIFFUSION_COEFFICIENTS')
            error_string = 'CHEMISTRY,AQUEOUS_DIFFUSION_COEFFICIENTS'
          case('GAS_DIFFUSION_COEFFICIENTS')
            error_string = 'CHEMISTRY,GAS_DIFFUSION_COEFFICIENTS'
        end select
        nullify(generic_list)
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)  
          call InputErrorMsg(input,option,'name',error_string)
          call InputReadDouble(input,option,temp_real)
          call InputErrorMsg(input,option, &
                             trim(word)//' coefficient',error_string)
          call InputReadAndConvertUnits(input,temp_real,'m^2/sec', &
                                        error_string//','//trim(name),option)
          call GenericParameterAddToList(generic_list, &
                 GenericParameterCreate(word,UNINITIALIZED_INTEGER,temp_real))
        enddo
        select case(card)
          case('AQUEOUS_DIFFUSION_COEFFICIENTS')
            reaction%aq_diffusion_coefficients => generic_list
          case('GAS_DIFFUSION_COEFFICIENTS')
            reaction%gas_diffusion_coefficients => generic_list
        end select
      case('ACTIVE_GAS_SPECIES')
        string = 'CHEMISTRY,ACTIVE_GAS_SPECIES'
        call RGasRead(reaction%gas%list,ACTIVE_GAS,string,input,option)
      !TODO(geh): remove GAS_SPECIES
      case('PASSIVE_GAS_SPECIES','GAS_SPECIES')
        string = 'CHEMISTRY,PASSIVE_GAS_SPECIES'
        call RGasRead(reaction%gas%list,PASSIVE_GAS,string,input,option)
      case('IMMOBILE_SPECIES')
        call ImmobileRead(reaction%immobile,input,option)
      case('IMMOBILE_DECAY_REACTION')
        call ImmobileDecayRxnRead(reaction%immobile,input,option)
      case('RADIOACTIVE_DECAY_REACTION')
        reaction%nradiodecay_rxn = reaction%nradiodecay_rxn + 1
        radioactive_decay_rxn => RadioactiveDecayRxnCreate()
        radioactive_decay_rxn%rate_constant = UNINITIALIZED_DOUBLE
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,RADIOACTIVE_DECAY_REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('REACTION')
              ! remainder of string should be the reaction equation
              radioactive_decay_rxn%reaction = trim(adjustl(input%buf))
              ! set flag for error message
              if (len_trim(radioactive_decay_rxn%reaction) < 2) input%ierr = 1
              call InputErrorMsg(input,option,'reaction', &
                               'CHEMISTRY,RADIOACTIVE_DECAY_REACTION,REACTION') 
            case('RATE_CONSTANT')
              call InputReadDouble(input,option, &
                                   radioactive_decay_rxn%rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                'CHEMISTRY,RADIOACTIVE_DECAY_REACTION,RATE_CONSTANT') 
              call InputReadAndConvertUnits(input, &
                     radioactive_decay_rxn%rate_constant,'unitless/sec', &
                    'CHEMISTRY,RADIOACTIVE_DECAY_REACTION,RATE_CONSTANT',option)
            case('HALF_LIFE')
              call InputReadDouble(input,option, &
                                   radioactive_decay_rxn%half_life)
              call InputErrorMsg(input,option,'half life', &
                'CHEMISTRY,RADIOACTIVE_DECAY_REACTION,HALF_LIFE') 
              call InputReadAndConvertUnits(input, &
                     radioactive_decay_rxn%half_life,'sec', &
                    'CHEMISTRY,RADIOACTIVE_DECAY_REACTION,HALF_LIFE',option)
              ! convert half life to rate constant
              radioactive_decay_rxn%rate_constant = &
                -1.d0*log(0.5d0)/radioactive_decay_rxn%half_life
            case default
              call InputKeywordUnrecognized(word, &
                'CHEMISTRY,IMMOBILE_DECAY_REACTION',option)
          end select
        enddo   
        if (Uninitialized(radioactive_decay_rxn%rate_constant)) then
          option%io_buffer = 'RATE_CONSTANT or HALF_LIFE must be set in ' // &
            'RADIOACTIVE_DECAY_REACTION.'
          call printErrMsg(option)
        endif
        if (.not.associated(reaction%radioactive_decay_rxn_list)) then
          reaction%radioactive_decay_rxn_list => radioactive_decay_rxn
          radioactive_decay_rxn%id = 1
        endif
        if (associated(prev_radioactive_decay_rxn)) then
          prev_radioactive_decay_rxn%next => radioactive_decay_rxn
          radioactive_decay_rxn%id = prev_radioactive_decay_rxn%id + 1
        endif
        prev_radioactive_decay_rxn => radioactive_decay_rxn
        nullify(radioactive_decay_rxn)
      case('GENERAL_REACTION')
        reaction%ngeneral_rxn = reaction%ngeneral_rxn + 1
        general_rxn => GeneralRxnCreate()
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,GENERAL_REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('REACTION')
              ! remainder of string should be the reaction equation
              general_rxn%reaction = trim(adjustl(input%buf))
              ! set flag for error message
              if (len_trim(general_rxn%reaction) < 2) input%ierr = 1
              call InputErrorMsg(input,option,'reaction', &
                                 'CHEMISTRY,GENERAL_REACTION,REACTION') 
! For now, the reactants are negative stoich, products positive in reaction equation - geh
#if 0
            case('FORWARD_SPECIES')
              nullify(prev_species)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                
                species => AqueousSpeciesCreate()
                call InputReadWord(input,option,species%name,PETSC_TRUE)  
                call InputErrorMsg(input,option,'keyword','CHEMISTRY, &
                                   GENERAL_REACTION,FORWARD_SPECIES')    
                if (.not.associated(general_rxn%forward_species_list)) then
                  general_rxn%forward_species_list => species
                  species%id = 1
                endif
                if (associated(prev_species)) then
                  prev_species%next => species
                  species%id = prev_species%id + 1
                endif
                prev_species => species
                nullify(species)
              enddo
              nullify(prev_species)
            case('BACKWARD_SPECIES')
              nullify(prev_species)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit
                
                species => AqueousSpeciesCreate()
                call InputReadWord(input,option,species%name,PETSC_TRUE)  
                call InputErrorMsg(input,option,'keyword','CHEMISTRY, &
                                   GENERAL_REACTION,BACKWARD_SPECIES')    
                if (.not.associated(general_rxn%backward_species_list)) then
                  general_rxn%forward_species_list => species
                  species%id = 1
                endif
                if (associated(prev_species)) then
                  prev_species%next => species
                  species%id = prev_species%id + 1
                endif
                prev_species => species
                nullify(species)
              enddo
              nullify(prev_species)
#endif                
            case('FORWARD_RATE')
              call InputReadDouble(input,option,general_rxn%forward_rate)  
              call InputErrorMsg(input,option,'forward rate', &
                                 'CHEMISTRY,GENERAL_REACTION') 
            case('BACKWARD_RATE')
              call InputReadDouble(input,option,general_rxn%backward_rate)  
              call InputErrorMsg(input,option,'backward rate', &
                                 'CHEMISTRY,GENERAL_REACTION') 
          end select
        enddo   
        if (.not.associated(reaction%general_rxn_list)) then
          reaction%general_rxn_list => general_rxn
          general_rxn%id = 1
        endif
        if (associated(prev_general_rxn)) then
          prev_general_rxn%next => general_rxn
          general_rxn%id = prev_general_rxn%id + 1
        endif
        prev_general_rxn => general_rxn
        nullify(general_rxn)

      case('REACTION_SANDBOX')
        call RSandboxRead(input,option)
        reaction_sandbox_read = PETSC_TRUE
      case('MINERALS')
        call MineralRead(reaction%mineral,input,option)
      case('MINERAL_KINETICS') ! mineral kinetics read on second round
        !geh: but we need to count the number of kinetic minerals this round
        temp_int = 0 ! used to count kinetic minerals
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,name,PETSC_TRUE)
          call InputErrorMsg(input,option,name,'CHEMISTRY,MINERAL_KINETICS')
          temp_int = temp_int + 1

          do
            call InputReadPflotranString(input,option)
            call InputReadStringErrorMsg(input,option,card)
            if (InputCheckExit(input,option)) exit
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword', &
                                    'CHEMISTRY,MINERAL_KINETICS')
            call StringToUpper(word)
            select case(word)
              case('PREFACTOR')
                do 
                  call InputReadPflotranString(input,option)
                  call InputReadStringErrorMsg(input,option,card)
                  if (InputCheckExit(input,option)) exit
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'keyword', &
                                      'CHEMISTRY,MINERAL_KINETICS,PREFACTOR')
                  call StringToUpper(word)
                  select case(word)
                    case('PREFACTOR_SPECIES')
                      call InputSkipToEnd(input,option,word)
                  end select
                enddo
            end select
          enddo
        enddo
        reaction%mineral%nkinmnrl = reaction%mineral%nkinmnrl + temp_int

      case('COLLOIDS')
        nullify(prev_colloid)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          
          reaction%ncoll = reaction%ncoll + 1
          
          colloid => ColloidCreate()
          call InputReadWord(input,option,colloid%name,PETSC_TRUE)  
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,COLLOIDS')    
          call InputReadDouble(input,option,colloid%mobile_fraction)  
          call InputDefaultMsg(input,option,'CHEMISTRY,COLLOIDS,MOBILE_FRACTION')          
          if (.not.associated(reaction%colloid_list)) then
            reaction%colloid_list => colloid
            colloid%id = 1
          endif
          if (associated(prev_colloid)) then
            prev_colloid%next => colloid
            colloid%id = prev_colloid%id + 1
          endif
          prev_colloid => colloid
          nullify(colloid)
        enddo
      case('SORPTION')
!geh        nullify(prev_srfcplx_rxn)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CHEMISTRY,SORPTION')
          call StringToUpper(word)   

          select case(trim(word))

            case('ISOTHERM_REACTIONS')
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                reaction%neqkdrxn = reaction%neqkdrxn + 1

                kd_rxn => KDRxnCreate()
                if (option%use_mc) then
                  sec_cont_kd_rxn => KDRxnCreate()
                endif
                ! first string is species name
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'species name', &
                                   'CHEMISTRY,ISOTHERM_REACTIONS')
                kd_rxn%species_name = trim(word)
                if (option%use_mc) then
                  sec_cont_kd_rxn%species_name = kd_rxn%species_name
                endif
                do 
                  call InputReadPflotranString(input,option)
                  if (InputError(input)) exit
                  if (InputCheckExit(input,option)) exit

                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'keyword', &
                                     'CHEMISTRY,ISOTHERM_REACTIONS')
                  call StringToUpper(word)
                  
                  ! default type is linear
                  kd_rxn%itype = SORPTION_LINEAR
                  select case(trim(word))
                    case('TYPE')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'type', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                      select case(word)
                        case('LINEAR')
                          kd_rxn%itype = SORPTION_LINEAR
                        case('LANGMUIR')
                          kd_rxn%itype = SORPTION_LANGMUIR
                        case('FREUNDLICH')
                          kd_rxn%itype = SORPTION_FREUNDLICH
                        case default
                          call InputKeywordUnrecognized(word, &
                                'CHEMISTRY,SORPTION,ISOTHERM_REACTIONS,TYPE', &
                                option)
                      end select
                      if (option%use_mc) then
                        sec_cont_kd_rxn%itype = kd_rxn%itype
                      endif
                    case('DISTRIBUTION_COEFFICIENT','KD')
                      call InputReadDouble(input,option,kd_rxn%Kd)
                      call InputErrorMsg(input,option, &
                                         'DISTRIBUTION_COEFFICIENT', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (input%ierr == 0) kd_units = trim(word)
                    ! S.Karra, 02/20/2014
                    case('SEC_CONT_DISTRIBUTION_COEFFICIENT', &
                         'SEC_CONT_KD')
                         if (.not.option%use_mc) then
                           option%io_buffer = 'Make sure MULTIPLE_CONTINUUM ' &
                                   // 'keyword is set, SECONDARY_CONTINUUM_KD.'
                           call printErrMsg(option)
                         else
                           call InputReadDouble(input,option,sec_cont_kd_rxn%Kd)
                           call InputErrorMsg(input,option, &
                             'SECONDARY_CONTINUUM_DISTRIBUTION_COEFFICIENT', &
                             'CHEMISTRY,ISOTHERM_REACTIONS')                    
                        endif
                    case('LANGMUIR_B')
                      call InputReadDouble(input,option,kd_rxn%Langmuir_B)
                      call InputErrorMsg(input,option,'Langmuir_B', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                      kd_rxn%itype = SORPTION_LANGMUIR
                    case('FREUNDLICH_N')
                      call InputReadDouble(input,option,kd_rxn%Freundlich_N)
                      call InputErrorMsg(input,option,'Freundlich_N', &
                                         'CHEMISTRY,ISOTHERM_REACTIONS')
                      kd_rxn%itype = SORPTION_FREUNDLICH
                    case('KD_MINERAL_NAME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'KD_MINERAL_NAME', &
                                         'ISOTHERM_REACTIONS,KD_MINERAL_NAME')
                      kd_rxn%kd_mineral_name = word                      
                    case default
                      call InputKeywordUnrecognized(word, &
                              'CHEMISTRY,SORPTION,ISOTHERM_REACTIONS',option)
                  end select
                enddo

                if (len_trim(kd_units) > 0) then
                  if (len_trim(kd_rxn%kd_mineral_name) > 0) then
                    internal_units = 'L/kg'
                    kd_rxn%Kd = kd_rxn%Kd * &
                      UnitsConvertToInternal(kd_units,internal_units,option)
                  else
                    internal_units = 'kg/m^3'
                    kd_rxn%Kd = kd_rxn%Kd * &
                      UnitsConvertToInternal(kd_units,internal_units,option)
                  endif
                endif

                ! add to list
                if (.not.associated(reaction%kd_rxn_list)) then
                  reaction%kd_rxn_list => kd_rxn
                  kd_rxn%id = 1
                endif
                if (associated(prev_kd_rxn)) then
                  prev_kd_rxn%next => kd_rxn
                  kd_rxn%id = prev_kd_rxn%id + 1
                endif
                prev_kd_rxn => kd_rxn
                nullify(kd_rxn)
                
                if (option%use_mc) then
                ! add to list
                  if (.not.associated(reaction%sec_cont_kd_rxn_list)) then
                    reaction%sec_cont_kd_rxn_list => sec_cont_kd_rxn
                    sec_cont_kd_rxn%id = 1
                  endif
                  if (associated(sec_cont_prev_kd_rxn)) then
                    sec_cont_prev_kd_rxn%next => sec_cont_kd_rxn
                    sec_cont_kd_rxn%id = sec_cont_prev_kd_rxn%id + 1
                  endif
                  sec_cont_prev_kd_rxn => sec_cont_kd_rxn
                  nullify(sec_cont_kd_rxn)
                endif
              enddo
            
            case('ION_EXCHANGE_RXN')
              ionx_rxn => IonExchangeRxnCreate()
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'keyword', &
                                   'CHEMISTRY,ION_EXCHANGE_RXN')
                call StringToUpper(word)
                
                select case(trim(word))
                  case('MINERAL')
                    call InputReadWord(input,option,ionx_rxn%mineral_name, &
                                       PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword', &
                      'CHEMISTRY,ION_EXCHANGE_RXN,MINERAL_NAME')
                  case('CEC')
                    call InputReadDouble(input,option,ionx_rxn%CEC)
                    call InputErrorMsg(input,option,'keyword', &
                                       'CHEMISTRY,ION_EXCHANGE_RXN,CEC')
                  case('CATIONS')
                    string = '' ! string denotes the reference cation 
                    nullify(prev_cation)
                    do
                      call InputReadPflotranString(input,option)
                      if (InputError(input)) exit
                      if (InputCheckExit(input,option)) exit
                      
                      cation => IonExchangeCationCreate()
                      reaction%neqionxcation = reaction%neqionxcation + 1
                      call InputReadWord(input,option,cation%name,PETSC_TRUE)
                      call InputErrorMsg(input,option,'keyword', &
                        'CHEMISTRY,ION_EXCHANGE_RXN,CATION_NAME')
                      call InputReadDouble(input,option,cation%k)
                      call InputErrorMsg(input,option,'keyword', &
                                         'CHEMISTRY,ION_EXCHANGE_RXN,K')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (input%ierr == 0) then
                        if (StringCompareIgnoreCase(word,'REFERENCE')) then
                          string = cation%name
                        else
                          call InputKeywordUnrecognized(word, &
                                  'CHEMISTRY,ION_EXCHANGE_RXN,CATIONS',option)
                        endif
                      endif
                      if (.not.associated(ionx_rxn%cation_list)) then
                        ionx_rxn%cation_list => cation
                      endif
                      if (associated(prev_cation)) then
                        prev_cation%next => cation
                      endif
                      prev_cation => cation
                      nullify(cation)
                    enddo
                    if (len_trim(string) == 0) then
                      option%io_buffer = &
                        'Reference cation missing in Ion Exchange reaction.'
                      call printErrMsg(option)
                    endif
                    cation => ionx_rxn%cation_list
                    nullify(prev_cation)
                    do
                      if (.not.associated(cation)) exit
                      if (StringCompare(cation%name,string)) then
                        if (dabs(cation%k - 1.d0) > 1.d-40) then
                          option%io_buffer = 'Reference cation "' // &
                            trim(cation%name) // '" must have k = 1.d0.'
                          call printErrMsg(option)
                        endif
                        ! move it to the front of the list
                        if (associated(prev_cation)) then
                          prev_cation%next => cation%next
                          cation%next => ionx_rxn%cation_list
                          ionx_rxn%cation_list => cation
                        else
                          ! nothing to do as it is at the front of the list
                        endif
                        exit
                      endif
                      prev_cation => cation
                      cation => cation%next
                    enddo
                  case default
                    call InputKeywordUnrecognized(word, &
                              'CHEMISTRY,ION_EXCHANGE_RXN',option)
                end select
              enddo
              if (.not.associated(reaction%ion_exchange_rxn_list)) then
                reaction%ion_exchange_rxn_list => ionx_rxn
                ionx_rxn%id = 1
              endif
              if (associated(prev_ionx_rxn)) then
                prev_ionx_rxn%next => ionx_rxn
                ionx_rxn%id = prev_ionx_rxn%id + 1
              endif
              prev_ionx_rxn => ionx_rxn

              reaction%neqionxrxn = ionx_rxn%id

              nullify(ionx_rxn)          
            case('JUMPSTART_KINETIC_SORPTION')
              option%transport%jumpstart_kinetic_sorption = PETSC_TRUE
              option%transport%no_restart_kinetic_sorption = PETSC_TRUE
            case('NO_CHECKPOINT_KINETIC_SORPTION')
              option%transport%no_checkpoint_kinetic_sorption = PETSC_TRUE
            case('NO_RESTART_KINETIC_SORPTION')
              option%transport%no_restart_kinetic_sorption = PETSC_TRUE
          end select
        enddo
      case('DATABASE')
        call InputReadNChars(input,option,reaction%database_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)  
        call InputErrorMsg(input,option,'keyword', &
                           'CHEMISTRY,DATABASE FILENAME')  
      case('LOG_FORMULATION')
        reaction%use_log_formulation = PETSC_TRUE
      case('TRUNCATE_CONCENTRATION')
        call InputReadDouble(input,option,reaction%truncated_concentration)
        call InputErrorMsg(input,option,'truncate_concentration','CHEMISTRY')
      case('GEOTHERMAL_HPT')
        reaction%use_geothermal_hpt = PETSC_TRUE           
      case('NO_CHECK_UPDATE')
        reaction%check_update = PETSC_FALSE       
      case('NO_RESTART_MINERAL_VOL_FRAC')
        option%transport%no_restart_mineral_vol_frac = PETSC_TRUE
      case('NO_CHECKPOINT_ACT_COEFS')
        reaction%checkpoint_activity_coefs = PETSC_FALSE
      case('ACTIVITY_COEFFICIENTS')
        reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_LAG        
        reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_TIMESTEP        
        do 
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr /= 0) exit
          select case(trim(word))
            case('OFF')
              reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
            case('LAG')
              reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_LAG    
            case('NEWTON')
              reaction%act_coef_update_algorithm = ACT_COEF_ALGORITHM_NEWTON       
            case('TIMESTEP')
              reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_TIMESTEP
            case('NEWTON_ITERATION')
              reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_NEWTON_ITER
            case default
              call InputKeywordUnrecognized(word, &
                        'CHEMISTRY,ACTIVITY_COEFFICIENTS',option)
          end select
        enddo
      case('NO_BDOT')
        reaction%act_coef_use_bdot = PETSC_FALSE
      case('UPDATE_POROSITY')
        reaction%update_porosity = PETSC_TRUE
        option%flow%transient_porosity = PETSC_TRUE
      case('UPDATE_TORTUOSITY')
        reaction%update_tortuosity = PETSC_TRUE
      case('UPDATE_PERMEABILITY')
        reaction%update_permeability = PETSC_TRUE
      case('UPDATE_MINERAL_SURFACE_AREA')
        reaction%update_mineral_surface_area = PETSC_TRUE
      case('UPDATE_MNRL_SURF_AREA_WITH_POR')
        reaction%update_mnrl_surf_with_porosity = PETSC_TRUE
      case('UPDATE_ARMOR_MINERAL_SURFACE')
        reaction%update_armor_mineral_surface = PETSC_TRUE
      case('UPDATE_ARMOR_MINERAL_SURFACE_FLAG')
        reaction%update_armor_mineral_surface = PETSC_TRUE
      case('MOLAL','MOLALITY')
        reaction%initialize_with_molality = PETSC_TRUE
      case('ACTIVITY_H2O','ACTIVITY_WATER')
        reaction%use_activity_h2o = PETSC_TRUE
      case('REDOX_SPECIES')
        call InputSkipToEnd(input,option,word)
      case('OUTPUT')
        call InputSkipToEnd(input,option,word)
      case('MAX_DLNC')
        call InputReadDouble(input,option,reaction%max_dlnC)
        call InputErrorMsg(input,option,trim(word),'CHEMISTRY')
      case('OPERATOR_SPLIT','OPERATOR_SPLITTING')
        option%io_buffer = 'OPERATOR_SPLIT functionality has not been ' // &
          'reimplemented in the refactored PFLOTRAN at this time.  Please ' // &
          'GLOBAL_IMPLICIT (remove OPERATOR_SPLIT(TING)) or ask for the ' // &
          'capability through pflotran-dev@googlegroups.com if you really ' // &
          'need it.'
        call printErrMsg(option)
        option%transport%reactive_transport_coupling = OPERATOR_SPLIT    
      case('EXPLICIT_ADVECTION')
        option%itranmode = EXPLICIT_ADVECTION
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          call StringToUpper(word)
          select case(word)
            !TODO(geh): fix these hardwired values.
            case('UPWIND')
              option%transport%tvd_flux_limiter = 1
            case('MINMOD')
              option%transport%tvd_flux_limiter = 3
            case('MC')
              option%transport%tvd_flux_limiter = 2
            case('SUPERBEE')
              option%transport%tvd_flux_limiter = 4
            case('VANLEER')
              option%transport%tvd_flux_limiter = 5
            case default
              call InputKeywordUnrecognized(word, &
                        'CHEMISTRY,EXPLICIT_ADVECTION',option)
          end select
          option%io_buffer = 'Flux Limiter: ' // trim(word)
          call printMsg(option)
        else
          call InputDefaultMsg(input,option,'TVD Flux Limiter')
        endif
      case('MAX_RELATIVE_CHANGE_TOLERANCE','REACTION_TOLERANCE')
        call InputReadDouble(input,option,reaction%max_relative_change_tolerance)
        call InputErrorMsg(input,option,'maximum relative change tolerance','CHEMISTRY')
      case('MAX_RESIDUAL_TOLERANCE')
        call InputReadDouble(input,option,reaction%max_residual_tolerance)
        call InputErrorMsg(input,option,'maximum residual tolerance','CHEMISTRY')
      case('MINIMUM_POROSITY')
        call InputReadDouble(input,option,reaction%minimum_porosity)
        call InputErrorMsg(input,option,'minimim porosity','CHEMISTRY')
      case('USE_FULL_GEOCHEMISTRY')
        reaction%use_full_geochemistry = PETSC_TRUE
      !-------------------------------------------------
      case('ONLY_VERTICAL_TRANSPORT')
        option%transport%only_vertical_tran = PETSC_TRUE
      !-------------------------------------------------
      case default
        call InputKeywordUnrecognized(word,'CHEMISTRY',option)
    end select
  enddo

  call GasSpeciesListMergeDuplicates(reaction%gas%list)
  
  reaction%neqsorb = reaction%neqionxrxn + &
                     reaction%neqkdrxn
  reaction%nsorb = reaction%neqsorb

    

  if (reaction%print_free_conc_type == 0) then
    if (reaction%initialize_with_molality) then
      reaction%print_free_conc_type = PRIMARY_MOLALITY
    else
      reaction%print_free_conc_type = PRIMARY_MOLARITY
    endif
  endif
  if (reaction%print_tot_conc_type == 0) then
    if (reaction%initialize_with_molality) then
      reaction%print_tot_conc_type = TOTAL_MOLALITY
    else
      reaction%print_tot_conc_type = TOTAL_MOLARITY
    endif
  endif
  if (reaction%print_secondary_conc_type == 0) then
    if (reaction%initialize_with_molality) then
      reaction%print_secondary_conc_type = SECONDARY_MOLALITY
    else
      reaction%print_secondary_conc_type = SECONDARY_MOLARITY
    endif
  endif
  if (reaction%neqcplx + reaction%nsorb + reaction%mineral%nmnrl + &
      reaction%ngeneral_rxn + &
      reaction%nradiodecay_rxn + reaction%immobile%nimmobile > 0 .or. &
      GasGetCount(reaction%gas%list,ACTIVE_AND_PASSIVE_GAS) > 0 .or. &
      reaction_clm_read .or. &
      reaction_sandbox_read) then
    reaction%use_full_geochemistry = PETSC_TRUE
  endif
      
  ! ensure that update porosity is ON if update of tortuosity, permeability or
  ! mineral surface area are ON
  if (.not.reaction%update_porosity .and. &
      (reaction%update_tortuosity .or. &
       reaction%update_permeability .or. &
       reaction%update_mnrl_surf_with_porosity)) then
    option%io_buffer = 'UPDATE_POROSITY must be listed under CHEMISTRY ' // &
      'card when UPDATE_TORTUOSITY, UPDATE_PERMEABILITY, or ' // &
      'UPDATE_MNRL_SURF_WITH_POR are listed.'
    call printErrMsg(option)
  endif
    
  if (len_trim(reaction%database_filename) < 2) &
    reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
  
end subroutine ReactionReadPass1

! ************************************************************************** !

subroutine ReactionReadPass2(reaction,input,option)
  ! 
  ! Reads chemistry on pass 2
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/03/13
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none

  type(reaction_type) :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
  
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'word','CHEMISTRY') 
    select case(trim(word))
      case('PRIMARY_SPECIES','SECONDARY_SPECIES','GAS_SPECIES', &
            'MINERALS','COLLOIDS','GENERAL_REACTION', &
            'IMMOBILE_SPECIES','RADIOACTIVE_DECAY_REACTION', &
            'IMMOBILE_DECAY_REACTION', &
            'ACTIVE_GAS_SPECIES','PASSIVE_GAS_SPECIES', &
            'AQUEOUS_DIFFUSION_COEFFICIENTS','GAS_DIFFUSION_COEFFICIENTS')
        call InputSkipToEND(input,option,card)
      case('REDOX_SPECIES')
        call ReactionReadRedoxSpecies(reaction,input,option)
      case('OUTPUT')
        call ReactionReadOutput(reaction,input,option)
      case('MINERAL_KINETICS')
        call MineralReadKinetics(reaction%mineral,input,option)
      case('REACTION_SANDBOX')
        call RSandboxSkipInput(input,option)
      case('SORPTION')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'SORPTION','CHEMISTRY') 
          select case(trim(word))
            case('ISOTHERM_REACTIONS')
              do
                call InputReadPflotranString(input,option)
                call InputReadStringErrorMsg(input,option,card)
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,word, &
                                    'CHEMISTRY,SORPTION,ISOTHERM_REACTIONS') 
                ! skip over remaining cards to end of each kd entry
                call InputSkipToEnd(input,option,word)
              enddo
            case('SURFACE_COMPLEXATION_RXN','ION_EXCHANGE_RXN')
              do
                call InputReadPflotranString(input,option)
                call InputReadStringErrorMsg(input,option,card)
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'SORPTION','CHEMISTRY')
                select case(trim(word))
                  case('COMPLEXES','CATIONS')
                    call InputSkipToEND(input,option,word)
                  case('COMPLEX_KINETICS')
                    do
                      call InputReadPflotranString(input,option)
                      call InputReadStringErrorMsg(input,option,card)
                      if (InputCheckExit(input,option)) exit
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,word, &
                             'CHEMISTRY,SURFACE_COMPLEXATION_RXN,KINETIC_RATES')
                      ! skip over remaining cards to end of each mineral entry
                      call InputSkipToEnd(input,option,word)
                    enddo
                end select 
              enddo
            case('NUM_THREADS')
            case('JUMPSTART_KINETIC_SORPTION')
            case('NO_CHECKPOINT_KINETIC_SORPTION')
            case('NO_RESTART_KINETIC_SORPTION')
              ! dummy placeholder
          end select
        enddo
      case('MICROBIAL_REACTION')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'MICROBIAL_REACTION','CHEMISTRY')
          select case(trim(word))
            case('INHIBITION','MONOD','BIOMASS')
              call InputSkipToEND(input,option,word)
          end select 
        enddo
      case('MOLAL','MOLALITY', &
            'UPDATE_POROSITY','UPDATE_TORTUOSITY', &
            'UPDATE_PERMEABILITY','UPDATE_MINERAL_SURFACE_AREA', &
            'NO_RESTART_MINERAL_VOL_FRAC','USE_FULL_GEOCHEMISTRY')
        ! dummy placeholder
    end select
  enddo  
  
end subroutine ReactionReadPass2

! ************************************************************************** !

subroutine ReactionReadRedoxSpecies(reaction,input,option)
  ! 
  ! Reads names of mineral species and sets flag
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/01/11
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Input_Aux_module
  use String_module  
  use Option_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: name
  
  type(aq_species_type), pointer :: cur_species

  input%ierr = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,name,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REDOX_SPECIES')
    
    cur_species => reaction%primary_species_list
    do 
      if (.not.associated(cur_species)) exit
      if (StringCompare(cur_species%name,name,MAXWORDLENGTH)) then
        cur_species%is_redox = PETSC_TRUE
        exit
      endif
      cur_species => cur_species%next
    enddo
    
    if (.not.associated(cur_species)) then
      option%io_buffer = 'Redox species "' // trim(name) // &
        '" not found among primary species.'
      call printErrMsg(option)
    endif
  enddo
  
end subroutine ReactionReadRedoxSpecies

! ************************************************************************** !

subroutine ReactionProcessConstraint(reaction,constraint_name, &
                                     aq_species_constraint, &
                                     free_ion_guess, &
                                     mineral_constraint, &
                                     colloid_constraint, &
                                     immobile_constraint, &
                                     option)
  ! 
  ! Initializes constraints based on primary
  ! species in system
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Transport_Constraint_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(guess_constraint_type), pointer :: free_ion_guess
  type(mineral_constraint_type), pointer :: mineral_constraint

  type(colloid_constraint_type), pointer :: colloid_constraint
  type(immobile_constraint_type), pointer :: immobile_constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: icomp, jcomp
  PetscInt :: icoll, jcoll
  PetscInt :: igas, imnrl
  PetscReal :: constraint_conc(reaction%naqcomp)
  PetscInt :: constraint_type(reaction%naqcomp)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(reaction%naqcomp)
  character(len=MAXWORDLENGTH) :: constraint_colloid_name(reaction%ncoll)
  PetscInt :: constraint_id(reaction%naqcomp)
  PetscBool :: external_dataset(reaction%naqcomp)
  
  constraint_id = 0
  constraint_aux_string = ''
  constraint_type = 0
  constraint_conc = 0.d0
  external_dataset = PETSC_FALSE
  
  ! aqueous species
  do icomp = 1, reaction%naqcomp
    found = PETSC_FALSE
    do jcomp = 1, reaction%naqcomp
      if (StringCompare(aq_species_constraint%names(icomp), &
                        reaction%primary_species_names(jcomp), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
               'Species ' // trim(aq_species_constraint%names(icomp)) // &
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among primary species.'
      call printErrMsg(option)
    else
      constraint_type(jcomp) = aq_species_constraint%constraint_type(icomp)
      constraint_aux_string(jcomp) = &
        aq_species_constraint%constraint_aux_string(icomp)
      constraint_conc(jcomp) = aq_species_constraint%constraint_conc(icomp)
      external_dataset(jcomp) = aq_species_constraint%external_dataset(icomp)
      
      ! link constraint species
      select case(constraint_type(jcomp))
        case(CONSTRAINT_MINERAL)
          found = PETSC_FALSE
          do imnrl = 1, reaction%mineral%nmnrl
            if (StringCompare(constraint_aux_string(jcomp), &
                              reaction%mineral%mineral_names(imnrl), &
                                MAXWORDLENGTH)) then
              constraint_id(jcomp) = imnrl
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            option%io_buffer = 'Constraint mineral: ' // &
                     trim(constraint_aux_string(jcomp)) // &
                     ' for aqueous species: ' // &
                     trim(reaction%primary_species_names(jcomp)) // &
                     ' in constraint: ' // &
                     trim(constraint_name) // ' not found.' 
            call printErrMsg(option)         
          endif
        case(CONSTRAINT_GAS, CONSTRAINT_SUPERCRIT_CO2)
          found = PETSC_FALSE
          do igas = 1, reaction%gas%npassive_gas
            if (StringCompare(constraint_aux_string(jcomp), &
                              reaction%gas%passive_names(igas), &
                              MAXWORDLENGTH)) then
              constraint_id(jcomp) = igas
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (.not.found) then
            option%io_buffer = 'Constraint gas: ' // &
                     trim(constraint_aux_string(jcomp)) // &
                     ' for aqueous species: ' // &
                     trim(reaction%primary_species_names(jcomp)) // &
                     ' in constraint: ' // &
                     trim(constraint_name) // ' not found.' 
            call printErrMsg(option)         
          endif
      end select

    endif
  enddo

  ! ensure that two aqueous species do not share the same mineral or gas
  ! constraint
  do icomp = 1, reaction%naqcomp
    do jcomp = icomp + 1, reaction%naqcomp
      if (constraint_type(icomp) == constraint_type(jcomp) .and. &
          constraint_id(icomp) == constraint_id(jcomp) .and. &
          constraint_id(icomp) > 0) then
        option%io_buffer = 'Two aqueous species (' // &
          trim(reaction%primary_species_names(icomp)) // ' and ' // &
          trim(reaction%primary_species_names(jcomp)) // ') are &
          &constrained by the same'
        select case(constraint_type(icomp))
          case(CONSTRAINT_MINERAL)
            option%io_buffer = trim(option%io_buffer) // ' mineral (' // &
              trim(reaction%mineral%mineral_names(constraint_id(icomp))) // &
              ')'
          case(CONSTRAINT_GAS, CONSTRAINT_SUPERCRIT_CO2)
            option%io_buffer = trim(option%io_buffer) // ' gas (' // &
              trim(reaction%gas%passive_names(constraint_id(icomp))) // &
              ')'
          case default
        end select
        option%io_buffer = trim(option%io_buffer) // ' in CONSTRAINT "' // &
          trim(constraint_name) // '".'
        call printErrMsg(option)
      endif
    enddo
  enddo
  
  ! place ordered constraint parameters back in original arrays
  aq_species_constraint%constraint_type = constraint_type
  aq_species_constraint%constraint_aux_string = constraint_aux_string
  aq_species_constraint%constraint_spec_id = constraint_id
  aq_species_constraint%constraint_conc = constraint_conc
  aq_species_constraint%external_dataset = external_dataset

  
  if (.not.reaction%use_full_geochemistry) return
  
  ! free ion guess
  if (associated(free_ion_guess)) then
    constraint_conc = 0.d0
    do icomp = 1, reaction%naqcomp
      found = PETSC_FALSE
      do jcomp = 1, reaction%naqcomp
        if (StringCompare(free_ion_guess%names(icomp), &
                          reaction%primary_species_names(jcomp), &
                          MAXWORDLENGTH)) then
          found = PETSC_TRUE
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = &
                 'Guess species ' // trim(free_ion_guess%names(icomp)) // &
                 ' from CONSTRAINT ' // trim(constraint_name) // &
                 ' not found among primary species.'
        call printErrMsg(option)
      else
        constraint_conc(jcomp) = free_ion_guess%conc(icomp)
      endif
    enddo
    ! place ordered concentrations back in original array
    free_ion_guess%conc = constraint_conc
  endif
  
  ! minerals
  call MineralProcessConstraint(reaction%mineral,constraint_name, &
                                mineral_constraint,option)

  ! microbial immobile
  call ImmobileProcessConstraint(reaction%immobile,constraint_name, &
                                 immobile_constraint,option)
  
end subroutine ReactionProcessConstraint

! ************************************************************************** !

subroutine ReactionEquilibrateConstraint(rt_auxvar,global_auxvar, &
                                         material_auxvar, &
                                         reaction,constraint_name, &
                                         aq_species_constraint, &
                                         free_ion_guess_constraint, &
                                         mineral_constraint, &
                                         colloid_constraint, &
                                         immobile_constraint, &
                                         num_iterations, &
                                         use_prev_soln_as_guess,option)
  ! 
  ! Equilibrates constraint concentrations
  ! with prescribed geochemistry
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/08
  !
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module  
  use Utility_module
  use Transport_Constraint_module
  use EOS_Water_module
  use Material_Aux_class

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type), pointer :: reaction
  character(len=MAXWORDLENGTH) :: constraint_name
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(guess_constraint_type), pointer :: free_ion_guess_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint

  type(colloid_constraint_type), pointer :: colloid_constraint
  type(immobile_constraint_type), pointer :: immobile_constraint
  PetscInt :: num_iterations
  
  PetscBool :: use_prev_soln_as_guess
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icomp, jcomp, kcomp
  PetscInt :: imnrl, jmnrl
  PetscInt :: icplx
  PetscInt :: irxn, isite, ncplx, k, ikinrxn
  PetscInt :: igas
  PetscReal :: conc(reaction%naqcomp)
  PetscInt :: constraint_type(reaction%naqcomp)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(reaction%naqcomp)

  type(mineral_type), pointer :: mineral_reaction

  PetscReal :: Res(reaction%naqcomp)
  PetscReal :: update(reaction%naqcomp)
  PetscReal :: total_conc(reaction%naqcomp)
  PetscReal :: free_conc(reaction%naqcomp)
  PetscReal :: Jac(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: indices(reaction%naqcomp)
  PetscReal :: norm
  PetscReal :: maximum_residual, maximum_relative_change
  PetscReal :: ratio, min_ratio
  PetscReal :: prev_molal(reaction%naqcomp)
  PetscBool :: compute_activity_coefs

  PetscInt :: constraint_id(reaction%naqcomp)
  PetscReal :: lnQK, QK
  PetscReal :: tempreal
  PetscReal :: pres, tc, xphico2, henry, m_na, m_cl, xmass 
  PetscInt :: comp_id
  PetscReal :: convert_molal_to_molar
  PetscReal :: convert_molar_to_molal
  
  PetscBool :: charge_balance_warning_flag = PETSC_FALSE
  PetscBool :: use_log_formulation
  
  PetscReal :: Jac_num(reaction%naqcomp)
  PetscReal :: Res_pert, pert, prev_value, coh0

  PetscInt :: iphase
  PetscInt :: idof
  PetscInt :: istartaq, iendaq
  PetscInt :: irate

  PetscInt :: num_it_act_coef_turned_on
  
  ! CO2-specific
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
  PetscInt :: iflag
  PetscErrorCode :: ierr

  mineral_reaction => reaction%mineral
    
  constraint_type = aq_species_constraint%constraint_type
  constraint_aux_string = aq_species_constraint%constraint_aux_string
  constraint_id = aq_species_constraint%constraint_spec_id
  conc = aq_species_constraint%constraint_conc

  istartaq = reaction%offset_aqueous
  iendaq = reaction%offset_aqueous + reaction%naqcomp    

  iphase = 1
  
  xmass = 1.d0
  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  
  if (reaction%initialize_with_molality) then
    convert_molal_to_molar = global_auxvar%den_kg(iphase)*xmass/1000.d0
    convert_molar_to_molal = 1.d0
  else
    convert_molal_to_molar = 1.d0
    convert_molar_to_molal = 1000.d0/global_auxvar%den_kg(iphase)/xmass
  endif

  !geh: We do need this setting of mineral volume fractions.  If not specified
  !     by a dataset, we must set the volume fraction and areas so that
  !     surface complexation is scaled correctly.  Yes, these will get 
  !     overwrittent with the same values when called from 
  !     CondControlAssignTranInitCond(), but big deal.
  if (associated(mineral_constraint)) then
    do imnrl = 1, mineral_reaction%nkinmnrl
      ! if read from a dataset, the mineral volume frac has already been set.
      if (.not.mineral_constraint%external_vol_frac_dataset(imnrl)) then
        rt_auxvar%mnrl_volfrac0(imnrl) = &
          mineral_constraint%constraint_vol_frac(imnrl)
        rt_auxvar%mnrl_volfrac(imnrl) = &
          mineral_constraint%constraint_vol_frac(imnrl)
      endif
      if (.not.mineral_constraint%external_area_dataset(imnrl)) then
        rt_auxvar%mnrl_area0(imnrl) = &
          mineral_constraint%constraint_area(imnrl)
        rt_auxvar%mnrl_area(imnrl) = &
          mineral_constraint%constraint_area(imnrl)
      endif
    enddo
  endif

  if (associated(colloid_constraint)) then      
    colloid_constraint%basis_conc_mob = colloid_constraint%constraint_conc_mob        
    colloid_constraint%basis_conc_imb = colloid_constraint%constraint_conc_imb        
    rt_auxvar%colloid%conc_mob = colloid_constraint%basis_conc_mob* &
                                 convert_molar_to_molal
    !TODO(geh): this can't be correct as immobile concentrations are mol/m^3
    rt_auxvar%colloid%conc_imb = colloid_constraint%basis_conc_imb* &
                                 convert_molar_to_molal
  endif  
  
  if (.not.reaction%use_full_geochemistry) then
    aq_species_constraint%basis_molarity = conc ! don't need to convert
    rt_auxvar%pri_molal = aq_species_constraint%basis_molarity* &
                          convert_molar_to_molal
    rt_auxvar%total(:,iphase) = aq_species_constraint%basis_molarity
    return
  endif
  
  if (.not.option%use_isothermal) then
    call RUpdateTempDependentCoefs(global_auxvar,reaction,PETSC_TRUE,option)
  endif
  
  if (use_prev_soln_as_guess) then
    free_conc = rt_auxvar%pri_molal
  else if (associated(free_ion_guess_constraint)) then
    free_conc = free_ion_guess_constraint%conc
  else
    free_conc = 1.d-9
  endif
  total_conc = 0.d0
  do icomp = 1, reaction%naqcomp
    select case(constraint_type(icomp))
      case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
        ! units = mol/L
        total_conc(icomp) = conc(icomp)*convert_molal_to_molar
        ! free_conc guess set above
      case(CONSTRAINT_TOTAL_SORB)
        ! units = mol/m^3 bulk
        total_conc(icomp) = conc(icomp)
      case(CONSTRAINT_FREE)
        free_conc(icomp) = conc(icomp)*convert_molar_to_molal
      case(CONSTRAINT_LOG)
        free_conc(icomp) = (10.d0**conc(icomp))*convert_molar_to_molal
      case(CONSTRAINT_CHARGE_BAL)
        if (.not.use_prev_soln_as_guess) then
          free_conc(icomp) = conc(icomp)*convert_molar_to_molal ! just a guess
        endif
      case(CONSTRAINT_PH)
        ! check if H+ id set
        if (associated(reaction%species_idx)) then
          if (reaction%species_idx%h_ion_id /= 0) then
            ! check if icomp is H+
            if (reaction%species_idx%h_ion_id /= icomp) then
              string = 'OH-'
              if (.not.StringCompare(reaction%primary_species_names(icomp), &
                                     string,MAXWORDLENGTH)) then
                option%io_buffer = &
                         'pH specified as constraint (constraint =' // &
                         trim(constraint_name) // &
                         ') for species other than H+ or OH-: ' // &
                         trim(reaction%primary_species_names(icomp))
                call printErrMsg(option)
              endif
            endif
            free_conc(icomp) = 10.d0**(-conc(icomp))
          else
            option%io_buffer = &
                     'pH specified as constraint (constraint =' // &
                     trim(constraint_name) // &
                     '), but H+ not found in chemical species.'
            call printErrMsg(option)
          endif
        endif        
      case(CONSTRAINT_MINERAL)
        if (.not.use_prev_soln_as_guess) then
          free_conc(icomp) = conc(icomp)*convert_molar_to_molal ! guess
        endif
      case(CONSTRAINT_GAS, CONSTRAINT_SUPERCRIT_CO2)
        if (conc(icomp) <= 0.d0) then ! log form
          conc(icomp) = 10.d0**conc(icomp) ! conc log10 partial pressure gas
        endif
        ! free_conc guess set above
    end select
  enddo

  rt_auxvar%pri_molal = free_conc

  num_iterations = 0
  num_it_act_coef_turned_on = 0

  ! if previous solution is provided as a guess, it should be close enough
  ! to use activity coefficients right away. - geh
  ! essentially the same as:
  !   compute_activity_coefficients = (use_prev_soln_as_guess == PETSC_TRUE)
  compute_activity_coefs = use_prev_soln_as_guess
  
  do

    ! there are cases (e.g. in a linear update) where components constrained by
    ! a free ion concentration can change a bit.  This reset hopefullly 
    ! eliminates this issue.
    do icomp = 1, reaction%naqcomp
      select case(constraint_type(icomp))
        case(CONSTRAINT_FREE,CONSTRAINT_LOG)
          rt_auxvar%pri_molal(icomp) = free_conc(icomp)
      end select
    enddo          
  
    if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF .and. &
        compute_activity_coefs) then
      call RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
    endif
    call RTotal(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    
    ! geh - for debugging
    !call RTPrintAuxVar(rt_auxvar,reaction,option)

    Jac = 0.d0

! for colloids later on    
!    if (reaction%ncoll > 0) then
!      do idof = istartcoll, iendcoll
!        Jac(idof,idof) = 1.d0
!      enddo
!    endif
        
    do icomp = 1, reaction%naqcomp

      select case(constraint_type(icomp))
      
        case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
        
          ! units = mol/L water
          Res(icomp) = rt_auxvar%total(icomp,1) - total_conc(icomp)
          ! dtotal units = kg water/L water

          ! Jac units = kg water/L water
          Jac(icomp,:) = rt_auxvar%aqueous%dtotal(icomp,:,1)
          
        case(CONSTRAINT_TOTAL_SORB)
        
          ! units = mol/m^3 bulk
          Res(icomp) = rt_auxvar%total_sorb_eq(icomp) - total_conc(icomp)
          ! dtotal_sorb units = kg water/m^3 bulk
          ! Jac units = kg water/m^3 bulk
          Jac(icomp,:) = rt_auxvar%dtotal_sorb_eq(icomp,:)

        case(CONSTRAINT_FREE,CONSTRAINT_LOG)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
!          Jac(:,icomp) = 0.d0
          Jac(icomp,icomp) = 1.d0
          
        case(CONSTRAINT_CHARGE_BAL)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
          do jcomp = 1, reaction%naqcomp
            Res(icomp) = Res(icomp) + reaction%primary_spec_Z(jcomp) * &
              rt_auxvar%total(jcomp,1)
            do kcomp = 1, reaction%naqcomp
              Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                reaction%primary_spec_Z(kcomp)*rt_auxvar%aqueous%dtotal(kcomp,jcomp,1)
            enddo
          enddo
          if (rt_auxvar%pri_molal(icomp) < 1.d-20 .and. &
              .not.charge_balance_warning_flag) then
            if ((Res(icomp) > 0.d0 .and. &
                 reaction%primary_spec_Z(icomp) > 0.d0) .or. &
                (Res(icomp) < 0.d0 .and. &
                 reaction%primary_spec_Z(icomp) < 0.d0)) then
              option%io_buffer = &
                       'Charge balance species ' // &
                       trim(reaction%primary_species_names(icomp)) // &
                       ' may not satisfy constraint ' // &
                       trim(constraint_name) // &
                       '.  Molality already below 1.e-20.'
              call printMsg(option)
              charge_balance_warning_flag = PETSC_TRUE
              rt_auxvar%pri_molal(icomp) = 1.e-3 ! reset guess
            endif
          endif
          
        case(CONSTRAINT_PH)
        
          Res(icomp) = 0.d0
          Jac(icomp,:) = 0.d0
          if (associated(reaction%species_idx)) then
            if (reaction%species_idx%h_ion_id > 0) then ! conc(icomp) = 10**-pH
              rt_auxvar%pri_molal(icomp) = 10.d0**(-conc(icomp)) / &
                                            rt_auxvar%pri_act_coef(icomp)
              Jac(icomp,icomp) = 1.d0
            else ! H+ is a complex
            
              icplx = abs(reaction%species_idx%h_ion_id)
              
              ! compute secondary species concentration
              ! *note that the sign was flipped below
              lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

              ! activity of water
              if (reaction%eqcplxh2oid(icplx) > 0) then
                lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
              endif

              do jcomp = 1, reaction%eqcplxspecid(0,icplx)
                comp_id = reaction%eqcplxspecid(jcomp,icplx)
                lnQK = lnQK + reaction%eqcplxstoich(jcomp,icplx)* &
                              log(rt_auxvar%pri_molal(comp_id)* &
                              rt_auxvar%pri_act_coef(comp_id))
              enddo
              lnQK = lnQK + conc(icomp)*LOG_TO_LN ! this is log activity H+
              QK = exp(lnQK)
              
              Res(icomp) = 1.d0 - QK

              do jcomp = 1,reaction%eqcplxspecid(0,icplx)
                comp_id = reaction%eqcplxspecid(jcomp,icplx)
                Jac(icomp,comp_id) = -QK/rt_auxvar%pri_molal(comp_id)* &
                                          reaction%eqcplxstoich(jcomp,icplx)
              enddo
            endif
          endif
                      
        case(CONSTRAINT_MINERAL)

          imnrl = constraint_id(icomp)
          lnQK = -mineral_reaction%mnrl_logK(imnrl)*LOG_TO_LN

          ! activity of water
          if (mineral_reaction%mnrlh2oid(imnrl) > 0) then
            lnQK = lnQK + mineral_reaction%mnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
          endif

          ! compute ion activity product
          do jcomp = 1, mineral_reaction%mnrlspecid(0,imnrl)
            comp_id = mineral_reaction%mnrlspecid(jcomp,imnrl)
            lnQK = lnQK + mineral_reaction%mnrlstoich(jcomp,imnrl)* &
                          log(rt_auxvar%pri_molal(comp_id)* &
                          rt_auxvar%pri_act_coef(comp_id))
          enddo
          
          Res(icomp) = lnQK

          do jcomp = 1,mineral_reaction%mnrlspecid(0,imnrl)
            comp_id = mineral_reaction%mnrlspecid(jcomp,imnrl)
            Jac(icomp,comp_id) = mineral_reaction%mnrlstoich(jcomp,imnrl)/ &
              rt_auxvar%pri_molal(comp_id)
          enddo
  
        case(CONSTRAINT_GAS)

          igas = constraint_id(icomp)
          lnQK = -reaction%gas%paseqlogK(igas)*LOG_TO_LN
 
          ! divide K by RT
          !lnQK = lnQK - log((auxvar%temp+273.15d0)*IDEAL_GAS_CONSTANT)
          
          ! activity of water
          if (reaction%gas%paseqh2oid(igas) > 0) then
            lnQK = lnQK + reaction%gas%paseqh2ostoich(igas)*rt_auxvar%ln_act_h2o
          endif

          ! compute ion activity product
          do jcomp = 1, reaction%gas%paseqspecid(0,igas)
            comp_id = reaction%gas%paseqspecid(jcomp,igas)
            lnQK = lnQK + reaction%gas%paseqstoich(jcomp,igas)* &
              log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
          enddo
          
!         QK = exp(lnQK)
          
!         Res(icomp) = QK - conc(icomp)
          Res(icomp) = lnQK - log(conc(icomp)) ! gas pressure
          Jac(icomp,:) = 0.d0
          do jcomp = 1,reaction%gas%paseqspecid(0,igas)
            comp_id = reaction%gas%paseqspecid(jcomp,igas)
!           Jac(icomp,comp_id) = QK/auxvar%primary_spec(comp_id)* &
!                                reaction%gas%paseqstoich(jcomp,igas)
            Jac(icomp,comp_id) = reaction%gas%paseqstoich(jcomp,igas)/ &
              rt_auxvar%pri_molal(comp_id)
          enddo

      end select
    enddo
    
    maximum_residual = maxval(abs(Res))
    
    if (reaction%use_log_formulation) then
      ! force at least 4 log updates, then cycle the next 5 updates between
      ! log/linear.  This improves convergence of linear problems or 
      ! primary components with no complexes, reactions, etc. (e.g. tracers)
      if (num_iterations > 3 .and. num_iterations < 9) then
        use_log_formulation = (mod(num_iterations,2) == 0)
      else
        use_log_formulation = PETSC_TRUE
      endif
    else
      use_log_formulation = PETSC_FALSE
    endif
      
    call RSolve(Res,Jac,rt_auxvar%pri_molal,update,reaction%naqcomp, &
                use_log_formulation)
    
    prev_molal = rt_auxvar%pri_molal

    if (use_log_formulation) then
      update = dsign(1.d0,update)*min(dabs(update),reaction%max_dlnC)
      rt_auxvar%pri_molal = rt_auxvar%pri_molal*exp(-update)    
    else ! linear update
      ! ensure non-negative concentration
      min_ratio = 1.d20 ! large number
      do icomp = 1, reaction%naqcomp
        if (prev_molal(icomp) <= update(icomp)) then
          ratio = abs(prev_molal(icomp)/update(icomp))
          if (ratio < min_ratio) min_ratio = ratio
        endif
      enddo
      if (min_ratio <= 1.d0) then
        ! scale by 0.99 to make the update slightly smaller than the min_ratio
        update = update*min_ratio*0.99d0
      endif
      rt_auxvar%pri_molal = prev_molal - update
      ! could use:
      ! rt_auxvar%pri_molal = prev_molal - update * minval(abs(prev_molal/update))
    endif
    
    ! check to ensure that minimum concentration is not less than or equal
    ! to zero
    tempreal = minval(rt_auxvar%pri_molal)
    if (tempreal <= 0.d0) then
      option%io_buffer = 'ERROR: Zero concentrations found in ' // &
        'constraint "' // trim(constraint_name) // '".'
      call printMsgByRank(option)
      ! now figure out which species have zero concentrations
      do idof = 1, reaction%naqcomp
        if (rt_auxvar%pri_molal(idof) <= 0.d0) then
          write(string,*) rt_auxvar%pri_molal(idof)
          option%io_buffer = '  Species "' // &
            trim(reaction%primary_species_names(idof)) // &
            '" has zero concentration (' // &
            trim(adjustl(string)) // ').'
          call printMsgByRank(option)
        endif
      enddo
      option%io_buffer = 'Free ion concentations RESULTING from ' // &
        'constraint concentrations must be positive.'
      call printErrMsgByRank(option)
    endif
    
#if 0
!geh cannot use this check as for many problems (e.g. Hanford 300 Area U), the 
!    concentrations temporarily go well above 100.

    ! check for excessively large maximum values, which likely indicates
    ! reaction going awry.
    tempreal = maxval(rt_auxvar%pri_molal)
    ! allow a few iterations; sometime charge balance constraint jumps
    ! during initial iterations
    if (tempreal > 100.d0 .and. num_iterations > 500) then
      !geh: for some reason, needs the array rank included in call to maxloc
      idof = maxloc(rt_auxvar%pri_molal,1)
      option%io_buffer = 'ERROR: Excessively large concentration for ' // &
        'species "' // trim(reaction%primary_species_names(idof)) // &
        '" in constraint "' // trim(constraint_name) // &
        '" in ReactionEquilibrateConstraint. Email input deck to ' // &
        'pflotran-dev@googlegroups.com.'
      call printErrMsg(option)
    endif
#endif
    
    maximum_relative_change = maxval(abs((rt_auxvar%pri_molal-prev_molal)/ &
                                         prev_molal))
    
    num_iterations = num_iterations + 1
    
    if (mod(num_iterations,1000) == 0) then
100   format('Constraint iteration count has exceeded: ',i5)
      write(option%io_buffer,100) num_iterations
      call printMsg(option)
      option%io_buffer = 'species_name  prev_free_ion  residual'
      call printMsg(option)
      do icomp=1,reaction%naqcomp
        write(option%io_buffer,200) reaction%primary_species_names(icomp), &
        prev_molal(icomp),Res(icomp)
        call printMsg(option)
      enddo
200   format(a12,1x,1p2e12.4)
      if (num_iterations >= 10000) then
        print *, 'cell id (natural):', option%iflag 
        print *, 'constraint:', conc
        print *, 'constraint type:', constraint_type
        print *, 'free_conc:', free_conc
        option%io_buffer = 'Equilibration of constraint "' // &
          trim(constraint_name) // &
          '" stopping due to excessive iteration count!'
        call printErrMsgByRank(option)
      endif
    endif
    
    ! check for convergence
    if (maximum_residual < reaction%max_residual_tolerance .and. &
        maximum_relative_change < reaction%max_relative_change_tolerance) then
      ! Need some sort of convergence before we kick in activities
      if (compute_activity_coefs .and. &
          ! With some constraints (e.g. pH), the total component concentration
          ! is not updated immediately after activity coefficients are turned
          ! on.  Therefore, we need at least two iterations to declare 
          ! convergence. - geh
          num_iterations - num_it_act_coef_turned_on > 1) exit
      if (.not. compute_activity_coefs) &
        num_it_act_coef_turned_on = num_iterations
      compute_activity_coefs = PETSC_TRUE
    endif

  enddo

  ! once equilibrated, compute sorbed concentrations
  if (reaction%nsorb > 0) then
    if (reaction%neqsorb > 0) then
      call RTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    endif
  endif

  ! do not scale by molal_to_molar since it could be 1.d0 if MOLAL flag set
  aq_species_constraint%basis_molarity = rt_auxvar%pri_molal* &
                                 global_auxvar%den_kg(option%liquid_phase)/ &
                                 1000.d0

#if 0
  call RCalculateCompression(global_auxvar,rt_auxvar,reaction,option)
#endif

! this is performed above
!  if (associated(colloid_constraint%colloids)) then                        
!    colloid_constraint%colloids%basis_conc_mob = rt_auxvar%colloid%conc_mob* &
!                           global_auxvar%den_kg(option%liquid_phase)/1000.d0
!    colloid_constraint%colloids%basis_conc_imb = rt_auxvar%colloid%conc_imb* &
!                           global_auxvar%den_kg(option%liquid_phase)/1000.d0
!  endif
    
!  write(option%io_buffer,111) trim(constraint_name),num_iterations
!  call printMsg(option)
!111 format(' Equilibrate Constraint: ',a30,i4)

end subroutine ReactionEquilibrateConstraint

! ************************************************************************** !

subroutine ReactionPrintConstraint(constraint_coupler,reaction,option)
  ! 
  ! Prints a constraint associated with reactive
  ! transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/28/08
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  use Transport_Constraint_module

  implicit none
  
  type(option_type) :: option
  type(tran_constraint_coupler_type) :: constraint_coupler
  type(reaction_type), pointer :: reaction
  
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint

  type(mineral_type), pointer :: mineral_reaction
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, icomp, irxn, j, jj, ncomp, ncplx, ieqrxn
  PetscInt :: icplx, icplx2
  PetscInt :: imnrl,igas
  PetscInt :: eqcplxsort(reaction%neqcplx+1)
  PetscInt :: eqcplxid(reaction%neqcplx+1)
  PetscInt :: eqminsort(reaction%mineral%nmnrl)
  PetscInt, allocatable :: eqsrfcplxsort(:)
  PetscBool :: finished, found
  PetscReal :: conc, conc2
  PetscReal :: lnQK(reaction%mineral%nmnrl), QK(reaction%mineral%nmnrl)
  PetscReal :: lnQKgas(reaction%gas%npassive_gas)
  PetscReal :: QKgas(reaction%gas%npassive_gas)
  PetscReal :: charge_balance, ionic_strength
  PetscReal :: percent(reaction%neqcplx+1)
  PetscReal :: totj, retardation, kd, ph
  PetscInt :: comp_id, jcomp
  PetscInt :: icount
  PetscInt :: iphase, ifo2
  PetscReal :: bulk_vol_to_fluid_vol, molar_to_molal, molal_to_molar
  PetscReal :: sum_molality, sum_mass, mole_fraction_h2o, mass_fraction_h2o, &
               mass_fraction_co2, mole_fraction_co2
  PetscReal :: ehfac,eh,pe,tk
  PetscReal :: affinity

  aq_species_constraint => constraint_coupler%aqueous_species
  mineral_constraint => constraint_coupler%minerals
  
  iphase = 1

  90 format(2x,76('-'))

  write(option%fid_out,'(/,''  Constraint: '',a)') &
    trim(constraint_coupler%constraint_name)

  rt_auxvar => constraint_coupler%rt_auxvar
  global_auxvar => constraint_coupler%global_auxvar

  mineral_reaction => reaction%mineral

  select case(option%iflowmode)
    case(NULL_MODE)
      global_auxvar%den_kg(iphase) = &
        option%reference_density(option%liquid_phase)
      global_auxvar%temp = option%reference_temperature
      global_auxvar%sat(iphase) = option%reference_saturation
  end select
        
  bulk_vol_to_fluid_vol = option%reference_porosity* &
                          global_auxvar%sat(iphase)*1000.d0

! compute mole and mass fractions of H2O
  if (reaction%use_full_geochemistry) then
    sum_molality = 0.d0
    do icomp = 1, reaction%naqcomp
      if (icomp /= reaction%species_idx%h2o_aq_id) then
        sum_molality = sum_molality + rt_auxvar%pri_molal(icomp)
      endif
    enddo
    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        sum_molality = sum_molality + rt_auxvar%sec_molal(i)
      enddo
    endif
    mole_fraction_h2o = 1.d0/(1.d0+FMWH2O*sum_molality*1.d-3)

    sum_mass = 0.d0
    do icomp = 1, reaction%naqcomp
      if (icomp /= reaction%species_idx%h2o_aq_id) then
        sum_mass = sum_mass + &
          reaction%primary_spec_molar_wt(icomp)*rt_auxvar%pri_molal(icomp)
      endif
    enddo
    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        sum_mass = sum_mass + reaction%eqcplx_molar_wt(i)*rt_auxvar%sec_molal(i)
      enddo
    endif
    mass_fraction_h2o = 1.d0/(1.d0 + sum_mass*1.d-3)
  endif
  
  molal_to_molar = global_auxvar%den_kg(iphase)/1000.d0
  molar_to_molal = 1.d0/molal_to_molar
    
  if (.not.reaction%use_full_geochemistry) then
    100 format(/,'  species       molality')  
    write(option%fid_out,100)
    101 format(2x,a12,es12.4)
    do icomp = 1, reaction%naqcomp
      write(option%fid_out,101) reaction%primary_species_names(icomp), &
                                rt_auxvar%pri_molal(icomp)
    enddo
  else

    if (.not.option%use_isothermal) then
      call RUpdateTempDependentCoefs(global_auxvar,reaction,PETSC_TRUE,option)
    endif
  

201 format(a20,i5)
203 format(a20,f8.4)
204 format(a20,es12.4)

    write(option%fid_out,90)
    write(option%fid_out,201) '      iterations: ', &
      constraint_coupler%num_iterations

    if (associated(reaction%species_idx)) then
!    output pH
      if (reaction%species_idx%h_ion_id > 0) then
        ph = &
          -log10(rt_auxvar%pri_molal(reaction%species_idx%h_ion_id)* &
                 rt_auxvar%pri_act_coef(reaction%species_idx%h_ion_id))
      else if (reaction%species_idx%h_ion_id < 0) then
        ph = &
          -log10(rt_auxvar%sec_molal(abs(reaction%species_idx%h_ion_id))* &
                 rt_auxvar%sec_act_coef(abs(reaction%species_idx%h_ion_id)))
      endif
      if (reaction%species_idx%h_ion_id > 0 .or. reaction%species_idx%h_ion_id < 0) &
      write(option%fid_out,203) '              pH: ',ph

!    output Eh and pe
      if (reaction%species_idx%o2_gas_id > 0 .and. (reaction%species_idx%h_ion_id > 0 &
          .or. reaction%species_idx%h_ion_id < 0)) then

        ifo2 = reaction%species_idx%o2_gas_id
      
      ! compute gas partial pressure
        lnQKgas(ifo2) = -reaction%gas%paseqlogK(ifo2)*LOG_TO_LN
      
      ! activity of water
        if (reaction%gas%paseqh2oid(ifo2) > 0) then
          lnQKgas(ifo2) = lnQKgas(ifo2) + reaction%gas%paseqh2ostoich(ifo2)*rt_auxvar%ln_act_h2o
        endif
        do jcomp = 1, reaction%gas%paseqspecid(0,ifo2)
          comp_id = reaction%gas%paseqspecid(jcomp,ifo2)
          lnQKgas(ifo2) = lnQKgas(ifo2) + reaction%gas%paseqstoich(jcomp,ifo2)* &
                      log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
        enddo

        tk = global_auxvar%temp+273.15d0
        ehfac = IDEAL_GAS_CONSTANT*tk*LOG_TO_LN/FARADAY
        eh = ehfac*(-4.d0*ph+lnQKgas(ifo2)*LN_TO_LOG+logKeh(tk))/4.d0
        pe = eh/ehfac

!       pe = (-4.d0*ph+lnQKgas(ifo2)*LN_TO_LOG+logKeh(tk))/4.d0
!       eh = pe*ehfac
        write(option%fid_out,203) '              pe: ',pe
        write(option%fid_out,203) '              Eh: ',eh
      endif
    endif
    
    ionic_strength = 0.d0
    charge_balance = 0.d0
    do icomp = 1, reaction%naqcomp      
      charge_balance = charge_balance + rt_auxvar%total(icomp,1)* &
                                        reaction%primary_spec_Z(icomp)
      ionic_strength = ionic_strength + rt_auxvar%pri_molal(icomp)* &
        reaction%primary_spec_Z(icomp)*reaction%primary_spec_Z(icomp)
    enddo
    
    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        ionic_strength = ionic_strength + rt_auxvar%sec_molal(i)* &
                                          reaction%eqcplx_Z(i)* &
                                          reaction%eqcplx_Z(i)
      enddo
    endif
    ionic_strength = 0.5d0 * ionic_strength
    
    write(option%fid_out,'(a20,es12.4,a8)') '  ionic strength: ', &
      ionic_strength,' [mol/L]'
    write(option%fid_out,204) '  charge balance: ', charge_balance
    
    write(option%fid_out,'(a20,1pe12.4,a5)') '        pressure: ', &
      global_auxvar%pres(1),' [Pa]'
    write(option%fid_out,'(a20,f8.2,a4)') '     temperature: ', &
      global_auxvar%temp,' [C]'
    write(option%fid_out,'(a20,f8.2,a9)') '     density H2O: ', &
      global_auxvar%den_kg(1),' [kg/m^3]'
    write(option%fid_out,'(a20,1p2e12.4,a9)') 'ln / activity H2O: ', &
      rt_auxvar%ln_act_h2o,exp(rt_auxvar%ln_act_h2o),' [---]'
    write(option%fid_out,'(a20,1pe12.4,a9)') 'mole fraction H2O: ', &
      mole_fraction_h2o,' [---]'
    write(option%fid_out,'(a20,1pe12.4,a9)') 'mass fraction H2O: ', &
      mass_fraction_h2o,' [---]'

    write(option%fid_out,90)

    102 format(/,'  primary               free        total')  
    103 format('  species               molal       molal       act coef     constraint')  
    write(option%fid_out,102)
    write(option%fid_out,103)
    write(option%fid_out,90)
  
    104 format(2x,a20,es12.4,es12.4,es12.4,4x,a)
    do icomp = 1, reaction%naqcomp
      select case(aq_species_constraint%constraint_type(icomp))
        case(CONSTRAINT_NULL,CONSTRAINT_TOTAL)
          string = 'total aq'
        case(CONSTRAINT_TOTAL_SORB)
          string = 'total sorb'
        case(CONSTRAINT_FREE)
          string = 'free'
        case(CONSTRAINT_CHARGE_BAL)
          string = 'chrg'
        case(CONSTRAINT_LOG)
          string = 'log'
        case(CONSTRAINT_PH)
          string = 'pH'
        case(CONSTRAINT_MINERAL,CONSTRAINT_GAS)
          string = aq_species_constraint%constraint_aux_string(icomp)
        case(CONSTRAINT_SUPERCRIT_CO2)
          string = 'SC ' // aq_species_constraint%constraint_aux_string(icomp)
      end select
      write(option%fid_out,104) reaction%primary_species_names(icomp), &
                                rt_auxvar%pri_molal(icomp), &
                                rt_auxvar%total(icomp,1)*molar_to_molal, &
                                rt_auxvar%pri_act_coef(icomp), &
                                trim(string)
    enddo 
  endif 
      
  if (reaction%neqcplx > 0) then    
    ! sort complex concentrations from largest to smallest
    do i = 1, reaction%neqcplx
      eqcplxsort(i) = i
    enddo
    do
      finished = PETSC_TRUE
      do i = 1, reaction%neqcplx-1
        icplx = eqcplxsort(i)
        icplx2 = eqcplxsort(i+1)
        if (rt_auxvar%sec_molal(icplx) < &
            rt_auxvar%sec_molal(icplx2)) then
          eqcplxsort(i) = icplx2
          eqcplxsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo
            
    110 format(/,'  complex               molality    act coef        logK')  
    write(option%fid_out,110)
    write(option%fid_out,90)
    111 format(2x,a20,es12.4,es12.4,2x,es12.4)
    do i = 1, reaction%neqcplx ! for each secondary species
      icplx = eqcplxsort(i)
      write(option%fid_out,111) reaction%secondary_species_names(icplx), &
                                rt_auxvar%sec_molal(icplx), &
                                rt_auxvar%sec_act_coef(icplx), &
                                reaction%eqcplx_logK(icplx)
    enddo 

    !print speciation precentages
    write(option%fid_out,92)
    92 format(/)
    134 format(2x,'complex species       percent   molality')
    135 format(2x,'primary species: ',a20,2x,' total conc: ',1pe12.4)
    136 format(2x,a20,2x,f6.2,2x,1pe12.4,1p2e12.4)
    do icomp = 1, reaction%naqcomp
    
      eqcplxsort = 0
      eqcplxid = 0
      percent = 0.d0
      totj = 0.d0
      
      icount = 0
      do icplx = 1, reaction%neqcplx
        found = PETSC_FALSE
        do i = 1, reaction%eqcplxspecid(0,icplx)
          if (reaction%eqcplxspecid(i,icplx) == icomp) then
            icount = icount + 1
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (found) then
          eqcplxid(icount) = icplx
          percent(icount) = dabs(rt_auxvar%sec_molal(icplx)* &
                                 reaction%eqcplxstoich(i,icplx))
          totj = totj + percent(icount)
        endif
      enddo
      icount = icount + 1
      eqcplxid(icount) = -icomp
      percent(icount) = rt_auxvar%pri_molal(icomp)
      totj = totj + percent(icount)
      percent = percent / totj
      
      eqcplxsort = 0
      do i = 1, icount
        eqcplxsort(i) = i
      enddo
      
      do
        finished = PETSC_TRUE
        do i = 1, icount-1
          icplx = eqcplxsort(i)
          icplx2 = eqcplxsort(i+1)
          if (percent(abs(icplx)) < percent(abs(icplx2))) then
            eqcplxsort(i) = icplx2
            eqcplxsort(i+1) = icplx
            finished = PETSC_FALSE
          endif
        enddo
        if (finished) exit
      enddo

      write(option%fid_out,90)
      write(option%fid_out,135) reaction%primary_species_names(icomp), &
                                rt_auxvar%total(icomp,iphase)
      write(option%fid_out,134)
      write(option%fid_out,90)
      do i = 1, icount
        j = eqcplxsort(i)
        if (percent(j) < 0.0001d0) cycle
        icplx = eqcplxid(j)
        if (icplx < 0) then
          icplx = abs(icplx)
          write(option%fid_out,136) reaction%primary_species_names(icplx), &
                                    percent(j)*100.d0, &
                                    rt_auxvar%pri_molal(icplx)
        else
          write(option%fid_out,136) reaction%secondary_species_names(icplx), &
                                    percent(j)*100.d0, &
                                    rt_auxvar%sec_molal(icplx)
        endif
      enddo
    enddo

  endif 
          

  ! Ion Exchange
  if (reaction%neqionxrxn > 0) then
    write(option%fid_out,125)
    do irxn = 1, reaction%neqionxrxn
      write(option%fid_out,90)
      write(option%fid_out,126) reaction%eqionx_rxn_CEC(irxn)
      write(option%fid_out,127)
      write(option%fid_out,90)
      ncomp = reaction%eqionx_rxn_cationid(0,irxn)
      do i = 1, ncomp
        icomp = reaction%eqionx_rxn_cationid(i,irxn)
        kd = rt_auxvar%eqionx_conc(i,irxn)/rt_auxvar%total(icomp,iphase) & 
                      /bulk_vol_to_fluid_vol
        write(option%fid_out,128) reaction%primary_species_names(icomp), &
          reaction%eqionx_rxn_k(i,irxn), & 
          rt_auxvar%eqionx_conc(i,irxn), &
          kd
      enddo
    enddo
    125 format(/,2x,'ion-exchange reactions')
    126 format(2x,'CEC = ',1pe12.4)
    127 format(2x,'cation    selectivity coef.   sorbed conc.     Kd',&
               /,30x,'[mol/m^3]')
    128 format(2x,a8,2x,1pe12.4,4x,1pe12.4,4x,1pe12.4,4x,1pe12.4)
  endif
  
! total retardation from ion exchange and equilibrium surface complexation
  if (reaction%neqsorb > 0) then
    write(option%fid_out,1128)
    write(option%fid_out,90)
    do jcomp = 1, reaction%naqcomp
      if (abs(rt_auxvar%total(jcomp,iphase)) > 0.d0) &
      retardation = 1.d0 + rt_auxvar%total_sorb_eq(jcomp)/bulk_vol_to_fluid_vol &
        /rt_auxvar%total(jcomp,iphase)
      totj = rt_auxvar%total(jcomp,iphase)+rt_auxvar%total_sorb_eq(jcomp)/bulk_vol_to_fluid_vol
      write(option%fid_out,129) reaction%primary_species_names(jcomp), &
        totj,retardation
    enddo
   1128 format(/,2x,'primary species     total(aq+sorbed)    total retardation', &
               /,25x,'[mol/L]',15x,'1+Kd')
    129 format(2x,a12,8x,1pe12.4,8x,1pe12.4)
  endif
  
  if (mineral_reaction%nmnrl > 0) then
  
    130 format(/,'  mineral                             log SI       Affinity     log K', &
           /,51x,'[kJ/mol]')
    131 format(2x,a30,2x,2f12.4,2x,1pe12.4)

    do imnrl = 1, mineral_reaction%nmnrl
      ! compute saturation
      lnQK(imnrl) = -mineral_reaction%mnrl_logK(imnrl)*LOG_TO_LN
      if (mineral_reaction%mnrlh2oid(imnrl) > 0) then
        lnQK(imnrl) = lnQK(imnrl) + mineral_reaction%mnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
      endif
      do jcomp = 1, mineral_reaction%mnrlspecid(0,imnrl)
        comp_id = mineral_reaction%mnrlspecid(jcomp,imnrl)
        lnQK(imnrl) = lnQK(imnrl) + mineral_reaction%mnrlstoich(jcomp,imnrl)* &
                      log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
      enddo
      QK(imnrl) = exp(lnQK(imnrl))    
    enddo

    ! sort mineral saturation indices from largest to smallest
    do i = 1, mineral_reaction%nmnrl
      eqminsort(i) = i
    enddo
    do
      finished = PETSC_TRUE
      do i = 1, mineral_reaction%nmnrl-1
        icplx = eqminsort(i)
        icplx2 = eqminsort(i+1)
        if (QK(icplx) < QK(icplx2)) then
          eqminsort(i) = icplx2
          eqminsort(i+1) = icplx
          finished = PETSC_FALSE
        endif
      enddo
      if (finished) exit
    enddo

    write(option%fid_out,130)
    write(option%fid_out,90)
  
    do imnrl = 1, mineral_reaction%nmnrl
      i = eqminsort(imnrl)
      affinity = -1.d0*IDEAL_GAS_CONSTANT*(global_auxvar%temp+273.15d0)*lnQK(i)
      write(option%fid_out,131) mineral_reaction%mineral_names(i), &
                                lnQK(i)*LN_TO_LOG, affinity, &
                                mineral_reaction%mnrl_logK(i)
    enddo
  endif
    
  if (reaction%gas%npassive_gas > 0) then
    
    132 format(/,'  gas        log part. press.  part. press. [bars]      log K')
    133 format(2x,a10,2x,1pe12.4,6x,1pe12.4,8x,1pe12.4)
    
    write(option%fid_out,132)
    write(option%fid_out,90)
    
    do igas = 1, reaction%gas%npassive_gas
      
      ! compute gas partial pressure
      lnQKgas(igas) = -reaction%gas%paseqlogK(igas)*LOG_TO_LN
      
      ! divide K by RT
      !lnQKgas = lnQKgas - log((auxvar%temp+273.15d0)*IDEAL_GAS_CONSTANT)
      
      ! activity of water
      if (reaction%gas%paseqh2oid(igas) > 0) then
        lnQKgas(igas) = lnQKgas(igas) + reaction%gas%paseqh2ostoich(igas)* &
                                        rt_auxvar%ln_act_h2o
      endif

      do jcomp = 1, reaction%gas%paseqspecid(0,igas)
        comp_id = reaction%gas%paseqspecid(jcomp,igas)
        lnQKgas(igas) = lnQKgas(igas) + reaction%gas%paseqstoich(jcomp,igas)* &
                      log(rt_auxvar%pri_molal(comp_id)*rt_auxvar%pri_act_coef(comp_id))
      enddo
      
      QKgas(igas) = exp(lnQKgas(igas))
          
      write(option%fid_out,133) reaction%gas%passive_names(igas),lnQKgas(igas)*LN_TO_LOG, &
      QKgas(igas),reaction%gas%paseqlogK(igas)
    enddo
  endif

#ifdef AMANZI_BGD
  ! output constraints for amanzi cfg formatted input
  if (OptionPrintToFile(option)) then
    string = trim(option%global_prefix) // '-' // &
             trim(constraint_coupler%constraint_name) // '.txt'
    open(unit=86,file=trim(string))

    write(86,'("# pflotran constraint preprocessing :")')
    !call date_and_time(date=word,time=word2)
    ! prints garbage? need to clear memory?
    !write(86,'("#        date : ",a,"   ",a)') trim(word), trim(word2)
    write(86,'("#       input : ",a)') trim(option%input_filename)
    write(86,'(/,"# Constraint: ",a)') &
         trim(constraint_coupler%constraint_name)

    write(86,'(/,"[total]")')
    do icomp = 1, reaction%naqcomp
      write(86,'(a," = ",1es13.6)') trim(reaction%primary_species_names(icomp)), &
                                rt_auxvar%pri_molal(icomp)
    enddo

    write(86,'(/,"[free_ion]")')
    do icomp = 1, reaction%naqcomp
      write(86,'(a," = ",1es13.6)') trim(reaction%primary_species_names(icomp)), &
                                rt_auxvar%total(icomp,1)*molar_to_molal
    enddo

    write(86,'(/,"[minerals]")')
    do imnrl = 1, mineral_reaction%nkinmnrl
      write(86,'(a," = ",f6.3)') trim(mineral_reaction%kinmnrl_names(imnrl)), &
                                mineral_reaction%mnrl_volfrac(imnrl)
    enddo

    if (associated(rt_auxvar%total_sorb_eq)) then
      write(86,'(/,"[total_sorbed]")')
      do icomp = 1, reaction%naqcomp
        write(86,'(a," = ",1es13.6)') trim(reaction%primary_species_names(icomp)), &
                                  rt_auxvar%total_sorb_eq(icomp)
      enddo
    endif

    write(86,'(/,"[ion_exchange]")')
    do icomp = 1, reaction%neqionxrxn
      write(86, '("X- = ",1es13.6)') reaction%eqionx_rxn_CEC(icomp)
    enddo
    close(86)
  endif
#endif
! end AMANZI_BGD

end subroutine ReactionPrintConstraint

! ************************************************************************** !

subroutine ReactionDoubleLayer(constraint_coupler,reaction,option)
  ! 
  ! Calculates double layer potential, surface charge, and
  ! sorbed surface complex concentrations
  ! 
  ! Author: Peter C. Lichtner
  ! Date: ???
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  use Transport_Constraint_module

  implicit none
  
  type(option_type) :: option
  type(tran_constraint_coupler_type) :: constraint_coupler
  type(reaction_type), pointer :: reaction

  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint

  type(mineral_type), pointer :: mineral_reaction

  PetscReal, parameter :: tk = 273.15d0
  PetscReal, parameter :: epsilon = 78.5d0
  PetscReal, parameter :: epsilon0 = 8.854187817d-12
  
  PetscReal :: fac, boltzmann, dbl_charge, surface_charge, ionic_strength, &
               charge_balance, potential, tempk, debye_length, &
               srfchrg_capacitance_model, capacitance
               
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)


  PetscReal :: free_site_conc
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total

  PetscInt :: iphase
  PetscInt :: i, j, icomp, icplx, irxn, ncomp, ncplx

  PetscReal :: site_density(2)
  PetscReal :: mobile_fraction
  PetscInt :: num_types_of_sites
  PetscInt :: isite

  PetscBool :: one_more


    rt_auxvar => constraint_coupler%rt_auxvar
    global_auxvar => constraint_coupler%global_auxvar

    iphase = 1
    global_auxvar%temp = option%reference_temperature
    tempk = tk + global_auxvar%temp
    
    potential = 0.1d0 ! initial guess
    boltzmann = exp(-FARADAY*potential/(IDEAL_GAS_CONSTANT*tempk))
        
    fac = sqrt(epsilon*epsilon0*IDEAL_GAS_CONSTANT*tempk)

    ionic_strength = 0.d0
    charge_balance = 0.d0
    dbl_charge = 0.d0
    do icomp = 1, reaction%naqcomp      
      charge_balance = charge_balance + reaction%primary_spec_Z(icomp)* &
                       rt_auxvar%total(icomp,1)
                                        
      ionic_strength = ionic_strength + reaction%primary_spec_Z(icomp)**2* &
                       rt_auxvar%pri_molal(icomp)
      dbl_charge = dbl_charge + rt_auxvar%pri_molal(icomp)* &
                   (boltzmann**reaction%primary_spec_Z(icomp) - 1.d0)
    enddo

    if (reaction%neqcplx > 0) then    
      do i = 1, reaction%neqcplx
        ionic_strength = ionic_strength + reaction%eqcplx_Z(i)**2* &
                         rt_auxvar%sec_molal(i)
        dbl_charge = dbl_charge + rt_auxvar%sec_molal(i)* &
                     (boltzmann**reaction%eqcplx_Z(i) - 1.d0)
      enddo
    endif
    ionic_strength = 0.5d0*ionic_strength
    if (dbl_charge > 0.d0) then
      dbl_charge = fac*sqrt(2.d0*dbl_charge)
    else
      print *,'neg. dbl_charge: ',dbl_charge
      dbl_charge = fac*sqrt(2.d0*(-dbl_charge))
    endif

    srfchrg_capacitance_model = FARADAY* &
      sqrt(2.d0*epsilon*epsilon0*ionic_strength*1.d3/(IDEAL_GAS_CONSTANT*tempk))

    surface_charge = 0.d0

    debye_length = sqrt(fac/(2.d0*ionic_strength*1.d3))/FARADAY
    capacitance = sqrt(2.d0*epsilon*epsilon0*ionic_strength*1.d3/ &
                    (IDEAL_GAS_CONSTANT*tempk)) * FARADAY
    
    print *,'========================='
    print *,'dbl: debye_length = ',debye_length
    print *,'surface charge = ',dbl_charge,surface_charge, &
      srfchrg_capacitance_model
    print *,'ionic strength = ',ionic_strength
    print *,'chrg bal. = ',charge_balance,' Tk = ',tempk,' Boltz. = ',boltzmann
    print *,'srfcmplx: ',rt_auxvar%eqsrfcplx_conc
    print *,'capacitance: ',capacitance
    print *,'========================='

!   compute surface complex concentrations  
    ln_conc = log(rt_auxvar%pri_molal)
    ln_act = ln_conc+log(rt_auxvar%pri_act_coef)


end subroutine ReactionDoubleLayer


! ************************************************************************** !

subroutine ReactionReadOutput(reaction,input,option)
  ! 
  ! Reads species to be printed in output
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/24/09
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Input_Aux_module
  use String_module  
  use Option_module
  use Variables_module, only : PRIMARY_MOLALITY, PRIMARY_MOLARITY, &
                               TOTAL_MOLALITY, TOTAL_MOLARITY  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  PetscBool :: found
  PetscInt :: temp_int

  type(aq_species_type), pointer :: cur_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_rxn_type), pointer :: cur_mineral
  type(immobile_species_type), pointer :: cur_immobile
  
  nullify(cur_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,OUTPUT,SPECIES_NAME')
    
    word = name
    call StringToUpper(word)
    select case(word)
      case('OFF')
        reaction%print_all_species = PETSC_FALSE
        reaction%print_all_primary_species = PETSC_FALSE
        reaction%print_all_secondary_species = PETSC_FALSE
        reaction%gas%print_all = PETSC_FALSE
        reaction%mineral%print_all = PETSC_FALSE
        reaction%print_pH = PETSC_FALSE
        reaction%print_Eh = PETSC_FALSE
        reaction%print_pe = PETSC_FALSE
        reaction%print_O2 = PETSC_FALSE
        reaction%print_kd = PETSC_FALSE
        reaction%print_total_sorb = PETSC_FALSE
        reaction%print_total_sorb_mobile = PETSC_FALSE
        reaction%print_colloid = PETSC_FALSE
        reaction%print_act_coefs = PETSC_FALSE
        reaction%print_total_component = PETSC_FALSE
        reaction%print_free_ion = PETSC_FALSE
      case('ALL')
        reaction%print_all_species = PETSC_TRUE
        reaction%print_all_primary_species = PETSC_TRUE
 !       reaction%print_all_secondary_species = PETSC_TRUE
 !       reaction%gas%print_all = PETSC_TRUE
        reaction%mineral%print_all = PETSC_TRUE
        reaction%immobile%print_all = PETSC_TRUE
!        reaction%print_pH = PETSC_TRUE
      case('PRIMARY_SPECIES')
        reaction%print_all_primary_species = PETSC_TRUE
!        reaction%print_pH = PETSC_TRUE
      case('SECONDARY_SPECIES')
        reaction%print_all_secondary_species = PETSC_TRUE
      case('GASES')
        reaction%gas%print_all = PETSC_TRUE
      case('MINERALS')
        reaction%mineral%print_all = PETSC_TRUE
      case('MINERAL_SATURATION_INDEX')
        reaction%mineral%print_saturation_index = PETSC_TRUE
      case('IMMOBILE')
        reaction%immobile%print_all = PETSC_TRUE
      case('PH')
        reaction%print_pH = PETSC_TRUE
      case('EH')
        reaction%print_Eh = PETSC_TRUE
      case('PE')
        reaction%print_pe = PETSC_TRUE
      case('O2')
        reaction%print_O2 = PETSC_TRUE
      case('KD')
        reaction%print_kd = PETSC_TRUE
      case('COLLOIDS')
        reaction%print_colloid = PETSC_TRUE
      case('TOTAL')
        reaction%print_total_component = PETSC_TRUE
      case('TOTAL_SORBED')
        reaction%print_total_sorb = PETSC_TRUE
      case('TOTAL_BULK')
        reaction%print_total_bulk = PETSC_TRUE
      case('TOTAL_SORBED_MOBILE')
        reaction%print_total_sorb_mobile = PETSC_TRUE
      case('FREE_ION')
        reaction%print_free_ion = PETSC_TRUE
      case('ACTIVITY_COEFFICIENTS')
        reaction%print_act_coefs = PETSC_TRUE
      case('MOLARITY')
        reaction%print_free_conc_type = PRIMARY_MOLARITY
        reaction%print_tot_conc_type = TOTAL_MOLARITY
      case('MOLALITY')
        reaction%print_free_conc_type = PRIMARY_MOLALITY
        reaction%print_tot_conc_type = TOTAL_MOLALITY
      case('AGE')
        reaction%print_age = PETSC_TRUE
        reaction%use_full_geochemistry = PETSC_TRUE
      case('AUXILIARY')
        reaction%print_auxiliary = PETSC_TRUE
      case ('SITE_DENSITY')
        call InputReadWord(input,option,name,PETSC_TRUE)  
        call InputErrorMsg(input,option,'Site Name', &
                           'CHEMISTRY,OUTPUT,SITE DENSITY')
      case default        
        found = PETSC_FALSE
        ! primary aqueous species
        if (.not.found) then
          cur_aq_spec => reaction%primary_species_list
          do
            if (.not.associated(cur_aq_spec)) exit
            if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
              cur_aq_spec%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_aq_spec => cur_aq_spec%next
          enddo
        endif
        ! secondary aqueous complex
        if (.not.found) then
          cur_aq_spec => reaction%secondary_species_list
          do
            if (.not.associated(cur_aq_spec)) exit
            if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
              cur_aq_spec%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_aq_spec => cur_aq_spec%next
          enddo  
        endif
        ! gas
        if (.not.found) then
          ! be sure to check both lists, not skipping the passive list
          ! if found in the active list
          cur_gas_spec => reaction%gas%list
          do
            if (.not.associated(cur_gas_spec)) exit
            if (StringCompare(name,cur_gas_spec%name,MAXWORDLENGTH)) then
              cur_gas_spec%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_gas_spec => cur_gas_spec%next
          enddo  
        endif
        ! minerals
        if (.not.found) then
          cur_mineral => reaction%mineral%mineral_list
          do
            if (.not.associated(cur_mineral)) exit
            if (StringCompare(name,cur_mineral%name,MAXWORDLENGTH)) then
              cur_mineral%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_mineral => cur_mineral%next
          enddo
        endif
        ! immobile
        if (.not.found) then
          cur_immobile => reaction%immobile%list
          do  
            if (.not.associated(cur_immobile)) exit
            if (StringCompare(name,cur_immobile%name,MAXWORDLENGTH)) then
              cur_immobile%print_me = PETSC_TRUE
              found = PETSC_TRUE
              exit
            endif
            cur_immobile => cur_immobile%next
          enddo
        endif 
        if (.not.found) then
          option%io_buffer = 'CHEMISTRY,OUTPUT species name: '//trim(name)// &
                             ' not found among chemical species'
          call printErrMsg(option)
        endif
    end select

  enddo

  ! check to ensure that the user has listed FREE_ION or TOTAL is a primary
  ! species is listed for output
  found = PETSC_FALSE
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%print_me) then
      found = PETSC_TRUE
      exit
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo

  if ((found .or. reaction%print_all_primary_species .or. &
       reaction%print_all_species) .and. &
      .not.(reaction%print_total_component .or. &
            reaction%print_free_ion)) then
    option%io_buffer = 'FREE_ION or TOTAL must be specified to print a ' // &
      'primary species.'
    call printErrMsg(option)
  endif

end subroutine ReactionReadOutput

! ************************************************************************** !

subroutine RJumpStartKineticSorption(rt_auxvar,global_auxvar, &
                                     material_auxvar,reaction,option)
  ! 
  ! Calculates the concentrations of species sorbing
  ! through kinetic sorption processes based
  ! on equilibrium with the aqueous phase.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/09
  ! 

  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscInt :: irate
  
  ! WARNING: below assumes site concentration multiplicative factor
  allocate(rt_auxvar%dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp))
  call RZeroSorb(rt_auxvar)
#if 0  
  !TODO(geh): sort this out
  do irate = 1, reaction%kinmr_nrate
    rt_auxvar%kinmr_total_sorb(:,irate) = reaction%kinmr_frac(irate) * &
                                          rt_auxvar%total_sorb_eq
  enddo
#endif  
  deallocate(rt_auxvar%dtotal_sorb_eq)
  nullify(rt_auxvar%dtotal_sorb_eq)

end subroutine RJumpStartKineticSorption

! ************************************************************************** !

subroutine RReact(rt_auxvar,global_auxvar,material_auxvar,tran_xx_p, &
                  num_iterations_,reaction,option)
  ! 
  ! Solves reaction portion of operator splitting using Newton-Raphson
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/04/10
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: tran_xx_p(reaction%ncomp)
  type(option_type) :: option
  PetscInt :: num_iterations_
  PetscReal :: sign_(reaction%ncomp)
  
  PetscReal :: residual(reaction%ncomp)
  PetscReal :: res(reaction%ncomp)
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  PetscReal :: one_over_dt
  PetscReal :: prev_solution(reaction%ncomp)
  PetscReal :: new_solution(reaction%ncomp)
  PetscReal :: update(reaction%ncomp)
  PetscReal :: maximum_relative_change
  PetscReal :: accumulation_coef
  PetscReal :: fixed_accum(reaction%ncomp)
  PetscInt :: num_iterations
  PetscInt :: icomp
  PetscInt :: immobile_start, immobile_end
  PetscReal :: ratio, min_ratio
  PetscReal :: scale
  
  PetscInt, parameter :: iphase = 1

  one_over_dt = 1.d0/option%tran_dt
  num_iterations = 0

  ! calculate fixed portion of accumulation term
  ! fixed_accum is overwritten in RTAccumulation
  ! Since RTAccumulation uses rt_auxvar%total, we must overwrite the 
  ! rt_auxvar total variables
  ! aqueous
  rt_auxvar%total(:,iphase) = tran_xx_p(1:reaction%naqcomp)
  
  if (reaction%ncoll > 0) then
    option%io_buffer = 'Colloids not set up for operator split mode.'
    call printErrMsg(option)
  endif

! skip chemistry if species nonreacting 
#if 1  
  if (.not.reaction%use_full_geochemistry) then
    rt_auxvar%pri_molal(:) = tran_xx_p(1:reaction%naqcomp) / &
                             global_auxvar%den_kg(iphase)*1.d3
    return
  endif
#endif  
  
  ! update immobile concentrations
  if (reaction%immobile%nimmobile > 0) then
    immobile_start = reaction%offset_immobile + 1
    immobile_end = reaction%offset_immobile + reaction%immobile%nimmobile
    rt_auxvar%immobile(1:reaction%immobile%nimmobile)  = &
      tran_xx_p(immobile_start:immobile_end)
  endif
  
  if (.not.option%use_isothermal) then
    call RUpdateTempDependentCoefs(global_auxvar,reaction,PETSC_FALSE,option)
  endif

  ! still need code to overwrite other phases
  call RTAccumulation(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                      option,fixed_accum)
  if (reaction%neqsorb > 0) then
    call RAccumulationSorb(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                           option,fixed_accum)  
  endif

  ! now update activity coefficients
  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    call RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
  endif
  
  do
  
    num_iterations = num_iterations + 1

    if (reaction%act_coef_update_frequency == &
        ACT_COEF_FREQUENCY_NEWTON_ITER) then
      call RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
    endif
    call RTAuxVarCompute(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
    
    ! Accumulation
    ! residual is overwritten in RTAccumulation()
    call RTAccumulation(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                        option,residual)
    ! J is overwritten in RTAccumulationDerivative()
    call RTAccumulationDerivative(rt_auxvar,global_auxvar,material_auxvar, &
                                  reaction,option,J)
    if (reaction%neqsorb > 0) then
      call RAccumulationSorb(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                             option,residual)
      call RAccumulationSorbDerivative(rt_auxvar,global_auxvar, &
                                       material_auxvar,reaction,option,J)
    endif
    ! must come after sorbed accumulation
    residual = (residual-fixed_accum) / option%tran_dt

                         ! derivative
    call RReaction(residual,J,PETSC_TRUE,rt_auxvar,global_auxvar, &
                   material_auxvar,reaction,option)
    
    if (maxval(abs(residual)) < reaction%max_residual_tolerance) exit

    call RSolve(residual,J,rt_auxvar%pri_molal,update,reaction%ncomp, &
                reaction%use_log_formulation)
    
    prev_solution(1:reaction%naqcomp) = rt_auxvar%pri_molal(1:reaction%naqcomp)
    if (reaction%immobile%nimmobile > 0) then
      prev_solution(immobile_start:immobile_end) = &
        rt_auxvar%immobile(1:reaction%immobile%nimmobile)
    endif

    if (reaction%use_log_formulation) then
      update = dsign(1.d0,update)*min(dabs(update),reaction%max_dlnC)
      new_solution = prev_solution*exp(-update)    
    else ! linear upage
      ! ensure non-negative concentration
      min_ratio = 1.d20 ! large number
      do icomp = 1, reaction%ncomp
        if (prev_solution(icomp) <= update(icomp)) then
          ratio = abs(prev_solution(icomp)/update(icomp))
          if (ratio < min_ratio) min_ratio = ratio
        endif
      enddo
      if (min_ratio < 1.d0) then
        ! scale by 0.99 to make the update slightly smaller than the min_ratio
        update = update*min_ratio*0.99d0
      endif
      new_solution = prev_solution - update
    endif

    maximum_relative_change = maxval(abs((new_solution-prev_solution)/ &
                                         prev_solution))
    
    if (maximum_relative_change < reaction%max_relative_change_tolerance) exit

    if (num_iterations > 50) then
      scale = 1.d0
      if (num_iterations > 50) then
        scale = 0.1d0
      else if (num_iterations > 100) then
        scale = 0.01d0
      else if (num_iterations > 150) then
        scale = 0.001d0
      else if (num_iterations > 500) then
        print *, 'Maximum iterations in RReact: stop: ',num_iterations
        print *, 'Maximum iterations in RReact: residual: ',residual
        print *, 'Maximum iterations in RReact: new solution: ',new_solution
        stop
      endif
      if (scale < 0.99d0) then
        ! apply scaling
        new_solution = scale*(new_solution-prev_solution(:))+prev_solution(:)
      endif
    endif
    
    rt_auxvar%pri_molal(1:reaction%naqcomp) = new_solution(1:reaction%naqcomp)
    if (reaction%immobile%nimmobile > 0) then
      rt_auxvar%immobile(1:reaction%immobile%nimmobile) = &
        new_solution(immobile_start:immobile_end)
    endif    
  
  enddo

  ! one last update
  call RTAuxVarCompute(rt_auxvar,global_auxvar,material_auxvar,reaction,option)

  num_iterations_ = num_iterations
  
end subroutine RReact

! ************************************************************************** !

subroutine RReaction(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                     material_auxvar,reaction,option)
  ! 
  ! Computes reactions
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/30/08
  ! 

  use Option_module

 
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscBool :: derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)

  if (reaction%mineral%nkinmnrl > 0) then
    call RKineticMineral(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                         material_auxvar,reaction,option)
  endif
  

  if (reaction%nradiodecay_rxn > 0) then
    call RRadioactiveDecay(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                           material_auxvar,reaction,option)
  endif
  
  if (reaction%ngeneral_rxn > 0) then
    call RGeneral(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                  material_auxvar,reaction,option)
  endif
  

  if (reaction%immobile%ndecay_rxn > 0) then
    call RImmobileDecay(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                        material_auxvar,reaction,option)
  endif
  
  if (associated(rxn_sandbox_list)) then
    call RSandbox(Res,Jac,derivative,rt_auxvar,global_auxvar, &
                  material_auxvar,reaction,option)
  endif
  
end subroutine RReaction

! ************************************************************************** !

subroutine RReactionDerivative(Res,Jac,rt_auxvar,global_auxvar, &
                               material_auxvar,reaction,option)
  ! 
  ! RReaction: Computes reactions
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/30/08
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  
  implicit none
  
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar 
  type(reactive_transport_auxvar_type) :: rt_auxvar_pert
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
   
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp, joffset
  PetscReal :: Jac_dummy(reaction%ncomp,reaction%ncomp)
  PetscReal :: pert
  PetscBool :: compute_derivative

  ! add new reactions in the 3 locations below

  Res = 0.d0
  Jac = 0.d0

  if (.not.option%transport%numerical_derivatives) then ! analytical derivative
    compute_derivative = PETSC_TRUE
    call RReaction(Res,Jac,compute_derivative,rt_auxvar, &
                   global_auxvar,material_auxvar,reaction,option)  

    ! add only in RReaction

  else ! numerical derivative
    compute_derivative = PETSC_FALSE
    Res_orig = 0.d0
    option%iflag = 0 ! be sure not to allocate mass_balance array
    call RTAuxVarInit(rt_auxvar_pert,reaction,option)
    call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)

    call RReaction(Res_orig,Jac_dummy,compute_derivative,rt_auxvar, &
                   global_auxvar,material_auxvar,reaction,option)     

    ! aqueous species
    do jcomp = 1, reaction%naqcomp
      Res_pert = 0.d0
      call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
      pert = rt_auxvar_pert%pri_molal(jcomp)*perturbation_tolerance
      rt_auxvar_pert%pri_molal(jcomp) = rt_auxvar_pert%pri_molal(jcomp) + pert
      
      call RTotal(rt_auxvar_pert,global_auxvar,material_auxvar,reaction,option)
      call RReaction(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                     global_auxvar,material_auxvar,reaction,option)    

      do icomp = 1, reaction%ncomp
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                           (Res_pert(icomp)-Res_orig(icomp))/pert
      enddo
    enddo
    ! immobile species
    do jcomp = 1, reaction%immobile%nimmobile
      Res_pert = 0.d0
      call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
      ! leave pri_molal, total, total sorbed as is; just copy
      pert = rt_auxvar_pert%immobile(jcomp)*perturbation_tolerance
      rt_auxvar_pert%immobile(jcomp) = rt_auxvar_pert%immobile(jcomp) + pert
      call RReaction(Res_pert,Jac_dummy,compute_derivative,rt_auxvar_pert, &
                     global_auxvar,material_auxvar,reaction,option)    

      ! j is the index in the residual vector and Jacobian
      joffset = reaction%offset_immobile + jcomp
      do icomp = 1, reaction%ncomp
        Jac(icomp,joffset) = Jac(icomp,joffset) + &
                           (Res_pert(icomp)-Res_orig(icomp))/pert
      enddo
    enddo
    
    ! zero small derivatives
    do icomp = 1, reaction%ncomp
      do jcomp = 1, reaction%ncomp
        if (dabs(Jac(icomp,jcomp)) < 1.d-40)  Jac(icomp,jcomp) = 0.d0
      enddo
    enddo
    call RTAuxVarStrip(rt_auxvar_pert)
  endif

end subroutine RReactionDerivative

! ************************************************************************** !

function RSumMoles(rt_auxvar,reaction,option)
  ! 
  ! Sums the total moles of primary and secondary aqueous species
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/01/14
  ! 

  use Option_module
  
  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscReal :: RSumMoles

  PetscInt :: i
  
  RSumMoles = 0.d0
  do i = 1, reaction%naqcomp
    RSumMoles = RSumMoles + rt_auxvar%pri_molal(i)
  enddo
  
  do i = 1, reaction%neqcplx
    RSumMoles = RSumMoles + rt_auxvar%sec_molal(i)
  enddo
  
end function RSumMoles

! ************************************************************************** !

function RCO2MoleFraction(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Sums the total moles of primary and secondary aqueous species
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/01/14
  ! 

  use Option_module
  
  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscReal :: RCO2MoleFraction

  PetscInt :: i
  PetscInt :: icplx
  PetscInt :: ico2
  PetscReal :: sum_co2, sum_mol
  
  if (associated(reaction%species_idx)) then
    ico2 = reaction%species_idx%co2_aq_id
  else
    option%io_buffer = 'reaction%species_idx not set in RCO2MoleFraction(). &
      &Have you defined CO2(aq) as a species in the input deck?'
    call printErrMsg(option)
  endif
  if (ico2 == 0) then
    option%io_buffer = 'CO2 is not set in RCO2MoleFraction(). Have you &
      &defined CO2(aq) as a species in the input deck?'
    call printErrMsg(option)
  endif
  
  sum_co2 = rt_auxvar%pri_molal(ico2)
  sum_mol = RSumMoles(rt_auxvar,reaction,option)
  ! sum_co2 and sum_mol are both in units mol/kg water
  ! FMWH2O is in units g/mol
  ! therefore, scale by 1.d-3 to convert from mol/kg water - g water/mol water
  ! to mol/mol water -- kg water / 1000g water
  RCO2MoleFraction = sum_co2 * FMWH2O * 1.d-3 / &
                     (1.d0 + FMWH2O * sum_mol * 1.d-3)
  
end function RCO2MoleFraction

! ************************************************************************** !

subroutine RActivityCoefficients(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Computes the ionic strength and activity coefficients
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/30/08
  ! 

  use Option_module
  
  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: icplx, icomp, it, j, jcomp, ncomp
  PetscReal :: I, sqrt_I, II, sqrt_II, f, fpri, didi, dcdi, den, dgamdi, &
    lnQK, sum, sum_pri_molal, sum_sec_molal
  PetscReal :: sum_molality
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: NaN

  if (reaction%use_activity_h2o) then
    sum_pri_molal = 0.d0
    do j = 1, reaction%naqcomp
      if (j /= reaction%species_idx%h2o_aq_id) then
        sum_pri_molal = sum_pri_molal + rt_auxvar%pri_molal(j)
      endif
    enddo
  endif

  if (reaction%act_coef_update_algorithm == ACT_COEF_ALGORITHM_NEWTON) then

    ln_conc = log(rt_auxvar%pri_molal)
    ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  
  ! compute primary species contribution to ionic strength
    fpri = 0.d0
    sum_molality = 0.d0
    do j = 1, reaction%naqcomp
      fpri = fpri + rt_auxvar%pri_molal(j)*reaction%primary_spec_Z(j)* &
                                         reaction%primary_spec_Z(j)
    enddo
  
    it = 0
    II = 0
    do
      it = it + 1
      
      if (it > 50) then
        write(option%io_buffer,*) &
          ' too many iterations in computing activity coefficients-stop',it,f,I, &
          ' setting all activity coefficients to NaNs to crash the code.'
        call printErrMsgNoStopByRank(option)
        NaN = 0.d0
        NaN = 1.d0/NaN
        NaN = 0.d0*NaN
        rt_auxvar%pri_molal = NaN
        rt_auxvar%pri_act_coef = NaN
        rt_auxvar%sec_act_coef = NaN
      endif
    
  ! add secondary species contribution to ionic strength
      I = fpri
      do icplx = 1, reaction%neqcplx ! for each secondary species
        I = I + rt_auxvar%sec_molal(icplx)*reaction%eqcplx_Z(icplx)* &
                                         reaction%eqcplx_Z(icplx)
      enddo
      I = 0.5d0*I
      f = I
    
      if (abs(I-II) < 1.d-6*I) exit
    
      if (reaction%neqcplx > 0) then
        didi = 0.d0
        sqrt_I = sqrt(I)
        do icplx = 1, reaction%neqcplx
          if (abs(reaction%eqcplx_Z(icplx)) > 0.d0) then
            sum = 0.5d0*reaction%debyeA*reaction%eqcplx_Z(icplx)* &
            reaction%eqcplx_Z(icplx) &
            /(sqrt_I*(1.d0+reaction%debyeB*reaction%eqcplx_a0(icplx)*sqrt_I)**2) &
            -reaction%debyeBdot
            ncomp = reaction%eqcplxspecid(0,icplx)
            do jcomp = 1, ncomp
              j = reaction%eqcplxspecid(jcomp,icplx)
              if (abs(reaction%primary_spec_Z(j)) > 0.d0) then
                dgamdi = -0.5d0*reaction%debyeA*reaction%primary_spec_Z(j)**2/(sqrt_I* &
                (1.d0+reaction%debyeB*reaction%primary_spec_a0(j)*sqrt_I)**2)+ &
                reaction%debyeBdot 
                sum = sum + reaction%eqcplxstoich(jcomp,icplx)*dgamdi
              endif
            enddo
            dcdi = rt_auxvar%sec_molal(icplx)*LOG_TO_LN*sum
            didi = didi+0.5d0*reaction%eqcplx_Z(icplx)*reaction%eqcplx_Z(icplx)*dcdi
          endif
        enddo
        den = 1.d0-didi
        if (abs(den) > 0.d0) then
          II = (f-I*didi)/den
        else
          II = f
        endif
      else
        II = f
      endif
    
      if (II < 0.d0) then
        write(option%io_buffer,*) 'ionic strength negative! it =',it, &
          ' I= ',I,II,den,didi,dcdi,sum
        call printErrMsgByRank(option)        
      endif
    
  ! compute activity coefficients
  ! primary species
      I = II
      sqrt_I = sqrt(I)
      do icomp = 1, reaction%naqcomp
        if (abs(reaction%primary_spec_Z(icomp)) > 0.d0) then
          rt_auxvar%pri_act_coef(icomp) = exp((-reaction%primary_spec_Z(icomp)* &
                                        reaction%primary_spec_Z(icomp)* &
                                        sqrt_I*reaction%debyeA/ &
                                        (1.d0+reaction%primary_spec_a0(icomp)* &
                                        reaction%debyeB*sqrt_I)+ &
                                        reaction%debyeBdot*I)* &
                                        LOG_TO_LN)
        else
          rt_auxvar%pri_act_coef(icomp) = 1.d0
        endif
      enddo
                
  ! secondary species
      sum_sec_molal = 0.d0
      do icplx = 1, reaction%neqcplx
        if (abs(reaction%eqcplx_Z(icplx)) > 0.d0) then
          rt_auxvar%sec_act_coef(icplx) = exp((-reaction%eqcplx_Z(icplx)* &
                                        reaction%eqcplx_Z(icplx)* &
                                        sqrt_I*reaction%debyeA/ &
                                        (1.d0+reaction%eqcplx_a0(icplx)* &
                                        reaction%debyeB*sqrt_I)+ &
                                        reaction%debyeBdot*I)* &
                                        LOG_TO_LN)
        else
          rt_auxvar%sec_act_coef(icplx) = 1.d0
        endif
    
    ! compute secondary species concentration
        lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

    ! activity of water
        if (reaction%eqcplxh2oid(icplx) > 0) then
          lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
        endif

        ncomp = reaction%eqcplxspecid(0,icplx)
        do jcomp = 1, ncomp
          icomp = reaction%eqcplxspecid(jcomp,icplx)
          lnQK = lnQK + reaction%eqcplxstoich(jcomp,icplx)*ln_act(icomp)
        enddo
        rt_auxvar%sec_molal(icplx) = exp(lnQK)/rt_auxvar%sec_act_coef(icplx)
        sum_sec_molal = sum_sec_molal + rt_auxvar%sec_molal(icplx)
      
      enddo
      
      if (reaction%use_activity_h2o) then
        rt_auxvar%ln_act_h2o = 1.d0-0.017d0*(sum_pri_molal+sum_sec_molal)
        if (rt_auxvar%ln_act_h2o > 0.d0) then
          rt_auxvar%ln_act_h2o = log(rt_auxvar%ln_act_h2o)
        else
          rt_auxvar%ln_act_h2o = 0.d0
          write(option%io_buffer,*) 'activity of H2O negative! ln act H2O =', &
            rt_auxvar%ln_act_h2o
          call printMsg(option)
        endif
      endif

    enddo
  
  else
  
  ! compute ionic strength
  ! primary species
    I = 0.d0
    do icomp = 1, reaction%naqcomp
      I = I + rt_auxvar%pri_molal(icomp)*reaction%primary_spec_Z(icomp)* &
                                       reaction%primary_spec_Z(icomp)
    enddo
  
  ! secondary species
    do icplx = 1, reaction%neqcplx ! for each secondary species
      I = I + rt_auxvar%sec_molal(icplx)*reaction%eqcplx_Z(icplx)* &
                                       reaction%eqcplx_Z(icplx)
    enddo
    I = 0.5d0*I
    sqrt_I = sqrt(I)
  
  ! compute activity coefficients
  ! primary species
    do icomp = 1, reaction%naqcomp
      if (abs(reaction%primary_spec_Z(icomp)) > 1.d-10) then
        rt_auxvar%pri_act_coef(icomp) = exp((-reaction%primary_spec_Z(icomp)* &
                                      reaction%primary_spec_Z(icomp)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%primary_spec_a0(icomp)* &
                                      reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                      LOG_TO_LN)
      else
        rt_auxvar%pri_act_coef(icomp) = 1.d0
      endif
    enddo
                
  ! secondary species
    sum_sec_molal = 0.d0
    do icplx = 1, reaction%neqcplx
      if (dabs(reaction%eqcplx_Z(icplx)) > 1.d-10) then
        rt_auxvar%sec_act_coef(icplx) = exp((-reaction%eqcplx_Z(icplx)* &
                                      reaction%eqcplx_Z(icplx)* &
                                      sqrt_I*reaction%debyeA/ &
                                      (1.d0+reaction%eqcplx_a0(icplx)* &
                                      reaction%debyeB*sqrt_I)+ &
                                      reaction%debyeBdot*I)* &
                                      LOG_TO_LN)
      else
        rt_auxvar%sec_act_coef(icplx) = 1.d0
      endif
      sum_sec_molal = sum_sec_molal + rt_auxvar%sec_molal(icplx)
    enddo
    
    if (reaction%use_activity_h2o) then
      rt_auxvar%ln_act_h2o = 1.d0-0.017d0*(sum_pri_molal+sum_sec_molal)
      if (rt_auxvar%ln_act_h2o > 0.d0) then
        rt_auxvar%ln_act_h2o = log(rt_auxvar%ln_act_h2o)
      else
        rt_auxvar%ln_act_h2o = 0.d0
      endif
    endif
  endif
  
end subroutine RActivityCoefficients

! ************************************************************************** !

subroutine RTotal(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Computes the total component concentrations and derivative with
  ! respect to free-ion for all phases (aqueous, sorbed, gas, etc.)
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/21/18
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
  if (reaction%neqsorb > 0) then
    call RTotalSorb(rt_auxvar,global_auxvar,material_auxvar, &
                    reaction,option)
  endif
  if (reaction%gas%nactive_gas > 0) then
    call RTotalGas(rt_auxvar,global_auxvar,reaction,option)
  endif



end subroutine RTotal

! ************************************************************************** !

subroutine RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Computes the total component concentrations and derivative with
  ! respect to free-ion for the aqueous phase
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/28/08
  ! 

  use Option_module

  use EOS_Water_module
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, icplx, icomp, jcomp, ncomp
  PetscInt, parameter :: iphase = 1
  PetscErrorCode :: ierr
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal
  PetscReal :: den_kg_per_L, xmass
  
  rt_auxvar%total = 0.d0 !debugging 
  
  xmass = 1.d0
  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  den_kg_per_L = global_auxvar%den_kg(iphase)*xmass*1.d-3

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  rt_auxvar%total(:,iphase) = rt_auxvar%pri_molal(:)
  
  ! initialize derivatives
  rt_auxvar%aqueous%dtotal = 0.d0
  do icomp = 1, reaction%naqcomp
    rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) = 1.d0
  enddo
   
  do icplx = 1, reaction%neqcplx ! for each secondary species
    ! compute secondary species concentration
    lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

    ! activity of water
    if (reaction%eqcplxh2oid(icplx) > 0) then
      lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
    endif

    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      lnQK = lnQK + reaction%eqcplxstoich(i,icplx)*ln_act(icomp)
    enddo
    rt_auxvar%sec_molal(icplx) = exp(lnQK)/rt_auxvar%sec_act_coef(icplx)
  
    ! add contribution to primary totals
    ! units of total here are mol/kg water (molality) since we are summing 
    ! molalities, but it will be coverted to molarity [mol/L] below     
    do i = 1, ncomp
      icomp = reaction%eqcplxspecid(i,icplx)
      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
                                      reaction%eqcplxstoich(i,icplx)* &
                                      rt_auxvar%sec_molal(icplx)
    enddo
    
    ! add contribution to derivatives of total with respect to free
    ! bear in mind that the water density portion is scaled below
    ! dtotal units here are [mol/kg water * kg water/mol] = [-]
    do j = 1, ncomp
      jcomp = reaction%eqcplxspecid(j,icplx)
      tempreal = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 rt_auxvar%sec_act_coef(icplx)
      do i = 1, ncomp
        icomp = reaction%eqcplxspecid(i,icplx)
        rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) = &
          rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) + &
          reaction%eqcplxstoich(i,icplx)*tempreal
      enddo
    enddo
  enddo

  ! convert molality -> molarity
  ! unit of total = mol/L water
  rt_auxvar%total(:,iphase) = rt_auxvar%total(:,iphase)*den_kg_per_L
  
  ! units of dtotal = kg water/L water
  rt_auxvar%aqueous%dtotal = rt_auxvar%aqueous%dtotal*den_kg_per_L

end subroutine RTotalAqueous

! ************************************************************************** !

subroutine RZeroSorb(rt_auxvar)
  ! 
  ! Zeros out arrays associated with sorption
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/12
  ! 

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  
  if (associated(rt_auxvar%total_sorb_eq)) rt_auxvar%total_sorb_eq = 0.d0
  if (associated(rt_auxvar%dtotal_sorb_eq)) rt_auxvar%dtotal_sorb_eq = 0.d0
  if (associated(rt_auxvar%eqsrfcplx_conc)) rt_auxvar%eqsrfcplx_conc = 0.d0
  
end subroutine RZeroSorb

! ************************************************************************** !

subroutine RTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Computes the total sorbed component concentrations and
  ! derivative with respect to free-ion
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/08
  ! 

  use Option_module
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  call RZeroSorb(rt_auxvar)

  if (reaction%neqionxrxn > 0) then
    call RTotalSorbEqIonx(rt_auxvar,global_auxvar,reaction,option)
  endif
  
  if (reaction%neqkdrxn > 0) then
    call RTotalSorbKD(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
  endif
  
end subroutine RTotalSorb

! ************************************************************************** !

subroutine RTotalSorbKD(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                        option)
  ! 
  ! Computes the total sorbed component concentrations and
  ! derivative with respect to free-ion for the linear
  ! K_D model
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/30/2010
  ! 

  use Option_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: irxn
  PetscInt :: icomp
  PetscReal :: res
  PetscReal :: dres_dc
  PetscReal :: molality
  PetscReal :: tempreal
  PetscReal :: one_over_n
  PetscReal :: molality_one_over_n
  PetscReal :: kd_kgw_m3b  
  PetscReal :: temp

  PetscInt, parameter :: iphase = 1

  do irxn = 1, reaction%neqkdrxn
    icomp = reaction%eqkdspecid(irxn)
    molality = rt_auxvar%pri_molal(icomp)
    if (reaction%eqkdmineral(irxn) > 0) then
      ! NOTE: mineral volume fraction here is solely a scaling factor.  It has 
      ! nothing to do with the soil volume; that is calculated through as a 
      ! function of porosity.
      temp = reaction%eqkddistcoef(irxn)
      temp = global_auxvar%den_kg(iphase)
      temp = (1.d0-material_auxvar%porosity)
      temp = material_auxvar%soil_particle_density
      temp = (rt_auxvar%mnrl_volfrac(reaction%eqkdmineral(irxn)))
      kd_kgw_m3b = reaction%eqkddistcoef(irxn) * & !KD units [mL water/g soil]
                   global_auxvar%den_kg(iphase) * &
                   (1.d0-material_auxvar%porosity) * &
                   material_auxvar%soil_particle_density * &
                   1.d-3 * & ! convert mL water/g soil to m^3 water/kg soil
                   (rt_auxvar%mnrl_volfrac(reaction%eqkdmineral(irxn)))
    else
      kd_kgw_m3b = reaction%eqkddistcoef(irxn)
    endif
    select case(reaction%eqkdtype(irxn))
      case(SORPTION_LINEAR)
        ! Csorb = Kd*Caq
        res = kd_kgw_m3b*molality
        dres_dc = kd_kgw_m3b
      case(SORPTION_LANGMUIR)
        ! Csorb = K*Caq*b/(1+K*Caq)
        tempreal = kd_kgw_m3b*molality
        res = tempreal*reaction%eqkdlangmuirb(irxn) / (1.d0 + tempreal)
        dres_dc = res/molality - &
                  res / (1.d0 + tempreal) * tempreal / molality
      case(SORPTION_FREUNDLICH)
        ! Csorb = Kd*Caq**(1/n)
        one_over_n = 1.d0/reaction%eqkdfreundlichn(irxn)
        molality_one_over_n = molality**one_over_n
        res = kd_kgw_m3b*molality**one_over_n
        dres_dc = res/molality*one_over_n
      case default
        res = 0.d0
        dres_dc = 0.d0
    end select
    rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + res
    rt_auxvar%dtotal_sorb_eq(icomp,icomp) = &
      rt_auxvar%dtotal_sorb_eq(icomp,icomp) + dres_dc 
  enddo

end subroutine RTotalSorbKD

! ************************************************************************** !

subroutine RTotalSorbEqIonx(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Computes the total sorbed component concentrations and
  ! derivative with respect to free-ion for equilibrium ion
  ! exchange
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/08; 05/26/09
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, iphase, ncomp, ncplx
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: tempreal, tempreal1, tempreal2, total
  PetscInt :: irxn
  PetscReal, parameter :: tol = 1.d-12
  PetscBool :: one_more
  PetscReal :: res
    
  PetscReal :: omega
  PetscReal :: ref_cation_X, ref_cation_conc, ref_cation_Z, ref_cation_k, &
               ref_cation_quotient
  PetscReal :: cation_X(reaction%naqcomp)
  PetscReal :: dres_dref_cation_X, dref_cation_X
  PetscReal :: sumZX, sumkm
    
  PetscReal :: total_pert, ref_cation_X_pert, pert
  PetscReal :: ref_cation_quotient_pert, dres_dref_cation_X_pert

  PetscReal :: KDj, dres_dKDj, delta_KDj
  PetscInt :: it

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
    
  ! Ion Exchange
  if (associated(rt_auxvar%eqionx_conc)) rt_auxvar%eqionx_conc = 0.d0
  do irxn = 1, reaction%neqionxrxn

    ncomp = reaction%eqionx_rxn_cationid(0,irxn)

    ! for now we assume that omega is equal to CEC.
    if (reaction%eqionx_rxn_to_surf(irxn) > 0) then
      ! if tied to a mineral vol frac
      omega = max(reaction%eqionx_rxn_CEC(irxn)* &
                  rt_auxvar%mnrl_volfrac(reaction%eqionx_rxn_to_surf(irxn)), &
                  1.d-40)
    else
      omega = reaction%eqionx_rxn_CEC(irxn)
    endif

    if (reaction%eqionx_rxn_Z_flag(irxn)) then ! Zi /= Zj for any i,j

      icomp = reaction%eqionx_rxn_cationid(1,irxn)
      ref_cation_conc = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)
      ref_cation_Z = reaction%primary_spec_Z(icomp)
      ref_cation_k = reaction%eqionx_rxn_k(1,irxn)
      ref_cation_X = ref_cation_Z* &
                     rt_auxvar%eqionx_ref_cation_sorbed_conc(irxn)/omega

      one_more = PETSC_FALSE
      cation_X = 0.d0
      KDj = ref_cation_X /(ref_cation_k*ref_cation_conc)
      it = 0
!geh: Change from 0 to 1 to run new implementation.
#if 1
      do
        it = it + 1
        if (it > 20000) then
          option%io_buffer = 'Too many Newton iterations in ion exchange.'
          call printErrMsgByRank(option)
        endif
        ref_cation_X = KDj*(ref_cation_k*ref_cation_conc)
        cation_X(1) = ref_cation_X
        total = ref_cation_X
        dres_dKDj = 0.d0
        do j = 2, ncomp
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          cation_X(j) = reaction%eqionx_rxn_k(j,irxn)* &
                        rt_auxvar%pri_molal(icomp)* &
                        rt_auxvar%pri_act_coef(icomp)* &
                        KDj**(reaction%primary_spec_Z(icomp)/ref_cation_Z)
          total = total + cation_X(j)
          dres_dKDj = dres_dKDj + cation_X(j)/KDj* &
                                  reaction%primary_spec_Z(icomp)
        enddo
        dres_dKDj = dres_dKDj/ref_cation_Z + (ref_cation_k*ref_cation_conc)
        res = 1.d0 - total

        if (one_more) exit

        ! no need to negate since res is subtracted above.
        delta_KDj = res/dres_dKDj
        KDj = KDj + delta_KDj
        KDj = max(KDj,1.d-40) ! prevent from going negative
        if (dabs(delta_KDj/KDj) < tol) then
          one_more = PETSC_TRUE
        endif
      enddo
#else 
      do
        if (ref_cation_X <= 0.d0) ref_cation_X = 1.d-8
        cation_X(1) = ref_cation_X
        ref_cation_quotient = ref_cation_X/(ref_cation_k*ref_cation_conc)
        total = ref_cation_X

        do j = 2, ncomp
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          cation_X(j) = rt_auxvar%pri_molal(icomp)* &
                        rt_auxvar%pri_act_coef(icomp)* &
                        reaction%eqionx_rxn_k(j,irxn)* &
                        ref_cation_quotient** &
                        (reaction%primary_spec_Z(icomp)/ref_cation_Z)
          total = total + cation_X(j)
        enddo
        
        if (one_more) exit
        
        res = 1.d0-total
          
        dres_dref_cation_X = 1.d0

#if 0
  ! test derivative
        pert = 1.d-6 * ref_cation_X
        ref_cation_X_pert = ref_cation_X + pert
        ref_cation_quotient_pert = ref_cation_X_pert/ &
        (ref_cation_k*ref_cation_conc)
        total_pert = ref_cation_X_pert

          do j = 2, ncomp
            icomp = reaction%eqionx_rxn_cationid(j,irxn)
            total_pert = total_pert + &
                         rt_auxvar%pri_molal(icomp)* &
                         rt_auxvar%pri_act_coef(icomp)* &
                         reaction%eqionx_rxn_k(j,irxn)* &
                         ref_cation_quotient_pert** &
                         (reaction%primary_spec_Z(icomp)/ref_cation_Z)
          enddo
        dres_dref_cation_X_pert = (1.d0-total_pert-res)/pert
  ! test
#endif

        do j = 2, ncomp
          icomp = reaction%eqionx_rxn_cationid(j,irxn)
          dres_dref_cation_X = dres_dref_cation_X + &
            (reaction%primary_spec_Z(icomp)/ref_cation_Z)* &
            cation_X(j)/ref_cation_X
        enddo

        dref_cation_X = res / (-dres_dref_cation_X)
!        dref_cation_X = res / dres_dref_cation_X_pert
        ref_cation_X = ref_cation_X - dref_cation_X
      
        if (dabs(dref_cation_X/ref_cation_X) < tol) then
          one_more = PETSC_TRUE
        endif
      
      enddo
#endif

      rt_auxvar%eqionx_ref_cation_sorbed_conc(irxn) = ref_cation_X*omega/ &
                                                      ref_cation_Z

    else ! Zi == Zj for all i,j
        
      sumkm = 0.d0
      cation_X = 0.d0
      
      do j = 1, ncomp  
        icomp = reaction%eqionx_rxn_cationid(j,irxn)
        cation_X(j) = rt_auxvar%pri_molal(icomp)* &
                      rt_auxvar%pri_act_coef(icomp)* &
                      reaction%eqionx_rxn_k(j,irxn)
        sumkm = sumkm + cation_X(j)
      enddo
          
      cation_X = cation_X / sumkm

    endif
                
    ! sum up charges
    sumZX = 0.d0
    do i = 1, ncomp
      icomp = reaction%eqionx_rxn_cationid(i,irxn)
      sumZX = sumZX + reaction%primary_spec_Z(icomp)*cation_X(i)
    enddo

    ! compute totals based on sorbed ions
    do i = 1, ncomp
      icomp = reaction%eqionx_rxn_cationid(i,irxn)
      tempreal1 = cation_X(i)*omega/reaction%primary_spec_Z(icomp)
      ! residual function entry
      
      rt_auxvar%eqionx_conc(i,irxn) = rt_auxvar%eqionx_conc(i,irxn) + tempreal1

      rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + &
                                       tempreal1

      tempreal2 = reaction%primary_spec_Z(icomp)/sumZX
      do j = 1, ncomp
        jcomp = reaction%eqionx_rxn_cationid(j,irxn)
        if (i == j) then
          rt_auxvar%dtotal_sorb_eq(icomp,jcomp) = &
            rt_auxvar%dtotal_sorb_eq(icomp,jcomp) + &
            tempreal1*(1.d0-(tempreal2*cation_X(j)))/ &
            rt_auxvar%pri_molal(jcomp)
        else
          rt_auxvar%dtotal_sorb_eq(icomp,jcomp) = &
            rt_auxvar%dtotal_sorb_eq(icomp,jcomp) + &
            (-tempreal1)*tempreal2*cation_X(j)/ &
            rt_auxvar%pri_molal(jcomp)
        endif
      enddo
    enddo    

  enddo
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
end subroutine RTotalSorbEqIonx

! ************************************************************************** !

subroutine RAccumulationSorb(rt_auxvar,global_auxvar,material_auxvar, &
                             reaction,option,Res)
  ! 
  ! Computes non-aqueous portion of the accumulation term in
  ! residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/26/09
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  
!  PetscReal :: v_t
  
  ! units = (mol solute/m^3 bulk)*(m^3 bulk)/(sec) = mol/sec
  ! all residual entries should be in mol/sec
!  v_t = material_auxvar%volume/option%tran_dt
  Res(1:reaction%naqcomp) = Res(1:reaction%naqcomp) + &
    rt_auxvar%total_sorb_eq(:)*material_auxvar%volume

end subroutine RAccumulationSorb

! ************************************************************************** !

subroutine RAccumulationSorbDerivative(rt_auxvar,global_auxvar, &
                                       material_auxvar,reaction,option,J)
  ! 
  ! Computes derivative of non-aqueous portion of
  ! the accumulation term in residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/26/09
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp
  PetscReal :: v_t
  
  ! units = (kg water/m^3 bulk)*(m^3 bulk)/(sec) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  v_t = material_auxvar%volume/option%tran_dt
  J(1:reaction%naqcomp,1:reaction%naqcomp) = &
    J(1:reaction%naqcomp,1:reaction%naqcomp) + &
    rt_auxvar%dtotal_sorb_eq(:,:)*v_t

end subroutine RAccumulationSorbDerivative

! ************************************************************************** !

subroutine RRadioactiveDecay(Res,Jac,compute_derivative,rt_auxvar, &
                             global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Computes radioactive decay with a single reactant
  ! (considering both the aqueous and sorbed phases) with
  ! the possibility of multiple daughter products
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/08/10, 01/07/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  PetscInt :: i, icomp, jcomp, irxn, ncomp
  PetscReal :: tempreal, L_water, sum, rate

  PetscInt, parameter :: iphase = 1

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! L water

  do irxn = 1, reaction%nradiodecay_rxn ! for each reaction
    
    ! units(kf): 1/sec
    
    ! we assume only one chemical component involved in decay reaction
    icomp = reaction%radiodecayforwardspecid(irxn)

    ! sum total moles of component in aqueous and sorbed phases
    sum = rt_auxvar%total(icomp,iphase)*L_water
    if (associated(rt_auxvar%total_sorb_eq)) then
      sum = sum + rt_auxvar%total_sorb_eq(icomp)*material_auxvar%volume
    endif
    
    rate = sum*reaction%radiodecay_kf(irxn)
    
    ! units(Res): mol/sec
    ncomp = reaction%radiodecayspecid(0,irxn)
    do i = 1, ncomp
      icomp = reaction%radiodecayspecid(i,irxn)
      ! units = mol/sec
      Res(icomp) = Res(icomp) - reaction%radiodecaystoich(i,irxn)*rate
    enddo    

    if (.not. compute_derivative) cycle   

    tempreal = -1.d0*reaction%radiodecay_kf(irxn)
    jcomp = reaction%radiodecayforwardspecid(irxn)
    if (associated(rt_auxvar%dtotal_sorb_eq)) then
      do i = 1, ncomp
        icomp = reaction%radiodecayspecid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,1:reaction%naqcomp) = Jac(icomp,1:reaction%naqcomp) + &
          tempreal * &
          reaction%radiodecaystoich(i,irxn) * &
          (rt_auxvar%aqueous%dtotal(jcomp,1:reaction%naqcomp,iphase)*L_water + &
           rt_auxvar%dtotal_sorb_eq(jcomp,1:reaction%naqcomp)* &
           material_auxvar%volume)
      enddo
    else ! no sorption
      do i = 1, ncomp
        icomp = reaction%radiodecayspecid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,1:reaction%naqcomp) = Jac(icomp,1:reaction%naqcomp) + &
          tempreal * &
          reaction%radiodecaystoich(i,irxn) * &
          rt_auxvar%aqueous%dtotal(jcomp,1:reaction%naqcomp,iphase)*L_water
      enddo
    endif
    
  enddo  ! loop over reactions
    
end subroutine RRadioactiveDecay

! ************************************************************************** !

subroutine RGeneral(Res,Jac,compute_derivative,rt_auxvar,global_auxvar, &
                    material_auxvar,reaction,option)
  ! 
  ! Computes the general reaction rates
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/08/10
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  
  PetscInt :: i, j, icomp, jcomp, irxn, ncomp
  PetscReal :: tempreal
  PetscReal :: kf, kr
  PetscReal :: Qkf, lnQkf
  PetscReal :: Qkr, lnQkr
  PetscReal :: por_den_sat_vol

  PetscInt, parameter :: iphase = 1

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  do irxn = 1, reaction%ngeneral_rxn ! for each mineral
    
    ! units
    ! for nth-order reaction
    ! kf/kr = kg^(n-1)/mol^(n-1)-sec
    ! thus for a 1st-order reaction, kf units = 1/sec
    
    kf = reaction%general_kf(irxn)
    kr = reaction%general_kr(irxn)

    if (kf > 0.d0) then
      ! compute ion activity product
      lnQkf = log(kf)

  ! currently not accommodating activity of water
      ! activity of water
  !    if (reaction%kinmnrlh2oid(irxn) > 0) then
  !      lnQkf = lnQkf + reaction%generalh2ostoich(irxn)*rt_auxvar%ln_act_h2o
  !    endif

      ncomp = reaction%generalforwardspecid(0,irxn)
      do i = 1, ncomp
        icomp = reaction%generalforwardspecid(i,irxn)
        lnQkf = lnQkf + reaction%generalforwardstoich(i,irxn)*ln_act(icomp)
      enddo
      Qkf = exp(lnQkf)
    else
      Qkf = 0.d0
    endif
    
    if (kr > 0.d0) then
      lnQkr = log(kr)

  ! currently not accommodating activity of water
      ! activity of water
  !    if (reaction%kinmnrlh2oid(irxn) > 0) then
  !      lnQkr = lnQkr + reaction%generalh2ostoich(irxn)*rt_auxvar%ln_act_h2o
  !    endif

      ncomp = reaction%generalbackwardspecid(0,irxn)
      do i = 1, ncomp
        icomp = reaction%generalbackwardspecid(i,irxn)
        lnQkr = lnQkr + reaction%generalbackwardstoich(i,irxn)*ln_act(icomp)
      enddo
      Qkr = exp(lnQkr)
    else
      Qkr = 0.d0
    endif

    ! Qkf/Qkr units are now mol/kg(water)-sec

    por_den_sat_vol = material_auxvar%porosity*global_auxvar%den_kg(iphase)* &
                      global_auxvar%sat(iphase)* &
                      material_auxvar%volume

    ncomp = reaction%generalspecid(0,irxn)
    do i = 1, ncomp
      icomp = reaction%generalspecid(i,irxn)
      ! units = mol/sec
      Res(icomp) = Res(icomp) - reaction%generalstoich(i,irxn)*(Qkf-Qkr)* &
                                por_den_sat_vol
    enddo 

    if (.not. compute_derivative) cycle   

    ! calculate derivatives of rate with respect to free
    ! units = mol/sec

    if (kf > 0.d0) then
      ! derivatives with respect to primary species in forward reaction
      do j = 1, reaction%generalforwardspecid(0,irxn)
        jcomp = reaction%generalforwardspecid(j,irxn)
        tempreal = -1.d0*reaction%generalforwardstoich(j,irxn)*exp(lnQkf-ln_conc(jcomp))* &
                   por_den_sat_vol
        do i = 1, reaction%generalspecid(0,irxn)
          icomp = reaction%generalspecid(i,irxn)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                             reaction%generalstoich(i,irxn)*tempreal
        enddo
      enddo
    endif
    
    if (kr > 0.d0) then
      ! derivatives with respect to primary species in forward reaction
      do j = 1, reaction%generalbackwardspecid(0,irxn)
        jcomp = reaction%generalbackwardspecid(j,irxn)
        tempreal = reaction%generalbackwardstoich(j,irxn)*exp(lnQkr-ln_conc(jcomp))* &
                   por_den_sat_vol
        do i = 1, reaction%generalspecid(0,irxn)
          icomp = reaction%generalspecid(i,irxn)
          ! units = (mol/sec)*(kg water/mol) = kg water/sec
          Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                             reaction%generalstoich(i,irxn)*tempreal
        enddo
      enddo
    endif

  enddo  ! loop over reactions
    
end subroutine RGeneral

! ************************************************************************** !

subroutine RSolve(Res,Jac,conc,update,ncomp,use_log_formulation)
  ! 
  ! Computes the kinetic mineral precipitation/dissolution
  ! rates
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/04/08
  ! 

  use Utility_module
  
  implicit none

  PetscInt :: ncomp
  PetscReal :: Res(ncomp)
  PetscReal :: Jac(ncomp,ncomp)
  PetscReal :: update(ncomp)
  PetscReal :: conc(ncomp)
  PetscBool :: use_log_formulation
  
  PetscInt :: indices(ncomp)
  PetscReal :: rhs(ncomp)
  PetscInt :: icomp
  PetscReal :: norm

  ! scale Jacobian
  do icomp = 1, ncomp
    norm = max(1.d0,maxval(abs(Jac(icomp,:))))
    norm = 1.d0/norm
    rhs(icomp) = Res(icomp)*norm
    Jac(icomp,:) = Jac(icomp,:)*norm
  enddo
    
  if (use_log_formulation) then
    ! for derivatives with respect to ln conc
    do icomp = 1, ncomp
      Jac(:,icomp) = Jac(:,icomp)*conc(icomp)
    enddo
  endif

  call ludcmp(Jac,ncomp,indices,icomp)
  call lubksb(Jac,ncomp,indices,rhs)
  
  update = rhs

end subroutine RSolve

! ************************************************************************** !

subroutine ReactionComputeKd(icomp,retardation,rt_auxvar,global_auxvar, &
                             material_auxvar,reaction,option)
  ! 
  ! RComputeKd: Computes the Kd for a given chemical component
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/14/09
  ! 

  use Option_module
  
  implicit none
  
  PetscInt :: icomp
  PetscReal :: retardation
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscReal :: bulk_vol_to_fluid_vol
  PetscInt :: i, j, jcomp, irxn, icplx, irate
  PetscInt, parameter :: iphase = 1

  retardation = 0.d0
  if (reaction%nsorb == 0) return
  
  bulk_vol_to_fluid_vol = material_auxvar%porosity* &
                          global_auxvar%sat(iphase)*1000.d0

  if (associated(rt_auxvar%total_sorb_eq)) then
    retardation = rt_auxvar%total_sorb_eq(icomp)
  endif

  if (dabs(rt_auxvar%total(icomp,iphase)) > 1.d-40) &
    retardation = retardation/bulk_vol_to_fluid_vol/ &
      rt_auxvar%total(icomp,iphase)

end subroutine ReactionComputeKd

! ************************************************************************** !

subroutine RAge(rt_auxvar,global_auxvar,material_auxvar,option,reaction,Res)
  ! 
  ! Computes the ages of the groundwater
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/10
  ! 

  use Option_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar 
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscInt, parameter :: iphase = 1
  
  Res(:) = 0.d0
  if (reaction%calculate_water_age) then
    Res(reaction%species_idx%water_age_id) = material_auxvar%porosity* &
                                             global_auxvar%sat(iphase)* &
                                             1000.d0 * material_auxvar%volume
  endif
  if (reaction%calculate_tracer_age) then
    Res(reaction%species_idx%tracer_age_id) = &
      -rt_auxvar%total(reaction%species_idx%tracer_aq_id,iphase)* &
      material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0* &
      material_auxvar%volume
  endif
end subroutine RAge

! ************************************************************************** !

subroutine RTAuxVarCompute(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                           option)
  ! 
  ! Computes secondary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/28/08
  ! 

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
#if 0  
  PetscReal :: Res_orig(reaction%ncomp)
  PetscReal :: Res_pert(reaction%ncomp)
  PetscInt :: icomp, jcomp, iphase
  PetscReal :: dtotal(reaction%naqcomp,reaction%naqcomp,2)
  PetscReal :: dtotalsorb(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: pert
  PetscInt :: nphase
  type(reactive_transport_auxvar_type) :: rt_auxvar_pert
#endif

  !already set  rt_auxvar%pri_molal = x
  call RTotal(rt_auxvar,global_auxvar,material_auxvar,reaction,option)

#if 0
! numerical check
  nphase = reaction%nphase
  Res_orig = 0.d0
  dtotal = 0.d0
  dtotalsorb = 0.d0
  option%iflag = 0 ! be sure not to allocate mass_balance array
  call RTAuxVarInit(rt_auxvar_pert,reaction,option)
  do jcomp = 1, reaction%naqcomp
    Res_pert = 0.d0
    call RTAuxVarCopy(rt_auxvar_pert,rt_auxvar,option)
    if (reaction%neqcplx > 0) then
      rt_auxvar%sec_molal = 0.d0
    endif
    if (reaction%gas%nactive_gas > 0) then
      rt_auxvar%gas_pp = 0.d0
    endif
    if (reaction%surface_complexation%neqsrfcplxrxn > 0) then
      rt_auxvar_pert%srfcplxrxn_free_site_conc = 1.d-9
      rt_auxvar_pert%eqsrfcplx_conc = 0.d0
    endif
    if (reaction%neqionxrxn > 0) then
      rt_auxvar%eqionx_ref_cation_sorbed_conc = 1.d-9
    endif
    pert = rt_auxvar_pert%pri_molal(jcomp)*perturbation_tolerance
    rt_auxvar_pert%pri_molal(jcomp) = rt_auxvar_pert%pri_molal(jcomp) + pert
    
    call RTotal(rt_auxvar_pert,global_auxvar,material_auxvar,reaction,option)
    dtotal(:,jcomp,1) = (rt_auxvar_pert%total(:,1) - &
                         rt_auxvar%total(:,1))/pert
    if (reaction%neqsorb > 0) then
      dtotalsorb(:,jcomp) = (rt_auxvar_pert%total_sorb_eq(:) - &
                             rt_auxvar%total_sorb_eq(:))/pert
    endif
    if (option%iflowmode == MPH_MODE .or. &
        option%iflowmode == IMS_MODE .or. &
        option%iflowmode == FLASH2_MODE .or. 
        reaction%gas%nactive_gas > 0) then
      dtotal(:,jcomp,2) = (rt_auxvar_pert%total(:,2) - &
                           rt_auxvar%total(:,2))/pert
    endif    
  enddo
  do iphase = 1, nphase
    do icomp = 1, reaction%naqcomp
      do jcomp = 1, reaction%naqcomp
        if (dabs(dtotal(icomp,jcomp,iphase)) < 1.d-16) &
          dtotal(icomp,jcomp,iphase) = 0.d0
      enddo
      if (reaction%neqsorb > 0) then
        if (dabs(dtotalsorb(icomp,jcomp)) < 1.d-16) &
          dtotalsorb(icomp,jcomp) = 0.d0
      endif
    enddo
  enddo
  rt_auxvar%aqueous%dtotal(:,:,:) = dtotal
  if (reaction%neqsorb > 0) rt_auxvar%dtotal_sorb_eq = dtotalsorb
  call RTAuxVarStrip(rt_auxvar_pert)
#endif
  
end subroutine RTAuxVarCompute

! ************************************************************************** !

subroutine RTAccumulation(rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option,Res)
  ! 
  ! Computes aqueous portion of the accumulation term in
  ! residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  
  PetscInt :: iphase
  PetscInt :: istart, iend
  PetscInt :: idof
  PetscInt :: iimob
  PetscInt :: icoll
  PetscInt :: icollcomp
  PetscInt :: iaqcomp
  PetscInt :: iimb
  PetscReal :: psv_t
  
  iphase = 1
  Res = 0.d0
  
  ! units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  !         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  ! 1000.d0 converts vol from m^3 -> L
  ! all residual entries should be in mol/sec
  psv_t = material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0* &
!          material_auxvar%volume / option%tran_dt  
          material_auxvar%volume
  istart = 1
  iend = reaction%naqcomp
  Res(istart:iend) = psv_t*rt_auxvar%total(:,iphase) 

  if (reaction%immobile%nimmobile > 0) then
    do iimob = 1, reaction%immobile%nimmobile
      idof = reaction%offset_immobile + iimob
      Res(idof) = Res(idof) + rt_auxvar%immobile(iimob)* &
!                              material_auxvar%volume/option%tran_dt 
                              material_auxvar%volume
    enddo
  endif
  if (reaction%ncoll > 0) then
    do icoll = 1, reaction%ncoll
      idof = reaction%offset_colloid + icoll
      Res(idof) = psv_t*rt_auxvar%colloid%conc_mob(icoll)
    enddo
  endif
  if (reaction%ncollcomp > 0) then
    do icollcomp = 1, reaction%ncollcomp
      iaqcomp = reaction%coll_spec_to_pri_spec(icollcomp)
      Res(iaqcomp) = Res(iaqcomp) + &
        psv_t*rt_auxvar%colloid%total_eq_mob(icollcomp)
    enddo
  endif
  if (reaction%gas%nactive_gas > 0) then
    iphase = 2
    psv_t = material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0* &
!            material_auxvar%volume / option%tran_dt  
            material_auxvar%volume
    istart = 1
    iend = reaction%naqcomp
    Res(istart:iend) = Res(istart:iend) + psv_t*rt_auxvar%total(:,iphase)     
  endif

end subroutine RTAccumulation

! ************************************************************************** !

subroutine RTAccumulationDerivative(rt_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    reaction,option,J)
  ! 
  ! Computes derivative of aqueous portion of the
  ! accumulation term in residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar  
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  
  PetscInt :: icomp, iphase
  PetscInt :: istart, iendaq
  PetscInt :: idof
  PetscInt :: icoll
  PetscInt :: iimob
  PetscReal :: psvd_t ! d is for density

  iphase = 1
  istart = 1
  iendaq = reaction%naqcomp 
  ! units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  !         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  ! all Jacobian entries should be in kg water/sec
  J = 0.d0
  if (associated(rt_auxvar%aqueous%dtotal)) then ! units of dtotal = kg water/L water
    psvd_t = material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0* &
             material_auxvar%volume/option%tran_dt
    J(istart:iendaq,istart:iendaq) = rt_auxvar%aqueous%dtotal(:,:,iphase)*psvd_t
  else
    psvd_t = material_auxvar%porosity*global_auxvar%sat(iphase)* &
             global_auxvar%den_kg(iphase)*material_auxvar%volume/ &
             option%tran_dt ! units of den = kg water/m^3 water
    do icomp=istart,iendaq
      J(icomp,icomp) = psvd_t
    enddo
  endif

  if (reaction%immobile%nimmobile > 0) then
    do iimob = 1, reaction%immobile%nimmobile
      idof = reaction%offset_immobile + iimob
      J(idof,idof) = material_auxvar%volume/option%tran_dt
    enddo
  endif
  if (reaction%ncoll > 0) then
    do icoll = 1, reaction%ncoll
      idof = reaction%offset_colloid + icoll
      ! shouldn't have to sum a this point
      J(idof,idof) = psvd_t
    enddo
  endif
  if (reaction%ncollcomp > 0) then
    ! dRj_dCj - mobile
    J(istart:iendaq,istart:iendaq) = J(istart:iendaq,istart:iendaq) + &
      rt_auxvar%colloid%dRj_dCj%dtotal(:,:,1)*psvd_t
    ! need the below
    ! dRj_dSic
    ! dRic_dCj                                 
  endif
  if (reaction%gas%nactive_gas > 0) then
    iphase = 2
    ! units of dtotal(:,:,2) = kg water / L gas
    psvd_t = material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0* &
             material_auxvar%volume/option%tran_dt
    J(istart:iendaq,istart:iendaq) = J(istart:iendaq,istart:iendaq) + &
                                     rt_auxvar%aqueous%dtotal(:,:,iphase) * &
                                     psvd_t
  endif

end subroutine RTAccumulationDerivative

! ************************************************************************** !

subroutine RCalculateCompression(global_auxvar,rt_auxvar,material_auxvar, &
                                 reaction,option)
  ! 
  ! Calculates the compression for the Jacobian block
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/12/10
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module

  implicit none
  
  type(reaction_type), pointer :: reaction
  type(option_type) :: option

  PetscInt :: dfill(reaction%ncomp,reaction%ncomp)
  PetscInt :: ofill(reaction%ncomp,reaction%ncomp)
  PetscReal :: J(reaction%ncomp,reaction%ncomp)
  PetscReal :: residual(reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  PetscInt :: i, jj
  PetscReal :: vol = 1.d0
  PetscReal :: por = 0.25d0
  PetscReal :: sum
  
  dfill = 0
  ofill = 0
  J = 0.d0
  residual = 0.d0

  call RTAuxVarCompute(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
  call RTAccumulationDerivative(rt_auxvar,global_auxvar, &
                                material_auxvar,reaction,option,J)
    
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (dabs(J(i,jj)) > 1.d-20) ofill(i,jj) = 1
    enddo
  enddo

  if (reaction%neqsorb > 0) then
    call RAccumulationSorbDerivative(rt_auxvar,global_auxvar, &
                                     material_auxvar, &
                                     reaction,option,J)
  endif

  call RReaction(residual,J,PETSC_TRUE,rt_auxvar,global_auxvar, &
                 material_auxvar,reaction,option)
 
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (dabs(J(i,jj)) > 1.d-20) dfill(i,jj) = 1
    enddo
  enddo

  sum = 0.d0
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (dfill(i,jj) == 1) sum = sum + 1.d0
    enddo
  enddo
  write(option%io_buffer,'(''Diagonal Fill (%): '',f6.2)') &
    sum / (reaction%ncomp*reaction%ncomp) * 100.d0
  call printMsg(option)
 
 
  sum = 0.d0
  do jj = 1, reaction%ncomp
    do i = 1, reaction%ncomp
      if (ofill(i,jj) == 1) sum = sum + 1.d0
    enddo
  enddo
  write(option%io_buffer,'(''Off-Diagonal Fill (%): '',f6.2)') &
    sum / (reaction%ncomp*reaction%ncomp) * 100.d0
  call printMsg(option)

end subroutine RCalculateCompression

! ************************************************************************** !

subroutine RUpdateKineticState(rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option)
  ! 
  ! Updates state variables such as mineral vol frac,
  ! etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/24/13
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar  
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: imnrl, iaqspec, ncomp, icomp
  PetscInt :: k, irate, irxn, icplx, ncplx, ikinrxn
  PetscReal :: kdt, one_plus_kdt, k_over_one_plus_kdt
  PetscReal :: delta_volfrac
  PetscReal :: res(reaction%ncomp)
  PetscReal :: jac(reaction%ncomp,reaction%ncomp)
    
  ! update mineral volume fractions
  if (reaction%mineral%nkinmnrl > 0) then
  
    ! Updates the mineral rates, res is not needed
    call RKineticMineral(res,jac,PETSC_FALSE,rt_auxvar,global_auxvar, &
                         material_auxvar,reaction,option)  
                    
    do imnrl = 1, reaction%mineral%nkinmnrl
      ! rate = mol/m^3/sec
      ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
      !                                mol_vol (m^3 mnrl/mol mnrl)
      delta_volfrac = rt_auxvar%mnrl_rate(imnrl)* &
                      reaction%mineral%kinmnrl_molar_vol(imnrl)* &
                      option%tran_dt
      rt_auxvar%mnrl_volfrac(imnrl) = rt_auxvar%mnrl_volfrac(imnrl) + &
                                      delta_volfrac
      if (rt_auxvar%mnrl_volfrac(imnrl) < 0.d0) &
        rt_auxvar%mnrl_volfrac(imnrl) = 0.d0

    enddo
  endif

  if (associated(rxn_sandbox_list)) then
    call RSandboxUpdateKineticState(rt_auxvar,global_auxvar, &
                                    material_auxvar,reaction,option)
  endif  

end subroutine RUpdateKineticState

! ************************************************************************** !

subroutine RUpdateTempDependentCoefs(global_auxvar,reaction, &
                                     update_mnrl,option)
  ! 
  ! Updates temperature dependent coefficients for
  ! anisothermal simulations
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/25/13
  ! 

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: global_auxvar  
  type(reaction_type) :: reaction
  PetscBool :: update_mnrl
  type(option_type) :: option
  
  PetscReal :: temp
  PetscReal :: pres
  
  PetscInt, parameter :: iphase = 1
  
  if (.not.reaction%use_geothermal_hpt)then
    temp = global_auxvar%temp
    pres = 0.d0
    if (associated(reaction%eqcplx_logKcoef)) then
      call ReactionInterpolateLogK(reaction%eqcplx_logKcoef, &
                                    reaction%eqcplx_logK, &
                                    temp, &
                                    reaction%neqcplx)
    endif
    if (associated(reaction%gas%acteqlogKcoef)) then
      call ReactionInterpolateLogK(reaction%gas%acteqlogKcoef, &
                                    reaction%gas%acteqlogK, &
                                    temp, &
                                    reaction%gas%npassive_gas)
    endif
    if (associated(reaction%gas%paseqlogKcoef)) then
      call ReactionInterpolateLogK(reaction%gas%paseqlogKcoef, &
                                    reaction%gas%paseqlogK, &
                                    temp, &
                                    reaction%gas%npassive_gas)
    endif
    call MineralUpdateTempDepCoefs(temp,pres,reaction%mineral, &
                                   reaction%use_geothermal_hpt, &
                                   update_mnrl, &
                                   option)
  else ! high pressure and temperature
    temp = global_auxvar%temp
    pres = global_auxvar%pres(iphase)
    if (associated(reaction%eqcplx_logKcoef)) then
      call ReactionInterpolateLogK_hpt(reaction%eqcplx_logKcoef, &
                                       reaction%eqcplx_logK, &
                                       temp, &
                                       pres, &
                                       reaction%neqcplx)
    endif
    if (associated(reaction%gas%acteqlogKcoef)) then
      call ReactionInterpolateLogK_hpt(reaction%gas%acteqlogKcoef, &
                                       reaction%gas%acteqlogK, &
                                       temp, &
                                       pres, &
                                       reaction%gas%npassive_gas)
    endif   
    if (associated(reaction%gas%paseqlogKcoef)) then
      call ReactionInterpolateLogK_hpt(reaction%gas%paseqlogKcoef, &
                                       reaction%gas%paseqlogK, &
                                       temp, &
                                       pres, &
                                       reaction%gas%npassive_gas)
    endif   
    call MineralUpdateTempDepCoefs(temp,pres,reaction%mineral, &
                                   reaction%use_geothermal_hpt, &
                                   update_mnrl, &
                                   option)    
  endif 
  
end subroutine RUpdateTempDependentCoefs

! ************************************************************************** !

subroutine RTPrintAuxVar(rt_auxvar,reaction,option)
  ! 
  ! PrintRTAuxVar: Prints data from RTAuxVar object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/18/2011
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  
  10 format(a20,':',10es19.11)
  20 format(a20,':',a20)
  30 format(/)

  if (OptionPrintToScreen(option)) write(*,30)
  if (OptionPrintToFile(option)) write(option%fid_out,30)

  if (OptionPrintToScreen(option)) &
    write(*,20) 'Primary', 'free molal., total molar., act. coef.'
  if (OptionPrintToFile(option)) &
    write(option%fid_out,20) 'Primary', 'free molal., total molar., act. coef.'
  do i = 1, reaction%naqcomp  
    if (OptionPrintToScreen(option)) &
      write(*,10) reaction%primary_species_names(i), &
        rt_auxvar%pri_molal(i), &
        rt_auxvar%total(i,1), &
        rt_auxvar%pri_act_coef(i)
    if (OptionPrintToFile(option)) &
      write(option%fid_out,10) reaction%primary_species_names(i), &
        rt_auxvar%pri_molal(i), &
        rt_auxvar%total(i,1), &
        rt_auxvar%pri_act_coef(i)
  enddo
  if (OptionPrintToScreen(option)) write(*,30)
  if (OptionPrintToFile(option)) write(option%fid_out,30)

  if (reaction%neqcplx > 0) then
    if (OptionPrintToScreen(option)) &
      write(*,20) 'Secondary Complex', 'molal., act. coef.'
    if (OptionPrintToFile(option)) &
      write(option%fid_out,20) 'Secondary Complex', 'molal., act. coef.'
    do i = 1, reaction%neqcplx  
      if (OptionPrintToScreen(option)) &
        write(*,10) reaction%secondary_species_names(i), &
          rt_auxvar%sec_molal(i), &
          rt_auxvar%sec_act_coef(i)
      if (OptionPrintToFile(option)) &
        write(option%fid_out,10) reaction%secondary_species_names(i), &
          rt_auxvar%sec_molal(i), &
          rt_auxvar%sec_act_coef(i)
    enddo
    if (OptionPrintToScreen(option)) write(*,30)
    if (OptionPrintToFile(option)) write(option%fid_out,30)
  endif

  if (reaction%neqsorb > 0) then  
    if (OptionPrintToScreen(option)) &
      write(*,20) 'Total Sorbed EQ', 'mol/m^3'
    if (OptionPrintToFile(option)) &
      write(option%fid_out,20) 'Total Sorbed EQ', 'mol/m^3'
    do i = 1, reaction%naqcomp  
      if (OptionPrintToScreen(option)) &
        write(*,10) reaction%primary_species_names(i), rt_auxvar%total_sorb_eq(i)
      if (OptionPrintToFile(option)) &
        write(option%fid_out,10) reaction%primary_species_names(i), &
          rt_auxvar%total_sorb_eq(i)
    enddo
    if (OptionPrintToScreen(option)) write(*,30)
    if (OptionPrintToFile(option)) write(option%fid_out,30)
  endif    


  if (reaction%mineral%nkinmnrl > 0) then
    if (OptionPrintToScreen(option)) &
      write(*,20) 'Kinetic Minerals', 'vol frac, area, rate'
    if (OptionPrintToFile(option)) &
      write(option%fid_out,20) 'Kinetic Minerals', 'vol frac, area, rate'
    do i = 1, reaction%mineral%nkinmnrl
      if (OptionPrintToScreen(option)) &
        write(*,10) reaction%mineral%kinmnrl_names(i), &
          rt_auxvar%mnrl_volfrac(i), &
          rt_auxvar%mnrl_area(i), &
          rt_auxvar%mnrl_rate(i)

      if (OptionPrintToFile(option)) &
        write(option%fid_out,10) reaction%mineral%kinmnrl_names(i), &
          rt_auxvar%mnrl_volfrac(i), &
          rt_auxvar%mnrl_area(i), &
          rt_auxvar%mnrl_rate(i)
    enddo
    if (OptionPrintToScreen(option)) write(*,30)
    if (OptionPrintToFile(option)) write(option%fid_out,30)
  endif

end subroutine RTPrintAuxVar

! ************************************************************************** !

subroutine RTSetPlotVariables(list,reaction,option,time_unit)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
  use Option_module
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(output_variable_list_type), pointer :: list
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: time_unit
  
  character(len=MAXWORDLENGTH) :: name,  units
  character(len=MAXSTRINGLENGTH) string
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscInt :: i
  
  if (reaction%print_free_conc_type == PRIMARY_MOLALITY) then
    free_mol_char = 'm'
  else
    free_mol_char = 'M'
  endif
  
  if (reaction%print_tot_conc_type == TOTAL_MOLALITY) then
    tot_mol_char = 'm'
  else
    tot_mol_char = 'M'
  endif
  
  if (reaction%print_secondary_conc_type == SECONDARY_MOLALITY) then
    sec_mol_char = 'm'
  else
    sec_mol_char = 'M'
  endif
  
  if (reaction%print_pH .and. associated(reaction%species_idx)) then
    if (reaction%species_idx%h_ion_id /= 0) then
      name = 'pH'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units,PH, &
                                   reaction%species_idx%h_ion_id)
    else
      option%io_buffer = 'pH may not be printed when H+ is not ' // &
        'defined as a species.'
      call printErrMsg(option)
    endif
  endif  
  
  if (reaction%print_EH .and. associated(reaction%species_idx)) then
    if (reaction%species_idx%o2_gas_id > 0) then
      name = 'Eh'
      units = 'V'
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units,EH, &
                                   reaction%species_idx%h_ion_id)
    else
      option%io_buffer = 'Eh may not be printed when O2(g) is not ' // &
        'defined as a species.'
      call printErrMsg(option)
    endif
  endif  
  
  if (reaction%print_pe .and. associated(reaction%species_idx)) then
    if (reaction%species_idx%o2_gas_id > 0) then
      name = 'pe'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units,PE, &
                                   reaction%species_idx%h_ion_id)
    else
      option%io_buffer = 'pe may not be printed when O2(g) is not ' // &
        'defined as a species.'
      call printErrMsg(option)   
    endif
  endif  
  
  if (reaction%print_O2 .and. associated(reaction%species_idx)) then
    if (reaction%species_idx%o2_gas_id > 0) then
      name = 'logfO2'
      units = 'bars'
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units,O2, &
                                   reaction%species_idx%o2_gas_id)
    else
      option%io_buffer = 'logfO2 may not be printed when O2(g) is not ' // &
        'defined as a species.'
      call printErrMsg(option)
    endif
  endif  
  
  if (reaction%print_total_component) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Total ' // trim(reaction%primary_species_names(i))
        units = trim(tot_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     reaction%print_tot_conc_type,i)
      endif
    enddo
  endif  
  
  if (reaction%print_free_ion) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Free ' // trim(reaction%primary_species_names(i)) 
        units = trim(free_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                      reaction%print_free_conc_type,i)
      endif
    enddo
  endif  
  
#if 0
!geh: not yet supported as we have no function set up to calculate the
!     passive gas partial pressures at this point.  They are calculated
!     in the CONSTRAINT output.
  do i=1,reaction%gas%npassive_gas
    if (reaction%gas%passive_print_me(i)) then
      name = 'Passive Gas ' // trim(reaction%gas%passive_names(i))
      units = trim('Bar')
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
          GAS_CONCENTRATION,i)
    endif
  enddo
#endif

  do i=1,reaction%gas%nactive_gas
    if (reaction%gas%active_print_me(i)) then
      name = 'Active Gas ' // trim(reaction%gas%active_names(i))
      units = trim('Bar')
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
          GAS_CONCENTRATION,i)
    endif
  enddo

  if (reaction%print_total_bulk) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Total Bulk ' // trim(reaction%primary_species_names(i))
        units = 'mol/m^3 bulk'
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     TOTAL_BULK,i)   
      endif
    enddo
  endif
  
  if (reaction%print_act_coefs) then
    do i=1,reaction%naqcomp
      if (reaction%primary_species_print(i)) then
        name = 'Gamma ' // trim(reaction%primary_species_names(i))
        units = ''
        call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                     PRIMARY_ACTIVITY_COEF,i) 
      endif
    enddo
  endif
  
  do i=1,reaction%neqcplx
    if (reaction%secondary_species_print(i)) then
      name = trim(reaction%secondary_species_names(i))
      units = trim(sec_mol_char)
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   reaction%print_secondary_conc_type,i) 
    endif
  enddo   

  do i=1,reaction%mineral%nkinmnrl
    if (reaction%mineral%kinmnrl_print(i)) then
      name = trim(reaction%mineral%kinmnrl_names(i)) // ' VF'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_VOLUME_FRACTION,units, &
                                   MINERAL_VOLUME_FRACTION,i)     
    endif
  enddo
  
  do i=1,reaction%mineral%nkinmnrl
    if (reaction%mineral%kinmnrl_print(i)) then
      name = trim(reaction%mineral%kinmnrl_names(i)) // ' Rate'
      units = 'mol/m^3/sec'
      call OutputVariableAddToList(list,name,OUTPUT_RATE,units, &
                                   MINERAL_RATE,i)      
    endif
  enddo  
  
  if (reaction%mineral%print_saturation_index) then
    do i=1,reaction%mineral%nmnrl
      if (reaction%mineral%mnrl_print(i)) then
        name = trim(reaction%mineral%mineral_names(i)) // ' SI'
        units = ''
        call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                     MINERAL_SATURATION_INDEX,i)    
      endif
    enddo
  endif
  
  do i=1,reaction%immobile%nimmobile
    if (reaction%immobile%print_me(i)) then
      name = trim(reaction%immobile%names(i)) 
      units = 'mol/m^3'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   IMMOBILE_SPECIES,i)
    endif
  enddo
  

  if (associated(reaction%kd_print)) then
    do i=1,reaction%naqcomp
      if (reaction%kd_print(i)) then
      name = trim(reaction%primary_species_names(i)) // ' KD'
      units = '-'
      call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                   PRIMARY_KD,i)
      endif
    enddo
  endif
  
  if (associated(reaction%total_sorb_print)) then
    do i=1,reaction%naqcomp
      if (reaction%total_sorb_print(i)) then
        name = 'Total Sorbed ' // trim(reaction%primary_species_names(i))
        units = 'mol/m^3'
        call  OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                      TOTAL_SORBED,i)        
      endif
    enddo
  endif
  
  if (associated(reaction%total_sorb_mobile_print)) then
    do i=1,reaction%ncollcomp
      if (reaction%total_sorb_mobile_print(i)) then
        name = 'Total Sorbed Mobile ' // &
               trim(reaction%colloid_species_names(i))
        units = trim(tot_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     TOTAL_SORBED_MOBILE,i)  
      endif
    enddo
  endif  
  
  if (reaction%print_colloid) then
    do i=1,reaction%ncoll
      if (reaction%colloid_print(i)) then
        name = 'Mobile Colloidal ' // trim(reaction%colloid_names(i)) 
        units = trim(tot_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     COLLOID_MOBILE,i)         
      endif
    enddo
    do i=1,reaction%ncoll
      if (reaction%colloid_print(i)) then
        name = 'Mobile Colloidal ' // trim(reaction%colloid_names(i))
        units = trim(tot_mol_char)
        call OutputVariableAddToList(list,name,OUTPUT_CONCENTRATION,units, &
                                     COLLOID_IMMOBILE,i)         
      endif
    enddo
  endif

  if (reaction%print_age) then
    if (reaction%species_idx%tracer_age_id > 0) then
      name = 'Tracer Age'
      units = trim(time_unit) // '-molar'
      call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                   AGE,reaction%species_idx%tracer_age_id, &
                                   reaction%species_idx%tracer_aq_id)       
    endif
  endif  

  if (reaction%print_auxiliary) then
    call RSandboxAuxiliaryPlotVariables(list,reaction,option)
  endif
  
end subroutine RTSetPlotVariables

end module Reaction_module
