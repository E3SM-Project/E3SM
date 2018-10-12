module EOS_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use EOS_Water_module
  use EOS_Gas_module
  use EOS_Oil_module
  use EOSData_module

  implicit none

  private


  public :: EOSInit, &
            EOSRead, &
            EOSReferenceDensity, &
            EOSProcess, &
            EOSInputRecord, &
            AllEOSDBaseDestroy

contains

! ************************************************************************** !

subroutine EOSInit()

  implicit none

  call EOSWaterInit()
  call EOSGasInit()
  call EOSOilInit()

  call EOSTableInitList()

end subroutine EOSInit

! ************************************************************************** !

subroutine EOSRead(input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word, subkeyword
  character(len=MAXWORDLENGTH) :: test_filename
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: tempreal, tempreal2
  PetscReal :: rks_tc = UNINITIALIZED_DOUBLE
  PetscReal :: rks_pc = UNINITIALIZED_DOUBLE
  PetscReal :: rks_acen = UNINITIALIZED_DOUBLE
  PetscReal :: rks_omegaa = UNINITIALIZED_DOUBLE
  PetscReal :: rks_omegab = UNINITIALIZED_DOUBLE
  PetscBool :: rks_hydrogen = PETSC_TRUE
  PetscBool :: rks_use_effective_properties = PETSC_TRUE
  PetscBool :: rks_use_cubic_root_solution = PETSC_FALSE
  PetscReal :: temparray(10)
  PetscReal :: test_t_high, test_t_low, test_p_high, test_p_low
  PetscInt :: test_n_temp, test_n_pres
  PetscBool :: test_uniform_temp, test_uniform_pres
  PetscErrorCode :: ierr

  input%ierr = 0

  call InputReadWord(input,option,keyword,PETSC_TRUE)
  call InputErrorMsg(input,option,'keyword','EOS')
  call StringToUpper(keyword)

  select case(trim(keyword))
    case('WATER')
      do
        temparray = 0.d0
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','EOS,WATER')
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('DENSITY')
            temparray = 0.d0
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DENSITY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                               'kg/m^3','EOS,WATER,DENSITY,CONSTANT',option)
              case('EXPONENTIAL','BRAGFLO')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'REFERENCE_DENSITY', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
                call InputReadDouble(input,option,temparray(2))
                call InputErrorMsg(input,option,'REFERENCE_PRESSURE', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
                call InputReadDouble(input,option,temparray(3))
                call InputErrorMsg(input,option,'WATER_COMPRESSIBILITY', &
                                   'EOS,WATER,DENSITY,EXPONENTIAL')
              case('LINEAR')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'REFERENCE_DENSITY', &
                                   'EOS,WATER,DENSITY,LINEAR')
                call InputReadDouble(input,option,temparray(2))
                call InputErrorMsg(input,option,'REFERENCE_PRESSURE', &
                                   'EOS,WATER,DENSITY,LINEAR')
                call InputReadDouble(input,option,temparray(3))
                call InputErrorMsg(input,option,'WATER_COMPRESSIBILITY', &
                                   'EOS,WATER,DENSITY,LINEAR')
              case('QUADRATIC')
                do
                  call InputReadPflotranString(input,option)
                  call InputReadStringErrorMsg(input,option, &
                                               'EOS,WATER,DENSITY,QUADRATIC')
                  if (InputCheckExit(input,option)) exit
                  if (InputError(input)) exit
                  call InputReadWord(input,option,subkeyword,PETSC_TRUE)
                  call InputErrorMsg(input,option,'subkeyword', &
                                       'EOS,WATER,DENSITY,QUADRATIC')
                  select case(trim(subkeyword))
                    case('REFERENCE_DENSITY')
                      call InputReadDouble(input,option,temparray(1))
                      call InputErrorMsg(input,option,'REFERENCE_DENSITY', &
                                         'EOS,WATER,DENSITY,QUADRATIC')
                    case('REFERENCE_PRESSURE')
                      call InputReadDouble(input,option,temparray(2))
                      call InputErrorMsg(input,option,'REFERENCE_PRESSURE', &
                                         'EOS,WATER,DENSITY,QUADRATIC')
                    case('WATER_COMPRESSIBILITY')
                      call InputReadDouble(input,option,temparray(3))
                      call InputErrorMsg(input,option,'WATER_COMPRESSIBILITY', &
                                         'EOS,WATER,DENSITY,QUADRATIC')
                    case default
                      call InputKeywordUnrecognized(subkeyword, &
                                'EOS,WATER,DENSITY,QUADRATIC',option)
                  end select
                enddo
              case('IFC67','DEFAULT','BATZLE_AND_WANG','TGDPB01','PLANAR', &
                   'TRANGENSTEIN')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,DENSITY',option)
            end select
            call EOSWaterSetDensity(word,temparray)
          case('ENTHALPY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'ENTHALPY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,ENTHALPY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                               'J/kmol','EOS,WATER,ENTHALPY,CONSTANT',option)
               case('IFC67','PAINTER','DEFAULT','PLANAR')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,ENTHALPY',option)
            end select
            call EOSWaterSetEnthalpy(word,temparray)
          case('VISCOSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'VISCOSITY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,VISCOSITY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                              'Pa-s','EOS,WATER,VISCOSITY,CONSTANT',option)
              case('DEFAULT','BATZLE_AND_WANG','GRABOWSKI')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,VISCOSITY', &
                                              option)
            end select
            call EOSWaterSetViscosity(word,temparray)
          case('STEAM_DENSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'STEAM_DENSITY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,STEAM_DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                           'kg/m^3','EOS,WATER,STEAM_DENSITY,CONSTANT',option)
              case('IFC67','DEFAULT','PLANAR')
              case default
                call InputKeywordUnrecognized(word,'EOS,WATER,STEAM_DENSITY', &
                                              option)
            end select
            call EOSWaterSetSteamDensity(keyword,temparray)
          case('STEAM_ENTHALPY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'STEAM_ENTHALPY','EOS,WATER')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,temparray(1))
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,WATER,STEAM_ENTHALPY,CONSTANT')
                call InputReadAndConvertUnits(input,temparray(1), &
                        'J/kmol','EOS,WATER,STEAM_ENTHALPY,CONSTANT',option)
              case('IFC67','DEFAULT','PLANAR')
              case default
                call InputKeywordUnrecognized(word, &
                       'EOS,WATER,STEAM_ENTHALPY',option)
            end select
            call EOSWaterSetSteamEnthalpy(keyword,temparray)
          case('TEST')
            if (option%global_rank == 0) then
              call InputReadDouble(input,option,test_t_low)
              call InputErrorMsg(input,option,'T_low', &
                                 'EOS,WATER,TEST')
              call InputReadDouble(input,option,test_t_high)
              call InputErrorMsg(input,option,'T_high', &
                                 'EOS,WATER,TEST')
              call InputReadDouble(input,option,test_p_low)
              call InputErrorMsg(input,option,'P_low', &
                                 'EOS,WATER,TEST')
              call InputReadDouble(input,option,test_p_high)
              call InputErrorMsg(input,option,'P_high', &
                                 'EOS,WATER,TEST')
              call InputReadInt(input,option,test_n_temp)
              call InputErrorMsg(input,option,'num_temperatures', &
                                 'EOS,WATER,TEST')
              call InputReadInt(input,option,test_n_pres)
              call InputErrorMsg(input,option,'num_pressures', &
                                 'EOS,WATER,TEST')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'temperature distribution type', &
                                 'EOS,WATER,TEST')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_temp = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_temp = PETSC_FALSE
              else
                option%io_buffer = 'Temperature distribution type "' // &
                  trim(word) // '" for EOS Water not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'pressure distribution type', &
                                 'EOS,WATER,TEST,')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_pres = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_pres = PETSC_FALSE
              else
                option%io_buffer = 'Pressure distribution type "' // &
                  trim(word) // '" for EOS Water not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              test_filename = ''
              if (input%ierr == 0) then
                test_filename = word
              endif
              call EOSWaterTest(test_t_low,test_t_high,test_p_low,test_p_high, &
                                test_n_temp, test_n_pres, &
                                test_uniform_temp, test_uniform_pres, &
                                test_filename)
            endif
          case default
            call InputKeywordUnrecognized(keyword,'EOS,WATER',option)
        end select
      enddo
      string = ''
      call EOSWaterVerify(ierr,string)
      if (ierr /= 0) then
        option%io_buffer = 'Error in Water EOS'
        if (len_trim(string) > 1) then
          option%io_buffer = trim(option%io_buffer) // ': ' // trim(string)
        endif
        call printErrMsg(option)
      endif
    case('GAS')
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','EOS,GAS')
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('REFERENCE_DENSITY','SURFACE_DENSITY','STANDARD_DENSITY')
            call InputReadDouble(input,option,tempreal)
            call InputErrorMsg(input,option,'VALUE', &
                               'EOS,GAS,REFERENCE_DENSITY')
            call InputReadAndConvertUnits(input,tempreal, &
                           'kg/m^3','EOS,GAS,REFERENCE_DENSITY',option)
            call EOSGasSetReferenceDensity(tempreal)
          case('DENSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DENSITY','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                            'kmol/m^3','EOS,GAS,DENSITY,CONSTANT',option)
                call EOSGasSetDensityConstant(tempreal)
              case('RKS')
                ! if nothing is entered, it will calculate as hydrogen gas
                  do
                    call InputReadPflotranString(input,option)
                    call InputReadStringErrorMsg(input,option, &
                                                 'EOS GAS,RKS')
                    if (InputCheckExit(input,option)) exit
                    if (InputError(input)) exit
                    call InputReadWord(input,option,word,PETSC_TRUE)
                    call InputErrorMsg(input,option,'keyword', &
                                       'EOS GAS, RKS')
                    select case(trim(word))
                      case('HYDROGEN')
                        rks_hydrogen = PETSC_TRUE
                      case('NON-HYDROGEN')
                        rks_hydrogen = PETSC_FALSE
                      case('USE_EFFECTIVE_PROPERTIES')
                        rks_use_effective_properties = PETSC_TRUE
                      case('DONT_USE_EFFECTIVE_PROPERTIES')
                        rks_use_effective_properties = PETSC_FALSE
                      case('USE_CUBIC_ROOT_SOLUTION')
                        rks_use_cubic_root_solution = PETSC_TRUE
                      case('DONT_USE_CUBIC_ROOT_SOLUTION')
                        rks_use_cubic_root_solution = PETSC_FALSE
                      case('CRITICAL_TEMPERATURE','TC')
                        call InputReadDouble(input,option,rks_tc)
                        call InputErrorMsg(input,option, &
                                            'critical temperature for RKS', &
                                            'EOS GAS,RKS')
                      case('CRITICAL_PRESSURE','PC')
                        call InputReadDouble(input,option,rks_pc)
                        call InputErrorMsg(input,option, &
                                            'critical pressure for RKS', &
                                            'EOS GAS,RKS')
                      case('ACENTRIC,ACENTRIC_FACTOR','ACEN','AC')
                        ! acentric factor is only used for non-hydrogen gas
                        call InputReadDouble(input,option,rks_acen)
                        call InputErrorMsg(input,option, &
                                            'accentric factor for RKS', &
                                            'EOS GAS,RKS')
                      case('OMEGAA','A')
                        call InputReadDouble(input,option,rks_omegaa)
                        call InputErrorMsg(input,option, &
                                        'omega_a factor for RKS', &
                                            'EOS GAS,RKS')
                      case('OMEGAB','B')
                        call InputReadDouble(input,option,rks_omegab)
                        call InputErrorMsg(input,option, &
                                        'omega_b factor for RKS', &
                                            'EOS GAS,RKS')
                      case default
                        call InputKeywordUnrecognized(word, &
                                'EOS GAS,RKS',option)
                    end select
                enddo
                call EOSGasSetDensityRKS(rks_hydrogen, &
                                         rks_use_effective_properties, &
                                         rks_use_cubic_root_solution, &
                                         rks_tc,rks_pc,rks_acen, &
                                         rks_omegaa,rks_omegab)
              case('PR_METHANE')
                call EOSGasSetDensityPRMethane()
              case('IDEAL','DEFAULT')
                call EOSGasSetDensityIdeal()
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,DENSITY',option)
            end select
          case('ENTHALPY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'ENTHALPY','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,ENTHALPY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                                 'J/kmol','EOS,GAS,ENTHALPY,CONSTANT',option)
                call EOSGasSetEnergyConstant(tempreal)
              case('IDEAL_METHANE')
                call EOSGasSetEnergyIdealMethane()
              case('IDEAL','DEFAULT')
                call EOSGasSetEnergyIdeal()
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,ENTHALPY',option)
            end select
          case('VISCOSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'VISCOSITY','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,VISCOSITY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                                 'Pa-s','EOS,GAS,VISCOSITY,CONSTANT',option)
                call EOSGasSetViscosityConstant(tempreal)
              case('DEFAULT')
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,VISCOSITY',option)
            end select
          case('HENRYS_CONSTANT')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'HENRYS_CONSTANT','EOS,GAS')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,GAS,HENRYS_CONSTANT,CONSTANT')
                call EOSGasSetHenryConstant(tempreal)
              case('DEFAULT')
                call EOSGasSetHenry()
              case default
                call InputKeywordUnrecognized(word,'EOS,GAS,HENRYS_CONSTANT', &
                                              option)
            end select
          case('TEST')
            if (option%global_rank == 0) then
              call InputReadDouble(input,option,test_t_low)
              call InputErrorMsg(input,option,'T_low', &
                                 'EOS,GAS,TEST')
              call InputReadDouble(input,option,test_t_high)
              call InputErrorMsg(input,option,'T_high', &
                                 'EOS,GAS,TEST')
              call InputReadDouble(input,option,test_p_low)
              call InputErrorMsg(input,option,'P_low', &
                                 'EOS,GAS,TEST')
              call InputReadDouble(input,option,test_p_high)
              call InputErrorMsg(input,option,'P_high', &
                                 'EOS,GAS,TEST')
              call InputReadInt(input,option,test_n_temp)
              call InputErrorMsg(input,option,'num_temperatures', &
                                 'EOS,GAS,TEST')
              call InputReadInt(input,option,test_n_pres)
              call InputErrorMsg(input,option,'num_pressures', &
                                 'EOS,GAS,TEST')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'temperature distribution type', &
                                 'EOS,GAS,TEST')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_temp = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_temp = PETSC_FALSE
              else
                option%io_buffer = 'Temperature distribution type "' // &
                  trim(word) // '" for EOS Gas not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'pressure distribution type', &
                                 'EOS,GAS,TEST,')
              if (StringCompareIgnoreCase(word,'uniform')) then
                test_uniform_pres = PETSC_TRUE
              else if (StringCompareIgnoreCase(word,'log')) then
                test_uniform_pres = PETSC_FALSE
              else
                option%io_buffer = 'Pressure distribution type "' // &
                  trim(word) // '" for EOS Gas not recognized.'
                call printErrMsg(option)
              endif
              call InputReadWord(input,option,word,PETSC_TRUE)
              test_filename = ''
              if (input%ierr == 0) then
                test_filename = word
              endif
              call EOSGasTest(test_t_low,test_t_high,test_p_low,test_p_high, &
                              test_n_temp, test_n_pres, &
                              test_uniform_temp, test_uniform_pres, &
                              test_filename)
            endif
          case('FORMULA_WEIGHT')
            call InputReadDouble(input,option,tempreal)
            call InputErrorMsg(input,option,'VALUE','EOS,GAS,FORMULA_WEIGHT')
            call InputReadAndConvertUnits(input,tempreal, &
                             'g/mol','EOS,GAS,FORMULA_WEIGHT',option)
            call EOSGasSetFMWConstant(tempreal)
          case('CO2_SPAN_WAGNER_DB')
            call EOSGasSetFMWConstant(FMWCO2)
            temparray = UNINITIALIZED_DOUBLE
            subkeyword =''
            do
              call InputReadPflotranString(input,option)
              call InputReadStringErrorMsg(input,option, &
                                           'EOS GAS,CO2_SPAN_WAGNER_DB')
              if (InputCheckExit(input,option)) exit
              if (InputError(input)) exit
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword', &
                                       'EOS GAS, CO2_SPANWAGNER_DB')
              select case(trim(word))
                case('PRESSURE_MIN')
                  call InputReadDouble(input,option,temparray(1))
                  call InputErrorMsg(input,option, &
                                    'min pressure - properties look up', &
                                    'EOS GAS,CO2_SPAN_WAGNER_DB')
                  call InputReadAndConvertUnits(input,temparray(1), &
                       'Pa','EOS,GAS,CO2_SPAN_WAGNER_DB,PRESSURE_MIN',option)
                case('PRESSURE_MAX')
                  call InputReadDouble(input,option,temparray(2))
                  call InputErrorMsg(input,option, &
                                    'MAX pressure - properties look up', &
                                    'EOS GAS,CO2_SPAN_WAGNER_DB')
                  call InputReadAndConvertUnits(input,temparray(2), &
                       'Pa','EOS,GAS,CO2_SPAN_WAGNER_DB,PRESSURE_MAX',option)
                case('PRESSURE_DELTA')
                  call InputReadDouble(input,option,temparray(3))
                  call InputErrorMsg(input,option, &
                                    'Delta pressure - properties look up', &
                                    'EOS GAS,CO2_SPAN_WAGNER_DB')
                  call InputReadAndConvertUnits(input,temparray(3), &
                      'Pa','EOS,GAS,CO2_SPAN_WAGNER_DB,PRESSURE_DELTA',option)
                case('TEMPERATURE_MIN')
                  call InputReadDouble(input,option,temparray(4))
                  call InputErrorMsg(input,option, &
                                    'min temperature - properties look up', &
                                    'EOS GAS,CO2_SPAN_WAGNER_DB')
                  call InputReadAndConvertUnits(input,temparray(4), &
                    'C','EOS,GAS,CO2_SPAN_WAGNER_DB,TEMPERATURE_MIN',option)
                case('TEMPERATURE_MAX')
                  call InputReadDouble(input,option,temparray(5))
                  call InputErrorMsg(input,option, &
                                    'MAX temperature - properties look up', &
                                    'EOS GAS,CO2_SPAN_WAGNER_DB')
                  call InputReadAndConvertUnits(input,temparray(5), &
                    'C','EOS,GAS,CO2_SPAN_WAGNER_DB,TEMPERATURE_MAX',option)
                case('TEMPERATURE_DELTA')
                  call InputReadDouble(input,option,temparray(6))
                  call InputErrorMsg(input,option, &
                                    'Delta temperature - properties look up', &
                                    'EOS GAS,CO2_SPAN_WAGNER_DB')
                  call InputReadAndConvertUnits(input,temparray(6), &
                    'C','EOS,GAS,CO2_SPAN_WAGNER_DB,TEMPERATURE_DELTA',option)
                case('DATABASE_FILE_NAME')
                  call InputReadWord(input,option,subkeyword,PETSC_TRUE)
                  call InputErrorMsg(input,option, &
                                     'databas file name',&
                                     'EOS,GAS,FORMULA_WEIGHT')
                case default
                  call InputKeywordUnrecognized(subkeyword,&
                                     'EOS,GAS,CO2_SPAN_WAGNER_DB',option)
              end select
            end do
            call MPI_Barrier(option%mycomm,ierr)
            call EOSGasSetEOSDBase(subkeyword,option)
          case('DATABASE')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'EOS,GAS','DATABASE filename')
            call EOSGasSetEOSDBase(word,option)
          case('PVDG')
            call EOSGasSetPVDG(input,option)
          case('CO2_DATABASE')
            call EOSGasSetFMWConstant(FMWCO2)
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'EOS,GAS','DATABASE filename')
            call EOSGasSetEOSDBase(word,option)
          case default
            call InputKeywordUnrecognized(keyword,'EOS,GAS',option)
        end select
      enddo
      string = ''
      call EOSGasVerify(ierr,string)
      if (ierr == 5) then
        option%io_buffer = 'set to default value for RKS hydrogen'
        if (len_trim(string) > 1) then
          option%io_buffer =  trim(string) // ': ' // trim(option%io_buffer)
        endif
        call printMsg(option)
      !else if (ierr == 6) then
      !  option%io_buffer = 'set as default value for gas fmw'
      !  if (len_trim(string) > 1) then
      !    option%io_buffer =  trim(string) // ': ' // trim(option%io_buffer)
      !  endif
      !  call printMsg(option)
      else if (ierr /= 0) then
        option%io_buffer = 'Error in Gas EOS'
        if (len_trim(string) > 1) then
          option%io_buffer = trim(option%io_buffer) // ': ' // trim(string)
        endif
        call printErrMsg(option)
      endif
    case('OIL')
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','EOS,OIL')
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('PVDO')
            call EOSOilSetPVDO(input,option)
          case('PVCO')
            call EOSOilSetPVCO(input,option)
          case('DATABASE')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'EOS,OIL','DATABASE filename')
            call EOSOilSetEOSDBase(word,option)
          case('VISCOSITY_LINLOG_INTERPOLATION')
            call EOSOilSetVisLinLogInterp(option)
          case('REFERENCE_DENSITY','SURFACE_DENSITY','STANDARD_DENSITY')
            call InputReadDouble(input,option,tempreal)
            call InputErrorMsg(input,option,'VALUE', &
                               'EOS,OIL,REFERENCE_DENSITY')
            call InputReadAndConvertUnits(input,tempreal, &
                             'kg/m^3','EOS,OIL,REFERENCE_DENSITY',option)
            call EOSOilSetReferenceDensity(tempreal)
          case('DENSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DENSITY','EOS,OIL')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,OIL,DENSITY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                                 'kg/m^3','EOS,OIL,DENSITY,CONSTANT',option)
                call EOSOilSetDensityConstant(tempreal)
              case('LINEAR','INVERSE_LINEAR')
                select case(trim(word))
                  case('LINEAR')
                    call EOSOilSetDensityLinear()
                  case('INVERSE_LINEAR')
                    call EOSOilSetDensityInverseLinear()
                end select
                do
                  call InputReadPflotranString(input,option)
                  if (InputCheckExit(input,option)) exit
                  call InputReadWord(input,option,subkeyword,PETSC_TRUE)
                  call InputErrorMsg(input,option,'subkeyword','EOS,OIL,VIS')
                  call StringToUpper(subkeyword)
                  select case(subkeyword)
                    case('REFERENCE_VALUE')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,DENSITY_LINEAR,REFERENCE_VALUE')
                      call EOSOilSetDenLinearRefDen(tempreal)
                    case('PRES_REF_VALUE')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,DENSITY_LINEAR,PRES_REF_VALUE')
                      call EOSOilSetDenLinearRefPres(tempreal)
                    case('TEMP_REF_VALUE')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,DENSITY_LINEAR,TEMP_REF_VALUE')
                      call EOSOilSetDenLinearRefTemp(tempreal)
                    case('COMPRESS_COEFF')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,DENSITY_LINEAR,COMPRESS_COEFF')
                      call EOSOilSetDenLinearComprCoef(tempreal)
                    case('THERMAL_EXPANSION_COEFF')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,DENSITY_LINEAR,THERMAL_EXPANSION_COEFF')
                      call EOSOilSetDenLinearExpanCoef(tempreal)
                    case default
                      call InputKeywordUnrecognized(subkeyword, &
                           'EOS,OIL,DENSITY_LINEAR',option)
                  end select
                end do
              case('DATABASE')
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'EOS,OIL','DEN DBASE filename')
                call EOSOilSetDenDBase(word,option)
              case default
                call InputKeywordUnrecognized(word,'EOS,OIL,DENSITY',option)
            end select
          case('ENTHALPY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'ENTHALPY','EOS,OIL')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,OIL,ENTHALPY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                                  'J/kmol','EOS,OIL,ENTHALPY,CONSTANT',option)
                call EOSOilSetEnthalpyConstant(tempreal)
              case('LINEAR_TEMP')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,OIL,ENTHALPY,LINEAR_TEMP')
                call EOSOilSetEnthalpyLinearTemp(tempreal)
              case('QUADRATIC_TEMP')
                call EOSOilSetEnthalpyQuadraticTemp()
                do
                  call InputReadPflotranString(input,option)
                  if (InputCheckExit(input,option)) exit
                  call InputReadWord(input,option,subkeyword,PETSC_TRUE)
                  call InputErrorMsg(input,option,'subkeyword',&
                                     'EOS,OIL,ENTHALPY')
                  call StringToUpper(subkeyword)
                  select case(subkeyword)
                    case('TEMP_REF_VALUES')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,ENTHALPY_QUAD,TEMP_REF_VALUES_1')
                      call InputReadDouble(input,option,tempreal2)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,ENTHALPY_QUAD,TEMP_REF_VALUES_2')
                      call EOSOilSetEntQuadRefTemp(tempreal,tempreal2)
                    case('TEMP_COEFFICIENTS')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,ENTHALPY_QUAD,TEMP_COEFF_1')
                      call InputReadDouble(input,option,tempreal2)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,ENTHALPY_QUAD,TEMP_COEFF_2')
                      call EOSOilSetEntQuadTempCoef(tempreal,tempreal2)
                  end select
                end do
              case('DATABASE')
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'EOS,OIL','ENT DBASE filename')
                call EOSOilSetEntDBase(word,option)
              case default
                call InputKeywordUnrecognized(word,'EOS,OIL,ENTHALPY',option)
            end select
          case('VISCOSITY')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'VISCOSITY','EOS,OIL')
            call StringToUpper(word)
            select case(trim(word))
              case('CONSTANT')
                call InputReadDouble(input,option,tempreal)
                call InputErrorMsg(input,option,'VALUE', &
                                   'EOS,OIL,VISCOSITY,CONSTANT')
                call InputReadAndConvertUnits(input,tempreal, &
                                 'Pa-s','EOS,OIL,VISCOSITY,CONSTANT',option)
                call EOSOilSetViscosityConstant(tempreal)
              case('QUADRATIC')
                call EOSOilSetViscosityQuad()
                do
                  call InputReadPflotranString(input,option)
                  if (InputCheckExit(input,option)) exit
                  call InputReadWord(input,option,subkeyword,PETSC_TRUE)
                  call InputErrorMsg(input,option,'subkeyword','EOS,OIL,VIS')
                  call StringToUpper(subkeyword)
                  select case(subkeyword)
                    case('REFERENCE_VALUE')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,REFERENCE_VALUE')
                      call EOSOilSetVisQuadRefVis(tempreal)
                    case('PRES_REF_VALUES')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,PRES_REF_VALUES_1')
                      call InputReadDouble(input,option,tempreal2)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,PRES_REF_VALUES_2')
                      call EOSOilSetVisQuadRefPres(tempreal,tempreal2)
                    case('TEMP_REF_VALUES')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,TEMP_REF_VALUES_1')
                      call InputReadDouble(input,option,tempreal2)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,TEMP_REF_VALUES_2')
                      call EOSOilSetVisQuadRefTemp(tempreal,tempreal2)
                    case('PRES_COEFFICIENTS')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,PRES_COEFF_1')
                      call InputReadDouble(input,option,tempreal2)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,PRES_COEFF_2')
                      call EOSOilSetVisQuadPresCoef(tempreal,tempreal2)
                    case('TEMP_COEFFICIENTS')
                      call InputReadDouble(input,option,tempreal)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,TEMP_COEFF_1')
                      call InputReadDouble(input,option,tempreal2)
                      call InputErrorMsg(input,option,'VALUE', &
                            'EOS,OIL,VISCOSITY_QUAD,TEMP_COEFF_2')
                      call EOSOilSetVisQuadTempCoef(tempreal,tempreal2)
                    case default
                      call InputKeywordUnrecognized(subkeyword, &
                           'EOS,OIL, VISCOSITY_QUAD',option)
                  end select
                end do
              case('DATABASE')
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'EOS,OIL','VIS DBASE filename')
                call EOSOilSetVisDBase(word,option)
              case default
                call InputKeywordUnrecognized(word,'EOS,OIL,VISCOSITY',option)
            end select
          case('FORMULA_WEIGHT')
            call InputReadDouble(input,option,tempreal)
            !call InputReadDouble(input,option,fmw_oil)
            call InputErrorMsg(input,option,'VALUE','EOS,OIL,FORMULA_WEIGHT')
            call InputReadAndConvertUnits(input,tempreal, &
                             'g/mol','EOS,OIL,FORMULA_WEIGHT',option)
            call EOSOilSetFMWConstant(tempreal)
          case default
            call InputKeywordUnrecognized(keyword,'EOS,OIL',option)
        end select
      enddo
      ! need to add verifying function - follow EOSgas template
      string = ''
      call EOSOilVerify(ierr,string)
      if (ierr /= 0) then
        option%io_buffer = 'Error in Oil EOS'
        if (len_trim(string) > 1) then
          option%io_buffer = trim(option%io_buffer) // ': ' // trim(string)
        endif
        call printErrMsg(option)
      endif
    case default
      call InputKeywordUnrecognized(keyword,'EOS',option)
  end select

end subroutine EOSRead

! **************************************************************************** !

subroutine EOSReferenceDensity(option)
  ! 
  ! Calculates reference densities for phases if not specified in input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/18
  ! 
  use Option_module
  
  implicit none
  
  type(option_type) :: option

  PetscReal :: vapor_saturation_pressure
  PetscReal :: vapor_density_kg
  PetscReal :: air_density_kmol
  PetscReal :: dum1, dum2
  PetscErrorCode :: ierr
  
  if (option%reference_density(option%liquid_phase) < 1.d-40) then
    call EOSWaterDensity(option%reference_temperature, &
                         option%reference_pressure, &
                         option%reference_density(option%liquid_phase), &
                         dum1,ierr)
  endif
  if (option%reference_density(option%gas_phase) < 1.d-40) then
    ! assume saturated vapor pressure
    call EOSWaterSaturationPressure(option%reference_temperature, &
                                    vapor_saturation_pressure,dum1,ierr)
    call EOSWaterSteamDensityEnthalpy(option%reference_temperature, &
                                      vapor_saturation_pressure, &
                                      vapor_density_kg,dum1,dum2,ierr)
    call EOSGasDensity(option%reference_temperature, &
                       option%reference_pressure-vapor_saturation_pressure, &
                       air_density_kmol,dum1,dum2,ierr)
    option%reference_density(option%gas_phase) = &
      vapor_density_kg + air_density_kmol*FMWAIR
  endif

end subroutine EOSReferenceDensity

! **************************************************************************** !

subroutine EOSProcess(option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/28/17
  !
  ! Process EOS data after all input data have been read

  use Option_module

  implicit none

  type(option_type) :: option

  call EOSOilTableProcess(option,EOSGasGetFMW(),EOSGasGetReferenceDensity())
  call EOSGasTableProcess(option)

  call EOSTableProcessList(option)

end subroutine EOSProcess

! **************************************************************************** !

subroutine EOSInputRecord()
  !
  ! Prints ingested equation of state information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 05/04/2016
  !
  use EOS_Water_module
  use EOS_Gas_module
  use EOS_Oil_module

  implicit none

  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'EQUATIONS OF STATE (EOS)'

  call EOSWaterInputRecord()
  call EOSGasInputRecord()
  call EOSOilInputRecord()

end subroutine EOSInputRecord

! ************************************************************************** !

subroutine AllEOSDBaseDestroy()
  !
  ! Author: Paolo Orsini
  ! Date: 12/14/15
  !

  implicit none

  call EOSTableDestroyList()

  call EOSOilDBaseDestroy()
  call EOSGasDBaseDestroy()

end subroutine AllEOSDBaseDestroy

! ************************************************************************** !

end module EOS_module
