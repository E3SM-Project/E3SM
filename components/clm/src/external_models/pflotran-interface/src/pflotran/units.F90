module Units_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
  
  public :: UnitsConvertToInternal, UnitsConvertToExternal
  
contains

! ************************************************************************** !

function UnitsConvertToInternal(units,internal_units,option)
  ! 
  ! Converts given units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! Notes: Updated/modified by Jenn Frederick 1/25/2016
  ! 

  use Option_module

  implicit none
  
  character(len=*) :: units
  character(len=MAXWORDLENGTH) :: internal_units
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: units_buff
  character(len=MAXSTRINGLENGTH) :: internal_units_buff1
  character(len=MAXSTRINGLENGTH) :: internal_units_buff2
  character(len=MAXSTRINGLENGTH) :: error_msg
  PetscBool :: multi_option, successful, error
  PetscInt :: length, ind_or, num_options
  PetscReal :: conversion_factor, UnitsConvertToInternal

  units_buff = trim(units)
  internal_units_buff1 = trim(internal_units)
  internal_units_buff2 = trim(internal_units)
  multi_option = PETSC_FALSE
  successful = PETSC_FALSE
  conversion_factor = 1.d0
  UnitsConvertToInternal = 1.d0
  num_options = 0
  ind_or = 1
  length = 0

  do while(ind_or /= 0)
    length = len_trim(internal_units_buff1)
    ind_or = index(trim(internal_units_buff1),"|")   
    if (ind_or == 0) then
      call UnitsConvertParse(units_buff,internal_units_buff1, &
                             conversion_factor,error,error_msg)
      if (.not.error) successful = PETSC_TRUE
      if (successful) exit
    else
      multi_option = PETSC_TRUE
      num_options = num_options + 1
      internal_units_buff2 = internal_units_buff1(1:(ind_or-1))
      call UnitsConvertParse(units_buff,internal_units_buff2, &
                             conversion_factor,error,error_msg)
      if (.not.error) successful = PETSC_TRUE
      if (successful) exit
      internal_units_buff1 = internal_units_buff1((ind_or+1):length)
    endif
  enddo
  
  if (.not.successful) then
    option%io_buffer = error_msg
    call printErrMsg(option)
  endif

  UnitsConvertToInternal = conversion_factor

end function UnitsConvertToInternal

! ************************************************************************** !

subroutine UnitsConvertParse(units,internal_units,units_conversion, &
                             error,error_msg)
  ! 
  ! Parses user's and internal units into numerator and denominator 
  ! (if they exist) and gets the conversion factor from user's units
  ! to internal Pflotran units.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/15/2016
  ! 

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: units, internal_units
  PetscReal :: units_conversion
  PetscBool :: error
  character(len=MAXSTRINGLENGTH) :: error_msg
  
  ! A maximum of 3 numerators/denominators are accepted 
  character(len=MAXSTRINGLENGTH) :: numerator_units
  character(len=MAXSTRINGLENGTH) :: denominator_units
  character(len=MAXSTRINGLENGTH) :: numerator_internal_units
  character(len=MAXSTRINGLENGTH) :: denominator_internal_units
  PetscReal :: numerator_conv_factor, denominator_conv_factor
  PetscBool :: set_units_denom, set_internal_units_denom
  PetscInt :: length, ind
  
  error = PETSC_FALSE
  set_units_denom = PETSC_FALSE
  set_internal_units_denom = PETSC_FALSE
  numerator_conv_factor = 1.d0
  denominator_conv_factor = 1.d0
  units_conversion = 1.d0
  numerator_units = ''
  denominator_units = ''
  numerator_internal_units = ''
  denominator_internal_units = ''

  length = len(units)
  ind = index(trim(units),"/")   
  if (ind == 0) then
  ! numerator(s) only
    numerator_units = trim(units)
  else
  ! denominator(s) also exists
    set_units_denom = PETSC_TRUE
    numerator_units = units(1:(ind-1))
    denominator_units = units((ind+1):length)
  endif

  length = len_trim(internal_units)
  ind = index(trim(internal_units),"/")
  if (ind == 0) then
  ! numerator(s) only
    numerator_internal_units = trim(internal_units)
  else
  ! denominator(s) also exists
    set_internal_units_denom = PETSC_TRUE
    numerator_internal_units = internal_units(1:(ind-1))
    denominator_internal_units = internal_units((ind+1):length)
  endif

  ! check if provided units structure matches internal units structure
  if (set_units_denom .neqv. set_internal_units_denom) then
    if ((.not.set_units_denom) .and. (set_internal_units_denom)) then
      error_msg = 'Units provided do not match the expected &
                   &unit structure (numerator/denominator). &
                   &A unit denominator was expected, but was not given: ' &
                   // trim(units)
      error = PETSC_TRUE
    endif
    if ((set_units_denom) .and. (.not.set_internal_units_denom)) then
      error_msg = 'Units provided do not match the expected &
                   &unit structure (numerator/denominator). &
                   &A unit denominator was not expected, but was given: ' &
                   // trim(units)
      error = PETSC_TRUE
    endif
  endif

  if (error) return
  call UnitsConvert(numerator_units,numerator_internal_units, &
                    numerator_conv_factor,error,error_msg)

  if (error) return
  if (set_internal_units_denom) then
    call UnitsConvert(denominator_units,denominator_internal_units, &
                      denominator_conv_factor,error,error_msg)
  endif

  if (error) return
  units_conversion = numerator_conv_factor/denominator_conv_factor

end subroutine UnitsConvertParse

! ************************************************************************** !

subroutine UnitsConvert(units_user,units_internal,units_conversion, &
                        error,error_msg)
  ! 
  ! Converts units given by the user to pflotran internal units
  ! 
  ! Author: Glenn Hammond; Jenn Frederick (updated)
  ! Date: 01/21/09; 2/15/2016 (updated)
  ! 

  use Option_module
  
  implicit none

  character(len=MAXWORDLENGTH) :: units_user
  character(len=MAXWORDLENGTH) :: units_internal
  PetscReal :: units_conversion
  PetscBool :: error
  character(len=MAXSTRINGLENGTH) :: error_msg
   
  ! a maximum of 3 unit categories are allowed
  character(len=MAXWORDLENGTH) :: unit_user(3)
  character(len=MAXWORDLENGTH) :: unit_user_cat(3)
  character(len=MAXSTRINGLENGTH) :: unit_user_buff
  character(len=MAXWORDLENGTH) :: unit_internal(3)
  character(len=MAXWORDLENGTH) :: unit_internal_cat(3)
  character(len=MAXSTRINGLENGTH) :: unit_internal_buff

  PetscInt :: length, ind_dash
  PetscInt :: k
  PetscInt :: num_user_units, num_internal_units
  PetscReal :: conv_user_to_SI, conv_internal_to_SI
  PetscReal :: conversion_user, conversion_internal

  error = PETSC_FALSE
  conv_user_to_SI = 1.d0
  conv_internal_to_SI = 1.d0
  conversion_user = 1.d0
  conversion_internal = 1.d0
  units_conversion = 1.d0 

! ======== SEPARATE OUT GIVEN UNITS ================================= !
  unit_user_buff = trim(units_user)
  unit_user(:) = 'not_assigned'
  unit_user_cat(:) = 'not_assigned'
  num_user_units = 0
  ind_dash = 1
  k = 1
  do while (ind_dash /= 0)
    length = len_trim(unit_user_buff)
    ind_dash = index(trim(unit_user_buff),"-")
    if (ind_dash == 0) then
      unit_user(k) = trim(unit_user_buff)
    else
      unit_user(k) = unit_user_buff(1:(ind_dash-1))
      unit_user_buff = unit_user_buff((ind_dash+1):length)
    endif
    k = k + 1
    if (k > 3) then
      error_msg = 'Maximum number of user units exceeded. &
                  &Unit numerators or denominators are limited to a &
                  &maximum of 3 units each.' 
      error = PETSC_TRUE
    endif
  enddo
  if (error) return
  num_user_units = k - 1

! ======== SEPARATE OUT INTERNAL UNITS ============================== !
  unit_internal_buff = trim(units_internal)
  unit_internal(:) = 'not_assigned'
  unit_internal_cat(:) = 'not_assigned'
  num_internal_units = 0
  k = 1
  ind_dash = 1
  do while (ind_dash /= 0)
    length = len_trim(unit_internal_buff)
    ind_dash = index(trim(unit_internal_buff),"-")
    if (ind_dash == 0) then
      unit_internal(k) = trim(unit_internal_buff)
    else
      unit_internal(k) = unit_internal_buff(1:(ind_dash-1))
      unit_internal_buff = unit_internal_buff((ind_dash+1):length)
    endif
    k = k + 1
    if (k > 3) then
      error_msg = 'Maximum number of internal units exceeded. &
                  &Unit numerators or denominators are limited to a &
                  &maximum of 3 units each.' 
      error = PETSC_TRUE
    endif
  enddo
  if (error) return
  num_internal_units = k - 1


! ======== GET UNIT CATEGORIES OF GIVEN AND INTERNAL UNITS ========== !
  call UnitsCategory(unit_user,unit_user_cat,error,error_msg)
  if (error) return
  call UnitsCategory(unit_internal,unit_internal_cat,error,error_msg)
  if (error) return

! ======== CHECK IF UNIT CATEGORIES ALIGN =========================== !
  call UnitsCategoryCheck(unit_user_cat,unit_internal_cat,error,error_msg)
  if (error) return

! ======== CONVERT INTERNAL UNITS TO SI ============================= !
  k = 1
  do while (k < (num_internal_units + 1))
    call UnitsConvertToSI(unit_internal(k),conversion_internal,error,error_msg)
    if (error) exit
    conv_internal_to_SI = conv_internal_to_SI * conversion_internal
    k = k + 1
  enddo
  if (error) return

! ======== CONVERT USER UNITS TO SI ================================= !
  k = 1
  do while (k < (num_user_units + 1))
    call UnitsConvertToSI(unit_user(k),conversion_user,error,error_msg)
    if (error) exit
    conv_user_to_SI = conv_user_to_SI * conversion_user
    k = k + 1
  enddo
  if (error) return

! ======== CONVERT USER UNITS TO INTERNAL UNITS ===================== !
  units_conversion = conv_user_to_SI / conv_internal_to_SI

end subroutine UnitsConvert

! ************************************************************************** !

subroutine UnitsCategory(unit,unit_category,error,error_msg)
  ! 
  ! Gives the unit category of the given units.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/15/2016
  ! 

  use Option_module
  
  implicit none
 
  ! a maximum of 3 unit categories are allowed
  character(len=MAXWORDLENGTH) :: unit(3)
  character(len=MAXWORDLENGTH) :: unit_category(3)
  PetscBool :: error
  character(len=MAXSTRINGLENGTH) :: error_msg

  PetscInt :: k 

  k = 1
  error = PETSC_FALSE

  do while (k < 4)
    select case(trim(unit(k)))
      case('cm^3','l','L','ml','mL','dm^3','m^3','gal','gallon',&
           'bbl','cf','Mcf')
        unit_category(k) = 'volume'
      case('cm^2','dm^2','m^2','km^2')
        unit_category(k) = 'area'
      case('km','m','met','meter','dm','cm','mm')
        unit_category(k) = 'length'
      case('s','sec','second','min','minute','h','hr','hour','d','day','w', &
           'week','mo','month','y','yr','year')
        unit_category(k) = 'time'
      case('Pa.s','cP','centiPoise','Poise','P')
        unit_category(k) = 'viscosity'
      case('J','kJ','MJ')
        unit_category(k) = 'energy'
      case('W','kW','MW')
        unit_category(k) = 'power'
      case('mol','mole','moles','kmol')
        unit_category(k) = 'molar_mass'
      case('ug','mg','g','kg')
        unit_category(k) = 'mass'
      case('C','Celcius')
        unit_category(k) = 'temperature'
      case('K','Kelvin')  
        unit_category(k) = 'temperature'
        error_msg = 'Kelvin temperature units are not supported. Use Celcius.' 
        error = PETSC_TRUE
      case('Pa','kPa','MPa','Bar','psi')
        unit_category(k) = 'pressure'
      case('M','mM')
        unit_category(k) = 'concentration'
      case('N')
        unit_category(k) = 'force'
      case('unitless','1')
        unit_category(k) = 'unitless'
      case('not_assigned')
        unit_category(k) = 'not_assigned'
      case default
        error_msg = 'The given units are unrecognized: ' // trim(unit(k))
        error = PETSC_TRUE
    end select
    k = k + 1
  enddo
  
end subroutine UnitsCategory

! ************************************************************************** !

subroutine UnitsConvertToSI(unit,conversion_factor,error,error_msg)
  ! 
  ! Converts a unit to SI units.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/21/2016
  ! 

  use Option_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: unit
  PetscReal :: conversion_factor
  PetscBool :: error
  character(len=MAXSTRINGLENGTH) :: error_msg

  conversion_factor = 1.d0
  error = PETSC_FALSE

  select case(trim(unit))
  !---> VOLUME ---> (meter^3)
    case('cm^3','ml','mL')
      conversion_factor = 1.d-6
    case('l','L','dm^3')
      conversion_factor = 1.d-3
    case('m^3')
      conversion_factor = 1.d0
    case('gal','gallon')
      conversion_factor = 3.785411784d-3
    case('cf')
      conversion_factor = 0.02831685
    case('Mcf') !In Field units M means 1000
      conversion_factor = 0.02831685d3
    case('bbl')
      conversion_factor = 0.1589873141
    case('')
  !---> AREA ---> (meter^2)
    case('cm^2')
      conversion_factor = 1.d-4
    case('dm^2')
      conversion_factor = 1.d-2
    case('m^2')
      conversion_factor = 1.d0
    case('km^2')
      conversion_factor = 1.d6
  ! ---> LENGTH ---> (meter)
    case('km')
      conversion_factor = 1000.d0
    case('m','met','meter')
      conversion_factor = 1.d0
    case('dm')
      conversion_factor = 1.d-1
    case('cm')
      conversion_factor = 1.d-2
    case('mm')
      conversion_factor = 1.d-3
  ! ---> TIME ---> (second)
    case('s','sec','second')
      conversion_factor = 1.d0
    case('min','minute')
      conversion_factor = 60.d0
    case('h','hr','hour')
      conversion_factor = 3600.d0
    case('d','day')
      conversion_factor = 24.d0*3600.d0 
    case('w','week')
      conversion_factor = 7.d0*24.d0*3600.d0 
    case('mo','month')
      conversion_factor = DAYS_PER_YEAR/12.d0*24.d0*3600.d0 
    case('y','yr','year')
      conversion_factor = DAYS_PER_YEAR*24.d0*3600.d0
    case('Pa.s')
      conversion_factor = 1.0
    case('cP','centiPoise')
      conversion_factor = 1.d-3
    case('P','Poise')
      conversion_factor = 1.d-1 
  ! ---> ENERGY ---> (Joule)
    case('J')   
      conversion_factor = 1.d0
    case('kJ')   
      conversion_factor = 1.d3
    case('MJ')   
      conversion_factor = 1.d6
  ! ---> ENERGY FLUX or POWER ---> (Watt)
    case('W')   
      conversion_factor = 1.d0
    case('kW')   
      conversion_factor = 1.d3
    case('MW')   
      conversion_factor = 1.d6
  ! ---> MOLAR MASS ---> (kilogram, mole)
    case('mol','mole','moles')
      conversion_factor = 1.d0
    case('kmol')
      conversion_factor = 1.d3
  ! ---> MASS ---> (kilogram, mole)
    case('ug')
      conversion_factor = 1.d-9
    case('mg')
      conversion_factor = 1.d-6
    case('g')
      conversion_factor = 1.d-3
    case('kg')
      conversion_factor = 1.d0
  ! ---> TEMPERATURE ---> (C)
    case('C','Celsius') 
      conversion_factor = 1.d0
  ! ---> PRESSURE ---> (Pascal)
    case('Pa') 
      conversion_factor = 1.d0
    case('kPa')
      conversion_factor = 1.d3
    case('MPa')
      conversion_factor = 1.d6
    case('Bar')
      conversion_factor = 1.d5
    case('psi')
      conversion_factor = 6894.757
  ! ---> CONCENTRATION ---> (M)
    case('M') 
      conversion_factor = 1.d0
    case('mM') 
      conversion_factor = 1.d-3
  ! ---> FORCE ---> (Newton)
    case('N') 
      conversion_factor = 1.d0
  ! ---> UNITLESS ---> (1)
    case('unitless','1')
      conversion_factor = 1.d0
  ! ---> NOT_ASSIGNED ---> (error)
    case('not_assigned','unknown')
      error_msg = 'Unit not assigned or unknown. Please e-mail &
                  &pflotran-dev@googlegroups.com (attn: jmfrede) &
                  &with your input file and screen output.'
      error = PETSC_TRUE
    case default
      error_msg = 'Unit [ ' // trim(unit) // ' ] not recognized when &
                  &converting to SI units. If this is not a spelling error, &
                  &then the unit is not supported.'
      error = PETSC_TRUE

  end select
  
end subroutine UnitsConvertToSI

! ************************************************************************** !

subroutine UnitsCategoryCheck(unit_user_cat,unit_internal_cat, &
                              error,error_msg)
  ! 
  ! Checks if the unit categories of the given and internal units align.
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/15/2016
  ! 

  use Option_module
  
  implicit none
 
  ! a maximum of 3 unit categories are allowed
  character(len=MAXWORDLENGTH) :: unit_user_cat(3)
  character(len=MAXWORDLENGTH) :: unit_internal_cat(3)
  PetscBool :: error
  character(len=MAXSTRINGLENGTH) :: error_msg

  PetscInt :: k, j
  PetscBool :: category_assigned(3), successful
  character(len=MAXWORDLENGTH) :: unit_cat

  category_assigned(:) = PETSC_FALSE
  error = PETSC_FALSE
  error_msg = ''
  k = 1

  do while (k < 4)
    successful = PETSC_FALSE
    unit_cat = trim(unit_user_cat(k))
    j = 1
    do while (j < 4)
      if ((trim(unit_internal_cat(j)) == unit_cat) .and. &
          (.not.category_assigned(j))) then
        category_assigned(j) = PETSC_TRUE
        successful = PETSC_TRUE
      endif
      if (successful) exit ! after first successful assignment
      j = j + 1
    enddo
    k = k + 1
  enddo

  k = 1
  do while (k < 4)
    if (.not.category_assigned(k)) then
      error_msg = 'Mismatch between the category of the given units &
                   &and the expected, internal units. Units of ' &
                   // trim(unit_internal_cat(k)) // ' were expected, but &
                   &units of ' // trim(unit_user_cat(k)) // ' were given.'
      error = PETSC_TRUE
    endif
    if (error) exit
    k = k + 1
  enddo
  
end subroutine UnitsCategoryCheck

! ************************************************************************** !

function UnitsConvertToExternal(units,units_category,option)
  ! 
  ! UnitsConvert: Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: units
  character(len=*) :: units_category
  type(option_type) :: option
  
  PetscReal :: UnitsConvertToExternal
  
  UnitsConvertToExternal = 1.d0/UnitsConvertToInternal(units, &
                                                       units_category,option)

end function UnitsConvertToExternal

end module Units_module
