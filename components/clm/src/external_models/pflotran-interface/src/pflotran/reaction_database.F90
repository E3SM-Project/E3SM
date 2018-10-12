module Reaction_Database_module

  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  public :: DatabaseRead, BasisInit
  
  public :: GetSpeciesBasisID, &
            BasisPrint

            
contains

! ************************************************************************** !

subroutine DatabaseRead(reaction,option)
  ! 
  ! Collects parameters from geochemical database
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/01/08
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Reaction_Mineral_Aux_module
  use Reaction_Mineral_module
  use Reaction_Immobile_Aux_module
  use Reaction_Immobile_module
  use Reaction_Gas_Aux_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec, cur_aq_spec2
  type(gas_species_type), pointer :: cur_gas_spec, cur_gas_spec2
  type(mineral_rxn_type), pointer :: cur_mineral, cur_mineral2
  type(immobile_species_type), pointer :: cur_immobile_spec
  type(colloid_type), pointer :: cur_colloid
  type(mineral_type), pointer :: mineral
  type(immobile_type), pointer :: immobile
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: null_name
  
  PetscBool :: flag, found, logK_error_flag
  PetscInt :: ispec, itemp, i
  PetscReal :: stoich
  PetscReal :: temp_real
  type(input_type), pointer :: input
  PetscInt :: iostat
  PetscInt :: num_nulls
  PetscInt :: num_logKs
  
  mineral => reaction%mineral
  immobile => reaction%immobile
  
  ! negate ids for use as flags
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    cur_aq_spec%id = -abs(cur_aq_spec%id)
    cur_aq_spec => cur_aq_spec%next
  enddo
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    cur_aq_spec%id = -abs(cur_aq_spec%id)
    cur_aq_spec => cur_aq_spec%next
  enddo  
  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    cur_gas_spec%id = -abs(cur_gas_spec%id)
    cur_gas_spec => cur_gas_spec%next
  enddo  
  cur_immobile_spec => immobile%list
  do
    if (.not.associated(cur_immobile_spec)) exit
    cur_immobile_spec%id = -abs(cur_immobile_spec%id)
    cur_immobile_spec => cur_immobile_spec%next
  enddo
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo
  
  if (len_trim(reaction%database_filename) < 2) then
    option%io_buffer = 'Database filename not included in input deck.'
    call printErrMsg(option)
  endif
#ifdef CLM_PFLOTRAN
  input => InputCreate(IUNIT_TEMP,trim(option%input_dir) // '/' //reaction%database_filename,option)
#else
  input => InputCreate(IUNIT_TEMP,reaction%database_filename,option)
#endif

  ! read temperatures
  call InputReadPflotranString(input,option)
  ! remove comment
  call InputReadQuotedWord(input,option,name,PETSC_TRUE)
  call InputReadInt(input,option,num_logKs)
  if (reaction%use_geothermal_hpt) then
    reaction%num_dbase_parameters = num_logKs
    call InputErrorMsg(input,option,'Number of database parameters','DATABASE')
  else
    reaction%num_dbase_temperatures = num_logKs
    call InputErrorMsg(input,option,'Number of database temperatures', &
                       'DATABASE')  
    allocate(reaction%dbase_temperatures(reaction%num_dbase_temperatures))
    reaction%dbase_temperatures = 0.d0 
   
    do itemp = 1, reaction%num_dbase_temperatures
      call InputReadDouble(input,option,reaction%dbase_temperatures(itemp))
      call InputErrorMsg(input,option,'Database temperatures','DATABASE')            
    enddo
  endif

  num_nulls = 0
  null_name = 'null'
  do ! loop over every entry in the database
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'DATABASE')

    call InputReadQuotedWord(input,option,name,PETSC_TRUE)
    ! 'null's mark the end of a section in the database.  We count these 
    ! to determine which species we are reading.
    ! --
    ! primary species
    ! null
    ! aq complexes
    ! null
    ! gases
    ! null
    ! minerals
    ! null
    ! surface complexes
    ! null
    ! --
    
    if (StringCompare(name,null_name,MAXWORDLENGTH)) then
      num_nulls = num_nulls + 1
      if (reaction%use_geothermal_hpt) then
        if (num_nulls >= 4) exit
        cycle
      else
        if (num_nulls >= 5) exit
        cycle
      endif
    endif
    
    select case(num_nulls)
      case(0,1) ! primary and secondary aq species and colloids
        cur_aq_spec => reaction%primary_species_list
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_aq_spec)) exit
          if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE
            ! change negative id to positive, indicating it was found in 
            ! database
            cur_aq_spec%id = abs(cur_aq_spec%id)
            exit
          endif
          cur_aq_spec => cur_aq_spec%next
        enddo
        if (.not.found) cur_aq_spec => reaction%secondary_species_list
        do
          if (found .or. .not.associated(cur_aq_spec)) exit
          if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
            ! change negative id to positive, indicating it was found in 
            ! database
            cur_aq_spec%id = abs(cur_aq_spec%id)
            exit
          endif
          cur_aq_spec => cur_aq_spec%next
        enddo
        ! check if a colloid
        if (.not.found) cur_colloid => reaction%colloid_list
        do
          if (found .or. .not.associated(cur_colloid)) exit
          if (StringCompare(name,cur_colloid%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
            ! change negative id to positive, indicating it was found in 
            ! database
            cur_colloid%id = abs(cur_colloid%id)

            ! skip the Debye-Huckel ion size parameter (a0)
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'Colloid skip a0','DATABASE')            
            ! skip the valence
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'Colloid skip Z','DATABASE')            
            ! read the molar weight
            call InputReadDouble(input,option,cur_colloid%molar_weight)
            call InputErrorMsg(input,option,'Colloid molar weight','DATABASE')
            
            cycle ! avoid the aqueous species parameters below
          endif
          cur_colloid => cur_colloid%next
        enddo
        ! check if immobile
        if (.not.found) cur_immobile_spec => immobile%list
        do
          if (found .or. .not.associated(cur_immobile_spec)) exit
          if (StringCompare(name,cur_immobile_spec%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
            ! change negative id to positive, indicating it was found in 
            ! database
            cur_immobile_spec%id = abs(cur_immobile_spec%id)

            ! skip the Debye-Huckel ion size parameter (a0)
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'Immobile skip a0','DATABASE')            
            ! skip the valence
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'Immobile skip Z','DATABASE')            
            ! read the molar weight
            call InputReadDouble(input,option,cur_immobile_spec%molar_weight)
            call InputErrorMsg(input,option,'Immobile molar weight','DATABASE')
            
            cycle ! avoid the aqueous species parameters below
          endif
          cur_immobile_spec => cur_immobile_spec%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        if (num_nulls > 0) then ! secondary species in database
          ! create aqueous equilibrium reaction
          if (.not.associated(cur_aq_spec%dbaserxn)) &
            cur_aq_spec%dbaserxn => DatabaseRxnCreate()
          ! read the number of primary species in secondary rxn
          call InputReadInt(input,option,cur_aq_spec%dbaserxn%nspec)
          call InputErrorMsg(input,option,'Number of species in aqueous ', &
                          'complex DATABASE')  
          ! allocate arrays for rxn
          allocate(cur_aq_spec%dbaserxn%spec_name(cur_aq_spec%dbaserxn%nspec))
          cur_aq_spec%dbaserxn%spec_name = ''
          allocate(cur_aq_spec%dbaserxn%stoich(cur_aq_spec%dbaserxn%nspec))
          cur_aq_spec%dbaserxn%stoich = 0.d0
          allocate(cur_aq_spec%dbaserxn%logK(num_logKs))
          cur_aq_spec%dbaserxn%logK = 0.d0
          ! read in species and stoichiometries
          do ispec = 1, cur_aq_spec%dbaserxn%nspec
            call InputReadDouble(input,option, &
                                 cur_aq_spec%dbaserxn%stoich(ispec))
            call InputErrorMsg(input,option,'EQRXN species stoichiometry', &
                               'DATABASE')            
            call InputReadQuotedWord(input,option, &
                            cur_aq_spec%dbaserxn%spec_name(ispec),PETSC_TRUE)
            call InputErrorMsg(input,option,'EQRXN species name','DATABASE')            
          enddo
          !note: logKs read are pK so that K is in the denominator (i.e. Q/K)
          do itemp = 1, num_logKs
            call InputReadDouble(input,option,cur_aq_spec%dbaserxn%logK(itemp))
            call InputErrorMsg(input,option,'EQRXN logKs','DATABASE') 
          enddo
        endif
        ! read the Debye-Huckel ion size parameter (a0)
        call InputReadDouble(input,option,cur_aq_spec%a0)
        call InputErrorMsg(input,option,'AQ Species a0','DATABASE')            
        ! read the valence
        call InputReadDouble(input,option,cur_aq_spec%Z)
        call InputErrorMsg(input,option,'AQ Species Z','DATABASE')            
        ! read the molar weight
        call InputReadDouble(input,option,cur_aq_spec%molar_weight)
        call InputErrorMsg(input,option,'AQ Species molar weight','DATABASE')
        
                    
      case(2) ! gas species
        cur_gas_spec => reaction%gas%list
        if (.not.associated(cur_gas_spec)) cycle
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_gas_spec)) exit
          if (StringCompare(name,cur_gas_spec%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
            ! change negative id to positive, indicating it was found in 
            ! database
            cur_gas_spec%id = abs(cur_gas_spec%id)
            exit
          endif
          cur_gas_spec => cur_gas_spec%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        ! read the molar volume
        call InputReadDouble(input,option,cur_gas_spec%molar_volume)
        call InputErrorMsg(input,option,'GAS molar volume','DATABASE')
        ! convert from cm^3/mol to m^3/mol
        cur_gas_spec%molar_volume = cur_gas_spec%molar_volume*1.d-6
        ! create aqueous equilibrium reaction
        if (.not.associated(cur_gas_spec%dbaserxn)) &
          cur_gas_spec%dbaserxn => DatabaseRxnCreate()
        ! read the number of aqueous species in secondary rxn
        call InputReadInt(input,option,cur_gas_spec%dbaserxn%nspec)
        call InputErrorMsg(input,option,'Number of species in gas reaction', &
                        'DATABASE')  
        ! allocate arrays for rxn
        allocate(cur_gas_spec%dbaserxn%spec_name(cur_gas_spec%dbaserxn%nspec))
        cur_gas_spec%dbaserxn%spec_name = ''
        allocate(cur_gas_spec%dbaserxn%stoich(cur_gas_spec%dbaserxn%nspec))
        cur_gas_spec%dbaserxn%stoich = 0.d0
        allocate(cur_gas_spec%dbaserxn%logK(num_logKs))
        cur_gas_spec%dbaserxn%logK = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_gas_spec%dbaserxn%nspec
          call InputReadDouble(input,option, &
                               cur_gas_spec%dbaserxn%stoich(ispec))
          call InputErrorMsg(input,option,'GAS species stoichiometry', &
                             'DATABASE')            
          call InputReadQuotedWord(input,option, &
                                   cur_gas_spec%dbaserxn%spec_name(ispec), &
                                   PETSC_TRUE)
          call InputErrorMsg(input,option,'GAS species name','DATABASE')            
        enddo
        !note: logKs read are pK so that K is in the denominator (i.e. Q/K)
        do itemp = 1, num_logKs
          call InputReadDouble(input,option,cur_gas_spec%dbaserxn%logK(itemp))
          call InputErrorMsg(input,option,'GAS logKs','DATABASE')            
        enddo
        ! read the molar weight
        call InputReadDouble(input,option,cur_gas_spec%molar_weight)
        call InputErrorMsg(input,option,'GAS molar weight','DATABASE')     
        
               
      case(3) ! minerals
        cur_mineral => mineral%mineral_list
        if (.not.associated(cur_mineral)) cycle
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_mineral)) exit
          if (StringCompare(name,cur_mineral%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
            ! change negative id to positive, indicating it was found in 
            ! database
            cur_mineral%id = abs(cur_mineral%id)
            exit
          endif
          cur_mineral => cur_mineral%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        call MineralReadFromDatabase(cur_mineral, &
                                     num_logKs,input, &
                                     option)
      
    end select
    
  enddo
  
  ! check for duplicate species
  flag = PETSC_FALSE
 
  ! aqueous primary species
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    
    ! aqueous primary species
    cur_aq_spec2 => cur_aq_spec%next
    do
      if (.not.associated(cur_aq_spec2)) exit
      if (cur_aq_spec%id /= cur_aq_spec2%id .and. &
          StringCompare(cur_aq_spec%name, &
                          cur_aq_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = &
                 'Aqueous primary species (' // trim(cur_aq_spec%name) // &
                 ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_aq_spec2 => cur_aq_spec2%next
    enddo

    cur_aq_spec2 => reaction%secondary_species_list
    do
      if (.not.associated(cur_aq_spec2)) exit
      if (StringCompare(cur_aq_spec%name, &
                          cur_aq_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous primary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated as secondary species in input file.'
        call printMsg(option)                          
      endif
      cur_aq_spec2 => cur_aq_spec2%next
    enddo

    cur_gas_spec2 => reaction%gas%list
    do
      if (.not.associated(cur_gas_spec2)) exit
      if (StringCompare(cur_aq_spec%name, &
                          cur_gas_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous primary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated as gas species in input file.'
        call printMsg(option)                          
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    cur_aq_spec => cur_aq_spec%next  
  enddo
  
  ! aqueous secondary species
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    
    ! already checked against primary
    ! aqueous secondary species
    cur_aq_spec2 => cur_aq_spec%next
    do
      if (.not.associated(cur_aq_spec2)) exit
      if (cur_aq_spec%id /= cur_aq_spec2%id .and. &
          StringCompare(cur_aq_spec%name, &
                          cur_aq_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous secondary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_aq_spec2 => cur_aq_spec2%next
    enddo

    cur_gas_spec2 => reaction%gas%list
    do
      if (.not.associated(cur_gas_spec2)) exit
      if (StringCompare(cur_aq_spec%name, &
                          cur_gas_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Aqueous secondary species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated as gas species in input file.'
        call printMsg(option)                          
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    cur_aq_spec => cur_aq_spec%next  
  enddo
  
  ! gas species
  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_aq_spec)) exit
    
    ! already checked against primary
    ! already checked against secondary
    ! gas species
    cur_gas_spec2 => cur_gas_spec%next
    do
      if (.not.associated(cur_gas_spec2)) exit
      if (cur_gas_spec%id /= cur_gas_spec2%id .and. &
          StringCompare(cur_aq_spec%name, &
                          cur_gas_spec2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Gas species (' // &
                           trim(cur_aq_spec%name) // &
                           ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_gas_spec2 => cur_gas_spec2%next
    enddo
    cur_aq_spec => cur_aq_spec%next  
  enddo
  
  ! minerals
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral2 => cur_mineral%next
    do
      if (.not.associated(cur_mineral2)) exit
      if (cur_mineral%id /= cur_mineral2%id .and. &
          StringCompare(cur_mineral%name, &
                          cur_mineral2%name,MAXWORDLENGTH)) then
        flag = PETSC_TRUE
        option%io_buffer = 'Mineral (' // &
                           trim(cur_mineral%name) // &
                           ') duplicated in input file.'
        call printMsg(option)                          
      endif
      cur_mineral2 => cur_mineral2%next
    enddo
    cur_mineral => cur_mineral%next
  enddo
  
  if (flag) call printErrMsg(option,'Species duplicated in input file.')

  ! check that all species, etc. were read
  ! also check whether legitimate logK values exist if non-isothermal and
  ! a database reaction exists
  flag = PETSC_FALSE
  logK_error_flag = PETSC_FALSE
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = 'Aqueous primary species (' // &
               trim(cur_aq_spec%name) // &
               ') not found in database.'
      call printMsg(option)
    endif
    if (.not.reaction%use_geothermal_hpt) then
      if (.not.DatabaseCheckLegitimateLogKs(cur_aq_spec%dbaserxn, &
                                            cur_aq_spec%name, &
                                            reaction%dbase_temperatures, &
                                            option)) then
        logK_error_flag = PETSC_TRUE
      endif
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = &
               'Aqueous secondary species (' // trim(cur_aq_spec%name) // &
               ') not found in database.'
      call printMsg(option)
    endif
    if (.not.reaction%use_geothermal_hpt) then
      if (.not.DatabaseCheckLegitimateLogKs(cur_aq_spec%dbaserxn, &
                                            cur_aq_spec%name, &
                                            reaction%dbase_temperatures, &
                                            option)) then
        logK_error_flag = PETSC_TRUE
      endif
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo  
  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    if (cur_gas_spec%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = 'Gas species (' // trim(cur_gas_spec%name) // &
                         ') not found in database.'
      call printMsg(option)
    endif
    if (.not.reaction%use_geothermal_hpt) then
      if (.not.DatabaseCheckLegitimateLogKs(cur_gas_spec%dbaserxn, &
                                            cur_gas_spec%name, &
                                            reaction%dbase_temperatures, &
                                            option)) then
        logK_error_flag = PETSC_TRUE
      endif
    endif
    cur_gas_spec => cur_gas_spec%next
  enddo  
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = 'Mineral (' // trim(cur_mineral%name) // &
               ') not found in database.'
      call printErrMsg(option)
    endif
    if (.not.reaction%use_geothermal_hpt) then
      if (.not.DatabaseCheckLegitimateLogKs(cur_mineral%dbaserxn, &
                                            cur_mineral%name, &
                                            reaction%dbase_temperatures, &
                                            option)) then
        logK_error_flag = PETSC_TRUE
      endif
    endif
    cur_mineral => cur_mineral%next
  enddo

  if (flag) call printErrMsg(option,'Species not found in database.')
  if (.not.option%use_isothermal) then
    !geh: only stop if running with temperature dependent log Ks.
    if (logK_error_flag) then
      option%io_buffer = 'Non-isothermal reactions not possible due to &
        &missing logKs in database.'
      call printErrMsg(option)
    endif
  endif

  call InputDestroy(input)
  
end subroutine DatabaseRead

! ************************************************************************** !

subroutine BasisInit(reaction,option)
  ! 
  ! Initializes the basis for geochemistry
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/01/08
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Utility_module
  use Input_Aux_module
  
  use Reaction_Mineral_Aux_module
  use Reaction_Immobile_Aux_module
  use Reaction_Gas_Aux_module
  
  use Reaction_Sandbox_module

  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec
  type(aq_species_type), pointer :: cur_pri_aq_spec
  type(aq_species_type), pointer :: cur_sec_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_rxn_type), pointer :: cur_mineral
  type(aq_species_type), pointer :: cur_sec_aq_spec1
  type(aq_species_type), pointer :: cur_sec_aq_spec2
  type(gas_species_type), pointer :: cur_gas_spec1
  type(gas_species_type), pointer :: cur_gas_spec2
  type(immobile_species_type), pointer :: cur_immobile_spec
  type(ion_exchange_rxn_type), pointer :: cur_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cur_cation
  type(general_rxn_type), pointer :: cur_general_rxn
  type(radioactive_decay_rxn_type), pointer :: cur_radiodecay_rxn

  type(immobile_decay_rxn_type), pointer :: cur_immobile_decay_rxn
  type(kd_rxn_type), pointer :: cur_kd_rxn, sec_cont_cur_kd_rxn
  type(colloid_type), pointer :: cur_colloid
  type(database_rxn_type), pointer :: dbaserxn
  type(transition_state_rxn_type), pointer :: tstrxn
  type(transition_state_prefactor_type), pointer :: cur_prefactor
  type(ts_prefactor_species_type), pointer :: cur_prefactor_species
  type(mineral_type), pointer :: mineral
  type(immobile_type), pointer :: immobile

  character(len=MAXWORDLENGTH), allocatable :: old_basis_names(:)
  character(len=MAXWORDLENGTH), allocatable :: new_basis_names(:)

  character(len=MAXWORDLENGTH), parameter :: h2oname = 'H2O'
  character(len=MAXWORDLENGTH) :: word, word2
  character(len=MAXSTRINGLENGTH) :: string, string2

  PetscInt, parameter :: h2o_id = 1

  PetscReal :: logK(reaction%num_dbase_temperatures)
  PetscReal, allocatable :: transformation(:,:), old_basis(:,:), new_basis(:,:)
  PetscReal, allocatable :: stoich_new(:), stoich_prev(:), logKvector(:,:)
  PetscInt, allocatable :: indices(:)
  
  PetscReal, allocatable :: pri_matrix(:,:), sec_matrix(:,:)
  PetscReal, allocatable :: sec_matrix_inverse(:,:)
  PetscReal, allocatable :: stoich_matrix(:,:)
  PetscReal, allocatable :: unit_vector(:)
  character(len=MAXWORDLENGTH), allocatable :: pri_names(:)
  character(len=MAXWORDLENGTH), allocatable :: sec_names(:)
  character(len=MAXWORDLENGTH), allocatable :: gas_names(:)
  PetscReal, allocatable :: logKvector_swapped(:,:)
  PetscBool, allocatable :: colloid_species_flag(:)
  PetscReal :: value
  
  PetscInt :: ispec, itemp
  PetscInt :: spec_id
  PetscInt :: ncomp_h2o, ncomp_secondary
  PetscInt :: icount_old, icount_new, icount, icount2, icount3
  PetscInt :: i, j, irow, icol
  PetscInt :: icomp, icplx, irxn, ieqrxn
  PetscInt :: ipri_spec, isec_spec, imnrl, ikinmnrl, icoll
  PetscInt :: i_old, i_new
  PetscInt :: isrfcplx
  PetscInt :: ication
  PetscInt :: temp_int
  PetscReal :: scale
  PetscReal :: temp_high, temp_low
  PetscInt :: itemp_high, itemp_low
  PetscInt :: species_count, max_species_count
  PetscInt :: max_monod_count, max_inhibition_count
  PetscInt :: monod_count, inhibition_count, activation_energy_count
  PetscInt :: forward_count, max_forward_count
  PetscInt :: backward_count, max_backward_count
  PetscInt :: max_aq_species
  PetscInt :: max_num_prefactors, max_num_prefactor_species
  
  PetscBool :: compute_new_basis
  PetscBool :: found
  PetscErrorCode :: ierr
  PetscInt :: num_logKs
  
  mineral => reaction%mineral
  immobile => reaction%immobile
    
  if (reaction%use_geothermal_hpt) then
    num_logKs = reaction%num_dbase_parameters
  else
    num_logKs = reaction%num_dbase_temperatures
  endif
    
! get database temperature based on REFERENCE_TEMPERATURE
  if (option%reference_temperature <= 0.01d0) then
    reaction%debyeA = 0.4939d0 
    reaction%debyeB = 0.3253d0 
    reaction%debyeBdot = 0.0374d0
  else if (option%reference_temperature > 0.d0 .and. &
           option%reference_temperature <= 25.d0) then
    temp_low = 0.d0
    temp_high = 25.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.5114d0,0.4939d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3288d0,0.3253d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0410d0,0.0374d0,reaction%debyeBdot)
  else if (option%reference_temperature > 25.d0 .and. &
           option%reference_temperature <= 60.d0) then
    temp_low = 25.d0
    temp_high = 60.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.5465d0,0.5114d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3346d0,0.3288d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0440d0,0.0410d0,reaction%debyeBdot)
  else if (option%reference_temperature > 60.d0 .and. &
           option%reference_temperature <= 100.d0) then
    temp_low = 60.d0
    temp_high = 100.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.5995d0,0.5465d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3421d0,0.3346d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0460d0,0.0440d0,reaction%debyeBdot)
  else if (option%reference_temperature > 100.d0 .and. &
           option%reference_temperature <= 150.d0) then
    temp_low = 100.d0
    temp_high = 150.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.6855d0,0.5995d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3525d0,0.3421d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0470d0,0.0460d0,reaction%debyeBdot)
  else if (option%reference_temperature > 150.d0 .and. &
           option%reference_temperature <= 200.d0) then
    temp_low = 150.d0
    temp_high = 200.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.7994d0,0.6855d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3639d0,0.3525d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0470d0,0.0470d0,reaction%debyeBdot)
  else if (option%reference_temperature > 200.d0 .and. &
           option%reference_temperature <= 250.d0) then
    temp_low = 200.d0
    temp_high = 250.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.9593d0,0.7994d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3766d0,0.3639d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0340d0,0.0470d0,reaction%debyeBdot)
  else if (option%reference_temperature > 250.d0 .and. &
           option%reference_temperature <= 300.d0) then
    temp_low = 250.d0
    temp_high = 300.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     1.2180d0,0.9593d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3925d0,0.3766d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0000d0,0.0340d0,reaction%debyeBdot)
  else if (option%reference_temperature > 300.d0 .and. &
           option%reference_temperature <= 350.d0) then
    temp_low = 300.d0
    temp_high = 350.d0
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     1.2180d0,1.2180d0,reaction%debyeA)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.3925d0,0.3925d0,reaction%debyeB)
    call Interpolate(temp_high,temp_low,option%reference_temperature, &
                     0.0000d0,0.0000d0,reaction%debyeBdot)
  else if (option%reference_temperature > 350.d0) then
    reaction%debyeA = 1.2180d0 
    reaction%debyeB = 0.3925d0 
    reaction%debyeBdot = 0.0000d0
  endif
  
  if (.not.reaction%act_coef_use_bdot) then
    reaction%debyeBdot = 0.d0
  endif

  if (.not. reaction%use_geothermal_hpt) then
    if (option%reference_temperature <= reaction%dbase_temperatures(1)) then
      itemp_low = 1
      itemp_high = 1
      temp_low = reaction%dbase_temperatures(itemp_low)
      temp_high = reaction%dbase_temperatures(itemp_high)
    else if (option%reference_temperature > &
             reaction%dbase_temperatures(reaction%num_dbase_temperatures)) then
      itemp_low = reaction%num_dbase_temperatures
      itemp_high = reaction%num_dbase_temperatures
      temp_low = reaction%dbase_temperatures(itemp_low)
      temp_high = reaction%dbase_temperatures(itemp_high)
    else
      do itemp = 1, reaction%num_dbase_temperatures-1
        itemp_low = itemp
        itemp_high = itemp+1
        temp_low = reaction%dbase_temperatures(itemp_low)
        temp_high = reaction%dbase_temperatures(itemp_high)
        if (option%reference_temperature > temp_low .and. &
            option%reference_temperature <= temp_high) then
          exit
        endif
      enddo
    endif
  endif
  
  reaction%naqcomp = GetPrimarySpeciesCount(reaction)
  reaction%neqcplx = GetSecondarySpeciesCount(reaction)
  reaction%gas%ngas = GasGetCount(reaction%gas%list,ACTIVE_AND_PASSIVE_GAS)
  reaction%nimcomp = GetImmobileCount(reaction)
  reaction%ncoll = GetColloidCount(reaction)
  ! set to naqcomp for now, will be adjusted later
  reaction%ncollcomp = reaction%naqcomp 
  
  reaction%offset_aqueous = 0
  reaction%offset_immobile = reaction%offset_aqueous + reaction%naqcomp
  reaction%offset_colloid = reaction%offset_immobile + reaction%nimcomp
  reaction%offset_collcomp = reaction%offset_colloid + reaction%ncoll

  ! account for H2O in the basis by adding 1
  ncomp_h2o = reaction%naqcomp+1
  
  allocate(old_basis_names(ncomp_h2o+reaction%neqcplx))
  allocate(new_basis_names(ncomp_h2o))
  old_basis_names = ''
  new_basis_names = ''
  
  call BasisPrint(reaction,'Initial Basis',option)
  
  !--------------------------------------------

  ! for now, remove equilibrium reactions from any primary species that are
  ! flagged as "redox"
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    if (cur_aq_spec%is_redox .and. associated(cur_aq_spec%dbaserxn)) then
      call DatabaseRxnDestroy(cur_aq_spec%dbaserxn)
    endif
    cur_aq_spec => cur_aq_spec%next
  enddo

  ncomp_secondary = reaction%neqcplx+reaction%gas%ngas
  
  ! check to ensure that the number of secondary aqueous and gas species
  ! (i.e. those with a reaction defined from the database) is equal to the 
  ! number of reactions read from the database.  If not, send an error 
  ! message.  
  
  icount = 0
  cur_pri_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    if (associated(cur_pri_aq_spec%dbaserxn)) then
      icount = icount + 1
    endif
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo

  cur_sec_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_sec_aq_spec)) exit
    if (associated(cur_sec_aq_spec%dbaserxn)) then
      icount = icount + 1
    endif
    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo

  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    if (associated(cur_gas_spec%dbaserxn)) then
      icount = icount + 1
    endif
    cur_gas_spec => cur_gas_spec%next
  enddo
  
  if (icount /= ncomp_secondary) then
    if (icount < ncomp_secondary) then
      option%io_buffer = 'Too few reactions read from database for & 
        &number of secondary species defined.'
    else
      option%io_buffer = 'Too many reactions read from database for &
        &number of secondary species defined.  Perhaps REDOX &
        &SPECIES need to be defined?'
    endif
    call printErrMsg(option)
  endif
  
  allocate(pri_matrix(ncomp_secondary,ncomp_h2o))
  pri_matrix = 0.d0
  allocate(pri_names(ncomp_h2o))
  pri_names = ''
  allocate(sec_matrix(ncomp_secondary,ncomp_secondary))
  sec_matrix = 0.d0
  allocate(sec_names(reaction%neqcplx))
  sec_names = ''
  allocate(gas_names(reaction%gas%ngas))
  gas_names = ''
  
  allocate(logKvector(num_logKs,ncomp_secondary))
  logKvector = 0.d0
  
  ! fill in names
  icount = 1
  pri_names(icount) = h2oname
  cur_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    icount = icount + 1
    pri_names(icount) = cur_aq_spec%name
    cur_aq_spec => cur_aq_spec%next
  enddo
  icount = 0
  cur_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_aq_spec)) exit
    icount = icount + 1
    sec_names(icount) = cur_aq_spec%name
    cur_aq_spec => cur_aq_spec%next
  enddo
  icount= 0
  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    icount = icount + 1
    gas_names(icount) = cur_gas_spec%name
    cur_gas_spec => cur_gas_spec%next
  enddo
  
  ! fill in matrices
  icount = 0
  cur_pri_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    if (associated(cur_pri_aq_spec%dbaserxn)) then
      icount = icount + 1
      logKvector(:,icount) = cur_pri_aq_spec%dbaserxn%logK
      i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                            cur_pri_aq_spec%name, &
                            cur_pri_aq_spec%name, &
                            pri_names,sec_names,gas_names)
      if (i < 0) then
        option%io_buffer = 'Primary species ' // &
                 trim(cur_pri_aq_spec%name) // &
                 ' found in secondary or gas list.'
        call printErrMsg(option)
      endif
      pri_matrix(icount,i) = -1.d0
      do ispec=1,cur_pri_aq_spec%dbaserxn%nspec
        i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                              cur_pri_aq_spec%name, &
                              cur_pri_aq_spec%dbaserxn%spec_name(ispec), &
                              pri_names,sec_names,gas_names)
        if (i > 0) then
          pri_matrix(icount,i) = cur_pri_aq_spec%dbaserxn%stoich(ispec)
        else
          sec_matrix(icount,-i) = cur_pri_aq_spec%dbaserxn%stoich(ispec)
        endif
      enddo
    endif
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo

  cur_sec_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_sec_aq_spec)) exit
    if (associated(cur_sec_aq_spec%dbaserxn)) then
      icount = icount + 1
      logKvector(:,icount) = cur_sec_aq_spec%dbaserxn%logK
      i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                            cur_sec_aq_spec%name, &
                            cur_sec_aq_spec%name, &
                            pri_names,sec_names,gas_names)
      if (i > 0) then
        option%io_buffer = 'Secondary aqueous species ' // &
                 trim(cur_sec_aq_spec%name) // &
                 ' found in primary species list.'
        call printErrMsg(option)
      endif
      sec_matrix(icount,-i) = -1.d0
      do ispec=1,cur_sec_aq_spec%dbaserxn%nspec
        i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                              cur_sec_aq_spec%name, &
                              cur_sec_aq_spec%dbaserxn%spec_name(ispec), &
                              pri_names,sec_names,gas_names)
        if (i > 0) then
          pri_matrix(icount,i) = cur_sec_aq_spec%dbaserxn%stoich(ispec)
        else
          sec_matrix(icount,-i) = cur_sec_aq_spec%dbaserxn%stoich(ispec)
        endif
      enddo
    endif
    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo

  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    if (associated(cur_gas_spec%dbaserxn)) then
      icount = icount + 1
      logKvector(:,icount) = cur_gas_spec%dbaserxn%logK
      i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                            cur_gas_spec%name, &
                            cur_gas_spec%name, &
                            pri_names,sec_names,gas_names)
      if (i > 0) then
        option%io_buffer = 'Gas species ' // &
                 trim(cur_gas_spec%name) // &
                 ' found in primary species list.'
        call printErrMsg(option)
      endif
      sec_matrix(icount,-i) = -1.d0
      do ispec=1,cur_gas_spec%dbaserxn%nspec
        i = GetSpeciesBasisID(reaction,option,ncomp_h2o, &
                              cur_gas_spec%name, &
                              cur_gas_spec%dbaserxn%spec_name(ispec), &
                              pri_names,sec_names,gas_names)
        if (i > 0) then
          pri_matrix(icount,i) = cur_gas_spec%dbaserxn%stoich(ispec)
        else
          sec_matrix(icount,-i) = cur_gas_spec%dbaserxn%stoich(ispec)
        endif
      enddo
    endif
    cur_gas_spec => cur_gas_spec%next
  enddo

  allocate(indices(ncomp_secondary))
  indices = 0
  allocate(unit_vector(ncomp_secondary))
  unit_vector = 0.d0
  allocate(sec_matrix_inverse(ncomp_secondary,ncomp_secondary))
  sec_matrix_inverse = 0.d0
 
  call ludcmp(sec_matrix,ncomp_secondary,indices,temp_int)
  do ispec = 1, ncomp_secondary
    unit_vector = 0.d0
    unit_vector(ispec) = 1.d0
    call lubksb(sec_matrix,ncomp_secondary,indices,unit_vector)
    sec_matrix_inverse(:,ispec) = unit_vector(:)
  enddo

  ! invert the secondary species matrix
  allocate(stoich_matrix(ncomp_secondary,ncomp_h2o))
  stoich_matrix = 0.d0
  do j = 1, ncomp_h2o
    do i = 1, ncomp_secondary
      do ispec = 1, ncomp_secondary
        stoich_matrix(i,j) = stoich_matrix(i,j) + &
          sec_matrix_inverse(i,ispec)*pri_matrix(ispec,j)
      enddo
    enddo
  enddo
  stoich_matrix = -1.d0*stoich_matrix

  allocate(logKvector_swapped(num_logKs,ncomp_secondary))
  logKvector_swapped = 0.d0
  
  do j = 1, ncomp_secondary
    do i = 1, num_logKs
      logKvector_swapped(i,j) = logKvector_swapped(i,j) - &
        dot_product(sec_matrix_inverse(j,1:ncomp_secondary), &
                    logKvector(i,1:ncomp_secondary))
    enddo
  enddo
    
  deallocate(pri_matrix)
  deallocate(sec_matrix)
  deallocate(indices)
  deallocate(unit_vector)
  deallocate(sec_matrix_inverse)
  deallocate(logKvector)
  
  cur_pri_aq_spec => reaction%primary_species_list
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    if (associated(cur_pri_aq_spec%dbaserxn)) then
      call DatabaseRxnDestroy(cur_pri_aq_spec%dbaserxn)
    endif
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo

  icount = 0
  cur_sec_aq_spec => reaction%secondary_species_list
  do
    if (.not.associated(cur_sec_aq_spec)) exit
    icount = icount + 1
    ! destory old reaction
    call DatabaseRxnDestroy(cur_sec_aq_spec%dbaserxn)
    ! allocate new
    cur_sec_aq_spec%dbaserxn => DatabaseRxnCreate()

    ! count # of species in reaction
    icount2 = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        cur_sec_aq_spec%dbaserxn%nspec = cur_sec_aq_spec%dbaserxn%nspec + 1
      endif
    enddo
    
    allocate(cur_sec_aq_spec%dbaserxn%stoich(cur_sec_aq_spec%dbaserxn%nspec))
    cur_sec_aq_spec%dbaserxn%stoich = 0.d0
    allocate(cur_sec_aq_spec%dbaserxn%spec_name(cur_sec_aq_spec%dbaserxn%nspec))
    cur_sec_aq_spec%dbaserxn%spec_name = ''
    allocate(cur_sec_aq_spec%dbaserxn%spec_ids(cur_sec_aq_spec%dbaserxn%nspec))
    cur_sec_aq_spec%dbaserxn%spec_ids = 0
    allocate(cur_sec_aq_spec%dbaserxn%logK(num_logKs))
    cur_sec_aq_spec%dbaserxn%logK = 0.d0

    ispec = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        ispec = ispec + 1
        cur_sec_aq_spec%dbaserxn%spec_name(ispec) = pri_names(icol)
        cur_sec_aq_spec%dbaserxn%stoich(ispec) = stoich_matrix(icount,icol)
        cur_sec_aq_spec%dbaserxn%spec_ids(ispec) = icol
      endif
    enddo

    cur_sec_aq_spec%dbaserxn%logK = logKvector_swapped(:,icount)

    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo

  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    icount = icount + 1
    ! destory old reaction
    call DatabaseRxnDestroy(cur_gas_spec%dbaserxn)
    ! allocate new
    cur_gas_spec%dbaserxn => DatabaseRxnCreate()

    ! count # of species in reaction
    icount2 = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        cur_gas_spec%dbaserxn%nspec = cur_gas_spec%dbaserxn%nspec + 1
      endif
    enddo
    
    allocate(cur_gas_spec%dbaserxn%stoich(cur_gas_spec%dbaserxn%nspec))
    cur_gas_spec%dbaserxn%stoich = 0.d0
    allocate(cur_gas_spec%dbaserxn%spec_name(cur_gas_spec%dbaserxn%nspec))
    cur_gas_spec%dbaserxn%spec_name = ''
    allocate(cur_gas_spec%dbaserxn%spec_ids(cur_gas_spec%dbaserxn%nspec))
    cur_gas_spec%dbaserxn%spec_ids = 0
    allocate(cur_gas_spec%dbaserxn%logK(num_logKs))
    cur_gas_spec%dbaserxn%logK = 0.d0

    ispec = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        ispec = ispec + 1
        cur_gas_spec%dbaserxn%spec_name(ispec) = pri_names(icol)
        cur_gas_spec%dbaserxn%stoich(ispec) = stoich_matrix(icount,icol)
        cur_gas_spec%dbaserxn%spec_ids(ispec) = icol
      endif
    enddo

    cur_gas_spec%dbaserxn%logK = logKvector_swapped(:,icount)

    cur_gas_spec => cur_gas_spec%next
  enddo

  new_basis_names = pri_names

  deallocate(stoich_matrix)
  deallocate(logKvector_swapped)

  deallocate(pri_names)
  deallocate(sec_names)
  deallocate(gas_names)

  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
    
  ! first off, lets remove all the secondary gases from all other reactions
  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    
    ! gases in mineral reactions
    cur_mineral => mineral%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%dbaserxn%nspec) exit
          if (StringCompare(cur_gas_spec%name, &
                              cur_mineral%dbaserxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn(cur_gas_spec%name, &
                                             cur_gas_spec%dbaserxn, &
                                             cur_mineral%dbaserxn, &
                                             scale)
!geh             cur_mineral%dbaserxn%logK = cur_mineral%dbaserxn%logK &
!geh                                       + scale*cur_gas_spec%dbaserxn%logK
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo
    nullify(cur_mineral)

    cur_gas_spec => cur_gas_spec%next
  enddo

  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)


  ! secondary aqueous species
  cur_sec_aq_spec => reaction%secondary_species_list
  do

    if (.not.associated(cur_sec_aq_spec)) exit
    
    ! secondary aqueous species in mineral reactions
    cur_mineral => mineral%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%dbaserxn%nspec) exit
          if (StringCompare(cur_sec_aq_spec%name, &
                              cur_mineral%dbaserxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn(cur_sec_aq_spec%name, &
                                             cur_sec_aq_spec%dbaserxn, &
                                             cur_mineral%dbaserxn, &
                                             scale)
!geh            cur_mineral%dbaserxn%logK = cur_mineral%dbaserxn%logK &
!geh                                      + scale*cur_sec_aq_spec%dbaserxn%logK
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo
    
    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo
  
  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)

  ! substitute new basis into mineral and surface complexation rxns,
  ! if necessary
  cur_mineral => mineral%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    if (.not.associated(cur_mineral%dbaserxn%spec_ids)) then
      allocate(cur_mineral%dbaserxn%spec_ids(cur_mineral%dbaserxn%nspec))
      cur_mineral%dbaserxn%spec_ids = 0
    endif

    call BasisAlignSpeciesInRxn(ncomp_h2o,new_basis_names, &
                                cur_mineral%dbaserxn%nspec, &
                                cur_mineral%dbaserxn%spec_name, &
                                cur_mineral%dbaserxn%stoich, &
                                cur_mineral%dbaserxn%spec_ids, &
                                cur_mineral%name,option)     
    cur_mineral => cur_mineral%next
  enddo  

  ! fill reaction arrays, swapping if necessary
  if (associated(reaction%primary_species_names)) &
    deallocate(reaction%primary_species_names)

  allocate(reaction%primary_species_names(reaction%naqcomp))
  reaction%primary_species_names = ''

  allocate(reaction%primary_species_print(reaction%naqcomp))
  reaction%primary_species_print = PETSC_FALSE

  allocate(reaction%primary_spec_Z(reaction%naqcomp))
  reaction%primary_spec_Z = 0.d0

  allocate(reaction%primary_spec_molar_wt(reaction%naqcomp))
  reaction%primary_spec_molar_wt = 0.d0

  allocate(reaction%primary_spec_a0(reaction%naqcomp))
  reaction%primary_spec_a0 = 0.d0

  allocate(reaction%kd_print(reaction%naqcomp))
  reaction%kd_print = PETSC_FALSE
  if (reaction%nsorb > 0) then
    allocate(reaction%total_sorb_print(reaction%naqcomp))
    reaction%total_sorb_print = PETSC_FALSE
  endif
  
    ! pack in reaction arrays
  cur_pri_aq_spec => reaction%primary_species_list
  ispec = 1
  do
    if (.not.associated(cur_pri_aq_spec)) exit
    reaction%primary_species_names(ispec) = cur_pri_aq_spec%name
    reaction%primary_spec_Z(ispec) = cur_pri_aq_spec%Z
    reaction%primary_spec_molar_wt(ispec) = cur_pri_aq_spec%molar_weight
    reaction%primary_spec_a0(ispec) = cur_pri_aq_spec%a0
    reaction%primary_species_print(ispec) = cur_pri_aq_spec%print_me .or. &
                                            reaction%print_all_primary_species
    reaction%kd_print(ispec) = (cur_pri_aq_spec%print_me .or. &
                                reaction%print_all_primary_species) .and. &
                                reaction%print_kd
    if (reaction%nsorb > 0) then
      reaction%total_sorb_print(ispec) = (cur_pri_aq_spec%print_me .or. &
                                  reaction%print_all_primary_species) .and. &
                                  reaction%print_total_sorb
    endif
    ispec = ispec + 1
    cur_pri_aq_spec => cur_pri_aq_spec%next
  enddo
  nullify(cur_pri_aq_spec)
  ispec = -1 ! to catch bugs
  
  ! secondary aqueous complexes
  reaction%neqcplx = GetSecondarySpeciesCount(reaction)
  
  if (reaction%neqcplx > 0) then
  
    ! get maximum # of aqueous species in a aqueous complexation reaction
    cur_sec_aq_spec => reaction%secondary_species_list
    max_aq_species = 0
    do
      if (.not.associated(cur_sec_aq_spec)) exit
      max_aq_species = max(cur_sec_aq_spec%dbaserxn%nspec,max_aq_species)
      cur_sec_aq_spec => cur_sec_aq_spec%next
    enddo
    
    allocate(reaction%secondary_species_names(reaction%neqcplx))
    reaction%secondary_species_names = ''

    allocate(reaction%secondary_species_print(reaction%neqcplx))
    reaction%secondary_species_print = PETSC_FALSE

    allocate(reaction%eqcplx_basis_names(max_aq_species,reaction%neqcplx))
    reaction%eqcplx_basis_names = ''

    allocate(reaction%eqcplxspecid(0:max_aq_species,reaction%neqcplx))
    reaction%eqcplxspecid = 0

    allocate(reaction%eqcplxstoich(0:max_aq_species,reaction%neqcplx))
    reaction%eqcplxstoich = 0.d0

    allocate(reaction%eqcplxh2oid(reaction%neqcplx))
    reaction%eqcplxh2oid = 0

    allocate(reaction%eqcplxh2ostoich(reaction%neqcplx))
    reaction%eqcplxh2ostoich = 0.d0

    allocate(reaction%eqcplx_logK(reaction%neqcplx))
    reaction%eqcplx_logK = 0.d0

    if (.not.reaction%use_geothermal_hpt) then
      if (option%use_isothermal) then
        allocate(reaction%eqcplx_logKcoef(reaction%num_dbase_temperatures, &
                                          reaction%neqcplx))
      else
        allocate(reaction%eqcplx_logKcoef(FIVE_INTEGER,reaction%neqcplx))
      endif
    else
      allocate(reaction%eqcplx_logKcoef(num_logKs,reaction%neqcplx))
    endif
    
    reaction%eqcplx_logKcoef = 0.d0

    allocate(reaction%eqcplx_Z(reaction%neqcplx))
    reaction%eqcplx_Z = 0.d0

    allocate(reaction%eqcplx_molar_wt(reaction%neqcplx))
    reaction%eqcplx_molar_wt = 0.d0

    allocate(reaction%eqcplx_a0(reaction%neqcplx))
    reaction%eqcplx_a0 = 0.d0

    ! pack in reaction arrays
    cur_sec_aq_spec => reaction%secondary_species_list
    isec_spec = 1
    do
      if (.not.associated(cur_sec_aq_spec)) exit

      reaction%secondary_species_names(isec_spec) = &
        cur_sec_aq_spec%name
      reaction%secondary_species_print(isec_spec) = &
        (cur_sec_aq_spec%print_me .or. reaction%print_all_secondary_species)
      ispec = 0
      do i = 1, cur_sec_aq_spec%dbaserxn%nspec
      
!       print *,'database: ',i,cur_sec_aq_spec%dbaserxn%spec_name(i)
        
        if (cur_sec_aq_spec%dbaserxn%spec_ids(i) /= h2o_id) then
          ispec = ispec + 1
          spec_id = cur_sec_aq_spec%dbaserxn%spec_ids(i)
          if (spec_id > h2o_id) spec_id = spec_id - 1
          reaction%eqcplxspecid(ispec,isec_spec) = spec_id
          reaction%eqcplx_basis_names(ispec,isec_spec) = &
            cur_sec_aq_spec%dbaserxn%spec_name(i)
          reaction%eqcplxstoich(ispec,isec_spec) = &
            cur_sec_aq_spec%dbaserxn%stoich(i)
            
        else ! fill in h2o id and stoich
          reaction%eqcplxh2oid(isec_spec) = h2o_id
          reaction%eqcplxh2ostoich(isec_spec) = &
            cur_sec_aq_spec%dbaserxn%stoich(i)
        endif
      enddo
      reaction%eqcplxspecid(0,isec_spec) = ispec

      if (.not.reaction%use_geothermal_hpt) then
        if (option%use_isothermal) then
          call Interpolate(temp_high,temp_low,option%reference_temperature, &
                      cur_sec_aq_spec%dbaserxn%logK(itemp_high), &
                      cur_sec_aq_spec%dbaserxn%logK(itemp_low), &
                      reaction%eqcplx_logK(isec_spec))
        else
          call ReactionFitLogKCoef(reaction%eqcplx_logKcoef(:,isec_spec), &
                                   cur_sec_aq_spec%dbaserxn%logK, &
                                   reaction%secondary_species_names(isec_spec), &
                                   option,reaction)
          call ReactionInitializeLogK(reaction%eqcplx_logKcoef(:,isec_spec), &
                                      cur_sec_aq_spec%dbaserxn%logK, &
                                      reaction%eqcplx_logK(isec_spec), &
                                      option,reaction)
        endif
      else
        reaction%eqcplx_logKcoef(:,isec_spec) = cur_sec_aq_spec%dbaserxn%logK
        call ReactionInitializeLogK_hpt(reaction%eqcplx_logKcoef(:,isec_spec), &
                                        reaction%eqcplx_logK(isec_spec), &
                                        option,reaction)        
    
      endif

      reaction%eqcplx_Z(isec_spec) = cur_sec_aq_spec%Z
      reaction%eqcplx_molar_wt(isec_spec) = cur_sec_aq_spec%molar_weight
      reaction%eqcplx_a0(isec_spec) = cur_sec_aq_spec%a0
  
      isec_spec = isec_spec + 1
      cur_sec_aq_spec => cur_sec_aq_spec%next
    enddo

  endif
  nullify(cur_sec_aq_spec)
  isec_spec = -1 ! to catch bugs

  ! gas complexes
  ! passive
  call ReactionDatabaseSetupGases(reaction,num_logKs,option,h2o_id, &
                                  temp_high,temp_low,itemp_high,itemp_low, &
                                  reaction%gas%list,PASSIVE_GAS, &
                                  reaction%gas%npassive_gas, &
                                  reaction%gas%passive_names, &
                                  reaction%gas%passive_print_me, &
                                  reaction%gas%paseqspecid, &
                                  reaction%gas%paseqstoich, &
                                  reaction%gas%paseqh2oid, &
                                  reaction%gas%paseqh2ostoich, &
                                  reaction%gas%paseqlogK, &
                                  reaction%gas%paseqlogKcoef)
  ! active
  call ReactionDatabaseSetupGases(reaction,num_logKs,option,h2o_id, &
                                  temp_high,temp_low,itemp_high,itemp_low, &
                                  reaction%gas%list,ACTIVE_GAS, &
                                  reaction%gas%nactive_gas, &
                                  reaction%gas%active_names, &
                                  reaction%gas%active_print_me, &
                                  reaction%gas%acteqspecid, &
                                  reaction%gas%acteqstoich, &
                                  reaction%gas%acteqh2oid, &
                                  reaction%gas%acteqh2ostoich, &
                                  reaction%gas%acteqlogK, &
                                  reaction%gas%acteqlogKcoef)
  
  
  ! immobile species
  immobile%nimmobile = ImmobileGetCount(immobile)
  if (immobile%nimmobile > 0) then
    allocate(immobile%names(immobile%nimmobile))
    immobile%names = ''
    allocate(immobile%print_me(immobile%nimmobile))
    immobile%print_me = PETSC_FALSE

    cur_immobile_spec => immobile%list
    temp_int = 0
    do
      if (.not.associated(cur_immobile_spec)) exit
      temp_int = temp_int + 1
      immobile%names(temp_int) = cur_immobile_spec%name
      immobile%print_me(temp_int) = cur_immobile_spec%print_me .or. &
                                   immobile%print_all
      cur_immobile_spec => cur_immobile_spec%next
    enddo
  endif
  
  ! minerals
  ! Count the number of kinetic mineral reactions, max number of prefactors
  ! in a tst reaction, and the maximum number or species in a prefactor
  temp_int = mineral%nkinmnrl !geh: store for check after processing
  mineral%nkinmnrl = 0
  max_num_prefactors = 0
  max_num_prefactor_species = 0
  cur_mineral => mineral%mineral_list
  !
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%itype == MINERAL_KINETIC .and. &
        associated(cur_mineral%tstrxn)) then
      ! increment number of kinetic minerals
      mineral%nkinmnrl = mineral%nkinmnrl + 1
      cur_prefactor => cur_mineral%tstrxn%prefactor
      ! zero number of prefactors
      i = 0
      do
        if (.not.associated(cur_prefactor)) exit
        i = i + 1
        cur_prefactor_species => cur_prefactor%species
        ! zero number of prefactor species
        j = 0
        do
          if (.not.associated(cur_prefactor_species)) exit
          j = j + 1
          cur_prefactor_species => cur_prefactor_species%next
        enddo
        if (j > max_num_prefactor_species) max_num_prefactor_species = j
        cur_prefactor => cur_prefactor%next
      enddo
      if (i > max_num_prefactors) max_num_prefactors = i
    endif
    cur_mineral => cur_mineral%next
  enddo
  
  if (mineral%nkinmnrl /= temp_int) then
    write(string,'(2i4)') temp_int, mineral%nkinmnrl
    option%io_buffer = 'Inconsistent number of kinetic minerals: ' // &
      trim(string)
    call printErrMsg(option)
  endif

  if (mineral%nmnrl > 0) then
  
    ! get maximum # of aqueous species in a mineral reaction
    cur_mineral => mineral%mineral_list
    max_aq_species = 0
    do
      if (.not.associated(cur_mineral)) exit
      max_aq_species = max(cur_mineral%dbaserxn%nspec,max_aq_species)
      cur_mineral => cur_mineral%next
    enddo
    
    allocate(mineral%mineral_names(mineral%nmnrl))
    mineral%mineral_names = ''
    allocate(mineral%mnrlspecid(0:max_aq_species,mineral%nmnrl))
    mineral%mnrlspecid = 0
    allocate(mineral%mnrlstoich(max_aq_species,mineral%nmnrl))
    mineral%mnrlstoich = 0.d0
    allocate(mineral%mnrlh2oid(mineral%nmnrl))
    mineral%mnrlh2oid = 0
    allocate(mineral%mnrlh2ostoich(mineral%nmnrl))
    mineral%mnrlh2ostoich = 0.d0
    allocate(mineral%mnrl_logK(mineral%nmnrl))
    mineral%mnrl_logK = 0.d0
    allocate(mineral%mnrl_print(mineral%nmnrl))
    mineral%mnrl_print = PETSC_FALSE
    if (.not.reaction%use_geothermal_hpt) then
      if (option%use_isothermal) then
        allocate(mineral%mnrl_logKcoef(reaction%num_dbase_temperatures, &
                                        mineral%nmnrl))
      else
        allocate(mineral%mnrl_logKcoef(FIVE_INTEGER,mineral%nmnrl))
      endif
    else
      allocate(mineral%mnrl_logKcoef(num_logKs,mineral%nmnrl))
    endif
    
    reaction%mineral%mnrl_logKcoef = 0.d0

    if (mineral%nkinmnrl > 0) then
    
      ! get maximum # of aqueous species in a mineral reaction
      cur_mineral => mineral%mineral_list
      max_aq_species = 0
      do
        if (.not.associated(cur_mineral)) exit
        if (associated(cur_mineral%tstrxn)) then ! reaction is kinetic
          max_aq_species = max(cur_mineral%dbaserxn%nspec,max_aq_species)
        endif
        cur_mineral => cur_mineral%next
      enddo
    
      allocate(mineral%kinmnrl_names(mineral%nkinmnrl))
      mineral%kinmnrl_names = ''
      allocate(mineral%kinmnrl_print(mineral%nkinmnrl))
      mineral%kinmnrl_print = PETSC_FALSE
      allocate(mineral%kinmnrlspecid(0:max_aq_species,mineral%nkinmnrl))
      mineral%kinmnrlspecid = 0
      allocate(mineral%kinmnrlstoich(max_aq_species,mineral%nkinmnrl))
      mineral%kinmnrlstoich = 0.d0
      allocate(mineral%kinmnrlh2oid(mineral%nkinmnrl))
      mineral%kinmnrlh2oid = 0
      allocate(mineral%kinmnrlh2ostoich(mineral%nkinmnrl))
      mineral%kinmnrlh2ostoich = 0.d0
      allocate(mineral%kinmnrl_logK(mineral%nkinmnrl))
      mineral%kinmnrl_logK = 0.d0
      if (.not.reaction%use_geothermal_hpt) then
        if (option%use_isothermal) then
          allocate(mineral%kinmnrl_logKcoef(reaction%num_dbase_temperatures, &
                                             mineral%nkinmnrl))
        else
          allocate(mineral%kinmnrl_logKcoef(FIVE_INTEGER,mineral%nkinmnrl))
        endif
      else
        allocate(mineral%kinmnrl_logKcoef(num_logKs,mineral%nkinmnrl))
      endif
      
      mineral%kinmnrl_logKcoef = 0.d0

      ! TST Rxn variables
      allocate(mineral%kinmnrl_affinity_threshold(mineral%nkinmnrl))
      mineral%kinmnrl_affinity_threshold = 0.d0
      allocate(mineral%kinmnrl_rate_limiter(mineral%nkinmnrl))
      mineral%kinmnrl_rate_limiter = 0.d0
      allocate(mineral%kinmnrl_irreversible(mineral%nkinmnrl))
      mineral%kinmnrl_irreversible = 0
      allocate(mineral%kinmnrl_rate_constant(mineral%nkinmnrl))
      mineral%kinmnrl_rate_constant = 0.d0
      allocate(mineral%kinmnrl_activation_energy(mineral%nkinmnrl))
      mineral%kinmnrl_activation_energy = 0.d0
      allocate(mineral%kinmnrl_molar_vol(mineral%nkinmnrl))
      mineral%kinmnrl_molar_vol = 0.d0
      allocate(mineral%kinmnrl_molar_wt(mineral%nkinmnrl))
      mineral%kinmnrl_molar_wt = 0.d0

      allocate(mineral%kinmnrl_armor_pwr(mineral%nkinmnrl))
      mineral%kinmnrl_armor_pwr = 0.d0

      allocate(mineral%kinmnrl_armor_crit_vol_frac(mineral%nkinmnrl))
      mineral%kinmnrl_armor_crit_vol_frac = 0.d0

      allocate(mineral%kinmnrl_armor_min_names(mineral%nkinmnrl))
      mineral%kinmnrl_armor_min_names = ''

      allocate(mineral%kinmnrl_num_prefactors(mineral%nkinmnrl))
      mineral%kinmnrl_num_prefactors = 0
      if (max_num_prefactors > 0) then
        allocate(mineral%kinmnrl_pref_rate(max_num_prefactors,mineral%nkinmnrl))
        mineral%kinmnrl_pref_rate = 0.d0
        allocate(mineral%kinmnrl_pref_activation_energy(max_num_prefactors, &
                                                         mineral%nkinmnrl))
        mineral%kinmnrl_pref_activation_energy = 0.d0
        allocate(mineral%kinmnrl_prefactor_id(0:max_num_prefactor_species, &
                                           max_num_prefactors,mineral%nkinmnrl))
        mineral%kinmnrl_prefactor_id = 0
        allocate(mineral%kinmnrl_pref_alpha(max_num_prefactor_species, &
                                           max_num_prefactors,mineral%nkinmnrl))
        mineral%kinmnrl_pref_alpha = 0.d0
        allocate(mineral%kinmnrl_pref_beta(max_num_prefactor_species, &
                                           max_num_prefactors,mineral%nkinmnrl))
        mineral%kinmnrl_pref_beta = 0.d0
        allocate(mineral%kinmnrl_pref_atten_coef(max_num_prefactor_species, &
                                           max_num_prefactors,mineral%nkinmnrl))
        mineral%kinmnrl_pref_atten_coef = 0.d0
      endif
    endif
    
    ! Determine whether mineral scale factor is used in any TST reactions
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (Initialized(cur_mineral%tstrxn%min_scale_factor)) then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (found) then
      allocate(mineral%kinmnrl_min_scale_factor(mineral%nkinmnrl))
      mineral%kinmnrl_min_scale_factor = 1.d0
    endif

    ! Determine whether Temkin's constant is used in any TST reactions
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (Initialized(cur_mineral%tstrxn%affinity_factor_sigma)) then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (found) then
      allocate(mineral%kinmnrl_Temkin_const(mineral%nkinmnrl))
      mineral%kinmnrl_Temkin_const = 1.d0
    endif

    ! Determine whether affinity factor has power
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (Initialized(cur_mineral%tstrxn%affinity_factor_beta)) then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (found) then
      allocate(mineral%kinmnrl_affinity_power(mineral%nkinmnrl))
      mineral%kinmnrl_affinity_power = 1.d0    
    endif

    ! Determine whether surface area volume fraction power defined
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (.not.Equal(cur_mineral%tstrxn%surf_area_vol_frac_pwr, &
                       0.d0)) then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (reaction%update_mineral_surface_area .or. found) then
      allocate(mineral%kinmnrl_surf_area_vol_frac_pwr(mineral%nkinmnrl))
      mineral%kinmnrl_surf_area_vol_frac_pwr = 0.d0    
    endif
    
    ! Determine whether surface area volume fraction power defined
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (.not.Equal(cur_mineral%tstrxn%surf_area_porosity_pwr, &
                       0.d0)) then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (found) then
      allocate(mineral%kinmnrl_surf_area_porosity_pwr(mineral%nkinmnrl))
      mineral%kinmnrl_surf_area_porosity_pwr = 0.d0    
    endif

#if 0
    ! Determine whether armor mineral name defined
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (.not. cur_mineral%tstrxn%armor_min_name == '') then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (found) then
      allocate(mineral%kinmnrl_armor_min_names(mineral%nkinmnrl))
      mineral%kinmnrl_armor_min_names = ''
    endif

    ! Determine whether armor mineral volume fraction power defined
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (.not.Equal(cur_mineral%tstrxn%armor_pwr,0.d0)) then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (found) then
      allocate(mineral%kinmnrl_armor_pwr(mineral%nkinmnrl))
      mineral%kinmnrl_armor_pwr = 0.d0
    endif
    
    ! Determine whether armor critical volume fraction defined
    cur_mineral => mineral%mineral_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_mineral)) exit
      if (associated(cur_mineral%tstrxn)) then 
        if (.not.Equal(cur_mineral%tstrxn%armor_crit_vol_frac, &
                       0.d0)) then
          found = PETSC_TRUE
          exit
        endif
      endif
      cur_mineral => cur_mineral%next
    enddo
    if (found) then
      allocate(mineral%kinmnrl_armor_crit_vol_frac(mineral%nkinmnrl))
      mineral%kinmnrl_armor_crit_vol_frac = 0.d0
    endif
#endif

    cur_mineral => mineral%mineral_list
    imnrl = 1
    ikinmnrl = 1
    do
      if (.not.associated(cur_mineral)) exit

      mineral%mineral_names(imnrl) = cur_mineral%name
      ispec = 0
      do i = 1, cur_mineral%dbaserxn%nspec
        if (cur_mineral%dbaserxn%spec_ids(i) /= h2o_id) then
          ispec = ispec + 1
          spec_id = cur_mineral%dbaserxn%spec_ids(i)
          if (spec_id > h2o_id) spec_id = spec_id - 1
          mineral%mnrlspecid(ispec,imnrl) = spec_id
          mineral%mnrlstoich(ispec,imnrl) = &
            cur_mineral%dbaserxn%stoich(i)
            
        else ! fill in h2o id and stoich
          mineral%mnrlh2oid(imnrl) = h2o_id
          mineral%mnrlh2ostoich(imnrl) = &
            cur_mineral%dbaserxn%stoich(i)
        endif
      enddo
      mineral%mnrlspecid(0,imnrl) = ispec

      if (.not.reaction%use_geothermal_hpt) then
        if (option%use_isothermal) then
          call Interpolate(temp_high,temp_low,option%reference_temperature, &
                           cur_mineral%dbaserxn%logK(itemp_high), &
                           cur_mineral%dbaserxn%logK(itemp_low), &
                           mineral%mnrl_logK(imnrl))
        else
          call ReactionFitLogKCoef(mineral%mnrl_logKcoef(:,imnrl), &
                                   cur_mineral%dbaserxn%logK, &
                                   mineral%mineral_names(imnrl), &
                                   option,reaction)
          call ReactionInitializeLogK(mineral%mnrl_logKcoef(:,imnrl), &
                                      cur_mineral%dbaserxn%logK, &
                                      mineral%mnrl_logK(imnrl), &
                                      option,reaction)
        endif
      else
        mineral%mnrl_logKcoef(:,imnrl) = cur_mineral%dbaserxn%logK
        call ReactionInitializeLogK_hpt(mineral%mnrl_logKcoef(:,imnrl), &
                                        mineral%mnrl_logK(imnrl), &
                                        option,reaction)      
      endif

      ! geh - for now, the user must specify they want each individual
      !       mineral printed for non-kinetic reactions (e.g. for SI).
      mineral%mnrl_print(imnrl) = cur_mineral%print_me .or. &
                                  reaction%mineral%print_all
      if (cur_mineral%itype == MINERAL_KINETIC) then
        mineral%kinmnrl_names(ikinmnrl) = mineral%mineral_names(imnrl)
        mineral%kinmnrl_print(ikinmnrl) = cur_mineral%print_me .or. &
                                           reaction%mineral%print_all
        mineral%kinmnrlspecid(:,ikinmnrl) = mineral%mnrlspecid(:,imnrl)
        mineral%kinmnrlstoich(:,ikinmnrl) = mineral%mnrlstoich(:,imnrl)
        mineral%kinmnrlh2oid(ikinmnrl) = mineral%mnrlh2oid(imnrl)
        mineral%kinmnrlh2ostoich(ikinmnrl) = mineral%mnrlh2ostoich(imnrl)

        if (.not.reaction%use_geothermal_hpt) then
          if (option%use_isothermal) then
            call Interpolate(temp_high,temp_low,option%reference_temperature, &
                             cur_mineral%dbaserxn%logK(itemp_high), &
                             cur_mineral%dbaserxn%logK(itemp_low), &
                             mineral%kinmnrl_logK(ikinmnrl))
          else
            call ReactionFitLogKCoef(mineral%kinmnrl_logKcoef(:,ikinmnrl), &
                                     cur_mineral%dbaserxn%logK, &
                                     mineral%kinmnrl_names(ikinmnrl), &
                                     option,reaction)
            call ReactionInitializeLogK(mineral%kinmnrl_logKcoef(:,ikinmnrl), &
                                        cur_mineral%dbaserxn%logK, &
                                        mineral%kinmnrl_logK(ikinmnrl), &
                                        option,reaction)
          endif
        else
          mineral%kinmnrl_logKcoef(:,ikinmnrl) = cur_mineral%dbaserxn%logK
          call ReactionInitializeLogK_hpt(mineral%kinmnrl_logKcoef(:,ikinmnrl), &
                                          mineral%kinmnrl_logK(ikinmnrl), &
                                          option,reaction)        
        endif

        tstrxn => cur_mineral%tstrxn
        if (associated(tstrxn)) then
          ! loop over transition state theory reactions/prefactors
          cur_prefactor => cur_mineral%tstrxn%prefactor
          i = 0
          do
            if (.not.associated(cur_prefactor)) exit
            ! ith prefactor
            i = i + 1

            mineral%kinmnrl_pref_rate(i,ikinmnrl) = cur_prefactor%rate
            mineral%kinmnrl_pref_activation_energy(i,ikinmnrl) = &
              cur_prefactor%activation_energy

            cur_prefactor_species => cur_prefactor%species
            j = 0
            do
              if (.not.associated(cur_prefactor_species)) exit
              ! jth prefactor species
              j = j + 1
              ! find the prefactor species
              do ispec = 1, reaction%naqcomp
                if (StringCompareIgnoreCase( &
                                    reaction%primary_species_names(ispec), &
                                    cur_prefactor_species%name)) then
                  cur_prefactor_species%id = ispec
                  exit
                endif
              enddo
              if (cur_prefactor_species%id == 0) then ! not found
                ! negative prefactor_species_id denotes a secondary species
                do ispec = 1, reaction%neqcplx
                  if (StringCompareIgnoreCase( &
                                   reaction%secondary_species_names(ispec), &
                                   cur_prefactor_species%name)) then
                    cur_prefactor_species%id = -ispec
                    exit
                  endif
                enddo
              endif
              if (cur_prefactor_species%id == 0) then
                option%io_buffer = 'Kinetic mineral prefactor species "' // &
                  trim(cur_prefactor_species%name) // &
                  '" not found among primary or secondary species.'
                call printErrMsg(option)
              endif
              mineral%kinmnrl_prefactor_id(j,i,ikinmnrl) = &
                cur_prefactor_species%id
              mineral%kinmnrl_pref_alpha(j,i,ikinmnrl) = &
                cur_prefactor_species%alpha
              mineral%kinmnrl_pref_beta(j,i,ikinmnrl) = &
                cur_prefactor_species%beta
              mineral%kinmnrl_pref_atten_coef(j,i,ikinmnrl) = &
                cur_prefactor_species%attenuation_coef
              cur_prefactor_species => cur_prefactor_species%next
            enddo
            ! store the number of species
            mineral%kinmnrl_prefactor_id(0,i,ikinmnrl) = j
            cur_prefactor => cur_prefactor%next
          enddo
          mineral%kinmnrl_num_prefactors(ikinmnrl) = i

          mineral%kinmnrl_affinity_threshold(ikinmnrl) = &
            tstrxn%affinity_threshold
          mineral%kinmnrl_rate_limiter(ikinmnrl) = tstrxn%rate_limiter
          mineral%kinmnrl_irreversible(ikinmnrl) = tstrxn%irreversible

          mineral%kinmnrl_armor_min_names(ikinmnrl) = tstrxn%armor_min_name
          mineral%kinmnrl_armor_pwr(ikinmnrl) = tstrxn%armor_pwr
          mineral%kinmnrl_armor_crit_vol_frac(ikinmnrl) = &
            tstrxn%armor_crit_vol_frac

          if (mineral%kinmnrl_num_prefactors(ikinmnrl) == 0) then
            ! no prefactors, rates stored in upper level
            mineral%kinmnrl_rate_constant(ikinmnrl) = tstrxn%rate
            mineral%kinmnrl_activation_energy(ikinmnrl) = &
              tstrxn%activation_energy
          endif
          if (Initialized(tstrxn%min_scale_factor)) then
            mineral%kinmnrl_min_scale_factor(ikinmnrl) = &
              tstrxn%min_scale_factor
          endif
          if (Initialized(tstrxn%affinity_factor_sigma)) then
            mineral%kinmnrl_Temkin_const(ikinmnrl) = &
              tstrxn%affinity_factor_sigma
          endif
          if (Initialized(tstrxn%affinity_factor_beta)) then
            mineral%kinmnrl_affinity_power(ikinmnrl) = &
              tstrxn%affinity_factor_beta
          endif
          if (.not.Equal(tstrxn%surf_area_vol_frac_pwr, &
                         0.d0)) then
            mineral%kinmnrl_surf_area_vol_frac_pwr(ikinmnrl) = &
              tstrxn%surf_area_vol_frac_pwr
          endif
          if (.not.Equal(tstrxn%surf_area_porosity_pwr, &
                         0.d0)) then
            mineral%kinmnrl_surf_area_porosity_pwr(ikinmnrl) = &
              tstrxn%surf_area_porosity_pwr
          endif
        endif ! associated(tstrxn)

        mineral%kinmnrl_molar_vol(ikinmnrl) = cur_mineral%molar_volume
        mineral%kinmnrl_molar_wt(ikinmnrl) = cur_mineral%molar_weight
        ikinmnrl = ikinmnrl + 1
      endif

      cur_mineral => cur_mineral%next
      imnrl = imnrl + 1
    enddo

  endif
  
  ! colloids
  ! already calculated above
  !reaction%ncoll = GetColloidCount(reaction)

  if (reaction%ncoll > 0) then
    allocate(reaction%colloid_names(reaction%ncoll))
    allocate(reaction%colloid_mobile_fraction(reaction%ncoll))
    allocate(reaction%colloid_print(reaction%ncoll))
    reaction%colloid_names = ''
    reaction%colloid_mobile_fraction = 0.d0
    reaction%colloid_print = PETSC_FALSE

    cur_colloid => reaction%colloid_list
    icoll = 1
    do
      if (.not.associated(cur_colloid)) exit

      reaction%colloid_names(icoll) = cur_colloid%name
      reaction%colloid_mobile_fraction(icoll) = cur_colloid%mobile_fraction
      reaction%colloid_print(icoll) = cur_colloid%print_me .or. &
                                      reaction%print_all_species
      cur_colloid => cur_colloid%next
      icoll = icoll + 1
    enddo
  endif
  
  ! use flags to determine whether a primary aqueous species is included
  ! in the list of colloid species
  allocate(colloid_species_flag(reaction%naqcomp))
  colloid_species_flag = PETSC_FALSE

  ! allocate colloids species names, mappings, etc.
  reaction%ncollcomp = 0
  icount = 0
  do i = 1, reaction%naqcomp
    if (colloid_species_flag(i)) then 
      icount = icount + 1
    endif
  enddo
  if (icount > 0) then
    allocate(reaction%pri_spec_to_coll_spec(reaction%naqcomp))
    allocate(reaction%colloid_species_names(icount))
    allocate(reaction%coll_spec_to_pri_spec(icount))
    reaction%pri_spec_to_coll_spec = UNINITIALIZED_INTEGER
    reaction%coll_spec_to_pri_spec = UNINITIALIZED_INTEGER
    reaction%colloid_species_names = ''
    reaction%ncollcomp = icount
    icount = 0
    do i = 1, reaction%naqcomp
      if (colloid_species_flag(i)) then
        icount = icount + 1
        reaction%colloid_species_names(icount) = &
          trim(reaction%primary_species_names(i))
        reaction%coll_spec_to_pri_spec(icount) = i
        reaction%pri_spec_to_coll_spec(i) = icount
      endif
    enddo
    if (minval(reaction%coll_spec_to_pri_spec) < 1) then
      option%io_buffer = 'Species colloid surface complexation reaction not &
                         & recognized among primary species'
      call printErrMsg(option)
    endif
    allocate(reaction%total_sorb_mobile_print(reaction%ncollcomp))
    reaction%total_sorb_mobile_print = PETSC_FALSE
    do i = 1, reaction%ncollcomp
      reaction%total_sorb_mobile_print(i) = &
        (reaction%primary_species_print(reaction%coll_spec_to_pri_spec(i)) .or. &
         reaction%print_all_species) .and. &
        reaction%print_total_sorb_mobile
    enddo
  endif
  deallocate(colloid_species_flag)

  if (reaction%neqionxrxn > 0) then

    ! determine max # cations for a given ionx exchange rxn
    icount = 0
    cur_ionx_rxn => reaction%ion_exchange_rxn_list
    do
      if (.not.associated(cur_ionx_rxn)) exit
      ication = 0
      cur_cation => cur_ionx_rxn%cation_list
      do
        if (.not.associated(cur_cation)) exit
        ication = ication + 1
        cur_cation => cur_cation%next
      enddo
      if (ication > icount) icount = ication
      cur_ionx_rxn => cur_ionx_rxn%next
    enddo
    nullify(cur_ionx_rxn)
    
    allocate(reaction%eqionx_rxn_cationid(0:icount,reaction%neqionxrxn))
    reaction%eqionx_rxn_cationid = 0
    allocate(reaction%eqionx_rxn_Z_flag(reaction%neqionxrxn))
    reaction%eqionx_rxn_Z_flag = PETSC_FALSE
    allocate(reaction%eqionx_rxn_cation_X_offset(reaction%neqionxrxn))
    reaction%eqionx_rxn_cation_X_offset = 0
    allocate(reaction%eqionx_rxn_CEC(reaction%neqionxrxn))
    reaction%eqionx_rxn_CEC = 0.d0
    allocate(reaction%eqionx_rxn_to_surf(reaction%neqionxrxn))
    reaction%eqionx_rxn_to_surf = 0
    allocate(reaction%eqionx_rxn_k(icount,reaction%neqionxrxn))
    reaction%eqionx_rxn_k = 0.d0

    irxn = 0
    icount = 0
    cur_ionx_rxn => reaction%ion_exchange_rxn_list
    do
      if (.not.associated(cur_ionx_rxn)) exit
      irxn = irxn + 1
      ication = 0
      reaction%eqionx_rxn_CEC(irxn) = cur_ionx_rxn%CEC
        ! compute the offset to the first cation in rxn
      reaction%eqionx_rxn_cation_X_offset(irxn) = icount
      if (len_trim(cur_ionx_rxn%mineral_name) > 1) then
        reaction%eqionx_rxn_to_surf(irxn) = &
          GetKineticMineralIDFromName(cur_ionx_rxn%mineral_name, &
                                      reaction%mineral,option)
        if (reaction%eqionx_rxn_to_surf(irxn) < 0) then
          option%io_buffer = 'Mineral ' // trim(cur_ionx_rxn%mineral_name) // &
            ' listed in ion exchange &reaction not found in mineral list'
          call printErrMsg(option)
        endif
      endif
      cur_cation => cur_ionx_rxn%cation_list
      do
        if (.not.associated(cur_cation)) exit
        ication = ication + 1
        icount = icount + 1
        reaction%eqionx_rxn_k(ication,irxn) = cur_cation%k

        found = PETSC_FALSE
        do i = 1, reaction%naqcomp
          if (StringCompare(cur_cation%name, &
                              reaction%primary_species_names(i), &
                              MAXWORDLENGTH)) then
            reaction%eqionx_rxn_cationid(ication,irxn) = i
            found = PETSC_TRUE        
          endif
        enddo
        if (.not.found) then
          option%io_buffer = 'Cation ' // trim(cur_cation%name) // &
                   ' in ion exchange reaction not found in swapped basis.'
          call printErrMsg(option)  
        endif
        cur_cation => cur_cation%next
      enddo
      reaction%eqionx_rxn_cationid(0,irxn) = ication
      ! Find any Zi /= Zj for all species i, j
      found = PETSC_FALSE
      do i = 1, reaction%eqionx_rxn_cationid(0,irxn)
        do j = 1, reaction%eqionx_rxn_cationid(0,irxn)
          if (abs( &
              reaction%primary_spec_Z(reaction%eqionx_rxn_cationid(i,irxn))- &
              reaction%primary_spec_Z(reaction%eqionx_rxn_cationid(j,irxn))) > &
              0.1d0) then
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (found) exit
      enddo
      reaction%eqionx_rxn_Z_flag(irxn) = found
      cur_ionx_rxn => cur_ionx_rxn%next
    enddo
    nullify(cur_ionx_rxn)

  endif

  ! radioactive decay reaction
  
  if (reaction%nradiodecay_rxn > 0) then
  
    ! process reaction equation into the database format
    cur_radiodecay_rxn => reaction%radioactive_decay_rxn_list
    do
      if (.not.associated(cur_radiodecay_rxn)) exit
      cur_radiodecay_rxn%dbaserxn => &
        DatabaseRxnCreateFromRxnString(cur_radiodecay_rxn%reaction, &
                                       reaction%naqcomp, &
                                       reaction%offset_aqueous, &
                                       reaction%primary_species_names, &
                                       reaction%nimcomp, &
                                       reaction%offset_immobile, &
                                       reaction%immobile%names, &
                                       PETSC_FALSE,option)
      cur_radiodecay_rxn => cur_radiodecay_rxn%next
    enddo
    nullify(cur_radiodecay_rxn)

    ! determine max # species for a given radiodecay rxn
    max_species_count = 0
    cur_radiodecay_rxn => reaction%radioactive_decay_rxn_list
    do
      if (.not.associated(cur_radiodecay_rxn)) exit

      ! zero count
      forward_count = 0

      ! max species in reaction
      species_count = cur_radiodecay_rxn%dbaserxn%nspec

      ! sum forward and reverse species
      dbaserxn => cur_radiodecay_rxn%dbaserxn
      do i = 1, dbaserxn%nspec
        if (dbaserxn%stoich(i) < 0.d0) then
          forward_count = forward_count + 1
        endif
      enddo
      
      if (forward_count > 1) then ! currently cannot have more than one species
        option%io_buffer = 'Cannot have more than one reactant in &
                           &radioactive decay reaction: (' // &
                           trim(cur_radiodecay_rxn%reaction) // ').'
        call printErrMsg(option)
      endif

      ! calculate maximum
      if (species_count > max_species_count) max_species_count = species_count

      cur_radiodecay_rxn => cur_radiodecay_rxn%next

    enddo
    nullify(cur_radiodecay_rxn)
    
    allocate(reaction%radiodecayspecid(0:max_species_count, &
                                       reaction%nradiodecay_rxn))
    reaction%radiodecayspecid = 0
    allocate(reaction%radiodecaystoich(max_species_count, &
                                       reaction%nradiodecay_rxn))
    reaction%radiodecaystoich = 0.d0
    allocate(reaction%radiodecayforwardspecid(reaction%nradiodecay_rxn))
    reaction%radiodecayforwardspecid = 0
    allocate(reaction%radiodecay_kf(reaction%nradiodecay_rxn))
    reaction%radiodecay_kf = 0.d0

    ! load the data into the compressed arrays
    irxn = 0
    cur_radiodecay_rxn => reaction%radioactive_decay_rxn_list
    do
      if (.not.associated(cur_radiodecay_rxn)) exit
      
      dbaserxn => cur_radiodecay_rxn%dbaserxn
      
      irxn = irxn + 1
     
      forward_count = 0
      backward_count = 0
      do i = 1, dbaserxn%nspec
        reaction%radiodecayspecid(i,irxn) = dbaserxn%spec_ids(i)
        reaction%radiodecaystoich(i,irxn) = dbaserxn%stoich(i)
        if (dbaserxn%stoich(i) < 0.d0) then
          reaction%radiodecayforwardspecid(irxn) = dbaserxn%spec_ids(i)
        endif
      enddo
      reaction%radiodecayspecid(0,irxn) = dbaserxn%nspec
      reaction%radiodecay_kf(irxn) = cur_radiodecay_rxn%rate_constant
      
      cur_radiodecay_rxn => cur_radiodecay_rxn%next
      
    enddo
              
  endif 
  
  ! general reaction
  
  if (reaction%ngeneral_rxn > 0) then
  
    ! process reaction equation into the database format
    cur_general_rxn => reaction%general_rxn_list
    do
      if (.not.associated(cur_general_rxn)) exit
      cur_general_rxn%dbaserxn => &
        DatabaseRxnCreateFromRxnString(cur_general_rxn%reaction, &
                                       reaction%naqcomp, &
                                       reaction%offset_aqueous, &
                                       reaction%primary_species_names, &
                                       reaction%nimcomp, &
                                       reaction%offset_immobile, &
                                       reaction%immobile%names, &
                                       PETSC_FALSE,option)
      cur_general_rxn => cur_general_rxn%next
    enddo
    nullify(cur_general_rxn)

    ! determine max # species, forward species and backward species
    !  for a given general rxn
    max_species_count = 0
    max_forward_count = 0
    max_backward_count = 0
    cur_general_rxn => reaction%general_rxn_list
    do
      if (.not.associated(cur_general_rxn)) exit

      ! zero count
      forward_count = 0
      backward_count = 0   

      ! max species in reaction
      species_count = cur_general_rxn%dbaserxn%nspec

      ! sum forward and reverse species
      dbaserxn => cur_general_rxn%dbaserxn
      do i = 1, dbaserxn%nspec
        if (dbaserxn%stoich(i) < 0.d0) then
          forward_count = forward_count + 1
        else if (dbaserxn%stoich(i) > 0.d0) then
          backward_count = backward_count + 1
        endif
      enddo

      ! calculate maximum
      if (forward_count > max_forward_count) max_forward_count = forward_count
      if (backward_count > max_backward_count) &
        max_backward_count = backward_count
      if (species_count > max_species_count) max_species_count = species_count

      cur_general_rxn => cur_general_rxn%next

    enddo
    nullify(cur_general_rxn)
    
    allocate(reaction%generalspecid(0:max_species_count,reaction%ngeneral_rxn))
    reaction%generalspecid = 0
    allocate(reaction%generalstoich(max_species_count,reaction%ngeneral_rxn))
    reaction%generalstoich = 0.d0
    allocate(reaction%generalforwardspecid(0:max_forward_count, &
                                           reaction%ngeneral_rxn))
    reaction%generalforwardspecid = 0
    allocate(reaction%generalforwardstoich(max_forward_count, &
                                           reaction%ngeneral_rxn))
    reaction%generalforwardstoich = 0.d0
    allocate(reaction%generalbackwardspecid(0:max_backward_count, &
                                            reaction%ngeneral_rxn))
    reaction%generalbackwardspecid = 0
    allocate(reaction%generalbackwardstoich(max_backward_count, &
                                            reaction%ngeneral_rxn))
    reaction%generalbackwardstoich = 0.d0
    allocate(reaction%generalh2oid(reaction%ngeneral_rxn))
    reaction%generalh2oid = 0
    allocate(reaction%generalh2ostoich(reaction%ngeneral_rxn))
    reaction%generalh2ostoich = 0.d0
    allocate(reaction%general_kf(reaction%ngeneral_rxn))
    reaction%general_kf = 0.d0
    allocate(reaction%general_kr(reaction%ngeneral_rxn))    
    reaction%general_kr = 0.d0

    ! load the data into the compressed arrays
    irxn = 0
    cur_general_rxn => reaction%general_rxn_list
    do
      if (.not.associated(cur_general_rxn)) exit
      
      dbaserxn => cur_general_rxn%dbaserxn
      
      irxn = irxn + 1
     
      forward_count = 0
      backward_count = 0
      do i = 1, dbaserxn%nspec
        reaction%generalspecid(i,irxn) = dbaserxn%spec_ids(i)
        reaction%generalstoich(i,irxn) = dbaserxn%stoich(i)
        if (dbaserxn%stoich(i) < 0.d0) then
          forward_count = forward_count + 1
          reaction%generalforwardspecid(forward_count,irxn) = &
            dbaserxn%spec_ids(i)
          ! ensure that forward stoich is positive for rate expression
          reaction%generalforwardstoich(forward_count,irxn) = &
            dabs(dbaserxn%stoich(i))
        else if (dbaserxn%stoich(i) > 0.d0) then
          backward_count = backward_count + 1
          reaction%generalbackwardspecid(backward_count,irxn) = &
            dbaserxn%spec_ids(i)
          reaction%generalbackwardstoich(backward_count,irxn) = &
            dbaserxn%stoich(i)
        endif
      enddo
      reaction%generalspecid(0,irxn) = dbaserxn%nspec
      reaction%generalforwardspecid(0,irxn) = forward_count
      reaction%generalbackwardspecid(0,irxn) = backward_count
      
      reaction%general_kf(irxn) = cur_general_rxn%forward_rate
      reaction%general_kr(irxn) = cur_general_rxn%backward_rate
      
      cur_general_rxn => cur_general_rxn%next
      
    enddo
              
  endif 
  
  ! immobile decay reaction
  
  if (reaction%immobile%ndecay_rxn > 0) then
  
    allocate(reaction%immobile%decayspecid(reaction%immobile%ndecay_rxn))
    allocate(reaction%immobile%decay_rate_constant(reaction%immobile%ndecay_rxn))
  
    cur_immobile_decay_rxn => reaction%immobile%decay_rxn_list
    irxn = 0
    do
      if (.not.associated(cur_immobile_decay_rxn)) exit
      
      irxn = irxn + 1

      found = PETSC_FALSE
      do i = 1, reaction%immobile%nimmobile
        if (StringCompare(cur_immobile_decay_rxn%species_name, &
                          reaction%immobile%names(i), &
                          MAXWORDLENGTH)) then
          reaction%immobile%decayspecid(irxn) = i
          found = PETSC_TRUE
          exit      
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'Species "' // &
        trim(cur_immobile_decay_rxn%species_name) // &
        '" in immobile decay reaction not found among immobile species.'
        call printErrMsg(option)
      endif
      reaction%immobile%decay_rate_constant(irxn) = &
        cur_immobile_decay_rxn%rate_constant
      cur_immobile_decay_rxn => cur_immobile_decay_rxn%next
    enddo
    nullify(cur_immobile_decay_rxn)

  endif 
  
  ! Kd reactions
  
  if (reaction%neqkdrxn > 0) then

    if (reaction%neqcplx > 0) then
      option%io_buffer = 'Isotherm reactions currently calculated as a &
                         &function of free-ion, not totals.  Contact Glenn!'
      call printErrMsg(option)
    endif
  
    ! allocate arrays
    allocate(reaction%eqkdspecid(reaction%neqkdrxn))
    reaction%eqkdspecid = 0
    allocate(reaction%eqkdtype(reaction%neqkdrxn))
    reaction%eqkdtype = 0
    allocate(reaction%eqkddistcoef(reaction%neqkdrxn))
    reaction%eqkddistcoef = 0.d0
    allocate(reaction%eqkdlangmuirb(reaction%neqkdrxn))
    reaction%eqkdlangmuirb = 0.d0
    allocate(reaction%eqkdfreundlichn(reaction%neqkdrxn))
    reaction%eqkdfreundlichn = 0.d0
    allocate(reaction%eqkdmineral(reaction%neqkdrxn))
    reaction%eqkdmineral = 0

    cur_kd_rxn => reaction%kd_rxn_list
    
    if (option%use_mc) then
      allocate(reaction%sec_cont_eqkdtype(reaction%neqkdrxn))
      reaction%sec_cont_eqkdtype = 0   
      allocate(reaction%sec_cont_eqkddistcoef(reaction%neqkdrxn))
      reaction%sec_cont_eqkddistcoef = 0.d0
      allocate(reaction%sec_cont_eqkdlangmuirb(reaction%neqkdrxn))
      reaction%sec_cont_eqkdlangmuirb = 0.d0
      allocate(reaction%sec_cont_eqkdfreundlichn(reaction%neqkdrxn))
      reaction%sec_cont_eqkdfreundlichn = 0.d0
      sec_cont_cur_kd_rxn => reaction%sec_cont_kd_rxn_list
    endif
    
    irxn = 0
    do  
      if (.not.associated(cur_kd_rxn)) exit

      irxn = irxn + 1

      found = PETSC_FALSE
      do i = 1, reaction%naqcomp
        if (StringCompare(cur_kd_rxn%species_name, &
                          reaction%primary_species_names(i), &
                          MAXWORDLENGTH)) then
          reaction%eqkdspecid(irxn) = i
          found = PETSC_TRUE
          exit      
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'Species ' // trim(cur_kd_rxn%species_name) // &
                 ' in kd reaction &
                 & not found among primary species list.'
        call printErrMsg(option)     
      endif
      reaction%eqkdtype(irxn) = cur_kd_rxn%itype
      ! associate mineral id
      if (len_trim(cur_kd_rxn%kd_mineral_name) > 1) then
        reaction%eqkdmineral(irxn) = &
          GetKineticMineralIDFromName(cur_kd_rxn%kd_mineral_name, &
                                      reaction%mineral,option)
        if (reaction%eqkdmineral(irxn) < 0) then
          option%io_buffer = 'Mineral ' // trim(cur_ionx_rxn%mineral_name) // &
                             ' listed in kd (linear sorption) &
                             &reaction not found in mineral list'
          call printErrMsg(option)
        endif
      endif      
      reaction%eqkddistcoef(irxn) = cur_kd_rxn%Kd
      reaction%eqkdlangmuirb(irxn) = cur_kd_rxn%Langmuir_b
      reaction%eqkdfreundlichn(irxn) = cur_kd_rxn%Freundlich_n
       
      cur_kd_rxn => cur_kd_rxn%next
      
      if (option%use_mc) then
        reaction%sec_cont_eqkdtype(irxn) = sec_cont_cur_kd_rxn%itype
        reaction%sec_cont_eqkddistcoef(irxn) = sec_cont_cur_kd_rxn%Kd
        reaction%sec_cont_eqkdlangmuirb(irxn) = sec_cont_cur_kd_rxn%Langmuir_b
        reaction%sec_cont_eqkdfreundlichn(irxn) = &
          sec_cont_cur_kd_rxn%Freundlich_n
        sec_cont_cur_kd_rxn => sec_cont_cur_kd_rxn%next
      endif
      
      
    enddo
  endif

  call BasisPrint(reaction,'Final Basis',option)

  ! locate specific species
  reaction%species_idx => SpeciesIndexCreate()
  do ispec = 1, reaction%naqcomp
    if (reaction%species_idx%h_ion_id == 0) then
      word = 'H+'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%h_ion_id = ispec
      endif
    endif
    if (reaction%species_idx%na_ion_id == 0) then
      word = 'Na+'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%na_ion_id = ispec
      endif
    endif
    if (reaction%species_idx%cl_ion_id == 0) then
      word = 'Cl-'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%cl_ion_id = ispec
      endif
    endif
    if (reaction%species_idx%co2_aq_id == 0) then
      word = 'CO2(aq)'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_aq_id = ispec
      endif
    endif
    if (reaction%species_idx%tracer_aq_id == 0) then
      word = 'Tracer'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%tracer_aq_id = ispec
      endif
    endif
    if (reaction%species_idx%h2o_aq_id == 0) then
      word = 'H2O'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%h2o_aq_id = ispec
      endif
    endif
    if (reaction%species_idx%tracer_age_id == 0) then
      word = 'Tracer_Age'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%tracer_age_id = ispec
        reaction%calculate_tracer_age = PETSC_TRUE
      endif
    endif
    if (reaction%species_idx%water_age_id == 0) then
      word = 'Water_Age'
      if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%water_age_id = ispec
        reaction%calculate_water_age = PETSC_TRUE
      endif
    endif
  enddo
  
  do ispec = 1, reaction%neqcplx
    if (reaction%species_idx%h_ion_id == 0) then
      word = 'H+'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%h_ion_id = -ispec
      endif
    endif
    if (reaction%species_idx%na_ion_id == 0) then
      word = 'Na+'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%na_ion_id = -ispec
      endif
    endif
    if (reaction%species_idx%cl_ion_id == 0) then
      word = 'Cl-'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%cl_ion_id = -ispec
      endif
    endif
    if (reaction%species_idx%co2_aq_id == 0) then
      word = 'CO2(aq)'
      if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_aq_id = -ispec
      endif
    endif
  enddo

  ! these are passive gases used in CONSTRAINTS
  do ispec = 1, reaction%gas%npassive_gas
    if (reaction%species_idx%o2_gas_id == 0) then
      word = 'O2(g)'
      if (StringCompareIgnoreCase(reaction%gas%passive_names(ispec), &
                                  word)) then
        reaction%species_idx%o2_gas_id = ispec
      endif
    endif
    if (reaction%species_idx%co2_gas_id == 0) then
      word = 'CO2(g)'
      if (StringCompareIgnoreCase(reaction%gas%passive_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_gas_id = ispec
      endif
      word = 'CO2(g)*'
      if (StringCompareIgnoreCase(reaction%gas%passive_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_gas_id = ispec
      endif

    endif

  enddo
  
  ! sandbox reactions
  call RSandboxSetup(reaction,option)
  
90 format(80('-'))
100 format(/,2x,i4,2x,a)
110 format(100(/,14x,3(a20,2x)))
120 format(/,a)

  if (OptionPrintToFile(option)) then
    write(option%fid_out,90)
    write(option%fid_out,100) reaction%naqcomp, 'Primary Species'
    write(option%fid_out,110) &
      (reaction%primary_species_names(i),i=1,reaction%naqcomp)
    
    write(option%fid_out,100) reaction%neqcplx, 'Secondary Complex Species'
    write(option%fid_out,110) &
      (reaction%secondary_species_names(i),i=1,reaction%neqcplx)
    
    write(option%fid_out,100) reaction%gas%nactive_gas, 'Active Gas Species'
    write(option%fid_out,110) (reaction%gas%active_names(i),i=1, &
                               reaction%gas%nactive_gas)
    
    write(option%fid_out,100) reaction%gas%npassive_gas, 'Passive Gas Species'
    write(option%fid_out,110) (reaction%gas%passive_names(i),i=1, &
                               reaction%gas%npassive_gas)
    
    write(option%fid_out,100) mineral%nmnrl, 'Reference Minerals'
    write(option%fid_out,110) (mineral%mineral_names(i),i=1,mineral%nmnrl)
    
    write(option%fid_out,100) mineral%nkinmnrl, 'Kinetic Mineral Reactions'
    write(option%fid_out,110) (mineral%kinmnrl_names(i),i=1,mineral%nkinmnrl)

    write(option%fid_out,100) reaction%neqionxrxn, 'Ion Exchange Reactions'
    write(option%fid_out,100) reaction%neqionxcation, 'Ion Exchange Cations'
    write(option%fid_out,90)
  endif
  
  if (allocated(new_basis)) deallocate(new_basis)
  if (allocated(old_basis)) deallocate(old_basis)
  if (allocated(transformation)) deallocate(transformation)
  if (allocated(stoich_prev)) deallocate(stoich_prev)
  if (allocated(stoich_new)) deallocate(stoich_new)
  if (allocated(logKvector)) deallocate(logKvector)
  if (allocated(indices)) deallocate(indices)

  if (allocated(new_basis_names)) deallocate(new_basis_names)
  if (allocated(old_basis_names)) deallocate(old_basis_names)
  
end subroutine BasisInit

! ************************************************************************** !

function GetSpeciesBasisID(reaction,option,ncomp_h2o,reaction_name, &
                           species_name, &
                           pri_names,sec_names,gas_names)
  ! 
  ! Reduces redundant coding above
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/08
  ! 

  use Option_module
  use String_module

  implicit none

  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscInt :: ncomp_h2o
  character(len=MAXWORDLENGTH) :: reaction_name
  character(len=MAXWORDLENGTH) :: species_name
  character(len=MAXWORDLENGTH) :: pri_names(:)
  character(len=MAXWORDLENGTH) :: sec_names(:)
  character(len=MAXWORDLENGTH) :: gas_names(:)

  PetscInt :: GetSpeciesBasisID
  PetscInt :: i

  GetSpeciesBasisID = 0
  do i=1,ncomp_h2o
    if (StringCompare(species_name, &
                        pri_names(i),MAXWORDLENGTH)) then
      GetSpeciesBasisID = i
      return
    endif
  enddo
  ! secondary aqueous and gas species denoted by negative id
  do i=1,reaction%neqcplx
    if (StringCompare(species_name, &
                        sec_names(i),MAXWORDLENGTH)) then
      GetSpeciesBasisID = -i
      return
    endif
  enddo
  do i=1,reaction%gas%ngas
    if (StringCompare(species_name, &
                      gas_names(i),MAXWORDLENGTH)) then
      GetSpeciesBasisID = -(reaction%neqcplx+i)
      return
    endif
  enddo
  
  option%io_buffer = 'Species ' // &
           trim(species_name) // &
           ' listed in reaction for ' // &
           trim(reaction_name) // &
           ' not found among primary, secondary, or gas species.'
  call printErrMsg(option)

end function GetSpeciesBasisID

! ************************************************************************** !

subroutine ReactionDatabaseSetupGases(reaction,num_logKs,option,h2o_id, &
                                      temp_high,temp_low, &
                                      itemp_high,itemp_low, &
                                      gas_species_list,gas_itype, &
                                      ngas,gas_names,gas_print, &
                                      eqspecid,eqstoich,eqh2oid,eqh2ostoich, &
                                      eqlogK,eqlogKcoef)
  ! 
  ! Sets up gas reactions (both active and passive).  Placing setup of both
  ! active and passive gases in a single subroutine removes redundancy
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/10/16
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Reaction_Gas_Aux_module
  use Utility_module
  
  implicit none
  
  type(reaction_type) :: reaction
  PetscInt :: num_logKs
  type(option_type) :: option
  PetscInt :: h2o_id
  PetscReal :: temp_high, temp_low
  PetscInt :: itemp_high, itemp_low
  type(gas_species_type), pointer :: gas_species_list
  PetscInt :: gas_itype
  PetscInt :: ngas
  character(len=MAXWORDLENGTH), pointer :: gas_names(:)
  PetscBool, pointer :: gas_print(:)
  PetscInt, pointer :: eqspecid(:,:)
  PetscReal, pointer :: eqstoich(:,:)
  PetscInt, pointer :: eqh2oid(:)
  PetscReal, pointer :: eqh2ostoich(:)
  PetscReal, pointer :: eqlogK(:)
  PetscReal, pointer :: eqlogKcoef(:,:)
  
  type(gas_species_type), pointer :: cur_gas_spec
  PetscInt :: max_aq_species
  PetscInt :: igas_spec
  PetscInt :: ispec
  PetscInt :: i
  PetscInt :: spec_id
    
  ngas = GasGetCount(gas_species_list,gas_itype)
  if (ngas > 0) then
  
    ! get maximum # of aqueous species in a gas reaction
    cur_gas_spec => gas_species_list
    max_aq_species = 0
    do
      if (.not.associated(cur_gas_spec)) exit
      if (cur_gas_spec%itype == gas_itype .or. &
          cur_gas_spec%itype == ACTIVE_AND_PASSIVE_GAS) then
        max_aq_species = max(cur_gas_spec%dbaserxn%nspec,max_aq_species)
      endif
      cur_gas_spec => cur_gas_spec%next
    enddo
    
    allocate(gas_names(ngas))
    gas_names = ''
    allocate(gas_print(ngas))
    gas_print = PETSC_FALSE
    allocate(eqspecid(0:max_aq_species,ngas))
    eqspecid = 0
    allocate(eqstoich(0:max_aq_species,ngas))
    eqstoich = 0.d0
    allocate(eqh2oid(ngas))
    eqh2oid = 0
    allocate(eqh2ostoich(ngas))
    eqh2ostoich = 0.d0
    allocate(eqlogK(ngas))
    eqlogK = 0.d0
    if (.not.reaction%use_geothermal_hpt) then
      if (option%use_isothermal) then
        allocate(eqlogKcoef(reaction%num_dbase_temperatures, &
                                         ngas))
      else
        allocate(eqlogKcoef(FIVE_INTEGER,ngas))
      endif
    else
      allocate(eqlogKcoef(num_logKs,ngas))
    endif
    
    eqlogKcoef = 0.d0

    ! pack in reaction arrays
    cur_gas_spec => reaction%gas%list
    igas_spec = 1
    do
      if (.not.associated(cur_gas_spec)) exit
      if (cur_gas_spec%itype == gas_itype .or. &
          cur_gas_spec%itype == ACTIVE_AND_PASSIVE_GAS) then      

        gas_names(igas_spec) = cur_gas_spec%name
        ! if the gas is active and passive, we only print the active side
        if (cur_gas_spec%print_me .or. reaction%gas%print_all) then
          if (.not.(cur_gas_spec%itype == ACTIVE_AND_PASSIVE_GAS .and. &
                    gas_itype == PASSIVE_GAS)) then
            gas_print(igas_spec) = PETSC_TRUE
          endif
        endif

        ispec = 0
        do i = 1, cur_gas_spec%dbaserxn%nspec
          if (cur_gas_spec%dbaserxn%spec_ids(i) /= h2o_id) then
            ispec = ispec + 1
            spec_id = cur_gas_spec%dbaserxn%spec_ids(i)
            if (spec_id > h2o_id) spec_id = spec_id - 1
            eqspecid(ispec,igas_spec) = spec_id
            eqstoich(ispec,igas_spec) = cur_gas_spec%dbaserxn%stoich(i)
            
          else ! fill in h2o id and stoich
            eqh2oid(igas_spec) = h2o_id
            eqh2ostoich(igas_spec) = cur_gas_spec%dbaserxn%stoich(i)
          endif
        enddo
        eqspecid(0,igas_spec) = ispec
      
        if (.not.reaction%use_geothermal_hpt) then
          if (option%use_isothermal) then
            eqlogKcoef(:,igas_spec) = cur_gas_spec%dbaserxn%logK
            call Interpolate(temp_high,temp_low,option%reference_temperature, &
                             cur_gas_spec%dbaserxn%logK(itemp_high), &
                             cur_gas_spec%dbaserxn%logK(itemp_low), &
                             eqlogK(igas_spec))
          else
            call ReactionFitLogKCoef(eqlogKcoef(:,igas_spec), &
                                     cur_gas_spec%dbaserxn%logK, &
                                     gas_names(igas_spec), &
                                     option,reaction)
            call ReactionInitializeLogK(eqlogKcoef(:,igas_spec), &
                                        cur_gas_spec%dbaserxn%logK, &
                                        eqlogK(igas_spec), &
                                        option,reaction)
          endif
        else
          eqlogKcoef(:,igas_spec) = cur_gas_spec%dbaserxn%logK
          call ReactionInitializeLogK_hpt(eqlogKcoef(:,igas_spec), &
                                          eqlogK(igas_spec), &
                                          option,reaction)  
        endif
        igas_spec = igas_spec + 1
      endif
      cur_gas_spec => cur_gas_spec%next
    enddo
  endif

end subroutine ReactionDatabaseSetupGases

! ************************************************************************** !

subroutine BasisPrint(reaction,title,option)
  ! 
  ! Prints the basis
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/01/08
  ! 

  use Option_module
  use Reaction_module
  use Reaction_Mineral_Aux_module
  use Reaction_Gas_Aux_module

  implicit none
  
  type(reaction_type) :: reaction
  character(len=*) :: title
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec
  type(gas_species_type), pointer :: cur_gas_spec
  type(mineral_rxn_type), pointer :: cur_mineral
  type(ion_exchange_rxn_type), pointer :: cur_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cur_cation
  
  character(len=MAXSTRINGLENGTH) :: reactant_string, product_string
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: fid

  PetscInt :: ispec, itemp

100 format(a)
110 format(a,f9.4,a)
120 format(a,f9.4,2x,a)
130 format(a,100f11.4)
140 format(a,f6.2)
150 format(a,es11.4,a)

  if (OptionPrintToFile(option)) then
    write(option%fid_out,*)
    write(option%fid_out,*) '! ********************************************&
                    &****************************** !'
    write(option%fid_out,*)
    write(option%fid_out,*) trim(title)
    write(option%fid_out,*)

    write(option%fid_out,*) 'Primary Species:'
    cur_aq_spec => reaction%primary_species_list
    do
      if (.not.associated(cur_aq_spec)) exit
      write(option%fid_out,100,advance='no') '  ' // trim(cur_aq_spec%name)
      if (cur_aq_spec%is_redox) then
        write(option%fid_out,100) '   (redox species)'
      else
        write(option%fid_out,100) ''
      endif
      write(option%fid_out,140) '    Charge: ', cur_aq_spec%Z
      write(option%fid_out,110) '    Molar Mass: ', &
        cur_aq_spec%molar_weight, ' [g/mol]'
      write(option%fid_out,110) '    Debye-Huckel a0: ', &
        cur_aq_spec%a0, ' [Angstrom]'
      if (associated(cur_aq_spec%dbaserxn)) then
        write(option%fid_out,100) '    Equilibrium Aqueous Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_aq_spec%name
        do ispec = 1, cur_aq_spec%dbaserxn%nspec
          write(option%fid_out,120) '      ', &
                          cur_aq_spec%dbaserxn%stoich(ispec), &
                          cur_aq_spec%dbaserxn%spec_name(ispec)
        enddo
        if (reaction%use_geothermal_hpt)then
          write(option%fid_out,130) '      logKCoeff(PT):', &
            (cur_aq_spec%dbaserxn%logK(itemp),&
             itemp=1, reaction%num_dbase_parameters)
        else
          write(option%fid_out,130) '      logK:', &
            (cur_aq_spec%dbaserxn%logK(itemp),itemp=1, &
             reaction%num_dbase_temperatures)
        endif
      endif
      write(option%fid_out,*)
      cur_aq_spec => cur_aq_spec%next
    enddo
    
    cur_aq_spec => reaction%secondary_species_list
    if (associated(cur_aq_spec)) then
      write(option%fid_out,*)
      write(option%fid_out,*) 'Secondary Species:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Secondary Species: None'
    endif
!#define WRITE_LATEX
#ifdef WRITE_LATEX
    fid = 86
    open(fid,file="rxns.txt",action="write")
#endif
    do
      if (.not.associated(cur_aq_spec)) exit
      write(option%fid_out,100) '  ' // trim(cur_aq_spec%name)
      write(option%fid_out,140) '    Charge: ', cur_aq_spec%Z
      write(option%fid_out,110) '    Molar Mass: ', &
        cur_aq_spec%molar_weight,' [g/mol]'
      write(option%fid_out,110) '    Debye-Huckel a0: ', &
        cur_aq_spec%a0, ' [Angstrom]'
      if (associated(cur_aq_spec%dbaserxn)) then
        write(option%fid_out,100) '    Equilibrium Aqueous Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_aq_spec%name
#ifdef WRITE_LATEX
        reactant_string = cur_aq_spec%name
        product_string = ''
#endif
        do ispec = 1, cur_aq_spec%dbaserxn%nspec
          write(option%fid_out,120) '      ', &
            cur_aq_spec%dbaserxn%stoich(ispec), &
            cur_aq_spec%dbaserxn%spec_name(ispec)
#ifdef WRITE_LATEX
          if (dabs(cur_aq_spec%dbaserxn%stoich(ispec)) > 1.d0) then
            write(word,160) &
              int(dabs(cur_aq_spec%dbaserxn%stoich(ispec))+1.e-10), &
              ' ' // trim(cur_aq_spec%dbaserxn%spec_name(ispec))
            word = adjustl(word)
          else
            word = cur_aq_spec%dbaserxn%spec_name(ispec)
          endif
          if (cur_aq_spec%dbaserxn%stoich(ispec) < 0.d0) then
            reactant_string = trim(reactant_string) // ' + ' // trim(word)
          else
            if (len_trim(product_string) > 0) then
              product_string = trim(product_string) // ' + ' // trim(word)
            else
              product_string = word
            endif
          endif
#endif          
        enddo
        if (reaction%use_geothermal_hpt)then
          write(option%fid_out,130) '      logKCoeff(PT):', &
            (cur_aq_spec%dbaserxn%logK(itemp),&
             itemp=1, reaction%num_dbase_parameters)
        else
          write(option%fid_out,130) '      logK:', &
            (cur_aq_spec%dbaserxn%logK(itemp),itemp=1, &
             reaction%num_dbase_temperatures)
        endif

      endif
#ifdef WRITE_LATEX
      write(word,130) '', cur_aq_spec%dbaserxn%logK(2)
      write(fid,*) trim(reactant_string) // ' $~\rightleftharpoons~$ ' // &
        trim(product_string) // ' & ' // trim(adjustl(word)) // ' \\'
#endif
      write(option%fid_out,*)
      cur_aq_spec => cur_aq_spec%next
    enddo
#ifdef WRITE_LATEX
    close(fid)
#endif
    
    cur_gas_spec => reaction%gas%list
    if (associated(cur_gas_spec)) then
      write(option%fid_out,*)
      write(option%fid_out,*) 'Gas Components:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Gas Components: None'
    endif
    do
      if (.not.associated(cur_gas_spec)) exit
      write(option%fid_out,100) '  ' // trim(cur_gas_spec%name)
      write(option%fid_out,110) '    Molar Mass: ', &
        cur_gas_spec%molar_weight,' [g/mol]'
      if (associated(cur_gas_spec%dbaserxn)) then
        write(option%fid_out,100) '    Gas Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_gas_spec%name
        do ispec = 1, cur_gas_spec%dbaserxn%nspec
          write(option%fid_out,120) '      ', &
            cur_gas_spec%dbaserxn%stoich(ispec), &
            cur_gas_spec%dbaserxn%spec_name(ispec)
        enddo
        if (reaction%use_geothermal_hpt)then
           write(option%fid_out,130) '      logKCoeff(PT):', &
             (cur_gas_spec%dbaserxn%logK(itemp),&
              itemp=1, reaction%num_dbase_parameters)
        else
        write(option%fid_out,130) '      logK:', &
          (cur_gas_spec%dbaserxn%logK(itemp),itemp=1, &
           reaction%num_dbase_temperatures)
        endif                                       
      endif
      write(option%fid_out,*)
      cur_gas_spec => cur_gas_spec%next
    enddo
    
    cur_mineral => reaction%mineral%mineral_list
    if (associated(cur_mineral)) then    
      write(option%fid_out,*)
      write(option%fid_out,*) 'Minerals:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Minerals: None'    
    endif
    do
      if (.not.associated(cur_mineral)) exit
      write(option%fid_out,100) '  ' // trim(cur_mineral%name)
      write(option%fid_out,110) '    Molar Mass: ', cur_mineral%molar_weight,' [g/mol]'
      write(option%fid_out,150) '    Molar Volume: ', cur_mineral%molar_volume,' [m^3/mol]'
      if (associated(cur_mineral%tstrxn)) then
        write(option%fid_out,100) '    Mineral Reaction: '
        write(option%fid_out,120) '      ', -1.d0, cur_mineral%name
        do ispec = 1, cur_mineral%dbaserxn%nspec
          write(option%fid_out,120) '      ', cur_mineral%dbaserxn%stoich(ispec), &
                          cur_mineral%dbaserxn%spec_name(ispec)
        enddo
        if (reaction%use_geothermal_hpt)then
          write(option%fid_out,130) '      logKCoeff(PT):', (cur_mineral%dbaserxn%logK(itemp),&
                                    itemp=1, reaction%num_dbase_parameters)
        else        
          write(option%fid_out,130) '      logK:', (cur_mineral%dbaserxn%logK(itemp),itemp=1, &
                                       reaction%num_dbase_temperatures)
        endif
      endif
      write(option%fid_out,*)
      cur_mineral => cur_mineral%next
    enddo
    

    cur_ionx_rxn => reaction%ion_exchange_rxn_list
    if (associated(cur_ionx_rxn)) then
      write(option%fid_out,*)
      write(option%fid_out,*) 'Ion Exchange Reactions:'
    else
      write(option%fid_out,*)
      write(option%fid_out,*) 'Ion Exchange Reactions: None'
    endif
    do
      if (.not.associated(cur_ionx_rxn)) exit
      write(option%fid_out,*) '  Mineral: ', trim(cur_ionx_rxn%mineral_name)
      write(option%fid_out,150) '      CEC: ', cur_ionx_rxn%CEC
      cur_cation => cur_ionx_rxn%cation_list
      if (associated(cur_cation)) then
        write(option%fid_out,*) '  Cations:'
      else
        write(option%fid_out,*) '  Cations: None'
      endif
      do
        if (.not.associated(cur_cation)) exit
        write(option%fid_out,150) '      ' // trim(cur_cation%name), &
          cur_cation%k
        cur_cation => cur_cation%next
      enddo
      write(option%fid_out,*)
      cur_ionx_rxn => cur_ionx_rxn%next
    enddo
    
    write(option%fid_out,*)
    write(option%fid_out,*) '! ********************************************&
                    &****************************** !'
    write(option%fid_out,*)
  endif

end subroutine BasisPrint

end module Reaction_Database_module
