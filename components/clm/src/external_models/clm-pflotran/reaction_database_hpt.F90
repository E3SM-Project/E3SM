module Reaction_Database_hpt_module

  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_module
  use Reaction_Database_Aux_module
  use Reaction_Surface_Complexation_Aux_module
  use Reaction_Mineral_Aux_module
  use Reaction_Gas_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  public :: DatabaseRead_hpt, BasisInit_hpt
            
contains

! ************************************************************************** !

subroutine DatabaseRead_hpt(reaction,option)
  ! 
  ! Collects parameters from geochemical database
  ! 
  ! Author: ???
  ! Date: ???
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  type(aq_species_type), pointer :: cur_aq_spec, cur_aq_spec2
  type(gas_species_type), pointer :: cur_gas_spec, cur_gas_spec2
  type(mineral_rxn_type), pointer :: cur_mineral, cur_mineral2
  type(colloid_type), pointer :: cur_colloid
  type(surface_complexation_rxn_type), pointer :: cur_srfcplx_rxn
  type(surface_complex_type), pointer :: cur_srfcplx, cur_srfcplx2
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: null_name
  
  PetscBool :: flag, found
  PetscInt :: ispec, itemp, i
  PetscReal :: stoich
  PetscReal :: temp_real
  type(input_type), pointer :: input
  PetscInt :: iostat
  PetscInt :: num_nulls
!TODO(geh)
#if 0  
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
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    cur_mineral%id = -abs(cur_mineral%id)
    cur_mineral => cur_mineral%next
  enddo
  cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_srfcplx_rxn)) exit
    cur_srfcplx => cur_srfcplx_rxn%complex_list
    do  
      if (.not.associated(cur_srfcplx)) exit
      cur_srfcplx%id = -abs(cur_srfcplx%id)
      cur_srfcplx => cur_srfcplx%next
    enddo
    cur_srfcplx_rxn => cur_srfcplx_rxn%next
  enddo  
  
  if (len_trim(reaction%database_filename) < 2) then
    option%io_buffer = 'Database filename not included in input deck.'
    call printErrMsg(option)
  endif
  input => InputCreate(IUNIT_TEMP,reaction%database_filename,option)

  ! read temperatures
  call InputReadPflotranString(input,option)
  ! remove comment
  call InputReadQuotedWord(input,option,name,PETSC_TRUE)
  call InputReadInt(input,option,reaction%num_dbase_parameters)
  call InputErrorMsg(input,option,'Number of database parameters','DATABASE')  
 ! allocate(reaction%dbase_temperatures(reaction%num_dbase_temperatures))
 ! reaction%dbase_temperatures = 0.d0  
!  do itemp = 1, reaction%num_dbase_temperatures
!    call InputReadDouble(input,option,reaction%dbase_temperatures(itemp))
!    call InputErrorMsg(input,option,'Database temperatures','DATABASE')            
!  enddo

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
      if (num_nulls >= 5) exit
      cycle
    endif
    
    select case(num_nulls)
      case(0,1) ! primary and secondary aq species and colloids
        cur_aq_spec => reaction%primary_species_list
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_aq_spec)) exit
          if (StringCompare(name,cur_aq_spec%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE
            ! change negative id to positive, indicating it was found in database
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
          ! change negative id to positive, indicating it was found in database
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
          ! change negative id to positive, indicating it was found in database
            cur_colloid%id = abs(cur_colloid%id)

            ! skip the Debye-Huckel ion size parameter (a0)
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'Colloid skip a0','DATABASE')            
            ! skipo the valence
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'Colloid skip Z','DATABASE')            
            ! read the molar weight
            call InputReadDouble(input,option,cur_colloid%molar_weight)
            call InputErrorMsg(input,option,'Colloid molar weight','DATABASE')
            
            cycle ! avoid the aqueous species parameters below
          endif
          cur_colloid => cur_colloid%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        if (num_nulls > 0) then ! secondary species in database
          ! create aqueous equilibrium reaction
          if (.not.associated(cur_aq_spec%dbaserxn)) &
            cur_aq_spec%dbaserxn => DatabaseRxnCreate()
          ! read the number of primary species in secondary rxn
          call InputReadInt(input,option,cur_aq_spec%dbaserxn%nspec)
          call InputErrorMsg(input,option,'Number of species in aqueous complex', &
                          'DATABASE')  
          ! allocate arrays for rxn
          allocate(cur_aq_spec%dbaserxn%spec_name(cur_aq_spec%dbaserxn%nspec))
          cur_aq_spec%dbaserxn%spec_name = ''
          allocate(cur_aq_spec%dbaserxn%stoich(cur_aq_spec%dbaserxn%nspec))
          cur_aq_spec%dbaserxn%stoich = 0.d0
          allocate(cur_aq_spec%dbaserxn%logKCoeff_hpt(reaction%num_dbase_parameters))
          cur_aq_spec%dbaserxn%logKCoeff_hpt = 0.d0
          ! read in species and stoichiometries
          do ispec = 1, cur_aq_spec%dbaserxn%nspec
            call InputReadDouble(input,option,cur_aq_spec%dbaserxn%stoich(ispec))
            call InputErrorMsg(input,option,'EQRXN species stoichiometry','DATABASE')            
            call InputReadQuotedWord(input,option,cur_aq_spec%dbaserxn%spec_name(ispec),PETSC_TRUE)
            call InputErrorMsg(input,option,'EQRXN species name','DATABASE')            
          enddo
          do itemp = 1, reaction%num_dbase_parameters
            call InputReadDouble(input,option,cur_aq_spec%dbaserxn%logKCoeff_hpt(itemp))
            call InputErrorMsg(input,option,'EQRXN logKs Coeff','DATABASE')            
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
          ! change negative id to positive, indicating it was found in database
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
        allocate(cur_gas_spec%dbaserxn%logKCoeff_hpt(reaction%num_dbase_parameters))
        cur_gas_spec%dbaserxn%logKCoeff_hpt = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_gas_spec%dbaserxn%nspec
          call InputReadDouble(input,option,cur_gas_spec%dbaserxn%stoich(ispec))
          call InputErrorMsg(input,option,'GAS species stoichiometry','DATABASE')            
          call InputReadQuotedWord(input,option,cur_gas_spec%dbaserxn%spec_name(ispec),PETSC_TRUE)
          call InputErrorMsg(input,option,'GAS species name','DATABASE')            
        enddo
        do itemp = 1, reaction%num_dbase_parameters
          call InputReadDouble(input,option,cur_gas_spec%dbaserxn%logKCoeff_hpt(itemp))
          call InputErrorMsg(input,option,'GAS logKs Coeff','DATABASE')            
        enddo
        ! read the molar weight
        call InputReadDouble(input,option,cur_gas_spec%molar_weight)
        call InputErrorMsg(input,option,'GAS molar weight','DATABASE')     
        
               
      case(3) ! minerals
        cur_mineral => reaction%mineral_list
        if (.not.associated(cur_mineral)) cycle
        found = PETSC_FALSE
        do
          if (found .or. .not.associated(cur_mineral)) exit
          if (StringCompare(name,cur_mineral%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE          
          ! change negative id to positive, indicating it was found in database
            cur_mineral%id = abs(cur_mineral%id)
            exit
          endif
          cur_mineral => cur_mineral%next
        enddo
        
        if (.not.found) cycle ! go to next line in database
        
        ! read the molar volume
        call InputReadDouble(input,option,cur_mineral%molar_volume)
        call InputErrorMsg(input,option,'MINERAL molar volume','DATABASE')            
        ! convert from cm^3/mol to m^3/mol
        cur_mineral%molar_volume = cur_mineral%molar_volume*1.d-6
        ! create mineral reaction
        if (.not.associated(cur_mineral%tstrxn)) then
          cur_mineral%tstrxn => TransitionStateTheoryRxnCreate()
        endif
        ! read the number of aqueous species in mineral rxn
        cur_mineral%dbaserxn => DatabaseRxnCreate()
        call InputReadInt(input,option,cur_mineral%dbaserxn%nspec)
        call InputErrorMsg(input,option,'Number of species in mineral reaction', &
                        'DATABASE')  
        ! allocate arrays for rxn
        allocate(cur_mineral%dbaserxn%spec_name(cur_mineral%dbaserxn%nspec))
        cur_mineral%dbaserxn%spec_name = ''
        allocate(cur_mineral%dbaserxn%stoich(cur_mineral%dbaserxn%nspec))
        cur_mineral%dbaserxn%stoich = 0.d0
        allocate(cur_mineral%dbaserxn%logKCoeff_hpt(reaction%num_dbase_parameters))
        cur_mineral%dbaserxn%logKCoeff_hpt = 0.d0
        ! read in species and stoichiometries
        do ispec = 1, cur_mineral%dbaserxn%nspec
          call InputReadDouble(input,option,cur_mineral%dbaserxn%stoich(ispec))
          call InputErrorMsg(input,option,'MINERAL species stoichiometry','DATABASE')            
          call InputReadQuotedWord(input,option,cur_mineral%dbaserxn% &
                                   spec_name(ispec),PETSC_TRUE)
          call InputErrorMsg(input,option,'MINERAL species name','DATABASE')            
        enddo
        do itemp = 1, reaction%num_dbase_parameters
          call InputReadDouble(input,option,cur_mineral%dbaserxn%logKCoeff_hpt(itemp))
          call InputErrorMsg(input,option,'MINERAL logKs','DATABASE')            
        enddo
        ! read the molar weight
        call InputReadDouble(input,option,cur_mineral%molar_weight)
        call InputErrorMsg(input,option,'MINERAL molar weight','DATABASE')            
        
        
      case(4) ! surface complexes
        cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
        found = PETSC_FALSE
        do
          if (.not.associated(cur_srfcplx_rxn)) exit
          cur_srfcplx => cur_srfcplx_rxn%complex_list
          do
            if (.not.associated(cur_srfcplx)) exit
            if (StringCompare(name,cur_srfcplx%name,MAXWORDLENGTH)) then
              found = PETSC_TRUE          
            ! change negative id to positive, indicating it was found in database
              cur_srfcplx%id = abs(cur_srfcplx%id)
              exit
            endif
            cur_srfcplx => cur_srfcplx%next
          enddo
          if (found) exit
          cur_srfcplx_rxn => cur_srfcplx_rxn%next
        enddo
        
        if (.not.found) cycle ! go to next line in database

        if (.not.associated(cur_srfcplx%dbaserxn)) &
          cur_srfcplx%dbaserxn => DatabaseRxnCreate()
            
        ! read the number of aqueous species in surface complexation rxn
        call InputReadInt(input,option,cur_srfcplx%dbaserxn%nspec)
        call InputErrorMsg(input,option,'Number of species in surface complexation reaction', &
                        'DATABASE')  
        ! decrement number of species since free site will not be included
        cur_srfcplx%dbaserxn%nspec = cur_srfcplx%dbaserxn%nspec - 1
        ! allocate arrays for rxn
        allocate(cur_srfcplx%dbaserxn%spec_name(cur_srfcplx%dbaserxn%nspec))
        cur_srfcplx%dbaserxn%spec_name = ''
        allocate(cur_srfcplx%dbaserxn%stoich(cur_srfcplx%dbaserxn%nspec))
        cur_srfcplx%dbaserxn%stoich = 0.d0
        allocate(cur_srfcplx%dbaserxn%logKCoeff_hpt(reaction%num_dbase_parameters))
        cur_srfcplx%dbaserxn%logKCoeff_hpt = 0.d0
        ! read in species and stoichiometries
        ispec = 0
        found = PETSC_FALSE
        do i = 1, cur_srfcplx%dbaserxn%nspec+1 ! recall that nspec was decremented above
          call InputReadDouble(input,option,stoich)
          call InputErrorMsg(input,option,'SURFACE COMPLEX species stoichiometry','DATABASE')            
          call InputReadQuotedWord(input,option,name,PETSC_TRUE)
          call InputErrorMsg(input,option,'SURFACE COMPLEX species name','DATABASE')            
          if (StringCompare(name,cur_srfcplx_rxn%free_site_name,MAXWORDLENGTH)) then
            found = PETSC_TRUE
            cur_srfcplx%free_site_stoich = stoich
          else
            ispec = ispec + 1
            cur_srfcplx%dbaserxn%stoich(ispec) = stoich
            cur_srfcplx%dbaserxn%spec_name(ispec) = name
          endif
        enddo
        if (.not.found) then
          option%io_buffer = 'Free site name: ' // &
                             trim(cur_srfcplx_rxn%free_site_name) // &
                             ' not found in surface complex:' // &
                             trim(cur_srfcplx%name)
          call printErrMsg(option)
        endif
        do itemp = 1, reaction%num_dbase_parameters
          call InputReadDouble(input,option,cur_srfcplx%dbaserxn%logKCoeff_hpt(itemp))
          call InputErrorMsg(input,option,'SURFACE COMPLEX logKs','DATABASE')            
        enddo
        ! read the valence
        call InputReadDouble(input,option,cur_srfcplx%Z)
        call InputErrorMsg(input,option,'Surface Complex Z','DATABASE')            

      
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
  cur_mineral => reaction%mineral_list
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
  
  ! surface complexes
  cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_srfcplx_rxn)) exit
    cur_srfcplx => cur_srfcplx_rxn%complex_list
    do
      if (.not.associated(cur_srfcplx)) exit
      cur_srfcplx2 => cur_srfcplx%next
      do
        if (.not.associated(cur_srfcplx2)) exit
        if (cur_srfcplx%id /= cur_srfcplx2%id .and. &
            StringCompare(cur_srfcplx%name, &
                            cur_srfcplx2%name,MAXWORDLENGTH)) then
          flag = PETSC_TRUE
          option%io_buffer = 'Surface complex (' // &
                             trim(cur_srfcplx2%name) // &
                      ') duplicated in input file surface complex reaction.'
          call printMsg(option)                          
        endif
        cur_srfcplx2 => cur_srfcplx2%next
      enddo
      cur_srfcplx => cur_srfcplx%next
    enddo
    cur_srfcplx_rxn => cur_srfcplx_rxn%next
  enddo  
  
  if (flag) call printErrMsg(option,'Species duplicated in input file.')

  ! check that all species, etc. were read
  flag = PETSC_FALSE
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
    cur_gas_spec => cur_gas_spec%next
  enddo  
  cur_mineral => reaction%mineral_list
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%id < 0) then
      flag = PETSC_TRUE
      option%io_buffer = 'Mineral (' // trim(cur_mineral%name) // &
               ') not found in database.'
      call printMsg(option)
    endif
    cur_mineral => cur_mineral%next
  enddo
  cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_srfcplx_rxn)) exit
    cur_srfcplx => cur_srfcplx_rxn%complex_list
    do
      if (.not.associated(cur_srfcplx)) exit
      if (cur_srfcplx%id < 0) then
        flag = PETSC_TRUE
        option%io_buffer = 'Surface species (' // trim(cur_srfcplx%name) // &
                 ') not found in database.'
        call printMsg(option)
      endif
      cur_srfcplx => cur_srfcplx%next
    enddo  
    cur_srfcplx_rxn => cur_srfcplx_rxn%next
  enddo 
    
  if (flag) call printErrMsg(option,'Species not found in database.')

  call InputDestroy(input)
!TODO(geh)
#endif
end subroutine DatabaseRead_hpt

! ************************************************************************** !

subroutine BasisInit_hpt(reaction,option)
  ! 
  ! Initializes the basis for geochemistry
  ! 
  ! Author: ???
  ! Date: ???
  ! 

  use Option_module
  use String_module
  use Utility_module
  use Input_Aux_module

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
  type(surface_complexation_rxn_type), pointer :: cur_srfcplx_rxn
  type(surface_complex_type), pointer :: cur_srfcplx
  type(surface_complex_type), pointer :: cur_srfcplx2
  type(ion_exchange_rxn_type), pointer :: cur_ionx_rxn
  type(ion_exchange_cation_type), pointer :: cur_cation
  type(general_rxn_type), pointer :: cur_general_rxn
  type(kd_rxn_type), pointer :: cur_kd_rxn
  type(colloid_type), pointer :: cur_colloid
  type(database_rxn_type), pointer :: dbaserxn
  type(transition_state_rxn_type), pointer :: tstrxn
  type(transition_state_prefactor_type), pointer :: cur_prefactor
  type(ts_prefactor_species_type), pointer :: cur_prefactor_species

  character(len=MAXWORDLENGTH), allocatable :: old_basis_names(:)
  character(len=MAXWORDLENGTH), allocatable :: new_basis_names(:)

  character(len=MAXWORDLENGTH), parameter :: h2oname = 'H2O'
  character(len=MAXWORDLENGTH) :: word, word2
  character(len=MAXSTRINGLENGTH) :: string, string2
  
  PetscInt, parameter :: h2o_id = 1

  PetscReal :: logK(reaction%num_dbase_temperatures)
  PetscReal, allocatable :: transformation(:,:), old_basis(:,:), new_basis(:,:)
  PetscReal, allocatable :: stoich_new(:), stoich_prev(:), logKCoeffvector(:,:)
  PetscInt, allocatable :: indices(:)
  
  PetscReal, allocatable :: pri_matrix(:,:), sec_matrix(:,:)
  PetscReal, allocatable :: sec_matrix_inverse(:,:)
  PetscReal, allocatable :: stoich_matrix(:,:)
  PetscReal, allocatable :: unit_vector(:)
  character(len=MAXWORDLENGTH), allocatable :: pri_names(:)
  character(len=MAXWORDLENGTH), allocatable :: sec_names(:)
  character(len=MAXWORDLENGTH), allocatable :: gas_names(:)
  PetscReal, allocatable :: logKCoeffvector_swapped(:,:)
  PetscBool, allocatable :: flags(:)
  PetscBool :: negative_flag
  PetscReal :: value
  
  PetscInt :: ispec, itemp
  PetscInt :: spec_id
  PetscInt :: ncomp_h2o, ncomp_secondary
  PetscInt :: icount_old, icount_new, icount, icount2
  PetscInt :: i, j, irow, icol
  PetscInt :: icomp, icplx, irxn
  PetscInt :: ipri_spec, isec_spec, imnrl, igas_spec, ikinmnrl, icoll
  PetscInt :: i_old, i_new
  PetscInt :: isrfcplx
  PetscInt :: ication
  PetscInt :: idum
  PetscReal :: temp_high, temp_low
  PetscInt :: itemp_high, itemp_low
  PetscInt :: species_count, max_species_count
  PetscInt :: forward_count, max_forward_count
  PetscInt :: backward_count, max_backward_count
  PetscInt :: midpoint
  PetscInt :: max_num_prefactors, max_num_prefactor_species
  
  PetscBool :: compute_new_basis
  PetscBool :: found
  PetscErrorCode :: ierr
!TODO(geh)
#if 0
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

  ! # of components sorbed to colloids
  reaction%naqcomp = GetPrimarySpeciesCount(reaction)
  reaction%ncoll = GetColloidCount(reaction)
  reaction%neqcplx = GetSecondarySpeciesCount(reaction)
  reaction%ngas = GetGasCount(reaction)

  reaction%ncollcomp = reaction%naqcomp ! set to naqcomp for now, will be adjusted later
  reaction%offset_aqueous = 0
  reaction%offset_colloid = reaction%offset_aqueous + reaction%naqcomp
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

  ncomp_secondary = reaction%neqcplx+reaction%ngas
  
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
      option%io_buffer = 'Too few reactions read from database for ' // & 
        'number of secondary species defined.'
    else
      option%io_buffer = 'Too many reactions read from database for ' // & 
        'number of secondary species defined.  Perhaps REDOX ' // &
        'SPECIES need to be defined?'
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
  allocate(gas_names(reaction%ngas))
  gas_names = ''
  
  allocate(logKCoeffvector(reaction%num_dbase_parameters,ncomp_secondary))
  logKCoeffvector = 0.d0
  
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
      logKCoeffvector(:,icount) = cur_pri_aq_spec%dbaserxn%logK
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
      logKCoeffvector(:,icount) = cur_sec_aq_spec%dbaserxn%logKCoeff_hpt
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
      logKCoeffvector(:,icount) = cur_gas_spec%dbaserxn%logKCoeff_hpt
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
 
  call ludcmp(sec_matrix,ncomp_secondary,indices,idum)
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

  allocate(logKCoeffvector_swapped(reaction%num_dbase_parameters,ncomp_secondary))
  logKCoeffvector_swapped = 0.d0
  
  do j = 1, ncomp_secondary
    do i = 1, reaction%num_dbase_parameters
      logKCoeffvector_swapped(i,j) = logKCoeffvector_swapped(i,j) - &
        dot_product(sec_matrix_inverse(j,1:ncomp_secondary), &
                    logKCoeffvector(i,1:ncomp_secondary))
    enddo
  enddo
    
  deallocate(pri_matrix)
  deallocate(sec_matrix)
  deallocate(indices)
  deallocate(unit_vector)
  deallocate(sec_matrix_inverse)
  deallocate(logKCoeffvector_swapped)
  
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
    allocate(cur_sec_aq_spec%dbaserxn%logKCoeff_hpt(reaction%num_dbase_parameters))
    cur_sec_aq_spec%dbaserxn%logKCoeff_hpt = 0.d0

    ispec = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        ispec = ispec + 1
        cur_sec_aq_spec%dbaserxn%spec_name(ispec) = pri_names(icol)
        cur_sec_aq_spec%dbaserxn%stoich(ispec) = stoich_matrix(icount,icol)
        cur_sec_aq_spec%dbaserxn%spec_ids(ispec) = icol
      endif
    enddo

    cur_sec_aq_spec%dbaserxn%logKCoeff_hpt = logKCoeffvector_swapped(:,icount)

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
    allocate(cur_gas_spec%dbaserxn%logKCoeff_hpt(reaction%num_dbase_parameters))
    cur_gas_spec%dbaserxn%logKCoeff_hpt = 0.d0

    ispec = 0
    do icol = 1, ncomp_h2o
      if (dabs(stoich_matrix(icount,icol)) > 1.d-40) then
        ispec = ispec + 1
        cur_gas_spec%dbaserxn%spec_name(ispec) = pri_names(icol)
        cur_gas_spec%dbaserxn%stoich(ispec) = stoich_matrix(icount,icol)
        cur_gas_spec%dbaserxn%spec_ids(ispec) = icol
      endif
    enddo

    cur_gas_spec%dbaserxn%logKCoeff_hpt = logKCoeffvector_swapped(:,icount)

    cur_gas_spec => cur_gas_spec%next
  enddo

  new_basis_names = pri_names

  deallocate(stoich_matrix)
  deallocate(logKCoeffvector_swapped)

  deallocate(pri_names)
  deallocate(sec_names)
  deallocate(gas_names)

  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_srfcplx_rxn)
  nullify(cur_srfcplx)
    
  ! first off, lets remove all the secondary gases from all other reactions
  cur_gas_spec => reaction%gas%list
  do
    if (.not.associated(cur_gas_spec)) exit
    
    ! gases in mineral reactions
    cur_mineral => reaction%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%dbaserxn%nspec) exit
          if (StringCompare(cur_gas_spec%name, &
                              cur_mineral%dbaserxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn_hpt(cur_gas_spec%name, &
                                                 cur_gas_spec%dbaserxn, &
                                                 cur_mineral%dbaserxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo
    nullify(cur_mineral)

    ! gases in surface complex reactions
    cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_srfcplx_rxn)) exit
      cur_srfcplx2 => cur_srfcplx_rxn%complex_list
      do
        if (.not.associated(cur_srfcplx2)) exit
        
        if (associated(cur_srfcplx2%dbaserxn)) then
          ispec = 1
          do
            if (ispec > cur_srfcplx2%dbaserxn%nspec) exit
            if (StringCompare(cur_gas_spec%name, &
                                cur_srfcplx2%dbaserxn%spec_name(ispec), &
                                MAXWORDLENGTH)) then
              call BasisSubSpecInGasOrSecRxn_hpt(cur_gas_spec%name, &
                                                cur_gas_spec%dbaserxn, &
                                                cur_srfcplx2%dbaserxn)
              ispec = 0
            endif
            ispec = ispec + 1
          enddo
        endif
        cur_srfcplx2 => cur_srfcplx2%next
      enddo
      nullify(cur_srfcplx2)
      cur_srfcplx_rxn => cur_srfcplx_rxn%next
    enddo
    nullify(cur_srfcplx_rxn)

    cur_gas_spec => cur_gas_spec%next
  enddo

  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_srfcplx_rxn)
  nullify(cur_srfcplx)

  ! secondary aqueous species
  cur_sec_aq_spec => reaction%secondary_species_list
  do

    if (.not.associated(cur_sec_aq_spec)) exit
    
    ! secondary aqueous species in mineral reactions
    cur_mineral => reaction%mineral_list
    do
      if (.not.associated(cur_mineral)) exit
      
      if (associated(cur_mineral%tstrxn)) then
        ispec = 1
        do
          if (ispec > cur_mineral%dbaserxn%nspec) exit
          if (StringCompare(cur_sec_aq_spec%name, &
                              cur_mineral%dbaserxn%spec_name(ispec), &
                              MAXWORDLENGTH)) then
            call BasisSubSpeciesInMineralRxn_hpt(cur_sec_aq_spec%name, &
                                                 cur_sec_aq_spec%dbaserxn, &
                                                 cur_mineral%dbaserxn)
            ispec = 0
          endif
          ispec = ispec + 1
        enddo
      endif
      cur_mineral => cur_mineral%next
    enddo

    ! secondary aqueous species in surface complex reactions
    cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_srfcplx_rxn)) exit
      cur_srfcplx2 => cur_srfcplx_rxn%complex_list
      do
        if (.not.associated(cur_srfcplx2)) exit
        
        if (associated(cur_srfcplx2%dbaserxn)) then
          ispec = 1
          do
            if (ispec > cur_srfcplx2%dbaserxn%nspec) exit
            if (StringCompare(cur_sec_aq_spec%name, &
                                cur_srfcplx2%dbaserxn%spec_name(ispec), &
                                MAXWORDLENGTH)) then
              call BasisSubSpecInGasOrSecRxn_hpt(cur_sec_aq_spec%name, &
                                                 cur_sec_aq_spec%dbaserxn, &
                                                 cur_srfcplx2%dbaserxn)
              ispec = 0
            endif
            ispec = ispec + 1
          enddo
        endif
        cur_srfcplx2 => cur_srfcplx2%next
      enddo
      nullify(cur_srfcplx2)
      cur_srfcplx_rxn => cur_srfcplx_rxn%next
    enddo
    nullify(cur_srfcplx_rxn)    
    
    cur_sec_aq_spec => cur_sec_aq_spec%next
  enddo
  
  nullify(cur_sec_aq_spec)
  nullify(cur_gas_spec)
  nullify(cur_mineral)
  nullify(cur_srfcplx_rxn)
  nullify(cur_srfcplx)

  ! substitute new basis into mineral and surface complexation rxns,
  ! if necessary
  cur_mineral => reaction%mineral_list
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

  cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
  do
    if (.not.associated(cur_srfcplx_rxn)) exit
    cur_srfcplx => cur_srfcplx_rxn%complex_list
    do
      if (.not.associated(cur_srfcplx)) exit
      if (.not.associated(cur_srfcplx%dbaserxn%spec_ids)) then
        allocate(cur_srfcplx%dbaserxn%spec_ids(cur_srfcplx%dbaserxn%nspec))
        cur_srfcplx%dbaserxn%spec_ids = 0
      endif
      call BasisAlignSpeciesInRxn(ncomp_h2o,new_basis_names, &
                                  cur_srfcplx%dbaserxn%nspec, &
                                  cur_srfcplx%dbaserxn%spec_name, &
                                  cur_srfcplx%dbaserxn%stoich, &
                                  cur_srfcplx%dbaserxn%spec_ids, &
                                  cur_srfcplx%name,option) 
      cur_srfcplx => cur_srfcplx%next
    enddo
    nullify(cur_srfcplx)
    cur_srfcplx_rxn => cur_srfcplx_rxn%next
  enddo
  nullify(cur_srfcplx_rxn)  

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
    allocate(reaction%secondary_species_names(reaction%neqcplx))
    reaction%secondary_species_names = ''

    allocate(reaction%secondary_species_print(reaction%neqcplx))
    reaction%secondary_species_print = PETSC_FALSE

    allocate(reaction%eqcplx_basis_names(reaction%naqcomp,reaction%neqcplx))
    reaction%eqcplx_basis_names = ''

    allocate(reaction%eqcplxspecid(0:reaction%naqcomp,reaction%neqcplx))
    reaction%eqcplxspecid = 0

    allocate(reaction%eqcplxstoich(0:reaction%naqcomp,reaction%neqcplx))
    reaction%eqcplxstoich = 0.d0

    allocate(reaction%eqcplxh2oid(reaction%neqcplx))
    reaction%eqcplxh2oid = 0

    allocate(reaction%eqcplxh2ostoich(reaction%neqcplx))
    reaction%eqcplxh2ostoich = 0.d0

    allocate(reaction%eqcplx_logK(reaction%neqcplx))
    reaction%eqcplx_logK = 0.d0

    allocate(reaction%eqcplx_logKcoef(reaction%num_dbase_parameters,reaction%neqcplx))
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
      reaction%secondary_species_print(isec_spec) = cur_sec_aq_spec%print_me .or. &
                                            reaction%print_all_secondary_species
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
          reaction%eqcplxstoich(ispec,isec_spec) = cur_sec_aq_spec%dbaserxn%stoich(i)
            
        else ! fill in h2o id and stoich
          reaction%eqcplxh2oid(isec_spec) = h2o_id
          reaction%eqcplxh2ostoich(isec_spec) = cur_sec_aq_spec%dbaserxn%stoich(i)
        endif
      enddo
      reaction%eqcplxspecid(0,isec_spec) = ispec
!#if 0
!     TODO(Peter): fix argument list
      call ReactionInitializeLogK_hpt(reaction%eqcplx_logKcoef(:,isec_spec), &
                                      reaction%eqcplx_logK(isec_spec), &
                                      option,reaction)
!#endif
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
  reaction%ngas = GetGasCount(reaction)
  
  if (reaction%ngas > 0) then
    allocate(reaction%gas_species_names(reaction%ngas))
    reaction%gas_species_names = ''
    allocate(reaction%gas_species_print(reaction%ngas))
    reaction%gas_species_print = PETSC_FALSE
    allocate(reaction%eqgasspecid(0:reaction%naqcomp,reaction%ngas))
    reaction%eqgasspecid = 0
    allocate(reaction%eqgasstoich(0:reaction%naqcomp,reaction%ngas))
    reaction%eqgasstoich = 0.d0
    allocate(reaction%eqgash2oid(reaction%ngas))
    reaction%eqgash2oid = 0
    allocate(reaction%eqgash2ostoich(reaction%ngas))
    reaction%eqgash2ostoich = 0.d0
    allocate(reaction%eqgas_logK(reaction%ngas))
    reaction%eqgas_logK = 0.d0
    allocate(reaction%eqgas_logKcoef(reaction%num_dbase_parameters, &
                                     reaction%ngas))
    reaction%eqgas_logKcoef = 0.d0

    ! pack in reaction arrays
    cur_gas_spec => reaction%gas%list
    igas_spec = 1
    do
      if (.not.associated(cur_gas_spec)) exit

      reaction%gas_species_names(igas_spec) = cur_gas_spec%name
      reaction%gas_species_print(igas_spec) = cur_gas_spec%print_me .or. &
                                              reaction%print_all_gas_species
      ispec = 0
      do i = 1, cur_gas_spec%dbaserxn%nspec
        if (cur_gas_spec%dbaserxn%spec_ids(i) /= h2o_id) then
          ispec = ispec + 1
          spec_id = cur_gas_spec%dbaserxn%spec_ids(i)
          if (spec_id > h2o_id) spec_id = spec_id - 1
          reaction%eqgasspecid(ispec,igas_spec) = spec_id
          reaction%eqgasstoich(ispec,igas_spec) = &
            cur_gas_spec%dbaserxn%stoich(i)
            
        else ! fill in h2o id and stoich
          reaction%eqgash2oid(igas_spec) = h2o_id
          reaction%eqgash2ostoich(igas_spec) = &
            cur_gas_spec%dbaserxn%stoich(i)
        endif
      enddo
      reaction%eqgasspecid(0,igas_spec) = ispec
!#if 0      
!     TODO(Peter): fix argument list
      call ReactionInitializeLogK_hpt(reaction%eqgas_logKcoef(:,igas_spec), &
                                   reaction%eqgas_logK(igas_spec), &
                                  option,reaction)
!#endif  
      igas_spec = igas_spec + 1
      cur_gas_spec => cur_gas_spec%next
    enddo

  endif
  nullify(cur_gas_spec)
  igas_spec = -1 ! to catch bugs

  ! minerals
  ! Count the number of kinetic mineral reactions, max number of prefactors in a
  !   tst reaction, and the maximum number or species in a prefactor
  reaction%nkinmnrl = 0
  max_num_prefactors = 0
  max_num_prefactor_species = 0
  cur_mineral => reaction%mineral_list
  !
  do
    if (.not.associated(cur_mineral)) exit
    if (cur_mineral%itype == MINERAL_KINETIC .and. &
        associated(cur_mineral%tstrxn)) then
      ! increment number of kinetic minerals
      reaction%nkinmnrl = reaction%nkinmnrl + 1
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

  if (reaction%nmnrl > 0) then
    allocate(reaction%mineral_names(reaction%nmnrl))
    reaction%mineral_names = ''
    allocate(reaction%mnrlspecid(0:reaction%naqcomp,reaction%nmnrl))
    reaction%mnrlspecid = 0
    allocate(reaction%mnrlstoich(reaction%naqcomp,reaction%nmnrl))
    reaction%mnrlstoich = 0.d0
    allocate(reaction%mnrlh2oid(reaction%nmnrl))
    reaction%mnrlh2oid = 0
    allocate(reaction%mnrlh2ostoich(reaction%nmnrl))
    reaction%mnrlh2ostoich = 0.d0
    allocate(reaction%mnrl_logK(reaction%nmnrl))
    reaction%mnrl_logK = 0.d0
    allocate(reaction%mnrl_print(reaction%nmnrl))
    reaction%mnrl_print = PETSC_FALSE

    allocate(reaction%mnrl_logKcoef(reaction%num_dbase_parameters, &
                                    reaction%nmnrl))
    reaction%mnrl_logKcoef = 0.d0

    if (reaction%nkinmnrl > 0) then
      allocate(reaction%kinmnrl_names(reaction%nkinmnrl))
      reaction%kinmnrl_names = ''
      allocate(reaction%kinmnrl_print(reaction%nkinmnrl))
      reaction%kinmnrl_print = PETSC_FALSE
      allocate(reaction%kinmnrlspecid(0:reaction%naqcomp,reaction%nkinmnrl))
      reaction%kinmnrlspecid = 0
      allocate(reaction%kinmnrlstoich(reaction%naqcomp,reaction%nkinmnrl))
      reaction%kinmnrlstoich = 0.d0
      allocate(reaction%kinmnrlh2oid(reaction%nkinmnrl))
      reaction%kinmnrlh2oid = 0
      allocate(reaction%kinmnrlh2ostoich(reaction%nkinmnrl))
      reaction%kinmnrlh2ostoich = 0.d0
      allocate(reaction%kinmnrl_logK(reaction%nkinmnrl))
      reaction%kinmnrl_logK = 0.d0

      allocate(reaction%kinmnrl_logKcoef(reaction%num_dbase_parameters, &
                                         reaction%nkinmnrl))
      reaction%kinmnrl_logKcoef = 0.d0

      ! TST Rxn variables
      allocate(reaction%kinmnrl_affinity_threshold(reaction%nkinmnrl))
      reaction%kinmnrl_affinity_threshold = 0.d0
      allocate(reaction%kinmnrl_rate_limiter(reaction%nkinmnrl))
      reaction%kinmnrl_rate_limiter = 0.d0
      allocate(reaction%kinmnrl_irreversible(reaction%nkinmnrl))
      reaction%kinmnrl_irreversible = 0
      allocate(reaction%kinmnrl_rate_constant(reaction%nkinmnrl))
      reaction%kinmnrl_rate_constant = 0.d0
      allocate(reaction%kinmnrl_activation_energy(reaction%nkinmnrl))
      reaction%kinmnrl_activation_energy = 0.d0
      allocate(reaction%kinmnrl_molar_vol(reaction%nkinmnrl))
      reaction%kinmnrl_molar_vol = 0.d0
      allocate(reaction%kinmnrl_molar_wt(reaction%nkinmnrl))
      reaction%kinmnrl_molar_wt = 0.d0

      allocate(reaction%kinmnrl_num_prefactors(reaction%nkinmnrl))
      reaction%kinmnrl_num_prefactors = 0
      if (max_num_prefactors > 0) then
        allocate(reaction%kinmnrl_pref_rate(max_num_prefactors,reaction%nkinmnrl))
        reaction%kinmnrl_pref_rate = 0.d0
        allocate(reaction%kinmnrl_pref_activation_energy(max_num_prefactors, &
                                                         reaction%nkinmnrl))
        reaction%kinmnrl_pref_activation_energy = 0.d0
        allocate(reaction%kinmnrl_prefactor_id(0:max_num_prefactor_species, &
                                             max_num_prefactors,reaction%nkinmnrl))
        reaction%kinmnrl_prefactor_id = 0
        allocate(reaction%kinmnrl_pref_alpha(max_num_prefactor_species, &
                                             max_num_prefactors,reaction%nkinmnrl))
        reaction%kinmnrl_pref_alpha = 0.d0
        allocate(reaction%kinmnrl_pref_beta(max_num_prefactor_species, &
                                             max_num_prefactors,reaction%nkinmnrl))
        reaction%kinmnrl_pref_beta = 0.d0
        allocate(reaction%kinmnrl_pref_atten_coef(max_num_prefactor_species, &
                                             max_num_prefactors,reaction%nkinmnrl))
        reaction%kinmnrl_pref_atten_coef = 0.d0
      endif
    endif
    
    cur_mineral => reaction%mineral_list
    imnrl = 1
    ikinmnrl = 1
    do
      if (.not.associated(cur_mineral)) exit

      reaction%mineral_names(imnrl) = cur_mineral%name
      ispec = 0
      do i = 1, cur_mineral%dbaserxn%nspec
        if (cur_mineral%dbaserxn%spec_ids(i) /= h2o_id) then
          ispec = ispec + 1
          spec_id = cur_mineral%dbaserxn%spec_ids(i)
          if (spec_id > h2o_id) spec_id = spec_id - 1
          reaction%mnrlspecid(ispec,imnrl) = spec_id
          reaction%mnrlstoich(ispec,imnrl) = &
            cur_mineral%dbaserxn%stoich(i)
            
        else ! fill in h2o id and stoich
          reaction%mnrlh2oid(imnrl) = h2o_id
          reaction%mnrlh2ostoich(imnrl) = &
            cur_mineral%dbaserxn%stoich(i)
        endif
      enddo
      reaction%mnrlspecid(0,imnrl) = ispec

      call ReactionInitializeLogK_hpt(reaction%mnrl_logKcoef(:,imnrl), &
                                  reaction%mnrl_logK(imnrl), &
                                  option,reaction)

      ! geh - for now, the user must specify they want each individual
      !       mineral printed for non-kinetic reactions (e.g. for SI).
      reaction%mnrl_print(imnrl) = cur_mineral%print_me
      if (cur_mineral%itype == MINERAL_KINETIC) then
        reaction%kinmnrl_names(ikinmnrl) = reaction%mineral_names(imnrl)
        reaction%kinmnrl_print(ikinmnrl) = cur_mineral%print_me .or. &
                                           reaction%mineral%print_all
        reaction%kinmnrlspecid(:,ikinmnrl) = reaction%mnrlspecid(:,imnrl)
        reaction%kinmnrlstoich(:,ikinmnrl) = reaction%mnrlstoich(:,imnrl)
        reaction%kinmnrlh2oid(ikinmnrl) = reaction%mnrlh2oid(imnrl)
        reaction%kinmnrlh2ostoich(ikinmnrl) = reaction%mnrlh2ostoich(imnrl)

        call ReactionInitializeLogK_hpt(reaction%kinmnrl_logKcoef(:,ikinmnrl), &
                                        reaction%kinmnrl_logK(ikinmnrl), &
                                        option,reaction)

        tstrxn => cur_mineral%tstrxn
        if (associated(tstrxn)) then
          ! loop over transition state theory reactions/prefactors
          cur_prefactor => cur_mineral%tstrxn%prefactor
          i = 0
          do
            if (.not.associated(cur_prefactor)) exit
            ! ith prefactor
            i = i + 1

            reaction%kinmnrl_pref_rate(i,ikinmnrl) = cur_prefactor%rate
            reaction%kinmnrl_pref_activation_energy(i,ikinmnrl) = &
              cur_prefactor%activation_energy

            cur_prefactor_species => cur_prefactor%species
            j = 0
            do
              if (.not.associated(cur_prefactor_species)) exit
              ! jth prefactor species
              j = j + 1
              ! find the prefactor species
              do ispec = 1, reaction%naqcomp
                if (StringCompareIgnoreCase(reaction%primary_species_names(ispec), &
                                            cur_prefactor_species%name)) then
                  cur_prefactor_species%id = ispec
                  exit
                endif
              enddo
              if (cur_prefactor_species%id == 0) then ! not found
                ! negative prefactor_species_id denotes a secondary species
                do ispec = 1, reaction%neqcplx
                  if (StringCompareIgnoreCase(reaction%secondary_species_names(ispec), &
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
              reaction%kinmnrl_prefactor_id(j,i,ikinmnrl) = cur_prefactor_species%id
              reaction%kinmnrl_pref_alpha(j,i,ikinmnrl) = cur_prefactor_species%alpha
              reaction%kinmnrl_pref_beta(j,i,ikinmnrl) = cur_prefactor_species%beta
              reaction%kinmnrl_pref_atten_coef(j,i,ikinmnrl) = &
                cur_prefactor_species%attenuation_coef
              cur_prefactor_species => cur_prefactor_species%next
            enddo
            ! store the number of species
            reaction%kinmnrl_prefactor_id(0,i,ikinmnrl) = j
            cur_prefactor => cur_prefactor%next
          enddo
          reaction%kinmnrl_num_prefactors(ikinmnrl) = i

          reaction%kinmnrl_affinity_threshold(ikinmnrl) = &
            tstrxn%affinity_threshold
          reaction%kinmnrl_rate_limiter(ikinmnrl) = tstrxn%rate_limiter
          reaction%kinmnrl_irreversible(ikinmnrl) = tstrxn%irreversible
          if (reaction%kinmnrl_num_prefactors(ikinmnrl) == 0) then
            ! no prefactors, rates stored in upper level
            reaction%kinmnrl_rate_constant(ikinmnrl) = tstrxn%rate
            reaction%kinmnrl_activation_energy(ikinmnrl) = &
              tstrxn%activation_energy
          endif
        endif ! associated(tstrxn)

        reaction%kinmnrl_molar_vol(ikinmnrl) = cur_mineral%molar_volume
        reaction%kinmnrl_molar_wt(ikinmnrl) = cur_mineral%molar_weight
        ikinmnrl = ikinmnrl + 1
      endif

      cur_mineral => cur_mineral%next
      imnrl = imnrl + 1
    enddo
  endif
  
  ! colloids
  reaction%ncoll = GetColloidCount(reaction)

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
  allocate(flags(reaction%naqcomp))
  flags = PETSC_FALSE

  if (reaction%neqsrfcplx > 0) then
  
    ! determine max # complexes for a given site
    icount = 0
    cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_srfcplx_rxn)) exit
      if (cur_srfcplx_rxn%itype == SRFCMPLX_RXN_EQUILIBRIUM .or. &
          cur_srfcplx_rxn%itype == SRFCMPLX_RXN_MULTIRATE_KINETIC) then
        isrfcplx = 0
        cur_srfcplx => cur_srfcplx_rxn%complex_list
        do
          if (.not.associated(cur_srfcplx)) exit
          isrfcplx = isrfcplx + 1
          cur_srfcplx => cur_srfcplx%next
        enddo
        if (isrfcplx > icount) icount = isrfcplx
      endif
      cur_srfcplx_rxn => cur_srfcplx_rxn%next
    enddo
    nullify(cur_srfcplx_rxn)  

    allocate(reaction%eqsrfcplx_rxn_to_surf(reaction%neqsrfcplxrxn))
    reaction%eqsrfcplx_rxn_to_surf = 0
    
    allocate(reaction%eqsrfcplx_rxn_surf_type(reaction%neqsrfcplxrxn))
    reaction%eqsrfcplx_rxn_surf_type = 0
    
    allocate(reaction%srfcplxrxn_to_complex(0:icount, &
                                                 reaction%neqsrfcplxrxn))
    reaction%srfcplxrxn_to_complex = 0
    
    allocate(reaction%eqsrfcplx_site_names(reaction%neqsrfcplxrxn))
    reaction%eqsrfcplx_site_names = ''
    
    allocate(reaction%eqsrfcplx_site_print(reaction%neqsrfcplxrxn))
    reaction%eqsrfcplx_site_print = PETSC_FALSE
    
    allocate(reaction%eqsrfcplx_site_density_print(reaction%neqsrfcplxrxn))
    reaction%eqsrfcplx_site_density_print = PETSC_FALSE
    
    allocate(reaction%eqsrfcplx_rxn_site_density(reaction%neqsrfcplxrxn))
    reaction%eqsrfcplx_rxn_site_density = 0.d0
    
    allocate(reaction%eqsrfcplx_rxn_stoich_flag(reaction%neqsrfcplxrxn))
    reaction%eqsrfcplx_rxn_stoich_flag = PETSC_FALSE
    
    allocate(reaction%eqsrfcplx_names(reaction%neqsrfcplx))
    reaction%eqsrfcplx_names = ''
    
    allocate(reaction%eqsrfcplx_print(reaction%neqsrfcplx))
    reaction%eqsrfcplx_print = PETSC_FALSE
    
    allocate(reaction%srfcplxspecid(0:reaction%naqcomp,reaction%neqsrfcplx))
    reaction%srfcplxspecid = 0
    
    allocate(reaction%eqsrfcplxstoich(reaction%naqcomp,reaction%neqsrfcplx))
    reaction%eqsrfcplxstoich = 0.d0
    
    allocate(reaction%eqsrfcplxh2oid(reaction%neqsrfcplx))
    reaction%eqsrfcplxh2oid = 0
    
    allocate(reaction%eqsrfcplxh2ostoich(reaction%neqsrfcplx))
    reaction%eqsrfcplxh2ostoich = 0.d0
    
    allocate(reaction%eqsrfcplx_free_site_id(reaction%neqsrfcplx))
    reaction%eqsrfcplx_free_site_id = 0
    
    allocate(reaction%eqsrfcplx_free_site_stoich(reaction%neqsrfcplx))
    reaction%eqsrfcplx_free_site_stoich = 0.d0
    
!    allocate(reaction%eqsrfcplx_mineral_id(reaction%neqsrfcplx))
!    reaction%eqsrfcplx_mineral_id = 0
    
    allocate(reaction%eqsrfcplx_logK(reaction%neqsrfcplx))
    reaction%eqsrfcplx_logK = 0.d0
    
    allocate(reaction%eqsrfcplx_logKcoef(reaction%num_dbase_parameters, &
                                           reaction%neqsrfcplx))
    reaction%eqsrfcplx_logKcoef = 0.d0

    allocate(reaction%eqsrfcplx_Z(reaction%neqsrfcplx))
    reaction%eqsrfcplx_Z = 0.d0

    isrfcplx = 0
    irxn = 0
    cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_srfcplx_rxn)) exit
      
      if (cur_srfcplx_rxn%itype == SRFCMPLX_RXN_EQUILIBRIUM .or. &
          cur_srfcplx_rxn%itype == SRFCMPLX_RXN_MULTIRATE_KINETIC) then

        irxn = irxn + 1
        reaction%eqsrfcplx_site_names(irxn) = cur_srfcplx_rxn%free_site_name
        reaction%eqsrfcplx_site_print(irxn) = cur_srfcplx_rxn%free_site_print_me .or. &
                                              reaction%print_all_species
        reaction%eqsrfcplx_site_density_print(irxn) = &
                                              cur_srfcplx_rxn%site_density_print_me .or. &
                                              reaction%print_all_species
        surface_complexation%srfcplxrxn_surf_type(irxn) = &
          cur_srfcplx_rxn%surface_itype
        select case(cur_srfcplx_rxn%surface_itype)
          case(ROCK_SURFACE)
            ! nothing to do here as the linkage to rick density is already set
          case(MINERAL_SURFACE)
            surface_complexation%srfcplxrxn_to_surf(irxn) = &
              GetKineticMineralIDFromName(reaction%mineral, &
                                          cur_srfcplx_rxn%surface_name)
            if (surface_complexation%srfcplxrxn_to_surf(irxn) < 0) then
              option%io_buffer = 'Mineral ' // &
                                  trim(cur_srfcplx_rxn%surface_name) // &
                                  ' listed in surface complexation ' // &
                                  'reaction not found in kinetic mineral list'
              call printErrMsg(option)
            endif
          case(COLLOID_SURFACE)
            surface_complexation%srfcplxrxn_to_surf(irxn) = &
              GetColloidIDFromName(reaction,cur_srfcplx_rxn%surface_name)
            if (surface_complexation%srfcplxrxn_to_surf(irxn) < 0) then
              option%io_buffer = 'Colloid ' // &
                                  trim(cur_srfcplx_rxn%surface_name) // &
                                  ' listed in surface complexation ' // &
                                  'reaction not found in colloid list'
              call printErrMsg(option)
            endif
            ! loop over primary species associated with colloid sorption and
            ! add to colloid species list, if not already listed
            cur_srfcplx_in_rxn => cur_srfcplx_rxn%complex_list
            do
              if (.not.associated(cur_srfcplx_in_rxn)) exit
              ! cur_srfcplx2%ptr is a pointer to complex in master list
              cur_srfcplx => cur_srfcplx_in_rxn%ptr 
              do i = 1, cur_srfcplx%dbaserxn%nspec
                if (cur_srfcplx%dbaserxn%spec_ids(i) == h2o_id) cycle
                spec_id = cur_srfcplx%dbaserxn%spec_ids(i)
                if (spec_id > h2o_id) spec_id = spec_id - 1              
                colloid_species_flag(spec_id) = PETSC_TRUE
              enddo
              nullify(cur_srfcplx)
              cur_srfcplx_in_rxn => cur_srfcplx_in_rxn%next
            enddo
          case(NULL_SURFACE)
            write(word,*) cur_srfcplx_rxn%id
            option%io_buffer = 'No mineral or colloid name specified ' // &
              'for equilibrium surface complexation reaction:' // &
              trim(adjustl(word))
            call printWrnMsg(option)
        end select
        reaction%eqsrfcplx_rxn_site_density(irxn) = cur_srfcplx_rxn%site_density
              
        cur_srfcplx => cur_srfcplx_rxn%complex_list
        do
          if (.not.associated(cur_srfcplx)) exit
          
          isrfcplx = isrfcplx + 1
          
          ! set up integer pointers from site to complexes
          ! increment count for site
          reaction%srfcplxrxn_to_complex(0,irxn) = &
            reaction%srfcplxrxn_to_complex(0,irxn) + 1
          reaction%srfcplxrxn_to_complex( &
            reaction%srfcplxrxn_to_complex(0,irxn),irxn) = isrfcplx 
          
          reaction%eqsrfcplx_names(isrfcplx) = cur_srfcplx%name
          reaction%eqsrfcplx_print(isrfcplx) = cur_srfcplx%print_me .or. &
                                              reaction%print_all_species
          reaction%eqsrfcplx_free_site_id(isrfcplx) = &
            cur_srfcplx_rxn%free_site_id
          reaction%eqsrfcplx_free_site_stoich(isrfcplx) =  &
            cur_srfcplx%free_site_stoich
            
          if (cur_srfcplx%free_site_stoich > 1.d0) then
            reaction%eqsrfcplx_rxn_stoich_flag(irxn) = PETSC_TRUE
          endif
   
          ispec = 0
          do i = 1, cur_srfcplx%dbaserxn%nspec
            if (cur_srfcplx%dbaserxn%spec_ids(i) /= h2o_id) then
              ispec = ispec + 1
              spec_id = cur_srfcplx%dbaserxn%spec_ids(i)
              if (spec_id > h2o_id) spec_id = spec_id - 1
              reaction%srfcplxspecid(ispec,isrfcplx) = spec_id
              reaction%eqsrfcplxstoich(ispec,isrfcplx) = &
                cur_srfcplx%dbaserxn%stoich(i)
              
            else ! fill in h2o id and stoich
              reaction%eqsrfcplxh2oid(isrfcplx) = h2o_id
              reaction%eqsrfcplxh2ostoich(isrfcplx) = &
                cur_srfcplx%dbaserxn%stoich(i)
            endif
          enddo
          reaction%srfcplxspecid(0,isrfcplx) = ispec
          call ReactionInitializeLogK_hpt(reaction%eqsrfcplx_logKcoef(:,isrfcplx), &
                                          reaction%eqsrfcplx_logK(isrfcplx), &
                                          option,reaction)

          reaction%eqsrfcplx_Z(isrfcplx) = cur_srfcplx%Z

          cur_srfcplx => cur_srfcplx%next
        enddo
        nullify(cur_srfcplx)
        
      endif
      cur_srfcplx_rxn => cur_srfcplx_rxn%next
    enddo
    nullify(cur_srfcplx_rxn)  
  
  endif
  
  if (reaction%nkinsrfcplxrxn > 0) then
  
    ! determine max # complexes for a given site
    icount = 0
    cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_srfcplx_rxn)) exit
      if (cur_srfcplx_rxn%itype == SRFCMPLX_RXN_KINETIC) then
        isrfcplx = 0
        cur_srfcplx => cur_srfcplx_rxn%complex_list
        do
          if (.not.associated(cur_srfcplx)) exit
          isrfcplx = isrfcplx + 1
          cur_srfcplx => cur_srfcplx%next
        enddo
        if (isrfcplx > icount) icount = isrfcplx
      endif
      cur_srfcplx_rxn => cur_srfcplx_rxn%next
    enddo
    nullify(cur_srfcplx_rxn) 
    allocate(reaction%kinsrfcplx_rxn_to_complex(0:icount, &
                                                 reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_rxn_to_complex = 0
    
    allocate(reaction%kinsrfcplx_rxn_to_site(reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_rxn_to_site = 0
    
    allocate(reaction%kinsrfcplx_rxn_to_surf(reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_rxn_to_surf = 0
    
    allocate(reaction%kinsrfcplx_rxn_surf_type(reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_rxn_surf_type = 0
    
    allocate(reaction%kinsrfcplx_site_names(reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_site_names = ''
    
    allocate(reaction%kinsrfcplx_site_print(reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_site_print = PETSC_FALSE
    
    allocate(reaction%kinsrfcplx_rxn_site_density(reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_rxn_site_density = 0.d0
    
    allocate(reaction%kinsrfcplx_rxn_stoich_flag(reaction%nkinsrfcplxrxn))
    reaction%kinsrfcplx_rxn_stoich_flag = PETSC_FALSE
    
    allocate(reaction%kinsrfcplx_names(reaction%nkinsrfcplx))
    reaction%kinsrfcplx_names = ''
    
    allocate(reaction%kinsrfcplx_print(reaction%nkinsrfcplx))
    reaction%kinsrfcplx_print = PETSC_FALSE
    
    allocate(reaction%kinsrfcplxspecid(0:reaction%naqcomp,reaction%nkinsrfcplx))
    reaction%kinsrfcplxspecid = 0
    
    allocate(reaction%kinsrfcplxstoich(reaction%naqcomp,reaction%nkinsrfcplx))
    reaction%kinsrfcplxstoich = 0.d0
    
    allocate(reaction%kinsrfcplxh2oid(reaction%nkinsrfcplx))
    reaction%kinsrfcplxh2oid = 0
    
    allocate(reaction%kinsrfcplxh2ostoich(reaction%nkinsrfcplx))
    reaction%kinsrfcplxh2ostoich = 0.d0
    
    allocate(reaction%kinsrfcplx_free_site_id(reaction%nkinsrfcplx))
    reaction%kinsrfcplx_free_site_id = 0
    
    allocate(reaction%kinsrfcplx_free_site_stoich(reaction%nkinsrfcplx))
    reaction%kinsrfcplx_free_site_stoich = 0.d0
    
    allocate(reaction%kinsrfcplx_forward_rate(reaction%nkinsrfcplx))
    reaction%kinsrfcplx_forward_rate = 0.d0
    
    allocate(reaction%kinsrfcplx_backward_rate(reaction%nkinsrfcplx))
    reaction%kinsrfcplx_backward_rate = 0.d0

!    allocate(reaction%kinsrfcplx_logK(reaction%nkinsrfcplx))
!    reaction%kinsrfcplx_logK = 0.d0
!#if TEMP_DEPENDENT_LOGK
!    allocate(reaction%kinsrfcplx_logKcoef(FIVE_INTEGER,reaction%nkinsrfcplx))
!    reaction%kinsrfcplx_logKcoef = 0.d0
!#else
!    allocate(reaction%kinsrfcplx_logKcoef(reaction%num_dbase_temperatures, &
!                                           reaction%nkinsrfcplx))
!    reaction%kinsrfcplx_logKcoef = 0.d0
!#endif
    allocate(reaction%kinsrfcplx_Z(reaction%nkinsrfcplx))
    reaction%kinsrfcplx_Z = 0.d0

    isrfcplx = 0
    irxn = 0
    cur_srfcplx_rxn => reaction%surface_complexation_rxn_list
    do
      if (.not.associated(cur_srfcplx_rxn)) exit
      
      if (cur_srfcplx_rxn%itype == SRFCMPLX_RXN_KINETIC) then

        irxn = irxn + 1
        
        reaction%kinsrfcplx_site_names(irxn) = cur_srfcplx_rxn%free_site_name
        reaction%kinsrfcplx_site_print(irxn) = cur_srfcplx_rxn%free_site_print_me .or. &
                                              reaction%print_all_species
        if (len_trim(cur_srfcplx_rxn%mineral_name) > 1) then
          reaction%kinsrfcplx_rxn_surf_type(irxn) = MINERAL_SURFACE
          reaction%kinsrfcplx_rxn_to_surf(irxn) = &
            GetKineticMineralIDFromName(reaction,cur_srfcplx_rxn%mineral_name)
          if (reaction%kinsrfcplx_rxn_to_surf(irxn) < 0) then
            option%io_buffer = 'Mineral ' // trim(cur_srfcplx_rxn%mineral_name) // &
                               'listed in kinetic surface complexation ' // &
                               'reaction not found in mineral list'
            call printErrMsg(option)
          endif
        else if (len_trim(cur_srfcplx_rxn%colloid_name) > 1) then
          reaction%kinsrfcplx_rxn_surf_type(irxn) = COLLOID_SURFACE
          reaction%kinsrfcplx_rxn_to_surf(irxn) = &
            GetColloidIDFromName(reaction,cur_srfcplx_rxn%colloid_name)
          if (reaction%kinsrfcplx_rxn_to_surf(irxn) < 0) then
            option%io_buffer = 'Colloid ' // trim(cur_srfcplx_rxn%colloid_name) // &
                               'listed in kinetic surface complexation ' // &
                               'reaction not found in colloid list'
            call printErrMsg(option)
          endif
          ! loop over primary species associated with colloid sorption and
          ! add to colloid species list, if not already listed
          cur_srfcplx => cur_srfcplx_rxn%complex_list
          do
            if (.not.associated(cur_srfcplx)) exit
            do i = 1, cur_srfcplx%dbaserxn%nspec
              if (cur_srfcplx%dbaserxn%spec_ids(i) == h2o_id) cycle
              spec_id = cur_srfcplx%dbaserxn%spec_ids(i)
              if (spec_id > h2o_id) spec_id = spec_id - 1              
              flags(spec_id) = PETSC_TRUE
            enddo
            cur_srfcplx => cur_srfcplx%next
          enddo

        else
          write(word,*) cur_srfcplx_rxn%id
          write(word,*) cur_srfcplx_rxn%id
          option%io_buffer = 'No mineral or colloid name specified for ' // &
            'kinetic surface complexation reaction:' // &
            trim(adjustl(word))
          call printWrnMsg(option)
          reaction%kinsrfcplx_rxn_surf_type(irxn) = NULL_SURFACE          
        endif
        reaction%kinsrfcplx_rxn_site_density(irxn) = cur_srfcplx_rxn%site_density
              
        cur_srfcplx => cur_srfcplx_rxn%complex_list
        do
          if (.not.associated(cur_srfcplx)) exit
          
          isrfcplx = isrfcplx + 1

          reaction%kinsrfcplx_forward_rate(isrfcplx) = cur_srfcplx%forward_rate
          if (Uninitialized(cur_srfcplx%backward_rate)) then
            ! backward rate will be calculated based on Kb = Kf * Keq
            call Interpolate(temp_high,temp_low,option%reference_temperature, &
                             cur_srfcplx%dbaserxn%logK(itemp_high), &
                             cur_srfcplx%dbaserxn%logK(itemp_low), &
                             value)
            reaction%kinsrfcplx_backward_rate(isrfcplx) = 10.d0**value * &
                                                          cur_srfcplx%forward_rate
          else
            reaction%kinsrfcplx_backward_rate(isrfcplx) = cur_srfcplx%backward_rate
          endif
          ! set up integer pointers from site to complexes
          ! increment count for site
          reaction%kinsrfcplx_rxn_to_complex(0,irxn) = &
            reaction%kinsrfcplx_rxn_to_complex(0,irxn) + 1
          reaction%kinsrfcplx_rxn_to_complex( &
            reaction%kinsrfcplx_rxn_to_complex(0,irxn),irxn) = isrfcplx 
          reaction%kinsrfcplx_rxn_to_site(irxn) = cur_srfcplx_rxn%free_site_id

          
          reaction%kinsrfcplx_names(isrfcplx) = cur_srfcplx%name
          reaction%kinsrfcplx_print(isrfcplx) = cur_srfcplx%print_me .or. &
                                              reaction%print_all_species
          reaction%kinsrfcplx_free_site_stoich(isrfcplx) =  &
            cur_srfcplx%free_site_stoich
            
          if (cur_srfcplx%free_site_stoich > 1.d0) then
            reaction%kinsrfcplx_rxn_stoich_flag(irxn) = PETSC_TRUE
          endif
   
          ispec = 0
          do i = 1, cur_srfcplx%dbaserxn%nspec
            if (cur_srfcplx%dbaserxn%spec_ids(i) /= h2o_id) then
              ispec = ispec + 1
              spec_id = cur_srfcplx%dbaserxn%spec_ids(i)
              if (spec_id > h2o_id) spec_id = spec_id - 1
              reaction%kinsrfcplxspecid(ispec,isrfcplx) = spec_id
              reaction%kinsrfcplxstoich(ispec,isrfcplx) = &
                cur_srfcplx%dbaserxn%stoich(i)
            else ! fill in h2o id and stoich
              reaction%kinsrfcplxh2oid(isrfcplx) = h2o_id
              reaction%kinsrfcplxh2ostoich(isrfcplx) = &
                cur_srfcplx%dbaserxn%stoich(i)
            endif
          enddo
          reaction%kinsrfcplxspecid(0,isrfcplx) = ispec
!  Chuan added, for surface complex reaction, the new data base 
          cur_srfcplx%dbaserxn%logK(:) = cur_srfcplx%dbaserxn%logKCoeff_hpt(:)
!  #if TEMP_DEPENDENT_LOGK
!        call ReactionFitLogKCoef(reaction%kinsrfcplx_logKcoef(:,isrfcplx),cur_srfcplx%dbaserxn%logK, &
!                                 reaction%kinsrfcplx_names(isrfcplx), &
!                                 option,reaction)
!        call ReactionInitializeLogK(reaction%kinsrfcplx_logKcoef(:,isrfcplx), &
!                                    cur_srfcplx%dbaserxn%logK, &
!                                    reaction%kinsrfcplx_logK(isrfcplx), &
!                                    option,reaction)
!  #else
!          call Interpolate(temp_high,temp_low,option%reference_temperature, &
!                           cur_srfcplx%dbaserxn%logK(itemp_high), &
!                           cur_srfcplx%dbaserxn%logK(itemp_low), &
!                           reaction%kinsrfcplx_logK(isrfcplx))
!          !reaction%kinsrfcplx_logK(isrfcplx) = cur_srfcplx%dbaserxn%logK(option%itemp_ref)
!  #endif
          reaction%kinsrfcplx_Z(isrfcplx) = cur_srfcplx%Z

          cur_srfcplx => cur_srfcplx%next
        enddo
        nullify(cur_srfcplx)
        
      endif
      
      cur_srfcplx_rxn => cur_srfcplx_rxn%next
    enddo
    nullify(cur_srfcplx_rxn)  
  
  endif


  ! allocate colloids species names, mappings, etc.
  reaction%ncollcomp = 0
  icount = 0
  do i = 1, reaction%naqcomp
    if (flags(i)) then 
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
      if (flags(i)) then
        icount = icount + 1
        reaction%colloid_species_names(icount) = &
          trim(reaction%primary_species_names(i))
        reaction%coll_spec_to_pri_spec(icount) = i
        reaction%pri_spec_to_coll_spec(i) = icount
      endif
    enddo
    if (minval(reaction%coll_spec_to_pri_spec) < 1) then
      option%io_buffer = 'Species colloid surface complexation reaction not' // &
                         ' recognized among primary species'
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
  deallocate(flags)

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
                   ' in ion exchange reaction' // &
                   ' not found in swapped basis.'
          call printErrMsg(option)     
        endif
        cur_cation => cur_cation%next
      enddo
      reaction%eqionx_rxn_cationid(0,irxn) = ication
      ! Find any Zi /= Zj for all species i, j
      found = PETSC_FALSE
      do i = 1, reaction%eqionx_rxn_cationid(0,irxn)
        do j = 1, reaction%eqionx_rxn_cationid(0,irxn)
          if (abs(reaction%primary_spec_Z(reaction%eqionx_rxn_cationid(i,irxn))- &
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
  
  ! general reaction
  
  if (reaction%ngeneral_rxn > 0) then
  
    ! process reaction equation into the database format
    cur_general_rxn => reaction%general_rxn_list
    do
      if (.not.associated(cur_general_rxn)) exit
      
      ! count # species
      icount = 0
      string = cur_general_rxn%reaction
      do
        ierr = 0
        call InputReadWord(string,word,PETSC_TRUE,ierr)
        if (InputError(ierr)) exit

        select case(word)
          case('+')
          case('-')
          case('=','<=>','<->')
          case default
          ! try reading as double precision
          string2 = word
          if (.not.StringStartsWithAlpha(string2)) then
            ! the word is the stoichiometry value
          else
            ! the word is the species name
            icount = icount + 1
          endif
        end select

      enddo
      
      ! load species into database format
      
      cur_general_rxn%dbaserxn => DatabaseRxnCreate()
      
      dbaserxn => DatabaseRxnCreate()
      dbaserxn%nspec = icount
      allocate(dbaserxn%spec_name(icount))
      dbaserxn%spec_name = ''
      allocate(dbaserxn%stoich(icount))
      dbaserxn%stoich = UNINITIALIZED_DOUBLE
      allocate(dbaserxn%spec_ids(icount))
      dbaserxn%spec_ids = 0

      string = cur_general_rxn%reaction
      icount = 1
      ! midpoint points to the first product species, as in
      ! reactant1 + reactant2 <-> product1 + product2
      midpoint = 0
      negative_flag = PETSC_FALSE
      do
        ierr = 0
        call InputReadWord(string,word,PETSC_TRUE,ierr)
        if (InputError(ierr)) exit

        select case(word)
          case('+')
          case('-')
            ! toggle negative flag
            if (negative_flag) then
              negative_flag = PETSC_FALSE
            else
              negative_flag = PETSC_TRUE
            endif
          case('=','<=>','<->')
            midpoint = icount
          case default
            ! try reading as double precision
            string2 = word
            if (.not.StringStartsWithAlpha(string2)) then
              ! negate if a product
              call InputReadDouble(string2,option,value,ierr)
              ! negate if negative stoichiometry
              if (negative_flag) value = -1.0*value
              dbaserxn%stoich(icount) = value
            else
              dbaserxn%spec_name(icount) = word
              if (negative_flag .and. &
                  (dbaserxn%stoich(icount) + 999.d0) < 1.d-10) then
                dbaserxn%stoich(icount) = -1.d0
              endif

              ! set the primary species id
              found = PETSC_FALSE
              do i = 1, reaction%naqcomp
                if (StringCompare(word, &
                                  reaction%primary_species_names(i), &
                                  MAXWORDLENGTH)) then
                  dbaserxn%spec_ids(icount) = i
                  found = PETSC_TRUE
                  exit      
                endif
              enddo
              ! check water
              word2 = 'H2O'
              if (StringCompareIgnoreCase(word,word2)) then
                ! don't increment icount
                exit
              endif              
              if (.not.found) then
                option%io_buffer = 'Species ' // trim(word) // &
                         ' in general reaction' // &
                         ' not found among primary species list.'
                call printErrMsg(option)     
              endif
              icount = icount + 1
            endif
            negative_flag = PETSC_FALSE
        end select

      enddo
      
      ! if no stoichiometry specified, default = 1.
      do i = 1, dbaserxn%nspec
        if ((dbaserxn%stoich(i) + 999.d0) < 1.d-10) dbaserxn%stoich(i) = 1.d0
      enddo
      ! negate stoichiometries after midpoint
      do i = midpoint, dbaserxn%nspec
        dbaserxn%stoich(i) = -1.d0*dbaserxn%stoich(i)
      enddo
      ! now negate all stoichiometries to have - for reactants; + for products
      do i = 1, dbaserxn%nspec
        dbaserxn%stoich(i) = -1.d0*dbaserxn%stoich(i)
      enddo
      ! reorder species ids in ascending order
      do i = 1, dbaserxn%nspec
        do j = i+1, dbaserxn%nspec
          if (dbaserxn%spec_ids(i) > dbaserxn%spec_ids(j)) then
            ! swap ids
            idum = dbaserxn%spec_ids(j)
            dbaserxn%spec_ids(j) = dbaserxn%spec_ids(i)
            dbaserxn%spec_ids(i) = idum
            ! swap stoichiometry
            value = dbaserxn%stoich(j)
            dbaserxn%stoich(j) = dbaserxn%stoich(i)
            dbaserxn%stoich(i) = value
            ! swap names
            word = dbaserxn%spec_name(j)
            dbaserxn%spec_name(j) = dbaserxn%spec_name(i)
            dbaserxn%spec_name(i) = word
          endif
        enddo
      enddo
      
      cur_general_rxn%dbaserxn => dbaserxn
      
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
      if (backward_count > max_backward_count) max_backward_count = backward_count
      if (species_count > max_species_count) max_species_count = species_count

      cur_general_rxn => cur_general_rxn%next

    enddo
    nullify(cur_general_rxn)
    
    allocate(reaction%generalspecid(0:max_species_count,reaction%ngeneral_rxn))
    reaction%generalspecid = 0
    allocate(reaction%generalstoich(max_species_count,reaction%ngeneral_rxn))
    reaction%generalstoich = 0.d0
    allocate(reaction%generalforwardspecid(0:max_forward_count,reaction%ngeneral_rxn))
    reaction%generalforwardspecid = 0
    allocate(reaction%generalforwardstoich(max_forward_count,reaction%ngeneral_rxn))
    reaction%generalforwardstoich = 0.d0
    allocate(reaction%generalbackwardspecid(0:max_backward_count,reaction%ngeneral_rxn))
    reaction%generalbackwardspecid = 0
    allocate(reaction%generalbackwardstoich(max_backward_count,reaction%ngeneral_rxn))
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
          reaction%generalforwardspecid(forward_count,irxn) = dbaserxn%spec_ids(i)
          ! ensure that forward stoich is positive for rate expression
          reaction%generalforwardstoich(forward_count,irxn) = dabs(dbaserxn%stoich(i))
        else if (dbaserxn%stoich(i) > 0.d0) then
          backward_count = backward_count + 1
          reaction%generalbackwardspecid(backward_count,irxn) = dbaserxn%spec_ids(i)
          reaction%generalbackwardstoich(backward_count,irxn) = dbaserxn%stoich(i)
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
  
  ! Kd reactions
  
  if (reaction%neqkdrxn > 0) then
  
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

    cur_kd_rxn => reaction%kd_rxn_list
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
        option%io_buffer = 'Species ' // trim(word) // &
                 ' in kd reaction' // &
                 ' not found among primary species list.'
        call printErrMsg(option)     
      endif
      reaction%eqkdtype(irxn) = cur_kd_rxn%itype
      reaction%eqkddistcoef(irxn) = cur_kd_rxn%Kd
      reaction%eqkdlangmuirb(irxn) = cur_kd_rxn%Langmuir_b
      reaction%eqkdfreundlichn(irxn) = cur_kd_rxn%Freundlich_n
       
      cur_kd_rxn => cur_kd_rxn%next
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

  do ispec = 1, reaction%ngas
    if (reaction%species_idx%o2_gas_id == 0) then
      word = 'O2(g)'
      if (StringCompareIgnoreCase(reaction%gas_species_names(ispec), &
                                  word)) then
        reaction%species_idx%o2_gas_id = ispec
      endif
    endif
    if (reaction%species_idx%co2_gas_id == 0) then
      word = 'CO2(g)'
      if (StringCompareIgnoreCase(reaction%gas_species_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_gas_id = ispec
      endif
      word = 'CO2(g)*'
      if (StringCompareIgnoreCase(reaction%gas_species_names(ispec), &
                                  word)) then
        reaction%species_idx%co2_gas_id = ispec
      endif

    endif

  enddo
  
90 format(80('-'))
100 format(/,2x,i4,2x,a)
110 format(100(/,14x,3(a20,2x)))

  if (OptionPrintToFile(option)) then
    write(option%fid_out,90)
    write(option%fid_out,100) reaction%naqcomp, 'Primary Species'
    write(option%fid_out,110) (reaction%primary_species_names(i),i=1,reaction%naqcomp)
    
    write(option%fid_out,100) reaction%neqcplx, 'Secondary Complex Species'
    write(option%fid_out,110) (reaction%secondary_species_names(i),i=1,reaction%neqcplx)
    
    write(option%fid_out,100) reaction%ngas, 'Gas Species'
    write(option%fid_out,110) (reaction%gas_species_names(i),i=1,reaction%ngas)
    
    write(option%fid_out,100) reaction%nmnrl, 'Reference Minerals'
    write(option%fid_out,110) (reaction%mineral_names(i),i=1,reaction%nmnrl)
    
    write(option%fid_out,100) reaction%nkinmnrl, 'Kinetic Mineral Reactions'
    write(option%fid_out,110) (reaction%kinmnrl_names(i),i=1,reaction%nkinmnrl)
    
    write(option%fid_out,100) reaction%neqsrfcplxrxn, 'Surface Complexation Reactions'
    write(option%fid_out,110) (reaction%eqsrfcplx_site_names(i),i=1,reaction%neqsrfcplxrxn)
    write(option%fid_out,100) reaction%neqsrfcplx, 'Surface Complexes'
    write(option%fid_out,110) (reaction%eqsrfcplx_names(i),i=1,reaction%neqsrfcplx)
    
    write(option%fid_out,100) reaction%neqionxrxn, 'Ion Exchange Reactions'
    write(option%fid_out,100) reaction%neqionxcation, 'Ion Exchange Cations'
    write(option%fid_out,90)
  endif
  
#if 0
  ! output for ASCEM reactions
  if (OptionPrintToFile(option)) then
    open(unit=86,file='reaction.dat')
    write(86,'(10i4)') reaction%naqcomp, reaction%neqcplx, reaction%ngeneral_rxn, & 
                       reaction%neqsrfcplxrxn, reaction%nkinmnrl
    do icomp = 1, reaction%naqcomp
      write(86,'(a12,f6.2,f6.2)') reaction%primary_species_names(icomp), &
                                  reaction%primary_spec_Z(icomp), &
                                  reaction%primary_spec_a0(icomp)
    enddo
    do icplx = 1, reaction%neqcplx
      write(86,'(a32,f6.2,f6.2)') reaction%secondary_species_names(icplx), &
                                  reaction%eqcplx_Z(icplx), &
                                  reaction%eqcplx_a0(icplx)
      write(86,'(40i4)') reaction%eqcplxspecid(:,icplx)
      write(86,'(40f6.2)') reaction%eqcplxstoich(:,icplx)
      write(86,'(i4)') reaction%eqcplxh2oid(icplx)
      write(86,'(f6.2)') reaction%eqcplxh2ostoich(icplx)
      write(86,'(1es13.5)') reaction%eqcplx_logK(icplx)
    enddo
    do irxn = 1, reaction%ngeneral_rxn
      write(86,'(40i4)') reaction%generalspecid(:,irxn)
      write(86,'(40f6.2)') reaction%generalstoich(:,irxn)
      write(86,'(40i4)') reaction%generalforwardspecid(:,irxn)
      write(86,'(40f6.2)') reaction%generalforwardstoich(:,irxn)
      write(86,'(40i4)') reaction%generalbackwardspecid(:,irxn)
      write(86,'(40f6.2)') reaction%generalbackwardstoich(:,irxn)
      write(86,'(f6.2)') reaction%generalh2ostoich(irxn)
      write(86,'(1es13.5)') reaction%general_kf(irxn)
      write(86,'(1es13.5)') reaction%general_kr(irxn)
    enddo
    do irxn = 1, reaction%neqsrfcplxrxn
      write(86,'(a32)')reaction%eqsrfcplx_site_names(irxn)
      write(86,'(1es13.5)') reaction%eqsrfcplx_rxn_site_density(irxn)
      write(86,'(i4)') reaction%srfcplxrxn_to_complex(0,irxn) ! # complexes
      do i = 1, reaction%srfcplxrxn_to_complex(0,irxn)
        icplx = reaction%srfcplxrxn_to_complex(i,irxn)
        write(86,'(a32,f6.2)') reaction%eqsrfcplx_names(icplx), &
                               reaction%eqsrfcplx_Z(icplx)
        write(86,'(40i4)') reaction%srfcplxspecid(:,icplx)
        write(86,'(40f6.2)') reaction%eqsrfcplxstoich(:,icplx)
        write(86,'(i4)') reaction%eqsrfcplxh2oid(icplx)
        write(86,'(f6.2)') reaction%eqsrfcplxh2ostoich(icplx)
        write(86,'(i4)') reaction%eqsrfcplx_free_site_id(icplx)
        write(86,'(f6.2)') reaction%eqsrfcplx_free_site_stoich(icplx)
        write(86,'(1es13.5)') reaction%eqsrfcplx_logK(icplx)

      enddo
    enddo
    do imnrl = 1, reaction%nkinmnrl
      write(86,'(a32)') reaction%kinmnrl_names(imnrl)
      write(86,'(40i4)') reaction%kinmnrlspecid(:,imnrl)
      write(86,'(40f6.2)') reaction%kinmnrlstoich(:,imnrl)
      write(86,'(i4)') reaction%kinmnrlh2oid(imnrl)
      write(86,'(f6.2)') reaction%kinmnrlh2ostoich(imnrl)
      write(86,'(1es13.5)') reaction%kinmnrl_logK(imnrl)
      write(86,'(1es13.5)') reaction%kinmnrl_molar_vol(imnrl)
      write(86,'(1es13.5)') reaction%kinmnrl_molar_wt(imnrl)
      write(86,'(1es13.5)') reaction%kinmnrl_rate_constant(1,imnrl)
      write(86,'(1es13.5)') 1.d0 ! specific surface area 1 cm^2 / cm^3
    enddo
        close(86)
  endif
#endif  
  
  if (allocated(new_basis)) deallocate(new_basis)
  if (allocated(old_basis)) deallocate(old_basis)
  if (allocated(transformation)) deallocate(transformation)
  if (allocated(stoich_prev)) deallocate(stoich_prev)
  if (allocated(stoich_new)) deallocate(stoich_new)
  if (allocated(logKCoeffvector)) deallocate(logKCoeffvector)
  if (allocated(indices)) deallocate(indices)

  if (allocated(new_basis_names)) deallocate(new_basis_names)
  if (allocated(old_basis_names)) deallocate(old_basis_names)
!TODO(geh)
#endif  
end subroutine BasisInit_hpt

! ************************************************************************** !

subroutine BasisSubSpecInGasOrSecRxn_hpt(name1,dbaserxn1,dbaserxn2)
  ! 
  ! Swaps out a chemical species in a chemical
  ! reaction, replacing it with the species in a
  ! secondary reaction (swaps 1 into 2)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/2012
  ! 

  use String_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: name1
  type(database_rxn_type) :: dbaserxn1
  type(database_rxn_type) :: dbaserxn2

  PetscReal :: scale
  
  call BasisSubSpeciesInGasOrSecRxn(name1,dbaserxn1,dbaserxn2,scale)
  dbaserxn2%logKCoeff_hpt = dbaserxn2%logKCoeff_hpt + &
    scale*dbaserxn1%logKCoeff_hpt

  end subroutine BasisSubSpecInGasOrSecRxn_hpt

! ************************************************************************** !

subroutine BasisSubSpeciesInMineralRxn_hpt(name,sec_dbaserxn,mnrl_dbaserxn)
  ! 
  ! Swaps out a chemical species in a chemical
  ! reaction, replacing it with the species in a
  ! secondary reaction (swaps 1 into 2)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/2012
  ! 

  use String_module
  use Reaction_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(database_rxn_type) :: sec_dbaserxn
  type(database_rxn_type) :: mnrl_dbaserxn
  
  PetscReal :: scale

  call BasisSubSpeciesInMineralRxn(name,sec_dbaserxn,mnrl_dbaserxn,scale)
  mnrl_dbaserxn%logKCoeff_hpt = mnrl_dbaserxn%logKCoeff_hpt + &
    scale*sec_dbaserxn%logKCoeff_hpt

end subroutine BasisSubSpeciesInMineralRxn_hpt

end module Reaction_Database_hpt_module
