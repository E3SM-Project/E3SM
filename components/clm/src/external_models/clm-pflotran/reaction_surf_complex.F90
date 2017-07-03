module Reaction_Surface_Complexation_module

  use Reaction_Surface_Complexation_Aux_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: SurfaceComplexationRead, &
            SrfCplxProcessConstraint, &
            RTotalSorbEqSurfCplx, &
            RMultiRateSorption, &
            RKineticSurfCplx, &
            RTotalSorbMultiRateAsEQ
            
contains

! ************************************************************************** !

subroutine SurfaceComplexationRead(reaction,input,option)
  ! 
  ! Reads chemical species
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/02/08
  ! 

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
  type(surface_complexation_type), pointer :: surface_complexation
  type(surface_complex_type), pointer :: srfcplx, cur_srfcplx, prev_srfcplx, &
                                         cur_srfcplx_in_rxn
  type(surface_complex_type), pointer :: rate_list, cur_srfcplx_rate, &
                                         prev_srfcplx_rate  
  type(surface_complexation_rxn_type), pointer :: srfcplx_rxn, cur_srfcplx_rxn
  PetscInt :: temp_srfcplx_count
  PetscBool :: found
  PetscReal :: tempreal
  PetscInt :: i
  PetscInt :: num_times_surface_type_set
  
  nullify(srfcplx_rxn)
  nullify(cur_srfcplx_rxn)
  nullify(cur_srfcplx)
  
  surface_complexation => reaction%surface_complexation
           
  srfcplx_rxn => SurfaceComplexationRxnCreate()
  ! default
  srfcplx_rxn%itype = SRFCMPLX_RXN_EQUILIBRIUM
  temp_srfcplx_count = 0
  num_times_surface_type_set = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                        'CHEMISTRY,SURFACE_COMPLEXATION_RXN')
    call StringToUpper(word)
                
    select case(trim(word))
      case('EQUILIBRIUM')
        srfcplx_rxn%itype = SRFCMPLX_RXN_EQUILIBRIUM
      case('MULTIRATE_KINETIC')
        srfcplx_rxn%itype = SRFCMPLX_RXN_MULTIRATE_KINETIC
      case('KINETIC')
        srfcplx_rxn%itype = SRFCMPLX_RXN_KINETIC
      case('COMPLEX_KINETICS')
        nullify(prev_srfcplx)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
                      
          srfcplx => SurfaceComplexCreate()
          call InputReadWord(input,option,srfcplx%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE')
                        
          do
            call InputReadPflotranString(input,option)
            call InputReadStringErrorMsg(input,option,card)
            if (InputCheckExit(input,option)) exit
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'word', &
                    'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE') 
            select case(trim(word))
              case('FORWARD_RATE_CONSTANT')
                call InputReadDouble(input,option,srfcplx%forward_rate)
                call InputErrorMsg(input,option,'forward_rate', &
                        'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE')
              case('BACKWARD_RATE_CONSTANT')
                call InputReadDouble(input,option,srfcplx%backward_rate)
                call InputErrorMsg(input,option,'backward_rate', &
                        'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE')
              case default
                call InputKeywordUnrecognized(word, &
                       'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_KINETIC_RATE',option)
            end select
          enddo
                                      
          if (.not.associated(rate_list)) then
            rate_list => srfcplx
          endif
          if (associated(prev_srfcplx)) then
            prev_srfcplx%next => srfcplx
          endif
          prev_srfcplx => srfcplx
          nullify(srfcplx)
        enddo
        nullify(prev_srfcplx)
      case('RATE','RATES') 
        srfcplx_rxn%itype = SRFCMPLX_RXN_MULTIRATE_KINETIC
        string = 'RATES inside SURFACE_COMPLEXATION_RXN'
        call UtilityReadArray(srfcplx_rxn%rates,NEG_ONE_INTEGER,string,input, &
                              option) 
      case('SITE_FRACTION') 
        string = 'SITE_FRACTION inside SURFACE_COMPLEXATION_RXN'
        call UtilityReadArray(srfcplx_rxn%site_fractions,NEG_ONE_INTEGER, &
                              string,input,option) 
      case('MULTIRATE_SCALE_FACTOR')
        call InputReadDouble(input,option,srfcplx_rxn%kinmr_scale_factor)
        call InputErrorMsg(input,option,'keyword', &
          'CHEMISTRY,SURFACE_COMPLEXATION_RXN,MULTIRATE_SCALE_FACTOR')
      case('MINERAL')
        srfcplx_rxn%surface_itype = MINERAL_SURFACE
        num_times_surface_type_set = num_times_surface_type_set + 1
        call InputReadWord(input,option,srfcplx_rxn%surface_name, &
          PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword', &
          'CHEMISTRY,SURFACE_COMPLEXATION_RXN,MINERAL_NAME')
      case('ROCK_DENSITY')
        srfcplx_rxn%surface_itype = ROCK_SURFACE
        num_times_surface_type_set = num_times_surface_type_set + 1
      case('COLLOID')
        srfcplx_rxn%surface_itype = COLLOID_SURFACE
        num_times_surface_type_set = num_times_surface_type_set + 1
        call InputReadWord(input,option,srfcplx_rxn%surface_name, &
          PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword', &
          'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COLLOID_NAME')
      case('SITE')
        call InputReadWord(input,option,srfcplx_rxn%free_site_name, &
          PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword', &
          'CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_NAME')
        ! site density in mol/m^3 bulk
        call InputReadDouble(input,option,srfcplx_rxn%site_density)
        call InputErrorMsg(input,option,'keyword', &
          'CHEMISTRY,SURFACE_COMPLEXATION_RXN,SITE_DENSITY')                   
      case('COMPLEXES')
        nullify(prev_srfcplx)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
                      
          temp_srfcplx_count = temp_srfcplx_count + 1
          srfcplx => SurfaceComplexCreate()
          srfcplx%id = temp_srfcplx_count
          call InputReadWord(input,option,srfcplx%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,SURFACE_COMPLEXATION_RXN,COMPLEX_NAME')
                
          if (.not.associated(srfcplx_rxn%complex_list)) then
            srfcplx_rxn%complex_list => srfcplx
          endif
          if (associated(prev_srfcplx)) then
            prev_srfcplx%next => srfcplx
          endif
          prev_srfcplx => srfcplx
          nullify(srfcplx)
        enddo
        nullify(prev_srfcplx)
      case default
        call InputKeywordUnrecognized(word, &
                'CHEMISTRY,SURFACE_COMPLEXATION_RXN',option)
    end select

  enddo
  
  if (num_times_surface_type_set > 1) then
    option%io_buffer = 'Surface site type (e.g. MINERAL, ROCK_DENSITY, ' // &
      'COLLOID) may only be set once under the SURFACE_COMPLEXATION_RXN card.'
    call printErrMsg(option)
  endif
  
  if (.not.associated(surface_complexation%rxn_list)) then
    surface_complexation%rxn_list => srfcplx_rxn
    srfcplx_rxn%id = 1
  else
    cur_srfcplx_rxn => surface_complexation%rxn_list
    do
      if (.not.associated(cur_srfcplx_rxn%next)) then
        cur_srfcplx_rxn%next => srfcplx_rxn
        srfcplx_rxn%id = cur_srfcplx_rxn%id + 1
        exit
      endif
      cur_srfcplx_rxn => cur_srfcplx_rxn%next
    enddo
    nullify(cur_srfcplx_rxn)
  endif
  
  ! Add surface complexes in reaction to master list, without duplicating.
  ! Set the id of the surface complex to the one in the master list, even
  ! if duplicated.
  cur_srfcplx_in_rxn => srfcplx_rxn%complex_list
  do
    if (.not.associated(cur_srfcplx_in_rxn)) exit
    cur_srfcplx => surface_complexation%complex_list
    found = PETSC_FALSE
    do
      if (.not.associated(cur_srfcplx)) exit
      if (StringCompare(cur_srfcplx_in_rxn%name,cur_srfcplx%name, &
                        MAXWORDLENGTH)) then
        ! set id to id in master list
        cur_srfcplx_in_rxn%id = cur_srfcplx%id
        cur_srfcplx_in_rxn%ptr => cur_srfcplx
        found = PETSC_TRUE
        exit
      endif
      prev_srfcplx => cur_srfcplx
      cur_srfcplx => cur_srfcplx%next
    enddo
    if (.not.found) then
      srfcplx => SurfaceComplexCreate()
      srfcplx%name = cur_srfcplx_in_rxn%name
      if (.not.associated(prev_srfcplx)) then
        surface_complexation%complex_list => srfcplx
        srfcplx%id = 1
      else
        prev_srfcplx%next => srfcplx
        srfcplx%id = prev_srfcplx%id + 1
      endif
    endif
    cur_srfcplx_in_rxn => cur_srfcplx_in_rxn%next
  enddo

  surface_complexation%nsrfcplxrxn = &
    surface_complexation%nsrfcplxrxn + 1
  select case(srfcplx_rxn%itype)
    ! default (NULL) to EQUILIBRIUM
    case(SRFCMPLX_RXN_NULL,SRFCMPLX_RXN_EQUILIBRIUM)
      surface_complexation%neqsrfcplx = surface_complexation%neqsrfcplx + &
        temp_srfcplx_count
      surface_complexation%neqsrfcplxrxn = &
        surface_complexation%neqsrfcplxrxn + 1
    case(SRFCMPLX_RXN_MULTIRATE_KINETIC)
      surface_complexation%nkinmrsrfcplx = &
        surface_complexation%nkinmrsrfcplx + temp_srfcplx_count
      surface_complexation%nkinmrsrfcplxrxn = &
        surface_complexation%nkinmrsrfcplxrxn + 1
      ! if site fractions is not specified, we can calculate this
      ! based on a uniform distribution and the number of rates
      if (.not.associated(srfcplx_rxn%site_fractions) .and. &
          associated(srfcplx_rxn%rates)) then
        allocate(srfcplx_rxn%site_fractions(size(srfcplx_rxn%rates)))
        ! it is possible to specify the number of site fractions
        srfcplx_rxn%site_fractions = 1.d0 / dble(size(srfcplx_rxn%rates))
      endif
      ! check to ensure that rates for multirate surface complexation 
      ! are aligned with surface fractions
      if (size(srfcplx_rxn%rates) /= size(srfcplx_rxn%site_fractions)) then
        write(word,*) size(srfcplx_rxn%rates)
        write(string,*) size(srfcplx_rxn%site_fractions)
        option%io_buffer = 'Number of kinetic rates (' // &
          trim(adjustl(word)) // &
          ') does not match the number of surface fractions (' // &
          trim(adjustl(string)) // ').'
        call printErrMsg(option)
      endif
      tempreal = 0.d0
      do i = 1, size(srfcplx_rxn%site_fractions)
        tempreal = tempreal + srfcplx_rxn%site_fractions(i)
        srfcplx_rxn%rates(i) = srfcplx_rxn%rates(i) * &
          srfcplx_rxn%kinmr_scale_factor
      enddo
    
      if (dabs(1.d0 - tempreal) > 1.d-6) then
        write(string,*) tempreal
        option%io_buffer = 'The sum of the surface site fractions for ' // &
          'multirate kinetic sorption does not add up to 1.d0 (' // &
          trim(adjustl(string)) // '.'
        call printErrMsg(option)
      endif
    case(SRFCMPLX_RXN_KINETIC)
      ! match up rates with their corresponding surface complex
      cur_srfcplx => srfcplx_rxn%complex_list
      do
        if (.not.associated(cur_srfcplx)) exit
        found = PETSC_FALSE
        nullify(prev_srfcplx_rate)
        cur_srfcplx_rate => rate_list
        do
          if (.not.associated(cur_srfcplx_rate)) exit
          ! check for same name
          if (StringCompare(cur_srfcplx_rate%name, &
                            cur_srfcplx%name, &
                            MAXWORDLENGTH)) then
            ! set rates
            cur_srfcplx%forward_rate = cur_srfcplx_rate%forward_rate
            cur_srfcplx%backward_rate = cur_srfcplx_rate%backward_rate
            ! remove srfcplx_rate from list of rates
            if (associated(prev_srfcplx_rate)) then
              prev_srfcplx_rate%next => cur_srfcplx_rate%next
            else
              rate_list => cur_srfcplx_rate%next
            endif
            ! destroy the object
            call SurfaceComplexDestroy(cur_srfcplx_rate)
            found = PETSC_TRUE
            exit
          endif
          prev_srfcplx_rate => cur_srfcplx_rate
          cur_srfcplx_rate => cur_srfcplx_rate%next
        enddo
        if (.not.found) then
          option%io_buffer = 'Rates for surface complex ' // &
            trim(cur_srfcplx%name) // ' not found in kinetic rate list'
          call printErrMsg(option)
        endif
        cur_srfcplx => cur_srfcplx%next
      enddo
      ! check to ensure that rates are matched
      if (associated(rate_list)) then
        option%io_buffer = '# of rates is greater than # of surface complexes'
        call printErrMsg(option)
      endif
      nullify(cur_srfcplx)
      nullify(prev_srfcplx)
      nullify(rate_list)
      nullify(cur_srfcplx_rate)
      nullify(prev_srfcplx_rate)                  
      surface_complexation%nkinsrfcplx = &
        surface_complexation%nkinsrfcplx + temp_srfcplx_count
      surface_complexation%nkinsrfcplxrxn = &
        surface_complexation%nkinsrfcplxrxn + 1
  end select
  srfcplx_rxn%free_site_id = srfcplx_rxn%id
              
  nullify(srfcplx_rxn)

end subroutine SurfaceComplexationRead

! ************************************************************************** !

subroutine SrfCplxProcessConstraint(surface_complexation,constraint_name, &
                                    constraint,option)
  ! 
  ! Initializes constraints based on surface complex
  ! species in system
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/07/13
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(surface_complexation_type), pointer :: surface_complexation
  character(len=MAXWORDLENGTH) :: constraint_name
  type(srfcplx_constraint_type), pointer :: constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: isrfcplx, jsrfcplx
  
  character(len=MAXWORDLENGTH) :: srfcplx_name(surface_complexation%nkinsrfcplx)
  PetscReal :: constraint_conc(surface_complexation%nkinsrfcplx)
  
  if (.not.associated(constraint)) return

  if (surface_complexation%nkinsrfcplx == 0) then
    option%io_buffer = 'Surface complexation specified in constraint "' // &
      trim(constraint_name) // '" requires that kinetic surface ' // &
      'complexation be defined in the CHEMISTRY section.'
    call printErrMsg(option)
  endif
  
  srfcplx_name = ''
  do isrfcplx = 1, surface_complexation%nkinsrfcplx
    found = PETSC_FALSE
    do jsrfcplx = 1, surface_complexation%nkinsrfcplx
      if (StringCompare(constraint%names(isrfcplx), &
                        surface_complexation%srfcplx_names(&
                          !TODO(geh): fix 0 index
                        surface_complexation%kinsrfcplx_to_name(jsrfcplx,0)), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
                'Surface complex ' // trim(constraint%names(isrfcplx)) // &
                'from CONSTRAINT ' // trim(constraint_name) // &
                ' not found among kinetic surface complexes.'
      call printErrMsg(option)
    else
      constraint_conc(jsrfcplx) = &
        constraint%constraint_conc(isrfcplx)
      srfcplx_name(jsrfcplx) = constraint%names(isrfcplx)
    endif  
  enddo
  constraint%names = srfcplx_name
  constraint%constraint_conc = constraint_conc
  
end subroutine SrfCplxProcessConstraint

! ************************************************************************** !

subroutine RTotalSorbEqSurfCplx(rt_auxvar,global_auxvar,material_auxvar, &
                                reaction,option)
  ! 
  ! Computes the total sorbed component concentrations and
  ! derivative with respect to free-ion for equilibrium
  ! surface complexation
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/08; 05/26/09
  ! 

  use Option_module
  use Matrix_Block_Aux_module
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: irxn, ieqrxn
  PetscReal, pointer :: colloid_array_ptr(:)
  PetscInt, parameter :: iphase = 1
  type(matrix_block_auxvar_type), pointer :: colloid_matrix_block_ptr
  type(surface_complexation_type), pointer :: surface_complexation

  surface_complexation => reaction%surface_complexation
  
  if (reaction%ncollcomp > 0) then  
    rt_auxvar%colloid%total_eq_mob = 0.d0
    rt_auxvar%colloid%dRj_dCj%dtotal = 0.d0
    colloid_array_ptr => rt_auxvar%colloid%total_eq_mob
    colloid_matrix_block_ptr => rt_auxvar%colloid%dRj_dCj
  else
    nullify(colloid_matrix_block_ptr)
    nullify(colloid_array_ptr)
  endif

  ! Surface Complexation
  do ieqrxn = 1, surface_complexation%neqsrfcplxrxn
  
    irxn = surface_complexation%eqsrfcplxrxn_to_srfcplxrxn(ieqrxn)
    
    !TODO(geh): clean up colloidpointers
    call RTotalSorbEqSurfCplx1(rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option, &
                               irxn, &
                               rt_auxvar%srfcplxrxn_free_site_conc(irxn), &
                               rt_auxvar%eqsrfcplx_conc, &
                               rt_auxvar%total_sorb_eq, &
                               rt_auxvar%dtotal_sorb_eq, &
                               colloid_array_ptr, &
                               colloid_matrix_block_ptr)    
  
  enddo ! irxn
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
end subroutine RTotalSorbEqSurfCplx

! ************************************************************************** !

subroutine RTotalSorbMultiRateAsEQ(rt_auxvar,global_auxvar,material_auxvar, &
                                   reaction,option)
  ! 
  ! Calculates the multirate surface complexation
  ! reaction as if it were equilibrium.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/19/12
  ! 

  use Option_module
  use Matrix_Block_Aux_module
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: irxn, ikinmrrxn
  PetscInt, parameter :: iphase = 1
  type(surface_complexation_type), pointer :: surface_complexation
  PetscReal :: total_sorb_eq(reaction%naqcomp)
  PetscReal :: dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp)
  PetscReal, pointer :: null_array_ptr(:)
  type(matrix_block_auxvar_type), pointer :: null_matrix_block

  surface_complexation => reaction%surface_complexation

  nullify(null_array_ptr)
  nullify(null_matrix_block)

  ! Surface Complexation
  do ikinmrrxn = 1, surface_complexation%nkinmrsrfcplxrxn
  
    irxn = surface_complexation%kinmrsrfcplxrxn_to_srfcplxrxn(ikinmrrxn)

    total_sorb_eq = 0.d0
    dtotal_sorb_eq = 0.d0  

    call RTotalSorbEqSurfCplx1(rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option, &
                               irxn, &
                               rt_auxvar%srfcplxrxn_free_site_conc(irxn), &
                               null_array_ptr, &
                               total_sorb_eq, &
                               dtotal_sorb_eq, &
                               null_array_ptr, &
                               null_matrix_block)    

    rt_auxvar%kinmr_total_sorb(:,0,ikinmrrxn) = total_sorb_eq(:)
  
  enddo ! irxn
  
end subroutine RTotalSorbMultiRateAsEQ

! ************************************************************************** !

subroutine RMultiRateSorption(Res,Jac,compute_derivative,rt_auxvar, &
                              global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Computes contribution to the accumualtion term due
  ! due to multirate sorption
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/20/09; 03/16/12
  ! 

  use Option_module
  use Matrix_Block_Aux_module
  
  implicit none

  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(option_type) :: option
  
  PetscInt :: irxn, ikinmrrxn
  type(surface_complexation_type), pointer :: surface_complexation
  
  PetscInt :: irate
  PetscInt, parameter :: iphase = 1
  PetscReal :: kdt, one_plus_kdt, k_over_one_plus_kdt
  PetscReal :: total_sorb_eq(reaction%naqcomp)
  PetscReal :: dtotal_sorb_eq(reaction%naqcomp,reaction%naqcomp)
  PetscReal, pointer :: null_array_ptr(:)
  type(matrix_block_auxvar_type), pointer :: null_matrix_block

  surface_complexation => reaction%surface_complexation

  nullify(null_array_ptr)
  nullify(null_matrix_block)

  ! only zero out the zero index.  The other indices hold values from 
  ! the previous time step
  rt_auxvar%kinmr_total_sorb(:,0,:) = 0.d0
  
  ! Surface Complexation
  do ikinmrrxn = 1, surface_complexation%nkinmrsrfcplxrxn
  
    irxn = surface_complexation%kinmrsrfcplxrxn_to_srfcplxrxn(ikinmrrxn)

    total_sorb_eq = 0.d0
    dtotal_sorb_eq = 0.d0  

    call RTotalSorbEqSurfCplx1(rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option, &
                               irxn, &
                               rt_auxvar%srfcplxrxn_free_site_conc(irxn), &
                               null_array_ptr, &
                               total_sorb_eq, &
                               dtotal_sorb_eq, &
                               null_array_ptr, &
                               null_matrix_block)

    ! WARNING: this assumes site fraction multiplicative factor 
    do irate = 1, surface_complexation%kinmr_nrate(ikinmrrxn)
      kdt = surface_complexation%kinmr_rate(irate,ikinmrrxn) * &
            option%tran_dt
      one_plus_kdt = 1.d0 + kdt
      k_over_one_plus_kdt = &
        surface_complexation%kinmr_rate(irate,ikinmrrxn)/one_plus_kdt
        
      ! this is the constribution to the accumulation term in the residual
      Res(:) = Res(:) + material_auxvar%volume * k_over_one_plus_kdt * &
        (surface_complexation%kinmr_frac(irate,ikinmrrxn)*total_sorb_eq(:) - &
         rt_auxvar%kinmr_total_sorb(:,irate,ikinmrrxn))
      
      if (compute_derivative) then
        Jac = Jac + material_auxvar%volume * k_over_one_plus_kdt * &
          surface_complexation%kinmr_frac(irate,ikinmrrxn) * dtotal_sorb_eq
      endif

    enddo

    ! store the target equilibrium concentration to update the sorbed 
    ! concentration at the end of the time step.
    rt_auxvar%kinmr_total_sorb(:,0,ikinmrrxn) = total_sorb_eq
  
  enddo ! ikinmrrxn
  
end subroutine RMultiRateSorption

! ************************************************************************** !

subroutine RTotalSorbEqSurfCplx1(rt_auxvar,global_auxvar,material_auxvar, &
                                 reaction,option, &
                                 irxn,external_free_site_conc, &
                                 external_srfcplx_conc, &
                                 external_total_sorb, &
                                 external_dtotal_sorb, &
                                 external_total_mob, &
                                 external_dRj_dCj)
  ! 
  ! Computes the total sorbed component concentrations and
  ! derivative with respect to free-ion for equilibrium
  ! surface complexation
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/08; 05/26/09; 03/16/12
  ! 

  use Option_module
  use Matrix_Block_Aux_module
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  PetscInt :: irxn
  PetscReal :: external_free_site_conc
  PetscReal, pointer :: external_srfcplx_conc(:)
  PetscReal :: external_total_sorb(reaction%naqcomp)
  PetscReal :: external_dtotal_sorb(reaction%naqcomp,reaction%naqcomp)
  PetscReal, pointer :: external_total_mob(:)
  type(matrix_block_auxvar_type), pointer :: external_dRj_dCj
  
  PetscInt :: i, j, k, icplx, icomp, jcomp, ncomp, ncplx
  PetscReal :: srfcplx_conc(reaction%surface_complexation%nsrfcplx)
  PetscReal :: dSx_dmi(reaction%naqcomp)
  PetscReal :: nui_Si_over_Sx
  PetscReal :: ln_free_site
  PetscReal :: lnQK, tempreal, tempreal1, tempreal2, total
  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: tol = 1.d-12
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscBool :: one_more
  PetscReal :: res, dres_dfree_site, dfree_site_conc
  PetscReal :: free_site_conc, rel_change_in_free_site_conc
  PetscReal :: site_density(2)
  PetscReal :: mobile_fraction
  PetscInt :: num_types_of_sites
  PetscInt :: isite
  PetscInt :: num_iterations
  PetscReal :: damping_factor
  type(surface_complexation_type), pointer :: surface_complexation

  surface_complexation => reaction%surface_complexation
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  
  ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
    
  free_site_conc = external_free_site_conc
  srfcplx_conc = 0.d0

  select case(surface_complexation%srfcplxrxn_surf_type(irxn))
    case(MINERAL_SURFACE)
      site_density(1) = surface_complexation%srfcplxrxn_site_density(irxn)* &
                rt_auxvar%mnrl_volfrac(surface_complexation% &
                                         srfcplxrxn_to_surf(irxn))
      num_types_of_sites = 1
    case(ROCK_SURFACE)
      site_density(1) = surface_complexation%srfcplxrxn_site_density(irxn)* &
                        material_auxvar%soil_particle_density * &
                        (1.d0-material_auxvar%porosity)
      num_types_of_sites = 1
    case(COLLOID_SURFACE)
      mobile_fraction = reaction%colloid_mobile_fraction( &
                          surface_complexation%srfcplxrxn_to_surf(irxn))
      site_density(1) = (1.d0-mobile_fraction)* &
                        surface_complexation%srfcplxrxn_site_density(irxn)
      site_density(2) = mobile_fraction* &
                        surface_complexation%srfcplxrxn_site_density(irxn)
      num_types_of_sites = 2 ! two types of sites (mobile and immobile) with
                             ! separate site densities
    case(NULL_SURFACE)
      site_density(1) = surface_complexation%srfcplxrxn_site_density(irxn)
      num_types_of_sites = 1
  end select
    
  do isite=1, num_types_of_sites
    ! isite == 1 - immobile (colloids, minerals, etc.)
    ! isite == 2 - mobile (colloids)
    
    if (site_density(isite) < 1.d-40) cycle
    
    ! get a pointer to the first complex (there will always be at least 1)
    ! in order to grab free site conc
    one_more = PETSC_FALSE
    num_iterations = 0
    damping_factor = 1.d0
    do

      num_iterations = num_iterations + 1
      total = free_site_conc
      ln_free_site = log(free_site_conc)
      do j = 1, ncplx
        icplx = surface_complexation%srfcplxrxn_to_complex(j,irxn)
        ! compute secondary species concentration
        lnQK = -surface_complexation%srfcplx_logK(icplx)*LOG_TO_LN

        ! activity of water
        if (surface_complexation%srfcplxh2oid(icplx) > 0) then
          lnQK = lnQK + surface_complexation%srfcplxh2ostoich(icplx)* &
                        rt_auxvar%ln_act_h2o
        endif

        lnQK = lnQK + surface_complexation%srfcplx_free_site_stoich(icplx)* &
                      ln_free_site
        
        ncomp = surface_complexation%srfcplxspecid(0,icplx)
        do i = 1, ncomp
          icomp = surface_complexation%srfcplxspecid(i,icplx)
          lnQK = lnQK + surface_complexation%srfcplxstoich(i,icplx)* &
                        ln_act(icomp)
        enddo
        srfcplx_conc(icplx) = exp(lnQK)
        total = total + surface_complexation%srfcplx_free_site_stoich(icplx)* &
                        srfcplx_conc(icplx) 
          
      enddo
        
      if (one_more) exit
        
      if (surface_complexation%srfcplxrxn_stoich_flag(irxn)) then 
        ! stoichiometry for free sites in one of reactions is not 1, thus must
        ! use nonlinear iteration to solve
        res = site_density(isite)-total
          
        dres_dfree_site = 1.d0

        do j = 1, ncplx
          icplx = surface_complexation%srfcplxrxn_to_complex(j,irxn)
          dres_dfree_site = dres_dfree_site + &
            surface_complexation%srfcplx_free_site_stoich(icplx)* &
            srfcplx_conc(icplx)/free_site_conc
        enddo

        dfree_site_conc = res / dres_dfree_site
        
        ! if number of Newton iterations is excessive, try damping the update
        if (num_iterations > 1000) then
          damping_factor = 0.5d0
        endif
        free_site_conc = free_site_conc + damping_factor*dfree_site_conc
        rel_change_in_free_site_conc = dabs(dfree_site_conc/free_site_conc)
        
        if (rel_change_in_free_site_conc < tol) then
          one_more = PETSC_TRUE
        endif
        
      else
        
        total = total / free_site_conc
        free_site_conc = site_density(isite) / total  
          
        one_more = PETSC_TRUE 
        
      endif

    enddo ! generic do
      
    external_free_site_conc = free_site_conc
   
!!!!!!!!!!!!
    ! 2.3-46

    ! Sx = free site
    ! mi = molality of component i
    dSx_dmi = 0.d0
    tempreal = 0.d0
    do j = 1, ncplx
      icplx = surface_complexation%srfcplxrxn_to_complex(j,irxn)
      ncomp = surface_complexation%srfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = surface_complexation%srfcplxspecid(i,icplx)
        ! sum of nu_li * nu_i * S_i
        dSx_dmi(icomp) = dSx_dmi(icomp) + &
          surface_complexation%srfcplxstoich(i,icplx)* &
          surface_complexation%srfcplx_free_site_stoich(icplx)* &
          srfcplx_conc(icplx)
      enddo
      ! sum of nu_i^2 * S_i
      tempreal = tempreal + &
        surface_complexation%srfcplx_free_site_stoich(icplx)* & 
        surface_complexation%srfcplx_free_site_stoich(icplx)* &
        srfcplx_conc(icplx)
    enddo 
    ! divide denominator by Sx
    tempreal = tempreal / free_site_conc
    ! add 1.d0 to denominator
    tempreal = tempreal + 1.d0
    ! divide numerator by denominator
    dSx_dmi = -dSx_dmi / tempreal
    ! convert from dlogm to dm
    dSx_dmi = dSx_dmi / rt_auxvar%pri_molal
!!!!!!!!!!!!
    
    if (isite == 1 .and. associated(external_srfcplx_conc)) then
      external_srfcplx_conc(:) = external_srfcplx_conc(:) + srfcplx_conc(:)
    endif
#if 0
!geh: for use if we decide to map surface complexes to master
    if (isite == 1 .and. associated(ext_srfcplx_to_rxn_srfcplx_map)) then
      do i = 1, surface_complexation%nsrfcplx
        icplx = ext_srfcplx_to_rxn_srfcplx_map(i)
        if (icplx > 0) then
          external_srfcplx_conc(icplx) = external_srfcplx_conc(icplx) + &
                                         srfcplx_conc(i)
        endif
      enddo
    endif    
#endif
   
    do k = 1, ncplx
      icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)

      ncomp = surface_complexation%srfcplxspecid(0,icplx)
      if (isite == 1) then ! immobile sites  
        do i = 1, ncomp
          icomp = surface_complexation%srfcplxspecid(i,icplx)
          external_total_sorb(icomp) = external_total_sorb(icomp) + &
            surface_complexation%srfcplxstoich(i,icplx)*srfcplx_conc(icplx)
        enddo
      else ! mobile sites
        do i = 1, ncomp
          icomp = reaction%pri_spec_to_coll_spec(surface_complexation% &
                                                 srfcplxspecid(i,icplx))
          external_total_mob(icomp) = external_total_mob(icomp) + &
            surface_complexation%srfcplxstoich(i,icplx)*srfcplx_conc(icplx)
        enddo
      endif
        
      ! for 2.3-47 which feeds into 2.3-50
      nui_Si_over_Sx = surface_complexation%srfcplx_free_site_stoich(icplx)* &
                        srfcplx_conc(icplx)/ &
                        free_site_conc

      do j = 1, ncomp
        jcomp = surface_complexation%srfcplxspecid(j,icplx)
        tempreal = surface_complexation%srfcplxstoich(j,icplx)* &
                   srfcplx_conc(icplx) / &
                   rt_auxvar%pri_molal(jcomp)+ &
                   nui_Si_over_Sx*dSx_dmi(jcomp)
        if (isite == 1) then ! immobile sites                  
          do i = 1, ncomp
            icomp = surface_complexation%srfcplxspecid(i,icplx)
            external_dtotal_sorb(icomp,jcomp) = &
                                external_dtotal_sorb(icomp,jcomp) + &
                                surface_complexation%srfcplxstoich(i,icplx)* &
                                tempreal
          enddo ! i
        else ! mobile sites
          do i = 1, ncomp
            icomp = surface_complexation%srfcplxspecid(i,icplx)
            external_dRj_dCj%dtotal(icomp,jcomp,1) = &
                                 external_dRj_dCj%dtotal(icomp,jcomp,1) + &
                                 surface_complexation%srfcplxstoich(i,icplx)* &
                                 tempreal
          enddo ! i
        endif
      enddo ! j
    enddo ! k
  enddo ! isite
  
end subroutine RTotalSorbEqSurfCplx1

! ************************************************************************** !

subroutine RKineticSurfCplx(Res,Jac,compute_derivative,rt_auxvar, &
                            global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Computes contribution to residual and jacobian for
  ! kinetic surface complexation reactions
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/07/09
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none
  
  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(option_type) :: option

  PetscInt :: i, j, k, l, icplx, icomp, jcomp, lcomp, ncomp, ncplx
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscInt :: irxn, isite, ikinrxn
  PetscReal :: dt
  type(surface_complexation_type), pointer :: surface_complexation

  PetscReal :: numerator_sum(reaction%surface_complexation%nkinsrfcplxrxn)
  PetscReal :: denominator_sum(reaction%surface_complexation%nkinsrfcplxrxn)

  PetscReal :: denominator
  PetscReal :: fac
  PetscReal :: fac_sum(reaction%naqcomp)
  PetscReal :: lnQ(reaction%surface_complexation%nkinsrfcplx)
  PetscReal :: Q(reaction%surface_complexation%nkinsrfcplx)
  PetscReal :: srfcplx_conc_k(maxval(reaction%surface_complexation%srfcplxrxn_to_complex(0,:)),1)
  PetscReal :: srfcplx_conc_kp1(maxval(reaction%surface_complexation%srfcplxrxn_to_complex(0,:)),1)
  
  surface_complexation => reaction%surface_complexation
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

! Members of the rt aux var object: mol/m^3
! PetscReal, pointer :: kinsrfcplx_conc(:)          ! S_{i\alpha}^k
! PetscReal, pointer :: kinsrfcplx_conc_kp1(:)      ! S_{i\alpha}^k+1
! PetscReal, pointer :: kinsrfcplx_free_site_conc(:) ! S_\alpha
  
! units
! k_f: dm^3/mol/sec
! k_b: 1/sec
! Res: mol/sec

  dt = option%tran_dt
  
! compute ion activity product and store: units mol/L
  lnQ = 0.d0
  do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
    irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
    ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)
      if (surface_complexation%srfcplxh2oid(icplx) > 0) then
        lnQ(icplx) = lnQ(icplx) + surface_complexation%srfcplxh2ostoich(icplx)* &
          rt_auxvar%ln_act_h2o
      endif
    
      ncomp = surface_complexation%srfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = surface_complexation%srfcplxspecid(i,icplx)
        lnQ(icplx) = lnQ(icplx) + surface_complexation%srfcplxstoich(i,icplx)* &
          ln_act(icomp)
      enddo
      Q(icplx) = exp(lnQ(icplx))
    enddo
  enddo
    
  ! compute summation in numerator of 5.1-29: units mol/m^3
  numerator_sum = 0.d0
  do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
    irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
    isite = surface_complexation%srfcplxrxn_to_surf(irxn)
    ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)
      numerator_sum(isite) = numerator_sum(isite) + &
                  rt_auxvar%kinsrfcplx_conc(icplx,ikinrxn)/ &
                  (1.d0+surface_complexation%kinsrfcplx_backward_rate(icplx,ikinrxn)*dt)
    enddo
  enddo

  do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
    irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
    isite = surface_complexation%srfcplxrxn_to_surf(irxn)
    numerator_sum(isite) = surface_complexation%srfcplxrxn_site_density(isite) - &
                           numerator_sum(isite)
  enddo
  
  ! compute summation in denominator of 5.1-29
  denominator_sum = 1.d0
  do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
    irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
    isite = surface_complexation%srfcplxrxn_to_surf(irxn)
    ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)
      denominator_sum(isite) = denominator_sum(isite) + &
                               (surface_complexation%kinsrfcplx_forward_rate(icplx,ikinrxn)*dt)/ &
                               (1.d0+surface_complexation%kinsrfcplx_backward_rate(icplx,ikinrxn)*dt)* &
                               Q(icplx)
    enddo
  enddo

! compute surface complex conc. at new time step (5.1-30)  
  do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
    irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
    isite = surface_complexation%srfcplxrxn_to_surf(irxn)
    ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)
      srfcplx_conc_k(icplx,ikinrxn) = rt_auxvar%kinsrfcplx_conc(icplx,ikinrxn)
      denominator = 1.d0 + surface_complexation%kinsrfcplx_backward_rate(icplx,ikinrxn)*dt
      srfcplx_conc_kp1(icplx,ikinrxn) = (srfcplx_conc_k(icplx,ikinrxn) + &
                                surface_complexation%kinsrfcplx_forward_rate(icplx,ikinrxn)*dt * &
                                numerator_sum(isite)/denominator_sum(isite)* &
                                Q(icplx))/denominator
      rt_auxvar%kinsrfcplx_conc_kp1(icplx,ikinrxn) = srfcplx_conc_kp1(icplx,ikinrxn)
    enddo
    rt_auxvar%kinsrfcplx_free_site_conc(isite) = numerator_sum(isite)/ &
                                               denominator_sum(isite)
  enddo

! compute residual (5.1-34)
  
  do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
    irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
    ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
    do k = 1, ncplx ! ncplx in rxn
      icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)
      ncomp = surface_complexation%srfcplxspecid(0,icplx)
      do i = 1, ncomp
        icomp = surface_complexation%srfcplxspecid(i,icplx)
        Res(icomp) = Res(icomp) + surface_complexation%srfcplxstoich(i,icplx)* &
                     (srfcplx_conc_kp1(icplx,ikinrxn)-rt_auxvar%kinsrfcplx_conc(icplx,ikinrxn))/ &
                     dt * material_auxvar%volume
      enddo
    enddo
  enddo

  if (compute_derivative) then
  ! compute jacobian (5.1-39)
    fac_sum = 0.d0
    do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
      irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
      ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
      do k = 1, ncplx ! ncplx in rxn
        icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)
        denominator = 1.d0 + surface_complexation%kinsrfcplx_backward_rate(icplx,ikinrxn)*dt
        fac = surface_complexation%kinsrfcplx_forward_rate(icplx,ikinrxn)/denominator
        ncomp = surface_complexation%srfcplxspecid(0,icplx)
        do j = 1, ncomp
          jcomp = surface_complexation%srfcplxspecid(j,icplx)
          fac_sum(jcomp) = fac_sum(jcomp) + surface_complexation%srfcplxstoich(j,icplx)* &
            fac * Q(icplx)
        enddo
      enddo
    enddo

    do ikinrxn = 1, surface_complexation%nkinsrfcplxrxn
      irxn = surface_complexation%kinsrfcplxrxn_to_srfcplxrxn(ikinrxn)
      isite = surface_complexation%srfcplxrxn_to_surf(irxn)
      ncplx = surface_complexation%srfcplxrxn_to_complex(0,irxn)
      do k = 1, ncplx ! ncplx in rxn
        icplx = surface_complexation%srfcplxrxn_to_complex(k,irxn)
        denominator = 1.d0 + surface_complexation%kinsrfcplx_backward_rate(icplx,ikinrxn)*dt
        fac = surface_complexation%kinsrfcplx_forward_rate(icplx,ikinrxn)/denominator
        ncomp = surface_complexation%srfcplxspecid(0,icplx)
        do j = 1, ncomp
          jcomp = surface_complexation%srfcplxspecid(j,icplx)
          do l = 1, ncomp
            lcomp = surface_complexation%srfcplxspecid(l,icplx)
            Jac(jcomp,lcomp) = Jac(jcomp,lcomp) + &
              (surface_complexation%srfcplxstoich(j,icplx) * fac * numerator_sum(isite) * &
              Q(icplx) * (surface_complexation%srfcplxstoich(l,icplx) - &
              dt * fac_sum(lcomp)/denominator_sum(isite)))/denominator_sum(isite) * &
              exp(-ln_conc(lcomp)) * material_auxvar%volume
          enddo
        enddo
      enddo
    enddo
  endif
  
  ! units of total_sorb = mol/m^3
  ! units of dtotal_sorb = kg water/m^3 bulk
  
end subroutine RKineticSurfCplx

end module Reaction_Surface_Complexation_module
