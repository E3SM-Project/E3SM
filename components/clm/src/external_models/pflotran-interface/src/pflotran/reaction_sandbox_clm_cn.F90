module Reaction_Sandbox_CLM_CN_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter :: CN_ratio_mass_to_mol = 1.16616d0 
  PetscInt, parameter :: CARBON_INDEX = 1
  PetscInt, parameter :: NITROGEN_INDEX = 2
  PetscInt, parameter :: SOM_INDEX = 1

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_clm_cn_type
    PetscInt :: nrxn
    PetscInt :: npool
    PetscReal, pointer :: CN_ratio(:)
    PetscReal, pointer :: rate_constant(:)
    PetscReal, pointer :: respiration_fraction(:)
    PetscReal, pointer :: inhibition_constant(:)
    PetscInt, pointer :: upstream_pool_id(:)
    PetscInt, pointer :: downstream_pool_id(:)
    PetscInt, pointer :: pool_id_to_species_id(:,:)
    PetscInt :: C_species_id
    PetscInt :: N_species_id
    type(pool_type), pointer :: pools
    type(clm_cn_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => CLM_CN_Read
    procedure, public :: Setup => CLM_CN_Setup
    procedure, public :: Evaluate => CLM_CN_React
    procedure, public :: Destroy => CLM_CN_Destroy
  end type reaction_sandbox_clm_cn_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: CN_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: clm_cn_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    character(len=MAXWORDLENGTH) :: downstream_pool_name
    PetscReal :: rate_constant
    PetscReal :: respiration_fraction
    PetscReal :: inhibition_constant
    type(clm_cn_reaction_type), pointer :: next
  end type clm_cn_reaction_type
  
  public :: CLM_CN_Create

contains

! ************************************************************************** !

function CLM_CN_Create()
  ! 
  ! Allocates CLM-CN reaction object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  implicit none
  
  type(reaction_sandbox_clm_cn_type), pointer :: CLM_CN_Create
  
  allocate(CLM_CN_Create)
  CLM_CN_Create%nrxn = 0
  CLM_CN_Create%npool = 0
  nullify(CLM_CN_Create%CN_ratio)
  nullify(CLM_CN_Create%rate_constant)
  nullify(CLM_CN_Create%respiration_fraction)
  nullify(CLM_CN_Create%inhibition_constant)
  nullify(CLM_CN_Create%pool_id_to_species_id)
  nullify(CLM_CN_Create%upstream_pool_id)
  nullify(CLM_CN_Create%downstream_pool_id)
  CLM_CN_Create%C_species_id = 0
  CLM_CN_Create%N_species_id = 0
  nullify(CLM_CN_Create%next)
  nullify(CLM_CN_Create%pools)
  nullify(CLM_CN_Create%reactions)

end function CLM_CN_Create

! ************************************************************************** !

subroutine CLM_CN_Read(this,input,option)
  ! 
  ! Reads input deck for reaction sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(clm_cn_reaction_type), pointer :: new_reaction, prev_reaction
  
  PetscReal :: rate_constant, turnover_time
  
  nullify(new_pool)
  nullify(prev_pool)
  
  nullify(new_reaction)
  nullify(prev_reaction)
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CLM-CN')
    call StringToUpper(word)   

    select case(trim(word))
      case('POOLS')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit   

          allocate(new_pool)
          new_pool%name = ''
          new_pool%CN_ratio = 0.d0
          nullify(new_pool%next)

          call InputReadWord(input,option,new_pool%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'pool name', &
            'CHEMISTRY,REACTION_SANDBOX,CLM-CN,POOLS')
          call InputReadDouble(input,option,new_pool%CN_ratio)
          if (InputError(input)) then
            new_pool%CN_ratio = UNINITIALIZED_DOUBLE
          else
            ! convert CN ratio from mass C/mass N to mol C/mol N
            new_pool%CN_ratio = new_pool%CN_ratio * CN_ratio_mass_to_mol
          endif
          if (associated(this%pools)) then
            prev_pool%next => new_pool
          else
            this%pools => new_pool
          endif
          prev_pool => new_pool
          nullify(new_pool)
        enddo
      case('REACTION')
      
        allocate(new_reaction)
        new_reaction%upstream_pool_name = ''
        new_reaction%downstream_pool_name = ''
        new_reaction%rate_constant = UNINITIALIZED_DOUBLE
        new_reaction%respiration_fraction = UNINITIALIZED_DOUBLE
        new_reaction%inhibition_constant = 0.d0
        nullify(new_reaction%next)
        
        ! need to set these temporarily in order to check that they
        ! are not both set.
        turnover_time = 0.d0
        rate_constant = 0.d0
        
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION')
            case('DOWNSTREAM_POOL')
              call InputReadWord(input,option, &
                                 new_reaction%downstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION')
            case('RATE_CONSTANT')
              call InputReadDouble(input,option,rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              internal_units = 'unitless/sec'
              if (InputError(input)) then
                input%err_buf = 'CLM-CN RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,internal_units,option)
              endif
            case('TURNOVER_TIME')
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              internal_units = 'sec'
              if (InputError(input)) then
                input%err_buf = 'CLM-CN TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,internal_units,option)
              endif
            case('RESPIRATION_FRACTION')
              call InputReadDouble(input,option, &
                                   new_reaction%respiration_fraction)
              call InputErrorMsg(input,option,'respiration fraction', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION')
            case('N_INHIBITION')
              call InputReadDouble(input,option, &
                                   new_reaction%inhibition_constant)
              call InputErrorMsg(input,option,'inhibition constant', &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION')
            case default
              call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN,REACTION',option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
            'be included in a CLM-CN reaction definition, but not both. ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        else if (turnover_time > 0.d0) then
          new_reaction%rate_constant = 1.d0 / turnover_time
        else
          new_reaction%rate_constant = rate_constant
        endif
        ! ensure that respiration fraction is 0-1.
        if (new_reaction%respiration_fraction < 0.d0 .or. &
                  new_reaction%respiration_fraction > 1.d0) then
          option%io_buffer = 'Respiratory fraction (rf) must be between ' // &
            'zero and one (i.e. 0. <= rf <= 1.) in a CLM-CN reaction ' // &
            'definition. See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        endif
        ! If no downstream pool exists, ensure that respiration fraction = 1
        if (len_trim(new_reaction%downstream_pool_name) < 1 .and. &
            (1.d0 - new_reaction%respiration_fraction) > 1.d-40) then
          option%io_buffer = 'Respiratory fraction (rf) must be set to ' // &
            '1.0 if no downstream pool is specified in a CLM-CN reaction ' // &
            'definition. See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        endif
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,CLM-CN',option)
    end select
  enddo
  
end subroutine CLM_CN_Read

! ************************************************************************** !

subroutine CLM_CN_Setup(this,reaction,option)
  ! 
  ! Sets up CLM-CN reaction after it has been read from input
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Reaction_Aux_module, only : reaction_type
  use Option_module
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(reaction_type) :: reaction  
  type(option_type) :: option
  
  call CLM_CN_Map(this,reaction,option)

end subroutine CLM_CN_Setup

! ************************************************************************** !

subroutine CLM_CN_Map(this,reaction,option)
  ! 
  ! Maps coefficients to primary dependent variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Aux_module, only : reaction_type
  use Option_module
  use String_module
  use Reaction_Immobile_Aux_module
  
  implicit none

  class(reaction_sandbox_clm_cn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH), allocatable :: pool_names(:)
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt :: icount

  type(pool_type), pointer :: cur_pool
  type(clm_cn_reaction_type), pointer :: cur_rxn
  
  ! count # pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    cur_pool => cur_pool%next
  enddo
  this%npool = icount
  
  ! count # reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    cur_rxn => cur_rxn%next
  enddo
  this%nrxn = icount
  
  ! allocate and initialize arrays
  allocate(this%CN_ratio(this%npool))
  allocate(this%pool_id_to_species_id(0:2,this%npool))
  allocate(this%upstream_pool_id(this%nrxn))
  allocate(this%downstream_pool_id(this%nrxn))
  allocate(this%rate_constant(this%nrxn))
  allocate(this%respiration_fraction(this%nrxn))
  allocate(this%inhibition_constant(this%nrxn))
  this%CN_ratio = 0.d0
  this%pool_id_to_species_id = 0
  this%upstream_pool_id = 0
  this%downstream_pool_id = 0
  this%rate_constant = 0.d0
  this%respiration_fraction = 0.d0
  this%inhibition_constant = 0.d0
  
  ! temporary array for mapping pools in reactions
  allocate(pool_names(this%npool))
  
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    this%CN_ratio(icount) = cur_pool%CN_ratio
    pool_names(icount) = cur_pool%name
    if (cur_pool%CN_ratio < 0.d0) then
      ! Since no CN ratio provided, must provide two species with the
      ! same name as the pool with C or N appended.
      word = trim(cur_pool%name) // 'C'
      this%pool_id_to_species_id(CARBON_INDEX,icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      word = trim(cur_pool%name) // 'N'
      this%pool_id_to_species_id(NITROGEN_INDEX,icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      this%pool_id_to_species_id(0,icount) = 2
      if (minval(this%pool_id_to_species_id(:,icount)) <= 0) then
        option%io_buffer = 'For CLM-CN pools with no CN ratio defined, ' // &
          'the user must define two immobile species with the same root ' // &
          'name as the pool with "C" or "N" appended, respectively.'
        call printErrMsg(option)
      endif
    else ! only one species (e.g. SOMX)
      this%pool_id_to_species_id(SOM_INDEX,icount) = &
        GetImmobileSpeciesIDFromName(cur_pool%name,reaction%immobile, &
                                     PETSC_TRUE,option)
      this%pool_id_to_species_id(0,icount) = 1
    endif
    cur_pool => cur_pool%next
  enddo
  
  ! map C and N species (solid phase for now)
  word = 'C'
  this%C_species_id = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                   PETSC_TRUE,option)
  word = 'N'
  this%N_species_id = &
      GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                   PETSC_TRUE,option)
  
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    this%upstream_pool_id(icount) = &
      StringFindEntryInList(cur_rxn%upstream_pool_name,pool_names)
    if (this%upstream_pool_id(icount) == 0) then
      option%io_buffer = 'Upstream pool "' // &
        trim(cur_rxn%upstream_pool_name) // &
        '" in reaction with downstream pool "' // &
        trim(cur_rxn%downstream_pool_name) // '" not found in list of pools.'
      call printErrMsg(option)
    endif
    if (len_trim(cur_rxn%downstream_pool_name) > 0) then
      this%downstream_pool_id(icount) = &
        StringFindEntryInList(cur_rxn%downstream_pool_name,pool_names)
      if (this%downstream_pool_id(icount) == 0) then
        option%io_buffer = 'Downstream pool "' // &
          trim(cur_rxn%downstream_pool_name) // &
          '" in reaction with upstream pool "' // &
          trim(cur_rxn%upstream_pool_name) // '" not found in list of pools.'
        call printErrMsg(option)
      endif
      if (this%CN_ratio(this%downstream_pool_id(icount)) < 0.d0) then
        option%io_buffer = 'For CLM-CN reactions, downstream pools ' // &
          'must have a constant C:N ratio (i.e. C and N are not tracked ' // &
          ' individually.  Therefore, pool "' // &
          trim(cur_rxn%downstream_pool_name) // &
          '" may not be used as a downstream pool.'
        call printErrMsg(option)
      endif
    endif
    this%rate_constant(icount) = cur_rxn%rate_constant
    this%respiration_fraction(icount) = cur_rxn%respiration_fraction
    this%inhibition_constant(icount) = cur_rxn%inhibition_constant
    cur_rxn => cur_rxn%next
  enddo 
  
  deallocate(pool_names)
  
end subroutine CLM_CN_Map

! ************************************************************************** !

subroutine CLM_CN_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                        global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Reaction_Aux_module, only : reaction_type
  use Material_Aux_class, only : material_auxvar_type
  
  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscInt :: ipool_up, ipool_down
  PetscInt :: ispec_pool_down
  PetscInt :: ispecC_pool_up, ispecN_pool_up
  PetscInt :: ires_pool_down, ires_C, ires_N
  PetscInt :: iresC_pool_up, iresN_pool_up
  PetscInt :: ispec_N
  PetscReal :: drate, scaled_rate_const, rate
  PetscInt :: irxn
  
  ! inhibition variables
  PetscReal :: F_t
  PetscReal :: F_theta
  PetscReal :: constant_inhibition
  PetscReal :: temp_K
  PetscReal, parameter :: one_over_71_02 = 1.408054069d-2
  PetscReal, parameter :: theta_min = 0.01d0     ! 1/nat log(0.01d0)
  PetscReal, parameter :: one_over_log_theta_min = -2.17147241d-1
  PetscReal, parameter :: twelve_over_14 = 0.857142857143d0

  PetscReal :: CN_ratio_up, CN_ratio_down
  PetscBool :: constant_CN_ratio_up
  PetscReal :: resp_frac
  PetscReal :: stoich_N
  PetscReal :: stoich_C
  PetscReal :: stoich_downstreamC_pool
  PetscReal :: stoich_upstreamC_pool, stoich_upstreamN_pool
  
  PetscReal :: N_inhibition, d_N_inhibition
  PetscReal :: drate_dN_inhibition
  PetscBool :: use_N_inhibition
  PetscReal :: temp_real
  
  PetscReal :: dCN_ratio_up_dC_pool_up, dCN_ratio_up_dN_pool_up
  PetscReal :: dstoich_upstreamN_pool_dC_pool_up
  PetscReal :: dstoich_upstreamN_pool_dN_pool_up
  PetscReal :: dstoichN_dC_pool_up
  PetscReal :: dstoichN_dN_pool_up
  
  PetscReal :: sumC
  PetscReal :: sumN
  
  ! inhibition due to temperature
  ! Equation: F_t = exp(308.56*(1/17.02 - 1/(T - 227.13)))
  temp_K = global_auxvar%temp + 273.15d0

  if (temp_K > 227.15d0) then
    F_t = exp(308.56d0*(one_over_71_02 - 1.d0/(temp_K - 227.13d0)))
  else
    F_t = 0.0
    return
  endif
  
  ! inhibition due to moisture content
  ! Equation: F_theta = log(theta_min/theta) / log(theta_min/theta_max)
  ! Assumptions: theta is saturation
  !              theta_min = 0.01, theta_max = 1.
  F_theta = log(theta_min/max(theta_min,global_auxvar%sat(1))) * one_over_log_theta_min 
  
  constant_inhibition = F_t * F_theta
  
  ! indices for C and N species
  ires_C = reaction%offset_immobile + this%C_species_id
  ispec_N = this%N_species_id
  ires_N = reaction%offset_immobile + ispec_N

  do irxn = 1, this%nrxn
  
    sumC = 0.d0
    sumN = 0.d0
  
    ! scaled_rate_const units: (m^3 bulk / s) = (1/s) * (m^3 bulk)
    scaled_rate_const = this%rate_constant(irxn)*material_auxvar%volume* &
                        constant_inhibition
    resp_frac = this%respiration_fraction(irxn)
    
    ipool_up = this%upstream_pool_id(irxn)
    constant_CN_ratio_up = (this%pool_id_to_species_id(0,ipool_up) == 1)

    if (.not. constant_CN_ratio_up) then
      ! upstream pool is Litter pool with two species (C,N)
      !
      !  a LitC_i + b LitN_i -> c SOM_i + d C + e N
      !
      ispecC_pool_up = this%pool_id_to_species_id(CARBON_INDEX,ipool_up)
      ispecN_pool_up = this%pool_id_to_species_id(NITROGEN_INDEX,ipool_up)
      CN_ratio_up = rt_auxvar%immobile(ispecC_pool_up) / &
                    rt_auxvar%immobile(ispecN_pool_up)
    else
      ! upstream pool is an SOM pool with one species
      !
      !  a SOM_i -> c SOM_(i+1) + d C + e N
      !
      ispecC_pool_up = this%pool_id_to_species_id(SOM_INDEX,ipool_up)
      CN_ratio_up = this%CN_ratio(ipool_up)
    endif

    ! a = fraction_C_up = 1.
    stoich_upstreamC_pool = 1.d0
    ! b = a / CN_ratio_up
    stoich_upstreamN_pool = stoich_upstreamC_pool / CN_ratio_up

    ! downstream pool
    ipool_down = this%downstream_pool_id(irxn)
    ! a downstream pool need not exist if last in succession and
    ! respiration fraction = 1.
    if (ipool_down > 0) then
      ! downstream pool is always SOM (i.e. not split between C and N
      ! as a litter would be).
      ispec_pool_down = this%pool_id_to_species_id(SOM_INDEX,ipool_down)
      CN_ratio_down = this%CN_ratio(ipool_down)
      ! c = (1-resp_frac) * a
      stoich_downstreamC_pool = (1.d0-resp_frac) * stoich_upstreamC_pool
    else    
      ispec_pool_down = 0
      stoich_downstreamC_pool = 0.d0
      CN_ratio_down = 1.d0 ! to prevent divide by zero below.
    endif
      
    ! d = resp_frac * a
    stoich_C = resp_frac * stoich_upstreamC_pool
    ! e = b - c / CN_ratio_dn
    stoich_N = stoich_upstreamN_pool - stoich_downstreamC_pool / CN_ratio_down
 
    ! Inhibition by nitrogen (inhibition concentration > 0 and N is a reactant)
    ! must be calculated here as the sign on the stoichiometry for N is 
    ! required.
    if (this%inhibition_constant(irxn) > 1.d-40 .and. stoich_N < 0.d0) then
      use_N_inhibition = PETSC_TRUE
      temp_real = rt_auxvar%immobile(ispec_N) + &
                  this%inhibition_constant(irxn)
      N_inhibition = rt_auxvar%immobile(ispec_N) / temp_real
      d_N_inhibition = this%inhibition_constant(irxn) / &
                             (temp_real * temp_real)
    else 
      use_N_inhibition = PETSC_FALSE
      N_inhibition = 1.d0
      d_N_inhibition = 0.d0
    endif
    
    ! residual units: (mol/sec) = (m^3 bulk/s) * (mol/m^3 bulk)
    rate = rt_auxvar%immobile(ispecC_pool_up) * &
           scaled_rate_const * N_inhibition

    ! calculation of residual
    
    ! carbon
    Residual(ires_C) = Residual(ires_C) - stoich_C * rate
    sumC = sumC + stoich_C * rate
    
    ! nitrogen
    Residual(ires_N) = Residual(ires_N) - stoich_N * rate
    sumN = sumN + stoich_N * rate

    ! C species in upstream pool (litter or SOM)
    iresC_pool_up = reaction%offset_immobile + ispecC_pool_up
    ! scaled by negative one since it is a reactant 
    Residual(iresC_pool_up) = Residual(iresC_pool_up) - &
      (-1.d0) * stoich_upstreamC_pool * rate
    sumC = sumC - stoich_upstreamC_pool * rate
    
    ! N species in upstream pool (litter only)
    if (.not.constant_CN_ratio_up) then
      ! N species in upstream pool
      iresN_pool_up = reaction%offset_immobile + ispecN_pool_up
      ! scaled by negative one since it is a reactant 
      Residual(iresN_pool_up) = Residual(iresN_pool_up) - &
        (-1.d0) * stoich_upstreamN_pool * rate
    endif
    sumN = sumN - stoich_upstreamN_pool * rate
    
    if (ispec_pool_down > 0) then
      ! downstream pool
      ires_pool_down = reaction%offset_immobile + ispec_pool_down
      Residual(ires_pool_down) = Residual(ires_pool_down) - &
        stoich_downstreamC_pool * rate
      sumC = sumC + stoich_downstreamC_pool * rate
      sumN = sumN + stoich_downstreamC_pool / CN_ratio_down * rate
    endif
    
    !for debugging
!    if (dabs(sumC) > 1.d-40 .or. dabs(sumN) > 1.d-40) then
!      print *, "sum C: ", sumC
!      print *, "sum N: ", sumN
!    endif
    
    if (compute_derivative) then
    
      drate = scaled_rate_const * N_inhibition
      
      ! upstream C pool
      Jacobian(iresC_pool_up,iresC_pool_up) = &
        Jacobian(iresC_pool_up,iresC_pool_up) - &
        ! scaled by negative one since it is a reactant 
        (-1.d0) * stoich_upstreamC_pool * drate
      if (use_N_inhibition) then
!       revision to avoid division by 0 when N -> 0 (N_inhibition -> 0)
        drate_dN_inhibition = rt_auxvar%immobile(ispecC_pool_up) * &
           scaled_rate_const * d_N_inhibition

        ! scaled by negative one since it is a reactant 
        Jacobian(iresC_pool_up,ires_N) = &
          Jacobian(iresC_pool_up,ires_N) - &
          (-1.d0) * stoich_upstreamC_pool * drate_dN_inhibition
      endif
      
      ! downstream pool
      if (ispec_pool_down > 0) then
        Jacobian(ires_pool_down,iresC_pool_up) = &
          Jacobian(ires_pool_down,iresC_pool_up) - &
          stoich_downstreamC_pool * drate
        if (use_N_inhibition) then
          Jacobian(ires_pool_down,ires_N) = &
            Jacobian(ires_pool_down,ires_N) - &
            stoich_downstreamC_pool * drate_dN_inhibition
        endif
      endif
      
      ! variable upstream N pool (for litter pools only!!!)
      if (.not.constant_CN_ratio_up) then

        ! derivative of upstream N pool with respect to upstream C pool
        ! scaled by negative one since it is a reactant 
        Jacobian(iresN_pool_up,iresC_pool_up) = &
          Jacobian(iresN_pool_up,iresC_pool_up) - &
          (-1.d0) * stoich_upstreamN_pool * drate
        if (use_N_inhibition) then
          ! scaled by negative one since it is a reactant 
          Jacobian(iresN_pool_up,ires_N) = &
            Jacobian(iresN_pool_up,ires_N) - &
            (-1.d0) * stoich_upstreamN_pool * drate_dN_inhibition
        endif

        ! Since CN_ratio is a function of unknowns, other derivatives to 
        ! calculate
        dCN_ratio_up_dC_pool_up = 1.d0 / rt_auxvar%immobile(ispecN_pool_up)
        dCN_ratio_up_dN_pool_up = -1.d0 * CN_ratio_up / &
                                  rt_auxvar%immobile(ispecN_pool_up)
        
        ! stoich_upstreamC_pool = 1.d0 (for upstream litter pool)
        ! dstoich_upstreamC_pool_dC_up = 0.
        ! stoich_upstreamN_pool = stoich_upstreamC_pool / CN_ratio_up
        temp_real = -1.d0 * stoich_upstreamN_pool / CN_ratio_up
        dstoich_upstreamN_pool_dC_pool_up = temp_real * dCN_ratio_up_dC_pool_up
        dstoich_upstreamN_pool_dN_pool_up = temp_real * dCN_ratio_up_dN_pool_up        

        Jacobian(iresN_pool_up,iresC_pool_up) = &
          Jacobian(iresN_pool_up,iresC_pool_up) - &
          ! scaled by negative one since it is a reactant 
!     revision to avoid division by 0 when upstream N -> 0 (line 727, 734, 735)
          (-1.d0) * (-1.d0) * rt_auxvar%immobile(ispecN_pool_up)/ &
                              rt_auxvar%immobile(ispecC_pool_up)* &
                              scaled_rate_const * N_inhibition

        Jacobian(iresN_pool_up,iresN_pool_up) = &
          Jacobian(iresN_pool_up,iresN_pool_up) - &
          ! scaled by negative one since it is a reactant 
!     revision to avoid division by 0 when upstream N -> 0 (line 727, 734, 735)
          (-1.d0) * scaled_rate_const * N_inhibition

        ! stoich_C = resp_frac * stoich_upstreamC_pool
        ! dstoichC_dC_pool_up = 0.
        ! dstoichC_dN_pool_up = 0

        ! nitrogen (stoichiometry a function of upstream C/N)
        ! stoich_N = stoich_upstreamN_pool - stoich_downstreamC_pool / &
        !                                    CN_ratio_down
        ! latter half is constant
        dstoichN_dC_pool_up = dstoich_upstreamN_pool_dC_pool_up
        dstoichN_dN_pool_up = dstoich_upstreamN_pool_dN_pool_up
        Jacobian(ires_N,iresC_pool_up) = Jacobian(ires_N,iresC_pool_up) - &
!     revision to avoid division by 0 when upstream N -> 0 (line 727, 734, 735)
          (-1.d0) * rt_auxvar%immobile(ispecN_pool_up)/ &
                    rt_auxvar%immobile(ispecC_pool_up)* &
                    scaled_rate_const * N_inhibition

        Jacobian(ires_N,iresN_pool_up) = Jacobian(ires_N,iresN_pool_up) - &
!     revision to avoid division by 0 when upstream N -> 0 (line 727, 734, 735)
          scaled_rate_const * N_inhibition
      endif
      
      ! carbon
      Jacobian(ires_C,iresC_pool_up) = Jacobian(ires_C,iresC_pool_up) - &
        stoich_C * drate
      ! nitrogen
      Jacobian(ires_N,iresC_pool_up) = Jacobian(ires_N,iresC_pool_up) - &
        stoich_N * drate
      if (use_N_inhibition) then
        Jacobian(ires_C,ires_N) = Jacobian(ires_C,ires_N) - & 
          stoich_C * drate_dN_inhibition
        Jacobian(ires_N,ires_N) = Jacobian(ires_N,ires_N) - &
          stoich_N * drate_dN_inhibition
      endif
    endif
  enddo
  
end subroutine CLM_CN_React

! ************************************************************************** !

subroutine CLM_CN_Destroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_clm_cn_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(clm_cn_reaction_type), pointer :: cur_reaction, prev_reaction
  
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    prev_pool => cur_pool
    cur_pool => cur_pool%next
    deallocate(prev_pool)
    nullify(prev_pool)
  enddo
  
  cur_reaction => this%reactions
  do
    if (.not.associated(cur_reaction)) exit
    prev_reaction => cur_reaction
    cur_reaction => cur_reaction%next
    deallocate(prev_reaction)
    nullify(prev_reaction)
  enddo
  
  call DeallocateArray(this%CN_ratio)
  call DeallocateArray(this%rate_constant)
  call DeallocateArray(this%respiration_fraction)
  call DeallocateArray(this%inhibition_constant)
  call DeallocateArray(this%upstream_pool_id)
  call DeallocateArray(this%downstream_pool_id)
  call DeallocateArray(this%pool_id_to_species_id)
  
end subroutine CLM_CN_Destroy

end module Reaction_Sandbox_CLM_CN_class
