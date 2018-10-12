module Reaction_Sandbox_Microbial_class
  ! this sandbox implements a microbial reaction with an electron donor, 
  ! an electron acceptor, a microbial mass, and a number of inhibitors in the 
  ! aqueous phase. Every reactant / product is assumed in the aqueous phase 
  ! except microbial mass, which is assumed to be immobile phase. 
  ! if a microbial mass is specified, the rate is of the first order with 
  ! respect to microbial mass
  ! if an electron donor or acceptor is specified, the rate is of Monod (C/(C+S))
  ! if an inhibitor is specified, the rate is inhibited by I/(C+I)
  ! inhibitors are neither reactant nor product

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Database_Aux_module

#ifdef CLM_PFLOTRAN
  use CLM_RspFuncs_module
#endif

  implicit none
  
  private
  
  type, public :: rate_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal                    :: value
    type(rate_type), pointer     :: next
  end type rate_type

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_microbial_type

    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscInt :: moisture_response_function
    PetscInt :: ph_response_function
    PetscReal :: fixed_ph

    character(len=MAXSTRINGLENGTH) :: str_reaction
    type(database_rxn_type), pointer ::dbase_rxn

    PetscReal :: rate_constant

    PetscReal :: residual_e_donor
    PetscReal :: residual_e_acceptor

    character(len=MAXWORDLENGTH) :: name_e_donor
    PetscReal :: half_saturation_e_donor
    PetscInt :: ispec_e_donor

    character(len=MAXWORDLENGTH) :: name_e_acceptor
    PetscReal :: half_saturation_e_acceptor
    PetscInt :: ispec_e_acceptor

    character(len=MAXWORDLENGTH) :: name_microbial_mass
    PetscInt :: ispec_microbial_mass

    ! allow multiple inhibition in terms of I/(C + I)
    PetscInt :: nInhibition

    type(rate_type), pointer :: Inhibition

    PetscInt, pointer :: ispec_inh(:)
    PetscReal, pointer :: inhibition_coef(:)

  contains
    procedure, public :: ReadInput => MicrobialRead
    procedure, public :: Setup => MicrobialSetup
    procedure, public :: Evaluate => MicrobialReact
    procedure, public :: Destroy => MicrobialDestroy
  end type reaction_sandbox_microbial_type

  public :: MicrobialCreate

contains

! ************************************************************************** !
!
! MicrobialCreate: Allocates Microbial reaction object.
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
function MicrobialCreate()

  implicit none
  
  class(reaction_sandbox_microbial_type), pointer :: MicrobialCreate

  allocate(MicrobialCreate)

#ifdef CLM_PFLOTRAN
  MicrobialCreate%temperature_response_function = -1
  MicrobialCreate%moisture_response_function = -1
  MicrobialCreate%ph_response_function = -1
#endif

  MicrobialCreate%Q10 = 1.5d0

  MicrobialCreate%str_reaction = ''

  nullify(MicrobialCreate%dbase_rxn)

  MicrobialCreate%rate_constant = 0.d0
  MicrobialCreate%residual_e_donor = 1.d-20
  MicrobialCreate%residual_e_acceptor = 1.d-20
  MicrobialCreate%fixed_ph = -1.0d0

  MicrobialCreate%name_e_donor = ''
  MicrobialCreate%half_saturation_e_donor = 1.0d-6
  MicrobialCreate%ispec_e_donor = -1

  MicrobialCreate%name_e_acceptor = ''
  MicrobialCreate%half_saturation_e_acceptor = 1.0d-6
  MicrobialCreate%ispec_e_acceptor = -1

  MicrobialCreate%name_microbial_mass = ''
  MicrobialCreate%ispec_microbial_mass = -1

  MicrobialCreate%nInhibition = 0

  nullify(MicrobialCreate%Inhibition)
  nullify(MicrobialCreate%ispec_inh)
  nullify(MicrobialCreate%inhibition_coef)
  nullify(MicrobialCreate%next)  

end function MicrobialCreate

function RateCreate()
  implicit none
  type(rate_type), pointer :: RateCreate
  allocate(RateCreate)  
  RateCreate%name = ''
  RateCreate%value = -1.0d0
  nullify(RateCreate%next)
end function RateCreate

recursive subroutine RateDestroy(rate)
  implicit none
  type(rate_type), pointer :: rate
  if (.not.associated(rate)) return
  call RateDestroy(rate%next)
  deallocate(rate)
  nullify(rate)
end subroutine RateDestroy

! ************************************************************************** !
!
! MicrobialRead:
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_microbial_type) :: this
  type(input_type), pointer              :: input
  type(option_type)                      :: option
  type(rate_type), pointer :: inhibition, inhibition_prev

  PetscReal :: tmp_real, rate_constant, turnover_time

  character(len=MAXWORDLENGTH) :: word

  nullify(inhibition_prev)

  turnover_time = 0.d0
  rate_constant = 0.d0

  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,Microbial')
    call StringToUpper(word)   

    select case(trim(word))
#ifdef CLM_PFLOTRAN
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,REACTION_SANDBOX,Microbial,' // &
            'TEMPERATURE RESPONSE FUNCTION.')
          call StringToUpper(word)   

          select case(trim(word))
            case('CLMCN')
              this%temperature_response_function = &
                TEMPERATURE_RESPONSE_FUNCTION_CLMCN
            case('Q10')
              this%temperature_response_function = &
                TEMPERATURE_RESPONSE_FUNCTION_Q10    
              call InputReadDouble(input,option,tmp_real)  
              call InputErrorMsg(input,option,'Q10', &
                'CHEMISTRY,REACTION_SANDBOX_Microbial,' // &
                'TEMPERATURE RESPONSE FUNCTION.')
                  this%Q10 = tmp_real
            case('DLEM')
              this%temperature_response_function = &
                TEMPERATURE_RESPONSE_FUNCTION_DLEM    
              call InputReadDouble(input,option,tmp_real)  
              call InputErrorMsg(input,option,'DLEM', &
                'CHEMISTRY,REACTION_SANDBOX_Microbial,' // &
                'TEMPERATURE RESPONSE FUNCTION')
              this%Q10 = tmp_real
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial,' // &
                'TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                trim(word) // ' not recognized - Valid keyword: "CLMCN","Q10" or "DLEM" '
              call printErrMsg(option)
          end select
        enddo 

      case('MOISTURE_RESPONSE_FUNCTION')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,REACTION_SANDBOX,Microbial,MOISTURE RESPONSE FUNCTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('CLMCN')
              this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_CLMCN
            case('DLEM')
              this%moisture_response_function = MOISTURE_RESPONSE_FUNCTION_DLEM    
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial,' // &
                'TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                trim(word) // ' not recognized  - Valid keyword: "CLMCN","DLEM"'
              call printErrMsg(option)
          end select
        enddo 

      case('PH_RESPONSE_FUNCTION')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,REACTION_SANDBOX,Microbial,PH RESPONSE FUNCTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('CENTURY')
              this%ph_response_function = PH_RESPONSE_FUNCTION_CENTURY    
            case('DLEM')
              this%ph_response_function = PH_RESPONSE_FUNCTION_DLEM    
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial,' // &
                'PH RESPONSE FUNCTION keyword: ' // &
                trim(word) // ' not recognized - Valid keyword: "CENTURY","DLEM".'
            call printErrMsg(option)
          end select
        enddo 

      case('FIXED_PH')
        call InputReadDouble(input,option,this%fixed_ph)
        call InputErrorMsg(input,option,'fixed ph', &
          'CHEMISTRY,REACTION_SANDBOX,Microbial,FIXED_PH.')
#endif

      case('REACTION')
        ! remainder of string should be the reaction equation
        this%str_reaction = trim(adjustl(input%buf))
        ! set flag for error message
        if (len_trim(this%str_reaction) < 2) input%ierr = 1
        call InputErrorMsg(input,option,'reaction string', &
          'CHEMISTRY,REACTION_SANDBOX_MICROBIAL_REACTION,REACTION')

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate constant', &
          'CHEMISTRY,REACTION_SANDBOX,Microbial,RATE_CONSTANT')

      case('ELECTRON_DONOR')
        call InputReadWord(input,option,this%name_e_donor,PETSC_TRUE)
        call InputErrorMsg(input,option,'Electron donor species name', &
          'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_DONOR')
        call InputReadDouble(input,option,this%half_saturation_e_donor)  
        call InputErrorMsg(input,option,'Electron donor half saturation coef', &
          'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_DONOR')

      case('ELECTRON_ACCEPTOR')
        call InputReadWord(input,option,this%name_e_acceptor,PETSC_TRUE)
        call InputErrorMsg(input,option,'Electron acceptor species name', &
          'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_ACCEPTOR')
        call InputReadDouble(input,option,this%half_saturation_e_acceptor)  
        call InputErrorMsg(input,option,'Electron acceptor half satur. coef', &
          'CHEMISTRY,REACTION_SANDBOX_Microbial,ELECTRON_ACCEPTOR')
      
      case('MICROBIAL_MASS')
        call InputReadWord(input,option,this%name_microbial_mass,PETSC_TRUE)
        call InputErrorMsg(input,option,'Microbial mass species name', &
          'CHEMISTRY,REACTION_SANDBOX_Microbial,MICROBIAL_MASS')
      case('RESIDUAL_ELECTRON_DONOR')
        call InputReadDouble(input,option,this%residual_e_donor)
        call InputErrorMsg(input,option,'e donor residual concentration', &
               'CHEMISTRY,REACTION_SANDBOX,Microbial,RESIDUAL_E_DONOR')
      case('RESIDUAL_ELECTRON_ACCEPTOR')
        call InputReadDouble(input,option,this%residual_e_acceptor)
        call InputErrorMsg(input,option,'e acceptor residual concentration', &
               'CHEMISTRY,REACTION_SANDBOX,Microbial,RESIDUAL_E_ACCEPTOR')
      case('INHIBITION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
          'CHEMISTRY,REACTION_SANDBOX_Microbial,INHIBITION')
        call InputReadDouble(input,option,tmp_real)  
        call InputErrorMsg(input,option,'inhibition constant', &
          'CHEMISTRY,REACTION_SANDBOX_Microbial,INHIBITION')
        inhibition => RateCreate()
        inhibition%name = word
        inhibition%value = tmp_real

        if (.not.associated(this%inhibition)) then
          this%inhibition => inhibition
        else
          inhibition_prev%next => inhibition
        endif
        inhibition_prev => inhibition
        nullify(inhibition)
        this%nInhibition = this%nInhibition + 1

      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial keyword: ' // &
          trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo

end subroutine MicrobialRead

! ************************************************************************** !
!
! MicrobialSetup: 
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialSetup(this,reaction,option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module
  use Reaction_Database_Aux_module

  implicit none

  class(reaction_sandbox_microbial_type) :: this
  type(reaction_type)                    :: reaction
  type(option_type)                      :: option
  type(rate_type), pointer               :: cur_rate

  PetscInt :: i, icount, ispec

  character(len=MAXWORDLENGTH) :: word

  ! parse reaction
  this%dbase_rxn => DatabaseRxnCreateFromRxnString(this%str_reaction,  &
                           reaction%naqcomp, reaction%offset_aqueous,  &
                           reaction%primary_species_names,             &
                           reaction%nimcomp, reaction%offset_immobile, &
                           reaction%immobile%names,                    &
                           PETSC_TRUE, option)

  ! electron donor
  if (trim(this%name_e_donor) /= '') then
    this%ispec_e_donor = GetPrimarySpeciesIDFromName(  &
      this%name_e_donor,reaction,PETSC_FALSE,option)
  endif

  if (trim(this%name_e_acceptor) /= '') then
    this%ispec_e_acceptor = GetPrimarySpeciesIDFromName(  &
      this%name_e_acceptor,reaction,PETSC_FALSE,option)
  endif

  if (trim(this%name_microbial_mass) /= '') then
    this%ispec_microbial_mass = GetImmobileSpeciesIDFromName( &
      this%name_microbial_mass,reaction%immobile,PETSC_FALSE,option) 
  endif

  ! inhibition rate terms
  if (this%nInhibition >= 1) then
    allocate(this%ispec_inh(this%nInhibition))
    allocate(this%inhibition_coef(this%nInhibition))
 
    cur_rate => this%Inhibition
    icount = 1
    do
      if (.not.associated(cur_rate)) exit
      ispec = GetPrimarySpeciesIDFromName(cur_rate%name,reaction,PETSC_FALSE, &
        option)
      if (ispec > 0) then
        this%ispec_inh(icount) = ispec
      else
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,Microbial check: ' // &
          'Inhibition term ' // trim(cur_rate%name)// &
           ' is not an aqueous species.'
        call printErrMsg(option)
      endif
      this%inhibition_coef(icount) = cur_rate%value
      cur_rate => cur_rate%next
      icount = icount + 1
    enddo 
  endif   

end subroutine MicrobialSetup

! ************************************************************************** !
!
! MicrobialReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)

  use Option_module
  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Material_Aux_class, only : material_auxvar_type
  use Utility_module, only : DeallocateArray
  
#ifdef CLM_PFLOTRAN
#include "petsc/finclude/petscvec.h"
  use petscvec
  use clm_pflotran_interface_data
#endif
  
  implicit none

  class(reaction_sandbox_microbial_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: Lwater
  PetscReal :: rate, drate, drate_ddonor, drate_dacceptor, drate_dbiomass
  PetscReal :: tmp_real, conc, activity 

  PetscReal :: f_t
  PetscReal :: f_w
  PetscReal :: f_ph
  PetscReal :: ph

#ifdef CLM_PFLOTRAN
  PetscReal :: tc
  PetscReal :: theta
#endif

  PetscReal, pointer :: drate_dinh(:)

  PetscInt :: local_id
  PetscInt :: i, j, icomp, ires
  PetscErrorCode :: ierr

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

! temperature response function 
#ifdef CLM_PFLOTRAN 
  tc = global_auxvar%temp
  f_t = GetTemperatureResponse(tc, this%temperature_response_function, this%Q10) 
#else
  f_t = 1.0d0
#endif

  ! moisture response function 
#ifdef CLM_PFLOTRAN
  local_id = option%iflag ! temporary measure suggested by Glenn
  theta = global_auxvar%sat(1) * porosity 
  f_w = GetMoistureResponse(theta, local_id, this%moisture_response_function)
#else
  f_w = 1.0d0
#endif

#ifdef CLM_PFLOTRAN
  if (this%ph_response_function > 0) then
    if (this%fixed_ph > 0.0d0) then
      ph = this%fixed_ph
    else
      if (reaction%species_idx%h_ion_id > 0) then
        ph = -log10(rt_auxvar%pri_molal(reaction%species_idx%h_ion_id) * &
              rt_auxvar%pri_act_coef(reaction%species_idx%h_ion_id))
      elseif (reaction%species_idx%h_ion_id < 0) then
        ph = &
          -log10(rt_auxvar%sec_molal(abs(reaction%species_idx%h_ion_id))* &
                 rt_auxvar%sec_act_coef(abs(reaction%species_idx%h_ion_id)))
      endif
    endif
    f_ph = GetpHResponse(ph, this%ph_response_function)  
  endif 
#else
  f_ph = 1.0d0
#endif 
  
  if (f_t < 1.0d-20 .or. f_w < 1.0d-20 .or. f_ph < 1.0d-20) then
     return
  endif

  ! residual
  Lwater = porosity * global_auxvar%sat(iphase) * volume * 1.d3

  rate = this%rate_constant * Lwater * f_t * f_w * f_ph

  if (compute_derivative) then 
    drate_ddonor    = rate
    drate_dacceptor = rate
    drate_dbiomass  = rate

    if (this%nInhibition >= 1) then
      allocate(drate_dinh(this%nInhibition))
    endif
 
    do i = 1, this%nInhibition
      drate_dinh(i) = rate
    enddo
  endif

  ! biomass
  icomp = this%ispec_microbial_mass
  if (icomp > 0) then
    conc = rt_auxvar%immobile(icomp)
    rate = rate * conc 
    if (compute_derivative) then 
      drate_ddonor = drate_ddonor * conc
      drate_dacceptor = drate_dacceptor * conc
      do i = 1, this%nInhibition
        drate_dinh(i) = drate_dinh(i) * conc
      enddo
    endif
  endif

  ! electron donor
  icomp = this%ispec_e_donor
  if (icomp > 0) then
    activity = (rt_auxvar%pri_molal(icomp) - this%residual_e_donor) &
             * rt_auxvar%pri_act_coef(icomp)
    tmp_real = activity / (activity + this%half_saturation_e_donor) 
    rate = rate * tmp_real
    if (compute_derivative) then 
      drate_dbiomass = drate_dbiomass * tmp_real
      drate_ddonor = drate_ddonor * this%half_saturation_e_donor / &
                     (activity + this%half_saturation_e_donor) / &
                     (activity + this%half_saturation_e_donor) * &
                     rt_auxvar%pri_act_coef(icomp)
      drate_dacceptor = drate_dacceptor * tmp_real
      do i = 1, this%nInhibition
        drate_dinh(i) = drate_dinh(i) * tmp_real
      enddo
    endif
  endif

  ! electron acceptor
  icomp = this%ispec_e_acceptor
  if (icomp > 0) then
    activity = (rt_auxvar%pri_molal(icomp) - this%residual_e_acceptor) &
             * rt_auxvar%pri_act_coef(icomp)
    tmp_real = activity /(activity + this%half_saturation_e_acceptor) 
    rate = rate * tmp_real
    if (compute_derivative) then 
      drate_dbiomass = drate_dbiomass * tmp_real
      drate_ddonor = drate_ddonor * tmp_real
      drate_dacceptor = drate_dacceptor * this%half_saturation_e_acceptor / &
                   (activity + this%half_saturation_e_acceptor) / &
                   (activity + this%half_saturation_e_acceptor) * &
                   rt_auxvar%pri_act_coef(icomp)
      do i = 1, this%nInhibition
        drate_dinh(i) = drate_dinh(i) * tmp_real
      enddo
    endif 
  endif

  ! inhibition term
  do i = 1, this%nInhibition
    icomp = this%ispec_inh(i)
    if (icomp <= 0) cycle
    activity = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp) 
    tmp_real = this%inhibition_coef(i)/(activity + this%inhibition_coef(i))  
    rate = rate * tmp_real  
    if (compute_derivative) then 
      drate_dbiomass = drate_dbiomass * tmp_real
      drate_ddonor = drate_ddonor * tmp_real
      drate_dacceptor = drate_dacceptor * tmp_real
      do j = 1, this%nInhibition
        if (i == j) then
          drate_dinh(j) = drate_dinh(j) * (-1.0d0) * this%inhibition_coef(j) &
                        / (activity + this%inhibition_coef(j)) &
                        / (activity + this%inhibition_coef(j)) &
                        * rt_auxvar%pri_act_coef(icomp)
        else 
          drate_dinh(j) = drate_dinh(j) * tmp_real
        endif
      enddo
    endif
  enddo
 
  ! Residual 
  do i = 1, this%dbase_rxn%nspec
    if (this%dbase_rxn%spec_ids(i) == this%ispec_microbial_mass) then
      ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
    else
      ires = this%dbase_rxn%spec_ids(i)
    endif
    Residual(ires) = Residual(ires) - this%dbase_rxn%stoich(i) * rate 
  enddo

  ! Jacobian
  if (.not.compute_derivative) return 

  ! with respect to biomass
  icomp = this%ispec_microbial_mass
  if (icomp > 0) then
    do i = 1, this%dbase_rxn%nspec
      if (this%dbase_rxn%spec_ids(i) == this%ispec_microbial_mass) then
        ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
      else
        ires = this%dbase_rxn%spec_ids(i)
      endif
      Jacobian(ires, icomp + reaction%offset_immobile) = &
        Jacobian(ires, icomp + reaction%offset_immobile) &
        - this%dbase_rxn%stoich(i) * drate_dbiomass 
    enddo
  endif

  ! with respect to electron donor     
  icomp = this%ispec_e_donor
  if (icomp > 0) then
    do i = 1, this%dbase_rxn%nspec
      if (this%dbase_rxn%spec_ids(i) == this%ispec_microbial_mass) then
        ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
      else
        ires = this%dbase_rxn%spec_ids(i)
      endif
      Jacobian(ires, icomp) = Jacobian(ires,icomp) &
                            - this%dbase_rxn%stoich(i) * drate_ddonor 
     enddo
  endif

  ! with respect to electron acceptor
  icomp = this%ispec_e_acceptor
  if (icomp > 0) then
    do i = 1, this%dbase_rxn%nspec
      if (this%dbase_rxn%spec_ids(i) == this%ispec_microbial_mass) then
        ires = this%dbase_rxn%spec_ids(i) + reaction%offset_immobile
      else
        ires = this%dbase_rxn%spec_ids(i)
      endif
      Jacobian(ires, icomp) = Jacobian(ires,icomp) &
                            - this%dbase_rxn%stoich(i) * drate_dacceptor 
    enddo
  endif

  ! with respect to inhibitors
  do i = 1, this%nInhibition
    icomp = this%ispec_inh(i)
    if (icomp <= 0) cycle
    do j = 1, this%dbase_rxn%nspec
      if (this%dbase_rxn%spec_ids(j) == this%ispec_microbial_mass) then
        ires = this%dbase_rxn%spec_ids(j) + reaction%offset_immobile
      else
        ires = this%dbase_rxn%spec_ids(j)
      endif
      Jacobian(ires, icomp) = Jacobian(ires,icomp) &
                            - this%dbase_rxn%stoich(j) * drate_dinh(i) 
    enddo
  enddo

  if (this%nInhibition >= 1) then
    call DeallocateArray(drate_dinh)
  endif

end subroutine MicrobialReact

! ************************************************************************** !
!
! MicrobialDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Guoping Tang
! date: 11/25/13
!
! ************************************************************************** !
subroutine MicrobialDestroy(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(reaction_sandbox_microbial_type) :: this  

  type(rate_type), pointer :: cur_rate, prev_rate

  call DatabaseRxnDestroy(this%dbase_rxn)    

  cur_rate => this%Inhibition
  do 
    if(.not.associated(cur_rate)) exit
       prev_rate => cur_rate
       cur_rate => cur_rate%next
       deallocate(prev_rate)
       nullify(prev_rate)
  enddo

  call DeallocateArray(this%ispec_inh) 

  call DeallocateArray(this%inhibition_coef) 
 
end subroutine MicrobialDestroy

end module Reaction_Sandbox_Microbial_class
