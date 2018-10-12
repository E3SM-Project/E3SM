module Global_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  type, public :: global_auxvar_type
    PetscInt :: istate
    PetscReal :: temp
    PetscReal, pointer :: pres(:)
    PetscReal, pointer :: pres_store(:,:)
    PetscReal, pointer :: temp_store(:)
    PetscReal, pointer :: sat(:)
    PetscReal, pointer :: sat_store(:,:)
    PetscReal, pointer :: den(:)  ! kmol/m^3
    PetscReal, pointer :: den_kg(:) ! kg/m^3
    PetscReal, pointer :: den_store(:,:)
    PetscReal, pointer :: den_kg_store(:,:)
    PetscReal, pointer :: fugacoeff(:)
    PetscReal, pointer :: fugacoeff_store(:,:)
    PetscReal, pointer :: m_nacl(:)
    PetscReal, pointer :: xmass(:)
    PetscReal, pointer :: mass_balance(:,:) ! kg
    PetscReal, pointer :: mass_balance_delta(:,:) ! kmol
    PetscReal, pointer :: reaction_rate(:)
    PetscReal, pointer :: reaction_rate_store(:)
    PetscReal, pointer :: dphi(:,:) !geh: why here?
!geh    PetscReal :: scco2_eq_logK ! SC CO2
  end type global_auxvar_type
  
  type, public :: global_type
    PetscReal :: time_t, time_tpdt
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(global_auxvar_type), pointer :: auxvars(:)
    type(global_auxvar_type), pointer :: auxvars_bc(:)
    type(global_auxvar_type), pointer :: auxvars_ss(:)
  end type global_type
  
  interface GlobalAuxVarDestroy
    module procedure GlobalAuxVarSingleDestroy
    module procedure GlobalAuxVarArrayDestroy
  end interface GlobalAuxVarDestroy
  
  public :: GlobalAuxCreate, GlobalAuxDestroy, &
            GlobalAuxVarInit, GlobalAuxVarCopy, &
            GlobalAuxVarDestroy, GlobalAuxVarStrip

contains

! ************************************************************************** !

function GlobalAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(global_type), pointer :: GlobalAuxCreate
  
  type(global_type), pointer :: aux

  allocate(aux) 
  aux%time_t = 0.d0
  aux%time_tpdt = 0.d0
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  GlobalAuxCreate => aux
  
end function GlobalAuxCreate

! ************************************************************************** !

subroutine GlobalAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: auxvar
  type(option_type) :: option

  PetscInt :: nphase
  
  auxvar%istate = 0
  auxvar%temp = 0.d0

  ! nullify everthing to begin with and allocate later
  nullify(auxvar%pres)
  nullify(auxvar%sat)
  nullify(auxvar%den)
  nullify(auxvar%den_kg)
  nullify(auxvar%pres_store)
  nullify(auxvar%temp_store)
  nullify(auxvar%sat_store)
  nullify(auxvar%den_store)
  nullify(auxvar%den_kg_store)
  nullify(auxvar%fugacoeff)
  nullify(auxvar%fugacoeff_store)
  nullify(auxvar%m_nacl)
  nullify(auxvar%xmass)
  nullify(auxvar%reaction_rate)
  nullify(auxvar%reaction_rate_store)
  nullify(auxvar%mass_balance)
  nullify(auxvar%mass_balance_delta)
  nullify(auxvar%dphi)

  nphase = max(option%nphase,option%transport%nphase)

  if (option%nflowdof > 0) then
    allocate(auxvar%den(nphase))
    auxvar%den = 0.d0
  endif
  allocate(auxvar%pres(nphase))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den_kg(nphase))
  auxvar%den_kg = 0.d0

  ! need these for reactive transport only if if flow if computed
  if (option%nflowdof > 0 .and. option%ntrandof > 0) then
    allocate(auxvar%sat_store(nphase,TWO_INTEGER))
    auxvar%sat_store = 0.d0
    allocate(auxvar%den_kg_store(nphase,TWO_INTEGER))
    auxvar%den_kg_store = 0.d0
  endif
 
  select case(option%iflowmode)
    case(TH_MODE)
      allocate(auxvar%pres_store(nphase,TWO_INTEGER))
      auxvar%pres_store = 0.d0
      allocate(auxvar%temp_store(TWO_INTEGER))
      auxvar%temp_store = 0.d0
      allocate(auxvar%den_store(nphase,TWO_INTEGER))
      auxvar%den_store = 0.d0
      nullify(auxvar%xmass)
      nullify(auxvar%fugacoeff)
      nullify(auxvar%fugacoeff_store)
      nullify(auxvar%m_nacl)
      nullify(auxvar%reaction_rate)
      nullify(auxvar%reaction_rate_store)  
    case default
  end select
  
  if (option%flow%density_depends_on_salinity) then
    allocate(auxvar%m_nacl(ONE_INTEGER))
    auxvar%m_nacl = 0.d0
  endif
  
  if (option%iflag /= 0 .and. option%compute_mass_balance_new) then
    allocate(auxvar%mass_balance(option%nflowspec,nphase))
    auxvar%mass_balance = 0.d0
    allocate(auxvar%mass_balance_delta(option%nflowspec,nphase))
    auxvar%mass_balance_delta = 0.d0
  endif
  
end subroutine GlobalAuxVarInit

! ************************************************************************** !

subroutine GlobalAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(global_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate = auxvar%istate
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
!  auxvar2%dphi = auxvar%dphi
  
  if (associated(auxvar2%reaction_rate)) then
    auxvar2%reaction_rate = auxvar%reaction_rate
  endif
  
  if (associated(auxvar2%m_nacl)) then
    auxvar2%m_nacl = auxvar%m_nacl
  endif

  if (associated(auxvar2%fugacoeff)) then
    auxvar2%fugacoeff = auxvar%fugacoeff  
  endif
  if (associated(auxvar2%xmass)) then
    auxvar2%xmass = auxvar%xmass  
  endif
  if (associated(auxvar2%pres_store)) then
    auxvar2%pres_store = auxvar%pres_store  
  endif
  if (associated(auxvar2%den_store)) then
    auxvar2%den_store = auxvar%den_store  
  endif
  if (associated(auxvar2%sat_store)) then
    auxvar2%sat_store = auxvar%sat_store  
  endif
  if (associated(auxvar2%den_kg_store)) then
    auxvar2%den_kg_store = auxvar%den_kg_store
  endif
  if (associated(auxvar2%temp_store)) then
    auxvar2%temp_store = auxvar%temp_store  
  endif
  if (associated(auxvar2%fugacoeff_store)) then
    auxvar2%fugacoeff_store = auxvar%fugacoeff_store  
  endif

  !geh: here we have to check on both as mass_balance often exists for bcs and
  !     src/sinks but not regular cells.
  if (associated(auxvar%mass_balance) .and. &
      associated(auxvar2%mass_balance)) then
    auxvar2%mass_balance = auxvar%mass_balance
    auxvar2%mass_balance_delta = auxvar%mass_balance_delta
  endif

end subroutine GlobalAuxVarCopy

! ************************************************************************** !

subroutine GlobalAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(global_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call GlobalAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine GlobalAuxVarSingleDestroy

! ************************************************************************** !

subroutine GlobalAuxVarArrayDestroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  implicit none

  type(global_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call GlobalAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine GlobalAuxVarArrayDestroy

! ************************************************************************** !

subroutine GlobalAuxVarStrip(auxvar)
  ! 
  ! Deallocates all members of single auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/10/12
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(global_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)
  call DeallocateArray(auxvar%sat)
  call DeallocateArray(auxvar%den)
  call DeallocateArray(auxvar%fugacoeff)
  call DeallocateArray(auxvar%den_kg)
  call DeallocateArray(auxvar%m_nacl)
  call DeallocateArray(auxvar%xmass)
  call DeallocateArray(auxvar%reaction_rate)
  call DeallocateArray(auxvar%dphi)

  call DeallocateArray(auxvar%pres_store)
  call DeallocateArray(auxvar%temp_store)
  call DeallocateArray(auxvar%fugacoeff_store)
  call DeallocateArray(auxvar%sat_store)
  call DeallocateArray(auxvar%den_kg_store)
  call DeallocateArray(auxvar%den_store)
  call DeallocateArray(auxvar%reaction_rate_store)
  
  call DeallocateArray(auxvar%mass_balance)
  call DeallocateArray(auxvar%mass_balance_delta)

end subroutine GlobalAuxVarStrip

! ************************************************************************** !

subroutine GlobalAuxDestroy(aux)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(global_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call GlobalAuxVarDestroy(aux%auxvars)
  call GlobalAuxVarDestroy(aux%auxvars_bc)
  call GlobalAuxVarDestroy(aux%auxvars_ss)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GlobalAuxDestroy

end module Global_Aux_module
