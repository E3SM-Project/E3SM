module Miscible_Aux_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

type, public :: Miscible_auxvar_elem_type
   PetscReal :: pres
    PetscReal :: temp
    PetscReal , pointer :: sat(:)
    PetscReal , pointer :: den(:)
    PetscReal , pointer :: avgmw(:)
    PetscReal , pointer :: vis(:)
    PetscReal , pointer :: h(:)
    PetscReal , pointer :: u(:)
    PetscReal , pointer :: pc(:)
    PetscReal , pointer :: kvr(:)
    PetscReal , pointer :: xmol(:)
    PetscReal , pointer :: diff(:)
    PetscReal , pointer :: hysdat(:)
    PetscReal :: zco2
 end type Miscible_auxvar_elem_type

  type, public :: Miscible_auxvar_type
    type(Miscible_auxvar_elem_type), pointer :: auxvar_elem(:) 
  end type Miscible_auxvar_type
  
  type, public :: Miscible_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckwet(:)
    PetscReal, pointer :: ckdry(:)
    PetscReal, pointer :: sir(:,:)
  end type Miscible_parameter_type
    
  type, public :: Miscible_type
     PetscInt :: n_zero_rows
     PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
     PetscBool :: auxvars_up_to_date
     PetscBool :: inactive_cells_exist
     PetscInt :: num_aux, num_aux_bc, num_aux_ss
     type(Miscible_parameter_type), pointer :: Miscible_parameter
     type(Miscible_auxvar_type), pointer :: auxvars(:)
     type(Miscible_auxvar_type), pointer :: auxvars_bc(:)
     type(Miscible_auxvar_type), pointer :: auxvars_ss(:)
     PetscReal, pointer :: Resold_AR(:,:)
     PetscReal, pointer :: Resold_BC(:,:)
     PetscReal, pointer :: Resold_FL(:,:)
     PetscReal, pointer :: delx(:,:)
  end type Miscible_type

  public :: MiscibleAuxCreate, MiscibleAuxDestroy, &
            MiscibleAuxVarCompute_NINC, MiscibleAuxVarCompute_WINC,&
            MiscibleAuxVarInit, MiscibleAuxVarCopy

contains

! ************************************************************************** !

function MiscibleAuxCreate()
  ! 
  ! MiscibleAuxVarCreate: Allocate and initialize auxiliary object
  ! 
  ! Author: Chuan Lu
  ! Date: 02/27/08
  ! 

  use Option_module

  implicit none
  
  type(Miscible_type), pointer :: MiscibleAuxCreate
  
  type(Miscible_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_zero_rows = 0
  allocate(aux%Miscible_parameter)
  nullify(aux%Miscible_parameter%sir)
  nullify(aux%Miscible_parameter%ckwet)
  nullify(aux%Miscible_parameter%dencpr)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
  nullify(aux%Resold_AR)
  nullify(aux%Resold_BC)
  nullify(aux%Resold_FL)
  nullify(aux%delx)

  MiscibleAuxCreate => aux
  
end function MiscibleAuxCreate

! ************************************************************************** !

subroutine MiscibleAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Chuan Lu
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(Miscible_auxvar_type) :: auxvar
  type(option_type) :: option

  PetscInt :: var_elem_size, var_node_size
  PetscInt :: nvar 

  allocate(auxvar%auxvar_elem(0 : option%nflowdof))
  allocate(auxvar%auxvar_elem(0)%hysdat(4))
 
  do nvar = 0, option%nflowdof
     allocate ( auxvar%auxvar_elem(nvar)%sat(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%den(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%avgmw(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%vis(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%h(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%u(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%pc(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%kvr(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%xmol(option%nphase*option%nflowspec))
     allocate ( auxvar%auxvar_elem(nvar)%diff(option%nphase*option%nflowspec))
     if (nvar>0)&
     auxvar%auxvar_elem(nvar)%hysdat => auxvar%auxvar_elem(0)%hysdat

     auxvar%auxvar_elem(nvar)%pres = 0.d0
     auxvar%auxvar_elem(nvar)%temp = 0.d0
     auxvar%auxvar_elem(nvar)%sat = 0.d0
     auxvar%auxvar_elem(nvar)%den = 0.d0
     auxvar%auxvar_elem(nvar)%avgmw = 0.d0
     auxvar%auxvar_elem(nvar)%vis = 0.d0
     auxvar%auxvar_elem(nvar)%h = 0.d0
     auxvar%auxvar_elem(nvar)%u = 0.d0
     auxvar%auxvar_elem(nvar)%pc = 0.d0
     auxvar%auxvar_elem(nvar)%kvr = 0.d0
     auxvar%auxvar_elem(nvar)%xmol = 0.d0
     auxvar%auxvar_elem(nvar)%diff = 0.d0
  enddo

end subroutine MiscibleAuxVarInit

! ************************************************************************** !

subroutine MiscibleAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Chuan Lu
  ! Date: 10/13/0
  ! 

  use Option_module

  implicit none
  
  type(Miscible_auxvar_elem_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%avgmw = auxvar%avgmw
  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
  auxvar2%kvr = auxvar%kvr
!  auxvar2%xmol = auxvar%xmol
!  auxvar2%diff = auxvar%diff

end subroutine MiscibleAuxVarCopy

! ************************************************************************** !

subroutine Water_glycol_density(y,p,dkg)
  ! 
  ! Computes water-propylene glycol mixture density
  ! 
  ! Author: Chuan Lu
  ! Date: 12/12/11
  ! 
  implicit none
  PetscReal y, p ! water mass fraction
  PetscReal dkg

  dkg = (((0.0806d0*y-0.203d0)*y + 0.0873d0)*y + 1.0341d0) * 1.d3
! dkg = (4.49758d-10*y +(1.d0-y)*5.d-10)*(p-1.01325d5) + dkg
! dkg = dkg * 1.d3  ! convert g/cm^3 to kg/m^3

  dkg = dkg * (1+(4.49758d-10*y +(1.d0-y)*5.d-10)*(p-1.01325d5))
 
end subroutine Water_glycol_density

! ************************************************************************** !

subroutine MiscibleAuxVarCompute_NINC(x,auxvar,global_auxvar, &
             fluid_properties,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! No increments
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Option_module
  use Global_Aux_module  
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(Miscible_auxvar_elem_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  
  PetscReal :: x(option%nflowdof)
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: p, t, temp, p2, err
  PetscReal :: henry,lngamco2
  PetscReal :: dg, dddp, dddt
  PetscReal :: fg, dfgdp, dfgdt, xphi
  PetscReal :: eng,hg, dhdp, dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase), u(option%nphase), kr(option%nphase)
  PetscReal :: tk, visw 
  PetscReal :: denw, yh2o, yppg
  PetscReal :: tmp 
  PetscInt :: iflag  
  
  auxvar%sat = 0.d0
  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%den = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%pc = 0.d0
  auxvar%kvr = 0.d0
  auxvar%xmol = 0.d0
! auxvar%diff = 0.d0
  kr = 0.d0

  auxvar%pres = x(1)

! auxvar%xmol(2) = x(2)
  auxvar%xmol(2) = exp(x(2))
! auxvar%xmol(2) = (atan(x(2))/3.1416*2+1)/2
  auxvar%xmol(1) = 1.D0 - auxvar%xmol(2)

! Glycol-Water mixture formula weight (kg/kmol)
  auxvar%avgmw(1) = auxvar%xmol(1)*FMWH2O + auxvar%xmol(2)*FMWGLYC
  
! Mass fraction water
  yh2o = auxvar%xmol(1)*FMWH2O/auxvar%avgmw(1)
  
  call Water_glycol_density(yh2o, auxvar%pres, denw)
  
  auxvar%den(1) = denw/auxvar%avgmw(1)
  
! Glycol-Water mixture viscosity (yh2o mass fraction water)
  yppg = 1.d0-yh2o
  visw = 10.d0**(1.6743d0*yppg-0.0758d0) * 1.0d-3 ! centipoise to Pa s.

  auxvar%vis(1) = visw
  
  auxvar%sat(1) = 1.d0
  auxvar%kvr(1) = 1.d0/visw
  auxvar%h(1) = denw*4.18d-3*global_auxvar%temp
  
! Glycol-Water mixture diffusivity (yh2o mass fraction water)
  auxvar%diff(2) = ((((-4.021d0*yh2o + 9.1181d0)*yh2o - 5.9703d0)*yh2o &
     + 0.4043d0)*yh2o + 0.5687d0) * 1.d-9 ! m^2/s
  auxvar%diff(1) = auxvar%diff(2)

! auxvar%diff(1:option%nflowspec) = fluid_properties%diffusion_coefficient

end subroutine MiscibleAuxVarCompute_NINC

! ************************************************************************** !

subroutine MiscibleAuxVarCompute_WINC(x,delx,auxvar,global_auxvar, &
                                    fluid_properties,option)

  use Option_module
  use Global_Aux_module
  
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(Miscible_auxvar_elem_type) :: auxvar(1:option%nflowdof)
  type(global_auxvar_type) :: global_auxvar

  PetscReal :: x(option%nflowdof), xx(option%nflowdof), delx(option%nflowdof)
  PetscInt :: idof 
  
  do idof = 1, option%nflowdof
    xx = x; xx(idof) = x(idof) + delx(idof)

!   print *,'Winc: ',idof,x(idof),xx(idof),delx(idof)

! ***   note: var_node here starts from 1 to option%flowdof ***
    call  MiscibleAuxVarCompute_NINC(xx,auxvar(idof),global_auxvar, &
      fluid_properties,option)
  enddo

end subroutine MiscibleAuxVarCompute_WINC

! ************************************************************************** !

subroutine MiscibleAuxVarElemDestroy(auxvar_elem)
  ! 
  ! Deallocates a mphase auxiliary elment object
  ! 
  implicit none

  type(miscible_auxvar_elem_type) :: auxvar_elem

  if (associated(auxvar_elem%xmol)) deallocate(auxvar_elem%xmol)
  nullify(auxvar_elem%xmol)
  if (associated(auxvar_elem%diff)) deallocate(auxvar_elem%diff)
  nullify(auxvar_elem%diff)
  if (associated(auxvar_elem%pc)) deallocate(auxvar_elem%pc)
  nullify(auxvar_elem%pc)
  if (associated(auxvar_elem%sat)) deallocate(auxvar_elem%sat)
  nullify(auxvar_elem%sat)
  if (associated(auxvar_elem%u)) deallocate(auxvar_elem%u)
  nullify(auxvar_elem%u)
  if (associated(auxvar_elem%h)) deallocate(auxvar_elem%h)
  nullify(auxvar_elem%h)
  if (associated(auxvar_elem%den)) deallocate(auxvar_elem%den)
  nullify(auxvar_elem%den)
  if (associated(auxvar_elem%den)) deallocate(auxvar_elem%vis)
  nullify(auxvar_elem%vis)
  if (associated(auxvar_elem%avgmw)) deallocate(auxvar_elem%avgmw)
  nullify(auxvar_elem%avgmw)

end subroutine MiscibleAuxVarElemDestroy

! ************************************************************************** !

subroutine MiscibleAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a miscible auxiliary object
  ! 
  implicit none

  type(miscible_auxvar_type) :: auxvar

  PetscInt :: ielem

  ! subtract 1 since indexing from 0
  if (associated(auxvar%auxvar_elem)) then
    do ielem = 0, size(auxvar%auxvar_elem) - 1 
      call MiscibleAuxVarElemDestroy(auxvar%auxvar_elem(ielem))
    enddo
    deallocate(auxvar%auxvar_elem)
    nullify(auxvar%auxvar_elem)
  endif

end subroutine MiscibleAuxVarDestroy

! ************************************************************************** !

subroutine MiscibleAuxDestroy(aux)
  ! 
  ! Deallocates a miscible auxiliary object
  ! 
  implicit none

  type(miscible_type), pointer :: aux

  PetscInt :: iaux
  
  if (.not.associated(aux)) return

  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call MiscibleAuxVarDestroy(aux%auxvars(iaux))
    enddo
    deallocate(aux%auxvars)
    nullify(aux%auxvars)
  endif
  
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call MiscibleAuxVarDestroy(aux%auxvars_bc(iaux))
    enddo
    deallocate(aux%auxvars_bc)
    nullify(aux%auxvars_bc)
  endif

#if 0
  if (associated(aux%auxvars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call MiscibleAuxVarDestroy(aux%auxvars_ss(iaux))
    enddo
    deallocate(aux%auxvars_ss)
    nullify(aux%auxvars_ss)
  endif
#endif

  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%miscible_parameter)) then
    if (associated(aux%miscible_parameter%dencpr)) deallocate(aux%miscible_parameter%dencpr)
    nullify(aux%miscible_parameter%dencpr)
    if (associated(aux%miscible_parameter%ckwet)) deallocate(aux%miscible_parameter%ckwet)
    nullify(aux%miscible_parameter%ckwet)
    if (associated(aux%miscible_parameter%sir)) deallocate(aux%miscible_parameter%sir)
    nullify(aux%miscible_parameter%sir)
    deallocate(aux%miscible_parameter)
  endif
  nullify(aux%miscible_parameter)
! if (associated(aux%res_old_AR)) deallocate(aux%res_old_AR)
! if (associated(aux%res_old_FL)) deallocate(aux%res_old_FL)
  if (associated(aux%delx)) deallocate(aux%delx)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine MiscibleAuxDestroy

end module Miscible_Aux_module
