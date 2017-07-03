module Immis_Aux_module

  use Mphase_pckr_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

type, public :: Immis_auxvar_elem_type
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
! PetscReal , pointer :: xmol(:)
! PetscReal , pointer :: diff(:)
  PetscReal , pointer :: hysdat(:)
  PetscReal :: zco2
!    PetscReal :: vis
!    PetscReal :: dvis_dp
!    PetscReal :: kr
!    PetscReal :: dkr_dp
 end type Immis_auxvar_elem_type

  type, public :: Immis_auxvar_type
    
    type(Immis_auxvar_elem_type), pointer :: auxvar_elem(:) 
#if 0
    PetscReal , pointer :: davgmw_dx(:)
    PetscReal , pointer :: dden_dp(:)
    PetscReal , pointer :: dden_dt(:)
    PetscReal , pointer :: dden_dx(:)
    PetscReal , pointer :: dkvr_dp(:)
    PetscReal , pointer :: dkvr_dt(:)
    PetscReal , pointer :: dkvr_ds(:)
    PetscReal , pointer :: dkvr_dx(:)
    PetscReal , pointer :: dh_dp(:)
    PetscReal , pointer :: dh_dt(:)
    PetscReal , pointer :: dh_dx(:)
    PetscReal , pointer :: du_dp(:)
    PetscReal , pointer :: du_dt(:)
    PetscReal , pointer :: du_dx(:)
#endif
  end type Immis_auxvar_type
  
  type, public :: Immis_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckwet(:)
    PetscReal, pointer :: ckdry(:)
    PetscReal, pointer :: sir(:,:)
  end type Immis_parameter_type
    
  type, public :: Immis_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(Immis_parameter_type), pointer :: immis_parameter
    type(Immis_auxvar_type), pointer :: auxvars(:)
    type(Immis_auxvar_type), pointer :: auxvars_bc(:)
    type(Immis_auxvar_type), pointer :: auxvars_ss(:)

    PetscReal, pointer :: res_old_AR(:,:), res_old_FL(:,:), delx(:,:)
  end type Immis_type

  

  public :: ImmisAuxCreate, ImmisAuxDestroy, &
            ImmisAuxVarCompute_NINC, ImmisAuxVarCompute_WINC,&
            ImmisAuxVarInit, ImmisAuxVarCopy

contains

! ************************************************************************** !

function ImmisAuxCreate()
  ! 
  ! ImmisAuxVarCreate: Allocate and initialize auxiliary object
  ! 
  ! Author: Chuan Lu
  ! Date: 02/27/08
  ! 

  use Option_module

  implicit none
  
  type(Immis_type), pointer :: ImmisAuxCreate
  
  type(Immis_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  aux%n_zero_rows = 0
  allocate(aux%immis_parameter)
  nullify(aux%immis_parameter%sir)
  nullify(aux%immis_parameter%ckwet)
  nullify(aux%immis_parameter%dencpr)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
  nullify(aux%res_old_FL)
  nullify(aux%res_old_AR)
  nullify(aux%delx)

  ImmisAuxCreate => aux
  
end function ImmisAuxCreate

! ************************************************************************** !

subroutine ImmisAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Chuan Lu
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(Immis_auxvar_type) :: auxvar
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
!   allocate ( auxvar%auxvar_elem(nvar)%xmol(option%nphase*option%nflowspec))
!   allocate ( auxvar%auxvar_elem(nvar)%diff(option%nphase*option%nflowspec))
    if (nvar>0) &
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
!   auxvar%auxvar_elem(nvar)%xmol = 0.d0
!   auxvar%auxvar_elem(nvar)%diff = 0.d0
#if 0
     auxvar%auxvar_elem(nvar)%dsat_dp = 0.d0
     auxvar%auxvar_elem(nvar)%dden_dp = 0.d0
     auxvar%auxvar_elem(nvar)%dkvr_dp = 0.d0
#endif
  enddo

end subroutine ImmisAuxVarInit

! ************************************************************************** !

subroutine ImmisAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Chuan Lu
  ! Date: 10/13/0
  ! 

  use Option_module

  implicit none
  
  type(Immis_auxvar_elem_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%avgmw = auxvar%avgmw
  auxvar2%vis = auxvar%vis
  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
!  auxvar2%kr = auxvar%kr
!  auxvar2%dkr_dp = auxvar%dkr_dp
!  auxvar2%vis = auxvar%vis
!  auxvar2%dvis_dp = auxvar%dvis_dp
  auxvar2%kvr = auxvar%kvr
#if 0
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dt = auxvar%dden_dt
  auxvar2%dkvr_dp = auxvar%dkvr_dp
  auxvar2%dkvr_dt = auxvar%dkvr_dt
  auxvar2%dh_dp = auxvar%dh_dp
  auxvar2%dh_dt = auxvar%dh_dt
  auxvar2%du_dp = auxvar%du_dp
  auxvar2%du_dt = auxvar%du_dt  
#endif
! auxvar2%xmol = auxvar%xmol
! auxvar2%diff = auxvar%diff

end subroutine ImmisAuxVarCopy

! ************************************************************************** !

subroutine ImmisAuxVarCompute_NINC(x,auxvar,saturation_function, &
                                   fluid_properties,option)
  ! 
  ! ImmisAuxVarCompute_NI: Computes auxiliary variables for each grid cell
  ! No increments
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Option_module
  use EOS_Water_module
  use Gas_EOS_module
  use co2eos_module
  use co2_span_wagner_module
  use co2_span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use Saturation_Function_module
  use Fluid_module
  use Mphase_pckr_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(Immis_auxvar_elem_type) :: auxvar
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: p,t,temp,p2,err
  PetscReal :: henry
  PetscReal :: dg,dddp,dddt
  PetscReal :: fg,dfgdp,dfgdt,xphi
  PetscReal :: eng,hg,dhdp,dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase),u(option%nphase),kr(option%nphase)
  PetscReal :: xm_nacl, x_nacl, vphi             
  PetscReal :: aux(1)
  PetscInt :: iflag
   
  
  
  auxvar%sat = 0.d0
  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%den = 0.d0
  auxvar%avgmw = 0.d0
! auxvar%vis = 0.d0
  auxvar%pc = 0.d0
  auxvar%kvr = 0.d0
! auxvar%xmol = 0.d0
! auxvar%diff = 0.d0
  kr = 0.d0
 
  auxvar%pres = x(1)  
  auxvar%temp = x(2)
  p = auxvar%pres
  t = auxvar%temp

  if (x(3)<0.D0)x(3) = 0.D0
  if (x(3)>1.D0)x(3) = 1.D0
  
  auxvar%sat(2) = x(3)
  if (auxvar%sat(2) < 0.D0) then
!   print *,'tran:',iphase, x(1:3)
    auxvar%sat(2) = 0.D0
  endif
! if (auxvar%sat(2) > 1.D0) print *,'tran:',iphase, x(1:3)
  auxvar%sat(1) = 1.D0 - auxvar%sat(2)
  auxvar%pc(:) = 0.D0
  temp = 1.D-2


! ********************* Gas phase properties ***********************
    call EOSWaterSaturationPressure(t, sat_pressure, ierr)
    err = 1.D0
    p2 = p

    if (p2 >= 5.d4) then
      if (option%co2eos == EOS_SPAN_WAGNER) then
! ************ Span-Wagner EOS ********************             
        select case(option%itable)  
          case(0,1,2,4,5)
            if (option%itable >=4) then
                ! print *,' interp', itable
              call co2_sw_interp(p2*1.D-6, t,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
            else
              iflag = 1
              call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag, &
                     option%itable)
            endif
            dg = dg / FMWCO2
            fg = fg*1.D6 
            hg = hg*FMWCO2
            xphi = fg/p2
            
! ************* Span-Wagner EOS with Bi-Cubic Spline interpolation ********
          case(3) 
            call sw_prop(t,p2*1D-6,dg,hg,eng,fg)
            call visco2(t,dg,visg)
            dg = dg / FMWCO2
            fg = fg*1.D6 
            hg = hg*FMWCO2
            xphi = fg/p2
        end select
          
      elseif (option%co2eos == EOS_MRK) then
! MRK eos [modified version from Kerrick and Jacobs (1981) and Weir et al. (1996).]     
        call CO2(t,p2,dg,fg,xphi,hg)
        call visco2(t,dg,visg)
        dg = dg / FMWCO2
        hg = hg*FMWCO2*option%scale
          !      print *, 'translator', p2, t, dg,hg,visg
      else
        call printErrMsg(option,'pflow Immis ERROR: Need specify CO2 EOS')
      endif
    else      
      call ideal_gaseos_noderiv(p2,t,dg,hg,eng)
      ! J/kmol -> whatever
      hg = hg * option%scale
      eng = eng * option%scale
      call visco2(t,dg*FMWCO2,visg)
      fg=p2
      xphi = 1.D0
    endif

!   call Henry_duan_sun(t,p2*1D-5,henry,lngamco2, &
!     option%m_nacl,option%m_nacl)
!   henry= 1D0 / (FMWH2O*1.D-3) / (henry*1.D-5)/xphi 
   
   pw = p
   auxvar%den(2) = dg
   auxvar%h(2) = hg  
   auxvar%u(2) = hg - p/dg * option%scale
   auxvar%pc(2) = 0.D0
   
!   auxvar%diff(option%nflowspec+1:option%nflowspec*2) = 2.13D-5
!       fluid_properties%diff_base(2)
! Note: not temperature dependent yet.       
   auxvar%zco2 = auxvar%den(2)/(p/IDEAL_GAS_CONSTANT/(t+273.15D0)*1D-3)
!***************  Liquid phase properties **************************
 
!  avgmw(1)= xmol(1)*FMWH2O + xmol(2)*FMWCO2 
  call EOSWaterDensity(t,pw,dw_kg,dw_mol,ierr) 
  call EOSWaterEnthalpy(t,pw,hw,ierr) 
  ! J/kmol -> whatever units
  hw = hw * option%scale

  auxvar%h(1) = hw
  auxvar%u(1) = auxvar%h(1) - pw /dw_mol*option%scale
!  auxvar%diff(1:option%nflowspec) = 1D-9
  ! fluid_properties%diff_base(1)

  xm_nacl = option%m_nacl*FMWNACL
  xm_nacl = xm_nacl /(1.D3 + xm_nacl)
  aux(1) = xm_nacl
  call EOSWaterDensityExt(t,p,aux,dw_kg,dw_mol,ierr)
  call EOSWaterViscosityNaCl(t,p,xm_nacl,visl)

!FEHM mixing ****************************
! den(1) = xmol(2)*dg + xmol(1)*dw_mol
! ideal mixing    
! den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

! Garcia mixing **************************
  x_nacl =  option%m_nacl/(option%m_nacl + 1D3/FMWH2O)
! **  xmol(1) = xh2o + xnacl

  auxvar%avgmw(1) = (1.D0 - x_nacl)*FMWH2O + x_nacl*FMWNACL
  auxvar%avgmw(2) = FMWCO2
  auxvar%den(1) = dw_kg/auxvar%avgmw(1)

 ! Hebach, J. Chem.Eng.Data 2004 (49),p950 ***********
 !  den(1) = 949.7109D0 + p * (0.559684D-6 - 0.00097D-12 * p) &  
 !      + (t+273.15)*(0.883148 - 0.00228*(t+273.15))  
 !  den(1) = dw_kg + (den(1)-dw_kg)*xmol(2)/p*henry
 !  den(1) = den(1)/avgmw(1)
!****************************** 2 phase S-Pc-kr relation *********************************
    auxvar%pc = 0.D0

    if (saturation_function%hysteresis_id <= 0.1D0) then 
      call pckrNH_noderiv(auxvar%sat,auxvar%pc,kr, &
                                   saturation_function, &
                                   option)
      pw=p !-pc(1)
     
    else
      call pckrHY_noderiv(auxvar%sat,auxvar%hysdat,auxvar%pc,kr, &
                                   saturation_function, &
                                   option)
    end if

!   call SaturationFunctionCompute(auxvar%pres,auxvar%sat,kr, &
!                                   ds_dp,dkr_dp, &
!                                   saturation_function, &
!                                   por,perm, &
!                                   option)
    auxvar%kvr(2) = kr(2)/visg     
    auxvar%kvr(1) = kr(1)/visl
    auxvar%vis(2) = visg     
    auxvar%vis(1) = visl
  
   
!  print *,'immis_aux: ',auxvar%den,auxvar%avgmw,auxvar%vis,auxvar%kvr

end subroutine ImmisAuxVarCompute_NINC

! ************************************************************************** !

subroutine ImmisAuxVarCompute_WINC(x, delx, auxvar,saturation_function, &
                                    fluid_properties,option)

  use Option_module
  
  use Saturation_Function_module
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof), xx(option%nflowdof), delx(option%nflowdof)
  type(Immis_auxvar_elem_type) :: auxvar(1:option%nflowdof)
 ! PetscInt :: iphase

  PetscInt :: n 
  
  do n=1, option%nflowdof
     xx=x;  xx(n)=x(n)+ delx(n)
! ***   note: var_node here starts from 1 to option%flowdof ***
    call  ImmisAuxVarCompute_NINC(xx,auxvar(n),saturation_function, &
                                   fluid_properties, option)
  enddo

end subroutine ImmisAuxVarCompute_WINC

! ************************************************************************** !

subroutine ImmisAuxVarDestroy(auxvar)
  ! 
  ! AuxVarDestroy: Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(Immis_auxvar_elem_type) :: auxvar
  
! if (associated(auxvar%xmol)) deallocate(auxvar%xmol)
! nullify(auxvar%xmol)
! if (associated(auxvar%diff))deallocate(auxvar%diff)
! nullify(auxvar%diff)
  if (associated(auxvar%pc))deallocate(auxvar%pc)
  nullify(auxvar%pc)
  if (associated(auxvar%sat))deallocate(auxvar%sat)
  nullify(auxvar%sat)
  if (associated(auxvar%u))deallocate(auxvar%u)
  nullify(auxvar%u)
  if (associated(auxvar%h))deallocate(auxvar%h)
  nullify(auxvar%h)
  if (associated(auxvar%den))deallocate(auxvar%den)
  nullify(auxvar%den)
  if (associated(auxvar%avgmw))deallocate(auxvar%avgmw)
  nullify(auxvar%avgmw)
  if (associated(auxvar%vis))deallocate(auxvar%vis)
  nullify(auxvar%vis)
end subroutine ImmisAuxVarDestroy

! ************************************************************************** !

subroutine ImmisAuxDestroy(aux, option)
  ! 
  ! RichardsAuxDestroy: Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module
  implicit none

  type(Immis_type), pointer :: aux
  type(option_type) :: option
  PetscInt :: iaux, ielem
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    do ielem= 0, option%nflowdof 
      call ImmisAuxVarDestroy(aux%auxvars(iaux)%auxvar_elem(ielem))
    enddo
    deallocate(aux%auxvars)
  enddo
  nullify(aux%auxvars)
  
  do iaux = 1, aux%num_aux_bc
    do ielem= 0, option%nflowdof 
      call ImmisAuxVarDestroy(aux%auxvars_bc(iaux)%auxvar_elem(ielem))
    enddo
    deallocate(aux%auxvars_bc)
  enddo
  nullify(aux%auxvars_bc)
  
  do iaux = 1, aux%num_aux_ss
    do ielem = 0, option%nflowdof 
      call ImmisAuxVarDestroy(aux%auxvars_ss(iaux)%auxvar_elem(ielem))
    enddo
    deallocate(aux%auxvars_ss)
  enddo
  nullify(aux%auxvars_ss)
  
  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) deallocate(aux%auxvars_bc)
  nullify(aux%auxvars_bc)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%res_old_AR)) deallocate(aux%res_old_AR)
  nullify(aux%res_old_AR)
  if (associated(aux%res_old_FL)) deallocate(aux%res_old_FL)
  nullify(aux%res_old_FL)
  if (associated(aux%delx)) deallocate(aux%delx)
  nullify(aux%delx)
  if (associated(aux%immis_parameter)) then
    if (associated(aux%immis_parameter%dencpr)) deallocate(aux%immis_parameter%dencpr)
    nullify(aux%immis_parameter%dencpr)
    if (associated(aux%immis_parameter%ckwet)) deallocate(aux%immis_parameter%ckwet)
    nullify(aux%immis_parameter%ckwet)
    if (associated(aux%immis_parameter%ckdry)) deallocate(aux%immis_parameter%ckdry)
    nullify(aux%immis_parameter%ckdry)
    if (associated(aux%immis_parameter%sir)) deallocate(aux%immis_parameter%sir)
    nullify(aux%immis_parameter%sir)
    deallocate(aux%immis_parameter)
  endif
  nullify(aux%immis_parameter%dencpr)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine ImmisAuxDestroy



end module Immis_Aux_module


