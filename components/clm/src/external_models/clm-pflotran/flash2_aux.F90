module Flash2_Aux_module
  use Mphase_pckr_module
  use PFLOTRAN_Constants_module

  implicit none
  
  private 
!#define GARCIA 1
#define DUANDEN 1

#include "petsc/finclude/petscsys.h"

  type, public :: Flash2_auxvar_elem_type
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
!    PetscReal :: dvis_dp
!    PetscReal :: kr
!    PetscReal :: dkr_dp
 end type Flash2_auxvar_elem_type

  type, public :: Flash2_auxvar_type
    
  type(Flash2_auxvar_elem_type), pointer :: auxvar_elem(:)
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
  end type Flash2_auxvar_type
  
  type, public :: Flash2_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckwet(:)
    PetscReal, pointer :: ckdry(:)
    PetscReal, pointer :: sir(:,:)
  end type Flash2_parameter_type
    
  type, public :: Flash2_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(Flash2_parameter_type), pointer :: Flash2_parameter
    type(Flash2_auxvar_type), pointer :: auxvars(:)
    type(Flash2_auxvar_type), pointer :: auxvars_bc(:)
    type(Flash2_auxvar_type), pointer :: auxvars_ss(:)
    PetscReal , pointer :: Resold_AR(:,:)
    PetscReal , pointer :: Resold_BC(:,:)
    PetscReal , pointer :: Resold_FL(:,:)
    PetscReal , pointer :: delx(:,:)
  end type Flash2_type

  

  public :: Flash2AuxCreate, Flash2AuxDestroy, &
            Flash2AuxVarCompute_NINC, Flash2AuxVarCompute_WINC,&
            Flash2AuxVarInit, Flash2AuxVarCopy

contains

! ************************************************************************** !

function Flash2AuxCreate()
  ! 
  ! Flash2AuxVarCreate: Allocate and initialize auxiliary object
  ! 
  ! Author: Chuan Lu
  ! Date: 02/27/08
  ! 

  use Option_module

  implicit none
  
  type(Flash2_type), pointer :: Flash2AuxCreate
  
  type(Flash2_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_zero_rows = 0
  allocate(aux%Flash2_parameter)
  nullify(aux%Flash2_parameter%sir)
  nullify(aux%Flash2_parameter%ckwet)
  nullify(aux%Flash2_parameter%dencpr)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  Flash2AuxCreate => aux
  
end function Flash2AuxCreate

! ************************************************************************** !

subroutine Flash2AuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Chuan Lu
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(Flash2_auxvar_type) :: auxvar
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
    if (nvar > 0) &
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
#if 0
    auxvar%auxvar_elem(nvar)%dsat_dp = 0.d0
    auxvar%auxvar_elem(nvar)%dden_dp = 0.d0
    auxvar%auxvar_elem(nvar)%dkvr_dp = 0.d0
#endif
  enddo

end subroutine Flash2AuxVarInit

! ************************************************************************** !

subroutine Flash2AuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Chuan Lu
  ! Date: 10/13/0
  ! 

  use Option_module

  implicit none
  
  type(Flash2_auxvar_elem_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%avgmw = auxvar%avgmw
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
!  auxvar2%xmol = auxvar%xmol
!  auxvar2%diff = auxvar%diff

end subroutine Flash2AuxVarCopy

! ************************************************************************** !

subroutine Flash2AuxVarCompute_NINC(x,auxvar,global_auxvar, &
             saturation_function,fluid_properties,option,xphico2)
  ! 
  ! Flash2AuxVarCompute_NI: Computes auxiliary variables for each grid cell
  ! No increments
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Option_module
  use Global_Aux_module  
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
  type(Flash2_auxvar_elem_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  PetscReal, optional :: xphico2

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal ::  p, t, temp, p2, err
  PetscReal :: henry,lngamco2
  PetscReal :: dg, dddp, dddt
  PetscReal :: fg, dfgdp, dfgdt, xphi
  PetscReal :: eng,hg, dhdp, dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase), u(option%nphase), kr(option%nphase)
  PetscReal :: m_na,m_cl,m_nacl, xm_nacl, x_nacl, y_nacl, vphi             
  PetscReal :: tk, xco2, pw_kg, x1, vphi_a1, vphi_a2 
  PetscReal :: Qkco2, mco2,xco2eq
  PetscReal :: tmp 
  PetscReal :: aux(1)
  PetscInt :: iflag  
  
  auxvar%sat = 0.d0
  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%den = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%pc = 0.d0
  auxvar%kvr = 0.d0
! auxvar%xmol = 0.d0
! auxvar%diff = 0.d0
  kr = 0.d0
 
  auxvar%pres = x(1)  
  auxvar%temp = x(2)

  p = auxvar%pres
  t = auxvar%temp

! ********************* Gas phase properties ***********************
    call EOSWaterSaturationPressure(t, sat_pressure, ierr)
    err=1.D0
    p2 = p

    if (p2 >= 5.d4) then
      if (option%co2eos == EOS_SPAN_WAGNER) then
! ************ Span-Wagner EOS ********************             
        select case(option%itable)  
          case(0,1,2,4,5)
            if (option%itable >= 4) then
                ! print *,' interp', itable
              call co2_sw_interp(p2*1.D-6, t,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
            else
              iflag = 1
              call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg,&
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag, &
                     option%itable)
            endif
            dg= dg / FMWCO2
            fg= fg * 1.D6 
            hg= hg * FMWCO2
            xphi = fg/p2
! ************* Span-Wagner EOS with Bi-Cubic Spline interpolation ********
          case(3) 
            call sw_prop(t,p2*1D-6,dg,hg, eng, fg)
            call visco2(t, dg, visg)
            dg= dg / FMWCO2
            fg= fg * 1.D6 
            hg= hg * FMWCO2
            xphi = fg/p2
        end select
      elseif (option%co2eos == EOS_MRK)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]     
        call CO2(t, p2,  dg,fg, xphi, hg)
        call visco2( t,dg,visg)
        dg = dg / FMWCO2
        hg = hg * FMWCO2 *option%scale
          !      print *, 'translator', p2, t, dg,hg,visg
      else
        call printErrMsg(option,'pflow Flash2 ERROR: Need specify CO2 EOS')
      endif
    else      
      call ideal_gaseos_noderiv(p2, t,dg,hg,eng)
      ! J/kmol -> whatever
      hg = hg * option%scale
      eng = eng * option%scale        
      call visco2(t,dg*FMWCO2,visg)
      fg=p2
      xphi = 1.D0
    endif
 
!*********** Get Salinity properties ***********************
    m_na=option%m_nacl; m_cl=m_na; m_nacl=m_na 
    if (option%ntrandof>0) then
      m_na = global_auxvar%m_nacl(1)
      m_cl = global_auxvar%m_nacl(2)
      m_nacl = m_na
      if (m_cl> m_na) m_nacl = m_cl
    endif  

!** Calculate solubility of CO2 in aqueous phase *************
    call Henry_duan_sun(t,p2*1.D-5,henry,lngamco2,m_na,m_cl)

    Qkco2 = henry*xphi  ! convert from bar to Pa
    henry = 1.D0 / (FMWH2O*1.D-3) / (henry*1.D-5) / xphi 
    if (present(xphico2)) xphico2 = xphi
   
    mco2 = (p - sat_pressure)*1D-5 * Qkco2
    xco2eq = mco2/(1D3/fmwh2o + mco2 + m_nacl) 
   
    tmp= Henry/p
    if (x(3) < xco2eq) then
      ! water only
      auxvar%xmol(2) = x(3)
      auxvar%xmol(1) = 1.D0 - auxvar%xmol(2)
      auxvar%xmol(4) = auxvar%xmol(2)*tmp
      auxvar%xmol(3) = 1.D0 - auxvar%xmol(4)
      auxvar%sat(1) = 1.D0
      auxvar%sat(2) = 0.D0
      iphase = 1
    elseif (x(3) > (1.D0-sat_pressure/p)) then
	    !gas only
      iphase = 2
      auxvar%xmol(4) = x(3)
      auxvar%xmol(3) = 1.D0 - auxvar%xmol(4)
      auxvar%xmol(2) = auxvar%xmol(4)/tmp
      auxvar%xmol(1) = 1.D0 - auxvar%xmol(2)
      auxvar%sat(1) = 0.D0 !1.D-8
      auxvar%sat(2) = 1.D0
    else 
      iphase = 3
      auxvar%xmol(1) = 1.D0 - xco2eq
      auxvar%xmol(2) = xco2eq
      auxvar%xmol(3) = sat_pressure/p*auxvar%xmol(1)
      auxvar%xmol(4) = 1.D0 - auxvar%xmol(3)
    endif 

! **************  Gas phase properties ********************
    auxvar%avgmw(2) = auxvar%xmol(3)*FMWH2O + auxvar%xmol(4)*FMWCO2
    pw = p
    call EOSWaterDensity(t,pw,dw_kg,dw_mol,ierr)
    call EOSWaterEnthalpy(t,pw,hw,ierr)
    hw = hw * option%scale ! J/kmol -> whatever units
    auxvar%den(2) = 1.D0/(auxvar%xmol(4)/dg + auxvar%xmol(3)/dw_mol)
    auxvar%h(2) = hg  
    auxvar%u(2) = hg - p/dg * option%scale
    
!   auxvar%diff(option%nflowspec+1:option%nflowspec*2) = 2.13D-5
    auxvar%diff(option%nflowspec+1:option%nflowspec*2) = &
      fluid_properties%gas_diffusion_coefficient &
      * 101325.d0 / p * ((t+273.15)/273.15d0)**1.8d0

!       fluid_properties%diff_base(2)

!  z factor    
    auxvar%zco2=auxvar%den(2)/(p/IDEAL_GAS_CONSTANT/(t+273.15D0)*1D-3)

 !***************  Liquid phase properties **************************
 
!    avgmw(1)= xmol(1)* FMWH2O + xmol(2) * FMWCO2 
    auxvar%h(1) = hw
    auxvar%u(1) = auxvar%h(1) - pw/dw_mol*option%scale
    
    auxvar%diff(1:option%nflowspec) = fluid_properties%diffusion_coefficient
  ! fluid_properties%diff_base(1) need more work here. Add temp. dependence.
  
    xm_nacl = m_nacl * FMWNACL
    xm_nacl = xm_nacl /(1.D3 + xm_nacl)
    aux(1) = xm_nacl
    call EOSWaterDensityExt(t,p,aux,dw_kg,dw_mol,ierr)
    call EOSWaterViscosityNaCl(t,p,xm_nacl,visl)
    
    y_nacl = m_nacl/( m_nacl + 1D3/FMWH2O)
!   y_nacl is the mole fraction
    auxvar%avgmw(1) = auxvar%xmol(1)*((1D0 - y_nacl) * FMWH2O &
       + y_nacl * FMWNACL) + auxvar%xmol(2) * FMWCO2

!duan mixing **************************
#ifdef DUANDEN
  call EOSWaterDuanMixture (t,p,auxvar%xmol(2),y_nacl, &
    auxvar%avgmw(1),dw_kg,auxvar%den(1))
#endif 

! Garcia mixing **************************
#ifdef GARCIA
  vphi = 1D-6*(37.51D0 + t&
       *(-9.585D-2 + t*(8.74D-4 - t*5.044D-7)))
  auxvar%den(1) = dw_kg/(1D0-(FMWCO2*1D-3-dw_kg*vphi)&
       *auxvar%xmol(2)/(auxvar%avgmw(1)*1D-3))
  auxvar%den(1) = auxvar%den(1)/auxvar%avgmw(1)
#endif  
       

!FEHM mixing ****************************
!  den(1) = xmol(2)*dg + xmol(1)*dw_mol
! ideal mixing    
  !den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

! Garcia mixing **************************
        
 ! Hebach, J. Chem.Eng.Data 2004 (49),p950 ***********
 !   den(1)= 949.7109D0 + p * (0.559684D-6 - 0.00097D-12 * p) &  
 !      + (t+273.15)*(0.883148 - 0.00228*(t+273.15))  
 !  den(1)=dw_kg + (den(1)-dw_kg)*xmol(2)/p*henry
 !  den(1)=den(1)/avgmw(1)

!****** calcultate phase splition for 2 phase coexist condition *******
  select case(iphase)
  case(1,2)
  case(3)
    auxvar%sat(2) = auxvar%den(1)* ( x(3) - auxvar%xmol(2))/&
      (auxvar%den(2) * (auxvar%xmol(4)-x(3)) - auxvar%den(1)*(auxvar%xmol(2)-x(3)))
    if (auxvar%sat(2) >1D0 .or. auxvar%sat(2) <0D0) print *,'z->s error: ',auxvar%sat(2)
    if (auxvar%sat(2) > 1D0) auxvar%sat(2) = 1D0
    if (auxvar%sat(2) < 0D0) auxvar%sat(2) = 0D0  
    auxvar%sat(1) = 1D0 - auxvar%sat(2)
  end select
 
!******************************** 2 phase S-Pc-kr relation ********************
    auxvar%pc =0.D0

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

end subroutine Flash2AuxVarCompute_NINC

! ************************************************************************** !

subroutine Flash2AuxVarCompute_WINC(x, delx, auxvar,global_auxvar,saturation_function, &
                                    fluid_properties,option)

  use Option_module
  use Global_Aux_module
  
  use Saturation_Function_module
  use Fluid_module
  
  implicit none

  type(option_type) :: option
  type(fluid_property_type) :: fluid_properties
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof), xx(option%nflowdof), delx(option%nflowdof)
  type(Flash2_auxvar_elem_type) :: auxvar(1:option%nflowdof)
  type(global_auxvar_type) :: global_auxvar
 ! PetscInt :: iphase

  PetscInt :: n 
  
  do n=1, option%nflowdof
     xx=x;  xx(n)=x(n)+ delx(n)
! ***   note: var_node here starts from 1 to option%flowdof ***
    call  Flash2AuxVarCompute_NINC(xx,auxvar(n),global_auxvar, &
      saturation_function,fluid_properties, option)
  enddo

end subroutine Flash2AuxVarCompute_WINC

! ************************************************************************** !

subroutine Flash2AuxVarDestroy(auxvar)
  ! 
  ! AuxVarDestroy: Deallocates a FLASH2 auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(Flash2_auxvar_elem_type) :: auxvar
  
!  if (associated(auxvar%xmol)) deallocate(auxvar%xmol)
!  nullify(auxvar%xmol)
!  if (associated(auxvar%diff))deallocate(auxvar%diff)
!  nullify(auxvar%diff)
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
  if (associated(auxvar%den))deallocate(auxvar%vis)
  nullify(auxvar%vis)
  if (associated(auxvar%avgmw))deallocate(auxvar%avgmw)
  nullify(auxvar%avgmw)
end subroutine Flash2AuxVarDestroy

! ************************************************************************** !

subroutine Flash2AuxDestroy(aux, option)
  ! 
  ! RichardsAuxDestroy: Deallocates a FLASH2 auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module
  implicit none

  type(Flash2_type), pointer :: aux
  type(option_type) :: option
  PetscInt :: iaux, ielem
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    do ielem= 0, option%nflowdof 
      call Flash2AuxVarDestroy(aux%auxvars(iaux)%auxvar_elem(ielem))
    enddo
  enddo

  do iaux = 1, aux%num_aux_bc
    do ielem= 0, option%nflowdof 
      call Flash2AuxVarDestroy(aux%auxvars_bc(iaux)%auxvar_elem(ielem))
    enddo
  enddo  

  do iaux = 1, aux%num_aux_ss
    do ielem= 0, option%nflowdof 
      call Flash2AuxVarDestroy(aux%auxvars_ss(iaux)%auxvar_elem(ielem))
    enddo
  enddo  

  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) deallocate(aux%auxvars_bc)
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) deallocate(aux%auxvars_ss)
  nullify(aux%auxvars_ss)

  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

  if (associated(aux%Flash2_parameter)) then
    if (associated(aux%Flash2_parameter%dencpr)) deallocate(aux%Flash2_parameter%dencpr)
    nullify(aux%Flash2_parameter%dencpr)
    if (associated(aux%Flash2_parameter%ckwet)) deallocate(aux%Flash2_parameter%ckwet)
    nullify(aux%Flash2_parameter%ckwet)
    if (associated(aux%Flash2_parameter%ckdry)) deallocate(aux%Flash2_parameter%ckdry)
    nullify(aux%Flash2_parameter%ckdry)
    if (associated(aux%Flash2_parameter%sir)) deallocate(aux%Flash2_parameter%sir)
    nullify(aux%Flash2_parameter%sir)
    deallocate(aux%Flash2_parameter)
  endif
  nullify(aux%Flash2_parameter%dencpr)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine Flash2AuxDestroy



end module Flash2_Aux_module


