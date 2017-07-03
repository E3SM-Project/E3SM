module Mphase_Aux_module
  
  use Mphase_pckr_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private

!#define GARCIA 1
#define DUANDEN 1

#include "petsc/finclude/petscsys.h"

  type, public :: mphase_auxvar_elem_type
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
  end type mphase_auxvar_elem_type

  type, public :: mphase_auxvar_type
    
    type(mphase_auxvar_elem_type), pointer :: auxvar_elem(:) 
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
  end type mphase_auxvar_type
  
  type, public :: mphase_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckwet(:)
    PetscReal, pointer :: sir(:,:)
  end type mphase_parameter_type
  
  type, public :: mphase_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss

    PetscReal, pointer :: res_old_AR(:,:)
    PetscReal, pointer :: res_old_FL(:,:)
    PetscReal, pointer :: delx(:,:)
  
    type(mphase_parameter_type), pointer :: mphase_parameter
    type(mphase_auxvar_type), pointer :: auxvars(:)
    type(mphase_auxvar_type), pointer :: auxvars_bc(:)
    type(mphase_auxvar_type), pointer :: auxvars_ss(:)
  end type mphase_type


  public :: MphaseAuxCreate, MphaseAuxDestroy, &
            MphaseAuxVarCompute_NINC, MphaseAuxVarCompute_WINC, &
            MphaseAuxVarInit, MphaseAuxVarCopy

contains

! ************************************************************************** !

function MphaseAuxCreate()
  ! 
  ! MphaseAuxVarCreate: Allocate and initialize auxiliary object
  ! Author: Chuan Lu
  ! 

  use Option_module

  implicit none
  
  type(mphase_type), pointer :: MphaseAuxCreate
  
  type(mphase_type), pointer :: aux

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
  allocate(aux%mphase_parameter)
  nullify(aux%mphase_parameter%sir)
  nullify(aux%mphase_parameter%ckwet)
  nullify(aux%mphase_parameter%dencpr)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
  nullify(aux%res_old_AR)
  nullify(aux%res_old_FL)
  nullify(aux%delx)
  
  MphaseAuxCreate => aux
  
end function MphaseAuxCreate

! ************************************************************************** !

subroutine MphaseAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! Author: Chuan Lu
  ! 

  use Option_module

  implicit none
  
  type(mphase_auxvar_type) :: auxvar
  type(option_type) :: option

  PetscInt :: var_elem_size, var_node_size
  PetscInt :: nvar 

  allocate(auxvar%auxvar_elem(0 : option%nflowdof))
  allocate(auxvar%auxvar_elem(0)%hysdat(4))
 
  do nvar = 0, option%nflowdof
     allocate ( auxvar%auxvar_elem(nvar)%sat(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%den(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%avgmw(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%h(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%u(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%pc(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%kvr(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%vis(option%nphase))
     allocate ( auxvar%auxvar_elem(nvar)%xmol(option%nphase*option%nflowspec))
     allocate ( auxvar%auxvar_elem(nvar)%diff(option%nphase*option%nflowspec))
     if (nvar>0) &
       auxvar%auxvar_elem(nvar)%hysdat => auxvar%auxvar_elem(0)%hysdat

     auxvar%auxvar_elem(nvar)%pres = 0.d0
     auxvar%auxvar_elem(nvar)%temp = 0.d0
     auxvar%auxvar_elem(nvar)%sat = 0.d0
     auxvar%auxvar_elem(nvar)%den = 0.d0
     auxvar%auxvar_elem(nvar)%avgmw = 0.d0
     auxvar%auxvar_elem(nvar)%h = 0.d0
     auxvar%auxvar_elem(nvar)%u = 0.d0
     auxvar%auxvar_elem(nvar)%pc = 0.d0
     auxvar%auxvar_elem(nvar)%kvr = 0.d0
     auxvar%auxvar_elem(nvar)%xmol = 0.d0
     auxvar%auxvar_elem(nvar)%diff = 0.d0
     auxvar%auxvar_elem(nvar)%vis = 0.d0
#if 0
     auxvar%auxvar_elem(nvar)%dsat_dp = 0.d0
     auxvar%auxvar_elem(nvar)%dden_dp = 0.d0
     auxvar%auxvar_elem(nvar)%dkvr_dp = 0.d0
#endif
  enddo

end subroutine MphaseAuxVarInit

! ************************************************************************** !

subroutine MphaseAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  use Option_module

  implicit none
  
  type(mphase_auxvar_elem_type) :: auxvar, auxvar2
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
  auxvar2%vis = auxvar%vis
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
  auxvar2%xmol = auxvar%xmol
  auxvar2%diff = auxvar%diff

end subroutine MphaseAuxVarCopy

! ************************************************************************** !

subroutine MphaseAuxVarCompute_NINC(x,auxvar,global_auxvar,iphase,saturation_function, &
                                   fluid_properties,option,xphico2)
  ! 
  ! MphaseAuxVarCompute_NI: Computes auxiliary variables for each grid cell
  ! No increments
  ! Author: Chuan Lu
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
  type(mphase_auxvar_elem_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  PetscReal, optional :: xphico2

  PetscErrorCode :: ierr
  PetscReal :: pw, dw_kg, dw_mol, hw, sat_pressure, visl
  PetscReal :: p, t, temp, p2, err
  PetscReal :: henry, lngamco2
  PetscReal :: dg, dddp, dddt, m_na, m_cl, m_nacl
  PetscReal :: fg, dfgdp, dfgdt, xphi
  PetscReal :: eng, hg, dhdp, dhdt
  PetscReal :: visg, dvdp, dvdt
  PetscReal :: h(option%nphase), u(option%nphase), kr(option%nphase)
  PetscReal :: xm_nacl, y_nacl, vphi             
  PetscReal :: tk, xco2, pw_kg, x1, vphi_a1, vphi_a2 
  PetscReal :: Qkco2, mco2, xco2eq
  PetscReal :: aux(1)
  PetscInt :: iflag
  
  auxvar%den = 0.d0
  auxvar%sat = 0.d0
  
  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%xmol = 0.d0
  auxvar%pc = 0.d0
  auxvar%kvr = 0.d0
  auxvar%diff = 0.d0
  kr = 0.d0
 
  auxvar%pres = x(1)  
  auxvar%temp = x(2)

  p = auxvar%pres
  t = auxvar%temp
  
  select case(iphase)
    case(1)
!******* aqueous phase exists ***********
      auxvar%xmol(2) = x(3)
!      if (auxvar%xmol(2) < 0.D0) print *,'tran:',iphase, x(1:3)
!      if (auxvar%xmol(2) > 1.D0) print *,'tran:',iphase, x(1:3)
!pcl  if (x(3) < 0.d0) then
!pcl    option%io_buffer = 'CO2 mole fraction below zero.  It is likely ' // &
!pcl      'that CO2 aqueous concentrations in transport are inconsistent with flow.'
!pcl    call printErrMsgByRank(option)
!pcl  endif
      auxvar%xmol(1) = 1.D0 - auxvar%xmol(2)
      auxvar%pc(:) = 0.D0
      auxvar%sat(1) = 1.D0
      auxvar%sat(2) = 0.D0
      kr(1)= 1.D0
      kr(2)= 0.D0
    case(2)
!******* gas phase exists ***********
      auxvar%xmol(4) = x(3)
!      if (auxvar%xmol(4) < 0.D0) print *,'tran:',iphase, x(1:3)
!      if (auxvar%xmol(4) > 1.D0) print *,'tran:',iphase, x(1:3)
      auxvar%xmol(3) = 1.D0 - auxvar%xmol(4)
      auxvar%pc(:) = 0.D0
      auxvar%sat(1) = 0.D0
      auxvar%sat(2) = 1.D0
      kr(1)= 0.D0
      kr(2)= 1.D0
    case(3)    
!******* 2-phase phase exists ***********
      auxvar%sat(2) = x(3)
      if (auxvar%sat(2) < 0.D0)then
!        print *,'tran:',iphase, x(1:3)
        auxvar%sat(2) = 0.D0
      endif
!      if (auxvar%sat(2)> 1.D0) print *,'tran:',iphase, x(1:3)
      auxvar%sat(1) = 1.D0 - auxvar%sat(2)
      auxvar%pc(:) = 0.D0
      temp = 1.D-2
      auxvar%xmol(1)=1.D0; auxvar%xmol(2)=0.D0
      auxvar%xmol(3)=temp; auxvar%xmol(4)=1.D0-auxvar%xmol(3)
   end select
! ********************* Gas phase properties ***********************
    call EOSWaterSaturationPressure(t, sat_pressure, ierr)
    err = 1.D0
    p2 = p

    if (p2 >= 5.d4) then
       
      if (option%co2eos == EOS_SPAN_WAGNER) then
! ************ Span-Wagner EOS ********************             
        select case(option%itable)  
          case(0,1,2,4,5)
            if (option%itable >= 4) then
                ! print *,' interp', itable
              call co2_sw_interp(p2*1.D-6,t,dg,dddt,dddp,fg, &
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
            else
              iflag = 1
              call co2_span_wagner(p2*1.D-6,t+273.15D0,dg,dddt,dddp,fg, &
                     dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,iflag, &
                     option%itable)
            endif

            dg = dg/FMWCO2
            fg = fg*1.D6 
            hg = hg*FMWCO2
            xphi = fg/p2
            
! ************* Span-Wagner EOS with Bi-Cubic Spline interpolation ********
          case(3) 
            call sw_prop(t,p2*1.D-6,dg,hg,eng,fg)
            call visco2(t,dg,visg)
            dg = dg/FMWCO2
            fg = fg*1.D6 
            hg = hg*FMWCO2
            xphi = fg/p2
          end select
          
       elseif (option%co2eos == EOS_MRK) then
       
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]     
          call CO2(t,p2,dg,fg,xphi,hg)
          call visco2(t,dg,visg)
          dg = dg/FMWCO2
          hg = hg*FMWCO2*option%scale
          !      print *, 'translator', p2,t,dg,hg,visg
       else
         call printErrMsg(option,'pflow mphase ERROR: Need specify CO2 EOS')
      endif
    else      
      call ideal_gaseos_noderiv(p2,t,dg,hg,eng)
      ! J/kmol -> whatever
      hg = hg * option%scale
      eng = eng * option%scale      
      call visco2(t,dg*FMWCO2,visg)
      fg = p2
      xphi = 1.D0
    endif

    m_na=option%m_nacl; m_cl=m_na; m_nacl=m_na 
    if (option%ntrandof > 0) then
      m_na = global_auxvar%m_nacl(1)
      m_cl = global_auxvar%m_nacl(2)
      m_nacl = m_na
      if (m_cl > m_na) m_nacl = m_cl
    endif  


    call Henry_duan_sun(t,p2*1.D-5,henry,lngamco2,m_na,m_cl)
    Qkco2 = henry*xphi  ! convert from bar to Pa
    henry = 1.D0/(FMWH2O*1.D-3)/(henry*1.D-5)/xphi 
    if (present(xphico2)) xphico2 = xphi
   
    mco2 = (p - sat_pressure)*1.D-5*Qkco2
    xco2eq = mco2/(1.D3/fmwh2o + mco2 + m_nacl)

!   print *,'mphase_duan_den: ',xco2eq,mco2,m_nacl,qkco2,p,p2,sat_pressure,t

!   question here: m_nacl or m_na+m_cl ?
   
    select case(iphase)     
    case(1)
      auxvar%xmol(4) = auxvar%xmol(2)*henry/p   
      auxvar%xmol(3) = 1.D0-auxvar%xmol(4)
      if (auxvar%xmol(3) < 0.D0) auxvar%xmol(3) = 0.D0
!     if (xmol(3) < 0.D0) xmol(3) = 0.D0
    case(2)   
      auxvar%xmol(2) = p*auxvar%xmol(4)/henry
      auxvar%xmol(1) = 1.D0 - auxvar%xmol(2)
    case(3)
      temp= sat_pressure/p
      auxvar%xmol(2) = xco2eq
      auxvar%xmol(1) = 1.D0 - xco2eq
      auxvar%xmol(3) = temp
      auxvar%xmol(4) = 1.D0 - temp            
    end select
    auxvar%avgmw(2) = auxvar%xmol(3)*FMWH2O + auxvar%xmol(4)*FMWCO2
    pw = p
    call EOSWaterDensity(t,pw,dw_kg,dw_mol,ierr) 
    call EOSWaterEnthalpy(t,pw,hw,ierr) 
    hw = hw * option%scale ! J/kmol -> whatever units
    auxvar%den(2) = 1.D0/(auxvar%xmol(4)/dg + auxvar%xmol(3)/dw_mol)
    auxvar%h(2) = hg  
    auxvar%u(2) = hg - p/dg*option%scale
    auxvar%pc(2) = 0.D0

!   auxvar%diff(option%nflowspec+1:option%nflowspec*2) = 2.13D-5
    auxvar%diff(option%nflowspec+1:option%nflowspec*2) = &
      fluid_properties%gas_diffusion_coefficient &
      * 101325.d0/p * ((t+273.15d0)/273.15d0)**1.8d0
!       fluid_properties%diff_base(2)

!   print *,'gas diff: ',fluid_properties%gas_diffusion_coefficient,p,t
       
!  z factor    
    auxvar%zco2=auxvar%den(2)/(p/IDEAL_GAS_CONSTANT/(t+273.15D0)*1.D-3)

!***************  Liquid phase properties **************************
 
!    avgmw(1)= xmol(1)* FMWH2O + xmol(2) * FMWCO2 
    auxvar%h(1) = hw
    auxvar%u(1) = hw - pw/dw_mol*option%scale

    auxvar%diff(1:option%nflowspec) = fluid_properties%diffusion_coefficient
  ! fluid_properties%diff_base(1)

  
    xm_nacl = m_nacl*FMWNACL
    xm_nacl = xm_nacl/(1.D3 + xm_nacl)
    aux(1) = xm_nacl
    call EOSWaterDensityExt(t,p,aux,dw_kg,dw_mol,ierr)
!   call EOSWaterViscosityNaCl(t,p,xm_nacl,visl)
    call EOSWaterViscosity(t,pw,sat_pressure,0.d0,visl,dvdt,dvdp,ierr)

!FEHM mixing ****************************
!  den(1) = xmol(2)*dg + xmol(1)*dw_mol
!  ideal mixing    
!  den(1) = 1.D0/(xmol(2)/dg + xmol(1)/dw_mol) !*c+(1-c)* 

!  m_nacl=option%m_nacl
!  if (reaction%species_idx%na_ion_id /= 0 .and. &
!     reaction%species_idx%cl_ion_id /= 0) then
!    m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
!    m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
!    m_nacl = m_na
!    if (m_cl > m_nacl) m_nacl=m_cl
!  endif  

    y_nacl = m_nacl/(m_nacl + 1.D3/FMWH2O)
! **  xmol(1) = xh2o + xnacl
    auxvar%avgmw(1) = auxvar%xmol(1)*((1.D0 - y_nacl)*FMWH2O &
       + y_nacl*FMWNACL) + auxvar%xmol(2)*FMWCO2

!duan mixing **************************
#ifdef DUANDEN
!                 units: t [C], p [MPa], dw_kg [kg/m^3]
  call EOSWaterDuanMixture (t,p,auxvar%xmol(2),y_nacl,auxvar%avgmw(1),dw_kg,auxvar%den(1))

#else
  auxvar%den(1) = dw_mol
#endif 
!print *,'mixture den: ',t,p,auxvar%xmol(2),y_nacl,auxvar%avgmw(1),dw_kg,auxvar%den(1),dw_mol

! Garcia mixing **************************
#ifdef GARCIA
  vphi = 1D-6*(37.51D0 + t &
       *(-9.585D-2 + t*(8.74D-4 - t*5.044D-7)))
  auxvar%den(1) = dw_kg/(1D0-(FMWCO2*1.D-3-dw_kg*vphi) &
       *auxvar%xmol(2)/(auxvar%avgmw(1)*1.D-3))
  auxvar%den(1) = auxvar%den(1)/auxvar%avgmw(1)
#endif  
       
 ! Hebach, J. Chem.Eng.Data 2004 (49),p950 ***********
 !   den(1) = 949.7109D0 + p*(0.559684D-6 - 0.00097D-12*p) &  
 !      + (t+273.15)*(0.883148 - 0.00228*(t+273.15))  
 !  den(1) = dw_kg + (den(1)-dw_kg)*xmol(2)/p*henry
 !  den(1) = den(1)/avgmw(1)

!****************************** 2 phase S-Pc-kr relation ************************
    if (option%nphase >= 2) then
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
    endif

!    call SaturationFunctionCompute(auxvar%pres,auxvar%sat,kr, &
!                                   ds_dp,dkr_dp, &
!                                   saturation_function, &
!                                   por,perm, &
!                                   option)
    auxvar%kvr(2) = kr(2)/visg     
    auxvar%kvr(1) = kr(1)/visl
    auxvar%vis(2) = visg     
    auxvar%vis(1) = visl
    select case(iphase)
      case(1)
        auxvar%pc =0.D0
      case(2)
        auxvar%pc =0.D0
    end select          

end subroutine MphaseAuxVarCompute_NINC

! ************************************************************************** !

subroutine MphaseAuxVarCompute_WINC(x,delx,auxvar,global_auxvar,iphase,saturation_function, &
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
  type(mphase_auxvar_elem_type) :: auxvar(1:option%nflowdof)
  type(global_auxvar_type) :: global_auxvar

  PetscInt :: iphase

  PetscInt :: n 
  
  do n=1, option%nflowdof
     xx=x;  xx(n)=x(n)+ delx(n)
! ***   note: var_node here starts from 1 to option%flowdof ***
    call  MphaseAuxVarCompute_NINC(xx,auxvar(n),global_auxvar,iphase, &
                                   saturation_function, &
                                   fluid_properties, option)
  enddo

end subroutine MphaseAuxVarCompute_WINC

! ************************************************************************** !

subroutine MphaseAuxVarElemDestroy(auxvar_elem)
  ! 
  ! Deallocates a mphase auxiliary elment object
  ! 
  implicit none

  type(mphase_auxvar_elem_type) :: auxvar_elem

  if (associated(auxvar_elem%xmol)) deallocate(auxvar_elem%xmol)
  nullify(auxvar_elem%xmol)
  if (associated(auxvar_elem%diff))deallocate(auxvar_elem%diff)
  nullify(auxvar_elem%diff)
  if (associated(auxvar_elem%pc))deallocate(auxvar_elem%pc)
  nullify(auxvar_elem%pc)
  if (associated(auxvar_elem%sat))deallocate(auxvar_elem%sat)
  nullify(auxvar_elem%sat)
  if (associated(auxvar_elem%u))deallocate(auxvar_elem%u)
  nullify(auxvar_elem%u)
  if (associated(auxvar_elem%h))deallocate(auxvar_elem%h)
  nullify(auxvar_elem%h)
  if (associated(auxvar_elem%den))deallocate(auxvar_elem%den)
  nullify(auxvar_elem%den)
  if (associated(auxvar_elem%vis))deallocate(auxvar_elem%vis)
   nullify(auxvar_elem%vis)
   if (associated(auxvar_elem%avgmw))deallocate(auxvar_elem%avgmw)
   nullify(auxvar_elem%avgmw)
  if (associated(auxvar_elem%kvr))deallocate(auxvar_elem%kvr)
  nullify(auxvar_elem%kvr)
  if (associated(auxvar_elem%vis))deallocate(auxvar_elem%vis)
  nullify(auxvar_elem%vis)

end subroutine MphaseAuxVarElemDestroy

! ************************************************************************** !

subroutine MphaseAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a mphase auxiliary object
  ! 
  implicit none

  type(mphase_auxvar_type) :: auxvar

  PetscInt :: ielem

  deallocate(auxvar%auxvar_elem(0)%hysdat)
  nullify(auxvar%auxvar_elem(0)%hysdat)

  ! subtract 1 since indexing from 0
  if (associated(auxvar%auxvar_elem)) then
    do ielem = 0, size(auxvar%auxvar_elem) - 1 
      call MphaseAuxVarElemDestroy(auxvar%auxvar_elem(ielem))
    enddo
    deallocate(auxvar%auxvar_elem)
    nullify(auxvar%auxvar_elem)
  endif

end subroutine MphaseAuxVarDestroy

! ************************************************************************** !

subroutine MphaseAuxDestroy(aux)
  ! 
  ! Deallocates a mphase auxiliary object
  ! 
  implicit none

  type(mphase_type), pointer :: aux

  PetscInt :: iaux
  
  if (.not.associated(aux)) return

  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call MphaseAuxVarDestroy(aux%auxvars(iaux))
    enddo
    deallocate(aux%auxvars)
    nullify(aux%auxvars)
  endif
  
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call MphaseAuxVarDestroy(aux%auxvars_bc(iaux))
    enddo
    deallocate(aux%auxvars_bc)
    nullify(aux%auxvars_bc)
  endif
  
  if (associated(aux%auxvars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call MphaseAuxVarDestroy(aux%auxvars_ss(iaux))
    enddo
    deallocate(aux%auxvars_ss)
    nullify(aux%auxvars_ss)
  endif
  
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%mphase_parameter)) then
    if (associated(aux%mphase_parameter%dencpr)) deallocate(aux%mphase_parameter%dencpr)
    nullify(aux%mphase_parameter%dencpr)
    if (associated(aux%mphase_parameter%ckwet)) deallocate(aux%mphase_parameter%ckwet)
    nullify(aux%mphase_parameter%ckwet)
    if (associated(aux%mphase_parameter%sir)) deallocate(aux%mphase_parameter%sir)
    nullify(aux%mphase_parameter%sir)
    deallocate(aux%mphase_parameter)
  endif
  nullify(aux%mphase_parameter)
  if (associated(aux%res_old_AR)) deallocate(aux%res_old_AR)
  if (associated(aux%res_old_FL)) deallocate(aux%res_old_FL)
  if (associated(aux%delx)) deallocate(aux%delx)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine MphaseAuxDestroy

end module Mphase_Aux_module


