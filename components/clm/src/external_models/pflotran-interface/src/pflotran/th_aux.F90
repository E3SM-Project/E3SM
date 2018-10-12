module TH_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: th_itol_scaled_res = 1.d-5
  PetscReal, public :: th_itol_rel_update = UNINITIALIZED_DOUBLE

  type, public :: TH_auxvar_type
    PetscReal :: avgmw
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
    PetscReal :: vis
    PetscReal :: kvr
    PetscReal :: dsat_dp
    PetscReal :: dsat_dt
    PetscReal :: dden_dp
    PetscReal :: dden_dt
    PetscReal :: dkvr_dp
    PetscReal :: dkvr_dt
    PetscReal :: dh_dp
    PetscReal :: dh_dt
    PetscReal :: du_dp
    PetscReal :: du_dt
    PetscReal :: transient_por
    PetscReal :: Ke
    PetscReal :: dKe_dp
    PetscReal :: dKe_dt
    PetscReal :: Dk_eff
    PetscReal :: dDk_eff_dp
    PetscReal :: dDK_eff_dt
    ! for ice
    type(th_ice_type), pointer :: ice
    ! For surface-flow
    type(th_surface_flow_type), pointer :: surface
  end type TH_auxvar_type

  type, public :: th_ice_type
    PetscReal :: Ke_fr
    PetscReal :: dKe_fr_dp
    PetscReal :: dKe_fr_dt
    ! ice
    PetscReal :: sat_ice
    PetscReal :: sat_gas
    PetscReal :: dsat_ice_dp
    PetscReal :: dsat_gas_dp
    PetscReal :: dsat_ice_dt
    PetscReal :: dsat_gas_dt
    PetscReal :: den_ice
    PetscReal :: dden_ice_dp
    PetscReal :: dden_ice_dt
    PetscReal :: u_ice
    PetscReal :: du_ice_dt
    PetscReal :: du_ice_dp
    PetscReal :: den_gas
    PetscReal :: dden_gas_dt
    PetscReal :: dden_gas_dp
    PetscReal :: u_gas
    PetscReal :: du_gas_dt
    PetscReal :: du_gas_dp
    PetscReal :: mol_gas
    PetscReal :: dmol_gas_dt
    PetscReal :: dmol_gas_dp
    ! For DallAmico model
    PetscReal :: pres_fh2o
    PetscReal :: dpres_fh2o_dp
    PetscReal :: dpres_fh2o_dt
  end type th_ice_type
  
  type, public :: th_surface_flow_type
    PetscBool :: surf_wat
    PetscReal :: P_min
    PetscReal :: P_max
    PetscReal :: coeff_for_cubic_approx(4)
    PetscReal :: coeff_for_deriv_cubic_approx(4)
    PetscReal :: range_for_linear_approx(4)
    PetscReal :: dlinear_slope_dT
    PetscBool :: bcflux_default_scheme
  end type th_surface_flow_type

  type, public :: TH_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckdry(:) ! Thermal conductivity (dry)
    PetscReal, pointer :: ckwet(:) ! Thermal conductivity (wet)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: ckfrozen(:) ! Thermal conductivity (frozen soil)
    PetscReal, pointer :: alpha_fr(:) ! exponent frozen
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
  end type TH_parameter_type
  
  type, public :: TH_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(TH_parameter_type), pointer :: TH_parameter
    type(TH_auxvar_type), pointer :: auxvars(:)
    type(TH_auxvar_type), pointer :: auxvars_bc(:)
    type(TH_auxvar_type), pointer :: auxvars_ss(:)
  end type TH_type

  PetscReal, parameter :: epsilon = 1.d-6

  public :: THAuxCreate, THAuxDestroy, &
            THAuxVarComputeNoFreezing, THAuxVarInit, &
            THAuxVarCopy, THAuxVarDestroy

  public :: THAuxVarComputeFreezing

contains

! ************************************************************************** !

function THAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(TH_type), pointer :: THAuxCreate
  
  type(TH_type), pointer :: aux

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

  allocate(aux%TH_parameter)
  nullify(aux%TH_parameter%dencpr)
  nullify(aux%TH_parameter%ckdry)
  nullify(aux%TH_parameter%ckwet)
  nullify(aux%TH_parameter%alpha)
  nullify(aux%TH_parameter%ckfrozen)
  nullify(aux%TH_parameter%alpha_fr)
  nullify(aux%TH_parameter%sir)
  nullify(aux%TH_parameter%diffusion_coefficient)
  nullify(aux%TH_parameter%diffusion_activation_energy)
  
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  allocate(aux%TH_parameter%diffusion_coefficient(option%nphase))
  allocate(aux%TH_parameter%diffusion_activation_energy(option%nphase))
  aux%TH_parameter%diffusion_coefficient = 1.d-9
  aux%TH_parameter%diffusion_activation_energy = 0.d0
 
  THAuxCreate => aux
  
end function THAuxCreate

! ************************************************************************** !

subroutine THAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none
  
  type(TH_auxvar_type) :: auxvar
  type(option_type) :: option
  
  PetscReal :: uninit_value
  uninit_value     = UNINITIALIZED_DOUBLE

  auxvar%avgmw     = uninit_value
  auxvar%h         = uninit_value
  auxvar%u         = uninit_value
  auxvar%pc        = uninit_value
  auxvar%vis       = uninit_value
  auxvar%kvr       = uninit_value
  auxvar%dsat_dp   = uninit_value
  auxvar%dsat_dt   = uninit_value
  auxvar%dden_dp   = uninit_value
  auxvar%dden_dt   = uninit_value
  auxvar%dkvr_dp   = uninit_value
  auxvar%dkvr_dt   = uninit_value
  auxvar%dh_dp     = uninit_value
  auxvar%dh_dt     = uninit_value
  auxvar%du_dp     = uninit_value
  auxvar%du_dt     = uninit_value    
  auxvar%transient_por = uninit_value
  auxvar%Dk_eff    = uninit_value
  auxvar%dDK_eff_dp= uninit_value
  auxvar%dDK_eff_dt= uninit_value
  auxvar%Ke        = uninit_value
  auxvar%dKe_dp    = uninit_value
  auxvar%dKe_dt    = uninit_value

  if (option%use_th_freezing) then
    allocate(auxvar%ice)
    auxvar%ice%Ke_fr     = uninit_value
    auxvar%ice%dKe_fr_dp = uninit_value
    auxvar%ice%dKe_fr_dt = uninit_value
    ! NOTE(bja, 2013-12) always initialize ice variables to zero, even if not used!
    auxvar%ice%sat_ice       = uninit_value
    auxvar%ice%sat_gas       = uninit_value
    auxvar%ice%dsat_ice_dp   = uninit_value
    auxvar%ice%dsat_gas_dp   = uninit_value
    auxvar%ice%dsat_ice_dt   = uninit_value
    auxvar%ice%dsat_gas_dt   = uninit_value
    auxvar%ice%den_ice       = uninit_value
    auxvar%ice%dden_ice_dp   = uninit_value
    auxvar%ice%dden_ice_dt   = uninit_value
    auxvar%ice%u_ice         = uninit_value
    auxvar%ice%du_ice_dt     = uninit_value
    auxvar%ice%du_ice_dp     = uninit_value
    auxvar%ice%den_gas       = uninit_value
    auxvar%ice%dden_gas_dt   = uninit_value
    auxvar%ice%dden_gas_dp   = uninit_value
    auxvar%ice%u_gas         = uninit_value
    auxvar%ice%du_gas_dt     = uninit_value
    auxvar%ice%du_gas_dp     = uninit_value
    auxvar%ice%mol_gas       = uninit_value
    auxvar%ice%dmol_gas_dt   = uninit_value
    auxvar%ice%dmol_gas_dp   = uninit_value
    auxvar%ice%pres_fh2o     = uninit_value
    auxvar%ice%dpres_fh2o_dp = uninit_value
    auxvar%ice%dpres_fh2o_dt = uninit_value
  else
    nullify(auxvar%ice)
  endif
  if (option%surf_flow_on) then
    allocate(auxvar%surface)
    auxvar%surface%surf_wat      = PETSC_FALSE
    auxvar%surface%P_min         = uninit_value
    auxvar%surface%P_max         = uninit_value
    auxvar%surface%coeff_for_cubic_approx(:)       = uninit_value
    auxvar%surface%coeff_for_deriv_cubic_approx(:) = uninit_value
    auxvar%surface%range_for_linear_approx(:)      = uninit_value
    auxvar%surface%dlinear_slope_dT                = uninit_value
    auxvar%surface%bcflux_default_scheme           = PETSC_FALSE
  else
    nullify(auxvar%surface)
  endif
  
end subroutine THAuxVarInit

! ************************************************************************** !

subroutine THAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(TH_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%avgmw = auxvar%avgmw
  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
  auxvar2%vis = auxvar%vis
  auxvar2%kvr = auxvar%kvr
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dsat_dt = auxvar%dsat_dt
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dt = auxvar%dden_dt
  auxvar2%dkvr_dp = auxvar%dkvr_dp
  auxvar2%dkvr_dt = auxvar%dkvr_dt
  auxvar2%dh_dp = auxvar%dh_dp
  auxvar2%dh_dt = auxvar%dh_dt
  auxvar2%du_dp = auxvar%du_dp
  auxvar2%du_dt = auxvar%du_dt  
  auxvar2%transient_por = auxvar%transient_por
  auxvar2%Dk_eff = auxvar%Dk_eff
  auxvar2%Ke = auxvar%Ke
  auxvar2%dKe_dp = auxvar%dKe_dp
  auxvar2%dKe_dt = auxvar%dKe_dt
  if (associated(auxvar%ice)) then
    auxvar2%ice%Ke_fr = auxvar%ice%Ke_fr
    auxvar2%ice%dKe_fr_dp = auxvar%ice%dKe_fr_dp
    auxvar2%ice%dKe_fr_dt = auxvar%ice%dKe_fr_dt
    auxvar2%ice%sat_ice = auxvar%ice%sat_ice 
    auxvar2%ice%sat_gas = auxvar%ice%sat_gas
    auxvar2%ice%dsat_ice_dp = auxvar%ice%dsat_ice_dp
    auxvar2%ice%dsat_gas_dp = auxvar%ice%dsat_gas_dp
    auxvar2%ice%dsat_ice_dt = auxvar%ice%dsat_ice_dt
    auxvar2%ice%dsat_gas_dt = auxvar%ice%dsat_gas_dt
    auxvar2%ice%den_ice = auxvar%ice%den_ice
    auxvar2%ice%dden_ice_dp = auxvar%ice%dden_ice_dp
    auxvar2%ice%dden_ice_dt = auxvar%ice%dden_ice_dt
    auxvar2%ice%u_ice = auxvar%ice%u_ice
    auxvar2%ice%du_ice_dt = auxvar%ice%du_ice_dt
    auxvar2%ice%du_ice_dp = auxvar%ice%du_ice_dp
    auxvar2%ice%pres_fh2o = auxvar%ice%pres_fh2o
    auxvar2%ice%dpres_fh2o_dp = auxvar%ice%dpres_fh2o_dp
    auxvar2%ice%dpres_fh2o_dt = auxvar%ice%dpres_fh2o_dt
    auxvar2%ice%den_gas = auxvar%ice%den_gas
    auxvar2%ice%dden_gas_dt = auxvar%ice%dden_gas_dt
    auxvar2%ice%dden_gas_dp = auxvar%ice%dden_gas_dp
    auxvar2%ice%u_gas = auxvar%ice%u_gas
    auxvar2%ice%du_gas_dt = auxvar%ice%du_gas_dt
    auxvar2%ice%du_gas_dp = auxvar%ice%du_gas_dp
    auxvar2%ice%mol_gas = auxvar%ice%mol_gas
    auxvar2%ice%dmol_gas_dt = auxvar%ice%dmol_gas_dt
    auxvar2%ice%dmol_gas_dp = auxvar%ice%dmol_gas_dp
  endif
  if (associated(auxvar%surface)) then
    auxvar2%surface%surf_wat = auxvar%surface%surf_wat
    auxvar2%surface%P_min = auxvar%surface%P_min
    auxvar2%surface%P_max = auxvar%surface%P_max
    auxvar2%surface%coeff_for_cubic_approx(:) = &
      auxvar%surface%coeff_for_cubic_approx(:)
    auxvar2%surface%coeff_for_deriv_cubic_approx(:) = &
      auxvar%surface%coeff_for_deriv_cubic_approx(:)
    auxvar2%surface%range_for_linear_approx(:) = &
      auxvar%surface%range_for_linear_approx(:)
    auxvar2%surface%dlinear_slope_dT = auxvar%surface%dlinear_slope_dT
    auxvar2%surface%bcflux_default_scheme = &
      auxvar%surface%bcflux_default_scheme
  endif

end subroutine THAuxVarCopy

! ************************************************************************** !

subroutine THAuxVarComputeNoFreezing(x,auxvar,global_auxvar, &
                                     material_auxvar, &
                                     iphase, &
                                     characteristic_curves,    &
                                     th_parameter, ithrm, &
                                     option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type), pointer :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
  class(material_auxvar_type) :: material_auxvar

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp, dkr_dsl
  PetscReal :: dvis_dt, dvis_dp
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: Ke
  PetscReal :: alpha
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: aux(1)

  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%kvr = 0.d0
  kr = 0.d0
 
  global_auxvar%pres = x(1)
  global_auxvar%temp = x(2)
 
! auxvar%pc = option%reference_pressure - auxvar%pres
  auxvar%pc = min(option%reference_pressure - global_auxvar%pres(1), &
                  characteristic_curves%saturation_function%pcmax)

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (auxvar%pc > 0.d0) then
    iphase = 3

    call characteristic_curves%saturation_function%Saturation(auxvar%pc, &
         global_auxvar%sat(1), ds_dp, option)

    call characteristic_curves%liq_rel_perm_function%RelativePermeability( &
         global_auxvar%sat(1), &
         kr, dkr_dsl, option)
    dkr_dp = dkr_dsl*ds_dp

    dpw_dp = 0.d0
  else
    iphase = 1
    auxvar%pc = 0.d0
    global_auxvar%sat(1) = 1.d0  
    kr = 1.d0    
    pw = global_auxvar%pres(1)
    dpw_dp = 1.d0
  endif  

  ! may need to compute dpsat_dt to pass to VISW
  call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,dpsat_dt,ierr)
  call EOSWaterEnthalpy(global_auxvar%temp,pw,hw,hw_dp,hw_dt,ierr)
  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,ierr)
    call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,dpsat_dt,visl, &
                           dvis_dt,dvis_dp,ierr)
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(global_auxvar%temp,pw,aux, &
                            dw_kg,dw_mol,dw_dp,dw_dt,ierr)
    call EOSWaterViscosityExt(global_auxvar%temp,pw,sat_pressure,dpsat_dt,aux, &
                              visl,dvis_dt,dvis_dp,ierr)
  endif
  ! J/kmol -> whatever units
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
  
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%kvr = kr/visl
  
  auxvar%vis = visl
  auxvar%dsat_dp = ds_dp
  auxvar%dsat_dt = 0.d0

  auxvar%dden_dt = dw_dt
  auxvar%dden_dp = dw_dp
  
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  if (iphase < 3) then !kludge since pw is constant in the unsat zone
    auxvar%dh_dp = hw_dp
    auxvar%du_dp = hw_dp - (dpw_dp/dw_mol-pw/(dw_mol*dw_mol)*dw_dp)*option%scale
  else
    auxvar%dh_dp = 0.d0
    auxvar%du_dp = 0.d0
  endif

  auxvar%dh_dt = hw_dt
  auxvar%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt
  
  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)

  !unfrozen soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  auxvar%Ke = Ke

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk_dry + (Dk - Dk_dry)*Ke

  ! Derivative of soil Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
  auxvar%dKe_dt = 0.d0

end subroutine THAuxVarComputeNoFreezing

! ************************************************************************** !

subroutine THAuxVarComputeFreezing(x, auxvar, global_auxvar, &
                                   material_auxvar,          &
                                   iphase,                   &
                                   characteristic_curves,    &
                                   th_parameter, ithrm,      &
                                   option)
  ! 
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 11/16/11
  ! 

!sk: Not sure if we need por, perm
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type), pointer :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw, dw_kg, dw_mol, hw, sat_pressure, visl
  PetscReal :: kr, dkr_dp, dkr_dt
  PetscReal :: dvis_dt, dvis_dp
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: liq_saturation, ice_saturation, gas_saturation
  PetscReal :: dsl_dp, dsl_dt
  PetscReal :: dsg_dp, dsg_dt
  PetscReal :: dsi_dp, dsi_dt
  PetscReal :: den_ice, dden_ice_dt, dden_ice_dp
  PetscReal :: u_ice, du_ice_dt
  PetscReal :: p_th

  PetscReal :: p_g, tk_g
  PetscReal :: p_sat
  PetscReal :: mol_g
  PetscReal :: C_g
  PetscReal :: dmolg_dt, dmolg_dp
  PetscReal, parameter :: C_a  = 1.86d-3   ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3  ! in MJ/kg/K

  PetscReal :: Ke
  PetscReal :: Ke_fr
  PetscReal :: alpha
  PetscReal :: alpha_fr
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: Dk_ice
  character(len=MAXSTRINGLENGTH) :: error_string

  PetscBool :: DTRUNC_FLAG = PETSC_TRUE  ! option for truncating deriatives to zero at bounds (default: TRUE)
  PetscReal :: tcmin, pcmax, pcmin
  PetscReal :: t_trunc, p_trunc, dt_trunc, dp_trunc


  !
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%kvr = 0.d0
   
  global_auxvar%pres(1) = x(1)
  global_auxvar%temp    = x(2)



!***************  P/T Bounds ***********************************************

  ! do the trunction of derivatives when assigning temporary values (local) to global variables (auxvar%, global_auxvar%)
  ! for P/T of out of bounds
  dt_trunc = 1.d0
  dp_trunc = 1.d0

  ! general bounds of P/T. We may need to further limit bounds for 3-phase water properties.
  if (global_auxvar%temp>=100.d0 .or. global_auxvar%temp<=-273.d0) dt_trunc = 0.d0
  global_auxvar%temp    = min(100.d0, max(-273.d0, global_auxvar%temp))
  if (global_auxvar%pres(1)>=16.54d6) dp_trunc = 0.d0           ! 16.54 MPa is upper limit for using IFC67 EOS.
  global_auxvar%pres(1) = min(16.54d6, global_auxvar%pres(1))
  
  ! Check if the capillary pressure is less than -100MPa, which also limit the lower limit of pres(1).
  pcmax = abs(characteristic_curves%saturation_function%pcmax)  ! non-negative
  if (global_auxvar%pres(1) - option%reference_pressure <= -pcmax) then
    dp_trunc = 0.d0
    global_auxvar%pres(1) = -pcmax + option%reference_pressure
  endif
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)  ! always non-negative



!***************  Characteristic Curves ***********************************
  
  ! using modules in 'characteristic_curves.F90' to calculate needed variables:
  ! saturations and derivatives for 3 phases
  call THAuxVarComputeCharacteristicCurves(global_auxvar%pres(1), global_auxvar%temp,   &
                                           characteristic_curves,                       &
                                           liq_saturation,  dsl_dp, dsl_dt,             &
                                           ice_saturation,  dsi_dp, dsi_dt,             &
                                           auxvar%ice%pres_fh2o,                        &
                                           auxvar%ice%dpres_fh2o_dp,                    &
                                           auxvar%ice%dpres_fh2o_dt,                    &
                                           gas_saturation,  dsg_dp, dsg_dt,             &
                                           kr,              dkr_dp, dkr_dt,             &
                                           option)



!***************  3-phase water properties **********************************************************

  ! ----- Liq. water ---------------------------------------------------------
  auxvar%avgmw = FMWH2O

  pw     = option%reference_pressure
  pcmin  = 0.d0  ! near-saturation PC zone (hint: a small positive value may be helpful ?)
  if (auxvar%pc > pcmin) then
    iphase    = 3
    dpw_dp    = 0.d0
  else
    iphase    = 1
    auxvar%pc = 0.d0
    pw        = global_auxvar%pres(1)
    dpw_dp    = 1.d0
  endif

  global_auxvar%sat(1) = liq_saturation
  auxvar%dsat_dp = dsl_dp
  auxvar%dsat_dt = dsl_dt
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    auxvar%dsat_dp = auxvar%dsat_dp * dp_trunc
    auxvar%dsat_dt = auxvar%dsat_dt * dt_trunc
  endif

  ! F.-M. Yuan (2016-07-10)
  ! Liq. water Enthalpy by IFC67 eq. may have a tiny offset of temperature as following
  ! i.e. @ t of ~-0.02404oC/1atm, hw ~ 0.d0; beyond that, hw is negative with a large negative derivative.
  ! So for a threshold of 0 ~ -0.02oC, it will give a small positive value of 'hw' for suppercooled liq. water
  tcmin = -0.02d0
  ! a note here: limit for p/t according to IFC67 p/t ranges.
  ! This may be helpful to solve small-timing issue

  call EOSWaterDensity(max(tcmin, global_auxvar%temp),      &
                        pw, dw_kg, dw_mol, dw_dp, dw_dt,ierr)
  if (DTRUNC_FLAG) dw_dp = dw_dp * dpw_dp  ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation
  if (DTRUNC_FLAG .and. global_auxvar%temp<tcmin) dw_dt = 0.d0

  call EOSWaterEnthalpy(max(tcmin, global_auxvar%temp),     &
                        pw, hw,hw_dp,hw_dt,ierr)
  if (DTRUNC_FLAG .and. global_auxvar%temp<tcmin) hw_dt = 0.d0

  ! J/kmol -> MJ/kmol
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
  if(DTRUNC_FLAG) hw_dp = hw_dp * dpw_dp  ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation

  ! A note here (F.-M. Yuan: 2017-01-17)
  ! The Viscosity Eq. shows that: temp< ~ -63oC (1atm), 'visl' sharply increases starting from ~ 1.e-2 order.
  ! 'visl' ~ 0. around -133oC, which produces 'inf' for 'kr'
  tcmin = -63.d0

  call EOSWaterSaturationPressure(max(tcmin,global_auxvar%temp), sat_pressure, dpsat_dt, ierr)
  ! the lowest Tk of 200 for vapor exists in EOS-h2o phase-diagram, but here make it consistent with 'Viscosity'
  if(DTRUNC_FLAG .and. global_auxvar%temp<=tcmin) dpsat_dt = 0.d0

  call EOSWaterViscosity(max(tcmin,global_auxvar%temp), pw,   &
                         sat_pressure, dpsat_dt,   &
                         visl, dvis_dt,dvis_dp, ierr)
  if(DTRUNC_FLAG .and. global_auxvar%temp<=tcmin) dvis_dt = 0.d0
  if(DTRUNC_FLAG) dvis_dp = dvis_dp*dpw_dp  ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation

  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%kvr = kr/visl
  auxvar%vis = visl
  auxvar%dden_dt = dw_dt
  auxvar%dden_dp = dw_dp
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity  
!  auxvar%dkvr_dt = -kr/(visl*visl)*(dvis_dt + dvis_dpsat*dpsat_dt) + dkr_dt/visl
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt + dkr_dt/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  auxvar%dh_dp = hw_dp
  auxvar%du_dp = hw_dp - (dpw_dp/dw_mol - pw/(dw_mol*dw_mol)*dw_dp)* option%scale
  auxvar%dh_dt = hw_dt
  auxvar%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt

  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    auxvar%dden_dt = auxvar%dden_dt * dt_trunc
    auxvar%dden_dp = auxvar%dden_dp * dp_trunc
    auxvar%dkvr_dt = auxvar%dkvr_dt * dt_trunc
    auxvar%dkvr_dp = auxvar%dkvr_dp * dp_trunc
    auxvar%dh_dt = auxvar%dh_dt * dt_trunc
    auxvar%dh_dp = auxvar%dh_dp * dp_trunc
    auxvar%du_dt = auxvar%du_dt * dt_trunc
    auxvar%du_dp = auxvar%du_dp * dp_trunc
  endif

  ! ----- ice water ---------------------------------------------------------

  auxvar%ice%sat_ice     = ice_saturation
  auxvar%ice%dsat_ice_dp = dsi_dp
  auxvar%ice%dsat_ice_dt = dsi_dt

  ! for ice Ih, Tk limit is ~127K (-146oC) in EOS-h2o phase-diagram
  tcmin = -146.d0

  call EOSWaterDensityIce(max(tcmin, global_auxvar%temp), pw,                    &
                          den_ice, dden_ice_dT, dden_ice_dp, ierr)
  if(DTRUNC_FLAG) dden_ice_dp = dden_ice_dp * dpw_dp  ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation
  if(DTRUNC_FLAG .and. global_auxvar%temp<=tcmin) dden_ice_dT = 0.d0

  call EOSWaterInternalEnergyIce(max(tcmin, global_auxvar%temp), u_ice, du_ice_dT)
  if(DTRUNC_FLAG .and. global_auxvar%temp<=tcmin) du_ice_dT = 0.d0

  auxvar%ice%den_ice     = den_ice
  auxvar%ice%dden_ice_dt = dden_ice_dT
  auxvar%ice%dden_ice_dp = dden_ice_dP
  auxvar%ice%u_ice     = u_ice*1.d-3              !kJ/kmol --> MJ/kmol
  auxvar%ice%du_ice_dt = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K 
  auxvar%ice%du_ice_dp = 0.d0

  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    auxvar%ice%dsat_ice_dp = auxvar%ice%dsat_ice_dp * dp_trunc
    auxvar%ice%dsat_ice_dt = auxvar%ice%dsat_ice_dt * dt_trunc
    auxvar%ice%dden_ice_dp = auxvar%ice%dden_ice_dp * dp_trunc
    auxvar%ice%dden_ice_dt = auxvar%ice%dden_ice_dt * dt_trunc
    auxvar%ice%du_ice_dp   = auxvar%ice%du_ice_dp * dp_trunc
    auxvar%ice%du_ice_dt   = auxvar%ice%du_ice_dt * dt_trunc
  endif

  ! ----- Air (incl. vapor) ---------------------------------------------------------

  auxvar%ice%sat_gas = gas_saturation
  auxvar%ice%dsat_gas_dp = dsg_dp
  auxvar%ice%dsat_gas_dt = dsg_dt

  ! Calculate the values and derivatives for vapor density and internal energy

  p_g      = pw
  ! the lowest Tk of ~200K for vapor exists in EOS-h2o phase-diagram
  tcmin    = -73.d0
  tk_g     = max(tcmin,global_auxvar%temp)+273.15d0

  auxvar%ice%den_gas     = p_g/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3                ! in kmol/m3 for all air-mixture
  auxvar%ice%dden_gas_dt = -p_g/(IDEAL_GAS_CONSTANT*tk_g**2)*1.d-3
  auxvar%ice%dden_gas_dp = 1.d0/(IDEAL_GAS_CONSTANT*tk_g)*1.d-3
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    ! the following is a MUST for reducing tiny-time step (NOT sure why).
    auxvar%ice%dden_gas_dp = auxvar%ice%dden_gas_dp*dpw_dp    ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation
    if (tk_g<=tcmin) auxvar%ice%dden_gas_dt = 0.d0
  endif

  call EOSWaterSaturationPressure(tk_g-273.15d0, p_sat, dpsat_dt, ierr)
  mol_g    = p_sat/p_g
  dmolg_dt = dpsat_dt/p_g
  dmolg_dp = -p_sat/p_g/p_g
  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    dmolg_dp = dmolg_dp * dpw_dp                              ! w.r.t from 'pw' to 'pres(1)' upon soil total saturation
    if (tk_g<=tcmin) then
      dpsat_dt = 0.d0
      dmolg_dt = 0.d0
    endif
  endif

  C_g                    = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR          ! in MJ/kmol/K
  auxvar%ice%u_gas       = C_g*tk_g                                               ! in MJ/kmol

#if 0
! 2017-01-30: tests show that the following may cause tiny-time recovering issue when re-freezing (so OFF now).
  ! NOTE: vapor 'mol_gas' is included in 'den_gas', 'u_gas'
  ! (because 'mol_g', fraction of vapor in air-mixture, going to be as multiplier in all 'gas' calculations in 'th.F90')
  ! so, if for air-mixture, may assign this to 1.0 ( appears very helpful to reduce time-step)
  !     if totally ignore gas (all), may assign this to 0.0
#ifdef NO_VAPOR_DIFFUSION
  auxvar%ice%mol_gas     = 0.d0   ! no gas (inc. vapor)
#else
  auxvar%ice%mol_gas     = 1.d0   ! air-mixture
#endif
  auxvar%ice%dmol_gas_dt = 0.d0
  auxvar%ice%dmol_gas_dp = 0.d0
#else
  ! vapor-only (because 'mol_g', fraction of vapor in air-mixture, going to be as multiplier in all 'gas' calculations in 'th.F90'
  auxvar%ice%mol_gas     = mol_g
  auxvar%ice%dmol_gas_dt = dmolg_dt
  auxvar%ice%dmol_gas_dp = dmolg_dp
#endif
  auxvar%ice%du_gas_dt   = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)*tk_g
  auxvar%ice%du_gas_dp   = (C_wv*FMWH2O-C_a*FMWAIR)*dmolg_dp*tk_g


  ! F.-M. Yuan (2017-01-18): truncating 'derivatives' at the bounds
  if(DTRUNC_FLAG) then
    ! something happened in the following, causing difficulties to re-freeze soils (TODO - checking: 2017-01-18)
    ! (maybe: except for DALL_AMICO model)
    ! 02-07-2017: It may be relevant to soil compressibility.
    !             This also causes LARGE temperature oscillation during F/T - with one layer supper HOT while other supper COLD.
    !if (option%ice_model == DALL_AMICO) &
    auxvar%ice%dsat_gas_dp = auxvar%ice%dsat_gas_dp * dp_trunc

    auxvar%ice%dsat_gas_dt = auxvar%ice%dsat_gas_dt * dt_trunc
    auxvar%ice%dden_gas_dp = auxvar%ice%dden_gas_dp * dp_trunc
    auxvar%ice%dden_gas_dt = auxvar%ice%dden_gas_dt * dt_trunc
    auxvar%ice%dmol_gas_dp = auxvar%ice%dmol_gas_dp * dp_trunc
    auxvar%ice%dmol_gas_dt = auxvar%ice%dmol_gas_dt * dt_trunc
    auxvar%ice%du_gas_dp = auxvar%ice%du_gas_dp * dp_trunc
    auxvar%ice%du_gas_dt = auxvar%ice%du_gas_dt * dt_trunc
  endif


  ! ----- Thermal conductivity ---------------------------------------------------------

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  alpha_fr = th_parameter%alpha_fr(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)
  Dk_ice = th_parameter%ckfrozen(ithrm)

  !Soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  Ke_fr = (auxvar%ice%sat_ice + epsilon)**(alpha_fr)
  auxvar%Ke = Ke
  auxvar%ice%Ke_fr = Ke_fr

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk*Ke + Dk_ice*Ke_fr + (1.d0 - Ke - Ke_fr)*Dk_dry

  ! Derivative of Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
  auxvar%dKe_dt = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dt
  auxvar%ice%dKe_fr_dt = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dt
  auxvar%ice%dKe_fr_dp = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dp

  ! derivative of 'Dk_eff': Dk_eff = Dk_dry + (Dk-Dk_dry)*Ke + (Dk_ice-Dk_dry)*Ke_fr
  auxvar%dDk_eff_dp = (Dk-Dk_dry)*auxvar%dKe_dp + (Dk_ice-Dk_dry)*auxvar%ice%dKe_fr_dp
  auxvar%dDk_eff_dt = (Dk-Dk_dry)*auxvar%dKe_dt + (Dk_ice-Dk_dry)*auxvar%ice%dKe_fr_dt

#if 0
  ! F.-M. Yuan (2017-1-30): comment the following out, so that DALL_AMICO model works as others.
  if (option%ice_model == DALL_AMICO) then
    auxvar%ice%den_ice = dw_mol
    auxvar%ice%dden_ice_dt = auxvar%dden_dt
    auxvar%ice%dden_ice_dp = auxvar%dden_dp
!    auxvar%ice%u_ice = auxvar%u  ! commented out by S.Karra 06/02/14. setting
!    internal energy of ice and water might not be correct.
!    auxvar%ice%du_ice_dt = auxvar%du_dt

    auxvar%ice%sat_gas       = 0.d0
    auxvar%ice%dsat_gas_dp   = 0.d0
    auxvar%ice%dsat_gas_dt   = 0.d0
    auxvar%ice%den_gas       = 0.d0
    auxvar%ice%dden_gas_dt   = 0.d0
    auxvar%ice%u_gas         = 0.d0
    auxvar%ice%du_gas_dt     = 0.d0
    auxvar%ice%mol_gas       = 0.d0
    auxvar%ice%dmol_gas_dt   = 0.d0
  endif
#endif

end subroutine THAuxVarComputeFreezing

! ************************************************************************** !
subroutine THAuxVarComputeCharacteristicCurves( pres_l,  tc,                &
                                    characteristic_curves,                  &
                                    sl,  dsl_dpl, dsl_dt,                   &
                                    si,  dsi_dpl, dsi_dt,                   &
                                    ice_presl, ice_presl_dpl, ice_presl_dt, &
                                    sg,  dsg_dpl, dsg_dt,                   &
                                    kr,  dkr_dpl, dkr_dt,                   &
                                    option)
  !
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  !
  ! Revised by fengming Yuan @03-08-2016/CCSI-ONRL
  ! NOTE: (1) ice_model 'PAINTER_EXPLICIT', or, 'PAINTER_KARRA_EXPLICIT',
  !             or, 'PAINTER_KARRA_EXPLICIT_SMOOTH', or, 'DALL_AMICO'
  !       (2) ANY saturation_function, from 'Characteristic_Curves_module'
  !       (3) ANY permissivity function, from 'Characteristci_Curves_module' as well

  use Option_module
  use Characteristic_Curves_module

  implicit none

  type(option_type) :: option
  PetscReal, intent(in) :: pres_l   ! unit: Pa (liq water-air interface pressure: -pc+reference_pressure)
  PetscReal, intent(in) :: tc       ! unit: oC
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal, intent(out) :: sl,  dsl_dpl, dsl_dt
  PetscReal, intent(out) :: si,  dsi_dpl, dsi_dt
  PetscReal, intent(out) :: ice_presl, ice_presl_dpl, ice_presl_dt
  PetscReal, intent(out) :: sg,  dsg_dpl, dsg_dt
  PetscReal, intent(out) :: kr,  dkr_dpl, dkr_dt

  ! local variables
  PetscErrorCode :: ierr

  PetscReal :: pc
  PetscReal :: sli, dsli_dpl, xplice, dxplice_dpl, dxplice_dt, slx, dslx_dx
  PetscReal :: dkr_dsl

  PetscReal :: se, dse_dpc, function_A, dfunc_A_dt, function_B, dfunc_B_dpl

  ! ----------------
  ! (0) inputs
  pc = max(0.d0, option%reference_pressure - pres_l)   ! always non-negative (0 = saturated)
  if (pc > abs(characteristic_curves%saturation_function%pcmax)) then
    pc = characteristic_curves%saturation_function%pcmax
  endif

  !
  ! (1) saturations
  sl      = 0.d0  !all init to zero
  dsl_dpl = 0.d0
  dsl_dt  = 0.d0
  si      = 0.d0
  dsi_dpl = 0.d0
  dsi_dt  = 0.d0
  sg      = 1.d0
  dsg_dpl = 0.d0
  dsg_dt  = 0.d0

  ice_presl    = pres_l
  ice_presl_dpl= 1.d0
  ice_presl_dt = 0.d0

  ! initial liq. saturation and derivatives
  call characteristic_curves%saturation_function%Saturation(pc, sli, dsli_dpl, option)  ! a note: in CC modules, already w.r.t 'pres_l' by dpc_dpres = -1.d0
  sl      = sli
  dsl_dpl = dsli_dpl
  dsl_dt  = 0.d0

  ! if ice module turns on, 2-phase saturation recalculated (liq. and ice) under common 'pres_l' and 'tc'
  xplice = pc   ! liq. water PC at ice-liq interface (positive, if no ice it's pc by default)
  if (option%use_th_freezing) then

    call characteristic_curves%saturation_function%IceCapillaryPressure(pres_l, tc, &
                                   xplice, dxplice_dpl, dxplice_dt, option) ! w.r.t 'pressure' already in %IceCapillaryPressure()

    select case (option%ice_model)
      case (PAINTER_EXPLICIT)

        ! NOTE: in 'saturation_function.F90': SatFuncComputeIcePExplicit(), 'se' and 'dse_dpc' are used
        ! while in 'characteristic_curves.F90': SF_VG_Saturaton(), the following are output:
        !   sl = this%Sr + (1.d0-this%Sr)*Se
        !   dsl_dpl = -(1.d0-this%Sr)*dSe_dpc

        function_B = 1.d0
        dfunc_B_dpl= 0.d0
        if (pc>0.d0) then
          se = (sl - characteristic_curves%saturation_function%Sr) &
               /(1.0d0 - characteristic_curves%saturation_function%Sr)
          dse_dpc = (-dsl_dpl)/(1.0d0 - characteristic_curves%saturation_function%Sr)

          function_B = 1.d0/se
          dfunc_B_dpl= 1.d0/(se*se)*dse_dpc
        endif
        !
        function_A = 1.d0
        dfunc_A_dt = 0.d0
        if(tc<0.d0) then
          call characteristic_curves%saturation_function%Saturation(xplice, slx, dslx_dx, option)
          se = (slx - characteristic_curves%saturation_function%Sr) &
              /(1.0d0 - characteristic_curves%saturation_function%Sr)
          dse_dpc = (-dslx_dx)/(1.0d0 - characteristic_curves%saturation_function%Sr)

          if (se>0.d0) then
            function_A = 1.d0/se
            ! dfunc_A_dt = dfuncA_dslx * dslx_dx * dxplice_dt   (note: dslx_dx = dslx_dxplice)
            dfunc_A_dt = 1.d0/(se*se)* dse_dpc
            dfunc_A_dt = dfunc_A_dt * (-dxplice_dt)
          endif
        endif

        sl = 1.d0/(function_A + function_B - 1.d0)
        sg = sl*(function_B - 1.d0)
        si = sl*(function_A - 1.d0)

        dsl_dpl = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_B_dpl)
        dsl_dt  = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_A_dt)

        dsg_dpl = dsl_dpl*(function_B - 1.d0) + sl*dfunc_B_dpl
        dsg_dt  = dsl_dt*(function_B - 1.d0)

        dsi_dpl = dsl_dpl*(function_A - 1.d0)
        dsi_dt  = dsl_dt*(function_A - 1.d0) + sl*dfunc_A_dt

      case (PAINTER_KARRA_EXPLICIT, PAINTER_KARRA_EXPLICIT_SMOOTH)

#if 0
        ! (TODO-checking) when coupled with CLM, the following not works well
        ! although it's theoretically better (because it's used for calculating 'dphi' for water Darcy flux)
        ice_presl    = (pres_l - pc) - xplice
        ice_presl_dpl= dxplice_dpl
        ice_presl_dt = dxplice_dt
#endif
        call characteristic_curves%saturation_function%Saturation(xplice, slx, dslx_dx, option)

        ! liq. saturation and its derivatives, with ice-adjusted capillary pressure
        sl      = slx
        dsl_dt  = -dslx_dx*dxplice_dt
        dsl_dpl = dslx_dx*dxplice_dpl

        ! ice satuation and its derivatives
        if (sli>0.d0) then
          si     = 1.d0 - sl/sli             ! P.-K. Eq.(19)
          dsi_dt = -1.d0/sli*dsl_dt          ! dsli_dt = 0 (see above)
          dsi_dpl= (sl*dsli_dpl-sli*dsl_dpl)/(sli**2)
        else
          si     = 0.d0
          dsi_dt = 0.d0
          dsi_dpl= 0.d0
        endif

        ! gas component (as difference)
        sg      = 1.d0 - sl - si
        dsg_dpl = -dsl_dpl - dsi_dpl
        dsg_dt  = -dsl_dt - dsi_dt

      case (DALL_AMICO)

        ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
        ! rewritten following 'saturation_function.F90:SatFuncComputeIceDallAmico()'
        ! NOTE: here calculate 'saturations and its derivatives'

        call characteristic_curves%saturation_function%Saturation(xplice, slx, dslx_dx, option)   ! Pc1 ---> S1, but '_dx' is w.r.t. '_dpres'

#if 0
        ! (TODO-checking) when coupled with CLM, the following not works well
        ! although it's theoretically better (because it's used for calculating 'dphi' for water Darcy flux)
        ice_presl    = (pres_l - pc) - xplice
        ice_presl_dpl= dxplice_dpl
        ice_presl_dt = dxplice_dt
#endif

        !
        ! liq. saturation and its derivatives, with ice-adjusted capillary pressure
        sl     = slx
        dsl_dpl= dslx_dx * dxplice_dpl
        dsl_dt = -dslx_dx * dxplice_dt

        ! ice satuation and its derivatives
        si     = sli - sl
        dsi_dpl= dsli_dpl - dsl_dpl
        dsi_dt = -dsl_dt                 ! dsli_dt = 0 (see above)

        ! gas phase
        sg      = 1.d0 - sl - si
        dsg_dpl = -dsl_dpl - dsi_dpl
        dsg_dt  = -dsl_dt - dsi_dt

      case default
        option%io_buffer = 'Ice module NOT recognized'
        call printErrMsg(option)

    end select
  endif ! 'option%use_th_freezing'

  ! Check for bounds on saturations
  if ((sl-1.d0)>1.d-15 .or. sl<-1.d-15) then
    print *, tc, pc, sl, si, sg, sli, xplice
    option%io_buffer = 'TH with ice mode: Liquid Saturation error: >1 or <0'
    call printErrMsg(option)
  endif
  if ((si-1.d0)>1.d-15 .or. si<-1.d-15) then
    print *, tc, pc, sl, si, sg, sli, xplice
    option%io_buffer = 'TH with ice mode: ICE Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if ((sg-1.d0)>1.d-15 .or. sg<-1.d-15) then
    print *, tc, pc, sl, si, sg, sli, xplice
    option%io_buffer = 'TH with ice mode: Gas Saturation error:  >1 or <0'
    call printErrMsg(option)
  endif
  if (abs((sl + si + sg)-1.d0)>1.d-10) then
    option%io_buffer = 'TH with ice mode: Saturation not summed to 1 '
    call printErrMsg(option)
  endif

  ! (2) relative permissivity of liq. water in multiple-phase mixture
  kr      = 0.d0  !all initialized to zero
  dkr_dsl = 0.d0
  dkr_dt  = 0.d0
  dkr_dpl = 0.d0

  call characteristic_curves%liq_rel_perm_function%RelativePermeability(sl, kr, dkr_dsl, option)
  dkr_dpl= dkr_dsl*dsl_dpl
  dkr_dt = dkr_dsl*dsl_dt

end subroutine THAuxVarComputeCharacteristicCurves

! ************************************************************************** !

subroutine THAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a TH auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(TH_auxvar_type) :: auxvar
  
  if (associated(auxvar%ice)) deallocate(auxvar%ice)
  nullify(auxvar%ice)
  if (associated(auxvar%surface)) deallocate(auxvar%surface)
  nullify(auxvar%surface)
  
end subroutine THAuxVarDestroy

! ************************************************************************** !

subroutine THAuxDestroy(aux)
  ! 
  ! Deallocates a TH auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(TH_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    call THAuxVarDestroy(aux%auxvars(iaux))
  enddo  
  do iaux = 1, aux%num_aux_bc
    call THAuxVarDestroy(aux%auxvars_bc(iaux))
  enddo  
  do iaux = 1, aux%num_aux_ss
    call THAuxVarDestroy(aux%auxvars_ss(iaux))
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
  if (associated(aux%TH_parameter)) then
    if (associated(aux%TH_parameter%diffusion_coefficient)) &
      deallocate(aux%TH_parameter%diffusion_coefficient)
    nullify(aux%TH_parameter%diffusion_coefficient)
    if (associated(aux%TH_parameter%diffusion_activation_energy)) &
      deallocate(aux%TH_parameter%diffusion_activation_energy)
    nullify(aux%TH_parameter%diffusion_activation_energy)
    if (associated(aux%TH_parameter%dencpr)) deallocate(aux%TH_parameter%dencpr)
    nullify(aux%TH_parameter%dencpr)
    if (associated(aux%TH_parameter%ckwet)) deallocate(aux%TH_parameter%ckwet)
    nullify(aux%TH_parameter%ckwet)
    if (associated(aux%TH_parameter%ckdry)) deallocate(aux%TH_parameter%ckdry)
    nullify(aux%TH_parameter%ckdry)
    if (associated(aux%TH_parameter%alpha)) deallocate(aux%TH_parameter%alpha)
    nullify(aux%TH_parameter%alpha)
    ! ice
    if (associated(aux%TH_parameter%ckfrozen)) deallocate(aux%TH_parameter%ckfrozen)
    nullify(aux%TH_parameter%ckfrozen)
    if (associated(aux%TH_parameter%alpha_fr)) deallocate(aux%TH_parameter%alpha_fr)
    nullify(aux%TH_parameter%alpha_fr)

    if (associated(aux%TH_parameter%sir)) deallocate(aux%TH_parameter%sir)
    nullify(aux%TH_parameter%sir)
  endif
  nullify(aux%TH_parameter)
  
  deallocate(aux)
  nullify(aux)  

  end subroutine THAuxDestroy

end module TH_Aux_module
