module General_Derivative_module

  
  use PFLOTRAN_Constants_module
  use Option_module
  use General_Aux_module
  use General_Common_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: GeneralDerivativeDriver
  
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

contains
  
! ************************************************************************** !

subroutine GeneralDerivativeDriver(option)
          
  use Characteristic_Curves_module

  implicit none
  
  type(option_type), pointer :: option
  
  PetscInt, parameter :: ACCUMULATION = 1
  PetscInt, parameter :: INTERIOR_FLUX = 2
  PetscInt, parameter :: BOUNDARY_FLUX = 3
  PetscInt, parameter :: SRCSINK = 4
  
  PetscInt :: itype

  PetscInt :: istate
  PetscReal :: xx(3)
  PetscReal :: pert(3)
  type(general_auxvar_type), pointer :: general_auxvar(:)
  type(global_auxvar_type), pointer :: global_auxvar(:)
  class(material_auxvar_type), pointer :: material_auxvar(:)

  PetscInt :: istate2
  PetscReal :: xx2(3)
  PetscReal :: pert2(3)
  PetscReal :: scale2
  type(general_auxvar_type), pointer :: general_auxvar2(:)
  type(global_auxvar_type), pointer :: global_auxvar2(:)
  class(material_auxvar_type), pointer :: material_auxvar2(:)
  
  PetscInt :: ibndtype(3)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)
  PetscReal :: auxvars(3) ! from aux_real_var array
  
  PetscInt :: flow_src_sink_type
  PetscReal :: qsrc(3)
  PetscReal :: srcsink_scale
  
  class(characteristic_curves_type), pointer :: characteristic_curves
  type(material_parameter_type), pointer :: material_parameter
  type(general_parameter_type), pointer :: general_parameter

  nullify(characteristic_curves)
  nullify(material_parameter)
  nullify(general_parameter)
  call GeneralDerivativeSetFlowMode(option)
  call GeneralDerivativeSetupEOS(option)
  
  call GeneralDerivativeSetup(general_parameter, &
                              characteristic_curves, &
                              material_parameter,option)  
  option%flow_dt = 1.d0
!  itype = ACCUMULATION
!  itype = INTERIOR_FLUX
!  itype = BOUNDARY_FLUX
  itype = SRCSINK
  
!  istate = LIQUID_STATE
!  istate = GAS_STATE
  istate = TWO_PHASE_STATE
  select case(istate)
    case(LIQUID_STATE)
      xx(1) = 1.d6
      xx(2) = 1.d-6
      xx(3) = 30.d0
!      xx(3) = 100.d0
    case(GAS_STATE)
!      xx(1) = 1.d4
!      xx(2) = 0.997d4
!      xx(3) = 15.d0
      xx(1) = 1.d6
      xx(2) = 0.997d6
      xx(3) = 30.d0
!      xx(1) = 1.d7
!      xx(2) = 0.997d7
!      xx(3) = 85.d0
    case(TWO_PHASE_STATE)
      xx(1) = 1.d6
      xx(2) = 0.5d0
      xx(3) = 30.d0
  end select

  call GeneralDerivativeSetupAuxVar(istate,xx,pert,general_auxvar, &
                                    global_auxvar,material_auxvar, &
                                    characteristic_curves, &
                                    option)  
  
  istate2 = istate
!  istate2 = LIQUID_STATE
!  istate2 = GAS_STATE
!  istate2 = TWO_PHASE_STATE
  ! scales must range (0.001 - 1.999d0)
  scale2 = 1.d0
!  scale2 = 1.d0 - 1.d-14
!  scale2 = 1.001d0
!  scale2 = 0.999d0
!  scale2 = 100.d0
!  scale2 = 0.1d0
  select case(istate2)
    case(LIQUID_STATE)
      xx2(1) = 1.d6*scale2
      xx2(2) = 1.d-6*scale2
      xx2(3) = 30.d0*scale2
!      xx2(3) = 100.d0
    case(GAS_STATE)
!      xx2(1) = 1.d4
!      xx2(2) = 0.997d4
!      xx2(3) = 15.d0
      xx2(1) = 1.d6*scale2
      xx2(2) = 0.997d6*scale2
!      xx2(2) = 0.996d6*scale2  ! to generate gradient in air mole fraction
      xx2(3) = 30.d0*scale2
!      xx2(1) = 1.d7
!      xx2(2) = 0.997d7
!      xx2(3) = 85.d0
    case(TWO_PHASE_STATE)
      xx2(1) = 1.d6*scale2
!      xx2(1) = 1.d8
      xx2(2) = 0.5d0*scale2
      xx2(3) = 30.d0*scale2
  end select    
  
  call GeneralDerivativeSetupAuxVar(istate2,xx2,pert2,general_auxvar2, &
                                    global_auxvar2,material_auxvar2, &
                                    characteristic_curves, &
                                    option)  
  
  call GeneralDerivativeAuxVar(pert,general_auxvar,global_auxvar, &
                               material_auxvar,option)

  select case(itype)
    case(ACCUMULATION)
      call GeneralDerivativeAccum(pert,general_auxvar,global_auxvar, &
                                  material_auxvar,material_parameter,option)
    case(INTERIOR_FLUX)
      call GeneralDerivativeAuxVar(pert2,general_auxvar2,global_auxvar2, &
                                   material_auxvar2,option)
      call GeneralDerivativeFlux(pert,general_auxvar,global_auxvar, &
                                 material_auxvar,characteristic_curves, &
                                 material_parameter,&
                                 pert2,general_auxvar2,global_auxvar2, &
                                 material_auxvar2,characteristic_curves, &
                                 material_parameter, &
                                 general_parameter,option)
    case(BOUNDARY_FLUX)
      ibndtype = DIRICHLET_BC     
    !  ibndtype = NEUMANN_BC
      auxvar_mapping(GENERAL_LIQUID_FLUX_INDEX) = 1
      auxvar_mapping(GENERAL_GAS_FLUX_INDEX) = 2
      auxvar_mapping(GENERAL_ENERGY_FLUX_INDEX) = 3
      auxvars(1) = -1.d-2
      auxvars(2) = -1.d-2
      auxvars(3) = -1.d-2
      ! everything downwind is XXX2, boundary is XXX
      call GeneralDerivativeAuxVar(pert2,general_auxvar2,global_auxvar2, &
                                   material_auxvar2,option)
      call GeneralDerivativeFluxBC(pert2, &
                                   ibndtype,auxvar_mapping,auxvars, &
                                   general_auxvar(ZERO_INTEGER), &
                                   global_auxvar(ZERO_INTEGER), &
                                   general_auxvar2,global_auxvar2, &
                                   material_auxvar2, &
                                   characteristic_curves, &
                                   material_parameter, &
                                   general_parameter,option)
    case(SRCSINK)
      flow_src_sink_type = VOLUMETRIC_RATE_SS
      qsrc = 0.d0
      qsrc(1) = 1.d-10
      srcsink_scale = 1.d0
      call GeneralDerivativeSrcSink(pert,qsrc,flow_src_sink_type, &
                                    general_auxvar,global_auxvar, &
                                    material_auxvar,srcsink_scale,option)
  end select
  
  ! Destroy objects
  call GeneralDerivativeDestroyAuxVar(general_auxvar,global_auxvar, &
                                      material_auxvar,option)  
  call GeneralDerivativeDestroyAuxVar(general_auxvar2,global_auxvar2, &
                                      material_auxvar2,option)
  call GeneralDerivativeDestroy(general_parameter, &
                                characteristic_curves, &
                                material_parameter,option)

end subroutine GeneralDerivativeDriver

! ************************************************************************** !

subroutine GeneralDerivativeSetup(general_parameter, &
                                  characteristic_curves, &
                                  material_parameter,option)
  use Characteristic_Curves_module
  use Material_Aux_class
  use Option_module
  
  implicit none

  type(general_parameter_type), pointer :: general_parameter
  class(characteristic_curves_type), pointer :: characteristic_curves
  type(material_parameter_type), pointer :: material_parameter
  type(option_type), pointer :: option
  
  class(sat_func_VG_type), pointer :: sf
  class(rpf_Mualem_VG_liq_type), pointer :: rpf_liq
  class(rpf_Mualem_VG_gas_type), pointer :: rpf_gas  
  
  if (.not.associated(general_parameter)) then
    allocate(general_parameter)
    allocate(general_parameter%diffusion_coefficient(2))
    general_parameter%diffusion_coefficient(1) = 1.d-9
    general_parameter%diffusion_coefficient(2) = 2.13d-5
  endif
  if (.not.associated(characteristic_curves)) then
    characteristic_curves => CharacteristicCurvesCreate()
    sf => SF_VG_Create()
    rpf_liq => RPF_Mualem_VG_Liq_Create()
    rpf_gas => RPF_Mualem_VG_Gas_Create()
    sf%m = 0.5d0
    sf%alpha = 1.d-4
    sf%Sr = 0.d0
    sf%pcmax = 1.d6
    characteristic_curves%saturation_function => sf
    rpf_liq%m = 0.5d0
    rpf_liq%Sr = 0.d0
    characteristic_curves%liq_rel_perm_function => rpf_liq
    rpf_gas%m = 0.5d0
    rpf_gas%Sr = 0.d0
    rpf_gas%Srg = 1.d-40
    characteristic_curves%gas_rel_perm_function => rpf_gas
  endif
  if (.not.associated(material_parameter)) then
    allocate(material_parameter)
    allocate(material_parameter%soil_residual_saturation(2,1))
    allocate(material_parameter%soil_heat_capacity(1))
    allocate(material_parameter%soil_thermal_conductivity(2,1))
    material_parameter%soil_residual_saturation(1,1) = rpf_liq%Sr
    material_parameter%soil_residual_saturation(2,1) = rpf_gas%Srg
    material_parameter%soil_heat_capacity(1) = 850.d0
    material_parameter%soil_thermal_conductivity(1,1) = 0.5d0
    material_parameter%soil_thermal_conductivity(2,1) = 2.d0
  endif
  
end subroutine GeneralDerivativeSetup

! ************************************************************************** !

subroutine GeneralDerivativeSetupAuxVar(istate,xx,pert,general_auxvar, &
                                        global_auxvar, &
                                        material_auxvar, &
                                        characteristic_curves,option)

  use Characteristic_Curves_module
  use Option_module
  use EOS_Gas_module
  use EOS_Water_module
  
  implicit none

  PetscInt :: istate
  PetscReal :: xx(3)
  PetscReal :: pert(3)
  type(general_auxvar_type), pointer :: general_auxvar(:)
  type(global_auxvar_type), pointer :: global_auxvar(:)
  class(material_auxvar_type), pointer :: material_auxvar(:)
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type), pointer :: option

  PetscReal :: xx_pert(3)
  PetscInt :: natural_id = 1
  PetscBool :: analytical_derivative = PETSC_TRUE
  PetscInt :: i
  
  allocate(general_auxvar(0:3))
  allocate(global_auxvar(0:3))
  allocate(material_auxvar(0:3))
  
  do i = 0, 3
    call GeneralAuxVarInit(general_auxvar(i),analytical_derivative,option)
    call GlobalAuxVarInit(global_auxvar(i),option)
    call MaterialAuxVarInit(material_auxvar(i),option)
    material_auxvar(i)%porosity_base = 0.25d0
    material_auxvar(i)%permeability = 1.d-12
    material_auxvar(i)%volume = 1.d0
    global_auxvar(i)%istate = istate
  enddo
  
  option%iflag = GENERAL_UPDATE_FOR_ACCUM
  call GeneralAuxVarCompute(xx,general_auxvar(0),global_auxvar(0), &
                            material_auxvar(0),characteristic_curves, &
                            natural_id,option)  

! default perturbation approach that does not consider global or material auxvar
#if 0
  call GeneralAuxVarPerturb(general_auxvar,global_auxvar(0), &
                            material_auxvar(0), &
                            characteristic_curves,natural_id, &
                            option)
#else
  do i = 1, 3
    option%iflag = GENERAL_UPDATE_FOR_DERIVATIVE
    xx_pert = xx
    pert(i) = perturbation_tolerance * xx(i)
    xx_pert(i) = xx(i) + pert(i)
    call GeneralAuxVarCompute(xx_pert,general_auxvar(i),global_auxvar(i), &
                              material_auxvar(i),characteristic_curves, &
                              natural_id,option)    
  enddo  
#endif

end subroutine GeneralDerivativeSetupAuxVar

! ************************************************************************** !

subroutine GeneralDerivativeSetupEOS(option)

  use Option_module
  use EOS_Gas_module
  use EOS_Water_module
  
  implicit none

  type(option_type), pointer :: option
  PetscReal :: tlow, thigh, plow, phigh
  PetscInt :: ntemp, npres  
  PetscReal :: aux(1)
  character(len=MAXWORDLENGTH) :: word
 
#if 1
  call EOSWaterSetDensity('PLANAR')
  call EOSWaterSetEnthalpy('PLANAR')
  call EOSWaterSetSteamDensity('PLANAR')
  call EOSWaterSetSteamEnthalpy('PLANAR')
#endif

  ! for ruling out density partial derivative
#if 0
  aux(1) = 996.000000000000d0
  call EOSWaterSetDensity('CONSTANT',aux)
  aux(1) = 2.27d6
  call EOSWaterSetEnthalpy('CONSTANT',aux)
#endif
#if 0
#if 1
  aux(1) = 9.899986173768605D-002
  call EOSWaterSetSteamDensity('CONSTANT',aux)  
  aux(1) = 45898000.0921749d0
  call EOSWaterSetSteamEnthalpy('CONSTANT',aux)
#endif
#if 1
!  call EOSGasSetDensityConstant(196.d0)
  call EOSGasSetDensityConstant(0.395063665868904d0)
  call EOSGasSetEnergyConstant(8841206.16464255d0)
#endif
#endif

  tlow = 1.d-1
  thigh = 350.d0
  plow = 1.d-1
  phigh = 16.d6
  ntemp = 100
  npres = 100
  word = ''
  call EOSGasTest(tlow,thigh,plow,phigh, &
                  ntemp, npres, &
                  PETSC_FALSE, PETSC_FALSE, &
                  word)  
  word = ''
  call EOSWaterTest(tlow,thigh,plow,phigh, &
                    ntemp, npres, &
                    PETSC_FALSE, PETSC_FALSE, &
                    word)  
  word = ''
  call EOSWaterSteamTest(tlow,thigh,plow,phigh, &
                         ntemp, npres, &
                         PETSC_FALSE, PETSC_FALSE, &
                         word) 
  
end subroutine GeneralDerivativeSetupEOS
  
! ************************************************************************** !

subroutine GeneralDerivativeAuxVar(pert,general_auxvar,global_auxvar, &
                                   material_auxvar,option)

  use Option_module
  
  implicit none

  PetscReal :: pert(3)
  type(general_auxvar_type) :: general_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar(0:)
  class(material_auxvar_type) :: material_auxvar(0:)
  type(option_type), pointer :: option

  PetscInt :: natural_id = 1
  PetscBool :: analytical_derivative = PETSC_TRUE
  character(len=MAXSTRINGLENGTH) :: strings(3,3)
  PetscInt :: i
  
  strings(1,1) = 'Liquid Pressure'
  strings(2,1) = 'Air Mole Fraction in Liquid'
  strings(3,1) = 'Temperature'
  strings(1,2) = 'Gas Pressure'
  strings(2,2) = 'Air Pressure'
  strings(3,2) = 'Temperature'
  strings(1,3) = 'Gas Pressure'
  strings(2,3) = 'Gas Saturation'
  strings(3,3) = 'Temperature'
  
  call GeneralPrintAuxVars(general_auxvar(0),global_auxvar(0),material_auxvar(0), &
                           natural_id,'',option)

  do i = 1, 3
    call GeneralAuxVarDiff(i,general_auxvar(0),global_auxvar(0), &
                           material_auxvar(0), &
                           general_auxvar(i),global_auxvar(i), &
                           material_auxvar(i), &
                           pert(i), &
                           strings(i,global_auxvar(i)%istate), &
                           analytical_derivative,option) 
  enddo

end subroutine GeneralDerivativeAuxVar

! ************************************************************************** !

subroutine GeneralDerivativeAccum(pert,general_auxvar,global_auxvar, &
                                  material_auxvar,material_parameter, &
                                  option)
  use Material_Aux_class
  use Option_module
  
  implicit none

  PetscReal :: pert(3)
  type(general_auxvar_type) :: general_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar(0:)
  class(material_auxvar_type) :: material_auxvar(0:)
  type(material_parameter_type) :: material_parameter
  type(option_type), pointer :: option

  PetscInt :: natural_id = 1
  PetscInt :: i
  PetscReal, parameter :: soil_heat_capacity = 850.d0
  
  PetscInt :: irow
  PetscReal :: res(3), res_pert(3,3)
  PetscReal :: jac_anal(3,3)
  PetscReal :: jac_num(3,3)
  PetscReal :: jac_dum(3,3)  
  
  call GeneralPrintAuxVars(general_auxvar(0),global_auxvar(0),material_auxvar(0), &
                           natural_id,'',option)

  call GeneralAccumulation(general_auxvar(ZERO_INTEGER), &
                           global_auxvar(ZERO_INTEGER), &
                           material_auxvar(ZERO_INTEGER), &
                           material_parameter%soil_heat_capacity(1), &
                           option, &
                           res,jac_anal,PETSC_TRUE,PETSC_FALSE)

  do i = 1, 3
    call GeneralAccumulation(general_auxvar(i), &
                             global_auxvar(i), &
                             material_auxvar(i), &
                             material_parameter%soil_heat_capacity(1), &
                             option, &
                             res_pert(:,i),jac_dum,PETSC_FALSE,PETSC_FALSE)
                           
    do irow = 1, option%nflowdof
      jac_num(irow,i) = (res_pert(irow,i)-res(irow))/pert(i)
    enddo !irow
  enddo
  
  call GeneralDiffJacobian('',jac_num,jac_anal,res,res_pert,pert, &
                           general_auxvar,option)
  
end subroutine GeneralDerivativeAccum

! ************************************************************************** !

subroutine GeneralDerivativeFlux(pert,general_auxvar,global_auxvar, &
                                 material_auxvar,characteristic_curves, &
                                 material_parameter,&
                                 pert2,general_auxvar2,global_auxvar2, &
                                 material_auxvar2,characteristic_curves2, &
                                 material_parameter2, &
                                 general_parameter,option)

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  PetscReal :: pert(3)
  type(general_auxvar_type) :: general_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar(0:)
  class(material_auxvar_type) :: material_auxvar(0:)
  class(characteristic_curves_type) :: characteristic_curves
  type(material_parameter_type) :: material_parameter
  PetscReal :: pert2(3)
  type(general_auxvar_type) :: general_auxvar2(0:)
  type(global_auxvar_type) :: global_auxvar2(0:)
  class(material_auxvar_type) :: material_auxvar2(0:)
  class(characteristic_curves_type) :: characteristic_curves2
  type(material_parameter_type) :: material_parameter2
  type(general_parameter_type) :: general_parameter
  type(option_type), pointer :: option

  PetscInt :: natural_id = 1
  PetscInt :: i
  PetscReal, parameter :: area = 1.d0
!  PetscReal, parameter :: dist(-1:3) = [0.5d0,1.d0,1.d0,0.d0,0.d0]
!  PetscReal, parameter :: dist(-1:3) = [0.5d0,1.d0,sqrt(2.d0/2.d0),0.d0,sqrt(2.d0/2.d0)]
  PetscReal, parameter :: dist(-1:3) = [0.5d0,1.d0,0.d0,0.d0,1.d0]
  
  PetscReal :: v_darcy(2)
  
  
  PetscInt :: irow
  PetscReal :: res(3)
  PetscReal :: res_pert(3,3)
  PetscReal :: jac_anal(3,3)
  PetscReal :: jac_num(3,3)
  PetscReal :: jac_dum(3,3)
  
  PetscReal :: res_pert2(3,3)
  PetscReal :: jac_anal2(3,3)
  PetscReal :: jac_num2(3,3)
  PetscReal :: jac_dum2(3,3)
  
  call GeneralPrintAuxVars(general_auxvar(0),global_auxvar(0),material_auxvar(0), &
                           natural_id,'upwind',option)
  call GeneralPrintAuxVars(general_auxvar2(0),global_auxvar2(0),material_auxvar2(0), &
                           natural_id,'downwind',option)

  call GeneralFlux(general_auxvar(ZERO_INTEGER), &
                   global_auxvar(ZERO_INTEGER), &
                   material_auxvar(ZERO_INTEGER), &
                   material_parameter%soil_thermal_conductivity(:,1), &
                   general_auxvar2(ZERO_INTEGER), &
                   global_auxvar2(ZERO_INTEGER), &
                   material_auxvar2(ZERO_INTEGER), &
                   material_parameter2%soil_thermal_conductivity(:,1), &
                   area, dist, general_parameter, &
                   option,v_darcy,res,jac_anal,jac_anal2, &
                   PETSC_TRUE,PETSC_FALSE)                           

  do i = 1, 3
    call GeneralFlux(general_auxvar(i), &
                     global_auxvar(i), &
                     material_auxvar(i), &
                     material_parameter%soil_thermal_conductivity(:,1), &
                     general_auxvar2(ZERO_INTEGER), &
                     global_auxvar2(ZERO_INTEGER), &
                     material_auxvar2(ZERO_INTEGER), &
                     material_parameter2%soil_thermal_conductivity(:,1), &
                     area, dist, general_parameter, &
                     option,v_darcy,res_pert(:,i),jac_dum,jac_dum2, &
                     PETSC_FALSE,PETSC_FALSE)                           
    do irow = 1, option%nflowdof
      jac_num(irow,i) = (res_pert(irow,i)-res(irow))/pert(i)
    enddo !irow
  enddo
  do i = 1, 3
    call GeneralFlux(general_auxvar(ZERO_INTEGER), &
                     global_auxvar(ZERO_INTEGER), &
                     material_auxvar(ZERO_INTEGER), &
                     material_parameter%soil_thermal_conductivity(:,1), &
                     general_auxvar2(i), &
                     global_auxvar2(i), &
                     material_auxvar2(i), &
                     material_parameter2%soil_thermal_conductivity(:,1), &
                     area, dist, general_parameter, &
                     option,v_darcy,res_pert2(:,i),jac_dum,jac_dum2, &
                     PETSC_FALSE,PETSC_FALSE)                           
    do irow = 1, option%nflowdof
      jac_num2(irow,i) = (res_pert2(irow,i)-res(irow))/pert2(i)
    enddo !irow
  enddo  
  
  call GeneralDiffJacobian('upwind',jac_num,jac_anal,res,res_pert,pert, &
                           general_auxvar,option)
  write(*,*) '-----------------------------------------------------------------'
  call GeneralDiffJacobian('downwind',jac_num2,jac_anal2,res,res_pert2,pert2, &
                           general_auxvar2,option)
  
end subroutine GeneralDerivativeFlux

! ************************************************************************** !

subroutine GeneralDerivativeFluxBC(pert, &
                                   ibndtype,auxvar_mapping,auxvars, &
                                   general_auxvar_bc,global_auxvar_bc, &
                                   general_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   characteristic_curves_dn, &
                                   material_parameter_dn, &
                                   general_parameter,option)

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type), pointer :: option
  PetscReal :: pert(3)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(GENERAL_MAX_INDEX)  
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(general_auxvar_type) :: general_auxvar_bc
  type(global_auxvar_type) :: global_auxvar_bc
  type(general_auxvar_type) :: general_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_dn(0:)
  class(material_auxvar_type) :: material_auxvar_dn(0:)
  class(characteristic_curves_type) :: characteristic_curves_dn
  type(material_parameter_type) :: material_parameter_dn
  type(general_parameter_type) :: general_parameter

  PetscInt :: natural_id = 1
  PetscInt :: i
  PetscReal, parameter :: area = 1.d0
  PetscReal, parameter :: dist(-1:3) = [0.5d0,1.d0,0.d0,0.d0,1.d0]
  
  PetscReal :: v_darcy(2)
  
  
  PetscInt :: irow
  PetscReal :: res(3)
  PetscReal :: res_pert(3,3)
  PetscReal :: jac_anal(3,3)
  PetscReal :: jac_num(3,3)
  PetscReal :: jac_dum(3,3)
  
  call GeneralPrintAuxVars(general_auxvar_bc,global_auxvar_bc,material_auxvar_dn(0), &
                           natural_id,'boundary',option)
  call GeneralPrintAuxVars(general_auxvar_dn(0),global_auxvar_dn(0),material_auxvar_dn(0), &
                           natural_id,'internal',option)

  call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     general_auxvar_bc,global_auxvar_bc, &
                     general_auxvar_dn(ZERO_INTEGER),global_auxvar_dn(ZERO_INTEGER), &
                     material_auxvar_dn(ZERO_INTEGER), &
                     material_parameter_dn%soil_thermal_conductivity(:,1), &
                     area,dist,general_parameter, &
                     option,v_darcy,res,jac_anal, &
                     PETSC_TRUE,PETSC_FALSE)                           
                           
  do i = 1, 3
    call GeneralBCFlux(ibndtype,auxvar_mapping,auxvars, &
                       general_auxvar_bc,global_auxvar_bc, &
                       general_auxvar_dn(i),global_auxvar_dn(i), &
                       material_auxvar_dn(i), &
                       material_parameter_dn%soil_thermal_conductivity(:,1), &
                       area,dist,general_parameter, &
                       option,v_darcy,res_pert(:,i),jac_dum, &
                       PETSC_FALSE,PETSC_FALSE)                           
    do irow = 1, option%nflowdof
      jac_num(irow,i) = (res_pert(irow,i)-res(irow))/pert(i)
    enddo !irow
  enddo
  
  write(*,*) '-----------------------------------------------------------------'
  call GeneralDiffJacobian('boundary',jac_num,jac_anal,res,res_pert,pert, &
                           general_auxvar_dn,option)
  
end subroutine GeneralDerivativeFluxBC

! ************************************************************************** !

subroutine GeneralDerivativeSrcSink(pert,qsrc,flow_src_sink_type, &
                                    general_auxvar,global_auxvar, &
                                    material_auxvar,scale,option)

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type), pointer :: option
  PetscReal :: pert(3)
  PetscReal :: qsrc(:)
  PetscInt :: flow_src_sink_type  
  type(general_auxvar_type) :: general_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar(0:)
  type(material_auxvar_type) :: material_auxvar(0:)
  PetscReal :: scale
  
  PetscReal :: ss_flow_vol_flux(option%nphase)  

  PetscInt :: natural_id = 1
  PetscInt :: i
  
  PetscInt :: irow
  PetscReal :: res(3)
  PetscReal :: res_pert(3,3)
  PetscReal :: jac_anal(3,3)
  PetscReal :: jac_num(3,3)
  PetscReal :: jac_dum(3,3)
  
  call GeneralPrintAuxVars(general_auxvar(0),global_auxvar(0),material_auxvar(0), &
                           natural_id,'srcsink',option)

  call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                      general_auxvar(ZERO_INTEGER),global_auxvar(ZERO_INTEGER), &
                      ss_flow_vol_flux, &
                      scale,res,jac_anal,PETSC_TRUE,PETSC_FALSE)                           
                           
  do i = 1, 3
    call GeneralSrcSink(option,qsrc,flow_src_sink_type, &
                        general_auxvar(i),global_auxvar(i), &
                        ss_flow_vol_flux, &
                        scale,res_pert(:,i),jac_dum,PETSC_FALSE,PETSC_FALSE)                           
    do irow = 1, option%nflowdof
      jac_num(irow,i) = (res_pert(irow,i)-res(irow))/pert(i)
    enddo !irow
  enddo
  
  write(*,*) '-----------------------------------------------------------------'
  call GeneralDiffJacobian('srcsink',jac_num,jac_anal,res,res_pert,pert, &
                           general_auxvar,option)
  
end subroutine GeneralDerivativeSrcSink

! ************************************************************************** !

subroutine GeneralDerivativeSetFlowMode(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  
  option%iflowmode = G_MODE
  option%nphase = 2
  option%liquid_phase = 1  ! liquid_pressure
  option%gas_phase = 2     ! gas_pressure

  option%air_pressure_id = 3
  option%capillary_pressure_id = 4
  option%vapor_pressure_id = 5
  option%saturation_pressure_id = 6

  option%water_id = 1
  option%air_id = 2
  option%energy_id = 3

  option%nflowdof = 3
  option%nflowspec = 2
  option%use_isothermal = PETSC_FALSE
      
end subroutine GeneralDerivativeSetFlowMode

! ************************************************************************** !

subroutine GeneralAuxVarDiff(idof,general_auxvar,global_auxvar, &
                             material_auxvar, &
                             general_auxvar_pert,global_auxvar_pert, &
                             material_auxvar_pert, &
                             pert,string,compare_analytical_derivative, &
                             option)

  use Option_module
  use General_Aux_module
  use Global_Aux_module
  use Material_Aux_class  

  implicit none
  
  type(option_type) :: option
  PetscInt :: idof
  type(general_auxvar_type) :: general_auxvar, general_auxvar_pert
  type(global_auxvar_type) :: global_auxvar, global_auxvar_pert
  class(material_auxvar_type) :: material_auxvar, material_auxvar_pert
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: compare_analytical_derivative
  PetscReal :: pert

  
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation
  PetscReal :: liquid_mass_pert, gas_mass_pert
  PetscReal :: liquid_density_pert, gas_density_pert
  PetscReal :: liquid_energy_pert, gas_energy_pert
  PetscReal :: liquid_saturation_pert, gas_saturation_pert
  
  PetscReal :: dpl 
  PetscReal :: dpg 
  PetscReal :: dpa 
  PetscReal :: dpc 
  PetscReal :: dpv 
  PetscReal :: dps 
  PetscReal :: dsatl
  PetscReal :: dsatg
  PetscReal :: ddenl  
  PetscReal :: ddeng  
  PetscReal :: ddenlkg
  PetscReal :: ddengkg
  PetscReal :: dUl 
  PetscReal :: dHl  
  PetscReal :: dUg  
  PetscReal :: dHg  
  PetscReal :: dUv
  PetscReal :: dHv  
  PetscReal :: dUa  
  PetscReal :: dHa  
  PetscReal :: dpsat  
  PetscReal :: dmobilityl  
  PetscReal :: dmobilityg  
  PetscReal :: dxmolwl
  PetscReal :: dxmolal
  PetscReal :: dxmolwg
  PetscReal :: dxmolag
  PetscReal :: denv
  PetscReal :: dena
  PetscReal :: dHc
  PetscReal :: dmug
  
  PetscReal, parameter :: uninitialized_value = -999.d0
  
  dpl = uninitialized_value
  dpg = uninitialized_value
  dpa = uninitialized_value
  dpc = uninitialized_value
  dpv = uninitialized_value
  dps = uninitialized_value
  dsatl = uninitialized_value
  dsatg = uninitialized_value
  ddenl = uninitialized_value
  ddeng = uninitialized_value
  ddenlkg = uninitialized_value
  ddengkg = uninitialized_value
  dUl = uninitialized_value
  dHl = uninitialized_value
  dUg = uninitialized_value
  dHg = uninitialized_value
  dUv = uninitialized_value
  dHv = uninitialized_value
  dUa = uninitialized_value
  dHa = uninitialized_value
  dpsat = uninitialized_value
  dmobilityl = uninitialized_value
  dmobilityg = uninitialized_value
  dxmolwl = uninitialized_value
  dxmolal = uninitialized_value
  dxmolwg = uninitialized_value
  dxmolag = uninitialized_value
  denv = uninitialized_value
  dena = uninitialized_value
  dHc = uninitialized_value
  dmug = uninitialized_value

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0
    
  if (compare_analytical_derivative) then
    select case(global_auxvar%istate)
      case(LIQUID_STATE)
        select case(idof)
          case(1) ! liquid pressure pl
            dpl = 1.d0
            dpv = 0.d0
            dpa = general_auxvar%d%Hc_p*general_auxvar%xmol(acid,lid)
            dps = 0.d0
            dpv = general_auxvar%d%pv_p
            ddenl = general_auxvar%d%denl_pl
            dUl = general_auxvar%d%Ul_pl
            dHl = general_auxvar%d%Hl_pl
            dmobilityl = general_auxvar%d%mobilityl_pl
          case(2) ! xmole air in liquid
            dxmolwl = -1.d0
            dxmolal = 1.d0
          case(3) ! temperature
            dpl = 0.d0 ! pl = pg - pc
            dpg = 0.d0 ! pg = pg
            dpa = general_auxvar%d%Hc_T*general_auxvar%xmol(acid,lid)
            dps = general_auxvar%d%psat_T
            dHc = general_auxvar%d%Hc_T            
            dsatl = 0.d0
            dsatg = 0.d0            
            ddenl = general_auxvar%d%denl_T
            ddeng = general_auxvar%d%deng_T
            ddenlkg = ddenl*fmw_comp(1)
            ddengkg = general_auxvar%d%dengkg_T
            dUl = general_auxvar%d%Ul_T
            dHl = general_auxvar%d%Hl_T
            
            dpv = general_auxvar%d%pv_T
            dmobilityl = general_auxvar%d%mobilityl_T
            dxmolwl = general_auxvar%d%xmol_T(wid,lid)
            dxmolal = general_auxvar%d%xmol_T(acid,lid)
        end select
      case(GAS_STATE)
        select case(idof)
          case(1)
            dpg = 1.d0 ! pg = pg
            dpv = general_auxvar%d%pv_p
            dps = 0.d0
            dHc = general_auxvar%d%Hc_p
            ddeng = general_auxvar%d%deng_pg
            ddengkg = general_auxvar%d%dengkg_pg
            dUg = general_auxvar%d%Ug_pg
            dHg = general_auxvar%d%Hg_pg
            
            dHv = general_auxvar%d%Hv_pg
            dUv = general_auxvar%d%Uv_pg
            dHa = general_auxvar%d%Ha_pg
            dUa = general_auxvar%d%Ua_pg
            
            dmug = general_auxvar%d%mug_pg
            dmobilityg = general_auxvar%d%mobilityg_pg
            dxmolwg = general_auxvar%d%xmol_p(wid,gid)
            dxmolag = general_auxvar%d%xmol_p(acid,gid)          
          case(2)
            dpg = 0.d0
            dpa = 1.d0
            ddeng = general_auxvar%d%deng_pa
            dpv = general_auxvar%d%pv_pa
            dUg = general_auxvar%d%Ug_pa
            dHg = general_auxvar%d%Hg_pa
            dHv = general_auxvar%d%Hv_pa
            dUv = general_auxvar%d%Uv_pa
            dHa = general_auxvar%d%Ha_pa
            dUa = general_auxvar%d%Ua_pa
            ! for gas state, derivative wrt air pressure is under lid
            dxmolwg = general_auxvar%d%xmol_p(wid,lid)
            dxmolag = general_auxvar%d%xmol_p(acid,lid)          
            dmobilityg = general_auxvar%d%mobilityg_pa
          case(3)
            dpg = 0.d0
            dpa = 0.d0
            dpv = 0.d0
            dps = general_auxvar%d%psat_T
            dHc = general_auxvar%d%Hc_T            
            ddeng = general_auxvar%d%deng_T
            ddengkg = general_auxvar%d%dengkg_T
            dUg = general_auxvar%d%Ug_T
            dHg = general_auxvar%d%Hg_T
            
            dHv = general_auxvar%d%Hv_T
            dUv = general_auxvar%d%Uv_T
            dHa = general_auxvar%d%Ha_T
            dUa = general_auxvar%d%Ua_T
            denv = general_auxvar%d%denv_T
            dena = general_auxvar%d%dena_T
            
            dmug = general_auxvar%d%mug_T
            dmobilityg = general_auxvar%d%mobilityg_T
            dxmolwg = general_auxvar%d%xmol_T(wid,gid)
            dxmolag = general_auxvar%d%xmol_T(acid,gid)          
        end select
      case(TWO_PHASE_STATE)
        select case(idof)
          case(1) ! gas pressure pg
            dpl = 1.d0 ! pl = pg - pc
            dpg = 1.d0 ! pg = pg
            dpa = 1.d0 ! pa = pg - pv
            dpv = 0.d0
            dps = 0.d0
            dsatl = 0.d0
            dsatg = 0.d0
            dHc = general_auxvar%d%Hc_p
            ddenl = general_auxvar%d%denl_pl*dpl
            ddeng = general_auxvar%d%deng_pg
            ddenlkg = ddenl*fmw_comp(1)
            ddengkg = general_auxvar%d%dengkg_pg
            dUl = general_auxvar%d%Ul_pl
            dHl = general_auxvar%d%Hl_pl
            dUg = general_auxvar%d%Ug_pg
            dHg = general_auxvar%d%Hg_pg
            
            denv = general_auxvar%d%denv_pg
            dena = general_auxvar%d%dena_pg
            dHv = general_auxvar%d%Hv_pg
            dUv = general_auxvar%d%Uv_pg
            dHa = general_auxvar%d%Ha_pg
            dUa = general_auxvar%d%Ua_pg
            
            dmug = general_auxvar%d%mug_pg
            dmobilityl = general_auxvar%d%mobilityl_pl
            dmobilityg = general_auxvar%d%mobilityg_pg
            dxmolwl = general_auxvar%d%xmol_p(wid,lid)
            dxmolal = general_auxvar%d%xmol_p(acid,lid)
            dxmolwg = general_auxvar%d%xmol_p(wid,gid)
            dxmolag = general_auxvar%d%xmol_p(acid,gid)
          case(2) ! gas saturation
            dpl = -1.d0*general_auxvar%d%pc_satg ! pl = pg - pc
            dpg = 0.d0
            dpa = 0.d0
            dpc = general_auxvar%d%pc_satg
            dpv = 0.d0 
            dps = 0.d0
            dsatl = -1.d0
            dsatg = 1.d0
            ddenl = 0.d0
            dmobilityl = -1.d0*general_auxvar%d%mobilityl_satg
            dmobilityg = general_auxvar%d%mobilityg_satg
          case(3) ! temperature
            dpl = 0.d0 ! pl = pg - pc
            dpg = 0.d0 ! pg = pg
            dpa = -1.d0*general_auxvar%d%psat_T ! pa = pg - pv
            dpv = general_auxvar%d%psat_T
            dps = general_auxvar%d%psat_T
            dHc = general_auxvar%d%Hc_T            
            dsatl = 0.d0
            dsatg = 0.d0            
            ddenl = general_auxvar%d%denl_T
            ddeng = general_auxvar%d%deng_T
            ddenlkg = ddenl*fmw_comp(1)
            ddengkg = general_auxvar%d%dengkg_T
            dUl = general_auxvar%d%Ul_T
            dHl = general_auxvar%d%Hl_T
            dUg = general_auxvar%d%Ug_T
            dHg = general_auxvar%d%Hg_T
            
            dHv = general_auxvar%d%Hv_T
            dUv = general_auxvar%d%Uv_T
            dHa = general_auxvar%d%Ha_T
            dUa = general_auxvar%d%Ua_T
            denv = general_auxvar%d%denv_T
            dena = general_auxvar%d%dena_T
            
            dmug = general_auxvar%d%mug_T
            dmobilityl = general_auxvar%d%mobilityl_T
            dmobilityg = general_auxvar%d%mobilityg_T
            dxmolwl = general_auxvar%d%xmol_T(wid,lid)
            dxmolal = general_auxvar%d%xmol_T(acid,lid)
            dxmolwg = general_auxvar%d%xmol_T(wid,gid)
            dxmolag = general_auxvar%d%xmol_T(acid,gid)
        end select
      end select
    endif

  print *, '--------------------------------------------------------'
  print *, 'Derivative with respect to ' // trim(string)
  select case(global_auxvar%istate)
    case(LIQUID_STATE)
      print *, '     Thermodynamic state: Liquid phase'
      liquid_density = general_auxvar%den(lid)
      liquid_energy = general_auxvar%U(lid)
      liquid_saturation = general_auxvar%sat(lid)
    case(GAS_STATE)
      print *, '     Thermodynamic state: Gas phase'
      gas_density = general_auxvar%den(gid)
      gas_energy = general_auxvar%U(gid)
      gas_saturation = general_auxvar%sat(gid)
    case(TWO_PHASE_STATE)
      print *, '     Thermodynamic state: Two phase'
      liquid_density = general_auxvar%den(lid)
      gas_density = general_auxvar%den(gid)
      liquid_energy = general_auxvar%U(lid)
      gas_energy = general_auxvar%U(gid)
      liquid_saturation = general_auxvar%sat(lid)
      gas_saturation = general_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*general_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*general_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 general_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*general_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*general_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              general_auxvar%effective_porosity*material_auxvar%volume
  select case(global_auxvar_pert%istate)
    case(LIQUID_STATE)
      print *, '     Thermodynamic state (pert): Liquid phase'
      liquid_density_pert = general_auxvar_pert%den(lid)
      liquid_energy_pert = general_auxvar_pert%U(lid)
      liquid_saturation_pert = general_auxvar_pert%sat(lid)
      gas_density_pert = 0.d0
      gas_energy_pert = 0.d0
      gas_saturation_pert = 0.d0
    case(GAS_STATE)
      print *, '     Thermodynamic state (pert): Gas phase'
      liquid_density_pert = 0.d0
      liquid_energy_pert = 0.d0
      liquid_saturation_pert = 0.d0
      gas_density_pert = general_auxvar_pert%den(gid)
      gas_energy_pert = general_auxvar_pert%U(gid)
      gas_saturation_pert = general_auxvar_pert%sat(gid)
    case(TWO_PHASE_STATE)
      print *, '     Thermodynamic state (pert): Two phase'
      liquid_density_pert = general_auxvar_pert%den(lid)
      gas_density_pert = general_auxvar_pert%den(gid)
      liquid_energy_pert = general_auxvar_pert%U(lid)
      gas_energy_pert = general_auxvar_pert%U(gid)
      liquid_saturation_pert = general_auxvar_pert%sat(lid)
      gas_saturation_pert = general_auxvar_pert%sat(gid)
  end select  
  liquid_mass_pert = (liquid_density_pert*general_auxvar_pert%xmol(lid,lid)* & 
                 liquid_saturation_pert+ &
                 gas_density_pert*general_auxvar_pert%xmol(lid,gid)* & 
                 gas_saturation_pert)* & 
                 general_auxvar_pert%effective_porosity*material_auxvar_pert%volume
  gas_mass_pert = (liquid_density_pert*general_auxvar_pert%xmol(gid,lid)* & 
              liquid_saturation_pert+ &
              gas_density_pert*general_auxvar_pert%xmol(gid,gid)* & 
              gas_saturation_pert)* & 
              general_auxvar_pert%effective_porosity*material_auxvar_pert%volume 
              
              
  call GeneralAuxVarPrintResult('tot liq comp mass [kmol]', &
                                (liquid_mass_pert-liquid_mass)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call GeneralAuxVarPrintResult('tot gas comp mass [kmol]', &
                                (gas_mass_pert-gas_mass)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call GeneralAuxVarPrintResult('             energy [MJ]', &
                                ((liquid_mass_pert*liquid_energy_pert + &
                                  gas_mass_pert*gas_energy_pert)- &
                                 (liquid_mass*liquid_energy + &
                                  gas_mass*gas_energy))/pert, &
                                uninitialized_value,uninitialized_value,option)
  call GeneralAuxVarPrintResult('         liquid pressure', &
                                (general_auxvar_pert%pres(lid)-general_auxvar%pres(lid))/pert, &
                                dpl,uninitialized_value,option)
  call GeneralAuxVarPrintResult('            gas pressure', &
                                (general_auxvar_pert%pres(gid)-general_auxvar%pres(gid))/pert, &
                                dpg,uninitialized_value,option)
  call GeneralAuxVarPrintResult('            air pressure', &
                                (general_auxvar_pert%pres(apid)-general_auxvar%pres(apid))/pert, &
                                dpa,uninitialized_value,option)
  call GeneralAuxVarPrintResult('      capillary pressure', &
                                (general_auxvar_pert%pres(cpid)-general_auxvar%pres(cpid))/pert, &
                                dpc,uninitialized_value,option)
  call GeneralAuxVarPrintResult('          vapor pressure', &
                                (general_auxvar_pert%pres(vpid)-general_auxvar%pres(vpid))/pert, &
                                dpv,uninitialized_value,option)
  call GeneralAuxVarPrintResult("        Henry's constant", &
                                (general_auxvar_pert%d%Hc-general_auxvar%d%Hc)/pert, &
                                dHc,uninitialized_value,option)
  call GeneralAuxVarPrintResult('     saturation pressure', &
                                (general_auxvar_pert%pres(spid)-general_auxvar%pres(spid))/pert, &
                                dps,uninitialized_value,option)
  call GeneralAuxVarPrintResult('       liquid saturation', &
                                (general_auxvar_pert%sat(lid)-general_auxvar%sat(lid))/pert, &
                                dsatl,uninitialized_value,option)
  call GeneralAuxVarPrintResult('          gas saturation', &
                                (general_auxvar_pert%sat(gid)-general_auxvar%sat(gid))/pert, &
                                dsatg,uninitialized_value,option)
  call GeneralAuxVarPrintResult('   liquid density [kmol]', &
                                (general_auxvar_pert%den(lid)-general_auxvar%den(lid))/pert, &
                                ddenl,uninitialized_value,option)
  call GeneralAuxVarPrintResult('      gas density [kmol]', &
                                (general_auxvar_pert%den(gid)-general_auxvar%den(gid))/pert, &
                                ddeng,uninitialized_value,option)
  call GeneralAuxVarPrintResult('     liquid density [kg]', &
                                (general_auxvar_pert%den_kg(lid)-general_auxvar%den_kg(lid))/pert, &
                                ddenlkg,uninitialized_value,option)
  call GeneralAuxVarPrintResult('        gas density [kg]', &
                                (general_auxvar_pert%den_kg(gid)-general_auxvar%den_kg(gid))/pert, &
                                ddengkg,uninitialized_value,option)
  call GeneralAuxVarPrintResult('         temperature [C]', &
                                (general_auxvar_pert%temp-general_auxvar%temp)/pert, &
                                uninitialized_value,uninitialized_value,option)
  call GeneralAuxVarPrintResult('      liquid H [MJ/kmol]', &
                                (general_auxvar_pert%H(lid)-general_auxvar%H(lid))/pert, &
                                dHl,uninitialized_value,option)
  call GeneralAuxVarPrintResult('         gas H [MJ/kmol]', &
                                (general_auxvar_pert%H(gid)-general_auxvar%H(gid))/pert, &
                                dHg,uninitialized_value,option)
  call GeneralAuxVarPrintResult('      liquid U [MJ/kmol]', &
                                (general_auxvar_pert%U(lid)-general_auxvar%U(lid))/pert, &
                                dUl,uninitialized_value,option)
  call GeneralAuxVarPrintResult('         gas U [MJ/kmol]', &
                                (general_auxvar_pert%U(gid)-general_auxvar%U(gid))/pert, &
                                dUg,uninitialized_value,option)
  !------------------------------
  call GeneralAuxVarPrintResult('       vapor H [MJ/kmol]', &
                                (general_auxvar_pert%d%Hv-general_auxvar%d%Hv)/pert, &
                                dHv,uninitialized_value,option)
  call GeneralAuxVarPrintResult('         air H [MJ/kmol]', &
                                (general_auxvar_pert%d%Ha-general_auxvar%d%Ha)/pert, &
                                dHa,uninitialized_value,option)
  call GeneralAuxVarPrintResult('       vapor U [MJ/kmol]', &
                                (general_auxvar_pert%d%Uv-general_auxvar%d%Uv)/pert, &
                                dUv,uninitialized_value,option)
  call GeneralAuxVarPrintResult('         air U [MJ/kmol]', &
                                (general_auxvar_pert%d%Ua-general_auxvar%d%Ua)/pert, &
                                dUa,uninitialized_value,option)
  call GeneralAuxVarPrintResult('    vapor density [kmol]', &
                                (general_auxvar_pert%d%denv-general_auxvar%d%denv)/pert, &
                                denv,uninitialized_value,option)
  call GeneralAuxVarPrintResult('      air density [kmol]', &
                                (general_auxvar_pert%d%dena-general_auxvar%d%dena)/pert, &
                                dena,uninitialized_value,option)
  !------------------------------                                
  call GeneralAuxVarPrintResult('     X (water in liquid)', &
                                (general_auxvar_pert%xmol(wid,lid)-general_auxvar%xmol(wid,lid))/pert, &
                                dxmolwl,uninitialized_value,option)
  call GeneralAuxVarPrintResult('       X (air in liquid)', &
                                (general_auxvar_pert%xmol(acid,lid)-general_auxvar%xmol(acid,lid))/pert, &
                                dxmolal,uninitialized_value,option)
  call GeneralAuxVarPrintResult('        X (water in gas)', &
                                (general_auxvar_pert%xmol(wid,gid)-general_auxvar%xmol(wid,gid))/pert, &
                                dxmolwg,uninitialized_value,option)
  call GeneralAuxVarPrintResult('          X (air in gas)', &
                                (general_auxvar_pert%xmol(acid,gid)-general_auxvar%xmol(acid,gid))/pert, &
                                dxmolag,uninitialized_value,option)
  call GeneralAuxVarPrintResult('         liquid mobility', &
                                (general_auxvar_pert%mobility(lid)-general_auxvar%mobility(lid))/pert, &
                                dmobilityl,uninitialized_value,option)
  call GeneralAuxVarPrintResult('            gas mobility', &
                                (general_auxvar_pert%mobility(gid)-general_auxvar%mobility(gid))/pert, &
                                dmobilityg,uninitialized_value,option)
  call GeneralAuxVarPrintResult('           gas viscosity', &
                                (general_auxvar_pert%d%mug-general_auxvar%d%mug)/pert, &
                                dmug,uninitialized_value,option)
  call GeneralAuxVarPrintResult('      effective porosity', &
                                (general_auxvar_pert%effective_porosity-general_auxvar%effective_porosity)/pert, &
                                uninitialized_value,uninitialized_value,option)
#if 0                                
100 format(a,2(es13.5),es16.8)  
  write(*,100) 'tot liq comp mass [kmol]: ', (liquid_mass_pert-liquid_mass)/pert
  write(*,100) 'tot gas comp mass [kmol]: ', (gas_mass_pert-gas_mass)/pert
  write(*,100) '             energy [MJ]: ', ((liquid_mass_pert*liquid_energy_pert + &
                                           gas_mass_pert*gas_energy_pert)- &
                                          (liquid_mass*liquid_energy + &
                                           gas_mass*gas_energy))/pert
  write(*,100) '         liquid pressure: ', (general_auxvar_pert%pres(lid)-general_auxvar%pres(lid))/pert,dpl
  write(*,100) '            gas pressure: ', (general_auxvar_pert%pres(gid)-general_auxvar%pres(gid))/pert,dpg
  write(*,100) '            air pressure: ', (general_auxvar_pert%pres(apid)-general_auxvar%pres(apid))/pert,dpa !,general_auxvar_pert%pres(apid)-general_auxvar%pres(apid)
  write(*,100) '      capillary pressure: ', (general_auxvar_pert%pres(cpid)-general_auxvar%pres(cpid))/pert,dpc
  write(*,100) '          vapor pressure: ', (general_auxvar_pert%pres(vpid)-general_auxvar%pres(vpid))/pert,dpv !,general_auxvar_pert%pres(vpid)-general_auxvar%pres(vpid)
  write(*,100) "        Henry's constant: ", (general_auxvar_pert%d%Hc-general_auxvar%d%Hc)/pert,dHc
  write(*,100) '     saturation pressure: ', (general_auxvar_pert%pres(spid)-general_auxvar%pres(spid))/pert,dps
  write(*,100) '       liquid saturation: ', (general_auxvar_pert%sat(lid)-general_auxvar%sat(lid))/pert,dsatl
  write(*,100) '          gas saturation: ', (general_auxvar_pert%sat(gid)-general_auxvar%sat(gid))/pert,dsatg
  write(*,100) '   liquid density [kmol]: ', (general_auxvar_pert%den(lid)-general_auxvar%den(lid))/pert,ddenl
  write(*,100) '      gas density [kmol]: ', (general_auxvar_pert%den(gid)-general_auxvar%den(gid))/pert,ddeng
  write(*,100) '     liquid density [kg]: ', (general_auxvar_pert%den_kg(lid)-general_auxvar%den_kg(lid))/pert,ddenl*fmw_comp(1)
  write(*,100) '        gas density [kg]: ', (general_auxvar_pert%den_kg(gid)-general_auxvar%den_kg(gid))/pert,ddengkg
  write(*,100) '         temperature [C]: ', (general_auxvar_pert%temp-general_auxvar%temp)/pert
  write(*,100) '      liquid H [MJ/kmol]: ', (general_auxvar_pert%H(lid)-general_auxvar%H(lid))/pert,dHl
  write(*,100) '         gas H [MJ/kmol]: ', (general_auxvar_pert%H(gid)-general_auxvar%H(gid))/pert,dHg
  write(*,100) '      liquid U [MJ/kmol]: ', (general_auxvar_pert%U(lid)-general_auxvar%U(lid))/pert,dUl
  write(*,100) '         gas U [MJ/kmol]: ', (general_auxvar_pert%U(gid)-general_auxvar%U(gid))/pert,dUg

  write(*,100) '       vapor H [MJ/kmol]: ', (general_auxvar_pert%d%Hv-general_auxvar%d%Hv)/pert,dHv
  write(*,100) '         air H [MJ/kmol]: ', (general_auxvar_pert%d%Ha-general_auxvar%d%Ha)/pert,dHa
  write(*,100) '       vapor U [MJ/kmol]: ', (general_auxvar_pert%d%Uv-general_auxvar%d%Uv)/pert,dUv

  write(*,100) '         air U [MJ/kmol]: ', (general_auxvar_pert%d%Ua-general_auxvar%d%Ua)/pert,dUa
  write(*,100) '    vapor density [kmol]: ', (general_auxvar_pert%d%denv-general_auxvar%d%denv)/pert,denv
  write(*,100) '      air density [kmol]: ', (general_auxvar_pert%d%dena-general_auxvar%d%dena)/pert,dena

  write(*,100) '     X (water in liquid): ', (general_auxvar_pert%xmol(wid,lid)-general_auxvar%xmol(wid,lid))/pert,dxmolwl
  write(*,100) '       X (air in liquid): ', (general_auxvar_pert%xmol(acid,lid)-general_auxvar%xmol(acid,lid))/pert,dxmolal
  write(*,100) '        X (water in gas): ', (general_auxvar_pert%xmol(wid,gid)-general_auxvar%xmol(wid,gid))/pert,dxmolwg
  write(*,100) '          X (air in gas): ', (general_auxvar_pert%xmol(acid,gid)-general_auxvar%xmol(acid,gid))/pert,dxmolag
  write(*,100) '         liquid mobility: ', (general_auxvar_pert%mobility(lid)-general_auxvar%mobility(lid))/pert,dmobilityl
  write(*,100) '            gas mobility: ', (general_auxvar_pert%mobility(gid)-general_auxvar%mobility(gid))/pert,dmobilityg
  write(*,100) '      effective porosity: ', (general_auxvar_pert%effective_porosity-general_auxvar%effective_porosity)/pert
#endif
  write(*,*) '--------------------------------------------------------'  
  
end subroutine GeneralAuxVarDiff

! ************************************************************************** !

subroutine GeneralAuxVarPrintResult(string,numerical,analytical, &
                                    uninitialized_value,option)

  use Option_module
  
  implicit none
  
  character(len=*) :: string
  PetscReal :: numerical
  PetscReal :: analytical
  PetscReal :: uninitialized_value
  type(option_type) :: option
  
  character(len=8) :: word
  PetscReal :: tempreal
  PetscReal, parameter :: tol = 1.d-5
  PetscInt :: precision
          
100 format(a24,': ',2(es13.5),2x,i2,x,a8,x,es16.8)

  precision = GeneralDigitsOfAccuracy(numerical,analytical)
  word = ''
  if (dabs(analytical-uninitialized_value) > 1.d-20) then
    if (dabs(analytical) > 0.d0) then
      tempreal = dabs((numerical-analytical)/analytical)
      if (tempreal < tol) then
        word = ' PASS'
      else
        word = '-FAIL-'
        if (tempreal < 1.d1*tol) word = trim(word) // ' *'
      endif
    else
      if (dabs(numerical) > 1.d-20) then
        word = '-FAIL-'
      else
        word = ' PASS'
      endif
    endif
    write(*,100) trim(string), numerical, analytical, precision, word
  else
    write(*,100) trim(string), numerical
  endif
  

end subroutine GeneralAuxVarPrintResult

! ************************************************************************** !

subroutine GeneralDiffJacobian(string,numerical_jacobian,analytical_jacobian, &
                               residual,residual_pert,perturbation, &
                               general_auxvar,option)

  use Option_module
  
  implicit none
  
  character(len=*) :: string
  PetscReal :: numerical_jacobian(3,3)
  PetscReal :: analytical_jacobian(3,3)
  PetscReal :: residual(3)
  PetscReal :: residual_pert(3,3)
  PetscReal :: perturbation(3)
  type(general_auxvar_type) :: general_auxvar(0:)
  type(option_type) :: option
  
  PetscInt :: irow, icol
  
100 format(2i2,2es13.5,i3,es16.8)

  if (len_trim(string) > 1) then
    write(*,'(x,a)') string
  endif
  write(*,'(" Perturbation tolerance: ",es12.4)') perturbation_tolerance
  write(*,'(" r c  numerical    analytical   digits of accuracy")')
  do icol = 1, 3
    do irow = 1, 3
      write(*,100) irow, icol, numerical_jacobian(irow,icol), &
                   analytical_jacobian(irow,icol), &
                   GeneralDigitsOfAccuracy(numerical_jacobian(irow,icol), &
                                           analytical_jacobian(irow,icol))
    enddo
  enddo

200 format(2es20.12)
300 format(a24,10es20.12)

#if 0
  do icol = 1, 3
    write(*,'(/," dof = ",i1,"  perturbation = ",es13.5)') icol, perturbation(icol)
!    write(*,300) 'density', general_auxvar(icol)%den(:), general_auxvar(0)%den(:)
!    write(*,300) 'energy', general_auxvar(icol)%U(:), general_auxvar(0)%U(:)
    write(*,'("  residual_pert       residual")')
    do irow = 1, 3
      write(*,200) residual_pert(irow,icol), residual(irow)
    enddo
  enddo
#endif  

  
end subroutine GeneralDiffJacobian

! ************************************************************************** !

function GeneralDigitsOfAccuracy(num1,num2)

  implicit none
  
  PetscReal :: num1
  PetscReal :: num2
  
  PetscInt :: GeneralDigitsOfAccuracy
  
  PetscReal :: tempreal
  PetscReal :: relative_difference
  
  GeneralDigitsOfAccuracy = 0
  if (dabs(num1) > 0.d0 .and. dabs(num2) > 0.d0) then
    relative_difference = dabs((num1-num2)/num2)
    if (relative_difference < 1.d-17) then
      ! accuracy is beyond double precision
      GeneralDigitsOfAccuracy = 99
    else
      tempreal = 1.d0 / relative_difference
      do
        if (tempreal < 10.d0) exit
        tempreal = tempreal / 10.d0
        GeneralDigitsOfAccuracy = GeneralDigitsOfAccuracy + 1
      enddo
    endif
  else
    ! change this value if you want to report something difference for
    ! double zeros.
    GeneralDigitsOfAccuracy = 0
  endif
    
end function GeneralDigitsOfAccuracy

! ************************************************************************** !

subroutine GeneralDerivativeDestroyAuxVar(general_auxvar,global_auxvar, &
                                          material_auxvar,option)

  use Option_module
  
  implicit none

  type(general_auxvar_type), pointer :: general_auxvar(:)
  type(global_auxvar_type), pointer :: global_auxvar(:)
  class(material_auxvar_type), pointer :: material_auxvar(:)
  type(option_type), pointer :: option
  
  PetscInt :: i

  do i = 0, 3
    call GeneralAuxVarStrip(general_auxvar(i))
    call GlobalAuxVarStrip(global_auxvar(i))
    call MaterialAuxVarStrip(material_auxvar(i))
  enddo
  deallocate(general_auxvar)
  deallocate(global_auxvar)
  deallocate(material_auxvar)

end subroutine GeneralDerivativeDestroyAuxVar

! ************************************************************************** !

subroutine GeneralDerivativeDestroy(general_parameter, &
                                  characteristic_curves, &
                                  material_parameter,option)
  use General_Aux_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use Option_module
  
  implicit none

  type(general_parameter_type), pointer :: general_parameter
  class(characteristic_curves_type), pointer :: characteristic_curves
  type(material_parameter_type), pointer :: material_parameter
  type(option_type), pointer :: option
  
  if (associated(general_parameter)) then
    deallocate(general_parameter%diffusion_coefficient)
    nullify(general_parameter%diffusion_coefficient)
    deallocate(general_parameter)
    nullify(general_parameter)
  endif
  call CharacteristicCurvesDestroy(characteristic_curves)
  if (associated(material_parameter)) then
    if (associated(material_parameter%soil_residual_saturation)) &
      deallocate(material_parameter%soil_residual_saturation)
    nullify(material_parameter%soil_residual_saturation)
    if (associated(material_parameter%soil_heat_capacity)) &
      deallocate(material_parameter%soil_heat_capacity)
    nullify(material_parameter%soil_heat_capacity)
    if (associated(material_parameter%soil_thermal_conductivity)) &
      deallocate(material_parameter%soil_thermal_conductivity)
    nullify(material_parameter%soil_thermal_conductivity)
    deallocate(material_parameter)
    nullify(material_parameter)    
  endif

end subroutine GeneralDerivativeDestroy

end module General_Derivative_module
