module PM_TOilIms_Aux_module

  use PFLOTRAN_Constants_module

  use PM_Base_Aux_module
  use AuxVars_TOilIms_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  !BEGINNING-Paramters to move into toil_ims_paramters 
  PetscReal, public :: toil_ims_window_epsilon = 1.d-4
  PetscReal, public :: toil_ims_fmw_comp(2) = & ! initialised after EOSread
                        [UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE]
  PetscReal, public :: toil_ims_max_pressure_change = 5.d4
  PetscInt, public :: toil_ims_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: toil_ims_damping_factor = 0.6d0
  PetscReal, public :: toil_ims_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: toil_ims_itol_scaled_res = 1.d-5
  PetscReal, public :: toil_ims_tgh2_itol_scld_res_e1(3) = 1.d-5
  PetscReal, public :: toil_ims_tgh2_itol_scld_res_e2 = 1.d0
  PetscBool, public :: toil_ims_tough2_conv_criteria = PETSC_FALSE
  PetscInt, public :: toil_ims_debug_cell_id = UNINITIALIZED_INTEGER
  !END-Paramters to move into toil_ims_paramters

  ! if needed must be specifi to toil_ims (e.g. TOIL_IMS_PREV_TS)
  !PetscInt, parameter, public :: TOIL_IMS_PREV_TS = 1
  !PetscInt, parameter, public :: TOIL_IMS_PREV_IT = 2

  ! Primary DOF indices 
  PetscInt, parameter, public :: TOIL_IMS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TOIL_IMS_SATURATION_DOF = 2
  PetscInt, parameter, public :: TOIL_IMS_ENERGY_DOF = 3

  ! Equation indices   
  PetscInt, parameter, public :: TOIL_IMS_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: TOIL_IMS_OIL_EQUATION_INDEX = 2
  PetscInt, parameter, public :: TOIL_IMS_ENERGY_EQUATION_INDEX = 3

  ! Indices used to map aux_real for condition values 
  PetscInt, parameter, public :: TOIL_IMS_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: TOIL_IMS_OIL_SATURATION_INDEX = 2
  PetscInt, parameter, public :: TOIL_IMS_TEMPERATURE_INDEX = 3
  PetscInt, parameter, public :: TOIL_IMS_LIQUID_FLUX_INDEX = 4
  PetscInt, parameter, public :: TOIL_IMS_OIL_FLUX_INDEX = 5
  PetscInt, parameter, public :: TOIL_IMS_ENERGY_FLUX_INDEX = 6
  PetscInt, parameter, public :: TOIL_IMS_LIQ_CONDUCTANCE_INDEX = 7
  PetscInt, parameter, public :: TOIL_IMS_OIL_CONDUCTANCE_INDEX = 8
  PetscInt, parameter, public :: TOIL_IMS_MAX_INDEX = 9

  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: TOIL_IMS_UPDATE_FOR_BOUNDARY = 2

  PetscReal, parameter, public :: toil_ims_pressure_scale = 1.d0
  
  ! these variables, which are global to general, can be modified
  PetscInt, public :: toil_ims_dof_to_primary_vars(3) ! only one 2ph state 
  PetscBool, public :: toil_ims_isothermal = PETSC_FALSE

  !phase mapping:
  !While LIQUID_PHASE = 1 alway, OIL_PHASE can be:
  ! 2 for brine/oil
  ! 3 for brine/gas/oil (can be 2 again, but then GAS_PHASE requires
  !                      a different index). Need mapping in any case
  PetscInt, parameter, public :: TOIL_IMS_OIL_PHASE = 2

  ! it might be required for thermal diffusion terms and tough conv criteria
  ! consider to pput here all 
  type, public :: toil_ims_parameter_type
  !  PetscReal, public :: window_epsilon 
  !  PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
  !  PetscReal :: newton_inf_scaled_res_tol
  !  PetscBool :: check_post_converged
  end type toil_ims_parameter_type

  !if required, could add other classes in the middle:
  ! pm_base_aux's doughters & pm_toil_ims_aux'parents 
  type, public, extends(pm_base_aux_type) :: pm_toil_ims_aux_type 
    type(toil_ims_parameter_type), pointer :: parameter
    !type(auxvar_toil_ims_type), pointer :: auxvars(:,:)
    !type(auxvar_toil_ims_type), pointer :: auxvars_bc(:)
    !type(auxvar_toil_ims_type), pointer :: auxvars_ss(:)
    class(auxvar_toil_ims_type), pointer :: auxvars(:,:)
    class(auxvar_toil_ims_type), pointer :: auxvars_bc(:)
    class(auxvar_toil_ims_type), pointer :: auxvars_ss(:)
  contains
    !add bound-procedure
    
    procedure, public :: Init => InitTOilImsAuxVars
    !procedure, public :: Perturb => PerturbTOilIms
  end type pm_toil_ims_aux_type

  interface TOilImsAuxVarStrip
    module procedure TOilImsAuxVarArray1Strip
    module procedure TOilImsAuxVarArray2Strip
  end interface TOilImsAuxVarStrip
 
  public :: TOilImsAuxCreate, TOilImsAuxVarCompute, &
            TOilImsAuxVarPerturb, TOilImsAuxDestroy, &
            TOilImsAuxVarStrip
  !          AuxDestroy
  ! create only the base part of the aux_vars   

contains

! ************************************************************************** !

function TOilImsAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/20/25
  ! 

  use Option_module
  use EOS_Oil_module

  implicit none

  type(option_type) :: option
    
  class(pm_toil_ims_aux_type), pointer :: TOilImsAuxCreate  

  class(pm_toil_ims_aux_type), pointer :: aux

  ! there is no variable switch, but this map can be used 
  ! to change the primary variable set 
  toil_ims_dof_to_primary_vars(1:3) = &
       [TOIL_IMS_PRESSURE_INDEX, TOIL_IMS_OIL_SATURATION_INDEX, &
        TOIL_IMS_TEMPERATURE_INDEX]

  toil_ims_fmw_comp(LIQUID_PHASE) = FMWH2O
  !toil_ims_fmw_comp(TOIL_IMS_OIL_PHASE) = fmw_oil ! this defines the dead oil
  toil_ims_fmw_comp(TOIL_IMS_OIL_PHASE) = EOSOilGetFMW() !defines the dead oil

    !allocate here to define this is a pm_toil_ims_aux_type
    allocate(aux) 

    call PMBaseAuxInit(aux) 

    nullify(aux%auxvars)
    nullify(aux%auxvars_bc)
    nullify(aux%auxvars_ss)

   !parameter not needed for now toil_ims
   allocate(aux%parameter)


  TOilImsAuxCreate => aux
  
end function TOilImsAuxCreate

! ************************************************************************** !

subroutine InitTOilImsAuxVars(this,grid,num_bc_connection, &
                              num_ss_connection,option)
  ! 
  ! Initialize pm_aux  
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 5/27/16
  ! 

  use Option_module
  use Grid_module 

  implicit none

  class(pm_toil_ims_aux_type) :: this
  PetscInt :: num_bc_connection
  PetscInt :: num_ss_connection  
  type(grid_type) :: grid
  type(option_type) :: option

  PetscInt :: ghosted_id, iconn, local_id
  PetscInt :: idof 
 
  allocate(this%auxvars(0:option%nflowdof,grid%ngmax)) 
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      !call TOilImsAuxVarInit(toil_auxvars(idof,ghosted_id),option)
      call this%auxvars(idof,ghosted_id)%Init(option)
    enddo
  enddo

  this%num_aux = grid%ngmax

  if (num_bc_connection > 0) then
    allocate(this%auxvars_bc(num_bc_connection))
    do iconn = 1, num_bc_connection
      call this%auxvars_bc(iconn)%Init(option)
    enddo
  endif
  this%num_aux_bc = num_bc_connection

  if (num_ss_connection > 0) then
    allocate(this%auxvars_ss(num_ss_connection))
    do iconn = 1, num_ss_connection
      call this%auxvars_ss(iconn)%Init(option)
    enddo
  endif
  this%num_aux_ss = num_ss_connection

  call PMBaseAuxSetup(this,grid,option)

end subroutine InitTOilImsAuxVars


! ************************************************************************** !
! this could (actually should) be moved to auxvar_toil_ims as done for Init

subroutine TOilImsAuxVarCompute(x,toil_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/21/15
  ! 

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Oil_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  !type(toil_ims_auxvar_type) :: toil_auxvar
  class(auxvar_toil_ims_type) :: toil_auxvar
  type(global_auxvar_type) :: global_auxvar ! passing this for salt conc.
                                            ! not currenty used  
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id !only for debugging/print out - currently not used 

  PetscInt :: lid, oid, cpid
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: krl, visl, dkrl_Se
  PetscReal :: kro, viso, dkro_Se
  PetscReal :: dummy
  PetscReal :: Uoil_J_kg, Hoil_J_kg
  PetscErrorCode :: ierr

  ! from SubsurfaceSetFlowMode
  ! option%liquid_phase = 1           ! liquid_pressure
  ! option%oil_phase = 2              ! oil_pressure
  ! option%capillary_pressure_id = 3  ! capillary pressure
 
  lid = option%liquid_phase
  oid = option%oil_phase
  cpid = option%capillary_pressure_id

  ! do we need this initializations? All values are overwritten below 
  toil_auxvar%H = 0.d0
  toil_auxvar%U = 0.d0
  toil_auxvar%pres = 0.d0
  toil_auxvar%sat = 0.d0
  toil_auxvar%den = 0.d0
  toil_auxvar%den_kg = 0.d0
  toil_auxvar%effective_porosity = 0.d0

  toil_auxvar%mobility = 0.d0

  !assing auxvars given by the solution variables
  toil_auxvar%pres(oid) = x(TOIL_IMS_PRESSURE_DOF)
  toil_auxvar%sat(oid) = x(TOIL_IMS_SATURATION_DOF)
  toil_auxvar%temp = x(TOIL_IMS_ENERGY_DOF)

  toil_auxvar%sat(lid) = 1.d0 - toil_auxvar%sat(oid)
      
  call characteristic_curves%saturation_function% &
             CapillaryPressure(toil_auxvar%sat(lid),toil_auxvar%pres(cpid), &
                               option)                             

  ! Testing zero capillary pressure
  !toil_auxvar%pres(cpid) = 0.d0

  toil_auxvar%pres(lid) = toil_auxvar%pres(oid) - &
                          toil_auxvar%pres(cpid)

  cell_pressure = max(toil_auxvar%pres(lid),toil_auxvar%pres(oid))


  ! calculate effective porosity as a function of pressure
  if (option%iflag /= TOIL_IMS_UPDATE_FOR_BOUNDARY) then
    toil_auxvar%effective_porosity = material_auxvar%porosity_base

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                toil_auxvar%effective_porosity,dummy)
    endif
    if (option%iflag /= TOIL_IMS_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = toil_auxvar%effective_porosity
    endif
  endif

  ! UPDATE THERMODYNAMIC PROPERTIES FOR BOTH PHASES!!!

  ! Liquid phase thermodynamic properties
  ! must use cell_pressure as the pressure, not %pres(lid)
  call EOSWaterDensity(toil_auxvar%temp,cell_pressure, &
                       toil_auxvar%den_kg(lid),toil_auxvar%den(lid),ierr)
  call EOSWaterEnthalpy(toil_auxvar%temp,cell_pressure,toil_auxvar%H(lid),ierr)
  toil_auxvar%H(lid) = toil_auxvar%H(lid) * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp
  toil_auxvar%U(lid) = toil_auxvar%H(lid) - &
                       ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                       (cell_pressure / toil_auxvar%den(lid) * &
                        1.d-6)

  ! ADD HERE BRINE dependency. Two options (see mphase)
  ! - salinity constant in space and time (passed in option%option%m_nacl)
  ! - salt can be trasnported by RT (sequential coupling) and passed 
  !   and passed with global_auxvar%m_nacl 
  !  ! Assign salinity 
  !  m_na=option%m_nacl; m_cl=m_na; m_nacl=m_na 
  !  if (option%ntrandof > 0) then
  !    m_na = global_auxvar%m_nacl(1)
  !    m_cl = global_auxvar%m_nacl(2)
  !    m_nacl = m_na
  !    if (m_cl > m_na) m_nacl = m_cl
  !  endif    
  !
  !  ! calculate density for pure water
  !  call EOSWaterDensityEnthalpy(t,pw,dw_kg,dw_mol,hw,ierr)
  !  !..................
  !  xm_nacl = m_nacl*FMWNACL
  !  xm_nacl = xm_nacl/(1.D3 + xm_nacl)
  !  ! corrects water densit previously calculated as pure water
  !  call EOSWaterDensityNaCl(t,p,xm_nacl,dw_kg)  
  !  ! water viscosity dependence on salt concetration, but no derivatives
  !  !  call EOSWaterViscosityNaCl(t,p,xm_nacl,visl)
  !  call EOSWaterViscosity(t,pw,sat_pressure,0.d0,visl,dvdt,dvdp,dvdps,ierr)

  call EOSOilDensityEnergy(toil_auxvar%temp,toil_auxvar%pres(oid),&
                           toil_auxvar%den(oid),toil_auxvar%H(oid), &
                           toil_auxvar%U(oid),ierr)

  !toil_auxvar%den_kg(oid) = toil_auxvar%den(oid) * fmw_oil 
  toil_auxvar%den_kg(oid) = toil_auxvar%den(oid) * EOSOilGetFMW()

  toil_auxvar%H(oid) = toil_auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  toil_auxvar%U(oid) = toil_auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! compute water mobility (rel. perm / viscostiy)
  call characteristic_curves%liq_rel_perm_function% &
         RelativePermeability(toil_auxvar%sat(lid),krl,dkrl_Se,option)
                            
  call EOSWaterSaturationPressure(toil_auxvar%temp, wat_sat_pres,ierr)                   

  ! use cell_pressure; cell_pressure - psat calculated internally
  call EOSWaterViscosity(toil_auxvar%temp,cell_pressure,wat_sat_pres,visl,ierr)

  toil_auxvar%mobility(lid) = krl/visl


  ! compute oil mobility (rel. perm / viscostiy)
  call characteristic_curves%oil_rel_perm_function% &
         RelativePermeability(toil_auxvar%sat(lid),kro,dkro_Se,option)

  call EOSOilViscosity(toil_auxvar%temp,toil_auxvar%pres(oid), &
                       toil_auxvar%den(oid), viso, ierr)

  toil_auxvar%mobility(oid) = kro/viso


end subroutine TOilImsAuxVarCompute

! ************************************************************************** !

#if 0
subroutine PerturbTOilIms(this,ghosted_id,global_auxvar,material_auxvar, &
                   characteristic_curves, natural_id,option)
  !
  ! Calculates auxiliary variables for perturbed system
  ! 
  ! Author: Paolo Orsini
  ! Date: 5/28/16
  ! 
  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  class(pm_toil_ims_aux_type) :: this
  !type(grid_type) :: grid
  PetscInt :: ghosted_id
  type(option_type) :: option
  PetscInt :: natural_id !only for debugging/print out - currently not used 
  !type(toil_ims_auxvar_type) :: toil_auxvar(0:)
  !class(auxvar_toil_ims_type) :: toil_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: tempreal

  ! IMS uses 1.d-8 for Press and 1.d-5 for Sat ad Temp
  ! MPHASE uses 1.d-8 for Press, Sat and Temp  
  ! GENERAL uses 1.d-8 for Press, Sat and Temp
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  !PetscReal, parameter :: perturbation_tolerance = 1.d-5

  !PetscReal, parameter :: min_mole_fraction_pert = 1.d-12
  ! From General (min perturb = 1.d-10)
  PetscReal, parameter :: min_perturbation = 1.d-10
  !PetscReal, parameter :: min_perturbation = 1.d-7
  PetscInt :: idof
  
  !in here we would need "select type" for auxvars so that
  ! the compiler could see the auxvar_toil_ims_type components
  ! no idea!!
  x(TOIL_IMS_PRESSURE_DOF) = &
     this%auxvars(ZERO_INTEGER,ghosted_id)%pres(option%oil_phase)
  x(TOIL_IMS_SATURATION_DOF) = &
     this%auxvars(ZERO_INTEGER,ghosted_id)%sat(option%oil_phase)

  x(TOIL_IMS_ENERGY_DOF) = this%auxvars(ZERO_INTEGER,ghosted_id)%temp

  pert(TOIL_IMS_ENERGY_DOF) = &
           perturbation_tolerance*x(TOIL_IMS_ENERGY_DOF)+min_perturbation

  pert(TOIL_IMS_PRESSURE_DOF) = &
         perturbation_tolerance*x(TOIL_IMS_PRESSURE_DOF)+min_perturbation

  ! always perturb toward 0.5 (0.9 is used in ims and mphase
  if (x(TOIL_IMS_SATURATION_DOF) > 0.5d0) then 
      pert(TOIL_IMS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
  else
     pert(TOIL_IMS_SATURATION_DOF) = perturbation_tolerance
  endif

! below constrains from ims and mphase - might not be needed is post-solve checks
!  if (pert(TOIL_IMS_SATURATION_DOF) < 1D-12 .and. pert(TOIL_IMS_SATURATION_DOF) >= 0.D0) &
!     pert(TOIL_IMS_SATURATION_DOF) = 1D-12
!  if (pert(TOIL_IMS_SATURATION_DOF) > -1D-12 .and. pert(TOIL_IMS_SATURATION_DOF) < 0.D0) &
!     pert(TOIL_IMS_SATURATION_DOF) = -1D-12
!       
!  if ( ( pert(TOIL_IMS_SATURATION_DOF)+ x(TOIL_IMS_SATURATION_DOF) )>1.D0) then
!    pert(TOIL_IMS_SATURATION_DOF) = (1.D0-pert(TOIL_IMS_SATURATION_DOF))*1D-6
!  endif
!  if (( pert(TOIL_IMS_SATURATION_DOF)+x(TOIL_IMS_SATURATION_DOF) ) < 0.D0)then
!    pert(TOIL_IMS_SATURATION_DOF) = x(TOIL_IMS_SATURATION_DOF)*1D-6
!  endif


  ! TOIL_IMS_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = TOIL_IMS_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    toil_auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    !call TOilImsAuxVarCompute(x_pert,toil_auxvar(idof),global_auxvar, &
    !                          material_auxvar, &
    !                          characteristic_curves,natural_id,option)
    call TOilImsAuxVarCompute(x_pert,this%auxvars(idof,ghosted_id),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id,option)
  enddo

  toil_auxvar(TOIL_IMS_PRESSURE_DOF)%pert = &
     toil_auxvar(TOIL_IMS_PRESSURE_DOF)%pert / toil_ims_pressure_scale


end subroutine PerturbTOilIms
! ************************************************************************** !
#endif

! ************************************************************************** !
! this could (actually should) be moved to auxvar_toil_ims as done for Init
! see example abive
subroutine TOilImsAuxVarPerturb(toil_auxvar,global_auxvar, &
                                material_auxvar, &
                                characteristic_curves,natural_id, &
                                option)
  ! 
  ! Calculates auxiliary variables for perturbed system
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/05/15
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id !only for debugging/print out - currently not used 
  !type(toil_ims_auxvar_type) :: toil_auxvar(0:)
  !class(auxvar_toil_ims_type) :: toil_auxvar(0:)
  type(auxvar_toil_ims_type) :: toil_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
     
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: tempreal

  ! IMS uses 1.d-8 for Press and 1.d-5 for Sat ad Temp
  ! MPHASE uses 1.d-8 for Press, Sat and Temp  
  ! GENERAL uses 1.d-8 for Press, Sat and Temp
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: sat_perturb_tolerance = 1.d-8
  !PetscReal, parameter :: perturbation_tolerance = 1.d-5

  !PetscReal, parameter :: min_mole_fraction_pert = 1.d-12
  ! From General (min perturb = 1.d-10)
  PetscReal, parameter :: min_perturbation = 1.d-10
  !PetscReal, parameter :: min_perturbation = 1.d-7
  PetscInt :: idof

  x(TOIL_IMS_PRESSURE_DOF) = &
     toil_auxvar(ZERO_INTEGER)%pres(option%oil_phase)
  x(TOIL_IMS_SATURATION_DOF) = &
     toil_auxvar(ZERO_INTEGER)%sat(option%oil_phase)

  x(TOIL_IMS_ENERGY_DOF) = toil_auxvar(ZERO_INTEGER)%temp

  pert(TOIL_IMS_ENERGY_DOF) = &
           perturbation_tolerance*x(TOIL_IMS_ENERGY_DOF)+min_perturbation

  ! below how the water pressure is perturbed in general
  ! pert(GENERAL_LIQUID_PRESSURE_DOF) = &
  !      perturbation_tolerance*x(GENERAL_LIQUID_PRESSURE_DOF) + &
  !      min_perturbation
  pert(TOIL_IMS_PRESSURE_DOF) = &
         perturbation_tolerance*x(TOIL_IMS_PRESSURE_DOF)+min_perturbation

  ! always perturb toward 0.5 (0.9 is used in ims and mphase
  if (x(TOIL_IMS_SATURATION_DOF) > 0.5d0) then 
      !pert(TOIL_IMS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      pert(TOIL_IMS_SATURATION_DOF) = -1.d0 * sat_perturb_tolerance
  else
     !pert(TOIL_IMS_SATURATION_DOF) = perturbation_tolerance
     pert(TOIL_IMS_SATURATION_DOF) = sat_perturb_tolerance
  endif

! below constrains from ims and mphase - might not be needed with post-solve checks
!  if (pert(TOIL_IMS_SATURATION_DOF) < 1D-12 .and. pert(TOIL_IMS_SATURATION_DOF) >= 0.D0) &
!     pert(TOIL_IMS_SATURATION_DOF) = 1D-12
!  if (pert(TOIL_IMS_SATURATION_DOF) > -1D-12 .and. pert(TOIL_IMS_SATURATION_DOF) < 0.D0) &
!     pert(TOIL_IMS_SATURATION_DOF) = -1D-12
!       
!  if ( ( pert(TOIL_IMS_SATURATION_DOF)+ x(TOIL_IMS_SATURATION_DOF) )>1.D0) then
!    pert(TOIL_IMS_SATURATION_DOF) = (1.D0-pert(TOIL_IMS_SATURATION_DOF))*1D-6
!  endif
!  if (( pert(TOIL_IMS_SATURATION_DOF)+x(TOIL_IMS_SATURATION_DOF) ) < 0.D0)then
!    pert(TOIL_IMS_SATURATION_DOF) = x(TOIL_IMS_SATURATION_DOF)*1D-6
!  endif


  ! TOIL_IMS_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = TOIL_IMS_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    toil_auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call TOilImsAuxVarCompute(x_pert,toil_auxvar(idof),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id,option)

  enddo

  toil_auxvar(TOIL_IMS_PRESSURE_DOF)%pert = &
     toil_auxvar(TOIL_IMS_PRESSURE_DOF)%pert / toil_ims_pressure_scale

  !write(*,"('within TOilImsAuxVarPerturb perturb = ',(4(e10.4,1x)))"), &
  !         toil_auxvar(0:3)%pert
 
end subroutine TOilImsAuxVarPerturb


! ************************************************************************** !

subroutine TOilImsAuxDestroy(aux)
  ! 
  ! Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/30/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  !type(TOil_ims_type), pointer :: aux
  class(pm_toil_ims_aux_type), pointer :: aux
  !class(pm_toil_ims_aux_type) :: aux
  PetscInt :: iaux, idof
  
  if (.not.associated(aux)) return
  
!debug testing 
#if 0  
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, size(aux%auxvars_bc(:))
      !print *, "i am in TOilImsAuxVarArray1Destroy iaux loop", iaux
      !call TOilImsAuxVarStrip(auxvars(iaux))
      call aux%auxvars_bc(iaux)%Strip()
    enddo  
    !print *, "i am in TOilImsAuxVarArray1Destroy after auxvars loop"
    deallocate(aux%auxvars_bc)
    !print *, "after deallocation 1D auxvars"
  endif
  nullify(aux%auxvars_bc)  
  !print *, "after 1D nullify auxvars"


  if (associated(aux%auxvars_ss)) then
    do iaux = 1, size(aux%auxvars_ss(:))
      !print *, "i am in TOilImsAuxVarArray1Destroy iaux loop", iaux
      !call TOilImsAuxVarStrip(auxvars(iaux))
      call aux%auxvars_ss(iaux)%Strip()
    enddo  
    !print *, "i am in TOilImsAuxVarArray1Destroy after auxvars loop"
    deallocate(aux%auxvars_ss)
    !print *, "after deallocation 1D auxvars"
  endif
  nullify(aux%auxvars_ss)  
  !print *, "after 1D nullify auxvars"


  if (associated(aux%auxvars)) then
    do iaux = 1, size(aux%auxvars,2)
      do idof = 1, size(aux%auxvars,1)
        !print *, "i am before TOilImsAuxVarStrip"
        !call TOilImsAuxVarStrip(aux%auxvars(idof-1,iaux))
        call aux%auxvars(idof-1,iaux)%Strip()
      enddo
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)  
#endif

  !debug printing  
  !print *, "den oil 01 int = ", aux%auxvars(0,1)%den(2)
  !print *, "den oil 11 int = ", aux%auxvars(1,1)%den(2)
  !print *, "den oil 21 int = ", aux%auxvars(2,1)%den(2)
  !print *, "den oil 31 int = ", aux%auxvars(3,1)%den(2)
  if (associated(aux%auxvars) ) then
    call TOilImsAuxVarStrip(aux%auxvars)
    deallocate(aux%auxvars)
  end if 
  nullify(aux%auxvars) 

  !print *, "den oil bc = ", aux%auxvars_bc(1)%den(2)
  if (associated(aux%auxvars_bc) ) then
    call TOilImsAuxVarStrip(aux%auxvars_bc)
    deallocate(aux%auxvars_bc)
  end if
  nullify(aux%auxvars_bc)

  !print *, "den oil bc = ", aux%auxvars_ss(1)%den(2)
  if ( associated(aux%auxvars_ss) ) then
    call TOilImsAuxVarStrip(aux%auxvars_ss)
    deallocate(aux%auxvars_ss)
  end if
  nullify(aux%auxvars_ss)

  call PMBaseAuxStrip(aux)

  if (associated(aux%parameter)) then
    deallocate(aux%parameter) 
  end if
  nullify(aux%parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine TOilImsAuxDestroy

! ************************************************************************** !

subroutine  TOilImsAuxVarArray1Strip(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes 
  ! using class(*) (unlimited polymorphic)
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/20/15
  ! 

  use AuxVars_TOilIms_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !here we can pass by pointer, we could destroy the array within the routine
  !but we don't to be consistent with TOilImsAuxVarArray2Strip 
  !class(auxvar_toil_ims_type), pointer :: auxvars(:)
  type(auxvar_toil_ims_type) :: auxvars(:)

  PetscInt :: iaux

  !print *, "den oil bc/ss pass = ", auxvars(1)%den(2)
  
  do iaux = 1, size(auxvars)
    call auxvars(iaux)%Strip
  enddo  

end subroutine TOilImsAuxVarArray1Strip

! ************************************************************************** !

subroutine TOilImsAuxVarArray2Strip(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes 
  ! using class(*) (unlimited polymorphic)
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/20/15
  ! 

  use AuxVars_TOilIms_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !class(auxvar_toil_ims_type) :: auxvars(0:,:)
  !cannot use type(...) with pointer attribute.
  !because the compiler does not allow to specify lower 0-bound in auxvar
  !type(auxvar_toil_ims_type), pointer :: auxvars(:,:)
  !class(auxvar_toil_ims_type) :: auxvars(0:,:)
  type(auxvar_toil_ims_type) :: auxvars(0:,:)

  PetscInt :: iaux, idof

  do iaux = 1, size(auxvars,2)
    do idof = 1, size(auxvars,1)
      !print *, "i am before TOilImsAuxVarStrip"
      !call TOilImsAuxVarStrip(auxvars(idof-1,iaux))
      call auxvars(idof-1,iaux)%Strip()
    enddo
  enddo  

end subroutine TOilImsAuxVarArray2Strip

! ************************************************************************** !

end module PM_TOilIms_Aux_module
