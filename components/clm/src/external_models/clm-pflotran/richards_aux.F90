module Richards_Aux_module

#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: richards_itol_scaled_res = 1.d-5
  PetscReal, public :: richards_itol_rel_update = UNINITIALIZED_DOUBLE

  type, public :: richards_auxvar_type
  
    PetscReal :: pc
#ifdef USE_ANISOTROPIC_MOBILITY
    PetscReal :: kvr_x
    PetscReal :: kvr_y
    PetscReal :: kvr_z
    PetscReal :: dkvr_x_dp
    PetscReal :: dkvr_y_dp
    PetscReal :: dkvr_z_dp
#else
    PetscReal :: kvr
    PetscReal :: dkvr_dp
#endif
    PetscReal :: dsat_dp
    PetscReal :: dden_dp
#ifdef CLM_PFLOTRAN
    PetscReal :: bc_alpha  ! Brooks Corey parameters: alpha
    PetscReal :: bc_lambda ! Brooks Corey parameters: lambda
    PetscReal :: bc_sr1    ! Brooks Corey parameters: sr(1)
#endif

    ! OLD-VAR-NAMES            = NEW-VAR
    ! ------------------------------------------------
    ! P_min                    = vars_for_sflow(1)
    ! P_max                    = vars_for_sflow(2)
    ! coeff_for_cubic_approx   = vars_for_sflow(3:6)
    ! range_for_linear_approx  = vars_for_sflow(7:10)
    ! bcflux_default_scheme    = vars_for_sflow(11)
    PetscReal, pointer :: vars_for_sflow(:)

  end type richards_auxvar_type
  
  type, public :: richards_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscBool :: auxvars_up_to_date
    PetscBool :: auxvars_cell_pressures_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(richards_auxvar_type), pointer :: auxvars(:)
    type(richards_auxvar_type), pointer :: auxvars_bc(:)
    type(richards_auxvar_type), pointer :: auxvars_ss(:)
#ifdef BUFFER_MATRIX
    type(matrix_buffer_type), pointer :: matrix_buffer
#endif
  end type richards_type

  public :: RichardsAuxCreate, RichardsAuxDestroy, &
            RichardsAuxVarCompute, RichardsAuxVarInit, &
            RichardsAuxVarCopy

contains

! ************************************************************************** !

function RichardsAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(richards_type), pointer :: RichardsAuxCreate
  
  type(richards_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%auxvars_cell_pressures_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
#ifdef BUFFER_MATRIX
  nullify(aux%matrix_buffer)
#endif

  RichardsAuxCreate => aux
  
end function RichardsAuxCreate

! ************************************************************************** !

subroutine RichardsAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%pc = 0.d0

#ifdef USE_ANISOTROPIC_MOBILITY
  auxvar%kvr_x = 0.d0
  auxvar%kvr_y = 0.d0
  auxvar%kvr_z = 0.d0
  auxvar%dkvr_x_dp = 0.d0
  auxvar%dkvr_y_dp = 0.d0
  auxvar%dkvr_z_dp = 0.d0
#else
  auxvar%kvr = 0.d0
  auxvar%dkvr_dp = 0.d0
#endif

  auxvar%dsat_dp = 0.d0
  auxvar%dden_dp = 0.d0

#ifdef CLM_PFLOTRAN
  auxvar%bc_alpha  = UNINITIALIZED_DOUBLE
  auxvar%bc_lambda = UNINITIALIZED_DOUBLE
  auxvar%bc_sr1    = UNINITIALIZED_DOUBLE
#endif

  if (option%surf_flow_on) then
    allocate(auxvar%vars_for_sflow(11))
    auxvar%vars_for_sflow(:) = 0.d0
  else
    nullify(auxvar%vars_for_sflow)
  endif
  

end subroutine RichardsAuxVarInit

! ************************************************************************** !

subroutine RichardsAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pc = auxvar%pc

#ifdef USE_ANISOTROPIC_MOBILITY
  auxvar2%kvr_x = auxvar%kvr_x 
  auxvar2%kvr_y = auxvar%kvr_y 
  auxvar2%kvr_z = auxvar%kvr_z 
  auxvar2%dkvr_x_dp = auxvar%dkvr_x_dp 
  auxvar2%dkvr_y_dp = auxvar%dkvr_y_dp 
  auxvar2%dkvr_z_dp = auxvar%dkvr_z_dp 
#else
  auxvar2%kvr = auxvar%kvr
  auxvar2%dkvr_dp = auxvar%dkvr_dp
#endif

  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dden_dp = auxvar%dden_dp
 
#ifdef CLM_PFLOTRAN
  auxvar2%bc_alpha  = auxvar%bc_alpha
  auxvar2%bc_lambda = auxvar%bc_lambda
  auxvar2%bc_sr1    = auxvar%bc_sr1
#endif

  if (option%surf_flow_on) &
    auxvar2%vars_for_sflow(:) = auxvar%vars_for_sflow(:)

end subroutine RichardsAuxVarCopy

! ************************************************************************** !

subroutine RichardsAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                 characteristic_curves,option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(richards_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  PetscInt :: i
  PetscBool :: saturated
  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: pert, pw_pert, dw_kg_pert
  PetscReal :: fs, ani_A, ani_B, ani_C, ani_n, ani_coef
  PetscReal :: dkr_sat
  PetscReal :: aux(1)
  PetscReal, parameter :: tol = 1.d-3
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscReal :: pcmax
  
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0
  
#ifdef USE_ANISOTROPIC_MOBILITY
  auxvar%kvr_x = 0.d0
  auxvar%kvr_y = 0.d0
  auxvar%kvr_z = 0.d0
#else
  auxvar%kvr = 0.d0
#endif

  kr = 0.d0
 
  global_auxvar%pres = x(1)
!if coupled with CLM, CLM will update temperature via the interface
#ifndef CLM_PFLOTRAN

  global_auxvar%temp = option%reference_temperature

#else
  ! Check if the capillary pressure is less than prescribed 'pcmax' (positive)
  pcmax = -characteristic_curves%saturation_function%pcmax + 1.d0
  if (global_auxvar%pres(1) - option%reference_pressure < pcmax) then
    global_auxvar%pres(1) = pcmax + option%reference_pressure
  endif

#endif
 
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)
  auxvar%pc = max(0.d0, auxvar%pc)   ! always non-negative (0 = saturated)
  
!***************  Liquid phase properties **************************
  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (auxvar%pc > 0.d0) then

#ifdef CLM_PFLOTRAN
  if(auxvar%bc_alpha /= UNINITIALIZED_DOUBLE .and. &
     auxvar%bc_lambda /= UNINITIALIZED_DOUBLE .and. &
     auxvar%bc_sr1 /= UNINITIALIZED_DOUBLE ) then

    ! fmy: the following only needs calling ONCE, but not yet figured out how
    ! because CLM's every single CELL has ONE set of SF/RPF parameters
    select type(sf => characteristic_curves%saturation_function)
      !class is(sat_func_VG_type)
        ! not-yet
      class is(sat_func_BC_type)
        sf%alpha  = auxvar%bc_alpha
        sf%lambda = auxvar%bc_lambda
        sf%Sr  = auxvar%bc_sr1
        ! needs to re-calculate some extra variables for 'saturation_function', if changed above
        if(associated(sf%sat_poly) .or. associated(sf%pres_poly)) then
          error_string = 'passing CLM characterisitc-curves parameters: sat_function'
          call sf%SetupPolynomials(option,error_string)
        endif

      class default
        option%io_buffer = 'Currently ONLY support Brooks_COREY saturation function type' // &
          ' when coupled with CLM.'
        call printErrMsg(option)
    end select

    select type(rpf => characteristic_curves%liq_rel_perm_function)
      !class is(rpf_Mualem_VG_liq_type)
      class is(rpf_Burdine_BC_liq_type)
        rpf%lambda = auxvar%bc_lambda
        rpf%Sr  = auxvar%bc_sr1

        ! Burdine_BC_liq RPF has no spline-smoothing (@ May-05-2016)
        if(associated(rpf%poly)) then
          error_string = 'passing CLM characterisitc-curves parameters: rpf_function'
          call rpf%SetupPolynomials(option,error_string)
        endif

      class default
        option%io_buffer = 'Currently ONLY support Brooks_COREY-Burdine liq. permissivity function type' // &
          ' when coupled with CLM.'
        call printErrMsg(option)
    end select

  endif
#endif

    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(auxvar%pc, &
                                          global_auxvar%sat(1), &
                                          ds_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (ds_dp < 1.d-40 .and. auxvar%pc<1.0d1) then
      saturated = PETSC_TRUE
    else
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(global_auxvar%sat(1), &
                                            kr,dkr_sat,option)
      dkr_dp = ds_dp * dkr_sat
    endif
  else
    saturated = PETSC_TRUE
  endif  
  
  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when 
  ! capillary pressure is very small
  if (saturated) then
    auxvar%pc = 0.d0
    global_auxvar%sat = 1.d0  
    kr = 1.d0    
    pw = global_auxvar%pres(1)
  endif

  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol, &
                         dw_dp,dw_dt,ierr)
    ! may need to compute dpsat_dt to pass to VISW
    call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,ierr)
    !geh: 0.d0 passed in for derivative of pressure w/respect to temp
    call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,0.d0, &
                           visl,dvis_dt,dvis_dp,ierr) 
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(global_auxvar%temp,pw,aux, &
                            dw_kg,dw_mol,dw_dp,dw_dt,ierr)
    call EOSWaterViscosityExt(global_auxvar%temp,pw,sat_pressure,0.d0,aux, &
                              visl,dvis_dt,dvis_dp,ierr) 
  endif
  if (.not.saturated) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dp = dw_dp
  auxvar%kvr = kr/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  
end subroutine RichardsAuxVarCompute

! ************************************************************************** !

subroutine AuxVarDestroy(auxvar)
  ! 
  ! Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(richards_auxvar_type) :: auxvar
  
end subroutine AuxVarDestroy

! ************************************************************************** !

subroutine RichardsAuxDestroy(aux)
  ! 
  ! Deallocates a richards auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(richards_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call AuxVarDestroy(aux%auxvars_bc(iaux))
    enddo  
    deallocate(aux%auxvars_bc)
  endif
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call AuxVarDestroy(aux%auxvars_ss(iaux))
    enddo  
    deallocate(aux%auxvars_ss)
  endif
  nullify(aux%auxvars_ss)
  
  call DeallocateArray(aux%zero_rows_local)
  call DeallocateArray(aux%zero_rows_local_ghosted)

#ifdef BUFFER_MATRIX
  if (associated(aux%matrix_buffer)) then
    call MatrixBufferDestroy(aux%matrix_buffer)
  endif
  nullify(aux%matrix_buffer)
#endif
  
  deallocate(aux)
  nullify(aux)
    
end subroutine RichardsAuxDestroy

end module Richards_Aux_module
