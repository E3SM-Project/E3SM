module TOilIms_module
! Brief description for the module
! Pimary variables for ToilIms: oil_pressure, oil_saturation, temperature

  !use TOilIms_Aux_module
  use PM_TOilIms_Aux_module
  use AuxVars_TOilIms_module 

  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petsclog.h"

#define TOIL_CONVECTION
#define TOIL_CONDUCTION

! Cutoff parameters - no public
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 1.d-24

  !pointing to default function
  procedure(TOilImsFluxDummy), pointer :: TOilImsFlux => null()

  abstract interface
    subroutine TOilImsFluxDummy(toil_auxvar_up,global_auxvar_up, &
                     material_auxvar_up,sir_up, thermal_conductivity_up, &
                     toil_auxvar_dn,global_auxvar_dn,material_auxvar_dn, &
                     sir_dn, thermal_conductivity_dn, area, dist, parameter, &
                     option,v_darcy,Res)
      
      use AuxVars_TOilIms_module
      use PM_TOilIms_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_class
      use Connection_module
 
      implicit none
  
      class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: sir_up(:), sir_dn(:)
      PetscReal :: v_darcy(option%nphase)
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(toil_ims_parameter_type) :: parameter
      PetscReal :: thermal_conductivity_dn(2)
      PetscReal :: thermal_conductivity_up(2)
      PetscReal :: Res(option%nflowdof)
                                     
    end subroutine TOilImsFluxDummy

  end interface

  public :: TOilImsSetup, &
            TOilImsFluxDipcSetup, &
            TOilImsDefaultSetup, & 
            TOilImsUpdateAuxVars, &
            TOilImsInitializeTimestep, &
            TOilImsComputeMassBalance, &
            TOilImsResidual, &
            ToilImsJacobian, &
            TOilImsUpdateSolution, &
            TOilImsTimeCut, &
            TOilImsMapBCAuxVarsToGlobal, &
            !TOilImsCheckUpdatePre, &
            !TOilImsCheckUpdatePost, &
            TOilImsDestroy

contains

! ************************************************************************** !

subroutine TOilImsSetup(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/20/15
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  !use Fluid_module
  use Material_Aux_class
  use Output_Aux_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter
  type(output_variable_list_type), pointer :: list

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: num_bc_connection, num_ss_connection
  PetscInt :: i, idof, count
  PetscBool :: error_found
  PetscInt :: flag(10)

  class(material_auxvar_type), pointer :: material_auxvars(:)
  !type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%TOil_ims => TOilImsAuxCreate(option)

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE

  if (minval(material_parameter%soil_residual_saturation(:,:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil residual saturation.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_heat_capacity(:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil heat capacity.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_thermal_conductivity(:,:)) < 0.d0) then
    option%io_buffer = 'Non-initialized soil thermal conductivity.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  
  material_auxvars => patch%aux%Material%auxvars
  flag = 0

  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'Non-initialized cell volume.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'Non-initialized porosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'Non-initialized tortuosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%soil_particle_density < 0.d0 .and. &
        flag(4) == 0) then
      flag(4) = 1
      option%io_buffer = 'Non-initialized soil particle density.'
      call printMsg(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'Non-initialized permeability.'
      call printMsg(option)
      flag(5) = 1
    endif
  enddo

  if (error_found .or. maxval(flag) > 0) then
    option%io_buffer = 'Material property errors found in TOilImsSetup.'
    call printErrMsg(option)
  endif

  num_bc_connection = &
               CouplerGetNumConnectionsInList(patch%boundary_condition_list)

  num_ss_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)

  call patch%aux%TOil_ims%Init(grid,num_bc_connection,num_ss_connection,option)

  list => realization%output_option%output_snap_variable_list
  call TOilImsSetPlotVariables(list)
  list => realization%output_option%output_obs_variable_list
  call TOilImsSetPlotVariables(list)
 
  ! covergence creteria to be chosen (can use TOUGH or general type) 
  !if (general_tough2_conv_criteria .and. &
  !    Initialized(option%flow%inf_scaled_res_tol)) then
  !  ! override what was set in OPTION block of GENERAL process model
  !  general_tough2_itol_scaled_res_e1 = option%flow%inf_scaled_res_tol
  !endif

end subroutine TOilImsSetup

! ************************************************************************** !
subroutine TOilImsDefaultSetup()
  ! 
  ! Setup Dipc Flux computation
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/10/16

  implicit none
  
  TOilImsFlux => TOilImsFluxPFL

end subroutine TOilImsDefaultSetup

! ************************************************************************** !
subroutine TOilImsFluxDipcSetup()
  ! 
  ! Setup Dipc Flux computation
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/10/16

  implicit none
  
  TOilImsFlux => TOilImsFluxDipc

end subroutine TOilImsFluxDipcSetup

! ************************************************************************** !

subroutine TOilImsInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/20/15
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call TOilImsUpdateFixedAccum(realization)
  

end subroutine TOilImsInitializeTimestep

! ************************************************************************** !

subroutine TOilImsSetPlotVariables(list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/20/15
  ! 
  
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(output_variable_list_type), pointer :: list
  
  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_type), pointer :: output_variable

  if (associated(list%first)) then
    return
  endif
  
  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               TEMPERATURE)

  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)

  name = 'Oil Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               OIL_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)
  
  name = 'Oil Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               OIL_SATURATION)
  
  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)
  
  name = 'Oil Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               OIL_DENSITY)
  
  name = 'Liquid Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)
  
  name = 'Oil Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               OIL_ENERGY)
  
 !name = 'Thermodynamic State'
 ! units = ''
 ! output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE,units,STATE)
 ! output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
 ! output_variable%iformat = 1 ! integer
 ! call OutputVariableAddToList(list,output_variable)   
  
end subroutine TOilImsSetPlotVariables

! ************************************************************************** !

subroutine TOilImsTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15
  ! 
  use Realization_Subsurface_class
  !use Option_module
  !use Field_module
  !use Patch_module
  !use Discretization_module
  !use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  ! if anything else to do when specific to cutting time step - 
  ! for TOIL IMS - add here

  call TOilImsInitializeTimestep(realization)

end subroutine TOilImsTimeCut

! ************************************************************************** !


! ************************************************************************** !

subroutine TOilImsUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with the TOilIms problem
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/21/15 - 28/05/2016
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_class
  !use EOS_Water_module 
  !use Saturation_Function_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  !commented for new auxvars data structure
  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:), toil_auxvars_bc(:)   
  !class(pm_toil_ims_aux_type), pointer :: toil_aux
 
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)  

  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end
  !PetscInt :: iphasebc, iphase
  PetscInt :: offset
  PetscInt :: istate

  !PetscReal :: gas_pressure, capillary_pressure, liquid_saturation
  !PetscReal :: saturation_pressure, temperature

  PetscInt :: real_index, variable
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)

  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  !idea to generalise pm_xxx_aux using three classes: 
  !  pm_base_aux,pm_subsurface_flow_aux,pm_xxx_aux
  ! pm_subsurface_flow_aux extension of pm_base_aux
  ! type, public, extends(pm_base_aux) :: pm_subsurface_flow_aux
  !   class(auxvar_base_type) :: auxvars_base  !pointer as in well
  !   ! class(auxvar_flow_energy_type) :: auxvars_flow_energy !pointer as in wells
  !   ....................
  ! end 
  ! 
  ! (from TOilImsUpdateAuxVars calls )
  ! PMSubSurfFlowUpdateAuxVars(grid,field,global_auxvars,characteristic_curves_array,..) 
  ! PMSubSurfFlowUpdateAuxVars as member function of pm_subsurface_flow_aux
  !
  ! ...............................
  ! do ghosted_id = 1, grid%ngmax
  !   pm_subsurface_flow_aux%auxvars_flow_energy(ZERO_INTEGER,ghosted_id)%compute(
  !                       grid,field,global_auxvars,characteristic_curves_array
  !                        option)
  ! end do

  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    ! TOIL_IMS_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = TOIL_IMS_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    call TOilImsAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id), & 
                       !toil_auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%characteristic_curves_array( &
                         patch%sat_func_id(ghosted_id))%ptr, &
                       natural_id, &
                       option)
   ! if TOilImsAuxVarCompute becomes a member function of auxvar_toil_ims
   ! call patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id)%compute( &
   !                   xx_loc_p(ghosted_start:ghosted_end), &
   !                   global_auxvars(ghosted_id), &
   !                   material_auxvars(ghosted_id), &
   !                   patch%characteristic_curves_array( &
   !                   patch%sat_func_id(ghosted_id))%ptr, &
   !                   natural_id, &
   !                   option)
   ! we will be forced to have a Compute member common to all modes
   ! where even in auxvar_base we must have the same list of dummy arguments
   ! in this case: xx_loc_p, global_auxvar, material_auxvar, characteristic_curves
   ! natural_id, option
   ! basically forced to reuse as much code as possible
   ! an example is auxvars()%Init already implemented
   
  enddo
  
  ! compute auxiliary variables for boundary cells
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      !negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id) 
      offset = (ghosted_id-1)*option%nflowdof
      if (patch%imat(ghosted_id) <= 0) cycle

      xxbc(:) = xx_loc_p(offset+1:offset+option%nflowdof)
      !istate = boundary_condition%flow_aux_int_var(GENERAL_STATE_INDEX,iconn)

      ! we do this for all BCs; Neumann bcs will be set later
      do idof = 1, option%nflowdof
        real_index = boundary_condition% &
                       flow_aux_mapping(toil_ims_dof_to_primary_vars(idof))
        if (real_index > 0) then
          xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
        else
          option%io_buffer = 'Error setting up boundary condition' // &
                             ' in TOilImsUpdateAuxVars'
          call printErrMsg(option)
        endif
      enddo
        
      ! state not required 
      !toil_auxvars_bc(sum_connection)%istate = istate
      ! TOIL_IMS_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = TOIL_IMS_UPDATE_FOR_BOUNDARY
      !call TOilImsAuxVarCompute(xxbc,toil_auxvars_bc(sum_connection), &
      call TOilImsAuxVarCompute(xxbc, &
                                patch%aux%TOil_ims%auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  !patch%aux%TOil_ims%auxvars_up_to_date = PETSC_TRUE
  patch%aux%TOil_ims%auxvars_up_to_date = PETSC_TRUE 

end subroutine TOilImsUpdateAuxVars

! ************************************************************************** !

subroutine TOilImsUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)

  type(global_auxvar_type), pointer :: global_auxvars(:)

  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  !new auxvar data strucutre
  !toil_auxvars => patch%aux%TOil_ims%auxvars

  global_auxvars => patch%aux%Global%auxvars

  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  !Tough2 conv. creteria: initialize accumulation term for every iteration
  if (toil_ims_tough2_conv_criteria) then
    call VecGetArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! TOIL_IMS_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = TOIL_IMS_UPDATE_FOR_FIXED_ACCUM ! not currently used
    call TOilImsAuxVarCompute(xx_p(local_start:local_end), &
                              !toil_auxvars(ZERO_INTEGER,ghosted_id), &
                              patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              ghosted_id, &
                              option)
    !call TOilImsAccumulation(toil_auxvars(ZERO_INTEGER,ghosted_id), &
    call TOilImsAccumulation( &
                          patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id), &
                          material_auxvars(ghosted_id), &
                          material_parameter%soil_heat_capacity(imat), &
                          option,accum_p(local_start:local_end) )
  enddo
  
  !Tough2 conv. creteria: initialize accumulation term for every iteration
  if (toil_ims_tough2_conv_criteria) then
    accum_p2 = accum_p
  endif
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  !Tough2 conv. creteria: initialize accumulation term for every iteration
  if (toil_ims_tough2_conv_criteria) then
    call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
end subroutine TOilImsUpdateFixedAccum

! ************************************************************************** !

subroutine TOilImsUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step: currently it updates only mass balance
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call TOilImsUpdateMassBalance(realization)
  endif
  
end subroutine TOilImsUpdateSolution

! ************************************************************************** !

subroutine TOilImsComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/12/15
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase !, icomp
  PetscReal :: vol_phase

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  !toil_auxvars => patch%aux%TOil_ims%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0
  ! note ::  only first column of mass_balance(1:2,1) is used

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = & 
        patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density

        mass_balance(iphase,1) = mass_balance(iphase,1) + &
          patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)* &
          toil_ims_fmw_comp(iphase)*vol_phase

    enddo
  enddo

end subroutine TOilImsComputeMassBalance

! ************************************************************************** !

subroutine TOilImsUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! Using existing data structure for two phase compositional
  ! For memory efficiency define new data structure 
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
  use EOS_Oil_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  
  PetscInt :: iconn
  PetscInt :: icomp

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  ! option%nflowspec = 2,
  ! two species (H2O,OIL): each present only in its own rich phase

  ! updating with mass balance, assuming molar quantity loaded in:
  ! mass_balance and mass_balance_delta

  do iconn = 1, patch%aux%TOil_ims%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        toil_ims_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%TOil_ims%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        toil_ims_fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine TOilImsUpdateMassBalance

! ************************************************************************** !

subroutine TOilImsZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%TOil_Ims%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%TOil_Ims%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine TOilImsZeroMassBalanceDelta

! ************************************************************************** !

subroutine TOilImsMapBCAuxVarsToGlobal(realization)
  ! 
  ! Maps toil_ims BC AuxVars to global auxvars
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  !type(toil_ims_auxvar_type), pointer :: toil_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  !NOTE: this is called only if cpoupling to RT
  if (option%ntrandof == 0) return ! no need to update

  !new auxvar data structure  
  !toil_auxvars_bc => patch%aux%TOil_ims%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        !toil_auxvars_bc(sum_connection)%sat
        patch%aux%TOil_ims%auxvars_bc(sum_connection)%sat 
      global_auxvars_bc(sum_connection)%den_kg = &
        !toil_auxvars_bc(sum_connection)%den_kg
        patch%aux%TOil_ims%auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        !toil_auxvars_bc(sum_connection)%temp
        patch%aux%TOil_ims%auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine TOilImsMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine TOilImsAccumulation(toil_auxvar,material_auxvar, &
                               soil_heat_capacity,option,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/23/15
  ! 

  use Option_module
  use Material_module
  use Material_Aux_class
  
  implicit none

  !type(toil_ims_auxvar_type) :: toil_auxvar
  class(auxvar_toil_ims_type) :: toil_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: soil_heat_capacity
  type(option_type) :: option
  PetscReal :: Res(option%nflowdof) 
  !PetscBool :: debug_cell
  
  PetscInt :: iphase, energy_id
  
  PetscReal :: porosity
  PetscReal :: v_over_t
  
  energy_id = option%energy_id
  
  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  v_over_t = material_auxvar%volume / option%flow_dt
  ! must use toil_auxvar%effective porosity here as it enables numerical 
  ! derivatives to be employed 
  porosity = toil_auxvar%effective_porosity
  
  ! flow accumulation term units = kmol/s
  Res = 0.d0
  do iphase = 1, option%nphase
    ! Res[kmol comp/m^3 void] = sat[m^3 phase/m^3 void] * 
    !                           den[kmol phase/m^3 phase] 
    ! molar balance formulation (kmol)
    Res(iphase) = Res(iphase) + toil_auxvar%sat(iphase) * &
                                toil_auxvar%den(iphase)  
  enddo

  ! scale by porosity * volume / dt
  ! Res[kmol/sec] = Res[kmol/m^3 void] * por[m^3 void/m^3 bulk] * 
  !                 vol[m^3 bulk] / dt[sec]
  Res(1:option%nphase) = Res(1:option%nphase) * &
                            porosity * v_over_t

  ! energy accumulation term units = MJ/s
  do iphase = 1, option%nphase
    ! Res[MJ/m^3 void] = sat[m^3 phase/m^3 void] *
    !                    den[kmol phase/m^3 phase] * U[MJ/kmol phase]
    Res(energy_id) = Res(energy_id) + toil_auxvar%sat(iphase) * &
                                      toil_auxvar%den(iphase) * &
                                      toil_auxvar%U(iphase)
  enddo
  ! Res[MJ/sec] = (Res[MJ/m^3 void] * por[m^3 void/m^3 bulk] + 
  !                (1-por)[m^3 rock/m^3 bulk] * 
  !                  dencpr[kg rock/m^3 rock * MJ/kg rock-K] * T[C]) &
  !               vol[m^3 bulk] / dt[sec]
  Res(energy_id) = (Res(energy_id) * porosity + &
                    (1.d0 - porosity) * &
                    material_auxvar%soil_particle_density * &
                    soil_heat_capacity * toil_auxvar%temp) * v_over_t
                    
end subroutine TOilImsAccumulation

! ************************************************************************** !

! ************************************************************************** !

subroutine TOilImsFluxPFL(toil_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       sir_up, &
                       thermal_conductivity_up, &
                       toil_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area, dist, parameter, &
                       option,v_darcy,Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/27/15
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
 
  ! no fractures considered for now
  ! use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  !type(toil_ims_auxvar_type) :: toil_auxvar_up, toil_auxvar_dn
  class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(toil_ims_parameter_type) :: parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  !PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscInt :: energy_id
  PetscInt :: iphase
 
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn           ! no mole fractions
  PetscReal :: delta_pressure, delta_temp !, delta_xmol,

  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux

  ! no diff fluxes - arrays used for debugging only
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)
  
  PetscReal :: dummy_perm_up, dummy_perm_dn

  energy_id = option%energy_id

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  
  ! Fracture permeability change only available for structured grid (Heeho)
  ! PO no fractures considered for now
  !if (associated(material_auxvar_up%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_up,perm_up,perm_up, &
  !                            dummy_perm_up,dist)
  !endif
  !if (associated(material_auxvar_dn%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_dn,perm_dn,perm_dn, &
  !                            dummy_perm_dn,dist)
  !endif
  
  !if (associated(klinkenberg)) then
  !  perm_ave_over_dist(1) = (perm_up * perm_dn) / &
  !                          (dist_up*perm_dn + dist_dn*perm_up)
  !  dummy_perm_up = klinkenberg%Evaluate(perm_up, &
  !                                       gen_auxvar_up%pres(option%gas_phase))
  !  dummy_perm_dn = klinkenberg%Evaluate(perm_dn, &
  !                                       gen_auxvar_dn%pres(option%gas_phase))
  !  perm_ave_over_dist(2) = (dummy_perm_up * dummy_perm_dn) / &
  !                          (dist_up*dummy_perm_dn + dist_dn*dummy_perm_up)
  !else
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
  !endif
      
  Res = 0.d0
  
  v_darcy = 0.d0

!#ifdef DEBUG_FLUXES  
!  adv_flux = 0.d0
!  diff_flux = 0.d0
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  debug_flux = 0.d0
!  debug_dphi = 0.d0
!#endif

#ifdef TOIL_CONVECTION
  do iphase = 1, option%nphase
 
    if (toil_auxvar_up%mobility(iphase) + &
        toil_auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif

    ! an alternative could be to avergae using oil_sat
    !density_kg_ave = 0.5d0* ( toil_auxvar_up%den_kg(iphase) + &
    !                          toil_auxvar_dn%den_kg(iphase) )
    density_kg_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                     toil_auxvar_dn%sat(iphase), &
                     toil_auxvar_up%den_kg(iphase), &
                     toil_auxvar_dn%den_kg(iphase))

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = toil_auxvar_up%pres(iphase) - &
                     toil_auxvar_dn%pres(iphase) + &
                     gravity_term

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_dphi(iphase) = delta_pressure
!#endif

    ! upwinding the mobilities and enthalpies
    if (delta_pressure >= 0.D0) then
      mobility = toil_auxvar_up%mobility(iphase)
      H_ave = toil_auxvar_up%H(iphase)
      uH = H_ave
      !density_ave = toil_auxvar_up%den(iphase)
    else
      mobility = toil_auxvar_dn%mobility(iphase)
      H_ave = toil_auxvar_dn%H(iphase)
      uH = H_ave
      !density_ave = toil_auxvar_dn%den(iphase)
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility * delta_pressure

      ! if comments below, use upwinding value
      !density_ave = 0.5d0*( toil_auxvar_up%den(iphase) + &
      !                      toil_auxvar_dn%den(iphase))

      density_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den(iphase), &
                           toil_auxvar_dn%den(iphase))       
 
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave
      ! Res[kmol total/sec]

      ! Res[kmol phase/sec] = mole_flux[kmol phase/sec]  
      Res(iphase) = Res(iphase) + mole_flux 

      !do icomp = 1, option%nflowspec
      !  ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
      !  !                      xmol[kmol comp/kmol phase]
      !  Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      !enddo

!#ifdef DEBUG_FLUXES  
!      do icomp = 1, option%nflowspec
!        adv_flux(icomp) = adv_flux(icomp) + mole_flux * xmol(icomp)
!      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      do icomp = 1, option%nflowspec
!        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
!#endif

      Res(energy_id) = Res(energy_id) + mole_flux * uH

!#ifdef DEBUG_FLUXES  
!      adv_flux(energy_id) = adv_flux(energy_id) + mole_flux * uH
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_dphi(iphase) = delta_pressure
!      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
!#endif

    endif  ! if mobility larger than given tolerance                 

  enddo
#endif 
! TOIL_CONVECTION

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then  
!    write(debug_unit,'(a,7es24.15)') 'delta pressure :', debug_dphi(:)
!    write(debug_unit,'(a,7es24.15)') 'adv flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'adv flux (gas):', debug_flux(:,2)
!  endif
!  debug_flux = 0.d0
!#endif                    

#ifdef TOIL_CONDUCTION
  ! model for liquid + gas
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
  !k_eff_up = thermal_conductivity_up(1) + &
  !           sqrt(gen_auxvar_up%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  !k_eff_dn = thermal_conductivity_dn(1) + &
  !           sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  !if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
  !  k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  !else
  !  k_eff_ave = 0.d0
  !endif
  ! considered the formation fully saturated in water for heat conduction 
  k_eff_up = thermal_conductivity_up(1)
  k_eff_dn = thermal_conductivity_dn(1)
  if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  else
    k_eff_ave = 0.d0
  endif

  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = toil_auxvar_up%temp - toil_auxvar_dn%temp
  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux
! CONDUCTION
#endif

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  debug_flux(energy_id,1) = debug_flux(energy_id,1) + heat_flux
!  if (debug_flag > 0) then  
!    write(debug_unit,'(a,7es24.15)') 'dif flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'dif flux (gas):', debug_flux(:,2)
!  endif
!#endif

end subroutine TOilImsFluxPFL

! ************************************************************************** !


! ************************************************************************** !
subroutine TOilImsFluxDipc(toil_auxvar_up,global_auxvar_up, &
                           material_auxvar_up, &
                           sir_up, &
                           thermal_conductivity_up, &
                           toil_auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           sir_dn, &
                           thermal_conductivity_dn, &
                           area, dist, parameter, &
                           option,v_darcy,Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! using a dip correction for distorted mesh
  ! The function assumes that z is the vertical direction  
  ! 
  ! Author: Paolo Orsini
  ! Date: 9/13/16
  ! 
  use Option_module
  use Material_Aux_class
  use Connection_module
  
  implicit none
  
  !type(toil_ims_auxvar_type) :: toil_auxvar_up, toil_auxvar_dn
  class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: v_darcy(option%nphase)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(toil_ims_parameter_type) :: parameter
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: Res(option%nflowdof)
  !PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscInt :: energy_id
  PetscInt :: iphase
 
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: uH
  PetscReal :: H_ave
  PetscReal :: perm_ave_over_dist(option%nphase)
  PetscReal :: perm_up, perm_dn           ! no mole fractions
  PetscReal :: delta_pressure, delta_temp !, delta_xmol,

  PetscReal :: pressure_ave
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: stpd_up, stpd_dn
  PetscReal :: sat_up, sat_dn, den_up, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, tempreal
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux

  ! no diff fluxes - arrays used for debugging only
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)

  PetscReal :: dummy_perm_up, dummy_perm_dn

  PetscBool :: horizontal_conn !flag if the connection is horizontal or not
  PetscReal :: dist_hrz_vec(2)
  PetscReal :: dist_hrz_projection_sq, dist_hrz_projection 
  PetscReal :: dist_hrz_up, dist_hrz_dn
  PetscReal :: dist_vrt_projection_sq !vertical projection: cell centres dh 
  PetscReal :: dip_corr

  energy_id = option%energy_id  

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  !determine if the connection is horixontal or not - 
  !TODO - PO: should this method work move this operation to pre-processing
  !           to be done only once  

  horizontal_conn = PETSC_FALSE
  if ( dabs(dist(3)) > dabs(dist(1)) .and. dabs(dist(3)) > dabs(dist(2)) ) then
    horizontal_conn = PETSC_TRUE
  end if 

  dip_corr = 1.0d0

  !if vertical connection, compute dip correction term and horizontal dists
  if (.not.horizontal_conn) then
    dist_hrz_vec(1:2) = dist(1:2) * dist(0) 
    dist_hrz_projection_sq = dot_product(dist_hrz_vec,dist_hrz_vec) 
    dist_vrt_projection_sq = ( dist(3) * dist(0) ) * (dist(3) * dist(0)) 

    dip_corr =  dist_hrz_projection_sq /  &
                (dist_hrz_projection_sq + dist_vrt_projection_sq)

    dist_hrz_projection = dsqrt(dist_hrz_projection_sq)
    dist_hrz_up = dist_hrz_projection * dist(-1)
    dist_hrz_dn = dist_hrz_projection - dist_hrz_up
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_hrz_up*perm_dn + dist_hrz_dn*perm_up)
    !perm_ave_over_dist(:) = (perm_up * perm_dn) / &
    !                        (dist_up*perm_dn + dist_dn*perm_up)
  else 
    perm_ave_over_dist(:) = (perm_up * perm_dn) / &
                            (dist_up*perm_dn + dist_dn*perm_up)
  end if

  !write(*,*) "dip_corr = ", dip_corr
      
  Res = 0.d0
  
  v_darcy = 0.d0

#ifdef TOIL_CONVECTION
  do iphase = 1, option%nphase
 
    if (toil_auxvar_up%mobility(iphase) + &
        toil_auxvar_dn%mobility(iphase) < eps) then
      cycle
    endif

    ! an alternative could be to avergae using oil_sat
    !density_kg_ave = 0.5d0* ( toil_auxvar_up%den_kg(iphase) + &
    !                          toil_auxvar_dn%den_kg(iphase) )
    density_kg_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                     toil_auxvar_dn%sat(iphase), &
                     toil_auxvar_up%den_kg(iphase), &
                     toil_auxvar_dn%den_kg(iphase))

    gravity_term = density_kg_ave * dist_gravity
    delta_pressure = toil_auxvar_up%pres(iphase) - &
                     toil_auxvar_dn%pres(iphase) + &
                     gravity_term

    ! upwinding the mobilities and enthalpies
    if (delta_pressure >= 0.D0) then
      mobility = toil_auxvar_up%mobility(iphase)
      H_ave = toil_auxvar_up%H(iphase)
      uH = H_ave
      !density_ave = toil_auxvar_up%den(iphase)
    else
      mobility = toil_auxvar_dn%mobility(iphase)
      H_ave = toil_auxvar_dn%H(iphase)
      uH = H_ave
      !density_ave = toil_auxvar_dn%den(iphase)
    endif      

    if (mobility > floweps) then
      ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
      !                    dP[Pa]]
      v_darcy(iphase) = perm_ave_over_dist(iphase) * mobility *  &
                        delta_pressure * dip_corr 

      ! if comments below, use upwinding value
      !density_ave = 0.5d0*( toil_auxvar_up%den(iphase) + &
      !                      toil_auxvar_dn%den(iphase))

      density_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den(iphase), &
                           toil_auxvar_dn%den(iphase))       

      !ovewrite area computed as for OLDTRAN 
      !0.5 below assumes uniform grid in x and y, i.e. 2*dist 
      !should compute the half volumes of entire cell hrz extensions DXi & DXj
      if (.not.horizontal_conn) then
        area = 0.5*(material_auxvar_up%volume + material_auxvar_dn%volume) / &
                dist_hrz_projection
      else 
        area = 0.5*(material_auxvar_up%volume + material_auxvar_dn%volume) / &
                (dist_up + dist_dn)  
      endif

      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area  
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                             density_ave[kmol phase/m^3 phase]        
      mole_flux = q*density_ave
      ! Res[kmol total/sec]

      ! Res[kmol phase/sec] = mole_flux[kmol phase/sec]  
      Res(iphase) = Res(iphase) + mole_flux 

      Res(energy_id) = Res(energy_id) + mole_flux * uH

    endif  ! if mobility larger than given tolerance                 

  enddo
#endif 
 TOIL_CONVECTION

#ifdef TOIL_CONDUCTION
  ! model for liquid + gas
  ! add heat conduction flux
  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
  !k_eff_up = thermal_conductivity_up(1) + &
  !           sqrt(gen_auxvar_up%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_up(2) - thermal_conductivity_up(1))
  !k_eff_dn = thermal_conductivity_dn(1) + &
  !           sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
  !           (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
  !if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
  !  k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  !else
  !  k_eff_ave = 0.d0
  !endif
  ! considered the formation fully saturated in water for heat conduction 
  k_eff_up = thermal_conductivity_up(1)
  k_eff_dn = thermal_conductivity_dn(1)
  if (k_eff_up > 0.d0 .or. k_eff_up > 0.d0) then
    k_eff_ave = (k_eff_up*k_eff_dn)/(k_eff_up*dist_dn+k_eff_dn*dist_up)
  else
    k_eff_ave = 0.d0
  endif

  ! units:
  ! k_eff = W/K-m = J/s/K-m
  ! delta_temp = K
  ! area = m^2
  ! heat_flux = k_eff * delta_temp * area = J/s
  delta_temp = toil_auxvar_up%temp - toil_auxvar_dn%temp
  heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! J/s -> MJ/s
  ! MJ/s
  Res(energy_id) = Res(energy_id) + heat_flux
! CONDUCTION
#endif

end subroutine TOilImsFluxDipc
! ************************************************************************** !

subroutine TOilImsBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         toil_auxvar_up,global_auxvar_up, &
                         toil_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         sir_dn, &
                         thermal_conductivity_dn, &
                         area,dist,toil_parameter, &
                         option,v_darcy,Res)
  ! 
  ! Computes the boundary flux terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/27/15
  ! 
  use Option_module                              
  use Material_Aux_class
  !use Fracture_module
  !use Klinkenberg_module
  
  implicit none
  
  !type(toil_ims_auxvar_type) :: toil_auxvar_up, toil_auxvar_dn
  class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn 
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: v_darcy(option%nphase), area
  type(toil_ims_parameter_type) :: toil_parameter
  PetscReal :: dist(-1:3)
  PetscReal :: Res(1:option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(TOIL_IMS_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)
  !PetscBool :: debug_connection
  
  PetscInt :: energy_id
  PetscInt :: iphase
  PetscInt :: bc_type
  PetscReal :: density_ave, density_kg_ave
  PetscReal :: H_ave, uH
  PetscReal :: perm_dn_adj(option%nphase)
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure, delta_temp
  PetscReal :: gravity_term
  PetscReal :: mobility, mole_flux, q
  PetscReal :: sat_dn, perm_dn, den_dn
  PetscReal :: temp_ave, stpd_ave_over_dist, pres_ave
  PetscReal :: k_eff_up, k_eff_dn, k_eff_ave, heat_flux
  ! for debugging only
  PetscReal :: adv_flux(3,2), diff_flux(2,2)
  PetscReal :: debug_flux(3,3), debug_dphi(2)

  PetscReal :: boundary_pressure
  PetscReal :: tempreal
  
  PetscInt :: idof
  PetscBool :: neumann_bc_present
  
  PetscReal :: dummy_perm_dn
  
  energy_id = option%energy_id

  Res = 0.d0
  v_darcy = 0.d0 
 
!#ifdef DEBUG_FLUXES    
!  adv_flux = 0.d0
!  diff_flux = 0.d0
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  debug_flux = 0.d0
!  debug_dphi = 0.d0
!#endif

  neumann_bc_present = PETSC_FALSE
  
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  ! currently no fractures considered 
  ! Fracture permeability change only available for structured grid (Heeho)
  !if (associated(material_auxvar_dn%fracture)) then
  !  call FracturePermEvaluate(material_auxvar_dn,perm_dn,perm_dn, &
  !                            dummy_perm_dn,dist)
  !endif  
  
  !if (associated(klinkenberg)) then
  !  perm_dn_adj(1) = perm_dn
  !                                        
  !  perm_dn_adj(2) = klinkenberg%Evaluate(perm_dn, &
  !                                        gen_auxvar_dn%pres(option%gas_phase))
  !else
    perm_dn_adj(:) = perm_dn
  !endif
  
#ifdef TOIL_CONVECTION  
  do iphase = 1, option%nphase
 
    bc_type = ibndtype(iphase) ! loop over equations 1.Liq and 2.Oil
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)

        ! dist(0) = scalar - magnitude of distance
        ! gravity = vector(3)
        ! dist(1:3) = vector(3) - unit vector
        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
      
        if (bc_type == CONDUCTANCE_BC) then !not implemented yet
          select case(iphase)
            case(LIQUID_PHASE)
              idof = auxvar_mapping(TOIL_IMS_LIQ_CONDUCTANCE_INDEX)
            case(TOIL_IMS_OIL_PHASE)
              idof = auxvar_mapping(TOIL_IMS_OIL_CONDUCTANCE_INDEX)
          end select        
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn_adj(iphase) / dist(0)
        endif
        
        ! PO need to check what values of saturations are assigned to the BC ghost cells  
        ! using residual saturation cannot be correct! - geh
        ! reusing sir_dn for bounary auxvar
!#define BAD_MOVE1 ! this works
!#ifndef BAD_MOVE1       
        if (toil_auxvar_up%sat(iphase) > sir_dn(iphase) .or. &
            toil_auxvar_dn%sat(iphase) > sir_dn(iphase)) then
!#endif
          boundary_pressure = toil_auxvar_up%pres(iphase)

          ! PO no free surfce boundaries considered  
          !if (iphase == LIQUID_PHASE .and. &
          !    global_auxvar_up%istate == GAS_STATE) then
          !  ! the idea here is to accommodate a free surface boundary
          !  ! face.  this will not work for an interior grid cell as
          !  ! there should be capillary pressure in force.
          !  boundary_pressure = gen_auxvar_up%pres(option%gas_phase)
          !endif

          !density_kg_ave = 0.5d0 * (toil_auxvar_up%den_kg(iphase) + &
          !                          toil_auxvar_dn%den_kg(iphase) )

          density_kg_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den_kg(iphase), &
                           toil_auxvar_dn%den_kg(iphase))

          gravity_term = density_kg_ave * dist_gravity
          delta_pressure = boundary_pressure - &
                           toil_auxvar_dn%pres(iphase) + &
                           gravity_term
          

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!          debug_dphi(iphase) = delta_pressure
!#endif
          ! PO CONDUCTANCE_BC and SEEPAGE_BC not implemented
          if (bc_type == SEEPAGE_BC .or. &
              bc_type == CONDUCTANCE_BC) then
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                toil_auxvar_up%pres(iphase) - &
                 option%reference_pressure < eps) then
              delta_pressure = 0.d0
            endif
          endif
          
          !upwinding mobilities and enthalpies   
          if (delta_pressure >= 0.D0) then
            mobility = toil_auxvar_up%mobility(iphase)
            uH = toil_auxvar_up%H(iphase)
            !density_ave = toil_auxvar_up%den(iphase)
          else
            mobility = toil_auxvar_dn%mobility(iphase)
            uH = toil_auxvar_dn%H(iphase)
            !density_ave = toil_auxvar_dn%den(iphase)
          endif      

          if (mobility > floweps) then
            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
            v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
            ! only need average density if velocity > 0.

            ! when this is commented - using upwinding value
            !density_ave = 0.5d0 * (toil_auxvar_up%den(iphase) + &
            !                       toil_auxvar_dn%den(iphase) )
            density_ave = TOilImsAverageDensity(toil_auxvar_up%sat(iphase), &
                           toil_auxvar_dn%sat(iphase), &
                           toil_auxvar_up%den(iphase), &
                           toil_auxvar_dn%den(iphase))
          endif
!#ifndef BAD_MOVE1        
        endif ! sat > eps
!#endif

      case(NEUMANN_BC)
        select case(iphase)
          case(LIQUID_PHASE)
            idof = auxvar_mapping(TOIL_IMS_LIQUID_FLUX_INDEX)
          case(TOIL_IMS_OIL_PHASE)
            idof = auxvar_mapping(TOIL_IMS_OIL_FLUX_INDEX)
        end select
        
        neumann_bc_present = PETSC_TRUE
        if (dabs(auxvars(idof)) > floweps) then
          v_darcy(iphase) = auxvars(idof)
          if (v_darcy(iphase) > 0.d0) then 
            density_ave = toil_auxvar_up%den(iphase)
            uH = toil_auxvar_up%H(iphase)
          else 
            density_ave = toil_auxvar_dn%den(iphase)
            uH = toil_auxvar_dn%H(iphase)
          endif 
        endif
      case default
        option%io_buffer = &
          'Boundary condition type not recognized in GeneralBCFlux phase loop.'
        call printErrMsg(option)
    end select

    if (dabs(v_darcy(iphase)) > 0.d0) then
      ! q[m^3 phase/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy(iphase) * area
      ! mole_flux[kmol phase/sec] = q[m^3 phase/sec] * 
      !                              density_ave[kmol phase/m^3 phase]
      mole_flux = q*density_ave       
  
      ! Res[kmol total/sec]
      Res(iphase) = Res(iphase) + mole_flux
 
      ! Res[kmol total/sec]
      !do icomp = 1, option%nflowspec
      !  ! Res[kmol comp/sec] = mole_flux[kmol phase/sec] * 
      !  !                      xmol[kmol comp/mol phase]
      !  Res(icomp) = Res(icomp) + mole_flux * xmol(icomp)
      !enddo
!#ifdef DEBUG_FLUXES  
!      do icomp = 1, option%nflowspec
!        adv_flux(icomp,iphase) = adv_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      do icomp = 1, option%nflowspec
!        debug_flux(icomp,iphase) = debug_flux(icomp,iphase) + mole_flux * xmol(icomp)
!      enddo
!#endif
      ! Res[MJ/sec] = mole_flux[kmol comp/sec] * H_ave[MJ/kmol comp]
      Res(energy_id) = Res(energy_id) + mole_flux * uH ! H_ave
!#ifdef DEBUG_FLUXES  
!      adv_flux(energy_id,iphase) = adv_flux(energy_id,iphase) + mole_flux * uH
!#endif
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!      debug_flux(energy_id,iphase) = debug_flux(energy_id,iphase) + mole_flux * uH
!#endif
    endif
  enddo
#endif 
! end of TOIL_CONVECTION
  
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then 
!    write(debug_unit,'(a,7es24.15)') 'bc delta pressure :', debug_dphi(:)  
!    write(debug_unit,'(a,7es24.15)') 'bc adv flux (liquid):', debug_flux(:,1)
!    write(debug_unit,'(a,7es24.15)') 'bc adv flux (gas):', debug_flux(:,2)
!  endif
!  debug_flux = 0.d0
!#endif  


#ifdef TOIL_CONDUCTION
  ! add heat conduction flux
  heat_flux = 0.d0
  select case (ibndtype(TOIL_IMS_ENERGY_EQUATION_INDEX))
    case (DIRICHLET_BC)
      ! based on Somerton et al., 1974:
      ! k_eff = k_dry + sqrt(s_l)*(k_sat-k_dry)
      !k_eff_dn = thermal_conductivity_dn(1) + &
      !           sqrt(gen_auxvar_dn%sat(option%liquid_phase)) * &
      !           (thermal_conductivity_dn(2) - thermal_conductivity_dn(1))
      ! considered the formation fully saturated in water for heat conduction
      k_eff_dn = thermal_conductivity_dn(1)
      ! units:
      ! k_eff = W/K/m/m = J/s/K/m/m
      ! delta_temp = K
      ! area = m^2
      ! heat_flux = J/s
      k_eff_ave = k_eff_dn / dist(0)
      delta_temp = toil_auxvar_up%temp - toil_auxvar_dn%temp
      heat_flux = k_eff_ave * delta_temp * area * 1.d-6 ! convert W -> MW
    case(NEUMANN_BC)
                  ! flux prescribed as MW/m^2
      heat_flux = auxvars(auxvar_mapping(TOIL_IMS_ENERGY_FLUX_INDEX)) * area
    case(ZERO_GRADIENT_BC)
      ! No contribution to heat_flux
    case default
      option%io_buffer = 'Boundary condition type not recognized in ' // &
        'TOilImsBCFlux heat conduction loop.'
      call printErrMsg(option)
  end select
  Res(energy_id) = Res(energy_id) + heat_flux ! MW
#endif 
! end of TOIL_CONDUCTION

end subroutine TOilImsBCFlux

! ************************************************************************** !

subroutine TOilImsSrcSink(option,src_sink_condition, toil_auxvar, &
                          global_auxvar,ss_flow_vol_flux,scale,Res)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/04/15
  ! 

  use Option_module
  use Condition_module  

  use EOS_Water_module
  use EOS_Oil_module

  implicit none

  type(option_type) :: option
  type(flow_toil_ims_condition_type), pointer :: src_sink_condition
  !type(toil_ims_auxvar_type) :: toil_auxvar
  class(auxvar_toil_ims_type) :: toil_auxvar
  type(global_auxvar_type) :: global_auxvar !keep global_auxvar for salinity
  PetscReal :: ss_flow_vol_flux(option%nphase)
  PetscReal :: scale  
  PetscReal :: Res(option%nflowdof)

  ! local parameter
  PetscInt, parameter :: SRC_TEMPERATURE = 1
  PetscInt, parameter :: SRC_ENTHALPY = 2 
  ! local variables
  PetscReal, pointer :: qsrc(:)
  PetscInt :: flow_src_sink_type    
  PetscReal :: qsrc_mol
  PetscReal :: den, den_kg, enthalpy, internal_energy, temperature
  PetscReal :: cell_pressure, dummy_pressure
  PetscInt :: iphase
  PetscInt :: energy_var
  PetscErrorCode :: ierr

  ! this can be removed when etxending to pressure condition
  if (.not.associated(src_sink_condition%rate) ) then
    option%io_buffer = 'TOilImsSrcSink fow condition rate not defined ' // &
    'rate is needed for a valid src/sink term'
    call printErrMsg(option)  
  end if

  !qsrc => src_sink_condition%rate%dataset%rarray(:)
  qsrc => src_sink_condition%rate%dataset%rarray

  energy_var = 0
  if ( associated(src_sink_condition%temperature) ) then
    energy_var = SRC_TEMPERATURE 
  else if ( associated(src_sink_condition%enthalpy) ) then
    energy_var = SRC_ENTHALPY
  end if

  flow_src_sink_type = src_sink_condition%rate%itype

 ! checks that qsrc(liquid_phase) and qsrc(oil_phase) 
 ! do not have different signs
  if ( (qsrc(option%liquid_phase)>0.0d0 .and. qsrc(option%oil_phase)<0.d0).or.&
      (qsrc(option%liquid_phase)<0.0d0 .and. qsrc(option%oil_phase)>0.d0)  & 
    ) then
    option%io_buffer = "TOilImsSrcSink error: " // &
      "src(wat) and src(oil) with opposite sign"
    call printErrMsg(option)
  end if

  ! approximates BHP with local pressure
  ! to compute BHP we need to solve an IPR equation
  if ( ( (flow_src_sink_type == VOLUMETRIC_RATE_SS) .or. &
         ( associated(src_sink_condition%temperature) ) &
       ) .and. &
       (  (qsrc(option%liquid_phase) > 0.d0).or. &
         (qsrc(option%oil_phase) > 0.d0) &
       ) & 
     ) then  
    cell_pressure = &
        maxval(toil_auxvar%pres(option%liquid_phase:option%oil_phase))
  end if

  ! if enthalpy is used to define enthelpy or energy rate is used  
  ! approximate bottom hole temperature (BHT) with local temp
  if ( energy_var == SRC_TEMPERATURE) then
    temperature = src_sink_condition%temperature%dataset%rarray(1)
  else   
    temperature = toil_auxvar%temp
  end if


  Res = 0.d0
  do iphase = 1, option%nphase
    qsrc_mol = 0.d0
    if ( qsrc(iphase) > 0.d0) then 
      select case(iphase)
        case(LIQUID_PHASE)
          call EOSWaterDensity(temperature,cell_pressure,den_kg,den,ierr)
        case(TOIL_IMS_OIL_PHASE)
            call EOSOilDensity(temperature,cell_pressure,den,ierr)
      end select 
    else
      den = toil_auxvar%den(iphase)
    end if

    select case(flow_src_sink_type)
      ! injection and production 
      case(MASS_RATE_SS)
        qsrc_mol = qsrc(iphase)/toil_ims_fmw_comp(iphase) ! kg/sec -> kmol/sec
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_mol = qsrc(iphase)/toil_ims_fmw_comp(iphase)*scale 
      case(VOLUMETRIC_RATE_SS)  ! assume local density for now 
                  ! qsrc(iphase) = m^3/sec  
        qsrc_mol = qsrc(iphase)*den ! den = kmol/m^3 
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
        ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_mol = qsrc(iphase)* den * scale
        !qsrc_mol = qsrc(iphase)*gen_auxvar%den(iphase)*scale 
    end select
    ss_flow_vol_flux(iphase) = qsrc_mol/ den
    Res(iphase) = qsrc_mol
  enddo

  ! when using scaled src/sinks, the rates (marr or vol) scaling 
  ! at this point the scale factor is already included in Res(iphase)

  ! Res(option%energy_id), energy units: MJ/sec

  if ( associated(src_sink_condition%temperature) .or. &
      associated(src_sink_condition%enthalpy) &
     ) then
    ! if injection compute local pressure that will be used as BHP
    ! approximation used to overcome the solution of an IPR
    !if ( qsrc(option%liquid_phase)>0.d0 .or. 
    !    qsrc(option%oil_phase)>0.d0 ) then
    !  cell_pressure = &
    !      maxval(toil_auxvar%pres(option%liquid_phase:option%oil_phase))
    !end if
    ! water injection 
    if (qsrc(option%liquid_phase) > 0.d0) then !implies qsrc(option%oil_phase)>=0
      if ( energy_var == SRC_TEMPERATURE ) then
        call EOSWaterDensity(src_sink_condition%temperature% &
                             dataset%rarray(1), cell_pressure, &
                             den_kg,den,ierr)
        call EOSWaterEnthalpy(src_sink_condition%temperature% &
                              dataset%rarray(1), cell_pressure, &
                              enthalpy,ierr)
        ! enthalpy = [J/kmol]
      else if ( energy_var == SRC_ENTHALPY ) then
        !input as J/kg
        enthalpy = src_sink_condition%enthalpy% &
                       dataset%rarray(option%liquid_phase)
                     ! J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * toil_ims_fmw_comp(option%liquid_phase) 
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! water component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * enthalpy
    end if
    ! oil injection 
    if (qsrc(option%oil_phase) > 0.d0) then !implies qsrc(option%liquid_phase)>=0
      if ( energy_var == SRC_TEMPERATURE ) then
        call EOSOilEnthalpy(src_sink_condition%temperature%dataset%rarray(1), &
                            cell_pressure, enthalpy, ierr)
        ! enthalpy = [J/kmol] 
      else if ( energy_var == SRC_ENTHALPY ) then
        enthalpy = src_sink_condition%enthalpy% &
                     dataset%rarray(option%oil_phase)
                      !J/kg * kg/kmol = J/kmol  
        enthalpy = enthalpy * toil_ims_fmw_comp(option%oil_phase)        
      end if
      enthalpy = enthalpy * 1.d-6 ! J/kmol -> whatever units
      ! enthalpy units: MJ/kmol ! oil component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * enthalpy
    end if
    ! water energy extraction due to water production
    if (qsrc(option%liquid_phase) < 0.d0) then !implies qsrc(option%oil_phase)<=0
      ! auxvar enthalpy units: MJ/kmol ! water component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%liquid_phase) * &
                              toil_auxvar%H(option%liquid_phase)
    end if
    !oil energy extraction due to oil production 
    if (qsrc(option%oil_phase) < 0.d0) then !implies qsrc(option%liquid_phase)<=0
      ! auxvar enthalpy units: MJ/kmol ! water component mass                     
      Res(option%energy_id) = Res(option%energy_id) + &
                              Res(option%oil_phase) * &
                              toil_auxvar%H(option%oil_phase)
    end if

  else !if not temp or enthalpy are given
    ! if energy rate is given, loaded in qsrc(3) in MJ/sec 
    Res(option%energy_id) = qsrc(THREE_INTEGER)* scale ! MJ/s
  end if


  nullify(qsrc)      
  
end subroutine TOilImsSrcSink

! ************************************************************************** !
subroutine TOilImsDerivativePassTest(toil_auxvar,option)

  use Option_module
  use AuxVars_Flow_module
  
  implicit none

  type(auxvar_flow_type) :: toil_auxvar(0:)
  !class(auxvar_flow_type) :: toil_auxvar(0:)
  type(option_type) :: option

  write(*,"('acc sat01 derivTest = ',e10.4)") toil_auxvar(0)%den(2)
  write(*,"('acc sat11 derivTest = ',e10.4)") toil_auxvar(1)%den(2) 
  write(*,"('acc sat21 derivTest = ',e10.4)") toil_auxvar(2)%den(2) 
  write(*,"('acc sat31 derivTest = ',e10.4)") toil_auxvar(3)%den(2) 
     

end subroutine TOilImsDerivativePassTest

! ************************************************************************** !

subroutine TOilImsAccumDerivative(toil_auxvar,material_auxvar, &
                                  soil_heat_capacity,option,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/06/15
  ! 

  use Option_module
  use Saturation_Function_module
  use Material_Aux_class
  
  implicit none

  !type(toil_ims_auxvar_type) :: toil_auxvar(0:)
  !class(auxvar_toil_ims_type) :: toil_auxvar(0:)
  type(auxvar_toil_ims_type) :: toil_auxvar(0:)
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: soil_heat_capacity
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  !print *, 'ToilImsAccumDerivative'

  !write(*,"('acc sat01 derivB = ',e10.4)"), toil_auxvar(0)%den(2)
  !write(*,"('acc sat11 derivB = ',e10.4)"), toil_auxvar(1)%den(2) 
  !write(*,"('acc sat21 derivB = ',e10.4)"), toil_auxvar(2)%den(2) 
  !write(*,"('acc sat31 derivB = ',e10.4)"), toil_auxvar(3)%den(2) 

  call TOilImsAccumulation(toil_auxvar(ZERO_INTEGER), &
                           material_auxvar,soil_heat_capacity,option,Res)

  !write(*,"('acc sat01 derivM = ',e10.4)"), toil_auxvar(0)%den(2)
  !write(*,"('acc sat11 derivM = ',e10.4)"), toil_auxvar(1)%den(2) 
  !write(*,"('acc sat21 derivM = ',e10.4)"), toil_auxvar(2)%den(2) 
  !write(*,"('acc sat31 derivM = ',e10.4)"), toil_auxvar(3)%den(2) 

  do idof = 1, option%nflowdof
    call TOilImsAccumulation(toil_auxvar(idof), &
                           material_auxvar,soil_heat_capacity,option,res_pert)
    do irow = 1, option%nflowdof
      J(irow,idof) = (res_pert(irow)-res(irow))/toil_auxvar(idof)%pert
      !print *, irow, idof, J(irow,idof), toil_auxvar(idof)%pert
    enddo !irow
  enddo ! idof


  !write(*,"('acc sat01 deriv = ',e10.4)"), toil_auxvar(0)%den(2)
  !write(*,"('acc sat11 deriv = ',e10.4)"), toil_auxvar(1)%den(2) 
  !write(*,"('acc sat21 deriv = ',e10.4)"), toil_auxvar(2)%den(2) 
  !write(*,"('acc sat31 deriv = ',e10.4)"), toil_auxvar(3)%den(2) 

!  write(*,"('acc sat derivative = ',(4(e10.4,1x)))"), &
!        toil_auxvar(0:3)%sat(1)

  if (toil_ims_isothermal) then
    J(TOIL_IMS_ENERGY_EQUATION_INDEX,:) = 0.d0
    J(:,TOIL_IMS_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then
!    write(debug_unit,'(a,10es24.15)') 'accum deriv:', J
!  endif
!#endif

end subroutine TOilImsAccumDerivative

! ************************************************************************** !

subroutine ToilImsFluxDerivative(toil_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 sir_up, &
                                 thermal_conductivity_up, &
                                 toil_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 sir_dn, &
                                 thermal_conductivity_dn, &
                                 area, dist, &
                                 toil_parameter, &
                                 option,Jup,Jdn)
  ! 
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/06/15
  ! 
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  !type(toil_ims_auxvar_type) :: toil_auxvar_up(0:), toil_auxvar_dn(0:)
  !class(auxvar_toil_ims_type) :: toil_auxvar_up(0:), toil_auxvar_dn(0:)
  type(auxvar_toil_ims_type) :: toil_auxvar_up(0:), toil_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: thermal_conductivity_dn(2)
  PetscReal :: thermal_conductivity_up(2)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(toil_ims_parameter_type) :: toil_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jup = 0.d0
  Jdn = 0.d0
  
  !geh:print *, 'ToilImsFluxDerivative'
  option%iflag = -2
  call ToilImsFlux(toil_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up,sir_up, &
                   thermal_conductivity_up, &
                   toil_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn,sir_dn, &
                   thermal_conductivity_dn, &
                   area,dist,toil_parameter, &
                   option,v_darcy,res)
                           
  ! upgradient derivatives
  do idof = 1, option%nflowdof
    call ToilImsFlux(toil_auxvar_up(idof),global_auxvar_up, &
                     material_auxvar_up,sir_up, &
                     thermal_conductivity_up, &
                     toil_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn,sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,toil_parameter, &
                     option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jup(irow,idof) = (res_pert(irow)-res(irow))/toil_auxvar_up(idof)%pert
      !print *, 'up: ', irow, idof, Jup(irow,idof), toil_auxvar_up(idof)%pert
    enddo !irow
  enddo ! idof

  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call ToilImsFlux(toil_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                     material_auxvar_up,sir_up, &
                     thermal_conductivity_up, &
                     toil_auxvar_dn(idof),global_auxvar_dn, &
                     material_auxvar_dn,sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,toil_parameter, &
                     option,v_darcy,res_pert)
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/toil_auxvar_dn(idof)%pert
!geh:print *, 'dn: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (toil_ims_isothermal) then
    Jup(TOIL_IMS_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jup(:,TOIL_IMS_ENERGY_EQUATION_INDEX) = 0.d0
    Jdn(TOIL_IMS_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,TOIL_IMS_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then
!    write(debug_unit,'(a,20es24.15)') 'flux deriv:', Jup, Jdn
!  endif
!#endif
  
end subroutine ToilImsFluxDerivative

! ************************************************************************** !

subroutine ToilImsBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   toil_auxvar_up, &
                                   global_auxvar_up, &
                                   toil_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   sir_dn, &
                                   thermal_conductivity_dn, &
                                   area,dist,toil_parameter, &
                                   option,Jdn)
  ! 
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/06/15
  ! 

  use Option_module 
  use Material_Aux_class
  
  implicit none

  PetscReal :: auxvars(:) ! from aux_real_var array
  !type(toil_ims_auxvar_type) :: toil_auxvar_up, toil_auxvar_dn(0:) 
  !class(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn(0:)
  type(auxvar_toil_ims_type) :: toil_auxvar_up, toil_auxvar_dn(0:) 
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn(:)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(toil_ims_parameter_type) :: toil_parameter
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscInt :: ibndtype(1:option%nflowdof)
  PetscInt :: auxvar_mapping(TOIL_IMS_MAX_INDEX)
  PetscReal :: thermal_conductivity_dn(2)

  PetscReal :: v_darcy(option%nphase)
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscInt :: idof, irow

  Jdn = 0.d0
!geh:print *, 'GeneralBCFluxDerivative'

  option%iflag = -2
  call ToilImsBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     toil_auxvar_up,global_auxvar_up, &
                     toil_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     sir_dn, &
                     thermal_conductivity_dn, &
                     area,dist,toil_parameter, &
                     option,v_darcy,res)                     
  ! downgradient derivatives
  do idof = 1, option%nflowdof
    call ToilImsBCFlux(ibndtype,auxvar_mapping,auxvars, &
                       toil_auxvar_up,global_auxvar_up, &
                       toil_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       sir_dn, &
                       thermal_conductivity_dn, &
                       area,dist,toil_parameter, &
                       option,v_darcy,res_pert)   
    do irow = 1, option%nflowdof
      Jdn(irow,idof) = (res_pert(irow)-res(irow))/toil_auxvar_dn(idof)%pert
!print *, 'bc: ', irow, idof, Jdn(irow,idof), gen_auxvar_dn(idof)%pert
    enddo !irow
  enddo ! idof

  if (toil_ims_isothermal) then
    Jdn(TOIL_IMS_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jdn(:,TOIL_IMS_ENERGY_EQUATION_INDEX) = 0.d0
  endif
  
 
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then
!    write(debug_unit,'(a,10es24.15)') 'bc flux deriv:', Jdn
!  endif
!#endif
  
end subroutine ToilImsBCFluxDerivative


! ************************************************************************** !

subroutine ToilImsSrcSinkDerivative(option,src_sink_condition, toil_auxvar, &
                                    global_auxvar,scale,Jac)
  ! 
  ! Computes the source/sink terms for the residual
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/06/15
  ! 

  use Option_module
  use Condition_module

  implicit none

  type(option_type) :: option
  type(flow_toil_ims_condition_type), pointer :: src_sink_condition
  !type(toil_ims_auxvar_type) :: toil_auxvar(0:)
  !class(auxvar_toil_ims_type) :: toil_auxvar(0:)
  type(auxvar_toil_ims_type) :: toil_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  PetscReal :: scale
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  option%iflag = -3

  call TOilImsSrcSink(option,src_sink_condition,toil_auxvar(ZERO_INTEGER), &
                          global_auxvar,dummy_real,scale,Res)

  ! downgradient derivatives
  do idof = 1, option%nflowdof

    call TOilImsSrcSink(option,src_sink_condition,toil_auxvar(idof), &
                        global_auxvar,dummy_real,scale,res_pert)
  
    do irow = 1, option%nflowdof
      Jac(irow,idof) = (res_pert(irow)-res(irow))/toil_auxvar(idof)%pert
    enddo !irow
  enddo ! idof
  
  if (toil_ims_isothermal) then
    Jac(TOIL_IMS_ENERGY_EQUATION_INDEX,:) = 0.d0
    Jac(:,TOIL_IMS_ENERGY_EQUATION_INDEX) = 0.d0
  endif
   
!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then
!    write(debug_unit,'(a,20es24.15)') 'src/sink deriv:', Jac
!  endif
!#endif

end subroutine ToilImsSrcSinkDerivative

! ************************************************************************** !

subroutine TOilImsResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/05/15
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module

  use Connection_module
  use Grid_module
  use Coupler_module  
  use Debug_module
  use Material_Aux_class

!#define DEBUG_WITH_TECPLOT
#ifdef DEBUG_WITH_TECPLOT
  use Output_Tecplot_module
#endif

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  Mat, parameter :: null_mat = 0
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  
  type(toil_ims_parameter_type), pointer :: toil_parameter

  ! new auxvar deta strucutre
  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:), toil_auxvars_bc(:)

  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn

  !PetscInt :: iphase
  
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)

  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn
  PetscInt, save :: iplot = 0

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word

  PetscInt :: icap_up, icap_dn
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)
 
  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  !New auxvar data structure
  !toil_auxvars => patch%aux%TOil_ims%auxvars
  !toil_auxvars_bc => patch%aux%TOil_ims%auxvars_bc
  !toil_parameter => patch%aux%Toil_ims%parameter
  toil_parameter => patch%aux%TOil_ims%parameter

  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  
  ! Communication -----------------------------------------
  ! These 3 must be called before GeneralUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  
  call TOilImsUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%TOil_ims%auxvars_up_to_date = PETSC_FALSE 

  ! always assume variables have been swapped; therefore, must copy back
  ! PO check when copied at the end of iteration - this copy might not
  ! be needed 
  call VecLockPop(xx,ierr); CHKERRQ(ierr) !unlock vector from writing
  call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                   NFLOWDOF)
  call VecLockPush(xx,ierr); CHKERRQ(ierr) ! block vector from writing 

  if (option%compute_mass_balance_new) then
    call TOilImsZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  r_p = -accum_p

  
  !Heeho dynamically update p+1 accumulation term
  if (toil_ims_tough2_conv_criteria) then
    call VecGetArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
  ! accumulation at t(k+1)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call TOilImsAccumulation( &
                          patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id), &
                          material_auxvars(ghosted_id), &
                          material_parameter%soil_heat_capacity(imat), &
                          option,Res) 
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    
    !TOUGH2 conv. creteria: update p+1 accumulation term
    if (toil_ims_tough2_conv_criteria) then
      accum_p2(local_start:local_end) = Res(:)
    endif
    
  enddo

  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  !TOUGH2 conv. creteria: update p+1 accumulation term
  if (toil_ims_tough2_conv_criteria) then
    call VecRestoreArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      imat_up = patch%imat(ghosted_id_up) 
      imat_dn = patch%imat(ghosted_id_dn) 
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)

      call TOilImsFlux(patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     material_parameter%soil_residual_saturation(:,icap_up), &
                     material_parameter%soil_thermal_conductivity(:,imat_up), &
                     patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     material_parameter%soil_residual_saturation(:,icap_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     toil_parameter,option,v_darcy,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy
      if (associated(patch%internal_flow_fluxes)) then
        patch%internal_flow_fluxes(:,sum_connection) = Res(:)
      endif
      
      if (local_id_up > 0) then
        local_end = local_id_up * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) + Res(:)
      endif
         
      if (local_id_dn > 0) then
        local_end = local_id_dn * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) - Res(:)
      endif
    enddo

    cur_connection_set => cur_connection_set%next
  enddo    

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id)

#ifdef WELL_DEBUG
  write(*,*) 'BC gh = ', ghosted_id
  write(*,"('BC cell press = ',e26.20)") &
    patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id)%pres(option%oil_phase)
      !this%flow_auxvars(dof,ghosted_id)%pres(i_ph)
  write(*,"('BC press = ',e26.20)") &
    patch%aux%TOil_ims%auxvars_bc(sum_connection)%pres(option%oil_phase)
  write(*,"('BC delta press = ',e26.20)") &
    patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id)%pres(option%oil_phase) - &
    patch%aux%TOil_ims%auxvars_bc(sum_connection)%pres(option%oil_phase)
#endif

      call TOilImsBCFlux(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     patch%aux%TOil_ims%auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     material_parameter%soil_residual_saturation(:,icap_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     toil_parameter,option, &
                     v_darcy,Res)

      patch%boundary_velocities(:,sum_connection) = v_darcy
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2) ! one-component phase molar fluxes 
      endif

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      r_p(local_start:local_end)= r_p(local_start:local_end) - Res(:)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set

#ifdef WELL_DEBUG
  print *,'my rank = ',option%myrank,' just before ExplUpdate'
#endif

#ifdef WELL_CLASS
    if ( associated(source_sink%well) ) then
      if (cur_connection_set%num_connections > 0 ) then
        call source_sink%well%ExplUpdate(grid,option)
      end if
    end if 
#endif

#ifdef WELL_DEBUG
  print *,'my rank = ',option%myrank,' just after ExplUpdate'
#endif
 
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1

      ! if the src/sink is not scaled, flow_aux_real_var is not allocated
      ! time varying rate loaded in flow_condition%toil_ims%rate%dataset%rarray(:)  
      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      !if ( associated(source_sink%well) ) then
         ! use the if well to decide if to call TOilImsSrcSink or the WellRes
         !call source_sink%well%PrintMsg(); 
      !end if
#ifdef WELL_CLASS
      if ( associated(source_sink%well) ) then
        call source_sink%well%ExplRes(iconn,ss_flow_vol_flux, &
                       toil_ims_isothermal,ghosted_id,ZERO_INTEGER,option,Res)
      else
#endif
        call TOilImsSrcSink(option,source_sink%flow_condition%toil_ims, &
                           patch%aux%TOil_ims%auxvars(ZERO_INTEGER,ghosted_id), &
                           global_auxvars(ghosted_id),ss_flow_vol_flux, &
                           scale,Res)
#ifdef WELL_CLASS
      end if
#endif

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif      
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif      
      if (option%compute_mass_balance_new) then
        ! src/sinks contribution
        global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

    enddo 
!#ifdef WELL_DEBUG    
    ! for debugging
#ifdef WELL_CLASS
    if ( associated(source_sink%well) ) then
      if (cur_connection_set%num_connections > 0 ) then
        call source_sink%well%DataOutput(grid,source_sink%name,option)
      end if
    end if
#endif
!#endif
    source_sink => source_sink%next
  enddo

  if (patch%aux%TOil_ims%inactive_cells_exist) then
    do i=1,patch%aux%TOil_ims%n_inactive_rows
      r_p(patch%aux%TOil_ims%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  !do not use sandbox
  !call GeneralSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
  !                      gen_auxvars,option)

  !if (Initialized(toil_ims_debug_cell_id)) then
  !  call VecGetArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
  !  do local_id = general_debug_cell_id-1, general_debug_cell_id+1
  !    write(*,'(''  residual   : '',i2,10es12.4)') local_id, &
  !      r_p((local_id-1)*option%nflowdof+1:(local_id-1)*option%nflowdof+2), &
  !      r_p(local_id*option%nflowdof)*1.d6
  !  enddo
  !  call VecRestoreArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
  !endif
  
  if (toil_ims_isothermal) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+TOIL_IMS_ENERGY_EQUATION_INDEX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif

  
  if (realization%debug%vecview_residual) then
    string = 'Gresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'Gxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  
end subroutine TOilImsResidual

! ************************************************************************** !

! ************************************************************************** !

subroutine TOilImsJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian for TOilIms Mode
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/05/15
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class

  use AuxVars_Flow_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscInt :: icap_up,icap_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id
  PetscInt :: irow
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = 0
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  PetscInt, pointer :: zeros(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(material_parameter_type), pointer :: material_parameter
  type(toil_ims_parameter_type), pointer :: toil_parameter
  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:), toil_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)

  class(auxvar_flow_type), pointer :: toil_auxvar_p(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  !New auxvar data structure
  !toil_auxvars => patch%aux%TOil_ims%auxvars
  !toil_auxvars_bc => patch%aux%TOil_ims%auxvars_bc
  !toil_parameter => patch%aux%TOil_ims%parameter
  toil_parameter => patch%aux%TOil_ims%parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  ! Perturb aux vars
  do ghosted_id = 1, grid%ngmax  ! For each local node do...
    if (patch%imat(ghosted_id) <= 0) cycle

    call TOilImsAuxVarPerturb(patch%aux%TOil_ims%auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                               patch%sat_func_id(ghosted_id))%ptr, &
                              ghosted_id,option)
!#ifdef WELL_DEBUG    
!  write(*,"('after perturb = ',(4(e10.4,1x)))"), patch%aux%TOil_ims%auxvars(0:3,ghosted_id)%pert
!#endif 
    ! alternative way to call TOilImsAuxVarPerturb using a type-bound procedure 
    !call patch%aux%TOil_ims%Perturb(ghosted_id, &
    !                       global_auxvars(ghosted_id), &
    !                       material_auxvars(ghosted_id), &
    !                       patch%characteristic_curves_array( &
    !                       patch%sat_func_id(ghosted_id))%ptr, &
    !                       ghosted_id,option)

  enddo

#ifdef WELL_DEBUG    
  write(*,"('AP p011 before = ',e10.4)") patch%aux%TOil_ims%auxvars(0,1)%pres(1)
  write(*,"('AP p111 before = ',e10.4)") patch%aux%TOil_ims%auxvars(1,1)%pres(1)
  write(*,"('AP p211 before = ',e10.4)") patch%aux%TOil_ims%auxvars(2,1)%pres(1)
  write(*,"('AP p311 before = ',e10.4)") patch%aux%TOil_ims%auxvars(3,1)%pres(1)
#endif 

  
!#ifdef DEBUG_GENERAL_LOCAL
!  call GeneralOutputAuxVars(gen_auxvars,global_auxvars,option)
!#endif 

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    !write(*,"('acc sat before = ',(4(e10.4,1x)))"), &
    !     patch%aux%TOil_ims%auxvars(0:3,ghosted_id)%sat(1)
    !write(*,"('acc sat01 before = ',e10.4)"), patch%aux%TOil_ims%auxvars(0,ghosted_id)%den(2)
    !write(*,"('acc sat11 before = ',e10.4)"), patch%aux%TOil_ims%auxvars(1,ghosted_id)%den(2) 
    !write(*,"('acc sat21 before = ',e10.4)"), patch%aux%TOil_ims%auxvars(2,ghosted_id)%den(2) 
    !write(*,"('acc sat31 before = ',e10.4)"), patch%aux%TOil_ims%auxvars(3,ghosted_id)%den(2) 

    call TOilImsAccumDerivative(patch%aux%TOil_ims%auxvars(:,ghosted_id), &
                                material_auxvars(ghosted_id), &
                                material_parameter%soil_heat_capacity(imat), & 
                                option,Jup)
    !toil_auxvar_p => patch%aux%TOil_ims%auxvars(:,ghosted_id)
    !call TOilImsDerivativePassTest(patch%aux%TOil_ims%auxvars(:,ghosted_id),option)
    !call TOilImsDerivativePassTest(toil_auxvar_p,option)

    !write(*,"('acc sat after = ',(4(e10.4,1x)))"), &
    !     patch%aux%TOil_ims%auxvars(0:3,ghosted_id)%sat(1)
    !write(*,"('acc sat01 after = ',e10.4)"), patch%aux%TOil_ims%auxvars(0,ghosted_id)%den(2)
    !write(*,"('acc sat11 after = ',e10.4)"), patch%aux%TOil_ims%auxvars(1,ghosted_id)%den(2) 
    !write(*,"('acc sat21 after = ',e10.4)"), patch%aux%TOil_ims%auxvars(2,ghosted_id)%den(2) 
    !write(*,"('acc sat31 after = ',e10.4)"), patch%aux%TOil_ims%auxvars(3,ghosted_id)%den(2) 

    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_accum'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif


  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      imat_up = patch%imat(ghosted_id_up)
      imat_dn = patch%imat(ghosted_id_dn)
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
      ! if issues in passing auxvars, pass the entire TOil_ims and 
      ! ghosted_id_up, ghosted_id_dn
      call TOilImsFluxDerivative(patch%aux%TOil_ims%auxvars(:,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       material_parameter%soil_residual_saturation(:,icap_up), &
                       material_parameter%soil_thermal_conductivity(:,imat_up), &
                       patch%aux%TOil_ims%auxvars(:,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       material_parameter%soil_residual_saturation(:,icap_dn), &
                       material_parameter%soil_thermal_conductivity(:,imat_dn), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       toil_parameter,option, &
                       Jup,Jdn)
     
      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_flux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icap_dn = patch%sat_func_id(ghosted_id)

      call TOilImsBCFluxDerivative(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     patch%aux%TOil_ims%auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     patch%aux%TOil_ims%auxvars(:,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     material_parameter%soil_residual_saturation(:,icap_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     toil_parameter,option, &
                     Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_bcflux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  ! Source/sinks
  source_sink => patch%source_sink_list%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
#ifdef WELL_CLASS      
      if (associated(source_sink%well) ) then
        !Jdn = 0.0d0
        call source_sink%well%ExplJDerivative(iconn,ghosted_id, &
                        toil_ims_isothermal,TOIL_IMS_ENERGY_EQUATION_INDEX, &
                         option,Jdn)
        Jdn = -Jdn  
        call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)

      else 
#endif
        !Jup = 0.d0
        call TOilImsSrcSinkDerivative(option, &
                          source_sink%flow_condition%toil_ims, &
                          patch%aux%TOil_ims%auxvars(:,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          scale,Jup)
        call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
#ifdef WELL_CLASS
      end if
#endif

    enddo
    source_sink => source_sink%next
  enddo
   
  ! SSSandBox not supported 
  !call GeneralSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
  !                      gen_auxvars,option)

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_srcsink'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out isothermal and inactive cells
  if (patch%aux%TOil_ims%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%TOil_ims%n_inactive_rows, &
                          patch%aux%TOil_ims%inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          ierr);CHKERRQ(ierr)
  endif

  if (toil_ims_isothermal) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => patch%aux%Toil_ims%row_zeroing_array
    ! zero energy residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        TOIL_IMS_ENERGY_EQUATION_INDEX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%matview_Jacobian) then
    string = 'Gjacobian'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option) 
  endif

!  call MatView(J,PETSC_VIEWER_STDOUT_WORLD,ierr)

!#if 0
!  imat = 1
!  if (imat == 1) then
!    call GeneralNumericalJacobianTest(xx,realization,J) 
!  endif
!#endif

!#ifdef DEBUG_GENERAL_FILEOUTPUT
!  if (debug_flag > 0) then
!    write(word,*) debug_timestep_count
!    string = 'jacobian_' // trim(adjustl(word))
!    write(word,*) debug_timestep_cut_count
!    string = trim(string) // '_' // trim(adjustl(word))
!    write(word,*) debug_iteration_count
!    string = trim(string) // '_' // trim(adjustl(word)) // '.out'
!    call PetscViewerASCIIOpen(realization%option%mycomm,trim(string), &
!                              viewer,ierr);CHKERRQ(ierr)
!    call MatView(J,viewer,ierr);CHKERRQ(ierr)
!    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
!    close(debug_unit)
!  endif
!#endif

end subroutine ToilImsJacobian

! ************************************************************************** !

subroutine TOilImsDestroy(realization)
  ! 
  ! Deallocates variables associated with TOilIms
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15
  ! 

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine TOilImsDestroy

! ************************************************************************** !

! ************************************************************************** !

function TOilImsAverageDensity(sat_up,sat_dn,density_up,density_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/28/15
  ! 

  implicit none

  PetscReal :: sat_up, sat_dn
  PetscReal :: density_up, density_dn

  PetscReal :: TOilImsAverageDensity

  if (sat_up < eps ) then
    TOilImsAverageDensity = density_dn
  else if (sat_dn < eps ) then 
    TOilImsAverageDensity = density_up
  else ! in here we could use an armonic average, 
       ! other idea sat weighted average but it needs truncation
    TOilImsAverageDensity = 0.5d0*(density_up+density_dn)
  end if

end function TOilImsAverageDensity

! ************************************************************************** !


end module TOilIms_module

