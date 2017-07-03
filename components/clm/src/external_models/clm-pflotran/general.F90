module General_module

  use General_Aux_module
  use General_Common_module
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
  
  public :: GeneralSetup, &
            GeneralInitializeTimestep, &
            GeneralUpdateSolution, &
            GeneralTimeCut,&
            GeneralUpdateAuxVars, &
            GeneralUpdateFixedAccum, &
            GeneralComputeMassBalance, &
            GeneralResidual, &
            GeneralJacobian, &
            GeneralGetTecplotHeader, &
            GeneralSetPlotVariables, &
            GeneralMapBCAuxVarsToGlobal, &
            GeneralDestroy

contains

! ************************************************************************** !

subroutine GeneralSetup(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Fluid_module
  use Material_Aux_class
  use Output_Aux_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: i, idof, count
  PetscBool :: error_found
  PetscInt :: flag(10)
                                                ! extra index for derivatives
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(general_auxvar_type), pointer :: gen_auxvars_bc(:)
  type(general_auxvar_type), pointer :: gen_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%General => GeneralAuxCreate(option)

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE
  if (minval(material_parameter%soil_residual_saturation(:,:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil residual saturation.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_heat_capacity(:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil heat capacity.'
    call printMsg(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_thermal_conductivity(:,:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil thermal conductivity.'
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
      option%io_buffer = 'ERROR: Non-initialized cell volume.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'ERROR: Non-initialized tortuosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%soil_particle_density < 0.d0 .and. &
        flag(4) == 0) then
      flag(4) = 1
      option%io_buffer = 'ERROR: Non-initialized soil particle density.'
      call printMsg(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'ERROR: Non-initialized permeability.'
      call printMsg(option)
      flag(5) = 1
    endif
  enddo
  
  if (error_found .or. maxval(flag) > 0) then
    option%io_buffer = 'Material property errors found in GeneralSetup.'
    call printErrMsg(option)
  endif
  
  ! allocate auxvar data structures for all grid cells  
  allocate(gen_auxvars(0:option%nflowdof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      call GeneralAuxVarInit(gen_auxvars(idof,ghosted_id), &
                             (general_analytical_derivatives .and. idof==0), &
                             option)
    enddo
  enddo
  patch%aux%General%auxvars => gen_auxvars
  patch%aux%General%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them 
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(gen_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_auxvars_bc(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%General%auxvars_bc => gen_auxvars_bc
  endif
  patch%aux%General%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(gen_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call GeneralAuxVarInit(gen_auxvars_ss(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%General%auxvars_ss => gen_auxvars_ss
  endif
  patch%aux%General%num_aux_ss = sum_connection

  ! create array for zeroing Jacobian entries if isothermal and/or no air
  allocate(patch%aux%General%row_zeroing_array(grid%nlmax))
  patch%aux%General%row_zeroing_array = 0
  
  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    patch%aux%General%general_parameter% &
      diffusion_coefficient(cur_fluid_property%phase_id) = &
        cur_fluid_property%diffusion_coefficient
    cur_fluid_property => cur_fluid_property%next
  enddo  
  ! check whether diffusion coefficients are initialized.
  if (Uninitialized(patch%aux%General%general_parameter% &
      diffusion_coefficient(LIQUID_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Liquid phase diffusion coefficient','')
    call printErrMsg(option)
  endif
  if (Uninitialized(patch%aux%General%general_parameter% &
      diffusion_coefficient(GAS_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Gas phase diffusion coefficient','')
    call printErrMsg(option)
  endif

  list => realization%output_option%output_snap_variable_list
  call GeneralSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call GeneralSetPlotVariables(realization,list)
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flag = 0
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = 0
  ! create new file
  open(debug_info_unit, file='debug_info.txt', action="write", &
       status="unknown")
  write(debug_info_unit,*) 'type timestep cut iteration'
  close(debug_info_unit)
#endif  

end subroutine GeneralSetup

! ************************************************************************** !

subroutine GeneralInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call GeneralUpdateFixedAccum(realization)
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_flag = 0
!  if (realization%option%time >= 35.6d0*3600d0*24.d0*365.d0 - 1.d-40) then
!  if (.false.) then
  if (.true.) then
    debug_iteration_count = 0
    debug_flag = 1
  endif
  debug_iteration_count = 0
#endif

end subroutine GeneralInitializeTimestep

! ************************************************************************** !

subroutine GeneralUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  gen_auxvars => patch%aux%General%auxvars  
  global_auxvars => patch%aux%Global%auxvars
  
  if (realization%option%compute_mass_balance_new) then
    call GeneralUpdateMassBalance(realization)
  endif
  
  ! update stored state
  do ghosted_id = 1, grid%ngmax
    gen_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS) = &
      global_auxvars(ghosted_id)%istate
  enddo
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_iteration_count = 0
  debug_timestep_cut_count = 0
  debug_timestep_count = debug_timestep_count + 1
#endif   
  
end subroutine GeneralUpdateSolution

! ************************************************************************** !

subroutine GeneralTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
  ! 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  gen_auxvars => patch%aux%General%auxvars

  ! restore stored state
  do ghosted_id = 1, grid%ngmax
    global_auxvars(ghosted_id)%istate = &
      gen_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS)
  enddo

#ifdef DEBUG_GENERAL_FILEOUTPUT
  debug_timestep_cut_count = debug_timestep_cut_count + 1
#endif 

  call GeneralInitializeTimestep(realization)  

end subroutine GeneralTimeCut

! ************************************************************************** !

subroutine GeneralNumericalJacobianTest(xx,realization,B)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/03/15
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_subsurface_type) :: realization
  Mat :: B

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  PetscReal :: derivative, perturbation
  PetscReal :: perturbation_tolerance = 1.d-6
  PetscInt, save :: icall = 0
  character(len=MAXWORDLENGTH) :: word, word2

  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  icall = icall + 1
  call VecDuplicate(xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res_pert,ierr);CHKERRQ(ierr)

  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof, &
                   grid%nlmax*option%nflowdof, &
                   ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(A,27,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)

  call VecZeroEntries(res,ierr);CHKERRQ(ierr)
  call GeneralResidual(PETSC_NULL_OBJECT,xx,res,realization,ierr)
#if 0
  word  = 'num_0.dat'
  call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
  call VecView(res,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call VecCopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      call VecZeroEntries(res_pert,ierr);CHKERRQ(ierr)
      call GeneralResidual(PETSC_NULL_OBJECT,xx_pert,res_pert,realization,ierr)
#if 0
      write(word,*) idof
      word  = 'num_' // trim(adjustl(word)) // '.dat'
      call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
      call VecView(res_pert,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
      call VecGetArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call MatSetValue(A,idof2-1,idof-1,derivative,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if 1
  write(word,*) icall
  word = 'numerical_jacobian-' // trim(adjustl(word)) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

!geh: uncomment to overwrite numerical Jacobian
!  call MatCopy(A,B,DIFFERENT_NONZERO_PATTERN,ierr)
  call MatDestroy(A,ierr);CHKERRQ(ierr)

  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)

end subroutine GeneralNumericalJacobianTest

! ************************************************************************** !

subroutine GeneralComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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
  type(general_auxvar_type), pointer :: general_auxvars(:,:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase, icomp
  PetscReal :: vol_phase

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  general_auxvars => patch%aux%General%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = &
        general_auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        general_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density
      do icomp = 1, option%nflowspec
        mass_balance(icomp,iphase) = mass_balance(icomp,iphase) + &
          general_auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)* &
          general_auxvars(ZERO_INTEGER,ghosted_id)%xmol(icomp,iphase) * &
          fmw_comp(icomp)*vol_phase
      enddo
    enddo
  enddo

end subroutine GeneralComputeMassBalance

! ************************************************************************** !

subroutine GeneralZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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

  do iconn = 1, patch%aux%General%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%General%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine GeneralZeroMassBalanceDelta

! ************************************************************************** !

subroutine GeneralUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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
  PetscInt :: icomp

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%General%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%General%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine GeneralUpdateMassBalance

! ************************************************************************** !

subroutine GeneralUpdateAuxVars(realization,update_state)
  ! 
  ! Updates the auxiliary variables associated with the General problem
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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
  use EOS_Water_module
  use Saturation_Function_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(general_auxvar_type), pointer :: gen_auxvars(:,:), gen_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end
  PetscInt :: iphasebc, iphase
  PetscInt :: offset
  PetscInt :: istate
  PetscReal :: gas_pressure, capillary_pressure, liquid_saturation
  PetscReal :: saturation_pressure, temperature
  PetscInt :: real_index, variable
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
!#define DEBUG_AUXVARS
#ifdef DEBUG_AUXVARS
  character(len=MAXWORDLENGTH) :: word
  PetscInt, save :: icall = 0
#endif
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  gen_auxvars => patch%aux%General%auxvars
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

#ifdef DEBUG_AUXVARS
  icall = icall + 1
  write(word,*) icall
  word = 'genaux' // trim(adjustl(word))
#endif
  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    ! GENERAL_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = GENERAL_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call GeneralAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       gen_auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%characteristic_curves_array( &
                         patch%sat_func_id(ghosted_id))%ptr, &
                       natural_id, &
                       option)
    if (update_state) then
      call GeneralAuxVarUpdateState(xx_loc_p(ghosted_start:ghosted_end), &
                                    gen_auxvars(ZERO_INTEGER,ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%sat_func_id(ghosted_id))%ptr, &
                                    natural_id, &  ! for debugging
                                    option)
    endif
#ifdef DEBUG_AUXVARS
!geh: for debugging
    call GeneralOutputAuxVars(gen_auxvars(0,ghosted_id), &
                              global_auxvars(ghosted_id),natural_id,word, &
                              PETSC_TRUE,option)
#endif
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(debug_unit,'(a,i5,i3,7es24.15)') 'auxvar:', natural_id, &
                        global_auxvars(ghosted_id)%istate, &
                        xx_loc_p(ghosted_start:ghosted_end)
  endif
#endif
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      !geh: negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id) 
      offset = (ghosted_id-1)*option%nflowdof
      if (patch%imat(ghosted_id) <= 0) cycle

      xxbc(:) = xx_loc_p(offset+1:offset+option%nflowdof)
      istate = boundary_condition%flow_aux_int_var(GENERAL_STATE_INDEX,iconn)
      if (istate == ANY_STATE) then
        istate = global_auxvars(ghosted_id)%istate
        select case(istate)
          case(LIQUID_STATE,GAS_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(DIRICHLET_BC,HYDROSTATIC_BC)
                  real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
              end select   
            enddo
          case(TWO_PHASE_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(HYDROSTATIC_BC)
                  real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                case(DIRICHLET_BC)
                  variable = dof_to_primary_variable(idof,istate)
                  select case(variable)
                    ! for gas pressure dof
                    case(GENERAL_GAS_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs gas pressure defined.'
                        call printErrMsg(option)
                      endif
                    ! for air pressure dof
                    case(GENERAL_AIR_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index == 0) then ! air pressure not found
                        ! if air pressure is not available, let's try temperature 
                        real_index = boundary_condition%flow_aux_mapping(GENERAL_TEMPERATURE_INDEX)
                        if (real_index /= 0) then
                          temperature = boundary_condition%flow_aux_real_var(real_index,iconn)
                          call EOSWaterSaturationPressure(temperature,saturation_pressure,ierr)
                          ! now verify whether gas pressure is provided through BC
                          if (boundary_condition%flow_bc_type(ONE_INTEGER) == NEUMANN_BC) then
                            gas_pressure = xxbc(ONE_INTEGER)
                          else
                            real_index = boundary_condition%flow_aux_mapping(GENERAL_GAS_PRESSURE_INDEX)
                            if (real_index /= 0) then
                              gas_pressure = boundary_condition%flow_aux_real_var(real_index,iconn)
                            else
                              option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                                trim(boundary_condition%flow_condition%name) // &
                                '" needs gas pressure defined to calculate air ' // &
                                'pressure from temperature.'
                              call printErrMsg(option)
                            endif
                          endif
                          xxbc(idof) = gas_pressure - saturation_pressure
                        else
                          option%io_buffer = 'Cannot find boundary constraint for air pressure.'
                          call printErrMsg(option)
                        endif
                      else
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      endif
                    ! for gas saturation dof
                    case(GENERAL_GAS_SATURATION_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
!geh: should be able to use the saturation within the cell
!                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
!                          trim(boundary_condition%flow_condition%name) // &
!                          '" needs saturation defined.'
!                        call printErrMsg(option)
                      endif
                    case(GENERAL_TEMPERATURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs temperature defined.'
                        call printErrMsg(option)
                      endif
                  end select
                case(NEUMANN_BC)
                case default
                  option%io_buffer = 'Unknown BC type in GeneralUpdateAuxVars().'
                  call printErrMsg(option)
              end select
            enddo  
        end select
      else
        ! we do this for all BCs; Neumann bcs will be set later
        do idof = 1, option%nflowdof
          real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
          if (real_index > 0) then
            xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
          else
            option%io_buffer = 'Error setting up boundary condition in GeneralUpdateAuxVars'
            call printErrMsg(option)
          endif
        enddo
      endif
          
      ! set this based on data given 
      global_auxvars_bc(sum_connection)%istate = istate
      ! GENERAL_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = GENERAL_UPDATE_FOR_BOUNDARY
      call GeneralAuxVarCompute(xxbc,gen_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
      ! update state and update aux var; this could result in two update to 
      ! the aux var as update state updates if the state changes
      call GeneralAuxVarUpdateState(xxbc,gen_auxvars_bc(sum_connection), &
                                    global_auxvars_bc(sum_connection), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%sat_func_id(ghosted_id))%ptr, &
                                    natural_id,option)
#ifdef DEBUG_GENERAL_FILEOUTPUT
      if (debug_flag > 0) then
        write(debug_unit,'(a,i5,i3,7es24.15)') 'bc_auxvar:', natural_id, &
                           global_auxvars_bc(ghosted_id)%istate, &
                            xxbc(:)
      endif
#endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%General%auxvars_up_to_date = PETSC_TRUE

end subroutine GeneralUpdateAuxVars

! ************************************************************************** !

subroutine GeneralUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/11
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
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  gen_auxvars => patch%aux%General%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  !Heeho initialize dynamic accumulation term for every p iteration step
  if (general_tough2_conv_criteria) then
    call VecGetArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! GENERAL_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = GENERAL_UPDATE_FOR_FIXED_ACCUM
    call GeneralAuxVarCompute(xx_p(local_start:local_end), &
                              gen_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              natural_id, &
                              option)
    call GeneralAccumulation(gen_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,accum_p(local_start:local_end), &
                             Jac_dummy,PETSC_FALSE, &
                             local_id == general_debug_cell_id) 
  enddo
  
  if (general_tough2_conv_criteria) then
    accum_p2 = accum_p
  endif
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
  !Heeho initialize dynamic accumulation term for every p iteration step
  if (general_tough2_conv_criteria) then
    call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  endif
  
end subroutine GeneralUpdateFixedAccum

! ************************************************************************** !

subroutine GeneralResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
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
  type(general_parameter_type), pointer :: general_parameter
  type(general_auxvar_type), pointer :: gen_auxvars(:,:), gen_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: iphase
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
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)
  
  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  gen_auxvars => patch%aux%General%auxvars
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  general_parameter => patch%aux%General%general_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  
#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    debug_iteration_count = debug_iteration_count + 1
    write(word,*) debug_timestep_count
    string = 'residual_debug_data_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    debug_filename = trim(string) // '_' // trim(adjustl(word)) // '.txt'
    open(debug_unit, file=debug_filename, action="write", status="unknown")
    open(debug_info_unit, file='debug_info.txt', action="write", &
         position="append", status="unknown")
    write(debug_info_unit,*) 'residual ', debug_timestep_count, &
      debug_timestep_cut_count, debug_iteration_count
    close(debug_info_unit)
  endif
#endif

  ! Communication -----------------------------------------
  ! These 3 must be called before GeneralUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  
                                             ! do update state
  call GeneralUpdateAuxVars(realization,PETSC_TRUE)

! for debugging a single grid cell
!  i = 6
!  call GeneralOutputAuxVars(gen_auxvars(0,i),global_auxvars(i),i,'genaux', &
!                            PETSC_TRUE,option)
#ifdef DEBUG_WITH_TECPLOT
! for debugging entire solution over a single SNES solve
  write(word,*) iplot
  iplot = iplot + 1
  realization%output_option%plot_name = 'general-ni-' // trim(adjustl(word))
  call OutputTecplotPoint(realization)
#endif

  ! override flags since they will soon be out of date
  patch%aux%General%auxvars_up_to_date = PETSC_FALSE 

  ! always assume variables have been swapped; therefore, must copy back
  call VecLockPop(xx,ierr); CHKERRQ(ierr)
  call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                   NFLOWDOF)
  call VecLockPush(xx,ierr); CHKERRQ(ierr)

  if (option%compute_mass_balance_new) then
    call GeneralZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  r_p = -accum_p

  
  !Heeho dynamically update p+1 accumulation term
  if (general_tough2_conv_criteria) then
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
    call GeneralAccumulation(gen_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,Res,Jac_dummy, &
                             general_analytical_derivatives, &
                             local_id == general_debug_cell_id) 
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    
    !Heeho dynamically update p+1 accumulation term
    if (general_tough2_conv_criteria) then
      accum_p2(local_start:local_end) = Res(:)
    endif
    
  enddo

  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  !Heeho dynamically update p+1 accumulation term
  if (general_tough2_conv_criteria) then
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

      call GeneralFlux(gen_auxvars(ZERO_INTEGER,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       material_parameter%soil_thermal_conductivity(:,imat_up), &
                       gen_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       material_parameter%soil_thermal_conductivity(:,imat_dn), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       general_parameter,option,v_darcy,Res, &
                       Jac_dummy,Jac_dummy, &
                       general_analytical_derivatives, &
                       (local_id_up == general_debug_cell_id .or. &
                        local_id_dn == general_debug_cell_id))

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

      call GeneralBCFlux(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     gen_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     gen_auxvars(ZERO_INTEGER,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     general_parameter,option, &
                     v_darcy,Res,Jac_dummy, &
                     general_analytical_derivatives, &
                     local_id == general_debug_cell_id)
      patch%boundary_velocities(:,sum_connection) = v_darcy
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
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
    
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      call GeneralSrcSink(option,source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                          source_sink%flow_condition%general%rate%itype, &
                          gen_auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          ss_flow_vol_flux, &
                          scale,Res,Jac_dummy, &
                          general_analytical_derivatives, &
                          local_id == general_debug_cell_id)

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif      
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif      
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%General%inactive_cells_exist) then
    do i=1,patch%aux%General%n_inactive_rows
      r_p(patch%aux%General%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  call GeneralSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                        gen_auxvars,option)

  ! Mass Transfer
  if (field%flow_mass_transfer /= 0) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif                      
                        
  if (Initialized(general_debug_cell_id)) then
    call VecGetArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
    do local_id = general_debug_cell_id-1, general_debug_cell_id+1
      write(*,'(''  residual   : '',i2,10es12.4)') local_id, &
        r_p((local_id-1)*option%nflowdof+1:(local_id-1)*option%nflowdof+2), &
        r_p(local_id*option%nflowdof)*1.d6
    enddo
    call VecRestoreArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  
  if (general_isothermal) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+GENERAL_ENERGY_EQUATION_INDEX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  if (general_no_air) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    ! zero energy residual
    do local_id = 1, grid%nlmax
      r_p((local_id-1)*option%nflowdof+GENERAL_GAS_EQUATION_INDEX) =  0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif  

#ifdef DEBUG_GENERAL_FILEOUTPUT
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'fixed residual:', local_id, &
      accum_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    write(debug_unit,'(a,i5,7es24.15)') 'residual:', local_id, &
      r_p((local_id-1)*option%nflowdof+1:local_id*option%nflowdof)
  enddo
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
#endif
  
  
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

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    close(debug_unit)
  endif
#endif
  
end subroutine GeneralResidual

! ************************************************************************** !

subroutine GeneralJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
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
  PetscInt :: local_id, ghosted_id, natural_id
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
  type(general_parameter_type), pointer :: general_parameter
  type(general_auxvar_type), pointer :: gen_auxvars(:,:), gen_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  gen_auxvars => patch%aux%General%auxvars
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  general_parameter => patch%aux%General%general_parameter
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

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(word,*) debug_timestep_count
    string = 'jacobian_debug_data_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    debug_filename = trim(string) // '_' // trim(adjustl(word)) // '.txt'
    open(debug_unit, file=debug_filename, action="write", status="unknown")
    open(debug_info_unit, file='debug_info.txt', action="write", &
         position="append", status="unknown")
    write(debug_info_unit,*) 'jacobian ', debug_timestep_count, &
      debug_timestep_cut_count, debug_iteration_count
    close(debug_info_unit)
  endif
#endif

  ! Perturb aux vars
  do ghosted_id = 1, grid%ngmax  ! For each local node do...
    if (patch%imat(ghosted_id) <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    call GeneralAuxVarPerturb(gen_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              natural_id,option)
  enddo
  
#ifdef DEBUG_GENERAL_LOCAL
  call GeneralOutputAuxVars(gen_auxvars,global_auxvars,option)
#endif 

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call GeneralAccumDerivative(gen_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              material_parameter%soil_heat_capacity(imat), &
                              option, &
                              Jup) 
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
                              
      call GeneralFluxDerivative(gen_auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     material_parameter%soil_thermal_conductivity(:,imat_up), &
                     gen_auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     general_parameter,option,&
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

      call GeneralBCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      gen_auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      gen_auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      material_parameter%soil_thermal_conductivity(:,imat_dn), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      general_parameter,option, &
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
      
      Jup = 0.d0
      call GeneralSrcSinkDerivative(option, &
                        source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                        source_sink%flow_condition%general%rate%itype, &
                        gen_auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo
  
  call GeneralSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
                        gen_auxvars,option)

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
  if (patch%aux%General%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%General%n_inactive_rows, &
                          patch%aux%General%inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          ierr);CHKERRQ(ierr)
  endif

  if (general_isothermal) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => patch%aux%General%row_zeroing_array
    ! zero energy residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        GENERAL_ENERGY_EQUATION_INDEX - 1 ! zero-based
    enddo
    call MatZeroRowsLocal(A,grid%nlmax,zeros,qsrc,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
  endif

  if (general_no_air) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    zeros => patch%aux%General%row_zeroing_array
    ! zero gas component mass balance residual
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      zeros(local_id) = (ghosted_id-1)*option%nflowdof+ &
                        GENERAL_GAS_EQUATION_INDEX - 1 ! zero-based
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

#if 0
  imat = 1
  if (imat == 1) then
    call GeneralNumericalJacobianTest(xx,realization,J) 
  endif
#endif

#ifdef DEBUG_GENERAL_FILEOUTPUT
  if (debug_flag > 0) then
    write(word,*) debug_timestep_count
    string = 'jacobian_' // trim(adjustl(word))
    write(word,*) debug_timestep_cut_count
    string = trim(string) // '_' // trim(adjustl(word))
    write(word,*) debug_iteration_count
    string = trim(string) // '_' // trim(adjustl(word)) // '.out'
    call PetscViewerASCIIOpen(realization%option%mycomm,trim(string), &
                              viewer,ierr);CHKERRQ(ierr)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    close(debug_unit)
  endif
#endif

end subroutine GeneralJacobian

! ************************************************************************** !

function GeneralGetTecplotHeader(realization,icolumn)
  ! 
  ! Returns General Lite contribution to
  ! Tecplot file header
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 
  
  use Realization_Subsurface_class
  use Option_module
  use Field_module
    
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: GeneralGetTecplotHeader
  type(realization_subsurface_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i

  option => realization%option
  field => realization%field
  
  string = ''
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-T [C]"'')') icolumn
  else
    write(string2,'('',"T [C]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-State"'')') icolumn
  else
    write(string2,'('',"State"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sat(l)"'')') icolumn
  else
    write(string2,'('',"Sat(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sat(g)"'')') icolumn
  else
    write(string2,'('',"Sat(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Rho(l)"'')') icolumn
  else
    write(string2,'('',"Rho(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Rho(g)"'')') icolumn
  else
    write(string2,'('',"Rho(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(l)"'')') icolumn
  else
    write(string2,'('',"U(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(g)"'')') icolumn
  else
    write(string2,'('',"U(g)"'')')
  endif
  string = trim(string) // trim(string2)
  
  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i2,'')"'')') icolumn, i
    else
      write(string2,'('',"Xl('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xg('',i2,'')"'')') icolumn, i
    else
      write(string2,'('',"Xg('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo
 
  GeneralGetTecplotHeader = string

end function GeneralGetTecplotHeader

! ************************************************************************** !

subroutine GeneralSetPlotVariables(realization,list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/13
  ! 
  
  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(realization_subsurface_type) :: realization
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

  name = 'Gas Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               GAS_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)
  
  name = 'Gas Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               GAS_SATURATION)
  
  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)
  
  name = 'Gas Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_DENSITY)
  
  name = 'X_g^l'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION, &
                               realization%option%air_id)
  
  name = 'X_l^l'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOLE_FRACTION, &
                               realization%option%water_id)
  
  name = 'X_g^g'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOLE_FRACTION, &
                               realization%option%air_id)
  
  name = 'X_l^g'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOLE_FRACTION, &
                               realization%option%water_id)
  
  name = 'Liquid Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)
  
  name = 'Gas Energy'
  units = 'MJ/kmol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_ENERGY)
  
  name = 'Thermodynamic State'
  units = ''
  output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE,units,STATE)
  output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList(list,output_variable)   
  
end subroutine GeneralSetPlotVariables

! ************************************************************************** !

function GeneralAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn,dden_up,dden_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/14
  ! 

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)
  PetscReal :: dden_up, dden_dn

  PetscReal :: GeneralAverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0
  if (iphase == LIQUID_PHASE) then
    if (istate_up == GAS_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0
    else if (istate_dn == GAS_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == LIQUID_STATE) then
      GeneralAverageDensity = density_dn(iphase)
      dden_dn = 1.d0      
    else if (istate_dn == LIQUID_STATE) then
      GeneralAverageDensity = density_up(iphase)
      dden_up = 1.d0      
    else
      GeneralAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0      
    endif
  endif

end function GeneralAverageDensity

! ************************************************************************** !

subroutine GeneralSSSandbox(residual,Jacobian,compute_derivative, &
                            grid,material_auxvars,general_auxvars,option)
  ! 
  ! Evaluates source/sink term storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Grid_module
  use Material_Aux_class, only: material_auxvar_type
  use SrcSink_Sandbox_module
  use SrcSink_Sandbox_Base_class
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"

  PetscBool :: compute_derivative
  Vec :: residual
  Mat :: Jacobian
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(general_auxvar_type), pointer :: general_auxvars(:,:)
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscReal, pointer :: r_p(:)
  PetscReal :: res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  PetscInt :: local_id, ghosted_id, istart, iend, irow, idof
  PetscReal :: res_pert(option%nflowdof)
  PetscReal :: aux_real(10)
  PetscErrorCode :: ierr
  
  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif
  
  cur_srcsink => ss_sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    aux_real = 0.d0
    local_id = cur_srcsink%local_cell_id
    ghosted_id = grid%nL2G(local_id)
    res = 0.d0
    Jac = 0.d0
    call GeneralSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                      general_auxvars(ZERO_INTEGER,ghosted_id),option)
    call cur_srcsink%Evaluate(res,Jac,PETSC_FALSE, &
                              material_auxvars(ghosted_id), &
                              aux_real,option)
    if (compute_derivative) then
      do idof = 1, option%nflowdof
        res_pert = 0.d0
        call GeneralSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                                    general_auxvars(idof,ghosted_id),option)
        call cur_srcsink%Evaluate(res_pert,Jac,PETSC_FALSE, &
                                  material_auxvars(ghosted_id), &
                                  aux_real,option)
        do irow = 1, option%nflowdof
          Jac(irow,idof) = (res_pert(irow)-res(irow)) / &
                            general_auxvars(idof,ghosted_id)%pert
        enddo
      enddo
      if (general_isothermal) then
        Jac(GENERAL_ENERGY_EQUATION_INDEX,:) = 0.d0
        Jac(:,GENERAL_ENERGY_EQUATION_INDEX) = 0.d0
      endif         
      if (general_no_air) then
        Jac(GENERAL_GAS_EQUATION_INDEX,:) = 0.d0
        Jac(:,GENERAL_GAS_EQUATION_INDEX) = 0.d0
      endif          
      call MatSetValuesBlockedLocal(Jacobian,1,ghosted_id-1,1, &
                                    ghosted_id-1,Jac,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
    else
      iend = local_id*option%nflowdof
      istart = iend - option%nflowdof + 1
      r_p(istart:iend) = r_p(istart:iend) - res
    endif
    cur_srcsink => cur_srcsink%next
  enddo
  
  if (.not.compute_derivative) then
    call VecRestoreArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine GeneralSSSandbox

! ************************************************************************** !

subroutine GeneralSSSandboxLoadAuxReal(srcsink,aux_real,gen_auxvar,option)

  use Option_module
  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_WIPP_Gas_class
  use SrcSink_Sandbox_WIPP_Well_class

  implicit none

  class(srcsink_sandbox_base_type) :: srcsink
  PetscReal :: aux_real(:)
  type(general_auxvar_type) gen_auxvar
  type(option_type) :: option
  
  aux_real = 0.d0
  select type(srcsink)
    class is(srcsink_sandbox_wipp_gas_type)
      aux_real(WIPP_GAS_WATER_SATURATION_INDEX) = &
        gen_auxvar%sat(option%liquid_phase)
      aux_real(WIPP_GAS_TEMPERATURE_INDEX) = &
        gen_auxvar%temp
    class is(srcsink_sandbox_wipp_well_type)
      aux_real(WIPP_WELL_LIQUID_MOBILITY) = &
        gen_auxvar%mobility(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_MOBILITY) = &
        gen_auxvar%mobility(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_PRESSURE) = &
        gen_auxvar%pres(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_PRESSURE) = &
        gen_auxvar%pres(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_ENTHALPY) = &
        gen_auxvar%H(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_ENTHALPY) = &
        gen_auxvar%H(option%gas_phase)
      aux_real(WIPP_WELL_XMOL_AIR_IN_LIQUID) = &
        gen_auxvar%xmol(option%air_id,option%liquid_phase)
      aux_real(WIPP_WELL_XMOL_WATER_IN_GAS) = &
        gen_auxvar%xmol(option%water_id,option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_DENSITY) = &
        gen_auxvar%den(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_DENSITY) = &
        gen_auxvar%den(option%gas_phase)
  end select
  
end subroutine GeneralSSSandboxLoadAuxReal

! ************************************************************************** !

subroutine GeneralMapBCAuxVarsToGlobal(realization)
  ! 
  ! Maps variables in general auxvar to global equivalent.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
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
  type(general_auxvar_type), pointer :: gen_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update
  
  gen_auxvars_bc => patch%aux%General%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        gen_auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        gen_auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        gen_auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine GeneralMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine GeneralDestroy(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/11
  ! 

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine GeneralDestroy

end module General_module
