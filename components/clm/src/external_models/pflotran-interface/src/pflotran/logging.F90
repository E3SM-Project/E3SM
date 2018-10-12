module Logging_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! stages
  PetscInt, parameter, public :: INIT_STAGE = 1
  PetscInt, parameter, public :: TS_STAGE = 2
  PetscInt, parameter, public :: OUTPUT_STAGE = 3
  PetscInt, parameter, public :: FINAL_STAGE = 4

  type, public :: logging_type 
  
    PetscInt :: stage_count
    PetscLogStage :: stage(10)
    
    PetscClassId :: class_pflotran
    
    PetscLogEvent :: event_init
    PetscLogEvent :: event_setup
   
    PetscLogEvent :: event_create_iogroups

    PetscLogEvent :: event_restart
    PetscLogEvent :: event_checkpoint

    PetscLogEvent :: event_flow_condition_read
    PetscLogEvent :: event_tran_condition_read
    PetscLogEvent :: event_tran_constraint_read
    PetscLogEvent :: event_flow_condition_read_values

    PetscLogEvent :: event_h5dread_f
    PetscLogEvent :: event_h5dwrite_f
    PetscLogEvent :: event_read_indices_hdf5
    PetscLogEvent :: event_map_indices_hdf5
    PetscLogEvent :: event_hash_create
    PetscLogEvent :: event_hash_map
    PetscLogEvent :: event_read_real_array_hdf5
    PetscLogEvent :: event_read_ndim_real_array_hdf5
    PetscLogEvent :: event_read_int_array_hdf5
    PetscLogEvent :: event_write_real_array_hdf5
    PetscLogEvent :: event_write_int_array_hdf5
    PetscLogEvent :: event_read_array_hdf5    
    PetscLogEvent :: event_read_xyz_dataset_hdf5    
    PetscLogEvent :: event_write_struct_dataset_hdf5
    PetscLogEvent :: event_region_read_hdf5
    PetscLogEvent :: event_region_read_ascii
    PetscLogEvent :: event_cell_indx_int_read_hdf5
    PetscLogEvent :: event_cell_indx_real_read_hdf5
    PetscLogEvent :: event_dataset_gridded_hdf5_read
    PetscLogEvent :: event_dataset_map_hdf5_read

    PetscLogEvent :: event_output_tecplot
    PetscLogEvent :: event_output_hdf5
    PetscLogEvent :: event_output_vtk
    PetscLogEvent :: event_output_grid_vtk
    PetscLogEvent :: event_output_write_vtk
    PetscLogEvent :: event_output_mad
    PetscLogEvent :: event_output_str_grid_tecplot
    PetscLogEvent :: event_output_write_tecplot
    PetscLogEvent :: event_output_write_flux_tecplot
    PetscLogEvent :: event_output_get_var_from_array
    PetscLogEvent :: event_output_get_cell_vel
    PetscLogEvent :: event_output_vec_tecplot
    PetscLogEvent :: event_output_observation
    PetscLogEvent :: event_output_coordinates_hdf5
    PetscLogEvent :: event_output_hydrograph
    PetscLogEvent :: event_output_secondary_tecplot

    PetscLogEvent :: event_r_residual
    PetscLogEvent :: event_r_jacobian
    PetscLogEvent :: event_r_auxvars
    PetscLogEvent :: event_r_auxvars_bc
    
    PetscLogEvent :: event_rt_residual
    PetscLogEvent :: event_rt_jacobian

    PetscLogEvent :: event_rt_jacobian_flux
    PetscLogEvent :: event_rt_jacobian_fluxbc
    PetscLogEvent :: event_rt_jacobian_accum
    PetscLogEvent :: event_rt_jacobian_zero_calc
    PetscLogEvent :: event_rt_jacobian_zero
    PetscLogEvent :: event_rt_jacobian_ss
    PetscLogEvent :: event_rt_jacobian1
    PetscLogEvent :: event_rt_jacobian2

    PetscLogEvent :: event_rt_res_reaction
    PetscLogEvent :: event_rt_jac_reaction
    PetscLogEvent :: event_rt_react
    PetscLogEvent :: event_rt_auxvars
    PetscLogEvent :: event_rt_auxvars_bc
    
    PetscLogEvent :: event_mass_balance

    PetscBool :: allow_new_stages

  end type logging_type
  
  type(logging_type), pointer, public :: logging
  
  public :: LoggingCreate, &
            LoggingCreateStage, &
            LoggingSetupComplete, &
            LoggingDestroy

contains

! ************************************************************************** !

subroutine LoggingCreate()
  ! 
  ! Allocates and initializes a new logging object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  PetscErrorCode :: ierr
  
  allocate(logging)

  logging%allow_new_stages = PETSC_TRUE
  logging%stage_count = FINAL_STAGE
  
  call PetscLogStageRegister('Init Stage',  & 
                             logging%stage(INIT_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogStageRegister('Time Step Stage', &
                             logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogStageRegister('Output Stage', &
                             logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogStageRegister('Finalization Stage', &
                             logging%stage(FINAL_STAGE),ierr);CHKERRQ(ierr)
                             
!!  call PetscCookieRegister('PFLOTRAN',logging%class_pflotran,ierr)
  call PetscClassIdRegister('PFLOTRAN',logging%class_pflotran, &
                            ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('Init', &
                             logging%class_pflotran, &
                             logging%event_init,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('Init,Setup', &
                             logging%class_pflotran, &
                             logging%event_setup,ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('Create_iogroups', &
                             logging%class_pflotran, &
                             logging%event_create_iogroups,ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('Restart', &
                             logging%class_pflotran, &
                             logging%event_restart,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('Checkpoint', &
                             logging%class_pflotran, &
                             logging%event_checkpoint,ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('FlowCondRead', &
                             logging%class_pflotran, &
                             logging%event_flow_condition_read, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('TranCondRead', &
                             logging%class_pflotran, &
                             logging%event_tran_condition_read, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('TranConstraintRd', &
                             logging%class_pflotran, &
                             logging%event_tran_constraint_read, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('FlowCondReadVals', &
                             logging%class_pflotran, &
                             logging%event_flow_condition_read_values, &
                             ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('H5DRead_F', &
                             logging%class_pflotran, &
                             logging%event_h5dread_f,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5DWrite_F', &
                             logging%class_pflotran, &
                             logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('HDF5ReadIndices', &
                             logging%class_pflotran, &
                             logging%event_read_indices_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5MapLoc2NatIndx', &
                             logging%class_pflotran, &
                             logging%event_map_indices_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('GrdCrNat2GhstHsh', &
                             logging%class_pflotran, &
                             logging%event_hash_create,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('GrdLocGhstIdHsh', &
                             logging%class_pflotran, &
                             logging%event_hash_map,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5ReadRealArray', &
                             logging%class_pflotran, &
                             logging%event_read_real_array_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5ReadNDimRealArray', &
                             logging%class_pflotran, &
                             logging%event_read_ndim_real_array_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5ReadIntArray', &
                             logging%class_pflotran, &
                             logging%event_read_int_array_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5ReadArray', &
                             logging%class_pflotran, &
                             logging%event_read_array_hdf5,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5WriteRealArray', &
                             logging%class_pflotran, &
                             logging%event_write_real_array_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5WriteIntArray', &
                             logging%class_pflotran, &
                             logging%event_write_int_array_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5WriteStrData', &
                             logging%class_pflotran, &
                             logging%event_write_struct_dataset_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5ReadRegFrmFile', &
                             logging%class_pflotran, &
                             logging%event_region_read_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RegReadFrmFileId', &
                             logging%class_pflotran, &
                             logging%event_region_read_ascii, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5RdCellIndxInt', &
                             logging%class_pflotran, &
                             logging%event_cell_indx_int_read_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('H5RdCellIndxReal', &
                             logging%class_pflotran, &
                             logging%event_cell_indx_real_read_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('DatasetGriddedHDF5Read', &
                             logging%class_pflotran, &
                             logging%event_dataset_gridded_hdf5_read, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('DatasetMapHDF5Read', &
                             logging%class_pflotran, &
                             logging%event_dataset_map_hdf5_read, &
                             ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('OutputTecplot', &
                             logging%class_pflotran, &
                             logging%event_output_tecplot,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputVTK', &
                             logging%class_pflotran, &
                             logging%event_output_vtk,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputMAD', &
                             logging%class_pflotran, &
                             logging%event_output_mad,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputHDF5', &
                             logging%class_pflotran, &
                             logging%event_output_hdf5,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputSecondaryTecplot', &
                             logging%class_pflotran, &
                             logging%event_output_secondary_tecplot, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('WriteTecStrGrid', &
                             logging%class_pflotran, &
                             logging%event_output_str_grid_tecplot, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('WriteVTKGrid', &
                             logging%class_pflotran, &
                             logging%event_output_grid_vtk,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('WriteTecDataSet', &
                             logging%class_pflotran, &
                             logging%event_output_write_tecplot, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('WriteVTKDataSet', &
                             logging%class_pflotran, &
                             logging%event_output_write_vtk, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputFluxVelTec', &
                             logging%class_pflotran, &
                             logging%event_output_write_flux_tecplot, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputGtVrFrmArr', &
                             logging%class_pflotran, &
                             logging%event_output_get_var_from_array, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('GetCellCentVel', &
                             logging%class_pflotran, &
                             logging%event_output_get_cell_vel, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputFluxVelTec', &
                             logging%class_pflotran, &
                             logging%event_output_vec_tecplot, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputBrkthuTec', &
                             logging%class_pflotran, &
                             logging%event_output_observation, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('WriteHDF5Coord', &
                             logging%class_pflotran, &
                             logging%event_output_coordinates_hdf5, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('OutputHydrograph', &
                             logging%class_pflotran, &
                             logging%event_output_hydrograph, &
                             ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('RResidual', &
                             logging%class_pflotran, &
                             logging%event_r_residual,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RJacobian', &
                             logging%class_pflotran, &
                             logging%event_r_jacobian,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RAuxVars', &
                             logging%class_pflotran, &
                             logging%event_r_auxvars,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RAuxVarsBC', &
                             logging%class_pflotran, &
                             logging%event_r_auxvars_bc,ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('RTResidual', &
                             logging%class_pflotran, &
                             logging%event_rt_residual,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacobian', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian,ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('RTJacobianFlux', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_flux, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacobianFluxBC', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_fluxbc, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacobianAccum', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_accum, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacobianZeroCalc', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_zero_calc, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacobianZero', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_zero, &
                             ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacobianSS', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian_ss,ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('RTJacobian1', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian1,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacobian2', &
                             logging%class_pflotran, &
                             logging%event_rt_jacobian2,ierr);CHKERRQ(ierr)

  call PetscLogEventRegister('RTResReaction', &
                             logging%class_pflotran, &
                             logging%event_rt_res_reaction,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTJacReaction', &
                             logging%class_pflotran, &
                             logging%event_rt_jac_reaction,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTReact', &
                             logging%class_pflotran, &
                             logging%event_rt_react,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTAuxVars', &
                             logging%class_pflotran, &
                             logging%event_rt_auxvars,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('RTAuxVarsBC', &
                             logging%class_pflotran, &
                             logging%event_rt_auxvars_bc,ierr);CHKERRQ(ierr)
                             
  call PetscLogEventRegister('MassBalance', &
                             logging%class_pflotran, &
                             logging%event_mass_balance,ierr);CHKERRQ(ierr)

end subroutine LoggingCreate

! ************************************************************************** !

subroutine LoggingCreateStage(stage_name,stage_id)
  ! 
  ! Creates a new custom stage
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: stage_name
  PetscInt :: stage_id
  
  character(len=MAXSTRINGLENGTH) :: full_stage_name
  character(len=MAXSTRINGLENGTH) :: temp_stage_name
  character(len=MAXWORDLENGTH) :: word
  PetscLogStage :: temp_stage_id
  PetscInt :: i
  PetscErrorCode :: ierr

  ! this conditional prevents duplicate stages that can be generated during
  ! multirealization simulations.
  if (.not. logging%allow_new_stages) return
  
  logging%stage_count = logging%stage_count + 1
  full_stage_name = trim(stage_name) // ' Stage'
  !TODO(geh): fix after PETSc fixes bug in implementation.  PetscLogStageGetId
  !           currently returns the number of stages, not -1, when one exists.
  ! No two stages can have the same name
#if 0
  i = 0
  temp_stage_name = full_stage_name
  do
    ! check if stage exists
    call PetscLogStageGetId(temp_stage_name,temp_stage_id,ierr);CHKERRQ(ierr)
    if (temp_stage_id > -1) then
      i = i + 1
      write(word,*) i
      ! append count
      temp_stage_name = trim(full_stage_name) // trim(adjustl(word))
    else
      full_stage_name = temp_stage_name
      exit
    endif
  enddo
#endif
  call PetscLogStageRegister(full_stage_name,stage_id,ierr);CHKERRQ(ierr)
  
  stage_id = logging%stage_count
  
end subroutine LoggingCreateStage
  
! ************************************************************************** !

subroutine LoggingSetupComplete()
  ! 
  ! Sets flag that indicates that setup is complete.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/30/14
  ! 

  implicit none

  logging%allow_new_stages = PETSC_FALSE
  
end subroutine LoggingSetupComplete

! ************************************************************************** !

subroutine LoggingDestroy()
  ! 
  ! Deallocates a logging object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  implicit none
  
  ! all kinds of stuff needs to be added here.
  
  deallocate(logging)
  nullify(logging)
  
end subroutine LoggingDestroy

end module Logging_module
