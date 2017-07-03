module Geomechanics_Force_module

  use Geomechanics_Global_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
!#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscts.h"

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: GeomechForceSetup, &
            GeomechForceUpdateAuxVars, &
            GeomechanicsForceInitialGuess, &
            GeomechForceResidual, &
            GeomechForceJacobian, &
            GeomechUpdateFromSubsurf, &
            GeomechUpdateSubsurfFromGeomech, &
            GeomechCreateGeomechSubsurfVec, &
            GeomechCreateSubsurfStressStrainVec, &
            GeomechUpdateSolution, &
            GeomechStoreInitialPressTemp, &
            GeomechStoreInitialDisp, &
            GeomechStoreInitialPorosity, &
            GeomechUpdateSubsurfPorosity, &
            GeomechForceJacobianLinearPart
 
contains

! ************************************************************************** !

subroutine GeomechForceSetup(geomech_realization)
  ! 
  ! Sets up the geomechanics calculations
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Output_Aux_module

  class(realization_geomech_type) :: geomech_realization
  type(output_variable_list_type), pointer :: list

  call GeomechForceSetupPatch(geomech_realization)

  list => geomech_realization%output_option%output_snap_variable_list
  call GeomechForceSetPlotVariables(list)
  list => geomech_realization%output_option%output_obs_variable_list
  call GeomechForceSetPlotVariables(list)
   
end subroutine GeomechForceSetup

! ************************************************************************** !

subroutine GeomechForceSetupPatch(geomech_realization)
  ! 
  ! Sets up the arrays for geomech parameters
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/11/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Option_module
 
  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch

  PetscInt :: i

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch

  allocate(patch%geomech_aux%GeomechParam%youngs_modulus &
    (size(geomech_realization%geomech_material_property_array)))
  allocate(patch%geomech_aux%GeomechParam%poissons_ratio &
    (size(geomech_realization%geomech_material_property_array)))
  allocate(patch%geomech_aux%GeomechParam%biot_coef &
    (size(geomech_realization%geomech_material_property_array)))
  allocate(patch%geomech_aux%GeomechParam%thermal_exp_coef &
    (size(geomech_realization%geomech_material_property_array)))
  allocate(patch%geomech_aux%GeomechParam%density &
    (size(geomech_realization%geomech_material_property_array)))

  do i = 1, size(geomech_realization%geomech_material_property_array)
    patch%geomech_aux%GeomechParam%youngs_modulus(geomech_realization% &
      geomech_material_property_array(i)%ptr%id) = geomech_realization% &
      geomech_material_property_array(i)%ptr%youngs_modulus
    patch%geomech_aux%GeomechParam%poissons_ratio(geomech_realization% &
      geomech_material_property_array(i)%ptr%id) = geomech_realization% &
      geomech_material_property_array(i)%ptr%poissons_ratio
    patch%geomech_aux%GeomechParam%density(geomech_realization% &
      geomech_material_property_array(i)%ptr%id) = geomech_realization% &
      geomech_material_property_array(i)%ptr%density
    patch%geomech_aux%GeomechParam%biot_coef(geomech_realization% &
      geomech_material_property_array(i)%ptr%id) = geomech_realization% &
      geomech_material_property_array(i)%ptr%biot_coeff
    patch%geomech_aux%GeomechParam%thermal_exp_coef(geomech_realization% &
      geomech_material_property_array(i)%ptr%id) = geomech_realization% &
      geomech_material_property_array(i)%ptr%thermal_exp_coeff
  enddo

end subroutine GeomechForceSetupPatch

! ************************************************************************** !

subroutine GeomechForceSetPlotVariables(list)
  ! 
  ! Set up of geomechanics plot variables
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 
  
  use Output_Aux_module
  use Variables_module
    
  implicit none

  type(output_variable_list_type), pointer :: list
  type(output_variable_type), pointer :: output_variable

  character(len=MAXWORDLENGTH) :: name, units
  
  if (associated(list%first)) then
    return
  endif

  name = 'displacement_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_X)
                               
  name = 'displacement_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_Y)
                               
  name = 'displacement_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_Z)

  units = ''
  name = 'Material ID'
  output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                          units,GEOMECH_MATERIAL_ID)
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList(list,output_variable)
                             
  name = 'strain_xx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_XX)
                               
  name = 'strain_yy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_YY)
                               
  name = 'strain_zz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_ZZ)
                               
  name = 'strain_xy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_XY)
                               
  name = 'strain_yz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_YZ)
                               
  name = 'strain_zx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_ZX)
                                                                              
  name = 'stress_xx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_XX)
                               
  name = 'stress_yy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_YY)
                               
  name = 'stress_zz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_ZZ)
                               
  name = 'stress_xy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_XY)
                               
  name = 'stress_yz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_YZ)
                               
  name = 'stress_zx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_ZX)
  name = 'relative_displacement_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_X)
                               
  name = 'relative_displacement_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_Y)
                               
  name = 'relative_displacement_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_Z)


end subroutine GeomechForceSetPlotVariables

! ************************************************************************** !

subroutine GeomechanicsForceInitialGuess(geomech_realization)
  ! 
  ! Sets up the inital guess for the solution
  ! The boundary conditions are set here
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/19/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Option_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Patch_module
  use Geomechanics_Coupler_module
  use Geomechanics_Region_module
  
  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: field
  type(geomech_patch_type), pointer :: patch
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_grid_type), pointer :: grid
  type(gm_region_type), pointer :: region
  
  PetscInt :: ghosted_id,local_id,total_verts,ivertex
  PetscReal, pointer :: xx_p(:)
  PetscErrorCode :: ierr
  
  option => geomech_realization%option
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  
  call GeomechGridVecGetArrayF90(grid,field%disp_xx,xx_p,ierr)
  
  boundary_condition => patch%geomech_boundary_condition_list%first
  total_verts = 0
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      total_verts = total_verts + 1
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif    
      
      ! X displacement 
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_X_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Y displacement 
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Y_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Z displacement      
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Z_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
    enddo
    boundary_condition => boundary_condition%next      
  enddo
  
  call GeomechGridVecRestoreArrayF90(grid,field%disp_xx,xx_p,ierr)

end subroutine GeomechanicsForceInitialGuess

! ************************************************************************** !

subroutine GeomechForceUpdateAuxVars(geomech_realization)
  ! 
  ! Updates the geomechanics variables
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/18/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Option_module
  use Geomechanics_Field_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Material_module
  use Geomechanics_Global_Aux_module
  use Geomechanics_Region_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
  type(geomech_grid_type), pointer :: grid
  type(geomech_field_type), pointer :: geomech_field
  type(gm_region_type), pointer :: region
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)

  PetscInt :: ghosted_id, local_id
  PetscReal, pointer :: xx_loc_p(:), xx_init_loc_p(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  geomech_field => geomech_realization%geomech_field

  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars
  
  call GeomechGridVecGetArrayF90(grid,geomech_field%disp_xx_loc,xx_loc_p,ierr)
  call GeomechGridVecGetArrayF90(grid,geomech_field%disp_xx_init_loc, &
                                 xx_init_loc_p,ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax_node
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_X_DOF) = &
      xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)
 
    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_X_DOF) = &
      xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)
 enddo
   
  call GeomechGridVecRestoreArrayF90(grid,geomech_field%disp_xx_loc, &
                                     xx_loc_p,ierr)
  call GeomechGridVecRestoreArrayF90(grid,geomech_field%disp_xx_init_loc, &
                                     xx_init_loc_p,ierr)


end subroutine GeomechForceUpdateAuxVars

! ************************************************************************** !

subroutine GeomechForceResidual(snes,xx,r,geomech_realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Satish Karra
  ! Date: 06/21/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(realization_geomech_type) :: geomech_realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: field
  type(option_type), pointer :: option
  
  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  option => geomech_realization%option

  ! Communication -----------------------------------------
  call GeomechDiscretizationGlobalToLocal(geomech_discretization,xx, &
                                          field%disp_xx_loc,NGEODOF)
  
  call GeomechForceResidualPatch(snes,xx,r,geomech_realization,ierr)

  if (geomech_realization%geomech_debug%vecview_residual) then
    call PetscViewerASCIIOpen(geomech_realization%option%mycomm, &
                              'Geomech_residual.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  endif
  
  if (geomech_realization%geomech_debug%vecview_solution) then
    call PetscViewerASCIIOpen(geomech_realization%option%mycomm, &
                              'Geomech_xx.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine GeomechForceResidual

! ************************************************************************** !

subroutine GeomechForceResidualPatch(snes,xx,r,geomech_realization,ierr)
  ! 
  ! Computes the residual equation on a patch
  ! 
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module
  use Geomechanics_Auxiliary_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(realization_geomech_type) :: geomech_realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_parameter_type), pointer :: GeomechParam

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, allocatable :: local_press(:), local_temp(:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: res_vec(:)
  PetscReal, pointer :: press(:), temp(:)
  PetscReal, pointer :: press_init(:), temp_init(:)
  PetscReal, allocatable :: beta_vec(:), alpha_vec(:)
  PetscReal, allocatable :: density_vec(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscInt :: ielem, ivertex 
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: petsc_id, local_id
  PetscReal :: error_H1_global, error_L2_global
  PetscReal :: error_L2, error_H1
  PetscReal, pointer :: imech_loc_p(:)      
  PetscInt :: size_elenodes
                    
  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars  
  GeomechParam => patch%geomech_aux%GeomechParam 

  call GeomechForceUpdateAuxVars(geomech_realization)
  ! Add flag for the update
  
  call VecSet(r,0.d0,ierr);CHKERRQ(ierr)
  
#if 0  
  error_H1_global = 0.d0
  error_L2_global = 0.d0
#endif

  ! Get pressure and temperature from subsurface
  call VecGetArrayF90(field%press_loc,press,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_loc,temp,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  ! Get initial pressure and temperature 
  call VecGetArrayF90(field%press_init_loc,press_init,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%temp_init_loc,temp_init,ierr);CHKERRQ(ierr)
 
  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes),option%ngeomechdof))
    allocate(local_press(size(elenodes)))
    allocate(local_temp(size(elenodes)))
    allocate(petsc_ids(size(elenodes)))
    allocate(ids(size(elenodes)*option%ngeomechdof))
    allocate(res_vec(size(elenodes)*option%ngeomechdof))
    allocate(beta_vec(size(elenodes)))
    allocate(alpha_vec(size(elenodes)))
    allocate(density_vec(size(elenodes)))
    allocate(youngs_vec(size(elenodes)))
    allocate(poissons_vec(size(elenodes)))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%EleType
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(ivertex,idof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
      local_press(ivertex) = press(ghosted_id) - press_init(ghosted_id)  ! p - p_0
      local_temp(ivertex) = temp(ghosted_id) - temp_init(ghosted_id)     ! T - T_0
      alpha_vec(ivertex) = &
        GeomechParam%thermal_exp_coef(int(imech_loc_p(ghosted_id))) 
      beta_vec(ivertex) = &
        GeomechParam%biot_coef(int(imech_loc_p(ghosted_id))) 
      density_vec(ivertex) = &
        GeomechParam%density(int(imech_loc_p(ghosted_id))) 
      youngs_vec(ivertex) = &
        GeomechParam%youngs_modulus(int(imech_loc_p(ghosted_id))) 
      poissons_vec(ivertex) = &
        GeomechParam%poissons_ratio(int(imech_loc_p(ghosted_id))) 
    enddo
    size_elenodes = size(elenodes)
    call GeomechForceLocalElemResidual(size_elenodes,local_coordinates, &
       local_disp,local_press,local_temp,youngs_vec,poissons_vec, &
       density_vec,beta_vec,alpha_vec,eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,res_vec,option)
    call VecSetValues(r,size(ids),ids,res_vec,ADD_VALUES,ierr);CHKERRQ(ierr)
#if 0
    call GeomechForceLocalElemError(size_elenodes,local_coordinates, &
                                    local_disp, &
                                    eletype,grid%gauss_node(ielem)%dim, &
                                    grid%gauss_node(ielem)%r, &
                                    grid%gauss_node(ielem)%w,error_L2, &
                                    error_H1,option)
    error_H1_global = error_H1_global + error_H1
    error_L2_global = error_L2_global + error_L2
#endif
    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(petsc_ids)
    deallocate(ids)
    deallocate(res_vec)
    deallocate(local_press)
    deallocate(local_temp)
    deallocate(beta_vec)
    deallocate(alpha_vec)
    deallocate(density_vec)
    deallocate(youngs_vec)
    deallocate(poissons_vec)
  enddo
      
  call VecRestoreArrayF90(field%press_loc,press,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%temp_loc,temp,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
    
  call VecRestoreArrayF90(field%press_init_loc,press_init,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%temp_init_loc,temp_init,ierr);CHKERRQ(ierr)

#if 0
  call MPI_Allreduce(error_H1_global,error_H1_global,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr)      
  call MPI_Allreduce(error_L2_global,error_L2_global,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr)   
                     
  if (option%myrank == option%io_rank) then                   
    print *, 'L2 error:', sqrt(error_L2_global)
    print *, 'H1 error:', sqrt(error_H1_global)
  endif
#endif      
      
      
  call VecAssemblyBegin(r,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(r,ierr);CHKERRQ(ierr)

  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  ! displacement boundary conditions
  boundary_condition => patch%geomech_boundary_condition_list%first
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif    
      
      ! X displacement 
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_X_DOF-1,0.d0,INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call printErrMsg(option)
        end select
      endif
      
      ! Y displacement 
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Y_DOF-1,0.d0,INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call printErrMsg(option)
        end select
      endif
      
      ! Z displacement      
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Z_DOF-1,0.d0,INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for displacement not available.'
            call printErrMsg(option)
        end select
      endif
      
    enddo
    boundary_condition => boundary_condition%next      
  enddo

  ! Need to assemby here since one cannot mix INSERT_VALUES
  ! and ADD_VALUES
  call VecAssemblyBegin(r,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(r,ierr);CHKERRQ(ierr)
  
  ! Force boundary conditions
  boundary_condition => patch%geomech_boundary_condition_list%first
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif    
       
      ! X force 
      if (associated(boundary_condition%geomech_condition%force_x)) then
        select case(boundary_condition%geomech_condition%force_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_X_DOF-1, &
              -boundary_condition%geomech_aux_real_var &
              (GEOMECH_DISP_X_DOF,ivertex),ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call printErrMsg(option)
        end select
      endif
      
       ! Y force 
      if (associated(boundary_condition%geomech_condition%force_y)) then
        select case(boundary_condition%geomech_condition%force_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Y_DOF-1, &
              -boundary_condition%geomech_aux_real_var &
              (GEOMECH_DISP_Y_DOF,ivertex),ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call printErrMsg(option)

        end select
      endif

       ! Z force 
      if (associated(boundary_condition%geomech_condition%force_z)) then
        select case(boundary_condition%geomech_condition%force_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(r,(petsc_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Z_DOF-1, &
              -boundary_condition%geomech_aux_real_var &
              (GEOMECH_DISP_Z_DOF,ivertex),ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            option%io_buffer = 'Neumann BC for force not available.'
            call printErrMsg(option)
        end select
      endif
 
    enddo
    boundary_condition => boundary_condition%next      
  enddo

  call VecAssemblyBegin(r,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(r,ierr);CHKERRQ(ierr)

end subroutine GeomechForceResidualPatch

! ************************************************************************** !

subroutine GeomechForceLocalElemResidual(size_elenodes,local_coordinates, &
                                         local_disp, &
                                         local_press,local_temp, &
                                         local_youngs,local_poissons, &
                                         local_density,local_beta, &
                                         local_alpha, &
                                         eletype,dim,r,w,res_vec,option)
  ! 
  ! Computes the residual for a local element
  ! 
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! 
                                         
  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module
  
  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: res_vec(:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, allocatable :: local_press(:)
  PetscReal, allocatable :: local_temp(:)
  PetscReal, allocatable :: local_youngs(:)
  PetscReal, allocatable :: local_poissons(:)
  PetscReal, allocatable :: local_density(:)
  PetscReal, allocatable :: local_beta(:)
  PetscReal, allocatable :: local_alpha(:)
      
  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j,d
  PetscReal :: eye_three(THREE_INTEGER)
  PetscInt :: indx(THREE_INTEGER)
  PetscInt :: dim
  PetscReal :: lambda, mu, beta, alpha
  PetscReal :: density, youngs_mod, poissons_ratio
  PetscInt :: load_type
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal, allocatable :: N(:,:)
  PetscReal, allocatable :: vecB_transpose(:,:)
  PetscReal, allocatable :: kron_B_eye(:,:)
  PetscReal, allocatable :: kron_B_transpose_eye(:,:)
  PetscReal, allocatable :: Trans(:,:)
  PetscReal, allocatable :: kron_eye_B_transpose(:,:)
  PetscReal, allocatable :: kron_N_eye(:,:)
  PetscReal, allocatable :: vec_local_disp(:,:)
  PetscReal, allocatable :: force(:), res_vec_mat(:,:)
  PetscInt :: size_elenodes
  
  allocate(B(size_elenodes,dim))
  allocate(Kmat(size_elenodes*option%ngeomechdof, &
                size_elenodes*option%ngeomechdof))  
  allocate(force(size_elenodes*option%ngeomechdof))
  allocate(res_vec_mat(size_elenodes*option%ngeomechdof,1))
  
  res_vec = 0.d0
  res_vec_mat = 0.d0
  Kmat = 0.d0
  force = 0.d0
  len_w = size(w)
  
  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo
  
  call Transposer(option%ngeomechdof,size_elenodes,Trans)
 
  do igpt = 1, len_w
    shapefunction%EleType = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    allocate(N(size(shapefunction%N),ONE_INTEGER))
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!' 
      call printErrMsg(option)        
    endif
    ! Find the inverse of J_map
    ! Set identity matrix
    call ludcmp(J_map,THREE_INTEGER,indx,d)
    do i = 1, THREE_INTEGER
      eye_three = 0.d0
      eye_three(i) = 1.d0
      call lubksb(J_map,THREE_INTEGER,indx,eye_three)
      inv_J_map(:,i) = eye_three
    enddo
    B = matmul(shapefunction%DN,inv_J_map)
    youngs_mod = dot_product(shapefunction%N,local_youngs)
    poissons_ratio = dot_product(shapefunction%N,local_poissons)
    alpha = dot_product(shapefunction%N,local_alpha)
    beta = dot_product(shapefunction%N,local_beta)
    density = dot_product(shapefunction%N,local_density) 
    call GeomechGetLambdaMu(lambda,mu,youngs_mod,poissons_ratio)
    call GeomechGetBodyForce(load_type,lambda,mu,x,bf,option) 
    call ConvertMatrixToVector(transpose(B),vecB_transpose)
    Kmat = Kmat + w(igpt)*lambda* &
      matmul(vecB_transpose,transpose(vecB_transpose))*detJ_map
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    N(:,1)= shapefunction%N    
    call Kron(N,identity,kron_N_eye)
    Kmat = Kmat + w(igpt)*mu*matmul(kron_B_eye,kron_B_transpose_eye)*detJ_map
    Kmat = Kmat + w(igpt)*mu* &
      matmul(matmul(kron_B_eye,kron_eye_B_transpose),Trans)*detJ_map
    force = force + w(igpt)*density*matmul(kron_N_eye,bf)*detJ_map
    force = force + w(igpt)*beta*dot_product(N(:,1),local_press)* &
      vecB_transpose(:,1)*detJ_map
    force = force + w(igpt)*alpha*(3*lambda+2*mu)* &
      dot_product(N(:,1),local_temp)*vecB_transpose(:,1)*detJ_map  
    call ShapeFunctionDestroy(shapefunction)
    deallocate(N)
    deallocate(vecB_transpose)
    deallocate(kron_B_eye)
    deallocate(kron_B_transpose_eye)
    deallocate(kron_eye_B_transpose)
    deallocate(kron_N_eye)
  enddo
  
  call ConvertMatrixToVector(transpose(local_disp),vec_local_disp)
  res_vec_mat = matmul(Kmat,vec_local_disp)
  res_vec = res_vec + res_vec_mat(:,1)
  res_vec = res_vec - force

  deallocate(B)
  deallocate(force)
  deallocate(Kmat)
  deallocate(res_vec_mat)
  deallocate(vec_local_disp)
  deallocate(Trans)

end subroutine GeomechForceLocalElemResidual

! ************************************************************************** !
!
! GeomechForceLocalElemError: Computes the error for a local element
! author: Satish Karra
! date: 07/08/13
!
! ************************************************************************** !
#if 0

! ************************************************************************** !

subroutine GeomechForceLocalElemError(size_elenodes,local_coordinates, &
                                      local_disp, &
                                      eletype,dim,r,w,error_L2,error_H1,option)
                                         
  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module
  
  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: res_vec(:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: u(THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j,d
  PetscReal :: eye_three(THREE_INTEGER)
  PetscInt :: indx(THREE_INTEGER)
  PetscInt :: dim
  PetscReal :: lambda, mu
  PetscInt :: load_type
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: den_rock
  PetscReal :: grad_u(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: grad_u_exact(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: u_exact(THREE_INTEGER)
  PetscReal :: error_H1, error_L2
  PetscReal :: trace_disp, trace_disp_grad
  
  allocate(B(size_elenodes,dim))
  
  error_H1 = 0.d0
  error_L2 = 0.d0

  len_w = size(w)
  
  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo
    
  do igpt = 1, len_w
    shapefunction%EleType = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    u = matmul(transpose(local_disp),shapefunction%N)
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!' 
      call printErrMsg(option)        
    endif
    ! Find the inverse of J_map
    ! Set identity matrix
    call ludcmp(J_map,THREE_INTEGER,indx,d)
    do i = 1, THREE_INTEGER
      eye_three = 0.d0
      eye_three(i) = 1.d0
      call lubksb(J_map,THREE_INTEGER,indx,eye_three)
      inv_J_map(:,i) = eye_three
    enddo
    B = matmul(shapefunction%DN,inv_J_map)
    grad_u = matmul(transpose(local_disp),B)
    call GeomechGetLambdaMu(lambda,mu,x)
    load_type = 2 ! Need to change
    call GeomechGetBodyForce(load_type,lambda,mu,x,bf,option) 
    call GetAnalytical(load_type,lambda,mu,x,u_exact,grad_u_exact)
    trace_disp = 0.d0
    do i = 1,3
        trace_disp = trace_disp + (u_exact(i) - u(i))**2
    enddo
    error_L2 = error_L2 + w(igpt)*trace_disp*detJ_map
    trace_disp_grad = 0.d0
    do i = 1,3
      do j = 1,3
        trace_disp_grad = trace_disp_grad + (grad_u_exact(i,j) - grad_u(i,j))**2
      enddo
    enddo
    error_H1 = error_H1 + w(igpt)*trace_disp_grad*detJ_map
    call ShapeFunctionDestroy(shapefunction)
  enddo

  deallocate(B)

end subroutine GeomechForceLocalElemError
#endif

! ************************************************************************** !

subroutine GetAnalytical(load_type,lambda,mu,coord,u,grad_u)
  ! 
  ! GeomechGetBodyForce: Gets the body force at a given position
  ! of the point
  ! 
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! 

  PetscReal :: lambda, mu
  PetscReal :: coord(THREE_INTEGER)
  PetscReal :: u(THREE_INTEGER)
  PetscReal :: grad_u(THREE_INTEGER,THREE_INTEGER)
  PetscInt :: load_type
  PetscReal :: x, y, z
  
  x = coord(1)
  y = coord(2)
  z = coord(3)
  
  select case(load_type)
    case(2)
      u(1) = 2*y*(x+y+z)
      u(2) = 4*x-y**2-z**2
      u(3) = sin(PI*x)*sin(PI*y)*sin(PI*z)
      grad_u(1,1) = 2*y
      grad_u(1,2) = 2*x + 4*y + 2*z
      grad_u(1,3) = 2*y
      grad_u(2,1) = 4
      grad_u(2,2) = (-2)*y
      grad_u(2,3) = (-2)*z
      grad_u(3,1) = PI*cos(PI*x)*sin(PI*y)*sin(PI*z)
      grad_u(3,2) = PI*cos(PI*y)*sin(PI*x)*sin(PI*z)
      grad_u(3,3) = PI*cos(PI*z)*sin(PI*x)*sin(PI*y)
    case default
  end select
  
end subroutine GetAnalytical

! ************************************************************************** !

subroutine GeomechForceLocalElemJacobian(size_elenodes,local_coordinates, &
                                         local_disp, &
                                         local_youngs,local_poissons, &
                                         eletype,dim,r,w,Kmat,option)
  ! 
  ! Computes the Jacobian for a local element
  ! 
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! 
                                         
  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module
  
  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: local_disp(:)
  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: igpt
  PetscInt :: len_w
  PetscInt :: eletype
  PetscReal :: x(THREE_INTEGER), J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: detJ_map
  PetscInt :: i,j,d
  PetscReal :: eye_three(THREE_INTEGER)
  PetscInt :: indx(THREE_INTEGER)
  PetscInt :: dim
  PetscReal :: lambda, mu
  PetscReal :: youngs_mod, poissons_ratio
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER)
  PetscReal, allocatable :: N(:,:)
  PetscReal, allocatable :: vecB_transpose(:,:)
  PetscReal, allocatable :: kron_B_eye(:,:)
  PetscReal, allocatable :: kron_B_transpose_eye(:,:)
  PetscReal, allocatable :: Trans(:,:)
  PetscReal, allocatable :: kron_eye_B_transpose(:,:)
  PetscReal, allocatable :: kron_N_eye(:,:)
  PetscReal, allocatable :: local_youngs(:)
  PetscReal, allocatable :: local_poissons(:)
  PetscInt :: size_elenodes

  allocate(B(size_elenodes,dim))
  
  Kmat = 0.d0
  len_w = size(w)
  
  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo
  
  call Transposer(option%ngeomechdof,size_elenodes,Trans)

  do igpt = 1, len_w
    shapefunction%EleType = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction)
    x = matmul(transpose(local_coordinates),shapefunction%N)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    allocate(N(size(shapefunction%N),ONE_INTEGER))
    call Determinant(J_map,detJ_map)
    if (detJ_map <= 0.d0) then
      option%io_buffer = 'GEOMECHANICS: Determinant of J_map has' // &
                         ' to be positive!' 
      call printErrMsg(option)        
    endif
    ! Find the inverse of J_map
    ! Set identity matrix
    call ludcmp(J_map,THREE_INTEGER,indx,d)
    do i = 1, THREE_INTEGER
      eye_three = 0.d0
      eye_three(i) = 1.d0
      call lubksb(J_map,THREE_INTEGER,indx,eye_three)
      inv_J_map(:,i) = eye_three
    enddo
    B = matmul(shapefunction%DN,inv_J_map)
    youngs_mod = dot_product(shapefunction%N,local_youngs)
    poissons_ratio = dot_product(shapefunction%N,local_poissons)
    call GeomechGetLambdaMu(lambda,mu,youngs_mod,poissons_ratio)
    call ConvertMatrixToVector(transpose(B),vecB_transpose)
    Kmat = Kmat + w(igpt)*lambda* &
      matmul(vecB_transpose,transpose(vecB_transpose))*detJ_map
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    N(:,1)= shapefunction%N    
    call Kron(N,identity,kron_N_eye)
    Kmat = Kmat + w(igpt)*mu* &
      matmul(kron_B_eye,kron_B_transpose_eye)*detJ_map
    Kmat = Kmat + w(igpt)*mu* &
      matmul(matmul(kron_B_eye,kron_eye_B_transpose),Trans)*detJ_map
    call ShapeFunctionDestroy(shapefunction)
    deallocate(N)
    deallocate(vecB_transpose)
    deallocate(kron_B_eye)
    deallocate(kron_B_transpose_eye)
    deallocate(kron_eye_B_transpose)
    deallocate(kron_N_eye)
  enddo
    
  deallocate(B)
  deallocate(Trans)

end subroutine GeomechForceLocalElemJacobian

! ************************************************************************** !

subroutine GeomechGetLambdaMu(lambda,mu,E,nu)
  ! 
  ! Gets the material properties given the position
  ! of the point
  ! 
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! 

  PetscReal :: lambda, mu
  PetscReal :: E, nu
  PetscReal :: coord(THREE_INTEGER)
 
  lambda = E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
  mu = E/2.d0/(1.d0+nu)


end subroutine GeomechGetLambdaMu

! ************************************************************************** !

subroutine GeomechGetBodyForce(load_type,lambda,mu,coord,bf,option)
  ! 
  ! Gets the body force at a given position
  ! of the point
  ! 
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! 

  use Option_module

  type(option_type) :: option

  PetscInt :: load_type
  PetscReal :: lambda, mu, den_rock
  PetscReal :: coord(THREE_INTEGER)
  PetscReal :: bf(THREE_INTEGER)
  PetscReal :: x, y, z
 
  bf = 0.d0
  
  x = coord(1)
  y = coord(2)
  z = coord(3)
    
  ! This subroutine needs major changes. For given position, it needs to give 
  ! out lambda, mu, also need to add density of rock
  
  
  select case(load_type)
    case default
      bf(GEOMECH_DISP_X_DOF) = option%geomech_gravity(X_DIRECTION)
      bf(GEOMECH_DISP_Y_DOF) = option%geomech_gravity(Y_DIRECTION)
      bf(GEOMECH_DISP_Z_DOF) = option%geomech_gravity(Z_DIRECTION)
  end select
  
end subroutine GeomechGetBodyForce

! ************************************************************************** !

subroutine GeomechForceJacobian(snes,xx,A,B,geomech_realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Satish Karra
  ! Date: 06/21/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  class(realization_geomech_type) :: geomech_realization
  PetscErrorCode :: ierr
  
  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(geomech_grid_type),  pointer :: grid
  type(option_type), pointer :: option
  PetscReal :: norm
  
  option => geomech_realization%option

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call GeomechForceJacobianPatch(snes,xx,J,J,geomech_realization,ierr)

  if (geomech_realization%geomech_debug%matview_Jacobian) then
    call PetscViewerASCIIOpen(geomech_realization%option%mycomm, &
                              'Geomech_jacobian.out', &
                              viewer,ierr);CHKERRQ(ierr)
   
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (geomech_realization%geomech_debug%norm_Jacobian) then
    option => geomech_realization%option
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

end subroutine GeomechForceJacobian

! ************************************************************************** !

subroutine GeomechForceJacobianPatch(snes,xx,A,B,geomech_realization,ierr)
  ! 
  ! Computes the nonlinear part of the Jacobian on a patch
  ! 
  ! Author: Satish Karra
  ! Date: 06/21/13
  ! Modified: 07/12/16
       
  use Geomechanics_Realization_class
      
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(inout) :: A
  Mat, intent(out) :: B
  PetscViewer :: viewer

  PetscErrorCode :: ierr
   
  class(realization_geomech_type) :: geomech_realization
  
  ! Do nothing here since Jacobian is always linear and is computed
  ! once at the setup of geomechanics realization

end subroutine GeomechForceJacobianPatch  

! ************************************************************************** !

subroutine GeomechForceJacobianLinearPart(A,geomech_realization)
  ! 
  ! Computes the Linear part of the Jacobian on a patch
  ! 
  ! Author: Satish Karra
  ! Date: 06/21/13
  ! Modified: 07/12/16
  ! 
       
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Field_module
  use Geomechanics_Debug_module
  use Geomechanics_Discretization_module
  use Option_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Auxiliary_module
      
  implicit none

  Mat :: A
  PetscViewer :: viewer

  PetscErrorCode :: ierr
   
  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_parameter_type), pointer :: GeomechParam

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:)
  PetscInt, allocatable :: ghosted_ids(:)
  PetscReal, allocatable :: Jac_full(:,:)
  PetscReal, allocatable :: Jac_sub_mat(:,:)
  PetscInt, allocatable :: rows(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscInt :: ielem,ivertex 
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: local_id, petsc_id
  PetscInt :: ghosted_id1, ghosted_id2
  PetscInt :: petsc_id1, petsc_id2
  PetscInt :: id1, id2, i, j, vertex_count, count
  PetscReal, pointer :: imech_loc_p(:)
  PetscInt :: size_elenodes
        
  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars  
  GeomechParam => patch%geomech_aux%GeomechParam 

  call MatZeroEntries(A,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes)*option%ngeomechdof))
    allocate(ghosted_ids(size(elenodes)))
    allocate(Jac_full(size(elenodes)*option%ngeomechdof, &
                      size(elenodes)*option%ngeomechdof))
    allocate(Jac_sub_mat(option%ngeomechdof,option%ngeomechdof))
    allocate(youngs_vec(size(elenodes)))
    allocate(poissons_vec(size(elenodes)))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%EleType
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      ghosted_ids(ivertex) = ghosted_id
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(idof + (ivertex-1)*option%ngeomechdof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
      enddo
      youngs_vec(ivertex) = &
        GeomechParam%youngs_modulus(int(imech_loc_p(ghosted_id))) 
      poissons_vec(ivertex) = &
        GeomechParam%poissons_ratio(int(imech_loc_p(ghosted_id))) 
    enddo
    size_elenodes = size(elenodes)
    call GeomechForceLocalElemJacobian(size_elenodes,local_coordinates, &
       local_disp,youngs_vec,poissons_vec,eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,Jac_full,option)
    do id1 = 1, size(ghosted_ids)
      ghosted_id1 = ghosted_ids(id1)
      petsc_id1 = grid%node_ids_ghosted_petsc(ghosted_id1)
      do id2 = 1, size(ghosted_ids)
        ghosted_id2 = ghosted_ids(id2)
        petsc_id2 = grid%node_ids_ghosted_petsc(ghosted_id2)
        Jac_sub_mat = 0.d0
        Jac_sub_mat =  &
          Jac_full(option%ngeomechdof*(id1-1)+GEOMECH_DISP_X_DOF: &
                   option%ngeomechdof*(id1-1)+GEOMECH_DISP_Z_DOF, &
                   option%ngeomechdof*(id2-1)+GEOMECH_DISP_X_DOF: &
                   option%ngeomechdof*(id2-1)+GEOMECH_DISP_Z_DOF) 
          
        call MatSetValuesBlocked(A,1,petsc_id1-1,1,petsc_id2-1, &
                                 Jac_sub_mat,ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
    enddo
   
    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(ghosted_ids)
    deallocate(Jac_full)
    deallocate(Jac_sub_mat)
    deallocate(youngs_vec)
    deallocate(poissons_vec)
  enddo
  
  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  
  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1
  
  ! Find the number of boundary vertices
  vertex_count = 0
  boundary_condition => patch%geomech_boundary_condition_list%first
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    vertex_count = vertex_count + region%num_verts
    boundary_condition => boundary_condition%next      
  enddo
  
  allocate(rows(vertex_count*option%ngeomechdof))
  count = 0
 
  boundary_condition => patch%geomech_boundary_condition_list%first
  do 
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif    
      
      ! X displacement 
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_X_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Y displacement 
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Y_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
      ! Z displacement      
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Z_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif
      
    enddo
    boundary_condition => boundary_condition%next      
  enddo
    
  call MatZeroRowsLocal(A,count,rows,1.d0, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)
  call MatStoreValues(A,ierr);CHKERRQ(ierr) ! Store the linear part of Jacobian
                    
  deallocate(rows)

end subroutine GeomechForceJacobianLinearPart  

! ************************************************************************** !

subroutine GeomechUpdateFromSubsurf(realization,geomech_realization)
  ! 
  ! The pressure/temperature from subsurface are
  ! mapped to geomech
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/10/13
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Option_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization
  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(geomech_field_type), pointer :: geomech_field
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr
  PetscReal, pointer :: vec_p(:), xx_loc_p(:)
  PetscInt :: local_id, ghosted_id

  option        => realization%option
  grid          => realization%discretization%grid
  field         => realization%field
  geomech_grid  => geomech_realization%geomech_discretization%grid
  geomech_field => geomech_realization%geomech_field

  ! use the subsurface output option parameters for geomechanics as well
  geomech_realization%output_option%tunit = realization%output_option%tunit
  geomech_realization%output_option%tconv = realization%output_option%tconv
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)


  ! pressure
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call GeomechGridVecGetArrayF90(geomech_grid, &
                                 geomech_field%subsurf_vec_1dof,vec_p,ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    vec_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id-1)+1) 
  enddo
  call GeomechGridVecRestoreArrayF90(geomech_grid, &
                                     geomech_field%subsurf_vec_1dof,vec_p,ierr)
  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  
  ! Scatter the data
  call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_field%subsurf_vec_1dof, &
                       geomech_field%press, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                     geomech_field%subsurf_vec_1dof, &
                     geomech_field%press, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
                     
  ! temperature
  if (option%nflowdof > 1) then
    call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
    call GeomechGridVecGetArrayF90(geomech_grid, &
                                   geomech_field%subsurf_vec_1dof,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      vec_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id-1)+2) 
    enddo
    call GeomechGridVecRestoreArrayF90(geomech_grid, &
                                       geomech_field%subsurf_vec_1dof, &
                                       vec_p,ierr)
    call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  
    ! Scatter the data
    call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                         geomech_field%subsurf_vec_1dof, &
                         geomech_field%temp, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_field%subsurf_vec_1dof, &
                       geomech_field%temp, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  endif                       
 
  call GeomechDiscretizationGlobalToLocal(&
                                geomech_realization%geomech_discretization, &
                                geomech_field%press, & 
                                geomech_field%press_loc,ONEDOF)
  
  if (option%nflowdof > 1) &
    call GeomechDiscretizationGlobalToLocal(&
                                geomech_realization%geomech_discretization, &
                                geomech_field%temp, &
                                geomech_field%temp_loc,ONEDOF)

end subroutine GeomechUpdateFromSubsurf

! ************************************************************************** !

subroutine GeomechUpdateSubsurfFromGeomech(realization,geomech_realization)
  ! 
  ! The stresses and strains from geomech
  ! are mapped to subsurf.
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 10/10/13
  ! 

  use Realization_Subsurface_class
  use Discretization_module
  use Grid_module
  use Field_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Option_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization
  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(geomech_field_type), pointer :: geomech_field
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr

  option        => realization%option
  grid          => realization%discretization%grid
  field         => realization%field
  geomech_grid  => geomech_realization%geomech_discretization%grid
  geomech_field => geomech_realization%geomech_field
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)
  
  ! Scatter the strains
  call VecScatterBegin(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%strain, &
                       geomech_field%strain_subsurf, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%strain, &
                       geomech_field%strain_subsurf, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
                       
  ! Scatter the stresses
  call VecScatterBegin(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%stress, &
                       geomech_field%stress_subsurf, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%stress, &
                       geomech_field%stress_subsurf, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
                       
  ! Scatter from global to local vectors
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   geomech_field%strain_subsurf, &
                                   geomech_field%strain_subsurf_loc, &
                                   NGEODOF)
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   geomech_field%stress_subsurf, &
                                   geomech_field%stress_subsurf_loc, &
                                   NGEODOF)
 
end subroutine GeomechUpdateSubsurfFromGeomech

! ************************************************************************** !

subroutine GeomechCreateGeomechSubsurfVec(realization,geomech_realization)
  ! 
  ! Creates the MPI vector that stores the
  ! variables from subsurface
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/10/13
  ! 

  use Grid_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module
  use String_module
  use Realization_Subsurface_class
  use Option_module

  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization

  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: geomech_field
  
  PetscErrorCode :: ierr
  
  option     => realization%option
  grid       => realization%discretization%grid
  geomech_field => geomech_realization%geomech_field
  
  call VecCreate(option%mycomm,geomech_field%subsurf_vec_1dof, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%subsurf_vec_1dof, &
                   grid%nlmax,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%subsurf_vec_1dof,ierr);CHKERRQ(ierr)
  
end subroutine GeomechCreateGeomechSubsurfVec

! ************************************************************************** !

subroutine GeomechCreateSubsurfStressStrainVec(realization,geomech_realization)
  ! 
  ! Creates the subsurface stress and strain
  ! MPI vectors to store information from geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 10/10/13
  ! 

  use Grid_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module
  use String_module
  use Realization_Subsurface_class
  use Option_module

  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization

  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: geomech_field
  
  PetscErrorCode :: ierr
  
  option     => realization%option
  grid       => realization%discretization%grid
  geomech_field => geomech_realization%geomech_field
  
  ! strain
  call VecCreate(option%mycomm,geomech_field%strain_subsurf, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%strain_subsurf, &
                   grid%nlmax*SIX_INTEGER,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%strain_subsurf,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%strain_subsurf,ierr);CHKERRQ(ierr)
  
  ! stress
  call VecCreate(option%mycomm,geomech_field%stress_subsurf, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%stress_subsurf, &
                   grid%nlmax*SIX_INTEGER,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%stress_subsurf,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%stress_subsurf,ierr);CHKERRQ(ierr)
  
  ! strain_loc
  call VecCreate(PETSC_COMM_SELF,geomech_field%strain_subsurf_loc, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%strain_subsurf_loc, &
                   grid%ngmax*SIX_INTEGER,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%strain_subsurf_loc,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%strain_subsurf_loc,ierr);CHKERRQ(ierr)
  
  ! stress_loc 
  call VecCreate(PETSC_COMM_SELF,geomech_field%stress_subsurf_loc, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%stress_subsurf_loc, &
                   grid%ngmax*SIX_INTEGER,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%stress_subsurf_loc,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%stress_subsurf_loc,ierr);CHKERRQ(ierr)
  
end subroutine GeomechCreateSubsurfStressStrainVec

! ************************************************************************** !

subroutine GeomechForceStressStrain(geomech_realization)
  ! 
  ! Computes the stress strain on a patch
  ! 
  ! Author: Satish Karra
  ! Date: 09/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module
  use Geomechanics_Auxiliary_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_parameter_type), pointer :: GeomechParam

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscReal, allocatable :: strain(:,:), stress(:,:)
  PetscInt, allocatable :: count(:)
  PetscInt :: ielem, ivertex 
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: petsc_id, local_id
  PetscInt :: size_elenodes
  PetscReal, pointer :: imech_loc_p(:)  
  PetscReal, pointer :: strain_loc_p(:)
  PetscReal, pointer :: stress_loc_p(:)
  PetscReal, pointer :: strain_p(:), stress_p(:)
  PetscReal, pointer :: no_elems_p(:)
  
  PetscErrorCode :: ierr
                  
  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%GeomechGlobal%aux_vars  
  GeomechParam => patch%geomech_aux%GeomechParam 

  call VecSet(field%strain,0.d0,ierr);CHKERRQ(ierr)
  call VecSet(field%stress,0.d0,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  
  strain_loc_p = 0.d0
  stress_loc_p = 0.d0
  
   ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    allocate(elenodes(grid%elem_nodes(0,ielem)))
    allocate(local_coordinates(size(elenodes),THREE_INTEGER))
    allocate(local_disp(size(elenodes),option%ngeomechdof))
    allocate(petsc_ids(size(elenodes)))
    allocate(ids(size(elenodes)*option%ngeomechdof))
    allocate(youngs_vec(size(elenodes)))
    allocate(poissons_vec(size(elenodes)))
    allocate(strain(size(elenodes),SIX_INTEGER))
    allocate(stress(size(elenodes),SIX_INTEGER))
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    eletype = grid%gauss_node(ielem)%EleType
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(ivertex,idof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
      youngs_vec(ivertex) = &
        GeomechParam%youngs_modulus(int(imech_loc_p(ghosted_id))) 
      poissons_vec(ivertex) = &
        GeomechParam%poissons_ratio(int(imech_loc_p(ghosted_id))) 
    enddo
    size_elenodes = size(elenodes)
    call GeomechForceLocalElemStressStrain(size_elenodes,local_coordinates, &
       local_disp,youngs_vec,poissons_vec, &
       eletype,grid%gauss_node(ielem)%dim,strain,stress,option)
 
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      do idof = 1, SIX_INTEGER
        strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + &
          strain(ivertex,idof)
        stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + & 
          stress(ivertex,idof)
      enddo
    enddo
   
    deallocate(elenodes)
    deallocate(local_coordinates)
    deallocate(local_disp)
    deallocate(petsc_ids)
    deallocate(ids)
    deallocate(youngs_vec)
    deallocate(poissons_vec)
    deallocate(strain)
    deallocate(stress)
  enddo

  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  
  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%strain_loc,field%strain, &
                                             SIX_INTEGER)
  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%stress_loc,field%stress, &
                                             SIX_INTEGER)
                                             
! Now take the average at each node for elements sharing the node
  call VecGetArrayF90(grid%no_elems_sharing_node,no_elems_p, &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%strain,strain_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress,stress_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax_node
    do idof = 1, SIX_INTEGER
      strain_p(idof + (local_id-1)*SIX_INTEGER) = &
        strain_p(idof + (local_id-1)*SIX_INTEGER)/int(no_elems_p(local_id))
      stress_p(idof + (local_id-1)*SIX_INTEGER) = &
        stress_p(idof + (local_id-1)*SIX_INTEGER)/int(no_elems_p(local_id))
    enddo
  enddo
  call VecRestoreArrayF90(field%stress,stress_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%strain,strain_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%no_elems_sharing_node,no_elems_p, &
                          ierr);CHKERRQ(ierr)

! Now scatter back to local domains
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%strain,field%strain_loc, &
                                          SIX_INTEGER)
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%stress,field%stress_loc, &
                                          SIX_INTEGER)
                                          
  call VecGetArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
! Copy them to global_aux_vars
  do ghosted_id = 1, grid%ngmax_node  
    do idof = 1, SIX_INTEGER
      geomech_global_aux_vars(ghosted_id)%strain(idof) = &
        strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
      geomech_global_aux_vars(ghosted_id)%stress(idof) = &
        stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
    enddo
  enddo
  call VecRestoreArrayF90(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)

end subroutine GeomechForceStressStrain

! ************************************************************************** !

subroutine GeomechForceLocalElemStressStrain(size_elenodes,local_coordinates, &
                                             local_disp, &
                                             local_youngs,local_poissons, &
                                             eletype,dim,strain,stress,option)
  ! 
  ! Computes the stress-strain for a local
  ! element
  ! 
  ! Author: Satish Karra
  ! Date: 09/17/13
  ! 
                                         
  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module
  
  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: B(:,:), Kmat(:,:)
  PetscReal, allocatable :: res_vec(:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscReal, allocatable :: local_youngs(:)
  PetscReal, allocatable :: local_poissons(:)
  PetscReal, allocatable :: strain(:,:)
  PetscReal, allocatable :: stress(:,:) 
  PetscReal :: strain_local(NINE_INTEGER,ONE_INTEGER)
  PetscReal :: stress_local(NINE_INTEGER,ONE_INTEGER) 
    
  PetscReal, pointer :: r(:,:), w(:)
  PetscInt :: ivertex
  PetscInt :: eletype
  PetscReal :: identity(THREE_INTEGER,THREE_INTEGER) 
  PetscInt :: indx(THREE_INTEGER)
  PetscInt :: dim
  PetscInt :: i, j, d
  PetscReal :: lambda, mu
  PetscReal :: youngs_mod, poissons_ratio
  PetscReal, allocatable :: kron_B_eye(:,:)
  PetscReal, allocatable :: kron_B_transpose_eye(:,:)
  PetscReal, allocatable :: Trans(:,:)
  PetscReal, allocatable :: kron_eye_B_transpose(:,:)
  PetscReal, allocatable :: vec_local_disp(:,:)
  PetscInt :: size_elenodes
  PetscReal :: J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: inv_J_map(THREE_INTEGER,THREE_INTEGER)
  PetscReal :: eye_three(THREE_INTEGER)
  PetscReal :: eye_vec(NINE_INTEGER,ONE_INTEGER) 
  
  allocate(B(size_elenodes,dim))
  
  call Transposer(option%ngeomechdof,size_elenodes,Trans)
  strain = 0.d0
  stress = 0.d0 

  call ConvertMatrixToVector(transpose(local_disp),vec_local_disp)

  identity = 0.d0
  do i = 1, THREE_INTEGER
    do j = 1, THREE_INTEGER
      if (i == j) identity(i,j) = 1.d0
    enddo
  enddo

  eye_vec = 0.d0
  eye_vec(1,1) = 1.d0
  eye_vec(5,1) = 1.d0
  eye_vec(9,1) = 1.d0
 
  do ivertex = 1, size_elenodes
    strain_local = 0.d0
    stress_local = 0.d0 
    shapefunction%EleType = eletype
    call ShapeFunctionInitialize(shapefunction)
    shapefunction%zeta = shapefunction%coord(ivertex,:)
    call ShapeFunctionCalculate(shapefunction)
    J_map = matmul(transpose(local_coordinates),shapefunction%DN)
    call ludcmp(J_map,THREE_INTEGER,indx,d)
    do i = 1, THREE_INTEGER
      eye_three = 0.d0
      eye_three(i) = 1.d0
      call lubksb(J_map,THREE_INTEGER,indx,eye_three)
      inv_J_map(:,i) = eye_three
    enddo
    B = matmul(shapefunction%DN,inv_J_map)
    youngs_mod = dot_product(shapefunction%N,local_youngs)
    poissons_ratio = dot_product(shapefunction%N,local_poissons)
    call GeomechGetLambdaMu(lambda,mu,youngs_mod,poissons_ratio)
    call Kron(B,identity,kron_B_eye)
    call Kron(transpose(B),identity,kron_B_transpose_eye)
    call Kron(identity,transpose(B),kron_eye_B_transpose)
    strain_local =  0.5*matmul((kron_B_transpose_eye + &
      matmul(kron_eye_B_transpose,Trans)),vec_local_disp)
    stress_local = lambda*(strain_local(1,1)+ &
                   strain_local(5,1)+strain_local(9,1))*eye_vec + &
                   2*mu*strain_local 
    call ShapeFunctionDestroy(shapefunction)
    deallocate(kron_B_eye)
    deallocate(kron_B_transpose_eye)
    deallocate(kron_eye_B_transpose)
    strain(ivertex,1) = strain_local(1,1)
    strain(ivertex,2) = strain_local(5,1)
    strain(ivertex,3) = strain_local(9,1)
    strain(ivertex,4) = strain_local(2,1)
    strain(ivertex,5) = strain_local(6,1)
    strain(ivertex,6) = strain_local(3,1)
    stress(ivertex,1) = stress_local(1,1)
    stress(ivertex,2) = stress_local(5,1)
    stress(ivertex,3) = stress_local(9,1)
    stress(ivertex,4) = stress_local(2,1)
    stress(ivertex,5) = stress_local(6,1)
    stress(ivertex,6) = stress_local(3,1)
  enddo
  
  deallocate(B)
  deallocate(vec_local_disp)
  deallocate(Trans)

end subroutine GeomechForceLocalElemStressStrain 

! ************************************************************************** !

subroutine GeomechUpdateSolution(geomech_realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  
  implicit none 
  
  class(realization_geomech_type) :: geomech_realization
  type(geomech_field_type), pointer :: field
  
  PetscErrorCode :: ierr 
  PetscViewer :: viewer
  
  field => geomech_realization%geomech_field

  call GeomechUpdateSolutionPatch(geomech_realization)

end subroutine GeomechUpdateSolution

! ************************************************************************** !

subroutine GeomechUpdateSolutionPatch(geomech_realization)
  ! 
  ! updates data in module after a successful time
  ! step
  ! 
  ! Author: satish karra, lanl
  ! Date: 09/17/13
  ! 

  use Geomechanics_Realization_class
    
  implicit none 
  
  class(realization_geomech_type) :: geomech_realization

  call GeomechForceStressStrain(geomech_realization)

end subroutine GeomechUpdateSolutionPatch

! ************************************************************************** !

subroutine GeomechStoreInitialPressTemp(geomech_realization)
  ! 
  ! Stores initial pressure and temperature from
  ! subsurface
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/24/13
  ! 

  use Geomechanics_Realization_class
    
  implicit none 
  
  class(realization_geomech_type) :: geomech_realization

  PetscErrorCode :: ierr

  call VecCopy(geomech_realization%geomech_field%press_loc, & 
               geomech_realization%geomech_field%press_init_loc, &
               ierr);CHKERRQ(ierr)
 
  call VecCopy(geomech_realization%geomech_field%temp_loc, & 
               geomech_realization%geomech_field%temp_init_loc, &
               ierr);CHKERRQ(ierr)
   
end subroutine GeomechStoreInitialPressTemp

! ************************************************************************** !

subroutine GeomechStoreInitialPorosity(realization,geomech_realization)
  ! 
  ! Stores initial porosity from
  ! subsurface
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 10/22/13
  ! 

  use Geomechanics_Realization_class
  use Realization_Subsurface_class
  use Discretization_module
    
  implicit none 
  
  class(realization_geomech_type) :: geomech_realization
  class(realization_subsurface_type) :: realization
  type(discretization_type) :: discretization

  PetscErrorCode :: ierr

  call DiscretizationDuplicateVector(discretization, &
                                     realization%field%work_loc, &
                                     geomech_realization%geomech_field% &
                                     porosity_init_loc)
   
end subroutine GeomechStoreInitialPorosity

! ************************************************************************** !

subroutine GeomechStoreInitialDisp(geomech_realization)
  ! 
  ! Stores initial displacement for calculating
  ! relative displacements
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/30/13
  ! 

  use Geomechanics_Realization_class
    
  implicit none 
  
  class(realization_geomech_type) :: geomech_realization

  PetscErrorCode :: ierr

  call VecCopy(geomech_realization%geomech_field%disp_xx_loc, & 
               geomech_realization%geomech_field%disp_xx_init_loc, &
               ierr);CHKERRQ(ierr)
   
end subroutine GeomechStoreInitialDisp

! ************************************************************************** !

subroutine GeomechUpdateSubsurfPorosity(realization,geomech_realization)
  ! 
  ! Updates the porosity in the subsurface based
  ! on the deformation in geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 10/08/13
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Discretization_module
  use Geomechanics_Field_module
  use Material_Aux_class
  use Material_module
  use Variables_module, only : POROSITY
  use Geomechanics_Realization_class

  implicit none
  
  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(geomech_field_type), pointer :: geomech_field
  type(grid_type), pointer :: grid
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscReal :: trace_epsilon
  PetscReal, pointer :: por0_loc_p(:), strain_loc_p(:)
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  geomech_field => geomech_realization%geomech_field
  material_auxvars => realization%patch%aux%Material%auxvars

  if (.not.associated(patch%imat)) then
    option%io_buffer = 'Materials IDs not present in run.  Material ' // &
      ' properties cannot be updated without material ids'
    call printErrMsg(option)
  endif
  
  call VecGetArrayF90(geomech_field%porosity_init_loc,por0_loc_p, &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(geomech_field%strain_subsurf_loc,strain_loc_p, &
                      ierr);CHKERRQ(ierr)
  
  do ghosted_id = 1, grid%ngmax
    trace_epsilon = strain_loc_p((ghosted_id-1)*SIX_INTEGER+ONE_INTEGER) + &
                    strain_loc_p((ghosted_id-1)*SIX_INTEGER+TWO_INTEGER) + &
                    strain_loc_p((ghosted_id-1)*SIX_INTEGER+THREE_INTEGER)
    material_auxvars(ghosted_id)%porosity = por0_loc_p(ghosted_id)/ &
      (1.d0 + (1.d0 - por0_loc_p(ghosted_id))*trace_epsilon)
  enddo
  
  call VecRestoreArrayF90(geomech_field%porosity_init_loc,por0_loc_p, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(geomech_field%strain_subsurf_loc,strain_loc_p, &
                          ierr);CHKERRQ(ierr)

  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,ZERO_INTEGER)
  call DiscretizationLocalToLocal(realization%discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,ZERO_INTEGER)

end subroutine GeomechUpdateSubsurfPorosity

end module Geomechanics_Force_module
