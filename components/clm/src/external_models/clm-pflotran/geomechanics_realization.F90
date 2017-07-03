module Geomechanics_Realization_class

  use Realization_Base_class
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Material_module
  use Geomechanics_Field_module
  use Geomechanics_Debug_module
  use Geomechanics_Region_module
  use Geomechanics_Condition_module
  use Input_Aux_module
  use Option_module
  use Output_Aux_module
  use Dataset_Base_class
  use PFLOTRAN_Constants_module

  implicit none

private

#include "petsc/finclude/petscsys.h"

  type, public, extends(realization_base_type) :: realization_geomech_type

    type(geomech_discretization_type), pointer :: geomech_discretization
    type(geomech_patch_type), pointer :: geomech_patch

    type(geomech_material_property_type), &
                           pointer :: geomech_material_properties
    type(geomech_material_property_ptr_type), &
                           pointer :: geomech_material_property_array(:)

    type(geomech_field_type), pointer :: geomech_field
    type(geomech_debug_type), pointer :: geomech_debug
    type(gm_region_list_type), pointer :: geomech_region_list
    type(geomech_condition_list_type),pointer :: geomech_conditions
    class(dataset_base_type), pointer :: geomech_datasets
    PetscReal :: dt_coupling

  end type realization_geomech_type

  public :: GeomechRealizCreate, &
            GeomechRealizDestroy, &
            GeomechRealizAddStrata, &
            GeomechRealizAddGeomechCoupler, &
            GeomechRealizLocalizeRegions, &
            GeomechRealizPassFieldPtrToPatch, &
            GeomechRealizProcessMatProp, &
            GeomechRealizProcessGeomechCouplers, &
            GeomechRealizCreateDiscretization, &
            GeomechRealizProcessGeomechConditions, &
            GeomechRealizInitAllCouplerAuxVars, &
            GeomechRealizPrintCouplers, &
            GeomechRealizAddWaypointsToList, &
            GeomechRealizGetDataset, &
            GeomechRealizLocalToLocalWithArray, &
            GeomechRealizMapSubsurfGeomechGrid, &
            GeomechGridElemSharedByNodes
contains

! ************************************************************************** !

function GeomechRealizCreate(option)
  ! 
  ! This subroutine creates realization for geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  ! 

  implicit none

  class(realization_geomech_type), pointer :: GeomechRealizCreate
  class(realization_geomech_type), pointer :: geomech_realization
  type(option_type), pointer :: option
  
  allocate(geomech_realization)
  geomech_realization%id = 0
  if (associated(option)) then
    geomech_realization%option => option
  else
    geomech_realization%option => OptionCreate()
  endif
  
  nullify(geomech_realization%input)
  geomech_realization%geomech_discretization => GeomechDiscretizationCreate()
  
  geomech_realization%geomech_field => GeomechFieldCreate()
  geomech_realization%output_option => OutputOptionCreate()
  geomech_realization%geomech_debug => GeomechDebugCreate()
  
  allocate(geomech_realization%geomech_region_list)
  call GeomechRegionInitList(geomech_realization%geomech_region_list)
  
  allocate(geomech_realization%geomech_conditions)
  call GeomechConditionInitList(geomech_realization%geomech_conditions)

  nullify(geomech_realization%geomech_material_properties)
  nullify(geomech_realization%geomech_material_property_array)
  
  nullify(geomech_realization%geomech_patch)
  geomech_realization%dt_coupling = 0.d0

  GeomechRealizCreate => geomech_realization
  
end function GeomechRealizCreate

! ************************************************************************** !

subroutine GeomechRealizAddStrata(geomech_realization,strata)
  ! 
  ! Adds strata to a list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  ! 

  use Geomechanics_Strata_module

  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  type(geomech_strata_type), pointer :: strata
  
  type(geomech_patch_type), pointer :: geomech_patch
  type(geomech_strata_type), pointer :: new_strata
  
  geomech_patch => geomech_realization%geomech_patch

  if (.not.associated(geomech_patch)) return
 
  new_strata => GeomechStrataCreate(strata)
  call GeomechStrataAddToList(new_strata,geomech_patch%geomech_strata_list)
  nullify(new_strata)
  
  call GeomechStrataDestroy(strata)
 
end subroutine GeomechRealizAddStrata

! ************************************************************************** !

subroutine GeomechRealizLocalizeRegions(geomech_realization)
  ! 
  ! This routine localizes geomechanics regions
  ! within each patch
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/07/13
  ! 

  use Option_module
  use String_module

  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  type(geomech_patch_type), pointer :: patch
  type(option_type), pointer :: option

  option => geomech_realization%option

  ! localize the regions on each patch
  patch => geomech_realization%geomech_patch
  call GeomechPatchLocalizeRegions(patch, &
                                   geomech_realization%geomech_region_list, &
                                   option)
                                   
end subroutine GeomechRealizLocalizeRegions

! ************************************************************************** !

subroutine GeomechRealizProcessMatProp(geomech_realization)
  ! 
  ! Setup for material properties
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  use String_module
  
  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  type(geomech_patch_type), pointer :: patch  
  type(option_type), pointer :: option

  
  option => geomech_realization%option
  
  ! organize lists
  call GeomechanicsMaterialPropConvertListToArray( &
                        geomech_realization%geomech_material_properties, &
                        geomech_realization%geomech_material_property_array, &
                        option)
  ! set up mirrored pointer arrays within patches to saturation functions
  ! and material properties
  patch => geomech_realization%geomech_patch
  patch%geomech_material_properties => geomech_realization% &
                                       geomech_material_properties
  call GeomechanicsMaterialPropConvertListToArray( &
                                    patch%geomech_material_properties, &
                                    patch%geomech_material_property_array, &
                                    option)
                                      
end subroutine GeomechRealizProcessMatProp

! ************************************************************************** !

subroutine GeomechRealizCreateDiscretization(geomech_realization)
  ! 
  ! Creates grid
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  ! 

  use Geomechanics_Grid_Aux_module
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: geomech_field
  type(gmdm_ptr_type), pointer :: dm_ptr
  PetscErrorCode :: ierr

  geomech_discretization => geomech_realization%geomech_discretization
  grid => geomech_discretization%grid
  option => geomech_realization%option
  geomech_field => geomech_realization%geomech_field
  
  call GeomechDiscretizationCreateDMs(geomech_discretization,option)
  
  ! n degree of freedom, global
  call GeomechDiscretizationCreateVector(geomech_discretization,NGEODOF, &
                                         geomech_field%disp_xx, &
                                         GLOBAL,option)
  call VecSet(geomech_field%disp_xx,0.d0,ierr);CHKERRQ(ierr)

  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%disp_xx, &
                                            geomech_field%disp_r)
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%disp_xx, &
                                            geomech_field%work)
  
  ! 1 degree of freedom, global                                                                                    
  call GeomechDiscretizationCreateVector(geomech_discretization,ONEDOF, &
                                         geomech_field%press, &
                                         GLOBAL,option)
  call VecSet(geomech_field%press,0.d0,ierr);CHKERRQ(ierr)
  
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%press, &
                                            geomech_field%temp)
                                            
  ! n degrees of freedom, local
  call GeomechDiscretizationCreateVector(geomech_discretization,NGEODOF, &
                                         geomech_field%disp_xx_loc, &
                                         LOCAL,option)
  call VecSet(geomech_field%disp_xx_loc,0.d0,ierr);CHKERRQ(ierr)
 
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%disp_xx_loc, &
                                            geomech_field%work_loc)

  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%disp_xx_loc, &
                                            geomech_field%disp_xx_init_loc)
                                            
  ! 1 degree of freedom, local
  call GeomechDiscretizationCreateVector(geomech_discretization,ONEDOF, &
                                         geomech_field%press_loc, &
                                         LOCAL,option)

  call VecSet(geomech_field%press_loc,0.d0,ierr);CHKERRQ(ierr)
  
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%press_loc, &
                                            geomech_field%temp_loc)

  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%press_loc, &
                                            geomech_field%press_init_loc)

  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%press_loc, &
                                            geomech_field%temp_init_loc)

  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%press_loc, &
                                            geomech_field%imech_loc)

  ! 6 dof for strain and stress
  call GeomechDiscretizationCreateVector(geomech_discretization,SIX_INTEGER, &
                                         geomech_field%strain_loc, &
                                         LOCAL,option)

  call VecSet(geomech_field%strain_loc,0.d0,ierr);CHKERRQ(ierr)
 
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%strain_loc, &
                                            geomech_field%stress_loc)

  call GeomechDiscretizationCreateVector(geomech_discretization,SIX_INTEGER, &
                                         geomech_field%strain, &
                                         GLOBAL,option)

  call VecSet(geomech_field%strain,0.d0,ierr);CHKERRQ(ierr)
 
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            geomech_field%strain, &
                                            geomech_field%stress) 

  grid => geomech_discretization%grid
  
  ! set up nG2L, NL2G, etc.
  call GMGridMapIndices(grid,geomech_discretization%dm_1dof%gmdm, &
                        grid%nG2L,grid%nL2G,grid%nG2A,option)
                        
  ! SK, Need to add a subroutine to ensure right hand rule
  ! SK, Need to add a subroutine equivalent to UGridComputeCoord                      

  
end subroutine GeomechRealizCreateDiscretization

! ************************************************************************** !

subroutine GeomechRealizMapSubsurfGeomechGrid(realization, &
                                              geomech_realization, &
                                              option)
  ! 
  ! This routine creates scatter contexts
  ! betweeen subsurface and geomech grids
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/09/13
  ! 

  use Option_module
  use Geomechanics_Grid_Aux_module
  use Realization_Subsurface_class
  use Grid_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscdm.h"  
#include "petsc/finclude/petscdm.h90"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscviewer.h"

  class(realization_subsurface_type), pointer :: realization
  class(realization_geomech_type), pointer :: geomech_realization
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type) :: option
  type(grid_type), pointer :: grid
  type(gmdm_type), pointer :: gmdm
  type(gmdm_ptr_type), pointer :: dm_ptr
  IS :: is_geomech, is_subsurf
  IS :: is_subsurf_natural
  IS :: is_subsurf_petsc
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  AO :: ao_geomech_to_subsurf_natural
  AO :: ao_subsurf_natual_to_petsc
  PetscInt, allocatable :: int_array(:)
  PetscInt :: local_id
  VecScatter :: scatter
  IS :: is_geomech_petsc
  PetscInt, pointer :: int_ptr(:)
  IS :: is_geomech_petsc_block
  IS :: is_subsurf_petsc_block

  geomech_grid => geomech_realization%geomech_discretization%grid
  grid => realization%discretization%grid
    
  ! Convert from 1-based to 0-based  
  ! Create IS for flow side cell ids
  call ISCreateGeneral(option%mycomm,geomech_grid%mapping_num_cells, &
                       geomech_grid%mapping_cell_ids_flow-1, &
                       PETSC_COPY_VALUES,is_subsurf,ierr);CHKERRQ(ierr)
                       
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_mapping_cell_ids_flow.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_subsurf,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif   

  ! Convert from 1-based to 0-based
  ! Create IS for geomech side vertex ids
  call ISCreateGeneral(option%mycomm,geomech_grid%mapping_num_cells, &
                       geomech_grid%mapping_vertex_ids_geomech-1, &
                       PETSC_COPY_VALUES,is_geomech,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_mapping_vertex_ids_geomech.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_geomech,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif   

  ! Create an application ordering between flow cell ids and geomech vertex ids
  call AOCreateMappingIS(is_geomech,is_subsurf,ao_geomech_to_subsurf_natural, &
                         ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_ao_geomech_to_subsurf_natural.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call AOView(ao_geomech_to_subsurf_natural,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  
  
  allocate(int_array(grid%nlmax))
  do local_id = 1, grid%nlmax
    int_array(local_id) = grid%nG2A(grid%nL2G(local_id)) - 1
  enddo

  ! Flow natural numbering IS
  call ISCreateGeneral(option%mycomm,grid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_subsurf_natural, &
                       ierr);CHKERRQ(ierr)
  deallocate(int_array)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_subsurf_natural.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_subsurf_natural,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  allocate(int_array(grid%nlmax))
  do local_id = 1, grid%nlmax
    int_array(local_id) = (local_id-1) + grid%global_offset
  enddo

  ! Flow petsc numbering IS
  call ISCreateGeneral(option%mycomm,grid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_subsurf_petsc, &
                       ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_subsurf_petsc.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_subsurf_natural,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  ! AO for flow natural to petsc numbering
  call AOCreateMappingIS(is_subsurf_natural,is_subsurf_petsc, &
                         ao_subsurf_natual_to_petsc, &
                         ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_ao_subsurf_natural_to_petsc.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call AOView(ao_subsurf_natual_to_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call AOApplicationToPetscIS(ao_subsurf_natual_to_petsc, &
                              is_subsurf,ierr);CHKERRQ(ierr)


  call ISDuplicate(is_geomech,is_geomech_petsc,ierr);CHKERRQ(ierr)
  call ISCopy(is_geomech,is_geomech_petsc,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_geomech_petsc.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_geomech_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  ! Now calculate the petsc ordering of the mapped geomech nodes
  call AOApplicationToPetscIS(geomech_grid%ao_natural_to_petsc_nodes, &
                              is_geomech_petsc,ierr);CHKERRQ(ierr)
                              
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_subsurf_petsc_geomech_petsc.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_geomech_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif                              

  ! Create scatter context between flow and geomech
  call VecScatterCreate(realization%field%porosity0,is_subsurf, &
                        geomech_realization%geomech_field%press, &
                        is_geomech_petsc,scatter,ierr);CHKERRQ(ierr)
                        
  if (ierr /= 0) then
    option%io_buffer = 'The number of cells specified in ' // &
                       'input file might not be same as the ' // &
                       'SUBSURF->GEOMECH mapping used.'
    call printErrMsg(option)
  endif

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_scatter_subsurf_to_geomech.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecScatterView(scatter,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)

  call VecScatterCopy(scatter,dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                      ierr);CHKERRQ(ierr)
  
  call VecScatterDestroy(scatter,ierr);CHKERRQ(ierr)
  
  ! Geomech to subsurf scatter
  
  allocate(int_array(grid%nlmax))
  call ISGetIndicesF90(is_geomech_petsc,int_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    int_array(local_id) = int_ptr(local_id)
  enddo  
  call ISRestoreIndicesF90(is_geomech_petsc,int_ptr,ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,SIX_INTEGER,grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,is_geomech_petsc_block, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_geomech_petsc_block.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_geomech_petsc_block,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif 

  allocate(int_array(grid%nlmax))
  call ISGetIndicesF90(is_subsurf_petsc,int_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    int_array(local_id) = int_ptr(local_id)
  enddo  
  call ISRestoreIndicesF90(is_subsurf_petsc,int_ptr,ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,SIX_INTEGER,grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,is_subsurf_petsc_block, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_subsurf_petsc_block.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_subsurf_petsc_block,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  
  
  call VecScatterCreate(geomech_realization%geomech_field%strain, &
                        is_geomech_petsc_block, &
                        geomech_realization%geomech_field%strain_subsurf, &
                        is_subsurf_petsc_block,scatter,ierr);CHKERRQ(ierr)
                        
  if (ierr /= 0) then
    option%io_buffer = 'The number of cells specified in ' // &
                       'input file might not be same as the ' // &
                       'GEOMECH->SUBSURF mapping used.'
    call printErrMsg(option)
  endif

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_scatter_geomech_to_subsurf_block.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecScatterView(scatter,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)

  call VecScatterCopy(scatter,dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                      ierr);CHKERRQ(ierr)
 
  call VecScatterDestroy(scatter,ierr);CHKERRQ(ierr) 
  call ISDestroy(is_geomech,ierr);CHKERRQ(ierr)
  call ISDestroy(is_subsurf,ierr);CHKERRQ(ierr)
  call ISDestroy(is_subsurf_natural,ierr);CHKERRQ(ierr)
  call ISDestroy(is_geomech_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(is_subsurf_petsc,ierr);CHKERRQ(ierr)
  call AODestroy(ao_geomech_to_subsurf_natural,ierr);CHKERRQ(ierr)
  call AODestroy(ao_subsurf_natual_to_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(is_subsurf_petsc_block,ierr);CHKERRQ(ierr)
  call ISDestroy(is_geomech_petsc_block,ierr);CHKERRQ(ierr)

end subroutine GeomechRealizMapSubsurfGeomechGrid

! ************************************************************************** !

subroutine GeomechGridElemSharedByNodes(geomech_realization,option)
  ! 
  ! GeomechGridElemsSharedByNodes: Calculates the number of elements common
  ! to a node (vertex)
  ! 
  ! Author: Satish Karra
  ! Date: 09/17/13
  ! 
  
  use Option_module
  use Geomechanics_Grid_Aux_module

  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_geomech_type) :: geomech_realization
  type(geomech_grid_type), pointer :: grid
  type(option_type) :: option
  
  PetscInt :: ielem
  PetscInt :: ivertex
  PetscInt :: ghosted_id
  PetscInt :: elenodes(10)
  PetscReal, pointer :: elem_sharing_node_loc_p(:)
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  
  grid => geomech_realization%geomech_discretization%grid
  
  call VecGetArrayF90(grid%no_elems_sharing_node_loc,elem_sharing_node_loc_p, &
                      ierr);CHKERRQ(ierr)
  
  ! Calculate the common elements to a node on a process
  do ielem = 1, grid%nlmax_elem
    elenodes(1:grid%elem_nodes(0,ielem)) = &
      grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex)
      elem_sharing_node_loc_p(ghosted_id) = &
        elem_sharing_node_loc_p(ghosted_id) + 1
    enddo
  enddo
    
  call VecRestoreArrayF90(grid%no_elems_sharing_node_loc, &
                          elem_sharing_node_loc_p,ierr);CHKERRQ(ierr)
                          
  ! Local to global scatter
  call GeomechDiscretizationLocalToGlobalAdd(&
                                geomech_realization%geomech_discretization, &
                                grid%no_elems_sharing_node_loc, &
                                grid%no_elems_sharing_node, &
                                ONEDOF)  
                                             
#if GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'no_elems_sharing_node_loc_' // trim(adjustl(string)) // '.out'

  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string), &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(grid%no_elems_sharing_node_loc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(option%mycomm,'no_elems_sharing_node.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(grid%no_elems_sharing_node,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

end subroutine GeomechGridElemSharedByNodes

! ************************************************************************** !

subroutine GeomechRealizInitAllCouplerAuxVars(geomech_realization)
  ! 
  ! This routine initializez coupler
  ! auxillary variables
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Option_module

  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  
  type(geomech_patch_type), pointer :: patch
  
  patch => geomech_realization%geomech_patch

  call GeomechPatchInitAllCouplerAuxVars(patch,geomech_realization%option)

end subroutine GeomechRealizInitAllCouplerAuxVars

! ************************************************************************** !

subroutine GeomechRealizLocalToLocalWithArray(geomech_realization,array_id)
  ! 
  ! This routine takes an F90 array that is
  ! ghosted and updates the ghosted values
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  PetscInt :: array_id
  
  type(geomech_patch_type), pointer :: patch
  type(geomech_grid_type), pointer :: grid
  type(geomech_field_type), pointer :: geomech_field

  geomech_field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid

  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GeomechGridCopyIntegerArrayToVec(grid,patch%imat, &
                                            geomech_field%work_loc, &
                                            grid%ngmax_node)
  end select
    
  call GeomechDiscretizationLocalToLocal(& 
                            geomech_realization%geomech_discretization, &
                            geomech_field%work_loc, &
                            geomech_field%work_loc,ONEDOF)
                                  
  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GeomechGridCopyVecToIntegerArray(grid,patch%imat, &
                                            geomech_field%work_loc, &
                                            grid%ngmax_node)
  end select
  
end subroutine GeomechRealizLocalToLocalWithArray

! ************************************************************************** !

subroutine GeomechRealizPrintCouplers(geomech_realization)
  ! 
  ! Print boundary data for geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Coupler_module
  
  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  
  type(geomech_patch_type), pointer :: patch
  type(geomech_coupler_type), pointer :: cur_coupler
  type(option_type), pointer :: option
 
  option => geomech_realization%option
 
  if (.not.OptionPrintToFile(option)) return
  
  patch => geomech_realization%geomech_patch
   
  cur_coupler => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(cur_coupler)) exit
    call GeomechRealizPrintCoupler(cur_coupler,option)    
    cur_coupler => cur_coupler%next
  enddo
     
  cur_coupler => patch%geomech_source_sink_list%first
  do
    if (.not.associated(cur_coupler)) exit
    call GeomechRealizPrintCoupler(cur_coupler,option)    
    cur_coupler => cur_coupler%next
  enddo
 
end subroutine GeomechRealizPrintCouplers

! ************************************************************************** !

subroutine GeomechRealizPrintCoupler(coupler,option)
  ! 
  ! Prints boundary condition coupler for geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Coupler_module
  
  implicit none
  
  type(geomech_coupler_type) :: coupler
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(geomech_condition_type), pointer :: geomech_condition
  type(gm_region_type), pointer :: region
   
98 format(40('=+'))
99 format(80('-'))
  
  geomech_condition => coupler%geomech_condition
  region => coupler%region

  write(option%fid_out,*)
  write(option%fid_out,98)


  select case(coupler%itype)
    case(GM_BOUNDARY_COUPLER_TYPE)
      string = 'Geomech Boundary Condition'
    case(GM_SRC_SINK_COUPLER_TYPE)
      string = 'Geomech Source Sink'
  end select
  write(option%fid_out,'(/,2x,a,/)') trim(string)

  write(option%fid_out,99)
101 format(5x,'     Geomech Condition: ',2x,a)
  if (associated(geomech_condition)) &
    write(option%fid_out,101) trim(geomech_condition%name)
102 format(5x,'             Region: ',2x,a)
  if (associated(region)) &
    write(option%fid_out,102) trim(region%name)
  write(option%fid_out,99)
  
  if (associated(geomech_condition)) then
    call GeomechConditionPrint(geomech_condition,option)
  endif
 
end subroutine GeomechRealizPrintCoupler

! ************************************************************************** !

subroutine GeomechRealizPassFieldPtrToPatch(geomech_realization)
  ! 
  ! This subroutine passes field to patch
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  use Option_module

  implicit none
  
  class(realization_geomech_type) :: geomech_realization

  type(geomech_patch_type), pointer :: patch

  patch => geomech_realization%geomech_patch
   
  patch%geomech_field => geomech_realization%geomech_field
  
end subroutine GeomechRealizPassFieldPtrToPatch

! ************************************************************************** !

subroutine GeomechRealizProcessGeomechCouplers(geomech_realization)
  ! 
  ! This subroutine sets up couplers in
  ! geomech realization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  implicit none
  
  class(realization_geomech_type) :: geomech_realization

  type(geomech_patch_type), pointer :: patch
  
  patch => geomech_realization%geomech_patch
  
  call GeomechPatchProcessGeomechCouplers(patch, &
                                   geomech_realization%geomech_conditions, &
                                   geomech_realization%option)
 
end subroutine GeomechRealizProcessGeomechCouplers

! ************************************************************************** !

subroutine GeomechRealizProcessGeomechConditions(geomech_realization)
  ! 
  ! This subroutine sets up condition in
  ! geomech realization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Dataset_Base_class
  use Dataset_module

  implicit none

  class(realization_geomech_type), pointer :: geomech_realization
  
  type(geomech_condition_type), pointer :: cur_geomech_condition
  type(geomech_sub_condition_type), pointer :: cur_geomech_sub_condition
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: dataset_name
  class(dataset_base_type), pointer :: dataset
  PetscInt :: i
  
  option => geomech_realization%option
  
  ! loop over geomech conditions looking for linkage to datasets
  cur_geomech_condition => geomech_realization%geomech_conditions%first
  do
    if (.not.associated(cur_geomech_condition)) exit
      do i = 1, size(cur_geomech_condition%sub_condition_ptr)
        ! find dataset
        call DatasetFindInList(geomech_realization%geomech_datasets, &
                 cur_geomech_condition%sub_condition_ptr(i)%ptr%dataset, &
                 cur_geomech_condition%default_time_storage, &
                 string,option)
      enddo
     cur_geomech_condition => cur_geomech_condition%next
  enddo
  
end subroutine GeomechRealizProcessGeomechConditions

! ************************************************************************** !

subroutine GeomechRealizGetDataset(geomech_realization,vec,ivar,isubvar, &
                                   isubvar1)
  ! 
  ! This routine extracts variables indexed by
  ! ivar and isubvar from geomechanics realization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/03/13
  ! 

  implicit none

  class(realization_geomech_type) :: geomech_realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1

  call GeomechPatchGetDataset(geomech_realization%geomech_patch, &
                              geomech_realization%geomech_field, &
                              geomech_realization%option, &
                              geomech_realization%output_option, &
                              vec,ivar,isubvar,isubvar1)

end subroutine GeomechRealizGetDataset

! ************************************************************************** !

subroutine GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
  ! 
  ! This subroutine addes a geomechanics
  ! coupler to a geomechanics realization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  use Geomechanics_Coupler_module

  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  type(geomech_coupler_type), pointer :: coupler
  
  type(geomech_patch_type), pointer :: patch
  type(geomech_coupler_type), pointer :: new_coupler
  
  patch => geomech_realization%geomech_patch
  
 ! only add to geomech list for now, since they will be split out later
  new_coupler => GeomechCouplerCreate(coupler)
  select case(coupler%itype)
    case(GM_BOUNDARY_COUPLER_TYPE)
      call GeomechCouplerAddToList(new_coupler, &
                                   patch%geomech_boundary_condition_list)
    case(GM_SRC_SINK_COUPLER_TYPE)
      call GeomechCouplerAddToList(new_coupler,patch%geomech_source_sink_list)
  end select
  nullify(new_coupler)
  
  call GeomechCouplerDestroy(coupler)
 
end subroutine GeomechRealizAddGeomechCoupler

! ************************************************************************** !

subroutine GeomechRealizAddWaypointsToList(geomech_realization,waypoint_list)
  ! 
  ! Adds waypoints from BCs and source/sink
  ! to waypoint list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Option_module
  use Waypoint_module
  use Time_Storage_module

  implicit none
  
  class(realization_geomech_type) :: geomech_realization
  type(waypoint_list_type), pointer :: waypoint_list

  type(geomech_condition_type), pointer :: cur_geomech_condition
  type(geomech_sub_condition_type), pointer :: sub_condition
  type(waypoint_type), pointer :: waypoint, cur_waypoint
  type(option_type), pointer :: option
  PetscInt :: itime, isub_condition
  PetscReal :: temp_real, final_time
  PetscReal, pointer :: times(:)

  option => geomech_realization%option
  nullify(times)
  
  ! set flag for final output
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final) then
      cur_waypoint%print_snap_output = &
        geomech_realization%output_option%print_final_snap
      exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  ! use final time in conditional below
  if (associated(cur_waypoint)) then
    final_time = cur_waypoint%time
  else
    option%io_buffer = 'Final time not found in GeomechRealizAddWaypointsToList'
    call printErrMsg(option)
  endif

  ! add update of geomech conditions
  cur_geomech_condition => geomech_realization%geomech_conditions%first
  do
    if (.not.associated(cur_geomech_condition)) exit
    if (cur_geomech_condition%sync_time_with_update) then
      do isub_condition = 1, cur_geomech_condition%num_sub_conditions
        sub_condition => cur_geomech_condition% &
                         sub_condition_ptr(isub_condition)%ptr
        call TimeStorageGetTimes(sub_condition%dataset%time_storage, option, &
                                final_time, times)
        if (associated(times)) then
          if (size(times) > 1000) then
            option%io_buffer = 'For geomech condition "' // &
              trim(cur_geomech_condition%name) // &
              '" dataset "' // trim(sub_condition%name) // &
              '", the number of times is excessive for synchronization ' // &
              'with waypoints.'
            call printErrMsg(option)
          endif
          do itime = 1, size(times)
            waypoint => WaypointCreate()
            waypoint%time = times(itime)
            waypoint%update_conditions = PETSC_TRUE
            call WaypointInsertInList(waypoint,waypoint_list)
          enddo
          deallocate(times)
          nullify(times)
        endif
      enddo
    endif
    cur_geomech_condition => cur_geomech_condition%next
  enddo
      
end subroutine GeomechRealizAddWaypointsToList

! ************************************************************************** !

subroutine GeomechRealizDestroy(geomech_realization)
  ! 
  ! This subroutine deallocates geomechanics realization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  ! 

  implicit none
  
  class(realization_geomech_type), pointer :: geomech_realization
  
  if (.not.associated(geomech_realization)) return
    
  call GeomechFieldDestroy(geomech_realization%geomech_field)

  call OutputOptionDestroy(geomech_realization%output_option)

  call GeomechRegionDestroyList(geomech_realization%geomech_region_list)

  call GeomechConditionDestroyList(geomech_realization%geomech_conditions)
  
  if (associated(geomech_realization%geomech_debug)) &
    deallocate(geomech_realization%geomech_debug)
  nullify(geomech_realization%geomech_debug)
  
  if (associated(geomech_realization%geomech_material_property_array)) &
    deallocate(geomech_realization%geomech_material_property_array)
  nullify(geomech_realization%geomech_material_property_array)
  if (associated(geomech_realization%geomech_patch)) &
    call GeomechanicsPatchDestroy(geomech_realization%geomech_patch)                                       
  call GeomechanicsMaterialPropertyDestroy(geomech_realization% &
                                           geomech_material_properties)
  call GeomechDiscretizationDestroy(geomech_realization%geomech_discretization)

  if (associated(geomech_realization%output_option)) &
    deallocate(geomech_realization%output_option)
  nullify(geomech_realization%output_option)

  if (associated(geomech_realization)) deallocate(geomech_realization)
  nullify(geomech_realization)
  
end subroutine GeomechRealizDestroy


end module Geomechanics_Realization_class
