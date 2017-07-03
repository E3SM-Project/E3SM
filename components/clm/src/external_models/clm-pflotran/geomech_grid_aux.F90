module Geomechanics_Grid_Aux_module

  use Grid_Unstructured_Cell_module
  use Gauss_module
  use PFLOTRAN_Constants_module
  use Geometry_module  

  implicit none

  private 
  
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#if defined(SCORPIO)
  include "scorpiof.h"
#endif
 
  type, public :: geomech_grid_type
    ! variables for geomechanics grid (unstructured) for finite element formulation
    ! The dofs (displacements, here) are solved at the nodes 
    ! Notation corresponding to finite volume: node implies cell vertex and element means cell
    PetscInt :: global_offset                ! offset in petsc ordering for the first vertex on a processor
    PetscInt :: global_offset_elem           ! offset in petsc ordering for the first element on a processor
    PetscInt :: nmax_elem                    ! Total number of elements in the global domain
    PetscInt :: nlmax_elem                   ! Total number of non-ghosted elements on a processor
    PetscInt :: nmax_node                    ! Total number of nodes in the global domain 
    PetscInt :: nlmax_node                   ! Total number of non-ghosted nodes on a processor
    PetscInt :: ngmax_node                   ! Total number of ghosted nodes on a processor
    PetscInt :: num_ghost_nodes              ! Number of ghost nodes on a processor
    PetscInt, pointer :: elem_ids_natural(:) ! Natural numbering of elements on a processor
    PetscInt, pointer :: elem_ids_petsc(:)   ! Petsc numbering of elements on a processor
    AO :: ao_natural_to_petsc                ! mapping of natural to Petsc ordering
    AO :: ao_natural_to_petsc_nodes          ! mapping of natural to Petsc ordering of vertices
    PetscInt :: max_ndual_per_elem           ! Max. number of dual elements connected to an element
    PetscInt :: max_nnode_per_elem           ! Max. number of nodes per element
    PetscInt :: max_elem_sharing_a_node      ! Max. number of elements sharing a common node
    PetscInt, pointer :: elem_type(:)        ! Type of element
    PetscInt, pointer :: elem_nodes(:,:)     ! Node number on each element
    type(point3d_type), pointer :: nodes(:)    ! Coordinates of the nodes
    PetscInt, pointer :: node_ids_ghosted_petsc(:)   ! Petsc ids of ghosted local nodes
    PetscInt, pointer :: node_ids_ghosted_natural(:) ! Natural ids of ghosted local nodes
    PetscInt, pointer :: node_ids_local_natural(:)   ! Natural ids of local nodes
    PetscInt, pointer :: ghosted_node_ids_natural(:) ! Natural ids of the ghost nodes only
    PetscInt, pointer :: ghosted_node_ids_petsc(:)   ! Petsc ids of the ghost nodes only
    PetscInt, pointer :: nL2G(:),nG2L(:),nG2A(:)
    type(gauss_type), pointer :: gauss_node(:)
    character(len=MAXSTRINGLENGTH) :: mapping_filename ! mapping between subsurf and geomech meshes
    PetscInt :: mapping_num_cells
    PetscInt, pointer :: mapping_cell_ids_flow(:)
    PetscInt, pointer :: mapping_vertex_ids_geomech(:)
    Vec :: no_elems_sharing_node_loc
    Vec :: no_elems_sharing_node
  end type geomech_grid_type
  

  type, public :: gmdm_type                  ! Geomech. DM type
    ! local: included both local (non-ghosted) and ghosted nodes
    ! global: includes only local (non-ghosted) nodes
    PetscInt :: ndof
    ! for the below
    ! ghosted = local (non-ghosted) and ghosted nodes
    ! local = local (non-ghosted) nodes
    IS :: is_ghosted_local                   ! IS for ghosted nodes with local on-processor numbering
    IS :: is_local_local                     ! IS for local nodes with local on-processor numbering
    IS :: is_ghosted_petsc                   ! IS for ghosted nodes with petsc numbering
    IS :: is_local_petsc                     ! IS for local nodes with petsc numbering
    IS :: is_ghosts_local                    ! IS for ghost nodes with local on-processor numbering
    IS :: is_ghosts_petsc                    ! IS for ghost nodes with petsc numbering
    IS :: is_local_natural                   ! IS for local nodes with natural (global) numbering
    VecScatter :: scatter_ltog               ! scatter context for local to global updates
    VecScatter :: scatter_gtol               ! scatter context for global to local updates
    VecScatter :: scatter_ltol               ! scatter context for local to local updates
    VecScatter :: scatter_gton               ! scatter context for global to natural updates
    VecScatter :: scatter_gton_elem          ! scatter context for global to natural updates for elements(cells)
    ISLocalToGlobalMapping :: mapping_ltog   ! petsc vec local to global mapping
!geh: deprecated in PETSc in spring 2014 
!    ISLocalToGlobalMapping :: mapping_ltogb  ! block form of mapping_ltog
    ISLocalToGlobalMapping :: mapping_ltog_elem   ! For elements
!geh: deprecated in PETSc in spring 2014    
!    ISLocalToGlobalMapping :: mapping_ltogb_elem  ! For elements
    Vec :: global_vec                        ! global vec (no ghost nodes), petsc-ordering
    Vec :: local_vec                         ! local vec (includes local and ghosted nodes), local ordering
    Vec :: global_vec_elem
    VecScatter :: scatter_subsurf_to_geomech_ndof ! scatter context between subsurface and
                                                  ! geomech grids for 1-DOF
    VecScatter :: scatter_geomech_to_subsurf_ndof ! scatter context between geomech and
                                                  ! subsurface grids for N-DOFs (Ndof = ngeomechdof)

  end type gmdm_type

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: GMGridCreate, &
            GMGridDestroy, &
            GMDMCreate, &
            GMDMDestroy, &
            GMCreateGMDM, &
            GMGridDMCreateVector, &
            GMGridDMCreateVectorElem, &
            GMGridDMCreateJacobian, &
            GMGridMapIndices
            
  
contains

! ************************************************************************** !
!
! GMDMCreate: Creates a geomech grid distributed mesh object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
function GMDMCreate()

  implicit none
  
  type(gmdm_type), pointer :: GMDMCreate
  type(gmdm_type), pointer :: gmdm

  allocate(gmdm)
  gmdm%is_ghosted_local = 0
  gmdm%is_local_local = 0
  gmdm%is_ghosted_petsc = 0
  gmdm%is_local_petsc = 0
  gmdm%is_ghosts_local = 0
  gmdm%is_ghosts_petsc = 0
  gmdm%is_local_natural = 0
  gmdm%scatter_ltog = 0
  gmdm%scatter_gtol = 0
  gmdm%scatter_ltol  = 0
  gmdm%scatter_gton = 0
  gmdm%scatter_gton_elem = 0
  gmdm%mapping_ltog = 0
  gmdm%mapping_ltog_elem = 0
  gmdm%global_vec = 0
  gmdm%local_vec = 0
  gmdm%global_vec_elem = 0
  gmdm%scatter_subsurf_to_geomech_ndof = 0
  gmdm%scatter_geomech_to_subsurf_ndof = 0

  GMDMCreate => gmdm

end function GMDMCreate

! ************************************************************************** !
!
! GMGridCreate: Creates a geomechanics grid object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
function GMGridCreate()

  implicit none
  
  type(geomech_grid_type), pointer :: GMGridCreate
  type(geomech_grid_type), pointer :: geomech_grid

  allocate(geomech_grid)

  ! variables for all unstructured grids
  
  geomech_grid%global_offset = 0
  geomech_grid%global_offset_elem = 0
  geomech_grid%nmax_elem = 0
  geomech_grid%nlmax_elem = 0
  geomech_grid%nmax_node = 0
  geomech_grid%nlmax_node = 0
  geomech_grid%num_ghost_nodes = 0
  nullify(geomech_grid%elem_ids_natural)
  nullify(geomech_grid%elem_ids_petsc)
  geomech_grid%ao_natural_to_petsc = 0
  geomech_grid%ao_natural_to_petsc_nodes = 0
  geomech_grid%max_ndual_per_elem = 0
  geomech_grid%max_nnode_per_elem = 0
  geomech_grid%max_elem_sharing_a_node = 0
  nullify(geomech_grid%elem_type)
  nullify(geomech_grid%elem_nodes)
  nullify(geomech_grid%nodes)
  nullify(geomech_grid%node_ids_ghosted_petsc)
  nullify(geomech_grid%node_ids_ghosted_natural)
  nullify(geomech_grid%ghosted_node_ids_natural)
  nullify(geomech_grid%ghosted_node_ids_petsc)
  nullify(geomech_grid%gauss_node)
  geomech_grid%no_elems_sharing_node_loc = 0
  geomech_grid%no_elems_sharing_node = 0

  GMGridCreate => geomech_grid
  
end function GMGridCreate

! ************************************************************************** !
!
! GMCreateGMDM: Mapping/scatter contexts are created for PETSc DM object
! author: Satish Karra, LANL
! date: 05/30/13
!
! ************************************************************************** !
subroutine GMCreateGMDM(geomech_grid,gmdm,ndof,option)

  use Option_module
  use Utility_module, only: reallocateIntArray
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscdm.h"  
#include "petsc/finclude/petscdm.h90"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscviewer.h"

  type(geomech_grid_type) :: geomech_grid
  type(option_type) :: option
  type(gmdm_type), pointer :: gmdm
  PetscInt :: ndof
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof
  IS :: is_tmp, is_tmp_petsc, is_tmp_natural
  Vec :: vec_tmp
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: ndof_word
  character(len=MAXSTRINGLENGTH) :: string
  PetscViewer :: viewer
  PetscInt, allocatable :: int_array(:), int_array2(:)
  
  gmdm => GMDMCreate()
  gmdm%ndof = ndof
  
#if GEOMECH_DEBUG
  write(ndof_word,*) ndof
  ndof_word = adjustl(ndof_word)
  ndof_word = '_' // trim(ndof_word)
  string = 'Vectors_nodes' // ndof_word
  call printMsg(option,string)
#endif

  ! create global vec
  call VecCreate(option%mycomm,gmdm%global_vec,ierr);CHKERRQ(ierr)
  call VecSetSizes(gmdm%global_vec,geomech_grid%nlmax_node*ndof,PETSC_DECIDE,&
                    ierr);CHKERRQ(ierr)
  call VecSetBlockSize(gmdm%global_vec,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(gmdm%global_vec,ierr);CHKERRQ(ierr)

  ! create local vec
  call VecCreate(PETSC_COMM_SELF,gmdm%local_vec,ierr);CHKERRQ(ierr)
  call VecSetSizes(gmdm%local_vec,geomech_grid%ngmax_node*ndof,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(gmdm%local_vec,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(gmdm%local_vec,ierr);CHKERRQ(ierr)
  
   
  ! IS for global numbering of local, non-ghosted vertices
  ! ISCreateBlock requires block ids, not indices.  Therefore, istart should be
  ! the offset of the block from the beginning of the vector.

#if GEOMECH_DEBUG
  string = 'Index_sets_nodes' // ndof_word
  call printMsg(option,string)
#endif  
  
  ! SK, Note: All the numbering are to be 0-based before creating IS
  
  allocate(int_array(geomech_grid%nlmax_node))
  do local_id = 1, geomech_grid%nlmax_node
    int_array(local_id) = (local_id-1) + geomech_grid%global_offset
  enddo

  ! arguments for ISCreateBlock():
  ! option%mycomm  - the MPI communicator
  ! ndof  - number of elements in each block
  ! geomech_grid%nlmax  - the length of the index set
  !                                      (the number of blocks
  ! int_array  - the list of integers, one for each block and count
  !              of block not indices
  ! PETSC_COPY_VALUES  - see PetscCopyMode, only PETSC_COPY_VALUES and
  !                      PETSC_OWN_POINTER are supported in this routine
  ! gmdm%is_local_petsc - the new index set
  ! ierr - PETScErrorCode
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%nlmax_node, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_local_petsc, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  string = 'geomech_is_local_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(gmdm%is_local_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  ! IS for local numbering of ghost nodes
  if (geomech_grid%num_ghost_nodes > 0) then
    allocate(int_array(geomech_grid%num_ghost_nodes))
    do ghosted_id = 1, geomech_grid%num_ghost_nodes
      int_array(ghosted_id) = (ghosted_id+geomech_grid%nlmax_node-1)
    enddo
  else
    allocate(int_array(1))
    int_array(1) = 0
  endif
  
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%num_ghost_nodes, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_ghosts_local, &
                     ierr);CHKERRQ(ierr)
                     
  if (allocated(int_array)) deallocate(int_array) 
                                       
#if GEOMECH_DEBUG
  string = 'geomech_is_ghosts_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(gmdm%is_ghosts_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif   

  ! IS for petsc numbering of ghost nodes
  if (geomech_grid%num_ghost_nodes > 0) then
    allocate(int_array(geomech_grid%num_ghost_nodes))
    do ghosted_id = 1, geomech_grid%num_ghost_nodes
      int_array(ghosted_id) = &
        (geomech_grid%ghosted_node_ids_petsc(ghosted_id)-1)
    enddo
  else
    allocate(int_array(1))
    int_array(1) = 0
  endif
  
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%num_ghost_nodes, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_ghosts_petsc, &
                     ierr);CHKERRQ(ierr)
                     
  if (allocated(int_array)) deallocate(int_array) 
                                       
#if GEOMECH_DEBUG
  string = 'geomech_is_ghosts_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(gmdm%is_ghosts_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif 

  ! IS for local numbering of local, non-ghosted cells
  allocate(int_array(geomech_grid%nlmax_node))
  do local_id = 1, geomech_grid%nlmax_node
    int_array(local_id) = (local_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%nlmax_node, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_local_local, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  string = 'geomech_is_local_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(gmdm%is_local_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  ! IS for ghosted numbering of local ghosted cells
  allocate(int_array(geomech_grid%ngmax_node))
  do ghosted_id = 1, geomech_grid%ngmax_node
    int_array(ghosted_id) = (ghosted_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%ngmax_node, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_ghosted_local, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  string = 'geomech_is_ghosted_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(gmdm%is_ghosted_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
                  
  ! IS for petsc numbering of local ghosted cells
  allocate(int_array(geomech_grid%ngmax_node))
  do local_id = 1, geomech_grid%nlmax_node
    int_array(local_id) = (local_id-1) + geomech_grid%global_offset
  enddo
  if (geomech_grid%num_ghost_nodes > 0) then
    do ghosted_id = 1,geomech_grid%num_ghost_nodes
      int_array(geomech_grid%nlmax_node+ghosted_id) = &
        (geomech_grid%ghosted_node_ids_petsc(ghosted_id)-1)
    enddo
  endif
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%ngmax_node, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_ghosted_petsc, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  string = 'geomech_is_ghosted_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(gmdm%is_ghosted_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif    

 ! create a local to global mapping
#if GEOMECH_DEBUG
  string = 'geomech_ISLocalToGlobalMapping' // ndof_word
  call printMsg(option,string)
#endif

  call ISLocalToGlobalMappingCreateIS(gmdm%is_ghosted_petsc, &
                                      gmdm%mapping_ltog,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  string = 'geomech_mapping_ltog' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISLocalToGlobalMappingView(gmdm%mapping_ltog,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
               
#if GEOMECH_DEBUG
  string = 'geomech_local to global' // ndof_word
  call printMsg(option,string)
#endif

  ! Create local to global scatter
  ! Note this will also use the data from ghosted nodes on local rank into
  ! a global vec
  call VecScatterCreate(gmdm%local_vec,gmdm%is_ghosted_local,gmdm%global_vec, &
                        gmdm%is_ghosted_petsc,gmdm%scatter_ltog, &
                        ierr);CHKERRQ(ierr)
                        
#if GEOMECH_DEBUG
  string = 'geomech_scatter_ltog' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(gmdm%scatter_ltog,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

#if GEOMECH_DEBUG
  string = 'geomech_global to local' // ndof_word
  call printMsg(option,string)
#endif

  ! Create global to local scatter
  call VecScatterCreate(gmdm%global_vec,gmdm%is_ghosted_petsc,gmdm%local_vec, &
                        gmdm%is_ghosted_local,gmdm%scatter_gtol, &
                        ierr);CHKERRQ(ierr)
                        
#if GEOMECH_DEBUG
  string = 'geomech_scatter_gtol' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(gmdm%scatter_gtol,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

#if GEOMECH_DEBUG
  string = 'geomech_local to local' // ndof_word
  call printMsg(option,string)
#endif
  
  ! Create local to local scatter.  Essentially remap the global to local as
  ! PETSc does in daltol.c
  call VecScatterCopy(gmdm%scatter_gtol,gmdm%scatter_ltol,ierr);CHKERRQ(ierr)
  call ISGetIndicesF90(gmdm%is_local_local,int_ptr,ierr);CHKERRQ(ierr)
  call VecScatterRemap(gmdm%scatter_ltol,int_ptr,PETSC_NULL_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call ISRestoreIndicesF90(gmdm%is_local_local,int_ptr,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  string = 'geomech_scatter_ltol' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(gmdm%scatter_ltol,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  allocate(int_array(geomech_grid%nlmax_node))
  do local_id = 1, geomech_grid%nlmax_node 
    int_array(local_id) = (local_id-1) + geomech_grid%global_offset
  enddo
  call ISCreateGeneral(option%mycomm,geomech_grid%nlmax_node, &
                       int_array,PETSC_COPY_VALUES,is_tmp,ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
  call AOPetscToApplicationIS(geomech_grid%ao_natural_to_petsc_nodes, &
                              is_tmp,ierr);CHKERRQ(ierr)

  allocate(int_array(geomech_grid%nlmax_node))
  call ISGetIndicesF90(is_tmp,int_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1, geomech_grid%nlmax_node
    int_array(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(is_tmp,int_ptr,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp,ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%nlmax_node, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_local_natural, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)

#if GEOMECH_DEBUG
  string = 'geomech_is_local_natural' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(gmdm%is_local_natural,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecCreate(option%mycomm,vec_tmp,ierr);CHKERRQ(ierr)
  call VecSetSizes(vec_tmp,geomech_grid%nlmax_node*ndof,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vec_tmp,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vec_tmp,ierr);CHKERRQ(ierr)
  call VecScatterCreate(gmdm%global_vec,gmdm%is_local_petsc,vec_tmp, &
                        gmdm%is_local_natural,gmdm%scatter_gton, &
                        ierr);CHKERRQ(ierr)
  call VecDestroy(vec_tmp,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  string = 'geomech_scatter_gton' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(gmdm%scatter_gton,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

#endif

  ! Now for elements. Need this for writing tecplot output

  allocate(int_array(geomech_grid%nlmax_elem))
  do local_id = 1, geomech_grid%nlmax_elem
    int_array(local_id) = (local_id-1) + geomech_grid%global_offset_elem
  enddo
  
  call ISCreateGeneral(option%mycomm,geomech_grid%nlmax_elem, &
                       int_array,PETSC_COPY_VALUES,is_tmp_petsc, &
                       ierr);CHKERRQ(ierr)
  
  call AOPetscToApplicationIS(geomech_grid%ao_natural_to_petsc, &
                              is_tmp_petsc,ierr);CHKERRQ(ierr)
  allocate(int_array2(geomech_grid%nlmax_elem))
  call ISGetIndicesF90(is_tmp_petsc,int_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1, geomech_grid%nlmax_elem
    int_array2(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(is_tmp_petsc,int_ptr,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp_petsc,ierr);CHKERRQ(ierr)
  
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%nlmax_elem, &
                     int_array,PETSC_COPY_VALUES,is_tmp_petsc, &
                     ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%nlmax_elem, &
                     int_array2,PETSC_COPY_VALUES,is_tmp_natural, &
                     ierr);CHKERRQ(ierr)

  ! create global vec
  call VecCreate(option%mycomm,gmdm%global_vec_elem,ierr);CHKERRQ(ierr)
  call VecSetSizes(gmdm%global_vec_elem,geomech_grid%nlmax_elem*ndof, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(gmdm%global_vec_elem,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(gmdm%global_vec_elem,ierr);CHKERRQ(ierr)

  call VecCreate(option%mycomm,vec_tmp,ierr);CHKERRQ(ierr)
  call VecSetSizes(vec_tmp,geomech_grid%nlmax_elem*ndof,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vec_tmp,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vec_tmp,ierr);CHKERRQ(ierr)
  call VecScatterCreate(gmdm%global_vec_elem,is_tmp_petsc,vec_tmp, &
                        is_tmp_natural,gmdm%scatter_gton_elem, &
                        ierr);CHKERRQ(ierr)
  call VecDestroy(vec_tmp,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp_natural,ierr);CHKERRQ(ierr)
  
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%nlmax_elem, &
                     int_array,PETSC_COPY_VALUES,is_tmp_petsc, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array) 
  deallocate(int_array2)

  ! create a local to global mapping
  call ISLocalToGlobalMappingCreateIS(is_tmp_petsc, &
                                      gmdm%mapping_ltog_elem, &
                                      ierr);CHKERRQ(ierr)
!  call ISDestroy(is_tmp_petsc,ierr);CHKERRQ(ierr)
              
end subroutine GMCreateGMDM

! ************************************************************************** !
!
! GMGridDMCreateJacobian: Creates a Jacobian matrix
! author: Satish Karra, LANL
! date: 06/05/13
!
! ************************************************************************** !
subroutine GMGridDMCreateJacobian(geomech_grid,gmdm,mat_type,J,option)

  use Option_module
  
  implicit none
  
  type(geomech_grid_type) :: geomech_grid
  type(gmdm_type) :: gmdm
  type(option_type) :: option 
  MatType :: mat_type
  Mat :: J
  IS :: is_tmp

  PetscInt, allocatable :: d_nnz(:), o_nnz(:)
  PetscInt :: local_id, ineighbor, neighbor_id
  PetscInt :: ndof_local
  PetscInt :: ielem, node_count
  PetscInt :: ivertex1, ivertex2
  PetscInt :: local_id1, local_id2
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  
  allocate(d_nnz(geomech_grid%nlmax_node))
  allocate(o_nnz(geomech_grid%nlmax_node))

  ! The following is an approximate estimate only.
  ! Need to come up with a more accurate way of calculating d_nnz, o_nnz
  ! Assuming max vertices based on hex. 
  d_nnz = min(27,geomech_grid%nlmax_node)
  o_nnz = min(27,geomech_grid%nmax_node - geomech_grid%nlmax_node)
   
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_d_nnz_jacobian' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do local_id1 = 1, geomech_grid%nlmax_node
    write(86,'(i5)') d_nnz(local_id1)
  enddo  
  close(86)
  write(string,*) option%myrank
  string = 'geomech_o_nnz_jacobian' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do local_id1 = 1, geomech_grid%nlmax_node
    write(86,'(i5)') o_nnz(local_id1)
  enddo  
  close(86) 
#endif   
    
  ndof_local = geomech_grid%nlmax_node*gmdm%ndof
  select case(mat_type)
    case(MATAIJ)
      d_nnz = d_nnz*gmdm%ndof
      o_nnz = o_nnz*gmdm%ndof
    case(MATBAIJ)
    case default
      option%io_buffer = 'MatType not recognized in GMGridDMCreateJacobian'
      call printErrMsg(option)
  end select 
  
  call MatCreate(option%mycomm,J,ierr);CHKERRQ(ierr)
  call MatSetType(J,mat_type,ierr);CHKERRQ(ierr)
  call MatSetSizes(J,ndof_local,ndof_local,PETSC_DETERMINE,PETSC_DETERMINE, &
                   ierr);CHKERRQ(ierr)
  call MatSetFromOptions(J,ierr);CHKERRQ(ierr)  
  call MatXAIJSetPreallocation(J,gmdm%ndof,d_nnz,o_nnz, &
                               PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                               ierr); CHKERRQ(ierr)
  call MatSetLocalToGlobalMapping(J,gmdm%mapping_ltog, &
                                  gmdm%mapping_ltog,ierr);CHKERRQ(ierr)
  
                        
  deallocate(d_nnz)
  deallocate(o_nnz)
  
end subroutine GMGridDMCreateJacobian

! ************************************************************************** !
!
! GMGridDMCreateVector: Creates a global vector with PETSc ordering
! author: Satish Karra, LANL
! date: 06/04/13
!
! ************************************************************************** !
subroutine GMGridDMCreateVector(geomech_grid,gmdm,vec,vec_type,option)

  use Option_module

  implicit none
  
  type(geomech_grid_type) :: geomech_grid
  type(gmdm_type) :: gmdm
  type(option_type) :: option 
  Vec :: vec
  PetscInt :: vec_type
  PetscErrorCode :: ierr
  
  select case(vec_type)
    case(GLOBAL)
      call VecCreate(option%mycomm,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,geomech_grid%nlmax_node*gmdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetLocalToGlobalMapping(vec,gmdm%mapping_ltog, &
                                      ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,gmdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
    case(LOCAL)
      call VecCreate(PETSC_COMM_SELF,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,geomech_grid%ngmax_node*gmdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,gmdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
    case(NATURAL)
      call VecCreate(option%mycomm,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,geomech_grid%nlmax_node*gmdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,gmdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
  end select
    
end subroutine GMGridDMCreateVector

! ************************************************************************** !
!
! GMGridDMCreateVectorElem: Creates a global vector with PETSc ordering
!                           with size of elements and not vertices
! author: Satish Karra, LANL
! date: 06/04/13
!
! ************************************************************************** !
subroutine GMGridDMCreateVectorElem(geomech_grid,gmdm,vec,vec_type,option)

  use Option_module

  implicit none
  
  type(geomech_grid_type) :: geomech_grid
  type(gmdm_type) :: gmdm
  type(option_type) :: option 
  Vec :: vec
  PetscInt :: vec_type
  PetscErrorCode :: ierr
  
  select case(vec_type)
    case(GLOBAL)
      call VecCreate(option%mycomm,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,geomech_grid%nlmax_elem*gmdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetLocalToGlobalMapping(vec,gmdm%mapping_ltog_elem, &
                                      ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,gmdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
    case(LOCAL)
      call VecCreate(PETSC_COMM_SELF,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,geomech_grid%nlmax_elem*gmdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,gmdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
    case(NATURAL)
      call VecCreate(option%mycomm,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,geomech_grid%nlmax_elem*gmdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,gmdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
  end select
    
end subroutine GMGridDMCreateVectorElem

! ************************************************************************** !
!
! GMGridMapIndices: Maps global, local and natural indices of nodes to
! each other for geomech grid.
! author: Satish Karra, LANL
! date: 06/04/13
!
! ************************************************************************** !
subroutine GMGridMapIndices(geomech_grid,gmdm,nG2L,nL2G,nG2A,option)

  use Option_module

  implicit none
  
  type(geomech_grid_type) :: geomech_grid
  type(gmdm_type) :: gmdm
  type(option_type) :: option
  PetscInt, pointer :: nG2L(:)
  PetscInt, pointer :: nL2G(:)
  PetscInt, pointer :: nG2A(:)
  PetscInt, pointer :: int_ptr(:)
  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id

  ! The index mapping arrays are the following:
  ! nL2G :  not collective, local processor: local  =>  ghosted local  
  ! nG2L :  not collective, local processor:  ghosted local => local  
  ! nG2A :  not collective, ghosted local => natural

  allocate(nG2L(geomech_grid%ngmax_node))
  allocate(nL2G(geomech_grid%nlmax_node))
  allocate(nG2A(geomech_grid%ngmax_node))
  
  nG2L = 0

  do local_id = 1, geomech_grid%nlmax_node
    nL2G(local_id) = local_id
    nG2L(local_id) = local_id
  enddo

  call ISGetIndicesF90(gmdm%is_ghosted_petsc,int_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, geomech_grid%ngmax_node
    nG2A(ghosted_id) = int_ptr(ghosted_id)
  enddo
  call ISRestoreIndicesF90(gmdm%is_ghosted_petsc,int_ptr,ierr);CHKERRQ(ierr)
  call AOPetscToApplication(geomech_grid%ao_natural_to_petsc_nodes, &
                            geomech_grid%ngmax_node, &
                            nG2A,ierr);CHKERRQ(ierr)
  nG2A = nG2A + 1 ! 1-based

end subroutine GMGridMapIndices

! ************************************************************************** !
!
! GMDMDestroy: Deallocates a geomechanics grid distributed mesh object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
subroutine GMGridDestroy(geomech_grid)

  use Utility_module, only : DeallocateArray

  implicit none
  
  type(geomech_grid_type), pointer :: geomech_grid
  
  PetscInt :: i
  PetscErrorCode :: ierr
  
  if (.not.associated(geomech_grid)) return
  
  call DeallocateArray(geomech_grid%elem_ids_natural)
  call DeallocateArray(geomech_grid%elem_ids_petsc)
  call DeallocateArray(geomech_grid%elem_type)
  call DeallocateArray(geomech_grid%elem_nodes)
  call DeallocateArray(geomech_grid%node_ids_ghosted_petsc)
  call DeallocateArray(geomech_grid%node_ids_ghosted_natural)
  call DeallocateArray(geomech_grid%node_ids_local_natural)
  call DeallocateArray(geomech_grid%ghosted_node_ids_natural)
  call DeallocateArray(geomech_grid%ghosted_node_ids_petsc)
  if (geomech_grid%ao_natural_to_petsc /= 0) &
    call AODestroy(geomech_grid%ao_natural_to_petsc,ierr) 
  if (geomech_grid%ao_natural_to_petsc_nodes /= 0) then
    call AODestroy(geomech_grid%ao_natural_to_petsc_nodes,ierr);CHKERRQ(ierr)
  endif
  
  if (associated(geomech_grid%nodes)) &
    deallocate(geomech_grid%nodes)
  nullify(geomech_grid%nodes)
  
  if (associated(geomech_grid%gauss_node)) then
    do i = 1, size(geomech_grid%gauss_node)
      call GaussDestroy(geomech_grid%gauss_node(i))
    enddo
    deallocate(geomech_grid%gauss_node)
  endif
  
  nullify(geomech_grid%gauss_node)
 
  if (geomech_grid%no_elems_sharing_node_loc /= 0) then
    call VecDestroy(geomech_grid%no_elems_sharing_node_loc,ierr);CHKERRQ(ierr)
  endif
  if ( geomech_grid%no_elems_sharing_node /= 0) then
    call VecDestroy(geomech_grid%no_elems_sharing_node,ierr);CHKERRQ(ierr)
  endif
 
  deallocate(geomech_grid%nG2L)
  deallocate(geomech_grid%nL2G)
  deallocate(geomech_grid%nG2A)
  deallocate(geomech_grid%mapping_cell_ids_flow)
  deallocate(geomech_grid%mapping_vertex_ids_geomech)

  deallocate(geomech_grid)
  nullify(geomech_grid)
  
end subroutine GMGridDestroy

! ************************************************************************** !
!
! GMGridDestroy: Deallocates a geomechanics grid object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
subroutine GMDMDestroy(gmdm)

  implicit none
  
  type(gmdm_type), pointer :: gmdm
  
  PetscErrorCode :: ierr
  
  if (.not.associated(gmdm)) return

  call ISDestroy(gmdm%is_ghosted_local,ierr);CHKERRQ(ierr)
  call ISDestroy(gmdm%is_local_local,ierr);CHKERRQ(ierr)
  call ISDestroy(gmdm%is_ghosted_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(gmdm%is_local_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(gmdm%is_ghosts_local,ierr);CHKERRQ(ierr)
  call ISDestroy(gmdm%is_ghosts_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(gmdm%is_local_natural,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(gmdm%scatter_ltog,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(gmdm%scatter_gtol,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(gmdm%scatter_ltol,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(gmdm%scatter_gton,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(gmdm%scatter_gton_elem,ierr);CHKERRQ(ierr)
  call ISLocalToGlobalMappingDestroy(gmdm%mapping_ltog,ierr);CHKERRQ(ierr)
  call ISLocalToGlobalMappingDestroy(gmdm%mapping_ltog_elem, &
                                     ierr);CHKERRQ(ierr)
  call VecDestroy(gmdm%global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(gmdm%local_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(gmdm%global_vec_elem,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(gmdm%scatter_subsurf_to_geomech_ndof, &
                         ierr);CHKERRQ(ierr)
  call VecScatterDestroy(gmdm%scatter_geomech_to_subsurf_ndof, &
                         ierr);CHKERRQ(ierr)
  deallocate(gmdm)
  nullify(gmdm)

end subroutine GMDMDestroy

end module Geomechanics_Grid_Aux_module
