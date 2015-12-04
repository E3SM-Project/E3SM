module domainLateralMod

#ifdef USE_PETSC_LIB
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: domainMod
  !
  ! !DESCRIPTION:
  ! Module containing:
  ! - Information regarding lateral connectivity of the land grid,
  ! - PETSc-based framework to exchange data across processors.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use clm_varctl  , only : iulog
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils  , only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"    
  !

  type, public :: ugdm_type
     PetscInt               :: ndof                      ! number of degree of freedoms per grid cells

     IS                     :: is_ghosted_lorder         ! IS for ghosted cells in "local" numbering order
     IS                     :: is_ghosted_porder         ! IS for ghosted cells in "PETSc" numerbing order
     IS                     :: is_local_lorder           ! IS for local cells in "local" numbering order
     IS                     :: is_local_porder           ! IS for local cells in "PETSc" numbering order
     IS                     :: is_local_norder           ! IS for local cells in "natural" numbering order

     VecScatter             :: scatter_l2g               ! scatter context for local to global
     VecScatter             :: scatter_g2l               ! scatter context for global to local
     VecScatter             :: scatter_l2l               ! scatter context for local to local
     VecScatter             :: scatter_g2n               ! scatter context for global to natural
     VecScatter             :: scatter_n2g               ! scatter context for natural to global

     ISLocalToGlobalMapping :: mapping_ltog              ! petsc vec local to global mapping

     Vec                    :: global_vec                ! global vec (no ghost cells), petsc-ordering
     Vec                    :: local_vec                 ! local vec (includes local and ghosted cells), local ordering

     AO                     :: ao_norder_to_porder       !
  end type ugdm_type

  type, public :: ugrid_type

     PetscInt,pointer         :: gridsOnGrid_local (:,:) ! grid cell connectivity information in local-order

     PetscInt,pointer         :: grid_id_norder(:)       ! grid cell ids in natural-order
     PetscInt,pointer         :: grid_id_porder(:)       ! grid cell ids in PETSc-order

     PetscInt                 :: ngrid_local             ! number of local grid cells
     PetscInt                 :: ngrid_ghost             ! number of ghost grid cells
     PetscInt                 :: ngrid_ghosted           ! total number of grid cells

     PetscInt                 :: maxEdges                ! max number of neighbors

     AO                       :: ao_norder_to_porder     ! mapping between natural and PETSc order

  end type ugrid_type
     
  type, public :: domainlateral_type

     type(ugrid_type), pointer :: ugrid                  ! unstructured grid object
     
     type(ugdm_type), pointer :: dm_1dof                 ! PETSc DM for 1 DOF
     type(ugdm_type), pointer :: dm_nlevgrnddof          !

  end type domainlateral_type

  type(domainlateral_type)    , public :: ldomain_lateral
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public domainlateral_init          ! allocates/nans domain types
  !
  !EOP
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: domainlateral_init
  !
  ! !INTERFACE:
  subroutine domainlateral_init(domain_l, cellsOnCell_old, ncells_loc_old, maxEdges)
    !
    ! !ARGUMENTS:
    implicit none
    !
    !
    type(domainlateral_type) :: domain_l             ! domain datatype
    integer , intent(in)     :: cellsOnCell_old(:,:) !
    integer , intent(in)     :: ncells_loc_old       !
    integer , intent(in)     :: maxEdges             !

    allocate(domain_l%ugrid)
    allocate(domain_l%dm_1dof)
    allocate(domain_l%dm_nlevgrnddof)

    call create_ugrid(domain_l%ugrid, cellsOnCell_old, ncells_loc_old, maxEdges)

    call create_ugdm(domain_l%ugrid, domain_l%dm_1dof, 1)

  end subroutine domainlateral_init

  !------------------------------------------------------------------------------
  subroutine create_ugrid(ugrid, cellsOnCell_old, ncells_loc_old, maxEdges)
    !
    ! !DESCRIPTION:
    ! This subroutine populates the unstructured grid object.
    ! - This subroutine is called after domain decomposition.
    ! - ldecomp%gdc2glo(:) has the local grid cell ids in natural order
    !   after domain decomposition.
    !
    use decompMod, only : ldecomp, get_proc_bounds
    ! !ARGUMENTS:
    implicit none
    type(ugrid_type), pointer :: ugrid                 ! unstructured grid object
    integer , intent(in)     :: cellsOnCell_old(:,:)   ! grid cell level connectivity information as read in from netcdf file
    integer , intent(in)     :: maxEdges               ! maximum number of grid neighbor [size(cellsOnCell_old,1)]
    integer , intent(in)     :: ncells_loc_old         ! no. of grid cell for which grid connectivity information was read [size(cellsOnCells_old,2)]
    !
    ! !LOCAL VARIABLES:
    integer                  :: beg, end               ! beg/end of grid bounds at processsor level
    integer                  :: ngrid_local            ! number of local grid cells
    integer                  :: n,i,m                  ! indices
    integer, pointer         :: int_array(:)           ! temporary
    integer, pointer         :: int_array2(:)          ! temporary
    integer, pointer         :: int_array3(:)          ! temporary
    integer, pointer         :: int_array4(:)          ! temporary
    integer, pointer         :: int_array5(:)          ! temporary
    
    PetscInt                 :: offset                 ! temporary
    PetscInt                 :: cell_stride            ! temporary
    PetscInt                 :: idx                    ! temporary
    PetscInt                 :: icell, iedge           ! indices
    PetscInt                 :: count                  ! temporary

    PetscReal, pointer       :: r_array(:)             ! temporary
    PetscReal, pointer       :: real_ptr(:)            ! temporary
    PetscReal, pointer       :: real_ptr2(:)           ! temporary

    IS                       :: is_from, is_to         ! temporary
    IS                       :: is_scatter             ! temporary

    Vec                      :: ids_petsc, ids_petsc_2 ! ids in PETSc order
    Vec                      :: conn_info_old          ! store connectivity information before domain decomposition [ids in natural order]
    Vec                      :: conn_info_natural      ! store connectivity information after domain decompoistion [ids in natural order]
    Vec                      :: conn_info_petsc        ! store connectivity information after domain decompoistion [ids in PETSc order]

    VecScatter               :: vec_scatter            ! temporary

    PetscErrorCode           :: ierr                   ! get error code from PETSc


    call get_proc_bounds(beg, end)
    ngrid_local = end-beg+1

    !
    ! Step-1: Create IS natural-to-petsc
    !
    offset = 0
    call MPI_Exscan(ngrid_local, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)

    call VecCreateMPI(mpicom, ngrid_local   , PETSC_DETERMINE, ids_petsc  , ierr); CHKERRQ(ierr)
    call VecCreateMPI(mpicom, ncells_loc_old, PETSC_DETERMINE, ids_petsc_2, ierr); CHKERRQ(ierr)

    allocate(int_array(ngrid_local))
    allocate(r_array(ngrid_local))

    do n = beg,end
       int_array(n-beg+1) = ldecomp%gdc2glo(n) - 1
       r_array(n-beg+1) = offset + n - beg + 1
    enddo

    call VecSetValues(ids_petsc, ngrid_local, int_array, r_array, INSERT_VALUES, ierr);
    CHKERRQ(ierr);

    call VecAssemblyBegin(ids_petsc, ierr);
    call VecAssemblyEnd(ids_petsc, ierr);
    deallocate(r_array)
    deallocate(int_array)
        

    offset = 0
    call MPI_Exscan(ncells_loc_old, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)
    allocate(int_array(ncells_loc_old))
    do n = 1,ncells_loc_old
       int_array(n) = n + offset -1
    enddo
    call ISCreateGeneral(mpicom, ncells_loc_old, int_array, PETSC_COPY_VALUES, is_from, ierr);
    CHKERRQ(ierr)
    deallocate(int_array)

    offset = 0
    call MPI_Exscan(ncells_loc_old, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)

    allocate(int_array(ncells_loc_old))
    do n = 1,ncells_loc_old
       int_array(n) = n + offset -1
    enddo
    call ISCreateGeneral(mpicom, ncells_loc_old, int_array, PETSC_COPY_VALUES, is_to, ierr);
    CHKERRQ(ierr)
    deallocate(int_array)
    
    call VecScatterCreate(ids_petsc ,is_from, ids_petsc_2, is_to, vec_scatter, ierr);
    CHKERRQ(ierr)

    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to, ierr); CHKERRQ(ierr)

    call VecScatterBegin(vec_scatter,ids_petsc, ids_petsc_2, INSERT_VALUES, SCATTER_FORWARD, ierr);
    CHKERRQ(ierr)
    call VecScatterEnd(vec_scatter,ids_petsc, ids_petsc_2, INSERT_VALUES, SCATTER_FORWARD, ierr);
    CHKERRQ(ierr)
    call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)
    call VecDestroy(ids_petsc, ierr)

    allocate(int_array(ncells_loc_old))
    call VecGetArrayF90(ids_petsc_2, real_ptr, ierr)
    int_array(:) = INT(real_ptr(:) - 1)
    call VecRestoreArrayF90(ids_petsc_2, real_ptr, ierr)
    call VecDestroy(ids_petsc_2, ierr)

    cell_stride = 1 + maxEdges

    call ISCreateBlock(mpicom, cell_stride, ncells_loc_old, int_array, &
         PETSC_COPY_VALUES, is_scatter, ierr);CHKERRQ(ierr);
         
    deallocate(int_array)
    

    !
    ! Step-2: Build connectivity information BEFORE domain decomposition
    !         with cell-ids in natural order.

    ! -cell_nat_id
    ! dual_1
    ! dual_2
    ! dual_3
    !
    ! daul_N

    offset = 0
    call MPI_Exscan(ncells_loc_old, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)

    call VecCreate(mpicom, conn_info_old, ierr)
    call VecSetSizes(conn_info_old, ncells_loc_old*cell_stride, PETSC_DECIDE, ierr);CHKERRQ(ierr)
    call VecSetFromOptions(conn_info_old, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(conn_info_old, real_ptr, ierr)
    count = 0
    do icell = 1,ncells_loc_old
       count = count + 1
       real_ptr(count) = -(icell + offset)
       do iedge = 1, maxEdges
          count = count + 1
          real_ptr(count) = cellsOnCell_old(iedge,icell)
       enddo
    enddo
    call VecRestoreArrayF90(conn_info_old, real_ptr, ierr)
    
    !
    ! Step-3: Build connectivity information AFTER domain decomposition
    !         with cell-ids in natural order.
    !
    call VecCreate(mpicom,conn_info_natural,ierr);CHKERRQ(ierr)
    call VecSetSizes(conn_info_natural, &
         cell_stride*ngrid_local, &
         PETSC_DECIDE,ierr);CHKERRQ(ierr)
    call VecSetFromOptions(conn_info_natural,ierr);CHKERRQ(ierr)

    call VecScatterCreate(conn_info_old,PETSC_NULL_OBJECT,conn_info_natural,is_scatter, &
         vec_scatter,ierr);CHKERRQ(ierr)
    call ISDestroy(is_scatter,ierr);CHKERRQ(ierr)

    call VecScatterBegin(vec_scatter, conn_info_old, conn_info_natural, &
         INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr)
    call VecScatterEnd(vec_scatter, conn_info_old, conn_info_natural, &
         INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr)
    call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

    !
    ! Step-4: Build mapping between natural and PETSc order
    !
    offset = 0
    call MPI_Exscan(ngrid_local, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)

    allocate(int_array(ngrid_local))
    allocate(int_array2(ngrid_local))
    do n = beg, end
       int_array(n-beg+1)     = n-beg + offset
       int_array2(n-beg+1) = ldecomp%gdc2glo(n) - 1
    enddo

    call AOCreateBasic(mpicom, ngrid_local, int_array2, int_array, &
         ugrid%ao_norder_to_porder, ierr); CHKERRQ(ierr);
    deallocate(int_array)
    deallocate(int_array2)

    !
    ! Step-5: Build connectivity information AFTER domain decomposition
    !         with cell-ids in PETSc order
    call VecGetArrayF90(conn_info_natural, real_ptr, ierr); CHKERRQ(ierr)
    count = 0
    idx = 0
    do icell = 1, ngrid_local
       idx = idx + 1
       count = count + 1
       do iedge = 1, maxEdges
          idx = idx + 1
          if (real_ptr(idx) > 0) count = count + 1
       enddo
    enddo

    allocate(int_array(count))
    count = 0
    idx = 0
    do icell = 1, ngrid_local
       idx = idx + 1
       count = count + 1
       int_array(count) = -INT(real_ptr(idx))
       do iedge = 1, maxEdges
          idx = idx + 1
          if (real_ptr(idx) > 0) then
             count = count + 1
             int_array(count) = INT(real_ptr(idx))
          endif
       enddo
    enddo    
    call VecRestoreArrayF90(conn_info_natural, real_ptr, ierr); CHKERRQ(ierr)

    int_array = int_array - 1
    call AOApplicationToPetsc(ugrid%ao_norder_to_porder, count, int_array, ierr)
    CHKERRQ(ierr)
    int_array = int_array + 1

    call VecDuplicate(conn_info_natural, conn_info_petsc, ierr);CHKERRQ(ierr)
    call VecCopy(conn_info_natural, conn_info_petsc, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(conn_info_petsc, real_ptr, ierr); CHKERRQ(ierr)
    idx = 0
    count = 0
    do icell = 1, ngrid_local
       idx = idx + 1
       count = count + 1
       real_ptr(idx) = -int_array(count)
       do iedge = 1, maxEdges
          idx = idx + 1
          if (real_ptr(idx) > 0) then
             count = count + 1
             real_ptr(idx) = int_array(count)
          endif
       enddo
    enddo
    call VecRestoreArrayF90(conn_info_petsc, real_ptr, ierr); CHKERRQ(ierr)

    !
    ! Step-6: Find ghost cell ids in PETSc order
    !
    call VecGetArrayF90(conn_info_petsc, real_ptr, ierr)
    count = 0
    idx = 0
    do icell = 1, ngrid_local
       idx = idx + 1
       do iedge = 1, maxEdges
          idx = idx + 1
          if (real_ptr(idx) > 0) then
             if (real_ptr(idx) <= offset .or. &
                  real_ptr(idx) >  offset + ngrid_local ) then
                count = count + 1
             end if
          end if
       enddo
    enddo

    allocate(int_array(count))
    count = 0
    idx = 0
    do icell = 1, ngrid_local
       idx = idx + 1
       do iedge = 1, maxEdges
          idx = idx + 1
          if (real_ptr(idx) > 0) then
             if (real_ptr(idx) <= offset .or. &
                  real_ptr(idx) >  offset + ngrid_local ) then
                count = count + 1
                int_array(count) = real_ptr(idx)
             end if
          end if
       enddo
    enddo
    call VecRestoreArrayF90(conn_info_petsc, real_ptr, ierr); CHKERRQ(ierr)

    if (count > 0) then

       allocate(int_array2(count))
       allocate(int_array3(count))
       allocate(int_array4(count))
       allocate(int_array5(count))

       do n = 1, count
          int_array2(n) = n
       enddo

       int_array  = int_array - 1
       int_array2 = int_array2 - 1

       call PetscSortIntWithPermutation(count, int_array, int_array2, ierr)
       CHKERRQ(ierr);

       int_array  = int_array + 1
       int_array2 = int_array2 + 1

       int_array3 = 0
       int_array4 = 0
       int_array5 = 0

       ugrid%ngrid_ghost = 1
       int_array3(1)            = int_array(int_array2(1))

       do n = 1, count
          if (int_array3(ugrid%ngrid_ghost) < &
               int_array(int_array2(n))) then
             ugrid%ngrid_ghost = ugrid%ngrid_ghost + 1
             int_array3(ugrid%ngrid_ghost) = int_array(int_array2(n))
          end if
          int_array5(int_array2(n)) = n
          int_array4(n) = ugrid%ngrid_ghost
       enddo

    else

       ugrid%ngrid_ghost = 0

    endif

    !
    ! Step-7: Begin populating the ugrid object
    !
    ugrid%ngrid_local   = ngrid_local
    ugrid%ngrid_ghosted = ngrid_local + ugrid%ngrid_ghost

    allocate(ugrid%grid_id_porder(1:ugrid%ngrid_ghosted))
    allocate(ugrid%grid_id_norder(1:ugrid%ngrid_ghosted))

    ugrid%grid_id_porder = 0
    ugrid%grid_id_norder = 0
    
    allocate(ugrid%gridsOnGrid_local(maxEdges, ugrid%ngrid_local))

    ugrid%gridsOnGrid_local = 0

    call VecGetArrayF90(conn_info_petsc  , real_ptr, ierr); CHKERRQ(ierr)    
    call VecGetArrayF90(conn_info_natural, real_ptr2, ierr); CHKERRQ(ierr)    

    idx = 0
    count = 0
    do icell = 1, ngrid_local
       idx = idx + 1
       ugrid%grid_id_porder(icell) = -INT(real_ptr( idx))
       ugrid%grid_id_norder(icell) = -INT(real_ptr2(idx))

       do iedge = 1, maxEdges
          idx = idx + 1
          if (real_ptr(idx) > 0) then
             if (real_ptr(idx) <= offset .or. &
                  real_ptr(idx) >  offset + ngrid_local ) then
                count = count + 1
                ugrid%gridsOnGrid_local(iedge, icell) = int_array4(count) + ngrid_local
             else
                ugrid%gridsOnGrid_local(iedge, icell) = INT(real_ptr(idx)) - offset
             end if
          end if
       enddo
    enddo
    call VecRestoreArrayF90(conn_info_petsc  , real_ptr , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(conn_info_natural, real_ptr2, ierr); CHKERRQ(ierr)

    ! save information about ghost grid cells
    do n = 1, ugrid%ngrid_ghost
       ugrid%grid_id_porder(n + ngrid_local) = int_array3(n)
    enddo
    
    if (count > 0) then
       deallocate(int_array2)
       deallocate(int_array3)
       deallocate(int_array4)
       deallocate(int_array5)
    endif
    deallocate(int_array)

    !
    ! Step-8: Find the ghost cell ids in natural order
    !
    allocate(int_array(ugrid%ngrid_ghost))
    do icell = ngrid_local+1, ngrid_local+ugrid%ngrid_ghost
       int_array(icell - ngrid_local) = ugrid%grid_id_porder(icell)
    enddo

    int_array = int_array - 1 
    call AOPetscToApplication(ugrid%ao_norder_to_porder, ugrid%ngrid_ghost, int_array, ierr)
    CHKERRQ(ierr)
    int_array = int_array + 1
    do icell = ngrid_local+1, ngrid_local+ugrid%ngrid_ghost
       ugrid%grid_id_norder(icell) = int_array(icell - ngrid_local)
    enddo
    deallocate(int_array)
    
  end subroutine create_ugrid

  !------------------------------------------------------------------------------
  subroutine create_ugdm(ugrid, ugdm, ndof)
    !
    ! !DESCRIPTION:
    ! This subroutine creates PETSc DM to manage data for an unstructured
    ! grid.
    !
    use decompMod, only : ldecomp, get_proc_bounds
    ! !ARGUMENTS:
    implicit none
    type(ugrid_type), pointer :: ugrid               !
    type(ugdm_type), pointer  :: ugdm                !
    integer , intent(in)      :: ndof                !
    !
    !
    ! !LOCAL VARIABLES:
    integer, pointer          :: int_array(:)        ! temporary

    PetscInt                  :: offset              !
    PetscInt                  :: icell, iedge, count !
    PetscInt                  :: nsize
    PetscInt, pointer         :: int_ptr(:)          ! temporary
    Vec                       :: natural_vec
    PetscErrorCode            :: ierr                ! get error code from PETSc


    ugdm%ndof = ndof

    ! 1) Create is_local*

    allocate(int_array(ugrid%ngrid_local))


    ! 1.1) is_local_norder
    do icell = 1,ugrid%ngrid_local
       int_array(icell) = ugrid%grid_id_norder(icell)
    enddo
    int_array = int_array - 1
    
    call ISCreateBlock(mpicom, ndof, ugrid%ngrid_local, int_array, &
         PETSC_COPY_VALUES, ugdm%is_local_norder, ierr);CHKERRQ(ierr);

    ! 1.2) is_local_porder
    do icell = 1,ugrid%ngrid_local
       int_array(icell) = ugrid%grid_id_porder(icell)
    enddo
    int_array = int_array - 1
    
    call ISCreateBlock(mpicom, ndof, ugrid%ngrid_local, int_array, &
         PETSC_COPY_VALUES, ugdm%is_local_porder, ierr);CHKERRQ(ierr);

    ! 1.3) is_local_lorder
    do icell = 1,ugrid%ngrid_local
       int_array(icell) = icell
    enddo
    int_array = int_array - 1
    
    call ISCreateBlock(mpicom, ndof, ugrid%ngrid_local, int_array, &
         PETSC_COPY_VALUES, ugdm%is_local_lorder, ierr);CHKERRQ(ierr);

    deallocate(int_array)

    ! 2) Create is_ghosted*
    allocate(int_array(ugrid%ngrid_ghosted))

    ! 2.1) is_ghosted_porder
    do icell = 1, ugrid%ngrid_ghosted
       int_array(icell) = ugrid%grid_id_porder(icell)
    enddo
    int_array = int_array - 1
    
    call ISCreateBlock(mpicom, ndof, ugrid%ngrid_ghosted, int_array, &
         PETSC_COPY_VALUES, ugdm%is_ghosted_porder, ierr);CHKERRQ(ierr);

    ! 2.2) is_ghosted_lorder
    do icell = 1, ugrid%ngrid_ghosted
       int_array(icell) = icell
    enddo
    int_array = int_array - 1
    
    call ISCreateBlock(mpicom, ndof, ugrid%ngrid_ghosted, int_array, &
         PETSC_COPY_VALUES, ugdm%is_ghosted_lorder, ierr);CHKERRQ(ierr);

    deallocate(int_array)
    
    ! Create mapping local-to-global indexing:
    !   Mapping of indices of all cells (local + ghost) from local-order to
    !   PETSc-order
    call ISLocalToGlobalMappingCreateIS(ugdm%is_ghosted_porder, &
         ugdm%mapping_ltog, ierr); CHKERRQ(ierr)

    ! Create Vectors
    call VecCreate(mpicom, ugdm%global_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(ugdm%global_vec, ugrid%ngrid_local*ndof, PETSC_DECIDE,ierr); CHKERRQ(ierr)
    call VecSetBlockSize(ugdm%global_vec, ndof, ierr); CHKERRQ(ierr)
    call VecSetFromOptions(ugdm%global_vec,ierr); CHKERRQ(ierr)    
  
    call VecCreate(PETSC_COMM_SELF, ugdm%local_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(ugdm%local_vec, ugrid%ngrid_ghosted*ndof, PETSC_DECIDE, ierr); CHKERRQ(ierr)
    call VecSetBlockSize(ugdm%local_vec, ndof, ierr); CHKERRQ(ierr)
    call VecSetFromOptions(ugdm%local_vec, ierr); CHKERRQ(ierr)
    
    call VecCreate(mpicom, natural_vec,ierr); CHKERRQ(ierr)
    call VecSetSizes(natural_vec, ugrid%ngrid_local*ndof,PETSC_DECIDE, &
         ierr); CHKERRQ(ierr)
    call VecSetBlockSize(natural_vec,ndof,ierr); CHKERRQ(ierr)
    call VecSetFromOptions(natural_vec,ierr); CHKERRQ(ierr)
    
    ! Create VecScatter

    call VecScatterCreate(ugdm%local_vec, ugdm%is_local_lorder, ugdm%global_vec, &
         ugdm%is_local_porder, ugdm%scatter_l2g, ierr); CHKERRQ(ierr)
    
    call VecScatterCreate(ugdm%global_vec, ugdm%is_ghosted_porder,ugdm%local_vec, &
         ugdm%is_ghosted_lorder, ugdm%scatter_g2l, ierr); CHKERRQ(ierr)
    
    call VecScatterCopy(ugdm%scatter_g2l, ugdm%scatter_l2l,ierr); CHKERRQ(ierr)
    call ISGetIndicesF90(ugdm%is_local_lorder,int_ptr,ierr); CHKERRQ(ierr)
    call VecScatterRemap(ugdm%scatter_l2l,int_ptr,PETSC_NULL_INTEGER, &
         ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(ugdm%is_local_lorder,int_ptr,ierr); CHKERRQ(ierr)
    
    call VecScatterCreate(ugdm%global_vec,ugdm%is_local_porder, natural_vec, &
         ugdm%is_local_norder,ugdm%scatter_g2n, ierr); CHKERRQ(ierr)
    call VecScatterCreate(ugdm%global_vec,ugdm%is_local_norder, natural_vec, &
         ugdm%is_local_porder,ugdm%scatter_n2g, ierr); CHKERRQ(ierr)

  end subroutine create_ugdm


#else

  !-----------------------------------------------------------------------
  ! This is a stub for the case when PETSc is unavailable
  !
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use clm_varctl  , only : iulog
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils  , only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
     
  type, public :: domainlateral_type
     integer :: dummy
  end type domainlateral_type

  type(domainlateral_type)    , public :: ldomain_lateral
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public domainlateral_init          ! allocates/nans domain types
  !
  !EOP
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: domainlateral_init
  !
  ! !INTERFACE:
  subroutine domainlateral_init(domain_l, cellsOnCell_old, ncells_loc_old, maxEdges)
    !
    ! !ARGUMENTS:
    implicit none
    !
    !
    type(domainlateral_type) :: domain_l             ! domain datatype
    integer , intent(in)     :: cellsOnCell_old(:,:) !
    integer , intent(in)     :: ncells_loc_old       !
    integer , intent(in)     :: maxEdges             !
    character(len=255)       :: subname = 'domainlateral_init'


    call endrun(msg='ERROR ' // trim(subname) //': Requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

  end subroutine domainlateral_init

#endif

end module domainLateralMod
