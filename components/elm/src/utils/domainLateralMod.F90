module domainLateralMod

#include "shr_assert.h"

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
#include <petsc/finclude/petsc.h>
  !
  ! !USES:
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscdm
  use petscdmda
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use elm_varctl  , only : iulog
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils  , only : endrun
  use UnstructuredGridType, only : ugdm_type, ugrid_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: domainlateral_type

     type(ugrid_type), pointer :: ugrid                  ! unstructured grid object
     
     type(ugdm_type), pointer :: dm_1dof                 ! PETSc DM for 1 DOF
     type(ugdm_type), pointer :: dm_nlevgrnddof          !

  end type domainlateral_type

  type(domainlateral_type)    , public :: ldomain_lateral
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public domainlateral_init          ! allocates/nans domain types
  public ExchangeColumnLevelGhostData
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
  subroutine domainlateral_init(domain_l, cellsOnCell_old, edgesOnCell_old, &
       nEdgesOnCell_old, areaCell_old, dcEdge_old, dvEdge_old, &
       nCells_loc_old, nEdges_loc_old, maxEdges)
    !
    use decompMod, only : ldecomp, get_proc_bounds
    use UnstructuredGridType, only : create_ugrid, create_ugdm
    ! !ARGUMENTS:
    implicit none
    !
    !
    type(domainlateral_type) :: domain_l                ! domain datatype
    integer , intent(in)     :: cellsOnCell_old(:,:)    ! grid cell level connectivity information as read in from netcdf file
    integer , intent(in)     :: edgesOnCell_old(:,:)    ! index to determine distance between neighbors from dcEdge [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdgesOnCell_old(:)     ! number of edged                                            [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dvEdge_old(:)           ! distance between neighbors                                 [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dcEdge_old(:)           ! distance between vertices                                  [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: areaCell_old(:)         ! area of grid cell                                          [in natural order prior to domain decomposition]
    integer , intent(in)     :: nCells_loc_old          ! number of local cell-to-cell connections                   [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdges_loc_old          ! number of edges                                            [in natural order prior to domain decomposition]
    integer , intent(in)     :: maxEdges                ! max number of edges/neighbors
    !
    integer :: begg, endg

    allocate(domain_l%ugrid)
    allocate(domain_l%dm_1dof)
    allocate(domain_l%dm_nlevgrnddof)

    call get_proc_bounds(begg, endg)

    call create_ugrid(domain_l%ugrid, mpicom, begg, endg, ldecomp%gdc2glo, &
         cellsOnCell_old, ncells_loc_old, maxEdges)

    call create_ugdm(domain_l%ugrid, domain_l%dm_1dof, 1)

    call save_geometric_attributes(edgesOnCell_old, &
         nEdgesOnCell_old, areaCell_old, dcEdge_old, dvEdge_old, &
         nCells_loc_old, nEdges_loc_old, maxEdges)

  end subroutine domainlateral_init

  !------------------------------------------------------------------------------
  subroutine save_geometric_attributes(edgesOnCell_old, &
       nEdgesOnCell_old, areaCell_old, dcEdge_old, dvEdge_old, &
       nCells_loc_old, nEdges_loc_old, maxEdges)
    !
    ! !DESCRIPTION:
    ! Save following geometric attributes:
    !  - grid cell area,
    !  - centroidal distance between neighboring grid cells, and
    !  - edge length between neighboring grid cells.
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer , intent(in)     :: edgesOnCell_old(:,:) ! index to determine distance between neighbors from dcEdge [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdgesOnCell_old(:)  ! number of edges                                           [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dcEdge_old(:)        ! distance between neighbors                                [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dvEdge_old(:)        ! distance between vertices                                 [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: areaCell_old(:)      ! area of grid cell                                         [in natural order prior to domain decomposition]
    integer , intent(in)     :: nCells_loc_old       ! number of local cell-to-cell connections                  [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdges_loc_old       ! number of edges                                           [in natural order prior to domain decomposition]
    integer , intent(in)     :: maxEdges             ! max number of edges/neighbors
    !
    ! !LOCAL VARIABLES:
    PetscInt                 :: ii                   ! temporary
    PetscInt                 :: icell, iedge         ! indices
    PetscInt                 :: dcdv_count           ! counter for non-zero dc/dv
    PetscInt                 :: count                ! temporary
    PetscInt                 :: nblocks              ! temporary
    PetscInt, pointer        :: int_array(:)         ! temporary
    IS                       :: is_from, is_to       ! temporary
    VecScatter               :: scatter              ! temprary
    Vec                      :: dcdvEdge_glb_vec     ! global vectors for dcEdge
    Vec                      :: dcdvEdge_loc_vec     ! global vectors for dcEdge
    Vec                      :: attr_glb_vec         ! temporary global vector
    Vec                      :: attr_loc_vec         ! temporary sequential vector
    PetscReal, pointer       :: real_ptr(:)          ! temporary
    PetscReal, pointer       :: dcOnCell_old(:,:)    ! temporary array to hold distance between neighbors
    PetscReal, pointer       :: dvOnCell_old(:,:)    ! temporary array to hold distance between vertices
    PetscErrorCode           :: ierr                 ! get error code from PETSc

    allocate(dcOnCell_old(maxEdges, nCells_loc_old))
    allocate(dvOnCell_old(maxEdges, nCells_loc_old))

    dcOnCell_old = 0._r8
    dvOnCell_old = 0._r8

    dcdv_count = 0
    do icell = 1, nCells_loc_old
       dcdv_count = dcdv_count + nEdgesOnCell_old(icell)
    enddo

    nblocks = 2
    call VecCreate(mpicom, dcdvEdge_glb_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(dcdvEdge_glb_vec, nEdges_loc_old*nblocks, PETSC_DECIDE, ierr);
    CHKERRQ(ierr)
    call VecSetBlockSize(dcdvEdge_glb_vec, nblocks, ierr);CHKERRQ(ierr)
    call VecSetFromOptions(dcdvEdge_glb_vec, ierr);CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_SELF, dcdvEdge_loc_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(dcdvEdge_loc_vec, dcdv_count*nblocks, PETSC_DECIDE, ierr);
    CHKERRQ(ierr)
    call VecSetBlockSize(dcdvEdge_loc_vec, nblocks, ierr);CHKERRQ(ierr)
    call VecSetFromOptions(dcdvEdge_loc_vec, ierr);CHKERRQ(ierr)

    call VecGetArrayF90(dcdvEdge_glb_vec, real_ptr, ierr); CHKERRQ(ierr)
    count = 0;
    do iedge = 1, nEdges_loc_old
       count = count + 1
       real_ptr(count) = dcEdge_old(iedge)
       count = count + 1
       real_ptr(count) = dvEdge_old(iedge)
    enddo
    call VecRestoreArrayF90(dcdvEdge_glb_vec, real_ptr, ierr); CHKERRQ(ierr)

    ! Populate dcOnCell_old
    allocate(int_array(dcdv_count))
    do ii = 1,dcdv_count
       int_array(ii) = ii-1
    enddo
    call ISCreateBlock(mpicom, nblocks, dcdv_count, int_array, &
         PETSC_COPY_VALUES, is_to, ierr);CHKERRQ(ierr);

    dcdv_count = 0
    do icell = 1, nCells_loc_old
       do iedge = 1, nEdgesOnCell_old(icell)
          dcdv_count            = dcdv_count + 1
          int_array(dcdv_count) = edgesOnCell_old(iedge, icell) - 1
       enddo
    enddo
    call ISCreateBlock(mpicom, nblocks, dcdv_count, int_array, &
         PETSC_COPY_VALUES, is_from, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    call VecScatterCreate(dcdvEdge_glb_vec, is_from, dcdvEdge_loc_vec, is_to, &
         scatter, ierr); CHKERRQ(ierr)
    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to, ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, dcdvEdge_glb_vec, dcdvEdge_loc_vec, INSERT_VALUES, SCATTER_FORWARD, ierr);
    CHKERRQ(ierr);
    call VecScatterEnd(scatter, dcdvEdge_glb_vec, dcdvEdge_loc_vec, INSERT_VALUES, SCATTER_FORWARD, ierr);
    CHKERRQ(ierr);
    call VecScatterDestroy(scatter, ierr)

    call VecGetArrayF90(dcdvEdge_loc_vec, real_ptr, ierr); CHKERRQ(ierr)
    count = 0
    do icell = 1, nCells_loc_old
       do iedge = 1, nEdgesOnCell_old(icell)
          count = count + 1;
          dcOnCell_old(iedge, icell) = real_ptr(count)
          count = count + 1
          dvOnCell_old(iedge, icell) = real_ptr(count)
       enddo
    enddo
    call VecRestoreArrayF90(dcdvEdge_loc_vec, real_ptr, ierr); CHKERRQ(ierr)
    call VecDestroy(dcdvEdge_loc_vec, ierr); CHKERRQ(ierr)

    ! Aggregate data to be sent
    nblocks = maxEdges*2
    call VecCreate(mpicom, attr_glb_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(attr_glb_vec, nCells_loc_old*nblocks, PETSC_DECIDE, ierr);
    CHKERRQ(ierr)
    call VecSetBlockSize(attr_glb_vec, nblocks, ierr);CHKERRQ(ierr)
    call VecSetFromOptions(attr_glb_vec, ierr);CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_SELF, attr_loc_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(attr_loc_vec, ldomain_lateral%ugrid%ngrid_local*nblocks, PETSC_DECIDE, ierr);
    CHKERRQ(ierr)
    call VecSetBlockSize(attr_loc_vec, nblocks, ierr);CHKERRQ(ierr)
    call VecSetFromOptions(attr_loc_vec, ierr);CHKERRQ(ierr)

    call VecGetArrayF90(attr_glb_vec, real_ptr, ierr); CHKERRQ(ierr)
    count = 0
    do icell = 1, nCells_loc_old
       do iedge = 1, maxEdges
          count           = count + 1;
          real_ptr(count) = dcOnCell_old(iedge, icell)
       enddo

       do iedge = 1, maxEdges
          count           = count + 1;
          real_ptr(count) = dvOnCell_old(iedge, icell)
       enddo
    enddo
    call VecRestoreArrayF90(attr_glb_vec, real_ptr, ierr); CHKERRQ(ierr)
    deallocate(dcOnCell_old)
    deallocate(dvOnCell_old)

    allocate(int_array(ldomain_lateral%ugrid%ngrid_local))
    do ii = 1, ldomain_lateral%ugrid%ngrid_local
       int_array(ii) = ldomain_lateral%ugrid%grid_id_norder(ii) - 1
    enddo

    call ISCreateBlock(mpicom, nblocks, ldomain_lateral%ugrid%ngrid_local, int_array, &
         PETSC_COPY_VALUES, is_from, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    allocate(int_array(ldomain_lateral%ugrid%ngrid_local))
    do ii = 1, ldomain_lateral%ugrid%ngrid_local
       int_array(ii) = ii - 1
    enddo

    call ISCreateBlock(mpicom, nblocks, ldomain_lateral%ugrid%ngrid_local, int_array, &
         PETSC_COPY_VALUES, is_to, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    call VecScatterCreate(attr_glb_vec, is_from, attr_loc_vec, is_to, &
         scatter, ierr); CHKERRQ(ierr)
    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to, ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, attr_glb_vec, attr_loc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr);
    call VecScatterEnd(scatter, attr_glb_vec, attr_loc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr);
    call VecScatterDestroy(scatter, ierr)

    allocate(ldomain_lateral%ugrid%dcOnGrid_local(maxEdges, ldomain_lateral%ugrid%ngrid_local))
    allocate(ldomain_lateral%ugrid%dvOnGrid_local(maxEdges, ldomain_lateral%ugrid%ngrid_local))

    call VecGetArrayF90(attr_loc_vec, real_ptr, ierr); CHKERRQ(ierr);
    count = 0
    do ii = 1, ldomain_lateral%ugrid%ngrid_local

       do iedge = 1, maxEdges
          count = count + 1
          ldomain_lateral%ugrid%dcOnGrid_local(iedge, ii) = real_ptr(count)
       enddo

       do iedge = 1, maxEdges
          count = count + 1
          ldomain_lateral%ugrid%dvOnGrid_local(iedge, ii) = real_ptr(count)
       enddo
    enddo
    call VecRestoreArrayF90(attr_loc_vec, real_ptr, ierr); CHKERRQ(ierr);

    call VecDestroy(attr_loc_vec, ierr); CHKERRQ(ierr);
    call VecDestroy(attr_glb_vec, ierr); CHKERRQ(ierr);

    !
    ! areaCell
    !
    allocate(ldomain_lateral%ugrid%areaGrid_ghosted(ldomain_lateral%ugrid%ngrid_ghosted))

    nblocks = 1
    call VecCreate(mpicom, attr_glb_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(attr_glb_vec, nCells_loc_old*nblocks, PETSC_DECIDE, ierr);
    CHKERRQ(ierr)
    call VecSetBlockSize(attr_glb_vec, nblocks, ierr);CHKERRQ(ierr)
    call VecSetFromOptions(attr_glb_vec, ierr);CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_SELF, attr_loc_vec, ierr); CHKERRQ(ierr)
    call VecSetSizes(attr_loc_vec, ldomain_lateral%ugrid%ngrid_ghosted*nblocks, PETSC_DECIDE, ierr);
    CHKERRQ(ierr)
    call VecSetBlockSize(attr_loc_vec, nblocks, ierr);CHKERRQ(ierr)
    call VecSetFromOptions(attr_loc_vec, ierr);CHKERRQ(ierr)

    allocate(int_array(ldomain_lateral%ugrid%ngrid_ghosted))
    do ii = 1,ldomain_lateral%ugrid%ngrid_ghosted
       int_array(ii) = ii-1
    enddo
    call ISCreateBlock(mpicom, nblocks, ldomain_lateral%ugrid%ngrid_ghosted, int_array, &
         PETSC_COPY_VALUES, is_to, ierr);CHKERRQ(ierr);

    do ii = 1, ldomain_lateral%ugrid%ngrid_ghosted
       int_array(ii) = ldomain_lateral%ugrid%grid_id_norder(ii) - 1
    enddo
    call ISCreateBlock(mpicom, nblocks, ldomain_lateral%ugrid%ngrid_ghosted, int_array, &
         PETSC_COPY_VALUES, is_from, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    call VecGetArrayF90(attr_glb_vec, real_ptr, ierr); CHKERRQ(ierr)
    do ii = 1, nCells_loc_old
       real_ptr(ii) = areaCell_old(ii)
    enddo
    call VecRestoreArrayF90(attr_glb_vec, real_ptr, ierr); CHKERRQ(ierr)

    call VecScatterCreate(attr_glb_vec, is_from, attr_loc_vec, is_to, &
         scatter, ierr); CHKERRQ(ierr)
    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to, ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, attr_glb_vec, attr_loc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr);
    call VecScatterEnd(scatter, attr_glb_vec, attr_loc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr);
    call VecScatterDestroy(scatter, ierr)

    call VecGetArrayF90(attr_loc_vec, real_ptr, ierr); CHKERRQ(ierr)
    do ii = 1, ldomain_lateral%ugrid%ngrid_ghosted
       ldomain_lateral%ugrid%areaGrid_ghosted(ii) = real_ptr(ii)
    enddo
    call VecRestoreArrayF90(attr_loc_vec, real_ptr, ierr); CHKERRQ(ierr)

    call VecDestroy(attr_loc_vec, ierr); CHKERRQ(ierr);
    call VecDestroy(attr_glb_vec, ierr); CHKERRQ(ierr);

  end subroutine save_geometric_attributes

  !-----------------------------------------------------------------------

  subroutine ExchangeColumnLevelGhostData(bounds_proc, nvals_per_col, &
       data_send_col, data_recv_col)
    !
    ! !DESCRIPTION:
    ! - Exchanges column level data between MPI tasks.
    ! - This subroutine must be called from OUTSIDE any loops over clumps
    !
    ! !USES:
    use landunit_varcon , only : max_lunit
    use decompMod       , only : bounds_type, BOUNDS_LEVEL_PROC
    use LandunitType    , only : lun_pp
    use ColumnType      , only : col_pp
    use UnstructuredGridType, only : ScatterDataG2L
    !
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)  :: bounds_proc       ! bound information at processor level
    integer , intent(in )          :: nvals_per_col     ! number of values per grid column
    real(r8), intent(in ), pointer :: data_send_col(:)  ! data to be send by each MPI task
    real(r8), intent(out), pointer :: data_recv_col(:)  ! data received by each MPI task
    !
    ! !LOCAL VARIABLES:
    integer             :: c,g,l,j                 ! indices
    integer             :: cidx, lidx              ! column/landunit index
    integer             :: ltype                   ! landunit type
    integer             :: col_ltype               ! landunit type of the column
    integer             :: ier                     ! error
    integer             :: ndata_send              ! number of data sent by local mpi rank
    integer             :: ndata_recv              ! number of data received by local mpi rank
    integer             :: max_ncol_local          ! maximum number of columns per grid cell for local mpi rank
    integer             :: max_ncol_global         ! maximum number of columns per grid cell across all mpi ranks
    integer             :: nblocks                 ! number of values per grid cell
    integer             :: nvals_col               ! number of values per subgrid category
    integer             :: nvals                   ! number of values per subgrid category + additional values
    integer             :: count                   ! temporary
    integer             :: beg_idx, end_idx        ! begin/end index for accessing values in data_send/data_recv
    integer             :: begc_idx, endc_idx      ! begin/end index for accessing values in data_send/data_recv
    integer, pointer    :: ncol(:)                 ! number of columns in grid cell
    integer, pointer    :: landunit_index(:,:)     ! index of the first landunit of a given landunit_itype within a grid cell
    real(r8) , pointer  :: data_send(:)            ! data sent by local mpi rank
    real(r8) , pointer  :: data_recv(:)            ! data received by local mpi rank
    real(r8) , pointer  :: lun_rank(:)             ! rank of a landunit in a given grid cell for a given landunit type
    real(r8) , pointer  :: grid_count(:)           ! temporary
    integer             :: l_rank                  ! rank of landunit
    integer             :: last_lun_type           ! temporary
    PetscErrorCode      :: ierr                    ! PETSc return

    character(len=*), parameter :: subname = 'ExchangeColumnLevelGhostData'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds_proc%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Compute index of the first landunit for a given landunit_itype within a grid cell
    allocate(landunit_index(bounds_proc%begg_all:bounds_proc%endg_all,max_lunit))
    landunit_index = 0

    do lidx = bounds_proc%begl_all,  bounds_proc%endl_all
       if (landunit_index(lun_pp%gridcell(lidx),lun_pp%itype(lidx)) == 0) then
          landunit_index(lun_pp%gridcell(lidx),lun_pp%itype(lidx)) = lidx
       endif
    enddo

    ! Compute number of columns for each grid cell
    allocate(ncol(bounds_proc%begg:bounds_proc%endg))
    ncol = 0

    max_ncol_local = 0
    do c = bounds_proc%begc, bounds_proc%endc
       g       = col_pp%gridcell(c)
       ncol(g) = ncol(g) + 1
       if (ncol(g) > max_ncol_local) max_ncol_local = ncol(g)
    enddo

    ! Determine the maximum number of columns for a grid cell
    call mpi_allreduce(max_ncol_local, max_ncol_global, 1, MPI_INTEGER, MPI_MAX, mpicom, ier)

    ! Determine the total number of data per subgrid category
    nvals     = nvals_per_col + 2

    ! Determine the number of data to be sent/received by
    ! local mpi rank and allocate memory
    nblocks = max_ncol_global * nvals

    ndata_send = nblocks*(bounds_proc%endg     - bounds_proc%begg     + 1)
    ndata_recv = nblocks*(bounds_proc%endg_all - bounds_proc%begg_all + 1)

    allocate(data_send(ndata_send))
    allocate(data_recv(ndata_recv))

    data_send = -9999.d0
    ! Determine the rank of first landunit for a given grid cell
    ! and given landunit type
    !
    ! NOTE: Assumption is that for subgrid category are contigously allocated
    !       for a given landunit type.
    !
    allocate(lun_rank  (bounds_proc%begl_all:bounds_proc%endl_all))
    allocate(grid_count(bounds_proc%begg_all:bounds_proc%endg_all))

    lun_rank(:)   = 0.d0
    grid_count(:) = 0.d0
    last_lun_type   = -1

    do l = bounds_proc%begl_all, bounds_proc%endl_all
       g             = lun_pp%gridcell(l)

       if (last_lun_type /= lun_pp%itype(l)) then
          grid_count(:) = 0.d0
          last_lun_type = lun_pp%itype(l)
       endif
       grid_count(g) = grid_count(g) + 1.d0
       lun_rank(l)   = grid_count(g)
    enddo

    ! Aggregate the data to send
    ncol = 0
    do c = bounds_proc%begc, bounds_proc%endc

       g = col_pp%gridcell(c)
       l = col_pp%landunit(c)

       beg_idx            = (g-bounds_proc%begg)*nblocks + ncol(g)*nvals + 1
       data_send(beg_idx) = real(lun_pp%itype(l))

       beg_idx            = beg_idx + 1
       data_send(beg_idx) = lun_rank(l)

       beg_idx = beg_idx + 1
       end_idx = beg_idx + nvals_per_col - 1

       begc_idx = (c-bounds_proc%begc)*nvals_per_col + 1
       endc_idx = begc_idx + nvals_per_col -1

       data_send(beg_idx:end_idx) = data_send_col(begc_idx:endc_idx)

       ncol(g) = ncol(g) + 1
    enddo

    ! Scatter: Global-to-Local
    call ScatterDataG2L(ldomain_lateral%ugrid, nblocks, ndata_send, data_send, ndata_recv, data_recv)

    ! Save data for ghost subgrid category
    c = bounds_proc%endc
    do ltype = 1, max_lunit
       do g = bounds_proc%endg + 1, bounds_proc%endg_all

          do cidx = 0, max_ncol_local-1

             beg_idx = (g-bounds_proc%begg)*nblocks + cidx*nvals + 1

             col_ltype = int(data_recv(beg_idx))

             beg_idx  = beg_idx + 1
             l_rank   = int(data_recv(beg_idx))

             if (col_ltype == ltype) then
                c       = c + 1

                beg_idx = beg_idx + 1
                end_idx = beg_idx + nvals_per_col - 1

                begc_idx = (c-bounds_proc%begc)*nvals_per_col + 1
                endc_idx = begc_idx + nvals_per_col -1

                data_recv_col(begc_idx:endc_idx) = data_recv(beg_idx:end_idx)

             endif
          enddo
       enddo
    enddo

    ! Free up memory
    deallocate(ncol           )
    deallocate(lun_rank       )
    deallocate(grid_count     )
    deallocate(landunit_index )
    deallocate(data_send      )
    deallocate(data_recv      )

  end subroutine ExchangeColumnLevelGhostData


#else

  !-----------------------------------------------------------------------
  ! This is a stub for the case when PETSc is unavailable
  !
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use elm_varctl  , only : iulog
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
  subroutine domainlateral_init(domain_l, cellsOnCell_old, edgesOnCell_old, &
       nEdgesOnCell_old, areaCell_old, dcEdge_old, dvEdge_old, &
       nCells_loc_old, nEdges_loc_old, maxEdges)
    !
    ! !ARGUMENTS:
    implicit none
    !
    !
    type(domainlateral_type) :: domain_l                     ! domain datatype
    integer , intent(in)     :: cellsOnCell_old(:,:)         ! grid cell level connectivity information
    integer , intent(in)     :: edgesOnCell_old(:,:)         ! index to determine distance between neighbors from dcEdge [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdgesOnCell_old(:)          ! number of edges                                           [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dcEdge_old(:)                ! distance between neighbors                                [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dvEdge_old(:)                ! distance between vertices                                 [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: areaCell_old(:)              ! area of grid cell                                         [in natural order prior to domain decomposition]
    integer , intent(in)     :: nCells_loc_old               ! number of local cell-to-cell connections                  [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdges_loc_old               ! number of edges                                           [in natural order prior to domain decomposition]
    integer , intent(in)     :: maxEdges                     ! max number of edges/neighbors

    character(len=*), parameter :: subname = 'domainlateral_init'

    call endrun(msg='ERROR ' // trim(subname) //': Requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

  end subroutine domainlateral_init

#endif

end module domainLateralMod
