#ifdef CLM_PFLOTRAN

module Mapping_module

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscviewer.h"
  use petscsys
  use petscvec
  use petscmat


  use PFLOTRAN_Constants_module

  implicit none


  private

  ! Note:
  !
  ! CLM has the following:
  !   (i) 3D subsurface grid (CLM_3DSUB);
  !   (ii) 2D top-cell grid (CLM_2DTOP).
  !   (iii) 2D bottom-cell grid (CLM_2DBOT)
  ! CLM decomposes the 3D subsurface grid across processors in a 2D (i.e.
  ! cells in Z are not split across processors). Thus, the surface cells of
  ! 3D subsurface grid are on the same processors as the 2D surface grid.
  !
  ! PFLOTRAN has the following:
  !   (i) 3D subsurface grid (PF_3DSUB);
  !   (ii) top-face control volumes of 3D subsurface grid (PF_2DTOP);
  !   (iii) bottom face control volumes of 3D subsurface grid (PF_2DBOT);
  ! In PFLOTRAN, control volumes in PF_2DSUB and PF_SRF may reside on different
  ! processors. PF_3DSUB and PF_2DSUB are derived from simulation%realization;
  ! while PF_SRF refers to simulation%surf_realization.

  PetscInt, parameter, public :: CLM_3DSUB_TO_PF_3DSUB       = 1 ! 3D --> 3D
  PetscInt, parameter, public :: PF_3DSUB_TO_CLM_3DSUB       = 2 ! 3D --> 3D

  PetscInt, parameter, public :: CLM_2DTOP_TO_PF_2DTOP       = 3 ! TOP face of 3D cell
  PetscInt, parameter, public :: PF_2DTOP_TO_CLM_2DTOP       = 4 ! TOP face of 3D cell
  PetscInt, parameter, public :: CLM_2DBOT_TO_PF_2DBOT       = 5 ! BOTTOM face of 3D cell
  PetscInt, parameter, public :: PF_2DBOT_TO_CLM_2DBOT       = 6 ! BOTTOM face of 3D cell

  type, public  :: mapping_type

    !
    ! Linear Mapping from Source mesh to Destination mesh is matrix-vector product and
    ! can be written as:
    !
    !                W * s = d                                 Eq[1]
    !
    ! where W - Weight matrix      (nd x ns)
    !       s - Source vector      (ns x 1)
    !       d - Destination vector (nd x 1)
    !
    ! In CLM-PFLOTRAN coupling, s and d vectors are decomposed over multiple processors.
    ! The decomposition of vectors need not be in a contiguous order.
    !
    ! Each processor perfoms a local matrix-vector product:
    !
    !                W_loc * s_dloc = d_loc                    Eq[2]
    !
    !   - Obtains from Global Vec 's', a subset of local 's_dloc'
    !   - Performs a local matrix-vector product 
    !

    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt                       :: id

    ! Note: IDs of source/destination mesh are 0-based

    ! Source mesh
    PetscInt           :: s_ncells_loc              ! # of local source mesh cells present
    PetscInt,pointer   :: s_ids_loc_nidx(:)         ! IDs of local source mesh cells
    PetscInt,pointer   :: s_locids_loc_nidx(:)      ! loal IDs of local source mesh cells

    ! Destination mesh
    PetscInt           :: d_ncells_loc              ! # of local destination mesh cells present
    PetscInt           :: d_ncells_gh               ! # of ghost destination mesh cells present
    PetscInt           :: d_ncells_ghd              ! local+ghost=ghosted destination mesh cells

    ! natuaral-index starting with 0
    PetscInt,pointer   :: d_ids_ghd_nidx(:)         ! IDs of ghosted destination mesh cells present
    PetscInt,pointer   :: d_ids_nidx_sor(:)         ! Sorted Ghosted IDs of destination mesh cells
    PetscInt,pointer   :: d_nGhd2Sor(:)             ! Ghosted to Sorted
    PetscInt,pointer   :: d_nSor2Ghd(:)             ! Sorted to Ghosted
    PetscInt,pointer   :: d_loc_or_gh(:)            ! To flag if a cell is local(1) or ghost(0)

    ! Mapping from Source-to-Destination mesh
    PetscInt           :: s2d_s_ncells              ! # of source cells for mapping
    PetscInt           :: s2d_s_ncells_dis          ! # of "distinct" source cells for mapping
    PetscInt,pointer   :: s2d_s_ids_nidx(:)         ! IDs of source cells for mapping
    PetscInt,pointer   :: s2d_s_ids_nidx_dis(:)     ! IDs of "distinct" source cells for mapping

    ! Compressed Sparse Row (CSR) Matrix
    PetscReal,pointer  :: s2d_wts(:)                ! Wts for mapping
    PetscInt           :: s2d_nwts                  ! Number of wts
    PetscInt,pointer   :: s2d_jcsr(:)               ! J-th entry for CSR
    PetscInt,pointer   :: s2d_icsr(:)               ! I-th entry for CSR
    PetscInt,pointer   :: s2d_nonzero_rcount_csr(:) ! Non-Zero entries within a row

    Mat                :: wts_mat                   ! Sparse matrix for linear mapping
    VecScatter         :: s2d_scat_s_gb2disloc      ! Vec-Scatter of source mesh:Global to "distinct" local
    
    Vec                :: s_disloc_vec              ! Sequential vector to save "distinct" local
                                                    ! component of source vector

    ! Header information about number of layers mapped
    PetscInt           :: clm_nlevsoi               ! Number of CLM nlevsoi
    PetscInt           :: clm_nlevgrnd              ! Number of CLM nlevgrnd
    PetscInt           :: clm_nlev_mapped           ! Number of CLM layers mapped
    PetscInt           :: pflotran_nlev             ! Number of PFLOTRAN layers
    PetscInt           :: pflotran_nlev_mapped      ! Number of PFLOTRAN layers mapped

  end type mapping_type

  public :: MappingCreate, &
            MappingSetSourceMeshCellIds, &
            MappingSetDestinationMeshCellIds, &
            MappingFromCLMGrids, &
            MappingReadTxtFile, &
            MappingReadHDF5, &
            MappingDecompose, &
            MappingFindDistinctSourceMeshCellIds, &
            MappingCreateWeightMatrix, &
            MappingCreateScatterOfSourceMesh, &
            MappingFreeNotNeeded, &
            MappingSourceToDestination, &
            MappingDestroy
contains

! ************************************************************************** !

  function MappingCreate()
  ! 
  ! This routine creates a mapping.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 

    implicit none

    type(mapping_type), pointer :: MappingCreate
    type(mapping_type), pointer :: map

    allocate(map)

    map%filename = ''
    map%id = -999
    map%s_ncells_loc = 0
    nullify(map%s_ids_loc_nidx)
    nullify(map%s_locids_loc_nidx)

    ! Destination mesh
    map%d_ncells_loc = 0
    map%d_ncells_gh  = 0
    map%d_ncells_ghd = 0

    ! natuaral-index starting with 0
    nullify(map%d_ids_ghd_nidx)
    nullify(map%d_ids_nidx_sor)
    nullify(map%d_nGhd2Sor)
    nullify(map%d_nSor2Ghd)
    nullify(map%d_loc_or_gh)

    ! Mapping from Source-to-Destination mesh
    map%s2d_s_ncells = 0
    map%s2d_s_ncells_dis = 0
    nullify(map%s2d_s_ids_nidx)
    nullify(map%s2d_s_ids_nidx_dis)

    ! Compressed Sparse Row (CSR) Matrix
    nullify(map%s2d_wts)
    map%s2d_nwts = 0
    nullify(map%s2d_jcsr)
    nullify(map%s2d_icsr)
    nullify(map%s2d_nonzero_rcount_csr)

    map%wts_mat              = PETSC_NULL_MAT
    map%s2d_scat_s_gb2disloc = PETSC_NULL_VECSCATTER
    map%s_disloc_vec         = PETSC_NULL_VEC

    map%clm_nlevsoi          = 0
    map%clm_nlevgrnd         = 0
    map%clm_nlev_mapped      = 0
    map%pflotran_nlev        = 0
    map%pflotran_nlev_mapped = 0
    
    MappingCreate => map

  end function MappingCreate

! ************************************************************************** !

  subroutine MappingFromCLMGrids(map,grid,if_dest_pf,option)
  !
  ! This routine directly obtains grids/soils from CLM.
  ! NOTE: only works with structured/cartsian type of pflotran mesh
  ! (modified from subroutine:: MappingReadTxtFile below
  ! Author: Fengming Yuan, ORNL
  ! Date: Jan. 2016
  !

    use Option_module
    use String_module
    use Grid_module
    use Grid_Structured_module
    use PFLOTRAN_Constants_module

    implicit none

    ! argument
    type(mapping_type), pointer     :: map
    type(grid_type), pointer        :: grid
    PetscBool                       :: if_dest_pf   ! True: PF mesh as destination mesh
    type(option_type), pointer      :: option

    ! local variables
    PetscInt                        :: temp_int
    PetscInt                        :: nwts
    PetscInt,pointer                :: wts_clmid(:), wts_pfid(:)
    PetscInt                        :: i, j, k, natural_id, grid_count, cell_count
    PetscInt,pointer                :: wts_row_tmp(:), wts_col_tmp(:)
    PetscReal,pointer               :: wts_tmp(:)

    PetscErrorCode                  :: ierr
    PetscMPIInt                     :: status_mpi(MPI_STATUS_SIZE)
    PetscInt                        :: nread, cum_nread, nwts_tmp, remainder, irank, ii

    ! Do the 'entire' domain of grid-mesh conversion/scaling through io_rank and communicate to other ranks

    select case(grid%itype)
      case(STRUCTURED_GRID)
        option%io_buffer = 'MAPPING INFO: CLM grid/column dimension will be mapped into ' // &
           'PF structured/cartesian grid mesh.'
        call printMsg(option)

        select case(grid%structured_grid%itype)
          case(CARTESIAN_GRID)
             !
          case default
            option%io_buffer = "MAPPING ERROR: Currently only works on structured CARTESIAN_GRID mesh."
             call printErrMsg(option)
        end select

      case default
        option%io_buffer = "MAPPING ERROR: Currently only works on structured grids."
        call printErrMsg(option)
    end select

    !
    nwts = -999
    if (map%clm_nlev_mapped == map%pflotran_nlev_mapped) then
      if (map%id == CLM_3DSUB_TO_PF_3DSUB .or. map%id == PF_3DSUB_TO_CLM_3DSUB) nwts = grid%structured_grid%nmax
      if (map%id == CLM_2DTOP_TO_PF_2DTOP .or. map%id == PF_2DTOP_TO_CLM_2DTOP) nwts = grid%structured_grid%nxy
      if (map%id == CLM_2DBOT_TO_PF_2DBOT .or. map%id == PF_2DBOT_TO_CLM_2DBOT) nwts = grid%structured_grid%nxy
    else
      option%io_buffer = 'not equal numbers of soil layers mapped btw CLM and PFLOTRAN! '
      call printErrMsg(option)
    end if

    !
    if(option%myrank == option%io_rank) then
      ! some information on CLM-PFLOTRAN mesh matching-up (passing)
      write(*,*) 'CLM nlevsoil mapped = ', map%clm_nlev_mapped
      write(*,*) 'PFLOTRAN nlevsoil mapped = ', map%pflotran_nlev_mapped
      write(*,*) 'CLM-PFLOTRAN ncell mapped = ', nwts

      ! obtain the global cell natural ids for both CLM and PFLTORAN
      allocate(wts_pfid(nwts))
      allocate(wts_clmid(nwts))
      cell_count = 0

      do k = 1, grid%structured_grid%nz
        do j = 1, grid%structured_grid%ny
          do i = 1, grid%structured_grid%nx

            ! PF global cell id in natural order (1-based)
            natural_id = i +                               &   ! x first
                         (j-1)*grid%structured_grid%nx  +  &   ! y second
                         (k-1)*grid%structured_grid%nxy        ! z third


            ! CLM global grid numering (1-based)
            grid_count = i +                          &        ! west-east (x: longitudal) direction third
                        (j-1)*grid%structured_grid%nx          ! south-north (y: latitudal) direction second

            ! cell ids globally
            if (map%id == CLM_3DSUB_TO_PF_3DSUB .or. map%id == PF_3DSUB_TO_CLM_3DSUB) then  ! 3D subsurface (soil) domain
              cell_count = cell_count + 1
              wts_pfid(cell_count)  = natural_id
              wts_clmid(cell_count) = (grid_count-1) * grid%structured_grid%nz + &  ! soil layer numbering first
                                      grid%structured_grid%nz-k+1                   ! reverse vertical-numbering

            elseif((map%id == CLM_2DTOP_TO_PF_2DTOP .or. map%id == PF_2DTOP_TO_CLM_2DTOP) &
               .and. k==grid%structured_grid%nz) then  ! 2D top face
              cell_count = cell_count + 1
              wts_pfid(cell_count)  = natural_id
              wts_clmid(cell_count) = (grid_count-1) * grid%structured_grid%nz + &  ! soil layer numbering first
                                      grid%structured_grid%nz-k+1                   ! reverse vertical-numbering

            elseif((map%id == CLM_2DBOT_TO_PF_2DBOT .or. map%id == PF_2DBOT_TO_CLM_2DBOT) &
               .and. k==1) then  ! 2D bottom face
              cell_count = cell_count + 1
              wts_pfid(cell_count)  = natural_id
              wts_clmid(cell_count) = (grid_count-1) * grid%structured_grid%nz + &  ! soil layer numbering first
                                      grid%structured_grid%nz-k+1                   ! reverse vertical-numbering

            endif

          end do   ! x
        end do   ! y
      end do   ! z
      ! just in case, do the following checking
      if (cell_count /= nwts) then
        option%io_buffer = 'not equal numbers of cell_counts mapped btw CLM and PFLOTRAN! '
        call printErrMsg(option)
      endif

      ! sort destination mesh id obtained above
      ! because mapping approach requires 'destination mesh' id (row-id in wts_matrix) in descending order
      if (if_dest_pf) then
        ! PF mesh as destination mesh (row-id in wts_matrix), while CLM mesh is source-mesh
        call PetscSortIntWithArray(nwts,wts_pfid,wts_clmid,ierr)

      else
        ! CLM mesh as destination mesh (row-id in wts_matrix), while PF mesh is source-mesh
        call PetscSortIntWithArray(nwts,wts_clmid,wts_pfid,ierr)
      end if


      ! calculting mapping matrix/vecs to be send to each processor equally (almost)
      nwts_tmp = nwts/option%mycommsize
      remainder= nwts - nwts_tmp*option%mycommsize

      allocate(wts_row_tmp(nwts_tmp+1))   ! row number is for destination mesh
      allocate(wts_col_tmp(nwts_tmp+1))
      allocate(wts_tmp(nwts_tmp+1))

      cum_nread = 0
      do irank = 0,option%mycommsize-1

        ! Determine the number of rows
        nread = nwts_tmp
        if(irank<remainder) nread = nread+1
        cum_nread = cum_nread + nread


        if (irank == option%myrank) then
          ! Save data locally, if io_rank
          map%s2d_nwts = nread
          allocate(map%s2d_icsr(map%s2d_nwts))
          allocate(map%s2d_jcsr(map%s2d_nwts))
          allocate(map%s2d_wts( map%s2d_nwts))

          do ii = 1, nread
            if (if_dest_pf) then
              map%s2d_icsr(ii) = wts_pfid(cum_nread-nread+ii) - 1        ! 1-based --> 0-based
              map%s2d_jcsr(ii) = wts_clmid(cum_nread-nread+ii) - 1
            else
              map%s2d_icsr(ii) = wts_clmid(cum_nread-nread+ii) - 1        ! 1-based --> 0-based
              map%s2d_jcsr(ii) = wts_pfid(cum_nread-nread+ii) - 1
            end if
           ! set 'wts' to 1.0, assuming 1:1 mapping btw CLM and PF
           ! i.e. no scaling or weighting of area/volume.
           ! SO, all data passing MUST be in unit of per area/volume or independent of area/volume.
            map%s2d_wts(ii)  = 1.d0

          enddo

        else
          ! Otherwise send data to other ranks

          do ii = 1, nread
            if (if_dest_pf) then
              wts_row_tmp(ii) = wts_pfid(cum_nread-nread+ii) - 1        ! 1-based --> 0-based
              wts_col_tmp(ii) = wts_clmid(cum_nread-nread+ii) - 1
            else
              wts_row_tmp(ii) = wts_clmid(cum_nread-nread+ii) - 1        ! 1-based --> 0-based
              wts_col_tmp(ii) = wts_pfid(cum_nread-nread+ii) - 1
            end if
           ! set 'wts' to 1.0, assuming 1:1 mapping btw CLM and PF
           ! i.e. no scaling or weighting of area/volume.
           ! SO, all data passing MUST be in unit of per area/volume or independent of area/volume.
            wts_tmp(ii)  = 1.d0

          enddo

          call MPI_Send(nread,1,MPI_INTEGER,irank,option%myrank,option%mycomm, &
                        ierr)
          call MPI_Send(wts_row_tmp,nread,MPI_INTEGER,irank,option%myrank, &
                        option%mycomm,ierr)
          call MPI_Send(wts_col_tmp,nread,MPI_INTEGER,irank,option%myrank, &
                        option%mycomm,ierr)
          call MPI_Send(wts_tmp,nread,MPI_DOUBLE_PRECISION,irank,option%myrank, &
                        option%mycomm,ierr)

        endif

      end do

      deallocate(wts_row_tmp)
      deallocate(wts_col_tmp)
      deallocate(wts_tmp)
      deallocate(wts_pfid)
      deallocate(wts_clmid)

    else
      ! Other ranks receive data from io_rank

      call MPI_Recv(nread,1,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)

      map%s2d_nwts = nread
      allocate(map%s2d_icsr(map%s2d_nwts))
      allocate(map%s2d_jcsr(map%s2d_nwts))
      allocate(map%s2d_wts( map%s2d_nwts))

      call MPI_Recv(map%s2d_icsr,nread,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
      call MPI_Recv(map%s2d_jcsr,nread,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
      call MPI_Recv(map%s2d_wts,nread,MPI_DOUBLE_PRECISION,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)

    end if !else of if(option%myrank == option%io_rank)

  end subroutine MappingFromCLMGrids

! ************************************************************************** !

  subroutine MappingReadTxtFile(map,map_filename,option)
  ! 
  ! This routine reads a ASCII mapping file.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Input_Aux_module
    use Option_module
    use String_module
    use PFLOTRAN_Constants_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer     :: map
    character(len=MAXSTRINGLENGTH)  :: map_filename
    type(option_type), pointer      :: option
    
    ! local variables
    type(input_type), pointer       :: input
    character(len=MAXSTRINGLENGTH)  :: card
    character(len=MAXSTRINGLENGTH)  :: string
    character(len=MAXWORDLENGTH)    :: word
    character(len=MAXSTRINGLENGTH)  :: hint
    PetscInt                        :: temp_int
    PetscInt                        :: temp_int_array(5)
    PetscInt                        :: nwts
    PetscInt                        :: nread
    PetscInt                        :: remainder, prev_row
    PetscInt                        :: irank, ii
    
    ! Only used by io_rank
    PetscInt                        :: nwts_tmp
    PetscInt,pointer                :: wts_row_tmp(:), wts_col_tmp(:)
    PetscReal,pointer               :: wts_tmp(:)
    PetscInt                        :: nheader
    
    PetscMPIInt                     :: status_mpi(MPI_STATUS_SIZE)
    PetscErrorCode                  :: ierr

    card     = 'MppingReadTxtFile'

    ! Read ASCII file through io_rank and communicate to other ranks
    if(option%myrank == option%io_rank) then

      input => InputCreate(IUNIT_TEMP,map_filename,option)

      nwts     = -1
      prev_row = -1

      ! Read first six entries in the mapping file.
      do nheader = 1, 6
        call InputReadPflotranString(input,option)
        call InputReadWord(input,option,card,PETSC_TRUE)
        call StringToLower(card)

        select case (trim(card))
          case('clm_nlevsoi')
            hint = 'clm_nlevsoi'
            call InputReadInt(input,option,map%clm_nlevsoi)
            call InputErrorMsg(input,option,'CLM nlevsoi',hint)
          case('clm_nlevgrnd')
            hint = 'clm_nlevgrnd'
            call InputReadInt(input,option,map%clm_nlevgrnd)
            call InputErrorMsg(input,option,'CLM nlevgrnd',hint)
          case('clm_nlev_mapped')
            hint = 'clm_nlev_mapped'
            call InputReadInt(input,option,map%clm_nlev_mapped)
            call InputErrorMsg(input,option,'CLM nlev mapped',hint)
          case('pflotran_nlev')
            hint = 'pflotran_nlev'
            call InputReadInt(input,option,map%pflotran_nlev)
            call InputErrorMsg(input,option,'PFLOTRAN nlev',hint)
          case('pflotran_nlev_mapped')
            hint = 'pflotran_nlev_mapped'
            call InputReadInt(input,option,map%pflotran_nlev_mapped)
            call InputErrorMsg(input,option,'PFLOTRAN nlev mapped',hint)
          case('num_weights')
            hint = 'num_weights'
            call InputReadInt(input,option,nwts)
            call InputErrorMsg(input,option,'Number of weights',hint)
            write(*,*) 'rank = ', option%myrank, 'nwts = ',nwts
          case default
            option%io_buffer = 'Unrecognized keyword "' // trim(card) // &
              '" in explicit grid file.'
            call printErrMsgByRank(option)
        end select
      enddo
      
      nwts_tmp = nwts/option%mycommsize
      remainder= nwts - nwts_tmp*option%mycommsize
      
      allocate(wts_row_tmp(nwts_tmp + 1))
      allocate(wts_col_tmp(nwts_tmp + 1))
      allocate(wts_tmp(    nwts_tmp + 1))

      do irank = 0,option%mycommsize-1
        
        ! Determine the number of row to be read
        nread = nwts_tmp
        if(irank<remainder) nread = nread+1
        
        ! Read the data
        do ii = 1,nread
          call InputReadPflotranString(input,option)
          call InputReadInt(input,option,wts_row_tmp(ii))
          call InputReadInt(input,option,wts_col_tmp(ii))
          call InputReadDouble(input,option,wts_tmp(ii))
          
          !Perform checks on the data read
          if(wts_row_tmp(ii) < 1) then
            write(*,string) 'Row entry for ii = ',ii,' less than 1'
            option%io_buffer = string
            call printErrMsg(option)
          endif
          
          if(wts_col_tmp(ii) < 1) then
            write(*,string) 'Col entry for ii = ',ii,' less than 1'
            option%io_buffer = string
            call printErrMsg(option)
          endif
          
          if((wts_tmp(ii) < 0.d0).or.(wts_tmp(ii) > 1.d0)) then
            write(*,string) 'Invalid wt value for ii = ',ii
            option%io_buffer = string
            call printErrMsg(option)
          endif
          
          ! ensure that row values in the data are stored in ascending order
          if(wts_row_tmp(ii) < prev_row) then
            write(*,string) 'Row value in the mapping data not store in ascending order: ii ',ii
            option%io_buffer = string
            call printErrMsg(option) 
          endif
          prev_row = wts_row_tmp(ii)
          
          ! Convert row/col values to 0-based
          wts_row_tmp(ii) = wts_row_tmp(ii) - 1
          wts_col_tmp(ii) = wts_col_tmp(ii) - 1
          
        enddo
        
        ! Save data locally
        if (irank == option%myrank) then
          
            map%s2d_nwts = nread
            allocate(map%s2d_icsr(map%s2d_nwts))
            allocate(map%s2d_jcsr(map%s2d_nwts))
            allocate(map%s2d_wts( map%s2d_nwts))
            
            do ii = 1,map%s2d_nwts
              map%s2d_icsr(ii) = wts_row_tmp(ii)
              map%s2d_jcsr(ii) = wts_col_tmp(ii)
              map%s2d_wts(ii)  = wts_tmp(ii)
            enddo
        
        else

          ! Otherwise communicate data to other ranks
          call MPI_Send(nread,1,MPI_INTEGER,irank,option%myrank,option%mycomm, &
                        ierr)
          call MPI_Send(wts_row_tmp,nread,MPI_INTEGER,irank,option%myrank, &
                        option%mycomm,ierr)
          call MPI_Send(wts_col_tmp,nread,MPI_INTEGER,irank,option%myrank, &
                        option%mycomm,ierr)
          call MPI_Send(wts_tmp,nread,MPI_DOUBLE_PRECISION,irank,option%myrank, &
                        option%mycomm,ierr)
          
        endif
        
      enddo

      deallocate(wts_row_tmp)
      deallocate(wts_col_tmp)
      deallocate(wts_tmp)
      call InputDestroy(input)
      
    else
      ! Other ranks receive data from io_rank
      
      ! Get the number of data
      call MPI_Recv(map%s2d_nwts,1,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)

      ! Allocate memory
      allocate(map%s2d_icsr(map%s2d_nwts))
      allocate(map%s2d_jcsr(map%s2d_nwts))
      allocate(map%s2d_wts( map%s2d_nwts))
      
      call MPI_Recv(map%s2d_icsr,map%s2d_nwts,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
      call MPI_Recv(map%s2d_jcsr,map%s2d_nwts,MPI_INTEGER,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
      call MPI_Recv(map%s2d_wts,map%s2d_nwts,MPI_DOUBLE_PRECISION,option%io_rank,MPI_ANY_TAG, &
                    option%mycomm,status_mpi,ierr)
                    
    endif

    ! Broadcast from root information regarding CLM/PFLOTRAN num soil layers
    temp_int_array(1) = map%clm_nlevsoi
    temp_int_array(2) = map%clm_nlevgrnd
    temp_int_array(3) = map%clm_nlev_mapped
    temp_int_array(4) = map%pflotran_nlev
    temp_int_array(5) = map%pflotran_nlev_mapped

    call MPI_Bcast(temp_int_array,FIVE_INTEGER,MPI_INTEGER,option%io_rank, &
                 option%mycomm,ierr)
    
    map%clm_nlevsoi = temp_int_array(1)
    map%clm_nlevgrnd = temp_int_array(2)
    map%clm_nlev_mapped = temp_int_array(3)
    map%pflotran_nlev = temp_int_array(4)
    map%pflotran_nlev_mapped = temp_int_array(5)

  end subroutine MappingReadTxtFile

! ************************************************************************** !

  subroutine MappingReadHDF5(map,map_filename,option)
  ! 
  ! This routine reads a mapping file in HDF5 format.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/8/2013
  ! 
  
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

    use Input_Aux_module
    use Option_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer :: map
    character(len=MAXSTRINGLENGTH) :: map_filename
    type(option_type), pointer :: option
    
    ! local
    PetscMPIInt       :: hdf5_err
    PetscMPIInt       :: rank_mpi
    PetscInt          :: ndims
    PetscInt          :: istart, iend, ii, jj
    PetscInt          :: num_cells_local
    PetscInt          :: num_cells_local_save
    PetscInt          :: num_vertices_local
    PetscInt          :: num_vertices_local_save
    PetscInt          :: remainder
    PetscInt,pointer  :: int_buffer(:)
    PetscReal,pointer :: double_buffer(:)
    PetscInt, parameter :: max_nvert_per_cell = 8  
    PetscErrorCode    :: ierr

    character(len=MAXSTRINGLENGTH) :: group_name
    character(len=MAXSTRINGLENGTH) :: dataset_name

#if defined(PETSC_HAVE_HDF5)
    integer(HID_T) :: file_id
    integer(HID_T) :: grp_id, grp_id2
    integer(HID_T) :: prop_id
    integer(HID_T) :: data_set_id
    integer(HID_T) :: file_space_id
    integer(HID_T) :: data_space_id
    integer(HID_T) :: memory_space_id
    integer(HSIZE_T) :: num_data_in_file
    integer(HSIZE_T) :: dims_h5(1), max_dims_h5(1)
    integer(HSIZE_T) :: offset(1), length(1), stride(1), block(1), dims(1)
#endif

    ! Initialize FORTRAN predefined datatypes
    call h5open_f(hdf5_err)

    ! Setup file access property with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)

#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif

    ! Open the file collectively
    call h5fopen_f(map_filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
    call h5pclose_f(prop_id,hdf5_err)
    
    !
    ! /row
    !
    
    ! Open group
    group_name = "/row"
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call printMsg(option)

    ! Open dataset
    call h5dopen_f(file_id,"row",data_set_id,hdf5_err)

    ! Get dataset's dataspace
    call h5dget_space_f(data_set_id,data_space_id,hdf5_err)
    
    ! Get number of dimensions and check
    call h5sget_simple_extent_ndims_f(data_space_id,ndims,hdf5_err)
    if (ndims /= 1) then
      option%io_buffer='Dimension of row dataset in ' // trim(map_filename) // &
            ' is not equal to 1.'
      call printErrMsg(option)
    endif

    ! Get dimensions of dataset
    call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5, &
                                     hdf5_err)
    
    ! Determine the number of cells each that will be saved on each processor
    map%s2d_nwts=INT(dims_h5(1))/option%mycommsize
    remainder=INT(dims_h5(1))-map%s2d_nwts*option%mycommsize
    if (option%myrank < remainder) map%s2d_nwts=map%s2d_nwts + 1
    
    ! allocate array to store vertices for each cell
    allocate(map%s2d_icsr(map%s2d_nwts))
    allocate(map%s2d_jcsr(map%s2d_nwts))
    allocate(map%s2d_wts( map%s2d_nwts))

    ! Find istart and iend
    istart = 0
    iend   = 0
    call MPI_Exscan(map%s2d_nwts,istart,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(map%s2d_nwts,iend,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    
    ! Determine the length and offset of data to be read by each processor
    length(1) = iend-istart
    offset(1) = istart
    
    !
    rank_mpi = 1
    memory_space_id = -1
    
    ! Create data space for dataset
    call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err)

    ! Select hyperslab
    call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
    call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                               hdf5_err)
    
    ! Initialize data buffer
    allocate(int_buffer(length(1)))
    
    ! Create property list
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
    
    ! Read the dataset collectively
    call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, int_buffer, &
                   dims_h5, hdf5_err, memory_space_id, data_space_id)
    
    ! Convert 1-based to 0-based
    map%s2d_icsr = int_buffer-1

    call h5dclose_f(data_set_id, hdf5_err)

    !
    ! /col
    !
    
    ! Open group
    group_name = "/col"
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call printMsg(option)

    ! Open dataset
    call h5dopen_f(file_id,"col",data_set_id,hdf5_err)

    ! Get dataset's dataspace
    call h5dget_space_f(data_set_id,data_space_id,hdf5_err)
    
    ! Get number of dimensions and check
    call h5sget_simple_extent_ndims_f(data_space_id,ndims,hdf5_err)
    if (ndims /= 1) then
      option%io_buffer='Dimension of row dataset in ' // trim(map_filename) // &
            ' is not equal to 1.'
      call printErrMsg(option)
    endif

    ! Get dimensions of dataset
    call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5, &
                                     hdf5_err)
    
    ! Determine the number of cells each that will be saved on each processor
    map%s2d_nwts=INT(dims_h5(1))/option%mycommsize
    remainder=INT(dims_h5(1))-map%s2d_nwts*option%mycommsize
    if (option%myrank < remainder) map%s2d_nwts=map%s2d_nwts + 1
    
    ! Find istart and iend
    istart = 0
    iend   = 0
    call MPI_Exscan(map%s2d_nwts,istart,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(map%s2d_nwts,iend,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    
    ! Determine the length and offset of data to be read by each processor
    length(1) = iend-istart
    offset(1) = istart
    
    !
    rank_mpi = 1
    memory_space_id = -1
    
    ! Create data space for dataset
    call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err)

    ! Select hyperslab
    call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
    call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                               hdf5_err)
    
    ! Create property list
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
    
    ! Read the dataset collectively
    call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, int_buffer, &
                   dims_h5, hdf5_err, memory_space_id, data_space_id)
    
    ! Convert 1-based to 0-based
    map%s2d_jcsr = int_buffer-1

    call h5dclose_f(data_set_id, hdf5_err)

    !
    ! /S
    !
    
    ! Open group
    group_name = "/S"
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call printMsg(option)

    ! Open dataset
    call h5dopen_f(file_id,"S",data_set_id,hdf5_err)

    ! Get dataset's dataspace
    call h5dget_space_f(data_set_id,data_space_id,hdf5_err)
    
    ! Get number of dimensions and check
    call h5sget_simple_extent_ndims_f(data_space_id,ndims,hdf5_err)
    if (ndims /= 1) then
      option%io_buffer='Dimension of row dataset in ' // trim(map_filename) // &
            ' is not equal to 1.'
      call printErrMsg(option)
    endif

    ! Get dimensions of dataset
    call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5, &
                                     hdf5_err)
    
    ! Determine the number of cells each that will be saved on each processor
    map%s2d_nwts=INT(dims_h5(1))/option%mycommsize
    remainder=INT(dims_h5(1))-map%s2d_nwts*option%mycommsize
    if (option%myrank < remainder) map%s2d_nwts=map%s2d_nwts + 1
    
    ! Find istart and iend
    istart = 0
    iend   = 0
    call MPI_Exscan(map%s2d_nwts,istart,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(map%s2d_nwts,iend,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
    
    ! Determine the length and offset of data to be read by each processor
    length(1) = iend-istart
    offset(1) = istart
    
    !
    rank_mpi = 1
    memory_space_id = -1
    
    ! Create data space for dataset
    call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err)

    ! Select hyperslab
    call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
    call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                               hdf5_err)
    
    ! Initialize data buffer
    allocate(double_buffer(length(1)))
    
    ! Create property list
    call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
    
    ! Read the dataset collectively
    call h5dread_f(data_set_id, H5T_NATIVE_DOUBLE, double_buffer, &
                   dims_h5, hdf5_err, memory_space_id, data_space_id)
    
    map%s2d_wts = double_buffer

    call h5dclose_f(data_set_id, hdf5_err)

    ! Close file
    call h5fclose_f(file_id, hdf5_err)
    call h5close_f(hdf5_err)

    ! Free memory
    deallocate(int_buffer)
    deallocate(double_buffer)

  end subroutine MappingReadHDF5

! ************************************************************************** !

  subroutine MappingSetSourceMeshCellIds( map,                           &
                                          option,                        &
                                          ncells_loc, ncells_gh,         &
                                          cell_ids_ghd,                  &
                                          local_ids_ghd                  &
                                        )
  !
  ! This routine sets cell ids source mesh.
  !
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  !

    use Option_module
    implicit none

    type(mapping_type), pointer :: map
    type(option_type),  pointer :: option

    PetscInt                    :: ncells_loc, ncells_gh, iloc, ii
    PetscInt,pointer            :: cell_ids_ghd(:), local_ids_ghd(:)

    map%s_ncells_loc = ncells_loc
    allocate(map%s_ids_loc_nidx(map%s_ncells_loc))
    allocate(map%s_locids_loc_nidx(map%s_ncells_loc))

    iloc = 0
    do ii = 1,ncells_loc+ncells_gh
      if (local_ids_ghd(ii)>=1) then
        iloc = iloc+1
        map%s_ids_loc_nidx(iloc)    = cell_ids_ghd(ii)
        map%s_locids_loc_nidx(iloc) = local_ids_ghd(ii)
      endif
    enddo

  end subroutine MappingSetSourceMeshCellIds

! ************************************************************************** !

  subroutine MappingSetDestinationMeshCellIds(map, &
                                              option, &
                                              ncells_loc, &
                                              ncells_gh, &
                                              cell_ids_ghd, &
                                              loc_or_gh &
                                              )
  !
  ! This routine sets the cell ids of destination mesh.
  !
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  !

    use Utility_module, only : DeallocateArray
    use Option_module
    implicit none

    type(mapping_type), pointer :: map
    type(option_type), pointer  :: option

    PetscInt                    :: ncells_loc, ncells_gh
    PetscInt,pointer            :: cell_ids_ghd(:), loc_or_gh(:)

    PetscInt                    :: ii
    PetscInt, pointer           :: index(:),rev_index(:)
    PetscErrorCode              :: ierr

    ! Initialize
    map%d_ncells_loc = ncells_loc
    map%d_ncells_gh  = ncells_gh
    map%d_ncells_ghd = ncells_loc + ncells_gh

    ! Allocate memory
    allocate(map%d_ids_ghd_nidx(map%d_ncells_ghd))
    allocate(map%d_ids_nidx_sor(map%d_ncells_ghd))
    allocate(map%d_loc_or_gh(   map%d_ncells_ghd))
    allocate(map%d_nGhd2Sor(    map%d_ncells_ghd))
    allocate(map%d_nSor2Ghd(    map%d_ncells_ghd))
    allocate(index(             map%d_ncells_ghd))
    allocate(rev_index(         map%d_ncells_ghd))

    do ii=1,ncells_loc+ncells_gh
      map%d_ids_ghd_nidx(ii) = cell_ids_ghd(ii)
      map%d_loc_or_gh(ii)    = loc_or_gh(ii)
      index(ii)              = ii
      rev_index(ii)          = ii
    enddo

    ! Sort cell_ids_ghd
    index = index - 1 ! Needs to be 0-based
    call PetscSortIntWithPermutation(map%d_ncells_ghd,cell_ids_ghd,index,ierr)
    index = index + 1

    do ii=1,ncells_loc+ncells_gh
      map%d_ids_nidx_sor(ii) = cell_ids_ghd(index(ii))
      map%d_nGhd2Sor(ii)     = index(ii)
    enddo

    ! Sort the index
    rev_index = rev_index - 1
    call PetscSortIntWithPermutation(map%d_ncells_ghd,index,rev_index,ierr)
    map%d_nSor2Ghd = rev_index

    ! free memory locally allocated
    deallocate(index)
    deallocate(rev_index)

  end subroutine MappingSetDestinationMeshCellIds

! ************************************************************************** !

  subroutine MappingDecompose(map,option)
  ! 
  ! This routine decomposes the mapping when running on more than 1 processor,
  ! while accounting for different domain decomposition of source and
  ! destination grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer   :: map
    type(option_type), pointer    :: option
    
    ! local variables
    Vec                           :: nonzero_row_count_vec           ! MPI
    Vec                           :: nonzero_row_count_loc_vec       ! Seq
    Vec                           :: cumsum_nonzero_row_count_vec    ! MPI
    Vec                           :: cumsum_nonzero_row_count_loc_vec! Seq
    
    Vec                           :: row_vec, col_vec, wts_vec             ! MPI
    Vec                           :: row_loc_vec, col_loc_vec, wts_loc_vec ! Seq
    
    IS                            :: is_from, is_to
    
    VecScatter                    :: vec_scat
    
    PetscViewer                   :: viewer

    PetscInt, pointer             :: row(:)
    PetscScalar,pointer           :: row_count(:)
    PetscInt, pointer             :: tmp_int_array(:)
    PetscInt                      :: ii,jj,kk
    PetscInt                      :: nrow,cumsum_start,count
    
    PetscScalar,pointer           :: vloc1(:),vloc2(:),vloc3(:),vloc4(:)
    PetscErrorCode                :: ierr
    character(len=MAXSTRINGLENGTH):: string

    ! 0) Save dataset related to wts in MPI-Vecs
    call VecCreateMPI(option%mycomm,map%s2d_nwts,PETSC_DECIDE,row_vec,ierr)
    call VecCreateMPI(option%mycomm,map%s2d_nwts,PETSC_DECIDE,col_vec,ierr)
    call VecCreateMPI(option%mycomm,map%s2d_nwts,PETSC_DECIDE,wts_vec,ierr)

    call VecGetArrayF90(row_vec,vloc1,ierr)
    call VecGetArrayF90(col_vec,vloc2,ierr)
    call VecGetArrayF90(wts_vec,vloc3,ierr)

    do ii = 1,map%s2d_nwts
       vloc1(ii) = map%s2d_icsr(ii)
       vloc2(ii) = map%s2d_jcsr(ii)
       vloc3(ii) = map%s2d_wts(ii)
    enddo

    call VecRestoreArrayF90(row_vec,vloc1,ierr)
    call VecRestoreArrayF90(col_vec,vloc2,ierr)
    call VecRestoreArrayF90(wts_vec,vloc3,ierr)

    ! 1) For each cell of destination mesh, find the number of source mesh cells
    !    overlapped.
    !                             OR
    !    Number of non-zero entries for each row of the global W matrix

    ! Create a MPI vector
    call VecCreateMPI(option%mycomm,map%d_ncells_loc,PETSC_DECIDE, &
                      nonzero_row_count_vec,ierr)

    call VecSet(nonzero_row_count_vec,0.d0,ierr)
    
    ! Find non-zero entries for each of the W matrix with the locally saved data
    allocate(row_count(map%s2d_nwts))
    allocate(row(map%s2d_nwts))

    ii = 1
    nrow = 1
    row(nrow)       = map%s2d_icsr(ii) 
    row_count(nrow) = 1.d0
    
    do ii = 2,map%s2d_nwts
      if (map%s2d_icsr(ii) == row(nrow)) then
        row_count(nrow) = row_count(nrow) + 1
      else
        nrow = nrow + 1
        row(nrow)       = map%s2d_icsr(ii) 
        row_count(nrow) = 1
      endif
    enddo

    ! Save values in the MPI vector
    call VecSetValues(nonzero_row_count_vec,nrow,row,row_count,ADD_VALUES,ierr)
    call VecAssemblyBegin(nonzero_row_count_vec,ierr)
    call VecAssemblyEnd(nonzero_row_count_vec,ierr)
    deallocate(row)
    deallocate(row_count)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_nonzero_row_count_vec.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(nonzero_row_count_vec, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    ! 2) Find cummulative sum of the nonzero_row_count_vec

    ! Create a MPI vector
    call VecCreateMPI(option%mycomm,map%d_ncells_loc,PETSC_DECIDE, &
                      cumsum_nonzero_row_count_vec,ierr)
    call VecGetArrayF90(nonzero_row_count_vec,vloc1,ierr)
    call VecGetArrayF90(cumsum_nonzero_row_count_vec,vloc2,ierr)
    
    ii = 1
    vloc2(ii) = vloc1(ii)
    do ii = 2,map%d_ncells_loc
      vloc2(ii) = vloc2(ii-1) + vloc1(ii)
    enddo

    cumsum_start = 0
    call MPI_Exscan(INT(vloc2(map%d_ncells_loc)),cumsum_start,ONE_INTEGER_MPI, &
                    MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

    do ii = 1,map%d_ncells_loc
      vloc2(ii) = vloc2(ii) + cumsum_start
    enddo
    
    call VecRestoreArrayF90(nonzero_row_count_vec,vloc1,ierr)
    call VecRestoreArrayF90(cumsum_nonzero_row_count_vec,vloc2,ierr)

#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_cumsum_nonzero_row_count_vec.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(cumsum_nonzero_row_count_vec, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    ! 3) On a given processor, in order to map a variable from source mesh to
    !    destination mesh, find 
    !    - the number source mesh cells required
    !    - cell ids of source mesh
    !
    !    Use VecScatter() to save portion of nonzero_row_count_vec and
    !    cumsum_nonzero_row_count_vec corresponding to
    !    ghosted (local+ghost) cell ids destination mesh on a given proc.
    
    !
    call VecCreateSeq(PETSC_COMM_SELF,map%d_ncells_ghd, &
                      nonzero_row_count_loc_vec,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,map%d_ncells_ghd, &
                      cumsum_nonzero_row_count_loc_vec,ierr)

    ! Create index sets (IS) for VecScatter()
    allocate(tmp_int_array(map%d_ncells_ghd))
    do ii = 1,map%d_ncells_ghd
      tmp_int_array(ii) = ii - 1
    enddo
    call ISCreateBlock(option%mycomm,1,map%d_ncells_ghd,tmp_int_array, &
                       PETSC_COPY_VALUES,is_to,ierr)
    deallocate(tmp_int_array)
    
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_to1.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_to,viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    call ISCreateBlock(option%mycomm,1,map%d_ncells_ghd, &
                      map%d_ids_ghd_nidx, PETSC_COPY_VALUES,is_from,ierr)

#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_from1.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_from, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Create VecScatter()
    call VecScatterCreate(nonzero_row_count_vec,is_from, &
                          nonzero_row_count_loc_vec,is_to, &
                          vec_scat,ierr)
    call ISDestroy(is_from,ierr)
    call ISDestroy(is_to,ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_vscat1a.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecScatterView(vec_scat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    ! Scatter the data
    call VecScatterBegin(vec_scat,nonzero_row_count_vec, &
                        nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                        ierr)
    call VecScatterEnd(vec_scat,nonzero_row_count_vec, &
                      nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                      ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_vscat1b.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(cumsum_nonzero_row_count_vec, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    call VecScatterBegin(vec_scat,cumsum_nonzero_row_count_vec, &
                        cumsum_nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                        ierr)
    call VecScatterEnd(vec_scat,cumsum_nonzero_row_count_vec, &
                      cumsum_nonzero_row_count_loc_vec,INSERT_VALUES,SCATTER_FORWARD, &
                      ierr)
    call VecScatterDestroy(vec_scat,ierr)

#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_vscat1c.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(cumsum_nonzero_row_count_vec, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif
    
    call VecGetArrayF90(nonzero_row_count_loc_vec,vloc1,ierr)
    call VecGetArrayF90(cumsum_nonzero_row_count_loc_vec,vloc2,ierr)
    
    allocate(map%s2d_nonzero_rcount_csr(map%d_ncells_ghd))
    
    ! For each destination, save the number of overlapped source cells
    count = 0
    do ii = 1,map%d_ncells_ghd
      count = count + INT(vloc1(ii))
      map%s2d_nonzero_rcount_csr(ii) = INT(vloc1(ii))
    enddo
    map%s2d_s_ncells = count

    ! Allocate memory
    if(associated(map%s2d_wts)) deallocate(map%s2d_wts)
    ! this free allocated memory in reading txt (Line 422/Line 461)
    ! or hdf mesh file (Line 601) above
    ! so that same pointer variable can be used but may be changed below

    allocate(map%s2d_s_ids_nidx(map%s2d_s_ncells))
    allocate(map%s2d_wts(map%s2d_s_ncells))
    allocate(tmp_int_array(map%s2d_s_ncells))
      
    ! For each cell in destination mesh, save indices of MPI Vectors, which
    ! contain data read from mapping file, for all overlapped cells of 
    ! of source mesh
    kk = 0
    do ii = 1,map%d_ncells_ghd
       do jj = 1,INT(vloc1(ii))
          kk = kk + 1
          tmp_int_array(kk) = INT(vloc2(ii)) - INT(vloc1(ii)) + jj - 1
       enddo
    enddo

    ! Create an index set to scatter from
    call ISCreateBlock(option%mycomm,1,map%s2d_s_ncells,tmp_int_array, &
         PETSC_COPY_VALUES,is_from,ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_from2.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_from,viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    do ii=1,map%s2d_s_ncells
       tmp_int_array(ii) = ii-1
    enddo

    call ISCreateBlock(option%mycomm,1,map%s2d_s_ncells,tmp_int_array, &
         PETSC_COPY_VALUES,is_to,ierr)
    deallocate(tmp_int_array)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_to2.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_to,viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Allocate memory
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_s_ncells,row_loc_vec,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_s_ncells,col_loc_vec,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_s_ncells,wts_loc_vec,ierr)

    ! Create scatter context
    call VecScatterCreate(row_vec,is_from,row_loc_vec,is_to,vec_scat,ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_vscat2.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecScatterView(vec_scat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Scatter the data
    call VecScatterBegin(vec_scat,col_vec,col_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)
    call VecScatterEnd(vec_scat,col_vec,col_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)
    call VecScatterBegin(vec_scat,wts_vec,wts_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)
    call VecScatterEnd(vec_scat,wts_vec,wts_loc_vec, INSERT_VALUES, &
         SCATTER_FORWARD,ierr)

    ! Attach to the local copy of the scatterd data
    call VecGetArrayF90(col_loc_vec,vloc3,ierr)
    call VecGetArrayF90(wts_loc_vec,vloc4,ierr)

    ! Save the scattered data
    do ii = 1,map%s2d_s_ncells
       map%s2d_s_ids_nidx(ii) = INT(vloc3(ii))
       map%s2d_wts(ii)        = vloc4(ii)
    enddo

    ! Restore data
    call VecRestoreArrayF90(col_loc_vec,vloc3,ierr)
    call VecRestoreArrayF90(wts_loc_vec,vloc4,ierr)

    ! Free memory
    call VecDestroy(row_loc_vec,ierr)
    call VecDestroy(col_loc_vec,ierr)
    call VecDestroy(wts_loc_vec,ierr)
    
    ! Restore data
    call VecRestoreArrayF90(nonzero_row_count_loc_vec,vloc1,ierr)
    call VecRestoreArrayF90(cumsum_nonzero_row_count_loc_vec,vloc2,ierr)

    ! Free memory
    call VecDestroy(nonzero_row_count_vec,ierr)
    call VecDestroy(cumsum_nonzero_row_count_vec,ierr)
    call VecDestroy(nonzero_row_count_loc_vec,ierr)
    call VecDestroy(cumsum_nonzero_row_count_loc_vec,ierr)

  end subroutine MappingDecompose

! ************************************************************************** !

  subroutine MappingFindDistinctSourceMeshCellIds(map,option)
  ! 
  ! This routine finds distinct cell ids of source mesh
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer   :: map
    type(option_type), pointer    :: option
    
    ! local variables
    PetscErrorCode               :: ierr
    PetscInt                     :: ii,jj,kk,count
    PetscInt, pointer            :: index(:)
    PetscInt,pointer             :: int_array(:),int_array2(:)
    PetscInt,pointer             :: int_array3(:),int_array4(:)

    Vec :: xx, yy
    
    ! No overlapped cells with Source Mesh, then return
    if(map%s2d_s_ncells == 0) then
       map%s2d_s_ncells_dis = 0
       call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells_dis, map%s_disloc_vec, ierr)
       return
    end if

    ! Allocate memory
    if(associated(map%s2d_jcsr)) deallocate(map%s2d_jcsr)
    ! it has been allocated before (when reading mesh file either in txt or in hdf5)
    allocate(map%s2d_jcsr( map%s2d_s_ncells))
    allocate(int_array (map%s2d_s_ncells))
    allocate(int_array2(map%s2d_s_ncells))
    allocate(int_array3(map%s2d_s_ncells))
    allocate(int_array4(map%s2d_s_ncells))

    !
    ! Follows Glenn's approach in unstructured code to remove duplicate
    !   vertices.
    !
    ! map%s2d_s_ids_loc_nidx - Contains source mesh cell ids with duplicate entries
    ! The algo is explained below: 
    !     int_array  : Cell-ids
    !     int_array2 : Sorted index
    !     int_array3 : Distinct values
    !     int_array4 : Indices w.r.t. new distinct value vector
    !
    !  ii  int_array  int_array2  int_array3  int_array4
    !   1     90         6           70          3
    !   2    100         3           80          4
    !   3     80         1           90          2
    !   4    100         2          100          4
    !   5    101         4          101          5
    !   6     70         5                       1
    !
  
    do ii = 1,map%s2d_s_ncells
      int_array(ii)  = map%s2d_s_ids_nidx(ii)
      int_array2(ii) = ii
    enddo

    int_array2 = int_array2 - 1
    call PetscSortIntWithPermutation(map%s2d_s_ncells,int_array,int_array2,ierr)
    int_array2 = int_array2 + 1
    
    int_array3 = 0
    int_array4 = 0
    count = 1
    int_array3(1)             = int_array(int_array2(1))
    int_array4(int_array2(1)) = count
    
    do ii=2,map%s2d_s_ncells
      jj = int_array(int_array2(ii))
      if (jj > int_array3(count)) then
        count = count + 1
        int_array3(count) = jj
      endif
      int_array4(int_array2(ii)) = count 
    enddo
    
    ! Change 1-based index to 0-based index
    int_array4 = int_array4 - 1
    
    map%s2d_s_ncells_dis = count

    ! Save the distinct ids
    allocate(map%s2d_s_ids_nidx_dis(map%s2d_s_ncells_dis))
    
    map%s2d_s_ids_nidx_dis(1:count) = int_array3(1:count)
    map%s2d_jcsr = int_array4

    ! Create a sequential vector
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_s_ncells_dis, map%s_disloc_vec, ierr)

    ! Free memory
    deallocate(int_array)
    deallocate(int_array2)
    deallocate(int_array3)
    deallocate(int_array4)

  end subroutine MappingFindDistinctSourceMeshCellIds

! ************************************************************************** !

  subroutine MappingCreateWeightMatrix(map,option)
  ! 
  ! This routine creates a weight matrix to map data from source to destination
  ! grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    
    implicit none
    
    ! argument
    type(mapping_type), pointer   :: map
    type(option_type), pointer    :: option
    
    ! local variables
    PetscInt, pointer            :: index(:)
    PetscInt                     :: ii,jj,kk
    PetscErrorCode               :: ierr
    character(len=MAXSTRINGLENGTH)     :: string, string1
    PetscViewer :: viewer

    allocate(index(map%s2d_s_ncells))

    kk = 0
    do ii = 1,map%d_ncells_ghd
       do jj = 1,map%s2d_nonzero_rcount_csr(ii)
          kk = kk + 1
          index(kk) = ii -1
       enddo
    enddo

    !
    ! size(wts_mat) = [d_ncells_ghd x s2d_s_ncells_dis]
    !
    call MatCreateSeqAIJ(PETSC_COMM_SELF, map%d_ncells_ghd, &
         map%s2d_s_ncells_dis, PETSC_DEFAULT_INTEGER, &    ! PETSC_NULL_INTEGER not works here
         map%s2d_nonzero_rcount_csr, map%wts_mat, ierr)

    do ii = 1,map%s2d_s_ncells
       call MatSetValues(map%wts_mat,1,index(ii),1,map%s2d_jcsr(ii), &
          map%s2d_wts(ii),INSERT_VALUES,ierr)
    enddo
    call MatAssemblyBegin(map%wts_mat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(  map%wts_mat,MAT_FINAL_ASSEMBLY,ierr)
  
    deallocate(index)

#ifdef MAP_DEBUG
    write(string,*) option%myrank
    write(string1,*) map%id
    string = 'mat_wts' // trim(adjustl(string)) // '_' // trim(adjustl(string1)) // '.out'
    call PetscViewerASCIIOpen(PETSC_COMM_SELF, trim(string), viewer, ierr)
    call MatView(map%wts_mat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

  end subroutine MappingCreateWeightMatrix

! ************************************************************************** !

  subroutine MappingCreateScatterOfSourceMesh(map,option)
  ! 
  ! This routine creates a vector scatter context from source to destination
  ! grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    implicit none
#include "petsc/finclude/petscviewer.h"

    ! argument
    type(mapping_type), pointer  :: map
    type(option_type), pointer   :: option

    ! local variables
    Vec :: pindex       ! PETSc index
    Vec :: nindex       ! Natural index
    Vec :: N2P          ! Natural to PETSc index
    Vec :: pindex_req   ! Required PETSc indices
    IS  :: is_from, is_to
    VecScatter :: vscat
    PetscViewer :: viewer

    PetscInt                     :: ii,jj,kk
    PetscInt                     :: istart, iend
    PetscInt,pointer             :: tmp_int_array(:)
    PetscScalar,pointer          :: v_loc(:)
    PetscErrorCode               :: ierr

    character(len=MAXSTRINGLENGTH)     :: string

    !
    ! Example:
    !
    ! GIVEN:
    ! source vector (MPI):
    !                       p0        |    p1         : processor id
    !                [ s3 s4 s5 s6 s7 | s2 s1 s0 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    !
    ! required source vector (Seq):
    !                     p0        |             p1         : processor id
    !                [ s0 s4 s6 ]   | [s0 s1 s2 s3 s4 s5 s6] : natural index
    !                   0  1  2     |   0  1  2  3  4  5  6  : PETSc order
    !
    ! End product of this subroutine is construction of a vector-scatter, which
    ! has indices of MPI vector needed for creation of sequential vector.
    !
    !                     p0        |             p1         : processor id
    !                [  7  1  3 ]   | [7  6  5  0  1  2  4]  : pindex_req 

    ! Create the vectors
    call VecCreateMPI(option%mycomm, map%s_ncells_loc, PETSC_DECIDE, pindex, ierr)
    call VecCreateMPI(option%mycomm, map%s_ncells_loc, PETSC_DECIDE, nindex, ierr)
    call VecCreateMPI(option%mycomm, map%s_ncells_loc, PETSC_DECIDE, N2P  , ierr)
    call VecCreateMPI(option%mycomm, map%s2d_s_ncells_dis, PETSC_DECIDE, pindex_req, ierr)


    ! STEP-1 -
    !
    ! GIVEN:
    ! source vector (MPI):
    !                       p0        |    p1         : processor id
    !                [ s3 s4 s5 s6 s7 | s2 s1 s0 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! source vector (MPI) sorted in asceding order:
    !                       p0        |    p1         : processor id
    !                [ s0 s1 s2 s3 s4 | s5 s6 s7 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! RESULT:
    ! Indices of MPI vector required to do the sorting (N2P)
    !                       p0        |    p1         : processor id
    !                [  7  6  5  0  1 |  2  3  4  ]   : N2P
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    
    ! Initialize 'nindex' vector
    call VecGetArrayF90(nindex,v_loc,ierr)
    do ii=1,map%s_ncells_loc
       v_loc(ii) = map%s_ids_loc_nidx(ii)
    enddo
    call VecRestoreArrayF90(nindex,v_loc,ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_nindex.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(nindex,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Initialize 'pindex' vector
    call VecGetOwnershipRange(pindex,istart,iend,ierr)
    call VecGetArrayF90(pindex,v_loc,ierr)
    do ii=1,map%s_ncells_loc
       v_loc(ii) = istart + ii - 1
    enddo
    call VecRestoreArrayF90(pindex,v_loc,ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_pindex.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(pindex,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Create 'is_from'
    allocate(tmp_int_array(map%s_ncells_loc))
    do ii=1,map%s_ncells_loc
       tmp_int_array(ii) = istart + ii - 1
    enddo
    call ISCreateBlock(option%mycomm, 1, map%s_ncells_loc, tmp_int_array, PETSC_COPY_VALUES, &
         is_from, ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_from3.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_from,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Create 'is_to'
    do ii=1,map%s_ncells_loc
       tmp_int_array(ii) = map%s_ids_loc_nidx(ii)
    enddo
    call ISCreateBlock(option%mycomm, 1, map%s_ncells_loc, tmp_int_array, PETSC_COPY_VALUES, &
         is_to, ierr)
    deallocate(tmp_int_array)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_to3.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_to,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Create 'vscat'
    call VecScatterCreate(nindex, is_from, N2P, is_to, vscat, ierr)
    call ISDestroy(is_to,ierr)
    call ISDestroy(is_from,ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_vscat3.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecScatterView(vscat, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Scatter data
    call VecScatterBegin(vscat, pindex, N2P, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vscat, pindex, N2P, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterDestroy(vscat,ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_N2P.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(N2P,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! STEP-2 -
    !
    ! GIVEN:
    ! source vector (MPI):
    !                       p0        |    p1         : processor id
    !                [ s3 s4 s5 s6 s7 | s2 s1 s0 ]    : natural index
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! Indices of MPI vector required to do the sorting (N2P)
    !                       p0        |    p1         : processor id
    !                [  7  6  5  0  1 |  2  3  4  ]   : N2P
    !                   0  1  2  3  4 |  5  6  7      : PETSc index
    !
    ! required source vector (Seq):
    !                     p0        |             p1         : processor id
    !                [ s0 s4 s6 ]   | [s0 s1 s2 s3 s4 s5 s6] : natural index
    !                   0  1  2     |   0  1  2  3  4  5  6  : PETSc order
    !
    ! RESULT:
    ! Indices of MPI vector needed for creation of sequential vector.
    !                     p0        |             p1         : processor id
    !                [  7  1  3 ]   | [7  6  5  0  1  2  4]  : pindex_req
    !

    ! Create 'is_to'
    call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells_dis, map%s2d_s_ids_nidx_dis, &
         PETSC_COPY_VALUES, is_to, ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_to4.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_to, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Create 'is_from'
    call VecGetOwnershipRange(pindex_req, istart, iend, ierr)
    allocate(tmp_int_array(map%s2d_s_ncells_dis))
    do ii = 1,map%s2d_s_ncells_dis
       tmp_int_array(ii) = istart + ii - 1
    enddo
    call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells_dis, tmp_int_array, &
         PETSC_COPY_VALUES, is_from, ierr)
    deallocate(tmp_int_array)

#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_is_from4.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call ISView(is_from, viewer,ierr)
    call PetscViewerDestroy(viewer, ierr)
#endif

    ! Create vector scatter
    call VecScatterCreate(N2P, is_to, pindex_req, is_from, vscat, ierr)

    call ISDestroy(is_to,ierr)
    call ISDestroy(is_from, ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_vscat4.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecScatterView(vscat,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Scatter data
    call VecScatterBegin(vscat, N2P, pindex_req, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(  vscat, N2P, pindex_req, INSERT_VALUES, SCATTER_FORWARD, ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_pindex_req.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecView(pindex_req,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    !
    ! Step-3 -
    ! Using 'pindex_req' create and save a vector-scatter from MPI source vector
    ! to sequential source vectors

    call VecGetArrayF90(pindex_req,v_loc,ierr)
    call ISCreateBlock(option%mycomm, 1, map%s2d_s_ncells_dis, INT(v_loc),&
         PETSC_COPY_VALUES, is_from, ierr)
    call VecRestoreArrayF90(pindex_req,v_loc,ierr)

    call VecScatterCreate(N2P, is_from, pindex_req, is_to, map%s2d_scat_s_gb2disloc, ierr)
#ifdef MAP_DEBUG
    write(string,*) map%id
    string = 'map' // trim(adjustl(string)) // '_s2d_scat_s_gb2disloc.out'
    call PetscViewerASCIIOpen(option%mycomm, trim(string), viewer, ierr)
    call VecScatterView(map%s2d_scat_s_gb2disloc,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

    ! Free-memory
    call VecDestroy(pindex, ierr)
    call VecDestroy(nindex, ierr)
    call VecDestroy(N2P  , ierr)
    call VecScatterDestroy(vscat,ierr)

  end subroutine MappingCreateScatterOfSourceMesh

! ************************************************************************** !
  ! what needed () in 'MappingSourceToDestination' after initializing mapping is as following 3,
  ! map%s_disloc_vec
  ! map%s2d_scat_s_gb2disloc
  ! map%wts_mat
  ! the rest map%??? (total 12 pointers) can be deallocated, if allocated memory before,
  ! so that free memory not needed any more (and avoid possible memory leak).
  subroutine MappingFreeNotNeeded(map)

    implicit none

    type(mapping_type), pointer :: map

    if(associated(map%s_ids_loc_nidx)) deallocate(map%s_ids_loc_nidx)
          ! allocated in L187, last used in L1409.

    if(associated(map%d_ids_ghd_nidx)) deallocate(map%d_ids_ghd_nidx)
          ! allocated in L234, last used in L968

    ! in 'MappingSetDestinationMeshCellIds', the following may not really used
    if(associated(map%d_ids_nidx_sor)) deallocate(map%d_ids_nidx_sor)   ! allocated in L235, last used in L255
    if(associated(map%d_loc_or_gh)) deallocate(map%d_loc_or_gh)         ! allocated in L236, last used in L244
    if(associated(map%d_nGhd2Sor)) deallocate(map%d_nGhd2Sor)           ! allocated in L237, last used in L256
    if(associated(map%d_nSor2Ghd)) deallocate(map%d_nSor2Ghd)           ! allocated in L238, last used in L262
    if(associated(map%d_nGhd2Sor)) deallocate(map%d_nGhd2Sor)           ! allocated in L237, last used in L256

    if(associated(map%s2d_icsr)) deallocate(map%s2d_icsr)
        ! allocated in L420/459 or L599 in reading mapping mesh file either in txt or in hdf5,
        ! used L888, and in pflotran_model's 'initMapFrom??To??' (btw 2D and 3D mesh)

    if(associated(map%s2d_jcsr)) deallocate(map%s2d_jcsr)
        ! allocated in L421/460 or L600 in reading mapping mesh file either in txt or in hdf5,
        ! used L855, and in pflotran_model's 'initMapFrom??To??' (btw 2D and 3D mesh)
        ! BUT, then re-allocated in L1149 and last used in L1264

    if(associated(map%s2d_wts)) deallocate(map%s2d_wts)
        ! allocated in L422/461 or L601 in reading mapping mesh file either in txt or in hdf5,
        ! last used L856,
        ! BUT, then re-allocated in L1022, last used in L1265

    if(associated(map%s2d_nonzero_rcount_csr)) deallocate(map%s2d_nonzero_rcount_csr)
          ! allocated in L1005, last used in L1261

  end subroutine MappingFreeNotNeeded

! ************************************************************************** !

  subroutine MappingSourceToDestination(map,option,s_vec,d_vec)
  ! 
  ! This routine maps the data from source to destination grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    use Option_module
    implicit none
    
    ! argument
    type(mapping_type), pointer :: map
    type(option_type), pointer  :: option
    Vec                         :: s_vec ! MPI
    Vec                         :: d_vec ! Seq
    
    ! local variables
    PetscErrorCode              :: ierr

    ! a note here:
    ! what needed after initializing mapping is as following,
    ! the rest map%??? can be deallocated, if allocated memory before,
    ! so that free memory not needed any more.
    ! map%s_disloc_vec
    ! map%s2d_scat_s_gb2disloc
    ! map%wts_mat

    if (map%s2d_s_ncells > 0) then  
       ! Initialize local vector
       call VecSet(map%s_disloc_vec, 0.d0, ierr);CHKERRQ(ierr)

       call VecSet(d_vec, 0.d0, ierr);CHKERRQ(ierr)

    end if

    ! Scatter the source vector
    call VecScatterBegin(map%s2d_scat_s_gb2disloc, s_vec, map%s_disloc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)

    call VecScatterEnd(map%s2d_scat_s_gb2disloc, s_vec, map%s_disloc_vec, &
         INSERT_VALUES, SCATTER_FORWARD, ierr)
    CHKERRQ(ierr)

    if (map%s2d_s_ncells > 0) then  
       ! Perform Matrix-Vector product
       call MatMult(map%wts_mat, map%s_disloc_vec, d_vec, ierr)
       CHKERRQ(ierr)
    end if

  end subroutine MappingSourceToDestination

! ************************************************************************** !

  subroutine MappingDestroy(map)
  !
  ! This routine frees up memoery
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/22/2014
  !

    implicit none

    ! argument
    type(mapping_type), pointer :: map

    if (associated(map%s_ids_loc_nidx)) deallocate(map%s_ids_loc_nidx)
    if (associated(map%s_locids_loc_nidx)) deallocate(map%s_locids_loc_nidx)
    if (associated(map%d_ids_ghd_nidx)) deallocate(map%d_ids_ghd_nidx)
    if (associated(map%d_ids_nidx_sor)) deallocate(map%d_ids_nidx_sor)
    if (associated(map%d_nGhd2Sor)) deallocate(map%d_nGhd2Sor)
    if (associated(map%d_nSor2Ghd)) deallocate(map%d_nSor2Ghd)
    if (associated(map%d_loc_or_gh)) deallocate(map%d_loc_or_gh)
    if (associated(map%s2d_s_ids_nidx)) deallocate(map%s2d_s_ids_nidx)
    if (associated(map%s2d_s_ids_nidx_dis)) deallocate(map%s2d_s_ids_nidx_dis)
    if (associated(map%s2d_wts)) deallocate(map%s2d_wts)
    if (associated(map%s2d_jcsr)) deallocate(map%s2d_jcsr)
    if (associated(map%s2d_icsr)) deallocate(map%s2d_icsr)
    if (associated(map%s2d_nonzero_rcount_csr)) deallocate(map%s2d_nonzero_rcount_csr)

    nullify(map%s_ids_loc_nidx)
    nullify(map%d_ids_ghd_nidx)
    nullify(map%d_ids_nidx_sor)
    nullify(map%d_nGhd2Sor)
    nullify(map%d_nSor2Ghd)
    nullify(map%d_loc_or_gh)
    nullify(map%s2d_s_ids_nidx)
    nullify(map%s2d_s_ids_nidx_dis)
    nullify(map%s2d_wts)
    nullify(map%s2d_jcsr)
    nullify(map%s2d_icsr)
    nullify(map%s2d_nonzero_rcount_csr)

  end subroutine MappingDestroy

end module Mapping_module

#endif
