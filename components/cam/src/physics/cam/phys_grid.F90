module phys_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of physics computational horizontal grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code.
! 
! Entry points:
!      phys_grid_init       initialize chunk'ed data structure
!      phys_grid_initialized    get physgrid_set flag
!
!      phys_grid_defaultopts   get default runtime options
!      phys_grid_setopts       set runtime options
!
!      get_chunk_indices_p get local chunk index range
!      get_ncols_p         get number of columns for a given chunk
!      get_xxx_all_p       get global indices, coordinates, or values
!                          for a given chunk
!      get_xxx_vec_p       get global indices, coordinates, or values
!                          for a subset of the columns in a chunk
!      get_xxx_p           get global indices, coordinates, or values
!                          for a single column
!      where xxx is
!       area               for column surface area (in radians squared)
!       gcol               for global column index
!       lat                for global latitude index
!       lon                for global longitude index
!       rlat               for latitude coordinate (in radians)
!       rlon               for longitude coordinate (in radians)
!       wght               for column integration weight
!
!      print_cost_p        print measured cost (for all chunks)
!      update_cost_p       add walltime to chunk cost for a given chunk
!
!      scatter_field_to_chunk
!                          distribute field
!                          to decomposed chunk data structure
!      gather_chunk_to_field
!                          reconstruct field
!                          from decomposed chunk data structure
!
!      read_chunk_from_field
!                          read and distribute field
!                          to decomposed chunk data structure
!      write_field_from_chunk
!                          write field
!                          from decomposed chunk data structure
!
!      block_to_chunk_send_pters
!                          return pointers into send buffer where data
!                          from decomposed fields should
!                          be copied to
!      block_to_chunk_recv_pters
!                          return pointers into receive buffer where data
!                          for decomposed chunk data structures should
!                          be copied from
!      transpose_block_to_chunk
!                          transpose buffer containing decomposed 
!                          fields to buffer
!                          containing decomposed chunk data structures
!
!      chunk_to_block_send_pters
!                          return pointers into send buffer where data
!                          from decomposed chunk data structures should
!                          be copied to
!      chunk_to_block_recv_pters
!                          return pointers into receive buffer where data
!                          for decomposed fields should
!                          be copied from
!      transpose_chunk_to_block
!                          transpose buffer containing decomposed
!                          chunk data structures to buffer
!                          containing decomposed fields
!
!      chunk_index         identify whether index is for a latitude or
!                          a chunk
!
! FOLLOWING ARE NO LONGER USED, AND ARE CURRENTLY COMMENTED OUT
!      get_gcol_owner_p    get owner of column
!                          for given global physics column index
!
!      buff_to_chunk       Copy from local buffer to local chunk data 
!                          structure. (Needed for cpl6.)
!
!      chunk_to_buff       Copy from local chunk data structure to 
!                          local buffer. (Needed for cpl6.)
!
! Author: Patrick Worley and John Drake
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,     only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use physconst,        only: pi
   use ppgrid,           only: pcols, pver, begchunk, endchunk
#if ( defined SPMD )
   use spmd_dyn,         only: block_buf_nrecs, chunk_buf_nrecs, &
                               local_dp_map
   use mpishorthand
#endif
   use spmd_utils,       only: iam, masterproc, npes, proc_smp_map, nsmps
   use m_MergeSorts,     only: IndexSet, IndexSort
   use cam_abortutils,   only: endrun
   use perf_mod
   use cam_logfile,      only: iulog
   use scamMod,          only: single_column, scmlat, scmlon, iop_mode
   use shr_const_mod,    only: SHR_CONST_PI
   use dycore,           only: dycore_is
   use units,            only: getunit, freeunit

   implicit none
   save

#if ( ! defined SPMD )
   integer, private :: block_buf_nrecs
   integer, private :: chunk_buf_nrecs
   logical, private :: local_dp_map=.true. 
#endif

! The identifier for the physics grid
   integer, parameter, public :: phys_decomp = 100

! dynamics field grid information
   integer, private :: hdim1_d, hdim2_d
                                       ! dimensions of rectangular horizontal grid
                                       ! data structure, If 1D data structure, then
                                       ! hdim2_d == 1.
   logical, private :: use_cost_d
                                       ! flag indicating that nontrivial colum cost 
                                       ! estimates are available for use in load 
                                       ! balancing
   real(r8), dimension(:), allocatable, private :: cost_d
                                       ! normalized estimated column computational 
                                       ! cost (from dynamics)

! physics field data structures
   integer         :: ngcols           ! global column count in physics grid (all)
   integer, public :: ngcols_p         ! global column count in physics grid 
                                       ! (without holes)

   integer, dimension(:), allocatable, public :: dyn_to_latlon_gcol_map
                                       ! map from unsorted (dynamics) to lat/lon sorted grid indices
   integer, dimension(:), allocatable, public :: latlon_to_dyn_gcol_map !now beingh used in RRTMG radiation.F90
                                       ! map from lat/lon sorted grid to unsorted (dynamics) indices
   integer, dimension(:), allocatable, private :: lonlat_to_dyn_gcol_map
                                       ! map from lon/lat sorted grid to unsorted (dynamics) indices

!   integer, private :: clat_p_tot ! number of unique latitudes
!   integer, private :: clon_p_tot ! number of unique longitudes
! these are public to support mozart chemistry in the short term
   integer, public :: clat_p_tot ! number of unique latitudes
   integer, public :: clon_p_tot ! number of unique longitudes

   integer, dimension(:), allocatable, private :: clat_p_cnt ! number of repeats for each latitude
   integer, dimension(:), allocatable, private :: clat_p_idx ! index in latlon ordering for first occurence
                                                             ! of latitude corresponding to given 
                                                             ! latitude index
   real(r8), dimension(:), allocatable :: clat_p  ! unique latitudes (radians, increasing)


   integer, dimension(:), allocatable, private :: clon_p_cnt ! number of repeats for each longitude
   integer, dimension(:), allocatable, private :: clon_p_idx ! index in lonlat ordering for first 
                                                             ! occurrence of longitude corresponding to 
                                                             ! given latitude index
   real(r8), dimension(:), allocatable :: clon_p  ! unique longitudes (radians, increasing)

   integer, dimension(:), allocatable, private :: lat_p      ! index into list of unique column latitudes
   integer, dimension(:), allocatable, private :: lon_p      ! index into list of unique column longitudes

! chunk data structures
   type chunk
     integer  :: ncols                 ! number of vertical columns
     integer, allocatable :: gcol(:)   ! global physics column indices
     integer, allocatable :: lon(:)    ! global longitude indices
     integer, allocatable :: lat(:)    ! global latitude indices
     integer  :: owner                 ! id of process where chunk assigned
     integer  :: lcid                  ! local chunk index
     real(r8) :: estcost               ! estimated computational cost (normalized)
   end type chunk

   integer :: nchunks                  ! global chunk count
   type (chunk), dimension(:), allocatable, public :: chunks  
                                       ! global computational grid

!!XXgoldyXX: v this should be private!
   integer, dimension(:), allocatable, public :: npchunks 
!   integer, dimension(:), allocatable, private :: npchunks 
!!XXgoldyXX: ^ this should be private
                                       ! number of chunks assigned to each process

   type lchunk
     integer  :: ncols                 ! number of vertical columns
     integer  :: cid                   ! global chunk index
     integer,  allocatable :: gcol(:)  ! global physics column indices
     real(r8), allocatable :: area(:)  ! column surface area (from dynamics)
     real(r8), allocatable :: wght(:)  ! column integration weight (from dynamics)
     real(r8) :: cost                  ! measured computational cost (seconds)
   end type lchunk

   integer, private :: nlchunks        ! local chunk count
   type (lchunk), dimension(:), allocatable, private :: lchunks  
                                       ! local chunks

   type knuhc
     integer  :: chunkid               ! chunk id
     integer  :: col                   ! column index in chunk
   end type knuhc

   type (knuhc), dimension(:), allocatable, public :: knuhcs !now beingh used in RRTMG radiation.F90
                                       ! map from global column indices
                                       ! to chunk'ed grid

! column mapping data structures
   type column_map
     integer  :: chunk                 ! global chunk index
     integer  :: ccol                  ! column ordering in chunk
   end type column_map

   integer, private :: nlcols           ! local column count
   type (column_map), dimension(:), allocatable, private :: pgcols
                                       ! ordered list of columns (for use in gather/scatter)
                                       ! NOTE: consistent with local ordering

! column remap data structures
   integer, dimension(:), allocatable, private :: gs_col_num
                                       ! number of columns scattered to each process in
                                       ! field_to_chunk scatter
   integer, dimension(:), allocatable, private :: gs_col_offset
                                       ! offset of columns (-1) in pgcols scattered to
                                       ! each process in field_to_chunk scatter

   integer, dimension(:), allocatable, private :: btofc_blk_num
                                       ! number of grid points scattered to each process in
                                       ! block_to_chunk alltoallv, and gathered from each
                                       ! process in chunk_to_block alltoallv

   integer, dimension(:), allocatable, private :: btofc_chk_num
                                       ! number of grid points gathered from each process in
                                       ! block_to_chunk alltoallv, and scattered to each
                                       ! process in chunk_to_block alltoallv

   type btofc_pters
     integer :: ncols                  ! number of columns in block
     integer :: nlvls                  ! number of levels in columns
     integer, dimension(:,:), pointer :: pter 
   end type btofc_pters
   type (btofc_pters), dimension(:), allocatable, private :: btofc_blk_offset
                                       ! offset in btoc send array (-1) where 
                                       ! (blockid, bcid, k) column should be packed in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob receive array (-1) from which
                                       ! (blockid, bcid, k) column should be unpacked in
                                       ! chunk_to_block alltoallv

   type (btofc_pters), dimension(:), allocatable, private :: btofc_chk_offset
                                       ! offset in btoc receive array (-1) from which
                                       ! (lcid, i, k) data should be unpacked in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob send array (-1) where
                                       ! (lcid, i, k) data should be packed in
                                       ! chunk_to_block alltoallv

! miscellaneous phys_grid data
   integer, private :: dp_coup_steps   ! number of swaps in transpose algorithm
   integer, dimension(:), private, allocatable :: dp_coup_proc
                                       ! swap partner in each step of 
                                       !  transpose algorithm
   logical :: physgrid_set = .false.   ! flag indicates physics grid has been set
   integer, private :: max_nproc_vsmp  ! maximum number of processes assigned to a
                                       !  single virtual SMP used to define physics 
                                       !  load balancing
   integer, private :: nproc_busy_d    ! number of processes active during the dynamics
                                       !  (assigned a dynamics block)

! Physics grid decomposition environment
   integer, private :: nlthreads       ! number of OpenMP threads available to this process
   integer, dimension(:), allocatable, private :: npthreads
                                       ! number of OpenMP threads available to each process
                                       !  (deallocated at end of phys_grid_init)

! Physics fields data structure (chunk) first dimension (pcols) option:
!  1 <= pcols_opt: pcols is set to pcols_opt
!  0 >= pcols_opt: calculate pcols based on lbal_opt, chunks_per_thread, number of
!       columns and threads in virtual SMP and relative costs per column (if provided), 
!       attempting to minimize wasted space and number of chunks, subject to:
!   0 <  pcols_max, then this is an upper bound on the calculated pcols.
!        Otherwise, pcols_max is ignored.
!   1 <  pcols_mult, then pcols is required to be a multiple of pcols_mult
!        (if pcols_max > 0 and pcols_mult <= pcols_max)
!        Otherwise, pcols_mult is ignored.
   integer, private, parameter :: def_pcols_opt  = 0              ! default 
   integer, private :: pcols_opt = def_pcols_opt

   integer, private, parameter :: min_pcols_max  = 1
   integer, private, parameter :: def_pcols_max  = -1             ! default 
   integer, private :: pcols_max = def_pcols_max

   integer, private, parameter :: min_pcols_mult = 1
   integer, private, parameter :: def_pcols_mult = 1              ! default 
   integer, private :: pcols_mult = def_pcols_mult

! Physics grid decomposition options:  
! -1: each chunk is a dynamics block
!  0: chunk definitions and assignments do not require interprocess comm.
!  1: chunk definitions and assignments do not require internode comm.
!  2: chunk definitions and assignments may require communication between all processes
!  3: chunk definitions and assignments only require communication with one other process
!  4: concatenated blocks, no load balancing, no interprocess communication
   integer, private, parameter :: min_lbal_opt = -1
   integer, private, parameter :: max_lbal_opt = 5
   integer, private, parameter :: def_lbal_opt = 2               ! default 
   integer, private :: lbal_opt = def_lbal_opt

! Physics grid load balancing options:  
!  0: assign columns to chunks as single columns, wrap mapped across chunks
!  1: use (day/night; north/south) twin algorithm to determine load-balanced pairs of 
!       columns and assign columns to chunks in pairs, wrap mapped
   integer, private, parameter :: min_twin_alg = 0
   integer, private, parameter :: max_twin_alg = 1
   integer, private, parameter :: def_twin_alg_lonlat = 1         ! default
   integer, private, parameter :: def_twin_alg_unstructured = 0
   integer, private :: twin_alg = def_twin_alg_lonlat

! Physics grid load balancing output options:  
!  T: write out both estimated (normalized) and actual (seconds) cost per chunk
!  F: do not write out costs
   logical, private, parameter :: def_output_chunk_costs = .true.
   logical, private :: output_chunk_costs = def_output_chunk_costs

! target number of chunks per thread
   integer, private, parameter :: min_chunks_per_thread = 1
   integer, private, parameter :: def_chunks_per_thread = &
                                    min_chunks_per_thread         ! default
   integer, private :: chunks_per_thread = def_chunks_per_thread

! Dynamics/physics transpose method for nonlocal load-balance:
! -1: use "0" if max_nproc_vsmp and nproc_busy_d are both > npes/2; otherwise use "1"
!  0: use mpi_alltoallv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
!  11-13: use mod_comm, choosing any of several methods internal to mod_comm.
!      The method within mod_comm (denoted mod_method) has possible values 0,1,2 and
!      is set according to mod_method = phys_alltoall - modmin_alltoall, where
!      modmin_alltoall is 11.
   integer, private, parameter :: min_alltoall = -1
   integer, private, parameter :: max_alltoall = 3
# if defined(MODCM_DP_TRANSPOSE)
   integer, private, parameter :: modmin_alltoall = 11
   integer, private, parameter :: modmax_alltoall = 13
# endif
   integer, private, parameter :: def_alltoall = 1                ! default
   integer, private :: phys_alltoall = def_alltoall

contains
!========================================================================
  integer function get_nlcols_p()
    get_nlcols_p = nlcols
  end function get_nlcols_p

  integer function get_clon_p_tot()
    get_clon_p_tot = clon_p_tot
  end function get_clon_p_tot
  integer function get_clat_p_tot()
    get_clat_p_tot = clat_p_tot
  end function get_clat_p_tot

  subroutine phys_grid_init( )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Physics mapping initialization routine:  
    ! 
    ! Method: 
    ! 
    ! Author: John Drake and Patrick Worley
    ! 
    !-----------------------------------------------------------------------
    use pmgrid,           only: plev
    use dycore,           only: dycore_is
    use dyn_grid,         only: get_block_bounds_d, get_block_gcol_d,     &
                                get_block_gcol_cnt_d, get_block_levels_d, &
                                get_block_lvl_cnt_d, get_block_owner_d,   &
                                get_gcol_block_d, get_gcol_block_cnt_d,   &
                                get_horiz_grid_dim_d, get_horiz_grid_d,   &
                                physgrid_copy_attributes_d
    use spmd_utils,       only: pair, ceil2
    use cam_grid_support, only: cam_grid_register, iMap
    use cam_grid_support, only: hcoord_len => max_hcoordname_len
    use cam_grid_support, only: horiz_coord_t, horiz_coord_create
    use cam_grid_support, only: cam_grid_attribute_copy
    !
    !------------------------------Arguments--------------------------------
    !
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: i, j, jb, k, p             ! loop indices
    integer :: pre_i                      ! earlier index in loop iteration
    integer :: clat_p_dex, clon_p_dex     ! indices into unique lat. and lon. arrays
    integer :: maxblksiz                  ! maximum number of columns in a dynamics block
    integer :: beg_dex, end_dex           ! index range
    integer :: cid, lcid                  ! global and local chunk ids
    integer :: max_ncols                  ! upper bound on number of columns in a block
    integer :: ncols                      ! number of columns in current chunk
    integer :: curgcol, curgcol_d         ! current global column index
    integer :: firstblock, lastblock      ! global block indices
    integer :: blksiz                     ! current block size
    integer :: glbcnt, curcnt             ! running grid point counts
    integer :: curp                       ! current process id
    integer :: block_cnt                  ! number of blocks containing data
    ! for a given vertical column
    integer :: numlvl                     ! number of vertical levels in block 
    integer :: levels(plev+1)             ! vertical level indices
    integer :: owner_d                    ! process owning given block column
    integer :: owner_p                    ! process owning given chunk column
    integer :: blockids(plev+1)           ! block indices
    integer :: bcids(plev+1)              ! block column indices
    real(r8), parameter :: deg2rad = SHR_CONST_PI/180.0
    real(r8), dimension(:), allocatable :: area_d     ! column surface area (from dynamics)
    real(r8), dimension(:), allocatable :: wght_d     ! column integration weight (from dynamics)
    integer,  dimension(:), allocatable :: pchunkid   ! chunk global ordering
    integer,  dimension(:), allocatable :: tmp_gcol   ! work array for global physics column indices

    ! permutation array used in physics column sorting;
    ! reused later as work space in (lbal_opt == -1) logic
    integer, dimension(:), allocatable :: cdex

    real(r8), dimension(:), allocatable :: clat_d   ! lat (radians) from dynamics columns
    real(r8), dimension(:), allocatable :: clon_d   ! lon (radians) from dynamics columns
    real(r8), dimension(:), allocatable :: lat_d    ! lat from dynamics columns
    real(r8), dimension(:), allocatable :: lon_d    ! lon from dynamics columns
    real(r8) :: clat_p_tmp
    real(r8) :: clon_p_tmp

    ! for calculation and output of physics decomposition statistics
    integer :: chunk_ncols                ! number of columns assigned to a chunk
    integer :: min_chunk_ncols            ! min number of columns assigned to a chunk
    integer :: max_chunk_ncols            ! max number of columns assigned to a chunk
    integer :: min_process_nthreads       ! min number of threads available to a process
    integer :: max_process_nthreads       ! max number of threads available to a process
    integer :: min_process_nchunks        ! min number of chunks assigned to a process
    integer :: max_process_nchunks        ! max number of chunks assigned to a process
    integer :: min_process_ncols          ! min number of columns assigned to a process
    integer :: max_process_ncols          ! max number of columns assigned to a process
    integer :: min_pcols                  ! min pcols over all processes
    integer :: max_pcols                  ! max pcols over all processes
    integer, dimension(:), allocatable :: pcols_proc
                                          ! pcols for all chunks assigned to each process
    integer, dimension(:), allocatable :: process_ncols ! number of columns per process
    integer, dimension(:), allocatable :: maxblksiz_proc
                                          ! maxblksiz for blocks assigned to each process

    ! Maps and values for physics grid
    real(r8),                   pointer :: lonvals(:)
    real(r8),                   pointer :: latvals(:)
    real(r8),               allocatable :: latdeg_p(:)
    real(r8),               allocatable :: londeg_p(:)
    integer(iMap),              pointer :: grid_map(:,:)
    integer(iMap),          allocatable :: coord_map(:)
    type(horiz_coord_t),        pointer :: lat_coord
    type(horiz_coord_t),        pointer :: lon_coord
    integer,                allocatable :: gcols(:)
    character(len=hcoord_len),  pointer :: copy_attributes(:)
    character(len=hcoord_len)           :: copy_gridname
    logical                             :: unstructured
    real(r8)                            :: lonmin, latmin

    integer :: nlthreads                  ! number of local OpenMP threads
#if ( defined _OPENMP )
    integer omp_get_max_threads
    external omp_get_max_threads
#endif

    nullify(lonvals)
    nullify(latvals)
    nullify(grid_map)
    nullify(lat_coord)
    nullify(lon_coord)

    call t_adj_detailf(-2)
    call t_startf("phys_grid_init")

    !-----------------------------------------------------------------------
    !
    ! Initialize physics grid, using dynamics grid
    ! a) column coordinates
    if (single_column .and. .not. iop_mode .and. dycore_is ('SE')) lbal_opt = -1
    call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
    if (single_column .and. .not. iop_mode .and. dycore_is('SE')) then
      ngcols = 1
    else
      ngcols = hdim1_d*hdim2_d
    endif
    allocate( clat_d(1:ngcols) )
    allocate( clon_d(1:ngcols) )
    allocate( lat_d(1:ngcols) )
    allocate( lon_d(1:ngcols) )
    allocate( cdex(1:ngcols) )
    clat_d = 100000.0_r8
    clon_d = 100000.0_r8
    if (single_column .and. dycore_is('SE')) then
      lat_d = scmlat
      lon_d = scmlon
      clat_d = scmlat * deg2rad
      clon_d = scmlon * deg2rad
    else
      call get_horiz_grid_d(ngcols, clat_d_out=clat_d, clon_d_out=clon_d, lat_d_out=lat_d, lon_d_out=lon_d)
    endif
    latmin = MINVAL(ABS(lat_d))
    lonmin = MINVAL(ABS(lon_d))

    ! Get estimated computational cost weight for each column (only from SE dynamics currently)
    allocate( cost_d (1:ngcols) )
    cost_d(:) = 1.0_r8
    use_cost_d = .false.
    if ((.not. single_column) .and. dycore_is('SE')) then
      call get_horiz_grid_d(ngcols, cost_d_out=cost_d)
      if (minval(cost_d) .ne. maxval(cost_d)) use_cost_d = .true.
    endif

!!XXgoldyXX: To do: replace collection above with local physics points

    ! count number of "real" column indices
    ngcols_p = 0
    do i=1,ngcols
       if (clon_d(i) < 100000.0_r8) then
          ngcols_p = ngcols_p + 1
       endif
    enddo

    ! sort over longitude and identify unique longitude coordinates
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clon_d,descend=.false.)
    clon_p_tmp = clon_d(cdex(1))
    clon_p_tot = 1

    do i=2,ngcols_p
       if (clon_d(cdex(i)) > clon_p_tmp) then
          clon_p_tot = clon_p_tot + 1
          clon_p_tmp = clon_d(cdex(i))
       endif
    enddo

    allocate( clon_p(1:clon_p_tot) )
    allocate( clon_p_cnt(1:clon_p_tot) )
    allocate( londeg_p(1:clon_p_tot) )

    pre_i = 1
    clon_p_tot = 1
    clon_p(1) = clon_d(cdex(1))
    londeg_p(1) = lon_d(cdex(1))
    do i=2,ngcols_p
       if (clon_d(cdex(i)) > clon_p(clon_p_tot)) then
          clon_p_cnt(clon_p_tot) = i-pre_i
          pre_i = i
          clon_p_tot = clon_p_tot + 1
          clon_p(clon_p_tot) = clon_d(cdex(i))
          londeg_p(clon_p_tot) = lon_d(cdex(i))
       endif
    enddo
    clon_p_cnt(clon_p_tot) = (ngcols_p+1)-pre_i

    ! sort over latitude and identify unique latitude coordinates
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clat_d,descend=.false.)
    clat_p_tmp = clat_d(cdex(1))
    clat_p_tot = 1
    do i=2,ngcols_p
       if (clat_d(cdex(i)) > clat_p_tmp) then
          clat_p_tot = clat_p_tot + 1
          clat_p_tmp = clat_d(cdex(i))
       endif
    enddo

    allocate( clat_p(1:clat_p_tot) )
    allocate( clat_p_cnt(1:clat_p_tot) )
    allocate( clat_p_idx(1:clat_p_tot) )
    allocate( latdeg_p(1:clat_p_tot) )

    pre_i = 1
    clat_p_tot = 1
    clat_p(1) = clat_d(cdex(1))
    latdeg_p(1) = lat_d(cdex(1))
    do i=2,ngcols_p
       if (clat_d(cdex(i)) > clat_p(clat_p_tot)) then
          clat_p_cnt(clat_p_tot) = i-pre_i
          pre_i = i
          clat_p_tot = clat_p_tot + 1
          clat_p(clat_p_tot) = clat_d(cdex(i))
          latdeg_p(clat_p_tot) = lat_d(cdex(i))
       endif
    enddo
    clat_p_cnt(clat_p_tot) = (ngcols_p+1)-pre_i

    clat_p_idx(1) = 1
    do j=2,clat_p_tot
       clat_p_idx(j) = clat_p_idx(j-1) + clat_p_cnt(j-1)
    enddo

    deallocate( lat_d )
    deallocate( lon_d )

    ! sort by longitude within latitudes
    end_dex = 0
    do j=1,clat_p_tot
       beg_dex = end_dex + 1
       end_dex = end_dex + clat_p_cnt(j)
       call IndexSort(cdex(beg_dex:end_dex),clon_d,descend=.false.)
    enddo

    ! Early clean-up, to minimize memory high water mark
    ! (not executing find_partner or find_twin)
    if (((twin_alg .ne. 1) .and. (lbal_opt .ne. 3)) .or. &
        (lbal_opt .eq. -1)) deallocate( clat_p_cnt)

    ! save "longitude within latitude" column ordering
    ! and determine mapping from unsorted global column index to 
    ! unique latitude/longitude indices
    allocate( lat_p(1:ngcols) )
    allocate( lon_p(1:ngcols) )
    allocate( dyn_to_latlon_gcol_map(1:ngcols) )
    if (lbal_opt .ne. -1) allocate( latlon_to_dyn_gcol_map(1:ngcols_p) )

    clat_p_dex = 1
    lat_p = -1
    dyn_to_latlon_gcol_map = -1
    do i=1,ngcols_p
       if (lbal_opt .ne. -1) latlon_to_dyn_gcol_map(i) = cdex(i)
       dyn_to_latlon_gcol_map(cdex(i)) = i

       do while ((clat_p(clat_p_dex) < clat_d(cdex(i))) .and. &
                 (clat_p_dex < clat_p_tot))
          clat_p_dex = clat_p_dex + 1
       enddo
       lat_p(cdex(i)) = clat_p_dex
    enddo

    ! sort by latitude within longitudes
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clon_d,descend=.false.)
    end_dex = 0
    do i=1,clon_p_tot
       beg_dex = end_dex + 1
       end_dex = end_dex + clon_p_cnt(i)
       call IndexSort(cdex(beg_dex:end_dex),clat_d,descend=.false.)
    enddo

    ! Early clean-up, to minimize memory high water mark
    ! (not executing find_twin)
    if ((twin_alg .ne. 1) .or. (lbal_opt .eq. -1)) deallocate( clon_p_cnt )

    ! save "latitude within longitude" column ordering
    ! (only need in find_twin)
    if ((twin_alg .eq. 1) .and. (lbal_opt .ne. -1)) &
       allocate( lonlat_to_dyn_gcol_map(1:ngcols_p) )

    clon_p_dex = 1
    lon_p = -1
    do i=1,ngcols_p
       if ((twin_alg .eq. 1) .and. (lbal_opt .ne. -1)) &
         lonlat_to_dyn_gcol_map(i) = cdex(i)
       do while ((clon_p(clon_p_dex) < clon_d(cdex(i))) .and. &
                 (clon_p_dex < clon_p_tot))
          clon_p_dex = clon_p_dex + 1
       enddo
       lon_p(cdex(i)) = clon_p_dex
    enddo

    ! Clean-up
    deallocate( clat_d )
    deallocate( clon_d )
    deallocate( cdex )

    !
    ! Determine number of threads per process, for use in create_chunks
    ! and in output of decomposition statistics
    !
    allocate( npthreads(0:npes-1) )
    npthreads(:) = 0

    nlthreads = 1
#if ( defined _OPENMP )
    nlthreads = OMP_GET_MAX_THREADS()
#endif
!
#if ( defined SPMD )
    call mpiallgatherint(nlthreads, 1, npthreads, 1, mpicom)
#else
    npthreads(0) = nlthreads
#endif

    !
    ! Determine block index bounds
    !
    call get_block_bounds_d(firstblock,lastblock)

    !
    ! Allocate storage to save number of chunks, pcols for all
    ! chunks, and columns assigned to each process during chunk
    ! creation and assignment.
    !
    allocate( npchunks(0:npes-1) )
    allocate( pcols_proc(0:npes-1) )
    allocate( gs_col_num(0:npes-1) )
    npchunks(:) = 0
    pcols_proc(:) = 0
    gs_col_num(:) = 0

    !
    ! Option -1: each dynamics block is a single chunk
    !          
    if (lbal_opt == -1) then

       ! Allocate storage to save maximum number of columns per block
       ! over blocks assigned to a given process
       allocate( maxblksiz_proc(0:npes-1) )

       !
       ! Calculate maximum block size for each process
       !
       if (single_column .and. .not. iop_mode .and. dycore_is('SE')) then
          maxblksiz_proc(:) = 1
       else
          maxblksiz_proc(:) = 0
          do j=firstblock,lastblock
             p = get_block_owner_d(j)
             blksiz = get_block_gcol_cnt_d(j)
             if (blksiz > maxblksiz_proc(p)) maxblksiz_proc(p) = blksiz
          enddo
       endif

       !
       ! Calculate pcols and check that pcols >= maxblksiz
       !
#ifdef PPCOLS
       ! use compile-time value
       pcols_proc(:) = pcols
#else
       if (pcols_opt > 0) then
          ! use runtime value
          pcols_proc(:) = pcols_opt
       else
          do p=0,npes-1
             ! calculated default is maxblksiz
             pcols_proc(p) = maxblksiz_proc(p)
             ! then increase (if necessary) to be a multiple of pcols_mult
             if (pcols_mult > 1) then
                pcols_proc(p) = pcols_mult*ceiling(real(pcols_proc(p),r8)/real(pcols_mult,r8))
             endif
             ! then decrease (if necessary) to be no greater than pcols_max
             if (pcols_max > 0) then
                pcols_proc(p) = min(pcols_proc(p), pcols_max)
             endif
          enddo
       endif
       pcols = pcols_proc(iam)
#endif
       max_pcols = maxval(pcols_proc)
       maxblksiz = maxval(maxblksiz_proc)

       if (pcols < maxblksiz_proc(iam)) then
	  write(iulog,*) 'pcols = ',pcols, ' maxblksiz=',maxblksiz_proc(iam)
          call endrun ('PHYS_GRID_INIT error: phys_loadbalance -1 specified but PCOLS < MAXBLKSIZ')
       endif

       !
       ! Determine total number of chunks
       !
       if (single_column .and. .not. iop_mode .and. dycore_is('SE')) then
         nchunks = 1
       else
	 nchunks = (lastblock-firstblock+1)
       endif

       !
       ! Set max virtual SMP node size
       !
       max_nproc_vsmp = 1

       !
       ! Allocate and initialize part of chunks data structure
       !
       allocate( cdex(1:maxblksiz) )
       allocate( chunks(1:nchunks) )
       do cid=1,nchunks
          allocate( chunks(cid)%gcol(max_pcols) )
       enddo
       chunks(:)%estcost = 0.0_r8

       do cid=1,nchunks
          ! get number of global column indices in block
          if (single_column .and. .not. iop_mode .and. dycore_is('SE')) then
            max_ncols = 1
          else
            max_ncols = get_block_gcol_cnt_d(cid+firstblock-1)
          endif
          ! fill cdex array with global indices from current block
          call get_block_gcol_d(cid+firstblock-1,max_ncols,cdex)

          ncols = 0
          do i=1,max_ncols
             ! check whether global index is for a column that dynamics
             ! intends to pass to the physics
             curgcol_d = cdex(i)
             if (dyn_to_latlon_gcol_map(curgcol_d) .ne. -1) then
                ! yes - then save the information
                ncols = ncols + 1
                chunks(cid)%gcol(ncols) = curgcol_d
                chunks(cid)%estcost = chunks(cid)%estcost + cost_d(curgcol_d)
             endif
          enddo
          chunks(cid)%ncols = ncols
       enddo

       ! Clean-up
       deallocate( cdex )
       deallocate( maxblksiz_proc )

       !
       ! Specify parallel decomposition 
       !
       do cid=1,nchunks
#if (defined SPMD)
          p = get_block_owner_d(cid+firstblock-1)
#else
          p = 0
#endif
          chunks(cid)%owner = p
          npchunks(p)       = npchunks(p) + 1
          gs_col_num(p)     = gs_col_num(p) + chunks(cid)%ncols
       enddo
       !
       ! Set flag indicating columns in physics and dynamics 
       ! decompositions reside on the same processes
       !
       local_dp_map = .true. 
       !
    else
       !
       ! Option == 0: split local blocks into chunks,
       !               while attempting to create load-balanced chunks.
       !               Does not work with vertically decomposed blocks.
       !               (default)
       ! Option == 1: split SMP-local blocks into chunks,
       !               while attempting to create load-balanced chunks.
       !               Does not work with vertically decomposed blocks.
       ! Option == 2: load balance chunks with respect to diurnal and
       !               seaonsal cycles and wth respect to latitude, 
       !               and assign chunks to processes
       !               in a way that attempts to minimize communication costs
       ! Option == 3: divide processes into pairs and split 
       !               blocks assigned to these pairs into 
       !               chunks, attempting to create load-balanced chunks.
       !               The process pairs are chosen to maximize load balancing
       !               opportunities.
       !               Does not work with vertically decomposed blocks.
       ! Option == 4: concatenate local blocks, then
       !               divide into chunks.
       !               Does not work with vertically decomposed blocks.
       ! Option == 5: split individual blocks into chunks,
       !               assigning columns using block ordering
       !
       !
       ! Allocate and initialize chunks data structure, then
       ! assign chunks to processes.
       !

       if  (twin_alg .eq. 1) then
          ! precompute clon_p_idx: index in lonlat ordering for first 
          ! occurrence of longitude corresponding to given latitude index,
          ! used in twin option in create_chunks
          allocate( clon_p_idx(1:clon_p_tot) )
          clon_p_idx(1) = 1
          do i=2,clon_p_tot
             clon_p_idx(i) = clon_p_idx(i-1) + clon_p_cnt(i-1)
          enddo
       endif

       call t_startf("create_chunks")
       call create_chunks(lbal_opt, chunks_per_thread, pcols_opt, &
                          pcols_max, pcols_mult, pcols_proc)
#ifndef PPCOLS
       pcols = pcols_proc(iam)
#endif
       call t_stopf("create_chunks")

       ! Early clean-up, to minimize memory high water mark
       !deallocate( latlon_to_dyn_gcol_map ) !do not deallocate as it is being used in RRTMG radiation.F90
       if  (twin_alg .eq. 1) deallocate( lonlat_to_dyn_gcol_map )
       if  (twin_alg .eq. 1) deallocate( clon_p_cnt )
       if  (twin_alg .eq. 1) deallocate( clon_p_idx )
       if ((twin_alg .eq. 1) .or. (lbal_opt .eq. 3)) deallocate( clat_p_cnt )

       !
       ! Determine whether dynamics and physics decompositions
       ! are colocated, not requiring any interprocess communication
       ! in the coupling.
       local_dp_map = .true.   
       do cid=1,nchunks
          do i=1,chunks(cid)%ncols
             curgcol_d = chunks(cid)%gcol(i)
             block_cnt = get_gcol_block_cnt_d(curgcol_d)
             call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
             do jb=1,block_cnt
                owner_d = get_block_owner_d(blockids(jb)) 
                if (owner_d .ne. chunks(cid)%owner) then
                   local_dp_map = .false.   
                endif
             enddo
          enddo
       enddo
    endif

    !
    ! Deallocate unneeded work space
    !
    deallocate( cost_d )

    !
    ! Resize gcol to dimension nlon to eliminate unused space, then
    ! allocate and set lat and lon member arrays
    !
    max_pcols = maxval(pcols_proc)
    allocate(tmp_gcol(max_pcols))
    do cid=1,nchunks

       ncols = chunks(cid)%ncols
       do i = 1, ncols
          tmp_gcol(i) = chunks(cid)%gcol(i)
       end do

       deallocate(chunks(cid)%gcol)
       allocate(chunks(cid)%gcol(ncols))

       allocate(chunks(cid)%lat(ncols))
       allocate(chunks(cid)%lon(ncols))

       do i = 1, ncols
          chunks(cid)%gcol(i) = tmp_gcol(i)
          chunks(cid)%lat(i)  = lat_p(tmp_gcol(i))
          chunks(cid)%lon(i)  = lon_p(tmp_gcol(i))
       end do

    enddo

    !
    ! Deallocate unneeded work space
    !
    deallocate( tmp_gcol )
    deallocate( lat_p )
    deallocate( lon_p )

    !
    ! Allocate and initialize data structures for gather/scatter
    !  
    allocate( pgcols(1:ngcols_p) )
    allocate( gs_col_offset(0:npes) )
    allocate( pchunkid(0:npes) )

    ! Initialize pchunkid and gs_col_offset by summing 
    ! number of chunks and columns per process, respectively
    pchunkid(0) = 0
    gs_col_offset(0) = 0
    do p=1,npes-1
       pchunkid(p)      = pchunkid(p-1)      + npchunks(p-1)
       gs_col_offset(p) = gs_col_offset(p-1) + gs_col_num(p-1)
    enddo
    
    ! Determine local ordering via "process id" bin sort
    do cid=1,nchunks
       p = chunks(cid)%owner
       pchunkid(p) = pchunkid(p) + 1

       chunks(cid)%lcid = pchunkid(p) + lastblock

       curgcol = gs_col_offset(p)
       do i=1,chunks(cid)%ncols
          curgcol = curgcol + 1
          pgcols(curgcol)%chunk = cid
          pgcols(curgcol)%ccol = i
       enddo
       gs_col_offset(p) = curgcol
    enddo

    ! Reinitialize pchunkid and gs_col_offset (for real)
    pchunkid(0) = 1
    gs_col_offset(0) = 1
    do p=1,npes-1
       pchunkid(p)      = pchunkid(p-1)      + npchunks(p-1)
       gs_col_offset(p) = gs_col_offset(p-1) + gs_col_num(p-1)
    enddo
    pchunkid(npes)      = pchunkid(npes-1)      + npchunks(npes-1)
    gs_col_offset(npes) = gs_col_offset(npes-1) + gs_col_num(npes-1)

    ! Save local information
    ! (Local chunk index range chosen so that it does not overlap 
    !  {begblock,...,endblock})
    ! 
    nlcols   = gs_col_num(iam)
    nlchunks = npchunks(iam)
    begchunk = pchunkid(iam)   + lastblock
    endchunk = pchunkid(iam+1) + lastblock - 1
    !
    allocate( lchunks(begchunk:endchunk) )
    do cid=1,nchunks
       if (chunks(cid)%owner == iam) then
          lcid = chunks(cid)%lcid
          lchunks(lcid)%ncols = chunks(cid)%ncols
          lchunks(lcid)%cid   = cid
          allocate( lchunks(lcid)%gcol(chunks(cid)%ncols) )
          allocate( lchunks(lcid)%area(chunks(cid)%ncols) )
          allocate( lchunks(lcid)%wght(chunks(cid)%ncols) )
          do i=1,chunks(cid)%ncols
             lchunks(lcid)%gcol(i) = chunks(cid)%gcol(i)
          enddo
       endif
    enddo
    lchunks(:)%cost = 0.0_r8

    deallocate( pchunkid )
    !deallocate( npchunks ) !do not deallocate as it is being used in RRTMG radiation.F90
    !
    !-----------------------------------------------------------------------
    !
    ! Initialize physics grid, using dynamics grid
    ! b) column area and integration weight

    allocate( area_d(1:ngcols) )
    allocate( wght_d(1:ngcols) )
    area_d = 0.0_r8
    wght_d = 0.0_r8

    if (single_column .and. .not. iop_mode .and. dycore_is('SE')) then
      area_d = 4.0_r8*pi
      wght_d = 4.0_r8*pi
    else
      call get_horiz_grid_d(ngcols, area_d_out=area_d, wght_d_out=wght_d)
    endif

    if ( abs(sum(area_d) - 4.0_r8*pi) > 1.e-10_r8 ) then
       write(iulog,*) ' ERROR: sum of areas on globe does not equal 4*pi'
       write(iulog,*) ' sum of areas = ', sum(area_d), sum(area_d)-4.0_r8*pi
       call endrun('phys_grid')
    end if

    if ( abs(sum(wght_d) - 4.0_r8*pi) > 1.e-10_r8 ) then
       write(iulog,*) ' ERROR: sum of integration weights on globe does not equal 4*pi'
       write(iulog,*) ' sum of weights = ', sum(wght_d), sum(wght_d)-4.0_r8*pi
       call endrun('phys_grid')
    end if

    do lcid=begchunk,endchunk
       do i=1,lchunks(lcid)%ncols
          lchunks(lcid)%area(i) = area_d(lchunks(lcid)%gcol(i))
          lchunks(lcid)%wght(i) = wght_d(lchunks(lcid)%gcol(i))
       enddo
    enddo

    deallocate( area_d )
    deallocate( wght_d )

    if (.not. local_dp_map) then

       !
       ! allocate and initialize data structures for transposes
       !  
       allocate( btofc_blk_num(0:npes-1) )
       btofc_blk_num = 0
       allocate( btofc_blk_offset(firstblock:lastblock) )
       do jb = firstblock,lastblock
          nullify( btofc_blk_offset(jb)%pter )
       enddo
       !
       glbcnt = 0
       curcnt = 0
       curp = 0
       do curgcol=1,ngcols_p
          cid = pgcols(curgcol)%chunk
          i   = pgcols(curgcol)%ccol
          owner_p   = chunks(cid)%owner
          do while (curp < owner_p)
             btofc_blk_num(curp) = curcnt
             curcnt = 0
             curp = curp + 1
          enddo
          curgcol_d = chunks(cid)%gcol(i)
          block_cnt = get_gcol_block_cnt_d(curgcol_d)
          call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
          do jb = 1,block_cnt
             owner_d = get_block_owner_d(blockids(jb))
             if (iam == owner_d) then
                if (.not. associated(btofc_blk_offset(blockids(jb))%pter)) then
                   blksiz = get_block_gcol_cnt_d(blockids(jb))
                   numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
                   btofc_blk_offset(blockids(jb))%ncols = blksiz
                   btofc_blk_offset(blockids(jb))%nlvls = numlvl
                   allocate( btofc_blk_offset(blockids(jb))%pter(blksiz,numlvl) )
                endif
                do k=1,btofc_blk_offset(blockids(jb))%nlvls
                   btofc_blk_offset(blockids(jb))%pter(bcids(jb),k) = glbcnt
                   curcnt = curcnt + 1
                   glbcnt = glbcnt + 1
                enddo
             endif
          enddo
       enddo
       
       btofc_blk_num(curp) = curcnt
       block_buf_nrecs = glbcnt
       !  
       allocate( btofc_chk_num(0:npes-1) )
       btofc_chk_num = 0
       allocate( btofc_chk_offset(begchunk:endchunk) )
       do lcid=begchunk,endchunk
          ncols = lchunks(lcid)%ncols
          btofc_chk_offset(lcid)%ncols = ncols
          btofc_chk_offset(lcid)%nlvls = pver+1
          allocate( btofc_chk_offset(lcid)%pter(ncols,pver+1) )
       enddo
       !
       curcnt = 0
       glbcnt = 0
       do p=0,npes-1
          do curgcol=gs_col_offset(iam),gs_col_offset(iam+1)-1
             cid  = pgcols(curgcol)%chunk
             owner_p  = chunks(cid)%owner
             if (iam == owner_p) then
                i    = pgcols(curgcol)%ccol
                lcid = chunks(cid)%lcid
                curgcol_d = chunks(cid)%gcol(i)
                block_cnt = get_gcol_block_cnt_d(curgcol_d)
                call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
                do jb = 1,block_cnt
                   owner_d = get_block_owner_d(blockids(jb))
                   if (p == owner_d) then
                      numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
                      call get_block_levels_d(blockids(jb),bcids(jb),numlvl,levels)
                      do k=1,numlvl
                         btofc_chk_offset(lcid)%pter(i,levels(k)+1) = glbcnt
                         curcnt = curcnt + 1
                         glbcnt = glbcnt + 1
                      enddo
                   endif
                enddo
             endif
          enddo
          btofc_chk_num(p) = curcnt
          curcnt = 0
       enddo
       chunk_buf_nrecs = glbcnt
       !
       ! Precompute swap partners and number of steps in point-to-point
       ! implementations of alltoall algorithm.
       ! First, determine number of swaps.
       !
       dp_coup_steps = 0
       do i=1,ceil2(npes)-1
          p = pair(npes,i,iam)
          if (p >= 0) then
             if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
                dp_coup_steps = dp_coup_steps + 1
             end if
          end if
       end do
       !
       ! Second, determine swap partners.
       !

       allocate( dp_coup_proc(dp_coup_steps) )
       dp_coup_steps = 0
       do i=1,ceil2(npes)-1
          p = pair(npes,i,iam)
          if (p >= 0) then
             if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
                dp_coup_steps = dp_coup_steps + 1
                dp_coup_proc(dp_coup_steps) = p
             end if
          end if
       end do
       !
    endif

    ! Final clean-up
    deallocate( gs_col_offset )
    ! (if eliminate get_lon_xxx, can also deallocate
    !  clat_p_idx, and grid_latlon?))

    ! Add physics-package grid to set of CAM grids
    ! physgrid always uses 'lat' and 'lon' as coordinate names; If dynamics
    !    grid is different, it will use different coordinate names

    ! First, create a map for the physics grid
    ! It's structure will depend on whether or not the physics grid is
    ! unstructured
    unstructured = dycore_is('UNSTRUCTURED')
    ! local chunks, so use pcols
    if (unstructured) then
      allocate(grid_map(3, pcols * (endchunk - begchunk + 1)))
    else
      allocate(grid_map(4, pcols * (endchunk - begchunk + 1)))
    end if
    grid_map = 0
    allocate( gcols(pcols) )
    allocate(latvals(size(grid_map, 2)))
    allocate(lonvals(size(grid_map, 2)))
    p = 0
    do lcid = begchunk, endchunk
      ncols = lchunks(lcid)%ncols
      call get_gcol_all_p(lcid, pcols, gcols)
      ! collect latvals and lonvals
      cid = lchunks(lcid)%cid
      do i = 1, chunks(cid)%ncols
        latvals(p + i) = latdeg_p(chunks(cid)%lat(i))
        lonvals(p + i) = londeg_p(chunks(cid)%lon(i))
      end do
      if (pcols > ncols) then
        ! Need to set these to detect unused columns
        latvals(p+ncols+1:p+pcols) = 1000.0_r8
        lonvals(p+ncols+1:p+pcols) = 1000.0_r8
      end if

      ! Set grid values for this chunk
      do i = 1, pcols
        p = p + 1
        grid_map(1, p) = i
        grid_map(2, p) = lcid
        if ((i <= ncols) .and. (gcols(i) > 0)) then
          if (unstructured) then
            grid_map(3, p) = gcols(i)
          else
            grid_map(3, p) = get_lon_p(lcid, i)
            grid_map(4, p) = get_lat_p(lcid, i)
          end if
        else
          if (i <= ncols) then
            call endrun("phys_grid_init: unmapped column")
          end if
        end if
      end do
    end do

    ! Note that if the dycore is using the same points as the physics grid,
    !      it will have already set up 'lat' and 'lon' axes for the physics grid
    !      However, these will be in the dynamics decomposition

    if (unstructured) then
      lon_coord => horiz_coord_create('lon', 'ncol', ngcols_p,                 &
                                      'longitude', 'degrees_east', 1,          &
                                      size(lonvals), lonvals, map=grid_map(3,:))
      lat_coord => horiz_coord_create('lat', 'ncol', ngcols_p,                 &
                                      'latitude', 'degrees_north', 1,          &
                                      size(latvals), latvals, map=grid_map(3,:))
    else

      allocate(coord_map(size(grid_map, 2)))

      ! Create a lon coord map which only writes from one of each unique lon
      where(latvals == latmin)
        coord_map(:) = grid_map(3, :)
      elsewhere
        coord_map(:) = 0_iMap
      end where
      lon_coord => horiz_coord_create('lon', 'lon', hdim1_d, 'longitude',     &
           'degrees_east', 1, size(lonvals), lonvals, map=coord_map)

      ! Create a lat coord map which only writes from one of each unique lat
      where(lonvals == lonmin)
        coord_map(:) = grid_map(4, :)
      elsewhere
        coord_map(:) = 0_iMap
      end where
      lat_coord => horiz_coord_create('lat', 'lat', hdim2_d, 'latitude',      &
           'degrees_north', 1, size(latvals), latvals, map=coord_map)

      deallocate(coord_map)

    end if ! unstructured
    
    call cam_grid_register('physgrid', phys_decomp, lat_coord, lon_coord,     &
         grid_map, unstruct=unstructured, block_indexed=.true.)
    ! Copy required attributes from the dynamics array
    nullify(copy_attributes)
    call physgrid_copy_attributes_d(copy_gridname, copy_attributes)
    do i = 1, size(copy_attributes)
      call cam_grid_attribute_copy(copy_gridname, 'physgrid', copy_attributes(i))
    end do
    ! Cleanup pointers (they belong to the grid now)
    nullify(grid_map)
    deallocate(latvals)
    nullify(latvals)
    deallocate(lonvals)
    nullify(lonvals)
    deallocate(gcols)
    ! Cleanup, we are responsible for copy attributes
    if (associated(copy_attributes)) then
      deallocate(copy_attributes)
      nullify(copy_attributes)
    end if

    !
    physgrid_set = .true.   ! Set flag indicating physics grid is now set
    !
    if (masterproc) then
!
! Determine number of threads per process
!
      nlthreads = 1
#if ( defined _OPENMP )
      nlthreads = OMP_GET_MAX_THREADS()
#endif
      allocate( process_ncols(0:npes-1) )
      process_ncols(:) = 0

      min_chunk_ncols = ngcols_p
      max_chunk_ncols = 0
      do cid = 1,nchunks
        chunk_ncols = chunks(cid)%ncols
        if (chunk_ncols < min_chunk_ncols) min_chunk_ncols = chunk_ncols
        if (chunk_ncols > max_chunk_ncols) max_chunk_ncols = chunk_ncols
        owner_p = chunks(cid)%owner
        process_ncols(owner_p) = process_ncols(owner_p) + chunk_ncols
      enddo

      min_process_nthreads = minval(npthreads)
      max_process_nthreads = maxval(npthreads)
      min_process_nchunks  = minval(npchunks)
      max_process_nchunks  = maxval(npchunks)
      min_process_ncols    = minval(process_ncols)
      max_process_ncols    = maxval(process_ncols)
      min_pcols            = minval(pcols_proc)
      deallocate(process_ncols)

      write(iulog,*) 'PHYS_GRID_INIT:  Using'
#ifdef PPCOLS
      write(iulog,*) '  PCOLS (compile-time parameter)=',pcols
#else
      write(iulog,*) '  PCOLS (masterproc)=            ',pcols
      write(iulog,*) '  phys_chnk_fdim=                ',pcols_opt
      if (pcols_opt <= 0) then
       write(iulog,*)'  phys_chnk_fdim_max=            ',pcols_max
       write(iulog,*)'  phys_chnk_fdim_mult=           ',pcols_mult
      endif
#endif
      write(iulog,*) '  phys_loadbalance=              ',lbal_opt
      write(iulog,*) '  phys_twin_algorithm=           ',twin_alg
      write(iulog,*) '  phys_alltoall=                 ',phys_alltoall
      write(iulog,*) '  chunks_per_thread=             ',chunks_per_thread
      write(iulog,*) '  num threads=                   ',nlthreads
      write(iulog,*) 'PHYS_GRID_INIT:  Decomposition Statistics:'
      write(iulog,*) '  total number of physics columns=   ',ngcols_p
      write(iulog,*) '  total number of chunks=            ',nchunks
      write(iulog,*) '  total number of physics processes= ',npes
      write(iulog,*) '  (min,max) # of threads per physics process: (',  &
                        min_process_nthreads,',',max_process_nthreads,')'
      write(iulog,*) '  (min,max) pcols per chunk:                  (',  &
                        min_pcols,',',max_pcols,')'
      write(iulog,*) '  (min,max) # of physics columns per chunk:   (',  &
                        min_chunk_ncols,',',max_chunk_ncols,')'
      write(iulog,*) '  (min,max) # of chunks per process:          (',  &
                        min_process_nchunks,',',max_process_nchunks,')'
      write(iulog,*) '  (min,max) # of physics columns per process: (',  &
                        min_process_ncols,',',max_process_ncols,')'
      write(iulog,*) ''
    endif

    ! Clean-up
    deallocate(pcols_proc)
    deallocate(npthreads)

    call t_stopf("phys_grid_init")
    call t_adj_detailf(+2)
    return
  end subroutine phys_grid_init

!========================================================================

subroutine phys_grid_find_col(lat, lon, owner, lcid, icol)

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: Find the global column closest to the point specified by lat
   !          and lon.  Return indices of owning process, local chunk, and 
   !          column.
   ! 
   ! Authors: Phil Rasch / Patrick Worley / B. Eaton
   ! 
   !-----------------------------------------------------------------------

   real(r8), intent(in) :: lat, lon    ! requested location in degrees
   integer, intent(out) :: owner       ! rank of chunk owner
   integer, intent(out) :: lcid      ! local chunk index
   integer, intent(out) :: icol        ! column index within the chunk

   ! local
   real(r8) dist2           ! the distance (in radians**2 from lat, lon)
   real(r8) distmin         ! the distance (in radians**2 from closest column)
   real(r8) latr, lonr      ! lat, lon (in radians) of requested location
   real(r8) clat, clon      ! lat, lon (in radians) of column being tested
   real(r8) const

   integer i
   integer cid
   !-----------------------------------------------------------------------

   ! Check that input lat and lon are in valid range
   if (lon < 0.0_r8 .or. lon >= 360._r8 .or. &
       lat < -90._r8 .or. lat > 90._r8) then
      if (masterproc) then
         write(iulog,*) &
            'phys_grid_find_col: ERROR: lon must satisfy 0.<=lon<360. and lat must satisfy -90<=lat<=90.'
         write(iulog,*) &
            'input lon=', lon, '  input lat=', lat
      endif
      call endrun('phys_grid_find_col: input ERROR')
   end if

   const = 180._r8/pi            ! degrees per radian
   latr = lat/const              ! to radians
   lonr = lon/const              ! to radians

   owner   = -999
   lcid  = -999
   icol    = -999
   distmin = 1.e10_r8

   ! scan all chunks for closest point to lat, lon
   do cid = 1, nchunks
      do i = 1, chunks(cid)%ncols
         clat = clat_p(chunks(cid)%lat(i))
         clon = clon_p(chunks(cid)%lon(i))
         dist2 = (clat-latr)**2 + (clon-lonr)**2
         if (dist2 < distmin ) then
            distmin = dist2
            owner = chunks(cid)%owner
            lcid = chunks(cid)%lcid
            icol = i
         endif
      enddo
   end do

end subroutine phys_grid_find_col

!========================================================================

subroutine phys_grid_find_cols(lat, lon, nclosest, owner, lcid, icol, distmin, mlats, mlons)

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: Find the global columns closest to the point specified by lat
   !          and lon.  Return indices of owning process, local chunk, and 
   !          column.
   ! 
   ! Authors: Phil Rasch / Patrick Worley / B. Eaton
   ! 
   !-----------------------------------------------------------------------
   use physconst,    only : rearth
   
   real(r8), intent(in) :: lat, lon            ! requested location in degrees
   integer, intent(in)  :: nclosest            ! number of closest points to find
   integer, intent(out) :: owner(nclosest)     ! rank of chunk owner
   integer, intent(out) :: lcid(nclosest)      ! local chunk index
   integer, intent(out) :: icol(nclosest)      ! column index within the chunk
   real(r8),intent(out) :: distmin(nclosest)   ! the distance (m) of the closest column(s)
   real(r8),intent(out) :: mlats(nclosest)     ! the latitude of the closest column(s)
   real(r8),intent(out) :: mlons(nclosest)     ! the longitude of the closest column(s)

   ! local
   real(r8) dist2           ! the distance (in radians**2 from lat, lon)
   real(r8) latr, lonr      ! lat, lon (in radians) of requested location
   real(r8) clat, clon      ! lat, lon (in radians) of column being tested
   real(r8) const

   integer i, j
   integer cid
   !-----------------------------------------------------------------------

   ! Check that input lat and lon are in valid range
   if (lon < 0.0_r8 .or. lon >= 360._r8 .or. &
       lat < -90._r8 .or. lat > 90._r8) then
      if (masterproc) then
         write(iulog,*) &
            'phys_grid_find_cols: ERROR: lon must satisfy 0.<=lon<360. and lat must satisfy -90<=lat<=90.'
         write(iulog,*) &
            'input lon=', lon, '  input lat=', lat
      endif
      call endrun('phys_grid_find_cols: input ERROR')
   end if

   const = 180._r8/pi            ! degrees per radian
   latr = lat/const              ! to radians
   lonr = lon/const              ! to radians

   owner(:)   = -999
   lcid(:)    = -999
   icol(:)    = -999
   mlats(:)   = -999
   mlons(:)   = -999
   distmin(:) = 1.e10_r8

   ! scan all chunks for closest point to lat, lon
   do cid = 1, nchunks
      do i = 1, chunks(cid)%ncols
         clat = clat_p(chunks(cid)%lat(i))
         clon = clon_p(chunks(cid)%lon(i))
         dist2 = acos(sin(latr) * sin(clat) + cos(latr) * cos(clat) * cos(clon - lonr)) * rearth       
         
         do j = nclosest, 1, -1
            if (dist2 < distmin(j)) then
            
               if (j < nclosest) then
                 distmin(j+1) = distmin(j)
                 owner(j+1)   = owner(j)
                 lcid(j+1)    = lcid(j)
                 icol(j+1)    = icol(j)
                 mlats(j+1)   = mlats(j)
                 mlons(j+1)    = mlons(j)
               end if
             
               distmin(j) = dist2
               owner(j)   = chunks(cid)%owner
               lcid(j)    = chunks(cid)%lcid
               icol(j)    = i
               mlats(j)   = clat * const
               mlons(j)   = clon * const
            else
               exit
            end if
         enddo
      enddo
   end do
   
end subroutine phys_grid_find_cols
!
!========================================================================

logical function phys_grid_initialized ()
!----------------------------------------------------------------------- 
! 
! Purpose: Identify whether phys_grid has been called yet or not
! 
! Method: Return physgrid_set
! 
! Author: Pat Worley
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   phys_grid_initialized = physgrid_set
!
   return
   end function phys_grid_initialized

!
!========================================================================
!
   subroutine phys_grid_defaultopts(phys_chnk_fdim_out, &
                                    phys_chnk_fdim_max_out, &
                                    phys_chnk_fdim_mult_out, &
                                    phys_loadbalance_out, &
                                    phys_twin_algorithm_out, &
                                    phys_alltoall_out, &
                                    phys_chnk_per_thd_out, &
                                    phys_chnk_cost_write_out )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
   use dycore, only: dycore_is
!------------------------------Arguments--------------------------------
     ! physics data structures declared first dimension option
     integer, intent(out), optional :: phys_chnk_fdim_out
     ! physics data structures declared first dimension upper bound
     integer, intent(out), optional :: phys_chnk_fdim_max_out
     ! physics data structures declared first dimension required factor
     integer, intent(out), optional :: phys_chnk_fdim_mult_out
     ! physics load balancing option
     integer, intent(out), optional :: phys_loadbalance_out
     ! algorithm to use when determining column pairs to assign to chunks
     integer, intent(out), optional :: phys_twin_algorithm_out
     ! alltoall option
     integer, intent(out), optional :: phys_alltoall_out
     ! number of chunks per thread
     integer, intent(out), optional :: phys_chnk_per_thd_out
     ! flag whether to write out estimated and actual cost per chunk
     logical, intent(out), optional :: phys_chnk_cost_write_out
!-----------------------------------------------------------------------
     if ( present(phys_chnk_fdim_out) ) then
       phys_chnk_fdim_out = def_pcols_opt
     endif
     if ( present(phys_chnk_fdim_max_out) ) then
       phys_chnk_fdim_max_out = def_pcols_max
     endif
     if ( present(phys_chnk_fdim_mult_out) ) then
       phys_chnk_fdim_mult_out = def_pcols_mult
     endif
     if ( present(phys_loadbalance_out) ) then
       phys_loadbalance_out = def_lbal_opt
     endif
     if ( present(phys_twin_algorithm_out) ) then
       if (dycore_is('UNSTRUCTURED')) then
          phys_twin_algorithm_out = def_twin_alg_unstructured
       else
          phys_twin_algorithm_out = def_twin_alg_lonlat
       endif
     endif
     if ( present(phys_alltoall_out) ) then
       phys_alltoall_out = def_alltoall
     endif
     if ( present(phys_chnk_per_thd_out) ) then
       phys_chnk_per_thd_out = def_chunks_per_thread
     endif
     if ( present(phys_chnk_cost_write_out) ) then
       phys_chnk_cost_write_out = def_output_chunk_costs
     endif
   end subroutine phys_grid_defaultopts
!
!========================================================================
!
   subroutine phys_grid_setopts(phys_chnk_fdim_in, &
                                phys_chnk_fdim_max_in, &
                                phys_chnk_fdim_mult_in, &
                                phys_loadbalance_in, &
                                phys_twin_algorithm_in, &
                                phys_alltoall_in,    &
                                phys_chnk_per_thd_in, &
                                phys_chnk_cost_write_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
   use spmd_utils, only: phys_mirror_decomp_req
#if defined(MODCM_DP_TRANSPOSE)
   use mod_comm, only: phys_transpose_mod
#endif
!------------------------------Arguments--------------------------------
     ! physics data structures declared first dimension option
     integer, intent(in), optional :: phys_chnk_fdim_in
     ! physics data structures declared first dimension upper bound
     integer, intent(in), optional :: phys_chnk_fdim_max_in
     ! physics data structures declared first dimension required factor
     integer, intent(in), optional :: phys_chnk_fdim_mult_in
     ! physics load balancing option
     integer, intent(in), optional :: phys_loadbalance_in
     ! option to use load balanced column pairs
     integer, intent(in), optional :: phys_twin_algorithm_in
     ! alltoall option
     integer, intent(in), optional :: phys_alltoall_in
     ! number of chunks per thread
     integer, intent(in), optional :: phys_chnk_per_thd_in
     ! flag whether to write out estimated and actual cost per chunk
     logical, intent(in), optional :: phys_chnk_cost_write_in
!-----------------------------------------------------------------------
#ifdef PPCOLS
     if ( present(phys_chnk_fdim_in) ) then
        if (phys_chnk_fdim_in /= pcols) then
           if (masterproc) then
              write(iulog,*)                                     &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_chnk_fdim=',  &
                 phys_chnk_fdim_in,                              &
                 '  differs from compile-time PCOLS parameter=', &
                 pcols,                                          &
                 '  .'
              write(iulog,*)                                     &
                 '  Must compile without -DPPCOLS to enable runtime'
              write(iulog,*)                                     &
                 '  option. Ignoring and using PCOLS parameter.'
           endif
        endif
     endif
#else
     if ( present(phys_chnk_fdim_in) ) then
        pcols_opt = phys_chnk_fdim_in
     endif
!
     if ( present(phys_chnk_fdim_max_in) ) then
        pcols_max = phys_chnk_fdim_max_in
        if (pcols_opt <= 0) then
           if (pcols_max < min_pcols_max) then
              if (masterproc) then
                 write(iulog,*)                                            &
                    'PHYS_GRID_SETOPTS:  ERROR:  phys_chnk_fdim_max=',     &
                    phys_chnk_fdim_max_in,                                 &
                    '  is out of range.  It must be at least as large as ',&
                    min_pcols_max,                                         &
                    '  . Ignoring.'
              endif
           endif
        endif
     endif
!
     if ( present(phys_chnk_fdim_in) ) then
        pcols_mult = phys_chnk_fdim_mult_in 
        if (pcols_opt <= 0) then
           if (pcols_mult < min_pcols_mult) then
              if (masterproc) then
                 write(iulog,*)                                            &
                    'PHYS_GRID_SETOPTS:  ERROR:  phys_chnk_fdim_mult=',    &
                    phys_chnk_fdim_mult_in,                                &
                    '  is out of range.  It must be at least as large as ',&
                    min_pcols_mult,                                        &
                    '  . Ignoring.'
              endif
           endif
        endif
     endif
#endif
!
     if ( present(phys_loadbalance_in) ) then
        lbal_opt = phys_loadbalance_in
        if ((lbal_opt < min_lbal_opt).or.(lbal_opt > max_lbal_opt)) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_loadbalance=', &
                 phys_loadbalance_in,                             &
                 '  is out of range.  It must be between ',       &
                 min_lbal_opt,' and ',max_lbal_opt
           endif
           call endrun
        endif
        if (lbal_opt .eq. 3) then
           phys_mirror_decomp_req = .true.
        else
           phys_mirror_decomp_req = .false.
        endif
     endif
!
     if ( present(phys_twin_algorithm_in) ) then
        twin_alg = phys_twin_algorithm_in
        if ((twin_alg < min_twin_alg).or.(twin_alg > max_twin_alg)) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_twin_algorithm=', &
                 phys_twin_algorithm_in,                             &
                 '  is out of range.  It must be between ',       &
                 min_twin_alg,' and ',max_twin_alg
           endif
           call endrun
        endif
     endif
!
     if ( present(phys_alltoall_in) ) then
        phys_alltoall = phys_alltoall_in
        if (((phys_alltoall .lt. min_alltoall) .or.    &
             (phys_alltoall .gt. max_alltoall))        &
# if defined(MODCM_DP_TRANSPOSE)
           .and.                                       &
            ((phys_alltoall .lt. modmin_alltoall) .or. &
             (phys_alltoall .gt. modmax_alltoall))     &
# endif
           ) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SET_OPTS:  ERROR:  phys_alltoall=',   &
                  phys_alltoall_in,                               &
                  '  is out of range.  It must be between ',      &
                  min_alltoall,' and ',max_alltoall
           endif
           call endrun
        endif
#if defined(SPMD)
# if defined(MODCM_DP_TRANSPOSE)
        phys_transpose_mod = phys_alltoall
# endif
#endif
     endif
!
     if ( present(phys_chnk_per_thd_in) ) then
        chunks_per_thread = phys_chnk_per_thd_in
        if (chunks_per_thread < min_chunks_per_thread) then
           if (masterproc) then
              write(iulog,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_chnk_per_thd=',&
                 phys_chnk_per_thd_in,                            &
                 ' is too small.  It must not be smaller than ',  &
                 min_chunks_per_thread
           endif
           call endrun
        endif
     endif
!
     if ( present(phys_chnk_cost_write_in) ) then
        output_chunk_costs = phys_chnk_cost_write_in
     endif
   end subroutine phys_grid_setopts
!
!========================================================================
!
   subroutine get_chunk_indices_p(index_beg, index_end)
!----------------------------------------------------------------------- 
! 
! Purpose: Return range of indices for local chunks
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(out) :: index_beg  ! first index used for local chunks
   integer, intent(out) :: index_end  ! last index used for local chunks
!-----------------------------------------------------------------------

   index_beg = begchunk
   index_end = endchunk

   return
   end subroutine get_chunk_indices_p
!
!========================================================================
!
   subroutine get_gcol_all_p(lcid, latdim, gcols)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all global column indices for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
     integer, intent(in)  :: lcid        ! local chunk id
     integer, intent(in)  :: latdim      ! declared size of output array

     integer, intent(out) :: gcols(:)    ! array of global latitude indices
!---------------------------Local workspace-----------------------------
     integer :: i                        ! loop index
     
!-----------------------------------------------------------------------
     gcols=-1
     do i=1,lchunks(lcid)%ncols
        gcols(i) = lchunks(lcid)%gcol(i)
     enddo
     return
   end subroutine get_gcol_all_p

!
!========================================================================
!
   integer function get_gcol_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global physics column index for chunk column
! 
! Method: 
! 
! Author: Jim Edwards / Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!-----------------------------------------------------------------------
   get_gcol_p = lchunks(lcid)%gcol(col)
   
   return
   end function get_gcol_p

!
!========================================================================

   subroutine get_gcol_vec_p(lcid, lth, cols, gcols)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global physics column indices for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid

!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   integer, intent(out) :: gcols(lth)    ! array of global physics 
                                         !  columns indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index

!-----------------------------------------------------------------------
   do i=1,lth
     gcols(i) = lchunks(lcid)%gcol(cols(i))
   enddo

   return
   end subroutine get_gcol_vec_p

!
!========================================================================
!
   integer function get_ncols_p(lcid)
!----------------------------------------------------------------------- 
! 
! Purpose: Return number of columns in chunk given the local chunk id.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid      ! local chunk id

!---------------------------Local workspace-----------------------------
   integer              :: cid       ! global chunk id

!-----------------------------------------------------------------------
   get_ncols_p = lchunks(lcid)%ncols

   return
   end function get_ncols_p
!
!========================================================================
!
   subroutine get_lat_all_p(lcid, latdim, lats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all global latitude indices for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: latdim        ! declared size of output array

   integer, intent(out) :: lats(latdim)  ! array of global latitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     lats(i) = chunks(cid)%lat(i)
   enddo

   return
   end subroutine get_lat_all_p
!
!========================================================================

   subroutine get_lat_vec_p(lcid, lth, cols, lats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global latitude indices for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid

!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   integer, intent(out) :: lats(lth)     ! array of global latitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,lth
     lats(i) = chunks(cid)%lat(cols(i))
   enddo

   return
   end subroutine get_lat_vec_p
!
!========================================================================

   integer function get_lat_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global latitude index for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   get_lat_p = chunks(cid)%lat(col)

   return
   end function get_lat_p
!
!========================================================================
!
   subroutine get_lon_all_p(lcid, londim, lons)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Was: Return all global longitude indices for chunk
!  Now: Return all longitude offsets (+1) for chunk. These are offsets
!       in ordered list of global columns from first
!       column with given latitude to column with given latitude
!       and longitude. This corresponds to the usual longitude indices
!       for full and reduced lon/lat grids.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: londim        ! declared size of output array

   integer, intent(out) :: lons(londim)  ! array of global longitude 
                                         !  indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: lat                        ! latitude index
   integer :: cid                        ! global chunk id
   integer :: gcol                       ! global column id in latlon 
                                         !  ordering

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     lat  = chunks(cid)%lat(i)
     gcol = dyn_to_latlon_gcol_map(chunks(cid)%gcol(i))
     lons(i) = (gcol - clat_p_idx(lat)) + 1
   enddo

   return
   end subroutine get_lon_all_p
!
!========================================================================

   subroutine get_lon_vec_p(lcid, lth, cols, lons)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Was: Return global longitude indices for set of chunk columns.
!  Now: Return longitude offsets (+1) for set of chunk columns. 
!       These are offsets in ordered list of global columns from first
!       column with given latitude to column with given latitude
!       and longitude. This corresponds to the usual longitude indices
!       for full and reduced lon/lat grids.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   integer, intent(out) :: lons(lth)     ! array of global longitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: lat                        ! latitude index
   integer :: cid                        ! global chunk id
   integer :: gcol                       ! global column id in latlon 
                                         !  ordering

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,lth
     lat = chunks(cid)%lat(cols(i))
     gcol = dyn_to_latlon_gcol_map(chunks(cid)%gcol(i))
     lons(i) = (gcol - clat_p_idx(lat)) + 1
   enddo

   return
   end subroutine get_lon_vec_p
!
!========================================================================

   integer function get_lon_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Was: Return global longitude index for chunk column.
!  Now: Return longitude offset (+1) for chunk column. This is the 
!       offset in ordered list of global columns from first
!       column with given latitude to column with given latitude
!       and longitude. This corresponds to the usual longitude index
!       for full and reduced lon/lat grids.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: cid                        ! global chunk id
   integer :: lat                        ! latitude index
   integer :: gcol                       ! global column id in latlon 
                                         !  ordering

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   lat = chunks(cid)%lat(col)
   gcol = dyn_to_latlon_gcol_map(chunks(cid)%gcol(col))
   get_lon_p = (gcol - clat_p_idx(lat)) + 1

   return
   end function get_lon_p
!
!========================================================================
!
   subroutine get_rlat_all_p(lcid, rlatdim, rlats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all latitudes (in radians) for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid           ! local chunk id
   integer, intent(in)  :: rlatdim        ! declared size of output array

   real(r8), intent(out) :: rlats(rlatdim)! array of latitudes

!---------------------------Local workspace-----------------------------
   integer :: i                           ! loop index
   integer :: cid                         ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     rlats(i) = clat_p(chunks(cid)%lat(i))
   enddo

   return
   end subroutine get_rlat_all_p
!
!========================================================================
!
   subroutine get_area_all_p(lcid, rdim, area)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all areas for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: rdim          ! declared size of output array

   real(r8), intent(out) :: area(rdim)   ! array of areas

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index

!-----------------------------------------------------------------------
   do i=1,lchunks(lcid)%ncols
     area(i) = lchunks(lcid)%area(i)
   enddo

   return
   end subroutine get_area_all_p
!
!========================================================================
!
   real(r8) function get_area_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return area for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!-----------------------------------------------------------------------
   get_area_p = lchunks(lcid)%area(col)

   return
   end function get_area_p
!
!========================================================================
!
   subroutine get_wght_all_p(lcid, rdim, wght)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all integration weights for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: rdim          ! declared size of output array

   real(r8), intent(out) :: wght(rdim)   ! array of integration weights

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index

!-----------------------------------------------------------------------
   do i=1,lchunks(lcid)%ncols
     wght(i) = lchunks(lcid)%wght(i)
   enddo

   return
   end subroutine get_wght_all_p
!
!========================================================================
!
   real(r8) function get_wght_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return integration weight for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!-----------------------------------------------------------------------
   get_wght_p = lchunks(lcid)%wght(col)

   return
   end function get_wght_p
!
!========================================================================
!
   subroutine get_rlat_vec_p(lcid, lth, cols, rlats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return latitudes (in radians) for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   real(r8), intent(out) :: rlats(lth)   ! array of latitudes

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,lth
     rlats(i) = clat_p(chunks(cid)%lat(cols(i)))
   enddo

   return
   end subroutine get_rlat_vec_p
!
!========================================================================

   real(r8) function get_rlat_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return latitude (in radians) for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   get_rlat_p = clat_p(chunks(cid)%lat(col))

   return
   end function get_rlat_p
!
!========================================================================
!
   subroutine get_rlon_all_p(lcid, rlondim, rlons)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all longitudes (in radians) for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid           ! local chunk id
   integer, intent(in)  :: rlondim        ! declared size of output array

   real(r8), intent(out) :: rlons(rlondim)! array of longitudes

!---------------------------Local workspace-----------------------------
   integer :: i                           ! loop index
   integer :: cid                         ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,chunks(cid)%ncols
     rlons(i) = clon_p(chunks(cid)%lon(i))
   enddo

   return
   end subroutine get_rlon_all_p
!
!========================================================================

   subroutine get_rlon_vec_p(lcid, lth, cols, rlons)
!----------------------------------------------------------------------- 
! 
! Purpose: Return longitudes (in radians) for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid         ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   real(r8), intent(out) :: rlons(lth)   ! array of longitudes

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   do i=1,lth
     rlons(i) = clon_p(chunks(cid)%lon(cols(i)))
   enddo

   return
   end subroutine get_rlon_vec_p
!
!========================================================================

   real(r8) function get_rlon_p(lcid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return longitude (in radians) for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: cid                        ! global chunk id

!-----------------------------------------------------------------------
   cid = lchunks(lcid)%cid
   get_rlon_p = clon_p(chunks(cid)%lon(col))

   return
   end function get_rlon_p
!
!========================================================================

   subroutine print_cost_p()
!----------------------------------------------------------------------- 
! 
! Purpose: Print walltime cost for all chunks.
! 
! Method: 
! 
! Author: P. Worley
! 
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
   integer  :: unitn                  ! file unit number
   integer  :: signal                 ! handshake variable
   integer  :: lcid                   ! local chunk id
   integer  :: cid                    ! global chunk id
   integer  :: owner                  ! chunk owner
   integer  :: ncols                  ! number of columns in chunk
   integer  :: ierr                   ! error return

   real(r8) :: cost_lsum, cost_gsum   ! local and global sums of 
                                      !  chunk costs
   real(r8) :: avg_cost               ! average chunk cost
   real(r8) :: avg_estcost            ! average estimated chunk cost
   real(r8) :: cost                   ! chunk cost
   real(r8) :: norm_cost              ! normalized chunk cost
   real(r8) :: norm_est_cost          ! normalized estimated chunk cost

   character(len=*), parameter :: fname = 'atm_chunk_costs.txt'

!-----------------------------------------------------------------------
   if (output_chunk_costs) then
      unitn = getunit()

      ! Calculate normalized measured cost
      cost_gsum = 0.0_r8
      cost_lsum = 0.0_r8
      do lcid = begchunk, endchunk
         cost_lsum = cost_lsum + lchunks(lcid)%cost
      enddo
#if ( defined SPMD )
      call MPI_ALLREDUCE(cost_lsum, cost_gsum, 1, MPI_REAL8, MPI_SUM, &
                         mpicom, ierr)
#else
      cost_gsum = cost_lsum
#endif
      if (abs(cost_gsum) > 0.000001_r8) then
         avg_cost = cost_gsum/nchunks
      else
         avg_cost = 1.0_r8
      endif

      ! Calculate normalized estimated cost
      cost_gsum = 0.0_r8
      do cid = 1,nchunks
         cost_gsum = cost_gsum + chunks(cid)%estcost
      enddo
      if (abs(cost_gsum) > 0.000001_r8) then
         avg_estcost = cost_gsum/nchunks
      else
         avg_estcost = 1.0_r8
      endif

      ! Take turns writing to fname.
      if (iam == 0) then
         open( unitn, file=trim(fname), status='REPLACE', &
                      form='FORMATTED', access='SEQUENTIAL' )

         write(unitn,'(a)') "ATM CHUNK COST"
         write(unitn,'(a)') &
            " owner   lcid    cid  pcols  ncols  estcost (norm)  cost (norm)  cost (seconds)"

         signal = 1
#if ( defined SPMD )
      else
         call mpirecv(signal, 1, mpiint, iam-1, iam, mpicom) 
         open( unitn, file=trim(fname), status='OLD', &
                      form='FORMATTED', access='SEQUENTIAL', position='APPEND' )
#endif
      endif

      do lcid = begchunk, endchunk
         cid   = lchunks(lcid)%cid
         owner = chunks(cid)%owner
         ncols = lchunks(lcid)%ncols
         cost  = lchunks(lcid)%cost
         norm_cost = cost/avg_cost
         norm_est_cost = chunks(cid)%estcost/avg_estcost
         write(unitn,'(i6,1x,i6,1x,i6,1x,i6,1x,i6,6x,f10.3,3x,f10.3,2x,e14.3)') &
               owner, lcid, cid, pcols, ncols, norm_est_cost, norm_cost, cost
      enddo

      close(unitn)

#if ( defined SPMD )
      if (iam+1 < npes) then
         call mpisend(signal, 1, mpiint, iam+1, iam+1, mpicom)
      endif
#endif

      call freeunit(unitn)
   endif

   return
   end subroutine print_cost_p
!
!========================================================================

   subroutine update_cost_p(lcid, cost_increment)
!----------------------------------------------------------------------- 
! 
! Purpose: Add cost_increment to walltime cost for chunk given the local
!          chunk id.
! 
! Method: 
! 
! Author: P. Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lcid          ! local chunk id
   real(r8), intent(in) :: cost_increment! increment to running walltime

!-----------------------------------------------------------------------
   lchunks(lcid)%cost = lchunks(lcid)%cost + cost_increment

   return
   end subroutine update_cost_p
!
!========================================================================
!
!  integer function get_gcol_owner_p(gcol)
!----------------------------------------------------------------------- 
! 
! Purpose: Return owner of physics column with indicate index
! 
! Method: 
! 
! Author: P. Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!  integer, intent(in)  :: gcol     ! physics column index
!
!-----------------------------------------------------------------------
!
!  get_gcol_owner_p = chunks(knuhcs(gcol)%chunkid)%owner
!
!  return
!  end function get_gcol_owner_p
!
!========================================================================

!  subroutine buff_to_chunk(fdim,mdim,lbuff,localchunks)
!-----------------------------------------------------------------------
!
! Purpose: Copy from local buffer 
!          to local chunk data structure.
!          Needed for cpl6.
!
! Method:
!
! Author: Pat Worley and Robert Jacob
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!  integer, intent(in) :: fdim      ! declared length of first lbuff dimension
!  integer, intent(in) :: mdim      ! declared length of middle lbuff dimension
!  real(r8), intent(in) :: lbuff(fdim, mdim) ! local lon/lat buffer
!
!  real(r8), intent(out):: localchunks(pcols,mdim,begchunk:endchunk) ! local chunks
!
!
!---------------------------Local workspace-----------------------------
!  integer :: i,j,m,n                      ! loop indices
!
!  integer, save :: numcols = 0
!  integer, allocatable, save :: columnid(:), chunkid(:)
!-----------------------------------------------------------------------
!
!  if (numcols .eq. 0) then
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!           endif
!        endif
!     enddo
!     allocate(columnid(1:n))
!     allocate(chunkid(1:n))
!
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!              columnid(n) = knuhcs(i)%col
!              chunkid(n)  = chunks(knuhcs(i)%chunkid)%lcid
!           endif
!        endif
!     end do
!
!     numcols = n
!  endif
!
!  if (numcols .gt. fdim) call endrun('buff_to_chunk')
!  do m=1,mdim
!#ifdef CPRCRAY
!!dir$ concurrent
!!dir$ prefervector, preferstream
!#endif
!     do n = 1, numcols
!        localchunks(columnid(n),m,chunkid(n)) = lbuff(n,m)
!     end do
!  end do
!
!  return
!  end subroutine buff_to_chunk
!
!========================================================================

   subroutine scatter_field_to_chunk(fdim,mdim,ldim, &
                                     hdim1d,globalfield,localchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute field
!          to decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 

!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: hdim1d    ! declared first horizontal index 
                                    ! dimension
   real(r8), intent(in) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
                                    ! global field

   real(r8), intent(out):: localchunks(fdim,pcols,mdim, &
                                       begchunk:endchunk,ldim) 
                                    ! local chunks

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local column index
   integer :: gcol                       ! global column index
   integer :: h1                         ! first horizontal dimension index
   integer :: h2                         ! second horizontal dimension index

#if ( defined SPMD )
   real(r8) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be scattered
   real(r8) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of scattered
                                         !  vector
   integer :: displs(0:npes-1)           ! scatter displacements
   integer :: sndcnts(0:npes-1)          ! scatter send counts
   integer :: recvcnt                    ! scatter receive count
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
   if (hdim1d < hdim1_d) then
      write(iulog,*) __FILE__,__LINE__,hdim1d,hdim1_d
      call endrun ('SCATTER_FIELD_TO_CHUNK error: hdim1d < hdim1_d')
   endif
   localchunks(:,:,:,:,:) = 0
#if ( defined SPMD )
   displs(0) = 0
   sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + sndcnts(p-1)
     sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   recvcnt = fdim*mdim*ldim*nlcols

   if (masterproc) then

! copy field into global (process-ordered) chunked data structure

      do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
         do i=1,ngcols_p
            cid  = pgcols(i)%chunk
            lid  = pgcols(i)%ccol
            gcol = chunks(cid)%gcol(lid)
            h2   = (gcol-1)/hdim1_d + 1
            h1   = mod((gcol-1),hdim1_d) + 1
            do m=1,mdim
               do f=1,fdim
                  gfield_p(f,m,l,i) = &
                     globalfield(f, h1, m, h2, l)
               end do
            end do
         end do
      end do
   endif

! scatter to other processes
! (pgcols ordering consistent with begchunk:endchunk 
! local ordering)

   call t_barrierf('sync_scat_ftoc', mpicom)
   call mpiscatterv(gfield_p, sndcnts, displs, mpir8, &
                    lfield_p, recvcnt, mpir8, 0, mpicom)

! copy into local chunked data structure

#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
   do i=1,nlcols
      cid = pgcols(beglcol+i)%chunk
      lcid = chunks(cid)%lcid
      lid = pgcols(beglcol+i)%ccol
      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                 lfield_p(f, m, l, i)
            end do
         end do
      end do
   end do
#else

! copy field into chunked data structure
! (pgcol ordering chosen to reflect begchunk:endchunk 
!  local ordering)

   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lcid = chunks(cid)%lcid
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                  globalfield(f, h1, m, h2, l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine scatter_field_to_chunk
!========================================================================

   subroutine scatter_field_to_chunk4(fdim,mdim,ldim, &
                                      hdim1d,globalfield,localchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute field
!          to decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: hdim1d    ! declared first horizontal index 
                                    ! dimension
   real(r4), intent(in) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
                                    ! global field

   real(r4), intent(out):: localchunks(fdim,pcols,mdim, &
                                       begchunk:endchunk,ldim) 
                                    ! local chunks

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local column index
   integer :: gcol                       ! global column index
   integer :: h1                         ! first horizontal dimension index
   integer :: h2                         ! second horizontal dimension index

#if ( defined SPMD )
   real(r4) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be scattered
   real(r4) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of scattered
                                         !  vector
   integer :: displs(0:npes-1)           ! scatter displacements
   integer :: sndcnts(0:npes-1)          ! scatter send counts
   integer :: recvcnt                    ! scatter receive count
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
   if (hdim1d < hdim1_d) then
      call endrun ('SCATTER_FIELD_TO_CHUNK4 error: hdim1d < hdim1_d')
   endif
#if ( defined SPMD )
   displs(0) = 0
   sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + sndcnts(p-1)
     sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   recvcnt = fdim*mdim*ldim*nlcols

   if (masterproc) then
      ! copy field into global (process-ordered) chunked data structure
      do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
         do i=1,ngcols_p
            cid  = pgcols(i)%chunk
            lid  = pgcols(i)%ccol
            gcol = chunks(cid)%gcol(lid)
            h2   = (gcol-1)/hdim1_d + 1
            h1   = mod((gcol-1),hdim1_d) + 1
            do m=1,mdim
               do f=1,fdim
                  gfield_p(f,m,l,i) = &
                     globalfield(f, h1, m, h2, l)
               end do
            end do
         end do
      end do
   endif

! scatter to other processes
! (pgcols ordering consistent with begchunk:endchunk 
!  local ordering)

   call t_barrierf('sync_scat_ftoc', mpicom)
   call mpiscatterv(gfield_p, sndcnts, displs, mpir4, &
                    lfield_p, recvcnt, mpir4, 0, mpicom)

! copy into local chunked data structure

#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
   do i=1,nlcols
      cid = pgcols(beglcol+i)%chunk
      lcid = chunks(cid)%lcid
      lid = pgcols(beglcol+i)%ccol
      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                 lfield_p(f, m, l, i)
            end do
         end do
      end do
   end do
#else

   ! copy field into chunked data structure
   ! (pgcol ordering chosen to reflect begchunk:endchunk 
   !  local ordering)
   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lcid = chunks(cid)%lcid
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                  globalfield(f, h1, m, h2, l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine scatter_field_to_chunk4
!========================================================================

   subroutine scatter_field_to_chunk_int(fdim,mdim,ldim, &
                                         hdim1d,globalfield,localchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute field
!          to decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: hdim1d    ! declared first horizontal index 
                                    ! dimension
   integer, intent(in) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
                                    ! global field

   integer, intent(out):: localchunks(fdim,pcols,mdim, &
                                       begchunk:endchunk,ldim) 
                                    ! local chunks

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local column index
   integer :: gcol                       ! global column index
   integer :: h1                         ! first horizontal dimension index
   integer :: h2                         ! second horizontal dimension index

#if ( defined SPMD )
   integer gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be scattered
   integer lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of scattered
                                         !  vector
   integer :: displs(0:npes-1)           ! scatter displacements
   integer :: sndcnts(0:npes-1)          ! scatter send counts
   integer :: recvcnt                    ! scatter receive count
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
   if (hdim1d < hdim1_d) then
      call endrun ('SCATTER_FIELD_TO_CHUNK_INT error: hdim1d < hdim1_d')
   endif
#if ( defined SPMD )
   displs(0) = 0
   sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + sndcnts(p-1)
     sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   recvcnt = fdim*mdim*ldim*nlcols

   if (masterproc) then

! copy field into global (process-ordered) chunked data structure

      do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
         do i=1,ngcols_p
            cid = pgcols(i)%chunk
            lid = pgcols(i)%ccol
            gcol = chunks(cid)%gcol(lid)
            h2   = (gcol-1)/hdim1_d + 1
            h1   = mod((gcol-1),hdim1_d) + 1
            do m=1,mdim
               do f=1,fdim
                  gfield_p(f,m,l,i) = &
                     globalfield(f, h1, m, h2, l)
               end do
            end do
         end do
      end do
   endif

! scatter to other processes
! (pgcols ordering consistent with begchunk:endchunk 
!  local ordering)

   call t_barrierf('sync_scat_ftoc', mpicom)
   call mpiscatterv(gfield_p, sndcnts, displs, mpiint, &
                    lfield_p, recvcnt, mpiint, 0, mpicom)

! copy into local chunked data structure

#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
   do i=1,nlcols
      cid = pgcols(beglcol+i)%chunk
      lcid = chunks(cid)%lcid
      lid = pgcols(beglcol+i)%ccol
      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                 lfield_p(f, m, l, i)
            end do
         end do
      end do
   end do
#else

! copy field into chunked data structure
! (pgcol ordering chosen to reflect begchunk:endchunk 
!  local ordering)
   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lcid = chunks(cid)%lcid
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                  globalfield(f, h1, m, h2, l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine scatter_field_to_chunk_int
!
!========================================================================
!
!  subroutine chunk_to_buff(fdim,mdim,localchunks,lbuff)
!
!-----------------------------------------------------------------------
!
! Purpose: Copy from local chunk data structure
!          to local buffer.  Needed for cpl6.
!          (local = assigned to same process)
!
! Method:
!
! Author: Pat Worley and Robert Jacob
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!  integer, intent(in) :: fdim      ! declared length of first lbuff dimension
!  integer, intent(in) :: mdim      ! declared length of middle lbuff dimension
!  real(r8), intent(in):: localchunks(pcols,mdim, begchunk:endchunk) ! local chunks
!
!  real(r8), intent(out) :: lbuff(fdim,mdim) ! local buff
!
!---------------------------Local workspace-----------------------------
!  integer :: i,j,m,n                  ! loop indices
!
!  integer, save :: numcols = 0
!  integer, allocatable, save :: columnid(:), chunkid(:)
!-----------------------------------------------------------------------
!
!  if (numcols .eq. 0) then
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!           endif
!        endif
!     enddo
!     allocate(columnid(1:n))
!     allocate(chunkid(1:n))
!
!     n = 0
!     do i=1,ngcols
!        if (dyn_to_latlon_gcol_map(i) .ne. -1) then
!           if(chunks(knuhcs(i)%chunkid)%owner .eq. iam) then
!              n = n + 1
!              columnid(n) = knuhcs(i)%col
!              chunkid(n)  = chunks(knuhcs(i)%chunkid)%lcid
!           endif
!        endif
!     end do
!
!     numcols = n
!  endif
!
!  if (numcols .gt. fdim) call endrun('chunk_to_buff')
!  do m=1,mdim
!#ifdef CPRCRAY
!!dir$ concurrent
!!dir$ prefervector, preferstream
!#endif
!     do n = 1, numcols
!        lbuff(n,m) = localchunks(columnid(n),m,chunkid(n))
!     end do
!  end do
!
!  return
!  end subroutine chunk_to_buff
!
!
!========================================================================
!
   subroutine gather_chunk_to_field(fdim,mdim,ldim, &
                                     hdim1d,localchunks,globalfield)

!----------------------------------------------------------------------- 
! 
! Purpose: Reconstruct field
!          from decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
   use spmd_utils,    only: fc_gatherv
#endif
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: hdim1d    ! declared first horizontal index 
                                    ! dimension
   real(r8), intent(in):: localchunks(fdim,pcols,mdim, &
                                      begchunk:endchunk,ldim) 
                                    ! local chunks

   real(r8), intent(out) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
                                    ! global field

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local column index
   integer :: gcol                       ! global column index
   integer :: h1                         ! first horizontal dimension index
   integer :: h2                         ! second horizontal dimension index

#if ( defined SPMD )
   real(r8) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be gathered
   real(r8) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of gather
                                         !  vector
   integer :: displs(0:npes-1)           ! gather displacements
   integer :: rcvcnts(0:npes-1)          ! gather receive count
   integer :: sendcnt                    ! gather send counts
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
   if (hdim1d < hdim1_d) then
      call endrun ('GATHER_CHUNK_TO_FIELD error: hdim1d < hdim1_d')
   endif
#if ( defined SPMD )
   displs(0) = 0
   rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + rcvcnts(p-1)
     rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   sendcnt = fdim*mdim*ldim*nlcols

! copy into local gather data structure

   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,nlcols
         cid = pgcols(beglcol+i)%chunk
         lcid = chunks(cid)%lcid
         lid = pgcols(beglcol+i)%ccol
         do m=1,mdim
            do f=1,fdim
               lfield_p(f, m, l, i) = &
                  localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

! gather from other processes

   call t_barrierf('sync_gath_ctof', mpicom)
   call fc_gatherv(lfield_p, sendcnt, mpir8, &
                   gfield_p, rcvcnts, displs, mpir8, 0, mpicom)

   if (masterproc) then

! copy gathered columns into lon/lat field

#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do l=1,ldim
            do m=1,mdim
               do f=1,fdim
                  globalfield(f, h1, m, h2, l)    &
                  = gfield_p(f,m,l,i)
               end do
            end do
         end do
      end do
   endif
   call mpibarrier(mpicom)
#else

   ! copy chunked data structure into dynamics field
   ! (pgcol ordering chosen to reflect begchunk:endchunk 
   !  local ordering)
   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lcid = chunks(cid)%lcid
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do m=1,mdim
            do f=1,fdim
               globalfield(f, h1, m, h2, l)    &
               = localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine gather_chunk_to_field

!
!========================================================================
!
   subroutine gather_chunk_to_field4 (fdim,mdim,ldim, &
                                      hdim1d,localchunks,globalfield)

!----------------------------------------------------------------------- 
! 
! Purpose: Reconstruct field
!          from decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
   use spmd_utils,    only: fc_gathervr4
#endif
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: hdim1d    ! declared first horizontal index 
                                    ! dimension
   real(r4), intent(in):: localchunks(fdim,pcols,mdim, &
                                      begchunk:endchunk,ldim) 
                                    ! local chunks

   real(r4), intent(out) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) 
                                    ! global field

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local column index
   integer :: gcol                       ! global column index
   integer :: h1                         ! first horizontal dimension index
   integer :: h2                         ! second horizontal dimension index

#if ( defined SPMD )
   real(r4) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be gathered
   real(r4) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of gather
                                         !  vector
   integer :: displs(0:npes-1)           ! gather displacements
   integer :: rcvcnts(0:npes-1)          ! gather receive count
   integer :: sendcnt                    ! gather send counts
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
   if (hdim1d < hdim1_d) then
      call endrun ('GATHER_CHUNK_TO_FIELD4 error: hdim1d < hdim1_d')
   endif
#if ( defined SPMD )
   displs(0) = 0
   rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + rcvcnts(p-1)
     rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   sendcnt = fdim*mdim*ldim*nlcols

! copy into local gather data structure

   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,nlcols
         cid = pgcols(beglcol+i)%chunk
         lcid = chunks(cid)%lcid
         lid = pgcols(beglcol+i)%ccol
         do m=1,mdim
            do f=1,fdim
               lfield_p(f, m, l, i) = &
                  localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

! gather from other processes

   call t_barrierf('sync_gath_ctof', mpicom)
   call fc_gathervr4(lfield_p, sendcnt, mpir4, &
                     gfield_p, rcvcnts, displs, mpir4, 0, mpicom)

   if (masterproc) then

! copy gathered columns into lon/lat field

#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do l=1,ldim
            do m=1,mdim
               do f=1,fdim
                  globalfield(f, h1, m, h2, l)    &
                  = gfield_p(f,m,l,i)
               end do
            end do
         end do
      end do
   endif

#else

! copy chunked data structure into dynamics field
! (pgcol ordering chosen to reflect begchunk:endchunk 
!  local ordering)

   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lcid = chunks(cid)%lcid
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do m=1,mdim
            do f=1,fdim
               globalfield(f, h1, m, h2, l)    &
               = localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine gather_chunk_to_field4

!
!========================================================================
!
   subroutine gather_chunk_to_field_int (fdim,mdim,ldim, &
                                         hdim1d,localchunks,globalfield)

!----------------------------------------------------------------------- 
! 
! Purpose: Reconstruct field
!          from decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
   use spmd_utils,    only: fc_gathervint
#endif
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: hdim1d    ! declared first horizontal index 
                                    ! dimension
   integer, intent(in):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks

   integer, intent(out) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim) ! global field

!---------------------------Local workspace-----------------------------

   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local column index
   integer :: gcol                       ! global column index
   integer :: h1                         ! first horizontal dimension index
   integer :: h2                         ! second horizontal dimension index

#if ( defined SPMD )
   integer gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be gathered
   integer lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of gather
                                         !  vector
   integer :: displs(0:npes-1)           ! gather displacements
   integer :: rcvcnts(0:npes-1)          ! gather receive count
   integer :: sendcnt                    ! gather send counts
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
   if (hdim1d < hdim1_d) then
      call endrun ('GATHER_CHUNK_TO_FIELD_INT error: hdim1d < hdim1_d')
   endif
#if ( defined SPMD )
   displs(0) = 0
   rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + rcvcnts(p-1)
     rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   sendcnt = fdim*mdim*ldim*nlcols

! copy into local gather data structure

   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,nlcols
         cid = pgcols(beglcol+i)%chunk
         lcid = chunks(cid)%lcid
         lid = pgcols(beglcol+i)%ccol
         do m=1,mdim
            do f=1,fdim
               lfield_p(f, m, l, i) = &
                  localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

! gather from other processes

   call t_barrierf('sync_gath_ctof', mpicom)
   call fc_gathervint(lfield_p, sendcnt, mpiint, &
                      gfield_p, rcvcnts, displs, mpiint, 0, mpicom)

   if (masterproc) then

! copy gathered columns into lon/lat field

#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do l=1,ldim
            do m=1,mdim
               do f=1,fdim
                  globalfield(f, h1, m, h2, l)    &
                  = gfield_p(f,m,l,i)
               end do
            end do
         end do
      end do
   endif

#else

   ! copy chunked data structure into lon/lat field
   ! (pgcol ordering chosen to reflect begchunk:endchunk 
   !  local ordering)
   do l=1,ldim
#ifdef CPRCRAY
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
#endif
      do i=1,ngcols_p
         cid  = pgcols(i)%chunk
         lcid = chunks(cid)%lcid
         lid  = pgcols(i)%ccol
         gcol = chunks(cid)%gcol(lid)
         h2   = (gcol-1)/hdim1_d + 1
         h1   = mod((gcol-1),hdim1_d) + 1
         do m=1,mdim
            do f=1,fdim
               globalfield(f, h1, m, h2, l)    &
               = localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine gather_chunk_to_field_int

!
!========================================================================
!
   subroutine write_field_from_chunk(iu,fdim,mdim,ldim,localchunks)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Write field from decomposed chunk data 
!          structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!------------------------------Arguments--------------------------------
   integer, intent(in) :: iu        ! logical unit
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   real(r8), intent(in):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks

!---------------------------Local workspace-----------------------------

   integer :: ioerr                 ! error return

   real(r8), allocatable :: globalfield(:,:,:,:,:)
                                    ! global field
!-----------------------------------------------------------------------

   allocate(globalfield(fdim,hdim1_d,mdim,hdim2_d,ldim))

   call gather_chunk_to_field (fdim,mdim,ldim,hdim1_d,localchunks,globalfield)
                               
   if (masterproc) then
      write (iu,iostat=ioerr) globalfield
      if (ioerr /= 0 ) then
         write(iulog,*) 'WRITE_FIELD_FROM_CHUNK ioerror ', ioerr,' on i/o unit = ',iu
         call endrun
      end if
   endif

   deallocate(globalfield)

   return
   end subroutine write_field_from_chunk

!
!========================================================================
!
   subroutine read_chunk_from_field(iu,fdim,mdim,ldim,localchunks)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Write field from decomposed chunk data 
!          structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!------------------------------Arguments--------------------------------
   integer, intent(in) :: iu        ! logical unit
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension

   real(r8), intent(out):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks

!---------------------------Local workspace-----------------------------

   integer :: ioerr                 ! error return

   real(r8), allocatable :: globalfield(:,:,:,:,:)
                                    ! global field
!-----------------------------------------------------------------------

   allocate(globalfield(fdim,hdim1_d,mdim,hdim2_d,ldim))

   if (masterproc) then
      read (iu,iostat=ioerr) globalfield
      if (ioerr /= 0 ) then
         write(iulog,*) 'READ_CHUNK_FROM_FIELD ioerror ', ioerr,' on i/o unit = ',iu
         call endrun
      end if
   endif

   call scatter_field_to_chunk (fdim,mdim,ldim,hdim1_d,globalfield,localchunks)

   deallocate(globalfield)

   return
   end subroutine read_chunk_from_field
!
!========================================================================

   subroutine transpose_block_to_chunk(record_size, block_buffer, &
                                       chunk_buffer, window)
                                       
!----------------------------------------------------------------------- 
! 
! Purpose: Transpose buffer containing decomposed 
!          fields to buffer
!          containing decomposed chunk data structures
! 
! Method: 
! 
! Author: Patrick Worley
! Modified: Art Mirin, Jan 04, to add support for mod_comm
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
# if defined(MODCM_DP_TRANSPOSE)
   use mod_comm, only: blockdescriptor, mp_sendirr, mp_recvirr,  &
                       get_partneroffset, max_nparcels
   use mpishorthand,  only : mpicom
# endif
   use spmd_utils,    only: altalltoallv
#endif
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 6000
!------------------------------Arguments--------------------------------
   integer, intent(in) :: record_size  ! per column amount of data 
   real(r8), intent(in) :: block_buffer(record_size*block_buf_nrecs)
                                       ! buffer of block data to be
                                       ! transposed
   real(r8), intent(out):: chunk_buffer(record_size*chunk_buf_nrecs)
                                       ! buffer of chunk data 
                                       ! transposed into
   integer, intent(in), optional :: window
                                       ! MPI-2 window id for
                                       ! chunk_buffer

!---------------------------Local workspace-----------------------------
#if ( defined SPMD )
   integer :: p                        ! loop indices
   integer :: bbuf_siz                 ! size of block_buffer
   integer :: cbuf_siz                 ! size of chunk_buffer
   integer :: lwindow                  ! placeholder for missing window
   integer :: lopt                     ! local copy of phys_alltoall
!
   logical, save :: first = .true.
   integer, allocatable, save :: sndcnts(:), sdispls(:)
   integer, allocatable, save :: rcvcnts(:), rdispls(:)
   integer, allocatable, save :: pdispls(:)
   integer, save :: prev_record_size = 0
# if defined(MODCM_DP_TRANSPOSE)
   type (blockdescriptor), allocatable, save :: sendbl(:), recvbl(:)
   integer ione, ierror, mod_method
# endif
!-----------------------------------------------------------------------
   if (first) then
! Compute send/recv/put counts and displacements
      allocate(sndcnts(0:npes-1))
      allocate(sdispls(0:npes-1))
      allocate(rcvcnts(0:npes-1))
      allocate(rdispls(0:npes-1))
      allocate(pdispls(0:npes-1))
!
# if defined(MODCM_DP_TRANSPOSE)
! This branch uses mod_comm. Admissable values of phys_alltoall are 
! 11,12 and 13. Each value corresponds to a different option 
! within mod_comm of implementing the communication. That option is expressed
! internally to mod_comm using the variable mod_method defined below; 
! mod_method will have values 0,1 or 2 and is defined as 
! phys_alltoall - modmin_alltoall, where modmin_alltoall equals 11.
! Also, sendbl and recvbl must have exactly npes elements, to match
! this size of the communicator, or the transpose will fail.
!
      if (phys_alltoall .ge. modmin_alltoall) then
         mod_method = phys_alltoall - modmin_alltoall
         ione = 1
         allocate( sendbl(0:npes-1) )
         allocate( recvbl(0:npes-1) )

         do p = 0,npes-1

            sendbl(p)%method = mod_method
            recvbl(p)%method = mod_method

            allocate( sendbl(p)%blocksizes(1) )
            allocate( sendbl(p)%displacements(1) )
            allocate( recvbl(p)%blocksizes(1) )
            allocate( recvbl(p)%displacements(1) )

         enddo

      endif
# endif

      first = .false.
   endif
!
   if (record_size .ne. prev_record_size) then
!
! Compute send/recv/put counts and displacements
      sdispls(0) = 0
      sndcnts(0) = record_size*btofc_blk_num(0)
      do p=1,npes-1
        sdispls(p) = sdispls(p-1) + sndcnts(p-1)
        sndcnts(p) = record_size*btofc_blk_num(p)
      enddo
!
      rdispls(0) = 0
      rcvcnts(0) = record_size*btofc_chk_num(0)
      do p=1,npes-1
         rdispls(p) = rdispls(p-1) + rcvcnts(p-1)
         rcvcnts(p) = record_size*btofc_chk_num(p)
      enddo
!
      call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
!
# if defined(MODCM_DP_TRANSPOSE)
      if (phys_alltoall .ge. modmin_alltoall) then
         do p = 0,npes-1

            sendbl(p)%type = MPI_DATATYPE_NULL
            if ( sndcnts(p) .ne. 0 ) then

               if (phys_alltoall .gt. modmin_alltoall) then
                  call MPI_TYPE_INDEXED(ione, sndcnts(p),   &
                       sdispls(p), mpir8, &
                       sendbl(p)%type, ierror)
                  call MPI_TYPE_COMMIT(sendbl(p)%type, ierror)
               endif

               sendbl(p)%blocksizes(1) = sndcnts(p)
               sendbl(p)%displacements(1) = sdispls(p)
               sendbl(p)%partneroffset = 0

            else

               sendbl(p)%blocksizes(1) = 0
               sendbl(p)%displacements(1) = 0
               sendbl(p)%partneroffset = 0

            endif
            sendbl(p)%nparcels = size(sendbl(p)%displacements)
            sendbl(p)%tot_size = sum(sendbl(p)%blocksizes)
            max_nparcels = max(max_nparcels, sendbl(p)%nparcels)

            recvbl(p)%type = MPI_DATATYPE_NULL
            if ( rcvcnts(p) .ne. 0) then

               if (phys_alltoall .gt. modmin_alltoall) then
                  call MPI_TYPE_INDEXED(ione, rcvcnts(p),   &
                       rdispls(p), mpir8, &
                       recvbl(p)%type, ierror)
                  call MPI_TYPE_COMMIT(recvbl(p)%type, ierror)
               endif

               recvbl(p)%blocksizes(1) = rcvcnts(p)
               recvbl(p)%displacements(1) = rdispls(p)
               recvbl(p)%partneroffset = 0 ! not properly initialized - do not use Mpi2
            else

               recvbl(p)%blocksizes(1) = 0
               recvbl(p)%displacements(1) = 0
               recvbl(p)%partneroffset = 0

            endif
            recvbl(p)%nparcels = size(recvbl(p)%displacements)
            recvbl(p)%tot_size = sum(recvbl(p)%blocksizes)
            max_nparcels = max(max_nparcels, recvbl(p)%nparcels)

         enddo

         call get_partneroffset(mpicom, sendbl, recvbl)

      endif
# endif
!
      prev_record_size = record_size
   endif
!
   call t_barrierf('sync_tran_btoc', mpicom)
   if (phys_alltoall < 0) then
      if ((max_nproc_vsmp > npes/2) .and. (nproc_busy_d > npes/2)) then
         lopt = 0
      else
         lopt = 1
      endif
   else
      lopt = phys_alltoall
      if ((lopt .eq. 2) .and. ( .not. present(window) )) lopt = 1
   endif
   if (lopt < 4) then
!
      bbuf_siz = record_size*block_buf_nrecs
      cbuf_siz = record_size*chunk_buf_nrecs
      if ( present(window) ) then
         call altalltoallv(lopt, iam, npes,    &
                           dp_coup_steps, dp_coup_proc, &
                           block_buffer, bbuf_siz, sndcnts, sdispls, mpir8, &
                           chunk_buffer, cbuf_siz, rcvcnts, rdispls, mpir8, &
                           msgtag, pdispls, mpir8, window, mpicom)
      else
         call altalltoallv(lopt, iam, npes,    &
                           dp_coup_steps, dp_coup_proc, &
                           block_buffer, bbuf_siz, sndcnts, sdispls, mpir8, &
                           chunk_buffer, cbuf_siz, rcvcnts, rdispls, mpir8, &
                           msgtag, pdispls, mpir8, lwindow, mpicom)
      endif
!
   else
!
# if defined(MODCM_DP_TRANSPOSE)
      call mp_sendirr(mpicom, sendbl, recvbl, block_buffer, chunk_buffer)
      call mp_recvirr(mpicom, sendbl, recvbl, block_buffer, chunk_buffer)
# else
      call mpialltoallv(block_buffer, sndcnts, sdispls, mpir8, &
                        chunk_buffer, rcvcnts, rdispls, mpir8, &
                        mpicom)
# endif
!
   endif
!
#endif
   return
   end subroutine transpose_block_to_chunk
!
!========================================================================

   subroutine block_to_chunk_send_pters(blockid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into send buffer where column from decomposed 
!          fields should be copied to
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid      ! block index
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offsets
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_blk_offset(blockid)%ncols > fdim) .or. &
       (btofc_blk_offset(blockid)%nlvls > ldim)) then
      write(iulog,*) "BLOCK_TO_CHUNK_SEND_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_blk_offset(blockid)%ncols,",", &
                  btofc_blk_offset(blockid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_blk_offset(blockid)%nlvls
      do i=1,btofc_blk_offset(blockid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_blk_offset(blockid)%pter(i,k))
      enddo
      do i=btofc_blk_offset(blockid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_blk_offset(blockid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine block_to_chunk_send_pters
!
!========================================================================

   subroutine block_to_chunk_recv_pters(lcid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into receive buffer where data for
!          decomposed chunk data structures should be copied from
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: lcid         ! local chunk id
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offset
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_chk_offset(lcid)%ncols > fdim) .or. &
       (btofc_chk_offset(lcid)%nlvls > ldim)) then
      write(iulog,*) "BLOCK_TO_CHUNK_RECV_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_chk_offset(lcid)%ncols,",", &
                  btofc_chk_offset(lcid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_chk_offset(lcid)%nlvls
      do i=1,btofc_chk_offset(lcid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_chk_offset(lcid)%pter(i,k))
      enddo
      do i=btofc_chk_offset(lcid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_chk_offset(lcid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine block_to_chunk_recv_pters
!
!========================================================================

   subroutine transpose_chunk_to_block(record_size, chunk_buffer, &
                                       block_buffer, window)
!----------------------------------------------------------------------- 
! 
! Purpose: Transpose buffer containing decomposed 
!          chunk data structures to buffer
!          containing decomposed fields 
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
# if defined(MODCM_DP_TRANSPOSE)
   use mod_comm, only: blockdescriptor, mp_sendirr, mp_recvirr,  &
                       get_partneroffset, max_nparcels
   use mpishorthand,  only : mpicom
# endif
   use spmd_utils,    only: altalltoallv
#endif
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 7000
!------------------------------Arguments--------------------------------
   integer, intent(in) :: record_size  ! per column amount of data 
   real(r8), intent(in):: chunk_buffer(record_size*chunk_buf_nrecs)
                                       ! buffer of chunk data to be
                                       ! transposed
   real(r8), intent(out) :: block_buffer(record_size*block_buf_nrecs)
                                       ! buffer of block data to
                                       ! transpose into
   integer, intent(in), optional :: window
                                       ! MPI-2 window id for
                                       ! chunk_buffer

!---------------------------Local workspace-----------------------------
#if ( defined SPMD )
   integer :: p                        ! loop indices
   integer :: bbuf_siz                 ! size of block_buffer
   integer :: cbuf_siz                 ! size of chunk_buffer
   integer :: lwindow                  ! placeholder for missing window
   integer :: lopt                     ! local copy of phys_alltoall
!
   logical, save :: first = .true.
   integer, allocatable, save :: sndcnts(:), sdispls(:)
   integer, allocatable, save :: rcvcnts(:), rdispls(:)
   integer, allocatable, save :: pdispls(:)
   integer, save :: prev_record_size = 0
# if defined(MODCM_DP_TRANSPOSE)
   type (blockdescriptor), allocatable, save :: sendbl(:), recvbl(:)
   integer ione, ierror, mod_method
# endif
!-----------------------------------------------------------------------
   if (first) then
! Compute send/recv/put counts and displacements
      allocate(sndcnts(0:npes-1))
      allocate(sdispls(0:npes-1))
      allocate(rcvcnts(0:npes-1))
      allocate(rdispls(0:npes-1))
      allocate(pdispls(0:npes-1))
!
# if defined(MODCM_DP_TRANSPOSE)
! This branch uses mod_comm. Admissable values of phys_alltoall are 
! 11,12 and 13. Each value corresponds to a differerent option 
! within mod_comm of implementing the communication. That option is expressed
! internally to mod_comm using the variable mod_method defined below; 
! mod_method will have values 0,1 or 2 and is defined as 
! phys_alltoall - modmin_alltoall, where modmin_alltoall equals 11.
! Also, sendbl and recvbl must have exactly npes elements, to match
! this size of the communicator, or the transpose will fail.
!
      if (phys_alltoall .ge. modmin_alltoall) then
         mod_method = phys_alltoall - modmin_alltoall
         ione = 1
         allocate( sendbl(0:npes-1) )
         allocate( recvbl(0:npes-1) )

         do p = 0,npes-1

            sendbl(p)%method = mod_method
            recvbl(p)%method = mod_method

            allocate( sendbl(p)%blocksizes(1) )
            allocate( sendbl(p)%displacements(1) )
            allocate( recvbl(p)%blocksizes(1) )
            allocate( recvbl(p)%displacements(1) )

         enddo

      endif
# endif
!
      first = .false.
   endif
!
   if (record_size .ne. prev_record_size) then
!
! Compute send/recv/put counts and displacements
      sdispls(0) = 0
      sndcnts(0) = record_size*btofc_chk_num(0)
      do p=1,npes-1
        sdispls(p) = sdispls(p-1) + sndcnts(p-1)
        sndcnts(p) = record_size*btofc_chk_num(p)
      enddo
!
      rdispls(0) = 0
      rcvcnts(0) = record_size*btofc_blk_num(0)
      do p=1,npes-1
         rdispls(p) = rdispls(p-1) + rcvcnts(p-1)
         rcvcnts(p) = record_size*btofc_blk_num(p)
      enddo
!
      call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
!
# if defined(MODCM_DP_TRANSPOSE)
      if (phys_alltoall .ge. modmin_alltoall) then
         do p = 0,npes-1

            sendbl(p)%type = MPI_DATATYPE_NULL
            if ( sndcnts(p) .ne. 0 ) then

               if (phys_alltoall .gt. modmin_alltoall) then
                  call MPI_TYPE_INDEXED(ione, sndcnts(p),   &
                       sdispls(p), mpir8, &
                       sendbl(p)%type, ierror)
                  call MPI_TYPE_COMMIT(sendbl(p)%type, ierror)
               endif

               sendbl(p)%blocksizes(1) = sndcnts(p)
               sendbl(p)%displacements(1) = sdispls(p)
               sendbl(p)%partneroffset = 0

            else

               sendbl(p)%blocksizes(1) = 0
               sendbl(p)%displacements(1) = 0
               sendbl(p)%partneroffset = 0

            endif
            sendbl(p)%nparcels = size(sendbl(p)%displacements)
            sendbl(p)%tot_size = sum(sendbl(p)%blocksizes)
            max_nparcels = max(max_nparcels, sendbl(p)%nparcels)

            recvbl(p)%type = MPI_DATATYPE_NULL
            if ( rcvcnts(p) .ne. 0) then

               if (phys_alltoall .gt. modmin_alltoall) then
                  call MPI_TYPE_INDEXED(ione, rcvcnts(p),   &
                       rdispls(p), mpir8, &
                       recvbl(p)%type, ierror)
                  call MPI_TYPE_COMMIT(recvbl(p)%type, ierror)
               endif

               recvbl(p)%blocksizes(1) = rcvcnts(p)
               recvbl(p)%displacements(1) = rdispls(p)
               recvbl(p)%partneroffset = 0 ! not properly initialized - do not use Mpi2
            else

               recvbl(p)%blocksizes(1) = 0
               recvbl(p)%displacements(1) = 0
               recvbl(p)%partneroffset = 0

            endif
            recvbl(p)%nparcels = size(recvbl(p)%displacements)
            recvbl(p)%tot_size = sum(recvbl(p)%blocksizes)
            max_nparcels = max(max_nparcels, recvbl(p)%nparcels)

         enddo

         call get_partneroffset(mpicom, sendbl, recvbl)

      endif
# endif
!
      prev_record_size = record_size
   endif
!
   call t_barrierf('sync_tran_ctob', mpicom)
   if (phys_alltoall < 0) then
      if ((max_nproc_vsmp > npes/2) .and. (nproc_busy_d > npes/2)) then
         lopt = 0
      else
         lopt = 1
      endif
   else
      lopt = phys_alltoall
      if ((lopt .eq. 2) .and. ( .not. present(window) )) lopt = 1
   endif
   if (lopt < 4) then
!
      bbuf_siz = record_size*block_buf_nrecs
      cbuf_siz = record_size*chunk_buf_nrecs
      if ( present(window) ) then
         call altalltoallv(lopt, iam, npes,    &
                           dp_coup_steps, dp_coup_proc, &
                           chunk_buffer, cbuf_siz, sndcnts, sdispls, mpir8, &
                           block_buffer, bbuf_siz, rcvcnts, rdispls, mpir8, &
                           msgtag, pdispls, mpir8, window, mpicom)
      else
         call altalltoallv(lopt, iam, npes,    &
                           dp_coup_steps, dp_coup_proc, &
                           chunk_buffer, cbuf_siz, sndcnts, sdispls, mpir8, &
                           block_buffer, bbuf_siz, rcvcnts, rdispls, mpir8, &
                           msgtag, pdispls, mpir8, lwindow, mpicom)
      endif
!
   else
# if defined(MODCM_DP_TRANSPOSE)
      call mp_sendirr(mpicom, sendbl, recvbl, block_buffer, chunk_buffer)
      call mp_recvirr(mpicom, sendbl, recvbl, block_buffer, chunk_buffer)
# else
      call mpialltoallv(chunk_buffer, sndcnts, sdispls, mpir8, &
                        block_buffer, rcvcnts, rdispls, mpir8, &
                        mpicom)
# endif
!
   endif
!
#endif

   return
   end subroutine transpose_chunk_to_block
!
!========================================================================

   subroutine chunk_to_block_send_pters(lcid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into send buffer where data for
!          decomposed chunk data structures should be copied to
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: lcid         ! local chunk id
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offset
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_chk_offset(lcid)%ncols > fdim) .or. &
       (btofc_chk_offset(lcid)%nlvls > ldim)) then
      write(iulog,*) "CHUNK_TO_BLOCK_SEND_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_chk_offset(lcid)%ncols,",", &
                  btofc_chk_offset(lcid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_chk_offset(lcid)%nlvls
      do i=1,btofc_chk_offset(lcid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_chk_offset(lcid)%pter(i,k))
      enddo
      do i=btofc_chk_offset(lcid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_chk_offset(lcid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine chunk_to_block_send_pters
!
!========================================================================

   subroutine chunk_to_block_recv_pters(blockid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into receive buffer where column from decomposed 
!          fields should be copied from
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid      ! block index
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offsets
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_blk_offset(blockid)%ncols > fdim) .or. &
       (btofc_blk_offset(blockid)%nlvls > ldim)) then
      write(iulog,*) "CHUNK_TO_BLOCK_RECV_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_blk_offset(blockid)%ncols,",", &
                  btofc_blk_offset(blockid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_blk_offset(blockid)%nlvls
      do i=1,btofc_blk_offset(blockid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_blk_offset(blockid)%pter(i,k))
      enddo
      do i=btofc_blk_offset(blockid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_blk_offset(blockid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine chunk_to_block_recv_pters
!
!========================================================================

   subroutine create_chunks(opt, chnks_per_thrd, cols_per_chnk, &
                            cols_per_chnk_max, cols_per_chnk_mult, &
                            pcols_proc)
!----------------------------------------------------------------------- 
! 
! Purpose: Decompose physics computational grid into chunks, for
!          improved serial efficiency and parallel load balance.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plev
   use dyn_grid, only: get_block_bounds_d, get_block_gcol_cnt_d, &
                       get_gcol_block_cnt_d, get_gcol_block_d, &
                       get_block_owner_d, get_block_gcol_d
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: opt           
      ! chunking option
      !  0: chunks may cross block boundaries, but retain same
      !     process mapping as blocks. If possible, columns assigned
      !     as day/night pairs. Columns (or pairs) are wrap-mapped.
      !     May not work with vertically decomposed blocks. (default)
      !  1: chunks may cross block boundaries, but retain same
      !     SMP-node mapping as blocks.  If possible, columns assigned
      !     as day/night pairs.  Columns (or pairs) are wrap-mapped.
      !     May not work with vertically decomposed blocks.
      !  2: 2-column day/night and season column pairs wrap-mapped
      !     to chunks to also balance assignment of polar, mid-latitude, 
      !     and equatorial columns across  chunks.
      !  3: same as 1 except that SMP defined to be pairs of consecutive
      !     processes
      !  4: chunks may cross block boundaries, but retain same
      !     process mapping as blocks. Columns assigned to chunks
      !     in block ordering.
      !     May not work with vertically decomposed blocks.
      !  5: Chunks do not cross block boundaries, and are block-mapped.
   integer, intent(in)  :: chnks_per_thrd 
      ! target number of chunks per thread
   integer, intent(in)  :: cols_per_chnk
      ! Physics fields data structure (chunk) first dimension option:
      ! >0: pcols is set to cols_per_chnk
      ! <1: calculate pcols based on opt, chnks_per_thrd,
      !      number columns and threads in virtual SMP, and relative
      !      costs per column (if provided), attempting to
      !      minimize wasted space and number of chunks, subject to
      !      additional requirements sepcified by cols_per_chnk_max
      !      and cols_per_chnk_mult
   integer, intent(in)  :: cols_per_chnk_max
      ! Physics fields data structure (chunk) first dimension upper bound:
      ! >0: if cols_per_chnk <= 0, then cols_per_chnk_max is an
      !      upper bound on the calculated pcols
      ! <1: ignore
   integer, intent(in)  :: cols_per_chnk_mult
      ! Physics fields data structure (chunk) first dimension factor:
      ! >0: if cols_per_chnk <= 0 and cols_per_chnk_max > 0 and
      !      cols_per_chnk_max >= cols_per_chnk_mult, then pcols is
      !      required to be a multiple of cols_per_chnk_mult. Otherwise
      !      it is ignored.
      ! <1: ignore
   integer, intent(out) :: pcols_proc(0:npes-1)
      ! pcols for all chunks assigned to each process
!---------------------------Local workspace-----------------------------
   integer :: i, j, p                    ! loop indices
   integer :: proc_vsmp_map(0:npes-1)    ! process/virtual SMP node map
   integer :: firstblock, lastblock      ! global block index bounds
   integer :: maxblksiz                  ! maximum number of columns in a dynamics block
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: nvsmp, nvsmp2              ! virtual SMP node counts and indices
   integer :: curgcol, twingcol          ! global physics and dynamics column indices
   integer :: smp                        ! SMP node index (both virtual and actual)
   integer :: tmp_chnks_per_thrd         ! temporary used to calculate number of
                                         !  chunks per thread
   integer :: max_pcols                  ! max pcols over all virtual SMPs
   integer :: lmax_pcols                 ! max pcols over chunks in one virtual SMP
   integer :: cid                        ! chunk id
   integer :: jb, ib                     ! global block and columns indices
   integer :: blksiz                     ! current block size
   integer :: ntmp1, ntmp2, nlchunks     ! work variables
   integer :: lcol                       ! chunk column index
   integer :: max_ncols                  ! upper bound on number of columns in a block
   integer :: ncols                      ! number of columns in current chunk

   logical :: error                      ! error flag 

   ! indices for dynamics columns in given block
   integer, dimension(:), allocatable :: cols_d

   ! number of MPI processes per virtual SMP node (0:nvsmp-1)
   integer, dimension(:), allocatable :: nsmpprocs      

   ! flag indicating whether a process is busy or idle during the dynamics (0:npes-1)
   logical, dimension(:), allocatable :: proc_busy_d

   ! flag indicating whether any of the processes assigned to an SMP node are busy 
   ! during the dynamics, or whether all of them are idle (0:nsmps-1)
   logical, dimension(:), allocatable :: smp_busy_d

   ! actual SMP node/virtual SMP node map (0:nsmps-1)    
   integer, dimension(:), allocatable :: smp_vsmp_map

   ! column/virtual SMP node map (ngcols)
   integer, dimension(:), allocatable :: col_vsmp_map

   ! number of columns assigned to a given virtual SMP node (0:nvsmp-1)
   integer, dimension(:), allocatable :: nvsmpcolumns

   ! number of OpenMP threads per virtual SMP node (0:nvsmp-1)
   integer, dimension(:), allocatable :: nvsmpthreads

   ! number of chunks assigned to a given virtual SMP node (0:nvsmp-1)
   integer, dimension(:), allocatable :: nvsmpchunks
                                         
   ! space allocated for columns assigned to a chunk in a given virtual SMP node (0:nvsmp-1)
   integer, dimension(:), allocatable :: pcols_vsmp
                                         
   ! maximum number of columns assigned to a chunk in a given virtual SMP node (0:nvsmp-1)
   integer, dimension(:), allocatable :: maxcol_chk
                                         
   ! number of chunks in given virtual SMP node receiving maximum number of columns 
   ! (0:nvsmp-1)
   integer, dimension(:), allocatable :: maxcol_chks

   ! maximum number of columns assigned to a dynamics block in given virtual SMP node
   ! (0:nsvmp-1)
   integer, dimension(:), allocatable :: maxblksiz_vsmp

   ! chunk id virtual offset (0:nvsmp-1)
   integer, dimension(:), allocatable :: cid_offset

   ! process-local chunk id (0:nvsmp-1)
   integer, dimension(:), allocatable :: local_cid

   ! permutation array used to sort columns by their computation cost
   integer, dimension(:), allocatable :: cdex

   ! array used to mark whether a column has been assigned when sorting
   ! by dynamics block
   logical, dimension(:), allocatable :: udex

   ! min heap array used to maintain chunks sorted by assigned computational
   ! cost. A separate heap is created for each virtual smp, but all
   ! heaps are implemented in this one 1-D array
   integer, dimension(:), allocatable :: heap
   integer, dimension(:), allocatable :: heap_len

!-----------------------------------------------------------------------

! Make sure that proc_smp_map is set appropriately even if not using MPI
#if ( ! defined SPMD )
   proc_smp_map(0) = 0
#endif

!
! Determine index range for dynamics blocks
!
   call get_block_bounds_d(firstblock,lastblock)

!
! Determine maximum number of columns in a block
!
   maxblksiz = 0
   do jb=firstblock,lastblock
      maxblksiz = max(maxblksiz,get_block_gcol_cnt_d(jb))
   enddo

!
! Determine which (and how many) processes are assigned
! dynamics blocks
!
   allocate( proc_busy_d(0:npes-1) )
   proc_busy_d = .false.
   nproc_busy_d = 0
   do jb=firstblock,lastblock
      p = get_block_owner_d(jb)
      if (.not. proc_busy_d(p) ) then
         proc_busy_d(p) = .true.
         nproc_busy_d = nproc_busy_d + 1
      endif
   enddo

!
! Determine virtual SMP count and processes/virtual SMP map.
!  If option 0 or >3, pretend that each SMP has only one process. 
!  If option 1, use SMP information.
!  If option 2, pretend that all processes are in one SMP node. 
!  If option 3, pretend that each SMP node is made up of two
!     processes, chosen to maximize load-balancing opportunities.
!
!  For all options < 5, if there are "idle" dynamics processes, 
!     assign them to the virtual SMP nodes in wrap fashion.
!     Communication between the active and idle dynamics 
!     processes is scatter/gather (no communications between 
!     idle dynamics processes) so there is no advantage to 
!     blocking the idle processes in these assignments.
!
   if ((opt <= 0) .or. (opt == 4)) then

      ! Assign active dynamics processes to virtual SMP nodes
      nvsmp = 0
      do p=0,npes-1
         if (proc_busy_d(p)) then
            proc_vsmp_map(p) = nvsmp
            nvsmp = nvsmp + 1
         endif
      enddo

      ! Assign idle dynamics processes to virtual SMP nodes (wrap map)
      nvsmp2 = 0
      do p=0,npes-1
         if (.not. proc_busy_d(p)) then
            proc_vsmp_map(p) = nvsmp2
            nvsmp2 = mod(nvsmp2+1,nvsmp)
         endif
      enddo

   elseif (opt == 1) then

      allocate( smp_busy_d(0:nsmps-1) )
      allocate( smp_vsmp_map(0:nsmps-1) )

      ! Determine SMP nodes assigned dynamics blocks
      smp_busy_d = .false.
      do p=0,npes-1
         if ( proc_busy_d(p) ) then
            smp = proc_smp_map(p)
            smp_busy_d(smp) = .true.
         endif
      enddo

      ! Determine number of SMP nodes assigned dynamics blocks
      nvsmp = 0
      do smp=0,nsmps-1
         if (smp_busy_d(smp)) then
            smp_vsmp_map(smp) = nvsmp
            nvsmp = nvsmp + 1
         endif
      enddo

      ! Assign processes in active dynamics SMP nodes to virtual SMP nodes
      do p=0,npes-1
         smp = proc_smp_map(p)
         if (smp_busy_d(smp)) then
            proc_vsmp_map(p) = smp_vsmp_map(smp)
         endif
      enddo

      ! Assign processes in idle dynamics SMP nodes to virtual SMP nodes (wrap map)
      nvsmp2 = 0
      do p=0,npes-1
         smp = proc_smp_map(p)
         if (.not. smp_busy_d(smp)) then
            proc_vsmp_map(p) = nvsmp2
            nvsmp2 = mod(nvsmp2+1,nvsmp)
         endif
      enddo

      deallocate( smp_busy_d )
      deallocate( smp_vsmp_map )

   elseif (opt == 2) then

      nvsmp = 1
      do p=0,npes-1
         proc_vsmp_map(p) = 0
      enddo

   elseif (opt == 3) then

      ! Find active process partners
      proc_vsmp_map = -1
      call find_partners(opt,proc_busy_d,nvsmp,proc_vsmp_map)

      ! Assign unassigned (idle dynamics) processes to virtual SMP nodes 
      ! (wrap map)
      nvsmp2 = 0
      do p=0,npes-1
         if (proc_vsmp_map(p) .eq. -1) then
            proc_vsmp_map(p) = nvsmp2
            nvsmp2 = mod(nvsmp2+1,nvsmp)
         endif
      enddo

   else

      nvsmp = npes
      do p=0,npes-1
         proc_vsmp_map(p) = p
      enddo

   endif

   deallocate( proc_busy_d )

!
! Determine maximum number of processes assigned to a single 
! virtual SMP node
!
   allocate( nsmpprocs(0:nvsmp-1) )

   nsmpprocs(:) = 0
   do p=0,npes-1
      smp = proc_vsmp_map(p)
      nsmpprocs(smp) = nsmpprocs(smp) + 1
   enddo
   max_nproc_vsmp = maxval(nsmpprocs)

   deallocate( nsmpprocs )   

!
! Determine number of columns assigned to each
! virtual SMP in block decomposition
!
   allocate( col_vsmp_map(ngcols) )

   col_vsmp_map(:) = -1
   error = .false.
   do i=1,ngcols
      if (dyn_to_latlon_gcol_map(i) .ne. -1) then
         block_cnt = get_gcol_block_cnt_d(i)
         call get_gcol_block_d(i,block_cnt,blockids,bcids)
         do jb=1,block_cnt
            p = get_block_owner_d(blockids(jb))
            if (col_vsmp_map(i) .eq. -1) then
               col_vsmp_map(i) = proc_vsmp_map(p)
            elseif (col_vsmp_map(i) .ne. proc_vsmp_map(p)) then
               error = .true.
            endif
         enddo
      endif
   enddo
   if (error) then
      write(iulog,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
               "but vertical decomposition not limited to virtual SMP"
      call endrun()
   endif  

   allocate( nvsmpcolumns(0:nvsmp-1) )

   nvsmpcolumns(:) = 0
   error = .false.
   do i=1,ngcols_p
      curgcol = latlon_to_dyn_gcol_map(i)
      smp = col_vsmp_map(curgcol)
      if (smp >= 0) then
         nvsmpcolumns(smp) = nvsmpcolumns(smp) + 1
      else
         error = .true.
         exit
     endif
   end do
   if (error) then
      write(iulog,*) "PHYS_GRID_INIT error: no SMP found for ", &
                     "column ", curgcol
      call endrun()
   endif

!
!  Allocate other work space
!
   allocate( nvsmpthreads  (0:nvsmp-1) )
   allocate( nvsmpchunks   (0:nvsmp-1) )
   allocate( pcols_vsmp    (0:nvsmp-1) )
   allocate( maxcol_chk    (0:nvsmp-1) )
   allocate( maxcol_chks   (0:nvsmp-1) )
   allocate( maxblksiz_vsmp(0:nvsmp-1) )
   allocate( cid_offset    (0:nvsmp-1) )
   allocate( local_cid     (0:nvsmp-1) )
   allocate( cols_d      (1:maxblksiz) )

!
! Calculate number of threads available in each virtual SMP node. 
!
   nvsmpthreads(:) = 0
   do p=0,npes-1
      smp = proc_vsmp_map(p)
      nvsmpthreads(smp) = nvsmpthreads(smp) + npthreads(p)
   enddo

!
! Calculate maximum block size for each virtual SMP
!
   maxblksiz_vsmp(:) = 0
   do j=firstblock,lastblock
      p = get_block_owner_d(j)
      smp = proc_vsmp_map(p)
      blksiz = get_block_gcol_cnt_d(j)
      if (blksiz > maxblksiz_vsmp(smp)) maxblksiz_vsmp(smp) = blksiz
   enddo

!
! Calculate pcols for chunks in each virtual SMP
!
#ifdef PPCOLS
   ! Use compile-time value
   pcols_vsmp(:) = pcols
#else
   if (cols_per_chnk > 0) then
      ! Use runtime value
      pcols_vsmp(:) = cols_per_chnk
   else
      do smp=0,nvsmp-1
         ! Pick pcols to minimize wasted space and the number of chunks per thread
         pcols_vsmp(smp) = &
            ceiling(real(nvsmpcolumns(smp),r8)/real(nvsmpthreads(smp)*chnks_per_thrd,r8))
         if (opt == 5) then
            ! pcols should not be (much) larger than the maximum block size
            if (maxblksiz_vsmp(smp) < pcols_vsmp(smp)) pcols_vsmp(smp) = maxblksiz_vsmp(smp)
         endif
         if (cols_per_chnk_mult > 1) then
            ! Then increase (if necessary) to be a multiple of cols_per_chnk_mult
            pcols_vsmp(smp) = &
               cols_per_chnk_mult*ceiling(real(pcols_vsmp(smp),r8)/real(cols_per_chnk_mult,r8))
         endif
         if (cols_per_chnk_max > 0) then
            ! If calculated pcols is too large, recalculate with more chunks per thread
            tmp_chnks_per_thrd = chnks_per_thrd
            do while (pcols_vsmp(smp) > cols_per_chnk_max)
               tmp_chnks_per_thrd = tmp_chnks_per_thrd + 1
               pcols_vsmp(smp) = &
                  ceiling(real(nvsmpcolumns(smp),r8)/real(nvsmpthreads(smp)*tmp_chnks_per_thrd,r8))
               if (opt == 5) then
                  ! pcols should not be (much) larger than the maximum block size
                  if (maxblksiz_vsmp(smp) < pcols_vsmp(smp)) pcols_vsmp(smp) = maxblksiz_vsmp(smp)
               endif
               if ((cols_per_chnk_mult > 1) .and. (cols_per_chnk_mult <= cols_per_chnk_max)) then
                  pcols_vsmp(smp) = &
                     cols_per_chnk_mult*ceiling(real(pcols_vsmp(smp),r8)/real(cols_per_chnk_mult,r8))
               endif
            enddo
         endif
      enddo
   endif
#endif
   max_pcols = maxval(pcols_vsmp)

!
! Options 0-3: split local dynamics blocks into chunks,
!              using wrap-map assignment of columns and
!              day/night and north/south column pairs
!              to chunks to improve load balance
!  Option 0: local is per process
!  Option 1: local is subset of`processes assigned to same SMP node
!  Option 2: local is global
!  Option 3: local is pair of processes chosen to maximize load-balance
!            wrt restriction that only communicate with one other
!            process.
! Option 4: split local dynamics blocks into chunks,
!           using block-map assignment of columns
!             
   if ((opt >= 0) .and. (opt <= 4)) then
!
! Determine number of chunks to keep all threads busy
!
      nchunks = 0
      do smp=0,nvsmp-1
         nvsmpchunks(smp) = nvsmpcolumns(smp)/pcols_vsmp(smp)
         if (mod(nvsmpcolumns(smp), pcols_vsmp(smp)) .ne. 0) then
            nvsmpchunks(smp) = nvsmpchunks(smp) + 1
         endif
         if (nvsmpchunks(smp) < chnks_per_thrd*nvsmpthreads(smp)) then
            nvsmpchunks(smp) = chnks_per_thrd*nvsmpthreads(smp)
         endif
         do while (mod(nvsmpchunks(smp), nvsmpthreads(smp)) .ne. 0)
            nvsmpchunks(smp) = nvsmpchunks(smp) + 1
         enddo
         if (nvsmpchunks(smp) > nvsmpcolumns(smp)) then
            nvsmpchunks(smp) = nvsmpcolumns(smp)
         endif
         nchunks = nchunks + nvsmpchunks(smp)
      enddo      
!
! Determine maximum number of columns to assign to chunks
! in a given SMP
!
      do smp=0,nvsmp-1
         if (nvsmpchunks(smp) /= 0) then
            ntmp1 = nvsmpcolumns(smp)/nvsmpchunks(smp)
            ntmp2 = mod(nvsmpcolumns(smp),nvsmpchunks(smp))
            if (ntmp2 > 0) then
               maxcol_chk(smp) = ntmp1 + 1
               maxcol_chks(smp) = ntmp2
            else
               maxcol_chk(smp) = ntmp1
               maxcol_chks(smp) = nvsmpchunks(smp)
            endif
         else
            maxcol_chk(smp) = 0
            maxcol_chks(smp) = 0
         endif
      enddo
#ifndef PPCOLS
!
! If column cost is provided and pcols is not specified (and opt /= 4),
! then can potentially get better load balanced chunks by picking 
! maxcol_chk, maxcol_chks, and pcols_vsmp as large as possible before
! chunk creation. After the chunks are created, can then set pcols_vsmp
! to the actual sizes generated (per virtual smp). There needs to be
! a bound on this initial max pcols value to guard against memory
! issues. Use cols_per_chnk_max if this is specified. Otherwise just
! stick with the previously calculated maxcol_chk, maxcol_chks, and
! pcols_vsmp and skip the following logic.
! Note that already used the previously calculated pcols_vsmp to determine
! the number of chunks, so this optimization is still constrained by this.
!
      if (cols_per_chnk <= 0) then
         if ((use_cost_d) .and. (cols_per_chnk_max > 0) .and. (opt /= 4)) then
            max_pcols = cols_per_chnk_max
            if ((cols_per_chnk_mult > 1) .and. (cols_per_chnk_mult <= cols_per_chnk_max)) then
!              Decrease (if necessary) to be a multiple of cols_per_chnk_mult
               max_pcols = &
                  cols_per_chnk_mult*floor(real(max_pcols,r8)/real(cols_per_chnk_mult,r8))
            endif
!           If the new max_pcols is larger than the previously calculated pcols_vsmp
!           for a given virtual smp, then set pcols_vsml(smp) to max_pcols.
            do smp=0,nvsmp-1
               if (max_pcols > pcols_vsmp(smp)) then
                  pcols_vsmp(smp) = max_pcols
               endif
            enddo
!           Reset maxcol_chk(smp) to the updated pcols_vsmp(smp), and reset 
!           maxcol_chks(smp) to be the number of processes assigned to the virtual smp.
            maxcol_chk(:) = pcols_vsmp(:)
            maxcol_chks(:) = nvsmpchunks(:)
!           Reset max_pcols. Note this could be larger than necessary (after resetting
!           pcols_vsmp based on the actual number of columns assigned to each chunk), 
!           but it is needed now in order to allocate the gcol array in the
!           chunks data structure. The gcol array will be reallocated after the return
!           from the create_chunks routine, after which there will be no wasted space.
            max_pcols = maxval(pcols_vsmp)
         endif
      endif
#endif
!
! Allocate chunks and knuhcs data structures
!
      allocate( chunks(1:nchunks) )
      do cid=1,nchunks
         allocate( chunks(cid)%gcol(max_pcols) )
      enddo
      allocate( knuhcs(1:ngcols) )
!
! Initialize chunks and knuhcs data structures
!
      chunks(:)%ncols = 0
      chunks(:)%estcost = 0.0_r8
      knuhcs(:)%chunkid = -1
      knuhcs(:)%col = -1
!
! Determine chunk id ranges for each SMP
!
      cid_offset(0) = 1
      local_cid(0) = 0
      do smp=1,nvsmp-1
         cid_offset(smp) = cid_offset(smp-1) + nvsmpchunks(smp-1)
         local_cid(smp) = 0
      enddo    

!
! Determine order in which to traverse columns for assignment to chunks
!
      allocate( cdex(1:ngcols) )

      if ((use_cost_d) .and. (opt < 4)) then
         ! If load balancing using column cost, then sort columns by cost
         ! first, maximum to minimum.
         call IndexSet(ngcols,cdex)
         call IndexSort(ngcols,cdex,cost_d,descend=.true.)
      else
         ! If not using column cost, then sort columns by block ordering,
         ! as done in the original algorithm.
         allocate( udex(1:ngcols) )
         udex(:) = .false.
         i = 0
         do jb=firstblock,lastblock
            blksiz = get_block_gcol_cnt_d(jb)
            call get_block_gcol_d(jb,blksiz,cols_d)

            do ib = 1,blksiz
               curgcol = cols_d(ib)

               ! Record column in cdex in block order if not already
               ! recorded
               if (.not. udex(curgcol)) then
                  i=i+1
                  if (i > ngcols) then
                     if (masterproc) then
                        write(iulog,*) &
                           "PHYS_GRID_INIT error: more dynamics ", &
                           "columns found in block traversal ", &
                           "than expected."
                     endif
                     call endrun()
                  endif
                  cdex(i) = curgcol
                  udex(curgcol) = .true.
               endif

            enddo

         enddo
         deallocate( udex )
      endif

!
! Allocate and initialize min heap data structure, for use in 
! maintaining list of chunks sorted by the assigned computional
! cost (sum of estimated cost for assigned columns)
!
      allocate( heap(1:nchunks) )
      do cid=1,nchunks
         heap(cid) = cid
      enddo

      allocate( heap_len(0:nvsmp-1) )
      do smp=0,nvsmp-1
         heap_len(smp) = nvsmpchunks(smp)
      enddo

!
! Assign columns to chunks
!
      do i=1,ngcols
         curgcol = cdex(i)
         smp = col_vsmp_map(i)

         ! Assign column to a chunk if not already assigned
         if ((dyn_to_latlon_gcol_map(curgcol) .ne. -1) .and. &
             (knuhcs(curgcol)%chunkid == -1)) then

            if ((use_cost_d) .and. (opt < 4)) then

               ! For opt==0,1,2,3 and when using column cost estimates,
               ! add column to chunk with lowest estimated cost chunk
               ! (and with space), i.e. to chunk at root of heap for
               ! current SMP
               if (heap_len(smp) > 0) then
                  cid = heap(cid_offset(smp))
               else
                  if (masterproc) then
                     write(iulog,*) &
                        "PHYS_GRID_INIT error: not enough chunks ", &
                        "or too many columns assigned to current SMP"
                  endif
                  call endrun()
               endif

            else
               ! For opt==4, find next chunk with space
               ! (maxcol_chks > 0 test necessary for opt==4 block map)
               cid = cid_offset(smp) + local_cid(smp)
               if (maxcol_chks(smp) > 0) then
                  do while (chunks(cid)%ncols >=  maxcol_chk(smp))
                     local_cid(smp) = mod(local_cid(smp)+1,nvsmpchunks(smp))
                     cid = cid_offset(smp) + local_cid(smp)
                  enddo
               else
                  do while (chunks(cid)%ncols >=  maxcol_chk(smp)-1)
                     local_cid(smp) = mod(local_cid(smp)+1,nvsmpchunks(smp))
                     cid = cid_offset(smp) + local_cid(smp)
                  enddo
               endif

            endif

            ! Update chunk with new column
            chunks(cid)%ncols = chunks(cid)%ncols + 1
            if (chunks(cid)%ncols .eq. maxcol_chk(smp)) &
               maxcol_chks(smp) = maxcol_chks(smp) - 1

            lcol = chunks(cid)%ncols
            chunks(cid)%gcol(lcol) = curgcol
            chunks(cid)%estcost = chunks(cid)%estcost + cost_d(curgcol)
            knuhcs(curgcol)%chunkid = cid
            knuhcs(curgcol)%col = lcol

            if (opt < 4) then

               ! If space available, look to assign a load-balancing
               ! "twin" to same chunk
               if ( (chunks(cid)%ncols <  maxcol_chk(smp)) .and. &
                    (maxcol_chks(smp) > 0) .and. (twin_alg > 0)) then

                  call find_twin(curgcol, smp, &
                                 proc_vsmp_map, twingcol)

                  if (twingcol > 0) then

                     ! Update chunk with twin column
                     chunks(cid)%ncols = chunks(cid)%ncols + 1
                     if (chunks(cid)%ncols .eq. maxcol_chk(smp)) &
                        maxcol_chks(smp) = maxcol_chks(smp) - 1

                      lcol = chunks(cid)%ncols
                      chunks(cid)%gcol(lcol) = twingcol
                      chunks(cid)%estcost = chunks(cid)%estcost + cost_d(twingcol)
                      knuhcs(twingcol)%chunkid = cid
                      knuhcs(twingcol)%col = lcol
                  endif

               endif

               if (use_cost_d) then

                  ! Re-heapify the min heap
                  call adjust_heap(nchunks, maxcol_chk(smp), &
                                   cid_offset(smp), heap_len(smp), heap)

               else

                  ! Move on to next chunk (wrap map)
                  local_cid(smp) = mod(local_cid(smp)+1,nvsmpchunks(smp))

               endif

            endif

         endif

      enddo

!
! Opt-specific clean up
!
      deallocate( heap_len )
      deallocate( heap     )
      deallocate( cdex     )

   else
!
! Option 5: split individual dynamics blocks into chunks,
!            assigning consecutive columns to the same chunk
!

!
! Determine total number of chunks and
! number of chunks in each "SMP node"
! (assuming no vertical decomposition)
!
      nchunks = 0
      nvsmpchunks(:) = 0
      do j=firstblock,lastblock
         p = get_block_owner_d(j)
         smp = proc_vsmp_map(p)
         blksiz = get_block_gcol_cnt_d(j)
         nlchunks = blksiz/pcols_vsmp(smp)
         if ((pcols_vsmp(smp)*nlchunks) /= blksiz) then
            nlchunks = nlchunks + 1
         endif
         nchunks = nchunks + nlchunks
         nvsmpchunks(smp) = nvsmpchunks(smp) + nlchunks
      enddo
!
! Determine chunk id ranges for each SMP
!
      cid_offset(0) = 1
      local_cid(0) = 0
      do smp=1,nvsmp-1
         cid_offset(smp) = cid_offset(smp-1) + nvsmpchunks(smp-1)
         local_cid(smp) = 0
      enddo
!
! Allocate chunks and knuhcs data structures
!
      allocate( chunks(1:nchunks) )
      do cid=1,nchunks
         allocate( chunks(cid)%gcol(max_pcols) )
      enddo
      allocate( knuhcs(1:ngcols) )
!
! Initialize chunks and knuhcs data structures
!
      chunks(:)%ncols = 0
      chunks(:)%estcost = 0.0_r8
      knuhcs(:)%chunkid = -1
      knuhcs(:)%col = -1
      cid = 0
      do jb=firstblock,lastblock
         p = get_block_owner_d(jb)
         smp = proc_vsmp_map(p)
         blksiz = get_block_gcol_cnt_d(jb)
         call get_block_gcol_d(jb,blksiz,cols_d)

         ib = 0
         do while (ib < blksiz)

            cid = cid_offset(smp) + local_cid(smp)
            max_ncols = min(pcols_vsmp(smp),blksiz-ib)

            ncols = 0
            do i=1,max_ncols
               ib = ib + 1
               ! Check whether global index is for a column that dynamics
               ! intends to pass to the physics
               curgcol = cols_d(ib)
               if (dyn_to_latlon_gcol_map(curgcol) .ne. -1) then
                  ! Yes - then save the information
                  ncols = ncols + 1
                  chunks(cid)%gcol(ncols) = curgcol
                  chunks(cid)%estcost = chunks(cid)%estcost + cost_d(curgcol)
                  knuhcs(curgcol)%chunkid = cid
                  knuhcs(curgcol)%col = ncols
               endif
            enddo
            chunks(cid)%ncols = ncols

            local_cid(smp) = local_cid(smp) + 1
         enddo
      enddo

   endif
!
! Assign chunks to processes.
!
   call assign_chunks(npthreads, nvsmp, proc_vsmp_map, &
                      nvsmpthreads, nvsmpchunks)		      
!
! Save pcols-per-process information, to use in setting pcols 
! (if necessary) and in reporting decomposition information
!
#ifdef PPCOLS
   pcols_proc(:) = pcols
#else
   if (cols_per_chnk > 0) then
      pcols_proc(:) = cols_per_chnk
   else
      ! If pcols is not specified, then calculate pcols_proc based
      ! on the actual size of the chunks that were created.
      pcols_proc(:) = 0
      do cid=1,nchunks
         p = chunks(cid)%owner
         ! Set pcols_proc(p) to the maximum number of columns over chunks
         ! assigned to the process
         if (chunks(cid)%ncols > pcols_proc(p)) then
            pcols_proc(p) = chunks(cid)%ncols
         endif
      enddo
      ! Modify pcols_proc(p), if possible, to be a multiple of
      ! cols_per_chnk_mult. If not already a multiple, increase until it is,
      ! subject to not exceeding cols_per_chnk_max. If this does exceed the
      ! max, then leave pcols_proc(p) unchanged.
      if (cols_per_chnk_mult > 1) then
         do p=0,npes-1
            lmax_pcols = &
               cols_per_chnk_mult*ceiling(real(pcols_proc(p),r8)/real(cols_per_chnk_mult,r8))
            if (cols_per_chnk_max > 0) then
               if (lmax_pcols <= cols_per_chnk_max) then
                  pcols_proc(p) = lmax_pcols
               endif
            else
               pcols_proc(p) = lmax_pcols
            endif
         enddo
      endif
   endif
#endif

!
! Clean up
!
   deallocate( col_vsmp_map   )
   deallocate( nvsmpcolumns   )
   deallocate( nvsmpthreads   )
   deallocate( nvsmpchunks    )
   deallocate( pcols_vsmp     )
   deallocate( maxcol_chk     )
   deallocate( maxcol_chks    )
   deallocate( maxblksiz_vsmp )
   deallocate( cid_offset     )
   deallocate( local_cid      )
   deallocate( cols_d         )
  !deallocate( knuhcs         ) !do not deallocate as it is being used in RRTMG radiation.F90

   return
   end subroutine create_chunks
!
!========================================================================

   subroutine find_partners(opt, proc_busy_d, nvsmp, proc_vsmp_map)
!----------------------------------------------------------------------- 
! 
! Purpose: Divide processes into pairs, attempting to maximize the
!          the number of columns in one process whose twins are in the 
!          other process.
! 
! Method: The day/night and north/south hemisphere complement is defined
!         to be the column twin.
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use dyn_grid, only: get_gcol_block_cnt_d, get_gcol_block_d, &
                       get_block_owner_d
   use pmgrid, only: plev
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: opt           ! chunking option
   logical, intent(in)  :: proc_busy_d(0:npes-1)
                                         ! active/idle dynamics process flags
   integer, intent(out) :: nvsmp         ! calculated number of virtual 
                                         !  SMP nodes
   integer, intent(out) :: proc_vsmp_map(0:npes-1)
                                         ! process/virtual smp map
!---------------------------Local workspace-----------------------------
   integer :: gcol_latlon                ! physics column index (latlon sorted)
   integer :: twingcol_latlon            ! physics column index (latlon sorted)
   integer :: gcol, twingcol             ! physics column indices
   integer :: lon, lat, twinlat          ! longitude and latitude indices
   integer :: twinlon_off                ! estimate as to offset of twinlon
                                         ! on a latitude line
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: jb                         ! block index
   integer :: p, twp                     ! process indices
   integer :: col_proc_map(ngcols)       ! location of columns in 
                                         !  dynamics decomposition
   integer :: twin_proc_map(ngcols)      ! location of column twins in 
                                         !  dynamics decomposition
   integer :: twin_cnt(0:npes-1)         ! for each process, number of twins 
                                         !  in each of the other processes
   logical :: assigned(0:npes-1)         ! flag indicating whether process
                                         !  assigned to an SMP node yet
   integer :: maxpartner, maxcnt         ! process with maximum number of 
                                         !  twins and this count

   logical :: error                      ! error flag 
!-----------------------------------------------------------------------
!
! Determine process location of column and its twin in dynamics decomposition
!
   col_proc_map(:) = -1
   twin_proc_map(:) = -1

   error = .false.
   do gcol_latlon=1,ngcols_p

      ! Assume latitude and longitude symmetries and that index manipulations
      ! are sufficient to find partners. (Will be true for lon/lat grids.)
      gcol = latlon_to_dyn_gcol_map(gcol_latlon)
      lat = lat_p(gcol)
      twinlat = clat_p_tot+1-lat
      lon = lon_p(gcol)
      twinlon_off = mod((lon-1)+(clat_p_cnt(twinlat)/2), clat_p_cnt(twinlat))
      twingcol_latlon = clat_p_idx(twinlat) + twinlon_off
      twingcol = latlon_to_dyn_gcol_map(twingcol_latlon)

      block_cnt = get_gcol_block_cnt_d(gcol)
      call get_gcol_block_d(gcol,block_cnt,blockids,bcids)
      do jb=1,block_cnt
         p = get_block_owner_d(blockids(jb)) 
         if (col_proc_map(gcol) .eq. -1) then
            col_proc_map(gcol) = p
         elseif (col_proc_map(gcol) .ne. p) then
            error = .true.
         endif
      enddo

      block_cnt = get_gcol_block_cnt_d(twingcol)
      call get_gcol_block_d(twingcol,block_cnt,blockids,bcids)
      do jb=1,block_cnt
         p = get_block_owner_d(blockids(jb)) 
         if (twin_proc_map(gcol) .eq. -1) then
            twin_proc_map(gcol) = p
         elseif (twin_proc_map(gcol) .ne. p) then
            error = .true.
         endif
      enddo

   end do

   if (error) then
      if (masterproc) then
         write(iulog,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
            "but vertical decomposition not limited to single process"
      endif
      call endrun()
   endif

!
! Assign process pairs to SMPs, attempting to maximize the number of column,twin
! pairs in same SMP.
!
   assigned(:) = .false.
   twin_cnt(:) = 0
   nvsmp = 0
   do p=0,npes-1
      if ((.not. assigned(p)) .and. (proc_busy_d(p))) then
!
! For each process, determine number of twins in each of the other processes
! (running over all columns multiple times to minimize memory requirements).
!
         do gcol_latlon=1,ngcols_p
            gcol = latlon_to_dyn_gcol_map(gcol_latlon)
            if (col_proc_map(gcol) .eq. p) then
               twin_cnt(twin_proc_map(gcol)) = &
                  twin_cnt(twin_proc_map(gcol)) + 1
            endif
         enddo
!
! Find process with maximum number of twins that has not yet been designated
! a partner.
!
         maxpartner = -1
         maxcnt = 0
         do twp=0,npes-1
            if ((.not. assigned(twp)) .and. (twp .ne. p)) then
               if (twin_cnt(twp) >= maxcnt) then
                  maxcnt = twin_cnt(twp)
                  maxpartner = twp
               endif
            endif
         enddo
!
! Assign p and twp to the same SMP node
!
         if (maxpartner .ne. -1) then
            assigned(p) = .true.
            assigned(maxpartner) = .true.
            proc_vsmp_map(p) = nvsmp
            proc_vsmp_map(maxpartner) = nvsmp
            nvsmp = nvsmp + 1
         else
            if (masterproc) then
               write(iulog,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
                  "but could not divide processes into pairs."
            endif
            call endrun()
         endif
!
      endif
!      
   enddo
!
   return
   end subroutine find_partners
!
!========================================================================

   subroutine find_twin(gcol, smp, proc_vsmp_map, twingcol_f)
!----------------------------------------------------------------------- 
! 
! Purpose: Find column that when paired with gcol in a chunk
!          balances the load. A column is a candidate to be paired with
!          gcol if it is in the same SMP node as gcol as defined
!          by proc_vsmp_map.
! 
! Method: The day/night and north/south hemisphere complement is
!         tried first. If it is not a candidate or if it has already been
!         assigned, then the day/night complement is tried next. If that
!         also is not available, then nothing is returned.
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use dyn_grid, only: get_gcol_block_d, get_block_owner_d

!------------------------------Arguments--------------------------------
   integer, intent(in)  :: gcol          ! global column index for column
                                         ! seeking a twin for
   integer, intent(in)  :: smp           ! index of SMP node 
                                         ! currently assigned to
   integer, intent(in)  :: proc_vsmp_map(0:npes-1)
                                         ! process/virtual smp map
   integer, intent(out) :: twingcol_f
                                         ! global column index for twin
!---------------------------Local workspace-----------------------------
   integer :: lon, lat                   ! global lon/lat indices for column
                                         ! seeking a twin for
   integer :: twinlon, twinlat           ! lon/lat indices of twin candidate
   integer :: twinlon_off                ! estimate as to offset of twinlon
                                         ! on a latitude line
   logical :: found                      ! found flag
   integer :: i                          ! loop index
   integer :: upper, lower               ! search temporaries
   integer :: twingcol_latlon            ! global physics column index (latlon sorted)
   integer :: twingcol_lonlat            ! global physics column index (lonlat sorted)
   integer :: twingcol                   ! global physics column indes
   integer :: diff, min_diff, min_i      ! search temporaries
   integer :: jbtwin(npes)               ! global block indices
   integer :: ibtwin(npes)               ! global column indices
   integer :: twinproc, twinsmp          ! process and smp ids

   real(r8):: twopi                      ! 2*pi
   real(r8):: clat, twinclat             ! latitude and twin
   real(r8):: clon, twinclon             ! longitude and twin

!-----------------------------------------------------------------------
   twingcol_f = -1

!
! Try day/night and north/south hemisphere complement first
!
   ! determine twin latitude
   lat = lat_p(gcol)
   clat = clat_p(lat)
   twinclat = -clat
   twinlat = clat_p_tot+1-lat
   if (clat_p(twinlat) .eq. twinclat) then
      found = .true.
   else
      found = .false.
      upper = twinlat
      lower = twinlat
      if (upper < clat_p_tot) upper = twinlat + 1
      if (lower > 1) lower = twinlat - 1
   endif
   do while (.not. found)
      if      ((abs(clat_p(upper)-twinclat) < abs(clat_p(twinlat)-twinclat)) .and. &
               (upper .ne. twinlat)) then
         twinlat = upper
         if (upper < clat_p_tot) then
            upper = twinlat + 1
         else
            found = .true.
         endif
      else if ((abs(clat_p(lower)-twinclat) < abs(clat_p(twinlat)-twinclat)) .and. &
               (lower .ne. twinlat))    then
         twinlat = lower
         if (lower > 1) then
            lower = twinlat - 1
         else
            found = .true.
         endif
      else
         found = .true.
      endif
    enddo

   ! determine twin longitude
   twopi = 2.0_r8*pi
   lon = lon_p(gcol)
   clon = clon_p(lon)
   twinclon = mod(clon+pi,twopi)
   twinlon = mod((lon-1)+(clon_p_tot/2), clon_p_tot) + 1
   if (clon_p(twinlon) .eq. twinclon) then
      found = .true.
   else
      found = .false.
      upper = twinlon
      lower = twinlon
      if (upper < clon_p_tot) upper = twinlon + 1
      if (lower > 1) lower = twinlon - 1
   endif
   do while (.not. found)
      if      ((abs(clon_p(upper)-twinclon) < abs(clon_p(twinlon)-twinclon)) .and. &
               (upper .ne. twinlon)) then
         twinlon = upper
         if (upper < clon_p_tot) then
            upper = twinlon + 1
         else
            found = .true.
         endif
      else if ((abs(clon_p(lower)-twinclon) < abs(clon_p(twinlon)-twinclon)) .and. &
               (lower .ne. twinlon))    then
         twinlon = lower
         if (lower > 1) then
            lower = twinlon - 1
         else
            found = .true.
         endif
      else
         found = .true.
      endif
   enddo

   ! first, look for an exact match (assuming latitude and longitude symmetries)
   twinlon_off = mod((lon-1)+(clat_p_cnt(twinlat)/2), clat_p_cnt(twinlat))
   twingcol_latlon = clat_p_idx(twinlat) + twinlon_off
   twingcol = latlon_to_dyn_gcol_map(twingcol_latlon)

   ! otherwise, look around for an approximate match using lonlat sorted indices
   if ((lon_p(twingcol) .ne. twinlon) .or. (lat_p(twingcol) .ne. twinlat)) then
      twingcol_lonlat = clon_p_idx(twinlon)
      twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
      min_diff = abs(lat_p(twingcol) - twinlat)
      min_i = 0
      do i = 1, clon_p_cnt(twinlon)-1
         twingcol_lonlat = clon_p_idx(twinlon)+i
         twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
         diff = abs(lat_p(twingcol) - twinlat)
         if (diff < min_diff) then
            min_diff = diff
            min_i = i
         endif
      enddo
      twingcol_lonlat = clon_p_idx(twinlon) + min_i
      twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
   endif

   ! Check whether twin and original are in same smp
   found = .false.
   call get_gcol_block_d(twingcol,npes,jbtwin,ibtwin)
   twinproc = get_block_owner_d(jbtwin(1))
   twinsmp  = proc_vsmp_map(twinproc)
!
   if ((twinsmp .eq. smp) .and. &
       (knuhcs(twingcol)%chunkid == -1)) then
      found = .true.
      twingcol_f = twingcol
   endif
!
! Try day/night complement next
   if (.not. found) then

      ! first, look for an exact match (assuming longitude symmetries)
      twinlon_off = mod((lon-1)+(clat_p_cnt(lat)/2), clat_p_cnt(lat))
      twingcol_latlon = clat_p_idx(lat) + twinlon_off
      twingcol = latlon_to_dyn_gcol_map(twingcol_latlon)

      ! otherwise, look around for an approximate match using lonlat
      ! column ordering
      if ((lon_p(twingcol) .ne. twinlon) .or. &
          (lat_p(twingcol) .ne. lat)) then
         twingcol_lonlat = clon_p_idx(twinlon)
         twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
         min_diff = abs(lat_p(twingcol) - lat)
         min_i = 0
         do i = 1, clon_p_cnt(twinlon)-1
            twingcol_lonlat = clon_p_idx(twinlon)+i
            twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
            diff = abs(lat_p(twingcol) - lat)
            if (diff < min_diff) then
               min_diff = diff
               min_i = i
            endif
         enddo
         twingcol_lonlat = clon_p_idx(twinlon) + min_i
         twingcol = lonlat_to_dyn_gcol_map(twingcol_lonlat)
      endif
!
      call get_gcol_block_d(twingcol,npes,jbtwin,ibtwin)
      twinproc = get_block_owner_d(jbtwin(1))
      twinsmp  = proc_vsmp_map(twinproc)
!
      if ((twinsmp .eq. smp) .and. &
          (knuhcs(twingcol)%chunkid == -1)) then
         found = .true.
         twingcol_f = twingcol
      endif
!
   endif
!
   return
   end subroutine find_twin
!
!========================================================================

   subroutine adjust_heap(nchunks, maxcol_chk, &
                          cid_offset, heap_len, heap)
!----------------------------------------------------------------------- 
! 
! Purpose: Adjust heap after adding columns (and updating the associated
!          computational cost) to the current root, restoring the min 
!          heap property
! 
! Method: Percolate the root down through the heap until find a legal
!         new location. This is a modified version of an algorithm 
!         described in Aho, Hopcroft and Ullman ("The Design and 
!         Analysis of Computer Algorithms") used in sorting. Here it
!         is customized for this slightly different application.
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

   ! size of min heap array (for all virtual SMPs)
   integer, intent(in)    :: nchunks

   ! maximum number of columns assigned to a chunk in current
   ! virtual SMP node
   integer, intent(in)    :: maxcol_chk

   ! beginning of min heap for current virtual SMP in heap array
   integer, intent(in)    :: cid_offset

   ! size of min heap for current virtual SMP (including only chunks
   ! that still have space to be assigned more columns)
   integer, intent(inout) :: heap_len

   ! min heap array used to maintain chunks sorted by assigned 
   ! computational cost
   integer, intent(inout) :: heap(1:nchunks)

!---------------------------Local workspace-----------------------------
   integer :: root                     ! first chunk id in heap (root)
   integer :: heap_last                ! index for last chunk id in heap
   integer :: last_nonleaf             ! index for last non-leaf in heap
   integer :: heap_i, heap_il, heap_ir ! indices used in navigating the 
                                       !  heap
   integer :: parent, lchild, rchild   ! column ids used in maintaining 
                                       !  the heap 
   logical :: done                     ! flag indicating whether heap
                                       !  adjustment is done or not
   real(r8):: min_cost                 ! minimum of two chunk estimated
                                       !  costs
!-----------------------------------------------------------------------
!
   root = heap(cid_offset)
!
   if (chunks(root)%ncols .eq. maxcol_chk) then
! Move chunk to the end of the heap (bringing end to the root)
! and decrement heap length by 1
      heap_last = cid_offset + heap_len - 1
      heap(cid_offset) = heap(heap_last)
      heap(heap_last) = root
      heap_len = heap_len - 1
!
   endif

! Percolate new or updated root to its proper location in the min heap.
! Swapping parent with child if chunk estimated cost is greater than or equal, 
! not just greater than (as in usual heap operations), so that increase
! likelihood of cycling through all chunks
   last_nonleaf = heap_len/2
   heap_i  = 1
   done = .false.
!
   do while (.not. done)
      if (heap_i > last_nonleaf) then
! If no children (a leaf), then done.
         done = .true.
!
      else if (2*heap_i == heap_len) then
! If only a left child, then only a single test.
         heap_il = 2*heap_i
         parent = heap(cid_offset+heap_i-1)
         lchild = heap(cid_offset+heap_il-1)
!
         if (chunks(parent)%estcost >= chunks(lchild)%estcost) then
! If larger or equal, then swap. 
            heap(cid_offset+heap_i-1)  = lchild
            heap(cid_offset+heap_il-1) = parent
            heap_i = heap_il
!
         else
! If smaller, then done.
            done = .true.
!
         endif
!
      else
! If both left and right children, then a number of different possibilities.
         heap_il = 2*heap_i
         heap_ir = heap_il+1
         parent = heap(cid_offset+heap_i-1)
         lchild = heap(cid_offset+heap_il-1)
         rchild = heap(cid_offset+heap_ir-1)
         min_cost = min(chunks(lchild)%estcost,chunks(rchild)%estcost)
!
         if (chunks(parent)%estcost < min_cost) then
! If smaller than both, then done.
            done = .true.
!
         else
            if (chunks(rchild)%estcost > chunks(lchild)%estcost) then
! If rchild has the larger cost, then parent cost must be larger than
! or equal to lchild, so swap with lchild. (For equal lchild and
! rchild values, go with right child since more likely to be a leaf,
! so a little less expensive.)
               heap(cid_offset+heap_i-1)  = lchild
               heap(cid_offset+heap_il-1) = parent
               heap_i = heap_il
!
            else
! If lchild has the larger or equal cost, then parent cost must be 
! larger than or equal to rchild, so swap with rchild. 
               heap(cid_offset+heap_i-1)  = rchild
               heap(cid_offset+heap_ir-1) = parent
               heap_i = heap_ir
!
            endif
!
         endif
!
      endif
!
   enddo

   return
   end subroutine adjust_heap
!
!========================================================================

   subroutine assign_chunks(npthreads, nvsmp, proc_vsmp_map, &
                            nvsmpthreads, nvsmpchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Assign chunks to processes, balancing the number of
!          chunks per thread and minimizing the communication costs
!          in dp_coupling subject to the restraint that columns
!          do not migrate outside of the current SMP node.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plev
   use dyn_grid, only: get_gcol_block_cnt_d, get_gcol_block_d,&
                       get_block_owner_d 
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: npthreads(0:npes-1)
                                         ! number of OpenMP threads per process
   integer, intent(in)  :: nvsmp         ! virtual smp count
   integer, intent(in)  :: proc_vsmp_map(0:npes-1)
                                         ! process/virtual smp map
   integer, intent(in)  :: nvsmpthreads(0:nvsmp-1)
                                         ! number of OpenMP threads 
                                         ! per virtual SMP
   integer, intent(in)  :: nvsmpchunks(0:nvsmp-1)
                                         ! number of chunks assigned 
                                         ! to a given virtual SMP
!---------------------------Local workspace-----------------------------
   integer :: i, jb, p                   ! loop indices
   integer :: cid                        ! chunk id
   integer :: smp                        ! SMP index
   integer :: curgcol                    ! global column index
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: ntsks_vsmp(0:nvsmp-1)      ! number of processes per virtual SMP
   integer :: vsmp_proc_map(max_nproc_vsmp,0:nvsmp-1)   
                                         ! virtual smp to process id map
   integer :: cid_offset(0:nvsmp)        ! chunk id virtual smp offset
   integer :: ntmp1_vsmp(0:nvsmp-1)      ! minimum number of chunks per thread
                                         !  in a virtual SMP
   integer :: ntmp2_vsmp(0:nvsmp-1)      ! number of extra chunks to be assigned
                                         !  in a virtual SMP
   integer :: ntmp3_vsmp(0:nvsmp-1)      ! number of processes in a virtual
                                         !  SMP that get more extra chunks
                                         !  than the others
   integer :: ntmp4_vsmp(0:nvsmp-1)      ! number of extra chunks per process
                                         !  in a virtual SMP
   integer :: ntmp1, ntmp2               ! work variables
!  integer :: npchunks(0:npes-1)         ! number of chunks to be assigned to
!                                        !  a given process
   integer :: cur_npchunks(0:npes-1)     ! current number of chunks assigned 
                                         !  to a given process
   integer :: column_count(0:npes-1)     ! number of columns from current chunk
                                         !  assigned to each process in dynamics
                                         !  decomposition
   integer :: first_nonfull              ! first process (in vsmp_proc_map 
                                         !  ordering) that has room to be assigned
                                         !  another chunk
   integer :: ndyn_task                  ! number of processes in the dynamics 
                                         !  decomposition that were assigned columns
                                         !  in the current chunk
   integer :: dyn_task(npes)             ! list of process ids that were assigned
                                         !  columns in the current chunk in the
                                         !  dynamics decomposition
!-----------------------------------------------------------------------
!
! Count number of processes per virtual SMP and determine virtual SMP
! to process id map
!
   ntsks_vsmp(:) = 0
   vsmp_proc_map(:,:) = -1
   do p=0,npes-1
      smp = proc_vsmp_map(p)
      ntsks_vsmp(smp) = ntsks_vsmp(smp) + 1
      vsmp_proc_map(ntsks_vsmp(smp),smp) = p
   enddo
!
! Determine chunk id ranges for each virtual SMP
!
   cid_offset(0) = 1
   do smp=1,nvsmp
      cid_offset(smp) = cid_offset(smp-1) + nvsmpchunks(smp-1)
   enddo
!
! Determine number of chunks to assign to each process
!
   do smp=0,nvsmp-1
!
! Minimum number of chunks per thread
      ntmp1_vsmp(smp) = nvsmpchunks(smp)/nvsmpthreads(smp)

! Number of extra chunks to be assigned
      ntmp2_vsmp(smp) = mod(nvsmpchunks(smp),nvsmpthreads(smp))

! Number of processes that get more extra chunks than the others
      ntmp3_vsmp(smp) = mod(ntmp2_vsmp(smp),ntsks_vsmp(smp))

! Number of extra chunks per process
      ntmp4_vsmp(smp) = ntmp2_vsmp(smp)/ntsks_vsmp(smp)
      if (ntmp3_vsmp(smp) > 0) then
         ntmp4_vsmp(smp) = ntmp4_vsmp(smp) + 1
      endif
   enddo

   do p=0,npes-1
      smp = proc_vsmp_map(p)

! Update number of extra chunks
      if (ntmp2_vsmp(smp) > ntmp4_vsmp(smp)) then
         ntmp2_vsmp(smp) = ntmp2_vsmp(smp) - ntmp4_vsmp(smp)
      else
         ntmp4_vsmp(smp) = ntmp2_vsmp(smp)
         ntmp2_vsmp(smp) = 0
         ntmp3_vsmp(smp) = 0
      endif

! Set number of chunks
      npchunks(p) = ntmp1_vsmp(smp)*npthreads(p) + ntmp4_vsmp(smp)

! Update extra chunk increment
      if (ntmp3_vsmp(smp) > 0) then
         ntmp3_vsmp(smp) = ntmp3_vsmp(smp) - 1
         if (ntmp3_vsmp(smp) .eq. 0) then
            ntmp4_vsmp(smp) = ntmp4_vsmp(smp) - 1
         endif
      endif
   enddo

!
! Assign chunks to processes: 
!
!  First, initialize number of chunks assigned to each process.
!  Then initialize number of chunk columns assigned to each process 
!  in the dynamics decomposition (column_count). column_count is
!  chunk-specific, and is reset for each chunk in the loop over
!  chunks below, except that '-1' values are retained, where '-1' indicates
!  that the process has already been assigned its quota of chunks.
!
   cur_npchunks(:) = 0
   column_count(:) = 0
!
   do smp=0,nvsmp-1
!
!  Initialize pointer to first process (in vsmp_proc_map ordering) that
!  has room to be assigned another chunk
      first_nonfull = 1
!
      do cid=cid_offset(smp),cid_offset(smp+1)-1
!
!  Determine number of chunk columns assigned to each process in
!  the dynamics decomposition (excepting processes that have already been
!  assigned their quota of chunks). Also build a list of these processes.
         ndyn_task = 0
         do i=1,chunks(cid)%ncols
            curgcol = chunks(cid)%gcol(i)
            block_cnt = get_gcol_block_cnt_d(curgcol)
            call get_gcol_block_d(curgcol,block_cnt,blockids,bcids)
            do jb=1,block_cnt
               p = get_block_owner_d(blockids(jb)) 
               if (column_count(p) > -1) then
                  column_count(p) = column_count(p) + 1
                  if (column_count(p) == 1) then
                     ndyn_task = ndyn_task + 1
                     dyn_task(ndyn_task) = p
                  endif
               endif
            enddo
         enddo
!                              
!  Identify process with most chunk columns in dynamics decomposition,
!  from among those that do not already have their assigned chunk quota 
         ntmp1 = -1
         ntmp2 = -1
         do i=1,ndyn_task
            p = dyn_task(i)
            if (column_count(p) > ntmp1) then
               ntmp1 = column_count(p)
               ntmp2 = p
            endif
         enddo
!
!  If no processes found that qualify, identify some other process that can
!  accept a new chunk
         if (ntmp1 == -1) then
            do i=first_nonfull,ntsks_vsmp(smp)
               p = vsmp_proc_map(i,smp)
               if (column_count(p) /= -1) then
                  ntmp2 = p
                  exit
               endif
            enddo
         endif
!
!  Assign chunk to indicated process, and check whether process now
!  has its quota of chunks
         chunks(cid)%owner   = ntmp2
         cur_npchunks(ntmp2) = cur_npchunks(ntmp2) + 1
         if (cur_npchunks(ntmp2) == npchunks(ntmp2)) then
            column_count(ntmp2) = -1
         endif
!
!  Update total number of columns assigned to this process
         gs_col_num(ntmp2)   = gs_col_num(ntmp2) + chunks(cid)%ncols
!
!  Zero out per process chunk columns counts, but retain -1 value
!  (indicating that process cannot accept any more chunks)
         do i=1,ndyn_task
            p = dyn_task(i)
            if (column_count(p) > 0) then
               column_count(p) = 0
            endif
         enddo
!
!  Update pointer to first nonfull process
         do i=first_nonfull,ntsks_vsmp(smp)
            p = vsmp_proc_map(i,smp)
            if (column_count(p) /= -1) then
               first_nonfull = i
               exit
            endif
         enddo
!
      enddo
!
   enddo
!
   return
   end subroutine assign_chunks
!
!========================================================================

!#######################################################################

end module phys_grid
