! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_fme_horiz_remap
!
!> \brief MPAS ocean FME horizontal remapping wrapper
!> \author FME Training Data Support
!> \date   April 2026
!> \details
!>  Wrapper around shr_horiz_remap_t for MPAS-Ocean. Handles MPAS cell
!>  decomposition mapping, send buffer packing, CRS SpMV remapping, and
!>  direct PIO output of remapped fields on a lat-lon grid.
!>
!>  The PIO output bypasses the MPAS stream manager because streams
!>  only support the native unstructured grid (nCells). Remapped output
!>  uses (nlon, nlat) dimensions directly.
!>
!>  Usage:
!>    call ocn_fme_horiz_remap_init(domain, mapfile, err)
!>    call ocn_fme_horiz_remap_init_time_epoch(domain, amName, err)
!>    call ocn_fme_horiz_remap_field(domain, field2d, numlev, fld_out, err)
!>    call ocn_fme_horiz_remap_open_file(domain, filename, err)
!>    call ocn_fme_horiz_remap_def_var(varname, units, long_name, numlev, err)
!>    call ocn_fme_horiz_remap_write_var(varname, fld_out, numlev, err)
!>    call ocn_fme_horiz_remap_write_time(time_lo, time_hi, err)
!>    call ocn_fme_horiz_remap_close_file(err)
!>    logical = ocn_fme_horiz_remap_active()
!>
!-----------------------------------------------------------------------

module ocn_fme_horiz_remap

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_kind_types
   use mpas_timekeeping
   use mpas_stream_manager, only: MPAS_stream_mgr_get_property
   use mpas_log, only: mpas_log_write

   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_horiz_remap_mod, only: shr_horiz_remap_t, SHR_FILL_VALUE
   use pio, only: file_desc_t, var_desc_t, io_desc_t, iosystem_desc_t, &
                  pio_openfile, pio_closefile, pio_createfile, &
                  pio_def_dim, pio_def_var, pio_enddef, &
                  pio_put_att, pio_put_var, pio_write_darray, pio_read_darray, &
                  pio_get_att, pio_get_var, &
                  pio_inq_dimid, pio_inq_varid, pio_inq_dimlen, &
                  pio_initdecomp, pio_freedecomp, pio_setframe, &
                  PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF4P, &
                  PIO_CLOBBER, PIO_WRITE, PIO_REAL, PIO_DOUBLE, PIO_INT, &
                  PIO_OFFSET_KIND, PIO_GLOBAL, PIO_NOERR, &
                  PIO_UNLIMITED

   implicit none
   private
   save

   public :: ocn_fme_horiz_remap_init
   public :: ocn_fme_horiz_remap_init_time_epoch
   public :: ocn_fme_horiz_remap_field
   public :: ocn_fme_horiz_remap_active
   public :: ocn_fme_horiz_remap_get_dims
   public :: ocn_fme_remap_file_t
   public :: ocn_fme_remap_build_filename
   public :: ocn_fme_horiz_remap_write_accum
   public :: ocn_fme_horiz_remap_read_accum

   type(shr_horiz_remap_t) :: remap_data

   ! Send-side cell mapping
   integer, allocatable :: send_local_cell(:)
   real(r8), allocatable :: ws_send_buf(:)
   logical :: is_initialized = .false.

   ! Cached at init from config_do_restart. Append-mode reopen of the
   ! `.remapped.nc` files (and accumulator-sidecar reads) only fire on
   ! warm restart. On a cold start we always clobber-create, even if an
   ! older `.remapped.nc` is still in the run dir from a prior CIME
   ! compare-copy that didn't get cleared.
   logical :: is_restart_run = .false.

   ! CF time-coordinate epoch shared across all FME output files for
   ! this component. Auto-discovered at AM init via init_time_epoch:
   !   1. If a `<amName>.fme_accum_restart.nc` sidecar exists (warm
   !      restart), read its `time_reference_date` global attribute.
   !   2. Otherwise capture `MPAS_START_TIME` from the MPAS clock. On
   !      cold start MPAS_START_TIME == CIME RUN_STARTDATE (verified at
   !      ocn_comp_mct.F:611-614 -- the driver sets MPAS_START_TIME from
   !      seq_timemgr_EClockGetData(EClock, ECurrTime=...) for
   !      runtype=='initial').
   ! On warm restart, sidecar persistence wins so the epoch never
   ! drifts to leg-start-time even though MPAS_START_TIME gets reset on
   ! each restart leg. Each *_write_accum call rewrites the current
   ! time_epoch_str into the sidecar so subsequent restarts can recover
   ! it. All FME output files write time(time) as days since this epoch
   ! and stamp time:units = "days since {time_epoch_cf}".
   type(MPAS_Time_Type) :: time_epoch
   character(len=64)    :: time_epoch_str = ''   ! `YYYY-MM-DD_HH:MM:SS` (MPAS form)
   character(len=64)    :: time_epoch_cf  = ''   ! `YYYY-MM-DD HH:MM:SS` (CF form, space)
   logical              :: time_epoch_set = .false.

   ! PIO decomposition cache (shared across all file instances)
   type(iosystem_desc_t), pointer :: pio_sys => null()
   type(io_desc_t) :: iodesc_2d        ! PIO_REAL  (output records)
   type(io_desc_t) :: iodesc_2d_dbl    ! PIO_DOUBLE (accumulator restart sidecars)
   logical :: iodesc_2d_valid = .false.

   ! MPAS-O fill value for input field detection (matches def_var _FillValue)
   real(r8), parameter :: MPAS_OCN_FILL = SHR_FILL_VALUE

   ! Per-file output state: each AM creates its own instance
   integer, parameter :: max_pio_vars = 512

   type :: ocn_fme_remap_file_t
      type(file_desc_t) :: pio_file
      logical :: file_is_open = .false.
      integer :: dim_lon_id = 0, dim_lat_id = 0, dim_time_id = 0
      integer :: dim_nbnd_id = 0
      integer :: time_record = 0
      integer :: n_pio_vars = 0
      character(len=64)  :: pio_var_names(max_pio_vars)
      type(var_desc_t)   :: pio_var_descs(max_pio_vars)
      ! Per-active-def is_static flag, parallel to pio_var_descs. True
      ! for vars defined with shape (lon, lat) only — no time dim, written
      ! once per file. Used by write_var to skip pio_setframe.
      logical            :: pio_var_is_static(max_pio_vars) = .false.
      type(var_desc_t)   :: pio_lon_var, pio_lat_var
      ! CF time variables (replace legacy MPAS xtime). See module-scope
      ! comment on time_epoch for epoch handling.
      type(var_desc_t)   :: pio_time_var, pio_time_bnds_var
      type(var_desc_t)   :: pio_date_var, pio_datesec_var
      ! Stored metadata for replaying var definitions after rotation
      character(len=64)  :: stored_units(max_pio_vars)
      character(len=128) :: stored_long_names(max_pio_vars)
      ! Per-registered-var is_static flag, replayed by check_rotate so
      ! every fresh file's def_var calls reproduce the original shape.
      logical            :: stored_is_static(max_pio_vars) = .false.
      integer :: n_registered_vars = 0
      ! Current open filename and stream name for rotation
      character(len=512) :: current_filename = ''
      character(len=128) :: stream_name = ''
      ! True when the file was reopened from disk (append mode, e.g. after
      ! restart) rather than freshly created. In this state lon/lat are
      ! already populated and var_descs are looked up via pio_inq_varid
      ! instead of being defined.
      logical :: existing_file_reopened = .false.
   contains
      procedure :: open_file    => remap_file_open
      procedure :: def_var      => remap_file_def_var
      procedure :: register_var => remap_file_register_var
      procedure :: write_var    => remap_file_write_var
      procedure :: write_time   => remap_file_write_time
      procedure :: close_file   => remap_file_close
      procedure :: check_rotate => remap_file_check_rotate
   end type ocn_fme_remap_file_t

!***********************************************************************

contains

!***********************************************************************

   subroutine ocn_fme_horiz_remap_init(domain, mapfile, err)!{{{
      use mpi, only: MPI_INTEGER, MPI_MAX

      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: mapfile
      integer, intent(out) :: err

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: meshPool
      integer, pointer :: nCellsSolve
      integer, dimension(:), pointer :: indexToCellID
      integer :: myrank, nprocs, comm, ierr
      integer :: iCell, gcol, map_col
      integer, allocatable :: gcol_to_rank(:), gcol_to_rank_global(:)
      integer, allocatable :: send_gcol_list(:)
      integer :: i
      integer :: N_total, ocean_count
      integer, allocatable :: is_ocean(:), gcol_to_ocean_col(:)
      integer, dimension(:), pointer :: maxLevelCell

      err = 0

      if (is_initialized) return

      if (len_trim(mapfile) == 0 .or. trim(mapfile) == 'none') then
         is_initialized = .false.
         return
      end if

      myrank = domain % dminfo % my_proc_id
      nprocs = domain % dminfo % nprocs
      comm   = domain % dminfo % comm

      ! Phase 1: Read map file
      call remap_data%read_mapfile(trim(mapfile), comm, myrank, nprocs, &
           domain % ioContext % pio_iosystem, ierr)
      if (ierr /= 0) then
         call mpas_log_write('ocn_fme_horiz_remap_init: failed to read map file: ' &
              // trim(mapfile), messageType=MPAS_LOG_ERR)
         err = ierr
         return
      end if

      ! Resolve cell ID mismatch: map files generated from masked SCRIP
      ! meshes use contiguous 1-based indexing over ocean-only cells
      ! (1..N_ocean), whereas MPAS indexToCellID covers all cells
      ! (1..N_total including land). Build a mapping from global cell ID
      ! to ocean-only index so gcol_to_rank uses the map file's numbering.

      ! Find global N_total (max indexToCellID across all ranks)
      N_total = 0
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
         do iCell = 1, nCellsSolve
            if (indexToCellID(iCell) > N_total) N_total = indexToCellID(iCell)
         end do
         block => block % next
      end do
      call mpi_allreduce(N_total, iCell, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
      N_total = iCell

      allocate(gcol_to_ocean_col(N_total))
      if (remap_data%n_a == N_total) then
         ! Map file covers all cells — identity mapping
         do i = 1, N_total
            gcol_to_ocean_col(i) = i
         end do
      else
         ! Map file was generated from masked mesh (ocean only).
         ! Build ocean mask via maxLevelCell and compute prefix-sum mapping.
         allocate(is_ocean(N_total))
         is_ocean(:) = 0
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
            do iCell = 1, nCellsSolve
               if (maxLevelCell(iCell) > 0) then
                  is_ocean(indexToCellID(iCell)) = 1
               end if
            end do
            block => block % next
         end do

         ! Gather ocean mask across all ranks
         allocate(gcol_to_rank_global(N_total))
         call mpi_allreduce(is_ocean, gcol_to_rank_global, N_total, &
              MPI_INTEGER, MPI_MAX, comm, ierr)
         is_ocean = gcol_to_rank_global
         deallocate(gcol_to_rank_global)

         ! Prefix sum: assign contiguous 1-based index to ocean cells
         ocean_count = 0
         do i = 1, N_total
            if (is_ocean(i) == 1) then
               ocean_count = ocean_count + 1
               gcol_to_ocean_col(i) = ocean_count
            else
               gcol_to_ocean_col(i) = 0
            end if
         end do
         deallocate(is_ocean)

         if (myrank == 0) then
            call mpas_log_write('ocn_fme_horiz_remap_init: N_total=$i, n_a=$i, ocean_count=$i', &
                 intArgs=(/N_total, remap_data%n_a, ocean_count/))
         end if
      end if

      ! Build gcol_to_rank mapping using map-file column indices
      allocate(gcol_to_rank(remap_data%n_a))
      gcol_to_rank(:) = -1

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

         do iCell = 1, nCellsSolve
            map_col = gcol_to_ocean_col(indexToCellID(iCell))
            if (map_col >= 1 .and. map_col <= remap_data%n_a) then
               gcol_to_rank(map_col) = myrank
            end if
         end do

         block => block % next
      end do

      ! Gather full gcol_to_rank across all ranks
      allocate(gcol_to_rank_global(remap_data%n_a))
      call mpi_allreduce(gcol_to_rank, gcol_to_rank_global, remap_data%n_a, &
           MPI_INTEGER, MPI_MAX, comm, ierr)
      gcol_to_rank = gcol_to_rank_global
      deallocate(gcol_to_rank_global)

      ! Phase 2: Build MPI communication pattern
      call remap_data%build_comm(gcol_to_rank, comm, myrank, nprocs, &
           send_gcol_list, ierr)
      deallocate(gcol_to_rank)

      if (ierr /= 0) then
         call mpas_log_write('ocn_fme_horiz_remap_init: failed to build comm pattern', &
              messageType=MPAS_LOG_ERR)
         err = ierr
         return
      end if

      ! Convert send_gcol_list to local cell indices
      allocate(send_local_cell(remap_data%n_send_total))
      send_local_cell(:) = 0
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

         do i = 1, remap_data%n_send_total
            if (send_local_cell(i) > 0) cycle  ! already found
            do iCell = 1, nCellsSolve
               map_col = gcol_to_ocean_col(indexToCellID(iCell))
               if (map_col == send_gcol_list(i)) then
                  send_local_cell(i) = iCell
                  exit
               end if
            end do
         end do

         block => block % next
      end do

      if (allocated(send_gcol_list)) deallocate(send_gcol_list)
      if (allocated(gcol_to_ocean_col)) deallocate(gcol_to_ocean_col)

      ! Validate that all needed source cells were found locally.
      ! A zero entry means the map file references a global cell not in any
      ! local block, which would cause out-of-bounds array access during remap.
      do i = 1, remap_data%n_send_total
         if (send_local_cell(i) == 0) then
            call mpas_log_write('ocn_fme_horiz_remap_init: FATAL - source cell $i ' // &
                 'not found in local blocks. Map file may be incompatible with mesh.', &
                 messageType=MPAS_LOG_CRIT, intArgs=(/i/))
            err = 1
            return
         end if
      end do

      ! Store PIO system reference for later file creation
      pio_sys => domain % ioContext % pio_iosystem

      ! Build 2D PIO decomposition for lat-lon output
      call build_iodesc_2d()

      ! Cache config_do_restart so append-mode logic doesn't need to
      ! re-query the config pool on every check_rotate call.
      block
         logical, pointer :: config_do_restart
         call mpas_pool_get_config(domain % configs, 'config_do_restart', &
              config_do_restart)
         if (associated(config_do_restart)) then
            is_restart_run = config_do_restart
         else
            is_restart_run = .false.
         end if
      end block

      is_initialized = .true.

      if (myrank == 0) then
         call mpas_log_write('ocn_fme_horiz_remap_init: initialized with ' &
              // 'n_a=$i, n_b=$i, nlat=$i, nlon=$i', &
              intArgs=(/remap_data%n_a, remap_data%n_b, &
                        remap_data%nlat, remap_data%nlon/))
      end if

   end subroutine ocn_fme_horiz_remap_init!}}}

!***********************************************************************

   subroutine build_iodesc_2d()!{{{
      integer(PIO_OFFSET_KIND), allocatable :: idof(:)
      integer :: i, global_row, ilon, ilat, nlat, nlon

      nlat = remap_data%nlat
      nlon = remap_data%nlon

      allocate(idof(remap_data%n_b_local))
      do i = 1, remap_data%n_b_local
         global_row = remap_data%row_start + i - 1
         ilon = mod(global_row - 1, nlon) + 1
         ilat = (global_row - 1) / nlon + 1
         idof(i) = int(ilon, PIO_OFFSET_KIND) + &
              int(nlon, PIO_OFFSET_KIND) * int(ilat - 1, PIO_OFFSET_KIND)
      end do

      call pio_initdecomp(pio_sys, PIO_REAL,   (/nlon, nlat/), idof, iodesc_2d)
      call pio_initdecomp(pio_sys, PIO_DOUBLE, (/nlon, nlat/), idof, iodesc_2d_dbl)
      deallocate(idof)
      iodesc_2d_valid = .true.

   end subroutine build_iodesc_2d!}}}

!***********************************************************************

   subroutine ocn_fme_horiz_remap_field(domain, field2d, numlev, fld_out, err, &
        mask)!{{{
   !> Remap a field from the MPAS unstructured grid to the lat-lon target grid.
   !> Builds a validity mask and uses apply_masked so that:
   !>   - land cells (mask==0) and fill values (MPAS_OCN_FILL) are excluded
   !>   - target cells with no valid sources get SHR_FILL_VALUE
   !>   - coastal/mixed cells are normalized by actual valid coverage
   !>
   !> field2d(numlev, *) is on the native MPAS grid indexed by cell.
   !> If mask(nCells) is provided, cells where mask(iCell)==0 are land.

      type(domain_type), intent(inout) :: domain
      integer, intent(in) :: numlev
      real(kind=RKIND), intent(in) :: field2d(numlev, *)
      real(r8), intent(out) :: fld_out(remap_data%n_b_local, numlev)
      integer, intent(out) :: err
      integer, intent(in), optional :: mask(*)

      integer :: i, k, needed, comm, nprocs, icell, nlev_packed
      real(r8) :: val
      logical :: cell_valid

      ! Fill detection threshold: 10% of MPAS fill value.
      ! Catches 1e34 fills and any floating-point variants while staying
      ! far above any physical value (geophysical quantities < 1e10).
      real(r8), parameter :: fill_detect = MPAS_OCN_FILL * 0.1_r8

      err = 0

      if (.not. is_initialized) then
         err = 1
         return
      end if

      nprocs = domain % dminfo % nprocs
      comm   = domain % dminfo % comm

      nlev_packed = numlev + 1

      needed = remap_data%n_send_total * nlev_packed
      if (.not. allocated(ws_send_buf) .or. size(ws_send_buf) < max(1, needed)) then
         if (allocated(ws_send_buf)) deallocate(ws_send_buf)
         allocate(ws_send_buf(max(1, needed)))
      end if

      if (remap_data%n_send_total > 0) then
         do i = 1, remap_data%n_send_total
            icell = send_local_cell(i)

            ! Determine cell validity from mask argument (maxLevelCell)
            if (present(mask)) then
               cell_valid = (mask(icell) >= 1)
            else
               cell_valid = .true.
            end if

            if (.not. cell_valid) then
               ! Land cell: zero field, mask = 0
               do k = 1, numlev
                  ws_send_buf((i-1)*nlev_packed + k) = 0.0_r8
               end do
               ws_send_buf((i-1)*nlev_packed + nlev_packed) = 0.0_r8
            else
               ! Ocean cell: pack field values, detect fills per level
               do k = 1, numlev
                  ! Below-bottom check (for multi-level fields)
                  if (present(mask) .and. mask(icell) < k) then
                     ws_send_buf((i-1)*nlev_packed + k) = 0.0_r8
                     cycle
                  end if
                  val = real(field2d(k, icell), r8)
                  if (abs(val) >= fill_detect) then
                     ws_send_buf((i-1)*nlev_packed + k) = 0.0_r8
                  else
                     ws_send_buf((i-1)*nlev_packed + k) = val
                  end if
               end do
               ! Mask = 1.0 for ocean cells
               ws_send_buf((i-1)*nlev_packed + nlev_packed) = 1.0_r8
            end if
         end do
      end if

      call remap_data%apply_masked(ws_send_buf, numlev, fld_out, comm, nprocs, err)

   end subroutine ocn_fme_horiz_remap_field!}}}

!***********************************************************************

   !-----------------------------------------------------------------
   ! CF time-coordinate epoch
   !-----------------------------------------------------------------

   subroutine ocn_fme_horiz_remap_init_time_epoch(domain, amName, err)!{{{
   !> Auto-discover the CF time-coordinate epoch and populate module state.
   !> Idempotent across the four FME AMs: only the first caller does work;
   !> later AMs see time_epoch_set=.true. and no-op.
   !>
   !> Discovery order (sidecar wins so warm-restart legs preserve the
   !> original cold-start epoch):
   !>   1. Try to read time_reference_date global attribute from
   !>      `<amName>.fme_accum_restart.nc`. Present iff the sidecar was
   !>      written by a prior leg with the new format.
   !>   2. Fall back to MPAS_START_TIME from the MPAS clock. On cold
   !>      start this equals CIME RUN_STARTDATE (driver sets it from
   !>      EClock for runtype='initial', see ocn_comp_mct.F:611-614).
   !>      On warm restart with no sidecar fallback (rare), this would
   !>      be the leg-start time -- the wrong value. Logged as WARN.
      type(domain_type),    intent(inout) :: domain
      character(len=*),     intent(in)    :: amName
      integer,              intent(out)   :: err
      character(len=64)    :: epoch_str
      character(len=StrKIND) :: epoch_str_buf  ! mpas_get_time requires StrKIND
      logical              :: from_sidecar
      type(MPAS_Time_Type) :: start_time
      integer              :: ierr

      err = 0
      if (time_epoch_set) return

      epoch_str = ''
      from_sidecar = .false.
      call try_read_epoch_from_sidecar(amName, epoch_str, from_sidecar)

      if (.not. from_sidecar) then
         start_time = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr)
         call mpas_get_time(start_time, dateTimeString=epoch_str_buf, ierr=ierr)
         epoch_str = trim(epoch_str_buf)
         if (ierr /= 0) then
            call mpas_log_write('ocn_fme_horiz_remap_init_time_epoch: failed ' // &
                 'to read MPAS_START_TIME from clock', messageType=MPAS_LOG_CRIT)
            err = ierr
            return
         end if
         if (is_restart_run) then
            call mpas_log_write('ocn_fme_horiz_remap_init_time_epoch: warm ' // &
                 'restart but no sidecar epoch found; falling back to ' // &
                 'MPAS_START_TIME = "' // trim(epoch_str) // '". This is ' // &
                 'the LEG start time, not the case start -- the resulting ' // &
                 '`time:units` may not align with EAM. Investigate why ' // &
                 'sidecar `' // trim(amName) // '.fme_accum_restart.nc` ' // &
                 'lacks `time_reference_date`.', messageType=MPAS_LOG_WARN)
         end if
      end if

      call set_time_epoch(epoch_str, err)

      if (err == 0) then
         if (from_sidecar) then
            call mpas_log_write('ocn_fme_horiz_remap_init_time_epoch: epoch ' // &
                 'restored from sidecar = "' // trim(time_epoch_cf) // '"')
         else
            call mpas_log_write('ocn_fme_horiz_remap_init_time_epoch: epoch ' // &
                 'captured from MPAS_START_TIME = "' // trim(time_epoch_cf) // '"')
         end if
      end if

   end subroutine ocn_fme_horiz_remap_init_time_epoch!}}}

   subroutine try_read_epoch_from_sidecar(amName, epoch_str, success)!{{{
   !> Read the time_reference_date global attribute from an existing
   !> accumulator sidecar file (if any). Returns success=.false. when
   !> the sidecar is absent or lacks the attribute.
      character(len=*), intent(in)  :: amName
      character(len=*), intent(out) :: epoch_str
      logical,          intent(out) :: success
      character(len=512)            :: filename
      type(file_desc_t)             :: pio_acc_file
      logical                       :: file_exists
      integer                       :: ierr

      success = .false.
      epoch_str = ''
      if (.not. is_initialized) return

      call fme_accum_filename(trim(amName), filename)
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) return

      ierr = pio_openfile(pio_sys, pio_acc_file, PIO_IOTYPE_PNETCDF, &
           trim(filename), PIO_WRITE)
      if (ierr /= PIO_NOERR) return

      ierr = pio_get_att(pio_acc_file, PIO_GLOBAL, 'time_reference_date', &
           epoch_str)
      call pio_closefile(pio_acc_file)
      if (ierr == PIO_NOERR .and. len_trim(epoch_str) > 0) then
         success = .true.
      end if

   end subroutine try_read_epoch_from_sidecar!}}}

   subroutine set_time_epoch(epoch_str, err)!{{{
   !> Internal helper. Stores the epoch in module state (string, MPAS
   !> Time_Type, and CF-form string with space separator). Validates
   !> against any previously-set value; mismatch is a hard error so two
   !> AMs can't disagree on the epoch.
      character(len=*), intent(in) :: epoch_str
      integer,          intent(out) :: err
      integer :: ierr, idx
      character(len=64) :: cf_str

      err = 0

      if (len_trim(epoch_str) == 0) then
         call mpas_log_write('set_time_epoch: empty epoch string', &
              messageType=MPAS_LOG_CRIT)
         err = 1; return
      end if

      if (time_epoch_set) then
         if (trim(epoch_str) /= trim(time_epoch_str)) then
            call mpas_log_write('set_time_epoch: epoch already set to "' // &
                 trim(time_epoch_str) // '"; rejecting reset to "' // &
                 trim(epoch_str) // '". Did sidecars get out of sync?', &
                 messageType=MPAS_LOG_CRIT)
            err = 2
         end if
         return
      end if

      call mpas_set_time(time_epoch, dateTimeString=trim(epoch_str), ierr=ierr)
      if (ierr /= 0) then
         call mpas_log_write('set_time_epoch: failed to parse "' // &
              trim(epoch_str) // '" as YYYY-MM-DD_HH:MM:SS', &
              messageType=MPAS_LOG_CRIT)
         err = ierr; return
      end if

      ! CF requires `days since YYYY-MM-DD HH:MM:SS` with a space, not '_'.
      cf_str = trim(epoch_str)
      idx = index(cf_str, '_')
      if (idx > 0) cf_str(idx:idx) = ' '

      time_epoch_str = trim(epoch_str)
      time_epoch_cf  = trim(cf_str)
      time_epoch_set = .true.
   end subroutine set_time_epoch!}}}

   subroutine time_to_days_since_epoch(t, days, err)!{{{
   !> Convert an MPAS_Time_Type to days since the configured epoch.
   !> Uses ESMF time-interval arithmetic via mpas_get_timeInterval so
   !> the calendar (noleap) is respected exactly.
      type(MPAS_Time_Type), intent(in) :: t
      real(r8), intent(out) :: days
      integer, intent(out) :: err
      type(MPAS_TimeInterval_Type) :: dt_interval
      real(kind=RKIND) :: dt_seconds
      integer :: ierr

      days = 0.0_r8
      err  = 0
      if (.not. time_epoch_set) then
         call mpas_log_write('time_to_days_since_epoch: epoch not set', &
              messageType=MPAS_LOG_ERR)
         err = 1
         return
      end if

      dt_interval = t - time_epoch
      call mpas_get_timeInterval(dt_interval, StartTimeIn=time_epoch, &
           dt=dt_seconds, ierr=ierr)
      if (ierr /= 0) then
         err = ierr
         return
      end if
      days = real(dt_seconds, r8) / 86400.0_r8

   end subroutine time_to_days_since_epoch!}}}

   !-----------------------------------------------------------------
   ! Type-bound methods for per-file PIO output
   !-----------------------------------------------------------------

   subroutine remap_file_open(self, filename, domain, err)!{{{
   !> Open the remapped output file. If a file with this name already
   !> exists on disk (warm restart reopening leg-1's monthly file),
   !> open it in append mode:
   !>   - skip dim/var define; look up existing IDs via pio_inq_*
   !>   - read existing numeric `time` values and seed time_record to
   !>     the count of frames whose `time < MPAS_NOW (in days)`. The
   !>     next write OVERWRITES the first frame whose time >= restart
   !>     time (rather than appending after it), which is the correct
   !>     semantics for ERS-style tests where leg-1 already wrote
   !>     records past the restart point: leg-2 stomps on those
   !>     stale-from-leg-2's-perspective records with its own values.
   !>     For BFB the overwritten value matches what was there.
   !> Otherwise create-new in PIO_CLOBBER mode (the cold-start path).
      class(ocn_fme_remap_file_t), intent(inout) :: self
      character(len=*), intent(in) :: filename
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      integer :: nlat, nlon, ierr, i, count_le_now
      integer(PIO_OFFSET_KIND) :: time_len
      logical :: file_exists
      type(MPAS_Time_Type) :: now_time
      character(len=StrKIND) :: now_str
      real(r8) :: now_days
      real(r8), allocatable :: existing_time(:)

      err = 0
      if (.not. is_initialized) then; err = 1; return; end if
      if (.not. time_epoch_set) then
         call mpas_log_write('ocn_fme_remap_file: time epoch not set; ' // &
              'AM init must call ocn_fme_horiz_remap_set_time_epoch first', &
              messageType=MPAS_LOG_CRIT)
         err = 1; return
      end if

      nlat = remap_data%nlat
      nlon = remap_data%nlon

      inquire(file=trim(filename), exist=file_exists)

      ! Append-mode reopen ONLY on warm restart. On cold start always
      ! clobber-create — protects against stale `.remapped.nc` files
      ! left behind in the run dir by a prior CIME compare-copy or
      ! aborted run.
      if (file_exists .and. is_restart_run) then
         ! --- Append mode: reopen existing file and continue writing -------
         ierr = pio_openfile(pio_sys, self%pio_file, PIO_IOTYPE_PNETCDF, &
              trim(filename), PIO_WRITE)
         if (ierr /= PIO_NOERR) then
            call mpas_log_write('ocn_fme_remap_file: failed to reopen ' &
                 // trim(filename), messageType=MPAS_LOG_ERR)
            err = ierr; return
         end if

         err = 0
         call check_pio('inq_dimid lon', &
              pio_inq_dimid(self%pio_file, 'lon', self%dim_lon_id), err)
         call check_pio('inq_dimid lat', &
              pio_inq_dimid(self%pio_file, 'lat', self%dim_lat_id), err)
         call check_pio('inq_dimid time', &
              pio_inq_dimid(self%pio_file, 'time', self%dim_time_id), err)
         call check_pio('inq_dimid nbnd', &
              pio_inq_dimid(self%pio_file, 'nbnd', self%dim_nbnd_id), err)

         call check_pio('inq_varid lon',       pio_inq_varid(self%pio_file, &
              'lon',       self%pio_lon_var),       err)
         call check_pio('inq_varid lat',       pio_inq_varid(self%pio_file, &
              'lat',       self%pio_lat_var),       err)
         call check_pio('inq_varid time',      pio_inq_varid(self%pio_file, &
              'time',      self%pio_time_var),      err)
         call check_pio('inq_varid time_bnds', pio_inq_varid(self%pio_file, &
              'time_bnds', self%pio_time_bnds_var), err)
         call check_pio('inq_varid date',      pio_inq_varid(self%pio_file, &
              'date',      self%pio_date_var),      err)
         call check_pio('inq_varid datesec',   pio_inq_varid(self%pio_file, &
              'datesec',   self%pio_datesec_var),   err)

         call check_pio('inq_dimlen time', &
              pio_inq_dimlen(self%pio_file, self%dim_time_id, time_len), err)

         ! Read existing numeric times and seed time_record at the largest
         ! frame index whose time is STRICTLY less than the current model
         ! time (in days since epoch). This makes leg-2's next write_time
         ! land on the FIRST frame whose time >= restart point — i.e., it
         ! overwrites stale records that leg-1 wrote past the restart.
         now_time = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
         call mpas_get_time(now_time, dateTimeString=now_str, ierr=ierr)
         call time_to_days_since_epoch(now_time, now_days, ierr)

         allocate(existing_time(int(time_len)))
         ierr = pio_get_var(self%pio_file, self%pio_time_var, existing_time)
         if (ierr /= PIO_NOERR) then
            call mpas_log_write('ocn_fme_remap_file: failed to read time ' // &
                 'from ' // trim(filename) // '; falling back to N_existing', &
                 messageType=MPAS_LOG_WARN)
            self%time_record = int(time_len)
         else
            count_le_now = 0
            do i = 1, int(time_len)
               if (existing_time(i) < now_days) then
                  count_le_now = count_le_now + 1
               end if
            end do
            self%time_record = count_le_now
         end if
         deallocate(existing_time)

         self%n_pio_vars = 0
         self%file_is_open = .true.
         self%existing_file_reopened = .true.

         if (err /= 0) then
            call mpas_log_write('ocn_fme_remap_file: append-mode setup ' // &
                 'failed for ' // trim(filename), messageType=MPAS_LOG_ERR)
         else
            call mpas_log_write('ocn_fme_remap_file: reopened ' // &
                 trim(filename) // ' in append mode; existing time=$i, ' // &
                 'next write at frame=$i+1 (MPAS_NOW=' // trim(now_str) // ')', &
                 intArgs=(/int(time_len), self%time_record/))
         end if
         return
      end if

      ! --- Fresh create -------------------------------------------------
      ierr = pio_createfile(pio_sys, self%pio_file, PIO_IOTYPE_PNETCDF, &
           trim(filename), PIO_CLOBBER)
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('ocn_fme_remap_file: failed to create ' // trim(filename), &
              messageType=MPAS_LOG_ERR)
         err = ierr; return
      end if

      err = 0
      call check_pio('def_dim lon',  pio_def_dim(self%pio_file, 'lon',  nlon, self%dim_lon_id), err)
      call check_pio('def_dim lat',  pio_def_dim(self%pio_file, 'lat',  nlat, self%dim_lat_id), err)
      call check_pio('def_dim time', pio_def_dim(self%pio_file, 'time', PIO_UNLIMITED, self%dim_time_id), err)
      call check_pio('def_dim nbnd', pio_def_dim(self%pio_file, 'nbnd', 2,    self%dim_nbnd_id), err)

      call check_pio('def_var lon', pio_def_var(self%pio_file, 'lon', PIO_DOUBLE, &
           (/self%dim_lon_id/), self%pio_lon_var), err)
      ierr = pio_put_att(self%pio_file, self%pio_lon_var, 'units', 'degrees_east')
      ierr = pio_put_att(self%pio_file, self%pio_lon_var, 'long_name', 'longitude')
      ierr = pio_put_att(self%pio_file, self%pio_lon_var, 'axis', 'X')

      call check_pio('def_var lat', pio_def_var(self%pio_file, 'lat', PIO_DOUBLE, &
           (/self%dim_lat_id/), self%pio_lat_var), err)
      ierr = pio_put_att(self%pio_file, self%pio_lat_var, 'units', 'degrees_north')
      ierr = pio_put_att(self%pio_file, self%pio_lat_var, 'long_name', 'latitude')
      ierr = pio_put_att(self%pio_file, self%pio_lat_var, 'axis', 'Y')

      ! CF time coordinate. Numeric `time` is the END of the averaging
      ! window for averaged records and the sample time for instantaneous
      ! records (matches EAM convention). `time_bnds` carries the window
      ! endpoints unambiguously. `date` (YYYYMMDD) and `datesec` (seconds
      ! of day) mirror EAM's parallel integer pair for ncdump-grade
      ! human readability.
      call check_pio('def_var time', pio_def_var(self%pio_file, 'time', PIO_DOUBLE, &
           (/self%dim_time_id/), self%pio_time_var), err)
      ierr = pio_put_att(self%pio_file, self%pio_time_var, 'long_name', 'time')
      ierr = pio_put_att(self%pio_file, self%pio_time_var, 'units', &
           'days since ' // trim(time_epoch_cf))
      ierr = pio_put_att(self%pio_file, self%pio_time_var, 'calendar', 'noleap')
      ierr = pio_put_att(self%pio_file, self%pio_time_var, 'bounds', 'time_bnds')
      ierr = pio_put_att(self%pio_file, self%pio_time_var, 'axis', 'T')

      call check_pio('def_var time_bnds', pio_def_var(self%pio_file, 'time_bnds', PIO_DOUBLE, &
           (/self%dim_nbnd_id, self%dim_time_id/), self%pio_time_bnds_var), err)
      ierr = pio_put_att(self%pio_file, self%pio_time_bnds_var, 'long_name', &
           'time interval endpoints')

      call check_pio('def_var date', pio_def_var(self%pio_file, 'date', PIO_INT, &
           (/self%dim_time_id/), self%pio_date_var), err)
      ierr = pio_put_att(self%pio_file, self%pio_date_var, 'long_name', &
           'current date (YYYYMMDD)')

      call check_pio('def_var datesec', pio_def_var(self%pio_file, 'datesec', PIO_INT, &
           (/self%dim_time_id/), self%pio_datesec_var), err)
      ierr = pio_put_att(self%pio_file, self%pio_datesec_var, 'long_name', &
           'current seconds of current date')

      ierr = pio_put_att(self%pio_file, PIO_GLOBAL, 'source', 'E3SM FME online remapping')
      ierr = pio_put_att(self%pio_file, PIO_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = pio_put_att(self%pio_file, PIO_GLOBAL, 'calendar', 'noleap')
      ierr = pio_put_att(self%pio_file, PIO_GLOBAL, 'time_reference_date', &
           trim(time_epoch_cf))

      if (err /= 0) then
         call mpas_log_write('ocn_fme_remap_file: define-mode error opening ' &
              // trim(filename), messageType=MPAS_LOG_ERR)
         return
      end if

      self%n_pio_vars = 0
      self%time_record = 0
      self%file_is_open = .true.
      self%existing_file_reopened = .false.

   end subroutine remap_file_open!}}}

   subroutine check_pio(label, ierr, accumulated_err)
      ! Aggregate first non-zero PIO error.  Lets us issue every dim/var
      ! definition without aborting on the first failure but still
      ! propagate the failure to the caller.
      character(len=*), intent(in) :: label
      integer, intent(in) :: ierr
      integer, intent(inout) :: accumulated_err
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('ocn_fme_remap_file: PIO error in ' // trim(label) // ' (code $i)', &
              messageType=MPAS_LOG_ERR, intArgs=(/ierr/))
         if (accumulated_err == 0) accumulated_err = ierr
      end if
   end subroutine check_pio

   subroutine remap_file_def_var(self, varname, units, long_name, numlev, err, &
        is_static)!{{{
   !> Register an output variable. On a freshly created file, defines the
   !> variable + attributes. On a reopened (append-mode) file, looks up
   !> the existing var_desc — definitions and attributes are already
   !> persisted on disk from the first create.
   !>
   !> When is_static=.true., the var is defined with (lon, lat) only (no
   !> Time dim) and gets no _FillValue attribute. Static vars are written
   !> once per file via write_var, which skips pio_setframe for these.
      class(ocn_fme_remap_file_t), intent(inout) :: self
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(in) :: numlev
      integer, intent(out) :: err
      logical, optional, intent(in) :: is_static
      integer :: ierr
      real(kind=4) :: fill_val
      logical :: is_static_local

      err = 0
      if (.not. self%file_is_open) then; err = 1; return; end if

      is_static_local = .false.
      if (present(is_static)) is_static_local = is_static

      self%n_pio_vars = self%n_pio_vars + 1
      if (self%n_pio_vars > max_pio_vars) then
         call mpas_log_write('ocn_fme_remap_file: exceeded max_pio_vars=$i for var: ' &
              // trim(varname), messageType=MPAS_LOG_ERR, &
              intArgs=(/max_pio_vars/))
         self%n_pio_vars = self%n_pio_vars - 1
         err = 2; return
      end if

      self%pio_var_names(self%n_pio_vars) = varname
      self%pio_var_is_static(self%n_pio_vars) = is_static_local

      if (self%existing_file_reopened) then
         ierr = pio_inq_varid(self%pio_file, trim(varname), &
              self%pio_var_descs(self%n_pio_vars))
         if (ierr /= PIO_NOERR) then
            call mpas_log_write('ocn_fme_remap_file: pio_inq_varid failed for ' &
                 // trim(varname) // ' (code $i)', messageType=MPAS_LOG_ERR, &
                 intArgs=(/ierr/))
            err = ierr; return
         end if
         return
      end if

      if (is_static_local) then
         ierr = pio_def_var(self%pio_file, trim(varname), PIO_REAL, &
              (/self%dim_lon_id, self%dim_lat_id/), &
              self%pio_var_descs(self%n_pio_vars))
      else
         ierr = pio_def_var(self%pio_file, trim(varname), PIO_REAL, &
              (/self%dim_lon_id, self%dim_lat_id, self%dim_time_id/), &
              self%pio_var_descs(self%n_pio_vars))
      end if
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('ocn_fme_remap_file: pio_def_var failed for ' &
              // trim(varname) // ' (code $i)', messageType=MPAS_LOG_ERR, &
              intArgs=(/ierr/))
         err = ierr; return
      end if
      ierr = pio_put_att(self%pio_file, self%pio_var_descs(self%n_pio_vars), &
           'units', trim(units))
      ierr = pio_put_att(self%pio_file, self%pio_var_descs(self%n_pio_vars), &
           'long_name', trim(long_name))
      ! Fill-value attributes only on time-varying vars; masks contain
      ! literal 0.0/1.0 with no missing data semantics.
      if (.not. is_static_local) then
         fill_val = real(SHR_FILL_VALUE, 4)
         ierr = pio_put_att(self%pio_file, self%pio_var_descs(self%n_pio_vars), &
              '_FillValue', fill_val)
         ierr = pio_put_att(self%pio_file, self%pio_var_descs(self%n_pio_vars), &
              'missing_value', fill_val)
      end if

   end subroutine remap_file_def_var!}}}

   subroutine remap_file_write_var(self, varname, fld_out, numlev, err)!{{{
      class(ocn_fme_remap_file_t), intent(inout) :: self
      character(len=*), intent(in) :: varname
      real(r8), intent(in) :: fld_out(:,:)
      integer, intent(in) :: numlev
      integer, intent(out) :: err
      integer :: iv, ierr
      ! Must match the PIO_REAL type used in pio_initdecomp for iodesc_2d.
      ! Passing real(r8) to a PIO_REAL decomposition garbles the output
      ! (PIO reads 4-byte chunks from 8-byte elements).
      real(4) :: fld_1d_r4(remap_data%n_b_local)

      err = 0
      if (.not. self%file_is_open) then; err = 1; return; end if

      do iv = 1, self%n_pio_vars
         if (trim(self%pio_var_names(iv)) == trim(varname)) then
            ! Static (no Time dim) vars: skip setframe — there's no
            ! UNLIMITED dim to advance, and pio_setframe on a 2D var
            ! errors on some PIO versions.
            if (.not. self%pio_var_is_static(iv)) then
               call pio_setframe(self%pio_file, self%pio_var_descs(iv), &
                    int(self%time_record, PIO_OFFSET_KIND))
            end if
            fld_1d_r4(1:remap_data%n_b_local) = &
                 real(fld_out(1:remap_data%n_b_local, 1), 4)
            call pio_write_darray(self%pio_file, self%pio_var_descs(iv), &
                 iodesc_2d, fld_1d_r4, ierr)
            err = ierr
            return
         end if
      end do
      err = 3

   end subroutine remap_file_write_var!}}}

   subroutine remap_file_write_time(self, time_lo, time_hi, err)!{{{
   !> Write the per-record CF time variables (time, time_bnds, date,
   !> datesec) for the next record and advance time_record.
   !>
   !> time_lo / time_hi are the start and end of the averaging window
   !> for averaged records (set time_lo == time_hi for instantaneous).
   !> Numeric `time` is stamped at time_hi (end of window) to match
   !> EAM convention; time_bnds carries [lo, hi] explicitly.
      class(ocn_fme_remap_file_t), intent(inout) :: self
      type(MPAS_Time_Type), intent(in) :: time_lo, time_hi
      integer, intent(out) :: err
      integer :: ierr, yr, mo, dy, hr, mn, sc
      real(r8) :: t_lo_days, t_hi_days
      real(r8) :: time_buf(1)
      real(r8) :: bnds_buf(2, 1)
      integer  :: date_buf(1), datesec_buf(1)

      err = 0
      if (.not. self%file_is_open) then; err = 1; return; end if

      ! Fresh-create files: leave define-mode and seed lon/lat once.
      ! Reopened (append-mode) files are already in data mode and lon/lat
      ! are already populated, so this whole block is skipped.
      if (self%time_record == 0 .and. .not. self%existing_file_reopened) then
         ierr = pio_enddef(self%pio_file)
         ierr = pio_put_var(self%pio_file, self%pio_lon_var, &
              remap_data%lon(1:remap_data%nlon))
         ierr = pio_put_var(self%pio_file, self%pio_lat_var, &
              remap_data%lat(1:remap_data%nlat))
      end if
      self%time_record = self%time_record + 1

      call time_to_days_since_epoch(time_lo, t_lo_days, ierr)
      call time_to_days_since_epoch(time_hi, t_hi_days, ierr)

      time_buf(1) = t_hi_days
      bnds_buf(1, 1) = t_lo_days
      bnds_buf(2, 1) = t_hi_days

      call mpas_get_time(time_hi, YYYY=yr, MM=mo, DD=dy, &
           H=hr, M=mn, S=sc, ierr=ierr)
      date_buf(1)    = yr * 10000 + mo * 100 + dy
      datesec_buf(1) = hr * 3600  + mn * 60  + sc

      ierr = pio_put_var(self%pio_file, self%pio_time_var, &
           (/self%time_record/), (/1/), time_buf)
      ierr = pio_put_var(self%pio_file, self%pio_time_bnds_var, &
           (/1, self%time_record/), (/2, 1/), bnds_buf)
      ierr = pio_put_var(self%pio_file, self%pio_date_var, &
           (/self%time_record/), (/1/), date_buf)
      ierr = pio_put_var(self%pio_file, self%pio_datesec_var, &
           (/self%time_record/), (/1/), datesec_buf)

   end subroutine remap_file_write_time!}}}

   subroutine remap_file_close(self, err)!{{{
      class(ocn_fme_remap_file_t), intent(inout) :: self
      integer, intent(out) :: err
      err = 0
      if (self%file_is_open) then
         call pio_closefile(self%pio_file)
         self%file_is_open = .false.
         self%n_pio_vars = 0
         self%time_record = 0
         self%current_filename = ''
         self%existing_file_reopened = .false.
      end if
   end subroutine remap_file_close!}}}

   subroutine remap_file_register_var(self, varname, units, long_name, err, &
        is_static)!{{{
   !> Register a variable for output. Stores metadata so that def_var can
   !> be replayed automatically when a file is opened or rotated.
   !> Silently skips if the variable is already registered (MPAS framework
   !> may call init twice: once for setup, once for compute_on_startup).
   !> When is_static=.true., the variable is defined with shape (lon, lat)
   !> only (no Time dim) and written once per file (used for masks).
      class(ocn_fme_remap_file_t), intent(inout) :: self
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(out) :: err
      logical, optional, intent(in) :: is_static
      integer :: iv

      err = 0

      ! Skip if already registered (guard against double-init)
      do iv = 1, self%n_registered_vars
         if (trim(self%pio_var_names(iv)) == trim(varname)) return
      end do

      self%n_registered_vars = self%n_registered_vars + 1
      if (self%n_registered_vars > max_pio_vars) then
         call mpas_log_write('ocn_fme_remap_file: exceeded max_pio_vars=$i for var: ' &
              // trim(varname), messageType=MPAS_LOG_ERR, &
              intArgs=(/max_pio_vars/))
         self%n_registered_vars = self%n_registered_vars - 1
         err = 2; return
      end if
      self%pio_var_names(self%n_registered_vars) = varname
      self%stored_units(self%n_registered_vars) = units
      self%stored_long_names(self%n_registered_vars) = long_name
      self%stored_is_static(self%n_registered_vars) = .false.
      if (present(is_static)) then
         self%stored_is_static(self%n_registered_vars) = is_static
      end if
   end subroutine remap_file_register_var!}}}

   subroutine remap_file_check_rotate(self, domain, err)!{{{
   !> Check whether the file needs rotation based on the current model date.
   !> Builds a filename from the stream template + current date. If different
   !> from the currently open file, closes the old file and opens a new one,
   !> replaying all registered variable definitions.
      class(ocn_fme_remap_file_t), intent(inout) :: self
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=512) :: new_filename
      integer :: ierr, iv

      err = 0
      if (self%n_registered_vars == 0) return
      if (len_trim(self%stream_name) == 0) return

      call ocn_fme_remap_build_filename(domain, trim(self%stream_name), &
           new_filename, ierr)
      if (ierr /= 0) then; err = ierr; return; end if

      ! First call or date boundary crossed: open (or rotate) file
      if (.not. self%file_is_open .or. &
           trim(new_filename) /= trim(self%current_filename)) then

         if (self%file_is_open) then
            call self%close_file(ierr)
         end if

         call self%open_file(trim(new_filename), domain, ierr)
         if (ierr /= 0) then; err = ierr; return; end if
         self%current_filename = new_filename

         ! Replay all registered variable definitions, preserving
         ! static-var shape across rotation.
         do iv = 1, self%n_registered_vars
            call self%def_var(trim(self%pio_var_names(iv)), &
                 trim(self%stored_units(iv)), &
                 trim(self%stored_long_names(iv)), 1, ierr, &
                 is_static=self%stored_is_static(iv))
         end do
      end if

   end subroutine remap_file_check_rotate!}}}

!***********************************************************************

   subroutine ocn_fme_horiz_remap_get_dims(nlat, nlon, n_b_local)!{{{
      integer, intent(out) :: nlat, nlon, n_b_local
      nlat = remap_data%nlat
      nlon = remap_data%nlon
      n_b_local = remap_data%n_b_local
   end subroutine ocn_fme_horiz_remap_get_dims!}}}

!***********************************************************************

   logical function ocn_fme_horiz_remap_active()!{{{
      ocn_fme_horiz_remap_active = is_initialized
   end function ocn_fme_horiz_remap_active!}}}

!***********************************************************************

   subroutine ocn_fme_remap_build_filename(domain, streamName, filename, err)!{{{
   !> Build a date-stamped remapped output filename by querying the stream
   !> manager for the native stream's filename_template and substituting
   !> the current date. The native stream's filename_template is set by
   !> buildnml and contains the casename, instance string, and $Y-$M-$D
   !> date placeholders (e.g., "mycase.mpaso.hist.am.fmeDerived.$Y-$M-$D.nc").
   !> We replace the ".nc" extension with ".remapped.nc" and substitute dates.

      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: streamName
      character(len=*), intent(out) :: filename
      integer, intent(out) :: err

      character(len=StrKIND) :: tmpl
      type(MPAS_Time_Type) :: currTime
      integer :: yr, mo, dy, ierr, idx
      character(len=4) :: yr_str
      character(len=2) :: mo_str, dy_str

      err = 0
      filename = ''

      ! Query the stream manager for the native stream's filename template
      call MPAS_stream_mgr_get_property(domain % streamManager, &
           trim(streamName), MPAS_STREAM_PROPERTY_FILENAME, tmpl, ierr=ierr)
      if (ierr /= 0) then
         call mpas_log_write('ocn_fme_remap_build_filename: could not query stream ' &
              // trim(streamName) // ', using fallback name', messageType=MPAS_LOG_WARN)
         ! Fallback: use stream name directly with date
         tmpl = trim(streamName) // '.$Y-$M-$D.nc'
      end if

      ! Insert ".remapped" before ".nc" extension
      idx = index(tmpl, '.nc', back=.true.)
      if (idx > 0) then
         tmpl = tmpl(1:idx-1) // '.remapped.nc'
      else
         tmpl = trim(tmpl) // '.remapped.nc'
      end if

      ! Substitute date placeholders from current model time
      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(currTime, YYYY=yr, MM=mo, DD=dy, ierr=ierr)

      write(yr_str, '(I4.4)') yr
      write(mo_str, '(I2.2)') mo
      write(dy_str, '(I2.2)') dy

      ! Replace $Y, $M, $D tokens
      filename = tmpl
      call replace_token(filename, '$Y', yr_str)
      call replace_token(filename, '$M', mo_str)
      call replace_token(filename, '$D', dy_str)

   end subroutine ocn_fme_remap_build_filename!}}}

!***********************************************************************

   subroutine replace_token(str, token, replacement)!{{{
      character(len=*), intent(inout) :: str
      character(len=*), intent(in) :: token, replacement
      integer :: idx
      character(len=len(str)) :: tmp

      idx = index(str, token)
      if (idx > 0) then
         tmp = str(1:idx-1) // trim(replacement) // str(idx+len(token):)
         str = tmp
      end if
   end subroutine replace_token!}}}

!***********************************************************************

   subroutine fme_accum_filename(amName, filename)!{{{
   !> Sidecar restart-state file for an FME analysis member's remap-space
   !> accumulators. Fixed name in the run dir; one file per AM, overwritten
   !> on every restart-write event. No date stamp because there is only
   !> ever one valid snapshot at a time (the most recent restart point).
      character(len=*), intent(in)  :: amName
      character(len=*), intent(out) :: filename
      filename = trim(amName) // '.fme_accum_restart.nc'
   end subroutine fme_accum_filename!}}}

!***********************************************************************

   subroutine ocn_fme_horiz_remap_write_accum(amName, accum_arr, &
        n_avg_fields, nAccum, avg_last_output_time_str, latest_inst_arr, &
        count_arr, err)!{{{
   !> Persist the remap-space accumulator state for one AM to its sidecar
   !> file. Called from the MPAS-O driver at the moment the coupler has
   !> requested a restart write (see `seq_timemgr_RestartAlarmIsOn` branch
   !> in ocn_comp_mct.F). Each call overwrites the previous sidecar.
   !>
   !> accum_arr has shape (n_b_local, n_avg_fields). latest_inst_arr is
   !> optional and matches accum_arr; pass it for AMs that maintain
   !> instantaneous companions (depth-coarsening) and omit otherwise.
   !> count_arr is the per-cell valid-sample counter parallel to
   !> accum_arr; passing it preserves the time-averaging mask across
   !> restart so warm-leg output continues to weight intensive
   !> fill-prone fields correctly (see AGENTS.md gotcha on
   !> time-averaging fill pollution).
      character(len=*), intent(in) :: amName
      integer,          intent(in) :: n_avg_fields
      real(r8),         intent(in) :: accum_arr(:,:)
      integer,          intent(in) :: nAccum
      character(len=*), intent(in) :: avg_last_output_time_str
      real(r8), optional, intent(in) :: latest_inst_arr(:,:)
      real(r8), optional, intent(in) :: count_arr(:,:)
      integer,          intent(out) :: err

      type(file_desc_t)        :: pio_acc_file
      type(var_desc_t)         :: var_accum, var_inst, var_count
      integer                  :: dim_lon, dim_lat, dim_field
      integer                  :: ierr, idx
      character(len=512)       :: filename
      real(r8), allocatable    :: slab(:)

      err = 0
      if (.not. is_initialized) return

      call fme_accum_filename(trim(amName), filename)

      ierr = pio_createfile(pio_sys, pio_acc_file, PIO_IOTYPE_PNETCDF, &
           trim(filename), PIO_CLOBBER)
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('ocn_fme_write_accum: failed to create ' &
              // trim(filename), messageType=MPAS_LOG_WARN)
         err = ierr; return
      end if

      call check_pio('accum def_dim lon', &
           pio_def_dim(pio_acc_file, 'lon', remap_data%nlon, dim_lon), err)
      call check_pio('accum def_dim lat', &
           pio_def_dim(pio_acc_file, 'lat', remap_data%nlat, dim_lat), err)
      ! n_avg_fields is the record dimension — must be UNLIMITED so that
      ! pio_setframe + pio_write_darray walks across slots.
      call check_pio('accum def_dim n_avg_fields', &
           pio_def_dim(pio_acc_file, 'n_avg_fields', PIO_UNLIMITED, dim_field), err)

      call check_pio('accum def_var remap_accum', &
           pio_def_var(pio_acc_file, 'remap_accum', PIO_DOUBLE, &
                       (/dim_lon, dim_lat, dim_field/), var_accum), err)
      if (present(latest_inst_arr)) then
         call check_pio('accum def_var remap_latest_inst', &
              pio_def_var(pio_acc_file, 'remap_latest_inst', PIO_DOUBLE, &
                          (/dim_lon, dim_lat, dim_field/), var_inst), err)
      end if
      if (present(count_arr)) then
         call check_pio('accum def_var remap_accum_count', &
              pio_def_var(pio_acc_file, 'remap_accum_count', PIO_DOUBLE, &
                          (/dim_lon, dim_lat, dim_field/), var_count), err)
      end if

      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, 'nAccum', nAccum)
      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, 'n_avg_fields', n_avg_fields)
      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, &
           'avg_last_output_time', trim(avg_last_output_time_str))
      ! Persist the CF time epoch so init_time_epoch can recover it on
      ! warm restart -- MPAS_START_TIME from the clock would otherwise
      ! return the leg-start time on restart, not the case start.
      if (time_epoch_set) then
         ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, &
              'time_reference_date', trim(time_epoch_str))
      end if
      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, 'source', &
           'E3SM FME accumulator restart sidecar')

      ierr = pio_enddef(pio_acc_file)
      if (err /= 0) then
         call pio_closefile(pio_acc_file); return
      end if

      allocate(slab(remap_data%n_b_local))
      do idx = 1, n_avg_fields
         call pio_setframe(pio_acc_file, var_accum, &
              int(idx, PIO_OFFSET_KIND))
         slab(:) = accum_arr(:, idx)
         call pio_write_darray(pio_acc_file, var_accum, iodesc_2d_dbl, &
              slab, ierr)
         if (ierr /= PIO_NOERR) err = ierr
         if (present(latest_inst_arr)) then
            call pio_setframe(pio_acc_file, var_inst, &
                 int(idx, PIO_OFFSET_KIND))
            slab(:) = latest_inst_arr(:, idx)
            call pio_write_darray(pio_acc_file, var_inst, iodesc_2d_dbl, &
                 slab, ierr)
            if (ierr /= PIO_NOERR) err = ierr
         end if
         if (present(count_arr)) then
            call pio_setframe(pio_acc_file, var_count, &
                 int(idx, PIO_OFFSET_KIND))
            slab(:) = count_arr(:, idx)
            call pio_write_darray(pio_acc_file, var_count, iodesc_2d_dbl, &
                 slab, ierr)
            if (ierr /= PIO_NOERR) err = ierr
         end if
      end do
      deallocate(slab)

      call pio_closefile(pio_acc_file)

   end subroutine ocn_fme_horiz_remap_write_accum!}}}

!***********************************************************************

   subroutine ocn_fme_horiz_remap_read_accum(amName, accum_arr, &
        n_avg_fields, nAccum, avg_last_output_time_str, latest_inst_arr, &
        count_arr, success)!{{{
   !> Restore an AM's remap-space accumulator state from its sidecar file.
   !> Called from the AM's *_init_* when config_do_restart=.true.
   !>
   !> Returns success=.true. only if the file exists, has matching
   !> n_avg_fields, and all reads succeed. On any failure (no file,
   !> dim mismatch, PIO error) returns success=.false. and leaves the
   !> caller's arrays untouched — the AM should treat this as a normal
   !> fresh-init (zero accumulators, MPAS_NOW baseline).
      character(len=*),   intent(in)    :: amName
      integer,            intent(in)    :: n_avg_fields
      real(r8),           intent(inout) :: accum_arr(:,:)
      integer,            intent(out)   :: nAccum
      character(len=*),   intent(out)   :: avg_last_output_time_str
      real(r8), optional, intent(inout) :: latest_inst_arr(:,:)
      real(r8), optional, intent(inout) :: count_arr(:,:)
      logical,            intent(out)   :: success

      type(file_desc_t)     :: pio_acc_file
      type(var_desc_t)      :: var_accum, var_inst, var_count
      integer               :: ierr, idx, n_avg_on_disk
      character(len=512)    :: filename
      logical               :: file_exists
      real(r8), allocatable :: slab(:)
      real(r8), allocatable :: tmp_accum(:,:), tmp_inst(:,:), tmp_count(:,:)

      success = .false.
      nAccum = 0
      avg_last_output_time_str = ''
      if (.not. is_initialized) return

      call fme_accum_filename(trim(amName), filename)
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) return

      ierr = pio_openfile(pio_sys, pio_acc_file, PIO_IOTYPE_PNETCDF, &
           trim(filename), PIO_WRITE)
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('ocn_fme_read_accum: pio_openfile failed for ' &
              // trim(filename), messageType=MPAS_LOG_WARN)
         return
      end if

      ierr = pio_get_att(pio_acc_file, PIO_GLOBAL, 'n_avg_fields', n_avg_on_disk)
      if (ierr /= PIO_NOERR .or. n_avg_on_disk /= n_avg_fields) then
         call mpas_log_write('ocn_fme_read_accum: n_avg_fields mismatch ' &
              // 'in ' // trim(filename) // ' (on-disk=$i, expected=$i)', &
              messageType=MPAS_LOG_WARN, &
              intArgs=(/n_avg_on_disk, n_avg_fields/))
         call pio_closefile(pio_acc_file); return
      end if

      ierr = pio_get_att(pio_acc_file, PIO_GLOBAL, 'nAccum', nAccum)
      if (ierr /= PIO_NOERR) nAccum = 0
      ierr = pio_get_att(pio_acc_file, PIO_GLOBAL, 'avg_last_output_time', &
           avg_last_output_time_str)
      if (ierr /= PIO_NOERR) avg_last_output_time_str = ''

      ierr = pio_inq_varid(pio_acc_file, 'remap_accum', var_accum)
      if (ierr /= PIO_NOERR) then
         call pio_closefile(pio_acc_file); return
      end if

      allocate(tmp_accum(remap_data%n_b_local, n_avg_fields))
      allocate(slab(remap_data%n_b_local))
      do idx = 1, n_avg_fields
         call pio_setframe(pio_acc_file, var_accum, &
              int(idx, PIO_OFFSET_KIND))
         call pio_read_darray(pio_acc_file, var_accum, iodesc_2d_dbl, &
              slab, ierr)
         if (ierr /= PIO_NOERR) then
            deallocate(tmp_accum, slab)
            call pio_closefile(pio_acc_file); return
         end if
         tmp_accum(:, idx) = slab(:)
      end do

      if (present(latest_inst_arr)) then
         ierr = pio_inq_varid(pio_acc_file, 'remap_latest_inst', var_inst)
         if (ierr == PIO_NOERR) then
            allocate(tmp_inst(remap_data%n_b_local, n_avg_fields))
            do idx = 1, n_avg_fields
               call pio_setframe(pio_acc_file, var_inst, &
                    int(idx, PIO_OFFSET_KIND))
               call pio_read_darray(pio_acc_file, var_inst, iodesc_2d_dbl, &
                    slab, ierr)
               if (ierr /= PIO_NOERR) then
                  deallocate(tmp_accum, tmp_inst, slab)
                  call pio_closefile(pio_acc_file); return
               end if
               tmp_inst(:, idx) = slab(:)
            end do
            latest_inst_arr(:,:) = tmp_inst(:,:)
            deallocate(tmp_inst)
         end if
      end if

      if (present(count_arr)) then
         ierr = pio_inq_varid(pio_acc_file, 'remap_accum_count', var_count)
         if (ierr == PIO_NOERR) then
            allocate(tmp_count(remap_data%n_b_local, n_avg_fields))
            do idx = 1, n_avg_fields
               call pio_setframe(pio_acc_file, var_count, &
                    int(idx, PIO_OFFSET_KIND))
               call pio_read_darray(pio_acc_file, var_count, iodesc_2d_dbl, &
                    slab, ierr)
               if (ierr /= PIO_NOERR) then
                  deallocate(tmp_accum, tmp_count, slab)
                  call pio_closefile(pio_acc_file); return
               end if
               tmp_count(:, idx) = slab(:)
            end do
            count_arr(:,:) = tmp_count(:,:)
            deallocate(tmp_count)
         end if
      end if

      accum_arr(:,:) = tmp_accum(:,:)
      deallocate(tmp_accum, slab)
      call pio_closefile(pio_acc_file)

      success = .true.
      call mpas_log_write('ocn_fme_read_accum: restored accumulators from ' &
           // trim(filename) // ' (nAccum=$i)', intArgs=(/nAccum/))

   end subroutine ocn_fme_horiz_remap_read_accum!}}}

end module ocn_fme_horiz_remap

! vim: foldmethod=marker
