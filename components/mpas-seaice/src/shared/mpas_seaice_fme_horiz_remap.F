! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_fme_horiz_remap
!
!> \brief MPAS sea ice FME horizontal remapping wrapper
!> \author FME Training Data Support
!> \date   April 2026
!> \details
!>  Wrapper around shr_horiz_remap_t for MPAS-Sea-Ice with PIO output.
!>  Same architecture as ocn_fme_horiz_remap.
!>
!-----------------------------------------------------------------------

module seaice_fme_horiz_remap

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
                  pio_openfile, pio_createfile, pio_closefile, &
                  pio_def_dim, pio_def_var, pio_enddef, &
                  pio_put_att, pio_put_var, pio_write_darray, pio_read_darray, &
                  pio_get_att, pio_get_var, &
                  pio_inq_dimid, pio_inq_varid, pio_inq_dimlen, &
                  pio_initdecomp, pio_setframe, &
                  PIO_IOTYPE_PNETCDF, PIO_CLOBBER, PIO_WRITE, &
                  PIO_REAL, PIO_DOUBLE, PIO_INT, &
                  PIO_OFFSET_KIND, PIO_GLOBAL, PIO_NOERR, &
                  PIO_UNLIMITED

   implicit none
   private
   save

   public :: seaice_fme_horiz_remap_init
   public :: seaice_fme_horiz_remap_init_time_epoch
   public :: seaice_fme_horiz_remap_set_provenance
   public :: seaice_fme_horiz_remap_field
   public :: seaice_fme_horiz_remap_active
   public :: seaice_fme_horiz_remap_open_file
   public :: seaice_fme_horiz_remap_def_var
   public :: seaice_fme_horiz_remap_write_var
   public :: seaice_fme_horiz_remap_write_time
   public :: seaice_fme_horiz_remap_close_file
   public :: seaice_fme_horiz_remap_get_dims
   public :: seaice_fme_remap_build_filename
   public :: seaice_fme_horiz_remap_register_var
   public :: seaice_fme_horiz_remap_check_rotate
   public :: seaice_fme_horiz_remap_write_accum
   public :: seaice_fme_horiz_remap_read_accum
   public :: remap_stream_name, file_is_open
   ! 5D-companion analogues -- parallel state-and-routine set so the
   ! AM can drive a second output file (longer cadence) alongside the
   ! primary stream without touching the primary path.
   public :: seaice_fme_horiz_remap_register_var_5D
   public :: seaice_fme_horiz_remap_check_rotate_5D
   public :: seaice_fme_horiz_remap_write_var_5D
   public :: seaice_fme_horiz_remap_write_time_5D
   public :: seaice_fme_horiz_remap_close_file_5D
   public :: seaice_fme_horiz_remap_write_accum_5D
   public :: seaice_fme_horiz_remap_read_accum_5D
   public :: remap_stream_name_5D, file_is_open_5D
   public :: seaice_fme_horiz_remap_set_rewind

   type(shr_horiz_remap_t) :: remap_data
   integer, allocatable :: send_local_cell(:)
   real(r8), allocatable :: ws_send_buf(:)
   logical :: is_initialized = .false.

   ! Cached at init from config_do_restart. Append-mode only fires on
   ! warm restart; cold start always clobber-creates, defending against
   ! stale `.remapped.nc` left behind by a prior CIME compare-copy.
   logical :: is_restart_run = .false.

   ! Set by the seaice AM's *_read_restart when it detects a branch/rewind
   ! restart (persisted avg_last_output_time later than MPAS_NOW). When
   ! true, the file-open paths clobber-create the re-crossed monthly files
   ! instead of append-reopening, so re-simulated months get fresh records.
   ! Stays true for the rest of the run; harmless once the clock catches up.
   logical :: rewind_restart = .false.

   ! MPAS-SI fill value for input field detection (matches def_var _FillValue)
   real(r8), parameter :: MPAS_SI_FILL = SHR_FILL_VALUE

   ! CF time-coordinate epoch shared across all FME output files for
   ! this component. Auto-discovered at AM init via init_time_epoch
   ! (sidecar lookup -> MPAS_START_TIME fallback). See the parallel
   ! block in the ocean horiz_remap module for the full design rationale.
   type(MPAS_Time_Type) :: time_epoch
   character(len=64)    :: time_epoch_str = ''
   character(len=64)    :: time_epoch_cf  = ''
   logical              :: time_epoch_set = .false.

   ! PIO output
   type(file_desc_t) :: pio_file
   type(iosystem_desc_t), pointer :: pio_sys => null()
   logical :: file_is_open = .false.
   integer :: dim_lon_id, dim_lat_id, dim_time_id
   integer :: dim_nbnd_id = 0
   integer :: time_record = 0
   type(io_desc_t) :: iodesc_2d        ! PIO_REAL  (output records)
   type(io_desc_t) :: iodesc_2d_dbl    ! PIO_DOUBLE (accumulator restart)
   logical :: iodesc_2d_valid = .false.
   ! True when the open file was reopened from disk (append after restart)
   ! rather than freshly created. Mirrors the ocean module's flag.
   logical :: existing_file_reopened = .false.
   ! Fire-once latch limiting restart append mode to the leg-start month's
   ! file. is_restart_run stays true for the whole run, so without this
   ! latch every NEW month rotated to during the leg would also try to
   ! reopen-append; reopening a stale/0-byte file trips PIO's "invalid IO
   ! type" assert. Mirrors the ocean module fix.
   logical :: restart_reopen_pending = .true.

   integer, parameter :: max_pio_vars = 512
   integer :: n_pio_vars = 0
   character(len=64) :: pio_var_names(max_pio_vars)
   type(var_desc_t) :: pio_var_descs(max_pio_vars)
   type(var_desc_t) :: pio_lon_var, pio_lat_var
   ! CF time variables (replace legacy MPAS xtime).
   type(var_desc_t) :: pio_time_var, pio_time_bnds_var
   type(var_desc_t) :: pio_date_var, pio_datesec_var

   ! Stored var metadata for replay on rotation
   integer :: n_registered_vars = 0
   character(len=64)  :: stored_units(max_pio_vars)
   character(len=128) :: stored_long_names(max_pio_vars)
   ! Optional CF cell_methods (e.g. 'time: mean' / 'time: point').
   ! Empty string => attribute is not written.
   character(len=64)  :: stored_cell_methods(max_pio_vars) = ''
   character(len=512) :: current_filename = ''
   character(len=128) :: remap_stream_name = ''

   ! CIME-side provenance, captured once via seaice_fme_horiz_remap_set_provenance
   ! and stamped as global attrs onto every fresh `.remapped.nc` file.
   character(len=StrKIND) :: prov_case          = ''
   character(len=StrKIND) :: prov_username      = ''
   character(len=StrKIND) :: prov_hostname      = ''
   character(len=StrKIND) :: prov_model_version = ''
   character(len=StrKIND) :: prov_history       = ''
   logical                :: provenance_set     = .false.

   ! ---- 5D companion state (parallel to the primary stream above) ----
   ! Independent PIO file context so the AM can drive two output files at
   ! different cadences. iodesc_2d/_dbl, remap_data, time_epoch_*, and
   ! provenance are all shared with the primary path (single mesh, single
   ! time epoch, single provenance — only the file handle and per-file
   ! var IDs differ).
   type(file_desc_t)  :: pio_file_5D
   logical            :: file_is_open_5D = .false.
   integer            :: dim_lon_id_5D, dim_lat_id_5D, dim_time_id_5D
   integer            :: dim_nbnd_id_5D = 0
   integer            :: time_record_5D = 0
   logical            :: existing_file_reopened_5D = .false.
   ! Fire-once restart-append latch for the 5D companion (see primary).
   logical            :: restart_reopen_pending_5D = .true.
   integer            :: n_pio_vars_5D = 0
   character(len=64)  :: pio_var_names_5D(max_pio_vars)
   type(var_desc_t)   :: pio_var_descs_5D(max_pio_vars)
   type(var_desc_t)   :: pio_lon_var_5D, pio_lat_var_5D
   type(var_desc_t)   :: pio_time_var_5D, pio_time_bnds_var_5D
   type(var_desc_t)   :: pio_date_var_5D, pio_datesec_var_5D
   integer            :: n_registered_vars_5D = 0
   character(len=64)  :: stored_units_5D(max_pio_vars)
   character(len=128) :: stored_long_names_5D(max_pio_vars)
   character(len=64)  :: stored_cell_methods_5D(max_pio_vars) = ''
   character(len=512) :: current_filename_5D = ''
   character(len=128) :: remap_stream_name_5D = ''

contains

   subroutine seaice_fme_horiz_remap_init(domain, mapfile, err)!{{{
      use mpi, only: MPI_INTEGER, MPI_MAX

      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: mapfile
      integer, intent(out) :: err

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: meshPool
      integer, pointer :: nCellsSolve
      integer, dimension(:), pointer :: indexToCellID
      integer :: myrank, nprocs, comm, ierr, iCell, gcol, i, map_col
      integer, allocatable :: gcol_to_rank(:), gcol_to_rank_global(:), send_gcol_list(:)
      integer(PIO_OFFSET_KIND), allocatable :: idof(:)
      integer :: N_total, ocean_count
      integer, allocatable :: is_ocean(:), gcol_to_ocean_col(:)
      integer, dimension(:), pointer :: cellMask

      err = 0
      if (is_initialized) return
      if (len_trim(mapfile) == 0 .or. trim(mapfile) == 'none') then
         is_initialized = .false.; return
      end if

      myrank = domain % dminfo % my_proc_id
      nprocs = domain % dminfo % nprocs
      comm   = domain % dminfo % comm

      call remap_data%read_mapfile(trim(mapfile), comm, myrank, nprocs, &
           domain % ioContext % pio_iosystem, ierr)
      if (ierr /= 0) then
         err = ierr; return
      end if

      ! The FME remap always uses the masked-renormalization path, which is
      ! only safe for non-negative maps. A negative-weight map (TempestRemap
      ! 2nd-order FV / "trfv2") silently corrupts partially-covered cells.
      ! Refuse it loudly at init rather than shipping a corrupt training tape.
      ! Use a conservative or bilinear map instead. See AGENTS.md.
      if (remap_data%has_negative_weights) then
         call mpas_log_write('seaice_fme_horiz_remap_init: map file ' // trim(mapfile) &
              // ' has NEGATIVE interpolation weights (min weight = $r). This is ' &
              // 'UNSAFE for the masked FME remap path and corrupts ' &
              // 'partially-covered cells (out-of-range values). Use a ' &
              // 'non-negative map (first-order conservative or bilinear).', &
              messageType=MPAS_LOG_CRIT, realArgs=(/ remap_data%min_weight /))
         err = 1
         return
      end if

      ! Resolve cell ID mismatch (same issue as MPAS-O): map files from
      ! masked SCRIP meshes use contiguous 1-based ocean-only indexing.

      ! Find global N_total
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
         do i = 1, N_total
            gcol_to_ocean_col(i) = i
         end do
      else
         ! Build ocean mask via cellMask (sea ice equivalent of maxLevelCell)
         allocate(is_ocean(N_total))
         is_ocean(:) = 0
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
            call mpas_pool_get_array(meshPool, 'cellMask', cellMask)
            do iCell = 1, nCellsSolve
               if (cellMask(iCell) > 0) then
                  is_ocean(indexToCellID(iCell)) = 1
               end if
            end do
            block => block % next
         end do

         allocate(gcol_to_rank_global(N_total))
         call mpi_allreduce(is_ocean, gcol_to_rank_global, N_total, &
              MPI_INTEGER, MPI_MAX, comm, ierr)
         is_ocean = gcol_to_rank_global
         deallocate(gcol_to_rank_global)

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
            call mpas_log_write('seaice_fme_horiz_remap_init: N_total=$i, n_a=$i, ocean_count=$i', &
                 intArgs=(/N_total, remap_data%n_a, ocean_count/))
         end if
      end if

      allocate(gcol_to_rank(remap_data%n_a))
      gcol_to_rank(:) = -1
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
         do iCell = 1, nCellsSolve
            map_col = gcol_to_ocean_col(indexToCellID(iCell))
            if (map_col >= 1 .and. map_col <= remap_data%n_a) gcol_to_rank(map_col) = myrank
         end do
         block => block % next
      end do

      ! Gather full gcol_to_rank across all ranks
      allocate(gcol_to_rank_global(remap_data%n_a))
      call mpi_allreduce(gcol_to_rank, gcol_to_rank_global, remap_data%n_a, &
           MPI_INTEGER, MPI_MAX, comm, ierr)
      gcol_to_rank = gcol_to_rank_global
      deallocate(gcol_to_rank_global)

      call remap_data%build_comm(gcol_to_rank, comm, myrank, nprocs, &
           send_gcol_list, ierr)
      deallocate(gcol_to_rank)
      if (ierr /= 0) then
         err = ierr; return
      end if

      allocate(send_local_cell(remap_data%n_send_total))
      send_local_cell(:) = 0
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
         do i = 1, remap_data%n_send_total
            if (send_local_cell(i) > 0) cycle
            do iCell = 1, nCellsSolve
               map_col = gcol_to_ocean_col(indexToCellID(iCell))
               if (map_col == send_gcol_list(i)) then
                  send_local_cell(i) = iCell; exit
               end if
            end do
         end do
         block => block % next
      end do
      if (allocated(send_gcol_list)) deallocate(send_gcol_list)
      if (allocated(gcol_to_ocean_col)) deallocate(gcol_to_ocean_col)

      ! Validate that all needed source cells were found locally.
      do i = 1, remap_data%n_send_total
         if (send_local_cell(i) == 0) then
            call mpas_log_write('seaice_fme_horiz_remap_init: FATAL - source cell $i ' // &
                 'not found in local blocks. Map file may be incompatible with mesh.', &
                 messageType=MPAS_LOG_CRIT, intArgs=(/i/))
            err = 1
            return
         end if
      end do

      pio_sys => domain % ioContext % pio_iosystem

      ! Build 2D PIO decomposition
      allocate(idof(remap_data%n_b_local))
      do i = 1, remap_data%n_b_local
         gcol = remap_data%row_start + i - 1
         idof(i) = int(mod(gcol-1, remap_data%nlon) + 1, PIO_OFFSET_KIND) + &
              int(remap_data%nlon, PIO_OFFSET_KIND) * &
              int((gcol-1) / remap_data%nlon, PIO_OFFSET_KIND)
      end do
      call pio_initdecomp(pio_sys, PIO_REAL, &
           (/remap_data%nlon, remap_data%nlat/), idof, iodesc_2d)
      call pio_initdecomp(pio_sys, PIO_DOUBLE, &
           (/remap_data%nlon, remap_data%nlat/), idof, iodesc_2d_dbl)
      deallocate(idof)
      iodesc_2d_valid = .true.

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
   end subroutine seaice_fme_horiz_remap_init!}}}

   subroutine seaice_fme_horiz_remap_field(domain, field1d, nCells, fld_out, err)!{{{
   !> nCells is informational only; field1d is assumed-size so the send list
   !> can safely access halo cells beyond nCellsSolve.
   !> Fill values (1.0E34) are detected and masked for correct interpolation.
   !> Target cells with no valid sources get SHR_FILL_VALUE.
      type(domain_type), intent(inout) :: domain
      integer, intent(in) :: nCells
      real(kind=RKIND), intent(in) :: field1d(*)
      real(r8), intent(out) :: fld_out(remap_data%n_b_local)
      integer, intent(out) :: err
      real(r8) :: fld_out_2d(remap_data%n_b_local, 1)
      integer :: i, nlev_packed, needed
      real(r8) :: val

      ! Fill detection: 10% of MPAS fill value (catches 1e34 and variants)
      real(r8), parameter :: fill_detect = MPAS_SI_FILL * 0.1_r8

      err = 0
      if (.not. is_initialized) then; err = 1; return; end if

      nlev_packed = 2  ! 1 field level + 1 mask level

      needed = remap_data%n_send_total * nlev_packed
      if (.not. allocated(ws_send_buf) .or. size(ws_send_buf) < max(1, needed)) then
         if (allocated(ws_send_buf)) deallocate(ws_send_buf)
         allocate(ws_send_buf(max(1, needed)))
      end if

      if (remap_data%n_send_total > 0) then
         do i = 1, remap_data%n_send_total
            val = real(field1d(send_local_cell(i)), r8)
            if (abs(val) >= fill_detect) then
               ! Fill value: zero field, mask = 0
               ws_send_buf((i-1)*nlev_packed + 1) = 0.0_r8
               ws_send_buf((i-1)*nlev_packed + 2) = 0.0_r8
            else
               ! Valid: field + mask = 1
               ws_send_buf((i-1)*nlev_packed + 1) = val
               ws_send_buf((i-1)*nlev_packed + 2) = 1.0_r8
            end if
         end do
      end if

      call remap_data%apply_masked(ws_send_buf, 1, fld_out_2d, &
           domain % dminfo % comm, domain % dminfo % nprocs, err)
      fld_out(1:remap_data%n_b_local) = fld_out_2d(1:remap_data%n_b_local, 1)
   end subroutine seaice_fme_horiz_remap_field!}}}

   subroutine seaice_fme_horiz_remap_init_time_epoch(domain, amName, err)!{{{
   !> Mirror of ocn_fme_horiz_remap_init_time_epoch for the seaice side.
   !> Sidecar lookup -> MPAS_START_TIME fallback. See the ocean
   !> module's docstring for the design rationale.
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
      call try_read_epoch_from_sidecar_si(amName, epoch_str, from_sidecar)

      if (.not. from_sidecar) then
         start_time = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr)
         call mpas_get_time(start_time, dateTimeString=epoch_str_buf, ierr=ierr)
         epoch_str = trim(epoch_str_buf)
         if (ierr /= 0) then
            call mpas_log_write('seaice_fme_horiz_remap_init_time_epoch: ' // &
                 'failed to read MPAS_START_TIME from clock', &
                 messageType=MPAS_LOG_CRIT)
            err = ierr
            return
         end if
         if (is_restart_run) then
            call mpas_log_write('seaice_fme_horiz_remap_init_time_epoch: ' // &
                 'warm restart but no sidecar epoch found; falling back ' // &
                 'to MPAS_START_TIME = "' // trim(epoch_str) // '". Likely ' // &
                 'the wrong epoch -- investigate sidecar.', &
                 messageType=MPAS_LOG_WARN)
         end if
      end if

      call set_time_epoch_si(epoch_str, err)

      if (err == 0) then
         if (from_sidecar) then
            call mpas_log_write('seaice_fme_horiz_remap_init_time_epoch: ' // &
                 'epoch restored from sidecar = "' // trim(time_epoch_cf) // '"')
         else
            call mpas_log_write('seaice_fme_horiz_remap_init_time_epoch: ' // &
                 'epoch captured from MPAS_START_TIME = "' // &
                 trim(time_epoch_cf) // '"')
         end if
      end if

   end subroutine seaice_fme_horiz_remap_init_time_epoch!}}}

   subroutine try_read_epoch_from_sidecar_si(amName, epoch_str, success)!{{{
   !> Read time_reference_date global attribute from an existing
   !> sidecar file (if any). success=.false. on missing file or attr.
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

      call seaice_fme_accum_filename(trim(amName), filename)
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

   end subroutine try_read_epoch_from_sidecar_si!}}}

   subroutine set_time_epoch_si(epoch_str, err)!{{{
   !> Internal helper. Stores the epoch in module state. Validates
   !> against any previously-set value.
      character(len=*), intent(in) :: epoch_str
      integer,          intent(out) :: err
      integer :: ierr, idx
      character(len=64) :: cf_str

      err = 0

      if (len_trim(epoch_str) == 0) then
         call mpas_log_write('set_time_epoch_si: empty epoch string', &
              messageType=MPAS_LOG_CRIT)
         err = 1; return
      end if

      if (time_epoch_set) then
         if (trim(epoch_str) /= trim(time_epoch_str)) then
            call mpas_log_write('set_time_epoch_si: epoch already set to "' &
                 // trim(time_epoch_str) // '"; rejecting reset to "' // &
                 trim(epoch_str) // '"', messageType=MPAS_LOG_CRIT)
            err = 2
         end if
         return
      end if

      call mpas_set_time(time_epoch, dateTimeString=trim(epoch_str), ierr=ierr)
      if (ierr /= 0) then
         call mpas_log_write('set_time_epoch_si: failed to parse "' // &
              trim(epoch_str) // '" as YYYY-MM-DD_HH:MM:SS', &
              messageType=MPAS_LOG_CRIT)
         err = ierr; return
      end if

      cf_str = trim(epoch_str)
      idx = index(cf_str, '_')
      if (idx > 0) cf_str(idx:idx) = ' '

      time_epoch_str = trim(epoch_str)
      time_epoch_cf  = trim(cf_str)
      time_epoch_set = .true.
   end subroutine set_time_epoch_si!}}}

   subroutine si_time_to_days_since_epoch(t, days, err)!{{{
   !> Convert an MPAS_Time_Type to days since the configured seaice epoch.
   !> Uses ESMF time-interval arithmetic via mpas_get_timeInterval so the
   !> calendar (noleap) is respected exactly.
      type(MPAS_Time_Type), intent(in) :: t
      real(r8), intent(out) :: days
      integer, intent(out) :: err
      type(MPAS_TimeInterval_Type) :: dt_interval
      real(kind=RKIND) :: dt_seconds
      integer :: ierr

      days = 0.0_r8
      err  = 0
      if (.not. time_epoch_set) then
         call mpas_log_write('si_time_to_days_since_epoch: epoch not set', &
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

   end subroutine si_time_to_days_since_epoch!}}}

   subroutine seaice_fme_horiz_remap_open_file(domain, filename, err)!{{{
   !> Open the remapped output file. If a file with this name already
   !> exists on disk (append-mode after warm restart), reopen it via
   !> pio_openfile and seed time_record from the existing time dim
   !> length. Otherwise create-new.
      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: filename
      integer, intent(out) :: err
      integer :: ierr, i, count_le_now
      integer(PIO_OFFSET_KIND) :: time_len
      logical :: file_exists, may_append
      type(MPAS_Time_Type) :: now_time
      character(len=StrKIND) :: now_str
      real(r8) :: now_days
      real(r8), allocatable :: existing_time(:)

      err = 0
      if (.not. time_epoch_set) then
         call mpas_log_write('seaice_fme_horiz_remap: time epoch not set; ' // &
              'AM init must call seaice_fme_horiz_remap_set_time_epoch first', &
              messageType=MPAS_LOG_CRIT)
         err = 1; return
      end if

      inquire(file=trim(filename), exist=file_exists)
      ! Append-mode reopen ONLY for the leg-start month's file on a warm
      ! restart, gated on a fire-once latch (is_restart_run stays true for
      ! the whole run). A new month, or a stale/0-byte leftover, must
      ! clobber-create instead. See the ocean module for full rationale.
      ! Branch/rewind restart must overwrite the re-crossed month, not
      ! append to it -- fall through to clobber-create (see rewind_restart).
      may_append = file_exists .and. is_restart_run .and. &
           restart_reopen_pending .and. .not. rewind_restart
      restart_reopen_pending = .false.
      if (may_append) then
         ierr = pio_openfile(pio_sys, pio_file, PIO_IOTYPE_PNETCDF, &
              trim(filename), PIO_WRITE)
         if (ierr /= PIO_NOERR) then
            ! Reopen failed (stale/corrupt file): fall through to a clean
            ! clobber-create rather than aborting the run.
            call mpas_log_write('seaice_fme_horiz_remap: failed to reopen ' &
                 // trim(filename) // '; clobber-creating instead', &
                 messageType=MPAS_LOG_WARN)
            pio_file%fh = -1
            go to 200
         end if

         call check_pio_si('inq_dimid lon', &
              pio_inq_dimid(pio_file, 'lon', dim_lon_id), err)
         call check_pio_si('inq_dimid lat', &
              pio_inq_dimid(pio_file, 'lat', dim_lat_id), err)
         call check_pio_si('inq_dimid time', &
              pio_inq_dimid(pio_file, 'time', dim_time_id), err)
         call check_pio_si('inq_dimid nbnd', &
              pio_inq_dimid(pio_file, 'nbnd', dim_nbnd_id), err)
         call check_pio_si('inq_varid lon',       pio_inq_varid(pio_file, &
              'lon',       pio_lon_var),       err)
         call check_pio_si('inq_varid lat',       pio_inq_varid(pio_file, &
              'lat',       pio_lat_var),       err)
         call check_pio_si('inq_varid time',      pio_inq_varid(pio_file, &
              'time',      pio_time_var),      err)
         call check_pio_si('inq_varid time_bnds', pio_inq_varid(pio_file, &
              'time_bnds', pio_time_bnds_var), err)
         call check_pio_si('inq_varid date',      pio_inq_varid(pio_file, &
              'date',      pio_date_var),      err)
         call check_pio_si('inq_varid datesec',   pio_inq_varid(pio_file, &
              'datesec',   pio_datesec_var),   err)
         call check_pio_si('inq_dimlen time', &
              pio_inq_dimlen(pio_file, dim_time_id, time_len), err)

         ! Seed time_record at largest frame whose `time < MPAS_NOW (in days)`
         ! so leg-2's next write overwrites the first frame whose time
         ! >= restart time.
         now_time = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
         call mpas_get_time(now_time, dateTimeString=now_str, ierr=ierr)
         call si_time_to_days_since_epoch(now_time, now_days, ierr)

         allocate(existing_time(int(time_len)))
         ierr = pio_get_var(pio_file, pio_time_var, existing_time)
         if (ierr /= PIO_NOERR) then
            call mpas_log_write('seaice_fme_horiz_remap: failed to read ' // &
                 'time; falling back to N_existing', &
                 messageType=MPAS_LOG_WARN)
            time_record = int(time_len)
         else
            count_le_now = 0
            do i = 1, int(time_len)
               if (existing_time(i) < now_days) then
                  count_le_now = count_le_now + 1
               end if
            end do
            time_record = count_le_now
         end if
         deallocate(existing_time)

         n_pio_vars = 0
         file_is_open = .true.
         existing_file_reopened = .true.
         if (err == 0) then
            call mpas_log_write('seaice_fme_horiz_remap: reopened ' // &
                 trim(filename) // ' in append mode; existing time=$i, ' // &
                 'next write at frame=$i+1 (MPAS_NOW=' // trim(now_str) // ')', &
                 intArgs=(/int(time_len), time_record/))
         end if
         return
      end if

200   continue
      ierr = pio_createfile(pio_sys, pio_file, PIO_IOTYPE_PNETCDF, trim(filename), PIO_CLOBBER)
      if (ierr /= PIO_NOERR) then; err = ierr; return; end if

      call check_pio_si('def_dim lon',  pio_def_dim(pio_file, 'lon',  remap_data%nlon, dim_lon_id), err)
      call check_pio_si('def_dim lat',  pio_def_dim(pio_file, 'lat',  remap_data%nlat, dim_lat_id), err)
      call check_pio_si('def_dim time', pio_def_dim(pio_file, 'time', PIO_UNLIMITED, dim_time_id), err)
      call check_pio_si('def_dim nbnd', pio_def_dim(pio_file, 'nbnd', 2,    dim_nbnd_id), err)

      ! Coordinate variables; descriptors saved at module scope for write_time.
      call check_pio_si('def_var lon', pio_def_var(pio_file, 'lon', PIO_DOUBLE, (/dim_lon_id/), pio_lon_var), err)
      ierr = pio_put_att(pio_file, pio_lon_var, 'units', 'degrees_east')
      ierr = pio_put_att(pio_file, pio_lon_var, 'long_name', 'longitude')
      ierr = pio_put_att(pio_file, pio_lon_var, 'axis', 'X')
      call check_pio_si('def_var lat', pio_def_var(pio_file, 'lat', PIO_DOUBLE, (/dim_lat_id/), pio_lat_var), err)
      ierr = pio_put_att(pio_file, pio_lat_var, 'units', 'degrees_north')
      ierr = pio_put_att(pio_file, pio_lat_var, 'long_name', 'latitude')
      ierr = pio_put_att(pio_file, pio_lat_var, 'axis', 'Y')

      ! CF time / time_bnds / date / datesec -- mirror of ocean module.
      call check_pio_si('def_var time', pio_def_var(pio_file, 'time', PIO_DOUBLE, &
           (/dim_time_id/), pio_time_var), err)
      ierr = pio_put_att(pio_file, pio_time_var, 'long_name', 'time')
      ierr = pio_put_att(pio_file, pio_time_var, 'units', &
           'days since ' // trim(time_epoch_cf))
      ierr = pio_put_att(pio_file, pio_time_var, 'calendar', 'noleap')
      ierr = pio_put_att(pio_file, pio_time_var, 'bounds', 'time_bnds')
      ierr = pio_put_att(pio_file, pio_time_var, 'axis', 'T')

      call check_pio_si('def_var time_bnds', pio_def_var(pio_file, 'time_bnds', PIO_DOUBLE, &
           (/dim_nbnd_id, dim_time_id/), pio_time_bnds_var), err)
      ierr = pio_put_att(pio_file, pio_time_bnds_var, 'long_name', 'time interval endpoints')

      call check_pio_si('def_var date', pio_def_var(pio_file, 'date', PIO_INT, &
           (/dim_time_id/), pio_date_var), err)
      ierr = pio_put_att(pio_file, pio_date_var, 'long_name', 'current date (YYYYMMDD)')

      call check_pio_si('def_var datesec', pio_def_var(pio_file, 'datesec', PIO_INT, &
           (/dim_time_id/), pio_datesec_var), err)
      ierr = pio_put_att(pio_file, pio_datesec_var, 'long_name', &
           'current seconds of current date')

      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'title', &
           'MPAS-Seaice FME training tape')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'source', 'E3SM FME online remapping')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'calendar', 'noleap')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'time_reference_date', trim(time_epoch_cf))
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'realm', 'seaIce')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'product', 'model-output')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'institution_id', 'E3SM-Project')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'contact', 'e3sm-data-support@llnl.gov')
      if (provenance_set) then
         if (len_trim(prov_case)          > 0) ierr = pio_put_att(pio_file, PIO_GLOBAL, &
              'case',          trim(prov_case))
         if (len_trim(prov_username)      > 0) ierr = pio_put_att(pio_file, PIO_GLOBAL, &
              'username',      trim(prov_username))
         if (len_trim(prov_hostname)      > 0) ierr = pio_put_att(pio_file, PIO_GLOBAL, &
              'hostname',      trim(prov_hostname))
         if (len_trim(prov_model_version) > 0) ierr = pio_put_att(pio_file, PIO_GLOBAL, &
              'git_version',   trim(prov_model_version))
         if (len_trim(prov_history)       > 0) ierr = pio_put_att(pio_file, PIO_GLOBAL, &
              'history',       trim(prov_history))
      end if

      if (err /= 0) then
         call mpas_log_write('seaice_fme_horiz_remap: define-mode error opening ' &
              // trim(filename), messageType=MPAS_LOG_ERR)
         ! Close the just-created file so its ncid is not orphaned in
         ! pio_file_list (a recycled orphan id later trips "invalid IO type").
         call pio_closefile(pio_file)
         pio_file%fh = -1
         file_is_open = .false.
         return
      end if

      n_pio_vars = 0
      time_record = 0
      file_is_open = .true.
      existing_file_reopened = .false.
   end subroutine seaice_fme_horiz_remap_open_file!}}}

   subroutine check_pio_si(label, ierr, accumulated_err)
      character(len=*), intent(in) :: label
      integer, intent(in) :: ierr
      integer, intent(inout) :: accumulated_err
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('seaice_fme_horiz_remap: PIO error in ' // trim(label) // ' (code $i)', &
              messageType=MPAS_LOG_ERR, intArgs=(/ierr/))
         if (accumulated_err == 0) accumulated_err = ierr
      end if
   end subroutine check_pio_si

   subroutine seaice_fme_horiz_remap_def_var(varname, units, long_name, numlev, err, &
        cell_methods)!{{{
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(in) :: numlev
      integer, intent(out) :: err
      character(len=*), optional, intent(in) :: cell_methods
      integer :: ierr
      real(kind=4) :: fill_val

      err = 0
      n_pio_vars = n_pio_vars + 1
      if (n_pio_vars > max_pio_vars) then
         call mpas_log_write('seaice_fme_horiz_remap: exceeded max_pio_vars=$i for var: ' &
              // trim(varname), messageType=MPAS_LOG_ERR, &
              intArgs=(/max_pio_vars/))
         n_pio_vars = n_pio_vars - 1
         err = 2; return
      end if

      pio_var_names(n_pio_vars) = varname

      if (existing_file_reopened) then
         ierr = pio_inq_varid(pio_file, trim(varname), pio_var_descs(n_pio_vars))
         if (ierr /= PIO_NOERR) then
            call mpas_log_write('seaice_fme_horiz_remap: pio_inq_varid failed for ' &
                 // trim(varname) // ' (code $i)', messageType=MPAS_LOG_ERR, &
                 intArgs=(/ierr/))
            err = ierr; return
         end if
         return
      end if

      ierr = pio_def_var(pio_file, trim(varname), PIO_REAL, &
           (/dim_lon_id, dim_lat_id, dim_time_id/), pio_var_descs(n_pio_vars))
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('seaice_fme_horiz_remap: pio_def_var failed for ' &
              // trim(varname) // ' (code $i)', messageType=MPAS_LOG_ERR, &
              intArgs=(/ierr/))
         err = ierr; return
      end if
      ierr = pio_put_att(pio_file, pio_var_descs(n_pio_vars), 'units', trim(units))
      ierr = pio_put_att(pio_file, pio_var_descs(n_pio_vars), 'long_name', trim(long_name))
      if (present(cell_methods)) then
         if (len_trim(cell_methods) > 0) then
            ierr = pio_put_att(pio_file, pio_var_descs(n_pio_vars), &
                 'cell_methods', trim(cell_methods))
         end if
      end if
      fill_val = real(SHR_FILL_VALUE, 4)
      ierr = pio_put_att(pio_file, pio_var_descs(n_pio_vars), '_FillValue', fill_val)
      ierr = pio_put_att(pio_file, pio_var_descs(n_pio_vars), 'missing_value', fill_val)
   end subroutine seaice_fme_horiz_remap_def_var!}}}

   subroutine seaice_fme_horiz_remap_write_var(varname, fld_out, numlev, err)!{{{
      character(len=*), intent(in) :: varname
      real(r8), intent(in) :: fld_out(:)
      integer, intent(in) :: numlev
      integer, intent(out) :: err
      integer :: iv, ierr
      ! Must match the PIO_REAL type used in pio_initdecomp for iodesc_2d.
      ! Passing real(r8) to a PIO_REAL decomposition garbles the output.
      real(4) :: fld_out_r4(remap_data%n_b_local)

      err = 0
      do iv = 1, n_pio_vars
         if (trim(pio_var_names(iv)) == trim(varname)) then
            call pio_setframe(pio_file, pio_var_descs(iv), &
                 int(time_record, PIO_OFFSET_KIND))
            fld_out_r4(1:remap_data%n_b_local) = &
                 real(fld_out(1:remap_data%n_b_local), 4)
            call pio_write_darray(pio_file, pio_var_descs(iv), iodesc_2d, &
                 fld_out_r4, ierr)
            err = ierr
            return
         end if
      end do
      err = 3
   end subroutine seaice_fme_horiz_remap_write_var!}}}

   subroutine seaice_fme_horiz_remap_write_time(time_lo, time_hi, err)!{{{
   !> Mirror of the ocean module's remap_file_write_time.
   !> See that subroutine's comment for the CF-time scheme.
      type(MPAS_Time_Type), intent(in) :: time_lo, time_hi
      integer, intent(out) :: err
      integer :: ierr, yr, mo, dy, hr, mn, sc
      real(r8) :: t_lo_days, t_hi_days
      real(r8) :: time_buf(1)
      real(r8) :: bnds_buf(2, 1)
      integer  :: date_buf(1), datesec_buf(1)

      err = 0
      if (.not. file_is_open) then; err = 1; return; end if

      ! Seed the lon/lat coordinate axes on the first record of a file.
      ! Fresh-create files must first leave define mode (pio_enddef);
      ! reopened (append) files are already in data mode.
      !
      ! time_record == 0 on a REOPENED file means a prior leg created the
      ! file's schema but closed it before its first flush, so lon/lat were
      ! never written (they live only in this block). Re-seed them so the
      ! reopened file does not ship fill-valued coordinate axes (the same
      ! created-but-empty-file case as the ocean static-mask self-heal --
      ! AGENTS.md gotcha #45). A normal warm-restart reopen of a file that
      ! already holds records seeds time_record >= 1 and skips this block,
      ! leaving its coordinates untouched (warm-restart BFB).
      if (time_record == 0) then
         if (.not. existing_file_reopened) then
            ierr = pio_enddef(pio_file)
         end if
         ierr = pio_put_var(pio_file, pio_lon_var, remap_data%lon(1:remap_data%nlon))
         ierr = pio_put_var(pio_file, pio_lat_var, remap_data%lat(1:remap_data%nlat))
      end if
      time_record = time_record + 1

      call si_time_to_days_since_epoch(time_lo, t_lo_days, ierr)
      call si_time_to_days_since_epoch(time_hi, t_hi_days, ierr)

      time_buf(1) = t_hi_days
      bnds_buf(1, 1) = t_lo_days
      bnds_buf(2, 1) = t_hi_days

      call mpas_get_time(time_hi, YYYY=yr, MM=mo, DD=dy, &
           H=hr, M=mn, S=sc, ierr=ierr)
      date_buf(1)    = yr * 10000 + mo * 100 + dy
      datesec_buf(1) = hr * 3600  + mn * 60  + sc

      ierr = pio_put_var(pio_file, pio_time_var, &
           (/time_record/), (/1/), time_buf)
      ierr = pio_put_var(pio_file, pio_time_bnds_var, &
           (/1, time_record/), (/2, 1/), bnds_buf)
      ierr = pio_put_var(pio_file, pio_date_var, &
           (/time_record/), (/1/), date_buf)
      ierr = pio_put_var(pio_file, pio_datesec_var, &
           (/time_record/), (/1/), datesec_buf)
   end subroutine seaice_fme_horiz_remap_write_time!}}}

   subroutine seaice_fme_horiz_remap_close_file(err)!{{{
      integer, intent(out) :: err
      err = 0
      if (file_is_open) then
         call pio_closefile(pio_file)
         file_is_open = .false.
         current_filename = ''
         existing_file_reopened = .false.
         ! Drop the dead ncid and stale dim ids so a leftover handle can't
         ! alias a recycled ncid on the next open.
         pio_file%fh  = -1
         dim_lon_id   = 0
         dim_lat_id   = 0
         dim_time_id  = 0
         dim_nbnd_id  = 0
         n_pio_vars   = 0
         time_record  = 0
      end if
   end subroutine seaice_fme_horiz_remap_close_file!}}}

   subroutine seaice_fme_horiz_remap_register_var(varname, units, long_name, err, &
        cell_methods)!{{{
   !> Silently skips if already registered (MPAS may call init twice).
   !> cell_methods is stored verbatim and replayed by check_rotate into
   !> def_var, which writes the CF `cell_methods` attribute.
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(out) :: err
      character(len=*), optional, intent(in) :: cell_methods
      integer :: iv

      err = 0

      ! Skip if already registered (guard against double-init)
      do iv = 1, n_registered_vars
         if (trim(pio_var_names(iv)) == trim(varname)) return
      end do

      n_registered_vars = n_registered_vars + 1
      if (n_registered_vars > max_pio_vars) then
         call mpas_log_write('seaice_fme_horiz_remap: exceeded max_pio_vars=$i for var: ' &
              // trim(varname), messageType=MPAS_LOG_ERR, &
              intArgs=(/max_pio_vars/))
         n_registered_vars = n_registered_vars - 1
         err = 2; return
      end if
      pio_var_names(n_registered_vars) = varname
      stored_units(n_registered_vars) = units
      stored_long_names(n_registered_vars) = long_name
      stored_cell_methods(n_registered_vars) = ''
      if (present(cell_methods)) stored_cell_methods(n_registered_vars) = cell_methods
   end subroutine seaice_fme_horiz_remap_register_var!}}}

   subroutine seaice_fme_horiz_remap_check_rotate(domain, err)!{{{
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=512) :: new_filename
      integer :: ierr, iv

      err = 0
      if (n_registered_vars == 0) return
      if (len_trim(remap_stream_name) == 0) return

      call seaice_fme_remap_build_filename(domain, &
           trim(remap_stream_name), new_filename, ierr)
      if (ierr /= 0) then; err = ierr; return; end if

      if (.not. file_is_open .or. &
           trim(new_filename) /= trim(current_filename)) then

         if (file_is_open) then
            call seaice_fme_horiz_remap_close_file(ierr)
         end if

         call seaice_fme_horiz_remap_open_file(domain, trim(new_filename), ierr)
         if (ierr /= 0) then; err = ierr; return; end if
         current_filename = new_filename

         do iv = 1, n_registered_vars
            call seaice_fme_horiz_remap_def_var(trim(pio_var_names(iv)), &
                 trim(stored_units(iv)), trim(stored_long_names(iv)), 1, ierr, &
                 cell_methods=trim(stored_cell_methods(iv)))
         end do
      end if

   end subroutine seaice_fme_horiz_remap_check_rotate!}}}

   subroutine seaice_fme_horiz_remap_get_dims(nlat, nlon, n_b_local)!{{{
      integer, intent(out) :: nlat, nlon, n_b_local
      nlat = remap_data%nlat; nlon = remap_data%nlon; n_b_local = remap_data%n_b_local
   end subroutine seaice_fme_horiz_remap_get_dims!}}}

   logical function seaice_fme_horiz_remap_active()!{{{
      seaice_fme_horiz_remap_active = is_initialized
   end function seaice_fme_horiz_remap_active!}}}

   subroutine seaice_fme_remap_build_filename(domain, streamName, filename, err)!{{{
   !> Build a date-stamped remapped output filename from the native stream's
   !> filename_template (set by buildnml, contains casename and date tokens).

      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: streamName
      character(len=*), intent(out) :: filename
      integer, intent(out) :: err

      character(len=StrKIND) :: tmpl
      type(MPAS_Time_Type) :: currTime
      integer :: yr, mo, dy, ierr, idx
      character(len=4) :: yr_str
      character(len=2) :: mo_str, dy_str
      character(len=len(filename)) :: tmp

      err = 0
      filename = ''

      call MPAS_stream_mgr_get_property(domain % streamManager, &
           trim(streamName), MPAS_STREAM_PROPERTY_FILENAME, tmpl, ierr=ierr)
      if (ierr /= 0) then
         tmpl = trim(streamName) // '.$Y-$M-$D.nc'
      end if

      ! Insert ".remapped" before ".nc"
      idx = index(tmpl, '.nc', back=.true.)
      if (idx > 0) then
         tmpl = tmpl(1:idx-1) // '.remapped.nc'
      else
         tmpl = trim(tmpl) // '.remapped.nc'
      end if

      ! Substitute date from current model time
      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(currTime, YYYY=yr, MM=mo, DD=dy, ierr=ierr)

      write(yr_str, '(I4.4)') yr
      write(mo_str, '(I2.2)') mo
      write(dy_str, '(I2.2)') dy

      filename = tmpl
      idx = index(filename, '$Y')
      if (idx > 0) then; tmp = filename(1:idx-1) // yr_str // filename(idx+2:); filename = tmp; end if
      idx = index(filename, '$M')
      if (idx > 0) then; tmp = filename(1:idx-1) // mo_str // filename(idx+2:); filename = tmp; end if
      idx = index(filename, '$D')
      if (idx > 0) then; tmp = filename(1:idx-1) // dy_str // filename(idx+2:); filename = tmp; end if

   end subroutine seaice_fme_remap_build_filename!}}}

!***********************************************************************

   subroutine seaice_fme_accum_filename(amName, filename)!{{{
      character(len=*), intent(in)  :: amName
      character(len=*), intent(out) :: filename
      filename = trim(amName) // '.fme_accum_restart.nc'
   end subroutine seaice_fme_accum_filename!}}}

!***********************************************************************

   subroutine seaice_fme_horiz_remap_write_accum(amName, accum_arr, &
        n_avg_fields, nAccum, avg_last_output_time_str, count_arr, err)!{{{
   !> Persist remap-space accumulator state to a sidecar file. count_arr
   !> is the per-cell valid-sample counter parallel to accum_arr — passing
   !> it preserves the time-averaging mask across restart so the warm-leg
   !> output continues to weight intensive fill-prone fields correctly
   !> (see AGENTS.md gotcha on time-averaging fill pollution).
      character(len=*), intent(in) :: amName
      integer,          intent(in) :: n_avg_fields
      real(r8),         intent(in) :: accum_arr(:,:)
      integer,          intent(in) :: nAccum
      character(len=*), intent(in) :: avg_last_output_time_str
      real(r8), optional, intent(in) :: count_arr(:,:)
      integer,          intent(out) :: err

      type(file_desc_t)        :: pio_acc_file
      type(var_desc_t)         :: var_accum, var_count
      integer                  :: dim_lon, dim_lat, dim_field
      integer                  :: ierr, idx
      character(len=512)       :: filename
      real(r8), allocatable    :: slab(:)

      err = 0
      if (.not. is_initialized) return

      call seaice_fme_accum_filename(trim(amName), filename)

      ierr = pio_createfile(pio_sys, pio_acc_file, PIO_IOTYPE_PNETCDF, &
           trim(filename), PIO_CLOBBER)
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('seaice_fme_write_accum: failed to create ' &
              // trim(filename), messageType=MPAS_LOG_WARN)
         err = ierr; return
      end if

      call check_pio_si('accum def_dim lon', &
           pio_def_dim(pio_acc_file, 'lon', remap_data%nlon, dim_lon), err)
      call check_pio_si('accum def_dim lat', &
           pio_def_dim(pio_acc_file, 'lat', remap_data%nlat, dim_lat), err)
      ! Record dim — must be UNLIMITED so pio_setframe is valid.
      call check_pio_si('accum def_dim n_avg_fields', &
           pio_def_dim(pio_acc_file, 'n_avg_fields', PIO_UNLIMITED, dim_field), err)
      call check_pio_si('accum def_var remap_accum', &
           pio_def_var(pio_acc_file, 'remap_accum', PIO_DOUBLE, &
                       (/dim_lon, dim_lat, dim_field/), var_accum), err)
      if (present(count_arr)) then
         call check_pio_si('accum def_var remap_accum_count', &
              pio_def_var(pio_acc_file, 'remap_accum_count', PIO_DOUBLE, &
                          (/dim_lon, dim_lat, dim_field/), var_count), err)
      end if

      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, 'nAccum', nAccum)
      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, 'n_avg_fields', n_avg_fields)
      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, &
           'avg_last_output_time', trim(avg_last_output_time_str))
      ! Persist the CF time epoch -- mirror of ocean module. Lets
      ! init_time_epoch recover the case-start time on warm restart
      ! (MPAS_START_TIME from the clock would return leg-start instead).
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
   end subroutine seaice_fme_horiz_remap_write_accum!}}}

!***********************************************************************

   subroutine seaice_fme_horiz_remap_set_rewind()!{{{
   !> Flag a branch/rewind restart so the file-open paths clobber-create
   !> the re-crossed monthly files instead of append-reopening them. Called
   !> by the seaice AM's *_read_restart when the persisted schedule time is
   !> later than MPAS_NOW.
      rewind_restart = .true.
      call mpas_log_write('seaice_fme_horiz_remap: rewind restart flagged; ' // &
           're-crossed monthly files will be clobber-created.')
   end subroutine seaice_fme_horiz_remap_set_rewind!}}}

!***********************************************************************

   subroutine seaice_fme_horiz_remap_read_accum(amName, accum_arr, &
        n_avg_fields, nAccum, avg_last_output_time_str, count_arr, success)!{{{
      character(len=*),   intent(in)    :: amName
      integer,            intent(in)    :: n_avg_fields
      real(r8),           intent(inout) :: accum_arr(:,:)
      integer,            intent(out)   :: nAccum
      character(len=*),   intent(out)   :: avg_last_output_time_str
      real(r8), optional, intent(inout) :: count_arr(:,:)
      logical,            intent(out)   :: success

      type(file_desc_t)     :: pio_acc_file
      type(var_desc_t)      :: var_accum, var_count
      integer               :: ierr, idx, n_avg_on_disk
      character(len=512)    :: filename
      logical               :: file_exists
      real(r8), allocatable :: slab(:)
      real(r8), allocatable :: tmp_accum(:,:), tmp_count(:,:)

      success = .false.
      nAccum = 0
      avg_last_output_time_str = ''
      if (.not. is_initialized) return

      call seaice_fme_accum_filename(trim(amName), filename)
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) return

      ierr = pio_openfile(pio_sys, pio_acc_file, PIO_IOTYPE_PNETCDF, &
           trim(filename), PIO_WRITE)
      if (ierr /= PIO_NOERR) then
         call mpas_log_write('seaice_fme_read_accum: pio_openfile failed for ' &
              // trim(filename), messageType=MPAS_LOG_WARN)
         return
      end if

      ierr = pio_get_att(pio_acc_file, PIO_GLOBAL, 'n_avg_fields', n_avg_on_disk)
      if (ierr /= PIO_NOERR .or. n_avg_on_disk /= n_avg_fields) then
         call mpas_log_write('seaice_fme_read_accum: n_avg_fields mismatch ' &
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
      call mpas_log_write('seaice_fme_read_accum: restored accumulators from ' &
           // trim(filename) // ' (nAccum=$i)', intArgs=(/nAccum/))

   end subroutine seaice_fme_horiz_remap_read_accum!}}}

!***********************************************************************

   subroutine seaice_fme_horiz_remap_set_provenance(case_name, username, hostname, &
        model_version, history_str)!{{{
   !> Capture CIME-side provenance once during ice_comp_mct init. The
   !> strings are written verbatim as global attributes onto every
   !> freshly-created `.remapped.nc` file (see seaice_fme_horiz_remap_open_file).
   !> No-op on subsequent calls so the FIRST values win.
      character(len=*), intent(in) :: case_name, username, hostname
      character(len=*), intent(in) :: model_version, history_str

      if (provenance_set) return
      prov_case          = case_name
      prov_username      = username
      prov_hostname      = hostname
      prov_model_version = model_version
      prov_history       = history_str
      provenance_set     = .true.
   end subroutine seaice_fme_horiz_remap_set_provenance!}}}

!***********************************************************************
!  5D companion procedures.
!
!  Parallel set of register/open/def/write/close procedures that operate
!  on the `*_5D` module-level state declared near the top of this module.
!  Mirrors the primary procedures line-for-line apart from the file
!  context. The accumulator sidecar routines (write_accum/read_accum)
!  are NOT duplicated -- they are amName-keyed and the AM passes a
!  distinct amName ('fmeSeaiceDerivedFields5D') for the 5D path.
!***********************************************************************

   subroutine seaice_fme_horiz_remap_register_var_5D(varname, units, long_name, &
        err, cell_methods)!{{{
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(out) :: err
      character(len=*), optional, intent(in) :: cell_methods
      integer :: iv

      err = 0

      do iv = 1, n_registered_vars_5D
         if (trim(pio_var_names_5D(iv)) == trim(varname)) return
      end do

      n_registered_vars_5D = n_registered_vars_5D + 1
      if (n_registered_vars_5D > max_pio_vars) then
         call mpas_log_write('seaice_fme_horiz_remap_5D: exceeded max_pio_vars=$i for var: ' &
              // trim(varname), messageType=MPAS_LOG_ERR, &
              intArgs=(/max_pio_vars/))
         n_registered_vars_5D = n_registered_vars_5D - 1
         err = 2; return
      end if
      pio_var_names_5D(n_registered_vars_5D) = varname
      stored_units_5D(n_registered_vars_5D) = units
      stored_long_names_5D(n_registered_vars_5D) = long_name
      stored_cell_methods_5D(n_registered_vars_5D) = ''
      if (present(cell_methods)) stored_cell_methods_5D(n_registered_vars_5D) = cell_methods
   end subroutine seaice_fme_horiz_remap_register_var_5D!}}}

   subroutine seaice_fme_horiz_remap_open_file_5D(domain, filename, err)!{{{
      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: filename
      integer, intent(out) :: err
      integer :: ierr, i, count_le_now
      integer(PIO_OFFSET_KIND) :: time_len
      logical :: file_exists, may_append
      type(MPAS_Time_Type) :: now_time
      character(len=StrKIND) :: now_str
      real(r8) :: now_days
      real(r8), allocatable :: existing_time(:)

      err = 0
      if (.not. time_epoch_set) then
         call mpas_log_write('seaice_fme_horiz_remap_5D: time epoch not set', &
              messageType=MPAS_LOG_CRIT)
         err = 1; return
      end if

      inquire(file=trim(filename), exist=file_exists)
      ! Fire-once restart-append latch (see primary open_file rationale).
      ! Branch/rewind restart must overwrite the re-crossed month, not
      ! append to it -- fall through to clobber-create (see rewind_restart).
      may_append = file_exists .and. is_restart_run .and. &
           restart_reopen_pending_5D .and. .not. rewind_restart
      restart_reopen_pending_5D = .false.
      if (may_append) then
         ierr = pio_openfile(pio_sys, pio_file_5D, PIO_IOTYPE_PNETCDF, &
              trim(filename), PIO_WRITE)
         if (ierr /= PIO_NOERR) then
            ! Reopen failed (stale/corrupt file): fall through to clobber.
            call mpas_log_write('seaice_fme_horiz_remap_5D: failed to reopen ' &
                 // trim(filename) // '; clobber-creating instead', &
                 messageType=MPAS_LOG_WARN)
            pio_file_5D%fh = -1
            go to 300
         end if

         call check_pio_si('5D inq_dimid lon', &
              pio_inq_dimid(pio_file_5D, 'lon', dim_lon_id_5D), err)
         call check_pio_si('5D inq_dimid lat', &
              pio_inq_dimid(pio_file_5D, 'lat', dim_lat_id_5D), err)
         call check_pio_si('5D inq_dimid time', &
              pio_inq_dimid(pio_file_5D, 'time', dim_time_id_5D), err)
         call check_pio_si('5D inq_dimid nbnd', &
              pio_inq_dimid(pio_file_5D, 'nbnd', dim_nbnd_id_5D), err)
         call check_pio_si('5D inq_varid lon', pio_inq_varid(pio_file_5D, &
              'lon', pio_lon_var_5D), err)
         call check_pio_si('5D inq_varid lat', pio_inq_varid(pio_file_5D, &
              'lat', pio_lat_var_5D), err)
         call check_pio_si('5D inq_varid time', pio_inq_varid(pio_file_5D, &
              'time', pio_time_var_5D), err)
         call check_pio_si('5D inq_varid time_bnds', pio_inq_varid(pio_file_5D, &
              'time_bnds', pio_time_bnds_var_5D), err)
         call check_pio_si('5D inq_varid date', pio_inq_varid(pio_file_5D, &
              'date', pio_date_var_5D), err)
         call check_pio_si('5D inq_varid datesec', pio_inq_varid(pio_file_5D, &
              'datesec', pio_datesec_var_5D), err)
         call check_pio_si('5D inq_dimlen time', &
              pio_inq_dimlen(pio_file_5D, dim_time_id_5D, time_len), err)

         now_time = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
         call mpas_get_time(now_time, dateTimeString=now_str, ierr=ierr)
         call si_time_to_days_since_epoch(now_time, now_days, ierr)

         allocate(existing_time(int(time_len)))
         ierr = pio_get_var(pio_file_5D, pio_time_var_5D, existing_time)
         if (ierr /= PIO_NOERR) then
            time_record_5D = int(time_len)
         else
            count_le_now = 0
            do i = 1, int(time_len)
               if (existing_time(i) < now_days) count_le_now = count_le_now + 1
            end do
            time_record_5D = count_le_now
         end if
         deallocate(existing_time)

         n_pio_vars_5D = 0
         file_is_open_5D = .true.
         existing_file_reopened_5D = .true.
         return
      end if

300   continue
      ierr = pio_createfile(pio_sys, pio_file_5D, PIO_IOTYPE_PNETCDF, &
           trim(filename), PIO_CLOBBER)
      if (ierr /= PIO_NOERR) then; err = ierr; return; end if

      call check_pio_si('5D def_dim lon',  pio_def_dim(pio_file_5D, 'lon',  remap_data%nlon, dim_lon_id_5D), err)
      call check_pio_si('5D def_dim lat',  pio_def_dim(pio_file_5D, 'lat',  remap_data%nlat, dim_lat_id_5D), err)
      call check_pio_si('5D def_dim time', pio_def_dim(pio_file_5D, 'time', PIO_UNLIMITED, dim_time_id_5D), err)
      call check_pio_si('5D def_dim nbnd', pio_def_dim(pio_file_5D, 'nbnd', 2, dim_nbnd_id_5D), err)

      call check_pio_si('5D def_var lon', pio_def_var(pio_file_5D, 'lon', PIO_DOUBLE, &
           (/dim_lon_id_5D/), pio_lon_var_5D), err)
      ierr = pio_put_att(pio_file_5D, pio_lon_var_5D, 'units', 'degrees_east')
      ierr = pio_put_att(pio_file_5D, pio_lon_var_5D, 'long_name', 'longitude')
      ierr = pio_put_att(pio_file_5D, pio_lon_var_5D, 'axis', 'X')
      call check_pio_si('5D def_var lat', pio_def_var(pio_file_5D, 'lat', PIO_DOUBLE, &
           (/dim_lat_id_5D/), pio_lat_var_5D), err)
      ierr = pio_put_att(pio_file_5D, pio_lat_var_5D, 'units', 'degrees_north')
      ierr = pio_put_att(pio_file_5D, pio_lat_var_5D, 'long_name', 'latitude')
      ierr = pio_put_att(pio_file_5D, pio_lat_var_5D, 'axis', 'Y')

      call check_pio_si('5D def_var time', pio_def_var(pio_file_5D, 'time', PIO_DOUBLE, &
           (/dim_time_id_5D/), pio_time_var_5D), err)
      ierr = pio_put_att(pio_file_5D, pio_time_var_5D, 'long_name', 'time')
      ierr = pio_put_att(pio_file_5D, pio_time_var_5D, 'units', &
           'days since ' // trim(time_epoch_cf))
      ierr = pio_put_att(pio_file_5D, pio_time_var_5D, 'calendar', 'noleap')
      ierr = pio_put_att(pio_file_5D, pio_time_var_5D, 'bounds', 'time_bnds')
      ierr = pio_put_att(pio_file_5D, pio_time_var_5D, 'axis', 'T')

      call check_pio_si('5D def_var time_bnds', pio_def_var(pio_file_5D, 'time_bnds', PIO_DOUBLE, &
           (/dim_nbnd_id_5D, dim_time_id_5D/), pio_time_bnds_var_5D), err)
      ierr = pio_put_att(pio_file_5D, pio_time_bnds_var_5D, 'long_name', 'time interval endpoints')

      call check_pio_si('5D def_var date', pio_def_var(pio_file_5D, 'date', PIO_INT, &
           (/dim_time_id_5D/), pio_date_var_5D), err)
      ierr = pio_put_att(pio_file_5D, pio_date_var_5D, 'long_name', 'current date (YYYYMMDD)')

      call check_pio_si('5D def_var datesec', pio_def_var(pio_file_5D, 'datesec', PIO_INT, &
           (/dim_time_id_5D/), pio_datesec_var_5D), err)
      ierr = pio_put_att(pio_file_5D, pio_datesec_var_5D, 'long_name', &
           'current seconds of current date')

      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'title', &
           'MPAS-Seaice FME training tape (5-day mean)')
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'source', 'E3SM FME online remapping')
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'calendar', 'noleap')
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'time_reference_date', trim(time_epoch_cf))
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'realm', 'seaIce')
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'product', 'model-output')
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'institution_id', 'E3SM-Project')
      ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, 'contact', 'e3sm-data-support@llnl.gov')
      if (provenance_set) then
         if (len_trim(prov_case)          > 0) ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, &
              'case',          trim(prov_case))
         if (len_trim(prov_username)      > 0) ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, &
              'username',      trim(prov_username))
         if (len_trim(prov_hostname)      > 0) ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, &
              'hostname',      trim(prov_hostname))
         if (len_trim(prov_model_version) > 0) ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, &
              'git_version',   trim(prov_model_version))
         if (len_trim(prov_history)       > 0) ierr = pio_put_att(pio_file_5D, PIO_GLOBAL, &
              'history',       trim(prov_history))
      end if

      if (err /= 0) then
         call mpas_log_write('seaice_fme_horiz_remap_5D: define-mode error opening ' &
              // trim(filename), messageType=MPAS_LOG_ERR)
         ! Close the just-created file so its ncid is not orphaned.
         call pio_closefile(pio_file_5D)
         pio_file_5D%fh = -1
         file_is_open_5D = .false.
         return
      end if

      n_pio_vars_5D = 0
      time_record_5D = 0
      file_is_open_5D = .true.
      existing_file_reopened_5D = .false.
   end subroutine seaice_fme_horiz_remap_open_file_5D!}}}

   subroutine seaice_fme_horiz_remap_def_var_5D(varname, units, long_name, numlev, err, &
        cell_methods)!{{{
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(in) :: numlev
      integer, intent(out) :: err
      character(len=*), optional, intent(in) :: cell_methods
      integer :: ierr
      real(kind=4) :: fill_val

      err = 0
      n_pio_vars_5D = n_pio_vars_5D + 1
      if (n_pio_vars_5D > max_pio_vars) then
         n_pio_vars_5D = n_pio_vars_5D - 1
         err = 2; return
      end if

      pio_var_names_5D(n_pio_vars_5D) = varname

      if (existing_file_reopened_5D) then
         ierr = pio_inq_varid(pio_file_5D, trim(varname), pio_var_descs_5D(n_pio_vars_5D))
         if (ierr /= PIO_NOERR) then; err = ierr; return; end if
         return
      end if

      ierr = pio_def_var(pio_file_5D, trim(varname), PIO_REAL, &
           (/dim_lon_id_5D, dim_lat_id_5D, dim_time_id_5D/), pio_var_descs_5D(n_pio_vars_5D))
      if (ierr /= PIO_NOERR) then; err = ierr; return; end if
      ierr = pio_put_att(pio_file_5D, pio_var_descs_5D(n_pio_vars_5D), 'units', trim(units))
      ierr = pio_put_att(pio_file_5D, pio_var_descs_5D(n_pio_vars_5D), 'long_name', trim(long_name))
      if (present(cell_methods)) then
         if (len_trim(cell_methods) > 0) then
            ierr = pio_put_att(pio_file_5D, pio_var_descs_5D(n_pio_vars_5D), &
                 'cell_methods', trim(cell_methods))
         end if
      end if
      fill_val = real(SHR_FILL_VALUE, 4)
      ierr = pio_put_att(pio_file_5D, pio_var_descs_5D(n_pio_vars_5D), '_FillValue', fill_val)
      ierr = pio_put_att(pio_file_5D, pio_var_descs_5D(n_pio_vars_5D), 'missing_value', fill_val)
   end subroutine seaice_fme_horiz_remap_def_var_5D!}}}

   subroutine seaice_fme_horiz_remap_write_var_5D(varname, fld_out, numlev, err)!{{{
      character(len=*), intent(in) :: varname
      real(r8), intent(in) :: fld_out(:)
      integer, intent(in) :: numlev
      integer, intent(out) :: err
      integer :: iv, ierr
      real(4) :: fld_out_r4(remap_data%n_b_local)

      err = 0
      do iv = 1, n_pio_vars_5D
         if (trim(pio_var_names_5D(iv)) == trim(varname)) then
            call pio_setframe(pio_file_5D, pio_var_descs_5D(iv), &
                 int(time_record_5D, PIO_OFFSET_KIND))
            fld_out_r4(1:remap_data%n_b_local) = &
                 real(fld_out(1:remap_data%n_b_local), 4)
            call pio_write_darray(pio_file_5D, pio_var_descs_5D(iv), iodesc_2d, &
                 fld_out_r4, ierr)
            err = ierr
            return
         end if
      end do
      err = 3
   end subroutine seaice_fme_horiz_remap_write_var_5D!}}}

   subroutine seaice_fme_horiz_remap_write_time_5D(time_lo, time_hi, err)!{{{
      type(MPAS_Time_Type), intent(in) :: time_lo, time_hi
      integer, intent(out) :: err
      integer :: ierr, yr, mo, dy, hr, mn, sc
      real(r8) :: t_lo_days, t_hi_days
      real(r8) :: time_buf(1)
      real(r8) :: bnds_buf(2, 1)
      integer  :: date_buf(1), datesec_buf(1)

      err = 0
      if (.not. file_is_open_5D) then; err = 1; return; end if

      ! Seed lon/lat on the first record. Re-seed on a reopened file with no
      ! records yet (prior leg created the schema but never flushed, so the
      ! coordinate axes were left unwritten -- see write_time above and
      ! AGENTS.md gotcha #45). A reopen of a file that already holds records
      ! seeds time_record_5D >= 1 and skips this block (warm-restart BFB).
      if (time_record_5D == 0) then
         if (.not. existing_file_reopened_5D) then
            ierr = pio_enddef(pio_file_5D)
         end if
         ierr = pio_put_var(pio_file_5D, pio_lon_var_5D, remap_data%lon(1:remap_data%nlon))
         ierr = pio_put_var(pio_file_5D, pio_lat_var_5D, remap_data%lat(1:remap_data%nlat))
      end if
      time_record_5D = time_record_5D + 1

      call si_time_to_days_since_epoch(time_lo, t_lo_days, ierr)
      call si_time_to_days_since_epoch(time_hi, t_hi_days, ierr)

      time_buf(1) = t_hi_days
      bnds_buf(1, 1) = t_lo_days
      bnds_buf(2, 1) = t_hi_days

      call mpas_get_time(time_hi, YYYY=yr, MM=mo, DD=dy, &
           H=hr, M=mn, S=sc, ierr=ierr)
      date_buf(1)    = yr * 10000 + mo * 100 + dy
      datesec_buf(1) = hr * 3600  + mn * 60  + sc

      ierr = pio_put_var(pio_file_5D, pio_time_var_5D, &
           (/time_record_5D/), (/1/), time_buf)
      ierr = pio_put_var(pio_file_5D, pio_time_bnds_var_5D, &
           (/1, time_record_5D/), (/2, 1/), bnds_buf)
      ierr = pio_put_var(pio_file_5D, pio_date_var_5D, &
           (/time_record_5D/), (/1/), date_buf)
      ierr = pio_put_var(pio_file_5D, pio_datesec_var_5D, &
           (/time_record_5D/), (/1/), datesec_buf)
   end subroutine seaice_fme_horiz_remap_write_time_5D!}}}

   subroutine seaice_fme_horiz_remap_close_file_5D(err)!{{{
      integer, intent(out) :: err
      err = 0
      if (file_is_open_5D) then
         call pio_closefile(pio_file_5D)
         file_is_open_5D = .false.
         current_filename_5D = ''
         existing_file_reopened_5D = .false.
         ! Drop dead ncid / stale dim ids (see primary close).
         pio_file_5D%fh = -1
         dim_lon_id_5D  = 0
         dim_lat_id_5D  = 0
         dim_time_id_5D = 0
         dim_nbnd_id_5D = 0
         n_pio_vars_5D  = 0
         time_record_5D = 0
      end if
   end subroutine seaice_fme_horiz_remap_close_file_5D!}}}

   subroutine seaice_fme_horiz_remap_check_rotate_5D(domain, err)!{{{
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err
      character(len=512) :: new_filename
      integer :: ierr, iv

      err = 0
      if (n_registered_vars_5D == 0) return
      if (len_trim(remap_stream_name_5D) == 0) return

      call seaice_fme_remap_build_filename(domain, &
           trim(remap_stream_name_5D), new_filename, ierr)
      if (ierr /= 0) then; err = ierr; return; end if

      if (.not. file_is_open_5D .or. &
           trim(new_filename) /= trim(current_filename_5D)) then

         if (file_is_open_5D) then
            call seaice_fme_horiz_remap_close_file_5D(ierr)
         end if

         call seaice_fme_horiz_remap_open_file_5D(domain, trim(new_filename), ierr)
         if (ierr /= 0) then; err = ierr; return; end if
         current_filename_5D = new_filename

         do iv = 1, n_registered_vars_5D
            call seaice_fme_horiz_remap_def_var_5D(trim(pio_var_names_5D(iv)), &
                 trim(stored_units_5D(iv)), trim(stored_long_names_5D(iv)), 1, ierr, &
                 cell_methods=trim(stored_cell_methods_5D(iv)))
         end do
      end if
   end subroutine seaice_fme_horiz_remap_check_rotate_5D!}}}

   subroutine seaice_fme_horiz_remap_write_accum_5D(amName, accum_arr, &
        n_avg_fields, nAccum, avg_last_output_time_str, count_arr, err)!{{{
   !> Thin alias around the primary write_accum -- the routine is
   !> amName-keyed and shares iodesc_2d_dbl across both channels, so
   !> all that's needed is a distinct sidecar filename via amName.
      character(len=*), intent(in) :: amName
      integer,          intent(in) :: n_avg_fields
      real(r8),         intent(in) :: accum_arr(:,:)
      integer,          intent(in) :: nAccum
      character(len=*), intent(in) :: avg_last_output_time_str
      real(r8), optional, intent(in) :: count_arr(:,:)
      integer,          intent(out) :: err
      if (present(count_arr)) then
         call seaice_fme_horiz_remap_write_accum(amName, accum_arr, n_avg_fields, &
              nAccum, avg_last_output_time_str, count_arr=count_arr, err=err)
      else
         call seaice_fme_horiz_remap_write_accum(amName, accum_arr, n_avg_fields, &
              nAccum, avg_last_output_time_str, err=err)
      end if
   end subroutine seaice_fme_horiz_remap_write_accum_5D!}}}

   subroutine seaice_fme_horiz_remap_read_accum_5D(amName, accum_arr, &
        n_avg_fields, nAccum, avg_last_output_time_str, count_arr, success)!{{{
      character(len=*),   intent(in)    :: amName
      integer,            intent(in)    :: n_avg_fields
      real(r8),           intent(inout) :: accum_arr(:,:)
      integer,            intent(out)   :: nAccum
      character(len=*),   intent(out)   :: avg_last_output_time_str
      real(r8), optional, intent(inout) :: count_arr(:,:)
      logical,            intent(out)   :: success
      if (present(count_arr)) then
         call seaice_fme_horiz_remap_read_accum(amName, accum_arr, n_avg_fields, &
              nAccum, avg_last_output_time_str, count_arr=count_arr, success=success)
      else
         call seaice_fme_horiz_remap_read_accum(amName, accum_arr, n_avg_fields, &
              nAccum, avg_last_output_time_str, success=success)
      end if
   end subroutine seaice_fme_horiz_remap_read_accum_5D!}}}

end module seaice_fme_horiz_remap
