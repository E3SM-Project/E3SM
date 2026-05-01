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
                  PIO_REAL, PIO_DOUBLE, &
                  PIO_CHAR, PIO_OFFSET_KIND, PIO_GLOBAL, PIO_NOERR, &
                  PIO_UNLIMITED

   implicit none
   private
   save

   public :: seaice_fme_horiz_remap_init
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

   type(shr_horiz_remap_t) :: remap_data
   integer, allocatable :: send_local_cell(:)
   real(r8), allocatable :: ws_send_buf(:)
   logical :: is_initialized = .false.

   ! Cached at init from config_do_restart. Append-mode only fires on
   ! warm restart; cold start always clobber-creates, defending against
   ! stale `.remapped.nc` left behind by a prior CIME compare-copy.
   logical :: is_restart_run = .false.

   ! MPAS-SI fill value for input field detection (matches def_var _FillValue)
   real(r8), parameter :: MPAS_SI_FILL = SHR_FILL_VALUE

   ! PIO output
   type(file_desc_t) :: pio_file
   type(iosystem_desc_t), pointer :: pio_sys => null()
   logical :: file_is_open = .false.
   integer :: dim_lon_id, dim_lat_id, dim_time_id
   integer :: time_record = 0
   type(io_desc_t) :: iodesc_2d        ! PIO_REAL  (output records)
   type(io_desc_t) :: iodesc_2d_dbl    ! PIO_DOUBLE (accumulator restart)
   logical :: iodesc_2d_valid = .false.
   ! True when the open file was reopened from disk (append after restart)
   ! rather than freshly created. Mirrors the ocean module's flag.
   logical :: existing_file_reopened = .false.

   integer, parameter :: max_pio_vars = 512
   integer :: n_pio_vars = 0
   character(len=64) :: pio_var_names(max_pio_vars)
   type(var_desc_t) :: pio_var_descs(max_pio_vars)
   type(var_desc_t) :: pio_lon_var, pio_lat_var, pio_xtime_var
   integer :: dim_strlen_id = 0

   ! Stored var metadata for replay on rotation
   integer :: n_registered_vars = 0
   character(len=64)  :: stored_units(max_pio_vars)
   character(len=128) :: stored_long_names(max_pio_vars)
   character(len=512) :: current_filename = ''
   character(len=128) :: remap_stream_name = ''

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

   subroutine seaice_fme_horiz_remap_open_file(domain, filename, err)!{{{
   !> Open the remapped output file. If a file with this name already
   !> exists on disk (append-mode after warm restart), reopen it via
   !> pio_openfile and seed time_record from the existing Time dim
   !> length. Otherwise create-new.
      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: filename
      integer, intent(out) :: err
      integer :: ierr, i, count_le_now
      integer(PIO_OFFSET_KIND) :: time_len
      logical :: file_exists
      type(MPAS_Time_Type) :: now_time
      character(len=64) :: now_str
      character(len=64), allocatable :: existing_xtime(:)

      err = 0

      inquire(file=trim(filename), exist=file_exists)
      ! Append-mode reopen ONLY on warm restart. Cold start always
      ! clobber-creates, defending against stale leftover files.
      if (file_exists .and. is_restart_run) then
         ierr = pio_openfile(pio_sys, pio_file, PIO_IOTYPE_PNETCDF, &
              trim(filename), PIO_WRITE)
         if (ierr /= PIO_NOERR) then
            call mpas_log_write('seaice_fme_horiz_remap: failed to reopen ' &
                 // trim(filename), messageType=MPAS_LOG_ERR)
            err = ierr; return
         end if

         call check_pio_si('inq_dimid lon', &
              pio_inq_dimid(pio_file, 'lon', dim_lon_id), err)
         call check_pio_si('inq_dimid lat', &
              pio_inq_dimid(pio_file, 'lat', dim_lat_id), err)
         call check_pio_si('inq_dimid Time', &
              pio_inq_dimid(pio_file, 'Time', dim_time_id), err)
         call check_pio_si('inq_dimid StrLen', &
              pio_inq_dimid(pio_file, 'StrLen', dim_strlen_id), err)
         call check_pio_si('inq_varid lon',   pio_inq_varid(pio_file, &
              'lon', pio_lon_var), err)
         call check_pio_si('inq_varid lat',   pio_inq_varid(pio_file, &
              'lat', pio_lat_var), err)
         call check_pio_si('inq_varid xtime', pio_inq_varid(pio_file, &
              'xtime', pio_xtime_var), err)
         call check_pio_si('inq_dimlen Time', &
              pio_inq_dimlen(pio_file, dim_time_id, time_len), err)

         ! Mirror of ocean: seed time_record at largest frame whose
         ! xtime <= MPAS_NOW so leg-2's next write overwrites the
         ! first frame whose xtime > restart time.
         now_time = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
         call mpas_get_time(now_time, dateTimeString=now_str, ierr=ierr)

         allocate(existing_xtime(int(time_len)))
         ierr = pio_get_var(pio_file, pio_xtime_var, existing_xtime)
         if (ierr /= PIO_NOERR) then
            call mpas_log_write('seaice_fme_horiz_remap: failed to read ' // &
                 'xtime; falling back to N_existing', &
                 messageType=MPAS_LOG_WARN)
            time_record = int(time_len)
         else
            ! Strict <: check_rotate is called from the trigger-fire
            ! branch, so a frame with xtime == MPAS_NOW is leg-1's
            ! leftover that we want to overwrite. Excluding it makes
            ! next write_time land on it.
            count_le_now = 0
            do i = 1, int(time_len)
               if (trim(existing_xtime(i)) < trim(now_str)) then
                  count_le_now = count_le_now + 1
               end if
            end do
            time_record = count_le_now
         end if
         deallocate(existing_xtime)

         n_pio_vars = 0
         file_is_open = .true.
         existing_file_reopened = .true.
         if (err == 0) then
            call mpas_log_write('seaice_fme_horiz_remap: reopened ' // &
                 trim(filename) // ' in append mode; existing Time=$i, ' // &
                 'next write at frame=$i+1 (MPAS_NOW=' // trim(now_str) // ')', &
                 intArgs=(/int(time_len), time_record/))
         end if
         return
      end if

      ierr = pio_createfile(pio_sys, pio_file, PIO_IOTYPE_PNETCDF, trim(filename), PIO_CLOBBER)
      if (ierr /= PIO_NOERR) then; err = ierr; return; end if

      call check_pio_si('def_dim lon', pio_def_dim(pio_file, 'lon', remap_data%nlon, dim_lon_id), err)
      call check_pio_si('def_dim lat', pio_def_dim(pio_file, 'lat', remap_data%nlat, dim_lat_id), err)
      call check_pio_si('def_dim Time', pio_def_dim(pio_file, 'Time', PIO_UNLIMITED, dim_time_id), err)
      call check_pio_si('def_dim StrLen', pio_def_dim(pio_file, 'StrLen', 64, dim_strlen_id), err)

      ! Define coordinate variables; descriptors saved at module scope for write_time
      call check_pio_si('def_var lon', pio_def_var(pio_file, 'lon', PIO_DOUBLE, (/dim_lon_id/), pio_lon_var), err)
      ierr = pio_put_att(pio_file, pio_lon_var, 'units', 'degrees_east')
      ierr = pio_put_att(pio_file, pio_lon_var, 'long_name', 'longitude')
      ierr = pio_put_att(pio_file, pio_lon_var, 'axis', 'X')
      call check_pio_si('def_var lat', pio_def_var(pio_file, 'lat', PIO_DOUBLE, (/dim_lat_id/), pio_lat_var), err)
      ierr = pio_put_att(pio_file, pio_lat_var, 'units', 'degrees_north')
      ierr = pio_put_att(pio_file, pio_lat_var, 'long_name', 'latitude')
      ierr = pio_put_att(pio_file, pio_lat_var, 'axis', 'Y')

      call check_pio_si('def_var xtime', pio_def_var(pio_file, 'xtime', PIO_CHAR, &
           (/dim_strlen_id, dim_time_id/), pio_xtime_var), err)
      ierr = pio_put_att(pio_file, pio_xtime_var, 'long_name', 'model time')

      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'source', 'E3SM FME online remapping')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'calendar', 'noleap')

      if (err /= 0) then
         call mpas_log_write('seaice_fme_horiz_remap: define-mode error opening ' &
              // trim(filename), messageType=MPAS_LOG_ERR)
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

   subroutine seaice_fme_horiz_remap_def_var(varname, units, long_name, numlev, err)!{{{
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(in) :: numlev
      integer, intent(out) :: err
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

   subroutine seaice_fme_horiz_remap_write_time(time_str, err)!{{{
      character(len=*), intent(in) :: time_str
      integer, intent(out) :: err
      integer :: ierr
      character(len=64) :: xtime_buf(1)

      err = 0
      ! Fresh-create files: leave define-mode and seed lon/lat once.
      ! Reopened (append) files are already in data mode and lon/lat are
      ! already populated, so this whole block is skipped.
      if (time_record == 0 .and. .not. existing_file_reopened) then
         ierr = pio_enddef(pio_file)
         ierr = pio_put_var(pio_file, pio_lon_var, remap_data%lon(1:remap_data%nlon))
         ierr = pio_put_var(pio_file, pio_lat_var, remap_data%lat(1:remap_data%nlat))
      end if
      time_record = time_record + 1

      ! Write xtime string for this record
      xtime_buf(1) = ' '
      xtime_buf(1)(1:min(len_trim(time_str), 64)) = &
           time_str(1:min(len_trim(time_str), 64))
      ierr = pio_put_var(pio_file, pio_xtime_var, &
           (/1, time_record/), (/64, 1/), xtime_buf)
   end subroutine seaice_fme_horiz_remap_write_time!}}}

   subroutine seaice_fme_horiz_remap_close_file(err)!{{{
      integer, intent(out) :: err
      err = 0
      if (file_is_open) then
         call pio_closefile(pio_file)
         file_is_open = .false.
         current_filename = ''
         existing_file_reopened = .false.
      end if
   end subroutine seaice_fme_horiz_remap_close_file!}}}

   subroutine seaice_fme_horiz_remap_register_var(varname, units, long_name, err)!{{{
   !> Silently skips if already registered (MPAS may call init twice).
      character(len=*), intent(in) :: varname, units, long_name
      integer, intent(out) :: err
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
                 trim(stored_units(iv)), trim(stored_long_names(iv)), 1, ierr)
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
        n_avg_fields, nAccum, avg_last_output_time_str, err)!{{{
   !> Persist remap-space accumulator state to a sidecar file. Mirror of
   !> ocn_fme_horiz_remap_write_accum; sea-ice has no instantaneous-companion
   !> arrays, so the optional latest_inst slot is omitted here.
      character(len=*), intent(in) :: amName
      integer,          intent(in) :: n_avg_fields
      real(r8),         intent(in) :: accum_arr(:,:)
      integer,          intent(in) :: nAccum
      character(len=*), intent(in) :: avg_last_output_time_str
      integer,          intent(out) :: err

      type(file_desc_t)        :: pio_acc_file
      type(var_desc_t)         :: var_accum
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

      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, 'nAccum', nAccum)
      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, 'n_avg_fields', n_avg_fields)
      ierr = pio_put_att(pio_acc_file, PIO_GLOBAL, &
           'avg_last_output_time', trim(avg_last_output_time_str))
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
      end do
      deallocate(slab)
      call pio_closefile(pio_acc_file)
   end subroutine seaice_fme_horiz_remap_write_accum!}}}

!***********************************************************************

   subroutine seaice_fme_horiz_remap_read_accum(amName, accum_arr, &
        n_avg_fields, nAccum, avg_last_output_time_str, success)!{{{
      character(len=*),   intent(in)    :: amName
      integer,            intent(in)    :: n_avg_fields
      real(r8),           intent(inout) :: accum_arr(:,:)
      integer,            intent(out)   :: nAccum
      character(len=*),   intent(out)   :: avg_last_output_time_str
      logical,            intent(out)   :: success

      type(file_desc_t)     :: pio_acc_file
      type(var_desc_t)      :: var_accum
      integer               :: ierr, idx, n_avg_on_disk
      character(len=512)    :: filename
      logical               :: file_exists
      real(r8), allocatable :: slab(:)
      real(r8), allocatable :: tmp_accum(:,:)

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

      accum_arr(:,:) = tmp_accum(:,:)
      deallocate(tmp_accum, slab)
      call pio_closefile(pio_acc_file)

      success = .true.
      call mpas_log_write('seaice_fme_read_accum: restored accumulators from ' &
           // trim(filename) // ' (nAccum=$i)', intArgs=(/nAccum/))

   end subroutine seaice_fme_horiz_remap_read_accum!}}}

end module seaice_fme_horiz_remap
