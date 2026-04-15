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
                  pio_createfile, pio_closefile, &
                  pio_def_dim, pio_def_var, pio_enddef, &
                  pio_put_att, pio_put_var, pio_write_darray, &
                  pio_initdecomp, pio_setframe, &
                  PIO_IOTYPE_PNETCDF, PIO_CLOBBER, PIO_REAL, PIO_DOUBLE, &
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
   public :: remap_stream_name, file_is_open

   type(shr_horiz_remap_t) :: remap_data
   integer, allocatable :: send_local_cell(:)
   real(r8), allocatable :: ws_send_buf(:)
   logical :: is_initialized = .false.

   ! MPAS-SI fill value for input field detection (matches def_var _FillValue)
   real(r8), parameter :: MPAS_SI_FILL = SHR_FILL_VALUE

   ! PIO output
   type(file_desc_t) :: pio_file
   type(iosystem_desc_t), pointer :: pio_sys => null()
   logical :: file_is_open = .false.
   integer :: dim_lon_id, dim_lat_id, dim_time_id
   integer :: time_record = 0
   type(io_desc_t) :: iodesc_2d
   logical :: iodesc_2d_valid = .false.

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
      deallocate(idof)
      iodesc_2d_valid = .true.

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
      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: filename
      integer, intent(out) :: err
      integer :: ierr

      err = 0
      ierr = pio_createfile(pio_sys, pio_file, PIO_IOTYPE_PNETCDF, trim(filename), PIO_CLOBBER)
      if (ierr /= PIO_NOERR) then; err = ierr; return; end if

      ierr = pio_def_dim(pio_file, 'lon', remap_data%nlon, dim_lon_id)
      ierr = pio_def_dim(pio_file, 'lat', remap_data%nlat, dim_lat_id)
      ierr = pio_def_dim(pio_file, 'Time', PIO_UNLIMITED, dim_time_id)
      ierr = pio_def_dim(pio_file, 'StrLen', 64, dim_strlen_id)

      ! Define coordinate variables; descriptors saved at module scope for write_time
      ierr = pio_def_var(pio_file, 'lon', PIO_DOUBLE, (/dim_lon_id/), pio_lon_var)
      ierr = pio_put_att(pio_file, pio_lon_var, 'units', 'degrees_east')
      ierr = pio_put_att(pio_file, pio_lon_var, 'long_name', 'longitude')
      ierr = pio_put_att(pio_file, pio_lon_var, 'axis', 'X')
      ierr = pio_def_var(pio_file, 'lat', PIO_DOUBLE, (/dim_lat_id/), pio_lat_var)
      ierr = pio_put_att(pio_file, pio_lat_var, 'units', 'degrees_north')
      ierr = pio_put_att(pio_file, pio_lat_var, 'long_name', 'latitude')
      ierr = pio_put_att(pio_file, pio_lat_var, 'axis', 'Y')

      ierr = pio_def_var(pio_file, 'xtime', PIO_CHAR, &
           (/dim_strlen_id, dim_time_id/), pio_xtime_var)
      ierr = pio_put_att(pio_file, pio_xtime_var, 'long_name', 'model time')

      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'source', 'E3SM FME online remapping')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'Conventions', 'CF-1.8')
      ierr = pio_put_att(pio_file, PIO_GLOBAL, 'calendar', 'noleap')

      n_pio_vars = 0
      time_record = 0
      file_is_open = .true.
   end subroutine seaice_fme_horiz_remap_open_file!}}}

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
      ierr = pio_def_var(pio_file, trim(varname), PIO_REAL, &
           (/dim_lon_id, dim_lat_id, dim_time_id/), pio_var_descs(n_pio_vars))
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
      if (time_record == 0) then
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

end module seaice_fme_horiz_remap
