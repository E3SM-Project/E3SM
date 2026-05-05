! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_fme_depth_coarsening
!
!> \brief MPAS ocean analysis member: FME depth coarsening
!> \author FME Training Data Support
!> \date   April 2026
!> \details
!>  Computes layer-thickness-weighted vertical averages of tracer and
!>  velocity fields onto a coarsened set of depth levels defined by
!>  user-specified depth boundaries. Uses shr_vcoarsen_avg from the
!>  shared vertical coarsening utility module.
!>
!-----------------------------------------------------------------------

module ocn_fme_depth_coarsening

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_kind_types

   use ocn_constants
   use ocn_diagnostics_variables
   use shr_vcoarsen_mod, only: shr_vcoarsen_avg
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_horiz_remap_mod, only: SHR_FILL_VALUE
   use mpas_log, only: mpas_log_write
   use ocn_fme_horiz_remap, only: ocn_fme_horiz_remap_init, &
        ocn_fme_horiz_remap_init_time_epoch, &
        ocn_fme_horiz_remap_field, ocn_fme_horiz_remap_active, &
        ocn_fme_horiz_remap_get_dims, ocn_fme_remap_file_t, &
        ocn_fme_remap_build_filename, &
        ocn_fme_horiz_remap_write_accum, ocn_fme_horiz_remap_read_accum

   implicit none
   private
   save

   public :: ocn_init_fme_depth_coarsening, &
             ocn_compute_fme_depth_coarsening, &
             ocn_restart_fme_depth_coarsening, &
             ocn_finalize_fme_depth_coarsening, &
             ocn_fme_depth_coarsening_write_restart

   integer :: nCoarsenLevels = 0
   real(r8), allocatable :: depthBounds(:)

   !> Number of fields requested via config_AM_fmeDepthCoarsening_fields
   integer :: nCoarsenFields = 0

   !> Field names parsed from comma-separated config string
   character(len=64), allocatable :: coarsenFieldNames(:)

   real(kind=RKIND), parameter :: fillValue = real(SHR_FILL_VALUE, RKIND)

   !> Horizontal remapping state
   logical :: do_horiz_remap = .false.
   type(ocn_fme_remap_file_t) :: remap_outfile
   integer :: remap_nlat = 0, remap_nlon = 0, remap_nb_local = 0

   !> Time-averaging state for remapped output
   logical :: time_avg_enabled = .false.
   integer :: nAccum = 0
   type(MPAS_TimeInterval_type) :: avg_output_interval
   type(MPAS_Time_type) :: avg_last_output_time
   integer :: n_avg_fields = 0       ! nCoarsenLevels * nCoarsenFields + nCoarsenLevels
   real(r8), allocatable :: remap_accum(:,:)        ! (remap_nb_local, n_avg_fields)
   ! Per-cell, per-(field,level) valid-sample counter parallel to
   ! remap_accum. Defensive: depth-coarsened fields use a per-layer
   ! bathymetry mask that is static, so under current usage every
   ! target cell is fill on every sample or on none. The per-cell
   ! count keeps the averaging path correct if that ever changes
   ! (e.g., layer thickness becoming time-varying via free-surface
   ! coupling). See AGENTS.md gotcha on time-averaging fill pollution.
   real(r8), allocatable :: remap_accum_count(:,:)  ! (remap_nb_local, n_avg_fields)

   !> Instantaneous companion: when time_averaging AND write_instantaneous_companion
   !> are both true, mirror the latest remapped value of each (field, level) into
   !> remap_latest_inst so it can be written alongside the averaged value at each
   !> output trigger.
   logical :: write_inst_companion = .false.
   real(r8), allocatable :: remap_latest_inst(:,:)  ! (remap_nb_local, n_avg_fields)

   ! Static per-level bathymetry masks in remap space, written as
   ! `mask_0` ... `mask_{nCoarsenLevels-1}` in every output file. 1.0 on
   ! cells where the depth-coarsened layer is between surface and seafloor,
   ! 0.0 elsewhere. Computed once lazily after layerThicknessCoarsened is
   ! filled. Matches the SamudrACE convention (compute_ocean_dataset_e3sm.py:
   ! mask_<level_index> for any *_<level_index> 3D variable).
   real(r8), allocatable :: mask_levels_remap(:,:)  ! (remap_nb_local, nCoarsenLevels)
   ! Single-slot reusable buffer for write_var calls; copied from
   ! mask_levels_remap(:, k) per level.
   real(r8), allocatable :: mask_write_buf(:,:)    ! (remap_nb_local, 1)

   ! Threshold for detecting fill in remapped output (10% of SHR_FILL_VALUE).
   real(r8), parameter :: MASK_FILL_THRESH = real(SHR_FILL_VALUE, r8) * 0.1_r8

!***********************************************************************

contains

!***********************************************************************

   subroutine ocn_init_fme_depth_coarsening(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND), pointer :: config_AM_fmeDepthCoarsening_depth_bounds
      character(len=StrKIND), pointer :: config_AM_fmeDepthCoarsening_fields
      character(len=StrKIND), pointer :: config_AM_fmeDepthCoarsening_horiz_remap_file
      logical, pointer :: config_time_averaging
      logical, pointer :: config_write_inst_companion
      character(len=StrKIND), pointer :: config_stream_output_interval
      character(len=StrKIND) :: tempStr
      character(len=128) :: varName, longName
      character(len=16)  :: cm   ! CF cell_methods choice
      integer :: nBounds, i, iStart, iEnd, strLen, ierr, iField, kOut

      err = 0

      ! --- Parse depth bounds ---
      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeDepthCoarsening_depth_bounds', &
           config_AM_fmeDepthCoarsening_depth_bounds)

      tempStr = trim(config_AM_fmeDepthCoarsening_depth_bounds)
      if (trim(tempStr) == 'none') then
         nCoarsenLevels = 0
         return
      end if

      ! Count commas to determine number of bounds
      nBounds = 1
      strLen = len_trim(tempStr)
      do i = 1, strLen
         if (tempStr(i:i) == ',') nBounds = nBounds + 1
      end do

      if (allocated(depthBounds)) deallocate(depthBounds)
      allocate(depthBounds(nBounds))

      ! Parse comma-separated values with error handling
      iStart = 1
      do i = 1, nBounds
         iEnd = index(tempStr(iStart:strLen), ',')
         if (iEnd == 0) then
            iEnd = strLen + 1
         else
            iEnd = iEnd + iStart - 1
         end if
         read(tempStr(iStart:iEnd-1), *, iostat=ierr) depthBounds(i)
         if (ierr /= 0) then
            call mpas_log_write('fmeDepthCoarsening: failed to parse depth bound ' // &
                 'at position $i from string: ' // trim(tempStr), &
                 messageType=MPAS_LOG_CRIT, intArgs=(/i/))
            err = 1; return
         end if
         iStart = iEnd + 1
      end do

      ! Validate: bounds must be non-negative and strictly increasing
      do i = 1, nBounds
         if (depthBounds(i) < 0.0_RKIND) then
            call mpas_log_write('fmeDepthCoarsening: depth bound $i is negative: $r', &
                 messageType=MPAS_LOG_CRIT, intArgs=(/i/), &
                 realArgs=(/depthBounds(i)/))
            err = 1; return
         end if
         if (i > 1) then
            if (depthBounds(i) <= depthBounds(i-1)) then
               call mpas_log_write('fmeDepthCoarsening: depth bounds not strictly ' // &
                    'increasing at position $i', &
                    messageType=MPAS_LOG_CRIT, intArgs=(/i/))
               err = 1; return
            end if
         end if
      end do

      nCoarsenLevels = nBounds - 1

      ! --- Parse field names for generic dispatch ---
      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeDepthCoarsening_fields', &
           config_AM_fmeDepthCoarsening_fields)

      tempStr = trim(config_AM_fmeDepthCoarsening_fields)
      if (allocated(coarsenFieldNames)) deallocate(coarsenFieldNames)
      if (trim(tempStr) == 'none' .or. len_trim(tempStr) == 0) then
         ! Default: compute all four fields
         nCoarsenFields = 4
         allocate(coarsenFieldNames(4))
         coarsenFieldNames(1) = 'temperature'
         coarsenFieldNames(2) = 'salinity'
         coarsenFieldNames(3) = 'velocityZonal'
         coarsenFieldNames(4) = 'velocityMeridional'
         call mpas_log_write('fmeDepthCoarsening: no fields specified, ' &
              // 'defaulting to temperature,salinity,velocityZonal,velocityMeridional')
      else
         ! Count comma-separated field names
         nCoarsenFields = 1
         strLen = len_trim(tempStr)
         do i = 1, strLen
            if (tempStr(i:i) == ',') nCoarsenFields = nCoarsenFields + 1
         end do

         allocate(coarsenFieldNames(nCoarsenFields))

         ! Parse comma-separated field names
         iStart = 1
         do i = 1, nCoarsenFields
            iEnd = index(tempStr(iStart:strLen), ',')
            if (iEnd == 0) then
               iEnd = strLen + 1
            else
               iEnd = iEnd + iStart - 1
            end if
            coarsenFieldNames(i) = adjustl(tempStr(iStart:iEnd-1))
            iStart = iEnd + 1
         end do

         call mpas_log_write('fmeDepthCoarsening: configured to coarsen $i fields', &
              intArgs=(/nCoarsenFields/))
      end if

      ! --- Initialise horizontal remapping (optional) ---
      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeDepthCoarsening_horiz_remap_file', &
           config_AM_fmeDepthCoarsening_horiz_remap_file)

      if (trim(config_AM_fmeDepthCoarsening_horiz_remap_file) /= 'none' .and. &
          len_trim(config_AM_fmeDepthCoarsening_horiz_remap_file) > 0) then

         call ocn_fme_horiz_remap_init(domain, &
              trim(config_AM_fmeDepthCoarsening_horiz_remap_file), ierr)
         if (ierr /= 0) then
            call mpas_log_write( &
                 'fmeDepthCoarsening: horiz remap init failed, disabling', &
                 messageType=MPAS_LOG_WARN)
            do_horiz_remap = .false.
         else
            ! Auto-discover the CF time epoch -- sidecar lookup first
            ! (warm restart), MPAS_START_TIME fallback (cold start).
            ! Idempotent across the four FME AMs; only the first caller
            ! does work, the rest see time_epoch_set=.true. and no-op.
            call ocn_fme_horiz_remap_init_time_epoch(domain, &
                 'fmeDepthCoarsening', ierr)
            if (ierr /= 0) then
               call mpas_log_write('fmeDepthCoarsening: time epoch setup ' // &
                    'failed, disabling horiz remap', &
                    messageType=MPAS_LOG_WARN)
               do_horiz_remap = .false.
               err = ierr
               return
            end if
            do_horiz_remap = ocn_fme_horiz_remap_active()
            if (do_horiz_remap) then
               call ocn_fme_horiz_remap_get_dims(remap_nlat, remap_nlon, remap_nb_local)
               call mpas_log_write('fmeDepthCoarsening: horiz remap active, ' &
                    // 'target grid $i x $i', &
                    intArgs=(/remap_nlon, remap_nlat/))
               remap_outfile%stream_name = 'fmeDepthCoarseningOutput'
               ! Pick CF cell_methods up-front based on runtime mode.
               ! Production: time_averaging=.true. ('time: mean'); inst-only
               ! tags 'time: point'. Depth-coarsening is a thickness-weighted
               ! mean and unit-preserving, so output units match the source
               ! field's MPAS-O Registry units.
               call mpas_pool_get_config(domain % configs, &
                    'config_AM_fmeDepthCoarsening_time_averaging', &
                    config_time_averaging)
               cm = 'time: point'
               if (associated(config_time_averaging)) then
                  if (config_time_averaging) cm = 'time: mean'
               end if
               do iField = 1, nCoarsenFields
                  do kOut = 0, nCoarsenLevels - 1
                     write(varName, '(A,A,I0)') &
                          trim(coarsenFieldNames(iField)), 'Coarsened_', kOut
                     write(longName, '(A,I0)') 'Depth-coarsened ' &
                          // trim(coarsenFieldNames(iField)) // ' level ', kOut
                     call remap_outfile%register_var(trim(varName), &
                          trim(coarsen_field_units(coarsenFieldNames(iField))), &
                          trim(longName), ierr, cell_methods=trim(cm))
                  end do
               end do
               do kOut = 0, nCoarsenLevels - 1
                  write(varName, '(A,I0)') 'layerThicknessCoarsened_', kOut
                  write(longName, '(A,I0)') &
                       'Depth-coarsened layer thickness level ', kOut
                  call remap_outfile%register_var(trim(varName), 'm', &
                       trim(longName), ierr, cell_methods=trim(cm))
               end do
               ! Static per-level wet masks: shape (lon, lat) only, 1=valid
               ! ocean cell at this depth, 0=below seafloor or land. Written
               ! once per file. Names match the SamudrACE convention used
               ! by compute_ocean_dataset_e3sm.py for any var ending `_<k>`.
               do kOut = 0, nCoarsenLevels - 1
                  write(varName, '(A,I0)') 'mask_', kOut
                  write(longName, '(A,I0,A)') &
                       'Wet mask for depth level ', kOut, &
                       ' (1=valid ocean, 0=below seafloor or land)'
                  call remap_outfile%register_var(trim(varName), '1', &
                       trim(longName), ierr, is_static=.true.)
               end do
               ! Per-layer-interface depth scalars idepth_0..idepth_N
               ! (N = nCoarsenLevels; one more than the number of layers).
               ! Layer k spans idepth_k -> idepth_{k+1}; thickness =
               ! idepth_{k+1} - idepth_k matches layerThicknessCoarsened_k.
               do kOut = 0, nCoarsenLevels
                  write(varName, '(A,I0)') 'idepth_', kOut
                  write(longName, '(A,I0,A)') &
                       'Depth at coarsened layer interface ', kOut, &
                       ' (m below sea surface)'
                  call remap_outfile%register_var(trim(varName), 'm', &
                       trim(longName), ierr, is_scalar=.true., &
                       scalar_value=depthBounds(kOut+1), &
                       positive='down', standard_name='depth')
               end do

               ! Time-averaging setup
               call mpas_pool_get_config(domain % configs, &
                    'config_AM_fmeDepthCoarsening_time_averaging', &
                    config_time_averaging)
               if (associated(config_time_averaging)) then
                  if (config_time_averaging) then
                     call mpas_pool_get_config(domain % configs, &
                          'config_AM_fmeDepthCoarsening_stream_output_interval', &
                          config_stream_output_interval)
                     call mpas_set_timeInterval(avg_output_interval, &
                          timeString=trim(config_stream_output_interval), &
                          ierr=ierr)
                     ! Default baseline is MPAS_NOW so a warm restart
                     ! without a sidecar uses the restart-time as the
                     ! averaging baseline (no catch-up firing). If a
                     ! sidecar accumulator file is present (config_do_restart
                     ! and read succeeds), the persisted state below
                     ! overrides this.
                     avg_last_output_time = mpas_get_clock_time( &
                          domain % clock, MPAS_NOW, ierr)
                     n_avg_fields = nCoarsenLevels * nCoarsenFields &
                          + nCoarsenLevels
                     if (.not. allocated(remap_accum)) then
                        allocate(remap_accum(remap_nb_local, n_avg_fields))
                     end if
                     if (.not. allocated(remap_accum_count)) then
                        allocate(remap_accum_count(remap_nb_local, n_avg_fields))
                     end if
                     remap_accum = 0.0_r8
                     remap_accum_count = 0.0_r8
                     nAccum = 0
                     time_avg_enabled = .true.
                     call mpas_log_write( &
                          'fmeDepthCoarsening: time averaging enabled')

                     ! Optional instantaneous companion alongside avg
                     call mpas_pool_get_config(domain % configs, &
                          'config_AM_fmeDepthCoarsening_write_instantaneous_companion', &
                          config_write_inst_companion)
                     if (associated(config_write_inst_companion)) then
                        if (config_write_inst_companion) then
                           write_inst_companion = .true.
                           if (.not. allocated(remap_latest_inst)) then
                              allocate(remap_latest_inst(remap_nb_local, n_avg_fields))
                           end if
                           remap_latest_inst = 0.0_r8
                           do iField = 1, nCoarsenFields
                              do kOut = 0, nCoarsenLevels - 1
                                 write(varName, '(A,A,I0,A)') &
                                      trim(coarsenFieldNames(iField)), &
                                      'Coarsened_', kOut, '_inst'
                                 write(longName, '(A,I0,A)') &
                                      'Instantaneous companion: depth-coarsened ' &
                                      // trim(coarsenFieldNames(iField)) // &
                                      ' level ', kOut, ' (snapshot at output trigger)'
                                 call remap_outfile%register_var(trim(varName), &
                                      trim(coarsen_field_units(coarsenFieldNames(iField))), &
                                      trim(longName), ierr, cell_methods='time: point')
                              end do
                           end do
                           do kOut = 0, nCoarsenLevels - 1
                              write(varName, '(A,I0,A)') &
                                   'layerThicknessCoarsened_', kOut, '_inst'
                              write(longName, '(A,I0,A)') &
                                   'Instantaneous companion: layer thickness level ', &
                                   kOut, ' (snapshot at output trigger)'
                              call remap_outfile%register_var(trim(varName), 'm', &
                                   trim(longName), ierr, cell_methods='time: point')
                           end do
                           call mpas_log_write( &
                                'fmeDepthCoarsening: instantaneous _inst companion ' // &
                                'variables enabled alongside time-averaged output')
                        end if
                     end if

                     ! Restore accumulator state on warm restart.
                     call ocn_fme_depth_coarsening_read_restart(domain, ierr)
                  end if
               end if
            end if
         end if
      end if

   end subroutine ocn_init_fme_depth_coarsening!}}}

!***********************************************************************

   subroutine ocn_fme_depth_coarsening_read_restart(domain, err)!{{{
   !> If config_do_restart is .true. and a sidecar accumulator file exists,
   !> restore remap_accum (+ remap_latest_inst when present) and reset the
   !> averaging baseline to the persisted avg_last_output_time. Silent
   !> no-op when not a restart run or when no sidecar is present (fresh
   !> phase-1 cold start).
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      logical, pointer :: config_do_restart
      logical          :: ok
      integer          :: ierr_local
      character(len=StrKIND) :: ts_str

      err = 0
      if (.not. time_avg_enabled) return

      call mpas_pool_get_config(domain % configs, 'config_do_restart', &
           config_do_restart)
      if (.not. associated(config_do_restart)) return
      if (.not. config_do_restart) return

      if (write_inst_companion) then
         call ocn_fme_horiz_remap_read_accum('fmeDepthCoarsening', &
              remap_accum, n_avg_fields, nAccum, ts_str, &
              latest_inst_arr=remap_latest_inst, &
              count_arr=remap_accum_count, success=ok)
      else
         call ocn_fme_horiz_remap_read_accum('fmeDepthCoarsening', &
              remap_accum, n_avg_fields, nAccum, ts_str, &
              count_arr=remap_accum_count, success=ok)
      end if

      if (ok .and. len_trim(ts_str) > 0) then
         call mpas_set_time(avg_last_output_time, &
              dateTimeString=trim(ts_str), ierr=ierr_local)
         if (ierr_local /= 0) then
            call mpas_log_write('fmeDepthCoarsening: failed to parse ' // &
                 'persisted avg_last_output_time="' // trim(ts_str) // &
                 '"; falling back to MPAS_NOW baseline.', &
                 messageType=MPAS_LOG_WARN)
            avg_last_output_time = mpas_get_clock_time( &
                 domain % clock, MPAS_NOW, ierr_local)
         end if
      end if

   end subroutine ocn_fme_depth_coarsening_read_restart!}}}

!***********************************************************************

   subroutine ocn_fme_depth_coarsening_write_restart(domain, err)!{{{
   !> Persist accumulator state to the AM's sidecar file. Called from
   !> the MPAS-O driver on the same restart-write trigger that writes
   !> all native MPAS restart streams.
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND) :: ts_str
      integer :: ierr_local

      err = 0
      if (.not. time_avg_enabled) return
      if (.not. allocated(remap_accum)) return

      call mpas_get_time(avg_last_output_time, &
           dateTimeString=ts_str, ierr=ierr_local)

      if (write_inst_companion .and. allocated(remap_latest_inst)) then
         call ocn_fme_horiz_remap_write_accum('fmeDepthCoarsening', &
              remap_accum, n_avg_fields, nAccum, trim(ts_str), &
              latest_inst_arr=remap_latest_inst, &
              count_arr=remap_accum_count, err=err)
      else
         call ocn_fme_horiz_remap_write_accum('fmeDepthCoarsening', &
              remap_accum, n_avg_fields, nAccum, trim(ts_str), &
              count_arr=remap_accum_count, err=err)
      end if

   end subroutine ocn_fme_depth_coarsening_write_restart!}}}

!***********************************************************************

   subroutine ocn_compute_fme_depth_coarsening(domain, timeLevel, err)!{{{

      integer, intent(in) :: timeLevel
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: verticalMeshPool
      type(mpas_pool_type), pointer :: tracersPool
      type(mpas_pool_type), pointer :: amPool

      integer, pointer :: nCellsSolve, nVertLevels
      integer, pointer :: index_temperature, index_salinity
      integer, pointer :: nFmeDepthLevels
      integer :: iCell, k, iField
      integer, dimension(:), pointer :: maxLevelCell

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers
      real(kind=RKIND), dimension(:, :), pointer :: restingThickness
      real(kind=RKIND), dimension(:, :), pointer :: temperatureCoarsened
      real(kind=RKIND), dimension(:, :), pointer :: salinityCoarsened
      real(kind=RKIND), dimension(:, :), pointer :: velocityZonalCoarsened
      real(kind=RKIND), dimension(:, :), pointer :: velocityMeridionalCoarsened
      real(kind=RKIND), dimension(:, :), pointer :: layerThicknessCoarsened

      ! Work arrays for shr_vcoarsen_avg (single column, r8 to match shared module)
      real(r8), allocatable :: coord_iface(:)
      real(r8), allocatable :: field_1d(:)
      real(r8), allocatable :: field_out(:)
      real(kind=RKIND) :: depthTop

      ! Flags for which fields to compute, driven by coarsenFieldNames
      logical :: doTemperature, doSalinity, doVelocityZonal, doVelocityMeridional

      ! Horizontal remap workspace
      integer, pointer :: nCells
      integer :: ifld, kOut, ierr
      real(r8), allocatable :: fld_remap(:,:)
      real(kind=RKIND), allocatable :: field_1lev(:,:)
      integer, allocatable :: levelMask(:)
      character(len=128) :: varName
      type(MPAS_Time_Type) :: currTime
      ! Pointers to coarsened arrays for generic field loop
      real(kind=RKIND), dimension(:,:), pointer :: coarsenedPtr

      err = 0

      if (nCoarsenLevels == 0) return

      ! Determine which fields to compute based on parsed config
      doTemperature = .false.
      doSalinity = .false.
      doVelocityZonal = .false.
      doVelocityMeridional = .false.

      do iField = 1, nCoarsenFields
         select case (trim(coarsenFieldNames(iField)))
         case ('temperature')
            doTemperature = .true.
         case ('salinity')
            doSalinity = .true.
         case ('velocityZonal')
            if (associated(velocityZonal)) then
               doVelocityZonal = .true.
            else
               call mpas_log_write('fmeDepthCoarsening: velocityZonal not available, skipping', &
                    messageType=MPAS_LOG_WARN)
            end if
         case ('velocityMeridional')
            if (associated(velocityMeridional)) then
               doVelocityMeridional = .true.
            else
               call mpas_log_write('fmeDepthCoarsening: velocityMeridional not available, skipping', &
                    messageType=MPAS_LOG_WARN)
            end if
         end select
      end do

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'fmeDepthCoarseningAM', amPool)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nFmeDepthLevels', nFmeDepthLevels)

         ! Guard: nCoarsenLevels must fit in the Registry dimension
         if (nCoarsenLevels > nFmeDepthLevels) then
            call mpas_log_write('fmeDepthCoarsening: nCoarsenLevels ($i) exceeds ' // &
                 'nFmeDepthLevels ($i) -- reduce depth_bounds entries', &
                 messageType=MPAS_LOG_CRIT, &
                 intArgs=(/nCoarsenLevels, nFmeDepthLevels/))
            err = 1; return
         end if

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
         call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(amPool, 'temperatureCoarsened', temperatureCoarsened)
         call mpas_pool_get_array(amPool, 'salinityCoarsened', salinityCoarsened)
         call mpas_pool_get_array(amPool, 'velocityZonalCoarsened', velocityZonalCoarsened)
         call mpas_pool_get_array(amPool, 'velocityMeridionalCoarsened', velocityMeridionalCoarsened)
         call mpas_pool_get_array(amPool, 'layerThicknessCoarsened', layerThicknessCoarsened)

         ! Initialize all output arrays to fillValue; only requested fields
         ! will be overwritten with computed values.
         temperatureCoarsened(:,:) = fillValue
         salinityCoarsened(:,:) = fillValue
         velocityZonalCoarsened(:,:) = fillValue
         velocityMeridionalCoarsened(:,:) = fillValue
         layerThicknessCoarsened(:,:) = fillValue

         allocate(coord_iface(nVertLevels + 1))
         allocate(field_1d(nVertLevels))
         allocate(field_out(nCoarsenLevels))

         do iCell = 1, nCellsSolve

            ! Build depth interface array from resting (reference) layer thicknesses
            coord_iface(1) = 0.0_r8
            do k = 1, maxLevelCell(iCell)
               coord_iface(k + 1) = coord_iface(k) + restingThickness(k, iCell)
            end do
            ! Fill remaining interfaces for unused levels
            do k = maxLevelCell(iCell) + 1, nVertLevels
               coord_iface(k + 1) = coord_iface(k)
            end do

            ! Temperature coarsening (conditional on config)
            if (doTemperature) then
               field_1d(1:nVertLevels) = tracers(index_temperature, 1:nVertLevels, iCell)
               call shr_vcoarsen_avg(field_1d, coord_iface, maxLevelCell(iCell), &
                    depthBounds, nCoarsenLevels, real(fillValue, r8), field_out)
               temperatureCoarsened(1:nCoarsenLevels, iCell) = field_out(1:nCoarsenLevels)
            end if

            ! Salinity coarsening (conditional on config)
            if (doSalinity) then
               field_1d(1:nVertLevels) = tracers(index_salinity, 1:nVertLevels, iCell)
               call shr_vcoarsen_avg(field_1d, coord_iface, maxLevelCell(iCell), &
                    depthBounds, nCoarsenLevels, real(fillValue, r8), field_out)
               salinityCoarsened(1:nCoarsenLevels, iCell) = field_out(1:nCoarsenLevels)
            end if

            ! Zonal velocity coarsening (conditional on config)
            if (doVelocityZonal) then
               field_1d(1:nVertLevels) = velocityZonal(1:nVertLevels, iCell)
               call shr_vcoarsen_avg(field_1d, coord_iface, maxLevelCell(iCell), &
                    depthBounds, nCoarsenLevels, real(fillValue, r8), field_out)
               velocityZonalCoarsened(1:nCoarsenLevels, iCell) = field_out(1:nCoarsenLevels)
            end if

            ! Meridional velocity coarsening (conditional on config)
            if (doVelocityMeridional) then
               field_1d(1:nVertLevels) = velocityMeridional(1:nVertLevels, iCell)
               call shr_vcoarsen_avg(field_1d, coord_iface, maxLevelCell(iCell), &
                    depthBounds, nCoarsenLevels, real(fillValue, r8), field_out)
               velocityMeridionalCoarsened(1:nCoarsenLevels, iCell) = field_out(1:nCoarsenLevels)
            end if

            ! Layer thickness coarsening is always computed (metadata
            ! needed by downstream consumers of coarsened fields)
            do k = 1, nCoarsenLevels
               layerThicknessCoarsened(k, iCell) = 0.0_RKIND
            end do
            depthTop = 0.0_RKIND
            do k = 1, maxLevelCell(iCell)
               call accumulate_thickness(depthTop, restingThickness(k, iCell), &
                    nCoarsenLevels, layerThicknessCoarsened(:, iCell))
               depthTop = depthTop + restingThickness(k, iCell)
            end do
            do k = 1, nCoarsenLevels
               if (layerThicknessCoarsened(k, iCell) == 0.0_RKIND) then
                  layerThicknessCoarsened(k, iCell) = fillValue
               end if
            end do

         end do

         deallocate(coord_iface, field_1d, field_out)

         ! --- Horizontal remapping to lat-lon grid ---
         if (do_horiz_remap) then
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

            allocate(fld_remap(remap_nb_local, 1))
            allocate(field_1lev(1, nCells))

            ! Per-level validity mask: 1=valid ocean cell at this depth, 0=fill
            allocate(levelMask(nCells))

            ! Lazy-compute static per-level wet masks in remap space.
            ! layerThicknessCoarsened was just filled above, so it's
            ! safe to derive each level's source-side mask from it.
            ! The result depends only on bathymetry (static), so one
            ! pass per AM lifetime suffices.
            if (.not. allocated(mask_levels_remap)) then
               call compute_mask_levels_remap()
            end if

            if (time_avg_enabled) then
               ! Accumulate path: remap into accumulators (no immediate write).
               do ifld = 1, nCoarsenFields
                  select case (trim(coarsenFieldNames(ifld)))
                  case ('temperature')
                     coarsenedPtr => temperatureCoarsened
                  case ('salinity')
                     coarsenedPtr => salinityCoarsened
                  case ('velocityZonal')
                     coarsenedPtr => velocityZonalCoarsened
                  case ('velocityMeridional')
                     coarsenedPtr => velocityMeridionalCoarsened
                  case default
                     cycle
                  end select

                  do kOut = 0, nCoarsenLevels - 1
                     call remap_1level(coarsenedPtr, kOut, &
                          (ifld - 1) * nCoarsenLevels + kOut + 1)
                  end do
               end do
               do kOut = 0, nCoarsenLevels - 1
                  call remap_1level(layerThicknessCoarsened, kOut, &
                       nCoarsenFields * nCoarsenLevels + kOut + 1)
               end do
               nAccum = nAccum + 1

               ! Check if output interval elapsed
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
               if (currTime - avg_last_output_time >= avg_output_interval) then
                  call remap_outfile%check_rotate(domain, ierr)
                  if (remap_outfile%file_is_open) then
                     ! time_bnds = [avg window start, avg window end]
                     call remap_outfile%write_time(avg_last_output_time, &
                          currTime, ierr)
                     call write_static_masks_if_needed()

                     do ifld = 1, nCoarsenFields
                        do kOut = 0, nCoarsenLevels - 1
                           write(varName, '(A,A,I0)') &
                                trim(coarsenFieldNames(ifld)), &
                                'Coarsened_', kOut
                           call write_avg_field(trim(varName), &
                                (ifld - 1) * nCoarsenLevels + kOut + 1)
                        end do
                     end do
                     do kOut = 0, nCoarsenLevels - 1
                        write(varName, '(A,I0)') &
                             'layerThicknessCoarsened_', kOut
                        call write_avg_field(trim(varName), &
                             nCoarsenFields * nCoarsenLevels + kOut + 1)
                     end do

                     ! Instantaneous companions at the same trigger moment
                     if (write_inst_companion) then
                        do ifld = 1, nCoarsenFields
                           do kOut = 0, nCoarsenLevels - 1
                              write(varName, '(A,A,I0,A)') &
                                   trim(coarsenFieldNames(ifld)), &
                                   'Coarsened_', kOut, '_inst'
                              call write_inst_companion_field(trim(varName), &
                                   (ifld - 1) * nCoarsenLevels + kOut + 1)
                           end do
                        end do
                        do kOut = 0, nCoarsenLevels - 1
                           write(varName, '(A,I0,A)') &
                                'layerThicknessCoarsened_', kOut, '_inst'
                           call write_inst_companion_field(trim(varName), &
                                nCoarsenFields * nCoarsenLevels + kOut + 1)
                        end do
                     end if
                  end if
                  remap_accum = 0.0_r8
                  remap_accum_count = 0.0_r8
                  nAccum = 0
                  ! Advance by exact interval to avoid drift.
                  avg_last_output_time = avg_last_output_time + avg_output_interval
               end if
            else
               ! Instantaneous: remap each field/level and write under its
               ! own name in the same iteration so fld_remap holds the
               ! correct field at every write call.
               call remap_outfile%check_rotate(domain, ierr)
               if (remap_outfile%file_is_open) then
                  currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
                  ! Instantaneous: time_bnds collapses to [t, t]
                  call remap_outfile%write_time(currTime, currTime, ierr)
                  call write_static_masks_if_needed()

                  do ifld = 1, nCoarsenFields
                     select case (trim(coarsenFieldNames(ifld)))
                     case ('temperature')
                        coarsenedPtr => temperatureCoarsened
                     case ('salinity')
                        coarsenedPtr => salinityCoarsened
                     case ('velocityZonal')
                        coarsenedPtr => velocityZonalCoarsened
                     case ('velocityMeridional')
                        coarsenedPtr => velocityMeridionalCoarsened
                     case default
                        cycle
                     end select

                     do kOut = 0, nCoarsenLevels - 1
                        write(varName, '(A,A,I0)') &
                             trim(coarsenFieldNames(ifld)), &
                             'Coarsened_', kOut
                        call remap_and_write_inst(coarsenedPtr, kOut, trim(varName))
                     end do
                  end do
                  do kOut = 0, nCoarsenLevels - 1
                     write(varName, '(A,I0)') &
                          'layerThicknessCoarsened_', kOut
                     call remap_and_write_inst(layerThicknessCoarsened, kOut, &
                          trim(varName))
                  end do
               end if
            end if

            deallocate(fld_remap, field_1lev, levelMask)
         end if

         block => block % next
      end do

   contains

      subroutine remap_1level(srcArray, kLevel, idx)
         ! Per-cell accumulation that excludes fill values from the running
         ! sum; see AGENTS.md gotcha on time-averaging fill pollution. The
         ! latest_inst snapshot is unmodified (it is overwritten in place
         ! at every call regardless of validity, matching the existing
         ! convention used by the verify dashboard).
         real(kind=RKIND), dimension(:,:), intent(in) :: srcArray
         integer, intent(in) :: kLevel, idx
         integer :: ib
         real(r8), parameter :: fill_thresh = real(SHR_FILL_VALUE, r8) * 0.1_r8
         do iCell = 1, nCells
            if (layerThicknessCoarsened(kLevel + 1, iCell) < fillValue * 0.5_RKIND) then
               levelMask(iCell) = 1
               field_1lev(1, iCell) = srcArray(kLevel + 1, iCell)
            else
               levelMask(iCell) = 0
               field_1lev(1, iCell) = 0.0_RKIND
            end if
         end do
         call ocn_fme_horiz_remap_field(domain, field_1lev, 1, &
              fld_remap, ierr, mask=levelMask)
         if (time_avg_enabled) then
            do ib = 1, remap_nb_local
               if (abs(fld_remap(ib, 1)) < fill_thresh) then
                  remap_accum(ib, idx)       = remap_accum(ib, idx)       + fld_remap(ib, 1)
                  remap_accum_count(ib, idx) = remap_accum_count(ib, idx) + 1.0_r8
               end if
            end do
            if (write_inst_companion) then
               remap_latest_inst(:, idx) = fld_remap(:, 1)
            end if
         end if
      end subroutine remap_1level

      subroutine write_avg_field(vname, idx)
         character(len=*), intent(in) :: vname
         integer, intent(in) :: idx
         integer :: ib
         do ib = 1, remap_nb_local
            if (remap_accum_count(ib, idx) > 0.0_r8) then
               fld_remap(ib, 1) = remap_accum(ib, idx) / remap_accum_count(ib, idx)
            else
               fld_remap(ib, 1) = real(SHR_FILL_VALUE, r8)
            end if
         end do
         call remap_outfile%write_var(vname, fld_remap, 1, ierr)
      end subroutine write_avg_field

      subroutine compute_mask_levels_remap()
         ! Build per-level static wet masks in remap space. The source-side
         ! mask comes from layerThicknessCoarsened (just filled above):
         ! cells where thickness != fillValue at level k are valid sources.
         ! Remap a constant-1 field through apply_masked once per level —
         ! result is 1.0 on target cells with any valid source, fill
         ! elsewhere. Threshold to a hard 0.0/1.0.
         integer :: ic, ib, kLvl
         allocate(mask_levels_remap(remap_nb_local, nCoarsenLevels))
         allocate(mask_write_buf(remap_nb_local, 1))
         do kLvl = 0, nCoarsenLevels - 1
            do ic = 1, nCells
               if (layerThicknessCoarsened(kLvl + 1, ic) < fillValue * 0.5_RKIND) then
                  levelMask(ic) = 1
                  field_1lev(1, ic) = 1.0_RKIND
               else
                  levelMask(ic) = 0
                  field_1lev(1, ic) = 0.0_RKIND
               end if
            end do
            call ocn_fme_horiz_remap_field(domain, field_1lev, 1, &
                 fld_remap, ierr, mask=levelMask)
            do ib = 1, remap_nb_local
               if (abs(fld_remap(ib, 1)) < MASK_FILL_THRESH) then
                  mask_levels_remap(ib, kLvl + 1) = 1.0_r8
               else
                  mask_levels_remap(ib, kLvl + 1) = 0.0_r8
               end if
            end do
         end do
      end subroutine compute_mask_levels_remap

      subroutine write_static_masks_if_needed()
         ! Static masks are written once per file at the first record on
         ! a freshly created file. Append-mode reopens skip this — masks
         ! are already in the file from the original create.
         integer :: kLvl
         character(len=128) :: maskName
         if (remap_outfile%existing_file_reopened) return
         if (remap_outfile%time_record /= 1) return
         if (.not. allocated(mask_levels_remap)) return
         do kLvl = 0, nCoarsenLevels - 1
            mask_write_buf(:, 1) = mask_levels_remap(:, kLvl + 1)
            write(maskName, '(A,I0)') 'mask_', kLvl
            call remap_outfile%write_var(trim(maskName), mask_write_buf, 1, ierr)
         end do
      end subroutine write_static_masks_if_needed

      subroutine remap_and_write_inst(srcArray, kLevel, vname)
         ! Remap one (level of one) field into fld_remap, then write under
         ! vname. Used in instantaneous mode to ensure each write gets the
         ! correct field's data (the previous bulk-remap-then-bulk-write
         ! pattern wrote the last remapped field under every name).
         real(kind=RKIND), dimension(:,:), intent(in) :: srcArray
         integer, intent(in) :: kLevel
         character(len=*), intent(in) :: vname
         do iCell = 1, nCells
            if (layerThicknessCoarsened(kLevel + 1, iCell) < fillValue * 0.5_RKIND) then
               levelMask(iCell) = 1
               field_1lev(1, iCell) = srcArray(kLevel + 1, iCell)
            else
               levelMask(iCell) = 0
               field_1lev(1, iCell) = 0.0_RKIND
            end if
         end do
         call ocn_fme_horiz_remap_field(domain, field_1lev, 1, &
              fld_remap, ierr, mask=levelMask)
         call remap_outfile%write_var(vname, fld_remap, 1, ierr)
      end subroutine remap_and_write_inst

      subroutine write_inst_companion_field(vname, idx)
         ! Snapshot of the latest instantaneous remapped value for slot idx.
         ! Reuses fld_remap as the write buffer so we don't allocate per-call.
         character(len=*), intent(in) :: vname
         integer, intent(in) :: idx
         fld_remap(:, 1) = remap_latest_inst(:, idx)
         call remap_outfile%write_var(vname, fld_remap, 1, ierr)
      end subroutine write_inst_companion_field

   end subroutine ocn_compute_fme_depth_coarsening!}}}

!***********************************************************************

   subroutine accumulate_thickness(depthTop, thickness, nLevs, thickOut)!{{{
      real(kind=RKIND), intent(in) :: depthTop, thickness
      integer, intent(in) :: nLevs
      real(kind=RKIND), intent(inout) :: thickOut(nLevs)

      real(kind=RKIND) :: depthBot, overlapTop, overlapBot, overlap
      integer :: kOut

      depthBot = depthTop + thickness
      do kOut = 1, nLevs
         overlapTop = max(depthTop, real(depthBounds(kOut), RKIND))
         overlapBot = min(depthBot, real(depthBounds(kOut + 1), RKIND))
         overlap = max(0.0_RKIND, overlapBot - overlapTop)
         thickOut(kOut) = thickOut(kOut) + overlap
      end do

   end subroutine accumulate_thickness!}}}

!***********************************************************************

   subroutine ocn_restart_fme_depth_coarsening(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      err = 0

   end subroutine ocn_restart_fme_depth_coarsening!}}}

!***********************************************************************

   subroutine ocn_finalize_fme_depth_coarsening(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err
      integer :: ierr

      err = 0

      call remap_outfile%close_file(ierr)

      if (allocated(depthBounds)) deallocate(depthBounds)
      if (allocated(coarsenFieldNames)) deallocate(coarsenFieldNames)
      if (allocated(remap_accum)) deallocate(remap_accum)
      if (allocated(remap_accum_count)) deallocate(remap_accum_count)
      if (allocated(remap_latest_inst)) deallocate(remap_latest_inst)
      if (allocated(mask_levels_remap)) deallocate(mask_levels_remap)
      if (allocated(mask_write_buf)) deallocate(mask_write_buf)

   end subroutine ocn_finalize_fme_depth_coarsening!}}}

!***********************************************************************

   function coarsen_field_units(fieldName) result(units_out)
   !> Map a depth-coarsen source field name to the unit string from the
   !> MPAS-O Registry. Depth coarsening is a thickness-weighted mean and
   !> therefore unit-preserving, so the output unit equals the source.
   !> Falls back to '1' (CF "dimensionless") with a one-line WARN if the
   !> field name isn't recognised — keeps the file CF-loadable rather
   !> than emitting empty unit strings.
      character(len=*), intent(in) :: fieldName
      character(len=32) :: units_out

      select case (trim(fieldName))
      case ('temperature')
         units_out = 'degC'
      case ('salinity')
         units_out = 'PSU'
      case ('velocityZonal', 'velocityMeridional')
         units_out = 'm s-1'
      case default
         units_out = '1'
         call mpas_log_write( &
              'fmeDepthCoarsening: WARN no units mapping for source ' // &
              'field "' // trim(fieldName) // '"; tagging units="1"', &
              messageType=MPAS_LOG_WARN)
      end select
   end function coarsen_field_units

end module ocn_fme_depth_coarsening

! vim: foldmethod=marker
