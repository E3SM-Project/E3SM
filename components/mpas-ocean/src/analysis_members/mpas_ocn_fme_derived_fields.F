! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_fme_derived_fields
!
!> \brief MPAS ocean analysis member: FME derived surface fields
!> \author FME Training Data Support
!> \date   April 2026
!> \details
!>  Computes derived surface fields for Full Model Emulation training data:
!>  - SST: sea surface temperature from top-level temperature tracer
!>  - SSS: sea surface salinity from top-level salinity tracer
!>  - surfaceHeatFluxTotal: combined SW + LW_down + latent + sensible
!>
!>  When config_AM_fmeDerivedFields_horiz_remap_file is set to a valid
!>  map file path, the derived fields are also remapped to a lat-lon grid
!>  and written to a separate NetCDF file (fmeDerivedFieldsRemapped.*.nc).
!>
!-----------------------------------------------------------------------

module ocn_fme_derived_fields

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_kind_types
   use mpas_log, only: mpas_log_write

   use ocn_constants
   use ocn_diagnostics_variables

   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_horiz_remap_mod, only: SHR_FILL_VALUE
   use ocn_fme_horiz_remap, only: ocn_fme_horiz_remap_init, &
        ocn_fme_horiz_remap_init_time_epoch, &
        ocn_fme_horiz_remap_field, ocn_fme_horiz_remap_active, &
        ocn_fme_horiz_remap_get_dims, ocn_fme_remap_file_t, &
        ocn_fme_remap_build_filename, &
        ocn_fme_horiz_remap_write_accum, ocn_fme_horiz_remap_read_accum

   implicit none
   private
   save

   public :: ocn_init_fme_derived_fields, &
             ocn_compute_fme_derived_fields, &
             ocn_restart_fme_derived_fields, &
             ocn_finalize_fme_derived_fields, &
             ocn_fme_derived_fields_write_restart

   ! Remapped output workspace
   logical :: do_horiz_remap = .false.
   integer :: remap_nlat = 0, remap_nlon = 0, remap_nb_local = 0
   type(ocn_fme_remap_file_t) :: remap_outfile

   ! Time-averaging state for remapped output
   logical :: time_avg_enabled = .false.
   integer :: nAccum = 0
   type(MPAS_TimeInterval_type) :: avg_output_interval
   type(MPAS_Time_type) :: avg_last_output_time
   ! Slot assignments for remap_accum (must stay in sync with init/compute):
   !   1  sst                          13  rainFlux
   !   2  sss                          14  snowFlux
   !   3  surfaceHeatFluxTotal         15  evaporationFlux
   !   4  ssh                          16  seaIceHeatFlux
   !   5  shortWaveHeatFlux            17  icebergHeatFlux
   !   6  longWaveHeatFluxDown         18  seaIceFreshWaterFlux
   !   7  latentHeatFlux               19  seaIceSalinityFlux
   !   8  sensibleHeatFlux             20  riverRunoffFlux
   !   9  windStressZonal              21  iceRunoffFlux
   !  10  windStressMeridional         22  snowFluxLf
   !  11  longWaveHeatFluxUp           23  iceRunoffFluxLf
   !  12  (unused, was totalFreshWaterTemperatureFlux until 2026-04-30,
   !       removed per SamudrACE spec p3ai/6154289880 -- struck through)
   !                                   24  riverRunoffFluxTemperature
   integer, parameter :: N_AVG_FIELDS = 24
   real(r8), allocatable :: remap_accum(:,:)        ! (n_b_local, N_AVG_FIELDS)
   ! Static surface ocean mask in remap space, written as `mask_2d` in
   ! every output file. 1.0 on ocean cells, 0.0 elsewhere. Computed once
   ! lazily in compute by remapping a constant-1 native field through
   ! apply_masked with the static maxLevelCell mask. The result is the
   ! canonical mask the SamudrACE pipeline expects (see compute_ocean
   ! _dataset_e3sm.py:53). Shape (n_b_local, 1) so the existing
   ! write_var(:, :) interface accepts it without reshape.
   real(r8), allocatable :: mask_2d_remap(:,:)
   ! Threshold for detecting fill in remapped output (10% of SHR_FILL_VALUE,
   ! same convention as ocn_fme_horiz_remap_field's fill_detect).
   real(r8), parameter :: MASK_FILL_THRESH = real(SHR_FILL_VALUE, r8) * 0.1_r8

   ! Per-cell, per-field valid-sample counter parallel to remap_accum.
   ! Ocean's land mask is static (maxLevelCell does not change at run
   ! time), so any given target cell is fill on every timestep or on
   ! none, and the legacy global-nAccum division would have averaged
   ! correctly here. The per-cell count is added defensively so that any
   ! future field whose validity flips in time (e.g., ice-shelf cavities,
   ! seasonal sea-ice influence) is handled the same way as MPAS-SI. See
   ! AGENTS.md gotcha on time-averaging fill pollution.
   real(r8), allocatable :: remap_accum_count(:,:)  ! (n_b_local, N_AVG_FIELDS)

!***********************************************************************

contains

!***********************************************************************

   subroutine ocn_init_fme_derived_fields(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND), pointer :: config_AM_fmeDerivedFields_horiz_remap_file
      character(len=StrKIND), pointer :: config_stream_output_interval
      logical, pointer :: config_time_averaging
      character(len=16) :: cm   ! CF cell_methods choice -- 'time: mean' / 'time: point'
      integer :: ierr

      err = 0

      ! Check if horizontal remapping is requested
      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeDerivedFields_horiz_remap_file', &
           config_AM_fmeDerivedFields_horiz_remap_file)

      if (associated(config_AM_fmeDerivedFields_horiz_remap_file)) then
         if (trim(config_AM_fmeDerivedFields_horiz_remap_file) /= 'none' .and. &
             len_trim(config_AM_fmeDerivedFields_horiz_remap_file) > 0) then

            call ocn_fme_horiz_remap_init(domain, &
                 trim(config_AM_fmeDerivedFields_horiz_remap_file), ierr)
            if (ierr /= 0) then
               call mpas_log_write( &
                    'fmeDerivedFields: horiz remap init failed, continuing without remap', &
                    messageType=MPAS_LOG_WARN)
               do_horiz_remap = .false.
            else
               ! Auto-discover the CF time epoch -- see fmeDepthCoarsening
               ! for the design. Idempotent across all four FME AMs.
               call ocn_fme_horiz_remap_init_time_epoch(domain, &
                    'fmeDerivedFields', ierr)
               if (ierr /= 0) then
                  call mpas_log_write('fmeDerivedFields: time epoch ' // &
                       'setup failed, disabling horiz remap', &
                       messageType=MPAS_LOG_WARN)
                  do_horiz_remap = .false.
                  err = ierr
                  return
               end if
               do_horiz_remap = ocn_fme_horiz_remap_active()
               if (do_horiz_remap) then
                  call ocn_fme_horiz_remap_get_dims(remap_nlat, remap_nlon, remap_nb_local)
                  call mpas_log_write('fmeDerivedFields: horiz remap active, ' &
                       // 'target grid $i x $i', &
                       intArgs=(/remap_nlon, remap_nlat/))
                  remap_outfile%stream_name = 'fmeDerivedFieldsOutput'
                  ! Pick CF cell_methods up-front based on the runtime
                  ! averaging mode. Production: time_averaging=.true.
                  ! ('time: mean'); inst-mode pathway tags 'time: point'.
                  call mpas_pool_get_config(domain % configs, &
                       'config_AM_fmeDerivedFields_time_averaging', &
                       config_time_averaging)
                  cm = 'time: point'
                  if (associated(config_time_averaging)) then
                     if (config_time_averaging) cm = 'time: mean'
                  end if
                  call remap_outfile%register_var('sst', 'degC', &
                       'Sea surface temperature', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('sss', 'PSU', &
                       'Sea surface salinity', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('surfaceHeatFluxTotal', 'W m^-2', &
                       'Total surface heat flux', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('ssh', 'm', &
                       'Sea surface height', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('shortWaveHeatFlux', 'W m^-2', &
                       'Shortwave heat flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('longWaveHeatFluxDown', 'W m^-2', &
                       'Downward longwave heat flux', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('latentHeatFlux', 'W m^-2', &
                       'Latent heat flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('sensibleHeatFlux', 'W m^-2', &
                       'Sensible heat flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('windStressZonal', 'N m^-2', &
                       'Zonal wind stress on ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('windStressMeridional', 'N m^-2', &
                       'Meridional wind stress on ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('longWaveHeatFluxUp', 'W m^-2', &
                       'Upward longwave heat flux from ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('rainFlux', 'kg m^-2 s^-1', &
                       'Rain mass flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('snowFlux', 'kg m^-2 s^-1', &
                       'Snow mass flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('evaporationFlux', 'kg m^-2 s^-1', &
                       'Evaporation mass flux from ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('seaIceHeatFlux', 'W m^-2', &
                       'Sea-ice heat flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('icebergHeatFlux', 'W m^-2', &
                       'Iceberg heat flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('seaIceFreshWaterFlux', 'kg m^-2 s^-1', &
                       'Sea-ice fresh water mass flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('seaIceSalinityFlux', 'kg m^-2 s^-1', &
                       'Sea-ice salt flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('riverRunoffFlux', 'kg m^-2 s^-1', &
                       'River runoff mass flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('iceRunoffFlux', 'kg m^-2 s^-1', &
                       'Land-ice runoff mass flux into ocean', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('snowFluxLf', 'W m^-2', &
                       'snowFlux times latent heat of fusion', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('iceRunoffFluxLf', 'W m^-2', &
                       'iceRunoffFlux times latent heat of fusion', ierr, cell_methods=trim(cm))
                  call remap_outfile%register_var('riverRunoffFluxTemperature', &
                       'kg degC m^-2 s^-1', &
                       'riverRunoffFlux times top-layer temperature', ierr, &
                       cell_methods=trim(cm))
                  ! Static surface ocean mask: shape (lon, lat), 1=ocean,
                  ! 0=land. Written once per file. Matches the SamudrACE
                  ! convention (compute_ocean_dataset_e3sm.py:53).
                  call remap_outfile%register_var('mask_2d', '1', &
                       'Surface ocean wet mask (1=ocean, 0=land/invalid)', ierr, &
                       is_static=.true.)

                  ! Time-averaging setup
                  call mpas_pool_get_config(domain % configs, &
                       'config_AM_fmeDerivedFields_time_averaging', config_time_averaging)
                  if (associated(config_time_averaging)) then
                     if (config_time_averaging) then
                        call mpas_pool_get_config(domain % configs, &
                             'config_AM_fmeDerivedFields_stream_output_interval', &
                             config_stream_output_interval)
                        call mpas_set_timeInterval(avg_output_interval, &
                             timeString=trim(config_stream_output_interval), ierr=ierr)
                        ! Default baseline is MPAS_NOW so a warm restart
                        ! without a sidecar uses the restart-time as the
                        ! averaging baseline. Persisted state below
                        ! overrides this when a sidecar is found.
                        avg_last_output_time = mpas_get_clock_time( &
                             domain % clock, MPAS_NOW, ierr)
                        if (.not. allocated(remap_accum)) then
                           allocate(remap_accum(remap_nb_local, N_AVG_FIELDS))
                        end if
                        if (.not. allocated(remap_accum_count)) then
                           allocate(remap_accum_count(remap_nb_local, N_AVG_FIELDS))
                        end if
                        remap_accum = 0.0_r8
                        remap_accum_count = 0.0_r8
                        nAccum = 0
                        time_avg_enabled = .true.
                        call mpas_log_write('fmeDerivedFields: time averaging enabled')
                        call ocn_fme_derived_fields_read_restart(domain, ierr)
                     end if
                  end if
               end if
            end if
         end if
      end if

   end subroutine ocn_init_fme_derived_fields!}}}

!***********************************************************************

   subroutine ocn_fme_derived_fields_read_restart(domain, err)!{{{
   !> Restore accumulator state on warm restart (silent no-op otherwise).
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

      call ocn_fme_horiz_remap_read_accum('fmeDerivedFields', remap_accum, &
           N_AVG_FIELDS, nAccum, ts_str, &
           count_arr=remap_accum_count, success=ok)

      if (ok .and. len_trim(ts_str) > 0) then
         call mpas_set_time(avg_last_output_time, &
              dateTimeString=trim(ts_str), ierr=ierr_local)
         if (ierr_local /= 0) then
            call mpas_log_write('fmeDerivedFields: failed to parse ' // &
                 'persisted avg_last_output_time="' // trim(ts_str) // &
                 '"; falling back to MPAS_NOW baseline.', &
                 messageType=MPAS_LOG_WARN)
            avg_last_output_time = mpas_get_clock_time( &
                 domain % clock, MPAS_NOW, ierr_local)
         end if
      end if

   end subroutine ocn_fme_derived_fields_read_restart!}}}

!***********************************************************************

   subroutine ocn_fme_derived_fields_write_restart(domain, err)!{{{
   !> Persist accumulator state to the AM's sidecar file.
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND) :: ts_str
      integer :: ierr_local

      err = 0
      if (.not. time_avg_enabled) return
      if (.not. allocated(remap_accum)) return

      call mpas_get_time(avg_last_output_time, &
           dateTimeString=ts_str, ierr=ierr_local)

      call ocn_fme_horiz_remap_write_accum('fmeDerivedFields', remap_accum, &
           N_AVG_FIELDS, nAccum, trim(ts_str), &
           count_arr=remap_accum_count, err=err)

   end subroutine ocn_fme_derived_fields_write_restart!}}}

!***********************************************************************

   subroutine ocn_compute_fme_derived_fields(domain, timeLevel, err)!{{{

      integer, intent(in) :: timeLevel
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: tracersPool
      type(mpas_pool_type), pointer :: forcingPool
      type(mpas_pool_type), pointer :: diagnosticsPool
      type(mpas_pool_type), pointer :: amPool

      integer, pointer :: nCellsSolve, nCells
      integer, pointer :: index_temperature, index_salinity
      integer, dimension(:), pointer :: maxLevelCell
      integer :: iCell, ierr

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers
      real(kind=RKIND), dimension(:), pointer :: sst, sss, surfaceHeatFluxTotal
      real(kind=RKIND), dimension(:), pointer :: snowFluxLf, iceRunoffFluxLf
      real(kind=RKIND), dimension(:), pointer :: riverRunoffFluxTemperature
      real(kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux
      real(kind=RKIND), dimension(:), pointer :: longWaveHeatFluxDown
      real(kind=RKIND), dimension(:), pointer :: longWaveHeatFluxUp
      real(kind=RKIND), dimension(:), pointer :: latentHeatFlux
      real(kind=RKIND), dimension(:), pointer :: sensibleHeatFlux
      real(kind=RKIND), dimension(:), pointer :: rainFluxPtr, snowFluxPtr
      real(kind=RKIND), dimension(:), pointer :: evaporationFluxPtr
      real(kind=RKIND), dimension(:), pointer :: seaIceHeatFluxPtr, icebergHeatFluxPtr
      real(kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFluxPtr
      real(kind=RKIND), dimension(:), pointer :: seaIceSalinityFluxPtr
      real(kind=RKIND), dimension(:), pointer :: riverRunoffFluxPtr, iceRunoffFluxPtr
      real(kind=RKIND), dimension(:), pointer :: ssh_ptr
      real(kind=RKIND), dimension(:), pointer :: windStressZonalPtr
      real(kind=RKIND), dimension(:), pointer :: windStressMeridionalPtr
      ! Fill value for land-masked cells
      real(kind=RKIND), parameter :: fillValue = real(SHR_FILL_VALUE, RKIND)
      ! Workspace for masked copies of forcing fields
      real(kind=RKIND), allocatable :: fld_masked(:)
      ! Remap workspace
      real(r8), allocatable :: fld_remap(:,:)
      real(kind=RKIND), allocatable :: field_1lev(:,:)
      type(MPAS_Time_Type) :: currTime

      err = 0

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'fmeDerivedFieldsAM', amPool)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(amPool, 'sst', sst)
         call mpas_pool_get_array(amPool, 'sss', sss)
         call mpas_pool_get_array(amPool, 'surfaceHeatFluxTotal', surfaceHeatFluxTotal)
         call mpas_pool_get_array(amPool, 'snowFluxLf', snowFluxLf)
         call mpas_pool_get_array(amPool, 'iceRunoffFluxLf', iceRunoffFluxLf)
         call mpas_pool_get_array(amPool, 'riverRunoffFluxTemperature', &
              riverRunoffFluxTemperature)

         call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
         call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
         call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
         call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
         call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
         call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonalPtr)
         call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridionalPtr)
         call mpas_pool_get_array(forcingPool, 'rainFlux', rainFluxPtr)
         call mpas_pool_get_array(forcingPool, 'snowFlux', snowFluxPtr)
         call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFluxPtr)
         call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFluxPtr)
         call mpas_pool_get_array(forcingPool, 'icebergHeatFlux', icebergHeatFluxPtr)
         call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', &
              seaIceFreshWaterFluxPtr)
         call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxPtr)
         call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFluxPtr)
         call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFluxPtr)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         do iCell = 1, nCellsSolve
            ! Land cells (maxLevelCell==0) get fillValue so the remap
            ! correctly excludes them from weighted averages.
            if (maxLevelCell(iCell) == 0) then
               sst(iCell) = fillValue
               sss(iCell) = fillValue
               surfaceHeatFluxTotal(iCell) = fillValue
               snowFluxLf(iCell) = fillValue
               iceRunoffFluxLf(iCell) = fillValue
               riverRunoffFluxTemperature(iCell) = fillValue
               cycle
            end if

            sst(iCell) = tracers(index_temperature, 1, iCell)
            sss(iCell) = tracers(index_salinity, 1, iCell)

            if (associated(shortWaveHeatFlux) .and. &
                associated(longWaveHeatFluxDown) .and. &
                associated(latentHeatFlux) .and. &
                associated(sensibleHeatFlux)) then
               ! Sum of coupler fluxes (raw, ice-fraction-weighted values
               ! from the coupler, matching timeSeriesStatsCustom behavior).
               surfaceHeatFluxTotal(iCell) = shortWaveHeatFlux(iCell) &
                    + longWaveHeatFluxDown(iCell) &
                    + latentHeatFlux(iCell) &
                    + sensibleHeatFlux(iCell)
            else
               surfaceHeatFluxTotal(iCell) = fillValue
            end if

            ! Derived heat-equivalent and temperature-flux products.
            ! latent_heat_fusion_mks is J/kg; mass-flux * J/kg = W/m^2.
            if (associated(snowFluxPtr)) then
               snowFluxLf(iCell) = snowFluxPtr(iCell) * latent_heat_fusion_mks
            else
               snowFluxLf(iCell) = fillValue
            end if
            if (associated(iceRunoffFluxPtr)) then
               iceRunoffFluxLf(iCell) = iceRunoffFluxPtr(iCell) * latent_heat_fusion_mks
            else
               iceRunoffFluxLf(iCell) = fillValue
            end if
            if (associated(riverRunoffFluxPtr)) then
               riverRunoffFluxTemperature(iCell) = riverRunoffFluxPtr(iCell) &
                    * tracers(index_temperature, 1, iCell)
            else
               riverRunoffFluxTemperature(iCell) = fillValue
            end if
         end do

         ! --- Horizontal remapping to lat-lon grid ---
         if (do_horiz_remap) then
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            allocate(fld_remap(remap_nb_local, 1))
            allocate(field_1lev(1, nCells))
            allocate(fld_masked(nCells))

            ! Lazy-compute the static surface ocean mask in remap space.
            ! Done once per AM lifetime — the source-side mask
            ! (maxLevelCell > 0) is static, so the remapped result is too.
            if (.not. allocated(mask_2d_remap)) then
               call compute_mask_2d_remap()
            end if

            ! Remap each field to lat-lon
            ! Use pressureAdjustedSSH (from diagnostics) instead of raw ssh
            ! (from state) to avoid extreme values under ice shelves.
            call mpas_pool_get_array(diagnosticsPool, 'pressureAdjustedSSH', ssh_ptr)

            if (time_avg_enabled) then
               ! Accumulate path: each remap adds to remap_accum; write when interval elapses.
               call remap_1field(sst, 1)
               call remap_1field(sss, 2)
               call remap_1field(surfaceHeatFluxTotal, 3)

               ! SSH: use pressure-adjusted version for physically meaningful values
               if (associated(ssh_ptr)) call remap_1field(ssh_ptr, 4)

               ! Coupler fluxes: land-mask and remap (raw coupler values,
               ! matching timeSeriesStatsCustom behavior).
               if (associated(shortWaveHeatFlux)) &
                    call mask_and_remap(shortWaveHeatFlux, 5)
               if (associated(longWaveHeatFluxDown)) &
                    call mask_and_remap(longWaveHeatFluxDown, 6)
               if (associated(latentHeatFlux)) &
                    call mask_and_remap(latentHeatFlux, 7)
               if (associated(sensibleHeatFlux)) &
                    call mask_and_remap(sensibleHeatFlux, 8)
               if (associated(windStressZonalPtr)) &
                    call mask_and_remap(windStressZonalPtr, 9)
               if (associated(windStressMeridionalPtr)) &
                    call mask_and_remap(windStressMeridionalPtr, 10)
               if (associated(longWaveHeatFluxUp)) &
                    call mask_and_remap(longWaveHeatFluxUp, 11)
               ! idx=12 unused (was totalFreshWaterTemperatureFlux, removed 2026-04-30)
               if (associated(rainFluxPtr)) &
                    call mask_and_remap(rainFluxPtr, 13)
               if (associated(snowFluxPtr)) &
                    call mask_and_remap(snowFluxPtr, 14)
               if (associated(evaporationFluxPtr)) &
                    call mask_and_remap(evaporationFluxPtr, 15)
               if (associated(seaIceHeatFluxPtr)) &
                    call mask_and_remap(seaIceHeatFluxPtr, 16)
               if (associated(icebergHeatFluxPtr)) &
                    call mask_and_remap(icebergHeatFluxPtr, 17)
               if (associated(seaIceFreshWaterFluxPtr)) &
                    call mask_and_remap(seaIceFreshWaterFluxPtr, 18)
               if (associated(seaIceSalinityFluxPtr)) &
                    call mask_and_remap(seaIceSalinityFluxPtr, 19)
               if (associated(riverRunoffFluxPtr)) &
                    call mask_and_remap(riverRunoffFluxPtr, 20)
               if (associated(iceRunoffFluxPtr)) &
                    call mask_and_remap(iceRunoffFluxPtr, 21)
               ! Derived AM products are already land-masked in the cell loop,
               ! so use remap_1field directly (no second masking pass).
               call remap_1field(snowFluxLf, 22)
               call remap_1field(iceRunoffFluxLf, 23)
               call remap_1field(riverRunoffFluxTemperature, 24)
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
                     call write_avg_field('sst', 1)
                     call write_avg_field('sss', 2)
                     call write_avg_field('surfaceHeatFluxTotal', 3)
                     call write_avg_field('ssh', 4)
                     call write_avg_field('shortWaveHeatFlux', 5)
                     call write_avg_field('longWaveHeatFluxDown', 6)
                     call write_avg_field('latentHeatFlux', 7)
                     call write_avg_field('sensibleHeatFlux', 8)
                     call write_avg_field('windStressZonal', 9)
                     call write_avg_field('windStressMeridional', 10)
                     call write_avg_field('longWaveHeatFluxUp', 11)
                     ! idx=12 unused (was totalFreshWaterTemperatureFlux)
                     call write_avg_field('rainFlux', 13)
                     call write_avg_field('snowFlux', 14)
                     call write_avg_field('evaporationFlux', 15)
                     call write_avg_field('seaIceHeatFlux', 16)
                     call write_avg_field('icebergHeatFlux', 17)
                     call write_avg_field('seaIceFreshWaterFlux', 18)
                     call write_avg_field('seaIceSalinityFlux', 19)
                     call write_avg_field('riverRunoffFlux', 20)
                     call write_avg_field('iceRunoffFlux', 21)
                     call write_avg_field('snowFluxLf', 22)
                     call write_avg_field('iceRunoffFluxLf', 23)
                     call write_avg_field('riverRunoffFluxTemperature', 24)
                  end if
                  remap_accum = 0.0_r8
                  remap_accum_count = 0.0_r8
                  nAccum = 0
                  ! Advance by exact interval to avoid drift.
                  avg_last_output_time = avg_last_output_time + avg_output_interval
               end if
            else
               ! Instantaneous: remap and write each field one at a time so
               ! fld_remap holds the correct field at every write call.
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
               call remap_outfile%check_rotate(domain, ierr)
               if (remap_outfile%file_is_open) then
                  ! Instantaneous: time_bnds collapses to [t, t]
                  call remap_outfile%write_time(currTime, currTime, ierr)
                  call write_static_masks_if_needed()
                  call inst_remap_and_write(sst, 'sst', .false.)
                  call inst_remap_and_write(sss, 'sss', .false.)
                  call inst_remap_and_write(surfaceHeatFluxTotal, 'surfaceHeatFluxTotal', .false.)
                  if (associated(ssh_ptr)) &
                       call inst_remap_and_write(ssh_ptr, 'ssh', .false.)
                  if (associated(shortWaveHeatFlux)) &
                       call inst_remap_and_write(shortWaveHeatFlux, 'shortWaveHeatFlux', .true.)
                  if (associated(longWaveHeatFluxDown)) &
                       call inst_remap_and_write(longWaveHeatFluxDown, 'longWaveHeatFluxDown', .true.)
                  if (associated(latentHeatFlux)) &
                       call inst_remap_and_write(latentHeatFlux, 'latentHeatFlux', .true.)
                  if (associated(sensibleHeatFlux)) &
                       call inst_remap_and_write(sensibleHeatFlux, 'sensibleHeatFlux', .true.)
                  if (associated(windStressZonalPtr)) &
                       call inst_remap_and_write(windStressZonalPtr, 'windStressZonal', .true.)
                  if (associated(windStressMeridionalPtr)) &
                       call inst_remap_and_write(windStressMeridionalPtr, 'windStressMeridional', .true.)
                  if (associated(longWaveHeatFluxUp)) &
                       call inst_remap_and_write(longWaveHeatFluxUp, 'longWaveHeatFluxUp', .true.)
                  if (associated(rainFluxPtr)) &
                       call inst_remap_and_write(rainFluxPtr, 'rainFlux', .true.)
                  if (associated(snowFluxPtr)) &
                       call inst_remap_and_write(snowFluxPtr, 'snowFlux', .true.)
                  if (associated(evaporationFluxPtr)) &
                       call inst_remap_and_write(evaporationFluxPtr, 'evaporationFlux', .true.)
                  if (associated(seaIceHeatFluxPtr)) &
                       call inst_remap_and_write(seaIceHeatFluxPtr, 'seaIceHeatFlux', .true.)
                  if (associated(icebergHeatFluxPtr)) &
                       call inst_remap_and_write(icebergHeatFluxPtr, 'icebergHeatFlux', .true.)
                  if (associated(seaIceFreshWaterFluxPtr)) &
                       call inst_remap_and_write(seaIceFreshWaterFluxPtr, &
                            'seaIceFreshWaterFlux', .true.)
                  if (associated(seaIceSalinityFluxPtr)) &
                       call inst_remap_and_write(seaIceSalinityFluxPtr, &
                            'seaIceSalinityFlux', .true.)
                  if (associated(riverRunoffFluxPtr)) &
                       call inst_remap_and_write(riverRunoffFluxPtr, 'riverRunoffFlux', .true.)
                  if (associated(iceRunoffFluxPtr)) &
                       call inst_remap_and_write(iceRunoffFluxPtr, 'iceRunoffFlux', .true.)
                  call inst_remap_and_write(snowFluxLf, 'snowFluxLf', .false.)
                  call inst_remap_and_write(iceRunoffFluxLf, 'iceRunoffFluxLf', .false.)
                  call inst_remap_and_write(riverRunoffFluxTemperature, &
                       'riverRunoffFluxTemperature', .false.)
               end if
            end if

            deallocate(fld_remap, field_1lev, fld_masked)
         end if

         block => block % next
      end do

   contains

      subroutine remap_1field(native_fld, idx)
         ! Per-cell accumulation that excludes fill values from the running
         ! sum. Defensive: ocean's land mask is static, but any field whose
         ! validity flips in time would otherwise pollute the average with
         ! fractional-of-fill values. See AGENTS.md gotcha on time-averaging
         ! fill pollution.
         real(kind=RKIND), dimension(:), intent(in) :: native_fld
         integer, intent(in) :: idx
         integer :: ib
         real(r8), parameter :: fill_thresh = real(SHR_FILL_VALUE, r8) * 0.1_r8
         field_1lev(1, 1:nCells) = native_fld(1:nCells)
         call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
              mask=maxLevelCell)
         if (time_avg_enabled) then
            do ib = 1, remap_nb_local
               if (abs(fld_remap(ib, 1)) < fill_thresh) then
                  remap_accum(ib, idx)       = remap_accum(ib, idx)       + fld_remap(ib, 1)
                  remap_accum_count(ib, idx) = remap_accum_count(ib, idx) + 1.0_r8
               end if
            end do
         end if
      end subroutine remap_1field

      subroutine mask_and_remap(native_fld, idx)
         ! Copy field with land masking, then remap.
         ! Land cells (maxLevelCell==0) get fillValue so the remap
         ! correctly excludes them from weighted averages.
         ! No ice-fraction correction is applied — raw coupler values
         ! are preserved, matching timeSeriesStatsCustom behavior.
         real(kind=RKIND), dimension(:), intent(in) :: native_fld
         integer, intent(in) :: idx
         integer :: ic
         fld_masked(1:nCells) = native_fld(1:nCells)
         do ic = 1, nCells
            if (maxLevelCell(ic) == 0) then
               fld_masked(ic) = fillValue
            end if
         end do
         call remap_1field(fld_masked, idx)
      end subroutine mask_and_remap

      subroutine write_avg_field(varname, idx)
         ! Divide by per-cell count rather than the global nAccum so cells
         ! that flip between valid and fill within the window are averaged
         ! over their valid samples only. Cells with count == 0 emit fill.
         character(len=*), intent(in) :: varname
         integer, intent(in) :: idx
         integer :: ib
         do ib = 1, remap_nb_local
            if (remap_accum_count(ib, idx) > 0.0_r8) then
               fld_remap(ib, 1) = remap_accum(ib, idx) / remap_accum_count(ib, idx)
            else
               fld_remap(ib, 1) = real(SHR_FILL_VALUE, r8)
            end if
         end do
         call remap_outfile%write_var(varname, fld_remap, 1, ierr)
      end subroutine write_avg_field

      subroutine compute_mask_2d_remap()
         ! Build the static surface ocean mask in remap space by remapping
         ! a constant-1 ocean field through apply_masked. Land cells
         ! (maxLevelCell == 0) have mask=0 in the remap, so target cells
         ! whose source-side coverage is entirely land come out as
         ! SHR_FILL_VALUE; everything else comes out as 1.0 (a weighted
         ! average of constant 1s = 1). Threshold to a hard 0.0/1.0.
         integer :: ic, ib
         allocate(mask_2d_remap(remap_nb_local, 1))
         do ic = 1, nCells
            if (maxLevelCell(ic) > 0) then
               field_1lev(1, ic) = 1.0_RKIND
            else
               field_1lev(1, ic) = 0.0_RKIND
            end if
         end do
         call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
              mask=maxLevelCell)
         do ib = 1, remap_nb_local
            if (abs(fld_remap(ib, 1)) < MASK_FILL_THRESH) then
               mask_2d_remap(ib, 1) = 1.0_r8
            else
               mask_2d_remap(ib, 1) = 0.0_r8
            end if
         end do
      end subroutine compute_mask_2d_remap

      subroutine write_static_masks_if_needed()
         ! Static masks are written once per file at the first record on
         ! a freshly created file. Append-mode reopens (warm restart) skip
         ! this — the masks are already in the file from leg 1.
         if (remap_outfile%existing_file_reopened) return
         if (remap_outfile%time_record /= 1) return
         if (.not. allocated(mask_2d_remap)) return
         call remap_outfile%write_var('mask_2d', mask_2d_remap, 1, ierr)
      end subroutine write_static_masks_if_needed

      subroutine inst_remap_and_write(native_fld, varname, do_mask)
         ! Remap one field into fld_remap, then write under varname.
         ! Used in instantaneous mode to ensure each varname gets the
         ! correct field's data (the previous bulk-remap-then-bulk-write
         ! pattern wrote the last remapped field under every name).
         real(kind=RKIND), dimension(:), intent(in) :: native_fld
         character(len=*), intent(in) :: varname
         logical, intent(in) :: do_mask
         if (do_mask) then
            fld_masked(1:nCells) = native_fld(1:nCells)
            do iCell = 1, nCells
               if (maxLevelCell(iCell) == 0) then
                  fld_masked(iCell) = fillValue
               end if
            end do
            field_1lev(1, 1:nCells) = fld_masked(1:nCells)
         else
            field_1lev(1, 1:nCells) = native_fld(1:nCells)
         end if
         call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
              mask=maxLevelCell)
         call remap_outfile%write_var(varname, fld_remap, 1, ierr)
      end subroutine inst_remap_and_write

   end subroutine ocn_compute_fme_derived_fields!}}}

!***********************************************************************

   subroutine ocn_restart_fme_derived_fields(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      err = 0

   end subroutine ocn_restart_fme_derived_fields!}}}

!***********************************************************************

   subroutine ocn_finalize_fme_derived_fields(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err
      integer :: ierr

      err = 0

      call remap_outfile%close_file(ierr)
      if (allocated(mask_2d_remap)) deallocate(mask_2d_remap)

   end subroutine ocn_finalize_fme_derived_fields!}}}

end module ocn_fme_derived_fields

! vim: foldmethod=marker
