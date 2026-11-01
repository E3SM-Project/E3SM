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
!>  - surfaceHeatFluxTotal: combined SW + LW_down + latent + sensible ONLY
!>  - surfaceHeatFluxNet: net surface heat flux (conservationCheck heat
!>    terms incl. LWup, sea ice, iceberg, -Lf*snow, -Lf*ice runoff, frazil)
!>  - frazilHeatFlux, massFluxHeatContent and the remaining freshwater
!>    budget terms (0301 output spec)
!>  - oceanWetFraction: static 1-deg open-ocean fraction (remapped file only)
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
        ocn_fme_horiz_remap_get_dims, ocn_fme_horiz_remap_get_coverage, &
        ocn_fme_remap_file_t, &
        ocn_fme_remap_build_filename, &
        ocn_fme_horiz_remap_write_accum, ocn_fme_horiz_remap_read_accum, &
        ocn_fme_horiz_remap_set_rewind

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
   !  25  surfaceHeatFluxNet           29  frazilIceFreshwaterFlux
   !  26  frazilHeatFlux               30  removedRiverRunoffFlux
   !  27  massFluxHeatContent          31  removedIceRunoffFlux
   !  28  icebergFreshWaterFlux
   ! Slots 25-31 added for E3SMv3 0301 (2026-11). Changing N_AVG_FIELDS
   ! changes the accumulator sidecar layout; read_restart aborts on a
   ! sidecar with a different slot count instead of silently restarting
   ! the averaging window (see ocn_fme_derived_fields_read_restart).
   integer, parameter :: N_AVG_FIELDS = 31
   real(r8), allocatable :: remap_accum(:,:)        ! (n_b_local, N_AVG_FIELDS)
   ! Static surface ocean mask in remap space, written as `mask_2d` in
   ! every output file. 1.0 on ocean cells, 0.0 elsewhere. Computed once
   ! lazily in compute by remapping a constant-1 native field through
   ! apply_masked with the static maxLevelCell mask. The result is the
   ! canonical mask the SamudrACE pipeline expects (see compute_ocean
   ! _dataset_e3sm.py:53). Shape (n_b_local, 1) so the existing
   ! write_var(:, :) interface accepts it without reshape.
   real(r8), allocatable :: mask_2d_remap(:,:)
   ! Static 1-deg ocean wet fraction, written as `oceanWetFraction` next
   ! to mask_2d: valid_frac = sum_j w_ij*wet_j (NOT divided by the row
   ! sum), wet = open-surface ocean: maxLevelCell > 0 AND not under land
   ! ice (landIceMask == 1). The cavity rule is the one ocn_comp_mct uses
   ! to build the coupler ocean domain mask/frac, so this is the
   ! MPAS-side counterpart of the coupler's ocean fraction. Because the
   ! map is built on the culled ocean-only mesh, its row sums are < 1 at
   ! coastlines; keeping the absolute sum is what makes coastal cells
   ! partial. Range [0, 1]; 0 where no MPAS cell contributes.
   real(r8), allocatable :: ocean_wet_frac_remap(:,:)
   ! Inputs for the derived budget terms, set once at init.
   real(kind=RKIND) :: fme_ocn_dt = 0.0_RKIND          ! model time step (s)
   real(kind=RKIND) :: fme_frazil_heat_of_fusion = 0.0_RKIND
   ! landIceHeatFlux joins surfaceHeatFluxNet only when the land-ice flux
   ! mode actually computes it (conservationCheck's landIceFreshwaterFluxesOn
   ! rule). Under 'off'/'pressure_only' (IcoswISC30E3r5 default) it is not
   ! allocated and not part of the net.
   logical :: fme_landice_heat_in_net = .false.
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

   !> Optional 5-day-average companion stream. Mirrors the primary
   !> averaging path with its own accumulator/output_file at a longer
   !> cadence. Disabled when
   !> config_AM_fmeDerivedFields_stream_output_interval_5D is empty
   !> or 'none' (or when time_averaging is .false.).
   type(ocn_fme_remap_file_t) :: remap_outfile_5D
   logical :: time_avg_5D_enabled = .false.
   integer :: nAccum_5D = 0
   type(MPAS_TimeInterval_type) :: avg_output_interval_5D
   type(MPAS_Time_type) :: avg_last_output_time_5D
   real(r8), allocatable :: remap_accum_5D(:,:)        ! (n_b_local, N_AVG_FIELDS)
   real(r8), allocatable :: remap_accum_count_5D(:,:)  ! (n_b_local, N_AVG_FIELDS)

!***********************************************************************

contains

!***********************************************************************

   subroutine ocn_init_fme_derived_fields(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND), pointer :: config_AM_fmeDerivedFields_horiz_remap_file
      character(len=StrKIND), pointer :: config_stream_output_interval
      character(len=StrKIND), pointer :: config_stream_output_interval_5D
      logical, pointer :: config_time_averaging
      character(len=16) :: cm   ! CF cell_methods choice -- 'time: mean' / 'time: point'
      integer :: ierr
      ! long_name of surfaceHeatFluxNet; names every term (<= 128 chars,
      ! the remap layer's stored_long_names length).
      character(len=128) :: net_long_name

      err = 0

      ! dt / frazil Lf / land-ice rule are needed by compute even when the
      ! remap is off (native AM fields), so set them unconditionally.
      call set_budget_config()

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
                       'SW + LWdown + LH + SH only (not the net surface heat flux)', &
                       ierr, cell_methods=trim(cm))
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
                  call register_budget_vars(remap_outfile, trim(cm))
                  ! Static surface ocean mask: shape (lon, lat), 1=ocean,
                  ! 0=land. Written once per file. Matches the SamudrACE
                  ! convention (compute_ocean_dataset_e3sm.py:53).
                  call remap_outfile%register_var('mask_2d', '1', &
                       'Surface ocean wet mask (1=ocean, 0=land/invalid)', ierr, &
                       is_static=.true.)
                  call register_wet_fraction(remap_outfile)

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

                        ! Optional 5-day-average companion stream. Same
                        ! field set as the primary stream; averaged-only.
                        call mpas_pool_get_config(domain % configs, &
                             'config_AM_fmeDerivedFields_stream_output_interval_5D', &
                             config_stream_output_interval_5D)
                        if (associated(config_stream_output_interval_5D)) then
                           if (len_trim(config_stream_output_interval_5D) > 0 .and. &
                               trim(config_stream_output_interval_5D) /= 'none') then
                              remap_outfile_5D%stream_name = 'fmeDerivedFieldsOutput5D'
                              call remap_outfile_5D%register_var('sst', 'degC', &
                                   'Sea surface temperature', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('sss', 'PSU', &
                                   'Sea surface salinity', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('surfaceHeatFluxTotal', 'W m^-2', &
                                   'SW + LWdown + LH + SH only (not the net surface heat flux)', &
                                   ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('ssh', 'm', &
                                   'Sea surface height', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('shortWaveHeatFlux', 'W m^-2', &
                                   'Shortwave heat flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('longWaveHeatFluxDown', 'W m^-2', &
                                   'Downward longwave heat flux', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('latentHeatFlux', 'W m^-2', &
                                   'Latent heat flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('sensibleHeatFlux', 'W m^-2', &
                                   'Sensible heat flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('windStressZonal', 'N m^-2', &
                                   'Zonal wind stress on ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('windStressMeridional', 'N m^-2', &
                                   'Meridional wind stress on ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('longWaveHeatFluxUp', 'W m^-2', &
                                   'Upward longwave heat flux from ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('rainFlux', 'kg m^-2 s^-1', &
                                   'Rain mass flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('snowFlux', 'kg m^-2 s^-1', &
                                   'Snow mass flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('evaporationFlux', 'kg m^-2 s^-1', &
                                   'Evaporation mass flux from ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('seaIceHeatFlux', 'W m^-2', &
                                   'Sea-ice heat flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('icebergHeatFlux', 'W m^-2', &
                                   'Iceberg heat flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('seaIceFreshWaterFlux', 'kg m^-2 s^-1', &
                                   'Sea-ice fresh water mass flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('seaIceSalinityFlux', 'kg m^-2 s^-1', &
                                   'Sea-ice salt flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('riverRunoffFlux', 'kg m^-2 s^-1', &
                                   'River runoff mass flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('iceRunoffFlux', 'kg m^-2 s^-1', &
                                   'Land-ice runoff mass flux into ocean', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('snowFluxLf', 'W m^-2', &
                                   'snowFlux times latent heat of fusion', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('iceRunoffFluxLf', 'W m^-2', &
                                   'iceRunoffFlux times latent heat of fusion', ierr, cell_methods='time: mean')
                              call remap_outfile_5D%register_var('riverRunoffFluxTemperature', &
                                   'kg degC m^-2 s^-1', &
                                   'riverRunoffFlux times top-layer temperature', ierr, &
                                   cell_methods='time: mean')
                              call register_budget_vars(remap_outfile_5D, 'time: mean')
                              call remap_outfile_5D%register_var('mask_2d', '1', &
                                   'Surface ocean wet mask (1=ocean, 0=land/invalid)', ierr, &
                                   is_static=.true.)
                              call register_wet_fraction(remap_outfile_5D)

                              call mpas_set_timeInterval(avg_output_interval_5D, &
                                   timeString=trim(config_stream_output_interval_5D), &
                                   ierr=ierr)
                              avg_last_output_time_5D = mpas_get_clock_time( &
                                   domain % clock, MPAS_NOW, ierr)
                              if (.not. allocated(remap_accum_5D)) then
                                 allocate(remap_accum_5D(remap_nb_local, N_AVG_FIELDS))
                              end if
                              if (.not. allocated(remap_accum_count_5D)) then
                                 allocate(remap_accum_count_5D(remap_nb_local, N_AVG_FIELDS))
                              end if
                              remap_accum_5D = 0.0_r8
                              remap_accum_count_5D = 0.0_r8
                              nAccum_5D = 0
                              time_avg_5D_enabled = .true.
                              call mpas_log_write( &
                                   'fmeDerivedFields: 5-day average companion ' // &
                                   'stream enabled (interval=' // &
                                   trim(config_stream_output_interval_5D) // ')')
                           end if
                        end if

                        call ocn_fme_derived_fields_read_restart(domain, ierr)
                     end if
                  end if
               end if
            end if
         end if
      end if

   contains

      subroutine set_budget_config()
         ! Time step and frazil latent heat for frazilHeatFlux, and whether
         ! landIceHeatFlux is part of the net (same rule as
         ! conservationCheck's landIceFreshwaterFluxesOn).
         character(len=StrKIND), pointer :: config_lif_mode
         real(kind=RKIND), pointer :: config_frazil_lf
         type(MPAS_TimeInterval_type) :: dt_interval
         integer :: ierr_c
         ! Clock time step, exactly as ocn_frazil_forcing_build_arrays gets
         ! the dt it divides the accumulated mass increment by.
         dt_interval = mpas_get_clock_timestep(domain % clock, ierr=ierr_c)
         call mpas_get_timeInterval(dt_interval, dt=fme_ocn_dt)
         call mpas_pool_get_config(domain % configs, 'config_frazil_heat_of_fusion', &
              config_frazil_lf)
         if (associated(config_frazil_lf)) fme_frazil_heat_of_fusion = config_frazil_lf
         fme_landice_heat_in_net = .false.
         call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', &
              config_lif_mode)
         if (associated(config_lif_mode)) then
            fme_landice_heat_in_net = (trim(config_lif_mode) == 'standalone' .or. &
                 trim(config_lif_mode) == 'coupled' .or. trim(config_lif_mode) == 'data')
         end if
         if (fme_landice_heat_in_net) then
            net_long_name = 'Net surface heat flux into ocean: SW+LWdown+LWup+LH+SH' // &
                 '+seaIce+iceberg-Lf*snow-Lf*iceRunoff+frazil+landIce (per ocean-cell area)'
         else
            net_long_name = 'Net surface heat flux into ocean: SW+LWdown+LWup+LH+SH' // &
                 '+seaIce+iceberg-Lf*snow-Lf*iceRunoff+frazil (per ocean-cell area)'
         end if
      end subroutine set_budget_config

      subroutine register_budget_vars(outfile, cmeth)
         ! Slots 25-31. Same set and attributes in the 1D and 5D files.
         ! All surface fluxes here are per unit ocean-cell area (the
         ! ice-fraction-weighted x2o values), positive into the ocean.
         type(ocn_fme_remap_file_t), intent(inout) :: outfile
         character(len=*), intent(in) :: cmeth
         call outfile%register_var('surfaceHeatFluxNet', 'W m^-2', &
              trim(net_long_name), ierr, cell_methods=cmeth)
         call outfile%register_var('frazilHeatFlux', 'W m^-2', &
              'Frazil latent heat released into ocean: Lf_frazil*d(accumulatedFrazilIceMass)/dt', &
              ierr, cell_methods=cmeth)
         call outfile%register_var('massFluxHeatContent', 'W m^-2', &
              'rho*cp*T*(rain+evap+seaIce+iceberg+river) flux, T rel. 0 degC, positive into ocean; not in surfaceHeatFluxNet', &
              ierr, cell_methods=cmeth)
         call outfile%register_var('icebergFreshWaterFlux', 'kg m^-2 s^-1', &
              'Iceberg melt fresh water mass flux into ocean', ierr, cell_methods=cmeth)
         call outfile%register_var('frazilIceFreshwaterFlux', 'kg m^-2 s^-1', &
              'Frazil ice mass flux into ocean (under land ice only; negative = leaving)', &
              ierr, cell_methods=cmeth)
         call outfile%register_var('removedRiverRunoffFlux', 'kg m^-2 s^-1', &
              'River runoff removed by config_remove_ais_river_runoff (not added to ocean)', &
              ierr, cell_methods=cmeth)
         call outfile%register_var('removedIceRunoffFlux', 'kg m^-2 s^-1', &
              'Ice runoff removed by config_remove_ais_ice_runoff (not added to ocean)', &
              ierr, cell_methods=cmeth)
      end subroutine register_budget_vars

      subroutine register_wet_fraction(outfile)
         type(ocn_fme_remap_file_t), intent(inout) :: outfile
         call outfile%register_var('oceanWetFraction', '1', &
              'Remap-weighted fraction of MPAS sources that are open ocean ' // &
              '(maxLevelCell>0, landIceMask/=1)', ierr, is_static=.true.)
      end subroutine register_wet_fraction

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
      type(MPAS_Time_type)   :: now_time

      err = 0
      if (.not. time_avg_enabled) return

      call mpas_pool_get_config(domain % configs, 'config_do_restart', &
           config_do_restart)
      if (.not. associated(config_do_restart)) return
      if (.not. config_do_restart) return

      now_time = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr_local)

      call ocn_fme_horiz_remap_read_accum('fmeDerivedFields', remap_accum, &
           N_AVG_FIELDS, nAccum, ts_str, &
           count_arr=remap_accum_count, success=ok)
      call abort_if_unreadable('fmeDerivedFields', ok)

      if (ok .and. len_trim(ts_str) > 0) then
         call mpas_set_time(avg_last_output_time, &
              dateTimeString=trim(ts_str), ierr=ierr_local)
         if (ierr_local /= 0) then
            call mpas_log_write('fmeDerivedFields: failed to parse ' // &
                 'persisted avg_last_output_time="' // trim(ts_str) // &
                 '"; falling back to MPAS_NOW baseline.', &
                 messageType=MPAS_LOG_WARN)
            avg_last_output_time = now_time
         end if
      end if
      call rewind_reset(avg_last_output_time, now_time, &
           remap_accum, remap_accum_count, nAccum, ts_str)

      if (time_avg_5D_enabled) then
         call ocn_fme_horiz_remap_read_accum('fmeDerivedFields5D', &
              remap_accum_5D, N_AVG_FIELDS, nAccum_5D, ts_str, &
              count_arr=remap_accum_count_5D, success=ok)
         call abort_if_unreadable('fmeDerivedFields5D', ok)
         if (ok .and. len_trim(ts_str) > 0) then
            call mpas_set_time(avg_last_output_time_5D, &
                 dateTimeString=trim(ts_str), ierr=ierr_local)
            if (ierr_local /= 0) then
               call mpas_log_write('fmeDerivedFields: failed to parse ' // &
                    'persisted avg_last_output_time_5D="' // trim(ts_str) // &
                    '"; falling back to MPAS_NOW baseline.', &
                    messageType=MPAS_LOG_WARN)
               avg_last_output_time_5D = now_time
            end if
         end if
         call rewind_reset(avg_last_output_time_5D, now_time, &
              remap_accum_5D, remap_accum_count_5D, nAccum_5D, ts_str)
      end if

   contains

      subroutine abort_if_unreadable(am_name, read_ok)
         !> A sidecar that exists but cannot be restored is almost always
         !> one written with a different N_AVG_FIELDS (e.g. a 24-slot
         !> 0201-era sidecar read by this 31-slot build); read_accum logs
         !> the on-disk/expected counts as a WARN just before this. The
         !> old behaviour -- silently zero the accumulators and restart the
         !> window at MPAS_NOW -- drops the pre-restart part of the window
         !> and, for a restart off a window boundary, shifts every later
         !> record off the 1850-anchored 5-day grid. Abort instead; to
         !> start the averages fresh on purpose, remove the sidecar.
         !> A missing sidecar is still a silent fresh start (cold branch).
         character(len=*), intent(in) :: am_name
         logical, intent(in) :: read_ok
         logical :: sidecar_exists
         if (read_ok) return
         inquire(file=trim(am_name) // '.fme_accum_restart.nc', exist=sidecar_exists)
         if (.not. sidecar_exists) return
         call mpas_log_write('fmeDerivedFields: accumulator sidecar ' // &
              trim(am_name) // '.fme_accum_restart.nc exists but could not be ' // &
              'restored (slot-count mismatch or read error; this build has ' // &
              'N_AVG_FIELDS=$i). Remove the file to restart the averaging ' // &
              'window from this restart time.', &
              messageType=MPAS_LOG_CRIT, intArgs=(/N_AVG_FIELDS/))
      end subroutine abort_if_unreadable

      subroutine rewind_reset(last_out, now_t, accum, accum_count, nacc, label)
         !> Branch/rewind restart handling (AGENTS.md). The accumulator
         !> sidecar restores avg_last_output_time verbatim, assuming a
         !> restart can only move time forward. On a rewind to an earlier
         !> date the restored time is in the FUTURE relative to MPAS_NOW, so
         !> the output gate (currTime - last_out >= interval) is negative for
         !> the whole re-run and never fires -- the re-crossed months get no
         !> output. Detect that here (strict '>' leaves an exact same-point
         !> warm restart untouched, preserving warm-restart BFB), reset the
         !> schedule to MPAS_NOW, discard the gap-spanning accumulator (it
         !> would otherwise average the whole gap into the first record --
         !> Bug 2), and tell the remap layer to clobber-create the
         !> re-crossed files so their static masks/history are rewritten.
         type(MPAS_Time_type), intent(inout) :: last_out
         type(MPAS_Time_type), intent(in)    :: now_t
         real(r8), intent(inout) :: accum(:,:), accum_count(:,:)
         integer,  intent(inout) :: nacc
         character(len=*), intent(in) :: label
         if (last_out > now_t) then
            call mpas_log_write('fmeDerivedFields: rewind restart detected ' // &
                 '(persisted ' // trim(label) // ' is later than MPAS_NOW); ' // &
                 'resetting schedule to MPAS_NOW and discarding stale ' // &
                 'accumulator.', messageType=MPAS_LOG_WARN)
            last_out    = now_t
            accum       = 0.0_r8
            accum_count = 0.0_r8
            nacc        = 0
            call ocn_fme_horiz_remap_set_rewind()
         end if
      end subroutine rewind_reset

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

      if (time_avg_5D_enabled .and. allocated(remap_accum_5D)) then
         call mpas_get_time(avg_last_output_time_5D, &
              dateTimeString=ts_str, ierr=ierr_local)
         call ocn_fme_horiz_remap_write_accum('fmeDerivedFields5D', &
              remap_accum_5D, N_AVG_FIELDS, nAccum_5D, trim(ts_str), &
              count_arr=remap_accum_count_5D, err=err)
      end if

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
      ! 0301 budget terms (slots 25-31)
      real(kind=RKIND), dimension(:), pointer :: surfaceHeatFluxNet, frazilHeatFlux
      real(kind=RKIND), dimension(:), pointer :: massFluxHeatContent
      real(kind=RKIND), dimension(:), pointer :: icebergFreshWaterFluxPtr
      real(kind=RKIND), dimension(:), pointer :: frazilIceFreshwaterFluxPtr
      real(kind=RKIND), dimension(:), pointer :: removedRiverRunoffFluxPtr
      real(kind=RKIND), dimension(:), pointer :: removedIceRunoffFluxPtr
      real(kind=RKIND), dimension(:), pointer :: landIceHeatFluxPtr
      real(kind=RKIND), dimension(:), pointer :: frazilMassLev2, frazilMassLev1
      real(kind=RKIND), dimension(:), pointer :: rainTemperatureFluxPtr
      real(kind=RKIND), dimension(:), pointer :: evapTemperatureFluxPtr
      real(kind=RKIND), dimension(:), pointer :: seaIceTemperatureFluxPtr
      real(kind=RKIND), dimension(:), pointer :: icebergTemperatureFluxPtr
      real(kind=RKIND), dimension(:,:), pointer :: runoffTracerFlux
      integer, dimension(:), pointer :: landIceMask
      integer, pointer :: index_temperature_sflux
      type(mpas_pool_type), pointer :: tracersSurfaceFluxPool
      logical :: have_net_terms
      real(kind=RKIND) :: tflux_sum
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

         call mpas_pool_get_array(amPool, 'surfaceHeatFluxNet', surfaceHeatFluxNet)
         call mpas_pool_get_array(amPool, 'frazilHeatFlux', frazilHeatFlux)
         call mpas_pool_get_array(amPool, 'massFluxHeatContent', massFluxHeatContent)
         call mpas_pool_get_array(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFluxPtr)
         call mpas_pool_get_array(forcingPool, 'frazilIceFreshwaterFlux', &
              frazilIceFreshwaterFluxPtr)
         call mpas_pool_get_array(forcingPool, 'removedRiverRunoffFlux', &
              removedRiverRunoffFluxPtr)
         call mpas_pool_get_array(forcingPool, 'removedIceRunoffFlux', removedIceRunoffFluxPtr)
         call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFluxPtr)
         call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
         call mpas_pool_get_array(forcingPool, 'rainTemperatureFlux', rainTemperatureFluxPtr)
         call mpas_pool_get_array(forcingPool, 'evapTemperatureFlux', evapTemperatureFluxPtr)
         call mpas_pool_get_array(forcingPool, 'seaIceTemperatureFlux', seaIceTemperatureFluxPtr)
         call mpas_pool_get_array(forcingPool, 'icebergTemperatureFlux', &
              icebergTemperatureFluxPtr)
         nullify(runoffTracerFlux, index_temperature_sflux)
         call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
         if (associated(tracersSurfaceFluxPool)) then
            call mpas_pool_get_dimension(tracersSurfaceFluxPool, &
                 'index_temperatureSurfaceFlux', index_temperature_sflux)
            call mpas_pool_get_array(tracersSurfaceFluxPool, &
                 'activeTracersSurfaceFluxRunoff', runoffTracerFlux)
            if (.not. associated(index_temperature_sflux)) nullify(runoffTracerFlux)
         end if
         ! Frazil mass. The AM runs AFTER mpas_pool_shift_time_levels, so
         ! time level 1 holds this step's updated mass and time level 2 the
         ! value the step started from (after the previous coupler send
         ! reset it, ocn_comp_mct export). Their difference is exactly the
         ! mass formed this step whether or not a reset happened in
         ! between. The formula below is conservationCheck's sumArray(10)
         ! verbatim: it reads the same levels post-shift (it calls level 2
         ! 'New' and level 1 'Old') and negates the level2-level1 difference.
         call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', frazilMassLev2, 2)
         call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', frazilMassLev1, 1)

         ! The net needs every one of the nine spec components; frazil and
         ! landIce enter when present (zero/absent otherwise).
         have_net_terms = associated(shortWaveHeatFlux) .and. &
              associated(longWaveHeatFluxDown) .and. associated(longWaveHeatFluxUp) .and. &
              associated(latentHeatFlux) .and. associated(sensibleHeatFlux) .and. &
              associated(seaIceHeatFluxPtr) .and. associated(icebergHeatFluxPtr) .and. &
              associated(snowFluxPtr) .and. associated(iceRunoffFluxPtr)

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
               surfaceHeatFluxNet(iCell) = fillValue
               frazilHeatFlux(iCell) = fillValue
               massFluxHeatContent(iCell) = fillValue
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

            ! Frazil heat (conservationCheck sumArray(10) per cell): mass
            ! formed this step times the frazil latent heat, positive into
            ! the ocean. Zero when frazil formation is off.
            frazilHeatFlux(iCell) = 0.0_RKIND
            if (associated(frazilMassLev2) .and. associated(frazilMassLev1) .and. &
                fme_ocn_dt > 0.0_RKIND) then
               frazilHeatFlux(iCell) = - fme_frazil_heat_of_fusion &
                    * (frazilMassLev2(iCell) - frazilMassLev1(iCell)) / fme_ocn_dt
            end if

            ! Net surface heat flux, the netEnergyFlux heat terms of
            ! conservationCheck without its rho*cp*temperature-flux terms
            ! (those go to massFluxHeatContent). Snow and ice runoff melt
            ! cool the ocean, hence -Lf*flux. All terms are per unit
            ! ocean-cell area (ice-fraction-weighted coupler values).
            if (have_net_terms) then
               surfaceHeatFluxNet(iCell) = shortWaveHeatFlux(iCell) &
                    + longWaveHeatFluxDown(iCell) &
                    + longWaveHeatFluxUp(iCell) &
                    + latentHeatFlux(iCell) &
                    + sensibleHeatFlux(iCell) &
                    + seaIceHeatFluxPtr(iCell) &
                    + icebergHeatFluxPtr(iCell) &
                    - snowFluxPtr(iCell) * latent_heat_fusion_mks &
                    - iceRunoffFluxPtr(iCell) * latent_heat_fusion_mks &
                    + frazilHeatFlux(iCell)
               if (fme_landice_heat_in_net .and. associated(landIceHeatFluxPtr)) &
                    surfaceHeatFluxNet(iCell) = surfaceHeatFluxNet(iCell) &
                    + landIceHeatFluxPtr(iCell)
            else
               surfaceHeatFluxNet(iCell) = fillValue
            end if

            ! Heat carried by the surface water fluxes relative to 0 degC
            ! (conservationCheck sumArray(11-15) * rho_sw*cp_sw). Kept out
            ! of the net: its value depends on the 0 degC reference, but
            ! the 0 degC-referenced ocean heat content budget needs it.
            tflux_sum = 0.0_RKIND
            if (associated(rainTemperatureFluxPtr)) &
                 tflux_sum = tflux_sum + rainTemperatureFluxPtr(iCell)
            if (associated(evapTemperatureFluxPtr)) &
                 tflux_sum = tflux_sum + evapTemperatureFluxPtr(iCell)
            if (associated(seaIceTemperatureFluxPtr)) &
                 tflux_sum = tflux_sum + seaIceTemperatureFluxPtr(iCell)
            if (associated(icebergTemperatureFluxPtr)) &
                 tflux_sum = tflux_sum + icebergTemperatureFluxPtr(iCell)
            if (associated(runoffTracerFlux)) &
                 tflux_sum = tflux_sum + runoffTracerFlux(index_temperature_sflux, iCell)
            massFluxHeatContent(iCell) = tflux_sum * rho_sw * cp_sw
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
               call remap_1field(surfaceHeatFluxNet, 25)
               call remap_1field(frazilHeatFlux, 26)
               call remap_1field(massFluxHeatContent, 27)
               if (associated(icebergFreshWaterFluxPtr)) &
                    call mask_and_remap(icebergFreshWaterFluxPtr, 28)
               if (associated(frazilIceFreshwaterFluxPtr)) &
                    call mask_and_remap(frazilIceFreshwaterFluxPtr, 29)
               if (associated(removedRiverRunoffFluxPtr)) &
                    call mask_and_remap(removedRiverRunoffFluxPtr, 30)
               if (associated(removedIceRunoffFluxPtr)) &
                    call mask_and_remap(removedIceRunoffFluxPtr, 31)
               nAccum = nAccum + 1
               if (time_avg_5D_enabled) nAccum_5D = nAccum_5D + 1

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
                     call write_avg_field('surfaceHeatFluxNet', 25)
                     call write_avg_field('frazilHeatFlux', 26)
                     call write_avg_field('massFluxHeatContent', 27)
                     call write_avg_field('icebergFreshWaterFlux', 28)
                     call write_avg_field('frazilIceFreshwaterFlux', 29)
                     call write_avg_field('removedRiverRunoffFlux', 30)
                     call write_avg_field('removedIceRunoffFlux', 31)
                  end if
                  remap_accum = 0.0_r8
                  remap_accum_count = 0.0_r8
                  nAccum = 0
                  ! Advance by exact interval to avoid drift.
                  avg_last_output_time = avg_last_output_time + avg_output_interval
               end if

               ! 5-day-average companion: parallel to the 1D path above,
               ! shares per-step samples but uses its own accumulator,
               ! last-output-time, output_interval, and output file.
               if (time_avg_5D_enabled .and. &
                    currTime - avg_last_output_time_5D >= avg_output_interval_5D) then
                  call remap_outfile_5D%check_rotate(domain, ierr)
                  if (remap_outfile_5D%file_is_open) then
                     call remap_outfile_5D%write_time(avg_last_output_time_5D, &
                          currTime, ierr)
                     call write_static_masks_5D_if_needed()
                     call write_avg_field_5D('sst', 1)
                     call write_avg_field_5D('sss', 2)
                     call write_avg_field_5D('surfaceHeatFluxTotal', 3)
                     call write_avg_field_5D('ssh', 4)
                     call write_avg_field_5D('shortWaveHeatFlux', 5)
                     call write_avg_field_5D('longWaveHeatFluxDown', 6)
                     call write_avg_field_5D('latentHeatFlux', 7)
                     call write_avg_field_5D('sensibleHeatFlux', 8)
                     call write_avg_field_5D('windStressZonal', 9)
                     call write_avg_field_5D('windStressMeridional', 10)
                     call write_avg_field_5D('longWaveHeatFluxUp', 11)
                     ! idx=12 unused (matches primary stream)
                     call write_avg_field_5D('rainFlux', 13)
                     call write_avg_field_5D('snowFlux', 14)
                     call write_avg_field_5D('evaporationFlux', 15)
                     call write_avg_field_5D('seaIceHeatFlux', 16)
                     call write_avg_field_5D('icebergHeatFlux', 17)
                     call write_avg_field_5D('seaIceFreshWaterFlux', 18)
                     call write_avg_field_5D('seaIceSalinityFlux', 19)
                     call write_avg_field_5D('riverRunoffFlux', 20)
                     call write_avg_field_5D('iceRunoffFlux', 21)
                     call write_avg_field_5D('snowFluxLf', 22)
                     call write_avg_field_5D('iceRunoffFluxLf', 23)
                     call write_avg_field_5D('riverRunoffFluxTemperature', 24)
                     call write_avg_field_5D('surfaceHeatFluxNet', 25)
                     call write_avg_field_5D('frazilHeatFlux', 26)
                     call write_avg_field_5D('massFluxHeatContent', 27)
                     call write_avg_field_5D('icebergFreshWaterFlux', 28)
                     call write_avg_field_5D('frazilIceFreshwaterFlux', 29)
                     call write_avg_field_5D('removedRiverRunoffFlux', 30)
                     call write_avg_field_5D('removedIceRunoffFlux', 31)
                  end if
                  remap_accum_5D = 0.0_r8
                  remap_accum_count_5D = 0.0_r8
                  nAccum_5D = 0
                  avg_last_output_time_5D = avg_last_output_time_5D + avg_output_interval_5D
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
                  call inst_remap_and_write(surfaceHeatFluxNet, 'surfaceHeatFluxNet', .false.)
                  call inst_remap_and_write(frazilHeatFlux, 'frazilHeatFlux', .false.)
                  call inst_remap_and_write(massFluxHeatContent, 'massFluxHeatContent', .false.)
                  if (associated(icebergFreshWaterFluxPtr)) &
                       call inst_remap_and_write(icebergFreshWaterFluxPtr, &
                            'icebergFreshWaterFlux', .true.)
                  if (associated(frazilIceFreshwaterFluxPtr)) &
                       call inst_remap_and_write(frazilIceFreshwaterFluxPtr, &
                            'frazilIceFreshwaterFlux', .true.)
                  if (associated(removedRiverRunoffFluxPtr)) &
                       call inst_remap_and_write(removedRiverRunoffFluxPtr, &
                            'removedRiverRunoffFlux', .true.)
                  if (associated(removedIceRunoffFluxPtr)) &
                       call inst_remap_and_write(removedIceRunoffFluxPtr, &
                            'removedIceRunoffFlux', .true.)
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
            if (time_avg_5D_enabled) then
               do ib = 1, remap_nb_local
                  if (abs(fld_remap(ib, 1)) < fill_thresh) then
                     remap_accum_5D(ib, idx)       = remap_accum_5D(ib, idx)       + fld_remap(ib, 1)
                     remap_accum_count_5D(ib, idx) = remap_accum_count_5D(ib, idx) + 1.0_r8
                  end if
               end do
            end if
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

      subroutine write_avg_field_5D(varname, idx)
         character(len=*), intent(in) :: varname
         integer, intent(in) :: idx
         integer :: ib
         do ib = 1, remap_nb_local
            if (remap_accum_count_5D(ib, idx) > 0.0_r8) then
               fld_remap(ib, 1) = remap_accum_5D(ib, idx) / remap_accum_count_5D(ib, idx)
            else
               fld_remap(ib, 1) = real(SHR_FILL_VALUE, r8)
            end if
         end do
         call remap_outfile_5D%write_var(varname, fld_remap, 1, ierr)
      end subroutine write_avg_field_5D

      subroutine compute_mask_2d_remap()
         ! Build the static surface ocean mask in remap space by remapping
         ! a constant-1 ocean field through apply_masked. Land cells
         ! (maxLevelCell == 0) have mask=0 in the remap, so target cells
         ! whose source-side coverage is entirely land come out as
         ! SHR_FILL_VALUE; everything else comes out as 1.0 (a weighted
         ! average of constant 1s = 1). Threshold to a hard 0.0/1.0.
         integer :: ic, ib, ierr_cov
         real(r8), allocatable :: coverage(:)
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

         ! oceanWetFraction = valid_frac = sum_j w_j*wet_j. Remap the
         ! open-ocean indicator with every MPAS cell valid (no mask
         ! argument): apply_masked then divides by the full row sum,
         ! giving sum w*wet / sum w, and multiplying back by the row sum
         ! (coverage) restores the absolute sum. Dividing by the row sum
         ! alone would set every coastal 1-deg cell the culled mesh only
         ! partly covers to 1 -- exactly the cells spec check 6 is about.
         ! Cavity cells (landIceMask==1) and any maxLevelCell==0 cell
         ! count as dry, matching the coupler ocean domain mask built in
         ! ocn_comp_mct. Target cells with no MPAS source come back as
         ! fill and are written as 0.
         allocate(ocean_wet_frac_remap(remap_nb_local, 1))
         allocate(coverage(remap_nb_local))
         call ocn_fme_horiz_remap_get_coverage(coverage, ierr_cov)
         if (ierr_cov /= 0) then
            call mpas_log_write('fmeDerivedFields: remap row sums unavailable ' // &
                 '(err=$i); oceanWetFraction not written', &
                 messageType=MPAS_LOG_WARN, intArgs=(/ierr_cov/))
         end if
         do ic = 1, nCells
            field_1lev(1, ic) = 0.0_RKIND
            if (maxLevelCell(ic) > 0) field_1lev(1, ic) = 1.0_RKIND
            if (associated(landIceMask)) then
               if (landIceMask(ic) == 1) field_1lev(1, ic) = 0.0_RKIND
            end if
         end do
         ! Collective (alltoallv): every rank calls it, whatever ierr_cov.
         call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr)
         do ib = 1, remap_nb_local
            if (abs(fld_remap(ib, 1)) < MASK_FILL_THRESH) then
               ! Clamp round-off only (row sums can exceed 1 by ~1e-15).
               ocean_wet_frac_remap(ib, 1) = &
                    min(1.0_r8, max(0.0_r8, fld_remap(ib, 1) * coverage(ib)))
            else
               ocean_wet_frac_remap(ib, 1) = 0.0_r8
            end if
         end do
         deallocate(coverage)
         ! No coverage => do not write a misleading all-zero field.
         if (ierr_cov /= 0) deallocate(ocean_wet_frac_remap)
      end subroutine compute_mask_2d_remap

      subroutine write_static_masks_if_needed()
         ! Write the static mask on the FIRST data record of the file,
         ! whether the file was freshly created OR append-reopened with no
         ! records yet. The latter case matters: if a prior leg created the
         ! monthly file (defining mask_2d) but stopped before its first
         ! flush, mask_2d was left defined-but-unwritten (= all-zero), and
         ! the reopening leg must fill it. A normal warm-restart reopen of a
         ! file that already holds records has time_record >= 2 here, so its
         ! correct mask is left untouched (warm-restart BFB preserved).
         ! mask_2d_remap is recomputed each leg from the static land mask,
         ! so the value written on reopen is correct. (Observed in
         ! production: fmeDerivedFields.0451-01 shipped an all-zero mask
         ! exactly this way — see AGENTS.md gotcha #45.)
         if (remap_outfile%time_record /= 1) return
         if (.not. allocated(mask_2d_remap)) return
         call remap_outfile%write_var('mask_2d', mask_2d_remap, 1, ierr)
         if (allocated(ocean_wet_frac_remap)) &
              call remap_outfile%write_var('oceanWetFraction', ocean_wet_frac_remap, 1, ierr)
      end subroutine write_static_masks_if_needed

      subroutine write_static_masks_5D_if_needed()
         if (remap_outfile_5D%time_record /= 1) return
         if (.not. allocated(mask_2d_remap)) return
         call remap_outfile_5D%write_var('mask_2d', mask_2d_remap, 1, ierr)
         if (allocated(ocean_wet_frac_remap)) &
              call remap_outfile_5D%write_var('oceanWetFraction', ocean_wet_frac_remap, 1, ierr)
      end subroutine write_static_masks_5D_if_needed

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
      if (time_avg_5D_enabled) call remap_outfile_5D%close_file(ierr)
      if (allocated(mask_2d_remap)) deallocate(mask_2d_remap)
      if (allocated(ocean_wet_frac_remap)) deallocate(ocean_wet_frac_remap)
      if (allocated(remap_accum_5D)) deallocate(remap_accum_5D)
      if (allocated(remap_accum_count_5D)) deallocate(remap_accum_count_5D)

   end subroutine ocn_finalize_fme_derived_fields!}}}

end module ocn_fme_derived_fields

! vim: foldmethod=marker
