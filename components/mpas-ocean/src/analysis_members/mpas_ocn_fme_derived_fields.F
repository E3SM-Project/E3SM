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
        ocn_fme_horiz_remap_field, ocn_fme_horiz_remap_active, &
        ocn_fme_horiz_remap_get_dims, ocn_fme_remap_file_t, &
        ocn_fme_remap_build_filename

   implicit none
   private
   save

   public :: ocn_init_fme_derived_fields, &
             ocn_compute_fme_derived_fields, &
             ocn_restart_fme_derived_fields, &
             ocn_finalize_fme_derived_fields

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
   real(r8), allocatable :: remap_accum(:,:)  ! (n_b_local, N_AVG_FIELDS)

!***********************************************************************

contains

!***********************************************************************

   subroutine ocn_init_fme_derived_fields(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND), pointer :: config_AM_fmeDerivedFields_horiz_remap_file
      character(len=StrKIND), pointer :: config_stream_output_interval
      logical, pointer :: config_time_averaging
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
               do_horiz_remap = ocn_fme_horiz_remap_active()
               if (do_horiz_remap) then
                  call ocn_fme_horiz_remap_get_dims(remap_nlat, remap_nlon, remap_nb_local)
                  call mpas_log_write('fmeDerivedFields: horiz remap active, ' &
                       // 'target grid $i x $i', &
                       intArgs=(/remap_nlon, remap_nlat/))
                  remap_outfile%stream_name = 'fmeDerivedFieldsOutput'
                  call remap_outfile%register_var('sst', 'degC', &
                       'Sea surface temperature', ierr)
                  call remap_outfile%register_var('sss', 'PSU', &
                       'Sea surface salinity', ierr)
                  call remap_outfile%register_var('surfaceHeatFluxTotal', 'W m^-2', &
                       'Total surface heat flux', ierr)
                  call remap_outfile%register_var('ssh', 'm', &
                       'Sea surface height', ierr)
                  call remap_outfile%register_var('shortWaveHeatFlux', 'W m^-2', &
                       'Shortwave heat flux into ocean', ierr)
                  call remap_outfile%register_var('longWaveHeatFluxDown', 'W m^-2', &
                       'Downward longwave heat flux', ierr)
                  call remap_outfile%register_var('latentHeatFlux', 'W m^-2', &
                       'Latent heat flux into ocean', ierr)
                  call remap_outfile%register_var('sensibleHeatFlux', 'W m^-2', &
                       'Sensible heat flux into ocean', ierr)
                  call remap_outfile%register_var('windStressZonal', 'N m^-2', &
                       'Zonal wind stress on ocean', ierr)
                  call remap_outfile%register_var('windStressMeridional', 'N m^-2', &
                       'Meridional wind stress on ocean', ierr)
                  call remap_outfile%register_var('longWaveHeatFluxUp', 'W m^-2', &
                       'Upward longwave heat flux from ocean', ierr)
                  call remap_outfile%register_var('rainFlux', 'kg m^-2 s^-1', &
                       'Rain mass flux into ocean', ierr)
                  call remap_outfile%register_var('snowFlux', 'kg m^-2 s^-1', &
                       'Snow mass flux into ocean', ierr)
                  call remap_outfile%register_var('evaporationFlux', 'kg m^-2 s^-1', &
                       'Evaporation mass flux from ocean', ierr)
                  call remap_outfile%register_var('seaIceHeatFlux', 'W m^-2', &
                       'Sea-ice heat flux into ocean', ierr)
                  call remap_outfile%register_var('icebergHeatFlux', 'W m^-2', &
                       'Iceberg heat flux into ocean', ierr)
                  call remap_outfile%register_var('seaIceFreshWaterFlux', 'kg m^-2 s^-1', &
                       'Sea-ice fresh water mass flux into ocean', ierr)
                  call remap_outfile%register_var('seaIceSalinityFlux', 'kg m^-2 s^-1', &
                       'Sea-ice salt flux into ocean', ierr)
                  call remap_outfile%register_var('riverRunoffFlux', 'kg m^-2 s^-1', &
                       'River runoff mass flux into ocean', ierr)
                  call remap_outfile%register_var('iceRunoffFlux', 'kg m^-2 s^-1', &
                       'Land-ice runoff mass flux into ocean', ierr)
                  call remap_outfile%register_var('snowFluxLf', 'W m^-2', &
                       'snowFlux times latent heat of fusion', ierr)
                  call remap_outfile%register_var('iceRunoffFluxLf', 'W m^-2', &
                       'iceRunoffFlux times latent heat of fusion', ierr)
                  call remap_outfile%register_var('riverRunoffFluxTemperature', &
                       'kg degC m^-2 s^-1', &
                       'riverRunoffFlux times top-layer temperature', ierr)

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
                        ! MPAS_NOW (not MPAS_START_TIME) so warm restart uses
                        ! the restart-time as the averaging baseline; otherwise
                        ! the trigger would fire every timestep until catching up.
                        ! Partial pre-restart samples are not persisted (the
                        ! accumulators live in remap-space and aren't in the
                        ! Registry); the post-restart window starts fresh.
                        avg_last_output_time = mpas_get_clock_time( &
                             domain % clock, MPAS_NOW, ierr)
                        if (.not. allocated(remap_accum)) then
                           allocate(remap_accum(remap_nb_local, N_AVG_FIELDS))
                        end if
                        remap_accum = 0.0_r8
                        nAccum = 0
                        time_avg_enabled = .true.
                        call mpas_log_write('fmeDerivedFields: time averaging enabled')
                     end if
                  end if
               end if
            end if
         end if
      end if

   end subroutine ocn_init_fme_derived_fields!}}}

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
      character(len=StrKIND) :: timeStr
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
                  call mpas_get_time(currTime, dateTimeString=timeStr)
                  call remap_outfile%check_rotate(domain, ierr)
                  if (remap_outfile%file_is_open) then
                     call remap_outfile%write_time(trim(timeStr), ierr)
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
                  nAccum = 0
                  ! Advance by exact interval to avoid drift.
                  avg_last_output_time = avg_last_output_time + avg_output_interval
               end if
            else
               ! Instantaneous: remap and write each field one at a time so
               ! fld_remap holds the correct field at every write call.
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
               call mpas_get_time(currTime, dateTimeString=timeStr)
               call remap_outfile%check_rotate(domain, ierr)
               if (remap_outfile%file_is_open) then
                  call remap_outfile%write_time(trim(timeStr), ierr)
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
         real(kind=RKIND), dimension(:), intent(in) :: native_fld
         integer, intent(in) :: idx
         field_1lev(1, 1:nCells) = native_fld(1:nCells)
         call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
              mask=maxLevelCell)
         if (time_avg_enabled) then
            remap_accum(:, idx) = remap_accum(:, idx) + fld_remap(:, 1)
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
         character(len=*), intent(in) :: varname
         integer, intent(in) :: idx
         if (nAccum > 0) then
            fld_remap(:, 1) = remap_accum(:, idx) / real(nAccum, r8)
         else
            fld_remap(:, 1) = real(SHR_FILL_VALUE, r8)
         end if
         call remap_outfile%write_var(varname, fld_remap, 1, ierr)
      end subroutine write_avg_field

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

   end subroutine ocn_finalize_fme_derived_fields!}}}

end module ocn_fme_derived_fields

! vim: foldmethod=marker
