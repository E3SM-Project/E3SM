!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_fme_derived_fields
!
!> \brief MPAS sea ice analysis mode member: fme_derived_fields
!> \author FME Development Team
!> \date   April 2026
!> \details
!>  MPAS sea ice analysis mode member: fme_derived_fields
!>  Computes derived fields for Full Model Emulation (FME) by
!>  aggregating category-level tracers into cell-level quantities.
!>
!-----------------------------------------------------------------------

module seaice_fme_derived_fields

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_kind_types
   use mpas_log, only: mpas_log_write

   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_horiz_remap_mod, only: SHR_FILL_VALUE
   use seaice_fme_horiz_remap, only: seaice_fme_horiz_remap_init, &
        seaice_fme_horiz_remap_init_time_epoch, &
        seaice_fme_horiz_remap_field, seaice_fme_horiz_remap_active, &
        seaice_fme_horiz_remap_write_var, seaice_fme_horiz_remap_write_time, &
        seaice_fme_horiz_remap_close_file, seaice_fme_horiz_remap_get_dims, &
        seaice_fme_horiz_remap_register_var, seaice_fme_horiz_remap_check_rotate, &
        seaice_fme_horiz_remap_write_accum, seaice_fme_horiz_remap_read_accum, &
        remap_stream_name, file_is_open

   implicit none
   private
   save

   ! Remap state
   logical :: do_horiz_remap = .false.
   integer :: remap_nlat = 0, remap_nlon = 0, remap_nb_local = 0

   ! Time-averaging state for remapped output
   logical :: time_avg_enabled = .false.
   integer :: nAccum = 0
   type(MPAS_TimeInterval_type) :: avg_output_interval
   type(MPAS_Time_type) :: avg_last_output_time
   ! Slot assignments for remap_accum (must stay in sync with init/compute):
   !   1 iceAreaTotal             6 airStressZonal
   !   2 iceVolumeTotal           7 airStressMeridional
   !   3 snowVolumeTotal          8 uVelocityGeoCell
   !   4 iceThicknessMean         9 vVelocityGeoCell
   !   5 surfaceTemperatureMean
   integer, parameter :: N_AVG_FIELDS = 9
   real(r8), allocatable :: remap_accum(:,:)        ! (remap_nb_local, N_AVG_FIELDS)
   ! Per-cell, per-field valid-sample counter parallel to remap_accum.
   ! Used so that ice-edge cells (where source-side ice presence flips
   ! between samples) average only over the timesteps where the remapped
   ! value was non-fill — without this, SHR_FILL_VALUE samples are summed
   ! into remap_accum and surface as fractional-of-fill pollution
   ! (k * 1e20 / nAccum) at the ice edge. See AGENTS.md gotcha on
   ! time-averaging fill pollution.
   real(r8), allocatable :: remap_accum_count(:,:)  ! (remap_nb_local, N_AVG_FIELDS)

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: seaice_bootstrap_fme_derived_fields, &
             seaice_init_fme_derived_fields, &
             seaice_precompute_fme_derived_fields, &
             seaice_compute_fme_derived_fields, &
             seaice_restart_fme_derived_fields, &
             seaice_finalize_fme_derived_fields, &
             seaice_fme_derived_fields_write_restart

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_fme_derived_fields
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  FME Development Team
!> \date    April 2026
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_fme_derived_fields(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_bootstrap_fme_derived_fields!}}}

!***********************************************************************
!
!  routine seaice_init_fme_derived_fields
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  FME Development Team
!> \date    April 2026
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_fme_derived_fields(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      character(len=StrKIND), pointer :: config_horiz_remap_file
      logical, pointer :: config_time_averaging
      character(len=StrKIND), pointer :: config_stream_output_interval
      integer :: ierr

      err = 0

      ! Initialize horizontal remapping if configured
      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeSeaiceDerivedFields_horiz_remap_file', &
           config_horiz_remap_file)

      if (associated(config_horiz_remap_file)) then
         if (trim(config_horiz_remap_file) /= 'none' .and. &
             len_trim(config_horiz_remap_file) > 0) then
            call seaice_fme_horiz_remap_init(domain, &
                 trim(config_horiz_remap_file), ierr)
            if (ierr == 0) then
               ! Auto-discover the CF time epoch.
               call seaice_fme_horiz_remap_init_time_epoch(domain, &
                    'fmeSeaiceDerivedFields', ierr)
               if (ierr /= 0) then
                  call mpas_log_write('fmeSeaiceDerivedFields: time epoch ' // &
                       'setup failed, disabling horiz remap', &
                       messageType=MPAS_LOG_WARN)
                  do_horiz_remap = .false.
                  err = ierr
                  return
               end if
               do_horiz_remap = seaice_fme_horiz_remap_active()
               if (do_horiz_remap) then
                  call seaice_fme_horiz_remap_get_dims( &
                       remap_nlat, remap_nlon, remap_nb_local)
                  call mpas_log_write( &
                       'fmeSeaiceDerivedFields: horiz remap active, target $i x $i', &
                       intArgs=(/remap_nlon, remap_nlat/))
                  remap_stream_name = 'fmeSeaiceDerivedFieldsOutput'
                  call seaice_fme_horiz_remap_register_var('iceAreaTotal', &
                       'unitless', 'Total ice area fraction', ierr)
                  call seaice_fme_horiz_remap_register_var('iceVolumeTotal', &
                       'm', 'Total ice volume per unit area', ierr)
                  call seaice_fme_horiz_remap_register_var('snowVolumeTotal', &
                       'm', 'Total snow volume per unit area', ierr)
                  call seaice_fme_horiz_remap_register_var('iceThicknessMean', &
                       'm', 'Area-weighted mean ice thickness', ierr)
                  call seaice_fme_horiz_remap_register_var('surfaceTemperatureMean', &
                       'K', 'Area-weighted mean surface temperature', ierr)
                  call seaice_fme_horiz_remap_register_var('airStressZonal', &
                       'N m^-2', 'Zonal air stress on ice', ierr)
                  call seaice_fme_horiz_remap_register_var('airStressMeridional', &
                       'N m^-2', 'Meridional air stress on ice', ierr)
                  call seaice_fme_horiz_remap_register_var('uVelocityGeoCell', &
                       'm s^-1', 'True eastwards ice velocity at cell centers', ierr)
                  call seaice_fme_horiz_remap_register_var('vVelocityGeoCell', &
                       'm s^-1', 'True northwards ice velocity at cell centers', ierr)

                  ! Time-averaging setup
                  call mpas_pool_get_config(domain % configs, &
                       'config_AM_fmeSeaiceDerivedFields_time_averaging', &
                       config_time_averaging)
                  if (associated(config_time_averaging)) then
                     if (config_time_averaging) then
                        call mpas_pool_get_config(domain % configs, &
                             'config_AM_fmeSeaiceDerivedFields_stream_output_interval', &
                             config_stream_output_interval)
                        call mpas_set_timeInterval(avg_output_interval, &
                             timeString=trim(config_stream_output_interval), &
                             ierr=ierr)
                        ! Default baseline is MPAS_NOW; persisted state
                        ! below overrides on warm restart with sidecar.
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
                        call mpas_log_write( &
                             'fmeSeaiceDerivedFields: time averaging enabled')
                        call seaice_fme_derived_fields_read_restart(domain, ierr)
                     end if
                  end if
               end if
            end if
         end if
      end if

   end subroutine seaice_init_fme_derived_fields!}}}

!***********************************************************************

   subroutine seaice_fme_derived_fields_read_restart(domain, err)!{{{
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

      call seaice_fme_horiz_remap_read_accum('fmeSeaiceDerivedFields', &
           remap_accum, N_AVG_FIELDS, nAccum, ts_str, &
           count_arr=remap_accum_count, success=ok)

      if (ok .and. len_trim(ts_str) > 0) then
         call mpas_set_time(avg_last_output_time, &
              dateTimeString=trim(ts_str), ierr=ierr_local)
         if (ierr_local /= 0) then
            call mpas_log_write('fmeSeaiceDerivedFields: failed to parse ' // &
                 'persisted avg_last_output_time="' // trim(ts_str) // &
                 '"; falling back to MPAS_NOW baseline.', &
                 messageType=MPAS_LOG_WARN)
            avg_last_output_time = mpas_get_clock_time( &
                 domain % clock, MPAS_NOW, ierr_local)
         end if
      end if

   end subroutine seaice_fme_derived_fields_read_restart!}}}

!***********************************************************************

   subroutine seaice_fme_derived_fields_write_restart(domain, err)!{{{
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND) :: ts_str
      integer :: ierr_local

      err = 0
      if (.not. time_avg_enabled) return
      if (.not. allocated(remap_accum)) return

      call mpas_get_time(avg_last_output_time, &
           dateTimeString=ts_str, ierr=ierr_local)
      call seaice_fme_horiz_remap_write_accum('fmeSeaiceDerivedFields', &
           remap_accum, N_AVG_FIELDS, nAccum, trim(ts_str), &
           count_arr=remap_accum_count, err=err)

   end subroutine seaice_fme_derived_fields_write_restart!}}}

!***********************************************************************
!
!  routine seaice_precompute_fme_derived_fields
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  FME Development Team
!> \date    April 2026
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_fme_derived_fields(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_precompute_fme_derived_fields!}}}

!***********************************************************************
!
!  routine seaice_compute_fme_derived_fields
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  FME Development Team
!> \date    April 2026
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member. It aggregates category-level sea ice
!>  tracers into cell-level derived fields for FME output.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_fme_derived_fields(domain, instance, timeLevel, err)!{{{

      use seaice_mesh, only: seaice_latlon_vector_rotation_backward

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (block_type), pointer :: &
           block

      type (mpas_pool_type), pointer :: &
           tracersPool, &
           meshPool, &
           velocitySolverPool, &
           amPool

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceAreaCategory, &
           iceVolumeCategory, &
           snowVolumeCategory, &
           surfaceTemperature

      real(kind=RKIND), dimension(:), pointer :: &
           iceAreaTotal, &
           iceVolumeTotal, &
           snowVolumeTotal, &
           iceThicknessMean, &
           surfaceTemperatureMean

      ! air stress rotation variables
      real(kind=RKIND), dimension(:), pointer :: &
           airStressCellU, &
           airStressCellV, &
           airStressZonal, &
           airStressMeridional, &
           latCell, &
           lonCell, &
           xCell, &
           yCell, &
           zCell

      ! ice velocity (cell-centered) and geographic-frame outputs
      real(kind=RKIND), dimension(:), pointer :: &
           uVelocityCell, &
           vVelocityCell, &
           uVelocityGeoCell, &
           vVelocityGeoCell

      real(kind=RKIND), pointer :: &
           sphereRadius

      logical, pointer :: &
           rotateCartesianGrid

      integer, pointer :: &
           nCellsSolve, &
           nCategories

      integer :: &
           iCell, &
           iCategory

      ! Remap workspace
      real(r8), allocatable :: fld_remap(:)
      integer :: ierr_remap
      type(MPAS_Time_Type) :: currTime

      real(kind=RKIND) :: &
           areaSum, &
           tempSum

      ! Fill value for intensive quantities at ice-free cells.
      ! Matches _FillValue attribute in remapped output files.
      real(kind=RKIND), parameter :: fillValue = real(SHR_FILL_VALUE, RKIND)

      err = 0

      ! get configs
      call MPAS_pool_get_config(domain % configs, 'config_rotate_cartesian_grid', rotateCartesianGrid)

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity_solver', velocitySolverPool)
         call mpas_pool_get_subpool(block % structs, 'fmeSeaiceDerivedFieldsAM', amPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCategories', nCategories)

         call mpas_pool_get_array(tracersPool, 'iceAreaCategory', iceAreaCategory, timeLevel)
         call mpas_pool_get_array(tracersPool, 'iceVolumeCategory', iceVolumeCategory, timeLevel)
         call mpas_pool_get_array(tracersPool, 'snowVolumeCategory', snowVolumeCategory, timeLevel)
         call mpas_pool_get_array(tracersPool, 'surfaceTemperature', surfaceTemperature, timeLevel)

         call mpas_pool_get_array(amPool, 'iceAreaTotal', iceAreaTotal)
         call mpas_pool_get_array(amPool, 'iceVolumeTotal', iceVolumeTotal)
         call mpas_pool_get_array(amPool, 'snowVolumeTotal', snowVolumeTotal)
         call mpas_pool_get_array(amPool, 'iceThicknessMean', iceThicknessMean)
         call mpas_pool_get_array(amPool, 'surfaceTemperatureMean', surfaceTemperatureMean)

         ! air stress rotation arrays
         call MPAS_pool_get_config(meshPool, 'sphere_radius', sphereRadius)
         call MPAS_pool_get_array(meshPool, 'latCell', latCell)
         call MPAS_pool_get_array(meshPool, 'lonCell', lonCell)
         call MPAS_pool_get_array(meshPool, 'xCell', xCell)
         call MPAS_pool_get_array(meshPool, 'yCell', yCell)
         call MPAS_pool_get_array(meshPool, 'zCell', zCell)

         call MPAS_pool_get_array(velocitySolverPool, 'airStressCellU', airStressCellU)
         call MPAS_pool_get_array(velocitySolverPool, 'airStressCellV', airStressCellV)
         call MPAS_pool_get_array(velocitySolverPool, 'uVelocityCell', uVelocityCell)
         call MPAS_pool_get_array(velocitySolverPool, 'vVelocityCell', vVelocityCell)

         call MPAS_pool_get_array(amPool, 'airStressZonal', airStressZonal)
         call MPAS_pool_get_array(amPool, 'airStressMeridional', airStressMeridional)
         call MPAS_pool_get_array(amPool, 'uVelocityGeoCell', uVelocityGeoCell)
         call MPAS_pool_get_array(amPool, 'vVelocityGeoCell', vVelocityGeoCell)

         do iCell = 1, nCellsSolve

            ! Sum over categories
            areaSum = 0.0_RKIND
            iceVolumeTotal(iCell) = 0.0_RKIND
            snowVolumeTotal(iCell) = 0.0_RKIND
            tempSum = 0.0_RKIND

            do iCategory = 1, nCategories
               areaSum = areaSum + iceAreaCategory(1, iCategory, iCell)
               iceVolumeTotal(iCell) = iceVolumeTotal(iCell) + iceVolumeCategory(1, iCategory, iCell)
               snowVolumeTotal(iCell) = snowVolumeTotal(iCell) + snowVolumeCategory(1, iCategory, iCell)
               tempSum = tempSum + iceAreaCategory(1, iCategory, iCell) * surfaceTemperature(1, iCategory, iCell)
            end do

            iceAreaTotal(iCell) = areaSum

            ! Area-weighted mean ice thickness and surface temperature.
            ! Intensive quantities are undefined for ice-free cells and must
            ! be set to fillValue so that the horizontal remap (apply_masked)
            ! correctly excludes them from weighted averages at the ice edge.
            if (areaSum > 0.0_RKIND) then
               iceThicknessMean(iCell) = iceVolumeTotal(iCell) / areaSum
               surfaceTemperatureMean(iCell) = tempSum / areaSum
            else
               iceThicknessMean(iCell) = fillValue
               surfaceTemperatureMean(iCell) = fillValue
            end if

            ! Rotate air stress from model coordinates to geographical coordinates.
            ! For ice-free cells, set to fillValue (stress on non-existent ice
            ! is undefined, not zero).
            if (areaSum > 0.0_RKIND) then
               call seaice_latlon_vector_rotation_backward(&
                    airStressZonal(iCell), &
                    airStressMeridional(iCell), &
                    airStressCellU(iCell), &
                    airStressCellV(iCell), &
                    latCell(iCell), &
                    lonCell(iCell), &
                    xCell(iCell), &
                    yCell(iCell), &
                    zCell(iCell), &
                    sphereRadius, &
                    rotateCartesianGrid)
            else
               airStressZonal(iCell) = fillValue
               airStressMeridional(iCell) = fillValue
            end if

            ! Rotate cell-centered ice velocity to geographic frame.
            ! Velocity on ice-free cells is undefined -> fillValue.
            if (areaSum > 0.0_RKIND) then
               call seaice_latlon_vector_rotation_backward(&
                    uVelocityGeoCell(iCell), &
                    vVelocityGeoCell(iCell), &
                    uVelocityCell(iCell), &
                    vVelocityCell(iCell), &
                    latCell(iCell), &
                    lonCell(iCell), &
                    xCell(iCell), &
                    yCell(iCell), &
                    zCell(iCell), &
                    sphereRadius, &
                    rotateCartesianGrid)
            else
               uVelocityGeoCell(iCell) = fillValue
               vVelocityGeoCell(iCell) = fillValue
            end if

         end do ! iCell

         ! --- Horizontal remapping to lat-lon grid ---
         if (do_horiz_remap) then
            allocate(fld_remap(remap_nb_local))

            if (time_avg_enabled) then
               ! Accumulate path: each remap adds to remap_accum.
               call remap_1field(iceAreaTotal, 1)
               call remap_1field(iceVolumeTotal, 2)
               call remap_1field(snowVolumeTotal, 3)
               call remap_1field(iceThicknessMean, 4)
               call remap_1field(surfaceTemperatureMean, 5)
               call remap_1field(airStressZonal, 6)
               call remap_1field(airStressMeridional, 7)
               call remap_1field(uVelocityGeoCell, 8)
               call remap_1field(vVelocityGeoCell, 9)
               nAccum = nAccum + 1

               ! Check if output interval elapsed
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr_remap)
               if (currTime - avg_last_output_time >= avg_output_interval) then
                  call seaice_fme_horiz_remap_check_rotate(domain, ierr_remap)
                  if (file_is_open) then
                     ! time_bnds = [avg window start, avg window end]
                     call seaice_fme_horiz_remap_write_time( &
                          avg_last_output_time, currTime, ierr_remap)
                     call write_avg_field('iceAreaTotal', 1)
                     call write_avg_field('iceVolumeTotal', 2)
                     call write_avg_field('snowVolumeTotal', 3)
                     call write_avg_field('iceThicknessMean', 4)
                     call write_avg_field('surfaceTemperatureMean', 5)
                     call write_avg_field('airStressZonal', 6)
                     call write_avg_field('airStressMeridional', 7)
                     call write_avg_field('uVelocityGeoCell', 8)
                     call write_avg_field('vVelocityGeoCell', 9)
                  end if
                  remap_accum = 0.0_r8
                  remap_accum_count = 0.0_r8
                  nAccum = 0
                  ! Advance by exact interval to avoid drift.
                  avg_last_output_time = avg_last_output_time + avg_output_interval
               end if
            else
               ! Instantaneous path: remap and write each field one at a time
               ! so fld_remap holds the correct field at each write.
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr_remap)
               call seaice_fme_horiz_remap_check_rotate(domain, ierr_remap)
               if (file_is_open) then
                  ! Instantaneous: time_bnds collapses to [t, t]
                  call seaice_fme_horiz_remap_write_time(currTime, currTime, &
                       ierr_remap)
                  call remap_and_write_inst('iceAreaTotal', iceAreaTotal)
                  call remap_and_write_inst('iceVolumeTotal', iceVolumeTotal)
                  call remap_and_write_inst('snowVolumeTotal', snowVolumeTotal)
                  call remap_and_write_inst('iceThicknessMean', iceThicknessMean)
                  call remap_and_write_inst('surfaceTemperatureMean', surfaceTemperatureMean)
                  call remap_and_write_inst('airStressZonal', airStressZonal)
                  call remap_and_write_inst('airStressMeridional', airStressMeridional)
                  call remap_and_write_inst('uVelocityGeoCell', uVelocityGeoCell)
                  call remap_and_write_inst('vVelocityGeoCell', vVelocityGeoCell)
               end if
            end if

            deallocate(fld_remap)
         end if

         block => block % next
      end do

   contains

      subroutine remap_1field(native_fld, idx)
         ! Per-cell accumulation that excludes fill values from the running
         ! sum. Intensive sea-ice fields (iceThicknessMean, surface temp,
         ! airStress*, u/vVelocityGeoCell) are SHR_FILL_VALUE on ice-free
         ! cells, and ice presence flips at the ice edge across the
         ! averaging window. Without per-cell masking, fill samples sum
         ! into remap_accum and surface as fractional-of-fill pollution
         ! after dividing by nAccum. See AGENTS.md gotcha on time-averaging
         ! fill pollution.
         real(kind=RKIND), dimension(:), intent(in) :: native_fld
         integer, intent(in) :: idx
         integer :: ib
         real(r8), parameter :: fill_thresh = real(SHR_FILL_VALUE, r8) * 0.1_r8
         call seaice_fme_horiz_remap_field(domain, native_fld, &
              nCellsSolve, fld_remap, ierr_remap)
         if (time_avg_enabled) then
            do ib = 1, remap_nb_local
               if (abs(fld_remap(ib)) < fill_thresh) then
                  remap_accum(ib, idx)       = remap_accum(ib, idx)       + fld_remap(ib)
                  remap_accum_count(ib, idx) = remap_accum_count(ib, idx) + 1.0_r8
               end if
            end do
         end if
      end subroutine remap_1field

      subroutine write_avg_field(varname, idx)
         ! Divide by per-cell count rather than the global nAccum so cells
         ! that were sometimes valid and sometimes fill within the window
         ! get the unbiased average over their valid samples. Cells that
         ! were fill on every sample (count == 0) emit SHR_FILL_VALUE.
         character(len=*), intent(in) :: varname
         integer, intent(in) :: idx
         integer :: ib
         do ib = 1, remap_nb_local
            if (remap_accum_count(ib, idx) > 0.0_r8) then
               fld_remap(ib) = remap_accum(ib, idx) / remap_accum_count(ib, idx)
            else
               fld_remap(ib) = real(SHR_FILL_VALUE, r8)
            end if
         end do
         call seaice_fme_horiz_remap_write_var(varname, fld_remap, 1, ierr_remap)
      end subroutine write_avg_field

      subroutine remap_and_write_inst(varname, native_fld)
         ! Remap one field to fld_remap, then write under varname.
         ! Used in instantaneous mode to ensure each varname gets the
         ! correct field's data (the previous bulk-remap-then-bulk-write
         ! pattern silently wrote the last remapped field under every name).
         character(len=*), intent(in) :: varname
         real(kind=RKIND), dimension(:), intent(in) :: native_fld
         call seaice_fme_horiz_remap_field(domain, native_fld, &
              nCellsSolve, fld_remap, ierr_remap)
         call seaice_fme_horiz_remap_write_var(varname, fld_remap, 1, ierr_remap)
      end subroutine remap_and_write_inst

   end subroutine seaice_compute_fme_derived_fields!}}}

!***********************************************************************
!
!  routine seaice_restart_fme_derived_fields
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  FME Development Team
!> \date    April 2026
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_fme_derived_fields(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_restart_fme_derived_fields!}}}

!***********************************************************************
!
!  routine seaice_finalize_fme_derived_fields
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  FME Development Team
!> \date    April 2026
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_fme_derived_fields(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: ierr

      err = 0

      call seaice_fme_horiz_remap_close_file(ierr)
      if (allocated(remap_accum)) deallocate(remap_accum)
      if (allocated(remap_accum_count)) deallocate(remap_accum_count)

   end subroutine seaice_finalize_fme_derived_fields!}}}

end module seaice_fme_derived_fields

! vim: foldmethod=marker
