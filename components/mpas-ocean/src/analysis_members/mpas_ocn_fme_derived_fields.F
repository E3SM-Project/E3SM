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
   integer, parameter :: N_AVG_FIELDS = 10
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
                        avg_last_output_time = mpas_get_clock_time( &
                             domain % clock, MPAS_START_TIME, ierr)
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
      type(mpas_pool_type), pointer :: amPool

      integer, pointer :: nCellsSolve, nCells
      integer, pointer :: index_temperature, index_salinity
      integer, dimension(:), pointer :: maxLevelCell
      integer :: iCell, ierr

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers
      real(kind=RKIND), dimension(:), pointer :: sst, sss, surfaceHeatFluxTotal
      real(kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux
      real(kind=RKIND), dimension(:), pointer :: longWaveHeatFluxDown
      real(kind=RKIND), dimension(:), pointer :: latentHeatFlux
      real(kind=RKIND), dimension(:), pointer :: sensibleHeatFlux
      real(kind=RKIND), dimension(:), pointer :: ssh_ptr
      real(kind=RKIND), dimension(:), pointer :: windStressZonalPtr
      real(kind=RKIND), dimension(:), pointer :: windStressMeridionalPtr
      real(kind=RKIND), dimension(:), pointer :: iceFractionPtr

      ! Fill value for land/ice-masked cells
      real(kind=RKIND), parameter :: fillValue = real(SHR_FILL_VALUE, RKIND)
      ! Threshold: mask coupler fluxes where ice covers > 99% of cell
      real(kind=RKIND), parameter :: iceMaskThresh = 0.99_RKIND
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
         call mpas_pool_get_subpool(block % structs, 'fmeDerivedFieldsAM', amPool)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(amPool, 'sst', sst)
         call mpas_pool_get_array(amPool, 'sss', sss)
         call mpas_pool_get_array(amPool, 'surfaceHeatFluxTotal', surfaceHeatFluxTotal)

         call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
         call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
         call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
         call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
         call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonalPtr)
         call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridionalPtr)
         call mpas_pool_get_array(forcingPool, 'iceFraction', iceFractionPtr)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         do iCell = 1, nCellsSolve
            ! Land cells (maxLevelCell==0) get fillValue so the remap
            ! correctly excludes them from weighted averages.
            if (maxLevelCell(iCell) == 0) then
               sst(iCell) = fillValue
               sss(iCell) = fillValue
               surfaceHeatFluxTotal(iCell) = fillValue
               cycle
            end if

            sst(iCell) = tracers(index_temperature, 1, iCell)
            sss(iCell) = tracers(index_salinity, 1, iCell)

            if (associated(shortWaveHeatFlux) .and. &
                associated(longWaveHeatFluxDown) .and. &
                associated(latentHeatFlux) .and. &
                associated(sensibleHeatFlux)) then
               ! Sum of coupler fluxes; undo ice-fraction weighting
               ! to recover atmospheric total (same correction as
               ! mask_and_remap applies to individual components).
               surfaceHeatFluxTotal(iCell) = shortWaveHeatFlux(iCell) &
                    + longWaveHeatFluxDown(iCell) &
                    + latentHeatFlux(iCell) &
                    + sensibleHeatFlux(iCell)
               if (associated(iceFractionPtr)) then
                  if (iceFractionPtr(iCell) >= iceMaskThresh) then
                     surfaceHeatFluxTotal(iCell) = fillValue
                  else if (iceFractionPtr(iCell) > 0.0_RKIND) then
                     surfaceHeatFluxTotal(iCell) = surfaceHeatFluxTotal(iCell) &
                          / (1.0_RKIND - iceFractionPtr(iCell))
                  end if
               end if
            else
               surfaceHeatFluxTotal(iCell) = fillValue
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
            call mpas_pool_get_array(statePool, 'ssh', ssh_ptr, 1)

            call remap_1field(sst, 1)
            call remap_1field(sss, 2)
            call remap_1field(surfaceHeatFluxTotal, 3)

            ! SSH: state variable, no ice-fraction correction needed
            if (associated(ssh_ptr)) call remap_1field(ssh_ptr, 4)

            ! Coupler fluxes: undo ice-fraction weighting to recover
            ! the coupler delivers near-zero flux to the ocean surface.
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

            if (time_avg_enabled) then
               ! Accumulate (remapped fields already added in remap_1field)
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
                  end if
                  remap_accum = 0.0_r8
                  nAccum = 0
                  avg_last_output_time = currTime
               end if
            else
               ! Instantaneous: write immediately
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
               call mpas_get_time(currTime, dateTimeString=timeStr)
               call remap_outfile%check_rotate(domain, ierr)
               if (remap_outfile%file_is_open) then
                  call remap_outfile%write_time(trim(timeStr), ierr)
                  call write_inst_field('sst', 1)
                  call write_inst_field('sss', 2)
                  call write_inst_field('surfaceHeatFluxTotal', 3)
                  if (associated(ssh_ptr)) call write_inst_field('ssh', 4)
                  if (associated(shortWaveHeatFlux)) call write_inst_field('shortWaveHeatFlux', 5)
                  if (associated(longWaveHeatFluxDown)) call write_inst_field('longWaveHeatFluxDown', 6)
                  if (associated(latentHeatFlux)) call write_inst_field('latentHeatFlux', 7)
                  if (associated(sensibleHeatFlux)) call write_inst_field('sensibleHeatFlux', 8)
                  if (associated(windStressZonalPtr)) call write_inst_field('windStressZonal', 9)
                  if (associated(windStressMeridionalPtr)) call write_inst_field('windStressMeridional', 10)
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
         ! Copy field, undo coupler ice-fraction weighting to recover the
         ! atmospheric flux, then remap.  Coupler delivers:
         !   flux_ocean = flux_atm * (1 - iceFraction)
         ! We recover flux_atm = flux_ocean / (1 - iceFraction).
         ! At nearly full ice coverage (>= iceMaskThresh), set fillValue
         ! to avoid division by near-zero.
         ! Land cells (maxLevelCell==0) also get fillValue.
         real(kind=RKIND), dimension(:), intent(in) :: native_fld
         integer, intent(in) :: idx
         integer :: ic
         real(kind=RKIND) :: openFrac
         fld_masked(1:nCells) = native_fld(1:nCells)
         do ic = 1, nCells
            if (maxLevelCell(ic) == 0) then
               fld_masked(ic) = fillValue
            else if (associated(iceFractionPtr)) then
               if (iceFractionPtr(ic) >= iceMaskThresh) then
                  fld_masked(ic) = fillValue
               else if (iceFractionPtr(ic) > 0.0_RKIND) then
                  openFrac = 1.0_RKIND - iceFractionPtr(ic)
                  fld_masked(ic) = native_fld(ic) / openFrac
               end if
            end if
         end do
         call remap_1field(fld_masked, idx)
      end subroutine mask_and_remap

      subroutine write_avg_field(varname, idx)
         character(len=*), intent(in) :: varname
         integer, intent(in) :: idx
         fld_remap(:, 1) = remap_accum(:, idx) / real(nAccum, r8)
         call remap_outfile%write_var(varname, fld_remap, 1, ierr)
      end subroutine write_avg_field

      subroutine write_inst_field(varname, idx)
         character(len=*), intent(in) :: varname
         integer, intent(in) :: idx
         call remap_outfile%write_var(varname, fld_remap, 1, ierr)
      end subroutine write_inst_field

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
