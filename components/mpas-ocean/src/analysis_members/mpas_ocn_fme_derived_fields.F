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

!***********************************************************************

contains

!***********************************************************************

   subroutine ocn_init_fme_derived_fields(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      character(len=StrKIND), pointer :: config_AM_fmeDerivedFields_horiz_remap_file
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

      ! Fill value for land cells (must match _FillValue in remapped output)
      real(kind=RKIND), parameter :: fillValue = 1.0E34_RKIND
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
               surfaceHeatFluxTotal(iCell) = shortWaveHeatFlux(iCell) &
                    + longWaveHeatFluxDown(iCell) &
                    + latentHeatFlux(iCell) &
                    + sensibleHeatFlux(iCell)
            else
               surfaceHeatFluxTotal(iCell) = fillValue
            end if
         end do

         ! --- Horizontal remapping to lat-lon grid ---
         if (do_horiz_remap) then
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            ! Open or rotate file based on current date
            call remap_outfile%check_rotate(domain, ierr)

            if (remap_outfile%file_is_open) then
               allocate(fld_remap(remap_nb_local, 1))
               allocate(field_1lev(1, nCells))

               call remap_outfile%write_time('', ierr)

               field_1lev(1, 1:nCells) = sst(1:nCells)
               call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                    mask=maxLevelCell)
               call remap_outfile%write_var('sst', fld_remap, 1, ierr)

               field_1lev(1, 1:nCells) = sss(1:nCells)
               call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                    mask=maxLevelCell)
               call remap_outfile%write_var('sss', fld_remap, 1, ierr)

               field_1lev(1, 1:nCells) = surfaceHeatFluxTotal(1:nCells)
               call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                    mask=maxLevelCell)
               call remap_outfile%write_var('surfaceHeatFluxTotal', fld_remap, 1, ierr)

               ! ssh (from state pool)
               call mpas_pool_get_array(statePool, 'ssh', ssh_ptr, 1)
               if (associated(ssh_ptr)) then
                  field_1lev(1, 1:nCells) = ssh_ptr(1:nCells)
                  call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                       mask=maxLevelCell)
                  call remap_outfile%write_var('ssh', fld_remap, 1, ierr)
               end if

               ! Flux fields (already fetched from forcingPool)
               if (associated(shortWaveHeatFlux)) then
                  field_1lev(1, 1:nCells) = shortWaveHeatFlux(1:nCells)
                  call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                       mask=maxLevelCell)
                  call remap_outfile%write_var('shortWaveHeatFlux', fld_remap, 1, ierr)
               end if

               if (associated(longWaveHeatFluxDown)) then
                  field_1lev(1, 1:nCells) = longWaveHeatFluxDown(1:nCells)
                  call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                       mask=maxLevelCell)
                  call remap_outfile%write_var('longWaveHeatFluxDown', fld_remap, 1, ierr)
               end if

               if (associated(latentHeatFlux)) then
                  field_1lev(1, 1:nCells) = latentHeatFlux(1:nCells)
                  call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                       mask=maxLevelCell)
                  call remap_outfile%write_var('latentHeatFlux', fld_remap, 1, ierr)
               end if

               if (associated(sensibleHeatFlux)) then
                  field_1lev(1, 1:nCells) = sensibleHeatFlux(1:nCells)
                  call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                       mask=maxLevelCell)
                  call remap_outfile%write_var('sensibleHeatFlux', fld_remap, 1, ierr)
               end if

               ! Wind stress (from forcing pool)
               if (associated(windStressZonalPtr)) then
                  field_1lev(1, 1:nCells) = windStressZonalPtr(1:nCells)
                  call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                       mask=maxLevelCell)
                  call remap_outfile%write_var('windStressZonal', fld_remap, 1, ierr)
               end if

               if (associated(windStressMeridionalPtr)) then
                  field_1lev(1, 1:nCells) = windStressMeridionalPtr(1:nCells)
                  call ocn_fme_horiz_remap_field(domain, field_1lev, 1, fld_remap, ierr, &
                       mask=maxLevelCell)
                  call remap_outfile%write_var('windStressMeridional', fld_remap, 1, ierr)
               end if

               deallocate(fld_remap, field_1lev)
            end if
         end if

         block => block % next
      end do

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
