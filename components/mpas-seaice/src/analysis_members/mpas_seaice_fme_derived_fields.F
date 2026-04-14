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
        seaice_fme_horiz_remap_field, seaice_fme_horiz_remap_active, &
        seaice_fme_horiz_remap_write_var, seaice_fme_horiz_remap_write_time, &
        seaice_fme_horiz_remap_close_file, seaice_fme_horiz_remap_get_dims, &
        seaice_fme_horiz_remap_register_var, seaice_fme_horiz_remap_check_rotate, &
        remap_stream_name, file_is_open

   implicit none
   private
   save

   ! Remap state
   logical :: do_horiz_remap = .false.
   integer :: remap_nlat = 0, remap_nlon = 0, remap_nb_local = 0

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
             seaice_finalize_fme_derived_fields

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
               end if
            end if
         end if
      end if

   end subroutine seaice_init_fme_derived_fields!}}}

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

         call MPAS_pool_get_array(amPool, 'airStressZonal', airStressZonal)
         call MPAS_pool_get_array(amPool, 'airStressMeridional', airStressMeridional)

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

         end do ! iCell

         ! --- Horizontal remapping to lat-lon grid ---
         if (do_horiz_remap) then
            call seaice_fme_horiz_remap_check_rotate(domain, ierr_remap)

            if (file_is_open) then
               allocate(fld_remap(remap_nb_local))

               call seaice_fme_horiz_remap_write_time('', ierr_remap)

               call seaice_fme_horiz_remap_field(domain, iceAreaTotal, &
                    nCellsSolve, fld_remap, ierr_remap)
               call seaice_fme_horiz_remap_write_var('iceAreaTotal', fld_remap, 1, ierr_remap)

               call seaice_fme_horiz_remap_field(domain, iceVolumeTotal, &
                    nCellsSolve, fld_remap, ierr_remap)
               call seaice_fme_horiz_remap_write_var('iceVolumeTotal', fld_remap, 1, ierr_remap)

               call seaice_fme_horiz_remap_field(domain, snowVolumeTotal, &
                    nCellsSolve, fld_remap, ierr_remap)
               call seaice_fme_horiz_remap_write_var('snowVolumeTotal', fld_remap, 1, ierr_remap)

               call seaice_fme_horiz_remap_field(domain, iceThicknessMean, &
                    nCellsSolve, fld_remap, ierr_remap)
               call seaice_fme_horiz_remap_write_var('iceThicknessMean', fld_remap, 1, ierr_remap)

               call seaice_fme_horiz_remap_field(domain, surfaceTemperatureMean, &
                    nCellsSolve, fld_remap, ierr_remap)
               call seaice_fme_horiz_remap_write_var('surfaceTemperatureMean', fld_remap, 1, ierr_remap)

               call seaice_fme_horiz_remap_field(domain, airStressZonal, &
                    nCellsSolve, fld_remap, ierr_remap)
               call seaice_fme_horiz_remap_write_var('airStressZonal', fld_remap, 1, ierr_remap)

               call seaice_fme_horiz_remap_field(domain, airStressMeridional, &
                    nCellsSolve, fld_remap, ierr_remap)
               call seaice_fme_horiz_remap_write_var('airStressMeridional', fld_remap, 1, ierr_remap)

               deallocate(fld_remap)
            end if
         end if

         block => block % next
      end do

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

   end subroutine seaice_finalize_fme_derived_fields!}}}

end module seaice_fme_derived_fields

! vim: foldmethod=marker
