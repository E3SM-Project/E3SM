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
        ocn_fme_horiz_remap_field, ocn_fme_horiz_remap_active, &
        ocn_fme_horiz_remap_get_dims, ocn_fme_remap_file_t, &
        ocn_fme_remap_build_filename

   implicit none
   private
   save

   public :: ocn_init_fme_depth_coarsening, &
             ocn_compute_fme_depth_coarsening, &
             ocn_restart_fme_depth_coarsening, &
             ocn_finalize_fme_depth_coarsening

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
   real(r8), allocatable :: remap_accum(:,:)  ! (remap_nb_local, n_avg_fields)

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
      character(len=StrKIND), pointer :: config_stream_output_interval
      character(len=StrKIND) :: tempStr
      character(len=128) :: varName, longName
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
            do_horiz_remap = ocn_fme_horiz_remap_active()
            if (do_horiz_remap) then
               call ocn_fme_horiz_remap_get_dims(remap_nlat, remap_nlon, remap_nb_local)
               call mpas_log_write('fmeDepthCoarsening: horiz remap active, ' &
                    // 'target grid $i x $i', &
                    intArgs=(/remap_nlon, remap_nlat/))
               remap_outfile%stream_name = 'fmeDepthCoarseningOutput'
               ! Register split-level 2D vars for each field and depth level
               do iField = 1, nCoarsenFields
                  do kOut = 0, nCoarsenLevels - 1
                     write(varName, '(A,A,I0)') &
                          trim(coarsenFieldNames(iField)), 'Coarsened_', kOut
                     write(longName, '(A,I0)') 'Depth-coarsened ' &
                          // trim(coarsenFieldNames(iField)) // ' level ', kOut
                     call remap_outfile%register_var(trim(varName), '', &
                          trim(longName), ierr)
                  end do
               end do
               do kOut = 0, nCoarsenLevels - 1
                  write(varName, '(A,I0)') 'layerThicknessCoarsened_', kOut
                  write(longName, '(A,I0)') &
                       'Depth-coarsened layer thickness level ', kOut
                  call remap_outfile%register_var(trim(varName), 'm', &
                       trim(longName), ierr)
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
                     avg_last_output_time = mpas_get_clock_time( &
                          domain % clock, MPAS_START_TIME, ierr)
                     n_avg_fields = nCoarsenLevels * nCoarsenFields &
                          + nCoarsenLevels
                     allocate(remap_accum(remap_nb_local, n_avg_fields))
                     remap_accum = 0.0_r8
                     nAccum = 0
                     time_avg_enabled = .true.
                     call mpas_log_write( &
                          'fmeDepthCoarsening: time averaging enabled')
                  end if
               end if
            end if
         end if
      end if

   end subroutine ocn_init_fme_depth_coarsening!}}}

!***********************************************************************

   subroutine ocn_compute_fme_depth_coarsening(domain, timeLevel, err)!{{{

      integer, intent(in) :: timeLevel
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: tracersPool
      type(mpas_pool_type), pointer :: amPool

      integer, pointer :: nCellsSolve, nVertLevels
      integer, pointer :: index_temperature, index_salinity
      integer, pointer :: nFmeDepthLevels
      integer :: iCell, k, iField
      integer, dimension(:), pointer :: maxLevelCell

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers
      real(kind=RKIND), dimension(:, :), pointer :: layerThickness
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
      character(len=StrKIND) :: timeStr
      character(len=128) :: varName, longName
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
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
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

            ! Build depth interface array from layer thicknesses
            coord_iface(1) = 0.0_r8
            do k = 1, maxLevelCell(iCell)
               coord_iface(k + 1) = coord_iface(k) + layerThickness(k, iCell)
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
               call accumulate_thickness(depthTop, layerThickness(k, iCell), &
                    nCoarsenLevels, layerThicknessCoarsened(:, iCell))
               depthTop = depthTop + layerThickness(k, iCell)
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

            ! Remap each coarsened field, one depth level at a time
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

            ! Remap layerThicknessCoarsened split levels
            do kOut = 0, nCoarsenLevels - 1
               call remap_1level(layerThicknessCoarsened, kOut, &
                    nCoarsenFields * nCoarsenLevels + kOut + 1)
            end do

            if (time_avg_enabled) then
               ! Accumulate (remapped fields already added in remap_1level)
               nAccum = nAccum + 1

               ! Check if output interval elapsed
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
               if (currTime - avg_last_output_time >= avg_output_interval) then
                  call mpas_get_time(currTime, dateTimeString=timeStr)
                  call remap_outfile%check_rotate(domain, ierr)
                  if (remap_outfile%file_is_open) then
                     call remap_outfile%write_time(trim(timeStr), ierr)

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
                  end if
                  remap_accum = 0.0_r8
                  nAccum = 0
                  avg_last_output_time = currTime
               end if
            else
               ! Instantaneous: write immediately
               call remap_outfile%check_rotate(domain, ierr)
               if (remap_outfile%file_is_open) then
                  currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
                  call mpas_get_time(currTime, dateTimeString=timeStr)
                  call remap_outfile%write_time(trim(timeStr), ierr)

                  do ifld = 1, nCoarsenFields
                     do kOut = 0, nCoarsenLevels - 1
                        write(varName, '(A,A,I0)') &
                             trim(coarsenFieldNames(ifld)), &
                             'Coarsened_', kOut
                        call write_inst_field(trim(varName))
                     end do
                  end do
                  do kOut = 0, nCoarsenLevels - 1
                     write(varName, '(A,I0)') &
                          'layerThicknessCoarsened_', kOut
                     call write_inst_field(trim(varName))
                  end do
               end if
            end if

            deallocate(fld_remap, field_1lev, levelMask)
         end if

         block => block % next
      end do

   contains

      subroutine remap_1level(srcArray, kLevel, idx)
         real(kind=RKIND), dimension(:,:), intent(in) :: srcArray
         integer, intent(in) :: kLevel, idx
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
            remap_accum(:, idx) = remap_accum(:, idx) + fld_remap(:, 1)
         end if
      end subroutine remap_1level

      subroutine write_avg_field(vname, idx)
         character(len=*), intent(in) :: vname
         integer, intent(in) :: idx
         fld_remap(:, 1) = remap_accum(:, idx) / real(nAccum, r8)
         call remap_outfile%write_var(vname, fld_remap, 1, ierr)
      end subroutine write_avg_field

      subroutine write_inst_field(vname)
         character(len=*), intent(in) :: vname
         call remap_outfile%write_var(vname, fld_remap, 1, ierr)
      end subroutine write_inst_field

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

   end subroutine ocn_finalize_fme_depth_coarsening!}}}

end module ocn_fme_depth_coarsening

! vim: foldmethod=marker
