! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_fme_vertical_reduce
!
!> \brief MPAS ocean analysis member: FME vertical reduction
!> \author FME Training Data Support
!> \date   April 2026
!> \details
!>  Computes column-integrated (vertical reduction) quantities for FME
!>  training data. Built-in fields:
!>    - oceanHeatContent:  sum(T * h * rho_sw * cp_sw)  [J/m^2]
!>    - freshwaterContent: sum((S_ref - S)/S_ref * h)    [m]
!>    - kineticEnergy:     sum(0.5*(u^2+v^2) * h * rho_sw) [J/m^2]
!>
!>  Also supports user-configurable fields via namelist with sum or
!>  thickness-weighted average reduction methods.
!>
!-----------------------------------------------------------------------

module ocn_fme_vertical_reduce

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_kind_types

   use ocn_constants
   use ocn_diagnostics_variables
   use mpas_log, only: mpas_log_write

   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_horiz_remap_mod, only: SHR_FILL_VALUE
   use ocn_fme_horiz_remap, only: ocn_fme_horiz_remap_init, &
        ocn_fme_horiz_remap_field, ocn_fme_horiz_remap_active, &
        ocn_fme_horiz_remap_get_dims, ocn_fme_remap_file_t, &
        ocn_fme_remap_build_filename

   implicit none
   private
   save

   public :: ocn_init_fme_vertical_reduce, &
             ocn_compute_fme_vertical_reduce, &
             ocn_restart_fme_vertical_reduce, &
             ocn_finalize_fme_vertical_reduce

   real(kind=RKIND) :: referenceSalinity = 34.8_RKIND

   ! Remapped output workspace
   logical :: do_horiz_remap = .false.
   type(ocn_fme_remap_file_t) :: remap_outfile
   integer :: remap_nlat = 0, remap_nlon = 0, remap_nb_local = 0

   ! Time-averaging state for remapped output
   logical :: time_avg_enabled = .false.
   integer :: nAccum = 0
   type(MPAS_TimeInterval_type) :: avg_output_interval
   type(MPAS_Time_type) :: avg_last_output_time
   integer, parameter :: N_AVG_FIELDS = 8  ! 3 built-in + up to 5 generic
   real(r8), allocatable :: remap_accum(:,:)  ! (remap_nb_local, N_AVG_FIELDS)

   !> Number of user-configured generic reduction fields (parsed from config)
   integer :: nReduceFields = 0

   !> Field names parsed from config_AM_fmeVerticalReduce_fields
   character(len=64), allocatable :: reduceFieldNames(:)

   !> Reduction methods parsed from config_AM_fmeVerticalReduce_methods
   !> Each entry is 'sum' or 'avg'
   character(len=64), allocatable :: reduceFieldMethods(:)

!***********************************************************************

contains

!***********************************************************************

   subroutine ocn_init_fme_vertical_reduce(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      real(kind=RKIND), pointer :: config_AM_fmeVerticalReduce_reference_salinity
      character(len=StrKIND), pointer :: config_AM_fmeVerticalReduce_fields
      character(len=StrKIND), pointer :: config_AM_fmeVerticalReduce_methods
      character(len=StrKIND), pointer :: config_AM_fmeVerticalReduce_horiz_remap_file
      logical, pointer :: config_time_averaging
      character(len=StrKIND), pointer :: config_stream_output_interval
      character(len=StrKIND) :: tempStr, methodStr
      character(len=32) :: slotName
      integer :: nFields, nMethods, i, iStart, iEnd, strLen
      integer :: ierr, iGeneric

      err = 0

      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeVerticalReduce_reference_salinity', &
           config_AM_fmeVerticalReduce_reference_salinity)

      referenceSalinity = config_AM_fmeVerticalReduce_reference_salinity

      ! --- Parse generic field names and methods ---
      ! The built-in fields (OHC, FWC, KE) are always computed when
      ! the AM is enabled. The generic fields/methods config allows
      ! users to request additional reductions. However, MPAS Registry
      ! requires compile-time variable definitions, so truly dynamic
      ! output fields are not supported without Registry changes.
      ! We parse the config here and log a warning if fields are
      ! specified, noting that output requires corresponding Registry
      ! variable definitions.

      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeVerticalReduce_fields', &
           config_AM_fmeVerticalReduce_fields)
      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeVerticalReduce_methods', &
           config_AM_fmeVerticalReduce_methods)

      tempStr = trim(config_AM_fmeVerticalReduce_fields)
      methodStr = trim(config_AM_fmeVerticalReduce_methods)

      if (trim(tempStr) /= 'none' .and. len_trim(tempStr) > 0) then
         ! Count comma-separated field names
         nFields = 1
         strLen = len_trim(tempStr)
         do i = 1, strLen
            if (tempStr(i:i) == ',') nFields = nFields + 1
         end do

         ! Count comma-separated methods
         nMethods = 1
         strLen = len_trim(methodStr)
         do i = 1, strLen
            if (methodStr(i:i) == ',') nMethods = nMethods + 1
         end do

         if (nFields /= nMethods) then
            call mpas_log_write( &
                 'fmeVerticalReduce WARNING: number of fields does not match ' &
                 // 'number of methods. Generic fields will not be computed.')
            nReduceFields = 0
         else
            nReduceFields = nFields
            if (allocated(reduceFieldNames)) deallocate(reduceFieldNames)
            if (allocated(reduceFieldMethods)) deallocate(reduceFieldMethods)
            allocate(reduceFieldNames(nReduceFields))
            allocate(reduceFieldMethods(nReduceFields))

            ! Parse field names
            iStart = 1
            strLen = len_trim(tempStr)
            do i = 1, nReduceFields
               iEnd = index(tempStr(iStart:strLen), ',')
               if (iEnd == 0) then
                  iEnd = strLen + 1
               else
                  iEnd = iEnd + iStart - 1
               end if
               reduceFieldNames(i) = adjustl(tempStr(iStart:iEnd-1))
               iStart = iEnd + 1
            end do

            ! Parse methods
            iStart = 1
            strLen = len_trim(methodStr)
            do i = 1, nReduceFields
               iEnd = index(methodStr(iStart:strLen), ',')
               if (iEnd == 0) then
                  iEnd = strLen + 1
               else
                  iEnd = iEnd + iStart - 1
               end if
               reduceFieldMethods(i) = adjustl(methodStr(iStart:iEnd-1))
               iStart = iEnd + 1
            end do

            if (nReduceFields > 5) then
               call mpas_log_write( &
                    'fmeVerticalReduce WARNING: $i fields requested but only 5' &
                    // ' generic output slots available (vertReduceField01..05).' &
                    // ' Extra fields will be ignored.', &
                    intArgs=(/nReduceFields/))
            else
               call mpas_log_write( &
                    'fmeVerticalReduce: $i generic field reductions configured' &
                    // ' (output in vertReduceField01..$i)', &
                    intArgs=(/nReduceFields, nReduceFields/))
            end if
         end if
      else
         nReduceFields = 0
      end if

      ! --- Initialize horizontal remapping if requested ---
      call mpas_pool_get_config(domain % configs, &
           'config_AM_fmeVerticalReduce_horiz_remap_file', &
           config_AM_fmeVerticalReduce_horiz_remap_file)

      if (associated(config_AM_fmeVerticalReduce_horiz_remap_file)) then
         if (trim(config_AM_fmeVerticalReduce_horiz_remap_file) /= 'none' .and. &
             len_trim(config_AM_fmeVerticalReduce_horiz_remap_file) > 0) then

            call ocn_fme_horiz_remap_init(domain, &
                 trim(config_AM_fmeVerticalReduce_horiz_remap_file), ierr)
            if (ierr /= 0) then
               call mpas_log_write( &
                    'fmeVerticalReduce: horiz remap init failed, continuing without remap', &
                    messageType=MPAS_LOG_WARN)
               do_horiz_remap = .false.
            else
               do_horiz_remap = ocn_fme_horiz_remap_active()
               if (do_horiz_remap) then
                  call ocn_fme_horiz_remap_get_dims(remap_nlat, remap_nlon, remap_nb_local)
                  call mpas_log_write('fmeVerticalReduce: horiz remap active, ' &
                       // 'target grid $i x $i', &
                       intArgs=(/remap_nlon, remap_nlat/))
                  remap_outfile%stream_name = 'fmeVerticalReduceOutput'
                  call remap_outfile%register_var('oceanHeatContent', 'J m^-2', &
                       'Column-integrated ocean heat content', ierr)
                  call remap_outfile%register_var('freshwaterContent', 'm', &
                       'Column-integrated freshwater content', ierr)
                  call remap_outfile%register_var('kineticEnergy', 'J m^-2', &
                       'Column-integrated kinetic energy', ierr)
                  do iGeneric = 1, min(nReduceFields, 5)
                     write(slotName, '(A,I2.2)') 'vertReduceField', iGeneric
                     call remap_outfile%register_var(trim(slotName), 'varies', &
                          'User-configured vertical reduction field', ierr)
                  end do

                  ! Time-averaging setup
                  call mpas_pool_get_config(domain % configs, &
                       'config_AM_fmeVerticalReduce_time_averaging', &
                       config_time_averaging)
                  if (associated(config_time_averaging)) then
                     if (config_time_averaging) then
                        call mpas_pool_get_config(domain % configs, &
                             'config_AM_fmeVerticalReduce_stream_output_interval', &
                             config_stream_output_interval)
                        call mpas_set_timeInterval(avg_output_interval, &
                             timeString=trim(config_stream_output_interval), &
                             ierr=ierr)
                        avg_last_output_time = mpas_get_clock_time( &
                             domain % clock, MPAS_START_TIME, ierr)
                        allocate(remap_accum(remap_nb_local, N_AVG_FIELDS))
                        remap_accum = 0.0_r8
                        nAccum = 0
                        time_avg_enabled = .true.
                        call mpas_log_write( &
                             'fmeVerticalReduce: time averaging enabled')
                     end if
                  end if
               end if
            end if
         end if
      end if

   end subroutine ocn_init_fme_vertical_reduce!}}}

!***********************************************************************

   subroutine ocn_compute_fme_vertical_reduce(domain, timeLevel, err)!{{{

      integer, intent(in) :: timeLevel
      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: tracersPool
      type(mpas_pool_type), pointer :: amPool

      integer, pointer :: nCellsSolve, nCells, nVertLevels
      integer, pointer :: index_temperature, index_salinity
      integer :: iCell, k, iGeneric, ierr
      integer, dimension(:), pointer :: maxLevelCell

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers
      real(kind=RKIND), dimension(:, :), pointer :: layerThickness
      real(kind=RKIND), dimension(:), pointer :: oceanHeatContent
      real(kind=RKIND), dimension(:), pointer :: freshwaterContent
      real(kind=RKIND), dimension(:), pointer :: kineticEnergy

      ! Generic reduction output slots
      real(kind=RKIND), dimension(:), pointer :: genSlot
      character(len=32) :: slotName
      integer, pointer :: index_tracer_generic

      real(kind=RKIND) :: ohc_sum, fwc_sum, ke_sum, h_k
      real(kind=RKIND) :: gen_sum, gen_wgt
      real(kind=RKIND), parameter :: fillValue = real(SHR_FILL_VALUE, RKIND)

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
         call mpas_pool_get_subpool(block % structs, 'fmeVerticalReduceAM', amPool)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(amPool, 'oceanHeatContent', oceanHeatContent)
         call mpas_pool_get_array(amPool, 'freshwaterContent', freshwaterContent)
         call mpas_pool_get_array(amPool, 'kineticEnergy', kineticEnergy)

         do iCell = 1, nCellsSolve
            ! Land cells (maxLevelCell==0) get fillValue so the remap
            ! correctly excludes them from weighted averages.
            if (maxLevelCell(iCell) == 0) then
               oceanHeatContent(iCell) = fillValue
               freshwaterContent(iCell) = fillValue
               kineticEnergy(iCell) = fillValue
               cycle
            end if

            ohc_sum = 0.0_RKIND
            fwc_sum = 0.0_RKIND
            ke_sum  = 0.0_RKIND

            do k = 1, maxLevelCell(iCell)
               h_k = layerThickness(k, iCell)

               ! Ocean heat content: T * h * rho_sw * cp_sw
               ohc_sum = ohc_sum + tracers(index_temperature, k, iCell) &
                    * h_k * rho_sw * cp_sw

               ! Freshwater content: (S_ref - S) / S_ref * h
               fwc_sum = fwc_sum + (referenceSalinity &
                    - tracers(index_salinity, k, iCell)) &
                    / referenceSalinity * h_k

               ! Kinetic energy: 0.5 * (u^2 + v^2) * h * rho_sw
               if (associated(velocityZonal) .and. associated(velocityMeridional)) then
                  ke_sum = ke_sum + 0.5_RKIND &
                       * (velocityZonal(k, iCell)**2 &
                        + velocityMeridional(k, iCell)**2) &
                       * h_k * rho_sw
               end if
            end do

            oceanHeatContent(iCell) = ohc_sum
            freshwaterContent(iCell) = fwc_sum
            kineticEnergy(iCell) = ke_sum
         end do

         ! --- Generic field reductions using output slots ---
         ! Up to 5 user-configured fields written to vertReduceField01..05
         do iGeneric = 1, min(nReduceFields, 5)
            write(slotName, '(A,I2.2)') 'vertReduceField', iGeneric
            call mpas_pool_get_array(amPool, trim(slotName), genSlot)
            if (.not. associated(genSlot)) cycle

            ! Look up tracer by name
            nullify(index_tracer_generic)
            select case (trim(reduceFieldNames(iGeneric)))
            case ('temperature')
               call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_tracer_generic)
            case ('salinity')
               call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_tracer_generic)
            case default
               ! Unknown field — fill with fillValue
               do iCell = 1, nCellsSolve
                  genSlot(iCell) = fillValue
               end do
               cycle
            end select

            do iCell = 1, nCellsSolve
               if (maxLevelCell(iCell) == 0) then
                  genSlot(iCell) = fillValue
                  cycle
               end if

               gen_sum = 0.0_RKIND
               gen_wgt = 0.0_RKIND
               do k = 1, maxLevelCell(iCell)
                  h_k = layerThickness(k, iCell)
                  gen_sum = gen_sum + tracers(index_tracer_generic, k, iCell) * h_k
                  gen_wgt = gen_wgt + h_k
               end do

               if (trim(reduceFieldMethods(iGeneric)) == 'avg') then
                  if (gen_wgt > 0.0_RKIND) then
                     genSlot(iCell) = gen_sum / gen_wgt
                  else
                     genSlot(iCell) = fillValue
                  end if
               else  ! 'sum' is default
                  genSlot(iCell) = gen_sum
               end if
            end do
         end do

         ! --- Horizontal remapping to lat-lon grid ---
         if (do_horiz_remap) then
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

            allocate(fld_remap(remap_nb_local, 1))
            allocate(field_1lev(1, nCells))

            ! Remap each field to lat-lon
            call remap_1field(oceanHeatContent, 1)
            call remap_1field(freshwaterContent, 2)
            call remap_1field(kineticEnergy, 3)

            ! Remap generic reduction slots
            do iGeneric = 1, min(nReduceFields, 5)
               write(slotName, '(A,I2.2)') 'vertReduceField', iGeneric
               call mpas_pool_get_array(amPool, trim(slotName), genSlot)
               if (associated(genSlot)) then
                  call remap_1field(genSlot, 3 + iGeneric)
               end if
            end do

            if (time_avg_enabled) then
               ! Accumulate (remapped fields already added in remap_1field)
               nAccum = nAccum + 1

               ! Check if output interval elapsed
               currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
               if (currTime - avg_last_output_time >= avg_output_interval) then
                  call remap_outfile%check_rotate(domain, ierr)
                  if (remap_outfile%file_is_open) then
                     call remap_outfile%write_time('', ierr)
                     call write_avg_field('oceanHeatContent', 1)
                     call write_avg_field('freshwaterContent', 2)
                     call write_avg_field('kineticEnergy', 3)
                     do iGeneric = 1, min(nReduceFields, 5)
                        write(slotName, '(A,I2.2)') 'vertReduceField', iGeneric
                        call write_avg_field(trim(slotName), 3 + iGeneric)
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

                  call write_inst_field('oceanHeatContent')
                  call write_inst_field('freshwaterContent')
                  call write_inst_field('kineticEnergy')
                  do iGeneric = 1, min(nReduceFields, 5)
                     write(slotName, '(A,I2.2)') 'vertReduceField', iGeneric
                     call write_inst_field(trim(slotName))
                  end do
               end if
            end if

            deallocate(fld_remap, field_1lev)
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

      subroutine write_avg_field(varname, idx)
         character(len=*), intent(in) :: varname
         integer, intent(in) :: idx
         fld_remap(:, 1) = remap_accum(:, idx) / real(nAccum, r8)
         call remap_outfile%write_var(varname, fld_remap, 1, ierr)
      end subroutine write_avg_field

      subroutine write_inst_field(varname)
         character(len=*), intent(in) :: varname
         call remap_outfile%write_var(varname, fld_remap, 1, ierr)
      end subroutine write_inst_field

   end subroutine ocn_compute_fme_vertical_reduce!}}}

!***********************************************************************

   subroutine ocn_restart_fme_vertical_reduce(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      err = 0

   end subroutine ocn_restart_fme_vertical_reduce!}}}

!***********************************************************************

   subroutine ocn_finalize_fme_vertical_reduce(domain, err)!{{{

      type(domain_type), intent(inout) :: domain
      integer, intent(out) :: err
      integer :: ierr

      err = 0

      call remap_outfile%close_file(ierr)

      if (allocated(reduceFieldNames)) deallocate(reduceFieldNames)
      if (allocated(reduceFieldMethods)) deallocate(reduceFieldMethods)
      if (allocated(remap_accum)) deallocate(remap_accum)

   end subroutine ocn_finalize_fme_vertical_reduce!}}}

end module ocn_fme_vertical_reduce

! vim: foldmethod=marker
